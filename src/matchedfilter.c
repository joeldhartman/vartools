/*     This file is part of VARTOOLS.                                       */
/*                                                                           */
/*     VARTOOLS is free software: you can redistribute it and/or modify      */
/*     it under the terms of the GNU General Public License as published by  */
/*     the Free Software Foundation, either version 3 of the License, or     */
/*     (at your option) any later version.                                   */
/*                                                                           */
/*     This program is distributed in the hope that it will be useful,       */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/*     GNU General Public License for more details.                          */
/*                                                                           */
/*     You should have received a copy of the GNU General Public License     */
/*     along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*                                                                           */
/*     Copyright 2026  Joel Hartman                                          */
/*                                                                           */
/*  Matched-filter transient/feature search for vartools.
 *
 *  This file was developed with assistance from Claude (Anthropic) by
 *  Joel Hartman.
 *
 *  At each trial centre tau the inverse-variance matched filter fits a
 *  scaled, time-shifted template g(t - tau) plus a local constant offset c
 *  to the light curve (so the absolute brightness of the LC doesn't have
 *  to be pre-subtracted; this is the textbook generalised MF used by BLS
 *  and other transient searches):
 *      y_i ~ a * g(t_i - tau) + c + noise_i,  w_i = 1/sigma_i^2
 *      a_hat(tau) = Cov_w(y, g) / Var_w(g)
 *      SNR(tau)   = a_hat * sqrt(sum_i w_i (g_i - <g>_w)^2)
 *  where g_i(tau) = g(t_i - tau) (zero outside the caller-supplied support
 *  window) and the weighted moments are taken over the support window.
 *  The SNR is signed -- positive for bumps that match the template
 *  orientation, negative for inverted matches (e.g. a transit when the
 *  template is a positive box).  Scale-invariant in g: the template need
 *  not be pre-normalised.
 *
 *  Template sources:
 *    - Named (kind = MF_TPL_EXP / DOUBLEEXP / FLARE / GAUSS / BOX /
 *      TRIANGLE / TRAP): closed-form g(s) with at most three scalar
 *      parameters (resolved per-LC from var/expr/fixed).
 *    - File (kind = MF_TPL_FILE): 2-col ASCII (t, amplitude) loaded
 *      once at parser time; the evaluator linearly interpolates
 *      between rows and returns 0 outside the file's t range.
 *
 *  Modes: window (exact, heteroscedastic sigma) and nfft (NFFT-batched,
 *  homoscedastic sigma; see matchedfilter_nfft.c).
 *
 *  Citations:
 *    - Davenport, J. R. A., Hawley, S. L., Hebb, L. et al. 2014, ApJ, 797,
 *        122 (ADS 2014ApJ...797..122D); empirical M-dwarf flare template
 *        used by the "flare" kind.  Equation (1) of that paper.
 */

#include "commands.h"
#include "programdata.h"
#include "functions.h"

#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MF_ERROR_SCORE       (-1.0e300)
#define MF_TINY              1.0e-32

/* Load a 2-column whitespace-separated ASCII file (t, amplitude) into
 * freshly-allocated arrays.  Skips lines whose first non-whitespace
 * character is '#' as well as fully-blank lines.  Sorts the rows by t
 * ascending and drops exact-duplicate t values (keeping the first
 * occurrence).  Returns 0 on success with (*t_out, *g_out, *N_out)
 * populated; -1 on file/open/parse error (no allocation leaked). */
int mf_load_template_file(const char *path,
                          int *N_out, double **t_out, double **g_out)
{
  FILE *fp;
  char line[4096];
  int  cap = 64, n = 0;
  double *tt = NULL, *gg = NULL;

  if (path == NULL || N_out == NULL || t_out == NULL || g_out == NULL)
    return -1;
  *N_out = 0; *t_out = NULL; *g_out = NULL;
  if ((fp = fopen(path, "r")) == NULL) {
    fprintf(stderr, "-matchedfilter: cannot open template file '%s'\n", path);
    return -1;
  }
  tt = (double *) malloc(cap * sizeof(double));
  gg = (double *) malloc(cap * sizeof(double));
  if (tt == NULL || gg == NULL) { if (tt) free(tt); if (gg) free(gg); fclose(fp); return -1; }

  while (fgets(line, sizeof(line), fp)) {
    char *p = line;
    while (*p == ' ' || *p == '\t') p++;
    if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0') continue;
    double tv, gv;
    if (sscanf(p, "%lf %lf", &tv, &gv) != 2) continue;
    if (n >= cap) {
      int newcap = cap * 2;
      double *nt = (double *) realloc(tt, newcap * sizeof(double));
      double *ng = (double *) realloc(gg, newcap * sizeof(double));
      if (nt == NULL || ng == NULL) {
        if (nt == NULL) free(tt); else tt = nt;
        if (ng == NULL) free(gg); else gg = ng;
        if (nt == NULL || ng == NULL) { fclose(fp); return -1; }
      }
      tt = nt; gg = ng; cap = newcap;
    }
    tt[n] = tv; gg[n] = gv; n++;
  }
  fclose(fp);

  if (n < 2) {
    fprintf(stderr, "-matchedfilter: template file '%s' has %d valid rows (need >= 2)\n",
            path, n);
    free(tt); free(gg);
    return -1;
  }

  /* Sort by t ascending (mysort2 carries gg along). */
  mysort2(n, tt, gg);

  /* Dedupe exact-duplicate t (keep first). */
  int w = 1, r;
  for (r = 1; r < n; r++) {
    if (tt[r] != tt[r - 1]) {
      tt[w] = tt[r];
      gg[w] = gg[r];
      w++;
    }
  }
  n = w;

  if (n < 2) {
    fprintf(stderr, "-matchedfilter: template file '%s' has %d unique-t rows after dedupe (need >= 2)\n",
            path, n);
    free(tt); free(gg);
    return -1;
  }

  *N_out = n; *t_out = tt; *g_out = gg;
  return 0;
}

/* ---- Named-template evaluator -----------------------------------------
 *
 * Returns g(s) at template-relative time s = t - tau.  Returns 0.0 outside
 * |s| > support (the outer truncation window) and outside each template's
 * intrinsic support.  Peak value of every template is 1.0 by construction
 * -- amplitude scaling is handled by the inverse-variance MF formula.
 *
 * For "flare", uses Davenport+2014 eq. (1) directly.  There is a small
 * (~1%) discontinuity at t_n=0 in the published piecewise expression:
 * the rise polynomial evaluates to 1.000 at t_n=0 while the decay sum
 * evaluates to 0.6890 + 0.3030 = 0.992; this is preserved as-is.
 *
 *   exp        param[0] = tau                    : 0..support
 *   doubleexp  param[0] = tau_rise, param[1] = tau_decay : 0..support
 *   flare      param[0] = tfwhm (FWHM in time)   : -tfwhm..support
 *   gauss      param[0] = sigma                  : -support..support
 *   box        param[0] = width                  : |s| <= width/2
 *   triangle   param[0] = width                  : |s| <= width/2
 *   trap       param[0..2] = rise, flat, fall    : centred, total = rise+flat+fall
 */
/* Linear interpolation of a file-loaded template g(s) sampled on the
 * sorted nodes ctx->t[0..ctx->N-1].  Returns 0 outside the file's t
 * range AND outside |s| > support. */
static double mf_template_eval_file(const _MFEvalCtx *ctx,
                                    double support, double s)
{
  if (fabs(s) > support || ctx == NULL || ctx->N < 2) return 0.0;
  if (s <= ctx->t[0] || s >= ctx->t[ctx->N - 1]) return 0.0;
  int lo = 0, hi = ctx->N - 1;
  while (hi - lo > 1) {
    int mid = (lo + hi) >> 1;
    if (ctx->t[mid] <= s) lo = mid; else hi = mid;
  }
  double dt = ctx->t[hi] - ctx->t[lo];
  if (!(dt > 0.0)) return ctx->g[lo];
  double frac = (s - ctx->t[lo]) / dt;
  return ctx->g[lo] + frac * (ctx->g[hi] - ctx->g[lo]);
}

/* Expression-template evaluator.  Sets the user-defined time-relative
 * variable to s, then calls EvaluateExpression.  Returns 0 outside
 * |s| > support. */
static double mf_template_eval_expr(const _MFEvalCtx *ctx,
                                    double support, double s)
{
  if (fabs(s) > support || ctx == NULL ||
      ctx->expr_var == NULL || ctx->expr == NULL) return 0.0;
  SetVariable_Value_Double(ctx->lcid, ctx->threadid, 0, ctx->expr_var, s);
  return EvaluateExpression(ctx->lcid, ctx->threadid, 0, ctx->expr);
}

/* Template evaluator dispatch.  The ctx pointer carries out-of-band
 * template state when needed: file rows for MF_TPL_FILE, or a stump
 * variable + parsed expression + (lcid, threadid) for MF_TPL_EXPR.
 * For the named-template kinds ctx is unused and may be NULL. */
static double mf_template_eval(int kind, const double *param,
                               const _MFEvalCtx *ctx,
                               double support, double s)
{
  if (kind == MF_TPL_FILE)
    return mf_template_eval_file(ctx, support, s);
  if (kind == MF_TPL_EXPR)
    return mf_template_eval_expr(ctx, support, s);
  if (fabs(s) > support) return 0.0;
  switch (kind) {
    case MF_TPL_EXP: {
      double tau = param[0];
      if (!(tau > 0.0) || s < 0.0) return 0.0;
      return exp(-s / tau);
    }
    case MF_TPL_DOUBLEEXP: {
      double tr = param[0], td = param[1];
      if (!(tr > 0.0) || !(td > 0.0) || s < 0.0) return 0.0;
      /* Profile (1 - exp(-s/tr)) * exp(-s/td) has peak amplitude
       *   A = (td/(tr+td)) * (tr/(tr+td))^(tr/td)
       * at s = tr*ln(1 + td/tr).  Normalise to peak=1 so the MF SNR maps
       * cleanly to "amplitude in LC units at the peak". */
      double tpk = tr * log(1.0 + td / tr);
      double Apk = (1.0 - exp(-tpk / tr)) * exp(-tpk / td);
      if (!(Apk > MF_TINY)) return 0.0;
      return ((1.0 - exp(-s / tr)) * exp(-s / td)) / Apk;
    }
    case MF_TPL_FLARE: {
      double tfwhm = param[0];
      if (!(tfwhm > 0.0)) return 0.0;
      double tn = s / tfwhm;
      if (tn < -1.0) return 0.0;
      if (tn <= 0.0) {
        /* Rise polynomial: F(t_n) = 1 + 1.941 t_n - 0.175 t_n^2
         *                          - 2.246 t_n^3 - 1.125 t_n^4 */
        double tn2 = tn * tn;
        return 1.0 + 1.941 * tn - 0.175 * tn2
               - 2.246 * tn2 * tn - 1.125 * tn2 * tn2;
      } else {
        /* Two-component exponential decay. */
        return 0.6890 * exp(-1.600 * tn) + 0.3030 * exp(-0.2783 * tn);
      }
    }
    case MF_TPL_GAUSS: {
      double sigma = param[0];
      if (!(sigma > 0.0)) return 0.0;
      double u = s / sigma;
      return exp(-0.5 * u * u);
    }
    case MF_TPL_BOX: {
      double w = param[0];
      if (!(w > 0.0)) return 0.0;
      return (fabs(s) <= 0.5 * w) ? 1.0 : 0.0;
    }
    case MF_TPL_TRIANGLE: {
      double w = param[0];
      if (!(w > 0.0)) return 0.0;
      double a = fabs(s);
      if (a > 0.5 * w) return 0.0;
      return 1.0 - 2.0 * a / w;
    }
    case MF_TPL_TRAP: {
      double rs = param[0], fl = param[1], fa = param[2];
      if (!(rs >= 0.0) || !(fl >= 0.0) || !(fa >= 0.0)) return 0.0;
      double tot = rs + fl + fa;
      if (!(tot > 0.0)) return 0.0;
      /* Centred span [-tot/2, +tot/2]; rise occupies [-tot/2, -tot/2+rs],
       * flat [-tot/2+rs, -tot/2+rs+fl], fall [-tot/2+rs+fl, +tot/2]. */
      double a = s + 0.5 * tot;        /* shift so rise starts at 0 */
      if (a < 0.0 || a > tot) return 0.0;
      if (a < rs) return (rs > 0.0) ? (a / rs) : 1.0;
      if (a < rs + fl) return 1.0;
      double t_in_fall = a - rs - fl;
      return (fa > 0.0) ? (1.0 - t_in_fall / fa) : 0.0;
    }
  }
  return 0.0;
}

/* Intrinsic half-support of the template (above and beyond the caller's
 * support_halfwidth).  Used only to choose a sensible lower bound when
 * the caller's support is unset -- the parser requires support, so this
 * is informational only. */

/* ---- Window-mode SNR(tau) and amplitude(tau) over a tau grid ----------
 *
 * Inputs (all length N for the data, all length Ntau for the trial grid):
 *   t, mag, w[]  -- sorted by t (caller guarantees via require_sort) and
 *                   already mask/NaN/sigma-filtered.
 *   t_tau[]      -- trial-centre array (= the LC's t[] in Phase A).
 * Outputs:
 *   snr[Ntau], amp[Ntau].  Both are MF_ERROR_SCORE at trial centres
 *   whose support window contains no valid template-evaluation points
 *   or whose template is constant in-window (centred denominator <= 0).
 *
 * Model:    y_i = a * g_i + c + noise_i,    noise_i ~ N(0, sigma_i^2)
 * with the constant offset c absorbed as a nuisance parameter (a local
 * baseline within the support window).  This is the textbook generalised
 * matched filter for "signal + DC + noise" and is what flare- /
 * transit-search codes use in practice (BLS, etc.) -- without c the MF
 * picks up the absolute magnitude of the LC instead of the local
 * deviation that the user actually cares about.
 *
 *   Let Sw   = sum w_i, Swg = sum w_i g_i, Swgg = sum w_i g_i^2,
 *       Swy  = sum w_i y_i, Swyg = sum w_i y_i g_i.
 *   Centred moments (the actual fit basis):
 *       Sgg_c = Swgg - Swg^2/Sw   = sum w_i (g_i - <g>_w)^2
 *       Syg_c = Swyg - Swg*Swy/Sw = sum w_i (y_i - <y>_w)(g_i - <g>_w)
 *   Best-fit amplitude and signed SNR (inverse-variance weighted):
 *       a_hat = Syg_c / Sgg_c
 *       SNR   = a_hat * sqrt(Sgg_c)
 *   Var(a_hat) = 1 / Sgg_c when w_i = 1/sigma_i^2, so SNR is the standard
 *   amplitude-over-amplitude-uncertainty.  Scale-invariant in g.
 *
 * Cost: for each tau, a two-cursor sweep finds the support window, then
 * O(n_window) template evals and accumulations.  The cursors move
 * monotonically as tau increases (t_tau is sorted in Phase A), so the
 * total work across all taus is bounded by O(N) plus O(sum n_window). */
static void mf_window_periodogram(int N, const double *t,
                                  const double *mag, const double *w,
                                  int kind, const double *params,
                                  const _MFEvalCtx *ctx, double support,
                                  int Ntau, const double *t_tau,
                                  double *snr_out, double *amp_out)
{
  int i, lo, hi;
  lo = 0; hi = 0;
  for (i = 0; i < Ntau; i++) {
    double tau = t_tau[i];
    double tlo = tau - support, thi = tau + support;
    while (lo < N && t[lo] < tlo) lo++;
    if (hi < lo) hi = lo;
    while (hi < N && t[hi] <= thi) hi++;
    double Sw = 0.0, Swg = 0.0, Swgg = 0.0, Swy = 0.0, Swyg = 0.0;
    int j;
    for (j = lo; j < hi; j++) {
      double g = mf_template_eval(kind, params, ctx, support, t[j] - tau);
      /* All points in the support window contribute to Sw / Swy (they
       * inform the local baseline c).  Points where g==0 still pin the
       * baseline; only Swg / Swgg / Swyg get the g-weighted forms. */
      double wi = w[j];
      Sw  += wi;
      Swy += wi * mag[j];
      if (g != 0.0) {
        double wg = wi * g;
        Swg  += wg;
        Swgg += wg * g;
        Swyg += wg * mag[j];
      }
    }
    if (Sw > MF_TINY) {
      double Sgg_c = Swgg - (Swg * Swg) / Sw;
      double Syg_c = Swyg - (Swg * Swy) / Sw;
      if (Sgg_c > MF_TINY) {
        double a = Syg_c / Sgg_c;
        amp_out[i] = a;
        snr_out[i] = a * sqrt(Sgg_c);
        continue;
      }
    }
    amp_out[i] = 0.0;
    snr_out[i] = MF_ERROR_SCORE;
  }
}

/* ---- Sign-filter helper: returns the value to compare for ranking,
 * or MF_ERROR_SCORE if the sample is excluded by the signs policy. */
static double mf_sign_rank(double snr, int signs)
{
  if (snr == MF_ERROR_SCORE) return MF_ERROR_SCORE;
  switch (signs) {
    case MF_SIGNS_POSITIVE: return (snr > 0.0) ? snr : MF_ERROR_SCORE;
    case MF_SIGNS_NEGATIVE: return (snr < 0.0) ? (-snr) : MF_ERROR_SCORE;
    case MF_SIGNS_BOTH:
    default:                return fabs(snr);
  }
}

/* ---- Mean / RMS of SNR(tau) over the sign-filtered trial grid ---------- */
static void mf_compute_mean_rms(int Ntau, const double *snr, int signs,
                                double *mean_out, double *rms_out)
{
  long double Sum = 0.0L, Sumsqr = 0.0L;
  int i, Ngood = 0;
  for (i = 0; i < Ntau; i++) {
    double rk = mf_sign_rank(snr[i], signs);
    if (rk == MF_ERROR_SCORE) continue;
    /* Mean / RMS are taken over the SIGNED SNR (so a negative-signs run
     * gives a negative mean if the data are dominated by below-zero dips).
     * mf_sign_rank only gates whether the sample contributes. */
    Sum    += (long double) snr[i];
    Sumsqr += (long double) snr[i] * (long double) snr[i];
    Ngood++;
  }
  if (Ngood > 0) {
    Sum /= Ngood; Sumsqr /= Ngood;
    *mean_out = (double) Sum;
    double v = (double)(Sumsqr - Sum * Sum);
    *rms_out  = (v > 0.0) ? sqrt(v) : 0.0;
  } else {
    *mean_out = 0.0;
    *rms_out  = 0.0;
  }
}

/* ---- Iterative peak finder with min_separation masking ----------------
 *
 * Largest |SNR| (after sign filter) wins; surrounding tau within
 * min_sep is masked out; repeat until Npeaks distinct peaks are found
 * or the grid is exhausted.  Output arrays are pre-allocated by the
 * caller; unused slots are filled with sentinels.
 */
static void mf_find_peaks(int Ntau, const double *t_tau,
                          const double *snr, const double *amp,
                          int signs, double min_sep, int Npeaks,
                          double *peakTime, double *peakSNR, double *peakAmp)
{
  unsigned char *masked = (unsigned char *) calloc(Ntau, sizeof(unsigned char));
  if (masked == NULL) vt_error(ERR_MEMALLOC);
  int k;
  for (k = 0; k < Npeaks; k++) {
    peakTime[k] = 0.0;
    peakSNR[k]  = 0.0;
    peakAmp[k]  = 0.0;
  }
  for (k = 0; k < Npeaks; k++) {
    int best = -1;
    double best_rk = -1.0;
    int i;
    for (i = 0; i < Ntau; i++) {
      if (masked[i]) continue;
      double rk = mf_sign_rank(snr[i], signs);
      if (rk == MF_ERROR_SCORE) continue;
      if (rk > best_rk) { best_rk = rk; best = i; }
    }
    if (best < 0) {
      /* no more accepted samples; remaining slots stay at sentinels */
      break;
    }
    peakTime[k] = t_tau[best];
    peakSNR[k]  = snr[best];
    peakAmp[k]  = amp[best];
    /* Mask everything within +/- min_sep of the best tau. */
    double tau_lo = peakTime[k] - min_sep;
    double tau_hi = peakTime[k] + min_sep;
    for (i = 0; i < Ntau; i++) {
      if (t_tau[i] >= tau_lo && t_tau[i] <= tau_hi) masked[i] = 1;
    }
  }
  free(masked);
}

/* ---- Subtract one fitted template from the LC (in-place on mag_loc) --- */
static void mf_subtract_template(int N, const double *t, double *mag_loc,
                                 int kind, const double *params,
                                 const _MFEvalCtx *ctx, double support,
                                 double tau, double amp)
{
  int j;
  for (j = 0; j < N; j++) {
    double s = t[j] - tau;
    if (fabs(s) > support) continue;
    double g = mf_template_eval(kind, params, ctx, support, s);
    if (g != 0.0) mag_loc[j] -= amp * g;
  }
}

/* ---- Aux file output (ASCII when p->ascii, else binary) --------------- */
static void mf_write_auxfile(const char *path, int ascii, int Ntau,
                             const double *t_tau, const double *snr,
                             const double *amp)
{
  FILE *fp = fopen(path, ascii ? "w" : "wb");
  if (fp == NULL) {
    fprintf(stderr, "Cannot write matched-filter aux file %s\n", path);
    exit(3);
  }
  if (ascii) {
    fprintf(fp, "# t  SNR  amplitude\n");
    int i;
    for (i = 0; i < Ntau; i++) {
      if (snr[i] == MF_ERROR_SCORE) continue;
      fprintf(fp, "%f %g %g\n", t_tau[i], snr[i], amp[i]);
    }
  } else {
    /* Header-less compact binary: int32 N, then double[N] tau, SNR, amp.
     * MF_ERROR_SCORE samples are written through as-is so the layout is
     * fixed and the caller can mask on the SNR sentinel. */
    int32_t n32 = (int32_t) Ntau;
    fwrite(&n32, sizeof(int32_t), 1, fp);
    fwrite(t_tau, sizeof(double), Ntau, fp);
    fwrite(snr,   sizeof(double), Ntau, fp);
    fwrite(amp,   sizeof(double), Ntau, fp);
  }
  fclose(fp);
}

/* Forward decl for the NFFT path (defined in matchedfilter_nfft.c).
 * Returns 0 on success (snr_out / amp_out populated), nonzero on
 * failure (caller falls back to window mode). */
extern int mf_nfft_periodogram(int N, const double *t, const double *mag,
                               const double *w_in,
                               int kind, const double *params,
                               const _MFEvalCtx *ctx, double support,
                               int Ntau, const double *t_tau,
                               double *snr_out, double *amp_out);

/* Dispatch between window-mode and NFFT-mode periodogram evaluators
 * based on mf->mode.  NFFT-mode falls through to window mode if the
 * NFFT routine returns nonzero (e.g. N_nfft cap exceeded). */
static void mf_periodogram_dispatch(const _MatchedFilter *mf,
                                    int N, const double *t, const double *mag,
                                    const double *w,
                                    const double *params,
                                    const _MFEvalCtx *ctx, double support,
                                    int Ntau, const double *t_tau,
                                    double *snr_out, double *amp_out,
                                    void (*windowfn)(int, const double *,
                                                     const double *, const double *,
                                                     int, const double *,
                                                     const _MFEvalCtx *, double,
                                                     int, const double *,
                                                     double *, double *))
{
  if (mf->mode == MF_MODE_NFFT) {
    int rc = mf_nfft_periodogram(N, t, mag, w, mf->kind, params, ctx, support,
                                 Ntau, t_tau, snr_out, amp_out);
    if (rc == 0) return;
    /* Fall through to window mode on failure. */
  }
  windowfn(N, t, mag, w, mf->kind, params, ctx, support, Ntau, t_tau,
           snr_out, amp_out);
}

/* ---- Resolve one _MFScalar to its per-LC scalar value ----------------- */
static double mf_resolve_scalar(_MFScalar *s, int lcnum, int lc_name_num)
{
  double v;
  switch (s->source) {
    case VARTOOLS_SOURCE_EVALEXPRESSION:
      v = EvaluateExpression(lc_name_num, lcnum, 0, s->expr); break;
    case VARTOOLS_SOURCE_EXISTINGVARIABLE:
      v = EvaluateVariable_Double(lc_name_num, lcnum, 0, s->var); break;
    case VARTOOLS_SOURCE_FIXED:
    default:
      v = s->fixed; break;
  }
  if (s->vals != NULL) s->vals[lcnum] = v;
  return v;
}

/* ---- Expression-template setup (MF_TPL_EXPR) -------------------------- */

/* Recurse through a parsed expression tree and return 1 if any
 * referenced variable has vectortype VARTOOLS_VECTORTYPE_LC.  Mirrors
 * fourierfilter.c::_fourierfilter_expr_has_lc_var; rejecting LC vectors
 * keeps the user from accidentally referencing the per-observation
 * t/mag/err arrays in a template expression that's supposed to be a
 * function of the template-relative time variable only. */
static int mf_expr_has_lc_var(_Expression *e, const char **badname)
{
  if (e == NULL) return 0;
  if (e->op1type == VARTOOLS_OPERANDTYPE_VARIABLE && e->op1_variable != NULL) {
    if (e->op1_variable->vectortype == VARTOOLS_VECTORTYPE_LC) {
      if (badname) *badname = e->op1_variable->varname;
      return 1;
    }
  }
  if (e->op2type == VARTOOLS_OPERANDTYPE_VARIABLE && e->op2_variable != NULL) {
    if (e->op2_variable->vectortype == VARTOOLS_VECTORTYPE_LC) {
      if (badname) *badname = e->op2_variable->varname;
      return 1;
    }
  }
  if (e->op1type == VARTOOLS_OPERANDTYPE_EXPRESSION
      && mf_expr_has_lc_var((_Expression *) e->op1_expression, badname))
    return 1;
  if (e->op2type == VARTOOLS_OPERANDTYPE_EXPRESSION
      && mf_expr_has_lc_var((_Expression *) e->op2_expression, badname))
    return 1;
  if (e->op1type == VARTOOLS_OPERANDTYPE_FUNCTION
      && e->op1_functioncall != NULL) {
    _FunctionCall *fc = (_FunctionCall *) e->op1_functioncall;
    int a;
    for (a = 0; a < fc->Nexpr; a++)
      if (mf_expr_has_lc_var(fc->arguments[a], badname)) return 1;
  }
  if (e->op2type == VARTOOLS_OPERANDTYPE_FUNCTION
      && e->op2_functioncall != NULL) {
    _FunctionCall *fc = (_FunctionCall *) e->op2_functioncall;
    int a;
    for (a = 0; a < fc->Nexpr; a++)
      if (mf_expr_has_lc_var(fc->arguments[a], badname)) return 1;
  }
  return 0;
}

/* Whole-word in-place substring replacement on src, writing into dst
 * (length dst_sz).  Only identifier characters (alnum or '_') on the
 * boundaries count as part of a word.  Returns 0 on success, -1 if dst
 * would overflow.  Mirrors fourierfilter.c's helper. */
static int mf_subst_whole_word(const char *src, const char *from,
                               const char *to, char *dst, size_t dst_sz)
{
  size_t flen = strlen(from), tlen = strlen(to);
  size_t wr = 0;
  const char *p = src;
  while (*p) {
    if (strncmp(p, from, flen) == 0) {
      int before_ok = (p == src) || !(isalnum((unsigned char) *(p-1)) || *(p-1) == '_');
      int after_ok  = (p[flen] == '\0') || !(isalnum((unsigned char) p[flen]) || p[flen] == '_');
      if (before_ok && after_ok) {
        if (wr + tlen >= dst_sz) return -1;
        memcpy(dst + wr, to, tlen);
        wr += tlen;
        p += flen;
        continue;
      }
    }
    if (wr + 1 >= dst_sz) return -1;
    dst[wr++] = *p++;
  }
  dst[wr] = '\0';
  return 0;
}

/* Called from analytic.c::CompileAllExpressions once per -matchedfilter
 * command with a template-expr string.  Allocates the stump
 * INTERNALSCALAR variable for the time-relative coordinate, substitutes
 * the user's name with the stump in the expression source, parses, and
 * validates that the expression depends only on per-star scalars,
 * constants, and the stump variable (no LC-vector references). */
void SetupMatchedFilterExpression(ProgramData *p, _MatchedFilter *mf, int cnum)
{
  char stump[64];
  char substituted[4 * MAXLEN];
  const char *badname = NULL;

  if (mf == NULL || mf->expr_string == NULL || mf->expr_varname == NULL) return;

  snprintf(stump, sizeof(stump), "__mf_s_%d", cnum);

  mf->expr_var = CreateVariable(p, stump, VARTOOLS_TYPE_DOUBLE,
                                VARTOOLS_VECTORTYPE_INTERNALSCALAR, NULL);
  RegisterDataFromLightCurve(p, mf->expr_var->dataptr, VARTOOLS_TYPE_DOUBLE,
                             0, 0, 0, 0, 0, NULL, mf->expr_var, -1, stump);

  if (mf_subst_whole_word(mf->expr_string, mf->expr_varname, stump,
                          substituted, sizeof(substituted))) {
    vt_error2(ERR_CODEERROR,
              "-matchedfilter template expr: substituted expression exceeds internal buffer");
  }

  mf->expr = ParseExpression(substituted, p);

  if (mf_expr_has_lc_var(mf->expr, &badname)) {
    fprintf(stderr,
            "\nError in -matchedfilter (command #%d): template expr may not "
            "reference light-curve-vector variables.\n"
            "'%s' has vectortype VECTORTYPE_LC.  Only per-star scalars, constants, "
            "and the template-relative variable (%s) are allowed.\n\n",
            cnum, badname ? badname : "<unknown>", mf->expr_varname);
    vt_error(ERR_CODEERROR);
  }
}

/* ---- Entry point invoked from processcommand.c ------------------------ */
void RunMatchedFilterCommand(ProgramData *p, Command *c, _MatchedFilter *mf,
                             int lcnum, int lc_name_num, int thisindex)
{
  int Nraw = p->NJD[lcnum];
  double *t_in = p->t[lcnum];
  double *mag_in = p->mag[lcnum];
  double *sig_in = p->sig[lcnum];
  double *t = NULL, *mag = NULL, *w = NULL, *mag_orig = NULL;
  double *snr = NULL, *amp = NULL;
  int Nused = 0;
  int k, j;
  char outname[MAXLEN];

  /* Sentinel-fill outputs in case we early-exit. */
  for (k = 0; k < mf->Npeaks; k++) {
    mf->peaktimes[lcnum][k] = 0.0;
    mf->peakSNR[lcnum][k]   = 0.0;
    mf->peakamps[lcnum][k]  = 0.0;
  }
  mf->mean_snr[lcnum] = 0.0;
  mf->rms_snr[lcnum]  = 0.0;

  if (Nraw <= 1) return;

  /* Resolve all var/expr/fixed scalars for this LC. */
  double params[3] = {0.0, 0.0, 0.0};
  for (k = 0; k < mf->nparams; k++)
    params[k] = mf_resolve_scalar(&mf->p[k], lcnum, lc_name_num);
  /* Evaluator context, NULL for named-template kinds.  Carries:
   *   - File rows (kind == MF_TPL_FILE; loaded once at parser time)
   *   - Stump variable + parsed expression + lcid/threadid for
   *     SetVariable_Value_Double + EvaluateExpression
   *     (kind == MF_TPL_EXPR; compiled in analytic.c). */
  _MFEvalCtx ctx_local;
  const _MFEvalCtx *ctx = NULL;
  memset(&ctx_local, 0, sizeof(ctx_local));
  if (mf->kind == MF_TPL_FILE && mf->file_t != NULL && mf->file_N >= 2) {
    ctx_local.N = mf->file_N;
    ctx_local.t = mf->file_t;
    ctx_local.g = mf->file_g;
    ctx = &ctx_local;
  } else if (mf->kind == MF_TPL_EXPR &&
             mf->expr_var != NULL && mf->expr != NULL) {
    ctx_local.expr_var = mf->expr_var;
    ctx_local.expr     = mf->expr;
    ctx_local.lcid     = lc_name_num;
    ctx_local.threadid = lcnum;
    ctx = &ctx_local;
  }
  double support = mf_resolve_scalar(&mf->support, lcnum, lc_name_num);
  double min_sep = mf->min_sep_given
                 ? mf_resolve_scalar(&mf->min_sep, lcnum, lc_name_num)
                 : support;
  if (!(support > 0.0)) return;
  if (!(min_sep >= 0.0)) min_sep = support;

  /* Build the filtered (t, mag, w) arrays. */
  t   = (double *) malloc(Nraw * sizeof(double));
  mag = (double *) malloc(Nraw * sizeof(double));
  w   = (double *) malloc(Nraw * sizeof(double));
  if (t == NULL || mag == NULL || w == NULL) vt_error(ERR_MEMALLOC);
  for (j = 0; j < Nraw; j++) {
    if (isnan(mag_in[j])) continue;
    if (mf->usemask && mf->maskvar != NULL) {
      if (!(EvaluateVariable_Double(lc_name_num, lcnum, j, mf->maskvar) > VARTOOLS_MASK_TINY))
        continue;
    }
    /* Weighting: 1/sigma^2 if sigma is positive, else uniform (matches
     * pyvartools' convention for non-noerr commands when sigma is missing).
     * Skip points with sigma<=0 or NaN to be safe with downstream math. */
    if (!(sig_in[j] > 0.0) || isnan(sig_in[j])) {
      w[Nused] = 1.0;
    } else {
      w[Nused] = 1.0 / (sig_in[j] * sig_in[j]);
    }
    t[Nused]   = t_in[j];
    mag[Nused] = mag_in[j];
    Nused++;
  }
  if (Nused < 2) {
    free(t); free(mag); free(w);
    return;
  }

  snr = (double *) malloc(Nused * sizeof(double));
  amp = (double *) malloc(Nused * sizeof(double));
  if (snr == NULL || amp == NULL) vt_error(ERR_MEMALLOC);

  if (mf->whiten) {
    /* Whiten branch: each peak iter recomputes SNR(tau) from the working
     * copy of mag, finds the single best tau (sign-filtered), subtracts
     * a_hat * g(t - tau) from the working copy, and proceeds.  Mean/RMS
     * are taken from cycle 0's SNR(tau) for parity with the non-whiten
     * branch (downstream tools that interpret the noise estimate as the
     * baseline noise floor get the same number in either branch). */
    mag_orig = (double *) malloc(Nused * sizeof(double));
    if (mag_orig == NULL) vt_error(ERR_MEMALLOC);
    memcpy(mag_orig, mag, Nused * sizeof(double));

    double cycle0_mean = 0.0, cycle0_rms = 0.0;
    int didcycle0 = 0;
    for (k = 0; k < mf->Npeaks; k++) {
      mf_periodogram_dispatch(mf, Nused, t, mag, w, params, ctx, support,
                              Nused, t, snr, amp, mf_window_periodogram);
      if (!didcycle0) {
        mf_compute_mean_rms(Nused, snr, mf->signs, &cycle0_mean, &cycle0_rms);
        didcycle0 = 1;
      }
      /* Find the single best (sign-filtered) tau on this cycle. */
      int best = -1;
      double best_rk = -1.0;
      for (j = 0; j < Nused; j++) {
        double rk = mf_sign_rank(snr[j], mf->signs);
        if (rk == MF_ERROR_SCORE) continue;
        if (rk > best_rk) { best_rk = rk; best = j; }
      }
      if (best < 0) break;
      mf->peaktimes[lcnum][k] = t[best];
      mf->peakSNR[lcnum][k]   = snr[best];
      mf->peakamps[lcnum][k]  = amp[best];
      /* Subtract this peak's contribution and continue. */
      mf_subtract_template(Nused, t, mag, mf->kind, params, ctx, support,
                           t[best], amp[best]);
    }
    /* Aux dump uses the cycle-0 (un-whitened) SNR(tau) so downstream
     * users see the baseline match-statistic surface. */
    if (mf->omatchfile) {
      memcpy(mag, mag_orig, Nused * sizeof(double));
      mf_periodogram_dispatch(mf, Nused, t, mag, w, params, ctx, support,
                              Nused, t, snr, amp, mf_window_periodogram);
      int i1, i2 = 0;
      for (i1 = 0; p->lcnames[lc_name_num][i1] != '\0'; i1++)
        if (p->lcnames[lc_name_num][i1] == '/') i2 = i1 + 1;
      sprintf(outname, "%s/%s%s", mf->outdir,
               &p->lcnames[lc_name_num][i2], mf->suffix);
      mf_write_auxfile(outname, p->ascii, Nused, t, snr, amp);
    }
    mf->mean_snr[lcnum] = cycle0_mean;
    mf->rms_snr[lcnum]  = cycle0_rms;
    free(mag_orig);
  } else {
    mf_periodogram_dispatch(mf, Nused, t, mag, w, params, ctx, support,
                            Nused, t, snr, amp, mf_window_periodogram);
    mf_compute_mean_rms(Nused, snr, mf->signs, &mf->mean_snr[lcnum],
                        &mf->rms_snr[lcnum]);
    mf_find_peaks(Nused, t, snr, amp, mf->signs, min_sep, mf->Npeaks,
                  mf->peaktimes[lcnum], mf->peakSNR[lcnum],
                  mf->peakamps[lcnum]);
    if (mf->omatchfile) {
      int i1, i2 = 0;
      for (i1 = 0; p->lcnames[lc_name_num][i1] != '\0'; i1++)
        if (p->lcnames[lc_name_num][i1] == '/') i2 = i1 + 1;
      sprintf(outname, "%s/%s%s", mf->outdir,
               &p->lcnames[lc_name_num][i2], mf->suffix);
      mf_write_auxfile(outname, p->ascii, Nused, t, snr, amp);
    }
  }

  free(t); free(mag); free(w);
  free(snr); free(amp);
}
