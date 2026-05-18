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
/*     Copyright 2026  Joel Hartman                                          */
/*                                                                           */
/*  NFFT-batched matched filter for vartools (-matchedfilter mode nfft).
 *
 *  Algorithm: at each trial centre tau (= the LC's time array in Phase A)
 *  the matched filter needs five weighted sums over the support window:
 *
 *      Sw(tau)   = sum_i w_i * B(t_i - tau)
 *      Wg(tau)   = sum_i w_i * g(t_i - tau)
 *      Wgg(tau)  = sum_i w_i * g(t_i - tau)^2
 *      Wy(tau)   = sum_i w_i * y_i * B(t_i - tau)
 *      Wyg(tau)  = sum_i w_i * y_i * g(t_i - tau)
 *
 *  where B(s) = (|s| <= support) is the box indicator of the support
 *  window and g(s) is the user-selected template kernel.  Each is a
 *  convolution f(tau) = sum_i alpha_i * K(t_i - tau) with two distinct
 *  alphas (w, w*y) and three distinct kernels (B, g, g^2).  We compute
 *  all five via the NFFT-batched scheme:
 *
 *    1. Embed the data times t_i in a periodic domain of length
 *       T_window = (tmax - tmin) + 2.5 * support so that no support
 *       window wraps to the opposite side.
 *    2. Adjoint NFFT of (alpha, t) -> rho_hat[k]  for alpha = w
 *       and again for alpha = w*y.   (2 adjoint NFFTs)
 *    3. Sample K_B, K_g, K_gg on the uniform wrap-around grid
 *       s_m = m * T_window/N_nfft (negative s wraps to s + T_window),
 *       FFT each via FFTW.   (3 forward FFTs)
 *    4. For each (alpha, K) pair: Q[k] = conj(K_hat[k]) * rho_hat[k];
 *       forward NFFT to tau_j (= t_j in Phase A) gives sum(tau_j).
 *       (5 forward NFFTs)
 *
 *  Cost: O(N_nfft * log(N_nfft) + M) per NFFT, with N_nfft chosen so
 *  the kernel grid step dt = T_window/N_nfft is at most support/16 to
 *  resolve the template shape.
 *
 *  Limitations vs window-mode:
 *
 *  - HOMOSCEDASTIC sigma: a single weight w = 1/sigma_median^2 is used
 *    for every point.  Mixing per-point weights would require one more
 *    adjoint NFFT per template kernel and is not worth the cost for
 *    Phase B.1.
 *  - SHARP-EDGED TEMPLATES (box, triangle, trap) develop a few-percent
 *    spectral-leakage artefact near the support boundary because the
 *    box indicator B has a discontinuity that no finite-N_nfft FFT can
 *    represent exactly.  Smooth templates (gauss, exp, doubleexp,
 *    flare) are unaffected.  Use mode = window for the cleanest
 *    semantics when sharp edges matter.
 *
 *  Returns 0 on success (caller's snr_out / amp_out filled), -1 on
 *  failure (caller should fall back to window mode).  Failure modes:
 *  insufficient data, support too wide vs baseline, NFFT plan setup
 *  error.
 */

#include "commands.h"
#include "programdata.h"
#include "functions.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_NFFT3
#include <fftw3.h>
#include <nfft3.h>
#endif

#define MF_NFFT_ERROR_SCORE   (-1.0e300)
#define MF_NFFT_TINY          1.0e-32

/* File-template linear-interp evaluator (mirrors matchedfilter.c's
 * mf_template_eval_file).  Returns 0 outside the file's t range AND
 * outside |s| > support. */
static double mf_nfft_template_file(const _MFEvalCtx *ctx,
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

/* Expression-template evaluator (mirrors matchedfilter.c's
 * mf_template_eval_expr).  Sets the user-defined time-relative variable
 * to s, then calls EvaluateExpression. */
static double mf_nfft_template_expr(const _MFEvalCtx *ctx,
                                    double support, double s)
{
  if (fabs(s) > support || ctx == NULL ||
      ctx->expr_var == NULL || ctx->expr == NULL) return 0.0;
  SetVariable_Value_Double(ctx->lcid, ctx->threadid, 0, ctx->expr_var, s);
  return EvaluateExpression(ctx->lcid, ctx->threadid, 0, ctx->expr);
}

/* Template evaluator (mirrors matchedfilter.c's mf_template_eval, kept
 * file-local here to avoid an inter-file include).  Dispatch by kind:
 * MF_TPL_FILE / MF_TPL_EXPR use the ctx out-of-band state, named kinds
 * use param[]. */
static double mf_nfft_template(int kind, const double *param,
                               const _MFEvalCtx *ctx,
                               double support, double s)
{
  if (kind == MF_TPL_FILE)
    return mf_nfft_template_file(ctx, support, s);
  if (kind == MF_TPL_EXPR)
    return mf_nfft_template_expr(ctx, support, s);
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
      double tpk = tr * log(1.0 + td / tr);
      double Apk = (1.0 - exp(-tpk / tr)) * exp(-tpk / td);
      if (!(Apk > MF_NFFT_TINY)) return 0.0;
      return ((1.0 - exp(-s / tr)) * exp(-s / td)) / Apk;
    }
    case MF_TPL_FLARE: {
      double tfwhm = param[0];
      if (!(tfwhm > 0.0)) return 0.0;
      double tn = s / tfwhm;
      if (tn < -1.0) return 0.0;
      if (tn <= 0.0) {
        double tn2 = tn * tn;
        return 1.0 + 1.941 * tn - 0.175 * tn2
               - 2.246 * tn2 * tn - 1.125 * tn2 * tn2;
      } else {
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
      double a = s + 0.5 * tot;
      if (a < 0.0 || a > tot) return 0.0;
      if (a < rs) return (rs > 0.0) ? (a / rs) : 1.0;
      if (a < rs + fl) return 1.0;
      double t_in_fall = a - rs - fl;
      return (fa > 0.0) ? (1.0 - t_in_fall / fa) : 0.0;
    }
  }
  return 0.0;
}

#ifdef HAVE_NFFT3

/* Round up to next power of 2 (or to the input if already a power of 2). */
static int mf_next_pow2(int n)
{
  int r = 1;
  while (r < n) r <<= 1;
  return r;
}

/* Sample a kernel K(s) on the FFT's wrap-around grid s_m = m * dt for
 * m = 0..N-1, with negative s_phys wrapping to s_phys + T_window.
 *
 * The output is multiplied by (-1)^m so that the subsequent FFTW3
 * FORWARD transform produces K_hat in NFFT3's CENTERED f_hat layout
 * (DC at index N/2, mode -N/2 at index 0, mode +N/2-1 at index N-1).
 * Without this prefactor, FFTW's output would be in "natural" layout
 * (DC at index 0) and would not align with NFFT3's rho_hat layout
 * when we multiply the two in Fourier space.
 *
 * Writes complex output into kbuf[N].  kernel_kind: 0 = B(s) (box
 * indicator), 1 = g(s), 2 = g(s)^2. */
static void mf_nfft_sample_kernel(int kind, const double *params,
                                  const _MFEvalCtx *ctx,
                                  double support, int kernel_kind,
                                  int N, double dt, double T_window,
                                  fftw_complex *kbuf)
{
  int m;
  for (m = 0; m < N; m++) {
    double s = (double) m * dt;
    if (s > 0.5 * T_window) s -= T_window;
    double K = 0.0;
    if (fabs(s) <= support) {
      switch (kernel_kind) {
        case 0: K = 1.0; break;
        case 1: K = mf_nfft_template(kind, params, ctx, support, s); break;
        case 2: {
          double g = mf_nfft_template(kind, params, ctx, support, s);
          K = g * g;
          break;
        }
      }
    }
    double sign = (m & 1) ? -1.0 : 1.0;   /* (-1)^m fftshift trick */
    kbuf[m][0] = sign * K;
    kbuf[m][1] = 0.0;
  }
}

/* In-place complex FFT (forward, FFTW3 sign convention exp(-2pi i k m / N)). */
static void mf_nfft_fftw_forward(int N, fftw_complex *buf)
{
  fftw_plan plan = fftw_plan_dft_1d(N, buf, buf, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
}

/* Build the per-tau sum f(tau_j) = sum_i alpha_i * K(t_i - tau_j) by
 * the kernel-FFT * NFFT scheme.
 *
 * Derivation: with NFFT3 conventions (adjoint uses exp(+2pi i k x_j),
 * trafo uses exp(-2pi i k x_j)), and FFTW3 FORWARD on the (-1)^m-
 * shifted kernel samples putting K_hat in CENTERED layout matching
 * NFFT3's f_hat, the convolution f(tau) = sum_i alpha_i K(t_i - tau)
 * becomes
 *     f(tau_j) = (1/N) * sum_k K_hat[k] * rho_hat[k] * exp(-2pi i k x_j)
 * which is exactly what nfft_trafo(plan) computes if we load
 * plan->f_hat[k] = K_hat[k] * rho_hat[k] / N.  No conjugation needed
 * (the sign convention is built into the centered layout).
 *
 * out_sum[] (length Ntau) is written. */
static void mf_nfft_evaluate(nfft_plan *plan, int N_nfft,
                             const fftw_complex *rho_hat,
                             const fftw_complex *K_hat,
                             double *out_sum, int Ntau)
{
  int k, j;
  double inv_N = 1.0 / (double) N_nfft;
  for (k = 0; k < N_nfft; k++) {
    double ar = K_hat[k][0],   ai = K_hat[k][1];
    double br = rho_hat[k][0], bi = rho_hat[k][1];
    plan->f_hat[k][0] = (ar * br - ai * bi) * inv_N;
    plan->f_hat[k][1] = (ar * bi + ai * br) * inv_N;
  }
  nfft_trafo(plan);
  for (j = 0; j < Ntau; j++) {
    out_sum[j] = plan->f[j][0];
  }
}

/* Entry point.  All arrays sized as documented; t_tau[] is the trial
 * grid (= t[] in Phase A) and must be sorted ascending.  Returns 0 on
 * success, -1 on failure (caller should fall back to window mode). */
int mf_nfft_periodogram(int N, const double *t, const double *mag,
                        const double *w_in,
                        int kind, const double *params,
                        const _MFEvalCtx *ctx, double support,
                        int Ntau, const double *t_tau,
                        double *snr_out, double *amp_out)
{
  int j, ret = -1;
  if (N < 2 || Ntau < 1 || !(support > 0.0)) return -1;

  /* Time range and t_mid. */
  double tmin = t[0], tmax = t[0];
  for (j = 1; j < N; j++) {
    if (t[j] < tmin) tmin = t[j];
    if (t[j] > tmax) tmax = t[j];
  }
  double baseline = tmax - tmin;
  if (!(baseline > 0.0)) return -1;
  /* T_window: enough headroom so support windows never wrap.  The 2.5x
   * factor gives 0.25 * support of margin between the data range and
   * the periodic boundary, well past any kernel support. */
  double T_window = baseline + 2.5 * support;
  double t_mid    = 0.5 * (tmin + tmax);

  /* NFFT bin count: dt_grid = T_window / N_nfft should be at most
   * support / 16 to resolve template shape with low aliasing. */
  double dt_target = support / 16.0;
  int N_nfft = mf_next_pow2((int) ceil(T_window / dt_target));
  if (N_nfft < 4096) N_nfft = 4096;
  /* Safety cap.  At N_nfft = 2^22 each complex buffer is 64 MB; we use
   * ~6 buffers plus the NFFT plan's internal scratch, so a 2^22 cap is
   * ~600 MB peak.  Beyond that the caller should use mode window. */
  if (N_nfft > (1 << 22)) return -1;
  double dt_grid = T_window / (double) N_nfft;

  /* Homoscedastic weight: median of (w_in if positive, skip otherwise).
   * Fall back to w = 1 if no positive weights exist. */
  double w_const = 1.0;
  {
    double *wpos = (double *) malloc(N * sizeof(double));
    if (wpos == NULL) return -1;
    int npos = 0;
    for (j = 0; j < N; j++) if (w_in[j] > 0.0) wpos[npos++] = w_in[j];
    if (npos > 0) w_const = median(npos, wpos);
    free(wpos);
    if (!(w_const > 0.0)) w_const = 1.0;
  }

  /* NFFT plan setup: M = N data nodes, F_HAT = N_nfft modes.  We set
   * the nodes once and reuse the same plan for both adjoint and forward
   * trafo, since tau_j == t_j in Phase A. */
  nfft_plan plan;
  nfft_init_1d(&plan, N_nfft, N);
  double inv_T = 1.0 / T_window;
  for (j = 0; j < N; j++) {
    double xj = (t[j] - t_mid) * inv_T;
    /* The 2.5*support margin guarantees |xj| < 0.5, but clamp defensively
     * in case of pathological data. */
    if (xj < -0.49999) xj = -0.49999;
    if (xj >  0.49999) xj =  0.49999;
    plan.x[j] = xj;
  }
  nfft_precompute_one_psi(&plan);

  /* Adjoint NFFTs: rho_hat_w[k] = sum_j w * exp(-2pi i k x_j),
   *               rho_hat_wy[k] = sum_j w * y_j * exp(-2pi i k x_j). */
  fftw_complex *rho_w  = (fftw_complex *) fftw_malloc(N_nfft * sizeof(fftw_complex));
  fftw_complex *rho_wy = (fftw_complex *) fftw_malloc(N_nfft * sizeof(fftw_complex));
  fftw_complex *K_B    = (fftw_complex *) fftw_malloc(N_nfft * sizeof(fftw_complex));
  fftw_complex *K_g    = (fftw_complex *) fftw_malloc(N_nfft * sizeof(fftw_complex));
  fftw_complex *K_gg   = (fftw_complex *) fftw_malloc(N_nfft * sizeof(fftw_complex));
  if (rho_w == NULL || rho_wy == NULL ||
      K_B == NULL || K_g == NULL || K_gg == NULL) {
    if (rho_w)  fftw_free(rho_w);
    if (rho_wy) fftw_free(rho_wy);
    if (K_B)    fftw_free(K_B);
    if (K_g)    fftw_free(K_g);
    if (K_gg)   fftw_free(K_gg);
    nfft_finalize(&plan);
    return -1;
  }

  /* Adjoint for alpha = w (constant; could simplify by w * sampling-FFT,
   * but the path is symmetric with the w*y case and the extra cost is
   * trivial). */
  for (j = 0; j < N; j++) { plan.f[j][0] = w_const;            plan.f[j][1] = 0.0; }
  nfft_adjoint(&plan);
  memcpy(rho_w, plan.f_hat, N_nfft * sizeof(fftw_complex));

  /* Adjoint for alpha = w * y. */
  for (j = 0; j < N; j++) { plan.f[j][0] = w_const * mag[j];   plan.f[j][1] = 0.0; }
  nfft_adjoint(&plan);
  memcpy(rho_wy, plan.f_hat, N_nfft * sizeof(fftw_complex));

  /* Sample and FFT the three kernels. */
  mf_nfft_sample_kernel(kind, params, ctx, support, 0, N_nfft, dt_grid, T_window, K_B);
  mf_nfft_fftw_forward(N_nfft, K_B);
  mf_nfft_sample_kernel(kind, params, ctx, support, 1, N_nfft, dt_grid, T_window, K_g);
  mf_nfft_fftw_forward(N_nfft, K_g);
  mf_nfft_sample_kernel(kind, params, ctx, support, 2, N_nfft, dt_grid, T_window, K_gg);
  mf_nfft_fftw_forward(N_nfft, K_gg);

  /* The five forward NFFTs.  Each writes into a length-N output array
   * indexed by data position; since tau_j = t_j in Phase A, the j-th
   * entry is the value at tau = t[j]. */
  double *Sw   = (double *) malloc(N * sizeof(double));
  double *Wg   = (double *) malloc(N * sizeof(double));
  double *Wgg  = (double *) malloc(N * sizeof(double));
  double *Wy   = (double *) malloc(N * sizeof(double));
  double *Wyg  = (double *) malloc(N * sizeof(double));
  if (Sw == NULL || Wg == NULL || Wgg == NULL || Wy == NULL || Wyg == NULL) {
    if (Sw)  free(Sw);
    if (Wg)  free(Wg);
    if (Wgg) free(Wgg);
    if (Wy)  free(Wy);
    if (Wyg) free(Wyg);
    fftw_free(rho_w); fftw_free(rho_wy);
    fftw_free(K_B); fftw_free(K_g); fftw_free(K_gg);
    nfft_finalize(&plan);
    return -1;
  }

  mf_nfft_evaluate(&plan, N_nfft, rho_w,  K_B,  Sw,  N);
  mf_nfft_evaluate(&plan, N_nfft, rho_w,  K_g,  Wg,  N);
  mf_nfft_evaluate(&plan, N_nfft, rho_w,  K_gg, Wgg, N);
  mf_nfft_evaluate(&plan, N_nfft, rho_wy, K_B,  Wy,  N);
  mf_nfft_evaluate(&plan, N_nfft, rho_wy, K_g,  Wyg, N);

  /* Combine: a_hat = (Wyg - Wg*Wy/Sw) / (Wgg - Wg*Wg/Sw),
   *          SNR   = a_hat * sqrt(Wgg - Wg*Wg/Sw). */
  for (j = 0; j < Ntau; j++) {
    double Sw_j  = Sw[j],  Wg_j  = Wg[j],  Wgg_j = Wgg[j];
    double Wy_j  = Wy[j],  Wyg_j = Wyg[j];
    if (Sw_j > MF_NFFT_TINY) {
      double Sgg_c = Wgg_j - (Wg_j * Wg_j) / Sw_j;
      double Syg_c = Wyg_j - (Wg_j * Wy_j) / Sw_j;
      if (Sgg_c > MF_NFFT_TINY) {
        double a = Syg_c / Sgg_c;
        amp_out[j] = a;
        snr_out[j] = a * sqrt(Sgg_c);
        continue;
      }
    }
    amp_out[j] = 0.0;
    snr_out[j] = MF_NFFT_ERROR_SCORE;
  }

  free(Sw); free(Wg); free(Wgg); free(Wy); free(Wyg);
  fftw_free(rho_w); fftw_free(rho_wy);
  fftw_free(K_B); fftw_free(K_g); fftw_free(K_gg);
  nfft_finalize(&plan);
  ret = 0;
  return ret;
}

#else /* !HAVE_NFFT3 */

/* Stub when libnfft3 / libfftw3 were not detected at configure time.
 * Caller must ensure mode != MF_MODE_NFFT in this build. */
int mf_nfft_periodogram(int N, const double *t, const double *mag,
                        const double *w_in,
                        int kind, const double *params,
                        const _MFEvalCtx *ctx, double support,
                        int Ntau, const double *t_tau,
                        double *snr_out, double *amp_out)
{
  (void) N; (void) t; (void) mag; (void) w_in;
  (void) kind; (void) params; (void) ctx; (void) support;
  (void) Ntau; (void) t_tau;
  (void) snr_out; (void) amp_out;
  return -1;
}

#endif /* HAVE_NFFT3 */
