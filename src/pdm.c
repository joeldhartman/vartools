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
/*  Phase Dispersion Minimization (PDM) period search for vartools.
 *
 *  This file was developed with assistance from Claude (Anthropic) by
 *  Joel Hartman.
 *
 *  References:
 *    Stellingwerf, R. F. 1978, ApJ, 224, 953
 *        Original PDM; proposed F-distribution for theta (later shown incorrect).
 *    Schwarzenberg-Czerny, A. 1989, MNRAS, 241, 153
 *        Pointed out that Stellingwerf's F-distribution claim is biased.
 *    Schwarzenberg-Czerny, A. 1997, ApJ, 489, 941
 *        "The Correct Probability Distribution for the PDM Periodogram":
 *        derived theta ~ Beta((N-M)/2, (M-1)/2) under H0 via Cochran's theorem
 *        on SS_total = SS_within + SS_between.  This is the result used here.
 *    Zalian, C., Chadid, M., Stellingwerf, R. F. 2014, MNRAS, 440, 68
 *        Modern restatement of the SCz 1997 Beta-distribution FAP.
 *    cuvarbase (https://github.com/johnh2o2/cuvarbase, GPL-3)
 *        Source for the linear-interpolated bin variant ("linterp").  The PDM
 *        module in cuvarbase was implemented by Attila Bodi.
 *
 *  Initial scope (commit 1 of the staged PDM rollout):
 *    - variants:    step  (Stellingwerf 1978 fixed bins)
 *                   linterp  (cuvarbase-style linear interpolation between bin means)
 *    - statistic:   theta = SS_within / SS_total  in [0,1].  Modern parametrisation
 *                   (matches cuvarbase); drops Stellingwerf's original (N-1)/(N-M)
 *                   prefactor.  Differs from Stellingwerf-original theta by at most
 *                   one part in N for typical light-curve sizes.  Lower theta =
 *                   stronger periodic signal.
 *    - weighting:   default 1/sigma^2 weighting if sigma column is meaningful;
 *                   `noerr` keyword falls back to unweighted accumulation.
 *    - peaks:       Npeaks via coarse + fine-tune sweep, period-multiple double-check.
 *    - FAP:         Always computed via the Schwarzenberg-Czerny 1997 incomplete-
 *                   beta formula.  Rigorous for `step` under the ANOVA assumptions
 *                   (Gaussian residuals, independent observations, bins populated
 *                   enough for the SS_total = SS_within + SS_between decomposition).
 *                   For `linterp` the same formula is applied as an approximation:
 *                   linear interpolation breaks the orthogonal SS decomposition that
 *                   the Beta result rests on, so the analytic FAP is no longer the
 *                   exact null distribution.  A planned `bootstrap` keyword will
 *                   recalibrate m_eff empirically (the trials factor inside the
 *                   FAP), not the distribution itself -- the analytic Beta-based
 *                   FAP is what we always emit.
 *
 *  Deferred to subsequent commits:
 *    - multicover (PDM2 / Schwarzenberg-Czerny Nb x Nc covers)
 *    - binless variants (tophat, gaussian)
 *    - bootstrap FAP
 *    - whiten, clip, fixperiodSNR, maskpoints (port from aov.c when needed)
 */

#include "commands.h"
#include "programdata.h"
#include "functions.h"

/* --- Constants --- */
#define PDM_ERROR_SCORE 100000.0
#define PDM_MAX_DOUBLE_CHECK_MULTIPLE 19
#define PDM_MAX_PERIOD_DIFF_MULTIPLE 5
#define PDM_DEFAULT_NBIN 8
#define PDM_MINVAR 1.e-32
#define PDM_CTMIN 1   /* minimum points per bin to contribute a per-bin variance */

/* PDM_KIND_* and PDM_FAP_* live in commands.h so the parser can use them. */

/* --- Per-bin accumulator block --- */
typedef struct {
  int Nbin;                /* number of bins currently in use */
  int sizehists;           /* allocated capacity */
  double *binsumw;         /* sum w_i within bin */
  double *binsumwy;        /* sum w_i y_i within bin */
  double *binsumwy2;       /* sum w_i y_i^2 within bin */
  unsigned long *binN;     /* unweighted point count per bin */
  double *binmean;         /* scratch: per-bin weighted mean (filled lazily) */
} _PDMHistType;

/* ---- Bin-array allocator / housekeeping ---- */

static void pdm_hist_alloc(_PDMHistType *h, int Nbin)
{
  if (Nbin < 2) Nbin = 2;
  if (h->sizehists >= Nbin) {
    h->Nbin = Nbin;
    return;
  }
  if (h->sizehists > 0) {
    free(h->binsumw);
    free(h->binsumwy);
    free(h->binsumwy2);
    free(h->binN);
    free(h->binmean);
  }
  h->binsumw   = (double *)        malloc(Nbin * sizeof(double));
  h->binsumwy  = (double *)        malloc(Nbin * sizeof(double));
  h->binsumwy2 = (double *)        malloc(Nbin * sizeof(double));
  h->binN      = (unsigned long *) malloc(Nbin * sizeof(unsigned long));
  h->binmean   = (double *)        malloc(Nbin * sizeof(double));
  if (h->binsumw == NULL || h->binsumwy == NULL || h->binsumwy2 == NULL ||
      h->binN == NULL || h->binmean == NULL)
    vt_error(ERR_MEMALLOC);
  h->sizehists = Nbin;
  h->Nbin = Nbin;
}

static void pdm_hist_free(_PDMHistType *h)
{
  if (h->sizehists > 0) {
    free(h->binsumw);
    free(h->binsumwy);
    free(h->binsumwy2);
    free(h->binN);
    free(h->binmean);
  }
  h->binsumw = h->binsumwy = h->binsumwy2 = NULL;
  h->binmean = NULL;
  h->binN = NULL;
  h->sizehists = 0;
  h->Nbin = 0;
}

static void pdm_hist_reset(_PDMHistType *h)
{
  memset(h->binsumw,   0, h->Nbin * sizeof(double));
  memset(h->binsumwy,  0, h->Nbin * sizeof(double));
  memset(h->binsumwy2, 0, h->Nbin * sizeof(double));
  memset(h->binN,      0, h->Nbin * sizeof(unsigned long));
}

/* ---- Period-similarity helpers (mirror aov.c) ---- */

static int pdm_isDifferentPeriods(double p1, double p2, double T)
{
  int a, b;
  double p1mul;
  if (T * fabs(p2 - p1) < p1 * p2)
    return 0;
  for (a = 1; a < PDM_MAX_PERIOD_DIFF_MULTIPLE; a++)
    for (b = a + 1; b <= PDM_MAX_PERIOD_DIFF_MULTIPLE; b++) {
      p1mul = p1 * b / a;
      if (T * fabs(p2 - p1mul) < p2 * p1mul)
        return 0;
    }
  return 1;
}

/* ---- Whitening helper ---- */

/* Subtract the (weighted) bin-mean phase model at `period` from mag in-place,
 * using `Nbin` step bins.  Used between whiten cycles in findPeaks_pdm.
 *
 * For non-step variants this is a deliberate simplification: -aov_harm follows
 * the same convention (always whiten with bin means, regardless of the
 * statistic used for periodogram peak detection). */
static void whitenlc_pdm(int N, double *t, double *mag, double *w,
                         double period, int Nbin)
{
  int i, ibin;
  double X, dph;
  double *binsumw, *binsumwy, *binmean;
  unsigned long *binN;

  if (Nbin < 2) Nbin = 2;
  binsumw  = (double *)        calloc(Nbin, sizeof(double));
  binsumwy = (double *)        calloc(Nbin, sizeof(double));
  binmean  = (double *)        calloc(Nbin, sizeof(double));
  binN     = (unsigned long *) calloc(Nbin, sizeof(unsigned long));
  if (binsumw == NULL || binsumwy == NULL || binmean == NULL || binN == NULL)
    vt_error(ERR_MEMALLOC);

  for (i = 0; i < N; i++) {
    if (isnan(mag[i])) continue;
    if (w[i] <= 0.0) continue;
    dph = t[i] / period;
    X = dph - floor(dph);
    ibin = (int)(Nbin * X);
    if (ibin >= Nbin) ibin = Nbin - 1;
    if (ibin < 0) ibin = 0;
    binsumw[ibin]  += w[i];
    binsumwy[ibin] += w[i] * mag[i];
    binN[ibin]++;
  }
  for (i = 0; i < Nbin; i++)
    binmean[i] = (binN[i] >= PDM_CTMIN && binsumw[i] > 0.0)
                 ? binsumwy[i] / binsumw[i] : 0.0;

  for (i = 0; i < N; i++) {
    if (isnan(mag[i])) continue;
    dph = t[i] / period;
    X = dph - floor(dph);
    ibin = (int)(Nbin * X);
    if (ibin >= Nbin) ibin = Nbin - 1;
    if (ibin < 0) ibin = 0;
    mag[i] -= binmean[ibin];
  }

  free(binsumw); free(binsumwy); free(binmean); free(binN);
}

/* ---- Core theta evaluator ---- */

/* Returns theta at a single test period.
 *
 * For step/linterp:  Nc=1, Nb=h->Nbin = number of phase bins.
 * For multicover:    h->Nbin = Nb*Nc fine bins; theta is the SCz/S78 pooled
 *                    SS_within (= average per-cover theta) over Nc covers.
 * For tophat/gauss:  no histogram; per-point model from a phase-window or
 *                    Gaussian kernel of half-width / sigma `dphi`.  O(N^2)
 *                    per trial period -- the per-point neighbour search is
 *                    inherent to the binless model.
 *
 * Single-pass O(N) binning into the (Nb*Nc)-element fine histogram; multicover
 * then combines Nc consecutive fine bins per wide bin per cover, so the total
 * cost is O(N + Nc^2*Nb) per trial period rather than O(Nc*N).
 *
 * Returns PDM_ERROR_SCORE on degenerate input. */
static double pdm_theta(int kind, int N, double *t, double *mag, double *w,
                        double mu_total, double var_total, double period,
                        int Nb, int Nc, double dphi, _PDMHistType *h)
{
  int i, ibin;
  double X, Y, wi, dph, alpha, mu_interp;
  double s_sq;
  int Nfilled;
  int total_bins = h->Nbin;             /* = Nb*Nc */

  if (period <= 0.0 || var_total < PDM_MINVAR)
    return PDM_ERROR_SCORE;

  /* Binless variants do not use the histogram.  For each point j the model is
   * a local weighted mean of phase-neighbours; theta = sum w_j (y_j - mu_j)^2
   * / var_total.  Phase distances are wrapped to (-0.5, 0.5]. */
  if (kind == PDM_KIND_TOPHAT || kind == PDM_KIND_GAUSS) {
    double s_sq_binless = 0.0;
    if (dphi <= 0.0) return PDM_ERROR_SCORE;
    if (kind == PDM_KIND_GAUSS) {
      double two_sigma2 = 2.0 * dphi * dphi;
      for (i = 0; i < N; i++) {
        double yi, wi_, num = 0.0, denom = 0.0, dph_jk, dph_i;
        int k;
        if (isnan(mag[i])) continue;
        wi_ = w[i];
        if (wi_ <= 0.0) continue;
        dph_i = t[i] / period;
        for (k = 0; k < N; k++) {
          double wk;
          if (isnan(mag[k])) continue;
          wk = w[k];
          if (wk <= 0.0) continue;
          dph_jk = (t[k] - t[i]) / period;
          dph_jk -= floor(dph_jk + 0.5);     /* wrap to (-0.5, 0.5] */
          double g = wk * exp(-0.5 * dph_jk * dph_jk / (dphi * dphi));
          num   += g * mag[k];
          denom += g;
        }
        if (denom <= 0.0) continue;
        yi = mag[i] - num / denom;
        s_sq_binless += wi_ * yi * yi;
      }
    } else {
      /* tophat: hard phase window of half-width dphi */
      for (i = 0; i < N; i++) {
        double yi, wi_, num = 0.0, denom = 0.0, dph_jk;
        int k;
        if (isnan(mag[i])) continue;
        wi_ = w[i];
        if (wi_ <= 0.0) continue;
        for (k = 0; k < N; k++) {
          double wk;
          if (isnan(mag[k])) continue;
          wk = w[k];
          if (wk <= 0.0) continue;
          dph_jk = (t[k] - t[i]) / period;
          dph_jk -= floor(dph_jk + 0.5);
          if (fabs(dph_jk) > dphi) continue;
          num   += wk * mag[k];
          denom += wk;
        }
        if (denom <= 0.0) continue;
        yi = mag[i] - num / denom;
        s_sq_binless += wi_ * yi * yi;
      }
    }
    if (s_sq_binless < 0.0) s_sq_binless = 0.0;
    return s_sq_binless / var_total;
  }

  pdm_hist_reset(h);

  for (i = 0; i < N; i++) {
    if (isnan(mag[i])) continue;
    wi = w[i];
    if (wi <= 0.0) continue;
    dph = t[i] / period;
    X = dph - floor(dph);                 /* phase in [0,1) */
    ibin = (int)(total_bins * X);
    if (ibin >= total_bins) ibin = total_bins - 1;
    if (ibin < 0) ibin = 0;
    Y = mag[i];
    h->binsumw[ibin]   += wi;
    h->binsumwy[ibin]  += wi * Y;
    h->binsumwy2[ibin] += wi * Y * Y;
    h->binN[ibin]++;
  }

  /* For step/linterp need at least two occupied bins; for multicover apply the
   * same check on the fine bins (looser than a per-cover check, but cheap). */
  Nfilled = 0;
  for (i = 0; i < total_bins; i++)
    if (h->binN[i] >= PDM_CTMIN) Nfilled++;
  if (Nfilled < 2)
    return PDM_ERROR_SCORE;

  if (kind == PDM_KIND_STEP) {
    /* Stellingwerf 1978: theta = (sum_j sum_{i in j} w_i (y_i - mu_j)^2) / var_total
     * Using single-pass identity:
     *     sum_{i in j} w_i (y_i - mu_j)^2 = binsumwy2[j] - binsumwy[j]^2 / binsumw[j] */
    s_sq = 0.0;
    for (i = 0; i < total_bins; i++) {
      if (h->binN[i] < PDM_CTMIN) continue;
      s_sq += h->binsumwy2[i] - (h->binsumwy[i] * h->binsumwy[i] / h->binsumw[i]);
    }
    if (s_sq < 0.0) s_sq = 0.0;            /* numerical guard */
    return s_sq / var_total;
  }
  else if (kind == PDM_KIND_LINTERP) {
    /* cuvarbase-style: bin means at bin centres; per-point model is linear
     * interpolation between adjacent bin centres.  Bin j has centre phase
     * (j + 0.5)/Nbin.  For a point at phase X with xn = X*Nbin and ibin =
     * floor(xn), alpha = xn - ibin - 0.5 is the signed distance from the
     * centre of ibin, in (-0.5, 0.5).  Sign of alpha tells us which adjacent
     * bin to interpolate against; the relative distance from the "lo" anchor
     * (the bin whose centre is to the left of the point) gives the weight
     * on the "hi" anchor. */
    int ibin_lo, ibin_hi;
    double w_hi;
    for (i = 0; i < total_bins; i++) {
      if (h->binN[i] >= PDM_CTMIN)
        h->binmean[i] = h->binsumwy[i] / h->binsumw[i];
      else
        h->binmean[i] = mu_total;          /* empty bin defaults to LC mean */
    }
    s_sq = 0.0;
    for (i = 0; i < N; i++) {
      if (isnan(mag[i])) continue;
      wi = w[i];
      if (wi <= 0.0) continue;
      dph = t[i] / period;
      X = dph - floor(dph);
      double xn = X * total_bins;
      ibin = (int)xn;
      if (ibin >= total_bins) ibin = total_bins - 1;
      if (ibin < 0) ibin = 0;
      alpha = xn - (double)ibin - 0.5;
      if (alpha < 0.0) {
        ibin_lo = ibin - 1; if (ibin_lo < 0)         ibin_lo += total_bins;
        ibin_hi = ibin;
        w_hi = alpha + 1.0;                /* xn - (lo_centre = ibin - 0.5) */
      } else {
        ibin_lo = ibin;
        ibin_hi = ibin + 1; if (ibin_hi >= total_bins)  ibin_hi -= total_bins;
        w_hi = alpha;                      /* xn - (lo_centre = ibin + 0.5) */
      }
      mu_interp = (1.0 - w_hi) * h->binmean[ibin_lo] + w_hi * h->binmean[ibin_hi];
      Y = mag[i] - mu_interp;
      s_sq += wi * Y * Y;
    }
    if (s_sq < 0.0) s_sq = 0.0;
    return s_sq / var_total;
  }
  else if (kind == PDM_KIND_MULTICOVER) {
    /* S78 multicover with the fine-binning optimisation.  Wide bin j of
     * cover c covers phase [(c + j*Nc)/(Nb*Nc), (c + j*Nc + Nc)/(Nb*Nc)),
     * i.e. Nc consecutive fine bins starting at (c + j*Nc) mod (Nb*Nc).
     *
     * SCz 1997 (pp. 943-944) explicitly notes that the analytic FAP is
     * unknown for Nc > 1 because the per-cover Theta values are correlated.
     * We still report theta_multicover = (sum_c SS_within^c) / (Nc*var_total),
     * which equals the simple per-cover average of theta_c.  The caller emits
     * an analytic FAP using the single-cover Beta((N-Nb)/2, (Nb-1)/2), which
     * is conservative; the help text documents this caveat. */
    int c, j, k, fb;
    double wsum, wysum, wy2sum;
    unsigned long Nsum;
    s_sq = 0.0;
    for (c = 0; c < Nc; c++) {
      for (j = 0; j < Nb; j++) {
        wsum = 0.0; wysum = 0.0; wy2sum = 0.0; Nsum = 0;
        for (k = 0; k < Nc; k++) {
          fb = (c + j * Nc + k) % total_bins;
          wsum   += h->binsumw[fb];
          wysum  += h->binsumwy[fb];
          wy2sum += h->binsumwy2[fb];
          Nsum   += h->binN[fb];
        }
        if (Nsum < PDM_CTMIN || wsum <= 0.0) continue;
        s_sq += wy2sum - wysum * wysum / wsum;
      }
    }
    if (s_sq < 0.0) s_sq = 0.0;
    return s_sq / ((double)Nc * var_total);
  }

  /* Unknown kind */
  return PDM_ERROR_SCORE;
}

/* Compute the periodogram (theta vs period) over a frequency sweep. */
static void pdm_periodogram(int kind, int N, double *t, double *mag, double *w,
                            double mu_total, double var_total,
                            int Nb, int Nc, double dphi,
                            int Nperiod, double *periods, double *periodogram,
                            _PDMHistType *h)
{
  int k;
  for (k = 0; k < Nperiod; k++)
    periodogram[k] = pdm_theta(kind, N, t, mag, w, mu_total, var_total,
                               periods[k], Nb, Nc, dphi, h);
}

/* ---- Bootstrap FAP calibration (mirrors -LS's bounded-statistic branch) ----
 *
 * Generates Nboot shuffled-LC realisations (sampling with replacement, like
 * -LS), computes a periodogram for each, records the minimum theta.  The
 * sorted log10(theta_min) values form an empirical reference distribution;
 * we fit a degree-1 polynomial to log10(FAP) vs log10(theta) in the most-
 * extreme 10% tail to extrapolate the FAP for peaks beyond the bootstrap
 * range.  This is the analog of -LS's polynomial extrapolation, with the
 * log-log transformation appropriate for a Beta-distributed statistic
 * (P(theta <= t) ~ t^((N-Nb)/2) for small theta under H_0). */
static void pdm_bootstrap_calibrate(
    int kind, int N, double *t, double *mag, double *w,
    int Nb, int Nc, double dphi,
    int Nperiod, double *periods, _PDMHistType *h,
    int Nboot,
    double *bootstrapdist /* [Nboot] -- on exit holds sorted log10(theta_min) */,
    double *fitcoeffs    /* [2]     -- log10_FAP = c0 + c1 * log10_theta tail fit */)
{
  int trial, i, j, k;
  double *mag_shuf, *w_shuf, *pgram, *bprobs;
  double sumw, mu_s, var_s;

  mag_shuf = (double *) malloc(N * sizeof(double));
  w_shuf   = (double *) malloc(N * sizeof(double));
  pgram    = (double *) malloc(Nperiod * sizeof(double));
  bprobs   = (double *) malloc(Nboot * sizeof(double));
  if (mag_shuf == NULL || w_shuf == NULL || pgram == NULL || bprobs == NULL)
    vt_error(ERR_MEMALLOC);

  for (trial = 0; trial < Nboot; trial++) {
    /* Shuffle with replacement (matches -LS's randlong-based sampler). */
    for (j = 0; j < N; j++) {
      long klong = randlong(N - 1);
      mag_shuf[j] = mag[klong];
      w_shuf[j]   = w[klong];
    }
    /* Renormalise weights and recompute the shuffled mu / var. */
    sumw = 0.0;
    for (j = 0; j < N; j++) sumw += w_shuf[j];
    if (!(sumw > 0.0)) { bootstrapdist[trial] = 1.0; continue; }
    for (j = 0; j < N; j++) w_shuf[j] /= sumw;
    mu_s = 0.0;
    for (j = 0; j < N; j++) mu_s += w_shuf[j] * mag_shuf[j];
    var_s = 0.0;
    for (j = 0; j < N; j++) {
      double dy = mag_shuf[j] - mu_s;
      var_s += w_shuf[j] * dy * dy;
    }
    if (var_s < PDM_MINVAR) { bootstrapdist[trial] = 1.0; continue; }

    pdm_periodogram(kind, N, t, mag_shuf, w_shuf, mu_s, var_s,
                    Nb, Nc, dphi, Nperiod, periods, pgram, h);

    /* Record the most extreme (smallest) theta over the search range. */
    double tmin = 1.0;
    for (i = 0; i < Nperiod; i++) {
      double v = pgram[i];
      if (v < PDM_ERROR_SCORE && v * 0.0 == 0.0)
        if (v < tmin) tmin = v;
    }
    bootstrapdist[trial] = tmin;
  }

  /* Sort ascending (small theta = front = most extreme). */
  mysort1(Nboot, bootstrapdist);

  /* log-transform; floor at PDM_MINVAR to keep log finite. */
  for (k = 0; k < Nboot; k++) {
    if (bootstrapdist[k] < PDM_MINVAR) bootstrapdist[k] = PDM_MINVAR;
    bootstrapdist[k] = log10(bootstrapdist[k]);
  }

  /* Empirical log10(FAP) at the k-th sorted entry = log10((k+1)/Nboot). */
  for (k = 0; k < Nboot; k++)
    bprobs[k] = log10((double)(k + 1) / (double)Nboot);

  /* Fit a degree-1 polynomial to the most-extreme 10% tail (smallest theta
   * = front of the sorted array).  log10_FAP = c0 + c1 * log10_theta. */
  {
    int nfit = (int) floor(0.1 * (double)Nboot);
    if (nfit < 2) nfit = 2;
    if (nfit > Nboot) nfit = Nboot;
    fitpoly(nfit, bootstrapdist, bprobs, NULL, 1, 0, fitcoeffs, NULL);
  }

  free(mag_shuf); free(w_shuf); free(pgram); free(bprobs);
}

/* Look up empirical log10(FAP) for theta_peak.  bootstrapdist is sorted
 * ascending log10(theta) values.  Returns log10(FAP) for that theta;
 * the caller converts to NEG_LN_FAP. */
static double pdm_bootstrap_log10fap(double theta_peak, int Nboot,
                                     double *bootstrapdist_log,
                                     double *fitcoeffs)
{
  double lt;
  int ind;
  if (theta_peak < PDM_MINVAR) theta_peak = PDM_MINVAR;
  if (theta_peak > 1.0) theta_peak = 1.0;
  lt = log10(theta_peak);
  ind = findX(bootstrapdist_log, lt, 0, Nboot);
  if (ind == 0) {
    /* theta more extreme than any bootstrap trial -- extrapolate via tail fit */
    return fitcoeffs[0] + fitcoeffs[1] * lt;
  }
  return log10((double) ind / (double) Nboot);
}

/* ---- Find Npeaks lowest-theta periods, with fine-tune + multiple check ---- */

void findPeaks_pdm(double *t_, double *mag_, double *sig_, int N,
                   int kind, int Nbin, int Nc, double dphi, int useerr,
                   double minP, double maxP, double subsample, double finetune,
                   int Npeaks,
                   double *perpeaks, double *thetapeaks, double *peakSNR,
                   double *peakFAP,
                   double clip, int clipiter,
                   double *ave_theta, double *rms_theta,
                   double *ave_theta_whiten, double *rms_theta_whiten,
                   int whiten, int bootstrap_Nboot,
                   int outflag, char *outname, int ascii,
                   int fixperiodSNR_on, double fixperiodSNR_period,
                   double *fixperiodSNR_value, double *fixperiodSNR_SNR,
                   double *fixperiodSNR_FAP,
                   int usemask, _Variable *maskvar,
                   int lcnum, int lc_name_num)
{
  int i, j, k, foundsofar, test, Nperiod, a, b, abest, bbest, ismultiple, m_eff;
  int Nused;
  double *t = NULL, *mag = NULL, *w = NULL;
  double *periods = NULL, *periodogram = NULL;
  double T, mu_total, var_total;
  double sumw, wnorm;
  double freq, minfreq, freqstep, smallfreqstep;
  double testperiod, bestscore, lastpoint, score;
  double ave_per, std_per, negln_m_eff;
  double a_, b_;
  long double Sum, Sumsqr;
  FILE *outfile = NULL;
  _PDMHistType h = {0};
  /* Bootstrap state (allocated only if bootstrap_Nboot > 0). */
  double *bootstrapdist = NULL;
  double bootstrap_fitcoeffs[2] = {0.0, 0.0};
  /* Whiten state (allocated only if whiten=1). */
  double **periodogram_whiten = NULL;

  /* Empty / trivial light curve -> fill zeros and bail. */
  if (N <= 1) {
    for (i = 0; i < Npeaks; i++) {
      perpeaks[i]   = -1.0;
      thetapeaks[i] = 1.0;
      peakSNR[i]    = 0.0;
      peakFAP[i]    = 0.0;
    }
    *ave_theta = 0.0;
    *rms_theta = 0.0;
    if (fixperiodSNR_on) {
      *fixperiodSNR_value = 1.0;
      *fixperiodSNR_SNR   = 0.0;
      *fixperiodSNR_FAP   = 0.0;
    }
    return;
  }

  /* Use the configured Nbin (or default if <=0) */
  if (Nbin <= 0) Nbin = PDM_DEFAULT_NBIN;
  if (Nc <= 0) Nc = 1;
  if (kind != PDM_KIND_MULTICOVER) Nc = 1;   /* covers are a multicover concept */
  pdm_hist_alloc(&h, Nbin * Nc);

  /* Copy t, mag, and build weight vector. */
  t   = (double *)malloc(N * sizeof(double));
  mag = (double *)malloc(N * sizeof(double));
  w   = (double *)malloc(N * sizeof(double));
  if (t == NULL || mag == NULL || w == NULL) vt_error(ERR_MEMALLOC);

  Nused = 0;
  for (i = 0; i < N; i++) {
    if (isnan(mag_[i])) continue;
    if (usemask && maskvar != NULL) {
      if (!(EvaluateVariable_Double(lc_name_num, lcnum, i, maskvar) > VARTOOLS_MASK_TINY))
        continue;
    }
    if (useerr) {
      if (!(sig_[i] > 0.0) || isnan(sig_[i])) continue;
      w[Nused] = 1.0 / (sig_[i] * sig_[i]);
    } else {
      w[Nused] = 1.0;
    }
    t[Nused]   = t_[i];
    mag[Nused] = mag_[i];
    Nused++;
  }
  if (Nused <= 1) {
    for (i = 0; i < Npeaks; i++) {
      perpeaks[i]   = -1.0;
      thetapeaks[i] = 1.0;
      peakSNR[i]    = 0.0;
      peakFAP[i]    = 0.0;
    }
    *ave_theta = 0.0;
    *rms_theta = 0.0;
    if (fixperiodSNR_on) {
      *fixperiodSNR_value = 1.0;
      *fixperiodSNR_SNR   = 0.0;
      *fixperiodSNR_FAP   = 0.0;
    }
    free(t); free(mag); free(w);
    pdm_hist_free(&h);
    return;
  }

  /* Subtract the time mean from the local t[] copy.  PDM theta at any fixed
   * trial period is invariant to constant time shifts modulo bin alignment,
   * but the SUB-RAYLEIGH structure of theta(f) across oversampled trial
   * frequencies is only correlated when |t|/T is small (concretely:
   * adjacent-trial phase shift dPhase = subsample * |t|/T must be << 1/Nbin).
   * Without this subtraction, large epochs (e.g. full JD ~ 2.46e6) decorrelate
   * adjacent oversampled trial frequencies and the periodogram looks noisy
   * around the true peak (Bodi 2026 priv. comm., comparing against cuvarbase's
   * mean-subtracted convention).  Subtraction is done ONCE here; the caller's
   * input array t_[] is unchanged. */
  {
    double mean_t = 0.0;
    for (i = 0; i < Nused; i++) mean_t += t[i];
    mean_t /= (double) Nused;
    for (i = 0; i < Nused; i++) t[i] -= mean_t;
  }

  /* Normalise weights to sum to 1 so theta is dimensionless. */
  sumw = 0.0;
  for (i = 0; i < Nused; i++) sumw += w[i];
  if (sumw <= 0.0) sumw = 1.0;
  wnorm = 1.0 / sumw;
  for (i = 0; i < Nused; i++) w[i] *= wnorm;

  /* Global weighted mean + variance (the denominator for theta). */
  mu_total = 0.0;
  for (i = 0; i < Nused; i++) mu_total += w[i] * mag[i];
  var_total = 0.0;
  for (i = 0; i < Nused; i++) {
    double dy = mag[i] - mu_total;
    var_total += w[i] * dy * dy;
  }
  if (var_total < PDM_MINVAR) {
    /* All points coincident -- nothing to find. */
    for (i = 0; i < Npeaks; i++) {
      perpeaks[i]   = -1.0;
      thetapeaks[i] = 1.0;
      peakSNR[i]    = 0.0;
      peakFAP[i]    = 0.0;
    }
    *ave_theta = 0.0;
    *rms_theta = 0.0;
    if (fixperiodSNR_on) {
      *fixperiodSNR_value = 1.0;
      *fixperiodSNR_SNR   = 0.0;
      *fixperiodSNR_FAP   = 0.0;
    }
    free(t); free(mag); free(w);
    pdm_hist_free(&h);
    return;
  }

  /* Frequency grid */
  T = t[Nused - 1] - t[0];
  if (T <= 0.0) T = 1.0;
  freqstep = subsample / T;

  Nperiod = 0;
  freq = 1.0 / minP;
  minfreq = 1.0 / maxP;
  while (freq >= minfreq) {
    freq -= freqstep;
    Nperiod++;
  }
  if (Nperiod <= 0) Nperiod = 1;

  periods     = (double *)malloc(Nperiod * sizeof(double));
  periodogram = (double *)malloc(Nperiod * sizeof(double));
  if (periods == NULL || periodogram == NULL) vt_error(ERR_MEMALLOC);

  freq = 1.0 / minP;
  for (i = 0; i < Nperiod && freq >= minfreq; i++) {
    periods[i] = 1.0 / freq;
    freq -= freqstep;
  }
  Nperiod = i;

  m_eff = (int)((1.0 / minP - 1.0 / maxP) * T);
  if (subsample > 1.0) m_eff = (int)((double)m_eff / subsample);
  if (m_eff < 1) m_eff = 1;
  negln_m_eff = -log((double)m_eff);

  /* Bootstrap calibration: run Nboot shuffled-LC trials BEFORE any whitening
   * (the H_0 distribution is determined by the original LC's time sampling
   * and noise structure).  Stays NULL when bootstrap is off. */
  if (bootstrap_Nboot > 0) {
    bootstrapdist = (double *) malloc(bootstrap_Nboot * sizeof(double));
    if (bootstrapdist == NULL) vt_error(ERR_MEMALLOC);
    pdm_bootstrap_calibrate(kind, Nused, t, mag, w, Nbin, Nc, dphi,
                            Nperiod, periods, &h, bootstrap_Nboot,
                            bootstrapdist, bootstrap_fitcoeffs);
  }

  a_ = 0.5 * ((double)(Nbin - 1));
  b_ = 0.5 * ((double)(Nused - Nbin));

  if (!whiten) {

  /* Compute periodogram */
  pdm_periodogram(kind, Nused, t, mag, w, mu_total, var_total,
                  Nbin, Nc, dphi, Nperiod, periods, periodogram, &h);

  /* Mean / RMS of the periodogram (for SNR).  Skip ERROR scores. */
  Sum = 0.0; Sumsqr = 0.0;
  int Ngood = 0;
  for (i = 0; i < Nperiod; i++) {
    double v = periodogram[i];
    if (v < PDM_ERROR_SCORE && v * 0.0 == 0.0) {
      Sum += (long double) v;
      Sumsqr += (long double)(v * v);
      Ngood++;
    }
  }
  if (Ngood > 0) {
    Sum /= Ngood;
    Sumsqr /= Ngood;
    ave_per = (double) Sum;
    std_per = sqrt((double)(Sumsqr - Sum*Sum));
    if (!(std_per > 0.0)) std_per = 1.0;
  } else {
    ave_per = 0.5;
    std_per = 1.0;
  }

  /* Iterative sigma-clipping of the noise estimate.  Without this the deep
   * signal peaks (at P, P/2, 2P, ...) contaminate ave_per / std_per and
   * depress the SNR.  For PDM the peaks lie below the noise floor, so we
   * keep periodogram values v > ave_per - clip * std_per.  Affects only
   * the SNR noise estimate -- peak finding still sees the full periodogram. */
  if (clip > 0.0 && Ngood > 0) {
    int nclippedthis = 0, nclippedlast;
    do {
      nclippedlast = nclippedthis;
      nclippedthis = 0;
      Sum = 0.0; Sumsqr = 0.0;
      int Ngood2 = 0;
      for (i = 0; i < Nperiod; i++) {
        double v = periodogram[i];
        if (v < PDM_ERROR_SCORE && v * 0.0 == 0.0) {
          if (v > ave_per - clip * std_per) {
            Sum += (long double) v;
            Sumsqr += (long double)(v * v);
            Ngood2++;
          } else {
            nclippedthis++;
          }
        }
      }
      if (Ngood2 > 0) {
        Sum /= Ngood2;
        Sumsqr /= Ngood2;
        ave_per = (double) Sum;
        std_per = sqrt((double)(Sumsqr - Sum*Sum));
        if (!(std_per > 0.0)) std_per = 1.0;
      }
    } while (clipiter && nclippedthis > nclippedlast);
  }

  *ave_theta = ave_per;
  *rms_theta = std_per;

  /* Optional periodogram dump (ASCII). */
  if (outflag) {
    if ((outfile = fopen(outname, "w")) == NULL) {
      fprintf(stderr, "Cannot write periodogram to %s\n", outname);
      exit(3);
    }
    if (ascii) {
      fprintf(outfile, "#Period PDM_Theta\n");
      for (i = 0; i < Nperiod; i++) {
        if (periodogram[i] < PDM_ERROR_SCORE && periodogram[i] * 0.0 == 0.0)
          fprintf(outfile, "%f %f\n", periods[i], periodogram[i]);
      }
    } else {
      fwrite(&Nperiod, 4, 1, outfile);
      fwrite(&ave_per, 8, 1, outfile);
      fwrite(&std_per, 8, 1, outfile);
      fwrite(periods,  8, Nperiod, outfile);
      fwrite(periodogram, 8, Nperiod, outfile);
    }
    fclose(outfile);
  }

  /* Reduce to local minima (theta dips below neighbours) */
  {
    int kk;
    lastpoint = periodogram[0] - 1.0;
    int kept = 0;
    for (kk = 0; kk < Nperiod - 1; kk++) {
      if (periodogram[kk] < lastpoint && periodogram[kk] < periodogram[kk + 1]) {
        lastpoint = periodogram[kk];
        periodogram[kept] = periodogram[kk];
        periods[kept] = periods[kk];
        kept++;
      } else {
        lastpoint = periodogram[kk];
      }
    }
    if (periodogram[Nperiod - 1] < lastpoint) {
      periodogram[kept] = periodogram[Nperiod - 1];
      periods[kept] = periods[Nperiod - 1];
      kept++;
    }
    Nperiod = kept;
  }

  /* Pick Npeaks distinct lowest-theta minima */
  foundsofar = 0;
  i = 0;
  while (foundsofar < Npeaks && i < Nperiod) {
    if (periodogram[i] < PDM_ERROR_SCORE && periodogram[i] * 0.0 == 0.0) {
      test = 1;
      for (j = 0; j < foundsofar; j++) {
        if (!pdm_isDifferentPeriods(perpeaks[j], periods[i], T)) {
          if (periodogram[i] < thetapeaks[j]) {
            perpeaks[j]   = periods[i];
            thetapeaks[j] = periodogram[i];
          }
          test = 0;
          break;
        }
      }
      if (test) {
        perpeaks[foundsofar]   = periods[i];
        thetapeaks[foundsofar] = periodogram[i];
        foundsofar++;
      }
    }
    i++;
  }
  /* Pad unused slots */
  for (k = foundsofar; k < Npeaks; k++) {
    perpeaks[k]   = -1.0;
    thetapeaks[k] = PDM_ERROR_SCORE + 1.0;
  }

  mysort2(Npeaks, thetapeaks, perpeaks);

  double worst_kept = thetapeaks[Npeaks - 1];
  for (; i < Nperiod; i++) {
    double v = periodogram[i];
    if (v < PDM_ERROR_SCORE && v * 0.0 == 0.0) {
      if (v < worst_kept) {
        test = 1;
        for (j = 0; j < Npeaks; j++) {
          if (!pdm_isDifferentPeriods(periods[i], perpeaks[j], T)) {
            if (v < thetapeaks[j]) {
              thetapeaks[j] = v;
              perpeaks[j]   = periods[i];
              mysort2(Npeaks, thetapeaks, perpeaks);
              worst_kept = thetapeaks[Npeaks - 1];
            }
            test = 0;
            break;
          }
        }
        if (test) {
          perpeaks[Npeaks - 1]   = periods[i];
          thetapeaks[Npeaks - 1] = v;
          mysort2(Npeaks, thetapeaks, perpeaks);
          worst_kept = thetapeaks[Npeaks - 1];
        }
      }
    }
  }

  /* High-resolution fine-tune around each peak */
  smallfreqstep = finetune / T;
  for (j = 0; j < Npeaks; j++) {
    if (thetapeaks[j] >= PDM_ERROR_SCORE) continue;
    freq    = dmin((1.0 / perpeaks[j]) + freqstep, 1.0 / minP);
    minfreq = dmax((1.0 / perpeaks[j]) - freqstep, 1.0 / maxP);
    while (freq >= minfreq) {
      testperiod = 1.0 / freq;
      score = pdm_theta(kind, Nused, t, mag, w, mu_total, var_total,
                        testperiod, Nbin, Nc, dphi, &h);
      if (score < PDM_ERROR_SCORE && score * 0.0 == 0.0)
        if (score < thetapeaks[j]) {
          thetapeaks[j] = score;
          perpeaks[j]   = testperiod;
        }
      freq -= smallfreqstep;
    }
    /* Double-check period multiples */
    bestscore = thetapeaks[j];
    ismultiple = 0;
    abest = bbest = 1;
    for (a = 1; a <= PDM_MAX_DOUBLE_CHECK_MULTIPLE; a++)
      for (b = 1; b <= PDM_MAX_DOUBLE_CHECK_MULTIPLE; b++)
        if (a != b) {
          testperiod = perpeaks[j] * a / b;
          if (testperiod > minP && testperiod < maxP) {
            score = pdm_theta(kind, Nused, t, mag, w, mu_total, var_total,
                              testperiod, Nbin, Nc, dphi, &h);
            if (score < PDM_ERROR_SCORE && score * 0.0 == 0.0)
              if (score < bestscore) {
                ismultiple = 1; abest = a; bbest = b; bestscore = score;
              }
          }
        }
    if (ismultiple) {
      perpeaks[j]   = perpeaks[j] * abest / bbest;
      thetapeaks[j] = bestscore;
    }
  }

  mysort2(Npeaks, thetapeaks, perpeaks);

  /* SNR and (optional) analytic FAP per peak */
  a_ = 0.5 * ((double)(Nbin - 1));
  b_ = 0.5 * ((double)(Nused - Nbin));
  for (j = 0; j < Npeaks; j++) {
    if (thetapeaks[j] < PDM_ERROR_SCORE) {
      /* For PDM the peak is a *minimum* of theta, so the SNR sign convention is
       * (mean - peak)/rms : positive when the peak dips below the mean. */
      peakSNR[j] = (ave_per - thetapeaks[j]) / std_per;
      if (bootstrap_Nboot > 0) {
        /* Empirical FAP from the bootstrap distribution (log-log polynomial
         * tail extrapolation; mirrors -LS for bounded statistics).  Converts
         * the log10 result to NEG_LN_FAP = -ln(FAP) = -ln(10) * log10(FAP). */
        double l10 = pdm_bootstrap_log10fap(thetapeaks[j], bootstrap_Nboot,
                                            bootstrapdist, bootstrap_fitcoeffs);
        peakFAP[j] = -log(10.0) * l10;
      }
      else if (b_ > 0.0 && a_ > 0.0) {
        double x = 1.0 - thetapeaks[j];
        if (x < 0.0) x = 0.0;
        if (x > 1.0) x = 1.0;
        /* Under H0 (no signal): theta ~ Beta(b,a) with a=(M-1)/2, b=(N-M)/2
         * (Schwarzenberg-Czerny 1997, ApJ 489, 941, eq. for PDM*; Cochran's
         * theorem on SS_total = SS_within + SS_between).  Then:
         *   FAP_single = I(theta; b, a) = 1 - I(1-theta; a, b)   [Pearson identity]
         * vartools' log1minusbetai(a,b,x) returns log(1 - I(a,b,x)), which equals
         * log(FAP_single) when called with x = 1 - theta.  We report NEG_LN_FAP
         * after an effective-trials-factor correction (negln_m_eff = -log m_eff),
         * i.e. NEG_LN_FAP = -log(m_eff * FAP_single). */
        peakFAP[j] = negln_m_eff - log1minusbetai(a_, b_, x);
      } else {
        peakFAP[j] = 0.0;
      }
    } else {
      perpeaks[j]   = -1.0;
      thetapeaks[j] = 1.0;
      peakSNR[j]    = 0.0;
      peakFAP[j]    = 0.0;
    }
  }

  } else {
    /* ---- whiten branch: iterate Npeaks times, whitening the LC between
     * cycles using the step-bin mean model.  Each cycle gets its own
     * mean/rms (with clipping) for SNR, stored in the per-cycle arrays
     * ave_theta_whiten[peakiter] / rms_theta_whiten[peakiter].
     *
     * To keep fixperiodSNR (which runs after this branch) honest, we save
     * the original mag and restore it at the end -- the fixed-period
     * evaluation runs against the ORIGINAL LC, not the progressively
     * whitened one. */
    int peakiter;
    int Nperiod_w = Nperiod;
    double *mag_orig    = (double *) malloc(Nused * sizeof(double));
    double mu_total_orig  = mu_total;
    double var_total_orig = var_total;
    if (mag_orig == NULL) vt_error(ERR_MEMALLOC);
    memcpy(mag_orig, mag, Nused * sizeof(double));

    periodogram_whiten = (double **) malloc(Npeaks * sizeof(double *));
    if (periodogram_whiten == NULL) vt_error(ERR_MEMALLOC);
    for (i = 0; i < Npeaks; i++) {
      periodogram_whiten[i] = (double *) malloc(Nperiod * sizeof(double));
      if (periodogram_whiten[i] == NULL) vt_error(ERR_MEMALLOC);
    }

    /* Initial periodogram (cycle 0) on the unwhitened LC. */
    pdm_periodogram(kind, Nused, t, mag, w, mu_total, var_total,
                    Nbin, Nc, dphi, Nperiod_w, periods, periodogram_whiten[0], &h);

    for (peakiter = 0; peakiter < Npeaks; peakiter++) {
      double cyc_ave = 0.5, cyc_std = 1.0;
      int Ngood_cyc = 0;
      long double Sum_cyc = 0.0, Sumsqr_cyc = 0.0;

      /* mean/rms of this cycle's periodogram, with the same iterative
       * sigma-clipping used in the single-pass path. */
      for (i = 0; i < Nperiod_w; i++) {
        double v = periodogram_whiten[peakiter][i];
        if (v < PDM_ERROR_SCORE && v * 0.0 == 0.0) {
          Sum_cyc    += (long double) v;
          Sumsqr_cyc += (long double)(v * v);
          Ngood_cyc++;
        }
      }
      if (Ngood_cyc > 0) {
        Sum_cyc /= Ngood_cyc; Sumsqr_cyc /= Ngood_cyc;
        cyc_ave = (double) Sum_cyc;
        cyc_std = sqrt((double)(Sumsqr_cyc - Sum_cyc * Sum_cyc));
        if (!(cyc_std > 0.0)) cyc_std = 1.0;
      }
      if (clip > 0.0 && Ngood_cyc > 0) {
        int nthis = 0, nlast;
        do {
          nlast = nthis; nthis = 0;
          Sum_cyc = 0.0; Sumsqr_cyc = 0.0;
          int Ngood2 = 0;
          for (i = 0; i < Nperiod_w; i++) {
            double v = periodogram_whiten[peakiter][i];
            if (v < PDM_ERROR_SCORE && v * 0.0 == 0.0) {
              if (v > cyc_ave - clip * cyc_std) {
                Sum_cyc += (long double) v;
                Sumsqr_cyc += (long double)(v * v);
                Ngood2++;
              } else {
                nthis++;
              }
            }
          }
          if (Ngood2 > 0) {
            Sum_cyc /= Ngood2; Sumsqr_cyc /= Ngood2;
            cyc_ave = (double) Sum_cyc;
            cyc_std = sqrt((double)(Sumsqr_cyc - Sum_cyc * Sum_cyc));
            if (!(cyc_std > 0.0)) cyc_std = 1.0;
          }
        } while (clipiter && nthis > nlast);
      }
      if (ave_theta_whiten != NULL) ave_theta_whiten[peakiter] = cyc_ave;
      if (rms_theta_whiten != NULL) rms_theta_whiten[peakiter] = cyc_std;

      /* Find the best (min-theta) period in this cycle. */
      double best = PDM_ERROR_SCORE;
      double best_p = -1.0;
      for (i = 0; i < Nperiod_w; i++) {
        double v = periodogram_whiten[peakiter][i];
        if (v < PDM_ERROR_SCORE && v * 0.0 == 0.0)
          if (v < best) { best = v; best_p = periods[i]; }
      }
      perpeaks[peakiter]   = best_p;
      thetapeaks[peakiter] = best;

      /* Fine-tune scan around the peak. */
      if (best_p > 0.0) {
        smallfreqstep = finetune / T;
        freq    = dmin((1.0 / best_p) + freqstep, 1.0 / minP);
        minfreq = dmax((1.0 / best_p) - freqstep, 1.0 / maxP);
        while (freq >= minfreq) {
          testperiod = 1.0 / freq;
          score = pdm_theta(kind, Nused, t, mag, w, mu_total, var_total,
                            testperiod, Nbin, Nc, dphi, &h);
          if (score < PDM_ERROR_SCORE && score * 0.0 == 0.0)
            if (score < thetapeaks[peakiter]) {
              thetapeaks[peakiter] = score;
              perpeaks[peakiter]   = testperiod;
            }
          freq -= smallfreqstep;
        }
        /* Period-multiple check. */
        bestscore = thetapeaks[peakiter]; ismultiple = 0; abest = bbest = 1;
        for (a = 1; a <= PDM_MAX_DOUBLE_CHECK_MULTIPLE; a++)
          for (b = 1; b <= PDM_MAX_DOUBLE_CHECK_MULTIPLE; b++)
            if (a != b) {
              testperiod = perpeaks[peakiter] * a / b;
              if (testperiod > minP && testperiod < maxP) {
                score = pdm_theta(kind, Nused, t, mag, w, mu_total, var_total,
                                  testperiod, Nbin, Nc, dphi, &h);
                if (score < PDM_ERROR_SCORE && score * 0.0 == 0.0)
                  if (score < bestscore) {
                    ismultiple = 1; abest = a; bbest = b; bestscore = score;
                  }
              }
            }
        if (ismultiple) {
          perpeaks[peakiter]   = perpeaks[peakiter] * abest / bbest;
          thetapeaks[peakiter] = bestscore;
        }
      }

      /* SNR + FAP for this peak using THIS cycle's noise estimate. */
      if (thetapeaks[peakiter] < PDM_ERROR_SCORE) {
        peakSNR[peakiter] = (cyc_ave - thetapeaks[peakiter]) / cyc_std;
        if (bootstrap_Nboot > 0) {
          double l10 = pdm_bootstrap_log10fap(thetapeaks[peakiter], bootstrap_Nboot,
                                              bootstrapdist, bootstrap_fitcoeffs);
          peakFAP[peakiter] = -log(10.0) * l10;
        }
        else if (b_ > 0.0 && a_ > 0.0) {
          double x = 1.0 - thetapeaks[peakiter];
          if (x < 0.0) x = 0.0;
          if (x > 1.0) x = 1.0;
          peakFAP[peakiter] = negln_m_eff - log1minusbetai(a_, b_, x);
        } else {
          peakFAP[peakiter] = 0.0;
        }
      } else {
        perpeaks[peakiter]   = -1.0;
        thetapeaks[peakiter] = 1.0;
        peakSNR[peakiter]    = 0.0;
        peakFAP[peakiter]    = 0.0;
      }

      /* Whiten the LC at the peak period; recompute periodogram for next cycle. */
      if (peakiter < Npeaks - 1 && perpeaks[peakiter] > 0.0) {
        whitenlc_pdm(Nused, t, mag, w, perpeaks[peakiter], Nbin);
        /* Recompute global mu / var on the whitened LC (mean is now ~0). */
        mu_total = 0.0;
        for (i = 0; i < Nused; i++) mu_total += w[i] * mag[i];
        var_total = 0.0;
        for (i = 0; i < Nused; i++) {
          double dy = mag[i] - mu_total;
          var_total += w[i] * dy * dy;
        }
        if (var_total < PDM_MINVAR) {
          /* LC is now flat; remaining cycles can't find anything. */
          for (i = 0; i < Nperiod_w; i++)
            periodogram_whiten[peakiter + 1][i] = 1.0;
        } else {
          pdm_periodogram(kind, Nused, t, mag, w, mu_total, var_total,
                          Nbin, Nc, dphi, Nperiod_w, periods, periodogram_whiten[peakiter + 1], &h);
        }
      }
    }

    /* Optional multi-column periodogram dump. */
    if (outflag) {
      if ((outfile = fopen(outname, "w")) == NULL) {
        fprintf(stderr, "Cannot write periodogram to %s\n", outname);
        exit(3);
      }
      if (ascii) {
        fprintf(outfile, "#Period");
        for (j = 0; j < Npeaks; j++) fprintf(outfile, " PDM_WhitenCycle_%d", j);
        fprintf(outfile, "\n");
        for (i = 0; i < Nperiod_w; i++) {
          double v0 = periodogram_whiten[0][i];
          if (!(v0 < PDM_ERROR_SCORE && v0 * 0.0 == 0.0)) continue;
          fprintf(outfile, "%f", periods[i]);
          for (j = 0; j < Npeaks; j++)
            fprintf(outfile, " %f", periodogram_whiten[j][i]);
          fprintf(outfile, "\n");
        }
      } else {
        fwrite(&Nperiod_w, 4, 1, outfile);
        fwrite(periods, 8, Nperiod_w, outfile);
        for (j = 0; j < Npeaks; j++) {
          double av = (ave_theta_whiten != NULL ? ave_theta_whiten[j] : 0.0);
          double rm = (rms_theta_whiten != NULL ? rms_theta_whiten[j] : 1.0);
          fwrite(&av, 8, 1, outfile);
          fwrite(&rm, 8, 1, outfile);
          fwrite(periodogram_whiten[j], 8, Nperiod_w, outfile);
        }
      }
      fclose(outfile);
    }

    *ave_theta = (ave_theta_whiten != NULL ? ave_theta_whiten[0] : 0.5);
    *rms_theta = (rms_theta_whiten != NULL ? rms_theta_whiten[0] : 1.0);

    /* Restore the original LC for fixperiodSNR; seed ave_per/std_per from
     * cycle 0 so the SNR_PeriodFix uses the original-LC noise estimate. */
    memcpy(mag, mag_orig, Nused * sizeof(double));
    mu_total  = mu_total_orig;
    var_total = var_total_orig;
    ave_per   = (ave_theta_whiten != NULL ? ave_theta_whiten[0] : 0.5);
    std_per   = (rms_theta_whiten != NULL ? rms_theta_whiten[0] : 1.0);
    if (!(std_per > 0.0)) std_per = 1.0;
    free(mag_orig);
  }

  /* fixperiodSNR: evaluate theta/SNR/FAP at a caller-supplied period (resolved
   * upstream from -aov/-ls/-pdm/-Injectharm/fix/list/fixcolumn).  Uses the
   * same ave/std-of-periodogram noise estimate (after clipping) as the peak
   * SNR; FAP uses the same single-cover Beta((N-Nb)/2, (Nb-1)/2) formula. */
  if (fixperiodSNR_on) {
    double t_fix;
    if (fixperiodSNR_period <= 0.0) {
      *fixperiodSNR_value = 1.0;
      *fixperiodSNR_SNR   = 0.0;
      *fixperiodSNR_FAP   = 0.0;
    } else {
      t_fix = pdm_theta(kind, Nused, t, mag, w, mu_total, var_total,
                        fixperiodSNR_period, Nbin, Nc, dphi, &h);
      if (t_fix < PDM_ERROR_SCORE && t_fix * 0.0 == 0.0) {
        *fixperiodSNR_value = t_fix;
        *fixperiodSNR_SNR = (ave_per - t_fix) / std_per;
        if (bootstrap_Nboot > 0) {
          double l10 = pdm_bootstrap_log10fap(t_fix, bootstrap_Nboot,
                                              bootstrapdist, bootstrap_fitcoeffs);
          *fixperiodSNR_FAP = -log(10.0) * l10;
        }
        else if (b_ > 0.0 && a_ > 0.0) {
          double x = 1.0 - t_fix;
          if (x < 0.0) x = 0.0;
          if (x > 1.0) x = 1.0;
          *fixperiodSNR_FAP = negln_m_eff - log1minusbetai(a_, b_, x);
        } else {
          *fixperiodSNR_FAP = 0.0;
        }
      } else {
        *fixperiodSNR_value = 1.0;
        *fixperiodSNR_SNR   = 0.0;
        *fixperiodSNR_FAP   = 0.0;
      }
    }
  }

  free(t); free(mag); free(w);
  free(periods); free(periodogram);
  pdm_hist_free(&h);
  if (bootstrapdist != NULL) free(bootstrapdist);
  if (periodogram_whiten != NULL) {
    for (i = 0; i < Npeaks; i++) free(periodogram_whiten[i]);
    free(periodogram_whiten);
  }
}

/* ---- Entry point invoked from processcommand.c ---- */

void RunPDMCommand(ProgramData *p, Command *c, _PDM *Pdm, int lcnum, int lc_name_num, int thisindex)
{
  char outname[MAXLEN];
  int i1, i2;
  int fix_on = 0;
  double fix_period = 1.0;
  double *fix_value_ptr = NULL, *fix_SNR_ptr = NULL, *fix_FAP_ptr = NULL;

  if (Pdm->operiodogram) {
    i1 = 0; i2 = 0;
    while (p->lcnames[lc_name_num][i1] != '\0') {
      if (p->lcnames[lc_name_num][i1] == '/') i2 = i1 + 1;
      i1++;
    }
    sprintf(outname, "%s/%s%s", Pdm->outdir, &p->lcnames[lc_name_num][i2], Pdm->suffix);
  }

  /* Resolve var/expr/fixed parameter sources to per-LC scalar values. */
  if (Pdm->minp_source == VARTOOLS_SOURCE_EVALEXPRESSION)
    Pdm->minp_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Pdm->minp_expr);
  else if (Pdm->minp_source == VARTOOLS_SOURCE_EXISTINGVARIABLE)
    Pdm->minp_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Pdm->minp_var);
  else
    Pdm->minp_vals[lcnum] = Pdm->minp;

  if (Pdm->maxp_source == VARTOOLS_SOURCE_EVALEXPRESSION)
    Pdm->maxp_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Pdm->maxp_expr);
  else if (Pdm->maxp_source == VARTOOLS_SOURCE_EXISTINGVARIABLE)
    Pdm->maxp_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Pdm->maxp_var);
  else
    Pdm->maxp_vals[lcnum] = Pdm->maxp;

  if (Pdm->subsample_source == VARTOOLS_SOURCE_EVALEXPRESSION)
    Pdm->subsample_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Pdm->subsample_expr);
  else if (Pdm->subsample_source == VARTOOLS_SOURCE_EXISTINGVARIABLE)
    Pdm->subsample_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Pdm->subsample_var);
  else
    Pdm->subsample_vals[lcnum] = Pdm->subsample;

  if (Pdm->finetune_source == VARTOOLS_SOURCE_EVALEXPRESSION)
    Pdm->finetune_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Pdm->finetune_expr);
  else if (Pdm->finetune_source == VARTOOLS_SOURCE_EXISTINGVARIABLE)
    Pdm->finetune_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Pdm->finetune_var);
  else
    Pdm->finetune_vals[lcnum] = Pdm->finetune;

  if (Pdm->Nbin_source == VARTOOLS_SOURCE_EVALEXPRESSION)
    Pdm->Nbin_vals[lcnum] = (int) ceil(EvaluateExpression(lc_name_num, lcnum, 0, Pdm->Nbin_expr));
  else if (Pdm->Nbin_source == VARTOOLS_SOURCE_EXISTINGVARIABLE)
    Pdm->Nbin_vals[lcnum] = (int) ceil(EvaluateVariable_Double(lc_name_num, lcnum, 0, Pdm->Nbin_var));
  else
    Pdm->Nbin_vals[lcnum] = Pdm->Nbin;

  if (Pdm->Nc_source == VARTOOLS_SOURCE_EVALEXPRESSION)
    Pdm->Nc_vals[lcnum] = (int) ceil(EvaluateExpression(lc_name_num, lcnum, 0, Pdm->Nc_expr));
  else if (Pdm->Nc_source == VARTOOLS_SOURCE_EXISTINGVARIABLE)
    Pdm->Nc_vals[lcnum] = (int) ceil(EvaluateVariable_Double(lc_name_num, lcnum, 0, Pdm->Nc_var));
  else
    Pdm->Nc_vals[lcnum] = Pdm->Nc;

  if (Pdm->dphi_source == VARTOOLS_SOURCE_EVALEXPRESSION)
    Pdm->dphi_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Pdm->dphi_expr);
  else if (Pdm->dphi_source == VARTOOLS_SOURCE_EXISTINGVARIABLE)
    Pdm->dphi_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Pdm->dphi_var);
  else
    Pdm->dphi_vals[lcnum] = Pdm->dphi;

  /* Resolve the fixperiodSNR period source (mirrors aov.c::RunAOVCommand). */
  if (Pdm->fixperiodSNR) {
    fix_on = 1;
    switch (Pdm->fixperiodSNR_pertype) {
      case PERTYPE_AOV:
        i1 = Pdm->fixperiodSNR_lastaovindex;
        if (c[i1 - thisindex].cnum == CNUM_AOV)
          Pdm->fixperiodSNR_periods[lcnum][0] = c[i1 - thisindex].Aov->peakperiods[lcnum][0];
        else if (c[i1 - thisindex].cnum == CNUM_HARMAOV)
          Pdm->fixperiodSNR_periods[lcnum][0] = c[i1 - thisindex].AovHarm->peakperiods[lcnum][0];
        break;
      case PERTYPE_LS:
        i1 = Pdm->fixperiodSNR_lastaovindex;
        Pdm->fixperiodSNR_periods[lcnum][0] = c[i1 - thisindex].Ls->peakperiods[lcnum][0];
        break;
      case PERTYPE_PDM:
        i1 = Pdm->fixperiodSNR_lastaovindex;
        Pdm->fixperiodSNR_periods[lcnum][0] = c[i1 - thisindex].Pdm->peakperiods[lcnum][0];
        break;
      case PERTYPE_INJECTHARM:
        i1 = Pdm->fixperiodSNR_lastaovindex;
        Pdm->fixperiodSNR_periods[lcnum][0] = c[i1 - thisindex].Injectharm->periodinject[lcnum];
        break;
      case PERTYPE_FIX:
        Pdm->fixperiodSNR_periods[lcnum][0] = Pdm->fixperiodSNR_fixedperiod;
        break;
      case PERTYPE_FIXCOLUMN:
        getoutcolumnvalue(Pdm->fixperiodSNR_linkedcolumn, lcnum, lc_name_num,
                          VARTOOLS_TYPE_DOUBLE, &(Pdm->fixperiodSNR_periods[lcnum][0]));
        break;
      default:
        /* PERTYPE_SPECIFIED (list-column): pre-loaded by inputlist parsing */
        break;
    }
    fix_period    = Pdm->fixperiodSNR_periods[lcnum][0];
    fix_value_ptr = &(Pdm->fixperiodSNR_peakvalues[lcnum]);
    fix_SNR_ptr   = &(Pdm->fixperiodSNR_peakSNR[lcnum]);
    fix_FAP_ptr   = &(Pdm->fixperiodSNR_peakFAP[lcnum]);
  }

  findPeaks_pdm(p->t[lcnum], p->mag[lcnum], p->sig[lcnum], p->NJD[lcnum],
                Pdm->kind, Pdm->Nbin_vals[lcnum], Pdm->Nc_vals[lcnum],
                Pdm->dphi_vals[lcnum], Pdm->useerr,
                Pdm->minp_vals[lcnum], Pdm->maxp_vals[lcnum],
                Pdm->subsample_vals[lcnum], Pdm->finetune_vals[lcnum],
                Pdm->Npeaks,
                Pdm->peakperiods[lcnum], Pdm->peakvalues[lcnum],
                Pdm->peakSNR[lcnum], Pdm->peakFAP[lcnum],
                Pdm->clip, Pdm->clipiter,
                &Pdm->avetheta[lcnum], &Pdm->rmstheta[lcnum],
                (Pdm->whiten ? Pdm->avetheta_whiten[lcnum] : NULL),
                (Pdm->whiten ? Pdm->rmstheta_whiten[lcnum] : NULL),
                Pdm->whiten, Pdm->bootstrap_Nboot,
                Pdm->operiodogram, outname, p->ascii,
                fix_on, fix_period, fix_value_ptr, fix_SNR_ptr, fix_FAP_ptr,
                Pdm->usemask, Pdm->maskvar, lcnum, lc_name_num);
}
