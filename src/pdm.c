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

/* ---- Core theta evaluator ---- */

/* Returns theta = s^2 / var_total at a single test period.
 * weights w[] sum to sumw (typically 1).  var_total is the weighted total
 * variance computed once by the caller (Sum w_i (y_i - mu_total)^2).
 * Returns PDM_ERROR_SCORE on degenerate input. */
static double pdm_theta(int kind, int N, double *t, double *mag, double *w,
                        double mu_total, double var_total, double period,
                        _PDMHistType *h)
{
  int i, ibin;
  double X, Y, wi, dph, alpha, mu_interp;
  double s_sq;
  int Nfilled;

  if (period <= 0.0 || var_total < PDM_MINVAR)
    return PDM_ERROR_SCORE;

  pdm_hist_reset(h);

  for (i = 0; i < N; i++) {
    if (isnan(mag[i])) continue;
    wi = w[i];
    if (wi <= 0.0) continue;
    dph = t[i] / period;
    X = dph - floor(dph);                 /* phase in [0,1) */
    ibin = (int)(h->Nbin * X);
    if (ibin >= h->Nbin) ibin = h->Nbin - 1;
    if (ibin < 0) ibin = 0;
    Y = mag[i];
    h->binsumw[ibin]   += wi;
    h->binsumwy[ibin]  += wi * Y;
    h->binsumwy2[ibin] += wi * Y * Y;
    h->binN[ibin]++;
  }

  /* Need at least two occupied bins, else the statistic is degenerate. */
  Nfilled = 0;
  for (i = 0; i < h->Nbin; i++)
    if (h->binN[i] >= PDM_CTMIN) Nfilled++;
  if (Nfilled < 2)
    return PDM_ERROR_SCORE;

  if (kind == PDM_KIND_STEP) {
    /* Stellingwerf 1978: theta = (sum_j sum_{i in j} w_i (y_i - mu_j)^2) / var_total
     * Using single-pass identity:
     *     sum_{i in j} w_i (y_i - mu_j)^2 = binsumwy2[j] - binsumwy[j]^2 / binsumw[j] */
    s_sq = 0.0;
    for (i = 0; i < h->Nbin; i++) {
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
    for (i = 0; i < h->Nbin; i++) {
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
      double xn = X * h->Nbin;
      ibin = (int)xn;
      if (ibin >= h->Nbin) ibin = h->Nbin - 1;
      if (ibin < 0) ibin = 0;
      alpha = xn - (double)ibin - 0.5;
      if (alpha < 0.0) {
        ibin_lo = ibin - 1; if (ibin_lo < 0)         ibin_lo += h->Nbin;
        ibin_hi = ibin;
        w_hi = alpha + 1.0;                /* xn - (lo_centre = ibin - 0.5) */
      } else {
        ibin_lo = ibin;
        ibin_hi = ibin + 1; if (ibin_hi >= h->Nbin)  ibin_hi -= h->Nbin;
        w_hi = alpha;                      /* xn - (lo_centre = ibin + 0.5) */
      }
      mu_interp = (1.0 - w_hi) * h->binmean[ibin_lo] + w_hi * h->binmean[ibin_hi];
      Y = mag[i] - mu_interp;
      s_sq += wi * Y * Y;
    }
    if (s_sq < 0.0) s_sq = 0.0;
    return s_sq / var_total;
  }

  /* Unsupported kind in this commit */
  return PDM_ERROR_SCORE;
}

/* Compute the periodogram (theta vs period) over a frequency sweep. */
static void pdm_periodogram(int kind, int N, double *t, double *mag, double *w,
                            double mu_total, double var_total,
                            int Nperiod, double *periods, double *periodogram,
                            _PDMHistType *h)
{
  int k;
  for (k = 0; k < Nperiod; k++)
    periodogram[k] = pdm_theta(kind, N, t, mag, w, mu_total, var_total,
                               periods[k], h);
}

/* ---- Find Npeaks lowest-theta periods, with fine-tune + multiple check ---- */

void findPeaks_pdm(double *t_, double *mag_, double *sig_, int N,
                   int kind, int Nbin, int useerr,
                   double minP, double maxP, double subsample, double finetune,
                   int Npeaks,
                   double *perpeaks, double *thetapeaks, double *peakSNR,
                   double *peakFAP,
                   double clip, int clipiter,
                   double *ave_theta, double *rms_theta,
                   int outflag, char *outname, int ascii)
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
    return;
  }

  /* Use the configured Nbin (or default if <=0) */
  if (Nbin <= 0) Nbin = PDM_DEFAULT_NBIN;
  pdm_hist_alloc(&h, Nbin);

  /* Copy t, mag, and build weight vector. */
  t   = (double *)malloc(N * sizeof(double));
  mag = (double *)malloc(N * sizeof(double));
  w   = (double *)malloc(N * sizeof(double));
  if (t == NULL || mag == NULL || w == NULL) vt_error(ERR_MEMALLOC);

  Nused = 0;
  for (i = 0; i < N; i++) {
    if (isnan(mag_[i])) continue;
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
    free(t); free(mag); free(w);
    pdm_hist_free(&h);
    return;
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

  /* Compute periodogram */
  pdm_periodogram(kind, Nused, t, mag, w, mu_total, var_total,
                  Nperiod, periods, periodogram, &h);

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
                        testperiod, &h);
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
                              testperiod, &h);
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
      if (b_ > 0.0 && a_ > 0.0) {
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

  free(t); free(mag); free(w);
  free(periods); free(periodogram);
  pdm_hist_free(&h);
}

/* ---- Entry point invoked from processcommand.c ---- */

void RunPDMCommand(ProgramData *p, Command *c, _PDM *Pdm, int lcnum, int lc_name_num, int thisindex)
{
  char outname[MAXLEN];
  int i1, i2;

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

  findPeaks_pdm(p->t[lcnum], p->mag[lcnum], p->sig[lcnum], p->NJD[lcnum],
                Pdm->kind, Pdm->Nbin_vals[lcnum], Pdm->useerr,
                Pdm->minp_vals[lcnum], Pdm->maxp_vals[lcnum],
                Pdm->subsample_vals[lcnum], Pdm->finetune_vals[lcnum],
                Pdm->Npeaks,
                Pdm->peakperiods[lcnum], Pdm->peakvalues[lcnum],
                Pdm->peakSNR[lcnum], Pdm->peakFAP[lcnum],
                Pdm->clip, Pdm->clipiter,
                &Pdm->avetheta[lcnum], &Pdm->rmstheta[lcnum],
                Pdm->operiodogram, outname, p->ascii);
}
