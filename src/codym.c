/*     This file is part of VARTOOLS version 1.40                            */
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
/*     Copyright 2026 Joel Hartman                                           */
/*                                                                           */
#include "commands.h"
#include "programdata.h"
#include "functions.h"

/* -CodyM: flux-asymmetry statistic M from Cody et al. 2014, AJ, 147, 82
 * (Eq. 7, Sec. 6.3):
 *
 *     M = (<d_10%> - d_med) / sigma_d
 *
 * computed on the long-term-detrended, outlier-filtered light curve, where
 *   d_med    = median of the curve,
 *   sigma_d  = rms (standard deviation) of the curve,
 *   <d_10%>  = mean of the combined faintest-decile + brightest-decile values.
 *
 * For magnitude-valued light curves M > 0.25 indicates a "dipping" curve
 * (asymmetric toward faint excursions), M < -0.25 a "bursting" curve, and
 * |M| < 0.25 a symmetric one.  The sign convention flips for flux input.
 *
 * Procedure (following the paper):
 *   1. Reject NaN magnitudes and, when "maskpoints" is in effect, points
 *      with maskvar <= VARTOOLS_MASK_TINY.
 *   2. Long-term detrend: subtract an inline boxcar-mean smooth of full
 *      width trendwindow (time-axis units).  M is only meaningful for
 *      asymmetry on timescales shorter than ~half the light-curve
 *      duration; the detrend removes longer trends.
 *   3. Identify sigma-clip outliers:
 *        - two-stage (paper-faithful): when "outlierwindow" is given, build
 *          a residual = mag - (boxcar-mean smooth of full width
 *          outlierwindow) and flag points deviating from the residual mean
 *          by more than sigclip * rms(residual).  Real variability is
 *          removed in this residual, so only genuine glitches are flagged
 *          -- the deep dips/bursts that M measures survive.
 *        - single-stage: when "outlierwindow" is absent, flag points on
 *          the trend-detrended curve directly.  Simpler, but a deep
 *          dip/burst can clip itself and bias M toward 0.
 *      sigclip <= 0 disables outlier rejection entirely.
 *   4. On the trend-detrended curve with the flagged points removed
 *      (N' points): d_med = median, sigma_d = stddev, and <d_10%> = mean
 *      of the n_dec faintest plus n_dec brightest values, with
 *      n_dec = round(0.1 * N').
 *   5. M = (<d_10%> - d_med) / sigma_d.
 *
 * trendwindow, outlierwindow and sigclip each accept a fixed value or a
 * per-light-curve "var"/"expr" source; the values are resolved here per LC.
 *
 * require_sort = 1 is set in the parser, so p->t[lcnum] is sorted in time
 * before this kernel runs -- the boxcar smoothers assume that ordering.
 *
 * Edge cases (M and the dependent diagnostics are reported as NaN):
 *   - fewer than 2 points survive filtering / outlier rejection;
 *   - a non-finite or non-positive trendwindow, or (two-stage) a
 *     non-finite or non-positive outlierwindow;
 *   - sigma_d == 0 (a flat curve);
 *   - the two deciles would overlap (2 * n_dec > N').
 */

/* Write NaN to every output cell for this light curve. */
static void codym_set_nan(_CodyM *cm, int lcnum)
{
  cm->M[lcnum]       = sqrt(-1.0);
  cm->d10[lcnum]     = sqrt(-1.0);
  cm->dmed[lcnum]    = sqrt(-1.0);
  cm->sigma_d[lcnum] = sqrt(-1.0);
  cm->Npoints[lcnum] = 0;
}

/* out[i] = mean of val[j] over all j with |t[j] - t[i]| <= halfwidth.
 * t must be sorted ascending.  O(N): both window edges advance
 * monotonically as i sweeps upward.  Point i is always within its own
 * window, so the divisor (hi - lo) is >= 1. */
static void codym_boxcar_mean(int N, const double *t, const double *val,
                              double halfwidth, double *out)
{
  int i, lo = 0, hi = 0;
  double sum = 0.0;
  for(i = 0; i < N; i++) {
    while(hi < N && t[hi] <= t[i] + halfwidth) { sum += val[hi]; hi++; }
    while(lo < hi && t[lo] < t[i] - halfwidth) { sum -= val[lo]; lo++; }
    out[i] = sum / (double)(hi - lo);
  }
}

/* Ascending comparator for qsort over doubles. */
static int codym_cmp_double(const void *a, const void *b)
{
  double da = *(const double *) a, db = *(const double *) b;
  if(da < db) return -1;
  if(da > db) return  1;
  return 0;
}

void RunCodyMCommand(ProgramData *p, _CodyM *cm, int lcnum, int lc_name_num)
{
  int i, N_in, N_filt, N_keep, n_dec;
  int do_reject, two_stage;
  double Wt, Wo, S;
  double *t_filt = NULL, *m_filt = NULL, *smooth = NULL;
  double *detrend = NULL, *keep = NULL;
  double resid_mean, resid_rms, d_med, sigma_d, d10, Mval;

  N_in = p->NJD[lcnum];

  /* Resolve the per-LC parameter values (fixed, "var" or "expr" source). */
  Wt = VT_EVAL_DOUBLE(cm, trendwindow, lc_name_num, lcnum);
  S  = VT_EVAL_DOUBLE(cm, sigclip,     lc_name_num, lcnum);
  Wo = cm->has_outlierwindow
       ? VT_EVAL_DOUBLE(cm, outlierwindow, lc_name_num, lcnum) : 0.0;

  do_reject = (S > 0.0);
  two_stage = (do_reject && cm->has_outlierwindow);

  /* The detrend window must be finite and positive; in the two-stage case
   * the outlier window must be too. */
  if(N_in < 2 || !(Wt > 0.0) || (two_stage && !(Wo > 0.0))) {
    codym_set_nan(cm, lcnum);
    return;
  }

  if((t_filt  = (double *) malloc(N_in * sizeof(double))) == NULL ||
     (m_filt  = (double *) malloc(N_in * sizeof(double))) == NULL ||
     (smooth  = (double *) malloc(N_in * sizeof(double))) == NULL ||
     (detrend = (double *) malloc(N_in * sizeof(double))) == NULL ||
     (keep    = (double *) malloc(N_in * sizeof(double))) == NULL)
    vt_error(ERR_MEMALLOC);

  /* Filter pass: drop NaN mags and (optionally) masked points. */
  N_filt = 0;
  for(i = 0; i < N_in; i++) {
    if(isnan(p->mag[lcnum][i])) continue;
    if(cm->usemask && cm->maskvar != NULL) {
      if(!(EvaluateVariable_Double(lc_name_num, lcnum, i, cm->maskvar) > VARTOOLS_MASK_TINY))
        continue;
    }
    t_filt[N_filt] = p->t[lcnum][i];
    m_filt[N_filt] = p->mag[lcnum][i];
    N_filt++;
  }

  if(N_filt < 2) {
    codym_set_nan(cm, lcnum);
    free(t_filt); free(m_filt); free(smooth); free(detrend); free(keep);
    return;
  }

  /* Long-term detrend: subtract a boxcar-mean smooth of full width Wt. */
  codym_boxcar_mean(N_filt, t_filt, m_filt, 0.5 * Wt, smooth);
  for(i = 0; i < N_filt; i++)
    detrend[i] = m_filt[i] - smooth[i];

  /* Build the curve of kept (non-outlier) detrended values. */
  if(do_reject) {
    /* Residual on which outliers are judged: a short-smooth residual for
     * the two-stage scheme, otherwise the trend-detrended curve itself. */
    if(two_stage) {
      codym_boxcar_mean(N_filt, t_filt, m_filt, 0.5 * Wo, smooth);
      for(i = 0; i < N_filt; i++)
        smooth[i] = m_filt[i] - smooth[i];
    } else {
      for(i = 0; i < N_filt; i++)
        smooth[i] = detrend[i];
    }
    resid_mean = getmean(N_filt, smooth);
    resid_rms  = stddev(N_filt, smooth);
    N_keep = 0;
    if(resid_rms > 0.0) {
      double thr = S * resid_rms;
      for(i = 0; i < N_filt; i++) {
        if(fabs(smooth[i] - resid_mean) <= thr)
          keep[N_keep++] = detrend[i];
      }
    } else {
      /* A zero-rms residual: no point can exceed a positive threshold. */
      for(i = 0; i < N_filt; i++)
        keep[N_keep++] = detrend[i];
    }
  } else {
    for(i = 0; i < N_filt; i++)
      keep[i] = detrend[i];
    N_keep = N_filt;
  }

  if(N_keep < 2) {
    codym_set_nan(cm, lcnum);
    free(t_filt); free(m_filt); free(smooth); free(detrend); free(keep);
    return;
  }

  /* sigma_d is order-independent; compute it before sorting keep[]. */
  sigma_d = stddev(N_keep, keep);

  n_dec = (int) (0.1 * (double) N_keep + 0.5);
  if(n_dec < 1) n_dec = 1;

  if(2 * n_dec > N_keep || !(sigma_d > 0.0)) {
    /* Overlapping deciles or a flat curve: M is undefined, but d_med /
     * sigma_d / Npoints remain meaningful and are still reported. */
    qsort(keep, (size_t) N_keep, sizeof(double), codym_cmp_double);
    cm->M[lcnum]       = sqrt(-1.0);
    cm->d10[lcnum]     = sqrt(-1.0);
    cm->dmed[lcnum]    = (N_keep % 2)
                         ? keep[N_keep/2]
                         : 0.5 * (keep[N_keep/2 - 1] + keep[N_keep/2]);
    cm->sigma_d[lcnum] = sigma_d;
    cm->Npoints[lcnum] = N_keep;
    free(t_filt); free(m_filt); free(smooth); free(detrend); free(keep);
    return;
  }

  /* Sort ascending: d_med from the middle, <d_10%> from the two ends. */
  qsort(keep, (size_t) N_keep, sizeof(double), codym_cmp_double);
  d_med = (N_keep % 2)
          ? keep[N_keep/2]
          : 0.5 * (keep[N_keep/2 - 1] + keep[N_keep/2]);

  d10 = 0.0;
  for(i = 0; i < n_dec; i++)
    d10 += keep[i] + keep[N_keep - 1 - i];
  d10 /= (double) (2 * n_dec);

  Mval = (d10 - d_med) / sigma_d;

  cm->M[lcnum]       = Mval;
  cm->d10[lcnum]     = d10;
  cm->dmed[lcnum]    = d_med;
  cm->sigma_d[lcnum] = sigma_d;
  cm->Npoints[lcnum] = N_keep;

  free(t_filt); free(m_filt); free(smooth); free(detrend); free(keep);
}
