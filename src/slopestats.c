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

/* Per-pair slope statistics over consecutive (after sorting in time, and
 * optionally after binning) light curve points:
 *
 *   s_i = (m[i+1] - m[i]) / (t[i+1] - t[i])
 *
 *   median_abs_dmdt = median(|s_i|)
 *   max_abs_dmdt    = max(|s_i|)
 *   mad_dmdt        = 1.483 * median(|s_i - median(s_i)|)
 *   sigma_s         = stddev(s_i)  (or MAD(s_i) if useMAD)
 *   frac_above_T    = #{ s_i : s_i - <s> >  T * sigma_s } / N_pairs
 *   frac_below_T    = #{ s_i : s_i - <s> < -T * sigma_s } / N_pairs
 *
 * <s> is the mean of the slopes; deviations are mean-centered regardless
 * of which sigma flavor is in use.
 *
 * The "bintime" keyword phase-bins the input (filtered) light curve in
 * time before forming pairs, with one (mean_t, mean_m) entry per occupied
 * bin.  Multiple bin sizes may be given; each one is run independently
 * and yields its own column set.  bintime values are interpreted as
 * minutes, assuming the LC times are in days; the kernel divides each
 * bintime by 1440 before binning.  When no binning is requested the
 * filtered LC is used directly (no extra copy).
 *
 * The "binshift" keyword shifts the start of the first bin by a
 * fraction of the binwidth: t0_bin = t[0] - binshift * binsize.
 * binshift is dimensionless, with canonical use 0 <= binshift
 * < 1 (a value of 0.5 produces a half-bin shift, etc.).  The shift is
 * applied per-binsize independently.  Note this deviates from -binlc's
 * literal `t0 - firstbinshift / binsize` formula in -binlc
 * (binlc.c:497), which appears to be a longstanding off-by-operator
 * bug: with `/` the shift has units of time^2 / time = 1/time, which is
 * dimensionally inconsistent with subtracting from t.  The intent here
 * is for -binlc to eventually grow a "binshift" keyword that matches
 * -slopestats' semantics, with "firstbinshift" retained (but
 * undocumented) for backward compatibility.
 *
 * The "maxgap" keyword drops pairs whose time separation exceeds maxgap.
 * maxgap is interpreted in days (the same units as the LC time axis).
 *
 * Cf. Richards et al. 2011, ApJ 733, 10 ("MaxSlope" = max_abs_dmdt;
 * "Beyond1Std" is the magnitude-domain analog of frac_above/below_T).
 *
 * NaN magnitudes and (when "maskpoints" is in effect) points with
 * maskvar <= VARTOOLS_MASK_TINY are filtered out before any statistic
 * is computed.
 *
 * Edge cases:
 *   N_pairs < 1                 -> NaN for all 5 stats (this binsize)
 *   N_pairs == 1 or sigma == 0  -> frac_above/below = 0; median/max/MAD
 *                                  still well-defined
 */
void RunSlopestatsCommand(ProgramData *p, _Slopestats *ss,
                          int lcnum, int lc_name_num)
{
  int i, b, t, idx;
  int N_in, N_filt, N_cur, N_pairs;
  int N_bin = ss->N_bin;
  int N_thresh = ss->N_thresh;
  double *t_filt = NULL, *m_filt = NULL;
  double *t_binned = NULL, *m_binned = NULL;
  double *slopes = NULL, *aux = NULL;
  double *cur_t, *cur_m;
  double sigval, mean_s;

  N_in = p->NJD[lcnum];

  /* Short-circuit empty input: write NaN everywhere. */
  if (N_in < 1) {
    for (b = 0; b < N_bin; b++) {
      ss->median_abs_dmdt[lcnum][b] = sqrt(-1.0);
      ss->max_abs_dmdt[lcnum][b]    = sqrt(-1.0);
      ss->mad_dmdt[lcnum][b]        = sqrt(-1.0);
      for (t = 0; t < N_thresh; t++) {
        idx = b * N_thresh + t;
        ss->frac_above[lcnum][idx] = sqrt(-1.0);
        ss->frac_below[lcnum][idx] = sqrt(-1.0);
      }
    }
    return;
  }

  if ((t_filt = (double *) malloc(N_in * sizeof(double))) == NULL ||
      (m_filt = (double *) malloc(N_in * sizeof(double))) == NULL ||
      (slopes = (double *) malloc(N_in * sizeof(double))) == NULL ||
      (aux    = (double *) malloc(N_in * sizeof(double))) == NULL)
    vt_error(ERR_MEMALLOC);

  if (ss->has_binning) {
    if ((t_binned = (double *) malloc(N_in * sizeof(double))) == NULL ||
        (m_binned = (double *) malloc(N_in * sizeof(double))) == NULL)
      vt_error(ERR_MEMALLOC);
  }

  /* Filter pass: drop NaN mags and (optionally) masked points. */
  N_filt = 0;
  for (i = 0; i < N_in; i++) {
    if (isnan(p->mag[lcnum][i])) continue;
    if (ss->usemask && ss->maskvar != NULL) {
      if (!(EvaluateVariable_Double(lc_name_num, lcnum, i, ss->maskvar) > VARTOOLS_MASK_TINY))
        continue;
    }
    t_filt[N_filt] = p->t[lcnum][i];
    m_filt[N_filt] = p->mag[lcnum][i];
    N_filt++;
  }

  /* Per-bintime loop. */
  for (b = 0; b < N_bin; b++) {

    if (ss->has_binning) {
      /* bintime is in minutes; convert to days for the binning math.
       * binshift is a unitless fraction of the binwidth. */
      double binsize = ss->bintimes[b] / 1440.0;
      double t0_bin  = (N_filt > 0) ? t_filt[0] : 0.0;
      long bin_idx_cur = 0;
      int started = 0;
      double sum_t = 0.0, sum_m = 0.0;
      int count = 0;
      int N_b = 0;

      if (ss->has_binshift)
        t0_bin -= ss->binshift * binsize;

      /* Single-pass binning: emit (mean_t, mean_m) as each bin closes.
       * Singleton bins are kept.  Empty bins are skipped (we only emit
       * when count > 0). */
      for (i = 0; i < N_filt; i++) {
        long bin_idx_i = (long) floor((t_filt[i] - t0_bin) / binsize);
        if (!started || bin_idx_i != bin_idx_cur) {
          if (count > 0) {
            t_binned[N_b] = sum_t / count;
            m_binned[N_b] = sum_m / count;
            N_b++;
          }
          bin_idx_cur = bin_idx_i;
          started = 1;
          sum_t = 0.0;
          sum_m = 0.0;
          count = 0;
        }
        sum_t += t_filt[i];
        sum_m += m_filt[i];
        count++;
      }
      if (count > 0) {
        t_binned[N_b] = sum_t / count;
        m_binned[N_b] = sum_m / count;
        N_b++;
      }

      cur_t = t_binned;
      cur_m = m_binned;
      N_cur = N_b;
    } else {
      cur_t = t_filt;
      cur_m = m_filt;
      N_cur = N_filt;
    }

    /* Form pairs s_i = dm/dt with inline maxgap filter. */
    N_pairs = 0;
    for (i = 0; i < N_cur - 1; i++) {
      double dt = cur_t[i+1] - cur_t[i];
      if (dt <= 0.0) continue;
      if (ss->has_maxgap && dt > ss->maxgap) continue;
      slopes[N_pairs++] = (cur_m[i+1] - cur_m[i]) / dt;
    }

    if (N_pairs < 1) {
      ss->median_abs_dmdt[lcnum][b] = sqrt(-1.0);
      ss->max_abs_dmdt[lcnum][b]    = sqrt(-1.0);
      ss->mad_dmdt[lcnum][b]        = sqrt(-1.0);
      for (t = 0; t < N_thresh; t++) {
        idx = b * N_thresh + t;
        ss->frac_above[lcnum][idx] = sqrt(-1.0);
        ss->frac_below[lcnum][idx] = sqrt(-1.0);
      }
      continue;
    }

    /* median_abs_dmdt and max_abs_dmdt: aux = |slopes|, then median_nocopy. */
    {
      double maxabs = 0.0;
      for (i = 0; i < N_pairs; i++) {
        aux[i] = fabs(slopes[i]);
        if (aux[i] > maxabs) maxabs = aux[i];
      }
      ss->max_abs_dmdt[lcnum][b]    = maxabs;
      ss->median_abs_dmdt[lcnum][b] = median_nocopy(N_pairs, aux);
    }

    /* mad_dmdt = 1.483 * median(|s - median(s)|), via MAD_nocopy on a copy. */
    memcpy(aux, slopes, N_pairs * sizeof(double));
    ss->mad_dmdt[lcnum][b] = MAD_nocopy(N_pairs, aux);

    /* Mean of slopes (used for centering deviations regardless of sigma flavor). */
    {
      double sum = 0.0;
      for (i = 0; i < N_pairs; i++) sum += slopes[i];
      mean_s = sum / (double) N_pairs;
    }

    /* sigma: stddev (n-1 normalized, matching vartools' stddev()) or MAD. */
    if (ss->useMAD) {
      sigval = ss->mad_dmdt[lcnum][b];
    } else if (N_pairs > 1) {
      double sumsq = 0.0;
      for (i = 0; i < N_pairs; i++) {
        double dev = slopes[i] - mean_s;
        sumsq += dev * dev;
      }
      sigval = sqrt(sumsq / (double) (N_pairs - 1));
    } else {
      sigval = 0.0;
    }

    if (sigval == 0.0) {
      for (t = 0; t < N_thresh; t++) {
        idx = b * N_thresh + t;
        ss->frac_above[lcnum][idx] = 0.0;
        ss->frac_below[lcnum][idx] = 0.0;
      }
      continue;
    }

    for (t = 0; t < N_thresh; t++) {
      idx = b * N_thresh + t;
      ss->frac_above[lcnum][idx] = 0.0;
      ss->frac_below[lcnum][idx] = 0.0;
    }
    for (i = 0; i < N_pairs; i++) {
      double dev = slopes[i] - mean_s;
      for (t = 0; t < N_thresh; t++) {
        double thr = ss->thresholds[t] * sigval;
        idx = b * N_thresh + t;
        if (dev >  thr) ss->frac_above[lcnum][idx] += 1.0;
        else if (dev < -thr) ss->frac_below[lcnum][idx] += 1.0;
      }
    }
    for (t = 0; t < N_thresh; t++) {
      idx = b * N_thresh + t;
      ss->frac_above[lcnum][idx] /= (double) N_pairs;
      ss->frac_below[lcnum][idx] /= (double) N_pairs;
    }
  }

  free(t_filt);
  free(m_filt);
  free(slopes);
  free(aux);
  if (t_binned) free(t_binned);
  if (m_binned) free(m_binned);
}
