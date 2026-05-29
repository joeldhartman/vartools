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

/* -CodyQ: quasi-periodicity statistic Q from Cody et al. 2014, AJ, 147, 82
 * (Eq. 6, Sec. 6.2):
 *
 *     Q = (rms_resid^2 - sigma^2) / (rms_raw^2 - sigma^2)
 *
 * where rms_raw is the rms (standard deviation) of the long-term-detrended
 * light curve, rms_resid is the rms of the same curve after subtracting a
 * phase-folded boxcar-smoothed model at the user-supplied period, and
 * sigma^2 is the mean of the per-point squared errors over the points
 * used.  Q ~ 0 for strictly periodic light curves, Q ~ 0.15-0.6 for
 * quasi-periodic ones, and Q ~ 0.6-1 (or with the denominator
 * approaching zero) for aperiodic ones.
 *
 * Procedure (following the paper):
 *   1. Reject NaN magnitudes, points with sig <= 0, and -- when
 *      "maskpoints" is in effect -- points with maskvar <= 0.  This
 *      matches the vartools rms() convention.
 *   2. Long-term detrend: subtract a boxcar-mean smooth of full width
 *      trendwindow (time-axis units).  Cody used 10 days for CoRoT.
 *   3. rms_raw = stddev of the detrended curve.
 *   4. sigma^2 = mean(sig_i^2) over the N' points used.
 *   5. Phase-fold each point to phi_i = ((t_i - t[0])/P) mod 1.  The
 *      phase-domain boxcar smoother below is invariant under a global
 *      phase shift, so the choice of T0 does not affect Q -- t[0] is
 *      used and is not exposed as a parameter.
 *   6. Build the smoothed phase model: model(phi_i) = boxcar mean of
 *      the detrended values of all points within circular phase
 *      distance <= phasesmooth/2 of phi_i.  Implemented as an
 *      O(N log N) two-pointer sliding-window sweep over a tripled
 *      phase array (the phases extended to [-1,0), [0,1), [1,2) so the
 *      circular window for any centre in [0,1) is a contiguous range).
 *   7. resid_i = detrended_i - model_i; rms_resid = stddev(resid).
 *   8. Q from Eq. 6.  When (rms_raw^2 - sigma^2) <= 0 (e.g. a quiet
 *      light curve whose total scatter is at or below the noise floor)
 *      Q is reported as NaN.
 *
 * The period is resolved in processcommand.c from the period-source
 * keyword (aov/ls/bls/pdm/ftp/injectharm/fix/fixcolumn/list/expr/var)
 * and is stored in cq->period[lcnum][0] before this kernel runs.
 *
 * trendwindow and phasesmooth each accept a fixed value or a per-LC
 * "var"/"expr" source; the values are resolved here per LC.
 *
 * require_sort = 1 is set in the parser, so p->t[lcnum] is sorted in
 * time before this kernel runs -- the time-domain boxcar assumes that
 * ordering.
 *
 * Edge cases (Q reported as NaN):
 *   - fewer than 2 points survive filtering;
 *   - a non-finite or non-positive trendwindow;
 *   - phasesmooth not in (0, 1];
 *   - a non-finite or non-positive resolved period (eg. a bad
 *     var/expr source);
 *   - the Q denominator (rms_raw^2 - sigma^2) is <= 0.
 * RMS_raw, RMS_resid, Sigma and Npoints remain meaningful and are still
 * reported when computable (eg. when only the Q denominator is bad).
 */

static void codyq_set_nan(_CodyQ *cq, int lcnum)
{
  cq->Q[lcnum]         = sqrt(-1.0);
  cq->RMS_raw[lcnum]   = sqrt(-1.0);
  cq->RMS_resid[lcnum] = sqrt(-1.0);
  cq->Sigma[lcnum]     = sqrt(-1.0);
  cq->Npoints[lcnum]   = 0;
}

/* out[i] = mean of val[j] over all j with |t[j] - t[i]| <= halfwidth.
 * t must be sorted ascending.  O(N) running-sum sweep. */
static void codyq_boxcar_time_mean(int N, const double *t, const double *val,
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

/* Sort record for the phase-domain smoother: stores the wrapped phase
 * of a point along with the detrended value and the original index, so
 * the sliding-window model values can be unshuffled back to time order
 * after the sweep. */
typedef struct {
  double phi;
  double val;
  int orig_idx;
} codyq_sortrec_t;

static int codyq_cmp_sortrec(const void *a, const void *b)
{
  double pa = ((const codyq_sortrec_t *) a)->phi;
  double pb = ((const codyq_sortrec_t *) b)->phi;
  if(pa < pb) return -1;
  if(pa > pb) return  1;
  return 0;
}

/* Phase-domain circular boxcar mean.  For each point i in the sorted
 * order, model[orig_idx[i]] = mean of val over all points whose phase
 * is within hw of ps[i] on the unit circle.  ps must already be in
 * [0, 1) and sorted ascending.  O(N log N) overall (the sort dominates;
 * the sliding-window sweep here is O(N)). */
static void codyq_boxcar_phase_mean(int N, const codyq_sortrec_t *recs,
                                    double hw, double *model)
{
  int i, k, lo, hi;
  double sum;
  double *phi3 = NULL, *val3 = NULL;

  if(hw >= 0.5) {
    /* Window covers the whole period -- model is the global mean. */
    double s = 0.0;
    for(i = 0; i < N; i++) s += recs[i].val;
    s /= (double) N;
    for(i = 0; i < N; i++) model[recs[i].orig_idx] = s;
    return;
  }

  if((phi3 = (double *) malloc(3 * N * sizeof(double))) == NULL ||
     (val3 = (double *) malloc(3 * N * sizeof(double))) == NULL)
    vt_error(ERR_MEMALLOC);

  /* Tripled-phase extension: copies at phi-1, phi, phi+1.  The middle
   * copy is the one whose model values we keep; the lower and upper
   * copies provide the wraparound on either side, so for any centre
   * phi in [0,1) and any hw < 0.5 the circular window [phi-hw, phi+hw]
   * is a contiguous range in phi3. */
  for(i = 0; i < N; i++) {
    phi3[i]         = recs[i].phi - 1.0;
    phi3[i + N]     = recs[i].phi;
    phi3[i + 2*N]   = recs[i].phi + 1.0;
    val3[i]         = recs[i].val;
    val3[i + N]     = recs[i].val;
    val3[i + 2*N]   = recs[i].val;
  }

  /* Two-pointer sweep over the tripled array.  We only retain model
   * values for the middle copy (k in [N, 2N)). */
  lo = 0; hi = 0; sum = 0.0;
  for(k = 0; k < 3 * N; k++) {
    while(hi < 3 * N && phi3[hi] <= phi3[k] + hw) { sum += val3[hi]; hi++; }
    while(lo < hi && phi3[lo] < phi3[k] - hw)    { sum -= val3[lo]; lo++; }
    if(k >= N && k < 2 * N) {
      int orig = recs[k - N].orig_idx;
      model[orig] = sum / (double)(hi - lo);
    }
  }

  free(phi3);
  free(val3);
}

void RunCodyQCommand(ProgramData *p, _CodyQ *cq, int lcnum, int lc_name_num)
{
  int i, N_in, N_filt;
  double Wt, PS, P, hw;
  double *t_filt = NULL, *m_filt = NULL, *s_filt = NULL;
  double *smooth = NULL, *detrend = NULL, *model = NULL, *resid = NULL;
  codyq_sortrec_t *recs = NULL;
  double rms_raw, rms_resid, sigma2, denom, Qval;

  N_in = p->NJD[lcnum];

  /* Per-LC parameter values (fixed, "var" or "expr" source). */
  Wt = VT_EVAL_DOUBLE(cq, trendwindow, lc_name_num, lcnum);
  PS = VT_EVAL_DOUBLE(cq, phasesmooth, lc_name_num, lcnum);
  P  = cq->period[lcnum][0];

  /* Input validation: any bad input -> all-NaN for this LC, no crash. */
  if(N_in < 2 || !(Wt > 0.0) || !(PS > 0.0) || PS > 1.0 || !(P > 0.0)) {
    codyq_set_nan(cq, lcnum);
    return;
  }

  if((t_filt  = (double *) malloc(N_in * sizeof(double))) == NULL ||
     (m_filt  = (double *) malloc(N_in * sizeof(double))) == NULL ||
     (s_filt  = (double *) malloc(N_in * sizeof(double))) == NULL ||
     (smooth  = (double *) malloc(N_in * sizeof(double))) == NULL ||
     (detrend = (double *) malloc(N_in * sizeof(double))) == NULL ||
     (model   = (double *) malloc(N_in * sizeof(double))) == NULL ||
     (resid   = (double *) malloc(N_in * sizeof(double))) == NULL ||
     (recs    = (codyq_sortrec_t *) malloc(N_in * sizeof(codyq_sortrec_t))) == NULL)
    vt_error(ERR_MEMALLOC);

  /* Filter pass: drop NaN mags, sig <= 0, and (optionally) masked points. */
  N_filt = 0;
  for(i = 0; i < N_in; i++) {
    if(isnan(p->mag[lcnum][i])) continue;
    if(!(p->sig[lcnum][i] > 0.0)) continue;
    if(cq->usemask && cq->maskvar != NULL) {
      if(!(EvaluateVariable_Double(lc_name_num, lcnum, i, cq->maskvar) > VARTOOLS_MASK_TINY))
        continue;
    }
    t_filt[N_filt] = p->t[lcnum][i];
    m_filt[N_filt] = p->mag[lcnum][i];
    s_filt[N_filt] = p->sig[lcnum][i];
    N_filt++;
  }

  if(N_filt < 2) {
    codyq_set_nan(cq, lcnum);
    goto cleanup;
  }

  /* Long-term detrend: subtract a boxcar-mean smooth of full width Wt. */
  codyq_boxcar_time_mean(N_filt, t_filt, m_filt, 0.5 * Wt, smooth);
  for(i = 0; i < N_filt; i++)
    detrend[i] = m_filt[i] - smooth[i];

  rms_raw = stddev(N_filt, detrend);

  /* sigma^2 = mean of per-point squared errors over the points used. */
  {
    double sumsq = 0.0;
    for(i = 0; i < N_filt; i++) sumsq += s_filt[i] * s_filt[i];
    sigma2 = sumsq / (double) N_filt;
  }

  /* Phase-fold to [0, 1). */
  for(i = 0; i < N_filt; i++) {
    double phi = (t_filt[i] - t_filt[0]) / P;
    phi = phi - floor(phi);
    recs[i].phi      = phi;
    recs[i].val      = detrend[i];
    recs[i].orig_idx = i;
  }
  qsort(recs, (size_t) N_filt, sizeof(codyq_sortrec_t), codyq_cmp_sortrec);

  /* Smoothed phase model (writes into model[] indexed by orig_idx). */
  hw = 0.5 * PS;
  codyq_boxcar_phase_mean(N_filt, recs, hw, model);

  /* Residuals and rms_resid. */
  for(i = 0; i < N_filt; i++)
    resid[i] = detrend[i] - model[i];
  rms_resid = stddev(N_filt, resid);

  denom = rms_raw * rms_raw - sigma2;
  if(!(denom > 0.0)) {
    /* Quiet light curve below noise floor: Q is undefined but the
     * diagnostics remain meaningful. */
    Qval = sqrt(-1.0);
  } else {
    Qval = (rms_resid * rms_resid - sigma2) / denom;
  }

  cq->Q[lcnum]         = Qval;
  cq->RMS_raw[lcnum]   = rms_raw;
  cq->RMS_resid[lcnum] = rms_resid;
  cq->Sigma[lcnum]     = sqrt(sigma2);
  cq->Npoints[lcnum]   = N_filt;

cleanup:
  free(t_filt); free(m_filt); free(s_filt);
  free(smooth); free(detrend); free(model); free(resid);
  free(recs);
}
