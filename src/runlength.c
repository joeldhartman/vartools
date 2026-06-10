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

/* -runlength: run-length statistics about the median and the +/-k*MAD band.
 *
 * For each light curve let
 *     m = median(x),   D = MAD(x) = 1.483 * median(|x_i - m|)
 * over the surviving points (NaN magnitudes and, when "maskpoints" is in
 * effect, points with maskvar <= VARTOOLS_MASK_TINY, are removed first).
 * D is the same 1.483-scaled estimator reported by -stats, so the band
 * +/-k*D is roughly +/-k Gaussian sigma.
 *
 * A "run" is a maximal block of consecutive points (in the stored, time-
 * sorted order; require_sort = 1 is set in the parser so the light curve is
 * sorted in time before this kernel runs) that all satisfy a predicate.
 * Four predicates are tracked, each on the current point order:
 *     above   : x_i > m
 *     below   : x_i < m
 *     outhigh : x_i - m >  k*D
 *     outlow  : x_i - m < -k*D
 * Strict comparisons: a point exactly at the median is "in band" and breaks
 * both the above and below runs.  The two outlier predicates are sign-
 * specific, so a band excursion that crosses from the high side to the low
 * side ends one run and starts another -- there is no combined sign-agnostic
 * outlier run.
 *
 * For each predicate a single O(N) scan accumulates {cur_len, max_len,
 * nruns, npoints}: when the predicate holds, cur_len++ (and nruns++ if
 * cur_len had been 0) and npoints++; otherwise cur_len resets to 0.  The
 * longest run is max_len; the mean run length is npoints/nruns (the run
 * lengths sum to npoints, so no separate sum accumulator is needed), and is
 * reported as NaN when nruns == 0.
 *
 * Complexity note: the run SCAN is a single O(N) pass, but computing m and D
 * requires median selection (O(N) average via quickselect inside median()
 * and MAD()), so the command as a whole is not a single streaming pass.
 *
 * Output (N = command number): RUNLENGTH_{ABOVE,BELOW,OUTHIGH,OUTLOW}_
 * {MAXLEN,NRUNS,MEANLEN}_N plus RUNLENGTH_MEDIAN_N, RUNLENGTH_MAD_N (the
 * 1.483-scaled D) and RUNLENGTH_K_N (the k used).  The run statistics are
 * 0 / 0 / NaN when no point survives filtering.
 */

/* Index labels for the four run categories. */
#define RL_ABOVE   0
#define RL_BELOW   1
#define RL_OUTHIGH 2
#define RL_OUTLOW  3
#define RL_NCAT    4

/* Empty result (no surviving points): zero runs, NaN means, NaN med/MAD.
 * k is still reported as the value that was evaluated for this light curve. */
static void runlength_set_empty(_Runlength *rl, int lcnum, double k)
{
  double nan = sqrt(-1.0);
  rl->above_maxlen[lcnum]   = 0; rl->above_nruns[lcnum]   = 0;
  rl->above_meanlen[lcnum]  = nan;
  rl->below_maxlen[lcnum]   = 0; rl->below_nruns[lcnum]   = 0;
  rl->below_meanlen[lcnum]  = nan;
  rl->outhigh_maxlen[lcnum] = 0; rl->outhigh_nruns[lcnum] = 0;
  rl->outhigh_meanlen[lcnum]= nan;
  rl->outlow_maxlen[lcnum]  = 0; rl->outlow_nruns[lcnum]  = 0;
  rl->outlow_meanlen[lcnum] = nan;
  rl->medval[lcnum] = nan;
  rl->madval[lcnum] = nan;
  rl->kval[lcnum]   = k;
}

void RunRunlengthCommand(ProgramData *p, _Runlength *rl, int lcnum,
                         int lc_name_num)
{
  int i, c, N_in, N;
  double k, med, mad, hi, lo, dev;
  double *xf = NULL;
  int pred;
  int cur[RL_NCAT], maxlen[RL_NCAT], nruns[RL_NCAT], npoints[RL_NCAT];

  N_in = p->NJD[lcnum];

  /* Outlier band half-width in MAD units (fixed, "var" or "expr" source). */
  k = VT_EVAL_DOUBLE(rl, k, lc_name_num, lcnum);

  if(N_in < 1) {
    runlength_set_empty(rl, lcnum, k);
    return;
  }

  if((xf = (double *) malloc(N_in * sizeof(double))) == NULL)
    vt_error(ERR_MEMALLOC);

  /* Filter pass: drop NaN mags and (optionally) masked points, preserving
   * the (time-sorted) point order. */
  N = 0;
  for(i = 0; i < N_in; i++) {
    if(isnan(p->mag[lcnum][i])) continue;
    if(rl->usemask && rl->maskvar != NULL) {
      if(!(EvaluateVariable_Double(lc_name_num, lcnum, i, rl->maskvar) > VARTOOLS_MASK_TINY))
        continue;
    }
    xf[N++] = p->mag[lcnum][i];
  }

  if(N < 1) {
    runlength_set_empty(rl, lcnum, k);
    free(xf);
    return;
  }

  /* median() and MAD() copy their input internally, so xf is preserved for
   * the run scan.  MAD() = 1.483 * median(|x - median|). */
  med = median(N, xf);
  mad = MAD(N, xf);
  hi  =  k * mad;   /* outhigh threshold on (x - med)            */
  lo  = -k * mad;   /* outlow  threshold on (x - med)            */

  for(c = 0; c < RL_NCAT; c++) {
    cur[c] = 0; maxlen[c] = 0; nruns[c] = 0; npoints[c] = 0;
  }

  for(i = 0; i < N; i++) {
    dev = xf[i] - med;
    for(c = 0; c < RL_NCAT; c++) {
      switch(c) {
      case RL_ABOVE:   pred = (xf[i] > med); break;
      case RL_BELOW:   pred = (xf[i] < med); break;
      case RL_OUTHIGH: pred = (dev > hi);    break;
      default:         pred = (dev < lo);    break; /* RL_OUTLOW */
      }
      if(pred) {
        if(cur[c] == 0) nruns[c]++;
        cur[c]++;
        npoints[c]++;
        if(cur[c] > maxlen[c]) maxlen[c] = cur[c];
      } else {
        cur[c] = 0;
      }
    }
  }

  rl->above_maxlen[lcnum]   = maxlen[RL_ABOVE];
  rl->above_nruns[lcnum]    = nruns[RL_ABOVE];
  rl->above_meanlen[lcnum]  = nruns[RL_ABOVE]
                              ? (double) npoints[RL_ABOVE] / nruns[RL_ABOVE]
                              : sqrt(-1.0);
  rl->below_maxlen[lcnum]   = maxlen[RL_BELOW];
  rl->below_nruns[lcnum]    = nruns[RL_BELOW];
  rl->below_meanlen[lcnum]  = nruns[RL_BELOW]
                              ? (double) npoints[RL_BELOW] / nruns[RL_BELOW]
                              : sqrt(-1.0);
  rl->outhigh_maxlen[lcnum] = maxlen[RL_OUTHIGH];
  rl->outhigh_nruns[lcnum]  = nruns[RL_OUTHIGH];
  rl->outhigh_meanlen[lcnum]= nruns[RL_OUTHIGH]
                              ? (double) npoints[RL_OUTHIGH] / nruns[RL_OUTHIGH]
                              : sqrt(-1.0);
  rl->outlow_maxlen[lcnum]  = maxlen[RL_OUTLOW];
  rl->outlow_nruns[lcnum]   = nruns[RL_OUTLOW];
  rl->outlow_meanlen[lcnum] = nruns[RL_OUTLOW]
                              ? (double) npoints[RL_OUTLOW] / nruns[RL_OUTLOW]
                              : sqrt(-1.0);
  rl->medval[lcnum] = med;
  rl->madval[lcnum] = mad;
  rl->kval[lcnum]   = k;

  free(xf);
}
