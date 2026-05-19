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

/* Robust scatter statistics from the magnitude distribution.
 *
 * For each pair (p, q) with 0 < p < q < 100 (canonicalized at parse time):
 *
 *   amp_p_q  = pct(q) - pct(p)
 *   asym_p_q = (pct(q) - median) / (median - pct(p))
 *
 * Plus one fixed scalar independent of the pair list:
 *
 *   medmeddev_over_stddev = median(|x - median(x)|) / stddev(x)
 *
 * For independent Gaussian noise medmeddev/stddev -> 0.6745 in the
 * large-N limit, and asym_p_q -> 1.0 for any symmetric distribution.
 *
 * Percentile interpolation matches the percentile() helper in
 * statistics.c (used by -stats), so PERCENTILERATIOS values are
 * directly comparable to the corresponding pct(p) columns from -stats.
 *
 * NaN magnitudes and (when the "maskpoints" keyword is in effect)
 * points with maskvar <= 0 are filtered out before any statistic is
 * computed.  Zero denominators (e.g. median == pct(p)) and degenerate
 * inputs (post-filter N < 2, or stddev == 0) produce NaN outputs.
 */
void RunPercentileratiosCommand(ProgramData *p, _Percentileratios *pr,
                                 int lcnum, int lc_name_num)
{
  int i, k, pp;
  int N_in, nrej;
  double *buf;
  double medval, sigval, madval, pv, qv, denom;

  N_in = p->NJD[lcnum];

  if (N_in < 2) {
    pr->medmeddev_over_stddev[lcnum] = sqrt(-1.0);
    for (pp = 0; pp < pr->Npairs; pp++) {
      pr->amp[lcnum][pp] = sqrt(-1.0);
      pr->asym[lcnum][pp] = sqrt(-1.0);
    }
    return;
  }

  if ((buf = (double *) malloc(N_in * sizeof(double))) == NULL)
    vt_error(ERR_MEMALLOC);

  nrej = 0;
  for (i = 0; i < N_in; i++) {
    if (isnan(p->mag[lcnum][i])) continue;
    if (pr->usemask && pr->maskvar != NULL) {
      if (!(EvaluateVariable_Double(lc_name_num, lcnum, i, pr->maskvar) > VARTOOLS_MASK_TINY))
        continue;
    }
    buf[nrej++] = p->mag[lcnum][i];
  }

  if (nrej < 2) {
    pr->medmeddev_over_stddev[lcnum] = sqrt(-1.0);
    for (pp = 0; pp < pr->Npairs; pp++) {
      pr->amp[lcnum][pp] = sqrt(-1.0);
      pr->asym[lcnum][pp] = sqrt(-1.0);
    }
    free(buf);
    return;
  }

  medval = median(nrej, buf);
  sigval = stddev(nrej, buf);
  madval = medmeddev(nrej, buf);

  pr->medmeddev_over_stddev[lcnum] =
    (sigval > 0.0) ? (madval / sigval) : sqrt(-1.0);

  for (pp = 0; pp < pr->Npairs; pp++) {
    pv = percentile(nrej, buf, pr->plow[pp]);
    qv = percentile(nrej, buf, pr->phigh[pp]);
    pr->amp[lcnum][pp] = qv - pv;
    denom = medval - pv;
    pr->asym[lcnum][pp] = (denom != 0.0)
                           ? ((qv - medval) / denom)
                           : sqrt(-1.0);
  }

  free(buf);
}
