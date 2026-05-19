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

/* Fraction of magnitudes more than N*sigma above / below the median,
 * for a user-supplied list of N values.
 *
 *   threshold(N) = N * sigma
 *   frac_above_N = #{ x : x > median + threshold } / N_rej
 *   frac_below_N = #{ x : x < median - threshold } / N_rej
 *
 * sigma is the sample standard deviation by default, or
 * 1.483 * median(|x - median(x)|) if the useMAD flag is given (vartools'
 * MAD(), the Gaussian-consistent calibration).
 *
 * Cf. the Beyond1Std feature of Nun et al. 2015, "FATS: Feature
 * Analysis for Time Series" (arXiv:1506.00010), which is the N=1
 * instance of this statistic.
 *
 * NaN magnitudes and (when the "maskpoints" keyword is in effect)
 * points with maskvar <= 0 are filtered out before any statistic is
 * computed.
 * Edge cases:
 *   N_rej < 2     -> NaN for all fractions
 *   sigma == 0    -> 0 for all fractions (degenerate but well-defined;
 *                    nothing strictly exceeds a zero threshold)
 *
 * Comparisons use strict > and <.
 */
void RunBeyondNsigmaCommand(ProgramData *p, _BeyondNsigma *bs,
                            int lcnum, int lc_name_num)
{
  int i, k;
  int N_in, nrej;
  double *buf;
  double medval, sigval, dev, thr;

  N_in = p->NJD[lcnum];

  if (N_in < 2) {
    for (k = 0; k < bs->NN; k++) {
      bs->frac_above[lcnum][k] = sqrt(-1.0);
      bs->frac_below[lcnum][k] = sqrt(-1.0);
    }
    return;
  }

  if ((buf = (double *) malloc(N_in * sizeof(double))) == NULL)
    vt_error(ERR_MEMALLOC);

  nrej = 0;
  for (i = 0; i < N_in; i++) {
    if (isnan(p->mag[lcnum][i])) continue;
    if (bs->usemask && bs->maskvar != NULL) {
      if (!(EvaluateVariable_Double(lc_name_num, lcnum, i, bs->maskvar) > VARTOOLS_MASK_TINY))
        continue;
    }
    buf[nrej++] = p->mag[lcnum][i];
  }

  if (nrej < 2) {
    for (k = 0; k < bs->NN; k++) {
      bs->frac_above[lcnum][k] = sqrt(-1.0);
      bs->frac_below[lcnum][k] = sqrt(-1.0);
    }
    free(buf);
    return;
  }

  medval = median(nrej, buf);
  sigval = bs->useMAD ? MAD(nrej, buf) : stddev(nrej, buf);

  if (sigval == 0.0) {
    for (k = 0; k < bs->NN; k++) {
      bs->frac_above[lcnum][k] = 0.0;
      bs->frac_below[lcnum][k] = 0.0;
    }
    free(buf);
    return;
  }

  for (k = 0; k < bs->NN; k++) {
    bs->frac_above[lcnum][k] = 0.0;
    bs->frac_below[lcnum][k] = 0.0;
  }

  for (i = 0; i < nrej; i++) {
    dev = buf[i] - medval;
    for (k = 0; k < bs->NN; k++) {
      thr = bs->Nvalues[k] * sigval;
      if (dev > thr) bs->frac_above[lcnum][k] += 1.0;
      else if (dev < -thr) bs->frac_below[lcnum][k] += 1.0;
    }
  }

  for (k = 0; k < bs->NN; k++) {
    bs->frac_above[lcnum][k] /= (double) nrej;
    bs->frac_below[lcnum][k] /= (double) nrej;
  }

  free(buf);
}
