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

/* von Neumann (1941) ratio
 *
 *   eta = delta^2 / s^2
 *     delta^2 = (1 / (N-1)) * sum_i (y[i+1] - y[i])^2
 *     s^2     = (1 / N)     * sum_i (y[i] - ybar)^2
 *
 * For independent Gaussian samples, E[eta] = 2 with variance ~ 4/N.  Smooth
 * (correlated) signals drive eta well below 2; anti-correlated alternating
 * signals push eta > 2.  Widely used as a variability indicator for sparse
 * uneven photometric time series (Sokolovsky et al. 2017, MNRAS 464, 274).
 *
 * The result depends on the ordering of points; the caller is responsible for
 * passing a time-sorted light curve (the parsecommandline.c branch for -vonNeumann
 * sets c[cn].require_sort = 1 to ensure this happens upstream).
 *
 * Weighted variant (chi-square inverse-variance weighting):
 *
 *   w_i        = 1 / sigma_i^2
 *   w_pair_i   = 1 / (sigma_i^2 + sigma_{i+1}^2)        (proper for pair difference)
 *   ybar_w     = sum w_i y_i / sum w_i
 *   eta_w      = (2N/(N-1)) * sum w_pair_i (y[i+1] - y[i])^2
 *                           / sum w_i (y_i - ybar_w)^2
 *
 * The (2N/(N-1)) prefactor in eta_w restores E[eta_w] ~ 2 for any sigma
 * distribution under white noise.  The naive ratio
 * sum_pair / Σw_pair / (sum_w / Σw) instead converges to ⟨w⟩/⟨w_pair⟩,
 * which is 2 only for homoscedastic sigma; for non-uniform sigma it drifts
 * (e.g. bimodal sigma in {0.5, 5.0} gives raw eta_w ~ 3.9, corrected ~ 2.0).
 * For homoscedastic sigma the corrected weighted form reduces exactly to the
 * unweighted form.
 *
 * Long-double accumulators throughout for the same numerical-stability reason
 * the alarm statistic uses them.
 */
double dovonneumann(int N_in, double *mag_in, double *err_in,
                     int weighted, int usemask,
                     _Variable *maskvar, int lcnum, int lc_name_num)
{
  int i, N;
  double *mag = NULL, *err = NULL;
  long double sumy = 0.0L, sumw = 0.0L, ybar;
  long double sumvar = 0.0L;
  long double sumdy = 0.0L;
  double delta_sq, s_sq;

  if (N_in < 2) return 0.0;

  /* Build a filtered local copy: drop NaN mag (and NaN / non-positive err when
   * weighting is requested) and points masked out via maskpoints. */
  mag = (double *) malloc(N_in * sizeof(double));
  err = (double *) malloc(N_in * sizeof(double));
  if (mag == NULL || err == NULL) vt_error(ERR_MEMALLOC);
  N = 0;
  for (i = 0; i < N_in; i++) {
    if (isnan(mag_in[i])) continue;
    if (weighted && (isnan(err_in[i]) || !(err_in[i] > 0.0))) continue;
    if (usemask && maskvar != NULL) {
      if (!(EvaluateVariable_Double(lc_name_num, lcnum, i, maskvar) > VARTOOLS_MASK_TINY))
        continue;
    }
    mag[N] = mag_in[i];
    err[N] = weighted ? err_in[i] : 1.0;
    N++;
  }
  if (N < 2) {
    free(mag); free(err);
    return 0.0;
  }

  /* Mean and variance. */
  if (weighted) {
    for (i = 0; i < N; i++) {
      long double w = 1.0L / ((long double) err[i] * (long double) err[i]);
      sumy += w * (long double) mag[i];
      sumw += w;
    }
    ybar = sumy / sumw;
    for (i = 0; i < N; i++) {
      long double w = 1.0L / ((long double) err[i] * (long double) err[i]);
      long double dy = (long double) mag[i] - ybar;
      sumvar += w * dy * dy;
    }
    /* Keep s_sq as the UNNORMALIZED Σw_i (y_i - ȳ_w)²; the (2N/(N-1))
     * factor at the end restores E[eta_w] ≈ 2 for any sigma distribution. */
    s_sq = (double) sumvar;
  } else {
    for (i = 0; i < N; i++) sumy += (long double) mag[i];
    ybar = sumy / (long double) N;
    for (i = 0; i < N; i++) {
      long double dy = (long double) mag[i] - ybar;
      sumvar += dy * dy;
    }
    s_sq = (double) (sumvar / (long double) N);
  }

  /* Mean-square successive difference. */
  if (weighted) {
    for (i = 0; i < N - 1; i++) {
      long double sig_pair_sq = (long double) err[i] * err[i]
                              + (long double) err[i+1] * err[i+1];
      long double w_pair = 1.0L / sig_pair_sq;
      long double dy = (long double) mag[i+1] - (long double) mag[i];
      sumdy += w_pair * dy * dy;
    }
    /* Unnormalized Σw_pair_i (Δy)²; see s_sq comment above. */
    delta_sq = (double) sumdy;
  } else {
    for (i = 0; i < N - 1; i++) {
      long double dy = (long double) mag[i+1] - (long double) mag[i];
      sumdy += dy * dy;
    }
    delta_sq = (double) (sumdy / (long double) (N - 1));
  }

  free(mag); free(err);

  if (!(s_sq > 0.0)) return 0.0;
  if (weighted) {
    return (2.0 * (double) N / (double) (N - 1)) * delta_sq / s_sq;
  }
  return delta_sq / s_sq;
}
