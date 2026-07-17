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

/* -structurefunction: ensemble (single-LC) structure function and optional
 * damped-random-walk model fit.
 *
 * Two estimators are supported (selected by the "estimator" keyword).
 * The labels "first-order" and "second-order" are avoided here as
 * descriptors of the two forms: those terms have opposite meanings in
 * the Rutman 1978 / Simonetti 1985 time-series literature (where they
 * label a polynomial-detrending differencing scheme) and in the AGN
 * variability literature (where they label the moment power of |dm|).
 * The two forms are distinguished here by which moment of |dm| is
 * averaged inside the bin.
 *
 *   SF_ESTIMATOR_SQUARED (default) -- squared-difference form with
 *   bin-averaged noise subtraction.  This is the D^(1)(tau) of
 *   Simonetti, Cordes & Heeschen 1985, ApJ, 296, 46 (their Eq. A1 with
 *   M=1; finite-sample estimator A8), with the noise correction of
 *   their Eq. A13 generalised from a single homoscedastic noise
 *   variance 2 sigma^2_delta-f to the heteroscedastic per-pair sum
 *   sigma_i^2 + sigma_j^2:
 *
 *      SF^2(dt) = < (m_j - m_i)^2 >_pairs in bin
 *               - < sigma_i^2 + sigma_j^2 >_pairs in bin,    SF = sqrt(SF^2)
 *
 *   Bins where the right-hand side is non-positive are noise-dominated
 *   and reported as NaN.  This form is unbiased in the variance at the
 *   bin-averaged level and pairs naturally with the DRW autocovariance
 *   derivation (Kelly et al. 2009).
 *
 *   SF_ESTIMATOR_MAD -- absolute-deviation form with per-pair noise
 *   subtraction (Hughes, Aller & Aller 1992, ApJ, 396, 469; canonical
 *   SDSS-quasar implementation per Schmidt et al. 2010, ApJ, 714, 1194,
 *   their Eq. 2):
 *
 *      V(dt) = < (sqrt(pi)/2) * |m_j - m_i| - sqrt(sigma_i^2 + sigma_j^2)
 *              >_pairs in bin
 *
 *   The sqrt(pi)/2 factor converts the mean absolute deviation of a
 *   Gaussian to its standard deviation.  This form is more robust to
 *   outliers and stays real-valued at the noise floor (per-pair
 *   subtraction in the radical form, biased low rather than NaN), at the
 *   cost of a slightly biased response to V approaching the noise floor
 *   (asymptotes V_estimator = sqrt(V + sigma_n^2) - sqrt(sigma_n^2),
 *   which is close to sqrt(V) only when V >> sigma_n^2).
 *
 * Both estimators recover the intrinsic per-lag standard deviation in
 * the well-determined regime and match the DRW analytic SF
 *      SF_DRW(dt) = SF_inf * sqrt(1 - exp(-dt / tau))
 * for a Gaussian DRW.  The fit code is the same for both.
 *
 * The damped-random-walk (Ornstein-Uhlenbeck) model has analytic SF
 * (Kelly et al. 2009, ApJ, 698, 895; MacLeod et al. 2010, ApJ, 721, 1014):
 *
 *    SF_DRW(dt) = SF_inf * sqrt(1 - exp(-dt / tau))
 *
 * The fitted scalar reported here is sigma_long = SF_inf / sqrt(2), the
 * long-term standard deviation of the OU process (MacLeod 2010 sigma).
 * Tau is the damping time-scale in the same units as the input time axis.
 *
 * Procedure:
 *   1. Reject NaN magnitudes and points with sig <= 0; if "maskpoints" is
 *      in effect also reject points with maskvar <= 0.
 *   2. Resolve the per-LC lagrange / sigma0 / tau0 if user-specified, or
 *      derive a default lagrange of [min(dt), max(dt)] over the surviving
 *      pairs.
 *   3. Sweep all unordered pairs (i,j) with i < j, accumulate per bin:
 *        sum_dm2[k]  += (m_j - m_i)^2
 *        sum_sig2[k] += sig_i^2 + sig_j^2
 *        n_pairs[k]  += 1
 *      With require_sort = 1 the time array is ascending, so the inner
 *      loop breaks as soon as t[j] - t[i] exceeds the largest bin edge.
 *   4. Per bin: SF^2 = sum_dm2/n - sum_sig2/n; SF = sqrt(SF^2) if positive
 *      else NaN; sigma_SF = SF / sqrt(2 * N_eff) where N_eff = min(n_pairs,
 *      N_obs/2) bounds the effective number of independent pairs.  At
 *      short lags N_eff = n_pairs reduces to the standard Gaussian-rms
 *      error SF/sqrt(2 n_pairs).  At long lags every observation
 *      participates in many pairs across many bins; since the maximum
 *      matching on N_obs vertices has N_obs/2 disjoint pairs, no bin
 *      can hold more than N_obs/2 fully independent pair-samples, and
 *      sigma_SF saturates at ~SF/sqrt(N_obs).  This recovers the correct
 *      ~SF/sqrt(N_obs) precision floor at long lags, where the LC's
 *      finite N_obs caps any estimator's information content.
 *   5. (Optional) DRW fit: amoeba on (log SF_inf, log tau) minimising
 *      chi^2 = sum_bins [(SF_obs - SF_model) / sigma_SF]^2 over the bins
 *      with a finite SF.  Initial guess: SF_inf from the mean of the
 *      top-quartile of well-determined bins, tau from the lag where
 *      SF_obs first reaches SF_inf/2 (or the geometric mean of the lag
 *      range as a fallback).  User-supplied sigma0 / tau0 override.
 *   6. (Optional) write a 4-column aux file (dt_center, SF, sigma_SF,
 *      n_pairs) per light curve.
 *
 * Edge cases (DRW scalar columns reported as NaN, CONVERGED = 0):
 *   - fewer than 2 finite-SF bins;
 *   - amoeba returns non-finite parameters or exceeds max iterations
 *     without converging to within ftol.
 *
 * require_sort = 1 is set in the parser; the kernel assumes p->t[lc] is
 * sorted ascending.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SF_AMOEBA_FTOL    1.0e-7
#define SF_AMOEBA_MAXEVAL 5000

/* User-parameter block handed through amoeba to the chi^2 callback.
 * amoeba calls funk(p, N, Nparams, t, mag, sig, userparams); we stash
 * the SF bin arrays (dt_centre / SF / sigma_SF over the well-determined
 * bins) in t / mag / sig and don't use userparams. */

static double drw_sf_model(double SF_inf, double tau, double dt)
{
  return SF_inf * sqrt(1.0 - exp(-dt / tau));
}

/* amoeba is fitting log(SF_inf) and log(tau) so the search is
 * scale-invariant and the simplex doesn't have to span orders of
 * magnitude. */
/* amoeba calls funk(params, ndim, N, t, mag, sig, userparams), so the
 * second argument is the parameter count (always 2 here) and the third
 * is the data-point count (the well-determined bin count nfit). */
static double drw_sf_chi2(double *params, int Nparams, int Nbins_fit,
                          double *dt_arr, double *sf_arr, double *err_arr,
                          void *userparams)
{
  double SF_inf = exp(params[0]);
  double tau    = exp(params[1]);
  double chi2   = 0.0;
  double model, resid;
  int    i;
  (void) Nparams;
  if(!(tau > 0.0) || !(SF_inf > 0.0))
    return 1.0e300;
  for(i = 0; i < Nbins_fit; i++) {
    model = drw_sf_model(SF_inf, tau, dt_arr[i]);
    resid = (sf_arr[i] - model) / err_arr[i];
    chi2 += resid * resid;
  }
  return chi2;
}

/* Write NaN to every output cell for this light curve (when fitDRW is on). */
static void sf_set_drw_nan(_StructureFunction *sf, int lcnum)
{
  sf->sigma_long[lcnum] = sqrt(-1.0);
  sf->tau[lcnum]        = sqrt(-1.0);
  sf->chi2[lcnum]       = sqrt(-1.0);
  sf->dof[lcnum]        = 0;
  sf->converged[lcnum]  = 0;
}

/* Write NaN / 0 to every per-edge output cell for this LC (when
 * reportsfvalsintable is on).  report_dt = NaN by default; the kernel
 * overrides it with the bin centre when the edge falls in a valid bin. */
static void sf_set_report_nan(_StructureFunction *sf, int lcnum)
{
  int ke;
  for(ke = 0; ke < sf->n_report_edges; ke++) {
    sf->report_dt[lcnum][ke]       = sqrt(-1.0);
    sf->report_sf[lcnum][ke]       = sqrt(-1.0);
    sf->report_sigma_sf[lcnum][ke] = sqrt(-1.0);
    sf->report_npairs[lcnum][ke]   = 0;
  }
}

/* Build the array of bin edges (length Nbins+1) and centres (length
 * Nbins) for the resolved (lagmin, lagmax).  For BINS_EDGES the user-
 * supplied edges are copied verbatim and lagmin/lagmax are ignored.
 * Centres are geometric for log binning, arithmetic otherwise. */
static int sf_build_bins(_StructureFunction *sf,
                         double lagmin_eff, double lagmax_eff,
                         double *edges, double *centres)
{
  int    k;
  double r, step;
  switch(sf->bin_mode) {
    case SF_BINS_LOG:
      if(!(lagmin_eff > 0.0) || !(lagmax_eff > lagmin_eff))
        return -1;
      r = log(lagmax_eff / lagmin_eff) / (double) sf->Nbins;
      for(k = 0; k <= sf->Nbins; k++)
        edges[k] = lagmin_eff * exp(r * (double) k);
      for(k = 0; k < sf->Nbins; k++)
        centres[k] = sqrt(edges[k] * edges[k+1]);
      break;
    case SF_BINS_LINEAR:
      if(!(lagmax_eff > lagmin_eff))
        return -1;
      step = (lagmax_eff - lagmin_eff) / (double) sf->Nbins;
      for(k = 0; k <= sf->Nbins; k++)
        edges[k] = lagmin_eff + step * (double) k;
      for(k = 0; k < sf->Nbins; k++)
        centres[k] = 0.5 * (edges[k] + edges[k+1]);
      break;
    case SF_BINS_EDGES:
      for(k = 0; k <= sf->Nbins; k++)
        edges[k] = sf->user_edges[k];
      for(k = 0; k < sf->Nbins; k++)
        centres[k] = 0.5 * (edges[k] + edges[k+1]);
      break;
    default:
      return -1;
  }
  return 0;
}

/* Find the bin index for lag dt.  Returns -1 if dt is outside the range.
 * Log and linear modes use O(1) arithmetic; explicit-edges uses bsearch.
 * The right edge of the last bin is inclusive; all other right edges are
 * exclusive (i.e. bin k covers [edges[k], edges[k+1]), and bin Nbins-1
 * covers [edges[Nbins-1], edges[Nbins]]). */
static int sf_bin_index(_StructureFunction *sf, double dt,
                        const double *edges, double lagmin_eff,
                        double inv_logstep, double inv_linstep)
{
  int    k, lo, hi, mid;
  if(dt < edges[0]) return -1;
  if(dt > edges[sf->Nbins]) return -1;
  switch(sf->bin_mode) {
    case SF_BINS_LOG:
      if(dt <= 0.0) return -1;
      k = (int) floor(log(dt / lagmin_eff) * inv_logstep);
      if(k < 0) k = 0;
      if(k >= sf->Nbins) k = sf->Nbins - 1;
      return k;
    case SF_BINS_LINEAR:
      k = (int) floor((dt - lagmin_eff) * inv_linstep);
      if(k < 0) k = 0;
      if(k >= sf->Nbins) k = sf->Nbins - 1;
      return k;
    case SF_BINS_EDGES:
      /* Binary search: locate the bin k with edges[k] <= dt < edges[k+1],
       * extending the last bin's right edge to be inclusive. */
      lo = 0;
      hi = sf->Nbins - 1;
      while(lo < hi) {
        mid = (lo + hi + 1) >> 1;
        if(edges[mid] <= dt) lo = mid;
        else hi = mid - 1;
      }
      return lo;
    default:
      return -1;
  }
}

/* Pick the DRW initial guess for (SF_inf, tau).  sigma0 / tau0 override
 * if user-supplied; otherwise SF_inf is the mean of the SF in the upper
 * quartile of well-determined bins, and tau is the centre lag at which
 * the SF first reaches SF_inf/2 (linearly interpolated in log-lag).
 * Falls back to the geometric mean of the lag range. */
static void sf_drw_initial_guess(int Nfit, const double *dt_fit,
                                 const double *sf_fit,
                                 int have_sigma0, double sigma0_in,
                                 int have_tau0, double tau0_in,
                                 double *SF_inf_init, double *tau_init)
{
  int    i, lo, n_top;
  double sum;
  if(have_sigma0 && sigma0_in > 0.0) {
    /* User supplied sigma_long; convert back to SF_inf. */
    *SF_inf_init = sqrt(2.0) * sigma0_in;
  } else {
    /* Mean SF over the top quartile of well-determined bins. */
    lo = (3 * Nfit) / 4;
    if(lo >= Nfit) lo = Nfit - 1;
    n_top = Nfit - lo;
    sum = 0.0;
    for(i = lo; i < Nfit; i++) sum += sf_fit[i];
    *SF_inf_init = (n_top > 0 ? sum / (double) n_top : sf_fit[Nfit - 1]);
    if(!(*SF_inf_init > 0.0)) *SF_inf_init = sf_fit[Nfit - 1];
  }
  if(have_tau0 && tau0_in > 0.0) {
    *tau_init = tau0_in;
  } else {
    /* First lag where SF reaches half SF_inf. */
    double half = 0.5 * (*SF_inf_init);
    *tau_init = -1.0;
    for(i = 0; i < Nfit; i++) {
      if(sf_fit[i] >= half) {
        *tau_init = dt_fit[i];
        break;
      }
    }
    if(!(*tau_init > 0.0))
      *tau_init = sqrt(dt_fit[0] * dt_fit[Nfit - 1]);
  }
}

/* amoeba wrapper for the DRW chi^2; returns chi2_min in *chi2_out and
 * the fitted (SF_inf, tau) in *SF_inf_out / *tau_out.  Convergence is
 * declared when amoeba returns within SF_AMOEBA_MAXEVAL function
 * evaluations and the returned chi^2 + parameters are finite.  Returns
 * 1 on convergence, 0 otherwise. */
static int sf_fit_drw(int Nfit, double *dt_fit, double *sf_fit,
                      double *err_fit, double SF_inf_init, double tau_init,
                      double *SF_inf_out, double *tau_out, double *chi2_out)
{
  double **simplex = NULL, *yvals = NULL;
  int    *ia = NULL;
  int    nfunk = 0, ndim = 2, i, j, ok;
  double base[2], step[2];

  if((simplex = (double **) malloc(3 * sizeof(double *))) == NULL ||
     (yvals   = (double *)  malloc(3 * sizeof(double)))   == NULL ||
     (ia      = (int *)     malloc(2 * sizeof(int)))      == NULL)
    vt_error(ERR_MEMALLOC);
  for(i = 0; i < 3; i++)
    if((simplex[i] = (double *) malloc(2 * sizeof(double))) == NULL)
      vt_error(ERR_MEMALLOC);

  base[0] = log(SF_inf_init);
  base[1] = log(tau_init);
  step[0] = log(2.0);   /* simplex spans a factor of 2 in each direction */
  step[1] = log(2.0);
  ia[0] = 1; ia[1] = 1;

  /* Build a 3-vertex simplex around the initial guess. */
  simplex[0][0] = base[0];          simplex[0][1] = base[1];
  simplex[1][0] = base[0] + step[0]; simplex[1][1] = base[1];
  simplex[2][0] = base[0];          simplex[2][1] = base[1] + step[1];
  for(i = 0; i < 3; i++)
    yvals[i] = drw_sf_chi2(simplex[i], ndim, Nfit,
                           dt_fit, sf_fit, err_fit, NULL);

  /* amoeba returns 0 when the simplex converges within ftol, 1 when the
   * max-eval budget is exhausted before convergence. */
  ok = (amoeba(simplex, yvals, ia, ndim, SF_AMOEBA_FTOL, drw_sf_chi2,
               &nfunk, SF_AMOEBA_MAXEVAL, Nfit,
               dt_fit, sf_fit, err_fit, NULL) == 0);

  /* Best vertex is the lowest-y one. */
  j = 0;
  for(i = 1; i < 3; i++)
    if(yvals[i] < yvals[j]) j = i;
  *SF_inf_out = exp(simplex[j][0]);
  *tau_out    = exp(simplex[j][1]);
  *chi2_out   = yvals[j];

  for(i = 0; i < 3; i++) free(simplex[i]);
  free(simplex); free(yvals); free(ia);

  if(!ok) return 0;
  if(!isfinite(*chi2_out)) return 0;
  if(!(*SF_inf_out > 0.0) || !(*tau_out > 0.0)) return 0;
  if(!isfinite(*SF_inf_out) || !isfinite(*tau_out)) return 0;
  return 1;
}

void RunStructureFunctionCommand(ProgramData *p, _StructureFunction *sf,
                                 int lcnum, int lc_name_num, char *outname)
{
  int    N_in, N_filt, i, j, k, nfit;
  int    bin_mode = sf->bin_mode;
  int    Nbins    = sf->Nbins;
  double lagmin_user = 0.0, lagmax_user = 0.0;
  double lagmin_eff, lagmax_eff;
  double inv_logstep = 0.0, inv_linstep = 0.0;
  double *t_filt = NULL, *m_filt = NULL, *s2_filt = NULL;
  double *edges = NULL, *centres = NULL;
  double *sum_dm2 = NULL, *sum_sig2 = NULL, *sum_mad = NULL;
  long   *n_pairs = NULL;
  double *SF_curve = NULL, *sigma_SF = NULL;
  double dt, dm, sf2, mean_dm2, mean_sig2;
  const double sqrt_pi_half = 1.2533141373155003;  /* sqrt(pi/2) */
  double sigma0_in = 0.0, tau0_in = 0.0;
  int    have_sigma0 = sf->have_sigma0, have_tau0 = sf->have_tau0;
  double SF_inf_init, tau_init, SF_inf_fit, tau_fit, chi2_fit;
  int    converged = 0;
  FILE   *fp;

  N_in = p->NJD[lcnum];

  /* If this LC produces a DRW scalar set, default to NaN; populated on
   * success below. */
  if(sf->do_fit_drw)
    sf_set_drw_nan(sf, lcnum);
  if(sf->do_report_in_table)
    sf_set_report_nan(sf, lcnum);

  /* Resolve per-LC parameter values. */
  if(sf->have_lagrange) {
    lagmin_user = VT_EVAL_DOUBLE(sf, lagmin, lc_name_num, lcnum);
    lagmax_user = VT_EVAL_DOUBLE(sf, lagmax, lc_name_num, lcnum);
    if(!(lagmin_user > 0.0) || !(lagmax_user > lagmin_user))
      return;   /* nothing to do; aux file unwritten, DRW already NaN */
  }
  if(sf->do_fit_drw) {
    if(have_sigma0) sigma0_in = VT_EVAL_DOUBLE(sf, sigma0, lc_name_num, lcnum);
    if(have_tau0)   tau0_in   = VT_EVAL_DOUBLE(sf, tau0,   lc_name_num, lcnum);
    if(have_sigma0 && !(sigma0_in > 0.0)) have_sigma0 = 0;
    if(have_tau0   && !(tau0_in   > 0.0)) have_tau0   = 0;
  }

  if(N_in < 2) return;

  if((t_filt  = (double *) malloc(N_in * sizeof(double))) == NULL ||
     (m_filt  = (double *) malloc(N_in * sizeof(double))) == NULL ||
     (s2_filt = (double *) malloc(N_in * sizeof(double))) == NULL)
    vt_error(ERR_MEMALLOC);

  /* Filter NaN mags / nonpositive errors / masked points.  Errors are
   * stored squared since every pair uses sig_i^2 + sig_j^2. */
  N_filt = 0;
  for(i = 0; i < N_in; i++) {
    if(isnan(p->mag[lcnum][i])) continue;
    if(!(p->sig[lcnum][i] > 0.0)) continue;
    if(sf->usemask && sf->maskvar != NULL) {
      if(!(EvaluateVariable_Double(lc_name_num, lcnum, i,
                                   sf->maskvar) > VARTOOLS_MASK_TINY))
        continue;
    }
    t_filt[N_filt]  = p->t[lcnum][i];
    m_filt[N_filt]  = p->mag[lcnum][i];
    s2_filt[N_filt] = p->sig[lcnum][i] * p->sig[lcnum][i];
    N_filt++;
  }

  if(N_filt < 2) {
    free(t_filt); free(m_filt); free(s2_filt);
    return;
  }

  /* Effective lag range. */
  if(sf->have_lagrange) {
    lagmin_eff = lagmin_user;
    lagmax_eff = lagmax_user;
  } else if(bin_mode == SF_BINS_EDGES) {
    /* Edges binning ignores lagrange anyway. */
    lagmin_eff = sf->user_edges[0];
    lagmax_eff = sf->user_edges[Nbins];
  } else {
    /* Auto: [min(dt), max(dt)] over the surviving points (sorted in
     * time, so min(dt) is the smallest consecutive spacing and
     * max(dt) is the full baseline). */
    lagmax_eff = t_filt[N_filt - 1] - t_filt[0];
    lagmin_eff = lagmax_eff;
    for(i = 1; i < N_filt; i++) {
      dt = t_filt[i] - t_filt[i-1];
      if(dt > 0.0 && dt < lagmin_eff) lagmin_eff = dt;
    }
    if(!(lagmin_eff > 0.0) || !(lagmax_eff > lagmin_eff)) {
      free(t_filt); free(m_filt); free(s2_filt);
      return;
    }
  }

  if((edges    = (double *) malloc((Nbins + 1) * sizeof(double))) == NULL ||
     (centres  = (double *) malloc(Nbins * sizeof(double)))       == NULL ||
     (sum_dm2  = (double *) malloc(Nbins * sizeof(double)))       == NULL ||
     (sum_sig2 = (double *) malloc(Nbins * sizeof(double)))       == NULL ||
     (sum_mad  = (double *) malloc(Nbins * sizeof(double)))       == NULL ||
     (n_pairs  = (long *)   malloc(Nbins * sizeof(long)))         == NULL ||
     (SF_curve = (double *) malloc(Nbins * sizeof(double)))       == NULL ||
     (sigma_SF = (double *) malloc(Nbins * sizeof(double)))       == NULL)
    vt_error(ERR_MEMALLOC);

  if(sf_build_bins(sf, lagmin_eff, lagmax_eff, edges, centres) != 0) {
    free(t_filt); free(m_filt); free(s2_filt);
    free(edges); free(centres); free(sum_dm2); free(sum_sig2); free(sum_mad);
    free(n_pairs); free(SF_curve); free(sigma_SF);
    return;
  }

  if(bin_mode == SF_BINS_LOG)
    inv_logstep = (double) Nbins / log(lagmax_eff / lagmin_eff);
  else if(bin_mode == SF_BINS_LINEAR)
    inv_linstep = (double) Nbins / (lagmax_eff - lagmin_eff);

  for(k = 0; k < Nbins; k++) {
    sum_dm2[k]  = 0.0;
    sum_sig2[k] = 0.0;
    sum_mad[k]  = 0.0;
    n_pairs[k]  = 0;
  }

  /* All-pairs sweep.  require_sort = 1 in the parser so t is ascending;
   * the inner loop breaks as soon as t[j] - t[i] exceeds the largest
   * bin's right edge.  For the SQUARED estimator we accumulate sum_dm2
   * and sum_sig2 separately; for the MAD estimator we accumulate
   * (sqrt(pi)/2)*|dm| - sqrt(sig_i^2+sig_j^2) per pair, since Schmidt
   * 2010 Eq. 2 subtracts noise inside the angle brackets. */
  if(sf->estimator == SF_ESTIMATOR_MAD) {
    for(i = 0; i < N_filt - 1; i++) {
      for(j = i + 1; j < N_filt; j++) {
        dt = t_filt[j] - t_filt[i];
        if(dt > edges[Nbins]) break;
        k = sf_bin_index(sf, dt, edges, lagmin_eff,
                         inv_logstep, inv_linstep);
        if(k < 0) continue;
        dm = m_filt[j] - m_filt[i];
        sum_mad[k] += sqrt_pi_half * fabs(dm)
                      - sqrt(s2_filt[i] + s2_filt[j]);
        n_pairs[k] += 1;
      }
    }
  } else {
    for(i = 0; i < N_filt - 1; i++) {
      for(j = i + 1; j < N_filt; j++) {
        dt = t_filt[j] - t_filt[i];
        if(dt > edges[Nbins]) break;
        k = sf_bin_index(sf, dt, edges, lagmin_eff,
                         inv_logstep, inv_linstep);
        if(k < 0) continue;
        dm = m_filt[j] - m_filt[i];
        sum_dm2[k]  += dm * dm;
        sum_sig2[k] += s2_filt[i] + s2_filt[j];
        n_pairs[k]  += 1;
      }
    }
  }

  /* Per-bin SF and sigma_SF.  N_eff = min(n_pairs, N_filt/2) caps the
   * effective independent-pair count: with N_filt observations the
   * maximum matching contains N_filt/2 disjoint pairs, so no bin can
   * carry more than N_filt/2 fully independent pair-samples.  This
   * gives the right ~SF/sqrt(N_obs) noise floor at long lags.
   *
   * The per-sample variance-to-sigma conversion is estimator-specific:
   *   - SQUARED: for a Gaussian (m_j-m_i)^2 has variance ~2 V^2 where
   *     V = E[(m_j-m_i)^2]; standard delta-method gives
   *     sigma_SF = SF / sqrt(2 N_eff).
   *   - MAD: (sqrt(pi)/2)*|m_j-m_i| has per-sample variance
   *     ((pi/2)-1) * V (where V = E[(m_j-m_i)^2]), so
   *     sigma_V = V * sqrt((pi/2)-1) / sqrt(N_eff)
   *            = V / sqrt(N_eff/((pi/2)-1)). */
  {
    long n_eff_cap = (long) N_filt / 2;
    const double mad_norm = 1.0 / (M_PI / 2.0 - 1.0);  /* (pi/2 - 1)^-1 */
    if(n_eff_cap < 1) n_eff_cap = 1;
    for(k = 0; k < Nbins; k++) {
      long n_eff;
      if(n_pairs[k] <= 0) {
        SF_curve[k] = sqrt(-1.0);
        sigma_SF[k] = sqrt(-1.0);
        continue;
      }
      n_eff = n_pairs[k] < n_eff_cap ? n_pairs[k] : n_eff_cap;
      if(sf->estimator == SF_ESTIMATOR_MAD) {
        double V_estimator = sum_mad[k] / (double) n_pairs[k];
        if(V_estimator > 0.0) {
          SF_curve[k] = V_estimator;
          sigma_SF[k] = V_estimator / sqrt((double) n_eff * mad_norm);
        } else {
          /* Per-pair noise-subtracted MAD has gone non-positive: the
           * bin is noise-dominated and the unbiased intrinsic
           * standard deviation is effectively zero / below the
           * detection floor.  Reported as NaN so the DRW fit ignores
           * the bin (consistent with the squared-form behaviour). */
          SF_curve[k] = sqrt(-1.0);
          sigma_SF[k] = sqrt(-1.0);
        }
      } else {
        mean_dm2  = sum_dm2[k]  / (double) n_pairs[k];
        mean_sig2 = sum_sig2[k] / (double) n_pairs[k];
        sf2 = mean_dm2 - mean_sig2;
        if(sf2 > 0.0) {
          SF_curve[k] = sqrt(sf2);
          sigma_SF[k] = SF_curve[k] / sqrt(2.0 * (double) n_eff);
        } else {
          SF_curve[k] = sqrt(-1.0);
          sigma_SF[k] = sqrt(-1.0);
        }
      }
    }
  }

  /* Optional in-table per-edge scalar output: for each user-supplied
   * lag, locate the bin containing it and copy its (centre, SF,
   * sigma_SF, n_pairs) into the output columns.  Out-of-range edges
   * leave the per-LC defaults (NaN / 0) untouched. */
  if(sf->do_report_in_table) {
    int ke;
    for(ke = 0; ke < sf->n_report_edges; ke++) {
      int kb = sf_bin_index(sf, sf->report_edges[ke], edges, lagmin_eff,
                            inv_logstep, inv_linstep);
      if(kb < 0) continue;
      sf->report_dt[lcnum][ke]       = centres[kb];
      sf->report_sf[lcnum][ke]       = SF_curve[kb];
      sf->report_sigma_sf[lcnum][ke] = sigma_SF[kb];
      sf->report_npairs[lcnum][ke]   = (int) n_pairs[kb];
    }
  }

  /* Optional aux-file write (always 4 columns, NaN rows preserved). */
  if(sf->do_save && outname != NULL && outname[0] != '\0') {
    if((fp = fopen(outname, "w")) == NULL) {
      fprintf(stderr,
        "Error: cannot open '%s' for writing the structure function\n",
        outname);
      exit(3);
    }
    fprintf(fp, "# dt_center  SF  sigma_SF  n_pairs\n");
    for(k = 0; k < Nbins; k++)
      fprintf(fp, "%.10g %.10g %.10g %ld\n",
              centres[k], SF_curve[k], sigma_SF[k], n_pairs[k]);
    fclose(fp);
  }

  /* Optional DRW fit on the well-determined bins. */
  if(sf->do_fit_drw) {
    double *dt_fit = NULL, *sf_fit = NULL, *err_fit = NULL;
    if((dt_fit  = (double *) malloc(Nbins * sizeof(double))) == NULL ||
       (sf_fit  = (double *) malloc(Nbins * sizeof(double))) == NULL ||
       (err_fit = (double *) malloc(Nbins * sizeof(double))) == NULL)
      vt_error(ERR_MEMALLOC);
    nfit = 0;
    for(k = 0; k < Nbins; k++) {
      if(isnan(SF_curve[k]) || isnan(sigma_SF[k])) continue;
      if(!(sigma_SF[k] > 0.0)) continue;
      dt_fit[nfit]  = centres[k];
      sf_fit[nfit]  = SF_curve[k];
      err_fit[nfit] = sigma_SF[k];
      nfit++;
    }
    if(nfit >= 3) {
      sf_drw_initial_guess(nfit, dt_fit, sf_fit,
                           have_sigma0, sigma0_in,
                           have_tau0,   tau0_in,
                           &SF_inf_init, &tau_init);
      if(SF_inf_init > 0.0 && tau_init > 0.0) {
        converged = sf_fit_drw(nfit, dt_fit, sf_fit, err_fit,
                               SF_inf_init, tau_init,
                               &SF_inf_fit, &tau_fit, &chi2_fit);
        if(converged) {
          sf->sigma_long[lcnum] = SF_inf_fit / sqrt(2.0);
          sf->tau[lcnum]        = tau_fit;
          sf->chi2[lcnum]       = chi2_fit;
          sf->dof[lcnum]        = nfit - 2;
          sf->converged[lcnum]  = 1;
        }
      }
    }
    free(dt_fit); free(sf_fit); free(err_fit);
  }

  free(t_filt); free(m_filt); free(s2_filt);
  free(edges); free(centres); free(sum_dm2); free(sum_sig2); free(sum_mad);
  free(n_pairs); free(SF_curve); free(sigma_SF);
}
