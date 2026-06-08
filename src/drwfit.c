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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* -drwfit: direct maximum-likelihood fit of a damped-random-walk (DRW /
 * CAR(1) / Ornstein-Uhlenbeck) model to a single light curve via the
 * Kelly, Bechtold & Siemiginowska 2009 (ApJ, 698, 895) state-space
 * recursion (their Equations 6-13).
 *
 * Model (Kelly Eq. 1):
 *
 *    dX(t) = -(1/tau) X(t) dt + sigma_K sqrt(dt) eps(t) + b dt
 *
 * with mean E[X] = b * tau = mu and long-term variance Var[X] = tau *
 * sigma_K^2 / 2.  We parameterise the recursion in (sigma_long, tau, mu)
 * with
 *
 *    sigma_long^2 = tau * sigma_K^2 / 2,
 *
 * i.e. sigma_long is the long-term magnitude standard deviation
 * (MacLeod et al. 2010, ApJ, 721, 1014 sigma; their SF_inf = sqrt(2) *
 * sigma_long).  This is the same convention as the optional fitDRW
 * output of -structurefunction, so the two commands are directly
 * comparable on the same light curve.  To convert to Kelly's SDE
 * driving-noise amplitude sigma_K (units mag * day^(-1/2)), use
 *
 *    sigma_K = sigma_long * sqrt(2 / tau)
 *
 * (MacLeod 2010 footnote 8).
 *
 * State-space recursion on time-sorted observations (x_i +/- sig_meas_i
 * at times t_i, i = 0..n-1):
 *
 *   x*_i      = x_i - mu                                   (Eq. 7)
 *   x_hat_0   = 0                                          (Eq. 8)
 *   Omega_0   = sigma_long^2                               (Eq. 9)
 *   a_i       = exp(-(t_i - t_{i-1}) / tau)                (Eq. 12)
 *   K_i       = Omega_{i-1} / (Omega_{i-1} + sig_meas_{i-1}^2)
 *   x_hat_i   = a_i * (x_hat_{i-1} + K_i * (x*_{i-1} - x_hat_{i-1}))   (Eq. 10)
 *   Omega_i   = sigma_long^2 * (1 - a_i^2)
 *               + a_i^2 * Omega_{i-1} * (1 - K_i)          (Eq. 11)
 *   -2 ln L   = sum_i [ ln(2 pi (Omega_i + sig_meas_i^2))
 *                       + (x*_i - x_hat_i)^2 / (Omega_i + sig_meas_i^2) ]
 *                                                          (Eq. 6, logged)
 *
 * Each likelihood evaluation is O(n) -- no matrix inversion needed.  The
 * recursion is exact for the irregularly-sampled CAR(1) process; it is
 * the same Kalman-filter form derived independently by Rybicki & Press
 * 1992 for the analogous Ornstein-Uhlenbeck process.
 *
 * The minimiser is Nelder-Mead (amoeba) over (log sigma_long, log tau)
 * when mu is held fixed or pre-subtracted, or over (log sigma_long,
 * log tau, mu) when mu is fit jointly.  Log parameterisation in
 * sigma_long and tau keeps the search scale-invariant and avoids
 * accidental excursions to negative values.
 *
 * Numerical care: 1 - a^2 is computed via expm1 to preserve precision
 * when consecutive Delta t is much smaller than tau (a known
 * floating-point cancellation pitfall when the sample spacing is small
 * compared to the damping time-scale).
 * Both terms of the Omega recursion are non-negative and the equilibrium
 * Omega = sigma_long^2 is a fixed point, so Omega stays in
 * [0, sigma_long^2] for all i and never overflows.
 *
 * Likelihood-ratio quality indicators (MacLeod et al. 2010 Sec. 3.1):
 *
 *   Delta L_noise = ln L_best - ln L_{sigma_long -> 0}
 *   Delta L_inf   = ln L_best - ln L_{tau         -> inf}
 *
 * Both are emitted unconditionally as scalar output columns:
 *   - L_{sigma_long -> 0} has a closed form: the Gaussian-noise-only
 *     likelihood at the analytic best-fit mu_weighted = weighted mean.
 *   - L_{tau -> inf} is evaluated by running the recursion at
 *     tau_eff = 1e6 * T_baseline with (sigma_long_best, mu_best).
 *     Holding sigma_long fixed as tau -> inf keeps the equilibrium
 *     variance finite (in contrast to holding sigma_K fixed, which
 *     gives an unbounded Brownian-motion limit), which is the regime
 *     these indicators are intended to flag.
 *
 * require_sort = 1 is set in the parser; the kernel assumes p->t[lcnum]
 * is sorted ascending.  All NaN magnitudes and points with sig_meas <= 0
 * are filtered out before the recursion; if maskpoints is in effect we
 * additionally reject points with maskvar <= 0.
 *
 * Edge cases (all scalars set to NaN, CONVERGED = 0):
 *   - fewer than 3 surviving points after the filter;
 *   - the entire surviving span is degenerate (T_baseline == 0);
 *   - amoeba returns non-finite parameters or exits max-eval without
 *     converging within ftol;
 *   - fitted sigma_long or tau is non-positive or non-finite.
 */

#define DRWFIT_AMOEBA_FTOL    1.0e-7
#define DRWFIT_AMOEBA_MAXEVAL 5000

/* Multiplier on T_baseline for the tau -> inf limit-likelihood.  At 1e6
 * the largest exp(-Delta t / tau_eff) over the LC is 1 - 1e-6, well below
 * the ftol threshold the amoeba converges to in the best-fit search, so
 * the kernel's response at tau_eff is numerically indistinguishable from
 * the actual tau -> inf limit. */
#define DRWFIT_TAU_INFINITY_MULTIPLIER 1.0e6

/* DRWFIT_MEAN_FIT / _FIX / _SUBTRACT defined in commands.h so the parser
 * in parsecommandline.c can also reference them. */

/* Userparams block passed through amoeba to the negative-log-likelihood
 * callback.  amoeba's funk signature is
 *   double funk(double *params, int Nparams, int Npoints,
 *               double *t, double *mag, double *sig, void *user);
 * We use t/mag/sig as the filtered time / magnitude / squared-sig arrays
 * and stash mean-handling mode + the fixed mu (when applicable) in
 * userparams. */
typedef struct {
  int    mean_mode;   /* DRWFIT_MEAN_FIT | _FIX | _SUBTRACT */
  double mu_fixed;    /* used when mean_mode != DRWFIT_MEAN_FIT */
} drw_userparams;

/* -2 ln L for the CAR(1) state-space recursion at (sigma_long, tau, mu).
 * Returns -2 ln L on success or +1e300 on numerical failure.
 *
 * Inputs:
 *   t          : observation times, ascending, length n
 *   x          : observed magnitudes, length n
 *   sig2       : squared 1-sigma photometric errors, length n (> 0)
 *   n          : number of observations (n >= 2)
 *   sigma_long : long-term magnitude standard deviation, > 0
 *   tau        : damping time-scale, time-axis units, > 0
 *   mu         : long-term mean of the process
 *
 * Optional outputs (each may be NULL; pass non-NULL only when the
 * forward filter trajectory is wanted, e.g. before the RTS smoother):
 *   x_hat_minus_out[i]  : prior mean at i (forecast, before observing x_i)
 *   Omega_minus_out[i]  : prior variance at i
 *   x_hat_plus_out[i]   : posterior mean at i (after Bayesian update with x_i)
 *   Omega_plus_out[i]   : posterior variance at i
 *
 * The recorded x_hat_*_out arrays are in CENTERED coordinates
 * (x* = x - mu); add mu when converting to a magnitude-scale prediction. */
static double kelly_neg2lnL(const double *t, const double *x,
                            const double *sig2, int n,
                            double sigma_long, double tau, double mu,
                            double *x_hat_minus_out, double *Omega_minus_out,
                            double *x_hat_plus_out,  double *Omega_plus_out)
{
  double Omega, x_hat, x_star_prev, K, a, one_minus_a2, var, resid, neg2lnL;
  double x_hat_plus_prev, Omega_plus_prev;
  double dt_over_tau, sl2;
  int    record_minus, record_plus;
  int    i;
  if(!(sigma_long > 0.0) || !(tau > 0.0)) return 1.0e300;
  record_minus = (x_hat_minus_out != NULL || Omega_minus_out != NULL);
  record_plus  = (x_hat_plus_out  != NULL || Omega_plus_out  != NULL);
  sl2     = sigma_long * sigma_long;
  Omega   = sl2;                         /* Eq. 9: Omega_0 = sigma_long^2 */
  x_hat   = 0.0;                          /* Eq. 8: x_hat_0 = 0 */
  neg2lnL = 0.0;
  for(i = 0; i < n; i++) {
    if(i > 0) {
      /* Update from i-1 -> i (Eqs. 10, 11, 12).  Use expm1 so 1 - a^2 is
       * accurate when Delta t / tau is small (a known floating-point
       * cancellation pitfall when the sample spacing is small compared
       * to the damping time-scale).  The update is logically two steps:
       *   posterior at i-1 = prior at i-1 + Bayesian update using x_{i-1}
       *   prior at i       = dynamics(posterior at i-1, dt_{i-1 -> i})
       * which lets the optional record_plus path capture the posterior. */
      dt_over_tau  = (t[i] - t[i-1]) / tau;
      one_minus_a2 = -expm1(-2.0 * dt_over_tau);   /* = 1 - exp(-2 dt/tau) = 1 - a^2 */
      a            = exp(-dt_over_tau);             /* needed for x_hat update */
      K            = Omega / (Omega + sig2[i-1]);
      x_star_prev  = x[i-1] - mu;
      x_hat_plus_prev = x_hat + K * (x_star_prev - x_hat);
      Omega_plus_prev = Omega * (1.0 - K);
      if(record_plus) {
        if(x_hat_plus_out) x_hat_plus_out[i-1] = x_hat_plus_prev;
        if(Omega_plus_out) Omega_plus_out[i-1] = Omega_plus_prev;
      }
      x_hat = a * x_hat_plus_prev;
      Omega = sl2 * one_minus_a2 + (1.0 - one_minus_a2) * Omega_plus_prev;
    }
    if(record_minus) {
      if(x_hat_minus_out) x_hat_minus_out[i] = x_hat;
      if(Omega_minus_out) Omega_minus_out[i] = Omega;
    }
    /* Likelihood contribution at observation i.  Eq. 6 in log form. */
    var = Omega + sig2[i];
    if(!(var > 0.0)) return 1.0e300;
    resid    = (x[i] - mu) - x_hat;
    neg2lnL += log(2.0 * M_PI * var) + resid * resid / var;
  }
  /* Posterior at the final step (not advanced into an iter n) -- only
   * the RTS smoother's initial conditions need this, so we only compute
   * it when the caller asked for the plus arrays. */
  if(record_plus && n > 0) {
    K = Omega / (Omega + sig2[n-1]);
    if(x_hat_plus_out) x_hat_plus_out[n-1] = x_hat + K * ((x[n-1] - mu) - x_hat);
    if(Omega_plus_out) Omega_plus_out[n-1] = Omega * (1.0 - K);
  }
  if(!isfinite(neg2lnL)) return 1.0e300;
  return neg2lnL;
}

/* Rauch-Tung-Striebel backward smoother for the CAR(1) Kalman filter.
 * Given the forward filter's prior/posterior moments at each step, fills
 * (x_smooth, Omega_smooth) with the smoothed posterior moments
 * conditional on the full observation sequence.
 *
 * Smoother gain G_i = a_{i+1} * Omega_plus[i] / Omega_minus[i+1].
 * Recurrence (running backward from i = n-1 to 0):
 *   x_smooth_i   = x_hat_plus_i + G_i * (x_smooth_{i+1} - x_hat_minus_{i+1})
 *   Omega_smooth_i = Omega_plus_i + G_i^2 * (Omega_smooth_{i+1} - Omega_minus_{i+1})
 *
 * Initial conditions at i = n-1: x_smooth = x_hat_plus, Omega_smooth = Omega_plus.
 * O(n) work, no allocation.  All x_*_out arrays in centered coordinates. */
static void kelly_rts_smoother(const double *t, int n, double tau,
                               const double *x_hat_minus,
                               const double *Omega_minus,
                               const double *x_hat_plus,
                               const double *Omega_plus,
                               double *x_smooth, double *Omega_smooth)
{
  double dt, a, G;
  int    i;
  if(n <= 0) return;
  x_smooth[n-1]     = x_hat_plus[n-1];
  Omega_smooth[n-1] = Omega_plus[n-1];
  for(i = n - 2; i >= 0; i--) {
    dt = t[i+1] - t[i];
    a  = exp(-dt / tau);
    /* Smoother gain.  Omega_minus[i+1] is strictly positive whenever
     * sigma_long > 0 (the prior variance after a forward step has the
     * sigma_long^2 (1 - a^2) term).  The kernel guards sigma_long > 0
     * before calling; defensive check just in case. */
    if(!(Omega_minus[i+1] > 0.0)) {
      G = 0.0;
    } else {
      G = a * Omega_plus[i] / Omega_minus[i+1];
    }
    x_smooth[i]     = x_hat_plus[i] + G * (x_smooth[i+1] - x_hat_minus[i+1]);
    Omega_smooth[i] = Omega_plus[i] + G * G * (Omega_smooth[i+1] - Omega_minus[i+1]);
  }
}

/* Closed-form -2 ln L_{sigma_long -> 0}: pure measurement noise.  Mu is
 * held at the analytic best-fit weighted mean for the noise-only model
 * (which is the argmax of L_noise(mu)).  Returns -2 ln L_noise. */
static double kelly_neg2lnL_noise_limit(const double *x, const double *sig2,
                                        int n)
{
  double w, wm, wsum = 0.0, wxsum = 0.0;
  double neg2lnL = 0.0;
  double resid;
  int    i;
  for(i = 0; i < n; i++) {
    w      = 1.0 / sig2[i];
    wsum  += w;
    wxsum += w * x[i];
  }
  if(!(wsum > 0.0)) return 1.0e300;
  wm = wxsum / wsum;
  for(i = 0; i < n; i++) {
    resid    = x[i] - wm;
    neg2lnL += log(2.0 * M_PI * sig2[i]) + resid * resid / sig2[i];
  }
  return neg2lnL;
}

/* -2 ln L in the tau -> inf limit: run the recursion at tau_eff =
 * T_baseline * DRWFIT_TAU_INFINITY_MULTIPLIER with the best-fit
 * (sigma_long, mu).  Holding sigma_long fixed as tau -> inf keeps the
 * equilibrium variance Omega = sigma_long^2 finite; the recursion
 * behaves as the slow-mixing limit, which is the regime MacLeod's cut
 * was designed to detect. */
static double kelly_neg2lnL_tau_infinity_limit(const double *t,
                                               const double *x,
                                               const double *sig2,
                                               int n,
                                               double sigma_long, double mu)
{
  double T_baseline, tau_eff;
  if(n < 2) return 1.0e300;
  T_baseline = t[n-1] - t[0];
  if(!(T_baseline > 0.0)) return 1.0e300;
  tau_eff = T_baseline * DRWFIT_TAU_INFINITY_MULTIPLIER;
  return kelly_neg2lnL(t, x, sig2, n, sigma_long, tau_eff, mu,
                       NULL, NULL, NULL, NULL);
}

/* amoeba callback.  params encode (log sigma_long, log tau) when
 * mean_mode != FIT, or (log sigma_long, log tau, mu) when mean_mode = FIT.
 * Returns -2 ln L, with 1e300 sentinel on numerical failure (non-finite
 * decoded parameters etc.). */
static double drw_amoeba_neg2lnL(double *params, int Nparams, int Npoints,
                                 double *t, double *x, double *sig2,
                                 void *userparams)
{
  drw_userparams *up = (drw_userparams *) userparams;
  double sigma_long, tau, mu;
  (void) Nparams;
  sigma_long = exp(params[0]);
  tau        = exp(params[1]);
  if(up->mean_mode == DRWFIT_MEAN_FIT) mu = params[2];
  else                                 mu = up->mu_fixed;
  if(!isfinite(sigma_long) || !isfinite(tau) || !isfinite(mu))
    return 1.0e300;
  return kelly_neg2lnL(t, x, sig2, Npoints, sigma_long, tau, mu,
                       NULL, NULL, NULL, NULL);
}

/* Set all DRWFit output scalars to NaN / 0 for this LC.  Called at entry
 * so any early-return path leaves the row well-defined. */
static void drwfit_set_nan(_DRWFit *drw, int lcnum)
{
  double nan_val = sqrt(-1.0);
  drw->sigma_long[lcnum] = nan_val;
  drw->tau[lcnum]        = nan_val;
  drw->mu[lcnum]         = nan_val;
  drw->lnL[lcnum]        = nan_val;
  drw->dlnL_noise[lcnum] = nan_val;
  drw->dlnL_inf[lcnum]   = nan_val;
  drw->converged[lcnum]  = 0;
}

/* Default initial guesses (only when not user-overridden).  Plain
 * heuristics; the amoeba does the real work but it needs a sane start.
 *
 *   sigma_long_init = sqrt(max(Var_LC - <sig_meas^2>, 0)),
 *                     floored at a small positive value;
 *   tau_init        = 0.1 * T_baseline  (the central decile of MacLeod's
 *                     empirical tau / T_LC distribution for S82 quasars);
 *   mu_init         = weighted mean of the LC. */
static void drwfit_initial_guess(int n, const double *t, const double *x,
                                 const double *sig2,
                                 int have_sigma0, double sigma0_in,
                                 int have_tau0,   double tau0_in,
                                 int have_mean0,  double mean0_in,
                                 double *sigma_long_init,
                                 double *tau_init,
                                 double *mu_init)
{
  double w, wsum = 0.0, wxsum = 0.0;
  double wm, sum_sq = 0.0, var_lc, mean_sig2 = 0.0;
  double T_baseline;
  int    i;

  /* Weighted mean. */
  for(i = 0; i < n; i++) {
    w      = 1.0 / sig2[i];
    wsum  += w;
    wxsum += w * x[i];
  }
  wm = (wsum > 0.0 ? wxsum / wsum : 0.0);

  /* Unweighted variance about wm minus mean photometric noise variance.
   * Floor avoids log(0) if the LC is noise-dominated. */
  for(i = 0; i < n; i++) {
    sum_sq    += (x[i] - wm) * (x[i] - wm);
    mean_sig2 += sig2[i];
  }
  var_lc    = sum_sq / (double) n;
  mean_sig2 = mean_sig2 / (double) n;
  var_lc   -= mean_sig2;
  if(!(var_lc > 0.0)) var_lc = 0.01 * mean_sig2 + 1.0e-12;

  T_baseline = t[n-1] - t[0];

  *sigma_long_init = have_sigma0 ? sigma0_in : sqrt(var_lc);
  *tau_init        = have_tau0   ? tau0_in   : 0.1 * T_baseline;
  *mu_init         = have_mean0  ? mean0_in  : wm;

  /* Sanity floor on positive-domain params. */
  if(!(*sigma_long_init > 0.0)) *sigma_long_init = 1.0e-6;
  if(!(*tau_init > 0.0))
    *tau_init = (T_baseline > 0.0 ? 0.1 * T_baseline : 1.0);
}

/* Wrap amoeba.  Returns 1 on convergence with a finite positive
 * (sigma_long, tau), 0 otherwise.  *sigma_long_out / *tau_out / *mu_out /
 * *neg2lnL_out are filled with the best-vertex values regardless. */
static int drwfit_run_amoeba(int n, double *t, double *x, double *sig2,
                             double sigma_long_init, double tau_init,
                             double mu_init,
                             int mean_mode, double mu_fixed,
                             double *sigma_long_out, double *tau_out,
                             double *mu_out, double *neg2lnL_out)
{
  drw_userparams up;
  int    ndim, nverts;
  int    *ia = NULL;
  double **simplex = NULL, *yvals = NULL;
  double base[3], step[3];
  int    nfunk = 0, i, j, ok;

  up.mean_mode = mean_mode;
  up.mu_fixed  = mu_fixed;

  ndim   = (mean_mode == DRWFIT_MEAN_FIT) ? 3 : 2;
  nverts = ndim + 1;

  if((simplex = (double **) malloc(nverts * sizeof(double *))) == NULL ||
     (yvals   = (double *)  malloc(nverts * sizeof(double)))   == NULL ||
     (ia      = (int *)     malloc(ndim   * sizeof(int)))      == NULL)
    vt_error(ERR_MEMALLOC);
  for(i = 0; i < nverts; i++)
    if((simplex[i] = (double *) malloc(ndim * sizeof(double))) == NULL)
      vt_error(ERR_MEMALLOC);

  base[0] = log(sigma_long_init);
  base[1] = log(tau_init);
  base[2] = mu_init;
  step[0] = log(2.0);     /* span factor of 2 around initial guess in log-sigma */
  step[1] = log(2.0);     /* span factor of 2 around initial guess in log-tau */
  step[2] = (sigma_long_init > 0.0 ? sigma_long_init : 1.0);
  for(i = 0; i < ndim; i++) ia[i] = 1;

  /* Vertex 0 = base; vertex k+1 = base + step[k] on dimension k. */
  for(j = 0; j < ndim; j++) simplex[0][j] = base[j];
  for(i = 1; i < nverts; i++) {
    for(j = 0; j < ndim; j++) simplex[i][j] = base[j];
    simplex[i][i-1] += step[i-1];
  }
  for(i = 0; i < nverts; i++)
    yvals[i] = drw_amoeba_neg2lnL(simplex[i], ndim, n,
                                  t, x, sig2, (void *) &up);

  /* amoeba returns 0 on convergence within ftol, non-zero on max-eval
   * exhaustion. */
  ok = (amoeba(simplex, yvals, ia, ndim, DRWFIT_AMOEBA_FTOL,
               drw_amoeba_neg2lnL, &nfunk, DRWFIT_AMOEBA_MAXEVAL,
               n, t, x, sig2, (void *) &up) == 0);

  /* Pick the lowest-y vertex. */
  j = 0;
  for(i = 1; i < nverts; i++)
    if(yvals[i] < yvals[j]) j = i;
  *sigma_long_out = exp(simplex[j][0]);
  *tau_out        = exp(simplex[j][1]);
  *mu_out         = (mean_mode == DRWFIT_MEAN_FIT) ? simplex[j][2] : mu_fixed;
  *neg2lnL_out    = yvals[j];

  for(i = 0; i < nverts; i++) free(simplex[i]);
  free(simplex); free(yvals); free(ia);

  if(!ok) return 0;
  if(!isfinite(*neg2lnL_out)) return 0;
  if(!(*sigma_long_out > 0.0) || !(*tau_out > 0.0)) return 0;
  if(!isfinite(*sigma_long_out) || !isfinite(*tau_out)) return 0;
  if(!isfinite(*mu_out)) return 0;
  return 1;
}

void RunDRWFitCommand(ProgramData *p, _DRWFit *drw, int lcnum,
                      int lc_name_num, char *outname)
{
  int    N_in, N_filt, i, j;
  double *t_filt = NULL, *m_filt = NULL, *s2_filt = NULL;
  int    *orig_idx = NULL;
  double sigma_long_init, tau_init, mu_init;
  double sigma_long_fit, tau_fit, mu_fit, neg2lnL_best;
  double neg2lnL_noise, neg2lnL_inf;
  double sigma0_in = 0.0, tau0_in = 0.0, mean0_in = 0.0, mu_fixed = 0.0;
  double wm_subtracted = 0.0;  /* mean=subtract: weighted mean removed */
  int    have_sigma0 = drw->have_sigma0;
  int    have_tau0   = drw->have_tau0;
  int    have_mean0  = drw->have_mean0;
  int    mean_mode   = drw->mean_mode;
  int    converged   = 0;
  int    need_forecast, need_smoothed;
  double *xhat_minus = NULL, *Omega_minus = NULL;
  double *xhat_plus  = NULL, *Omega_plus  = NULL;
  double *x_smooth   = NULL, *Omega_smooth = NULL;
  double *model_smoothed = NULL, *model_forecast = NULL;
  double  nan_val = sqrt(-1.0);

  /* Default-NaN every scalar output column for this LC; updated on success. */
  drwfit_set_nan(drw, lcnum);

  /* Default-NaN the modelvar at every original LC index.  Done early so a
   * non-converged fit leaves the variable well-defined for downstream
   * commands rather than holding stale data from a previous LC. */
  N_in = p->NJD[lcnum];
  if(drw->do_modelvar && drw->modelvar != NULL) {
    for(i = 0; i < N_in; i++)
      (*((double ***) drw->modelvar->dataptr))[lcnum][i] = nan_val;
  }

  if(N_in < 3) return;

  /* Resolve per-LC parameter values from var / expr sources. */
  if(mean_mode == DRWFIT_MEAN_FIX)
    mu_fixed = VT_EVAL_DOUBLE(drw, mu_fix, lc_name_num, lcnum);
  if(have_sigma0) sigma0_in = VT_EVAL_DOUBLE(drw, sigma0, lc_name_num, lcnum);
  if(have_tau0)   tau0_in   = VT_EVAL_DOUBLE(drw, tau0,   lc_name_num, lcnum);
  if(have_mean0)  mean0_in  = VT_EVAL_DOUBLE(drw, mean0,  lc_name_num, lcnum);

  if(have_sigma0 && !(sigma0_in > 0.0)) have_sigma0 = 0;
  if(have_tau0   && !(tau0_in   > 0.0)) have_tau0   = 0;

  if((t_filt   = (double *) malloc(N_in * sizeof(double))) == NULL ||
     (m_filt   = (double *) malloc(N_in * sizeof(double))) == NULL ||
     (s2_filt  = (double *) malloc(N_in * sizeof(double))) == NULL ||
     (orig_idx = (int *)    malloc(N_in * sizeof(int)))    == NULL)
    vt_error(ERR_MEMALLOC);

  /* Filter NaN mags / non-positive errors / masked points.  require_sort
   * = 1 in the parser, so t_filt[] inherits the ascending order.  We
   * track orig_idx[i] = original index of the i-th surviving point so the
   * smoother/forecast can be written back to the right LC positions. */
  N_filt = 0;
  for(i = 0; i < N_in; i++) {
    if(isnan(p->mag[lcnum][i])) continue;
    if(!(p->sig[lcnum][i] > 0.0)) continue;
    if(drw->usemask && drw->maskvar != NULL) {
      if(!(EvaluateVariable_Double(lc_name_num, lcnum, i,
                                   drw->maskvar) > VARTOOLS_MASK_TINY))
        continue;
    }
    t_filt[N_filt]   = p->t[lcnum][i];
    m_filt[N_filt]   = p->mag[lcnum][i];
    s2_filt[N_filt]  = p->sig[lcnum][i] * p->sig[lcnum][i];
    orig_idx[N_filt] = i;
    N_filt++;
  }

  if(N_filt < 3) {
    free(t_filt); free(m_filt); free(s2_filt); free(orig_idx);
    return;
  }

  /* mean=subtract: pre-subtract the weighted mean from the LC so the
   * recursion sees a mean-zero series.  Remember wm_subtracted so we
   * can put the smoother / forecast back on the original magnitude
   * scale when emitting models. */
  if(mean_mode == DRWFIT_MEAN_SUBTRACT) {
    double w, wsum = 0.0, wxsum = 0.0;
    for(i = 0; i < N_filt; i++) {
      w      = 1.0 / s2_filt[i];
      wsum  += w;
      wxsum += w * m_filt[i];
    }
    wm_subtracted = (wsum > 0.0 ? wxsum / wsum : 0.0);
    for(i = 0; i < N_filt; i++) m_filt[i] -= wm_subtracted;
    mu_fixed = 0.0;
  }

  drwfit_initial_guess(N_filt, t_filt, m_filt, s2_filt,
                       have_sigma0, sigma0_in,
                       have_tau0,   tau0_in,
                       have_mean0,  mean0_in,
                       &sigma_long_init, &tau_init, &mu_init);

  if(mean_mode == DRWFIT_MEAN_FIX)      mu_init = mu_fixed;
  if(mean_mode == DRWFIT_MEAN_SUBTRACT) mu_init = 0.0;

  converged = drwfit_run_amoeba(N_filt, t_filt, m_filt, s2_filt,
                                sigma_long_init, tau_init, mu_init,
                                mean_mode, mu_fixed,
                                &sigma_long_fit, &tau_fit, &mu_fit,
                                &neg2lnL_best);

  if(converged) {
    /* MacLeod 2010 quality indicators.  ln L = -(1/2) * (-2 ln L), so
     *   Delta L = ln L_best - ln L_alt = (1/2) (neg2lnL_alt - neg2lnL_best).
     * Both quantities are positive when L_best > L_alt (the data prefer
     * the DRW over the alternative). */
    neg2lnL_noise = kelly_neg2lnL_noise_limit(m_filt, s2_filt, N_filt);
    neg2lnL_inf   = kelly_neg2lnL_tau_infinity_limit(t_filt, m_filt, s2_filt,
                                                      N_filt,
                                                      sigma_long_fit, mu_fit);

    drw->sigma_long[lcnum] = sigma_long_fit;
    drw->tau[lcnum]        = tau_fit;
    drw->mu[lcnum]         = (mean_mode == DRWFIT_MEAN_SUBTRACT)
                             ? nan_val : mu_fit;
    drw->lnL[lcnum]        = -0.5 * neg2lnL_best;
    drw->dlnL_noise[lcnum] = 0.5 * (neg2lnL_noise - neg2lnL_best);
    drw->dlnL_inf[lcnum]   = 0.5 * (neg2lnL_inf   - neg2lnL_best);
    drw->converged[lcnum]  = 1;
  }

  /* --- Phase B: filter trajectory + smoother + correctlc/modelvar/save -- */

  /* Do we need to run the filter+smoother at all?  Only when there's a
   * downstream consumer (save, correctlc, modelvar) AND the fit converged. */
  need_forecast = (drw->do_save ||
                   (drw->do_correctlc && drw->correctlc_kind == DRWFIT_MODEL_FORECAST) ||
                   (drw->do_modelvar  && drw->modelvar_kind  == DRWFIT_MODEL_FORECAST));
  need_smoothed = (drw->do_save ||
                   (drw->do_correctlc && drw->correctlc_kind == DRWFIT_MODEL_SMOOTHED) ||
                   (drw->do_modelvar  && drw->modelvar_kind  == DRWFIT_MODEL_SMOOTHED));

  if(converged && (need_forecast || need_smoothed)) {
    double center_offset;
    /* center_offset is what we add to a centered (x*-scale) prediction
     * to put it on the same magnitude scale as p->mag[lcnum][i].  For
     * mean=fit/fix it's the fitted/fixed mu; for mean=subtract the
     * recursion saw a mean-zero LC so we add back the weighted-mean
     * offset that was subtracted from m_filt at the top. */
    center_offset = (mean_mode == DRWFIT_MEAN_SUBTRACT)
                    ? wm_subtracted : mu_fit;

    if((xhat_minus  = (double *) malloc(N_filt * sizeof(double))) == NULL ||
       (Omega_minus = (double *) malloc(N_filt * sizeof(double))) == NULL ||
       (xhat_plus   = (double *) malloc(N_filt * sizeof(double))) == NULL ||
       (Omega_plus  = (double *) malloc(N_filt * sizeof(double))) == NULL)
      vt_error(ERR_MEMALLOC);

    /* One final forward pass at the best-fit (sigma, tau, mu) with full
     * recording.  Likelihood return value discarded (it equals
     * neg2lnL_best up to amoeba tolerance). */
    (void) kelly_neg2lnL(t_filt, m_filt, s2_filt, N_filt,
                         sigma_long_fit, tau_fit, mu_fit,
                         xhat_minus, Omega_minus,
                         xhat_plus,  Omega_plus);

    if(need_smoothed) {
      if((x_smooth    = (double *) malloc(N_filt * sizeof(double))) == NULL ||
         (Omega_smooth = (double *) malloc(N_filt * sizeof(double))) == NULL)
        vt_error(ERR_MEMALLOC);
      kelly_rts_smoother(t_filt, N_filt, tau_fit,
                         xhat_minus, Omega_minus,
                         xhat_plus,  Omega_plus,
                         x_smooth, Omega_smooth);
    }

    /* Build per-original-LC-index model arrays (magnitude scale).  NaN
     * at indices that were filtered out at the top. */
    if((model_smoothed = (double *) malloc(N_in * sizeof(double))) == NULL ||
       (model_forecast = (double *) malloc(N_in * sizeof(double))) == NULL)
      vt_error(ERR_MEMALLOC);
    for(i = 0; i < N_in; i++) {
      model_smoothed[i] = nan_val;
      model_forecast[i] = nan_val;
    }
    /* xhat_minus and x_smooth are in centered (x* = x - mu) coordinates,
     * so adding center_offset converts to the magnitude scale matching
     * p->mag[lcnum][orig_idx[j]]. */
    for(j = 0; j < N_filt; j++) {
      int oi = orig_idx[j];
      if(need_forecast) model_forecast[oi] = center_offset + xhat_minus[j];
      if(need_smoothed) model_smoothed[oi] = center_offset + x_smooth[j];
    }

    /* "save" aux file: 8 columns over all N_in original indices.
     * Filtered-out rows have NaN for the model columns. */
    if(drw->do_save && outname != NULL && outname[0] != '\0') {
      FILE *fp = fopen(outname, "w");
      if(fp == NULL) {
        fprintf(stderr,
          "Error: cannot open '%s' for writing the drwfit aux file\n",
          outname);
        exit(3);
      }
      fprintf(fp, "# t x sig_meas x_hat_fwd Omega_fwd chi_fwd x_smoothed Omega_smoothed\n");
      /* Build an inverse map index_in_filtered[orig_index] = j or -1. */
      int *inv = (int *) malloc(N_in * sizeof(int));
      if(inv == NULL) vt_error(ERR_MEMALLOC);
      for(i = 0; i < N_in; i++) inv[i] = -1;
      for(j = 0; j < N_filt; j++) inv[orig_idx[j]] = j;
      for(i = 0; i < N_in; i++) {
        double xv = p->mag[lcnum][i];
        double sv = p->sig[lcnum][i];
        if(inv[i] < 0) {
          fprintf(fp, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                  p->t[lcnum][i], xv, sv,
                  nan_val, nan_val, nan_val, nan_val, nan_val);
        } else {
          j = inv[i];
          double xhat_fwd_mag    = center_offset + xhat_minus[j];
          double xsmooth_mag     = center_offset + x_smooth[j];
          double var_fwd         = Omega_minus[j] + s2_filt[j];
          /* The recursion's residual is (x*_i - x_hat_minus_i) where
           * x*_i = x_i - mu.  m_filt has wm_subtracted removed up-front
           * in subtract mode (mu_fit = 0); in fit/fix mode m_filt is the
           * original x and we still subtract mu_fit.  Both reduce to
           * "centered observation minus centered prior". */
          double x_star          = m_filt[j] - mu_fit;
          double chi_fwd         = (x_star - xhat_minus[j]) / sqrt(var_fwd);
          fprintf(fp, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                  p->t[lcnum][i], xv, sv,
                  xhat_fwd_mag, Omega_minus[j], chi_fwd,
                  xsmooth_mag, Omega_smooth[j]);
        }
      }
      free(inv);
      fclose(fp);
    }

    /* "modelvar": copy chosen model into the per-LC variable.  NaN was
     * pre-written for every LC index at the top of this function, so
     * filtered-out points stay NaN automatically. */
    if(drw->do_modelvar && drw->modelvar != NULL) {
      double *src = (drw->modelvar_kind == DRWFIT_MODEL_FORECAST)
                    ? model_forecast : model_smoothed;
      for(i = 0; i < N_in; i++)
        (*((double ***) drw->modelvar->dataptr))[lcnum][i] = src[i];
    }

    /* "correctlc": subtract chosen model from p->mag[lcnum] in place.
     * Filtered-out points (NaN model) left untouched. */
    if(drw->do_correctlc) {
      double *src = (drw->correctlc_kind == DRWFIT_MODEL_FORECAST)
                    ? model_forecast : model_smoothed;
      for(i = 0; i < N_in; i++) {
        if(!isnan(src[i])) p->mag[lcnum][i] -= src[i];
      }
    }
  }

  if(xhat_minus)     free(xhat_minus);
  if(Omega_minus)    free(Omega_minus);
  if(xhat_plus)      free(xhat_plus);
  if(Omega_plus)     free(Omega_plus);
  if(x_smooth)       free(x_smooth);
  if(Omega_smooth)   free(Omega_smooth);
  if(model_smoothed) free(model_smoothed);
  if(model_forecast) free(model_forecast);
  free(t_filt); free(m_filt); free(s2_filt); free(orig_idx);
}
