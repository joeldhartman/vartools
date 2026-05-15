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
/*  Fast Template Periodogram (FTP) period search for vartools.
 *
 *  This file was developed with assistance from Claude (Anthropic) by
 *  Joel Hartman.
 *
 *  References:
 *    Hoffman, J., VanderPlas, J., Hartman, J. D., & Bakos, G. A. 2021,
 *        arXiv:2101.12348, "A Fast Template Periodogram for Detecting
 *        Non-sinusoidal Fixed-shape Signals in Irregularly Sampled
 *        Time Series".
 *    Zechmeister, M., Kuerster, M. 2009, A&A, 496, 577
 *        Generalized Lomb-Scargle (GLS); FTP is the template-fit extension.
 *    Sesar, B., Hernitschek, N., Mitrovic, S. et al. 2017, AJ, 153, 204
 *        RR Lyrae template-fitting motivation for the FTP method.
 *    pyftp (https://github.com/PrincetonUniversity/FastTemplatePeriodogram)
 *        Reference Python implementation by Hoffman & VanderPlas.
 *
 *  Algorithm: At each trial frequency omega, fit
 *      y_i ~= theta_1 * M(omega t_i - theta_2) + theta_3
 *  where M is a known periodic template approximated by a truncated Fourier
 *  series of length H:
 *      M(phi) = sum_{n=1..H} [c_n cos(n phi) + s_n sin(n phi)].
 *  The best-fit (theta_1, theta_2, theta_3) minimises the weighted sum of
 *  squared residuals; the periodogram value is
 *      P(omega) = (YM)^2 / (YY * MM)         in [0, 1]
 *  where YM = Cov(M_{theta_2}, y), MM = Var(M_{theta_2}), YY = Var(y).
 *  P -> 1 means the template fits the data exactly at this frequency;
 *  P -> 0 means no improvement over a constant fit.
 *
 *  Initial scope (Phase A of the FTP rollout):
 *    - Template input: single-file, two-column whitespace-separated ASCII
 *                      with rows (c_n s_n), n=1..H.  H inferred from the
 *                      file's row count.
 *    - Per-frequency optimisation: brute-force scan of theta_2 in [0, pi)
 *                      (the +/- sign of sin(theta_2) is enumerated
 *                      separately) at fine resolution + Brent refinement.
 *                      Slower than the polynomial root-finding fast path
 *                      from the paper, but algorithmically correct to
 *                      O(1e-12) by construction.  The polynomial fast path
 *                      is deferred to a later phase, with this brute-force
 *                      implementation as the validator.
 *    - Per-frequency summations: direct loops over data points (no NFFT).
 *                      O(N*H) per trial frequency.
 *    - Negative-amplitude policy: theta_1 < 0 corresponds to a flipped
 *                      template, NOT a phase shift, for non-sinusoidal
 *                      templates.  By default we report the best fit
 *                      regardless of theta_1's sign, with FTP_NegAmp_N_M
 *                      = 1 when the best fit had theta_1 < 0.  The
 *                      "posamponly" keyword forces the search to skip
 *                      negative-amplitude fits.
 *    - Weighting:      default 1/sigma^2 if the input LC has meaningful
 *                      errors; "noerr" keyword forces uniform weighting.
 *    - Peaks:          Npeaks via coarse + fine-tune sweep + period-multiple
 *                      double-check, mirroring -aov / -PDM.
 *    - SNR:            (peak - clipped_mean) / clipped_rms over the
 *                      periodogram, with a hardcoded 5-sigma single-pass
 *                      iterative clip that excludes the bright peaks.
 *                      Sign convention is positive when the peak rises
 *                      ABOVE the noise (opposite of PDM).
 *    - FAP:            Not reported in Phase A.  The Hoffman et al. paper
 *                      does not derive an analytic null distribution for
 *                      P(omega), and ad-hoc GLS-style approximations are
 *                      biased by the polynomial-amplitude search.  A
 *                      bootstrap FAP is the natural follow-up (Phase D).
 *
 *  Deferred to subsequent commits:
 *    - Polynomial root-finding fast path (Phase B).  Cross-validate
 *      against this brute-force baseline.
 *    - Filelist + inline template input modes (Phase B).
 *    - clip / noerr / maskpoints in their canonical strict-order positions
 *      (Phase C).
 *    - whiten / fixperiodSNR / bootstrap (Phase D).
 */

#include "commands.h"
#include "programdata.h"
#include "functions.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* --- Constants --- */
#define FTP_ERROR_SCORE      (-1.0)            /* error sentinel for P(omega); valid P in [0,1] */
#define FTP_MAX_DOUBLE_CHECK_MULTIPLE 19
#define FTP_MAX_PERIOD_DIFF_MULTIPLE  5
#define FTP_MINVAR           1.0e-32
#define FTP_DEFAULT_NTHETA   720               /* 0.5 deg phase resolution per (b, sgn) branch */
#define FTP_BRENT_TOL        1.0e-12
#define FTP_BRENT_MAX_ITER   100
#define FTP_DEFAULT_CLIP     5.0
#define FTP_DEFAULT_CLIPITER 1

/* ---- Template loader ---- */

/* Read a 2-column whitespace-separated ASCII file with H rows, columns
 * (c_n s_n) for n=1..H.  Lines beginning with '#' or empty lines are skipped.
 * Allocates *cn_out and *sn_out (caller frees) and returns H.  Returns 0 on
 * parse error or empty file. */
int ftp_load_template_file(const char *path, double **cn_out, double **sn_out)
{
  FILE *fp;
  char *line = NULL;
  size_t linelen = 0;
  int H_alloc = 16, H = 0;
  double *cn = NULL, *sn = NULL;
  ssize_t got;

  *cn_out = *sn_out = NULL;

  if ((fp = fopen(path, "r")) == NULL) {
    fprintf(stderr, "-FTP: cannot open template file '%s' for reading\n", path);
    return 0;
  }

  cn = (double *) malloc(H_alloc * sizeof(double));
  sn = (double *) malloc(H_alloc * sizeof(double));
  if (cn == NULL || sn == NULL) {
    fclose(fp);
    if (cn != NULL) free(cn);
    if (sn != NULL) free(sn);
    return 0;
  }

  while ((got = getline(&line, &linelen, fp)) != -1) {
    char *p = line;
    double c_val, s_val;
    int nread;

    /* Skip leading whitespace */
    while (*p && isspace((unsigned char)*p)) p++;
    /* Skip blank lines and comments */
    if (*p == '\0' || *p == '#') continue;

    nread = sscanf(p, "%lf %lf", &c_val, &s_val);
    if (nread != 2) {
      fprintf(stderr, "-FTP: malformed template line in '%s' "
                      "(expected two whitespace-separated doubles): %s",
              path, line);
      free(line); free(cn); free(sn); fclose(fp);
      return 0;
    }

    if (H >= H_alloc) {
      H_alloc *= 2;
      cn = (double *) realloc(cn, H_alloc * sizeof(double));
      sn = (double *) realloc(sn, H_alloc * sizeof(double));
      if (cn == NULL || sn == NULL) {
        fprintf(stderr, "-FTP: out of memory loading template '%s'\n", path);
        if (line != NULL) free(line);
        if (cn != NULL) free(cn);
        if (sn != NULL) free(sn);
        fclose(fp);
        return 0;
      }
    }
    cn[H] = c_val;
    sn[H] = s_val;
    H++;
  }

  if (line != NULL) free(line);
  fclose(fp);

  if (H == 0) {
    fprintf(stderr, "-FTP: template file '%s' has no usable rows\n", path);
    free(cn); free(sn);
    return 0;
  }

  *cn_out = cn;
  *sn_out = sn;
  return H;
}

/* ---- Period-similarity helper (mirror aov.c / pdm.c) ---- */

static int ftp_isDifferentPeriods(double p1, double p2, double T)
{
  int a, b;
  double p1mul;
  if (T * fabs(p2 - p1) < p1 * p2)
    return 0;
  for (a = 1; a < FTP_MAX_PERIOD_DIFF_MULTIPLE; a++)
    for (b = a + 1; b <= FTP_MAX_PERIOD_DIFF_MULTIPLE; b++) {
      p1mul = p1 * b / a;
      if (T * fabs(p2 - p1mul) < p2 * p1mul)
        return 0;
    }
  return 1;
}

/* ---- Per-frequency direct summations ----
 *
 * Computes
 *   C[n]  = sum_i w_i  cos(n omega t_i)        for n = 1..twoH
 *   S[n]  = sum_i w_i  sin(n omega t_i)
 *   YC[n] = sum_i w_i  (y_i - ybar) cos(n omega t_i)   for n = 1..H
 *   YS[n] = sum_i w_i  (y_i - ybar) sin(n omega t_i)
 *
 * (twoH = 2*H is needed for the bilinear forms CC/CS/SS, below.)
 * Storage: C[0..twoH-1], S[0..twoH-1], YC[0..H-1], YS[0..H-1]
 * (i.e. C[k] holds the n=k+1 sum; zero index unused for the n=0 trivial
 * sum which is 1 = sum w_i by construction).
 *
 * Cost: O(N * twoH) per call.  For each point i we compute cos(omega t_i)
 * and sin(omega t_i) once, then iterate the multi-angle recurrence
 *   cos((n+1) x) = 2 cos(x) cos(n x) - cos((n-1) x)
 *   sin((n+1) x) = 2 cos(x) sin(n x) - sin((n-1) x).
 */
static void ftp_compute_sums(int N, const double *t, const double *y_centered,
                              const double *w, double omega, int H,
                              double *C, double *S, double *YC, double *YS)
{
  int twoH = 2 * H;
  int i, n;

  for (n = 0; n < twoH; n++) C[n] = S[n] = 0.0;
  for (n = 0; n < H; n++)    YC[n] = YS[n] = 0.0;

  for (i = 0; i < N; i++) {
    double wi = w[i];
    double yi = y_centered[i];
    double x  = omega * t[i];
    double cn_prev, sn_prev, cn_curr, sn_curr, cn_next, sn_next;
    double two_cos_x;

    /* n = 1 directly */
    cn_curr = cos(x);
    sn_curr = sin(x);
    two_cos_x = 2.0 * cn_curr;

    cn_prev = 1.0;          /* cos(0) */
    sn_prev = 0.0;          /* sin(0) */

    /* n = 1 contributions */
    C[0] += wi * cn_curr;
    S[0] += wi * sn_curr;
    YC[0] += wi * yi * cn_curr;
    YS[0] += wi * yi * sn_curr;

    /* n = 2..twoH via Chebyshev-style recurrence */
    for (n = 2; n <= twoH; n++) {
      cn_next = two_cos_x * cn_curr - cn_prev;
      sn_next = two_cos_x * sn_curr - sn_prev;
      C[n - 1] += wi * cn_next;
      S[n - 1] += wi * sn_next;
      if (n <= H) {
        YC[n - 1] += wi * yi * cn_next;
        YS[n - 1] += wi * yi * sn_next;
      }
      cn_prev = cn_curr;
      sn_prev = sn_curr;
      cn_curr = cn_next;
      sn_curr = sn_next;
    }
  }
}

/* ---- Bilinear sums CC, CS, SS ----
 *
 * Computed from C[1..2H], S[1..2H] (and the trivial C[0]=1, S[0]=0):
 *   CC[n,m] = sum_i w_i cos(n omega t_i) cos(m omega t_i) - C_n C_m
 *           = 0.5 * (Cf(|n-m|) + Cf(n+m)) - C_n C_m
 *   CS[n,m] = sum_i w_i cos(n omega t_i) sin(m omega t_i) - C_n S_m
 *           = 0.5 * (sgn(m-n) Sf(|m-n|) + Sf(n+m)) - C_n S_m
 *   SS[n,m] = sum_i w_i sin(n omega t_i) sin(m omega t_i) - S_n S_m
 *           = 0.5 * (Cf(|n-m|) - Cf(n+m)) - S_n S_m
 *
 * where Cf, Sf are the "full" sums: Cf(0) = 1, Sf(0) = 0, Cf(k) = C[k-1]
 * for k=1..2H, Sf(k) = S[k-1].
 *
 * CC and SS are symmetric.  CS is NOT symmetric: CS[n,m] != CS[m,n] in
 * general.
 *
 * Storage: row-major, CC[n*H + m] for n,m in [0,H-1] (i.e. n=index+1, m=index+1).
 */
static void ftp_compute_bilinears(int H, const double *C, const double *S,
                                   double *CC, double *CS, double *SS)
{
  int n, m;
  /* Indexing helper: Cf(k) and Sf(k) for k in [0, 2H]. */
  #define CF(k) ((k) == 0 ? 1.0 : C[(k) - 1])
  #define SF(k) ((k) == 0 ? 0.0 : S[(k) - 1])
  for (n = 1; n <= H; n++) {
    for (m = 1; m <= H; m++) {
      int diff = n - m;
      int sum  = n + m;
      int absd = diff < 0 ? -diff : diff;
      double Cnm_diff = CF(absd);
      /* CS uses sgn(m - n); diff = n - m so sgn(m-n) = -sgn(diff). */
      double sgn_mn = (diff > 0) ? -1.0 : (diff < 0 ? 1.0 : 0.0);
      double Snm_diff_csign = sgn_mn * SF(absd);
      double Cn = CF(n), Cm = CF(m), Sn = SF(n), Sm = SF(m);
      CC[(n - 1) * H + (m - 1)] = 0.5 * (Cnm_diff + CF(sum)) - Cn * Cm;
      CS[(n - 1) * H + (m - 1)] = 0.5 * (Snm_diff_csign + SF(sum)) - Cn * Sm;
      SS[(n - 1) * H + (m - 1)] = 0.5 * (Cnm_diff - CF(sum)) - Sn * Sm;
    }
  }
  #undef CF
  #undef SF
}

/* ---- P(omega, theta_2) evaluator ----
 *
 * Given a phase shift theta_2 (radians) and the precomputed summation
 * vectors / matrices at the current frequency, compute
 *
 *   A_n = c_n cos(n theta_2) - s_n sin(n theta_2)
 *   B_n = s_n cos(n theta_2) + c_n sin(n theta_2)
 *
 * (verified against pyftp's Avec / Bvec via the Chebyshev identities
 *  T_n(cos t) = cos(n t),  U_{n-1}(cos t) sin t = sin(n t))
 *
 *   YM = sum_n [A_n YC_n + B_n YS_n]
 *   MM = sum_{n,m} [A_n A_m CC_nm + (A_n B_m + A_m B_n) CS_nm + B_n B_m SS_nm]
 *
 * CS is not symmetric, so the cross term uses (A_n B_m + A_m B_n).
 * P(omega, theta_2) = (YM)^2 / (YY * MM), bounded in [0, 1].
 *
 * Caller provides scratch arrays scratch_A, scratch_B (each of length >= H)
 * to avoid per-call allocation.  Writes the sign of YM (= sign of optimal
 * theta_1) into *sign_out when non-NULL; flag is +1 for non-negative
 * amplitudes and -1 for negative amplitudes (flipped template).
 */
static double ftp_eval_P(int H, double theta_2,
                          const double *cn, const double *sn,
                          const double *YC, const double *YS,
                          const double *CC, const double *CS, const double *SS,
                          double YY,
                          double *scratch_A, double *scratch_B,
                          int *sign_out)
{
  int n, m;
  double YM, MM;
  double c1, s1, cprev, sprev, ccurr, scurr;

  /* A_n, B_n by trigonometric recurrence on n*theta_2.  Standard double-angle
   * recurrence: cos((n+1)t) = 2 cos(t) cos(nt) - cos((n-1)t),  sin similarly. */
  c1 = cos(theta_2);
  s1 = sin(theta_2);
  cprev = 1.0; sprev = 0.0;
  ccurr = c1;  scurr = s1;
  scratch_A[0] = cn[0] * ccurr - sn[0] * scurr;
  scratch_B[0] = sn[0] * ccurr + cn[0] * scurr;
  for (n = 1; n < H; n++) {
    double cnext = 2.0 * c1 * ccurr - cprev;
    double snext = 2.0 * c1 * scurr - sprev;
    scratch_A[n] = cn[n] * cnext - sn[n] * snext;
    scratch_B[n] = sn[n] * cnext + cn[n] * snext;
    cprev = ccurr; sprev = scurr;
    ccurr = cnext; scurr = snext;
  }

  YM = 0.0;
  for (n = 0; n < H; n++) YM += scratch_A[n] * YC[n] + scratch_B[n] * YS[n];

  /* MM = Var(M_{theta_2}) = sum_{n,m} [
   *        A_n A_m CC_{nm}              (CC is symmetric)
   *      + A_n B_m CS_{nm} + B_n A_m CS_{mn}     (CS is NOT symmetric)
   *      + B_n B_m SS_{nm}              (SS is symmetric)
   *      ].
   * The CS pair simplifies to 2 sum_{n,m} A_n B_m CS_{nm} after renaming
   * (n,m) <-> (m,n) in the second sum -- the factor 2 is the right
   * coefficient, NOT the (A_n B_m + A_m B_n) symmetrisation that would be
   * correct only if CS were symmetric. */
  MM = 0.0;
  for (n = 0; n < H; n++) {
    double An = scratch_A[n], Bn = scratch_B[n];
    for (m = 0; m < H; m++) {
      int idx = n * H + m;
      double Am = scratch_A[m], Bm = scratch_B[m];
      MM += An * Am * CC[idx]
          + 2.0 * An * Bm * CS[idx]
          + Bn * Bm * SS[idx];
    }
  }

  /* MM and YY are variances and must be non-negative.  Numerical cancellation
   * in the bilinear sums can produce tiny-negative MM at frequencies where the
   * template is ill-conditioned for the data; treat these as degenerate. */
  if (!(MM > FTP_MINVAR) || !(YY > FTP_MINVAR)) {
    if (sign_out != NULL) *sign_out = 1;
    return 0.0;
  }

  {
    double P = (YM * YM) / (YY * MM);
    /* By Cauchy-Schwarz |YM|^2 <= YY * MM, so P <= 1 mathematically.  Numerical
     * cancellation in MM can produce P slightly > 1; clamp. */
    if (P > 1.0) P = 1.0;
    if (P < 0.0) P = 0.0;
    if (sign_out != NULL) *sign_out = (YM >= 0.0) ? 1 : -1;
    return P;
  }
}

/* ---- Brent's-method-style refinement around a coarse-scan maximum ----
 *
 * Given a bracket (theta_lo, theta_mid, theta_hi) with f(theta_mid) >
 * f(theta_lo) and f(theta_mid) > f(theta_hi), refine theta_mid to local
 * max via golden-section.  Returns the refined theta_2 and writes the
 * refined P into *P_out.
 *
 * Plain golden-section is sufficient here: f is smooth in theta_2
 * (analytic in (cos(n theta_2), sin(n theta_2))) and we already have a
 * tight bracket from the coarse scan.  ~30 iterations bring the bracket
 * to ~1e-10 in theta_2.
 */
static double ftp_golden_refine(double theta_lo, double theta_mid, double theta_hi,
                                  double f_mid,
                                  int H, const double *cn, const double *sn,
                                  const double *YC, const double *YS,
                                  const double *CC, const double *CS,
                                  const double *SS, double YY,
                                  double *scratch_A, double *scratch_B,
                                  double *P_out, int *sign_out)
{
  const double phi_inv = 0.6180339887498949;     /* 1/phi */
  double a = theta_lo, b = theta_hi;
  double c = b - phi_inv * (b - a);
  double d = a + phi_inv * (b - a);
  double fc, fd;
  int iter;

  fc = ftp_eval_P(H, c, cn, sn, YC, YS, CC, CS, SS, YY, scratch_A, scratch_B, NULL);
  fd = ftp_eval_P(H, d, cn, sn, YC, YS, CC, CS, SS, YY, scratch_A, scratch_B, NULL);

  for (iter = 0; iter < FTP_BRENT_MAX_ITER; iter++) {
    if (fabs(b - a) < FTP_BRENT_TOL * (fabs(c) + fabs(d) + 1.0)) break;
    if (fc > fd) {
      b = d;
      d = c;
      fd = fc;
      c = b - phi_inv * (b - a);
      fc = ftp_eval_P(H, c, cn, sn, YC, YS, CC, CS, SS, YY, scratch_A, scratch_B, NULL);
    } else {
      a = c;
      c = d;
      fc = fd;
      d = a + phi_inv * (b - a);
      fd = ftp_eval_P(H, d, cn, sn, YC, YS, CC, CS, SS, YY, scratch_A, scratch_B, NULL);
    }
  }

  {
    double theta_best = (fc > fd) ? c : d;
    double f_best = (fc > fd) ? fc : fd;
    if (f_best < f_mid) {
      /* Coarse-scan midpoint won out -- keep it. */
      *P_out = f_mid;
      if (sign_out != NULL)
        ftp_eval_P(H, theta_mid, cn, sn, YC, YS, CC, CS, SS, YY,
                    scratch_A, scratch_B, sign_out);
      return theta_mid;
    }
    *P_out = f_best;
    if (sign_out != NULL)
      ftp_eval_P(H, theta_best, cn, sn, YC, YS, CC, CS, SS, YY,
                  scratch_A, scratch_B, sign_out);
    return theta_best;
  }
}

/* ---- Per-frequency P(omega) maximised over theta_2 ----
 *
 * Brute-force scan of theta_2 in [0, 2pi) at Ntheta uniform samples,
 * locating the global maximum of P(omega, theta_2).  Refine the best
 * sample with golden-section.  If allow_neg_amp == 0, samples whose YM
 * is negative are excluded (theta_1 = YM/MM must be non-negative).
 *
 * Returns peak P; writes optimal theta_2 to *theta_out and the sign of
 * theta_1 (i.e. sign(YM)) to *sign_out.
 *
 * Cost: Ntheta * O(H^2) per call (the H^2 is the MM bilinear form).
 */
static double ftp_max_P_over_theta(int H, int Ntheta, int allow_neg_amp,
                                     const double *cn, const double *sn,
                                     const double *YC, const double *YS,
                                     const double *CC, const double *CS,
                                     const double *SS, double YY,
                                     double *scratch_A, double *scratch_B,
                                     double *Pgrid, int *signgrid,
                                     double *theta_out, int *sign_out)
{
  int i;
  double dtheta = 2.0 * M_PI / Ntheta;
  double P_best = 0.0;
  double theta_best = 0.0;
  int sign_best = 1;
  int idx_best = -1;

  for (i = 0; i < Ntheta; i++) {
    double th = i * dtheta;
    int sgn;
    double P = ftp_eval_P(H, th, cn, sn, YC, YS, CC, CS, SS, YY,
                           scratch_A, scratch_B, &sgn);
    Pgrid[i] = P;
    signgrid[i] = sgn;
    if (!allow_neg_amp && sgn < 0) continue;
    if (P > P_best) {
      P_best = P;
      theta_best = th;
      sign_best = sgn;
      idx_best = i;
    }
  }

  /* Refine around idx_best if both neighbours bracket it (a clean local
   * maximum) and pass the sign filter. */
  if (idx_best >= 0 && Ntheta >= 3) {
    int ilo = (idx_best - 1 + Ntheta) % Ntheta;
    int ihi = (idx_best + 1) % Ntheta;
    int ok_lo = allow_neg_amp || signgrid[ilo] >= 0;
    int ok_hi = allow_neg_amp || signgrid[ihi] >= 0;
    if (ok_lo && ok_hi && Pgrid[idx_best] > Pgrid[ilo] && Pgrid[idx_best] > Pgrid[ihi]) {
      double th_lo = theta_best - dtheta;
      double th_hi = theta_best + dtheta;
      int sgn_refined;
      double P_refined;
      theta_best = ftp_golden_refine(th_lo, theta_best, th_hi, P_best,
                                      H, cn, sn, YC, YS, CC, CS, SS, YY,
                                      scratch_A, scratch_B,
                                      &P_refined, &sgn_refined);
      if (allow_neg_amp || sgn_refined >= 0) {
        if (P_refined > P_best) {
          P_best = P_refined;
          sign_best = sgn_refined;
        }
      }
    }
  }

  /* Wrap theta into [0, 2pi). */
  while (theta_best < 0.0)        theta_best += 2.0 * M_PI;
  while (theta_best >= 2.0 * M_PI) theta_best -= 2.0 * M_PI;

  if (theta_out != NULL) *theta_out = theta_best;
  if (sign_out  != NULL) *sign_out  = sign_best;
  return P_best;
}

/* ---- Per-LC scratch workspace ----
 *
 * Bundles the per-frequency working memory so callers allocate once per
 * LC and pass scratch down through ftp_periodogram and ftp_one_period.
 */
typedef struct {
  int H;
  double *C;          /* size 2H */
  double *S;          /* size 2H */
  double *YC;         /* size H  */
  double *YS;         /* size H  */
  double *CC;         /* size H*H */
  double *CS;         /* size H*H */
  double *SS;         /* size H*H */
  double *A;          /* size H, used by ftp_eval_P */
  double *B;          /* size H */
  double *Pgrid;      /* size FTP_DEFAULT_NTHETA */
  int    *signgrid;   /* size FTP_DEFAULT_NTHETA */
} _FTPScratch;

static void ftp_scratch_alloc(_FTPScratch *s, int H, int Ntheta)
{
  int twoH = 2 * H;
  s->H = H;
  s->C  = (double *) malloc(twoH * sizeof(double));
  s->S  = (double *) malloc(twoH * sizeof(double));
  s->YC = (double *) malloc(H    * sizeof(double));
  s->YS = (double *) malloc(H    * sizeof(double));
  s->CC = (double *) malloc(H*H  * sizeof(double));
  s->CS = (double *) malloc(H*H  * sizeof(double));
  s->SS = (double *) malloc(H*H  * sizeof(double));
  s->A  = (double *) malloc(H    * sizeof(double));
  s->B  = (double *) malloc(H    * sizeof(double));
  s->Pgrid    = (double *) malloc(Ntheta * sizeof(double));
  s->signgrid = (int *)    malloc(Ntheta * sizeof(int));
  if (s->C == NULL || s->S == NULL || s->YC == NULL || s->YS == NULL ||
      s->CC == NULL || s->CS == NULL || s->SS == NULL ||
      s->A == NULL || s->B == NULL || s->Pgrid == NULL || s->signgrid == NULL)
    vt_error(ERR_MEMALLOC);
}

static void ftp_scratch_free(_FTPScratch *s)
{
  free(s->C); free(s->S); free(s->YC); free(s->YS);
  free(s->CC); free(s->CS); free(s->SS);
  free(s->A); free(s->B);
  free(s->Pgrid); free(s->signgrid);
  memset(s, 0, sizeof(_FTPScratch));
}

/* Compute periodogram across periods[0..Nperiod-1]. */
static void ftp_periodogram(int N, double *t, double *y_centered, double *w,
                             double YY, int H, double *cn, double *sn,
                             int allow_neg_amp, _FTPScratch *sc,
                             int Nperiod, double *periods, double *periodogram,
                             int *negamp_grid, double *theta_grid)
{
  int i;
  for (i = 0; i < Nperiod; i++) {
    double omega = 2.0 * M_PI / periods[i];
    double th = 0.0;
    int sgn = 1;
    ftp_compute_sums(N, t, y_centered, w, omega, H, sc->C, sc->S, sc->YC, sc->YS);
    ftp_compute_bilinears(H, sc->C, sc->S, sc->CC, sc->CS, sc->SS);
    periodogram[i] = ftp_max_P_over_theta(H, FTP_DEFAULT_NTHETA, allow_neg_amp,
                                           cn, sn, sc->YC, sc->YS,
                                           sc->CC, sc->CS, sc->SS, YY,
                                           sc->A, sc->B, sc->Pgrid, sc->signgrid,
                                           &th, &sgn);
    if (negamp_grid != NULL) negamp_grid[i] = (sgn < 0) ? 1 : 0;
    if (theta_grid  != NULL) theta_grid[i]  = th;
  }
}

/* One-shot P(omega) at a single test period (used in fine-tune + multiple checks). */
static double ftp_one_period(int N, double *t, double *y_centered, double *w,
                              double YY, int H, double *cn, double *sn,
                              int allow_neg_amp, _FTPScratch *sc, double period,
                              double *theta_out, int *sign_out)
{
  double omega, P;

  if (period <= 0.0) {
    if (theta_out != NULL) *theta_out = 0.0;
    if (sign_out  != NULL) *sign_out  = 1;
    return FTP_ERROR_SCORE;
  }
  omega = 2.0 * M_PI / period;
  ftp_compute_sums(N, t, y_centered, w, omega, H, sc->C, sc->S, sc->YC, sc->YS);
  ftp_compute_bilinears(H, sc->C, sc->S, sc->CC, sc->CS, sc->SS);
  P = ftp_max_P_over_theta(H, FTP_DEFAULT_NTHETA, allow_neg_amp,
                            cn, sn, sc->YC, sc->YS, sc->CC, sc->CS, sc->SS, YY,
                            sc->A, sc->B, sc->Pgrid, sc->signgrid,
                            theta_out, sign_out);
  return P;
}

/* ---- Per-LC peak finder ---- */

void findPeaks_ftp(double *t_, double *mag_, double *sig_, int N,
                    int H, double *cn, double *sn,
                    int useerr, int allow_neg_amp,
                    double minP, double maxP, double subsample, double finetune,
                    int Npeaks,
                    double *perpeaks, double *Ppeaks,
                    int *peakNegAmp,
                    double *peakSNR,
                    double *peakTheta,
                    double *avePower, double *rmsPower,
                    int outflag, char *outname, int ascii,
                    int lcnum, int lc_name_num)
{
  int i, j, k, foundsofar, test, Nperiod, a, b, abest, bbest, ismultiple;
  int Nused;
  double *t = NULL, *mag = NULL, *w = NULL, *y_centered = NULL;
  double *periods = NULL, *periodogram = NULL;
  int *negamp_grid = NULL;
  double *theta_grid = NULL;
  double T, mu_total, var_total, YY;
  double sumw, wnorm;
  double freq, minfreq, freqstep, smallfreqstep;
  double testperiod, bestscore, lastpoint;
  double ave_per, std_per;
  double clip_val = FTP_DEFAULT_CLIP;
  int    clipiter = FTP_DEFAULT_CLIPITER;
  long double Sum, Sumsqr;
  FILE *outfile = NULL;
  _FTPScratch sc = {0};

  /* Avoid -Wunused warnings on parameters not used in Phase A scope. */
  (void) lcnum; (void) lc_name_num;

  if (N <= 1 || H <= 0) {
    for (i = 0; i < Npeaks; i++) {
      perpeaks[i]   = -1.0;
      Ppeaks[i]     = 0.0;
      peakSNR[i]    = 0.0;
      peakNegAmp[i] = 0;
      peakTheta[i]  = 0.0;
    }
    *avePower = 0.0;
    *rmsPower = 1.0;
    return;
  }

  /* Copy + filter LC. */
  t          = (double *) malloc(N * sizeof(double));
  mag        = (double *) malloc(N * sizeof(double));
  w          = (double *) malloc(N * sizeof(double));
  y_centered = (double *) malloc(N * sizeof(double));
  if (t == NULL || mag == NULL || w == NULL || y_centered == NULL) vt_error(ERR_MEMALLOC);

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
      Ppeaks[i]     = 0.0;
      peakSNR[i]    = 0.0;
      peakNegAmp[i] = 0;
      peakTheta[i]  = 0.0;
    }
    *avePower = 0.0;
    *rmsPower = 1.0;
    free(t); free(mag); free(w); free(y_centered);
    return;
  }

  /* Normalise weights to sum to 1. */
  sumw = 0.0;
  for (i = 0; i < Nused; i++) sumw += w[i];
  if (sumw <= 0.0) sumw = 1.0;
  wnorm = 1.0 / sumw;
  for (i = 0; i < Nused; i++) w[i] *= wnorm;

  /* Weighted mean and variance. */
  mu_total = 0.0;
  for (i = 0; i < Nused; i++) mu_total += w[i] * mag[i];
  YY = 0.0;
  for (i = 0; i < Nused; i++) {
    double dy = mag[i] - mu_total;
    y_centered[i] = dy;
    YY += w[i] * dy * dy;
  }
  var_total = YY;
  if (var_total < FTP_MINVAR) {
    for (i = 0; i < Npeaks; i++) {
      perpeaks[i]   = -1.0;
      Ppeaks[i]     = 0.0;
      peakSNR[i]    = 0.0;
      peakNegAmp[i] = 0;
      peakTheta[i]  = 0.0;
    }
    *avePower = 0.0;
    *rmsPower = 1.0;
    free(t); free(mag); free(w); free(y_centered);
    return;
  }

  /* Frequency grid. */
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

  periods     = (double *) malloc(Nperiod * sizeof(double));
  periodogram = (double *) malloc(Nperiod * sizeof(double));
  negamp_grid = (int *)    malloc(Nperiod * sizeof(int));
  theta_grid  = (double *) malloc(Nperiod * sizeof(double));
  if (periods == NULL || periodogram == NULL || negamp_grid == NULL || theta_grid == NULL)
    vt_error(ERR_MEMALLOC);

  freq = 1.0 / minP;
  for (i = 0; i < Nperiod && freq >= minfreq; i++) {
    periods[i] = 1.0 / freq;
    freq -= freqstep;
  }
  Nperiod = i;

  ftp_scratch_alloc(&sc, H, FTP_DEFAULT_NTHETA);

  /* Compute periodogram. */
  ftp_periodogram(Nused, t, y_centered, w, YY, H, cn, sn,
                  allow_neg_amp, &sc, Nperiod, periods, periodogram,
                  negamp_grid, theta_grid);

  /* Mean / RMS over the periodogram (for SNR). */
  Sum = 0.0; Sumsqr = 0.0;
  int Ngood = 0;
  for (i = 0; i < Nperiod; i++) {
    double v = periodogram[i];
    if (v >= 0.0 && v * 0.0 == 0.0) {
      Sum    += (long double) v;
      Sumsqr += (long double)(v * v);
      Ngood++;
    }
  }
  if (Ngood > 0) {
    Sum /= Ngood;
    Sumsqr /= Ngood;
    ave_per = (double) Sum;
    std_per = sqrt((double)(Sumsqr - Sum * Sum));
    if (!(std_per > 0.0)) std_per = 1.0;
  } else {
    ave_per = 0.0;
    std_per = 1.0;
  }

  /* Iterative sigma-clip the peaks (which lie ABOVE the noise) for SNR purposes. */
  if (clip_val > 0.0 && Ngood > 0) {
    int nclippedthis = 0, nclippedlast;
    do {
      nclippedlast = nclippedthis;
      nclippedthis = 0;
      Sum = 0.0; Sumsqr = 0.0;
      int Ngood2 = 0;
      for (i = 0; i < Nperiod; i++) {
        double v = periodogram[i];
        if (v >= 0.0 && v * 0.0 == 0.0) {
          if (v < ave_per + clip_val * std_per) {
            Sum    += (long double) v;
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
        std_per = sqrt((double)(Sumsqr - Sum * Sum));
        if (!(std_per > 0.0)) std_per = 1.0;
      }
    } while (clipiter && nclippedthis > nclippedlast);
  }

  *avePower = ave_per;
  *rmsPower = std_per;

  /* Optional periodogram dump. */
  if (outflag) {
    if ((outfile = fopen(outname, "w")) == NULL) {
      fprintf(stderr, "Cannot write periodogram to %s\n", outname);
      exit(3);
    }
    if (ascii) {
      fprintf(outfile, "#Period FTP_Power FTP_Theta FTP_NegAmp\n");
      for (i = 0; i < Nperiod; i++) {
        if (periodogram[i] >= 0.0 && periodogram[i] * 0.0 == 0.0)
          fprintf(outfile, "%f %f %f %d\n",
                  periods[i], periodogram[i], theta_grid[i], negamp_grid[i]);
      }
    } else {
      fwrite(&Nperiod, 4, 1, outfile);
      fwrite(&ave_per, 8, 1, outfile);
      fwrite(&std_per, 8, 1, outfile);
      fwrite(periods,     8, Nperiod, outfile);
      fwrite(periodogram, 8, Nperiod, outfile);
    }
    fclose(outfile);
  }

  /* Reduce to local maxima (P spikes ABOVE neighbours). */
  {
    int kk;
    int kept = 0;
    lastpoint = periodogram[0] + 1.0;
    for (kk = 0; kk < Nperiod - 1; kk++) {
      if (periodogram[kk] > lastpoint && periodogram[kk] > periodogram[kk + 1]) {
        lastpoint = periodogram[kk];
        periodogram[kept] = periodogram[kk];
        periods[kept] = periods[kk];
        negamp_grid[kept] = negamp_grid[kk];
        theta_grid[kept]  = theta_grid[kk];
        kept++;
      } else {
        lastpoint = periodogram[kk];
      }
    }
    if (periodogram[Nperiod - 1] > lastpoint) {
      periodogram[kept] = periodogram[Nperiod - 1];
      periods[kept]     = periods[Nperiod - 1];
      negamp_grid[kept] = negamp_grid[Nperiod - 1];
      theta_grid[kept]  = theta_grid[Nperiod - 1];
      kept++;
    }
    Nperiod = kept;
  }

  /* Pick Npeaks highest-P distinct maxima (mirrors -aov / -PDM logic, with
   * P-MAXIMA semantics).  Two passes:
   *   (1) walk local maxima in period order and keep the first Npeaks distinct
   *       periods (a "different" period means not a multiple p*a/b within
   *       resolution; see ftp_isDifferentPeriods);
   *   (2) sort what we have descending in P, then continue the walk for any
   *       remaining period that beats the worst-kept slot.
   *
   * Without pass 2, a weak low-period local max can shadow the actual best
   * peak at higher period.  Pass 2 fixes that. */
  foundsofar = 0;
  i = 0;
  while (foundsofar < Npeaks && i < Nperiod) {
    if (periodogram[i] >= 0.0 && periodogram[i] * 0.0 == 0.0) {
      test = 1;
      for (j = 0; j < foundsofar; j++) {
        if (!ftp_isDifferentPeriods(perpeaks[j], periods[i], T)) {
          if (periodogram[i] > Ppeaks[j]) {
            perpeaks[j]   = periods[i];
            Ppeaks[j]     = periodogram[i];
            peakNegAmp[j] = negamp_grid[i];
            peakTheta[j]  = theta_grid[i];
          }
          test = 0;
          break;
        }
      }
      if (test) {
        perpeaks[foundsofar]   = periods[i];
        Ppeaks[foundsofar]     = periodogram[i];
        peakNegAmp[foundsofar] = negamp_grid[i];
        peakTheta[foundsofar]  = theta_grid[i];
        foundsofar++;
      }
    }
    i++;
  }
  for (k = foundsofar; k < Npeaks; k++) {
    perpeaks[k]   = -1.0;
    Ppeaks[k]     = -1.0;
    peakNegAmp[k] = 0;
    peakTheta[k]  = 0.0;
  }

  /* Sort descending in P (selection sort over Npeaks; lockstep across the
   * four parallel arrays). */
  {
    int n, m;
    for (n = 0; n < Npeaks - 1; n++) {
      int best_idx = n;
      for (m = n + 1; m < Npeaks; m++)
        if (Ppeaks[m] > Ppeaks[best_idx]) best_idx = m;
      if (best_idx != n) {
        double tp = perpeaks[n], tP = Ppeaks[n], tT = peakTheta[n];
        int    tA = peakNegAmp[n];
        perpeaks[n]   = perpeaks[best_idx];
        Ppeaks[n]     = Ppeaks[best_idx];
        peakTheta[n]  = peakTheta[best_idx];
        peakNegAmp[n] = peakNegAmp[best_idx];
        perpeaks[best_idx]   = tp;
        Ppeaks[best_idx]     = tP;
        peakTheta[best_idx]  = tT;
        peakNegAmp[best_idx] = tA;
      }
    }
  }

  /* Pass 2: continue walking the local-maxima list; any periodogram value
   * beating the worst-kept slot replaces it (subject to the distinctness
   * check), then re-sort. */
  {
    double worst_kept = (Ppeaks[Npeaks - 1] >= 0.0 ? Ppeaks[Npeaks - 1] : -1.0);
    for (; i < Nperiod; i++) {
      double v = periodogram[i];
      if (!(v >= 0.0 && v * 0.0 == 0.0)) continue;
      if (v <= worst_kept) continue;
      test = 1;
      for (j = 0; j < Npeaks; j++) {
        if (perpeaks[j] > 0.0 && !ftp_isDifferentPeriods(periods[i], perpeaks[j], T)) {
          if (v > Ppeaks[j]) {
            Ppeaks[j]     = v;
            perpeaks[j]   = periods[i];
            peakNegAmp[j] = negamp_grid[i];
            peakTheta[j]  = theta_grid[i];
          }
          test = 0;
          break;
        }
      }
      if (test) {
        /* Replace the worst-kept slot. */
        perpeaks[Npeaks - 1]   = periods[i];
        Ppeaks[Npeaks - 1]     = v;
        peakNegAmp[Npeaks - 1] = negamp_grid[i];
        peakTheta[Npeaks - 1]  = theta_grid[i];
      }
      /* Re-sort descending in P (lockstep). */
      {
        int n, m;
        for (n = 0; n < Npeaks - 1; n++) {
          int best_idx = n;
          for (m = n + 1; m < Npeaks; m++)
            if (Ppeaks[m] > Ppeaks[best_idx]) best_idx = m;
          if (best_idx != n) {
            double tp = perpeaks[n], tP = Ppeaks[n], tT = peakTheta[n];
            int    tA = peakNegAmp[n];
            perpeaks[n]   = perpeaks[best_idx];
            Ppeaks[n]     = Ppeaks[best_idx];
            peakTheta[n]  = peakTheta[best_idx];
            peakNegAmp[n] = peakNegAmp[best_idx];
            perpeaks[best_idx]   = tp;
            Ppeaks[best_idx]     = tP;
            peakTheta[best_idx]  = tT;
            peakNegAmp[best_idx] = tA;
          }
        }
      }
      worst_kept = (Ppeaks[Npeaks - 1] >= 0.0 ? Ppeaks[Npeaks - 1] : -1.0);
    }
  }

  /* High-resolution fine-tune around each peak. */
  smallfreqstep = finetune / T;
  for (j = 0; j < Npeaks; j++) {
    if (perpeaks[j] <= 0.0) continue;
    freq    = dmin((1.0 / perpeaks[j]) + freqstep, 1.0 / minP);
    minfreq = dmax((1.0 / perpeaks[j]) - freqstep, 1.0 / maxP);
    while (freq >= minfreq) {
      double th, P_test;
      int sgn;
      testperiod = 1.0 / freq;
      P_test = ftp_one_period(Nused, t, y_centered, w, YY, H, cn, sn,
                                allow_neg_amp, &sc, testperiod, &th, &sgn);
      if (P_test >= 0.0 && P_test * 0.0 == 0.0)
        if (P_test > Ppeaks[j]) {
          Ppeaks[j]     = P_test;
          perpeaks[j]   = testperiod;
          peakNegAmp[j] = (sgn < 0) ? 1 : 0;
          peakTheta[j]  = th;
        }
      freq -= smallfreqstep;
    }
    /* Period-multiple double-check. */
    bestscore = Ppeaks[j];
    ismultiple = 0;
    abest = bbest = 1;
    for (a = 1; a <= FTP_MAX_DOUBLE_CHECK_MULTIPLE; a++)
      for (b = 1; b <= FTP_MAX_DOUBLE_CHECK_MULTIPLE; b++)
        if (a != b) {
          double th, P_test;
          int sgn;
          testperiod = perpeaks[j] * a / b;
          if (testperiod > minP && testperiod < maxP) {
            P_test = ftp_one_period(Nused, t, y_centered, w, YY, H, cn, sn,
                                      allow_neg_amp, &sc, testperiod, &th, &sgn);
            if (P_test >= 0.0 && P_test * 0.0 == 0.0)
              if (P_test > bestscore) {
                ismultiple = 1; abest = a; bbest = b; bestscore = P_test;
              }
          }
        }
    if (ismultiple) {
      double th, P_test;
      int sgn;
      perpeaks[j] = perpeaks[j] * abest / bbest;
      P_test = ftp_one_period(Nused, t, y_centered, w, YY, H, cn, sn,
                                allow_neg_amp, &sc, perpeaks[j], &th, &sgn);
      Ppeaks[j]     = P_test;
      peakNegAmp[j] = (sgn < 0) ? 1 : 0;
      peakTheta[j]  = th;
    }
  }

  /* Final descending sort on Ppeaks. */
  {
    int n, m;
    /* O(Npeaks^2) selection sort, since Npeaks is small (~5) and we need
     * to permute multiple parallel arrays in lockstep. */
    for (n = 0; n < Npeaks - 1; n++) {
      int best_idx = n;
      for (m = n + 1; m < Npeaks; m++) {
        if (Ppeaks[m] > Ppeaks[best_idx]) best_idx = m;
      }
      if (best_idx != n) {
        double tp = perpeaks[n], tP = Ppeaks[n], tT = peakTheta[n];
        int    tA = peakNegAmp[n];
        perpeaks[n]   = perpeaks[best_idx];
        Ppeaks[n]     = Ppeaks[best_idx];
        peakTheta[n]  = peakTheta[best_idx];
        peakNegAmp[n] = peakNegAmp[best_idx];
        perpeaks[best_idx]   = tp;
        Ppeaks[best_idx]     = tP;
        peakTheta[best_idx]  = tT;
        peakNegAmp[best_idx] = tA;
      }
    }
  }

  /* SNR per peak (positive when peak rises above noise). */
  for (j = 0; j < Npeaks; j++) {
    if (perpeaks[j] > 0.0)
      peakSNR[j] = (Ppeaks[j] - ave_per) / std_per;
    else {
      perpeaks[j]   = -1.0;
      Ppeaks[j]     = 0.0;
      peakSNR[j]    = 0.0;
      peakNegAmp[j] = 0;
      peakTheta[j]  = 0.0;
    }
  }

  free(t); free(mag); free(w); free(y_centered);
  free(periods); free(periodogram);
  free(negamp_grid); free(theta_grid);
  ftp_scratch_free(&sc);
}

/* ---- Entry point invoked from processcommand.c ---- */

void RunFTPCommand(ProgramData *p, Command *c, _FTP *Ftp, int lcnum, int lc_name_num, int thisindex)
{
  char outname[MAXLEN];
  int i1, i2;
  (void) c; (void) thisindex;   /* not yet used; reserved for fixperiodSNR backref */

  if (Ftp->operiodogram) {
    i1 = 0; i2 = 0;
    while (p->lcnames[lc_name_num][i1] != '\0') {
      if (p->lcnames[lc_name_num][i1] == '/') i2 = i1 + 1;
      i1++;
    }
    sprintf(outname, "%s/%s%s", Ftp->outdir, &p->lcnames[lc_name_num][i2], Ftp->suffix);
  }

  /* Resolve var/expr/fixed parameter sources to per-LC scalar values. */
  if (Ftp->minp_source == VARTOOLS_SOURCE_EVALEXPRESSION)
    Ftp->minp_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Ftp->minp_expr);
  else if (Ftp->minp_source == VARTOOLS_SOURCE_EXISTINGVARIABLE)
    Ftp->minp_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Ftp->minp_var);
  else
    Ftp->minp_vals[lcnum] = Ftp->minp;

  if (Ftp->maxp_source == VARTOOLS_SOURCE_EVALEXPRESSION)
    Ftp->maxp_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Ftp->maxp_expr);
  else if (Ftp->maxp_source == VARTOOLS_SOURCE_EXISTINGVARIABLE)
    Ftp->maxp_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Ftp->maxp_var);
  else
    Ftp->maxp_vals[lcnum] = Ftp->maxp;

  if (Ftp->subsample_source == VARTOOLS_SOURCE_EVALEXPRESSION)
    Ftp->subsample_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Ftp->subsample_expr);
  else if (Ftp->subsample_source == VARTOOLS_SOURCE_EXISTINGVARIABLE)
    Ftp->subsample_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Ftp->subsample_var);
  else
    Ftp->subsample_vals[lcnum] = Ftp->subsample;

  if (Ftp->finetune_source == VARTOOLS_SOURCE_EVALEXPRESSION)
    Ftp->finetune_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Ftp->finetune_expr);
  else if (Ftp->finetune_source == VARTOOLS_SOURCE_EXISTINGVARIABLE)
    Ftp->finetune_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Ftp->finetune_var);
  else
    Ftp->finetune_vals[lcnum] = Ftp->finetune;

  findPeaks_ftp(p->t[lcnum], p->mag[lcnum], p->sig[lcnum], p->NJD[lcnum],
                Ftp->H, Ftp->cn, Ftp->sn,
                Ftp->useerr, Ftp->allow_neg_amp,
                Ftp->minp_vals[lcnum], Ftp->maxp_vals[lcnum],
                Ftp->subsample_vals[lcnum], Ftp->finetune_vals[lcnum],
                Ftp->Npeaks,
                Ftp->peakperiods[lcnum], Ftp->peakvalues[lcnum],
                Ftp->peakNegAmp[lcnum],
                Ftp->peakSNR[lcnum],
                Ftp->peakTheta[lcnum],
                &Ftp->avepower[lcnum], &Ftp->rmspower[lcnum],
                Ftp->operiodogram, outname, p->ascii,
                lcnum, lc_name_num);
}
