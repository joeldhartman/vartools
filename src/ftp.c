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

/* ---- ASCII light-curve reader (for the "fitlc" template-source mode) ----
 *
 * Reads three columns by 1-indexed position from a whitespace-separated
 * ASCII file.  Lines beginning with '#' and blank lines are skipped.
 *
 *   t_col, mag_col:  required (>= 1).
 *   err_col:         0 means "no error column, do an unweighted fit".
 *
 * Returns N on success (caller frees *t_out, *mag_out, *err_out); returns 0
 * on parse error or empty file. */
static int ftp_read_lc_ascii(const char *path,
                              int t_col, int mag_col, int err_col,
                              double **t_out, double **mag_out,
                              double **err_out, int *N_out)
{
  FILE *fp;
  char *line = NULL;
  size_t linelen = 0;
  ssize_t got;
  int N_alloc = 1024, N = 0;
  double *t = NULL, *mag = NULL, *err = NULL;
  int max_col;

  *t_out = *mag_out = *err_out = NULL;
  *N_out = 0;

  if (t_col < 1 || mag_col < 1) {
    fprintf(stderr, "-FTP fitlc ascii: t_col and mag_col must be >= 1 "
                    "(1-indexed); got t_col=%d, mag_col=%d\n", t_col, mag_col);
    return 0;
  }
  if (err_col < 0) err_col = 0;
  max_col = t_col;
  if (mag_col > max_col) max_col = mag_col;
  if (err_col > max_col) max_col = err_col;

  if ((fp = fopen(path, "r")) == NULL) {
    fprintf(stderr, "-FTP fitlc ascii: cannot open '%s' for reading\n", path);
    return 0;
  }

  t   = (double *) malloc(N_alloc * sizeof(double));
  mag = (double *) malloc(N_alloc * sizeof(double));
  if (err_col > 0) err = (double *) malloc(N_alloc * sizeof(double));
  if (t == NULL || mag == NULL || (err_col > 0 && err == NULL)) {
    fclose(fp);
    if (t) free(t); if (mag) free(mag); if (err) free(err);
    return 0;
  }

  while ((got = getline(&line, &linelen, fp)) != -1) {
    char *tok, *saveptr = NULL;
    int col = 1;
    double t_val = 0.0, mag_val = 0.0, err_val = 1.0;
    int got_t = 0, got_mag = 0, got_err = (err_col == 0 ? 1 : 0);
    char *q = line;
    while (*q && isspace((unsigned char)*q)) q++;
    if (*q == '\0' || *q == '#') continue;

    tok = strtok_r(line, " \t\n\r", &saveptr);
    while (tok != NULL && col <= max_col) {
      if (col == t_col)   { t_val   = atof(tok); got_t = 1; }
      if (col == mag_col) { mag_val = atof(tok); got_mag = 1; }
      if (err_col > 0 && col == err_col) { err_val = atof(tok); got_err = 1; }
      tok = strtok_r(NULL, " \t\n\r", &saveptr);
      col++;
    }
    if (!got_t || !got_mag || !got_err) continue;

    if (N >= N_alloc) {
      N_alloc *= 2;
      t   = (double *) realloc(t,   N_alloc * sizeof(double));
      mag = (double *) realloc(mag, N_alloc * sizeof(double));
      if (err) err = (double *) realloc(err, N_alloc * sizeof(double));
      if (t == NULL || mag == NULL || (err_col > 0 && err == NULL)) {
        if (line) free(line); fclose(fp);
        if (t) free(t); if (mag) free(mag); if (err) free(err);
        return 0;
      }
    }
    t[N]   = t_val;
    mag[N] = mag_val;
    if (err) err[N] = err_val;
    N++;
  }
  if (line) free(line);
  fclose(fp);

  if (N == 0) {
    fprintf(stderr, "-FTP fitlc ascii: no usable rows in '%s'\n", path);
    free(t); free(mag); if (err) free(err);
    return 0;
  }

  *t_out   = t;
  *mag_out = mag;
  *err_out = err;     /* NULL when err_col == 0 */
  *N_out   = N;
  return N;
}

/* ---- FITS light-curve reader (for the "fitlc" template-source mode) ----
 *
 * Reads three columns by name from a FITS binary-table extension.  Behaves
 * the same as ftp_read_lc_ascii on success.  When err_name is NULL, the
 * empty string, or the literal string "none" (case-insensitive), no error
 * column is read and *err_out is set to NULL (unweighted fit).
 *
 * Compiled out when vartools was built without cfitsio support; the
 * stub returns 0 with a clear error so the parser fails fast. */
#ifdef USECFITSIO
#include "fitsio.h"
#endif
static int ftp_read_lc_fits(const char *path,
                              const char *t_name, const char *mag_name,
                              const char *err_name,
                              double **t_out, double **mag_out,
                              double **err_out, int *N_out)
{
#ifdef USECFITSIO
  fitsfile *infile = NULL;
  int status = 0, hdutype = 0, hdunum = 0;
  long nrows = 0;
  int t_colnum = 0, mag_colnum = 0, err_colnum = 0;
  double *t = NULL, *mag = NULL, *err = NULL;
  int has_err;
  int anynul = 0;

  *t_out = *mag_out = *err_out = NULL;
  *N_out = 0;

  has_err = (err_name != NULL && err_name[0] != '\0' &&
             strcasecmp(err_name, "none") != 0);

  if (fits_open_file(&infile, path, READONLY, &status)) {
    fprintf(stderr, "-FTP fitlc fits: cannot open '%s' (cfitsio status %d)\n",
            path, status);
    return 0;
  }

  /* Move to the first table HDU if currently at the primary HDU. */
  if (fits_get_hdu_num(infile, &hdunum) == 1) {
    fits_movabs_hdu(infile, 2, &hdutype, &status);
  } else {
    fits_get_hdu_type(infile, &hdutype, &status);
  }
  if (hdutype == IMAGE_HDU) {
    fprintf(stderr, "-FTP fitlc fits: '%s' has an image HDU where a table was expected\n", path);
    fits_close_file(infile, &status);
    return 0;
  }

  fits_get_num_rows(infile, &nrows, &status);
  if (nrows <= 0) {
    fprintf(stderr, "-FTP fitlc fits: '%s' has zero rows\n", path);
    fits_close_file(infile, &status);
    return 0;
  }

  fits_get_colnum(infile, 0, (char *) t_name, &t_colnum, &status);
  if (status == COL_NOT_FOUND) {
    fprintf(stderr, "-FTP fitlc fits: column '%s' not found in '%s'\n", t_name, path);
    fits_close_file(infile, &status);
    return 0;
  }
  status = 0;
  fits_get_colnum(infile, 0, (char *) mag_name, &mag_colnum, &status);
  if (status == COL_NOT_FOUND) {
    fprintf(stderr, "-FTP fitlc fits: column '%s' not found in '%s'\n", mag_name, path);
    fits_close_file(infile, &status);
    return 0;
  }
  if (has_err) {
    status = 0;
    fits_get_colnum(infile, 0, (char *) err_name, &err_colnum, &status);
    if (status == COL_NOT_FOUND) {
      fprintf(stderr, "-FTP fitlc fits: column '%s' not found in '%s'\n", err_name, path);
      fits_close_file(infile, &status);
      return 0;
    }
  }
  status = 0;

  t   = (double *) malloc(nrows * sizeof(double));
  mag = (double *) malloc(nrows * sizeof(double));
  if (has_err) err = (double *) malloc(nrows * sizeof(double));
  if (t == NULL || mag == NULL || (has_err && err == NULL)) {
    if (t) free(t); if (mag) free(mag); if (err) free(err);
    fits_close_file(infile, &status);
    return 0;
  }

  fits_read_col(infile, TDOUBLE, t_colnum,   1, 1, nrows, NULL, t,   &anynul, &status);
  fits_read_col(infile, TDOUBLE, mag_colnum, 1, 1, nrows, NULL, mag, &anynul, &status);
  if (has_err)
    fits_read_col(infile, TDOUBLE, err_colnum, 1, 1, nrows, NULL, err, &anynul, &status);

  if (status) {
    fprintf(stderr, "-FTP fitlc fits: cfitsio read error (status=%d) on '%s'\n", status, path);
    free(t); free(mag); if (err) free(err);
    fits_close_file(infile, &status);
    return 0;
  }

  fits_close_file(infile, &status);

  *t_out   = t;
  *mag_out = mag;
  *err_out = err;
  *N_out   = (int) nrows;
  return (int) nrows;
#else
  (void) path; (void) t_name; (void) mag_name; (void) err_name;
  (void) t_out; (void) mag_out; (void) err_out; (void) N_out;
  fprintf(stderr, "-FTP fitlc fits: vartools was built without cfitsio "
                  "support; cannot read FITS files\n");
  return 0;
#endif
}

/* ---- Linear least-squares fit of a Fourier series to an LC ----
 *
 * Fits the model
 *     mag(t) ~= c0 + sum_{n=1..H_total} [c_n cos(n omega t) + s_n sin(n omega t)]
 * with omega = 2*pi/period and H_total = Nharm + 1, by solving the
 * (2*H_total + 1) x (2*H_total + 1) normal-equations system M^T W M theta
 * = M^T W y via gaussj.  Weights are w_i = 1/sigma_i^2 if err is non-NULL
 * and sigma_i > 0; uniform otherwise.
 *
 * Convention: "Nharm" follows -harmonicfilter / -Injectharm and counts
 * the harmonics ABOVE the fundamental.  Nharm = 0 fits only the
 * fundamental (one sin + one cos); Nharm = N fits the fundamental plus N
 * extra harmonics (N+1 sin + N+1 cos).  Total sinusoidal terms = 2*(Nharm+1).
 *
 * The constant term c0 is fitted (so the cosine/sine amplitudes are not
 * biased by the LC mean) but is then dropped: the FTP template is zero-mean
 * by construction.  After extraction, cn[] and sn[] are normalised so that
 * the half peak-to-peak amplitude of the model
 *     M(phi) = sum_n [c_n cos(2 pi n phi) + s_n sin(2 pi n phi)]
 * is 1, evaluated on a 1024-point uniform phase grid.
 *
 * Returns H_total = Nharm + 1 on success (writes *cn_out, *sn_out as
 * fresh allocations of length H_total); returns 0 on degenerate input or
 * solver failure. */
#define FTP_FITLC_NPHASE 1024
static int ftp_fit_harmonics(const double *t, const double *mag, const double *err,
                               int N, int Nharm, double period,
                               double **cn_out, double **sn_out)
{
  int p, i, j, k;
  int Nused = 0;
  int H_total;
  double omega;
  double *Adata = NULL, **A = NULL;
  double *Bdata = NULL, **B = NULL;
  double *row = NULL;
  double *cn = NULL, *sn = NULL;
  double model_min = 0.0, model_max = 0.0, half_amp;

  *cn_out = *sn_out = NULL;

  if (Nharm < 0) {
    fprintf(stderr, "-FTP fitlc: Nharm must be >= 0 (Nharm counts harmonics above the fundamental; "
                    "Nharm=0 fits only the fundamental); got %d\n", Nharm);
    return 0;
  }
  if (!(period > 0.0)) {
    fprintf(stderr, "-FTP fitlc: period must be > 0; got %g\n", period);
    return 0;
  }
  H_total = Nharm + 1;
  p = 2 * H_total + 1;
  omega = 2.0 * M_PI / period;

  /* Allocate normal-equations system, contiguous backing for cache locality. */
  Adata = (double *) calloc(p * p, sizeof(double));
  A     = (double **) malloc(p * sizeof(double *));
  Bdata = (double *) calloc(p * 1, sizeof(double));
  B     = (double **) malloc(p * sizeof(double *));
  row   = (double *) malloc(p * sizeof(double));
  if (Adata == NULL || A == NULL || Bdata == NULL || B == NULL || row == NULL)
    vt_error(ERR_MEMALLOC);
  for (i = 0; i < p; i++) { A[i] = Adata + i * p; B[i] = Bdata + i; }

  /* Accumulate. */
  for (i = 0; i < N; i++) {
    double wi, yi, x, ccurr, scurr, cprev, sprev, two_cos_x;
    if (isnan(mag[i])) continue;
    if (err != NULL) {
      if (!(err[i] > 0.0) || isnan(err[i])) continue;
      wi = 1.0 / (err[i] * err[i]);
    } else {
      wi = 1.0;
    }
    yi = mag[i];
    x  = omega * t[i];

    row[0] = 1.0;
    ccurr = cos(x); scurr = sin(x);
    row[1] = ccurr;
    row[2] = scurr;
    cprev = 1.0; sprev = 0.0;
    two_cos_x = 2.0 * ccurr;
    for (k = 2; k <= H_total; k++) {
      double cnext = two_cos_x * ccurr - cprev;
      double snext = two_cos_x * scurr - sprev;
      row[2*k - 1] = cnext;
      row[2*k]     = snext;
      cprev = ccurr; sprev = scurr;
      ccurr = cnext; scurr = snext;
    }

    for (j = 0; j < p; j++) {
      A[j][j] += wi * row[j] * row[j];
      for (k = j + 1; k < p; k++) {
        double v = wi * row[j] * row[k];
        A[j][k] += v;
        A[k][j] += v;
      }
      B[j][0] += wi * yi * row[j];
    }
    Nused++;
  }
  free(row);

  if (Nused < p) {
    fprintf(stderr, "-FTP fitlc: %d usable points but %d are needed for a "
                    "Nharm=%d fit (fundamental + %d extra harmonics, %d free parameters; "
                    "reduce Nharm or use a longer LC)\n",
            Nused, p, Nharm, Nharm, p);
    free(Adata); free(A); free(Bdata); free(B);
    return 0;
  }

  if (gaussj(A, p, B, 1)) {
    fprintf(stderr, "-FTP fitlc: linear-system solver (gaussj) failed; "
                    "the design matrix is singular -- check for "
                    "duplicate timestamps or an LC that does not span the period\n");
    free(Adata); free(A); free(Bdata); free(B);
    return 0;
  }

  cn = (double *) malloc(H_total * sizeof(double));
  sn = (double *) malloc(H_total * sizeof(double));
  if (cn == NULL || sn == NULL) vt_error(ERR_MEMALLOC);
  /* Drop c0 = B[0][0]; the FTP template is zero-mean by construction. */
  for (k = 0; k < H_total; k++) {
    cn[k] = B[2*k + 1][0];
    sn[k] = B[2*k + 2][0];
  }
  free(Adata); free(A); free(Bdata); free(B);

  /* Normalise so that the half peak-to-peak amplitude of the model is 1. */
  for (i = 0; i < FTP_FITLC_NPHASE; i++) {
    double phi = (double) i / (double) FTP_FITLC_NPHASE;
    double m_val = 0.0;
    double tau   = 2.0 * M_PI * phi;
    for (k = 1; k <= H_total; k++) {
      m_val += cn[k - 1] * cos(k * tau) + sn[k - 1] * sin(k * tau);
    }
    if (i == 0 || m_val < model_min) model_min = m_val;
    if (i == 0 || m_val > model_max) model_max = m_val;
  }
  half_amp = 0.5 * (model_max - model_min);
  if (!(half_amp > 0.0)) {
    fprintf(stderr, "-FTP fitlc: fitted model has zero peak-to-peak amplitude "
                    "(degenerate fit -- LC has no power at the specified period?)\n");
    free(cn); free(sn);
    return 0;
  }
  for (k = 0; k < H_total; k++) {
    cn[k] /= half_amp;
    sn[k] /= half_amp;
  }

  *cn_out = cn;
  *sn_out = sn;
  return H_total;
}

/* ---- Top-level fitlc orchestrator ----
 *
 * Loads an LC from disk, fits a Fourier series at the given period, drops
 * the mean and amplitude-normalises the result, and returns the cn/sn
 * arrays ready for use as an FTP template.  Format token is
 * "ascii" (1-indexed integer column numbers) or "fits" (string column
 * names; cfitsio backend).  Mirrors ftp_load_template_file's contract. */
int ftp_load_template_fitlc(const char *lc_path, const char *format,
                              const char *t_col_str,
                              const char *mag_col_str,
                              const char *err_col_str,
                              int Nharm, double period,
                              double **cn_out, double **sn_out)
{
  double *t = NULL, *mag = NULL, *err = NULL;
  int N = 0;
  int H_out;

  *cn_out = *sn_out = NULL;

  if (!strcmp(format, "ascii")) {
    int t_col   = atoi(t_col_str);
    int mag_col = atoi(mag_col_str);
    int err_col = atoi(err_col_str);     /* 0 means no err column */
    if (!ftp_read_lc_ascii(lc_path, t_col, mag_col, err_col,
                            &t, &mag, &err, &N))
      return 0;
  } else if (!strcmp(format, "fits")) {
    if (!ftp_read_lc_fits(lc_path, t_col_str, mag_col_str, err_col_str,
                           &t, &mag, &err, &N))
      return 0;
  } else {
    fprintf(stderr, "-FTP fitlc: format must be 'ascii' or 'fits'; got '%s'\n", format);
    return 0;
  }

  H_out = ftp_fit_harmonics(t, mag, err, N, Nharm, period, cn_out, sn_out);

  free(t); free(mag); if (err) free(err);
  return H_out;
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

/* =====================================================================
 *  Polynomial root-finding fast path (Phase B.4)
 *  ---------------------------------------------------------------------
 *  An alternative to ftp_max_P_over_theta's brute-force scan: at each
 *  frequency, build the optimality polynomial 2*MM*d(YM)/dtheta_2
 *  - YM*d(MM)/dtheta_2 = 0 as a polynomial in b = cos(theta_2), find
 *  its roots in [-1, 1] via gsl_poly_complex_solve, and evaluate P at
 *  each candidate (b, sgn(sin theta_2)) pair to pick the global max.
 *
 *  This is a port of the algorithm from Hoffman et al. (2021,
 *  arXiv:2101.12348), following the PseudoPolynomial scaffolding
 *  introduced by pyftp (https://github.com/PrincetonUniversity/FastTemplatePeriodogram).
 *  IMPORTANT: pyftp's compute_zeros and PseudoPolynomial.root_finding_poly
 *  compute (1-b^2)*p^2 - q^2 -- with p and q swapped relative to the
 *  correct squared form of "p(b) + sqrt(1-b^2)*q(b) = 0".  The right
 *  squared form is  p^2(b) - (1-b^2)*q^2(b) = 0,  which is what this
 *  code uses.  An empirical verification on H=1..5 synthetic LCs
 *  shows the corrected polynomial agrees with the brute-force theta_2
 *  scan to ~1e-7 (essentially the brute-force scan's finite-resolution
 *  noise floor), where pyftp's polynomial gives a 4-5 order-of-
 *  magnitude worse agreement (~1e-3 to 1e-2).  See the project root's
 *  FTP-rollout notes for the derivation and the reproduction script.
 *
 *  A PseudoPolynomial PP(b) = (1-b^2)^r * (p(b) + sqrt(1-b^2)*q(b))
 *  is represented here as separate p[], q[] coefficient arrays in
 *  ascending order, plus the integer power r (r <= 0).  Compared to
 *  pyftp's full class hierarchy this is intentionally minimal: we
 *  only need the operations multiply, derivative, and "squared root
 *  polynomial" -- so they're inlined as simple coefficient-array
 *  primitives rather than wrapped in a class.
 *  ===================================================================== */

/* --- Polynomial primitives over ascending-order coefficient arrays --- */

/* out[k] += in[k], up to min(out_cap, in_len). */
static void ftp_poly_addto(double *out, int out_cap, const double *in, int in_len)
{
  int k;
  int n = in_len < out_cap ? in_len : out_cap;
  for (k = 0; k < n; k++) out[k] += in[k];
}

/* out[k] -= in[k]. */
static void ftp_poly_subfrom(double *out, int out_cap, const double *in, int in_len)
{
  int k;
  int n = in_len < out_cap ? in_len : out_cap;
  for (k = 0; k < n; k++) out[k] -= in[k];
}

/* out = a * b.  out must have capacity >= a_len + b_len - 1.  out is
 * zeroed first.  a and b may have length 0 (returns zero polynomial). */
static void ftp_poly_mul(double *out, int out_cap,
                          const double *a, int a_len,
                          const double *b, int b_len)
{
  int i, j;
  int out_len = (a_len > 0 && b_len > 0) ? (a_len + b_len - 1) : 0;
  if (out_len > out_cap) out_len = out_cap;
  for (i = 0; i < out_len; i++) out[i] = 0.0;
  for (i = 0; i < a_len; i++) {
    double ai = a[i];
    if (ai == 0.0) continue;
    int jmax = b_len < out_cap - i ? b_len : out_cap - i;
    for (j = 0; j < jmax; j++) out[i + j] += ai * b[j];
  }
}

/* out = (1 - x^2) * in.  out must have capacity >= in_len + 2.
 * out is zeroed first. */
static void ftp_poly_mul_oneMinusXsq(double *out, int out_cap,
                                      const double *in, int in_len)
{
  int k;
  int out_len = in_len + 2;
  if (out_len > out_cap) out_len = out_cap;
  for (k = 0; k < out_len; k++) out[k] = 0.0;
  for (k = 0; k < in_len && k < out_cap; k++) out[k] += in[k];
  for (k = 0; k < in_len && k + 2 < out_cap; k++) out[k + 2] -= in[k];
}

/* out = x * in.  out must have capacity >= in_len + 1.  out is zeroed first. */
static void ftp_poly_mul_x(double *out, int out_cap, const double *in, int in_len)
{
  int k;
  int out_len = in_len + 1;
  if (out_len > out_cap) out_len = out_cap;
  out[0] = 0.0;
  for (k = 0; k < in_len && k + 1 < out_cap; k++) out[k + 1] = in[k];
}

/* out = d(in)/dx.  out_len = max(0, in_len - 1). */
static void ftp_poly_deriv(double *out, int out_cap, const double *in, int in_len)
{
  int k;
  int out_len = in_len > 0 ? in_len - 1 : 0;
  if (out_len > out_cap) out_len = out_cap;
  for (k = 0; k < out_len; k++) out[k] = (double)(k + 1) * in[k + 1];
}

/* --- PseudoPolynomial-aware helpers ---
 *
 * Each PP is a triple (p[], q[], r) stored as two coefficient arrays and
 * an integer.  No allocator: caller manages buffer capacity. */

/* Product of two PPs:
 *   out_p = a_p * b_p + (1-x^2) * a_q * b_q
 *   out_q = a_p * b_q + a_q * b_p
 *   out_r = a_r + b_r
 *
 * Buffer requirements:
 *   *out_p_len returns the new p length.
 *   *out_q_len returns the new q length.
 *   out_p_cap must be >= a_p_len + b_p_len - 1  AND  >= a_q_len + b_q_len + 1.
 *   out_q_cap must be >= max(a_p_len + b_q_len - 1, a_q_len + b_p_len - 1).
 * Scratch arrays of capacity >= max(a_p_len + b_p_len, a_q_len + b_q_len + 2)
 *   needed to avoid aliasing with out_p / out_q. */
static void ftp_pp_mul(double *out_p, int out_p_cap, int *out_p_len,
                        double *out_q, int out_q_cap, int *out_q_len,
                        int *out_r,
                        const double *a_p, int a_p_len,
                        const double *a_q, int a_q_len, int a_r,
                        const double *b_p, int b_p_len,
                        const double *b_q, int b_q_len, int b_r,
                        double *scratch1, int scratch1_cap,
                        double *scratch2, int scratch2_cap)
{
  int i;
  int pp_len, qq_len_raw, qq_len, pq_len, qp_len;
  /* out_p = a_p*b_p + (1-x^2)*(a_q*b_q) */
  ftp_poly_mul(out_p, out_p_cap, a_p, a_p_len, b_p, b_p_len);
  pp_len = (a_p_len > 0 && b_p_len > 0) ? (a_p_len + b_p_len - 1) : 0;
  if (pp_len > out_p_cap) pp_len = out_p_cap;
  /* scratch1 = a_q * b_q */
  ftp_poly_mul(scratch1, scratch1_cap, a_q, a_q_len, b_q, b_q_len);
  qq_len_raw = (a_q_len > 0 && b_q_len > 0) ? (a_q_len + b_q_len - 1) : 0;
  if (qq_len_raw > scratch1_cap) qq_len_raw = scratch1_cap;
  /* scratch2 = (1-x^2) * scratch1 */
  ftp_poly_mul_oneMinusXsq(scratch2, scratch2_cap, scratch1, qq_len_raw);
  qq_len = qq_len_raw + 2;
  if (qq_len > scratch2_cap) qq_len = scratch2_cap;
  /* out_p += scratch2 */
  ftp_poly_addto(out_p, out_p_cap, scratch2, qq_len);
  *out_p_len = pp_len > qq_len ? pp_len : qq_len;
  if (*out_p_len > out_p_cap) *out_p_len = out_p_cap;

  /* out_q = a_p*b_q + a_q*b_p */
  ftp_poly_mul(out_q, out_q_cap, a_p, a_p_len, b_q, b_q_len);
  pq_len = (a_p_len > 0 && b_q_len > 0) ? (a_p_len + b_q_len - 1) : 0;
  if (pq_len > out_q_cap) pq_len = out_q_cap;
  ftp_poly_mul(scratch1, scratch1_cap, a_q, a_q_len, b_p, b_p_len);
  qp_len = (a_q_len > 0 && b_p_len > 0) ? (a_q_len + b_p_len - 1) : 0;
  if (qp_len > scratch1_cap) qp_len = scratch1_cap;
  ftp_poly_addto(out_q, out_q_cap, scratch1, qp_len);
  *out_q_len = pq_len > qp_len ? pq_len : qp_len;
  if (*out_q_len > out_q_cap) *out_q_len = out_q_cap;

  *out_r = a_r + b_r;
  (void) i;
}

/* Derivative of a PP:
 *   out_p = (1-x^2) * dp/dx - 2r * x * p
 *   out_q = (1-x^2) * dq/dx - (2r+1) * x * q
 *   out_r = r - 1
 *
 * Buffer requirements:
 *   out_p_cap >= max(in_p_len + 1, in_p_len + 1)  ==> in_p_len + 1 (since dp has deg in_p_len-2)
 *     (1-x^2)*dp has degree in_p_len, so length in_p_len + 1.
 *     x*p has degree in_p_len, so length in_p_len + 1.
 *     out_p has length max of these = in_p_len + 1.
 *   out_q_cap >= in_q_len + 1.
 * Scratch needs capacity >= max(in_p_len + 1, in_q_len + 1). */
static void ftp_pp_deriv(double *out_p, int out_p_cap, int *out_p_len,
                          double *out_q, int out_q_cap, int *out_q_len,
                          int *out_r,
                          const double *in_p, int in_p_len,
                          const double *in_q, int in_q_len, int in_r,
                          double *scratch, int scratch_cap)
{
  int i;
  int dp_len, ddp_len, xp_len;
  int dq_len, ddq_len, xq_len;
  for (i = 0; i < out_p_cap; i++) out_p[i] = 0.0;
  for (i = 0; i < out_q_cap; i++) out_q[i] = 0.0;

  /* out_p part 1: (1-x^2) * dp/dx */
  ftp_poly_deriv(scratch, scratch_cap, in_p, in_p_len);
  dp_len = in_p_len > 0 ? in_p_len - 1 : 0;
  if (dp_len > scratch_cap) dp_len = scratch_cap;
  ftp_poly_mul_oneMinusXsq(out_p, out_p_cap, scratch, dp_len);
  ddp_len = dp_len + 2;
  if (ddp_len > out_p_cap) ddp_len = out_p_cap;
  /* out_p part 2: -2r * x * p */
  if (in_r != 0) {
    ftp_poly_mul_x(scratch, scratch_cap, in_p, in_p_len);
    xp_len = in_p_len + 1;
    if (xp_len > scratch_cap) xp_len = scratch_cap;
    for (i = 0; i < xp_len && i < out_p_cap; i++)
      out_p[i] -= 2.0 * (double) in_r * scratch[i];
    *out_p_len = ddp_len > xp_len ? ddp_len : xp_len;
  } else {
    *out_p_len = ddp_len;
  }
  if (*out_p_len > out_p_cap) *out_p_len = out_p_cap;

  /* out_q part 1: (1-x^2) * dq/dx */
  ftp_poly_deriv(scratch, scratch_cap, in_q, in_q_len);
  dq_len = in_q_len > 0 ? in_q_len - 1 : 0;
  if (dq_len > scratch_cap) dq_len = scratch_cap;
  ftp_poly_mul_oneMinusXsq(out_q, out_q_cap, scratch, dq_len);
  ddq_len = dq_len + 2;
  if (ddq_len > out_q_cap) ddq_len = out_q_cap;
  /* out_q part 2: -(2r+1) * x * q */
  {
    double coeff = -(double)(2 * in_r + 1);
    if (coeff != 0.0) {
      ftp_poly_mul_x(scratch, scratch_cap, in_q, in_q_len);
      xq_len = in_q_len + 1;
      if (xq_len > scratch_cap) xq_len = scratch_cap;
      for (i = 0; i < xq_len && i < out_q_cap; i++)
        out_q[i] += coeff * scratch[i];
      *out_q_len = ddq_len > xq_len ? ddq_len : xq_len;
    } else {
      *out_q_len = ddq_len;
    }
  }
  if (*out_q_len > out_q_cap) *out_q_len = out_q_cap;

  *out_r = in_r - 1;
}

/* ---- Build A_n / B_n / dA_n / dB_n as PseudoPolynomials from cn, sn ----
 *
 * For n = 1..H (stored at index n-1):
 *   A_n(b) = c_n * T_n(b) - sgn * s_n * U_{n-1}(b) * sqrt(1-b^2)
 *          = (c_n*T_n)[p] + (-sgn*s_n*U_{n-1})[q]
 *   B_n(b) = s_n * T_n(b) + sgn * c_n * U_{n-1}(b) * sqrt(1-b^2)
 *          = (s_n*T_n)[p] + ( sgn*c_n*U_{n-1})[q]
 *
 * Convention: this code builds the PPs for sgn = +1.  The sgn = -1 case
 * has the q components negated, which we handle at evaluation time by
 * flipping sin(theta_2)'s sign in ftp_eval_P.  Storing only the sgn=+1
 * tensors halves the memory and the polynomial-coefficient work per
 * frequency (mirrors pyftp's convention; the polynomial roots are the
 * same up to a sign flip in q anyway, which doesn't affect b values).
 *
 * Each output PP is stored in caller-provided buffers (p_buf, q_buf)
 * of capacity >= H + 2 per PP.  Returns the per-n p_len and q_len in
 * the per-n length arrays.  All output buffers are owned by the
 * caller. */
static void ftp_chebyshev_T_coefs(int n, double *out, int out_cap, int *out_len)
{
  /* Chebyshev T_n(x) coefficients in ascending order.
   * T_0=1, T_1=x, T_{n+1}=2x*T_n - T_{n-1}. */
  int i, k;
  double *Tprev, *Tcurr, *Tnext, *tmp;
  static double *workA = NULL;
  static double *workB = NULL;
  static double *workC = NULL;
  static int work_cap = 0;
  if (n < 0) { *out_len = 0; return; }
  if (work_cap < n + 2) {
    free(workA); free(workB); free(workC);
    work_cap = n + 8;
    workA = (double *) calloc(work_cap, sizeof(double));
    workB = (double *) calloc(work_cap, sizeof(double));
    workC = (double *) calloc(work_cap, sizeof(double));
    if (workA == NULL || workB == NULL || workC == NULL) vt_error(ERR_MEMALLOC);
  }
  Tprev = workA; Tcurr = workB; Tnext = workC;
  for (i = 0; i < work_cap; i++) { Tprev[i] = 0.0; Tcurr[i] = 0.0; }
  Tprev[0] = 1.0;            /* T_0 = 1 */
  if (n == 0) {
    *out_len = (out_cap >= 1) ? 1 : 0;
    if (*out_len) out[0] = 1.0;
    return;
  }
  Tcurr[1] = 1.0;            /* T_1 = x */
  if (n == 1) {
    *out_len = (out_cap >= 2) ? 2 : out_cap;
    for (i = 0; i < *out_len; i++) out[i] = Tcurr[i];
    return;
  }
  for (k = 2; k <= n; k++) {
    /* T_k = 2x*T_{k-1} - T_{k-2} */
    for (i = 0; i < work_cap; i++) Tnext[i] = -Tprev[i];
    for (i = 0; i + 1 < work_cap; i++) Tnext[i + 1] += 2.0 * Tcurr[i];
    /* rotate: prev <- curr, curr <- next */
    tmp = Tprev; Tprev = Tcurr; Tcurr = Tnext; Tnext = tmp;
  }
  *out_len = (n + 1 < out_cap) ? (n + 1) : out_cap;
  for (i = 0; i < *out_len; i++) out[i] = Tcurr[i];
}

static void ftp_chebyshev_U_coefs(int n, double *out, int out_cap, int *out_len)
{
  /* Chebyshev U_n(x) coefficients in ascending order.
   * U_0=1, U_1=2x, U_{n+1}=2x*U_n - U_{n-1}.  Same recurrence shape as T. */
  int i, k;
  double *Uprev, *Ucurr, *Unext, *tmp;
  static double *workA = NULL;
  static double *workB = NULL;
  static double *workC = NULL;
  static int work_cap = 0;
  if (n < 0) { *out_len = 0; return; }
  if (work_cap < n + 2) {
    free(workA); free(workB); free(workC);
    work_cap = n + 8;
    workA = (double *) calloc(work_cap, sizeof(double));
    workB = (double *) calloc(work_cap, sizeof(double));
    workC = (double *) calloc(work_cap, sizeof(double));
    if (workA == NULL || workB == NULL || workC == NULL) vt_error(ERR_MEMALLOC);
  }
  Uprev = workA; Ucurr = workB; Unext = workC;
  for (i = 0; i < work_cap; i++) { Uprev[i] = 0.0; Ucurr[i] = 0.0; }
  Uprev[0] = 1.0;            /* U_0 = 1 */
  if (n == 0) {
    *out_len = (out_cap >= 1) ? 1 : 0;
    if (*out_len) out[0] = 1.0;
    return;
  }
  Ucurr[1] = 2.0;            /* U_1 = 2x */
  if (n == 1) {
    *out_len = (out_cap >= 2) ? 2 : out_cap;
    for (i = 0; i < *out_len; i++) out[i] = Ucurr[i];
    return;
  }
  for (k = 2; k <= n; k++) {
    for (i = 0; i < work_cap; i++) Unext[i] = -Uprev[i];
    for (i = 0; i + 1 < work_cap; i++) Unext[i + 1] += 2.0 * Ucurr[i];
    tmp = Uprev; Uprev = Ucurr; Ucurr = Unext; Unext = tmp;
  }
  *out_len = (n + 1 < out_cap) ? (n + 1) : out_cap;
  for (i = 0; i < *out_len; i++) out[i] = Ucurr[i];
}

/* ---- Per-LC polynomial path: build per-frequency Pp[L], Pq[L], find
 * roots of Pp^2 - (1-b^2)*Pq^2 = 0, evaluate P at each root + sgn=+/-1,
 * return max. ----
 *
 * Caller pre-allocates and reuses scratch via _FTPScratch + the
 * polynomial-specific extension below.  For Phase B.4.a the polynomial
 * path is used only by ftp_periodogram_method when method == FTP_METHOD_POLY.
 */

#include <gsl/gsl_poly.h>
#include <gsl/gsl_errno.h>

#ifdef HAVE_NFFT3
#include <nfft3.h>
#endif

#define FTP_METHOD_BRUTE  0
#define FTP_METHOD_POLY   1
#define FTP_METHOD_VERIFY 2     /* compute both, compare, return poly result */
#define FTP_VERIFY_TOL    1.0e-3 /* max allowed |P_poly - P_brute| before warning */

/* Polynomial-mode scratch.  Per-frequency this gets filled with the
 * polynomial coefficient arrays, then root-found.  All buffers are sized
 * for max H given at allocation; H must not exceed it on subsequent
 * frequencies. */
typedef struct {
  int H;
  int p_cap;          /* per-PP-buffer capacity, sized for one A_i*A_j*dA_k product */
  int L_cap;          /* root-polynomial length capacity = 6H + 4 (with margin) */

  /* Per-template (built once per template): A_n, B_n, dA_n, dB_n.
   * Each is an array of H PseudoPolynomials.  Stored as flat buffers
   * with per-n p_len/q_len recording the active length. */
  double *A_p, *A_q;          /* H rows of length p_cap each; row n at offset n*p_cap */
  double *B_p, *B_q;
  double *dA_p, *dA_q;
  double *dB_p, *dB_q;
  int    *A_p_len, *A_q_len;  /* size H */
  int    *B_p_len, *B_q_len;
  int    *dA_p_len, *dA_q_len;
  int    *dB_p_len, *dB_q_len;
  /* r values: A/B have r=0, dA/dB have r=-1 (built that way). */

  /* Per-frequency scratch for products and accumulators. */
  double *prod_p, *prod_q;    /* intermediate (X*Y).p, .q */
  int     prod_p_len, prod_q_len, prod_r;
  double *triple_p, *triple_q; /* (X*Y*Z).p, .q -- the three-PP product */
  int     triple_p_len, triple_q_len, triple_r;
  double *Pp, *Pq;            /* accumulated Pp[L], Pq[L] -- size L_cap */
  int     Pp_len, Pq_len;
  double *root_poly;          /* Pp^2 - (1-b^2)*Pq^2, size L_cap */
  int     root_poly_len;
  double *scratch1, *scratch2;
  int     sc_cap;

  /* Precomputed template tensors -- built once at template-load time so
   * the per-frequency cost drops from O(H^5) (with per-frequency PP
   * triple-products) to O(H^4) (with K-tensor builds + contractions).
   * Each tensor is shape (H, H, H, L_T), indexed
   *    T[i*H*H*L_T + j*H*L_T + k*L_T + l]
   * with i, j, k in [0..H-1] and l in [0..L_T-1].  L_T = 3H+4 (covers
   * max triple-product length 3H+2 with some margin).
   *
   *   T_AAdAp[i,j,k,l]: l-th coef of (A_i * A_j * dA_k).p
   *   T_AAdAq[i,j,k,l]: l-th coef of (A_i * A_j * dA_k).q
   *   (12 tensors total: AAdA, AAdB, ABdA, ABdB, BBdA, BBdB; each with p, q.)
   */
  int     L_T;                /* per-tensor polynomial length capacity */
  double *T_AAdAp; double *T_AAdAq;
  double *T_AAdBp; double *T_AAdBq;
  double *T_ABdAp; double *T_ABdAq;
  double *T_ABdBp; double *T_ABdBq;
  double *T_BBdAp; double *T_BBdAq;
  double *T_BBdBp; double *T_BBdBq;
  /* Per-frequency K-tensor scratch (H*H*H each). */
  double *K_aada; double *K_aadb; double *K_abda;
  double *K_abdb; double *K_bbda; double *K_bbdb;

  /* Root tracking across frequencies.  At nearby frequencies, the roots
   * of the optimality polynomial move smoothly.  Instead of running a
   * full gsl_poly_complex_solve at every frequency, we Newton-refine
   * the previous frequency's roots against the current polynomial.  A
   * full eigensolve runs at frequency 0 and periodically thereafter
   * (every full_every frequencies) to recover any roots that drifted
   * out of the [-1, 1] interval or that emerged from complex collisions. */
  int     prev_roots_valid;    /* 1 if prev_real_roots has been populated */
  int     prev_n_real;         /* number of real roots in [-1, 1] from last freq */
  double *prev_real_roots;     /* size L_cap, real-root b values */
  int     freqs_since_full;    /* counter; reset to 0 when we full-eigensolve */
  int     full_every;          /* full eigensolve every N frequencies */
  double *Fprime_coefs;        /* derivative of root_poly, size L_cap */

  /* gsl_poly_complex_solve output buffer: 2*(L_cap-1) doubles for re/im pairs. */
  double *roots_packed;
  gsl_poly_complex_workspace *gsl_ws;
  int     gsl_ws_cap;
} _FTPPolyState;

static void ftp_poly_state_alloc(_FTPPolyState *ps, int H)
{
  int p_cap, L_cap, sc_cap;
  /* For a triple product A_i*A_j*dA_k of PPs of degrees ~H, the result
   * has p of degree up to 3H+1 and q of degree up to 3H.  Allocate
   * generously: p_cap = 3H+8 covers any intermediate. */
  p_cap = 3 * H + 8;
  L_cap = 6 * H + 8;
  sc_cap = L_cap;

  memset(ps, 0, sizeof(*ps));
  ps->H = H;
  ps->p_cap = p_cap;
  ps->L_cap = L_cap;
  ps->sc_cap = sc_cap;

  ps->A_p    = (double *) calloc(H * p_cap, sizeof(double));
  ps->A_q    = (double *) calloc(H * p_cap, sizeof(double));
  ps->B_p    = (double *) calloc(H * p_cap, sizeof(double));
  ps->B_q    = (double *) calloc(H * p_cap, sizeof(double));
  ps->dA_p   = (double *) calloc(H * p_cap, sizeof(double));
  ps->dA_q   = (double *) calloc(H * p_cap, sizeof(double));
  ps->dB_p   = (double *) calloc(H * p_cap, sizeof(double));
  ps->dB_q   = (double *) calloc(H * p_cap, sizeof(double));
  ps->A_p_len = (int *) calloc(H, sizeof(int));
  ps->A_q_len = (int *) calloc(H, sizeof(int));
  ps->B_p_len = (int *) calloc(H, sizeof(int));
  ps->B_q_len = (int *) calloc(H, sizeof(int));
  ps->dA_p_len = (int *) calloc(H, sizeof(int));
  ps->dA_q_len = (int *) calloc(H, sizeof(int));
  ps->dB_p_len = (int *) calloc(H, sizeof(int));
  ps->dB_q_len = (int *) calloc(H, sizeof(int));

  ps->prod_p   = (double *) calloc(p_cap, sizeof(double));
  ps->prod_q   = (double *) calloc(p_cap, sizeof(double));
  ps->triple_p = (double *) calloc(p_cap, sizeof(double));
  ps->triple_q = (double *) calloc(p_cap, sizeof(double));
  ps->Pp       = (double *) calloc(L_cap, sizeof(double));
  ps->Pq       = (double *) calloc(L_cap, sizeof(double));
  ps->root_poly = (double *) calloc(L_cap, sizeof(double));
  ps->scratch1 = (double *) calloc(sc_cap, sizeof(double));
  ps->scratch2 = (double *) calloc(sc_cap, sizeof(double));

  ps->gsl_ws_cap = L_cap;
  ps->roots_packed = (double *) calloc(2 * (L_cap - 1), sizeof(double));
  ps->gsl_ws = gsl_poly_complex_workspace_alloc((size_t) L_cap);
  if (ps->gsl_ws == NULL) vt_error(ERR_MEMALLOC);

  /* Precomputed template tensors -- 12 of shape (H, H, H, L_T). */
  ps->L_T = 3 * H + 4;
  {
    int tsize = H * H * H * ps->L_T;
    ps->T_AAdAp = (double *) calloc(tsize, sizeof(double));
    ps->T_AAdAq = (double *) calloc(tsize, sizeof(double));
    ps->T_AAdBp = (double *) calloc(tsize, sizeof(double));
    ps->T_AAdBq = (double *) calloc(tsize, sizeof(double));
    ps->T_ABdAp = (double *) calloc(tsize, sizeof(double));
    ps->T_ABdAq = (double *) calloc(tsize, sizeof(double));
    ps->T_ABdBp = (double *) calloc(tsize, sizeof(double));
    ps->T_ABdBq = (double *) calloc(tsize, sizeof(double));
    ps->T_BBdAp = (double *) calloc(tsize, sizeof(double));
    ps->T_BBdAq = (double *) calloc(tsize, sizeof(double));
    ps->T_BBdBp = (double *) calloc(tsize, sizeof(double));
    ps->T_BBdBq = (double *) calloc(tsize, sizeof(double));
    if (ps->T_AAdAp == NULL || ps->T_AAdAq == NULL ||
        ps->T_AAdBp == NULL || ps->T_AAdBq == NULL ||
        ps->T_ABdAp == NULL || ps->T_ABdAq == NULL ||
        ps->T_ABdBp == NULL || ps->T_ABdBq == NULL ||
        ps->T_BBdAp == NULL || ps->T_BBdAq == NULL ||
        ps->T_BBdBp == NULL || ps->T_BBdBq == NULL)
      vt_error(ERR_MEMALLOC);
  }
  {
    int Ksize = H * H * H;
    ps->K_aada = (double *) calloc(Ksize, sizeof(double));
    ps->K_aadb = (double *) calloc(Ksize, sizeof(double));
    ps->K_abda = (double *) calloc(Ksize, sizeof(double));
    ps->K_abdb = (double *) calloc(Ksize, sizeof(double));
    ps->K_bbda = (double *) calloc(Ksize, sizeof(double));
    ps->K_bbdb = (double *) calloc(Ksize, sizeof(double));
    if (ps->K_aada == NULL || ps->K_aadb == NULL ||
        ps->K_abda == NULL || ps->K_abdb == NULL ||
        ps->K_bbda == NULL || ps->K_bbdb == NULL)
      vt_error(ERR_MEMALLOC);
  }

  /* Root-tracking state. */
  ps->prev_roots_valid = 0;
  ps->prev_n_real = 0;
  ps->prev_real_roots = (double *) calloc(L_cap, sizeof(double));
  ps->Fprime_coefs    = (double *) calloc(L_cap, sizeof(double));
  ps->freqs_since_full = 0;
  ps->full_every       = 50;        /* re-eigensolve every 50 frequencies */
  if (ps->prev_real_roots == NULL || ps->Fprime_coefs == NULL)
    vt_error(ERR_MEMALLOC);
}

static void ftp_poly_state_free(_FTPPolyState *ps)
{
  free(ps->A_p); free(ps->A_q); free(ps->B_p); free(ps->B_q);
  free(ps->dA_p); free(ps->dA_q); free(ps->dB_p); free(ps->dB_q);
  free(ps->A_p_len); free(ps->A_q_len); free(ps->B_p_len); free(ps->B_q_len);
  free(ps->dA_p_len); free(ps->dA_q_len); free(ps->dB_p_len); free(ps->dB_q_len);
  free(ps->prod_p); free(ps->prod_q);
  free(ps->triple_p); free(ps->triple_q);
  free(ps->Pp); free(ps->Pq);
  free(ps->root_poly);
  free(ps->scratch1); free(ps->scratch2);
  free(ps->roots_packed);
  if (ps->gsl_ws) gsl_poly_complex_workspace_free(ps->gsl_ws);
  free(ps->T_AAdAp); free(ps->T_AAdAq);
  free(ps->T_AAdBp); free(ps->T_AAdBq);
  free(ps->T_ABdAp); free(ps->T_ABdAq);
  free(ps->T_ABdBp); free(ps->T_ABdBq);
  free(ps->T_BBdAp); free(ps->T_BBdAq);
  free(ps->T_BBdBp); free(ps->T_BBdBq);
  free(ps->K_aada); free(ps->K_aadb); free(ps->K_abda);
  free(ps->K_abdb); free(ps->K_bbda); free(ps->K_bbdb);
  free(ps->prev_real_roots); free(ps->Fprime_coefs);
  memset(ps, 0, sizeof(*ps));
}

/* Build all A_n, B_n PPs (r = 0) and their derivatives dA_n, dB_n (r = -1)
 * for n = 1..H, sgn = +1 fixed.  Done once per template. */
static void ftp_poly_build_templates(_FTPPolyState *ps,
                                       const double *cn, const double *sn)
{
  int n, i;
  int H = ps->H;
  int p_cap = ps->p_cap;
  double *Tn_buf, *Un_buf;
  int Tn_len, Un_len;

  Tn_buf = (double *) calloc(H + 2, sizeof(double));
  Un_buf = (double *) calloc(H + 2, sizeof(double));
  if (Tn_buf == NULL || Un_buf == NULL) vt_error(ERR_MEMALLOC);

  /* A_n / B_n */
  for (n = 1; n <= H; n++) {
    ftp_chebyshev_T_coefs(n,     Tn_buf, H + 2, &Tn_len);
    ftp_chebyshev_U_coefs(n - 1, Un_buf, H + 2, &Un_len);
    /* A_n: p = c_n * T_n, q = -s_n * U_{n-1}  (sgn = +1) */
    for (i = 0; i < p_cap; i++) ps->A_p[(n - 1) * p_cap + i] = 0.0;
    for (i = 0; i < p_cap; i++) ps->A_q[(n - 1) * p_cap + i] = 0.0;
    for (i = 0; i < Tn_len; i++) ps->A_p[(n - 1) * p_cap + i] = cn[n - 1] * Tn_buf[i];
    for (i = 0; i < Un_len; i++) ps->A_q[(n - 1) * p_cap + i] = -sn[n - 1] * Un_buf[i];
    ps->A_p_len[n - 1] = Tn_len;
    ps->A_q_len[n - 1] = Un_len;
    /* B_n: p = s_n * T_n, q =  c_n * U_{n-1}  (sgn = +1) */
    for (i = 0; i < p_cap; i++) ps->B_p[(n - 1) * p_cap + i] = 0.0;
    for (i = 0; i < p_cap; i++) ps->B_q[(n - 1) * p_cap + i] = 0.0;
    for (i = 0; i < Tn_len; i++) ps->B_p[(n - 1) * p_cap + i] = sn[n - 1] * Tn_buf[i];
    for (i = 0; i < Un_len; i++) ps->B_q[(n - 1) * p_cap + i] =  cn[n - 1] * Un_buf[i];
    ps->B_p_len[n - 1] = Tn_len;
    ps->B_q_len[n - 1] = Un_len;
  }

  /* dA_n, dB_n */
  for (n = 1; n <= H; n++) {
    int out_r;
    ftp_pp_deriv(&ps->dA_p[(n - 1) * p_cap], p_cap, &ps->dA_p_len[n - 1],
                  &ps->dA_q[(n - 1) * p_cap], p_cap, &ps->dA_q_len[n - 1],
                  &out_r,
                  &ps->A_p[(n - 1) * p_cap], ps->A_p_len[n - 1],
                  &ps->A_q[(n - 1) * p_cap], ps->A_q_len[n - 1], 0,
                  ps->scratch1, ps->sc_cap);
    ftp_pp_deriv(&ps->dB_p[(n - 1) * p_cap], p_cap, &ps->dB_p_len[n - 1],
                  &ps->dB_q[(n - 1) * p_cap], p_cap, &ps->dB_q_len[n - 1],
                  &out_r,
                  &ps->B_p[(n - 1) * p_cap], ps->B_p_len[n - 1],
                  &ps->B_q[(n - 1) * p_cap], ps->B_q_len[n - 1], 0,
                  ps->scratch1, ps->sc_cap);
    /* dA/dB have r = -1 (verified by out_r). */
  }

  free(Tn_buf); free(Un_buf);

  /* Fill the 12 template tensors: T_XYdZp[i,j,k,l] and T_XYdZq[i,j,k,l]
   * are the p and q coefficients of the PP triple product X_i*Y_j*dZ_k
   * for each (X,Y,Z) ∈ {(A,A,A), (A,A,B), (A,B,A), (A,B,B), (B,B,A), (B,B,B)}. */
  {
    struct conf { double *Xp, *Xq, *Yp, *Yq, *Zp, *Zq;
                  int *Xpl, *Xql, *Ypl, *Yql, *Zpl, *Zql;
                  int Xr, Yr, Zr;
                  double *Tp, *Tq; };
    struct conf cfgs[6] = {
      /* AAdA */ { ps->A_p, ps->A_q, ps->A_p, ps->A_q, ps->dA_p, ps->dA_q,
                   ps->A_p_len, ps->A_q_len, ps->A_p_len, ps->A_q_len,
                   ps->dA_p_len, ps->dA_q_len,
                   0, 0, -1, ps->T_AAdAp, ps->T_AAdAq },
      /* AAdB */ { ps->A_p, ps->A_q, ps->A_p, ps->A_q, ps->dB_p, ps->dB_q,
                   ps->A_p_len, ps->A_q_len, ps->A_p_len, ps->A_q_len,
                   ps->dB_p_len, ps->dB_q_len,
                   0, 0, -1, ps->T_AAdBp, ps->T_AAdBq },
      /* ABdA */ { ps->A_p, ps->A_q, ps->B_p, ps->B_q, ps->dA_p, ps->dA_q,
                   ps->A_p_len, ps->A_q_len, ps->B_p_len, ps->B_q_len,
                   ps->dA_p_len, ps->dA_q_len,
                   0, 0, -1, ps->T_ABdAp, ps->T_ABdAq },
      /* ABdB */ { ps->A_p, ps->A_q, ps->B_p, ps->B_q, ps->dB_p, ps->dB_q,
                   ps->A_p_len, ps->A_q_len, ps->B_p_len, ps->B_q_len,
                   ps->dB_p_len, ps->dB_q_len,
                   0, 0, -1, ps->T_ABdBp, ps->T_ABdBq },
      /* BBdA */ { ps->B_p, ps->B_q, ps->B_p, ps->B_q, ps->dA_p, ps->dA_q,
                   ps->B_p_len, ps->B_q_len, ps->B_p_len, ps->B_q_len,
                   ps->dA_p_len, ps->dA_q_len,
                   0, 0, -1, ps->T_BBdAp, ps->T_BBdAq },
      /* BBdB */ { ps->B_p, ps->B_q, ps->B_p, ps->B_q, ps->dB_p, ps->dB_q,
                   ps->B_p_len, ps->B_q_len, ps->B_p_len, ps->B_q_len,
                   ps->dB_p_len, ps->dB_q_len,
                   0, 0, -1, ps->T_BBdBp, ps->T_BBdBq },
    };
    int H = ps->H, p_cap = ps->p_cap, L_T = ps->L_T;
    int i, j, k, l, t;
    for (t = 0; t < 6; t++) {
      struct conf *c = &cfgs[t];
      for (i = 0; i < H; i++) {
        const double *Xi_p = &c->Xp[i * p_cap];
        const double *Xi_q = &c->Xq[i * p_cap];
        int Xi_pl = c->Xpl[i], Xi_ql = c->Xql[i];
        for (j = 0; j < H; j++) {
          const double *Yj_p = &c->Yp[j * p_cap];
          const double *Yj_q = &c->Yq[j * p_cap];
          int Yj_pl = c->Ypl[j], Yj_ql = c->Yql[j];
          int prod_r;
          /* prod = X_i * Y_j */
          ftp_pp_mul(ps->prod_p, p_cap, &ps->prod_p_len,
                      ps->prod_q, p_cap, &ps->prod_q_len,
                      &prod_r,
                      Xi_p, Xi_pl, Xi_q, Xi_ql, c->Xr,
                      Yj_p, Yj_pl, Yj_q, Yj_ql, c->Yr,
                      ps->scratch1, ps->sc_cap,
                      ps->scratch2, ps->sc_cap);
          for (k = 0; k < H; k++) {
            const double *Zk_p = &c->Zp[k * p_cap];
            const double *Zk_q = &c->Zq[k * p_cap];
            int Zk_pl = c->Zpl[k], Zk_ql = c->Zql[k];
            int triple_r;
            int offset = (i * H * H + j * H + k) * L_T;
            /* triple = prod * dZ_k */
            ftp_pp_mul(ps->triple_p, p_cap, &ps->triple_p_len,
                        ps->triple_q, p_cap, &ps->triple_q_len,
                        &triple_r,
                        ps->prod_p, ps->prod_p_len,
                        ps->prod_q, ps->prod_q_len, prod_r,
                        Zk_p, Zk_pl, Zk_q, Zk_ql, c->Zr,
                        ps->scratch1, ps->sc_cap,
                        ps->scratch2, ps->sc_cap);
            for (l = 0; l < ps->triple_p_len && l < L_T; l++)
              c->Tp[offset + l] = ps->triple_p[l];
            for (l = ps->triple_p_len; l < L_T; l++) c->Tp[offset + l] = 0.0;
            for (l = 0; l < ps->triple_q_len && l < L_T; l++)
              c->Tq[offset + l] = ps->triple_q[l];
            for (l = ps->triple_q_len; l < L_T; l++) c->Tq[offset + l] = 0.0;
          }
        }
      }
    }
  }
}

/* ---- Per-frequency polynomial path ----
 *
 * Given the bilinear data sums and the precomputed template PPs, build
 *   Pp(b) + sqrt(1-b^2) * Pq(b) = 0
 * by accumulating six PP triple-products weighted by K-tensor entries,
 * find the real roots of Pp^2 - (1-b^2)*Pq^2 in [-1, 1], and evaluate
 * P(omega, theta_2) at each candidate.  Returns the maximum P found. */
static double ftp_poly_P_at_omega(int H, _FTPPolyState *ps,
                                    int allow_neg_amp,
                                    const double *cn, const double *sn,
                                    const double *YC, const double *YS,
                                    const double *CC, const double *CS,
                                    const double *SS, double YY,
                                    double *scratch_A, double *scratch_B,
                                    double *theta_out, int *sign_out)
{
  int i, j, k, l;
  int p_cap = ps->p_cap;
  int L_cap = ps->L_cap;

  /* Zero accumulator polynomials. */
  for (l = 0; l < L_cap; l++) { ps->Pp[l] = 0.0; ps->Pq[l] = 0.0; }
  ps->Pp_len = 0; ps->Pq_len = 0;

  /* Build the six K tensors from this frequency's data sums.  Formulas
   * mirror pyftp::compute_zeros lines 432-437.
   *   K_aada[i,j,k] = YC[i]*CC[j,k] - YC[k]*CC[i,j]   (AAdA)
   *   K_aadb[i,j,k] = YC[i]*CS[j,k] - YS[k]*CC[i,j]   (AAdB)
   *   K_abda[i,j,k] = YC[i]*CS[k,j] + YS[j]*CC[i,k]   (ABdA)
   *   K_abdb[i,j,k] = YC[i]*SS[j,k] + YS[j]*CS[i,k]   (ABdB)
   *   K_bbda[i,j,k] = YS[i]*CS[k,j] - YC[k]*SS[i,j]   (BBdA)
   *   K_bbdb[i,j,k] = YS[i]*SS[j,k] - YS[k]*SS[i,j]   (BBdB)
   * (CC, CS, SS are row-major [n,m] = n*H + m for n,m in [0..H-1].)
   */
  #define IDX2(n, m) ((n) * H + (m))
  for (i = 0; i < H; i++) {
    for (j = 0; j < H; j++) {
      for (k = 0; k < H; k++) {
        int kidx = (i * H + j) * H + k;
        ps->K_aada[kidx] = YC[i] * CC[IDX2(j,k)] - YC[k] * CC[IDX2(i,j)];
        ps->K_aadb[kidx] = YC[i] * CS[IDX2(j,k)] - YS[k] * CC[IDX2(i,j)];
        ps->K_abda[kidx] = YC[i] * CS[IDX2(k,j)] + YS[j] * CC[IDX2(i,k)];
        ps->K_abdb[kidx] = YC[i] * SS[IDX2(j,k)] + YS[j] * CS[IDX2(i,k)];
        ps->K_bbda[kidx] = YS[i] * CS[IDX2(k,j)] - YC[k] * SS[IDX2(i,j)];
        ps->K_bbdb[kidx] = YS[i] * SS[IDX2(j,k)] - YS[k] * SS[IDX2(i,j)];
      }
    }
  }
  #undef IDX2

  /* Contract each of the 12 precomputed tensors with its K tensor and
   * accumulate into Pp / Pq.  Each contraction is
   *     Pp_or_Pq[l] += sum_{i,j,k} T[i,j,k,l] * K[i,j,k]
   * which costs O(H^3 * L_T).  Six p-contractions and six q-contractions
   * with H^3 * L_T flops each gives a per-frequency total of O(H^4) for
   * the polynomial coefficient construction.
   *
   * Attempting to fuse the 12 contractions into a single outer-ijk loop
   * with 12 fan-in streams turned out to be slightly slower (register
   * pressure prevents the compiler from vectorising the inner ll loop
   * as effectively as the original 1-in/1-out pattern), so we keep the
   * 12 separate passes. */
  {
    int L_T = ps->L_T;
    int H3 = H * H * H;
    int triplet;
    struct { double *T; double *K; double *target; } contracts[12] = {
      { ps->T_AAdAp, ps->K_aada, ps->Pp },
      { ps->T_AAdBp, ps->K_aadb, ps->Pp },
      { ps->T_ABdAp, ps->K_abda, ps->Pp },
      { ps->T_ABdBp, ps->K_abdb, ps->Pp },
      { ps->T_BBdAp, ps->K_bbda, ps->Pp },
      { ps->T_BBdBp, ps->K_bbdb, ps->Pp },
      { ps->T_AAdAq, ps->K_aada, ps->Pq },
      { ps->T_AAdBq, ps->K_aadb, ps->Pq },
      { ps->T_ABdAq, ps->K_abda, ps->Pq },
      { ps->T_ABdBq, ps->K_abdb, ps->Pq },
      { ps->T_BBdAq, ps->K_bbda, ps->Pq },
      { ps->T_BBdBq, ps->K_bbdb, ps->Pq },
    };
    /* Each contraction is target[l] += sum_{ijk} T[ijk, l] * K[ijk].
     * Tried cblas_dgemv (matrix-vector product); on systems linked
     * against the reference libgslcblas (Debian/Ubuntu default) it's
     * SLOWER than the hand-rolled gcc-vectorised loop because the
     * reference cblas isn't SIMD-tuned.  Stick with the hand-rolled
     * loop; users wanting more can swap in an optimised BLAS at link
     * time. */
    for (triplet = 0; triplet < 12; triplet++) {
      const double * __restrict__ T = contracts[triplet].T;
      const double * __restrict__ K = contracts[triplet].K;
      double       * __restrict__ target = contracts[triplet].target;
      int ijk;
      for (ijk = 0; ijk < H3; ijk++) {
        double Kijk = K[ijk];
        const double *Tijk = &T[ijk * L_T];
        int ll;
        if (Kijk == 0.0) continue;
        for (ll = 0; ll < L_T; ll++) target[ll] += Kijk * Tijk[ll];
      }
    }
    if (ps->Pp_len < L_T) ps->Pp_len = L_T;
    if (ps->Pq_len < L_T) ps->Pq_len = L_T;
    if (ps->Pp_len > L_cap) ps->Pp_len = L_cap;
    if (ps->Pq_len > L_cap) ps->Pq_len = L_cap;
  }

  /* Build root polynomial:  F(b) = Pp^2(b) - (1-b^2) * Pq^2(b).
   * (CORRECTED form -- pyftp uses (1-b^2)*Pp^2 - Pq^2, which is wrong.) */
  {
    int F_len;
    /* scratch1 = Pp^2 */
    ftp_poly_mul(ps->scratch1, ps->sc_cap,
                  ps->Pp, ps->Pp_len, ps->Pp, ps->Pp_len);
    int pp2_len = (ps->Pp_len > 0) ? (2 * ps->Pp_len - 1) : 0;
    if (pp2_len > ps->sc_cap) pp2_len = ps->sc_cap;
    /* scratch2 = Pq^2 */
    ftp_poly_mul(ps->scratch2, ps->sc_cap,
                  ps->Pq, ps->Pq_len, ps->Pq, ps->Pq_len);
    int qq2_len = (ps->Pq_len > 0) ? (2 * ps->Pq_len - 1) : 0;
    if (qq2_len > ps->sc_cap) qq2_len = ps->sc_cap;
    /* root_poly = (1-b^2) * scratch2 */
    ftp_poly_mul_oneMinusXsq(ps->root_poly, L_cap, ps->scratch2, qq2_len);
    int omxqq2_len = qq2_len + 2;
    if (omxqq2_len > L_cap) omxqq2_len = L_cap;
    /* root_poly = scratch1 - root_poly, i.e. Pp^2 - (1-b^2)*Pq^2 */
    for (l = 0; l < L_cap; l++) ps->root_poly[l] = -ps->root_poly[l];
    ftp_poly_addto(ps->root_poly, L_cap, ps->scratch1, pp2_len);
    F_len = pp2_len > omxqq2_len ? pp2_len : omxqq2_len;
    if (F_len > L_cap) F_len = L_cap;
    /* Trim trailing zeros. */
    while (F_len > 1 && fabs(ps->root_poly[F_len - 1]) < 1.0e-30) F_len--;
    ps->root_poly_len = F_len;
  }

  /* Find roots: try Newton refinement from the previous frequency's roots
   * if available; fall back to a full eigensolve at frequency 0, every
   * full_every frequencies, or when Newton fails. */
  {
    int nroots = ps->root_poly_len - 1;
    int r;
    double best_P = 0.0;
    double best_theta = 0.0;
    int best_sign = 1;
    int do_full = 0;
    int n_real_now = 0;
    double *real_roots_now = ps->scratch1;     /* reuse scratch1 buffer */

    if (nroots < 1) {
      if (theta_out) *theta_out = 0.0;
      if (sign_out)  *sign_out = 1;
      return 0.0;
    }
    if (fabs(ps->root_poly[nroots]) < 1.0e-30) {
      if (theta_out) *theta_out = 0.0;
      if (sign_out)  *sign_out = 1;
      return 0.0;
    }

    /* Build derivative of root_poly for Newton iterations. */
    {
      int Fp_len = ftp_poly_deriv ? 1 : 0;       /* (suppress warning) */
      (void) Fp_len;
      ftp_poly_deriv(ps->Fprime_coefs, L_cap, ps->root_poly, ps->root_poly_len);
    }

    /* Decide whether to do a full eigensolve. */
    if (!ps->prev_roots_valid || ps->freqs_since_full >= ps->full_every) {
      do_full = 1;
    } else {
      /* Newton refinement from previous frequency's roots. */
      int Fp_eff_len = (ps->root_poly_len > 0) ? (ps->root_poly_len - 1) : 0;
      int max_iters = 15;
      double tol = 1.0e-13;
      for (r = 0; r < ps->prev_n_real; r++) {
        double b = ps->prev_real_roots[r];
        int iter, converged = 0;
        for (iter = 0; iter < max_iters; iter++) {
          double F_val  = gsl_poly_eval(ps->root_poly,    ps->root_poly_len, b);
          double Fp_val = gsl_poly_eval(ps->Fprime_coefs, Fp_eff_len,         b);
          double delta;
          if (fabs(Fp_val) < 1.0e-200) { break; }
          delta = F_val / Fp_val;
          b -= delta;
          if (b < -1.5 || b > 1.5) { break; }    /* escaped feasible region */
          if (fabs(delta) < tol)   { converged = 1; break; }
        }
        if (!converged || b < -1.0 - 1.0e-6 || b > 1.0 + 1.0e-6) {
          /* This root went bad -- fall back to full eigensolve. */
          do_full = 1;
          break;
        }
        if (b > 1.0)  b = 1.0;
        if (b < -1.0) b = -1.0;
        real_roots_now[n_real_now++] = b;
      }
    }

    if (do_full) {
      int gsl_status;
      if (nroots + 1 != ps->gsl_ws_cap) {
        gsl_poly_complex_workspace_free(ps->gsl_ws);
        ps->gsl_ws = gsl_poly_complex_workspace_alloc((size_t) (nroots + 1));
        ps->gsl_ws_cap = nroots + 1;
        if (ps->gsl_ws == NULL) vt_error(ERR_MEMALLOC);
      }
      gsl_status = gsl_poly_complex_solve(ps->root_poly, (size_t) (nroots + 1),
                                           ps->gsl_ws, ps->roots_packed);
      if (gsl_status != GSL_SUCCESS) {
        if (theta_out) *theta_out = 0.0;
        if (sign_out)  *sign_out  = 1;
        return -1.0;
      }
      n_real_now = 0;
      for (r = 0; r < nroots; r++) {
        double br = ps->roots_packed[2 * r];
        double bi = ps->roots_packed[2 * r + 1];
        if (fabs(bi) > 1.0e-5) continue;
        if (br > 1.0)  br = 1.0;
        if (br < -1.0) br = -1.0;
        real_roots_now[n_real_now++] = br;
      }
      ps->freqs_since_full = 0;
    } else {
      ps->freqs_since_full++;
    }

    /* Evaluate P at each real root and each endpoint, both sgns; pick max. */
    {
      double bz[2] = { -1.0, 1.0 };
      int bi2, rr, s;
      for (rr = 0; rr < n_real_now; rr++) {
        double b = real_roots_now[rr];
        for (s = -1; s <= 1; s += 2) {
          double theta_2 = acos(b);
          int sign_eff;
          double P_cand;
          if (s < 0) theta_2 = 2.0 * M_PI - theta_2;
          P_cand = ftp_eval_P(H, theta_2, cn, sn, YC, YS, CC, CS, SS, YY,
                               scratch_A, scratch_B, &sign_eff);
          if (!allow_neg_amp && sign_eff < 0) continue;
          if (P_cand > best_P) {
            best_P = P_cand;
            best_theta = theta_2;
            best_sign = sign_eff;
          }
        }
      }
      for (bi2 = 0; bi2 < 2; bi2++) {
        for (s = -1; s <= 1; s += 2) {
          double theta_2 = acos(bz[bi2]);
          int sign_eff;
          double P_cand;
          if (s < 0) theta_2 = 2.0 * M_PI - theta_2;
          P_cand = ftp_eval_P(H, theta_2, cn, sn, YC, YS, CC, CS, SS, YY,
                               scratch_A, scratch_B, &sign_eff);
          if (!allow_neg_amp && sign_eff < 0) continue;
          if (P_cand > best_P) {
            best_P = P_cand;
            best_theta = theta_2;
            best_sign = sign_eff;
          }
        }
      }
    }

    /* Save real roots for next frequency's Newton refinement. */
    {
      int rr;
      for (rr = 0; rr < n_real_now; rr++)
        ps->prev_real_roots[rr] = real_roots_now[rr];
      ps->prev_n_real = n_real_now;
      ps->prev_roots_valid = 1;
    }

    if (theta_out) *theta_out = best_theta;
    if (sign_out)  *sign_out  = best_sign;
    return best_P;
  }
}

/* =====================================================================
 *  End of polynomial fast path
 *  ===================================================================== */


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

/* Per-frequency dispatcher: brute vs poly vs verify (both, compare). */
typedef struct {
  long n_verified;       /* number of frequencies where both methods ran */
  double max_abs_diff;   /* worst-case |P_poly - P_brute| seen */
  double sum_abs_diff;   /* for mean */
  long n_warned;         /* number of frequencies that exceeded FTP_VERIFY_TOL */
} _FTPVerifyStats;

static double ftp_one_freq_dispatch(int method, int H, _FTPPolyState *ps,
                                     int allow_neg_amp,
                                     const double *cn, const double *sn,
                                     const double *YC, const double *YS,
                                     const double *CC, const double *CS,
                                     const double *SS, double YY,
                                     double *scratch_A, double *scratch_B,
                                     double *Pgrid, int *signgrid,
                                     double *theta_out, int *sign_out,
                                     _FTPVerifyStats *vs)
{
  double P_brute = 0.0, P_poly = 0.0;
  double theta_b = 0.0, theta_p = 0.0;
  int sign_b = 1, sign_p = 1;
  if (method == FTP_METHOD_BRUTE || method == FTP_METHOD_VERIFY) {
    P_brute = ftp_max_P_over_theta(H, FTP_DEFAULT_NTHETA, allow_neg_amp,
                                     cn, sn, YC, YS, CC, CS, SS, YY,
                                     scratch_A, scratch_B, Pgrid, signgrid,
                                     &theta_b, &sign_b);
  }
  if ((method == FTP_METHOD_POLY || method == FTP_METHOD_VERIFY) && ps != NULL) {
    P_poly = ftp_poly_P_at_omega(H, ps, allow_neg_amp,
                                   cn, sn, YC, YS, CC, CS, SS, YY,
                                   scratch_A, scratch_B,
                                   &theta_p, &sign_p);
    if (P_poly < 0.0) {           /* poly sentinel: fall back to brute */
      P_poly = 0.0;
      if (method == FTP_METHOD_POLY) {
        /* No brute computed; trigger brute as fallback. */
        P_brute = ftp_max_P_over_theta(H, FTP_DEFAULT_NTHETA, allow_neg_amp,
                                         cn, sn, YC, YS, CC, CS, SS, YY,
                                         scratch_A, scratch_B, Pgrid, signgrid,
                                         &theta_b, &sign_b);
        if (theta_out) *theta_out = theta_b;
        if (sign_out)  *sign_out  = sign_b;
        return P_brute;
      }
      /* In verify mode, just record P_brute. */
    }
  }
  if (method == FTP_METHOD_VERIFY && vs != NULL && ps != NULL) {
    double d = fabs(P_poly - P_brute);
    vs->n_verified++;
    vs->sum_abs_diff += d;
    if (d > vs->max_abs_diff) vs->max_abs_diff = d;
    if (d > FTP_VERIFY_TOL) vs->n_warned++;
  }
  if (method == FTP_METHOD_POLY) {
    if (theta_out) *theta_out = theta_p;
    if (sign_out)  *sign_out  = sign_p;
    return P_poly;
  }
  /* BRUTE or VERIFY: return brute (it's the gold standard).  VERIFY's
   * purpose is to compare, not to switch behaviour. */
  if (theta_out) *theta_out = theta_b;
  if (sign_out)  *sign_out  = sign_b;
  return P_brute;
}

/* =====================================================================
 *  NFFT-batched sums (Phase B.5)
 *  ---------------------------------------------------------------------
 *  For a uniform frequency grid f_i = k_i * df (with k_i integers), the
 *  per-frequency sums C_n[f_i], S_n[f_i], YC_n[f_i], YS_n[f_i] for
 *  n=1..2H can be computed in a single batched O(M * log(NFFT_size) +
 *  NFFT_size) pass via the non-equispaced FFT (Keiner et al. 2009).
 *  Compared to the per-frequency direct loops (O(M * H * Nfreq)), this
 *  is a substantial speedup at large M (number of data points) -- a 4x
 *  win at M=10000 in our use case.
 *
 *  The C/S sums need NFFT bins up to k = 2H * k_max (one plan).  The
 *  YC/YS sums need bins up to k = H * k_max (second plan, smaller).
 *  Phase correction handles the time shift x_j = (t_j - t_mid) * df.
 *
 *  Returns 1 on success (sums filled), 0 if the frequency grid is not
 *  NFFT-compatible (non-uniform spacing, df * T_baseline >= 1, etc.).
 *  On failure the caller should fall back to direct per-frequency sums.
 *  ===================================================================== */
#ifdef HAVE_NFFT3
static int ftp_compute_sums_batch_nfft(int N, const double *t,
                                          const double *y_centered,
                                          const double *w, int H,
                                          const double *periods, int Nperiod,
                                          double *C_all, double *S_all,
                                          double *YC_all, double *YS_all)
{
  nfft_plan plan_w, plan_u;
  int N_nfft_w = 0, N_nfft_u = 0;
  int twoH = 2 * H;
  int k_max;
  double df, t_mid, tmin, tmax;
  double f0, f1, dfi;
  int i, j, n;
  int ret = 0;

  if (Nperiod < 2 || N < 2 || H < 1) return 0;

  /* Determine df from the first two periods.  Periods are stored in
   * decreasing-frequency order (highest f first), so f0 = 1/periods[0]
   * is the maximum frequency. */
  f0 = 1.0 / periods[0];
  f1 = 1.0 / periods[1];
  df = f0 - f1;
  if (!(df > 0.0)) return 0;

  /* Verify uniformity. */
  for (i = 0; i < Nperiod - 1; i++) {
    dfi = (1.0 / periods[i]) - (1.0 / periods[i + 1]);
    if (fabs(dfi - df) > 1.0e-9 * df) return 0;
  }

  /* k_max = highest NFFT bin we'll touch = round(f0 / df). */
  k_max = (int) (f0 / df + 0.5);
  if (k_max < 1) return 0;

  /* Centre the time series at t_mid for the NFFT node range. */
  tmin = t[0]; tmax = t[0];
  for (j = 1; j < N; j++) {
    if (t[j] < tmin) tmin = t[j];
    if (t[j] > tmax) tmax = t[j];
  }
  t_mid = 0.5 * (tmin + tmax);

  /* NFFT node range check: max |x_j| = (T/2) * df < 1/2  iff  df * T < 1.
   * The frequency grid is df = subsample/T_baseline (where subsample < 1
   * for typical use), so this is comfortably satisfied. */
  if ((tmax - tmin) * df >= 1.0) return 0;

  /* Choose NFFT plan sizes: need bins up to ±(2H * k_max) for plan_w and
   * ±(H * k_max) for plan_u.  Round up to next power of 2 for efficient
   * FFTW behaviour. */
  {
    int need_w = 2 * (twoH * k_max + 1);
    int need_u = 2 * (H * k_max + 1);
    N_nfft_w = 2;
    while (N_nfft_w < need_w) N_nfft_w *= 2;
    N_nfft_u = 2;
    while (N_nfft_u < need_u) N_nfft_u *= 2;
  }

  /* Init plans. */
  nfft_init_1d(&plan_w, N_nfft_w, N);
  nfft_init_1d(&plan_u, N_nfft_u, N);

  /* Fill nodes and function values.  plan.f is fftw_complex (double[2]):
   * index [0] is the real part, [1] is the imag part. */
  for (j = 0; j < N; j++) {
    double xj = (t[j] - t_mid) * df;
    plan_w.x[j] = xj;
    plan_u.x[j] = xj;
    plan_w.f[j][0] = w[j];
    plan_w.f[j][1] = 0.0;
    plan_u.f[j][0] = w[j] * y_centered[j];
    plan_u.f[j][1] = 0.0;
  }

  /* Precompute interpolation weights and run adjoint. */
  nfft_precompute_one_psi(&plan_w);
  nfft_precompute_one_psi(&plan_u);
  nfft_adjoint(&plan_w);
  nfft_adjoint(&plan_u);

  /* Extract per-(i, n) sums with phase correction.
   * desired[k] = NFFT_result[k] * exp(+2*pi*i*k*df*t_mid). */
  for (i = 0; i < Nperiod; i++) {
    int k_i = (int) ((1.0 / periods[i]) / df + 0.5);
    for (n = 1; n <= twoH; n++) {
      int k_eff = n * k_i;
      int idx_w = k_eff + N_nfft_w / 2;
      double phase = 2.0 * M_PI * (double) k_eff * df * t_mid;
      double cp = cos(phase), sp = sin(phase);
      double re = plan_w.f_hat[idx_w][0];
      double im = plan_w.f_hat[idx_w][1];
      double re_c = re * cp - im * sp;
      double im_c = re * sp + im * cp;
      C_all[i * twoH + (n - 1)] = re_c;
      S_all[i * twoH + (n - 1)] = im_c;
      if (n <= H) {
        int idx_u = k_eff + N_nfft_u / 2;
        re = plan_u.f_hat[idx_u][0];
        im = plan_u.f_hat[idx_u][1];
        re_c = re * cp - im * sp;
        im_c = re * sp + im * cp;
        YC_all[i * H + (n - 1)] = re_c;
        YS_all[i * H + (n - 1)] = im_c;
      }
    }
  }

  nfft_finalize(&plan_w);
  nfft_finalize(&plan_u);
  ret = 1;
  return ret;
}
#endif /* HAVE_NFFT3 */

#define FTP_SUMS_DIRECT 0
#define FTP_SUMS_NFFT   1
#define FTP_SUMS_AUTO   2   /* default: NFFT if compiled, else direct */

/* Compute periodogram across periods[0..Nperiod-1].  When sums_mode == NFFT
 * (or AUTO with HAVE_NFFT3), pre-compute all C/S/YC/YS sums via NFFT in a
 * single batched pass and then dispatch per-frequency using the cached
 * sums; otherwise fall back to the direct per-frequency summation path. */
static void ftp_periodogram(int N, double *t, double *y_centered, double *w,
                             double YY, int H, double *cn, double *sn,
                             int allow_neg_amp, _FTPScratch *sc,
                             int method, _FTPPolyState *ps, _FTPVerifyStats *vs,
                             int sums_mode,
                             int Nperiod, double *periods, double *periodogram,
                             int *negamp_grid, double *theta_grid)
{
  int i;
  int twoH = 2 * H;
  double *C_all = NULL, *S_all = NULL, *YC_all = NULL, *YS_all = NULL;
  int use_nfft = 0;
  (void) sums_mode;

#ifdef HAVE_NFFT3
  if (sums_mode == FTP_SUMS_NFFT || sums_mode == FTP_SUMS_AUTO) {
    C_all  = (double *) malloc(Nperiod * twoH * sizeof(double));
    S_all  = (double *) malloc(Nperiod * twoH * sizeof(double));
    YC_all = (double *) malloc(Nperiod * H    * sizeof(double));
    YS_all = (double *) malloc(Nperiod * H    * sizeof(double));
    if (C_all == NULL || S_all == NULL || YC_all == NULL || YS_all == NULL)
      vt_error(ERR_MEMALLOC);
    use_nfft = ftp_compute_sums_batch_nfft(N, t, y_centered, w, H,
                                            periods, Nperiod,
                                            C_all, S_all, YC_all, YS_all);
    if (!use_nfft) {
      /* NFFT couldn't be used (non-uniform grid, etc.).  Fall back. */
      free(C_all); free(S_all); free(YC_all); free(YS_all);
      C_all = S_all = YC_all = YS_all = NULL;
    }
  }
#else
  if (sums_mode == FTP_SUMS_NFFT) {
    fprintf(stderr, "-FTP method ... 'sums nfft': vartools was built without "
                    "NFFT support (configure with --with-nfft); falling back to direct sums\n");
  }
#endif

  for (i = 0; i < Nperiod; i++) {
    double omega = 2.0 * M_PI / periods[i];
    double th = 0.0;
    int sgn = 1;
    if (use_nfft) {
      /* Copy cached sums into the per-frequency scratch buffers. */
      memcpy(sc->C,  &C_all[i * twoH], twoH * sizeof(double));
      memcpy(sc->S,  &S_all[i * twoH], twoH * sizeof(double));
      memcpy(sc->YC, &YC_all[i * H],   H    * sizeof(double));
      memcpy(sc->YS, &YS_all[i * H],   H    * sizeof(double));
    } else {
      ftp_compute_sums(N, t, y_centered, w, omega, H,
                        sc->C, sc->S, sc->YC, sc->YS);
    }
    ftp_compute_bilinears(H, sc->C, sc->S, sc->CC, sc->CS, sc->SS);
    periodogram[i] = ftp_one_freq_dispatch(method, H, ps, allow_neg_amp,
                                             cn, sn, sc->YC, sc->YS,
                                             sc->CC, sc->CS, sc->SS, YY,
                                             sc->A, sc->B, sc->Pgrid, sc->signgrid,
                                             &th, &sgn, vs);
    if (negamp_grid != NULL) negamp_grid[i] = (sgn < 0) ? 1 : 0;
    if (theta_grid  != NULL) theta_grid[i]  = th;
  }

  if (C_all) { free(C_all); free(S_all); free(YC_all); free(YS_all); }
}

/* One-shot P(omega) at a single test period (used in fine-tune + multiple checks). */
static double ftp_one_period(int N, double *t, double *y_centered, double *w,
                              double YY, int H, double *cn, double *sn,
                              int allow_neg_amp, _FTPScratch *sc,
                              int method, _FTPPolyState *ps,
                              double period,
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
  P = ftp_one_freq_dispatch(method, H, ps, allow_neg_amp,
                              cn, sn, sc->YC, sc->YS, sc->CC, sc->CS, sc->SS, YY,
                              sc->A, sc->B, sc->Pgrid, sc->signgrid,
                              theta_out, sign_out, NULL);
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
                    int method, int sums_mode,
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
  _FTPPolyState ps = {0};
  _FTPVerifyStats vs = {0};
  _FTPPolyState *ps_ptr = NULL;
  _FTPVerifyStats *vs_ptr = (method == FTP_METHOD_VERIFY) ? &vs : NULL;
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

  /* Allocate polynomial state if the chosen method needs it.  Built once
   * per LC; reused across the N_f-frequency sweep, the fine-tune sweep,
   * and the period-multiple double-check. */
  if (method == FTP_METHOD_POLY || method == FTP_METHOD_VERIFY) {
    ftp_poly_state_alloc(&ps, H);
    ftp_poly_build_templates(&ps, cn, sn);
    ps_ptr = &ps;
  }

  /* Compute periodogram. */
  ftp_periodogram(Nused, t, y_centered, w, YY, H, cn, sn,
                  allow_neg_amp, &sc, method, ps_ptr, vs_ptr,
                  sums_mode,
                  Nperiod, periods, periodogram,
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
                                allow_neg_amp, &sc, method, ps_ptr, testperiod, &th, &sgn);
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
                                      allow_neg_amp, &sc, method, ps_ptr, testperiod, &th, &sgn);
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
                                allow_neg_amp, &sc, method, ps_ptr, perpeaks[j], &th, &sgn);
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

  /* Verify-mode summary: emit stats so the user knows the poly path is
   * agreeing with brute-force to within tolerance.  Per-LC summary line
   * so multi-LC runs don't drown the output. */
  if (method == FTP_METHOD_VERIFY && vs.n_verified > 0) {
    fprintf(stderr,
            "-FTP method verify: LC %d: %ld freqs compared, max|P_poly - P_brute| = %.3e, "
            "mean = %.3e, %ld over tol %.1e\n",
            lc_name_num, vs.n_verified, vs.max_abs_diff,
            vs.sum_abs_diff / (double) vs.n_verified, vs.n_warned, FTP_VERIFY_TOL);
  }

  free(t); free(mag); free(w); free(y_centered);
  free(periods); free(periodogram);
  free(negamp_grid); free(theta_grid);
  ftp_scratch_free(&sc);
  if (ps_ptr != NULL) ftp_poly_state_free(&ps);
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

  /* Filelist mode: each LC has its own template file (path stored in
   * Ftp->template_filenames[lc_name_num]).  Load it now, replacing the
   * previous LC's cn/sn.  H can vary across LCs; the per-LC scratch in
   * findPeaks_ftp is allocated fresh per call so it follows automatically. */
  if (Ftp->filelist_mode) {
    int idx;
    if (Ftp->cn != NULL) { free(Ftp->cn); Ftp->cn = NULL; }
    if (Ftp->sn != NULL) { free(Ftp->sn); Ftp->sn = NULL; }
    Ftp->H = ftp_load_template_file(Ftp->template_filenames[lc_name_num],
                                     &Ftp->cn, &Ftp->sn);
    if (Ftp->H <= 0) {
      fprintf(stderr, "-FTP filelist: skipping LC '%s' due to failed template "
                      "load from '%s'\n",
              p->lcnames[lc_name_num], Ftp->template_filenames[lc_name_num]);
      for (idx = 0; idx < Ftp->Npeaks; idx++) {
        Ftp->peakperiods[lcnum][idx] = -1.0;
        Ftp->peakvalues[lcnum][idx]  = 0.0;
        Ftp->peakSNR[lcnum][idx]     = 0.0;
        Ftp->peakNegAmp[lcnum][idx]  = 0;
        Ftp->peakTheta[lcnum][idx]   = 0.0;
      }
      Ftp->avepower[lcnum] = 0.0;
      Ftp->rmspower[lcnum] = 1.0;
      return;
    }
  }

  /* Inline mode: resolve each c_k / s_k from its per-slot source (literal,
   * variable, or expression) into the shared Ftp->cn / Ftp->sn buffers.
   * For file / fitlc modes this loop is skipped (inline_mode == 0) and the
   * buffers already contain the parser-time-loaded template. */
  if (Ftp->inline_mode) {
    int kk;
    for (kk = 0; kk < Ftp->H; kk++) {
      switch (Ftp->cn_source[kk]) {
        case VARTOOLS_SOURCE_EVALEXPRESSION:
          Ftp->cn[kk] = EvaluateExpression(lc_name_num, lcnum, 0, Ftp->cn_expr[kk]);
          break;
        case VARTOOLS_SOURCE_EXISTINGVARIABLE:
          Ftp->cn[kk] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Ftp->cn_var[kk]);
          break;
        default:
          Ftp->cn[kk] = Ftp->cn_lit[kk];
          break;
      }
      switch (Ftp->sn_source[kk]) {
        case VARTOOLS_SOURCE_EVALEXPRESSION:
          Ftp->sn[kk] = EvaluateExpression(lc_name_num, lcnum, 0, Ftp->sn_expr[kk]);
          break;
        case VARTOOLS_SOURCE_EXISTINGVARIABLE:
          Ftp->sn[kk] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Ftp->sn_var[kk]);
          break;
        default:
          Ftp->sn[kk] = Ftp->sn_lit[kk];
          break;
      }
    }
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
                Ftp->method, Ftp->sums_mode,
                lcnum, lc_name_num);
}
