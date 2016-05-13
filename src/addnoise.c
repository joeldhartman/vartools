/*     This file is part of VARTOOLS version 1.31                      */
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
/*     Copyright 2007, 2008, 2009  Joel Hartman                              */
/*                                                                           */
#include "commands.h"
#include "programdata.h"
#include "functions.h"

#ifdef _HAVE_GSL
#include <gsl/gsl_sort.h>
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_sf.h>

#define MAXSIZE_NOISESIM 20
#endif

#ifndef VARTOOLS_LOG_COVARIANCE_TINY
#define VARTOOLS_LOG_COVARIANCE_TINY -23
#endif

/* This file contains the function to add time-correlated noise to a 
   light curve */

/* Follows the procedure to generate 1/f noise via wavelets described
   by McCoy & Walden, 1996, Journal of Computational and Graphical
   Statistics, Vol. 5, No. 1, p. 26-56 (specifically section 4.2 in
   that paper) */

double ran1(void)
{
  return (double) (rand() / (double) RAND_MAX);
}

double gasdev(void)
{
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;

  if(iset == 0) {
    do {
      v1 = 2.0*ran1() - 1.0;
      v2 = 2.0*ran1() - 1.0;
      rsq = v1*v1 + v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }
}

#ifdef _HAVE_GSL

void addnoise_wavelet(ProgramData *p, double gamval, double sig_r, double sig_w, int lc)
{
  int i, k, Nin, N, n0, M, m, j, j1, j2;

  double *datain, *data, minsep, sep, r1, r2, *dataout;

  double d, sig2tot, sig2sum, sigbm, sigsm;

  size_t *p_;

  gsl_wavelet *w;
  gsl_wavelet_workspace *work;

  double val;

  if((datain = (double *) malloc(p->NJD[lc] * sizeof(double))) == NULL ||
     (dataout = (double *) malloc(p->NJD[lc] * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  Nin = p->NJD[lc];

  for(i=0; i < Nin; i++) {
    datain[i] = p->t[lc][i];
  }

  if(gamval >= 1.)
    gamval = 0.999;
  if(gamval <= -1.)
    gamval = -0.999;

  d = gamval / 2.;

  /* Standard Deviation of full time series is twice sig_r when sig_w
     is set to 0. Divide sig_r by 2 to correct for this. */
  sig_r = sig_r / 2.;

  n0 = 2;

  if((p_ = malloc (Nin * sizeof (size_t))) == NULL)
    error(ERR_MEMALLOC);
  
  gsl_sort_index (p_, datain, 1, Nin);
  
  if(Nin <= 4) {
    fprintf(stderr,"Input data stream to short.\n");
    exit(1);
  }
  minsep = datain[p_[1]] - datain[p_[0]];
  for(i=2; i < Nin; i++) {
    if(datain[p_[i]] - datain[p_[i-1]] < minsep) {
      minsep = datain[p_[i]] - datain[p_[i-1]];
    }
  }
  if(minsep > 0) {
    N = ceil(log((double) Nin)/log(2.));
    if(N > MAXSIZE_NOISESIM) {
      N = MAXSIZE_NOISESIM;
    }
    M = N;
    N = 2 << N;
  } else {
    M = MAXSIZE_NOISESIM;
    N = 2 << MAXSIZE_NOISESIM;
  }
  
  w = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4);
  work = gsl_wavelet_workspace_alloc (N);
  
  if((data = (double *) malloc(N * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  //g = 1./(2.*log(2.));

  j = 1;
  sig2sum = 0.;
  for(m=M; m >= 1; m--) {
    for(k=1; k <= 2<<(M-m); k++) {
      if(k == 1) {
	sigbm = 2.*pow((2*M_PI),-2.*d)*sig_r*sig_r*pow(2.,(-m*(1-2.*d)))*(1 - pow(2.,(2.*d - 1.)))/(1. - 2.*d);
	sig2sum += sigbm;
	sigsm = (2 << m)*sigbm;
	sigsm = sqrt(sigsm);
      }
      data[j] = sigsm*gasdev();
      j++;
    }
  }

  sig2tot = sig_r*sig_r*gsl_sf_gamma((1. - 2.*d))/gsl_sf_gamma((1. - d))/gsl_sf_gamma((1.-d));
  sigbm = sig2tot - sig2sum;
  sigsm = (2 << M)*sigbm;
  sigsm = sqrt(sigsm);
  data[0] = sigsm*gasdev();
  
  gsl_wavelet_transform_inverse (w, data, 1, N, work);

  sep = (datain[p_[Nin-1]] - datain[p_[0]])/ ((double) N);
  
  for (i=0; i < Nin; i++) {
    j1 = floor((datain[p_[i]] - datain[p_[0]])/sep);
    j2 = j1+1;
    r1 = datain[p_[0]] + j1*sep;
    r2 = datain[p_[0]] + j2*sep;
    dataout[i] = data[j1] + (datain[p_[i]] - r1)*(data[j2]-data[j1])/sep;
    dataout[i] = dataout[i] + sig_w*gasdev();
  }

  val = 0.;
  for(i=0; i < Nin; i++) {
    val += dataout[i];
  }

  val = val / Nin;

  for(i=0; i < Nin; i++) {

    p->mag[lc][p_[i]] += (dataout[i] - val);
  }

  
  gsl_wavelet_free (w);
  gsl_wavelet_workspace_free (work);
  
  free (data);
  free (datain);
  free (dataout);
  free (p_);
  return;
}

#endif

void addnoise_covar_squareexp(ProgramData *p, double cov_rho, double sig_r, double sig_w, int lc, double bintime)
{
  int i, j, i0, k;
  double **Cov;
  int *Nvec;
  int NJD;
  double *t, *mag, *err, term, *b, *p_chol;

  NJD = p->NJD[lc];
  t = p->t[lc];
  mag = p->mag[lc];
  err = p->sig[lc];

  if(NJD <= 0) return;
  if((Cov = (double **) malloc(NJD * sizeof(double *))) == NULL ||
     (Nvec = (int *) malloc(NJD * sizeof(int))) == NULL ||
     (b = (double *) malloc(NJD * sizeof(double))) == NULL ||
     (p_chol = (double *) malloc(NJD * sizeof(double))) == NULL) {
    error(ERR_MEMALLOC);
  }
  for(i=0; i < NJD; i++) {
    if((Cov[i] = (double *) malloc(NJD * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
  }
  if(cov_rho <= 0.0) {
    error(ERR_NEGATIVE_COVAR_PARAM_ADDNOISE);
  }
  if(sig_r <= 0.0) {
    error(ERR_NEGATIVE_COVAR_PARAM_ADDNOISE);
  }
  cov_rho = cov_rho*cov_rho;
  sig_r = sig_r*sig_r;
  sig_w = sig_w*sig_w;
  if(bintime <= 0) {
    for(i=0; i < NJD; i++) {
      Cov[i][i] = sig_w + sig_r;
      b[i] = gasdev();
      for(j = i+1; j < NJD; j++) {
	term = -(t[i]-t[j])*(t[i]-t[j])/2.0/cov_rho;
	if(term < VARTOOLS_LOG_COVARIANCE_TINY) {
	  Cov[i][j] = 0.0;
	  break;
	}
	Cov[i][j] = sig_r*exp(term);
      }
    }
    choldc_sparse_neardiag(Cov, NJD, Nvec, p_chol);
    cholmult_sparse_neardiag(Cov, NJD, Nvec, p_chol, b, b);
    for(i=0; i < NJD; i++) {
      mag[i] = mag[i] + b[i];
    }
  } else {
    i = 0;
    do {
      i0 = i;
      for(; i < NJD ? (t[i] - t[i0] < bintime) : 0; i++) {
	k = i - i0;
	Cov[k][k] = sig_w + sig_r;
	b[k] = gasdev();
	for(j = i+1; j < NJD ? (t[j] - t[i0] < bintime) : 0; j++) {
	  term = -(t[i]-t[j])*(t[i]-t[j])/2.0/cov_rho;
	  if(term < VARTOOLS_LOG_COVARIANCE_TINY) {
	    Cov[k][j-i0] = 0.0;
	    break;
	  }
	  Cov[k][j-i0] = sig_r*exp(term);
	}
      }
      choldc_sparse_neardiag(Cov, i-i0, Nvec, p_chol);
      cholmult_sparse_neardiag(Cov, i-i0, Nvec, p_chol, b, b);
      for(k=i0; k < i; k++) {
	mag[k] = mag[k] + b[k-i0];
      }
    } while(i < NJD);
  }
  free(Nvec);
  free(p_chol);
  free(b);
  for(i=0; i < NJD; i++)
    free(Cov[i]);
  free(Cov);
}

void addnoise_covar_exp(ProgramData *p, double cov_rho, double sig_r, double sig_w, int lc, double bintime)
{
  int i, j, i0, k;
  double **Cov;
  int *Nvec;
  int NJD;
  double *t, *mag, *err, term, *b, *p_chol;

  NJD = p->NJD[lc];
  t = p->t[lc];
  mag = p->mag[lc];
  err = p->sig[lc];

  if(NJD <= 0) return;
  if((Cov = (double **) malloc(NJD * sizeof(double *))) == NULL ||
     (Nvec = (int *) malloc(NJD * sizeof(int))) == NULL ||
     (b = (double *) malloc(NJD * sizeof(double))) == NULL ||
     (p_chol = (double *) malloc(NJD * sizeof(double))) == NULL) {
    error(ERR_MEMALLOC);
  }
  for(i=0; i < NJD; i++) {
    if((Cov[i] = (double *) malloc(NJD * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
  }
  if(cov_rho <= 0.0) {
    error(ERR_NEGATIVE_COVAR_PARAM_ADDNOISE);
  }
  if(sig_r <= 0.0) {
    error(ERR_NEGATIVE_COVAR_PARAM_ADDNOISE);
  }
  sig_r = sig_r*sig_r;
  sig_w = sig_w*sig_w;
  if(bintime <= 0) {
    for(i=0; i < NJD; i++) {
      Cov[i][i] = sig_w + sig_r;
      b[i] = gasdev();
      for(j = i+1; j < NJD; j++) {
	term = -fabs((t[i]-t[j]))/cov_rho;
	if(term < VARTOOLS_LOG_COVARIANCE_TINY) {
	  Cov[i][j] = 0.0;
	  break;
	}
	Cov[i][j] = sig_r*exp(term);
      }
    }
    choldc_sparse_neardiag(Cov, NJD, Nvec, p_chol);
    cholmult_sparse_neardiag(Cov, NJD, Nvec, p_chol, b, b);
    for(i=0; i < NJD; i++) {
      mag[i] = mag[i] + b[i];
    }
  } else {
    i = 0;
    do {
      i0 = i;
      for(; i < NJD ? (t[i] - t[i0] < bintime) : 0; i++) {
	k = i - i0;
	Cov[k][k] = sig_w + sig_r;
	b[k] = gasdev();
	for(j = i+1; j < NJD ? (t[j] - t[i0] < bintime) : 0; j++) {
	  term = -fabs((t[i]-t[j]))/cov_rho;
	  if(term < VARTOOLS_LOG_COVARIANCE_TINY) {
	    Cov[k][j-i0] = 0.0;
	    break;
	  }
	  Cov[k][j-i0] = sig_r*exp(term);
	}
      }
      choldc_sparse_neardiag(Cov, i-i0, Nvec, p_chol);
      cholmult_sparse_neardiag(Cov, i-i0, Nvec, p_chol, b, b);
      for(k=i0; k < i; k++) {
	mag[k] = mag[k] + b[k-i0];
      }
    } while(i < NJD);
  }
  free(Nvec);
  free(p_chol);
  free(b);
  for(i=0; i < NJD; i++)
    free(Cov[i]);
  free(Cov);
}


void addnoise_covar_matern(ProgramData *p, double cov_nu, double cov_rho, double sig_r, double sig_w, int lc, double bintime)
{
  int i, j, i0, k;
  double **Cov;
  int *Nvec;
  int NJD;
  double *t, *mag, *err, term, *b, *p_chol;
  double x1, constterm, ri, rk, rip, rkp;

  NJD = p->NJD[lc];
  t = p->t[lc];
  mag = p->mag[lc];
  err = p->sig[lc];

  if(NJD <= 0) return;
  if((Cov = (double **) malloc(NJD * sizeof(double *))) == NULL ||
     (Nvec = (int *) malloc(NJD * sizeof(int))) == NULL ||
     (b = (double *) malloc(NJD * sizeof(double))) == NULL ||
     (p_chol = (double *) malloc(NJD * sizeof(double))) == NULL) {
    error(ERR_MEMALLOC);
  }
  for(i=0; i < NJD; i++) {
    if((Cov[i] = (double *) malloc(NJD * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
  }
  if(cov_nu <= 0.0) {
    error(ERR_NEGATIVE_COVAR_PARAM_ADDNOISE);
  }
  if(cov_rho <= 0.0) {
    error(ERR_NEGATIVE_COVAR_PARAM_ADDNOISE);
  }
  if(sig_r <= 0.0) {
    error(ERR_NEGATIVE_COVAR_PARAM_ADDNOISE);
  }
  sig_r = sig_r*sig_r;
  sig_w = sig_w*sig_w;
  x1 = sqrt(2.0*cov_nu)/cov_rho;
  constterm = sig_r/(exp(gammln(cov_nu))*pow(2.0,(cov_nu-1.0)));
  if(bintime <= 0) {
    for(i=0; i < NJD; i++) {
      Cov[i][i] = sig_w + sig_r;
      b[i] = gasdev();
      for(j = i+1; j < NJD; j++) {
	term = x1*fabs(t[i]-t[j]);
	if(-term < VARTOOLS_LOG_COVARIANCE_TINY) {
	  Cov[i][j] = 0.0;
	  break;
	}
	bessik(term, cov_nu, &ri, &rk, &rip, &rkp);
	Cov[i][j] = constterm*(pow(term,cov_nu))*rk;
      }
    }
    choldc_sparse_neardiag(Cov, NJD, Nvec, p_chol);
    cholmult_sparse_neardiag(Cov, NJD, Nvec, p_chol, b, b);
    for(i=0; i < NJD; i++) {
      mag[i] = mag[i] + b[i];
    }
  }
  else {
    i = 0;
    do {
      i0 = i;
      for(; i < NJD ? (t[i] - t[i0] < bintime) : 0; i++) {
	k = i - i0;
	Cov[k][k] = sig_w + sig_r;
	b[k] = gasdev();
	for(j = i+1; j < NJD ? (t[j] - t[i0] < bintime) : 0; j++) {
	  term = x1*fabs(t[i]-t[j]);
	  if(term < VARTOOLS_LOG_COVARIANCE_TINY) {
	    Cov[k][j-i0] = 0.0;
	    break;
	  }
	  bessik(term, cov_nu, &ri, &rk, &rip, &rkp);
	  Cov[k][j-i0] = constterm*(pow(term,cov_nu))*rk;
	}
      }
      choldc_sparse_neardiag(Cov, i-i0, Nvec, p_chol);
      cholmult_sparse_neardiag(Cov, i-i0, Nvec, p_chol, b, b);
      for(k=i0; k < i; k++) {
	mag[k] = mag[k] + b[k-i0];
      }
    } while(i < NJD);
  }
  free(Nvec);
  free(p_chol);
  free(b);
  for(i=0; i < NJD; i++)
    free(Cov[i]);
  free(Cov);
}

void addnoise_white(ProgramData *p, double sig_w, int lc)
{
  int i, j;
  int NJD;
  double *t, *mag, *err;

  NJD = p->NJD[lc];
  t = p->t[lc];
  mag = p->mag[lc];
  err = p->sig[lc];

  if(NJD <= 0) return;

  if(sig_w < 0.0) {
    error(ERR_NEGATIVE_COVAR_PARAM_ADDNOISE);
  }

  for(i=0; i < NJD; i++) {
    mag[i] = mag[i] + sig_w*gasdev();
  }
}

void addnoise(ProgramData *p, _AddNoise *c, int threadid, int lcid) {
  double d0, d1, d2, d3, d4;
  if(c->noise_type == VARTOOLS_ADDNOISE_COVAR_SQUAREDEXPONENTIAL) {
    sortlcbytime(p->NJD[threadid], p->t[threadid], threadid, p);
    if(c->rho_r_type == PERTYPE_SPECIFIED)
      d1 = c->rho_r[lcid][0];
    else
      d1 = c->rho_r_fix;
    if(c->sig_r_type == PERTYPE_SPECIFIED)
      d2 = c->sig_r[lcid][0];
    else
      d2 = c->sig_r_fix;
    if(c->sig_w_type == PERTYPE_SPECIFIED)
      d3 = c->sig_w[lcid][0];
    else
      d3 = c->sig_w_fix;
    if(c->bintime_type == PERTYPE_SPECIFIED)
      d4 = c->bintime[lcid][0];
    else
      d4 = c->bintime_fix;
    addnoise_covar_squareexp(p, d1, d2, d3, threadid, d4);
  }
  if(c->noise_type == VARTOOLS_ADDNOISE_COVAR_EXPONENTIAL) {
    sortlcbytime(p->NJD[threadid], p->t[threadid], threadid, p);
    if(c->rho_r_type == PERTYPE_SPECIFIED)
      d1 = c->rho_r[lcid][0];
    else
      d1 = c->rho_r_fix;
    if(c->sig_r_type == PERTYPE_SPECIFIED)
      d2 = c->sig_r[lcid][0];
    else
      d2 = c->sig_r_fix;
    if(c->sig_w_type == PERTYPE_SPECIFIED)
      d3 = c->sig_w[lcid][0];
    else
      d3 = c->sig_w_fix;
    if(c->bintime_type == PERTYPE_SPECIFIED)
      d4 = c->bintime[lcid][0];
    else
      d4 = c->bintime_fix;
    addnoise_covar_exp(p, d1, d2, d3, threadid, d4);
  }
  else if(c->noise_type == VARTOOLS_ADDNOISE_COVAR_MATERN) {
    sortlcbytime(p->NJD[threadid], p->t[threadid], threadid, p);
    if(c->nu_r_type == PERTYPE_SPECIFIED)
      d0 = c->nu_r[lcid][0];
    else
      d0 = c->nu_r_fix;
    if(c->rho_r_type == PERTYPE_SPECIFIED)
      d1 = c->rho_r[lcid][0];
    else
      d1 = c->rho_r_fix;
    if(c->sig_r_type == PERTYPE_SPECIFIED)
      d2 = c->sig_r[lcid][0];
    else
      d2 = c->sig_r_fix;
    if(c->sig_w_type == PERTYPE_SPECIFIED)
      d3 = c->sig_w[lcid][0];
    else
      d3 = c->sig_w_fix;
    if(c->bintime_type == PERTYPE_SPECIFIED)
      d4 = c->bintime[lcid][0];
    else
      d4 = c->bintime_fix;
    addnoise_covar_matern(p, d0, d1, d2, d3, threadid, d4);
  }
  else if(c->noise_type == VARTOOLS_ADDNOISE_WHITE) {
    if(c->sig_w_type == PERTYPE_SPECIFIED)
      d3 = c->sig_w[lcid][0];
    else
      d3 = c->sig_w_fix;
    addnoise_white(p, d3, threadid);
  }    
#ifdef _HAVE_GSL
  else if(c->noise_type == VARTOOLS_ADDNOISE_WAVELET) {
    if(sortlcbytime(p->NJD[threadid], p->t[threadid], threadid, p))
      mergeequallctimes(p, threadid);
    if(c->gammaval_type == PERTYPE_SPECIFIED)
      d1 = c->gammaval[lcid][0];
    else
      d1 = c->gammaval_fix;
    if(c->sig_r_type == PERTYPE_SPECIFIED)
      d2 = c->sig_r[lcid][0];
    else
      d2 = c->sig_r_fix;
    if(c->sig_w_type == PERTYPE_SPECIFIED)
      d3 = c->sig_w[lcid][0];
    else
      d3 = c->sig_w_fix;
    addnoise_wavelet(p, d1, d2, d3, threadid);
  }
#endif
}

