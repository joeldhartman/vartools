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
/*     This file is part of VARTOOLS version 1.152                      */
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

/* This file contains functions for computing the Lomb-Scargle periodogram for the program vartools by J. Hartman */

/* This function calculates a fast Lomb-Scargle periodogram using the Press & Rybicki method (1989, ApJ, 338, 277) as implemented by Press et al. 1992 */

#define MOD(a,b) while((a) >= (b)) (a) -= (b)
#define MACC 4

void LS_getfixedperiodSNR(double *f, double *p, int N, double ave, double rms, double Pin, double *peak, double *SNR)
{
  int i;
  double fin;

  fin = 1./Pin;
  i = findX(f, fin, 0, N);
  if(i > 0 && i < N - 1)
    {
      (*peak) = (p[i]*(fin - f[i]) + p[i+1]*(f[i+1] - fin))/(f[i+1] - f[i]);
    }
  else if(i == 0)
    {
      if(fin > f[0])
	{
	  (*peak) = (p[i]*(fin - f[i]) + p[i+1]*(f[i+1] - fin))/(f[i+1] - f[i]);
	}
      else
	(*peak) = p[0];
    }
  else if(i == N)
    (*peak) = p[N];
  (*SNR) = ((*peak) - ave)/rms;
}

int fasper(double *x, double *y, int n, double ofac, double hifac, double *wk1, double *wk2, int nwk, int *nout)
{
  /* Given n data points with abscissas x[0...n-1] (which need not be equally spaced) and ordinates y[0...n-1], and given a desired minimum/maximum frequencies and a frequency step, this routine fills array wk1[0...nwk-1] with a sequence of nout increasing frequencies (not angular frequencies) up to hifac times the average Nyquist frequency, and fills array wk2[0...nwk-1] with the values of the Lomb normalized periodogram at those frequencies. The arrays x and y are not altered. nwk, the dimension of wk1 and wk2, must be large enough for intermediate work space, or an error results. The routine also returns jmax such that wk2[jmax] is the maximum element in wk2, and prob, and estimate of the significance of that maximum against the hypothesis of random noise. A small value of prob indicates that a significant periodic signal is present. */
  void avevar(double *data, int n, double *ave, double *var);
  void realft(double *data, int n, int isign);
  int spread(double y, double *yy, int n, double x, int m);
  int j, k, ndim, nfreq, nfreqt, ngood;
  double ave, ck, ckk, cterm, cwt, den, df, fac, fndim, hc2wt;
  double hs2wt, hypo, pmax, sterm, swt, var, xdif, xmax, xmin;
  //freq = minfreq;
  //*nout = 0;
  //while(freq <= maxfreq)
  //  {
  //    (*nout)++;
  //    freq += freqstep;
  //  }
  *nout = 0.5*ofac*hifac*n;
  nfreqt=ofac*hifac*n*MACC;
  //nfreqt=2.0*MACC*(*nout);
  nfreq=64;
  while (nfreq < nfreqt) nfreq <<= 1;
  ndim = nfreq << 1;
  if (ndim > nwk) {
    return(1);
  }
  avevar(y,n,&ave,&var);
  if (var == 0.0) {
    return(1);
  }
  j = 0;
  ngood = 0;
  while(isnan(x[j]))
    j++;
  xmin=x[j];
  ngood = 1;
  xmax = xmin;
  for (j = j + 1; j < n; j++) {
    if(!isnan(x[j]))
      {
	ngood++;
	if (x[j] < xmin) xmin=x[j];
	if (x[j] > xmax) xmax = x[j];
      }
  }
  xdif = xmax - xmin;
  for (j=0;j<ndim;j++) wk1[j]=wk2[j]=0.0;
  fac=ndim/(xdif*ofac);
  fndim=ndim;
  for (j=0;j<n;j++) {
    if(!isnan(x[j]))
      {
	ck = (x[j] - xmin)*fac;
	MOD(ck,fndim);
	ckk=2.0*(ck++);
	MOD(ckk,fndim);
	++ckk;
	if(spread(y[j]-ave,wk1,ndim,ck,MACC))
	  return(1);
	if(spread(1.0,wk2,ndim,ckk,MACC))
	  return(1);
      }
  }
  realft(wk1,ndim,1);
  realft(wk2,ndim,1);
  df=1.0/(xdif*ofac);
  pmax = -1.0;
  for (k=3,j=1;j<=(*nout);j++,k+=2) {
    hypo=sqrt(wk2[k-1]*wk2[k-1]+wk2[k]*wk2[k]);
    hc2wt=0.5*wk2[k-1]/hypo;
    hs2wt=0.5*wk2[k]/hypo;
    cwt=sqrt(0.5+hc2wt);
    swt=SIGN(sqrt(0.5-hc2wt),hs2wt);
    den=0.5*n+hc2wt*wk2[k-1]+hs2wt*wk2[k];
    cterm=SQR(cwt*wk1[k-1]+swt*wk1[k])/den;
    sterm=SQR(cwt*wk1[k]-swt*wk1[k-1])/(n-den);
    wk1[j-1]=j*df;
    wk2[j-1]=(cterm+sterm)/(2.0*var);
  }
  return(0);
}

int gfasper(double *x, double *y, double *err, int n, double ofac, double hifac, double *wk1, double *wk2, int nwk, int *nout, int *ngood_out)
{
  /* Similar to fasper, here we compute the generalized lomb-scargle periodogram due to Zechmeister and K\"urster 2009, A&A, 496, 577. */
  void avevar(double *data, int n, double *ave, double *var);
  void realft(double *data, int n, int isign);
  int spread(double y, double *yy, int n, double x, int m);
  int j, k, ndim, nfreq, nfreqt, ngood;
  double ave, ck, ckk, cterm, cwt, den, df, fac, fndim, hc2wt;
  double hs2wt, hypo1, hypo2, hypofrac;
  double pmax, sterm, swt, var, xdif, xmax, xmin, wsum;
  double ysum, yhatsum, YY, YC, YS, CC, SS, CS, D, Y, C, S, YYhat, YChat;
  double YShat, CChat, SShat, CShat;
  double whc2wt, whs2wt, wcsqrwt, wcswt, wcwt, wswt, wycwt, wyswt;

  double *w = NULL;
  double *wk3 = NULL;

  if(n > 0) {
    if((w = (double *) malloc(n * sizeof(double))) == NULL) error(ERR_MEMALLOC);
  }


  //freq = minfreq;
  //*nout = 0;
  //while(freq <= maxfreq)
  //  {
  //    (*nout)++;
  //    freq += freqstep;
  //  }
  *nout = 0.5*ofac*hifac*n;
  *ngood_out = 0;
  nfreqt=ofac*hifac*n*MACC;
  //nfreqt=2.0*MACC*(*nout);
  nfreq=64;
  while (nfreq < nfreqt) nfreq <<= 1;
  ndim = nfreq << 1;
  if (ndim > nwk) {
    if(w != NULL)
      free(w);
    return(1);
  }

  if((wk3 = (double *) malloc(nwk * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  j = 0;
  ngood = 0;
  wsum = 0.;
  while(isnan(x[j]) || isnan(err[j]) || err[j] <= 0.) {
    j++;
  }
  if(j >= n) {
    if(w != NULL)
      free(w);
    return(1);
  }
  w[j] = 1./err[j]/err[j];
  wsum = w[j];
  xmin=x[j];
  ngood = 1;
  xmax = xmin;
  for (j = j + 1; j < n; j++) {
    if(!isnan(x[j]) && !isnan(err[j]) && err[j] > 0.)
      {
	w[j] = 1./err[j]/err[j];
	wsum += w[j];
	ngood++;
	if (x[j] < xmin) xmin=x[j];
	if (x[j] > xmax) xmax = x[j];
      }
  }
  *ngood_out = ngood;
  xdif = xmax - xmin;
  for (j=0;j<ndim;j++) wk1[j]=wk2[j]=wk3[j]=0.0;
  fac=ndim/(xdif*ofac);
  fndim=ndim;
  Y = 0.;
  YYhat = 0.;
  for (j=0;j<n;j++) {
    if(!isnan(x[j]) && !isnan(err[j]) && err[j] > 0.)
      {
	w[j] = w[j]/wsum;
	Y += (y[j]*w[j]);
	YYhat += (y[j]*y[j]*w[j]);
	ck = (x[j] - xmin)*fac;
	MOD(ck,fndim);
	ckk=2.0*(ck++);
	MOD(ckk,fndim);
	++ckk;
	if(spread(w[j]*y[j],wk1,ndim,ck,MACC)) {
	  if(w != NULL)
	    free(w);
	  if(wk3 != NULL)
	    free(wk3);
	  return(1);
	}
	if(spread(w[j],wk2,ndim,ck,MACC)) {
	  if(w != NULL)
	    free(w);
	  if(wk3 != NULL)
	    free(wk3);
	  return(1);
	}
	if(spread(w[j],wk3,ndim,ckk,MACC)) {
	  if(w != NULL)
	    free(w);
	  if(wk3 != NULL)
	    free(wk3);
	  return(1);
	}
      }
  }
  YY = YYhat - (Y*Y);
  if(YY <= 0.) {
    if(w != NULL)
      free(w);
    if(wk3 != NULL)
      free(wk3);
    return(1);
  }
  realft(wk1,ndim,1);
  realft(wk2,ndim,1);
  realft(wk3,ndim,1);
  df=1.0/(xdif*ofac);
  pmax = -1.0;
  for (k=3,j=1;j<=(*nout);j++,k+=2) {
    //hypo1=sqrt(wk2[k-1]*wk2[k-1]+wk2[k]*wk2[k]);
    //hypo2=sqrt(wk3[k-1]*wk3[k-1]+wk3[k]*wk3[k]);
    //hypofrac = hypo1/hypo2;
    C = wk2[k-1];
    S = wk2[k];
    YChat = wk1[k-1];
    YShat = wk1[k];
    CChat = 0.5*(1.0 + wk3[k-1]);
    SShat = 0.5*(1.0 - wk3[k-1]);
    CShat = 0.5*wk3[k];
    YC = YChat - Y*C;
    YS = YShat - Y*S;
    CC = CChat - C*C;
    SS = SShat - S*S;
    CS = CShat - C*S;
    D = CC*SS - CS*CS;
    wk1[j-1] = j*df;
    wk2[j-1] = (SS*YC*YC + CC*YS*YS - 2.0*CS*YC*YS)/(YY*D);
  }
  if(w != NULL)
    free(w);
  if(wk3 != NULL)
    free(wk3);
  return(0);
}

#undef SIGN

int spread(double y, double *yy, int n, double x, int m)
{
  /* Given an array yy[0...n-1], extirpolate (spread) a value y into m actual array elements that best approximate the "fictional" (i.e., possibly noninteger) array element number x. The weights used are coefficients of the Lagrange interpolating polynomial.*/
  int ihi, ilo, ix, j, nden;
  static int nfac[11]={0,1,1,2,6,24,120,720,5040,40320,362880};
  double fac;
  if (m > 10)
    return(1);
  ix = (int)x;
  if (x == (double)ix) yy[ix-1] += y;
  else {
    ilo=MIN_(MAX_((x-0.5*m+1.0),1),n-m+1);
    ihi=ilo+m-1;
    nden=nfac[m];
    fac=x-ilo;
    for (j=ilo+1;j<=ihi;j++) fac *= (x-j);
    yy[ihi-1] += y*fac/(nden*(x-ihi));
    for (j=ihi-1;j>=ilo;j--) {
      nden=(nden/(j+1-ilo))*(j-ihi);
      yy[j-1] += y*fac/(nden*(x-j));
    }
  }
  return(0);
}

void realft(double *data, int n, int isign)
{
  /* Calculates the Fourier transform of a set of n real-valued data points. Replaces this data (which is stored in array data[0...n-1]) by the positive frequency half of its complex Fourier transform. The real-valued first and last components of the complex transform are returned as elements data[0] and data[1], respectively. n must be a power of 2. This routine also calculates the inverse transform of a complex data array if it is the transform of real data (Result in this case must be multipled by 2/n.)*/
  void four1(double *data, int nn, int isign);
  int i, i1, i2, i3, i4, np3;
  double c1=0.5, c2, h1r, h1i, h2r, h2i;
  double wr, wi, wpr, wpi, wtemp, theta;

  theta = 3.141592653589793/(double) (n>>1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data,n>>1,1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np3=n+3;
  for(i=2;i<=(n>>2);i++){
    i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
    h1r=c1*(data[i1-1]+data[i3-1]);
    h1i=c1*(data[i2-1]-data[i4-1]);
    h2r = -c2*(data[i2-1]+data[i4-1]);
    h2i=c2*(data[i1-1]-data[i3-1]);
    data[i1-1]=h1r+wr*h2r-wi*h2i;
    data[i2-1]=h1i+wr*h2i+wi*h2r;
    data[i3-1]=h1r-wr*h2r+wi*h2i;
    data[i4-1]=-h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[0] = (h1r=data[1])+data[1];
    data[1] = h1r-data[1];
  } else {
    data[0] = c1*((h1r=data[0])+data[1]);
    data[1]=c1*(h1r-data[1]);
    four1(data,n>>1,-1);
  }
}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(double *data, int nn, int isign)
{
  /* Replaces data[0...(2*nn)-1] by its discrete Fourier transform, if isign is input as 1; or replaces data[0....(2*nn)-1] by nn times its inverse discrete Fourier transform, if isign is input as -1. data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn MUST be an integer power of 2 (this is not checked for!). */
  int n, mmax, m, j, istep, i;
  double wtemp, wr, wpr, wpi, wi, theta;
  double tempr, tempi;

  n=nn << 1;
  j=1;
  for(i=1;i<n;i+=2) {
    if(j > i) {
      SWAP(data[j-1],data[i-1]);
      SWAP(data[j],data[i]);
    }
    m=nn;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep = mmax << 1;
    theta = isign*(6.28318530717959/mmax);
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for(m=1;m<mmax;m+=2) {
      for(i=m;i<=n;i+=istep) {
	j = i + mmax;
	tempr = wr*data[j-1]-wi*data[j];
	tempi = wr*data[j]+wi*data[j-1];
	data[j-1] = data[i-1]-tempr;
	data[j] = data[i] - tempi;
	data[i-1] += tempr;
	data[i] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}

void avevar(double *data, int n, double *ave, double *var)
{
  /* Given array data[0...n-1], returns its mean as ave and its variance as var. */
  int j, n_;
  double s,ep;
  for (*ave=0.0,j=0, n_=0;j<n;j++) {
    if(!isnan(data[j]))
      {
	n_++;
	*ave += data[j];
      }
  }
  *ave /= n_;
  *var=ep=0.0;
  for (j=0;j<n;j++) {
    if(!isnan(data[j]))
      {
	s=data[j]-(*ave);
	ep += s;
	*var += s*s;
      }
  }
  *var = (*var - ep*ep/n_)/(n_-1);
}

/* Routine that will run LS on a given light curve */
void Lombscargle (int N_in, double *t_in, double *mag_in, double *sig_in, double minper, double maxper, double subsample, int Npeaks, double *periods, double *peaks, double *probs, double *SNR, int outputflag, char *outfile, int ascii, int whiten, double clip, int clipiter, int fixperiodSNR, double fixperiodSNR_period, double *fixperiodSNR_FAPvalues, double *fixperiodSNR_SNRvalues, double *fixperiodSNR_peakvalues, int use_orig_ls, int dobootstrapfap, int Nbootstrap, int usemask, _Variable *maskvar, int lcindex, int threadindex)
{
  int nf, nfreqt, nfreq, ndim, outval, nout, i, j, k, foundflag, peakiter, test, nclippedlast, nclippedthis, ngood, bootstrapind;
  double T, freq, freqstep, minfreq, *wk1, *wk2, maxfreq, trueminfreq, expy, effm, ofac, hifac, testperiod, lastpoint, tmpprob;
  double **wk2_whiten, *t_cpy, *mag_cpy, *sig_cpy, fundA, fundB, meanval, amp, zval;
  double *pbest_whiten;
  double *t_cpy2 = NULL, *mag_cpy2 = NULL, *sig_cpy2 = NULL;
  double *t_mask = NULL, *mag_mask = NULL, *sig_mask = NULL;
  double *bootstrapdist = NULL, *bootstrapprobs = NULL, *fitcoeffs = NULL;
  long klong;

  long double Sum, Sumsqr;
  double ave, rms, val, pbest;
  int nfit;

  FILE *outf;
  int N;
  double *t, *mag, *sig;

  if(!usemask) {
    N = N_in;
    t = t_in;
    mag = mag_in;
    sig = sig_in;
  } else {
    if((t_mask = (double *) malloc(N_in*sizeof(double))) == NULL ||
       (mag_mask = (double *) malloc(N_in*sizeof(double))) == NULL ||
       (sig_mask = (double *) malloc(N_in*sizeof(double))) == NULL) {
      error(ERR_MEMALLOC);
    }
    N = 0;
    for(i = 0; i < N_in; i++) {
      if(!isnan(mag_in[i]) && EvaluateVariable_Double(lcindex, threadindex, i, maskvar) > VARTOOLS_MASK_TINY) {
	t_mask[N] = t_in[i];
	mag_mask[N] = mag_in[i];
	sig_mask[N] = sig_in[i];
	N++;
      }
    }
    if(N <= 1) {
      if(t_mask != NULL) free(t_mask);
      if(mag_mask != NULL) free(mag_mask);
      if(sig_mask != NULL) free(sig_mask);
      return;
    }
    t = t_mask;
    mag = mag_mask;
    sig = sig_mask;
  }
	

  T = t[N-1] - t[0];

  nf = 0;
  freq = 1./minper;
  maxfreq = freq;
  trueminfreq = 1./maxper;
  minfreq = subsample / T;
  freqstep = subsample / T;

  ofac = 1./subsample;
  hifac = maxfreq * (2.0  * T / N);

  nout = 0.5*ofac*hifac*N;
  nfreqt=ofac*hifac*N*MACC;
  nfreq=64;
  while (nfreq < nfreqt) nfreq <<= 1;
  ndim = nfreq << 1;

  /* What we do depends on whether or not the light curve will be whitened between peaks */
  if(!whiten)
    {
      if((wk1 = (double *) malloc(ndim * sizeof(double))) == NULL ||
	 (wk2 = (double *) malloc(ndim * sizeof(double))) == NULL)
	{
	  fprintf(stderr,"Memory Allocation Error\n");
	  exit(3);
	}

      if(dobootstrapfap) {
	if((bootstrapdist = (double *) malloc(Nbootstrap * sizeof(double))) == NULL ||
	   (bootstrapprobs = (double *) malloc(Nbootstrap * sizeof(double))) == NULL ||
	   (fitcoeffs = (double *) malloc(2 * sizeof(double))) == NULL ||
	   (t_cpy2 = (double *) malloc(N * sizeof(double))) == NULL ||
	   (mag_cpy2 = (double *) malloc(N * sizeof(double))) == NULL ||
	   (sig_cpy2 = (double *) malloc(N * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
	
	memcpy(t_cpy2,t,N*sizeof(double));
	for(i=0; i < Nbootstrap; i++) {
	  for(j=0; j < N; j++) {
	    klong = randlong(N-1);
	    mag_cpy2[j] = mag[klong];
	    sig_cpy2[j] = sig[klong];
	  }
	  if(use_orig_ls) {
	    outval = fasper(t_cpy2, mag_cpy2, N, ofac, hifac, wk1, wk2, ndim, &nout);
	  } else {
	    outval = gfasper(t_cpy2, mag_cpy2, sig_cpy2, N, ofac, hifac, wk1, wk2, ndim, &nout, &ngood);
	  }
	  bootstrapdist[i] = -1.;
	  for(k=0;k<nout;k++)
	    {
	      if(wk1[k] >= minfreq)
		{
		  if(wk2[k] > bootstrapdist[i]) bootstrapdist[i] = wk2[k];
		}
	    }
	}
	mysort1(Nbootstrap, bootstrapdist);
	for(k=0; k < Nbootstrap; k++) {
	  bootstrapprobs[k] = log10(((Nbootstrap - k)/((double) Nbootstrap)));
	}
	if(use_orig_ls) {
	  nfit = floor(0.1*Nbootstrap);
	  if(nfit < 1) nfit = 1;
	  fitpoly(nfit, &(bootstrapdist[Nbootstrap-nfit]), &(bootstrapprobs[Nbootstrap-nfit]), NULL, 1, 0, fitcoeffs, NULL);
	}
	else {
	  for(k = 0; k < Nbootstrap; k++) {
	    bootstrapdist[k] = log10(1.0 - bootstrapdist[k]);
	  }
	  nfit = floor(0.1*Nbootstrap);
	  if(nfit < 1) nfit = 1;
	  fitpoly(nfit, &(bootstrapdist[Nbootstrap-nfit]), &(bootstrapprobs[Nbootstrap-nfit]), NULL, 1, 0, fitcoeffs, NULL);
	  for(k = 0; k < Nbootstrap; k++) {
	    bootstrapdist[k] = 1.0 - pow(10.0,bootstrapdist[k]);
	  }
	}
      }


      if(use_orig_ls) {
	outval = fasper(t, mag, N, ofac, hifac, wk1, wk2, ndim, &nout);
      } else {
	outval = gfasper(t, mag, sig, N, ofac, hifac, wk1, wk2, ndim, &nout, &ngood);
      }
      if(!use_orig_ls) {
	pbest = -1.;
	for(i=0;i<nout;i++)
	  {
	    if(wk1[i] >= minfreq)
	      {
		if(wk2[i] > pbest) pbest = wk2[i];
	      }
	  }
      }
      if(outval)
	{
	  for(j=0;j<Npeaks;j++)
	    {
	      periods[j] = -1.;
	      peaks[j] = -1.;
	      probs[j] = -1.;
	      SNR[j] = -1.;
	    }
	}
      else
	{
	  minfreq = 1./maxper;
	  if(outputflag)
	    {
	      if((outf = fopen(outfile,"w")) == NULL)
		{
		  fprintf(stderr,"Error: Cannot Write to File %s\n",outfile);
		  exit(5);
		}
	      if(ascii)
		{
		  if(use_orig_ls) {
		    fprintf(outf,"# Column 1 = Frequency in cycles per input light curve time unit.\n");
		    fprintf(outf,"# Column 2 = Normalized L-S periodogram (equation 13.8.4 from\n");
		    fprintf(outf,"#            Press et al. 1992).\n");
		    fprintf(outf,"# Column 3 = Logarithm of the false alarm probability.\n");
		  }
		  else {
		    fprintf(outf,"# Column 1 = Frequency in cycles per input light curve time unit.\n");
		    fprintf(outf,"# Column 2 = Unnormalized P(omega) (equation 5 of Zechmeister &\n");
		    fprintf(outf,"#            K\\\"urster 2009, A&A, 496, 577).\n");
		    fprintf(outf,"# Column 3 = Logarithm of the false alarm probability.\n");
		  }
		  if(use_orig_ls) {
		    effm = log10(2.0 * maxfreq * T);
		    for(i=0;i<nout;i++)
		      {
			if(wk1[i] >= minfreq)
			  {
			    if(!dobootstrapfap) {
			      expy = -wk2[i]*0.434294482;
			      tmpprob=effm + expy;
			      if(tmpprob > -2) {
				expy = exp(-wk2[i]);
				tmpprob=1.0-pow(1.0-expy,2.0 * maxfreq * T);
				tmpprob = log10(tmpprob);
			      }
			    } else {
			      bootstrapind = findX(bootstrapdist,wk2[i],0,Nbootstrap);
			      if(bootstrapind == Nbootstrap)
				tmpprob = fitcoeffs[0] + fitcoeffs[1]*wk2[i];
			      else {
				tmpprob = (Nbootstrap - bootstrapind)/((double) Nbootstrap);
				tmpprob = log10(tmpprob);
			      }
			    }
			    fprintf(outf,"%.17g %.17g %.17g\n",wk1[i],wk2[i], tmpprob);
			  }
		      }
		  } else {
		    effm = log10(2.0 * maxfreq * T);
		    for(i=0;i<nout;i++)
		      {
			if(wk1[i] >= minfreq)
			  {
			    if(!dobootstrapfap) {
			      zval = (ngood - 3)*(wk2[i]/(1.0 - pbest))/2.0;
			      expy = -(ngood - 3)*log10((1.0 + (2.0*zval)/(ngood - 3)))/2.0;
			      tmpprob=effm + expy;
			      if(tmpprob > -2) {
				expy = pow(((1.0 + (2.0*zval)/(ngood - 3))),(-(ngood-3)/2.0)); 
				tmpprob=1.0-pow(1.0-expy,2.0 * maxfreq * T);
				tmpprob = log10(tmpprob);
			      }
			    } else {
			      bootstrapind = findX(bootstrapdist,wk2[i],0,Nbootstrap);
			      if(bootstrapind == Nbootstrap)
				tmpprob = fitcoeffs[0] + fitcoeffs[1]*log10(1.0 - wk2[i]);
			      else {
				tmpprob = (Nbootstrap - bootstrapind)/((double) Nbootstrap);
				tmpprob = log10(tmpprob);
			      }
			      
			    }

			    fprintf(outf,"%.17g %.17g %.17g\n",wk1[i],wk2[i], tmpprob);
			  }
		      }
		  }
		}
	      else
		{
		  fwrite(&nout,4,1,outf);
		  fwrite(wk1,8,nout,outf);
		  fwrite(wk2,8,nout,outf);
		}
	      fclose(outf);
	    }
	  Sum = 0.;
	  Sumsqr = 0.;
	  for(i=0;i<nout;i++)
	    {
	      Sum += (long double) wk2[i];
	      Sumsqr += (long double) (wk2[i]*wk2[i]);
	    }
	  ave = (double) (Sum / ((long double) nout));
	  rms = sqrt(((double) (Sumsqr / ((long double) nout))) - (ave*ave));
	  nclippedthis = 0;
	  do {
	    nclippedlast = nclippedthis;
	    nclippedthis = 0;
	    Sum = 0.;
	    Sumsqr = 0.;
	    for(i=0;i<nout;i++)
	      {
		if(wk2[i] < ave + clip*rms)
		  {
		    Sum += (long double) wk2[i];
		    Sumsqr += (long double) (wk2[i]*wk2[i]);
		  }
		else
		  nclippedthis++;
	      }
	    ave = (double) (Sum / ((long double) (nout - nclippedthis)));
	    rms = sqrt(((double) (Sumsqr / ((long double) (nout - nclippedthis)))) - (ave*ave));
	  } while(clipiter && nclippedthis > nclippedlast);
	  /* Get the FAP and SNR for the fixed period if we're doing that */
	  if(fixperiodSNR)
	    {
	      LS_getfixedperiodSNR(wk1,wk2,nout,ave,rms,fixperiodSNR_period,fixperiodSNR_peakvalues,fixperiodSNR_SNRvalues);
	    }
	  /* Replace the periodogram with only points that are local maxima */
	  for(k=0; k < nout; k++) {if(wk1[k] > trueminfreq) break;}
	  if(k == nout) k = nout-1;
	  lastpoint = wk2[k] - 1.;
	  for(i=0;k<nout-1;k++)
	    {
	      if(wk2[k] > lastpoint && wk2[k] > wk2[k + 1])
		{
		  lastpoint = wk2[k];
		  wk2[i] = wk2[k];
		  wk1[i] = wk1[k];
		  i++;
		}
	      else
		{
		  lastpoint = wk2[k];
		}
	    }
	  if(wk2[k] > lastpoint)
	    {
	      wk2[i] = wk2[k];
	      wk1[i] = wk1[k];
	      i++;
	    }
	  nout = i;
	  mysort2(nout,wk2,wk1);
	  effm = log10(2.0 * maxfreq * T);
	  if(fixperiodSNR)
	    {
	      val = *fixperiodSNR_peakvalues;
	      if(use_orig_ls) {
		expy = -val*0.434294482;
		*fixperiodSNR_FAPvalues = effm + expy;
		if((*fixperiodSNR_FAPvalues) > -2.) {
		  expy = exp(-val);
		  *fixperiodSNR_FAPvalues = 1.0 - pow(1.0 - expy, 2.0 * maxfreq * T);
		  *fixperiodSNR_FAPvalues = log10(*fixperiodSNR_FAPvalues);
		}
	      } else {
		zval = (ngood - 3)*(val/(1.0 - pbest))/2.0;
		expy = -(ngood - 3)*log10((1.0 + (2.0*zval)/(ngood - 3)))/2.0;
		*fixperiodSNR_FAPvalues = effm + expy;
		if((*fixperiodSNR_FAPvalues) > -2.) {
		  expy = pow(((1.0 + (2.0*zval)/(ngood - 3))),(-(ngood-3)/2.0));
		  *fixperiodSNR_FAPvalues = 1.0 - pow(1.0 - expy, 2.0 * maxfreq * T);
		  *fixperiodSNR_FAPvalues = log10(*fixperiodSNR_FAPvalues);
		}		
	      }
	    }
	  for(j=0,i=nout-1;i>=0 && j < Npeaks;i--)
	    {
	      if(wk1[i] > trueminfreq)
		{
		  if(!j)
		    {
		      periods[j] = 1./wk1[i];
		      //expy=exp(-wk2[i]);
		      //effm=2.0 * maxfreq * T;
		      if(use_orig_ls) {
			expy = -wk2[i]*0.434294482;
			peaks[j] = wk2[i];
			SNR[j] = (wk2[i] - ave)/rms;
			if(!dobootstrapfap) {
			  probs[j]=effm + expy;
			  if (probs[j] > -2) {
			    expy = exp(-wk2[i]);
			    probs[j]=1.0-pow(1.0-expy,2.0 * maxfreq * T);
			    probs[j] = log10(probs[j]);
			  }
			} else {
			  bootstrapind = findX(bootstrapdist,wk2[i],0,Nbootstrap);
			  if(bootstrapind == Nbootstrap)
			    probs[j] = fitcoeffs[0] + fitcoeffs[1]*wk2[i];
			  else {
			    probs[j] = (Nbootstrap - bootstrapind)/((double) Nbootstrap);
			    probs[j] = log10(probs[j]);
			  }
			}
		      } else {
			zval = (ngood - 3)*(wk2[i]/(1.0 - pbest))/2.0;
			expy = -(ngood - 3)*log10((1.0 + (2.0*zval)/(ngood - 3)))/2.0;
			peaks[j] = wk2[i];
			SNR[j] = (wk2[i] - ave)/rms;
			if(!dobootstrapfap) {
			  probs[j]=effm + expy;
			  if(probs[j] > -2.) {
			    expy = pow(((1.0 + (2.0*zval)/(ngood - 3))),(-(ngood-3)/2.0));
			    probs[j]=1.0-pow(1.0-expy,2.0 * maxfreq * T);
			    probs[j] = log10(probs[j]);
			  }
			} else {
			  bootstrapind = findX(bootstrapdist,wk2[i],0,Nbootstrap);
			  if(bootstrapind == Nbootstrap)
			    probs[j] = fitcoeffs[0] + fitcoeffs[1]*log10(1.0 - wk2[i]);
			  else {
			    probs[j] = (Nbootstrap - bootstrapind)/((double) Nbootstrap);
			    probs[j] = log10(probs[j]);
			  }

			}
		      }
		      j++;
		    }
		  else
		    {
		      foundflag = 0;
		      testperiod = 1./wk1[i];
		      for(k=0;k<j && !foundflag;k++)
			{
			  if(!isDifferentPeriods(periods[k],testperiod,T))
			    foundflag = 1;
			}
		      if(!foundflag)
			{
			  periods[j] = testperiod;
			  if(use_orig_ls) {
			    expy = -wk2[i]*0.434294482;
			    peaks[j] = wk2[i];
			    SNR[j] = (wk2[i] - ave)/rms;
			    if(!dobootstrapfap) {
			      probs[j]=effm + expy;
			      if (probs[j] > -2) {
				expy = exp(-wk2[i]);
				probs[j]=1.0-pow(1.0-expy,2.0 * maxfreq * T);
				probs[j] = log10(probs[j]);
			      }
			    } else {
			      bootstrapind = findX(bootstrapdist,wk2[i],0,Nbootstrap);
			      if(bootstrapind == Nbootstrap)
				probs[j] = fitcoeffs[0] + fitcoeffs[1]*wk2[i];
			      else {
				probs[j] = (Nbootstrap - bootstrapind)/((double) Nbootstrap);
				probs[j] = log10(probs[j]);
			      }
			    }
			  } else {
			    zval = (ngood - 3)*(wk2[i]/(1.0 - pbest))/2.0;
			    expy = -(ngood - 3)*log10((1.0 + (2.0*zval)/(ngood - 3)))/2.0;
			    peaks[j] = wk2[i];
			    SNR[j] = (wk2[i] - ave)/rms;
			    if(!dobootstrapfap) {
			      probs[j]=effm + expy;
			      if(probs[j] > -2.) {
				expy = pow(((1.0 + (2.0*zval)/(ngood - 3))),(-(ngood-3)/2.0));
				probs[j]=1.0-pow(1.0-expy,2.0 * maxfreq * T);
				probs[j] = log10(probs[j]);
			      }
			    } else {
			      bootstrapind = findX(bootstrapdist,wk2[i],0,Nbootstrap);
			      if(bootstrapind == Nbootstrap)
				probs[j] = fitcoeffs[0] + fitcoeffs[1]*log10(1.0 - wk2[i]);
			      else {
				probs[j] = (Nbootstrap - bootstrapind)/((double) Nbootstrap);
				probs[j] = log10(probs[j]);
			      }
			    }
			  }
			  j++;
			}
		    }
		}
	    }
	  for(;j<Npeaks;j++)
	    {
	      periods[j] = -1.;
	      peaks[j] = -1.;
	      probs[j] = 1.;
	      SNR[j] = -1.;
	    }

	}
      free(wk1);
      free(wk2);
    }
  else
    {
      /* We are whitening the light curve and re-computing the periodogram after finding each peak */
      if((wk1 = (double *) malloc(ndim * sizeof(double))) == NULL ||
	 (wk2_whiten = (double **) malloc((Npeaks + 1) * sizeof(double *))) == NULL ||
	 (pbest_whiten = (double *) malloc((Npeaks + 1) * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0;i<Npeaks+1;i++)
	{
	  if((wk2_whiten[i] = (double *) malloc(ndim * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}


      if(dobootstrapfap) {
	if((bootstrapdist = (double *) malloc(Nbootstrap * sizeof(double))) == NULL ||
	   (bootstrapprobs = (double *) malloc(Nbootstrap * sizeof(double))) == NULL ||
	   (fitcoeffs = (double *) malloc(2 * sizeof(double))) == NULL ||
	   (t_cpy2 = (double *) malloc(N * sizeof(double))) == NULL ||
	   (mag_cpy2 = (double *) malloc(N * sizeof(double))) == NULL ||
	   (sig_cpy2 = (double *) malloc(N * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
	
	memcpy(t_cpy2,t,N*sizeof(double));
	for(i=0; i < Nbootstrap; i++) {
	  for(j=0; j < N; j++) {
	    klong = randlong(N-1);
	    mag_cpy2[j] = mag[klong];
	    sig_cpy2[j] = sig[klong];
	  }
	  if(use_orig_ls) {
	    outval = fasper(t_cpy2, mag_cpy2, N, ofac, hifac, wk1, wk2_whiten[0], ndim, &nout);
	  } else {
	    outval = gfasper(t_cpy2, mag_cpy2, sig_cpy2, N, ofac, hifac, wk1, wk2_whiten[0], ndim, &nout, &ngood);
	  }
	  bootstrapdist[i] = -1.;
	  for(k=0;k<nout;k++)
	    {
	      if(wk1[k] >= minfreq)
		{
		  if(wk2_whiten[0][k] > bootstrapdist[i]) bootstrapdist[i] = wk2_whiten[0][k];
		}
	    }
	}
	mysort1(Nbootstrap, bootstrapdist);
	for(k=0; k < Nbootstrap; k++) {
	  bootstrapprobs[k] = log10(((Nbootstrap - k)/((double) Nbootstrap)));
	}
	if(use_orig_ls) {
	  nfit = floor(0.1*Nbootstrap);
	  if(nfit < 1) nfit = 1;
	  fitpoly(nfit, &(bootstrapdist[Nbootstrap-nfit]), &(bootstrapprobs[Nbootstrap-nfit]), NULL, 1, 0, fitcoeffs, NULL);
	}
	else {
	  for(k = 0; k < Nbootstrap; k++) {
	    bootstrapdist[k] = log10(1.0 - bootstrapdist[k]);
	  }
	  nfit = floor(0.1*Nbootstrap);
	  if(nfit < 1) nfit = 1;
	  fitpoly(nfit, &(bootstrapdist[Nbootstrap-nfit]), &(bootstrapprobs[Nbootstrap-nfit]), NULL, 1, 0, fitcoeffs, NULL);
	  for(k = 0; k < Nbootstrap; k++) {
	    bootstrapdist[k] = 1.0 - pow(10.0,bootstrapdist[k]);
	  }
	}
	  
      }

      /* Make a copy of the light curve */
      if((t_cpy = (double *) malloc(N * sizeof(double))) == NULL ||
	 (mag_cpy = (double *) malloc(N * sizeof(double))) == NULL ||
	 (sig_cpy = (double *) malloc(N * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      memcpy(t_cpy,t,N*sizeof(double));
      memcpy(mag_cpy,mag,N*sizeof(double));
      memcpy(sig_cpy,sig,N*sizeof(double));

      effm = log10(2.0 * maxfreq * T);
      /* Compute the periodogram */
      if(use_orig_ls) {
	outval = fasper(t_cpy, mag_cpy, N, ofac, hifac, wk1, wk2_whiten[0], ndim, &nout);
      } else {
	outval = gfasper(t_cpy, mag_cpy, sig_cpy, N, ofac, hifac, wk1, wk2_whiten[0], ndim, &nout, &ngood);
      }
      if(!use_orig_ls) {
	pbest_whiten[0] = -1.;
	for(i=0;i<nout;i++)
	  {
	    if(wk1[i] >= minfreq)
	      {
		if(wk2_whiten[0][i] > pbest_whiten[0]) pbest_whiten[0] = wk2_whiten[0][i];
	      }
	  }
      }
      if(outval)
	{
	  for(j=0;j<Npeaks;j++)
	    {
	      periods[j] = -1.;
	      peaks[j] = -1.;
	      probs[j] = -1.;
	      SNR[j] = -1.;
	      free(wk2_whiten[i]);
	    }
	  free(wk2_whiten[i]);
	  free(wk1);
	  free(wk2_whiten);
	  free(t_cpy);
	  free(mag_cpy);
	  free(sig_cpy);
	  free(pbest_whiten);
	  return;
	}
      for(peakiter=0;peakiter < Npeaks;peakiter++)
	{
      	  /* Find the peak period */
	  minfreq = 1./maxper;

	  /* Find the peak period */
	  periods[peakiter] = -1.;
	  peaks[peakiter] = -1.;
	  probs[peakiter] = -1.;
	  SNR[peakiter] = -1.;
	  for(i=0; i < nout; i++)
	    {
	      if(wk1[i] > trueminfreq)
		{
		  if(wk2_whiten[peakiter][i] > probs[peakiter])
		    {
		      test = 1;
		      testperiod = 1./wk1[i];
		      for(j=0;j<peakiter;j++)
			{
			  if(periods[j] > 0.)
			    {
			      if(!isDifferentPeriods(periods[j],testperiod,T))
				{
				  test = 0;
				  break;
				}
			    }
			}
		      if(test)
			{
			  periods[peakiter] = testperiod;
			  probs[peakiter] = wk2_whiten[peakiter][i];
			}
		    }
		}
	    }
	  /* Whiten the light curve at this period */
	  dokillharms(N, t_cpy, mag_cpy, sig_cpy, 1, &(periods[peakiter]), 0, 0, NULL, NULL, NULL, NULL, &fundA, &fundB, &meanval, 0, NULL, &amp, 0, KILLHARM_OUTTYPE_DEFAULT, -1.);

	  /* Compute the periodogram of the whitened light curve */
	  if(use_orig_ls) {
	    outval = fasper(t_cpy, mag_cpy, N, ofac, hifac, wk1, wk2_whiten[peakiter + 1], ndim, &nout);
	  } else {
	    outval = gfasper(t_cpy, mag_cpy, sig_cpy, N, ofac, hifac, wk1, wk2_whiten[peakiter + 1], ndim, &nout, &ngood);
	  }
	  if(!use_orig_ls) {
	    pbest_whiten[peakiter+1] = -1.;
	    for(i=0;i<nout;i++)
	      {
		if(wk1[i] >= minfreq)
		  {
		    if(wk2_whiten[peakiter+1][i] > pbest_whiten[peakiter+1]) pbest_whiten[peakiter+1] = wk2_whiten[peakiter+1][i];
		  }
	      }
	  }
	  Sum = 0.;
	  Sumsqr = 0.;
	  for(i=0;i<nout;i++)
	    {
	      Sum += (long double) wk2_whiten[peakiter+1][i];
	      Sumsqr += (long double) (wk2_whiten[peakiter+1][i]*wk2_whiten[peakiter+1][i]);
	    }
	  ave = (double) (Sum / ((long double) nout));
	  rms = sqrt(((double) (Sumsqr / ((long double) nout))) - (ave*ave));
	  nclippedthis = 0;
	  do {
	    nclippedlast = nclippedthis;
	    nclippedthis = 0;
	    Sum = 0.;
	    Sumsqr = 0.;
	    for(i=0;i<nout;i++)
	      {
		if(wk2_whiten[peakiter+1][i] < ave + clip*rms)
		  {
		    Sum += (long double) wk2_whiten[peakiter+1][i];
		    Sumsqr += (long double) (wk2_whiten[peakiter+1][i]*wk2_whiten[peakiter+1][i]);
		  }
		else
		  nclippedthis++;
	      }
	    ave = (double) (Sum / ((long double) (nout - nclippedthis)));
	    rms = sqrt(((double) (Sumsqr / ((long double) (nout - nclippedthis)))) - (ave*ave));
	  } while(clipiter && nclippedthis > nclippedlast);


	  /* Get the FAP and SNR for the fixed period if we're doing that */
	  if(fixperiodSNR && !peakiter)
	    {
	      LS_getfixedperiodSNR(wk1,wk2_whiten[peakiter],nout,ave,rms,fixperiodSNR_period,fixperiodSNR_peakvalues,fixperiodSNR_SNRvalues);
	      val = *fixperiodSNR_peakvalues;
	      if(use_orig_ls) {
		expy = -val*0.434294482;
		*fixperiodSNR_FAPvalues = effm + expy;
		if((*fixperiodSNR_FAPvalues) > -2.) {
		  expy = exp(-val);
		  *fixperiodSNR_FAPvalues = 1.0 - pow(1.0 - expy, 2.0 * maxfreq * T);
		  *fixperiodSNR_FAPvalues = log10(*fixperiodSNR_FAPvalues);
		}
	      } else {
		zval = (ngood - 3)*(val/(1.0 - pbest_whiten[peakiter]))/2.0;
		expy = -(ngood - 3)*log10((1.0 + (2.0*zval)/(ngood - 3)))/2.0;
		*fixperiodSNR_FAPvalues=effm + expy;
		if(*fixperiodSNR_FAPvalues > -2.) {
		  expy = pow(((1.0 + (2.0*zval)/(ngood - 3))),(-(ngood-3)/2.0));
		  *fixperiodSNR_FAPvalues=1.0-pow(1.0-expy,2.0 * maxfreq * T);
		  *fixperiodSNR_FAPvalues = log10(*fixperiodSNR_FAPvalues);
		}
	      }
	    }

	  val = probs[peakiter];
	  peaks[peakiter] = val;
	  SNR[peakiter] = (val - ave)/rms;
	  if(use_orig_ls) {
	    expy = -val*0.434294482;
	    if(!dobootstrapfap) {
	      probs[peakiter]=effm + expy;
	      if (probs[peakiter] > -2) {
		expy = exp(-val);
		probs[peakiter]=1.0-pow(1.0-expy,2.0 * maxfreq * T);
		probs[peakiter] = log10(probs[peakiter]);
	      }
	    } else {
	      bootstrapind = findX(bootstrapdist,val,0,Nbootstrap);
	      if(bootstrapind == Nbootstrap)
		probs[peakiter] = fitcoeffs[0] + fitcoeffs[1]*val;
	      else {
		probs[peakiter] = (Nbootstrap - bootstrapind)/((double) Nbootstrap);
		probs[peakiter] = log10(probs[peakiter]);
	      }
	    }
	  } else {
	    zval = (ngood - 3)*(val/(1.0 - pbest_whiten[peakiter]))/2.0;
	    expy = -(ngood - 3)*log10((1.0 + (2.0*zval)/(ngood - 3)))/2.0;
	    if(!dobootstrapfap) {
	      probs[peakiter]=effm + expy;
	      if(probs[peakiter] > -2.) {
		expy = pow(((1.0 + (2.0*zval)/(ngood - 3))),(-(ngood-3)/2.0));
		probs[peakiter]=1.0-pow(1.0-expy,2.0 * maxfreq * T);
		probs[peakiter] = log10(probs[peakiter]);
	      }
	    } else {
	      bootstrapind = findX(bootstrapdist,val,0,Nbootstrap);
	      if(bootstrapind == Nbootstrap)
		probs[peakiter] = fitcoeffs[0] + fitcoeffs[1]*log10(1.0 - val);
	      else {
		probs[peakiter] = (Nbootstrap - bootstrapind)/((double) Nbootstrap);
		probs[peakiter] = log10(probs[peakiter]);
	      }
	    }
	    
	  }
	}

      if(outputflag)
	{
	  if((outf = fopen(outfile,"w")) == NULL)
	    {
	      fprintf(stderr,"Error: Cannot Write to File %s\n",outfile);
	      exit(5);
	    }
	  if(ascii)
	    {
	      if(use_orig_ls) {
		fprintf(outf,"# Column 1 = Frequency in cycles per input light curve time unit.\n");
		for(j=0; j < Npeaks; j++) {
		  fprintf(outf,"# Column %d = Normalized L-S periodogram (equation 13.8.4 from\n", 2+2*j);
		  fprintf(outf,"#            Press et al. 1992). Whitening Cycle %d\n", j);
		  fprintf(outf,"# Column %d = Log10 of the false alarm probability.\n", 3+2*j);
		  fprintf(outf,"#             Whitening Cycle %d.\n", j);
		}
	      }
	      else {
		fprintf(outf,"# Column 1 = Frequency in cycles per input light curve time unit.\n");
		for(j=0; j < Npeaks; j++) {
		  fprintf(outf,"# Column %d = Unnormalized P(omega) (equation 5 of Zechmeister &\n", 2+2*j);
		  fprintf(outf,"#            K\\\"urster 2009, A&A, 496, 577). Whitening Cycle %d.\n", j);
		  fprintf(outf,"# Column %d = Log10 of the false alarm probability.\n", 3+2*j);
		  fprintf(outf,"#             Whitening Cycle %d.\n", j);
		}
	      }
	      if(use_orig_ls) {
		effm = log10(2.0 * maxfreq * T);
		for(i=0;i<nout;i++)
		  {
		    if(wk1[i] >= minfreq)
		      {
			fprintf(outf,"%.17g",wk1[i]);
			for(j=0;j < Npeaks; j++) {
			  if(!dobootstrapfap) {
			    expy = -wk2_whiten[j][i]*0.434294482;
			    tmpprob=effm + expy;
			    if(tmpprob > -2) {
			      expy = exp(-wk2_whiten[j][i]);
			      tmpprob=1.0-pow(1.0-expy,2.0 * maxfreq * T);
			      tmpprob = log10(tmpprob);
			    }
			  } else {
			    bootstrapind = findX(bootstrapdist,wk2_whiten[j][i],0,Nbootstrap);
			    if(bootstrapind == Nbootstrap)
			      tmpprob = fitcoeffs[0] + fitcoeffs[1]*wk2_whiten[j][i];
			    else {
			      tmpprob = (Nbootstrap - bootstrapind)/((double) Nbootstrap);
			      tmpprob = log10(tmpprob);
			    }
			  }

			  fprintf(outf," %.17g %.17g\n",wk2_whiten[j][i],tmpprob);
			}
			fprintf(outf,"\n");
		      }
		  }
	      } else {
		effm = log10(2.0 * maxfreq * T);
		for(i=0;i<nout;i++)
		  {
		    if(wk1[i] >= minfreq)
		      {
			fprintf(outf,"%.17g",wk1[i]);
			for(j=0; j < Npeaks; j++) {
			  if(!dobootstrapfap) {
			    zval = (ngood - 3)*(wk2_whiten[j][i]/(1.0 - pbest_whiten[j]))/2.0;
			    expy = -(ngood - 3)*log10((1.0 + (2.0*zval)/(ngood - 3)))/2.0;
			    tmpprob=effm + expy;
			    if(tmpprob > -2) {
			      expy = pow(((1.0 + (2.0*zval)/(ngood - 3))),(-(ngood-3)/2.0)); 
			      tmpprob=1.0-pow(1.0-expy,2.0 * maxfreq * T);
			      tmpprob = log10(tmpprob);
			    }
			  } else {
			    bootstrapind = findX(bootstrapdist,wk2_whiten[j][i],0,Nbootstrap);
			    if(bootstrapind == Nbootstrap)
			      tmpprob = fitcoeffs[0] + fitcoeffs[1]*log10(1.0 - wk2_whiten[j][i]);
			    else {
			      tmpprob = (Nbootstrap - bootstrapind)/((double) Nbootstrap);
			      tmpprob = log10(tmpprob);
			    }
			  }

			  fprintf(outf," %.17g %.17g",wk2_whiten[j][i], tmpprob);
			}
			fprintf(outf,"\n");
		      }
		  }
	      }
	    }
	  else
	    {
	      fwrite(&nout,4,1,outf);
	      fwrite(&Npeaks,4,1,outf);
	      fwrite(wk1,8,nout,outf);
	      for(j=0;j<Npeaks;j++)
		fwrite(wk2_whiten[j],8,nout,outf);
	    }
	  fclose(outf);
	}
      for(i=0;i<=Npeaks;i++)
	free(wk2_whiten[i]);
      free(wk1);
      free(wk2_whiten);
      free(pbest_whiten);
      free(t_cpy);
      free(mag_cpy);
      free(sig_cpy);

    }

  if(t_mask != NULL) free(t_mask);
  if(mag_mask != NULL) free(mag_mask);
  if(sig_mask != NULL) free(sig_mask);

  if(t_cpy2 != NULL) free(t_cpy2);
  if(mag_cpy2 != NULL) free(mag_cpy2);
  if(sig_cpy2 != NULL) free(sig_cpy2);
  if(bootstrapdist != NULL) free(bootstrapdist);
  if(bootstrapprobs != NULL) free(bootstrapprobs);
  if(fitcoeffs != NULL) free(fitcoeffs);
}

void RunLombScargleCommand(ProgramData *p, _Ls *Ls, Command *c, int lcnum, int lc_name_num, int thisindex)
{
  int i1, i2;
  double d1;
  double *d1ptr, *d2ptr, *d3ptr;
  char outname[MAXLEN];
  if(Ls->operiodogram)
    {
      i1 = 0;
      i2 = 0;
      while(p->lcnames[lc_name_num][i1] != '\0')
	{
	  if(p->lcnames[lc_name_num][i1] == '/')
	    i2 = i1 + 1;
	  i1++;
	}
      sprintf(outname,"%s/%s%s",Ls->outdir,&p->lcnames[lc_name_num][i2],Ls->suffix);
    }
  if(Ls->minp_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    Ls->minp_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Ls->minp_expr);
  }
  else if(Ls->minp_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
    Ls->minp_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Ls->minp_var);
  }
  else {
    Ls->minp_vals[lcnum] = Ls->minp;
  }

  if(Ls->maxp_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    Ls->maxp_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Ls->maxp_expr);
  }
  else if(Ls->maxp_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
    Ls->maxp_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Ls->maxp_var);
  }
  else {
    Ls->maxp_vals[lcnum] = Ls->maxp;
  }

  if(Ls->subsample_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    Ls->subsample_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Ls->subsample_expr);
  }
  else if(Ls->subsample_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
    Ls->subsample_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Ls->subsample_var);
  }
  else {
    Ls->subsample_vals[lcnum] = Ls->subsample;
  }

  if(Ls->fixperiodSNR)
    {
      if(Ls->fixperiodSNR_pertype == PERTYPE_AOV)
	{
	  i1=Ls->fixperiodSNR_lastaovindex;
	  if(c[i1-thisindex].cnum == CNUM_AOV)
	    Ls->fixperiodSNR_periods[lcnum][0] = c[i1-thisindex].Aov->peakperiods[lcnum][0];
	  else if(c[i1-thisindex].cnum == CNUM_HARMAOV)
	    Ls->fixperiodSNR_periods[lcnum][0] = c[i1-thisindex].AovHarm->peakperiods[lcnum][0];
	  
	}
      else if(Ls->fixperiodSNR_pertype == PERTYPE_LS)
	{
	  i1 = Ls->fixperiodSNR_lastaovindex;
	  Ls->fixperiodSNR_periods[lcnum][0] = c[i1-thisindex].Ls->peakperiods[lcnum][0];
	}
      else if(Ls->fixperiodSNR_pertype == PERTYPE_INJECTHARM)
	{
	  i1 = Ls->fixperiodSNR_lastaovindex;
	  Ls->fixperiodSNR_periods[lcnum][0] = c[i1-thisindex].Injectharm->periodinject[lcnum];
	}
      else if(Ls->fixperiodSNR_pertype == PERTYPE_FIX)
	{
	  Ls->fixperiodSNR_periods[lcnum][0] = Ls->fixperiodSNR_fixedperiod;
	}
      else if(Ls->fixperiodSNR_pertype == PERTYPE_FIXCOLUMN)
	{
	  getoutcolumnvalue(Ls->fixperiodSNR_linkedcolumn, lcnum, lc_name_num, VARTOOLS_TYPE_DOUBLE, &(Ls->fixperiodSNR_periods[lcnum][0]));
	}
      if(Ls->fixperiodSNR_pertype != PERTYPE_SPECIFIED)
	d1 = Ls->fixperiodSNR_periods[lcnum][0];
      else
	d1 = Ls->fixperiodSNR_periods[lc_name_num][0];
      d1ptr = &(Ls->fixperiodSNR_FAPvalues[lcnum]);
      d2ptr = &(Ls->fixperiodSNR_SNRvalues[lcnum]);
      d3ptr = &(Ls->fixperiodSNR_peakvalues[lcnum]);
    }
  else
    {
      d1 = 1.;
      d1ptr = NULL;
      d2ptr = NULL;
      d3ptr = NULL;
    }
  if(p->NJD[lcnum] > 1) {
    Lombscargle (p->NJD[lcnum], p->t[lcnum], p->mag[lcnum], p->sig[lcnum], Ls->minp_vals[lcnum], Ls->maxp_vals[lcnum], Ls->subsample_vals[lcnum], Ls->Npeaks, Ls->peakperiods[lcnum], Ls->peakvalues[lcnum], Ls->peakFAP[lcnum], Ls->SNRvalues[lcnum],Ls->operiodogram, outname,p->ascii,Ls->whiten,Ls->clip,Ls->clipiter,Ls->fixperiodSNR,d1,d1ptr,d2ptr,d3ptr,Ls->use_orig_ls,Ls->dobootstrapfap,Ls->Nbootstrap, Ls->usemask, Ls->maskvar, lc_name_num, lcnum);
  }
}
