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

/* This file contains functions for computing the minimum detectable amplitude for a given signal using lomb-scargle for the program vartools by J. Hartman */

#define TWOPID 6.2831853071795865

double ls_oneperiod(double *x, double *y, int n, double period)
{
  /* This routine is based on the "period" function in Numerical Recipes in C by Press et al., it skips the loop over periods though */
  void avevar(double *data, int n, double *ave, double *var);
  int j;
  double ave, c, cc, cwtau, pnow, pymax, s, ss, sumc, sumcy, sums, sumsh, sumsy, swtau, var, wtau, xave, xdif, xmax, xmin, yy;
  double arg, *wi, *wpi, *wpr, *wr;

  wi = (double *) malloc(n*sizeof(double));
  wpi = (double *) malloc(n * sizeof(double));
  wpr = (double *) malloc(n * sizeof(double));
  wr = (double *) malloc(n * sizeof(double));
  avevar(y,n,&ave,&var);
  xmax = xmin=x[0];
  for (j=0;j<n;j++) {
    if(x[j] > xmax) xmax = x[j];
    if(x[j] < xmin) xmin = x[j];
  }
  xdif = xmax - xmin;
  xave = 0.5*(xmax + xmin);
  pymax = 0.0;
  pnow = 1.0/(period);
  for (j=0;j<n;j++) {
    arg = TWOPID*((x[j]-xave)*pnow);
    wpr[j] = -2.0*SQR(sin(0.5*arg));
    wpi[j] = sin(arg);
    wr[j] = cos(arg);
    wi[j] = wpi[j];
  }

  pnow = 1./period;
  sumsh = sumc = 0.0;
  for(j=0;j<n;j++) {
    c = wr[j];
    s = wi[j];
    sumsh += s*c;
    sumc += (c-s)*(c+s);
  }
  wtau = 0.5*atan2(2.0*sumsh,sumc);
  swtau = sin(wtau);
  cwtau = cos(wtau);
  sums=sumc=sumsy=sumcy=0.0;
  for(j=0;j<n;j++) {
    s = wi[j];
    c = wr[j];
    ss = s*cwtau-c*swtau;
    cc = c*cwtau+s*swtau;
    sums += ss*ss;
    sumc += cc*cc;
    yy = y[j] - ave;
    sumsy += yy*ss;
    sumcy += yy*cc;
  }
  pymax=0.5*(sumcy*sumcy/sumc+sumsy*sumsy/sums)/var;

  free(wi);
  free(wpi);
  free(wpr);
  free(wr);

  return pymax;
}

double gls_oneperiod(double *x, double *y, double *err, int n, double period)
{
  /* Similar to ls_oneperiod, this is for the generalized Lomb-Scargle periodogram */
  void avevar(double *data, int n, double *ave, double *var);
  int j;
  double ave, c, cc, cwtau, pnow, pymax, s, ss, sumc, sumcy, sums, sumsh, sumsy, swtau, var, wtau, xave, xdif, xmax, xmin, yy;
  double arg, *w, wsum;

  double Y, C, S, YY, YYhat, YC, YChat, YS, YShat, CC, CChat, SS, SShat;
  double CS, CShat, D, cwt, swt;
  int i;

  Y = 0.; C = 0.; S = 0.;
  YY = 0.; YYhat = 0.; YC = 0.; YChat = 0.; YS = 0.; YShat = 0.;
  CC = 0.; CChat = 0.; SS = 0.; SShat = 0.; CS = 0.; CShat = 0.;
  D = 0.; wsum = 0.;

  if(n <= 0) return -1.;
  if((w = (double *) malloc(n * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  for(i=0; i < n; i ++) {
    if(!isnan(y[i]) && !isnan(err[i]) && err[i] > 0.) {
      w[i] = 1./err[i]/err[i];
      wsum += w[i];
    }
  }
  if(wsum <= 0) {
    free(w);
    return -1.;
  }
  xmax = xmin=x[0];
  for (j=0;j<n;j++) {
    if(x[j] > xmax) xmax = x[j];
    if(x[j] < xmin) xmin = x[j];
  }
  xdif = xmax - xmin;
  xave = 0.5*(xmax + xmin);
  pnow = 1.0/(period);
  for(i=0; i < n; i ++) {
    if(!isnan(y[i]) && !isnan(err[i]) && err[i] > 0.) {
      w[i] = w[i]/wsum;
      Y += w[i]*y[i];
      YYhat += w[i]*y[i]*y[i];
      arg = TWOPID*((x[i]-xave)*pnow);
      cwt = cos(arg);
      swt = sin(arg);
      C += w[i]*cwt;
      S += w[i]*swt;
      YChat += w[i]*y[i]*cwt;
      YShat += w[i]*y[i]*swt;
      CChat += w[i]*cwt*cwt;
      SShat += w[i]*swt*swt;
      CShat += w[i]*cwt*swt;
    }
  }
  YY = YYhat - Y*Y;
  YC = YChat - Y*C;
  YS = YShat - Y*S;
  CC = CChat - C*C;
  SS = SShat - S*S;
  CS = CShat - C*S;
  D = CC*SS - CS*CS;
  pymax = (SS*YC*YC + CC*YS*YS - 2.0*CS*YC*YS)/(YY*D);
  free(w);

  if(pymax < 0.0) pymax = 0.0;
  else if(pymax > 1.0) pymax = 1.0;
  return pymax;
}

#define ITMAX 100
#define EPS 3.0e-8

void getlsampthresh(int N, double *t, double *mag, double *sig, double period, int harm_specsigflag, FILE *signalfile, int Nsubharm, int Nharm, double minPer, double thresh, double *ampthresh_scale, double *amp, int use_orig_ls)
{
  /* This routine takes a light curve and a period, it fits a fourier series to the light curve to determine the signal, which is subtracted from the light curve, it then scales the signal and adds it back to the light curve and calculates the lomb-scargle probability at the specified period. It repeats, adjusting the scale until the LS probability is above the thresh-hold at which point the threshhold scaling factor will be returned. If the original light curve does not pass the signal, a negative factor will be returned.

Embedded in this routine is the zbrent algorithm (see Numerical Recipes in C) to find the zero of LS-Prob = thresh.
*/
  double getlsampthresh_func(double amp_scale, int N, double *t, double *mag_orig, double *sig, double *mag_signal, double *mag_tmp, double period, int Nsubharm, int Nharm, double *subharmA, double *subharmB, double *harmA, double *harmB, double fundA, double fundB, double minPer, double thresh, int use_orig_ls);
  int i;
  double *subharmA, *subharmB, *harmA, *harmB, fundA, fundB, meanval, dum1, dum2;
  double *t_tmp, *mag_tmp, *sig_tmp, *mag_signal, *mag_orig, val0, val1;
  int iter;
  double a, b, c, d, e, min1, min2, fa, fb, fc, p, q, r, s, tol1, xm, tol;
  char *line;
  size_t line_size = MAXLEN;
  
  if((t_tmp = (double *) malloc(N * sizeof(double))) == NULL ||
     (mag_tmp = (double *) malloc(N * sizeof(double))) == NULL ||
     (sig_tmp = (double *) malloc(N * sizeof(double))) == NULL ||
     (mag_signal = (double *) malloc(N * sizeof(double))) == NULL ||
     (mag_orig = (double *) malloc(N * sizeof(double))) == NULL ||
     (subharmA = (double *) malloc((Nsubharm + 1) * sizeof(double))) == NULL ||
     (subharmB = (double *) malloc((Nsubharm + 1) * sizeof(double))) == NULL ||
     (harmA = (double *) malloc((Nharm + 1) * sizeof(double))) == NULL ||
     (harmB = (double *) malloc((Nharm + 1) * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  /* First Get the signal */
  /* Read the signal in from the file if we're doing that */
  if(harm_specsigflag)
    {
      line = malloc(line_size);
      rewind(signalfile);
      val0 = 0.;
      for(i=0;i<N;i++)
	{
	  gnu_getline(&line,&line_size,signalfile);
	  sscanf(line,"%lf %lf %lf",&dum1,&dum2,&mag_signal[i]);
	  val0 += mag[i];
	}
      val0 = val0 / (double) N;
      for(i=0;i<N;i++)
	{
	  mag_signal[i] -= val0;
	  mag_orig[i] = mag[i] - mag_signal[i];
	}

      fclose(signalfile);
      free(line);
    }
  /* Otherwise fit a fourier series at the period to compute the signal */
  else
    {
      for(i=0;i<N;i++)
	{
	  t_tmp[i] = t[i]; mag_tmp[i] = mag[i]; sig_tmp[i] = sig[i];
	}
      dokillharms(N, t_tmp, mag_tmp, sig_tmp, 1, &period, Nsubharm, Nharm, &subharmA, &subharmB, &harmA, &harmB, &fundA, &fundB, &meanval, 0, NULL, amp, 0, KILLHARM_OUTTYPE_DEFAULT, -1.);
      for(i=0;i<N;i++)
	{
	  mag_signal[i] = mag[i] - mag_tmp[i];
	  mag_orig[i] = mag_tmp[i];
	}
    }

  /* Check to see if no signal still lies above the threshhold, or if the original signal lies below the threshhold */
  val0 = getlsampthresh_func(0., N, t, mag_orig, sig, mag_signal, mag_tmp, period, Nsubharm, Nharm, subharmA, subharmB, harmA, harmB, fundA, fundB, minPer, thresh, use_orig_ls);
  val1 = getlsampthresh_func(1., N, t, mag_orig, sig, mag_signal, mag_tmp, period, Nsubharm, Nharm, subharmA, subharmB, harmA, harmB, fundA, fundB, minPer, thresh, use_orig_ls);
  if(val0 < 0)
    {
      *ampthresh_scale = 0.;
      *amp = (*amp)*(*ampthresh_scale);
      free(t_tmp);
      free(mag_tmp);
      free(sig_tmp);
      free(mag_signal);
      free(mag_orig);
      free(subharmA);
      free(subharmB);
      free(harmA);
      free(harmB);
      return;
    }
  if(val1 > 0)
    {
      *ampthresh_scale = -1.;
      *amp = (*amp)*(*ampthresh_scale);
      free(t_tmp);
      free(mag_tmp);
      free(sig_tmp);
      free(mag_signal);
      free(mag_orig);
      free(subharmA);
      free(subharmB);
      free(harmA);
      free(harmB);
      return;
    }
  a = 0.;
  b = 1.;
  tol = 1.0e-5;

  /* Now start the zbrent loop */
  fa = val0; fb = val1;
  fc = fb;
  for (iter=1; iter<=ITMAX; iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c = a;
      fc = fa;
      e=d=b-a;
    }
    if(fabs(fc) < fabs(fb)) {
      a = b; b = c; c = a; fa = fb; fb = fc; fc = fa;
    }
    tol1 = 2.0*EPS*fabs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0)
      {
	*ampthresh_scale = b;
	*amp = (*amp)*(*ampthresh_scale);
	free(t_tmp);
	free(mag_tmp);
	free(sig_tmp);
	free(mag_signal);
	free(mag_orig);
	free(subharmA);
	free(subharmB);
	free(harmA);
	free(harmB);
	return;
      }
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s = fb/fa;
      if (a == c) {
	p = 2.0*xm*s;
	q = 1.0 - s;
      } else {
	q = fa/fc;
	r = fb/fc;
	p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q = (q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p = fabs(p);
      min1 = 3.0*xm*q - fabs(tol1*q);
      min2 = fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
	e=d;
	d=p/q;
      } else {
	d=xm;
	e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1,xm);
    fb = getlsampthresh_func(b, N, t, mag_orig, sig, mag_signal, mag_tmp, period, Nsubharm, Nharm, subharmA, subharmB, harmA, harmB, fundA, fundB, minPer, thresh, use_orig_ls);
  }
  /* Too many iterations - set ampthresh_scale to -2. */
  *ampthresh_scale = -2.;
  *amp = (*amp)*(*ampthresh_scale);
  free(t_tmp);
  free(mag_tmp);
  free(sig_tmp);
  free(mag_signal);
  free(mag_orig);
  free(subharmA);
  free(subharmB);
  free(harmA);
  free(harmB);
  return;
}


double getlsampthresh_func(double amp_scale, int N, double *t, double *mag_orig, double *sig, double *mag_signal, double *mag_tmp, double period, int Nsubharm, int Nharm, double *subharmA, double *subharmB, double *harmA, double *harmB, double fundA, double fundB, double minPer, double thresh, int use_orig_ls)
{
  /* Add the scaled signal to the light curve without the signal, then call ls_oneperiod */
  int i, ngood = 0;
  double val, T, effm, expy, outprob, zval;
  for(i=0;i<N;i++)
    {
      mag_tmp[i] = mag_orig[i] + amp_scale*mag_signal[i];
      if(!isnan(mag_tmp[i]) && !isnan(sig[i]) && sig[i] > 0)
	ngood++;
    }
  if(use_orig_ls) {
    val = ls_oneperiod(t, mag_tmp, N, period);
  } else {
    val = gls_oneperiod(t, mag_tmp, sig, N, period);
  }
  T = t[N-1] - t[0];
  effm = log10(2.0 * T / minPer);
  if(use_orig_ls) {
    expy = -val*0.434294482;
    outprob = effm + expy;
    if(outprob > -2.0) {
      expy = exp(-val);
      outprob = 1.0-pow(1.0-expy,2.0 * T/minPer);
      outprob = log10(outprob);
    }
  } else {
    if(val == 0.) outprob = 0.0;
    else {
      zval = (ngood - 3)*(val/(1.0 - val))/2.0;
      expy = -(ngood - 3)*log10((1.0 + (2.0*zval)/(ngood - 3)))/2.0;
      outprob = effm + expy;
      if(outprob > -2.0) {
	expy = pow(((1.0 + (2.0*zval)/(ngood - 3))),(-(ngood-3)/2.0));
	outprob = 1.0 - pow(1.0-expy,2.0 * T/minPer);
	outprob = log10(outprob);
      }
    }
  }
  return(outprob - thresh);
}









