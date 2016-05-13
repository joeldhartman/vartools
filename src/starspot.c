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

/* This file contains functions for fitting a starspot model to a light curve with a known period

The starspot model comes from Dorren, 1987, ApJ, 320, 756.
*/

#define TWOPI 6.28318530717958647692528676656
#define PIOVERTWO 1.570796327
#define _PI_ 3.14159265358979323846264338328
#define PIOVERTHREE 1.047197551
#define TWOTHIRDSPI 2.094395102
#define LNTEN 2.302585093

double chisqstarspot(double *a_, int ma, int N, double *t, double *mag, double *sig, void *userparams)
{
  void starspot_(int N, double *t, double *dmag, double a, double b, double alpha, double i, double chi, double psi0, double P);
  int j;
  double chisq, *newmag, temp;
  chisq = 0.;
  if((newmag = (double *) malloc(N * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  starspot_(N, t, newmag, a_[1], a_[2], a_[3], a_[4], a_[5], a_[6], a_[0]);
  for(j=0;j<N;j++)
    if(!isnan(mag[j]))
      {
	temp = (mag[j] - newmag[j] - a_[7]);
	chisq += temp*temp/(sig[j]*sig[j]);
      }
  free(newmag);
  return(chisq);
}

/*
#ifdef TINY
#undef TINY
#endif
#define TINY 1.0e-10
#define NMAX 5000

#ifdef SWAP
#undef SWAP
#endif

#define SWAP(a,b) swap=(a);(a)=(b);(b)=swap

int amoeba(double **p, double *y, int *ia, int ndim, double ftol, double (*funk)(double *, int, int, double *, double *, double *), int *nfunk, int N, double *t, double *mag, double *sig)
{
  double amotry(double **p, double *y, int *ia, double *psum, int ndim, int ndimused, double (*funk)(double *, int, int, double *, double *, double *), int ihi, double fac, int N, double *t, double *mag, double *sig);
  int i, ihi, ilo, inhi, j, mpts, ndimused;
  double rtol, sum, swap, ysave, ytry, *psum;

  ndimused = 0;
  for(j=0;j<ndim;j++)
    if(ia[j])
      ndimused++;

  if(!ndimused)
    return 1;

  mpts = ndimused+1;

  if((psum = (double *) malloc(ndim * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  *nfunk = 0;
  for (j=0;j<ndim;j++) {
    for(sum=0.0,i=0;i<mpts;i++) sum += p[i][j];
    psum[j]=sum;
  }
  for(;;) {
    ilo = 1;
    ihi = y[0] > y[1] ? (inhi=1,0) : (inhi=0,1);
    for (i=0;i<mpts;i++) {
      if(y[i] <= y[ilo]) ilo = i;
      if(y[i] > y[ihi]) {
	inhi = ihi;
	ihi = i;
      } else if(y[i] > y[inhi] && i != ihi) inhi = i;
    }
    rtol = 2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
    if(rtol < ftol) {
      SWAP(y[0],y[ilo]);
      for(i=0;i<ndim;i++) {SWAP(p[0][i],p[ilo][i]);}
      break;
    }
    if (*nfunk >= NMAX)
      return 1;
    //error(ERR_TOOMANYAMOEBAITERATIONS);

    *nfunk += 2;
    ytry = amotry(p,y,ia,psum,ndim,ndimused,funk,ihi,-1.0, N, t, mag, sig);
    if(ytry <= y[ilo])
      ytry = amotry(p,y,ia,psum,ndim,ndimused,funk,ihi,2.0,N,t,mag,sig);
    else if(ytry >= y[inhi])
      {
	ysave=y[ihi];
	ytry = amotry(p,y,ia,psum,ndim,ndimused,funk,ihi,0.5,N,t,mag,sig);
	if (ytry >= ysave) {
	  for (i = 0; i<mpts;i++) {
	    if (i != ilo) {
	      for(j=0;j<ndim;j++)
		{
		  if(ia[j])
		    {
		      p[i][j] = psum[j] = 0.5*(p[i][j] + p[ilo][j]);
		    }
		  else
		    {
		      psum[j] = p[i][j];
		    }
		}
	      y[i]=(*funk)(psum,ndim,N,t,mag,sig);
	    }
	  }
	  *nfunk += ndim;
	  for (j=0;j<ndim;j++) {
	    for(sum=0.0,i=0;i<mpts;i++) sum += p[i][j];
	    psum[j]=sum;
	  }
	}
      } else --(*nfunk);
  }
  free(psum);
  return 0;
}

double amotry(double **p, double *y, int *ia, double *psum, int ndim, int ndimused, double (*funk)(double *, int, int, double *, double *, double *), int ihi, double fac, int N, double *t, double *mag, double *sig)
{
  int j;
  double fac1, fac2, ytry, *ptry;
  if((ptry = (double *) malloc(ndim * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  fac1 = (1.0 - fac)/ndimused;
  fac2 = fac1 - fac;
  for(j=0;j<ndim;j++)
    {
      if(ia[j])
	ptry[j]=psum[j]*fac1 - p[ihi][j]*fac2;
      else
	ptry[j]=p[ihi][j];
    }
  ytry=(*funk)(ptry,ndim,N,t,mag,sig);
  if (ytry < y[ihi]) {
    y[ihi]=ytry;
    for (j = 0; j< ndim;j++) {
      if(ia[j])
	{
	  psum[j] += ptry[j] - p[ihi][j];
	  p[ihi][j] = ptry[j];
	}
    }
  }
  free(ptry);
  return ytry;
}
*/

void starspot(double *t, double *a_, double *yfit, double **dyda, int ma, int N, void *userparams)
{
/*
Compute a Dorren starspot model for a given set of parameters:

The function is in the format used by the mrqmin non-linear minimization routine from Press et al. 1992.

The input parameters are stored in the vector a_:
They are:
a_[0] = P, a_[1] = a, a_[2] = b, a_[3] = alpha, a_[4] = i, a_[5] = chi, a_[6] = psi0, a_[7] = mconst
see Dorren, 1987 for description
t is the input times of observations
mconst is the offset of the light curve.

dyda is a [0....n-1][0....ma-1] matrix giving the derivative of the light curve w.r.t. the parameters at each input time.

yfit is the output light curve
ma is the number of parameters
N is the number of data points
*/
  int j;
  double psi, T, delta, zeta, beta, cbeta, sbeta, calpha, salpha, s2alpha, ci, si, cchi, schi, cpsi, spsi, sdelta, cdelta, szeta, czeta, A, B, c3alpha, P, a, b, alpha, i, chi, psi0, mconst, dAdbeta, dBdbeta, dzetadbeta, dTdbeta, ddeltadbeta, dbetadi, dbetadchi, dbetadpsi, denom, dAdalpha, dBdalpha, dzetadalpha, dTdalpha, ddeltadalpha, cT;

  if(ma != 8)
    error(ERR_CODEERROR);
  P = a_[0];
  a = a_[1];
  b = a_[2];
  alpha = a_[3];
  i = a_[4];
  chi = a_[5];
  psi0 = a_[6];
  mconst = a_[7];
  ci = cos(i);
  schi = sin(chi);
  si = sin(i);
  cchi = cos(chi);
  salpha = sin(alpha);
  s2alpha = salpha*salpha;
  calpha = cos(alpha);
  c3alpha = calpha*calpha*calpha;
  for(j=0;j<N;j++)
    {
      psi = TWOPI * ((t[j] - t[0]) / P) + psi0;
      cpsi = cos(psi);
      cbeta = (ci*schi + si*cchi*cpsi);

      /*Check which case we're in */
      beta = acos(cbeta);
      if(beta - alpha >= PIOVERTWO)
	{
	  A = 0.;
	  B = 0.;
	}
      else if(beta + alpha <= PIOVERTWO)
	{
	  A = _PI_ * cbeta * s2alpha;
	  sbeta = sin(beta);
	  B = PIOVERTHREE*(-2.0*c3alpha - 3.0*sbeta*sbeta*calpha*s2alpha) + TWOTHIRDSPI;
	}
      else
	{
	  sbeta = sin(beta);
	  cdelta = (calpha * cbeta) / (salpha * sbeta);
	  delta = acos(cdelta);
	  sdelta = sin(delta);
	  szeta = sdelta * salpha;
	  czeta = calpha / sbeta;
	  if(beta <= PIOVERTWO)
	    T = atan(szeta*sbeta/cbeta);
	  else
	    T = _PI_ - atan(-szeta*sbeta/cbeta);

	  zeta = atan2(szeta,czeta);
	  A = zeta + (_PI_ - delta)*cbeta*s2alpha - szeta*sbeta*calpha;
	  B = (PIOVERTHREE - (delta/3.0))*(-2.0*c3alpha - 3.0*sbeta*sbeta*calpha*s2alpha) + (TWOTHIRDSPI - (2.0*T/3.0)) + szeta*2.0*sbeta*cbeta*(2.0 - 3.0*calpha*calpha)/6.0;
	}

      yfit[j] = -2.5*log10(1.0 - a*A - b*B) + mconst;

      /* Evaluate the derivatives */
      //mconst
      dyda[j][7] = 1.;

      denom = LNTEN * (1.0 - a*A - b*B);
      //a
      dyda[j][1] = 2.5*A / denom;
      //b
      dyda[j][2] = 2.5*B / denom;
      /* Now for the more complicated derivatives */
      if(beta - alpha >= PIOVERTWO)
	{
	  //P
	  dyda[j][0] = 0.;
	  //alpha
	  dyda[j][3] = 0.;
	  //i
	  dyda[j][4] = 0.;
	  //chi
	  dyda[j][5] = 0.;
	  //psi0
	  dyda[j][6] = 0.;
	}
      else if(beta + alpha <= PIOVERTWO)
	{
	  dAdbeta = - _PI_ * sbeta * s2alpha;
	  dBdbeta = - TWOPI * sbeta * cbeta *calpha * s2alpha;
	  spsi = sin(psi);
	  dbetadpsi = si * cchi * spsi / sbeta;
	  dbetadi = (si*schi - ci*cchi*cpsi)/sbeta;
	  dbetadchi = (-ci*cchi + si*schi*cpsi)/sbeta;
	  //P
	  dyda[j][0] = 2.5*(a*dAdbeta + b*dBdbeta)*(dbetadpsi*(TWOPI*(t[0] - t[j])/(P*P)))/denom;
	  //i
	  dyda[j][4] = 2.5*(a*dAdbeta + b*dBdbeta)*dbetadi/denom;
	  //chi
	  dyda[j][5] = 2.5*(a*dAdbeta + b*dBdbeta)*dbetadchi/denom;
	  //psi0
	  dyda[j][6] = 2.5*(a*dAdbeta + b*dBdbeta)*dbetadpsi/denom;

	  //alpha
	  dAdalpha = TWOPI * cbeta * salpha * calpha;
	  dBdalpha = TWOPI * cbeta * cbeta * calpha * calpha * salpha + _PI_ * sbeta * sbeta * salpha * s2alpha;
	  dyda[j][3] = 2.5*(a*dAdalpha + b*dBdalpha)/denom;
	}
      else
	{
	  ddeltadbeta = calpha /(salpha * sdelta * sbeta * sbeta);
	  dzetadbeta = calpha * cbeta / (sbeta * sbeta * szeta);
	  if(beta <= PIOVERTWO)
	    {
	      cT = cos(T);
	      dTdbeta = cT*cT*(czeta*sbeta*dzetadbeta/(cbeta) + szeta/(cbeta*cbeta));
	    }
	  else
	    {
	      cT = cos(_PI_ - T);
	      dTdbeta = cT*cT*(czeta*sbeta*dzetadbeta/(cbeta) + szeta/(cbeta*cbeta));
	    }
	  dBdbeta = -ddeltadbeta * (-2.0*c3alpha - 3.0*sbeta*sbeta*calpha*s2alpha)/3.0 + (PIOVERTHREE - (delta/3.0))*(-6.0*sbeta*cbeta*calpha*s2alpha) - 2.0*dTdbeta/3.0 + czeta*sbeta*cbeta*(2.0-3.0*calpha*calpha)*dzetadbeta/3.0 + szeta*(cbeta*cbeta - sbeta*sbeta)*(2.0-3.0*calpha*calpha)/3.0;
	  dAdbeta = -(_PI_ - delta)*sbeta*s2alpha - czeta*sbeta*calpha*dzetadbeta - szeta*cbeta*calpha + dzetadbeta - ddeltadbeta*cbeta*s2alpha;

	  spsi = sin(psi);
	  dbetadpsi = si * cchi * spsi / sbeta;
	  dbetadi = (si*schi - ci*cchi*cpsi)/sbeta;
	  dbetadchi = (-ci*cchi + si*schi*cpsi)/sbeta;
	  //P
	  dyda[j][0] = 2.5*(a*dAdbeta + b*dBdbeta)*(dbetadpsi*(TWOPI*(t[0] - t[j])/(P*P)))/denom;
	  //i
	  dyda[j][4] = 2.5*(a*dAdbeta + b*dBdbeta)*dbetadi/denom;
	  //chi
	  dyda[j][5] = 2.5*(a*dAdbeta + b*dBdbeta)*dbetadchi/denom;
	  //psi0
	  dyda[j][6] = 2.5*(a*dAdbeta + b*dBdbeta)*dbetadpsi/denom;

	  //alpha
	  dzetadalpha = salpha/(sbeta*szeta);
	  ddeltadalpha = cbeta/(sdelta*sbeta*s2alpha);
	  if(beta <= PIOVERTWO)
	    {
	      cT = cos(T);
	      dTdalpha = cT*cT*czeta*sbeta*dzetadalpha/cbeta;
	    }
	  else
	    {
	      cT = cos(_PI_ - T);
	      dTdalpha = cT*cT*czeta*sbeta*dzetadalpha/cbeta;
	    }
	  dAdalpha = dzetadalpha - ddeltadalpha*cbeta*s2alpha + (_PI_ - delta)*2.0*cbeta*salpha*calpha - czeta*sbeta*calpha*dzetadalpha + szeta*sbeta*salpha;
	  dBdalpha = -ddeltadalpha*(-2.0*c3alpha - 3.0*sbeta*sbeta*calpha*s2alpha)/3.0 + (_PI_ - delta)*(2.0*calpha*calpha*salpha + sbeta*sbeta*s2alpha*salpha - 2.0*sbeta*sbeta*calpha*calpha*salpha) - 2.0*dTdalpha/3.0 + czeta*sbeta*cbeta*(2.0-3.0*calpha*calpha)*dzetadalpha/3.0 + 2.0*szeta*sbeta*cbeta*calpha*salpha;
	  dyda[j][3] = 2.5*(a*dAdalpha + b*dBdalpha)/denom;
	}

    }
}



void starspot_(int N, double *t, double *dmag, double a, double b, double alpha, double i, double chi, double psi0, double P)
{
  int j;
  double psi, T, delta, zeta, beta, cbeta, sbeta, calpha, salpha, s2alpha, ci, si, cchi, schi, cpsi, sdelta, cdelta, szeta, czeta, A, B, c3alpha;
  alpha = alpha * _PI_ / 180.;
  i = i * _PI_ / 180.;
  chi = chi * _PI_ / 180.;
  psi0 = psi0 * _PI_ / 180.;
  ci = cos(i);
  schi = sin(chi);
  si = sin(i);
  cchi = cos(chi);
  salpha = sin(alpha);
  s2alpha = salpha*salpha;
  calpha = cos(alpha);
  c3alpha = calpha*calpha*calpha;
  for(j=0;j<N;j++)
    {
      psi = TWOPI * ((t[j] - t[0]) / P) + psi0;
      cpsi = cos(psi);
      cbeta = (ci*schi + si*cchi*cpsi);

      /*Check which case we're in */
      beta = acos(cbeta);
      if(beta - alpha >= PIOVERTWO)
	{
	  A = 0.;
	  B = 0.;
	}
      else if(beta + alpha <= PIOVERTWO)
	{
	  A = _PI_ * cbeta * s2alpha;
	  sbeta = sin(beta);
	  B = PIOVERTHREE*(-2.0*c3alpha - 3.0*sbeta*sbeta*calpha*s2alpha) + TWOTHIRDSPI;
	}
      else
	{
	  sbeta = sin(beta);
	  cdelta = (calpha * cbeta) / (salpha * sbeta);
	  delta = acos(cdelta);
	  sdelta = sin(delta);
	  szeta = sdelta * salpha;
	  czeta = calpha / sbeta;
	  if(beta <= PIOVERTWO)
	    T = atan(szeta*sbeta/cbeta);
	  else
	    T = _PI_ - atan(-szeta*sbeta/cbeta);

	  zeta = atan2(szeta,czeta);
	  A = zeta + (_PI_ - delta)*cbeta*s2alpha - szeta*sbeta*calpha;
	  B = (PIOVERTHREE - (delta/3.0))*(-2.0*c3alpha - 3.0*sbeta*sbeta*calpha*s2alpha) + (TWOTHIRDSPI - (2.0*T/3.0)) + szeta*2.0*sbeta*cbeta*(2.0 - 3.0*calpha*calpha)/6.0;
	}
      dmag[j] = -2.5*log10(1.0 - a*A - b*B);
    }
}

#define CONVERGENCELIMIT 0.0001
#define INITIALSTEP 0.1

void fitstarspot_amoeba(int N, double *t, double *mag, double *sig, double *P, double *a, double *b, double *alpha, double *i, double *chi, double *psi0, double *mconst, int fitP, int fita, int fitb, int fitalpha, int fiti, int fitchi, int fitpsi0, int fitmconst, double *chi2_, int correctlc, int omodel, char *modelname)
{
  /* This is a function that fits a single starspot model to a light curve with t, mag, and sig. The parameters are:
P = period, a, b such that:
delta_mag = -2.5*log(1 - a*A - b*B) and A and B are as in Dorren, 1987.
i = inclination in degrees,
chi = spot latitude in degrees,
alpha = spot angular radius in degrees,
psi0 = spot initial longitude in degrees,
mconst = magnitude of star without spot.
fitP, fita etc are flags denoting whether or not a parameter is to be fit for.
eP, ea... etc are the output uncertainties on each of the parameters.
chi2_ will be the output chi2_ per dof after fitting the light curve.
correctlc is a flag denoting whether or not to subtract off the spot model from the light curve */

  int l, j, k, ma, ngood, nfunk, nvar;
  double **p, *y, *delmag, ftol, meanval1, meanval2;
  int *ia, amoeba_val;
  double (*funcs)(double *, int , int,  double *, double *, double *, void *);
  FILE *outfile;
  funcs = &chisqstarspot;

  ma = 8;

  if((*P <= 0))
    {
      *chi2_ = 9999999.;
      return;
    }

  if((ia = (int *) malloc(ma * sizeof(int))) == NULL)
    error(ERR_MEMALLOC);

  if(fitP) ia[0] = 1; else ia[0] = 0;
  if(fita) ia[1] = 1; else ia[1] = 0;
  if(fitb) ia[2] = 1; else ia[2] = 0;
  if(fitalpha) ia[3] = 1; else ia[3] = 0;
  if(fiti) ia[4] = 1; else ia[4] = 0;
  if(fitchi) ia[5] = 1; else ia[5] = 0;
  if(fitpsi0) ia[6] = 1; else ia[6] = 0;
  if(fitmconst) ia[7] = 1; else ia[7] = 0;
  for(nvar=0,j=0;j<ma;j++)
    if(ia[j])
      nvar++;

  if((p = (double **) malloc((nvar + 1)*sizeof(double *))) == NULL ||
     (y = (double *) malloc((nvar + 1)*sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  for(j=0;j<nvar+1;j++)
    if((p[j] = (double *) malloc(ma * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);

  /* estimate the minimum mag if mconst is less than zero */
  if(*mconst < 0)
    {
      j = 0;
      while(isnan(mag[j]))
	j++;
      *mconst = mag[j];
      for(j++;j<N;j++)
	if(mag[j] < *mconst)
	  *mconst = mag[j];
    }

  if(correctlc || omodel)
    {
      for(j=0,meanval1=0.0,meanval2=0.0;j<N;j++)
	if(!isnan(mag[j]))
	  {
	    meanval1 += mag[j]/(sig[j]*sig[j]);
	    meanval2 += 1./(sig[j]*sig[j]);
	  }
      meanval1 /= meanval2;
    }

  /* Initialize the trial guesses */
  p[0][0] = *P; p[0][1] = *a; p[0][2] = *b; p[0][3] = *alpha; p[0][4] = *i; p[0][5] = *chi; p[0][6] = *psi0; p[0][7] = *mconst;


  for(k=1;k<=nvar;k++)
    for(l=0,j=0;j<ma;j++)
      {
	if(ia[j])
	  {
	    if(l == k-1)
	      p[k][j] = p[0][j]*(1 + INITIALSTEP);
	    else
	      p[k][j] = p[0][j];
	    l++;
	  }
	else
	  p[k][j] = p[0][j];
      }
  for(k=0;k<nvar+1;k++)
    y[k] = (*funcs)(p[k],ma,N,t,mag,sig, NULL);

  ftol = CONVERGENCELIMIT;

  /* Run amoeba */
  amoeba_val = amoeba(p, y, ia, ma, ftol, funcs, &nfunk, 0, N, t, mag, sig, NULL);

  /* If amoeba didn't converge, then write out a garbage model */
  if(amoeba_val && nvar > 0)
    {
      *P = -1.; *a = -1.; *b = -1; *alpha = -1.; *i = -1.; *chi = -1.; *psi0 = -1.; *mconst = -1.;
    }
  else if(nvar > 0)
    {
      /* Find the minimum among the vertices */
      k = 0;
      for(j=0;j<=nvar;j++)
	if(y[j] < y[k])
	  k=j;

      *P = p[k][0]; *a = p[k][1]; *b = p[k][2]; *alpha = p[k][3]; *i = p[k][4]; *chi = p[k][5]; *psi0 = p[k][6]; *mconst = p[k][7];

    }
  /* Remove the signal if we're doing that */
  if(correctlc || omodel)
    {
      if((delmag = (double *) malloc(N * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      if(!amoeba_val || nvar == 0)
	starspot_(N, t, delmag, *a, *b, *alpha, *i, *chi, *psi0, *P);
      if(omodel)
	{
	  if((outfile = fopen(modelname,"w")) == NULL)
	    error2(ERR_CANNOTWRITE,modelname);
	  if(!amoeba_val || nvar == 0)
	    {
	      for(j=0;j<N;j++)
		if(!isnan(mag[j]))
		  {
		    fprintf(outfile,"%f %f %f %f\n",t[j],mag[j],delmag[j] + *mconst,sig[j]);
		  }
	    }
	  else
	    {
	      for(j=0;j<N;j++)
		if(!isnan(mag[j]))
		  fprintf(outfile,"%f %f 0. %f\n",t[j],mag[j],sig[j]);
	    }
	  fclose(outfile);
	}
      if(correctlc && (!amoeba_val || nvar == 0))
	{
	  for(j=0;j<N;j++)
	    if(!isnan(mag[j]))
	      { mag[j] = mag[j] - delmag[j] - *mconst + meanval1; }
	}
    }

  ngood = 0;
  for(j=0;j<N;j++)
    if(!isnan(mag[j]))
      ngood++;

  if(ngood > 9)
    *chi2_ = (y[k] / (ngood - nvar));
  else
    *chi2_ = -1.;
  if(correctlc || omodel)
    free(delmag);
  free(ia);
  free(y);
  for(j=0;j<=nvar;j++)
    {
      free(p[j]);
    }
  free(p);
}



