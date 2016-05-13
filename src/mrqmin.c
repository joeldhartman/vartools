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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "commands.h"
#include "programdata.h"
#include "functions.h"

void mrqmin(double *x, double *y, double *sig, int ndata, double *a, int *ia, int ma, double **covar, double **alpha, double *chisq, double *alamda, void (*funcs)(double *, double *, double *, double **, int, int, int, double *, double **, double *, double *, int *, void *), int Nlin_coeff, double **Design_Matrix, double *lin_coeffs, int *varylin_coeffs,
#ifdef PARALLEL
	    int mfit, double *ochisq, double *atry, double *beta, double *da, double **oneda,
#endif
	    void *userparams)
{
  /*This function is taken from Press et al., 1992 -
    Levenberg-Marquardt method, attempting to reduce the value of chi^2 of a fit between a set of data points x[0....ndata-1], y[0.....ndata-1] with individual standard deviations sig[0.....ndata-1], and a nonlinear function dependent on ma coefficients a[0...ma-1]. The input array ia[0....ma-1] indicates by nonzero entries those components of a that should be fitted for, and by zero entries, those components that should be held fixed at their input values. The program returns current best-fit values for the parameters a[0...ma-1], and chi^2. The arrays covar[0...ma-1][0...ma-1], alpha[0....ma-1][0....ma-1] are used as working space during most iterations. Supply a routine funcs(x,a,yfit,dyda,ma) that evaluates the fitting fuction yfit, and its derivatives dyda[0....ma-1] with respect to the fitting parameters a at x. On the first call provide an initial guess for the parameters a, and set alamda<0 for initialization (which then sets alamda=0.001). If a step succeeds chisq becomes smaller and alamda decreases by a factor of 10. If a step fails alamda grows by a factor of 10. You must call this routine repeatedly until convergence is achieved. Then make one final call with alamda=0, so that covar[0....ma-1][0....ma-1] returns the covariance matrix, and alpha the curvature matrix. (Parameters help fixed will return zero covariances). */
  void covsrt(double **covar, int ma, int *ia, int mfit);
  int gaussj(double **a, int n, double **b, int m);
  void mrqcof(double *x, double *y, double *sig, int ndata, double *a, int *ia, int ma, double **alpha, double *beta, double *chisq, void (*funcs)(double *, double *, double *, double **, int, int, int, double *, double **, double *, double *, int *, void *), int, double **, double *, int *, void *);
  int j, k, l;
#ifndef PARALLEL
  static int mfit;
  static double ochisq[1], *atry, *beta, *da, **oneda;
#endif

  if(*alamda < 0.0) {
#ifndef PARALLEL
    if((atry = (double *) malloc(ma * sizeof(double))) == NULL ||
       (beta = (double *) malloc(ma * sizeof(double))) == NULL ||
       (da = (double *) malloc(ma * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
    for(mfit=0,j=0;j<ma;j++)
      if(ia[j]) mfit++;
    if((oneda = (double **) malloc(mfit * sizeof(double *))) == NULL)
      error(ERR_MEMALLOC);
    for(j=0;j<mfit;j++)
      if((oneda[j] = (double *) malloc(sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
#endif
    *alamda = 0.001;
    mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs, Nlin_coeff, Design_Matrix, lin_coeffs, varylin_coeffs, userparams);
    ochisq[0]=(*chisq);
    for(j=0;j<ma;j++)
      atry[j] = a[j];
  }
  for(j=0;j<mfit;j++) {
    for(k=0;k<mfit;k++) covar[j][k] = alpha[j][k];
    covar[j][j]=alpha[j][j]*(1.0+(*alamda));
    oneda[j][0] = beta[j];
  }
  if(gaussj(covar,mfit,oneda,1))
    {
      *alamda = 0;
    }
  for(j=0;j<mfit;j++) da[j] = oneda[j][0];
  if (*alamda == 0.0) {
    covsrt(covar,ma,ia,mfit);
    covsrt(alpha,ma,ia,mfit);
#ifndef PARALLEL
    for(j=0;j<mfit;j++)
      free(oneda[j]);
    free(oneda);
    free(da);
    free(beta);
    free(atry);
#endif
    return;
  }
  for(j=0,l=0;l<ma;l++)
    if(ia[l]) atry[l]=a[l]+da[j++];
  mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs,Nlin_coeff,Design_Matrix,lin_coeffs,varylin_coeffs, userparams);
  if (*chisq < ochisq[0]) {
    *alamda *= 0.1;
    ochisq[0] = (*chisq);
    for (j=0;j<mfit;j++) {
      for(k=0;k<mfit;k++) alpha[j][k] = covar[j][k];
      beta[j] = da[j];
    }
    for(l=0;l<ma;l++) a[l]=atry[l];
  } else {
    *alamda *= 10.0;
    *chisq=ochisq[0];
  }
}

void mrqcof(double *x, double *y, double *sig, int ndata, double *a, int *ia, int ma, double **alpha, double *beta, double *chisq, void (*funcs)(double *, double *, double *, double **, int, int,int,double *, double **, double *, double *, int *, void *), int Nlin_coeff, double **Design_Matrix, double *lin_coeffs, int *varylin_coeffs, void *userparams)
{
  int i, j, k, l, m, mfit = 0;
  double *ymod, wt, sig2i, dy, **dyda;
  if((dyda = (double **) malloc(ndata * sizeof(double *))) == NULL ||
     (ymod = (double *) malloc(ndata * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  for(j=0;j<ndata;j++)
    if((dyda[j] = (double *) malloc(ma * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
  for(j=0;j<ma;j++)
    if(ia[j]) mfit++;
  for(j=0;j<mfit;j++) {
    for(k=0;k<j;k++) alpha[j][k] = 0.;
    beta[j] = 0.;
  }
  *chisq = 0.;

  /* Evaluate the function */
  (*funcs)(x,a,ymod,dyda,ma,ndata,Nlin_coeff,lin_coeffs,Design_Matrix,y,sig,varylin_coeffs, userparams);

  for(i=0;i<ndata;i++)
    {
      if(!isnan(y[i]) && sig[i] > 0.)
	{
	  sig2i = 1.0/(sig[i]*sig[i]);
	  dy = y[i] - ymod[i];
	  for (j=0,l=0;l<ma;l++)
	    {
	      if(ia[l]) {
		wt=dyda[i][l]*sig2i;
		for(j++,k=0,m=0;m<=l;m++)
		  if(ia[m]) alpha[j-1][k++] += wt*dyda[i][m];
		beta[j-1] += dy*wt;
	      }
	    }
	  *chisq += dy*dy*sig2i;
	}
    }
  for(j=1;j<mfit;j++)
    for(k=0;k<j;k++)
      alpha[k][j]=alpha[j][k];
  free(ymod);
  for(j=0;j<ndata;j++)
    free(dyda[j]);
  free(dyda);
}

#ifdef SWAP
#undef SWAP
#endif

#define SWAP(a,b) swap=(a);(a)=(b);(b)=swap

void covsrt(double **covar, int ma, int *ia, int mfit)
{
  /* Expand in storage the covariance matrix covar, so as to take into account parameters that are being held fixed. (For the latter, return zero covariances.) */
  int i, j, k;
  double swap;
  for (i=mfit;i<ma;i++)
    for(j=0;j<i;j++)
      covar[i][j]=covar[j][i]=0.0;
  k=mfit-1;
  for(j=ma-1;j>=0;j--)
    {
      if(ia[j]) {
	for(i=0;i<ma;i++) { SWAP(covar[i][k],covar[i][j]);}
	for(i=0;i<ma;i++) { SWAP(covar[k][i],covar[j][i]);}
	k--;
      }
    }
}

#undef SWAP
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp

int gaussj(double **a, int n, double **b, int m)
{
  /* Linear equation solution by Gauss-Jordan elimination. a[0...n-1][0....n-1] is the input matrix. b[0....n-1][0...m-1] is input containing the m right-hand side vectors. On output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of solution vectors.

This algorithm is taken from Press et al. 1992, it has been modified slightly to put it into the program vartools by J. Hartman

*/
  int *indxc, *indxr, *ipiv;
  int i, icol, irow, j, k, l, ll;
  double big, dum, pivinv, temp;

  if((indxc = (int *) malloc(n * sizeof(int))) == NULL ||
     (indxr = (int *) malloc(n * sizeof(int))) == NULL ||
     (ipiv = (int *) malloc(n * sizeof(int))) == NULL)
    error(ERR_MEMALLOC);

  for(j=0;j<n;j++) ipiv[j] = 0;
  for(i=0;i<n;i++) {
    big = 0.0;
    for(j=0;j<n;j++)
      if(ipiv[j] != 1)
	for(k=0;k<n;k++) {
	  if(ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big = fabs(a[j][k]);
	      irow = j;
	      icol = k;
	    }
	  }
	}
    ++(ipiv[icol]);
    if(irow != icol) {
      for(l=0;l<n;l++) { SWAP(a[irow][l],a[icol][l]); }
      for(l=0;l<m;l++) { SWAP(b[irow][l],b[icol][l]); }
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if(a[icol][icol] == 0.0)
      {
	free(ipiv);
	free(indxr);
	free(indxc);
	return(1);
      }
    pivinv=1.0/a[icol][icol];
    a[icol][icol] = 1.0;
    for(l=0;l<n;l++) a[icol][l] *= pivinv;
    for(l=0;l<m;l++) b[icol][l] *= pivinv;
    for(ll=0;ll<n;ll++)
      if(ll != icol) {
	dum=a[ll][icol];
	a[ll][icol] = 0.0;
	for(l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
	for(l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }

  for(l=n-1;l>=0;l--) {
    if(indxr[l] != indxc[l])
      for(k=0;k<n;k++){ SWAP(a[k][indxr[l]],a[k][indxc[l]]); }
  }
  free(ipiv);
  free(indxr);
  free(indxc);
  return(0);
}
