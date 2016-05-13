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
/*                                                                           */
/*     The routines in this file were adapted from Numerical Recipes in C,   */
/*     second edition by Press, Teukolsky, Vetterling, and Flannery.         */
#include "commands.h"
#include "programdata.h"
#include "functions.h"

#ifdef TINY
#undef TINY
#endif
#define TINY 1.0e-10
#define VARTOOLS_AMOEBA_DEFAULT_MAXSTEPS 5000

#ifdef SWAP
#undef SWAP
#endif

#define SWAP(a,b) swap=(a);(a)=(b);(b)=swap

int amoeba(double **p, double *y, int *ia, int ndim, double ftol, double (*funk)(double *, int, int, double *, double *, double *, void *), int *nfunk, int maxeval, int N, double *t, double *mag, double *sig, void *userparams)
/* This function runs the amoeba downhill simplex algorithm.
   p - a 2-d array storing the simplex
   y - a vector storing the chi2 values for each row in the simplex
   ia - an integer array storing flags that indicate whether or
        not the parameter is to be varied. The function will allocate memory for
        the ia array, and update it, for each parameter added. This array will
        be passed to amoeba.
   ndim - the number of parameters in the simplex.
   ftol - the fractional in chi2 to decide when to stop the minimization.
   funk - a pointer to a function which provides chi2.
      funk has the expected format:
      double funk(double *param, int ndim, int N, double *t, double *mag,
                  double *sig, void *userparams)
           where param is a set of parameter values to calculate chi2 for,
           and the function should return chi2.
      The other input parameters to funk are as described below.
   nfunk - on return, this stores the number of function calls executed by
           amoeba.
   maxeval - maximum number of function evaluations to perform before giving
             up. If this is <= 0 then the default number is used.
   The following terms are passed directly to funk and not actually used
   elsewhere within amoeba, they can, but do not have to, take on the context
   given below:
    N - number of points in the time series that is being fit.
    t - times of observation
    mag - set of observed magnitude values being fit
    sig - set of magnitude uncertainties
    userparams - structure containing other user-defined parameters.
         must be cast to the appropriate type by funk.
*/
{
  double amotry(double **p, double *y, int *ia, double *psum, int ndim, int ndimused, double (*funk)(double *, int, int, double *, double *, double *, void *), int ihi, double fac, int N, double *t, double *mag, double *sig, void *userparams);
  int i, ihi, ilo, inhi, j, mpts, ndimused;
  double rtol, sum, swap, ysave, ytry, *psum;

  if(maxeval <= 0) maxeval = VARTOOLS_AMOEBA_DEFAULT_MAXSTEPS;

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
    if (*nfunk >= maxeval) {
      free(psum);
      return 1;
    }
    //error(ERR_TOOMANYAMOEBAITERATIONS);

    *nfunk += 2;
    ytry = amotry(p,y,ia,psum,ndim,ndimused,funk,ihi,-1.0, N, t, mag, sig, userparams);
    
    if(!isnan(ytry) && (isnan(y[ilo]) || ytry <= y[ilo]))
      ytry = amotry(p,y,ia,psum,ndim,ndimused,funk,ihi,2.0,N,t,mag,sig, userparams);
    else if(isnan(ytry) || ytry >= y[inhi])
      {
	ysave=y[ihi];
	ytry = amotry(p,y,ia,psum,ndim,ndimused,funk,ihi,0.5,N,t,mag,sig, userparams);
	if (isnan(ytry) || ytry >= ysave) {
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
	      y[i]=(*funk)(psum,ndim,N,t,mag,sig, userparams);
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

double amotry(double **p, double *y, int *ia, double *psum, int ndim, int ndimused, double (*funk)(double *, int, int, double *, double *, double *, void *), int ihi, double fac, int N, double *t, double *mag, double *sig, void *userparams)
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
  ytry=(*funk)(ptry,ndim,N,t,mag,sig,userparams);
  if (!isnan(ytry) && (ytry < y[ihi])) {
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

void amoeba_initializesimplexchi2(int Nparameters, int Ntovary, double **p, double **chi2vals, double (*funk)(double *, int, int, double *, double *, double *, void *), int N, double *t, double *mag, double *err, void * userparam)
/* This function sets up the initial chi2 vector for the amoeba simplex.
   This is called after the simplex has been setup using 
   incrementparameters_foramoeba.
   Nparameters - the number of parameters in the simplex.
   Ntovary - the total number of parameters that are free to vary.
   p - A 2-d array storing the simplex.
   chi2vals - A pointer to a 1-d vector which will store the chi2 values for
     each point on the simplex. The program will allocate memory for this 
     vector.
   funk - a pointer to a function which provides chi2.
      funk has the expected format:
      double funk(double *param, int Nparameters, int N, double *t, double *mag,
                  double *sig, void *userparams)
           where param is a set of parameter values to calculate chi2 for,
           and the function should return chi2.
      The other input parameters to funk are as described below.
   The following terms are passed directly to funk and not actually used
   elsewhere within this function, they can, but do not have to, 
   take on the context given below:
    N - number of points in the time series that is being fit.
    t - times of observation
    mag - set of observed magnitude values being fit
    sig - set of magnitude uncertainties
    userparams - structure containing other user-defined parameters.
         must be cast to the appropriate type by funk.
*/
{
  int i;
  if(((*chi2vals) = (double *) malloc((Ntovary + 1)*sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0; i < Ntovary + 1; i++) {
    (*chi2vals)[i] = (*funk)(p[i], Nparameters, N, t, mag, err, userparam);
  }
}

#define DEFAULT_INITIAL_PVECSIZE 10
void incrementparameters_foramoeba(int *Nparameters, int *Ntovary, double ***p, int **ia, int varyparam, double initialval, double stepval)
/* This function sets up the initial simplex for amoeba.
   This is called once for each parameter that one wishes to include in the
   simplex.
   Nparameters - the number of parameters in the simplex before adding the new
                 parameter. On output this will be the number of parameters
		 including the new parameter. Nparameters should be zero on
		 input for the first parameter added to the simplex.
   Ntovary - the total number of parameters that are free to vary. This is
             incremented by the function if varyparam == 1.
   p - pointer to a 2-d array storing the simplex. The function will allocate
       memory for p, and update it, for each parameter added. This array will
       then be passed to amoeba.
   ia - pointer to an integer array storing flags that will indicate whether or
        not the parameter is to be varied. The function will allocate memory for
        the ia array, and update it, for each parameter added. This array will
        be passed to amoeba.
   varyparam - a flag indicating whether or not the new parameter is to be
        varied. 1 - yes, 0 - no.
   initialval - initial value to adopt for the parameter.
   stepval - initial step-size to use for this parameter in defining the
             simplex.
*/
{
  int i,j,k;
  int N;
  N = *Nparameters;
  if(N % DEFAULT_INITIAL_PVECSIZE == 0)
    {
      if(!N)
	{
	  *p = (double **) malloc((DEFAULT_INITIAL_PVECSIZE + 1) * sizeof(double *));
	  *ia = (int *) malloc((DEFAULT_INITIAL_PVECSIZE) * sizeof(int));
	  for(i=0;i< (DEFAULT_INITIAL_PVECSIZE + 1); i++)
	    (*p)[i] = (double *) malloc((DEFAULT_INITIAL_PVECSIZE) * sizeof(double));
	}
      else
	{
	  *p = (double **) realloc((*p), (N + DEFAULT_INITIAL_PVECSIZE + 1) * sizeof(double *));
	  *ia = (int *) realloc((*ia), (N + DEFAULT_INITIAL_PVECSIZE) * sizeof(int));
	  for(i=0;i< (N + 1); i++)
	    (*p)[i] = (double *) realloc((*p)[i], (N + DEFAULT_INITIAL_PVECSIZE) * sizeof(double));
	  for(i=(N+1); i < (N + DEFAULT_INITIAL_PVECSIZE + 1); i++)
	    (*p)[i] = (double *) malloc((N + DEFAULT_INITIAL_PVECSIZE) * sizeof(double));
	}
    }
  if(varyparam)
    (*ia)[N] = 1;
  else
    (*ia)[N] = 0;
  for(i=0;i<(*Ntovary)+1;i++) {
    (*p)[i][N] = initialval;
  }
  if(varyparam) {
    k = (*Ntovary)+1;
    for(j=0;j<N;j++)
      (*p)[k][j] = (*p)[0][j];
    (*p)[k][N] = initialval + stepval;
    (*Ntovary)++;
  }
  (*Nparameters)++;
}

void amoeba_cleanup(int *Nparameters, int *Ntovary, double ***p, int **ia, double **chi2vals)
{
  int i, k;
  int N, Ns;
  N = *Nparameters;
  if(N <= 0)
    return;
  Ns = ((N - 1) / DEFAULT_INITIAL_PVECSIZE);
  Ns *= DEFAULT_INITIAL_PVECSIZE;
  if(*p != NULL) {
    k = DEFAULT_INITIAL_PVECSIZE + Ns + 1;
    for(i=0; i < k; i++) {
      if((*p)[i] != NULL) {
	free((*p)[i]);
      }
    }
    free((*p));
  }
  if(ia != NULL) {
    if(*ia != NULL)
      free((*ia));
  }
  if(chi2vals != NULL) {
    if(*chi2vals != NULL)
      free((*chi2vals));
  }
}
