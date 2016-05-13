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

/* This file contains functions to linearly decorrelate light curves against various signals for the program vartools by J. Hartman. This version using singular value decomposition to do the decorrelation. */

#define TOLERANCE 1.0E-9

void docorr(double *mag, double *err, int Npoints, int ndecorr, double **decorr, int *order, double *Avector, double *A_errvector, double mag_ave, int zeropoint)
{
  int ncomp, nused, i, j, k, l, s;
  double **Amatrix, *Bvector, **v, *w, *wti, tmp, term, wmax, thresh, sum;

  if(zeropoint)
    ncomp = 1;
  else
    ncomp = 0;
  for(i=0;i<ndecorr;i++)
    ncomp += order[i];
  if(ncomp > 0)
    {
      if((Amatrix = (double **) malloc(Npoints * sizeof(double *))) == NULL ||
	 (Bvector = (double *) malloc(Npoints * sizeof(double))) == NULL ||
	 (v = (double **) malloc(ncomp * sizeof(double *))) == NULL ||
	 (w = (double *) malloc(ncomp * sizeof(double))) == NULL ||
	 (wti = (double *) malloc(ncomp * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0;i<Npoints;i++)
	if((Amatrix[i] = (double *) malloc(ncomp * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
      for(i=0;i<ncomp;i++)
	if((v[i] = (double *) malloc(ncomp * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);

      /* Fill out the design matrix */
      s = 0;
      for(i = 0; i<Npoints; i++)
	{
	  if(!isnan(mag[i]) && !isinf(mag[i]) && err[i] > 0.)
	    {
	      tmp = 1./err[i];
	      Bvector[s] = mag[i]*tmp;
	      if(zeropoint)
		{
		  Amatrix[s][0] = tmp;
		  j=1;
		}
	      else
		j = 0;
	      for(k=0;k<ndecorr;k++)
		{
		  if(!isnan(decorr[i][k]) && !isinf(decorr[i][k]))
		    {
		      term = 1.;
		      for(l=0;l<order[k];l++)
			{
			  term *= decorr[i][k];
			  Amatrix[s][j] = term*tmp;
			  j++;
			}
		    }
		  else
		    {
		      for(l=0;l<order[k];l++)
			{
			  Amatrix[s][j] = 0.;
			  j++;
			}
		    }
		}
	      s++;
	    }
	}
      nused = s;

      /* Run singular value decomposition */
      svdcmp(Amatrix,nused,ncomp,w,v);
      wmax = 0.;
      for(j=0;j<ncomp;j++)
	if(w[j] > wmax) wmax = w[j];
      thresh = TOLERANCE * wmax;
      for(j=0;j<ncomp;j++)
	if(w[j] < thresh) w[j] = 0.0;

      /* Do back substitution */
      svbksb(Amatrix,w,v,nused,ncomp,Bvector,Avector);

      /* Now get the parameter errors */
      for(i=0;i<ncomp;i++)
	{
	  wti[i] = 0.;
	  if(w[i]) wti[i] = 1./(w[i]*w[i]);
	}
      for(i=0;i<ncomp;i++)
	{
	  for(sum=0.,k=0;k<ncomp;k++)
	    sum += v[i][k]*v[i][k]*wti[k];
	  A_errvector[i] = sqrt(sum);
	}

      /* free the allocated memory */
      free(Bvector);
      free(wti);
      free(w);
      for(i=0;i<ncomp;i++)
	free(v[i]);
      free(v);
      for(i=0;i<Npoints;i++)
	free(Amatrix[i]);
      free(Amatrix);
    }
}

void  magcorr(void *t,int ttype, double *mag,double *sig,int Npoints,int ndecorr,double **decorr,int *order,double *coefficients,double *chi2val, double *mean, double mag_ave, int omodel, char *modelname, int zeropoint)
{
  int i, j, k, l, n;
  double term, modelterm;
  double val1, val2;
  FILE *outfile;

  if(omodel)
    {
      if((outfile = fopen(modelname,"w")) == NULL)
	error2(ERR_CANNOTWRITE,modelname);
    }

  for(i=0;i<Npoints;i++)
    {
      if(!isnan(mag[i]) && !isinf(mag[i]))
	{
	  if(zeropoint)
	    {
	      modelterm = (double) coefficients[0];
	      //mag[i] -= (double) coefficients[0];
	      //mag[i] += mag_ave;
	      l = 1;
	    }
	  else
	    {
	      modelterm = 0.;
	      l = 0;
	    }
	  for(j=0;j<ndecorr;j++)
	    {
	      if(!isnan(decorr[i][j]) && !isinf(decorr[i][j]))
		{
		  term = 1.;
		  for(k=0;k<order[j];k++)
		    {
		      term *= decorr[i][j];
		      modelterm += term*((double) coefficients[l]);
		      //mag[i] -= term*((double) coefficients[l]);
		      l++;
		    }
		}
	      else
		{
		  l += order[j];
		}
	    }
	  if(omodel)
	    {
	      switch(ttype)
		{
		case VARTOOLS_TYPE_DOUBLE:
		  fprintf(outfile,"%f %f %f %f\n",((double *) t)[i],mag[i],modelterm,sig[i]);
		  break;
		case VARTOOLS_TYPE_STRING:
		  fprintf(outfile,"%s %f %f %f\n",((char **) t)[i],mag[i],modelterm,sig[i]);
		  break;
		default:
		  error(ERR_BADTYPE);
		}
	    }
	  mag[i] = mag[i] - modelterm + mag_ave;
	}
    }
  if(omodel)
    fclose(outfile);
  val1 = 0.;
  val2 = 0.;
  n = 0;
  for(i=0;i<Npoints;i++)
    {
      if(!isnan(mag[i]))
	{
	  n++;
	  val1 += mag[i];
	  val2 += mag[i]*mag[i];
	}
    }
  val1 /= n;
  val2 /= n;
  *chi2val = 0.;
  *mean = val1;
  for(i=0;i<Npoints;i++)
    {
      if(!isnan(mag[i]) && sig[i] > 0.)
	(*chi2val) += SQR((mag[i] - val1)/sig[i]);
    }
}

void  magcorr_chi2only(double *t,double *mag,double *sig,int Npoints,int ndecorr,double **decorr,int *order,double *coefficients,double *chi2val, double *mean, double mag_ave, int omodel, char *modelname, int zeropoint)
{
  int i, j, k, l, n;
  double term;
  double val1, val2, magtemp;
  FILE *outfile;
  val1 = 0.;
  val2 = 0.;
  if(omodel)
    {
      if((outfile = fopen(modelname,"w")) == NULL)
	error2(ERR_CANNOTWRITE,modelname);
    }
  for(i=0;i<Npoints;i++)
    {
      if(!isnan(mag[i]))
	{
	  magtemp = mag[i];
	  if(zeropoint)
	    {
	      magtemp -= (double) coefficients[0];
	      magtemp += mag_ave;
	      l = 1;
	    }
	  else
	    {
	      magtemp += mag_ave;
	      l = 0;
	    }
	  for(j=0;j<ndecorr;j++)
	    {
	      if(!isnan(decorr[i][j]))
		{
		  term = 1.;
		  for(k=0;k<order[j];k++)
		    {
		      term *= decorr[i][j];
		      magtemp -= term*((double) coefficients[l]);
		      l++;
		    }
		}
	      else
		{
		  l += order[j];
		}
	    }
	  if(omodel)
	    fprintf(outfile,"%f %f %f %f\n",t[i],mag[i],magtemp,sig[i]);
	  val1 += magtemp/(sig[i]*sig[i]);
	  val2 += 1./(sig[i]*sig[i]);
	}
    }
  if(omodel)
    fclose(outfile);
  *mean = val1/val2;
  *chi2val = 0.;
  n = 0;
  for(i=0;i<Npoints;i++)
    {
      if(!isnan(mag[i]))
	{
	  magtemp = mag[i];
	  if(zeropoint)
	    {
	      magtemp -= (double) coefficients[0];
	      magtemp += mag_ave;
	      l = 1;
	    }
	  else
	    {
	      magtemp += mag_ave;
	      l = 0;
	    }
	  for(j=0;j<ndecorr;j++)
	    {
	      if(!isnan(decorr[i][j]))
		{
		  term = 1.;
		  for(k=0;k<order[j];k++)
		    {
		      term *= decorr[i][j];
		      magtemp -= term*((double) coefficients[l]);
		      l++;
		    }
		}
	      else
		{
		  l += order[j];
		}
	    }
	  (*chi2val) += SQR((magtemp - *mean)/sig[i]);
	}
    }
}

/* These ludcmp routines are not used, but saved here in case we decide to switch to them */

void ludcmp(long double **a, int n, int *indx, long double *d)
{
  int i,imax,j,k;
  long double big,dum,sum,temp;
  long double *vv;
  if ((vv = (long double *) malloc(n * sizeof(long double))) == NULL)
    {
      fprintf(stderr,"Memory Allocation Error\n");
      exit(4);
    }
  *d = 1.0;
  for(i=0;i<n;i++)
    {
      big = 0.0;
      for(j=0;j<n;j++)
	if ((temp = ABS_(a[i][j])) > big) big = temp;
      if(big == 0.0)
	{
	  fprintf(stderr,"Singular Matrix in row %d, Aborting!\n",i);
	  exit(4);
	}
      vv[i] = 1.0/big;
    }
  for(j=0;j<n;j++) {
    for(i=0;i<j;i++) {
      sum=a[i][j];
      for(k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for(i=j;i<n;i++) {
      sum=a[i][j];
      for(k=0;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if ( (dum=vv[i]*ABS_(sum)) >= big) {
	big = dum;
	imax = i;
      }
    }
    if ( j != imax) {
      for(k=0;k<n;k++) {
	dum = a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if(a[j][j] == 0.0) a[j][j] = TINY;
    if (j != n) {
      dum = 1.0/(a[j][j]);
      for(i=j+1;i<n;i++) a[i][j] *= dum;
    }
  }
  free(vv);
}

void lubksb(long double **a, int n, int *indx, long double *b)
{
  int i,ii=-1,ip,j;
  long double sum;
  for (i=0; i<n; i++) {
    ip = indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii >= 0)
      for(j=ii;j<=i - 1;j++) sum -= a[i][j]*b[j];
    else if(sum) ii = i;
    b[i] = sum;
  }
  for (i=n-1;i>=0;i--) {
    sum=b[i];
    for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

double fitpoly(int N, double *x, double *y, double *sig, int order,
	       int subtractfit, double *fitparams, double *paramerrs)
/* This function provides easy access to the docorr and magcorr functions
   to fit a polynomial of the form y = sum_i=0^i=order{a_i*x^i}.
   It returns chi2 (or the sum of residuals squared if sig=NULL).

   N - the number of points in the x and y vectors.
   x - vector of independent variables.
   y - vector of dependent variables.
   sig - optional uncertainty vector. If this is NULL then uniform weighting
         of points will be assumed.
   order - order of the polynomial to fit. Must be >= 0.
   subtractfit - 1 to subtract the best-fit model from y (the contents of the
      vector y will be modified in this case). 0 to not subtract it.
   fitparams - The fitted coefficients a_0 through a_N. This must be a vector
         of size order+1, or it can be NULL (in which case the parameters will
         not be returned.)
   paramerrs - The formal uncertainties on the fitted coefficients a_0 through
         a_N. This must also be a vector of size order+1, or it can be NULL (in
         which case the parameters will not be returned.)
*/
{
  int i;
  int Nparam;
  int ordervec;
  double **decorr;
  double *err;
  double *Avector;
  double *Aerrvector;
  double out_chi2;
  double out_mean;

  /* Set up arrays needed for call to docorr */

  if((decorr = (double **) malloc(N * sizeof(double *))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0; i < N; i++) {
    if((decorr[i] = (double *) malloc(sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
    decorr[i][0] = x[i];
  }

  Nparam = order + 1;
  if(sig == NULL) {
    if((err = (double *) malloc(N * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
    for(i=0; i < N; i++) {
      err[i] = 1.;
    }
  } else {
    err = sig;
  }

  if(fitparams == NULL) {
    if((Avector = (double *) malloc(Nparam * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
  } else {
    Avector = fitparams;
  }

  if(paramerrs == NULL) {
    if((Aerrvector = (double *) malloc(Nparam * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
  } else {
    Aerrvector = paramerrs;
  }


  /* Do the fit */
  docorr(y, err, N, 1, decorr, &order, Avector, Aerrvector, 0., 1);

  /* Calculate chi2, and subtract the model if requested */
  if(subtractfit) {
    magcorr(x, 0, y, err, N, 1, decorr, &order, Avector, &out_chi2, &out_mean, 0., 0, NULL, 1);
  } else {
    magcorr_chi2only(x, y, err, N, 1, decorr, &order, Avector, &out_chi2, &out_mean, 0., 0, NULL, 1);
  }

  for(i=0; i < N; i++)
    free(decorr[i]);
  free(decorr);

  if(sig == NULL)
    free(err);
  if(fitparams == NULL)
    free(Avector);
  if(paramerrs == NULL)
    free(paramerrs);

  return(out_chi2);

}
