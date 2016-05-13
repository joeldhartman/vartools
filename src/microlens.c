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

/* Below are functions that we will use for fitting only if the GNU scientific library is available */
#ifdef _HAVE_GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

struct _gslfitdata {
  size_t n;
  double * y;
  double * t;
  double * sigma;
};

int chi2microlensmodel_f (const gsl_vector * x, void *data, gsl_vector *f)
{
  size_t n, i;
  double *y, *t, *sigma;
  double F0, F1, u0, t0, tmax;
  double u02, term, u2, Yi;

  n = ((struct _gslfitdata *)data)->n;
  y = ((struct _gslfitdata *)data)->y;
  t = ((struct _gslfitdata *)data)->t;
  sigma = ((struct _gslfitdata *) data)->sigma;

  F0 = gsl_vector_get (x, 0);
  F1 = gsl_vector_get (x, 1);
  u0 = gsl_vector_get (x, 2);
  t0 = gsl_vector_get (x, 3);
  tmax = gsl_vector_get (x, 4);

  u02 = u0*u0;

  for (i = 0; i < n; i++)
    {
      term = (t[i] - tmax) / t0;
      u2 = u02 + (term*term);
      Yi = F0 + F1*(((u2 + 2.)/sqrt(u2)/sqrt(u2 + 4.)) - 1.);
      gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
    }
  return GSL_SUCCESS;
}

int chi2microlensmodel_df (const gsl_vector * x, void *data, gsl_matrix *J)
{

  size_t n, i, s;
  double *y, *t, *sigma;
  double F0, F1, u0, t0, tmax, A, dF0, dF1, du0, dt0, dtmax, dA;
  double u02, term, u2, Yi;

  n = ((struct _gslfitdata *)data)->n;
  y = ((struct _gslfitdata *)data)->y;
  t = ((struct _gslfitdata *)data)->t;
  sigma = ((struct _gslfitdata *) data)->sigma;

  F0 = gsl_vector_get (x, 0);
  F1 = gsl_vector_get (x, 1);
  u0 = gsl_vector_get (x, 2);
  t0 = gsl_vector_get (x, 3);
  tmax = gsl_vector_get (x, 4);

  u02 = u0*u0;

  for (i = 0; i < n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dxj */
      /* where fi = (Yi - yi)/sigma[i];     */
      /*       Yi = microlens function      */
      /* and the xj are the parameters      */
      term = (t[i] - tmax) / t0;
      u2 = u02 + (term*term);
      A = ((u2 + 2.)/sqrt(u2)/sqrt(u2 + 4.));
      s = sigma[i];
      dF0 = 1./s;
      dF1 = (A - 1.)/s;
      dA = (2.*u2*sqrt(u2 + 4.) - (u2 + 2.)*(sqrt(u2 + 4.) + u2/(sqrt(u2 + 4.))))/(u2*(u2 + 4.));
      du0 = F1 * dA * u0 / sqrt(u2) / s;
      dt0 = -F1 * dA * term * term / t0 / sqrt(u2) / s;
      dtmax = -F1 * dA * term / t0 / sqrt(u2) / s;
      gsl_matrix_set (J, i, 0, dF0);
      gsl_matrix_set (J, i, 1, dF1);
      gsl_matrix_set (J, i, 2, du0);
      gsl_matrix_set (J, i, 3, dt0);
      gsl_matrix_set (J, i, 4, dtmax);
    }
}

#endif

void microlensmodel(int N, double *t, double *newlc, double F0, double F1, double u0, double t0, double tmax)
{
  int i;
  double u02, term, u2;
  u02 = u0*u0;
  for(i=0;i<N;i++)
    {
      term = (t[i] - tmax)/t0;
      u2 = u02 + (term*term);
      newlc[i] = F0 + F1*(((u2 + 2.)/sqrt(u2)/sqrt(u2 + 4.)) - 1.);
    }
}

double chi2microlensmodel(double *p, int ma, int N, double *t, double *mag, double *sig, void *userparams)
{
#ifdef PARALLEL
  double *newlc = NULL;
  int size_newlc = 0;
#else
  static double *newlc;
  static int size_newlc = 0;
#endif
  double F0, F1, u0, t0, tmax, chisqval;
  int i;

  F0 = p[0]; F1 = p[1]; u0 = p[2]; t0 = p[3]; tmax = p[4];
  if(N > size_newlc)
    {
      if(!size_newlc)
	{
	  if((newlc = (double *) malloc(N * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      else
	{
	  if((newlc = (double *) realloc(newlc, N * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      size_newlc = N;
    }

  microlensmodel(N, t, newlc, F0, F1, u0, t0, tmax);

  chisqval = 0.;
  for(i=0; i<N; i++)
    {
      newlc[i] = -2.5*log(newlc[i])/log(10.0);
      if(!isnan(mag[i]))
	chisqval += (mag[i] - newlc[i])*(mag[i] - newlc[i])/sig[i]/sig[i];
    }
#ifdef PARALLEL
  if(newlc != NULL) free(newlc);
#endif
  return chisqval;
}

void microlens_initialparams(double *t, double *mag, double *sig, int N, double *F0, double *F1, double *u0, double *t0, double *tmax, int fixf0, int fixf1, int fixu0, int fixt0, int fixtmax)
{
  int i;
  double maxval, minval, tmax_, term;

  /* if tmax is not fixed then find the time of maximum */
  maxval = pow(10.,-0.4*mag[0]);
  minval = pow(10.,-0.4*mag[0]);
  tmax_ = t[0];
  for(i=1;i<N;i++)
    {
      term = pow(10.,-0.4*mag[i]);
      if(term > maxval)
	{
	  maxval = term;
	  tmax_ = t[i];
	}
      if(term < minval)
	minval = term;
    }
  if(!fixtmax)
    *tmax = tmax_;

  /* if f0 and f1 are not fixed, then find the appropriate values */
  if(!fixf0)
    *F0 = minval;
  if(!fixf1)
    *F1 = maxval;

  if(!fixu0)
    *u0 = 1.;

  if(!fixt0)
    {
      *t0 = (t[N-1] - t[0])/10.;
    }
}

#define CONVERGENCELIMIT 0.0001
#define INITIALSTEP 0.1
#define INITIALSTEP_tmax 1


/* This is the function to fit a microlens model to a light curve using downhill simplex. We use the functional form of a microlensing model given by Wozniak, P.R. 2001, AcA, 51, 175 */
void microlens(double *t, double *mag, double *sig, int N, int lc, _MicroLens *m, char *outname, double *f0_out, double *f1_out, double *u0_out, double *t0_out, double *tmax_out, double *chi2_)
{
  int i, j, l, ma, ngood, k, nfunk;
  int *ia, amoeba_val;
  double (*func)(double *, int, int, double *, double *, double *, void *);

  double **p, *y, *initialstep, meanval1, meanval2, ftol;
  double *model, f0, f1, u0, t0, tmax;
  int nvar, Nautofind = 0;

  double *mag_cpy;

  FILE *outfile;

  func = &chi2microlensmodel;

  /* make space for the fitting simplex */
  ma = 5;

  if((ia = (int *) malloc(ma * sizeof(int))) == NULL ||
     (initialstep = (double *) malloc(ma * sizeof(double))) == NULL ||
     (mag_cpy = (double *) malloc(N * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  for(j=0,meanval1=0.0,meanval2 = 0.0; j<N;j++)
    {
      if(!isnan(mag[j]))
	{
	  meanval1 += mag[j]/(sig[j]*sig[j]);
	  meanval2 += 1./(sig[j]*sig[j]);
	}
      meanval1 /= meanval2;
    }

  for(i=0;i<N;i++)
    mag_cpy[i] = mag[i] - meanval1;

  if(m->fitf0) ia[0] = 1; else ia[0] = 0;
  if(m->fitf1) ia[1] = 1; else ia[1] = 0;
  if(m->fitu0) ia[2] = 1; else ia[2] = 0;
  if(m->fitt0) ia[3] = 1; else ia[3] = 0;
  if(m->fittmax) ia[4] = 1; else ia[4] = 0;


  for(nvar=0,j=0;j<ma;j++)
    {
      if(ia[j])
	nvar++;
    }

  if((p = (double **) malloc((nvar + 1) * sizeof(double *))) == NULL ||
     (y = (double *) malloc((nvar + 1) * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  for(j=0;j<nvar+1;j++)
    {
      if((p[j] = (double *) malloc(ma * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }

  /* Get the initial parameters */
  switch(m->f0_source)
    {
    case PERTYPE_SPECIFIED:
      f0 = m->f00[lc][0];
      break;
    case PERTYPE_FIXCOLUMN:
      getoutcolumnvalue(m->f0_linkedcolumn, 0, lc, VARTOOLS_TYPE_DOUBLE, &f0);
      break;
    case PERTYPE_FIX:
      f0 = m->f00_fix;
      break;
    case PERTYPE_AUTOFIND:
      Nautofind = 1;
    }
  switch(m->f1_source)
    {
    case PERTYPE_SPECIFIED:
      f1 = m->f10[lc][0];
      break;
    case PERTYPE_FIXCOLUMN:
      getoutcolumnvalue(m->f1_linkedcolumn, 0, lc, VARTOOLS_TYPE_DOUBLE, &f1);
      break;
    case PERTYPE_FIX:
      f1 = m->f10_fix;
      break;
    case PERTYPE_AUTOFIND:
      Nautofind = 1;
    }
  switch(m->u0_source)
    {
    case PERTYPE_SPECIFIED:
      u0 = m->u00[lc][0];
      break;
    case PERTYPE_FIXCOLUMN:
      getoutcolumnvalue(m->u0_linkedcolumn, 0, lc, VARTOOLS_TYPE_DOUBLE, &u0);
      break;
    case PERTYPE_FIX:
      u0 = m->u00_fix;
      break;
    case PERTYPE_AUTOFIND:
      Nautofind = 1;
    }
  switch(m->t0_source)
    {
    case PERTYPE_SPECIFIED:
      t0 = m->t00[lc][0];
      break;
    case PERTYPE_FIXCOLUMN:
      getoutcolumnvalue(m->t0_linkedcolumn, 0, lc, VARTOOLS_TYPE_DOUBLE, &t0);
      break;
    case PERTYPE_FIX:
      t0 = m->t00_fix;
      break;
    case PERTYPE_AUTOFIND:
      Nautofind = 1;
    }
  switch(m->tmax_source)
    {
    case PERTYPE_SPECIFIED:
      tmax = m->tmax0[lc][0];
      break;
    case PERTYPE_FIXCOLUMN:
      getoutcolumnvalue(m->tmax_linkedcolumn, 0, lc, VARTOOLS_TYPE_DOUBLE, &tmax);
      break;
    case PERTYPE_FIX:
      tmax = m->tmax0_fix;
      break;
    case PERTYPE_AUTOFIND:
      Nautofind = 1;
    }

  if(Nautofind)
    {
      microlens_initialparams(t, mag_cpy, sig, N, &f0, &f1, &u0, &t0, &tmax, !(m->f0_source == PERTYPE_AUTOFIND), !(m->f1_source == PERTYPE_AUTOFIND), !(m->u0_source == PERTYPE_AUTOFIND), !(m->t0_source == PERTYPE_AUTOFIND), !(m->tmax_source == PERTYPE_AUTOFIND));
    }


  /* Initialize the simplex */
  p[0][0] = f0; p[0][1] = f1; p[0][2] = u0; p[0][3] = t0; p[0][4] = tmax;

  /* Set the initial step size */
  for(i=0;i<ma-1;i++)
    {
      if(p[0][i] != 0.)
	initialstep[i] = p[0][i]*INITIALSTEP;
      else
	initialstep[i] = INITIALSTEP;
    }
  initialstep[4] = INITIALSTEP_tmax;

  if(m->f0_initialstep)
    initialstep[0] = m->f0_initialstepval;
  if(m->f1_initialstep)
    initialstep[1] = m->f1_initialstepval;
  if(m->u0_initialstep)
    initialstep[2] = m->u0_initialstepval;
  if(m->t0_initialstep)
    initialstep[3] = m->t0_initialstepval;
  if(m->tmax_initialstep)
    initialstep[4] = m->tmax_initialstepval;

  for(k=1;k<=nvar;k++)
    {
      for(l=0,j=0;j<ma;j++)
	{
	  if(ia[j])
	    {
	      if(l == k - 1)
		{
		  p[k][j] = p[0][j] + initialstep[j];
		}
	      else
		p[k][j] = p[0][j];
	      l++;
	    }
	  else
	    p[k][j] = p[0][j];
	}
    }

  for(k=0;k<nvar+1;k++)
    y[k] = (*func)(p[k],ma,N,t,mag_cpy,sig, NULL);

  ftol = CONVERGENCELIMIT;

  /* Run amoeba if there are any variable terms */
  if(nvar > 0)
    amoeba_val = amoeba(p, y, ia, ma, ftol, func, &nfunk, 0, N, t, mag_cpy, sig, NULL);
  else
    amoeba_val = 0;

  /* If amoeba didn't converge, then write out a garbage model */
  if(amoeba_val)
    {
      *f0_out = -1.; *f1_out = -1.; *u0_out = -1.; *t0_out = -1.; *tmax_out = -1.;
    }
  else
    {
      /* Find the minimum along the vertices */
      k = 0;
      for(j=0;j<=nvar;j++)
	if(y[j] < y[k])
	  k=j;

      *f0_out = p[k][0]; *f1_out = p[k][1]; *u0_out = p[k][2]; *t0_out = p[k][3], *tmax_out = p[k][4];

    }

  /* Remove the signal if we're doing that */
  if(m->correctlc || m->omodel)
    {
      if((model = (double *) malloc(N * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      if(!amoeba_val)
	{
	  microlensmodel(N, t, model, *f0_out, *f1_out, *u0_out, *t0_out, *tmax_out);
	  for(j=0;j<N;j++)
	    model[j] = -2.5*log(model[j])/log(10.0) + meanval1;
	}

      if(m->omodel)
	{
	  if(!strncmp(outname,"-",1) && strlen(outname) == 1)
	    {
	      outfile = stdout;
	    }
	  else
	    {
	      if((outfile = fopen(outname,"w")) == NULL)
		error2(ERR_CANNOTWRITE,outname);
	    }
	  if(!amoeba_val)
	    {
	      for(j=0;j<N;j++)
		{
		  if(!isnan(mag[j]))
		    {
		      fprintf(outfile,"%f %f %f %f\n",t[j],mag[j],model[j], sig[j]);
		    }
		}
	    }
	  else
	    {
	      for(j=0;j<N;j++)
		{
		  if(!isnan(mag[j]))
		    {
		      fprintf(outfile,"%f %f 0. %f\n",t[j],mag[j],sig[j]);
		    }
		}
	    }
	  if(outfile != stdout)
	    fclose(outfile);
	}
      if(m->correctlc && !amoeba_val)
	{
	  for(j=0;j<N;j++)
	    if(!isnan(mag[j]))
	      {
		mag[j] = mag[j] - model[j] + meanval1;
	      }
	}
    }

  ngood = 0;
  for(j=0; j<N; j++)
    if(!isnan(mag[j]))
      ngood++;

  if(ngood > 5)
    *chi2_ = (y[k] / (ngood - nvar));
  else
    *chi2_ = -1.;

  if(m->correctlc || m->omodel)
    free(model);

  free(ia);
  free(y);
  free(initialstep);
  for(j=0;j<=nvar;j++)
    free(p[j]);
  free(p);
}
