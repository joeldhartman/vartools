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

/* This file contains functions for fitting a softened transit model to a light curve, the model comes from

Protopapas, P., Jimenez, R., \& Alcock, C. 2005, MNRAS, 362, 460

*/

#define PIVAL 3.1415926535897932384626433832795
#define TWOPI 6.2831853

#define TOLERANCE 1.0E-9

void softened_transit_(int N, double *t, double *mag, double P, double T0, double eta, double cval, double delta, double mconst, int Nlin_coeffs, double *lin_coeffs, double **Design_Matrix)
{
  int i, j;
  double tanhval1, tanhval2, cosval, sinval, tp, val1;
  for(i=0;i<N;i++)
    {
      val1 = PIVAL*(t[i] - T0)/P;
      sinval = sin(val1);
      cosval = cos(val1);
      tp = P*sinval / PIVAL / eta;
      tanhval1 = tanh(cval*(tp + 0.5));
      tanhval2 = tanh(cval*(tp - 0.5));
      mag[i] = 0.5*delta*(tanhval1 - tanhval2) + mconst;

      for(j=2;j<Nlin_coeffs;j++)
	mag[i] += Design_Matrix[i][j]*lin_coeffs[j];

    }
}

void softened_transit(double *t, double *a_, double *yfit, double **dyda, int ma, int N, int Nlin_coeffs, double *lin_coeffs, double **Design_Matrix, double *y, double *sig, int *varylin_coeffs,void *userparams)
{
  /*
Compute a Protopapas et al. softened transit model for a given set of parameters:

The function is in the format used by the mrqmin non-linear minimization routine from Press et al. 1992.

The input parameters are stored in the vector a_:
They are:
a_[0] = T0, a_[1] = T1, a_[2] = eta, a_[3] = c, a_[4] = N0, a_[5] = N1;

We take P = (T1 - T0)/(N1 - N0), N1 and N0 should not be varied.

see Protopapas' paper for the definition of eta and c.
t is the input times of observations

dyda is a [0.....n-1][0....ma-1] matrix giving the derivative of the model with respect to the parameters at each input time.

yfit is the output light curve
ma is the number of parameters
N is the number of data points

Nlin_coeffs is the number of linear coefficients (2 + 2*(1 + nsubharm + nharm))

lin_coeffs is a vector holding the values of the linear coefficients.

lin_coeffs[0] = delta
lin_coeffs[1] = mconst

Design_Matrix is the N x Nlin_coeffs design matrix used to fit for the linear parameters

y is the N column vector of observed values

sig is the N column vector of uncertainties

varylin_coeffs is a Nlin_coeffs column vector of 1 and 0 denoting whether or not a linear coefficient is free to vary.

userparams is a dummy variable
  */

  int i,j, k, Nlinvary;
  double sechval1, sechval2, tanhval1, tanhval2, dudtp, cosval, sinval, tp, val1, **Amatrix, *Avector, *Bvector, *w, **v, wmax, thresh, dtpdt0, dtpdt1;
  double P, T0;

  /* First find the number of linear coefficients that can vary */
  for(i=0, Nlinvary = 0;i<Nlin_coeffs;i++)
    if(varylin_coeffs[i])
      Nlinvary++;

  /* Allocate memory for the matrix and vectors used to solve the linear problem */
  if(Nlinvary > 0)
    {
      if((Amatrix = (double **) malloc(N * sizeof(double *))) == NULL ||
	 (Avector = (double *) malloc(Nlinvary * sizeof(double))) == NULL ||
	 (Bvector = (double *) malloc(N * sizeof(double))) == NULL ||
	 (v = (double **) malloc(Nlinvary * sizeof(double *))) == NULL ||
	 (w = (double *) malloc(Nlinvary * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0;i<N;i++)
	if((Amatrix[i] = (double *) malloc(Nlinvary * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
      for(i=0;i<Nlinvary;i++)
	if((v[i] = (double *) malloc(Nlinvary * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);

      /* Fill out the fitting design matrix */
      for(i=0;i<N;i++)
	Bvector[i] = y[i];

      for(k=0,j=0;k<Nlin_coeffs;k++)
	{
	  if(varylin_coeffs[k] && k < 1)
	    {
	      j++;
	    }
	  else if(varylin_coeffs[k] && k == 1)
	    {
	      for(i=0;i<N;i++)
		{
		  Amatrix[i][j] = 1./sig[i];
		}
	      j++;
	    }
	  else if(varylin_coeffs[k] && k > 1)
	    {
	      for(i=0;i<N;i++)
		{
		  Amatrix[i][j] = Design_Matrix[i][k]/sig[i];
		}
	      j++;
	    }
	  else if(!varylin_coeffs[k])
	    {
	      for(i=0;i<N;i++)
		Bvector[i] -= lin_coeffs[k]*Design_Matrix[i][k];
	    }
	}
    }

  P = (a_[1] - a_[0])/(a_[5] - a_[4]);
  T0 = a_[0];
  /* Do the Transit Model */
  for(i=0;i<N;i++)
    {
      val1 = PIVAL*(t[i] - T0)/P;
      sinval = sin(val1);
      cosval = cos(val1);
      tp = P*sinval / PIVAL / a_[2];
      tanhval1 = tanh(a_[3]*(tp + 0.5));
      tanhval2 = tanh(a_[3]*(tp - 0.5));
      yfit[i] = 0.5*(tanhval1 - tanhval2);
      sechval1 = cosh(a_[3]*(tp + 0.5));
      sechval1 = 1./sechval1/sechval1;
      sechval2 = cosh(a_[3]*(tp - 0.5));
      sechval2 = 1./sechval2/sechval2;
      dudtp = 0.5*(a_[3]*(sechval1 - sechval2));
      dtpdt0 = -tp/(a_[5] - a_[4])/P + cosval*((t[i] - T0)/P/(a_[5] - a_[4]) - 1.)/a_[2];
      dtpdt1 = tp/(a_[5] - a_[4])/P - cosval*(t[i] - T0)/P/(a_[5] - a_[4])/a_[2];
      //dyda[i][0] = dudtp * (sinval - (t[i] - a_[1])*cosval / a_[0]) / PIVAL / a_[2];
      //dyda[i][1] = -dudtp * cosval / a_[2];
      dyda[i][0] = dudtp*dtpdt0;
      dyda[i][1] = dudtp*dtpdt1;
      dyda[i][2] = -dudtp * P * sinval / PIVAL / a_[2] / a_[2];
      dyda[i][3] = 0.5*((tp + 0.5)*sechval1 - (tp - 0.5)*sechval2);
      //dyda[i][4] = yfit[i] / a_[4];
      //yfit[i] += a_[5];
      //dyda[i][5] = 1.;
    }
  /*Add the delta term */
  if(varylin_coeffs[0])
    {
      for(i=0;i<N;i++)
	Amatrix[i][0] = yfit[i]/sig[i];
    }
  else if(Nlinvary > 0)
    {
      for(i=0;i<N;i++)
	Bvector[i] -= yfit[i]*lin_coeffs[0];
    }
  if(Nlinvary > 0)
    {
      for(i=0;i<N;i++)
	Bvector[i] /= sig[i];

      /* Run singular value decomposition */
      svdcmp(Amatrix,N,Nlinvary,w,v);
      wmax = 0.;
      for(j=0;j<Nlinvary;j++)
	if(w[j] > wmax) wmax = w[j];
      thresh = TOLERANCE * wmax;
      for(j=0;j<Nlinvary;j++)
	if(w[j] < thresh) w[j] = 0.0;
      /* Do back substitution */
      svbksb(Amatrix,w,v,N,Nlinvary,Bvector,Avector);

      /* Update the lin_coeffs vector with the values from Avector */
      for(i=0,j=0;i<Nlin_coeffs;i++)
	{
	  if(varylin_coeffs[i])
	    {
	      lin_coeffs[i] = Avector[j];
	      j++;
	    }
	}
    }
  /* Now adjust the yfit and dyda values to take into account the linear results */
  for(i=0;i<N;i++)
    {
      yfit[i] = yfit[i]*lin_coeffs[0] + lin_coeffs[1];
      for(j=0;j<4;j++)
	dyda[i][j] = dyda[i][j]*lin_coeffs[0];
    }
  if(Nlin_coeffs > 2)
    {
      for(j=2;j<Nlin_coeffs;j++)
	{
	  for(i=0;i<N;i++)
	    yfit[i] += lin_coeffs[j]*Design_Matrix[i][j];
	}
    }

  /* Free the temporary storage arrays */
  if(Nlinvary > 0)
    {
      free(Avector);
      free(Bvector);
      free(w);
      for(i=0;i<N;i++)
	free(Amatrix[i]);
      for(i=0;i<Nlinvary;i++)
	free(v[i]);
      free(Amatrix);
      free(v);
    }
}

#define CONVERGENCELIMIT 0.0001
#define INITIALSTEP 0.1

void fitsoftened_transit(int N, double *t, double *mag, double *sig, double *P, double *T0, double *eta, double *cval, double *delta, double *mconst, int fitephem, int fiteta, int fitcval, int fitdelta, int fitmconst, double *chi2_, int correctlc, int omodel, char *modelname, int dokillharm, int nharm, int nsubharm, double perharm, double *subharmA, double *subharmB, double *harmA, double *harmB, double *fundA, double *fundB)
{
  /* This function is called by processcommand, it fits the Protopapas softened transit model to a light curve.
   */

  int i, j, k, ma, ngood, nvar, iter, Nlin_coeff, *varylin_coeffs, freqnum;
  double alamda, *a_, **covar, **alpha_, chisq, ochisq, deltachi, *delmag, oalamda, **Design_Matrix, *lin_coeffs, freq, ct, st, st_, ct_, d1, d2, meanval1, meanval2, phase, N0, N1, Tbaseline;
  int *ia;
  FILE *outfile;
  void (*funcs)(double *, double *, double *, double **, int, int, int, double *, double **, double *, double *, int *, void *);

#ifdef PARALLEL
  int mfit;
  double ochisq_;
  double *atry, *beta, *da, **oneda;
#endif

  funcs = &softened_transit;

  if(*P <= 0.)
    {
      return;
    }

  ma = 6;

  Nlin_coeff = 2;

  if(dokillharm)
    Nlin_coeff += 2*(1 + nharm + nsubharm);

  if((ia = (int *) malloc(ma * sizeof(int))) == NULL ||
     (covar = (double **) malloc(ma * sizeof(double *))) == NULL ||
     (alpha_ = (double **) malloc(ma * sizeof(double *))) == NULL ||
     (a_ = (double *) malloc(ma * sizeof(double))) == NULL ||
     (Design_Matrix = (double **) malloc(N * sizeof(double *))) == NULL ||
     (varylin_coeffs = (int *) malloc(Nlin_coeff * sizeof(int))) == NULL ||
     (lin_coeffs = (double *) malloc(Nlin_coeff * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  for(j=0;j<ma;j++)
    if((covar[j] = (double *) malloc(ma * sizeof(double))) == NULL ||
       (alpha_[j] = (double *) malloc(ma * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);

  for(i=0;i<N;i++)
    if((Design_Matrix[i] = (double *) malloc(Nlin_coeff * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);

  alamda = -1.;
  a_[0] = *T0;
  N0 = 0.;
  a_[4] = N0;
  Tbaseline = t[N-1] - t[0];
  N1 = floor(Tbaseline / (*P));
  if(N1 <= 0.)
    N1 = 1.;
  a_[5] = N1;
  a_[1] = a_[0] + (*P)*N1;

  //a_[0] = *P; a_[1] = *T0;

  a_[2] = *eta; a_[3] = *cval;

  lin_coeffs[0] = *delta;
  lin_coeffs[1] = *mconst;

  if(fitephem) ia[0] = 1; else ia[0] = 0;
  if(fitephem) ia[1] = 1; else ia[1] = 0;
  if(fiteta) ia[2] = 1; else ia[2] = 0;
  if(fitcval) ia[3] = 1; else ia[3] = 0;
  ia[4] = 0; ia[5] = 0;

  if(fitdelta) varylin_coeffs[0] = 1; else varylin_coeffs[0] = 0;
  if(fitmconst) varylin_coeffs[1] = 1; else varylin_coeffs[1] = 0;

  for(nvar=0,j=0;j<ma;j++)
    if(ia[j])
      nvar++;

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

      lin_coeffs[1] = *mconst;
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

  /* Fill out the Design Matrix if dokillharm is set */
  if(dokillharm)
    {
      for(i=0;i<N;i++)
	{
	  j = 2;
	  freq = TWOPI / perharm;
	  for(k=0;k<nsubharm;k++)
	    {
	      freqnum = (nsubharm - k) + 1;
	      st = sin(t[i] * (freq / (double) freqnum));
	      ct = cos(t[i] * (freq / (double) freqnum));
	      Design_Matrix[i][j] = st;
	      varylin_coeffs[j] = 1;
	      j++;
	      Design_Matrix[i][j] = ct;
	      varylin_coeffs[j] = 1;
	      j++;
	    }
	  st = sin(t[i] * freq);
	  ct = cos(t[i] * freq);
	  st_ = st;
	  ct_ = ct;
	  Design_Matrix[i][j] = st;
	  varylin_coeffs[j] = 1;
	  j++;
	  Design_Matrix[i][j] = ct;
	  varylin_coeffs[j] = 1;
	  j++;
	  for(k=0;k<nharm;k++)
	    {
	      d1 = ct_*st + st_*ct;
	      d2 = ct_*ct - st_*st;
	      st_ = d1;
	      ct_ = d2;
	      Design_Matrix[i][j] = st_;
	      varylin_coeffs[j] = 1;
	      j++;
	      Design_Matrix[i][j] = ct_;
	      varylin_coeffs[j] = 1;
	      j++;
	    }
	}
    }

  /* Initialize mrqmin */
#ifdef PARALLEL
  if((atry = (double *) malloc(ma * sizeof(double))) == NULL ||
     (beta = (double *) malloc(ma * sizeof(double))) == NULL ||
     (da = (double *) malloc(ma * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  for(mfit=0,k=0;k<ma;k++)
    if(ia[k]) mfit++;
  if((oneda = (double **) malloc(mfit * sizeof(double *))) == NULL)
    error(ERR_MEMALLOC);
  for(k=0;k<mfit;k++)
    if((oneda[k] = (double *) malloc(sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
  mrqmin(t, mag, sig, N, a_, ia, ma, covar, alpha_, &chisq, &alamda, funcs, Nlin_coeff,Design_Matrix,lin_coeffs,varylin_coeffs,mfit,&ochisq_,atry,beta,da,oneda,NULL);
#else
  mrqmin(t, mag, sig, N, a_, ia, ma, covar, alpha_, &chisq, &alamda, funcs, Nlin_coeff,Design_Matrix,lin_coeffs,varylin_coeffs,NULL);
#endif

  ochisq = chisq;

  iter = 0;

  do
    {
      oalamda = alamda;
#ifdef PARALLEL
      mrqmin(t, mag, sig, N, a_, ia, ma, covar, alpha_, &chisq, &alamda, funcs, Nlin_coeff,Design_Matrix,lin_coeffs,varylin_coeffs,mfit,&ochisq_,atry,beta,da,oneda,NULL);
#else
      mrqmin(t, mag, sig, N, a_, ia, ma, covar, alpha_, &chisq, &alamda, funcs, Nlin_coeff,Design_Matrix,lin_coeffs,varylin_coeffs,NULL);
#endif
      deltachi = (ochisq - chisq);
      ochisq = chisq;
      iter++;
    }
  while((deltachi < 0. || deltachi > CONVERGENCELIMIT || iter < 10 || alamda > oalamda) && alamda != 0.0);

  /* Run a final pass to get the covariance matrix */
  alamda = 0.0;
#ifdef PARALLEL
  mrqmin(t, mag, sig, N, a_, ia, ma, covar, alpha_, &chisq, &alamda, funcs, Nlin_coeff,Design_Matrix,lin_coeffs,varylin_coeffs,mfit,&ochisq_,atry,beta,da,oneda,NULL);
  for(k=0;k<mfit;k++)
    free(oneda[k]);
  free(oneda);
  free(da);
  free(beta);
  free(atry);
#else
  mrqmin(t, mag, sig, N, a_, ia, ma, covar, alpha_, &chisq, &alamda, funcs, Nlin_coeff,Design_Matrix,lin_coeffs,varylin_coeffs,NULL);
#endif

  /* Update the parameters */
  *T0 = a_[0]; *P = (a_[1] - a_[0])/(N1 - N0);
  //*P = a_[0]; *T0 = a_[1];
  *eta = a_[2]; *cval = a_[3];
  *delta = lin_coeffs[0]; *mconst = lin_coeffs[1];

  if(dokillharm)
    {
      j = 2;
      for(k = 0; k<nsubharm;k++)
	{
	  subharmA[k] = lin_coeffs[j];
	  j++;
	  subharmB[k] = lin_coeffs[j];
	  j++;
	}
      *fundA = lin_coeffs[j];
      j++;
      *fundB = lin_coeffs[j];
      j++;
      for(k = 0;k<nharm;k++)
	{
	  harmA[k] = lin_coeffs[j];
	  j++;
	  harmB[k] = lin_coeffs[j];
	  j++;
	}
    }

  /* Remove the signal if we're doing that */
  if(correctlc || omodel)
    {
      if((delmag = (double *) malloc(N * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      softened_transit_(N, t, delmag, *P, *T0, *eta, *cval, *delta, *mconst, Nlin_coeff,lin_coeffs,Design_Matrix);
      if(omodel)
	{
	  if((outfile = fopen(modelname,"w")) == NULL)
	    error2(ERR_CANNOTWRITE,modelname);
	  fprintf(outfile,"#Time Mag ModelMag eMag Phase\n");
	  for(j=0;j<N;j++)
	    if(!isnan(mag[j]))
	      {
		phase = (t[j] - *T0)/(*P);
		phase -= (double) floor(phase);
		fprintf(outfile,"%f %f %f %f %f\n",t[j],mag[j],delmag[j],sig[j],phase);
	      }
	  fclose(outfile);
	}
      if(correctlc)
	{
	  for(j=0;j<N;j++)
	    if(!isnan(mag[j]))
	      { mag[j] = mag[j] - delmag[j]  + meanval1; }
	}
    }

  ngood = 0;
  for(j=0;j<N;j++)
    if(!isnan(mag[j]))
      ngood++;

  if(ngood > 7)
    *chi2_ = (chisq / (ngood - nvar));
  else
    *chi2_ = -1.;
  if(correctlc || omodel)
    free(delmag);
  free(ia);
  free(a_);
  for(j=0;j<ma;j++)
    {
      free(covar[j]);
      free(alpha_[j]);
    }
  free(covar);
  free(alpha_);
  for(i=0;i<N;i++)
    free(Design_Matrix[i]);
  free(Design_Matrix);
  free(lin_coeffs);
  free(varylin_coeffs);
}


