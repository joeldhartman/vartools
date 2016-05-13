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

/* Routines to compute chi2 and rms for the program vartools by J. Hartman */

/* Function to calculate binned-chi^2 for a light curve */
int findX(double *x, double xval, int i1, int N)
{
  int u;
  if(x[i1] >= xval)
    return(i1);
  if(x[N-1] < xval)
    return(N);
  while(1)
    {
      u = (N+i1)/2;
      if(u == N || u == i1)
	return(u);
      if(x[u] >= xval)
	{
	  if(x[u-1] < xval)
	    return(u);
	  N = u;
	}
      if(x[u] <= xval)
	{
	  if(u == N-1)
	    return(N);
	  if(x[u+1] >= xval)
	    return(u+1);
	  i1 = u;
	}
    }
}


double binnedchi2(int N, double *t, double *mag, double *sig, double bintime, double *aveval, int *ngood)
{
  int i, jmin, jmax, jminold, jmaxold, foundgood;
  double avesum1, avesum2, mindt, maxdt, chi2val, v;
  double *sumval1, *sumval2, *sumval3, *binmag, *binsig, ave;

  *aveval = -1.;

  *ngood = 0;
  if(N > 0)
    {
      if((sumval1 = (double *) malloc(N * sizeof(double))) == NULL ||
	 (sumval2 = (double *) malloc(N * sizeof(double))) == NULL ||
	 (sumval3 = (double *) malloc(N * sizeof(double))) == NULL ||
	 (binmag = (double *) malloc(N * sizeof(double))) == NULL ||
	 (binsig = (double *) malloc(N * sizeof(double))) == NULL)
	{
	  fprintf(stderr,"Memory Allocation Error\n");
	  exit(2);
	}

      /* First get the average binned magnitude */

      if(!isnan(mag[0]) && sig[0] > 0.)
	{
	  (*ngood)++;
	  sumval1[0] = mag[0] / (sig[0] * sig[0]);
	  sumval2[0] = 1. / (sig[0] * sig[0]);
	  sumval3[0] = sig[0] * sig[0];
	}
      else
	{
	  sumval1[0] = 0.;
	  sumval2[0] = 0.;
	  sumval3[0] = 0.;
	}
      for(i=1;i<N;i++)
	{
	  if(!isnan(mag[i]) && sig[i] > 0.)
	    {
	      (*ngood)++;
	      sumval1[i] = sumval1[i-1] + (mag[i] / (sig[i] * sig[i]));
	      sumval2[i] = sumval2[i-1] + (1. / (sig[i]*sig[i]) );
	      sumval3[i] = sumval3[i-1] + (sig[i] * sig[i]);
	    }
	}

      /* Go through the list find the minimum and maximum times to include via bisection and compute the binned average magnitude and error */
      foundgood = 0;
      jminold = 0;
      jmaxold = 0;
      for(i=0;i<N;i++)
	{
	  mindt = dmax(t[i] - bintime, t[0]);
	  maxdt = dmin(t[i] + bintime, t[N-1]);

	  /*Determine the id of the first point above the minimum t and the first point above the maximum t*/
	  jmin = findX(t,mindt,jminold,N);
	  jmax = findX(t,maxdt,jmaxold,N);

	  jmax = (jmax < N ? (t[jmax] > maxdt ? jmax - 1 : jmax) : N - 1);
	  jminold = jmin;
	  jmaxold = jmax;

	  if(jmin < N && jmin > 0)
	    {
	      if((v = sumval2[jmax] - sumval2[jmin-1]) > 0)
		{
		  binmag[i] = (sumval1[jmax] - sumval1[jmin-1]) / v;
		  binsig[i] = sqrt((sumval3[jmax] - sumval3[jmin-1]) / ((double) (jmax - jmin + 1) * (jmax - jmin + 1)));
		}
	      else
		{
		  binmag[i] = 0.;
		  binsig[i] = 0.;
		}
	    }
	  else if(jmin == 0 && sumval2[jmax] > 0.)
	    {
	      binmag[i] = sumval1[jmax] / sumval2[jmax];
	      binsig[i] = sqrt(sumval3[jmax] / ((double) (jmax + 1) * (jmax + 1)));
	    }
	  else
	    {
	      binmag[i] = 0.;
	      binsig[i] = 0.;
	    }
	}
      avesum1 = 0.;
      avesum2 = 0.;
      for(i=0;i<N;i++)
	{
	  if(binsig[i] > 0.)
	    {
	      avesum1 += binmag[i] / (binsig[i]*binsig[i]);
	      avesum2 += 1. / (binsig[i] * binsig[i]);
	    }
	}
      if(avesum2 > 0)
	{
	  ave = avesum1 / avesum2;
	  *aveval = ave;
	  chi2val = 0.;
	  for(i=0;i<N;i++)
	    {
	      if(binsig[i] > 0.)
		{
		  v = (binmag[i] - ave);
		  v = v*v;
		  chi2val += v / (binsig[i] * binsig[i]);
		}
	    }
	}
      else
	{
	  *aveval = -1.;
	  chi2val = -1.;
	}
      free(sumval1);
      free(sumval2);
      free(sumval3);
      free(binmag);
      free(binsig);
    }
  else
    {
      *aveval = -1.;
      chi2val = -1.;
    }

  return(chi2val);
}

/* Function to calculate un-binned chi2 for a light curve */
double chi2(int N, double *t, double *mag, double *sig, double *aveval, int *ngood)
{
  double avesum1, avesum2, ave, chi2val, v;
  int i;
  *aveval = -1.;
  *ngood = 0;
  if(N > 0)
    {
      avesum1 = 0.;
      avesum2 = 0.;
      for(i=0; i<N; i++)
	{
	  if(!isnan(mag[i]) && sig[i] > 0.)
	    {
	      (*ngood)++;
	      v = sig[i] * sig[i];
	      avesum1 += mag[i] / v;
	      avesum2 += 1. / v;
	    }
	}
      if(avesum2 > 0.)
	{
	  ave = (avesum1 / avesum2);
	  *aveval = ave;
	  chi2val = 0.;
	  for(i=0; i<N; i++)
	    {
	      if(!isnan(mag[i]) && sig[i] > 0.)
		{
		  v = (mag[i] - ave);
		  v = v * v;
		  chi2val += v / (sig[i] * sig[i]);
		}
	    }
	}
      else
	chi2val = -1.;
    }
  else
    chi2val = -1.;
  return(chi2val);
}

/* Function to calculate binned-rms for a light curve */
double binnedrms(int N, double *t, double *mag, double *sig, double bintime, double *aveval, double *rmsthy, int *ngood)
{
  int i, n, jmin, jmax, jminold, jmaxold, foundgood;
  double avesum1, avesum2, avesum3, mindt, maxdt, rmsval, v;
  double *sumval1, *sumval2, *sumval3, *binmag, *binsig, ave;

  *aveval = -1.;
  if(N > 0)
    {
      if((sumval1 = (double *) malloc(N * sizeof(double))) == NULL ||
	 (sumval2 = (double *) malloc(N * sizeof(double))) == NULL ||
	 (sumval3 = (double *) malloc(N * sizeof(double))) == NULL ||
	 (binmag = (double *) malloc(N * sizeof(double))) == NULL ||
	 (binsig = (double *) malloc(N * sizeof(double))) == NULL)
	{
	  fprintf(stderr,"Memory Allocation Error\n");
	  exit(2);
	}

      /* First get the average binned magnitude */

      if(!isnan(mag[0]) && sig[0] > 0.)
	{
	  sumval1[0] = mag[0] / (sig[0] * sig[0]);
	  sumval2[0] = 1. / (sig[0] * sig[0]);
	  sumval3[0] = sig[0] * sig[0];
	}
      else
	{
	  sumval1[0] = 0.;
	  sumval2[0] = 0.;
	  sumval3[0] = 0.;
	}
      for(i=1;i<N;i++)
	{
	  if(!isnan(mag[i]) && sig[i] > 0.)
	    {
	      sumval1[i] = sumval1[i-1] + (mag[i] / (sig[i] * sig[i]));
	      sumval2[i] = sumval2[i-1] + (1. / (sig[i]*sig[i]) );
	      sumval3[i] = sumval3[i-1] + (sig[i] * sig[i]);
	    }
	}

      /* Go through the list find the minimum and maximum times to include via bisection and compute the binned average magnitude and error */
      foundgood = 0;
      jminold = 0;
      jmaxold = 0;
      for(i=0;i<N;i++)
	{
	  mindt = dmax(t[i] - bintime, t[0]);
	  maxdt = dmin(t[i] + bintime, t[N-1]);

	  /*Determine the id of the first point above the minimum t and the first point above the maximum t*/
	  jmin = findX(t,mindt,jminold,N);
	  jmax = findX(t,maxdt,jmaxold,N);

	  jmax = (jmax < N ? (t[jmax] > maxdt ? jmax - 1 : jmax) : N - 1);
	  jminold = jmin;
	  jmaxold = jmax;

	  if(jmin < N && jmin > 0)
	    {
	      if((v = sumval2[jmax] - sumval2[jmin-1]) > 0)
		{
		  binmag[i] = (sumval1[jmax] - sumval1[jmin-1]) / v;
		  binsig[i] = sqrt((sumval3[jmax] - sumval3[jmin-1]) / ((double) (jmax - jmin + 1) * (jmax - jmin + 1)));
		}
	      else
		{
		  binmag[i] = 0.;
		  binsig[i] = 0.;
		}
	    }
	  else if(jmin == 0 && sumval2[jmax] > 0.)
	    {
	      binmag[i] = sumval1[jmax] / sumval2[jmax];
	      binsig[i] = sqrt(sumval3[jmax] / ((double) (jmax + 1)*(jmax + 1)));
	    }
	  else
	    {
	      binmag[i] = 0.;
	      binsig[i] = 0.;
	    }
	}
      avesum1 = 0.;
      avesum2 = 0.;
      avesum3 = 0.;
      n = 0;
      for(i=0;i<N;i++)
	{
	  if(binsig[i] > 0.)
	    {
	      avesum1 += binmag[i];
	      avesum2 += (binmag[i] * binmag[i]);
	      avesum3 += binsig[i] * binsig[i];
	      n++;
	    }
	}
      if(n > 0)
	{
	  ave = avesum1 / (double) n;
	  *aveval = ave;
	  *ngood = n;
	  *rmsthy = sqrt(avesum3 / (double) n);
	  rmsval = sqrt((avesum2 / (double) n) - (ave * ave));
	}
      else
	{
	  *ngood = 0;
	  *rmsthy = -1.;
	  rmsval = -1.;
	}
      free(sumval1);
      free(sumval2);
      free(sumval3);
      free(binmag);
      free(binsig);
    }
  else
    {
      *ngood = 0;
      *rmsthy = -1.;
      rmsval = -1.;
    }

  return(rmsval);
}

/* Function to calculate un-binned rms for a light curve */
double rms(int N, double *t, double *mag, double *sig, double *aveval, double *rmsthy, int *ngood)
{
  double avesum1, avesum2, avesum3, ave, rmsval, n;
  int i;
  *aveval = -1.;
  if(N > 0)
    {
      avesum1 = 0.;
      avesum2 = 0.;
      avesum3 = 0.;
      n = 0;
      for(i=0; i<N; i++)
	{
	  if(!isnan(mag[i]) && sig[i] > 0.)
	    {
	      avesum1 += mag[i];
	      avesum2 += mag[i]*mag[i];
	      avesum3 += sig[i]*sig[i];
	      n++;
	    }
	}
      if(n > 0.)
	{
	  ave = (avesum1 / (double) n);
	  *aveval = ave;
	  *ngood = n;
	  *rmsthy = sqrt(avesum3 / (double) n);
	  rmsval = sqrt((avesum2 / (double) n) - (ave*ave));
	}
      else
	{
	  *ngood = 0;
	  *rmsthy = -1.;
	  rmsval = -1.;
	}
    }
  else
    {
      *ngood = 0;
      *rmsthy = -1.;
      rmsval = -1.;
    }
  return(rmsval);
}


double changeerror(int N, double *t, double *mag, double *sig, double *aveval, int *ngood)
{
  double avesum1, avesum2, avesum3, ave, rmsval, n;
  int i;
  *aveval = -1.;
  if(N > 0)
    {
      avesum1 = 0.;
      avesum2 = 0.;
      avesum3 = 0.;
      n = 0;
      for(i=0; i<N; i++)
	{
	  if(!isnan(mag[i]) && sig[i] > 0.)
	    {
	      avesum1 += mag[i];
	      avesum2 += mag[i]*mag[i];
	      avesum3 += sig[i]*sig[i];
	      n++;
	    }
	}
      if(n > 0.)
	{
	  ave = (avesum1 / (double) n);
	  *aveval = ave;
	  *ngood = n;
	  rmsval = sqrt((avesum2 / (double) n) - (ave*ave));
	  for(i=0; i<N; i++)
	    {
	      sig[i] = rmsval;
	    }
	}
      else
	{
	  *ngood = 0;
	  rmsval = -1.;
	}
    }
  else
    {
      *ngood = 0;
      rmsval = -1.;
    }
  return(rmsval);
}
