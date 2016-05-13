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

/* Routines to perform sigmarescaling and ensemble sigmarescaling */

/* Sigma-rescaling function */

void rescalesigma_linear(int N, double *sig, double a, double b)
{
  // Routine to perform the sigma rescaling when doing ensemble rescaling//
  int i;
  for(i=0;i<N;i++)
    {
      if(sig[i] > 0. && !isnan(sig[i]))
	sig[i] = sqrt(sig[i]*sig[i]*a + b);
    }
}

void rescalesigma_chi2(int N, double *sig, double chi2val)
{
  // Routine to perform the sigma rescaling when forcing chi2 to be unity for individual light curves//
  int i;
  double factor;
  factor = sqrt(chi2val);
  for(i=0; i<N; i++)
    {
      if(sig[i] > 0. && !isnan(sig[i]))
	sig[i] = sig[i]*factor;
    }
}

/* Linear fitting functions - used for determining ensemble sigma rescaling */
long double mean(int n, double *data, double *ysig){
  int i;
  int num = 0;
  long double meanval = 0;
  for (i=0;i<n;i++){
    if(!isnan(data[i]) && !isnan(ysig[i]))
      {
	meanval += (long double) (data[i]/ (ysig[i]*ysig[i]));
	num++;
      }
  }
  meanval /= (double) num;

  return meanval;
}

long double mean1(int n, double *ysig){
  int i;
  int num = 0;
  long double meanval = 0;
  for(i=0;i<n;i++){
    if(!isnan(ysig[i]))
      {
	meanval += (long double) (1.0/(ysig[i]*ysig[i]));
	num++;
      }
  }
  meanval /= (double) num;
  return meanval;
}

long double mean2(int n, double *data1, double *data2, double *ysig){
  int i;
  int num = 0;
  long double meanval = 0;
  for (i=0;i<n;i++){
    if(!isnan(data1[i]) && !isnan(data2[i]) && !isnan(ysig[i]))
      {
	meanval += (long double) (data1[i]*data2[i]/ (ysig[i]*ysig[i]));
	num++;
      }
  }
  meanval /= (double) num;
  return meanval;
}

int purge_bad(int *k, double *x, double *y, double *ysig, double a, double b, double max_sig, int use_sig)
{
  int i, j;
  int numbad = 0;
  double stdval = 0.0;
  double *sigsqr;
  sigsqr = (double *) malloc((*k) * sizeof(double));
  if(!use_sig)
    {
      for(i = 0; i < *k; i++)
	sigsqr[i] = (y[i] - (a*x[i]) - b)*(y[i] - (a*x[i]) - b);
      stdval = sqrt(mean(*k,sigsqr,ysig));
      for(i = 0; i < *k; i++)
	sigsqr[i] = sqrt(sigsqr[i])/stdval;
      mysort3(*k, sigsqr, x, y);
    }
  else
    {
      for(i = 0; i < *k; i++)
	sigsqr[i] = sqrt((y[i]-(a*x[i]) - b)*(y[i] - (a*x[i]) - b))/ysig[i];
      mysort4(*k, sigsqr, x, y, ysig);
    }
  j = 0;
  while(sigsqr[j] < max_sig)
    {
      j++;
      if(j == *k) break;
    }
  numbad = *k - j;
  *k = j;
  free(sigsqr);
  return(numbad);
}








