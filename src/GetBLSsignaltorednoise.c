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
/* This file contains functions to calculate the white noise, the red noise, the depth, the number of in transit points (nt), the number of distinct transits sampled (Nt) and the signal to pink noise (delta^2 / ((sigma_w ^2 / nt) + (sigma_r ^2 / Nt))). for a given BLS run */

#include "commands.h"
#include "programdata.h"
#include "functions.h"

#define BLS_NTV_FIELD_LENGTH 5

double getrms(int N, double *t, double *mag, double *sig, double *aveval, double *rmsthy, int *ngood)
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

double getweightedrms(int N, double *t, double *mag, double *sig, double *aveval, double *rmsthy, int *ngood)
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
	      avesum1 += mag[i]/sig[i]/sig[i];
	      avesum2 += 1.0/sig[i]/sig[i];
	      avesum3 += sig[i]*sig[i];
	      n++;
	    }
	}
      if(n > 0.)
	{
	  ave = (avesum1)/(avesum2);
	  avesum1 = 0.0;
	  for(i=0; i < N; i++)
	    {
	      if(!isnan(mag[i]) && sig[i] > 0.)
		{
		  avesum1 += ((mag[i]-ave)*(mag[i]-ave))/sig[i]/sig[i];
		}
	    }
	  *aveval = ave;
	  if(n > 1.) {
	    rmsval = sqrt( (avesum1)/((n - 1)*avesum2/n));
	    *rmsthy = sqrt(avesum3 / (double) n);
	  } else {
	    rmsval = -1.;
	    *rmsthy = -1.;
	  }
	  *ngood = n;
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


void redwhitenoise(int N, double *t, double *mag, double *sig, double timespan, double *rednoise, double *whitenoise)
{
  double rmsbinval, rmsbinthy, binaveval, rmsval, rmsthy, aveval, expectedbinrms;
  int ngood, ngoodbin;
  
  rmsval = getweightedrms(N, t, mag, sig, &aveval, &rmsthy, &ngood);
  rmsbinval = binnedrms(N, t, mag, sig, timespan, &binaveval, &rmsbinthy, &ngoodbin);
  
  expectedbinrms = rmsval * rmsbinthy / rmsthy;
  
  *whitenoise = rmsval;
  *rednoise = rmsbinval*rmsbinval - expectedbinrms*expectedbinrms;
  if((*rednoise) > 0)
    *rednoise = sqrt((rmsbinval*rmsbinval) - (expectedbinrms*expectedbinrms));
  else
    *rednoise = 0.;
}

void subtractbls(int N, double *t, double *mag, double *sig, double P, double q, double depth, double ph1, int *nt, int *Nt, int *nbefore, int *nafter, double qingress, double OOTmag, int *ntv)
{

  double t1, t0, ph, f0, phb1, phb2, pha1, pha2, ph2;
  int maxNt, i, j, k, len, val;
#ifdef PARALLEL
  int *Npint = NULL;
  int sizeNpint = 0;
#else
  static int *Npint;
  static int sizeNpint = 0;
#endif

  maxNt = ceil((t[N-1]-t[0])/P);
  if(sizeNpint == 0)
    {
      sizeNpint = maxNt + 1;
      if((Npint = (int *) malloc(sizeNpint * sizeof(int))) == NULL)
	error(ERR_MEMALLOC);
    }
  else if(maxNt + 1 > sizeNpint)
    {
      sizeNpint = maxNt + 1;
      if((Npint = (int *) realloc(Npint, sizeNpint * sizeof(int))) == NULL)
	error(ERR_MEMALLOC);
    }
  for(i = 0; i<maxNt; i++)
    Npint[i] = 0;
  t0 = t[0];
  t1 = t[0] + ph1*P;
  if(ph1 > 1. - q)
    t1 = t1 - P;
  f0 = 1./P;
  *nt = 0;
  *Nt = 0;
  *nbefore = 0;
  *nafter = 0;
  //ph1 = (double) in1 / (double) nb;
  //ph2 = (double) in2 / (double) nb;
  
  ph2 = q;
  phb1 = qingress*q;
  phb2 = q - qingress*q;

  for(i=0; i<N; i++)
    {
      ph = (t[i] - t1)*f0;
      j = (int) ((t[i] - t1)*f0);
      ph = ph - floor(ph);
      if(ph < q) {
	Npint[j]++;
	(*nt) += 1;
	if(ph >= phb1 && ph <= phb2) {
	  mag[i] = mag[i] - depth;
	} else if(ph < phb1) {
	  mag[i] = mag[i] - ph*depth/q/qingress;
	} else {
	  mag[i] = mag[i] - (q - ph)*depth/q/qingress;
	}
      }
      else if(ph < (2.0*q)) {
	(*nafter)++;
      }
      else if(ph > (1. - q)) {
	(*nbefore)++;
      }
    }
      

  /*
  if(in2 > in1)
    {
      phb2 = ph1; phb1 = phb2 - (ph2 - ph1);
      pha1 = ph2; pha2 = pha1 + (ph2 - ph1);
      for(i=0;i<N;i++)
	{
	  ph = (t[i]-t1)*f0;
	  j = (int) ph;
	  ph -= (int) ph;
	  k = (int) (nb * ph);
	  if(k >= in1 && k <= in2)
	    {
	      mag[i] = mag[i] - depth;
	      (*nt) += 1;
	      Npint[j]++;
	    }
	  else if((ph > pha1 && ph < pha2) || (ph - 1 > pha1 && ph -1 < pha2))
	    {
	      (*nafter)++;
	    }
	  else if((ph > phb1 && ph < phb2) || (ph + 1 > phb1 && ph + 1 < phb2))
	    {
	      (*nbefore)++;
	    }
	}
    }
  else
    {
      phb2 = ph1; phb1 = phb2 - (ph2 + 1. - ph1);
      pha1 = ph2; pha2 = pha1 + (ph2 + 1. - ph1);
      for(i=0;i<N;i++)
	{
	  ph = (t[i]-t1)*f0;
	  j = (int) ph;
	  ph -= (int) ph;
	  k = (int) (nb * ph);
	  if(k <= in2)
	    {
	      mag[i] = mag[i] - depth;
	      (*nt) += 1;
	      Npint[j]++;
	    }
	  else if( k >= in1)
	    { 
	      mag[i] = mag[i] - depth;
	      (*nt) += 1;
	      Npint[j+1]++;
	    }
	  else if((ph > pha1 && ph < pha2))
	    {
	      (*nafter)++;
	    }
	  else if((ph > phb1 && ph < phb2))
	    {
	      (*nbefore)++;
	    }
	}
    }
  */
  for(i=0;i<maxNt;i++)
    {
      if(Npint[i])
	(*Nt)++;
    }

  if(ntv != NULL) {
    mysort1int(maxNt,Npint);
    len = maxNt;
    (*ntv) = 0;
    for(i=0; i < BLS_NTV_FIELD_LENGTH; i++) {
      if(0 < (len-(i+2))) {
	val = 10*(((float)(Npint[len-(i+2)]))/Npint[len-1]);
	if(val > 10) val = 9;
	(*ntv) += pow(10,BLS_NTV_FIELD_LENGTH-(i+1))*val;
      }
    }
  }
#ifdef PARALLEL
  if(Npint != NULL) free(Npint);
#endif

}

void getsignaltopinknoiseforgivenblsmodel(int N, double *t, double *mag, double *sig, double P, double q, double depth, double in1_ph, int *nt, int *Nt, int *Nbefore, int *Nafter, double *rn, double *wn, double *sigtopink, double qingress, double OOTmag, int *ntv)
{
#ifdef PARALLEL
  int size_vec = 0;
  double *tstore = NULL, *magstore = NULL, *sigstore = NULL;
#else
  static int size_vec = 0;
  static double *tstore, *magstore, *sigstore;
#endif

  double duration;



  /* Allocate memory for the temporary light curve vectors if necessary - note that we define these as static so we don't have to reallocate memory for them every time we call this function - there may be a memory leak though as a result if a user includes an arbitrarily large number of -BLS commands on the command-line */

  if(size_vec == 0)
    {
      size_vec = N;
      if((tstore = (double *) malloc(N * sizeof(double))) == NULL ||
	 (magstore = (double *) malloc(N * sizeof(double))) == NULL ||
	 (sigstore = (double *) malloc(N * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
  else if(N > size_vec)
    {
      size_vec = N;
      if((tstore = (double *) realloc(tstore, N * sizeof(double))) == NULL ||
	 (magstore = (double *) realloc(magstore, N * sizeof(double))) == NULL ||
	 (sigstore = (double *) realloc(sigstore, N * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
  
  /* Copy the light curve to the temporary light curve */
  memcpy(tstore,t,(N * sizeof(double)));
  memcpy(magstore,mag,(N * sizeof(double)));
  memcpy(sigstore,sig,(N * sizeof(double)));
  
  subtractbls(N, tstore, magstore, sigstore, P, q, depth, in1_ph, nt, Nt, Nbefore, Nafter, qingress, OOTmag, ntv);

  duration = q*P;
	      
  redwhitenoise(N, tstore, magstore, sigstore, duration, rn, wn);

  *sigtopink = depth*depth/(((*wn)*(*wn)/(double)(*nt)) + ((*rn)*(*rn)/(double)(*Nt)));
  *sigtopink = sqrt((*sigtopink));
#ifdef PARALLEL
  if(tstore != NULL) free(tstore);
  if(magstore != NULL) free(magstore);
  if(sigstore != NULL) free(sigstore);
#endif
}

	    
