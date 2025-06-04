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

#define MAX_NHARM_AUTO 50

void AOV_getfixedperiodSNR(double *p, double *pr, int N, double ave, double rms, double Pin, double *peak, double *SNR)
{
  int i;
  i = findX(p, Pin, 0, N);
  if(i > 0 && i < N-1)
    {
      (*peak) = (pr[i]*(Pin - p[i]) + pr[i+1]*(p[i+1] - Pin))/(p[i+1] - p[i]);
    }
  else if(i == 0)
    {
      if(Pin > p[0])
	{
	  (*peak) = (pr[i]*(Pin - p[i]) + pr[i+1]*(p[i+1] - Pin))/(p[i+1] - p[i]);
	}
      else
	(*peak) = pr[0];
    }
  else if(i == N)
    (*peak) = pr[N - 1];
  (*SNR) = ((*peak) - ave)/rms;
}

/* This file contains functions for computing the AoV periodogram and period peaks using harmonics (Schwarzenberg-Czerny, 1996) for the program vartools by J. Hartman */



#define TWOPID 6.2831853071795865
#define ERROR_SCORE -100000.0  // assumed to be less than any valid score (>> 1)
#define MAX_DOUBLE_CHECK_MULTIPLE 19
#define MAX_PERIOD_DIFF_MULTIPLE 5

double TestPeriod_aov_harm(int N, double *t, double *m, double *sig, int Nharm, double testperiod, double *m_noave, double *t_nostart, double *weight, double lcvariance, int *Nharm_used)
{
  /* This is an implementation of the multi-harmonic AoV periodogram (Schwarzenberg-Czerny, 1996, ApJ, 460, L107.

This function returns AoV at a specific period

The input values are:
N - the number of data points
t - vector of times
m - vector of magnitudes
m_noave - a storage vector of size N
sig - vector of magnitude uncertainties
Nharm - number of harmonics to include (if Nharm is < 1 then it will be automatically varied until the FAP is minimized).
testperiod - period to calculate the aov value at

m_noave, t_nostart, weight, lcvariance are vectors created by aov_harm.

c_r, c_i are storage vectors of size 2*Nharm
psi_r and psi_i are storage vectors of size N
z_r and z_i are storage vectors of size N
phi_r and phi_i are storage vectors of size N

Note that memory for periodogram should be allocated before calling this function.
  */
  int i, n, n_, Nharmtrial;
  double freq, Nharm_orig_d, alpha_n_r, alpha_n_i;
  double var1, var1_r, var1_i, var2_r, var2_i, var3_r, var3_i, var4_r, var4_i, var5;
  double Ntimesfreq, tmp1, tmp2, th_coeff1, th_coeff2, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, periodogram, periodogram_new, a, b;

#ifdef PARALLEL
  int sizeNvecs = 0, sizeNharmvecs = 0;
  double *c_r = NULL, *c_i = NULL, *psi_r = NULL, *psi_i = NULL, 
    *z_r = NULL, *z_i = NULL, *phi_r = NULL, *phi_i = NULL, *zn_r = NULL, 
    *zn_i = NULL;
#else
  static int sizeNvecs = 0, sizeNharmvecs = 0;
  static double *c_r = NULL, *c_i = NULL, *psi_r = NULL, *psi_i = NULL, 
    *z_r = NULL, *z_i = NULL, *phi_r = NULL, *phi_i = NULL, *zn_r = NULL, 
    *zn_i = NULL;
#endif

  /* change the frequency to angular frequency for the computation, also double Nharm to account for the complex polynomials being of double the order of the real space polynomials */

  freq = TWOPID / testperiod;
  
  if(Nharm > 0)
    {
      Nharm_orig_d = (double) Nharm;
      Nharm = Nharm * 2;
      
      th_coeff1 = (double) (N - Nharm - 1);
      th_coeff2 = (double) Nharm;
    }

  /* Initialize memory for the storage vectors if this is our first pass, otherwise if the old vectors are not long enough allocation additional memory */
  if(!sizeNvecs)
    {
      sizeNvecs = N;
      if((psi_r = (double *) malloc(N * sizeof(double))) == NULL ||
	 (psi_i = (double *) malloc(N * sizeof(double))) == NULL ||
	 (z_r = (double *) malloc(N * sizeof(double))) == NULL ||
	 (z_i = (double *) malloc(N * sizeof(double))) == NULL ||
	 (zn_r = (double *) malloc(N * sizeof(double))) == NULL ||
	 (zn_i = (double *) malloc(N * sizeof(double))) == NULL ||
	 (phi_r = (double *) malloc(N * sizeof(double))) == NULL ||
	 (phi_i = (double *) malloc(N * sizeof(double))) == NULL)
	//error(ERR_MEMALLOC);
	exit(3);
    }
  else if(sizeNvecs < N)
    {
      sizeNvecs = N;
      if((psi_r = (double *) realloc(psi_r, N * sizeof(double))) == NULL ||
	 (psi_i = (double *) realloc(psi_i, N * sizeof(double))) == NULL ||
	 (z_r = (double *) realloc(z_r, N * sizeof(double))) == NULL ||
	 (z_i = (double *) realloc(z_i, N * sizeof(double))) == NULL ||
	 (zn_r = (double *) realloc(zn_r, N * sizeof(double))) == NULL ||
	 (zn_i = (double *) realloc(zn_i, N * sizeof(double))) == NULL ||
	 (phi_r = (double *) realloc(phi_r, N * sizeof(double))) == NULL ||
	 (phi_i = (double *) realloc(phi_i, N * sizeof(double))) == NULL)
	//error(ERR_MEMALLOC);
	exit(3);
    }
  if(!sizeNharmvecs)
    {
      if(Nharm > 0)
	{
	  sizeNharmvecs = Nharm + 1;
	  if((c_r = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL ||
	     (c_i = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL)
	    //error(ERR_MEMALLOC);
	    exit(3);
	}
      else
	{
	  sizeNharmvecs = MAX_NHARM_AUTO + 1;
	  if((c_r = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL ||
	     (c_i = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL)
	    exit(3);
	}
    }
  else if(sizeNharmvecs < Nharm + 1)
    {
      sizeNharmvecs = Nharm + 1;
      if((c_r = (double *) realloc(c_r, sizeNharmvecs * sizeof(double))) == NULL ||
	 (c_i = (double *) realloc(c_i, sizeNharmvecs * sizeof(double))) == NULL)
	//error(ERR_MEMALLOC);
	exit(3);
    }

  /* Calculate aov at the test frequency */

  if(Nharm > 0)
    {
      /* Do the calculation for fixed Nharm */
      /* For the trial frequency calculate the z, phi and psi values for all points */
      Ntimesfreq = freq*Nharm_orig_d;
      for(i=0;i<N;i++)
	{
	  z_r[i] = cos(freq*t[i]);
	  z_i[i] = sin(freq*t[i]);
	  zn_r[i] = 1.;
	  zn_i[i] = 0.;
	  psi_r[i] = m_noave[i]*(cos(Ntimesfreq*t[i]));
	  psi_i[i] = m_noave[i]*(sin(Ntimesfreq*t[i]));
	  phi_r[i] = 1.;
	  phi_i[i] = 0.;
	}
      
      /* Now get the cn values using the recurrence algorithm */
      
      /* First get the n = 0 term */
      var1_r = 0.; var1_i = 0.;
      var2_r = 0.; var2_i = 0.;
      var3_r = 0.; var3_i = 0.;
      var4_r = 0.; var4_i = 0.;
      for(i=0;i<N;i++)
	{
	  s1 = z_r[i]*phi_r[i];
	  s2 = z_i[i]*phi_i[i];
	  s3 = z_r[i]*phi_i[i];
	  s4 = z_i[i]*phi_r[i];
	  s5 = zn_r[i]*phi_r[i];
	  s6 = zn_i[i]*phi_i[i];
	  s7 = zn_r[i]*phi_i[i];
	  s8 = zn_i[i]*phi_r[i];
	  s9 = s1 - s2;
	  s10 = s4 + s3;
	  s11 = s5 + s6;
	  s12 = s8 - s7;
	  /*
	  var1_r += (double) weight[i]*(s9*phi_r[i] + s10*phi_i[i]);
	  var1_i += (double) weight[i]*(-s9*phi_i[i] + s10*phi_r[i]);
	  var2_r += (double) weight[i]*(s11*phi_r[i] + s12*phi_i[i]);
	  var2_i += (double) weight[i]*(-s11*phi_i[i] + s12*phi_r[i]);*/
	  var1_r += (double) weight[i]*s9;
	  var1_i += (double) weight[i]*s10;
	  var2_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
	  var3_r += (double) weight[i]*(psi_r[i]*phi_r[i] + psi_i[i]*phi_i[i]);
	  var3_i += (double) weight[i]*(-psi_r[i]*phi_i[i] + psi_i[i]*phi_r[i]);
	  var4_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
	}
      var1 = (double) sqrt((double) var4_r);
      c_r[0] = (double) var3_r / var1;
      c_i[0] = (double) var3_i / var1;
      /*var5 = var2_r*var2_r + var2_i*var2_i;

      alpha_n_r = (double) ((var1_r*var2_r + var1_i*var2_i)/var5);
      alpha_n_i = (double) ((-var1_r*var2_i + var1_i*var2_r)/var5);
      */
      /** Using revised expressions for alpha_n from
	  Schwarzenberg-Czerny 2012, IAU, 285, p 81; these avoid a singularity
	  that was present in the previous expression for even sampling **/
      alpha_n_r = (double) (var1_r / (double) var2_r);
      alpha_n_i = (double) (var1_i / (double) var2_r);
      
      /* Get the rest of the harmonics */
      for(n=1;n<=Nharm;n++)
	{
	  for(i=0;i<N;i++)
	    {
	      /* Get the new phi values */
	      s1 = zn_r[i]*phi_r[i];
	      s2 = zn_i[i]*phi_i[i];
	      s3 = zn_r[i]*phi_i[i];
	      s4 = zn_i[i]*phi_r[i];
	      s5 = s1 + s2;
	      s6 = s4 - s3;
	      tmp1 = z_r[i]*phi_r[i] - z_i[i]*phi_i[i] - (alpha_n_r*s5 - alpha_n_i*s6);
	      tmp2 = z_r[i]*phi_i[i] + z_i[i]*phi_r[i] - (alpha_n_r*s6 + alpha_n_i*s5);
	      phi_r[i] = tmp1;
	      phi_i[i] = tmp2;
	      /* Update zn */
	      tmp1 = zn_r[i]*z_r[i] - zn_i[i]*z_i[i];
	      tmp2 = zn_r[i]*z_i[i] + zn_i[i]*z_r[i];
	      zn_r[i] = tmp1;
	      zn_i[i] = tmp2;
	    }
	  /* Calculate the new alpha_n and c values */
	  var1_r = 0.; var1_i = 0.;
	  var2_r = 0.; var2_i = 0.;
	  var3_r = 0.; var3_i = 0.;
	  var4_r = 0.; var4_i = 0.;
	  for(i=0;i<N;i++)
	    {
	      s1 = z_r[i]*phi_r[i];
	      s2 = z_i[i]*phi_i[i];
	      s3 = z_r[i]*phi_i[i];
	      s4 = z_i[i]*phi_r[i];
	      s5 = zn_r[i]*phi_r[i];
	      s6 = zn_i[i]*phi_i[i];
	      s7 = zn_r[i]*phi_i[i];
	      s8 = zn_i[i]*phi_r[i];
	      s9 = s1 - s2;
	      s10 = s4 + s3;
	      s11 = s5 + s6;
	      s12 = s8 - s7;
	      /*var1_r += (double) weight[i]*(s9*phi_r[i] + s10*phi_i[i]);
	      var1_i += (double) weight[i]*(-s9*phi_i[i] + s10*phi_r[i]);
	      var2_r += (double) weight[i]*(s11*phi_r[i] + s12*phi_i[i]);
	      var2_i += (double) weight[i]*(-s11*phi_i[i] + s12*phi_r[i]);*/
	      var1_r += (double) weight[i]*s9;
	      var1_i += (double) weight[i]*s10;
	      var2_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
	      var3_r += (double) weight[i]*(psi_r[i]*phi_r[i] + psi_i[i]*phi_i[i]);
	      var3_i += (double) weight[i]*(-psi_r[i]*phi_i[i] + psi_i[i]*phi_r[i]);
	      var4_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
	    }
	  var1 = (double) sqrt((double) var4_r);
	  c_r[n] = (double) var3_r / var1;
	  c_i[n] = (double) var3_i / var1;
	  /*var5 = var2_r*var2_r + var2_i*var2_i;
	  alpha_n_r = (double) ((var1_r*var2_r + var1_i*var2_i)/var5);
	  alpha_n_i = (double) ((-var1_r*var2_i + var1_i*var2_r)/var5);*/
	  /** Using revised expressions for alpha_n from
	  Schwarzenberg-Czerny 2012, IAU, 285, p 81; these avoid a singularity
	  that was present in the previous expression for even sampling **/
	  alpha_n_r = (double) (var1_r / (double) var2_r);
	  alpha_n_i = (double) (var1_i / (double) var2_r);
	}
      /* Finally calculate the ANOVA statistic */
      var1 = 0;
      for(n=0;n<=Nharm;n++)
	var1 += (double) (c_r[n]*c_r[n] + c_i[n]*c_i[n]);
      periodogram = (th_coeff1 * (double) var1)/(th_coeff2 * (lcvariance - (double) var1));
      
      *Nharm_used = (Nharm / 2);
#ifdef PARALLEL
      if(c_r != NULL) free(c_r);
      if(c_i != NULL) free(c_i);
      if(psi_r != NULL) free(psi_r);
      if(psi_i != NULL) free(psi_i);
      if(z_r != NULL) free(z_r);
      if(z_i != NULL) free(z_i);
      if(phi_r != NULL) free(phi_r);
      if(phi_i != NULL) free(phi_i);
      if(zn_r != NULL) free(zn_r);
      if(zn_i != NULL) free(zn_i);
#endif
      return(periodogram);
    }
  else
    {
      /* Do the calculation for optimized Nharm */
      Nharmtrial = 2;
      periodogram_new = 0.;
      do
	{
	  th_coeff1 = (double) (N - Nharmtrial - 2);
	  th_coeff2 = (double) (Nharmtrial + 1);
	  periodogram = periodogram_new;
	  Ntimesfreq = freq*((double) Nharmtrial / 2);
	  for(i=0;i<N;i++)
	    {
	      z_r[i] = cos(freq*t[i]);
	      z_i[i] = sin(freq*t[i]);
	      zn_r[i] = 1.;
	      zn_i[i] = 0.;
	      psi_r[i] = m_noave[i]*(cos(Ntimesfreq*t[i]));
	      psi_i[i] = m_noave[i]*(sin(Ntimesfreq*t[i]));
	      phi_r[i] = 1.;
	      phi_i[i] = 0.;
	    }
	  
	  /* Now get the cn values using the recurrence algorithm */
	  
	  /* First get the n = 0 term */
	  var1_r = 0.; var1_i = 0.;
	  var2_r = 0.; var2_i = 0.;
	  var3_r = 0.; var3_i = 0.;
	  var4_r = 0.; var4_i = 0.;
	  for(i=0;i<N;i++)
	    {
	      s1 = z_r[i]*phi_r[i];
	      s2 = z_i[i]*phi_i[i];
	      s3 = z_r[i]*phi_i[i];
	      s4 = z_i[i]*phi_r[i];
	      s5 = zn_r[i]*phi_r[i];
	      s6 = zn_i[i]*phi_i[i];
	      s7 = zn_r[i]*phi_i[i];
	      s8 = zn_i[i]*phi_r[i];
	      s9 = s1 - s2;
	      s10 = s4 + s3;
	      s11 = s5 + s6;
	      s12 = s8 - s7;
	      /*var1_r += (double) weight[i]*(s9*phi_r[i] + s10*phi_i[i]);
	      var1_i += (double) weight[i]*(-s9*phi_i[i] + s10*phi_r[i]);
	      var2_r += (double) weight[i]*(s11*phi_r[i] + s12*phi_i[i]);
	      var2_i += (double) weight[i]*(-s11*phi_i[i] + s12*phi_r[i]);*/
	      var1_r += (double) weight[i]*s9;
	      var1_i += (double) weight[i]*s10;
	      var2_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
	      var3_r += (double) weight[i]*(psi_r[i]*phi_r[i] + psi_i[i]*phi_i[i]);
	      var3_i += (double) weight[i]*(-psi_r[i]*phi_i[i] + psi_i[i]*phi_r[i]);
	      var4_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
	    }
	  var1 = (double) sqrt((double) var4_r);
	  c_r[0] = (double) var3_r / var1;
	  c_i[0] = (double) var3_i / var1;
	  /*var5 = var2_r*var2_r + var2_i*var2_i;
	  alpha_n_r = (double) ((var1_r*var2_r + var1_i*var2_i)/var5);
	  alpha_n_i = (double) ((-var1_r*var2_i + var1_i*var2_r)/var5);*/
	  /** Using revised expressions for alpha_n from
	  Schwarzenberg-Czerny 2012, IAU, 285, p 81; these avoid a singularity
	  that was present in the previous expression for even sampling **/
	  alpha_n_r = (double) (var1_r / (double) var2_r);
	  alpha_n_i = (double) (var1_i / (double) var2_r);
	  
	  /* Get the rest of the harmonics */
	  for(n=1;n<=Nharmtrial;n++)
	    {
	      for(i=0;i<N;i++)
		{
		  /* Get the new phi values */
		  s1 = zn_r[i]*phi_r[i];
		  s2 = zn_i[i]*phi_i[i];
		  s3 = zn_r[i]*phi_i[i];
		  s4 = zn_i[i]*phi_r[i];
		  s5 = s1 + s2;
		  s6 = s4 - s3;
		  tmp1 = z_r[i]*phi_r[i] - z_i[i]*phi_i[i] - (alpha_n_r*s5 - alpha_n_i*s6);
		  tmp2 = z_r[i]*phi_i[i] + z_i[i]*phi_r[i] - (alpha_n_r*s6 + alpha_n_i*s5);
		  phi_r[i] = tmp1;
		  phi_i[i] = tmp2;
		  /* Update zn */
		  tmp1 = zn_r[i]*z_r[i] - zn_i[i]*z_i[i];
		  tmp2 = zn_r[i]*z_i[i] + zn_i[i]*z_r[i];
		  zn_r[i] = tmp1;
		  zn_i[i] = tmp2;
		}
	      /* Calculate the new alpha_n and c values */
	      var1_r = 0.; var1_i = 0.;
	      var2_r = 0.; var2_i = 0.;
	      var3_r = 0.; var3_i = 0.;
	      var4_r = 0.; var4_i = 0.;
	      for(i=0;i<N;i++)
		{
		  s1 = z_r[i]*phi_r[i];
		  s2 = z_i[i]*phi_i[i];
		  s3 = z_r[i]*phi_i[i];
		  s4 = z_i[i]*phi_r[i];
		  s5 = zn_r[i]*phi_r[i];
		  s6 = zn_i[i]*phi_i[i];
		  s7 = zn_r[i]*phi_i[i];
		  s8 = zn_i[i]*phi_r[i];
		  s9 = s1 - s2;
		  s10 = s4 + s3;
		  s11 = s5 + s6;
		  s12 = s8 - s7;
		  /*var1_r += (double) weight[i]*(s9*phi_r[i] + s10*phi_i[i]);
		  var1_i += (double) weight[i]*(-s9*phi_i[i] + s10*phi_r[i]);
		  var2_r += (double) weight[i]*(s11*phi_r[i] + s12*phi_i[i]);
		  var2_i += (double) weight[i]*(-s11*phi_i[i] + s12*phi_r[i]);*/
		  var1_r += (double) weight[i]*s9;
		  var1_i += (double) weight[i]*s10;
		  var2_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
		  var3_r += (double) weight[i]*(psi_r[i]*phi_r[i] + psi_i[i]*phi_i[i]);
		  var3_i += (double) weight[i]*(-psi_r[i]*phi_i[i] + psi_i[i]*phi_r[i]);
		  var4_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
		}
	      var1 = (double) sqrt((double) var4_r);
	      c_r[n] = (double) var3_r / var1;
	      c_i[n] = (double) var3_i / var1;
	      /*var5 = var2_r*var2_r + var2_i*var2_i;
	      alpha_n_r = (double) ((var1_r*var2_r + var1_i*var2_i)/var5);
	      alpha_n_i = (double) ((-var1_r*var2_i + var1_i*var2_r)/var5);*/
	      /** Using revised expressions for alpha_n from
	      Schwarzenberg-Czerny 2012, IAU, 285, p 81; these avoid a singularity
	      that was present in the previous expression for even sampling **/
	      alpha_n_r = (double) (var1_r / (double) var2_r);
	      alpha_n_i = (double) (var1_i / (double) var2_r);
	    }
	  /* Finally calculate the ANOVA statistic */
	  var1 = 0;
	  for(n=0;n<=Nharmtrial;n++)
	    var1 += (double) (c_r[n]*c_r[n] + c_i[n]*c_i[n]);
	  periodogram_new = (th_coeff1 * (double) var1)/(th_coeff2 * (lcvariance - (double) var1));
	  a = (double) th_coeff2 * 0.5;
	  b = th_coeff1 * 0.5;
	  periodogram_new = -log1minusbetai(a, b, ((2. * a * periodogram_new)/(2. * (b + a*periodogram_new))));
	  //periodogram_new = betai((2. * ((double) Nharmtrial)), (((double) N - 2.*((double) Nharmtrial) - 1.)), periodogram_new);
	  Nharmtrial += 2;
	} while ((Nharmtrial < 6 || periodogram_new > periodogram) && Nharmtrial < MAX_NHARM_AUTO);
      
      *Nharm_used = ((Nharmtrial - 2)/2);
#ifdef PARALLEL
      if(c_r != NULL) free(c_r);
      if(c_i != NULL) free(c_i);
      if(psi_r != NULL) free(psi_r);
      if(psi_i != NULL) free(psi_i);
      if(z_r != NULL) free(z_r);
      if(z_i != NULL) free(z_i);
      if(phi_r != NULL) free(phi_r);
      if(phi_i != NULL) free(phi_i);
      if(zn_r != NULL) free(zn_r);
      if(zn_i != NULL) free(zn_i);
#endif
      return(periodogram);
    }
}
			   

void aov_harm(int N, double *t, double *m, double *sig, int Nharm, int Nfreq, double freqmin, double freqstep, double *periodogram, double **out_m_noave, double **out_t_nostart, double **out_weight, double *out_lcvariance)
{
  /* This is an implementation of the multi-harmonic AoV periodogram (Schwarzenberg-Czerny, 1996, ApJ, 460, L107.

The input values are:
N - the number of data points
t - vector of times
m - vector of magnitudes
m_noave - a storage vector of size N
sig - vector of magnitude uncertainties
Nharm - number of harmonics to include (if it is negative, then it will be varied for each frequency until the false alarm probability at that frequency is minimized).
Nfreq - number of frequencies to compute the periodogram for
freqmin - minimum frequency to begin computation at
freqstep - the size of the frequency step
periodogram - output periodogram vector (The AoV theta statistic if input Nharm is positive, the false alarm probability if input Nharm is negative)

out_m_noave, out_t_nostart, out_weight, out_lcvariance - these are pointers to variables allocated by this routine. They are output here and then passed to TestPeriod_aov_harm above later so that they don't have to be recomputed.

c_r, c_i are storage vectors of size 2*Nharm
psi_r and psi_i are storage vectors of size N
z_r and z_i are storage vectors of size N
phi_r and phi_i are storage vectors of size N

Note that memory for periodogram should be allocated before calling this function.
  */
  int i, n, fnum;
  double freq, Nharm_orig_d, alpha_n_r, alpha_n_i;
  double var1, var2, var1_r, var1_i, var2_r, var2_i, var3_r, var3_i, var4_r, var4_i, var5;
  double lcave, lcvariance, Ntimesfreq, tmp1, tmp2, th_coeff1, th_coeff2, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, periodogram_new, a, b;

  int Nharmtrial;
#ifdef PARALLEL
  int sizeNvecs = 0, sizeNharmvecs = 0;
  double *m_noave = NULL, *t_nostart = NULL, *weight = NULL, 
    *c_r = NULL, *c_i = NULL, *psi_r = NULL, *psi_i = NULL, *z_r = NULL, 
    *z_i = NULL, *phi_r = NULL, *phi_i = NULL, *zn_r = NULL, *zn_i = NULL;
#else
  static int sizeNvecs = 0, sizeNharmvecs = 0;
  static double *m_noave, *t_nostart, *weight, *c_r, *c_i, *psi_r, *psi_i, *z_r, *z_i, *phi_r, *phi_i, *zn_r, *zn_i;
#endif

  /* change the frequency to angular frequency for the computation, also double Nharm to account for the complex polynomials being of double the order of the real space polynomials */
  freqmin = TWOPID * freqmin;
  freqstep = TWOPID * freqstep;
  Nharm_orig_d = (double) Nharm;
  Nharm = Nharm * 2;

  th_coeff1 = (double) (N - Nharm - 1);
  th_coeff2 = (double) Nharm;

  /* Initialize memory for the storage vectors if this is our first pass, otherwise if the old vectors are not long enough allocation additional memory */
  if(!sizeNvecs)
    {
      sizeNvecs = N;
      if((m_noave = (double *) malloc(N * sizeof(double))) == NULL ||
	 (t_nostart = (double *) malloc(N * sizeof(double))) == NULL ||
	 (weight = (double *) malloc(N * sizeof(double))) == NULL ||
	 (psi_r = (double *) malloc(N * sizeof(double))) == NULL ||
	 (psi_i = (double *) malloc(N * sizeof(double))) == NULL ||
	 (z_r = (double *) malloc(N * sizeof(double))) == NULL ||
	 (z_i = (double *) malloc(N * sizeof(double))) == NULL ||
	 (zn_r = (double *) malloc(N * sizeof(double))) == NULL ||
	 (zn_i = (double *) malloc(N * sizeof(double))) == NULL ||
	 (phi_r = (double *) malloc(N * sizeof(double))) == NULL ||
	 (phi_i = (double *) malloc(N * sizeof(double))) == NULL)
	//error(ERR_MEMALLOC);
	exit(3);
    }
  else if(sizeNvecs < N)
    {
      sizeNvecs = N;
      if((m_noave = (double *) realloc(m_noave, N * sizeof(double))) == NULL ||
	 (t_nostart = (double *) realloc(t_nostart, N * sizeof(double))) == NULL ||
	 (weight = (double *) realloc(weight, N * sizeof(double))) == NULL ||
	 (psi_r = (double *) realloc(psi_r, N * sizeof(double))) == NULL ||
	 (psi_i = (double *) realloc(psi_i, N * sizeof(double))) == NULL ||
	 (z_r = (double *) realloc(z_r, N * sizeof(double))) == NULL ||
	 (z_i = (double *) realloc(z_i, N * sizeof(double))) == NULL ||
	 (zn_r = (double *) realloc(zn_r, N * sizeof(double))) == NULL ||
	 (zn_i = (double *) realloc(zn_i, N * sizeof(double))) == NULL ||
	 (phi_r = (double *) realloc(phi_r, N * sizeof(double))) == NULL ||
	 (phi_i = (double *) realloc(phi_i, N * sizeof(double))) == NULL)
	//error(ERR_MEMALLOC);
	exit(3);
    }
  if(!sizeNharmvecs)
    {
      if(Nharm > 0)
	{
	  sizeNharmvecs = Nharm + 1;
	  if((c_r = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL ||
	     (c_i = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL)
	    //error(ERR_MEMALLOC);
	    exit(3);
	}
      else
	{
	  sizeNharmvecs = MAX_NHARM_AUTO + 1;
	  if((c_r = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL ||
	     (c_i = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL)
	    exit(3);
	}
    }
  else if(sizeNharmvecs < Nharm + 1)
    {
      sizeNharmvecs = Nharm + 1;
      if((c_r = (double *) realloc(c_r, sizeNharmvecs * sizeof(double))) == NULL ||
	 (c_i = (double *) realloc(c_i, sizeNharmvecs * sizeof(double))) == NULL)
	//error(ERR_MEMALLOC);
	exit(3);
    }

  *out_m_noave = m_noave;
  *out_t_nostart = t_nostart;
  *out_weight = weight;

  /* subtract the average from the light curve and the zero time from the light curve, finally compute the light curve variance */
  var1 = 0.;
  var2 = 0.;
  for(i=0;i<N;i++)
    {
      weight[i] = 1./sig[i]/sig[i];
    }
  for(i=0;i<N;i++)
    {
      t_nostart[i] = t[i] - t[0];
      var1 += (double) m[i]*weight[i];
      var2 += (double) weight[i];
    }
  lcave = (double) (var1 / var2);
  for(i=0;i<N;i++)
    {
      m_noave[i] = m[i] - lcave;
    }
  lcvariance = 0.;
  for(i=0;i<N;i++)
    {
      lcvariance += weight[i]*m_noave[i]*m_noave[i];
    }

  *out_lcvariance = lcvariance;
      
  if(Nharm > 0)
    {
      /* Do the calculation for fixed Nharm */
      /* Step through frequency space, for each frequency fit the harmonic series (set of periodic orthogonal polynomials) using the fast recurrence algorithm. Then compute the AoV statistic for that fit. */
      
      for(fnum=0, freq = (freqmin + freqstep*((double)(Nfreq - 1))); fnum < Nfreq; fnum++, freq -= freqstep)
	{
	  /* For the trial frequency calculate the z, phi and psi values for all points */
	  Ntimesfreq = freq*Nharm_orig_d;
	  for(i=0;i<N;i++)
	    {
	      z_r[i] = cos(freq*t[i]);
	      z_i[i] = sin(freq*t[i]);
	      zn_r[i] = 1.;
	      zn_i[i] = 0.;
	      psi_r[i] = m_noave[i]*(cos(Ntimesfreq*t[i]));
	      psi_i[i] = m_noave[i]*(sin(Ntimesfreq*t[i]));
	      phi_r[i] = 1.;
	      phi_i[i] = 0.;
	    }
	  
	  /* Now get the cn values using the recurrence algorithm */
	  
	  /* First get the n = 0 term */
	  var1_r = 0.; var1_i = 0.;
	  var2_r = 0.; var2_i = 0.;
	  var3_r = 0.; var3_i = 0.;
	  var4_r = 0.; var4_i = 0.;
	  for(i=0;i<N;i++)
	    {
	      s1 = z_r[i]*phi_r[i];
	      s2 = z_i[i]*phi_i[i];
	      s3 = z_r[i]*phi_i[i];
	      s4 = z_i[i]*phi_r[i];
	      s5 = zn_r[i]*phi_r[i];
	      s6 = zn_i[i]*phi_i[i];
	      s7 = zn_r[i]*phi_i[i];
	      s8 = zn_i[i]*phi_r[i];
	      s9 = s1 - s2;
	      s10 = s4 + s3;
	      s11 = s5 + s6;
	      s12 = s8 - s7;
	      /*var1_r += (double) weight[i]*(s9*phi_r[i] + s10*phi_i[i]);
	      var1_i += (double) weight[i]*(-s9*phi_i[i] + s10*phi_r[i]);
	      var2_r += (double) weight[i]*(s11*phi_r[i] + s12*phi_i[i]);
	      var2_i += (double) weight[i]*(-s11*phi_i[i] + s12*phi_r[i]);*/
	      var1_r += (double) weight[i]*s9;
	      var1_i += (double) weight[i]*s10;
	      var2_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
	      var3_r += (double) weight[i]*(psi_r[i]*phi_r[i] + psi_i[i]*phi_i[i]);
	      var3_i += (double) weight[i]*(-psi_r[i]*phi_i[i] + psi_i[i]*phi_r[i]);
	      var4_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
	    }
	  var1 = (double) sqrt((double) var4_r);
	  c_r[0] = (double) var3_r / var1;
	  c_i[0] = (double) var3_i / var1;
	  /*var5 = var2_r*var2_r + var2_i*var2_i;
	  alpha_n_r = (double) ((var1_r*var2_r + var1_i*var2_i)/var5);
	  alpha_n_i = (double) ((-var1_r*var2_i + var1_i*var2_r)/var5);*/
	  /** Using revised expressions for alpha_n from
	      Schwarzenberg-Czerny 2012, IAU, 285, p 81; these avoid a singularity
	      that was present in the previous expression for even sampling **/
	  alpha_n_r = (double) (var1_r / (double) var2_r);
	  alpha_n_i = (double) (var1_i / (double) var2_r);
	  
	  /* Get the rest of the harmonics */
	  for(n=1;n<=Nharm;n++)
	    {
	      for(i=0;i<N;i++)
		{
		  /* Get the new phi values */
		  s1 = zn_r[i]*phi_r[i];
		  s2 = zn_i[i]*phi_i[i];
		  s3 = zn_r[i]*phi_i[i];
		  s4 = zn_i[i]*phi_r[i];
		  s5 = s1 + s2;
		  s6 = s4 - s3;
		  tmp1 = z_r[i]*phi_r[i] - z_i[i]*phi_i[i] - (alpha_n_r*s5 - alpha_n_i*s6);
		  tmp2 = z_r[i]*phi_i[i] + z_i[i]*phi_r[i] - (alpha_n_r*s6 + alpha_n_i*s5);
		  phi_r[i] = tmp1;
		  phi_i[i] = tmp2;
		  /* Update zn */
		  tmp1 = zn_r[i]*z_r[i] - zn_i[i]*z_i[i];
		  tmp2 = zn_r[i]*z_i[i] + zn_i[i]*z_r[i];
		  zn_r[i] = tmp1;
		  zn_i[i] = tmp2;
		}
	      /* Calculate the new alpha_n and c values */
	      var1_r = 0.; var1_i = 0.;
	      var2_r = 0.; var2_i = 0.;
	      var3_r = 0.; var3_i = 0.;
	      var4_r = 0.; var4_i = 0.;
	      for(i=0;i<N;i++)
		{
		  s1 = z_r[i]*phi_r[i];
		  s2 = z_i[i]*phi_i[i];
		  s3 = z_r[i]*phi_i[i];
		  s4 = z_i[i]*phi_r[i];
		  s5 = zn_r[i]*phi_r[i];
		  s6 = zn_i[i]*phi_i[i];
		  s7 = zn_r[i]*phi_i[i];
		  s8 = zn_i[i]*phi_r[i];
		  s9 = s1 - s2;
		  s10 = s4 + s3;
		  s11 = s5 + s6;
		  s12 = s8 - s7;
		  /*var1_r += (double) weight[i]*(s9*phi_r[i] + s10*phi_i[i]);
		  var1_i += (double) weight[i]*(-s9*phi_i[i] + s10*phi_r[i]);
		  var2_r += (double) weight[i]*(s11*phi_r[i] + s12*phi_i[i]);
		  var2_i += (double) weight[i]*(-s11*phi_i[i] + s12*phi_r[i]);*/
		  var1_r += (double) weight[i]*s9;
		  var1_i += (double) weight[i]*s10;
		  var2_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
		  var3_r += (double) weight[i]*(psi_r[i]*phi_r[i] + psi_i[i]*phi_i[i]);
		  var3_i += (double) weight[i]*(-psi_r[i]*phi_i[i] + psi_i[i]*phi_r[i]);
		  var4_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
		}
	      var1 = (double) sqrt((double) var4_r);
	      c_r[n] = (double) var3_r / var1;
	      c_i[n] = (double) var3_i / var1;
	      /*var5 = var2_r*var2_r + var2_i*var2_i;
	      alpha_n_r = (double) ((var1_r*var2_r + var1_i*var2_i)/var5);
	      alpha_n_i = (double) ((-var1_r*var2_i + var1_i*var2_r)/var5);*/
	      /** Using revised expressions for alpha_n from
	      Schwarzenberg-Czerny 2012, IAU, 285, p 81; these avoid a singularity
	      that was present in the previous expression for even sampling **/
	      alpha_n_r = (double) (var1_r / (double) var2_r);
	      alpha_n_i = (double) (var1_i / (double) var2_r);
	    }
	  /* Finally calculate the ANOVA statistic */
	  var1 = 0;
	  for(n=0;n<=Nharm;n++)
	    var1 += (double) (c_r[n]*c_r[n] + c_i[n]*c_i[n]);
	  periodogram[fnum] = (th_coeff1 * (double) var1)/(th_coeff2 * (lcvariance - (double) var1));
	}
    }
  else
    {
      /* Do the calculation for optimal Nharm */
      for(fnum=0, freq = (freqmin + freqstep*((double)(Nfreq - 1))); fnum < Nfreq; fnum++, freq -= freqstep)
	{
	  Nharmtrial = 2;
	  periodogram_new = 0.;
	  do
	    {
	      th_coeff1 = (double) (N - Nharmtrial - 2);
	      th_coeff2 = (double) (Nharmtrial + 1);
	      periodogram[fnum] = periodogram_new;
	      /* For the trial frequency calculate the z, phi and psi values for all points */
	      Ntimesfreq = freq*((double) Nharmtrial / 2);
	      for(i=0;i<N;i++)
		{
		  z_r[i] = cos(freq*t[i]);
		  z_i[i] = sin(freq*t[i]);
		  zn_r[i] = 1.;
		  zn_i[i] = 0.;
		  psi_r[i] = m_noave[i]*(cos(Ntimesfreq*t[i]));
		  psi_i[i] = m_noave[i]*(sin(Ntimesfreq*t[i]));
		  phi_r[i] = 1.;
		  phi_i[i] = 0.;
		}
	      
	      /* Now get the cn values using the recurrence algorithm */
	      
	      /* First get the n = 0 term */
	      var1_r = 0.; var1_i = 0.;
	      var2_r = 0.; var2_i = 0.;
	      var3_r = 0.; var3_i = 0.;
	      var4_r = 0.; var4_i = 0.;
	      for(i=0;i<N;i++)
		{
		  s1 = z_r[i]*phi_r[i];
		  s2 = z_i[i]*phi_i[i];
		  s3 = z_r[i]*phi_i[i];
		  s4 = z_i[i]*phi_r[i];
		  s5 = zn_r[i]*phi_r[i];
		  s6 = zn_i[i]*phi_i[i];
		  s7 = zn_r[i]*phi_i[i];
		  s8 = zn_i[i]*phi_r[i];
		  s9 = s1 - s2;
		  s10 = s4 + s3;
		  s11 = s5 + s6;
		  s12 = s8 - s7;
		  /*var1_r += (double) weight[i]*(s9*phi_r[i] + s10*phi_i[i]);
		  var1_i += (double) weight[i]*(-s9*phi_i[i] + s10*phi_r[i]);
		  var2_r += (double) weight[i]*(s11*phi_r[i] + s12*phi_i[i]);
		  var2_i += (double) weight[i]*(-s11*phi_i[i] + s12*phi_r[i]);*/
		  var1_r += (double) weight[i]*s9;
		  var1_i += (double) weight[i]*s10;
		  var2_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
		  var3_r += (double) weight[i]*(psi_r[i]*phi_r[i] + psi_i[i]*phi_i[i]);
		  var3_i += (double) weight[i]*(-psi_r[i]*phi_i[i] + psi_i[i]*phi_r[i]);
		  var4_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
		}
	      var1 = (double) sqrt((double) var4_r);
	      c_r[0] = (double) var3_r / var1;
	      c_i[0] = (double) var3_i / var1;
	      /*var5 = var2_r*var2_r + var2_i*var2_i;
	      alpha_n_r = (double) ((var1_r*var2_r + var1_i*var2_i)/var5);
	      alpha_n_i = (double) ((-var1_r*var2_i + var1_i*var2_r)/var5);*/
	      /** Using revised expressions for alpha_n from
	      Schwarzenberg-Czerny 2012, IAU, 285, p 81; these avoid a singularity
	      that was present in the previous expression for even sampling **/
	      alpha_n_r = (double) (var1_r / (double) var2_r);
	      alpha_n_i = (double) (var1_i / (double) var2_r);
	      
	      /* Get the rest of the harmonics */
	      for(n=1;n<=Nharmtrial;n++)
		{
		  for(i=0;i<N;i++)
		    {
		      /* Get the new phi values */
		      s1 = zn_r[i]*phi_r[i];
		      s2 = zn_i[i]*phi_i[i];
		      s3 = zn_r[i]*phi_i[i];
		      s4 = zn_i[i]*phi_r[i];
		      s5 = s1 + s2;
		      s6 = s4 - s3;
		      tmp1 = z_r[i]*phi_r[i] - z_i[i]*phi_i[i] - (alpha_n_r*s5 - alpha_n_i*s6);
		      tmp2 = z_r[i]*phi_i[i] + z_i[i]*phi_r[i] - (alpha_n_r*s6 + alpha_n_i*s5);
		      phi_r[i] = tmp1;
		      phi_i[i] = tmp2;
		      /* Update zn */
		      tmp1 = zn_r[i]*z_r[i] - zn_i[i]*z_i[i];
		      tmp2 = zn_r[i]*z_i[i] + zn_i[i]*z_r[i];
		      zn_r[i] = tmp1;
		      zn_i[i] = tmp2;
		    }
		  /* Calculate the new alpha_n and c values */
		  var1_r = 0.; var1_i = 0.;
		  var2_r = 0.; var2_i = 0.;
		  var3_r = 0.; var3_i = 0.;
		  var4_r = 0.; var4_i = 0.;
		  for(i=0;i<N;i++)
		    {
		      s1 = z_r[i]*phi_r[i];
		      s2 = z_i[i]*phi_i[i];
		      s3 = z_r[i]*phi_i[i];
		      s4 = z_i[i]*phi_r[i];
		      s5 = zn_r[i]*phi_r[i];
		      s6 = zn_i[i]*phi_i[i];
		      s7 = zn_r[i]*phi_i[i];
		      s8 = zn_i[i]*phi_r[i];
		      s9 = s1 - s2;
		      s10 = s4 + s3;
		      s11 = s5 + s6;
		      s12 = s8 - s7;
		      /*var1_r += (double) weight[i]*(s9*phi_r[i] + s10*phi_i[i]);
		      var1_i += (double) weight[i]*(-s9*phi_i[i] + s10*phi_r[i]);
		      var2_r += (double) weight[i]*(s11*phi_r[i] + s12*phi_i[i]);
		      var2_i += (double) weight[i]*(-s11*phi_i[i] + s12*phi_r[i]);*/
		      var1_r += (double) weight[i]*s9;
		      var1_i += (double) weight[i]*s10;
		      var2_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
		      var3_r += (double) weight[i]*(psi_r[i]*phi_r[i] + psi_i[i]*phi_i[i]);
		      var3_i += (double) weight[i]*(-psi_r[i]*phi_i[i] + psi_i[i]*phi_r[i]);
		      var4_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
		    }
		  var1 = (double) sqrt((double) var4_r);
		  c_r[n] = (double) var3_r / var1;
		  c_i[n] = (double) var3_i / var1;
		  /*var5 = var2_r*var2_r + var2_i*var2_i;
		  alpha_n_r = (double) ((var1_r*var2_r + var1_i*var2_i)/var5);
		  alpha_n_i = (double) ((-var1_r*var2_i + var1_i*var2_r)/var5);*/
		  /** Using revised expressions for alpha_n from
		      Schwarzenberg-Czerny 2012, IAU, 285, p 81; these avoid a singularity
		      that was present in the previous expression for even sampling **/
		  alpha_n_r = (double) (var1_r / (double) var2_r);
		  alpha_n_i = (double) (var1_i / (double) var2_r);
		}
	      /* Finally calculate the ANOVA statistic */
	      var1 = 0;
	      for(n=0;n<=Nharmtrial;n++)
		var1 += (double) (c_r[n]*c_r[n] + c_i[n]*c_i[n]);
	      periodogram_new = (th_coeff1 * (double) var1)/(th_coeff2 * (lcvariance - (double) var1));
	      a = (double) th_coeff2 * 0.5;
	      b = th_coeff1 * 0.5;
	      periodogram_new = -log1minusbetai(a, b, ((2. * a * periodogram_new)/(2. * (b + a*periodogram_new))));
	      //periodogram_new = betai((2. * ((double) Nharmtrial)), (((double) N - 2.*((double) Nharmtrial) - 1.)), periodogram_new);
	      Nharmtrial += 2;
	    } while ((Nharmtrial < 6 || periodogram_new > periodogram[fnum]) && Nharmtrial < MAX_NHARM_AUTO);
	  periodogram[fnum] = periodogram[fnum];
	}
    }
#ifdef PARALLEL
  if(c_r != NULL) free(c_r);
  if(c_i != NULL) free(c_i);
  if(psi_r != NULL) free(psi_r);
  if(psi_i != NULL) free(psi_i);
  if(z_r != NULL) free(z_r);
  if(z_i != NULL) free(z_i);
  if(phi_r != NULL) free(phi_r);
  if(phi_i != NULL) free(phi_i);
  if(zn_r != NULL) free(zn_r);
  if(zn_i != NULL) free(zn_i);
#endif
}

/* Given a light curve, this function will compute an AOV Multiharmonic periodogram and find the top Npeaks peaks */
void findPeaks_aovharm(double *t_in, double *mag_in, double *sig_in, int N_in, double *perpeaks, double *aovpeaks, double *aovSNR, double *aovFAP, int *Nharm_used, int Npeaks, double minP, double maxP, double subsample, double fine_tune, int outflag, char *outname, double *aveaov, double *stddevaov, double *aveaov_whiten, double *stddevaov_whiten, int ascii, int Nharm, int whiten, double clip, int clipiter, int fixperiodSNR, double fixperiodSNR_period, double *fixperiodSNR_value, double *fixperiodSNR_SNR, double *fixperiodSNR_FAP, int lcnum, int lclistnum, int usemask, _Variable *maskvar)
{
  int i, j, k, foundsofar, test, Nperiod, a, b, abest, bbest, ismultiple, Ngood,l1, l2, peakiter, nclippedthis, nclippedlast, testNharm, m_eff;

  double testperiod, bestscore;
  double T, freq, minfreq, minbest, freqstep, smallfreqstep, score, aveper, stdper, *m_noave = NULL, *t_nostart = NULL, *weight = NULL, lcvariance, lastpoint, negln_m_eff;

  long double Sum, Sumsqr;

  FILE *outfile;

  char outstring[MAXLEN];

  double *mag_cpy, *sig_cpy, *t_cpy, a_, b_;

  int bestNharmused;

#ifdef PARALLEL
  int sizeperiodvec = 0;
  int sizeperiodvec_whiten = 0;
  double *periods = NULL, *periodogram = NULL, **periodogram_whiten = NULL;

  double *aveper_whiten = NULL, *stdper_whiten = NULL, *harmA = NULL, *harmB = NULL;
  int size_aveper_whiten = 0;
#else
  static int sizeperiodvec = 0;
  static int sizeperiodvec_whiten = 0;
  static double *periods, *periodogram, **periodogram_whiten;

  static double *aveper_whiten, *stdper_whiten, *harmA, *harmB;
  static int size_aveper_whiten = 0;
#endif

  double fundA, fundB, meanval, amp;

  int N, N_mask;
  double *t_mask = NULL, *mag_mask = NULL, *sig_mask = NULL;
  double *t, *mag, *sig;

  if(N_in <= 1) {
    for(i = 0; i < Npeaks; i++) {
      perpeaks[i] = -1.0;
      aovpeaks[i] = 0.0;
      aovSNR[i] = 0.0;
      aovFAP[i] = 0.0;
      Nharm_used[i] = -1;
    }
    *aveaov = 0.0;
    *stddevaov = 0.0;
    if(whiten) {
      for(i=0; i < Npeaks; i++) {
	aveaov_whiten[i] = 0.0;
	stddevaov_whiten[i] = 0.0;
      }
    }
    if(fixperiodSNR) {
      *fixperiodSNR_value = 0.0;
      *fixperiodSNR_SNR = 0.0;
      *fixperiodSNR_FAP = 0.0;
    }
    return;
  }

  if(!usemask) {
    N = N_in;
    t = t_in;
    mag = mag_in;
    sig = sig_in;
  } else {
    if((t_mask = (double *) malloc(N_in * sizeof(double))) == NULL ||
       (mag_mask = (double *) malloc(N_in * sizeof(double))) == NULL ||
       (sig_mask = (double *) malloc(N_in * sizeof(double))) == NULL) {
      error(ERR_MEMALLOC);
    }
    N = 0;
    for(i = 0; i < N_in; i++) {
      if(!isnan(mag_in[i]) && sig_in[i] > 0. && EvaluateVariable_Double(lclistnum, lcnum, i, maskvar) > VARTOOLS_MASK_TINY) {
	t_mask[N] = t_in[i];
	mag_mask[N] = mag_in[i];
	sig_mask[N] = sig_in[i];
	N++;
      }
    }
    t = t_mask;
    mag = mag_mask;
    sig = sig_mask;
    if(N <= 1) {
      for(i = 0; i < Npeaks; i++) {
	perpeaks[i] = -1.0;
	aovpeaks[i] = 0.0;
	aovSNR[i] = 0.0;
	aovFAP[i] = 0.0;
	Nharm_used[i] = -1;
      }
      *aveaov = 0.0;
      *stddevaov = 0.0;
      if(whiten) {
	for(i=0; i < Npeaks; i++) {
	  aveaov_whiten[i] = 0.0;
	  stddevaov_whiten[i] = 0.0;
	}
      }
      if(fixperiodSNR) {
	*fixperiodSNR_value = 0.0;
	*fixperiodSNR_SNR = 0.0;
	*fixperiodSNR_FAP = 0.0;
      }
      free(t_mask);
      free(mag_mask);
      free(sig_mask);
      return;
    }
  }
      

  if(!size_aveper_whiten)
    {
      aveper_whiten = (double *) malloc((Npeaks + 1) * sizeof(double));
      stdper_whiten = (double *) malloc((Npeaks + 1) * sizeof(double));
      if(Nharm > 0)
	{
	  harmA = (double *) malloc((Nharm + 1) * sizeof(double));
	  harmB = (double *) malloc((Nharm + 1) * sizeof(double));
	}
      else
	{
	  harmA = (double *) malloc((MAX_NHARM_AUTO + 1) * sizeof(double));
	  harmB = (double *) malloc((MAX_NHARM_AUTO + 1) * sizeof(double));
	}
      size_aveper_whiten = Npeaks + 1;
    }

  /* Find the number of periods we need to compute, and then allocate memory for the periodogram if we need to */

  /* initialize some of the period search variables*/
  T = t[N - 1] - t[0];

  // Initialize the periodogram if it hasn't already been initialized

  Nperiod = 0;
  freq = 1./minP;
  minfreq = 1./maxP;
  freqstep = subsample / T;
  while(freq >= minfreq)
    {
      freq -= freqstep;
      Nperiod++;
    }

  /* Set the estimate for the effective number of independent frequency samples */
  m_eff = (int) ((1./minP - 1./maxP) * T);
  if(subsample > 1.)
    m_eff = (int) ((double) m_eff / subsample);

  negln_m_eff = -log((double) m_eff);

  if(!sizeperiodvec && !whiten)
    {
      if(!sizeperiodvec && !sizeperiodvec_whiten)
	{
	  if((periods = (double *) malloc(Nperiod * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      else if(Nperiod > MAX_(sizeperiodvec_whiten,sizeperiodvec))
	{
	  if((periods = (double *) realloc(periods, Nperiod * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      sizeperiodvec = Nperiod;
      if((periodogram = (double *) malloc((Nperiod) * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
  else if(sizeperiodvec < Nperiod && !whiten)
    {
      if(Nperiod > MAX_(sizeperiodvec_whiten,sizeperiodvec))
	{
	  if((periods = (double *) realloc(periods, Nperiod * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}    
      sizeperiodvec = Nperiod;
      if((periodogram = (double *) realloc(periodogram, Nperiod * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
  if(!sizeperiodvec_whiten && whiten)
    {
      if(!sizeperiodvec && !sizeperiodvec_whiten)
	{
	  if((periods = (double *) malloc(Nperiod * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      else if(Nperiod > MAX_(sizeperiodvec_whiten,sizeperiodvec))
	{
	  if((periods = (double *) realloc(periods, Nperiod * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      sizeperiodvec_whiten = Nperiod;
      if((periodogram_whiten = (double **) malloc((Npeaks + 1) * sizeof(double *))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0;i<Npeaks+1;i++)
	{
	  if((periodogram_whiten[i] = (double *) malloc(Nperiod * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
    }
  else if(sizeperiodvec_whiten < Nperiod && whiten)
    {
      if(Nperiod > MAX_(sizeperiodvec_whiten,sizeperiodvec))
	{
	  if((periods = (double *) realloc(periods, Nperiod * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}    
      sizeperiodvec_whiten = Nperiod;
      for(i=0;i<Npeaks+1;i++)
	{
	  if((periodogram_whiten[i] = (double *) realloc(periodogram_whiten[i], Nperiod * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
    }
  freq = 1./minP;
  minfreq = 1./maxP;
  i = 0;
  while(freq >= minfreq)
    {
      periods[i] = 1./freq;
      freq -= freqstep;
      i++;
    }

  /* Here we do things differently depending on whether or not we're going to whiten the signal between finding peaks */
  if(!whiten)
    {
      /* Compute the periodogram */
#ifdef PARALLEL
      if(m_noave != NULL) { free(m_noave); m_noave = NULL;}
      if(t_nostart != NULL) { free(t_nostart); t_nostart = NULL;}
      if(weight != NULL) { free(weight); weight = NULL;}
#endif
      aov_harm(N, t, mag, sig, Nharm, Nperiod, minfreq, freqstep, periodogram, &m_noave, &t_nostart, &weight, &lcvariance);

      l1= (int) (log(maxP)/log(10.0));
      l1 += 2;
      if(l1 < 0)
	l1 = 1;
      l2= (int) (log(freqstep*minP*minP)/log(10.0));
      l2--;
      if(l2 < 0)
	l1 = l1 - l2 + 1;
      else
	l1 += 2;
      if(l2 < 0)
	sprintf(outstring,"%%%d.%df %%f\n",l1,-l2);
      else
	sprintf(outstring,"%%%d.1f %%f\n",l1);
      /* Write out the periodogram if asked to */
      if(outflag)
	{
	  if((outfile = fopen(outname,"w")) == NULL)
	    {
	      fprintf(stderr,"Cannot Write periodogram to %s\n",outname);
	      exit(3);
	    }
	  if(ascii)
	    {
	      fprintf(outfile,"#Period theta_AOV_SDE\n");
	      for(i=0;i<Nperiod;i++)
		{
		  fprintf(outfile,outstring,periods[i], periodogram[i]);
		}
	    }
	  else
	    {
	      fwrite(&Nperiod,4,1,outfile);
	      fwrite(&aveper,8,1,outfile);
	      fwrite(&stdper,8,1,outfile);
	      fwrite(periods,8,Nperiod,outfile);
	      fwrite(periodogram,8,Nperiod,outfile);
	    }
	  
	  fclose(outfile);
	}
      
      /* Compute the average and std-dev of the periodogram */
      Sum = 0.;
      Sumsqr = 0.;
      Ngood = 0;
      for(i=0;i<Nperiod;i++)
	{
	  if(periodogram[i] > ERROR_SCORE && !isnan(periodogram[i]) && periodogram[i]*0.0 == 0.0)
	    {
	      Sum += (long double) periodogram[i];
	      Sumsqr += (long double) (periodogram[i] * periodogram[i]);
	      Ngood++;
	    }
	}
      
      if(Ngood > 0) {
	Sum /= Ngood;
	Sumsqr /= Ngood;
	aveper = (double) Sum;
	stdper = sqrt((double)(Sumsqr - (Sum*Sum)));
      } else {
	aveper = 0.;
	stdper = 0.;
      }
      nclippedthis = 0;
      do {
	nclippedlast = nclippedthis;
	nclippedthis = 0;
	Sum = 0.;
	Sumsqr = 0.;
	Ngood = 0;
	for(i=0;i<Nperiod;i++)
	  {
	    if(periodogram[i] > ERROR_SCORE && !isnan(periodogram[i]) && periodogram[i]*0.0 == 0.0)
	      {
		if(periodogram[i] < aveper + clip*stdper)
		  {
		    Sum += (long double) periodogram[i];
		    Sumsqr += (long double) (periodogram[i]*periodogram[i]);
		    Ngood++;
		  }
		else
		  nclippedthis++;
	      }
	  }
	if(Ngood > 0) {
	  Sum /= Ngood;
	  Sumsqr /= Ngood;
	  aveper = (double) Sum;
	  stdper = sqrt((double)(Sumsqr - (Sum*Sum)));
	} else {
	  aveper = 0.;
	  stdper = 0.;
	  break;
	}
      } while(clipiter && nclippedthis > nclippedlast);
      
      *aveaov = aveper;
      *stddevaov = stdper;
      
      /* Get the peak value and SNR for the fixed period if we're doing that */
      if(fixperiodSNR)
	{
	  AOV_getfixedperiodSNR(periods,periodogram,Nperiod,aveper,stdper,fixperiodSNR_period,fixperiodSNR_value,fixperiodSNR_SNR);
	  if(Nharm > 0)
	    {
	      a_ = (double) Nharm;
	      b_ = 0.5*(((double) N - 2. *((double) Nharm) - 1.));
	      if(!isnan((*fixperiodSNR_value)) && (*fixperiodSNR_value > ERROR_SCORE)) {
		*fixperiodSNR_FAP = -log1minusbetai(a_, b_, ((2. * a_ * (*fixperiodSNR_value))/(2. * (b_ + a_*(*fixperiodSNR_value)))));
	      } else {
		*fixperiodSNR_FAP = 0.;
	      }
	    }
	}
      
      /* Replace the periodogram with only points that are local maxima */
      i = 0;
      while(i < Nperiod && !(periodogram[i] > ERROR_SCORE && !isnan(periodogram[i]) && periodogram[i]*0.0 == 0.0))
	i++;
      if(i < Nperiod) {
	lastpoint = periodogram[i] - 1.;
	for(k=i;k<Nperiod-1;k++)
	  {
	    if(periodogram[k] > ERROR_SCORE && !isnan(periodogram[k]) && periodogram[k]*0.0 == 0.0 && periodogram[k+1] > ERROR_SCORE && !isnan(periodogram[k+1]) && periodogram[k+1]*0.0 == 0.0) {
	      if(periodogram[k] > lastpoint && periodogram[k] > periodogram[k + 1])
		{
		  lastpoint = periodogram[k];
		  periodogram[i] = periodogram[k];
		  periods[i] = periods[k];
		  i++;
		}
	      else
		{
		  lastpoint = periodogram[k];
		}
	    }
	  }
	if(periodogram[k] > ERROR_SCORE && !isnan(periodogram[k]) && periodogram[k]*0.0 == 0.0) {
	  if(periodogram[k] > lastpoint)
	    {
	      periodogram[i] = periodogram[k];
	      periods[i] = periods[k];
	      i++;
	    }
	}
      }
      Nperiod = i;
      
      /* Search through the periodogram to identify the best Npeaks periods */
      
      foundsofar = 0;
      i = 0;
      while(foundsofar < Npeaks && i < Nperiod)
	{
	  test = 1;
	  for(j=0;j<foundsofar;j++)
	    {
	      if(!isDifferentPeriods(perpeaks[j],periods[i],T))
		{
		  if(periodogram[i] > aovpeaks[j])
		    {
		      perpeaks[j] = periods[i];
		      aovpeaks[j] = periodogram[i];
		    }
		  test = 0;
		  break;
		}
	    }
	  if(test)
	    {
	      perpeaks[foundsofar] = periods[i];
	      aovpeaks[foundsofar] = periodogram[i];
	      foundsofar++;
	    }
	  i++;
	}
      
      if(foundsofar < Npeaks)
	{
	  for(k=foundsofar;k<Npeaks;k++)
	    {
	      perpeaks[k] = ERROR_SCORE - 1.;
	      aovpeaks[k] = ERROR_SCORE - 1.;
	    }
	}
      
      mysort2(Npeaks,aovpeaks,perpeaks);
      
      minbest = aovpeaks[0];
      
      for(;i<Nperiod;i++)
	{
	  if(periodogram[i] > minbest)
	    {
	      test = 1;
	      for(j=0;j<Npeaks;j++)
		{
		  if(!isDifferentPeriods(periods[i],perpeaks[j],T))
		    {
		      if(periodogram[i] > aovpeaks[j])
			{
			  aovpeaks[j] = periodogram[i];
			  perpeaks[j] = periods[i];
			  mysort2(Npeaks,aovpeaks,perpeaks);
			  minbest = aovpeaks[0];
			}
		      test = 0;
		      break;
		    }
		}
	      if(test)
		{
		  perpeaks[0] = periods[i];
		  aovpeaks[0] = periodogram[i];
		  mysort2(Npeaks,aovpeaks,perpeaks);
		  minbest = aovpeaks[0];
		}
	    }
	}
      
      /* Now perform the high-resolution period scan on the peaks */
      
      smallfreqstep = fine_tune/T;
      
      for(j=0;j<Npeaks;j++)
	{
	  if(!isnan(perpeaks[j]) && perpeaks[j] > ERROR_SCORE)
	    {
	      freq = dmin((1./perpeaks[j]) + freqstep,(1./minP));
	      minfreq = dmax((1./perpeaks[j]) - freqstep,(1./maxP));
	      while (freq >= minfreq)
		{
		  testperiod = 1./freq;
		  score = TestPeriod_aov_harm(N, t, mag, sig, Nharm, testperiod, m_noave, t_nostart, weight, lcvariance, &testNharm);
		  if(score > aovpeaks[j])
		    {
		      aovpeaks[j] = score;
		      perpeaks[j] = testperiod;
		      Nharm_used[j] = testNharm;
		    }
		  freq -= smallfreqstep;
		}
	      // double-check period multiples
	      bestscore = aovpeaks[j];
	      ismultiple = 0;
	      for(a=1 ; a <= MAX_DOUBLE_CHECK_MULTIPLE ; a++)
		for(b=1 ; b <= MAX_DOUBLE_CHECK_MULTIPLE ; b++)
		  if(a != b)
		    {
		      testperiod = perpeaks[j] * a / b;
		      if((testperiod > minP) && (testperiod < maxP))
			{
			  score = TestPeriod_aov_harm(N, t, mag, sig, Nharm, testperiod, m_noave, t_nostart, weight, lcvariance, &testNharm);
			  if (score > bestscore)
			    {
			      ismultiple = 1;
			      abest = a;
			      bbest = b;
			      bestscore = score;
			      bestNharmused = testNharm;
			    }
			}
		    }
	      if(ismultiple)
		{
		  perpeaks[j] = perpeaks[j] * abest / bbest;
		  aovpeaks[j] = bestscore;
		  Nharm_used[j] = bestNharmused;
		}
	    }
	}
      
      /* Sort it so that the higher aov values come first in the vector */
      for(k=0;k<Npeaks;k++)
	aovpeaks[k] = -aovpeaks[k];
      mysort3_int(Npeaks,aovpeaks,perpeaks,Nharm_used);
      for(k=0;k<Npeaks;k++)
	{
	  aovpeaks[k] = -aovpeaks[k];
	  aovSNR[k] = (aovpeaks[k] - aveper)/stdper;
	  if(Nharm > 0)
	    {
	      a_ = (double) Nharm;
	      b_ = 0.5*(((double) N - 2. *((double) Nharm) - 1.));
	      aovFAP[k] = -log1minusbetai(a_, b_, ((2. * a_ * aovpeaks[k])/(2. * (b_ + a_*aovpeaks[k])))) + negln_m_eff;
	    }
	}

      /* Set ERROR period peaks to 1. so that other routines that take the period from aov don't have trouble with the negative period. */
      for(k=0;k<Npeaks;k++)
	{
	  if(perpeaks[k] < ERROR_SCORE || isnan(perpeaks[k])) {
	    perpeaks[k] = 1.;
	    aovSNR[k] = ERROR_SCORE-1;
	    aovFAP[k] = ERROR_SCORE-1;
	  }
	}
    }
  else
    {
      /* We need to whiten the light curve and re-compute the periodogram after finding each peak */

      /* Make a copy of the light curve */
      if((t_cpy = (double *) malloc(N * sizeof(double))) == NULL ||
	 (mag_cpy = (double *) malloc(N * sizeof(double))) == NULL ||
	 (sig_cpy = (double *) malloc(N * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      memcpy(t_cpy,t,N*sizeof(double));
      memcpy(mag_cpy,mag,N*sizeof(double));
      memcpy(sig_cpy,sig,N*sizeof(double));

      /* Compute the periodogram */
#ifdef PARALLEL
      if(m_noave != NULL) { free(m_noave); m_noave = NULL;}
      if(t_nostart != NULL) { free(t_nostart); t_nostart = NULL;}
      if(weight != NULL) { free(weight); weight = NULL;}
#endif
      aov_harm(N, t_cpy, mag_cpy, sig_cpy, Nharm, Nperiod, minfreq, freqstep, periodogram_whiten[0], &m_noave, &t_nostart, &weight, &lcvariance);
      
      for(peakiter=0;peakiter < Npeaks; peakiter++)
	{
	  /* Find the peak period */
	  perpeaks[peakiter] = 1.0;
	  aovpeaks[peakiter] = ERROR_SCORE-1;
	  for(i=0; i < Nperiod; i++)
	    {
	      if(periodogram_whiten[peakiter][i] > ERROR_SCORE && periodogram_whiten[peakiter][i]*0.0 == 0.0 && !isnan(periodogram_whiten[peakiter][i]))
		{
		  if(periodogram_whiten[peakiter][i] > aovpeaks[peakiter])
		    {
		      test = 1;
		      for(j=0;j<peakiter;j++)
			{
			  if(!isDifferentPeriods(periods[i],perpeaks[j],T))
			    {
			      test = 0;
			      break;
			    }
			}
		      if(test)
			{
			  perpeaks[peakiter] = periods[i];
			  aovpeaks[peakiter] = periodogram_whiten[peakiter][i];
			}
		    }
		}
	    }

	  /* Now perform the high-resolution period scan on the peak */
	  
	  smallfreqstep = fine_tune/T;
	  j = peakiter;	  

	  if(!isnan(aovpeaks[j]) && aovpeaks[j] > ERROR_SCORE && !isnan(perpeaks[j]))
	    {
	      freq = dmin((1./perpeaks[j]) + freqstep,(1./minP));
	      minfreq = dmax((1./perpeaks[j]) - freqstep,(1./maxP));
	      while (freq >= minfreq)
		{
		  testperiod = 1./freq;
		  score = TestPeriod_aov_harm(N, t_cpy, mag_cpy, sig_cpy, Nharm, testperiod, m_noave, t_nostart, weight, lcvariance, &testNharm);
		  if(score > aovpeaks[j])
		    {
		      aovpeaks[j] = score;
		      perpeaks[j] = testperiod;
		      Nharm_used[j] = testNharm;
		    }
		  freq -= smallfreqstep;
		}
	      // double-check period multiples
	      bestscore = aovpeaks[j];
	      ismultiple = 0;
	      for(a=1 ; a <= MAX_DOUBLE_CHECK_MULTIPLE ; a++)
		for(b=1 ; b <= MAX_DOUBLE_CHECK_MULTIPLE ; b++)
		  if(a != b)
		    {
		      testperiod = perpeaks[j] * a / b;
		      if((testperiod > minP) && (testperiod < maxP))
			{
			  score = TestPeriod_aov_harm(N, t_cpy, mag_cpy, sig_cpy, Nharm, testperiod, m_noave, t_nostart, weight, lcvariance, &testNharm);
			  if (score > bestscore)
			    {
			      ismultiple = 1;
			      abest = a;
			      bbest = b;
			      bestscore = score;
			      bestNharmused = testNharm;
			    }
			}
		    }
	      if(ismultiple)
		{
		  perpeaks[j] = perpeaks[j] * abest / bbest;
		  aovpeaks[j] = bestscore;
		  Nharm_used[j] = bestNharmused;
		}
	    }

	  if(!isnan(aovpeaks[peakiter]) && aovpeaks[peakiter] > ERROR_SCORE)
	    {
	      /* Whiten the light curve at this period */
	      if(Nharm > 0)
		{
		  dokillharms(N, t_cpy, mag_cpy, sig_cpy, 1, &(perpeaks[peakiter]), 0, (Nharm-1), NULL, NULL, &harmA, &harmB, &fundA, &fundB, &meanval, 0, NULL, &amp, 0, KILLHARM_OUTTYPE_DEFAULT, -1.);
		}
	      else
		{
		  dokillharms(N, t_cpy, mag_cpy, sig_cpy, 1, &(perpeaks[peakiter]), 0, (Nharm_used[peakiter]-1), NULL, NULL, &harmA, &harmB, &fundA, &fundB, &meanval, 0, NULL, &amp, 0, KILLHARM_OUTTYPE_DEFAULT, -1.);
		}
	      /* Compute the periodogram of the whitened light curve */
#ifdef PARALLEL
	      if(m_noave != NULL) { free(m_noave); m_noave = NULL;}
	      if(t_nostart != NULL) { free(t_nostart); t_nostart = NULL;}
	      if(weight != NULL) { free(weight); weight = NULL;}
#endif
	      aov_harm(N, t_cpy, mag_cpy, sig_cpy, Nharm, Nperiod, minfreq, freqstep, periodogram_whiten[peakiter+1], &m_noave, &t_nostart, &weight, &lcvariance);

	      /* Compute the average and std-dev of the periodogram */
	      Sum = 0.;
	      Sumsqr = 0.;
	      Ngood = 0;
	      for(i=0;i<Nperiod;i++)
		{
		  if(periodogram_whiten[peakiter+1][i] > ERROR_SCORE && !isnan(periodogram_whiten[peakiter+1][i]) && periodogram_whiten[peakiter+1][i]*0.0 == 0.0)
		    {
		      Sum += (long double) periodogram_whiten[peakiter+1][i];
		      Sumsqr += (long double) (periodogram_whiten[peakiter+1][i] * periodogram_whiten[peakiter+1][i]);
		      Ngood++;
		    }
		}
	      
	      Sum /= Ngood;
	      Sumsqr /= Ngood;
	      
	      aveper_whiten[peakiter] = (double) Sum;
	      stdper_whiten[peakiter] = sqrt((double)(Sumsqr - (Sum*Sum)));
	      
	      nclippedthis = 0;
	      do {
		nclippedlast = nclippedthis;
		nclippedthis = 0;
		Sum = 0.;
		Sumsqr = 0.;
		Ngood = 0;
		for(i=0;i<Nperiod;i++)
		  {
		    if(periodogram_whiten[peakiter+1][i] > ERROR_SCORE && !isnan(periodogram_whiten[peakiter+1][i]) && periodogram_whiten[peakiter+1][i]*0.0 == 0.0)
		      {
			if(periodogram_whiten[peakiter+1][i] < aveper_whiten[peakiter] + clip*stdper_whiten[peakiter])
			  {
			    Sum += (long double) periodogram_whiten[peakiter+1][i];
			    Sumsqr += (long double) (periodogram_whiten[peakiter+1][i]*periodogram_whiten[peakiter+1][i]);
			    Ngood++;
			  }
			else
			  nclippedthis++;
		      }
		  }
		if(Ngood > 0) {
		  Sum /= Ngood;
		  Sumsqr /= Ngood;
		  aveper_whiten[peakiter] = (double) Sum;
		  stdper_whiten[peakiter] = sqrt((double)(Sumsqr - (Sum*Sum)));
		} else {
		  aveper_whiten[peakiter] = ERROR_SCORE - 1;
		  stdper_whiten[peakiter] = ERROR_SCORE - 1;
		  break;
		}
	      } while(clipiter && nclippedthis > nclippedlast);
	      
	      if(Ngood > 0) {
		if(fixperiodSNR && !peakiter)
		  {
		    AOV_getfixedperiodSNR(periods,periodogram_whiten[peakiter],Nperiod,aveper_whiten[peakiter],stdper_whiten[peakiter],fixperiodSNR_period,fixperiodSNR_value,fixperiodSNR_SNR);
		    if(Nharm > 0)
		      {
			a_ = (double) Nharm;
			b_ = 0.5*(((double) N - 2. *((double) Nharm) - 1.));
			*fixperiodSNR_FAP = -log1minusbetai(a_, b_, ((2. * a_ * (*fixperiodSNR_value))/(2. * (b_ + a_*(*fixperiodSNR_value)))));
		      }
		  }
		aveaov_whiten[peakiter] = aveper_whiten[peakiter];
		stddevaov_whiten[peakiter] = stdper_whiten[peakiter];
		aovSNR[peakiter] = (aovpeaks[peakiter] - aveper_whiten[peakiter])/stdper_whiten[peakiter];
		if(Nharm > 0)
		  {
		    a_ = (double) Nharm;
		    b_ = 0.5*(((double) N - 2. *((double) Nharm) - 1.));
		    aovFAP[peakiter] = -log1minusbetai(a_, b_, ((2. * a_ * aovpeaks[peakiter])/(2. * (b_ + a_*aovpeaks[peakiter])))) + negln_m_eff;
		  }
	      } else {
		aovpeaks[peakiter] = ERROR_SCORE - 1;
		perpeaks[peakiter] = 1.;
		aovFAP[peakiter] = ERROR_SCORE - 1;
		aovSNR[peakiter] = ERROR_SCORE - 1;
		aveaov_whiten[peakiter] = 0.;
		stddevaov_whiten[peakiter] = 0.;
	      }
	    }
	  else
	    {
	      for(;peakiter < Npeaks; peakiter++)
		{
		  aovpeaks[peakiter] = ERROR_SCORE-1;
		  perpeaks[peakiter] = 1.;
		  aovSNR[peakiter] = ERROR_SCORE-1;
		  aovFAP[peakiter] = ERROR_SCORE-1;
		  aveaov_whiten[peakiter] = 0.;
		  stddevaov_whiten[peakiter] = 0.;
		}
	    }     
	}

      l1= (int) (log(maxP)/log(10.0));
      l1 += 2;
      if(l1 < 0)
	l1 = 1;
      l2= (int) (log(freqstep*minP*minP)/log(10.0));
      l2--;
      if(l2 < 0)
	l1 = l1 - l2 + 1;
      else
	l1 += 2;
      if(l2 < 0)
	sprintf(outstring,"%%%d.%df",l1,-l2);
      else
	sprintf(outstring,"%%%d.1f",l1);
      /* Write out the periodogram if asked to */
      if(outflag)
	{
	  if((outfile = fopen(outname,"w")) == NULL)
	    {
	      fprintf(stderr,"Cannot Write periodogram to %s\n",outname);
	      exit(3);
	    }
	  if(ascii)
	    {
	      fprintf(outfile,"#Period theta_AOV_SDE\n");
	      for(i=0;i<Nperiod;i++)
		{
		  fprintf(outfile,outstring,periods[i]);
		  for(j=0;j<Npeaks;j++)
		    fprintf(outfile," %f",periodogram_whiten[j][i]);
		  fprintf(outfile,"\n");
		}
	    }
	  else
	    {
	      fwrite(&Nperiod,4,1,outfile);
	      fwrite(periods,8,Nperiod,outfile);
	      for(j=0;j<Npeaks;j++)
		{
		  fwrite(&(aveper_whiten[j]),8,1,outfile);
		  fwrite(&(stdper_whiten[j]),8,1,outfile);
		  fwrite(periodogram_whiten[j],8,Nperiod,outfile);
		}
	    }
	  
	  fclose(outfile);
	}
      
      /* Set ERROR period peaks to 1. so that other routines that take the period from aov don't have trouble with the negative period. */
      for(k=0;k<Npeaks;k++)
	{
	  if(perpeaks[k] < ERROR_SCORE || isnan(perpeaks[k]))
	    perpeaks[k] = 1.;
	}
      
      free(t_cpy);
      free(mag_cpy);
      free(sig_cpy);
    }
#ifdef PARALLEL
  if(m_noave != NULL) { free(m_noave); m_noave = NULL;}
  if(t_nostart != NULL) { free(t_nostart); t_nostart = NULL;}
  if(weight != NULL) { free(weight); weight = NULL;}

  if(periods != NULL) free(periods);
  if(periodogram != NULL) free(periodogram);
  if(periodogram_whiten != NULL) {
    for(k=0; k <Npeaks+1; k++) {
      free(periodogram_whiten[k]);
    }
    free(periodogram_whiten);
  }
  if(aveper_whiten != NULL) free(aveper_whiten);
  if(stdper_whiten != NULL) free(stdper_whiten);
  if(harmA != NULL) free(harmA);
  if(harmB != NULL) free(harmB);
  if(t_mask != NULL) free(t_mask);
  if(mag_mask != NULL) free(mag_mask);
  if(sig_mask != NULL) free(sig_mask);

#endif
}

void RunAOVHarmCommand(ProgramData *p, Command *c, _AovHarm *AovHarm, int lcnum, int lc_name_num, int thisindex)
{
  int i1, i2;
  double d1;
  double *d1ptr, *d2ptr, *d3ptr;
  char outname[MAXLEN];

  if(AovHarm->operiodogram)
    {
      i1 = 0;
      i2 = 0;
      while(p->lcnames[lc_name_num][i1] != '\0')
	{
	  if(p->lcnames[lc_name_num][i1] == '/')
	    i2 = i1 + 1;
	  i1++;
	}
      sprintf(outname,"%s/%s%s",AovHarm->outdir,&p->lcnames[lc_name_num][i2],AovHarm->suffix);
    }

  if(AovHarm->Nharm_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    AovHarm->Nharm_vals[lcnum] = ceil(EvaluateExpression(lc_name_num, lcnum, 0, AovHarm->Nharm_expr));
  }
  else if(AovHarm->Nharm_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
    AovHarm->Nharm_vals[lcnum] = ceil(EvaluateVariable_Double(lc_name_num, lcnum, 0, AovHarm->Nharm_var));
  }
  else {
    AovHarm->Nharm_vals[lcnum] = AovHarm->Nharm;
  }

  if(AovHarm->minp_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    AovHarm->minp_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, AovHarm->minp_expr);
  }
  else if(AovHarm->minp_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
    AovHarm->minp_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, AovHarm->minp_var);
  }
  else {
    AovHarm->minp_vals[lcnum] = AovHarm->minp;
  }

  if(AovHarm->maxp_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    AovHarm->maxp_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, AovHarm->maxp_expr);
  }
  else if(AovHarm->maxp_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
    AovHarm->maxp_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, AovHarm->maxp_var);
  }
  else {
    AovHarm->maxp_vals[lcnum] = AovHarm->maxp;
  }

  if(AovHarm->subsample_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    AovHarm->subsample_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, AovHarm->subsample_expr);
  }
  else if(AovHarm->subsample_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
    AovHarm->subsample_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, AovHarm->subsample_var);
  }
  else {
    AovHarm->subsample_vals[lcnum] = AovHarm->subsample;
  }

  if(AovHarm->finetune_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    AovHarm->finetune_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, AovHarm->finetune_expr);
  }
  else if(AovHarm->finetune_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
    AovHarm->finetune_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, AovHarm->finetune_var);
  }
  else {
    AovHarm->finetune_vals[lcnum] = AovHarm->finetune;
  }


  if(AovHarm->fixperiodSNR)
    {
      if(AovHarm->fixperiodSNR_pertype == PERTYPE_AOV)
	{
	  i1=AovHarm->fixperiodSNR_lastaovindex;
	  if(c[i1-thisindex].cnum == CNUM_AOV)
	    AovHarm->fixperiodSNR_periods[lcnum][0] = c[i1-thisindex].Aov->peakperiods[lcnum][0];
	  else if(c[i1-thisindex].cnum == CNUM_HARMAOV)
	    AovHarm->fixperiodSNR_periods[lcnum][0] = c[i1-thisindex].AovHarm->peakperiods[lcnum][0];
	  
	}
      else if(AovHarm->fixperiodSNR_pertype == PERTYPE_LS)
	{
	  i1 = AovHarm->fixperiodSNR_lastaovindex;
	  AovHarm->fixperiodSNR_periods[lcnum][0] = c[i1-thisindex].AovHarm->peakperiods[lcnum][0];
	}
      else if(AovHarm->fixperiodSNR_pertype == PERTYPE_INJECTHARM)
	{
	  i1 = AovHarm->fixperiodSNR_lastaovindex;
	  AovHarm->fixperiodSNR_periods[lcnum][0] = c[i1-thisindex].Injectharm->periodinject[lcnum];
	}
      else if(AovHarm->fixperiodSNR_pertype == PERTYPE_FIX)
	{
	  AovHarm->fixperiodSNR_periods[lcnum][0] = AovHarm->fixperiodSNR_fixedperiod;
	}
      else if(AovHarm->fixperiodSNR_pertype == PERTYPE_FIXCOLUMN)
	{
	  getoutcolumnvalue(AovHarm->fixperiodSNR_linkedcolumn, lcnum, lc_name_num, VARTOOLS_TYPE_DOUBLE, &(AovHarm->fixperiodSNR_periods[lcnum][0]));
	}
      d1 = AovHarm->fixperiodSNR_periods[lcnum][0];
      d1ptr = &(AovHarm->fixperiodSNR_peakvalues[lcnum]);
      d2ptr = &(AovHarm->fixperiodSNR_peakSNR[lcnum]);
      d3ptr = &(AovHarm->fixperiodSNR_peakFAP[lcnum]);
    }
  else
    {
      d1 = 1.;
      d1ptr = NULL;
      d2ptr = NULL;
      d3ptr = NULL;
    }
  if(p->NJD[lcnum] > 1) {
    findPeaks_aovharm(p->t[lcnum], p->mag[lcnum], p->sig[lcnum], p->NJD[lcnum], AovHarm->peakperiods[lcnum], AovHarm->peakvalues[lcnum], AovHarm->peakSNR[lcnum], AovHarm->peakFAP[lcnum], AovHarm->peakNharm[lcnum], AovHarm->Npeaks, AovHarm->minp_vals[lcnum], AovHarm->maxp_vals[lcnum], AovHarm->subsample_vals[lcnum], AovHarm->finetune_vals[lcnum], AovHarm->operiodogram, outname, &AovHarm->aveaov[lcnum], &AovHarm->rmsaov[lcnum],AovHarm->aveaov_whiten[lcnum],AovHarm->rmsaov_whiten[lcnum],p->ascii, AovHarm->Nharm_vals[lcnum],AovHarm->whiten, AovHarm->clip, AovHarm->clipiter, AovHarm->fixperiodSNR, d1, d1ptr, d2ptr, d3ptr, lcnum, lc_name_num, AovHarm->usemask, AovHarm->maskvar);
  }
}
