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

double chisqtraptransit(double *a_, int ma, int N, double *t, double *mag, double *sig, void *userparams)
{
  double ph;
  double P;
  double q;
  double qingress;
  double in1_ph;
  double t0;
  double depth;
  double OOTmag;
  int i;
  double chisq, ph1, ph2;
  double mc, temp;

  P = a_[0]; q = a_[1]; qingress = a_[2]; in1_ph = a_[3];
  depth = a_[4]; OOTmag = a_[5];

  if(qingress > 0.5)
    qingress = 0.5;
  if(qingress < 0.)
    qingress = 0.;
  if(q > 1.)
    q = 1.;
  if(q < 0.)
    q = 0.;

  t0 = t[0] + in1_ph*P;

  chisq = 0.;
  ph1 = qingress*q;
  ph2 = q - qingress*q;
  for(i=0; i < N; i++) {
    if(!isnan(mag[i])) {
      ph = ((t[i] - t0)/P);
      ph = ph - floor(ph);
      if(ph > q) mc = OOTmag;
      else if(ph > ph1 && ph < ph2) mc = OOTmag + depth;
      else {
	if(ph < ph1) {
	  mc = OOTmag + depth*ph/ph1;
	}
	else {
	  mc = OOTmag + depth*(q - ph)/ph1;
	}
      }
      temp = (mag[i] - mc)/sig[i];
      chisq += temp*temp;
    }
  }
  return(chisq);
}

#define CONVERGENCELIMIT 0.0001
#define INITIALSTEP 0.1

void dofittrap_amoeba(int N, double *t, double *mag, double *sig, double P, double *q, double *qingress, double *in1_ph, double *in2_ph, double *depth, double *OOTmag)
{
  /* This function uses amoeba to fit a trapezoidal transit to a light curve.
     The input parameters are:
     N - number of data points
     t - times of observation
     mag - magnitude values
     sig - errors
     P - period of the transit
     q - duration / Period
     qingress - duration of ingress / duration of transit
     in1_ph - phase of transit start relative to the first time in the lc
     in2_ph - phase of transit stop relative to the first time in the lc
     depth - depth in magnitude of the transit
     OOTmag - out-of-transit magnitude level.

     The free parameters in the fit are q, qingress, in1_ph, depth and OOTmag.
     These parameters, together with in2_ph, are output.
  */

  int ma, i, j, nfunk;
  int ia[6], nvar;
  double y[6];
  double **p;

  double ftol;
  int amoeba_val;

  p = (double **) malloc(6*sizeof(double *));
  for(i=0; i < 6; i++) {
    p[i] = (double *) malloc(6 * sizeof(double));
  }
  ma = 6;
  ia[0] = 0; ia[1] = 1; ia[2] = 1; ia[3] = 1; ia[4] = 1; ia[5] = 1;
  nvar = 5;

  p[0][0] = P; p[0][1] = *q; p[0][2] = *qingress; p[0][3] = *in1_ph;
  p[0][4] = *depth; p[0][5] = *OOTmag;

  for(i=1; i <= 5; i++) {
    p[i][0] = P;
    for(j=1; j <= 5; j++) {
      if(i == j) {
	p[i][j] = p[0][j]*(1. + INITIALSTEP);
      } else
	p[i][j] = p[0][j];
    }
  }

  for(i=0; i < 6; i++)
    y[i] = chisqtraptransit(p[i], ma, N, t, mag, sig, NULL);

  ftol = CONVERGENCELIMIT;

  amoeba_val = amoeba(p, y, ia, ma, ftol, &chisqtraptransit, &nfunk, 0, N, t, mag, sig, NULL);

  if(!amoeba_val && nvar > 0)
    {
      j = 0;
      for(i=0; i <= nvar; i++){
	if(y[i] < y[j])
	  j=i;
      }
      *q = p[j][1]; *qingress = p[j][2]; *in1_ph = p[j][3];
      *depth = p[j][4]; *OOTmag = p[j][5];
      *in2_ph = *in1_ph + *q;
    }
  else
    {
      *qingress = 0.;
    }

  for(i=0; i < 6; i++)
    free(p[i]);
  free(p);

}

void dofittrap_amoeba_fixdur(int N, double *t, double *mag, double *sig, double P, double q, double *qingress, double in1_ph, double in2_ph, double *depth, double *OOTmag)
{
  /* This function uses amoeba to fit a trapezoidal transit to a light curve.
     The input parameters are:
     N - number of data points
     t - times of observation
     mag - magnitude values
     sig - errors
     P - period of the transit
     q - duration / Period
     qingress - duration of ingress / duration of transit
     in1_ph - phase of transit start relative to the first time in the lc
     in2_ph - phase of transit stop relative to the first time in the lc
     depth - depth in magnitude of the transit
     OOTmag - out-of-transit magnitude level.

     The free parameters in the fit are q, qingress, in1_ph, depth and OOTmag.
     These parameters, together with in2_ph, are output.
  */

  int ma, i, j, nfunk;
  int ia[6], nvar;
  double y[6];
  double **p;

  double ftol;
  int amoeba_val;

  p = (double **) malloc(4*sizeof(double *));
  for(i=0; i < 4; i++) {
    p[i] = (double *) malloc(6 * sizeof(double));
  }
  ma = 6;
  ia[0] = 0; ia[1] = 0; ia[2] = 1; ia[3] = 0; ia[4] = 1; ia[5] = 1;
  nvar = 3;

  for(i=0; i <= 3; i++) {
    p[i][0] = P;
    p[i][1] = q;
    p[i][2] = *qingress;
    p[i][3] = in1_ph;
    p[i][4] = *depth;
    p[i][5] = *OOTmag;
  }
  p[1][2] = p[1][2]*(1. + INITIALSTEP);
  p[2][4] = p[2][4]*(1. + INITIALSTEP);
  p[3][5] = p[3][5]*(1. + INITIALSTEP);
  
  for(i=0; i < 4; i++)
    y[i] = chisqtraptransit(p[i], ma, N, t, mag, sig, NULL);

  ftol = CONVERGENCELIMIT;

  amoeba_val = amoeba(p, y, ia, ma, ftol, &chisqtraptransit, &nfunk, 0, N, t, mag, sig, NULL);

  if(!amoeba_val && nvar > 0)
    {
      j = 0;
      for(i=0; i <= nvar; i++){
	if(y[i] < y[j])
	  j=i;
      }
      *qingress = p[j][2];
      *depth = p[j][4]; *OOTmag = p[j][5];
    }
  else
    {
      *qingress = 0.;
    }
  
  for(i=0; i < 4; i++)
    free(p[i]);
  free(p);

}

