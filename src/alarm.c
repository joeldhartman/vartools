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

/* This file contains functions for computing the Alarm variability statistic for the program vartools by J. Hartman */

#ifdef SIGN
#undef SIGN
#endif

#define SIGN(A) ((A) >= 0 ? (1) : (-1))
#define ROUND(A) ((A) - floor((double) A) > 0.5 ? ((A) - floor((double) A) + 1) : ((A) - floor((double) A)));
#define DEFAULT_ALARMVAL 2.2732395
double doalarm(int N_in, double *resid_in, double *err_in, int lcnum, int lclistnum, int usemask, _Variable *maskvar)
{
  double a, ave;
  long double sum1, sum2;
  int sign, i;
  int N;
  double *resid_mask = NULL, *err_mask = NULL;

  if(N_in < 1)
    return 0.;

  if((resid_mask = (double *) malloc(N_in * sizeof(double))) == NULL ||
       (err_mask = (double *) malloc(N_in * sizeof(double))) == NULL)
      vt_error(ERR_MEMALLOC);

  /* Always create a filtered copy, regardless a mask is used. */
  N = 0;
  for(i=0; i < N_in; i++) {
    if(usemask && !(EvaluateVariable_Double(lclistnum, lcnum, i, maskvar) > VARTOOLS_MASK_TINY))
      continue;

    if(!isfinite(resid_in[i]) || !isfinite(err_in[i]) || err_in[i] <= 0.0)
      continue;

    resid_mask[N] = resid_in[i];
    err_mask[N] = err_in[i];
    N++;
  }

  if(N < 2) {
    free(resid_mask);
    free(err_mask);
    return(0.);
  }

  /* Calculate the inverse-variance-weighted mean. */
  sum1 = 0.;
  sum2 = 0.;
  for(i=0;i<N;i++)
    {
	  sum1 += (long double) (resid_mask[i] / (err_mask[i]*err_mask[i]));
	  sum2 += (long double) (1. / (err_mask[i]*err_mask[i]));
    }

  if(!(sum2 > 0.0L) || !isfinite(sum1) || !isfinite(sum2)) {
    free(resid_mask);
    free(err_mask);
    return(0.);
  }

  ave = (double) (sum1 / sum2);

  if(!isfinite(ave)) {
    free(resid_mask);
    free(err_mask);
    return 0.0;
  }

  a = 0.;
  resid_mask[0] -= ave;
  sum1 = resid_mask[0]*resid_mask[0] / (err_mask[0]*err_mask[0]);
  sum2 = resid_mask[0] / err_mask[0];
  sign = SIGN(resid_mask[0]);
  for(i=1;i<N;i++)
  {
	resid_mask[i] -= ave;
    if(SIGN(resid_mask[i]) != sign) {
	  a += (double) (sum2 * sum2);
	  sum2 = resid_mask[i] / err_mask[i];
	  sign = SIGN(resid_mask[i]);
      }
    else {
	  sum2 += resid_mask[i] / err_mask[i];
      }
	sum1 += resid_mask[i]*resid_mask[i] / (err_mask[i]*err_mask[i]);
  }
  a += (double) (sum2 * sum2);
  a = (double) (a / sum1) - DEFAULT_ALARMVAL;

  free(resid_mask);
  free(err_mask);

  return(a);
}

