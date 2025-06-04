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
  int sign, i, n;
  int N;
  double *resid, *err;
  double *resid_mask = NULL, *err_mask = NULL;

  if(N_in < 1)
    return 0.;

  if(!usemask) {
    N = N_in;
    resid = resid_in;
    err = err_in;
  } else {
    if((resid_mask = (double *) malloc(N_in * sizeof(double))) == NULL ||
       (err_mask = (double *) malloc(N_in * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
    N = 0;
    for(i=0; i < N_in; i++) {
      if(!isnan(resid_in[i]) && !isnan(err_in[i]) && err_in[i] > 0. &&
	 EvaluateVariable_Double(lclistnum, lcnum, i, maskvar) > VARTOOLS_MASK_TINY) {
	resid_mask[N] = resid_in[i];
	err_mask[N] = err_in[i];
	N++;
      }
    }
    resid = resid_mask;
    err = err_mask;
    if(N < 1) {
      free(resid_mask);
      free(err_mask);
      return(0.);
    }
  }
      

  sum1 = 0.;
  sum2 = 0.;
  for(i=0;i<N;i++)
    {
      if(!isnan(resid[i]) && !isnan(err[i]) && err[i] > 0.)
	{
	  sum1 += (long double) (resid[i] / (err[i]*err[i]));
	  sum2 += (long double) (1. / (err[i]*err[i]));
	}
    }

  a = 0.;
  if(N > 0)
    {
      ave = sum1 / sum2;
      resid[0] -= ave;
      sum1 = resid[0]*resid[0] / (err[0]*err[0]);
      sum2 = resid[0] / err[0];
      n = 1;
      sign = SIGN(resid[0]);
      for(i=1;i<N;i++)
	{
	  if(!isnan(resid[i]))
	    {
	      resid[i] -= ave;
	      if(SIGN(resid[i]) != sign)
		{
		  a += (sum2 * sum2);
		  sum2 = resid[i] / err[i];
		  n = 1;
		  sign = SIGN(resid[i]);
		}
	      else
		{
		  sum2 += resid[i] / err[i];
		  n++;
		}
	      sum1 += resid[i]*resid[i] / (err[i]*err[i]);
	    }
	}
      a += (sum2 * sum2);
      a = (a / sum1) - DEFAULT_ALARMVAL;
    }
  if(resid_mask != NULL) free(resid_mask);
  if(err_mask != NULL) free(err_mask);
  return(a);
}

