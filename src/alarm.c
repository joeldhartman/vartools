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
double doalarm(int N, double *resid, double *err)
{
  double a, ave;
  long double sum1, sum2;
  int sign, i, n;

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
  return(a);
}

