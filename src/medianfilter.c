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

/* This file contains the function to median filter a light curve */

/* Public/ABI-stable entry point (7 arguments, exported to user libraries as
   VARTOOLS_medianfilter).  Masking support is provided by domedianfilter;
   keep this signature unchanged so the userlib ABI is not broken. */
void medianfilter(int N, double *t, double *mag, double *sig, double timesize, int meanflag, int replace)
{
  domedianfilter(N, t, mag, sig, timesize, meanflag, replace, 0, NULL, 0, 0);
}

void domedianfilter(int N, double *t, double *mag, double *sig, double timesize, int meanflag, int replace, int usemask, _Variable *maskvar, int lcindex, int threadindex)
{
  /* N, t, mag and sig define the input light curve, the output light curve is stored in these variables as well.

We assume the input light curve is sorted in time.

meanflag - 0 means median, 1 means average, 2 means weighted average.

timesize - the median/average is calculated for all points within timesize of the observation

replace - 0 means the median is subtracted from each point.
          1 means each point is replaced with the running median.

usemask - if non-zero, only points j with maskvar > VARTOOLS_MASK_TINY
          contribute to the window statistic.  The filter value is still
          applied to every point; a point whose window contains no unmasked
          points is left unchanged.
  */

  int i, j, jstart, jstop;
  double *vals, *wvals = NULL, *magout, filtval;
  int Nvals;

  if(N <= 0)
    return;

  if((magout = (double *) malloc(N * sizeof(double))) == NULL ||
     (vals = (double *) malloc(N * sizeof(double))) == NULL)
    vt_error(ERR_MEMALLOC);
  if(meanflag == 2) {
    if((wvals = (double *) malloc(N * sizeof(double))) == NULL)
      vt_error(ERR_MEMALLOC);
  }

  jstart = 0;
  jstop = 0;
  for(i=0;i<N;i++)
    {
      while(t[i] - t[jstart] > timesize)
	jstart++;
      if(jstop < i)
	jstop = i;
      while(jstop < N ? (t[jstop] - t[i] < timesize) : 0)
	jstop++;
      jstop--;

      /* Gather the points in the window [jstart, jstop], keeping only the
	 unmasked ones when a mask is in use.  With no mask this collects the
	 full contiguous window, so the result is identical to the original. */
      Nvals = 0;
      for(j = jstart; j <= jstop; j++)
	{
	  if(usemask && !(EvaluateVariable_Double(lcindex, threadindex, j, maskvar) > VARTOOLS_MASK_TINY))
	    continue;
	  vals[Nvals] = mag[j];
	  if(meanflag == 2)
	    wvals[Nvals] = sig[j];
	  Nvals++;
	}

      if(Nvals <= 0)
	{
	  /* No unmasked points in the window; leave this point unchanged. */
	  magout[i] = mag[i];
	  continue;
	}

      if(!meanflag)
	filtval = median(Nvals, vals);
      else if(meanflag == 1)
	filtval = getmean(Nvals, vals);
      else
	filtval = getweightedmean(Nvals, vals, wvals);

      if(!replace)
	magout[i] = mag[i] - filtval;
      else
	magout[i] = filtval;
    }
  memcpy(mag,magout,N*sizeof(double));

  free(magout);
  free(vals);
  if(wvals != NULL)
    free(wvals);
}
