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

void medianfilter(int N, double *t, double *mag, double *sig, double timesize, int meanflag, int replace)
{
  /* N, t, mag and sig define the input light curve, the output light curve is stored in these variables as well.

We assume the input light curve is sorted in time.

meanflag - 0 means median, 1 means average, 2 means weighted average.

timesize - the median/average is calculated for all points within timesize of the observation

replace - 0 means the median is subtracted from each point.
          1 means each point is replaced with the running median.
  */

  int i, j, k, jstart, jstop;
  double *vals, *magout, medval, meanval, val1, val2;
  int Nvals;

  if(N <= 0)
    return;

  if((magout = (double *) malloc(N * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  if(!meanflag)
    {
      if((vals = (double *) malloc(N * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);

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
	  medval = median((jstop - jstart + 1), &mag[jstart]);
	  if(!replace)
	    magout[i] = mag[i] - medval;
	  else
	    magout[i] = medval;
	}
      memcpy(mag,magout,N*sizeof(double));
      free(vals);
    }
  else if(meanflag == 1)
    {
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
	  meanval = getmean((jstop - jstart + 1), &mag[jstart]);
	  if(!replace)
	    magout[i] = mag[i] - meanval;
	  else
	    magout[i] = meanval;
	}
      memcpy(mag,magout,N*sizeof(double));
    }
  else if(meanflag == 2)
    {
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
	  meanval = getweightedmean((jstop - jstart + 1), &mag[jstart], &sig[jstart]);
	  if(!replace)
	    magout[i] = mag[i] - meanval;
	  else
	    magout[i] = meanval;
	}
      memcpy(mag,magout,N*sizeof(double));
    }

  free(magout);
}
