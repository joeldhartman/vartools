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

/* This file contains the function to bin a light curve */

#define GROWBINVECTOR(indx,typ,vec) do { \
    if(!Nptr[(indx)]) { \
        if((ptrindx[(indx)] = (int *) malloc(sizeof(int))) == NULL || \
           ((vec) = malloc(sizeof(typ))) == NULL) \
        error(ERR_MEMALLOC); \
    } else { \
        if((ptrindx[(indx)] = (int *) realloc(ptrindx[(indx)],(Nptr[(indx)]+1)*sizeof(int))) == NULL || \
           ((vec) = realloc(vec, (Nptr[(indx)]+1)*sizeof(typ))) == NULL) \
        error(ERR_MEMALLOC); \
    } \
    ptrindx[(indx)][Nptr[(indx)]] = k; \
    Nptr[(indx)]++; \
  } while(0)

#define GROWPTRVECTOR(indx) do { \
    if(!Nptr[(indx)]) { \
      if((ptrindx[(indx)] = (int *) malloc(sizeof(int))) == NULL)	\
        error(ERR_MEMALLOC); \
    } else { \
      if((ptrindx[(indx)] = (int *) realloc(ptrindx[(indx)],(Nptr[(indx)]+1)*sizeof(int))) == NULL) \
        error(ERR_MEMALLOC);						\
    } \
    ptrindx[(indx)][Nptr[(indx)]] = k; \
    Nptr[(indx)]++; \
  } while(0)


void binlc(ProgramData *p, int lcnum, int medflag, int binsize_Nbins_flag, double binsize, int Nbins, int firstbinflag, double firstbin, int tflag)
{
  /* N, t, mag and sig define the input light curve, the output binned light curve is stored in these variables as well.

We assume the input light curve is sorted in time.

medflag - 0 means average bin, 1 means median bin, 2 means weighted average bin

binsize_Nbins_flag - 0 means that the binsize (in units of time) is specified, 1 means the number of bins is specified.

firstbinflag - 0 means that the first bin starts at t[0], 1 means the first bin starts at t[0] - firstbin/binsize

tflag - 0 means the output time for each bin is the time at the center of the bin, 1 means to take the average of the times of points that fall within the bin, 2 means to take the median.
  */

  int *N;
  double *t, *mag, *sig;

  int i, j, k, n, u, otherdata = 0;
#ifdef PARALLEL
  double **bin_mag = NULL, **bin_sig = NULL, **bin_time = NULL;
  int size1 = 0, size2 = 0, *nbin = NULL;
  int *startindx = NULL;
#else
  static double **bin_mag, **bin_sig, **bin_time;
  static int *startindx;
  static int size1 = 0, size2 = 0, *nbin;
#endif
  double ***bin_dblptr = NULL, ****bin_dbl2ptr = NULL;
  short ***bin_shortptr = NULL, ****bin_short2ptr = NULL;
  int ***bin_intptr = NULL, ****bin_int2ptr = NULL;
  char ***bin_charptr = NULL, ****bin_char2ptr = NULL;
  char ****bin_stringptr = NULL, *****bin_string2ptr = NULL;
  float ***bin_floatptr = NULL, ****bin_float2ptr = NULL;
  long ***bin_longptr = NULL, ****bin_long2ptr = NULL;

  double t0, t1, var1, var2;
  int binnum, Nc;   

  double ***dblptr;
  double ****dbl2ptr;
  short ***shortptr;
  short ****short2ptr;
  int ***intptr;
  int ****int2ptr;
  char ***charptr;
  char ****char2ptr;
  char ****stringptr;
  char *****string2ptr;
  float ***floatptr;
  float ****float2ptr;
  long ***longptr;
  long ****long2ptr;

  int indx;

  int Nptr[14];
  int* ptrindx[14];
  int ptr_datasize1[14];
  int ptr_datasize2[14];
  int ptr_datasize3[14];
  
  _DataFromLightCurve *d;

  

  for(i=0; i < 14; i++) Nptr[i] = 0;

  N = &p->NJD[lcnum];
  t = p->t[lcnum];
  mag = p->mag[lcnum];
  sig = p->sig[lcnum];

  n = *N;
  
  if(n <= 0)
    return;

  /* First determine the number of bins that we need to allocate memory for */
  t0 = t[0];
  t1 = t[n - 1];  

  if(t1 < t0)
    error(ERR_UNSORTEDLIGHTCURVE);

  if(binsize_Nbins_flag)
    {
      binsize = (t1 - t0)/(double) Nbins;
    }
  else
    {
      Nbins = ceil((t1 - t0)/binsize);
    }
  
  /* adjust the initial time if needed */
  if(firstbinflag)
    {
      t0 = t0 - firstbin/binsize;
    }

  /* Only run if Nbins > 0 */
  if(Nbins <= 0) {
    *N = 0;
    return;
  }

  /* Figure out what light-curve vectors exist in addition to those
     given above, and prepare storage for sorting them into bins */
  for(k=0; k < p->NDataFromLightCurve; k++) {
    d = &(p->DataFromLightCurve[k]);
    Nc = d->Ncolumns;
    /* skip vectors corresponding to t, mag or sig */
    if(Nc == 0 && d->datatype == VARTOOLS_TYPE_DOUBLE) {
      if(*((double ***) d->dataptr) == p->t ||
	 *((double ***) d->dataptr) == p->mag ||
	 *((double ***) d->dataptr) == p->sig)
	continue;
    }
    otherdata = 1;
    switch(d->datatype)
      {
      case VARTOOLS_TYPE_DOUBLE:
	if(Nc == 0)
	  GROWBINVECTOR(0,double **,bin_dblptr);
	else 
	  GROWBINVECTOR(1, double ***, bin_dbl2ptr);
	break;
      case VARTOOLS_TYPE_CONVERTJD:
	if(Nc == 0)
	  GROWBINVECTOR(0,double **,bin_dblptr);
	else
	  GROWBINVECTOR(1,double ***,bin_dbl2ptr);
	break;
      case VARTOOLS_TYPE_SHORT:
	if(Nc == 0) {
	  GROWBINVECTOR(2,short **,bin_shortptr);
	} else {
	  GROWBINVECTOR(3,short ***,bin_short2ptr);
	}
	break;
      case VARTOOLS_TYPE_INT:
	if(Nc == 0) {
	  GROWBINVECTOR(4,int **,bin_intptr);
	} else {
	  GROWBINVECTOR(5,int ***,bin_int2ptr);
	}
	break;
      case VARTOOLS_TYPE_CHAR:
	
	if(Nc == 0) {
	  GROWPTRVECTOR(6);
	} else {
	  GROWPTRVECTOR(7);
	}
	break;
      case VARTOOLS_TYPE_STRING:
	if(Nc == 0) {
	  GROWPTRVECTOR(8);
	} else {
	  GROWPTRVECTOR(9);
	}
	break;
      case VARTOOLS_TYPE_FLOAT:
	if(Nc == 0) {
	  GROWBINVECTOR(10,float **,bin_floatptr);
	} else {
	  GROWBINVECTOR(11,float ***,bin_float2ptr);
	}
	break;
      case VARTOOLS_TYPE_LONG:
	if(Nc == 0) {
	  GROWBINVECTOR(12,long **,bin_longptr);
	} else {
	  GROWBINVECTOR(13,long ***,bin_long2ptr);
	}
	break;
      default:
	error(ERR_BADTYPE);
	break;
      }
  }
    

  /* Only run this process if there are more than 0 bins to make */
  if(Nbins > 0)
    {
      
      if(Nbins > size1)
	{
	  if(size1 == 0)
	    {
	      size1 = Nbins;
	      size2 = n;
	      if((bin_mag = (double **) malloc(Nbins * sizeof(double *))) == NULL ||
		 (bin_sig = (double **) malloc(Nbins * sizeof(double *))) == NULL ||
		 (bin_time = (double **) malloc(Nbins * sizeof(double *))) == NULL ||
		 (nbin = (int *) malloc(Nbins * sizeof(int))) == NULL ||
		 (startindx = (int *) malloc(Nbins * sizeof(int))) == NULL)
		error(ERR_MEMALLOC);
	      for(i=0;i<Nbins;i++)
		{
		  if((bin_mag[i] = (double *) malloc(n * sizeof(double))) == NULL ||
		     (bin_sig[i] = (double *) malloc(n * sizeof(double))) == NULL ||
		     (bin_time[i] = (double *) malloc(n * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	    }
	  else
	    {
	      if((bin_mag = (double **) realloc(bin_mag, Nbins * sizeof(double *))) == NULL ||
		 (bin_sig = (double **) realloc(bin_sig, Nbins * sizeof(double *))) == NULL ||
		 (bin_time = (double **) realloc(bin_time, Nbins * sizeof(double *))) == NULL ||
		 (nbin = (int *) realloc(nbin, Nbins * sizeof(int))) == NULL ||
		 (startindx = (int *) realloc(startindx, Nbins * sizeof(int))) == NULL)
		error(ERR_MEMALLOC);
	      if(n > size2)
		{
		  size2 = n;
		  for(i=0;i<size1;i++)
		    {
		      if((bin_mag[i] = (double *) realloc(bin_mag[i], n * sizeof(double))) == NULL ||
			 (bin_sig[i] = (double *) realloc(bin_sig[i], n * sizeof(double))) == NULL ||
			 (bin_time[i] = (double *) realloc(bin_time[i], n * sizeof(double))) == NULL)
			error(ERR_MEMALLOC);
		    }
		  for(;i<Nbins;i++)
		    {
		      if((bin_mag[i] = (double *) malloc(n * sizeof(double))) == NULL ||
			 (bin_sig[i] = (double *) malloc(n * sizeof(double))) == NULL ||
			 (bin_time[i] = (double *) malloc(n * sizeof(double))) == NULL)
			error(ERR_MEMALLOC);
		    }
		}
	      else
		{
		  for(i=size1;i<Nbins;i++)
		    {
		      if((bin_mag[i] = (double *) malloc(size2 * sizeof(double))) == NULL ||
			 (bin_sig[i] = (double *) malloc(size2 * sizeof(double))) == NULL ||
			 (bin_time[i] = (double *) malloc(size2 * sizeof(double))) == NULL)
			error(ERR_MEMALLOC);
		    }
		}
	      size1 = Nbins;
	    }
	}
      if(n > size2)
	{
	  size2 = n;
	  for(i=0;i<size1;i++)
	    {
	      if((bin_mag[i] = (double *) realloc(bin_mag[i], n * sizeof(double))) == NULL ||
		 (bin_sig[i] = (double *) realloc(bin_sig[i], n * sizeof(double))) == NULL ||
		 (bin_time[i] = (double *) realloc(bin_time[i], n * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	}

      /* Allocate memory for the other data vectors */
      if(otherdata) {
	for(k=0;k<Nptr[0];k++) {
	  if((bin_dblptr[k] = (double **) malloc(Nbins * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(i=0; i < Nbins; i++) {
	    if((bin_dblptr[k][i] = (double *) malloc(n * sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	}
	for(k=0;k<Nptr[1];k++) {
	  d = &(p->DataFromLightCurve[ptrindx[1][k]]);
	  if((bin_dbl2ptr[k] = (double ***) malloc(d->Ncolumns*sizeof(double **))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < d->Ncolumns; j++) {
	    if((bin_dbl2ptr[k][j] = (double **) malloc(Nbins * sizeof(double *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(i=0; i < Nbins; i++) {
	      if((bin_dbl2ptr[k][j][i] = (double *) malloc(n * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  }
	}
	for(k=0;k<Nptr[2];k++) {
	  if((bin_shortptr[k] = (short **) malloc(Nbins * sizeof(short *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(i=0; i < Nbins; i++) {
	    if((bin_shortptr[k][i] = (short *) malloc(n * sizeof(short))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	}
	for(k=0;k<Nptr[3];k++) {
	  d = &(p->DataFromLightCurve[ptrindx[3][k]]);
	  if((bin_short2ptr[k] = (short ***) malloc(d->Ncolumns*sizeof(short **))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < d->Ncolumns; j++) {
	    if((bin_short2ptr[k][j] = (short **) malloc(Nbins * sizeof(short *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(i=0; i < Nbins; i++) {
	      if((bin_short2ptr[k][j][i] = (short *) malloc(n * sizeof(short))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  }
	}
	for(k=0;k<Nptr[4];k++) {
	  if((bin_intptr[k] = (int **) malloc(Nbins * sizeof(int *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(i=0; i < Nbins; i++) {
	    if((bin_intptr[k][i] = (int *) malloc(n * sizeof(int))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	}
	for(k=0;k<Nptr[5];k++) {
	  d = &(p->DataFromLightCurve[ptrindx[5][k]]);
	  if((bin_int2ptr[k] = (int ***) malloc(d->Ncolumns*sizeof(int **))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < d->Ncolumns; j++) {
	    if((bin_int2ptr[k][j] = (int **) malloc(Nbins * sizeof(int *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(i=0; i < Nbins; i++) {
	      if((bin_int2ptr[k][j][i] = (int *) malloc(n * sizeof(int))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  }
	}
	for(k=0;k<Nptr[10];k++) {
	  if((bin_floatptr[k] = (float **) malloc(Nbins * sizeof(float *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(i=0; i < Nbins; i++) {
	    if((bin_floatptr[k][i] = (float *) malloc(n * sizeof(float))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	}
	for(k=0;k<Nptr[11];k++) {
	  d = &(p->DataFromLightCurve[ptrindx[11][k]]);
	  if((bin_float2ptr[k] = (float ***) malloc(d->Ncolumns*sizeof(float **))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < d->Ncolumns; j++) {
	    if((bin_float2ptr[k][j] = (float **) malloc(Nbins * sizeof(float *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(i=0; i < Nbins; i++) {
	      if((bin_float2ptr[k][j][i] = (float *) malloc(n * sizeof(float))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  }
	}
	for(k=0;k<Nptr[12];k++) {
	  if((bin_longptr[k] = (long **) malloc(Nbins * sizeof(long *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(i=0; i < Nbins; i++) {
	    if((bin_longptr[k][i] = (long *) malloc(n * sizeof(long))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	}
	for(k=0;k<Nptr[13];k++) {
	  d = &(p->DataFromLightCurve[ptrindx[13][k]]);
	  if((bin_long2ptr[k] = (long ***) malloc(d->Ncolumns*sizeof(long **))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < d->Ncolumns; j++) {
	    if((bin_long2ptr[k][j] = (long **) malloc(Nbins * sizeof(long *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(i=0; i < Nbins; i++) {
	      if((bin_long2ptr[k][j][i] = (long *) malloc(n * sizeof(long))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  }
	}
      }

      for(i=0;i<Nbins;i++)
	nbin[i] = 0;
      
      /* Sort the data into the bin vectors */
      for(i=0; i<n; i++)
	{
	  binnum = floor((t[i] - t0)/binsize);
	  if(binnum >= Nbins)
	    binnum = Nbins - 1;
	  if(binnum < 0)
	    binnum = 0;
	  
	  j = nbin[binnum];
	  bin_mag[binnum][j] = mag[i];
	  bin_time[binnum][j] = t[i];
	  bin_sig[binnum][j] = sig[i];
	  if(!j) {
	    startindx[binnum] = i;
	  }
	  if(otherdata) {
	    for(k=0;k<Nptr[0];k++) {
	      d = &(p->DataFromLightCurve[ptrindx[0][k]]);
	      dblptr = (double ***) d->dataptr;
	      bin_dblptr[k][binnum][j] = (*dblptr)[lcnum][i];
	    }
	    for(k=0;k<Nptr[1];k++) {
	      d = &(p->DataFromLightCurve[ptrindx[1][k]]);
	      dbl2ptr = (double ****) d->dataptr;
	      for(u=0; u < d->Ncolumns; u++) {
		bin_dbl2ptr[k][u][binnum][j] = (*dbl2ptr)[lcnum][u][i];
	      }
	    }
	    for(k=0;k<Nptr[2];k++) {
	      d = &(p->DataFromLightCurve[ptrindx[2][k]]);
	      shortptr = (short ***) d->dataptr;
	      bin_shortptr[k][binnum][j] = (*shortptr)[lcnum][i];
	    }
	    for(k=0;k<Nptr[3];k++) {
	      d = &(p->DataFromLightCurve[ptrindx[3][k]]);
	      short2ptr = (short ****) d->dataptr;
	      for(u=0; u < d->Ncolumns; u++) {
		bin_short2ptr[k][u][binnum][j] = (*short2ptr)[lcnum][u][i];
	      }
	    }
	    for(k=0;k<Nptr[4];k++) {
	      d = &(p->DataFromLightCurve[ptrindx[4][k]]);
	      intptr = (int ***) d->dataptr;
	      bin_intptr[k][binnum][j] = (*intptr)[lcnum][i];
	    }
	    for(k=0;k<Nptr[5];k++) {
	      d = &(p->DataFromLightCurve[ptrindx[5][k]]);
	      int2ptr = (int ****) d->dataptr;
	      for(u=0; u < d->Ncolumns; u++) {
		bin_int2ptr[k][u][binnum][j] = (*int2ptr)[lcnum][u][i];
	      }
	    }
	    for(k=0;k<Nptr[10];k++) {
	      d = &(p->DataFromLightCurve[ptrindx[10][k]]);
	      floatptr = (float ***) d->dataptr;
	      bin_floatptr[k][binnum][j] = (*floatptr)[lcnum][i];
	    }
	    for(k=0;k<Nptr[11];k++) {
	      d = &(p->DataFromLightCurve[ptrindx[11][k]]);
	      float2ptr = (float ****) d->dataptr;
	      for(u=0; u < d->Ncolumns; u++) {
		bin_float2ptr[k][u][binnum][j] = (*float2ptr)[lcnum][u][i];
	      }
	    }
	    for(k=0;k<Nptr[12];k++) {
	      d = &(p->DataFromLightCurve[ptrindx[12][k]]);
	      longptr = (long ***) d->dataptr;
	      bin_longptr[k][binnum][j] = (*longptr)[lcnum][i];
	    }
	    for(k=0;k<Nptr[13];k++) {
	      d = &(p->DataFromLightCurve[ptrindx[13][k]]);
	      long2ptr = (long ****) d->dataptr;
	      for(u=0; u < d->Ncolumns; u++) {
		bin_long2ptr[k][u][binnum][j] = (*long2ptr)[lcnum][u][i];
	      }
	    }
	  }
	  
	      
	  nbin[binnum]++;
	}
      
      /* Now fill out the new light curve */
      for(i=0, j=0; i < Nbins; i++)
	{
	  if(nbin[i] > 0)
	    {
	      if(!medflag)
		{
		  mag[j] = getmean(nbin[i],bin_mag[i]);
		  var1 = 0.;
		  for(k=0;k<nbin[i];k++)
		    {
		      var1 += bin_sig[i][k]*bin_sig[i][k];
		    }
		  sig[j] = sqrt(var1 / (double) (nbin[i]*nbin[i]));
		  if(otherdata) {
		    for(k=0;k<Nptr[0];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[0][k]]);
		      dblptr = (double ***) d->dataptr;
		      (*dblptr)[lcnum][j] = getmean(nbin[i],bin_dblptr[k][i]);
		    }
		    for(k=0;k<Nptr[1];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[1][k]]);
		      dbl2ptr = (double ****) d->dataptr;
		      for(u=0; u < d->Ncolumns; u++) {
			(*dbl2ptr)[lcnum][u][j] = getmean(nbin[i],bin_dbl2ptr[k][u][i]);
		      }
		    }
		    for(k=0;k<Nptr[2];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[2][k]]);
		      shortptr = (short ***) d->dataptr;
		      (*shortptr)[lcnum][j] = getmean_short(nbin[i],bin_shortptr[k][i]);
		    }
		    for(k=0;k<Nptr[3];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[3][k]]);
		      short2ptr = (short ****) d->dataptr;
		      for(u=0; u < d->Ncolumns; u++) {
			(*short2ptr)[lcnum][u][j] = getmean_short(nbin[i],bin_short2ptr[k][u][i]);
		      }
		    }
		    for(k=0;k<Nptr[4];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[4][k]]);
		      intptr = (int ***) d->dataptr;
		      (*intptr)[lcnum][j] = getmean_int(nbin[i],bin_intptr[k][i]);
		    }
		    for(k=0;k<Nptr[5];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[5][k]]);
		      int2ptr = (int ****) d->dataptr;
		      for(u=0; u < d->Ncolumns; u++) {
			(*int2ptr)[lcnum][u][j] = getmean_int(nbin[i],bin_int2ptr[k][u][i]);
		      }
		    }
		    for(k=0;k<Nptr[10];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[10][k]]);
		      floatptr = (float ***) d->dataptr;
		      (*floatptr)[lcnum][j] = getmean_float(nbin[i],bin_floatptr[k][i]);
		    }
		    for(k=0;k<Nptr[11];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[11][k]]);
		      float2ptr = (float ****) d->dataptr;
		      for(u=0; u < d->Ncolumns; u++) {
			(*float2ptr)[lcnum][u][j] = getmean_float(nbin[i],bin_float2ptr[k][u][i]);
		      }
		    }
		    for(k=0;k<Nptr[12];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[12][k]]);
		      longptr = (long ***) d->dataptr;
		      (*longptr)[lcnum][j] = getmean_long(nbin[i],bin_longptr[k][i]);
		    }
		    for(k=0;k<Nptr[13];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[13][k]]);
		      long2ptr = (long ****) d->dataptr;
		      for(u=0; u < d->Ncolumns; u++) {
			(*long2ptr)[lcnum][u][j] = getmean_long(nbin[i],bin_long2ptr[k][u][i]);
		      }
		    }
		  }
		    
		}
	      else if(medflag == 1)
		{
		  mag[j] = median(nbin[i],bin_mag[i]);
		  var1 = 0.;
		  for(k=0;k<nbin[i];k++)
		    {
		      var1 += bin_sig[i][k]*bin_sig[i][k];
		    }
		  sig[j] = 1.253 * sqrt(var1 / (double) (nbin[i]*nbin[i]));
		  if(otherdata) {
		    for(k=0;k<Nptr[0];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[0][k]]);
		      dblptr = (double ***) d->dataptr;
		      (*dblptr)[lcnum][j] = median(nbin[i],bin_dblptr[k][i]);
		    }
		    for(k=0;k<Nptr[1];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[1][k]]);
		      dbl2ptr = (double ****) d->dataptr;
		      for(u=0; u < d->Ncolumns; u++) {
			(*dbl2ptr)[lcnum][u][j] = median(nbin[i],bin_dbl2ptr[k][u][i]);
		      }
		    }
		    for(k=0;k<Nptr[2];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[2][k]]);
		      shortptr = (short ***) d->dataptr;
		      (*shortptr)[lcnum][j] = median_short(nbin[i],bin_shortptr[k][i]);
		    }
		    for(k=0;k<Nptr[3];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[3][k]]);
		      short2ptr = (short ****) d->dataptr;
		      for(u=0; u < d->Ncolumns; u++) {
			(*short2ptr)[lcnum][u][j] = median_short(nbin[i],bin_short2ptr[k][u][i]);
		      }
		    }
		    for(k=0;k<Nptr[4];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[4][k]]);
		      intptr = (int ***) d->dataptr;
		      (*intptr)[lcnum][j] = median_int(nbin[i],bin_intptr[k][i]);
		    }
		    for(k=0;k<Nptr[5];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[5][k]]);
		      int2ptr = (int ****) d->dataptr;
		      for(u=0; u < d->Ncolumns; u++) {
			(*int2ptr)[lcnum][u][j] = median_int(nbin[i],bin_int2ptr[k][u][i]);
		      }
		    }
		    for(k=0;k<Nptr[10];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[10][k]]);
		      floatptr = (float ***) d->dataptr;
		      (*floatptr)[lcnum][j] = median_float(nbin[i],bin_floatptr[k][i]);
		    }
		    for(k=0;k<Nptr[11];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[11][k]]);
		      float2ptr = (float ****) d->dataptr;
		      for(u=0; u < d->Ncolumns; u++) {
			(*float2ptr)[lcnum][u][j] = median_float(nbin[i],bin_float2ptr[k][u][i]);
		      }
		    }
		    for(k=0;k<Nptr[12];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[12][k]]);
		      longptr = (long ***) d->dataptr;
		      (*longptr)[lcnum][j] = median_long(nbin[i],bin_longptr[k][i]);
		    }
		    for(k=0;k<Nptr[13];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[13][k]]);
		      long2ptr = (long ****) d->dataptr;
		      for(u=0; u < d->Ncolumns; u++) {
			(*long2ptr)[lcnum][u][j] = median_long(nbin[i],bin_long2ptr[k][u][i]);
		      }
		    }
		  }
		}
	      else if(medflag == 2)
		{
		  if(otherdata) {
		    for(k=0;k<Nptr[0];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[0][k]]);
		      dblptr = (double ***) d->dataptr;
		      (*dblptr)[lcnum][j] = getweightedmean(nbin[i],bin_dblptr[k][i],bin_sig[i]);
		    }
		    for(k=0;k<Nptr[1];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[1][k]]);
		      dbl2ptr = (double ****) d->dataptr;
		      for(u=0; u < d->Ncolumns; u++) {
			(*dbl2ptr)[lcnum][u][j] = getweightedmean(nbin[i],bin_dbl2ptr[k][u][i],bin_sig[i]);
		      }
		    }
		    for(k=0;k<Nptr[2];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[2][k]]);
		      shortptr = (short ***) d->dataptr;
		      (*shortptr)[lcnum][j] = getweightedmean_short(nbin[i],bin_shortptr[k][i],bin_sig[i]);
		    }
		    for(k=0;k<Nptr[3];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[3][k]]);
		      short2ptr = (short ****) d->dataptr;
		      for(u=0; u < d->Ncolumns; u++) {
			(*short2ptr)[lcnum][u][j] = getweightedmean_short(nbin[i],bin_short2ptr[k][u][i],bin_sig[i]);
		      }
		    }
		    for(k=0;k<Nptr[4];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[4][k]]);
		      intptr = (int ***) d->dataptr;
		      (*intptr)[lcnum][j] = getweightedmean_int(nbin[i],bin_intptr[k][i],bin_sig[i]);
		    }
		    for(k=0;k<Nptr[5];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[5][k]]);
		      int2ptr = (int ****) d->dataptr;
		      for(u=0; u < d->Ncolumns; u++) {
			(*int2ptr)[lcnum][u][j] = getweightedmean_int(nbin[i],bin_int2ptr[k][u][i],bin_sig[i]);
		      }
		    }
		    for(k=0;k<Nptr[10];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[10][k]]);
		      floatptr = (float ***) d->dataptr;
		      (*floatptr)[lcnum][j] = getweightedmean_float(nbin[i],bin_floatptr[k][i],bin_sig[i]);
		    }
		    for(k=0;k<Nptr[11];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[11][k]]);
		      float2ptr = (float ****) d->dataptr;
		      for(u=0; u < d->Ncolumns; u++) {
			(*float2ptr)[lcnum][u][j] = getweightedmean_float(nbin[i],bin_float2ptr[k][u][i],bin_sig[i]);
		      }
		    }
		    for(k=0;k<Nptr[12];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[12][k]]);
		      longptr = (long ***) d->dataptr;
		      (*longptr)[lcnum][j] = getweightedmean_long(nbin[i],bin_longptr[k][i],bin_sig[i]);
		    }
		    for(k=0;k<Nptr[13];k++) {
		      d = &(p->DataFromLightCurve[ptrindx[13][k]]);
		      long2ptr = (long ****) d->dataptr;
		      for(u=0; u < d->Ncolumns; u++) {
			(*long2ptr)[lcnum][u][j] = getweightedmean_long(nbin[i],bin_long2ptr[k][u][i],bin_sig[i]);
		      }
		    }
		  }

		  var1 = 0.;
		  var2 = 0.;
		  for(k=0;k<nbin[i];k++)
		    {
		      var1 += bin_mag[i][k]/bin_sig[i][k]/bin_sig[i][k];
		      var2 += 1./bin_sig[i][k]/bin_sig[i][k];
		    }
		  mag[j] = var1 / var2;
		  sig[j] = sqrt(1./var2);
		}
	      if(otherdata && (Nptr[6] || Nptr[7] || Nptr[8] || Nptr[9])) {
		indx = startindx[i] + nbin[i]/2;
		for(k=0;k<Nptr[6];k++) {
		  d = &(p->DataFromLightCurve[ptrindx[6][k]]);
		  charptr = (char ***) d->dataptr;
		  (*charptr)[lcnum][j] = (*charptr)[lcnum][indx];
		}
		for(k=0;k<Nptr[7];k++) {
		  d = &(p->DataFromLightCurve[ptrindx[7][k]]);
		  char2ptr = (char ****) d->dataptr;
		  for(u=0; u < d->Ncolumns; u++) {
		    (*char2ptr)[lcnum][u][j] = (*char2ptr)[lcnum][u][indx];
		  }
		}
		for(k=0;k<Nptr[8];k++) {
		  d = &(p->DataFromLightCurve[ptrindx[8][k]]);
		  stringptr = (char ****) d->dataptr;
		  sprintf((*stringptr)[lcnum][j],"%s",(*stringptr)[lcnum][indx]);
		}
		for(k=0;k<Nptr[9];k++) {
		  d = &(p->DataFromLightCurve[ptrindx[9][k]]);
		  string2ptr = (char *****) d->dataptr;
		  for(u=0; u < d->Ncolumns; u++) {
		    sprintf((*string2ptr)[lcnum][u][j],"%s",(*string2ptr)[lcnum][u][indx]);
		  }
		}
	      }
	      if(!tflag)
		{
		  t[j] = t0 + (((double)i) + 0.5)*binsize;
		}
	      else if(tflag == 1)
		{
		  t[j] = getmean(nbin[i], bin_time[i]);
		}
	      else if(tflag == 2)
		{
		  t[j] = median(nbin[i], bin_time[i]);
		}
	      j++;
	    }
	}
      *N = j;
    }
  else
    *N = 0;
  if(otherdata) {

    for(k=0;k<Nptr[0];k++) {
      for(i=0; i < size1; i++) {
	if(n > 0)
	  free(bin_dblptr[k][i]);
      }
      if(Nptr[0] > 0)
	free(bin_dblptr[k]);
    }
    if(bin_dblptr != NULL)
      free(bin_dblptr);
    
    for(k=0;k<Nptr[1];k++) {
      d = &(p->DataFromLightCurve[ptrindx[1][k]]);
      for(u=0; u < d->Ncolumns; u++) {
	for(i=0; i < size1; i++) {
	  if(n > 0)
	    free(bin_dbl2ptr[k][u][i]);
	}
	if(size1 > 0)
	  free(bin_dbl2ptr[k][u]);
      }
      if(d->Ncolumns > 0)
	free(bin_dbl2ptr[k]);
    }
    if(bin_dbl2ptr != NULL)
      free(bin_dbl2ptr);
    
    for(k=0;k<Nptr[2];k++) {
      for(i=0; i < size1; i++) {
	if(n > 0)
	  free(bin_shortptr[k][i]);
      }
      if(size1 > 0)
	free(bin_shortptr[k]);
    }
    if(bin_shortptr != NULL)
      free(bin_shortptr);
    
    for(k=0;k<Nptr[3];k++) {
      d = &(p->DataFromLightCurve[ptrindx[3][k]]);
      for(u=0; u < d->Ncolumns; u++) {
	for(i=0; i < size1; i++) {
	  if(n > 0)
	    free(bin_short2ptr[k][u][i]);
	}
	if(size1 > 0)
	  free(bin_short2ptr[k][u]);
      }
      if(d->Ncolumns > 0)
	free(bin_short2ptr[k]);
    }
    if(bin_short2ptr != NULL)
      free(bin_short2ptr);

    for(k=0;k<Nptr[4];k++) {
      for(i=0; i < size1; i++) {
	if(n > 0)
	  free(bin_intptr[k][i]);
      }
      if(size1 > 0)
	free(bin_intptr[k]);
    }
    if(bin_intptr != NULL)
      free(bin_intptr);
    
    for(k=0;k<Nptr[5];k++) {
      d = &(p->DataFromLightCurve[ptrindx[5][k]]);
      for(u=0; u < d->Ncolumns; u++) {
	for(i=0; i < size1; i++) {
	  if(n > 0)
	    free(bin_int2ptr[k][u][i]);
	}
	if(size1 > 0)
	  free(bin_int2ptr[k][u]);
      }
      if(d->Ncolumns > 0)
	free(bin_int2ptr[k]);
    }
    if(bin_int2ptr != NULL)
      free(bin_int2ptr);
    
    for(k=0;k<Nptr[10];k++) {
      for(i=0; i < size1; i++) {
	if(n > 0)
	  free(bin_floatptr[k][i]);
      }
      if(size1 > 0)
	free(bin_floatptr[k]);
    }
    if(bin_floatptr != NULL)
      free(bin_floatptr);
    
    for(k=0;k<Nptr[11];k++) {
      d = &(p->DataFromLightCurve[ptrindx[11][k]]);
      for(u=0; u < d->Ncolumns; u++) {
	for(i=0; i < size1; i++) {
	  if(n > 0)
	    free(bin_float2ptr[k][u][i]);
	}
	if(size1 > 0)
	  free(bin_float2ptr[k][u]);
      }
      if(d->Ncolumns > 0)
	free(bin_float2ptr[k]);
    }
    if(bin_float2ptr != NULL)
      free(bin_float2ptr);

    for(k=0;k<Nptr[12];k++) {
      for(i=0; i < size1; i++) {
	if(n > 0)
	  free(bin_longptr[k][i]);
      }
      if(size1 > 0)
	free(bin_longptr[k]);
    }
    if(bin_longptr != NULL)
      free(bin_longptr);
    
    for(k=0;k<Nptr[13];k++) {
      d = &(p->DataFromLightCurve[ptrindx[13][k]]);
      for(u=0; u < d->Ncolumns; u++) {
	for(i=0; i < size1; i++) {
	  if(n > 0)
	    free(bin_long2ptr[k][u][i]);
	}
	if(size1 > 0)
	  free(bin_long2ptr[k][u]);
      }
      if(d->Ncolumns > 0)
	free(bin_long2ptr[k]);
    }
    if(bin_long2ptr != NULL)
      free(bin_long2ptr);

    for(i=0; i < 14; i++) {
      if(Nptr[i])
	free(ptrindx[i]);
    }
  }
#ifdef PARALLEL
  if(bin_mag != NULL) {
    for(i=0; i < size1; i++)
      free(bin_mag[i]);
    free(bin_mag);
  }
  if(bin_sig != NULL) {
    for(i=0; i < size1; i++)
      free(bin_sig[i]);
    free(bin_sig);
  }
  if(bin_time != NULL) {
    for(i=0; i < size1; i++)
      free(bin_time[i]);
    free(bin_time);
  }
  if(nbin != NULL) free(nbin);
  if(startindx != NULL) free(startindx);
#endif


}


	    

