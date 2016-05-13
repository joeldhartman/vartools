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

/* This file contains the function to perform the sigma-clipping command on light curves for the program vartools by J. Hartman */


void sigclip_copyterms(int i,int j,ProgramData *p,int lc)
{
  int k, Nc, u;
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
  _DataFromLightCurve *d;
  for(k=0;k<p->NDataFromLightCurve;k++)
    {
      d = &(p->DataFromLightCurve[k]);
      Nc = d->Ncolumns;
      if(Nc == 0) {
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = (double ***) d->dataptr;
	  (*dblptr)[lc][j] = (*dblptr)[lc][i];
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  dblptr = (double ***) d->dataptr;
	  (*dblptr)[lc][j] = (*dblptr)[lc][i];
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = (char ****) d->dataptr;
	  sprintf(((*stringptr)[lc][j]),"%s",((*stringptr)[lc][i]));
	  /*memcpy(((*stringptr)[lc][j]),((*stringptr)[lc][i]),d->maxstringlength);*/
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = (int ***) d->dataptr;
	  (*intptr)[lc][j] = (*intptr)[lc][i];
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = (short ***) d->dataptr;
	  (*shortptr)[lc][j] = (*shortptr)[lc][i];
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = (long ***) d->dataptr;
	  (*longptr)[lc][j] = (*longptr)[lc][i];
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = (float ***) d->dataptr;
	  (*floatptr)[lc][j] = (*floatptr)[lc][i];
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = (char ***) d->dataptr;
	  (*charptr)[lc][j] = (*charptr)[lc][i];
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      } else if(Nc > 0) {
	for(u=0; u < Nc; u++) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dbl2ptr = (double ****) d->dataptr;
	    (*dbl2ptr)[lc][u][j] = (*dbl2ptr)[lc][u][i];
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    dbl2ptr = (double ****) d->dataptr;
	    (*dbl2ptr)[lc][u][j] = (*dbl2ptr)[lc][u][i];
	    break;
	  case VARTOOLS_TYPE_STRING:
	    string2ptr = (char *****) d->dataptr;
	    sprintf(((*string2ptr)[lc][u][j]),"%s",((*string2ptr)[lc][u][i]));
	    /*memcpy(((*string2ptr)[lc][u][j]),((*string2ptr)[lc][u][i]),d->maxstringlength);*/
	    break;
	  case VARTOOLS_TYPE_INT:
	    int2ptr = (int ****) d->dataptr;
	    (*int2ptr)[lc][u][j] = (*int2ptr)[lc][u][i];
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    short2ptr = (short ****) d->dataptr;
	    (*short2ptr)[lc][u][j] = (*short2ptr)[lc][u][i];
	    break;
	  case VARTOOLS_TYPE_LONG:
	    long2ptr = (long ****) d->dataptr;
	    (*long2ptr)[lc][u][j] = (*long2ptr)[lc][u][i];
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    float2ptr = (float ****) d->dataptr;
	    (*float2ptr)[lc][u][j] = (*float2ptr)[lc][u][i];
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    char2ptr = (char ****) d->dataptr;
	    (*char2ptr)[lc][u][j] = (*char2ptr)[lc][u][i];
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      } else {
	for(u=0; u < (-Nc); u++) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dbl2ptr = (double ****) d->dataptr;
	    (*dbl2ptr)[lc][j][u] = (*dbl2ptr)[lc][i][u];
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    dbl2ptr = (double ****) d->dataptr;
	    (*dbl2ptr)[lc][j][u] = (*dbl2ptr)[lc][i][u];
	    break;
	  case VARTOOLS_TYPE_STRING:
	    string2ptr = (char *****) d->dataptr;
	    sprintf(((*string2ptr)[lc][j][u]),"%s",((*string2ptr)[lc][i][u]));
	    /*memcpy(((*string2ptr)[lc][j][u]),((*string2ptr)[lc][i][u]),d->maxstringlength);*/
	    break;
	  case VARTOOLS_TYPE_INT:
	    int2ptr = (int ****) d->dataptr;
	    (*int2ptr)[lc][j][u] = (*int2ptr)[lc][i][u];
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    short2ptr = (short ****) d->dataptr;
	    (*short2ptr)[lc][j][u] = (*short2ptr)[lc][i][u];
	    break;
	  case VARTOOLS_TYPE_LONG:
	    long2ptr = (long ****) d->dataptr;
	    (*long2ptr)[lc][j][u] = (*long2ptr)[lc][i][u];
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    float2ptr = (float ****) d->dataptr;
	    (*float2ptr)[lc][j][u] = (*float2ptr)[lc][i][u];
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    char2ptr = (char ****) d->dataptr;
	    (*char2ptr)[lc][j][u] = (*char2ptr)[lc][i][u];
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      }
    }
}


/* Function that performs sigma clipping on a light curve and returns the reduced size of the light curve */
int sigclip (int size, double *t, double *mag, double *sig, double *ave, double *stddev, double *rmsthy, int *ngood, double clip, int iter, int lc, ProgramData *p, int niter, int usemedian)
{
  int i, j, k, clipped, nclipped, lastnclipped;
  double avein, stddevin;
  double *magcpy = NULL;
  int sizemagcpy = 0;

  nclipped = 0;
  if(size > 0 && clip > 0)
    {
      avein = 0.;
      *ngood = 0;
      for(i=0 ; i<size; i++)
	if(sig[i] > 0. && !isnan(mag[i]))
	  {
	    avein += mag[i];
	    (*ngood)++;
	  }
      if((*ngood) > 0)
	{
	  if(usemedian) {
	    if(*ngood == size) {
	      avein = median(size, mag);
	    } else {
	      if(!sizemagcpy) {
		if((magcpy = (double *) malloc((*ngood)*sizeof(double))) == NULL)
		  error(ERR_MEMALLOC);
		sizemagcpy = *ngood;
	      }
	      else if(sizemagcpy < (*ngood)) {
		if((magcpy = (double *) realloc(magcpy, (*ngood)*sizeof(double))) == NULL)
		  error(ERR_MEMALLOC);
		sizemagcpy = *ngood;
	      }
	      for(i=0, k=0; i < size; i++) {
		if(sig[i] > 0. && !isnan(mag[i])) {
		  magcpy[k] = mag[i];
		  k++;
		}
	      }
	      avein = median(k, magcpy);
	    }
	  } else {
	    avein /= (double) (*ngood);
	  }
	  stddevin = 0.;
	  for(i = 0; i < size; i++)
	    if(sig[i] > 0. && !isnan(mag[i]))
	      stddevin += SQR(mag[i] - avein);
	  stddevin = sqrt(stddevin / ((double) (*ngood)));


	  *ave = avein;
	  *stddev = stddevin;

	  lastnclipped = 0;
	  do
	    {
	      i = 0;
	      clipped = 0;
	      *ngood = 0;
	      nclipped = 0;
	      for(i=0;i<size;i++)
		{
		  if(isnan(mag[i]) || sig[i] <= 0. || (!isnan(mag[i]) && (fabs(mag[i] - (*ave)) > ((*stddev) * clip) || sig[i] <= 0.)))
		    {
		      nclipped++;
		      mag[i] = sqrt(-1);
		    }
		  else
		    (*ngood)++;
		}

	      if((*ngood) > 0)
		{
		  i = 0;
		  while(isnan(mag[i]) || sig[i] <= 0.)
		    i++;
		  *ave = mag[i];
		  *rmsthy = sig[i]*sig[i];
		  for(i = i + 1; i < size ; i++)
		    {
		      if(!isnan(mag[i]) && sig[i] > 0.)
			{
			  *ave += mag[i];
			  *rmsthy += sig[i]*sig[i];
			}
		    }
		  *rmsthy = sqrt((*rmsthy) / (*ngood));
		  if(usemedian) {
		    if(*ngood == size) {
		      avein = median(size, mag);
		    } else {
		      if(!sizemagcpy) {
			if((magcpy = (double *) malloc((*ngood)*sizeof(double))) == NULL)
			  error(ERR_MEMALLOC);
			sizemagcpy = *ngood;
		      }
		      else if(sizemagcpy < (*ngood)) {
			if((magcpy = (double *) realloc(magcpy, (*ngood)*sizeof(double))) == NULL)
			  error(ERR_MEMALLOC);
			sizemagcpy = *ngood;
		      }
		      for(i=0, k=0; i < size; i++) {
			if(sig[i] > 0. && !isnan(mag[i])) {
			  magcpy[k] = mag[i];
			  k++;
			}
		      }
		      *ave = median(k, magcpy);
		    }
		  } else {
		    *ave /= (*ngood);
		  }
		  i = 0;
		  while(isnan(mag[i]) || sig[i] <= 0.)
		    i++;
		  *stddev = SQR(mag[i] - *ave);
		  for ( i = i + 1; i < size; i++)
		    {
		      if(!isnan(mag[i]) && sig[i] > 0.)
			{
			  *stddev += SQR(mag[i] - *ave);
			}
		    }
		  *stddev = sqrt((*stddev) / (*ngood));
		}
	      else
		{
		  *ave = -1.;
		  *stddev = -1.;
		  *rmsthy = -1.;
		}
	      if(niter)
		niter--;
	      if(nclipped > lastnclipped)
		clipped = 1;
	      lastnclipped = nclipped;
	    }
	  while (clipped && (iter || niter > 0)) ;
	}
      else
	{
	  size = 0;
	  *ave = -1.;
	  *stddev = -1.;
	  *rmsthy = -1.;
	}
    }
  else if(size > 0) {
    nclipped = 0;
    avein = 0.;
    *ngood = 0;
    for(i=0 ; i<size; i++) {
      if(sig[i] > 0. && !isnan(mag[i]))
	{
	  avein += mag[i];
	  (*ngood)++;
	}
      else
	{
	  nclipped++;
	  mag[i] = sqrt(-1);
	}
    }
    if((*ngood) > 0)
      {
	if(usemedian) {
	  if(*ngood == size) {
	    avein = median(size, mag);
	  } else {
	    if(!sizemagcpy) {
	      if((magcpy = (double *) malloc((*ngood)*sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	      sizemagcpy = *ngood;
	    }
	    else if(sizemagcpy < (*ngood)) {
	      if((magcpy = (double *) realloc(magcpy, (*ngood)*sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	      sizemagcpy = *ngood;
	    }
	    for(i=0, k=0; i < size; i++) {
	      if(sig[i] > 0. && !isnan(mag[i])) {
		magcpy[k] = mag[i];
		k++;
	      }
	    }
	    avein = median(k, magcpy);
	  }
	} else {
	  avein /= (double) (*ngood);
	}
	stddevin = 0.;
	for(i = 0; i < size; i++)
	  if(sig[i] > 0. && !isnan(mag[i]))
	    stddevin += SQR(mag[i] - avein);
	stddevin = sqrt(stddevin / ((double) (*ngood)));

	*ave = avein;
	*stddev = stddevin;
      }
    else {
      *ave = -1.;
      *stddev = -1.;
    }
  }
  else
    {
      *ave = -1.;
      *stddev = -1.;
    }
  /* Now remove the isnan points from the light curve */
  if(nclipped)
    {
      j = 0;
      for(i=0;i<size;i++)
	{
	  if(!isnan(mag[i]))
	    {
	      if(i != j)
		sigclip_copyterms(i,j,p,lc);
	      j++;
	    }
	}
      p->NJD[lc] = j;
      if(p->readimagestring)
	{
	  for(i=0;i<p->NJD[lc];i++)
	    p->stringid_idx[0][i] = i;
	  mysortstringint(p->NJD[lc], MAXIDSTRINGLENGTH, p->stringid[lc], p->stringid_idx[lc]);
	}
    }
  if(magcpy != NULL)
    free(magcpy);
  return(nclipped);
}


