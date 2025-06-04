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
           ((vec) = malloc(sizeof(typ))) == NULL || \
	   (stats_type[(indx)] = (int *) malloc(sizeof(int))) == NULL || \
           (pctval[(indx)] = (double *) malloc(sizeof(double))) == NULL) \
        error(ERR_MEMALLOC); \
    } else { \
        if((ptrindx[(indx)] = (int *) realloc(ptrindx[(indx)],(Nptr[(indx)]+1)*sizeof(int))) == NULL || \
           ((vec) = realloc(vec, (Nptr[(indx)]+1)*sizeof(typ))) == NULL || \
	   (stats_type[(indx)] = (int *) realloc(stats_type[(indx)],(Nptr[(indx)]+1)*sizeof(int))) == NULL || \
           (pctval[(indx)] = (double *) realloc(pctval[(indx)],(Nptr[(indx)]+1)*sizeof(double))) == NULL ) \
        error(ERR_MEMALLOC); \
    } \
    ptrindx[(indx)][Nptr[(indx)]] = k; \
    stats_type[(indx)][Nptr[(indx)]] = default_stats_type; \
    Nptr[(indx)]++; \
  } while(0)

#define GROWPTRVECTOR(indx) do { \
    if(!Nptr[(indx)]) { \
      if((ptrindx[(indx)] = (int *) malloc(sizeof(int))) == NULL || \
	 (stats_type[(indx)] = (int *) malloc(sizeof(int))) == NULL || \
	 (pctval[(indx)] = (double *) malloc(sizeof(double))) == NULL) \
        error(ERR_MEMALLOC); \
    } else { \
      if((ptrindx[(indx)] = (int *) realloc(ptrindx[(indx)],(Nptr[(indx)]+1)*sizeof(int))) == NULL || \
	 (stats_type[(indx)] = (int *) realloc(stats_type[(indx)],(Nptr[(indx)]+1)*sizeof(int))) == NULL || \
         (pctval[(indx)] = (double *) realloc(pctval[(indx)],(Nptr[(indx)]+1)*sizeof(double))) == NULL) \
        error(ERR_MEMALLOC);						\
    } \
    ptrindx[(indx)][Nptr[(indx)]] = k; \
    stats_type[(indx)][Nptr[(indx)]] = default_stats_type; \
    Nptr[(indx)]++; \
  } while(0)


#define GETBINSTATS1(indx,typ,vec1,vec2) do { \
  d = &(p->DataFromLightCurve[ptrindx[indx][k]]); \
  vec1 = (typ ***) d->dataptr; \
  switch(stats_type[indx][k]) { \
    case VARTOOLS_STATSTYPE_MEAN: \
      (*vec1)[lcnum][j] = getmean_##typ(nbin[i],vec2[k][i]); \
      break; \
    case VARTOOLS_STATSTYPE_WEIGHTEDMEAN: \
      (*vec1)[lcnum][j] = getweightedmean_##typ(nbin[i],vec2[k][i],bin_sig[i]); \
      break; \
    case VARTOOLS_STATSTYPE_MEDIAN: \
      (*vec1)[lcnum][j] = median_##typ(nbin[i],vec2[k][i]); \
      break; \
    case VARTOOLS_STATSTYPE_MEDIAN_WEIGHT: \
      (*vec1)[lcnum][j] = median_weight_##typ(nbin[i],vec2[k][i],bin_sig[i]); \
      break; \
    case VARTOOLS_STATSTYPE_STDDEV: \
      (*vec1)[lcnum][j] = stddev_##typ(nbin[i],vec2[k][i]); \
      break; \
    case VARTOOLS_STATSTYPE_MEDDEV: \
      (*vec1)[lcnum][j] = meddev_##typ(nbin[i],vec2[k][i]); \
      break; \
    case VARTOOLS_STATSTYPE_MEDMEDDEV: \
      (*vec1)[lcnum][j] = medmeddev_##typ(nbin[i],vec2[k][i]); \
      break; \
    case VARTOOLS_STATSTYPE_MAD: \
      (*vec1)[lcnum][j] = MAD_##typ(nbin[i],vec2[k][i]); \
      break; \
    case VARTOOLS_STATSTYPE_KURTOSIS: \
      (*vec1)[lcnum][j] = kurtosis_##typ(nbin[i],vec2[k][i]); \
      break; \
    case VARTOOLS_STATSTYPE_SKEWNESS: \
      (*vec1)[lcnum][j] = skewness_##typ(nbin[i],vec2[k][i]); \
      break; \
    case VARTOOLS_STATSTYPE_PERCENTILE: \
      (*vec1)[lcnum][j] = percentile_##typ(nbin[i],vec2[k][i],pctval[indx][k]); \
      break; \
    case VARTOOLS_STATSTYPE_PERCENTILE_WEIGHT: \
      (*vec1)[lcnum][j] = percentile_weight_##typ(nbin[i],vec2[k][i],bin_sig[i],pctval[indx][k]); \
      break; \
    case VARTOOLS_STATSTYPE_MAXIMUM: \
      (*vec1)[lcnum][j] = getmaximum_##typ(nbin[i],vec2[k][i]); \
      break; \
    case VARTOOLS_STATSTYPE_MINIMUM: \
      (*vec1)[lcnum][j] = getminimum_##typ(nbin[i],vec2[k][i]); \
      break; \
    case VARTOOLS_STATSTYPE_SUM: \
      (*vec1)[lcnum][j] = getsum_##typ(nbin[i],vec2[k][i]); \
      break; \
    default: \
      error(ERR_CODEERROR); \
    } \
    if(tflag == VARTOOLS_BINLC_TIMETYPE_NOSHRINK) { \
        for(j2=1; j2 < nbin[i]; j2++) { \
	  (*vec1)[lcnum][j+j2] = (*vec1)[lcnum][j]; \
        } \
    } \
    } while(0)


#define GETBINSTATS2(indx,typ,vec1,vec2) do { \
  d = &(p->DataFromLightCurve[ptrindx[indx][k]]); \
  vec1 = (typ ****) d->dataptr; \
  switch(stats_type[indx][k]) { \
    case VARTOOLS_STATSTYPE_MEAN: \
      for(u=0; u < d->Ncolumns; u++) { \
        (*vec1)[lcnum][u][j] = getmean_##typ (nbin[i],vec2[k][u][i]); \
      } \
      break; \
    case VARTOOLS_STATSTYPE_WEIGHTEDMEAN: \
      for(u=0; u < d->Ncolumns; u++) { \
        (*vec1)[lcnum][u][j] = getweightedmean_##typ (nbin[i],vec2[k][u][i], bin_sig[i]); \
      } \
      break; \
    case VARTOOLS_STATSTYPE_MEDIAN: \
      for(u=0; u < d->Ncolumns; u++) { \
        (*vec1)[lcnum][u][j] = median_##typ (nbin[i],vec2[k][u][i]); \
      } \
      break; \
    case VARTOOLS_STATSTYPE_MEDIAN_WEIGHT: \
      for(u=0; u < d->Ncolumns; u++) { \
        (*vec1)[lcnum][u][j] = median_weight_##typ (nbin[i],vec2[k][u][i],bin_sig[i]); \
      } \
      break; \
    case VARTOOLS_STATSTYPE_STDDEV: \
      for(u=0; u < d->Ncolumns; u++) { \
        (*vec1)[lcnum][u][j] = stddev_##typ (nbin[i],vec2[k][u][i]); \
      } \
      break; \
    case VARTOOLS_STATSTYPE_MEDDEV: \
      for(u=0; u < d->Ncolumns; u++) { \
        (*vec1)[lcnum][u][j] = meddev_##typ (nbin[i],vec2[k][u][i]); \
      } \
      break; \
    case VARTOOLS_STATSTYPE_MEDMEDDEV: \
      for(u=0; u < d->Ncolumns; u++) { \
        (*vec1)[lcnum][u][j] = medmeddev_##typ (nbin[i],vec2[k][u][i]); \
      } \
      break; \
    case VARTOOLS_STATSTYPE_MAD: \
      for(u=0; u < d->Ncolumns; u++) { \
        (*vec1)[lcnum][u][j] = MAD_##typ (nbin[i],vec2[k][u][i]); \
      } \
      break; \
    case VARTOOLS_STATSTYPE_KURTOSIS: \
      for(u=0; u < d->Ncolumns; u++) { \
        (*vec1)[lcnum][u][j] = kurtosis_##typ (nbin[i],vec2[k][u][i]); \
      } \
      break; \
    case VARTOOLS_STATSTYPE_SKEWNESS: \
      for(u=0; u < d->Ncolumns; u++) { \
        (*vec1)[lcnum][u][j] = skewness_##typ (nbin[i],vec2[k][u][i]); \
      } \
      break; \
    case VARTOOLS_STATSTYPE_PERCENTILE: \
      for(u=0; u < d->Ncolumns; u++) { \
        (*vec1)[lcnum][u][j] = percentile_##typ (nbin[i],vec2[k][u][i],pctval[indx][k]); \
      } \
      break; \
    case VARTOOLS_STATSTYPE_PERCENTILE_WEIGHT: \
      for(u=0; u < d->Ncolumns; u++) { \
        (*vec1)[lcnum][u][j] = percentile_weight_##typ (nbin[i],vec2[k][u][i],bin_sig[i],pctval[indx][k]); \
      } \
      break; \
    case VARTOOLS_STATSTYPE_MAXIMUM: \
      for(u=0; u < d->Ncolumns; u++) { \
        (*vec1)[lcnum][u][j] = getmaximum_##typ (nbin[i],vec2[k][u][i]); \
      } \
      break; \
    case VARTOOLS_STATSTYPE_MINIMUM: \
      for(u=0; u < d->Ncolumns; u++) { \
        (*vec1)[lcnum][u][j] = getminimum_##typ (nbin[i],vec2[k][u][i]); \
      } \
      break; \
    case VARTOOLS_STATSTYPE_SUM: \
      for(u=0; u < d->Ncolumns; u++) { \
        (*vec1)[lcnum][u][j] = getsum_##typ (nbin[i],vec2[k][u][i]); \
      } \
      break; \
    default: \
      error(ERR_CODEERROR); \
    } \
    if(tflag == VARTOOLS_BINLC_TIMETYPE_NOSHRINK) { \
      for(u=0; u < d->Ncolumns; u++) { \
        for(j2=1; j2 < nbin[i]; j2++) { \
	  (*vec1)[lcnum][u][j+j2] = (*vec1)[lcnum][u][j]; \
        } \
      } \
    } \
    } while(0)

int binlc_parsevarstring(_Binlc *c) {
  int j, i1, i2, k;
  char *tmpstring = NULL;
  int lentmpstring = 0;
  double pctval;
  c->Nvar = 1;
  j = 0;
  while(c->binvarstring[j] != '\0') {
    if(c->binvarstring[j] == ',')
      c->Nvar += 1;
    j++;
  }
  if((c->binvarnames = (char **) malloc(c->Nvar * sizeof(char *))) == NULL ||
     (c->binstats = (int *) malloc(c->Nvar * sizeof(int))) == NULL ||
     (c->binvars = (_Variable **) malloc(c->Nvar * sizeof(_Variable *))) == NULL ||
     (c->pctval = (double *) malloc(c->Nvar * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  i1 = 0;
  i2 = 0;
  for(k=0; k < c->Nvar; k++) {
    while(c->binvarstring[i2] != '\0' && c->binvarstring[i2] != ',' && c->binvarstring[i2] != ':') {
      i2++;
    }
    if(i2 <= i1) {
      return 1;
    }
    if((c->binvarnames[k] = (char *) malloc((i2 - i1 + 1)*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
    for(j=i1; j < i2; j++) {
      c->binvarnames[k][j-i1] = c->binvarstring[j];
    }
    c->binvarnames[k][j-i1] = '\0';
    i1 = i2+1;
    if(c->binvarstring[i2] == ':') {
      i2 = i2+1;
      while(c->binvarstring[i2] != '\0' && c->binvarstring[i2] != ',') {
	i2++;
      }
      if(i2 <= i1) {
	return 1;
      }
      if(i2 - i1 + 1 > lentmpstring) {
	if(!lentmpstring) {
	  lentmpstring = i2 - i1 + 1;
	  if((tmpstring = (char *) malloc(lentmpstring * sizeof(char))) == NULL)
	    error(ERR_MEMALLOC);
	} else {
	  lentmpstring = i2 - i1 + 1;
	  if((tmpstring = (char *) realloc(tmpstring, lentmpstring * sizeof(char))) == NULL)
	    error(ERR_MEMALLOC);
	}
      }
      for(j=i1; j < i2; j++) {
	tmpstring[j-i1] = c->binvarstring[j];
      }
      tmpstring[j-i1] = '\0';
      i1 = i2+1;
      i2 = i2+1;
      if(!strcmp(tmpstring,"mean")) {
	c->binstats[k] = VARTOOLS_STATSTYPE_MEAN;
      }
      else if(!strcmp(tmpstring,"weightedmean")) {
	c->binstats[k] = VARTOOLS_STATSTYPE_WEIGHTEDMEAN;
      }
      else if(!strcmp(tmpstring,"median")) {
	c->binstats[k] = VARTOOLS_STATSTYPE_MEDIAN;
      }
      else if(!strcmp(tmpstring,"wmedian")) {
	c->binstats[k] = VARTOOLS_STATSTYPE_MEDIAN_WEIGHT;
      }
      else if(!strcmp(tmpstring,"stddev")) {
	c->binstats[k] = VARTOOLS_STATSTYPE_STDDEV;
      }
      else if(!strcmp(tmpstring,"meddev")) {
	c->binstats[k] = VARTOOLS_STATSTYPE_MEDDEV;
      }
      else if(!strcmp(tmpstring,"medmeddev")) {
	c->binstats[k] = VARTOOLS_STATSTYPE_MEDMEDDEV;
      }
      else if(!strcmp(tmpstring,"MAD")) {
	c->binstats[k] = VARTOOLS_STATSTYPE_MAD;
      }
      else if(!strcmp(tmpstring,"kurtosis")) {
	c->binstats[k] = VARTOOLS_STATSTYPE_KURTOSIS;
      }
      else if(!strcmp(tmpstring,"skewness")) {
	c->binstats[k] = VARTOOLS_STATSTYPE_SKEWNESS;
      }
      else if(sscanf(tmpstring,"wpct%lf",&pctval) == 1) {
	c->binstats[k] = VARTOOLS_STATSTYPE_PERCENTILE_WEIGHT;
	c->pctval[k] = pctval;
      }
      else if(sscanf(tmpstring,"pct%lf",&pctval) == 1) {
	c->binstats[k] = VARTOOLS_STATSTYPE_PERCENTILE;
	c->pctval[k] = pctval;
      }
      else if(!strcmp(tmpstring,"max")) {
	c->binstats[k] = VARTOOLS_STATSTYPE_MAXIMUM;
      }
      else if(!strcmp(tmpstring,"min")) {
	c->binstats[k] = VARTOOLS_STATSTYPE_MINIMUM;
      }
      else if(!strcmp(tmpstring,"sum")) {
	c->binstats[k] = VARTOOLS_STATSTYPE_SUM;
      }
      else {
	error2(ERR_INVALIDSTATISTIC,tmpstring);
      }
      
    }
    else {
      if(c->medflag == VARTOOLS_BINLC_BINTYPE_AVERAGE) {
	c->binstats[k] = VARTOOLS_STATSTYPE_MEAN;
      }
      else if(c->medflag == VARTOOLS_BINLC_BINTYPE_MEDIAN) {
	c->binstats[k] = VARTOOLS_STATSTYPE_MEDIAN;
      }
      else if(c->medflag == VARTOOLS_BINLC_BINTYPE_WEIGHTEDAVERAGE) {
	c->binstats[k] = VARTOOLS_STATSTYPE_WEIGHTEDMEAN;
      }
      else
	error(ERR_CODEERROR);
      i2 = i2+1;
    }
  }
  if(tmpstring != NULL)
    free(tmpstring);
  return(0);
}

void binlc(ProgramData *p, _Binlc *c, int lcnum, int lclistnum)
{
  /* N, t, mag and sig define the input light curve, the output binned light curve is stored in these variables as well.

We assume the input light curve is sorted in time.

medflag - 0 means average bin, 1 means median bin, 2 means weighted average bin

binsize_Nbins_flag - 0 means that the binsize (in units of time) is specified, 1 means the number of bins is specified.

firstbinflag - 0 means that the first bin starts at t[0], 1 means the first bin starts at t[0] - firstbin/binsize; This is only relevant if T0source = -1

T0source - -1 means that firstbinflag settings apply (bin starts at t[0] or t[0] - firstbin/binsize);
           PERTYPE_FIXCOLUMN - it comes from a previously computed output column
           PERTYPE_FIX - it is specified on the command line
           PERTYPE_SPECIFIED - it comes from the input list
           PERTYPE_EXPR - it is computed from an analytic expression

tflag - 0 means the output time for each bin is the time at the center of the bin, 1 means to take the average of the times of points that fall within the bin, 2 means to take the median.
  */

  int *N;
  double *t, *mag, *sig;

  int i, j, j2, k, n, u, otherdata = 0;
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

  double t0, t1, var1, var2, t0lc;
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
  int* stats_type[14];
  double *pctval[14];
  int ptr_datasize1[14];
  int ptr_datasize2[14];
  int ptr_datasize3[14];
  
  _DataFromLightCurve *d;

  int medflag, binsize_Nbins_flag, Nbins, firstbinflag, tflag;
  int T0source;
  double binsize, firstbin;
  int default_stats_type;
  int indx2;

  int testflag;

  medflag = c->medflag; 
  binsize_Nbins_flag = c->binsize_Nbins_flag;
  Nbins = c->Nbins;
  firstbinflag = c->firstbinflag;
  tflag = c->tflag;
  binsize = c->binsize;
  firstbin = c->firstbin;
  T0source = c->T0source;

  switch(medflag) {
  case VARTOOLS_BINLC_BINTYPE_AVERAGE:
    default_stats_type = VARTOOLS_STATSTYPE_MEAN;
    break;
  case VARTOOLS_BINLC_BINTYPE_MEDIAN:
    default_stats_type = VARTOOLS_STATSTYPE_MEDIAN;
    break;
  case VARTOOLS_BINLC_BINTYPE_WEIGHTEDAVERAGE:
    default_stats_type = VARTOOLS_STATSTYPE_WEIGHTEDMEAN;
    break;
  default:
    default_stats_type = VARTOOLS_STATSTYPE_MEAN;
    break;
  }

  for(i=0; i < 14; i++) Nptr[i] = 0;

  N = &p->NJD[lcnum];
  t = p->t[lcnum];
  mag = p->mag[lcnum];
  sig = p->sig[lcnum];

  n = *N;
  
  if(n <= 0)
    return;

  /* First determine the number of bins that we need to allocate memory for */
  t1 = t[n - 1];  
  t0lc = t[0];

  if(T0source == -1) {
    
    t0 = t0lc;
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
  } else {
    if(T0source == PERTYPE_FIXCOLUMN) {
      getoutcolumnvalue(c->t0_linkedcolumn, lcnum, lclistnum, VARTOOLS_TYPE_DOUBLE, &t0);
    } else if(T0source == PERTYPE_FIX) {
      t0 = c->t0fixval;
    } else if(T0source == PERTYPE_SPECIFIED) {
      t0 = c->t0listval[lclistnum];
    } else if(T0source == PERTYPE_EXPR) {
      t0 = EvaluateExpression(lclistnum, lcnum, 0, c->t0expr);
    }
    if(t1 < t0) {
      *N = 0;
      return;
    }

    if(binsize_Nbins_flag)
      {
	binsize = (t1 - t0)/(double) Nbins;
      }
    else
      {
	if(t0 < t0lc - binsize) {
	  t0 = t0 + binsize*floor((t0lc - t0)/binsize);
	}
	Nbins = ceil((t1 - t0)/binsize);
      }
  }

  /* Only run if Nbins > 0 */
  if(Nbins <= 0) {
    *N = 0;
    return;
  }

  /* Figure out what light-curve vectors exist in addition to those
     given above, and prepare storage for sorting them into bins, also
     determine which statistic to calculate */
  for(k=0; k < p->NDataFromLightCurve; k++) {
    d = &(p->DataFromLightCurve[k]);
    Nc = d->Ncolumns;
    /* skip vectors corresponding to t, mag or sig */
    if(Nc == 0 && (d->datatype == VARTOOLS_TYPE_DOUBLE || 
		   d->datatype == VARTOOLS_TYPE_CONVERTJD) 
       && !c->only_bin_columns) {
      if((*((double ***) d->dataptr))[lcnum] == p->t[lcnum] ||
	 (*((double ***) d->dataptr))[lcnum] == p->mag[lcnum] ||
	 (*((double ***) d->dataptr))[lcnum] == p->sig[lcnum])
	continue;
    }
    if(c->only_bin_columns) {
      testflag = 0;
      if(d->variable != NULL) {
	for(i=0; i < c->Nvar; i++) {
	  if(c->binvars[i] == d->variable) {
	    testflag = 1;
	    break;
	  }
	}
      }
      if(!testflag)
	continue;
    }
    otherdata = 1;
    switch(d->datatype)
      {
      case VARTOOLS_TYPE_DOUBLE:
	if(Nc == 0) {
	  indx2 = 0;
	  GROWBINVECTOR(0,double **,bin_dblptr);
	} else {
	  indx2 = 1;
	  GROWBINVECTOR(1, double ***, bin_dbl2ptr);
	}
	break;
      case VARTOOLS_TYPE_CONVERTJD:
	if(Nc == 0) {
	  indx2 = 0;
	  GROWBINVECTOR(0,double **,bin_dblptr);
	} else {
	  indx2 = 1;
	  GROWBINVECTOR(1,double ***,bin_dbl2ptr);
	}
	break;
      case VARTOOLS_TYPE_SHORT:
	if(Nc == 0) {
	  indx2 = 2;
	  GROWBINVECTOR(2,short **,bin_shortptr);
	} else {
	  indx2 = 3;
	  GROWBINVECTOR(3,short ***,bin_short2ptr);
	}
	break;
      case VARTOOLS_TYPE_INT:
	if(Nc == 0) {
	  indx2 = 4;
	  GROWBINVECTOR(4,int **,bin_intptr);
	} else {
	  indx2 = 5;
	  GROWBINVECTOR(5,int ***,bin_int2ptr);
	}
	break;
      case VARTOOLS_TYPE_CHAR:
	
	if(Nc == 0) {
	  indx2 = 6;
	  GROWPTRVECTOR(6);
	} else {
	  indx2 = 7;
	  GROWPTRVECTOR(7);
	}
	break;
      case VARTOOLS_TYPE_STRING:
	if(Nc == 0) {
	  indx2 = 8;
	  GROWPTRVECTOR(8);
	} else {
	  indx2 = 9;
	  GROWPTRVECTOR(9);
	}
	break;
      case VARTOOLS_TYPE_FLOAT:
	if(Nc == 0) {
	  indx2 = 10;
	  GROWBINVECTOR(10,float **,bin_floatptr);
	} else {
	  indx2 = 11;
	  GROWBINVECTOR(11,float ***,bin_float2ptr);
	}
	break;
      case VARTOOLS_TYPE_LONG:
	if(Nc == 0) {
	  indx2 = 12;
	  GROWBINVECTOR(12,long **,bin_longptr);
	} else {
	  indx2 = 13;
	  GROWBINVECTOR(13,long ***,bin_long2ptr);
	}
	break;
      default:
	error(ERR_BADTYPE);
	break;
      }
    if(d->variable != NULL) {
      for(i=0; i < c->Nvar; i++) {
	if(c->binvars[i] == d->variable) {
	  stats_type[indx2][Nptr[indx2]-1] = c->binstats[i];
	  if(c->binstats[i] == VARTOOLS_STATSTYPE_PERCENTILE_WEIGHT ||
	     c->binstats[i] == VARTOOLS_STATSTYPE_PERCENTILE) {
	    pctval[indx2][Nptr[indx2]-1] = c->pctval[i];
	  }
	  break;
	}
      }
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
      if(!c->usemask) {
	/* We are not masking any points to fill into the bins */
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
      } else {
	/* We are masking points to fill into the bins */
	for(i=0; i<n; i++)
	  {
	    if(EvaluateVariable_Double(lclistnum, lcnum, i, c->maskvar) <= VARTOOLS_MASK_TINY)
	      continue;
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
      }
      
      /* Now fill out the new light curve */
      for(i=0, j=0; i < Nbins; i++)
	{
	  if(nbin[i] > 0)
	    {
	      if(otherdata) {
		for(k=0;k<Nptr[0];k++) {
		  GETBINSTATS1(0,double,dblptr,bin_dblptr);
		}
		for(k=0;k<Nptr[1];k++) {
		  GETBINSTATS2(1,double,dbl2ptr,bin_dbl2ptr);
		}
		for(k=0;k<Nptr[2];k++) {
		  GETBINSTATS1(2,short,shortptr,bin_shortptr);
		}
		for(k=0;k<Nptr[3];k++) {
		  GETBINSTATS2(3,short,short2ptr,bin_short2ptr);
		}
		for(k=0;k<Nptr[4];k++) {
		  GETBINSTATS1(4,int,intptr,bin_intptr);
		}
		for(k=0;k<Nptr[5];k++) {
		  GETBINSTATS2(5,int,int2ptr,bin_int2ptr);
		}
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
		for(k=0;k<Nptr[10];k++) {
		  GETBINSTATS1(10,float,floatptr,bin_floatptr);
		}
		for(k=0;k<Nptr[11];k++) {
		  GETBINSTATS2(11,float,float2ptr,bin_float2ptr);
		}
		for(k=0;k<Nptr[12];k++) {
		  GETBINSTATS1(12,long,longptr,bin_longptr);
		}
		for(k=0;k<Nptr[13];k++) {
		  GETBINSTATS2(13,long,long2ptr,bin_long2ptr);
		}
	      }
	      if(medflag == VARTOOLS_BINLC_BINTYPE_AVERAGE &&
		 !c->only_bin_columns)
		{
		  mag[j] = getmean(nbin[i],bin_mag[i]);
		  var1 = 0.;
		  for(k=0;k<nbin[i];k++)
		    {
		      var1 += bin_sig[i][k]*bin_sig[i][k];
		    }
		  sig[j] = sqrt(var1 / (double) (nbin[i]*nbin[i]));
		}
	      else if(medflag == VARTOOLS_BINLC_BINTYPE_MEDIAN &&
		      !c->only_bin_columns)
		{
		  mag[j] = median(nbin[i],bin_mag[i]);
		  var1 = 0.;
		  for(k=0;k<nbin[i];k++)
		    {
		      var1 += bin_sig[i][k]*bin_sig[i][k];
		    }
		  sig[j] = 1.253 * sqrt(var1 / (double) (nbin[i]*nbin[i]));
		}
	      else if(medflag == VARTOOLS_BINLC_BINTYPE_WEIGHTEDAVERAGE &&
		      !c->only_bin_columns)
		{
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
	      if(tflag == VARTOOLS_BINLC_TIMETYPE_CENTER)
		{
		  t[j] = t0 + (((double)i) + 0.5)*binsize;
		}
	      else if(tflag == VARTOOLS_BINLC_TIMETYPE_AVERAGE)
		{
		  t[j] = getmean(nbin[i], bin_time[i]);
		}
	      else if(tflag == VARTOOLS_BINLC_TIMETYPE_MEDIAN)
		{
		  t[j] = median(nbin[i], bin_time[i]);
		}
	      if(tflag == VARTOOLS_BINLC_TIMETYPE_NOSHRINK &&
		 !c->only_bin_columns) {
		for(j2=1; j2 < nbin[i]; j2++) {
		  mag[j+j2] = mag[j];
		  sig[j+j2] = sig[j];
		}
	      }
	      if(tflag == VARTOOLS_BINLC_TIMETYPE_NOSHRINK) {
		j += nbin[i];
	      }
	      else {
		j++;
	      }
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
      if(Nptr[i]) {
	free(ptrindx[i]);
	free(stats_type[i]);
	free(pctval[i]);
      }
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

