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
/* This is an implementation of the Trend Filtering Algorithm due to Kovacs, Bakos and Noyes 2005.  It takes in a list of light curves and a dates file. We compute the fitting matrix and its singular value decomposition. The template list for each light curve is all other light curves that are a specified number of pixels away. For each light curve that we read in we subtract the average value and remove outliers. We then correct the fitting matrix and recompute its singular value decomposition if the light curve was also used in defining the trends, we then calculate the coefficients for each trend and write out the corrected light curves and statistics.

Compile line:

gcc -o tfa tfa.c -lm

*/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include "commands.h"
#include "programdata.h"
#include "functions.h"

#ifdef USECFITSIO
#include "fitsio.h"
#endif

#define MAXSTRINGLENGTH 2048
/*#define JDTOL 0.00001*/
#define DEFAULT_CLIPPING_VALUE 5.0

#define SIGN_(A,B) ((B) >= 0.0 ? ABS_(A) : (-(ABS_(A))))
#define MAX_(A,B) ((A) > (B) ? (A) : (B))
#define MIN_(A,B) ((A) < (B) ? (A) : (B))
#define TINY 1.0e-20

/* Number of cycles to find the preliminary magnitude corrections */

void zero_lc_averages(int Njd, double *b, double clipping)
{
  int i, n;
  double ave1, ave2, std;
  long double val1, val2;
  val1 = 0.;
  val2 = 0.;
  n = 0;
  for(i=0;i<Njd;i++)
    {
      if(!isnan(b[i]))
	{
	  val1 += b[i]*b[i];
	  val2 += b[i];
	  n++;
	}
    }
  if(n > 0)
    {
      ave1 = (double) (val2 / (long double) n);
      std = sqrt((double)((val1 / (long double) n) - (val2 * val2 / (long double) (n * n))));
    }
  else
    {
      ave1 = 0.;
      std = 0.;
    }
  val1 = 0.;
  val2 = 0.;
  n = 0;
  for(i=0;i<Njd;i++)
    {
      if(!isnan(b[i]))
	if(ABS_(b[i] - ave1) < clipping * std)
	  {
	    val1 += b[i]*b[i];
	    val2 += b[i];
	    n++;
	  }
    }
  if(n > 0)
    ave2 = (double) (val2 / (long double) n);
  else
    ave2 = 0.;
  for(i=0;i<Njd;i++)
    {
      if(isnan(b[i]) || ABS_(b[i] - ave1) >= clipping * std)
	{
	  b[i] = sqrt(-1.);
	}
      else
	{
	  b[i] -= ave2;
	}
    }
}


void zero_trend_averages(int Njd, int Ntrends, double **trends, double clipping)
{
  int i, j, n;
  double ave1, ave2, std;
  long double val1, val2;
  for(j=0;j<Ntrends;j++)
    {
      val1 = 0.;
      val2 = 0.;
      n = 0;
      for(i=0;i<Njd;i++)
	{
	  if(!isnan(trends[i][j]))
	    {
	      val1 += trends[i][j]*trends[i][j];
	      val2 += trends[i][j];
	      n++;
	    }
	}
      if(n > 0)
	{
	  ave1 = (double) (val2 / (long double) n);
	  std = sqrt((double)((val1 / (long double) n) - (val2 * val2 / (long double) (n * n))));
	}
      else
	{
	  ave1 = 0.;
	  std = 0.;
	}

      val1 = 0.;
      val2 = 0.;
      n = 0;
      for(i=0;i<Njd;i++)
	{
	  if(!isnan(trends[i][j]))
	    if(ABS_(trends[i][j] - ave1) < clipping * std)
	      {
		val1 += trends[i][j]*trends[i][j];
		val2 += trends[i][j];
		n++;
	      }
	}
      if(n > 0)
	ave2 = (double) (val2 / (long double) n);
      else
	ave2 = 0.;
      for(i=0;i<Njd;i++)
	{
	  if(isnan(trends[i][j]) || ABS_(trends[i][j] - ave1) >= clipping * std)
	    {
	      trends[i][j] = 0.;
	    }
	  else
	    {
	      trends[i][j] -= ave2;
	    }
	}
    }
}

#ifdef USECFITSIO
int ReadFitsTFATemplate(ProgramData *p, int jdcol_isfromheader, char *jdcol_headername, int *JDcol_trend_ptr, int magcol_isfromheader, char *magcol_headername, int *magcol_trend_ptr, char *trendname, int *sizevec, int *lengthvec, double **magin_tmp, double **jdin_tmp, char ***stringid_tmp, int **stringid_tmpidx)
{
  fitsfile *infile;
  int status = 0, hdunum, ncols, hdutype, anynulallcolumns;
  long nrows;
  int j, jold, k, i, l, N, oldsizesinglelc, anynul;
  int JDcol_trend;
  int magcol_trend;
  char *nullarray, *nullarraystore;
  char tmpchar;


  int Nc, u;

  *lengthvec = 0;

#ifdef PARALLEL
  if(p->Nproc_allow > 1) {
    while(pthread_mutex_trylock(&(p->cfitsio_mutex)));
  }
#endif

  hdutype = 0; status = 0; nrows = 0; ncols = 0;

  if((fits_open_file(&infile,trendname,READONLY,&status)))
    {
      error2_noexit(ERR_CANNOTOPEN,trendname);
      return(ERR_CANNOTOPEN);
    }
  if(fits_get_hdu_num(infile, &hdunum) == 1)
    {
      fits_movabs_hdu(infile, 2, &hdutype, &status);
    }
  else
    fits_get_hdu_type(infile, &hdutype, &status);
  
  if(hdutype == IMAGE_HDU) {
    error2_noexit(ERR_IMAGEHDU,trendname);
    return(ERR_IMAGEHDU);
  }

  fits_get_num_rows(infile, &nrows, &status);
  fits_get_num_cols(infile, &ncols, &status);

  if((nullarray = (char *) calloc(nrows, sizeof(char))) == NULL)
    error(ERR_MEMALLOC);
  if((nullarraystore = (char *) calloc(nrows, sizeof(char))) == NULL)
    error(ERR_MEMALLOC);

  if(jdcol_isfromheader) {
    fits_get_colnum(infile, 0, jdcol_headername, JDcol_trend_ptr, &status);
    if(status == COL_NOT_FOUND) {
      error2_noexit(ERR_MISSING_FITSLC_HEADERNAME, jdcol_headername);
      return(ERR_MISSING_FITSLC_HEADERNAME);
    }
  }
  JDcol_trend = *JDcol_trend_ptr;
  if(magcol_isfromheader) {
    fits_get_colnum(infile, 0, magcol_headername, magcol_trend_ptr, &status);
    if(status == COL_NOT_FOUND) {
      error2_noexit(ERR_MISSING_FITSLC_HEADERNAME, magcol_headername);
      return(ERR_MISSING_FITSLC_HEADERNAME);
    }
  }
  magcol_trend = *magcol_trend_ptr;

  /* Increase the memory for the data vectors if needed */
  if(nrows >= *sizevec) {
    if(*sizevec == 0) {
      *sizevec = nrows;
      if((*magin_tmp = (double *) malloc(*sizevec * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      if(p->matchstringid) {
	if((*stringid_tmp = (char **) malloc(*sizevec * sizeof(char *))) == NULL ||
	   (*stringid_tmpidx = (int *) malloc(*sizevec * sizeof(int))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0;j<*sizevec;j++)
	  {
	    if(((*stringid_tmp)[j] = (char *) malloc(MAXIDSTRINGLENGTH * sizeof(char))) == NULL)
	      error(ERR_MEMALLOC);
	  }
      } else {
	if((*jdin_tmp = (double *) malloc(*sizevec * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
      }
    } else {
      if((*magin_tmp = (double *) realloc(*magin_tmp, nrows*sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      if(p->matchstringid) {
	if((*stringid_tmp = (char **) realloc(*stringid_tmp, nrows * sizeof(char *))) == NULL ||
	   (*stringid_tmpidx = (int *) realloc(*stringid_tmpidx, nrows * sizeof(int))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=*sizevec; j<nrows; j++)
	  {
	    if(((*stringid_tmp)[j] = (char *) malloc(MAXIDSTRINGLENGTH * sizeof(char))) == NULL)
	      error(ERR_MEMALLOC);
	  }
      } else {
	if((*jdin_tmp = (double *) realloc(*jdin_tmp, nrows*sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
      }
      *sizevec = nrows;
    }
  }

  anynulallcolumns = 0;

  /* Read in the Columns */
  status = 0;
  anynul = 0;

  fits_read_colnull(infile, TDOUBLE, magcol_trend, 1, 1, nrows, *magin_tmp, nullarray, &anynul,&status);

  if(p->matchstringid) {
    for(k=1; k<=nrows && !status; k++) {
      fits_read_col_str(infile, JDcol_trend, k, 1, 1, 0, &((*stringid_tmp)[k-1]), &anynul, &status);
      if(anynul)
	{
	  nullarray[k-1] = 1;
	  anynulallcolumns = 1;
	  anynul = 0;
	}
    }
  } else {
    fits_read_colnull(infile, TDOUBLE, JDcol_trend, 1, 1, nrows, *jdin_tmp, nullarray, &anynul,&status);
  }
  if(anynul) {
    anynulallcolumns = 1;
    for(k=0; k < nrows; k++) {
      nullarraystore[k] = nullarraystore[k] || nullarray[k];
    }
  }
  if(status) {
    fits_report_error(stderr, status);
    error(ERR_FITSERROR);
  }
  
  fits_close_file(infile, &status);
  if(status) {
    fits_report_error(stderr, status);
    error(ERR_FITSERROR);
  }
  
#ifdef PARALLEL
  if(p->Nproc_allow > 1) {
    pthread_mutex_unlock(&(p->cfitsio_mutex));
  }
#endif

  if(anynulallcolumns) {
    for(i=0,j=0; i < nrows; i++) {
        if(!nullarraystore[i])
	{
	  if(i != j) {
	    (*magin_tmp)[j] = (*magin_tmp)[i];
	    if(p->matchstringid) {
	      sprintf((*stringid_tmp)[j],"%s",(*stringid_tmp)[i]);
	    } else {
	      (*jdin_tmp)[j] = (*jdin_tmp)[i];
	    }
	  }
	  j++;
	}
    }
    nrows = j;
  }


  if(p->matchstringid)
    {
      for(i=0;i<nrows;i++)
	(*stringid_tmpidx)[i] = i;
    }
  
  *lengthvec = nrows;
  
  free(nullarray);
  free(nullarraystore);
  return 0;
} 
#endif


void initialize_tfa(_TFA *tfa, ProgramData *p)
{
  FILE *dates, *trendin, *trend_list;
  double *magin_tmp = NULL;
  double *jdin_tmp = NULL;
  char **stringid_tmp = NULL;
  int *stringid_tmpidx = NULL;
  int lengthtmp, sizetmpvec = 0;
  double jdin, jdin_last, magin;
  char stringidin[MAXIDSTRINGLENGTH];
  void *ptr1, *ptr2;
  int i, j, k, lcline, lccol, lcindx, lcindx_old, col1, col2, type1, type2, ii, jj, kk;
  void error2(int, char *);
  char trend_name[MAXSTRINGLENGTH], dums[MAXSTRINGLENGTH];
  char *line;
  size_t line_size = MAXSTRINGLENGTH;
  int Nthreads;

  line = malloc(line_size);

#ifdef PARALLEL
  Nthreads = p->Nproc_allow;
#else
  Nthreads = 1;
#endif

  if(tfa->JDcol_trend <= tfa->magcol_trend)
    {
      col1 = tfa->JDcol_trend;
      col2 = tfa->magcol_trend;
      if(p->matchstringid)
	{
	  ptr1 = (void *) stringidin;
	  type1 = VARTOOLS_TYPE_STRING;
	}
      else
	{
	  ptr1 = (void *) &jdin;
	  type1 = VARTOOLS_TYPE_DOUBLE;
	}
      ptr2 = &magin;
      type2 = VARTOOLS_TYPE_DOUBLE;
    }
  else
    {
      col1 = tfa->magcol_trend;
      col2 = tfa->JDcol_trend;
      ptr1 = (void *) &magin;
      type1 = VARTOOLS_TYPE_DOUBLE;
      if(p->matchstringid)
	{
	  ptr2 = (void *) stringidin;
	  type2 = VARTOOLS_TYPE_STRING;
	}
      else
	{
	  ptr2 = (void *) &jdin;
	  type2 = VARTOOLS_TYPE_DOUBLE;
	}
    }

  if((dates = fopen(tfa->dates_name,"r")) == NULL)
    {
      error2(ERR_FILENOTFOUND,tfa->dates_name);
    }
  if((trend_list = fopen(tfa->trend_list_name,"r")) == NULL)
    {
      error2(ERR_FILENOTFOUND,tfa->trend_list_name);
    }
  tfa->Njd = 0;
  while(gnu_getline(&line,&line_size,dates) >= 0)
    {
      tfa->Njd++;
    }
  rewind(dates);

  if(p->matchstringid)
    {
      if((tfa->stringid = (char **) malloc(tfa->Njd * sizeof(char *))) == NULL ||
	 (tfa->stringid_idx = (int *) malloc(tfa->Njd * sizeof(int))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0;i<tfa->Njd;i++)
	{
	  if((tfa->stringid[i] = (char *) malloc(MAXIDSTRINGLENGTH * sizeof(char))) == NULL)
	    error(ERR_MEMALLOC);
	  tfa->stringid_idx[i] = i;
	}
      tfa->Njd = 0;
      while(gnu_getline(&line,&line_size,dates) >= 0)
	{
	  sscanf(line,"%s ",tfa->stringid[tfa->Njd]);
	  tfa->Njd++;
	}
      mysortstringint(tfa->Njd, MAXIDSTRINGLENGTH, tfa->stringid, tfa->stringid_idx);
    }
  else
    {
      if((tfa->JD = (double *) malloc(tfa->Njd * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);

      tfa->Njd = 0;
      while(gnu_getline(&line,&line_size,dates) >= 0)
	{
	  sscanf(line,"%s %lf",dums,&(tfa->JD[tfa->Njd]));
	  tfa->Njd++;
	}
      mysort1(tfa->Njd,tfa->JD);
    }

  fclose(dates);

  tfa->clipping = DEFAULT_CLIPPING_VALUE;

  if((tfa->Njd_mout = (int *) malloc(Nthreads * sizeof(int))) == NULL)
    error(ERR_MEMALLOC);

  for(i=0; i < Nthreads; i++)
    tfa->Njd_mout[i] = tfa->Njd;
  /*
  if((tfa->m_out = (double *) malloc(tfa->Njd_mout * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  */
  if((tfa->m_out = (double **) malloc(Nthreads * sizeof(double **))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0; i < Nthreads; i++) {
    if((tfa->m_out[i] = (double *) malloc(tfa->Njd_mout[i] * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
  }

  if((tfa->trends = (double **) malloc(tfa->Njd * sizeof(double *))) == NULL)
    error(ERR_MEMALLOC);

  tfa->Ntrends = 0;
  while(gnu_getline(&line,&line_size,trend_list) >= 0)
    tfa->Ntrends++;

  if((tfa->trendx = (double *) malloc(tfa->Ntrends * sizeof(double))) == NULL ||
     (tfa->trendy = (double *) malloc(tfa->Ntrends * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0;i<tfa->Njd;i++)
    {
      if((tfa->trends[i] = (double *) malloc(tfa->Ntrends * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
  if((tfa->trend_names = (char **) malloc(tfa->Ntrends * sizeof(char *))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0;i<tfa->Ntrends;i++)
    {
      if((tfa->trend_names[i] = (char *) malloc(MAXSTRINGLENGTH * sizeof(char))) == NULL)
	error(ERR_MEMALLOC);
    }
  rewind(trend_list);
  i=0;
  while(gnu_getline(&line,&line_size,trend_list) >= 0)
    {
      sscanf(line,"%s %lf %lf",trend_name, &tfa->trendx[i], &tfa->trendy[i]);
      sprintf(tfa->trend_names[i],"%s",trend_name);
#ifdef USECFITSIO
      jj = strlen(tfa->trend_names[i]);
      ii = jj - 5;
      if(ii > 0) {
	if(!strcmp(&(tfa->trend_names[i][ii]),".fits")) {
	  /* This is a fits light curve */
	  if(ReadFitsTFATemplate(p,tfa->jdcol_isfromheader,tfa->jdcol_headername,&(tfa->JDcol_trend),tfa->magcol_isfromheader,tfa->magcol_headername,&(tfa->magcol_trend),tfa->trend_names[i],&sizetmpvec,&lengthtmp,&magin_tmp,&jdin_tmp,&stringid_tmp,&stringid_tmpidx)) {
	    error2(ERR_TFATEMPLATEFITSREADERROR,trend_name);
	  }
	  if(!p->matchstringid) {
	    j = 0;
	    kk = 0;
	    jdin_last = 0.;
	    while(kk < lengthtmp)
	      {
		if(!j) {
		  while((j < tfa->Njd ? (jdin_tmp[kk] > tfa->JD[j] + JDTOL) : 0))
		    {
		      tfa->trends[j][i] = sqrt(-1);
		      j++;
		    }
		  if((j < tfa->Njd ? (jdin_tmp[kk] < tfa->JD[j] + JDTOL && jdin_tmp[kk] > tfa->JD[j] - JDTOL) : 0))
		    {
		      tfa->trends[j][i] = magin_tmp[kk];
		      j++;
		    }
		}
		else {
		  if(jdin_tmp[kk] > jdin_last) {
		    while((j < tfa->Njd ? (jdin_tmp[kk] > tfa->JD[j] + JDTOL) : 0))
		      {
			tfa->trends[j][i] = sqrt(-1);
			j++;
		      }
		    if((j < tfa->Njd ? (jdin_tmp[kk] < tfa->JD[j] + JDTOL && jdin_tmp[kk] > tfa->JD[j] - JDTOL) : 0))
		      {
			tfa->trends[j][i] = magin_tmp[kk];
			j++;
		      }
		  }
		  else {
		    j = 0;
		    while((j < tfa->Njd ? (jdin_tmp[kk] > tfa->JD[j] + JDTOL) : 0))
		      {
			tfa->trends[j][i] = sqrt(-1);
			j++;
		      }
		    if((j < tfa->Njd ? (jdin_tmp[kk] < tfa->JD[j] + JDTOL && jdin_tmp[kk] > tfa->JD[j] - JDTOL) : 0))
		      {
			tfa->trends[j][i] = magin_tmp[kk];
			j++;
		      }
		  }
		}
		jdin_last = jdin_tmp[kk];
		kk++;
	      }
	    for(;j < tfa->Njd; j++)
	      tfa->trends[j][i] = sqrt(-1);
	  } else {
	    j = 0;
	    kk = 0;
	    jdin_last = 0.;
	    mysortstringint(lengthtmp, MAXIDSTRINGLENGTH, stringid_tmp, stringid_tmpidx);
	    j=0;
	    k=0;
	    while(k < lengthtmp && j < tfa->Njd)
	      {
		while((j < tfa->Njd ? strncmp(stringid_tmp[stringid_tmpidx[k]],tfa->stringid[tfa->stringid_idx[j]],MAXIDSTRINGLENGTH) > 0 : 0))
		  {
		    tfa->trends[tfa->stringid_idx[j]][i] = sqrt(-1);
		    j++;
		  }
		if( j < tfa->Njd ? !strncmp(stringid_tmp[stringid_tmpidx[k]],tfa->stringid[tfa->stringid_idx[j]],MAXIDSTRINGLENGTH) : 0)
		  {
		    tfa->trends[tfa->stringid_idx[j]][i] = magin_tmp[stringid_tmpidx[k]];
		    j++;
		  }
		k++;
	      }
	    for(;j < tfa->Njd; j++)
	      tfa->trends[tfa->stringid_idx[j]][i] = sqrt(-1);
	  }
	  i++;
	  continue;
	}
      }
#endif
      if((trendin = fopen(trend_name,"r")) == NULL)
	error2(ERR_FILENOTFOUND,trend_name);
      if(!p->matchstringid)
	{
	  /* We are not using string-ids for the matching, we don't need to read the template light curve into a separate data vector */
	  j = 0;
	  lcline = 0;
	  while((lcline < tfa->Nskip_trend ? gnu_getline(&line,&line_size,trendin) >= 0 : 0))
	    lcline++;
	  jdin_last = 0.;
	  while(gnu_getline(&line,&line_size,trendin) >= 0)
	    {
	      if(line[0] != '#')
		{
		  lccol = 1;
		  lcindx = 0;
		  while (lccol < col1)
		    {
		      lcindx += skipone(&line[lcindx]);
		      lccol++;
		    }
		  if(line[lcindx] != '\0' && line[lcindx] != '\n')
		    {
		      lcindx_old = lcindx;
		      lcindx += parseone(&line[lcindx], ptr1, type1);
		      lccol++;
		    }
		  else
		    error2(ERR_INPUTMISSINGCOLUMN,trend_name);
		  if(lccol > col2)
		    {
		      lcindx = lcindx_old;
		      lcindx += parseone(&line[lcindx], ptr2, type2);
		      lccol++;
		    }
		  else
		    {
		      while (lccol < col2)
			{
			  lcindx += skipone(&line[lcindx]);
			  lccol++;
			}
		      if(line[lcindx] != '\0' && line[lcindx] != '\n')
			{
			  lcindx += parseone(&line[lcindx], ptr2, type2);
			  lccol++;
			}
		      else
			error2(ERR_INPUTMISSINGCOLUMN,trend_name);
		    }

		  if(!j) {
		    while((j < tfa->Njd ? (jdin > tfa->JD[j] + JDTOL) : 0))
		      {
			tfa->trends[j][i] = sqrt(-1);
			j++;
		      }
		    if((j < tfa->Njd ? (jdin < tfa->JD[j] + JDTOL && jdin > tfa->JD[j] - JDTOL) : 0))
		      {
			tfa->trends[j][i] = magin;
			j++;
		      }
		  }
		  else {
		    if(jdin > jdin_last) {
		      while((j < tfa->Njd ? (jdin > tfa->JD[j] + JDTOL) : 0))
			{
			  tfa->trends[j][i] = sqrt(-1);
			  j++;
			}
		      if((j < tfa->Njd ? (jdin < tfa->JD[j] + JDTOL && jdin > tfa->JD[j] - JDTOL) : 0))
			{
			  tfa->trends[j][i] = magin;
			  j++;
			}
		    }
		    else {
		      j = 0;
		      while((j < tfa->Njd ? (jdin > tfa->JD[j] + JDTOL) : 0))
			{
			  tfa->trends[j][i] = sqrt(-1);
			  j++;
			}
		      if((j < tfa->Njd ? (jdin < tfa->JD[j] + JDTOL && jdin > tfa->JD[j] - JDTOL) : 0))
			{
			  tfa->trends[j][i] = magin;
			  j++;
			}
		    }
		  }
		  jdin_last = jdin;
		}
	    }
	  for(;j < tfa->Njd; j++)
	    tfa->trends[j][i] = sqrt(-1);
	}
      else
	{
	  /* We are doing string-id matching */
	  j = 0;
	  lcline = 0;
	  while((lcline < tfa->Nskip_trend ? gnu_getline(&line,&line_size,trendin) >= 0 : 0))
	    lcline++;
	  /* Get the length of the trendin file */
	  lengthtmp = 0;
	  while(gnu_getline(&line,&line_size,trendin) >= 0)
	    {
	      if(line[0] != '#')
		lengthtmp++;
	    }
	  /* Add space to the temporary vectors if necessary */
	  if(lengthtmp > sizetmpvec)
	    {
	      if(!sizetmpvec)
		{
		  sizetmpvec = lengthtmp;
		  if((magin_tmp = (double *) malloc(sizetmpvec * sizeof(double))) == NULL ||
		     (stringid_tmp = (char **) malloc(sizetmpvec * sizeof(char *))) == NULL ||
		     (stringid_tmpidx = (int *) malloc(sizetmpvec * sizeof(int))) == NULL)
		    error(ERR_MEMALLOC);
		  for(j=0;j<sizetmpvec;j++)
		    {
		      if((stringid_tmp[j] = (char *) malloc(MAXIDSTRINGLENGTH * sizeof(char))) == NULL)
			error(ERR_MEMALLOC);
		    }
		}
	      else
		{
		  if((magin_tmp = (double *) realloc(magin_tmp, lengthtmp * sizeof(double))) == NULL ||
		     (stringid_tmp = (char **) realloc(stringid_tmp, lengthtmp * sizeof(char *))) == NULL ||
		     (stringid_tmpidx = (int *) realloc(stringid_tmpidx, lengthtmp * sizeof(int))) == NULL)
		    error(ERR_MEMALLOC);
		  for(j=sizetmpvec;j<lengthtmp;j++)
		    {
		      if((stringid_tmp[j] = (char *) malloc(MAXIDSTRINGLENGTH * sizeof(char))) == NULL)
			error(ERR_MEMALLOC);
		    }
		  sizetmpvec = lengthtmp;
		}
	    }
	  rewind(trendin);

	  lengthtmp = 0;
	  lcline = 0;
	  while((lcline < tfa->Nskip_trend ? gnu_getline(&line,&line_size,trendin) >= 0 : 0))
	    lcline++;
	  /* Read in the trend light curve */
	  while(gnu_getline(&line,&line_size,trendin) >= 0)
	    {
	      if(line[0] != '#')
		{
		  lccol = 1;
		  lcindx = 0;
		  while (lccol < col1)
		    {
		      lcindx += skipone(&line[lcindx]);
		      lccol++;
		    }
		  if(line[lcindx] != '\0' && line[lcindx] != '\n')
		    {
		      lcindx_old = lcindx;
		      lcindx += parseone(&line[lcindx], ptr1, type1);
		      lccol++;
		    }
		  else
		    error2(ERR_INPUTMISSINGCOLUMN,trend_name);
		  if(lccol > col2)
		    {
		      lcindx = lcindx_old;
		      lcindx += parseone(&line[lcindx], ptr2, type2);
		      lccol++;
		    }
		  else
		    {
		      while (lccol < col2)
			{
			  lcindx += skipone(&line[lcindx]);
			  lccol++;
			}
		      if(line[lcindx] != '\0' && line[lcindx] != '\n')
			{
			  lcindx += parseone(&line[lcindx], ptr2, type2);
			  lccol++;
			}
		      else
			error2(ERR_INPUTMISSINGCOLUMN,trend_name);
		    }

		  strncpy(stringid_tmp[lengthtmp],stringidin,MAXIDSTRINGLENGTH);
		  magin_tmp[lengthtmp] = magin;
		  stringid_tmpidx[lengthtmp] = lengthtmp;
		  lengthtmp++;
		}
	    }
	  mysortstringint(lengthtmp, MAXIDSTRINGLENGTH, stringid_tmp, stringid_tmpidx);
	  j=0;
	  k=0;
	  while(k < lengthtmp && j < tfa->Njd)
	    {
	      while((j < tfa->Njd ? strncmp(stringid_tmp[stringid_tmpidx[k]],tfa->stringid[tfa->stringid_idx[j]],MAXIDSTRINGLENGTH) > 0 : 0))
		{
		  tfa->trends[tfa->stringid_idx[j]][i] = sqrt(-1);
		  j++;
		}
	      if( j < tfa->Njd ? !strncmp(stringid_tmp[stringid_tmpidx[k]],tfa->stringid[tfa->stringid_idx[j]],MAXIDSTRINGLENGTH) : 0)
		{
		  tfa->trends[tfa->stringid_idx[j]][i] = magin_tmp[stringid_tmpidx[k]];
		  j++;
		}
	      k++;
	    }
	  for(;j < tfa->Njd; j++)
	    tfa->trends[tfa->stringid_idx[j]][i] = sqrt(-1);
	}
      fclose(trendin);
      i++;
    }
  if(jdin_tmp != NULL)
    free(jdin_tmp);
  if(magin_tmp != NULL)
    free(magin_tmp);
  if(stringid_tmpidx != NULL)
    free(stringid_tmpidx);
  if(stringid_tmp != NULL)
    {
      for(j=0;j<sizetmpvec;j++)
	free(stringid_tmp[j]);
      free(stringid_tmp);
    }
  fclose(trend_list);

  free(line);

  zero_trend_averages(tfa->Njd,tfa->Ntrends,tfa->trends,tfa->clipping);

  /* Compute the Singular Value Decomposition of the trends matrix */
  if((tfa->u = (double **) malloc(tfa->Njd * sizeof(double *))) == NULL ||
     (tfa->v = (double **) malloc(tfa->Ntrends * sizeof(double *))) == NULL ||
     (tfa->w1 = (double *) malloc(tfa->Ntrends * sizeof(double))) == NULL)
    {
      fprintf(stderr,"Memory Allocation Error in Function main\n");
      exit(3);
    }

  for(i=0;i<tfa->Njd;i++)
    if((tfa->u[i] = (double *) malloc(tfa->Ntrends * sizeof(double))) == NULL)
      {
	error(ERR_MEMALLOC);
      }
  for(i=0;i<tfa->Ntrends;i++)
    if((tfa->v[i] = (double *) malloc(tfa->Ntrends * sizeof(double))) == NULL)
      {
	error(ERR_MEMALLOC);
      }

  if((tfa->u2 = (double ***) malloc(Nthreads * sizeof(double **))) == NULL ||
     (tfa->v2 = (double ***) malloc(Nthreads * sizeof(double **))) == NULL ||
     (tfa->w2 = (double **) malloc(Nthreads * sizeof(double *))) == NULL ||
     (tfa->a = (double **) malloc(Nthreads * sizeof(double *))) == NULL ||
     (tfa->b = (double **) malloc(Nthreads * sizeof(double *))) == NULL)
    error(ERR_MEMALLOC);

  for(i=0; i < Nthreads; i++) {
    if((tfa->u2[i] = (double **) malloc(tfa->Njd * sizeof(double *))) == NULL ||
       (tfa->v2[i] = (double **) malloc(tfa->Ntrends * sizeof(double *))) == NULL ||
       (tfa->w2[i] = (double *) malloc(tfa->Ntrends * sizeof(double))) == NULL ||
       (tfa->a[i] = (double *) malloc(tfa->Ntrends * sizeof(double))) == NULL ||
       (tfa->b[i] = (double *) malloc(tfa->Njd * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);

    for(j=0; j < tfa->Njd; j++) {
      if((tfa->u2[i][j] = (double *) malloc(tfa->Ntrends * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
    for(j=0; j < tfa->Ntrends; j++) {
      if((tfa->v2[i][j] = (double *) malloc(tfa->Ntrends * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
  }

  /* Compute the singular value decomposition for all trends in the trend-list right now */
  for(i=0;i<tfa->Njd;i++)
    for(j=0;j<tfa->Ntrends;j++)
      tfa->u[i][j] = tfa->trends[i][j];
  svdcmp(tfa->u,tfa->Njd,tfa->Ntrends,tfa->w1,tfa->v);
}

void detrend_tfa(_TFA *tfa, int N, double *t, double *m, double *e, double lcx, double lcy, char *lc_name, char *coeff_file_name, int coeff_flag, int correctlc, int outlc, char *lc_out_name, double *ave_out, double *rms_out, int matchstringid, char **stringid, int *stringid_idx, int threadid)
{
  int isatrend_flag, k, i, l, h, n;
  FILE *coeff_file, *lcout;
  double val1, val2;

  /* Increase the output vector length if necessary */
  if(N > tfa->Njd_mout[threadid])
    {
      tfa->Njd_mout[threadid] = N;
      if((tfa->m_out[threadid] = (double *) realloc((void *) tfa->m_out[threadid],N * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }

  /* set-up the b vector for running tfa */
  if(!matchstringid)
    {
      for(i=0,k=0; i < N; i++)
	{
	  while((k < tfa->Njd ? (t[i] > tfa->JD[k] + JDTOL) : 0))
	    {
	      tfa->b[threadid][k] = sqrt(-1);
	      k++;
	    }
	  if ( k < tfa->Njd ? (t[i] < tfa->JD[k] + JDTOL && t[i] > tfa->JD[k] - JDTOL) : 0)
	    {
	      tfa->b[threadid][k] = m[i];
	      k++;
	    }
	}
      for(;k<tfa->Njd;k++)
	tfa->b[threadid][k] = sqrt(-1);
    }
  else
    {
      for(i=0,k=0; i < N; i++)
	{
	  while((k < tfa->Njd ? strncmp(stringid[stringid_idx[i]],tfa->stringid[tfa->stringid_idx[k]],MAXIDSTRINGLENGTH) > 0 : 0))
	    {
	      tfa->b[threadid][tfa->stringid_idx[k]] = sqrt(-1);
	      k++;
	    }
	  if ( k < tfa->Njd ? !strncmp(stringid[stringid_idx[i]],tfa->stringid[tfa->stringid_idx[k]],MAXIDSTRINGLENGTH) : 0)
	    {
	      tfa->b[threadid][tfa->stringid_idx[k]] = m[stringid_idx[i]];
	      k++;
	    }
	}
      for(;k<tfa->Njd;k++)
	tfa->b[threadid][tfa->stringid_idx[k]] = sqrt(-1);
    }
  zero_lc_averages(tfa->Njd,tfa->b[threadid],tfa->clipping);

  isatrend_flag = 0;
  for(i=0;i<tfa->Ntrends;i++)
    {
      if((!strncmp(lc_name,tfa->trend_names[i],strlen(tfa->trend_names[i])) && strlen(lc_name) == strlen(tfa->trend_names[i])) || sqrt((((lcx - tfa->trendx[i])*(lcx - tfa->trendx[i])) + ((lcy - tfa->trendy[i])*(lcy - tfa->trendy[i])))) < tfa->pixelsep)
	isatrend_flag = 1;
    }
  if(!isatrend_flag)
    svbksb(tfa->u,tfa->w1,tfa->v,tfa->Njd,tfa->Ntrends,tfa->b[threadid],tfa->a[threadid]);
  else
    {
      /* Remove all near-by stars from the trend-list */
      l = 0;
      for(i=0;i<tfa->Ntrends;i++)
	{
	  if(!(!strncmp(lc_name,tfa->trend_names[i],strlen(tfa->trend_names[i])) && strlen(lc_name) == strlen(tfa->trend_names[i])) && !(sqrt((((lcx - tfa->trendx[i])*(lcx - tfa->trendx[i])) + ((lcy - tfa->trendy[i])*(lcy - tfa->trendy[i])))) < tfa->pixelsep))
	    {
	      for(k=0;k<tfa->Njd;k++)
		{
		  tfa->u2[threadid][k][l] = tfa->trends[k][i];
		}
	      l++;
	    }
	}
      svdcmp(tfa->u2[threadid],tfa->Njd,l,tfa->w2[threadid],tfa->v2[threadid]);
      svbksb(tfa->u2[threadid],tfa->w2[threadid],tfa->v2[threadid],tfa->Njd,l,tfa->b[threadid],tfa->a[threadid]);
    }
  /* Write out the coefficients if we're doing that sort of thing */
  if(coeff_flag)
    {
      if((coeff_file = fopen(coeff_file_name,"w")) == NULL)
	{
	  error2(ERR_CANNOTWRITE,coeff_file_name);
	}
      l=0;
      for(i=0; i<tfa->Ntrends;i++)
	{
	  if(!(!strncmp(lc_name,tfa->trend_names[i],strlen(tfa->trend_names[i])) && strlen(lc_name) == strlen(tfa->trend_names[i])) && !(sqrt((((lcx - tfa->trendx[i])*(lcx - tfa->trendx[i])) + ((lcy - tfa->trendy[i])*(lcy - tfa->trendy[i])))) < tfa->pixelsep))
	    {
	      fprintf(coeff_file,"%s %f\n",tfa->trend_names[i],tfa->a[threadid][l]);
	      l++;
	    }
	}
      fclose(coeff_file);
    }
  /* Open the output light curve, correct the light curve and write it out */
  if(outlc)
    {
      if((lcout = fopen(lc_out_name,"w")) == NULL)
	error2(ERR_CANNOTWRITE,lc_out_name);
    }
  l=0;
  for(i=0;i<N;i++)
    tfa->m_out[threadid][i] = 0.;
  if(correctlc || outlc)
    {
      if(!matchstringid)
	{
	  for(i=0;i<tfa->Ntrends;i++)
	    {
	      if(!(!strncmp(lc_name,tfa->trend_names[i],strlen(tfa->trend_names[i])) && strlen(lc_name) == strlen(tfa->trend_names[i])) && !(sqrt((((lcx - tfa->trendx[i])*(lcx - tfa->trendx[i])) + ((lcy - tfa->trendy[i])*(lcy - tfa->trendy[i])))) < tfa->pixelsep))
		{
		  k = 0;
		  for(h = 0; h<N; h++)
		    {
		      while((k < tfa->Njd ? (t[h] > tfa->JD[k] + JDTOL) : 0)) k++;
		      if(k < tfa->Njd ? (t[h] < tfa->JD[k] + JDTOL && t[h] > tfa->JD[k] - JDTOL) : 0)
			{
			  tfa->m_out[threadid][h] += tfa->a[threadid][l]*tfa->trends[k][i];
			  k++;
			}
		    }
		  l++;
		}
	    }
	}
      else
	{
	  for(i=0;i<tfa->Ntrends;i++)
	    {
	      if(!(!strncmp(lc_name,tfa->trend_names[i],strlen(tfa->trend_names[i])) && strlen(lc_name) == strlen(tfa->trend_names[i])) && !(sqrt((((lcx - tfa->trendx[i])*(lcx - tfa->trendx[i])) + ((lcy - tfa->trendy[i])*(lcy - tfa->trendy[i])))) < tfa->pixelsep))
		{
		  k = 0;
		  for(h = 0; h<N; h++)
		    {
		      while((k < tfa->Njd ? strncmp(stringid[stringid_idx[h]],tfa->stringid[tfa->stringid_idx[k]],MAXIDSTRINGLENGTH) > 0 : 0)) k++;
		      if(k < tfa->Njd ? !strncmp(stringid[stringid_idx[h]],tfa->stringid[tfa->stringid_idx[k]],MAXIDSTRINGLENGTH) : 0)
			{
			  tfa->m_out[threadid][stringid_idx[h]] += tfa->a[threadid][l]*tfa->trends[tfa->stringid_idx[k]][i];
			  k++;
			}
		    }
		  l++;
		}
	    }
	}
      for(i=0;i<N;i++)
	{
	  if(outlc)
	    {
	      fprintf(lcout,"%f %f %f %f\n",t[i],m[i],tfa->m_out[threadid][i],e[i]);
	    }
	  if(correctlc)
	    {
	      m[i] = m[i] - tfa->m_out[threadid][i];
	    }
	}
      if(outlc)
	fclose(lcout);
    }
  if(correctlc)
    {
      n = 0;
      val1 = 0.;
      val2 = 0.;
      for(i=0;i<N;i++)
	{
	  n++;
	  val1 += m[i]*m[i];
	  val2 += m[i];
	}
    }
  else
    {
      n = 0;
      val1 = 0.;
      val2 = 0.;
      for(i = 0; i<N;i++)
	{
	  n++;
	  val1 += (m[i] - tfa->m_out[threadid][i])*(m[i] - tfa->m_out[threadid][i]);
	  val2 += (m[i] - tfa->m_out[threadid][i]);
	}
    }
  *ave_out = (double) (val2 / (double) n);
  *rms_out = sqrt((double) ((val1 / (double) n) - (val2 * val2 / (double) (n*n))));
}

void initialize_tfa_sr(_TFA_SR *tfa, int Nlcs, ProgramData *p)
{
  FILE *dates, *trendin, *trend_list, *listfile;
  double *magin_tmp = NULL;
  double *jdin_tmp = NULL;
  char **stringid_tmp = NULL;
  int *stringid_tmpidx = NULL;
  int lengthtmp, sizetmpvec = 0;
  char stringidin[MAXIDSTRINGLENGTH];
  double jdin, magin, jdin_last;
  void *ptr1, *ptr2;
  int i, j, k, i_, nbins, lcline, lccol, lcindx, lcindx_old, col1, col2, type1, type2, ii, jj, kk;
  void error2(int, char *);
  char trend_name[MAXSTRINGLENGTH], dums[MAXSTRINGLENGTH];
  char *line;
  size_t line_size = MAXSTRINGLENGTH;
  int use_harm, Nharm, Nsubharm;

  int Nthreads;

  line = malloc(line_size);

#ifdef PARALLEL
  Nthreads = p->Nproc_allow;
#else
  Nthreads = 1;
#endif

  use_harm = tfa->use_harm;
  Nharm = tfa->Nharm;
  Nsubharm = tfa->Nsubharm;

  if(tfa->JDcol_trend <= tfa->magcol_trend)
    {
      col1 = tfa->JDcol_trend;
      col2 = tfa->magcol_trend;
      if(p->matchstringid)
	{
	  ptr1 = (void *) stringidin;
	  type1 = VARTOOLS_TYPE_STRING;
	}
      else
	{
	  ptr1 = (void *) &jdin;
	  type1 = VARTOOLS_TYPE_DOUBLE;
	}
      ptr2 = &magin;
      type2 = VARTOOLS_TYPE_DOUBLE;
    }
  else
    {
      col1 = tfa->magcol_trend;
      col2 = tfa->JDcol_trend;
      ptr1 = (void *) &magin;
      type1 = VARTOOLS_TYPE_DOUBLE;
      if(p->matchstringid)
	{
	  ptr2 = (void *) stringidin;
	  type2 = VARTOOLS_TYPE_STRING;
	}
      else
	{
	  ptr2 = (void *) &jdin;
	  type2 = VARTOOLS_TYPE_DOUBLE;
	}
    }

  nbins = tfa->nbins;

  if(tfa->use_bin)
    {
      if((tfa->signal_bin = (double **) malloc(Nthreads * sizeof(double *))) == NULL ||
	 (tfa->signal_bin_N = (int **) malloc(Nthreads * sizeof(int *))) == NULL)
	error(ERR_MEMALLOC);

      for(i_=0; i_ < Nthreads; i_++) {
	if((tfa->signal_bin[i_] = (double *) malloc(nbins * sizeof(double))) == NULL ||
	   (tfa->signal_bin_N[i_] = (int *) malloc(nbins * sizeof(int))) == NULL)
	  error(ERR_MEMALLOC);
      }
    }

  /* Read in the signal list if we're doing that */
  if(!tfa->use_bin && !tfa->use_harm)
    {
      if((listfile = fopen(tfa->signal_listname,"r")) == NULL)
	error2(ERR_FILENOTFOUND,tfa->signal_listname);
      i=0;
      while(gnu_getline(&line,&line_size,listfile) >= 0)
	i++;
      if(i != Nlcs)
	error2(ERR_GETLSAMPTHRESH_FILETOSHORT, tfa->signal_listname);
      rewind(listfile);
      if((tfa->signalfilenames = (char **) malloc(Nlcs * sizeof(char *))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0;i<Nlcs;i++)
	{
	  if((tfa->signalfilenames[i] = (char *) malloc(MAXLEN * sizeof(char))) == NULL)
	    error(ERR_MEMALLOC);
	}
      i = 0;
      while(gnu_getline(&line,&line_size,listfile) >= 0)
	{
	  sscanf(line,"%s",tfa->signalfilenames[i]);
	  i++;
	}
      fclose(listfile);
    }

  if((dates = fopen(tfa->dates_name,"r")) == NULL)
    {
      error2(ERR_FILENOTFOUND,tfa->dates_name);
    }
  if((trend_list = fopen(tfa->trend_list_name,"r")) == NULL)
    {
      error2(ERR_FILENOTFOUND,tfa->trend_list_name);
    }
  tfa->Njd = 0;
  while(gnu_getline(&line,&line_size,dates) >= 0)
    {
      tfa->Njd++;
    }
  rewind(dates);

  if(p->matchstringid)
    {
      if((tfa->stringid = (char **) malloc(tfa->Njd * sizeof(char *))) == NULL ||
	 (tfa->stringid_idx = (int *) malloc(tfa->Njd * sizeof(int))) == NULL)
	error(ERR_MEMALLOC);
      if(tfa->use_harm)
	{
	  if((tfa->JD = (double *) malloc(tfa->Njd * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      for(i=0;i<tfa->Njd;i++)
	{
	  if((tfa->stringid[i] = (char *) malloc(MAXIDSTRINGLENGTH * sizeof(char))) == NULL)
	    error(ERR_MEMALLOC);
	  tfa->stringid_idx[i] = i;
	}
      tfa->Njd = 0;
      if(!tfa->use_harm)
	{
	  while(gnu_getline(&line,&line_size,dates) >= 0)
	    {
	      sscanf(line,"%s ",tfa->stringid[tfa->Njd]);
	      tfa->Njd++;
	    }
	}
      else
	{
	  while(gnu_getline(&line,&line_size,dates) >= 0)
	    {
	      sscanf(line,"%s %lf",tfa->stringid[tfa->Njd],&(tfa->JD[tfa->Njd]));
	      tfa->Njd++;
	    }
	}
      mysortstringint(tfa->Njd, MAXIDSTRINGLENGTH, tfa->stringid, tfa->stringid_idx);
    }
  else
    {
      if((tfa->JD = (double *) malloc(tfa->Njd * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);

      tfa->Njd = 0;
      while(gnu_getline(&line,&line_size,dates) >= 0)
	{
	  sscanf(line,"%s %lf",dums,&(tfa->JD[tfa->Njd]));
	  tfa->Njd++;
	}
      mysort1(tfa->Njd, tfa->JD);
    }

  fclose(dates);

  tfa->clipping = DEFAULT_CLIPPING_VALUE;

  if((tfa->Njd_mout = (int *) malloc(Nthreads * sizeof(int))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0; i < Nthreads; i++)
    tfa->Njd_mout[i] = tfa->Njd;

  if((tfa->m_out = (double **) malloc(Nthreads * sizeof(double *))) == NULL ||
     (tfa->signal = (double **) malloc(Nthreads * sizeof(double *))) == NULL ||
     (tfa->signal_bin_ids = (int **) malloc(Nthreads * sizeof(int *))) == NULL ||
     (tfa->inputsignal = (double **) malloc(Nthreads * sizeof(double *))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0; i < Nthreads; i++) {
    if((tfa->m_out[i] = (double *) malloc(tfa->Njd_mout[i] * sizeof(double))) == NULL ||
       (tfa->signal[i] = (double *) malloc(tfa->Njd_mout[i] * sizeof(double))) == NULL ||
       (tfa->signal_bin_ids[i] = (int *) malloc(tfa->Njd_mout[i] * sizeof(int))) == NULL||
       (tfa->inputsignal[i] = (double *) malloc(tfa->Njd_mout[i] * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
  }
  if(tfa->use_harm)
    {
      if((tfa->harmterm = (double ***) malloc(Nthreads * sizeof(double **))) == NULL)
	error(ERR_MEMALLOC);
      for(i_ = 0; i_ < Nthreads; i_++) {
	if((tfa->harmterm[i_] = (double **) malloc(tfa->Njd_mout[i_] * sizeof(double *))) == NULL)
	  error(ERR_MEMALLOC);
	for(i=0;i<tfa->Njd_mout[i_];i++)
	  {
	    if((tfa->harmterm[i_][i] = (double *) malloc((2*(1 + tfa->Nharm + tfa->Nsubharm))*sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	  }
      }
    }

  if((tfa->trends = (double **) malloc(tfa->Njd * sizeof(double *))) == NULL)
    error(ERR_MEMALLOC);

  tfa->Ntrends = 0;
  while(gnu_getline(&line,&line_size,trend_list) >= 0)
    tfa->Ntrends++;

  if((tfa->trendx = (double *) malloc(tfa->Ntrends * sizeof(double))) == NULL ||
     (tfa->trendy = (double *) malloc(tfa->Ntrends * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0;i<tfa->Njd;i++)
    {
      if((tfa->trends[i] = (double *) malloc(tfa->Ntrends * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
  if((tfa->trend_names = (char **) malloc(tfa->Ntrends * sizeof(char *))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0;i<tfa->Ntrends;i++)
    {
      if((tfa->trend_names[i] = (char *) malloc(MAXSTRINGLENGTH * sizeof(char))) == NULL)
	error(ERR_MEMALLOC);
    }
  rewind(trend_list);
  i=0;
  while(gnu_getline(&line,&line_size,trend_list) >= 0)
    {
      sscanf(line,"%s %lf %lf",trend_name, &tfa->trendx[i], &tfa->trendy[i]);
      sprintf(tfa->trend_names[i],"%s",trend_name);
#ifdef USECFITSIO
      jj = strlen(tfa->trend_names[i]);
      ii = jj - 5;
      if(ii > 0) {
	if(!strcmp(&(tfa->trend_names[i][ii]),".fits")) {
	  /* This is a fits light curve */
	  if(ReadFitsTFATemplate(p,tfa->jdcol_isfromheader,tfa->jdcol_headername,&(tfa->JDcol_trend),tfa->magcol_isfromheader,tfa->magcol_headername,&(tfa->magcol_trend),tfa->trend_names[i],&sizetmpvec,&lengthtmp,&magin_tmp,&jdin_tmp,&stringid_tmp,&stringid_tmpidx)) {
	    error2(ERR_TFATEMPLATEFITSREADERROR,trend_name);
	  }
	  if(!p->matchstringid) {
	    j = 0;
	    kk = 0;
	    jdin_last = 0.;
	    while(kk < lengthtmp)
	      {
		if(!j) {
		  while((j < tfa->Njd ? (jdin_tmp[kk] > tfa->JD[j] + JDTOL) : 0))
		    {
		      tfa->trends[j][i] = sqrt(-1);
		      j++;
		    }
		  if((j < tfa->Njd ? (jdin_tmp[kk] < tfa->JD[j] + JDTOL && jdin_tmp[kk] > tfa->JD[j] - JDTOL) : 0))
		    {
		      tfa->trends[j][i] = magin_tmp[kk];
		      j++;
		    }
		}
		else {
		  if(jdin_tmp[kk] > jdin_last) {
		    while((j < tfa->Njd ? (jdin_tmp[kk] > tfa->JD[j] + JDTOL) : 0))
		      {
			tfa->trends[j][i] = sqrt(-1);
			j++;
		      }
		    if((j < tfa->Njd ? (jdin_tmp[kk] < tfa->JD[j] + JDTOL && jdin_tmp[kk] > tfa->JD[j] - JDTOL) : 0))
		      {
			tfa->trends[j][i] = magin_tmp[kk];
			j++;
		      }
		  }
		  else {
		    j = 0;
		    while((j < tfa->Njd ? (jdin_tmp[kk] > tfa->JD[j] + JDTOL) : 0))
		      {
			tfa->trends[j][i] = sqrt(-1);
			j++;
		      }
		    if((j < tfa->Njd ? (jdin_tmp[kk] < tfa->JD[j] + JDTOL && jdin_tmp[kk] > tfa->JD[j] - JDTOL) : 0))
		      {
			tfa->trends[j][i] = magin_tmp[kk];
			j++;
		      }
		  }
		}
		jdin_last = jdin_tmp[kk];
		kk++;
	      }
	    for(;j < tfa->Njd; j++)
	      tfa->trends[j][i] = sqrt(-1);
	  } else {
	    j = 0;
	    kk = 0;
	    jdin_last = 0.;
	    mysortstringint(lengthtmp, MAXIDSTRINGLENGTH, stringid_tmp, stringid_tmpidx);
	    j=0;
	    k=0;
	    while(k < lengthtmp && j < tfa->Njd)
	      {
		while((j < tfa->Njd ? strncmp(stringid_tmp[stringid_tmpidx[k]],tfa->stringid[tfa->stringid_idx[j]],MAXIDSTRINGLENGTH) > 0 : 0))
		  {
		    tfa->trends[tfa->stringid_idx[j]][i] = sqrt(-1);
		    j++;
		  }
		if( j < tfa->Njd ? !strncmp(stringid_tmp[stringid_tmpidx[k]],tfa->stringid[tfa->stringid_idx[j]],MAXIDSTRINGLENGTH) : 0)
		  {
		    tfa->trends[tfa->stringid_idx[j]][i] = magin_tmp[stringid_tmpidx[k]];
		    j++;
		  }
		k++;
	      }
	    for(;j < tfa->Njd; j++)
	      tfa->trends[tfa->stringid_idx[j]][i] = sqrt(-1);
	  }
	  i++;
	  continue;
	}
      }
#endif
      if((trendin = fopen(trend_name,"r")) == NULL)
	error2(ERR_FILENOTFOUND,trend_name);
      if(!p->matchstringid)
	{
	  /* We are not using string-ids for the matching, we don't need to read the template light curve into a separate data vector */
	  j = 0;
	  lcline = 0;
	  while((lcline < tfa->Nskip_trend ? gnu_getline(&line,&line_size,trendin) >= 0 : 0))
	    lcline++;
	  while(gnu_getline(&line,&line_size,trendin) >= 0)
	    {
	      if(line[0] != '#')
		{
		  lccol = 1;
		  lcindx = 0;
		  while (lccol < col1)
		    {
		      lcindx += skipone(&line[lcindx]);
		      lccol++;
		    }
		  if(line[lcindx] != '\0' && line[lcindx] != '\n')
		    {
		      lcindx_old = lcindx;
		      lcindx += parseone(&line[lcindx], ptr1, type1);
		      lccol++;
		    }
		  else
		    error2(ERR_INPUTMISSINGCOLUMN,trend_name);
		  if(lccol > col2)
		    {
		      lcindx = lcindx_old;
		      lcindx += parseone(&line[lcindx], ptr2, type2);
		      lccol++;
		    }
		  else
		    {
		      while (lccol < col2)
			{
			  lcindx += skipone(&line[lcindx]);
			  lccol++;
			}
		      if(line[lcindx] != '\0' && line[lcindx] != '\n')
			{
			  lcindx += parseone(&line[lcindx], ptr2, type2);
			  lccol++;
			}
		      else
			error2(ERR_INPUTMISSINGCOLUMN,trend_name);
		    }

		  while((j < tfa->Njd ? (jdin > tfa->JD[j] + JDTOL) : 0))
		    {
		      tfa->trends[j][i] = sqrt(-1);
		      j++;
		    }
		  if((j < tfa->Njd ? (jdin < tfa->JD[j] + JDTOL && jdin > tfa->JD[j] - JDTOL) : 0))
		    {
		      tfa->trends[j][i] = magin;
		      j++;
		    }
		}
	    }
	  for(;j < tfa->Njd; j++)
	    tfa->trends[j][i] = sqrt(-1);
	}
      else
	{
	  /* We are doing string-id matching */
	  j = 0;
	  lcline = 0;
	  while((lcline < tfa->Nskip_trend ? gnu_getline(&line,&line_size,trendin) >= 0 : 0))
	    lcline++;
	  /* Get the length of the trendin file */
	  lengthtmp = 0;
	  while(gnu_getline(&line,&line_size,trendin) >= 0)
	    {
	      if(line[0] != '#')
		lengthtmp++;
	    }
	  /* Add space to the temporary vectors if necessary */
	  if(lengthtmp > sizetmpvec)
	    {
	      if(!sizetmpvec)
		{
		  sizetmpvec = lengthtmp;
		  if((magin_tmp = (double *) malloc(sizetmpvec * sizeof(double))) == NULL ||
		     (stringid_tmp = (char **) malloc(sizetmpvec * sizeof(char *))) == NULL ||
		     (stringid_tmpidx = (int *) malloc(sizetmpvec * sizeof(int))) == NULL)
		    error(ERR_MEMALLOC);
		  for(j=0;j<sizetmpvec;j++)
		    {
		      if((stringid_tmp[j] = (char *) malloc(MAXIDSTRINGLENGTH * sizeof(char))) == NULL)
			error(ERR_MEMALLOC);
		    }
		}
	      else
		{
		  if((magin_tmp = (double *) realloc(magin_tmp, lengthtmp * sizeof(double))) == NULL ||
		     (stringid_tmp = (char **) realloc(stringid_tmp, lengthtmp * sizeof(char *))) == NULL ||
		     (stringid_tmpidx = (int *) realloc(stringid_tmpidx, lengthtmp * sizeof(int))) == NULL)
		    error(ERR_MEMALLOC);
		  for(j=sizetmpvec;j<lengthtmp;j++)
		    {
		      if((stringid_tmp[j] = (char *) malloc(MAXIDSTRINGLENGTH * sizeof(char))) == NULL)
			error(ERR_MEMALLOC);
		    }
		  sizetmpvec = lengthtmp;
		}
	    }
	  rewind(trendin);

	  lengthtmp = 0;
	  lcline = 0;
	  while((lcline < tfa->Nskip_trend ? gnu_getline(&line,&line_size,trendin) >= 0 : 0))
	    lcline++;
	  /* Read in the trend light curve */
	  while(gnu_getline(&line,&line_size,trendin) >= 0)
	    {
	      if(line[0] != '#')
		{
		  lccol = 1;
		  lcindx = 0;
		  while (lccol < col1)
		    {
		      lcindx += skipone(&line[lcindx]);
		      lccol++;
		    }
		  if(line[lcindx] != '\0' && line[lcindx] != '\n')
		    {
		      lcindx_old = lcindx;
		      lcindx += parseone(&line[lcindx], ptr1, type1);
		      lccol++;
		    }
		  else
		    error2(ERR_INPUTMISSINGCOLUMN,trend_name);
		  if(lccol > col2)
		    {
		      lcindx = lcindx_old;
		      lcindx += parseone(&line[lcindx], ptr2, type2);
		      lccol++;
		    }
		  else
		    {
		      while (lccol < col2)
			{
			  lcindx += skipone(&line[lcindx]);
			  lccol++;
			}
		      if(line[lcindx] != '\0' && line[lcindx] != '\n')
			{
			  lcindx += parseone(&line[lcindx], ptr2, type2);
			  lccol++;
			}
		      else
			error2(ERR_INPUTMISSINGCOLUMN,trend_name);
		    }

		  strncpy(stringid_tmp[lengthtmp],stringidin,MAXIDSTRINGLENGTH);
		  magin_tmp[lengthtmp] = magin;
		  stringid_tmpidx[lengthtmp] = lengthtmp;
		  lengthtmp++;
		}
	    }
	  mysortstringint(lengthtmp, MAXIDSTRINGLENGTH, stringid_tmp, stringid_tmpidx);
	  j=0;
	  k=0;
	  while(k < lengthtmp && j < tfa->Njd)
	    {
	      while((j < tfa->Njd ? strncmp(stringid_tmp[stringid_tmpidx[k]],tfa->stringid[tfa->stringid_idx[j]],MAXIDSTRINGLENGTH) > 0 : 0))
		{
		  tfa->trends[tfa->stringid_idx[j]][i] = sqrt(-1);
		  j++;
		}
	      if( j < tfa->Njd ? !strncmp(stringid_tmp[stringid_tmpidx[k]],tfa->stringid[tfa->stringid_idx[j]],MAXIDSTRINGLENGTH) : 0)
		{
		  tfa->trends[tfa->stringid_idx[j]][i] = magin_tmp[stringid_tmpidx[k]];
		  j++;
		}
	      k++;
	    }
	  for(;j < tfa->Njd; j++)
	    tfa->trends[tfa->stringid_idx[j]][i] = sqrt(-1);
	}
      fclose(trendin);
      i++;
    }
  if(magin_tmp != NULL)
    free(magin_tmp);
  if(jdin_tmp != NULL)
    free(jdin_tmp);
  if(stringid_tmpidx != NULL)
    free(stringid_tmpidx);
  if(stringid_tmp != NULL)
    {
      for(j=0;j<sizetmpvec;j++)
	free(stringid_tmp[j]);
      free(stringid_tmp);
    }
  fclose(trend_list);
  zero_trend_averages(tfa->Njd,tfa->Ntrends,tfa->trends,tfa->clipping);

  /* Make adjustments for the type of external parameter decorrelation */
  tfa->Ntfatot = tfa->Ntrends;
  if(tfa->decorrflag)
    {
      tfa->Ndecorr = 0;
      for(i=0;i<tfa->decorr_Nlcterms;i++)
	{
	  tfa->Ndecorr += tfa->decorr_lc_order[i];
	}
      if(0)
      //      if(!tfa->decorr_iterate)
	{
	  tfa->Ntfatot += tfa->Ndecorr;
	}
      else
	{
	  if(!tfa->decorr_iterate)
	    tfa->Ntfatot += tfa->Ndecorr;
	  if((tfa->decorr_trends = (double ***) malloc(Nthreads * sizeof(double **))) == NULL ||
	     (tfa->u_decorr = (double ***) malloc(Nthreads * sizeof(double **))) == NULL ||
	     (tfa->v_decorr = (double ***) malloc(Nthreads * sizeof(double **))) == NULL ||
	     (tfa->w1_decorr = (double **) malloc(Nthreads * sizeof(double *))) == NULL ||
	     (tfa->a_decorr = (double **) malloc(Nthreads * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);

	  for(i_ = 0; i_ < Nthreads; i_++) {

	    if((tfa->decorr_trends[i_] = (double **) malloc(tfa->Njd * sizeof(double *))) == NULL ||
	       (tfa->u_decorr[i_] = (double **) malloc(tfa->Njd * sizeof(double *))) == NULL ||
	       (tfa->v_decorr[i_] = (double **) malloc(tfa->Ndecorr * sizeof(double *))) == NULL ||
	       (tfa->w1_decorr[i_] = (double *) malloc(tfa->Ndecorr * sizeof(double))) == NULL ||
	       (tfa->a_decorr[i_] = (double *) malloc(tfa->Ndecorr * sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);

	    for(i=0;i<tfa->Njd;i++)
	      if((tfa->decorr_trends[i_][i] = (double *) malloc(tfa->Ndecorr * sizeof(double))) == NULL ||
		 (tfa->u_decorr[i_][i] = (double *) malloc(tfa->Ndecorr * sizeof(double))) == NULL)
		{
		  error(ERR_MEMALLOC);
		}
	    for(i=0;i<tfa->Ndecorr;i++)
	      if((tfa->v_decorr[i_][i] = (double *) malloc(tfa->Ndecorr * sizeof(double))) == NULL)
		{
		  error(ERR_MEMALLOC);
		}

	  }
	}
    }
  if(use_harm)
    tfa->Ntfatot += 2*(1 + Nharm + Nsubharm);

  free(line);

  /* Compute the Singular Value Decomposition of the trends matrix */
  if((tfa->u = (double **) malloc(tfa->Njd * sizeof(double *))) == NULL ||
     (tfa->v = (double **) malloc(tfa->Ntfatot * sizeof(double *))) == NULL ||
     (tfa->w1 = (double *) malloc(tfa->Ntfatot * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  for(i=0;i<tfa->Njd;i++)
    if((tfa->u[i] = (double *) malloc(tfa->Ntfatot * sizeof(double))) == NULL)
      {
	error(ERR_MEMALLOC);
      }
  for(i=0;i<tfa->Ntfatot;i++)
    if((tfa->v[i] = (double *) malloc(tfa->Ntfatot * sizeof(double))) == NULL)
      {
	error(ERR_MEMALLOC);
      }

  if((tfa->u2 = (double ***) malloc(Nthreads * sizeof(double **))) == NULL ||
     (tfa->v2 = (double ***) malloc(Nthreads * sizeof(double **))) == NULL ||
     (tfa->w2 = (double **) malloc(Nthreads * sizeof(double *))) == NULL ||
     (tfa->a = (double **) malloc(Nthreads * sizeof(double *))) == NULL ||
     (tfa->b = (double **) malloc(Nthreads * sizeof(double *))) == NULL ||
     (tfa->bstore = (double **) malloc(Nthreads * sizeof(double *))) == NULL)
    error(ERR_MEMALLOC);

  for(i_=0;i_<Nthreads;i_++) {
    if((tfa->u2[i_] = (double **) malloc(tfa->Njd * sizeof(double *))) == NULL ||
       (tfa->v2[i_] = (double **) malloc(tfa->Ntfatot * sizeof(double *))) == NULL ||
       (tfa->w2[i_] = (double *) malloc(tfa->Ntfatot * sizeof(double))) == NULL ||
       (tfa->a[i_] = (double *) malloc(tfa->Ntfatot * sizeof(double))) == NULL ||
       (tfa->b[i_] = (double *) malloc(tfa->Njd * sizeof(double))) == NULL ||
       (tfa->bstore[i_] = (double *) malloc(tfa->Njd * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);

    for(i=0;i<tfa->Njd;i++)
      if((tfa->u2[i_][i] = (double *) malloc(tfa->Ntfatot * sizeof(double))) == NULL)
	{
	  error(ERR_MEMALLOC);
	}
    for(i=0;i<tfa->Ntfatot;i++)
      if((tfa->v2[i_][i] = (double *) malloc(tfa->Ntfatot * sizeof(double))) == NULL)
	{
	  error(ERR_MEMALLOC);
	}
  }


  /* Compute the singular value decomposition for all trends in the trend-list right now */
  if(tfa->Ntfatot == tfa->Ntrends)
    {
      for(i=0;i<tfa->Njd;i++)
	for(j=0;j<tfa->Ntrends;j++)
	  tfa->u[i][j] = tfa->trends[i][j];
      svdcmp(tfa->u,tfa->Njd,tfa->Ntrends,tfa->w1,tfa->v);
    }
}

/* This function runs tfa iteratively in signal reconstruction mode */
void detrend_tfa_sr(_TFA_SR *tfa, int N, double *t, double *m, double *e, double lcx, double lcy, char *lc_name, char *coeff_file_name, int coeff_flag, int correctlc, int outlc, char *lc_out_name, double *ave_out, double *rms_out, double period, char *signalfilename, int matchstringid, char **stringid, int *stringid_idx, int lcindex, int threadid)
{
  int isatrend_flag, j, k, i, l, lstart, h, n, Ntrends_, iter;
  FILE *coeff_file, *lcout, *signalfile;
  double val1, val2, val3, val4, a, b, ph, **u_, **v_, *w_, rmsnew, rmsold;
  int dotfafirst, use_bin, nbins, use_period, maxiter, use_harm, Nharm, Nsubharm;
  double tfathresh, term;
  char *line;
  size_t line_size = MAXLEN;

  line = malloc(line_size);

  dotfafirst = tfa->dotfafirst;
  use_bin = tfa->use_bin;
  nbins = tfa->nbins;
  use_period = tfa->use_period;
  maxiter = tfa->maxiter;
  tfathresh = tfa->tfathresh;
  use_harm = tfa->use_harm;
  Nharm = tfa->Nharm;
  Nsubharm = tfa->Nsubharm;

  /* Increase the output vector length if necessary */
  if(N > tfa->Njd_mout[threadid])
    {
      tfa->Njd_mout[threadid] = N;
      if((tfa->m_out[threadid] = (double *) realloc((void *) tfa->m_out[threadid],N * sizeof(double))) == NULL ||
	 (tfa->signal[threadid] = (double *) realloc((void *) tfa->signal[threadid], N * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      if(use_harm)
	{
	  if((tfa->harmterm[threadid] = (double **) realloc(tfa->harmterm[threadid], N * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(i=0;i<N;i++)
	    {
	      if((tfa->harmterm[threadid][i] = (double *) malloc((2*(1 + Nharm + Nsubharm))*sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	}
      if(use_bin)
	{
	  if((tfa->signal_bin_ids[threadid] = (int *) realloc((void *) tfa->signal_bin_ids[threadid], N * sizeof(int))) == NULL ||
	     (tfa->inputsignal[threadid] = (double *) realloc((void *) tfa->inputsignal[threadid], N * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
    }

  if(!use_bin && !use_harm)
    {
      if((signalfile = fopen(signalfilename,"r")) == NULL)
	error2(ERR_FILENOTFOUND,signalfilename);
      i=0;
      while(gnu_getline(&line,&line_size,signalfile) >= 0)
	i++;
      if(i != N)
	error2(ERR_SIGFILEWRONGLENGTH,signalfilename);
      rewind(signalfile);
      i=0;
      while(gnu_getline(&line,&line_size,signalfile) >= 0)
	{
	  sscanf(line,"%lf %lf",&val1,&tfa->inputsignal[threadid][i]);
	  i++;
	}
      fclose(signalfile);
    }

  free(line);

  /* set-up the storage b vector for running tfa */
  val1 = 0.;
  val2 = 0.;
  l = 0;
  if(!matchstringid)
    {
      for(i=0,k=0; i < N; i++)
	{
	  while((k < tfa->Njd ? (t[i] > tfa->JD[k] + JDTOL) : 0))
	    {
	      tfa->bstore[threadid][k] = sqrt(-1);
	      k++;
	    }
	  if ( k < tfa->Njd ? (t[i] < tfa->JD[k] + JDTOL && t[i] > tfa->JD[k] - JDTOL) : 0)
	    {
	      tfa->bstore[threadid][k] = m[i];
	      val1 += m[i];
	      val2 += m[i]*m[i];
	      l++;
	      k++;
	    }
	}
      for(;k<tfa->Njd;k++)
	tfa->bstore[threadid][k] = sqrt(-1);
    }
  else
    {
      for(i=0,k=0; i < N; i++)
	{
	  while((k < tfa->Njd ? (strncmp(stringid[stringid_idx[i]], tfa->stringid[tfa->stringid_idx[k]],MAXIDSTRINGLENGTH) > 0) : 0))
	    {
	      tfa->bstore[threadid][tfa->stringid_idx[k]] = sqrt(-1);
	      k++;
	    }
	  if ( k < tfa->Njd ? !strncmp(stringid[stringid_idx[i]], tfa->stringid[tfa->stringid_idx[k]],MAXIDSTRINGLENGTH) : 0)
	    {
	      tfa->bstore[threadid][tfa->stringid_idx[k]] = m[stringid_idx[i]];
	      val1 += m[stringid_idx[i]];
	      val2 += m[stringid_idx[i]]*m[stringid_idx[i]];
	      l++;
	      k++;
	    }
	}
      for(;k<tfa->Njd;k++)
	tfa->bstore[threadid][tfa->stringid_idx[k]] = sqrt(-1);
    }

  rmsnew = sqrt(val2/l - val1*val1/l/l);

  /* Now if we're not doing TFA first, then subtract the signal from the light curve as well */
  if(!dotfafirst && !use_harm)
    {
      /* Get the signal */
      if(!use_bin)
	{
	  for(i=0;i<N;i++)
	    tfa->signal[threadid][i] = tfa->inputsignal[threadid][i];
	}
      else
	{
	  for(i=0;i<nbins;i++)
	    {
	      tfa->signal_bin_N[threadid][i] = 0;
	      tfa->signal_bin[threadid][i] = 0.;
	    }
	  if(use_period)
	    {
	      for(i=0;i<N;i++)
		{
		  ph = (t[i] - t[0]) / period;
		  ph -= floor(ph);
		  k = floor(ph * nbins);
		  tfa->signal_bin_ids[threadid][i] = k;
		  tfa->signal_bin_N[threadid][k]++;
		  tfa->signal_bin[threadid][k] += m[i];
		}
	    }
	  else
	    {
	      for(i=0;i<N;i++)
		{
		  ph = (t[i] - t[0]) / (t[N-1] + JDTOL - t[0]);
		  k = floor(ph * nbins);
		  tfa->signal_bin_ids[threadid][i] = k;
		  tfa->signal_bin_N[threadid][k]++;
		  tfa->signal_bin[threadid][k] += m[i];
		}
	    }
	  for(k=0;k<nbins;k++)
	    {
	      if(tfa->signal_bin_N[threadid][k] > 0)
		tfa->signal_bin[threadid][k] /= tfa->signal_bin_N[threadid][k];
	    }
	  for(i=0;i<N;i++)
	    tfa->signal[threadid][i] = tfa->signal_bin[threadid][tfa->signal_bin_ids[threadid][i]];
	}

      l = 0;
      val1 = 0.;
      val2 = 0.;
      if(!matchstringid)
	{
	  for(i=0,k=0; i < N; i++)
	    {
	      while((k < tfa->Njd ? (t[i] > tfa->JD[k] + JDTOL) : 0))
		{
		  tfa->b[threadid][k] = tfa->bstore[threadid][k];
		  k++;
		}
	      if ( k < tfa->Njd ? (t[i] < tfa->JD[k] + JDTOL && t[i] > tfa->JD[k] - JDTOL) : 0)
		{
		  tfa->b[threadid][k] = tfa->bstore[threadid][k] - tfa->signal[threadid][i];
		  val1 += tfa->b[threadid][k];
		  val2 += tfa->b[threadid][k]*tfa->b[threadid][k];
		  l++;
		  k++;
		}
	    }
	  for(;k<tfa->Njd;k++)
	    tfa->b[threadid][k] = tfa->bstore[threadid][k];
	}
      else
	{
	  for(i=0,k=0; i < N; i++)
	    {
	      while((k < tfa->Njd ? strncmp(stringid[stringid_idx[i]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) > 0 : 0))
		{
		  tfa->b[threadid][tfa->stringid_idx[k]] = tfa->bstore[threadid][tfa->stringid_idx[k]];
		  k++;
		}
	      if ( k < tfa->Njd ? !strncmp(stringid[stringid_idx[i]], tfa->stringid[tfa->stringid_idx[k]],MAXIDSTRINGLENGTH) : 0)
		{
		  tfa->b[threadid][tfa->stringid_idx[k]] = tfa->bstore[threadid][tfa->stringid_idx[k]] - tfa->signal[threadid][stringid_idx[i]];
		  val1 += tfa->b[threadid][tfa->stringid_idx[k]];
		  val2 += tfa->b[threadid][tfa->stringid_idx[k]]*tfa->b[threadid][tfa->stringid_idx[k]];
		  l++;
		  k++;
		}
	    }
	  for(;k<tfa->Njd;k++)
	    tfa->b[threadid][tfa->stringid_idx[k]] = tfa->bstore[threadid][tfa->stringid_idx[k]];
	}
      rmsnew = sqrt((val2 / l) - (val1*val1/l/l));
    }
  else if(!use_harm) {
      /* Get the signal */
      if(use_bin)
	{
	  for(i=0;i<nbins;i++)
	    {
	      tfa->signal_bin_N[threadid][i] = 0;
	      tfa->signal_bin[threadid][i] = 0.;
	    }
	  if(use_period)
	    {
	      for(i=0;i<N;i++)
		{
		  ph = (t[i] - t[0]) / period;
		  ph -= floor(ph);
		  k = floor(ph * nbins);
		  tfa->signal_bin_ids[threadid][i] = k;
		  tfa->signal_bin_N[threadid][k]++;
		  tfa->signal_bin[threadid][k] = 0.;
		}
	    }
	  else
	    {
	      for(i=0;i<N;i++)
		{
		  ph = (t[i] - t[0]) / (t[N-1] + JDTOL - t[0]);
		  k = floor(ph * nbins);
		  tfa->signal_bin_ids[threadid][i] = k;
		  tfa->signal_bin_N[threadid][k]++;
		  tfa->signal_bin[threadid][k] = 0.;
		}
	    }
	}
      for(i=0;i<tfa->Njd;i++)
	tfa->b[threadid][i] = tfa->bstore[threadid][i];
  }
  else
    {
      for(i=0;i<tfa->Njd;i++)
	tfa->b[threadid][i] = tfa->bstore[threadid][i];
    }

  zero_lc_averages(tfa->Njd,tfa->b[threadid],tfa->clipping);

  /* Prepare the decorr_trends matrix if needed */
  if(tfa->decorrflag)
    {
      l=0;
      lstart = 0;
      for(j=0;j<tfa->decorr_Nlcterms;j++)
	{
	  lstart = l;
	  if(!matchstringid)
	    {
	      for(i=0,k=0; i < N; i++)
		{
		  while((k < tfa->Njd ? (t[i] > tfa->JD[k] + JDTOL) : 0))
		    {
		      l = lstart;
		      tfa->decorr_trends[threadid][k][l] = 0.;
		      l++;
		      for(h=2; h <= tfa->decorr_lc_order[j]; h++)
			{
			  tfa->decorr_trends[threadid][k][l] = 0.;
			  l++;
			}
		      k++;
		    }
		  if ( k < tfa->Njd ? (t[i] < tfa->JD[k] + JDTOL && t[i] > tfa->JD[k] - JDTOL) : 0)
		    {
		      l = lstart;
		      tfa->decorr_trends[threadid][k][l] = tfa->lcdecorr_terms_in[j][lcindex][i];
		      term = tfa->decorr_trends[threadid][k][l];
		      l++;
		      for(h=2; h <= tfa->decorr_lc_order[j]; h++)
			{
			  tfa->decorr_trends[threadid][k][l] = tfa->decorr_trends[threadid][k][l-1]*term;
			  l++;
			}
		      k++;
		    }
		}
	      for(;k<tfa->Njd;k++)
		{
		  l = lstart;
		  tfa->decorr_trends[threadid][k][l] = 0.;
		  l++;
		  for(h=2; h <= tfa->decorr_lc_order[j]; h++)
		    {
		      tfa->decorr_trends[threadid][k][l] = 0.;
		      l++;
		    }
		}
	    }
	  else
	    {
	      for(i=0,k=0; i < N; i++)
		{
		  while((k < tfa->Njd ? strncmp(stringid[stringid_idx[i]],tfa->stringid[tfa->stringid_idx[k]],MAXIDSTRINGLENGTH) > 0 : 0))
		    {
		      l = lstart;
		      tfa->decorr_trends[threadid][tfa->stringid_idx[k]][l] = 0.;
		      l++;
		      for(h=2; h <= tfa->decorr_lc_order[j]; h++)
			{
			  tfa->decorr_trends[threadid][tfa->stringid_idx[k]][l] = 0.;
			  l++;
			}
		      k++;
		    }
		  if ( k < tfa->Njd ? !strncmp(stringid[stringid_idx[i]],tfa->stringid[tfa->stringid_idx[k]],MAXIDSTRINGLENGTH) : 0)
		    {
		      l = lstart;
		      tfa->decorr_trends[threadid][tfa->stringid_idx[k]][l] = tfa->lcdecorr_terms_in[j][lcindex][stringid_idx[i]];
		      term = tfa->decorr_trends[threadid][tfa->stringid_idx[k]][l];
		      l++;
		      for(h=2; h <= tfa->decorr_lc_order[j]; h++)
			{
			  tfa->decorr_trends[threadid][tfa->stringid_idx[k]][l] = tfa->decorr_trends[threadid][tfa->stringid_idx[k]][l-1]*term;
			  l++;
			}
		      k++;
		    }
		}
	      for(;k<tfa->Njd;k++)
		{
		  l = lstart;
		  tfa->decorr_trends[threadid][tfa->stringid_idx[k]][l] = 0.;
		  l++;
		  for(h=2; h <= tfa->decorr_lc_order[j]; h++)
		    {
		      tfa->decorr_trends[threadid][tfa->stringid_idx[k]][l] = 0.;
		      l++;
		    }
		}
	    }
	}
      zero_trend_averages(tfa->Njd, tfa->Ndecorr, tfa->decorr_trends[threadid], tfa->clipping);
    }

  isatrend_flag = 0;
  for(i=0;i<tfa->Ntrends;i++)
    {
      if((!strncmp(lc_name,tfa->trend_names[i],strlen(tfa->trend_names[i])) && strlen(lc_name) == strlen(tfa->trend_names[i])) || sqrt((((lcx - tfa->trendx[i])*(lcx - tfa->trendx[i])) + ((lcy - tfa->trendy[i])*(lcy - tfa->trendy[i])))) < tfa->pixelsep)
	isatrend_flag = 1;
    }
  if(!isatrend_flag && (tfa->Ntrends == tfa->Ntfatot))
    {
      Ntrends_ = tfa->Ntrends;
      u_ = tfa->u;
      v_ = tfa->v;
      w_ = tfa->w1;
    }
  else
    {
      l = 0;
      if(isatrend_flag)
	{
	  for(i=0;i<tfa->Ntrends;i++)
	    {
	      if(!(!strncmp(lc_name,tfa->trend_names[i],strlen(tfa->trend_names[i])) && strlen(lc_name) == strlen(tfa->trend_names[i])) && !(sqrt((((lcx - tfa->trendx[i])*(lcx - tfa->trendx[i])) + ((lcy - tfa->trendy[i])*(lcy - tfa->trendy[i])))) < tfa->pixelsep))
		{
		  for(k=0;k<tfa->Njd;k++)
		    {
		      tfa->u2[threadid][k][l] = tfa->trends[k][i];
		    }
		  l++;
		}
	    }
	}
      else
	{
	  for(i=0;i<tfa->Ntrends;i++)
	    {
	      for(k=0;k<tfa->Njd;k++)
		{
		  tfa->u2[threadid][k][l] = tfa->trends[k][i];
		}
	      l++;
	    }
	}
      if(tfa->Ntfatot > tfa->Ntrends)
	{
	  /* We're doing simultaneous decorrelation */
	  if(tfa->decorrflag ? !tfa->decorr_iterate : 0)
	    {
	      for(i=0;i<tfa->Ndecorr;i++)
		{
		  for(k=0;k<tfa->Njd;k++)
		    {
		      tfa->u2[threadid][k][l] = tfa->decorr_trends[threadid][k][i];
		    }
		  l++;
		}
	    }
	  /* We're fitting a fourier series to the light curve */
	  if(use_harm)
	    {
	      for(k=0;k<tfa->Njd;k++)
		{
		  if(!isnan(tfa->b[threadid][k])) {
		    tfa->u2[threadid][k][l] = cos((tfa->JD[k])* 2.0 * M_PI / period);
		    tfa->u2[threadid][k][l+1] = sin((tfa->JD[k])* 2.0 * M_PI / period);
		  } else {
		    tfa->u2[threadid][k][l] = 0.;
		    tfa->u2[threadid][k][l+1] = 0.;
		  }
		  tfa->harmterm[threadid][k][0] = tfa->u2[threadid][k][l];
		  tfa->harmterm[threadid][k][1] = tfa->u2[threadid][k][l+1];
		}
	      l += 2;
	      for(i=0;i<Nharm;i++)
		{
		  for(k=0;k<tfa->Njd;k++)
		    {
		      if(!isnan(tfa->b[threadid][k])) {
			tfa->u2[threadid][k][l] = cos((tfa->JD[k])* 2.0 * M_PI * ((double) (i+2))/period);
			tfa->u2[threadid][k][l+1] = sin((tfa->JD[k])* 2.0 * M_PI * ((double) (i+2))/period);
		      } else {
			tfa->u2[threadid][k][l] = 0.;
			tfa->u2[threadid][k][l+1] = 0.;
		      }
		      tfa->harmterm[threadid][k][2*(i + 1)] = tfa->u2[threadid][k][l];
		      tfa->harmterm[threadid][k][2*(i + 1) + 1] = tfa->u2[threadid][k][l+1];
		    }
		  l += 2;
		}
	      for(i=0;i<Nsubharm;i++)
		{
		  for(k=0;k<tfa->Njd;k++)
		    {
		      if(!isnan(tfa->b[threadid][k])) {
			tfa->u2[threadid][k][l] = cos((tfa->JD[k])* 2.0 * M_PI / period / ((double) (i + 2)));
			tfa->u2[threadid][k][l+1] = sin((tfa->JD[k])* 2.0 * M_PI / period / ((double) (i + 2)));
		      } else {
			tfa->u2[threadid][k][l] = 0.;
			tfa->u2[threadid][k][l+1] = 0.;
		      }
		      tfa->harmterm[threadid][k][2*(i + 1 + Nharm)] = tfa->u2[threadid][k][l];
		      tfa->harmterm[threadid][k][2*(i + 1 + Nharm) + 1] = tfa->u2[threadid][k][l+1];
		    }
		  l += 2;
		}
	    }
	}

      svdcmp(tfa->u2[threadid],tfa->Njd,l,tfa->w2[threadid],tfa->v2[threadid]);
      Ntrends_ = l;
      u_ = tfa->u2[threadid];
      v_ = tfa->v2[threadid];
      w_ = tfa->w2[threadid];
    }

  if(tfa->decorrflag ? tfa->decorr_iterate : 0)
    {
      /* We're doing iterative decorrelation */
      for(i=0;i<tfa->Ndecorr;i++)
	{
	  for(k=0;k<tfa->Njd;k++)
	    {
	      tfa->u_decorr[threadid][k][i] = tfa->decorr_trends[threadid][k][i];
	    }
	}
      svdcmp(tfa->u_decorr[threadid],tfa->Njd,tfa->Ndecorr,tfa->w1_decorr[threadid],tfa->v_decorr[threadid]);
    }

  iter = 0;
  /* The main TFA-SR iteration loop */
  do
    {
      rmsold = rmsnew;

      for(i=0;i<N;i++)
	tfa->m_out[threadid][i] = 0.;
      /* Apply the per-light curve decorrelation first if we're doing that */
      if(tfa->decorrflag ? tfa->decorr_iterate : 0)
	{
	  /* Do the back-substitution */
	  svbksb(tfa->u_decorr[threadid],tfa->w1_decorr[threadid],tfa->v_decorr[threadid],tfa->Njd,tfa->Ndecorr,tfa->b[threadid],tfa->a_decorr[threadid]);

	  /* Get the new b vector */
	  l=0;
	  for(i=0;i<tfa->Ndecorr;i++)
	    {
	      k = 0;
	      if(!matchstringid)
		{
		  for(h = 0; h<N; h++)
		    {
		      while((k < tfa->Njd ? (t[h] > tfa->JD[k] + JDTOL) : 0)) k++;
		      if(k < tfa->Njd ? (t[h] < tfa->JD[k] + JDTOL && t[h] > tfa->JD[k] - JDTOL) : 0)
			{
			  tfa->m_out[threadid][h] += tfa->a_decorr[threadid][i]*tfa->decorr_trends[threadid][k][i];
			  k++;
			}
		    }
		}
	      else
		{
		  for(h = 0; h<N; h++)
		    {
		      while((k < tfa->Njd ? strncmp(stringid[stringid_idx[h]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) > 0 : 0)) k++;
		      if(k < tfa->Njd ? !strncmp(stringid[stringid_idx[h]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) : 0)
			{
			  tfa->m_out[threadid][stringid_idx[h]] += tfa->a_decorr[threadid][i]*tfa->decorr_trends[threadid][tfa->stringid_idx[k]][i];
			  k++;
			}
		    }
		}
	      l++;
	    }

	  if(!matchstringid)
	    {
	      for(i=0,k=0; i < N; i++)
		{
		  while((k < tfa->Njd ? (t[i] > tfa->JD[k] + JDTOL) : 0))
		    {
		      k++;
		    }
		  if ( k < tfa->Njd ? (t[i] < tfa->JD[k] + JDTOL && t[i] > tfa->JD[k] - JDTOL) : 0)
		    {
		      tfa->b[threadid][k] = tfa->b[threadid][k] - tfa->m_out[threadid][i];
		      k++;
		    }
		}
	    }
	  else
	    {
	      for(i=0,k=0; i < N; i++)
		{
		  while((k < tfa->Njd ? strncmp(stringid[stringid_idx[i]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) > 0 : 0))
		    {
		      k++;
		    }
		  if ( k < tfa->Njd ? !strncmp(stringid[stringid_idx[i]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) : 0)
		    {
		      tfa->b[threadid][tfa->stringid_idx[k]] = tfa->b[threadid][tfa->stringid_idx[k]] - tfa->m_out[threadid][stringid_idx[i]];
		      k++;
		    }
		}
	    }
	}

      svbksb(u_,w_,v_,tfa->Njd,Ntrends_,tfa->b[threadid],tfa->a[threadid]);

      l = 0;
      /* Correct the light curve */
      for(i=0;i<tfa->Ntrends;i++)
	{
	  if(!(!strncmp(lc_name,tfa->trend_names[i],strlen(tfa->trend_names[i])) && strlen(lc_name) == strlen(tfa->trend_names[i])) && !(sqrt((((lcx - tfa->trendx[i])*(lcx - tfa->trendx[i])) + ((lcy - tfa->trendy[i])*(lcy - tfa->trendy[i])))) < tfa->pixelsep))
	    {
	      k = 0;
	      if(!matchstringid)
		{
		  for(h = 0; h<N; h++)
		    {
		      while((k < tfa->Njd ? (t[h] > tfa->JD[k] + JDTOL) : 0)) k++;
		      if(k < tfa->Njd ? (t[h] < tfa->JD[k] + JDTOL && t[h] > tfa->JD[k] - JDTOL) : 0)
			{
			  tfa->m_out[threadid][h] += tfa->a[threadid][l]*tfa->trends[k][i];
			  k++;
			}
		    }
		}
	      else
		{
		  for(h = 0; h<N; h++)
		    {
		      while((k < tfa->Njd ? strncmp(stringid[stringid_idx[h]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) > 0 : 0)) k++;
		      if(k < tfa->Njd ? !strncmp(stringid[stringid_idx[h]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) : 0)
			{
			  tfa->m_out[threadid][stringid_idx[h]] += tfa->a[threadid][l]*tfa->trends[tfa->stringid_idx[k]][i];
			  k++;
			}
		    }
		}
	      l++;
	    }
	}
      /* Do the simultaneous trends */
      i=l;
      l=0;
      if(tfa->decorrflag ? !tfa->decorr_iterate : 0)
	{
	  for(l=0;l<tfa->Ndecorr;l++,i++)
	    {
	      k = 0;
	      if(!matchstringid)
		{
		  for(h = 0; h<N; h++)
		    {
		      while((k < tfa->Njd ? (t[h] > tfa->JD[k] + JDTOL) : 0)) k++;
		      if(k < tfa->Njd ? (t[h] < tfa->JD[k] + JDTOL && t[h] > tfa->JD[k] - JDTOL) : 0)
			{
			  tfa->m_out[threadid][h] += tfa->a[threadid][i]*tfa->decorr_trends[threadid][k][l];
			  k++;
			}
		    }
		}
	      else
		{
		  for(h = 0; h<N; h++)
		    {
		      while((k < tfa->Njd ? strncmp(stringid[stringid_idx[h]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) > 0 : 0)) k++;
		      if(k < tfa->Njd ? !strncmp(stringid[stringid_idx[h]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) : 0)
			{
			  tfa->m_out[threadid][stringid_idx[h]] += tfa->a[threadid][i]*tfa->decorr_trends[threadid][tfa->stringid_idx[k]][l];
			  k++;
			}
		    }
		}
	    }
	}
      /* Do the harmonic terms */
      if(use_harm)
	{
	  k = 0;
	  if(!matchstringid)
	    {
	      for(h = 0; h<N; h++)
		{
		  while((k < tfa->Njd ? (t[h] > tfa->JD[k] + JDTOL) : 0)) k++;
		  if(k < tfa->Njd ? (t[h] < tfa->JD[k] + JDTOL && t[h] > tfa->JD[k] - JDTOL) : 0)
		    {
		      for(l=0;l<2*(1 + Nharm + Nsubharm);l++)
			{
			  tfa->m_out[threadid][h] += tfa->a[threadid][i+l]*tfa->harmterm[threadid][k][l];
			}
		      k++;
		    }
		}
	    }
	  else
	    {
	      for(h = 0; h<N; h++)
		{
		  while((k < tfa->Njd ? strncmp(stringid[stringid_idx[h]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) > 0 : 0)) k++;
		  if(k < tfa->Njd ? !strncmp(stringid[stringid_idx[h]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) : 0)
		    {
		      for(l=0;l<2*(1 + Nharm + Nsubharm);l++)
			{
			  tfa->m_out[threadid][stringid_idx[h]] += tfa->a[threadid][i+l]*tfa->harmterm[threadid][tfa->stringid_idx[k]][l];
			}
		      k++;
		    }
		}
	    }
	}

      for(i=0; i<N; i++)
	tfa->m_out[threadid][i] = m[i] - tfa->m_out[threadid][i];

      /* Get the signal */
      if(!use_bin && !use_harm)
	{
	  val1 = 0.;
	  val2 = 0.;
	  val3 = 0.;
	  val4 = 0.;
	  l = 0;
	  for(i=0;i<N;i++)
	    {
	      if(tfa->inputsignal[threadid][i] != 0.)
		{
		  val1 += tfa->inputsignal[threadid][i];
		  val2 += tfa->inputsignal[threadid][i]*tfa->inputsignal[threadid][i];
		  val3 += tfa->inputsignal[threadid][i]*tfa->m_out[threadid][i];
		  val4 += tfa->m_out[threadid][i];
		  l++;
		}
	    }
	  if(l > 0)
	    {
	      b = (val2 - val1*val1/l)/(val3 - val4*val1/l);
	      a = val4/l - b*val1/l;
	    }
	  for(i=0;i<N;i++)
	    tfa->signal[threadid][i] = a + b*tfa->inputsignal[threadid][i];
	}
      else if(!use_harm)
	{
	  for(k=0;k<nbins;k++)
	    {
	      tfa->signal_bin[threadid][k] = 0.;
	    }
	  for(i=0;i<N;i++)
	    {
	      k = tfa->signal_bin_ids[threadid][i];
	      tfa->signal_bin[threadid][k] += tfa->m_out[threadid][i];
	    }
	  for(k=0;k<nbins;k++)
	    {
	      if(tfa->signal_bin_N[threadid][k] > 0)
		tfa->signal_bin[threadid][k] /= tfa->signal_bin_N[threadid][k];
	    }
	  for(i=0;i<N;i++)
	    tfa->signal[threadid][i] = tfa->signal_bin[threadid][tfa->signal_bin_ids[threadid][i]];
	}

      if(!use_harm)
	{
	  l = 0;
	  val1 = 0.;
	  val2 = 0.;
	  if(!matchstringid)
	    {
	      for(i=0,k=0; i < N; i++)
		{
		  while((k < tfa->Njd ? (t[i] > tfa->JD[k] + JDTOL) : 0))
		    {
		      tfa->b[threadid][k] = tfa->bstore[threadid][k];
		      k++;
		    }
		  if ( k < tfa->Njd ? (t[i] < tfa->JD[k] + JDTOL && t[i] > tfa->JD[k] - JDTOL) : 0)
		    {
		      tfa->b[threadid][k] = tfa->bstore[threadid][k] - tfa->signal[threadid][i];
		      val1 += tfa->b[threadid][k];
		      val2 += tfa->b[threadid][k]*tfa->b[threadid][k];
		      l++;
		      k++;
		    }
		}
	      for(;k<tfa->Njd;k++)
		tfa->b[threadid][k] = tfa->bstore[threadid][k];
	    }
	  else
	    {
	      for(i=0,k=0; i < N; i++)
		{
		  while((k < tfa->Njd ? strncmp(stringid[stringid_idx[i]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) > 0 : 0))
		    {
		      tfa->b[threadid][tfa->stringid_idx[k]] = tfa->bstore[threadid][tfa->stringid_idx[k]];
		      k++;
		    }
		  if ( k < tfa->Njd ? !strncmp(stringid[stringid_idx[i]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) : 0)
		    {
		      tfa->b[threadid][tfa->stringid_idx[k]] = tfa->bstore[threadid][tfa->stringid_idx[k]] - tfa->signal[threadid][stringid_idx[i]];
		      val1 += tfa->b[threadid][tfa->stringid_idx[k]];
		      val2 += tfa->b[threadid][tfa->stringid_idx[k]]*tfa->b[threadid][tfa->stringid_idx[k]];
		      l++;
		      k++;
		    }
		}
	      for(;k<tfa->Njd;k++)
		tfa->b[threadid][tfa->stringid_idx[k]] = tfa->bstore[threadid][tfa->stringid_idx[k]];
	    }
	  rmsnew = sqrt((val2 / l) - (val1*val1/l/l));
	  zero_lc_averages(tfa->Njd,tfa->b[threadid],tfa->clipping);
	  iter++;
	}
    } while(!use_harm ? fabs(rmsold - rmsnew) < tfathresh*rmsold && iter < maxiter : 0);


  /* Write out the coefficients if we're doing that sort of thing */
  if(coeff_flag)
    {
      if((coeff_file = fopen(coeff_file_name,"w")) == NULL)
	{
	  error2(ERR_CANNOTWRITE,coeff_file_name);
	}
      l=0;
      for(i=0; i<tfa->Ntrends;i++)
	{
	  if(!(!strncmp(lc_name,tfa->trend_names[i],strlen(tfa->trend_names[i])) && strlen(lc_name) == strlen(tfa->trend_names[i])) && !(sqrt((((lcx - tfa->trendx[i])*(lcx - tfa->trendx[i])) + ((lcy - tfa->trendy[i])*(lcy - tfa->trendy[i])))) < tfa->pixelsep))
	    {
	      fprintf(coeff_file,"%s %f\n",tfa->trend_names[i],tfa->a[threadid][l]);
	      l++;
	    }
	}
      i=l;
      l=0;
      for(;i < Ntrends_;i++)
	{
	  fprintf(coeff_file,"decorr_term_%d %f\n",l,tfa->a[threadid][i]);
	  l++;
	}
      if(tfa->decorrflag ? tfa->decorr_iterate : 0)
	{
	  for(i=0;i<tfa->Ndecorr;i++)
	    {
	      fprintf(coeff_file,"decorr_term_%d %f\n",i,tfa->a_decorr[threadid][i]);
	    }
	}
      fclose(coeff_file);
    }
  /* Open the output light curve, correct the light curve and write it out */
  if(outlc)
    {
      if((lcout = fopen(lc_out_name,"w")) == NULL)
	error2(ERR_CANNOTWRITE,lc_out_name);
    }
  l=0;
  for(i=0;i<N;i++)
    tfa->m_out[threadid][i] = 0.;
  if(use_harm)
    {
      for(i=0;i<N;i++)
	tfa->signal[threadid][i] = 0.;
    }

  if(correctlc || outlc)
    {
      if(tfa->decorrflag ? tfa->decorr_iterate : 0)
	{
	  for(i=0;i<tfa->Ndecorr;i++)
	    {
	      k = 0;
	      if(!matchstringid)
		{
		  for(h = 0; h<N; h++)
		    {
		      while((k < tfa->Njd ? (t[h] > tfa->JD[k] + JDTOL) : 0)) k++;
		      if(k < tfa->Njd ? (t[h] < tfa->JD[k] + JDTOL && t[h] > tfa->JD[k] - JDTOL) : 0)
			{
			  tfa->m_out[threadid][h] += tfa->a_decorr[threadid][i]*tfa->decorr_trends[threadid][k][i];
			  k++;
			}
		    }
		}
	      else
		{
		  for(h = 0; h<N; h++)
		    {
		      while((k < tfa->Njd ? strncmp(stringid[stringid_idx[h]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) > 0 : 0)) k++;
		      if(k < tfa->Njd ? !strncmp(stringid[stringid_idx[h]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) : 0)
			{
			  tfa->m_out[threadid][stringid_idx[h]] += tfa->a_decorr[threadid][i]*tfa->decorr_trends[threadid][tfa->stringid_idx[k]][i];
			  k++;
			}
		    }
		}
	    }
	}

      l=0;
      for(i=0;i<tfa->Ntrends;i++)
	{
	  if(!(!strncmp(lc_name,tfa->trend_names[i],strlen(tfa->trend_names[i])) && strlen(lc_name) == strlen(tfa->trend_names[i])) && !(sqrt((((lcx - tfa->trendx[i])*(lcx - tfa->trendx[i])) + ((lcy - tfa->trendy[i])*(lcy - tfa->trendy[i])))) < tfa->pixelsep))
	    {
	      k = 0;
	      if(!matchstringid)
		{
		  for(h = 0; h<N; h++)
		    {
		      while((k < tfa->Njd ? (t[h] > tfa->JD[k] + JDTOL) : 0)) k++;
		      if(k < tfa->Njd ? (t[h] < tfa->JD[k] + JDTOL && t[h] > tfa->JD[k] - JDTOL) : 0)
			{
			  tfa->m_out[threadid][h] += tfa->a[threadid][l]*tfa->trends[k][i];
			  k++;
			}
		    }
		}
	      else
		{
		  for(h = 0; h<N; h++)
		    {
		      while((k < tfa->Njd ? strncmp(stringid[stringid_idx[h]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) > 0 : 0)) k++;
		      if(k < tfa->Njd ? !strncmp(stringid[stringid_idx[h]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) : 0)
			{
			  tfa->m_out[threadid][stringid_idx[h]] += tfa->a[threadid][l]*tfa->trends[tfa->stringid_idx[k]][i];
			  k++;
			}
		    }
		}
	      l++;
	    }
	}

      /* Do the simultaneous trends */
      i=l;
      l=0;
      if(tfa->decorrflag ? !tfa->decorr_iterate : 0)
	{
	  for(;i<tfa->Ntrends + tfa->Ndecorr;i++)
	    {
	      k = 0;
	      if(!matchstringid)
		{
		  for(h = 0; h<N; h++)
		    {
		      while((k < tfa->Njd ? (t[h] > tfa->JD[k] + JDTOL) : 0)) k++;
		      if(k < tfa->Njd ? (t[h] < tfa->JD[k] + JDTOL && t[h] > tfa->JD[k] - JDTOL) : 0)
			{
			  tfa->m_out[threadid][h] += tfa->a[threadid][i]*tfa->decorr_trends[threadid][k][l];
			  k++;
			}
		    }
		}
	      else
		{
		  for(h = 0; h<N; h++)
		    {
		      while((k < tfa->Njd ? strncmp(stringid[stringid_idx[h]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) > 0 : 0)) k++;
		      if(k < tfa->Njd ? !strncmp(stringid[stringid_idx[h]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) : 0)
			{
			  tfa->m_out[threadid][stringid_idx[h]] += tfa->a[threadid][i]*tfa->decorr_trends[threadid][tfa->stringid_idx[k]][l];
			  k++;
			}
		    }
		}
	      l++;
	    }
	}
      if(use_harm)
	{
	  /* Do the harmonic terms, we'll put these in the signal vector rather than include them in the TFA model for output. */
	  k = 0;
	  if(!matchstringid)
	    {
	      for(h = 0; h<N; h++)
		{
		  while((k < tfa->Njd ? (t[h] > tfa->JD[k] + JDTOL) : 0)) k++;
		  if(k < tfa->Njd ? (t[h] < tfa->JD[k] + JDTOL && t[h] > tfa->JD[k] - JDTOL) : 0)
		    {
		      for(l=0;l<2*(1 + Nharm + Nsubharm);l++)
			{
			  tfa->signal[threadid][h] += tfa->a[threadid][i+l]*tfa->harmterm[threadid][k][l];
			}
		      k++;
		    }
		}
	    }
	  else
	    {
	      for(h = 0; h<N; h++)
		{
		  while((k < tfa->Njd ? strncmp(stringid[stringid_idx[h]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) > 0 : 0)) k++;
		  if(k < tfa->Njd ? !strncmp(stringid[stringid_idx[h]], tfa->stringid[tfa->stringid_idx[k]], MAXIDSTRINGLENGTH) : 0)
		    {
		      for(l=0;l<2*(1 + Nharm + Nsubharm);l++)
			{
			  tfa->signal[threadid][stringid_idx[h]] += tfa->a[threadid][i+l]*tfa->harmterm[threadid][tfa->stringid_idx[k]][l];
			}
		      k++;
		    }
		}
	    }
	}


      for(i=0;i<N;i++)
	{
	  if(outlc)
	    {
	      fprintf(lcout,"%f %f %f %f %f\n",t[i],m[i],tfa->m_out[threadid][i],tfa->signal[threadid][i],e[i]);
	    }
	  if(correctlc)
	    {
	      m[i] = m[i] - tfa->m_out[threadid][i];
	    }
	}
      if(outlc)
	fclose(lcout);
    }
  if(correctlc)
    {
      n = 0;
      val1 = 0.;
      val2 = 0.;
      for(i=0;i<N;i++)
	{
	  n++;
	  val1 += m[i]*m[i];
	  val2 += m[i];
	}
    }
  else
    {
      n = 0;
      val1 = 0.;
      val2 = 0.;
      for(i = 0; i<N;i++)
	{
	  n++;
	  val1 += (m[i] - tfa->m_out[threadid][i])*(m[i] - tfa->m_out[threadid][i]);
	  val2 += (m[i] - tfa->m_out[threadid][i]);
	}
    }
  *ave_out = (double) (val2 / (double) n);
  *rms_out = sqrt((double) ((val1 / (double) n) - (val2 * val2 / (double) (n*n))));
}
