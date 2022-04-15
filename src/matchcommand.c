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
#include "commands.h"
#include "programdata.h"
#include "functions.h"

#ifdef USECFITSIO
#include "fitsio.h"
#endif


#define VARTOOLS_MATCHCOMMAND_MISSINGMETHOD_CULLMISSING 0
#define VARTOOLS_MATCHCOMMAND_MISSINGMETHOD_NANMISSING 1
#define VARTOOLS_MATCHCOMMAND_MISSINGMETHOD_MISSINGVAL 2

void FreeMatchDataVectors(_MatchCommand *m, _MatchData *md) {
  int i, j;
  double **dblptr;
  short **shortptr;
  int **intptr;
  char **charptr;
  char ***stringptr;
  float **floatptr;
  long **longptr;
  if(md == NULL)
    return;
  for(i=0; i < m->Naddvars+1; i++) {
    if(md[i].dataptr == NULL)
      continue;
    switch(md[i].datatype) {
    case VARTOOLS_TYPE_DOUBLE:
    case VARTOOLS_TYPE_CONVERTJD:
      dblptr = (double **) md[i].dataptr;
      free(*dblptr);
      free(dblptr);
      md[i].dataptr = NULL;
      break;
    case VARTOOLS_TYPE_STRING:
      stringptr = (char ***) md[i].dataptr;
      for(j=0; j < md[i].sizevec; j++) {
	if(stringptr[0][j] != NULL)
	  free(stringptr[0][j]);
	stringptr[0][j] = NULL;
      }
      free(*stringptr);
      free(stringptr);
      md[i].dataptr = NULL;
      break;
    case VARTOOLS_TYPE_INT:
      intptr = (int **) md[i].dataptr;
      free(*intptr);
      free(intptr);
      md[i].dataptr = NULL;
      break;
    case VARTOOLS_TYPE_SHORT:
      shortptr = (short **) md[i].dataptr;
      free(*shortptr);
      free(shortptr);
      md[i].dataptr = NULL;
      break;
    case VARTOOLS_TYPE_FLOAT:
      floatptr = (float **) md[i].dataptr;
      free(*floatptr);
      free(floatptr);
      md[i].dataptr = NULL;
      break;
    case VARTOOLS_TYPE_LONG:
      longptr = (long **) md[i].dataptr;
      free(*longptr);
      free(longptr);
      md[i].dataptr = NULL;
      break;
    case VARTOOLS_TYPE_CHAR:
      charptr = (char **) md[i].dataptr;
      free(*charptr);
      free(charptr);
      md[i].dataptr = NULL;
      break;
    default:
      error(ERR_BADTYPE);
    }
    md[i].sizevec = 0;
    md[i].Npoints = 0;
  }
}

void FreeMatchData(_MatchCommand *m, _MatchData *md) {
  FreeMatchDataVectors(m, md);
  free(md);
}
  
void MemAllocMatchData(_MatchData *md, int Nv, int Np) {
  int i, j;

  double **dblptr;
  short **shortptr;
  int **intptr;
  char **charptr;
  char ***stringptr;
  float **floatptr;
  long **longptr;

  if(Np <= 0) return;
  for(i=0; i < Nv; i++) {
    if(Np <= md[i].sizevec)
      continue;
    if(!md[i].sizevec) {
      md[i].sizevec = 1024;
      while(Np > md[i].sizevec) 
	md[i].sizevec = md[i].sizevec*2;
      switch(md[i].datatype) {
      case VARTOOLS_TYPE_DOUBLE:
      case VARTOOLS_TYPE_CONVERTJD:
	dblptr = (double **) malloc(sizeof(double *));
	if((dblptr[0] = (double *) malloc(md[i].sizevec*sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
	md[i].dataptr = (void *) dblptr;
	break;
      case VARTOOLS_TYPE_STRING:
	stringptr = (char ***) malloc(sizeof(double **));
	if((stringptr[0] = (char **) malloc(md[i].sizevec*sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < md[i].sizevec; j++) {
	  if((stringptr[0][j] = (char *) malloc(MAXLEN*sizeof(char))) == NULL)
	    error(ERR_MEMALLOC);
	}
	md[i].dataptr = (void *) stringptr;
	break;
      case VARTOOLS_TYPE_INT:
	intptr = (int **) malloc(sizeof(int *));
	if((intptr[0] = (int *) malloc(md[i].sizevec*sizeof(int))) == NULL)
	  error(ERR_MEMALLOC);
	md[i].dataptr = (void *) intptr;
	break;
      case VARTOOLS_TYPE_SHORT:
	shortptr = (short **) malloc(sizeof(short *));
	if((shortptr[0] = (short *) malloc(md[i].sizevec*sizeof(short))) == NULL)
	  error(ERR_MEMALLOC);
	md[i].dataptr = (void *) shortptr;
	break;
      case VARTOOLS_TYPE_FLOAT:
	floatptr = (float **) malloc(sizeof(float *));
	if((floatptr[0] = (float *) malloc(md[i].sizevec*sizeof(float))) == NULL)
	  error(ERR_MEMALLOC);
	md[i].dataptr = (void *) floatptr;
	break;
      case VARTOOLS_TYPE_LONG:
	longptr = (long **) malloc(sizeof(long *));
	if((longptr[0] = (long *) malloc(md[i].sizevec*sizeof(long))) == NULL)
	  error(ERR_MEMALLOC);
	md[i].dataptr = (void *) longptr;
	break;
      case VARTOOLS_TYPE_CHAR:
	charptr = (char **) malloc(sizeof(char *));
	if((charptr[0] = (char *) malloc(md[i].sizevec*sizeof(char))) == NULL)
	  error(ERR_MEMALLOC);
	md[i].dataptr = (void *) charptr;
      default:
	error(ERR_BADTYPE);
      }
    } else {
      while(Np > md[i].sizevec)
	md[i].sizevec = md[i].sizevec*2;
      switch(md[i].datatype) {
      case VARTOOLS_TYPE_DOUBLE:
      case VARTOOLS_TYPE_CONVERTJD:
	dblptr = (double **) md[i].dataptr;
	if((dblptr[0] = (double *) realloc(dblptr[0], md[i].sizevec*sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
	md[i].dataptr = (void *) dblptr;
	break;
      case VARTOOLS_TYPE_STRING:
	stringptr = (char ***) md[i].dataptr;
	if((stringptr[0] = (char **) realloc(stringptr[0],md[i].sizevec*sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=md[i].sizevec/2; j < md[i].sizevec; j++) {
	  if((stringptr[0][j] = (char *) malloc(MAXLEN*sizeof(char))) == NULL)
	    error(ERR_MEMALLOC);
	}
	md[i].dataptr = (void *) stringptr;
	break;
      case VARTOOLS_TYPE_INT:
	intptr = (int **) md[i].dataptr;
	if((intptr[0] = (int *) realloc(intptr[0],md[i].sizevec*sizeof(int))) == NULL)
	  error(ERR_MEMALLOC);
	md[i].dataptr = (void *) intptr;
	break;
      case VARTOOLS_TYPE_SHORT:
	shortptr = (short **) md[i].dataptr;
	if((shortptr[0] = (short *) realloc(shortptr[0],md[i].sizevec*sizeof(short))) == NULL)
	  error(ERR_MEMALLOC);
	md[i].dataptr = (void *) shortptr;
	break;
      case VARTOOLS_TYPE_FLOAT:
	floatptr = (float **) md[i].dataptr;
	if((floatptr[0] = (float *) realloc(floatptr[0],md[i].sizevec*sizeof(float))) == NULL)
	  error(ERR_MEMALLOC);
	md[i].dataptr = (void *) floatptr;
	break;
      case VARTOOLS_TYPE_LONG:
	longptr = (long **) md[i].dataptr;
	if((longptr[0] = (long *) realloc(longptr[0],md[i].sizevec*sizeof(long))) == NULL)
	  error(ERR_MEMALLOC);
	md[i].dataptr = (void *) longptr;
	break;
      case VARTOOLS_TYPE_CHAR:
	charptr = (char **) md[i].dataptr;
	if((charptr[0] = (char *) realloc(charptr[0],md[i].sizevec*sizeof(char))) == NULL)
	  error(ERR_MEMALLOC);
	md[i].dataptr = (void *) charptr;
      default:
	error(ERR_BADTYPE);
      }
    }
  }
}

FILE *ExecMatchOpenCommand(char *open_command_str, char *infilename)
/* This function uses popen to carry out a requested shell command in
   opening a file */
{
  char *execcommand;
  int size_execcommand = 2048;
  int i, i2;
  int lc_in_name_length;
  FILE *return_pipe;
  lc_in_name_length = strlen(infilename);
  if((execcommand = (char *) malloc((size_execcommand+1))) == NULL)
    error(ERR_MEMALLOC);
  i = 0; i2 = 0;

  /* Parse the command string, substituting the input file name as needed */
  while(open_command_str[i2] != '\0')
    {
      if(open_command_str[i2] != '%')
	{
	  if(i >= (size_execcommand)) {
	    size_execcommand *= 2;
	    if((execcommand = (char *) realloc(execcommand, (size_execcommand+1))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  execcommand[i] = open_command_str[i2];
	  i++;
	  execcommand[i] = '\0';
	  i2++;
	}
      else
	{
	  i2++;
	  if(open_command_str[i2] == 's')
	    {
	      i2++;
	      if((i + lc_in_name_length) >= (size_execcommand)) {
		while((i + lc_in_name_length) >= (size_execcommand))
		  size_execcommand *= 2;
		if((execcommand = (char *) realloc(execcommand, (size_execcommand+1))) == NULL)
		  error(ERR_MEMALLOC);
	      }
	      sprintf(&execcommand[i],"%s",infilename);
	      i = strlen(execcommand);
	    }
	  else if(open_command_str[i2] == '%')
	    {
	      i2++;
	      if(i >= (size_execcommand)) {
		size_execcommand *= 2;
		if((execcommand = (char *) realloc(execcommand, (size_execcommand+1))) == NULL)
		  error(ERR_MEMALLOC);
	      }
	      execcommand[i] = '%';
	      i++;
	      execcommand[i] = '\0';
	    }
	  else
	    error(ERR_INVALIDEXECCOMMANDSTRFORMAT);
	}
    }
    
  /* Open the pipe to the light curve command, check for errors, and
     return the handle to the pipe */

  if((return_pipe = popen(execcommand,"r")) == NULL) {
    error2(ERR_CANNOTOPEN,execcommand);
  }
  return return_pipe;
}

#ifdef USECFITSIO

/* Read the data from the match file, assuming it is a FITS binary table. This is similar to the ReadFitsLightCurve function in parselc.c */
_MatchData *ReadFitsMatchFile(_MatchCommand *m, ProgramData *p, int lc, int threadid, char *infilename) {
  _MatchData *matchdata;
  fitsfile *infile;
  int status, hdunum, ncols, hdutype, anynulallcolumns;
  long nrows;
  int j, jold, k, i, l, N, oldsizesinglelc, anynul;
  char *nullarray, *nullarraystore;
  char tmpchar;

  int colnum, datatype;

  double **dblptr;
  short **shortptr;
  int **intptr;
  char **charptr;
  char ***stringptr;
  float **floatptr;
  long **longptr;

  if((matchdata = (_MatchData *) malloc((m->Naddvars + 1)*sizeof(_MatchData))) == NULL)
    error(ERR_MEMALLOC);

#ifdef PARALLEL
  if(p->Nproc_allow > 1) {
    while(pthread_mutex_trylock(&(p->cfitsio_mutex)));
  }
#endif

  hdutype = 0; status = 0; nrows = 0; ncols = 0;

  if((fits_open_file(&infile,infilename,READONLY,&status)))
    {
      error2(ERR_CANNOTOPEN,infilename);
    }
  if(fits_get_hdu_num(infile, &hdunum) == 1)
    {
      fits_movabs_hdu(infile, 2, &hdutype, &status);
    }
  else
    fits_get_hdu_type(infile, &hdutype, &status);
  
  if(hdutype == IMAGE_HDU) {
    error2(ERR_IMAGEHDU,infilename);
  }

  fits_get_num_rows(infile, &nrows, &status);
  fits_get_num_cols(infile, &ncols, &status);

  if((nullarray = (char *) calloc(nrows, sizeof(char))) == NULL)
    error(ERR_MEMALLOC);
  if((nullarraystore = (char *) calloc(nrows, sizeof(char))) == NULL)
    error(ERR_MEMALLOC);

  anynulallcolumns = 0;
  
  for(i=0; i < m->Naddvars+1; i++) {
    matchdata[i].sizevec = 0;
    matchdata[i].Npoints = 0;
    if(!i) {
      if(m->matchvar != NULL) {
	matchdata[i].datatype = m->matchvar->datatype;
      } else {
	if(p->matchstringid) {
	  matchdata[i].datatype = VARTOOLS_TYPE_STRING;
	} else {
	  matchdata[i].datatype = VARTOOLS_TYPE_DOUBLE;
	}
      }
      matchdata[i].format = NULL;
      if(!m->getmatchcolumnfromheader)
	matchdata[i].incol = m->matchcolumn;
      else {
	fits_get_colnum(infile, 0, m->matchcolumn_header_name, &(matchdata[i].incol), &status);
	if(status == COL_NOT_FOUND)
	  error2(ERR_MISSING_MATCHFILE_HEADERNAME,
		 m->matchcolumn_header_name);
      }
    } else {
      matchdata[i].datatype = m->addvar_datatypes[i-1];
      if(m->addvar_incolumn_header_names[i-1] == NULL) {
	matchdata[i].incol = m->addvar_columns[i-1];
      } else {
	fits_get_colnum(infile, 0, m->addvar_incolumn_header_names[i-1], &(matchdata[i].incol), &status);
	if(status == COL_NOT_FOUND)
	  error2(ERR_MISSING_MATCHFILE_HEADERNAME,
		 m->addvar_incolumn_header_names[i-1]);
      }
      if(m->addvar_formats[i-1][0] != '\0')
	matchdata[i].format = m->addvar_formats[i-1];
      else
	matchdata[i].format = NULL;
    }
    matchdata[i].dataptr = NULL;
  }

  MemAllocMatchData(matchdata, m->Naddvars+1, nrows);
  
  /* Read in the columns */
  for(j = 0; j < m->Naddvars+1; j++) {
    status = 0;
    anynul = 0;
    switch(matchdata[j].datatype) {
    case VARTOOLS_TYPE_DOUBLE:
    case VARTOOLS_TYPE_CONVERTJD:
      dblptr = (double **) matchdata[j].dataptr;
      fits_read_colnull(infile, TDOUBLE, matchdata[j].incol, 1, 1, nrows, &((*dblptr)[0]), nullarray, &anynul, &status);
      break;
    case VARTOOLS_TYPE_STRING:
      stringptr = (char ***) matchdata[j].dataptr;
      for(k=1; k<= nrows && !status; k++) {
	fits_read_col_str(infile, matchdata[j].incol, k, 1, 1, 0, &(((*stringptr)[k-1])), &anynul, &status);
	if(anynul)
	  {
	    nullarray[k-1] = 1;
	    anynulallcolumns = 1;
	    anynul = 0;
	  }
      }
      break;
    case VARTOOLS_TYPE_INT:
      intptr = (int **) matchdata[j].dataptr;
      fits_read_colnull(infile, TINT, matchdata[j].incol, 1, 1, nrows, &((*intptr)[0]), nullarray, &anynul, &status);
      break;
    case VARTOOLS_TYPE_FLOAT:
      floatptr = (float **) matchdata[j].dataptr;
      fits_read_colnull(infile, TFLOAT, matchdata[j].incol, 1, 1, nrows, &((*floatptr)[0]), nullarray, &anynul, &status);
      break;
    case VARTOOLS_TYPE_LONG:
      longptr = (long **) matchdata[j].dataptr;
      fits_read_colnull(infile, TLONG, matchdata[j].incol, 1, 1, nrows, &((*longptr)[0]), nullarray, &anynul, &status);
      break;
    case VARTOOLS_TYPE_SHORT:
      shortptr = (short **) matchdata[j].dataptr;
      fits_read_colnull(infile, TSHORT, matchdata[j].incol, 1, 1, nrows, &((*shortptr)[0]), nullarray, &anynul, &status);
      break;
    case VARTOOLS_TYPE_CHAR:
      charptr = (char **) matchdata[j].dataptr;
      fits_read_colnull(infile, TBYTE, matchdata[j].incol, 1, 1, nrows, &((*charptr)[0]), nullarray, &anynul, &status);
      break;
    default:
      error(ERR_BADTYPE);
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
    matchdata[j].Npoints = nrows;
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
	    for(k=0; k < m->Naddvars+1; k++) {
	      switch(matchdata[k].datatype) {
	      case VARTOOLS_TYPE_DOUBLE:
	      case VARTOOLS_TYPE_CONVERTJD:
		dblptr = (double **) matchdata[k].dataptr;
		(*dblptr)[j] = (*dblptr)[i];
		break;
	      case VARTOOLS_TYPE_STRING:
		stringptr = (char ***) matchdata[k].dataptr;
		sprintf((*stringptr)[j],"%s",(*stringptr)[i]);
		break;
	      case VARTOOLS_TYPE_INT:
		intptr = (int **) matchdata[k].dataptr;
		(*intptr)[j] = (*intptr)[i];
		break;
	      case VARTOOLS_TYPE_FLOAT:
		floatptr = (float **) matchdata[k].dataptr;
		(*floatptr)[j] = (*floatptr)[i];
		break;
	      case VARTOOLS_TYPE_LONG:
		longptr = (long **) matchdata[k].dataptr;
		(*longptr)[j] = (*longptr)[i];
		break;
	      case VARTOOLS_TYPE_SHORT:
		shortptr = (short **) matchdata[k].dataptr;
		(*shortptr)[j] = (*shortptr)[i];
		break;
	      case VARTOOLS_TYPE_CHAR:
		charptr = (char **) matchdata[k].dataptr;
		(*charptr)[j] = (*charptr)[i];
		break;
	      default:
		error(ERR_BADTYPE);
	      }
	    }
	  }
	  j++;
	}
      }
    for(k=0; k < m->Naddvars+1; k++) {
      matchdata[k].Npoints = j;
    }
  }

  free(nullarray);
  free(nullarraystore);
  return matchdata;
}

#endif

/* Read the data from the match file, assuming it is an ascii file. This is similar to the ReadSingleLightCurve function in parselc.c */
_MatchData *ReadAsciiMatchFile(_MatchCommand *m, ProgramData *p, int lc, int threadid, char *infilename) {
  _MatchData *matchdata;
  char *line;
  FILE *infile;
  char **incols;
  int u, j, jold, k, i, l, N, oldsizesinglelc, colmax, testskip, Nc;

  double **dblptr;
  short **shortptr;
  int **intptr;
  char **charptr;
  char ***stringptr;
  float **floatptr;
  long **longptr;

  int datatype;

  size_t line_size = MAXLEN;

  if((matchdata = (_MatchData *) malloc((m->Naddvars + 1)*sizeof(_MatchData))) == NULL)
    error(ERR_MEMALLOC);

  for(i=0; i < m->Naddvars+1; i++) {
    matchdata[i].sizevec = 0;
    matchdata[i].Npoints = 0;
    if(!i) {
      if(m->matchvar != NULL) {
	matchdata[i].datatype = m->matchvar->datatype;
      } else {
	if(p->matchstringid) {
	  matchdata[i].datatype = VARTOOLS_TYPE_STRING;
	} else {
	  matchdata[i].datatype = VARTOOLS_TYPE_DOUBLE;
	}
      }
      matchdata[i].format = NULL;
      matchdata[i].incol = m->matchcolumn - 1; 
    } else {
      matchdata[i].datatype = m->addvar_datatypes[i-1];
      matchdata[i].incol = m->addvar_columns[i-1] - 1;
      if(m->addvar_formats[i-1][0] != '\0')
	matchdata[i].format = m->addvar_formats[i-1];
      else
	matchdata[i].format = NULL;
    }
    matchdata[i].dataptr = NULL;
  }

  line = malloc(line_size);
  if(m->opencommand != NULL) {
    infile = ExecMatchOpenCommand(m->opencommand, infilename);
    if(infile == NULL) error2(ERR_FILENOTFOUND, infilename);
  } else {
    if((infile = fopen(infilename,"r")) == NULL) {
      error2(ERR_FILENOTFOUND,infilename);
    }
  }

  colmax = m->matchcolumn;
  for(i=0; i < m->Naddvars; i++) {
    if(m->addvar_columns[i] > colmax)
      colmax = m->addvar_columns[i];
  }
  
  if(colmax > 0) {
    if((incols = (char **) malloc(colmax * sizeof(char *))) == NULL)
      error(ERR_MEMALLOC);
    for(i=0; i < colmax; i++) {
      if((incols[i] = (char *) malloc(MAXLEN)) == NULL)
	error(ERR_MEMALLOC);
    }
  }

  j = 0;
  l = m->Nskip;
  N = 0;
  while(j < m->Nskip ? gnu_getline(&line,&line_size,infile) >= 0 : 0)
    j++;
  while(gnu_getline(&line,&line_size,infile) >= 0)
    {
      l++;
      if(m->delimtype == VARTOOLS_LC_DELIMTYPE_WHITESPACE) {
	testskip = ParseLineToColumns_testskip(line, incols, colmax, m->Nskipchar, m->skipchars);
      } else if(m->delimtype == VARTOOLS_LC_DELIMTYPE_CHAR) {
	testskip = ParseLineToColumnsDelimChar_testskip(line, incols, colmax, m->Nskipchar, m->skipchars, m->delimchar);
      } else if(m->delimtype == VARTOOLS_LC_DELIMTYPE_STRING) {
	testskip = ParseLineToColumnsDelimString_testskip(line, incols, colmax, m->Nskipchar, m->skipchars, m->delimstring);
      }
      if(testskip)
	continue;
      MemAllocMatchData(matchdata, m->Naddvars+1, (N+1));

      for(i=0; i < m->Naddvars+1; i++) {
	k = matchdata[i].incol;
	switch(matchdata[i].datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	case VARTOOLS_TYPE_CONVERTJD:
	  dblptr = (double **) matchdata[i].dataptr;
	  if(matchdata[i].format == NULL)
	    (*dblptr)[N] = atof(incols[k]);
	  else
	    sscanf(incols[k],matchdata[i].format,&((*dblptr)[N]));
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = (char ***) matchdata[i].dataptr;
	  if(matchdata[i].format == NULL)
	    sprintf(((*stringptr)[N]),"%s",incols[k]);
	  else
	    sscanf(incols[k],matchdata[i].format,((*stringptr)[N]));
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = (int **) matchdata[i].dataptr;
	  if(matchdata[i].format == NULL)
	    (*intptr)[N] = atoi(incols[k]);
	  else
	    sscanf(incols[k],matchdata[i].format,&((*intptr)[N]));
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = (short **) matchdata[i].dataptr;
	  if(matchdata[i].format == NULL)
	    (*shortptr)[N]= atoi(incols[k]);
	  else
	    sscanf(incols[k],matchdata[i].format,&((*shortptr)[N]));
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = (float **) matchdata[i].dataptr;
	  if(matchdata[i].format == NULL)
	    (*floatptr)[N] = atof(incols[k]);
	  else
	    sscanf(incols[k],matchdata[i].format,&((*floatptr)[N]));
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = (long **) matchdata[i].dataptr;
	  if(matchdata[i].format == NULL)
	    (*longptr)[N] = atol(incols[k]);
	  else
	    sscanf(incols[k],matchdata[i].format,&((*longptr)[N]));
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = (char **) matchdata[i].dataptr;
	  if(matchdata[i].format == NULL)
	    (*charptr)[N] = incols[k][0];
	  else
	    sscanf(incols[k],matchdata[i].format,&((*charptr)[N]));
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      }
      N++;
    }

  for(i=0; i < m->Naddvars+1; i++) {
    matchdata[i].Npoints = N;
  }

  if(colmax > 0) {
    for(i=0; i < colmax; i++)
      free(incols[i]);
    free(incols);
  }
  if(m->opencommand != NULL)
    pclose(infile);
  else
    fclose(infile);

  free(line);

  return(matchdata);
}
      

/* Read the data from the match file */
_MatchData *ReadMatchFile(_MatchCommand *m, ProgramData *p, int lc, int threadid)
{
  char *infilename;
  int j, i;


  if(m->inputfilename != NULL) {
    infilename = m->inputfilename;
  } else {
    infilename = m->inputfilenamelist[lc];
  }
  
  /* Check if the end of the file is .fits, if so assume it is a binary fits
     table. */

#ifdef USECFITSIO

  j = strlen(infilename);
  i = j - 5;
  if(i >= 0) {
    if(!strcmp(&(infilename[i]),".fits")) {
      return(ReadFitsMatchFile(m, p, lc, threadid, infilename));
    }
  }

#endif

  /* If not, then assume it is in ascii format */
  return(ReadAsciiMatchFile(m, p, lc, threadid, infilename));
}


/* This is similar to the CheckCreateCommandOutputLCVariable function
   in analytic.c, but here we need to allow the variables to have
   other types besides DOUBLE */
void SetupMatchCommandVariables(_MatchCommand *m, ProgramData *p)
{
  int j, k;
  _Variable *v;
  _Variable **omodelvar;
  char *varname;

  if(m->matchcolumnvarname != NULL) {
    for(j=0; j < p->NDefinedVariables; j++) {
      v = p->DefinedVariables[j];
      if(!strcmp(m->matchcolumnvarname, v->varname)) {
	if(v->vectortype != VARTOOLS_VECTORTYPE_LC) {
	  error2(ERR_INVALIDVARIABLELCVARIABLE, m->matchcolumnvarname);
	}
	m->matchvar = v;
	break;
      }
    }
    if(j == p->NDefinedVariables) {
      /* We cannot match on a variable that has not been set already for the light curve */
      error2(ERR_MATCHCOMMAND_BADMATCHVARIABLE,m->matchcolumnvarname);
    }
  }

  for(k=0; k < m->Naddvars; k++) {
    varname = m->addvar_varnames[k];
    if(varname != NULL) {
      if(varname[0] != '\0') {
	for(j=0; j < p->NDefinedVariables; j++) {
	  v = p->DefinedVariables[j];
	  if(!strcmp(varname,v->varname)) {
	    /* This is an existing variable, make sure it is the correct type */
	    if(v->vectortype != VARTOOLS_VECTORTYPE_LC ||
	       v->datatype != m->addvar_datatypes[k]) {
	      error2(ERR_INVALIDVARIABLELCVARIABLE,varname);
	    }
	    m->addvars[k] = v;
	    break;
	  }
	}
	if(j == p->NDefinedVariables) {
	  /* This is a new variable, create it */
	  v = CreateVariable(p, varname, m->addvar_datatypes[k],
			     VARTOOLS_VECTORTYPE_LC, NULL);
	  RegisterDataFromLightCurve(p,
				     v->dataptr,
				     m->addvar_datatypes[k],
				     (m->addvar_datatypes[k] == 
				      VARTOOLS_TYPE_STRING ? MAXLEN : 0), 
				     0, 0, 0, 0, NULL,
				     v,
				     -1, varname);
	  m->addvars[k] = v;
	}
      }
    }
  }
}

/* This function parses a string of the form:
   var1:col1[:type1[:fmt][,var2:col2[:type2[:fmt2]],...]] for the
   -match command */
int ParseMatchInputFormatString(char *argv, _MatchCommand *m)
{
  int i, j, termscan, startindex;
  char *varname;
  int column;
  int datatype;
  char *format;
  char *parsecopy;
  int sizestring;
  char *incolumn_header_name = NULL;

  sizestring = MAXLEN;
  format = (char *) malloc(sizestring);
  parsecopy = (char *) malloc(sizestring);
  varname = (char *) malloc(sizestring);

  varname[0] = '\0'; column = -1; datatype = -1; format[0] = '\0';
  termscan = 0;

  m->getcolumnsfromheader = 0;

  m->Naddvars = 0;

  i = 0, j=0;
  startindex = 0;
  do {
    if(argv[i] == ',' || (argv[i] == ':' && termscan != 3) || argv[i] == '\0') {
      if(startindex == i) {
	free(parsecopy);
	free(format);
	free(varname);
	return 1;
      }
      if(j >= sizestring) {
	sizestring *= 2;
	if((parsecopy = (char *) realloc(parsecopy, sizestring)) == NULL ||
	   (format = (char *) realloc(format, sizestring)) == NULL ||
	   (varname = (char *) realloc(varname, sizestring)) == NULL)
	  error(ERR_MEMALLOC);
      }
      parsecopy[j] = '\0';
      if(termscan == 0) {
	sprintf(varname,"%s",parsecopy);
      }
      else if(termscan == 1) {
#ifdef _USEBINARY_LC
	if(!isstringint(parsecopy)) {
	  incolumn_header_name = malloc(strlen(parsecopy)+1);
	  sprintf(incolumn_header_name,"%s",parsecopy);
	  column = 1;
	  m->getcolumnsfromheader = 1;
	}
	else {
	  column = atoi(parsecopy);
	}
#elif USECFITSIO
	if(!isstringint(parsecopy)) {
	  incolumn_header_name = malloc(strlen(parsecopy)+1);
	  sprintf(incolumn_header_name,"%s",parsecopy);
	  column = 1;
	  m->getcolumnsfromheader = 1;
	}
	else {
	  column = atoi(parsecopy);
	}
#else
	column = atoi(parsecopy);
#endif
      } 
      else if(termscan == 2) {
	if(!strcmp(parsecopy,"double")) {
	  datatype = VARTOOLS_TYPE_DOUBLE;
	}
	else if(!strcmp(parsecopy,"float")) {
	  datatype = VARTOOLS_TYPE_FLOAT;
	}
	else if(!strcmp(parsecopy,"int")) {
	  datatype = VARTOOLS_TYPE_INT;
	}
	else if(!strcmp(parsecopy,"long")) {
	  datatype = VARTOOLS_TYPE_LONG;
	}
	else if(!strcmp(parsecopy,"short")) {
	  datatype = VARTOOLS_TYPE_SHORT;
	}
	else if(!strcmp(parsecopy,"char")) {
	  datatype = VARTOOLS_TYPE_CHAR;
	}
	else if(!strcmp(parsecopy,"utc")) {
	  datatype = VARTOOLS_TYPE_CONVERTJD;
	}
	else if(!strcmp(parsecopy,"string")) {
	  datatype = VARTOOLS_TYPE_STRING;
	}
	else {
	  free(parsecopy);
	  free(format);
	  free(varname);
	  return 1;
	}
      }
      else if(termscan == 3) {
	sprintf(format, "%s", parsecopy);
      }
      else {
	free(parsecopy);
	free(format);
	free(varname);
	return 1;
      }
      j = 0;
      startindex = i+1;
      if(argv[i] == ':') {
	termscan++;
      }
      else {
	if(termscan < 1 || column < 0) {
	  free(parsecopy);
	  free(format);
	  free(varname);
	  return 1;
	}
	if(datatype < 0) {
	  if(!strcmp(varname,"id")) datatype = VARTOOLS_TYPE_STRING;
	  else datatype = VARTOOLS_TYPE_DOUBLE;
	}
	if(!m->Naddvars) {
	  if((m->addvars = (_Variable **) malloc(sizeof(_Variable *))) == NULL ||
	     (m->addvar_varnames = (char **) malloc(sizeof(char *))) == NULL ||
	     (m->addvar_columns = (int *) malloc(sizeof(int))) == NULL ||
	     (m->addvar_datatypes = (int *) malloc(sizeof(int))) == NULL ||
	     (m->addvar_formats = (char **) malloc(sizeof(char *))) == NULL ||
	     (m->addvar_incolumn_header_names = (char **) malloc(sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	} else {
	  if((m->addvars = (_Variable **) realloc(m->addvars, (m->Naddvars+1)*sizeof(_Variable *))) == NULL ||
	     (m->addvar_varnames = (char **) realloc(m->addvar_varnames, (m->Naddvars+1)*sizeof(char *))) == NULL ||
	     (m->addvar_columns = (int *) realloc(m->addvar_columns, (m->Naddvars+1)*sizeof(int))) == NULL ||
	     (m->addvar_datatypes = (int *) realloc(m->addvar_datatypes, (m->Naddvars+1)*sizeof(int))) == NULL ||
	     (m->addvar_formats = (char **) realloc(m->addvar_formats, (m->Naddvars+1)*sizeof(char *))) == NULL ||
	     (m->addvar_incolumn_header_names = (char **) realloc(m->addvar_incolumn_header_names, (m->Naddvars+1)*sizeof(char *))) == NULL) 
	    error(ERR_MEMALLOC);
	}
	if((m->addvar_varnames[m->Naddvars] = (char *) malloc((strlen(varname)+1)*sizeof(char))) == NULL ||
	   (m->addvar_formats[m->Naddvars] = (char *) malloc((strlen(format)+1)*sizeof(char))) == NULL)
	  error(ERR_MEMALLOC);
	sprintf(m->addvar_varnames[m->Naddvars],"%s",varname);
	m->addvar_columns[m->Naddvars] = column;
	m->addvar_datatypes[m->Naddvars] = datatype;
	sprintf(m->addvar_formats[m->Naddvars],"%s",format);
	if(incolumn_header_name != NULL) {
	   if((m->addvar_incolumn_header_names[m->Naddvars] = (char *) malloc((strlen(incolumn_header_name)+1)*sizeof(char))) == NULL)
	     error(ERR_MEMALLOC);
	   sprintf(m->addvar_incolumn_header_names[m->Naddvars],"%s",incolumn_header_name);
	} else
	  m->addvar_incolumn_header_names[m->Naddvars] = NULL;
	m->Naddvars += 1;
	varname[0] = '\0';
	column = -1;
	datatype = -1;
	format[0] = '\0';
	termscan = 0;
	if(incolumn_header_name != NULL) {
	  free(incolumn_header_name);
	  incolumn_header_name = NULL;
	}
	if(argv[i] == '\0')
	  break;
      }
    }
    else {
      parsecopy[j] = argv[i];
      j++;
    }
    i++;
  } while(1);

  free(format);
  free(parsecopy);
  free(varname);
  return 0;
}
  

/* -match [\"file\" filename | \"inlist\" inlistcolumn] [\"opencommand\" command] [\"skipnum\" Nskip] [\"skipchar\" <skipchar1[,skipchar2,...]>] [\"delimiter\" delimiter] <\"matchcolumn\" <varname:colnum | colnum>> <\"addcolumns\" varname1:colnum1[:coltype1[:colformat1]][,varname2:colnum2[:coltype2[:colformat2]],...]> <\"cullmissing\" | \"nanmissing\" | \"missingval\" value> */
int ParseMatchCommand(int *iret, int argc, char **argv, ProgramData *p, _MatchCommand *m, int cn)
{
  int i;
  int j;
  m->inputfilename = NULL;
  m->inlistcolumn = -1;
  m->Ninlist = 0;
  m->inputfilenamelist = NULL;
  m->opencommand = NULL;

  m->Nskip = 0;
  
  m->delimtype = VARTOOLS_LC_DELIMTYPE_WHITESPACE; 

  m->matchvar = NULL;

  m->matchcolumn = -1;
  m->matchcolumnvarname = NULL;
  
  m->missingmethod = -1;

  m->missingvalue = 0.0;

  m->getmatchcolumnfromheader = 0;

  i = *iret;
  if(i >= argc)
    return 1;
  if(!strcmp(argv[i],"file")) {
    i++;
    if(i >= argc)
      return 1;
    if((m->inputfilename = (char *) malloc(strlen(argv[i])+1)) == NULL)
      error(ERR_MEMALLOC);
    sprintf(m->inputfilename, "%s", argv[i]);
  }
  else if(!strcmp(argv[i],"inlist")) {
    i++;
    if(i >= argc)
      return 1;
    m->inlistcolumn = atoi(argv[i]);
    if(m->inlistcolumn <= 0) {
      error(ERR_MATCHCOMMAND_INLISTCOLUMNNOTPOSITIVE);
    }
    RegisterDataFromInputList(p,
			      (void *) (&(m->inputfilenamelist)),
			      VARTOOLS_TYPE_STRING,
			      0, cn, 0, 0, NULL, m->inlistcolumn,
			      "MATCH_FILENAME");
  }
  else
    return 1;

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"opencommand")) {
      i++;
      if(i >= argc)
	return 1;
      if((m->opencommand = (char *) malloc(strlen(argv[i])+1)) == NULL)
	error(ERR_MEMALLOC);
      sprintf(m->opencommand,"%s",argv[i]);
    }
    else
      i--;
  } else
    i--;

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"skipnum")) {
      i++;
      if(i >= argc)
	return 1;
      m->Nskip = atoi(argv[i]);
    } else
      i--;
  } else
    i--;

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"skipchar")) {
      i++;
      if(i >= argc)
	return 1;
      ParseSkipCharString(argv[i],&m->Nskipchar,&m->skipchars);
    } else {
      m->Nskipchar = 1;
      m->skipchars = malloc(1);
      m->skipchars[0] = '#';
      i--;
    }
  } else {
    m->Nskipchar = 1;
    m->skipchars = malloc(1);
    m->skipchars[0] = '#';
    i--;
  }

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"delimiter")) {
      i++;
      if(i >= argc)
	return 1;
      if(strlen(argv[i]) == 1) {
	m->delimtype = VARTOOLS_LC_DELIMTYPE_CHAR;
	m->delimchar = argv[i][0];
      } else {
	m->delimtype = VARTOOLS_LC_DELIMTYPE_STRING;
	m->delimstring = malloc(strlen(argv[i])+1);
	sprintf(m->delimstring,"%s",argv[i]);
      }
    } else
      i--;
  } else
    i--;

  i++;
  if(i >= argc)
    return 1;
  if(!strcmp(argv[i],"matchcolumn")) {
    i++;
    if(i >= argc)
      return 1;
    j = 0;
    while(argv[i][j] != '\0') {
      if(argv[i][j] == ':')
	break;
      j++;
    }
    if(argv[i][j] == ':') {
      if(!j)
	return 1;
      if(argv[i][j+1] == '\0')
	return 1;
      if((m->matchcolumnvarname = (char *) malloc(j+2)) == NULL)
	error(ERR_MEMALLOC);
      memcpy(m->matchcolumnvarname, argv[i], (size_t) (j*sizeof(char)));
      m->matchcolumnvarname[j] = '\0';
      if(!isstringint(&(argv[i][j+1]))) {
	m->matchcolumn_header_name = malloc(strlen(&(argv[i][j+1]))+1);
	sprintf(m->matchcolumn_header_name,"%s",&(argv[i][j+1]));
	m->matchcolumn = 1;
	m->getmatchcolumnfromheader = 1;
      } else {
	m->matchcolumn = atoi(&(argv[i][j+1]));
	if(m->matchcolumn <= 0) 
	  error(ERR_MATCHCOMMAND_MATCHCOLUMNNOTPOSITIVE);
      }
    } else {
      if(!j)
	return 1;
      if(!isstringint(argv[i])) {
	m->matchcolumn_header_name = malloc(strlen(argv[i])+1);
	sprintf(m->matchcolumn_header_name,"%s",argv[i]);
	m->matchcolumn = 1;
	m->getmatchcolumnfromheader = 1;
      } else {
	m->matchcolumn = atoi(argv[i]);
	if(m->matchcolumn <= 0)
	  error(ERR_MATCHCOMMAND_MATCHCOLUMNNOTPOSITIVE);
      }
    }
  } else 
    return 1;

  i++;
  if(i >= argc)
    return 1;
  if(!strcmp(argv[i],"addcolumns")) {
    i++;
    if(i >= argc)
      return 1;
    if(ParseMatchInputFormatString(argv[i], m))
      return 1;
  } else
    return 1;

  i++;
  if(i >= argc)
    return 1;
  if(!strcmp(argv[i],"cullmissing")) {
    m->missingmethod = VARTOOLS_MATCHCOMMAND_MISSINGMETHOD_CULLMISSING;
  }
  else if(!strcmp(argv[i],"nanmissing")) {
    m->missingmethod = VARTOOLS_MATCHCOMMAND_MISSINGMETHOD_NANMISSING;
  }
  else if(!strcmp(argv[i],"missingval")) {
    m->missingmethod = VARTOOLS_MATCHCOMMAND_MISSINGMETHOD_MISSINGVAL;
    i++;
    if(i >= argc)
      return 1;
    m->missingvalue = atof(argv[i]);
  }
  else
    return 1;

  *iret = i;
  return 0;
}

/* Fill out an index which sorts the column to match on; 
   if *indxptr is NULL on input, the memory will be allocated for
   it, otherwise it will be assumed to be the correct size already */
void SortMatchFileColumn(_MatchData *md, int **indxptr) {
  int *indx;
  int i;
  void *dataptr;
  
  if(md->Npoints <= 0)
    return;

  if((*indxptr) == NULL) {
    if((indx = (int *) malloc(md->Npoints * sizeof(int))) == NULL)
      error(ERR_MEMALLOC);
    *indxptr = indx;
  } else {
    indx = *indxptr;
  }
  for(i = 0; i < md->Npoints; i++)
    indx[i] = i;
  switch(md->datatype) {
  case VARTOOLS_TYPE_DOUBLE:
  case VARTOOLS_TYPE_CONVERTJD:
    dataptr = (void *) (((double **) md->dataptr)[0]);
    break;
  case VARTOOLS_TYPE_STRING:
    dataptr = (void *) (((char ***) md->dataptr)[0]);
    break;
  case VARTOOLS_TYPE_INT:
    dataptr = (void *) (((int **) md->dataptr)[0]);
    break;
  case VARTOOLS_TYPE_FLOAT:
    dataptr = (void *) (((float **) md->dataptr)[0]);
    break;
  case VARTOOLS_TYPE_LONG:
    dataptr = (void *) (((long **) md->dataptr)[0]);
    break;
  case VARTOOLS_TYPE_CHAR:
    dataptr = (void *) (((char **) md->dataptr)[0]);
    break;
  case VARTOOLS_TYPE_SHORT:
    dataptr = (void *) (((short **) md->dataptr)[0]);
    break;
  default:
    error(ERR_BADTYPE);
  }

  if(md->datatype == VARTOOLS_TYPE_STRING)
    sort_generic(md->Npoints, 0, indx, 1, md->datatype, dataptr, 1, MAXLEN);
  else
    sort_generic(md->Npoints, 0, indx, 1, md->datatype, dataptr, 1);
}

/* Read in the match file if a global match file is being used */
void InitMatchCommand(ProgramData *p, _MatchCommand *m, int lc, int threadid)
{
  
  if(m->inputfilename == NULL)
    return;

  /* Read any global match file, and sort it on the column for matching */
  m->matchdata = ReadMatchFile(m, p, lc, threadid);
  m->matchdata_sortindx = NULL;
  SortMatchFileColumn(m->matchdata, &(m->matchdata_sortindx));
}

/* Run the match command */
void RunMatchCommand(ProgramData *p, _MatchCommand *m, int lc, int threadid)
{
  _MatchData *md;
  int *md_indx = NULL;
  int *lc_indx = NULL;
  int lc_indx_isalloc = 0;
  int md_indx_isalloc = 0;
  int *match_indx = NULL;
  int i, j, k, kk, isclip = 0;
  void *dataptr;
  int lcmatch_datatype;


  /* Nothing to do, just return */
  if(!p->NJD[threadid])
    return;

  if(m->inputfilename == NULL) {
    /* We are not using a global match file, so read it in for this
       light curve, and sort on the match column */
    md = ReadMatchFile(m, p, lc, threadid);
    SortMatchFileColumn(md, &(md_indx));
    md_indx_isalloc = 1;
  } else {
    /* We are using the global match file, which should already have been read in */
    md = m->matchdata;
    md_indx = m->matchdata_sortindx;
  }

    if((match_indx = (int *) malloc(p->NJD[threadid]*sizeof(int))) == NULL)
      error(ERR_MEMALLOC);

    for(i=0; i < p->NJD[threadid]; i++) {
      match_indx[i] = -1;
    }


  /* Prepare the index for matching; 
     The process for matching will depend on whether we are matching
     by time, by ID, or by another variable */

  if(m->matchvar != NULL) {
    /* We are matching on a variable specified by the USER */

    if((lc_indx = (int *) malloc(p->NJD[threadid]*sizeof(int))) == NULL)
      error(ERR_MEMALLOC);
    
    for(i=0; i < p->NJD[threadid]; i++) {
      lc_indx[i] = i;
    }

    lc_indx_isalloc = 1;

    lcmatch_datatype = m->matchvar->datatype;

    /* Get a pointer to the variable data */
    switch(m->matchvar->datatype) {
    case VARTOOLS_TYPE_DOUBLE:
    case VARTOOLS_TYPE_CONVERTJD:
      dataptr = (void *) ((*((double ***) m->matchvar->dataptr))[threadid]);
      break;
    case VARTOOLS_TYPE_FLOAT:
      dataptr = (void *) ((*((float ***) m->matchvar->dataptr))[threadid]);
      break;
    case VARTOOLS_TYPE_INT:
      dataptr = (void *) ((*((int ***) m->matchvar->dataptr))[threadid]);
      break;
    case VARTOOLS_TYPE_SHORT:
      dataptr = (void *) ((*((short ***) m->matchvar->dataptr))[threadid]);
      break;
    case VARTOOLS_TYPE_LONG:
      dataptr = (void *) ((*((long ***) m->matchvar->dataptr))[threadid]);
      break;
    case VARTOOLS_TYPE_CHAR:
      dataptr = (void *) ((*((char ***) m->matchvar->dataptr))[threadid]);
      break;
    case VARTOOLS_TYPE_STRING:
      dataptr = (void *) ((*((char ****) m->matchvar->dataptr))[threadid]);
      break;
    default:
      error(ERR_BADTYPE);
      break;
    }

    /* Sort the light curve vector into the index */
    if(m->matchvar->datatype == VARTOOLS_TYPE_STRING)
      sort_generic(p->NJD[threadid], 0, lc_indx, 1, m->matchvar->datatype, dataptr, 1, MAXLEN);
    else
      sort_generic(p->NJD[threadid], 0, lc_indx, 1, m->matchvar->datatype, dataptr, 1);
    
  } 
  else if(p->matchstringid)  {
    /* We are matching on the string image ID; 
       A sorted index is already available, we do not need to 
       do this again */
    dataptr = (void *) (p->stringid[threadid]);
    lc_indx = p->stringid_idx[threadid];

    lcmatch_datatype = VARTOOLS_TYPE_STRING;

  } else {

    /* We are matching on the time variable;
       The light curve is already sorted by time, we just need to fill out
       the index */
    dataptr = (void *) (p->t[threadid]);

    if((lc_indx = (int *) malloc(p->NJD[threadid]*sizeof(int))) == NULL)
      error(ERR_MEMALLOC);
    
    for(i=0; i < p->NJD[threadid]; i++) {
      lc_indx[i] = i;
    }

    lc_indx_isalloc = 1;

    lcmatch_datatype = VARTOOLS_TYPE_DOUBLE;

  }

  /* Perform the match, handling each datatype individually */
  if(lcmatch_datatype == VARTOOLS_TYPE_STRING) {
    char **stringvec_md, **stringvec_lc;
    
    stringvec_lc = (char **) dataptr;
    if(md[0].Npoints > 0)
      stringvec_md = *((char ***) md[0].dataptr);

    j = 0;
    k = 0;
    kk = 0;

    while(k < md[0].Npoints && j < p->NJD[threadid]) {
      while((j < p->NJD[threadid] ? strncmp(stringvec_md[md_indx[k]], stringvec_lc[lc_indx[j]],MAXIDSTRINGLENGTH) > 0 : 0)) 
	j++;
      while ( j < p->NJD[threadid] ? !strncmp(stringvec_md[md_indx[k]], stringvec_lc[lc_indx[j]], MAXIDSTRINGLENGTH) : 0) {
	match_indx[lc_indx[j]] = md_indx[k];
	j++;
      }
      kk = k;
      k++;
      while(k < md[0].Npoints ? !strncmp(stringvec_md[md_indx[kk]],stringvec_md[md_indx[k]], MAXIDSTRINGLENGTH) : 0)
	k++;
    }
  }
  else if(lcmatch_datatype == VARTOOLS_TYPE_DOUBLE || 
	  lcmatch_datatype == VARTOOLS_TYPE_CONVERTJD) {
    double *dblvec_md, *dblvec_lc;
    
    dblvec_lc = (double *) dataptr;
    if(md[0].Npoints > 0)
      dblvec_md = *((double **) md[0].dataptr);

    j = 0;
    k = 0;

    while(k < md[0].Npoints && j < p->NJD[threadid]) {
      while((j < p->NJD[threadid] ? dblvec_lc[lc_indx[j]] < dblvec_md[md_indx[k]] - JDTOL : 0))
	j++;
      while ( j < p->NJD[threadid] ? dblvec_lc[lc_indx[j]] >= dblvec_md[md_indx[k]] - JDTOL && dblvec_lc[lc_indx[j]] <= dblvec_md[md_indx[k]] + JDTOL : 0) {
	match_indx[lc_indx[j]] = md_indx[k];
	j++;
      }
      kk = k;
      k++;
      while(k < md[0].Npoints ? dblvec_md[md_indx[k]] >= dblvec_md[md_indx[kk]] - JDTOL && dblvec_md[md_indx[k]] < dblvec_md[md_indx[kk]] + JDTOL : 0)
	k++;
    }
  }
  else if(lcmatch_datatype == VARTOOLS_TYPE_FLOAT) {
    float *fltvec_md, *fltvec_lc;
    
    fltvec_lc = (float *) dataptr;
    if(md[0].Npoints > 0)
      fltvec_md = *((float **) md[0].dataptr);

    j = 0;
    k = 0;

    while(k < md[0].Npoints && j < p->NJD[threadid]) {
      while((j < p->NJD[threadid] ? fltvec_lc[lc_indx[j]] < fltvec_md[md_indx[k]] - JDTOL : 0))
	j++;
      while ( j < p->NJD[threadid] ? fltvec_lc[lc_indx[j]] >= fltvec_md[md_indx[k]] - JDTOL && fltvec_lc[lc_indx[j]] <= fltvec_md[md_indx[k]] + JDTOL : 0) {
	match_indx[lc_indx[j]] = md_indx[k];
	j++;
      }
      kk = k;
      k++;
      while(k < md[0].Npoints ? fltvec_md[md_indx[k]] >= fltvec_md[md_indx[kk]] - JDTOL && fltvec_md[md_indx[k]] < fltvec_md[md_indx[kk]] + JDTOL : 0)
	k++;
    }
  }
  else if(lcmatch_datatype == VARTOOLS_TYPE_INT) {
    int *intvec_md, *intvec_lc;
    
    intvec_lc = (int *) dataptr;
    if(md[0].Npoints > 0)
      intvec_md = *((int **) md[0].dataptr);

    j = 0;
    k = 0;

    while(k < md[0].Npoints && j < p->NJD[threadid]) {
      while((j < p->NJD[threadid] ? intvec_lc[lc_indx[j]] < intvec_md[md_indx[k]] : 0))
	j++;
      while ( j < p->NJD[threadid] ? intvec_lc[lc_indx[j]] == intvec_md[md_indx[k]] : 0) {
	match_indx[lc_indx[j]] = md_indx[k];
	j++;
      }
      kk = k;
      k++;
      while(k < md[0].Npoints ? intvec_md[md_indx[k]] == intvec_md[md_indx[kk]] : 0)
	k++;
    }
  }
  else if(lcmatch_datatype == VARTOOLS_TYPE_LONG) {
    long *longvec_md, *longvec_lc;
    
    longvec_lc = (long *) dataptr;
    if(md[0].Npoints > 0)
      longvec_md = *((long **) md[0].dataptr);

    j = 0;
    k = 0;

    while(k < md[0].Npoints && j < p->NJD[threadid]) {
      while((j < p->NJD[threadid] ? longvec_lc[lc_indx[j]] < longvec_md[md_indx[k]] : 0))
	j++;
      while ( j < p->NJD[threadid] ? longvec_lc[lc_indx[j]] == longvec_md[md_indx[k]] : 0) {
	match_indx[lc_indx[j]] = md_indx[k];
	j++;
      }
      kk = k;
      k++;
      while(k < md[0].Npoints ? longvec_md[md_indx[k]] == longvec_md[md_indx[kk]] : 0)
	k++;
    }
  }
  else if(lcmatch_datatype == VARTOOLS_TYPE_SHORT) {
    short *shortvec_md, *shortvec_lc;
    
    shortvec_lc = (short *) dataptr;
    if(md[0].Npoints > 0)
      shortvec_md = *((short **) md[0].dataptr);

    j = 0;
    k = 0;

    while(k < md[0].Npoints && j < p->NJD[threadid]) {
      while((j < p->NJD[threadid] ? shortvec_lc[lc_indx[j]] < shortvec_md[md_indx[k]] : 0))
	j++;
      while ( j < p->NJD[threadid] ? shortvec_lc[lc_indx[j]] == shortvec_md[md_indx[k]] : 0) {
	match_indx[lc_indx[j]] = md_indx[k];
	j++;
      }
      kk = k;
      k++;
      while(k < md[0].Npoints ? shortvec_md[md_indx[k]] == shortvec_md[md_indx[kk]] : 0)
	k++;
    }
  }
  else if(lcmatch_datatype == VARTOOLS_TYPE_CHAR) {
    char *charvec_md, *charvec_lc;
    
    charvec_lc = (char *) dataptr;
    if(md[0].Npoints > 0)
      charvec_md = *((char **) md[0].dataptr);

    j = 0;
    k = 0;

    while(k < md[0].Npoints && j < p->NJD[threadid]) {
      while((j < p->NJD[threadid] ? charvec_lc[lc_indx[j]] < charvec_md[md_indx[k]] : 0))
	j++;
      while ( j < p->NJD[threadid] ? charvec_lc[lc_indx[j]] == charvec_md[md_indx[k]] : 0) {
	match_indx[lc_indx[j]] = md_indx[k];
	j++;
      }
      kk = k;
      k++;
      while(k < md[0].Npoints ? charvec_md[md_indx[k]] == charvec_md[md_indx[kk]] : 0)
	k++;
    }
  }
  else
    error(ERR_MEMALLOC);
    
  /* Now copy over the match data from the match file into the light curve */
  for(j=0; j < p->NJD[threadid]; j++) {
    if(match_indx[j] >= 0) {
      for(i=0; i < m->Naddvars; i++) {
	switch(m->addvars[i]->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	case VARTOOLS_TYPE_CONVERTJD:
	  (*((double ***) m->addvars[i]->dataptr))[threadid][j] = (*((double **) md[i+1].dataptr))[match_indx[j]];
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  (*((float ***) m->addvars[i]->dataptr))[threadid][j] = (*((float **) md[i+1].dataptr))[match_indx[j]];
	  break;
	case VARTOOLS_TYPE_INT:
	  (*((int ***) m->addvars[i]->dataptr))[threadid][j] = (*((int **) md[i+1].dataptr))[match_indx[j]];
	  break;
	case VARTOOLS_TYPE_LONG:
	  (*((long ***) m->addvars[i]->dataptr))[threadid][j] = (*((long **) md[i+1].dataptr))[match_indx[j]];
	  break;
	case VARTOOLS_TYPE_SHORT:
	  (*((short ***) m->addvars[i]->dataptr))[threadid][j] = (*((short **) md[i+1].dataptr))[match_indx[j]];
	  break;
	case VARTOOLS_TYPE_CHAR:
	  (*((char ***) m->addvars[i]->dataptr))[threadid][j] = (*((char **) md[i+1].dataptr))[match_indx[j]];
	  break;
	case VARTOOLS_TYPE_STRING:
	  sprintf((*((char ****) m->addvars[i]->dataptr))[threadid][j],"%s",(*((char ***) md[i+1].dataptr))[match_indx[j]]);
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      }
    } else {
      /* No match was found; What we do depends on the missingmethod */
      if(m->missingmethod == VARTOOLS_MATCHCOMMAND_MISSINGMETHOD_CULLMISSING) {
	/* We will deal with culling the missing points later if needed */
	isclip = 1;
      }
      else if(m->missingmethod == VARTOOLS_MATCHCOMMAND_MISSINGMETHOD_NANMISSING) {
	for(i=0; i < m->Naddvars; i++) {
	  switch(m->addvars[i]->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	  case VARTOOLS_TYPE_CONVERTJD:
	    (*((double ***) m->addvars[i]->dataptr))[threadid][j] = sqrt(-1.0);
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    (*((float ***) m->addvars[i]->dataptr))[threadid][j] = sqrt(-1.0);
	    break;
	  case VARTOOLS_TYPE_INT:
	    (*((int ***) m->addvars[i]->dataptr))[threadid][j] = 0;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    (*((long ***) m->addvars[i]->dataptr))[threadid][j] = 0;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    (*((short ***) m->addvars[i]->dataptr))[threadid][j] = 0;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    (*((char ***) m->addvars[i]->dataptr))[threadid][j] = '0';
	    break;
	  case VARTOOLS_TYPE_STRING:
	    sprintf((*((char ****) m->addvars[i]->dataptr))[threadid][j],"NaN");
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      }
      else if(m->missingmethod == VARTOOLS_MATCHCOMMAND_MISSINGMETHOD_MISSINGVAL) {
	for(i=0; i < m->Naddvars; i++) {
	  switch(m->addvars[i]->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	  case VARTOOLS_TYPE_CONVERTJD:
	    (*((double ***) m->addvars[i]->dataptr))[threadid][j] = m->missingvalue;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    (*((float ***) m->addvars[i]->dataptr))[threadid][j] = (float) m->missingvalue;
	    break;
	  case VARTOOLS_TYPE_INT:
	    (*((int ***) m->addvars[i]->dataptr))[threadid][j] = (int) m->missingvalue;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    (*((long ***) m->addvars[i]->dataptr))[threadid][j] = (long) m->missingvalue;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    (*((short ***) m->addvars[i]->dataptr))[threadid][j] = (short) m->missingvalue;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    (*((char ***) m->addvars[i]->dataptr))[threadid][j] = (char) m->missingvalue;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    sprintf((*((char ****) m->addvars[i]->dataptr))[threadid][j],"%.17g",m->missingvalue);
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      }
    }
  }

  /* Now deal with the case where we are clipping points */
  if(isclip) {
    j = 0;
    for(k = 0; k < p->NJD[threadid]; k++) {
      if(match_indx[k] >= 0) {
	if(j != k)
	  sigclip_copyterms(k, j, p, threadid);
	j++;
      }
    }
    p->NJD[threadid] = j;
    if(p->readimagestring) {
      for(i=0; i < p->NJD[threadid]; i++) 
	p->stringid_idx[threadid][i] = i;
      mysortstringint(p->NJD[threadid], MAXIDSTRINGLENGTH, p->stringid[threadid], p->stringid_idx[threadid]);
    }
  }

  if(lc_indx_isalloc && lc_indx != NULL)
    free(lc_indx);

  if(md_indx_isalloc && md_indx != NULL)
    free(md_indx);

  if(match_indx != NULL)
    free(match_indx);

  if(m->inputfilename == NULL)
    FreeMatchData(m, md);
}
