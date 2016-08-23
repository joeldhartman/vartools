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

/* This function parses a string of the form:
   var1:col1[:type1[:fmt1][,var2:col2[:type2[:fmt2]],...]]
   and creates the relevant _DataFromInputList and _Variable objects */
int ParseInListVarsString(char *argv, ProgramData *p)
{
  int i, j, termscan, startindex;
  char *varname;
  int column;
  int datatype;
  char *format;
  char *parsecopy;
  int sizestring;

  sizestring = MAXLEN;
  format = (char *) malloc(sizestring);
  parsecopy = (char *) malloc(sizestring);
  varname = (char *) malloc(sizestring);

  varname[0] = '\0'; column = -1; datatype = -1; format[0] = '\0';
  termscan = 0;

  i = 0, j=0;
  startindex = 0;
  do {
    if(argv[i] == ',' || argv[i] == ':' || argv[i] == '\0') {
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
	column = atoi(parsecopy);
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
	else if(!strcmp(parsecopy,"string")) {
	  datatype = VARTOOLS_TYPE_STRING;
	}
	else if(!strcmp(parsecopy,"utc")) {
	  datatype = VARTOOLS_TYPE_CONVERTJD;
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
	  datatype = VARTOOLS_TYPE_DOUBLE;
	}
	if(!strcmp(varname,"id") ||
	   !strcmp(varname,"t") ||
	   !strcmp(varname,"err") ||
	   !strcmp(varname,"mag")) {
	  error2(ERR_RESERVEDVARIABLENAME,varname);
	}
	SetupInListVariable(p, varname, column, datatype, format);
	varname[0] = '\0';
	column = -1;
	datatype = -1;
	format[0] = '\0';
	termscan = 0;
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

void SetupInListVariable(ProgramData *p, char *varname, int column, int datatype, char *format)
{
  _Variable *variable;

  variable = CreateVariable(p, varname, (char) datatype, VARTOOLS_VECTORTYPE_INLIST, NULL);

  RegisterDataFromInputList(p,
			     p->DefinedVariables[p->NDefinedVariables - 1]->dataptr,
			     datatype,
			     0, 0, 0, 0, format,
			     (column == 0 ? column - 1 : column),
			     varname);
}


void RegisterDataFromInputList(ProgramData *p, void *dataptr, int datatype,
			       int Ncolumns, int cnum, int disjointcolumns,
			       int Nonuniformnames, char *scanformat, ...)
/* dataptr = ptr to the vector or array where the data should be stored.
      Each row in the vector or array corresponds to a light curve.

   datatype = VARTOOLS_TYPE_DOUBLE, VARTOOLS_TYPE_STRING ...

   Ncolumns = The number of columns in the array. If this is 0 then dataptr
      is a vector (e.g. (double *)), if it is > 0 then it is an array
      (e.g. (double **)).

   cnum = command number (this is appended to the column name, if it is < 0
        then it is not appended.)

   disjointcolumns = 1 if you will be specifying the list columns for each of
       the columns in the array separately, 0 if they are sequential.

   Nonuniformnames = 1 if you will be specifying separate names for each
       column in the array, 0 to give them the same base name, with the column
       index appended.

   scanformat = scanf format string to use in parsing the column, give NULL to
       use the default format.

   After scanformat, give the column(s) in the list to read in
       (== 0 the program will take the next column in the
       list, < 0 the data will not be read from the list, but memory will
       be allocated.)

   Then give the base column name(s).
*/
{
  va_list varlist;
  char *columnname;
  int incolumn;
  int i, Nmalloc;
  _DataFromInputList *d;

  if(Ncolumns < 0)
    error(ERR_CODEERROR);

  va_start(varlist, scanformat);

  if(!p->NDataFromInputList) {
    if((p->DataFromInputList = (_DataFromInputList *) malloc(sizeof(_DataFromInputList))) == NULL)
      error(ERR_MEMALLOC);
  }
  else {
    if((p->DataFromInputList = (_DataFromInputList *) realloc(p->DataFromInputList, (p->NDataFromInputList + 1)*sizeof(_DataFromInputList))) == NULL)
      error(ERR_MEMALLOC);
  }

  d = &(p->DataFromInputList[p->NDataFromInputList]);

  p->NDataFromInputList += 1;

  d->datatype = datatype;
  d->dataptr = dataptr;
  d->Ncolumns = Ncolumns;
  d->expression = NULL;
  if(scanformat == NULL ? 1 : scanformat[0] == '\0') d->scanformat = NULL;
  else {
    d->scanformat = (char *) malloc((strlen(scanformat)+1));
    sprintf(d->scanformat,"%s",scanformat);
  }

  Nmalloc = Ncolumns;
  if(Nmalloc == 0)
    Nmalloc = 1;

  if((d->incolumns = (int *) malloc(Nmalloc*sizeof(int))) == NULL)
    error(ERR_MEMALLOC);

  if((d->incolumn_names = (char **) malloc(Nmalloc*sizeof(char *))) == NULL)
    error(ERR_MEMALLOC);

  for(i=0; i < Nmalloc; i++) {
    if((d->incolumn_names[i] = (char *) malloc(MAXLEN * sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
  }

  if(disjointcolumns) {
    for(i=0; i < Ncolumns; i++) {
      incolumn = va_arg(varlist, int);
      if(incolumn <= 0) {
	d->incolumns[i] = p->inputcolumn_iter_index + 1;
	p->inputcolumn_iter_index += 1;
      }
      else
	d->incolumns[i] = incolumn;
      if(d->incolumns[i] > p->maxinputcolumn)
	p->maxinputcolumn = d->incolumns[i];
    }
  } else {
    incolumn = va_arg(varlist, int);
    if(incolumn == 0) {
      d->incolumns[0] = p->inputcolumn_iter_index + 1;
      if(Ncolumns == 0)
	p->inputcolumn_iter_index += 1;
      else
	p->inputcolumn_iter_index += Ncolumns;
    } else {
      d->incolumns[0] = incolumn;
    }
    if(d->incolumns[0] >= 0) {
      for(i=1; i < Ncolumns; i++) {
	d->incolumns[i] = d->incolumns[i-1] + 1;
      }
      if(d->incolumns[i-1] > p->maxinputcolumn)
	p->maxinputcolumn = d->incolumns[i-1];
    }
  }

  if(Nonuniformnames) {
    if(!Ncolumns) {
      columnname = va_arg(varlist, char *);
      if(cnum >= 0) {
	sprintf(d->incolumn_names[i],"%s_%d", columnname, cnum);
      } else {
	sprintf(d->incolumn_names[i],"%s", columnname);
      }
    }
    else {
      for(i=0; i < Ncolumns; i++) {
	columnname = va_arg(varlist, char *);
	if(cnum >= 0) {
	  sprintf(d->incolumn_names[i],"%s_%d", columnname, cnum);
	} else {
	  sprintf(d->incolumn_names[i],"%s", columnname);
	}
      }
    }
  } else {
    columnname = va_arg(varlist, char *);
    if(Ncolumns <= 1) {
      if(cnum >= 0) {
	sprintf(d->incolumn_names[0],"%s_%d", columnname, cnum);
      } else {
	sprintf(d->incolumn_names[0],"%s", columnname);
      }
    } else {
      for(i=0; i < Ncolumns; i++) {
	if(cnum >= 0) {
	  sprintf(d->incolumn_names[i],"%s_%d_%d", columnname, (i+1), cnum);
	} else {
	  sprintf(d->incolumn_names[i],"%s_%d", columnname, (i+1));
	}
      }
    }
  }
}

void MemAllocDataFromInputList(ProgramData *p, int Nlc) {
  _DataFromInputList *d;
  int i, j, k, Nc;
  double **dblptr;
  double ***dbl2ptr;
  short **shortptr;
  short ***short2ptr;
  int **intptr;
  int ***int2ptr;
  char **charptr;
  char ***char2ptr;
  char ***stringptr;
  char ****string2ptr;
  float **floatptr;
  float ***float2ptr;
  long **longptr;
  long ***long2ptr;
  for(i=0; i < p->NDataFromInputList; i++) {
    d = &(p->DataFromInputList[i]);
    Nc = d->Ncolumns;
    if(Nc == 0) {
      switch(d->datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dblptr = (double **) d->dataptr;
	if(((*dblptr) = (double *) malloc(Nlc * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
	for(k=0; k < Nlc; k++) (*dblptr)[k] = 0.;
	break;
      case VARTOOLS_TYPE_STRING:
	stringptr = (char ***) d->dataptr;
	if(((*stringptr) = (char **) malloc(Nlc * sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nlc; j++) {
	  if((((*stringptr)[j]) = (char *) malloc(MAXLEN)) == NULL)
	    error(ERR_MEMALLOC);
	  (*stringptr)[j][0] = '\0';
	}
	break;
      case VARTOOLS_TYPE_INT:
	intptr = (int **) d->dataptr;
	if(((*intptr) = (int *) malloc(Nlc * sizeof(int))) == NULL)
	  error(ERR_MEMALLOC);
	for(k=0; k < Nlc; k++) (*intptr)[k] = 0;
	break;
      case VARTOOLS_TYPE_SHORT:
	shortptr = (short **) d->dataptr;
	if(((*shortptr) = (short *) malloc(Nlc * sizeof(short))) == NULL)
	  error(ERR_MEMALLOC);
	for(k=0; k < Nlc; k++) (*shortptr)[k] = 0;
	break;
      case VARTOOLS_TYPE_FLOAT:
	floatptr = (float **) d->dataptr;
	if(((*floatptr) = (float *) malloc(Nlc * sizeof(float))) == NULL)
	  error(ERR_MEMALLOC);
	for(k=0; k < Nlc; k++) (*floatptr)[k] = 0.;
	break;
      case VARTOOLS_TYPE_LONG:
	longptr = (long **) d->dataptr;
	if(((*longptr) = (long *) malloc(Nlc * sizeof(long))) == NULL)
	  error(ERR_MEMALLOC);
	for(k=0; k < Nlc; k++) (*longptr)[k] = 0;
	break;
      case VARTOOLS_TYPE_CHAR:
	charptr = (char **) d->dataptr;
	if(((*charptr) = (char *) malloc(Nlc)) == NULL)
	  error(ERR_MEMALLOC);
	for(k=0; k < Nlc; k++) (*charptr)[k] = 0;
	break;
      default:
	error(ERR_BADTYPE);
      }
    } else {
      switch(d->datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dbl2ptr = (double ***) d->dataptr;
	if(((*dbl2ptr) = (double **) malloc(Nlc * sizeof(double *))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j<Nlc; j++) {
	  if((((*dbl2ptr)[j]) = (double *) malloc(Nc * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  for(k=0; k < Nc; k++) (*dbl2ptr)[j][k] = 0.;
	}
	break;
      case VARTOOLS_TYPE_STRING:
	string2ptr = (char ****) d->dataptr;
	if(((*string2ptr) = (char ***) malloc(Nlc * sizeof(char **))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nlc; j++) {
	  if((((*string2ptr)[j]) = (char **) malloc(Nc * sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(k=0; k < Nc; k++) {
	    if((((*string2ptr)[j][k]) = (char *) malloc(MAXLEN)) == NULL)
	      error(ERR_MEMALLOC);
	    (*string2ptr)[j][k][0] = '\0';
	  }
	}
	break;
      case VARTOOLS_TYPE_INT:
	int2ptr = (int ***) d->dataptr;
	if(((*int2ptr) = (int **) malloc(Nlc * sizeof(int *))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nlc; j++) {
	  if((((*int2ptr)[j]) = (int *) malloc(Nc * sizeof(int))) == NULL)
	    error(ERR_MEMALLOC);
	  for(k=0; k < Nc; k++) (*int2ptr)[j][k] = 0;
	}
	break;
      case VARTOOLS_TYPE_SHORT:
	short2ptr = (short ***) d->dataptr;
	if(((*short2ptr) = (short **) malloc(Nlc * sizeof(short *))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nlc; j++) {
	  if((((*short2ptr)[j]) = (short *) malloc(Nc * sizeof(short))) == NULL)
	    error(ERR_MEMALLOC);
	  for(k=0; k < Nc; k++) (*short2ptr)[j][k] = 0;
	}
	break;
      case VARTOOLS_TYPE_FLOAT:
	float2ptr = (float ***) d->dataptr;
	if(((*float2ptr) = (float **) malloc(Nlc * sizeof(float *))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nlc; j++) {
	  if((((*float2ptr)[j]) = (float *) malloc(Nc * sizeof(float))) == NULL)
	    error(ERR_MEMALLOC);
	  for(k=0; k < Nc; k++) (*float2ptr)[j][k] = 0.;
	}
	break;
      case VARTOOLS_TYPE_LONG:
	long2ptr = (long ***) d->dataptr;
	if(((*long2ptr) = (long **) malloc(Nlc * sizeof(long *))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nlc; j++) {
	  if((((*long2ptr)[j]) = (long *) malloc(Nc * sizeof(long))) == NULL)
	    error(ERR_MEMALLOC);
	  for(k=0; k < Nc; k++) (*long2ptr)[j][k] = 0;
	}
	break;
      case VARTOOLS_TYPE_CHAR:
	char2ptr = (char ***) d->dataptr;
	if(((*char2ptr) = (char **) malloc(Nlc * sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nlc; j++) {
	  if((((*char2ptr)[j]) = (char *) malloc(Nc)) == NULL)
	    error(ERR_MEMALLOC);
	  for(k=0; k < Nc; k++) (*char2ptr)[j][k] = 0;
	}
	break;
      default:
	error(ERR_BADTYPE);

      }
    }
  }
}

void ParseLineToColumns(char *line, char **cols, int maxcols)
{
  int j, i;
  j = 0;
  i = 0;
  while(i < maxcols && line[j] != '\0' && line[j] != '\n') {
    j += parseone(&line[j],(void *) ((cols[i])), VARTOOLS_TYPE_STRING);
    i++;
  }
  if(i < maxcols) {
    error2(ERR_INPUTMISSINGCOLUMN,"Input List");
  }
}

/* Parse the input list */
void ParseInputList(ProgramData *p, char **inputlines, int Nlcs)
{
  int Np, Nc, i, j, k, u, in_indx;
  char **incols;
  _DataFromInputList *d;

  double **dblptr;
  double ***dbl2ptr;
  int **intptr;
  int ***int2ptr;
  short **shortptr;
  short ***short2ptr;
  char **charptr;
  char ***char2ptr;
  char ***stringptr;
  char ****string2ptr;
  float **floatptr;
  float ***float2ptr;
  long **longptr;
  long ***long2ptr;

  int Nskip;

  if(p->Ncopycommands > 0) {
    Nskip = p->Ncopiestotal;
  } else {
    Nskip = 1;
  }

  Np = p->maxinputcolumn;

  if(Np <= 0)
    return;

  if((incols = (char **) malloc(Np * sizeof(char *))) == NULL)
    error(ERR_MEMALLOC);
  for(i = 0; i < Np; i++) {
    if((incols[i] = (char *) malloc(MAXLEN)) == NULL)
      error(ERR_MEMALLOC);
  }

  for(i=0, in_indx = 0; i < Nlcs; i += Nskip, in_indx += 1) {
    ParseLineToColumns(inputlines[in_indx],incols,Np);
    for(j = 0; j < p->NDataFromInputList; j++) {
      d = &(p->DataFromInputList[j]);
      Nc = d->Ncolumns;
      if(Nc == 0) {
	if(d->incolumns[0] <= 0)
	  continue;
	k = d->incolumns[0] - 1;
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = (double **) d->dataptr;
	  (*dblptr)[i] = atof(incols[k]);
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = (char ***) d->dataptr;
	  sprintf(((*stringptr)[i]),"%s",incols[k]);
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = (int **) d->dataptr;
	  (*intptr)[i] = atoi(incols[k]);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = (short **) d->dataptr;
	  (*shortptr)[i] = (short) atoi(incols[k]);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = (float **) d->dataptr;
	  (*floatptr)[i] = atof(incols[k]);
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = (long **) d->dataptr;
	  (*longptr)[i] = atol(incols[k]);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = (char **) d->dataptr;
	  (*charptr)[i] = incols[k][0];
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      } else {
	for(u = 0; u < Nc; u++) {
	  if(d->incolumns[u] <= 0)
	    continue;
	  k = d->incolumns[u] - 1;
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dbl2ptr = (double ***) d->dataptr;
	    (*dbl2ptr)[i][u] = atof(incols[k]);
	    break;
	  case VARTOOLS_TYPE_STRING:
	    string2ptr = (char ****) d->dataptr;
	    sprintf(((*string2ptr)[i][u]),"%s",incols[k]);
	    break;
	  case VARTOOLS_TYPE_INT:
	    int2ptr = (int ***) d->dataptr;
	    (*int2ptr)[i][u] = atoi(incols[k]);
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    short2ptr = (short ***) d->dataptr;
	    (*short2ptr)[i][u] = (short) atoi(incols[k]);
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    float2ptr = (float ***) d->dataptr;
	    (*float2ptr)[i][u] = atof(incols[k]);
	    break;
	  case VARTOOLS_TYPE_LONG:
	    long2ptr = (long ***) d->dataptr;
	    (*long2ptr)[i][u] = atol(incols[k]);
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    char2ptr = (char ***) d->dataptr;
	    (*char2ptr)[i][u] = incols[k][0];
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      }
    }
  }

  /* Fill in any variables calculated from analytic expressions */
  for(j = 0; j < p->NDataFromInputList; j++) {
    d = &(p->DataFromInputList[j]);
    if(d->incolumns[0] <= 0)
      {
	if(d->scanformat != NULL) {
	  Nc = d->Ncolumns;
	  if(Nc != 0)
	    error(ERR_CODEERROR);
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    for(i=0; i < Nlcs; i += Nskip) {
	      dblptr = (double **) d->dataptr;
	      (*dblptr)[i] = EvaluateExpression(i, 0, 0, d->expression);
	    }
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    for(i=0; i < Nlcs; i += Nskip) {
	      floatptr = (float **) d->dataptr;
	      (*floatptr)[i] = EvaluateExpression(i, 0, 0, d->expression);
	    }
	    break;
	  case VARTOOLS_TYPE_INT:
	    for(i=0; i < Nlcs; i += Nskip) {
	      intptr = (int **) d->dataptr;
	      (*intptr)[i] = EvaluateExpression(i, 0, 0, d->expression);
	    }
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    for(i=0; i < Nlcs; i += Nskip) {
	      shortptr = (short **) d->dataptr;
	      (*shortptr)[i] = EvaluateExpression(i, 0, 0, d->expression);
	    }
	    break;
	  case VARTOOLS_TYPE_LONG:
	    for(i=0; i < Nlcs; i += Nskip) {
	      longptr = (long **) d->dataptr;
	      (*longptr)[i] = EvaluateExpression(i, 0, 0, d->expression);
	    }
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	} else {
	  Nc = d->Ncolumns;
	  if(Nc == 0) {
	    switch(d->datatype) {
	    case VARTOOLS_TYPE_DOUBLE:
	      for(i=0; i < Nlcs; i += Nskip) {
		dblptr = (double **) d->dataptr;
		(*dblptr)[i] = 0.0;
	      }
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      for(i=0; i < Nlcs; i += Nskip) {
		floatptr = (float **) d->dataptr;
		(*floatptr)[i] = 0.0;
	      }
	      break;
	    case VARTOOLS_TYPE_INT:
	      for(i=0; i < Nlcs; i += Nskip) {
		intptr = (int **) d->dataptr;
		(*intptr)[i] = 0;
	      }
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      for(i=0; i < Nlcs; i += Nskip) {
		shortptr = (short **) d->dataptr;
		(*shortptr)[i] = 0;
	      }
	      break;
	    case VARTOOLS_TYPE_LONG:
	      for(i=0; i < Nlcs; i += Nskip) {
		longptr = (long **) d->dataptr;
		(*longptr)[i] = 0;
	      }
	      break;
	    default:
	      error(ERR_BADTYPE);
	    }
	  }
	  else {
	    for(u = 0; u < Nc; u++) {
	      if(d->incolumns[u] > 0)
		continue;
	      switch(d->datatype) {
	      case VARTOOLS_TYPE_DOUBLE:
		dbl2ptr = (double ***) d->dataptr;
		for(i=0; i < Nlcs; i += Nskip) {
		  (*dbl2ptr)[i][u] = 0.0;
		}
		break;
	      case VARTOOLS_TYPE_STRING:
 		string2ptr = (char ****) d->dataptr;
		for(i=0; i < Nlcs; i += Nskip) {
		  sprintf(((*string2ptr)[i][u]),"");
		}
		break;
	      case VARTOOLS_TYPE_INT:
		int2ptr = (int ***) d->dataptr;
		for(i=0; i < Nlcs; i += Nskip) {
		  (*int2ptr)[i][u] = 0;
		}
		break;
	      case VARTOOLS_TYPE_SHORT:
		short2ptr = (short ***) d->dataptr;
		for(i=0; i < Nlcs; i += Nskip) {
		  (*short2ptr)[i][u] = 0;
		}
		break;
	      case VARTOOLS_TYPE_FLOAT:
		float2ptr = (float ***) d->dataptr;
		for(i=0; i < Nlcs; i += Nskip) {
		  (*float2ptr)[i][u] = 0.0;
		}
		break;
	      case VARTOOLS_TYPE_LONG:
		long2ptr = (long ***) d->dataptr;
		for(i=0; i < Nlcs; i += Nskip) {
		  (*long2ptr)[i][u] = 0;
		}
		break;
	      case VARTOOLS_TYPE_CHAR:
		char2ptr = (char ***) d->dataptr;
		for(i=0; i < Nlcs; i += Nskip) {
		  (*char2ptr)[i][u] = '\0';
		}
		break;
	      default:
		error(ERR_BADTYPE);
	      }
	    }
	  }
	}
      }
  }

  for(i = 0; i < Np; i++)
    free(incols[i]);
  free(incols);

}

void MoveInputListData(ProgramData *p, int isrc, int idest) {
  _DataFromInputList *d;
  int i, j, k, Nc;
  double **dblptr;
  double ***dbl2ptr;
  short **shortptr;
  short ***short2ptr;
  int **intptr;
  int ***int2ptr;
  char **charptr;
  char ***char2ptr;
  char ***stringptr;
  char ****string2ptr;
  float **floatptr;
  float ***float2ptr;
  long **longptr;
  long ***long2ptr;
  for(i=0; i < p->NDataFromInputList; i++) {
    d = &(p->DataFromInputList[i]);
    Nc = d->Ncolumns;
    if(Nc == 0) {
      switch(d->datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dblptr = (double **) d->dataptr;
	(*dblptr)[idest] = (*dblptr)[isrc];
	break;
      case VARTOOLS_TYPE_STRING:
	stringptr = (char ***) d->dataptr;
	(*stringptr)[idest] = (*stringptr)[isrc];
	break;
      case VARTOOLS_TYPE_INT:
	intptr = (int **) d->dataptr;
	(*intptr)[idest] = (*intptr)[isrc];
	break;
      case VARTOOLS_TYPE_SHORT:
	shortptr = (short **) d->dataptr;
	(*shortptr)[idest] = (*shortptr)[isrc];
	break;
      case VARTOOLS_TYPE_FLOAT:
	floatptr = (float **) d->dataptr;
	(*floatptr)[idest] = (*floatptr)[isrc];
	break;
      case VARTOOLS_TYPE_LONG:
	longptr = (long **) d->dataptr;
	(*longptr)[idest] = (*longptr)[isrc];
	break;
      case VARTOOLS_TYPE_CHAR:
	charptr = (char **) d->dataptr;
	(*charptr)[idest] = (*charptr)[isrc];
	break;
      default:
	error(ERR_BADTYPE);
      }
    } else {
      switch(d->datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dbl2ptr = (double ***) d->dataptr;
	(*dbl2ptr)[idest] = (*dbl2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_STRING:
	string2ptr = (char ****) d->dataptr;
	(*string2ptr)[idest] = (*string2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_INT:
	int2ptr = (int ***) d->dataptr;
	(*int2ptr)[idest] = (*int2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_SHORT:
	short2ptr = (short ***) d->dataptr;
	(*short2ptr)[idest] = (*short2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_FLOAT:
	float2ptr = (float ***) d->dataptr;
	(*float2ptr)[idest] = (*float2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_LONG:
	long2ptr = (long ***) d->dataptr;
	(*long2ptr)[idest] = (*long2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_CHAR:
	char2ptr = (char ***) d->dataptr;
	(*char2ptr)[idest] = (*char2ptr)[isrc];
	break;
      default:
	error(ERR_BADTYPE);

      }
    }
  }
}


/* Show the expected format of the input list */
void printinputlistformat(ProgramData *p, FILE *outfile)
{
  int i, j, Np, Nc, u;
  _DataFromInputList *d;

  fprintf(outfile,
	  "#Expected Format Of Input List\n");
  fprintf(outfile,
	  "#-----------------------------\n");

  Np = p->maxinputcolumn;

  for(i=1; i <= Np; i++) {
    for(j=0; j < p->NDataFromInputList; j++) {
      d = &(p->DataFromInputList[j]);
      Nc = d->Ncolumns;
      if(Nc == 0) {
	if(d->incolumns[0] <= 0)
	  continue;
	if(d->incolumns[0] == i) {
	  fprintf(outfile,
		  "Column %d: %s\n",
		  i, d->incolumn_names[0]);
	}
      } else {
	for(u=0; u < Nc; u++) {
	  if(d->incolumns[u] <= 0)
	    continue;
	  if(d->incolumns[u] == i) {
	    fprintf(outfile,
		    "Column %d: %s\n",
		    i, d->incolumn_names[u]);
	  }
	}
      }
    }
  }
}
