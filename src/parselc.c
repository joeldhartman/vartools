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

void ParseSkipCharString(char *argv, int *Nskip, char **skipchars)
{
  int i, j, nfound, lastindex;
  if(!strcmp(argv,"None") || !strcmp(argv,"none")) {
    *Nskip = 0;
    return;
  }
  *Nskip = 0;
  lastindex = -1;
  i = 0;
  while(argv[i] != '\0') {
    if(argv[i] != ' ' && argv[i] != '\t' && argv[i] != '\n') {
      if(argv[i] != ',') {
	if(i == 0 ? 1 : (argv[i-1] == ',' || argv[i-1] == ' ' || argv[i-1] == '\t' || argv[i-1] == '\n')) {
	  if((*Nskip) == 0) {
	    *Nskip = 1;
	    if(((*skipchars) = (char *) malloc(1)) == NULL)
	      error(ERR_MEMALLOC);
	  } else {
	    (*Nskip) = (*Nskip) + 1;
	    if(((*skipchars) = (char *) realloc((*skipchars),(*Nskip))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  (*skipchars)[(*Nskip)-1] = argv[i];
	}
      }
    }
    i++;
  }
}

void SetupInputVariable(ProgramData *p, char *varname, int column, int datatype, char *format, char *incolumn_header_name)
{
  _Variable *variable;

  variable = CreateVariable(p, varname, (char) datatype, VARTOOLS_VECTORTYPE_LC, NULL);

  /* Check the variable name for special cases */
  if(!strcmp(varname,"t")) {
    if(datatype != VARTOOLS_TYPE_DOUBLE && 
       datatype != VARTOOLS_TYPE_CONVERTJD) {
      error(ERR_WRONGTYPEFORSPECIALVARIABLE);
    }
    p->tvar = variable;
  }
  else if(!strcmp(varname,"mag")) {
    if(datatype != VARTOOLS_TYPE_DOUBLE) {
      error(ERR_WRONGTYPEFORSPECIALVARIABLE);
    }
    p->magvar = variable;
  }
  else if(!strcmp(varname,"err")) {
    if(datatype != VARTOOLS_TYPE_DOUBLE) {
      error(ERR_WRONGTYPEFORSPECIALVARIABLE);
    }
    p->sigvar = variable;
  }
  else if(!strcmp(varname,"id")) {
    if(datatype != VARTOOLS_TYPE_STRING) {
      error(ERR_WRONGTYPEFORSPECIALVARIABLE);
    }
    p->idvar = variable;
  }

#ifdef _USEBINARY_LC
  if(incolumn_header_name == NULL ? 1 : incolumn_header_name[0] == '\0') {
#endif
    RegisterDataFromLightCurve(p, 
			     p->DefinedVariables[p->NDefinedVariables-1]->dataptr, 
			     datatype, 
			     (datatype == VARTOOLS_TYPE_STRING ? MAXLEN : 0), 
			     0, 0, 0, 0, format, variable,
			     (column == 0 ? column - 1 : column), 
			     varname);
#ifdef _USEBINARY_LC
  } else {
    RegisterDataFromLightCurve(p, 
			       p->DefinedVariables[p->NDefinedVariables-1]->dataptr, 
			       datatype, 
			       (datatype == VARTOOLS_TYPE_STRING ? MAXLEN : 0), 
			       0, 0, 0, 0, format, variable,
			       INCOLUMN_HEADER_INDICATOR,
			       varname);
    if((p->DataFromLightCurve[p->NDataFromLightCurve-1].incolumn_header_names = (char **) malloc(sizeof(char *))) == NULL)
      error(ERR_MEMALLOC);
    if((p->DataFromLightCurve[p->NDataFromLightCurve-1].incolumn_header_names[0] = (char *) malloc(MAXLEN)) == NULL)
      error(ERR_MEMALLOC);
    sprintf(p->DataFromLightCurve[p->NDataFromLightCurve-1].incolumn_header_names[0], "%s", incolumn_header_name);
  }
#endif
}

char isstringint(char *s)
{
  int i;
  if(s == NULL)
    return 0;
  if(s[0] == '\0')
    return 0;
  i = 0;
  while(s[i] != '\0') {
    if(s[i] < '0' || s[i] > '9')
      return 0;
    i++;
  }
  return 1;
}

/* This function parses a string of the form:
   var1:col1[:type1[:fmt1][,var2:col2[:type2[:fmt2]],...]]
   and creates the relevant _DataFromLightCurve and _Variable objects */
int ParseInputLCFormatString(char *argv, ProgramData *p)
{
  int i, j, termscan, startindex;
  char *varname;
  int column;
  int datatype;
  char *format;
  char *parsecopy;
  char *incolumn_header_name = NULL;
  int sizestring;

  p->coljd = -1;
  p->colmag = -1;
  p->colsig = -1;

  sizestring = MAXLEN;
  format = (char *) malloc(sizestring);
  parsecopy = (char *) malloc(sizestring);
  varname = (char *) malloc(sizestring);

  varname[0] = '\0'; column = -1; datatype = -1; format[0] = '\0';
  termscan = 0;

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
	  p->lc_getcolumnsfromheader = 1;
	}
	else {
	  column = atoi(parsecopy);
	}
#elif USECFITSIO
	if(!isstringint(parsecopy)) {
	  incolumn_header_name = malloc(strlen(parsecopy)+1);
	  sprintf(incolumn_header_name,"%s",parsecopy);
	  column = 1;
	  p->lc_getcolumnsfromheader = 1;
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
	if(!strcmp(varname,"id")) {
	  p->readimagestring = 1;
	}
	if(!strcmp(varname,"t")) {
	  p->coljd = column;
	  if(column == 0 && format[0] == '\0') {
	    sprintf(format,"NR");
	  }
	}
	else if(!strcmp(varname,"err")) {
	  p->colsig = column;
	  if(column == 0 && format[0] == '\0') {
	    sprintf(format,"1.0");
	  }
	}
	else if(!strcmp(varname,"mag")) {
	  p->colmag = column;
	}
	SetupInputVariable(p, varname, column, datatype, format, incolumn_header_name);
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

int ParseLineToColumns_testskip(char *line, char **cols, int maxcols, int Nskipchar, char *skipchars)
{
  int j, i, k;
  j = 0;
  i = 0;
  /* Find the first non-white space character on the line and check if it is
     a character to reject */
  while((line[j] == ' ' || line[j] == '\t') && line[j] != '\0' && line[j] != '\n')
    j++;
  /* Return 1 if the line is blank, or if the first character is
   a comment */
  for(k=0; k < Nskipchar; k++) {
    if(line[j] == skipchars[k]) return 1;
  }
  if((line[j] == '\n' || line[j] == '\0') && maxcols > 0)
    return 1;

  while(i < maxcols && line[j] != '\0' && line[j] != '\n') {
    j += parseone(&line[j],(void *) ((cols[i])), VARTOOLS_TYPE_STRING);
    i++;
  }
  if(i < maxcols) {
    error2(ERR_INPUTMISSINGCOLUMN,"Input List");
  }
  return 0;
}

int ParseLineToColumnsDelimString_testskip(char *line, char **cols, int maxcols, int Nskipchar, char *skipchars, char *delim)
{
  int j, i, k;
  j = 0;
  i = 0;
  /* Return 1 if the line is blank, or if the first character is
   a comment */
  for(k=0; k < Nskipchar; k++) {
    if(line[j] == skipchars[k]) return 1;
  }
  if((line[j] == '\n' || line[j] == '\0') && maxcols > 0)
    return 1;

  while(i < maxcols && line[j] != '\0' && line[j] != '\n') {
    j += parseonedelimstring(&line[j],(void *) ((cols[i])), VARTOOLS_TYPE_STRING, delim);
    i++;
  }
  if(i < maxcols) {
    error2(ERR_INPUTMISSINGCOLUMN,"Input List");
  }
  return 0;
}

int ParseLineToColumnsDelimChar_testskip(char *line, char **cols, int maxcols, int Nskipchar, char *skipchars, char delim)
{
  int j, i, k;
  j = 0;
  i = 0;
  /* Return 1 if the line is blank, or if the first character is
   a comment */
  for(k=0; k < Nskipchar; k++) {
    if(line[j] == skipchars[k]) return 1;
  }
  if((line[j] == '\n' || line[j] == '\0') && maxcols > 0)
    return 1;

  while(i < maxcols && line[j] != '\0' && line[j] != '\n') {
    j += parseonedelimchar(&line[j],(void *) ((cols[i])), VARTOOLS_TYPE_STRING, delim);
    i++;
  }
  if(i < maxcols) {
    error2(ERR_INPUTMISSINGCOLUMN,"Input List");
  }
  return 0;
}

void RegisterDataFromLightCurve(ProgramData *p, void *dataptr, int datatype,
				int maxstringlength, int Ncolumns, int cnum, 
				int disjointcolumns,
				int Nonuniformnames, char *scanformat, 
				_Variable *variable, ...)
/* dataptr = ptr to an array of vectors, or array of 2-d matrices,
      where the data should be stored. The array will be indexed as
      follows:  dataptr[Nthread][NJD], dataptr[Nthread][Nc][NJD]
         or dataptr[Nthread][NJD][Nc] for Ncolumns = 0, Ncolumns > 0 and
         Ncolumns < 0, respectively.
      Here Nthread is the number of light curves that will be
      simultaneously read-in, NJD is the number of points in a light
      curve, and Nc is the number of columns in the matrix. VARTOOLS
      will determine Nthread, and NJD, and allocate memory. One should
      declare a variable of type ** (if it will store a vector), or
      *** (if it will store a matrix) and then pass a pointer to it
      cast to (void *).

   datatype = VARTOOLS_TYPE_DOUBLE, VARTOOLS_TYPE_STRING ...

   Ncolumns = The number of columns in the matrix. If this is 0 then
      dataptr is an array of vectors (e.g. (double **)), if it is > 0 or < 0
      then it is an array of matrices (e.g. (double ***)). If it is > 0, then
      dataptr will be indexed as dataptr[Nthread][Nc][NJD], if it is < 0, then
      dataptr will be indexed as dataptr[Nthread][NJD][Nc].

   cnum = command number (this is appended to the column name, if it is < 0
      then it is not appended.)

   disjointcolumns = 1 if you will be specifying the list columns for
       each of the columns in the array as separate arguments to the
       function, 2 if you will be supplying them in a vector, or 0 if
       they are sequential.

   Nonuniformnames = 1 if you will be specifying separate names for
       each column in the array as separate arguments to the function,
       2 if you will be supplying them in an array of strings, 0 to
       give them the same base name, with the column index appended.

   scanformat = A scanf format string to use in reading in the 
       light curve. If it is NULL, then the full input column string
       will be converted directly to the appropriate data type.  If 
       the column number is < 0 (data not read-in from the light curve) and
       the format is not NULL, then it will be treated as an analytic 
       expression which the vector is to be initialized to.

   variable = a pointer to an existing variable object to associate with this
       column. This can be NULL.

   After variable, give the column(s) in the list to read in
       (== 0 the program will take the next column in the list.; < 0
       the program will not read-in the data from the light curve, but
       it will allocate memory to store the data; If disjointcolumns
       == 0, you only need to provide a single number; if == 1 you
       should provide one number per column as separate arguments to
       the function; == 2 you should pass an int * vector with the
       list)

   Then give the base column name(s) (Nonuniformnames == 0 just give 1
       name; Nonuniformnames == 1 give the name for each column as a
       separate argument to the function; Nonuniformnames == 2, you
       should pass a char ** array with the list).
*/
{
  va_list varlist;
  char *columnname, **columnnamevec;
  int incolumn, *incolumnvec;
  int i, Nmalloc;
  _DataFromLightCurve *d;
  
  va_start(varlist, variable);
  
  if(!p->NDataFromLightCurve) {
    if((p->DataFromLightCurve = (_DataFromLightCurve *) malloc(sizeof(_DataFromLightCurve))) == NULL)
      error(ERR_MEMALLOC);
  }
  else {
    if((p->DataFromLightCurve = (_DataFromLightCurve *) realloc(p->DataFromLightCurve, (p->NDataFromLightCurve + 1)*sizeof(_DataFromLightCurve))) == NULL)
      error(ERR_MEMALLOC);
  }
  
  d = &(p->DataFromLightCurve[p->NDataFromLightCurve]);
  d->maxstringlength = maxstringlength;
  p->NDataFromLightCurve += 1;

  d->datatype = datatype;
  d->dataptr = dataptr;
  d->Ncolumns = Ncolumns;
  d->expression = NULL;
  d->incolumn_header_names = NULL;

  if(scanformat == NULL ? 1 : scanformat[0] == '\0') d->scanformat = NULL;
  else {
    d->scanformat = (char *) malloc((strlen(scanformat)+1));
    sprintf(d->scanformat,"%s",scanformat);
  }
  d->variable = variable;

  if(datatype == VARTOOLS_TYPE_CONVERTJD) {
    checkUTCFormat(d->scanformat, d->UTCindex);
  }

  Nmalloc = abs(Ncolumns);
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
  
  if(disjointcolumns == 1) {
    for(i=0; i < abs(Ncolumns); i++) {
      incolumn = va_arg(varlist, int);
      if(incolumn == 0) {
	d->incolumns[i] = p->inputlccolumn_iter_index + 1;
	p->inputlccolumn_iter_index += 1;
      }
      else
	d->incolumns[i] = incolumn;
      if(d->incolumns[i] > p->maxinputlccolumn)
	p->maxinputlccolumn = d->incolumns[i];
    }
  }  else if(disjointcolumns == 2) {
    incolumnvec = va_arg(varlist, int *);
    for(i=0; i < abs(Ncolumns); i++) {
      if(incolumnvec[i] == 0) {
	d->incolumns[i] = p->inputlccolumn_iter_index + 1;
	p->inputlccolumn_iter_index += 1;
      }
      else
	d->incolumns[i] = incolumnvec[i];
      if(d->incolumns[i] > p->maxinputlccolumn)
	p->maxinputlccolumn = d->incolumns[i];
    }
  }
  else {
    incolumn = va_arg(varlist, int);
    if(incolumn == 0) {
      d->incolumns[0] = p->inputlccolumn_iter_index + 1;
      if(Ncolumns == 0)
	p->inputlccolumn_iter_index += 1;
      else
	p->inputlccolumn_iter_index += abs(Ncolumns);
    } else {
      d->incolumns[0] = incolumn;
    }
    for(i=1; i < abs(Ncolumns); i++) {
      d->incolumns[i] = d->incolumns[i-1] + 1;
    }
    if(d->incolumns[i-1] > p->maxinputlccolumn)
      p->maxinputlccolumn = d->incolumns[i-1];
  }

  if(incolumn < 0)
    return;

  if(Nonuniformnames == 1) {
    if(!Ncolumns) {
      columnname = va_arg(varlist, char *);
      if(cnum >= 0) {
	sprintf(d->incolumn_names[i],"%s_%d", columnname, cnum);
      } else {
	sprintf(d->incolumn_names[i],"%s", columnname);
      }
    }
    else {
      for(i=0; i < abs(Ncolumns); i++) {
	columnname = va_arg(varlist, char *);
	if(cnum >= 0) {
	  sprintf(d->incolumn_names[i],"%s_%d", columnname, cnum);
	} else {
	  sprintf(d->incolumn_names[i],"%s", columnname);
	}
      }
    }
  } 
  else if(Nonuniformnames == 2) {
    if(!Ncolumns) {
      columnnamevec = va_arg(varlist, char **);
      if(cnum >= 0) {
	sprintf(d->incolumn_names[i],"%s_%d", columnnamevec[0], cnum);
      } else {
	sprintf(d->incolumn_names[i],"%s", columnnamevec[0]);
      }
    }
    else {
      columnnamevec = va_arg(varlist, char **);
      for(i=0; i < abs(Ncolumns); i++) {
	if(cnum >= 0) {
	  sprintf(d->incolumn_names[i],"%s_%d", columnnamevec[i], cnum);
	} else {
	  sprintf(d->incolumn_names[i],"%s", columnnamevec[i]);
	}
      }
    }
  }
  else {
    columnname = va_arg(varlist, char *);
    if(abs(Ncolumns) <= 1) {
      if(cnum >= 0) {
	sprintf(d->incolumn_names[0],"%s_%d", columnname, cnum);
      } else {
	sprintf(d->incolumn_names[0],"%s", columnname);
      }
    } else {
      for(i=0; i < abs(Ncolumns); i++) {
	if(cnum >= 0) {
	  sprintf(d->incolumn_names[i],"%s_%d_%d", columnname, (i+1), cnum);
	} else {
	  sprintf(d->incolumn_names[i],"%s_%d", columnname, (i+1));
	}
      }
    }
  }
}

void InitializeMemAllocDataFromLightCurve(ProgramData *p, int Nthread) {
  _DataFromLightCurve *d;
  int i, j, k, Nc;
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
  p->Nthread = Nthread;
  if(p->readimagestring) {
    if((p->stringid_idx = (int **) malloc(Nthread * sizeof(int *))) == NULL ||
       (p->stringid = (char ***) malloc(Nthread * sizeof(char **))) == NULL)
      error(ERR_MEMALLOC);
  }
  if((p->t = (double **) malloc(Nthread * sizeof(double *))) == NULL ||
     (p->mag = (double **) malloc(Nthread * sizeof(double *))) == NULL ||
     (p->sig = (double **) malloc(Nthread * sizeof(double *))) == NULL)
    error(ERR_MEMALLOC);
#ifdef _USEBINARY_LC
  if(p->binarylcinput) {
    if((p->binlc = (binarylightcurve *) malloc(Nthread * sizeof(binarylightcurve))) == NULL)
      error(ERR_MEMALLOC);
    for(i=0; i < Nthread; i++) {
      p->binlc[i].memsize_lc_header = 0;
      p->binlc[i].memsize_lc_columns = 0;
    }
  }
#endif
  for(i=0; i < p->NDataFromLightCurve; i++) {
    d = &(p->DataFromLightCurve[i]);
    Nc = d->Ncolumns;
    if(Nc == 0) {
      switch(d->datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dblptr = (double ***) d->dataptr;
	if(((*dblptr) = (double **) malloc(Nthread * sizeof(double *))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_CONVERTJD:
	dblptr = (double ***) d->dataptr;
	if(((*dblptr) = (double **) malloc(Nthread * sizeof(double *))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_STRING:
	stringptr = (char ****) d->dataptr;
	if(((*stringptr) = (char ***) malloc(Nthread * sizeof(char **))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_INT:
	intptr = (int ***) d->dataptr;
	if(((*intptr) = (int **) malloc(Nthread * sizeof(int *))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_SHORT:
	shortptr = (short ***) d->dataptr;
	if(((*shortptr) = (short **) malloc(Nthread * sizeof(short *))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_FLOAT:
	floatptr = (float ***) d->dataptr;
	if(((*floatptr) = (float **) malloc(Nthread * sizeof(float *))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_LONG:
	longptr = (long ***) d->dataptr;
	if(((*longptr) = (long **) malloc(Nthread * sizeof(long *))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_CHAR:
	charptr = (char ***) d->dataptr;
	if(((*charptr) = (char **) malloc(Nthread * sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      default:
	error(ERR_BADTYPE);
      }
    } else if(Nc > 0) {
      switch(d->datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dbl2ptr = (double ****) d->dataptr;
	if(((*dbl2ptr) = (double ***) malloc(Nthread * sizeof(double **))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j<Nthread; j++) {
	  if((((*dbl2ptr)[j]) = (double **) malloc(Nc * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	}
	break;
      case VARTOOLS_TYPE_CONVERTJD:
	dbl2ptr = (double ****) d->dataptr;
	if(((*dbl2ptr) = (double ***) malloc(Nthread * sizeof(double **))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j<Nthread; j++) {
	  if((((*dbl2ptr)[j]) = (double **) malloc(Nc * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	}
	break;
      case VARTOOLS_TYPE_STRING:
	string2ptr = (char *****) d->dataptr;
	if(((*string2ptr) = (char ****) malloc(Nthread * sizeof(char ***))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nthread; j++) {
	  if((((*string2ptr)[j]) = (char ***) malloc(Nc * sizeof(char **))) == NULL)
	    error(ERR_MEMALLOC);
	}
	break;
      case VARTOOLS_TYPE_INT:
	int2ptr = (int ****) d->dataptr;
	if(((*int2ptr) = (int ***) malloc(Nthread * sizeof(int **))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nthread; j++) {
	  if((((*int2ptr)[j]) = (int **) malloc(Nc * sizeof(int *))) == NULL)
	    error(ERR_MEMALLOC);
	}
	break;
      case VARTOOLS_TYPE_SHORT:
	short2ptr = (short ****) d->dataptr;
	if(((*short2ptr) = (short ***) malloc(Nthread * sizeof(short **))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nthread; j++) {
	  if((((*short2ptr)[j]) = (short **) malloc(Nc * sizeof(short *))) == NULL)
	    error(ERR_MEMALLOC);
	}
	break;
      case VARTOOLS_TYPE_FLOAT:
	float2ptr = (float ****) d->dataptr;
	if(((*float2ptr) = (float ***) malloc(Nthread * sizeof(float **))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nthread; j++) {
	  if((((*float2ptr)[j]) = (float **) malloc(Nc * sizeof(float *))) == NULL)
	    error(ERR_MEMALLOC);
	}
	break;
      case VARTOOLS_TYPE_LONG:
	long2ptr = (long ****) d->dataptr;
	if(((*long2ptr) = (long ***) malloc(Nthread * sizeof(long **))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nthread; j++) {
	  if((((*long2ptr)[j]) = (long **) malloc(Nc * sizeof(long *))) == NULL)
	    error(ERR_MEMALLOC);
	}
	break;
      case VARTOOLS_TYPE_CHAR:
	char2ptr = (char ****) d->dataptr;
	if(((*char2ptr) = (char ***) malloc(Nthread * sizeof(char **))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nthread; j++) {
	  if((((*char2ptr)[j]) = (char **) malloc(Nc * sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	}
	break;
      default:
	error(ERR_BADTYPE);
	
      }
    } else {
      switch(d->datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dbl2ptr = (double ****) d->dataptr;
	if(((*dbl2ptr) = (double ***) malloc(Nthread * sizeof(double **))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_CONVERTJD:
	dbl2ptr = (double ****) d->dataptr;
	if(((*dbl2ptr) = (double ***) malloc(Nthread * sizeof(double **))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_STRING:
	string2ptr = (char *****) d->dataptr;
	if(((*string2ptr) = (char ****) malloc(Nthread * sizeof(char ***))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_INT:
	int2ptr = (int ****) d->dataptr;
	if(((*int2ptr) = (int ***) malloc(Nthread * sizeof(int **))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_SHORT:
	short2ptr = (short ****) d->dataptr;
	if(((*short2ptr) = (short ***) malloc(Nthread * sizeof(short **))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_FLOAT:
	float2ptr = (float ****) d->dataptr;
	if(((*float2ptr) = (float ***) malloc(Nthread * sizeof(float **))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_LONG:
	long2ptr = (long ****) d->dataptr;
	if(((*long2ptr) = (long ***) malloc(Nthread * sizeof(long **))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_CHAR:
	char2ptr = (char ****) d->dataptr;
	if(((*char2ptr) = (char ***) malloc(Nthread * sizeof(char **))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      default:
	error(ERR_BADTYPE);
	
      }
    }
  }
}

void MemAllocDataFromLightCurve(ProgramData *p, int threadid, int Nterm) {
  _DataFromLightCurve *d;
  int i, j, k, Nc;
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
  int s, oldsizesinglelc;

  if(Nterm <= p->sizesinglelc[threadid]) return;

  if(p->sizesinglelc[threadid] == 0)
    {
      p->sizesinglelc[threadid] = 1024;
      while(p->sizesinglelc[threadid] < Nterm)
	p->sizesinglelc[threadid] *= 2;
      s = p->sizesinglelc[threadid];
      if(p->readimagestring) {
	if((p->stringid_idx[threadid] = (int *) malloc(s * sizeof(int))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < s; j++)
	  p->stringid_idx[threadid][j] = 0;
      }
      for(i=0; i < p->NDataFromLightCurve; i++) {
	d = &(p->DataFromLightCurve[i]);
	Nc = d->Ncolumns;
	if(Nc == 0) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dblptr = (double ***) d->dataptr;
	    if(((*dblptr)[threadid] = (double *) malloc(s * sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < s; j++) (*dblptr)[threadid][j] = 0.;
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    dblptr = (double ***) d->dataptr;
	    if(((*dblptr)[threadid] = (double *) malloc(s * sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < s; j++) (*dblptr)[threadid][j] = 0.;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    stringptr = (char ****) d->dataptr;
	    if(((*stringptr)[threadid] = (char **) malloc(s * sizeof(char *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < s; j++) {
	      if(((*stringptr)[threadid][j] = (char *) malloc(d->maxstringlength)) == NULL)
		error(ERR_MEMALLOC);
	      (*stringptr)[threadid][j][0] = '\0';
	    }
	    break;
	  case VARTOOLS_TYPE_INT:
	    intptr = (int ***) d->dataptr;
	    if(((*intptr)[threadid] = (int *) malloc(s * sizeof(int))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < s; j++) (*intptr)[threadid][j] = 0;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    shortptr = (short ***) d->dataptr;
	    if(((*shortptr)[threadid] = (short *) malloc(s * sizeof(short))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < s; j++) (*shortptr)[threadid][j] = 0;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    floatptr = (float ***) d->dataptr;
	    if(((*floatptr)[threadid] = (float *) malloc(s * sizeof(float))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < s; j++) (*floatptr)[threadid][j] = 0.;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    longptr = (long ***) d->dataptr;
	    if(((*longptr)[threadid] = (long *) malloc(s * sizeof(long))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < s; j++) (*longptr)[threadid][j] = 0;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    charptr = (char ***) d->dataptr;
	    if(((*charptr)[threadid] = (char *) malloc(s * sizeof(char))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < s; j++) (*charptr)[threadid][j] = 0;
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	} else if(Nc > 0) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dbl2ptr = (double ****) d->dataptr;
	    for(j=0; j < Nc; j++) {
	      if(((*dbl2ptr)[threadid][j] = (double *) malloc(s * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0; k < s; k++) (*dbl2ptr)[threadid][j][k] = 0.;
	    }
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    dbl2ptr = (double ****) d->dataptr;
	    for(j=0; j < Nc; j++) {
	      if(((*dbl2ptr)[threadid][j] = (double *) malloc(s * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0; k < s; k++) (*dbl2ptr)[threadid][j][k] = 0.;
	    }
	    break;
	  case VARTOOLS_TYPE_STRING:
	    string2ptr = (char *****) d->dataptr;
	    for(j=0; j < Nc; j++) {
	      if(((*string2ptr)[threadid][j] = (char **) malloc(s * sizeof(char *))) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0; k < s; k++) {
		if(((*string2ptr)[threadid][j][k] = (char *) malloc(d->maxstringlength)) == NULL)
		  error(ERR_MEMALLOC);
		(*string2ptr)[threadid][j][k][0] = '\0';
	      }
	    }
	    break;
	  case VARTOOLS_TYPE_INT:
	    int2ptr = (int ****) d->dataptr;
	    for(j=0; j < Nc; j++) {
	      if(((*int2ptr)[threadid][j] = (int *) malloc(s * sizeof(int))) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0; k < s; k++) (*int2ptr)[threadid][j][k] = 0;
	    }
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    short2ptr = (short ****) d->dataptr;
	    for(j=0; j < Nc; j++) {
	      if(((*short2ptr)[threadid][j] = (short *) malloc(s * sizeof(short))) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0; k < s; k++) (*short2ptr)[threadid][j][k] = 0;
	    }
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    float2ptr = (float ****) d->dataptr;
	    for(j=0; j < Nc; j++) {
	      if(((*float2ptr)[threadid][j] = (float *) malloc(s * sizeof(float))) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0; k < s; k++) (*float2ptr)[threadid][j][k] = 0.;
	    }
	    break;
	  case VARTOOLS_TYPE_LONG:
	    long2ptr = (long ****) d->dataptr;
	    for(j=0; j < Nc; j++) {
	      if(((*long2ptr)[threadid][j] = (long *) malloc(s * sizeof(long))) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0; k < s; k++) (*long2ptr)[threadid][j][k] = 0;
	    }
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    char2ptr = (char ****) d->dataptr;
	    for(j=0; j < Nc; j++) {
	      if(((*char2ptr)[threadid][j] = (char *) malloc(s)) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0; k < s; k++) (*char2ptr)[threadid][j][k] = 0;
	    }
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	} else {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dbl2ptr = (double ****) d->dataptr;
	    if(((*dbl2ptr)[threadid] = (double **) malloc(s * sizeof(double *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < s; j++) {
	      if(((*dbl2ptr)[threadid][j] = (double *) malloc((-Nc) * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0; k < -Nc; k++) (*dbl2ptr)[threadid][j][k] = 0.;
	    }
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    dbl2ptr = (double ****) d->dataptr;
	    if(((*dbl2ptr)[threadid] = (double **) malloc(s * sizeof(double *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < s; j++) {
	      if(((*dbl2ptr)[threadid][j] = (double *) malloc((-Nc) * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0; k < -Nc; k++) (*dbl2ptr)[threadid][j][k] = 0.;
	    }
	    break;
	  case VARTOOLS_TYPE_STRING:
	    string2ptr = (char *****) d->dataptr;
	    if(((*string2ptr)[threadid] = (char ***) malloc(s * sizeof(char **))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < s; j++) {
	      if(((*string2ptr)[threadid][j] = (char **) malloc((-Nc) * sizeof(char *))) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0; k < -Nc; k++) {
		if(((*string2ptr)[threadid][j][k] = (char *) malloc(d->maxstringlength)) == NULL)
		  error(ERR_MEMALLOC);
		(*string2ptr)[threadid][j][k][0] = '\0';
	      }
	    }
	    break;
	  case VARTOOLS_TYPE_INT:
	    int2ptr = (int ****) d->dataptr;
	    if(((*int2ptr)[threadid] = (int **) malloc(s * sizeof(int *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < s; j++) {
	      if(((*int2ptr)[threadid][j] = (int *) malloc((-Nc) * sizeof(int))) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0; k < -Nc; k++) (*int2ptr)[threadid][j][k] = 0;
	    }
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    short2ptr = (short ****) d->dataptr;
	    if(((*short2ptr)[threadid] = (short **) malloc(s * sizeof(short *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < s; j++) {
	      if(((*short2ptr)[threadid][j] = (short *) malloc((-Nc) * sizeof(short))) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0; k < -Nc; k++) (*short2ptr)[threadid][j][k] = 0;
	    }
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    float2ptr = (float ****) d->dataptr;
	    if(((*float2ptr)[threadid] = (float **) malloc(s * sizeof(float *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < s; j++) {
	      if(((*float2ptr)[threadid][j] = (float *) malloc((-Nc) * sizeof(float))) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0; k < -Nc; k++) (*float2ptr)[threadid][j][k] = 0.;
	    }
	    break;
	  case VARTOOLS_TYPE_LONG:
	    long2ptr = (long ****) d->dataptr;
	    if(((*long2ptr)[threadid] = (long **) malloc(s * sizeof(long *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < s; j++) {
	      if(((*long2ptr)[threadid][j] = (long *) malloc((-Nc) * sizeof(long))) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0; k < -Nc; k++) (*long2ptr)[threadid][j][k] = 0;
	    }
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    char2ptr = (char ****) d->dataptr;
	    if(((*char2ptr)[threadid] = (char **) malloc(s * sizeof(char *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < s; j++) {
	      if(((*char2ptr)[threadid][j] = (char *) malloc((-Nc) * sizeof(char))) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0; k < -Nc; k++) (*char2ptr)[threadid][j][k] = 0;
	    }
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      }
    }
  else {
    oldsizesinglelc = p->sizesinglelc[threadid];
    while(p->sizesinglelc[threadid] < Nterm)
      p->sizesinglelc[threadid] *= 2;
    s = p->sizesinglelc[threadid];
    if(p->readimagestring) {
      if((p->stringid_idx[threadid] = (int *) realloc(p->stringid_idx[threadid], s * sizeof(int))) == NULL)
	error(ERR_MEMALLOC);
      for(k=oldsizesinglelc; k < s; k++) p->stringid_idx[threadid][k] = 0;
    }
    for(i=0; i < p->NDataFromLightCurve; i++) {
      d = &(p->DataFromLightCurve[i]);
      Nc = d->Ncolumns;
      if(Nc == 0) {
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = (double ***) d->dataptr;
	  if(((*dblptr)[threadid] = (double *) realloc((*dblptr)[threadid], s * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  for(k=oldsizesinglelc; k < s; k++) (*dblptr)[threadid][k] = 0.;
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  dblptr = (double ***) d->dataptr;
	  if(((*dblptr)[threadid] = (double *) realloc((*dblptr)[threadid], s * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  for(k=oldsizesinglelc; k < s; k++) (*dblptr)[threadid][k] = 0.;
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = (char ****) d->dataptr;
	  if(((*stringptr)[threadid] = (char **) realloc((*stringptr)[threadid], s * sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=oldsizesinglelc; j < s; j++) {
	    if(((*stringptr)[threadid][j] = (char *) malloc(d->maxstringlength)) == NULL)
	      error(ERR_MEMALLOC);
	    (*stringptr)[threadid][j][0] = '\0';
	  }
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = (int ***) d->dataptr;
	  if(((*intptr)[threadid] = (int *) realloc((*intptr)[threadid], s * sizeof(int))) == NULL)
	    error(ERR_MEMALLOC);
	  for(k=oldsizesinglelc; k < s; k++) (*intptr)[threadid][k] = 0;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = (short ***) d->dataptr;
	  if(((*shortptr)[threadid] = (short *) realloc((*shortptr)[threadid], s * sizeof(short))) == NULL)
	    error(ERR_MEMALLOC);
	  for(k=oldsizesinglelc; k < s; k++) (*shortptr)[threadid][k] = 0;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = (float ***) d->dataptr;
	  if(((*floatptr)[threadid] = (float *) realloc((*floatptr)[threadid], s * sizeof(float))) == NULL)
	    error(ERR_MEMALLOC);
	  for(k=oldsizesinglelc; k < s; k++) (*floatptr)[threadid][k] = 0.;
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = (long ***) d->dataptr;
	  if(((*longptr)[threadid] = (long *) realloc((*longptr)[threadid], s * sizeof(long))) == NULL)
	    error(ERR_MEMALLOC);
	  for(k=oldsizesinglelc; k < s; k++) (*longptr)[threadid][k] = 0;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = (char ***) d->dataptr;
	  if(((*charptr)[threadid] = (char *) realloc((*charptr)[threadid], s * sizeof(char))) == NULL)
	    error(ERR_MEMALLOC);
	  for(k=oldsizesinglelc; k < s; k++) (*charptr)[threadid][k] = 0;
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      } else if (Nc > 0) {
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dbl2ptr = (double ****) d->dataptr;
	  for(j=0; j < Nc; j++) {
	    if(((*dbl2ptr)[threadid][j] = (double *) realloc((*dbl2ptr)[threadid][j], s * sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=oldsizesinglelc; k < s; k++) (*dbl2ptr)[threadid][j][k] = 0.;
	  }
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  dbl2ptr = (double ****) d->dataptr;
	  for(j=0; j < Nc; j++) {
	    if(((*dbl2ptr)[threadid][j] = (double *) realloc((*dbl2ptr)[threadid][j], s * sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=oldsizesinglelc; k < s; k++) (*dbl2ptr)[threadid][j][k] = 0.;
	  }
	  break;
	case VARTOOLS_TYPE_STRING:
	  string2ptr = (char *****) d->dataptr;
	  for(j=0; j < Nc; j++) {
	    if(((*string2ptr)[threadid][j] = (char **) realloc((*string2ptr)[threadid][j], s * sizeof(char *))) == NULL)
	      error(ERR_MEMALLOC);
	      for(k=oldsizesinglelc; k < s; k++) {
		if(((*string2ptr)[threadid][j][k] = (char *) malloc(d->maxstringlength)) == NULL)
		  error(ERR_MEMALLOC);
		(*string2ptr)[threadid][j][k][0] = '\0';
	      }
	  }
	  break;
	case VARTOOLS_TYPE_INT:
	  int2ptr = (int ****) d->dataptr;
	  for(j=0; j < Nc; j++) {
	    if(((*int2ptr)[threadid][j] = (int *) realloc((*int2ptr)[threadid][j], s * sizeof(int))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=oldsizesinglelc; k < s; k++) (*int2ptr)[threadid][j][k] = 0;
	  }
	  break;
	case VARTOOLS_TYPE_SHORT:
	  short2ptr = (short ****) d->dataptr;
	  for(j=0; j < Nc; j++) {
	    if(((*short2ptr)[threadid][j] = (short *) realloc((*short2ptr)[threadid][j], s * sizeof(short))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=oldsizesinglelc; k < s; k++) (*short2ptr)[threadid][j][k] = 0;
	  }
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  float2ptr = (float ****) d->dataptr;
	  for(j=0; j < Nc; j++) {
	    if(((*float2ptr)[threadid][j] = (float *) realloc((*float2ptr)[threadid][j], s * sizeof(float))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=oldsizesinglelc; k < s; k++) (*float2ptr)[threadid][j][k] = 0.;
	  }
	  break;
	case VARTOOLS_TYPE_LONG:
	  long2ptr = (long ****) d->dataptr;
	  for(j=0; j < Nc; j++) {
	    if(((*long2ptr)[threadid][j] = (long *) realloc((*long2ptr)[threadid][j], s * sizeof(long))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=oldsizesinglelc; k < s; k++) (*long2ptr)[threadid][j][k] = 0;
	  }
	  break;
	case VARTOOLS_TYPE_CHAR:
	  char2ptr = (char ****) d->dataptr;
	  for(j=0; j < Nc; j++) {
	    if(((*char2ptr)[threadid][j] = (char *) realloc((*char2ptr)[threadid][j], s)) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=oldsizesinglelc; k < s; k++) (*char2ptr)[threadid][j][k] = 0;
	  }
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      } else {
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dbl2ptr = (double ****) d->dataptr;
	  if(((*dbl2ptr)[threadid] = (double **) realloc((*dbl2ptr)[threadid], s * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=oldsizesinglelc; j < s; j++) {
	    if(((*dbl2ptr)[threadid][j] = (double *) malloc((-Nc) * sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=0; k < -Nc; k++) (*dbl2ptr)[threadid][j][k] = 0.;
	  }
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  dbl2ptr = (double ****) d->dataptr;
	  if(((*dbl2ptr)[threadid] = (double **) realloc((*dbl2ptr)[threadid], s * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=oldsizesinglelc; j < s; j++) {
	    if(((*dbl2ptr)[threadid][j] = (double *) malloc((-Nc) * sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=0; k < -Nc; k++) (*dbl2ptr)[threadid][j][k] = 0.;
	  }
	  break;
	case VARTOOLS_TYPE_STRING:
	  string2ptr = (char *****) d->dataptr;
	  if(((*string2ptr)[threadid] = (char ***) realloc((*string2ptr)[threadid], s * sizeof(char **))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=oldsizesinglelc; j < s; j++) {
	    if(((*string2ptr)[threadid][j] = (char **) malloc((-Nc)*sizeof(char *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=0; k < (-Nc); k++) {
	      if(((*string2ptr)[threadid][j][k] = (char *) malloc(d->maxstringlength)) == NULL)
		error(ERR_MEMALLOC);
	      (*string2ptr)[threadid][j][k][0] = '\0';
	    }
	  }
	  break;
	case VARTOOLS_TYPE_INT:
	  int2ptr = (int ****) d->dataptr;
	  if(((*int2ptr)[threadid] = (int **) realloc((*int2ptr)[threadid], s * sizeof(int *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=oldsizesinglelc; j < s; j++) {
	    if(((*int2ptr)[threadid][j] = (int *) malloc((-Nc) * sizeof(int))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=0; k < -Nc; k++) (*int2ptr)[threadid][j][k] = 0;
	  }
	  break;
	case VARTOOLS_TYPE_SHORT:
	  short2ptr = (short ****) d->dataptr;
	  if(((*short2ptr)[threadid] = (short **) realloc((*short2ptr)[threadid], s * sizeof(short *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=oldsizesinglelc; j < s; j++) {
	    if(((*short2ptr)[threadid][j] = (short *) malloc((-Nc) * sizeof(short))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=0; k < -Nc; k++) (*short2ptr)[threadid][j][k] = 0.;
	  }
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  float2ptr = (float ****) d->dataptr;
	  if(((*float2ptr)[threadid] = (float **) realloc((*float2ptr)[threadid], s * sizeof(float *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=oldsizesinglelc; j < s; j++) {
	    if(((*float2ptr)[threadid][j] = (float *) malloc((-Nc) * sizeof(float))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=0; k < -Nc; k++) (*float2ptr)[threadid][j][k] = 0.;
	  }
	  break;
	case VARTOOLS_TYPE_LONG:
	  long2ptr = (long ****) d->dataptr;
	  if(((*long2ptr)[threadid] = (long **) realloc((*long2ptr)[threadid], s * sizeof(long *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=oldsizesinglelc; j < s; j++) {
	    if(((*long2ptr)[threadid][j] = (long *) malloc((-Nc) * sizeof(long))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=0; k < -Nc; k++) (*long2ptr)[threadid][j][k] = 0;
	  }
	  break;
	case VARTOOLS_TYPE_CHAR:
	  char2ptr = (char ****) d->dataptr;
	  if(((*char2ptr)[threadid] = (char **) realloc((*char2ptr)[threadid], s * sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=oldsizesinglelc; j < s; j++) {
	    if(((*char2ptr)[threadid][j] = (char *) malloc((-Nc) * sizeof(char))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=0; k < -Nc; k++) (*char2ptr)[threadid][j][k] = 0;
	  }
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      }
    }
  }
  
}

/* The version of MemAllocDataFromLightCurve below should be called if additional space needs to be allocated for the light curves within a processing thread; This version takes care of re-aligning the pointers to t, mag, err and id. */
void MemAllocDataFromLightCurveMidProcess(ProgramData *p, int threadid, int Nterm) {
  int i;
  _Variable *tvar = NULL;
  _Variable *magvar = NULL;
  _Variable *sigvar = NULL;
  _Variable *stringidvar = NULL;

  for(i=0; i < p->NDefinedVariables; i++) {
    if(p->DefinedVariables[i]->vectortype == VARTOOLS_VECTORTYPE_LC && (
       p->DefinedVariables[i]->datatype == VARTOOLS_TYPE_DOUBLE ||
       p->DefinedVariables[i]->datatype == VARTOOLS_TYPE_CONVERTJD)) {
      if(p->t[threadid] == (*((double ***) p->DefinedVariables[i]->dataptr))[threadid] && tvar == NULL) {
	tvar = p->DefinedVariables[i];
      }
      if(p->mag[threadid] == (*((double ***) p->DefinedVariables[i]->dataptr))[threadid] && magvar == NULL) {
	magvar = p->DefinedVariables[i];
      }
      if(p->sig[threadid] == (*((double ***) p->DefinedVariables[i]->dataptr))[threadid] && sigvar == NULL) {
	sigvar = p->DefinedVariables[i];
      }
    } 
    else if(p->DefinedVariables[i]->vectortype == VARTOOLS_VECTORTYPE_LC &&
	    p->DefinedVariables[i]->datatype == VARTOOLS_TYPE_STRING) {
      if(p->stringid[threadid] == (*((char ****) p->DefinedVariables[i]->dataptr))[threadid] && stringidvar == NULL)
	stringidvar = p->DefinedVariables[i];
    }
  }
  MemAllocDataFromLightCurve(p, threadid, Nterm);
  if(tvar != NULL)
    p->t[threadid] = (*((double ***) tvar->dataptr))[threadid];
  if(magvar != NULL)
    p->mag[threadid] = (*((double ***) magvar->dataptr))[threadid];
  if(sigvar != NULL)
    p->sig[threadid] = (*((double ***) sigvar->dataptr))[threadid];
  if(stringidvar != NULL)
    p->stringid[threadid] = (*((char ****) stringidvar->dataptr))[threadid];
}

void DoChangeVariable(ProgramData *p, _Changevariable *c, int threadid)
{
  int i;
  switch(c->changevar) {
  case VARTOOLS_CHANGEVAR_TIME:
    p->t[threadid] = (*((double ***) c->newvar->dataptr))[threadid];
    sortlcbytime(p->NJD[threadid],p->t[threadid],threadid,p);
    break;
  case VARTOOLS_CHANGEVAR_MAG:
    p->mag[threadid] = (*((double ***) c->newvar->dataptr))[threadid];
    break;
  case VARTOOLS_CHANGEVAR_ERR:
    p->sig[threadid] = (*((double ***) c->newvar->dataptr))[threadid];
    break;
  case VARTOOLS_CHANGEVAR_ID:
    p->stringid[threadid] = (*((char ****) c->newvar->dataptr))[threadid];
    for(i=0;i<p->NJD[threadid];i++)
      p->stringid_idx[threadid][i] = i;
    mysortstringint(p->NJD[threadid], MAXIDSTRINGLENGTH, p->stringid[threadid], p->stringid_idx[threadid]);
    break;
  default:
    error(ERR_CODEERROR);
  }
}

void SetTimeMagSigPointers(ProgramData *p, int threadid)
{
  p->mag[threadid] = (*((double ***) p->magvar->dataptr))[threadid];
  p->sig[threadid] = (*((double ***) p->sigvar->dataptr))[threadid];
  p->t[threadid] = (*((double ***) p->tvar->dataptr))[threadid];
  if(p->readimagestring) {
    p->stringid[threadid] = (*((char ****) p->idvar->dataptr))[threadid];
    mysortstringint(p->NJD[threadid], MAXIDSTRINGLENGTH, p->stringid[threadid], p->stringid_idx[threadid]);
  }
}

void DetermineColumns(ProgramData *p, Command *c)
{
  /* This function registers some of the data which is read in from the light curves, including the time, mag, and err data, and data used by the decorr and other legacy commands from the pre-parselc days. Newer commands should register the data vectors when the command line is parsed.  */
  int Ncol, i, j, k, *incols, sizeincols = 0;
  char format[MAXLEN];
  Ncol = 0;

  _Variable *variable;

  /* Register the vector to store the light curve time data */
  if(!p->inputlcformatused) {
    variable = CreateVariable(p, "t", (char) VARTOOLS_TYPE_DOUBLE, VARTOOLS_VECTORTYPE_LC, NULL);
    p->tvar = variable;

    if(p->inputUTC)
      RegisterDataFromLightCurve(p, variable->dataptr, VARTOOLS_TYPE_CONVERTJD,
				 0, 0, 0, 0, 0, p->UTCformat, variable,
				 (p->coljd == 0 ? -1 : p->coljd), "Time_UTC");
    else {
      if(p->coljd == 0) {
	RegisterDataFromLightCurve(p, variable->dataptr, VARTOOLS_TYPE_DOUBLE,
				   0, 0, 0, 0, 0, "NR", variable,
				   (p->coljd == 0 ? -1 : p->coljd), "Time");
      } else {
	RegisterDataFromLightCurve(p, variable->dataptr, VARTOOLS_TYPE_DOUBLE,
				   0, 0, 0, 0, 0, NULL, variable,
				   (p->coljd), "Time");
      }
    }
    /* Register the vector to store the light curve magnitude data */
    variable = CreateVariable(p, "mag", (char) VARTOOLS_TYPE_DOUBLE, VARTOOLS_VECTORTYPE_LC, NULL);
    p->magvar = variable;
    RegisterDataFromLightCurve(p, variable->dataptr, VARTOOLS_TYPE_DOUBLE,
			       0, 0, 0, 0, 0, NULL, variable,
			       (p->colmag == 0 ? -1 : p->colmag), "Mag");
    
    /* Register the vector to store the light curve magnitude uncertainty data */
    variable = CreateVariable(p, "err", (char) VARTOOLS_TYPE_DOUBLE, VARTOOLS_VECTORTYPE_LC, NULL);
    p->sigvar = variable;
    if(p->colsig == 0) {
      RegisterDataFromLightCurve(p, variable->dataptr, VARTOOLS_TYPE_DOUBLE,
				 0, 0, 0, 0, 0, "1", variable,
				 (p->colsig == 0 ? -1 : p->colsig), "errMag");
    } else {
      RegisterDataFromLightCurve(p, variable->dataptr, VARTOOLS_TYPE_DOUBLE,
				 0, 0, 0, 0, 0, NULL, variable,
				 p->colsig, "errMag");
    }
    
    /* Register the vector to store the light curve stringids if used */
    if(p->readimagestring)
      {
	variable = CreateVariable(p, "id", (char) VARTOOLS_TYPE_STRING, VARTOOLS_VECTORTYPE_LC, NULL);
	p->idvar = variable;
	RegisterDataFromLightCurve(p, variable->dataptr, VARTOOLS_TYPE_STRING,
				   MAXIDSTRINGLENGTH, 0, 0, 0, 0, NULL, NULL, p->colstringid,
				   "ID_String");
      }
  } else {
    if(p->coljd < 0) {
      variable = CreateVariable(p, "t", (char) VARTOOLS_TYPE_DOUBLE, VARTOOLS_VECTORTYPE_LC, NULL);
      p->tvar = variable;
      RegisterDataFromLightCurve(p, variable->dataptr, VARTOOLS_TYPE_DOUBLE,
				 0, 0, 0, 0, 0, "NR", variable, -1, "Time");
    }
    if(p->colmag < 0) {
      variable = CreateVariable(p, "mag", (char) VARTOOLS_TYPE_DOUBLE, VARTOOLS_VECTORTYPE_LC, NULL);
      p->magvar = variable;
      RegisterDataFromLightCurve(p, variable->dataptr, VARTOOLS_TYPE_DOUBLE,
				 0, 0, 0, 0, 0, NULL, variable, -1, "Mag");
    }
    if(p->colsig < 0) {
      variable = CreateVariable(p, "err", (char) VARTOOLS_TYPE_DOUBLE, VARTOOLS_VECTORTYPE_LC, NULL);
      p->sigvar = variable;
      RegisterDataFromLightCurve(p, variable->dataptr, VARTOOLS_TYPE_DOUBLE,
				 0, 0, 0, 0, 0, "1", variable, -1, "errMag");
    }
  }

  /* Register column vectors to use with decorr or tfa_sr commands */
  for(i=0;i < p->Ncommands;i++)
    {
      if(c[i].cnum == CNUM_DECORR)
	{
	  if(c[i].Decorr->N_lcterms > 0) {
	    for(k=0;k<c[i].Decorr->N_lcterms;k++)
	      {
		RegisterDataFromLightCurve(p, 
					   (void *) &(c[i].Decorr->lcdecorr_terms_in[k]),
					   VARTOOLS_TYPE_DOUBLE,
					   0, 0, i, 0, 0, NULL, NULL,
					   c[i].Decorr->lc_columns[k],
					   "LC_Decorr_Incolumn");
	      }
	  }
	  if(c[i].Decorr->N_decorrterms > 0) {
	    RegisterDataFromLightCurve(p,
				       (void *) &(c[i].Decorr->decorr_terms),
				       VARTOOLS_TYPE_DOUBLE,
				       0, -c[i].Decorr->N_decorrterms, i, 0, 0,
				       NULL, NULL, -1, "DecorrData");
	  }
	}
      else if(c[i].cnum == CNUM_TFA_SR ? c[i].TFA_SR->decorrflag : 0)
	{
	  for(k=0;k<c[i].TFA_SR->decorr_Nlcterms;k++)
	    {
	      RegisterDataFromLightCurve(p,
					 (void *) &(c[i].TFA_SR->lcdecorr_terms_in[k]),
					 VARTOOLS_TYPE_DOUBLE,
					 0, 0, i, 0, 0, NULL, NULL,
					 c[i].TFA_SR->decorr_lc_columns[k],
					 "TFA_SR_LC_Decorr_Incolumn");
	    }
	}
    }
}
	  
#ifdef USECFITSIO

/* Find the column numbers from the binary lc header keywords */
void get_fitslc_header_columns(fitsfile *infile,
			       ProgramData *p)
{
  _DataFromLightCurve *d;
  int j, Nc, u, Ncabs, i, status = 0;
  /* Check if the program needs to do this at all */
  if(p->lc_getcolumnsfromheader) {
    /* Check if this has not yet been done. Note that if parallelization
       is enabled, we will need to use locks to ensure this is only done
       once */
    if(p->lc_getcolumnsfromheader_notyetset) {
#ifdef PARALLEL
      if(p->Nproc_allow > 1) {
	while(pthread_mutex_trylock(&(p->lc_getcolumnsfromheader_mutex)));
	if(!p->lc_getcolumnsfromheader_notyetset) {
	  pthread_mutex_unlock(&(p->lc_getcolumnsfromheader_mutex));
	  return;
	}
      } else {
	if(!p->lc_getcolumnsfromheader_notyetset) {
	  return;
	}
      }
#endif
      for(j=0; j < p->NDataFromLightCurve; j++) {
	d = &(p->DataFromLightCurve[j]);
	if(d->incolumn_header_names == NULL)
	  continue;
	Nc = d->Ncolumns;
	Ncabs = Nc < 0 ? -Nc : Nc;
	if(Nc == 0) Ncabs += 1;
	for(u = 0; u < Ncabs; u++) {
	  if(d->incolumn_header_names[u] != NULL) {
	    fits_get_colnum(infile, 0, d->incolumn_header_names[0],
			    &(d->incolumns[u]), &status);
	    if(status == COL_NOT_FOUND)
	      error2(ERR_MISSING_FITSLC_HEADERNAME,
		     d->incolumn_header_names[u]);
	    if(d->variable != NULL) {
	      if(d->variable == p->tvar)
		p->coljd = d->incolumns[u];
	      else if(d->variable == p->magvar)
		p->colmag = d->incolumns[u];
	      else if(d->variable == p->sigvar)
		p->colsig = d->incolumns[u];
	    }
	  }
	}
      }
      p->lc_getcolumnsfromheader_notyetset = 0;
#ifdef PARALLEL
      if(p->Nproc_allow > 1) {
	pthread_mutex_unlock(&(p->lc_getcolumnsfromheader_mutex));
      }
#endif
    }
  }
}

#endif

FILE *ExecLCOpenCommand(ProgramData *p, Command *c, int lc, int lc2)
/* This function uses popen to carry out a requested shell command in
   opening a light curve */
{
  char *execcommand;
  int size_execcommand = 2048;
  int i, i2;
  int lc_in_name_length;
  FILE *return_pipe;
  lc_in_name_length = strlen(p->lcnames[lc2]);
  if((execcommand = (char *) malloc((size_execcommand+1))) == NULL)
    error(ERR_MEMALLOC);
  i = 0; i2 = 0;

  /* Parse the command string, substituting the input file name as needed */
  while(p->lc_open_exec_command_str[i2] != '\0')
    {
      if(p->lc_open_exec_command_str[i2] != '%')
	{
	  if(i >= (size_execcommand)) {
	    size_execcommand *= 2;
	    if((execcommand = (char *) realloc(execcommand, (size_execcommand+1))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  execcommand[i] = p->lc_open_exec_command_str[i2];
	  i++;
	  execcommand[i] = '\0';
	  i2++;
	}
      else
	{
	  i2++;
	  if(p->lc_open_exec_command_str[i2] == 's')
	    {
	      i2++;
	      if((i + lc_in_name_length) >= (size_execcommand)) {
		while((i + lc_in_name_length) >= (size_execcommand))
		  size_execcommand *= 2;
		if((execcommand = (char *) realloc(execcommand, (size_execcommand+1))) == NULL)
		  error(ERR_MEMALLOC);
	      }
	      sprintf(&execcommand[i],"%s",p->lcnames[lc]);
	      i = strlen(execcommand);
	    }
	  else if(p->lc_open_exec_command_str[i2] == '%')
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
    if(p->skipmissing) {
      error2_noexit(ERR_CANNOTOPEN,execcommand);
    }
    else
      error2(ERR_CANNOTOPEN,execcommand);
  }
  return return_pipe;
}

#ifdef USECFITSIO

int ReadFitsLightCurve(ProgramData *p, Command *c, int lc, int lc2)
{
  fitsfile *infile;
  int status, hdunum, ncols, hdutype, anynulallcolumns;
  long nrows;
  int j, jold, k, i, l, N, oldsizesinglelc, anynul;
  char *nullarray, *nullarraystore;
  char **tmpstring;
  char tmpchar;

  _DataFromLightCurve *d;
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

  int Nc, u;

  p->NJD[lc2] = 0;

#ifdef PARALLEL
  if(p->Nproc_allow > 1) {
    while(pthread_mutex_trylock(&(p->cfitsio_mutex)));
  }
#endif

  tmpstring = (char **) malloc(sizeof(char *));
  tmpstring[0] = (char *) malloc(MAXLEN);

  hdutype = 0; status = 0; nrows = 0; ncols = 0;

  if((fits_open_file(&infile,p->lcnames[lc],READONLY,&status)))
    {
      if(p->skipmissing) {
	error2_noexit(ERR_CANNOTOPEN,p->lcnames[lc]);
	return(ERR_CANNOTOPEN);
      }
      else
	error2(ERR_CANNOTOPEN,p->lcnames[lc]);
    }
  if(fits_get_hdu_num(infile, &hdunum) == 1)
    {
      fits_movabs_hdu(infile, 2, &hdutype, &status);
    }
  else
    fits_get_hdu_type(infile, &hdutype, &status);
  
  if(hdutype == IMAGE_HDU) {
    if(p->skipmissing) {
      error2_noexit(ERR_IMAGEHDU,p->lcnames[lc]);
      return(ERR_IMAGEHDU);
    }
    else
      error2(ERR_IMAGEHDU,p->lcnames[lc]);
  }

  fits_get_num_rows(infile, &nrows, &status);
  fits_get_num_cols(infile, &ncols, &status);

  if((nullarray = (char *) calloc(nrows, sizeof(char))) == NULL)
    error(ERR_MEMALLOC);
  if((nullarraystore = (char *) calloc(nrows, sizeof(char))) == NULL)
    error(ERR_MEMALLOC);

  if(p->lc_getcolumnsfromheader) {
    if(p->lc_getcolumnsfromheader_notyetset) {
      get_fitslc_header_columns(infile,p);
    }
  }

  /* Increase the memory for the data vectors if needed */
  MemAllocDataFromLightCurve(p, lc2, nrows);

  anynulallcolumns = 0;

  /* Read in the Columns */
  for(j=0; j < p->NDataFromLightCurve; j++) {
    d = &(p->DataFromLightCurve[j]);
    if(d->incolumns[0] <= 0)
      continue;
    Nc = d->Ncolumns;
    status = 0;
    anynul = 0;
    if(Nc == 0) {
      switch(d->datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dblptr = (double ***) d->dataptr;
	fits_read_colnull(infile, TDOUBLE, d->incolumns[0], 1, 1, nrows, &((*dblptr)[lc2][0]), nullarray, &anynul,&status);
	break;
      case VARTOOLS_TYPE_CONVERTJD:
	dblptr = (double ***) d->dataptr;
	for(k=1; k<=nrows && !status; k++) {
	  fits_read_col_str(infile, d->incolumns[0], k, 1, 1, 0, tmpstring, &anynul, &status);
	  if(anynul)
	    {
	      nullarray[k-1] = 1;
	      anynulallcolumns = 1;
	      anynul = 0;
	    } 
	  else 
	    {
	      convertUTCtoJD(tmpstring[0],d->scanformat,d->UTCindex,(double *) (&((*dblptr)[lc2][k-1])));
	    }
	}
	break;
      case VARTOOLS_TYPE_STRING:
	stringptr = (char ****) d->dataptr;
	for(k=1; k<=nrows && !status; k++) {
	  fits_read_col_str(infile, d->incolumns[0], k, 1, 1, 0, &(((*stringptr)[lc2][k-1])), &anynul, &status);
	  if(anynul)
	    {
	      nullarray[k-1] = 1;
	      anynulallcolumns = 1;
	      anynul = 0;
	    }
	}
	break;
      case VARTOOLS_TYPE_INT:
	intptr = (int ***) d->dataptr;
	fits_read_colnull(infile, TINT, d->incolumns[0], 1, 1, nrows, &((*intptr)[lc2][0]), nullarray, &anynul,&status);
	break;
      case VARTOOLS_TYPE_FLOAT:
	floatptr = (float ***) d->dataptr;
	fits_read_colnull(infile, TFLOAT, d->incolumns[0], 1, 1, nrows, &((*floatptr)[lc2][0]), nullarray, &anynul,&status);
	break;
      case VARTOOLS_TYPE_LONG:
	longptr = (long ***) d->dataptr;
	fits_read_colnull(infile, TLONG, d->incolumns[0], 1, 1, nrows, &((*longptr)[lc2][0]), nullarray, &anynul,&status);
	break;
      case VARTOOLS_TYPE_SHORT:
	shortptr = (short ***) d->dataptr;
	fits_read_colnull(infile, TSHORT, d->incolumns[0], 1, 1, nrows, &((*shortptr)[lc2][0]), nullarray, &anynul,&status);
	break;
      case VARTOOLS_TYPE_CHAR:
	charptr = (char ***) d->dataptr;
	fits_read_colnull(infile, TBYTE, d->incolumns[0], 1, 1, nrows, &((*charptr)[lc2][0]), nullarray, &anynul,&status);
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
    } else if (Nc > 0) {
      for(u=0; u < Nc; u++) {
	anynul = 0;
	status = 0;
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dbl2ptr = (double ****) d->dataptr;
	  fits_read_colnull(infile, TDOUBLE, d->incolumns[u], 1, 1, nrows, &((*dbl2ptr)[lc2][u][0]), nullarray, &anynul,&status);
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  dbl2ptr = (double ****) d->dataptr;
	  for(k=1; k<=nrows && !status; k++) {
	    fits_read_col_str(infile, d->incolumns[u], k, 1, 1, 0, tmpstring, &anynul, &status);
	    if(anynul)
	      {
		nullarray[k-1] = 1;
		anynulallcolumns = 1;
		anynul = 0;
	      } 
	    else 
	      {
		convertUTCtoJD(tmpstring[0],d->scanformat,d->UTCindex,(double *) (&((*dbl2ptr)[lc2][u][k-1])));
	      }
	  }
	  break;
	case VARTOOLS_TYPE_STRING:
	  string2ptr = (char *****) d->dataptr;
	  for(k=1; k<=nrows && !status; k++) {
	    fits_read_col_str(infile, d->incolumns[u], k, 1, 1, 0, &(((*string2ptr)[lc2][u][k-1])), &anynul, &status);
	    if(anynul)
	      {
		nullarray[k-1] = 1;
		anynulallcolumns = 1;
		anynul = 0;
	      }
	  }
	  break;
	case VARTOOLS_TYPE_INT:
	  int2ptr = (int ****) d->dataptr;
	  fits_read_colnull(infile, TINT, d->incolumns[u], 1, 1, nrows, &((*int2ptr)[lc2][u][0]), nullarray, &anynul,&status);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  float2ptr = (float ****) d->dataptr;
	  fits_read_colnull(infile, TFLOAT, d->incolumns[u], 1, 1, nrows, &((*float2ptr)[lc2][u][0]), nullarray, &anynul,&status);
	  break;
	case VARTOOLS_TYPE_LONG:
	  long2ptr = (long ****) d->dataptr;
	  fits_read_colnull(infile, TLONG, d->incolumns[u], 1, 1, nrows, &((*long2ptr)[lc2][u][0]), nullarray, &anynul,&status);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  short2ptr = (short ****) d->dataptr;
	  fits_read_colnull(infile, TSHORT, d->incolumns[u], 1, 1, nrows, &((*short2ptr)[lc2][u][0]), nullarray, &anynul,&status);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  char2ptr = (char ****) d->dataptr;
	  fits_read_colnull(infile, TBYTE, d->incolumns[u], 1, 1, nrows, &((*char2ptr)[lc2][u][0]), nullarray, &anynul,&status);
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
      }
    } else {
      for(u=0; u < (-Nc); u++) {
	anynul = 0;
	status = 0;
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dbl2ptr = (double ****) d->dataptr;
	  for(k=1; k<=nrows && !status; k++) {
	    fits_read_colnull(infile, TDOUBLE, d->incolumns[u], k, 1, 1, &((*dbl2ptr)[lc2][k-1][u]), &tmpchar, &anynul,&status);
	    if(anynul)
	      {
		nullarraystore[k-1] = 1;
		anynulallcolumns = 1;
		anynul = 0;
	      }
	  }
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  dbl2ptr = (double ****) d->dataptr;
	  for(k=1; k<=nrows && !status; k++) {
	    fits_read_col_str(infile, d->incolumns[u], k, 1, 1, 0, tmpstring, &anynul, &status);
	    if(anynul)
	      {
		nullarraystore[k-1] = 1;
		anynulallcolumns = 1;
		anynul = 0;
	      } 
	    else 
	      {
		convertUTCtoJD(tmpstring[0],d->scanformat,d->UTCindex,(double *) (&((*dbl2ptr)[lc2][k-1][u])));
	      }
	  }
	  break;
	case VARTOOLS_TYPE_STRING:
	  string2ptr = (char *****) d->dataptr;
	  for(k=1; k<=nrows && !status; k++) {
	    fits_read_col_str(infile, d->incolumns[u], k, 1, 1, 0, &(((*string2ptr)[lc2][k-1][u])), &anynul, &status);
	    if(anynul)
	      {
		nullarraystore[k-1] = 1;
		anynulallcolumns = 1;
		anynul = 0;
	      }
	  }
	  break;
	case VARTOOLS_TYPE_INT:
	  int2ptr = (int ****) d->dataptr;
	  for(k=1; k<=nrows && !status; k++) {
	    fits_read_colnull(infile, TINT, d->incolumns[u], k, 1, 1, &((*int2ptr)[lc2][k-1][u]), &tmpchar, &anynul,&status);
	    if(anynul)
	      {
		nullarraystore[k-1] = 1;
		anynulallcolumns = 1;
		anynul = 0;
	      }
	  }
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  float2ptr = (float ****) d->dataptr;
	  for(k=1; k<=nrows && !status; k++) {
	    fits_read_colnull(infile, TFLOAT, d->incolumns[u], k, 1, 1, &((*float2ptr)[lc2][k-1][u]), &tmpchar, &anynul,&status);
	    if(anynul)
	      {
		nullarraystore[k-1] = 1;
		anynulallcolumns = 1;
		anynul = 0;
	      }
	  }
	  break;
	case VARTOOLS_TYPE_LONG:
	  long2ptr = (long ****) d->dataptr;
	  for(k=1; k<=nrows && !status; k++) {
	    fits_read_colnull(infile, TLONG, d->incolumns[u], k, 1, 1, &((*long2ptr)[lc2][k-1][u]), &tmpchar, &anynul,&status);
	    if(anynul)
	      {
		nullarraystore[k-1] = 1;
		anynulallcolumns = 1;
		anynul = 0;
	      }
	  }
	  break;
	case VARTOOLS_TYPE_SHORT:
	  short2ptr = (short ****) d->dataptr;
	  for(k=1; k<=nrows && !status; k++) {
	    fits_read_colnull(infile, TSHORT, d->incolumns[u], k, 1, 1, &((*short2ptr)[lc2][k-1][u]), &tmpchar, &anynul,&status);
	    if(anynul)
	      {
		nullarraystore[k-1] = 1;
		anynulallcolumns = 1;
		anynul = 0;
	      }
	  }
	  break;
	case VARTOOLS_TYPE_CHAR:
	  char2ptr = (char ****) d->dataptr;
	  for(k=1; k<=nrows && !status; k++) {
	    fits_read_colnull(infile, TBYTE, d->incolumns[u], k, 1, 1, &((*char2ptr)[lc2][k-1][u]), &tmpchar, &anynul,&status);
	    if(anynul)
	      {
		nullarraystore[k-1] = 1;
		anynulallcolumns = 1;
		anynul = 0;
	      }
	  }
	  break;
	default:
	  error(ERR_BADTYPE);
	}
	if(anynul)
	  anynulallcolumns = 1;
	if(status) {
	  fits_report_error(stderr, status);
	  error(ERR_FITSERROR);
	}
      }
    }
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
	    for(k=0; k < p->NDataFromLightCurve; k++) {
	      d = &(p->DataFromLightCurve[k]);
	      if(d->incolumns[0] <= 0)
		continue;
	      Nc = d->Ncolumns;
	      if(Nc == 0) {
		switch(d->datatype) {
		case VARTOOLS_TYPE_DOUBLE:
		  dblptr = (double ***) d->dataptr;
		  (*dblptr)[lc2][j] = (*dblptr)[lc2][i];
		  break;
		case VARTOOLS_TYPE_CONVERTJD:
		  dblptr = (double ***) d->dataptr;
		  (*dblptr)[lc2][j] = (*dblptr)[lc2][i];
		  break;
		case VARTOOLS_TYPE_STRING:
		  stringptr = (char ****) d->dataptr;
		  sprintf((*stringptr)[lc2][j],"%s",(*stringptr)[lc2][i]);
		  break;
		case VARTOOLS_TYPE_INT:
		  intptr = (int ***) d->dataptr;
		  (*intptr)[lc2][j] = (*intptr)[lc2][i];
		  break;
		case VARTOOLS_TYPE_FLOAT:
		  floatptr = (float ***) d->dataptr;
		  (*floatptr)[lc2][j] = (*floatptr)[lc2][i];
		  break;
		case VARTOOLS_TYPE_LONG:
		  longptr = (long ***) d->dataptr;
		  (*longptr)[lc2][j] = (*longptr)[lc2][i];
		  break;
		case VARTOOLS_TYPE_SHORT:
		  shortptr = (short ***) d->dataptr;
		  (*shortptr)[lc2][j] = (*shortptr)[lc2][i];
		  break;
		case VARTOOLS_TYPE_CHAR:
		  charptr = (char ***) d->dataptr;
		  (*charptr)[lc2][j] = (*charptr)[lc2][i];
		  break;
		default:
		  error(ERR_BADTYPE);
		}
	      }
	      else if(Nc > 0) {
		for(u=0; u < Nc; u++) {
		  switch(d->datatype) {
		  case VARTOOLS_TYPE_DOUBLE:
		    dbl2ptr = (double ****) d->dataptr;
		    (*dbl2ptr)[lc2][u][j] = (*dbl2ptr)[lc2][u][i];
		    break;
		  case VARTOOLS_TYPE_CONVERTJD:
		    dbl2ptr = (double ****) d->dataptr;
		    (*dbl2ptr)[lc2][u][j] = (*dbl2ptr)[lc2][u][i];
		    break;
		  case VARTOOLS_TYPE_STRING:
		    string2ptr = (char *****) d->dataptr;
		    sprintf((*string2ptr)[lc2][u][j],"%s",(*string2ptr)[lc2][u][i]);
		    break;
		  case VARTOOLS_TYPE_INT:
		    int2ptr = (int ****) d->dataptr;
		    (*int2ptr)[lc2][u][j] = (*int2ptr)[lc2][u][i];
		    break;
		  case VARTOOLS_TYPE_FLOAT:
		    float2ptr = (float ****) d->dataptr;
		    (*float2ptr)[lc2][u][j] = (*float2ptr)[lc2][u][i];
		    break;
		  case VARTOOLS_TYPE_LONG:
		    long2ptr = (long ****) d->dataptr;
		    (*long2ptr)[lc2][u][j] = (*long2ptr)[lc2][u][i];
		    break;
		  case VARTOOLS_TYPE_SHORT:
		    short2ptr = (short ****) d->dataptr;
		    (*short2ptr)[lc2][u][j] = (*short2ptr)[lc2][u][i];
		    break;
		  case VARTOOLS_TYPE_CHAR:
		    char2ptr = (char ****) d->dataptr;
		    (*char2ptr)[lc2][u][j] = (*char2ptr)[lc2][u][i];
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
		    (*dbl2ptr)[lc2][j][u] = (*dbl2ptr)[lc2][i][u];
		    break;
		  case VARTOOLS_TYPE_CONVERTJD:
		    dbl2ptr = (double ****) d->dataptr;
		    (*dbl2ptr)[lc2][j][u] = (*dbl2ptr)[lc2][i][u];
		    break;
		  case VARTOOLS_TYPE_STRING:
		    string2ptr = (char *****) d->dataptr;
		    sprintf((*string2ptr)[lc2][j][u],"%s",(*string2ptr)[lc2][i][u]);
		    break;
		  case VARTOOLS_TYPE_INT:
		    int2ptr = (int ****) d->dataptr;
		    (*int2ptr)[lc2][j][u] = (*int2ptr)[lc2][i][u];
		    break;
		  case VARTOOLS_TYPE_FLOAT:
		    float2ptr = (float ****) d->dataptr;
		    (*float2ptr)[lc2][j][u] = (*float2ptr)[lc2][i][u];
		    break;
		  case VARTOOLS_TYPE_LONG:
		    long2ptr = (long ****) d->dataptr;
		    (*long2ptr)[lc2][j][u] = (*long2ptr)[lc2][i][u];
		    break;
		  case VARTOOLS_TYPE_SHORT:
		    short2ptr = (short ****) d->dataptr;
		    (*short2ptr)[lc2][j][u] = (*short2ptr)[lc2][i][u];
		    break;
		  case VARTOOLS_TYPE_CHAR:
		    char2ptr = (char ****) d->dataptr;
		    (*char2ptr)[lc2][j][u] = (*char2ptr)[lc2][i][u];
		    break;
		  default:
		  error(ERR_BADTYPE);
		  }
		}
	      }
	    }
	  }
	  j++;
	}
    }
    p->NJD[lc2] = j;
  }
  else
    p->NJD[lc2] = nrows;
  if(p->readimagestring)
    {
      for(i=0;i<p->NJD[lc2];i++)
	p->stringid_idx[lc2][i] = i;
    }

  /* Fill in any variables calculated from analytic expressions */
  for(j=0; j < p->NDataFromLightCurve; j++) {
    d = &(p->DataFromLightCurve[j]);
    if(d->incolumns[0] <= 0)
      {
	if(d->scanformat != NULL) {
	  Nc = d->Ncolumns;
	  if(Nc != 0)
	    error(ERR_CODEERROR);
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    for(i=0; i < p->NJD[lc2]; i++) {
	      dblptr = (double ***) d->dataptr;
	      (*dblptr)[lc2][i] = 
		EvaluateExpression(lc, lc2, i, d->expression);
	    }
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    for(i=0; i < p->NJD[lc2]; i++) {
	      floatptr = (float ***) d->dataptr;
	      (*floatptr)[lc2][i] = 
		(float) EvaluateExpression(lc, lc2, i, d->expression);
	    }
	    break;
	  case VARTOOLS_TYPE_INT:
	    for(i=0; i < p->NJD[lc2]; i++) {
	      intptr = (int ***) d->dataptr;
	      (*intptr)[lc2][i] = 
		(int) EvaluateExpression(lc, lc2, i, d->expression);
	    }
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    for(i=0; i < p->NJD[lc2]; i++) {
	      shortptr = (short ***) d->dataptr;
	      (*shortptr)[lc2][i] = 
		(short) EvaluateExpression(lc, lc2, i, d->expression);
	    }
	    break;
	  case VARTOOLS_TYPE_LONG:
	    for(i=0; i < p->NJD[lc2]; i++) {
	      longptr = (long ***) d->dataptr;
	      (*longptr)[lc2][i] = 
		(long) EvaluateExpression(lc, lc2, i, d->expression);
	    }
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
	else {
	  Nc = d->Ncolumns;
	  if(Nc != 0)
	    error(ERR_CODEERROR);
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    for(i=0; i < p->NJD[lc2]; i++) {
	      dblptr = (double ***) d->dataptr;
	      (*dblptr)[lc2][i] = 0.;
	    }
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    for(i=0; i < p->NJD[lc2]; i++) {
	      floatptr = (float ***) d->dataptr;
	      (*floatptr)[lc2][i] = 0.;
	    }
	    break;
	  case VARTOOLS_TYPE_INT:
	    for(i=0; i < p->NJD[lc2]; i++) {
	      intptr = (int ***) d->dataptr;
	      (*intptr)[lc2][i] = 0;
	    }
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    for(i=0; i < p->NJD[lc2]; i++) {
	      shortptr = (short ***) d->dataptr;
	      (*shortptr)[lc2][i] = 0;
	    }
	    break;
	  case VARTOOLS_TYPE_LONG:
	    for(i=0; i < p->NJD[lc2]; i++) {
	      longptr = (long ***) d->dataptr;
	      (*longptr)[lc2][i] = 0;
	    }
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      }
  }
  
  /* Fill the decorr matrix if we are decorrelating */
  if(p->decorrflag)
    Filldecorr_matrix(p,c,lc2);
  
  free(nullarray);
  free(nullarraystore);
  free(tmpstring[0]);
  free(tmpstring);
  return 0;
} 
#endif

#ifdef _USEBINARY_LC

/* Find the column numbers from the binary lc header keywords */
void get_binarylc_header_columns(int num_columns, 
				 char **lc_column_names,
				 ProgramData *p)
{
  _DataFromLightCurve *d;
  int j, Nc, u, Ncabs, i;
  /* Check if the program needs to do this at all */
  if(p->lc_getcolumnsfromheader) {
    /* Check if this has not yet been done. Note that if parallelization
       is enabled, we will need to use locks to ensure this is only done
       once */
    if(p->lc_getcolumnsfromheader_notyetset) {
#ifdef PARALLEL
      if(p->Nproc_allow > 1) {
	while(pthread_mutex_trylock(&(p->lc_getcolumnsfromheader_mutex)));
	if(!p->lc_getcolumnsfromheader_notyetset) {
	  pthread_mutex_unlock(&(p->lc_getcolumnsfromheader_mutex));
	  return;
	}
      } else {
	if(!p->lc_getcolumnsfromheader_notyetset) {
	  return;
	}
      }
#endif
      for(j=0; j < p->NDataFromLightCurve; j++) {
	d = &(p->DataFromLightCurve[j]);
	if(d->incolumn_header_names == NULL)
	  continue;
	Nc = d->Ncolumns;
	Ncabs = Nc < 0 ? -Nc : Nc;
	if(Nc == 0) Ncabs += 1;
	for(u = 0; u < Ncabs; u++) {
	  if(d->incolumn_header_names[u] != NULL) {
	    for(i=0; i < num_columns; i++) {
	      if(!strcmp(lc_column_names[i], d->incolumn_header_names[u])) {
		d->incolumns[u] = i+1;
		if(d->variable != NULL) {
		  if(d->variable == p->tvar)
		    p->coljd = d->incolumns[u];
		  else if(d->variable == p->magvar)
		    p->colmag = d->incolumns[u];
		  else if(d->variable == p->sigvar)
		    p->colsig = d->incolumns[u];
		}
		break;
	      }
	    }
	    if(i >= num_columns) {
	      error2(ERR_MISSING_BINARYLC_HEADERNAME,d->incolumn_header_names[u]);
	    }
	  }
	}
      }
      p->lc_getcolumnsfromheader_notyetset = 0;
#ifdef PARALLEL
      if(p->Nproc_allow > 1) {
	pthread_mutex_unlock(&(p->lc_getcolumnsfromheader_mutex));
      }
#endif
    }
  }
}


/* stopped working here.  need to implement the binary lc i/o option in 
   parsecommandline.c. Also need to allow output of binary lcs. */

/* Reads the relevant portions of the header in order to fill in the part of
 * the lightcurve that depends only on the header and skips to the binary
 * part. Also sets the light curve format. If lc==NULL simply skips over the
 * header part of the file. */
int read_binarylightcurve_header(FILE *f, int *num_apertures, int *num_columns,
			   int *num_pts, int *hdr_size,
			   int *lc_record_size, char ****lc_header,
			   int *memsize_lc_header,
			   char ***lc_column_names,
			   BinLC_OutputFormat **lc_column_format,
			   int *memsize_lc_columns)
{
  int hdr_rec_size = MAXLEN;
  char *keyword, *value, *comment;
  char missing_apertures=0;
  int line, kwsize, col_ind=0, indx;
  int min_headersize=100;

  *num_pts = -1; *hdr_size = -1; *num_columns = -1; *num_apertures = -1;
  *lc_record_size = 0;

  keyword = malloc(MAXLEN);

  if(!(*memsize_lc_header)) {
    *memsize_lc_header = min_headersize;
    if(((*lc_header) = (char ***) malloc((*memsize_lc_header)*sizeof(char **))) == NULL)
      error(ERR_MEMALLOC);
    for(indx = 0; indx < (*memsize_lc_header); indx++) {
      if(((*lc_header)[indx] = (char **) malloc(3*sizeof(char *))) == NULL)
	error(ERR_MEMALLOC);
      if(((*lc_header)[indx][0] = (char *) malloc(MAXLEN)) == NULL ||
	 ((*lc_header)[indx][1] = (char *) malloc(MAXLEN)) == NULL ||
	 ((*lc_header)[indx][2] = (char *) malloc(MAXLEN)) == NULL)
	error(ERR_MEMALLOC);
    }
  }

  for(line=0; (*hdr_size)<0 || line<(*hdr_size) ; line++) {
    if(!fgets(keyword, hdr_rec_size, f)) {
      free(keyword); 
      return -1;
    }
    bin_lc_parse_header_record(keyword, hdr_rec_size, 
			       &value, &comment, 
			       &kwsize);
    if((*num_pts)<0 && strcmp(keyword, "NUMPOINTS")==0) {
      (*num_pts)=atol(value);
    } 
    else if((*num_columns)<0 && 
	    strcmp(keyword, "NUMCOLUMNS")==0) {
      (*num_columns)=atoi(value);
      if((*num_columns) > (*memsize_lc_columns)) {
	if(!(*memsize_lc_columns)) {
	  if(((*lc_column_names) = (char **) malloc((*num_columns)*sizeof(char *))) == NULL ||
	     ((*lc_column_format) = (BinLC_OutputFormat *) malloc((*num_columns)*sizeof(BinLC_OutputFormat))) == NULL)
	    error(ERR_MEMALLOC);
	  for(indx=0; indx < (*num_columns); indx++) {
	    if(((*lc_column_names)[indx] = (char *) malloc(MAXIDSTRINGLENGTH)) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  (*memsize_lc_columns) = (*num_columns);
	}
	else {
	  if(((*lc_column_names) = (char **) realloc((*lc_column_names), (*num_columns)*sizeof(char *))) == NULL ||
	     ((*lc_column_format) = (BinLC_OutputFormat *) realloc((*lc_column_names), (*num_columns)*sizeof(BinLC_OutputFormat))) == NULL)
	    error(ERR_MEMALLOC);
	  for(indx=(*memsize_lc_columns); indx < (*num_columns); indx++) {
	    if(((*lc_column_names)[indx] = (char *) malloc(MAXIDSTRINGLENGTH)) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  (*memsize_lc_columns) = (*num_columns);
	}
      }
    } else if((*num_apertures)<0 && 
	      strcmp(keyword, "NUMAPERTURES")==0) {
      (*num_apertures)=atoi(value);
      missing_apertures=(*num_apertures);
    } else if((*hdr_size)<0 && strcmp(keyword, "HEADERSIZE")==0) {
      (*hdr_size)=atoi(value);
      if((*hdr_size) > 0) {
	if(!(*memsize_lc_header)) {
	  if(((*lc_header) = (char ***) malloc((*hdr_size)*sizeof(char **))) == NULL)
	    error(ERR_MEMALLOC);
	  for(indx = 0; indx < (*hdr_size); indx++) {
	    if(((*lc_header)[indx] = (char **) malloc(3*sizeof(char *))) == NULL)
	      error(ERR_MEMALLOC);
	    if(((*lc_header)[indx][0] = (char *) malloc(MAXLEN)) == NULL ||
	       ((*lc_header)[indx][1] = (char *) malloc(MAXLEN)) == NULL ||
	       ((*lc_header)[indx][2] = (char *) malloc(MAXLEN)) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  *memsize_lc_header = *hdr_size;
	}
	else if((*hdr_size) > (*memsize_lc_header)){
	  if(((*lc_header) = (char ***) realloc((*lc_header), (*hdr_size)*sizeof(char **))) == NULL)
	    error(ERR_MEMALLOC);
	  for(indx = (*memsize_lc_header); indx < (*hdr_size); indx++) {
	    if(((*lc_header)[indx] = (char **) malloc(3*sizeof(char *))) == NULL)
	      error(ERR_MEMALLOC);
	    if(((*lc_header)[indx][0] = (char *) malloc(MAXLEN)) == NULL ||
	       ((*lc_header)[indx][1] = (char *) malloc(MAXLEN)) == NULL ||
	       ((*lc_header)[indx][2] = (char *) malloc(MAXLEN)) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  *memsize_lc_header = *hdr_size;
	}
      }
    } else if(missing_apertures && 
	      strncmp(keyword, "APERTURE", 8)==0) missing_apertures--;
    
    else if(!missing_apertures && (*num_columns)>=0 && 
	    col_ind<(*num_columns)) {
      kwsize=strlen(keyword);
      if(kwsize >= MAXIDSTRINGLENGTH) {
	error(ERR_BINARYLC_KEYWORDTOOLONG);
      }
      strcpy((*lc_column_names)[col_ind], keyword);
      if(parse_format_string(value, &((*lc_column_format)[col_ind]), NULL)) {
	free(keyword);
	return -1;
      }
      (*lc_record_size) += (*lc_column_format)[col_ind].numbits;
      if(col_ind == 0)
	(*lc_column_format)[col_ind].bytesfromstart = 0;
      else
	(*lc_column_format)[col_ind].bytesfromstart = 
	  (*lc_column_format)[col_ind-1].bytesfromstart +
	  (*lc_column_format)[col_ind-1].numbits/8;
      col_ind++;
    }
    if(keyword != NULL)
      strcpy((*lc_header)[line][0],keyword);
    else
      (*lc_header)[line][0][0] = '\0';
    if(value != NULL)
      strcpy((*lc_header)[line][1],value);
    else
      (*lc_header)[line][1][0] = '\0';
    if(comment != NULL)
      strcpy((*lc_header)[line][2],comment);
    else
      (*lc_header)[line][2][0] = '\0';
  }
  (*lc_record_size)=((*lc_record_size)+7)/8;
  if((*num_pts)<0 || (*num_columns)<0 || (*num_apertures)<0 || 
     missing_apertures || col_ind<(*num_columns)) {
    free(keyword);
    return -1;
  }
  free(keyword);
  return (*num_pts);
}

/* This program prints out the header info for a binary light curve to stdout */
void show_binarylc_header(int hdr_size, char ***lc_header)
{
  int line;
  for(line=0; line < hdr_size; line++) {
    printf("%s %s # %s\n", lc_header[line][0], lc_header[line][1], 
	   lc_header[line][2]);
  }
}



void push_binary_lc_value(char *dest, BinLC_OutputFormat *fmt,
			  void *source, int source_type)
{
  long out_long;
  double out_double;
  if(fmt->tp == INT) {
   if(fmt->has_sign) {
      switch(source_type) {
      case VARTOOLS_TYPE_DOUBLE:
	out_long = (long) (*((double *) source));
	break;
      case VARTOOLS_TYPE_CONVERTJD:
	error(ERR_BADTYPE_BINARYLC);
	break;
      case VARTOOLS_TYPE_STRING:
	error(ERR_BADTYPE_BINARYLC);
	break;
      case VARTOOLS_TYPE_INT:
	out_long = (long) (*((int *) source));
	break;
      case VARTOOLS_TYPE_FLOAT:
	out_long = (long) (*((float *) source));
	break;
      case VARTOOLS_TYPE_LONG:
	out_long = (long) (*((long *) source));
	break;
      case VARTOOLS_TYPE_SHORT:
	out_long = (long) (*((short *) source));
	break;
      case VARTOOLS_TYPE_CHAR:
	out_long = (long) (*((char *) source));
	break;
      default:
	error(ERR_BADTYPE);
      }
      push_signed_int(dest, fmt->numbits, out_long);
    }
    else {
      switch(source_type) {
      case VARTOOLS_TYPE_DOUBLE:
	out_long = (unsigned long) (*((double *) source));
	break;
      case VARTOOLS_TYPE_CONVERTJD:
	error(ERR_BADTYPE_BINARYLC);
	break;
      case VARTOOLS_TYPE_STRING:
	error(ERR_BADTYPE_BINARYLC);
	break;
      case VARTOOLS_TYPE_INT:
	out_long = (unsigned long) (*((int *) source));
	break;
      case VARTOOLS_TYPE_FLOAT:
	out_long = (unsigned long) (*((float *) source));
	break;
      case VARTOOLS_TYPE_LONG:
	out_long = (unsigned long) (*((long *) source));
	break;
      case VARTOOLS_TYPE_SHORT:
	out_long = (unsigned long) (*((short *) source));
	break;
      case VARTOOLS_TYPE_CHAR:
	out_long = (unsigned long) (*((char *) source));
	break;
      default:
	error(ERR_BADTYPE);
      }
      push_unsigned_int(dest, fmt->numbits, (unsigned long) out_long);      
    }
  }
  else {
    switch(source_type) {
    case VARTOOLS_TYPE_DOUBLE:
      out_double = (double) (*((double *) source));
      break;
    case VARTOOLS_TYPE_CONVERTJD:
      error(ERR_BADTYPE_BINARYLC);
      break;
    case VARTOOLS_TYPE_STRING:
      error(ERR_BADTYPE_BINARYLC);
      break;
    case VARTOOLS_TYPE_INT:
      out_double = (double) (*((int *) source));
      break;
    case VARTOOLS_TYPE_FLOAT:
      out_double = (double) (*((float *) source));
      break;
    case VARTOOLS_TYPE_LONG:
      out_double = (double) (*((long *) source));
      break;
    case VARTOOLS_TYPE_SHORT:
      out_double = (double) (*((short *) source));
      break;
    case VARTOOLS_TYPE_CHAR:
      out_double = (double) (*((char *) source));
      break;
    default:
      error(ERR_BADTYPE);
    }
    if(fmt->tp == PACKED_FLOAT) {
      pack_float(dest, fmt->numbits, fmt->precision,
		 fmt->offset, fmt->has_sign, out_double,
		 fmt->nan_repr);
    }
    else {
      push_float(dest, fmt->numbits, out_double);
    }
  }
}

void pop_binary_lc_value(char *source, BinLC_OutputFormat *fmt,
			 int column,
			 void * dest,
			 int dest_type) {
  long inp_long;
  double inp_double;
  if(fmt[column].tp == INT) {
    if(fmt[column].has_sign) {
      inp_long = pop_signed_int(&(source[fmt[column].bytesfromstart]),
				fmt[column].numbits) + fmt[column].offset;
    } else {
      inp_long = fmt[column].offset + pop_unsigned_int(&(source[fmt[column].bytesfromstart]),fmt[column].numbits);
    }
    switch(dest_type) {
    case VARTOOLS_TYPE_DOUBLE:
      *((double *) dest) = (double) inp_long;
      break;
    case VARTOOLS_TYPE_CONVERTJD:
      error(ERR_BADTYPE_BINARYLC);
      break;
    case VARTOOLS_TYPE_STRING:
      error(ERR_BADTYPE_BINARYLC);
      break;
    case VARTOOLS_TYPE_INT:
      *((int *) dest) = (int) inp_long;
      break;
    case VARTOOLS_TYPE_FLOAT:
      *((float *) dest) = (float) inp_long;
      break;
    case VARTOOLS_TYPE_LONG:
      *((long *) dest) = inp_long;
      break;
    case VARTOOLS_TYPE_SHORT:
      *((short *) dest) = (short) inp_long;
      break;
    case VARTOOLS_TYPE_CHAR:
      *((char *) dest) = (char) inp_long;
      break;
    }
  }
  else {
    if(fmt[column].tp == PACKED_FLOAT) {
      inp_double = unpack_float(&(source[fmt[column].bytesfromstart]),
				fmt[column].numbits,
				fmt[column].precision,
				fmt[column].offset,
				fmt[column].has_sign,
				fmt[column].nan_repr);
    }
    else {
      inp_double = pop_float(&(source[fmt[column].bytesfromstart]),
			     fmt[column].numbits);
    }
    switch(dest_type) {
    case VARTOOLS_TYPE_DOUBLE:
      *((double *) dest) = inp_double;
      break;
    case VARTOOLS_TYPE_CONVERTJD:
      error(ERR_BADTYPE_BINARYLC);
      break;
    case VARTOOLS_TYPE_STRING:
      error(ERR_BADTYPE_BINARYLC);
      break;
    case VARTOOLS_TYPE_INT:
      *((int *) dest) = (int) inp_double;
      break;
    case VARTOOLS_TYPE_FLOAT:
      *((float *) dest) = (float) inp_double;
      break;
    case VARTOOLS_TYPE_LONG:
      *((long *) dest) = (long) inp_double;
      break;
    case VARTOOLS_TYPE_SHORT:
      *((short *) dest) = (short) inp_double;
      break;
    case VARTOOLS_TYPE_CHAR:
      *((char *) dest) = (char) inp_double;
      break;
    }
  }
}

/* This function reads in a light curve that is stored in Kalo's binary
   light curve format */
int ReadBinaryLightCurve(ProgramData *p, Command *c, int lc, int threadid)
{
  FILE *infile;
  int status, hdunum, ncols, hdutype, anynulallcolumns;
  long nrows;
  int j, jold, k, i, l, N, oldsizesinglelc, anynul;
  char *nullarray, *nullarraystore;
  char **tmpstring;
  char tmpchar;

  _DataFromLightCurve *d;
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
  char *rec;

  int Nc, u, rec_ind, col_ind;

  p->NJD[threadid] = 0;

  if(p->use_lc_open_exec_command) {
    infile = ExecLCOpenCommand(p, c, lc, threadid);
    if(infile == NULL) return(ERR_FILENOTFOUND);
  }
  else if(p->readfromstdinflag == 1 && p->fileflag == 1)
    infile = stdin;
  else
    {
      if((infile = fopen(p->lcnames[lc],"r")) == NULL) {
	if(p->skipmissing) {
	  error2_noexit(ERR_FILENOTFOUND,p->lcnames[lc]);
	  return(ERR_FILENOTFOUND);
	}
	else
	  error2(ERR_FILENOTFOUND,p->lcnames[lc]);
      }
    }
  
  /* Read in the header for the binary light curve file */
  if(read_binarylightcurve_header(infile, 
				  &(p->binlc[threadid].num_apertures),
				  &(p->binlc[threadid].num_columns),
				  &(p->binlc[threadid].num_pts),
				  &(p->binlc[threadid].hdr_size),
				  &(p->binlc[threadid].lc_record_size),
				  &(p->binlc[threadid].lc_header),
				  &(p->binlc[threadid].memsize_lc_header),
				  &(p->binlc[threadid].lc_column_names),
				  &(p->binlc[threadid].lc_column_format),
				  &(p->binlc[threadid].memsize_lc_columns)) <= 0) {
    if(p->skipmissing) {
      error2_noexit(ERR_BINARYLIGHTCURVE_INVALIDFORMAT,p->lcnames[lc]);
      return(ERR_BINARYLIGHTCURVE_INVALIDFORMAT);
    }
    else
      error2(ERR_BINARYLIGHTCURVE_INVALIDFORMAT,p->lcnames[lc]);
  }
  if((rec = (char *) malloc(p->binlc[threadid].lc_record_size)) == NULL)
    error(ERR_MEMALLOC);

  if(p->lc_getcolumnsfromheader) {
    get_binarylc_header_columns(p->binlc[threadid].num_columns,
				p->binlc[threadid].lc_column_names,
				p);
  }


  p->NJD[threadid] = p->binlc[threadid].num_pts;  
  MemAllocDataFromLightCurve(p, threadid, p->binlc[threadid].num_pts);
  /* Read in the light curve data */
  for(rec_ind=0; rec_ind < p->binlc[threadid].num_pts; rec_ind++) {
    if(fread(rec, p->binlc[threadid].lc_record_size, 1, infile)!=1) {
      error2(ERR_BINARYLIGHTCURVE_MISSINGDATA,p->lcnames[lc]);
    }
    col_ind = 0;
    for(j=0; j < p->NDataFromLightCurve; j++) {
      d = &(p->DataFromLightCurve[j]);
      if(d->incolumns[0] <= 0)
	continue;
      Nc = d->Ncolumns;
      if(Nc == 0) {
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = (double ***) d->dataptr;
	  pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
			      (d->incolumns[0]-1), 
			      &((*dblptr)[threadid][rec_ind]),
			      d->datatype);
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  error(ERR_BADTYPE_BINARYLC);
	  break;
	case VARTOOLS_TYPE_STRING:
	  error(ERR_BADTYPE_BINARYLC);
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = (int ***) d->dataptr;
	  pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
			      (d->incolumns[0]-1), 
			      &((*intptr)[threadid][rec_ind]),
			      d->datatype);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = (float ***) d->dataptr;
	  pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
			      (d->incolumns[0]-1), 
			      &((*floatptr)[threadid][rec_ind]),
			      d->datatype);
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = (long ***) d->dataptr;
	  pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
			      (d->incolumns[0]-1), 
			      &((*longptr)[threadid][rec_ind]),
			      d->datatype);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = (short ***) d->dataptr;
	  pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
			      (d->incolumns[0]-1), 
			      &((*shortptr)[threadid][rec_ind]),
			      d->datatype);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = (char ***) d->dataptr;
	  pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
			      (d->incolumns[0]-1), 
			      &((*charptr)[threadid][rec_ind]),
			      d->datatype);
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      } else if (Nc > 0) {
	for(u=0; u < Nc; u++) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dbl2ptr = (double ****) d->dataptr;
	    pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
				(d->incolumns[u]-1), 
				&((*dbl2ptr)[threadid][u][rec_ind]),
 				d->datatype);
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    error(ERR_BADTYPE_BINARYLC);
	    break;
	  case VARTOOLS_TYPE_STRING:
	    error(ERR_BADTYPE_BINARYLC);
	    break;
	  case VARTOOLS_TYPE_INT:
	    int2ptr = (int ****) d->dataptr;
	    pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
				(d->incolumns[u]-1), 
				&((*int2ptr)[threadid][u][rec_ind]),
 				d->datatype);
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    float2ptr = (float ****) d->dataptr;
	    pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
				(d->incolumns[u]-1), 
				&((*float2ptr)[threadid][u][rec_ind]),
 				d->datatype);
	    break;
	  case VARTOOLS_TYPE_LONG:
	    long2ptr = (long ****) d->dataptr;
	    pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
				(d->incolumns[u]-1), 
				&((*long2ptr)[threadid][u][rec_ind]),
 				d->datatype);
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    short2ptr = (short ****) d->dataptr;
	    pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
				(d->incolumns[u]-1), 
				&((*short2ptr)[threadid][u][rec_ind]),
 				d->datatype);
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    char2ptr = (char ****) d->dataptr;
	    pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
				(d->incolumns[u]-1), 
				&((*char2ptr)[threadid][u][rec_ind]),
 				d->datatype);
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
	    pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
				(d->incolumns[u]-1), 
				&((*dbl2ptr)[threadid][rec_ind][u]),
				d->datatype);
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    error(ERR_BADTYPE_BINARYLC);
	    break;
	  case VARTOOLS_TYPE_STRING:
	    error(ERR_BADTYPE_BINARYLC);
	    break;
	  case VARTOOLS_TYPE_INT:
	    int2ptr = (int ****) d->dataptr;
	    pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
				(d->incolumns[u]-1), 
				&((*int2ptr)[threadid][rec_ind][u]),
				d->datatype);
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    float2ptr = (float ****) d->dataptr;
	    pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
				(d->incolumns[u]-1), 
				&((*float2ptr)[threadid][rec_ind][u]),
				d->datatype);
	    break;
	  case VARTOOLS_TYPE_LONG:
	    long2ptr = (long ****) d->dataptr;
	    pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
				(d->incolumns[u]-1), 
				&((*long2ptr)[threadid][rec_ind][u]),
				d->datatype);
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    short2ptr = (short ****) d->dataptr;
	    pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
				(d->incolumns[u]-1), 
				&((*short2ptr)[threadid][rec_ind][u]),
				d->datatype);
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    char2ptr = (char ****) d->dataptr;
	    pop_binary_lc_value(rec, p->binlc[threadid].lc_column_format,
				(d->incolumns[u]-1), 
				&((*char2ptr)[threadid][rec_ind][u]),
				d->datatype);
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      }
    }
  }
  free(rec);
  if(p->use_lc_open_exec_command) {
    pclose(infile);
  }
  else if(infile != stdin)
    fclose(infile);
}
#endif
      
int ReadSingleLightCurve(ProgramData *p, Command *c, int lc, int threadid)
{
  FILE *infile;
  char *line;
  char inputUTC[MAXLEN], **incols;
  int u, j, jold, k, i, l, N, oldsizesinglelc, colmax, testskip, Nc;

  int yr, mo, day, hr, min;
  double sec, jdout;

  _DataFromLightCurve *d;
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

  size_t line_size = MAXLEN;

  p->NJD[threadid] = 0;

#ifdef _USEBINARY_LC
  /* Read in the light curve from Kalo's binary format if requested */
  if(p->binarylcinput)
    return ReadBinaryLightCurve(p, c, lc, threadid);
#endif

#ifdef USECFITSIO

  if(!(p->readfromstdinflag && p->fileflag)) {
    /* Check if the end of the file is .fits */
    /* If it is, then assume the input is a binary fits table */
    j = strlen(p->lcnames[lc]);
    i = j - 5;
    if(i >= 0) {
      if(!strcmp(&(p->lcnames[lc][i]),".fits")) {
	p->is_inputlc_fits[lc] = 1;
	return ReadFitsLightCurve(p, c, lc, threadid);
      }
    }
  }

#endif

  line = malloc(line_size);

  if(p->use_lc_open_exec_command) {
    infile = ExecLCOpenCommand(p, c, lc, threadid);
    if(infile == NULL) return(ERR_FILENOTFOUND);
  }
  else if(p->readfromstdinflag == 1 && p->fileflag == 1)
    infile = stdin;
  else
    {
      if((infile = fopen(p->lcnames[lc],"r")) == NULL) {
	if(p->skipmissing) {
	  error2_noexit(ERR_FILENOTFOUND,p->lcnames[lc]);
	  return(ERR_FILENOTFOUND);
	}
	else
	  error2(ERR_FILENOTFOUND,p->lcnames[lc]);
      }
    }

  colmax = p->maxinputlccolumn;
  if(colmax > 0) {
    if((incols = (char **) malloc(colmax * sizeof(char *))) == NULL)
      error(ERR_MEMALLOC);
    for(i=0; i < colmax; i++) {
      if((incols[i] = (char *) malloc(MAXLEN)) == NULL)
	error(ERR_MEMALLOC);
    }
  }
  
  j = 0;
  l = p->Nskip;
  N = 0;
  while(j < p->Nskip ? gnu_getline(&line,&line_size,infile) >= 0 : 0)
    j++;
  while(gnu_getline(&line,&line_size,infile) >= 0)
    {
      l++;
      if(p->lcdelimtype == VARTOOLS_LC_DELIMTYPE_WHITESPACE) {
	testskip = ParseLineToColumns_testskip(line, incols, colmax, p->Nskipchar, p->skipchars);
      } else if(p->lcdelimtype == VARTOOLS_LC_DELIMTYPE_CHAR) {
	testskip = ParseLineToColumnsDelimChar_testskip(line, incols, colmax, p->Nskipchar, p->skipchars, p->delimchar);
      } else if(p->lcdelimtype == VARTOOLS_LC_DELIMTYPE_STRING) {
	testskip = ParseLineToColumnsDelimString_testskip(line, incols, colmax, p->Nskipchar, p->skipchars, p->delimstring);
      }
      if(testskip)
	continue;
      MemAllocDataFromLightCurve(p, threadid, (N+1));
      for(j=0; j < p->NDataFromLightCurve; j++) {
	d = &(p->DataFromLightCurve[j]);
	if(d->incolumns[0] <= 0)
	  continue;
	Nc = d->Ncolumns;
	if(Nc == 0) {
	  k = d->incolumns[0] - 1;
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dblptr = (double ***) d->dataptr;
	    if(d->scanformat == NULL)
	      (*dblptr)[threadid][N] = atof(incols[k]);
	    else
	      sscanf(incols[k],d->scanformat,&((*dblptr)[threadid][N]));
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    dblptr = (double ***) d->dataptr;
	    convertUTCtoJD(incols[k],d->scanformat,d->UTCindex,(double *) (&((*dblptr)[threadid][N])));
	    break;
	  case VARTOOLS_TYPE_STRING:
	    stringptr = (char ****) d->dataptr;
	    if(d->scanformat == NULL)
	      sprintf(((*stringptr)[threadid][N]),"%s",incols[k]);
	    else
	      sscanf(incols[k],d->scanformat,((*stringptr)[threadid][N]));
	    break;
	  case VARTOOLS_TYPE_INT:
	    intptr = (int ***) d->dataptr;
	    if(d->scanformat == NULL)
	      (*intptr)[threadid][N] = atoi(incols[k]);
	    else
	      sscanf(incols[k],d->scanformat,&((*intptr)[threadid][N]));
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    shortptr = (short ***) d->dataptr;
	    if(d->scanformat == NULL)
	      (*shortptr)[threadid][N]= atoi(incols[k]);
	    else
	      sscanf(incols[k],d->scanformat,&((*shortptr)[threadid][N]));
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    floatptr = (float ***) d->dataptr;
	    if(d->scanformat == NULL)
	      (*floatptr)[threadid][N] = atof(incols[k]);
	    else
	      sscanf(incols[k],d->scanformat,&((*floatptr)[threadid][N]));
	    break;
	  case VARTOOLS_TYPE_LONG:
	    longptr = (long ***) d->dataptr;
	    if(d->scanformat == NULL)
	      (*longptr)[threadid][N] = atol(incols[k]);
	    else
	      sscanf(incols[k],d->scanformat,&((*longptr)[threadid][N]));
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    charptr = (char ***) d->dataptr;
	    if(d->scanformat == NULL)
	      (*charptr)[threadid][N] = incols[k][0];
	    else
	      sscanf(incols[k],d->scanformat,&((*charptr)[threadid][N]));
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	} else if(Nc > 0) {
	  for(u = 0; u < Nc; u++) {
	    k = d->incolumns[u] - 1;
	    switch(d->datatype) {
	    case VARTOOLS_TYPE_DOUBLE:
	      dbl2ptr = (double ****) d->dataptr;
	      if(d->scanformat == NULL)
		(*dbl2ptr)[threadid][u][N] = atof(incols[k]);
	      else
		sscanf(incols[k],d->scanformat,&((*dbl2ptr)[threadid][u][N]));
	      break;
	    case VARTOOLS_TYPE_CONVERTJD:
	      dbl2ptr = (double ****) d->dataptr;
	      convertUTCtoJD(incols[k],d->scanformat,d->UTCindex,(double *) (&((*dbl2ptr)[threadid][u][N])));
	      break;
	    case VARTOOLS_TYPE_STRING:
	      string2ptr = (char *****) d->dataptr;
	      if(d->scanformat == NULL)
		sprintf(((*string2ptr)[threadid][u][N]),"%s",incols[k]);
	      else
		sscanf(incols[k],d->scanformat,((*string2ptr)[threadid][u][N]));
	      break;
	    case VARTOOLS_TYPE_INT:
	      int2ptr = (int ****) d->dataptr;
	      if(d->scanformat == NULL)
		(*int2ptr)[threadid][u][N] = atoi(incols[k]);
	      else
		sscanf(incols[k],d->scanformat,&((*int2ptr)[threadid][u][N]));
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      short2ptr = (short ****) d->dataptr;
	      if(d->scanformat == NULL)
		(*short2ptr)[threadid][u][N]= atoi(incols[k]);
	      else
		sscanf(incols[k],d->scanformat,&((*short2ptr)[threadid][u][N]));
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      float2ptr = (float ****) d->dataptr;
	      if(d->scanformat == NULL)
		(*float2ptr)[threadid][u][N] = atof(incols[k]);
	      else
		sscanf(incols[k],d->scanformat,&((*float2ptr)[threadid][u][N]));
	      break;
	    case VARTOOLS_TYPE_LONG:
	      long2ptr = (long ****) d->dataptr;
	      if(d->scanformat == NULL)
		(*long2ptr)[threadid][u][N] = atol(incols[k]);
	      else
		sscanf(incols[k],d->scanformat,&((*long2ptr)[threadid][u][N]));
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      char2ptr = (char ****) d->dataptr;
	      if(d->scanformat == NULL)
		(*char2ptr)[threadid][u][N] = incols[k][0];
	      else
		sscanf(incols[k],d->scanformat,&((*char2ptr)[threadid][u][N]));
	      break;
	    default:
		error(ERR_BADTYPE);
	    }
	  }
	} else {
	  for(u = 0; u < (-Nc); u++) {
	    k = d->incolumns[u] - 1;
	    switch(d->datatype) {
	    case VARTOOLS_TYPE_DOUBLE:
	      dbl2ptr = (double ****) d->dataptr;
	      if(d->scanformat == NULL)
		(*dbl2ptr)[threadid][N][u] = atof(incols[k]);
	      else
		sscanf(incols[k],d->scanformat,&((*dbl2ptr)[threadid][N][u]));
	      break;
	    case VARTOOLS_TYPE_CONVERTJD:
	      dbl2ptr = (double ****) d->dataptr;
	      convertUTCtoJD(incols[k],d->scanformat,d->UTCindex,(double *) (&((*dbl2ptr)[threadid][N][u])));
	      break;
	    case VARTOOLS_TYPE_STRING:
	      string2ptr = (char *****) d->dataptr;
	      if(d->scanformat == NULL)
		sprintf(((*string2ptr)[threadid][N][u]),"%s",incols[k]);
	      else
		sscanf(incols[k],d->scanformat,((*string2ptr)[threadid][N][u]));
	      break;
	    case VARTOOLS_TYPE_INT:
	      int2ptr = (int ****) d->dataptr;
	      if(d->scanformat == NULL)
		(*int2ptr)[threadid][N][u] = atoi(incols[k]);
	      else
		sscanf(incols[k],d->scanformat,&((*int2ptr)[threadid][N][u]));
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      short2ptr = (short ****) d->dataptr;
	      if(d->scanformat == NULL)
		(*short2ptr)[threadid][N][u] = atoi(incols[k]);
	      else
		sscanf(incols[k],d->scanformat,&((*short2ptr)[threadid][N][u]));
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      float2ptr = (float ****) d->dataptr;
	      if(d->scanformat == NULL)
		(*float2ptr)[threadid][N][u] = atof(incols[k]);
	      else
		sscanf(incols[k],d->scanformat,&((*float2ptr)[threadid][N][u]));
	      break;
	    case VARTOOLS_TYPE_LONG:
	      long2ptr = (long ****) d->dataptr;
	      if(d->scanformat == NULL)
		(*long2ptr)[threadid][N][u] = atol(incols[k]);
	      else
		sscanf(incols[k],d->scanformat,&((*long2ptr)[threadid][N][u]));
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      char2ptr = (char ****) d->dataptr;
	      if(d->scanformat == NULL)
		(*char2ptr)[threadid][N][u] = incols[k][0];
	      else
		sscanf(incols[k],d->scanformat,&((*char2ptr)[threadid][N][u]));
	      break;
	    default:
	      error(ERR_BADTYPE);
	    }
	  }
	}
      }
      N++;
    }

  if(colmax > 0) {
    for(i=0; i < colmax; i++)
      free(incols[i]);
    free(incols);
  }

  p->NJD[threadid] = N;
  if(p->readimagestring)
    {
      for(i=0;i<N;i++)
	p->stringid_idx[threadid][i] = i;
    }
  if(p->use_lc_open_exec_command)
    pclose(infile);
  else if(infile != stdin)
    fclose(infile);

  /* Fill in any variables calculated from analytic expressions */
  for(j=0; j < p->NDataFromLightCurve; j++) {
    d = &(p->DataFromLightCurve[j]);
    if(d->incolumns[0] <= 0 && d->expression != NULL)
      {
	if(d->scanformat != NULL) {
	  Nc = d->Ncolumns;
	  if(Nc != 0)
	    error(ERR_CODEERROR);
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    for(i=0; i < p->NJD[threadid]; i++) {
	      dblptr = (double ***) d->dataptr;
	      (*dblptr)[threadid][i] = 
		EvaluateExpression(lc, threadid, i, d->expression);
	    }
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    for(i=0; i < p->NJD[threadid]; i++) {
	      floatptr = (float ***) d->dataptr;
	      (*floatptr)[threadid][i] = 
		(float) EvaluateExpression(lc, threadid, i, d->expression);
	    }
	    break;
	  case VARTOOLS_TYPE_INT:
	    for(i=0; i < p->NJD[threadid]; i++) {
	      intptr = (int ***) d->dataptr;
	      (*intptr)[threadid][i] = 
		(int) EvaluateExpression(lc, threadid, i, d->expression);
	    }
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    for(i=0; i < p->NJD[threadid]; i++) {
	      shortptr = (short ***) d->dataptr;
	      (*shortptr)[threadid][i] = 
		(short) EvaluateExpression(lc, threadid, i, d->expression);
	    }
	    break;
	  case VARTOOLS_TYPE_LONG:
	    for(i=0; i < p->NJD[threadid]; i++) {
	      longptr = (long ***) d->dataptr;
	      (*longptr)[threadid][i] = 
		(long) EvaluateExpression(lc, threadid, i, d->expression);
	    }
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
	else {
	  Nc = d->Ncolumns;
	  if(Nc != 0)
	    error(ERR_CODEERROR);
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    for(i=0; i < p->NJD[threadid]; i++) {
	      dblptr = (double ***) d->dataptr;
	      (*dblptr)[threadid][i] = 0.;
	    }
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    for(i=0; i < p->NJD[threadid]; i++) {
	      floatptr = (float ***) d->dataptr;
	      (*floatptr)[threadid][i] = 0.;
	    }
	    break;
	  case VARTOOLS_TYPE_INT:
	    for(i=0; i < p->NJD[threadid]; i++) {
	      intptr = (int ***) d->dataptr;
	      (*intptr)[threadid][i] = 0;
	    }
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    for(i=0; i < p->NJD[threadid]; i++) {
	      shortptr = (short ***) d->dataptr;
	      (*shortptr)[threadid][i] = 0;
	    }
	    break;
	  case VARTOOLS_TYPE_LONG:
	    for(i=0; i < p->NJD[threadid]; i++) {
	      longptr = (long ***) d->dataptr;
	      (*longptr)[threadid][i] = 0;
	    }
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      }
  }

  /* Fill the decorr matrix if we are decorrelating */
  if(p->decorrflag)
    Filldecorr_matrix(p,c,threadid);
  
  free(line);

  return 0;
}


int ReadAllLightCurves(ProgramData *p, Command *c)
{
  FILE *infile;
  char line[MAXLEN];
  char inputUTC[MAXLEN];
  int j, k, i, N, l, lc, jold, isempty = 0;

  for(lc=0;lc < p->Nlcs ;lc++)
    {

#ifdef USECFITSIO

      /* Check if the end of the file is .fits */
      /* If it is, then assume the input is a binary fits table */
      j = strlen(p->lcnames[lc]);
      i = j - 5;
      if(i >= 0) {
	if(!strcmp(&(p->lcnames[lc][i]),".fits")) {
	  p->is_inputlc_fits[lc] = 1;
	  if(ReadFitsLightCurve(p, c, lc, lc)) {
	    isempty++;
	  }
	  continue;
	}
      }
#endif

      if(ReadSingleLightCurve(p, c, lc, lc))
	isempty++;
    }
  return isempty;
}

/* Moves all the light curve data associated with index isrc to
   the index idest */
void MoveLCData(ProgramData *p, Command *c, int isrc, int idest) {
  int i, Nc;
  _DataFromLightCurve *d;
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
  if(p->readimagestring) {
    p->stringid_idx[idest] = p->stringid_idx[isrc];
    p->stringid[idest] = p->stringid[isrc];
  }
  p->t[idest] = p->t[isrc];
  p->mag[idest] = p->mag[isrc];
  p->sig[idest] = p->sig[isrc];
  p->NJD[idest] = p->NJD[isrc];
  for(i=0; i < p->NDataFromLightCurve; i++) {
    d = &(p->DataFromLightCurve[i]);
    Nc = d->Ncolumns;
    if(Nc == 0) {
      switch(d->datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dblptr = (double ***) d->dataptr;
	(*dblptr)[idest] = (*dblptr)[isrc];
	break;
      case VARTOOLS_TYPE_CONVERTJD:
	dblptr = (double ***) d->dataptr;
	(*dblptr)[idest] = (*dblptr)[isrc];
      case VARTOOLS_TYPE_STRING:
	stringptr = (char ****) d->dataptr;
	(*stringptr)[idest] = (*stringptr)[isrc];
	break;
      case VARTOOLS_TYPE_INT:
	intptr = (int ***) d->dataptr;
	(*intptr)[idest] = (*intptr)[isrc];
	break;
      case VARTOOLS_TYPE_SHORT:
	shortptr = (short ***) d->dataptr;
	(*shortptr)[idest] = (*shortptr)[isrc];
	break;
      case VARTOOLS_TYPE_FLOAT:
	floatptr = (float ***) d->dataptr;
	(*floatptr)[idest] = (*floatptr)[isrc];
	break;
      case VARTOOLS_TYPE_LONG:
	longptr = (long ***) d->dataptr;
	(*longptr)[idest] = (*longptr)[isrc];
	break;
      case VARTOOLS_TYPE_CHAR:
	charptr = (char ***) d->dataptr;
	(*charptr)[idest] = (*charptr)[isrc];
	break;
      default:
	error(ERR_BADTYPE);
      }
    } else if(Nc > 0) {
      switch(d->datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dbl2ptr = (double ****) d->dataptr;
	(*dbl2ptr)[idest] = (*dbl2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_CONVERTJD:
	dbl2ptr = (double ****) d->dataptr;
	(*dbl2ptr)[idest] = (*dbl2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_STRING:
	string2ptr = (char *****) d->dataptr;
	(*string2ptr)[idest] = (*string2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_INT:
	int2ptr = (int ****) d->dataptr;
	(*int2ptr)[idest] = (*int2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_SHORT:
	short2ptr = (short ****) d->dataptr;
	(*short2ptr)[idest] = (*short2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_FLOAT:
	float2ptr = (float ****) d->dataptr;
	(*float2ptr)[idest] = (*float2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_LONG:
	long2ptr = (long ****) d->dataptr;
	(*long2ptr)[idest] = (*long2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_CHAR:
	char2ptr = (char ****) d->dataptr;
	(*char2ptr)[idest] = (*char2ptr)[isrc];
	break;
      default:
	error(ERR_BADTYPE);
	
      }
    } else {
      switch(d->datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dbl2ptr = (double ****) d->dataptr;
	(*dbl2ptr)[idest] = (*dbl2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_CONVERTJD:
	dbl2ptr = (double ****) d->dataptr;
	(*dbl2ptr)[idest] = (*dbl2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_STRING:
	string2ptr = (char *****) d->dataptr;
	(*string2ptr)[idest] = (*string2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_INT:
	int2ptr = (int ****) d->dataptr;
	(*int2ptr)[idest] = (*int2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_SHORT:
	short2ptr = (short ****) d->dataptr;
	(*short2ptr)[idest] = (*short2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_FLOAT:
	float2ptr = (float ****) d->dataptr;
	(*float2ptr)[idest] = (*float2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_LONG:
	long2ptr = (long ****) d->dataptr;
	(*long2ptr)[idest] = (*long2ptr)[isrc];
	break;
      case VARTOOLS_TYPE_CHAR:
	char2ptr = (char ****) d->dataptr;
	(*char2ptr)[idest] = (*char2ptr)[isrc];
	break;
      default:
	error(ERR_BADTYPE);
	
      }
    }
  }
  if(p->decorrflag)
    Filldecorr_matrix(p,c,idest);
}

/* This function removes empty light curves from the data set,
   compacting all vectors */
void RemoveEmptyLightCurves(ProgramData *p, Command *c) {
  int lc_orig, lc_good;
  lc_good = 0;
  for(lc_orig=0; lc_orig < p->Nlcs; lc_orig++) {
    if(p->NJD[lc_orig] != 0) {
      if(lc_good != lc_orig) {
	MoveInputListData(p, lc_orig, lc_good);
	MoveLCData(p, c, lc_orig, lc_good);
      }
      lc_good++;
    }
  }
  p->Nlcs = lc_good;
}

/* Show the expected format of the input light curve(s) */
void printinputlcformat(ProgramData *p, FILE *outfile)
{
  int i, j, Np, Nc, u;
  _DataFromLightCurve *d;

  fprintf(outfile,
	  "#Expected Format Of Input Light Curves\n");
  fprintf(outfile,
	  "#-------------------------------------\n");

  Np = p->maxinputlccolumn;

  for(i=1; i <= Np; i++) {
    for(j=0; j < p->NDataFromLightCurve; j++) {
      d = &(p->DataFromLightCurve[j]);
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

