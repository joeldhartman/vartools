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
#include <Rinternals.h>
#include <Rembedded.h>
#include <R_ext/Parse.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/socket.h>
#include <sys/types.h>

#include "runR.h"


/* Some functions we will need from the other files */
void SetVariable_Value_Double(int lcindex, int threadindex, int jdindex, _Variable *var, double val);
double EvaluateVariable_Double(int lcindex, int threadindex, int jdindex, _Variable *var);
void MemAllocDataFromLightCurveMidProcess(ProgramData *p, int threadid, int Nterm);

#ifndef DO_ERROR_MEMALLOC
#define DO_ERROR_MEMALLOC do {fprintf(stderr,"Memory Allocation Error\n"); exit(1);} while(1)
#endif

void simpleprinttostring(OutText *text, const char *stoadd)
{
  int l, lold, j, k, k1, k2;
  l = strlen(stoadd);
  lold = text->len_s;
  if(!text->space)
    {
      text->Nchar_cur_line = 0;
      text->space = MAXLEN;
      if((text->s = malloc(MAXLEN)) == NULL)
	DO_ERROR_MEMALLOC;
      text->s[0] = '\0';
    }

  k = 0;
  j = lold;
  do {
    k1 = k;
    while(stoadd[k] != '\0')
      k++;
    while(text->len_s + (k-k1) + 2 >= text->space) {
      text->space = text->space * 2;
      if((text->s = realloc(text->s, text->space)) == NULL)
	DO_ERROR_MEMALLOC;
    }
    while(k1 < k) {
      text->s[j] = stoadd[k1];
      text->Nchar_cur_line += 1;
      j++;
      k1++;
    }
    text->s[j] = '\0';
    text->len_s = j;
  } while(stoadd[k] != '\0');
}

void simpleprinttostring_tabindent(OutText *text, const char *stoadd)
{
  int l, lold, j, k, k1, k2, nnewtab;
  l = strlen(stoadd);
  lold = text->len_s;
  if(!text->space)
    {
      text->Nchar_cur_line = 0;
      text->space = MAXLEN;
      if((text->s = malloc(MAXLEN)) == NULL)
	DO_ERROR_MEMALLOC;
      text->s[0] = '\0';
    }

  k = 0;
  j = lold;
  k1 = k;
  nnewtab = 0;
  while(stoadd[k] != '\0') {
    if(stoadd[k] != '\n')
      k++;
    else {
      nnewtab++;
      k++;
    }
  }
  while(text->len_s + (k-k1) + 2 + nnewtab >= text->space) {
    text->space = text->space * 2;
    if((text->s = realloc(text->s, text->space)) == NULL)
      DO_ERROR_MEMALLOC;
  }
  while(k1 < k) {
    text->s[j] = stoadd[k1];
    text->Nchar_cur_line += 1;
    if(stoadd[k1] == '\n') {
      j++;
      text->s[j] = '\t';
      text->Nchar_cur_line = 1;
    }
    j++;
    k1++;
  }
  text->s[j] = '\0';
  text->len_s = j;
}



void CreateVartoolsRUserFunctionString(ProgramData *p, _RCommand *cparent, OutText *text, int cid)
{
  int i, k;
  int isterm;
  char tmpstr[MAXLEN];
  char tmpstrtype[MAXLEN];
  char astypecpy[] = ", copy=False)";
  char astypencpy[] = ")";
  char *astypever;
  _RCommand *c;

  if(!cid) c = cparent;
  else {
    c = ((_RCommand **) cparent->childcommandptrs)[cid-1];
  }
  /* The function will be named
     _VARTOOLS_R_USERFUNCTION_%d  where %d is the command number.
     It will take as input a List of variables of the name
     _VARTOOLS_VARIABLES_LIST. */
  simpleprinttostring(text,"VARTOOLS_R_USERFUNCTION_");
  sprintf(tmpstr,"%d",c->cnum);
  simpleprinttostring(text,tmpstr);
  simpleprinttostring(text," <- function(VARTOOLS_VARIABLES_LIST) {\n");

  k = -1;

  /* We have to transpose the data if we are receiving all of it from VARTOOLS */
  if(c->RequireReadAll) {
    simpleprinttostring(text,"\tVARTOOLS_INPUT_NLCS <- length(VARTOOLS_VARIABLES_LIST)\n");
    for(i = 0; i < c->Nvars; i++) {
      simpleprinttostring(text,"\t");
      simpleprinttostring(text,c->vars[i]->varname);
      simpleprinttostring(text," <- list()\n");
      simpleprinttostring(text,"\tfor (VARTOOLS_ITER_VAR in 1:VARTOOLS_INPUT_NLCS) {\n");
      simpleprinttostring(text,"\t\t");
      simpleprinttostring(text,c->vars[i]->varname);
      simpleprinttostring(text,"[[VARTOOLS_ITER_VAR]] <- ");
      simpleprinttostring(text,"VARTOOLS_VARIABLES_LIST[[VARTOOLS_ITER_VAR]][[");
      sprintf(tmpstr,"%d",i+1);
      simpleprinttostring(text,tmpstr);
      simpleprinttostring(text,"]]\n");
      simpleprinttostring(text,"\t}\n");
    }
  } else {

    /* Set all of the variables to the appropriate values from the List */
    for(i = 0; i < c->Nvars; i++) {
      simpleprinttostring(text,"\t");
      simpleprinttostring(text,c->vars[i]->varname);
      simpleprinttostring(text," <- VARTOOLS_VARIABLES_LIST[[");
      sprintf(tmpstr,"%d",i+1);
      simpleprinttostring(text,tmpstr);
      simpleprinttostring(text,"]]\n");
    }
  }

  /* Copy over the code supplied by the user; note that because we are in a
     function here, add a leading tab */
  simpleprinttostring(text,"\t");
  simpleprinttostring_tabindent(text,c->Rcommandstring);
  
  if(text->s[text->len_s-1] != '\t')
    simpleprinttostring_tabindent(text,"\n");
  simpleprinttostring(text,"\n");


  if(c->RequireReadAll) {
    /* Create any output variables that don't exist, these need to be lists */
    for(i=0; i < c->Nvars_outonly; i++) {
      simpleprinttostring(text,"\tif(!exists(\"");
      simpleprinttostring(text,c->outonlyvars[i]->varname);
      simpleprinttostring(text,"\")) {\n");
      simpleprinttostring(text,"\t\t");
      simpleprinttostring(text,c->outonlyvars[i]->varname);
      simpleprinttostring(text," <- list()\n");
      simpleprinttostring(text,"\t\tfor (VARTOOLS_ITER_VAR in 1:VARTOOLS_INPUT_NLCS) {\n");
      simpleprinttostring(text,"\t\t\t");
      simpleprinttostring(text,c->outonlyvars[i]->varname);
      simpleprinttostring(text,"[[VARTOOLS_ITER_VAR]] <- c(");
      switch(c->outonlyvars[i]->datatype) {
      case VARTOOLS_TYPE_CONVERTJD:
      case VARTOOLS_TYPE_DOUBLE:
      case VARTOOLS_TYPE_FLOAT:
	simpleprinttostring(text,"0.0");
	break;
      case VARTOOLS_TYPE_INT:
      case VARTOOLS_TYPE_LONG:
      case VARTOOLS_TYPE_SHORT:
	simpleprinttostring(text,"0");
	break;
      case VARTOOLS_TYPE_STRING:
      case VARTOOLS_TYPE_CHAR:
	simpleprinttostring(text,"\"\"");
	break;
      default:
	simpleprinttostring(text,"0.0");
	break;
      }
      simpleprinttostring(text,")\n");
      simpleprinttostring(text,"\t\t}\n\t}\n");
    }

    simpleprinttostring(text,"\tVARTOOLS_VARIABLES_OUTPUTLIST <- list()\n");
    simpleprinttostring(text,"\tfor (VARTOOLS_ITER_VAR in 1:VARTOOLS_INPUT_NLCS) {\n");
    simpleprinttostring(text,"\t\tVARTOOLS_VARIABLES_OUTPUTLIST_TMP <- list(");
    isterm = 0;
    for(i=0; i < c->Nvars; i++) {
      if(i > 0) {
	simpleprinttostring(text,", ");
      }
      simpleprinttostring(text,c->vars[i]->varname);
      simpleprinttostring(text,"[[VARTOOLS_ITER_VAR]]");
      isterm = 1;
    }
    /* Append any variables that are created by this process to the list */
    for(i=0; i < c->Nvars_outonly; i++) {
      if(isterm)
	simpleprinttostring(text,", ");
      simpleprinttostring(text,c->outonlyvars[i]->varname);
      simpleprinttostring(text,"[[VARTOOLS_ITER_VAR]]");
    }
    simpleprinttostring(text,")\n");
    simpleprinttostring(text,"\t\tVARTOOLS_VARIABLES_OUTPUTLIST[[VARTOOLS_ITER_VAR]] <- VARTOOLS_VARIABLES_OUTPUTLIST_TMP\n");
    simpleprinttostring(text,"\t}\n");
  } else {
    /* Create any output variables that don't exist */
    for(i=0; i < c->Nvars_outonly; i++) {
      simpleprinttostring(text,"\tif(!exists(\"");
      simpleprinttostring(text,c->outonlyvars[i]->varname);
      simpleprinttostring(text,"\")) {");
      simpleprinttostring(text,"\t\t");
      simpleprinttostring(text,c->outonlyvars[i]->varname);
      simpleprinttostring(text," <- c(");
      switch(c->outonlyvars[i]->datatype) {
      case VARTOOLS_TYPE_CONVERTJD:
      case VARTOOLS_TYPE_DOUBLE:
      case VARTOOLS_TYPE_FLOAT:
	simpleprinttostring(text,"0.0");
	break;
      case VARTOOLS_TYPE_INT:
      case VARTOOLS_TYPE_LONG:
      case VARTOOLS_TYPE_SHORT:
	simpleprinttostring(text,"0");
	break;
      case VARTOOLS_TYPE_STRING:
      case VARTOOLS_TYPE_CHAR:
	simpleprinttostring(text,"\"\"");
	break;
      default:
	simpleprinttostring(text,"0.0");
	break;
      }
      simpleprinttostring(text,")\n");
      simpleprinttostring(text,"\t}\n");
    }

    simpleprinttostring(text,"\tVARTOOLS_VARIABLES_OUTPUTLIST <- list(");
    
    isterm = 0;
    for(i=0; i < c->Nvars; i++) {
      if(i > 0) {
	simpleprinttostring(text,", ");
      }
      simpleprinttostring(text,c->vars[i]->varname);
      isterm = 1;
    }
    /* Append any variables that are created by this process to the list */
    for(i=0; i < c->Nvars_outonly; i++) {
      if(isterm)
	simpleprinttostring(text,", ");
      simpleprinttostring(text,c->outonlyvars[i]->varname);
    }
    simpleprinttostring(text,")\n");
  }
  simpleprinttostring(text,"\treturn(VARTOOLS_VARIABLES_OUTPUTLIST)\n");
  simpleprinttostring(text,"}\n");

  //fprintf(stderr,"%s", text->s);
}

void CleanUpRObjectContainerVariables(_RObjectContainer *Rcontainer) {
  int i;
  if(Rcontainer->data != NULL) {
    for(i=0; i < Rcontainer->Nvars; i++) {
      if(Rcontainer->data[i].dataptr != NULL)
	free(Rcontainer->data[i].dataptr);
    }
    free(Rcontainer->data);
  }
  Rcontainer->data = NULL;
  Rcontainer->Nvars = 0;
}

void CleanUpRObjectContainer(_RObjectContainer **Rcontainer) {
  int i;
  if(Rcontainer == NULL) return;
  if(*Rcontainer == NULL) return;
  CleanUpRObjectContainerVariables(*Rcontainer);
  free((*Rcontainer));
  *Rcontainer = NULL;
}

_RObjectContainer *CreateRObjectContainer(void) {
  _RObjectContainer *Rcontainer;
  if((Rcontainer = (_RObjectContainer *) malloc(sizeof(_RObjectContainer))) == NULL)
    DO_ERROR_MEMALLOC;
  Rcontainer->Nvars = 0;
  Rcontainer->data = NULL;
  Rcontainer->Nfunc = 0;
  return Rcontainer;
}



int InitializeR(ProgramData *p, _RCommand *c, int threadindex)
{
  OutText usercodetext;
  _RObjectContainer *Rcontainer;
  char tmpstr[MAXLEN];
  int i, j, testneedreadall;

  R_len_t ii;

  int r_argc = 2;
  char *r_argv[] = {"R", "--slave"};

  SEXP cmdSexp;
  SEXP cmdexp;
  SEXP ans;
  ParseStatus status;

  usercodetext.s = NULL;
  usercodetext.space = 0;
  usercodetext.len_s = 0;
  usercodetext.Nchar_cur_line = 0;

#ifdef RHOME
  if(!getenv("R_HOME")) {
    char r_home_var[MAXLEN];
    sprintf(r_home_var,"R_HOME=" RHOME);
    putenv(r_home_var);
  }
#endif

  Rf_initEmbeddedR(r_argc, r_argv);

  if(c->len_Rinitializationtextstring > 0) {
    if(c->Rinitializationtext[0] != '\0') {
      simpleprinttostring(&usercodetext, c->Rinitializationtext);
    }
  }
  simpleprinttostring(&usercodetext,"\n");
	
  Rcontainer = CreateRObjectContainer();

  ((_RObjectContainer **) c->Robjects)[threadindex] = Rcontainer;

  testneedreadall = 0;
  if(c->RequireReadAll) {
    testneedreadall = 1;
  }
  if(!testneedreadall && c->Nchildren > 0) {
    for(j=0; j < c->Nchildren; j++) {
      if(((_RCommand **)(c->childcommandptrs))[j]->RequireReadAll) {
	testneedreadall = 1;
	break;
      }
    }
  }
  if(testneedreadall) {
    for(j=1; j < p->Nlcs; j++) {
      ((_RObjectContainer **) c->Robjects)[j] = CreateRObjectContainer();
    }
  }

  Rcontainer->Nfunc = c->Nchildren + 1;


  /* Create the functions to run the R code supplied by the user */
  for(i=0; i < c->Nchildren + 1; i++) {
    CreateVartoolsRUserFunctionString(p, c, &usercodetext, i);
  }

  /* Parse the R code and evaluate it */

  cmdSexp = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT(cmdSexp, 0, mkChar(usercodetext.s));
  cmdexp = PROTECT(R_ParseVector(cmdSexp, -1, &(status), R_NilValue));
  if(status != PARSE_OK) {
    UNPROTECT(2);
    fprintf(stderr,"Error parsing the following R code from vartools:\n\n%s", usercodetext.s);
    exit(1);
  }
  for(i=0; i < length(cmdexp); i++) {
    ans = eval(VECTOR_ELT(cmdexp, i), R_GlobalEnv);
  }
  UNPROTECT(2);
 
  return 0;

  /* Everything should be setup now to call this function through the
     RunRProcessingLoop */
}


#define _EXIT_READ_VARIABLES do { if(tmpinpstr != NULL) free(tmpinpstr); if(tmpinpstr2 != NULL) free(tmpinpstr2); UNPROTECT(Nprotected); *statusval = 1; return VariableList; } while(0)

SEXP ReadVariablesFromSocketIntoR(ProgramData *p, _RCommand *c, 
				      int threadindex, int isall, int lcnum, int *statusval) {
  _RObjectContainer *Rcontainer = NULL;

  int Nvars, i, j, lenvec, k;
  char datatype;
  size_t databytesize;
  char *tmpinpstr = NULL;
  int sizetmpinpstr = 0;
  char *tmpinpstr2 = NULL;
  int sizetmpinpstr2 = 0;
  int Nprotected = 0;
  SEXP VariableList;
  SEXP InputVariable;
  size_t sizedataread;
  size_t sizedatareadtot;

  *statusval = 0;

  if(!isall) {
    Rcontainer = ((_RObjectContainer **) c->Robjects)[threadindex];

    if(Rcontainer->data != NULL)
      CleanUpRObjectContainerVariables(Rcontainer);
  } else {
    Rcontainer = ((_RObjectContainer **) c->Robjects)[lcnum];

    if(Rcontainer->data != NULL)
      CleanUpRObjectContainerVariables(Rcontainer);
  }

  if(read(c->sockets[threadindex][1], &Nvars, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Error reading the number of variables to read in the parent in the ReadVariablesFromSocketIntoR function.\n");
    _EXIT_READ_VARIABLES;
  }


  /* The only argument is a list of input variables; First create the list */
  VariableList = PROTECT(allocVector(VECSXP, Nvars));
  Nprotected++;
  
  Rcontainer->Nvars = Nvars;
  if(Nvars >= 0) {
    if((Rcontainer->data = (_RArrayData *) malloc(Nvars*sizeof(_RArrayData))) == NULL) {
      fprintf(stderr,"Memory Allocation Error in VARTOOLS R sub-process.\n");
      _EXIT_READ_VARIABLES;
    }
    for(i=0; i < Nvars; i++) Rcontainer->data[i].dataptr = NULL;
    for(i=0; i < Nvars; i++) {
      if(read(c->sockets[threadindex][1], &datatype, sizeof(char)) < sizeof(char)) {
	fprintf(stderr,"Error reading datatype for expected variable index %d in the ReadVariablesFromSocketIntoR function.\n", i);
	_EXIT_READ_VARIABLES;
      }
      if(datatype != VARTOOLS_TYPE_STRING && datatype != VARTOOLS_TYPE_CHAR) {
	switch(datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  databytesize = sizeof(double);
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  databytesize = sizeof(double);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  databytesize = sizeof(float);
	  break;
	case VARTOOLS_TYPE_INT:
	  databytesize = sizeof(int);
	  break;
	case VARTOOLS_TYPE_LONG:
	  databytesize = sizeof(long);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  databytesize = sizeof(short);
	  break;
	default:
	  fprintf(stderr,"Error: invalid datatype received for variable index %d in the ReadVariablesFromSocketIntoR function\n", i);
	  _EXIT_READ_VARIABLES;
	}
	if(read(c->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	  fprintf(stderr,"Error reading length of vector for expected variable index %d in the ReadVariablesFromSocketIntoR function.\n", i);
	  _EXIT_READ_VARIABLES;
	}
	Rcontainer->data[i].datatype = datatype;
	if(lenvec <= 0) {
	  continue;
	}
	if((Rcontainer->data[i].dataptr = (void *) malloc(((size_t) lenvec)*databytesize)) == NULL) {
	  fprintf(stderr,"Memory Allocation Error in VARTOOLS R sub-process.\n");
	  _EXIT_READ_VARIABLES;
	}
	
	sizedatareadtot = 0;
	while(sizedatareadtot < (((size_t) lenvec)*databytesize)) {
	  sizedataread = read(c->sockets[threadindex][1], (void *) (((char *) (Rcontainer->data[i].dataptr)) + sizedatareadtot), (((size_t) lenvec)*databytesize) - sizedatareadtot);
	  sizedatareadtot += sizedataread;
	}
	if(sizedatareadtot < (((size_t) lenvec)*databytesize)) {
	  fprintf(stderr,"Error reading data for variable index %d in the ReadVariablesFromSocketIntoR function.\n", i);
	  _EXIT_READ_VARIABLES;
	}
	switch(datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  InputVariable = PROTECT(allocVector(REALSXP, lenvec));
	  Nprotected++;
	  memcpy(REAL(InputVariable), Rcontainer->data[i].dataptr, (((size_t) lenvec)*databytesize));
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  InputVariable = PROTECT(allocVector(REALSXP, lenvec));
	  Nprotected++;
	  memcpy(REAL(InputVariable), Rcontainer->data[i].dataptr, (((size_t) lenvec)*databytesize));
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  InputVariable = PROTECT(allocVector(REALSXP, lenvec));
	  Nprotected++;
	  for(j=0; j < lenvec; j++) {
	    REAL(InputVariable)[j] = ((float *) Rcontainer->data[i].dataptr)[j];
	  }
	  break;
	case VARTOOLS_TYPE_INT:
	  InputVariable = PROTECT(allocVector(INTSXP, lenvec));
	  Nprotected++;
	  memcpy(INTEGER(InputVariable), Rcontainer->data[i].dataptr, (((size_t) lenvec)*databytesize));
	  break;
	case VARTOOLS_TYPE_LONG:
	  InputVariable = PROTECT(allocVector(INTSXP, lenvec));
	  Nprotected++;
	  for(j=0; j < lenvec; j++) {
	    INTEGER(InputVariable)[j] = ((long *) Rcontainer->data[i].dataptr)[j];
	  }
	  break;
	case VARTOOLS_TYPE_SHORT:
	  InputVariable = PROTECT(allocVector(INTSXP, lenvec));
	  Nprotected++;
	  for(j=0; j < lenvec; j++) {
	    INTEGER(InputVariable)[j] = ((short *) Rcontainer->data[i].dataptr)[j];
	  }
	  break;
	default:
	  fprintf(stderr,"Error: invalid datatype received for variable index %d in the ReadVariablesFromSocketIntoR function\n", i);
	  _EXIT_READ_VARIABLES;
	}
      }
      else if(datatype == VARTOOLS_TYPE_CHAR) {
	if(read(c->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	  fprintf(stderr,"Error reading length of vector for expected variable index %d in the ReadVariablesFromSocketIntoR function.\n", i);
	  _EXIT_READ_VARIABLES;
	}
	if(lenvec <= 0) {
	  continue;
	}
	if(lenvec > sizetmpinpstr) {
	  if(!sizetmpinpstr) {
	    if((tmpinpstr = (char *) malloc(lenvec * sizeof(char))) == NULL) {
	      fprintf(stderr,"Memory Allocation Error in VARTOOLS R sub-process.\n");
	      _EXIT_READ_VARIABLES;
	    }
	    sizetmpinpstr = lenvec;
	  } else {
	    if((tmpinpstr = (char *) realloc(tmpinpstr, lenvec * sizeof(char))) == NULL) {
	      fprintf(stderr,"Memory Allocation Error in VARTOOLS R sub-process.\n");
	      _EXIT_READ_VARIABLES;
	    }
	    sizetmpinpstr = lenvec;
	  }
	}
	databytesize = sizeof(char);
	sizedatareadtot = 0;
	while(sizedatareadtot < (((size_t) lenvec)*databytesize)) {
	  sizedataread = read(c->sockets[threadindex][1], (void *) (((char *) (tmpinpstr)) + sizedatareadtot), (((size_t) lenvec)*databytesize) - sizedatareadtot);
	  sizedatareadtot += sizedataread;
	}
	if(sizedatareadtot < (((size_t) lenvec)*databytesize)) {
	  fprintf(stderr,"Error reading data for variable index %d in the ReadVariablesFromSocketIntoR function.\n", i);
	  _EXIT_READ_VARIABLES;
	}
	if(sizetmpinpstr2 < 2) {
	  if(!sizetmpinpstr2) {
	    if((tmpinpstr2 = (char *) malloc(2*sizeof(char))) == NULL) {
	      fprintf(stderr,"Memory Allocation Error in VARTOOLS R sub-process.\n");
	      _EXIT_READ_VARIABLES;
	    }
	  } else {
	    if((tmpinpstr2 = (char *) realloc(tmpinpstr2,2*sizeof(char))) == NULL) {
	      fprintf(stderr,"Memory Allocation Error in VARTOOLS R sub-process.\n");
	      _EXIT_READ_VARIABLES;
	    }
	  }
	}
	Rcontainer->data[i].datatype = datatype;
	InputVariable = PROTECT(allocVector(STRSXP, lenvec));
	Nprotected++;
	for(j=0; j < lenvec; j++) {
	  tmpinpstr2[0] = tmpinpstr[j];
	  tmpinpstr2[1] = '\0';
	  SET_STRING_ELT(InputVariable, j, mkChar(tmpinpstr2));
	}
      } else if(datatype == VARTOOLS_TYPE_STRING) {
	if(read(c->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	  fprintf(stderr,"Error reading length of vector for expected variable index %d in the ReadVariablesFromSocketIntoR function.\n", i);
	  _EXIT_READ_VARIABLES;
	}
	if(lenvec <= 0) {
	  continue;
	}
	Rcontainer->data[i].datatype = datatype;
	InputVariable = PROTECT(allocVector(STRSXP, lenvec));
	Nprotected++;
	for(j=0; j < lenvec; j++) {
	  if(read(c->sockets[threadindex][1], &k, sizeof(int)) < sizeof(int)) {
	    fprintf(stderr,"Error reading length of string for vector %d, item %d in the ReadVariablesFromSocketIntoR function.\n", i, j);
	    _EXIT_READ_VARIABLES;
	  }
	  if(k < 0) k = 0;
	  if(k + 1 > sizetmpinpstr) {
	    if(!sizetmpinpstr) {
	      if((tmpinpstr = (char *) malloc((k+1) * sizeof(char))) == NULL) {
		fprintf(stderr,"Memory Allocation Error in VARTOOLS R sub-process.\n");
		_EXIT_READ_VARIABLES;
	      }
	      sizetmpinpstr = k+1;
	    } else {
	      if((tmpinpstr = (char *) realloc(tmpinpstr, (k+1) * sizeof(char))) == NULL) {
		fprintf(stderr,"Memory Allocation Error in VARTOOLS R sub-process.\n");
		_EXIT_READ_VARIABLES;
	      }
	      sizetmpinpstr = k+1;
	    }
	  }
	  if(k > 0) {

	    sizedatareadtot = 0;
	    while(sizedatareadtot < (((size_t) k)*sizeof(char))) {
	      sizedataread = read(c->sockets[threadindex][1], (void *) (((char *) (tmpinpstr)) + sizedatareadtot), (((size_t) k)*sizeof(char)) - sizedatareadtot);
	      sizedatareadtot += sizedataread;
	    }
	    if(sizedatareadtot < (((size_t)k)*sizeof(char))) {
	      fprintf(stderr,"Error reading length of string for vector %d, item %d in the ReadVariablesFromSocketIntoR function.\n", i, j);
	      _EXIT_READ_VARIABLES;
	    }
	  }
	  tmpinpstr[k] = '\0';
	  SET_STRING_ELT(InputVariable, j, mkChar(tmpinpstr));
	}
      }
      /* Transfer the variable to the list */
      SET_VECTOR_ELT(VariableList, i, InputVariable);
      UNPROTECT(1);
      Nprotected--;
    }
  }
  if(tmpinpstr != NULL)
    free(tmpinpstr);
  if(tmpinpstr2 != NULL)
    free(tmpinpstr2);
  if(Nprotected > 0)
    UNPROTECT(Nprotected);
  return VariableList;
}

#undef _EXIT_READ_VARIABLES

SEXP ReadVariablesFromSocketIntoR_all_lcs(ProgramData *p, 
					      _RCommand *c, 
					      int threadindex,
					  int *statusval) {
  int Nlcs, i;

  _RObjectContainer *Rcontainer = NULL;
  int Nprotected = 0;

  SEXP FullList;
  SEXP VariableList;
  int status2;

  *statusval = 0;

  if(read(c->sockets[threadindex][1], &Nlcs, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Error reading the number of light curves in the ReadVariablesFromSocketIntoR_all_lcs function.\n");
    *statusval = 1;
    return FullList;
  }

  /* The only argument is a list of input variables; First create the list */
  FullList = PROTECT(allocVector(VECSXP, Nlcs));
  Nprotected++;
  
  for(i=0; i < Nlcs; i++) {
    VariableList = PROTECT(ReadVariablesFromSocketIntoR(p, c, threadindex, 1, i, &status2));
    Nprotected++;
    if(status2) {
      UNPROTECT(Nprotected);
      *statusval = 1;
      return FullList;
    }
    SET_VECTOR_ELT(FullList,  i, VariableList);
  }
  
  UNPROTECT(Nprotected);
  return FullList;

}

int WriteVariablesFromRToSocket(ProgramData *p, _RCommand *cparent, 
				int threadindex, int cid, SEXP VariableListOut) {
  _RObjectContainer *Rcontainer = NULL;
  int i, j, k, ii;
  SEXP tmparray;
  SEXP tmplist;
  SEXP tmpcopyarray;

  int outdims;

  double tmpdblout;
  float tmpfloatout;
  int tmpintout;
  long tmplongout;
  short tmpshortout;
  int lenvec;

  char *tmpstr;
  char *tmpcharvec = NULL;
  int sizetmpcharvec = 0;

  double *dblptrout;
  float *floatptrout;
  int *intptrout;
  long *longptrout;
  short *shortptrout;
  
  _Variable *v;

  char *outstr;
  char outchar;
  int outstrlen;

  int tmpindx = 0;

  _RCommand *c;
  int Nprotected = 0;

  if(!cid) c = cparent;
  else c = ((_RCommand **)cparent->childcommandptrs)[cid-1];

  Rcontainer = ((_RObjectContainer **) cparent->Robjects)[threadindex];

  i = 0;
  j = 0;
  ii = 0;
  while(!(j == 1 && i >= c->Nvars_outonly)) {
    if(!j && i >= c->Nvars) { 
      j = 1;
      i = 0;
      if(i >= c->Nvars_outonly) break;
    }
    if(!j) {
      if(!c->isvaroutput[i]) {
	i++; ii++;
	continue;
      }
      v = c->vars[i];
    } else {
      v = c->outonlyvars[i];
    }
    if(v->datatype != VARTOOLS_TYPE_CHAR &&
       v->datatype != VARTOOLS_TYPE_STRING) {
      tmparray = PROTECT(VECTOR_ELT(VariableListOut,ii));
      Nprotected++;

      outdims = length(tmparray);
      if(v->vectortype != VARTOOLS_VECTORTYPE_LC) {
	if(outdims != 1) {
	  fprintf(stderr,"Error the variable %s has the wrong dimension upon output from R. This variable is expected to be a scalar, but was instead found to have dimension %d\n",v->varname,(int) (outdims));
	  if(tmpcharvec != NULL) free(tmpcharvec);
	  UNPROTECT(Nprotected);
	  return 1;
	}
	switch(v->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	case VARTOOLS_TYPE_CONVERTJD:
	  lenvec = 1;
	  if(write(cparent->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    UNPROTECT(Nprotected);
	    return 1;
	  }
	  tmpdblout = REAL(tmparray)[0];
	  if(write(cparent->sockets[threadindex][1], &tmpdblout, sizeof(double)) < sizeof(double)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    UNPROTECT(Nprotected);
	    return 1;
	  }
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  lenvec = 1;
	  if(write(cparent->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    UNPROTECT(Nprotected);
	    return 1;
	  }
	  tmpfloatout = (float) REAL(tmparray)[0];
	  if(write(cparent->sockets[threadindex][1], &tmpfloatout, sizeof(float)) < sizeof(float)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    UNPROTECT(Nprotected);
	    return 1;
	  }
	  break;
	case VARTOOLS_TYPE_INT:
	  lenvec = 1;
	  if(write(cparent->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    UNPROTECT(Nprotected);
	    return 1;
	  }
	  tmpintout = INTEGER(tmparray)[0];
	  if(write(cparent->sockets[threadindex][1], &tmpintout, sizeof(int)) < sizeof(int)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    UNPROTECT(Nprotected);
	    return 1;
	  }
	  break;
	case VARTOOLS_TYPE_LONG:
	  lenvec = 1;
	  if(write(cparent->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    UNPROTECT(Nprotected);
	    return 1;
	  }
	  tmplongout = (long) INTEGER(tmparray)[0];
	  if(write(cparent->sockets[threadindex][1], &tmplongout, sizeof(long)) < sizeof(long)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    UNPROTECT(Nprotected);
	    return 1;
	  }
	  break;
	case VARTOOLS_TYPE_SHORT:
	  lenvec = 1;
	  if(write(cparent->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    UNPROTECT(Nprotected);
	    return 1;
	  }
	  tmpshortout = (short) INTEGER(tmparray)[0];
	  if(write(cparent->sockets[threadindex][1], &tmpshortout, sizeof(short)) < sizeof(short)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    UNPROTECT(Nprotected);
	    return 1;
	  }
	  break;
	}
      } else {
	lenvec = outdims;
	if(write(cparent->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	  if(tmpcharvec != NULL) free(tmpcharvec);
	  UNPROTECT(Nprotected);
	  return 1;
	}
	if(lenvec <= 0) continue;
	switch(v->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	case VARTOOLS_TYPE_CONVERTJD:
	  dblptrout = REAL(tmparray);
	  if(write(cparent->sockets[threadindex][1], (void *) dblptrout, (((size_t) lenvec)*(sizeof(double)))) < (((size_t) lenvec)*(sizeof(double)))) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    UNPROTECT(Nprotected);
	    return 1;
	  }
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  if((floatptrout = malloc(lenvec*sizeof(float))) == NULL) {
	    DO_ERROR_MEMALLOC;
	  }
	  for(k=0; k < lenvec; k++) {
	    floatptrout[k] = (float) REAL(tmparray)[k];
	  }
	  if(write(cparent->sockets[threadindex][1], (void *) floatptrout, (((size_t) lenvec)*(sizeof(float)))) < (((size_t) lenvec)*(sizeof(float)))) {
	    free(floatptrout);
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    UNPROTECT(Nprotected);
	    return 1;
	  }
	  free(floatptrout);
	  break;
	case VARTOOLS_TYPE_INT:
	  intptrout = INTEGER(tmparray);
	  if(write(cparent->sockets[threadindex][1], (void *) intptrout, (((size_t) lenvec)*(sizeof(int)))) < (((size_t) lenvec)*(sizeof(int)))) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    UNPROTECT(Nprotected);
	    return 1;
	  }
	  break;
	case VARTOOLS_TYPE_LONG:
	  if((longptrout = malloc(lenvec*sizeof(long))) == NULL) {
	    DO_ERROR_MEMALLOC;
	  }
	  for(k=0; k < lenvec; k++) {
	    longptrout[k] = (long) INTEGER(tmparray)[k];
	  }
	  if(write(cparent->sockets[threadindex][1], (void *) longptrout, (((size_t) lenvec)*(sizeof(long)))) < (((size_t) lenvec)*(sizeof(long)))) {
	    free(longptrout);
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    UNPROTECT(Nprotected);
	    return 1;
	  }
	  free(longptrout);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  if((shortptrout = malloc(lenvec*sizeof(short))) == NULL) {
	    DO_ERROR_MEMALLOC;
	  }
	  for(k=0; k < lenvec; k++) {
	    shortptrout[k] = (short) INTEGER(tmparray)[k];
	  }
	  if(write(cparent->sockets[threadindex][1], (void *) shortptrout, (((size_t) lenvec)*(sizeof(short)))) < (((size_t) lenvec)*(sizeof(short)))) {
	    free(shortptrout);
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    UNPROTECT(Nprotected);
	    return 1;
	  }
	  free(shortptrout);
	  break;
	}
      }
      UNPROTECT(1);
      Nprotected--;
    } else {
      tmparray = PROTECT(VECTOR_ELT(VariableListOut,ii));
      Nprotected++;
      outdims = length(tmparray);
      lenvec = outdims;
      if(v->vectortype != VARTOOLS_VECTORTYPE_LC && lenvec != 1) {
	if(lenvec != 1) {
	  fprintf(stderr,"Error the variable %s has the wrong dimension upon output from R. This variable is expected to be a scalar, but was instead found to have dimension %d\n",v->varname,lenvec);
	  if(tmpcharvec != NULL) free(tmpcharvec);
	  UNPROTECT(Nprotected);
	  return 1;
	}
      }
      if(write(cparent->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	if(tmpcharvec != NULL) free(tmpcharvec);
	UNPROTECT(Nprotected);
	return 1;
      }
      for(k=0; k < lenvec; k++) {
	tmpstr = CHAR(STRING_ELT(tmparray, k));
	if(tmpstr == NULL) {
	  fprintf(stderr,"Error: NULL value returned from R for string or ASCII char variable %s\n", v->varname);
	  if(tmpcharvec != NULL) free(tmpcharvec);
	  UNPROTECT(Nprotected);
	  return 1;
	}
	if(v->datatype == VARTOOLS_TYPE_CHAR) {
	  if(lenvec <= 1) {
	    outchar = tmpstr[0];
	    if(write(cparent->sockets[threadindex][1], &outchar, sizeof(char)) < sizeof(char)) {
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      UNPROTECT(Nprotected);
	      return 1;
	    }
	  } else {
	    if(!k) {
	      if(lenvec > sizetmpcharvec) {
		if(!sizetmpcharvec) {
		  if((tmpcharvec = (char *) malloc(lenvec * sizeof(char))) == NULL)
		    DO_ERROR_MEMALLOC;
		} else {
		  if((tmpcharvec = (char *) realloc(tmpcharvec, lenvec * sizeof(char))) == NULL)
		    DO_ERROR_MEMALLOC;
		}
		sizetmpcharvec = lenvec;
	      }
	    }
	    tmpcharvec[k] = tmpstr[0];
	    if(k == lenvec - 1) {
	      if(write(cparent->sockets[threadindex][1], tmpcharvec, (((size_t) lenvec)*sizeof(char))) < (((size_t) lenvec)*sizeof(char))) {
		if(tmpcharvec != NULL) free(tmpcharvec);
		UNPROTECT(Nprotected);
		return 1;
	      }
	    }
	  }
	} else {
	  outstrlen = strlen(tmpstr);
	  if(write(cparent->sockets[threadindex][1], &outstrlen, sizeof(int)) < sizeof(int)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    UNPROTECT(Nprotected);
	    return 1;
	  }
	  if(outstrlen > 0) {
	    if(write(cparent->sockets[threadindex][1], tmpstr, (((size_t) (outstrlen))*sizeof(char))) < (((size_t) (outstrlen))*sizeof(char))) {
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      UNPROTECT(Nprotected);
	      return 1;
	    }
	  }
	}
      }
      UNPROTECT(1);
      Nprotected--;
    }
    i++; ii++;
  }
  /* The data has all been transferred back to the parent process; we can
     clean-up the R structures associated with this particular light
     curve */
  CleanUpRObjectContainerVariables(Rcontainer);
  

  if(Nprotected > 0) UNPROTECT(Nprotected);
  if(tmpcharvec != NULL) free(tmpcharvec);
  return 0;
}

int WriteVariablesFromRToSocket_all_lcs(ProgramData *p, 
					     _RCommand *cparent, 
					     int threadindex, int cid, SEXP FullList) {
  int Nlcs, i;
  _RObjectContainer *Rcontainer = NULL;

  _RCommand *c;
  int Nprotected = 0;
  SEXP VariableListOut;

  if(!cid) c = cparent;
  else c = ((_RCommand **)cparent->childcommandptrs)[cid-1];
  Rcontainer = ((_RObjectContainer **) cparent->Robjects)[threadindex];

  Nlcs = length(FullList);

  if(write(cparent->sockets[threadindex][1], &Nlcs, sizeof(int)) < sizeof(int))
    return 1;

  for(i=0; i < Nlcs; i++) {
    VariableListOut = PROTECT(VECTOR_ELT(FullList,i));
    if(WriteVariablesFromRToSocket(p, cparent, threadindex, cid, VariableListOut)) {
      UNPROTECT(1);
      return 1;
    }
    UNPROTECT(1);
  }
  
  return 0;
}

SEXP RunRUserFunctionOnLightCurve(ProgramData *p, _RCommand *c, 
				  int threadindex, int cid, int isall, int *statusval, SEXP ArgumentList) {
  int i, cnumval;
  _RObjectContainer *Rcontainer;
  int Nprotected = 0;
  SEXP function_call;
  SEXP OutputList;
  
  char functionnametext[256];

  *statusval = 0;

  if(!cid) 
    cnumval = c->cnum;
  else
    cnumval = ((_RCommand **)c->childcommandptrs)[cid-1]->cnum;

  sprintf(functionnametext,"VARTOOLS_R_USERFUNCTION_%d",cnumval);

  function_call = PROTECT(lang2(install(functionnametext), ArgumentList));


  OutputList = PROTECT(R_tryEval(function_call, R_GlobalEnv, statusval));
  
  if(*statusval) {
    fprintf(stderr,"Error running R command.\n");
  }

  UNPROTECT(2);

  return(OutputList);
}

void TerminateR(ProgramData *p, _RCommand *c, int threadindex)
{
  if(!c->IsRRunning[threadindex]) return;
  c->IsRRunning[threadindex] = 0;
  if(!c->iscontinueprocess) {
    CleanUpRObjectContainer(&(((_RObjectContainer **) c->Robjects)[threadindex]));
    Rf_endEmbeddedR(0);
  } else {
    TerminateR(p, (_RCommand *) c->continueprocesscommandptr, threadindex);
  }
}

void RunRProcessingLoop(ProgramData *p, _RCommand *c, int threadindex)
{
  int msg, retval = 0, newcnum, cid;
  _RCommand *ccheck;
  SEXP InputVariableList;
  SEXP OutputVariableList;
  int Nprotected = 0;

  while(1) {
    /* Read a message from the parent indicating what we will be doing */
    /* If no message is received, terminate the loop */
    if(read(c->sockets[threadindex][1], &msg, sizeof(int)) < sizeof(int))
      return;
  
    switch(msg) {
    case VARTOOLS_R_MESSAGE_ENDPROCESS:
      return;
    case VARTOOLS_R_MESSAGE_READDATA:

      if(read(c->sockets[threadindex][1], &cid, sizeof(int)) < sizeof(int))
	return;

      if(cid > 0 && (cid - 1 > c->Nchildren))
	return;

      if(!cid) ccheck = c;
      else ccheck = ((_RCommand **)c->childcommandptrs)[cid-1];


      if(ccheck->RequireReadAll) {
	InputVariableList = PROTECT(ReadVariablesFromSocketIntoR_all_lcs(p, c, threadindex, &retval));
	if(retval) {
	  UNPROTECT(1);
	  if(write(c->sockets[threadindex][1], &retval, sizeof(int)) < sizeof(int))
	    return;
	  break;
	}
	OutputVariableList = PROTECT(RunRUserFunctionOnLightCurve(p, c, threadindex, cid, 1, &retval, InputVariableList));
	if(retval) {
	  UNPROTECT(2);
	  if(write(c->sockets[threadindex][1], &retval, sizeof(int)) < sizeof(int))
	    return;
	  break;
	}
      } else {
	InputVariableList = PROTECT(ReadVariablesFromSocketIntoR(p, c, threadindex, 0, 0, &retval));
	if(retval) {
	  UNPROTECT(1);
	  if(write(c->sockets[threadindex][1], &retval, sizeof(int)) < sizeof(int))
	    return;
	  break;
	}
	OutputVariableList = PROTECT(RunRUserFunctionOnLightCurve(p, c, threadindex, cid, 0, &retval, InputVariableList));
	if(retval) {
	  UNPROTECT(2);
	  if(write(c->sockets[threadindex][1], &retval, sizeof(int)) < sizeof(int))
	    return;
	  break;
	}
      }
      if(write(c->sockets[threadindex][1], &retval, sizeof(int)) < sizeof(int)) {
	UNPROTECT(2);
	return;
      }
      if(ccheck->RequireReadAll) {
	if(WriteVariablesFromRToSocket_all_lcs(p, c, threadindex, cid, OutputVariableList)) {
	  UNPROTECT(2);
	  return;
	}
      } else {
	if(WriteVariablesFromRToSocket(p, c, threadindex, cid, OutputVariableList)) {
	  UNPROTECT(2);
	  return;
	}
      }
      UNPROTECT(2);
      break;
    default:
      fprintf(stderr,"Error: Invalid message received by R sub-process from main vartools process. Shutting down the R sub-process now.\n");
      return;
    }
  }
}


void StartRProcess(ProgramData *p, _RCommand *c, int threadindex)
{
  pid_t pid;

  /* If R is already running on this thread, just return */
  if(c->IsRRunning[threadindex]) return;
  
  c->IsRRunning[threadindex] = 1;

  /* Create a socket pair */
  if(socketpair(AF_UNIX, SOCK_STREAM, 0, c->sockets[threadindex]) < 0) {
    perror("Error Opening socket pair in starting R process\n");
    exit(1);
  }

  /* Create the sub-process */
  pid = fork();
  if(pid == -1) {
    perror("Error creating a sub-process initiating a call to the vartools R command.\n");
    exit(1);
  }
  if(!pid) {
    /* This is the child sub-process which will be used for running R */

    /* Close the parent end of the socket */
    close(c->sockets[threadindex][0]);    

    /* Redirect stdout to stderr */
    dup2(2, 1);

    /* Start R, load all of the
       initialization codes and create the function that will be used for
       processing individual light curves */
    InitializeR(p, c, threadindex);

    /* The loop that will receive light curve data from the parent process,
       run the R routine, and write the data back to the parent */
    RunRProcessingLoop(p, c, threadindex);

    /* If the loop is done running we can terminate this process */
    TerminateR(p, c, threadindex);

    /* Close the socket, then terminate the process */
    close(c->sockets[threadindex][1]);
    exit(0);
  } else {

    /* This is the main parent loop. In this case we close the child end of 
       the socket, and then continue with the processing */
    close(c->sockets[threadindex][1]);
  }
}


int SendVariablesToChildRProcess(ProgramData *p, int lcindex,
				      int threadindex, int Rthreadindex, _RCommand *c)
{
  
  int i, k;
  size_t sizedatasent;
  size_t databytesize;
  int lenvec;
  int tmpindex;
  void *ptrtosend;
  double tmpdblval;
  float tmpfloatval;
  int tmpintval;
  int strlentosend;
  long tmplongval;
  short tmpshortval;
  char tmpcharval;
  char *tmpstringval = NULL;
  int sizetmpstringval = 0;

  if(write(c->sockets[Rthreadindex][0], &c->Nvars, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Error sending the number of variables to read to the child R process.\n");
    if(tmpstringval != NULL) free(tmpstringval);
    return 1;
  }

  for(i=0; i < c->Nvars; i++) {
    if(write(c->sockets[Rthreadindex][0], &(c->vars[i]->datatype), sizeof(char)) < sizeof(char)) {
      fprintf(stderr,"Error sending the datatype for variable %s to the child R process\n", c->vars[i]->varname);
      if(tmpstringval != NULL) free(tmpstringval);
      return 1;
    }
    if(c->vars[i]->datatype != VARTOOLS_TYPE_STRING) {
      switch(c->vars[i]->datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	databytesize = sizeof(double);
	break;
      case VARTOOLS_TYPE_CONVERTJD:
	databytesize = sizeof(double);
	break;
      case VARTOOLS_TYPE_FLOAT:
	databytesize = sizeof(float);
	break;
      case VARTOOLS_TYPE_INT:
	databytesize = sizeof(int);
	break;
      case VARTOOLS_TYPE_LONG:
	databytesize = sizeof(long);
	break;
      case VARTOOLS_TYPE_SHORT:
	databytesize = sizeof(short);
	break;
      case VARTOOLS_TYPE_CHAR:
	databytesize = sizeof(char);
	break;
      default:
	fprintf(stderr,"Error: invalid datatype received for variable %s in the SendVariablesToChildRProcess function\n", c->vars[i]->varname);
	if(tmpstringval != NULL) free(tmpstringval);
	return 1;
      }
      if(c->vars[i]->vectortype != VARTOOLS_VECTORTYPE_LC)
	lenvec = 1;
      else
	lenvec = p->NJD[threadindex];
      
      if(write(c->sockets[Rthreadindex][0], &lenvec, sizeof(int)) < sizeof(int)) {
	fprintf(stderr,"Error sending the length of vector for variable %s in the SendVariablesToChildRProcess function.\n", c->vars[i]->varname);
	if(tmpstringval != NULL) free(tmpstringval);
	return 1;
      }

      if(lenvec <= 0) {
	continue;
      }

      switch(c->vars[i]->vectortype) {
      case VARTOOLS_VECTORTYPE_CONSTANT:
	ptrtosend = c->vars[i]->dataptr;
	break;
      case VARTOOLS_VECTORTYPE_SCALAR:
      case VARTOOLS_VECTORTYPE_INLIST:
	if(c->vars[i]->vectortype == VARTOOLS_VECTORTYPE_SCALAR)
	  tmpindex = threadindex;
	else
	  tmpindex = lcindex;
	switch(c->vars[i]->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	case VARTOOLS_TYPE_CONVERTJD:
	  ptrtosend = (void *) &((*((double **) c->vars[i]->dataptr))[tmpindex]);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  ptrtosend = (void *) &((*((float **) c->vars[i]->dataptr))[tmpindex]);
	  break;
	case VARTOOLS_TYPE_INT:
	  ptrtosend = (void *) &((*((int **) c->vars[i]->dataptr))[tmpindex]);
	  break;
	case VARTOOLS_TYPE_LONG:
	  ptrtosend = (void *) &((*((long **) c->vars[i]->dataptr))[tmpindex]);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  ptrtosend = (void *) &((*((short **) c->vars[i]->dataptr))[tmpindex]);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  ptrtosend = (void *) &((*((char **) c->vars[i]->dataptr))[tmpindex]);
	  break;
	default:
	  ptrtosend = NULL;
	  break;
	}
	break;
      case VARTOOLS_VECTORTYPE_OUTCOLUMN:
	switch(c->vars[i]->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	case VARTOOLS_TYPE_CONVERTJD:
	  getoutcolumnvalue(c->vars[i]->outc, threadindex, lcindex, VARTOOLS_TYPE_DOUBLE, &tmpdblval);
	  ptrtosend = (void *) &tmpdblval;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  getoutcolumnvalue(c->vars[i]->outc, threadindex, lcindex, VARTOOLS_TYPE_FLOAT, &tmpfloatval);
	  ptrtosend = (void *) &tmpfloatval;
	  break;
	case VARTOOLS_TYPE_INT:
	  getoutcolumnvalue(c->vars[i]->outc, threadindex, lcindex, VARTOOLS_TYPE_INT, &tmpintval);
	  ptrtosend = (void *) &tmpintval;
	  break;
	case VARTOOLS_TYPE_LONG:
	  getoutcolumnvalue(c->vars[i]->outc, threadindex, lcindex, VARTOOLS_TYPE_LONG, &tmplongval);
	  ptrtosend = (void *) &tmplongval;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  getoutcolumnvalue(c->vars[i]->outc, threadindex, lcindex, VARTOOLS_TYPE_SHORT, &tmpshortval);
	  ptrtosend = (void *) &tmpshortval;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  getoutcolumnvalue(c->vars[i]->outc, threadindex, lcindex, VARTOOLS_TYPE_CHAR, &tmpcharval);
	  ptrtosend = (void *) &tmpcharval;
	  break;
	default:
	  ptrtosend = NULL;
	  break;
	}
	break;
      case VARTOOLS_VECTORTYPE_LC:
	switch(c->vars[i]->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	case VARTOOLS_TYPE_CONVERTJD:
	  ptrtosend = (void *) &((*((double ***) c->vars[i]->dataptr))[threadindex][0]);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  ptrtosend = (void *) &((*((float ***) c->vars[i]->dataptr))[threadindex][0]);
	  break;
	case VARTOOLS_TYPE_INT:
	  ptrtosend = (void *) &((*((int ***) c->vars[i]->dataptr))[threadindex][0]);
	  break;
	case VARTOOLS_TYPE_LONG:
	  ptrtosend = (void *) &((*((long ***) c->vars[i]->dataptr))[threadindex][0]);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  ptrtosend = (void *) &((*((short ***) c->vars[i]->dataptr))[threadindex][0]);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  ptrtosend = (void *) &((*((char ***) c->vars[i]->dataptr))[threadindex][0]);
	  break;
	default:
	  ptrtosend = NULL;
	  break;
	}
	break;
      default:
	ptrtosend = NULL;
	break;
      }

      sizedatasent = write(c->sockets[Rthreadindex][0], ptrtosend, (((size_t) lenvec)*databytesize));
      if(sizedatasent < (((size_t) lenvec)*databytesize)) {
	fprintf(stderr,"Error sending data for variable %s in the SendVariablesToChildRProcess function.\n", c->vars[i]->varname);
	if(tmpstringval != NULL) free(tmpstringval);
	return 1;
      }
    } else {
      databytesize = sizeof(char);
      if(c->vars[i]->vectortype != VARTOOLS_VECTORTYPE_LC)
	lenvec = 1;
      else
	lenvec = p->NJD[threadindex];
      
      if(write(c->sockets[Rthreadindex][0], &lenvec, sizeof(int)) < sizeof(int)) {
	fprintf(stderr,"Error sending the length of vector for variable %s in the SendVariablesToChildRProcess function.\n", c->vars[i]->varname);
	if(tmpstringval != NULL) free(tmpstringval);
	return 1;
      }
      if(lenvec <= 0) {
	continue;
      }
      if(c->vars[i]->vectortype != VARTOOLS_VECTORTYPE_LC) {
	switch(c->vars[i]->vectortype) {
	case VARTOOLS_VECTORTYPE_CONSTANT:
	  ptrtosend = (void *) &((*((char **)c->vars[i]->dataptr))[0]);
	  strlentosend = strlen((char *)ptrtosend);
	  break;
	case VARTOOLS_VECTORTYPE_SCALAR:
	case VARTOOLS_VECTORTYPE_INLIST:
	  if(c->vars[i]->vectortype == VARTOOLS_VECTORTYPE_SCALAR)
	    tmpindex = threadindex;
	  else
	    tmpindex = lcindex;
	  ptrtosend = (void *) &((*((char ***) c->vars[i]->dataptr))[tmpindex][0]);
	  strlentosend = strlen((char *)ptrtosend);
	  break;
	case VARTOOLS_VECTORTYPE_OUTCOLUMN:
	  if(c->vars[i]->outc->stringsize > sizetmpstringval) {
	    if(!sizetmpstringval) {
	      if((tmpstringval = (char *) malloc(c->vars[i]->outc->stringsize)) == NULL)
		DO_ERROR_MEMALLOC;
	    } else {
	      if((tmpstringval = (char *) realloc(tmpstringval, c->vars[i]->outc->stringsize)) == NULL)
		DO_ERROR_MEMALLOC;
	    }
	    sizetmpstringval = c->vars[i]->outc->stringsize;
	  }
	  getoutcolumnvalue(c->vars[i]->outc, threadindex, lcindex, VARTOOLS_TYPE_STRING, tmpstringval,c->vars[i]->outc->stringsize);
	  ptrtosend = tmpstringval;
	  strlentosend = strlen((char *)ptrtosend);
	  break;
	default:
	  ptrtosend = NULL;
	  strlentosend = 0;
	}
	if(write(c->sockets[Rthreadindex][0], &strlentosend, sizeof(int)) < sizeof(int)) {
	  fprintf(stderr,"Error sending the string length for variable %s in the SendVariablesToChildRProcess function.\n", c->vars[i]->varname);
	  if(tmpstringval != NULL) free(tmpstringval);
	  return 1;
	}
	if(strlentosend <= 0) continue;
	if(write(c->sockets[Rthreadindex][0], ptrtosend, (((size_t) strlentosend)*sizeof(char))) < (((size_t) strlentosend)*sizeof(char))) {
	  fprintf(stderr,"Error sending the data for variable %s in the SendVariablesToChildRProcess function.\n", c->vars[i]->varname);
	  if(tmpstringval != NULL) free(tmpstringval);
	  return 1;
	}
      } else {
	for(k=0; k < lenvec; k++) {
	  ptrtosend = (void *) &((*((char ****) c->vars[i]->dataptr))[threadindex][k][0]);
	  strlentosend = strlen((char *) ptrtosend);
	  if(write(c->sockets[Rthreadindex][0], &strlentosend, sizeof(int)) < sizeof(int)) {
	    fprintf(stderr,"Error sending the string length for variable %s in the SendVariablesToChildRProcess function.\n", c->vars[i]->varname);
	    if(tmpstringval != NULL) free(tmpstringval);
	    return 1;
	  }
	  if(strlentosend <= 0) continue;
	  if(write(c->sockets[Rthreadindex][0], ptrtosend, (((size_t) strlentosend)*sizeof(char))) < (((size_t) strlentosend)*sizeof(char))) {
	    fprintf(stderr,"Error sending the data for variable %s in the SendVariablesToChildRProcess function.\n", c->vars[i]->varname);
	    if(tmpstringval != NULL) free(tmpstringval);
	    return 1;
	  }
	}
      }
    }
  }
  if(tmpstringval != NULL) free(tmpstringval);
  return 0;
}

int SendVariablesToChildRProcess_all_lcs(ProgramData *p, _RCommand *c)
{
  int j, Nlcs;

  Nlcs = p->Nlcs;
  if(write(c->sockets[0][0], &Nlcs, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Error sending the number of light curves to read to the child R process.\n");
    return 1;
  }
  for(j=0; j < Nlcs; j++) {
    if(SendVariablesToChildRProcess(p, j, j, 0, c))
      return 1;
  }
  return 0;
}

int ReadVariablesFromChildRProcess(ProgramData *p, int lcindex, 
					int threadindex, int Rthreadindex, 
					_RCommand *c)
{
  int i, j, k, ii, tmpindex;

  double tmpdblout;
  float tmpfloatout;
  int tmpintout;
  long tmplongout;
  short tmpshortout;
  int lenvec;

  double *dblptrout;
  float *floatptrout;
  int *intptrout;
  long *longptrout;
  short *shortptrout;
  
  _Variable *v;

  int *outlcvecptr;

  void *ptrtoget;

  char *outstr;
  char outchar;
  int instrlen;

  size_t databytesize;

  char *tmpinstr = NULL;
  int lentmpinstr = 0;

  int islclengthchange = 0;
  int maxlcoutlength = 0;
  double oldNJD;

  size_t sizedatareadtot;
  size_t sizedataread;

  i = 0;
  j = 0;
  while(!(j == 1 && i >= c->Nvars_outonly)) {
    if(!j && i >= c->Nvars) { 
      j = 1;
      i = 0;
      if(i >= c->Nvars_outonly) break;
    }
    if(!j) {
      if(!c->isvaroutput[i]) {
	i++;
	continue;
      }
      v = c->vars[i];
      outlcvecptr = &(c->outlcvecs_invars[threadindex][i]);
    } else {
      v = c->outonlyvars[i];
      outlcvecptr = &(c->outlcvecs_outonlyvars[threadindex][i]);
    }

    if(v->vectortype == VARTOOLS_VECTORTYPE_CONSTANT ||
       v->vectortype == VARTOOLS_VECTORTYPE_OUTCOLUMN) {
      fprintf(stderr,"Error: somehow a constant or outcolumn variable is set to receive data from a R command. This should not happen and indicates a bug in vartools.\n");
      exit(1);
    }
    
    if(read(c->sockets[Rthreadindex][0], &lenvec, sizeof(int)) < sizeof(int)) {
      fprintf(stderr,"Error receiving the vector length for variable %s in the ReadVariablesFromChildRProcess function.\n", v->varname);
      if(tmpinstr != NULL) free(tmpinstr);
      return 1;
    }
    
    if(v->datatype != VARTOOLS_TYPE_STRING) {
      
      if(v->vectortype != VARTOOLS_VECTORTYPE_LC) {
	if(lenvec != 1) {
	  fprintf(stderr,"Error the vector length received for variable %s has the wrong dimension in the ReadVariablesFromChildRProcess function. Expected only 1 value, received a length of %d\n",v->varname,lenvec);
	  if(tmpinstr != NULL) free(tmpinstr);
	  return 1;
	}
	if(v->vectortype == VARTOOLS_VECTORTYPE_SCALAR)
	  tmpindex = threadindex;
	else
	  tmpindex = lcindex;
	switch(v->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	case VARTOOLS_TYPE_CONVERTJD:
	  databytesize = sizeof(double);
	  ptrtoget = (void *) &((*((double **) v->dataptr))[tmpindex]);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  databytesize = sizeof(float);
	  ptrtoget = (void *) &((*((float **) v->dataptr))[tmpindex]);
	  break;
	case VARTOOLS_TYPE_INT:
	  databytesize = sizeof(int);
	  ptrtoget = (void *) &((*((int **) v->dataptr))[tmpindex]);
	  break;
	case VARTOOLS_TYPE_LONG:
	  databytesize = sizeof(long);
	  ptrtoget = (void *) &((*((long **) v->dataptr))[tmpindex]);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  databytesize = sizeof(short);
	  ptrtoget = (void *) &((*((short **) v->dataptr))[tmpindex]);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  databytesize = sizeof(char);
	  ptrtoget = (void *) &((*((char **) v->dataptr))[tmpindex]);
	  break;
	default:
	  ptrtoget = NULL;
	  break;
	}
	if(read(c->sockets[Rthreadindex][0], ptrtoget, databytesize) < databytesize) {
	  fprintf(stderr,"Error receiving the data for variable %s in the ReadVariablesFromChildRProcess function.\n",v->varname);
	  if(tmpinstr != NULL) free(tmpinstr);
	  return 1;
	}
      } else {
	*outlcvecptr = lenvec;
	if(lenvec > p->sizesinglelc[threadindex]) {
	  /* Not enough space to store the new data in the existing
	     light curve vectors; We will need to grow the vectors; */
	  MemAllocDataFromLightCurveMidProcess(p, threadindex, lenvec);
	}
	if(lenvec > maxlcoutlength) maxlcoutlength = lenvec;
	if(lenvec != p->NJD[threadindex]) {
	  islclengthchange = 1;
	}
	switch(v->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	case VARTOOLS_TYPE_CONVERTJD:
	  databytesize = sizeof(double);
	  ptrtoget = (void *) &((*((double ***) v->dataptr))[threadindex][0]);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  databytesize = sizeof(float);
	  ptrtoget = (void *) &((*((float ***) v->dataptr))[threadindex][0]);
	  break;
	case VARTOOLS_TYPE_INT:
	  databytesize = sizeof(int);
	  ptrtoget = (void *) &((*((int ***) v->dataptr))[threadindex][0]);
	  break;
	case VARTOOLS_TYPE_LONG:
	  databytesize = sizeof(long);
	  ptrtoget = (void *) &((*((long ***) v->dataptr))[threadindex][0]);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  databytesize = sizeof(short);
	  ptrtoget = (void *) &((*((short ***) v->dataptr))[threadindex][0]);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  databytesize = sizeof(char);
	  ptrtoget = (void *) &((*((char ***) v->dataptr))[threadindex][0]);
	  break;
	default:
	  ptrtoget = NULL;
	  break;
	}
	if(lenvec <= 0) continue;
	sizedatareadtot = 0;
	while(sizedatareadtot < (((size_t) lenvec)*databytesize)) {
	  sizedataread = read(c->sockets[Rthreadindex][0], (void *) (((char *) ptrtoget) + sizedatareadtot), (((size_t)lenvec)*databytesize)-sizedatareadtot);
	  sizedatareadtot += sizedataread;
	}
	if(sizedatareadtot < (((size_t) lenvec)*databytesize)) {
	  fprintf(stderr,"Error receiving the data for variable %s in the ReadVariablesFromChildRProcess function.\n",v->varname);
	  if(tmpinstr != NULL) free(tmpinstr);
	  return 1;
	}
      }
    } else {
      if(v->vectortype != VARTOOLS_VECTORTYPE_LC) {
	if(lenvec != 1) {
	  fprintf(stderr,"Error the vector length received for variable %s has the wrong dimension in the ReadVariablesFromChildRProcess function. Expected only 1 value, received a length of %d\n",v->varname,lenvec);
	  if(tmpinstr != NULL) free(tmpinstr);
	  return 1;
	}
	if(read(c->sockets[Rthreadindex][0], &instrlen, sizeof(int)) < sizeof(int)) {
	  fprintf(stderr,"Error receiving the string length for variable %s in the ReadVariablesFromChildRProcess function.\n",v->varname);
	  if(tmpinstr != NULL) free(tmpinstr);
	  return 1;
	}

	if(v->vectortype == VARTOOLS_VECTORTYPE_SCALAR)
	  tmpindex = threadindex;
	else
	  tmpindex = lcindex;

	ptrtoget = (void *) &((*((char ***) v->dataptr))[tmpindex][0]);

	if(instrlen > 0) {
	  if(instrlen+1 >= MAXLEN) {
	    if(instrlen+1 > lentmpinstr) {
	      if(!lentmpinstr) {
		if((tmpinstr = (char *) malloc((instrlen+1)*sizeof(char))) == NULL)
		  DO_ERROR_MEMALLOC;
	      } else {
		if((tmpinstr = (char *) realloc(tmpinstr, (instrlen+1)*sizeof(char))) == NULL)
		  DO_ERROR_MEMALLOC;
	      }
	      lentmpinstr = instrlen+1;
	    }
	    sizedatareadtot = 0;
	    while(sizedatareadtot < (((size_t) instrlen)*sizeof(char))) {
	      sizedataread = read(c->sockets[Rthreadindex][0], (void *) (((char *) tmpinstr) + sizedatareadtot), (((size_t)instrlen)*sizeof(char))-sizedatareadtot);
	      sizedatareadtot += sizedataread;
	    }
	    if(sizedatareadtot < (((size_t) instrlen)*sizeof(char))) {
	      fprintf(stderr,"Error receiving data for variable %s in the ReadVariablesFromChildRProcess function.\n",v->varname);
	      if(tmpinstr != NULL) free(tmpinstr);
	      return 1;
	    }
	    for(ii=0; ii < MAXLEN-1; ii++)
	      ((char *)ptrtoget)[ii] = tmpinstr[ii];
	    ((char *)ptrtoget)[ii] = '\0';
	  } else {
	    sizedatareadtot = 0;
	    while(sizedatareadtot < (((size_t) instrlen)*sizeof(char))) {
	      sizedataread = read(c->sockets[Rthreadindex][0], (void *) (((char *) ptrtoget) + sizedatareadtot), (((size_t)instrlen)*sizeof(char))-sizedatareadtot);
	      sizedatareadtot += sizedataread;
	    }
	    if(sizedatareadtot < (((size_t) instrlen)*sizeof(char))) {
	      fprintf(stderr,"Error receiving data for variable %s in the ReadVariablesFromChildRProcess function.\n",v->varname);
	      if(tmpinstr != NULL) free(tmpinstr);
	      return 1;
	    }
	    ((char *)ptrtoget)[instrlen] = '\0';
	  }
	} else {
	  ((char *)ptrtoget)[0] = '\0';
	}
      } else {
	for(k=0; k < lenvec; k++) {
	  ptrtoget = (void *) &((*((char ****) v->dataptr))[threadindex][k][0]);
	  if(read(c->sockets[Rthreadindex][0], &instrlen, sizeof(int)) < sizeof(int)) {
	    fprintf(stderr,"Error receiving the string length for variable %s in the ReadVariablesFromChildRProcess function.\n",v->varname);
	    if(tmpinstr != NULL) free(tmpinstr);
	    return 1;
	  }
	  if(instrlen > 0) {
	    if(instrlen+1 >= MAXLEN) {
	      if(instrlen+1 > lentmpinstr) {
		if(!lentmpinstr) {
		  if((tmpinstr = (char *) malloc((instrlen+1)*sizeof(char))) == NULL)
		    DO_ERROR_MEMALLOC;
		} else {
		  if((tmpinstr = (char *) realloc(tmpinstr, (instrlen+1)*sizeof(char))) == NULL)
		    DO_ERROR_MEMALLOC;
		}
		lentmpinstr = instrlen+1;
	      }
	      sizedatareadtot = 0;
	      while(sizedatareadtot < (((size_t) instrlen)*sizeof(char))) {
		sizedataread = read(c->sockets[Rthreadindex][0], (void *) (((char *) tmpinstr) + sizedatareadtot), (((size_t)instrlen)*sizeof(char))-sizedatareadtot);
		sizedatareadtot += sizedataread;
	      }
	      if(sizedatareadtot < (((size_t) instrlen)*sizeof(char))) {
		fprintf(stderr,"Error receiving data for variable %s in the ReadVariablesFromChildRProcess function.\n",v->varname);
		if(tmpinstr != NULL) free(tmpinstr);
		return 1;
	      }
	      for(ii=0; ii < MAXLEN-1; ii++)
		((char *)ptrtoget)[ii] = tmpinstr[ii];
	      ((char *)ptrtoget)[ii] = '\0';
	    } else {
	      sizedatareadtot = 0;
	      while(sizedatareadtot < (((size_t) instrlen)*sizeof(char))) {
		sizedataread = read(c->sockets[Rthreadindex][0], (void *) (((char *) ptrtoget) + sizedatareadtot), (((size_t)instrlen)*sizeof(char))-sizedatareadtot);
		sizedatareadtot += sizedataread;
	      }
	      if(sizedatareadtot < (((size_t) instrlen)*sizeof(char))) {
		fprintf(stderr,"Error receiving data for variable %s in the ReadVariablesFromChildRProcess function.\n",v->varname);
		if(tmpinstr != NULL) free(tmpinstr);
		return 1;
	      }
	      ((char *)ptrtoget)[instrlen] = '\0';
	    }
	  } else {
	    ((char *)ptrtoget)[0] = '\0';
	  }
	}
      }
    }
    i++;
  }

  /* Check if the lengths of any of the lcvecs has changed */
  if(islclengthchange) {
    /* If any output lc vector is longer than the previous light curve size;
       set the new light curve size to this new length and then append 0 to
       any vector that is shorter than this length */
    if(maxlcoutlength > p->NJD[threadindex]) {
      oldNJD = p->NJD[threadindex];
      p->NJD[threadindex] = maxlcoutlength;
      /* Grow all vectors with a length less than this limit */
      for(i=0; i < c->Nvars; i++) {
	if(c->isvaroutput[i]) {
	  if(c->vars[i]->vectortype == VARTOOLS_VECTORTYPE_LC) {
	    if(c->outlcvecs_invars[threadindex][i] < maxlcoutlength) {
	      if(c->vars[i]->datatype != VARTOOLS_TYPE_STRING &&
		 c->vars[i]->datatype != VARTOOLS_TYPE_CHAR) {
		for(j=c->outlcvecs_invars[threadindex][i]; j < maxlcoutlength; j++) {
		  SetVariable_Value_Double(lcindex, threadindex, j, c->vars[i], 0.0);
		}
	      } else if(c->vars[i]->datatype == VARTOOLS_TYPE_CHAR) {
		for(j=c->outlcvecs_invars[threadindex][i]; j < maxlcoutlength; j++) {
		  (*((char ***) c->vars[i]->dataptr))[threadindex][j] = '0';
		}
	      } else if(c->vars[i]->datatype == VARTOOLS_TYPE_STRING) {
		for(j=c->outlcvecs_invars[threadindex][i]; j < maxlcoutlength; j++) {
		  sprintf((*((char ****) c->vars[i]->dataptr))[threadindex][j], "0");
		}
	      }
	    }
	  }
	}
      }
      for(i=0; i < c->Nvars_outonly; i++) {
	if(c->outonlyvars[i]->vectortype == VARTOOLS_VECTORTYPE_LC) {
	  if(c->outlcvecs_outonlyvars[threadindex][i] < maxlcoutlength) {
	    if(c->outonlyvars[i]->datatype != VARTOOLS_TYPE_STRING &&
	       c->outonlyvars[i]->datatype != VARTOOLS_TYPE_CHAR) {
	      for(j=c->outlcvecs_outonlyvars[threadindex][i]; j < maxlcoutlength; j++) {
		SetVariable_Value_Double(lcindex, threadindex, j, c->outonlyvars[i], 0.0);
	      }
	    } else if(c->outonlyvars[i]->datatype == VARTOOLS_TYPE_CHAR) {
	      for(j=c->outlcvecs_outonlyvars[threadindex][i]; j < maxlcoutlength; j++) {
		(*((char ***) c->outonlyvars[i]->dataptr))[threadindex][j] = '0';
	      }
	    } else if(c->outonlyvars[i]->datatype == VARTOOLS_TYPE_STRING) {
	      for(j=c->outlcvecs_outonlyvars[threadindex][i]; j < maxlcoutlength; j++) {
		sprintf((*((char ****) c->outonlyvars[i]->dataptr))[threadindex][j], "0");
	      }
	    }
	  }
	}
      }
      for(i=0; i < c->Nlcvars_nonupdate; i++) {
	if(c->lcvars_nonupdate[i]->vectortype == VARTOOLS_VECTORTYPE_LC) {
	  if(c->lcvars_nonupdate[i]->datatype != VARTOOLS_TYPE_STRING &&
	     c->lcvars_nonupdate[i]->datatype != VARTOOLS_TYPE_CHAR) {
	    for(j=oldNJD; j < maxlcoutlength; j++) {
	      SetVariable_Value_Double(lcindex, threadindex, j, c->lcvars_nonupdate[i], 0.0);
	    }
	  } else if(c->lcvars_nonupdate[i]->datatype == VARTOOLS_TYPE_CHAR) {
	    for(j=oldNJD; j < maxlcoutlength; j++) {
	      (*((char ***) c->lcvars_nonupdate[i]->dataptr))[threadindex][j] = '0';
	    }
	  } else if(c->lcvars_nonupdate[i]->datatype == VARTOOLS_TYPE_STRING) {
	    for(j=oldNJD; j < maxlcoutlength; j++) {
	      sprintf((*((char ****) c->lcvars_nonupdate[i]->dataptr))[threadindex][j], "0");
	    }
	  }
	}
      }
    } else {
      /* If none of the output lc vectors are longer than the previous
	 light curve size; then if all lc vectors are updated by this
	 command, set the new lc size to the largest of the output lc
	 vectors, otherwise keep the light curves at their previous
	 size. Append zeros to any vectors that are shorter than the
	 light curve size. */
      if(!c->Nlcvars_nonupdate) {
	p->NJD[threadindex] = maxlcoutlength;
	/* Grow all vectors with a length less than this limit */
	for(i=0; i < c->Nvars; i++) {
	  if(c->isvaroutput[i]) {
	    if(c->vars[i]->vectortype == VARTOOLS_VECTORTYPE_LC) {
	      if(c->outlcvecs_invars[threadindex][i] < maxlcoutlength) {
		if(c->vars[i]->datatype != VARTOOLS_TYPE_STRING &&
		   c->vars[i]->datatype != VARTOOLS_TYPE_CHAR) {
		  for(j=c->outlcvecs_invars[threadindex][i]; j < maxlcoutlength; j++) {
		    SetVariable_Value_Double(lcindex, threadindex, j, c->vars[i], 0.0);
		  }
		} else if(c->vars[i]->datatype == VARTOOLS_TYPE_CHAR) {
		  for(j=c->outlcvecs_invars[threadindex][i]; j < maxlcoutlength; j++) {
		    (*((char ***) c->vars[i]->dataptr))[threadindex][j] = '0';
		  }
		} else if(c->vars[i]->datatype == VARTOOLS_TYPE_STRING) {
		  for(j=c->outlcvecs_invars[threadindex][i]; j < maxlcoutlength; j++) {
		    sprintf((*((char ****) c->vars[i]->dataptr))[threadindex][j], "0");
		  }
		}
	      }
	    }
	  }
	}
	for(i=0; i < c->Nvars_outonly; i++) {
	  if(c->outonlyvars[i]->vectortype == VARTOOLS_VECTORTYPE_LC) {
	    if(c->outlcvecs_outonlyvars[threadindex][i] < maxlcoutlength) {
	      if(c->outonlyvars[i]->datatype != VARTOOLS_TYPE_STRING &&
		 c->outonlyvars[i]->datatype != VARTOOLS_TYPE_CHAR) {
		for(j=c->outlcvecs_outonlyvars[threadindex][i]; j < maxlcoutlength; j++) {
		  SetVariable_Value_Double(lcindex, threadindex, j, c->outonlyvars[i], 0.0);
		}
	      } else if(c->outonlyvars[i]->datatype == VARTOOLS_TYPE_CHAR) {
		for(j=c->outlcvecs_outonlyvars[threadindex][i]; j < maxlcoutlength; j++) {
		  (*((char ***) c->outonlyvars[i]->dataptr))[threadindex][j] = '0';
		}
	      } else if(c->outonlyvars[i]->datatype == VARTOOLS_TYPE_STRING) {
		for(j=c->outlcvecs_outonlyvars[threadindex][i]; j < maxlcoutlength; j++) {
		  sprintf((*((char ****) c->outonlyvars[i]->dataptr))[threadindex][j], "0");
		}
	      }
	    }
	  }
	}
      } else {
	/* We are keeping the old size, grow any vectors that are too short */
	for(i=0; i < c->Nvars; i++) {
	  if(c->isvaroutput[i]) {
	    if(c->vars[i]->vectortype == VARTOOLS_VECTORTYPE_LC) {
	      if(c->outlcvecs_invars[threadindex][i] < p->NJD[threadindex]) {
		if(c->vars[i]->datatype != VARTOOLS_TYPE_STRING &&
		   c->vars[i]->datatype != VARTOOLS_TYPE_CHAR) {
		  for(j=c->outlcvecs_invars[threadindex][i]; j < p->NJD[threadindex]; j++) {
		    SetVariable_Value_Double(lcindex, threadindex, j, c->vars[i], 0.0);
		  }
		} else if(c->vars[i]->datatype == VARTOOLS_TYPE_CHAR) {
		  for(j=c->outlcvecs_invars[threadindex][i]; j < p->NJD[threadindex]; j++) {
		    (*((char ***) c->vars[i]->dataptr))[threadindex][j] = '0';
		  }
		} else if(c->vars[i]->datatype == VARTOOLS_TYPE_STRING) {
		  for(j=c->outlcvecs_invars[threadindex][i]; j < p->NJD[threadindex]; j++) {
		    sprintf((*((char ****) c->vars[i]->dataptr))[threadindex][j], "0");
		  }
		}
	      }
	    }
	  }
	}
	for(i=0; i < c->Nvars_outonly; i++) {
	  if(c->outonlyvars[i]->vectortype == VARTOOLS_VECTORTYPE_LC) {
	    if(c->outlcvecs_outonlyvars[threadindex][i] < maxlcoutlength) {
	      if(c->outonlyvars[i]->datatype != VARTOOLS_TYPE_STRING &&
		 c->outonlyvars[i]->datatype != VARTOOLS_TYPE_CHAR) {
		for(j=c->outlcvecs_outonlyvars[threadindex][i]; j < p->NJD[threadindex]; j++) {
		  SetVariable_Value_Double(lcindex, threadindex, j, c->outonlyvars[i], 0.0);
		}
	      } else if(c->outonlyvars[i]->datatype == VARTOOLS_TYPE_CHAR) {
		for(j=c->outlcvecs_outonlyvars[threadindex][i]; j < p->NJD[threadindex]; j++) {
		  (*((char ***) c->outonlyvars[i]->dataptr))[threadindex][j] = '0';
		}
	      } else if(c->outonlyvars[i]->datatype == VARTOOLS_TYPE_STRING) {
		for(j=c->outlcvecs_outonlyvars[threadindex][i]; j < p->NJD[threadindex]; j++) {
		  sprintf((*((char ****) c->outonlyvars[i]->dataptr))[threadindex][j], "0");
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  if(tmpinstr != NULL) free(tmpinstr);
  return 0;
}


int ReadVariablesFromChildRProcess_all_lcs(ProgramData *p, 
						_RCommand *c)
{
  int i, j, k, ii, tmpindex;
  int Nlcsout;
  int Nlcs = p->Nlcs;

  if(read(c->sockets[0][0], &Nlcsout, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Error receiving the number of light curves to read to the child R process.\n");
    return 1;
  }
  if(Nlcsout != p->Nlcs) {
    fprintf(stderr,"Error: the number of light curves implied by the length of the variables returned by R does not match the number of light curves sent to R.\n");
    return 1;
  }
  for(j=0; j < Nlcs; j++) {
    if(ReadVariablesFromChildRProcess(p, j, j, 0, c))
      return 1;
  }
  return 0;
}


void StopRunningRCommand(ProgramData *p, int threadindex, _RCommand *c)
{
  int msg = VARTOOLS_R_MESSAGE_ENDPROCESS;
  if(!c->iscontinueprocess) {
    if(c->IsRRunning[threadindex]) {
      write(c->sockets[threadindex][0], &msg, sizeof(int));
      close(c->sockets[threadindex][0]);
    }
  } else {
    if(!((_RCommand *)c->continueprocesscommandptr)->IsRRunning[threadindex]) {
      write(c->sockets[threadindex][0], &msg, sizeof(int));
      close(c->sockets[threadindex][0]);
    }
  }
  return;
}

void GetRCommandOutputColumnValues(ProgramData *p, int lcindex, int threadindex, _RCommand *c)
{
  int i;
  for(i=0; i < c->Noutcolumnvars; i++) {
    c->outcolumndata[threadindex][i] = EvaluateVariable_Double(lcindex, threadindex, 0, c->outcolumnvars[i]);
  }
}

void GetRCommandOutputColumnValues_all_lcs(ProgramData *p, _RCommand *c)
{
  int i, j;
  for(j=0; j < p->Nlcs; j++) {
    for(i=0; i < c->Noutcolumnvars; i++) {
      c->outcolumndata[j][i] = EvaluateVariable_Double(j, j, 0, c->outcolumnvars[i]);
    }
  }
}


void RunRCommand(ProgramData *p, int lcindex, int threadindex, int Rthreadindex, _RCommand *c)
{
  /* Check if the R process for this command is running for this
     thread; if not start it */
  int msg, retval;

  if(!c->iscontinueprocess) {

    if(!c->IsRRunning[Rthreadindex]) 
      StartRProcess(p, c, Rthreadindex);
    
  } else {
    if(!((_RCommand *)c->continueprocesscommandptr)->IsRRunning[Rthreadindex]) {
      StartRProcess(p, ((_RCommand *)c->continueprocesscommandptr), Rthreadindex);
    }
    
    if(!c->IsRRunning[Rthreadindex]) {
      c->IsRRunning[Rthreadindex] = 1;
      c->sockets[Rthreadindex][0] = ((_RCommand *)c->continueprocesscommandptr)->sockets[Rthreadindex][0];
    }
    
  }

  msg = VARTOOLS_R_MESSAGE_READDATA;
  if(write(c->sockets[Rthreadindex][0], &msg, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Lost communication with R subprocess\n");
    close(c->sockets[Rthreadindex][0]);
    exit(1);
  }

  if(write(c->sockets[Rthreadindex][0], &c->cid, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Lost communication with R subprocess\n");
    close(c->sockets[Rthreadindex][0]);
    exit(1);
  }
  
  if(SendVariablesToChildRProcess(p, lcindex, threadindex, Rthreadindex, c)) {
    /* The process failed; terminate with an error */
    fprintf(stderr,"Failed to send variables to R subprocess\n");
    close(c->sockets[Rthreadindex][0]);
    exit(1);
  }

  if(read(c->sockets[Rthreadindex][0], &retval, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Lost communication with R subprocess\n");
    close(c->sockets[Rthreadindex][0]);
    exit(1);
  }
  
  if(retval) {
    fprintf(stderr,"R function %d returned an error for light curve number %d\n", c->cnum, lcindex);
    exit(1);
  }
  
  if(ReadVariablesFromChildRProcess(p, lcindex, threadindex, Rthreadindex, c)) {
    /* The process failed; terminate with an error */
    fprintf(stderr,"Failed to receive output variables from the R subprocess\n");
    close(c->sockets[Rthreadindex][0]);
    exit(1);
  }
  GetRCommandOutputColumnValues(p, lcindex, threadindex, c);
  
}


void RunRCommand_all_lcs(ProgramData *p, _RCommand *c)
{
  /* Check if the R process for this command is running for this
     thread; if not start it */
  int msg, retval;

  if(!c->iscontinueprocess) {

    if(!c->IsRRunning[0]) 
      StartRProcess(p, c, 0);
    
  } else {
    if(!((_RCommand *)c->continueprocesscommandptr)->IsRRunning[0]) {
      StartRProcess(p, ((_RCommand *)c->continueprocesscommandptr), 0);
    }
    
    if(!c->IsRRunning[0]) {
      c->IsRRunning[0] = 1;
      c->sockets[0][0] = ((_RCommand *)c->continueprocesscommandptr)->sockets[0][0];
    }
    
  }

  msg = VARTOOLS_R_MESSAGE_READDATA;
  if(write(c->sockets[0][0], &msg, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Lost communication with R subprocess\n");
    close(c->sockets[0][0]);
    exit(1);
  }

  if(write(c->sockets[0][0], &c->cid, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Lost communication with R subprocess\n");
    close(c->sockets[0][0]);
    exit(1);
  }
  
  if(SendVariablesToChildRProcess_all_lcs(p, c)) {
    /* The process failed; terminate with an error */
    fprintf(stderr,"Failed to send variables to R subprocess\n");
    close(c->sockets[0][0]);
    exit(1);
  }

  if(read(c->sockets[0][0], &retval, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Lost communication with R subprocess\n");
    close(c->sockets[0][0]);
    exit(1);
  }
  
  if(retval) {
    fprintf(stderr,"R function %d returned an error\n", c->cnum);
    exit(1);
  }
  
  if(ReadVariablesFromChildRProcess_all_lcs(p, c)) {
    /* The process failed; terminate with an error */
    fprintf(stderr,"Failed to receive output variables from the R subprocess\n");
    close(c->sockets[0][0]);
    exit(1);
  }
  GetRCommandOutputColumnValues_all_lcs(p, c);
  
}


void InitRCommand(ProgramData *p, _RCommand *c, int Nlcs)
{
  int j;
  if(!p->readallflag) {
    if((c->Robjects = (void *)(((_RObjectContainer **) malloc(Nlcs * sizeof(_RObjectContainer *))))) == NULL ||
       (c->outcolumndata = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
      DO_ERROR_MEMALLOC;
    if(c->Noutcolumnvars > 0) {
      for(j=0; j < Nlcs; j++) {
	if((c->outcolumndata[j] = (double *) malloc(c->Noutcolumnvars * sizeof(double))) == NULL)
	  DO_ERROR_MEMALLOC;
      }
    }
#ifdef PARALLEL
    if(p->Nproc_allow > 1) {
      for(j = 0; j < p->Nproc_allow; j++)
	c->IsRRunning[j] = 0;
    } else {
#endif
      c->IsRRunning[0] = 0;
#ifdef PARALLEL
    }
#endif
  } else {
    if((c->Robjects = (void *)(((_RObjectContainer **) malloc(sizeof(_RObjectContainer *))))) == NULL ||
       (c->outcolumndata = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
      DO_ERROR_MEMALLOC;
    if(c->Noutcolumnvars > 0) {
      for(j=0; j < Nlcs; j++) {
	if((c->outcolumndata[j] = (double *) malloc(c->Noutcolumnvars * sizeof(double))) == NULL)
	  DO_ERROR_MEMALLOC;
      }
    }
    c->IsRRunning[0] = 0;
  }

}


int ParseRCommand(int *iret, int argc, char **argv, ProgramData *p, 
			_RCommand *c, Command *allcommands, int cnum)
/* Parse the command line for the "-R" command; the expected syntax is:
   -R < "fromfile" commandfile | commandstring >
           [ "init" < "file" initializationfile | initializationstring > |
             "continueprocess" prior_R_command_number ]
           [ "vars" variablelist | 
             [ "invars" inputvariablelist ] [ "outvars" outputvariablelist ] ]
           [ "outputcolumns" variablelist ] [ "process_all_lcs" ]
*/
{
  int i, j, k, ii;
  FILE *infile;
  char *line = NULL;
  size_t line_size = 0;
  char oldval;
  int s;

  c->cid = 0;
  c->cnum = cnum;

  i = *iret;
  if(i >= argc) {*iret = i; return 1;}

  if(!strcmp(argv[i],"fromfile")) {
    i++;
    if(i >= argc) {*iret = i; return 1;}
    if((infile = fopen(argv[i],"r")) == NULL) {
      fprintf(stderr,"Cannot Open %s\n", argv[i]);
      exit(1);
      *iret = i; return 1;
    }
    if (!fseek(infile, 0L, SEEK_END)) {
      c->len_Rcommandstring = ftell(infile) + 1;
      if(c->len_Rcommandstring <= 0) {
	fprintf(stderr,"Cannot Read %s\n", argv[i]);
	fclose(infile);
	exit(1);
	*iret = i; return 1;
      }
      if((c->Rcommandstring = (char *) malloc((c->len_Rcommandstring)*sizeof(char))) == NULL) {
	DO_ERROR_MEMALLOC;
	*iret = i; return 1;
      }
      if(fseek(infile, 0L, SEEK_SET)) {
	fprintf(stderr,"Error reading %s\n", argv[i]);
	fclose(infile);
	exit(1);
	*iret = i; return 1;
      }
      if(fread(c->Rcommandstring, sizeof(char),
	       (c->len_Rcommandstring-1),
	       infile) < (c->len_Rcommandstring-1)) {
	fprintf(stderr,"Error reading %s\n", argv[i]);
	fclose(infile);
	exit(1);
	*iret = i; return 1;
      }
      c->Rcommandstring[c->len_Rcommandstring-1] = '\0';
    } else {
      fprintf(stderr,"Error reading %s\n", argv[i]);
      fclose(infile);
      exit(1);
      *iret = i; return 1;
    }
  } else {
    c->len_Rcommandstring = strlen(argv[i])+1;
    if((c->Rcommandstring = (char *) malloc((c->len_Rcommandstring)*sizeof(char))) == NULL) {
      DO_ERROR_MEMALLOC;
      *iret = i; return 1;
    }
    sprintf(c->Rcommandstring,"%s",argv[i]);
  }
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"init")) {
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if(!strcmp(argv[i],"file")) {
	i++;
	if(i >= argc) {*iret = i; return 1;}
	if((infile = fopen(argv[i],"r")) == NULL) {
	  fprintf(stderr,"Error reading %s\n", argv[i]);
	  exit(1);
	  *iret = i; return 1;
	}
	if (!fseek(infile, 0L, SEEK_END)) {
	  c->len_Rinitializationtextstring = ftell(infile) + 1;
	  if(c->len_Rinitializationtextstring <= 0) {
	    fprintf(stderr,"Error reading %s\n", argv[i]);
	    fclose(infile);
	    exit(1);
	    *iret = i; return 1;
	  }
	  if((c->Rinitializationtext = (char *) malloc((c->len_Rinitializationtextstring)*sizeof(char))) == NULL) {
	    DO_ERROR_MEMALLOC;
	    *iret = i; return 1;
	  }
	  if(fseek(infile, 0L, SEEK_SET)) {
	    fprintf(stderr,"Error reading %s\n", argv[i]);
	    fclose(infile);
	    exit(1);
	    *iret = i; return 1;
	  }
	  if(fread(c->Rinitializationtext, sizeof(char),
		   (c->len_Rinitializationtextstring-1),
		   infile) < (c->len_Rinitializationtextstring-1)) {
	    fprintf(stderr,"Error reading %s\n", argv[i]);
	    fclose(infile);
	    exit(1);
	    *iret = i; return 1;
	  }
	  c->Rinitializationtext[c->len_Rinitializationtextstring-1] = '\0';
	} else {
	  fprintf(stderr,"Error reading %s\n", argv[i]);
	  fclose(infile);
	  exit(1);
	  *iret = i; return 1;
	}
      } else {
	c->len_Rinitializationtextstring = strlen(argv[i])+1;
	if((c->Rinitializationtext = (char *) malloc((c->len_Rinitializationtextstring)*sizeof(char))) == NULL) {
	  DO_ERROR_MEMALLOC;
	  *iret = i; return 1;
	}
	sprintf(c->Rinitializationtext,"%s",argv[i]);
      }
    }
    else if(!strcmp(argv[i],"continueprocess")) {
      c->iscontinueprocess = 1;
      i++;
      if(i >= argc) {*iret = i; return 1;}
      j = atoi(argv[i]);
      if(j <= 0) {
	fprintf(stderr,"Error parsing -R command. The value of prior_R_command_number must be >= 1.\n\n");
	*iret = i; return 1;
      }
      ii = 0;
      /* Find the previous process */
      for(k = 0; k < cnum; k++) {
	if(allcommands[k].cnum == CNUM_R) {
	  ii++;
	  if(ii == j) {
	    if(allcommands[k].RCommand->iscontinueprocess) {
	      c->continueprocesscommandptr = (void *) allcommands[k].RCommand->continueprocesscommandptr;
	    } else {
	      c->continueprocesscommandptr = (void *) allcommands[k].RCommand;
	    }
	    if(!(((_RCommand *) c->continueprocesscommandptr)->Nchildren)) {
	      if((((_RCommand *) c->continueprocesscommandptr)->childcommandptrs = (void *) malloc(sizeof(_RCommand *))) == NULL ||
		 (((_RCommand *) c->continueprocesscommandptr)->childcnumvals = (int *) malloc(sizeof(int))) == NULL) {
		DO_ERROR_MEMALLOC; *iret = i; return 1;
	      }
	    } else {
	      if((((_RCommand *) c->continueprocesscommandptr)->childcommandptrs = (void *) realloc((_RCommand **)(((_RCommand *) c->continueprocesscommandptr)->childcommandptrs), (((_RCommand *) c->continueprocesscommandptr)->Nchildren + 1)*sizeof(_RCommand *))) == NULL ||
		 (((_RCommand *) c->continueprocesscommandptr)->childcnumvals = (int *) realloc(((_RCommand *) c->continueprocesscommandptr)->childcnumvals, (((_RCommand *) c->continueprocesscommandptr)->Nchildren + 1)*sizeof(int))) == NULL) {
		DO_ERROR_MEMALLOC; *iret = i; return 1;
	      }
	    }
	    ((_RCommand **) (((_RCommand *) c->continueprocesscommandptr)->childcommandptrs))[((_RCommand *) c->continueprocesscommandptr)->Nchildren] = c;
	    ((_RCommand *) c->continueprocesscommandptr)->childcnumvals[((_RCommand *) c->continueprocesscommandptr)->Nchildren] = cnum;
	    ((_RCommand *) c->continueprocesscommandptr)->Nchildren += 1;
	    c->cid = ((_RCommand *) c->continueprocesscommandptr)->Nchildren;
	    break;
	  }
	}
      }
      if(ii < j) {
	fprintf(stderr,"Error parsing -R command. Attempting to continue processing from R command number %d from command number %d, but only %d previous instances of -R have been given on the command line.\n", j, cnum, ii);
	*iret = i; return 1;
      }
    } else
      i--;
  } else
    i--;

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"vars")) {
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if((c->inoutvarliststring = (char *) malloc((strlen(argv[i])+1)*sizeof(char))) == NULL) {
	DO_ERROR_MEMALLOC; *iret = i; return 1;
      }
      sprintf(c->inoutvarliststring,"%s",argv[i]);
    } else
      i--;
  } else
    i--;

  if(c->inoutvarliststring == NULL) {
    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"invars")) {
	i++;
	if(i >= argc) {*iret = i; return 1;}
	if((c->invarliststring = (char *) malloc((strlen(argv[i])+1)*sizeof(char))) == NULL) {
	  DO_ERROR_MEMALLOC; *iret = i; return 1;
	}
	sprintf(c->invarliststring,"%s",argv[i]);
      } else
	i--;
    } else
      i--;

    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"outvars")) {
	i++;
	if(i >= argc) {*iret = i; return 1;}
	if((c->outvarliststring = (char *) malloc((strlen(argv[i])+1)*sizeof(char))) == NULL) {
	  DO_ERROR_MEMALLOC; *iret = i; return 1;
	}
	sprintf(c->outvarliststring,"%s",argv[i]);
      } else
	i--;
    } else
      i--;
  }

  c->Noutcolumnvars = 0;
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"outputcolumns")) {
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if((c->outcolumnliststring = (char *) malloc((strlen(argv[i])+1)*sizeof(char))) == NULL) {
	DO_ERROR_MEMALLOC; *iret = i; return 1;
      }
      sprintf(c->outcolumnliststring,"%s",argv[i]);
      
      /* Now parse the column list */
      k = 0; j = 0; c->Noutcolumnvars = 0;
      do {
	if(c->outcolumnliststring[k] == '\0' || c->outcolumnliststring[k] == ',') {
	  oldval = c->outcolumnliststring[k];
	  c->outcolumnliststring[k] = '\0';
	  if(!c->Noutcolumnvars) {
	    if((c->outcolumnnames = (char **) malloc(sizeof(char *))) == NULL)
	      DO_ERROR_MEMALLOC;
	  } else {
	    if((c->outcolumnnames = (char **) realloc(c->outcolumnnames, (c->Noutcolumnvars + 1)*sizeof(char *))) == NULL)
	      DO_ERROR_MEMALLOC;
	  }
	  if((s = strlen(&(c->outcolumnliststring[j]))) == 0) {
	    fprintf(stderr,"Bad variable name \"\" given to the vartools analytic interpreter.\n");
	    exit(1);
	  }
	  if((c->outcolumnnames[c->Noutcolumnvars] = (char *) malloc((s+1))) == NULL)
	    DO_ERROR_MEMALLOC;
	  sprintf(c->outcolumnnames[c->Noutcolumnvars],"%s",&(c->outcolumnliststring[j]));
	  j = k+1;
	  c->outcolumnliststring[k] = oldval;
	  c->Noutcolumnvars += 1;
	}
	k++;
      } while(c->outcolumnliststring[k-1] != '\0');
    } else
      i--;
  } else
    i--;

  c->RequireReadAll = 0;
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"process_all_lcs")) {
      c->RequireReadAll = 1;
      p->readallflag = 1;
    } else
      i--;
  } else
    i--;
  
  if(c->invarliststring == NULL && c->outvarliststring == NULL &&
     c->inoutvarliststring == NULL)
    c->processallvariables = 1;
  else
    c->processallvariables = 0;

  if(c->inoutvarliststring != NULL) {
    k = 0; j = 0; c->Ninoutvarnames = 0;
    do {
      if(c->inoutvarliststring[k] == '\0' || c->inoutvarliststring[k] == ',') {
	oldval = c->inoutvarliststring[k];
	c->inoutvarliststring[k] = '\0';
	if(!c->Ninoutvarnames) {
	  if((c->inoutvarnames = (char **) malloc(sizeof(char *))) == NULL)
	    DO_ERROR_MEMALLOC;
	} else {
	  if((c->inoutvarnames = (char **) realloc(c->inoutvarnames, (c->Ninoutvarnames + 1)*sizeof(char *))) == NULL)
	    DO_ERROR_MEMALLOC;
	}
	if((s = strlen(&(c->inoutvarliststring[j]))) == 0) {
	  fprintf(stderr,"Bad variable name \"\" given to the vartools analytic expression evaluator\n");
	  exit(1);
	}
	if((c->inoutvarnames[c->Ninoutvarnames] = (char *) malloc((s+1))) == NULL)
	  DO_ERROR_MEMALLOC;
	sprintf(c->inoutvarnames[c->Ninoutvarnames],"%s",&(c->inoutvarliststring[j]));
	j = k+1;
	c->inoutvarliststring[k] = oldval;
	c->Ninoutvarnames += 1;
      }
      k++;
    } while(c->inoutvarliststring[k-1] != '\0');
  } else
    c->Ninoutvarnames = 0;

  if(c->invarliststring != NULL) {
    k = 0; j = 0; c->Ninvarnames = 0;
    do {
      if(c->invarliststring[k] == '\0' || c->invarliststring[k] == ',') {
	oldval = c->invarliststring[k];
	c->invarliststring[k] = '\0';
	if(!c->Ninvarnames) {
	  if((c->invarnames = (char **) malloc(sizeof(char *))) == NULL)
	    DO_ERROR_MEMALLOC;
	} else {
	  if((c->invarnames = (char **) realloc(c->invarnames, (c->Ninvarnames + 1)*sizeof(char *))) == NULL)
	    DO_ERROR_MEMALLOC;
	}
	if((s = strlen(&(c->invarliststring[j]))) == 0) {
	  fprintf(stderr,"Bad variable name \"\" given to the vartools analytic expression evaluator.\n");
	  exit(1);
	}
	if((c->invarnames[c->Ninvarnames] = (char *) malloc((s+1))) == NULL)
	  DO_ERROR_MEMALLOC;
	sprintf(c->invarnames[c->Ninvarnames],"%s",&(c->invarliststring[j]));
	j = k+1;
	c->invarliststring[k] = oldval;
	c->Ninvarnames += 1;
      }
      k++;
    } while(c->invarliststring[k-1] != '\0');
  } else
    c->Ninvarnames = 0;

  if(c->outvarliststring != NULL) {
    k = 0; j = 0; c->Noutvarnames = 0;
    do {
      if(c->outvarliststring[k] == '\0' || c->outvarliststring[k] == ',') {
	oldval = c->outvarliststring[k];
	c->outvarliststring[k] = '\0';
	if(!c->Noutvarnames) {
	  if((c->outvarnames = (char **) malloc(sizeof(char *))) == NULL)
	    DO_ERROR_MEMALLOC;
	} else {
	  if((c->outvarnames = (char **) realloc(c->outvarnames, (c->Noutvarnames + 1)*sizeof(char *))) == NULL)
	    DO_ERROR_MEMALLOC;
	}
	if((s = strlen(&(c->outvarliststring[j]))) == 0) {
	  fprintf(stderr,"Bad variable name \"\" given to the R analytic expression evaluator.\n");
	  exit(1);
	}
	if((c->outvarnames[c->Noutvarnames] = (char *) malloc((s+1))) == NULL)
	  DO_ERROR_MEMALLOC;
	sprintf(c->outvarnames[c->Noutvarnames],"%s",&(c->outvarliststring[j]));
	j = k+1;
	c->outvarliststring[k] = oldval;
	c->Noutvarnames += 1;
      }
      k++;
    } while(c->outvarliststring[k-1] != '\0');
  } else
    c->Noutvarnames = 0;  

  
  *iret = i; return 0;
}

void KillAllRProcessesOneThread(ProgramData *p, Command *allcommands, int threadindex) {
  int i;
  for(i=0; i < p->Ncommands; i++) {
    if(allcommands[i].cnum == CNUM_R) {
      StopRunningRCommand(p, threadindex, allcommands[i].RCommand);
    }
  }
}

void KillAllRProcesses(ProgramData *p, Command *allcommands) {
  int i;
#ifdef PARALLEL
  if(p->Nproc_allow > 1) {
    for(i = 0; i < p->Nproc_allow; i++)
      KillAllRProcessesOneThread(p, allcommands, i);
  } else {
#endif
    KillAllRProcessesOneThread(p, allcommands, 0);
#ifdef PARALLEL
  }
#endif
}

void StartAllRProcessesOneThread(ProgramData *p, Command *allcommands, int Rthreadindex) {
  int i;
  _RCommand *c;
  for(i=0; i < p->Ncommands; i++) {
    if(allcommands[i].cnum == CNUM_R) {
      c = allcommands[i].RCommand;
      if(!c->iscontinueprocess) {
	if(!c->IsRRunning[Rthreadindex]) 
	  StartRProcess(p, c, Rthreadindex);
      } else {
	if(!((_RCommand *)c->continueprocesscommandptr)->IsRRunning[Rthreadindex]) {
	  StartRProcess(p, ((_RCommand *)c->continueprocesscommandptr), Rthreadindex);
	}
	if(!c->IsRRunning[Rthreadindex]) {
	  c->IsRRunning[Rthreadindex] = 1;
	  c->sockets[Rthreadindex][0] = ((_RCommand *)c->continueprocesscommandptr)->sockets[Rthreadindex][0];
	}
      }
    }
  }
}


void StartAllRProcesses(ProgramData *p, Command *allcommands) {
  int i;
#ifdef PARALLEL
  if(p->Nproc_allow > 1) {
    for(i = 0; i < p->Nproc_allow; i++)
      StartAllRProcessesOneThread(p, allcommands, i);
  } else {
#endif
    StartAllRProcessesOneThread(p, allcommands, 0);
#ifdef PARALLEL
  }
#endif
}

