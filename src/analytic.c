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

/*void SetupLinfitExpression(ProgramData *p, _Linfit *l) {
  return;
  }*/

void CompileAllExpressions(ProgramData *p, Command *c)
/* This function generates _Expression objects for all of the 
   analytic expressions to be evaluated by the program */
{
  int i, j, k, test;
  _DataFromLightCurve *d;
  _DataFromInputList *d2;
  /* Setup any Variables associated with the output columns */
  for(i=0; i < p->Ncolumns; i++) {
    if(p->outcolumns[i].type != VARTOOLS_TYPE_STRING &&
       p->outcolumns[i].type != VARTOOLS_TYPE_CHAR)
      CreateVariable(p, p->outcolumns[i].columnname, p->outcolumns[i].type,
		     VARTOOLS_VECTORTYPE_OUTCOLUMN, p->outcolumns[i].ptr,
		     &(p->outcolumns[i]));
  }

  /* Look for -expr, -linfit, -nonlinfit, -if -elif or -else commands and set up their expressions */
  for(i=0; i < p->Ncommands; i++) {
    if(c[i].cnum == CNUM_EXPRESSION) {
      /* Check if a new variable is being defined here */
      test = 0;
      for(j=0; j < p->NDefinedVariables; j++) {
	if(!strcmp(p->DefinedVariables[j]->varname,
		   c[i].ExpressionCommand->lhsstring)) {
	  test = 1;
	  c[i].ExpressionCommand->outputvar = p->DefinedVariables[j];
	  break;
	}
      }
      if(!test) {
	c[i].ExpressionCommand->outputvar = CreateVariable(p, c[i].ExpressionCommand->lhsstring, VARTOOLS_TYPE_DOUBLE, VARTOOLS_VECTORTYPE_LC, NULL);
	RegisterDataFromLightCurve(p,
				   c[i].ExpressionCommand->outputvar->dataptr,
				   VARTOOLS_TYPE_DOUBLE,
				   0, 0, 0, 0, 0, NULL, 
				   c[i].ExpressionCommand->outputvar,
				   -1, c[i].ExpressionCommand->lhsstring);
      }
      /* Parse the expression */
      c[i].ExpressionCommand->expression = 
	ParseExpression(c[i].ExpressionCommand->rhsstring, p);
    }
    else if(c[i].cnum == CNUM_LINFIT) {
      InitLinfit(p, c[i].Linfit, i);
      SetupLinfitExpression(p, c[i].Linfit);
    }
    else if(c[i].cnum == CNUM_NONLINFIT) {
      InitNonlinfit(p, c[i].Nonlinfit, i);
      SetupNonlinfitExpression(p, c[i].Nonlinfit);
    }
    else if(c[i].cnum == CNUM_RESAMPLE) {
      SetupResampleExpression(p, c[i].Resample);
    }
    else if(c[i].cnum == CNUM_IF) {
      if(c[i].IfCommand->iftype != VARTOOLS_IFTYPE_FI) {
	c[i].IfCommand->ifs->expressions[c[i].IfCommand->ifindex] = ParseExpression(c[i].IfCommand->exprstring, p);
      }
    }
    else if(c[i].cnum == CNUM_RESTRICTTIMES) {
      if(c[i].RestrictTimes->restricttype == VARTOOLS_RESTRICTTIMES_JDRANGE) {
	if(c[i].RestrictTimes->minJDtype == PERTYPE_EXPR) {
	  c[i].RestrictTimes->minJDexpr = ParseExpression(c[i].RestrictTimes->minJDexprstring, p);
	}
	if(c[i].RestrictTimes->maxJDtype == PERTYPE_EXPR) {
	  c[i].RestrictTimes->maxJDexpr = ParseExpression(c[i].RestrictTimes->maxJDexprstring, p);
	}
      }
    }
    else if(c[i].cnum == CNUM_PHASE) {
      CheckCreateCommandOutputLCVariable(c[i].Phase->phasevarname,&(c[i].Phase->phasevar),p);
    }
    else if(c[i].cnum == CNUM_MANDELAGOLTRANSIT) {
      CheckCreateCommandOutputLCVariable(c[i].MandelAgolTransit->modelvarname,&(c[i].MandelAgolTransit->modelvar),p);
    }
#ifdef DYNAMICLIB
    else if(c[i].cnum == CNUM_USERCOMMAND) {
      for(j=0; j < c[i].UserCommand->Nexpr; j++) {
	*(c[i].UserCommand->UserDataExpressions[j]) = ParseExpression(c[i].UserCommand->expr_strings[j], p);
      }
    }
#endif
    /* Do any evaluations that are set through the "expr" option for
       built-in commands */
    for(j=0; j < c[i].N_setparam_expr; j++) {
      *(c[i].setparam_EvalExpressions[j]) = ParseExpression(c[i].setparam_EvalExprStrings[j], p);
    }
  }

  /* Now parse any expressions which are evaluated upon lc input */
  for(i=0; i < p->NDataFromLightCurve; i++) {
    d = &(p->DataFromLightCurve[i]);
    if(d->Ncolumns == 0 && d->incolumns[0] <= 0 && 
       (d->scanformat != NULL ? d->scanformat[0] != '\0' : 0)) {
      d->expression = ParseExpression(d->scanformat, p);
    }
  }

  /* Now parse any expressions which are evaluated upon reading the input list
   */
  for(i=0; i < p->NDataFromInputList; i++) {
    d2 = &(p->DataFromInputList[i]);
    if(d2->Ncolumns == 0 && d2->incolumns[0] <= 0 && 
       (d2->scanformat != NULL ? d2->scanformat[0] != '\0' : 0)) {
      d2->expression = ParseExpression(d2->scanformat, p);
    }
  }

  /* Now setup any Variable pointers which are required by the -o, -changevariable, or -stats commands */
  for(i=0; i < p->Ncommands; i++) {
    if(c[i].cnum == CNUM_OUTPUTLCS) {
      if(c[i].Outputlcs->usecolumnformat) {
	for(k=0; k < c[i].Outputlcs->Nvar; k++) {
	  for(j=0; j < p->NDefinedVariables; j++) {
	    if(!strcmp(c[i].Outputlcs->varnames[k],
		       p->DefinedVariables[j]->varname)) {
	      c[i].Outputlcs->variables[k] = p->DefinedVariables[j];
	      break;
	    }
	  }
	  if(j == p->NDefinedVariables) {
	    error2(ERR_UNDEFINEDVARIABLE,c[i].Outputlcs->varnames[k]);
	  }
	}
      }
    }
    else if(c[i].cnum == CNUM_CHANGEVARIABLE) {
      for(j=0; j < p->NDefinedVariables; j++) {
	if(!strcmp(c[i].Changevariable->newvarname,p->DefinedVariables[j]->varname)) {
	  if(p->DefinedVariables[j]->vectortype != VARTOOLS_VECTORTYPE_LC) {
	    error(ERR_INVALIDVARIABLEFORCHANGEVAR);
	  }
	  if(c[i].Changevariable->changevar == VARTOOLS_CHANGEVAR_ID) {
	    if(p->DefinedVariables[j]->datatype != VARTOOLS_TYPE_STRING) {
	      error(ERR_INVALIDVARIABLEFORCHANGEVAR);
	    }
	  } else {
	    if(p->DefinedVariables[j]->datatype != VARTOOLS_TYPE_DOUBLE) {
	      error(ERR_INVALIDVARIABLEFORCHANGEVAR);
	    }
	  }
	  c[i].Changevariable->newvar = p->DefinedVariables[j];
	  break;
	}
      }
      if(j == p->NDefinedVariables) {
	error2(ERR_UNDEFINEDVARIABLE,c[i].Changevariable->newvarname);
      }
    }
    else if(c[i].cnum == CNUM_STATS) {
      for(k = 0; k < c[i].Stats->Nvar; k++) {
	for(j=0; j < p->NDefinedVariables; j++) {
	  if(!strcmp(c[i].Stats->varnames[k],
		     p->DefinedVariables[j]->varname)) {
	    if(p->DefinedVariables[j]->vectortype != VARTOOLS_VECTORTYPE_LC) {
	      error(ERR_BADVARIABLETYPE_STATSCOMMAND);
	    }
	    c[i].Stats->vars[k] = p->DefinedVariables[j];
	    break;
	  }
	}
	if(j == p->NDefinedVariables) {
	  error2(ERR_UNDEFINEDVARIABLE,c[i].Stats->varnames[k]);
	}
      }
    }
    else if(c[i].cnum == CNUM_BINLC) {
      for(k = 0; k < c[i].Binlc->Nvar; k++) {
	for(j=0; j < p->NDefinedVariables; j++) {
	  if(!strcmp(c[i].Binlc->binvarnames[k],
		     p->DefinedVariables[j]->varname)) {
	    if(p->DefinedVariables[j]->vectortype != VARTOOLS_VECTORTYPE_LC) {
	      error(ERR_BADVARIABLETYPE_STATSCOMMAND);
	    }
	    c[i].Binlc->binvars[k] = p->DefinedVariables[j];
	    break;
	  }
	}
	if(j == p->NDefinedVariables) {
	  error2(ERR_UNDEFINEDVARIABLE,c[i].Binlc->binvarnames[k]);
	}
      }
    }
  }
}

void RunExpressionCommand(int lcindex, int threadindex, 
			  ProgramData *p, _ExpressionCommand *c)
/* Execute the -expr command, this is the function called from
   processcommand.c */
{
  int i, j;

  double dblval;

  switch(c->outputvar->vectortype) {
  case VARTOOLS_VECTORTYPE_CONSTANT:
    dblval = EvaluateExpression(lcindex, threadindex, 0, c->expression);
    switch(c->outputvar->datatype) {
    case VARTOOLS_TYPE_DOUBLE:
      *((double *) c->outputvar->dataptr) = dblval;
      break;
    case VARTOOLS_TYPE_FLOAT:
      *((float *) c->outputvar->dataptr) = dblval;
      break;
    case VARTOOLS_TYPE_INT:
      *((int *) c->outputvar->dataptr) = dblval;
      break;
    case VARTOOLS_TYPE_LONG:
      *((long *) c->outputvar->dataptr) = dblval;
      break;
    case VARTOOLS_TYPE_SHORT:
      *((short *) c->outputvar->dataptr) = dblval;
      break;
    default:
      error(ERR_BADTYPE);
    }
    break;
  case VARTOOLS_VECTORTYPE_SCALAR:
    dblval = EvaluateExpression(lcindex, threadindex, 0, c->expression);
    switch(c->outputvar->datatype) {
    case VARTOOLS_TYPE_DOUBLE:
      (*((double **) c->outputvar->dataptr))[threadindex] = dblval;
      break;
    case VARTOOLS_TYPE_FLOAT:
      (*((float **) c->outputvar->dataptr))[threadindex] = dblval;
      break;
    case VARTOOLS_TYPE_INT:
      (*((int **) c->outputvar->dataptr))[threadindex] = dblval;
      break;
    case VARTOOLS_TYPE_LONG:
      (*((long **) c->outputvar->dataptr))[threadindex] = dblval;
      break;
    case VARTOOLS_TYPE_SHORT:
      (*((short **) c->outputvar->dataptr))[threadindex] = dblval;
      break;
    default:
      error(ERR_BADTYPE);
    }
    break;
  case VARTOOLS_VECTORTYPE_INLIST:
    dblval = EvaluateExpression(lcindex, threadindex, 0, c->expression);
    switch(c->outputvar->datatype) {
    case VARTOOLS_TYPE_DOUBLE:
      (*((double **) c->outputvar->dataptr))[lcindex] = dblval;
      break;
    case VARTOOLS_TYPE_FLOAT:
      (*((float **) c->outputvar->dataptr))[lcindex] = dblval;
      break;
    case VARTOOLS_TYPE_INT:
      (*((int **) c->outputvar->dataptr))[lcindex] = dblval;
      break;
    case VARTOOLS_TYPE_LONG:
      (*((long **) c->outputvar->dataptr))[lcindex] = dblval;
      break;
    case VARTOOLS_TYPE_SHORT:
      (*((short **) c->outputvar->dataptr))[lcindex] = dblval;
      break;
    default:
      error(ERR_BADTYPE);
    }
    break;
  case VARTOOLS_VECTORTYPE_LC:
    switch(c->outputvar->datatype) {
    case VARTOOLS_TYPE_DOUBLE:
      for(j=0; j < p->NJD[threadindex]; j++) {
	dblval = EvaluateExpression(lcindex, threadindex, j, c->expression);
	(*((double ***) c->outputvar->dataptr))[threadindex][j] = dblval;
      }
      break;
    case VARTOOLS_TYPE_FLOAT:
      for(j=0; j < p->NJD[threadindex]; j++) {
	dblval = EvaluateExpression(lcindex, threadindex, j, c->expression);
	(*((float ***) c->outputvar->dataptr))[threadindex][j] = dblval;
      }
      break;
    case VARTOOLS_TYPE_INT:
      for(j=0; j < p->NJD[threadindex]; j++) {
	dblval = EvaluateExpression(lcindex, threadindex, j, c->expression);
	(*((int ***) c->outputvar->dataptr))[threadindex][j] = dblval;
      }
      break;
    case VARTOOLS_TYPE_SHORT:
      for(j=0; j < p->NJD[threadindex]; j++) {
	dblval = EvaluateExpression(lcindex, threadindex, j, c->expression);
	(*((short ***) c->outputvar->dataptr))[threadindex][j] = dblval;
      }
      break;
    case VARTOOLS_TYPE_LONG:
      for(j=0; j < p->NJD[threadindex]; j++) {
	dblval = EvaluateExpression(lcindex, threadindex, j, c->expression);
	(*((long ***) c->outputvar->dataptr))[threadindex][j] = dblval;
      }
      break;
    default:
      error(ERR_BADTYPE);
    }
    break;
  case VARTOOLS_VECTORTYPE_OUTCOLUMN:
    error(ERR_OUTCOLUMN_ON_LHS);
  default:
    error(ERR_CODEERROR);
  }    
}

_ExpressionCommand* CreateExpressionCommand(ProgramData *p, char *argvin){
  /* Create the ExpressionCommand object when the user issues a -expr
     command. argv is the argument to -expr, it should be an equality. */
  int i, j, jstart2, k;
  _ExpressionCommand *c;
  char *argv;
  argv = (char *) malloc((strlen(argvin)+1)*sizeof(char));
  sprintf(argv,"%s",argvin);
  /* Find the left and right hand sides of the equality */
  i=0, j=0;
  while(argv[i] != '=' && argv[i] != '\0') {
    if(argv[i] != ' ') {
      argv[j] = argv[i];
      j++;
    }
    i++;
  }
  if(argv[i] != '=' || j == 0) {
    error2(ERR_INVALIDARGUMENTTOEXPR,argv);
  }
  argv[j] = '\0';
  i++;
  j++;
  jstart2 = j;
  while(argv[i] != '\0') {
    if(argv[i] != ' ') {
      argv[j] = argv[i];
      j++;
    }
    i++;
  }
  argv[j] = '\0';
  
  if((c = (_ExpressionCommand *) malloc(sizeof(_ExpressionCommand))) == NULL)
    error(ERR_MEMALLOC);

  i = strlen(argv);
  j = strlen(&(argv[jstart2]));
  if(!i || !j) {
    argv[jstart2 - 1] = '=';
    error2(ERR_INVALIDARGUMENTTOEXPR,argv);
  }
  
  if((c->lhsstring = (char *) malloc((i+1))) == NULL)
    error(ERR_MEMALLOC);
  if((c->rhsstring = (char *) malloc((j+1))) == NULL)
    error(ERR_MEMALLOC);
  
  sprintf(c->lhsstring,"%s",argv);
  sprintf(c->rhsstring,"%s",&(argv[jstart2]));

  free(argv);
  return c;
}


void ParseOutputColumnFormat(_Outputlcs *o)
{
  /* This function takes a string given with the \"columnformat\" keyword to
     the -o command, and parses it into variable names and printf format
     strings */
  int i, j, k, invarname, informat;
  char copystring[MAXLEN];
  o->Nvar = 0;
  j = 0;
  i = 0;
  invarname = 1;
  informat = 0;
  do {
    if(o->columnformat[i] == ',' || o->columnformat[i] == ':' ||
       o->columnformat[i] == '\0') {
      copystring[j] = '\0';
      if((o->columnformat[i] == ':' && informat == 1) || !j) {
	error2(ERR_BADCOLUMNFORMATSTRING,o->columnformat);
      }
      if(invarname == 1) {
	if(!o->Nvar) {
	  if((o->variables = (_Variable **) malloc(sizeof(_Variable *))) == NULL ||
	     (o->printfformats = (char **) malloc(sizeof(char *))) == NULL ||
	     (o->varnames = (char **) malloc(sizeof(char *))) == NULL) {
	    error(ERR_MEMALLOC);
	  }
	} else {
	  if((o->variables = (_Variable **) realloc(o->variables, (o->Nvar + 1)*sizeof(_Variable *))) == NULL ||
	     (o->printfformats = (char **) realloc(o->printfformats, (o->Nvar + 1)*sizeof(char *))) == NULL ||
	     (o->varnames = (char **) realloc(o->varnames, (o->Nvar + 1)*sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	}
	if((o->printfformats[o->Nvar] = (char *) malloc(MAXLEN)) == NULL ||
	   (o->varnames[o->Nvar] = (char *) malloc(MAXLEN)) == NULL)
	  error(ERR_MEMALLOC);
	o->printfformats[o->Nvar][0] = '\0';
	sprintf(o->varnames[o->Nvar],"%s",copystring);
	if(o->columnformat[i] == ':') {
	  informat = 1; invarname = 0;
	}
      }
      else {
	sprintf(o->printfformats[o->Nvar],"%s",copystring);
	informat = 0; invarname = 1;
      }
      j = 0;
      if(o->columnformat[i] == ',' || o->columnformat[i] == '\0')
	o->Nvar = o->Nvar + 1;
    } else {
      copystring[j] = o->columnformat[i];
      j++;
    }
    i++;
  } while(o->columnformat[i-1] != '\0');
}

int CheckFunctionArgVariableNameNotAcceptable(char *varname, ProgramData *p) {
  int i, j, userindx;

  /* Check if the variable name is not legitimate */
  if(varname[0] != '_' && !(varname[0] >= 'A' && varname[0] <= 'Z')
     && !(varname[0] >= 'a' && varname[0] <= 'z')) {
    return(ERR_BADVARIABLENAME);
  }
  i = 1;
  while(varname[i] != '\0') {
    if(varname[i] != '_' && !(varname[i] >= '0' && varname[i] <= '9')
       && !(varname[i] >= 'A' && varname[i] <= 'Z')
       && !(varname[i] >= 'a' && varname[i] <= 'z'))
      return(ERR_BADVARIABLENAME);
    i++;
  }
  return 0;
}

int CheckVariableNameNotAcceptable(char *varname, ProgramData *p) {
  int i, j, userindx;
  /* Check if this variable name is already defined */
  if(!strcmp(varname,"NR") || !strcmp(varname,"NF"))
    return(ERR_VARIABLEALREADYINUSE);
  for(i = 0; i < p->NDefinedVariables; i++) {
    if(!strcmp(p->DefinedVariables[i]->varname,varname)) {
      return(ERR_VARIABLEALREADYINUSE);
    }
  }

  /* Check if it is the name of a function that has been defined already */
  if(!strcmp(varname,"exp") ||
     !strcmp(varname,"log") ||
     !strcmp(varname,"log10") ||
     !strcmp(varname,"sqrt") ||
     !strcmp(varname,"abs") ||
     !strcmp(varname,"max") ||
     !strcmp(varname,"min") ||
     !strcmp(varname,"hypot") ||
     !strcmp(varname,"sin") ||
     !strcmp(varname,"sindegr") ||
     !strcmp(varname,"cos") ||
     !strcmp(varname,"cosdegr") ||
     !strcmp(varname,"tan") ||
     !strcmp(varname,"tandegr") ||
     !strcmp(varname,"asin") ||
     !strcmp(varname,"asindegr") ||
     !strcmp(varname,"acos") ||
     !strcmp(varname,"acosdegr") ||
     !strcmp(varname,"atan2") ||
     !strcmp(varname,"atan2degr") ||
     !strcmp(varname,"ceil") ||
     !strcmp(varname,"floor") ||
     !strcmp(varname,"cosh") ||
     !strcmp(varname,"sinh") ||
     !strcmp(varname,"tanh") ||
     !strcmp(varname,"erf") ||
     !strcmp(varname,"erfc") ||
     !strcmp(varname,"lgamma") ||
     !strcmp(varname,"gamma") ||
     !strcmp(varname,"round") ||
     !strcmp(varname,"theta") ||
     !strcmp(varname,"acosh") ||
     !strcmp(varname,"asinh") ||
     !strcmp(varname,"atanh") ||
     !strcmp(varname,"rand") ||
     !strcmp(varname,"gauss")) {
      return(ERR_VARIABLEALREADYINUSE);
  }

#ifdef DYNAMICLIB
  if(p->NUserFunc > 0) {
    for(userindx=0; userindx < p->NUserFunc; userindx++) {
      if(!strcmp(varname,p->UserFunc[userindx].funcname)) {
	return(ERR_VARIABLEALREADYINUSE);
	break;
      }
    }
  }
  if(p->NAnalyticUserFunc > 0) {
    for(userindx=0; userindx < p->NAnalyticUserFunc; userindx++) {
      if(!strcmp(varname,p->AnalyticUserFunc[userindx].funcname)) {
	return(ERR_VARIABLEALREADYINUSE);
	break;
      }
    }
  }
#endif


  /* Check if the variable name is not legitimate */
  if(varname[0] != '_' && !(varname[0] >= 'A' && varname[0] <= 'Z')
     && !(varname[0] >= 'a' && varname[0] <= 'z')) {
    return(ERR_BADVARIABLENAME);
  }
  i = 1;
  while(varname[i] != '\0') {
    if(varname[i] != '_' && !(varname[i] >= '0' && varname[i] <= '9')
       && !(varname[i] >= 'A' && varname[i] <= 'Z')
       && !(varname[i] >= 'a' && varname[i] <= 'z'))
      return(ERR_BADVARIABLENAME);
    i++;
  }
  return 0;
}


_Variable* CreateVariable(ProgramData *p, char *varname, char datatype, char vectortype, void *vptrinput, ...)
{
  va_list varlist;
  double *dblptr = NULL;
  float *floatptr = NULL;
  int *intptr = NULL;
  long *longptr = NULL;
  short *shortptr = NULL;
  char *charptr = NULL;
  char **stringptr = NULL;

  double **dbl2ptr = NULL;
  float **float2ptr = NULL;
  int **int2ptr = NULL;
  long **long2ptr = NULL;
  short **short2ptr = NULL;
  char **char2ptr = NULL;
  char ***string2ptr = NULL;
  
  double ***dbl3ptr = NULL;
  float ***float3ptr = NULL;
  int ***int3ptr = NULL;
  long ***long3ptr = NULL;
  short ***short3ptr = NULL;
  char ***char3ptr = NULL;
  char ****string3ptr = NULL;
  
  _Variable *v;
  int j, i, len;
  void *vptr;
  int testval;

  if(vectortype == VARTOOLS_VECTORTYPE_OUTCOLUMN && vptrinput == NULL)
    error(ERR_CREATEVARIABLE_OUTCOLUMN_NEEDS_VPTRINPUT);    


  if(vectortype == VARTOOLS_VECTORTYPE_OUTCOLUMN) {
    va_start(varlist,vptrinput);
    i=0;
    while(varname[i] != '\0' ? (varname[i] >= '0' && varname[i] <= '9') || varname[i] == '_' : 0) i++;
    varname += i;
  }

  testval = CheckVariableNameNotAcceptable(varname, p);
  if(testval)
    error2(testval, varname);

  j = p->NDefinedVariables;

  if(!p->NDefinedVariables) {
    p->NDefinedVariables = 1;
    if((p->DefinedVariables = (_Variable **) malloc(sizeof(_Variable *))) == NULL)
      error(ERR_MEMALLOC);
  } else {
    p->NDefinedVariables++;
    if((p->DefinedVariables = (_Variable **) realloc(p->DefinedVariables, p->NDefinedVariables * sizeof(_Variable *))) == NULL)
      error(ERR_MEMALLOC);
  }

  if((p->DefinedVariables[j] = (_Variable *) malloc(sizeof(_Variable))) == NULL)
    error(ERR_MEMALLOC);

  v = p->DefinedVariables[j];

  if(vectortype == VARTOOLS_VECTORTYPE_OUTCOLUMN) {
    v->outc = va_arg(varlist,OutColumn *);
  }
  
  len = strlen(varname)+1;
  if((v->varname = (char *) malloc(len)) == NULL)
    error(ERR_MEMALLOC);
  sprintf(v->varname,"%s",varname);
  
  v->vectortype = vectortype;

  if(datatype == VARTOOLS_TYPE_CONVERTJD) {
    v->datatype = VARTOOLS_TYPE_DOUBLE;
  } else {
    v->datatype = datatype;
  }

  /* Allocate memory. If this is a constant we allocate space for the
     value of the variable, otherwise we simply allocate a pointer
     that will store the vector or array. The actual vector or array
     should be separately registered using RegisterDataFromLightCurve
     or RegisterDataFromInputList so that space for it will be
     automatically allocated. */
     
  if(vptrinput == NULL) {

    switch(vectortype) {
    case VARTOOLS_VECTORTYPE_CONSTANT:
      switch(datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dblptr = (double *) malloc(sizeof(double));
	vptr = (void *) dblptr;
	break;
      case VARTOOLS_TYPE_CONVERTJD:
	dblptr = (double *) malloc(sizeof(double));
	vptr = (void *) dblptr;
	break;
      case VARTOOLS_TYPE_FLOAT:
	floatptr = (float *) malloc(sizeof(float));
	vptr = (void *) floatptr;
	break;
      case VARTOOLS_TYPE_INT:
	intptr = (int *) malloc(sizeof(int));
	vptr = (void *) intptr;
	break;
      case VARTOOLS_TYPE_LONG:
	longptr = (long *) malloc(sizeof(long));
	vptr = (void *) longptr;
	break;
      case VARTOOLS_TYPE_SHORT:
	shortptr = (short *) malloc(sizeof(short));
	vptr = (void *) shortptr;
	break;
      case VARTOOLS_TYPE_CHAR:
	charptr = (char *) malloc(sizeof(char));
	vptr = (void *) charptr;
	break;
      case VARTOOLS_TYPE_STRING:
	stringptr = (char **) malloc(sizeof(char *));
	*stringptr = (char *) malloc(MAXLEN);
	vptr = (void *) stringptr;
	break;
      default:
	error(ERR_BADTYPE);
      }
      break;
    case VARTOOLS_VECTORTYPE_SCALAR:
      switch(datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dbl2ptr = (double **) malloc(sizeof(double *));
	vptr = (void *) dbl2ptr;
	break;
      case VARTOOLS_TYPE_CONVERTJD:
	dbl2ptr = (double **) malloc(sizeof(double *));
	vptr = (void *) dbl2ptr;
	break;
      case VARTOOLS_TYPE_FLOAT:
	float2ptr = (float **) malloc(sizeof(float *));
	vptr = (void *) float2ptr;
	break;
      case VARTOOLS_TYPE_INT:
	int2ptr = (int **) malloc(sizeof(int *));
	vptr = (void *) int2ptr;
	break;
      case VARTOOLS_TYPE_LONG:
	long2ptr = (long **) malloc(sizeof(long *));
	vptr = (void *) long2ptr;
	break;
      case VARTOOLS_TYPE_SHORT:
	short2ptr = (short **) malloc(sizeof(short *));
	vptr = (void *) short2ptr;
	break;
      case VARTOOLS_TYPE_CHAR:
	char2ptr = (char **) malloc(sizeof(char *));
	vptr = (void *) char2ptr;
	break;
      case VARTOOLS_TYPE_STRING:
	string2ptr = (char ***) malloc(sizeof(char **));
	vptr = (void *) string2ptr;
	break;
      default:
	error(ERR_BADTYPE);
      }
      break;
    case VARTOOLS_VECTORTYPE_INLIST:
      switch(datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dbl2ptr = (double **) malloc(sizeof(double *));
	vptr = (void *) dbl2ptr;
	break;
      case VARTOOLS_TYPE_CONVERTJD:
	dbl2ptr = (double **) malloc(sizeof(double *));
	vptr = (void *) dbl2ptr;
	break;
      case VARTOOLS_TYPE_FLOAT:
	float2ptr = (float **) malloc(sizeof(float *));
	vptr = (void *) float2ptr;
	break;
      case VARTOOLS_TYPE_INT:
	int2ptr = (int **) malloc(sizeof(int *));
	vptr = (void *) int2ptr;
	break;
      case VARTOOLS_TYPE_LONG:
	long2ptr = (long **) malloc(sizeof(long *));
	vptr = (void *) long2ptr;
	break;
      case VARTOOLS_TYPE_SHORT:
	short2ptr = (short **) malloc(sizeof(short *));
	vptr = (void *) short2ptr;
	break;
      case VARTOOLS_TYPE_CHAR:
	char2ptr = (char **) malloc(sizeof(char *));
	vptr = (void *) char2ptr;
	break;
      case VARTOOLS_TYPE_STRING:
	string2ptr = (char ***) malloc(sizeof(char **));
	vptr = (void *) string2ptr;
	break;
      default:
	error(ERR_BADTYPE);
      }
      break;
    case VARTOOLS_VECTORTYPE_LC:
      switch(datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dbl3ptr = (double ***) malloc(sizeof(double **));
	vptr = (void *) dbl3ptr;
	break;
      case VARTOOLS_TYPE_CONVERTJD:
	dbl3ptr = (double ***) malloc(sizeof(double **));
	vptr = (void *) dbl3ptr;
	break;
      case VARTOOLS_TYPE_FLOAT:
	float3ptr = (float ***) malloc(sizeof(float **));
	vptr = (void *) float3ptr;
	break;
      case VARTOOLS_TYPE_INT:
	int3ptr = (int ***) malloc(sizeof(int **));
	vptr = (void *) int3ptr;
	break;
      case VARTOOLS_TYPE_LONG:
	long3ptr = (long ***) malloc(sizeof(long **));
	vptr = (void *) long3ptr;
	break;
      case VARTOOLS_TYPE_SHORT:
	short3ptr = (short ***) malloc(sizeof(short **));
	vptr = (void *) short3ptr;
	break;
      case VARTOOLS_TYPE_CHAR:
	char3ptr = (char ***) malloc(sizeof(char **));
	vptr = (void *) char3ptr;
	break;
      case VARTOOLS_TYPE_STRING:
	string3ptr = (char ****) malloc(sizeof(char ***));
	vptr = (void *) string3ptr;
	break;
      default:
	error(ERR_BADTYPE);
      }
      break;
    default:
      error(ERR_CODEERROR);
    }
  } else {
    vptr = vptrinput;
  }
  
  v->dataptr = vptr;

  return v;
  
}

double EvaluateExpression(int lcindex, int threadindex, int jdindex, _Expression *expression)
{
  double val1, val2;
  switch(expression->op1type) {
  case VARTOOLS_OPERANDTYPE_CONSTANT:
    val1 = expression->op1_constant;
    break;
  case VARTOOLS_OPERANDTYPE_VARIABLE:
    val1 = EvaluateVariable_Double(lcindex, threadindex, jdindex, expression->op1_variable);
    break;
  case VARTOOLS_OPERANDTYPE_EXPRESSION:
    val1 = EvaluateExpression(lcindex, threadindex, jdindex, (_Expression *)(expression->op1_expression));
    break;
  case VARTOOLS_OPERANDTYPE_FUNCTION:
    val1 = EvaluateFunctionCall(lcindex, threadindex, jdindex, (_FunctionCall *)(expression->op1_functioncall));
    break;
  case VARTOOLS_OPERANDTYPE_ITERATORNR:
    val1 = (double) jdindex;
    break;
  case VARTOOLS_OPERANDTYPE_ITERATORNF:
    val1 = (double) lcindex;
    break;
  default:
    error(ERR_CODEERROR);
  }
				
  if(expression->operatortype == VARTOOLS_OPERATORTYPE_CONSTANT)
    return val1;
  else if(expression->operatortype == VARTOOLS_OPERATORTYPE_NOT)
    return ((double) (!val1));

  switch(expression->op2type) {
  case VARTOOLS_OPERANDTYPE_CONSTANT:
    val2 = expression->op2_constant;
    break;
  case VARTOOLS_OPERANDTYPE_VARIABLE:
    val2 = EvaluateVariable_Double(lcindex, threadindex, jdindex, expression->op2_variable);
    break;
  case VARTOOLS_OPERANDTYPE_EXPRESSION:
    val2 = EvaluateExpression(lcindex, threadindex, jdindex, (_Expression *)(expression->op2_expression));
    break;
  case VARTOOLS_OPERANDTYPE_FUNCTION:
    val2 = EvaluateFunctionCall(lcindex, threadindex, jdindex, (_FunctionCall *)(expression->op2_functioncall));
    break;
  case VARTOOLS_OPERANDTYPE_ITERATORNR:
    val2 = (double) jdindex;
    break;
  case VARTOOLS_OPERANDTYPE_ITERATORNF:
    val2 = (double) lcindex;
    break;
  default:
    error(ERR_CODEERROR);
  }
  
  switch(expression->operatortype) {
  case VARTOOLS_OPERATORTYPE_ADD:
    return (val1+val2);
    break;
  case VARTOOLS_OPERATORTYPE_SUBTRACT:
    return (val1-val2);
    break;
  case VARTOOLS_OPERATORTYPE_MULTIPLY:
    return (val1*val2);
    break;
  case VARTOOLS_OPERATORTYPE_DIVIDE:
    return (val1/val2);
    break;
  case VARTOOLS_OPERATORTYPE_MODULO:
    return fmod(val1, val2);
    break;
  case VARTOOLS_OPERATORTYPE_POWER:
    return pow(val1,val2);
    break;
  case VARTOOLS_OPERATORTYPE_GREATERTHAN:
    return ((double) (val1 > val2));
    break;
  case VARTOOLS_OPERATORTYPE_GREATERTHANEQUAL:
    return ((double) (val1 >= val2));
    break;
  case VARTOOLS_OPERATORTYPE_LESSTHAN:
    return ((double) (val1 < val2));
    break;
  case VARTOOLS_OPERATORTYPE_LESSTHANEQUAL:
    return ((double) (val1 <= val2));
    break;
  case VARTOOLS_OPERATORTYPE_ISEQUAL:
    return ((double) (val1 == val2));
    break;
  case VARTOOLS_OPERATORTYPE_NOTEQUAL:
    return ((double) (val1 != val2));
    break;
  case VARTOOLS_OPERATORTYPE_NOT:
    error(ERR_CODEERROR);
  case VARTOOLS_OPERATORTYPE_AND:
    return ((double) (val1 && val2));
    break;
  case VARTOOLS_OPERATORTYPE_OR:
    return ((double) (val1 || val2));
    break;
  default:
    error(ERR_CODEERROR);
  }
  
  return 0.0;
}

void SetVariable_Value_Double(int lcindex, int threadindex, int jdindex, _Variable *var, double val)
{
  switch(var->vectortype) {
  case VARTOOLS_VECTORTYPE_CONSTANT:
    switch(var->datatype) {
    case VARTOOLS_TYPE_DOUBLE:
      *((double *) var->dataptr) = val;
      break;
    case VARTOOLS_TYPE_FLOAT:
      *((float *) var->dataptr) = (float) val;
      break;
    case VARTOOLS_TYPE_INT:
      *((int *) var->dataptr) = (int) val;
      break;
    case VARTOOLS_TYPE_SHORT:
      *((short *) var->dataptr) = (short) val;
      break;
    case VARTOOLS_TYPE_LONG:
      *((long *) var->dataptr) = (long) val;
      break;
    case VARTOOLS_TYPE_STRING:
      error(ERR_BADTYPE);
      break;
    case VARTOOLS_TYPE_CHAR:
      error(ERR_BADTYPE);
      break;
    default:
      error(ERR_BADTYPE);
      break;
    }
    break;
  case VARTOOLS_VECTORTYPE_SCALAR:
    switch(var->datatype) {
    case VARTOOLS_TYPE_DOUBLE:
      (*((double **) var->dataptr))[threadindex] = val;
      break;
    case VARTOOLS_TYPE_FLOAT:
      (*((float **) var->dataptr))[threadindex] = (float) val;
      break;
    case VARTOOLS_TYPE_INT:
      (*((int **) var->dataptr))[threadindex] = (int) val;
      break;
    case VARTOOLS_TYPE_SHORT:
      (*((short **) var->dataptr))[threadindex] = (short) val;
      break;
    case VARTOOLS_TYPE_LONG:
      (*((long **) var->dataptr))[threadindex] = (long) val;
      break;
    case VARTOOLS_TYPE_STRING:
      error(ERR_BADTYPE);
      break;
    case VARTOOLS_TYPE_CHAR:
      error(ERR_BADTYPE);
      break;
    default:
      error(ERR_BADTYPE);
      break;
    }
    break;
  case VARTOOLS_VECTORTYPE_INLIST:
    switch(var->datatype) {
    case VARTOOLS_TYPE_DOUBLE:
      (*((double **) var->dataptr))[lcindex] = val;
      break;
    case VARTOOLS_TYPE_FLOAT:
      (*((float **) var->dataptr))[lcindex] = (float) val;
      break;
    case VARTOOLS_TYPE_INT:
      (*((int **) var->dataptr))[lcindex] = (int) val;
      break;
    case VARTOOLS_TYPE_SHORT:
      (*((short **) var->dataptr))[lcindex] = (short) val;
      break;
    case VARTOOLS_TYPE_LONG:
      (*((long **) var->dataptr))[lcindex] = (long) val;
      break;
    case VARTOOLS_TYPE_STRING:
      error(ERR_BADTYPE);
      break;
    case VARTOOLS_TYPE_CHAR:
      error(ERR_BADTYPE);
      break;
    default:
      error(ERR_BADTYPE);
    }
    break;
  case VARTOOLS_VECTORTYPE_LC:
    switch(var->datatype) {
    case VARTOOLS_TYPE_DOUBLE:
      (*((double ***) var->dataptr))[threadindex][jdindex] = val;
      break;
    case VARTOOLS_TYPE_FLOAT:
      (*((float ***) var->dataptr))[threadindex][jdindex] = (float) val;
      break;
    case VARTOOLS_TYPE_INT:
      (*((int ***) var->dataptr))[threadindex][jdindex] = (int) val;
      break;
    case VARTOOLS_TYPE_SHORT:
      (*((short ***) var->dataptr))[threadindex][jdindex] = (short) val;
      break;
    case VARTOOLS_TYPE_LONG:
      (*((long ***) var->dataptr))[threadindex][jdindex] = (long) val;
      break;
    case VARTOOLS_TYPE_STRING:
      error(ERR_BADTYPE);
      break;
    case VARTOOLS_TYPE_CHAR:
      error(ERR_BADTYPE);
      break;
    default:
      error(ERR_BADTYPE);
    }
    break;
  case VARTOOLS_VECTORTYPE_OUTCOLUMN:
    setoutcolumnvalue(var->outc, threadindex, lcindex, VARTOOLS_TYPE_DOUBLE, &val);
    break;
  default:
    error(ERR_CODEERROR);
  }
}

double EvaluateVariable_Double(int lcindex, int threadindex, int jdindex, _Variable *var)
{
  double val;
  switch(var->vectortype) {
  case VARTOOLS_VECTORTYPE_CONSTANT:
    switch(var->datatype) {
    case VARTOOLS_TYPE_DOUBLE:
      val = *((double *) var->dataptr);
      break;
    case VARTOOLS_TYPE_FLOAT:
      val = (double) (*((float *) var->dataptr));
      break;
    case VARTOOLS_TYPE_INT:
      val = (double) (*((int *) var->dataptr));
      break;
    case VARTOOLS_TYPE_SHORT:
      val = (double) (*((short *) var->dataptr));
      break;
    case VARTOOLS_TYPE_LONG:
      val = (double) (*((long *) var->dataptr));
      break;
    case VARTOOLS_TYPE_STRING:
      error(ERR_BADTYPE);
      break;
    case VARTOOLS_TYPE_CHAR:
      error(ERR_BADTYPE);
      break;
    default:
      error(ERR_BADTYPE);
      break;
    }
    return val;
    break;
  case VARTOOLS_VECTORTYPE_SCALAR:
    switch(var->datatype) {
    case VARTOOLS_TYPE_DOUBLE:
      val = (*((double **) var->dataptr))[threadindex];
      break;
    case VARTOOLS_TYPE_FLOAT:
      val = (double) ((*((float **) var->dataptr))[threadindex]);
      break;
    case VARTOOLS_TYPE_INT:
      val = (double) ((*((int **) var->dataptr))[threadindex]);
      break;
    case VARTOOLS_TYPE_SHORT:
      val = (double) ((*((short **) var->dataptr))[threadindex]);
      break;
    case VARTOOLS_TYPE_LONG:
      val = (double) ((*((long **) var->dataptr))[threadindex]);
      break;
    case VARTOOLS_TYPE_STRING:
      error(ERR_BADTYPE);
      break;
    case VARTOOLS_TYPE_CHAR:
      error(ERR_BADTYPE);
      break;
    default:
      error(ERR_BADTYPE);
    }
    break;
  case VARTOOLS_VECTORTYPE_INLIST:
    switch(var->datatype) {
    case VARTOOLS_TYPE_DOUBLE:
      val = (*((double **) var->dataptr))[lcindex];
      break;
    case VARTOOLS_TYPE_FLOAT:
      val = (double) ((*((float **) var->dataptr))[lcindex]);
      break;
    case VARTOOLS_TYPE_INT:
      val = (double) ((*((int **) var->dataptr))[lcindex]);
      break;
    case VARTOOLS_TYPE_SHORT:
      val = (double) ((*((short **) var->dataptr))[lcindex]);
      break;
    case VARTOOLS_TYPE_LONG:
      val = (double) ((*((long **) var->dataptr))[lcindex]);
      break;
    case VARTOOLS_TYPE_STRING:
      error(ERR_BADTYPE);
      break;
    case VARTOOLS_TYPE_CHAR:
      error(ERR_BADTYPE);
      break;
    default:
      error(ERR_BADTYPE);
    }
    break;
  case VARTOOLS_VECTORTYPE_LC:
    switch(var->datatype) {
    case VARTOOLS_TYPE_DOUBLE:
      val = (*((double ***) var->dataptr))[threadindex][jdindex];
      break;
    case VARTOOLS_TYPE_FLOAT:
      val = (double) (*((float ***) var->dataptr))[threadindex][jdindex];
      break;
    case VARTOOLS_TYPE_INT:
      val = (double) (*((int ***) var->dataptr))[threadindex][jdindex];
      break;
    case VARTOOLS_TYPE_SHORT:
      val = (double) (*((short ***) var->dataptr))[threadindex][jdindex];
      break;
    case VARTOOLS_TYPE_LONG:
      val = (double) (*((long ***) var->dataptr))[threadindex][jdindex];
      break;
    case VARTOOLS_TYPE_STRING:
      error(ERR_BADTYPE);
      break;
    case VARTOOLS_TYPE_CHAR:
      error(ERR_BADTYPE);
      break;
    default:
      error(ERR_BADTYPE);
    }
    break;
  case VARTOOLS_VECTORTYPE_OUTCOLUMN:
    getoutcolumnvalue(var->outc, threadindex, lcindex, VARTOOLS_TYPE_DOUBLE, &val);
    break;
  default:
    error(ERR_CODEERROR);
  }
  return val;
}

double EvaluateFunctionCall(int lcindex, int threadindex, int jdindex, _FunctionCall *call) {
  double *val = NULL;
  double outval;
  int i, indx;
  if(call->Nexpr > 0) {
    val = malloc(call->Nexpr * sizeof(double));
    for(i=0; i < call->Nexpr; i++) {
      val[i] = EvaluateExpression(lcindex, threadindex, jdindex, call->arguments[i]);
    }
  }

  switch(call->functionid) {
  case VARTOOLS_FUNCTIONCALL_EXP:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "exp");
    outval = exp(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_LOG:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "exp");
    outval = log(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_LOG10:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "log10");
    outval = log(val[0])/M_LN10;
    break;
  case VARTOOLS_FUNCTIONCALL_SQRT:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "sqrt");
    outval = sqrt(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_ABS:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "abs");
    outval = fabs(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_MAX:
    if(call->Nexpr != 2)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "max");
    outval = (val[0] > val[1] ? val[0] : val[1]);
    break;
  case VARTOOLS_FUNCTIONCALL_MIN:
    if(call->Nexpr != 2)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "min");
    outval = (val[0] < val[1] ? val[0] : val[1]);
    break;
  case VARTOOLS_FUNCTIONCALL_HYPOT:
    if(call->Nexpr != 2)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "hypot");
    outval = hypot(val[0],val[1]);
    break;
  case VARTOOLS_FUNCTIONCALL_SIN:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "sin");
    outval = sin(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_SINDEGR:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "sindegr");
    outval = sin(M_PI*val[0]/180.0);
    break;
  case VARTOOLS_FUNCTIONCALL_COS:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "cos");
    outval = cos(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_COSDEGR:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "cosdegr");
    outval = cos(M_PI*val[0]/180.0);
    break;
  case VARTOOLS_FUNCTIONCALL_TAN:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "tan");
    outval = tan(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_TANDEGR:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "tandegr");
    outval = tan(M_PI*val[0]/180.0);
    break;
  case VARTOOLS_FUNCTIONCALL_ASIN:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "asin");
    outval = asin(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_ASINDEGR:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "asindegr");
    outval = 180.0*asin(val[0])/M_PI;
    break;
  case VARTOOLS_FUNCTIONCALL_ACOS:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "acos");
    outval = acos(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_ACOSDEGR:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "acosdegr");
    outval = 180.0*acos(val[0])/M_PI;
    break;
  case VARTOOLS_FUNCTIONCALL_ATAN2:
    if(call->Nexpr != 2)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "atan2");
    outval = atan2(val[0],val[1]);
    break;
  case VARTOOLS_FUNCTIONCALL_ATAN2DEGR:
    if(call->Nexpr != 2)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "atan2degr");
    outval = 180.0*atan2(val[0],val[1])/M_PI;
    break;
  case VARTOOLS_FUNCTIONCALL_CEIL:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "ceil");
    outval = ceil(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_FLOOR:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "floor");
    outval = floor(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_COSH:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "cosh");
    outval = cosh(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_SINH:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "sinh");
    outval = sinh(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_TANH:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "tanh");
    outval = tanh(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_ERF:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "erf");
    outval = erf(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_ERFC:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "erfc");
    outval = erfc(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_LGAMMA:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "lgamma");
    outval = lgamma(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_GAMMA:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "gamma");
    outval = tgamma(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_ROUND:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "round");
    outval = round(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_THETA:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "theta");
    outval = (val[0] < 0.0 ? 0.0 : 1.0);
    break;
  case VARTOOLS_FUNCTIONCALL_ACOSH:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "acosh");
    outval = acosh(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_ASINH:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "asinh");
    outval = asinh(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_ATANH:
    if(call->Nexpr != 1)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "atanh");
    outval = atanh(val[0]);
    break;
  case VARTOOLS_FUNCTIONCALL_RAND:
    if(call->Nexpr != 0)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "rand");
    outval = (double) (rand() / (double) RAND_MAX);
    break;
  case VARTOOLS_FUNCTIONCALL_GAUSS:
    if(call->Nexpr != 0)
      error2(ERR_FUNCTIONCALL_INVALIDNEXPR, "gauss");
    outval = gasdev();
    break;
  default:
#ifdef DYNAMICLIB
    if(call->UserFunc != NULL) {
      if(call->Nexpr != call->UserFunc->Nargs)
	error2(ERR_FUNCTIONCALL_INVALIDNEXPR, call->UserFunc->funcname);
      outval = call->UserFunc->EvalFunction_ptr(val);
      break;
    }
    else if(call->AnalyticUserFunc != NULL) {
      if(call->Nexpr != call->AnalyticUserFunc->Nargs)
	error2(ERR_FUNCTIONCALL_INVALIDNEXPR, call->AnalyticUserFunc->funcname);
      for(i=0; i < call->Nexpr; i++) {
	SetVariable_Value_Double(lcindex, threadindex, jdindex, call->AnalyticUserFunc->input_argvars[i], val[i]);
      }
      outval = EvaluateExpression(lcindex, threadindex, jdindex, call->AnalyticUserFunc->func_expression);
      break;
    }
    else {
#endif
      error(ERR_CODEERROR);
#ifdef DYNAMICLIB
    }
#endif
  }

  if(val != NULL)
    free(val);

  return(outval);
}

/* This function takes a string and checks if it is the name of a function,
   a previously defined variable, or a constant.

   Input variable:
       term -- string to check
       p -- pointer to the ProgramData structure which stores the array of
            defined variables.
   Output variables:
       functionid -- ID number of the matched function.
       constval -- constant value.
       varptr -- pointer to the matched variable object.
   Output values:
       0 -- neither a function nor a variable.
       1 -- is a function.
       2 -- is a constant.
       3 -- is a variable.
       4 -- the string is a compound expression (operator symbols or parentheses
            are present, and it isn't a complete function call).
       5 -- the string is the NR iterator.
       6 -- the string is the NF iterator.
*/
int CheckIsFunctionConstantVariableExpression(char *term, ProgramData *p, char *functionid, double *constval, _Variable **varptr)
{
  int retval = 0;
  int i, j, userindx, test, ndec, numeterms, firstparen;
  char *term2 = NULL;
  int sizeterm2=0, Nparen;

  /* Check if term is a number */
  i=0; test = 1; ndec = 0; numeterms = 0;
  while(term[i] != '\0') {
    if(term[i] == '.') {
      ndec++;
      if(numeterms > 0) {
	/* Decimal cannot appear in an exponent, and decimal shouldn't
           appear in function or variable names */
	return 0;
      }
      if(ndec > 1) {
	/* Can't have more than one decimal in a number, and decimal
           shouldn't appear in function or variable names */
	return 0;
      }
    } 
    else if(term[i] == 'e' || term[i] == 'E') {
      /* Deal with the case of an exponential. It cannot be the first term
         in the number, it cannot appear more than once, and it must be
	 followed by a number or by a sign and a number */
      numeterms++;
      if(numeterms > 1) {test = 0; break;}
      if(i == 0) {test = 0; break;}
      if(term[i+1] == '\0') {test = 0; break;}
      if(term[i+1] >= '0' && term[i+1] <= '9') {
	i = i + 1;
      } else {
	if(term[i+1] != '+' && term[i+1] != '-') {test = 0; break;}
	if(term[i+2] == '\0') {test = 0; break;}
	if(term[i+2] < '0' || term[i+2] > '9') {return 0;}
	i = i + 2;
      }
    }
    else if(term[i] < '0' || term[i] > '9') {
      if(!((term[i] == '+' || term[i] == '-') && i == 0)) {
	test = 0;
	break;
      }
    }
    i++;
  }
  if(test) {
    *constval = atof(term);
    return 2;
  }

  /* Check if term is NR or NF */
  if(!strcmp(term,"NR")) {
    return 5;
  }
  if(!strcmp(term,"NF")) {
    return 6;
  }

  /* Check if term is a function call */
  i = 0; test = 0;
  Nparen = 0;
  firstparen = -1;
  while(term[i] != '\0' && term[i] != '(') {
    i++;
  }
  if(term[i] == '(' && i > 0) {
    firstparen = i;
    Nparen = 1;
    i++;
    while(term[i] != '\0' && !(term[i] == ')' && Nparen == 1)) {
      if(term[i] == '(') Nparen++;
      else if(term[i] == ')') Nparen--;
      i++;
    }
    if(term[i] == ')' && Nparen == 1) {
      if(term[i+1] == '\0') {
	/* This has the right form to be a function call */
	if((term2 = malloc((firstparen+1))) == NULL)
	  error(ERR_MEMALLOC);
	for(i=0; i < firstparen; i++) {
	  term2[i] = term[i];
	}
	term2[i] = '\0';
	test = 0;
	/* Check if term is the name of a function */
	if(!strcmp(term2,"exp")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_EXP;
	  test = 1;
	}
	else if(!strcmp(term2,"log")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_LOG;
	  test = 1;
	}
	else if(!strcmp(term2,"log10")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_LOG10;
	  test = 1;
	}
	else if(!strcmp(term2,"sqrt")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_SQRT;
	  test = 1;
	}
	else if(!strcmp(term2,"abs")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_ABS;
	  test = 1;
	}
	else if(!strcmp(term2,"max")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_MAX;
	  test = 1;
	}
	else if(!strcmp(term2,"min")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_MIN;
	  test = 1;
	}
	else if(!strcmp(term2,"hypot")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_HYPOT;
	  test = 1;
	}
	else if(!strcmp(term2,"sin")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_SIN;
	  test = 1;
	}
	else if(!strcmp(term2,"sindegr")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_SINDEGR;
	  test = 1;
	}
	else if(!strcmp(term2,"cos")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_COS;
	  test = 1;
	}
	else if(!strcmp(term2,"cosdegr")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_COSDEGR;
	  test = 1;
	}
	else if(!strcmp(term2,"tan")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_TAN;
	  test = 1;
	}
	else if(!strcmp(term2,"tandegr")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_TANDEGR;
	  test = 1;
	}
	else if(!strcmp(term2,"asin")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_ASIN;
	  test = 1;
	}
	else if(!strcmp(term2,"asindegr")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_ASINDEGR;
	  test = 1;
	}
	else if(!strcmp(term2,"acos")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_ACOS;
	  test = 1;
	}
	else if(!strcmp(term2,"acosdegr")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_ACOSDEGR;
	  test = 1;
	}
	else if(!strcmp(term2,"atan2")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_ATAN2;
	  test = 1;
	}
	else if(!strcmp(term2,"atan2degr")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_ATAN2DEGR;
	  test = 1;
	}
	else if(!strcmp(term2,"ceil")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_CEIL;
	  test = 1;
	}
	else if(!strcmp(term2,"floor")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_FLOOR;
	  test = 1;
	}
	else if(!strcmp(term2,"cosh")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_COSH;
	  test = 1;
	}
	else if(!strcmp(term2,"sinh")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_SINH;
	  test = 1;
	}
	else if(!strcmp(term2,"tanh")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_TANH;
	  test = 1;
	}
	else if(!strcmp(term2,"erf")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_ERF;
	  test = 1;
	}
	else if(!strcmp(term2,"erfc")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_ERFC;
	  test = 1;
	}
	else if(!strcmp(term2,"lgamma")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_LGAMMA;
	  test = 1;
	}
	else if(!strcmp(term2,"gamma")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_GAMMA;
	  test = 1;
	}
	else if(!strcmp(term2,"round")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_ROUND;
	  test = 1;
	}
	else if(!strcmp(term2,"theta")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_THETA;
	  test = 1;
	}
	else if(!strcmp(term2,"acosh")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_ACOSH;
	  test = 1;
	}
	else if(!strcmp(term2,"asinh")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_ASINH;
	  test = 1;
	}
	else if(!strcmp(term2,"atanh")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_ATANH;
	  test = 1;
	}
	else if(!strcmp(term2,"rand")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_RAND;
	  test = 1;
	}
	else if(!strcmp(term2,"gauss")) {
	  *functionid = VARTOOLS_FUNCTIONCALL_GAUSS;
	  test = 1;
	}
#ifdef DYNAMICLIB
        else {
	  if(p->NUserFunc > 0) {
	    for(userindx=0; userindx < p->NUserFunc; userindx++) {
	      if(!strcmp(term2,p->UserFunc[userindx].funcname)) {
		*functionid = VARTOOLS_FUNCTIONCALL_USERFUNC + userindx;
		test = 1;
		break;
	      }
	    }
	  }
	  if(!test && p->NAnalyticUserFunc > 0) {
	    for(userindx=0; userindx < p->NAnalyticUserFunc; userindx++) {
	      if(!strcmp(term2,p->AnalyticUserFunc[userindx].funcname)) {
		*functionid = VARTOOLS_FUNCTIONCALL_USERFUNC + p->NUserFunc + userindx;
		test = 1;
		break;
	      }
	    }
	  }
	}
#endif
	if(test == 1) {
	  free(term2);
	  return 1;
	}
	else {
	  /* Check if term2 has an operator in it */
	  test = 0;
	  for(j=0; j < strlen(term2); j++) {
	    if(term2[j] == '|' || term2[j] == '&' ||
	       term2[j] == '!' || term2[j] == '=' ||
	       term2[j] == '<' || term2[j] == '>' ||
	       term2[j] == '+' || term2[j] == '-' ||
	       term2[j] == '*' || term2[j] == '/' ||
	       term2[j] == '%' || term2[j] == '^') {
	      test = 1;
	      break;
	    }
	  }
	  if(!test) {
	    free(term2);
	    error2(ERR_ANALYTICPARSE,term);
	  }
	}
      }
    }
  }

  /* Check if term is an expression */
  i = 0; test = 0;
  while(term[i] != '\0') {
    if(term[i] == '+' ||
       term[i] == '-' || term[i] == '*' || term[i] == '/' ||
       term[i] == '%' || term[i] == '^' || term[i] == '>' ||
       (term[i] == '=' && term[i+1] == '=') || term[i] == '<' || 
       term[i] == '!' ||
       (term[i] == '&' && term[i+1] == '&') || 
       (term[i] == '|' && term[i+1] == '|'))
      return 4;
    i++;
  }

  /* Check for parentheses without a string in front */
  if(term[0] == '(' && term[i-1] == ')')
    return 4;
  

  /* Check if term is the name of a special constant */
  if(!strcmp(term,"pi")) {
    *constval = M_PI;
    return 2;
  }
  else if(!strcmp(term,"e")) {
    *constval = M_E;
    return 2;
  }
  
  /* Check if term is the name of a previously defined variable */
  for(i=0; i < p->NDefinedVariables; i++) {
    if(!strcmp(term,p->DefinedVariables[i]->varname)) {
      *varptr = p->DefinedVariables[i];
      return 3;
    }
  }

  /* term doesn't match any of the checks */
  return 0;

}

_FunctionCall* ParseFunctionCall(char *term, ProgramData *p, char functionid)
/* This function parses the string term and returns a pointer to a 
   "FunctionCall" structure which it creates */
{
  /* the format of the function call should be:
     XXXXX ( XXX, XXX, XXX, ... )
  */
  int i, ilast, j, k, Nexpr, Nparen, indx;
  _FunctionCall *retval;
  char *term2 = NULL;
  int sizeterm2 = 0;

  if((retval = (_FunctionCall *) malloc(sizeof(_FunctionCall))) == NULL) {
    error(ERR_MEMALLOC);
  }

#ifdef DYNAMICLIB
  retval->UserFunc = NULL;
  retval->AnalyticUserFunc = NULL;
#endif

  retval->functionid = functionid;

  i=0;
  Nexpr = 0;
  while(term[i] != '\0' && term[i] != '(') 
    i++;
  if(i == 0 || term[i] != '(') {
    error2(ERR_ANALYTICPARSE,term);
  }
  i++;
  
  Nparen = 0;
  ilast = i;
  while(term[i] != '\0' && !(term[i] == ')' && Nparen == 0))  {
    if(term[i] == '(') {
      Nparen++;
    }
    else if(term[i] == ')') {
      Nparen--;
    }
    else if(term[i] == ',' && Nparen == 0) {
      if(i == ilast) {
	error2(ERR_ANALYTICPARSE,term);
      }
      if(i-ilast+1 > sizeterm2) {
	if(!sizeterm2) {
	  sizeterm2 = i-ilast+1;
	  if((term2 = (char *) malloc(sizeterm2)) == NULL)
	    error(ERR_MEMALLOC);
	}
	else {
	  sizeterm2 = i-ilast+1;
	  if((term2 = (char *) realloc(term2,sizeterm2)) == NULL)
	    error(ERR_MEMALLOC);
	}
      }
      for(k=ilast,j=0;k<i;k++,j++) {
	term2[j] = term[k];
      }
      term2[j] = '\0';
      Nexpr++;
      if(Nexpr == 1) {
	if((retval->arguments = (_Expression **) malloc(sizeof(_Expression *))) == NULL)
	  error(ERR_MEMALLOC);
      } else {
	if((retval->arguments = (_Expression **) realloc(retval->arguments, Nexpr*sizeof(_Expression *))) == NULL)
	  error(ERR_MEMALLOC);
      }
      retval->arguments[Nexpr-1] = ParseExpression(term2,p);
      if(retval->arguments[Nexpr-1] == NULL)
	error2(ERR_ANALYTICPARSE,term);
      ilast = i+1;
    }
    i++;
  }
  if(!(term[i] == ')' && Nparen == 0)) {
    error2(ERR_ANALYTICPARSE,term);
  }
  if(i == ilast && Nexpr > 0) {
    error2(ERR_ANALYTICPARSE,term);
  }
  else if(i != ilast) {
    if(i-ilast+1 > sizeterm2) {
      if(!sizeterm2) {
	sizeterm2 = i-ilast+1;
	if((term2 = (char *) malloc(sizeterm2)) == NULL)
	  error(ERR_MEMALLOC);
      }
      else {
	sizeterm2 = i-ilast+1;
	if((term2 = (char *) realloc(term2,sizeterm2)) == NULL)
	  error(ERR_MEMALLOC);
      }
    }
    for(k=ilast,j=0;k<i;k++,j++) {
      term2[j] = term[k];
    }
    term2[j] = '\0';
    Nexpr++;
    if(Nexpr == 1) {
      if((retval->arguments = (_Expression **) malloc(sizeof(_Expression *))) == NULL)
	error(ERR_MEMALLOC);
    } else {
      if((retval->arguments = (_Expression **) realloc(retval->arguments, Nexpr*sizeof(_Expression *))) == NULL)
	error(ERR_MEMALLOC);
    }
    retval->arguments[Nexpr-1] = ParseExpression(term2,p);
    if(retval->arguments[Nexpr-1] == NULL)
      error2(ERR_ANALYTICPARSE,term);
    ilast = i+1;
  }

  if(term[i+1] != '\0')
    error2(ERR_ANALYTICPARSE,term);
  
  retval->Nexpr = Nexpr;

  /* Now check that the number of terms is correct for this function */
  switch(retval->functionid) {
  case VARTOOLS_FUNCTIONCALL_MAX:
    if(Nexpr != 2)
      error2(ERR_ANALYTICPARSE, term);
    break;
  case VARTOOLS_FUNCTIONCALL_MIN:
    if(Nexpr != 2)
      error2(ERR_ANALYTICPARSE, term);
    break;
  case VARTOOLS_FUNCTIONCALL_HYPOT:
    if(Nexpr != 2)
      error2(ERR_ANALYTICPARSE, term);
    break;
  case VARTOOLS_FUNCTIONCALL_ATAN2:
    if(Nexpr != 2)
      error2(ERR_ANALYTICPARSE, term);
    break;
  case VARTOOLS_FUNCTIONCALL_ATAN2DEGR:
    if(Nexpr != 2)
      error2(ERR_ANALYTICPARSE, term);
    break;
  case VARTOOLS_FUNCTIONCALL_RAND:
    if(Nexpr != 0)
      error2(ERR_ANALYTICPARSE, term);
    break;
  case VARTOOLS_FUNCTIONCALL_GAUSS:
    if(Nexpr != 0)
      error2(ERR_ANALYTICPARSE, term);
    break;
  default:
#ifdef DYNAMICLIB
    if(retval->functionid >= VARTOOLS_FUNCTIONCALL_USERFUNC && 
       retval->functionid < VARTOOLS_FUNCTIONCALL_USERFUNC + p->NUserFunc) {
      indx = retval->functionid - VARTOOLS_FUNCTIONCALL_USERFUNC;
      if(Nexpr != p->UserFunc[indx].Nargs)
	error2(ERR_ANALYTICPARSE, term);
      retval->UserFunc = &(p->UserFunc[indx]);
      break;
    }
    else if(retval->functionid >= VARTOOLS_FUNCTIONCALL_USERFUNC + p->NUserFunc && 
	    retval->functionid < VARTOOLS_FUNCTIONCALL_USERFUNC + p->NUserFunc + p->NAnalyticUserFunc) {
      indx = retval->functionid - VARTOOLS_FUNCTIONCALL_USERFUNC - p->NUserFunc;
      if(Nexpr != p->AnalyticUserFunc[indx].Nargs)
	error2(ERR_ANALYTICPARSE, term);
      retval->AnalyticUserFunc = &(p->AnalyticUserFunc[indx]);
      break;
    }
#endif
    if(Nexpr != 1)
      error2(ERR_ANALYTICPARSE, term);
    break;
  }

  if(sizeterm2 > 0)
    free(term2);

  return retval;
}

_Expression* SplitExpression(char *term, char operatortype, int i1, int sizeterm, ProgramData *p){
/* This function splits a string about an operator and returns a pointer to an
   Expression structure to evaluate the operation. 
   i1 = index of first character in term where the operator appears.
   sizeterm = total length of the string term.
   operatortype = type of operator that we are splitting on.
*/
  char *term1 = NULL, *term2 = NULL;
  _Expression *retval;
  _Variable *varptr;
  double constval;
  char functionid;
  int len = 1, i, j, expressiontype;
  if((retval = (_Expression *) malloc(sizeof(_Expression))) == NULL)
    error(ERR_MEMALLOC);
  retval->operatortype = operatortype;
  switch(operatortype) {
  case VARTOOLS_OPERATORTYPE_GREATERTHANEQUAL:
    len = 2;
    break;
  case VARTOOLS_OPERATORTYPE_LESSTHANEQUAL:
    len = 2;
    break;
  case VARTOOLS_OPERATORTYPE_ISEQUAL:
    len = 2;
    break;
  case VARTOOLS_OPERATORTYPE_NOTEQUAL:
    len = 2;
    break;
  case VARTOOLS_OPERATORTYPE_AND:
    len = 2;
    break;
  case VARTOOLS_OPERATORTYPE_OR:
    len = 2;
    break;
  default:
    len = 1;
  }
  if((i1 + len) == sizeterm) {
    error2(ERR_ANALYTICPARSE, term);
  }
  if(i1 > 0) {
    if((term1 = malloc(i1+1)) == NULL) {
      error(ERR_MEMALLOC);
    }
  }
  if((term2 = malloc(sizeterm - (i1+len) + 1)) == NULL) {
    error(ERR_MEMALLOC);
  }
  if(term1 != NULL) {
    for(i=0; i < i1; i++) {
      term1[i] = term[i];
    }
    term1[i] = '\0';
  }
  for(i=i1+len,j=0;i<sizeterm;i++,j++) {
    term2[j] = term[i];
  }
  term2[j] = '\0';

  if(term1 != NULL) {
    
    expressiontype = 
      CheckIsFunctionConstantVariableExpression(term2, p, 
						&functionid, &constval,
						&varptr);
    if(!expressiontype) {
          error2(ERR_ANALYTICPARSE, term);
    }
    else if(expressiontype == 1) {
      /***** It is a function call ****/
      retval->op2type = VARTOOLS_OPERANDTYPE_FUNCTION;
      retval->op2_functioncall = (void *) (ParseFunctionCall(term2, p, 
							     functionid));
    }
    else if(expressiontype == 2) {
      /****** It is a constant *****/
      retval->op2type = VARTOOLS_OPERANDTYPE_CONSTANT;
      retval->op2_constant = constval;
    }
    else if(expressiontype == 3) {
      /****** It is a variable *****/
      retval->op2type = VARTOOLS_OPERANDTYPE_VARIABLE;
      retval->op2_variable = varptr;
    }
    else if(expressiontype == 4) {
      /****** It is a compound expression ****/
      retval->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
      retval->op2_expression = (void *) (ParseExpression(term2, p));
      if(retval->op2_expression == NULL) {
	error2(ERR_ANALYTICPARSE, term);
      }
    }
    else if(expressiontype == 5) {
      /****** It is the NR iterator *****/
      retval->op2type = VARTOOLS_OPERANDTYPE_ITERATORNR;
    }
    else if(expressiontype == 6) {
      /****** It is the NF iterator *****/
      retval->op2type = VARTOOLS_OPERANDTYPE_ITERATORNF;
    }
  }
  else {
    term1 = term2;
    term2 = NULL;
  }
  expressiontype = 
    CheckIsFunctionConstantVariableExpression(term1, p, 
					      &functionid, &constval,
					      &varptr);
  if(!expressiontype) {
          error2(ERR_ANALYTICPARSE, term);
  }
  else if(expressiontype == 1) {
    /***** It is a function call ****/
    retval->op1type = VARTOOLS_OPERANDTYPE_FUNCTION;
    retval->op1_functioncall = (void *) (ParseFunctionCall(term1, p, 
							   functionid));
  }
  else if(expressiontype == 2) {
    /****** It is a constant *****/
    retval->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
    retval->op1_constant = constval;
  }
  else if(expressiontype == 3) {
    /****** It is a variable *****/
    retval->op1type = VARTOOLS_OPERANDTYPE_VARIABLE;
    retval->op1_variable = varptr;
  }
  else if(expressiontype == 4) {
    /****** It is a compound expression ****/
    retval->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
    retval->op1_expression = (void *) (ParseExpression(term1, p));
    if(retval->op1_expression == NULL) {
      error2(ERR_ANALYTICPARSE, term);
    }
  }
  else if(expressiontype == 5) {
    /****** It is the NR iterator *****/
    retval->op1type = VARTOOLS_OPERANDTYPE_ITERATORNR;
  }
  else if(expressiontype == 6) {
    /****** It is the NF iterator *****/
    retval->op1type = VARTOOLS_OPERANDTYPE_ITERATORNF;
  }
  if(term1 != NULL)
    free(term1);
  if(term2 != NULL)
    free(term2);
  return retval;
}
  
_Expression* ParseExpression(char *term, ProgramData *p){
/* This function parses a string called 'term' into an expression
   structure which can be evaluated with the EvaluateExpression
   function. A pointer to the created expression is returned. The
   function will exit vartools with an error if it is unable to parse
   the expression. A NULL pointer is returned if there is no expression
   to evaluate (e.g. a set of empty parentheses, or pure white space).
*/
  int starti, i, j, k, k1, k2, iendparen, istartparen, Nparen;
  int expressiontype, checknonnum;
  char functionid;
  double constval;
  _Variable *varptr;
  _Expression *retval = NULL;
  int Nsubexpressions = 0;
  char *term2 = NULL;
  int sizeterm = 0;
  int sizeterm2 = 0;

  /* First remove any white space from the expression */
  i = 0, j = 0;
  while(term[i] != '\0') {
    if(term[i] != ' ' && term[i] != '\t' && term[i] != '\n') {
      term[j] = term[i];
      j++;
    }
    i++;
  }
  term[j] = term[i];
  sizeterm = j;

  /* We parse the expression from the left-hand side. */
  
  /* We attempt to split the expression into two terms about an
     operator, starting from the lowest priority operators and proceeding
     to the highest priority operators. From lowest to highest priority, the
     operators are:
        Logical OR, Logical AND, Not equal to, Equal to, 
        greater than or equal to, greater than, less than or equal to,
        less than, addition, subtraction, modulo, division, multiplication,
        logical not, unary minus, unary plus, exponentiation, function call 
  */

  /* Check for || */
  Nparen = 0;
  i = 0;
  while(term[i] != '\0') {
    if(term[i] == '(') Nparen++;
    else if(term[i] == ')') Nparen--;
    else if(term[i] == '|' && term[i+1] == '|' && Nparen == 0) {
      if(!i) {
	error2(ERR_ANALYTICPARSE, term);
      }
      retval = SplitExpression(term, VARTOOLS_OPERATORTYPE_OR, i, sizeterm, p);
      return retval;
    }
    i++;
  }
  if(Nparen != 0)
    error2(ERR_ANALYTICPARSE, term);

  /* Check for && */
  Nparen = 0;
  i = 0;
  while(term[i] != '\0') {
    if(term[i] == '(') Nparen++;
    else if(term[i] == ')') Nparen--;
    else if(term[i] == '&' && term[i+1] == '&' && Nparen == 0) {
      if(!i) {
	error2(ERR_ANALYTICPARSE, term);
      }
      retval = SplitExpression(term, VARTOOLS_OPERATORTYPE_AND, i, sizeterm, p);
      return retval;
    }
    i++;
  }

  /* Check for != or == */
  Nparen = 0;
  i = 0;
  while(term[i] != '\0') {
    if(term[i] == '(') Nparen++;
    else if(term[i] == ')') Nparen--;
    else if(term[i] == '!' && term[i+1] == '=' && Nparen == 0) {
      if(!i) {
	error2(ERR_ANALYTICPARSE, term);
      }
      retval = SplitExpression(term, VARTOOLS_OPERATORTYPE_NOTEQUAL, i, sizeterm, p);
      return retval;
    }
    else if(term[i] == '=' && term[i+1] == '=' && Nparen == 0) {
      if(!i) {
	error2(ERR_ANALYTICPARSE, term);
      }
      retval = SplitExpression(term, VARTOOLS_OPERATORTYPE_NOTEQUAL, i, sizeterm, p);
      return retval;
    }
    i++;
  }

  /* Check for >= , >, <= or < */
  Nparen = 0;
  i = 0;
  while(term[i] != '\0') {
    if(term[i] == '(') Nparen++;
    else if(term[i] == ')') Nparen--;
    else if(term[i] == '>' && term[i+1] == '=' && Nparen == 0) {
      if(!i) {
	error2(ERR_ANALYTICPARSE, term);
      }
      retval = SplitExpression(term, VARTOOLS_OPERATORTYPE_GREATERTHANEQUAL, i, sizeterm, p);
      return retval;
    }
    else if(term[i] == '>' && term[i+1] != '=' && Nparen == 0) {
      if(!i) {
	error2(ERR_ANALYTICPARSE, term);
      }
      retval = SplitExpression(term, VARTOOLS_OPERATORTYPE_GREATERTHAN, i, sizeterm, p);
      return retval;
    }
    else if(term[i] == '<' && term[i+1] == '=' && Nparen == 0) {
      if(!i) {
	error2(ERR_ANALYTICPARSE, term);
      }
      retval = SplitExpression(term, VARTOOLS_OPERATORTYPE_LESSTHANEQUAL, i, sizeterm, p);
      return retval;
    }
    else if(term[i] == '<' && term[i+1] != '=' && Nparen == 0) {
      if(!i) {
	error2(ERR_ANALYTICPARSE, term);
      }
      retval = SplitExpression(term, VARTOOLS_OPERATORTYPE_LESSTHAN, i, sizeterm, p);
      return retval;
    }
    i++;
  }

  /* Check for + or -, note here we ignore leading terms which should be
     read as part of a number,  we also want to start from the RHS in this
     case to prevent an expression like 5-4+3 from being parsed as
     5-(4+3); finally we need to check for exponential notation. e.g.
     1.234536e-05, in which case we do not split the expression.
  */
  Nparen = 0;
  i = strlen(term);
  k2 = i;
  while(i >= 0) {
    if(term[i] == ')') Nparen++;
    else if(term[i] == '(') Nparen--;
    else if(i > 0 && term[i] == '+' && Nparen == 0) {
      checknonnum = 1;
      if(term[i-1] == 'e' || term[i-1] == 'E') {
	checknonnum = 0;
	if(i+1 < k2) {
	  if(term[i+1] < '0' || term[i+1] > '9')
	    checknonnum = 1;
	}
	else
	  checknonnum = 1;
	if(!checknonnum) {
	  k1 = i-2;
	  if(k1 < 0) checknonnum = 1;
	  else {
	    while(k1 >= 0 ? (term[k1] != '+' && term[k1] != '-' && term[k1] != '!' && term[k1] != '*' && term[k1] != '/' && term[k1] != '|' && term[k1] != '&' && term[k1] != '^' && term[k1] != '%' && term[k1] != ')' && term[k1] != '(' && term[k1] != '=' && term[k1] != '>' && term[k1] != '<') : 0) {
	      if(!(term[k1] == '.' || (term[k1] >= '0' && term[k1] <= '9'))) {
		checknonnum = 1;
		break;
	      }
	      k1--;
	    }
	  }
	}
      }
      if(checknonnum) {
	retval = SplitExpression(term, VARTOOLS_OPERATORTYPE_ADD, i, sizeterm, p);
	return retval;
      }
    }
    else if(i > 0 && term[i] == '-' && Nparen == 0) {
      checknonnum = 1;
      if(term[i-1] == 'e' || term[i-1] == 'E') {
	checknonnum = 0;
	if(i+1 < k2) {
	  if(term[i+1] < '0' || term[i+1] > '9')
	    checknonnum = 1;
	}
	else
	  checknonnum = 1;
	if(!checknonnum) {
	  k1 = i-2;
	  if(k1 < 0) checknonnum = 1;
	  else {
	    while(k1 >= 0 ? (term[k1] != '+' && term[k1] != '-' && term[k1] != '!' && term[k1] != '*' && term[k1] != '/' && term[k1] != '|' && term[k1] != '&' && term[k1] != '^' && term[k1] != '%' && term[k1] != ')' && term[k1] != '(' && term[k1] != '=' && term[k1] != '>' && term[k1] != '<') : 0) {
	      if(!(term[k1] == '.' || (term[k1] >= '0' && term[k1] <= '9'))) {
		checknonnum = 1;
		break;
	      }
	      k1--;
	    }
	  }
	}
      }
      if(checknonnum) {
	retval = SplitExpression(term, VARTOOLS_OPERATORTYPE_SUBTRACT, i, sizeterm, p);
	return retval;
      }
    }
    i--;
  }

  /* Check for *, / or % */
  Nparen = 0;
  i = strlen(term);
  while(i >= 0) {
    if(term[i] == ')') Nparen++;
    else if(term[i] == '(') Nparen--;
    else if(term[i] == '*' && Nparen == 0) {
      if(!i) {
	error2(ERR_ANALYTICPARSE, term);
      }
      retval = SplitExpression(term, VARTOOLS_OPERATORTYPE_MULTIPLY, i, sizeterm, p);
      return retval;
    }
    else if(term[i] == '/' && Nparen == 0) {
      if(!i) {
	error2(ERR_ANALYTICPARSE, term);
      }
      retval = SplitExpression(term, VARTOOLS_OPERATORTYPE_DIVIDE, i, sizeterm, p);
      return retval;
    }
    else if(term[i] == '%' && Nparen == 0) {
      if(!i) {
	error2(ERR_ANALYTICPARSE, term);
      }
      retval = SplitExpression(term, VARTOOLS_OPERATORTYPE_MODULO, i, sizeterm, p);
      return retval;
    }
    i--;
  }

  /* Check for !, only values leading an expression should be allowed. */
  Nparen = 0;
  i = 0;
  if(term[0] == '!')
    {
      retval = SplitExpression(term, VARTOOLS_OPERATORTYPE_NOT, 0, sizeterm, p);
      return retval;
    }

  /* Check for a leading '-' or a leading '+' */
  if(term[0] == '-') {
    if(sizeterm + 2 > sizeterm2) {
      if(!sizeterm2) {
	sizeterm2 = sizeterm + 2;
	if((term2 = (char *) malloc((sizeterm2))) == NULL) {
	  error(ERR_MEMALLOC);
	}
      } else {
	sizeterm2 = sizeterm + 2;
	if((term2 = (char *) realloc(term2, sizeterm2)) == NULL) {
	  error(ERR_MEMALLOC);
	}
      }
    }
    sprintf(term2,"0-%s",&(term[1]));
    retval = SplitExpression(term2, VARTOOLS_OPERATORTYPE_SUBTRACT, 1, sizeterm+1, p);
    free(term2);
    return retval;
  }
  
  if(term[0] == '+') {
    if(sizeterm + 2 > sizeterm2) {
      if(!sizeterm2) {
	sizeterm2 = sizeterm + 2;
	if((term2 = (char *) malloc((sizeterm2))) == NULL) {
	  error(ERR_MEMALLOC);
	}
      } else {
	sizeterm2 = sizeterm + 2;
	if((term2 = (char *) realloc(term2, sizeterm2)) == NULL) {
	  error(ERR_MEMALLOC);
	}
      }
    }
    sprintf(term2,"1+%s",&(term[1]));
    retval = SplitExpression(term2, VARTOOLS_OPERATORTYPE_ADD, 1, sizeterm+1, p);
    free(term2);
    return retval;
  }


  /* Check for ^ */
  Nparen = 0;
  i = 0;
  while(term[i] != '\0') {
    if(term[i] == '(') Nparen++;
    else if(term[i] == ')') Nparen--;
    else if(term[i] == '^' && Nparen == 0) {
      if(!i) {
	error2(ERR_ANALYTICPARSE, term);
      }
      retval = SplitExpression(term, VARTOOLS_OPERATORTYPE_POWER, i, sizeterm, p);
      return retval;
    }
    i++;
  }
  
  /* If the first character is a '(', this is a compound expression, remove
     the trailing ')' and pass it back */
  if(term[0] == '(') {
    if(term[sizeterm-1] != ')') {
      error2(ERR_ANALYTICPARSE, term);
    }
    if(sizeterm - 2 <= 1)
      return NULL;
    if(sizeterm - 1 > sizeterm2) {
      if(!sizeterm2) {
	sizeterm2 = sizeterm - 1;
	if((term2 = (char *) malloc((sizeterm2))) == NULL) {
	  error(ERR_MEMALLOC);
	}
      } else {
	sizeterm2 = sizeterm - 1;
	if((term2 = (char *) realloc(term2, sizeterm2)) == NULL) {
	  error(ERR_MEMALLOC);
	}
      }
    }
    for(i=1,j=0; i < sizeterm - 1; i++, j++) {
      term2[j] = term[i];
    }
    term2[j] = '\0';
    retval = ParseExpression(term2, p);
    free(term2);
    return retval;
  }
  
  /* Check if this is a variable, or a function call */
  if((retval = (_Expression *) malloc(sizeof(_Expression))) == NULL)
    error(ERR_MEMALLOC);

  expressiontype = CheckIsFunctionConstantVariableExpression(term, p, 
							     &functionid, 
							     &constval,
							     &varptr);

  if(!expressiontype)
    error2(ERR_ANALYTICPARSE, term);

  retval->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
  
  if(expressiontype == 1) {
    retval->op1type = VARTOOLS_OPERANDTYPE_FUNCTION;
    retval->op1_functioncall = (void *) (ParseFunctionCall(term, p, 
							   functionid));
  }
  else if(expressiontype == 2) {
    retval->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
    retval->op1_constant = constval;
  }
  else if(expressiontype == 3) {
    retval->op1type = VARTOOLS_OPERANDTYPE_VARIABLE;
    retval->op1_variable = varptr;
  }
  else if(expressiontype == 5) {
    retval->op1type = VARTOOLS_OPERANDTYPE_ITERATORNR;
  }
  else if(expressiontype == 6) {
    retval->op1type = VARTOOLS_OPERANDTYPE_ITERATORNF;
  }
  else {
    error2(ERR_ANALYTICPARSE, term);
  }

  return retval;
}

void PrintVartoolsFunctionList(void)
{
  OutText s;
  
  s.s = NULL;
  s.space = 0;
  s.len_s = 0;
  s.Nchar_cur_line = 0;

  printtostring(&s,"\n");
  printtostring(&s,"The following terms are understood by the VARTOOLS analytic expression evaluator:\n\n");
  printtostring(&s,"Operators:\n\n");
  printtostring(&s,"----------\n\n");
  printtostring(&s,"a+b\t\t- Addition\n\n");
  printtostring(&s,"a-b\t\t- Subtraction\n\n");
  printtostring(&s,"a*b\t\t- Multiplication\n\n");
  printtostring(&s,"a/b\t\t- Division\n\n");
  printtostring(&s,"a%b\t\t- Floating point reminder (fmod function in c)\n\n");
  printtostring(&s,"a^b\t\t- Exponentiation\n\n");
  printtostring(&s,"a>b\t\t- Greater than comparison\n\n");
  printtostring(&s,"a>=b\t\t- Greater than or equal to comparison\n\n");
  printtostring(&s,"a<b\t\t- Less than comparison\n\n");
  printtostring(&s,"a<=b\t\t- Less than or equal to comparison\n\n");
  printtostring(&s,"a==b\t\t- Logical equals\n\n");
  printtostring(&s,"a!=b\t\t- Logical not equal\n\n");
  printtostring(&s,"a&&b\t\t- Logical \"and\" comparison\n\n");
  printtostring(&s,"a||b\t\t- Logical \"or\" comparison\n\n");
  printtostring(&s,"!a\t\t- Logical \"not\"\n\n");

  printtostring(&s,"Functions:\n\n");
  printtostring(&s,"----------\n\n");
  printtostring(&s,"exp(x)\t\t- exponential of x.\n\n");
  printtostring(&s,"log(x)\t\t- natural logarithm of x.\n\n");
  printtostring(&s,"log10(x)\t\t- base 10 logarithm of x.\n\n");
  printtostring(&s,"sqrt(x)\t\t- square root of x.\n\n");
  printtostring(&s,"abs(x)\t\t- absolute value of x.\n\n");
  printtostring(&s,"max(x,y)\t\t- the larger of x or y.\n\n");
  printtostring(&s,"min(x,y)\t\t- the smaller of x or y.\n\n");
  printtostring(&s,"hypot(x,y)\t\t- sqrt(x*x+y*y).\n\n");
  printtostring(&s,"sin(x)\t\t- trigonometric sine of x. Input in radians.\n\n");
  printtostring(&s,"cos(x)\t\t- trigonometric cosine of x. Input in radians.\n\n");
  printtostring(&s,"tan(x)\t\t- trigonometric tangent of x. Input in radians.\n\n");
  printtostring(&s,"sindegr(x)\t\t- trigonometric sine of x. Input in degrees.\n\n");
  printtostring(&s,"cosdegr(x)\t\t- trigonometric cosine of x. Input in degrees.\n\n");
  printtostring(&s,"tandegr(x)\t\t- trigonometric tangent of x. Input in degrees.\n\n");
  printtostring(&s,"asin(x)\t\t- inverse sine of x. Output in radians.\n\n");
  printtostring(&s,"acos(x)\t\t- inverse cosine of x. Output in radians.\n\n");
  printtostring(&s,"atan2(y,x)\t\t- 4 quadrant inverse tangent of y/x. Output in radians.\n\n");
  printtostring(&s,"asindegr(x)\t\t- inverse sine of x. Output in degrees.\n\n");
  printtostring(&s,"acosdegr(x)\t\t- inverse cosine of x. Output in degrees.\n\n");
  printtostring(&s,"atan2degr(y,x)\t\t- 4 quadrant inverse tangent of y/x. Output in degrees.\n\n");
  printtostring(&s,"sinh(x)\t\t- hyperbolic sine of x.\n\n");
  printtostring(&s,"cosh(x)\t\t- hyperbolic cosine of x.\n\n");
  printtostring(&s,"tanh(x)\t\t- hyperbolic tangent of x.\n\n");
  printtostring(&s,"asinh(x)\t\t- inverse hyperbolic sine of x.\n\n");
  printtostring(&s,"acosh(x)\t\t- inverse hyperbolic cosine of x.\n\n");
  printtostring(&s,"atanh(x)\t\t- inverse hyperbolic tangent of x.\n\n");
  printtostring(&s,"erf(x)\t\t- error function of x.\n\n");
  printtostring(&s,"erfc(x)\t\t- complimentary error function of x.\n\n");
  printtostring(&s,"lgamma(x)\t\t- natural logarithm of the gamma function of x.\n\n");
  printtostring(&s,"gamma(x)\t\t- gamma function of x.\n\n");
  printtostring(&s,"theta(x)\t\t- 1 for x >= 0, 0 for x < 0.\n\n");
  printtostring(&s,"round(x)\t\t- round x to the nearest integer.\n\n");
  printtostring(&s,"ceil(x)\t\t- smallest integer that is greater than or equal to x.\n\n");
  printtostring(&s,"floor(x)\t\t- largest integer that is less than or equal to x.\n\n");
  printtostring(&s,"rand()\t\t- random number drawn from a uniform distribution between 0 and 1.\n\n");
  printtostring(&s,"gauss()\t\t- random number drawn from a normal distribution with 0 mean and unit variance.\n\n");
  printtostring(&s,"Constants:\n");
  printtostring(&s,"----------\n\n");
  printtostring(&s,"pi\n\n");
  printtostring(&s,"e\n\n");
  printtostring(&s,"Special Variables:\n");
  printtostring(&s,"------------------\n\n");
  printtostring(&s,"NR\t\t- image index in the light curve starting from 0\n\n");
  printtostring(&s,"NF\t\t- light curve index starting from 0\n\n");

  if(s.s != NULL) {
    printf(s.s);
    free(s.s);
  }
  exit(0);
}

char * GenerateInternalVariableName(ProgramData *p) {
  char *ret;
  int len;
  len = (p->NInternalVars / 10) + 10;
  if((ret = (char *) malloc(len)) == NULL)
    error(ERR_MEMALLOC);
  sprintf(ret,"InT_vAr_%d", p->NInternalVars);
  p->NInternalVars += 1;
  return ret;
}

void CheckCreateCommandOutputLCVariable(char *varname, _Variable **omodelvar, ProgramData *p) {
  int j;
  _Variable *v;
  if(varname != NULL) {
    for(j=0; j < p->NDefinedVariables; j++) {
      v = p->DefinedVariables[j];
      if(!strcmp(varname,v->varname)) {
	/* This is an existing variable, make sure it is the correct type */
	if(v->vectortype != VARTOOLS_VECTORTYPE_LC ||
	   v->datatype != VARTOOLS_TYPE_DOUBLE) {
	  error2(ERR_INVALIDVARIABLELCVARIABLE,varname);
	}
	*(omodelvar) = v;
	break;
      }
    }
    if(j == p->NDefinedVariables) {
      /* This is a new variable, create it */
      v = CreateVariable(p, varname, VARTOOLS_TYPE_DOUBLE,
			 VARTOOLS_VECTORTYPE_LC, NULL);
      RegisterDataFromLightCurve(p,
				v->dataptr,
				VARTOOLS_TYPE_DOUBLE,
				 0, 0, 0, 0, 0, NULL,
				 v,
				-1, varname);
      *(omodelvar) = v;
    }
  }
}

void ParseDefineAnalyticUserFunction(ProgramData *p, char *argv)
{
#ifdef DYNAMICLIB
  int i, iterm, k, j, m, m2, testval, sizeterm;
  _AnalyticUserFunc *f;
  char *funcname = NULL, **outargvars = NULL, **inargvars = NULL, *funcexprstr = NULL;
  int Narg, *inlen = NULL, *outlen = NULL, lenexpr, sizeexpr;
  if(!p->NAnalyticUserFunc) {
    if((p->AnalyticUserFunc = malloc(sizeof(_AnalyticUserFunc))) == NULL)
      error(ERR_MEMALLOC);
  } else {
    if((p->AnalyticUserFunc = realloc(p->AnalyticUserFunc,(p->NAnalyticUserFunc+1)*sizeof(_AnalyticUserFunc))) == NULL)
      error(ERR_MEMALLOC);
  }

  /* Remove any white space from the expression */
  i = 0, j = 0;
  while(argv[i] != '\0') {
    if(argv[i] != ' ' && argv[i] != '\t' && argv[i] != '\n') {
      argv[j] = argv[i];
      j++;
    }
    i++;
  }
  argv[j] = argv[i];
  sizeterm = j;


  /* First find the name of the function, and make sure that it is an acceptable name */
  i = 0;
  while(argv[i] != '\0' && argv[i] != '(') i++;
  if(argv[i] != '(' || i == 0)
    error2(ERR_INVALIDANALYTICFUNCTIONDEFINITION, argv);

  if((funcname = malloc((i+1)*sizeof(char))) == NULL)
    error(ERR_MEMALLOC);
  for(j=0; j < i; j++)
    funcname[j] = argv[j];
  funcname[j] = '\0';

  testval = CheckVariableNameNotAcceptable(funcname, p);
  if(testval)
    error2(testval, funcname);

  if(strlen(funcname) >= MAXFUNCNAMELENGTH) {
    error2(ERR_FUNCNAMETOOLONG, funcname);
  }

  /* determine the number of arguments */

  Narg = 0;
  i++;
  if(argv[i] != ')') Narg = 1;
  k = i;
  while(argv[i] != '\0' && argv[i] != ')') {
    if(argv[i] == ',') Narg++;
    i++;
  }
  if(argv[i] != ')' || i == k)
    error2(ERR_INVALIDANALYTICFUNCTIONDEFINITION, argv);
  
  if(Narg > 0) {
    if((inargvars = (char **) malloc(Narg*sizeof(char *))) == NULL ||
       (outargvars = (char **) malloc(Narg*sizeof(char *))) == NULL ||
       (inlen = (int *) malloc(Narg*sizeof(int))) == NULL ||
       (outlen = (int *) malloc(Narg*sizeof(int))) == NULL) {
      error(ERR_MEMALLOC);
    }
  }
  for(j=0; j < Narg; j++) {
    if((inargvars[j] = (char *) malloc((i-k + 2)*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
  }

  /* Get the names of the input argument variables, and make sure they are
   legitimate */

  i = k;
  j = 0;
  m = 0;
  while(argv[i] != '\0' && argv[i] != ')') {
    if(argv[i] == ',') {
      inargvars[m][j] = '\0';
      m++; j = 0;
    } else {
      inargvars[m][j] = argv[i];
      j++;
    }
    i++;
  }
  inargvars[m][j] = '\0';

  for(j=0; j < Narg; j++) {
    testval = CheckFunctionArgVariableNameNotAcceptable(inargvars[j],p);
    if(testval)
      error2(testval,inargvars[j]);
    inlen[j] = strlen(inargvars[j]);
    for(m=j+1; m < Narg; m++) {
      if(!strcmp(inargvars[j],inargvars[m])) {
	error2(ERR_ANALYTICFUNCTIONDUPLICATEINPUTARG,inargvars[j]);
      }
    }
  }

  /* Generate a set of internal variables to substitute in place of these
     external variables */
  for(j=0; j < Narg; j++) {
    outargvars[j] = GenerateInternalVariableName(p);
    outlen[j] = strlen(outargvars[j]);
  }

  /* Make sure the next term is an equal sign */
  i++;
  if(argv[i] != '=')
    error2(ERR_INVALIDANALYTICFUNCTIONDEFINITION, argv);

  sizeexpr = 2048;
  if((funcexprstr = (char *) malloc(sizeexpr*sizeof(char))) == NULL)
    error(ERR_MEMALLOC);

  i++;

  lenexpr = 0;
  /* Copy the remaining expression to funcexprstr, substituting internal variable names as needed */
  while(argv[i] != '\0') {
    /* Check if this is the start of a variable name */
    if(argv[i] == '_' || (argv[i] >= 'A' && argv[i] <= 'Z')
       || (argv[i] >= 'a' && argv[i] <= 'z')) {
      testval = 0;
      for(j=0; j < Narg; j++) {
	if(argv[i] == inargvars[j][0]) {
	  for(m = 1; m < inlen[j]; m++) {
	    if(argv[i+m] != inargvars[j][m])
	      break;
	  }
	  if(m == inlen[j]) {
	    /* If there is still a variable name underway, break */
	    if(argv[i+m] == '_' || (argv[i+m] >= '0' && argv[i+m] <= '9')
	       || (argv[i+m] >= 'A' && argv[i+m] <= 'Z')
	       || (argv[i+m] >= 'a' && argv[i+m] <= 'z'))
	      continue;
	  
	    /* Otherwise copy out the substitution */
	    for(k=0; k < outlen[j]; k++) {
	      if(lenexpr >= sizeexpr) {
		sizeexpr *= 2;
		if((funcexprstr = (char *) realloc(funcexprstr, sizeexpr*sizeof(char))) == NULL)
		  error(ERR_MEMALLOC);
	      }
	      funcexprstr[lenexpr] = outargvars[j][k];
	      lenexpr++;
	    }
	    i = i+m;
	    testval = 1;
	    break;
	  }
	}
      }
      if(!testval) {
	/* This is not one of the input arguments copy the variable expression to the string */
	if(lenexpr >= sizeexpr) {
	  sizeexpr *= 2;
	  if((funcexprstr = (char *) realloc(funcexprstr, sizeexpr*sizeof(char))) == NULL)
	    error(ERR_MEMALLOC);
	}
	funcexprstr[lenexpr] = argv[i];
	lenexpr++;
	i++;
	while(argv[i] == '_' || (argv[i] >= '0' && argv[i] <= '9') ||
	      (argv[i] >= 'A' && argv[i] <= 'Z') ||
	      (argv[i] >= 'a' && argv[i] <= 'z')) {
	  if(lenexpr >= sizeexpr) {
	    sizeexpr *= 2;
	    if((funcexprstr = (char *) realloc(funcexprstr, sizeexpr*sizeof(char))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  funcexprstr[lenexpr] = argv[i];
	  lenexpr++;
	  i++;
	}
      }
    }
    else {
      /* This is not a variable, just copy it */
      if(lenexpr >= sizeexpr) {
	sizeexpr *= 2;
	if((funcexprstr = (char *) realloc(funcexprstr, sizeexpr*sizeof(char))) == NULL)
	  error(ERR_MEMALLOC);
      }
      funcexprstr[lenexpr] = argv[i];
      lenexpr++;
      i++;
    }
  }
  if(lenexpr >= sizeexpr) {
    sizeexpr *= 2;
    if((funcexprstr = (char *) realloc(funcexprstr, sizeexpr*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
  }
  funcexprstr[lenexpr] = '\0';

  /* Copy everything over to the userfunction */
  iterm = p->NAnalyticUserFunc;
  p->NAnalyticUserFunc += 1;

  f = &(p->AnalyticUserFunc[iterm]);

  sprintf(f->funcname,"%s",funcname);
  f->Nargs = Narg;
  if(Narg > 0) {
    if((f->input_argvars = (_Variable **) malloc(Narg * sizeof(_Variable *))) == NULL)
      error(ERR_MEMALLOC);
    for(m=0; m < Narg; m++) {
      f->input_argvars[m] = CreateVariable(p, outargvars[m], 
					   VARTOOLS_TYPE_DOUBLE,
					   VARTOOLS_VECTORTYPE_CONSTANT, NULL);
    }
  }

  /* Create the expression */
  f->func_expression = ParseExpression(funcexprstr, p);

  /* Clean up any temporary storage */
  if(funcname != NULL) free(funcname);
  if(outargvars != NULL) {
    for(i=0; i < Narg; i++)
      free(outargvars[i]);
    free(outargvars);
  }
  if(inargvars != NULL) {
    for(i=0; i < Narg; i++)
      free(inargvars[i]);
    free(inargvars);
  }
  if(funcexprstr != NULL)
    free(funcexprstr);
  if(inlen != NULL)
    free(inlen);
  if(outlen != NULL)
    free(outlen);
  return;
#else
  return;
#endif
}
