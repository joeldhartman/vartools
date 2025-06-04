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

_Expression *GetLinfitVariableExpression(_Variable *v, ProgramData *p,
					 _Linfit *c, 
					 _Expression *Fullexpression);

_Expression *GetLinfitConstantExpression(ProgramData *p,
					 _Linfit *c, 
					 _Expression *Fullexpression);

int CheckFunctionCallForVariables(_FunctionCall *f, _Variable **vars, int Nvar);

int CheckExpressionForVariables(_Expression *expression, int Nvar, 
				_Variable **vars);

void DoLinfit(ProgramData *p, _Linfit *c, int threadid, int lcid) {
  /* This function executes the linear fit command. The command needs
     to first be initialized (one time) with the InitLinfit function */
  
  int Nparam, i, j, k, Njd, i1, i2;
  double **decorr;
  int *order;
  int doiterate;
  int iternum;
  int numrej;
  int numrejlast;
  double *A_errvector;
  double *constantterm;
  double *magtofit;
  double *tmpmag = NULL, *tmpmag2 = NULL;
  double rmsval;
  double modelval;
  char outname[MAXLEN];
  FILE *outfile;

  Nparam = c->Nparams;
  Njd = p->NJD[threadid];
  

  /* We will make use of the docorr function to carry out the fit.  We
   need to create the input matrices and vectors expected by this
   function. */

  if((decorr = (double **) malloc(Njd * sizeof(double *))) == NULL ||
     (magtofit = (double *) malloc(Njd * sizeof(double))) == NULL ||
     (constantterm = (double *) malloc(Njd * sizeof(double))) == NULL ||
     (order = (int *) malloc(Nparam * sizeof(int))) == NULL) {
    error(ERR_MEMALLOC);
  }
  for(j=0; j < Nparam; j++)
    order[j] = 1;
  for(i=0; i < Njd; i++) {
    if((decorr[i] = (double *) malloc(Nparam * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
    for(j=0; j < Nparam; j++) {
      decorr[i][j] = EvaluateExpression(lcid, threadid, i, c->expressions[j]);
    }
    if(c->constantexpression != NULL) {
      constantterm[i] = EvaluateExpression(lcid, threadid, i, c->constantexpression);
    } else {
      constantterm[i] = 0.;
    }
    if(c->usemask) {
      if(EvaluateVariable_Double(lcid, threadid, i, c->maskvar) > VARTOOLS_MASK_TINY)
	magtofit[i] = p->mag[threadid][i] - constantterm[i];
      else
	magtofit[i] = sqrt(-1.0);
    }
    else {
      magtofit[i] = p->mag[threadid][i] - constantterm[i];
    }
  }
  
  /* Do the fit */
  iternum = -1;
  numrej = 0;
  numrejlast = 0;
  do {
    doiterate = 0;
    docorr(magtofit, p->sig[threadid], Njd, Nparam, decorr, order, 
	   c->param_outvals[threadid], c->param_uncertainties[threadid], 0., 0,
	   0, NULL, 0, 0);

    /* Store the resulting parameter values in the appropriate variables */
    for(j=0; j < Nparam; j++) {
      (*((double **) c->params[j]->dataptr))[threadid] = c->param_outvals[threadid][j];
    }
    iternum++;
    if(c->rejectoutliers) {
      if((!c->rejiterate && iternum == 1) || 
	 (c->rejiterate && !c->rejfixnum) ||
	 (c->rejiterate && c->rejfixnum && iternum < c->rejiternum)) {
	numrejlast = numrej;
	if(tmpmag == NULL) {
	  if((tmpmag = malloc(Njd * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
	if(tmpmag2 == NULL) {
	  if((tmpmag2 = malloc(Njd * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
	k = 0;
	for(i=0; i < Njd; i++) {
	  modelval = constantterm[i];
	  for(j=0; j < Nparam; j++) {
	    modelval += c->param_outvals[threadid][j]*decorr[i][j];
	  }
	  tmpmag[i] = p->mag[threadid][i] - modelval;
	  if(!isnan(tmpmag[i]) && !isinf(tmpmag[i])) {
	    if(c->usemask) {
	      if(EvaluateVariable_Double(lcid, threadid, i, c->maskvar) <= VARTOOLS_MASK_TINY)
		continue;
	    }
	    tmpmag2[k] = tmpmag[i];
	    k++;
	  }
	}
	if(!c->rejuseMAD) {
	  rmsval = stddev(k, tmpmag2);
	} else {
	  rmsval = MAD(k, tmpmag2);
	}
	numrej=0;
	for(i=0; i < Njd; i++) {
	  if(c->usemask) {
	    if(EvaluateVariable_Double(lcid, threadid, i, c->maskvar) <= VARTOOLS_MASK_TINY)
	      continue;
	  }
	  if(!isnan(p->mag[threadid][i]) && !isinf(p->mag[threadid][i])) {
	    if(tmpmag[i] > rmsval*c->rejsigclip || tmpmag[i] < -rmsval*c->rejsigclip) {
	      magtofit[i] = sqrt(-1.0);
	      numrej++;
	    }
	    else {
	      magtofit[i] = p->mag[threadid][i] - constantterm[i];
	    }
	  }
	}
	if(!c->rejiterate && iternum == 1) {
	  if(numrej > 0) doiterate = 1;
	}
	else if(c->rejiterate && !c->rejfixnum) {
	  if(numrej > numrejlast) doiterate = 1;
	}
	else if(c->rejiterate && c->rejfixnum && iternum < c->rejiternum) {
	  if(numrej > numrejlast) doiterate = 1;
	}
      }
    }
  } while(doiterate);

  c->numrej[threadid] = numrej;
  c->iternum[threadid] = iternum;

  if(tmpmag != NULL) free(tmpmag);
  if(tmpmag2 != NULL) free(tmpmag2);
  
  /* generate the model, correct the light curve, and output the result,
     as requested */
  if(c->correctlc || c->omodel || c->modelvarname != NULL || c->calcchi2out) {
    if(c->omodel) {
      GetOutputFilename(outname,p->lcnames[lcid],c->outdir,c->outfile_extension,
			c->outfilename_format, lcid);
      /*      i1 = 0; i2 = 0;
	      while(p->lcnames[lcid][i1] != '\0')
	      {
	      if(p->lcnames[lcid][i1] == '/')
	      i2 = i1 + 1;
	      i1++;
	      }
	      sprintf(outname,"%s/%s.linfit.model", c->outdir, &p->lcnames[lcid][i2]);*/
      if((outfile = fopen(outname,"w")) == NULL)
	error2(ERR_CANNOTWRITE,outname);
    }
    if(c->calcchi2out) {
      c->chi2out[threadid] = 0.;
    }
    if(c->omodel) {
      fprintf(outfile,"#Time Mag_lc Err Mag_model\n");
    }
    for(i=0; i < Njd; i++) {
      modelval = constantterm[i];
      for(j=0; j < Nparam; j++) {
	modelval += c->param_outvals[threadid][j]*decorr[i][j];
      }
      if(c->calcchi2out) {
	if(!isnan(p->mag[threadid][i]) && p->sig[threadid][i] > 0.) {
	  c->chi2out[threadid] += ((modelval - p->mag[threadid][i])*(modelval - p->mag[threadid][i]))/p->sig[threadid][i]/p->sig[threadid][i];
	}
      }
      if(c->omodel) {
	fprintf(outfile,"%17.9f %f %f %f\n",p->t[threadid][i],
		p->mag[threadid][i],p->sig[threadid][i],modelval);
      }
      if(c->correctlc) {
	p->mag[threadid][i] -= modelval;
      }
      if(c->modelvarname != NULL) {
	(*((double ***) c->modelvar->dataptr))[threadid][i] = modelval;
      }
    }
    if(c->omodel)
      fclose(outfile);
  }
	
  free(order);
  for(i=0; i < Njd; i++)
    free(decorr[i]);
  free(decorr);
  free(magtofit);
  free(constantterm);

}

void InitLinfit(ProgramData *p, _Linfit *c, int cnum) {
  /* Initialize the Linfit command */

  int nparam, s;
  char **paramnames;
  char oldval;
  int i, j;
  _Variable *v;

  /* Create any new variables, and 
     make sure that old variables are of the correct type */
  nparam = c->Nparams;

  for(i=0; i < nparam; i++) {
    for(j=0; j < p->NDefinedVariables; j++) {
      v = p->DefinedVariables[j];
      if(!strcmp(c->paramnames[i],v->varname)) {
	/* This is an existing variable, make sure it is the correct type */
	if((v->vectortype != VARTOOLS_VECTORTYPE_INLIST 
	    && v->vectortype != VARTOOLS_VECTORTYPE_SCALAR) ||
	   v->datatype != VARTOOLS_TYPE_DOUBLE) {
	  error2(ERR_INVALIDVARIABLEFORLINFIT,v->varname);
	}
	c->params[i] = v;
	break;
      }
    }
    if(j == p->NDefinedVariables) {
      /* This is a new variable, create it */
      c->params[i] = CreateVariable(p, c->paramnames[i], VARTOOLS_TYPE_DOUBLE,
				    VARTOOLS_VECTORTYPE_INLIST, NULL);
      RegisterDataFromInputList(p,
				c->params[i]->dataptr,
				VARTOOLS_TYPE_DOUBLE,
				0, cnum, 0, 0, NULL,
				-1, c->paramnames[i]);
    }
  }

  /*
  if(nparam > 0) {
    RegisterDataFromInputList(p,
			      &(c->param_uncertainties),
			      VARTOOLS_TYPE_DOUBLE,
			      nparam, cnum, 0, 0, -nparam, "Linfit_Param_Uncertainty");
  }
  */

  if((c->expressions = (_Expression **) malloc(nparam * sizeof(_Expression *))) == NULL)
    error(ERR_MEMALLOC);

  for(i=0; i < nparam; i++)
    RegisterDataFromInputList(p,
			      c->params[i]->dataptr,
			      VARTOOLS_TYPE_DOUBLE,
			      0, cnum, 0, 0, NULL,
			      -1, "err");

  /* Set-up the model variable if it is being using */
  CheckCreateCommandOutputLCVariable(c->modelvarname,&(c->modelvar),p);

  if(c->calcchi2out) {
    RegisterScalarData(p, (void *) (&(c->chi2out)), VARTOOLS_TYPE_DOUBLE, 0);
  }


  CheckCreateCommandOutputLCVariable(c->maskvarname,&(c->maskvar),p);

}

_Expression *GetLinfitVariableExpression(_Variable *v, ProgramData *p,
					 _Linfit *c, 
					 _Expression *Fullexpression) {
  int foundit = 0;
  int i;
  _Expression *Newexpression;
  _Expression *Newexpression2;
  _Expression *Newexpression3;
  switch(Fullexpression->operatortype) {
  case VARTOOLS_OPERATORTYPE_CONSTANT:
    if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_VARIABLE ?
       Fullexpression->op1_variable == v : 0) {
      if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	error(ERR_MEMALLOC);
      Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
      Newexpression->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
      Newexpression->op1_constant = 1.0;
      return Newexpression;
    } else {
      return NULL;
    }
    break;
  case VARTOOLS_OPERATORTYPE_ADD:
    switch(Fullexpression->op1type) {
    case VARTOOLS_OPERANDTYPE_VARIABLE:
      if(Fullexpression->op1_variable == v) {
	if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	  error(ERR_MEMALLOC);
	Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	Newexpression->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	Newexpression->op1_constant = 1.0;
	foundit = 1;
      }
      break;
    case VARTOOLS_OPERANDTYPE_EXPRESSION:
      Newexpression = 
	GetLinfitVariableExpression(v, p, c, 
				    (_Expression *) Fullexpression->op1_expression);
      if(Newexpression != NULL) {
	foundit = 1;
      }
      break;
    case VARTOOLS_OPERANDTYPE_FUNCTION:
      if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op1_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
      }
      break;
    default:
      break;
    }
    switch(Fullexpression->op2type) {
    case VARTOOLS_OPERANDTYPE_VARIABLE:
      if(Fullexpression->op2_variable == v) {
	if(foundit)
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	  error(ERR_MEMALLOC);
	Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	Newexpression->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	Newexpression->op1_constant = 1.0;
	return Newexpression;
      } else {
	if(foundit)
	  return Newexpression;
	else
	  return NULL;
      }
      break;
    case VARTOOLS_OPERANDTYPE_EXPRESSION:
      if(foundit) {
	if(GetLinfitVariableExpression(v, p, c, 
				       (_Expression *) Fullexpression->op2_expression) != NULL) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	} else {
	  return Newexpression;
	}
      }
      return(GetLinfitVariableExpression(v, p, c, 
					 (_Expression *) Fullexpression->op2_expression));
      break;
    case VARTOOLS_OPERANDTYPE_FUNCTION:
      if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
      }
      if(foundit)
	return Newexpression;
      else
	return NULL;
      break;
    default:
      if(foundit)
	return Newexpression;
      else
	return NULL;
    }
    break;
  case VARTOOLS_OPERATORTYPE_SUBTRACT:
    switch(Fullexpression->op1type) {
    case VARTOOLS_OPERANDTYPE_VARIABLE:
      if(Fullexpression->op1_variable == v) {
	if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	  error(ERR_MEMALLOC);
	Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	Newexpression->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	Newexpression->op1_constant = 1.0;
	foundit = 1;
      }
      break;
    case VARTOOLS_OPERANDTYPE_EXPRESSION:
      Newexpression = 
	GetLinfitVariableExpression(v, p, c, 
				    (_Expression *) Fullexpression->op1_expression);
      if(Newexpression != NULL) {
	foundit = 1;
      }
      break;
    case VARTOOLS_OPERANDTYPE_FUNCTION:
      if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op1_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
      }
      break;
    default:
      break;
    }
    switch(Fullexpression->op2type) {
    case VARTOOLS_OPERANDTYPE_VARIABLE:
      if(Fullexpression->op2_variable == v) {
	if(foundit)
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	  error(ERR_MEMALLOC);
	Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	Newexpression->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	Newexpression->op1_constant = -1.0;
	return Newexpression;
      } else {
	if(foundit)
	  return Newexpression;
	else
	  return NULL;
      }
      break;
    case VARTOOLS_OPERANDTYPE_EXPRESSION:
      if(foundit) {
	if(GetLinfitVariableExpression(v, p, c, 
				       (_Expression *) Fullexpression->op2_expression)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	} else {
	  return Newexpression;
	}
      }
      Newexpression = GetLinfitVariableExpression(v, p, c, 
						  (_Expression *) Fullexpression->op2_expression);
      if(Newexpression == NULL) return NULL;
      if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	error(ERR_MEMALLOC);
      Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
      Newexpression2->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
      Newexpression2->op1_constant = -1.0;
      Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
      Newexpression2->op2_expression = (void *) Newexpression;
      return Newexpression2;
      break;
    case VARTOOLS_OPERANDTYPE_FUNCTION:
      if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
      }
      if(foundit)
	return Newexpression;
      else
	return NULL;
      break;
    default:
      if(foundit)
	return Newexpression;
      else
	return NULL;
    }
    break;
  case VARTOOLS_OPERATORTYPE_MULTIPLY:
    switch(Fullexpression->op1type) {
    case VARTOOLS_OPERANDTYPE_VARIABLE:
      if(Fullexpression->op1_variable == v) {
	if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_EXPRESSION) {
	  Newexpression = (_Expression *) Fullexpression->op2_expression;
	  if(CheckExpressionForVariables(Newexpression, c->Nparams, c->params)){
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	  }
	  return Newexpression;
	} else {
	  if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  Newexpression->op1type = Fullexpression->op2type;
	  switch(Newexpression->op1type) {
	  case VARTOOLS_OPERANDTYPE_CONSTANT:
	    Newexpression->op1_constant = Fullexpression->op2_constant;
	    break;
	  case VARTOOLS_OPERANDTYPE_VARIABLE:
	    for(i=0; i < c->Nparams; i++) {
	      if(Fullexpression->op2_variable == c->params[i])
		error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	    }
	    Newexpression->op1_variable = Fullexpression->op2_variable;
	    break;
	  case VARTOOLS_OPERANDTYPE_FUNCTION:
	    if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	      error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	    }
	    Newexpression->op1_functioncall = (void *) Fullexpression->op2_functioncall;
	    break;
	  default:
	    break;
	  }
	  return Newexpression;
	  
	}
	foundit = 1;
      } else {
	switch(Fullexpression->op2type) {
	case VARTOOLS_OPERANDTYPE_VARIABLE:
	  if(Fullexpression->op2_variable == v) {
	    for(i=0; i < c->Nparams; i++) {
	      if(Fullexpression->op1_variable == c->params[i])
		error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	    }
	    if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	      error(ERR_MEMALLOC);
	    Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	    Newexpression->op1type = VARTOOLS_OPERANDTYPE_VARIABLE;
	    Newexpression->op1_variable = Fullexpression->op1_variable;
	    return Newexpression;
	  } else {
	    return NULL;
	  }
	  break;
	case VARTOOLS_OPERANDTYPE_EXPRESSION:
	  foundit = 0;
	  for(i=0; i < c->Nparams; i++) {
	    if(Fullexpression->op1_variable == c->params[i])
	      foundit = 1;
	  }
	  if(foundit) {
	    if(CheckExpressionForVariables((_Expression *) Fullexpression->op2_expression, c->Nparams, c->params)){
	      error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	    }
	    else 
	      return NULL;
	  }
	  else {
	    Newexpression = GetLinfitVariableExpression(v, p, c, 
							(_Expression *) Fullexpression->op2_expression);
	    if(Newexpression != NULL) {
	      if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
		error(ERR_MEMALLOC);
	      Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	      Newexpression2->op1type = VARTOOLS_OPERANDTYPE_VARIABLE;
	      Newexpression2->op1_variable = Fullexpression->op1_variable;
	      Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	      Newexpression2->op2_expression = (void *) Newexpression;
	      return Newexpression2;
	    } else {
	      return NULL;
	    }
	  }
	  break;
	case VARTOOLS_OPERANDTYPE_FUNCTION:
	  if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	  }
	  return NULL;
	  break;
	default:
	  return NULL;
	  break;
	}
      }
      break;
    case VARTOOLS_OPERANDTYPE_CONSTANT:
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_CONSTANT:
	return NULL;
	break;
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	if(Fullexpression->op2_variable == v) {
	  if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  Newexpression->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	  Newexpression->op1_constant = Fullexpression->op1_constant;
	  return Newexpression;
	} else {
	  return NULL;
	}
	break;
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	
	Newexpression = GetLinfitVariableExpression(v, p, c, 
						    (_Expression *) Fullexpression->op2_expression);
	if(Newexpression != NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	  Newexpression2->op1_constant = Fullexpression->op1_constant;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op2_expression = (void *) Newexpression;
	  return Newexpression2;
	} else {
	  return NULL;
	}
	break;
      case VARTOOLS_OPERANDTYPE_FUNCTION:
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	return NULL;
	break;
      default:
	return NULL;
	break;
      }
      break;
    case VARTOOLS_OPERANDTYPE_EXPRESSION:
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_CONSTANT:
	Newexpression = GetLinfitVariableExpression(v, p, c, 
						    (_Expression *) Fullexpression->op1_expression);
	if(Newexpression != NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	  Newexpression2->op1_constant = Fullexpression->op2_constant;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op2_expression = (void *) Newexpression;
	  return Newexpression2;
	} else {
	  return NULL;
	}
	break;
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	foundit = 0;
	if(Fullexpression->op2_variable == v) {
	  foundit = 2;
	} else {
	  for(i=0; i < c->Nparams; i++) {
	    if(Fullexpression->op2_variable == c->params[i])
	      foundit = 1;
	  }
	}
	if(foundit == 2) {
	  if(CheckExpressionForVariables((_Expression *) Fullexpression->op1_expression, c->Nparams, c->params)){
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	  }
	  else {
	    return(Fullexpression->op1_expression);
	  }
	} else if(foundit == 1) {
	  if(CheckExpressionForVariables((_Expression *) Fullexpression->op1_expression, c->Nparams, c->params)){
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	  }
	  else
	    return NULL;
	}
	else {
	  Newexpression = GetLinfitVariableExpression(v, p, c, 
						      (_Expression *) Fullexpression->op1_expression);
	  if(Newexpression != NULL) {
	    if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	      error(ERR_MEMALLOC);
	    Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	    Newexpression2->op1type = VARTOOLS_OPERANDTYPE_VARIABLE;
	    Newexpression2->op1_variable = Fullexpression->op2_variable;
	    Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	    Newexpression2->op2_expression = Newexpression;
	    return Newexpression2;
	  } else {
	    return NULL;
	  }
	}
	break;

      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	Newexpression = GetLinfitVariableExpression(v, p, c, 
						    (_Expression *) Fullexpression->op1_expression);
	if(Newexpression != NULL) {
	  if(CheckExpressionForVariables((_Expression *) Fullexpression->op2_expression, c->Nparams, c->params)){
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	  }
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op1_expression = Newexpression;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op2_expression = Fullexpression->op2_expression;
	  return Newexpression2;
	} else {
	  Newexpression = GetLinfitVariableExpression(v, p, c, 
						      (_Expression *) Fullexpression->op2_expression);
	  if(Newexpression != NULL) {
	    if(CheckExpressionForVariables((_Expression *) Fullexpression->op1_expression, c->Nparams, c->params)){
	      error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	    }
	    if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	      error(ERR_MEMALLOC);
	    Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	    Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	    Newexpression2->op1_expression = Newexpression;
	    Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	    Newexpression2->op2_expression = Fullexpression->op1_expression;
	    return Newexpression2;
	  } else {
	    return NULL;
	  }
	}
	break;
      case VARTOOLS_OPERANDTYPE_FUNCTION:
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	Newexpression = GetLinfitVariableExpression(v, p, c, 
						    (_Expression *) Fullexpression->op1_expression);
	if(Newexpression != NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op1_expression = Newexpression;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_FUNCTION;
	  Newexpression2->op2_functioncall = Fullexpression->op2_functioncall;
	  return Newexpression2;
	} else {
	  return NULL;
	}
	break;
      default:
	Newexpression = GetLinfitVariableExpression(v, p, c, 
						    (_Expression *) Fullexpression->op1_expression);
	if(Newexpression != NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op1_expression = Newexpression;
	  Newexpression2->op2type = Fullexpression->op2type;
	  return Newexpression2;
	} else {
	  return NULL;
	}
	break;
      }
    case VARTOOLS_OPERANDTYPE_FUNCTION:
      if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op1_functioncall, c->params, c->Nparams)) {
	error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
      }
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	if(Fullexpression->op2_variable == v) {
	  if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  Newexpression->op1type = VARTOOLS_OPERANDTYPE_FUNCTION;
	  Newexpression->op1_functioncall = Fullexpression->op1_functioncall;
	  return Newexpression;
	} else {
	  return NULL;
	}
	break;
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	Newexpression = GetLinfitVariableExpression(v, p, c, 
						    (_Expression *) Fullexpression->op2_expression);
	if(Newexpression != NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_FUNCTION;
	  Newexpression2->op1_functioncall = Fullexpression->op1_functioncall;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op2_expression = Newexpression;
	  return Newexpression2;
	} else {
	  return NULL;
	}
	break;
      case VARTOOLS_OPERANDTYPE_FUNCTION:
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	return NULL;
	break;
      default:
	return NULL;
	break;
      }
      break;
    default:
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	if(Fullexpression->op2_variable == v) {
	  if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  Newexpression->op1type = Fullexpression->op1type;
	  return Newexpression;
	} else {
	  return NULL;
	}
	break;
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	Newexpression = GetLinfitVariableExpression(v, p, c, 
						    (_Expression *) Fullexpression->op2_expression);
	if(Newexpression != NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression2->op1type = Fullexpression->op1type;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op2_expression = Newexpression;
	  return Newexpression2;
	} else {
	  return NULL;
	}
	break;
      case VARTOOLS_OPERANDTYPE_FUNCTION:
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	return NULL;
	break;
      default:
	return NULL;
	break;
      }
      break;
    }
  case VARTOOLS_OPERATORTYPE_DIVIDE:
    if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_VARIABLE) {
      for(i=0; i < c->Nparams; i++) {
	if(Fullexpression->op2_variable == c->params[i])
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
      }
    }
    else if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_EXPRESSION) {
      if(CheckExpressionForVariables((_Expression *) Fullexpression->op2_expression, c->Nparams, c->params)){
	error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
      }
    }
    else if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_FUNCTION) {
      if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
      }
    }
    if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_CONSTANT ||
       Fullexpression->op1type == VARTOOLS_OPERANDTYPE_ITERATORNR ||
       Fullexpression->op1type == VARTOOLS_OPERANDTYPE_ITERATORNF) {
      return NULL;
    }
    if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
      error(ERR_MEMALLOC);
    Newexpression->operatortype = VARTOOLS_OPERATORTYPE_DIVIDE;
    Newexpression->op2type = Fullexpression->op2type;
    if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_VARIABLE) {
      Newexpression->op2_variable = Fullexpression->op2_variable;
    }
    else if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_EXPRESSION) {
      Newexpression->op2_expression = Fullexpression->op2_expression;
    }
    else if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_FUNCTION) {
      Newexpression->op2_functioncall = Fullexpression->op2_functioncall;
    }
    if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_VARIABLE) {
      if(Fullexpression->op1_variable == v) {
	Newexpression->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	Newexpression->op1_constant = 1.0;
	return Newexpression;
      } else {
	free(Newexpression);
	return NULL;
      }
    }
    else if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_EXPRESSION) {
      Newexpression2 = GetLinfitVariableExpression(v, p, c, 
						   (_Expression *) Fullexpression->op1_expression);
      if(Newexpression2 != NULL) {
	Newexpression->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	Newexpression->op1_expression = Newexpression2;
	return Newexpression;
      } else {
	free(Newexpression);
	return NULL;
      }
    }
    else if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_FUNCTION) {
      if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op1_functioncall, c->params, c->Nparams)) {
	error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
      }
      free(Newexpression);
      return NULL;
    }
    return NULL;
    break;
  default:
    if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_VARIABLE){
      for(i=0; i < c->Nparams; i++) {
	if(Fullexpression->op1_variable == c->params[i])
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
      }
    }
    else if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_EXPRESSION){
      if(CheckExpressionForVariables((_Expression *) Fullexpression->op1_expression, c->Nparams, c->params)){
	error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
      }
    }
    else if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_FUNCTION){
      if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op1_functioncall, c->params, c->Nparams)) {
	error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
      }
    }
    if(Fullexpression->operatortype != VARTOOLS_OPERATORTYPE_NOT) {
      if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_VARIABLE){
	for(i=0; i < c->Nparams; i++) {
	  if(Fullexpression->op2_variable == c->params[i])
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
      }
      else if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_EXPRESSION){
	if(CheckExpressionForVariables((_Expression *) Fullexpression->op2_expression, c->Nparams, c->params)){
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
      }
      else if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_FUNCTION){
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
      }
    }
    return NULL;
    break;
  }
}


/* Get the expression for the constant term in a linear fit */
_Expression *GetLinfitConstantExpression(ProgramData *p,
					 _Linfit *c, 
					 _Expression *Fullexpression) {
  int foundit = 0;
  int i;
  _Expression *Newexpression;
  _Expression *Newexpression2;
  _Expression *Newexpression3;
  switch(Fullexpression->operatortype) {
  case VARTOOLS_OPERATORTYPE_CONSTANT:
    switch(Fullexpression->op1type) {
    case VARTOOLS_OPERANDTYPE_VARIABLE:
      foundit = 0;
      for(i=0; i < c->Nparams; i++) {
	if(Fullexpression->op1_variable == c->params[i])
	  foundit = 1;
      }
      if(!foundit)
	return Fullexpression;
      else
	return NULL;
      break;
    case VARTOOLS_OPERANDTYPE_EXPRESSION:
      Newexpression = GetLinfitConstantExpression(p, c, (_Expression *) Fullexpression->op1_expression);
      return Newexpression;
      break;
    case VARTOOLS_OPERANDTYPE_FUNCTION:
      if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op1_functioncall, c->params, c->Nparams)) {
	error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
      }
      return Fullexpression;
      break;
    default:
      return Fullexpression;
      break;
    }
    break;
  case VARTOOLS_OPERATORTYPE_ADD:
    switch(Fullexpression->op1type) {
    case VARTOOLS_OPERANDTYPE_VARIABLE:
      foundit = 0;
      for(i=0; i < c->Nparams; i++) {
	if(Fullexpression->op1_variable == c->params[i])
	  foundit = 1;
      }
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_CONSTANT:
	if(!foundit)
	  return Fullexpression;
	else {
	  if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  Newexpression->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	  Newexpression->op1_constant = Fullexpression->op2_constant;
	  return Newexpression;
	}
	break;
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	for(i=0; i < c->Nparams; i++) {
	  if(Fullexpression->op2_variable == c->params[i]) {
	    if(foundit)
	      return NULL;
	    else {
	      if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
		error(ERR_MEMALLOC);
	      Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	      Newexpression->op1type = VARTOOLS_OPERANDTYPE_VARIABLE;
	      Newexpression->op1_variable = Fullexpression->op1_variable;
	      return Newexpression;
	    }
	  }
	}
	if(!foundit)
	  return Fullexpression;
	else {
	  if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  Newexpression->op1type = VARTOOLS_OPERANDTYPE_VARIABLE;
	  Newexpression->op1_variable = Fullexpression->op2_variable;
	  return Newexpression;
	}
	break;
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	Newexpression = GetLinfitConstantExpression(p, c, (_Expression *) Fullexpression->op2_expression);
	if(Newexpression == NULL) {
	  if(foundit)
	    return NULL;
	  else {
	    if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	      error(ERR_MEMALLOC);
	    Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	    Newexpression2->op1type = VARTOOLS_OPERANDTYPE_VARIABLE;
	    Newexpression2->op1_variable = Fullexpression->op1_variable;
	    return Newexpression2;
	  }
	} 
	else {
	  if(foundit) {
	    return Newexpression;
	  }
	  else {
	    if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	      error(ERR_MEMALLOC);
	    Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_ADD;
	    Newexpression2->op1type = VARTOOLS_OPERANDTYPE_VARIABLE;
	    Newexpression2->op1_variable = Fullexpression->op1_variable;
	    Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	    Newexpression2->op1_expression = (void *) Newexpression;
	    return Newexpression2;
	  }
	}
	break;
      case VARTOOLS_OPERANDTYPE_FUNCTION:
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	if(!foundit)
	  return Fullexpression;
	else {
	  if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  Newexpression->op1type = VARTOOLS_OPERANDTYPE_FUNCTION;
	  Newexpression->op1_functioncall = Fullexpression->op2_functioncall;
	  return Newexpression;
	}
	break;
      default:
	if(!foundit)
	  return Fullexpression;
	else {
	  if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  Newexpression->op1type = Fullexpression->op2type;
	  return Newexpression;
	}
	break;
      }
      break;
    case VARTOOLS_OPERANDTYPE_EXPRESSION:
      Newexpression = GetLinfitConstantExpression(p, c, (_Expression *) Fullexpression->op1_expression);
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_CONSTANT:
	if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	  error(ERR_MEMALLOC);
	if(Newexpression != NULL) {
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_ADD;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op1_expression = (void *) Newexpression;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_CONSTANT;
	  Newexpression2->op2_constant = Fullexpression->op2_constant;
	  return Newexpression2;
	} else {
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	  Newexpression2->op1_constant = Fullexpression->op1_constant;
	  return Newexpression2;
	}
	break;
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	foundit = 0;
	for(i=0; i < c->Nparams; i++) {
	  if(Fullexpression->op2_variable == c->params[i])
	    foundit = 1;
	}
	if(foundit) {
	  return Newexpression;
	}
	else {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  if(Newexpression == NULL) {
	    Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	    Newexpression2->op1type = VARTOOLS_OPERANDTYPE_VARIABLE;
	    Newexpression2->op1_variable = Fullexpression->op2_variable;
	    return Newexpression2;
	  } else {
	    Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_ADD;
	    Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	    Newexpression2->op1_expression = (void *) Newexpression;
	    Newexpression2->op2type = VARTOOLS_OPERANDTYPE_VARIABLE;
	    Newexpression2->op2_variable = Fullexpression->op2_variable;
	    return Newexpression2;
	  }
	}
	break;
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	Newexpression2 = GetLinfitConstantExpression(p, c, (_Expression *) Fullexpression->op2_expression); 
	if(Newexpression == NULL)
	  return Newexpression2;
	if(Newexpression2 == NULL)
	  return Newexpression;
	if((Newexpression3 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	  error(ERR_MEMALLOC);
	Newexpression3->operatortype = VARTOOLS_OPERATORTYPE_ADD;
	Newexpression3->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	Newexpression3->op1_expression = (void *) Newexpression;
	Newexpression3->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	Newexpression3->op2_expression = (void *) Newexpression2;
	return Newexpression3;
	break;
      case VARTOOLS_OPERANDTYPE_FUNCTION:
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	if(Newexpression == NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_FUNCTION;
	  Newexpression2->op1_functioncall = Fullexpression->op2_functioncall;
	  return Newexpression2;
	} 
	else {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_ADD;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op1_expression = (void *) Newexpression;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_FUNCTION;
	  Newexpression2->op2_functioncall = Fullexpression->op2_functioncall;
	  return Newexpression2;
	}
	break;
      default:
	if(Newexpression == NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  Newexpression2->op1type = Fullexpression->op2type;
	  return Newexpression2;
	}
	else {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_ADD;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op1_expression = (void *) Newexpression;
	  Newexpression2->op2type = Fullexpression->op2type;
	  return Newexpression2;
	}
	break;
      }
    default:
      if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_FUNCTION) {
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op1_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
      }
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	foundit = 0;
	for(i=0; i < c->Nparams; i++) {
	  if(Fullexpression->op2_variable == c->params[i])
	    foundit = 1;
	}
	if(!foundit) {
	  return Fullexpression;
	} else {
	  if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  Newexpression->op1type = Fullexpression->op1type;
	  if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_CONSTANT) {
	    Newexpression->op1_constant = Fullexpression->op1_constant;
	  }
	  else if(Fullexpression->op1type = VARTOOLS_OPERANDTYPE_FUNCTION) {
	    Newexpression->op1_functioncall = Fullexpression->op1_functioncall;
	  }
	  return Newexpression;
	}
	break;
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	Newexpression = GetLinfitConstantExpression(p, c, (_Expression *) Fullexpression->op2_expression); 
	if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	  error(ERR_MEMALLOC);
	Newexpression2->op1type = Fullexpression->op1type;
	if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_CONSTANT) {
	  Newexpression2->op1_constant = Fullexpression->op1_constant;
	}
	else if(Fullexpression->op1type = VARTOOLS_OPERANDTYPE_FUNCTION) {
	  Newexpression2->op1_functioncall = Fullexpression->op1_functioncall;
	}
	if(Newexpression == NULL) {
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  return Newexpression2;
	}
	else {
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_ADD;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op2_expression = Newexpression;
	  return Newexpression2;
	}
	break;
      default:
	if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_FUNCTION) {
	  if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	  }
	}
	return Fullexpression;
	break;
      }
      break;
    }


  case VARTOOLS_OPERATORTYPE_SUBTRACT:
    switch(Fullexpression->op1type) {
    case VARTOOLS_OPERANDTYPE_VARIABLE:
      foundit = 0;
      for(i=0; i < c->Nparams; i++) {
	if(Fullexpression->op1_variable == c->params[i])
	  foundit = 1;
      }
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_CONSTANT:
	if(!foundit)
	  return Fullexpression;
	else {
	  if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  Newexpression->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	  Newexpression->op1_constant = -(Fullexpression->op2_constant);
	  return Newexpression;
	}
	break;
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	for(i=0; i < c->Nparams; i++) {
	  if(Fullexpression->op2_variable == c->params[i]) {
	    if(foundit)
	      return NULL;
	    else {
	      if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
		error(ERR_MEMALLOC);
	      Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	      Newexpression->op1type = VARTOOLS_OPERANDTYPE_VARIABLE;
	      Newexpression->op1_variable = Fullexpression->op1_variable;
	      return Newexpression;
	    }
	  }
	}
	if(!foundit)
	  return Fullexpression;
	else {
	  if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	  Newexpression->op1_constant = -1.0;
	  Newexpression->op2type = VARTOOLS_OPERANDTYPE_VARIABLE;
	  Newexpression->op2_variable = Fullexpression->op2_variable;
	  return Newexpression;
	}
	break;
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	Newexpression = GetLinfitConstantExpression(p, c, (_Expression *) Fullexpression->op2_expression);
	if(Newexpression == NULL) {
	  if(foundit)
	    return NULL;
	  else {
	    if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	      error(ERR_MEMALLOC);
	    Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	    Newexpression2->op1type = VARTOOLS_OPERANDTYPE_VARIABLE;
	    Newexpression2->op1_variable = Fullexpression->op1_variable;
	    return Newexpression2;
	  }
	} 
	else {
	  if(foundit) {
	    if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	      error(ERR_MEMALLOC);
	    Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	    Newexpression2->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	    Newexpression2->op1_constant = -1.0;
	    Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	    Newexpression2->op2_expression = (void *) Newexpression;
	    return Newexpression2;
	  }
	  else {
	    if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	      error(ERR_MEMALLOC);
	    Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_SUBTRACT;
	    Newexpression2->op1type = VARTOOLS_OPERANDTYPE_VARIABLE;
	    Newexpression2->op1_variable = Fullexpression->op1_variable;
	    Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	    Newexpression2->op1_expression = (void *) Newexpression;
	    return Newexpression2;
	  }
	}
	break;
      case VARTOOLS_OPERANDTYPE_FUNCTION:
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	if(!foundit)
	  return Fullexpression;
	else {
	  if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	  Newexpression->op1_constant = -1.0;
	  Newexpression->op2type = VARTOOLS_OPERANDTYPE_FUNCTION;
	  Newexpression->op2_functioncall = Fullexpression->op2_functioncall;
	  return Newexpression;
	}
	break;
      default:
	if(!foundit)
	  return Fullexpression;
	else {
	  if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	  Newexpression->op1_constant = -1.0;
	  Newexpression->op2type = Fullexpression->op2type;
	  return Newexpression;
	}
	break;
      }
      break;
    case VARTOOLS_OPERANDTYPE_EXPRESSION:
      Newexpression = GetLinfitConstantExpression(p, c, (_Expression *) Fullexpression->op1_expression);
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_CONSTANT:
	if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	  error(ERR_MEMALLOC);
	if(Newexpression != NULL) {
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_SUBTRACT;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op1_expression = (void *) Newexpression;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_CONSTANT;
	  Newexpression2->op2_constant = Fullexpression->op2_constant;
	  return Newexpression2;
	} else {
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	  Newexpression2->op1_constant = Fullexpression->op1_constant;
	  return Newexpression2;
	}
	break;
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	foundit = 0;
	for(i=0; i < c->Nparams; i++) {
	  if(Fullexpression->op2_variable == c->params[i])
	    foundit = 1;
	}
	if(foundit) {
	  return Newexpression;
	}
	else {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  if(Newexpression == NULL) {
	    Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	    Newexpression2->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	    Newexpression2->op1_constant = -1.0;
	    Newexpression2->op2type = VARTOOLS_OPERANDTYPE_VARIABLE;
	    Newexpression2->op2_variable = Fullexpression->op2_variable;
	    return Newexpression2;
	  } else {
	    Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_SUBTRACT;
	    Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	    Newexpression2->op1_expression = (void *) Newexpression;
	    Newexpression2->op2type = VARTOOLS_OPERANDTYPE_VARIABLE;
	    Newexpression2->op2_variable = Fullexpression->op2_variable;
	    return Newexpression2;
	  }
	}
	break;
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	Newexpression2 = GetLinfitConstantExpression(p, c, (_Expression *) Fullexpression->op2_expression); 
	if(Newexpression == NULL) {
	  if(Newexpression2 == NULL)
	    return NULL;
	  if((Newexpression3 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression3->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression3->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	  Newexpression3->op1_constant = -1.0;
	  Newexpression3->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression3->op2_expression = (void *) Newexpression2;
	  return Newexpression3;
	}
	if(Newexpression2 == NULL)
	  return Newexpression;
	if((Newexpression3 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	  error(ERR_MEMALLOC);
	Newexpression3->operatortype = VARTOOLS_OPERATORTYPE_SUBTRACT;
	Newexpression3->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	Newexpression3->op1_expression = (void *) Newexpression;
	Newexpression3->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	Newexpression3->op2_expression = (void *) Newexpression2;
	return Newexpression3;
	break;
      case VARTOOLS_OPERANDTYPE_FUNCTION:
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	if(Newexpression == NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	  Newexpression2->op1_constant = -1.0;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_FUNCTION;
	  Newexpression2->op2_functioncall = Fullexpression->op2_functioncall;
	  return Newexpression2;
	} 
	else {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_SUBTRACT;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op1_expression = (void *) Newexpression;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_FUNCTION;
	  Newexpression2->op2_functioncall = Fullexpression->op2_functioncall;
	  return Newexpression2;
	}
	break;
      default:
	if(Newexpression == NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	  Newexpression2->op1_constant = -1.0;
	  Newexpression2->op2type = Fullexpression->op2type;
	  return Newexpression2;
	}
	else {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_SUBTRACT;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op1_expression = (void *) Newexpression;
	  Newexpression2->op2type = Fullexpression->op2type;
	  return Newexpression2;
	}
	break;
      }
    default:
      if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_FUNCTION) {
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op1_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
      }
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	foundit = 0;
	for(i=0; i < c->Nparams; i++) {
	  if(Fullexpression->op2_variable == c->params[i])
	    foundit = 1;
	}
	if(!foundit) {
	  return Fullexpression;
	} else {
	  if((Newexpression = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  Newexpression->op1type = Fullexpression->op1type;
	  if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_CONSTANT) {
	    Newexpression->op1_constant = Fullexpression->op1_constant;
	  }
	  else if(Fullexpression->op1type = VARTOOLS_OPERANDTYPE_FUNCTION) {
	    Newexpression->op1_functioncall = Fullexpression->op1_functioncall;
	  }
	  return Newexpression;
	}
	break;
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	Newexpression = GetLinfitConstantExpression(p, c, (_Expression *) Fullexpression->op2_expression); 
	if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	  error(ERR_MEMALLOC);
	Newexpression2->op1type = Fullexpression->op1type;
	if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_CONSTANT) {
	  Newexpression2->op1_constant = Fullexpression->op1_constant;
	}
	else if(Fullexpression->op1type = VARTOOLS_OPERANDTYPE_FUNCTION) {
	  Newexpression2->op1_functioncall = Fullexpression->op1_functioncall;
	}
	if(Newexpression == NULL) {
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_CONSTANT;
	  return Newexpression2;
	}
	else {
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_SUBTRACT;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op2_expression = Newexpression;
	  return Newexpression2;
	}
	break;
      default:
	if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_FUNCTION) {
	  if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	  }
	}
	return Fullexpression;
	break;
      }
      break;
    }
    break;
  case VARTOOLS_OPERATORTYPE_MULTIPLY:
    switch(Fullexpression->op1type) {
    case VARTOOLS_OPERANDTYPE_CONSTANT:
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	for(i=0; i < c->Nparams; i++) {
	  if(Fullexpression->op2_variable == c->params[i])
	    return NULL;
	}
	return Fullexpression;
	break;
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	Newexpression = GetLinfitConstantExpression(p, c, (_Expression *) Fullexpression->op2_expression); 
	if(Newexpression == NULL) {
	  return NULL;
	}
	else {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_CONSTANT;
	  Newexpression2->op1_constant = Fullexpression->op1_constant;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op2_expression = (void *) Newexpression;
	  return Newexpression2;
	}
	break;
      default:
	if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_FUNCTION) {
	  if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	  }
	}
	return Fullexpression;
	break;
      }
      break;
    case VARTOOLS_OPERANDTYPE_VARIABLE:
      foundit = 0;
      for(i=0; i < c->Nparams; i++) {
	if(Fullexpression->op1_variable == c->params[i])
	  foundit = 1;
      }
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	for(i=0; i < c->Nparams; i++) {
	  if(Fullexpression->op2_variable == c->params[i]) {
	    if(foundit)
	      error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	    else
	      return NULL;
	  }
	}
	if(foundit)
	  return NULL;
	else
	  return Fullexpression;
	break;
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	if(foundit) {
	  if(CheckExpressionForVariables((_Expression *) Fullexpression->op2_expression, c->Nparams, c->params)){
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	  }
	  return NULL;
	} else {
	  Newexpression = GetLinfitConstantExpression(p, c, (_Expression *) Fullexpression->op2_expression); 
	  if(Newexpression == NULL)
	    return NULL;
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_VARIABLE;
	  Newexpression2->op1_variable = Fullexpression->op1_variable;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op2_expression = Newexpression;
	  return Newexpression2;
	}
	break;
      case VARTOOLS_OPERANDTYPE_FUNCTION:
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	if(foundit)
	  return NULL;
	else
	  return Fullexpression;
	break;
      default:
	if(foundit)
	  return NULL;
	else
	  return Fullexpression;
	break;
      }
      break;
    case VARTOOLS_OPERANDTYPE_EXPRESSION:
      Newexpression = GetLinfitConstantExpression(p, c, (_Expression *) Fullexpression->op1_expression); 
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	foundit = CheckExpressionForVariables((_Expression *) Fullexpression->op1_expression, c->Nparams, c->params);
	if(foundit) {
	  if(CheckExpressionForVariables((_Expression *) Fullexpression->op2_expression, c->Nparams, c->params))
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	  if(Newexpression != NULL) {
	    if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	      error(ERR_MEMALLOC);
	    Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	    Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	    Newexpression2->op1_expression = (void *) Newexpression;
	    Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	    Newexpression2->op2_expression = (void *) Fullexpression->op2_expression;
	    return Newexpression2;
	  } else {
	    return NULL;
	  }
	} else {
	  if(Newexpression != NULL) {
	    Newexpression2 = GetLinfitConstantExpression(p, c, (_Expression *) Fullexpression->op2_expression);
	    if(Newexpression2 != NULL) {
	      if((Newexpression3 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
		error(ERR_MEMALLOC);
	      Newexpression3->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	      Newexpression3->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	      Newexpression3->op1_expression = (void *) Newexpression;
	      Newexpression3->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	      Newexpression3->op2_expression = (void *) Newexpression2;
	      return Newexpression3;
	    } else {
	      return NULL;
	    }
	  } else {
	    return NULL;
	  }
	}
	break;
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	foundit = 0;
	for(i=0; i < c->Nparams; i++) {
	  if(Fullexpression->op2_variable == c->params[i]) {
	    if(CheckExpressionForVariables((_Expression *) Fullexpression->op1_expression, c->Nparams, c->params))
	      error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	    else
	      return NULL;
	  }
	}
	if(Newexpression != NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op1_expression = (void *) Newexpression;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_VARIABLE;
	  Newexpression2->op2_variable = Fullexpression->op2_variable;
	  return Newexpression2;
	  break;
	}
	else
	  return NULL;
	break;
      case VARTOOLS_OPERANDTYPE_FUNCTION:
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	if(Newexpression != NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op1_expression = (void *) Newexpression;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_FUNCTION;
	  Newexpression2->op2_functioncall = Fullexpression->op2_functioncall;
	  return Newexpression2;
	} else {
	  return NULL;
	}
	break;
      default:
	if(Newexpression != NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op1_expression = (void *) Newexpression;
	  Newexpression2->op2type = Fullexpression->op2type;
	  if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_CONSTANT) {
	    Newexpression2->op2_constant = Fullexpression->op2_constant;
	  }
	  return Newexpression2;
	} else {
	  return NULL;
	}
	break;
      }
      break;
    default:
      if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_FUNCTION) {
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op1_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
      }
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	for(i=0; i < c->Nparams; i++) {
	  if(Fullexpression->op1_variable == c->params[i])
	    return NULL;
	}
	return Fullexpression;
	break;
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	Newexpression = GetLinfitConstantExpression(p, c, (_Expression *) Fullexpression->op2_expression);
	if(Newexpression == NULL)
	  return NULL;
	if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	  error(ERR_MEMALLOC);
	Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_MULTIPLY;
	Newexpression2->op1type = Fullexpression->op1type;
	if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_FUNCTION)
	  Newexpression2->op1_functioncall = Fullexpression->op1_functioncall;
	else if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_CONSTANT)
	  Newexpression2->op1_constant = Fullexpression->op1_constant;
	Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	Newexpression2->op2_expression = (void *) Newexpression;
	return Newexpression2;
	break;
      case VARTOOLS_OPERANDTYPE_FUNCTION:
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	return Fullexpression;
	break;
      default:
	return Fullexpression;
	break;
      }
      break;
    }
    break;

  case VARTOOLS_OPERATORTYPE_DIVIDE:
    switch(Fullexpression->op1type) {
    case VARTOOLS_OPERANDTYPE_CONSTANT:
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	for(i=0; i < c->Nparams; i++) {
	  if(Fullexpression->op2_variable == c->params[i])
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	return Fullexpression;
	break;
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	if(CheckExpressionForVariables((_Expression *) Fullexpression->op2_expression, c->Nparams, c->params))
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	return Fullexpression;
	break;
      default:
	if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_FUNCTION) {
	  if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	  }
	}
	return Fullexpression;
	break;
      }
      break;
    case VARTOOLS_OPERANDTYPE_VARIABLE:
      foundit = 0;
      for(i=0; i < c->Nparams; i++) {
	if(Fullexpression->op1_variable == c->params[i])
	  foundit = 1;
      }
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	for(i=0; i < c->Nparams; i++) {
	  if(Fullexpression->op2_variable == c->params[i]) {
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	  }
	}
	if(foundit)
	  return NULL;
	else
	  return Fullexpression;
	break;
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	if(CheckExpressionForVariables((_Expression *) Fullexpression->op2_expression, c->Nparams, c->params)){
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	if(foundit) {
	  return NULL;
	} else {
	  return Fullexpression;
	}
	break;
      case VARTOOLS_OPERANDTYPE_FUNCTION:
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	if(foundit)
	  return NULL;
	else
	  return Fullexpression;
	break;
      default:
	if(foundit)
	  return NULL;
	else
	  return Fullexpression;
	break;
      }
      break;
    case VARTOOLS_OPERANDTYPE_EXPRESSION:
      Newexpression = GetLinfitConstantExpression(p, c, (_Expression *) Fullexpression->op1_expression); 
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	if(CheckExpressionForVariables((_Expression *) Fullexpression->op2_expression, c->Nparams, c->params))
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	if(Newexpression != NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_DIVIDE;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op1_expression = (void *) Newexpression;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op2_expression = (void *) Fullexpression->op2_expression;
	  return Newexpression2;
	} else {
	  return NULL;
	}
	break;
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	for(i=0; i < c->Nparams; i++) {
	  if(Fullexpression->op2_variable == c->params[i])
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	if(Newexpression != NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_DIVIDE;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op1_expression = (void *) Newexpression;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_VARIABLE;
	  Newexpression2->op2_variable = Fullexpression->op2_variable;
	  return Newexpression2;
	  break;
	}
	else
	  return NULL;
	break;
      case VARTOOLS_OPERANDTYPE_FUNCTION:
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	if(Newexpression != NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_DIVIDE;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op1_expression = (void *) Newexpression;
	  Newexpression2->op2type = VARTOOLS_OPERANDTYPE_FUNCTION;
	  Newexpression2->op2_functioncall = Fullexpression->op2_functioncall;
	  return Newexpression2;
	} else {
	  return NULL;
	}
	break;
      default:
	if(Newexpression != NULL) {
	  if((Newexpression2 = (_Expression *) malloc(sizeof(_Expression))) == NULL)
	    error(ERR_MEMALLOC);
	  Newexpression2->operatortype = VARTOOLS_OPERATORTYPE_DIVIDE;
	  Newexpression2->op1type = VARTOOLS_OPERANDTYPE_EXPRESSION;
	  Newexpression2->op1_expression = (void *) Newexpression;
	  Newexpression2->op2type = Fullexpression->op2type;
	  if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_CONSTANT) {
	    Newexpression2->op2_constant = Fullexpression->op2_constant;
	  }
	  return Newexpression2;
	}
	else 
	  return NULL;
	break;
      }
      break;
    default:
      if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_FUNCTION) {
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op1_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
      }
      switch(Fullexpression->op2type) {
      case VARTOOLS_OPERANDTYPE_VARIABLE:
	for(i=0; i < c->Nparams; i++) {
	  if(Fullexpression->op1_variable == c->params[i])
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	return Fullexpression;
	break;
      case VARTOOLS_OPERANDTYPE_EXPRESSION:
	if(CheckExpressionForVariables((_Expression *) Fullexpression->op2_expression, c->Nparams, c->params))
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	return Fullexpression;
	break;
      case VARTOOLS_OPERANDTYPE_FUNCTION:
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
	return Fullexpression;
	break;
      default:
	return Fullexpression;
	break;
      }
      break;
    }
    break;
    
  default:
    if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_VARIABLE){
      for(i=0; i < c->Nparams; i++) {
	if(Fullexpression->op1_variable == c->params[i])
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
      }
    }
    else if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_EXPRESSION){
      if(CheckExpressionForVariables((_Expression *) Fullexpression->op1_expression, c->Nparams, c->params)){
	error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
      }
    }
    else if(Fullexpression->op1type == VARTOOLS_OPERANDTYPE_FUNCTION){
      if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op1_functioncall, c->params, c->Nparams)) {
	error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
      }
    }
    if(Fullexpression->operatortype != VARTOOLS_OPERATORTYPE_NOT) {
      if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_VARIABLE){
	for(i=0; i < c->Nparams; i++) {
	  if(Fullexpression->op2_variable == c->params[i])
	    error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
      }
      else if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_EXPRESSION){
	if(CheckExpressionForVariables((_Expression *) Fullexpression->op2_expression, c->Nparams, c->params)){
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
      }
      else if(Fullexpression->op2type == VARTOOLS_OPERANDTYPE_FUNCTION){
	if(CheckFunctionCallForVariables((_FunctionCall *) Fullexpression->op2_functioncall, c->params, c->Nparams)) {
	  error2(ERR_INVALIDFUNCTIONLINFIT,c->functionstring);
	}
      }
    }
    return Fullexpression;
    break;
  }
}

void SetupLinfitExpression(ProgramData *p, _Linfit *c)
{
  _Expression *Fullexpression;
  int i = 0;

  Fullexpression = ParseExpression(c->functionstring, p);
  
  for(i=0; i < c->Nparams; i++)
    {
      c->expressions[i] = GetLinfitVariableExpression(c->params[i],p,c, 
						      Fullexpression);
      if(c->expressions[i] == NULL) {
	error2(ERR_LINFITMISSINGPARAM, c->params[i]->varname);
      }
    }
  c->constantexpression = GetLinfitConstantExpression(p, c, Fullexpression);

}

int CheckFunctionCallForVariables(_FunctionCall *f, _Variable **vars, int Nvar) {
  int N, i;
  N = f->Nexpr;
  for(i = 0 ; i < N; i++) {
    if(CheckExpressionForVariables(f->arguments[i], Nvar, vars))
      return 1;
  }
  return 0;
}

int CheckExpressionForVariables(_Expression *expression, int Nvar, 
				_Variable **vars) {
  int i;
  if(expression->operatortype == VARTOOLS_OPERATORTYPE_CONSTANT) {
    switch(expression->op1type) {
    case VARTOOLS_OPERANDTYPE_CONSTANT:
      return 0;
      break;
    case VARTOOLS_OPERANDTYPE_VARIABLE:
      for(i=0; i < Nvar; i++) {
	if(expression->op1_variable == vars[i])
	  return 1;
      }
      return 0;
      break;
    case VARTOOLS_OPERANDTYPE_EXPRESSION:
      return(CheckExpressionForVariables((_Expression *) expression->op1_expression, Nvar, vars));
      break;
    case VARTOOLS_OPERANDTYPE_FUNCTION:
      if(CheckFunctionCallForVariables((_FunctionCall *) expression->op1_functioncall, vars, Nvar))
	return 1;
      return 0;
      break;
    default:
      return 0;
    }
  } else {
    switch(expression->op1type) {
    case VARTOOLS_OPERANDTYPE_CONSTANT:
      break;
    case VARTOOLS_OPERANDTYPE_VARIABLE:
      for(i=0; i < Nvar; i++) {
	if(expression->op1_variable == vars[i])
	  return 1;
      }
      break;
    case VARTOOLS_OPERANDTYPE_EXPRESSION:
      if(CheckExpressionForVariables((_Expression *) expression->op1_expression, Nvar, vars))
	return 1;
      break;
    case VARTOOLS_OPERANDTYPE_FUNCTION:
      if(CheckFunctionCallForVariables((_FunctionCall *) expression->op1_functioncall, vars, Nvar))
	return 1;
      break;
    default:
      break;
    }
    switch(expression->op2type) {
    case VARTOOLS_OPERANDTYPE_CONSTANT:
      return 0;
      break;
    case VARTOOLS_OPERANDTYPE_VARIABLE:
      for(i=0; i < Nvar; i++) {
	if(expression->op2_variable == vars[i])
	  return 1;
      }
      return 0;
      break;
    case VARTOOLS_OPERANDTYPE_EXPRESSION:
      return(CheckExpressionForVariables((_Expression *) expression->op2_expression, Nvar, vars));
      break;
    case VARTOOLS_OPERANDTYPE_FUNCTION:
      if(CheckFunctionCallForVariables((_FunctionCall *) expression->op2_functioncall, vars, Nvar))
	return 1;
      return 0;
      break;
    default:
      return 0;
    }
  }
}

int ParseLinfitCommand(int *iret, int argc, char **argv, ProgramData *p,
		       _Linfit *c)
{
  int i, k, j, nparam, s;
  char oldval;

  c->outfile_extension = (char *) malloc(20);
  c->outfilename_format = (char *) malloc(1);
  sprintf(c->outfile_extension,"linfit.model");
  c->outfilename_format[0] = '\0';
  
  i = *iret;
  
  if(i >= argc)
    return(1);

  if((c->functionstring = (char *) malloc((strlen(argv[i])+1))) == NULL)
    error(ERR_MEMALLOC);
  sprintf(c->functionstring,"%s",argv[i]);
  i++;
  if(i >= argc) {*iret = i; return 1;}

  if((c->paramliststring = (char *) malloc((strlen(argv[i])+1))) == NULL)
    error(ERR_MEMALLOC);
  sprintf(c->paramliststring,"%s",argv[i]);
  c->modelvarname = NULL;
  c->calcchi2out = 0;
  c->correctlc = 0;
  c->omodel = 0;
  c->rejectoutliers = 0;
  c->rejsigclip = 0.0;
  c->rejuseMAD = 0;
  c->rejiterate = 0;
  c->rejfixnum = 0;
  c->rejiternum = 0;
  c->usemask = 0;
  c->maskvarname = NULL;
  c->maskvar = NULL;
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"modelvar")) {
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if((c->modelvarname= (char *) malloc((strlen(argv[i])+1))) == NULL)
	error(ERR_MEMALLOC);
      sprintf(c->modelvarname,"%s",argv[i]);
    } else
      i--;
  }
  else
    i--;
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"reject")) {
      c->rejectoutliers = 1;
      i++;
      if(i >= argc) {*iret = i; return 1;}
      c->rejsigclip = atof(argv[i]);
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"useMAD")) {
	  c->rejuseMAD = 1;
	} else
	  i--;
      }
      else
	i--;
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"iter")) {
	  c->rejiterate = 1;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"fixednum")) {
	      c->rejfixnum = 1;
	      i++;
	      if(i >= argc) {*iret = i; return 1;}
	      c->rejiternum = atoi(argv[i]);
	    } else
	      i--;
	  }
	  else
	    i--;
	} else
	  i--;
      }
      else
	i--;
    } else
      i--;
  }
  else
    i--;
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"correctlc")) {
      c->correctlc = 1;
    }
    else
      i--;
  }
  else
    i--;
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"omodel")) {
      c->omodel = 1;
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if((c->outdir= (char *) malloc((strlen(argv[i])+1))) == NULL)
	error(ERR_MEMALLOC);
      sprintf(c->outdir,"%s",argv[i]);
      /* Check if the user gave the "format" keyword */
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"format")) {
	  i++;
	  if(i < argc) {
	    c->outfilename_format = (char *) realloc(c->outfilename_format, (strlen(argv[i])+1)*sizeof(char));
	    sprintf(c->outfilename_format,"%s",argv[i]);
	    i++;
	  } else {
	    *iret = i; return 1;
	  }
	} else
	  i--;
      } else
	i--;
    } else
      i--;
  }
  else
    i--;

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"fitmask")) {
      c->usemask = 1;
      i++;
      if(i >= argc) {
	*iret = i; return 1;
      }
      if((c->maskvarname = (char *) malloc(strlen(argv[i]+1))) == NULL)
	error(ERR_MEMALLOC);
      sprintf(c->maskvarname,"%s",argv[i]);
    }
    else
      i--;
  }
  else
    i--;
    

  /* Parse the parameter string to get the number of output parameters,
     this is needed for the CreateOutputColumns function. Later
     CompileAllExpressions will create the variables. */
  k = 0; j = 0; nparam = 0;
  do {
    if(c->paramliststring[k] == '\0' || c->paramliststring[k] == ',') {
      oldval = c->paramliststring[k];
      c->paramliststring[k] = '\0';
      if(!nparam) {
	if((c->paramnames = (char **) malloc(sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
      } else {
	if((c->paramnames = (char **) realloc(c->paramnames, (nparam + 1)*sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
      }
      if((s = strlen(&(c->paramliststring[j]))) == 0)
	error2(ERR_BADVARIABLENAME,"");
      if((c->paramnames[nparam] = (char *) malloc((s+1))) == NULL)
	error(ERR_MEMALLOC);
      sprintf(c->paramnames[nparam],"%s",&(c->paramliststring[j]));
      j = k+1;
      c->paramliststring[k] = oldval;
      nparam++;
    }
    k++;
  } while(c->paramliststring[k-1] != '\0');
  
  if((c->params = (_Variable **) malloc(nparam * sizeof(_Variable *))) == NULL)
    error(ERR_MEMALLOC);
  c->Nparams = nparam;

  *iret = i;
  return 0;
}
