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
#include "vartools_functionpointers.h"

//#ifdef DYNAMICLIB  
//#include <ltdl.h>
//#endif



//void Set_Function_Pointers_Callback(_VARTOOLS_FUNCTION_POINTER_STRUCT *fptr);

_PythonCommand *CreatePythonCommandStruct(ProgramData *p, char *argv0) {
#ifdef _HAVE_PYTHON
  _PythonCommand *c;
  if((c = (_PythonCommand *) malloc(sizeof(_PythonCommand))) == NULL)
    error(ERR_MEMALLOC);

  c->Nvars = 0;
  c->vars = NULL;
  c->isvaroutput = NULL;
  c->Nvars_outonly = 0;
  c->outonlyvars = NULL;
  c->Nlcvars_nonupdate = 0;
  c->outlcvecs_invars = NULL;
  c->outlcvecs_outonlyvars = NULL;
  c->lcvars_nonupdate = NULL;

  c->IsPythonRunning = NULL;
  RegisterScalarData(p, (void *) &(c->IsPythonRunning), VARTOOLS_TYPE_INT, 0);

  c->sockets = NULL;
  RegisterScalarData(p, (void *) &(c->sockets), VARTOOLS_TYPE_INT, 2);

  if((c->progname = (char *) malloc((strlen(argv0)+1)*sizeof(char))) == NULL)
    error(ERR_MEMALLOC);
  sprintf(c->progname, argv0);

  c->pythoninitializationtext = NULL;
  c->len_pythoninitializationtextstring = 0;
  c->pythoncommandstring = NULL;
  c->len_pythoncommandstring = 0;
  c->inputpythonfilename = NULL;
  c->pythonobjects = NULL;
  c->iscontinueprocess = 0;
  c->continueprocesscommandptr = NULL;

  c->Nchildren = 0;
  c->childcommandptrs = NULL;
  c->childcnumvals = NULL;

  c->cid = 0;

  c->inoutvarliststring = NULL;
  c->inoutvarnames = NULL;
  c->Ninoutvarnames = 0;

  c->invarliststring = NULL;
  c->outvarliststring = NULL;

  c->invarnames = NULL;
  c->outvarnames = NULL;

  c->Ninvarnames = 0;
  c->Noutvarnames = 0;

  c->outcolumnliststring = NULL;
  c->outcolumnnames = NULL;
  c->Noutcolumnvars = 0;
  c->outcolumnvars = NULL;

  c->processallvariables = 1;

  c->RequireReadAll = 0;

  c->skipfail = 0;

  c->cnum = -1;

  return c;
#else
  return NULL;
#endif
}

void SetupRunPythonVariables(_PythonCommand *c, ProgramData *p) {
  /* Organize the input/output variables into the c->vars and
     c->outonlyvars vectors; also setup the c->isvaroutput, outlcvecs_invars
     and outlcvecs_outonlyvars, and lcvars_nonupdate vectors */
#ifdef _HAVE_PYTHON
  int i, ii, j, k, jvarout;
  int testval;

  if(c->processallvariables) {
    c->Nvars_outonly = 0;
    jvarout = 0;
    for(i=0; i < p->NDefinedVariables; i++) {
      if(p->DefinedVariables[i]->vectortype == VARTOOLS_VECTORTYPE_OUTCOLUMN) {
	if(p->DefinedVariables[i]->outc->cnum < c->cnum) jvarout++;
      }
      else jvarout++;
    }
    c->Nvars = jvarout;
    if(jvarout > 0) {
      if((c->vars = (_Variable **) malloc(c->Nvars * sizeof(_Variable *))) == NULL ||
	 (c->isvaroutput = (int *) malloc(c->Nvars * sizeof(int))) == NULL)
	error(ERR_MEMALLOC);
      jvarout = 0;
      for(i=0; i < p->NDefinedVariables; i++) {
	if(p->DefinedVariables[i]->vectortype == VARTOOLS_VECTORTYPE_OUTCOLUMN) {
	  if(p->DefinedVariables[i]->outc->cnum >= c->cnum) continue;
	}
	c->vars[jvarout] = p->DefinedVariables[i];
	if(p->DefinedVariables[i]->vectortype == VARTOOLS_VECTORTYPE_CONSTANT ||
	   p->DefinedVariables[i]->vectortype == VARTOOLS_VECTORTYPE_OUTCOLUMN) {
	  c->isvaroutput[jvarout] = 0;
	} else {
	  c->isvaroutput[jvarout] = 1;
	}
	jvarout++;
      }
    }
  } else if(c->Ninoutvarnames > 0) {
    c->Nvars = c->Ninoutvarnames;
    c->Nvars_outonly = 0;
    if((c->vars = (_Variable **) malloc(c->Nvars * sizeof(_Variable *))) == NULL ||
       (c->isvaroutput = (int *) malloc(c->Nvars * sizeof(int))) == NULL)
      error(ERR_MEMALLOC);
    for(i=0; i < c->Nvars; i++) {
      for(j=0; j < p->NDefinedVariables; j++) {
	if(!strcmp(c->inoutvarnames[i],p->DefinedVariables[j]->varname)) {
	  c->vars[i] = p->DefinedVariables[j];
	  if(c->vars[i]->vectortype == VARTOOLS_VECTORTYPE_OUTCOLUMN) {
	    if(c->vars[i]->outc->cnum >= c->cnum) {
	      fprintf(stderr,"Error: attempting to pass the output column variable %s to a -python command before that column is created.\n\n",p->DefinedVariables[j]->varname);
	      exit(1);
	    }
	  }
	  if(c->vars[i]->vectortype == VARTOOLS_VECTORTYPE_CONSTANT ||
	     c->vars[i]->vectortype == VARTOOLS_VECTORTYPE_OUTCOLUMN) {
	    c->isvaroutput[i] = 0;
	  } else {
	    c->isvaroutput[i] = 1;
	  }
	  break;
	}
      }
      if(j >= p->NDefinedVariables) {
	/* Check if this variable is listed in one of the output columns; then it will need to be a VARTOOLS_VECTORTYPE_SCALAR to variable */
	for(k=0; k < c->Noutcolumnvars; k++) {
	  if(!strcmp(c->inoutvarnames[i],c->outcolumnnames[k])) {
	    c->vars[i] = SetupScalarVariable(p, c->inoutvarnames[i], VARTOOLS_TYPE_DOUBLE);
	    break;
	  }
	}
	/* Otherwise we will assume that its a light curve vector */
	if(k >= c->Noutcolumnvars) {
	  c->vars[i] = CreateVariable(p, c->inoutvarnames[i], VARTOOLS_TYPE_DOUBLE, VARTOOLS_VECTORTYPE_LC, NULL);
	  RegisterDataFromLightCurve(p,
				     c->vars[i]->dataptr,
				     VARTOOLS_TYPE_DOUBLE,
				     0, 0, 0, 0, 0, NULL, 
				     c->vars[i],
				     -1, c->inoutvarnames[i]);
	}
	/* Either way, we will allow this to be output */
	c->isvaroutput[i] = 1;
      }
    }
  } else {
    c->Nvars =  c->Ninvarnames;
    c->Nvars_outonly = 0;
    if(c->Nvars > 0) {
      if((c->vars = (_Variable **) malloc(c->Nvars * sizeof(_Variable *))) == NULL ||
	 (c->isvaroutput = (int *) malloc(c->Nvars * sizeof(int))) == NULL)
	error(ERR_MEMALLOC);
    }
    for(i=0; i < c->Ninvarnames; i++) {
      c->isvaroutput[i] = 0;
      for(j=0; j < p->NDefinedVariables; j++) {
	if(!strcmp(c->invarnames[i],p->DefinedVariables[j]->varname)) {
	  c->vars[i] = p->DefinedVariables[j];
	  if(c->vars[i]->vectortype == VARTOOLS_VECTORTYPE_OUTCOLUMN) {
	    if(c->vars[i]->outc->cnum >= c->cnum) {
	      fprintf(stderr,"Error: attempting to pass the output column variable %s to a -python command before that column is created.\n\n",p->DefinedVariables[j]->varname);
	      exit(1);
	    }
	  }
	  break;
	}
      }
      if(j >= p->NDefinedVariables) {
	/* Check if this variable is listed in one of the output columns; then it will need to be a VARTOOLS_VECTORTYPE_SCALAR to variable */
	for(k=0; k < c->Noutcolumnvars; k++) {
	  if(!strcmp(c->invarnames[i],c->outcolumnnames[k])) {
	    c->vars[i] = SetupScalarVariable(p, c->invarnames[i], VARTOOLS_TYPE_DOUBLE);
	    break;
	  }
	}
	/* Otherwise we will assume that its a light curve vector */
	if(k >= c->Noutcolumnvars) {
	  c->vars[i] = CreateVariable(p, c->invarnames[i], VARTOOLS_TYPE_DOUBLE, VARTOOLS_VECTORTYPE_LC, NULL);
	  RegisterDataFromLightCurve(p,
				     c->vars[i]->dataptr,
				     VARTOOLS_TYPE_DOUBLE,
				     0, 0, 0, 0, 0, NULL, 
				     c->vars[i],
				     -1, c->invarnames[i]);
	}
      }
    }

    ii = 0;
    for(i=0; i < c->Noutvarnames; i++) {
      for(j=0; j < c->Nvars; j++) {
	if(!strcmp(c->outvarnames[i],c->invarnames[j])) {
	  c->isvaroutput[j] = 1;
	  break;
	}
      }
      if(j < c->Nvars) continue;
      else c->Nvars_outonly += 1;
    }
    if(c->Nvars_outonly > 0) {
      if((c->outonlyvars = (_Variable **) malloc(c->Nvars_outonly * sizeof(_Variable *))) == NULL)
	error(ERR_MEMALLOC);

      ii = 0;
      for(i=0; i < c->Noutvarnames; i++) {
	for(j=0; j < c->Nvars; j++) {
	  if(!strcmp(c->outvarnames[i],c->invarnames[j])) {
	    break;
	  }
	}
	if(j < c->Nvars) continue;

	for(j=0; j < p->NDefinedVariables; j++) {
	  if(!strcmp(c->outvarnames[i],p->DefinedVariables[j]->varname)) {
	    c->outonlyvars[ii] = p->DefinedVariables[j];
	    if(c->outonlyvars[ii]->vectortype == VARTOOLS_VECTORTYPE_OUTCOLUMN) {
	      fprintf(stderr,"Error: the column variable %s cannot be included in the outvars list in a -python command.\n\n",c->outonlyvars[ii]->varname);
	      exit(1);
	    }
	    if(c->outonlyvars[ii]->vectortype == VARTOOLS_VECTORTYPE_CONSTANT) {
	      fprintf(stderr,"Error: the constant variable %s cannot be included in the outvars list in a -python command.\n\n",c->outonlyvars[ii]->varname);
	      exit(1);
	    }
	    break;
	  }
	}
	if(j >= p->NDefinedVariables) {
	  /* Check if this variable is listed in one of the output columns; then it will need to be a VARTOOLS_VECTORTYPE_SCALAR to variable */
	  for(k=0; k < c->Noutcolumnvars; k++) {
	    if(!strcmp(c->outvarnames[i],c->outcolumnnames[k])) {
	      c->outonlyvars[ii] = SetupScalarVariable(p, c->outvarnames[i], VARTOOLS_TYPE_DOUBLE);
	      break;
	    }
	  }

	  /* Otherwise we will assume that its a light curve vector */
	  if(k >= c->Noutcolumnvars) {
	    c->outonlyvars[ii] = CreateVariable(p, c->outvarnames[i], VARTOOLS_TYPE_DOUBLE, VARTOOLS_VECTORTYPE_LC, NULL);
	    RegisterDataFromLightCurve(p,
				       c->outonlyvars[ii]->dataptr,
				       VARTOOLS_TYPE_DOUBLE,
				       0, 0, 0, 0, 0, NULL, 
				       c->outonlyvars[ii],
				       -1, c->outvarnames[i]);
	  }
	}
	ii++;
      }
    }
  }

  if(c->Nvars > 0)
    RegisterScalarData(p, (void *) &(c->outlcvecs_invars), VARTOOLS_TYPE_INT, c->Nvars);
  if(c->Nvars_outonly > 0)
    RegisterScalarData(p, (void *) &(c->outlcvecs_outonlyvars), VARTOOLS_TYPE_INT, c->Nvars_outonly);

  if(c->Noutcolumnvars > 0) {
    if((c->outcolumnvars = (_Variable **) malloc(c->Noutcolumnvars * sizeof(_Variable *))) == NULL)
      error(ERR_MEMALLOC);
    for(i=0; i < c->Noutcolumnvars; i++) {
      for(j=0; j < p->NDefinedVariables; j++) {
	if(!strcmp(c->outcolumnnames[i],p->DefinedVariables[j]->varname)) {
	  if(p->DefinedVariables[j]->vectortype == VARTOOLS_VECTORTYPE_LC ||
	     p->DefinedVariables[j]->datatype == VARTOOLS_TYPE_STRING ||
	     p->DefinedVariables[j]->datatype == VARTOOLS_TYPE_CHAR) {
	    error2(ERR_BADVECTORTYPEFOROUTPUTCOLUMNVARIABLE,c->outcolumnnames[i]);
	  }
	  c->outcolumnvars[i] = p->DefinedVariables[j];
	  break;
	}
      }
      if(j >= p->NDefinedVariables) {
	error2(ERR_PYTHONOUTPUTUNDEFINEDVARIABLE,c->outcolumnnames[i]);
      }
    }
  }

  /* Find any lc vectors that are not updated by this command */
  c->Nlcvars_nonupdate = 0;
  for(i=0; i < p->NDefinedVariables; i++) {
    if(p->DefinedVariables[i]->vectortype == VARTOOLS_VECTORTYPE_LC) {
      testval = 0;
      for(j=0; j < c->Nvars; j++) {
	if(c->isvaroutput[j]) {
	  if(!strcmp(c->vars[j]->varname,p->DefinedVariables[i]->varname)) {
	    testval = 1;
	    break;
	  }
	}
      }
      if(!testval) {
	for(j=0; j < c->Nvars_outonly; j++) {
	  if(!strcmp(c->outonlyvars[j]->varname,p->DefinedVariables[i]->varname)) {
	    testval = 1;
	    break;
	  }
	}
      }
      if(!testval) c->Nlcvars_nonupdate += 1;
    }
  }
  if(c->Nlcvars_nonupdate > 0) {
    if((c->lcvars_nonupdate = (_Variable **) malloc(c->Nlcvars_nonupdate * sizeof(_Variable *))) == NULL)
      error(ERR_MEMALLOC);
    c->Nlcvars_nonupdate = 0;
    for(i=0; i < p->NDefinedVariables; i++) {
      if(p->DefinedVariables[i]->vectortype == VARTOOLS_VECTORTYPE_LC) {
	testval = 0;
	for(j=0; j < c->Nvars; j++) {
	  if(c->isvaroutput[j]) {
	    if(!strcmp(c->vars[j]->varname,p->DefinedVariables[i]->varname)) {
	      testval = 1;
	      break;
	    }
	  }
	}
	if(!testval) {
	  for(j=0; j < c->Nvars_outonly; j++) {
	    if(!strcmp(c->outonlyvars[j]->varname,p->DefinedVariables[i]->varname)) {
	      testval = 1;
	      break;
	    }
	  }
	}
	if(!testval) {
	  c->lcvars_nonupdate[c->Nlcvars_nonupdate] = p->DefinedVariables[i];
	  c->Nlcvars_nonupdate += 1;
	}
      }
    }
  }
#else
  return;
#endif
}

/*
void LoadVartoolsRunPythonLibrary(ProgramData *p) {
  
#ifdef DYNAMICLIB  
  void *func;

  lt_dlhandle lib;
  char tmpstring[MAXLEN];
  char libbasename[MAXLEN];
  char libname[] = "vartoolsrunpython";
  int i1, i2, i3, i, j;

  lib = lt_dlopenext(libname);

  if(lib == NULL) {
    error2(ERR_OPEN_LIBRARY,libname);
  }

  func = lt_dlsym(lib,"tmpinitialize");

  ((void (*)(void)) func)();

  func = lt_dlsym(lib,"VARTOOLS_Set_Function_Pointers");

  if(func == NULL) {
    error2(ERR_LIBRARY_MISSING_FUNCTION,"VARTOOLS_Set_Function_Pointers");
  }
  
  ((void (*)(void (*)(_VARTOOLS_FUNCTION_POINTER_STRUCT *))) 
   func)(&Set_Function_Pointers_Callback);

  func = lt_dlsym(lib,"ParsePythonCommand");

  if(func == NULL) {
    error2(ERR_LIBRARY_MISSING_FUNCTION,"ParsePythonCommand");
  }
  p->VartoolsPythonLib.ParsePythonCommand_ptr = func;

  func = lt_dlsym(lib,"RunPythonCommand");

  if(func == NULL) {
    error2(ERR_LIBRARY_MISSING_FUNCTION,"RunPythonCommand");
  }
  p->VartoolsPythonLib.RunPythonCommand_ptr = func;

  func = lt_dlsym(lib,"InitPythonCommand");
  
  if(func == NULL) {
    error2(ERR_LIBRARY_MISSING_FUNCTION,"InitPythonCommand");
  }
  p->VartoolsPythonLib.InitPythonCommand_ptr = func;

  p->pythonlibraryloaded = 1;
  return;
#else

  return;
#endif
}*/

/*
int ParsePythonCommand(int *inum, int argc, char **argv, ProgramData *p, 
			_PythonCommand *c, Command *allcommands, int cnum)
{
#ifdef _HAVE_PYTHON
  return ((int (*)(int *, int, char **, ProgramData *, _PythonCommand *, Command *, int))p->VartoolsPythonLib.ParsePythonCommand_ptr)(inum, argc, argv, p, c, allcommands, cnum);
#else
  
  return 1;
#endif
}

void RunPythonCommand(ProgramData *p, int lcindex, int threadindex, _PythonCommand *c)
{
#ifdef _HAVE_PYTHON
  ((void (*)(ProgramData *,int, int, _PythonCommand *))p->VartoolsPythonLib.RunPythonCommand_ptr)(p, lcindex, threadindex, c);
  return;
#else
  return;
#endif
}

void InitPythonCommand(ProgramData *p, _PythonCommand *c, int Nlcs)
{
#ifdef _HAVE_PYTHON
  ((void (*)(ProgramData *, _PythonCommand *, int))p->VartoolsPythonLib.InitPythonCommand_ptr)(p, c, Nlcs);
  return;
#else
  return;
#endif
}
*/
