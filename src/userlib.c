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
/*                                                                           */

/* This file contains functions to dynamically load and execute user-compiled */
/* libraries.                                                                 */

#include "commands.h"
#include "programdata.h"
#include "functions.h"
#include "vartools_functionpointers.h"

#include <ctype.h>
#include <ltdl.h>

/*
#ifdef DYNAMICLIB
#ifndef ISWINDOWS
#include <dlfcn.h>
#include <ctype.h>
#else
#include <windows.h>
#endif
#endif 
*/

typedef	void (*initialize_usertype_funcP)(int, void *, void *);


int ParseParameter(ProgramData *p, Command *c, int *iret, char **argv,
		   int argc, const char *keyword, 
		   int Nvec, _ParseFixSpecFixcolumnStruct *s,
		   _ParseParameter_InitializeStruct *ppstruct);

int ParseConstantParameter(ProgramData *p, Command *c, int *iret, char **argv,
			   int argc, const char *keyword, char datatype,
			   void *dataptr, int Ncolumns);

int ParseFixSpecFixcolumn(ProgramData *p, Command *c, int *iret, char **argv,
			  int argc, int Nvec, _ParseFixSpecFixcolumnStruct *s);

void vRegisterDataVector(ProgramData *p, Command *c, void *dataptr, 
			char datatype, int Ncolumns, int source, 
			int output, char *outname, va_list varlist);

void Set_Function_Pointers_Callback(_VARTOOLS_FUNCTION_POINTER_STRUCT *fptr);

void vAdd_Keyword_To_OutputLC_FitsHeader(ProgramData *p, int lcnum, char *keyname,
					char *comment, int hdutouse, int updateexisting,
					int dtype, va_list argp);

int load_user_library(char *libname, ProgramData *p, int islib, ...) {


#ifdef DYNAMICLIB  
/*
#ifndef ISWINDOWS
  void *lib;
#else
  HINSTANCE lib;
#endif
*/

  void *func;

  char tmpstring[MAXLEN];
  char libbasename[MAXLEN];

  int i1, i2, i3, i, j;

  lt_dlhandle lib;

  va_list varlist;

  if(islib) {
    va_start(varlist, islib);
    lib = ((lt_dlhandle) va_arg(varlist,lt_dlhandle));
    va_end(varlist);
  }


  /*
#ifndef ISWINDOWS
  lib = dlopen(libname, RTLD_LAZY);
#else
  lib = LoadLibrary(libname);
#endif
  */
  if(!islib) {
    lib = lt_dlopenext(libname);

    if(lib == NULL) {
      error2(ERR_OPEN_LIBRARY,libname);
    }  
  }

  /* Create space to store the pointers to the various library routines */

  if(!p->NUserLib) {
    if((p->UserLib = (_UserLib *) malloc(sizeof(_UserLib))) == NULL)
      error(ERR_MEMALLOC);
    p->NUserLib = 1;
  }
  else {
    if((p->UserLib = (_UserLib *) realloc(p->UserLib, (p->NUserLib + 1)*sizeof(_UserLib))) == NULL)
      error(ERR_MEMALLOC);
    p->NUserLib += 1;
  }

  j = p->NUserLib - 1;

  /* Strip off the basename of the library to get the names of the functions
     to be called */

  i1 = 0;
  i2 = 0;
  while(libname[i1] != '\n' && libname[i1] != '\0') {
    if(libname[i1] == '/')
      i2 = i1+1;
    i1++;
  }
  i1--;
  i3 = i1;
  while(i1 > i2) {
    if(libname[i1] == '.')
      i3 = i1 - 1;
    i1--;
  }
  for(i=i2; i <= i3; i++) {
    libbasename[i-i2] = libname[i];
  }
  libbasename[i-i2] = '\0';

  /* A library called "USERLIBS/examplelibrary.so" will be expected to
     contain the functions:
        examplelibrary_ParseCL
        examplelibrary_RunCommand
        examplelibrary_ShowSyntax
        examplelibrary_ShowHelp
	examplelibrary_Initialize

     it may also contain the optional functions:
        examplelibrary_ShowExample
	examplelibrary_CloseCommand

     as well as any other functions needed to execute the command

  */

  sprintf(tmpstring,"%s_ParseCL",libbasename);
  func = lt_dlsym(lib,tmpstring);
  /*
#ifndef ISWINDOWS
  func = dlsym(lib,tmpstring);
#else
  func = GetProcAddress(lib,tmpstring);
#endif
  */
  if(func == NULL) {
    error2(ERR_LIBRARY_MISSING_FUNCTION,tmpstring);
  }
  p->UserLib[j].ParseCL_function_ptr = func;

  sprintf(tmpstring,"%s_RunCommand",libbasename);
  func = lt_dlsym(lib,tmpstring);
  /*
#ifndef ISWINDOWS
  func = dlsym(lib,tmpstring);
#else
  func = GetProcAddress(lib,tmpstring);
#endif
  */
  if(func == NULL) {
    error2(ERR_LIBRARY_MISSING_FUNCTION,tmpstring);
  }
  p->UserLib[j].RunCommand_function_ptr = func;

  sprintf(tmpstring,"%s_ShowSyntax",libbasename);
  func = lt_dlsym(lib,tmpstring);
  /*
#ifndef ISWINDOWS
  func = dlsym(lib,tmpstring);
#else
  func = GetProcAddress(lib,tmpstring);
#endif
  */
  if(func == NULL) {
    error2(ERR_LIBRARY_MISSING_FUNCTION,tmpstring);
  }
  p->UserLib[j].ShowSyntax_function = ((ShowSyntax_function_type)(func));

  sprintf(tmpstring,"%s_ShowHelp",libbasename);
  func = lt_dlsym(lib,tmpstring);
  /*
#ifndef ISWINDOWS
  func = dlsym(lib,tmpstring);
#else
  func = GetProcAddress(lib,tmpstring);
#endif
  */
  if(func == NULL) {
    error2(ERR_LIBRARY_MISSING_FUNCTION,tmpstring);
  }
  p->UserLib[j].ShowHelp_function = ((ShowHelp_function_type)(func));

  sprintf(tmpstring,"%s_Initialize",libbasename);
  func = lt_dlsym(lib,tmpstring);
  /*
#ifndef ISWINDOWS
  func = dlsym(lib,tmpstring);
#else
  func = GetProcAddress(lib,tmpstring);
#endif
  */
  if(func == NULL) {
    error2(ERR_LIBRARY_MISSING_FUNCTION,tmpstring);
  }
  p->UserLib[j].Initialize_function = ((Initialize_function_type)(func));

  sprintf(tmpstring,"%s_ShowExample",libbasename);
  func = lt_dlsym(lib,tmpstring);
  /*
#ifndef ISWINDOWS
  func = dlsym(lib,tmpstring);
#else
  func = GetProcAddress(lib,tmpstring);
#endif
  */
  if(func != NULL) {
    p->UserLib[j].ShowExample_function = ((ShowExample_function_type)(func));
  }
  else {
    p->UserLib[j].ShowExample_function = NULL;
  }

  sprintf(tmpstring,"%s_CloseCommand",libbasename);
  func = lt_dlsym(lib,tmpstring);
  /*
#ifndef ISWINDOWS
  func = dlsym(lib,tmpstring);
#else
  func = GetProcAddress(lib,tmpstring);
#endif
  */
  if(func != NULL) {
    p->UserLib[j].CloseCommand_function_ptr = func;
  }
  else {
    p->UserLib[j].CloseCommand_function_ptr = NULL;
  }

  /* Run the initialization routine which gives the name of the
     command, a flag indicating whether or not the routine requires
     all lcs to be read in, a flag indicating whether or not the
     command expects the light curves to be sorted in time, a flag
     indicating whether or not the command expects the times in each
     light curve to be unique, and the size of the struct holding the
     data for this command */

  p->UserLib[j].Initialize_function(p->UserLib[j].commandname,
				    &(p->UserLib[j].RequireReadAll),
				    &(p->UserLib[j].RequireSortLC),
				    &(p->UserLib[j].RequireDistinctTimes),
				    &(p->UserLib[j].sizeuserdata));
     
  /* Share pointers to vartools library functions with the user library */
  func = lt_dlsym(lib,"VARTOOLS_Set_Function_Pointers");
  /*
#ifndef ISWINDOWS
  func = dlsym(lib,"VARTOOLS_Set_Function_Pointers");
#else
  func = GetProcAddress(lib,"VARTOOLS_Set_Function_Pointers");
#endif
  */
  if(func == NULL) {
    /* If this function is missing, it means the user did not link to the
       libvartools.a library in making their own library. */
    error2(ERR_LIBRARY_MISSING_FUNCTION,"VARTOOLS_Set_Function_Pointers");
  }
  
  ((void (*)(void (*)(_VARTOOLS_FUNCTION_POINTER_STRUCT *))) 
   func)(&Set_Function_Pointers_Callback);
  
  /*
  ((void (*)(int (*)(ProgramData *, Command *, int *, char **, int, int, _ParseFixSpecFixcolumnStruct *), 
	     int (*)(double **, double *, int *, int, double, double (*)(double *, int, int, double *, double *, double *, void *), int *, int, double *, double *, double *, void *), 
	     void (*)(ProgramData *, Command *, void *, int, int, int, int, char *, va_list varlist),
	     void (*)(char *, char *, char *, char *, char *, int),
	     void (*)(int *, int *, double ***, int **, int, double, double),
	     void (*)(double, int, double *, double *, int, double *, double, double, double, double, double, int),
	     void (*)(int, double *, double *, int, double *, double,
		      double, double, double, double),
	     void (*)(double *,double *,int,double,double,double *,double *),
	     void (*)(double *,double *,double *,int,double,double *),
	     void (*)(int, double *, double *, double *),
	     double (*)(int, double *, double *, double *, double),
	     void (*)(int, double *, double *, double *, double, int, int))
    ) func)(&ParseFixSpecFixcolumn,&amoeba,&vRegisterDataVector,
	    &GetOutputFilename, &incrementparameters_foramoeba,
	    &integratemandelagoltransitmodel,
	    &mandelagoltransitmodel,
	    &spline, &splint, &spline_monotonic, &splint_monotonic,
	    &medianfilter);
  */

  return 0;
#else

  return 1;
#endif
}

int load_userfunction_library(char *libname, ProgramData *p) {

#ifdef DYNAMICLIB  
  /*
#ifndef ISWINDOWS
  void *lib;
#else
  HINSTANCE lib;
#endif
  */
  void *func;

  lt_dlhandle lib;
  char tmpstring[MAXLEN];
  char libbasename[MAXLEN];

  int i1, i2, i3, i, j;

  lib = lt_dlopenext(libname);
  /*
#ifndef ISWINDOWS
  lib = dlopen(libname, RTLD_LAZY);
#else
  lib = LoadLibrary(libname);
#endif
  */

  if(lib == NULL) {
    error2(ERR_OPEN_LIBRARY,libname);
  }

  /* Call the initialization routine */
  i1 = 0;
  i2 = 0;
  while(libname[i1] != '\n' && libname[i1] != '\0') {
    if(libname[i1] == '/')
      i2 = i1+1;
    i1++;
  }
  i1--;
  i3 = i1;
  while(i1 > i2) {
    if(libname[i1] == '.')
      i3 = i1 - 1;
    i1--;
  }
  for(i=i2; i <= i3; i++) {
    libbasename[i-i2] = libname[i];
  }
  libbasename[i-i2] = '\0';

  /* Share pointers to vartools library functions with the user library */
  func = lt_dlsym(lib,"VARTOOLS_Set_Function_Pointers");
  /*
#ifndef ISWINDOWS
  func = dlsym(lib,"VARTOOLS_Set_Function_Pointers");
#else
  func = GetProcAddress(lib,"VARTOOLS_Set_Function_Pointers");
#endif
  */
  if(func == NULL) {
    /* If this function is missing, it means the user did not link to the
       libvartools.a library in making their own library. */
    error2(ERR_LIBRARY_MISSING_FUNCTION,"VARTOOLS_Set_Function_Pointers");
  }
  
  ((void (*)(void (*)(_VARTOOLS_FUNCTION_POINTER_STRUCT *))) 
   func)(&Set_Function_Pointers_Callback);

  sprintf(tmpstring,"%s_Initialize",libbasename);
  func = lt_dlsym(lib,tmpstring);
  /*
#ifndef ISWINDOWS
  func = dlsym(lib,tmpstring);
#else
  func = GetProcAddress(lib,tmpstring);
#endif
  */
  if(func == NULL) {
    error2(ERR_LIBRARY_MISSING_FUNCTION,tmpstring);
  }

  ((void (*)(ProgramData *))func)(p);
  return 0;
#else

  return 1;
#endif
}

int vRegisterUserFunction(ProgramData *p, char *funcname, int Nexpr, double (*func)(double *), int ishelp, va_list varlist) {
#ifdef DYNAMICLIB
  int i, j;
  char *instr;
  if(p->NUserFunc == 0) {
    p->NUserFunc = 1;
    if((p->UserFunc = (_UserFunc *) malloc(sizeof(_UserFunc))) == NULL)
      error(ERR_MEMALLOC);
  } else {
    p->NUserFunc += 1;
    if((p->UserFunc = (_UserFunc *) realloc(p->UserFunc, (p->NUserFunc*sizeof(_UserFunc)))) == NULL)
      error(ERR_MEMALLOC);
  }
  sprintf(p->UserFunc[p->NUserFunc - 1].funcname, "%s", funcname);
  p->UserFunc[p->NUserFunc - 1].Nargs = Nexpr;
  p->UserFunc[p->NUserFunc - 1].EvalFunction_ptr = func;
  p->UserFunc[p->NUserFunc - 1].ishelp = ishelp;
  if(ishelp) {
    InitOutTextStruct(&(p->UserFunc[p->NUserFunc - 1].helptext));
    instr = va_arg(varlist, char *);
    printtostring_nowrap(&(p->UserFunc[p->NUserFunc - 1].helptext), instr);
    if(Nexpr > 0) {
      if((p->UserFunc[p->NUserFunc - 1].argnames = (OutText *) malloc(Nexpr * sizeof(OutText))) == NULL ||
	 (p->UserFunc[p->NUserFunc - 1].argsummaries = (OutText *) malloc(Nexpr * sizeof(OutText))) == NULL)
	error(ERR_MEMALLOC);
    }
    for(i = 0; i < Nexpr; i++) {
      InitOutTextStruct(&(p->UserFunc[p->NUserFunc - 1].argnames[i]));
      InitOutTextStruct(&(p->UserFunc[p->NUserFunc - 1].argsummaries[i]));
      instr = va_arg(varlist, char *);
      printtostring_nowrap(&(p->UserFunc[p->NUserFunc - 1].argnames[i]), instr);
      instr = va_arg(varlist, char *);
      printtostring_nowrap(&(p->UserFunc[p->NUserFunc - 1].argsummaries[i]), instr);
    }
  }
  return 0;
#else
  return 1;
#endif
}

int RegisterUserFunction(ProgramData *p, char *funcname, int Nexpr, double (*func)(double *), int ishelp, ...) {
#ifdef DYNAMICLIB
  va_list varlist;
  int retval;
  va_start(varlist, ishelp);
  retval = vRegisterUserFunction(p, funcname, Nexpr, func, ishelp, varlist);
  va_end(varlist);
  return retval;
#else
  return 1;
#endif
}

int CheckIfUserCommandIsCalled(ProgramData *p, Command *c, int cn, char *argv) {
#ifdef DYNAMICLIB
  int i;
  int check = 0;
  /* Cycle through all the libraries which have been loaded so far */
  for(i=0; i < p->NUserLib; i++) {
    if(!strcmp(p->UserLib[i].commandname,argv)) {
      /* We Have a match */
      check = 1;
      increaseNcommands(p,&c);
      c[cn].cnum = CNUM_USERCOMMAND;
      if((c[cn].UserCommand = (_UserCommand *) malloc(sizeof(_UserCommand))) == NULL)
	error(ERR_MEMALLOC);
      /*c[cn].UserCommand->lib = &(p->UserLib[i]);*/
      c[cn].UserCommand->libnum = i;
      if((c[cn].UserCommand->userdata = (void *) malloc(p->UserLib[i].sizeuserdata)) == NULL)
	error(ERR_MEMALLOC);

      /*c[cn].require_sort = c[cn].UserCommand->lib->RequireSortLC;*/
      c[cn].require_sort = p->UserLib[i].RequireSortLC;
      /*c[cn].require_distinct = c[cn].UserCommand->lib->RequireDistinctTimes;*/
      c[cn].require_distinct = p->UserLib[i].RequireDistinctTimes;
      c[cn].UserCommand->Nfix = 0;
      c[cn].UserCommand->Ninlist = 0;
      c[cn].UserCommand->Ninlc = 0;
      c[cn].UserCommand->Ncomputed = 0;
      c[cn].UserCommand->Nprior = 0;
      c[cn].UserCommand->Noutput = 0;
      c[cn].UserCommand->Nptrs = 0;
      c[cn].UserCommand->Nexpr = 0;
      
      c[cn].UserCommand->UserDataPointers = NULL;
      c[cn].UserCommand->FixValues = NULL;
      c[cn].UserCommand->OutputData = NULL;
      c[cn].UserCommand->UserDataExpressions = NULL;

      c[cn].UserCommand->expr_strings = NULL;

      if(p->UserLib[i].RequireReadAll)
	p->readallflag = 1;
      break;
    }
  }
  if(!check) {
    /* Try to load a library with a name given by the command */
    if(strlen(argv) <= 1) return check;
    lt_dlhandle lib = lt_dlopenext(&(argv[1]));
    
    if(lib != NULL) {
      load_user_library(&(argv[1]), p, 1, lib);
      i = p->NUserLib-1;
      if(!strcmp(p->UserLib[i].commandname,argv)) {
	/* We Have a match */
	check = 1;
	increaseNcommands(p,&c);
	c[cn].cnum = CNUM_USERCOMMAND;
	if((c[cn].UserCommand = (_UserCommand *) malloc(sizeof(_UserCommand))) == NULL)
	  error(ERR_MEMALLOC);
	/*c[cn].UserCommand->lib = &(p->UserLib[i]);*/
	c[cn].UserCommand->libnum = i;
	if((c[cn].UserCommand->userdata = (void *) malloc(p->UserLib[i].sizeuserdata)) == NULL)
	  error(ERR_MEMALLOC);
	
	/*c[cn].require_sort = c[cn].UserCommand->lib->RequireSortLC;*/
	c[cn].require_sort = p->UserLib[i].RequireSortLC;
	/*c[cn].require_distinct = c[cn].UserCommand->lib->RequireDistinctTimes;*/
	c[cn].require_distinct = p->UserLib[i].RequireDistinctTimes;
	c[cn].UserCommand->Nfix = 0;
	c[cn].UserCommand->Ninlist = 0;
	c[cn].UserCommand->Ninlc = 0;
	c[cn].UserCommand->Ncomputed = 0;
	c[cn].UserCommand->Nprior = 0;
	c[cn].UserCommand->Noutput = 0;
	c[cn].UserCommand->Nptrs = 0;
	c[cn].UserCommand->Nexpr = 0;
	
	c[cn].UserCommand->UserDataPointers = NULL;
	c[cn].UserCommand->FixValues = NULL;
	c[cn].UserCommand->OutputData = NULL;
	c[cn].UserCommand->UserDataExpressions = NULL;
	
	c[cn].UserCommand->expr_strings = NULL;

	if(p->UserLib[i].RequireReadAll)
	  p->readallflag = 1;
      }
    }
  }


    
  return check;
#else
  return 0;
#endif
}


int CheckIfUserCommandHelpIsCalled(ProgramData *p, int all, char *argv) {
#ifdef DYNAMICLIB
  int i;
  int check = 0;
  /* Cycle through all the libraries which have been loaded so far */
  for(i=0; i < p->NUserLib; i++) {
    if(all == 1 || (!strcmp(p->UserLib[i].commandname,argv))) {
      /* We Have a match */
      check = 1;
      p->UserLib[i].ShowSyntax_function(stderr);
      fprintf(stderr,"\n");
      p->UserLib[i].ShowHelp_function(stderr);
      fprintf(stderr,"\n");
    }
  }
  if(!check) {
    /* Try to load a library with a name given by the command */
    if(strlen(argv) <= 1) return check;
    lt_dlhandle lib = lt_dlopenext(&(argv[1]));
    
    if(lib != NULL) {
      load_user_library(&(argv[1]), p, 1, lib);
      i = p->NUserLib-1;
      if(!strcmp(p->UserLib[i].commandname,argv)) {
	/* We Have a match */
	check = 1;
	p->UserLib[i].ShowSyntax_function(stderr);
	fprintf(stderr,"\n");
	p->UserLib[i].ShowHelp_function(stderr);
	fprintf(stderr,"\n");
      }
    }
  }
  return check;
#else
  return 0;
#endif
}


int CheckIfUserCommandExampleIsCalled(ProgramData *p, char *argv) {
#ifdef DYNAMICLIB
  int i;
  int check = 0;
  /* Cycle through all the libraries which have been loaded so far */
  for(i=0; i < p->NUserLib; i++) {
    if((!strcmp(p->UserLib[i].commandname,argv))) {
      /* We Have a match */
      check = 1;
      if(p->UserLib[i].ShowExample_function != NULL) {
	p->UserLib[i].ShowExample_function(stderr);
	fprintf(stderr,"\n");
      } else {
	fprintf(stderr,"No example is available for the user command: %s\n", argv);
      }
    }
  }
  if(!check) {
    /* Try to load a library with a name given by the command */
    if(strlen(argv) <= 1) return check;
    lt_dlhandle lib = lt_dlopenext(&(argv[1]));
    
    if(lib != NULL) {
      load_user_library(&(argv[1]), p, 1, lib);
      i = p->NUserLib-1;
      if(!strcmp(p->UserLib[i].commandname,argv)) {
	/* We Have a match */
	check = 1;
	if(p->UserLib[i].ShowExample_function != NULL) {
	  p->UserLib[i].ShowExample_function(stderr);
	  fprintf(stderr,"\n");
	} else {
	  fprintf(stderr,"No example is available for the user command: %s\n", argv);
	}
      }
    }
  }
  return check;
#else
  return 0;
#endif
}


/* This is an internal function called from the CreateOutputColumns function
   the user should not need to call this function */
void CreateOutputColumns_UserCommand(ProgramData *p, Command *c, int cnum) {
#ifdef DYNAMICLIB
  int i, j, k;
  _UserDataPointer *d;
  _UserCommand *co;
  _UserLib *lib;
  char tmpstring[MAXLEN];
  char fmt[MAXLEN];
  char firstlet;

  double *dblptr;
  float *fltptr;
  long *longptr;
  int *intptr;
  short *shortptr;
  char *charptr;
  char **stringptr;

  double **dbl2ptr;
  float **flt2ptr;
  long **long2ptr;
  int **int2ptr;
  short **short2ptr;
  char **char2ptr;
  char ***string2ptr;

  void *voidptr;

  co = c[cnum].UserCommand;
  lib = &(p->UserLib[co->libnum]);
  for(i=0; i < co->Noutput; i++)
    {
      d = &(co->OutputData[i]);
      switch(d->datatype){
      case VARTOOLS_TYPE_DOUBLE:
	sprintf(fmt,"%%.17g");
	break;
      case VARTOOLS_TYPE_FLOAT:
	sprintf(fmt,"%%.9g");
	break;
      case VARTOOLS_TYPE_STRING:
	sprintf(fmt,"%%s");
	break;
      case VARTOOLS_TYPE_CHAR:
	sprintf(fmt,"%%c");
	break;
      case VARTOOLS_TYPE_INT:
	sprintf(fmt,"%%d");
	break;
      case VARTOOLS_TYPE_SHORT:
	sprintf(fmt,"%%d");
	break;
      case VARTOOLS_TYPE_LONG:
	sprintf(fmt,"%%d");
	break;
      default:
	error(ERR_BADTYPE);
      }
      if(d->Ncolumns <= 0) {
	/*
	switch(d->datatype){
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = *((double **) d->dataptr);
	  voidptr = (void *) dblptr;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  fltptr = *((float **) d->dataptr);
	  voidptr = (void *) fltptr;
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = *((char ***) d->dataptr);
	  voidptr = (void *) stringptr;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = *((char **) d->dataptr);
	  voidptr = (void *) charptr;
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = *((int **) d->dataptr);
	  voidptr = (void *) intptr;
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = *((long **) d->dataptr);
	  voidptr = (void *) longptr;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = *((short **) d->dataptr);
	  voidptr = (void *) shortptr;
	  break;
	default:
	  error(ERR_BADTYPE);
	  }*/
	k = 0;
	while(lib->commandname[k] != '\0' && lib->commandname[k] == '-')
	  k++;
	if(lib->commandname[k] != '\0') {
	  firstlet = (char) toupper((int) lib->commandname[k]);
	  sprintf(tmpstring,"%c%s_%s_%%d",firstlet,&(lib->commandname[k+1]),d->outname);
	} else {
	  sprintf(tmpstring,"%s_%s_%%d",&(lib->commandname[k]),d->outname);
	}
	addcolumn(p, c, cnum, d->datatype, 0, d->dataptr, fmt, 1, 0, 0, 0, tmpstring, cnum);
      } else {
	/*
	switch(d->datatype){
	case VARTOOLS_TYPE_DOUBLE:
	  dbl2ptr = *((double ***) d->dataptr);
	  voidptr = (void *) dbl2ptr;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  flt2ptr = *((float ***) d->dataptr);
	  voidptr = (void *) flt2ptr;
	  break;
	case VARTOOLS_TYPE_STRING:
	  string2ptr = *((char ****) d->dataptr);
	  voidptr = (void *) string2ptr;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  char2ptr = *((char ***) d->dataptr);
	  voidptr = (void *) char2ptr;
	  break;
	case VARTOOLS_TYPE_INT:
	  int2ptr = *((int ***) d->dataptr);
	  voidptr = (void *) int2ptr;
	  break;
	case VARTOOLS_TYPE_LONG:
	  long2ptr = *((long ***) d->dataptr);
	  voidptr = (void *) long2ptr;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  short2ptr = *((short ***) d->dataptr);
	  voidptr = (void *) short2ptr;
	  break;
	default:
	  error(ERR_BADTYPE);
	  }*/
	k = 0;
	while(lib->commandname[k] != '\0' && lib->commandname[k] == '-')
	  k++;
	if(lib->commandname[k] != '\0') {
	  firstlet = (char) toupper((int) lib->commandname[k]);
	  sprintf(tmpstring,"%c%s_%s_%%d_%%d",firstlet,&(lib->commandname[k+1]),d->outname);
	} else {
	  sprintf(tmpstring,"%s_%s_%%d_%%d",&(lib->commandname[k]),d->outname);
	}
	for(j=0; j < d->Ncolumns; j++) {
	  addcolumn(p, c, cnum, d->datatype, 0, d->dataptr, fmt, 2, 0, 0, 0, j, tmpstring, j+1, cnum);
	}
      }
    }
#endif
  return;
}

void SetLinkedColumns_UserCommand(ProgramData *p, Command *c, int cnum) {
#ifdef DYNAMICLIB
  int i, j;
  _UserDataPointer *d;
  _UserCommand *co;

  co = c->UserCommand;

  for(i=0; i < co->Nptrs; i++) {
    d = &(co->UserDataPointers[i]);
    if(d->source == VARTOOLS_SOURCE_PRIORCOLUMN) {
      if(d->Ncolumns <= 0) {
	increaselinkedcols(p, &(d->linkedcolumn[0]), d->priorname[0], cnum);
      } else {
	for(j=0; j < d->Ncolumns; j++) {
	  increaselinkedcols(p, &(d->linkedcolumn[j]), d->priorname[j], cnum);
	}
      }
    }
  }
#else
  return;
#endif
}

void RunUserCommand(ProgramData *p, Command *c, int lc_list_num, int lc_num) {
#ifdef DYNAMICLIB
  int i, k;
  int Nc;
  double *dblptr;
  float *floatptr;
  int *intptr;
  long *longptr;
  short *shortptr;
  char *charptr;
  char **stringptr, **lstringptr;
  void *voidptr;
  double eval_result;

  _UserDataPointer *d;
  _UserDataPointer *fix;
  _UserCommand *co;
  _UserLib *lib;

  void (*RunCommand_Function)(ProgramData *, void *, int, int);

  co=c->UserCommand;
  lib = &(p->UserLib[co->libnum]);
  /* Fill in values which come from previously executed commands, or which
     are fixed */
  for(i=0; i < co->Nptrs; i++) {
    d = &(co->UserDataPointers[i]);
    if(d->source == VARTOOLS_SOURCE_PRIORCOLUMN ||
       (d->source == VARTOOLS_SOURCE_EXISTINGVARIABLE ? 
	(d->expectedvectortype == VARTOOLS_VECTORTYPE_PERSTARDATA &&
	(*(d->existingvariable))->vectortype == VARTOOLS_VECTORTYPE_OUTCOLUMN) : 0)) {
      Nc = d->Ncolumns;
      if(Nc <= 0) {
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = &((*((double **) d->dataptr))[lc_num]);
	  voidptr = (void *) dblptr;
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = &((*((char ***) d->dataptr))[lc_num]);
	  voidptr = (void *) stringptr;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = &((*((char **) d->dataptr))[lc_num]);
	  voidptr = (void *) charptr;
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = &((*((int **) d->dataptr))[lc_num]);
	  voidptr = (int *) intptr;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = &((*((short **) d->dataptr))[lc_num]);
	  voidptr = (short *) shortptr;
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = &((*((long **) d->dataptr))[lc_num]);
	  voidptr = (long *) longptr;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = &((*((float **) d->dataptr))[lc_num]);
	  voidptr = (float *) floatptr;
	  break;
	default:
	  error(ERR_BADTYPE);
	}
	if(d->source == VARTOOLS_SOURCE_PRIORCOLUMN)
	  getoutcolumnvalue(d->linkedcolumn[0], lc_num, 
			    lc_list_num, d->datatype, voidptr);
	else {
	  getoutcolumnvalue((*(d->existingvariable))->outc, lc_num,
			    lc_list_num, d->datatype, voidptr);
	}
      } else {
	for(k=0; k < Nc; k++) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dblptr = &((*((double ***) d->dataptr))[lc_num][k]);
	    voidptr = (void *) dblptr;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    stringptr = &((*((char ****) d->dataptr))[lc_num][k]);
	    voidptr = (void *) stringptr;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    charptr = &((*((char ***) d->dataptr))[lc_num][k]);
	    voidptr = (void *) charptr;
	    break;
	  case VARTOOLS_TYPE_INT:
	    intptr = &((*((int ***) d->dataptr))[lc_num][k]);
	    voidptr = (int *) intptr;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    shortptr = &((*((short ***) d->dataptr))[lc_num][k]);
	    voidptr = (short *) shortptr;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    longptr = &((*((long ***) d->dataptr))[lc_num][k]);
	    voidptr = (long *) longptr;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    floatptr = &((*((float ***) d->dataptr))[lc_num][k]);
	    voidptr = (float *) floatptr;
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	  getoutcolumnvalue(d->linkedcolumn[k], lc_num, lc_list_num,
			    d->datatype, voidptr);
	}
      }
    } else if(d->source == VARTOOLS_SOURCE_FIXED) {
      fix = &(co->FixValues[d->Nfixptr]);
      Nc = d->Ncolumns;
      if(Nc <= 0) {
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = &((*((double **) d->dataptr))[lc_num]);
	  *dblptr = *((double *) fix->dataptr);
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = &((*((char ***) d->dataptr))[lc_num]);
	  sprintf(*stringptr, "%s", *((char **) fix->dataptr));
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = &((*((char **) d->dataptr))[lc_num]);
	  *charptr = *((char *) fix->dataptr);
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = &((*((int **) d->dataptr))[lc_num]);
	  *intptr = *((int *) fix->dataptr);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = &((*((short **) d->dataptr))[lc_num]);
	  *shortptr = *((short *) fix->dataptr);
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = &((*((long **) d->dataptr))[lc_num]);
	  *longptr = *((long *) fix->dataptr);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = &((*((float **) d->dataptr))[lc_num]);
	  *floatptr = *((float *) fix->dataptr);
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      } else {
	for(k=0; k < Nc; k++) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dblptr = &((*((double ***) d->dataptr))[lc_num][k]);
	    *dblptr = (*((double **) fix->dataptr))[k];
	    break;
	  case VARTOOLS_TYPE_STRING:
	    stringptr = &((*((char ****) d->dataptr))[lc_num][k]);
	    sprintf(*stringptr, "%s", (*((char ***) fix->dataptr))[k]);
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    charptr = &((*((char ***) d->dataptr))[lc_num][k]);
	    *charptr = (*((char **) fix->dataptr))[k];
	    break;
	  case VARTOOLS_TYPE_INT:
	    intptr = &((*((int ***) d->dataptr))[lc_num][k]);
	    *intptr = (*((int **) fix->dataptr))[k];
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    shortptr = &((*((short ***) d->dataptr))[lc_num][k]);
	    *shortptr = (*((short **) fix->dataptr))[k];
	    break;
	  case VARTOOLS_TYPE_LONG:
	    longptr = &((*((long ***) d->dataptr))[lc_num][k]);
	    *longptr = (*((long **) fix->dataptr))[k];
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    floatptr = &((*((float ***) d->dataptr))[lc_num][k]);
	    *floatptr = (*((float **) fix->dataptr))[k];
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      }
    } else if(d->source == VARTOOLS_SOURCE_INLIST) {
      Nc = d->Ncolumns;
      if(Nc <= 0) {
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = &((*((double **) d->dataptr))[lc_num]);
	  *dblptr = ((*((double **) d->inlistdataptr))[lc_list_num]);
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = &((*((char ***) d->dataptr))[lc_num]);
	  lstringptr = &((*((char ***) d->inlistdataptr))[lc_list_num]);
	  sprintf(*stringptr, "%s", *lstringptr);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = &((*((char **) d->dataptr))[lc_num]);
	  *charptr = ((*((char **) d->inlistdataptr))[lc_list_num]);
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = &((*((int **) d->dataptr))[lc_num]);
	  *intptr = ((*((int **) d->inlistdataptr))[lc_list_num]);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = &((*((short **) d->dataptr))[lc_num]);
	  *shortptr = ((*((short **) d->inlistdataptr))[lc_list_num]);
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = &((*((long **) d->dataptr))[lc_num]);
	  *longptr = ((*((long **) d->inlistdataptr))[lc_list_num]);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = &((*((float **) d->dataptr))[lc_num]);
	  *floatptr = ((*((float **) d->inlistdataptr))[lc_list_num]);
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      } else {
	for(k=0; k < Nc; k++) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dblptr = &((*((double ***) d->dataptr))[lc_num][k]);
	    *dblptr = ((*((double ***) d->inlistdataptr))[lc_list_num][k]);
	    break;
	  case VARTOOLS_TYPE_STRING:
	    stringptr = &((*((char ****) d->dataptr))[lc_num][k]);
	    lstringptr = &((*((char ****) d->inlistdataptr))[lc_list_num][k]);
	    sprintf(*stringptr, "%s", *lstringptr);
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    charptr = &((*((char ***) d->dataptr))[lc_num][k]);
	    *charptr = ((*((char ***) d->inlistdataptr))[lc_list_num][k]);
	    break;
	  case VARTOOLS_TYPE_INT:
	    intptr = &((*((int ***) d->dataptr))[lc_num][k]);
	    *intptr = ((*((int ***) d->inlistdataptr))[lc_list_num][k]);
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    shortptr = &((*((short ***) d->dataptr))[lc_num][k]);
	    *shortptr = ((*((short ***) d->inlistdataptr))[lc_list_num][k]);
	    break;
	  case VARTOOLS_TYPE_LONG:
	    longptr = &((*((long ***) d->dataptr))[lc_num][k]);
	    *longptr = ((*((long ***) d->inlistdataptr))[lc_list_num][k]);
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    floatptr = &((*((float ***) d->dataptr))[lc_num][k]);
	    *floatptr = ((*((float ***) d->inlistdataptr))[lc_list_num][k]);
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      }
    } else if(d->source == VARTOOLS_SOURCE_EVALEXPRESSION) {
      eval_result = EvaluateExpression(lc_list_num, lc_num, 0, d->evalexpression);
      switch(d->datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dblptr = &((*((double **) d->dataptr))[lc_num]);
	*dblptr = (double) eval_result;
	break;
      case VARTOOLS_TYPE_STRING:
	stringptr = &((*((char ***) d->dataptr))[lc_num]);
	sprintf(*stringptr, "%f", eval_result);
	break;
      case VARTOOLS_TYPE_CHAR:
	charptr = &((*((char **) d->dataptr))[lc_num]);
	*charptr = (char) eval_result;
	break;
      case VARTOOLS_TYPE_INT:
	intptr = &((*((int **) d->dataptr))[lc_num]);
	*intptr = (int) eval_result;
	break;
      case VARTOOLS_TYPE_SHORT:
	shortptr = &((*((short **) d->dataptr))[lc_num]);
	*shortptr = (short) eval_result;
	break;
      case VARTOOLS_TYPE_LONG:
	longptr = &((*((long **) d->dataptr))[lc_num]);
	*longptr = (long) eval_result;
	break;
      case VARTOOLS_TYPE_FLOAT:
	floatptr = &((*((float **) d->dataptr))[lc_num]);
	*floatptr = (float) eval_result;
	break;
      default:
	error(ERR_BADTYPE);
      }
    } else if(d->source == VARTOOLS_SOURCE_EVALEXPRESSION_LC) {
      for(k=0; k < p->NJD[lc_num]; k++) {
	eval_result = EvaluateExpression(lc_list_num, lc_num, k, d->evalexpression);
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = &((*((double ***) d->dataptr))[lc_num][k]);
	  *dblptr = (double) eval_result;
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = &((*((char ****) d->dataptr))[lc_num][k]);
	  sprintf(*stringptr, "%f", eval_result);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = &((*((char ***) d->dataptr))[lc_num][k]);
	  *charptr = (char) eval_result;
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = &((*((int ***) d->dataptr))[lc_num][k]);
	  *intptr = (int) eval_result;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = &((*((short ***) d->dataptr))[lc_num][k]);
	  *shortptr = (short) eval_result;
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = &((*((long ***) d->dataptr))[lc_num][k]);
	  *longptr = (long) eval_result;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = &((*((float ***) d->dataptr))[lc_num][k]);
	  *floatptr = (float) eval_result;
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      }
    } else if(d->source == VARTOOLS_SOURCE_EXISTINGVARIABLE ? 
	      (d->expectedvectortype == VARTOOLS_VECTORTYPE_PERSTARDATA &&
	       (*(d->existingvariable))->vectortype == VARTOOLS_VECTORTYPE_INLIST) : 0) {
      Nc = d->Ncolumns;
      if(Nc <= 0) {
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = &((*((double **) d->dataptr))[lc_num]);
	  *dblptr = ((*((double **) (*(d->existingvariable))->dataptr))[lc_list_num]);
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = &((*((char ***) d->dataptr))[lc_num]);
	  lstringptr = &((*((char ***) (*(d->existingvariable))->dataptr))[lc_list_num]);
	  sprintf(*stringptr, "%s", *lstringptr);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = &((*((char **) d->dataptr))[lc_num]);
	  *charptr = ((*((char **) (*(d->existingvariable))->dataptr))[lc_list_num]);
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = &((*((int **) d->dataptr))[lc_num]);
	  *intptr = ((*((int **) (*(d->existingvariable))->dataptr))[lc_list_num]);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = &((*((short **) d->dataptr))[lc_num]);
	  *shortptr = ((*((short **) (*(d->existingvariable))->dataptr))[lc_list_num]);
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = &((*((long **) d->dataptr))[lc_num]);
	  *longptr = ((*((long **) (*(d->existingvariable))->dataptr))[lc_list_num]);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = &((*((float **) d->dataptr))[lc_num]);
	  *floatptr = ((*((float **) (*(d->existingvariable))->dataptr))[lc_list_num]);
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      } else {
	error(ERR_CODEERROR);
      }
    }  else if(d->source == VARTOOLS_SOURCE_EXISTINGVARIABLE ? 
	       (d->expectedvectortype == VARTOOLS_VECTORTYPE_PERSTARDATA &&
		(*(d->existingvariable))->vectortype == VARTOOLS_VECTORTYPE_CONSTANT) : 0) {
      Nc = d->Ncolumns;
      if(Nc <= 0) {
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = &((*((double **) d->dataptr))[lc_num]);
	  *dblptr = ((*((double *) (*(d->existingvariable))->dataptr)));
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = &((*((char ***) d->dataptr))[lc_num]);
	  lstringptr = &((*((char **) (*(d->existingvariable))->dataptr)));
	  sprintf(*stringptr, "%s", *lstringptr);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = &((*((char **) d->dataptr))[lc_num]);
	  *charptr = ((*((char *) (*(d->existingvariable))->dataptr)));
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = &((*((int **) d->dataptr))[lc_num]);
	  *intptr = ((*((int *) (*(d->existingvariable))->dataptr)));
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = &((*((short **) d->dataptr))[lc_num]);
	  *shortptr = ((*((short *) (*(d->existingvariable))->dataptr)));
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = &((*((long **) d->dataptr))[lc_num]);
	  *longptr = ((*((long *) (*(d->existingvariable))->dataptr)));
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = &((*((float **) d->dataptr))[lc_num]);
	  *floatptr = ((*((float *) (*(d->existingvariable))->dataptr)));
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      } else {
	error(ERR_CODEERROR);
      }
    } else if(d->source == VARTOOLS_SOURCE_EXISTINGVARIABLE ? 
	      (d->expectedvectortype == VARTOOLS_VECTORTYPE_PERSTARDATA &&
	       (*(d->existingvariable))->vectortype == VARTOOLS_VECTORTYPE_SCALAR) : 0) {
      Nc = d->Ncolumns;
      if(Nc <= 0) {
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = &((*((double **) d->dataptr))[lc_num]);
	  *dblptr = ((*((double **) (*(d->existingvariable))->dataptr))[lc_num]);
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = &((*((char ***) d->dataptr))[lc_num]);
	  lstringptr = &((*((char ***) (*(d->existingvariable))->dataptr))[lc_num]);
	  sprintf(*stringptr, "%s", *lstringptr);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = &((*((char **) d->dataptr))[lc_num]);
	  *charptr = ((*((char **) (*(d->existingvariable))->dataptr))[lc_num]);
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = &((*((int **) d->dataptr))[lc_num]);
	  *intptr = ((*((int **) (*(d->existingvariable))->dataptr))[lc_num]);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = &((*((short **) d->dataptr))[lc_num]);
	  *shortptr = ((*((short **) (*(d->existingvariable))->dataptr))[lc_num]);
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = &((*((long **) d->dataptr))[lc_num]);
	  *longptr = ((*((long **) (*(d->existingvariable))->dataptr))[lc_num]);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = &((*((float **) d->dataptr))[lc_num]);
	  *floatptr = ((*((float **) (*(d->existingvariable))->dataptr))[lc_num]);
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      } else {
	error(ERR_CODEERROR);
      }
    }
  }

  RunCommand_Function = (void (*)(ProgramData *, void *, int, int)) 
    lib->RunCommand_function_ptr;

  RunCommand_Function(p, c->UserCommand->userdata,
		      lc_list_num, lc_num);
#else
  return;
#endif
}

void RunUserCommand_all_lcs(ProgramData *p, Command *c) {
#ifdef DYNAMICLIB
  int i, k, lc_num;
  int Nc;
  double *dblptr;
  float *floatptr;
  int *intptr;
  long *longptr;
  short *shortptr;
  char *charptr;
  char **stringptr, **lstringptr;
  void *voidptr;
  double eval_result;
  
  _UserDataPointer *d;
  _UserDataPointer *fix;
  _UserCommand *co;
  _UserLib *lib;

  void (*RunCommand_Function)(ProgramData *, void *, int, int);

  co=c->UserCommand;
  lib = &(p->UserLib[co->libnum]);

  /* Fill in values which come from previously executed commands, or which
     are fixed */
  for(lc_num=0; lc_num < p->Nlcs; lc_num++) {
    for(i=0; i < co->Nptrs; i++) {
      d = &(co->UserDataPointers[i]);
      if(d->source == VARTOOLS_SOURCE_PRIORCOLUMN || 
	 (d->source == VARTOOLS_SOURCE_EXISTINGVARIABLE ? 
	  (d->expectedvectortype == VARTOOLS_VECTORTYPE_PERSTARDATA &&
	   (*(d->existingvariable))->vectortype == VARTOOLS_VECTORTYPE_OUTCOLUMN) : 0)) {
	Nc = d->Ncolumns;
	if(Nc <= 0) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dblptr = &((*((double **) d->dataptr))[lc_num]);
	    voidptr = (void *) dblptr;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    stringptr = &((*((char ***) d->dataptr))[lc_num]);
	    voidptr = (void *) stringptr;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    charptr = &((*((char **) d->dataptr))[lc_num]);
	    voidptr = (void *) charptr;
	    break;
	  case VARTOOLS_TYPE_INT:
	    intptr = &((*((int **) d->dataptr))[lc_num]);
	    voidptr = (int *) intptr;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    shortptr = &((*((short **) d->dataptr))[lc_num]);
	    voidptr = (short *) shortptr;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    longptr = &((*((long **) d->dataptr))[lc_num]);
	    voidptr = (long *) longptr;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    floatptr = &((*((float **) d->dataptr))[lc_num]);
	    voidptr = (float *) floatptr;
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	  if(d->source == VARTOOLS_SOURCE_PRIORCOLUMN)
	    getoutcolumnvalue(d->linkedcolumn[0], lc_num, 
			      lc_num, d->datatype, voidptr);
	  else {
	    getoutcolumnvalue((*(d->existingvariable))->outc, lc_num,
			      lc_num, d->datatype, voidptr);
	  }
	  
	} else {
	  for(k=0; k < Nc; k++) {
	    switch(d->datatype) {
	    case VARTOOLS_TYPE_DOUBLE:
	      dblptr = &((*((double ***) d->dataptr))[lc_num][k]);
	      voidptr = (void *) dblptr;
	      break;
	    case VARTOOLS_TYPE_STRING:
	      stringptr = &((*((char ****) d->dataptr))[lc_num][k]);
	      voidptr = (void *) stringptr;
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      charptr = &((*((char ***) d->dataptr))[lc_num][k]);
	      voidptr = (void *) charptr;
	      break;
	    case VARTOOLS_TYPE_INT:
	      intptr = &((*((int ***) d->dataptr))[lc_num][k]);
	      voidptr = (int *) intptr;
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      shortptr = &((*((short ***) d->dataptr))[lc_num][k]);
	      voidptr = (short *) shortptr;
	      break;
	    case VARTOOLS_TYPE_LONG:
	      longptr = &((*((long ***) d->dataptr))[lc_num][k]);
	      voidptr = (long *) longptr;
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      floatptr = &((*((float ***) d->dataptr))[lc_num][k]);
	      voidptr = (float *) floatptr;
	      break;
	    default:
	      error(ERR_BADTYPE);
	    }
	    getoutcolumnvalue(d->linkedcolumn[k], lc_num, lc_num,
			      d->datatype, voidptr);
	  }
	}
      } else if(d->source == VARTOOLS_SOURCE_FIXED) {
	fix = &(co->FixValues[d->Nfixptr]);
	Nc = d->Ncolumns;
	if(Nc <= 0) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dblptr = &((*((double **) d->dataptr))[lc_num]);
	    *dblptr = *((double *) fix->dataptr);
	    break;
	  case VARTOOLS_TYPE_STRING:
	    stringptr = &((*((char ***) d->dataptr))[lc_num]);
	    sprintf(*stringptr, "%s", *((char **) fix->dataptr));
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    charptr = &((*((char **) d->dataptr))[lc_num]);
	    *charptr = *((char *) fix->dataptr);
	    break;
	  case VARTOOLS_TYPE_INT:
	    intptr = &((*((int **) d->dataptr))[lc_num]);
	    *intptr = *((int *) fix->dataptr);
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    shortptr = &((*((short **) d->dataptr))[lc_num]);
	    *shortptr = *((short *) fix->dataptr);
	    break;
	  case VARTOOLS_TYPE_LONG:
	    longptr = &((*((long **) d->dataptr))[lc_num]);
	    *longptr = *((long *) fix->dataptr);
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    floatptr = &((*((float **) d->dataptr))[lc_num]);
	    *floatptr = *((float *) fix->dataptr);
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	} else {
	  for(k=0; k < Nc; k++) {
	    switch(d->datatype) {
	    case VARTOOLS_TYPE_DOUBLE:
	      dblptr = &((*((double ***) d->dataptr))[lc_num][k]);
	      *dblptr = (*((double **) fix->dataptr))[k];
	      break;
	    case VARTOOLS_TYPE_STRING:
	      stringptr = &((*((char ****) d->dataptr))[lc_num][k]);
	      sprintf(*stringptr, "%s", (*((char ***) fix->dataptr))[k]);
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      charptr = &((*((char ***) d->dataptr))[lc_num][k]);
	      *charptr = (*((char **) fix->dataptr))[k];
	      break;
	    case VARTOOLS_TYPE_INT:
	      intptr = &((*((int ***) d->dataptr))[lc_num][k]);
	      *intptr = (*((int **) fix->dataptr))[k];
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      shortptr = &((*((short ***) d->dataptr))[lc_num][k]);
	      *shortptr = (*((short **) fix->dataptr))[k];
	      break;
	    case VARTOOLS_TYPE_LONG:
	      longptr = &((*((long ***) d->dataptr))[lc_num][k]);
	      *longptr = (*((long **) fix->dataptr))[k];
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      floatptr = &((*((float ***) d->dataptr))[lc_num][k]);
	      *floatptr = (*((float **) fix->dataptr))[k];
	      break;
	    default:
	      error(ERR_BADTYPE);
	    }
	  }
	}
      } else if(d->source == VARTOOLS_SOURCE_INLIST) {
	Nc = d->Ncolumns;
	if(Nc <= 0) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dblptr = &((*((double **) d->dataptr))[lc_num]);
	    *dblptr = ((*((double **) d->inlistdataptr))[lc_num]);
	    break;
	  case VARTOOLS_TYPE_STRING:
	    stringptr = &((*((char ***) d->dataptr))[lc_num]);
	    lstringptr = &((*((char ***) d->inlistdataptr))[lc_num]);
	    sprintf(*stringptr, "%s", *lstringptr);
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    charptr = &((*((char **) d->dataptr))[lc_num]);
	    *charptr = ((*((char **) d->inlistdataptr))[lc_num]);
	    break;
	  case VARTOOLS_TYPE_INT:
	    intptr = &((*((int **) d->dataptr))[lc_num]);
	    *intptr = ((*((int **) d->inlistdataptr))[lc_num]);
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    shortptr = &((*((short **) d->dataptr))[lc_num]);
	    *shortptr = ((*((short **) d->inlistdataptr))[lc_num]);
	    break;
	  case VARTOOLS_TYPE_LONG:
	    longptr = &((*((long **) d->dataptr))[lc_num]);
	    *longptr = ((*((long **) d->inlistdataptr))[lc_num]);
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    floatptr = &((*((float **) d->dataptr))[lc_num]);
	    *floatptr = ((*((float **) d->inlistdataptr))[lc_num]);
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	} else {
	  for(k=0; k < Nc; k++) {
	    switch(d->datatype) {
	    case VARTOOLS_TYPE_DOUBLE:
	      dblptr = &((*((double ***) d->dataptr))[lc_num][k]);
	      *dblptr = ((*((double ***) d->inlistdataptr))[lc_num][k]);
	      break;
	    case VARTOOLS_TYPE_STRING:
	      stringptr = &((*((char ****) d->dataptr))[lc_num][k]);
	      lstringptr = &((*((char ****) d->inlistdataptr))[lc_num][k]);
	      sprintf(*stringptr, "%s", *lstringptr);
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      charptr = &((*((char ***) d->dataptr))[lc_num][k]);
	      *charptr = ((*((char ***) d->inlistdataptr))[lc_num][k]);
	      break;
	    case VARTOOLS_TYPE_INT:
	      intptr = &((*((int ***) d->dataptr))[lc_num][k]);
	      *intptr = ((*((int ***) d->inlistdataptr))[lc_num][k]);
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      shortptr = &((*((short ***) d->dataptr))[lc_num][k]);
	      *shortptr = ((*((short ***) d->inlistdataptr))[lc_num][k]);
	      break;
	    case VARTOOLS_TYPE_LONG:
	      longptr = &((*((long ***) d->dataptr))[lc_num][k]);
	      *longptr = ((*((long ***) d->inlistdataptr))[lc_num][k]);
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      floatptr = &((*((float ***) d->dataptr))[lc_num][k]);
	      *floatptr = ((*((float ***) d->inlistdataptr))[lc_num][k]);
	      break;
	    default:
	      error(ERR_BADTYPE);
	    }
	  }
	}
      } else if(d->source == VARTOOLS_SOURCE_EVALEXPRESSION) {
	eval_result = EvaluateExpression(lc_num, lc_num, 0, d->evalexpression);
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = &((*((double **) d->dataptr))[lc_num]);
	  *dblptr = (double) eval_result;
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = &((*((char ***) d->dataptr))[lc_num]);
	  sprintf(*stringptr, "%f", eval_result);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = &((*((char **) d->dataptr))[lc_num]);
	  *charptr = (char) eval_result;
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = &((*((int **) d->dataptr))[lc_num]);
	  *intptr = (int) eval_result;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = &((*((short **) d->dataptr))[lc_num]);
	  *shortptr = (short) eval_result;
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = &((*((long **) d->dataptr))[lc_num]);
	  *longptr = (long) eval_result;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = &((*((float **) d->dataptr))[lc_num]);
	  *floatptr = (float) eval_result;
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      } else if(d->source == VARTOOLS_SOURCE_EVALEXPRESSION_LC) {
	for(k=0; k < p->NJD[lc_num]; k++) {
	  eval_result = EvaluateExpression(lc_num, lc_num, k, d->evalexpression);
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dblptr = &((*((double ***) d->dataptr))[lc_num][k]);
	    *dblptr = (double) eval_result;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    stringptr = &((*((char ****) d->dataptr))[lc_num][k]);
	    sprintf(*stringptr, "%f", eval_result);
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    charptr = &((*((char ***) d->dataptr))[lc_num][k]);
	    *charptr = (char) eval_result;
	    break;
	  case VARTOOLS_TYPE_INT:
	    intptr = &((*((int ***) d->dataptr))[lc_num][k]);
	    *intptr = (int) eval_result;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    shortptr = &((*((short ***) d->dataptr))[lc_num][k]);
	    *shortptr = (short) eval_result;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    longptr = &((*((long ***) d->dataptr))[lc_num][k]);
	    *longptr = (long) eval_result;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    floatptr = &((*((float ***) d->dataptr))[lc_num][k]);
	    *floatptr = (float) eval_result;
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      } else if(d->source == VARTOOLS_SOURCE_EXISTINGVARIABLE ? 
		(d->expectedvectortype == VARTOOLS_VECTORTYPE_PERSTARDATA &&
		 (*(d->existingvariable))->vectortype == VARTOOLS_VECTORTYPE_INLIST) : 0) {
	Nc = d->Ncolumns;
	if(Nc <= 0) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dblptr = &((*((double **) d->dataptr))[lc_num]);
	    *dblptr = ((*((double **) (*(d->existingvariable))->dataptr))[lc_num]);
	    break;
	  case VARTOOLS_TYPE_STRING:
	    stringptr = &((*((char ***) d->dataptr))[lc_num]);
	    lstringptr = &((*((char ***) (*(d->existingvariable))->dataptr))[lc_num]);
	    sprintf(*stringptr, "%s", *lstringptr);
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    charptr = &((*((char **) d->dataptr))[lc_num]);
	    *charptr = ((*((char **) (*(d->existingvariable))->dataptr))[lc_num]);
	    break;
	  case VARTOOLS_TYPE_INT:
	    intptr = &((*((int **) d->dataptr))[lc_num]);
	    *intptr = ((*((int **) (*(d->existingvariable))->dataptr))[lc_num]);
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    shortptr = &((*((short **) d->dataptr))[lc_num]);
	    *shortptr = ((*((short **) (*(d->existingvariable))->dataptr))[lc_num]);
	    break;
	  case VARTOOLS_TYPE_LONG:
	    longptr = &((*((long **) d->dataptr))[lc_num]);
	    *longptr = ((*((long **) (*(d->existingvariable))->dataptr))[lc_num]);
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    floatptr = &((*((float **) d->dataptr))[lc_num]);
	    *floatptr = ((*((float **) (*(d->existingvariable))->dataptr))[lc_num]);
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	} else {
	  error(ERR_CODEERROR);
	}
      }  else if(d->source == VARTOOLS_SOURCE_EXISTINGVARIABLE ? 
		 (d->expectedvectortype == VARTOOLS_VECTORTYPE_PERSTARDATA &&
		  (*(d->existingvariable))->vectortype == VARTOOLS_VECTORTYPE_CONSTANT) : 0) {
	Nc = d->Ncolumns;
	if(Nc <= 0) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dblptr = &((*((double **) d->dataptr))[lc_num]);
	    *dblptr = ((*((double *) (*(d->existingvariable))->dataptr)));
	    break;
	  case VARTOOLS_TYPE_STRING:
	    stringptr = &((*((char ***) d->dataptr))[lc_num]);
	    lstringptr = &((*((char **) (*(d->existingvariable))->dataptr)));
	    sprintf(*stringptr, "%s", *lstringptr);
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    charptr = &((*((char **) d->dataptr))[lc_num]);
	    *charptr = ((*((char *) (*(d->existingvariable))->dataptr)));
	    break;
	  case VARTOOLS_TYPE_INT:
	    intptr = &((*((int **) d->dataptr))[lc_num]);
	    *intptr = ((*((int *) (*(d->existingvariable))->dataptr)));
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    shortptr = &((*((short **) d->dataptr))[lc_num]);
	    *shortptr = ((*((short *) (*(d->existingvariable))->dataptr)));
	    break;
	  case VARTOOLS_TYPE_LONG:
	    longptr = &((*((long **) d->dataptr))[lc_num]);
	    *longptr = ((*((long *) (*(d->existingvariable))->dataptr)));
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    floatptr = &((*((float **) d->dataptr))[lc_num]);
	    *floatptr = ((*((float *) (*(d->existingvariable))->dataptr)));
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	} else {
	  error(ERR_CODEERROR);
	}
      } else if(d->source == VARTOOLS_SOURCE_EXISTINGVARIABLE ? 
		(d->expectedvectortype == VARTOOLS_VECTORTYPE_PERSTARDATA &&
		 (*(d->existingvariable))->vectortype == VARTOOLS_VECTORTYPE_SCALAR) : 0) {
	Nc = d->Ncolumns;
	if(Nc <= 0) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dblptr = &((*((double **) d->dataptr))[lc_num]);
	    *dblptr = ((*((double **) (*(d->existingvariable))->dataptr))[lc_num]);
	    break;
	  case VARTOOLS_TYPE_STRING:
	    stringptr = &((*((char ***) d->dataptr))[lc_num]);
	    lstringptr = &((*((char ***) (*(d->existingvariable))->dataptr))[lc_num]);
	    sprintf(*stringptr, "%s", *lstringptr);
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    charptr = &((*((char **) d->dataptr))[lc_num]);
	    *charptr = ((*((char **) (*(d->existingvariable))->dataptr))[lc_num]);
	    break;
	  case VARTOOLS_TYPE_INT:
	    intptr = &((*((int **) d->dataptr))[lc_num]);
	    *intptr = ((*((int **) (*(d->existingvariable))->dataptr))[lc_num]);
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    shortptr = &((*((short **) d->dataptr))[lc_num]);
	    *shortptr = ((*((short **) (*(d->existingvariable))->dataptr))[lc_num]);
	    break;
	  case VARTOOLS_TYPE_LONG:
	    longptr = &((*((long **) d->dataptr))[lc_num]);
	    *longptr = ((*((long **) (*(d->existingvariable))->dataptr))[lc_num]);
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    floatptr = &((*((float **) d->dataptr))[lc_num]);
	    *floatptr = ((*((float **) (*(d->existingvariable))->dataptr))[lc_num]);
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	} else {
	  error(ERR_CODEERROR);
	}
      }


    }
  }

  RunCommand_Function = (void (*)(ProgramData *, void *, int, int)) 
    lib->RunCommand_function_ptr;

  RunCommand_Function(p, c->UserCommand->userdata,
		      -1, -1);
#else
  return;
#endif
}


void CloseUserCommand(ProgramData *p, Command *c) {
#ifdef DYNAMICLIB

  _UserCommand *co;
  _UserLib *lib;

  void (*CloseCommand_Function)(ProgramData *, void *);

  co=c->UserCommand;
  lib = &(p->UserLib[co->libnum]);

  if(lib->CloseCommand_function_ptr != NULL) {
    CloseCommand_Function = (void (*)(ProgramData *, void *))
      lib->CloseCommand_function_ptr;
    
    CloseCommand_Function(p, c->UserCommand->userdata);
  }

#else
  return;
#endif
}

int ParseCL_UserCommand(ProgramData *p, Command *c, int *iterm, char **argv, int argc) {
#ifdef DYNAMICLIB
  int retval;
  int itermnew;

  int idiff;

  _UserLib *lib;

  int (*ParseCL_function)(ProgramData *, Command *, void *, int *, char **, int);

  lib = &(p->UserLib[c->UserCommand->libnum]);

  itermnew = 0;

  ParseCL_function = (int (*)(ProgramData *, Command *, void *, int *, char **, int))
    lib->ParseCL_function_ptr;

  idiff = argc - (*iterm);

  retval = (int) ParseCL_function(p, c,
			    c->UserCommand->userdata,
			    &itermnew,
			    &(argv[*iterm]),
			    idiff);
  if(!retval) {
    *iterm = (*iterm) + (itermnew - 1);
    return 0;
  }
  else {
    lib->ShowSyntax_function(stderr);
    exit(ERR_USAGE);
  }
#else
  return 0;
#endif
}

void GetColumnNameForRecentCommand(ProgramData *p, Command *c, int cnum,
				   const char *paramname, char *outcolname)
/* Use this function to find the column name of the last command of type
   cnum, and the parameter paramname. The result is returned in outcolname.
   Options for cnum and paramname are:
        CNUM_BLS  -   Period
                      Tc
                      Q
        CNUM_AOV  -   Period
        CNUM_LS   -   Period
        CNUM_HARMAOV - Period
*/
{
  int k, i;
  switch(cnum) {
  case CNUM_BLS:
    k = p->Ncommands - 1;
    i = -1;
    while(k > 0) {
      if(c[i].cnum == CNUM_BLS)
	break;
      i--;
      k--;
    }
    if(k <= 0) {
      error2(ERR_NO_PREVIOUS_COMMAND,"BLS");
    }
    if(!strcmp(paramname,"Period")) {
      sprintf(outcolname,"BLS_Period_1_%d",(k-1));
      return;
    }
    else if(!strcmp(paramname,"Tc")) {
      sprintf(outcolname,"BLS_Tc_1_%d",(k-1));
      return;
    }
    else if(!strcmp(paramname,"Q")) {
      sprintf(outcolname,"BLS_Qtran_1_%d",(k-1));
      return;
    }
    else {
      error(ERR_INVALID_CALL_GETCOLUMNNAMEFORRECENTCOMMAND);
    }
    break;
  case CNUM_AOV:
    k = p->Ncommands - 1;
    i = -1;
    while(k > 0) {
      if(c[i].cnum == CNUM_AOV)
	break;
      i--;
      k--;
    }
    if(k <= 0) {
      error2(ERR_NO_PREVIOUS_COMMAND,"aov");
    }
    if(!strcmp(paramname,"Period")) {
      sprintf(outcolname,"Period_1_%d",(k-1));
      return;
    }
    else {
      error(ERR_INVALID_CALL_GETCOLUMNNAMEFORRECENTCOMMAND);
    }
    break;
  case CNUM_HARMAOV:
    k = p->Ncommands - 1;
    i = -1;
    while(k > 0) {
      if(c[i].cnum == CNUM_HARMAOV)
	break;
      i--;
      k--;
    }
    if(k <= 0) {
      error2(ERR_NO_PREVIOUS_COMMAND,"aov_harm");
    }
    if(!strcmp(paramname,"Period")) {
      sprintf(outcolname,"Period_1_%d",(k-1));
      return;
    }
    else {
      error(ERR_INVALID_CALL_GETCOLUMNNAMEFORRECENTCOMMAND);
    }
    break;
  case CNUM_LS:
    k = p->Ncommands - 1;
    i = -1;
    while(k > 0) {
      if(c[i].cnum == CNUM_LS)
	break;
      i--;
      k--;
    }
    if(k <= 0) {
      error2(ERR_NO_PREVIOUS_COMMAND,"LS");
    }
    if(!strcmp(paramname,"Period")) {
      sprintf(outcolname,"LS_Period_1_%d",(k-1));
      return;
    }
    else {
      error(ERR_INVALID_CALL_GETCOLUMNNAMEFORRECENTCOMMAND);
    }
    break;
  default:
    error(ERR_INVALID_CALL_GETCOLUMNNAMEFORRECENTCOMMAND);
  }
}

void GetDoubleParameterValue(int threadid, int lcid,
			     double *outparam, int source,
			     double fixvalue, double *inlistvec,
			     OutColumn *column, _Expression *exprsn)
/* Fills *outparam with the correct value depending on the value of source.
   source = VARTOOLS_SOURCE_FIXED: *outparam will be set to fixvalue
   source = VARTOOLS_SOURCE_INLIST: *outparam will be set to inlistvec[lcid]
   source = VARTOOLS_SOURCE_PRIORCOLUMN: *outparam will be set based on the 
                             linked column.
   source = VARTOOLS_SOURCE_PRIORCOLUMN: *outparam will be set by evaluating
                             the expression. 

   If any of the pointers are NULL, that option will not be supported and *outparam will be unchanged.
*/
{
  switch(source){
  case VARTOOLS_SOURCE_FIXED:
    *outparam = fixvalue;
    break;
  case VARTOOLS_SOURCE_INLIST:
    if(inlistvec != NULL)
      *outparam = inlistvec[lcid];
    break;
  case VARTOOLS_SOURCE_PRIORCOLUMN:
    if(column != NULL)
      getoutcolumnvalue(column, threadid, lcid, VARTOOLS_TYPE_DOUBLE, outparam);
    break;
  case VARTOOLS_SOURCE_EVALEXPRESSION:
    if(exprsn != NULL)
      *outparam = EvaluateExpression(lcid, threadid, 0, exprsn);
    break;
  }
  return;
}

void RegisterDataVector(ProgramData *p, Command *c, void *dataptr, 
			char datatype, int Ncolumns, int source, 
			int output, char *outname, ...) 
/* Use this function to register data for a user command. The data is
   a vector, or 2-d array, for which each row corresponds to a
   different light curve, or to a different observation within a light curve. 
   In general, data which holds input or output parameters which are
   light curve dependent, or which are read-in from the light curves,
   should be registered with this function.

   The program will take care of allocating memory for registered
   vectors, and will align the vector with its appropriate data-source
   if necessary (input list, a previously executed command, or a light curve)

   p = pointer to the main ProgramData structure.

   c = pointer to the command structure for the given user command.

   dataptr = pointer to array storing the data. For example, if the
             user wishes to register a vector of doubles, they should
             define it as a pointer (e.g. double* mydoubles) and then
             pass a pointer to it to this function (&mydoubles). This can
             be NULL if source = VARTOOLS_SOURCE_EXISTINGVARIABLE

   datatype = VARTOOLS_TYPE_DOUBLE, VARTOOLS_TYPE_INT, VARTOOLS_TYPE_SHORT,
              VARTOOLS_TYPE_FLOAT, VARTOOLS_TYPE_LONG, VARTOOLS_TYPE_CHAR,
              VARTOOLS_TYPE_STRING, VARTOOLS_TYPE_USERDEF.

   Ncolumns = Number of columns in the array. If this is 0 and the
              source is not VARTOOLS_SOURCE_LC,
              VARTOOLS_SOURCE_EVALEXPRESSION_LC, or 
              VARTOOLS_SOURCE_EXISTINGVARIABLE, then dataptr is a
              pointer to a vector (e.g. double **), if it is > 0 then
              it is a pointer to an array (e.g. double ***). If the
              source is VARTOOLS_SOURCE_LC, or
              VARTOOLS_SOURCE_EVALEXPRESSION_LC, then if this is 0 it
              is a pointer to an array (e.g. it should be type double
              ***), if it is > 0 or < 0 it is a pointer to a 3-d array
              (e.g. type double ****). For
              VARTOOLS_SOURCE_EVALEXPRESSION_LC, and 
              VARTOOLS_SOURCE_EXISTINGVARIABLE only Nc = 0 is
              permitted. For VARTOOLS_SOURCE_LC and
              VARTOOLS_SOURCE_EVALEXPRESSION_LC, Nc=0 implies that the
              array will be indexed as (*dataptr)[Nthread][NJD]. Nc >
              0 means it is indexed as
              (*dataptr)[Nthread][Nc][NJD]. Nc < 0 means it is indexed
              as (*dataptr)[Nthread][NJD][Nc].

   source = VARTOOLS_SOURCE_INLIST -- data to be read in from an input list.
            VARTOOLS_SOURCE_COMPUTED -- data that will be computed by the
                                        command.
            VARTOOLS_SOURCE_FIXED -- data that will be fixed for all light
                                    curves (typically read from the command
                                    line).
            VARTOOLS_SOURCE_PRIORCOLUMN -- data to be taken from the output
                                          of a previous command.
            VARTOOLS_SOURCE_LC -- data to be read in as a column from each
                                  light curve.
            VARTOOLS_SOURCE_RECENTCOMMAND -- similar to 
                                  VARTOOLS_SOURCE_PRIORCOLUMN, in this case
                                  the data is to be taken from the most recent
                                  instance of a command. This will determine
                                  the appropriate column name for the desired
                                  parameter. Ncolumns must be <= 0.
            VARTOOLS_SOURCE_EVALEXPRESSION -- data that will be determined by
                                  evaluating an expression.
            VARTOOLS_SOURCE_EVALEXPRESSION_LC -- a light curve vector
                                  whose data will be determined by
                                  evaluating an expression.
            VARTOOLS_SOURCE_EXISTINGVARIABLE -- in this case the data will come
                                  from an existing variable that should have
                                  already been defined by an earlier vartools
                                  command. 

   output = 1 if the data should be included in the output ascii
            table, 0 if it will not be included. Note, if the source
            is VARTOOLS_SOURCE_LC or
            VARTOOLS_SOURCE_EVALEXPRESSION_LC, then the data will not
            be included in the ascii table, no matter what output is
            set equal to. If output = 1 and dataptr = NULL, and error
            will be given. Also if output = 1 and source =
            VARTOOLS_SOURCE_EXISTINGVARIABLE, an error will be given.


   outname = root name of the corresponding column in the output ascii
             table (can be NULL if output == 0). The name of the
             library will be added as prefix to the name, and the
             column number (if Ncolumns >= 1) and command numbers will
             be appended to the name.

   ... --- additional variables depend on the source and possibly the type:

   VARTOOLS_TYPE_USERDEF:

        In this case the data type is some user-defined structure.
	The following additional parameters are then given in the
	function call. Note that this datatype can only
        be used for VARTOOLS_SOURCE_COMPUTED.

        size_t sizeusertype - Size of the user-defined data type.

	void (*initialize_usertype)(int lc_num, 
	                    void *userdata, 
			    void *userdata2) - 
                             
                             Optional function to call
                             after allocating memory for this data.
                             lc_num will be the index as in
                               the *_RunCommand functions, userdata is a pointer
                               to the data structure in question, and 
                               userdata2 is an optional pointer to additional
                               data used in the initialization procedure.

        void *userdata2 - An optional pointer to additional data that will be
                          passed to the initialize_usertype function. 


   VARTOOLS_SOURCE_INLIST:
        const char *name - string giving the name to associate with this 
                           input column (used with -inputlistformat option).

        int colnum - Column number in the input list from which to read this
                     data. (== 0, the program will take the next column in
                     the list, < 0, the data will not be read-in from the list
		     but memory will be allocated for it (this has the same
		     effect as VARTOOLS_SOURCE_COMPUTED)).

   VARTOOLS_SOURCE_COMPUTED: No additional variables.

   VARTOOLS_SOURCE_FIXED:

        void *fixptr - pointer to a variable storing the fixed value
                       to assign. If Ncolumns >= 1, then fixptr should
                       be a pointer to a vector with Ncolumns rows.

   VARTOOLS_SOURCE_PRIORCOLUMN: 

        char *colname - If Ncolumns <= 0: string giving the name or
                        number of the prior column to associate with
                        this data vector.

        char **colnames - If Ncolumns >= 1: array of strings giving the name
                          or number of the prior columns to associate with
                          each column in this data array.

   VARTOOLS_SOURCE_LC:

        const char *name - string giving the base name to associate with this
                           input column (use with -showinputlcformat option,
                           which is in development).

        int colnum - Column number in the input light curve from which
                     to read this data. (== 0, the program will take
                     the next column in the light curve, it is not
                     suggested to use this option as the number of the
                     next available column can not easily be
                     determined by the user; < 0, the program will not
                     read-in the data from the light curve, but it
                     will allocate memory to store the data; If Nc >
                     0, then this specified the column number for the
                     first column in the data array, subsequent
                     columns in the array are assumed to come from
                     sequential columns in the light curve).

	char *scanformat - A scanf format string to use in reading in
                     the light curve. If it is NULL, then the full
                     input column string will be converted directly to
                     the appropriate data type. If the colnum < 0
                     (data not read-in from the light curve) and the
                     format is not NULL, then it will be treated as an
                     analytic expression which the vector is to be
                     initialized to. Use this, for example, to
                     initialize a vector to the contents of an
                     existing variable rather than reading it in from
                     the light curve directly. An error will be given
                     if Ncolumns != 0, colnum < 0 and this is not
                     NULL.

	char *varnameout - For Ncolumns=0, Optionally associate this
                     data with a (possibly new) variable. This creates
                     a variable that the user can use to access this
                     data from subsequent commands.

   VARTOOLS_SOURCE_RECENTCOMMAND: 

        Note - Ncolumns must be <= 0 for this option.

        int cnum - the command-type to take data from:
                   Allowed values are:
                   CNUM_BLS -- A -BLS command
                   CNUM_AOV -- A -aov command
                   CNUM_LS  -- A -LS command
                   CNUM_HARMAOV -- A -aov_harm command
        
        char *paramname - the type of parameter to take the data from.
                          Allowed values depend on cnum. For any cnum
                          you can provide "Period", for CNUM_BLS this
                          can also be "Tc" or "Q".

    VARTOOLS_SOURCE_EVALEXPRESSION:

        Note - Ncolumns must be == 0 for this option.

        char *exprstring - a string with the expression to be parsed.

    VARTOOLS_SOURCE_EVALEXPRESSION_LC:

        Note - Ncolumns must be == 0 for this option.

        char *exprstring - a string with the expression to be parsed.


	char *varnameout - Optionally associate this
                     data with a (possibly new) variable. This creates
                     a variable that
 the user can use to access this
                     data from subsequent commands.

    VARTOOLS_SOURCE_EXISTINGVARIABLE:

        Note - Ncolumns must be == 0 for this option. Use this option if the
               variable has 1 value per light curve.

        char *varname - name of the existing variable to match to.

	char vectortype - the vectortype that is expected for this variable
                          VARTOOLS_VECTORTYPE_CONSTANT
                          VARTOOLS_VECTORTYPE_SCALAR
                          VARTOOLS_VECTORTYPE_INLIST
                          VARTOOLS_VECTORTYPE_LC
                          VARTOOLS_VECTORTYPE_OUTCOLUMN
			  VARTOOLS_VECTORTYPE_PERSTARDATA - use this to allow
                                constant, scalar, inlist or outcolumn vector
                                types for the source.

        _Variable **var - a pointer to a (_Variable *) that will point to the
                          matched variable.

*/
{
  va_list varlist;
  va_start(varlist, outname);
  vRegisterDataVector(p, c, dataptr, datatype, Ncolumns, source,
		      output, outname, varlist);
  va_end(varlist);
}

void vRegisterDataVector(ProgramData *p, Command *c, void *dataptr, 
			char datatype, int Ncolumns, int source, 
			int output, char *outname, va_list varlist) 
/* This function does the actual work of the RegisterDataVector function
   described below */
{
#ifdef DYNAMICLIB
  size_t sizeval;
  _UserCommand *co;
  char *tmpstring, **tmpstring2, *scanformat, *varnameout;
  int colnum, k, j;
  _UserDataPointer *ptr;
  
  double **dblptr, ***dbl2ptr, ****dbl3ptr;
  float **fltptr, ***flt2ptr, ****flt3ptr;
  int **intptr, ***int2ptr, ****int3ptr;
  long **longptr, ***long2ptr, ****long3ptr;
  short **shortptr, ***short2ptr, ****short3ptr;
  char **charptr, ***char2ptr, ****char3ptr;
  char ***stringptr, ****string2ptr, *****string3ptr;
  
  double *dblptr1;
  float *fltptr1;
  int *intptr1;
  long *longptr1;
  short *shortptr1;
  char *charptr1;
  char **stringptr1;

  void *voidptr;

  initialize_usertype_funcP tmpfptr;
  void *userdata2;

  char *exprstring;

  co = c->UserCommand;

  ptr = co->UserDataPointers;

  /* Make sure that VARTOOLS_TYPE_USERDEF is only
     given for VARTOOLS_SOURCE_COMPUTED, and that if it
     is given, then it is a vector, and it won't be included
     in the output ascii table.
  */
  if(datatype == VARTOOLS_TYPE_USERDEF &&
     (source != VARTOOLS_SOURCE_COMPUTED ||
      Ncolumns != 0 ||
      output != 0)
     )
    error(ERR_BADTYPE);

  if((source == VARTOOLS_SOURCE_EVALEXPRESSION && Ncolumns != 0) ||
     (source == VARTOOLS_SOURCE_EVALEXPRESSION_LC && Ncolumns != 0) ||
     (source == VARTOOLS_SOURCE_EXISTINGVARIABLE && Ncolumns != 0))
    error(ERR_BADTYPE);

  if(dataptr == NULL && output) {
    error(ERR_BADTYPE);
  }

  if(source == VARTOOLS_SOURCE_EXISTINGVARIABLE && output) {
    error(ERR_BADTYPE);
  }

  switch(datatype)
    {
    case VARTOOLS_TYPE_DOUBLE:
      sizeval = sizeof(double);
      break;
    case VARTOOLS_TYPE_INT:
      sizeval = sizeof(int);
      break;
    case VARTOOLS_TYPE_SHORT:
      sizeval = sizeof(short);
      break;
    case VARTOOLS_TYPE_FLOAT:
      sizeval = sizeof(float);
      break;
    case VARTOOLS_TYPE_STRING:
      sizeval = MAXLEN*sizeof(char);
      break;
    case VARTOOLS_TYPE_CHAR:
      sizeval = sizeof(char);
      break;
    case VARTOOLS_TYPE_LONG:
      sizeval = sizeof(long);
      break;
    case VARTOOLS_TYPE_USERDEF:
      sizeval = va_arg(varlist,size_t);
      tmpfptr = va_arg(varlist, initialize_usertype_funcP);
      userdata2 = va_arg(varlist, void *);
      break;
    default:
      break;
    }

  if(!co->Nptrs) {
    if((co->UserDataPointers = (_UserDataPointer *) malloc(sizeof(_UserDataPointer))) == NULL)
      error(ERR_MEMALLOC);
    ptr = co->UserDataPointers;
  } else {
    if((co->UserDataPointers = (_UserDataPointer *) realloc(co->UserDataPointers, (co->Nptrs + 1)*sizeof(_UserDataPointer))) == NULL)
      error(ERR_MEMALLOC);
    ptr = co->UserDataPointers;
    k = 0;
    for(j=0; j < co->Nptrs; j++) {
      if(ptr[j].source == VARTOOLS_SOURCE_EVALEXPRESSION ||
	 ptr[j].source == VARTOOLS_SOURCE_EVALEXPRESSION_LC) {
	co->UserDataExpressions[k] = &(ptr[j].evalexpression);
	k++;
      }
    }
  }

  ptr[co->Nptrs].expectedvectortype = -1;

  ptr[co->Nptrs].datatype = datatype;
  ptr[co->Nptrs].source = source;
  if(source == VARTOOLS_SOURCE_RECENTCOMMAND) {
    ptr[co->Nptrs].source = VARTOOLS_SOURCE_PRIORCOLUMN;
  }
  ptr[co->Nptrs].size_element = sizeval;
  ptr[co->Nptrs].Ncolumns = Ncolumns;
  ptr[co->Nptrs].dataptr = dataptr;
  ptr[co->Nptrs].inlistdataptr = NULL;
  ptr[co->Nptrs].inlcdataptr = NULL;
  if(datatype == VARTOOLS_TYPE_USERDEF) {
    ptr[co->Nptrs].initialize_usertype_ptr = tmpfptr;
    ptr[co->Nptrs].extra_user_data = userdata2;
  } else {
    ptr[co->Nptrs].initialize_usertype_ptr = NULL;
    ptr[co->Nptrs].extra_user_data = NULL;
  }

  if(output) {
    if(!co->Noutput) {
      if((co->OutputData = (_UserDataPointer *) malloc(sizeof(_UserDataPointer))) == NULL)
	error(ERR_MEMALLOC);
    } else {
      if((co->OutputData = (_UserDataPointer *) realloc(co->OutputData, (co->Noutput + 1)*sizeof(_UserDataPointer))) == NULL)
	error(ERR_MEMALLOC);
    }
    
    if((ptr[co->Nptrs].outname = (char *) malloc(MAXLEN)) == NULL)
      error(ERR_MEMALLOC);
    if(outname != NULL)
      sprintf(ptr[co->Nptrs].outname,"%s",outname);
    else
      ptr[co->Nptrs].outname[0] = '\0';

    memcpy(&(co->OutputData[co->Noutput]), &(ptr[co->Nptrs]), sizeof(_UserDataPointer));
    /*co->OutputData[co->Noutput] = &(ptr[co->Nptrs]);*/
    co->Noutput += 1;
  }

  switch(source)
    {
    case VARTOOLS_SOURCE_INLIST:
      /* This vector stores data which is read in from the input list */
      tmpstring = va_arg(varlist,char *);
      colnum = va_arg(varlist,int);
      if(Ncolumns <= 0) {
	switch(datatype){
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = (double **) malloc(sizeof(double *));
	  voidptr = (void *) dblptr;
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = (int **) malloc(sizeof(int *));
	  voidptr = (void *) intptr;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = (short **) malloc(sizeof(short *));
	  voidptr = (void *) shortptr;
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = (long **) malloc(sizeof(long *));
	  voidptr = (void *) longptr;
	case VARTOOLS_TYPE_FLOAT:
	  fltptr = (float **) malloc(sizeof(float *));
	  voidptr = (void *) fltptr;
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = (char ***) malloc(sizeof(char **));
	  voidptr = (void *) stringptr;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = (char **) malloc(sizeof(char *));
	  voidptr = (void *) charptr;
	  break;
	case VARTOOLS_TYPE_USERDEF:
	  voidptr = (void *) malloc(sizeof(void *));
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      }
      else {
	switch(datatype){
	case VARTOOLS_TYPE_DOUBLE:
	  dbl2ptr = (double ***) malloc(sizeof(double **));
	  voidptr = (void *) dbl2ptr;
	  break;
	case VARTOOLS_TYPE_INT:
	  int2ptr = (int ***) malloc(sizeof(int **));
	  voidptr = (void *) int2ptr;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  short2ptr = (short ***) malloc(sizeof(short **));
	  voidptr = (void *) short2ptr;
	  break;
	case VARTOOLS_TYPE_LONG:
	  long2ptr = (long ***) malloc(sizeof(long **));
	  voidptr = (void *) long2ptr;
	case VARTOOLS_TYPE_FLOAT:
	  flt2ptr = (float ***) malloc(sizeof(float **));
	  voidptr = (void *) flt2ptr;
	  break;
	case VARTOOLS_TYPE_STRING:
	  string2ptr = (char ****) malloc(sizeof(char ***));
	  voidptr = (void *) string2ptr;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  char2ptr = (char ***) malloc(sizeof(char **));
	  voidptr = (void *) char2ptr;
	  break;
	case VARTOOLS_TYPE_USERDEF:
	  voidptr = (void *) malloc(sizeof(void *));
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      }
      RegisterDataFromInputList(p, voidptr, datatype, Ncolumns, -1, 0, 0, 
				NULL, colnum, tmpstring);
      ptr[co->Nptrs].inlistdataptr = voidptr;
      co->Ninlist += 1;
      break;

    case VARTOOLS_SOURCE_COMPUTED:
      /* The array stores data which is computed by the command */
      co->Ncomputed += 1;
      break;

    case VARTOOLS_SOURCE_FIXED:
      /* The array stores data which is fixed for all light curves
	 the value is typically given on the command line */
      if(!co->Nfix) {
	if((co->FixValues = (_UserDataPointer *) malloc(sizeof(_UserDataPointer))) == NULL)
	  error(ERR_MEMALLOC);
      } else {
	if((co->FixValues = (_UserDataPointer *) realloc(co->FixValues, (co->Nfix + 1)*sizeof(_UserDataPointer))) == NULL)
	  error(ERR_MEMALLOC);
      }
      co->FixValues[co->Nfix].datatype = datatype;
      voidptr = va_arg(varlist,void *);

      /* Although the user passes a pointer to a variable storing the 
	 fixed value, we cannot assume that the memory allocated to this
	 variable will be preserved (if this is called explicitly from the
	 UserLib_ParseCL function it may be stored in a local variable only).
	 Therefore we will allocate the memory for the fix value and copy
	 over the data passed from the user to it. */
      if(Ncolumns < 1) {
	switch(datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr1 = (double *) malloc(sizeof(double));
	  *dblptr1 = *((double *) voidptr);
	  voidptr = (void *) dblptr1;
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr1 = (int *) malloc(sizeof(int));
	  *intptr1 = *((int *) voidptr);
	  voidptr = (void *) intptr1;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr1 = (short *) malloc(sizeof(short));
	  *shortptr1 = *((short *) voidptr);
	  voidptr = (void *) shortptr1;
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr1 = (long *) malloc(sizeof(long));
	  *longptr1 = *((long *) voidptr);
	  voidptr = (void *) longptr1;
	case VARTOOLS_TYPE_FLOAT:
	  fltptr1 = (float *) malloc(sizeof(float));
	  *fltptr1 = *((float *) voidptr);
	  voidptr = (void *) fltptr1;
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr1 = (char **) malloc(sizeof(char *));
	  *(stringptr1) = (char *) malloc(MAXLEN*sizeof(char));
	  sprintf((*(stringptr1)),"%s",(*((char **) voidptr)));
	  voidptr = (void *) stringptr1;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr1 = (char *) malloc(sizeof(char));
	  *charptr1 = *((char *) voidptr);
	  voidptr = (void *) charptr1;
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      } else {
	switch(datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = (double **) malloc(sizeof(double *));
	  *dblptr = (double *) malloc(Ncolumns*sizeof(double));
	  for(j=0; j < Ncolumns; j++) {
	    (*dblptr)[j] = (*((double **) voidptr))[j];
	  }
	  voidptr = (void *) dblptr;
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = (int **) malloc(sizeof(int *));
	  *intptr = (int *) malloc(Ncolumns*sizeof(int));
	  for(j=0; j < Ncolumns; j++) {
	    (*intptr)[j] = (*((int **) voidptr))[j];
	  }
	  voidptr = (void *) intptr;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = (short **) malloc(sizeof(short *));
	  *shortptr = (short *) malloc(Ncolumns*sizeof(short));
	  for(j=0; j < Ncolumns; j++) {
	    (*shortptr)[j] = (*((short **) voidptr))[j];
	  }
	  voidptr = (void *) shortptr;
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = (long **) malloc(sizeof(long *));
	  *longptr = (long *) malloc(Ncolumns*sizeof(long));
	  for(j=0; j < Ncolumns; j++) {
	    (*longptr)[j] = (*((long **) voidptr))[j];
	  }
	  voidptr = (void *) longptr;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  fltptr = (float **) malloc(sizeof(float *));
	  *fltptr = (float *) malloc(Ncolumns*sizeof(float));
	  for(j=0; j < Ncolumns; j++) {
	    (*fltptr)[j] = (*((float **) voidptr))[j];
	  }
	  voidptr = (void *) fltptr;
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = (char ***) malloc(sizeof(char **));
	  *stringptr = (char **) malloc(Ncolumns*sizeof(char *));
	  for(j=0; j < Ncolumns; j++) {
	    (*stringptr)[j] = (char *) malloc(MAXLEN*sizeof(char));
	    sprintf((*stringptr)[j],"%s",(*((char ***) voidptr))[j]);
	  }
	  voidptr = (void *) stringptr;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = (char **) malloc(sizeof(char *));
	  *charptr = (char *) malloc(Ncolumns*sizeof(char));
	  for(j=0; j < Ncolumns; j++) {
	    (*charptr)[j] = (*((char **) voidptr))[j];
	  }
	  voidptr = (void *) charptr;
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      }
      co->FixValues[co->Nfix].dataptr = voidptr;
      co->FixValues[co->Nfix].Ncolumns = Ncolumns;
      co->FixValues[co->Nfix].size_element = sizeval;

      ptr[co->Nptrs].Nfixptr = co->Nfix;
      co->Nfix += 1;
      break;

    case VARTOOLS_SOURCE_PRIORCOLUMN:
      /* This vector stores data which comes from a prior command */

      if(Ncolumns >= 1) {
	if((ptr[co->Nptrs].linkedcolumn = (OutColumn **) malloc(Ncolumns*sizeof(OutColumn *))) == NULL ||
	   (ptr[co->Nptrs].priorname = (char **) malloc(Ncolumns*sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
	for(k=0; k < Ncolumns; k++) {
	  if((ptr[co->Nptrs].priorname[k] = (char *) malloc(MAXLEN)) == NULL)
	    error(ERR_MEMALLOC);
	}
      } else {
	if((ptr[co->Nptrs].linkedcolumn = (OutColumn **) malloc(sizeof(OutColumn *))) == NULL ||
	   (ptr[co->Nptrs].priorname = (char **) malloc(sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
	if((ptr[co->Nptrs].priorname[0] = (char *) malloc(MAXLEN)) == NULL)
	  error(ERR_MEMALLOC);
      }

      if(Ncolumns <= 0) {
	tmpstring = va_arg(varlist,char *);

	sprintf(ptr[co->Nptrs].priorname[0],"%s",tmpstring);
      } else {
	tmpstring2 = va_arg(varlist,char **);
	for(k=0; k < Ncolumns; k++) {
	  sprintf(ptr[co->Nptrs].priorname[k],"%s",tmpstring2[k]);
	}
      }

      co->Nprior += 1;
      break;

    case VARTOOLS_SOURCE_LC:
      /* This vector stores data which is read in from the light curves */
      tmpstring = va_arg(varlist,char *);
      colnum = va_arg(varlist,int);
      scanformat = va_arg(varlist,char *);
      varnameout = va_arg(varlist,char *);
      /*if(Ncolumns == 0) {
	switch(datatype){
	case VARTOOLS_TYPE_DOUBLE:
	  dbl2ptr = (double ***) malloc(sizeof(double **));
	  voidptr = (void *) dbl2ptr;
	  break;
	case VARTOOLS_TYPE_INT:
	  int2ptr = (int ***) malloc(sizeof(int **));
	  voidptr = (void *) int2ptr;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  short2ptr = (short ***) malloc(sizeof(short **));
	  voidptr = (void *) short2ptr;
	  break;
	case VARTOOLS_TYPE_LONG:
	  long2ptr = (long ***) malloc(sizeof(long **));
	  voidptr = (void *) long2ptr;
	case VARTOOLS_TYPE_FLOAT:
	  flt2ptr = (float ***) malloc(sizeof(float **));
	  voidptr = (void *) flt2ptr;
	  break;
	case VARTOOLS_TYPE_STRING:
	  string2ptr = (char ****) malloc(sizeof(char ***));
	  voidptr = (void *) string2ptr;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  char2ptr = (char ***) malloc(sizeof(char **));
	  voidptr = (void *) char2ptr;
	  break;
	case VARTOOLS_TYPE_USERDEF:
	  voidptr = (void *) malloc(sizeof(void *));
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      }
      else {
	switch(datatype){
	case VARTOOLS_TYPE_DOUBLE:
	  dbl3ptr = (double ****) malloc(sizeof(double ***));
	  voidptr = (void *) dbl3ptr;
	  break;
	case VARTOOLS_TYPE_INT:
	  int3ptr = (int ****) malloc(sizeof(int ***));
	  voidptr = (void *) int3ptr;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  short3ptr = (short ****) malloc(sizeof(short ***));
	  voidptr = (void *) short3ptr;
	  break;
	case VARTOOLS_TYPE_LONG:
	  long3ptr = (long ****) malloc(sizeof(long ***));
	  voidptr = (void *) long3ptr;
	case VARTOOLS_TYPE_FLOAT:
	  flt3ptr = (float ****) malloc(sizeof(float ***));
	  voidptr = (void *) flt3ptr;
	  break;
	case VARTOOLS_TYPE_STRING:
	  string3ptr = (char *****) malloc(sizeof(char ****));
	  voidptr = (void *) string3ptr;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  char3ptr = (char ****) malloc(sizeof(char ***));
	  voidptr = (void *) char3ptr;
	  break;
	case VARTOOLS_TYPE_USERDEF:
	  voidptr = (void *) malloc(sizeof(void *));
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      }
      RegisterDataFromLightCurve(p, voidptr, datatype, MAXLEN, Ncolumns, -1, 
      0, 0, scanformat, NULL, colnum, tmpstring);*/
      RegisterDataFromLightCurve(p, dataptr, datatype, MAXLEN, Ncolumns, -1, 
				 0, 0, scanformat, NULL, colnum, tmpstring);
      ptr[co->Nptrs].inlcdataptr = voidptr;
      if(varnameout != NULL) {
	CreateVariable(p, varnameout, (char) datatype, VARTOOLS_VECTORTYPE_LC, dataptr);
      }
      co->Ninlc += 1;
      break;
      
    case VARTOOLS_SOURCE_RECENTCOMMAND: 
      /* This vector stores data which comes from a prior command, we need
       to look up the column name for the desired command */

      if(Ncolumns >= 1) {
	fprintf(stderr,"Error: Bad call to RegisterDataVector, Ncolumns must be <= 0 to use VARTOOLS_SOURCE_RECENTCOMMAND.\n");
	exit(ERR_CODEERROR);
      }

      if((ptr[co->Nptrs].linkedcolumn = (OutColumn **) malloc(sizeof(OutColumn *))) == NULL ||
	 (ptr[co->Nptrs].priorname = (char **) malloc(sizeof(char *))) == NULL)
	error(ERR_MEMALLOC);
      if((ptr[co->Nptrs].priorname[0] = (char *) malloc(MAXLEN)) == NULL)
	error(ERR_MEMALLOC);

      colnum = va_arg(varlist,int);
      charptr1 = va_arg(varlist,char *);

      tmpstring = malloc(MAXLEN);
      GetColumnNameForRecentCommand(p, c, colnum, charptr1, 
				    ptr[co->Nptrs].priorname[0]);
      sprintf(ptr[co->Nptrs].priorname[0],"%s",tmpstring);
      
      co->Nprior += 1;
      break;

    case VARTOOLS_SOURCE_EVALEXPRESSION:

      if(co->Nexpr == 0) {
	if((co->expr_strings = (char **) malloc(sizeof(char *))) == NULL ||
	   (co->UserDataExpressions = (_Expression ***) malloc(sizeof(_Expression **))) == NULL)
	  error(ERR_MEMALLOC);
      } else {
	if((co->expr_strings = (char **) realloc(co->expr_strings, (co->Nexpr + 1)*sizeof(char *))) == NULL ||
	   (co->UserDataExpressions = (_Expression ***) realloc(co->UserDataExpressions, (co->Nexpr + 1)*sizeof(_Expression **))) == NULL)
	  error(ERR_MEMALLOC);
      }
      co->UserDataExpressions[co->Nexpr] = &(ptr[co->Nptrs].evalexpression);
      exprstring = va_arg(varlist,char *);
      if((co->expr_strings[co->Nexpr] = malloc(sizeof(exprstring)+1)) == NULL)
	error(ERR_MEMALLOC);
      sprintf(co->expr_strings[co->Nexpr],"%s",exprstring);

      co->Nexpr += 1;
      break;

    case VARTOOLS_SOURCE_EVALEXPRESSION_LC:

      if(co->Nexpr == 0) {
	if((co->expr_strings = (char **) malloc(sizeof(char *))) == NULL ||
	   (co->UserDataExpressions = (_Expression ***) malloc(sizeof(_Expression **))) == NULL)
	  error(ERR_MEMALLOC);
      } else {
	if((co->expr_strings = (char **) realloc(co->expr_strings, (co->Nexpr + 1)*sizeof(char *))) == NULL ||
	   (co->UserDataExpressions = (_Expression ***) realloc(co->UserDataExpressions, (co->Nexpr + 1)*sizeof(_Expression **))) == NULL)
	  error(ERR_MEMALLOC);
      }
      co->UserDataExpressions[co->Nexpr] = &(ptr[co->Nptrs].evalexpression);
      exprstring = va_arg(varlist,char *);
      if((co->expr_strings[co->Nexpr] = malloc(sizeof(exprstring)+1)) == NULL)
	error(ERR_MEMALLOC);
      sprintf(co->expr_strings[co->Nexpr],"%s",exprstring);

      co->Nexpr += 1;

      varnameout = va_arg(varlist,char *);
      RegisterDataFromLightCurve(p, dataptr, datatype, MAXLEN, Ncolumns, -1, 
				 0, 0, NULL, NULL, -1, NULL);
      ptr[co->Nptrs].inlcdataptr = voidptr;
      if(varnameout != NULL) {
	CreateVariable(p, varnameout, (char) datatype, VARTOOLS_VECTORTYPE_LC, dataptr);
      }
      break;

    case VARTOOLS_SOURCE_EXISTINGVARIABLE:
      
      if(c->N_prior_vars == 0) {
	if((c->prior_var_datatypes = (char *) malloc(sizeof(char))) == NULL ||
	   (c->prior_var_vectortypes = (char *) malloc(sizeof(char))) == NULL ||
	   (c->prior_var_names = (char **) malloc(sizeof(char *))) == NULL ||
	   (c->prior_vars = (_Variable ***) malloc(sizeof(_Variable **))) == NULL)
	  error(ERR_MEMALLOC);
      } else {
	if((c->prior_var_datatypes = (char *) realloc(c->prior_var_datatypes, (c->N_prior_vars + 1)*sizeof(char))) == NULL ||
	   (c->prior_var_vectortypes = (char *) realloc(c->prior_var_vectortypes, (c->N_prior_vars + 1)*sizeof(char))) == NULL ||
	   (c->prior_var_names = (char **) realloc(c->prior_var_names, (c->N_prior_vars + 1)*sizeof(char *))) == NULL ||
	   (c->prior_vars = (_Variable ***) realloc(c->prior_vars, (c->N_prior_vars + 1)*sizeof(_Variable **))) == NULL)
	  error(ERR_MEMALLOC);
      }
      exprstring = va_arg(varlist,char *);
      if((c->prior_var_names[c->N_prior_vars] = (char *) malloc(strlen(exprstring)+1)) == NULL)
	error(ERR_MEMALLOC);
      sprintf(c->prior_var_names[c->N_prior_vars],"%s",exprstring);
      c->prior_var_vectortypes[c->N_prior_vars] = (char) va_arg(varlist,int);
      c->prior_vars[c->N_prior_vars] = va_arg(varlist,_Variable **);
      ptr[co->Nptrs].existingvariable = c->prior_vars[c->N_prior_vars];
      ptr[co->Nptrs].expectedvectortype = c->prior_var_vectortypes[c->N_prior_vars];
      c->prior_var_datatypes[c->N_prior_vars] = datatype;
      if(c->prior_var_vectortypes[c->N_prior_vars] == VARTOOLS_VECTORTYPE_LC &&
	 output) {
	error(ERR_BADTYPE);
      }

      c->N_prior_vars += 1;
      break;

    default:
      error(ERR_BADTYPE);
      break;
    }

  co->Nptrs += 1;  
#else
  return;
#endif
}  


int ParseFixSpecFixcolumn(ProgramData *p, Command *c, int *iret, char **argv,
			   int argc, int Nvec, _ParseFixSpecFixcolumnStruct *s)
/* Use this function to parse a command of the form:
   <"fix" value | "list" [\"column\" col] | "fixcolumn" <colname | colnum> | 
    "expr" expression>

   p - pointer to the ProgramData struct

   c - pointer to the Command struct associated with this UserCommand

   iret - pointer to a variable storing the current command line index
          number.  On input this should be the term storing "fix",
          "list", or "fixcolumn". This is updated on output.

   argv - Standard array holding the command-line strings.

   argc - Number of terms on the CL.

   Nvec - Number of terms to parse

   .... - For each term to parse, the following additional arguments must be
          provided:

         datatype - VARTOOLS_TYPE_DOUBLE, VARTOOLS_TYPE_FLOAT,
                    VARTOOLS_TYPE_INT, VARTOOLS_TYPE_SHORT,
                    VARTOOLS_TYPE_LONG, VARTOOLS_TYPE_STRING,
                    VARTOOLS_TYPE_CHAR

         dataptr - Pointer to the vector or array that will store the input
                   data.

         Ncolumns - Number of columns in the array, if it is <= 0 then
                    dataptr is taken to be a vector rather than an array
                    (e.g. dataptr is a pointer to an variable of type double*
                     rather than double**).

         output - 1 if this data should be included in the output ascii table,
                  0 if it should not be.

         name - root name of the data vector for display in the output
             ascii table and for the input table (if the user gives
             the \"list\" keyword on the command line, and issues the
             -inputlistformat command). The name of the library will
             be added as prefix to the name, and the column number (if
             Ncolumns >= 1) and command numbers will be appended to
             the name.  */
{
#ifdef DYNAMICLIB
  int datatype;
  void *dataptr;
  int Ncolumns;
  int output;
  char *name;
  char **priornames;
  int sizepriornames;
  int i, j, k;
  int Nc, incol;

  void *fixptr;
  double *dblptr;
  float *fltptr;
  int *intptr;
  short *shortptr;
  long *longptr;
  char *charptr;
  char **stringptr;

  double **dbl2ptr;
  float **flt2ptr;
  int **int2ptr;
  short **short2ptr;
  long **long2ptr;
  char **char2ptr;
  char ***string2ptr;

  double dblval;
  float fltval;
  int intval;
  short shortval;
  long longval;
  char charval;
  char stringval[MAXLEN];

  sizepriornames = 0;

  if(Nvec <= 0)
    error(ERR_INVALIDUSEOFPARSEFIXSPECFIXCOLUMN);
 
  i = *iret;

  for(j=0; j < Nvec; j++) {

    datatype = s[j].datatype;
    dataptr = s[j].dataptr;
    Ncolumns = s[j].Ncolumns;
    output = s[j].output;
    name = s[j].name;

    if(i < argc) {
      if(!strcmp(argv[i],"fix")) {
	if(Ncolumns <= 0) {
	  i++;
	  if(i < argc) {
	    switch(datatype) {
	    case VARTOOLS_TYPE_DOUBLE:
	      dblval = atof(argv[i]);
	      fixptr = (void *) &dblval;
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      fltval = (float) atof(argv[i]);
	      fixptr = (void *) &fltval;
	      break;
	    case VARTOOLS_TYPE_INT:
	      intval = atoi(argv[i]);
	      fixptr = (void *) &intval;
	      break;
	    case VARTOOLS_TYPE_LONG:
	      longval = atol(argv[i]);
	      fixptr = (void *) &longval;
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      shortval = (short) atoi(argv[i]);
	      fixptr = (void *) &shortval;
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      charval = argv[i][0];
	      fixptr = (void *) &charval;
	      break;
	    case VARTOOLS_TYPE_STRING:
	      sprintf(stringval,"%s",argv[i]);
	      fixptr = (void *) &stringval;
	      break;
	    default:
	      error(ERR_BADTYPE);
	      break;
	    }
	    i++;
	  } else {
	    return 1;
	  }
	} else {
	  switch(datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    if((dbl2ptr = (double **) malloc(sizeof(double *))) == NULL)
	      error(ERR_MEMALLOC);
	    if((*dbl2ptr = (double *) malloc(Ncolumns*sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	    fixptr = dbl2ptr;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    if((flt2ptr = (float **) malloc(sizeof(float *))) == NULL)
	      error(ERR_MEMALLOC);
	    if((*flt2ptr = (float *) malloc(Ncolumns*sizeof(float))) == NULL)
	      error(ERR_MEMALLOC);
	    fixptr = flt2ptr;
	    break;
	  case VARTOOLS_TYPE_INT:
	    if((int2ptr = (int **) malloc(sizeof(int *))) == NULL)
	      error(ERR_MEMALLOC);
	    if((*int2ptr = (int *) malloc(Ncolumns*sizeof(int))) == NULL)
	      error(ERR_MEMALLOC);
	    fixptr = int2ptr;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    if((long2ptr = (long **) malloc(sizeof(long *))) == NULL)
	      error(ERR_MEMALLOC);
	    if((*long2ptr = (long *) malloc(Ncolumns*sizeof(long))) == NULL)
	      error(ERR_MEMALLOC);
	    fixptr = long2ptr;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    if((short2ptr = (short **) malloc(sizeof(short *))) == NULL)
	      error(ERR_MEMALLOC);
	    if((*short2ptr = (short *) malloc(Ncolumns*sizeof(short))) == NULL)
	      error(ERR_MEMALLOC);
	    fixptr = short2ptr;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    if((char2ptr = (char **) malloc(sizeof(char *))) == NULL)
	      error(ERR_MEMALLOC);
	    if((*char2ptr = (char *) malloc(Ncolumns)) == NULL)
	      error(ERR_MEMALLOC);
	    fixptr = char2ptr;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    if((string2ptr = (char ***) malloc(sizeof(char **))) == NULL)
	      error(ERR_MEMALLOC);
	    if((*string2ptr = (char **) malloc(Ncolumns*sizeof(char *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=0; k < Ncolumns; k++) {
	      if(((*string2ptr)[k] = (char *) malloc(MAXLEN)) == NULL)
		error(ERR_MEMALLOC);
	    }
	    fixptr = string2ptr;
	    break;
	  default:
	    error(ERR_BADTYPE);
	    break;
	  }
	  for(k=0; k < Ncolumns; k++) {
	    i++;
	    if(i < argc) {
	      switch(datatype){
	      case VARTOOLS_TYPE_DOUBLE:
		(*dbl2ptr)[k] = atof(argv[i]);
		break;
	      case VARTOOLS_TYPE_FLOAT:
		(*flt2ptr)[k] = (float) atof(argv[i]);
		break;
	      case VARTOOLS_TYPE_INT:
		(*int2ptr)[k] = atoi(argv[i]);
		break;
	      case VARTOOLS_TYPE_LONG:
		(*long2ptr)[k] = atol(argv[i]);
		break;
	      case VARTOOLS_TYPE_SHORT:
		(*short2ptr)[k] = (short) atoi(argv[i]);
		break;
	      case VARTOOLS_TYPE_CHAR:
		(*char2ptr)[k] = argv[i][0];
		break;
	      case VARTOOLS_TYPE_STRING:
		sprintf((*string2ptr)[k],"%s",argv[i]);
		break;
	      default:
		error(ERR_BADTYPE);
		break;
	      }
	    } else {
	      return 1;
	    }
	  }
	  i++;
	}
	RegisterDataVector(p, c, dataptr, datatype, Ncolumns,
			   VARTOOLS_SOURCE_FIXED, output, name,
			   fixptr);
      }
      else if(!strcmp(argv[i],"expr")) {
	i++;
	if(i < argc) {
	  RegisterDataVector(p, c, dataptr, datatype, Ncolumns,
			     VARTOOLS_SOURCE_EVALEXPRESSION, output, name,
			     argv[i]);
	  i++;
	}
	else
	  return 1;
      }
      else if(!strcmp(argv[i],"list")) {
	i++;
	if(i < argc) {
	  if(!strcmp(argv[i],"column")) {
	    i++;
	    if(i < argc) {
	      incol = atoi(argv[i]);
	      RegisterDataVector(p, c, dataptr, datatype, Ncolumns,
				 VARTOOLS_SOURCE_INLIST, output, name,
				 name, incol);
	      i++;
	    } else 
	      return 1;
	  } 
	  else {
	    RegisterDataVector(p, c, dataptr, datatype, Ncolumns,
			       VARTOOLS_SOURCE_INLIST, output, name,
			       name, 0);
	  }
	} else {
	  RegisterDataVector(p, c, dataptr, datatype, Ncolumns,
			     VARTOOLS_SOURCE_INLIST, output, name,
			     name, 0);
	}
      } else if(!strcmp(argv[i],"fixcolumn")) {
	if(Ncolumns <= 0) {
	  if(!sizepriornames) {
	    if((priornames = (char **) malloc(sizeof(char *))) == NULL)
	      error(ERR_MEMALLOC);
	    if((priornames[0] = (char *) malloc(MAXLEN)) == NULL)
	      error(ERR_MEMALLOC);
	    sizepriornames = 1;
	  }
	} else {
	  if(!sizepriornames) {
	    if((priornames = (char **) malloc(Ncolumns*sizeof(char *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=0; k < Ncolumns; k++) {
	      if((priornames[k] = (char *) malloc(MAXLEN)) == NULL)
		error(ERR_MEMALLOC);
	    }
	    sizepriornames = Ncolumns;
	  } 
	  else if(sizepriornames > Ncolumns) {
	    if((priornames = (char **) realloc(priornames, Ncolumns*sizeof(char *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=sizepriornames; k < Ncolumns; k++) {
	      if((priornames[k] = (char *) malloc(MAXLEN)) == NULL)
		error(ERR_MEMALLOC);
	    }
	    sizepriornames = Ncolumns;
	  }
	}
	for(k=0; k < (Ncolumns <= 0 ? 1 : Ncolumns); k++) {
	  i++;
	  if(i < argc) {
	    sprintf(priornames[k],"%s",argv[i]);
	  } else {
	    for(k=0; k < sizepriornames; k++)
	      free(priornames[k]);
	    free(priornames);
	    return 1;
	  }
	}
	i++;
	if(Ncolumns <= 0) {
	  RegisterDataVector(p, c, dataptr, datatype, Ncolumns,
			     VARTOOLS_SOURCE_PRIORCOLUMN, output, name,
			     priornames[0]);
	} else {
	  RegisterDataVector(p, c, dataptr, datatype, Ncolumns,
			     VARTOOLS_SOURCE_PRIORCOLUMN, output, name,
			     priornames);
	}
      }
      else
	return 1;
    }
    else
      return 1;
  }

  if(sizepriornames > 0) {
    for(k=0; k < sizepriornames; k++)
      free(priornames[k]);
    free(priornames);
  }
    
  *iret = i;
    
  return 0;
#else
  return 1;
#endif
}

int ParseOutNameKeyword(ProgramData *p, Command *c, int *iret, char **argv,
			int argc, const char *keyword, int *outputflag,
			char **outdir, int *formatflag, char **format)
/* Use this function to parse a command of the form:
   "keyword" outdir ["nameformat" format]

   p - pointer to the ProgramData struct

   c - pointer to the Command struct associated with this UserCommand

   iret - pointer to a variable storing the current command line index
          number.  On input this should be the term storing
          "keyword". This is updated on output.

   argv - Standard array holding the command-line strings.

   argc - Number of terms on the CL.

   keyword - the keyword to check. If this is NULL, then the keyword will
             not be checked.

   outputflag - integer flag set to 1 if the command is parsed or 0 if not.

   outdir - pointer to a char* variable to store the name of the output
            directory. Memory will be allocated to store the string by this
            command. This will be returned as NULL if the keyword is not
            given.

   formatflag - integer flag set to 1 if the user gave the "nameformat" keyword,
           or 0 if not.

   format - pointer to a char* variable to store the format string. Memory
            will be allocated to store the string by this command. This will
            be returned as NULL if the keyword is not given.

Return values:
      0 - command is successfully parsed.
      1 - "keyword" is not given, command is not parsed.
      2 - "keyword" is given, outdir is not given. Or "nameformat" is given but
          the output format is not given.
*/
{
  *outputflag = 0; *formatflag = 0;
  *outdir = NULL; *format = NULL;
  if(keyword != NULL) {
    if((*iret >= argc))
      return 1;
    if(strcmp(argv[*iret],keyword))
      return 1;
    else
      (*iret)++;
  }
  if(*iret >= argc)
    return 2;
  if(((*outdir) = (char *) malloc((strlen(argv[*iret])+1)*sizeof(char))) == NULL)
    error(ERR_MEMALLOC);
  sprintf((*outdir),"%s",argv[*iret]);
  *outputflag = 1;
  (*iret)++;
  if(*iret >= argc)
    return 0;
  if(strcmp(argv[*iret],"nameformat"))
    return 0;
  (*iret)++;
  if(*iret >= argc)
    return 2;
  if(((*format) = (char *) malloc((strlen(argv[*iret])+1)*sizeof(char))) == NULL)
    error(ERR_MEMALLOC);
  sprintf((*format),"%s",argv[*iret]);
  *formatflag = 1;
  return 0;
}		        

int ParseConstantParameter(ProgramData *p, Command *c, int *iret, char **argv,
			   int argc, const char *keyword, char datatype,
			   void *dataptr, int Ncolumns)
/* Use this function to parse a command of the form:
   "keyword" value [value2 value3 ... valueNvec]

   or of the form:
   value [value2 value3 ... valueNvec]

   p - pointer to the ProgramData struct

   c - pointer to the Command struct associated with this UserCommand

   iret - pointer to a variable storing the current command line index
          number.  On input this should be the term storing "keyword"
          or value if the keyword is not used. This is updated on
          output.

   argv - Standard array holding the command-line strings.

   argc - Number of terms on the CL.

   keyword - the keyword to check. If this is NULL, then the keyword will
             not be checked.

   datatype - VARTOOLS_TYPE_DOUBLE, VARTOOLS_TYPE_FLOAT,
              VARTOOLS_TYPE_INT, VARTOOLS_TYPE_SHORT,
              VARTOOLS_TYPE_LONG, VARTOOLS_TYPE_STRING,
              VARTOOLS_TYPE_CHAR

   dataptr - Pointer to the variable or vector that will store the
             input data. If Ncolumns <= 0, you should declare a
             variable like "double dblval", and then pass &dblval to
             this function. If Ncolumns > 0, you should declare
             "double *dblval", and then pass &dblval to this
             function. Note that memory will be allocated for the
             vector by this routine if the command is successfully
             parsed.  If it is not parsed, then memory will not be
             allocated. If the datatype is a string, then memory to
             store the string will be allocated whether the data is a
             scalar or a vector.  In otherwords, if Ncolumn <= 0, and
             you are reading a string into a variable stringvar. You
             would declare stringvar as "char *stringvar;" rather than
             as "char stringvar[sizestring];", then pass &stringvar to
             this function for dataptr. If Ncolumn > 0, you would
             declare stringvar as "char **stringvar", and pass
             &stringvar to this function for dataptr.

   Ncolumns - Number of parameters to read in. If it is <= 0 then one
              parameter is read in, and dataptr is a pointer to a
              scalar variable. If it is > 0, then Ncolumns terms are
              read in, and stored in the vector pointed to be
              dataptr. 

   Return values:
      0 - command is successfully parsed.
      1 - "keyword" is not given, command is not parsed.
      2 - "keyword" is given, but the rest of the command is not parsed
          successfully.
*/
{
  double *dblptr;
  float *fltptr;
  int *intptr;
  long *longptr;
  short *shortptr;
  char *charptr;
  char **stringptr;
  
  int i;

  if(keyword != NULL) {
    if(*iret >= argc)
      return 1;
    if(strcmp(argv[*iret],keyword))
      return 1;
    (*iret)++;
  }
  
  if((*iret) + Ncolumns > argc)
    return 2;

  /* Allocate memory for the data-vector if the parameter is not a scalar */
  if(Ncolumns > 0) {
    switch(datatype) {
    case VARTOOLS_TYPE_DOUBLE:
      if((dblptr = (double *) malloc(Ncolumns * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0; i < Ncolumns; i++) {
	dblptr[i] = atof(argv[*iret]);
	(*iret)++;
      }
      *((double **) dataptr) = dblptr;
      break;
    case VARTOOLS_TYPE_FLOAT:
      if((fltptr = (float *) malloc(Ncolumns * sizeof(float))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0; i < Ncolumns; i++) {
	fltptr[i] = (float) atof(argv[*iret]);
	(*iret)++;
      }
      *((float **) dataptr) = fltptr;
      break;
    case VARTOOLS_TYPE_STRING:
      if((stringptr = (char **) malloc(Ncolumns * sizeof(char *))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0; i < Ncolumns; i++) {
	if((stringptr[i] = (char *) malloc((strlen(argv[*iret])+1)*sizeof(char))) == NULL)
	  error(ERR_MEMALLOC);
	sprintf(stringptr[i],"%s",argv[*iret]);
	(*iret)++;
      }
      *((char ***) dataptr) = stringptr;
      break;
    case VARTOOLS_TYPE_CHAR:
      if((charptr = (char *) malloc(Ncolumns * sizeof(char))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0; i < Ncolumns; i++) {
	charptr[i] = argv[*iret][0];
	(*iret)++;
      }
      *((char **) dataptr) = charptr;
      break;
    case VARTOOLS_TYPE_INT:
      if((intptr = (int *) malloc(Ncolumns * sizeof(int))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0; i < Ncolumns; i++) {
	intptr[i] = atoi(argv[*iret]);
	(*iret)++;
      }
      *((int **) dataptr) = intptr;
      break;
    case VARTOOLS_TYPE_SHORT:
      if((shortptr = (short *) malloc(Ncolumns * sizeof(short))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0; i < Ncolumns; i++) {
	shortptr[i] = (short) atoi(argv[*iret]);
	(*iret)++;
      }
      *((short **) dataptr) = shortptr;
      break;
    case VARTOOLS_TYPE_LONG:
      if((longptr = (long *) malloc(Ncolumns * sizeof(long))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0; i < Ncolumns; i++) {
	longptr[i] = atol(argv[*iret]);
	(*iret)++;
      }
      *((long **) dataptr) = longptr;
      break;
    default:
      error(ERR_BADTYPE);
    }
  } else {
    switch(datatype) {
    case VARTOOLS_TYPE_DOUBLE:
      *((double *) dataptr) = atof(argv[*iret]);
      break;
    case VARTOOLS_TYPE_FLOAT:
      *((float *) dataptr) = (float) atof(argv[*iret]);
      break;
    case VARTOOLS_TYPE_STRING:
      if((charptr = (char *) malloc((strlen(argv[*iret])+1)*sizeof(char))) == NULL)
	error(ERR_MEMALLOC);
      sprintf(charptr,"%s",argv[*iret]);
      *((char **) dataptr) = charptr;
      break;
    case VARTOOLS_TYPE_CHAR:
      *((char *) dataptr) = argv[*iret][0];
      break;
    case VARTOOLS_TYPE_INT:
      *((int *) dataptr) = atoi(argv[*iret]);
      break;
    case VARTOOLS_TYPE_SHORT:
      *((short *) dataptr) = (short) atoi(argv[*iret]);
      break;
    case VARTOOLS_TYPE_LONG:
      *((long *) dataptr) = atol(argv[*iret]);
      break;
    default:
      error(ERR_BADTYPE);
    }
    (*iret)++;
  }
  return 0;
  
}

int ParseParameterBuiltInCommand(ProgramData *p, int cnum, 
				 int *iret, char **argv,
				 int argc, const char *keyword, int Nvec,
				 ...)
/* Use this function for built-in VARTOOLS commands to parse a command of the
   form:
   "keyword" <"fix" value | "list" [\"column\" col] | "fixcolumn" <colname | colnum> | "expr" expression>

   p - pointer to the ProgramData struct

   cnum - the command number.

   iret - pointer to a variable storing the current command line index
          number.  On input this should be the term storing "keyword".
          This is updated on output.

   argv - Standard array holding the command-line strings.

   argc - Number of terms on the CL.

   keyword - the keyword to check. If this is NULL then the check for the
          keyword is skipped (the return value is 0 or 2 in this case).

   Nvec - Number of terms to parse

   .... - For each term to parse, the following additional arguments must be
          provided:

       int datatype - VARTOOLS_TYPE_DOUBLE, VARTOOLS_TYPE_FLOAT,
                  VARTOOLS_TYPE_INT, VARTOOLS_TYPE_SHORT,
                  VARTOOLS_TYPE_LONG, VARTOOLS_TYPE_STRING,
                  VARTOOLS_TYPE_CHAR

       int *source - Pointer to a variable which will store the user's 
                selection. On output it will have the values:
                VARTOOLS_SOURCE_FIXED, VARTOOLS_SOURCE_PRIORCOLUMN,
                VARTOOLS_SOURCE_INLIST, VARTOOLS_SOURCE_EVALEXPRESSION

       int Ncolumns - Number of columns in the array, if it is <= 0 then
                      listptr is taken to be a vector rather than an array
                      (e.g. listptr is a pointer to an variable of type double*
                      rather than double**).

       void *fixval - Pointer to a variable where the fixed value should 
                be stored if selected by the user. If Ncolumns >= 1, this 
                should be a pointer to a vector.

       void *dataptr - If Ncolumns >= 1 this should be a pointer to an array
                       otherwise it should be a pointer to a vector.

       void *outcolumn - If Ncolumns >= 1 this should be a pointer to a
                         variable of type **OutColumn, otherwise it should
                         be a pointer to a variable of type *OutColumn.

       void *expr_char - If Ncolumns >= 1 this should be a pointer to a
                         variable of type **char, otherwise it should be
                         a pointer to a variable of type *char.

       name - root name of the data vector for display in the input table 
             (if the user gives the \"list\" keyword on the command line, 
              and issues the -inputlistformat command).

       initialize - 1 if the data vector/array should be created and fixed
                    to a default value when the "keyword" is not given on 
                    the command-line. 
                    0 if the data vector/array should not be created.

       defaultval - the default value to initialize the array to.


   Return values:
      0 - command is successfully parsed.
      1 - "keyword" is not given, command is not parsed.
      2 - "keyword" is given, but the rest of the command is not parsed.
*/
{

  va_list varlist;
  int i, j, k;
  int datatype, *source, Ncolumns;
  void *fixval, *dataptr, *outcolumn, *expr_char;
  char *name;
  int initialize;
  void *defaultval;
  int incol;

  char **priornames = NULL;
  int sizepriornames = 0;

  if(keyword != NULL) {
    if(*iret >= argc) return 1;
    if(strcmp(argv[*iret],keyword)) {
      va_start(varlist, Nvec);
      for(i=0; i < Nvec; i++) {
	datatype = ((int) va_arg(varlist,int));
	source = ((int *) va_arg(varlist,int *));
	Ncolumns = ((int) va_arg(varlist,int));
	fixval = ((void *) va_arg(varlist,void *));
	dataptr = ((void *) va_arg(varlist,void *));
	outcolumn = ((void *) va_arg(varlist,void *));
	expr_char = ((void *) va_arg(varlist,void *));
	name = ((char *) va_arg(varlist,char *));
	initialize = ((int) va_arg(varlist,int));
	defaultval = ((void *) va_arg(varlist,void *));
	
	if(initialize) {
	  *source = VARTOOLS_SOURCE_FIXED;
	  if(Ncolumns <= 0) {
	    switch(datatype) {
	    case VARTOOLS_TYPE_DOUBLE:
	      *((double *) fixval) = *((double *) defaultval);
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      *((float *) fixval) = *((float *) defaultval);
	      break;
	    case VARTOOLS_TYPE_INT:
	      *((int *) fixval) = *((int *) defaultval);
	      break;
	    case VARTOOLS_TYPE_LONG:
	      *((long *) fixval) = *((long *) defaultval);
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      *((short *) fixval) = *((short *) defaultval);
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      *((char *) fixval) = *((char *) defaultval);
	      break;
	    case VARTOOLS_TYPE_STRING:
	      if(((*((char **) fixval)) = (char *) malloc(strlen(*((char **) defaultval))+1)) == NULL)
		error(ERR_MEMALLOC);
	      sprintf(*((char **) fixval),"%s", (*((char **)defaultval)));
	      break;
	    default:
	      error(ERR_BADTYPE);
	      break;
	    }
	  } else {
	    switch(datatype) {
	    case VARTOOLS_TYPE_DOUBLE:
	      if((*((double **) fixval) = (double *) malloc(Ncolumns*sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=0; j < Ncolumns; j++)
		(*((double **) fixval))[j] = (*((double **) defaultval))[j];
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      if((*((float **) fixval) = (float *) malloc(Ncolumns*sizeof(float))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=0; j < Ncolumns; j++)
		(*((float **) fixval))[j] = (*((float **) defaultval))[j];
	      break;
	    case VARTOOLS_TYPE_INT:
	      if((*((int **) fixval) = (int *) malloc(Ncolumns*sizeof(int))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=0; j < Ncolumns; j++)
		(*((int **) fixval))[j] = (*((int **) defaultval))[j];
	      break;
	    case VARTOOLS_TYPE_LONG:
	      if((*((long **) fixval) = (long *) malloc(Ncolumns*sizeof(long))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=0; j < Ncolumns; j++)
		(*((long **) fixval))[j] = (*((long **) defaultval))[j];
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      if((*((short **) fixval) = (short *) malloc(Ncolumns*sizeof(short))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=0; j < Ncolumns; j++)
		(*((short **) fixval))[j] = (*((short **) defaultval))[j];
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      if((*((char **) fixval) = (char *) malloc(Ncolumns*sizeof(char))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=0; j < Ncolumns; j++)
		(*((char **) fixval))[j] = (*((char **) defaultval))[j];
	      break;
	    case VARTOOLS_TYPE_STRING:
	      if((*((char ***) fixval) = (char **) malloc(Ncolumns*sizeof(char))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=0; j < Ncolumns; j++) {
		if(((*((char ***) fixval))[j] = (char *) malloc(strlen((*((char ***) defaultval))[j])+1)) == NULL)
		  error(ERR_MEMALLOC);
		sprintf((*((char ***) fixval))[j],"%s", (*((char ***)defaultval))[j]);
		break;
	      default:
		error(ERR_BADTYPE);
		break;
	      }
	    }
	  }
	}
      }
      va_end(varlist);
      return 1;
    }
    *iret = *iret + 1;
  }
  if(Nvec <= 0)
    return 0;

  va_start(varlist, Nvec);
  for(i=0; i < Nvec; i++) {
    datatype = ((int) va_arg(varlist,int));
    source = ((int *) va_arg(varlist,int *));
    Ncolumns = ((int) va_arg(varlist,int));
    fixval = ((void *) va_arg(varlist,void *));
    dataptr = ((void *) va_arg(varlist,void *));
    outcolumn = ((void *) va_arg(varlist,void *));
    expr_char = ((void *) va_arg(varlist,void *));
    name = ((char *) va_arg(varlist,char *));
    initialize = ((int) va_arg(varlist,int));
    defaultval = ((void *) va_arg(varlist,void *));

    if(*iret < argc) {
      if(!strcmp(argv[*iret],"fix")) {
	*source = VARTOOLS_SOURCE_FIXED;
	if(Ncolumns <= 0) {
	  *iret = *iret + 1;
	  if(*iret < argc) {
	    switch(datatype) {
	    case VARTOOLS_TYPE_DOUBLE:
	      *((double *) fixval) = atof(argv[*iret]);
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      *((float *) fixval)  = (float) atof(argv[*iret]);
	      break;
	    case VARTOOLS_TYPE_INT:
	      *((int *) fixval) = atoi(argv[*iret]);
	      break;
	    case VARTOOLS_TYPE_LONG:
	      *((long *) fixval) = atol(argv[*iret]);
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      *((short *) fixval) = (short) atoi(argv[*iret]);
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      *((char *) fixval) = argv[*iret][0];
	      break;
	    case VARTOOLS_TYPE_STRING:
	      if((*((char **) fixval) = (char *) malloc(strlen(argv[*iret])+1)) == NULL)
		error(ERR_MEMALLOC);
	      sprintf(*((char **) fixval),"%s",argv[*iret]);
	      break;
	    default:
	      error(ERR_BADTYPE);
	      break;
	    }
	    *iret = *iret + 1;
	  } else {
	    va_end(varlist);
	    return 1;
	  }
	} else {
	  switch(datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    if((*((double **) fixval) = (double *) malloc(Ncolumns*sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    if((*((float **) fixval) = (float *) malloc(Ncolumns*sizeof(float))) == NULL)
	      error(ERR_MEMALLOC);
	    break;
	  case VARTOOLS_TYPE_INT:
	    if((*((int **) fixval) = (int *) malloc(Ncolumns*sizeof(int))) == NULL)
	      error(ERR_MEMALLOC);
	    break;
	  case VARTOOLS_TYPE_LONG:
	    if((*((long **) fixval) = (long *) malloc(Ncolumns*sizeof(long))) == NULL)
	      error(ERR_MEMALLOC);
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    if((*((short **) fixval) = (short *) malloc(Ncolumns*sizeof(short))) == NULL)
	      error(ERR_MEMALLOC);
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    if((*((char **) fixval) = (char *) malloc(Ncolumns*sizeof(char))) == NULL)
	      error(ERR_MEMALLOC);
	    break;
	  case VARTOOLS_TYPE_STRING:
	    if((*((char ***) fixval) = (char **) malloc(Ncolumns*sizeof(char *))) == NULL)
	      error(ERR_MEMALLOC);
	    break;
	  default:
	    error(ERR_BADTYPE);
	    break;
	  }
	  for(k=0; k < Ncolumns; k++) {
	    *iret = *iret + 1;
	    if(*iret < argc) {
	      switch(datatype){
	      case VARTOOLS_TYPE_DOUBLE:
		(*((double **) fixval))[k] = atof(argv[*iret]);
		break;
	      case VARTOOLS_TYPE_FLOAT:
		(*((float **) fixval))[k] = (float) atof(argv[*iret]);
		break;
	      case VARTOOLS_TYPE_INT:
		(*((int **) fixval))[k] = (int) atof(argv[*iret]);
		break;
	      case VARTOOLS_TYPE_LONG:
		(*((long **) fixval))[k] = (long) atof(argv[*iret]);
		break;
	      case VARTOOLS_TYPE_SHORT:
		(*((short **) fixval))[k] = (short) atof(argv[*iret]);
		break;
	      case VARTOOLS_TYPE_CHAR:
		(*((char **) fixval))[k] = argv[*iret][0];
		break;
	      case VARTOOLS_TYPE_STRING:
		if(((*((char ***) fixval))[j] = (char *) malloc(strlen(argv[*iret])+1)) == NULL)
		  error(ERR_MEMALLOC);
		sprintf((*((char ***) fixval))[k],"%s", argv[*iret]);
		break;
	      default:
		error(ERR_BADTYPE);
		break;
	      }
	    } else {
	      va_end(varlist);
	      return 1;
	    }
	  }
	  *iret = *iret + 1;
	}
      }
      else if(!strcmp(argv[*iret],"expr")) {
	*source = VARTOOLS_SOURCE_EVALEXPRESSION;
	if(Ncolumns <= 0) {
	  *iret = *iret + 1;
	  if(*iret >= argc) { va_end(varlist); return 1; }
	  if((*((char **) expr_char) = (char *) malloc(strlen(argv[*iret])+1)) == NULL)
	    error(ERR_MEMALLOC);
	  sprintf((*((char **) expr_char)),"%s",argv[*iret]);
	} else {
	  if((*((char ***) expr_char) = (char **) malloc(Ncolumns*sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < Ncolumns; j++) {
	    *iret = *iret + 1;
	    if(*iret >= argc) { va_end(varlist); return 1; }
	    if(((*((char ***) expr_char))[j] = (char *) malloc(strlen(argv[*iret])+1)) == NULL)
	      error(ERR_MEMALLOC);
	    sprintf((*((char ***) expr_char))[j],"%s",argv[*iret]);
	  }
	}
	*iret = *iret + 1;
      }
      else if(!strcmp(argv[*iret],"list")) {
	*source = VARTOOLS_SOURCE_INLIST;
	incol = 0;
	*iret = *iret + 1;
	if(*iret < argc) {
	  if(!strcmp(argv[*iret],"column")) {
	    *iret = *iret + 1;
	    if(*iret >= argc) { va_end(varlist); return 1; }
	    incol = atoi(argv[*iret]);
	  } else
	    *iret = *iret - 1;
	} else
	  *iret = *iret - 1;
	RegisterDataFromInputList(p,
				  dataptr, datatype, Ncolumns, cnum, 0, 0,
				  NULL, incol, name);
	*iret = *iret + 1;
      } else if(!strcmp(argv[*iret],"fixcolumn")) {
	*source = VARTOOLS_SOURCE_PRIORCOLUMN;
	if(Ncolumns <= 0) {
	  if(!sizepriornames) {
	    if((priornames = (char **) malloc(sizeof(char *))) == NULL)
	      error(ERR_MEMALLOC);
	    if((priornames[0] = (char *) malloc(MAXLEN)) == NULL)
	      error(ERR_MEMALLOC);
	    sizepriornames = 1;
	  }
	} else {
	  if((*((OutColumn ***) outcolumn) = (OutColumn **) malloc(Ncolumns*sizeof(OutColumn *))) == NULL)
	    error(ERR_MEMALLOC);
	  if(!sizepriornames) {
	    if((priornames = (char **) malloc(Ncolumns*sizeof(char *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=0; k < Ncolumns; k++) {
	      if((priornames[k] = (char *) malloc(MAXLEN)) == NULL)
		error(ERR_MEMALLOC);
	    }
	    sizepriornames = Ncolumns;
	  } 
	  else if(sizepriornames > Ncolumns) {
	    if((priornames = (char **) realloc(priornames, Ncolumns*sizeof(char *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(k=sizepriornames; k < Ncolumns; k++) {
	      if((priornames[k] = (char *) malloc(MAXLEN)) == NULL)
		error(ERR_MEMALLOC);
	    }
	    sizepriornames = Ncolumns;
	  }
	}
	for(k=0; k < (Ncolumns <= 0 ? 1 : Ncolumns); k++) {
	  *iret = *iret + 1;
	  if(*iret < argc) {
	    sprintf(priornames[k],"%s",argv[*iret]);
	  } else {
	    for(k=0; k < sizepriornames; k++)
	      free(priornames[k]);
	    free(priornames);
	    va_end(varlist);
	    return 1;
	  }
	}
	*iret = *iret + 1;
	if(Ncolumns <= 0) {
	  increaselinkedcols(p, (OutColumn **) outcolumn, priornames[0], cnum);
	} else {
	  for(k=0; k < Ncolumns; k++) {
	    increaselinkedcols(p, &((*((OutColumn ***) outcolumn))[k]), priornames[k], cnum);
	  }
	}
      }
      else {
	va_end(varlist);
	return 1;
      }
    }
    else {
      va_end(varlist);
      return 1;
    }
  }

  if(sizepriornames > 0) {
    for(k=0; k < sizepriornames; k++)
      free(priornames[k]);
    free(priornames);
  }
    
  va_end(varlist);
  return 0;
}
    
int ParseParameter(ProgramData *p, Command *c, int *iret, char **argv,
		   int argc, const char *keyword, 
		   int Nvec, _ParseFixSpecFixcolumnStruct *s,
		   _ParseParameter_InitializeStruct *ppstruct)
/* Use this function to parse a command of the form:
   "keyword" <"fix" value | "list" [\"column\" col] | "fixcolumn" <colname | colnum> | "expr" expression>

   p - pointer to the ProgramData struct

   c - pointer to the Command struct associated with this UserCommand

   iret - pointer to a variable storing the current command line index
          number.  On input this should be the term storing "keyword".
          This is updated on output.

   argv - Standard array holding the command-line strings.

   argc - Number of terms on the CL.

   keyword - the keyword to check.

   Nvec - Number of terms to parse

   .... - For each term to parse, the following additional arguments must be
          provided:

         datatype - VARTOOLS_TYPE_DOUBLE, VARTOOLS_TYPE_FLOAT,
                    VARTOOLS_TYPE_INT, VARTOOLS_TYPE_SHORT,
                    VARTOOLS_TYPE_LONG, VARTOOLS_TYPE_STRING,
                    VARTOOLS_TYPE_CHAR

         dataptr - Pointer to the vector or array that will store the input
                   data.

         Ncolumns - Number of columns in the array, if it is <= 0 then
                    dataptr is taken to be a vector rather than an array
                    (e.g. dataptr is a pointer to an variable of type double*
                     rather than double**).

         output - 1 if this data should be included in the output ascii table,
                  0 if it should not be.

         name - root name of the data vector for display in the output
             ascii table and for the input table (if the user gives
             the \"list\" keyword on the command line, and issues the
             -inputlistformat command). The name of the library will
             be added as prefix to the name, and the column number (if
             Ncolumns >= 1) and command numbers will be appended to
             the name.  

	 initialize - 1 if the data vector/array should be created and fixed
                        to a default value when the "keyword" is not given on 
                        the command-line. 
                      0 if the data vector/array should not be created.

         defaultval - the default value to initialize the array to.

    Return values:
       0 - command is successfully parsed.
       1 - "keyword" is not given, command is not parsed.
       2 - "keyword" is given, but the rest of the command is not parsed.
*/
{
  int check, istart, i;
  void *fixptr;
  istart = *iret;
  if(*iret >= argc ? 1 : strcmp(argv[*iret],keyword)) {
    for(i=0; i < Nvec; i++) {
      if(ppstruct[i].doinitialize) {
	fixptr = ppstruct[i].fixptr;
	RegisterDataVector(p, c, s[i].dataptr, s[i].datatype, s[i].Ncolumns,
			   VARTOOLS_SOURCE_FIXED, s[i].output, s[i].name,
			   fixptr);
      }
    }
    return 1;
  }
  (*iret) += 1;
  check = ParseFixSpecFixcolumn(p, c, iret, argv, argc, Nvec, s);
  if(check) {
    *iret = istart;
    return 2;
  }
  else
    return 0;
}

/* This function allocates memory for the registered data arrays.

   We do not need to allocate memory for the data which is read from the
   input list or input lc, that is taken care of by MemAllocDataFromInputList,
   and MemAllocDataFromLightCurve.

   This function does not need to be called by the user.
*/
void MemAllocDataForUserCommand(Command *c, int Nlc)
{
#ifdef DYNAMICLIB
  int i, j, k;
  _UserCommand *co;
  double **dblptr;
  double ***dbl2ptr;
  double ****dbl3ptr;
  int **intptr;
  int ***int2ptr;
  int ****int3ptr;
  short **shortptr;
  short ***short2ptr;
  short ****short3ptr;
  char **charptr;
  char ***char2ptr;
  char ****char3ptr;
  char ***stringptr;
  char ****string2ptr;
  char *****string3ptr;
  float **floatptr;
  float ***float2ptr;
  float ****float3ptr;
  long **longptr;
  long ***long2ptr;
  long ****long3ptr;
  _UserDataPointer *d;
  int Nc;

  co = c->UserCommand;

  for(i=0; i < co->Nptrs; i++) {
    d = &(co->UserDataPointers[i]);
    Nc = d->Ncolumns;
    if((d->source != VARTOOLS_SOURCE_LC && d->source != VARTOOLS_SOURCE_EVALEXPRESSION_LC && d->source != VARTOOLS_SOURCE_EXISTINGVARIABLE) ||
       (d->source == VARTOOLS_SOURCE_EXISTINGVARIABLE ? 
	d->expectedvectortype == VARTOOLS_VECTORTYPE_PERSTARDATA : 0)) {
      if(Nc <= 0) {
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = (double **) d->dataptr;
	  if(((*dblptr) = (double *) malloc(Nlc * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = (char ***) d->dataptr;
	  if(((*stringptr) = (char **) malloc(Nlc * sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < Nlc; j++) {
	    if((((*stringptr)[j]) = (char *) malloc(MAXLEN)) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = (int **) d->dataptr;
	  if(((*intptr) = (int *) malloc(Nlc * sizeof(int))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = (short **) d->dataptr;
	  if(((*shortptr) = (short *) malloc(Nlc * sizeof(short))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = (float **) d->dataptr;
	  if(((*floatptr) = (float *) malloc(Nlc * sizeof(float))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = (long **) d->dataptr;
	  if(((*longptr) = (long *) malloc(Nlc * sizeof(long))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = (char **) d->dataptr;
	  if(((*charptr) = (char *) malloc(Nlc)) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case VARTOOLS_TYPE_USERDEF:
	  if(((*((char **)d->dataptr)) = 
	      (char *) malloc(Nlc * d->size_element)) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < Nlc; j++) {
	    d->initialize_usertype_ptr(Nlc, (void *) (&(*((char **)d->dataptr))[j*d->size_element]), d->extra_user_data);
	  }
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
	  }
	  break;
	case VARTOOLS_TYPE_SHORT:
	  short2ptr = (short ***) d->dataptr;
	  if(((*short2ptr) = (short **) malloc(Nlc * sizeof(short *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < Nlc; j++) {
	    if((((*short2ptr)[j]) = (short *) malloc(Nc * sizeof(short))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  float2ptr = (float ***) d->dataptr;
	  if(((*float2ptr) = (float **) malloc(Nlc * sizeof(float *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < Nlc; j++) {
	    if((((*float2ptr)[j]) = (float *) malloc(Nc * sizeof(float))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  break;
	case VARTOOLS_TYPE_LONG:
	  long2ptr = (long ***) d->dataptr;
	  if(((*long2ptr) = (long **) malloc(Nlc * sizeof(long *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < Nlc; j++) {
	    if((((*long2ptr)[j]) = (long *) malloc(Nc * sizeof(long))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  break;
	case VARTOOLS_TYPE_CHAR:
	  char2ptr = (char ***) d->dataptr;
	  if(((*char2ptr) = (char **) malloc(Nlc * sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < Nlc; j++) {
	    if((((*char2ptr)[j]) = (char *) malloc(Nc)) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      }
    }
  }
#else 
  return;
#endif
}

void RegisterTrackedOpenFile(ProgramData *p, FILE *f){
  if(!p->N_tracked_open_files) {
    if((p->tracked_open_files = (FILE **) malloc(sizeof(FILE *))) == NULL)
      error(ERR_MEMALLOC);
  } else {
    if((p->tracked_open_files = (FILE **) realloc(p->tracked_open_files, (p->N_tracked_open_files+1)*sizeof(FILE *))) == NULL)
      error(ERR_MEMALLOC);
  }
  p->tracked_open_files[p->N_tracked_open_files] = f;
  p->N_tracked_open_files += 1;
  return;
}

void CloseTrackedOpenFiles(ProgramData *p) {
  int i;
  for(i = 0; i < p->N_tracked_open_files; i++) {
    fclose(p->tracked_open_files[i]);
  }
  p->N_tracked_open_files = 0;
  if(p->tracked_open_files != NULL)
    free(p->tracked_open_files);
  p->tracked_open_files = NULL;
  return;
}

void Set_Function_Pointers_Callback(_VARTOOLS_FUNCTION_POINTER_STRUCT *fptr){
#ifdef DYNAMICLIB
  fptr->ParseParameter = &ParseParameter;
  fptr->ParseOutNameKeyword = &ParseOutNameKeyword;
  fptr->ParseConstantParameter = &ParseConstantParameter;
  fptr->ParseFixSpecFixcolumn = &ParseFixSpecFixcolumn;
  fptr->amoeba = &amoeba;
  fptr->mrqmin = &mrqmin;
  fptr->vRegisterDataVector = &vRegisterDataVector;
  fptr->GetOutputFilename = &GetOutputFilename;
  fptr->incrementparameters_foramoeba = &incrementparameters_foramoeba;
  fptr->amoeba_initializesimplexchi2 = &amoeba_initializesimplexchi2;
  fptr->amoeba_cleanup = &amoeba_cleanup;
  fptr->integratemandelagoltransitmodel = &integratemandelagoltransitmodel;
  fptr->mandelagoltransitmodel = &mandelagoltransitmodel;
  fptr->spline = &spline;
  fptr->splint = &splint;
  fptr->spline_monotonic = &spline_monotonic;
  fptr->splint_monotonic = &splint_monotonic;
  fptr->medianfilter = &medianfilter;
  fptr->getweightedmean = &getweightedmean;
  fptr->getmean = &getmean;
  fptr->median = &median;
  fptr->MAD = &MAD;
  fptr->stddev = &stddev;
  fptr->kurtosis = &kurtosis;
  fptr->skewness = &skewness;
  fptr->percentile = &percentile;
  fptr->error = &error;
  fptr->error2 = &error2;
  fptr->fitpoly = &fitpoly;
  fptr->chi2 = &chi2;
  fptr->isDifferentPeriods = &isDifferentPeriods;
  fptr->vsort_generic = &vsort_generic;
  fptr->sortvec_double = &mysort1;
  fptr->vRegisterUserFunction = &vRegisterUserFunction;
  fptr->occultquad = &occultquad;
  fptr->occultnl = &occultnl;
  fptr->memallocdatafromlightcurve = &MemAllocDataFromLightCurve;
  fptr->memallocdatafromlightcurvemidprocess = &MemAllocDataFromLightCurveMidProcess;
  fptr->gnu_getline = &gnu_getline;
  fptr->mysortstringint = &mysortstringint;
  fptr->docorr = &docorr;
  fptr->vAdd_Keyword_To_OutputLC_FitsHeader = &vAdd_Keyword_To_OutputLC_FitsHeader;
  fptr->findX = &findX;
  fptr->findX_string = &findX_string;
  fptr->RegisterTrackedOpenFile = &RegisterTrackedOpenFile;
  fptr->parseone = &parseone;
  fptr->printtostring = &printtostring;
#else
  return;
#endif
}
