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
#include <Python.h>
#include NUMPY_HEADER_FILE
#include <stdio.h>
#include <stdlib.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>

#include "runpython.h"


/* Some functions we will need from the other files */
void SetVariable_Value_Double(int lcindex, int threadindex, int jdindex, _Variable *var, double val);
double EvaluateVariable_Double(int lcindex, int threadindex, int jdindex, _Variable *var);
void MemAllocDataFromLightCurve(ProgramData *p, int threadid, int Nterm);

#define DO_ERROR_MEMALLOC do {fprintf(stderr,"Memory Allocation Error\n"); exit(1);} while(1)


void simpleprinttostring(OutText *text, const char *stoadd)
{
  int lold, j, k, k1;
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
  int lold, j, k, k1, nnewtab;
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



void CreateVartoolsPythonUserFunctionString(ProgramData *p, _PythonCommand *cparent, OutText *text, int cid)
{
  int i, k;
  char tmpstr[MAXLEN];
  char tmpstrtype[MAXLEN];
  char astypecpy[] = ", copy=False)";
  char astypencpy[] = ")";
  char *astypever;
  _PythonCommand *c;

#ifdef HAVE_NUMPY_ASTYPE_COPY
  astypever = astypecpy;
#else
  astypever = astypencpy;
#endif

  if(!cid) c = cparent;
  else {
    c = ((_PythonCommand **) cparent->childcommandptrs)[cid-1];
  }
  /* The function will be named
     _VARTOOLS_PYTHON_USERFUNCTION_%d  where %d is the command number.
     It will take as input a List of variables of the name
     _VARTOOLS_VARIABLES_LIST. */
  simpleprinttostring(text,"def _VARTOOLS_PYTHON_USERFUNCTION_");
  sprintf(tmpstr,"%d",c->cnum);
  simpleprinttostring(text,tmpstr);
  simpleprinttostring(text,"(_VARTOOLS_VARIABLES_LIST):\n");

  k = -1;

  if(c->RequireReadAll) {
    simpleprinttostring(text,"\t_VARTOOLS_INPUT_NLCS = len(_VARTOOLS_VARIABLES_LIST)\n");
  }
  
  /* Set all of the variables to the appropriate values from the List */
  for(i = 0; i < c->Nvars; i++) {
    simpleprinttostring(text,"\t");
    simpleprinttostring(text,c->vars[i]->varname);
    sprintf(tmpstr,"%d",i);
    if(c->RequireReadAll) {
      /* This is an inefficient/slow way of transforming the data, but was the
	 easiest way to code up using the functions written for communicating
	 the single light curve data between vartools and python.  TBD is to
	 write new functions for the case of transfering all of the light curves
	 to python to speed this up so we dont have to do this data
         transposition. */
      if(c->vars[i]->vectortype == VARTOOLS_VECTORTYPE_LC) {
	if(k < 0) k = i;
	simpleprinttostring(text,"=[_VARTOOLS_TMP_ITEM[");
	simpleprinttostring(text,tmpstr);
	simpleprinttostring(text,"]");
	simpleprinttostring(text," for _VARTOOLS_TMP_ITEM in _VARTOOLS_VARIABLES_LIST]\n");
      } else {
	if(c->vars[i]->datatype == VARTOOLS_TYPE_CHAR ||
	   c->vars[i]->datatype == VARTOOLS_TYPE_STRING) {
	  simpleprinttostring(text,"=[_VARTOOLS_TMP_ITEM[");
	  simpleprinttostring(text,tmpstr);
	  simpleprinttostring(text,"][0] for _VARTOOLS_TMP_ITEM in _VARTOOLS_VARIABLES_LIST]\n");
	}
	else {
	  simpleprinttostring(text,"=numpy.core.multiarray.array([_VARTOOLS_TMP_ITEM[");
	  simpleprinttostring(text,tmpstr);
	  simpleprinttostring(text,"][0] for _VARTOOLS_TMP_ITEM in _VARTOOLS_VARIABLES_LIST])\n");
	}
      }
    } else {
      simpleprinttostring(text,"=_VARTOOLS_VARIABLES_LIST[");
      simpleprinttostring(text,tmpstr);
      if(c->vars[i]->vectortype == VARTOOLS_VECTORTYPE_LC) {
	if(k < 0) k = i;
	simpleprinttostring(text,"]\n");
      } else {
	simpleprinttostring(text,"][0]\n");
      }
    }
  }

  /* Copy over the code supplied by the user; note that because we are in a
     function here, add a leading tab */
  simpleprinttostring(text,"\t");
  simpleprinttostring_tabindent(text,c->pythoncommandstring);
  
  if(text->s[text->len_s-1] != '\t')
    simpleprinttostring_tabindent(text,"\n");
  simpleprinttostring(text,"\n");

  simpleprinttostring(text,"\t_VARTOOLS_VARIABLES_OUTPUTLIST = []\n");

  if(!c->RequireReadAll) {

    /* Update the values in the list for all variables */
    for(i=0; i < c->Nvars; i++) {
      if(c->vars[i]->datatype != VARTOOLS_TYPE_STRING &&
	 c->vars[i]->datatype != VARTOOLS_TYPE_CHAR) {
	switch(c->vars[i]->datatype) {
	case VARTOOLS_TYPE_CONVERTJD:
	case VARTOOLS_TYPE_DOUBLE:
	case VARTOOLS_TYPE_FLOAT:
	  sprintf(tmpstrtype,".astype(float%s",astypever);
	  break;
	case VARTOOLS_TYPE_INT:
	case VARTOOLS_TYPE_LONG:
	case VARTOOLS_TYPE_SHORT:
	  sprintf(tmpstrtype,".astype(int%s",astypever);
	  break;
	default:
	  tmpstrtype[0] = '\0';
	  break;
	}
	simpleprinttostring(text,"\tif '");
	simpleprinttostring(text,c->vars[i]->varname);
	simpleprinttostring(text,"' in locals():\n");
	simpleprinttostring(text,"\t\tif type(");
	simpleprinttostring(text,c->vars[i]->varname);
	simpleprinttostring(text,") == type(numpy.core.multiarray.zeros(0)):\n");
	simpleprinttostring(text,"\t\t\t_VARTOOLS_VARIABLES_OUTPUTLIST.append(");
	simpleprinttostring(text,c->vars[i]->varname);
	simpleprinttostring(text,tmpstrtype);
	simpleprinttostring(text,")\n");
	simpleprinttostring(text,"\t\telse:\n");
	simpleprinttostring(text,"\t\t\t_VARTOOLS_VARIABLES_OUTPUTLIST.append(numpy.core.multiarray.array([");
	simpleprinttostring(text,c->vars[i]->varname);
	simpleprinttostring(text,"])");
	simpleprinttostring(text,tmpstrtype);
	simpleprinttostring(text,")\n");
	simpleprinttostring(text,"\telse:\n");
	if(c->vars[i]->vectortype != VARTOOLS_VECTORTYPE_LC) {
	  simpleprinttostring(text,"\t\t_VARTOOLS_VARIABLES_OUTPUTLIST.append(numpy.core.multiarray.zeros(1)");
	  simpleprinttostring(text,tmpstrtype);
	  simpleprinttostring(text,")\n");
	} else {
	  simpleprinttostring(text,"\t\t_VARTOOLS_VARIABLES_OUTPUTLIST.append(numpy.core.multiarray.zeros(1)");
	  simpleprinttostring(text,tmpstrtype);
	  simpleprinttostring(text,")\n");
	}
      } else {
	simpleprinttostring(text,"\tif '");
	simpleprinttostring(text,c->vars[i]->varname);
	simpleprinttostring(text,"' in locals():\n");
	simpleprinttostring(text,"\t\tif type(");
	simpleprinttostring(text,c->vars[i]->varname);
	simpleprinttostring(text,") is list:\n");
	simpleprinttostring(text,"\t\t\t_VARTOOLS_VARIABLES_OUTPUTLIST.append(");
	simpleprinttostring(text,c->vars[i]->varname);
	simpleprinttostring(text,")\n");
	simpleprinttostring(text,"\t\telse:\n");
	simpleprinttostring(text,"\t\t\t_VARTOOLS_VARIABLES_OUTPUTLIST.append([");
	simpleprinttostring(text,c->vars[i]->varname);
	simpleprinttostring(text,"])\n");
	simpleprinttostring(text,"\telse:\n");
	simpleprinttostring(text,"\t\t_VARTOOLS_VARIABLES_OUTPUTLIST.append([\"0\"])\n");
      }
    }


    /* Append any variables that are created by this process to the list; if
       they are not defined, set them to zero */
    for(i=0; i < c->Nvars_outonly; i++) {
      if(c->outonlyvars[i]->datatype != VARTOOLS_TYPE_CHAR &&
	 c->outonlyvars[i]->datatype != VARTOOLS_TYPE_STRING) {
	switch(c->outonlyvars[i]->datatype) {
	case VARTOOLS_TYPE_CONVERTJD:
	case VARTOOLS_TYPE_DOUBLE:
	case VARTOOLS_TYPE_FLOAT:
	  sprintf(tmpstrtype,".astype(float%s",astypever);
	  break;
	case VARTOOLS_TYPE_INT:
	case VARTOOLS_TYPE_LONG:
	case VARTOOLS_TYPE_SHORT:
	  sprintf(tmpstrtype,".astype(int%s",astypever);
	  break;
	case VARTOOLS_TYPE_STRING:
	case VARTOOLS_TYPE_CHAR:
	  sprintf(tmpstrtype,".astype(string%s",astypever);
	  break;
	default:
	  tmpstrtype[0] = '\0';
	  break;
	}
	simpleprinttostring(text,"\tif '");
	simpleprinttostring(text,c->outonlyvars[i]->varname);
	simpleprinttostring(text,"' in locals():\n");
	simpleprinttostring(text,"\t\tif type(");
	simpleprinttostring(text,c->outonlyvars[i]->varname);
	simpleprinttostring(text,") == type(numpy.core.multiarray.zeros(0)):\n");
	simpleprinttostring(text,"\t\t\t_VARTOOLS_VARIABLES_OUTPUTLIST.append(");
	simpleprinttostring(text,c->outonlyvars[i]->varname);
	simpleprinttostring(text,tmpstrtype);
	simpleprinttostring(text,")\n");
	simpleprinttostring(text,"\t\telse:\n");
	simpleprinttostring(text,"\t\t\t_VARTOOLS_VARIABLES_OUTPUTLIST.append(numpy.core.multiarray.array([");
	simpleprinttostring(text,c->outonlyvars[i]->varname);
	simpleprinttostring(text,"])");
	simpleprinttostring(text,tmpstrtype);
	simpleprinttostring(text,")\n");
	simpleprinttostring(text,"\telse:\n");
	if(c->outonlyvars[i]->vectortype != VARTOOLS_VECTORTYPE_LC) {
	  simpleprinttostring(text,"\t\t_VARTOOLS_VARIABLES_OUTPUTLIST.append(numpy.core.multiarray.zeros(1)");
	  simpleprinttostring(text,tmpstrtype);
	  simpleprinttostring(text,")\n");
	} else {
	  simpleprinttostring(text,"\t\t_VARTOOLS_VARIABLES_OUTPUTLIST.append(numpy.core.multiarray.zeros(1)");
	  simpleprinttostring(text,tmpstrtype);
	  simpleprinttostring(text,")\n");
	}
      } else {
	simpleprinttostring(text,"\tif '");
	simpleprinttostring(text,c->outonlyvars[i]->varname);
	simpleprinttostring(text,"' in locals():\n");
	simpleprinttostring(text,"\t\tif type(");
	simpleprinttostring(text,c->outonlyvars[i]->varname);
	simpleprinttostring(text,") is list:\n");
	simpleprinttostring(text,"\t\t\t_VARTOOLS_VARIABLES_OUTPUTLIST.append(");
	simpleprinttostring(text,c->outonlyvars[i]->varname);
	simpleprinttostring(text,")\n");
	simpleprinttostring(text,"\t\telse:\n");
	simpleprinttostring(text,"\t\t\t_VARTOOLS_VARIABLES_OUTPUTLIST.append([");
	simpleprinttostring(text,c->outonlyvars[i]->varname);
	simpleprinttostring(text,"])\n");
	simpleprinttostring(text,"\telse:\n");
	simpleprinttostring(text,"\t\t_VARTOOLS_VARIABLES_OUTPUTLIST.append([\"0\"])\n");
      }
    }
  } else {
    /* Again for the case of transferring all of the data to python this is
       inefficient but was easiest to code up quickly. No error checking is
       performed here, anticipating that this will need to be rewritten when
       a better method for transferring the data is implemented. */
    simpleprinttostring(text,"\tfor _VARTOOLS_TMP_ITEM in range(0,_VARTOOLS_INPUT_NLCS):\n");
    simpleprinttostring(text,"\t\t_VARTOOLS_VARIABLES_TMPOUTLIST = []\n");
    for(i=0; i < c->Nvars; i++) {
      if(c->vars[i]->vectortype == VARTOOLS_VECTORTYPE_LC) {
	simpleprinttostring(text,"\t\t_VARTOOLS_VARIABLES_TMPOUTLIST.append(");
	simpleprinttostring(text,c->vars[i]->varname);
	simpleprinttostring(text,"[_VARTOOLS_TMP_ITEM])\n");
      } else {
	if(c->vars[i]->datatype != VARTOOLS_TYPE_STRING &&
	   c->vars[i]->datatype != VARTOOLS_TYPE_CHAR) {
	  simpleprinttostring(text,"\t\t_VARTOOLS_VARIABLES_TMPOUTLIST.append(numpy.core.multiarray.array([");
	  simpleprinttostring(text,c->vars[i]->varname);
	  simpleprinttostring(text,"[_VARTOOLS_TMP_ITEM]]))\n");
	} else {
	  simpleprinttostring(text,"\t\t_VARTOOLS_VARIABLES_TMPOUTLIST.append([");
	  simpleprinttostring(text,c->vars[i]->varname);
	  simpleprinttostring(text,"[_VARTOOLS_TMP_ITEM]])\n");
	}
      }
    }
    for(i=0; i < c->Nvars_outonly; i++) {
      if(c->outonlyvars[i]->vectortype == VARTOOLS_VECTORTYPE_LC) {
	simpleprinttostring(text,"\t\t_VARTOOLS_VARIABLES_TMPOUTLIST.append(");
	simpleprinttostring(text,c->outonlyvars[i]->varname);
	simpleprinttostring(text,"[_VARTOOLS_TMP_ITEM])\n");
      } else {
	if(c->outonlyvars[i]->datatype != VARTOOLS_TYPE_STRING &&
	   c->outonlyvars[i]->datatype != VARTOOLS_TYPE_CHAR) {
	  simpleprinttostring(text,"\t\t_VARTOOLS_VARIABLES_TMPOUTLIST.append(numpy.core.multiarray.array([");
	  simpleprinttostring(text,c->outonlyvars[i]->varname);
	  simpleprinttostring(text,"[_VARTOOLS_TMP_ITEM]]))\n");
	} else {
	  simpleprinttostring(text,"\t\t_VARTOOLS_VARIABLES_TMPOUTLIST.append([");
	  simpleprinttostring(text,c->outonlyvars[i]->varname);
	  simpleprinttostring(text,"[_VARTOOLS_TMP_ITEM]])\n");
	}
      }
    }
    simpleprinttostring(text,"\t\t_VARTOOLS_VARIABLES_OUTPUTLIST.append(_VARTOOLS_VARIABLES_TMPOUTLIST)\n");
  }
  //simpleprinttostring(text,"\tprint _VARTOOLS_VARIABLES_OUTPUTLIST\n");
  simpleprinttostring(text,"\treturn _VARTOOLS_VARIABLES_OUTPUTLIST\n\n");  

  //fprintf(stderr,"%s", text->s);
}

void CleanUpPythonObjectContainerVariables(_PythonObjectContainer *py) {
  int i;
  if(py->Variables == NULL) {
    if(py->VariableList != NULL) {
      if(py->VariableList->ob_refcnt > 0)
	Py_CLEAR(py->VariableList);
      py->VariableList = NULL;
    }
    if(py->data != NULL) {
      for(i=0; i < py->Nvars; i++) {
	if(py->data[i].dataptr != NULL)
	  free(py->data[i].dataptr);
      }
      free(py->data);
    }
    py->data = NULL;
    py->Nvars = 0;
    return;
  }
  if(py->VariableList != NULL) {
    if(py->VariableList->ob_refcnt > 0)
      Py_CLEAR(py->VariableList);
    py->VariableList = NULL;
  }
  for(i=0; i < py->Nvars; i++) {
    if(py->Variables[i] != NULL) {
      if(py->Variables[i]->ob_refcnt > 0)
	Py_CLEAR(py->Variables[i]);
    }
    py->Variables[i] = NULL;
  }
  if(py->data != NULL) {
    for(i=0; i < py->Nvars; i++) {
      if(py->data[i].dataptr != NULL)
	free(py->data[i].dataptr);
    }
    free(py->data);
  }
  py->data = NULL;
  if(py->Variables != NULL)
    free(py->Variables);
  py->Variables = NULL;
  py->Nvars = 0;
}

void CleanUpPythonObjectContainer(_PythonObjectContainer **py) {
  if(py == NULL) return;
  if(*py == NULL) return;
  CleanUpPythonObjectContainerVariables(*py);
  free((*py));
  *py = NULL;
}

_PythonObjectContainer *CreatePythonObjectContainer(void) {
  _PythonObjectContainer *py;
  if((py = (_PythonObjectContainer *) malloc(sizeof(_PythonObjectContainer))) == NULL)
    DO_ERROR_MEMALLOC;
  py->CompiledUserCode = NULL;
  py->UserModule = NULL;
  py->UserFunctionToRun = NULL;
  py->Variables = NULL;
  py->VariableList = NULL;
  py->Nvars = 0;
  py->data = NULL;
  py->Nfunc = 0;
  py->FullList = NULL;
  py->FullListOut = NULL;
  return py;
}



int InitializePython(ProgramData *p, _PythonCommand *c, int threadindex)
{
  OutText usercodetext;
  _PythonObjectContainer *py;
  char tmpstr[MAXLEN];
  int i, j, testneedreadall;
#ifdef HAVE_PYTHON_UNICODE
  wchar_t *tmp1, *tmp2, *tmp3;
#endif
  
  usercodetext.s = NULL;
  usercodetext.space = 0;
  usercodetext.len_s = 0;
  usercodetext.Nchar_cur_line = 0;

#ifdef PYTHONHOMEPATH
#ifdef HAVE_PYTHON_UNICODE
  tmp1 = Py_DecodeLocale(PYTHONHOMEPATH,NULL);
  Py_SetPythonHome(tmp1);
  tmp3 = Py_DecodeLocale(PYTHONSEARCHPATH,NULL);
  Py_SetPath(tmp3);
#else
  Py_SetPythonHome(PYTHONHOMEPATH);
  Py_SetPath(PYTHONSEARCHPATH);
#endif
#endif

#ifdef HAVE_PYTHON_UNICODE
  tmp1 = Py_DecodeLocale(c->progname,NULL);
  Py_SetProgramName(tmp1);
#else
  Py_SetProgramName(c->progname);
#endif
  Py_Initialize();
  import_array1(1);

  simpleprinttostring(&usercodetext,"import numpy\n");

  if(c->len_pythoninitializationtextstring > 0) {
    if(c->pythoninitializationtext[0] != '\0') {
      simpleprinttostring(&usercodetext, c->pythoninitializationtext);
    }
  }
  simpleprinttostring(&usercodetext,"\n");
	
  py = CreatePythonObjectContainer();


  ((_PythonObjectContainer **) c->pythonobjects)[threadindex] = py;

  testneedreadall = 0;
  if(c->RequireReadAll) {
    testneedreadall = 1;
  }
  if(!testneedreadall && c->Nchildren > 0) {
    for(j=0; j < c->Nchildren; j++) {
      if(((_PythonCommand **)(c->childcommandptrs))[j]->RequireReadAll) {
	testneedreadall = 1;
	break;
      }
    }
  }
  if(testneedreadall) {
    for(j=1; j < p->Nlcs; j++) {
      ((_PythonObjectContainer **) c->pythonobjects)[j] = CreatePythonObjectContainer();
    }
  }

  if((py->UserFunctionToRun = (PyObject **) malloc((c->Nchildren + 1)*sizeof(PyObject *))) == NULL) {
    DO_ERROR_MEMALLOC;
  }
  py->Nfunc = c->Nchildren + 1;

  /* Create the functions to run the python code supplied by the user */
  for(i=0; i < c->Nchildren + 1; i++) {
    CreateVartoolsPythonUserFunctionString(p, c, &usercodetext, i);
  }

  /* Compile all of the code supplied by the user */
  py->CompiledUserCode = Py_CompileString(usercodetext.s,"", Py_file_input);
  
  if(py->CompiledUserCode == NULL) {
    fprintf(stderr,"Error compiling the following python code from vartools:\n\n%s\n", usercodetext.s);
    if(PyErr_Occurred() != NULL) {
      fprintf(stderr,"The error message from the python compiler is as follows:\n\n");
      PyErr_Print();
    }
    else {
      fprintf(stderr,"No error message from python\n\n");
    }
    exit(1);
  }

  /* Import this code as a module */
  py->UserModule = PyImport_ExecCodeModule( "_VARTOOLS_PYTHON_USERFUNCTION_MODULE",
					    py->CompiledUserCode);
    
  if(py->UserModule == NULL) {
    fprintf(stderr,"Error creating a module to run the following python code from vartools:\n\n%s\n", usercodetext.s);
    if(PyErr_Occurred() != NULL) {
      fprintf(stderr,"The error message from python is as follows:\n\n");
      PyErr_Print();
    }
    else {
      fprintf(stderr,"No error message from python\n\n");
    }
    exit(1);
  }
  
  /* Get pointers to the executable user function code from the module */
  for(i=0; i < c->Nchildren + 1; i++) {
    sprintf(tmpstr,"_VARTOOLS_PYTHON_USERFUNCTION_");
    if(!i) {
      sprintf(tmpstr,"%s%d",tmpstr,c->cnum);
    } else {
      sprintf(tmpstr,"%s%d",tmpstr,((_PythonCommand **)c->childcommandptrs)[i-1]->cnum);
    }
    
    py->UserFunctionToRun[i] = PyObject_GetAttrString( py->UserModule, tmpstr);
    if(py->UserFunctionToRun[i] == NULL) {
      fprintf(stderr,"Error loading the following python code from vartools:\n\n%s\n", usercodetext.s);
      if(PyErr_Occurred() != NULL) {
	fprintf(stderr,"The error message from python is as follows:\n\n");
	PyErr_Print();
      }
      else {
	fprintf(stderr,"No error message from python\n\n");
      }
      exit(1);
    }
  }

  return 0;

  /* Everything should be setup now to call this function through the
     RunPythonProcessingLoop */
}


#define _EXIT_READ_VARIABLES do { if(tmpinpstr != NULL) free(tmpinpstr); if(tmpinpstr2 != NULL) free(tmpinpstr2); return 1; } while(0)

int ReadVariablesFromSocketIntoPython(ProgramData *p, _PythonCommand *c, 
				      int threadindex, int isall, int lcnum) {
  _PythonObjectContainer *py = NULL;

  PyObject *tmpobj = NULL;

  int Nvars, i, j, lenvec, k;
  char datatype;
  npy_intp dims[1];
  size_t databytesize;
  char *tmpinpstr = NULL;
  int sizetmpinpstr = 0;
  char *tmpinpstr2 = NULL;
  int sizetmpinpstr2 = 0;

  int datatype_npy;

  size_t sizedataread;
  size_t sizedatareadtot;


  if(!isall) {
    py = ((_PythonObjectContainer **) c->pythonobjects)[threadindex];

    if(py->Variables != NULL)
      CleanUpPythonObjectContainerVariables(py);
  } else {
    py = ((_PythonObjectContainer **) c->pythonobjects)[lcnum];

    if(py->Variables != NULL)
      CleanUpPythonObjectContainerVariables(py);
  }

  if(read(c->sockets[threadindex][1], &Nvars, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Error reading the number of variables to read in the parent in the ReadVariablesFromSocketIntoPython function.\n");
    _EXIT_READ_VARIABLES;
  }


  /* The only argument is a list of input variables; First create the list */
  py->VariableList = PyList_New(Nvars);
  Py_INCREF(py->VariableList);
  
  py->Nvars = Nvars;
  if(Nvars >= 0) {
    if((py->Variables = (PyObject **) malloc(Nvars*sizeof(PyObject *))) == NULL) {
      fprintf(stderr,"Memory Allocation Error in VARTOOLS Python sub-process.\n");
      _EXIT_READ_VARIABLES;
    }
    if((py->data = (_PythonArrayData *) malloc(Nvars*sizeof(_PythonArrayData))) == NULL) {
      fprintf(stderr,"Memory Allocation Error in VARTOOLS Python sub-process.\n");
      _EXIT_READ_VARIABLES;
    }
    for(i=0; i < Nvars; i++) py->Variables[i] = NULL;
    for(i=0; i < Nvars; i++) py->data[i].dataptr = NULL;
    for(i=0; i < Nvars; i++) {
      if(read(c->sockets[threadindex][1], &datatype, sizeof(char)) < sizeof(char)) {
	fprintf(stderr,"Error reading datatype for expected variable index %d in the ReadVariablesFromSocketIntoPython function.\n", i);
	_EXIT_READ_VARIABLES;
      }
      if(datatype != VARTOOLS_TYPE_STRING && datatype != VARTOOLS_TYPE_CHAR) {
	switch(datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  databytesize = sizeof(double);
	  datatype_npy = NPY_DOUBLE;
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  databytesize = sizeof(double);
	  datatype_npy = NPY_DOUBLE;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  databytesize = sizeof(float);
	  datatype_npy = NPY_FLOAT;
	  break;
	case VARTOOLS_TYPE_INT:
	  databytesize = sizeof(int);
	  datatype_npy = NPY_INT;
	  break;
	case VARTOOLS_TYPE_LONG:
	  databytesize = sizeof(long);
	  datatype_npy = NPY_LONG;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  databytesize = sizeof(short);
	  datatype_npy = NPY_SHORT;
	  break;
	default:
	  fprintf(stderr,"Error: invalid datatype received for variable index %d in the ReadVariablesFromSocketIntoPython function\n", i);
	  _EXIT_READ_VARIABLES;
	}
	if(read(c->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	  fprintf(stderr,"Error reading length of vector for expected variable index %d in the ReadVariablesFromSocketIntoPython function.\n", i);
	  _EXIT_READ_VARIABLES;
	}
	py->data[i].datatype = datatype;
	if(lenvec <= 0) {
	  continue;
	}
	if((py->data[i].dataptr = (void *) malloc(((size_t) lenvec)*databytesize)) == NULL) {
	  fprintf(stderr,"Memory Allocation Error in VARTOOLS Python sub-process.\n");
	  _EXIT_READ_VARIABLES;
	}
	sizedatareadtot = 0;
	while(sizedatareadtot < (((size_t) lenvec)*databytesize)) {
	  sizedataread = read(c->sockets[threadindex][1], (void *) (((char *) (py->data[i].dataptr)) + sizedatareadtot), (((size_t) lenvec)*databytesize) - sizedatareadtot);
	  sizedatareadtot += sizedataread;
	}
	if(sizedatareadtot < (((size_t) lenvec)*databytesize)) {
	  fprintf(stderr,"Error reading data for variable index %d in the ReadVariablesFromSocketIntoPython function.\n", i);
	  _EXIT_READ_VARIABLES;
	}
	dims[0] = lenvec;
	py->Variables[i] = PyArray_SimpleNewFromData(1, &(dims[0]), datatype_npy, py->data[i].dataptr);
	Py_INCREF(py->Variables[i]);
	if(py->Variables[i] == NULL) {
	  fprintf(stderr,"Error creating numpy array variable for variable index %d in the ReadVariablesFromSocketIntoPython function.\n", i);
	  _EXIT_READ_VARIABLES;
	}
      }
      else if(datatype == VARTOOLS_TYPE_CHAR) {
	if(read(c->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	  fprintf(stderr,"Error reading length of vector for expected variable index %d in the ReadVariablesFromSocketIntoPython function.\n", i);
	  _EXIT_READ_VARIABLES;
	}
	if(lenvec <= 0) {
	  continue;
	}
	if(lenvec > sizetmpinpstr) {
	  if(!sizetmpinpstr) {
	    if((tmpinpstr = (char *) malloc(lenvec * sizeof(char))) == NULL) {
	      fprintf(stderr,"Memory Allocation Error in VARTOOLS Python sub-process.\n");
	      _EXIT_READ_VARIABLES;
	    }
	    sizetmpinpstr = lenvec;
	  } else {
	    if((tmpinpstr = (char *) realloc(tmpinpstr, lenvec * sizeof(char))) == NULL) {
	      fprintf(stderr,"Memory Allocation Error in VARTOOLS Python sub-process.\n");
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
	  fprintf(stderr,"Error reading data for variable index %d in the ReadVariablesFromSocketIntoPython function.\n", i);
	  _EXIT_READ_VARIABLES;
	}
	if(sizetmpinpstr2 < 2) {
	  if(!sizetmpinpstr2) {
	    if((tmpinpstr2 = (char *) malloc(2*sizeof(char))) == NULL) {
	      fprintf(stderr,"Memory Allocation Error in VARTOOLS Python sub-process.\n");
	      _EXIT_READ_VARIABLES;
	    }
	  } else {
	    if((tmpinpstr2 = (char *) realloc(tmpinpstr2,2*sizeof(char))) == NULL) {
	      fprintf(stderr,"Memory Allocation Error in VARTOOLS Python sub-process.\n");
	      _EXIT_READ_VARIABLES;
	    }
	  }
	}
	py->Variables[i] = PyList_New(lenvec);
	Py_INCREF(py->Variables[i]);
	if(py->Variables[i] == NULL) {
	  fprintf(stderr,"Error creating list variable for variable index %d in the ReadVariablesFromSocketIntoPython function.\n", i);
	  _EXIT_READ_VARIABLES;
	}
	for(j=0; j < lenvec; j++) {
	  tmpinpstr2[0] = tmpinpstr[j];
	  tmpinpstr2[1] = '\0';
#ifdef HAVE_PYTHON_UNICODE
	  tmpobj = PyUnicode_FromString(tmpinpstr2);
#else
	  tmpobj = PyString_FromString(tmpinpstr2);
#endif
	  if(PyList_SetItem(py->Variables[i],j,tmpobj)) {
	    fprintf(stderr,"Error adding a string variable to the list for variable index %d in the ReadVariablesFromSocketIntoPython function.\n", i);
	    _EXIT_READ_VARIABLES;
	  }
	}
      } else if(datatype == VARTOOLS_TYPE_STRING) {
	if(read(c->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	  fprintf(stderr,"Error reading length of vector for expected variable index %d in the ReadVariablesFromSocketIntoPython function.\n", i);
	  _EXIT_READ_VARIABLES;
	}
	if(lenvec <= 0) {
	  continue;
	}
	py->Variables[i] = PyList_New(lenvec);
	if(py->Variables[i] == NULL) {
	  fprintf(stderr,"Error creating list variable for variable index %d in the ReadVariablesFromSocketIntoPython function.\n", i);
	  _EXIT_READ_VARIABLES;
	}
	for(j=0; j < lenvec; j++) {
	  if(read(c->sockets[threadindex][1], &k, sizeof(int)) < sizeof(int)) {
	    fprintf(stderr,"Error reading length of string for vector %d, item %d in the ReadVariablesFromSocketIntoPython function.\n", i, j);
	    _EXIT_READ_VARIABLES;
	  }
	  if(k < 0) k = 0;
	  if(k + 1 > sizetmpinpstr) {
	    if(!sizetmpinpstr) {
	      if((tmpinpstr = (char *) malloc((k+1) * sizeof(char))) == NULL) {
		fprintf(stderr,"Memory Allocation Error in VARTOOLS Python sub-process.\n");
		_EXIT_READ_VARIABLES;
	      }
	      sizetmpinpstr = k+1;
	    } else {
	      if((tmpinpstr = (char *) realloc(tmpinpstr, (k+1) * sizeof(char))) == NULL) {
		fprintf(stderr,"Memory Allocation Error in VARTOOLS Python sub-process.\n");
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
	      fprintf(stderr,"Error reading length of string for vector %d, item %d in the ReadVariablesFromSocketIntoPython function.\n", i, j);
	      _EXIT_READ_VARIABLES;
	    }
	  }
	  tmpinpstr[k] = '\0';
#ifdef HAVE_PYTHON_UNICODE
	  tmpobj = PyUnicode_FromString(tmpinpstr);
#else
	  tmpobj = PyString_FromString(tmpinpstr);
#endif
	  if(PyList_SetItem(py->Variables[i],j,tmpobj)) {
	    fprintf(stderr,"Error adding a string variable to the list for variable index %d in the ReadVariablesFromSocketIntoPython function.\n", i);
	    _EXIT_READ_VARIABLES;
	  }
	}
      }
      /* Transfer the variable to the list */
      if(PyList_SetItem(py->VariableList, i, py->Variables[i])) {
	fprintf(stderr,"Error running python command.\n");
	if(PyErr_Occurred() != NULL) {
	  fprintf(stderr,"The error message from python is as follows:\n\n");
	  PyErr_Print();
	} else {
	  fprintf(stderr,"No error reported by python");
	}
	_EXIT_READ_VARIABLES;
      }
    }
  }
  if(tmpinpstr != NULL)
    free(tmpinpstr);
  if(tmpinpstr2 != NULL)
    free(tmpinpstr2);
  return 0;
}

#undef _EXIT_READ_VARIABLES

int ReadVariablesFromSocketIntoPython_all_lcs(ProgramData *p, 
					      _PythonCommand *c, 
					      int threadindex) {
  int Nlcs, i;

  PyObject *FullList;

  _PythonObjectContainer *py = NULL;

  if(read(c->sockets[threadindex][1], &Nlcs, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Error reading the number of light curves in the ReadVariablesFromSocketIntoPython_all_lcs function.\n");
    return 1;
  }

  FullList = PyList_New(Nlcs);  
  Py_INCREF(FullList);

  for(i=0; i < Nlcs; i++) {
    if(ReadVariablesFromSocketIntoPython(p, c, threadindex, 1, i)) {
      return 1;
    }
    py = ((_PythonObjectContainer **) c->pythonobjects)[i];
    PyList_SetItem(FullList, i, py->VariableList);
  }
  py = ((_PythonObjectContainer **) c->pythonobjects)[threadindex];
  py->FullList = FullList;
  return 0;

}

int WriteVariablesFromPythonToSocket(ProgramData *p, _PythonCommand *cparent, 
				     int threadindex, int cid) {
  _PythonObjectContainer *py = NULL;
  int i, j, k, ii;
  PyObject *tmparray = NULL;
  PyObject *tmplist = NULL;
  PyObject *tmpcopyarray = NULL;
  PyArray_Descr *descr;
  npy_intp *outdims;

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

  npy_intp tmpindx = 0;

  PyObject *tmpobj = NULL;

  _PythonCommand *c;

  if(!cid) c = cparent;
  else c = ((_PythonCommand **)cparent->childcommandptrs)[cid-1];

  py = ((_PythonObjectContainer **) cparent->pythonobjects)[threadindex];

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
      tmparray = PyList_GetItem(py->VariableListOut,ii);
      if(tmparray == NULL) {
	fprintf(stderr,"Error in WriteVariablesFromPythonToSocket. The output variable list from python is shorter than expected.\n");
	if(tmpcharvec != NULL) free(tmpcharvec);
	return 1;
      }
      if(PyArray_NDIM(tmparray) != 1) {
	fprintf(stderr,"Error the variable %s has the wrong number of array dimensions upon output from python. It should either be a scalar or a one dimensional array, but instead has dimension %d.\n", v->varname, (int) PyArray_NDIM(tmparray));
	if(tmpcharvec != NULL) free(tmpcharvec);
	return 1;
      }
      outdims = PyArray_DIMS(tmparray);
      if(v->vectortype != VARTOOLS_VECTORTYPE_LC) {
	if(outdims[0] != 1) {
	  fprintf(stderr,"Error the variable %s has the wrong dimension upon output from python. This variable is expected to be a scalar, but was instead found to have dimension %d\n",v->varname,(int) (outdims[0]));
	  if(tmpcharvec != NULL) free(tmpcharvec);
	  return 1;
	}
	switch(v->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	case VARTOOLS_TYPE_CONVERTJD:
	  lenvec = 1;
	  if(write(cparent->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    return 1;
	  }
	  tmpobj = PyArray_GETITEM((PyArrayObject *) tmparray,PyArray_GetPtr((PyArrayObject *) tmparray, &tmpindx));
	  tmpdblout = PyFloat_AS_DOUBLE(tmpobj);
	  if(write(cparent->sockets[threadindex][1], &tmpdblout, sizeof(double)) < sizeof(double)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    return 1;
	  }
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  lenvec = 1;
	  if(write(cparent->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    return 1;
	  }
	  tmpobj = PyArray_GETITEM((PyArrayObject *) tmparray,PyArray_GetPtr((PyArrayObject *) tmparray, &tmpindx));
	  tmpfloatout = (float) PyFloat_AS_DOUBLE(tmpobj);
	  if(write(cparent->sockets[threadindex][1], &tmpfloatout, sizeof(float)) < sizeof(float)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    return 1;
	  }
	  break;
	case VARTOOLS_TYPE_INT:
	  lenvec = 1;
	  if(write(cparent->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    return 1;
	  }
	  tmpobj = PyArray_GETITEM((PyArrayObject *) tmparray,PyArray_GetPtr((PyArrayObject *) tmparray, &tmpindx));
#ifdef HAVE_PYTHON_UNICODE
	  tmpintout = (int) PyLong_AS_LONG(tmpobj);
#else
	  tmpintout = (int) PyInt_AS_LONG(tmpobj);
#endif
	  if(write(cparent->sockets[threadindex][1], &tmpintout, sizeof(int)) < sizeof(int)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    return 1;
	  }
	  break;
	case VARTOOLS_TYPE_LONG:
	  lenvec = 1;
	  if(write(cparent->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    return 1;
	  }
	  tmpobj = PyArray_GETITEM((PyArrayObject *) tmparray,PyArray_GetPtr((PyArrayObject *) tmparray, &tmpindx));
	  tmplongout = (long) PyLong_AsLong(tmpobj);
	  if(write(cparent->sockets[threadindex][1], &tmplongout, sizeof(long)) < sizeof(long)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    return 1;
	  }
	  break;
	case VARTOOLS_TYPE_SHORT:
	  lenvec = 1;
	  if(write(cparent->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    return 1;
	  }
	  tmpobj = PyArray_GETITEM((PyArrayObject *) tmparray,PyArray_GetPtr((PyArrayObject *) tmparray, &tmpindx));
#ifdef HAVE_PYTHON_UNICODE
	  tmpshortout = (short) PyLong_AS_LONG(tmpobj);
#else
	  tmpshortout = (short) PyInt_AS_LONG(tmpobj);
#endif
	  if(write(cparent->sockets[threadindex][1], &tmpshortout, sizeof(short)) < sizeof(short)) {
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    return 1;
	  }
	  break;
	}
      } else {
	lenvec = outdims[0];
	if(write(cparent->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	  if(tmpcharvec != NULL) free(tmpcharvec);
	  return 1;
	}
	if(lenvec <= 0) continue;
	descr = PyArray_DESCR(tmparray);
	switch(v->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	case VARTOOLS_TYPE_CONVERTJD:
	  if(descr->type_num != NPY_DOUBLE) {
	    tmpcopyarray = PyArray_SimpleNew(1, outdims, NPY_DOUBLE);
	    if(PyArray_CopyInto((PyArrayObject *) tmpcopyarray, (PyArrayObject *) tmparray)) {
	      fprintf(stderr,"Error converting vector %s to the expected double datatype in output from python.\n",v->varname);
	      Py_CLEAR(tmpcopyarray);
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      return 1;
	    }
	    dblptrout = (double *) PyArray_DATA(tmpcopyarray);
	    if(write(cparent->sockets[threadindex][1], (void *) dblptrout, (((size_t) lenvec)*(sizeof(double)))) < (((size_t) lenvec)*(sizeof(double)))) {
	      Py_CLEAR(tmpcopyarray);
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      return 1;
	    }
	    Py_CLEAR(tmpcopyarray);
	  } else {
	    dblptrout = (double *) PyArray_DATA(tmparray);
	    if(write(cparent->sockets[threadindex][1], (void *) dblptrout, (((size_t) lenvec)*(sizeof(double)))) < (((size_t) lenvec)*(sizeof(double)))) {
	      Py_CLEAR(tmpcopyarray);
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      return 1;
	    }
	  }
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  if(descr->type_num != NPY_FLOAT) {
	    tmpcopyarray = PyArray_SimpleNew(1, outdims, NPY_FLOAT);
	    if(PyArray_CopyInto((PyArrayObject *) tmpcopyarray, (PyArrayObject *) tmparray)) {
	      fprintf(stderr,"Error converting vector %s to the expected float datatype in output from python.\n",v->varname);
	      Py_CLEAR(tmpcopyarray);
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      return 1;
	    }
	    floatptrout = (float *) PyArray_DATA(tmpcopyarray);
 	    if(write(cparent->sockets[threadindex][1], (void *) floatptrout, (((size_t) lenvec)*(sizeof(float)))) < (((size_t) lenvec)*(sizeof(float)))) {
	      Py_CLEAR(tmpcopyarray);
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      return 1;
	    }
	    Py_CLEAR(tmpcopyarray);
	  } else {
	    floatptrout = (float *) PyArray_DATA(tmparray);
	    if(write(cparent->sockets[threadindex][1], (void *) dblptrout, (((size_t) lenvec)*(sizeof(float)))) < (((size_t) lenvec)*(sizeof(float)))) {
	      Py_CLEAR(tmpcopyarray);
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      return 1;
	    }
	  }
	  break;
	case VARTOOLS_TYPE_INT:
	  if(descr->type_num != NPY_INT) {
	    tmpcopyarray = PyArray_SimpleNew(1, outdims, NPY_INT);
	    if(PyArray_CopyInto((PyArrayObject *) tmpcopyarray, (PyArrayObject *) tmparray)) {
	      fprintf(stderr,"Error converting vector %s to the expected int datatype in output from python.\n",v->varname);
	      Py_CLEAR(tmpcopyarray);
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      return 1;
	    }
	    intptrout = (int *) PyArray_DATA(tmpcopyarray);
	    if(write(cparent->sockets[threadindex][1], (void *) intptrout, (((size_t) lenvec)*(sizeof(int)))) < (((size_t) lenvec)*(sizeof(int)))) {
	      Py_CLEAR(tmpcopyarray);
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      return 1;
	    }
	    Py_CLEAR(tmpcopyarray);
	  } else {
	    intptrout = (int *) PyArray_DATA(tmparray);
	    if(write(cparent->sockets[threadindex][1], (void *) intptrout, (((size_t) lenvec)*(sizeof(int)))) < (((size_t) lenvec)*(sizeof(int)))) {
	      Py_CLEAR(tmpcopyarray);
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      return 1;
	    }
	  }
	  break;
	case VARTOOLS_TYPE_LONG:
	  if(descr->type_num != NPY_LONG) {
	    tmpcopyarray = PyArray_SimpleNew(1, outdims, NPY_LONG);
	    if(PyArray_CopyInto((PyArrayObject *) tmpcopyarray, (PyArrayObject *) tmparray)) {
	      fprintf(stderr,"Error converting vector %s to the expected long datatype in output from python.\n",v->varname);
	      Py_CLEAR(tmpcopyarray);
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      return 1;
	    }
	    longptrout = (long *) PyArray_DATA(tmpcopyarray);
	    if(write(cparent->sockets[threadindex][1], (void *) longptrout, (((size_t) lenvec)*(sizeof(int)))) < (((size_t) lenvec)*(sizeof(long)))) {
	      Py_CLEAR(tmpcopyarray);
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      return 1;
	    }
	    Py_CLEAR(tmpcopyarray);
	  } else {
	    longptrout = (long *) PyArray_DATA(tmparray);
	    if(write(cparent->sockets[threadindex][1], (void *) longptrout, (((size_t) lenvec)*(sizeof(int)))) < (((size_t) lenvec)*(sizeof(long)))) {
	      Py_CLEAR(tmpcopyarray);
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      return 1;
	    }
	  }
	  break;
	case VARTOOLS_TYPE_SHORT:
	  if(descr->type_num != NPY_SHORT) {
	    tmpcopyarray = PyArray_SimpleNew(1, outdims, NPY_SHORT);
	    if(PyArray_CopyInto((PyArrayObject *) tmpcopyarray, (PyArrayObject *) tmparray)) {
	      fprintf(stderr,"Error converting vector %s to the expected short datatype in output from python.\n",v->varname);
	      Py_CLEAR(tmpcopyarray);
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      return 1;
	    }
	    shortptrout = (short *) PyArray_DATA(tmpcopyarray);
	    if(write(cparent->sockets[threadindex][1], (void *) longptrout, (((size_t) lenvec)*(sizeof(short)))) < (((size_t) lenvec)*(sizeof(short)))) {
	      Py_CLEAR(tmpcopyarray);
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      return 1;
	    }
	    Py_CLEAR(tmpcopyarray);
	  } else {
	    shortptrout = (short *) PyArray_DATA(tmparray);
	    if(write(cparent->sockets[threadindex][1], (void *) shortptrout, (((size_t) lenvec)*(sizeof(short)))) < (((size_t) lenvec)*(sizeof(short)))) {
	      Py_CLEAR(tmpcopyarray);
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      return 1;
	    }
	  }
	  break;
	}
      }
    } else {
      tmplist = PyList_GetItem(py->VariableListOut,ii);
      if(!PyList_Check(tmplist)) {
	fprintf(stderr,"Error: expected a List from python storing the output value(s) for the vartools string or character variable '%s', but received a different type object.\n", v->varname);
	if(tmpcharvec != NULL) free(tmpcharvec);
	return 1;
      }
      lenvec = (int) PyList_Size(tmplist);
      if(v->vectortype != VARTOOLS_VECTORTYPE_LC && lenvec != 1) {
	if(lenvec != 1) {
	  fprintf(stderr,"Error the variable %s has the wrong dimension upon output from python. This variable is expected to be a scalar, but was instead found to have dimension %d\n",v->varname,lenvec);
	  if(tmpcharvec != NULL) free(tmpcharvec);
	  return 1;
	}
      }
      if(write(cparent->sockets[threadindex][1], &lenvec, sizeof(int)) < sizeof(int)) {
	if(tmpcharvec != NULL) free(tmpcharvec);
	return 1;
      }
      for(k=0; k < lenvec; k++) {
#ifdef HAVE_PYTHON_UNICODE
	if(!PyUnicode_Check(PyList_GetItem(tmplist,k))) {
#else
	  if(!PyString_Check(PyList_GetItem(tmplist,k))) {
#endif
	    /* Will need to convert the result to a string */
	  fprintf(stderr,"Error: the variable %s is expected to be a string or ASCII char, but a non-string value was returned from python.\n", v->varname);
	  if(tmpcharvec != NULL) free(tmpcharvec);
	  return 1;
	}
	else {
#ifdef HAVE_PYTHON_UNICODE
	  tmpstr = PyBytes_AS_STRING(PyUnicode_AsEncodedString(PyList_GetItem(tmplist,k), "UTF-8", "strict"));
#else
	  tmpstr = PyString_AsString(PyList_GetItem(tmplist,k));
#endif
	  if(tmpstr == NULL) {
	    fprintf(stderr,"Error: NULL value returned from python for string or ASCII char variable %s\n", v->varname);
	    if(tmpcharvec != NULL) free(tmpcharvec);
	    return 1;
	  }
	  if(v->datatype == VARTOOLS_TYPE_CHAR) {
	    if(lenvec <= 1) {
	      outchar = tmpstr[0];
	      if(write(cparent->sockets[threadindex][1], &outchar, sizeof(char)) < sizeof(char)) {
		if(tmpcharvec != NULL) free(tmpcharvec);
		return 1;
	      }
	    } else {
	      if(!k) {
		if(lenvec > sizetmpcharvec) {
		  if(!sizetmpcharvec) {
		    if((tmpcharvec = (char *) malloc(lenvec * sizeof(char))) == NULL)
		      error(ERR_MEMALLOC);
		  } else {
		    if((tmpcharvec = (char *) realloc(tmpcharvec, lenvec * sizeof(char))) == NULL)
		      error(ERR_MEMALLOC);
		  }
		  sizetmpcharvec = lenvec;
		}
	      }
	      tmpcharvec[k] = tmpstr[0];
	      if(k == lenvec - 1) {
		if(write(cparent->sockets[threadindex][1], tmpcharvec, (((size_t) lenvec)*sizeof(char))) < (((size_t) lenvec)*sizeof(char))) {
		  if(tmpcharvec != NULL) free(tmpcharvec);
		  return 1;
		}
	      }
	    }
	  } else {
	    outstrlen = strlen(tmpstr);
	    if(write(cparent->sockets[threadindex][1], &outstrlen, sizeof(int)) < sizeof(int)) {
	      if(tmpcharvec != NULL) free(tmpcharvec);
	      return 1;
	    }
	    if(outstrlen > 0) {
	      if(write(cparent->sockets[threadindex][1], tmpstr, (((size_t) (outstrlen))*sizeof(char))) < (((size_t) (outstrlen))*sizeof(char))) {
		if(tmpcharvec != NULL) free(tmpcharvec);
		return 1;
	      }
	    }
	  }
	}
      }
    }
    i++; ii++;
  }
  /* The data has all been transferred back to the parent process; we can
     clean-up the Python structures associated with this particular light
     curve */
  CleanUpPythonObjectContainerVariables(py);
  
  if(py->VariableListOut != NULL) {
    if(py->VariableListOut->ob_refcnt > 0)
      Py_CLEAR(py->VariableListOut);
    py->VariableListOut = NULL;
  }
  

  if(tmpcharvec != NULL) free(tmpcharvec);
  return 0;
}

int WriteVariablesFromPythonToSocket_all_lcs(ProgramData *p, 
					     _PythonCommand *cparent, 
					     int threadindex, int cid) {
  int Nlcs, i;
  PyObject* FullList;
  _PythonObjectContainer *py = NULL;

  _PythonCommand *c;

  if(!cid) c = cparent;
  else c = ((_PythonCommand **)cparent->childcommandptrs)[cid-1];
  py = ((_PythonObjectContainer **) cparent->pythonobjects)[threadindex];

  FullList =  py->FullListOut;
  Nlcs = (int) PyList_Size(FullList);

  if(write(cparent->sockets[threadindex][1], &Nlcs, sizeof(int)) < sizeof(int))
    return 1;

  for(i=0; i < Nlcs; i++) {
    py->VariableListOut = PyList_GetItem(FullList,i);
    Py_INCREF(py->VariableListOut);
    if(WriteVariablesFromPythonToSocket(p, cparent, threadindex, cid))
      return 1;
  }
  
  Py_CLEAR(FullList);
  return 0;
}

int RunPythonUserFunctionOnLightCurve(ProgramData *p, _PythonCommand *c, 
				      int threadindex, int cid, int isall) {
  PyObject *ArgTuple;
  int i;
  _PythonObjectContainer *py;

  py = ((_PythonObjectContainer **) c->pythonobjects)[threadindex];

  /* Create the argument tuple which will be passed to the user function */
  if(!isall)
    ArgTuple = PyTuple_Pack(1, py->VariableList);
  else
    ArgTuple = PyTuple_Pack(1, py->FullList);
  //Py_DECREF(py->VariableList);

  /* Set the list as the single item in the argument tuple */
  /*if(PyTuple_SetItem(ArgTuple, 0, py->VariableList)) {
      fprintf(stderr,"Error running python command.\n");
      if(PyErr_Occurred() != NULL) {
	fprintf(stderr,"The error message from python is as follows:\n\n");
	PyErr_Print();
      } else {
	fprintf(stderr,"No error reported by python");
      }
      Py_CLEAR(ArgTuple);
      return 1;
      }*/

  /* Call the User Function */
  if(!isall) {
    py->VariableListOut = PyObject_CallObject(py->UserFunctionToRun[cid], ArgTuple);
    if(py->VariableListOut == NULL) {
      fprintf(stderr,"Error running python command.\n");
      if(PyErr_Occurred() != NULL) {
	fprintf(stderr,"The error message from python is as follows:\n\n");
	PyErr_Print();
      } else {
	fprintf(stderr,"No error reported by python");
      }
      Py_CLEAR(ArgTuple);
      return 1;
    }
    Py_INCREF(py->VariableListOut);
  }
  else {
    py->FullListOut = PyObject_CallObject(py->UserFunctionToRun[cid], ArgTuple);
    if(py->FullListOut == NULL) {
      fprintf(stderr,"Error running python command.\n");
      if(PyErr_Occurred() != NULL) {
	fprintf(stderr,"The error message from python is as follows:\n\n");
	PyErr_Print();
      } else {
	fprintf(stderr,"No error reported by python");
      }
      Py_CLEAR(ArgTuple);
      return 1;
    }
    Py_INCREF(py->FullListOut);
  }

  /* Clean up the temporary variables that were allocated */
  Py_CLEAR(ArgTuple);
  return 0;
}

void TerminatePython(ProgramData *p, _PythonCommand *c, int threadindex)
{
  c->IsPythonRunning[threadindex] = 0;
  if(!c->iscontinueprocess) {
   /* PyThreadState_Swap((((_PythonObjectContainer **) c->pythonobjects)[threadindex])->mainstate);*/
    CleanUpPythonObjectContainer(&(((_PythonObjectContainer **) c->pythonobjects)[threadindex]));
    Py_Finalize();
  } else {
    TerminatePython(p, (_PythonCommand *) c->continueprocesscommandptr, threadindex);
  }
}

void RunPythonProcessingLoop(ProgramData *p, _PythonCommand *c, int threadindex)
{
  int msg, retval, newcnum, cid;
  _PythonCommand *ccheck;

  while(1) {
    /* Read a message from the parent indicating what we will be doing */
    /* If no message is received, terminate the loop */
    if(read(c->sockets[threadindex][1], &msg, sizeof(int)) < sizeof(int))
      return;
  
    switch(msg) {
    case VARTOOLS_PYTHON_MESSAGE_ENDPROCESS:
      return;
    case VARTOOLS_PYTHON_MESSAGE_READDATA:

      if(read(c->sockets[threadindex][1], &cid, sizeof(int)) < sizeof(int))
	return;

      if(cid > 0 && (cid - 1 > c->Nchildren))
	return;

      if(!cid) ccheck = c;
      else ccheck = ((_PythonCommand **)c->childcommandptrs)[cid-1];

      if(ccheck->RequireReadAll) {
	if((retval = ReadVariablesFromSocketIntoPython_all_lcs(p, c, threadindex))) {
	  if(write(c->sockets[threadindex][1], &retval, sizeof(int)) < sizeof(int))
	    return;
	  break;
	}
	if((retval = RunPythonUserFunctionOnLightCurve(p, c, threadindex, cid, 1))) {
	  if(write(c->sockets[threadindex][1], &retval, sizeof(int)) < sizeof(int))
	    return;
	  break;
	}
      } else {
	if((retval = ReadVariablesFromSocketIntoPython(p, c, threadindex, 0, 0))) {
	  if(write(c->sockets[threadindex][1], &retval, sizeof(int)) < sizeof(int))
	    return;
	  break;
	}
	if((retval = RunPythonUserFunctionOnLightCurve(p, c, threadindex, cid, 0))) {
	  if(write(c->sockets[threadindex][1], &retval, sizeof(int)) < sizeof(int))
	    return;
	  break;
	}
      }
      if(write(c->sockets[threadindex][1], &retval, sizeof(int)) < sizeof(int))
	return;
      if(ccheck->RequireReadAll) {
	if(WriteVariablesFromPythonToSocket_all_lcs(p, c, threadindex, cid))
	  return;
      } else {
	if(WriteVariablesFromPythonToSocket(p, c, threadindex, cid))
	  return;
      }
      break;
    default:
      fprintf(stderr,"Error: Invalid message received by python sub-process from main vartools process. Shutting down the python sub-process now.\n");
      return;
    }
  }
}


void StartPythonProcess(ProgramData *p, _PythonCommand *c, int threadindex)
{
  pid_t pid;

  /* If Python is already running on this thread, just return */
  if(c->IsPythonRunning[threadindex]) return;
  
  c->IsPythonRunning[threadindex] = 1;

  /* Create a socket pair */
  if(socketpair(AF_UNIX, SOCK_STREAM, 0, c->sockets[threadindex]) < 0) {
    perror("Error Opening socket pair in starting python process\n");
    exit(1);
  }

  /* Create the sub-process */
  pid = fork();
  if(pid == -1) {
    perror("Error creating a sub-process initiating a call to the vartools python command.\n");
    exit(1);
  }
  if(!pid) {
    /* This is the child sub-process which will be used for running python */

    /* Close the parent end of the socket */
    close(c->sockets[threadindex][0]);    

    /* Redirect stdout to stderr */
    dup2(2, 1);

    /* Start python, load all of the
       initialization codes and create the function that will be used for
       processing individual light curves */
    InitializePython(p, c, threadindex);

    /* The loop that will receive light curve data from the parent process,
       run the python routine, and write the data back to the parent */

    RunPythonProcessingLoop(p, c, threadindex);

    /* If the loop is done running we can terminate this process */
    TerminatePython(p, c, threadindex);

    /* Close the socket, then terminate the process */
    close(c->sockets[threadindex][1]);
    exit(0);
  } else {

    /* This is the main parent loop. In this case we close the child end of 
       the socket, and then continue with the processing */
    close(c->sockets[threadindex][1]);
  }
}


int SendVariablesToChildPythonProcess(ProgramData *p, int lcindex,
				      int threadindex, int pythreadindex, _PythonCommand *c)
{
  
  int i, k;
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

  if(write(c->sockets[pythreadindex][0], &c->Nvars, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Error sending the number of variables to read to the child python process.\n");
    if(tmpstringval != NULL) free(tmpstringval);
    return 1;
  }

  for(i=0; i < c->Nvars; i++) {
    if(write(c->sockets[pythreadindex][0], &(c->vars[i]->datatype), sizeof(char)) < sizeof(char)) {
      fprintf(stderr,"Error sending the datatype for variable %s to the child python process\n", c->vars[i]->varname);
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
	fprintf(stderr,"Error: invalid datatype received for variable %s in the SendVariablesToChildPythonProcess function\n", c->vars[i]->varname);
	if(tmpstringval != NULL) free(tmpstringval);
	return 1;
      }
      if(c->vars[i]->vectortype != VARTOOLS_VECTORTYPE_LC)
	lenvec = 1;
      else
	lenvec = p->NJD[threadindex];
      
      if(write(c->sockets[pythreadindex][0], &lenvec, sizeof(int)) < sizeof(int)) {
	fprintf(stderr,"Error sending the length of vector for variable %s in the SendVariablesToChildPythonProcess function.\n", c->vars[i]->varname);
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

      if(write(c->sockets[pythreadindex][0], ptrtosend, (((size_t) lenvec)*databytesize)) < (((size_t) lenvec)*databytesize)) {
	fprintf(stderr,"Error sending data for variable %s in the SendVariablesToChildPythonProcess function.\n", c->vars[i]->varname);
	if(tmpstringval != NULL) free(tmpstringval);
	return 1;
      }
    } else {
      databytesize = sizeof(char);
      if(c->vars[i]->vectortype != VARTOOLS_VECTORTYPE_LC)
	lenvec = 1;
      else
	lenvec = p->NJD[threadindex];
      
      if(write(c->sockets[pythreadindex][0], &lenvec, sizeof(int)) < sizeof(int)) {
	fprintf(stderr,"Error sending the length of vector for variable %s in the SendVariablesToChildPythonProcess function.\n", c->vars[i]->varname);
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
		error(ERR_MEMALLOC);
	    } else {
	      if((tmpstringval = (char *) realloc(tmpstringval, c->vars[i]->outc->stringsize)) == NULL)
		error(ERR_MEMALLOC);
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
	if(write(c->sockets[pythreadindex][0], &strlentosend, sizeof(int)) < sizeof(int)) {
	  fprintf(stderr,"Error sending the string length for variable %s in the SendVariablesToChildPythonProcess function.\n", c->vars[i]->varname);
	  if(tmpstringval != NULL) free(tmpstringval);
	  return 1;
	}
	if(strlentosend <= 0) continue;
	if(write(c->sockets[pythreadindex][0], ptrtosend, (((size_t) strlentosend)*sizeof(char))) < (((size_t) strlentosend)*sizeof(char))) {
	  fprintf(stderr,"Error sending the data for variable %s in the SendVariablesToChildPythonProcess function.\n", c->vars[i]->varname);
	  if(tmpstringval != NULL) free(tmpstringval);
	  return 1;
	}
      } else {
	for(k=0; k < lenvec; k++) {
	  ptrtosend = (void *) &((*((char ****) c->vars[i]->dataptr))[threadindex][k][0]);
	  strlentosend = strlen((char *) ptrtosend);
	  if(write(c->sockets[pythreadindex][0], &strlentosend, sizeof(int)) < sizeof(int)) {
	    fprintf(stderr,"Error sending the string length for variable %s in the SendVariablesToChildPythonProcess function.\n", c->vars[i]->varname);
	    if(tmpstringval != NULL) free(tmpstringval);
	    return 1;
	  }
	  if(strlentosend <= 0) continue;
	  if(write(c->sockets[pythreadindex][0], ptrtosend, (((size_t) strlentosend)*sizeof(char))) < (((size_t) strlentosend)*sizeof(char))) {
	    fprintf(stderr,"Error sending the data for variable %s in the SendVariablesToChildPythonProcess function.\n", c->vars[i]->varname);
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

int SendVariablesToChildPythonProcess_all_lcs(ProgramData *p, _PythonCommand *c)
{
  int j, Nlcs;

  Nlcs = p->Nlcs;
  if(write(c->sockets[0][0], &Nlcs, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Error sending the number of light curves to read to the child python process.\n");
    return 1;
  }
  for(j=0; j < Nlcs; j++) {
    if(SendVariablesToChildPythonProcess(p, j, j, 0, c))
      return 1;
  }
  return 0;
}

int ReadVariablesFromChildPythonProcess(ProgramData *p, int lcindex, 
					int threadindex, int pythreadindex, 
					_PythonCommand *c)
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
      fprintf(stderr,"Error: somehow a constant or outcolumn variable is set to receive data from a python command. This should not happen and indicates a bug in vartools.\n");
      exit(1);
    }
    
    if(read(c->sockets[pythreadindex][0], &lenvec, sizeof(int)) < sizeof(int)) {
      fprintf(stderr,"Error receiving the vector length for variable %s in the ReadVariablesFromChildPythonProcess function.\n", v->varname);
      if(tmpinstr != NULL) free(tmpinstr);
      return 1;
    }
    
    if(v->datatype != VARTOOLS_TYPE_STRING) {
      
      if(v->vectortype != VARTOOLS_VECTORTYPE_LC) {
	if(lenvec != 1) {
	  fprintf(stderr,"Error the vector length received for variable %s has the wrong dimension in the ReadVariablesFromChildPythonProcess function. Expected only 1 value, received a length of %d\n",v->varname,lenvec);
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
	if(read(c->sockets[pythreadindex][0], ptrtoget, databytesize) < databytesize) {
	  fprintf(stderr,"Error receiving the data for variable %s in the ReadVariablesFromChildPythonProcess function.\n",v->varname);
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
	  sizedataread = read(c->sockets[pythreadindex][0], (void *) (((char *) ptrtoget) + sizedatareadtot), (((size_t)lenvec)*databytesize)-sizedatareadtot);
	  sizedatareadtot += sizedataread;
	}
	if(sizedatareadtot < (((size_t)lenvec)*databytesize)) {
	  fprintf(stderr,"Error receiving the data for variable %s in the ReadVariablesFromChildPythonProcess function.\n",v->varname);
	  if(tmpinstr != NULL) free(tmpinstr);
	  return 1;
	}
      }
    } else {
      if(v->vectortype != VARTOOLS_VECTORTYPE_LC) {
	if(lenvec != 1) {
	  fprintf(stderr,"Error the vector length received for variable %s has the wrong dimension in the ReadVariablesFromChildPythonProcess function. Expected only 1 value, received a length of %d\n",v->varname,lenvec);
	  if(tmpinstr != NULL) free(tmpinstr);
	  return 1;
	}
	if(read(c->sockets[pythreadindex][0], &instrlen, sizeof(int)) < sizeof(int)) {
	  fprintf(stderr,"Error receiving the string length for variable %s in the ReadVariablesFromChildPythonProcess function.\n",v->varname);
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
	      sizedataread = read(c->sockets[pythreadindex][0], (void *) (((char *) tmpinstr) + sizedatareadtot), (((size_t) instrlen)*sizeof(char))-sizedatareadtot);
	      sizedatareadtot += sizedataread;
	    }
	    if(sizedatareadtot < (((size_t) instrlen)*sizeof(char))) {
	      fprintf(stderr,"Error receiving data for variable %s in the ReadVariablesFromChildPythonProcess function.\n",v->varname);
	      if(tmpinstr != NULL) free(tmpinstr);
	      return 1;
	    }
	    for(ii=0; ii < MAXLEN-1; ii++)
	      ((char *)ptrtoget)[ii] = tmpinstr[ii];
	    ((char *)ptrtoget)[ii] = '\0';
	  } else {
	    sizedatareadtot = 0;
	    while(sizedatareadtot < (((size_t) instrlen)*sizeof(char))) {
	      sizedataread = read(c->sockets[pythreadindex][0], (void *) (((char *) ptrtoget) + sizedatareadtot), (((size_t)instrlen)*sizeof(char))-sizedatareadtot);
	      sizedatareadtot += sizedataread;
	    }
	    if(sizedatareadtot < (((size_t) instrlen)*sizeof(char))) {
	      fprintf(stderr,"Error receiving data for variable %s in the ReadVariablesFromChildPythonProcess function.\n",v->varname);
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
	  if(read(c->sockets[pythreadindex][0], &instrlen, sizeof(int)) < sizeof(int)) {
	    fprintf(stderr,"Error receiving the string length for variable %s in the ReadVariablesFromChildPythonProcess function.\n",v->varname);
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
		sizedataread = read(c->sockets[pythreadindex][0], (void *) (((char *) tmpinstr) + sizedatareadtot), (((size_t)instrlen)*sizeof(char))-sizedatareadtot);
		sizedatareadtot += sizedataread;
	      }
	      if(sizedatareadtot < (((size_t) instrlen)*sizeof(char))) {
		fprintf(stderr,"Error receiving data for variable %s in the ReadVariablesFromChildPythonProcess function.\n",v->varname);
		if(tmpinstr != NULL) free(tmpinstr);
		return 1;
	      }
	      for(ii=0; ii < MAXLEN-1; ii++)
		((char *)ptrtoget)[ii] = tmpinstr[ii];
	      ((char *)ptrtoget)[ii] = '\0';
	    } else {
	      sizedatareadtot = 0;
	      while(sizedatareadtot < (((size_t) instrlen)*sizeof(char))) {
		sizedataread = read(c->sockets[pythreadindex][0], (void *) (((char *) ptrtoget) + sizedatareadtot), (((size_t)instrlen)*sizeof(char))-sizedatareadtot);
		sizedatareadtot += sizedataread;
	      }
	      if(sizedatareadtot < (((size_t) instrlen)*sizeof(char))) {
		fprintf(stderr,"Error receiving data for variable %s in the ReadVariablesFromChildPythonProcess function.\n",v->varname);
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


int ReadVariablesFromChildPythonProcess_all_lcs(ProgramData *p, 
						_PythonCommand *c)
{
  int i, j, k, ii, tmpindex;
  int Nlcsout;
  int Nlcs = p->Nlcs;

  if(read(c->sockets[0][0], &Nlcsout, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Error receiving the number of light curves to read to the child python process.\n");
    return 1;
  }
  if(Nlcsout != p->Nlcs) {
    fprintf(stderr,"Error: the number of light curves implied by the length of the variables returned by python does not match the number of light curves sent to python.\n");
    return 1;
  }
  for(j=0; j < Nlcs; j++) {
    if(ReadVariablesFromChildPythonProcess(p, j, j, 0, c))
      return 1;
  }
  return 0;
}


void StopRunningPythonCommand(ProgramData *p, int threadindex, _PythonCommand *c)
{
  int msg = VARTOOLS_PYTHON_MESSAGE_ENDPROCESS;
  if(!c->iscontinueprocess) {
    if(c->IsPythonRunning[threadindex]) {
      write(c->sockets[threadindex][0], &msg, sizeof(int));
      close(c->sockets[threadindex][0]);
    }
  } else {
    if(!((_PythonCommand *)c->continueprocesscommandptr)->IsPythonRunning[threadindex]) {
      write(c->sockets[threadindex][0], &msg, sizeof(int));
      close(c->sockets[threadindex][0]);
    }
  }
  return;
}

void GetPythonCommandOutputColumnValues(ProgramData *p, int lcindex, int threadindex, _PythonCommand *c)
{
  int i;
  for(i=0; i < c->Noutcolumnvars; i++) {
    c->outcolumndata[threadindex][i] = EvaluateVariable_Double(lcindex, threadindex, 0, c->outcolumnvars[i]);
  }
}

void GetPythonCommandOutputColumnValues_all_lcs(ProgramData *p, _PythonCommand *c)
{
  int i, j;
  for(j=0; j < p->Nlcs; j++) {
    for(i=0; i < c->Noutcolumnvars; i++) {
      c->outcolumndata[j][i] = EvaluateVariable_Double(j, j, 0, c->outcolumnvars[i]);
    }
  }
}


void RunPythonCommand(ProgramData *p, int lcindex, int threadindex, int pythreadindex, _PythonCommand *c)
{
  /* Check if the python process for this command is running for this
     thread; if not start it */
  int msg, retval;

  if(!c->iscontinueprocess) {

    if(!c->IsPythonRunning[pythreadindex]) 
      StartPythonProcess(p, c, pythreadindex);
    
  } else {
    if(!((_PythonCommand *)c->continueprocesscommandptr)->IsPythonRunning[pythreadindex]) {
      StartPythonProcess(p, ((_PythonCommand *)c->continueprocesscommandptr), pythreadindex);
    }
    
    if(!c->IsPythonRunning[pythreadindex]) {
      c->IsPythonRunning[pythreadindex] = 1;
      c->sockets[pythreadindex][0] = ((_PythonCommand *)c->continueprocesscommandptr)->sockets[pythreadindex][0];
    }
    
  }

  msg = VARTOOLS_PYTHON_MESSAGE_READDATA;
  if(write(c->sockets[pythreadindex][0], &msg, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Lost communication with python subprocess\n");
    close(c->sockets[pythreadindex][0]);
    exit(1);
  }

  if(write(c->sockets[pythreadindex][0], &c->cid, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Lost communication with python subprocess\n");
    close(c->sockets[pythreadindex][0]);
    exit(1);
  }
  
  if(SendVariablesToChildPythonProcess(p, lcindex, threadindex, pythreadindex, c)) {
    /* The process failed; terminate with an error */
    fprintf(stderr,"Failed to send variables to python subprocess\n");
    close(c->sockets[pythreadindex][0]);
    exit(1);
  }

  if(read(c->sockets[pythreadindex][0], &retval, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Lost communication with python subprocess\n");
    close(c->sockets[pythreadindex][0]);
    exit(1);
  }
  
  if(retval) {
    fprintf(stderr,"Python function %d returned an error for light curve number %d\n", c->cnum, lcindex);
    exit(1);
  }
  
  if(ReadVariablesFromChildPythonProcess(p, lcindex, threadindex, pythreadindex, c)) {
    /* The process failed; terminate with an error */
    fprintf(stderr,"Failed to receive output variables from the python subprocess\n");
    close(c->sockets[pythreadindex][0]);
    exit(1);
  }
  GetPythonCommandOutputColumnValues(p, lcindex, threadindex, c);
  
}


void RunPythonCommand_all_lcs(ProgramData *p, _PythonCommand *c)
{
  /* Check if the python process for this command is running for this
     thread; if not start it */
  int msg, retval;

  if(!c->iscontinueprocess) {

    if(!c->IsPythonRunning[0]) 
      StartPythonProcess(p, c, 0);
    
  } else {
    if(!((_PythonCommand *)c->continueprocesscommandptr)->IsPythonRunning[0]) {
      StartPythonProcess(p, ((_PythonCommand *)c->continueprocesscommandptr), 0);
    }
    
    if(!c->IsPythonRunning[0]) {
      c->IsPythonRunning[0] = 1;
      c->sockets[0][0] = ((_PythonCommand *)c->continueprocesscommandptr)->sockets[0][0];
    }
    
  }

  msg = VARTOOLS_PYTHON_MESSAGE_READDATA;
  if(write(c->sockets[0][0], &msg, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Lost communication with python subprocess\n");
    close(c->sockets[0][0]);
    exit(1);
  }

  if(write(c->sockets[0][0], &c->cid, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Lost communication with python subprocess\n");
    close(c->sockets[0][0]);
    exit(1);
  }
  
  if(SendVariablesToChildPythonProcess_all_lcs(p, c)) {
    /* The process failed; terminate with an error */
    fprintf(stderr,"Failed to send variables to python subprocess\n");
    close(c->sockets[0][0]);
    exit(1);
  }

  if(read(c->sockets[0][0], &retval, sizeof(int)) < sizeof(int)) {
    fprintf(stderr,"Lost communication with python subprocess\n");
    close(c->sockets[0][0]);
    exit(1);
  }
  
  if(retval) {
    fprintf(stderr,"Python function %d returned an error\n", c->cnum);
    exit(1);
  }
  
  if(ReadVariablesFromChildPythonProcess_all_lcs(p, c)) {
    /* The process failed; terminate with an error */
    fprintf(stderr,"Failed to receive output variables from the python subprocess\n");
    close(c->sockets[0][0]);
    exit(1);
  }
  GetPythonCommandOutputColumnValues_all_lcs(p, c);
  
}


void InitPythonCommand(ProgramData *p, _PythonCommand *c, int Nlcs)
{
  int j;
  if(!p->readallflag) {
    if((c->pythonobjects = (void *)(((_PythonObjectContainer **) malloc(Nlcs * sizeof(_PythonObjectContainer *))))) == NULL ||
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
	c->IsPythonRunning[j] = 0;
    } else {
#endif
      c->IsPythonRunning[0] = 0;
#ifdef PARALLEL
    }
#endif
  } else {
    if((c->pythonobjects = (void *)(((_PythonObjectContainer **) malloc(sizeof(_PythonObjectContainer *))))) == NULL ||
       (c->outcolumndata = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
      DO_ERROR_MEMALLOC;
    if(c->Noutcolumnvars > 0) {
      for(j=0; j < Nlcs; j++) {
	if((c->outcolumndata[j] = (double *) malloc(c->Noutcolumnvars * sizeof(double))) == NULL)
	  DO_ERROR_MEMALLOC;
      }
    }
    c->IsPythonRunning[0] = 0;
  }

}


int ParsePythonCommand(int *iret, int argc, char **argv, ProgramData *p, 
			_PythonCommand *c, Command *allcommands, int cnum)
/* Parse the command line for the "-python" command; the expected syntax is:
   -python < "fromfile" commandfile | commandstring >
           [ "init" < "file" initializationfile | initializationstring > |
             "continueprocess" prior_python_command_number ]
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
  _PythonCommand *cparent;

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
      c->len_pythoncommandstring = ftell(infile) + 1;
      if(c->len_pythoncommandstring <= 0) {
	fprintf(stderr,"Cannot Read %s\n", argv[i]);
	fclose(infile);
	exit(1);
	*iret = i; return 1;
      }
      if((c->pythoncommandstring = (char *) malloc((c->len_pythoncommandstring)*sizeof(char))) == NULL) {
	DO_ERROR_MEMALLOC;
	*iret = i; return 1;
      }
      if(fseek(infile, 0L, SEEK_SET)) {
	fprintf(stderr,"Error reading %s\n", argv[i]);
	fclose(infile);
	exit(1);
	*iret = i; return 1;
      }
      if(fread(c->pythoncommandstring, sizeof(char),
	       (c->len_pythoncommandstring-1),
	       infile) < (c->len_pythoncommandstring-1)) {
	fprintf(stderr,"Error reading %s\n", argv[i]);
	fclose(infile);
	exit(1);
	*iret = i; return 1;
      }
      c->pythoncommandstring[c->len_pythoncommandstring-1] = '\0';
    } else {
      fprintf(stderr,"Error reading %s\n", argv[i]);
      fclose(infile);
      exit(1);
      *iret = i; return 1;
    }
  } else {
    c->len_pythoncommandstring = strlen(argv[i])+1;
    if((c->pythoncommandstring = (char *) malloc((c->len_pythoncommandstring)*sizeof(char))) == NULL) {
      DO_ERROR_MEMALLOC;
      *iret = i; return 1;
    }
    sprintf(c->pythoncommandstring,"%s",argv[i]);
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
	  c->len_pythoninitializationtextstring = ftell(infile) + 1;
	  if(c->len_pythoninitializationtextstring <= 0) {
	    fprintf(stderr,"Error reading %s\n", argv[i]);
	    fclose(infile);
	    exit(1);
	    *iret = i; return 1;
	  }
	  if((c->pythoninitializationtext = (char *) malloc((c->len_pythoninitializationtextstring)*sizeof(char))) == NULL) {
	    DO_ERROR_MEMALLOC;
	    *iret = i; return 1;
	  }
	  if(fseek(infile, 0L, SEEK_SET)) {
	    fprintf(stderr,"Error reading %s\n", argv[i]);
	    fclose(infile);
	    exit(1);
	    *iret = i; return 1;
	  }
	  if(fread(c->pythoninitializationtext, sizeof(char),
		   (c->len_pythoninitializationtextstring-1),
		   infile) < (c->len_pythoninitializationtextstring-1)) {
	    fprintf(stderr,"Error reading %s\n", argv[i]);
	    fclose(infile);
	    exit(1);
	    *iret = i; return 1;
	  }
	  c->pythoninitializationtext[c->len_pythoninitializationtextstring-1] = '\0';
	} else {
	  fprintf(stderr,"Error reading %s\n", argv[i]);
	  fclose(infile);
	  exit(1);
	  *iret = i; return 1;
	}
      } else {
	c->len_pythoninitializationtextstring = strlen(argv[i])+1;
	if((c->pythoninitializationtext = (char *) malloc((c->len_pythoninitializationtextstring)*sizeof(char))) == NULL) {
	  DO_ERROR_MEMALLOC;
	  *iret = i; return 1;
	}
	sprintf(c->pythoninitializationtext,"%s",argv[i]);
      }
    }
    else if(!strcmp(argv[i],"continueprocess")) {
      c->iscontinueprocess = 1;
      i++;
      if(i >= argc) {*iret = i; return 1;}
      j = atoi(argv[i]);
      if(j <= 0) {
	fprintf(stderr,"Error parsing -python command. The value of prior_python_command_number must be >= 1.\n\n");
	*iret = i; return 1;
      }
      ii = 0;
      /* Find the previous process */
      for(k = 0; k < cnum; k++) {
	if(allcommands[k].cnum == CNUM_PYTHON) {
	  ii++;
	  if(ii == j) {
	    if(allcommands[k].PythonCommand->iscontinueprocess) {
	      c->continueprocesscommandptr = (void *) allcommands[k].PythonCommand->continueprocesscommandptr;
	    } else {
	      c->continueprocesscommandptr = (void *) allcommands[k].PythonCommand;
	    }
	    cparent = (_PythonCommand *) c->continueprocesscommandptr;
	    if(!cparent->Nchildren) {
	      if((cparent->childcommandptrs = (void *) malloc(sizeof(_PythonCommand *))) == NULL ||
		 (cparent->childcnumvals = (int *) malloc(sizeof(int))) == NULL) {
		DO_ERROR_MEMALLOC; *iret = i; return 1;
	      }
	    } else {
	      if((cparent->childcommandptrs = (void *) realloc((_PythonCommand **)cparent->childcommandptrs, (cparent->Nchildren + 1)*sizeof(_PythonCommand *))) == NULL ||
		 (cparent->childcnumvals = (int *) realloc(cparent->childcnumvals, (cparent->Nchildren + 1)*sizeof(int))) == NULL) {
		DO_ERROR_MEMALLOC; *iret = i; return 1;
	      }
	    }
	    ((_PythonCommand **) cparent->childcommandptrs)[cparent->Nchildren] = c;
	    cparent->childcnumvals[cparent->Nchildren] = cnum;
	    cparent->Nchildren += 1;
	    c->cid = cparent->Nchildren;
	    break;
	  }
	}
      }
      if(ii < j) {
	fprintf(stderr,"Error parsing -python command. Attempting to continue processing from python command number %d from command number %d, but only %d previous instances of -python have been given on the command line.\n", j, cnum, ii);
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
	  fprintf(stderr,"Bad variable name \"\" given to the python analytic expression evaluator.\n");
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

void KillAllPythonProcessesOneThread(ProgramData *p, Command *allcommands, int threadindex) {
  int i;
  for(i=0; i < p->Ncommands; i++) {
    if(allcommands[i].cnum == CNUM_PYTHON) {
      StopRunningPythonCommand(p, threadindex, allcommands[i].PythonCommand);
    }
  }
}

void KillAllPythonProcesses(ProgramData *p, Command *allcommands) {
  int i;
#ifdef PARALLEL
  if(p->Nproc_allow > 1) {
    for(i = 0; i < p->Nproc_allow; i++)
      KillAllPythonProcessesOneThread(p, allcommands, i);
  } else {
#endif
    KillAllPythonProcessesOneThread(p, allcommands, 0);
#ifdef PARALLEL
  }
#endif
}
