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

_Variable *SetupScalarVariable(ProgramData *p, char *varname, int datatype)
{
  _Variable *variable;
  variable = CreateVariable(p, varname, (char) datatype, VARTOOLS_VECTORTYPE_SCALAR, NULL);
  RegisterScalarData(p, variable->dataptr, datatype, 0);
  return variable;
}

void RegisterScalarData(ProgramData *p, void *dataptr, int datatype, int Ncolumns)
{
  if(!p->NScalarData) {
    p->ScalarData = (_ScalarData *) malloc(sizeof(_ScalarData));
  } else {
    p->ScalarData = (_ScalarData *) realloc((p->ScalarData),(p->NScalarData+1)*sizeof(_ScalarData));
  }
  p->ScalarData[p->NScalarData].dataptr = dataptr;
  p->ScalarData[p->NScalarData].datatype = datatype;
  p->ScalarData[p->NScalarData].Ncolumns = Ncolumns;
  p->NScalarData += 1;
}

void MemAllocScalarData(ProgramData *p, int Nthreads)
{
  int i, j, k;
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
  _ScalarData *s;
  int Nc;

  for(i=0; i < p->NScalarData; i++) {
    s = &(p->ScalarData[i]);
    Nc = s->Ncolumns;
    if(Nc <= 0) {
      switch(s->datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dblptr = (double **) s->dataptr;
	if(((*dblptr) = (double *) malloc(Nthreads * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_STRING:
	stringptr = (char ***) s->dataptr;
	if(((*stringptr) = (char **) malloc(Nthreads * sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nthreads; j++) {
	  if((((*stringptr)[j]) = (char *) malloc(MAXLEN)) == NULL)
	    error(ERR_MEMALLOC);
	}
	break;
      case VARTOOLS_TYPE_INT:
	intptr = (int **) s->dataptr;
	if(((*intptr) = (int *) malloc(Nthreads * sizeof(int))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_SHORT:
	shortptr = (short **) s->dataptr;
	if(((*shortptr) = (short *) malloc(Nthreads * sizeof(short))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_FLOAT:
	floatptr = (float **) s->dataptr;
	if(((*floatptr) = (float *) malloc(Nthreads * sizeof(float))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_LONG:
	longptr = (long **) s->dataptr;
	if(((*longptr) = (long *) malloc(Nthreads * sizeof(long))) == NULL)
	  error(ERR_MEMALLOC);
	break;
      case VARTOOLS_TYPE_CHAR:
	charptr = (char **) s->dataptr;
	if(((*charptr) = (char *) malloc(Nthreads)) == NULL)
	  error(ERR_MEMALLOC);
	break;
      default:
	error(ERR_BADTYPE);
      }
    } else {
      switch(s->datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dbl2ptr = (double ***) s->dataptr;
	if(((*dbl2ptr) = (double **) malloc(Nthreads * sizeof(double *))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j<Nthreads; j++) {
	  if((((*dbl2ptr)[j]) = (double *) malloc(Nc * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
	break;
      case VARTOOLS_TYPE_STRING:
	string2ptr = (char ****) s->dataptr;
	if(((*string2ptr) = (char ***) malloc(Nthreads * sizeof(char **))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nthreads; j++) {
	  if((((*string2ptr)[j]) = (char **) malloc(Nc * sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(k=0; k < Nc; k++) {
	    if((((*string2ptr)[j][k]) = (char *) malloc(MAXLEN)) == NULL)
	      error(ERR_MEMALLOC);
	  }
	}
	break;
      case VARTOOLS_TYPE_INT:
	int2ptr = (int ***) s->dataptr;
	if(((*int2ptr) = (int **) malloc(Nthreads * sizeof(int *))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nthreads; j++) {
	  if((((*int2ptr)[j]) = (int *) malloc(Nc * sizeof(int))) == NULL)
	    error(ERR_MEMALLOC);
	}
	break;
      case VARTOOLS_TYPE_SHORT:
	short2ptr = (short ***) s->dataptr;
	if(((*short2ptr) = (short **) malloc(Nthreads * sizeof(short *))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nthreads; j++) {
	  if((((*short2ptr)[j]) = (short *) malloc(Nc * sizeof(short))) == NULL)
	    error(ERR_MEMALLOC);
	}
	break;
      case VARTOOLS_TYPE_FLOAT:
	float2ptr = (float ***) s->dataptr;
	if(((*float2ptr) = (float **) malloc(Nthreads * sizeof(float *))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nthreads; j++) {
	  if((((*float2ptr)[j]) = (float *) malloc(Nc * sizeof(float))) == NULL)
	    error(ERR_MEMALLOC);
	}
	break;
      case VARTOOLS_TYPE_LONG:
	long2ptr = (long ***) s->dataptr;
	if(((*long2ptr) = (long **) malloc(Nthreads * sizeof(long *))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nthreads; j++) {
	  if((((*long2ptr)[j]) = (long *) malloc(Nc * sizeof(long))) == NULL)
	    error(ERR_MEMALLOC);
	}
	break;
      case VARTOOLS_TYPE_CHAR:
	char2ptr = (char ***) s->dataptr;
	if(((*char2ptr) = (char **) malloc(Nthreads * sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
	for(j=0; j < Nthreads; j++) {
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

void InitCommands(ProgramData *p, Command *c)
{

  int i, j, k, Nlcs, Ncommands;
  FILE *inlist;
  char **inputlistlines, *teststring;
  size_t *size_inputlistline;

  int Nkillharmterms;
  int sizeperptrs;
  int *Npers;
  double ***perptrs;
  int sizeinputlistvec = 1000;

  Ncommands = p->Ncommands;

  /* Read in the light curve list.  allocate memory for all the various arrays within the commands. The amount of memory to allocate depends on whether or not we're reading in all the light curves. */

  if(p->listflag && !p->headeronly && !p->inputlistformat &&
     !p->showinputlcformat)
    {
      if(!p->readfromstdinflag)
	{
	  if((inlist = fopen(p->lclistname,"r")) == NULL)
	    error2(ERR_FILENOTFOUND,p->lclistname);
	}
      else
	inlist = stdin;

      /* Read the file first into a buffer */
      if((inputlistlines = (char **) malloc(sizeinputlistvec * sizeof(char *))) == NULL ||
	 (size_inputlistline = (size_t *) malloc(sizeinputlistvec * sizeof(size_t))) == NULL)
	error(ERR_MEMALLOC);
      for(j=0;j<sizeinputlistvec;j++) {
	size_inputlistline[j] = MAXLEN;
	if((inputlistlines[j] = (char *) malloc(MAXLEN * sizeof(char))) == NULL)
	  error(ERR_MEMALLOC);
      }

      p->Nlcs = 0;
      while(gnu_getline(&(inputlistlines[p->Nlcs]), &(size_inputlistline[p->Nlcs]), inlist) >= 0)
	{
	  if(inputlistlines[p->Nlcs][0] != '#')
	    p->Nlcs++;
	  /* Increase the size of the buffer if we need more space for the next line */
	  if(p->Nlcs >= sizeinputlistvec)
	    {
	      sizeinputlistvec += 1000;
	      if((inputlistlines = (char **) realloc(inputlistlines, sizeinputlistvec * sizeof(char *))) == NULL ||
		 (size_inputlistline = (size_t *) realloc(size_inputlistline, sizeinputlistvec * sizeof(size_t))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=sizeinputlistvec - 1000; j < sizeinputlistvec; j++)
		{
		  size_inputlistline[j] = MAXLEN;
		  if((inputlistlines[j] = (char *) malloc(MAXLEN * sizeof(char))) == NULL)
		    error(ERR_MEMALLOC);
		}
	    }
	}

      Nlcs = p->Nlcs;
    }
  else {
    p->Nlcs = 1;
    Nlcs = 1;
  }

  if(p->Ncopycommands > 0) {
    p->Nlcs = p->Nlcs * p->Ncopiestotal;
    Nlcs = p->Nlcs;
    if(p->fileflag) {
      if((p->lcnames = (char **) realloc(p->lcnames,Nlcs * sizeof(char *))) == NULL ||
	 (p->NJD = (int *) realloc(p->NJD,Nlcs * sizeof(int))) == NULL)
	error(ERR_MEMALLOC);
      for(j=1; j < p->Nlcs; j++) {
	if((p->lcnames[j] = (char *) malloc(MAXLEN * sizeof(char))) == NULL)
	  error(ERR_MEMALLOC);
      }
    }
  }

  /* Determine the input lc columns */
  DetermineColumns(p, c);

  /* Prepare the output columns, and link them to command inputs if needed */
  p->Ncolumns = 0;
  CreateOutputColumns(p, c, p->Ncommands);
  linkcolumns(p);


 /* Compile any analytic expressions which will need to be evaluated */
  CompileAllExpressions(p, c);

  /* Allocate memory for the data to be read in from the input list */
  MemAllocDataFromInputList(p, p->Nlcs);

    /* Parse the input list */
  if(p->listflag && !p->headeronly && !p->inputlistformat &&
     !p->showinputlcformat)
    {
      ParseInputList(p, inputlistlines, Nlcs);

      if(inlist != stdin)
	fclose(inlist);

      for(i=0; i < sizeinputlistvec; i++)
	free(inputlistlines[i]);
      free(inputlistlines);
      free(size_inputlistline);
    }

  if(p->Ncopycommands > 0)
    SetupLCCopies(p, c);

# ifdef PARALLEL
  /* Prepare the thread data */
  //  if(p->Nproc_allow > 1) {
    if((p->threadsinuse = (int *) malloc(p->Nproc_allow * sizeof(int))) == NULL ||
       (p->pth = (pthread_t *) malloc(p->Nproc_allow * sizeof(pthread_t))) == NULL ||
       (p->pth_init = (char *) malloc(p->Nproc_allow * sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
    for(i=0; i < p->Nproc_allow; i++)
      p->threadsinuse[i] = 0;
    pthread_mutex_init(&(p->Nproc_mutex),NULL);
    pthread_mutex_init(&(p->outfile_mutex),NULL);
    pthread_mutex_init(&(p->cfitsio_mutex),NULL);
    pthread_mutex_init(&(p->is_lc_ready_mutex),NULL);
    pthread_mutex_init(&(p->lc_getcolumnsfromheader_mutex),NULL);
    //  }
#endif

  /* If we're reading all the light curves take Nlcs = the real number of light curves, otherwise set it equal to 1 or to the number of threads. Then allocate the memory. */
  if(p->readallflag && !p->headeronly && !p->inputlistformat &&
     !p->showinputlcformat)
    Nlcs = p->Nlcs;
  else
    {
#ifdef PARALLEL
      Nlcs = p->Nproc_allow;
#else
      Nlcs = 1;
#endif
    }

  MemAllocScalarData(p, Nlcs);

  if((p->sizesinglelc = (int *) malloc(Nlcs * sizeof(int))) == NULL)
    error(ERR_MEMALLOC);

  for(i=0; i < Nlcs; i++)
    p->sizesinglelc[i] = 0;

  if(p->isifcommands) {
    if((p->IfStack = (_IfStack **) malloc(Nlcs * sizeof(_IfStack *))) == NULL)
      error(ERR_MEMALLOC);
    for(i=0; i < Nlcs; i++) {
      p->IfStack[i] = CreateIfStack();
    }
  }

  if(p->NJD == NULL) {
    if((p->NJD = (int *) malloc(Nlcs * sizeof(int))) == NULL)
      error(ERR_MEMALLOC);
  }
  for(i=0;i<Ncommands;i++)
    {
      switch(c[i].cnum)
	{
#ifdef DYNAMICLIB
	case CNUM_USERCOMMAND:
	  MemAllocDataForUserCommand(&(c[i]), Nlcs);
	  break;
#endif
	case CNUM_CLIP:
	  if((c[i].Clip->Nclip = (int *) malloc(Nlcs * sizeof(int))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case CNUM_ENSEMBLERESCALESIG:
	  if((c[i].Ensemblerescalesig->rescalefactor = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Ensemblerescalesig->chi2_old = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Ensemblerescalesig->chi2_new = (double *) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case CNUM_RESCALESIG:
	  if((c[i].Rescalesig->rescalefactor = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Rescalesig->chi2_old = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Rescalesig->chi2_new = (double *) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case CNUM_COPYLC:
	  if((c[i].CopyLC->SaveListData->lclistindx = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].CopyLC->SaveListData->runyet = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].CopyLC->SaveListData->Ndblterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].CopyLC->SaveListData->Nsterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].CopyLC->SaveListData->Nshterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].CopyLC->SaveListData->Nlterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].CopyLC->SaveListData->Niterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].CopyLC->SaveListData->Nfterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].CopyLC->SaveListData->Ncterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].CopyLC->SaveListData->dblterms = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].CopyLC->SaveListData->sterms = (char **) malloc(Nlcs * sizeof(char *))) == NULL ||
	     (c[i].CopyLC->SaveListData->iterms = (int **) malloc(Nlcs * sizeof(int *))) == NULL ||
	     (c[i].CopyLC->SaveListData->lterms = (long **) malloc(Nlcs * sizeof(long *))) == NULL ||
	     (c[i].CopyLC->SaveListData->shterms = (short **) malloc(Nlcs * sizeof(short *))) == NULL ||
	     (c[i].CopyLC->SaveListData->fterms = (float **) malloc(Nlcs * sizeof(float *))) == NULL ||
	     (c[i].CopyLC->SaveListData->cterms = (char **) malloc(Nlcs * sizeof(char *))) == NULL ||
	     (c[i].CopyLC->sizearray_IfStruct_wasfoundtrue_copy = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].CopyLC->IfStruct_wasfoundtrue_copy = (char **) malloc(Nlcs * sizeof(char *))) == NULL ||
	     (c[i].CopyLC->IfStack = (_IfStack **) malloc(Nlcs * sizeof(_IfStack *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<Nlcs;j++)
	    {
	      c[i].CopyLC->sizearray_IfStruct_wasfoundtrue_copy[j] = 0;
	      c[i].CopyLC->IfStack[j] = NULL;
	      c[i].CopyLC->SaveListData->runyet[j] = 0;
	    }
	  /* Note that the copylc command has an associated savelc command
             so we can use the same initialization code for its SaveLC 
             structure, hence no break. */
	case CNUM_SAVELC:
	  if((c[i].Savelc->NJD = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].Savelc->lclistindx = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].Savelc->dblterms = (double ***) malloc(Nlcs * sizeof(double **))) == NULL ||
	     (c[i].Savelc->sterms = (char ***) malloc(Nlcs * sizeof(char ***))) == NULL ||
	     (c[i].Savelc->iterms = (int ***) malloc(Nlcs * sizeof(int ***))) == NULL ||
	     (c[i].Savelc->lterms = (long ***) malloc(Nlcs * sizeof(long ***))) == NULL ||
	     (c[i].Savelc->shterms = (short ***) malloc(Nlcs * sizeof(short ***))) == NULL ||
	     (c[i].Savelc->fterms = (float ***) malloc(Nlcs * sizeof(float ***))) == NULL ||
	     (c[i].Savelc->cterms = (char ***) malloc(Nlcs * sizeof(char ***))) == NULL ||
	     (c[i].Savelc->Nshterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].Savelc->Nlterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].Savelc->Niterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].Savelc->Nfterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].Savelc->Ncterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].Savelc->stringid_idx = (int **) malloc(Nlcs * sizeof(int *))) == NULL ||
	     (c[i].Savelc->sizestringid_idxvecs = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
             (c[i].Savelc->sizevecs = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
             (c[i].Savelc->sizesvecs = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].Savelc->Ndblterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].Savelc->Nsterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].Savelc->runyet = (int *) malloc(Nlcs * sizeof(int))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<Nlcs;j++)
	    {
	      c[i].Savelc->runyet[j] = 0;
	      c[i].Savelc->sizevecs[j] = 0;
	      c[i].Savelc->sizesvecs[j] = 0;
	      c[i].Savelc->dblterms[j] = NULL;
	      c[i].Savelc->sterms[j] = NULL;
	      c[i].Savelc->stringid_idx[j] = NULL;
	      c[i].Savelc->iterms[j] = NULL;
	      c[i].Savelc->lterms[j] = NULL;
	      c[i].Savelc->shterms[j] = NULL;
	      c[i].Savelc->fterms[j] = NULL;
	      c[i].Savelc->cterms[j] = NULL;
	    }
	  break;
	case CNUM_RESTRICTTIMES:
	  if(c[i].RestrictTimes->restricttype == VARTOOLS_RESTRICTTIMES_JDRANGE) {
	    if(c[i].RestrictTimes->minJDtype != PERTYPE_SPECIFIED){
	      if((c[i].RestrictTimes->minJD = (double *) malloc(Nlcs * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	    if(c[i].RestrictTimes->maxJDtype != PERTYPE_SPECIFIED){
	      if((c[i].RestrictTimes->maxJD = (double *) malloc(Nlcs * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  }
	  if(c[i].RestrictTimes->saveexcludedpoints) {
	    if((c[i].RestrictTimes->s = (_Savelc *) malloc(sizeof(_Savelc))) == NULL)
	      error(ERR_MEMALLOC);
	    if((c[i].RestrictTimes->s->NJD = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	       (c[i].RestrictTimes->s->lclistindx = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	       (c[i].RestrictTimes->s->dblterms = (double ***) malloc(Nlcs * sizeof(double **))) == NULL ||
	       (c[i].RestrictTimes->s->sterms = (char ***) malloc(Nlcs * sizeof(char ***))) == NULL ||
	       (c[i].RestrictTimes->s->iterms = (int ***) malloc(Nlcs * sizeof(int ***))) == NULL ||
	       (c[i].RestrictTimes->s->lterms = (long ***) malloc(Nlcs * sizeof(long ***))) == NULL ||
	       (c[i].RestrictTimes->s->shterms = (short ***) malloc(Nlcs * sizeof(short ***))) == NULL ||
	       (c[i].RestrictTimes->s->fterms = (float ***) malloc(Nlcs * sizeof(float ***))) == NULL ||
	       (c[i].RestrictTimes->s->cterms = (char ***) malloc(Nlcs * sizeof(char ***))) == NULL ||
	       (c[i].RestrictTimes->s->Nshterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	       (c[i].RestrictTimes->s->Nlterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	       (c[i].RestrictTimes->s->Niterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	       (c[i].RestrictTimes->s->Nfterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	       (c[i].RestrictTimes->s->Ncterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	       (c[i].RestrictTimes->s->stringid_idx = (int **) malloc(Nlcs * sizeof(int *))) == NULL ||
	       (c[i].RestrictTimes->s->sizestringid_idxvecs = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	       (c[i].RestrictTimes->s->sizevecs = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	       (c[i].RestrictTimes->s->sizesvecs = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	       (c[i].RestrictTimes->s->Ndblterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	       (c[i].RestrictTimes->s->Nsterms = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	       (c[i].RestrictTimes->s->runyet = (int *) malloc(Nlcs * sizeof(int))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0;j<Nlcs;j++)
	      {
		c[i].RestrictTimes->s->runyet[j] = 0;
		c[i].RestrictTimes->s->sizevecs[j] = 0;
		c[i].RestrictTimes->s->sizesvecs[j] = 0;
		c[i].RestrictTimes->s->dblterms[j] = NULL;
		c[i].RestrictTimes->s->sterms[j] = NULL;
		c[i].RestrictTimes->s->stringid_idx[j] = NULL;
		c[i].RestrictTimes->s->iterms[j] = NULL;
		c[i].RestrictTimes->s->lterms[j] = NULL;
		c[i].RestrictTimes->s->shterms[j] = NULL;
		c[i].RestrictTimes->s->fterms[j] = NULL;
		c[i].RestrictTimes->s->cterms[j] = NULL;
	      }
	  }
	  break;
        case CNUM_CHI2_NOBIN:
	  if((c[i].Chi2_NoBin->chi2val = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Chi2_NoBin->wtave = (double *) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case CNUM_CHI2_BIN:
	  if((c[i].Chi2_Bin->chi2binvals = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Chi2_Bin->wtavebin = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<Nlcs;j++)
	    if((c[i].Chi2_Bin->chi2binvals[j] = (double *) malloc(c[i].Chi2_Bin->Nbin * sizeof(double))) == NULL ||
	       (c[i].Chi2_Bin->wtavebin[j] = (double *) malloc(c[i].Chi2_Bin->Nbin * sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	  break;
	case CNUM_CHANGEERROR:
	  if((c[i].Changeerror->rmsval = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Changeerror->ave = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Changeerror->ngood = (int *) malloc(Nlcs * sizeof(int))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case CNUM_RMS_NOBIN:
	  if((c[i].RMS_NoBin->rmsval = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].RMS_NoBin->ave = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].RMS_NoBin->rmsthy = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].RMS_NoBin->ngood = (int *) malloc(Nlcs * sizeof(int))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case CNUM_RMS_BIN:
	  if((c[i].RMS_Bin->rmsbinvals = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].RMS_Bin->rmsthybin = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<Nlcs;j++)
	    if((c[i].RMS_Bin->rmsbinvals[j] = (double *) malloc(c[i].RMS_Bin->Nbin * sizeof(double))) == NULL ||
	       (c[i].RMS_Bin->rmsthybin[j] = (double *) malloc(c[i].RMS_Bin->Nbin * sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	  break;
	case CNUM_JSTET:
	  if((c[i].Jstet->jst = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Jstet->kur = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Jstet->lst = (double *) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case CNUM_ALARM:
	  if((c[i].Alarm->alarmvals = (double *) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case CNUM_AOV:
	  if((c[i].Aov->aveaov = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Aov->rmsaov = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Aov->aveaov_whiten = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Aov->rmsaov_whiten = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Aov->peakperiods = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Aov->peakvalues = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Aov->peakSNR = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Aov->peakFAP = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<Nlcs;j++)
	    if((c[i].Aov->peakperiods[j] = (double *) malloc(c[i].Aov->Npeaks * sizeof(double))) == NULL ||
	       (c[i].Aov->peakvalues[j] = (double *) malloc(c[i].Aov->Npeaks * sizeof(double))) == NULL ||
	       (c[i].Aov->peakSNR[j] = (double *) malloc(c[i].Aov->Npeaks * sizeof(double))) == NULL ||
	       (c[i].Aov->peakFAP[j] = (double *) malloc(c[i].Aov->Npeaks * sizeof(double))) == NULL ||
	       (c[i].Aov->aveaov_whiten[j] = (double *) malloc(c[i].Aov->Npeaks * sizeof(double))) == NULL ||
	       (c[i].Aov->rmsaov_whiten[j] = (double *) malloc(c[i].Aov->Npeaks * sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	  if(c[i].Aov->fixperiodSNR)
	    {
	      if((c[i].Aov->fixperiodSNR_peakvalues = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
		 (c[i].Aov->fixperiodSNR_peakSNR = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
		 (c[i].Aov->fixperiodSNR_peakFAP = (double *) malloc(Nlcs * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	      if(c[i].Aov->fixperiodSNR_pertype != PERTYPE_SPECIFIED)
		{
		  if((c[i].Aov->fixperiodSNR_periods = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
		    error(ERR_MEMALLOC);
		  for(j=0;j<Nlcs;j++)
		    {
		      if((c[i].Aov->fixperiodSNR_periods[j] = (double *) malloc(sizeof(double))) == NULL)
			error(ERR_MEMALLOC);
		    }
		}
	    }
	  break;
	case CNUM_HARMAOV:
	  if((c[i].AovHarm->aveaov = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].AovHarm->rmsaov = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].AovHarm->aveaov_whiten = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].AovHarm->rmsaov_whiten = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].AovHarm->peakperiods = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].AovHarm->peakvalues = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].AovHarm->peakSNR = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].AovHarm->peakFAP = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].AovHarm->peakNharm = (int **) malloc(Nlcs * sizeof(int *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<Nlcs;j++)
	    if((c[i].AovHarm->peakperiods[j] = (double *) malloc(c[i].AovHarm->Npeaks * sizeof(double))) == NULL ||
	       (c[i].AovHarm->peakvalues[j] = (double *) malloc(c[i].AovHarm->Npeaks * sizeof(double))) == NULL ||
	       (c[i].AovHarm->aveaov_whiten[j] = (double *) malloc(c[i].AovHarm->Npeaks * sizeof(double))) == NULL ||
	       (c[i].AovHarm->rmsaov_whiten[j] = (double *) malloc(c[i].AovHarm->Npeaks * sizeof(double))) == NULL ||
	       (c[i].AovHarm->peakSNR[j] = (double *) malloc(c[i].AovHarm->Npeaks * sizeof(double))) == NULL ||
	       (c[i].AovHarm->peakFAP[j] = (double *) malloc(c[i].AovHarm->Npeaks * sizeof(double))) == NULL ||
	       (c[i].AovHarm->peakNharm[j] = (int *) malloc(c[i].AovHarm->Npeaks * sizeof(int))) == NULL)
	      error(ERR_MEMALLOC);
	  if(c[i].AovHarm->fixperiodSNR)
	    {
	      if((c[i].AovHarm->fixperiodSNR_peakvalues = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
		 (c[i].AovHarm->fixperiodSNR_peakSNR = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
		 (c[i].AovHarm->fixperiodSNR_peakFAP = (double *) malloc(Nlcs * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	      if(c[i].AovHarm->fixperiodSNR_pertype != PERTYPE_SPECIFIED)
		{
		  if((c[i].AovHarm->fixperiodSNR_periods = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
		    error(ERR_MEMALLOC);
		  for(j=0;j<Nlcs;j++)
		    {
		      if((c[i].AovHarm->fixperiodSNR_periods[j] = (double *) malloc(sizeof(double))) == NULL)
			error(ERR_MEMALLOC);
		    }
		}
	    }
	  break;
	case CNUM_LS:
	  if((c[i].Ls->peakperiods = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Ls->peakvalues = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Ls->peakFAP = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Ls->SNRvalues = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<Nlcs;j++)
	    if((c[i].Ls->peakperiods[j] = (double *) malloc(c[i].Ls->Npeaks * sizeof(double))) == NULL ||
	       (c[i].Ls->peakvalues[j] = (double *) malloc(c[i].Ls->Npeaks * sizeof(double))) == NULL ||
	       (c[i].Ls->peakFAP[j] = (double *) malloc(c[i].Ls->Npeaks * sizeof(double))) == NULL ||
	       (c[i].Ls->SNRvalues[j] = (double *) malloc(c[i].Ls->Npeaks * sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	  if(c[i].Ls->fixperiodSNR)
	    {
	      if((c[i].Ls->fixperiodSNR_peakvalues = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
		 (c[i].Ls->fixperiodSNR_SNRvalues = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
		 (c[i].Ls->fixperiodSNR_FAPvalues = (double *) malloc(Nlcs * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	      if(c[i].Ls->fixperiodSNR_pertype != PERTYPE_SPECIFIED)
		{
		  if((c[i].Ls->fixperiodSNR_periods = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
		    error(ERR_MEMALLOC);
		  for(j=0;j<Nlcs;j++)
		    {
		      if((c[i].Ls->fixperiodSNR_periods[j] = (double *) malloc(sizeof(double))) == NULL)
			error(ERR_MEMALLOC);
		    }
		}
	    }
	  break;
	case CNUM_DECORR:
	  if((c[i].Decorr->b = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Decorr->b_err = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Decorr->chi2val = (double *) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<Nlcs;j++)
	    if(c[i].Decorr->N_decorrterms_total > 0)
	      if((c[i].Decorr->b[j] = (double *) malloc(c[i].Decorr->N_decorrterms_total * sizeof(double))) == NULL ||
		 (c[i].Decorr->b_err[j] = (double *) malloc(c[i].Decorr->N_decorrterms_total * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	  break;
	case CNUM_GETLSAMPTHRESH:
	  if(c[i].GetLSAmpThresh->pertype != PERTYPE_SPECIFIED)
	    {
	      if((c[i].GetLSAmpThresh->period = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=0;j<Nlcs;j++)
		if((c[i].GetLSAmpThresh->period[j] = (double *) malloc(sizeof(double))) == NULL)
		  error(ERR_MEMALLOC);
	    }
	  if((c[i].GetLSAmpThresh->ampthresh_scale = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].GetLSAmpThresh->amp = (double *) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case CNUM_KILLHARM:
	  if(c[i].Killharm->pertype != PERTYPE_SPECIFIED)
	    {
	      if((c[i].Killharm->periods = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=0;j<Nlcs;j++)
		if((c[i].Killharm->periods[j] = (double *) malloc(c[i].Killharm->Nper * sizeof(double))) == NULL)
		  error(ERR_MEMALLOC);
	    }
	  if((c[i].Killharm->subharmA = (double ***) malloc(Nlcs * sizeof(double **))) == NULL ||
	     (c[i].Killharm->subharmB = (double ***) malloc(Nlcs * sizeof(double **))) == NULL ||
	     (c[i].Killharm->harmA = (double ***) malloc(Nlcs * sizeof(double **))) == NULL ||
	     (c[i].Killharm->harmB = (double ***) malloc(Nlcs * sizeof(double **))) == NULL ||
	     (c[i].Killharm->fundA = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Killharm->fundB = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Killharm->mean = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Killharm->amp = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<Nlcs;j++)
	    {
	      if((c[i].Killharm->subharmA[j] = (double **) malloc(c[i].Killharm->Nper * sizeof(double *))) == NULL ||
		 (c[i].Killharm->subharmB[j] = (double **) malloc(c[i].Killharm->Nper * sizeof(double *))) == NULL ||
		 (c[i].Killharm->harmA[j] = (double **) malloc(c[i].Killharm->Nper * sizeof(double *))) == NULL ||
		 (c[i].Killharm->harmB[j] = (double **) malloc(c[i].Killharm->Nper * sizeof(double *))) == NULL ||
		 (c[i].Killharm->fundA[j] = (double *) malloc(c[i].Killharm->Nper * sizeof(double))) == NULL ||
		 (c[i].Killharm->fundB[j] = (double *) malloc(c[i].Killharm->Nper * sizeof(double))) == NULL ||
		 (c[i].Killharm->amp[j] = (double *) malloc(c[i].Killharm->Nper * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	      for(k=0;k<c[i].Killharm->Nper;k++)
		if((c[i].Killharm->subharmA[j][k] = (double *) malloc(c[i].Killharm->Nsubharm * sizeof(double))) == NULL ||
		   (c[i].Killharm->subharmB[j][k] = (double *) malloc(c[i].Killharm->Nsubharm * sizeof(double))) == NULL ||
		   (c[i].Killharm->harmA[j][k] = (double *) malloc(c[i].Killharm->Nharm * sizeof(double))) == NULL ||
		   (c[i].Killharm->harmB[j][k] = (double *) malloc(c[i].Killharm->Nharm * sizeof(double))) == NULL)
		  error(ERR_MEMALLOC);
		}
	  break;
	case CNUM_INJECTHARM:
	  if(c[i].Injectharm->pertype != PERTYPE_SPECIFIED)
	    {
	      if((c[i].Injectharm->periods = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=0;j<Nlcs;j++)
		if((c[i].Injectharm->periods[j] = (double *) malloc(sizeof(double))) == NULL)
		  error(ERR_MEMALLOC);
	    }
	  if((c[i].Injectharm->periodinject = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Injectharm->harm_amp = (double **) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Injectharm->harm_phase = (double **) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Injectharm->subharm_amp = (double **) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Injectharm->subharm_phase = (double **) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<Nlcs;j++)
	    {
	      if((c[i].Injectharm->harm_amp[j] = (double *) malloc((c[i].Injectharm->Nharm + 1) * sizeof(double))) == NULL ||
		 (c[i].Injectharm->harm_phase[j] = (double *) malloc((c[i].Injectharm->Nharm + 1) * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	      if(c[i].Injectharm->Nsubharm > 0)
		{
		  if((c[i].Injectharm->subharm_amp[j] = (double *) malloc((c[i].Injectharm->Nsubharm) * sizeof(double))) == NULL ||
		     (c[i].Injectharm->subharm_phase[j] = (double *) malloc((c[i].Injectharm->Nsubharm) * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	    }
	  break;
	case CNUM_INJECTTRANSIT:
	  for(j=0;j<c[i].Injecttransit->Nparam;j++)
	    {
	      if((c[i].Injecttransit->paraminject[j] = (double *) malloc(Nlcs * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  break;
	case CNUM_STARSPOT:
	  if(c[i].Starspot->pertype != PERTYPE_SPECIFIED)
	    {
	      if((c[i].Starspot->period = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=0;j<Nlcs;j++)
		if((c[i].Starspot->period[j] = (double *) malloc(sizeof(double))) == NULL)
		  error(ERR_MEMALLOC);
	    }
	  if((c[i].Starspot->a = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Starspot->b = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Starspot->chi = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Starspot->inclination = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Starspot->alpha = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Starspot->psi0 = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Starspot->mconst = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Starspot->chisq = (double *) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case CNUM_BLS:
#ifdef PARALLEL
	  if((c[i].Bls->p = (double **) malloc(p->Nproc_allow*sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < p->Nproc_allow; j++) {
	    if((c[i].Bls->p[j] = (double *) malloc((c[i].Bls->nf+1)*sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	  }
#endif
	  if((c[i].Bls->sizeuv = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].Bls->u = (double **) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Bls->v = (double **) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Bls->nf2 = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].Bls->fmin = (double *) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < Nlcs; j++)
	    c[i].Bls->sizeuv[j] = 0;
	  if((c[i].Bls->bper = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Bls->bt0 = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Bls->bpow = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Bls->sde = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Bls->snval = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Bls->depth = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Bls->qtran = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Bls->chisqrplus = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Bls->chisqrminus = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Bls->bperpos = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Bls->meanmagval = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Bls->fraconenight = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Bls->i1 = (int **) malloc(Nlcs * sizeof(int *))) == NULL ||
	     (c[i].Bls->i2 = (int **) malloc(Nlcs * sizeof(int *))) == NULL ||
	     (c[i].Bls->i1_ph = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Bls->i2_ph = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Bls->nt = (int **) malloc(Nlcs * sizeof(int *))) == NULL ||
	     (c[i].Bls->Nt = (int **) malloc(Nlcs * sizeof(int *))) == NULL ||
	     (c[i].Bls->Nbefore = (int **) malloc(Nlcs * sizeof(int *))) == NULL ||
	     (c[i].Bls->Nafter = (int **) malloc(Nlcs * sizeof(int *))) == NULL ||
	     (c[i].Bls->rednoise = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Bls->whitenoise = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Bls->sigtopink = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Bls->qingress = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Bls->OOTmag = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<Nlcs;j++)
	    {
	      if((c[i].Bls->bper[j] = (double *) malloc(c[i].Bls->Npeak * sizeof(double))) == NULL ||
		 (c[i].Bls->bt0[j] = (double *) malloc(c[i].Bls->Npeak * sizeof(double))) == NULL ||
		 (c[i].Bls->bpow[j] = (double *) malloc(c[i].Bls->Npeak * sizeof(double))) == NULL ||
		 (c[i].Bls->sde[j] = (double *) malloc(c[i].Bls->Npeak * sizeof(double))) == NULL ||
		 (c[i].Bls->snval[j] = (double *) malloc(c[i].Bls->Npeak * sizeof(double))) == NULL ||
		 (c[i].Bls->depth[j] = (double *) malloc(c[i].Bls->Npeak * sizeof(double))) == NULL ||
		 (c[i].Bls->qtran[j] = (double *) malloc(c[i].Bls->Npeak * sizeof(double))) == NULL ||
		 (c[i].Bls->chisqrplus[j] = (double *) malloc(c[i].Bls->Npeak * sizeof(double))) == NULL ||
		 (c[i].Bls->fraconenight[j] = (double *) malloc(c[i].Bls->Npeak * sizeof(double))) == NULL ||
		 (c[i].Bls->i1[j] = (int *) malloc(c[i].Bls->Npeak * sizeof(int))) == NULL ||
		 (c[i].Bls->i2[j] = (int *) malloc(c[i].Bls->Npeak * sizeof(int))) == NULL ||
		 (c[i].Bls->i1_ph[j] = (double *) malloc(c[i].Bls->Npeak * sizeof(double))) == NULL ||
		 (c[i].Bls->i2_ph[j] = (double *) malloc(c[i].Bls->Npeak * sizeof(double))) == NULL ||
		 (c[i].Bls->nt[j] = (int *) malloc(c[i].Bls->Npeak * sizeof(int))) == NULL ||
		 (c[i].Bls->Nt[j] = (int *) malloc(c[i].Bls->Npeak * sizeof(int))) == NULL ||
		 (c[i].Bls->Nbefore[j] = (int *) malloc(c[i].Bls->Npeak * sizeof(int))) == NULL ||
		 (c[i].Bls->Nafter[j] = (int *) malloc(c[i].Bls->Npeak * sizeof(int))) == NULL ||
		 (c[i].Bls->rednoise[j] = (double *) malloc(c[i].Bls->Npeak * sizeof(double))) == NULL ||
		 (c[i].Bls->whitenoise[j] = (double *) malloc(c[i].Bls->Npeak * sizeof(double))) == NULL ||
		 (c[i].Bls->sigtopink[j] = (double *) malloc(c[i].Bls->Npeak * sizeof(double))) == NULL ||
		 (c[i].Bls->qingress[j] = (double *) malloc(c[i].Bls->Npeak * sizeof(double))) == NULL ||
		 (c[i].Bls->OOTmag[j] = (double *) malloc(c[i].Bls->Npeak * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  break;
	case CNUM_FINDBLENDS:

	  if((c[i].FindBlends->varx = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].FindBlends->vary = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].FindBlends->varnames = (char **) malloc(Nlcs * sizeof(char *))) == NULL ||
	     (c[i].FindBlends->varblendnames = (char **) malloc(Nlcs * sizeof(char *))) == NULL ||
	     (c[i].FindBlends->blendamps = (double *) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<Nlcs;j++)
	    {
	      if((c[i].FindBlends->varnames[j] = (char *) malloc(MAXLEN * sizeof(char))) == NULL ||
		 (c[i].FindBlends->varblendnames[j] = (char *) malloc(MAXLEN * sizeof(char))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  if(c[i].FindBlends->pertype != PERTYPE_SPECIFIED) {
	    if((c[i].FindBlends->periods = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0;j<Nlcs;j++) {
	      if((c[i].FindBlends->periods[j] = (double *) malloc(sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  }
	  break;
	case CNUM_FIXPERBLS:
	  if(c[i].BlsFixPer->pertype != PERTYPE_SPECIFIED)
	    {
	      if((c[i].BlsFixPer->period = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=0;j<Nlcs;j++)
		if((c[i].BlsFixPer->period[j] = (double *) malloc(sizeof(double))) == NULL)
		  error(ERR_MEMALLOC);
	    }
	  if((c[i].BlsFixPer->sizeuv = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].BlsFixPer->u = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].BlsFixPer->v = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j<Nlcs; j++)
	    c[i].BlsFixPer->sizeuv[j] = 0;
	  if((c[i].BlsFixPer->bpow = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixPer->bt0 = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixPer->sde = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixPer->snval = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixPer->depth = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixPer->qtran = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixPer->chisqrplus = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixPer->chisqrminus = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixPer->meanmagval = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixPer->fraconenight = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixPer->i1 = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].BlsFixPer->i2 = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].BlsFixPer->i1_ph = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixPer->i2_ph = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixPer->nt = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].BlsFixPer->Nt = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].BlsFixPer->Nbefore = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].BlsFixPer->Nafter = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].BlsFixPer->rednoise = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixPer->whitenoise = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixPer->sigtopink = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixPer->qingress = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixPer->OOTmag = (double *) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case CNUM_BLSFIXDURTC:
	  if(c[i].BlsFixDurTc->durtype != PERTYPE_SPECIFIED)
	    {
	      if((c[i].BlsFixDurTc->inputdur = (double *) malloc(Nlcs * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  if(c[i].BlsFixDurTc->TCtype != PERTYPE_SPECIFIED)
	    {
	      if((c[i].BlsFixDurTc->inputTC = (double *) malloc(Nlcs * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  if(c[i].BlsFixDurTc->fixdepth) {
	    if(c[i].BlsFixDurTc->depthtype != PERTYPE_SPECIFIED)
	      {
		if((c[i].BlsFixDurTc->inputdepth = (double *) malloc(Nlcs * sizeof(double))) == NULL)
		  error(ERR_MEMALLOC);
	      }
	    if(c[i].BlsFixDurTc->qgresstype != PERTYPE_SPECIFIED)
	      {
		if((c[i].BlsFixDurTc->inputqgress = (double *) malloc(Nlcs * sizeof(double))) == NULL)
		  error(ERR_MEMALLOC);
	      }
	  }
#ifdef PARALLEL
	  if((c[i].BlsFixDurTc->p = (double **) malloc(p->Nproc_allow*sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < p->Nproc_allow; j++) {
	    if((c[i].BlsFixDurTc->p[j] = (double *) malloc((c[i].BlsFixDurTc->nf+1)*sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	  }
#endif
	  if((c[i].BlsFixDurTc->sizeuv = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].BlsFixDurTc->u = (double **) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixDurTc->v = (double **) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixDurTc->nf2 = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	     (c[i].BlsFixDurTc->fmin = (double *) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < Nlcs; j++)
	    c[i].BlsFixDurTc->sizeuv[j] = 0;
	  if((c[i].BlsFixDurTc->bper = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].BlsFixDurTc->bt0 = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].BlsFixDurTc->bpow = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].BlsFixDurTc->sde = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].BlsFixDurTc->snval = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].BlsFixDurTc->depth = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].BlsFixDurTc->qtran = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].BlsFixDurTc->chisqrplus = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].BlsFixDurTc->chisqrminus = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixDurTc->bperpos = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixDurTc->meanmagval = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].BlsFixDurTc->fraconenight = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].BlsFixDurTc->nt = (int **) malloc(Nlcs * sizeof(int *))) == NULL ||
	     (c[i].BlsFixDurTc->Nt = (int **) malloc(Nlcs * sizeof(int *))) == NULL ||
	     (c[i].BlsFixDurTc->Nbefore = (int **) malloc(Nlcs * sizeof(int *))) == NULL ||
	     (c[i].BlsFixDurTc->Nafter = (int **) malloc(Nlcs * sizeof(int *))) == NULL ||
	     (c[i].BlsFixDurTc->rednoise = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].BlsFixDurTc->whitenoise = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].BlsFixDurTc->sigtopink = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].BlsFixDurTc->qingress = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].BlsFixDurTc->OOTmag = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<Nlcs;j++)
	    {
	      if((c[i].BlsFixDurTc->bper[j] = (double *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(double))) == NULL ||
		 (c[i].BlsFixDurTc->bt0[j] = (double *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(double))) == NULL ||
		 (c[i].BlsFixDurTc->bpow[j] = (double *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(double))) == NULL ||
		 (c[i].BlsFixDurTc->sde[j] = (double *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(double))) == NULL ||
		 (c[i].BlsFixDurTc->snval[j] = (double *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(double))) == NULL ||
		 (c[i].BlsFixDurTc->depth[j] = (double *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(double))) == NULL ||
		 (c[i].BlsFixDurTc->qtran[j] = (double *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(double))) == NULL ||
		 (c[i].BlsFixDurTc->chisqrplus[j] = (double *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(double))) == NULL ||
		 (c[i].BlsFixDurTc->fraconenight[j] = (double *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(double))) == NULL ||
		 (c[i].BlsFixDurTc->nt[j] = (int *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(int))) == NULL ||
		 (c[i].BlsFixDurTc->Nt[j] = (int *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(int))) == NULL ||
		 (c[i].BlsFixDurTc->Nbefore[j] = (int *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(int))) == NULL ||
		 (c[i].BlsFixDurTc->Nafter[j] = (int *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(int))) == NULL ||
		 (c[i].BlsFixDurTc->rednoise[j] = (double *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(double))) == NULL ||
		 (c[i].BlsFixDurTc->whitenoise[j] = (double *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(double))) == NULL ||
		 (c[i].BlsFixDurTc->sigtopink[j] = (double *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(double))) == NULL ||
		 (c[i].BlsFixDurTc->qingress[j] = (double *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(double))) == NULL ||
		 (c[i].BlsFixDurTc->OOTmag[j] = (double *) malloc(c[i].BlsFixDurTc->Npeak * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  break;
	case CNUM_SOFTENEDTRANSIT:
	  if((c[i].SoftenedTransit->period = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].SoftenedTransit->T0 = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].SoftenedTransit->eta = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].SoftenedTransit->cval = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].SoftenedTransit->delta = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].SoftenedTransit->mconst = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].SoftenedTransit->chisq = (double *) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  if(c[i].SoftenedTransit->dokillharm)
	    {
	      if((c[i].SoftenedTransit->subharmA = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
		 (c[i].SoftenedTransit->subharmB = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
		 (c[i].SoftenedTransit->harmA = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
		 (c[i].SoftenedTransit->harmB = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
		 (c[i].SoftenedTransit->fundA = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
		 (c[i].SoftenedTransit->fundB = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
		 (c[i].SoftenedTransit->per_harm_out = (double *) malloc(Nlcs * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=0;j<Nlcs;j++)
		{
		  if((c[i].SoftenedTransit->subharmA[j] = (double *) malloc(c[i].SoftenedTransit->nsubharm * sizeof(double))) == NULL ||
		     (c[i].SoftenedTransit->subharmB[j] = (double *) malloc(c[i].SoftenedTransit->nsubharm * sizeof(double))) == NULL ||
		     (c[i].SoftenedTransit->harmA[j] = (double *) malloc(c[i].SoftenedTransit->nharm * sizeof(double))) == NULL ||
		     (c[i].SoftenedTransit->harmB[j] = (double *) malloc(c[i].SoftenedTransit->nharm * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	    }
	  break;
	case CNUM_MANDELAGOLTRANSIT:
	  if((c[i].MandelAgolTransit->period = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MandelAgolTransit->T0 = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MandelAgolTransit->r = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MandelAgolTransit->a = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MandelAgolTransit->sin_i = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MandelAgolTransit->inc = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MandelAgolTransit->bimpact = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MandelAgolTransit->e = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MandelAgolTransit->omega = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MandelAgolTransit->K = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MandelAgolTransit->gamma = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MandelAgolTransit->mconst = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MandelAgolTransit->ldcoeffs = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].MandelAgolTransit->chisq = (double *) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<Nlcs;j++)
	    if((c[i].MandelAgolTransit->ldcoeffs[j] = (double *) malloc(4*sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	  break;
	case CNUM_MICROLENS:
	  if((c[i].MicroLens->f0 = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MicroLens->f1 = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MicroLens->u0 = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MicroLens->t0 = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MicroLens->tmax = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].MicroLens->chi2_ = (double *) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  break;
	case CNUM_SYSREM:
	  /* We'll initialize Sysrem here, including reading in the initial airmass terms */
	  if(!p->headeronly && !p->inputlistformat &&
	     !p->showinputlcformat)
	    initialize_sysrem(c[i].Sysrem,Nlcs,p->matchstringid);
	  break;
	case CNUM_TFA:
	  /* Note we'll initialize the tfa here, including reading in the trend light curves and inverting the tfa design matrix - this will be slow depending on the number of trends included to filter */
	  if(!p->headeronly && !p->inputlistformat &&
	     !p->showinputlcformat)
	    {
	      initialize_tfa(c[i].TFA,p);
	      if((c[i].TFA->ave_out = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
		 (c[i].TFA->rms_out = (double *) malloc(Nlcs * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  break;

	case CNUM_TFA_SR:
	  /* Note we'll initialize the tfa here, including reading in the trend light curves and inverting the tfa design matrix - this will be slow depending on the number of trends included to filter */
	  if(!p->headeronly && !p->inputlistformat &&
	     !p->showinputlcformat)
	    {
	      initialize_tfa_sr(c[i].TFA_SR, Nlcs, p);
	      if((c[i].TFA_SR->ave_out = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
		 (c[i].TFA_SR->rms_out = (double *) malloc(Nlcs * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	      if((c[i].TFA_SR->use_bin || c[i].TFA_SR->use_harm) && c[i].TFA_SR->use_period && c[i].TFA_SR->pertype != PERTYPE_SPECIFIED && c[i].TFA_SR->pertype != PERTYPE_FIX)
		{
		  if((c[i].TFA_SR->periods = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
		    error(ERR_MEMALLOC);
		  for(j=0;j<Nlcs;j++)
		    if((c[i].TFA_SR->periods[j] = (double *) malloc(sizeof(double))) == NULL)
		      error(ERR_MEMALLOC);
		}
	    }
	  break;

	case CNUM_PHASE:
	  if(c[i].Phase->pertype != PERTYPE_SPECIFIED)
	    {
	      if((c[i].Phase->period = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=0;j<Nlcs;j++)
		if((c[i].Phase->period[j] = (double *) malloc(sizeof(double))) == NULL)
		  error(ERR_MEMALLOC);
	    }
	  break;
	case CNUM_DFTCLEAN:
	  if((c[i].Dftclean->peakfreqs_dirty = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Dftclean->peakpows_dirty = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Dftclean->SNR_dirty = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Dftclean->peakfreqs_clean = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Dftclean->peakpows_clean = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Dftclean->SNR_clean = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Dftclean->aveper_dirty = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Dftclean->stdper_dirty = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Dftclean->aveper_noclip_dirty = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Dftclean->stdper_noclip_dirty = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Dftclean->aveper_clean = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Dftclean->stdper_clean = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Dftclean->aveper_noclip_clean = (double *) malloc(Nlcs * sizeof(double))) == NULL ||
	     (c[i].Dftclean->stdper_noclip_clean = (double *) malloc(Nlcs * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  if(c[i].Dftclean->Npeaks_dirty > 0)
	    {
	      for(j=0;j<Nlcs;j++)
		{
		  if((c[i].Dftclean->peakfreqs_dirty[j] = (double *) malloc(c[i].Dftclean->Npeaks_dirty * sizeof(double))) == NULL ||
		     (c[i].Dftclean->peakpows_dirty[j] = (double *) malloc(c[i].Dftclean->Npeaks_dirty * sizeof(double))) == NULL ||
		     (c[i].Dftclean->SNR_dirty[j] = (double *) malloc(c[i].Dftclean->Npeaks_dirty * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	    }
	  if(c[i].Dftclean->Npeaks_clean > 0)
	    {
	      for(j=0; j<Nlcs;j++)
		{
		  if((c[i].Dftclean->peakfreqs_clean[j] = (double *) malloc(c[i].Dftclean->Npeaks_clean * sizeof(double))) == NULL ||
		     (c[i].Dftclean->peakpows_clean[j] = (double *) malloc(c[i].Dftclean->Npeaks_clean * sizeof(double))) == NULL ||
		     (c[i].Dftclean->SNR_clean[j] = (double *) malloc(c[i].Dftclean->Npeaks_clean * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	    }
	  break;
	case CNUM_LINFIT:
	  if((c[i].Linfit->param_outvals = (double **) malloc(Nlcs * sizeof(double *))) == NULL ||
	     (c[i].Linfit->param_uncertainties = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j = 0; j < Nlcs; j++)
	    {
	      if((c[i].Linfit->param_outvals[j] = (double *) malloc(c[i].Linfit->Nparams * sizeof(double))) == NULL ||
		 (c[i].Linfit->param_uncertainties[j] = (double *) malloc(c[i].Linfit->Nparams * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  break;
	case CNUM_NONLINFIT:
	  if(c[i].Nonlinfit->use_covar) {
	    if((c[i].Nonlinfit->Corr_sizemat = (int *) malloc(Nlcs * sizeof(int))) == NULL ||
	       (c[i].Nonlinfit->Corr_mat1 = (double ***) malloc(Nlcs * sizeof(double **))) == NULL ||
	       (c[i].Nonlinfit->Corr_mat2 = (double ***) malloc(Nlcs * sizeof(double **))) == NULL ||
	       (c[i].Nonlinfit->Corr_Nvec = (int **) malloc(Nlcs * sizeof(int *))) == NULL ||
	       (c[i].Nonlinfit->Corr_store_NJD = (int *) malloc(Nlcs * sizeof(int))) == NULL) 
	      error(ERR_MEMALLOC);
	    for(j=0; j < Nlcs; j++) {
	      c[i].Nonlinfit->Corr_sizemat[j] = 0;
	      c[i].Nonlinfit->Corr_store_NJD[j] = 0;
	      c[i].Nonlinfit->Corr_Nvec[j] = NULL;
	      c[i].Nonlinfit->Corr_mat1[j] = NULL;
	      c[i].Nonlinfit->Corr_mat2[j] = NULL;
	    }
	  }
	  break;
	case CNUM_WWZ:
	  break;
	case CNUM_STATS:
	  if((c[i].Stats->statsout = (double **) malloc(Nlcs * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0; j < Nlcs; j++)
	    {
	      if((c[i].Stats->statsout[j] = (double *) malloc(c[i].Stats->Nstatstot * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  break;
	case CNUM_IF:
	  if(c[i].IfCommand->iftype == VARTOOLS_IFTYPE_IF) {
	    if((c[i].IfCommand->ifs->wasfoundtrue = (char *) malloc(Nlcs * sizeof(char))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  break;
#ifdef _HAVE_PYTHON
	case CNUM_PYTHON:
	  InitPythonCommand(p, c[i].PythonCommand, Nlcs);
	  break;
#endif
	default:
	  break;
	}
    }



}
