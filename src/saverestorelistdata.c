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
/* This file contains functions to save and restore list data for the VARTOOLS light curve analysis program written by J. Hartman. */

#include "commands.h"
#include "programdata.h"
#include "functions.h"

void dosavelistdata(ProgramData *p, _SaveListData *s, int threadid, int lcid)
{
  int i, j, idbl, istring, iint, ilong, ishort, ichar, ifloat, k, u, idx;
  _Variable *v;
  _DataFromInputList *d;
  int Nc;
  void *lcnamesptr;

  lcnamesptr = (void *) (&p->lcnames);

  s->lclistindx[threadid] = lcid;
  if(s->runyet[threadid] == 0)
    {
      s->Ndblterms[threadid] = 0;
      s->Nsterms[threadid] = 0;
      s->Nshterms[threadid] = 0;
      s->Nlterms[threadid] = 0;
      s->Niterms[threadid] = 0;
      s->Nfterms[threadid] = 0;
      s->Ncterms[threadid] = 0;

      for(i=0; i < p->NDefinedVariables; i++) {
	if(p->DefinedVariables[i]->vectortype == VARTOOLS_VECTORTYPE_INLIST ||
	   p->DefinedVariables[i]->vectortype == VARTOOLS_VECTORTYPE_OUTCOLUMN)
	  {
	    v = p->DefinedVariables[i];
	    if(v->dataptr == lcnamesptr)
	      continue;
	    switch(v->datatype)
	      {
	      case VARTOOLS_TYPE_DOUBLE:
		s->Ndblterms[threadid]++;
		break;
	      case VARTOOLS_TYPE_CONVERTJD:
		s->Ndblterms[threadid]++;
		break;
	      case VARTOOLS_TYPE_STRING:
		s->Nsterms[threadid] += MAXLEN;
		break;
	      case VARTOOLS_TYPE_INT:
		s->Niterms[threadid]++;
		break;
	      case VARTOOLS_TYPE_SHORT:
		s->Nshterms[threadid]++;
		break;
	      case VARTOOLS_TYPE_FLOAT:
		s->Nfterms[threadid]++;
		break;
	      case VARTOOLS_TYPE_LONG:
		s->Nlterms[threadid]++;
		break;
	      case VARTOOLS_TYPE_CHAR:
		s->Ncterms[threadid]++;
		break;
	      }
	  }
      }
      for(i=0;i<p->NDataFromInputList;i++)
	{
	  d = &(p->DataFromInputList[i]);
	  if(d->dataptr == lcnamesptr)
	    continue;
	  Nc = d->Ncolumns;
	  switch(d->datatype)
	    {
	    case VARTOOLS_TYPE_DOUBLE:
	      if(Nc == 0)
		s->Ndblterms[threadid]++;
	      else
		s->Ndblterms[threadid] += abs(Nc);
	      break;
	    case VARTOOLS_TYPE_CONVERTJD:
	      if(Nc == 0)
		s->Ndblterms[threadid]++;
	      else
		s->Ndblterms[threadid] += abs(Nc);
	      break;
	    case VARTOOLS_TYPE_STRING:
	      if(Nc == 0)
		s->Nsterms[threadid] += MAXLEN;
	      else
		s->Nsterms[threadid] += abs(Nc)*MAXLEN;
	      break;
	    case VARTOOLS_TYPE_INT:
	      if(Nc == 0)
		s->Niterms[threadid]++;
	      else
		s->Niterms[threadid] += abs(Nc);
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      if(Nc == 0)
		s->Nshterms[threadid]++;
	      else
		s->Nshterms[threadid] += abs(Nc);
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      if(Nc == 0)
		s->Nfterms[threadid]++;
	      else
		s->Nfterms[threadid] += abs(Nc);
	      break;
	    case VARTOOLS_TYPE_LONG:
	      if(Nc == 0)
		s->Nlterms[threadid]++;
	      else
		s->Nlterms[threadid] += abs(Nc);
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      if(Nc == 0)
		s->Ncterms[threadid]++;
	      else
		s->Ncterms[threadid] += abs(Nc);
	      break;
	    }
	}
      /* Allocate space to store the copies */
      if(s->Ndblterms[threadid] > 0)
	{
	  if((s->dblterms[threadid] = (double *) malloc(s->Ndblterms[threadid] * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      if(s->Nsterms[threadid] > 0)
	{
	  if((s->sterms[threadid] = (char *) malloc(s->Nsterms[threadid])) == NULL)
	    error(ERR_MEMALLOC);
	}
      if(s->Niterms[threadid] > 0)
	{
	  if((s->iterms[threadid] = (int *) malloc(s->Niterms[threadid] * sizeof(int))) == NULL)
	    error(ERR_MEMALLOC);
	}
      if(s->Nlterms[threadid] > 0)
	{
	  if((s->lterms[threadid] = (long *) malloc(s->Nlterms[threadid] * sizeof(long))) == NULL)
	    error(ERR_MEMALLOC);
	}
      if(s->Nshterms[threadid] > 0)
	{
	  if((s->shterms[threadid] = (short *) malloc(s->Nshterms[threadid] * sizeof(short))) == NULL)
	    error(ERR_MEMALLOC);
	}
      if(s->Nfterms[threadid] > 0)
	{
	  if((s->fterms[threadid] = (float *) malloc(s->Nfterms[threadid] * sizeof(float))) == NULL)
	    error(ERR_MEMALLOC);
	}
      if(s->Ncterms[threadid] > 0)
	{
	  if((s->cterms[threadid] = (char *) malloc(s->Ncterms[threadid] * sizeof(char))) == NULL)
	    error(ERR_MEMALLOC);
	}

      s->runyet[threadid] = 1;
    }
  
  /* copy the data */
  for(idbl=0,istring=0,ishort=0,ilong=0,ifloat=0,ichar=0,iint=0,i=0;
      i<p->NDefinedVariables;i++)
    {
      if(p->DefinedVariables[i]->vectortype == VARTOOLS_VECTORTYPE_INLIST) {
	v = p->DefinedVariables[i];
	if(v->dataptr == lcnamesptr)
	  continue;
	switch(v->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  s->dblterms[threadid][idbl] = (*((double **) v->dataptr))[lcid];
	  idbl++;
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  s->dblterms[threadid][idbl] = (*((double **) v->dataptr))[lcid];
	  idbl++;
	  break;
	case VARTOOLS_TYPE_STRING:
	  memcpy(&(s->sterms[threadid][istring]),(*((char ***) v->dataptr))[lcid],MAXLEN);
	  istring += MAXLEN;
	  break;
	case VARTOOLS_TYPE_INT:
	  s->iterms[threadid][iint] = (*((int **) v->dataptr))[lcid];
	  iint++;
	  break;
	case VARTOOLS_TYPE_LONG:
	  s->lterms[threadid][ilong] = (*((long **) v->dataptr))[lcid];
	  ilong++;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  s->fterms[threadid][ifloat] = (*((float **) v->dataptr))[lcid];
	  ifloat++;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  s->shterms[threadid][ishort] = (*((short **) v->dataptr))[lcid];
	  ishort++;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  s->cterms[threadid][ichar] = (*((char **) v->dataptr))[lcid];
	  ichar++;
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      }
      else if (p->DefinedVariables[i]->vectortype == VARTOOLS_VECTORTYPE_OUTCOLUMN) {
	v = p->DefinedVariables[i];
	if(v->dataptr == lcnamesptr)
	  continue;
	switch(v->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  getoutcolumnvalue(v->outc, threadid, lcid, VARTOOLS_TYPE_DOUBLE, &s->dblterms[threadid][idbl]);
	  idbl++;
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  getoutcolumnvalue(v->outc, threadid, lcid, VARTOOLS_TYPE_DOUBLE, &s->dblterms[threadid][idbl]);
	  idbl++;
	  break;
	case VARTOOLS_TYPE_STRING:
	  getoutcolumnvalue(v->outc, threadid, lcid, VARTOOLS_TYPE_STRING, &s->sterms[threadid][istring]);
	  istring += MAXLEN;
	  break;
	case VARTOOLS_TYPE_INT:
	  getoutcolumnvalue(v->outc, threadid, lcid, VARTOOLS_TYPE_INT, &s->iterms[threadid][iint]);
	  iint++;
	  break;
	case VARTOOLS_TYPE_LONG:
	  getoutcolumnvalue(v->outc, threadid, lcid, VARTOOLS_TYPE_LONG, &s->lterms[threadid][ilong]);
	  ilong++;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  getoutcolumnvalue(v->outc, threadid, lcid, VARTOOLS_TYPE_FLOAT, &s->fterms[threadid][ifloat]);
	  ifloat++;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  getoutcolumnvalue(v->outc, threadid, lcid, VARTOOLS_TYPE_SHORT, &s->shterms[threadid][ishort]);
	  ishort++;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  getoutcolumnvalue(v->outc, threadid, lcid, VARTOOLS_TYPE_CHAR, &s->cterms[threadid][ichar]);
	  ichar++;
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      }
    }

  for(idbl=0,istring=0,ishort=0,ilong=0,ifloat=0,ichar=0,iint=0,i=0;
      i<p->NDataFromInputList;i++)
    {
      d = &(p->DataFromInputList[i]);
      if(d->dataptr == lcnamesptr)
	continue;
      Nc = d->Ncolumns;
      if(Nc == 0) {
	switch(d->datatype)
	  {
	  case VARTOOLS_TYPE_DOUBLE:
	    s->dblterms[threadid][idbl] = (*((double **) d->dataptr))[lcid];
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    s->dblterms[threadid][idbl] = (*((double **) d->dataptr))[lcid];
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    memcpy(&(s->sterms[threadid][istring]),((*((char ***)d->dataptr))[lcid]),MAXLEN);
	    istring += MAXLEN;
	    break;
	  case VARTOOLS_TYPE_INT:
	    s->iterms[threadid][iint] = (*((int **) d->dataptr))[lcid];
	    iint++;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    s->lterms[threadid][ilong] = (*((long **) d->dataptr))[lcid];
	    ilong++;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    s->fterms[threadid][ifloat] = (*((float **) d->dataptr))[lcid];
	    ifloat++;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    s->shterms[threadid][ishort] = (*((short **) d->dataptr))[lcid];
	    ishort++;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    s->cterms[threadid][ichar] = (*((char **) d->dataptr))[lcid];
	    ichar++;
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
      } else if(Nc != 0) {
	for(u=0; u < abs(Nc); u++) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    s->dblterms[threadid][idbl] = (*((double ***) d->dataptr))[lcid][u];
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    s->dblterms[threadid][idbl] = (*((double ***) d->dataptr))[lcid][u];
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    memcpy(&(s->sterms[threadid][istring]),((*((char ****)d->dataptr))[lcid][u]),MAXLEN);
	    istring += MAXLEN;
	    break;
	  case VARTOOLS_TYPE_INT:
	    s->iterms[threadid][iint] = (*((int ***) d->dataptr))[lcid][u];
	    iint++;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    s->lterms[threadid][ilong] = (*((long ***) d->dataptr))[lcid][u];
	    ilong++;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    s->fterms[threadid][ifloat] = (*((float ***) d->dataptr))[lcid][u];
	    ifloat++;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    s->shterms[threadid][ishort] = (*((short ***) d->dataptr))[lcid][u];
	    ishort++;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    s->cterms[threadid][ichar] = (*((char ***) d->dataptr))[lcid][u];
	    ichar++;
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      }
    }
}


void dorestorelistdata(ProgramData *p, _SaveListData *s, int sthreadid, int rthreadid, int lcid)
{
  int i, j, idbl, istring, iint, ilong, ishort, ichar, ifloat, k, u, idx;
  _Variable *v;
  _DataFromInputList *d;
  int Nc;
  
  void *lcnamesptr;
  
  lcnamesptr = (void *) (&p->lcnames);

  /* copy the data */
  for(idbl=0,istring=0,ishort=0,ilong=0,ifloat=0,ichar=0,iint=0,i=0;
      i<p->NDefinedVariables;i++)
    {
      if(p->DefinedVariables[i]->vectortype == VARTOOLS_VECTORTYPE_INLIST) {
	v = p->DefinedVariables[i];
	if(v->dataptr == lcnamesptr)
	  continue;
	switch(v->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  (*((double **) v->dataptr))[lcid] = s->dblterms[sthreadid][idbl];
	  idbl++;
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  (*((double **) v->dataptr))[lcid] = s->dblterms[sthreadid][idbl];
	  idbl++;
	  break;
	case VARTOOLS_TYPE_STRING:
	  memcpy((*((char ***) v->dataptr))[lcid],&(s->sterms[sthreadid][istring]),MAXLEN);
	  istring += MAXLEN;
	  break;
	case VARTOOLS_TYPE_INT:
	  (*((int **) v->dataptr))[lcid] = s->iterms[sthreadid][iint];
	  iint++;
	  break;
	case VARTOOLS_TYPE_LONG:
	  (*((long **) v->dataptr))[lcid] = s->lterms[sthreadid][ilong];
	  ilong++;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  (*((float **) v->dataptr))[lcid] = s->fterms[sthreadid][ifloat];
	  ifloat++;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  (*((short **) v->dataptr))[lcid] = s->shterms[sthreadid][ishort];
	  ishort++;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  (*((char **) v->dataptr))[lcid] = s->cterms[sthreadid][ichar];
	  ichar++;
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      }
      else if (p->DefinedVariables[i]->vectortype == VARTOOLS_VECTORTYPE_OUTCOLUMN) {
	v = p->DefinedVariables[i];
	if(v->dataptr == lcnamesptr)
	  continue;
	switch(v->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  setoutcolumnvalue(v->outc, rthreadid, lcid, VARTOOLS_TYPE_DOUBLE, &s->dblterms[sthreadid][idbl]);
	  idbl++;
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  setoutcolumnvalue(v->outc, rthreadid, lcid, VARTOOLS_TYPE_DOUBLE, &s->dblterms[sthreadid][idbl]);
	  idbl++;
	  break;
	case VARTOOLS_TYPE_STRING:
	  setoutcolumnvalue(v->outc, rthreadid, lcid, VARTOOLS_TYPE_STRING, &s->sterms[sthreadid][istring]);
	  istring += MAXLEN;
	  break;
	case VARTOOLS_TYPE_INT:
	  setoutcolumnvalue(v->outc, rthreadid, lcid, VARTOOLS_TYPE_INT, &s->iterms[sthreadid][iint]);
	  iint++;
	  break;
	case VARTOOLS_TYPE_LONG:
	  setoutcolumnvalue(v->outc, rthreadid, lcid, VARTOOLS_TYPE_LONG, &s->lterms[sthreadid][ilong]);
	  ilong++;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  setoutcolumnvalue(v->outc, rthreadid, lcid, VARTOOLS_TYPE_FLOAT, &s->fterms[sthreadid][ifloat]);
	  ifloat++;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  setoutcolumnvalue(v->outc, rthreadid, lcid, VARTOOLS_TYPE_SHORT, &s->shterms[sthreadid][ishort]);
	  ishort++;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  setoutcolumnvalue(v->outc, rthreadid, lcid, VARTOOLS_TYPE_CHAR, &s->cterms[sthreadid][ichar]);
	  ichar++;
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      }
    }

  for(idbl=0,istring=0,ishort=0,ilong=0,ifloat=0,ichar=0,iint=0,i=0;
      i<p->NDataFromInputList;i++)
    {
      d = &(p->DataFromInputList[i]);
      if(d->dataptr == lcnamesptr)
	continue;
      Nc = d->Ncolumns;
      if(Nc == 0) {
	switch(d->datatype)
	  {
	  case VARTOOLS_TYPE_DOUBLE:
	    (*((double **) d->dataptr))[lcid] = s->dblterms[sthreadid][idbl];
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    (*((double **) d->dataptr))[lcid] = s->dblterms[sthreadid][idbl];
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    memcpy(((*((char ***)d->dataptr))[lcid]),&(s->sterms[sthreadid][istring]),MAXLEN);
	    istring += MAXLEN;
	    break;
	  case VARTOOLS_TYPE_INT:
	    (*((int **) d->dataptr))[lcid] = s->iterms[sthreadid][iint];
	    iint++;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    (*((long **) d->dataptr))[lcid] = s->lterms[sthreadid][ilong];
	    ilong++;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    (*((float **) d->dataptr))[lcid] = s->fterms[sthreadid][ifloat];
	    ifloat++;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    (*((short **) d->dataptr))[lcid] = s->shterms[sthreadid][ishort];
	    ishort++;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    (*((char **) d->dataptr))[lcid] = s->cterms[sthreadid][ichar];
	    ichar++;
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
      } else if(Nc != 0) {
	for(u=0; u < abs(Nc); u++) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    (*((double ***) d->dataptr))[lcid][u] = s->dblterms[sthreadid][idbl];
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    (*((double ***) d->dataptr))[lcid][u] = s->dblterms[sthreadid][idbl];
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    memcpy(((*((char ****)d->dataptr))[lcid][u]),&(s->sterms[sthreadid][istring]),MAXLEN);
	    istring += MAXLEN;
	    break;
	  case VARTOOLS_TYPE_INT:
	    (*((int ***) d->dataptr))[lcid][u] = s->iterms[sthreadid][iint];
	    iint++;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    (*((long ***) d->dataptr))[lcid][u] = s->lterms[sthreadid][ilong];
	    ilong++;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    (*((float ***) d->dataptr))[lcid][u] = s->fterms[sthreadid][ifloat];
	    ifloat++;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    (*((short ***) d->dataptr))[lcid][u] = s->shterms[sthreadid][ishort];
	    ishort++;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    (*((char ***) d->dataptr))[lcid][u] = s->cterms[sthreadid][ichar];
	    ichar++;
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      }
    }
}
