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
/* This file contains functions to save and restore a light curve for the VARTOOLS light curve analysis program written by J. Hartman. */

#include "commands.h"
#include "programdata.h"
#include "functions.h"

void dosavelc(ProgramData *p, _Savelc *s, int threadid, int lcid)
{
  int i, j, idbl, istring, iint, ilong, ishort, ichar, ifloat, k, u;
  _DataFromLightCurve *d;
  int Nc;

  s->NJD[threadid] = p->NJD[threadid];
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
      s->maxstring = 0;
      for(i=0;i<p->NDataFromLightCurve;i++)
	{
	  d = &(p->DataFromLightCurve[i]);
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
		s->Nsterms[threadid] += d->maxstringlength;
	      else
		s->Nsterms[threadid] += abs(Nc)*d->maxstringlength;
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
      s->runyet[threadid] = 1;
    }

  /* Make sure we have enough space to store the copy of the light curve */
  if(s->Ndblterms[threadid] > 0)
    {
      if(s->dblterms[threadid] == NULL)
	{
	  if((s->dblterms[threadid] = (double **) malloc(s->Ndblterms[threadid] * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  if(p->NJD[threadid] > 0)
	    {
	      for(i=0;i<s->Ndblterms[threadid];i++)
		{
		  if((s->dblterms[threadid][i] = (double *) malloc(p->NJD[threadid] * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	}
      else
	{
	  if(p->NJD[threadid] > 0 && s->sizevecs[threadid] == 0)
	    {
	      for(i=0;i<s->Ndblterms[threadid];i++)
		{
		  if((s->dblterms[threadid][i] = (double *) malloc(p->NJD[threadid] * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	  else if(p->NJD[threadid] > s->sizevecs[threadid])
	    {
	      for(i=0;i<s->Ndblterms[threadid];i++)
		{
		  if((s->dblterms[threadid][i] = (double *) realloc(s->dblterms[threadid][i],p->NJD[threadid] * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	}
    }
  if(s->Nsterms[threadid] > 0)
    {
      if(s->sterms[threadid] == NULL)
	{
	  if(p->NJD[threadid] > 0)
	    {
	      if((s->sterms[threadid] = (char **) malloc(p->NJD[threadid] * sizeof(char *))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=0;j<p->NJD[threadid];j++)
		if((s->sterms[threadid][j] = (char *) malloc(s->Nsterms[threadid])) == NULL)
		  error(ERR_MEMALLOC);
	    }
	  s->sizesvecs[threadid] = p->NJD[threadid];
	}
      else
	{
	  if(p->NJD[threadid] > s->sizesvecs[threadid])
	    {
	      if((s->sterms[threadid] = (char **) realloc(s->sterms[threadid],p->NJD[threadid] * sizeof(char *))) == NULL)
		error(ERR_MEMALLOC);
	      for(j=s->sizesvecs[threadid];j<p->NJD[threadid];j++)
		{
		  if((s->sterms[threadid][j] = (char *) malloc(s->Nsterms[threadid])) == NULL)
		    error(ERR_MEMALLOC);
		}
	    }
	  s->sizesvecs[threadid] = p->NJD[threadid];
	}
    }
  if(s->Niterms[threadid] > 0)
    {
      if(s->iterms[threadid] == NULL)
	{
	  if((s->iterms[threadid] = (int **) malloc(s->Niterms[threadid] * sizeof(int *))) == NULL)
	    error(ERR_MEMALLOC);
	  if(p->NJD[threadid] > 0)
	    {
	      for(i=0;i<s->Niterms[threadid];i++)
		{
		  if((s->iterms[threadid][i] = (int *) malloc(p->NJD[threadid] * sizeof(int))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	}
      else
	{
	  if(p->NJD[threadid] > 0 && s->sizevecs[threadid] == 0)
	    {
	      for(i=0;i<s->Niterms[threadid];i++)
		{
		  if((s->iterms[threadid][i] = (int *) malloc(p->NJD[threadid] * sizeof(int))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	  else if(p->NJD[threadid] > s->sizevecs[threadid])
	    {
	      for(i=0;i<s->Niterms[threadid];i++)
		{
		  if((s->iterms[threadid][i] = (int *) realloc(s->iterms[threadid][i],p->NJD[threadid] * sizeof(int))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	}
    }
  if(s->Nlterms[threadid] > 0)
    {
      if(s->lterms[threadid] == NULL)
	{
	  if((s->lterms[threadid] = (long **) malloc(s->Nlterms[threadid] * sizeof(long *))) == NULL)
	    error(ERR_MEMALLOC);
	  if(p->NJD[threadid] > 0)
	    {
	      for(i=0;i<s->Nlterms[threadid];i++)
		{
		  if((s->lterms[threadid][i] = (long *) malloc(p->NJD[threadid] * sizeof(long))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	}
      else
	{
	  if(p->NJD[threadid] > 0 && s->sizevecs[threadid] == 0)
	    {
	      for(i=0;i<s->Nlterms[threadid];i++)
		{
		  if((s->lterms[threadid][i] = (long *) malloc(p->NJD[threadid] * sizeof(long))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	  else if(p->NJD[threadid] > s->sizevecs[threadid])
	    {
	      for(i=0;i<s->Nlterms[threadid];i++)
		{
		  if((s->lterms[threadid][i] = (long *) realloc(s->lterms[threadid][i],p->NJD[threadid] * sizeof(long))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	}
    }
  if(s->Nshterms[threadid] > 0)
    {
      if(s->shterms[threadid] == NULL)
	{
	  if((s->shterms[threadid] = (short **) malloc(s->Nshterms[threadid] * sizeof(short *))) == NULL)
	    error(ERR_MEMALLOC);
	  if(p->NJD[threadid] > 0)
	    {
	      for(i=0;i<s->Nshterms[threadid];i++)
		{
		  if((s->shterms[threadid][i] = (short *) malloc(p->NJD[threadid] * sizeof(short))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	}
      else
	{
	  if(p->NJD[threadid] > 0 && s->sizevecs[threadid] == 0)
	    {
	      for(i=0;i<s->Nshterms[threadid];i++)
		{
		  if((s->shterms[threadid][i] = (short *) malloc(p->NJD[threadid] * sizeof(short))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	  else if(p->NJD[threadid] > s->sizevecs[threadid])
	    {
	      for(i=0;i<s->Nshterms[threadid];i++)
		{
		  if((s->shterms[threadid][i] = (short *) realloc(s->shterms[threadid][i],p->NJD[threadid] * sizeof(short))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	}
    }
  if(s->Nfterms[threadid] > 0)
    {
      if(s->fterms[threadid] == NULL)
	{
	  if((s->fterms[threadid] = (float **) malloc(s->Nfterms[threadid] * sizeof(float *))) == NULL)
	    error(ERR_MEMALLOC);
	  if(p->NJD[threadid] > 0)
	    {
	      for(i=0;i<s->Nfterms[threadid];i++)
		{
		  if((s->fterms[threadid][i] = (float *) malloc(p->NJD[threadid] * sizeof(float))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	}
      else
	{
	  if(p->NJD[threadid] > 0 && s->sizevecs[threadid] == 0)
	    {
	      for(i=0;i<s->Nfterms[threadid];i++)
		{
		  if((s->fterms[threadid][i] = (float *) malloc(p->NJD[threadid] * sizeof(float))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	  else if(p->NJD[threadid] > s->sizevecs[threadid])
	    {
	      for(i=0;i<s->Nfterms[threadid];i++)
		{
		  if((s->fterms[threadid][i] = (float *) realloc(s->fterms[threadid][i],p->NJD[threadid] * sizeof(float))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	}
    }
  if(s->Ncterms[threadid] > 0)
    {
      if(s->cterms[threadid] == NULL)
	{
	  if((s->cterms[threadid] = (char **) malloc(s->Ncterms[threadid] * sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	  if(p->NJD[threadid] > 0)
	    {
	      for(i=0;i<s->Ncterms[threadid];i++)
		{
		  if((s->cterms[threadid][i] = (char *) malloc(p->NJD[threadid] * sizeof(char))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	}
      else
	{
	  if(p->NJD[threadid] > 0 && s->sizevecs[threadid] == 0)
	    {
	      for(i=0;i<s->Ncterms[threadid];i++)
		{
		  if((s->cterms[threadid][i] = (char *) malloc(p->NJD[threadid] * sizeof(char))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	  else if(p->NJD[threadid] > s->sizevecs[threadid])
	    {
	      for(i=0;i<s->Ncterms[threadid];i++)
		{
		  if((s->cterms[threadid][i] = (char *) realloc(s->cterms[threadid][i],p->NJD[threadid] * sizeof(char))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      s->sizevecs[threadid] = p->NJD[threadid];
	    }
	}
    }

  if(p->readimagestring)
    {
      if(s->stringid_idx[threadid] == NULL)
	{
	  if(p->NJD[threadid] > 0)
	    {
	      if((s->stringid_idx[threadid] = (int *) malloc(p->NJD[threadid] * sizeof(int))) == NULL)
		error(ERR_MEMALLOC);
	      s->sizestringid_idxvecs[threadid] = p->NJD[threadid];
	    }
	}
      else if(p->NJD[threadid] > s->sizestringid_idxvecs[threadid])
	{
	  if((s->stringid_idx[threadid] = (int *) realloc(s->stringid_idx[threadid], p->NJD[threadid] * sizeof(int))) == NULL)
	    error(ERR_MEMALLOC);
	  s->sizestringid_idxvecs[threadid] = p->NJD[threadid];
	}
    }
  s->NJD[threadid] = p->NJD[threadid];

  /* copy the light curve */
  for(idbl=0,istring=0,ishort=0,ilong=0,ifloat=0,ichar=0,iint=0,i=0;
      i<p->NDataFromLightCurve;i++)
    {
      d = &(p->DataFromLightCurve[i]);
      Nc = d->Ncolumns;
      if(Nc == 0) {
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  memcpy(s->dblterms[threadid][idbl],&((*((double ***)d->dataptr))[threadid][0]),(p->NJD[threadid]*sizeof(double)));
	  idbl++;
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  memcpy(s->dblterms[threadid][idbl],&((*((double ***)d->dataptr))[threadid][0]),(p->NJD[threadid]*sizeof(double)));
	  idbl++;
	  break;
	case VARTOOLS_TYPE_STRING:
	  for(k=0;k<p->NJD[threadid];k++) {
	    memcpy(&(s->sterms[threadid][k][istring]),((*((char ****)d->dataptr))[threadid][k]),d->maxstringlength);
	  }
	  istring += d->maxstringlength;
	  break;
	case VARTOOLS_TYPE_INT:
	  memcpy(s->iterms[threadid][iint],&((*((int ***)d->dataptr))[threadid][0]),(p->NJD[threadid]*sizeof(int)));
	  iint++;
	  break;
	case VARTOOLS_TYPE_LONG:
	  memcpy(s->lterms[threadid][ilong],&((*((long ***)d->dataptr))[threadid][0]),(p->NJD[threadid]*sizeof(long)));
	  ilong++;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  memcpy(s->fterms[threadid][ifloat],&((*((float ***)d->dataptr))[threadid][0]),(p->NJD[threadid]*sizeof(float)));
	  ifloat++;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  memcpy(s->shterms[threadid][ishort],&((*((short ***)d->dataptr))[threadid][0]),(p->NJD[threadid]*sizeof(short)));
	  ishort++;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  memcpy(s->cterms[threadid][ichar],&((*((char ***)d->dataptr))[threadid][0]),(p->NJD[threadid]*sizeof(char)));
	  ichar++;
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      } else if(Nc > 0) {
	for(u=0; u < Nc; u++) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    memcpy(s->dblterms[threadid][idbl],&((*((double ****)d->dataptr))[threadid][u][0]),(p->NJD[threadid]*sizeof(double)));
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    memcpy(s->dblterms[threadid][idbl],&((*((double ****)d->dataptr))[threadid][u][0]),(p->NJD[threadid]*sizeof(double)));
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    for(k=0;k<p->NJD[threadid];k++) {
	      memcpy(&(s->sterms[threadid][k][istring]),((*((char *****)d->dataptr))[threadid][u][k]),d->maxstringlength);
	    }
	    istring += d->maxstringlength;
	    break;
	  case VARTOOLS_TYPE_INT:
	    memcpy(s->iterms[threadid][iint],&((*((int ****)d->dataptr))[threadid][u][0]),(p->NJD[threadid]*sizeof(int)));
	    iint++;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    memcpy(s->lterms[threadid][ilong],&((*((long ****)d->dataptr))[threadid][u][0]),(p->NJD[threadid]*sizeof(long)));
	    ilong++;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    memcpy(s->fterms[threadid][ifloat],&((*((float ****)d->dataptr))[threadid][u][0]),(p->NJD[threadid]*sizeof(float)));
	    ifloat++;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    memcpy(s->shterms[threadid][ishort],&((*((short ****)d->dataptr))[threadid][u][0]),(p->NJD[threadid]*sizeof(short)));
	    ishort++;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    memcpy(s->iterms[threadid][ichar],&((*((char ****)d->dataptr))[threadid][u][0]),(p->NJD[threadid]*sizeof(char)));
	    ichar++;
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      } else {
	for(u=0; u < (-Nc); u++) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    for(k=0;k<p->NJD[threadid];k++) {
	      memcpy(&(s->dblterms[threadid][idbl][k]),&((*((double ****)d->dataptr))[threadid][k][u]),(sizeof(double)));
	    }
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    for(k=0;k<p->NJD[threadid];k++) {
	      memcpy(&(s->dblterms[threadid][idbl][k]),&((*((double ****)d->dataptr))[threadid][k][u]),(sizeof(double)));
	    }
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    for(k=0;k<p->NJD[threadid];k++) {
	      memcpy(&(s->sterms[threadid][k][istring]),((*((char *****)d->dataptr))[threadid][k][u]),d->maxstringlength);
	    }
	    istring += d->maxstringlength;
	    break;
	  case VARTOOLS_TYPE_INT:
	    for(k=0;k<p->NJD[threadid];k++) {
	      memcpy(&(s->iterms[threadid][iint][k]),&((*((int ****)d->dataptr))[threadid][k][u]),(sizeof(int)));
	    }
	    iint++;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    for(k=0;k<p->NJD[threadid];k++) {
	      memcpy(&(s->lterms[threadid][ilong][k]),&((*((long ****)d->dataptr))[threadid][k][u]),(sizeof(long)));
	    }
	    ilong++;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    for(k=0;k<p->NJD[threadid];k++) {
	      memcpy(&(s->fterms[threadid][ifloat][k]),&((*((float ****)d->dataptr))[threadid][k][u]),(sizeof(float)));
	    }
	    ifloat++;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    for(k=0;k<p->NJD[threadid];k++) {
	      memcpy(&(s->shterms[threadid][ishort][k]),&((*((short ****)d->dataptr))[threadid][k][u]),(sizeof(short)));
	    }
	    ishort++;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    for(k=0;k<p->NJD[threadid];k++) {
	      memcpy(&(s->cterms[threadid][ichar][k]),&((*((char ****)d->dataptr))[threadid][k][u]),(sizeof(char)));
	    }
	    ichar++;
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      }
    }
  if(p->readimagestring)
    {
      memcpy(s->stringid_idx[threadid],p->stringid_idx[threadid],(p->NJD[threadid]*sizeof(int)));
    }
}

void dorestorelc(ProgramData *p, _Savelc *s, _Restorelc *r, int sthreadid, int rthreadid, int lcid)
{
  int idbl, istring, iint, ishort, ilong, ifloat, ichar, i, k, u, ll;
  _DataFromLightCurve *d;
  int Nc;

  if(r == NULL ? 1 : !r->ispartialrestore) {
    if(s->NJD[sthreadid] > p->sizesinglelc[rthreadid])
      MemAllocDataFromLightCurveMidProcess(p, rthreadid, s->NJD[sthreadid]);
    
    for(idbl=0,istring=0,ishort=0,ilong=0,ifloat=0,ichar=0,iint=0,i=0;
	i<p->NDataFromLightCurve;i++)
      {
	d = &(p->DataFromLightCurve[i]);
	Nc = d->Ncolumns;
	if(Nc == 0) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    memcpy(&((*((double ***)d->dataptr))[rthreadid][0]),s->dblterms[sthreadid][idbl],(s->NJD[sthreadid]*sizeof(double)));
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    memcpy(&((*((double ***)d->dataptr))[rthreadid][0]),s->dblterms[sthreadid][idbl],(s->NJD[sthreadid]*sizeof(double)));
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    for(k=0;k<s->NJD[sthreadid];k++) {
	      memcpy(((*((char ****)d->dataptr))[rthreadid][k]),&(s->sterms[sthreadid][k][istring]),d->maxstringlength);
	    }
	    istring += d->maxstringlength;
	    break;
	  case VARTOOLS_TYPE_INT:
	    memcpy(&((*((int ***)d->dataptr))[rthreadid][0]),s->iterms[sthreadid][iint],(s->NJD[sthreadid]*sizeof(int)));
	    iint++;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    memcpy(&((*((long ***)d->dataptr))[rthreadid][0]),s->lterms[sthreadid][ilong],(s->NJD[sthreadid]*sizeof(long)));
	    ilong++;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    memcpy(&((*((float ***)d->dataptr))[rthreadid][0]),s->fterms[sthreadid][ifloat],(s->NJD[sthreadid]*sizeof(float)));
	    ifloat++;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    memcpy(&((*((short ***)d->dataptr))[rthreadid][0]),s->shterms[sthreadid][ishort],(s->NJD[sthreadid]*sizeof(short)));
	    ishort++;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    memcpy(&((*((char ***)d->dataptr))[rthreadid][0]),s->iterms[sthreadid][ichar],(s->NJD[sthreadid]*sizeof(char)));
	    ichar++;
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	} else if(Nc > 0) {
	  for(u=0; u < Nc; u++) {
	    switch(d->datatype) {
	    case VARTOOLS_TYPE_DOUBLE:
	      memcpy(&((*((double ****)d->dataptr))[rthreadid][u][0]),s->dblterms[sthreadid][idbl],(s->NJD[sthreadid]*sizeof(double)));
	      idbl++;
	      break;
	    case VARTOOLS_TYPE_CONVERTJD:
	      memcpy(&((*((double ****)d->dataptr))[rthreadid][u][0]),s->dblterms[sthreadid][idbl],(s->NJD[sthreadid]*sizeof(double)));
	      idbl++;
	      break;
	    case VARTOOLS_TYPE_STRING:
	      for(k=0;k<s->NJD[sthreadid];k++) {
		memcpy(((*((char *****)d->dataptr))[rthreadid][u][k]),&(s->sterms[sthreadid][k][istring]),d->maxstringlength);
	      }
	      istring += d->maxstringlength;
	      break;
	    case VARTOOLS_TYPE_INT:
	      memcpy(&((*((int ****)d->dataptr))[rthreadid][u][0]),s->iterms[sthreadid][iint],(s->NJD[sthreadid]*sizeof(int)));
	      iint++;
	      break;
	    case VARTOOLS_TYPE_LONG:
	      memcpy(&((*((long ****)d->dataptr))[rthreadid][u][0]),s->lterms[sthreadid][ilong],(s->NJD[sthreadid]*sizeof(long)));
	      ilong++;
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      memcpy(&((*((float ****)d->dataptr))[rthreadid][u][0]),s->fterms[sthreadid][ifloat],(s->NJD[sthreadid]*sizeof(float)));
	      ifloat++;
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      memcpy(&((*((short ****)d->dataptr))[rthreadid][u][0]),s->shterms[sthreadid][ishort],(s->NJD[sthreadid]*sizeof(short)));
	      ishort++;
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      memcpy(&((*((char ****)d->dataptr))[rthreadid][u][0]),s->iterms[sthreadid][ichar],(s->NJD[sthreadid]*sizeof(char)));
	      ichar++;
	      break;
	    default:
	      error(ERR_BADTYPE);
	    }
	  }
	} else {
	  for(u=0; u < (-Nc); u++) {
	    switch(d->datatype) {
	    case VARTOOLS_TYPE_DOUBLE:
	      for(k=0;k<s->NJD[sthreadid];k++) {
		memcpy(&((*((double ****)d->dataptr))[rthreadid][k][u]),&(s->dblterms[sthreadid][idbl][k]),(sizeof(double)));
	      }
	      idbl++;
	      break;
	    case VARTOOLS_TYPE_CONVERTJD:
	      for(k=0;k<s->NJD[sthreadid];k++) {
		memcpy(&((*((double ****)d->dataptr))[rthreadid][k][u]),&(s->dblterms[sthreadid][idbl][k]),(sizeof(double)));
	      }
	      idbl++;
	      break;
	    case VARTOOLS_TYPE_STRING:
	      for(k=0;k<s->NJD[sthreadid];k++) {
		memcpy(((*((char *****)d->dataptr))[rthreadid][k][u]),&(s->sterms[sthreadid][k][istring]),d->maxstringlength);
	      }
	      istring += d->maxstringlength;
	      break;
	    case VARTOOLS_TYPE_INT:
	      for(k=0;k<s->NJD[sthreadid];k++) {
		memcpy(&((*((int ****)d->dataptr))[rthreadid][k][u]),&(s->iterms[sthreadid][iint][k]),(sizeof(int)));
	      }
	      iint++;
	      break;
	    case VARTOOLS_TYPE_LONG:
	      for(k=0;k<s->NJD[sthreadid];k++) {
		memcpy(&((*((long ****)d->dataptr))[rthreadid][k][u]),&(s->lterms[sthreadid][ilong][k]),(sizeof(long)));
	      }
	      ilong++;
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      for(k=0;k<s->NJD[sthreadid];k++) {
		memcpy(&((*((float ****)d->dataptr))[rthreadid][k][u]),&(s->fterms[sthreadid][ifloat][k]),(sizeof(float)));
	      }
	      ifloat++;
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      for(k=0;k<s->NJD[sthreadid];k++) {
		memcpy(&((*((short ****)d->dataptr))[rthreadid][k][u]),&(s->shterms[sthreadid][ishort][k]),(sizeof(short)));
	      }
	      ishort++;
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      for(k=0;k<s->NJD[sthreadid];k++) {
		memcpy(&((*((char ****)d->dataptr))[rthreadid][k][u]),&(s->cterms[sthreadid][ichar][k]),(sizeof(char)));
	      }
	      ichar++;
	      break;
	    default:
	      error(ERR_BADTYPE);
	    }
	  }
	}
      }
    if(p->readimagestring)
      {
	memcpy(p->stringid_idx[rthreadid],s->stringid_idx[sthreadid],(s->NJD[sthreadid]*sizeof(int)));
      }

    p->NJD[rthreadid] = s->NJD[sthreadid];
  } else {
    
    /* Handle the case where the list of variables to restore is explicitly given */
    if(s->NJD[sthreadid] != p->NJD[rthreadid]) {
      /* The light curve vectors are different sizes! */
      fprintf(stderr,"Warning: the saved light curve vectors are a different size from the current vectors for light curve number %d. A partial restoration of light curve vectors is not supported when the vectors are different lengths. The light curve is being skipped.\n", lcid);
      return;
    }

    for(idbl=0,istring=0,ishort=0,ilong=0,ifloat=0,ichar=0,iint=0,i=0;
	i<p->NDataFromLightCurve;i++)
      {
	d = &(p->DataFromLightCurve[i]);
	Nc = d->Ncolumns;
	if(Nc == 0) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    for(ll=0; ll < r->Nrestorevars; ll++) {
	      if((*((double ***)d->dataptr))[rthreadid] == (*((double ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		memcpy(&((*((double ***)d->dataptr))[rthreadid][0]),s->dblterms[sthreadid][idbl],(s->NJD[sthreadid]*sizeof(double)));
	      }
	    }
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    for(ll=0; ll < r->Nrestorevars; ll++) {
	      if((*((double ***)d->dataptr))[rthreadid] == (*((double ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		memcpy(&((*((double ***)d->dataptr))[rthreadid][0]),s->dblterms[sthreadid][idbl],(s->NJD[sthreadid]*sizeof(double)));
	      }
	    }
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    for(ll=0; ll < r->Nrestorevars; ll++) {
	      if((*((char ****)d->dataptr))[rthreadid] == (*((char ****)r->restorevars[ll]->dataptr))[rthreadid]) {
		for(k=0;k<s->NJD[sthreadid];k++) {
		  memcpy(((*((char ****)d->dataptr))[rthreadid][k]),&(s->sterms[sthreadid][k][istring]),d->maxstringlength);
		}
	      }
	    }
	    istring += d->maxstringlength;
	    break;
	  case VARTOOLS_TYPE_INT:
	    for(ll=0; ll < r->Nrestorevars; ll++) {
	      if((*((int ***)d->dataptr))[rthreadid] == (*((int ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		memcpy(&((*((int ***)d->dataptr))[rthreadid][0]),s->iterms[sthreadid][iint],(s->NJD[sthreadid]*sizeof(int)));
	      }
	    }
	    iint++;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    for(ll=0; ll < r->Nrestorevars; ll++) {
	      if((*((long ***)d->dataptr))[rthreadid] == (*((long ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		memcpy(&((*((long ***)d->dataptr))[rthreadid][0]),s->lterms[sthreadid][ilong],(s->NJD[sthreadid]*sizeof(long)));
	      }
	    }
	    ilong++;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    for(ll=0; ll < r->Nrestorevars; ll++) {
	      if((*((float ***)d->dataptr))[rthreadid] == (*((float ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		memcpy(&((*((float ***)d->dataptr))[rthreadid][0]),s->fterms[sthreadid][ifloat],(s->NJD[sthreadid]*sizeof(float)));
	      }
	    }
	    ifloat++;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    for(ll=0; ll < r->Nrestorevars; ll++) {
	      if((*((short ***)d->dataptr))[rthreadid] == (*((short ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		memcpy(&((*((short ***)d->dataptr))[rthreadid][0]),s->shterms[sthreadid][ishort],(s->NJD[sthreadid]*sizeof(short)));
	      }
	    }
	    ishort++;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    for(ll=0; ll < r->Nrestorevars; ll++) {
	      if((*((char ***)d->dataptr))[rthreadid] == (*((char ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		memcpy(&((*((char ***)d->dataptr))[rthreadid][0]),s->iterms[sthreadid][ichar],(s->NJD[sthreadid]*sizeof(char)));
	      }
	    }
	    ichar++;
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	} else if(Nc > 0) {
	  for(u=0; u < Nc; u++) {
	    switch(d->datatype) {
	    case VARTOOLS_TYPE_DOUBLE:
	      for(ll=0; ll < r->Nrestorevars; ll++) {
		if((*((double ****)d->dataptr))[rthreadid][u] == (*((double ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		  memcpy(&((*((double ****)d->dataptr))[rthreadid][u][0]),s->dblterms[sthreadid][idbl],(s->NJD[sthreadid]*sizeof(double)));
		}
	      }
	      idbl++;
	      break;
	    case VARTOOLS_TYPE_CONVERTJD:
	      for(ll=0; ll < r->Nrestorevars; ll++) {
		if((*((double ****)d->dataptr))[rthreadid][u] == (*((double ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		  memcpy(&((*((double ****)d->dataptr))[rthreadid][u][0]),s->dblterms[sthreadid][idbl],(s->NJD[sthreadid]*sizeof(double)));
		}
	      }
	      idbl++;
	      break;
	    case VARTOOLS_TYPE_STRING:
	      for(ll=0; ll < r->Nrestorevars; ll++) {
		if((*((char *****)d->dataptr))[rthreadid][u] == (*((char ****)r->restorevars[ll]->dataptr))[rthreadid]) {
		  for(k=0;k<s->NJD[sthreadid];k++) {
		    memcpy(((*((char *****)d->dataptr))[rthreadid][u][k]),&(s->sterms[sthreadid][k][istring]),d->maxstringlength);
		  }
		}
	      }
	      istring += d->maxstringlength;
	      break;
	    case VARTOOLS_TYPE_INT:
	      for(ll=0; ll < r->Nrestorevars; ll++) {
		if((*((int ****)d->dataptr))[rthreadid][u] == (*((int ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		  memcpy(&((*((int ****)d->dataptr))[rthreadid][u][0]),s->iterms[sthreadid][iint],(s->NJD[sthreadid]*sizeof(int)));
		}
	      }
	      iint++;
	      break;
	    case VARTOOLS_TYPE_LONG:
	      for(ll=0; ll < r->Nrestorevars; ll++) {
		if((*((long ****)d->dataptr))[rthreadid][u] == (*((long ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		  memcpy(&((*((long ****)d->dataptr))[rthreadid][u][0]),s->lterms[sthreadid][ilong],(s->NJD[sthreadid]*sizeof(long)));
		}
	      }
	      ilong++;
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      for(ll=0; ll < r->Nrestorevars; ll++) {
		if((*((float ****)d->dataptr))[rthreadid][u] == (*((float ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		  memcpy(&((*((float ****)d->dataptr))[rthreadid][u][0]),s->fterms[sthreadid][ifloat],(s->NJD[sthreadid]*sizeof(float)));
		}
	      }
	      ifloat++;
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      for(ll=0; ll < r->Nrestorevars; ll++) {
		if((*((short ****)d->dataptr))[rthreadid][u] == (*((short ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		  memcpy(&((*((short ****)d->dataptr))[rthreadid][u][0]),s->shterms[sthreadid][ishort],(s->NJD[sthreadid]*sizeof(short)));
		}
	      }
	      ishort++;
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      for(ll=0; ll < r->Nrestorevars; ll++) {
		if((*((char ****)d->dataptr))[rthreadid][u] == (*((char ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		  memcpy(&((*((char ****)d->dataptr))[rthreadid][u][0]),s->iterms[sthreadid][ichar],(s->NJD[sthreadid]*sizeof(char)));
		}
	      }
	      ichar++;
	      break;
	    default:
	      error(ERR_BADTYPE);
	    }
	  }
	} else {
	  for(u=0; u < (-Nc); u++) {
	    switch(d->datatype) {
	    case VARTOOLS_TYPE_DOUBLE:
	      for(ll=0; ll < r->Nrestorevars; ll++) {
		if(&(*((double ****)d->dataptr))[rthreadid][0][u] == (*((double ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		  for(k=0;k<s->NJD[sthreadid];k++) {
		    memcpy(&((*((double ****)d->dataptr))[rthreadid][k][u]),&(s->dblterms[sthreadid][idbl][k]),(sizeof(double)));
		  }
		}
	      }
	      idbl++;
	      break;
	    case VARTOOLS_TYPE_CONVERTJD:
	      for(ll=0; ll < r->Nrestorevars; ll++) {
		if(&(*((double ****)d->dataptr))[rthreadid][0][u] == (*((double ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		  for(k=0;k<s->NJD[sthreadid];k++) {
		    memcpy(&((*((double ****)d->dataptr))[rthreadid][k][u]),&(s->dblterms[sthreadid][idbl][k]),(sizeof(double)));
		  }
		}
	      }
	      idbl++;
	      break;
	    case VARTOOLS_TYPE_STRING:
	      for(ll=0; ll < r->Nrestorevars; ll++) {
		if(&(*((char *****)d->dataptr))[rthreadid][0][u] == (*((char ****)r->restorevars[ll]->dataptr))[rthreadid]) {
		  for(k=0;k<s->NJD[sthreadid];k++) {
		    memcpy(((*((char *****)d->dataptr))[rthreadid][k][u]),&(s->sterms[sthreadid][k][istring]),d->maxstringlength);
		  }
		}
	      }
	      istring += d->maxstringlength;
	      break;
	    case VARTOOLS_TYPE_INT:
	      for(ll=0; ll < r->Nrestorevars; ll++) {
		if(&(*((int ****)d->dataptr))[rthreadid][0][u] == (*((int ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		  for(k=0;k<s->NJD[sthreadid];k++) {
		    memcpy(&((*((int ****)d->dataptr))[rthreadid][k][u]),&(s->iterms[sthreadid][iint][k]),(sizeof(int)));
		  }
		}
	      }
	      iint++;
	      break;
	    case VARTOOLS_TYPE_LONG:
	      for(ll=0; ll < r->Nrestorevars; ll++) {
		if(&(*((long ****)d->dataptr))[rthreadid][0][u] == (*((long ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		  for(k=0;k<s->NJD[sthreadid];k++) {
		    memcpy(&((*((long ****)d->dataptr))[rthreadid][k][u]),&(s->lterms[sthreadid][ilong][k]),(sizeof(long)));
		  }
		}
	      }
	      ilong++;
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      for(ll=0; ll < r->Nrestorevars; ll++) {
		if(&(*((float ****)d->dataptr))[rthreadid][0][u] == (*((float ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		  for(k=0;k<s->NJD[sthreadid];k++) {
		    memcpy(&((*((float ****)d->dataptr))[rthreadid][k][u]),&(s->fterms[sthreadid][ifloat][k]),(sizeof(float)));
		  }
		}
	      }
	      ifloat++;
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      for(ll=0; ll < r->Nrestorevars; ll++) {
		if(&(*((short ****)d->dataptr))[rthreadid][0][u] == (*((short ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		  for(k=0;k<s->NJD[sthreadid];k++) {
		    memcpy(&((*((short ****)d->dataptr))[rthreadid][k][u]),&(s->shterms[sthreadid][ishort][k]),(sizeof(short)));
		  }
		}
	      }
	      ishort++;
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      for(ll=0; ll < r->Nrestorevars; ll++) {
		if(&(*((char ****)d->dataptr))[rthreadid][0][u] == (*((char ***)r->restorevars[ll]->dataptr))[rthreadid]) {
		  for(k=0;k<s->NJD[sthreadid];k++) {
		    memcpy(&((*((char ****)d->dataptr))[rthreadid][k][u]),&(s->cterms[sthreadid][ichar][k]),(sizeof(char)));
		  }
		}
	      }
	      ichar++;
	      break;
	    default:
	      error(ERR_BADTYPE);
	    }
	  }
	}
      }
    if(p->readimagestring)
      {
	for(ll=0; ll < r->Nrestorevars; ll++) {
	  if(p->stringid[rthreadid] == (*((char ****) r->restorevars[ll]->dataptr))[rthreadid]) {
	    memcpy(p->stringid_idx[rthreadid],s->stringid_idx[sthreadid],(s->NJD[sthreadid]*sizeof(int)));
	  }
	}
      }

  }


}


