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

/* This file contains functions to perform the -restricttimes command for the
   vartools program by J. Hartman */

void RestrictTimes_InitStoreTimes(_RestrictTimes *c, ProgramData *p, int threadid){
  int i, j, idbl, istring, iint, ilong, ishort, ichar, ifloat, k, u, jdidsave;
  _DataFromLightCurve *d;
  int Nc;
  _Savelc *s = c->s;

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
  s->NJD[threadid] = 0;
}

void RestrictTimes_StoreTime(_RestrictTimes *c, ProgramData *p, int threadid, int jdidfull){
  int i, j, idbl, istring, iint, ilong, ishort, ichar, ifloat, k, u, jdidsave;
  _DataFromLightCurve *d;
  int Nc;
  _Savelc *s = c->s;

  jdidsave = s->NJD[threadid];
  s->NJD[threadid] = s->NJD[threadid] + 1;

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

  /* copy the light curve point */
  for(idbl=0,istring=0,ishort=0,ilong=0,ifloat=0,ichar=0,iint=0,i=0;
      i<p->NDataFromLightCurve;i++)
    {
      d = &(p->DataFromLightCurve[i]);
      Nc = d->Ncolumns;
      if(Nc == 0) {
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  memcpy(&(s->dblterms[threadid][idbl][jdidsave]),&((*((double ***)d->dataptr))[threadid][jdidfull]),sizeof(double));
	  idbl++;
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  memcpy(&(s->dblterms[threadid][idbl][jdidsave]),&((*((double ***)d->dataptr))[threadid][jdidfull]),sizeof(double));
	  idbl++;
	  break;
	case VARTOOLS_TYPE_STRING:
	  memcpy(&(s->sterms[threadid][jdidsave][istring]),((*((char ****)d->dataptr))[threadid][jdidfull]),d->maxstringlength);
	  istring += d->maxstringlength;
	  break;
	case VARTOOLS_TYPE_INT:
	  memcpy(&(s->iterms[threadid][iint][jdidsave]),&((*((int ***)d->dataptr))[threadid][jdidfull]),sizeof(int));
	  iint++;
	  break;
	case VARTOOLS_TYPE_LONG:
	  memcpy(&(s->lterms[threadid][ilong][jdidsave]),&((*((long ***)d->dataptr))[threadid][jdidfull]),sizeof(long));
	  ilong++;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  memcpy(&(s->fterms[threadid][ifloat][jdidsave]),&((*((float ***)d->dataptr))[threadid][jdidfull]),sizeof(float));
	  ifloat++;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  memcpy(&(s->shterms[threadid][ishort][jdidsave]),&((*((short ***)d->dataptr))[threadid][jdidfull]),sizeof(short));
	  ishort++;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  memcpy(&(s->cterms[threadid][ichar][jdidsave]),&((*((char ***)d->dataptr))[threadid][jdidfull]),sizeof(char));
	  ichar++;
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      } else if(Nc > 0) {
	for(u=0; u < Nc; u++) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    memcpy(&(s->dblterms[threadid][idbl][jdidsave]),&((*((double ****)d->dataptr))[threadid][u][jdidfull]),sizeof(double));
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    memcpy(&(s->dblterms[threadid][idbl][jdidsave]),&((*((double ****)d->dataptr))[threadid][u][jdidfull]),sizeof(double));
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    memcpy(&(s->sterms[threadid][jdidsave][istring]),((*((char *****)d->dataptr))[threadid][u][jdidfull]),d->maxstringlength);
	    istring += d->maxstringlength;
	    break;
	  case VARTOOLS_TYPE_INT:
	    memcpy(&(s->iterms[threadid][iint][jdidsave]),&((*((int ****)d->dataptr))[threadid][u][jdidfull]),sizeof(int));
	    iint++;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    memcpy(&(s->lterms[threadid][ilong][jdidsave]),&((*((long ****)d->dataptr))[threadid][u][jdidfull]),sizeof(long));
	    ilong++;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    memcpy(&(s->fterms[threadid][ifloat][jdidsave]),&((*((float ****)d->dataptr))[threadid][u][jdidfull]),sizeof(float));
	    ifloat++;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    memcpy(&(s->shterms[threadid][ishort][jdidsave]),&((*((short ****)d->dataptr))[threadid][u][jdidfull]),sizeof(short));
	    ishort++;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    memcpy(&(s->iterms[threadid][ichar][jdidsave]),&((*((char ****)d->dataptr))[threadid][u][jdidfull]),sizeof(char));
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
	    memcpy(&(s->dblterms[threadid][idbl][jdidsave]),&((*((double ****)d->dataptr))[threadid][jdidfull][u]),(sizeof(double)));
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    memcpy(&(s->dblterms[threadid][idbl][jdidsave]),&((*((double ****)d->dataptr))[threadid][jdidfull][u]),(sizeof(double)));
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    memcpy(&(s->sterms[threadid][jdidsave][istring]),((*((char *****)d->dataptr))[threadid][jdidfull][u]),d->maxstringlength);
	    istring += d->maxstringlength;
	    break;
	  case VARTOOLS_TYPE_INT:
	    memcpy(&(s->iterms[threadid][iint][jdidsave]),&((*((int ****)d->dataptr))[threadid][jdidfull][u]),(sizeof(int)));
	    iint++;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    memcpy(&(s->lterms[threadid][ilong][jdidsave]),&((*((long ****)d->dataptr))[threadid][jdidfull][u]),(sizeof(long)));
	    ilong++;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    memcpy(&(s->fterms[threadid][ifloat][jdidsave]),&((*((float ****)d->dataptr))[threadid][jdidfull][u]),(sizeof(float)));
	    ifloat++;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    memcpy(&(s->shterms[threadid][ishort][jdidsave]),&((*((short ****)d->dataptr))[threadid][jdidfull][u]),(sizeof(short)));
	    ishort++;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    memcpy(&(s->cterms[threadid][ichar][jdidsave]),&((*((char ****)d->dataptr))[threadid][jdidfull][u]),(sizeof(char)));
	    ichar++;
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      }
    }
}

/* This is the function that is run when the -restoretimes command is called to
   add excluded times back into the light curve */
void RestoreTimes(ProgramData *p, _RestoreTimes *RestoreTimes, int sthreadid, int rthreadid){
  _Savelc *s;
  int idbl, istring, iint, ishort, ilong, ifloat, ichar, i, k, u;
  int copyindx;
  _DataFromLightCurve *d;
  int Nc;

  s = RestoreTimes->RestrictTimes->s;

  /* Return if there are no points to add back to the light curve */
  if(!s->NJD[sthreadid]) return;

  copyindx = p->NJD[rthreadid];

  if(s->NJD[sthreadid] + p->NJD[rthreadid] > p->sizesinglelc[rthreadid])
    MemAllocDataFromLightCurve(p, rthreadid, (s->NJD[sthreadid] + p->NJD[rthreadid]));

  for(idbl=0,istring=0,ishort=0,ilong=0,ifloat=0,ichar=0,iint=0,i=0;
      i<p->NDataFromLightCurve;i++)
    {
      d = &(p->DataFromLightCurve[i]);
      Nc = d->Ncolumns;
      if(Nc == 0) {
	switch(d->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	  memcpy(&((*((double ***)d->dataptr))[rthreadid][copyindx]),s->dblterms[sthreadid][idbl],(s->NJD[sthreadid]*sizeof(double)));
	  idbl++;
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  memcpy(&((*((double ***)d->dataptr))[rthreadid][copyindx]),s->dblterms[sthreadid][idbl],(s->NJD[sthreadid]*sizeof(double)));
	  idbl++;
	  break;
	case VARTOOLS_TYPE_STRING:
	  for(k=0;k<s->NJD[sthreadid];k++) {
	    memcpy(((*((char ****)d->dataptr))[rthreadid][copyindx+k]),&(s->sterms[sthreadid][k][istring]),d->maxstringlength);
	  }
	  istring += d->maxstringlength;
	  break;
	case VARTOOLS_TYPE_INT:
	  memcpy(&((*((int ***)d->dataptr))[rthreadid][copyindx]),s->iterms[sthreadid][iint],(s->NJD[sthreadid]*sizeof(int)));
	  iint++;
	  break;
	case VARTOOLS_TYPE_LONG:
	  memcpy(&((*((long ***)d->dataptr))[rthreadid][copyindx]),s->lterms[sthreadid][ilong],(s->NJD[sthreadid]*sizeof(long)));
	  ilong++;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  memcpy(&((*((float ***)d->dataptr))[rthreadid][copyindx]),s->fterms[sthreadid][ifloat],(s->NJD[sthreadid]*sizeof(float)));
	  ifloat++;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  memcpy(&((*((short ***)d->dataptr))[rthreadid][copyindx]),s->shterms[sthreadid][ishort],(s->NJD[sthreadid]*sizeof(short)));
	  ishort++;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  memcpy(&((*((char ***)d->dataptr))[rthreadid][copyindx]),s->iterms[sthreadid][ichar],(s->NJD[sthreadid]*sizeof(char)));
	  ichar++;
	  break;
	default:
	  error(ERR_BADTYPE);
	}
      } else if(Nc > 0) {
	for(u=0; u < Nc; u++) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    memcpy(&((*((double ****)d->dataptr))[rthreadid][u][copyindx]),s->dblterms[sthreadid][idbl],(s->NJD[sthreadid]*sizeof(double)));
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    memcpy(&((*((double ****)d->dataptr))[rthreadid][u][copyindx]),s->dblterms[sthreadid][idbl],(s->NJD[sthreadid]*sizeof(double)));
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    for(k=0;k<s->NJD[sthreadid];k++) {
	      memcpy(((*((char *****)d->dataptr))[rthreadid][u][copyindx+k]),&(s->sterms[sthreadid][k][istring]),d->maxstringlength);
	    }
	    istring += d->maxstringlength;
	    break;
	  case VARTOOLS_TYPE_INT:
	    memcpy(&((*((int ****)d->dataptr))[rthreadid][u][copyindx]),s->iterms[sthreadid][iint],(s->NJD[sthreadid]*sizeof(int)));
	    iint++;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    memcpy(&((*((long ****)d->dataptr))[rthreadid][u][copyindx]),s->lterms[sthreadid][ilong],(s->NJD[sthreadid]*sizeof(long)));
	    ilong++;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    memcpy(&((*((float ****)d->dataptr))[rthreadid][u][copyindx]),s->fterms[sthreadid][ifloat],(s->NJD[sthreadid]*sizeof(float)));
	    ifloat++;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    memcpy(&((*((short ****)d->dataptr))[rthreadid][u][copyindx]),s->shterms[sthreadid][ishort],(s->NJD[sthreadid]*sizeof(short)));
	    ishort++;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    memcpy(&((*((char ****)d->dataptr))[rthreadid][u][copyindx]),s->iterms[sthreadid][ichar],(s->NJD[sthreadid]*sizeof(char)));
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
	      memcpy(&((*((double ****)d->dataptr))[rthreadid][copyindx+k][u]),&(s->dblterms[sthreadid][idbl][k]),(sizeof(double)));
	    }
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    for(k=0;k<s->NJD[sthreadid];k++) {
	      memcpy(&((*((double ****)d->dataptr))[rthreadid][copyindx+k][u]),&(s->dblterms[sthreadid][idbl][k]),(sizeof(double)));
	    }
	    idbl++;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    for(k=0;k<s->NJD[sthreadid];k++) {
	      memcpy(((*((char *****)d->dataptr))[rthreadid][copyindx+k][u]),&(s->sterms[sthreadid][k][istring]),d->maxstringlength);
	    }
	    istring += d->maxstringlength;
	    break;
	  case VARTOOLS_TYPE_INT:
	    for(k=0;k<s->NJD[sthreadid];k++) {
	      memcpy(&((*((int ****)d->dataptr))[rthreadid][copyindx+k][u]),&(s->iterms[sthreadid][iint][k]),(sizeof(int)));
	    }
	    iint++;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    for(k=0;k<s->NJD[sthreadid];k++) {
	      memcpy(&((*((long ****)d->dataptr))[rthreadid][copyindx+k][u]),&(s->lterms[sthreadid][ilong][k]),(sizeof(long)));
	    }
	    ilong++;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    for(k=0;k<s->NJD[sthreadid];k++) {
	      memcpy(&((*((float ****)d->dataptr))[rthreadid][copyindx+k][u]),&(s->fterms[sthreadid][ifloat][k]),(sizeof(float)));
	    }
	    ifloat++;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    for(k=0;k<s->NJD[sthreadid];k++) {
	      memcpy(&((*((short ****)d->dataptr))[rthreadid][copyindx+k][u]),&(s->shterms[sthreadid][ishort][k]),(sizeof(short)));
	    }
	    ishort++;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    for(k=0;k<s->NJD[sthreadid];k++) {
	      memcpy(&((*((char ****)d->dataptr))[rthreadid][copyindx+k][u]),&(s->cterms[sthreadid][ichar][k]),(sizeof(char)));
	    }
	    ichar++;
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
      }
    }


  p->NJD[rthreadid] = p->NJD[rthreadid] + s->NJD[sthreadid];
  if(p->readimagestring)
    {
      for(i=0;i<p->NJD[rthreadid];i++)
	p->stringid_idx[rthreadid][i] = i;
      mysortstringint(p->NJD[rthreadid],MAXIDSTRINGLENGTH, p->stringid[rthreadid], p->stringid_idx[rthreadid]);
    }
  sortlcbytime(p->NJD[rthreadid], p->t[rthreadid], rthreadid, p);
}

void RestrictTimes_ParseExpr(int *iret, int argc, char **argv, ProgramData *p, _RestrictTimes *RestrictTimes, char min_or_max){
  if(!min_or_max) {
    if((RestrictTimes->minJDexpr = (_Expression *) malloc(sizeof(_Expression))) == NULL)
      error(ERR_MEMALLOC);
    if((RestrictTimes->minJDexprstring = (char *) malloc((strlen(argv[*iret])+1)*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
    sprintf(RestrictTimes->minJDexprstring,"%s",argv[*iret]);
  }
  else {
    if((RestrictTimes->maxJDexpr = (_Expression *) malloc(sizeof(_Expression))) == NULL)
      error(ERR_MEMALLOC);
    if((RestrictTimes->maxJDexprstring = (char *) malloc((strlen(argv[*iret])+1)*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
    sprintf(RestrictTimes->maxJDexprstring,"%s",argv[*iret]);
  }
}

void RestrictTimes_readJDlist(char *filename, double **JDlist, int *Nlist) {
  FILE *infile;
  char *line;
  size_t line_size = MAXLEN;
  double inJD;
  int sizeJDlist = 1026;
  int N = 0;
  line = malloc(line_size);
  if((*JDlist = (double *) malloc(sizeJDlist * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  if(!strcmp(filename,"-")) {
    infile = stdin;
  }
  else if((infile = fopen(filename,"r")) == NULL) {
    error2(ERR_FILENOTFOUND,filename);
  }
  while(gnu_getline(&line,&line_size,infile) >= 0) {
    if(line[0] != '#') {
      if(N >= sizeJDlist) {
	sizeJDlist *= 2;
	if((*JDlist = (double *) realloc((*JDlist),sizeJDlist * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
      }
      sscanf(line,"%lf",&((*JDlist)[N]));
      N++;
    }
  }
  if(infile != stdin)
    fclose(infile);
  mysort1(N, *JDlist);
  *Nlist = N;
  free(line);
}

void RestrictTimes_readimagelist(char *filename, char ***imagelist, int **imagelist_indx, int *Nlist) {
  FILE *infile;
  char *line;
  size_t line_size = MAXLEN;
  double inJD;
  int sizeJDlist = 1026;
  int N = 0, i;
  line = malloc(line_size);
  if((*imagelist = (char **) malloc(sizeJDlist * sizeof(char *))) == NULL)
    error(ERR_MEMALLOC);
  if(!strcmp(filename,"-")) {
    infile = stdin;
  }
  else if((infile = fopen(filename,"r")) == NULL) {
    error2(ERR_FILENOTFOUND,filename);
  }
  while(gnu_getline(&line,&line_size,infile) >= 0) {
    if(line[0] != '#') {
      if(N >= sizeJDlist) {
	sizeJDlist *= 2;
	if((*imagelist = (char **) realloc((*imagelist),sizeJDlist * sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
      }
      if(((*imagelist)[N] = (char *) malloc(MAXIDSTRINGLENGTH * sizeof(char))) == NULL)
	error(ERR_MEMALLOC);
      sscanf(line,"%s ",((*imagelist)[N]));
      N++;
    }
  }
  if(infile != stdin)
    fclose(infile);
  if(N > 0) {
    if(((*imagelist_indx) = (int *) malloc(N * sizeof(int))) == NULL) 
      error(ERR_MEMALLOC);
    for(i=0; i < N; i++)
      (*imagelist_indx)[i] = i;
    mysortstringint(N, MAXIDSTRINGLENGTH, *imagelist, *imagelist_indx);
  }
  *Nlist = N;
  free(line);
}

/* Apply the JD range restriction to a light curve; we will use the
 sigclip_copyterms function to remove points in a way that is 
 compatible with the savelc, restorelc, and multi-column support */
int RestrictTimes_JDrange_apply(int N, double *t,
				int lc, ProgramData *p,
				_RestrictTimes *c,
				double JDmin, double JDmax,
				char exclude)
{
  int i, j, test;
  if(c->saveexcludedpoints)
    RestrictTimes_InitStoreTimes(c, p, lc);
  j = 0;
  for(i=0; i < N; i++) {
    test = t[i] >= JDmin && t[i] <= JDmax;
    if((test && !exclude) || (!test && exclude)) {
      if(i != j) {
	sigclip_copyterms(i,j,p,lc);
      }
      j++;
    } else if(c->saveexcludedpoints) {
      RestrictTimes_StoreTime(c, p, lc, i);
    }
  }
  p->NJD[lc] = j;
  if(p->readimagestring)
    {
      for(i=0;i<p->NJD[lc];i++)
	p->stringid_idx[lc][i] = i;
      mysortstringint(p->NJD[lc],MAXIDSTRINGLENGTH, p->stringid[lc], p->stringid_idx[lc]);
    }
}

/* Apply the JDlist restriction to a light curve; we will use the
 sigclip_copyterms function to remove points in a way that is 
 compatible with the savelc, restorelc, and multi-column support */
int RestrictTimes_JDlist_apply(int N, double *t,
			       int lc, ProgramData *p,
			       _RestrictTimes *c,
			       double *JDlist, int Nlist,
			       char exclude)
{
  int i, j, k, test;
  if(c->saveexcludedpoints)
    RestrictTimes_InitStoreTimes(c, p, lc);
  j = 0;
  k = 0;
  for(i=0; i < N; i++) {
    while((k < Nlist ? (t[i] > JDlist[k] + JDTOL) : 0))
      k++;
    test = k < Nlist ? (t[i] < JDlist[k] + JDTOL && t[i] > JDlist[k] - JDTOL) : 0;
    if((!exclude && test) || (exclude && !test)){
      sigclip_copyterms(i,j,p,lc);
      j++;
    } else if(c->saveexcludedpoints) {
      RestrictTimes_StoreTime(c, p, lc, i);
    }
  }
  p->NJD[lc] = j;
  if(p->readimagestring)
    {
      for(i=0;i<p->NJD[lc];i++)
	p->stringid_idx[lc][i] = i;
      mysortstringint(p->NJD[lc],MAXIDSTRINGLENGTH, p->stringid[lc], p->stringid_idx[lc]);
    }
}

/* Apply the imagelist restriction to a light curve; we will use the
 sigclip_copyterms function to remove points in a way that is 
 compatible with the savelc, restorelc, and multi-column support */
void RestrictTimes_imagelist_apply(int N, char **stringID, int *stringID_indx, 
				  int lc, ProgramData *p,
				   _RestrictTimes *c,
				  char **imagelist, int *imagelist_indx, 
				  int Nlist, char exclude)
{
  int i, j, k, test;
  int *good;
  if(c->saveexcludedpoints)
    RestrictTimes_InitStoreTimes(c, p, lc);
  j = 0;
  k = 0;
  if(N <= 0)
    return;
  if((good = (int *) malloc(N * sizeof(int))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0; i < N; i++) {
    while((k < Nlist ? (strcmp(stringID[stringID_indx[i]],imagelist[imagelist_indx[k]]) > 0) : 0))
      k++;
    test = k < Nlist ? (!strcmp(stringID[stringID_indx[i]],imagelist[imagelist_indx[k]])) : 0;
    if((!exclude && test) || (exclude && !test)) {
      /* Record the index of selected points */
      good[j] = stringID_indx[i];
      j++;
    } else if(c->saveexcludedpoints) {
      RestrictTimes_StoreTime(c, p, lc, i);
    }
  }
  /* Sort the index and then copy out only those elements that have
     been selected */
  mysort1int(j, good);
  p->NJD[lc] = j;
  for(i=0; i < j; i++) {
    sigclip_copyterms(good[i],i,p,lc);
  }
  if(p->readimagestring)
    {
      for(i=0;i<p->NJD[lc];i++)
	p->stringid_idx[lc][i] = i;
      mysortstringint(p->NJD[lc],MAXIDSTRINGLENGTH, p->stringid[lc], p->stringid_idx[lc]);
    }
  free(good);
}

void RestrictTimes_expr_apply(ProgramData *p, _RestrictTimes *c, int threadindex, int lcindex) {
  int i, j, Ntest, test;
  double testdbl;
  if(c->saveexcludedpoints)
    RestrictTimes_InitStoreTimes(c, p, threadindex);
  j = 0;
  Ntest = p->NJD[threadindex];
  for(i=0; i < Ntest; i++) {
    testdbl = EvaluateExpression(lcindex, threadindex, i, c->restrictexpr);
    if(((testdbl > 0) && !c->exclude) || ((testdbl <= 0) &&  c->exclude)) {
      if(i != j) {
	sigclip_copyterms(i,j,p,threadindex);
      }
      j++;
    } else if(c->saveexcludedpoints) {
      RestrictTimes_StoreTime(c, p, threadindex, i);
    }
  }
  p->NJD[threadindex] = j;
  if(p->readimagestring)
    {
      for(i=0;i<p->NJD[threadindex];i++)
	p->stringid_idx[threadindex][i] = i;
      mysortstringint(p->NJD[threadindex],MAXIDSTRINGLENGTH, p->stringid[threadindex], p->stringid_idx[threadindex]);
    }
}
