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

_CopyLC * CreateCopyLCCommand(ProgramData *p, char *argv, int cnum)
{
  _CopyLC *ret;
  if((ret = (_CopyLC *) malloc(sizeof(_CopyLC))) == NULL)
    error(ERR_MEMALLOC);
  if((ret->s = (_Savelc *) malloc(sizeof(_Savelc))) == NULL)
    error(ERR_MEMALLOC);
  if((ret->SaveListData = (_SaveListData *) malloc(sizeof(_SaveListData))) == NULL)
    error(ERR_MEMALLOC);
  ret->Ncopies = atoi(argv);
  if(ret->Ncopies <= 0) {
    error2(ERR_INVALID_PARAMETERVALUE,"-copylc, Ncopies must be > 0");
  }
  ret->cnum = cnum;
  ret->copycommand_index = p->Ncopycommands;
  ret->priorcopies = p->Ncopiestotal;
  if(!p->Ncopycommands) {
    if((p->Ncopies = (int *) malloc(sizeof(int))) == NULL ||
       (p->copy_cnum = (int *) malloc(sizeof(int))) == NULL)
      error(ERR_MEMALLOC);
    p->Ncopiestotal = 1;
  }
  else {
    if((p->Ncopies = (int *) realloc(p->Ncopies, (p->Ncopycommands+1)*sizeof(int))) == NULL ||
       (p->copy_cnum = (int *) realloc(p->copy_cnum, (p->Ncopycommands+1)*sizeof(int))) == NULL)
      error(ERR_MEMALLOC);
  }
  p->Ncopies[p->Ncopycommands] = ret->Ncopies;
  p->copy_cnum[p->Ncopycommands] = cnum;
  p->Ncopiestotal = p->Ncopiestotal * (1 + ret->Ncopies);
  p->Ncopycommands += 1;
  return ret;
}

void MakeCopyNames(ProgramData *p, int *lcindx, int copyindx) {
  int lcindx_base;
  int j, k;
  if(copyindx >= p->Ncopycommands) return;
  lcindx_base = *lcindx;
  k = *lcindx;
  /* Make the next set of copies for the base index */
  MakeCopyNames(p, &k, (copyindx+1));
  for(j=0; j < p->Ncopies[copyindx]; j++) {
    k++;
    p->copycommand_index[k] = copyindx;
    p->start_cnum[k] = p->copy_cnum[copyindx]+1;
    p->is_lc_ready[k] = 0;
    p->copy_origlc_index[k] = lcindx_base;
    sprintf(p->lcnames[k],"%s_copy%d.%d",p->lcnames[lcindx_base],(copyindx+1),(j+1));
    MakeCopyNames(p, &k, (copyindx+1));
  }
  *lcindx = k;
  return;
}

void SetupLCCopies(ProgramData *p, Command *c)
{
  int i, j;
  char *strtmp;
  j = 0;
  if((p->copycommand_index = (int *) malloc(p->Nlcs * sizeof(int))) == NULL ||
     (p->is_lc_ready = (int *) malloc(p->Nlcs * sizeof(int))) == NULL ||
     (p->copy_origlc_index = (int *) malloc(p->Nlcs * sizeof(int))) == NULL ||
     (p->start_cnum = (int *) malloc(p->Nlcs * sizeof(int))) == NULL)
    error(ERR_MEMALLOC);
  while(j < p->Nlcs) {
    p->is_lc_ready[j] = 1;
    p->copycommand_index[j] = -1;
    p->copy_origlc_index[j] = j;
    p->start_cnum[j] = 0;
    MakeCopyNames(p, &j, 0);
    j++;
  }
  for(i=0; i < p->Ncopycommands; i++) {
    if((c[p->copy_cnum[i]].CopyLC->lcid_tothreadid = (int *) malloc(p->Nlcs*sizeof(int))) == NULL)
      error(ERR_MEMALLOC);
    for(j=0; j < p->Nlcs; j++) {
      c[p->copy_cnum[i]].CopyLC->lcid_tothreadid[j] = -1;
    }
  }
}

void turnoffcopies_onecommand(ProgramData *p, Command *c, int threadid, int lcid)
{
  int copyindex;
  int i, j;
  int Nskip;
  if(c->cnum == CNUM_COPYLC) {
    copyindex = c->CopyLC->copycommand_index;
    Nskip = p->Ncopiestotal;
    if(c->CopyLC->priorcopies > 0)
      Nskip = Nskip/c->CopyLC->priorcopies;
    Nskip = Nskip / (c->CopyLC->Ncopies + 1);
  }
#ifdef PARALLEL
  while(pthread_mutex_trylock(&(p->is_lc_ready_mutex)));
#endif
  i = lcid+Nskip;
  for(j=0; j < c->CopyLC->Ncopies; j++) {
    p->is_lc_ready[i] = -1;
#ifdef PARALLEL
    pthread_mutex_unlock(&(p->is_lc_ready_mutex));
#endif
    turnoffcopies(p, &(c[0-p->copy_cnum[copyindex]]), p->copy_cnum[copyindex]+1, threadid, i);
#ifdef PARALLEL
    while(pthread_mutex_trylock(&(p->is_lc_ready_mutex)));
#endif
    i += (Nskip);
  }
#ifdef PARALLEL
  pthread_mutex_unlock(&(p->is_lc_ready_mutex));
#endif
}

void turnoffcopies(ProgramData *p, Command *c, int cnum_start, int threadid, int lcid)
{
  int copyindex;
  int i, j, k;
  int Nskip;
  for(k = cnum_start; k < p->Ncommands; k++) {
    if(c[k].cnum == CNUM_COPYLC) {
      copyindex = c[k].CopyLC->copycommand_index;
      Nskip = p->Ncopiestotal;
      if(c[k].CopyLC->priorcopies > 0)
	Nskip = Nskip/c[k].CopyLC->priorcopies;
      Nskip = Nskip / (c[k].CopyLC->Ncopies + 1);
#ifdef PARALLEL
      while(pthread_mutex_trylock(&(p->is_lc_ready_mutex)));
#endif
      i = lcid+Nskip;
      for(j=0; j < c[k].CopyLC->Ncopies; j++) {
	p->is_lc_ready[i] = -1;
#ifdef PARALLEL
	pthread_mutex_unlock(&(p->is_lc_ready_mutex));
#endif
	turnoffcopies(p, c, k+1, threadid, i);
#ifdef PARALLEL
	while(pthread_mutex_trylock(&(p->is_lc_ready_mutex)));
#endif
	i += (Nskip);
      }
#ifdef PARALLEL
      pthread_mutex_unlock(&(p->is_lc_ready_mutex));
#endif
      }
  }
}
    
void docopylccommand(ProgramData *p, _CopyLC *c, int threadid, int lcid)
{
  _DataFromLightCurve *d;
  int Nc;
  int copyindex;
  int i, j;
  int Nskip;

  copyindex = c->copycommand_index;
  Nskip = p->Ncopiestotal;
  if(c->priorcopies > 0)
    Nskip = Nskip/c->priorcopies;
  Nskip = Nskip / (c->Ncopies + 1);

  dosavelc(p, c->s, threadid, lcid);
  dosavelistdata(p, c->SaveListData, threadid, lcid);
  if(p->isifcommands) {
    dosaveifstackcopy(p, c, threadid);
  }

  c->lcid_tothreadid[lcid] = threadid;

#ifdef PARALLEL
  while(pthread_mutex_trylock(&(p->is_lc_ready_mutex)));
#endif
  i = lcid+Nskip;
  for(j=0; j < c->Ncopies; j++) {
    p->is_lc_ready[i] = 1;
    i += (Nskip);
  }
#ifdef PARALLEL
  pthread_mutex_unlock(&(p->is_lc_ready_mutex));
#endif
}

void getlccopy(ProgramData *p, Command *c, int threadid, int lcid) {
  int copy_cnum;

  int copycommand_index;
  int copy_origlc_index;
  int j;
  
  copy_cnum = p->copy_cnum[p->copycommand_index[lcid]];
  dorestorelc(p, c[copy_cnum].Savelc, c[copy_cnum].CopyLC->lcid_tothreadid[p->copy_origlc_index[lcid]], threadid);
  dorestorelistdata(p, c[copy_cnum].CopyLC->SaveListData, c[copy_cnum].CopyLC->lcid_tothreadid[p->copy_origlc_index[lcid]], threadid, lcid);
  if(p->isifcommands) {
    dorestoreifstackcopy(p, c[copy_cnum].CopyLC, c[copy_cnum].CopyLC->lcid_tothreadid[p->copy_origlc_index[lcid]], threadid);
  }
  
  /* Execute the most recent change variable request */
  for(j = copy_cnum; j >= 0; j--) {
    if(c[j].cnum == CNUM_CHANGEVARIABLE) {
      DoChangeVariable(p, c[j].Changevariable, threadid);
      break;
    }
  }
  /* Otherwise set the default vectors */
  if(j < 0) {
    SetTimeMagSigPointers(p, threadid);
  }

}
