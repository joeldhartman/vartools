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
				double JDmin, double JDmax,
				char exclude)
{
  int i, j, test;
  j = 0;
  for(i=0; i < N; i++) {
    test = t[i] >= JDmin && t[i] <= JDmax;
    if((test && !exclude) || (!test && exclude)) {
      if(i != j) {
	sigclip_copyterms(i,j,p,lc);
      }
      j++;
    }
  }
  p->NJD[lc] = j;
  if(p->readimagestring)
    {
      for(i=0;i<p->NJD[lc];i++)
	p->stringid_idx[0][i] = i;
      mysortstringint(p->NJD[lc],MAXIDSTRINGLENGTH, p->stringid[lc], p->stringid_idx[lc]);
    }
}

/* Apply the JDlist restriction to a light curve; we will use the
 sigclip_copyterms function to remove points in a way that is 
 compatible with the savelc, restorelc, and multi-column support */
int RestrictTimes_JDlist_apply(int N, double *t,
			       int lc, ProgramData *p,
			       double *JDlist, int Nlist,
			       char exclude)
{
  int i, j, k, test;
  j = 0;
  k = 0;
  for(i=0; i < N; i++) {
    while((k < Nlist ? (t[i] > JDlist[k] + JDTOL) : 0))
      k++;
    test = k < Nlist ? (t[i] < JDlist[k] + JDTOL && t[i] > JDlist[k] - JDTOL) : 0;
    if((!exclude && test) || (exclude && !test)){
      sigclip_copyterms(i,j,p,lc);
      j++;
      k++;
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
				  char **imagelist, int *imagelist_indx, 
				  int Nlist, char exclude)
{
  int i, j, k, test;
  int *good;
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
      k++;
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
