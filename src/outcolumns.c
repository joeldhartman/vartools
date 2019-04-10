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

void increaseNcolumns(ProgramData *p, int cnum)
{
  int i, j, k;
  if(!p->Ncolumns) {
    p->outcolumns = malloc(sizeof(OutColumn));
  }
  else {
    p->outcolumns = realloc(p->outcolumns,(p->Ncolumns + 1)*sizeof(OutColumn));
  }

  if(cnum < 0) {
    (p->Ncolumns)++;
    return;
  }

  if(cnum >= p->size_colcommandvec) {
    if(p->size_colcommandvec < 0) {
      p->size_colcommandvec = 2*(cnum + 1);
      j = p->size_colcommandvec;
      p->col_commandstart = malloc(j * sizeof(int));
      p->col_commandstop = malloc(j * sizeof(int));
      for(i=0; i < j; i++) {
	p->col_commandstart[i] = -1;
	p->col_commandstop[i] = -1;
      }
    } else {
      k = p->size_colcommandvec;
      p->size_colcommandvec = 2*(cnum + 1);
      j = p->size_colcommandvec;
      p->col_commandstart = realloc(p->col_commandstart, j * sizeof(int));
      p->col_commandstop = realloc(p->col_commandstop, j * sizeof(int));
      for(i=k; i < j; i++) {
	p->col_commandstart[i] = -1;
	p->col_commandstop[i] = -1;
      }
    }
  }
  if(cnum > p->max_colcommand)
    p->max_colcommand = cnum;
  if(p->col_commandstart[cnum] < 0) {
    p->col_commandstart[cnum] = p->Ncolumns;
    p->col_commandstop[cnum] = p->Ncolumns;
  }
  else {
    p->col_commandstop[cnum] = p->Ncolumns;
  }

  (p->Ncolumns)++;
}

void addcolumn(ProgramData *p, int cnum, int type, int stringsize, void *ptr, char *outputformat, int Ndereference, int usereallc, int lcdereferencecol, ...)
{
  va_list varlist;
  OutColumn *c;
  char *columnnamefmt;
  char colnamefmtcpy[MAX_COLNAME_LENGTH];
  int i;
  increaseNcolumns(p, cnum);
  c = &(p->outcolumns[p->Ncolumns-1]);
  va_start(varlist, lcdereferencecol);
  if((c->dereference = (int *) malloc(Ndereference * sizeof(int))) == NULL)
    error(ERR_MEMALLOC);
  c->type = type;
  c->stringsize = stringsize;
  c->ptr = ptr;
  sprintf(c->outputformat,"%s",outputformat);
  c->Ndereference = Ndereference;
  c->usereallc = usereallc;
  c->lcdereferencecol = lcdereferencecol;
  for(i=0;i<Ndereference;i++)
    c->dereference[i] = va_arg(varlist,int);
  columnnamefmt=va_arg(varlist,char *);
  if(p->numbercolumns)
    {
      sprintf(colnamefmtcpy,"%d_%s",p->Ncolumns,columnnamefmt);
    }
  else
    sprintf(colnamefmtcpy,"%s",columnnamefmt);
  vsprintf(c->columnname,colnamefmtcpy,varlist);
  c->length_columnname = strlen(c->columnname);
  c->cnum = cnum;
}

/* Add a new command to the list of commands that have input read from a previous output column */
void increaselinkedcols(ProgramData *p, OutColumn **c, char *s, int cmdidx)
{
  if(!p->Ncolstolink)
    {
      p->Ncolstolink = 1;
      if((p->colnamestolink = (char **) malloc(sizeof(char *))) == NULL ||
	 (p->outcolumnstolink = (OutColumn ***) malloc(sizeof(OutColumn **))) == NULL ||
	 (p->columnstolink_cmdidx = (int *) malloc(sizeof(int))) == NULL)
	error(ERR_MEMALLOC);
      if((p->colnamestolink[0] = (char *) malloc(MAX_COLNAME_LENGTH)) == NULL)
	error(ERR_MEMALLOC);
    }
  else
    {
      p->Ncolstolink++;
      if((p->colnamestolink = (char **) realloc(p->colnamestolink, p->Ncolstolink * sizeof(char *))) == NULL ||
	 (p->outcolumnstolink = (OutColumn ***) realloc(p->outcolumnstolink, p->Ncolstolink * sizeof(OutColumn **))) == NULL ||
	 (p->columnstolink_cmdidx = (int *) realloc(p->columnstolink_cmdidx, p->Ncolstolink * sizeof(int))) == NULL)
	error(ERR_MEMALLOC);
      if((p->colnamestolink[p->Ncolstolink - 1] = (char *) malloc(MAX_COLNAME_LENGTH)) == NULL)
	error(ERR_MEMALLOC);
    }
  sprintf(p->colnamestolink[p->Ncolstolink - 1],"%s",s);
  p->outcolumnstolink[p->Ncolstolink - 1] = c;
  p->columnstolink_cmdidx[p->Ncolstolink - 1] = cmdidx;
}

/* This function returns a pointer to the column whose number or name is stored in the input string s */
OutColumn *linktocolumn(ProgramData *p, char *s, int commandidx)
{
  int i, k;
  int testint, colid, l1, l2, nonumber;
  /* First check if the string is an integer */
  for(testint = 1, i =0; i < strlen(s); i++)
    {
      if(s[i] < '0' || s[i] > '9')
	{
	  testint = 0;
	  break;
	}
    }

  /* If it's an integer get the column value */
  if(testint)
    {
      colid = atoi(s);
      colid--;
      if(colid < 0 || colid >= p->Ncolumns)
	error2(ERR_NOCOLUMN, s);
    }
  else
    {
      /* Otherwise search for a column with that name */
      colid = -1;

      nonumber = 0;
      if(p->numbercolumns) {
	if(s[0] < '0' || s[0] > '9')
	  nonumber = 1;
      }

      if(!nonumber) {
	l1 = strlen(s);
	for(i=0; i<p->Ncolumns; i++)
	  {
	    l2 = p->outcolumns[i].length_columnname;
	    if(l1 == l2 && !strncmp(s,p->outcolumns[i].columnname,l1))
	      {
		colid = i;
		break;
	      }
	  }
      }
      else {
	l1 = strlen(s);
	for(i=0; i<p->Ncolumns; i++)
	  {
	    l2 = p->outcolumns[i].length_columnname;
	    k = 0;
	    while(p->outcolumns[i].columnname[k] >= '0' &&
		  p->outcolumns[i].columnname[k] <= '9')
	      k++;
	    while(p->outcolumns[i].columnname[k] == '_')
	      k++;
	    if(l1 == (l2-k) && !strncmp(s,&(p->outcolumns[i].columnname[k]),l1))
	      {
		colid = i;
		break;
	      }
	  }
      }
      if(colid < 0 || colid >= p->Ncolumns)
	error2(ERR_NOCOLUMN, s);
    }

  /* Get the command index for the linked column */
  if(colid > 0)
    {
      l2 = p->outcolumns[colid].length_columnname;
      l2--;
      while(l2 >= 0 ? p->outcolumns[colid].columnname[l2] >= '0' && p->outcolumns[colid].columnname[l2] <= '9' : 0)
	l2--;
      l2++;
      i = atoi(&(p->outcolumns[colid].columnname[l2]));
      if(i >= commandidx)
	error2(ERR_BADCOLUMNLINK, s);
    }
  return &(p->outcolumns[colid]);
}

void linkcolumns(ProgramData *p)
{
  int i;
  for(i=0; i < p->Ncolstolink; i++)
    {
      *(p->outcolumnstolink[i]) = linktocolumn(p, p->colnamestolink[i], p->columnstolink_cmdidx[i]);
    }
}

void getoutcolumnvalue(OutColumn *c, int lc, int reallc, int outtype, void *outvalue, ...)
{
  va_list varlist;
  void *ptr;
  double *dptr;
  char **cptr;
  int *iptr;
  float *fptr;
  short *shortptr;
  long *longptr;
  char *charptr;
  int i;
  int size_out_string;
  ptr = (void *) (*((char ***) c->ptr));
  for(i=0;i<c->Ndereference-1;i++)
    {
      if(i != c->lcdereferencecol)
	ptr = (void *) (*(((char **)ptr) + c->dereference[i]));
      else if(c->usereallc)
	{
	  ptr = (void *) (*(((char **)ptr) + reallc));
	}
      else
	{
	  ptr = (void *) (*(((char **)ptr) + lc));
	}
    }
  if(i == c->Ndereference-1)
    {
      if(i == c->lcdereferencecol)
	{
	  if(c->usereallc)
	    {
	      switch(c->type)
		{
		case VARTOOLS_TYPE_DOUBLE:
		  dptr = (double *) ptr + reallc;
		  break;
		case VARTOOLS_TYPE_FLOAT:
		  fptr = (float *) ptr + reallc;
		  break;
		case VARTOOLS_TYPE_INT:
		  iptr = (int *) ptr + reallc;
		  break;
		case VARTOOLS_TYPE_SHORT:
		  shortptr = (short *) ptr + reallc;
		  break;
		case VARTOOLS_TYPE_LONG:
		  longptr = (long *) ptr + reallc;
		  break;
		case VARTOOLS_TYPE_CHAR:
		  charptr = (char *) ptr + reallc;
		  break;
		case VARTOOLS_TYPE_STRING:
		  cptr = (char **) ptr + reallc;
		  break;
		}
	    }
	  else
	    {
	      switch(c->type)
		{
		case VARTOOLS_TYPE_DOUBLE:
		  dptr = (double *) ptr + lc;
		  break;
		case VARTOOLS_TYPE_FLOAT:
		  fptr = (float *) ptr + lc;
		  break;
		case VARTOOLS_TYPE_INT:
		  iptr = (int *) ptr + lc;
		  break;
		case VARTOOLS_TYPE_SHORT:
		  shortptr = (short *) ptr + lc;
		  break;
		case VARTOOLS_TYPE_LONG:
		  longptr = (long *) ptr + lc;
		  break;
		case VARTOOLS_TYPE_CHAR:
		  charptr = (char *) ptr + lc;
		  break;
		case VARTOOLS_TYPE_STRING:
		  cptr = (char **) ptr + lc;
		  break;
		}
	    }
	}
      else
	{
	  switch(c->type)
	    {
	    case VARTOOLS_TYPE_DOUBLE:
	      dptr = (double *) ptr + c->dereference[i];
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      fptr = (float *) ptr + c->dereference[i];
	      break;
	    case VARTOOLS_TYPE_INT:
	      iptr = (int *) ptr + c->dereference[i];
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      shortptr = (short *) ptr + c->dereference[i];
	      break;
	    case VARTOOLS_TYPE_LONG:
	      longptr = (long *) ptr + c->dereference[i];
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      charptr = (char *) ptr + c->dereference[i];
	      break;
	    case VARTOOLS_TYPE_STRING:
	      cptr = (char **) ptr + c->dereference[i];
	      break;
	    }
	}
    }
  switch(outtype)
    {
    case VARTOOLS_TYPE_DOUBLE:
      switch(c->type)
	{
	case VARTOOLS_TYPE_DOUBLE:
	  *((double *) outvalue) = *(dptr);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  *((double *) outvalue) = (double) (*(fptr));
	  break;
	case VARTOOLS_TYPE_INT:
	  *((double *) outvalue) = (double) (*(iptr));
	  break;
	case VARTOOLS_TYPE_SHORT:
	  *((double *) outvalue) = (double) (*(shortptr));
	  break;
	case VARTOOLS_TYPE_LONG:
	  *((double *) outvalue) = (double) (*(longptr));
	  break;
	case VARTOOLS_TYPE_CHAR:
	  *((double *) outvalue) = (double) (*(charptr));
	  break;
	case VARTOOLS_TYPE_STRING:
	  sscanf(*cptr,"%lf",(double *) outvalue);
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      break;
    case VARTOOLS_TYPE_STRING:
      va_start(varlist,outvalue);
      size_out_string = va_arg(varlist,int);
      switch(c->type)
	{
	case VARTOOLS_TYPE_STRING:
	  memcpy(outvalue,(*cptr),MIN_(size_out_string,c->stringsize));
	  break;
	case VARTOOLS_TYPE_DOUBLE:
	  sprintf(((char *) outvalue),"%f",(*(dptr)));
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  sprintf(((char *) outvalue),"%f",(*(fptr)));
	  break;
	case VARTOOLS_TYPE_INT:
	  sprintf(((char *) outvalue),"%d",(*(iptr)));
	  break;
	case VARTOOLS_TYPE_SHORT:
	  sprintf(((char *) outvalue),"%d",(*(shortptr)));
	  break;
	case VARTOOLS_TYPE_LONG:
	  sprintf(((char *) outvalue),"%ld",(*(longptr)));
	  break;
	case VARTOOLS_TYPE_CHAR:
	  sprintf(((char *) outvalue),"%c",(*(charptr)));
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      break;
    case VARTOOLS_TYPE_INT:
      switch(c->type)
	{
	case VARTOOLS_TYPE_INT:
	  *((int *) outvalue) = *(iptr);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  *((int *) outvalue) = (int) (*(shortptr));
	  break;
	case VARTOOLS_TYPE_LONG:
	  *((int *) outvalue) = (int) (*(longptr));
	  break;
	case VARTOOLS_TYPE_CHAR:
	  *((int *) outvalue) = (int) (*(charptr));
	  break;
	case VARTOOLS_TYPE_DOUBLE:
	  *((int *) outvalue) = (int) (*(dptr));
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  *((int *) outvalue) = (int) (*(fptr));
	  break;
	case VARTOOLS_TYPE_STRING:
	  sscanf((*cptr),"%d",((int *) outvalue));
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      break;
    case VARTOOLS_TYPE_FLOAT:
      switch(c->type)
	{
	case VARTOOLS_TYPE_FLOAT:
	  *((float *) outvalue) = *(fptr);
	  break;
	case VARTOOLS_TYPE_DOUBLE:
	  *((float *) outvalue) = (float) (*(dptr));
	  break;
	case VARTOOLS_TYPE_INT:
	  *((float *) outvalue) = (float) (*(iptr));
	  break;
	case VARTOOLS_TYPE_SHORT:
	  *((float *) outvalue) = (float) (*(shortptr));
	  break;
	case VARTOOLS_TYPE_LONG:
	  *((float *) outvalue) = (float) (*(longptr));
	  break;
	case VARTOOLS_TYPE_CHAR:
	  *((float *) outvalue) = (float) (*(charptr));
	  break;
	case VARTOOLS_TYPE_STRING:
	  sscanf((*cptr),"%f",((float *) outvalue));
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      break;
    case VARTOOLS_TYPE_SHORT:
      switch(c->type)
	{
	case VARTOOLS_TYPE_SHORT:
	  *((short *) outvalue) = *(shortptr);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  *((short *) outvalue) = (short) (*(fptr));
	  break;
	case VARTOOLS_TYPE_DOUBLE:
	  *((short *) outvalue) = (short) (*(dptr));
	  break;
	case VARTOOLS_TYPE_INT:
	  *((short *) outvalue) = (short) (*(iptr));
	  break;
	case VARTOOLS_TYPE_LONG:
	  *((short *) outvalue) = (short) (*(longptr));
	  break;
	case VARTOOLS_TYPE_CHAR:
	  *((short *) outvalue) = (short) (*(charptr));
	  break;
	case VARTOOLS_TYPE_STRING:
	  sscanf((*cptr),"%hd",((short *) outvalue));
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      break;
    case VARTOOLS_TYPE_LONG:
      switch(c->type)
	{
	case VARTOOLS_TYPE_LONG:
	  *((long *) outvalue) = *(longptr);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  *((long *) outvalue) = (long) (*(fptr));
	  break;
	case VARTOOLS_TYPE_DOUBLE:
	  *((long *) outvalue) = (long) (*(dptr));
	  break;
	case VARTOOLS_TYPE_INT:
	  *((long *) outvalue) = (long) (*(iptr));
	  break;
	case VARTOOLS_TYPE_SHORT:
	  *((long *) outvalue) = (long) (*(shortptr));
	  break;
	case VARTOOLS_TYPE_CHAR:
	  *((long *) outvalue) = (long) (*(charptr));
	  break;
	case VARTOOLS_TYPE_STRING:
	  sscanf((*cptr),"%ld",((long *) outvalue));
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      break;
    case VARTOOLS_TYPE_CHAR:
      switch(c->type)
	{
	case VARTOOLS_TYPE_CHAR:
	  *((char *) outvalue) = *(charptr);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  *((char *) outvalue) = (char) (*(fptr));
	  break;
	case VARTOOLS_TYPE_DOUBLE:
	  *((char *) outvalue) = (char) (*(dptr));
	  break;
	case VARTOOLS_TYPE_INT:
	  *((char *) outvalue) = (char) (*(iptr));
	  break;
	case VARTOOLS_TYPE_SHORT:
	  *((char *) outvalue) = (char) (*(shortptr));
	  break;
	case VARTOOLS_TYPE_LONG:
	  *((char *) outvalue) = (char) (*(longptr));
	  break;
	case VARTOOLS_TYPE_STRING:
	  sscanf((*cptr),"%c",((char *) outvalue));
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      break;
    default:
      error(ERR_BADTYPE);
      break;
    }
}

void setoutcolumnvalue(OutColumn *c, int lc, int reallc, int intype, void *invalue, ...)
{
  va_list varlist;
  void *ptr;
  double *dptr;
  char **cptr;
  int *iptr;
  float *fptr;
  short *shortptr;
  long *longptr;
  char *charptr;
  double *indptr;
  char **incptr;
  int *iniptr;
  float *infptr;
  short *inshortptr;
  long *inlongptr;
  char *incharptr;
  int i;
  int size_out_string;
  ptr = (void *) (*((char ***) c->ptr));
  for(i=0;i<c->Ndereference-1;i++)
    {
      if(i != c->lcdereferencecol)
	ptr = (void *) (*(((char **)ptr) + c->dereference[i]));
      else if(c->usereallc)
	{
	  ptr = (void *) (*(((char **)ptr) + reallc));
	}
      else
	{
	  ptr = (void *) (*(((char **)ptr) + lc));
	}
    }
  if(i == c->Ndereference-1)
    {
      if(i == c->lcdereferencecol)
	{
	  if(c->usereallc)
	    {
	      switch(c->type)
		{
		case VARTOOLS_TYPE_DOUBLE:
		  dptr = (double *) ptr + reallc;
		  break;
		case VARTOOLS_TYPE_FLOAT:
		  fptr = (float *) ptr + reallc;
		  break;
		case VARTOOLS_TYPE_INT:
		  iptr = (int *) ptr + reallc;
		  break;
		case VARTOOLS_TYPE_SHORT:
		  shortptr = (short *) ptr + reallc;
		  break;
		case VARTOOLS_TYPE_LONG:
		  longptr = (long *) ptr + reallc;
		  break;
		case VARTOOLS_TYPE_CHAR:
		  charptr = (char *) ptr + reallc;
		  break;
		case VARTOOLS_TYPE_STRING:
		  cptr = (char **) ptr + reallc;
		  break;
		}
	    }
	  else
	    {
	      switch(c->type)
		{
		case VARTOOLS_TYPE_DOUBLE:
		  dptr = (double *) ptr + lc;
		  break;
		case VARTOOLS_TYPE_FLOAT:
		  fptr = (float *) ptr + lc;
		  break;
		case VARTOOLS_TYPE_INT:
		  iptr = (int *) ptr + lc;
		  break;
		case VARTOOLS_TYPE_SHORT:
		  shortptr = (short *) ptr + lc;
		  break;
		case VARTOOLS_TYPE_LONG:
		  longptr = (long *) ptr + lc;
		  break;
		case VARTOOLS_TYPE_CHAR:
		  charptr = (char *) ptr + lc;
		  break;
		case VARTOOLS_TYPE_STRING:
		  cptr = (char **) ptr + lc;
		  break;
		}
	    }
	}
      else
	{
	  switch(c->type)
	    {
	    case VARTOOLS_TYPE_DOUBLE:
	      dptr = (double *) ptr + c->dereference[i];
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      fptr = (float *) ptr + c->dereference[i];
	      break;
	    case VARTOOLS_TYPE_INT:
	      iptr = (int *) ptr + c->dereference[i];
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      shortptr = (short *) ptr + c->dereference[i];
	      break;
	    case VARTOOLS_TYPE_LONG:
	      longptr = (long *) ptr + c->dereference[i];
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      charptr = (char *) ptr + c->dereference[i];
	      break;
	    case VARTOOLS_TYPE_STRING:
	      cptr = (char **) ptr + c->dereference[i];
	      break;
	    }
	}
    }
  switch(intype)
    {
    case VARTOOLS_TYPE_DOUBLE:
      switch(c->type)
	{
	case VARTOOLS_TYPE_DOUBLE:
	  *(dptr) = *((double *) invalue);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  *(fptr) = (float) *((double *) invalue);
	  break;
	case VARTOOLS_TYPE_INT:
	  *(iptr) = (int) *((double *) invalue);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  *(shortptr) = (short) *((double *) invalue);
	  break;
	case VARTOOLS_TYPE_LONG:
	  *(longptr) = (long) *((double *) invalue);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  *(charptr) = (char) *((double *) invalue);
	  break;
	case VARTOOLS_TYPE_STRING:
	  sprintf(*cptr,"%lf",*((double *) invalue));
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      break;
    case VARTOOLS_TYPE_STRING:
      va_start(varlist,invalue);
      size_out_string = va_arg(varlist,int);
      switch(c->type)
	{
	case VARTOOLS_TYPE_STRING:
	  memcpy((*cptr),((char *) invalue),MIN_(size_out_string,c->stringsize));
	  break;
	case VARTOOLS_TYPE_DOUBLE:
	  *(dptr) = atof((char *) invalue);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  *(fptr) = (float) atof((char *) invalue);
	  break;
	case VARTOOLS_TYPE_INT:
	  *(iptr) = atoi((char *) invalue);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  *(shortptr) = (short) atoi((char *) invalue);
	  break;
	case VARTOOLS_TYPE_LONG:
	  *(longptr) = atol((char *) invalue);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  *(charptr) = ((char *) invalue)[0];
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      break;
    case VARTOOLS_TYPE_INT:
      switch(c->type)
	{
	case VARTOOLS_TYPE_INT:
	  *(iptr) = *((int *) invalue);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  *(shortptr) = (short) *((int *) invalue);
	  break;
	case VARTOOLS_TYPE_LONG:
	  *(longptr) = (long) *((int *) invalue);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  *(charptr) = (char) *((int *) invalue);
	  break;
	case VARTOOLS_TYPE_DOUBLE:
	  *(dptr) = (double) *((int *) invalue);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  *(fptr) = (float) *((int *) invalue);
	  break;
	case VARTOOLS_TYPE_STRING:
	  sprintf((*cptr),"%d",*((int *) invalue));
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      break;
    case VARTOOLS_TYPE_FLOAT:
      switch(c->type)
	{
	case VARTOOLS_TYPE_FLOAT:
	  *(fptr) = *((float *) invalue);
	  break;
	case VARTOOLS_TYPE_DOUBLE:
	  *(dptr) = (double) *((float *) invalue);
	  break;
	case VARTOOLS_TYPE_INT:
	  *(iptr) = (int) *((float *) invalue);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  *(shortptr) = (short) *((float *) invalue);
	  break;
	case VARTOOLS_TYPE_LONG:
	  *(longptr) = (long) *((float *) invalue);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  *(charptr) = (char) *((float *) invalue);
	  break;
	case VARTOOLS_TYPE_STRING:
	  sprintf((*cptr),"%f",*((float *) invalue));
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      break;
    case VARTOOLS_TYPE_SHORT:
      switch(c->type)
	{
	case VARTOOLS_TYPE_SHORT:
	  *(shortptr) = *((short *) invalue);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  *(fptr) = (float) *((short *) invalue);
	  break;
	case VARTOOLS_TYPE_DOUBLE:
	  *(dptr) = (double) *((short *) invalue);
	  break;
	case VARTOOLS_TYPE_INT:
	  *(iptr) = (int) *((short *) invalue);
	  break;
	case VARTOOLS_TYPE_LONG:
	  *(longptr) = (long) *((short *) invalue);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  *(charptr) = (char) *((short *) invalue);
	  break;
	case VARTOOLS_TYPE_STRING:
	  sprintf((*cptr),"%d",*((short *) invalue));
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      break;
    case VARTOOLS_TYPE_LONG:
      switch(c->type)
	{
	case VARTOOLS_TYPE_LONG:
	  *(longptr) = *((long *) invalue);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  *(fptr) = (float) *((long *) invalue);
	  break;
	case VARTOOLS_TYPE_DOUBLE:
	  *(dptr) = (double) *((long *) invalue);
	  break;
	case VARTOOLS_TYPE_INT:
	  *(iptr) = (int) *((long *) invalue);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  *(shortptr) = (short) *((long *) invalue);
	  break;
	case VARTOOLS_TYPE_CHAR:
	  *(charptr) = (char) *((long *) invalue);
	  break;
	case VARTOOLS_TYPE_STRING:
	  sprintf((*cptr),"%ld",*((long *) invalue));
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      break;
    case VARTOOLS_TYPE_CHAR:
      switch(c->type)
	{
	case VARTOOLS_TYPE_CHAR:
	  *(charptr) = *((char *) invalue);
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  *(fptr) = (float) *((char *) invalue);
	  break;
	case VARTOOLS_TYPE_DOUBLE:
	  *(dptr) = (double) *((char *) invalue);
	  break;
	case VARTOOLS_TYPE_INT:
	  *(iptr) = (int) *((char *) invalue);
	  break;
	case VARTOOLS_TYPE_SHORT:
	  *(shortptr) = (short) *((char *) invalue);
	  break;
	case VARTOOLS_TYPE_LONG:
	  *(longptr) = (long) *((char *) invalue);
	  break;
	case VARTOOLS_TYPE_STRING:
	  sprintf((*cptr),"%c",*((char *) invalue));
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
      break;
    default:
      error(ERR_BADTYPE);
      break;
    }
}

void zerooutcolumnvalue(OutColumn *c, int lc, int reallc)
{
  void *ptr;
  double *dptr;
  char **cptr;
  int *iptr;
  float *fptr;
  short *shortptr;
  long *longptr;
  char *charptr;
  int i;
  int size_out_string;
  ptr = (void *) (*((char ***) c->ptr));
  for(i=0;i<c->Ndereference-1;i++)
    {
      if(i != c->lcdereferencecol)
	ptr = (void *) (*(((char **)ptr) + c->dereference[i]));
      else if(c->usereallc)
	{
	  ptr = (void *) (*(((char **)ptr) + reallc));
	}
      else
	{
	  ptr = (void *) (*(((char **)ptr) + lc));
	}
    }
  if(i == c->Ndereference-1)
    {
      if(i == c->lcdereferencecol)
	{
	  if(c->usereallc)
	    {
	      switch(c->type)
		{
		case VARTOOLS_TYPE_DOUBLE:
		  dptr = (double *) ptr + reallc;
		  *dptr = 0.0;
		  break;
		case VARTOOLS_TYPE_FLOAT:
		  fptr = (float *) ptr + reallc;
		  *fptr = 0.0;
		  break;
		case VARTOOLS_TYPE_INT:
		  iptr = (int *) ptr + reallc;
		  *iptr = 0;
		  break;
		case VARTOOLS_TYPE_SHORT:
		  shortptr = (short *) ptr + reallc;
		  *shortptr = 0;
		  break;
		case VARTOOLS_TYPE_LONG:
		  longptr = (long *) ptr + reallc;
		  *longptr = 0;
		  break;
		case VARTOOLS_TYPE_CHAR:
		  charptr = (char *) ptr + reallc;
		  *charptr = '0';
		  break;
		case VARTOOLS_TYPE_STRING:
		  cptr = (char **) ptr + reallc;
		  sprintf(*cptr, "0");
		  break;
		}
	    }
	  else
	    {
	      switch(c->type)
		{
		case VARTOOLS_TYPE_DOUBLE:
		  dptr = (double *) ptr + lc;
		  *dptr = 0.0;
		  break;
		case VARTOOLS_TYPE_FLOAT:
		  fptr = (float *) ptr + lc;
		  *fptr = 0.0;
		  break;
		case VARTOOLS_TYPE_INT:
		  iptr = (int *) ptr + lc;
		  *iptr = 0;
		  break;
		case VARTOOLS_TYPE_SHORT:
		  shortptr = (short *) ptr + lc;
		  *shortptr = 0;
		  break;
		case VARTOOLS_TYPE_LONG:
		  longptr = (long *) ptr + lc;
		  *longptr = 0;
		  break;
		case VARTOOLS_TYPE_CHAR:
		  charptr = (char *) ptr + lc;
		  *charptr = '0';
		  break;
		case VARTOOLS_TYPE_STRING:
		  cptr = (char **) ptr + lc;
		  sprintf(*cptr,"0");
		  break;
		}
	    }
	}
      else
	{
	  switch(c->type)
	    {
	    case VARTOOLS_TYPE_DOUBLE:
	      dptr = (double *) ptr + c->dereference[i];
	      *dptr = 0.0;
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      fptr = (float *) ptr + c->dereference[i];
	      *fptr = 0.0;
	      break;
	    case VARTOOLS_TYPE_INT:
	      iptr = (int *) ptr + c->dereference[i];
	      *iptr = 0;
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      shortptr = (short *) ptr + c->dereference[i];
	      *shortptr = 0;
	      break;
	    case VARTOOLS_TYPE_LONG:
	      longptr = (long *) ptr + c->dereference[i];
	      *longptr = 0;
	      break;
	    case VARTOOLS_TYPE_CHAR:
	      charptr = (char *) ptr + c->dereference[i];
	      *charptr = '0';
	      break;
	    case VARTOOLS_TYPE_STRING:
	      cptr = (char **) ptr + c->dereference[i];
	      sprintf(*cptr,"0");
	      break;
	    }
	}
    }
}


void CreateOutputColumns(ProgramData *p, Command *c, int Ncommands)
{
  int i, j, k, l, m;
  char tmpstring[MAXLEN];
  addcolumn(p, -1, VARTOOLS_TYPE_STRING, MAXLEN, &(p->lcnames), "%s", 1, 1, 0, 0, "Name");
  for(l=0;l<Ncommands;l++)
    {
      switch(c[l].cnum)
	{
#ifdef DYNAMICLIB
	case CNUM_USERCOMMAND:
	  CreateOutputColumns_UserCommand(p, c, l);
	  break;
#endif
#ifdef _HAVE_GSL
	case CNUM_ADDNOISE:
	  if(c[l].AddNoise->gammaval_type == PERTYPE_SPECIFIED)
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AddNoise->gammaval), "%9.5f", 2, 1, 0, 0, 0, "AddNoise_Gamma_%d_%d", 1, l);
	  if(c[l].AddNoise->sig_r_type == PERTYPE_SPECIFIED)
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AddNoise->sig_r), "%9.5f", 2, 1, 0, 0, 0, "AddNoise_Sig_Red_%d_%d", 1, l);
	  if(c[l].AddNoise->sig_w_type == PERTYPE_SPECIFIED)
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AddNoise->sig_w), "%9.5f", 2, 1, 0, 0, 0, "AddNoise_Sig_White_%d_%d", 1, l);
	  if(c[l].AddNoise->rho_r_type == PERTYPE_SPECIFIED)
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AddNoise->rho_r), "%9.5f", 2, 1, 0, 0, 0, "AddNoise_Rho_%d_%d", 1, l);
	  if(c[l].AddNoise->nu_r_type == PERTYPE_SPECIFIED)
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AddNoise->nu_r), "%9.5f", 2, 1, 0, 0, 0, "AddNoise_Nu_%d_%d", 1, l);
	  break;
#endif
	case CNUM_CLIP:
	  addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].Clip->Nclip), "%5d", 1, 0, 0, 0, "Nclip_%d", l);
	  break;
	case CNUM_ENSEMBLERESCALESIG:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Ensemblerescalesig->rescalefactor), "%9.5f", 1, 0, 0, 0, "SigmaRescaleFactor_%d", l);
	  break;
	case CNUM_RESCALESIG:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Rescalesig->rescalefactor), "%9.5f", 1, 0, 0, 0, "SigmaRescaleFactor_%d", l);
	  break;
	case CNUM_RESTRICTTIMES:
	  if(c[l].RestrictTimes->restricttype == VARTOOLS_RESTRICTTIMES_JDRANGE) {
	    if(c[l].RestrictTimes->minJDtype != PERTYPE_SPECIFIED) {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].RestrictTimes->minJD), "%.17g", 1, 0, 0, 0, "RestrictTimes_MinJD_%d", l);
	    } 
	    else {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].RestrictTimes->minJD), "%.17g", 1, 1, 0, 0, "RestrictTimes_MinJD_%d", l);
	    } 
	    if(c[l].RestrictTimes->maxJDtype != PERTYPE_SPECIFIED) {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].RestrictTimes->maxJD), "%.17g", 1, 0, 0, 0, "RestrictTimes_MaxJD_%d", l);
	    }
	    else {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].RestrictTimes->maxJD), "%.17g", 1, 1, 0, 0, "RestrictTimes_MaxJD_%d", l);
	    }
	  }
	  break;
	case CNUM_CHI2_NOBIN:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Chi2_NoBin->chi2val), "%12.5f", 1, 0, 0, 0, "Chi2_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Chi2_NoBin->wtave), "%9.5f", 1, 0, 0, 0, "Weighted_Mean_Mag_%d", l);
	  break;
	case CNUM_CHI2_BIN:
	  for(i=0;i<c[l].Chi2_Bin->Nbin;i++)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Chi2_Bin->chi2binvals), "%12.5f", 2, 0, 0, 0, i, "Chi2Bin_%d.%d_%d", (int) (c[l].Chi2_Bin->bintimes[i]*MINUTESPERDAY), (int) (100.0 * (c[l].Chi2_Bin->bintimes[i] * MINUTESPERDAY - (int) (c[l].Chi2_Bin->bintimes[i] * MINUTESPERDAY))), l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Chi2_Bin->wtavebin), "%9.5f", 2, 0, 0, 0, i, "Weight_Mean_Mag_Bin_%d.%d_%d", (int) (c[l].Chi2_Bin->bintimes[i]*MINUTESPERDAY), (int) (100.0 * (c[l].Chi2_Bin->bintimes[i] * MINUTESPERDAY - (int) (c[l].Chi2_Bin->bintimes[i] * MINUTESPERDAY))),l);
	    }
	  break;
	case CNUM_CHANGEERROR:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Changeerror->ave), "%9.5f", 1, 0, 0, 0, "Mean_Mag_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Changeerror->rmsval), "%9.5f", 1, 0, 0, 0, "RMS_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].Changeerror->ngood), "%5d", 1, 0, 0, 0, "Npoints_%d", l);
	  break;
	case CNUM_RMS_NOBIN:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].RMS_NoBin->ave), "%9.5f", 1, 0, 0, 0, "Mean_Mag_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].RMS_NoBin->rmsval), "%9.5f", 1, 0, 0, 0, "RMS_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].RMS_NoBin->rmsthy), "%9.5f", 1, 0, 0, 0, "Expected_RMS_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].RMS_NoBin->ngood), "%5d", 1, 0, 0, 0, "Npoints_%d", l);
	  break;
	case CNUM_RMS_BIN:
	  for(i=0;i<c[l].RMS_Bin->Nbin;i++)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].RMS_Bin->rmsbinvals), "%9.5f", 2, 0, 0, 0, i, "RMSBin_%d.%d_%d", (int) (c[l].RMS_Bin->bintimes[i] * MINUTESPERDAY), (int) (100.0 * (c[l].RMS_Bin->bintimes[i] * MINUTESPERDAY - (int) (c[l].RMS_Bin->bintimes[i] * MINUTESPERDAY))),l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].RMS_Bin->rmsthybin), "%9.5f", 2, 0, 0, 0, i, "Expected_RMS_Bin_%d.%d_%d", (int) (c[l].RMS_Bin->bintimes[i] * MINUTESPERDAY), (int) (100.0 * (c[l].RMS_Bin->bintimes[i] * MINUTESPERDAY - (int) (c[l].RMS_Bin->bintimes[i] * MINUTESPERDAY))),l);
	    }
	  break;
	case CNUM_JSTET:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Jstet->jst), "%9.5f", 1, 0, 0, 0, "Jstet_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Jstet->kur), "%9.5f", 1, 0, 0, 0, "Kurtosis_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Jstet->lst), "%9.5f", 1, 0, 0, 0, "Lstet_%d", l);
	  break;
	case CNUM_ALARM:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Alarm->alarmvals), "%9.5f", 1, 0, 0, 0, "Alarm_%d", l);
	  break;
	case CNUM_FINDBLENDS:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].FindBlends->periods), "%14.8f", 2, 0, 0, 0, 0, "FindBlends_Period_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_STRING, MAXLEN, &(c[l].FindBlends->varblendnames), "%s", 1, 0, 0, 0, "FindBlends_LCname_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].FindBlends->blendamps), "%9.5g", 1, 0, 0, 0, "FindBlends_FluxAmp_%d",l);
	  break;
	case CNUM_AOV:
	  for(i=1;i<=c[l].Aov->Npeaks;i++)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Aov->peakperiods), "%14.8f", 2, 0, 0, 0, i-1, "Period_%d_%d",i,l);
	      if(c[l].Aov->uselog)
		addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Aov->peakvalues), "%9.5f", 2, 0, 0, 0, i-1, "AOV_LOGSNR_%d_%d",i,l);
	      else
		{
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Aov->peakvalues), "%9.5f", 2, 0, 0, 0, i-1, "AOV_%d_%d",i,l);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Aov->peakSNR), "%9.5f", 2, 0, 0, 0, i-1, "AOV_SNR_%d_%d",i,l);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Aov->peakFAP), "%9.5f", 2, 0, 0, 0, i-1, "AOV_NEG_LN_FAP_%d_%d",i,l);
		}
	      if(c[l].Aov->whiten && c[l].Aov->uselog)
		{
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Aov->aveaov_whiten), "%9.5f", 2, 0, 0, 0, i-1, "Mean_lnAOV_%d_%d",i,l);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Aov->rmsaov_whiten), "%9.5f", 2, 0, 0, 0, i-1, "RMS_lnAOV_%d_%d",i,l);
		}
	    }
	  if(c[l].Aov->fixperiodSNR)
	    {
	      if(c[l].Aov->uselog)
		{
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Aov->fixperiodSNR_periods), "%14.8f", 2, 0, 0, 0, 0, "PeriodFix_%d",l);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Aov->fixperiodSNR_peakvalues), "%10.5f", 1, 0, 0, 0, "AOV_LOGSNR_PeriodFix_%d",l);
		}
	      else
		{
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Aov->fixperiodSNR_periods), "%14.8f", 2, 0, 0, 0, 0, "PeriodFix_%d",l);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Aov->fixperiodSNR_peakvalues), "%10.5f", 1, 0, 0, 0, "AOV_PeriodFix_%d",l);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Aov->fixperiodSNR_peakSNR), "%10.5f", 1, 0, 0, 0, "AOV_SNR_PeriodFix_%d",l);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Aov->fixperiodSNR_peakFAP), "%10.5f", 1, 0, 0, 0, "AOV_NEG_LN_FAP_PeriodFix_%d", l);
		}
	    }
	  if(!c[l].Aov->whiten && c[l].Aov->uselog)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Aov->aveaov), "%9.5f", 1, 0, 0, 0, "Mean_lnAOV_%d", l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Aov->rmsaov), "%9.5f", 1, 0, 0, 0, "RMS_lnAOV_%d", l);
	    }
	  break;
	case CNUM_HARMAOV:
	  for(i=1;i<=c[l].AovHarm->Npeaks;i++)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AovHarm->peakperiods), "%14.8f", 2, 0, 0, 0, i-1, "Period_%d_%d",i,l);
	      if(c[l].AovHarm->Nharm > 0)
		addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AovHarm->peakvalues), "%9g", 2, 0, 0, 0, i-1, "AOV_HARM_%d_%d",i,l);
	      else
		{
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AovHarm->peakvalues), "%9g", 2, 0, 0, 0, i-1, "AOV_HARM_NEG_LOG_FAP_%d_%d",i,l);
		  addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].AovHarm->peakNharm), "%2d", 2, 0, 0, 0, i-1, "AOV_HARM_NHARM_%d_%d",i,l);
		}
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AovHarm->peakSNR), "%9g", 2, 0, 0, 0, i-1, "AOV_HARM_SNR_%d_%d",i,l);
	      if(c[l].AovHarm->Nharm > 0)
		addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AovHarm->peakFAP), "%9g", 2, 0, 0, 0, i-1, "AOV_HARM_NEG_LOG_FAP_%d_%d",i,l);
	      if(c[l].AovHarm->whiten)
		{
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AovHarm->aveaov_whiten), "%9g", 2, 0, 0, 0, i-1, "Mean_AOV_HARM_%d_%d",i,l);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AovHarm->rmsaov_whiten), "%9g", 2, 0, 0, 0, i-1, "RMS_AOV_HARM_%d_%d",i,l);
		}
	    }
	  if(c[l].AovHarm->fixperiodSNR)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AovHarm->fixperiodSNR_periods), "%14.8f", 2, 0, 0, 0, 0, "PeriodFix_%d",l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AovHarm->fixperiodSNR_peakvalues), "%10.5f", 1, 0, 0, 0, "AOV_HARM_PeriodFix_%d",l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AovHarm->fixperiodSNR_peakSNR), "%10.5f", 1, 0, 0, 0, "AOV_HARM_SNR_PeriodFix_%d",l);
	      if(c[l].AovHarm->Nharm > 0)
		addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AovHarm->fixperiodSNR_peakFAP), "%10.5f", 1, 0, 0, 0, "AOV_HARM_NEG_LN_FAP_PeriodFix_%d", l);
	    }
	  if(!c[l].AovHarm->whiten)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AovHarm->aveaov), "%9g", 1, 0, 0, 0, "Mean_AOV_HARM_%d", l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].AovHarm->rmsaov), "%9g", 1, 0, 0, 0, "RMS_AOV_HARM_%d", l);
	    }
	  break;
	case CNUM_LS:
	  for(i=1;i<=c[l].Ls->Npeaks;i++)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Ls->peakperiods), "%14.8f", 2, 0, 0, 0, i-1, "LS_Period_%d_%d",i,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Ls->peakFAP), "%10.5f", 2, 0, 0, 0, i-1, "Log10_LS_Prob_%d_%d",i,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Ls->peakvalues), "%10.5f", 2, 0, 0, 0, i-1, "LS_Periodogram_Value_%d_%d",i,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Ls->SNRvalues), "%10.5f", 2, 0, 0, 0, i-1, "LS_SNR_%d_%d",i,l);
	    }
	  if(c[l].Ls->fixperiodSNR)
	    {
	      if(c[l].Ls->fixperiodSNR_pertype != PERTYPE_SPECIFIED)
		addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Ls->fixperiodSNR_periods), "%14.8f", 2, 0, 0, 0, 0, "LS_PeriodFix_%d",l);
	      else
		addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Ls->fixperiodSNR_periods), "%14.8f", 2, 1, 0, 0, 0, "LS_PeriodFix_%d",l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Ls->fixperiodSNR_FAPvalues), "%10.5f", 1, 0, 0, 0, "Log10_LS_Prob_PeriodFix_%d",l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Ls->fixperiodSNR_peakvalues), "%10.5f", 1, 0, 0, 0, "LS_Periodogram_Value_PeriodFix_%d",l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Ls->fixperiodSNR_SNRvalues), "%10.5f", 1, 0, 0, 0, "LS_SNR_PeriodFix_%d",l);
	    }
	  break;
	case CNUM_GETLSAMPTHRESH:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].GetLSAmpThresh->ampthresh_scale), "%9.5f", 1, 0, 0, 0, "LS_AmplitudeScaleFactor_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].GetLSAmpThresh->amp), "%9.5f", 1, 0, 0, 0, "LS_MinimumAmplitude_%d",l);
	  break;
	case CNUM_DECORR:
	  k = 0;
	  if(c[l].Decorr->zeropointterm)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Decorr->b), "%20.13f", 2, 0, 0, 0, 0, "Decorr_constant_term_%d", l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Decorr->b_err), "%20.13f", 2, 0, 0, 0, 0, "Decorr_constant_term_err_%d", l);
	      m = 1;
	    }
	  else
	    m = 0;
	  for(i=1;i<=c[l].Decorr->N_globalterms;i++)
	    {
	      for(j=1;j<=c[l].Decorr->order[k];j++)
		{
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Decorr->b), "%20.13f", 2, 0, 0, 0, m, "Global_%d_coeff_%d_%d",i,j,l);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Decorr->b_err), "%20.13f", 2, 0, 0, 0, m, "Global_%d_coeff_err_%d_%d",i,j,l);
		  m++;
		}
	      k++;
	    }
	  for(i=0;i<c[l].Decorr->N_lcterms;i++)
	    {
	      for(j=1;j<=c[l].Decorr->order[k];j++)
		{
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Decorr->b), "%20.13f", 2, 0, 0, 0, m, "LCColumn_%d_coeff_%d_%d",c[l].Decorr->lc_columns[i],j,l);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Decorr->b_err), "%20.13f", 2, 0, 0, 0, m, "LCColumn_%d_coeff_err_%d_%d",c[l].Decorr->lc_columns[i],j,l);
		  m++;
		}
	      k++;
	    }
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Decorr->chi2val), "%12.5f", 1, 0, 0, 0, "Decorr_chi2_%d", l);
	  break;
	case CNUM_LINFIT:
	  for(j=1;j<=c[l].Linfit->Nparams;j++)
	    {
	      sprintf(tmpstring,"Linfit_%s_%%d", c[l].Linfit->paramnames[j-1]);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Linfit->param_outvals), "%.17g", 2, 0, 0, 0, j-1, tmpstring, l, 0);

	      sprintf(tmpstring,"Linfit_err%s_%%d", c[l].Linfit->paramnames[j-1]);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Linfit->param_uncertainties), "%.17g", 2, 0, 0, 0, j-1, tmpstring, l, 0);
	    }
	  break;
	case CNUM_NONLINFIT:
	  if(c[l].Nonlinfit->fittype == VARTOOLS_NONLINFIT_FITTYPE_AMOEBA) {
	    addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].Nonlinfit->amoeba_isconverged), "%d", 1, 0, 0, 0, "Nonlinfit_HasConverged_%d", l);
	  }
	  for(j=1;j<=c[l].Nonlinfit->Nparams;j++)
	    {
	      sprintf(tmpstring,"Nonlinfit_%s_BestFit_%%d", c[l].Nonlinfit->paramnames[j-1]);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Nonlinfit->param_outvals), "%.17g", 2, 0, 0, 0, j-1, tmpstring, l, 0);
	      if(c[l].Nonlinfit->fittype == VARTOOLS_NONLINFIT_FITTYPE_AMOEBA) {
		sprintf(tmpstring,"Nonlinfit_%s_Err_%%d", c[l].Nonlinfit->paramnames[j-1]);
		addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Nonlinfit->param_uncertainties), "%.17g", 2, 0, 0, 0, j-1, tmpstring, l, 0);
	      }
	    }
	  if(c[l].Nonlinfit->uselinfit)
	    {
	      for(m=0; m < c[l].Nonlinfit->linfit->Nparams;m++) {
		sprintf(tmpstring,"Nonlinfit_%s_BestFit_%%d", c[l].Nonlinfit->linfit->paramnames[m]);
		addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Nonlinfit->linfit->param_outvals), "%.17g", 2, 0, 0, 0, m, tmpstring, l, 0);
		if(c[l].Nonlinfit->fittype == VARTOOLS_NONLINFIT_FITTYPE_AMOEBA) {
		  sprintf(tmpstring,"Nonlinfit_%s_Err_%%d", c[l].Nonlinfit->linfit->paramnames[m]);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Nonlinfit->linfit->param_uncertainties), "%.17g", 2, 0, 0, 0, m, tmpstring, l, 0);
		}

	      }
	    }
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Nonlinfit->chi2out),
		    "%.17g", 1, 0, 0, 0, "Nonlinfit_BestFit_Chi2_%d", l);
	  for(k=0,i=0; i < c[l].Nonlinfit->N_mcmc_chain_expressions; i++)
	    {
	      for(j=0,m=0; j < c[l].Nonlinfit->N_mcmc_chain_stats; j++, k++)
		{
		  sprintf(tmpstring,"Nonlinfit_%s_",c[l].Nonlinfit->mcmc_chain_expr_strings[i]);
		  switch(c[l].Nonlinfit->mcmc_statstocalc[j]) {
		  case VARTOOLS_STATSTYPE_MEAN:
		    sprintf(tmpstring,"%sMEAN_%%d", tmpstring);
		    break;
		  case VARTOOLS_STATSTYPE_WEIGHTEDMEAN:
		    sprintf(tmpstring,"%sWEIGHTEDMEAN_%%d",tmpstring);
		    break;
		  case VARTOOLS_STATSTYPE_MEDIAN:
		    sprintf(tmpstring,"%sMEDIAN_%%d",tmpstring);
		    break;
		  case VARTOOLS_STATSTYPE_MEDIAN_WEIGHT:
		    sprintf(tmpstring,"%sWEIGHTEDMEDIAN_%%d",tmpstring);
		    break;
		  case VARTOOLS_STATSTYPE_STDDEV:
		    sprintf(tmpstring,"%sSTDDEV_%%d",tmpstring);
		    break;
		  case VARTOOLS_STATSTYPE_MEDDEV:
		    sprintf(tmpstring,"%sMEDDEV_%%d",tmpstring);
		    break;
		  case VARTOOLS_STATSTYPE_MEDMEDDEV:
		    sprintf(tmpstring,"%sMEDMEDDEV_%%d",tmpstring);
		    break;
		  case VARTOOLS_STATSTYPE_MAD:
		    sprintf(tmpstring,"%sMAD_%%d",tmpstring);
		    break;
		  case VARTOOLS_STATSTYPE_KURTOSIS:
		    sprintf(tmpstring,"%sKURTOSIS_%%d",tmpstring);
		    break;
		  case VARTOOLS_STATSTYPE_SKEWNESS:
		    sprintf(tmpstring,"%sSKEWNESS_%%d",tmpstring);
		    break;
		  case VARTOOLS_STATSTYPE_PERCENTILE:
		    sprintf(tmpstring,"%sPCT%.2f_%%d",tmpstring,c[l].Nonlinfit->pctval[m]);
		    m++;
		    break;
		  case VARTOOLS_STATSTYPE_PERCENTILE_WEIGHT:
		    sprintf(tmpstring,"%sWPCT%.2f_%%d",tmpstring,c[l].Nonlinfit->pctval[m]);
		    m++;
		    break;
		  case VARTOOLS_STATSTYPE_MAXIMUM:
		    sprintf(tmpstring,"%sMAX_%%d",tmpstring);
		    break;
		  case VARTOOLS_STATSTYPE_MINIMUM:
		    sprintf(tmpstring,"%sMIN_%%d",tmpstring);
		    break;
		  case VARTOOLS_STATSTYPE_SUM:
		    sprintf(tmpstring,"%sSUM_%%d",tmpstring);
		    break;
		  default:
		    break;
		  }
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Nonlinfit->mcmc_statsout), "%.17g", 2, 0, 0, 0, k, tmpstring, l);
		}
	    }
	  break;
	case CNUM_KILLHARM:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->mean), "%9.5f", 1, 0, 0, 0, "Killharm_Mean_Mag_%d", l);
	  for(j=1;j<=c[l].Killharm->Nper;j++)
	    {
	      if(c[l].Killharm->pertype == PERTYPE_SPECIFIED)
		addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->periods), "%14.8f", 2, 1, 0, 0, j-1, "Killharm_Period_%d_%d",j,l);
	      else
		addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->periods), "%14.8f", 2, 0, 0, 0, j-1, "Killharm_Period_%d_%d",j,l);
	      for(i=2;i<=c[l].Killharm->Nsubharm+1;i++)
		{
		  if(c[l].Killharm->outtype == KILLHARM_OUTTYPE_DEFAULT)
		    {
		      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->subharmA), "%9.5f", 3, 0, 0, 0, j-1, i-2,"Killharm_Per%d_Subharm_%d_Sincoeff_%d",j,i,l);
		      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->subharmB), "%9.5f", 3, 0, 0, 0, j-1, i-2,"Killharm_Per%d_Subharm_%d_Coscoeff_%d",j,i,l);
		    }
		  else if(c[l].Killharm->outtype == KILLHARM_OUTTYPE_AMPPHASE || c[l].Killharm->outtype == KILLHARM_OUTTYPE_AMPRADPHASE)
		    {
		      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->subharmA), "%9.5f", 3, 0, 0, 0, j-1, i-2,"Killharm_Per%d_Subharm_Amp_%d_%d",j,i,l);
		      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->subharmB), "%9.5f", 3, 0, 0, 0, j-1, i-2,"Killharm_Per%d_Subharm_Phi_%d_%d",j,i,l);
		    }
		  else if(c[l].Killharm->outtype == KILLHARM_OUTTYPE_RPHI || c[l].Killharm->outtype == KILLHARM_OUTTYPE_RRADPHI)
		    {
		      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->subharmA), "%9.5f", 3, 0, 0, 0, j-1, i-2,"Killharm_Per%d_Subharm_R_%d_1_%d",j,i,l);
		      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->subharmB), "%9.5f", 3, 0, 0, 0, j-1, i-2,"Killharm_Per%d_Subharm_Phi_%d_1_%d",j,i,l);
		    }
		}
	      if(c[l].Killharm->outtype == KILLHARM_OUTTYPE_DEFAULT)
		{
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->fundA), "%9.5f", 2, 0, 0, 0, j-1, "Killharm_Per%d_Fundamental_Sincoeff_%d",j,l);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->fundB), "%9.5f", 2, 0, 0, 0, j-1, "Killharm_Per%d_Fundamental_Coscoeff_%d",j,l);
		}
	      else
		{
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->fundA), "%9.5f", 2, 0, 0, 0, j-1, "Killharm_Per%d_Fundamental_Amp_%d",j,l);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->fundB), "%9.5f", 2, 0, 0, 0, j-1, "Killharm_Per%d_Fundamental_Phi_%d",j,l);
		}
	      for(i=2;i<=c[l].Killharm->Nharm+1;i++)
		{
		  if(c[l].Killharm->outtype == KILLHARM_OUTTYPE_DEFAULT)
		    {
		      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->harmA), "%9.5f", 3, 0, 0, 0, j-1, i-2,"Killharm_Per%d_Harm_%d_Sincoeff_%d",j,i,l);
		      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->harmB), "%9.5f", 3, 0, 0, 0, j-1, i-2,"Killharm_Per%d_Harm_%d_Coscoeff_%d",j,i,l);
		    }
		  else if(c[l].Killharm->outtype == KILLHARM_OUTTYPE_AMPPHASE || c[l].Killharm->outtype == KILLHARM_OUTTYPE_AMPRADPHASE)
		    {
		      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->harmA), "%9.5f", 3, 0, 0, 0, j-1, i-2,"Killharm_Per%d_Harm_Amp_%d_%d",j,i,l);
		      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->harmB), "%9.5f", 3, 0, 0, 0, j-1, i-2,"Killharm_Per%d_Harm_Phi_%d_%d",j,i,l);
		    }
		  else if(c[l].Killharm->outtype == KILLHARM_OUTTYPE_RPHI || c[l].Killharm->outtype == KILLHARM_OUTTYPE_RRADPHI)
		    {
		      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->harmA), "%9.5f", 3, 0, 0, 0, j-1, i-2,"Killharm_Per%d_Harm_R_%d_1_%d",j,i,l);
		      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->harmB), "%9.5f", 3, 0, 0, 0, j-1, i-2,"Killharm_Per%d_Harm_Phi_%d_1_%d",j,i,l);
		    }
		}
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Killharm->amp), "%9.5f", 2, 0, 0, 0, j-1, "Killharm_Per%d_Amplitude_%d",j,l);
	    }
	  break;
	case CNUM_INJECTHARM:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injectharm->periodinject), "%14.8f", 1, 0, 0, 0, "Injectharm_Period_%d", l);
	  for(i=2;i<=c[l].Injectharm->Nsubharm+1;i++)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injectharm->subharm_amp), "%9.5f", 2, 0, 0, 0, i-2, "Injectharm_Subharm_%d_Amp_%d",i,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injectharm->subharm_phase), "%9.5f", 2, 0, 0, 0, i-2, "Injectharm_Subharm_%d_Phase_%d",i,l);
	    }
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injectharm->harm_amp), "%9.5f", 2, 0, 0, 0, 0, "Injectharm_Fundamental_Amp_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injectharm->harm_phase), "%9.5f", 2, 0, 0, 0, 0, "Injectharm_Fundamental_Phase_%d",l);
	  for(i=2;i<=c[l].Injectharm->Nharm+1;i++)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injectharm->harm_amp), "%9.5f", 2, 0, 0, 0, i-1, "Injectharm_Harm_%d_Amp_%d",i,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injectharm->harm_phase), "%9.5f", 2, 0, 0, 0, i-1, "Injectharm_Harm_%d_Phase_%d",i,l);
	    }
	  break;
	case CNUM_INJECTTRANSIT:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injecttransit->paraminject[INJECTTR_IDX_PERIOD]), "%14.8f", 1, 0, 0, 0, "Injecttransit_Period_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injecttransit->paraminject[INJECTTR_IDX_RP]), "%9.5f", 1, 0, 0, 0, "Injecttransit_Rp_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injecttransit->paraminject[INJECTTR_IDX_MP]), "%9.5f", 1, 0, 0, 0, "Injecttransit_Mp_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injecttransit->paraminject[INJECTTR_IDX_PHASE]), "%9.5f", 1, 0, 0, 0, "Injecttransit_phase_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injecttransit->paraminject[INJECTTR_IDX_SINI]), "%9.5f", 1, 0, 0, 0, "Injecttransit_sin_i_%d",l);
	  if(c[l].Injecttransit->eomegatype)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injecttransit->paraminject[INJECTTR_IDX_E]), "%9.5f", 1, 0, 0, 0, "Injecttransit_e_%d",l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injecttransit->paraminject[INJECTTR_IDX_OMEGA]), "%9.5f", 1, 0, 0, 0, "Injecttransit_omega_%d",l);
	    }
	  else
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injecttransit->paraminject[INJECTTR_IDX_H]), "%9.5f", 1, 0, 0, 0, "Injecttransit_h_%d",l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injecttransit->paraminject[INJECTTR_IDX_K]), "%9.5f", 1, 0, 0, 0, "Injecttransit_k_%d",l);
	    }
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injecttransit->paraminject[INJECTTR_IDX_MSTAR]), "%9.5f", 1, 0, 0, 0, "Injecttransit_Mstar_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injecttransit->paraminject[INJECTTR_IDX_RSTAR]), "%9.5f", 1, 0, 0, 0, "Injecttransit_Rstar_%d",l);
	  for(i=1;i<=c[l].Injecttransit->Nld;i++)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Injecttransit->paraminject[INJECTTR_IDX_LD+i-1]), "%9.5f", 1, 0, 0, 0, "Injecttransit_ld_%d_%d",i,l);
	    }
	  break;
	case CNUM_STARSPOT:
	  if(c[l].Starspot->pertype == PERTYPE_SPECIFIED)
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Starspot->period), "%14.8f", 2, 1, 0, 0, 0, "Starspot_Period_%d",l);
	  else
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Starspot->period), "%14.8f", 2, 0, 0, 0, 0, "Starspot_Period_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Starspot->a), "%9.5f", 1, 0, 0, 0, "Starspot_a_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Starspot->b), "%9.5f", 1, 0, 0, 0, "Starspot_b_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Starspot->alpha), "%9.5f", 1, 0, 0, 0, "Starspot_alpha_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Starspot->inclination), "%9.5f", 1, 0, 0, 0, "Starspot_inclination_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Starspot->chi), "%9.5f", 1, 0, 0, 0, "Starspot_chi_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Starspot->psi0), "%9.5f", 1, 0, 0, 0, "Starspot_psi0_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Starspot->mconst), "%9.5f", 1, 0, 0, 0, "Starspot_mconst_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Starspot->chisq), "%9.5f", 1, 0, 0, 0, "Starspot_chi2perdof_%d",l);
	  break;
	case CNUM_BLS:
	  for(j=1;j<=c[l].Bls->Npeak;j++)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->bper), "%14.8f", 2, 0, 0, 0, j-1, "BLS_Period_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->bt0), "%.17g", 2, 0, 0, 0, j-1, "BLS_Tc_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->snval), "%9.5f", 2, 0, 0, 0, j-1, "BLS_SN_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->bpow), "%9.5f", 2, 0, 0, 0, j-1, "BLS_SR_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->sde), "%9.5f", 2, 0, 0, 0, j-1, "BLS_SDE_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->depth), "%9.5f", 2, 0, 0, 0, j-1, "BLS_Depth_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->qtran), "%9.5f", 2, 0, 0, 0, j-1, "BLS_Qtran_%d_%d",j,l);
	      if(c[l].Bls->fittrap)
		{
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->qingress), "%9.5f", 2, 0, 0, 0, j-1, "BLS_Qingress_%d_%d",j,l);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->OOTmag), "%9.5f", 2, 0, 0, 0, j-1, "BLS_OOTmag_%d_%d",j,l);
		}
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->i1_ph), "%9.5f", 2, 0, 0, 0, j-1, "BLS_i1_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->i2_ph), "%9.5f", 2, 0, 0, 0, j-1, "BLS_i2_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->chisqrplus), "%9.5f", 2, 0, 0, 0, j-1, "BLS_deltaChi2_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->fraconenight), "%9.5f", 2, 0, 0, 0, j-1, "BLS_fraconenight_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].Bls->nt), "%5d", 2, 0, 0, 0, j-1, "BLS_Npointsintransit_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].Bls->Nt), "%5d", 2, 0, 0, 0, j-1, "BLS_Ntransits_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].Bls->Nbefore), "%5d", 2, 0, 0, 0, j-1, "BLS_Npointsbeforetransit_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].Bls->Nafter), "%5d", 2, 0, 0, 0, j-1, "BLS_Npointsaftertransit_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->rednoise), "%9.5f", 2, 0, 0, 0, j-1, "BLS_Rednoise_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->whitenoise), "%9.5f", 2, 0, 0, 0, j-1, "BLS_Whitenoise_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->sigtopink), "%9.5f", 2, 0, 0, 0, j-1, "BLS_SignaltoPinknoise_%d_%d",j,l);
	    }
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->bperpos), "%14.8f", 1, 0, 0, 0, "BLS_Period_invtransit_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->chisqrminus), "%9.5f", 1, 0, 0, 0, "BLS_deltaChi2_invtransit_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Bls->meanmagval), "%9.5f", 1, 0, 0, 0, "BLS_MeanMag_%d",l);
	  break;
	case CNUM_FIXPERBLS:
	  if(c[l].BlsFixPer->pertype != PERTYPE_SPECIFIED)
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->period), "%14.8f", 2, 0, 0, 0, 0, "BLSFixPer_Period_%d",l);
	  else
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->period), "%14.8f", 2, 1, 0, 0, 0, "BLSFixPer_Period_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->bt0), "%.17g", 1, 0, 0, 0, "BLSFixPer_Tc_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->bpow), "%9.5f", 1, 0, 0, 0, "BLSFixPer_SR_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->depth), "%9.5f", 1, 0, 0, 0, "BLSFixPer_Depth_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->qtran), "%9.5f", 1, 0, 0, 0, "BLSFixPer_Qtran_%d",l);
	  if(c[l].BlsFixPer->fittrap) {
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->qingress), "%9.5f", 1, 0, 0, 0, "BLSFixPer_Qingress_%d",l);
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->OOTmag), "%9.5f", 1, 0, 0, 0, "BLSFixPer_OOTmag_%d",l);
	  }
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->i1_ph), "%9.5f", 1, 0, 0, 0, "BLSFixPer_i1_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->i2_ph), "%9.5f", 1, 0, 0, 0, "BLSFixPer_i2_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->chisqrplus), "%9.5f", 1, 0, 0, 0, "BLSFixPer_deltaChi2_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->fraconenight), "%9.5f", 1, 0, 0, 0, "BLSFixPer_fraconenight_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].BlsFixPer->nt), "%5d", 1, 0, 0, 0, "BLSFixPer_Npointsintransit_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].BlsFixPer->Nt), "%5d", 1, 0, 0, 0, "BLSFixPer_Ntransits_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].BlsFixPer->Nbefore), "%5d", 1, 0, 0, 0, "BLSFixPer_Npointsbeforetransit_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].BlsFixPer->Nafter), "%5d", 1, 0, 0, 0, "BLSFixPer_Npointsaftertransit_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->rednoise), "%9.5f", 1, 0, 0, 0, "BLSFixPer_Rednoise_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->whitenoise), "%9.5f", 1, 0, 0, 0, "BLSFixPer_Whitenoise_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->sigtopink), "%9.5f", 1, 0, 0, 0, "BLSFixPer_SignaltoPinknoise_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->chisqrminus), "%9.5f", 1, 0, 0, 0, "BLSFixPer_deltaChi2_invtransit_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixPer->meanmagval), "%9.5f", 1, 0, 0, 0, "BLSFixPer_MeanMag_%d",l);
	  break;
	case CNUM_BLSFIXDURTC:
	  if(c[l].BlsFixDurTc->durtype != PERTYPE_SPECIFIED)
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->inputdur), "%9.5f", 1, 0, 0, 0, "BLSFixDurTc_Duration_%d",l);
	  else
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->inputdur), "%9.5f", 1, 1, 0, 0, "BLSFixDurTc_Duration_%d",l);	  
	  if(c[l].BlsFixDurTc->TCtype != PERTYPE_SPECIFIED)
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->inputTC), "%.17g", 1, 0, 0, 0, "BLSFixDurTc_Tc_%d",l);
	  else
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->inputTC), "%.17g", 1, 1, 0, 0, "BLSFixDurTc_Tc_%d",l);
	  if(c[l].BlsFixDurTc->fixdepth) {
	    if(c[l].BlsFixDurTc->depthtype != PERTYPE_SPECIFIED)
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->inputdepth), "%9.5f", 1, 0, 0, 0, "BLSFixDurTc_Depth_%d",l);
	    else
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->inputdepth), "%9.5f", 1, 1, 0, 0, "BLSFixDurTc_Depth_%d",l);
	    if(c[l].BlsFixDurTc->qgresstype != PERTYPE_SPECIFIED)
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->inputqgress), "%9.5f", 1, 0, 0, 0, "BLSFixDurTc_Qingress_%d",l);
	    else
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->inputqgress), "%9.5f", 1, 1, 0, 0, "BLSFixDurTc_Qingress_%d",l);
	  }
	  for(j=1;j<=c[l].BlsFixDurTc->Npeak;j++)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->bper), "%14.8f", 2, 0, 0, 0, j-1, "BLSFixDurTc_Period_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->snval), "%9.5f", 2, 0, 0, 0, j-1, "BLSFixDurTc_SN_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->bpow), "%9.5f", 2, 0, 0, 0, j-1, "BLSFixDurTc_SR_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->sde), "%9.5f", 2, 0, 0, 0, j-1, "BLSFixDurTc_SDE_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->depth), "%9.5f", 2, 0, 0, 0, j-1, "BLSFixDurTc_Depth_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->qtran), "%9.5f", 2, 0, 0, 0, j-1, "BLSFixDurTc_Qtran_%d_%d",j,l);
	      if(c[l].BlsFixDurTc->fittrap && !c[l].BlsFixDurTc->fixdepth)
		{
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->qingress), "%9.5f", 2, 0, 0, 0, j-1, "BLSFixDurTc_Qingress_%d_%d",j,l);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->OOTmag), "%9.5f", 2, 0, 0, 0, j-1, "BLSFixDurTc_OOTmag_%d_%d",j,l);
		}
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->chisqrplus), "%9.5f", 2, 0, 0, 0, j-1, "BLSFixDurTc_deltaChi2_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->fraconenight), "%9.5f", 2, 0, 0, 0, j-1, "BLSFixDurTc_fraconenight_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].BlsFixDurTc->nt), "%5d", 2, 0, 0, 0, j-1, "BLSFixDurTc_Npointsintransit_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].BlsFixDurTc->Nt), "%5d", 2, 0, 0, 0, j-1, "BLSFixDurTc_Ntransits_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].BlsFixDurTc->Nbefore), "%5d", 2, 0, 0, 0, j-1, "BLSFixDurTc_Npointsbeforetransit_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_INT, 0, &(c[l].BlsFixDurTc->Nafter), "%5d", 2, 0, 0, 0, j-1, "BLSFixDurTc_Npointsaftertransit_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->rednoise), "%9.5f", 2, 0, 0, 0, j-1, "BLSFixDurTc_Rednoise_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->whitenoise), "%9.5f", 2, 0, 0, 0, j-1, "BLSFixDurTc_Whitenoise_%d_%d",j,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->sigtopink), "%9.5f", 2, 0, 0, 0, j-1, "BLSFixDurTc_SignaltoPinknoise_%d_%d",j,l);
	    }
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->bperpos), "%14.8f", 1, 0, 0, 0, "BLSFixDurTc_Period_invtransit_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->chisqrminus), "%9.5f", 1, 0, 0, 0, "BLSFixDurTc_deltaChi2_invtransit_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].BlsFixDurTc->meanmagval), "%9.5f", 1, 0, 0, 0, "BLSFixDurTc_MeanMag_%d",l);
	  break;
	case CNUM_SOFTENEDTRANSIT:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].SoftenedTransit->period), "%14.8f", 1, 0, 0, 0, "SoftenedTransit_Period_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].SoftenedTransit->T0), "%14.8f", 1, 0, 0, 0, "SoftenedTransit_T0_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].SoftenedTransit->eta), "%14.8f", 1, 0, 0, 0, "SoftenedTransit_eta_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].SoftenedTransit->cval), "%14.8f", 1, 0, 0, 0, "SoftenedTransit_cval_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].SoftenedTransit->delta), "%14.8f", 1, 0, 0, 0, "SoftenedTransit_delta_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].SoftenedTransit->mconst), "%14.8f", 1, 0, 0, 0, "SoftenedTransit_mconst_%d",l);
	  if(c[l].SoftenedTransit->dokillharm)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].SoftenedTransit->per_harm_out), "%14.8f", 1, 0, 0, 0, "SoftenedTransit_perharm_%d",l);
	      for(i=2;i<=c[l].SoftenedTransit->nsubharm+1;i++)
		{
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].SoftenedTransit->subharmA), "%9.5f", 2, 0, 0, 0, i-2, "SoftenedTransit_Subharm_%d_Sincoeff_%d",i,l);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].SoftenedTransit->subharmB), "%9.5f", 2, 0, 0, 0, i-2, "SoftenedTransit_Subharm_%d_Coscoeff_%d",i,l);
		}
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].SoftenedTransit->fundA), "%9.5f", 1, 0, 0, 0, "SoftenedTransit_Fundamental_Sincoeff_%d",l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].SoftenedTransit->fundB), "%9.5f", 1, 0, 0, 0, "SoftenedTransit_Fundamental_Coscoeff_%d",l);
	      for(i=2;i<=c[l].SoftenedTransit->nharm+1;i++)
		{
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].SoftenedTransit->harmA), "%9.5f", 2, 0, 0, 0, i-2, "SoftenedTransit_Harm_%d_Sincoeff_%d",i,l);
		  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].SoftenedTransit->harmB), "%9.5f", 2, 0, 0, 0, i-2, "SoftenedTransit_Harm_%d_Coscoeff_%d",i,l);
		}
	    }
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].SoftenedTransit->chisq), "%14.8f", 1, 0, 0, 0, "SoftenedTransit_chi2perdof_%d",l);
	  break;
	case CNUM_MANDELAGOLTRANSIT:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MandelAgolTransit->period), "%14.8f", 1, 0, 0, 0, "MandelAgolTransit_Period_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MandelAgolTransit->T0), "%14.8f", 1, 0, 0, 0, "MandelAgolTransit_T0_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MandelAgolTransit->r), "%9.5f", 1, 0, 0, 0, "MandelAgolTransit_r_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MandelAgolTransit->a), "%9.5f", 1, 0, 0, 0, "MandelAgolTransit_a_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MandelAgolTransit->bimpact), "%9.5f", 1, 0, 0, 0, "MandelAgolTransit_bimpact_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MandelAgolTransit->inc), "%9.5f", 1, 0, 0, 0, "MandelAgolTransit_inc_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MandelAgolTransit->e), "%9.5f", 1, 0, 0, 0, "MandelAgolTransit_e_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MandelAgolTransit->omega), "%9.5f", 1, 0, 0, 0, "MandelAgolTransit_omega_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MandelAgolTransit->mconst), "%9.5f", 1, 0, 0, 0, "MandelAgolTransit_mconst_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MandelAgolTransit->ldcoeffs), "%9.5f", 2, 0, 0, 0, 0, "MandelAgolTransit_ldcoeff1_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MandelAgolTransit->ldcoeffs), "%9.5f", 2, 0, 0, 0, 1, "MandelAgolTransit_ldcoeff2_%d",l);
	  if(c[l].MandelAgolTransit->type)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MandelAgolTransit->ldcoeffs), "%9.5f", 2, 0, 0, 0, 2, "MandelAgolTransit_ldcoeff3_%d",l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MandelAgolTransit->ldcoeffs), "%9.5f", 2, 0, 0, 0, 3, "MandelAgolTransit_ldcoeff4_%d",l);
	    }
	  if(c[l].MandelAgolTransit->fitRV)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MandelAgolTransit->K), "%9.5f", 1, 0, 0, 0, "MandelAgolTransit_K_%d",l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MandelAgolTransit->gamma), "%9.5f", 1, 0, 0, 0, "MandelAgolTransit_gamma_%d",l);
	    }
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MandelAgolTransit->chisq), "%9.5f", 1, 0, 0, 0, "MandelAgolTransit_chi2_%d",l);
	  break;
	case CNUM_MICROLENS:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MicroLens->f0), "%.14g", 1, 0, 0, 0, "Microlens_f0_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MicroLens->f1), "%.14g", 1, 0, 0, 0, "Microlens_f1_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MicroLens->f0), "%.14g", 1, 0, 0, 0, "Microlens_u0_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MicroLens->t0), "%.14g", 1, 0, 0, 0, "Microlens_t0_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MicroLens->tmax), "%.14g", 1, 0, 0, 0, "Microlens_tmax_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].MicroLens->chi2_), "%.14g", 1, 0, 0, 0, "Microlens_chi2perdof_%d",l);
	  break;
	case CNUM_SYSREM:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Sysrem->mag_ave), "%9.5f", 1, 0, 0, 0, "SYSREM_MeanMag_%d",l);
	  for(i=0;i<c[l].Sysrem->Nsysrem_total;i++)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Sysrem->final_colors[i]), "%9.5f", 1, 0, 0, 0, "SYSREM_Trend_%d_Coeff_%d",i,l);
	    }
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Sysrem->rms_out), "%9.5f", 1, 0, 0, 0, "SYSREM_RMS_%d",l);
	  break;
	case CNUM_TFA:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].TFA->ave_out), "%9.5f", 1, 0, 0, 0, "TFA_MeanMag_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].TFA->rms_out), "%9.5f", 1, 0, 0, 0, "TFA_RMS_%d",l);
	  break;
	case CNUM_TFA_SR:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].TFA_SR->ave_out), "%9.5f", 1, 0, 0, 0, "TFA_SR_MeanMag_%d",l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].TFA_SR->rms_out), "%9.5f", 1, 0, 0, 0, "TFA_SR_RMS_%d",l);
	  break;
	case CNUM_WWZ:
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].WWZ->max_z), "%.17g", 1, 0, 0, 0, "MaxWWZ_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].WWZ->max_freq), "%.17g", 1, 0, 0, 0, "MaxWWZ_Freq_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].WWZ->max_tau), "%.17g", 1, 0, 0, 0, "MaxWWZ_TShift_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].WWZ->max_pow), "%.17g", 1, 0, 0, 0, "MaxWWZ_Power_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].WWZ->max_amp), "%.17g", 1, 0, 0, 0, "MaxWWZ_Amplitude_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].WWZ->max_neff), "%.17g", 1, 0, 0, 0, "MaxWWZ_Neffective_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].WWZ->max_con), "%.17g", 1, 0, 0, 0, "MaxWWZ_AverageMag_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].WWZ->med_z), "%.17g", 1, 0, 0, 0, "Med_WWZ_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].WWZ->med_freq), "%.17g", 1, 0, 0, 0, "Med_Freq_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].WWZ->med_pow), "%.17g", 1, 0, 0, 0, "Med_Power_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].WWZ->med_amp), "%.17g", 1, 0, 0, 0, "Med_Amplitude_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].WWZ->med_neff), "%.17g", 1, 0, 0, 0, "Med_Neffective_%d", l);
	  addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].WWZ->med_con), "%.17g", 1, 0, 0, 0, "Med_AverageMag_%d", l);
	  break;
	case CNUM_DFTCLEAN:
	  for(i=0;i<c[l].Dftclean->Npeaks_dirty;i++)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Dftclean->peakfreqs_dirty), "%14.8f", 2, 0, 0, 0, i, "DFTCLEAN_DSPEC_PEAK_FREQ_%d_%d",i,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Dftclean->peakpows_dirty), "%9g", 2, 0, 0, 0, i, "DFTCLEAN_DSPEC_PEAK_POW_%d_%d",i,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Dftclean->SNR_dirty), "%9g", 2, 0, 0, 0, i, "DFTCLEAN_DSPEC_PEAK_SNR_%d_%d",i,l);
	    }
	  if(c[l].Dftclean->verboseout && c[l].Dftclean->Npeaks_dirty > 0) {
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Dftclean->aveper_dirty), "%9g", 1, 0, 0, 0, "DFTCLEAN_DSPEC_AVESPEC_%d",l);
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Dftclean->stdper_dirty), "%9g", 1, 0, 0, 0, "DFTCLEAN_DSPEC_STDSPEC_%d",l);
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Dftclean->aveper_noclip_dirty), "%9g", 1, 0, 0, 0, "DFTCLEAN_DSPEC_AVESPEC_NOCLIP_%d", l);
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Dftclean->stdper_noclip_dirty), "%9g", 1, 0, 0, 0, "DFTCLEAN_DSPEC_STDSPEC_NOCLIP_%d", l);
          }
	  for(i=0;i<c[l].Dftclean->Npeaks_clean;i++)
	    {
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Dftclean->peakfreqs_clean), "%14.8f", 2, 0, 0, 0, i, "DFTCLEAN_CSPEC_PEAK_FREQ_%d_%d",i,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Dftclean->peakpows_clean), "%9g", 2, 0, 0, 0, i, "DFTCLEAN_CSPEC_PEAK_POW_%d_%d",i,l);
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Dftclean->SNR_clean), "%9g", 2, 0, 0, 0, i, "DFTCLEAN_CSPEC_PEAK_SNR_%d_%d",i,l);
	    }
	  if(c[l].Dftclean->verboseout && c[l].Dftclean->Npeaks_clean > 0) {
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Dftclean->aveper_clean), "%9g", 1, 0, 0, 0, "DFTCLEAN_CSPEC_AVESPEC_%d",l);
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Dftclean->stdper_clean), "%9g", 1, 0, 0, 0, "DFTCLEAN_CSPEC_STDSPEC_%d",l);
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Dftclean->aveper_noclip_clean), "%9g", 1, 0, 0, 0, "DFTCLEAN_CSPEC_AVESPEC_NOCLIP_%d", l);
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Dftclean->stdper_noclip_clean), "%9g", 1, 0, 0, 0, "DFTCLEAN_CSPEC_STDSPEC_NOCLIP_%d", l);
          }
	  break;
	case CNUM_STATS:
	  for(k=0,i=0; i < c[l].Stats->Nvar; i++) {
	    for(j=0,m=0; j < c[l].Stats->Nstats; j++, k++) {
	      sprintf(tmpstring,"STATS_%s_",c[l].Stats->varnames[i]);
	      switch(c[l].Stats->statstocalc[j]) {
	      case VARTOOLS_STATSTYPE_MEAN:
		sprintf(tmpstring,"%sMEAN_%%d", tmpstring);
		break;
	      case VARTOOLS_STATSTYPE_WEIGHTEDMEAN:
		sprintf(tmpstring,"%sWEIGHTEDMEAN_%%d",tmpstring);
		break;
	      case VARTOOLS_STATSTYPE_MEDIAN:
		sprintf(tmpstring,"%sMEDIAN_%%d",tmpstring);
		break;
	      case VARTOOLS_STATSTYPE_MEDIAN_WEIGHT:
		sprintf(tmpstring,"%sWEIGHTEDMEDIAN_%%d",tmpstring);
		break;
	      case VARTOOLS_STATSTYPE_STDDEV:
		sprintf(tmpstring,"%sSTDDEV_%%d",tmpstring);
		break;
	      case VARTOOLS_STATSTYPE_MEDDEV:
		sprintf(tmpstring,"%sMEDDEV_%%d",tmpstring);
		break;
	      case VARTOOLS_STATSTYPE_MEDMEDDEV:
		sprintf(tmpstring,"%sMEDMEDDEV_%%d",tmpstring);
		break;
	      case VARTOOLS_STATSTYPE_MAD:
		sprintf(tmpstring,"%sMAD_%%d",tmpstring);
		break;
	      case VARTOOLS_STATSTYPE_KURTOSIS:
		sprintf(tmpstring,"%sKURTOSIS_%%d",tmpstring);
		break;
	      case VARTOOLS_STATSTYPE_SKEWNESS:
		sprintf(tmpstring,"%sSKEWNESS_%%d",tmpstring);
		break;
	      case VARTOOLS_STATSTYPE_PERCENTILE:
		sprintf(tmpstring,"%sPCT%.2f_%%d",tmpstring,c[l].Stats->pctval[m]);
		m++;
		break;
	      case VARTOOLS_STATSTYPE_PERCENTILE_WEIGHT:
		sprintf(tmpstring,"%sWPCT%.2f_%%d",tmpstring,c[l].Stats->pctval[m]);
		m++;
		break;
	      case VARTOOLS_STATSTYPE_MAXIMUM:
		sprintf(tmpstring,"%sMAX_%%d",tmpstring);
		break;
	      case VARTOOLS_STATSTYPE_MINIMUM:
		sprintf(tmpstring,"%sMIN_%%d",tmpstring);
		break;
	      case VARTOOLS_STATSTYPE_SUM:
		sprintf(tmpstring,"%sSUM_%%d",tmpstring);
		break;
	      default:
		break;
	      }
	      addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].Stats->statsout), "%.17g", 2, 0, 0, 0, k, tmpstring, l);
	    }
	  }
	  break;
#ifdef _HAVE_PYTHON
	case CNUM_PYTHON:
	  for(i=0; i < c[l].PythonCommand->Noutcolumnvars; i++) {
	    sprintf(tmpstring,"PYTHON_%s_%%d",c[l].PythonCommand->outcolumnnames[i]);
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].PythonCommand->outcolumndata), "%.17g", 2, 0, 0, 0, i, tmpstring, l);
	  }
	  break;
#endif
#ifdef _HAVE_R
	case CNUM_R:
	  for(i=0; i < c[l].RCommand->Noutcolumnvars; i++) {
	    sprintf(tmpstring,"R_%s_%%d",c[l].RCommand->outcolumnnames[i]);
	    addcolumn(p, l, VARTOOLS_TYPE_DOUBLE, 0, &(c[l].RCommand->outcolumndata), "%.17g", 2, 0, 0, 0, i, tmpstring, l);
	  }
	  break;
#endif
	}
    }
}

void printheader_new(ProgramData *p, FILE *outfile)
{
  int i, j, k, l;
  if(!p->tabflag)
    fprintf(outfile,"#");
  for(i=0;i<p->Ncolumns;i++)
    {
      if(i)
	dotab(outfile,p->tabflag);
      fprintf(outfile,p->outcolumns[i].columnname);
    }
  fprintf(outfile,"\n");
  if(p->tabflag)
    {
      for(i=0;i<p->Ncolumns;i++)
	{
	  if(i)
	    fprintf(outfile,"\t");
	  for(j=0;j<p->outcolumns[i].length_columnname;j++)
	    fprintf(outfile,"-");
	}
      fprintf(outfile,"\n");
    }
}

void printresults_new(ProgramData *p, int lc, int reallc, FILE *outfile)
{
  int i, j, outi, maxlength;
  double outd;
  float outf;
  char outc;
  short outshort;
  long outlong;
  static char *outs;
  static int sizeouts = 0;

  if(!p->oneline) {
    for(i=0;i<p->Ncolumns;i++)
      {
	if(i)
	  dotab(outfile,p->tabflag);
	switch(p->outcolumns[i].type)
	  {
	  case VARTOOLS_TYPE_DOUBLE:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outd);
	    fprintf(outfile,p->outcolumns[i].outputformat,outd);
	    break;
	  case VARTOOLS_TYPE_INT:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outi);
	    fprintf(outfile,p->outcolumns[i].outputformat,outi);
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outf);
	    fprintf(outfile,p->outcolumns[i].outputformat,outf);
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outshort);
	    fprintf(outfile,p->outcolumns[i].outputformat,outshort);
	    break;
	  case VARTOOLS_TYPE_LONG:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outlong);
	    fprintf(outfile,p->outcolumns[i].outputformat,outlong);
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outc);
	    fprintf(outfile,p->outcolumns[i].outputformat,outc);
	    break;
	  case VARTOOLS_TYPE_STRING:
	    if(p->outcolumns[i].stringsize > sizeouts)
	      {
		if(!sizeouts)
		  {
		    if((outs = (char *) malloc(p->outcolumns[i].stringsize)) == NULL)
		      error(ERR_MEMALLOC);
		  }
		else
		  {
		    if((outs = (char *) realloc(outs, p->outcolumns[i].stringsize)) == NULL)
		      error(ERR_MEMALLOC);
		  }
		sizeouts = p->outcolumns[i].stringsize;
	      }
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, outs, p->outcolumns[i].stringsize);
	    fprintf(outfile,p->outcolumns[i].outputformat,outs);
	    break;
	  }
      }
  }
  else {
    maxlength = 0;
    for(i=0;i<p->Ncolumns;i++)
      {
	if(p->outcolumns[i].length_columnname > maxlength)
	  maxlength = p->outcolumns[i].length_columnname;
      }
    for(i=0;i<p->Ncolumns;i++)
      {
	fprintf(outfile,p->outcolumns[i].columnname);
	for(j=p->outcolumns[i].length_columnname; j < maxlength; j++)
	  fprintf(outfile," ");
	fprintf(outfile," = ");
	switch(p->outcolumns[i].type)
	  {
	  case VARTOOLS_TYPE_DOUBLE:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outd);
	    fprintf(outfile,p->outcolumns[i].outputformat,outd);
	    break;
	  case VARTOOLS_TYPE_INT:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outi);
	    fprintf(outfile,p->outcolumns[i].outputformat,outi);
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outf);
	    fprintf(outfile,p->outcolumns[i].outputformat,outf);
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outshort);
	    fprintf(outfile,p->outcolumns[i].outputformat,outshort);
	    break;
	  case VARTOOLS_TYPE_LONG:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outlong);
	    fprintf(outfile,p->outcolumns[i].outputformat,outlong);
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outc);
	    fprintf(outfile,p->outcolumns[i].outputformat,outc);
	    break;
	  case VARTOOLS_TYPE_STRING:
	    if(p->outcolumns[i].stringsize > sizeouts)
	      {
		if(!sizeouts)
		  {
		    if((outs = (char *) malloc(p->outcolumns[i].stringsize)) == NULL)
		      error(ERR_MEMALLOC);
		  }
		else
		  {
		    if((outs = (char *) realloc(outs, p->outcolumns[i].stringsize)) == NULL)
		      error(ERR_MEMALLOC);
		  }
		sizeouts = p->outcolumns[i].stringsize;
	      }
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, outs, p->outcolumns[i].stringsize);
	    fprintf(outfile,p->outcolumns[i].outputformat,outs);
	    break;
	  }
	fprintf(outfile,"\n");
      }
  }
  fprintf(outfile,"\n");
}

void growbuffer(char **buf, int *sizebuf, int newlen)
{
  if(newlen + 1 >= (*sizebuf)) {
    (*sizebuf) = 2*(*sizebuf);
    if(((*buf) = (char *) realloc(*buf, (*sizebuf)*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
  }
}

void printresults_buffer_new(ProgramData *p, int lc, int reallc, char **buf, int *sizebuf)
{
  int i, j, outi, maxlength;
  double outd;
  float outf;
  char outc;
  short outshort;
  long outlong;
  char *outs;
  int sizeouts = 0;
  int len;

  if((*sizebuf) == 0) {
    *sizebuf = 256;
    if((*buf = (char *) malloc((*sizebuf)*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
  }
  int indx=0;
  (*buf)[indx] = '\0';

  if(!p->oneline) {
    for(i=0;i<p->Ncolumns;i++)
      {
	if(i) {
	  growbuffer(buf, sizebuf, indx+1);
	  dotab_buffer(&((*buf)[indx]),p->tabflag);
	  indx++;
	}
	switch(p->outcolumns[i].type)
	  {
	  case VARTOOLS_TYPE_DOUBLE:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outd);
	    if((len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outd)) >= ((*sizebuf)-indx)) {
	      growbuffer(buf, sizebuf, indx+len+1);
	      len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outd);
	    }
	    indx += len;
	    break;
	  case VARTOOLS_TYPE_INT:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outi);
	    if((len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outi)) >= ((*sizebuf)-indx)) {
	      growbuffer(buf, sizebuf, indx+len+1);
	      len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outi);
	    }
	    indx += len;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outf);
	    if((len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outf)) >= ((*sizebuf)-indx)) {
	      growbuffer(buf, sizebuf, indx+len+1);
	      len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outf);
	    }
	    indx += len;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outshort);
	    if((len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outshort)) >= ((*sizebuf)-indx)) {
	      growbuffer(buf, sizebuf, indx+len+1);
	      len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outshort);
	    }
	    indx += len;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outlong);
	    if((len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outlong)) >= ((*sizebuf)-indx)) {
	      growbuffer(buf, sizebuf, indx+len+1);
	      len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outlong);
	    }
	    indx += len;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outc);
	    if((len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outc)) >= ((*sizebuf)-indx)) {
	      growbuffer(buf, sizebuf, indx+len+1);
	      len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outc);
	    }
	    indx += len;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    if(p->outcolumns[i].stringsize > sizeouts)
	      {
		if(!sizeouts)
		  {
		    if((outs = (char *) malloc(p->outcolumns[i].stringsize)) == NULL)
		      error(ERR_MEMALLOC);
		  }
		else
		  {
		    if((outs = (char *) realloc(outs, p->outcolumns[i].stringsize)) == NULL)
		      error(ERR_MEMALLOC);
		  }
		sizeouts = p->outcolumns[i].stringsize;
	      }
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, outs, p->outcolumns[i].stringsize);
	    if((len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outs)) >= ((*sizebuf)-indx)) {
	      growbuffer(buf, sizebuf, indx+len+1);
	      len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outs);
	    }
	    indx += len;
	    break;
	  }
      }
  }
  else {
    maxlength = 0;
    for(i=0;i<p->Ncolumns;i++)
      {
	if(p->outcolumns[i].length_columnname > maxlength)
	  maxlength = p->outcolumns[i].length_columnname;
      }
    for(i=0;i<p->Ncolumns;i++)
      {
	if((len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].columnname)) >= ((*sizebuf)-indx)) {
	  growbuffer(buf, sizebuf, indx+len+1);
	  len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].columnname);
	}
	indx += len;
	for(j=p->outcolumns[i].length_columnname; j < maxlength; j++) {
	  growbuffer(buf, sizebuf, indx+1);
	  sprintf(&((*buf)[indx])," ");
	  indx++;
	}
	growbuffer(buf, sizebuf, indx+3);
	sprintf(&((*buf)[indx])," = ");
	indx = indx+3;
	switch(p->outcolumns[i].type)
	  {
	  case VARTOOLS_TYPE_DOUBLE:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outd);
	    if((len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outd)) >= ((*sizebuf)-indx)) {
	      growbuffer(buf, sizebuf, indx+len+1);
	      len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outd);
	    }
	    indx += len;
	    break;
	  case VARTOOLS_TYPE_INT:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outi);
	    if((len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outi)) >= ((*sizebuf)-indx)) {
	      growbuffer(buf, sizebuf, indx+len+1);
	      len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outi);
	    }
	    indx += len;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outf);
	    if((len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outf)) >= ((*sizebuf)-indx)) {
	      growbuffer(buf, sizebuf, indx+len+1);
	      len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outf);
	    }
	    indx += len;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outshort);
	    if((len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outshort)) >= ((*sizebuf)-indx)) {
	      growbuffer(buf, sizebuf, indx+len+1);
	      len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outshort);
	    }
	    indx += len;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outlong);
	    if((len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outlong)) >= ((*sizebuf)-indx)) {
	      growbuffer(buf, sizebuf, indx+len+1);
	      len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outlong);
	    }
	    indx += len;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, &outc);
	    if((len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outc)) >= ((*sizebuf)-indx)) {
	      growbuffer(buf, sizebuf, indx+len+1);
	      len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outc);
	    }
	    indx += len;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    if(p->outcolumns[i].stringsize > sizeouts)
	      {
		if(!sizeouts)
		  {
		    if((outs = (char *) malloc(p->outcolumns[i].stringsize)) == NULL)
		      error(ERR_MEMALLOC);
		  }
		else
		  {
		    if((outs = (char *) realloc(outs, p->outcolumns[i].stringsize)) == NULL)
		      error(ERR_MEMALLOC);
		  }
		sizeouts = p->outcolumns[i].stringsize;
	      }
	    getoutcolumnvalue(&(p->outcolumns[i]), lc, reallc, p->outcolumns[i].type, outs, p->outcolumns[i].stringsize);
	    if((len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outs)) >= ((*sizebuf)-indx)) {
	      growbuffer(buf, sizebuf, indx+len+1);
	      len = snprintf(&((*buf)[indx]),(*sizebuf)-indx,p->outcolumns[i].outputformat,outs);
	    }
	    indx += len;
	    break;
	  }
	growbuffer(buf, sizebuf, indx+1);
	sprintf(&((*buf)[indx]),"\n");
	indx++;
      }
  }
  growbuffer(buf, sizebuf, indx+1);
  sprintf(&((*buf)[indx]),"\n");
  indx++;
}

void emptyresults_buffer(ProgramData *p, FILE *outfile) {
  _StringBuffer *buf;
  int i;
  for(i=0; i < p->Nbuffs_full; i++) {
    buf = p->full_buffer_stack[i];
    fprintf(outfile,buf->data);
    pushFreeBuffer(p,buf);
  }
  p->Nbuffs_full = 0;
}

_StringBuffer *popFreeBuffer(ProgramData *p) {
  int i;
  _StringBuffer *buff;
  if(!p->Nbuffs_free) return NULL;
  if(p->Nbuffs_free_stack > 0) {
    buff = p->free_buffer_stack[p->Nbuffs_free_stack-1];
    p->Nbuffs_free_stack--;
    p->Nbuffs_free--;
    return(buff);
  } else {
    buff = (_StringBuffer *) malloc(sizeof(_StringBuffer));
    buff->len = 0;
    p->Nbuffs_free--;
    return(buff);
  }
}

void pushFullBuffer(ProgramData *p, _StringBuffer *buf) {
  if(!p->full_buffer_stack_alloclen) {
    if((p->full_buffer_stack = (_StringBuffer **) malloc(sizeof(_StringBuffer *))) == NULL)
      error(ERR_MEMALLOC);
    p->full_buffer_stack_alloclen = 1;
  } else if(p->Nbuffs_full >= p->full_buffer_stack_alloclen-1) {
    if((p->full_buffer_stack = (_StringBuffer **) realloc(p->full_buffer_stack, (p->Nbuffs_full+1)*sizeof(_StringBuffer *))) == NULL)
      error(ERR_MEMALLOC);
    p->full_buffer_stack_alloclen = p->Nbuffs_full+1;
  }
  p->full_buffer_stack[p->Nbuffs_full] = buf;
  p->Nbuffs_full++;
}

void pushFreeBuffer(ProgramData *p, _StringBuffer *buf) {
  if(!p->free_buffer_stack_alloclen) {
    if((p->free_buffer_stack = (_StringBuffer **) malloc(sizeof(_StringBuffer *))) == NULL)
      error(ERR_MEMALLOC);
    p->free_buffer_stack_alloclen = 1;
  } else if(p->Nbuffs_free_stack >= p->free_buffer_stack_alloclen-1) {
    if((p->free_buffer_stack = (_StringBuffer **) realloc(p->free_buffer_stack, (p->Nbuffs_free_stack+1)*sizeof(_StringBuffer *))) == NULL)
      error(ERR_MEMALLOC);
    p->free_buffer_stack_alloclen = p->Nbuffs_free_stack+1;
  }
  p->free_buffer_stack[p->Nbuffs_free_stack] = buf;
  p->Nbuffs_free_stack++;
  p->Nbuffs_free++;
}
    
void InitializeOutputBufferStacks(ProgramData *p) {
  p->free_buffer_stack_alloclen = 0;
  p->full_buffer_stack_alloclen = 0;
  p->Nbuffs_free_stack = 0;
  p->Nbuffs_full = 0;
  p->free_buffer_stack = NULL;
  p->full_buffer_stack = NULL;
}
