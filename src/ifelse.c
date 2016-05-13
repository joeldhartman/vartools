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

_IfStack *CreateIfStack(void) {
  int j;
  _IfStack *ifs;
  if((ifs = (_IfStack *) malloc(sizeof(_IfStack))) == NULL)
    error(ERR_MEMALLOC);
  ifs->sizearray = 4;
  ifs->curpos = 0;
  if((ifs->IfPtrs = (_IfStruct **) malloc(ifs->sizearray * sizeof(_IfStruct *))) == NULL ||
     (ifs->istrue = (char *) malloc(ifs->sizearray * sizeof(char))) == NULL)
    error(ERR_MEMALLOC);
  for(j = 0; j < ifs->sizearray; j++) {
    ifs->IfPtrs[j] = NULL;
  }
  return(ifs);
}


void pushIfStack(_IfStack *stack, _IfStruct *IfStruct) {
  int j, iold;
  if(stack->curpos + 1 >= stack->sizearray) {
    iold = stack->sizearray;
    stack->sizearray *= 2;
    if((stack->IfPtrs = (_IfStruct **) realloc(stack->IfPtrs, stack->sizearray * sizeof(_IfStruct *))) == NULL ||
       (stack->istrue = (char *) realloc(stack->istrue, stack->sizearray*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
    for(j=iold; j < stack->sizearray; j++) 
      stack->IfPtrs[j] = NULL;
  }
  stack->IfPtrs[stack->curpos] = IfStruct;
  stack->curpos += 1;
}

_IfStruct * popIfStack(_IfStack *stack) {
  if(stack->curpos <= 0)
    return NULL;
  stack->curpos -= 1;
  return stack->IfPtrs[stack->curpos];
}

char TestIf(_IfStack *stack, ProgramData *p, Command *c, int lcindex, int threadindex)
/* Takes as input the stack of if structures for this command,
   pointers to the program data structure, and the command structure
   for the current command, and returns 1 if execution of the command
   is not blocked by an if statement, or 0 if it is blocked. */
{
  /* Take care of the case where the current command is not part of an
     if statement. */
  if(c->cnum != CNUM_IF) {
    if(!stack->curpos)
      return 1;
    else {
      if(!stack->istrue[stack->curpos - 1]) {
	if(c->cnum == CNUM_COPYLC) {
	  turnoffcopies_onecommand(p, c, threadindex, lcindex);
	}
      }
      return(stack->istrue[stack->curpos - 1]);
    }
  }

  switch(c->IfCommand->iftype) {
  case VARTOOLS_IFTYPE_IF:
    if(!stack->curpos ? 1 : stack->istrue[stack->curpos - 1]) {
      pushIfStack(stack, c->IfCommand->ifs);
      if(EvaluateExpression(lcindex, threadindex, 0, c->IfCommand->ifs->expressions[c->IfCommand->ifindex])) {
	c->IfCommand->ifs->wasfoundtrue[threadindex] = 1;
	stack->istrue[stack->curpos - 1] = 1;
	return 1;
      }
      else {
	c->IfCommand->ifs->wasfoundtrue[threadindex] = 0;
	stack->istrue[stack->curpos - 1] = 0;
	return 0;
      }
    } else {
      return 0;
    }
    break;
  case VARTOOLS_IFTYPE_ELIF:
  case VARTOOLS_IFTYPE_ELSE:
    if(c->IfCommand->ifs != stack->IfPtrs[stack->curpos - 1])
      return 0;
    if(c->IfCommand->ifs->wasfoundtrue[threadindex]) {
      stack->istrue[stack->curpos - 1] = 0;
      return 0;
    }
    else {
      if(EvaluateExpression(lcindex, threadindex, 0, c->IfCommand->ifs->expressions[c->IfCommand->ifindex])) {
	c->IfCommand->ifs->wasfoundtrue[threadindex] = 1;
	stack->istrue[stack->curpos - 1] = 1;
	return 1;
      }
      else {
	stack->istrue[stack->curpos - 1] = 0;
	return 0;
      }
    }
    break;
  case VARTOOLS_IFTYPE_FI:
    if(c->IfCommand->ifs != stack->IfPtrs[stack->curpos - 1])
      return 0;
    popIfStack(stack);
    if(stack->curpos <= 0)
      return 1;
    else
      return stack->istrue[stack->curpos - 1];
    break;
  default:
    error(ERR_CODEERROR);
  }
}

int ParseIfCommand(int *iret, int argc, char **argv, int cn, ProgramData *p, Command *c)
{
  int i, cn2, ifcounter;
  i = *iret;
  if((c[cn].IfCommand = (_IfCommand *) malloc(sizeof(_IfCommand))) == NULL)
    error(ERR_MEMALLOC);
  if(!strcmp(argv[i],"-if")) {
    p->isifcommands = 1;
    c[cn].IfCommand->iftype = VARTOOLS_IFTYPE_IF;
    c[cn].IfCommand->ifindex = 0;
    if((c[cn].IfCommand->ifs = (_IfStruct *) malloc(sizeof(_IfStruct))) == NULL)
      error(ERR_MEMALLOC);
    c[cn].IfCommand->ifs->sizearray = 4;
    if((c[cn].IfCommand->ifs->expressions = (_Expression **) malloc(c[cn].IfCommand->ifs->sizearray * sizeof(_Expression *))) == NULL)
      error(ERR_MEMALLOC);
    c[cn].IfCommand->ifs->nterms = 1;
    i++;
    if(i >= argc) {*iret = i; return 1;}
    if((c[cn].IfCommand->exprstring = (char *) malloc((strlen(argv[i])+1)*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
    sprintf(c[cn].IfCommand->exprstring,"%s",argv[i]);
  }
  else if(!strcmp(argv[i],"-elif")) {
    c[cn].IfCommand->iftype = VARTOOLS_IFTYPE_ELIF;
    ifcounter = 0;
    i++;
    if(i >= argc) {*iret = i; return 1;}
    if((c[cn].IfCommand->exprstring = (char *) malloc((strlen(argv[i])+1)*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
    sprintf(c[cn].IfCommand->exprstring,"%s",argv[i]);
    /* Find the if command corresponding to this object */
    for(cn2 = cn-1; cn2 >= 0; cn2--) {
      if(c[cn2].cnum == CNUM_IF) {
	if(c[cn2].IfCommand->iftype == VARTOOLS_IFTYPE_IF) {
	  if(!ifcounter) {
	    c[cn].IfCommand->ifindex = c[cn2].IfCommand->ifs->nterms;
	    c[cn].IfCommand->ifs = c[cn2].IfCommand->ifs;
	    if(c[cn].IfCommand->ifs->nterms >= c[cn].IfCommand->ifs->sizearray) {
	      c[cn].IfCommand->ifs->sizearray *= 2;
	      if((c[cn].IfCommand->ifs->expressions = (_Expression **) realloc(c[cn].IfCommand->ifs->expressions, c[cn].IfCommand->ifs->sizearray * sizeof(_Expression *))) == NULL)
		error(ERR_MEMALLOC);
	    }
	    c[cn].IfCommand->ifs->nterms += 1;
	    break;
	  }
	  else ifcounter--;
	}
	else if(c[cn2].IfCommand->iftype == VARTOOLS_IFTYPE_ELSE) {
	  if(!ifcounter) {
	    error(ERR_BADIFTHENELSE);
	  }
	}
	else if(c[cn2].IfCommand->iftype == VARTOOLS_IFTYPE_FI) ifcounter++;
      }
    }
    if(cn2 < 0)
      error(ERR_BADIFTHENELSE);
  }
  else if(!strcmp(argv[i],"-else")) {
    c[cn].IfCommand->iftype = VARTOOLS_IFTYPE_ELSE;
    ifcounter = 0;
    if((c[cn].IfCommand->exprstring = (char *) malloc(2*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
    sprintf(c[cn].IfCommand->exprstring,"1");
    /* Find the if command corresponding to this object */
    for(cn2 = cn-1; cn2 >= 0; cn2--) {
      if(c[cn2].cnum == CNUM_IF) {
	if(c[cn2].IfCommand->iftype == VARTOOLS_IFTYPE_IF) {
	  if(!ifcounter) {
	    c[cn].IfCommand->ifindex = c[cn2].IfCommand->ifs->nterms;
	    c[cn].IfCommand->ifs = c[cn2].IfCommand->ifs;
	    if(c[cn].IfCommand->ifs->nterms >= c[cn].IfCommand->ifs->sizearray) {
	      c[cn].IfCommand->ifs->sizearray *= 2;
	      if((c[cn].IfCommand->ifs->expressions = (_Expression **) realloc(c[cn].IfCommand->ifs->expressions, c[cn].IfCommand->ifs->sizearray * sizeof(_Expression *))) == NULL)
		error(ERR_MEMALLOC);
	    }
	    c[cn].IfCommand->ifs->nterms += 1;
	    break;
	  }
	  else ifcounter--;
	}
	else if(c[cn2].IfCommand->iftype == VARTOOLS_IFTYPE_ELSE) {
	  if(!ifcounter) {
	    error(ERR_BADIFTHENELSE);
	  }
	}
	else if(c[cn2].IfCommand->iftype == VARTOOLS_IFTYPE_FI) ifcounter++;
      }
    }
    if(cn2 < 0)
      error(ERR_BADIFTHENELSE);
  }
  else if(!strcmp(argv[i],"-fi")) {
    c[cn].IfCommand->iftype = VARTOOLS_IFTYPE_FI;
    ifcounter = 0;
    for(cn2 = cn-1; cn2 >= 0; cn2--) {
      if(c[cn2].cnum == CNUM_IF) {
	if(c[cn2].IfCommand->iftype == VARTOOLS_IFTYPE_IF) {
	  if(!ifcounter) {
	    c[cn].IfCommand->ifs = c[cn2].IfCommand->ifs;
	    break;
	  }
	  else ifcounter--;
	}
	else if(c[cn2].IfCommand->iftype == VARTOOLS_IFTYPE_FI) ifcounter++;
      }
    }
    if(cn2 < 0)
      error(ERR_BADIFTHENELSE);
  }
  *iret = i;
  return 0;
}

void dosaveifstackcopy(ProgramData *p, _CopyLC *c, int threadid) {
  _IfStruct *IfPtr;
  int j;
  if(c->IfStack[threadid] == NULL) {
    c->IfStack[threadid] = CreateIfStack();
  }
  while(c->IfStack[threadid]->curpos > 0)
    popIfStack(c->IfStack[threadid]);
  for(j=0; j < p->IfStack[threadid]->curpos; j++) {
    pushIfStack(c->IfStack[threadid],p->IfStack[threadid]->IfPtrs[j]);
    c->IfStack[threadid]->istrue[j] = p->IfStack[threadid]->istrue[j];
  }

  if(c->sizearray_IfStruct_wasfoundtrue_copy[threadid] < c->IfStack[threadid]->sizearray) {
    if(!c->sizearray_IfStruct_wasfoundtrue_copy[threadid]) {
      if((c->IfStruct_wasfoundtrue_copy[threadid] = (char *) malloc(c->IfStack[threadid]->sizearray*sizeof(char))) == NULL)
	error(ERR_MEMALLOC);
    } else {
      if((c->IfStruct_wasfoundtrue_copy[threadid] = (char *) realloc(c->IfStruct_wasfoundtrue_copy[threadid], c->IfStack[threadid]->sizearray*sizeof(char))) == NULL)
	c->sizearray_IfStruct_wasfoundtrue_copy[threadid] = c->IfStack[threadid]->sizearray;
	error(ERR_MEMALLOC);
    }
  }
  for(j=0; j < p->IfStack[threadid]->curpos; j++) {
    c->IfStruct_wasfoundtrue_copy[threadid][j] = p->IfStack[threadid]->IfPtrs[j]->wasfoundtrue[threadid];
  }

  return;
}


void dorestoreifstackcopy(ProgramData *p, _CopyLC *c, int sthreadid, int rthreadid) {
  _IfStruct *IfPtr;
  int j;
  while(p->IfStack[rthreadid]->curpos > 0)
    popIfStack(p->IfStack[rthreadid]);
  for(j=0; j < c->IfStack[sthreadid]->curpos; j++) {
    pushIfStack(p->IfStack[rthreadid],c->IfStack[sthreadid]->IfPtrs[j]);
    p->IfStack[rthreadid]->istrue[j] = c->IfStack[sthreadid]->istrue[j];
    p->IfStack[rthreadid]->IfPtrs[j]->wasfoundtrue[rthreadid] = 
      c->IfStruct_wasfoundtrue_copy[sthreadid][j];
  }

  return;
}
