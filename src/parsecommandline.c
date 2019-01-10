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

/* Routines for parsing the command line - these are part of the vartools program by J. Hartman */

void increaseNcommands(ProgramData *p, Command **c)
{
  p->Ncommands++;
  int i;
  if(p->Ncommands > p->sizecommandvector)
    {
      p->sizecommandvector += INITCOMMANDSIZE;
      *c = (Command *) realloc(*c,p->sizecommandvector * sizeof(Command));
      for(i=p->Ncommands; i < p->sizecommandvector; i++) {
	(*c)[i].require_sort = 0;
	(*c)[i].require_distinct = 0;
	(*c)[i].N_setparam_expr = 0;
      }
    }
}

void parse_setparam_expr(Command *c, char *exprstr, _Expression **exprptr)
{
  if(!c->N_setparam_expr) {
    if((c->setparam_EvalExprStrings = (char **) malloc(sizeof(char *))) == NULL ||
       (c->setparam_EvalExpressions = (_Expression ***) malloc(sizeof(_Expression **))) == NULL)
      error(ERR_MEMALLOC);
  } else {
    if((c->setparam_EvalExprStrings = (char **) realloc(c->setparam_EvalExprStrings, (c->N_setparam_expr + 1)*sizeof(char *))) == NULL ||
       (c->setparam_EvalExpressions = (_Expression ***) realloc(c->setparam_EvalExpressions, (c->N_setparam_expr + 1)*sizeof(_Expression **))) == NULL)
      error(ERR_MEMALLOC);
  }
  c->setparam_EvalExpressions[c->N_setparam_expr] = exprptr;
  if((c->setparam_EvalExprStrings[c->N_setparam_expr] = malloc(sizeof(exprstr)+1)) == NULL)
    error(ERR_MEMALLOC);
  sprintf(c->setparam_EvalExprStrings[c->N_setparam_expr],"%s",exprstr);
  c->N_setparam_expr += 1;
}

#define DEFAULTSIZEKILLHARMTERMS 100

void parsecommandline(int argc, char **argv, ProgramData *p, Command **cptr)
{
  int i, iterm, j, k, cn, l, m, Nlcs, Ncommands;

  FILE *inlist;

  char **inputlistlines, *teststring;

  int Nkillharmterms;
  int sizeperptrs;
  int *Npers;
  double ***perptrs;
  int sizeinputlistvec = 1000;

  Command *c;

  c = *cptr;

  p->Ncopycommands = 0;
  p->Ncopiestotal = 0;
  p->Ncopies = NULL;
  p->readformatused = 0;
  p->inputlcformatused = 0;
  p->inlistvars = 0;
  p->t = NULL;
  p->mag = NULL;
  p->sig = NULL;
  p->NJD = NULL;
  p->lcnames = NULL;
  p->stringid = NULL;
  p->stringid_idx = NULL;

  p->Nskipchar = 0;
  p->skipchars = NULL;

  cn = 0;
  for(i=1;i<argc;i++)
    {
      if(strlen(argv[i]) > 1 && (argv[i][0] == '-' && argv[i][1] == '-')) {
	sprintf(argv[i],"%s",&(argv[i][1]));
      }
      /* -i lcname : A single input light curve - allocate memory to store the names and the light curves */
      if(!strncmp(argv[i],"-i",2) && strlen(argv[i]) == 2)
	{
	  iterm = i;
	  if(p->fileflag)
	    help(argv[iterm],p);
	  p->fileflag = 1;
	  if(p->listflag)
	    help(argv[iterm],p);
	  p->Nlcs = 1;
	  if((p->lcnames = (char **) malloc(p->Nlcs * sizeof(char *))) == NULL ||
	     (p->NJD = (int *) malloc(p->Nlcs * sizeof(int))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<p->Nlcs;j++)
	    {
	      if((p->lcnames[j] = (char *) malloc(MAXLEN * sizeof(char))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"-",1) && strlen(argv[i]) == 1)
		{
		  sprintf(p->lcnames[0],"stdin");
		  p->readfromstdinflag = 1;
		}
	      else
		{
		  sprintf(p->lcnames[0],"%s",argv[i]);
		  p->readfromstdinflag = 0;
		}
	    }
	  else
	    help(argv[iterm],p);
#ifdef _USEBINARY_LC
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"binary")) {
	      p->binarylcinput = 1;
	    }
	    else
	      i--;
	  }
	  else
	    i--;
#endif	      
	    
	}

      /* -l lclist : A list of light curves, just store the list-name for now */
      else if(!strncmp(argv[i],"-l",2) && strlen(argv[i]) == 2)
	{
	  iterm = i;
	  p->listflag = 1;
	  if(p->fileflag)
	    help(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"-",1) && strlen(argv[i]) == 1)
		{
		  sprintf(p->lclistname,"stdin");
		  p->readfromstdinflag = 1;
		}
	      else
		{
		  sprintf(p->lclistname,"%s",argv[i]);
		  p->readfromstdinflag = 0;
		}
	    }
	  else
	    help(argv[iterm],p);
	  k = 1;
	  i++;
	  if(i < argc)
	    {
	      if(!strcmp(argv[i],"column"))
		{
		  i++;
		  if(i < argc)
		    {
		      k = atoi(argv[i]);
		    }
		  else
		    help(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  RegisterDataFromInputList(p, 
				    (void *)(&p->lcnames), 
				    VARTOOLS_TYPE_STRING,
				    0, -1, 0, 0, NULL,
				    k, "LC_Name");
				    
#ifdef _USEBINARY_LC
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"binary")) {
	      p->binarylcinput = 1;
	    }
	    else
	      i--;
	  }
	  else
	    i--;
#endif

	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"opencommand")) {
	      p->use_lc_open_exec_command = 1;
	      i++;
	      if(i < argc)
		{
		  if((p->lc_open_exec_command_str = (char *) malloc((strlen(argv[i])+1))) == NULL)
		    error(ERR_MEMALLOC);
		  sprintf(p->lc_open_exec_command_str,"%s",argv[i]);
		}
	      else
		help(argv[iterm],p);
	    }
	    else
	      i--;
	  }
	  else
	    i--;

	}

      /* -copylc Ncopies */
      else if(!strcmp(argv[i],"-copylc"))
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_COPYLC;
	  i++;
	  if(i >= argc)
	    listcommands(argv[iterm],p);
	  c[cn].CopyLC = CreateCopyLCCommand(p,argv[i],cn);
	  c[cn].Savelc = c[cn].CopyLC->s;
	  cn++;
	}

      /* -expr var\"=\"expression */
      else if(!strcmp(argv[i],"-expr"))
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_EXPRESSION;
	  i++;
	  if(i >= argc)
	    listcommands(argv[iterm],p);
	  c[cn].ExpressionCommand = CreateExpressionCommand(p,argv[i]);
	  cn++;
	}

      /* -findblends matchrad [\"radec\"] [\"xycol\" colx coly] <\"fix\" period | \"list\" [\"column\" col] | \"fixcolumn\" <colname | colnum>> [\"starlist\" starlistfile] [\"zeromag\" zeromagval] [\"nofluxconvert\"] [\"Nharm\" Nharm] [\"omatches\" outputmatchfile] */
      else if(!strncmp(argv[i],"-findblends",11) && strlen(argv[i]) == 11)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_FINDBLENDS;
	  if((c[cn].FindBlends = (_FindBlends *) malloc(sizeof(_FindBlends))) == NULL)
	    error(ERR_MEMALLOC);
	  /* Set the default values */
	  increaselinkedcols(p, &(c[cn].FindBlends->linkedcolumn_varname), "1", cn);
	  c[cn].FindBlends->zeromag = 25.0;
	  c[cn].FindBlends->sepstarlist = 0;
	  c[cn].FindBlends->converttoflux = 1;
	  c[cn].FindBlends->radec = 0;
	  c[cn].FindBlends->outputmatches = 0;
	  c[cn].FindBlends->Nharm = 2;
	  i++;
	  if(i < argc)
	    c[cn].FindBlends->matchrad = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"radec",5) && strlen(argv[i]) == 5)
		c[cn].FindBlends->radec = 1;
	      else
		i--;
	    }
	  else
	    i--;
	  k = 0;
	  l = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"xycol",5) && strlen(argv[i]) == 5) {
		i++;
		if(i < argc)
		  k = atoi(argv[i]);
		else
		  listcommands(argv[iterm],p);
		i++;
		if(i < argc)
		  l = atoi(argv[i]);
		else
		  listcommands(argv[iterm],p);
	      } else i--;
	    } else i--;
	  RegisterDataFromInputList(p,
				    (void *) (&(c[cn].FindBlends->varxyin)),
				    VARTOOLS_TYPE_DOUBLE,
				    2, cn, 1, 1, NULL, k, l,
				    "FINDBLENDS_STAR_XPOS",
				    "FINDBLENDS_STAR_YPOS");
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
		{
		  c[cn].FindBlends->pertype = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    {
		      c[cn].FindBlends->fixperiod = atof(argv[i]);
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"fixcolumn",9) && strlen(argv[i]) == 9)
		{
		  c[cn].FindBlends->pertype = PERTYPE_FIXCOLUMN;
		  i++;
		  if(i < argc)
		    increaselinkedcols(p, &(c[cn].FindBlends->linkedcolumn), argv[i], cn);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
		{
		  c[cn].FindBlends->pertype = PERTYPE_SPECIFIED;
		  k = 0;
		  i++;
		  if(i < argc) {
		    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
		      i++;
		      if(i < argc)
			k = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		    } else i--;
		  } else i--;
		  RegisterDataFromInputList(p,
					    (void *) (&(c[cn].FindBlends->periods)),
					    VARTOOLS_TYPE_DOUBLE,
					    1, cn, 0, 0, NULL, k,
					    "FINDBLENDS_PERIOD");
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"starlist",8) && strlen(argv[i]) == 8)
		{
		  c[cn].FindBlends->sepstarlist = 1;
		  i++;
		  if(i < argc)
		    {
		      strcpy(c[cn].FindBlends->starlistname,argv[i]);
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"zeromag",7) && strlen(argv[i]) == 7)
		{
		  i++;
		  if(i < argc)
		    c[cn].FindBlends->zeromag = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"nofluxconvert",13) && strlen(argv[i]) == 13)
		c[cn].FindBlends->converttoflux = 0;
	      else
		i--;
	    }
	  else
	    i--;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"Nharm",5) && strlen(argv[i]) == 5)
		{
		  i++;
		  if(i < argc)
		    c[cn].FindBlends->Nharm = atoi(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"omatches",8) && strlen(argv[i]) == 8)
		{
		  c[cn].FindBlends->outputmatches = 1;
		  i++;
		  if(i < argc)
		    strcpy(c[cn].FindBlends->outmatchesfilename,argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  if(! c[cn].FindBlends->sepstarlist)
	    p->readallflag = 1;
	  cn++;
	}

      /* -fluxtomag mag_constant offset : converting isis differential flux light curves into magnitudes */
      else if(!strncmp(argv[i],"-difffluxtomag",14) && strlen(argv[i]) == 14)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_DIFFFLUXTOMAG;
	  if((c[cn].DiffFluxtomag = (_DiffFluxtomag *) malloc(sizeof(_DiffFluxtomag))) == NULL)
	    error(ERR_MEMALLOC);
	  
	  i++;
	  if(i < argc)
	    c[cn].DiffFluxtomag->mag_constant1 = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  
	  i++;
	  if(i < argc)
	    c[cn].DiffFluxtomag->offset = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  k = 0;
	  i++;
	  if(i < argc) {
	    if(!strncmp(argv[i],"magcolumn",9) && strlen(argv[i]) == 9) {
	      i++;
	      if(i < argc)
		k = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
	    } else i--;
	  } else i--;
	  RegisterDataFromInputList(p,
				    (void *) (&(c[cn].DiffFluxtomag->magstar)),
				    VARTOOLS_TYPE_DOUBLE,
				    1, cn, 0, 0, NULL, k,
				    "DIFFFLUXTOMAG_MAGSTAR");
	  cn++;
	}

      else if(!strncmp(argv[i],"-fluxtomag",10) && strlen(argv[i]) == 10)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_FLUXTOMAG;
	  if((c[cn].Fluxtomag = (_Fluxtomag *) malloc(sizeof(_Fluxtomag))) == NULL)
	    error(ERR_MEMALLOC);
	  
	  i++;
	  if(i < argc)
	    c[cn].Fluxtomag->mag_constant1 = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  
	  i++;
	  if(i < argc)
	    c[cn].Fluxtomag->offset = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  cn++;
	}

      /* -clip clipfactor iter ["niter" num] ["median"] : light curve clipping */
      else if(!strncmp(argv[i],"-clip",5) && strlen(argv[i]) == 5)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_CLIP;
	  if((c[cn].Clip = (_Clip *) malloc(sizeof(_Clip))) == NULL)
	    error(ERR_MEMALLOC);

	  c[cn].Clip->niter = 0;
	  c[cn].Clip->usemedian = 0;

	  i++;
	  if(i < argc)
	    c[cn].Clip->sigclip = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  
	  i++;
	  if(i < argc)
	    {
	      c[cn].Clip->iter = atoi(argv[i]);
	      if(c[cn].Clip->iter != 0 && c[cn].Clip->iter != 1)
		error(ERR_WRONGITER);
	    }
	  else
	    listcommands(argv[iterm],p);
	  
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"niter")) {
	      i++;
	      if(i >= argc) listcommands(argv[iterm],p);
	      c[cn].Clip->niter = atoi(argv[i]);
	    } else
	      i--;
	  }
	  else
	    i--;
	  
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"median")) {
	      c[cn].Clip->usemedian = 1;
	    } else
	      i--;
	  }
	  else
	    i--;

	  cn++;
	}

      /* -rescalesig : simple sigma rescaling */
      else if(!strncmp(argv[i],"-rescalesig",11) && strlen(argv[i]) == 11)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_RESCALESIG;
	  if((c[cn].Rescalesig = (_Rescalesig *) malloc(sizeof(_Rescalesig))) == NULL)
	    error(ERR_MEMALLOC);
	  cn++;
	}

      /* -ensemblerescalesig erssigclip : ensemble sigma rescaling */
      else if(!strncmp(argv[i],"-ensemblerescalesig",19) && strlen(argv[i]) == 19)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  p->readallflag = 1;
	  c[cn].cnum = CNUM_ENSEMBLERESCALESIG;
	  if((c[cn].Ensemblerescalesig = (_Ensemblerescalesig *) malloc(sizeof(_Ensemblerescalesig))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    c[cn].Ensemblerescalesig->erssigclip = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  cn++;
	}

      /* -decorr Nglobalterms globalfile1 order1 ... Nlcterms lccolumn1 lcorder1 ... */
      else if(!strncmp(argv[i],"-decorr",7) && strlen(argv[i]) == 7)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_DECORR;
	  p->decorrflag = 1;
	  if((c[cn].Decorr = (_Decorr *) malloc(sizeof(_Decorr))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    c[cn].Decorr->correctlc = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Decorr->zeropointterm = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Decorr->subtractfirstterm = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Decorr->N_globalterms = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if((c[cn].Decorr->global_file_names = (char **) malloc(c[cn].Decorr->N_globalterms * sizeof(char *))) == NULL ||
	     (c[cn].Decorr->globalfile_order = (int *) malloc(c[cn].Decorr->N_globalterms * sizeof(int))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<c[cn].Decorr->N_globalterms;j++)
	    if((c[cn].Decorr->global_file_names[j] = (char *) malloc(MAXLEN * sizeof(char))) == NULL)
	      error(ERR_MEMALLOC);
	  for(j=0;j<c[cn].Decorr->N_globalterms;j++)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].Decorr->global_file_names[j],"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		c[cn].Decorr->globalfile_order[j] = atoi(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      if(c[cn].Decorr->globalfile_order[j] <= 0)
		error(ERR_WRONGORDER);
	    }
	  i++;
	  if(i < argc)
	    c[cn].Decorr->N_lcterms = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if((c[cn].Decorr->lc_order = (int *) malloc(c[cn].Decorr->N_lcterms * sizeof(int))) == NULL ||
	     (c[cn].Decorr->lc_columns = (int *) malloc(c[cn].Decorr->N_lcterms * sizeof(int))) == NULL)
	    error(ERR_MEMALLOC);
	  
	  if(c[cn].Decorr->N_lcterms > 0) {
	    if((c[cn].Decorr->lcdecorr_terms_in = (double ***) malloc(c[cn].Decorr->N_lcterms * sizeof(double **))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  
	  for(j=0;j<c[cn].Decorr->N_lcterms;j++)
	    {
	      i++;
	      if(i < argc)
		c[cn].Decorr->lc_columns[j] = atoi(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		c[cn].Decorr->lc_order[j] = atoi(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      if(c[cn].Decorr->lc_order[j] <= 0)
		error(ERR_WRONGORDER);
	    }
	  c[cn].Decorr->N_decorrterms = c[cn].Decorr->N_globalterms + c[cn].Decorr->N_lcterms;
	  if(c[cn].Decorr->zeropointterm)
	    c[cn].Decorr->N_decorrterms_total = 1;
	  else
	    c[cn].Decorr->N_decorrterms_total = 0;
	  k = 0;
	  if(c[cn].Decorr->N_decorrterms > 0)
	    if((c[cn].Decorr->order = (int *) malloc(c[cn].Decorr->N_decorrterms * sizeof(int))) == NULL)
	      error(ERR_MEMALLOC);
	  for(j=0;j<c[cn].Decorr->N_globalterms;j++)
	    {
	      c[cn].Decorr->order[k] = c[cn].Decorr->globalfile_order[j];
	      c[cn].Decorr->N_decorrterms_total += c[cn].Decorr->order[k];
	      k++;
	    }
	  for(j=0;j<c[cn].Decorr->N_lcterms;j++)
	    {
	      c[cn].Decorr->order[k] = c[cn].Decorr->lc_order[j];
	      c[cn].Decorr->N_decorrterms_total += c[cn].Decorr->order[k];
	      k++;
	    }
	  i++;
	  if(i < argc)
	    c[cn].Decorr->omodel = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].Decorr->omodel)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].Decorr->modeloutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].Decorr->modelsuffix,".decorr.model");
	    }
	  cn++;
	}

      /*-dftclean nbeam [\"maxfreq\" maxf] [\"outdspec\" dspec_outdir] [\"finddirtypeaks\" Npeaks [\"clip\" clip clipiter]] [\"outwfunc\" wfunc_outdir] [\"clean\" gain SNlimit [\"outcbeam\" cbeam_outdir] [\"outcspec\" cspec_outdir] [\"findcleanpeaks\" Npeaks [\"clip\" clip clipiter]]] [\"useampspec\"] [\"verboseout\"] */
      else if(!strncmp(argv[i],"-dftclean",9) && strlen(argv[i]) == 9)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].require_sort = 1;
	  c[cn].require_distinct = 1;
	  c[cn].cnum = CNUM_DFTCLEAN;
	  if((c[cn].Dftclean = (_Dftclean *) malloc(sizeof(_Dftclean))) == NULL)
	    error(ERR_MEMALLOC);

	  i++;
	  if(i < argc)
	    c[cn].Dftclean->nbeam = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  c[cn].Dftclean->maxfreq = -1.;
	  c[cn].Dftclean->outdspec = 0;
	  c[cn].Dftclean->finddirtypeaks = 0;
	  c[cn].Dftclean->outwspec = 0;
	  c[cn].Dftclean->runclean = 0;
	  c[cn].Dftclean->outcbeam = 0;
	  c[cn].Dftclean->outcspec = 0;
	  c[cn].Dftclean->findcleanpeaks = 0;
	  c[cn].Dftclean->useampspec = 0;
	  c[cn].Dftclean->verboseout = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"maxfreq",7) && strlen(argv[i]) == 7)
		{
		  i++;
		  if(i < argc)
		    c[cn].Dftclean->maxfreq = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"outdspec",8) && strlen(argv[i]) == 8)
		{
		  c[cn].Dftclean->outdspec = 1;
		  i++;
		  if(i < argc)
		    {
		      sprintf(c[cn].Dftclean->dirtyspec_outdir,"%s",argv[i]);
		      sprintf(c[cn].Dftclean->dirtyspec_suffix,".dftclean.dspec");
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"finddirtypeaks",14) && strlen(argv[i]) == 14)
		{
		  c[cn].Dftclean->finddirtypeaks = 1;
		  i++;
		  if(i < argc)
		    {
		      c[cn].Dftclean->Npeaks_dirty = atoi(argv[i]);
		    }
		  else
		    listcommands(argv[iterm],p);
		  c[cn].Dftclean->clip_dirty = 5.0;
		  c[cn].Dftclean->clipiter_dirty = 1;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"clip",4) && strlen(argv[i]) == 4)
			{
			  i++;
			  if(i < argc) 
			    c[cn].Dftclean->clip_dirty = atof(argv[i]);
			  else
			    listcommands(argv[iterm],p);
			  i++;
			  if(i < argc)
			    c[cn].Dftclean->clipiter_dirty = atoi(argv[i]);
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			i--;
		    }
		  else 
		    i--;
		}
	      else
		i--;
	    }
	  else
	    i--;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"outwfunc",8) && strlen(argv[i]) == 8)
		{
		  c[cn].Dftclean->outwspec = 1;
		  i++;
		  if(i < argc)
		    {
		      sprintf(c[cn].Dftclean->wspec_outdir,"%s",argv[i]);
		      sprintf(c[cn].Dftclean->wspec_suffix,".dftclean.wfunc");
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"clean",5) && strlen(argv[i]) == 5)
		{
		  c[cn].Dftclean->runclean = 1;
		  i++;
		  if(i < argc)
		    {
		      c[cn].Dftclean->gain = atof(argv[i]);
		    }
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Dftclean->SNlimit = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"outcbeam",8) && strlen(argv[i]) == 8)
			{
			  c[cn].Dftclean->outcbeam = 1;
			  i++;
			  if(i < argc)
			    {
			      sprintf(c[cn].Dftclean->cbeam_outdir,"%s",argv[i]);
			      sprintf(c[cn].Dftclean->cbeam_suffix,".dftclean.cbeam");
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			i--;
		    }
		  else
		    i--;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"outcspec",8) && strlen(argv[i]) == 8)
			{
			  c[cn].Dftclean->outcspec = 1;
			  i++;
			  if(i < argc)
			    {
			      sprintf(c[cn].Dftclean->cspec_outdir,"%s",argv[i]);
			      sprintf(c[cn].Dftclean->cspec_suffix,".dftclean.cspec");
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			i--;
		    }
		  else
		    i--;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"findcleanpeaks",14) && strlen(argv[i]) == 14)
			{
			  c[cn].Dftclean->findcleanpeaks = 1;
			  i++;
			  if(i < argc)
			    c[cn].Dftclean->Npeaks_clean = atoi(argv[i]);
			  else
			    listcommands(argv[iterm],p);
			  c[cn].Dftclean->clip_clean = 5.0;
			  c[cn].Dftclean->clipiter_clean = 1;
			  i++;
			  if(i < argc)
			    {
			      if(!strncmp(argv[i],"clip",4) && strlen(argv[i]) == 4)
				{
				  i++;
				  if(i < argc) 
				    c[cn].Dftclean->clip_clean = atof(argv[i]);
				  else
				    listcommands(argv[iterm],p);
				  i++;
				  if(i < argc)
				    c[cn].Dftclean->clipiter_clean = atoi(argv[i]);
				  else
				    listcommands(argv[iterm],p);
				}
			      else
				i--;
			    }
			  else 
			    i--;
			}
		      else
			i--;
		    }
		  else
		    i--;
		}
	      else
		i--;
	    }
	  else
	    i--;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"useampspec",10) && strlen(argv[i]) == 10)
		c[cn].Dftclean->useampspec = 1;
	      else
		i--;
	    }
	  else
	    i--;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"verboseout",10) && strlen(argv[i]) == 10)
		c[cn].Dftclean->verboseout = 1;
	      else
		i--;
	    }
	  else
	    i--;
		  
	  cn++;
	}
      
      /* -chi2 : Calculate un-binned chi2 */
      else if(!strncmp(argv[i],"-chi2",5) && strlen(argv[i]) == 5)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_CHI2_NOBIN;
	  if((c[cn].Chi2_NoBin = (_Chi2_NoBin *) malloc(sizeof(_Chi2_NoBin))) == NULL)
	    error(ERR_MEMALLOC);
	  cn++;
	}
      
      /* -chi2bin Nbin bin1 ... binn : Calculate chi2 in a moving mean filter */
      else if(!strncmp(argv[i],"-chi2bin",8) && strlen(argv[i]) == 8)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_CHI2_BIN;
	  if((c[cn].Chi2_Bin = (_Chi2_Bin *) malloc(sizeof(_Chi2_Bin))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    c[cn].Chi2_Bin->Nbin = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);

	  if((c[cn].Chi2_Bin->bintimes = (double *) malloc(c[cn].Chi2_Bin->Nbin * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<c[cn].Chi2_Bin->Nbin;j++)
	    {
	      i++;
	      if(i < argc)
		c[cn].Chi2_Bin->bintimes[j] = atof(argv[i]) / MINUTESPERDAY;
	      else
		listcommands(argv[iterm],p);
	    }
	  cn++;
	}

      /* -changeerror : Replace the formal errors in a light curve with the light curve RMS */
      else if(!strncmp(argv[i],"-changeerror",12) && strlen(argv[i]) == 12)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_CHANGEERROR;
	  if((c[cn].Changeerror = (_Changeerror *) malloc(sizeof(_Changeerror))) == NULL)
	    error(ERR_MEMALLOC);
	  cn++;
	}

      /* -changevariable <\"t\" | \"mag\" | \"err\" | \"id\"> var */
      else if(!strcmp(argv[i],"-changevariable"))
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_CHANGEVARIABLE;
	  if((c[cn].Changevariable = (_Changevariable *) malloc(sizeof(_Changevariable))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i >= argc)
	    listcommands(argv[iterm],p);
	  if(!strcmp(argv[i],"t")) {
	    c[cn].Changevariable->changevar = VARTOOLS_CHANGEVAR_TIME;
	  }
	  else if(!strcmp(argv[i],"mag")) {
	    c[cn].Changevariable->changevar = VARTOOLS_CHANGEVAR_MAG;
	  }
	  else if(!strcmp(argv[i],"err")) {
	    c[cn].Changevariable->changevar = VARTOOLS_CHANGEVAR_ERR;
	  }
	  else if(!strcmp(argv[i],"id")) {
	    c[cn].Changevariable->changevar = VARTOOLS_CHANGEVAR_ID;
	  }
	  else {
	    listcommands(argv[iterm],p);
	  }
	  i++;
	  if(i >= argc)
	    listcommands(argv[iterm],p);
	  sprintf(c[cn].Changevariable->newvarname,"%s",argv[i]);
	  cn++;
	}
      /* -restorelc savenumber */
      else if(!strncmp(argv[i],"-restorelc",10) && strlen(argv[i]) == 10)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_RESTORELC;
	  if((c[cn].Restorelc = (_Restorelc *) malloc(sizeof(_Restorelc))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    {
	      c[cn].Restorelc->savenumber = atoi(argv[i]);
	      /* Get the -savelc command corresponding to this savenumber */
	      m = 0;
	      /* Get the index for the last aov command */
	      for(l = 0; l < cn; l++)
		{
		  if(c[l].cnum == CNUM_SAVELC)
		    m++;
		  if(m == c[cn].Restorelc->savenumber)
		    break;
		}
	      if(m == c[cn].Restorelc->savenumber)
		{
		  c[cn].Restorelc->saveindex = l;
		}
	      else
		error2(ERR_MISSINGSAVELC,argv[i]);
	    }
	  else
	    listcommands(argv[iterm],p);	      
	  cn++;
	}

      /* -savelc */
      else if(!strncmp(argv[i],"-savelc",7) && strlen(argv[i]) == 7)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_SAVELC;
	  if((c[cn].Savelc = (_Savelc *) malloc(sizeof(_Savelc))) == NULL)
	    error(ERR_MEMALLOC);
	  cn++;
	}

      /* -harmonicfilter <"full" |
                          "highpass"
                             <"minfreq" <"fix" value | 
                                         "list" [\"column\" col] | 
                                         "fixcolumn" <colname | colnum> | 
                                         "expr" expression>> |
                          "lowpass"
                             <"maxfreq" <"fix" value | 
                                         "list" [\"column\" col] | 
                                         "fixcolumn" <colname | colnum> | 
                                         "expr" expression>> |
                          "bandpass"
                             <"minfreq" <"fix" value | 
                                         "list" [\"column\" col] | 
                                         "fixcolumn" <colname | colnum> | 
                                         "expr" expression>>
                             <"maxfreq" <"fix" value | 
                                         "list" [\"column\" col] | 
                                         "fixcolumn" <colname | colnum> | 
                                         "expr" expression>> |
                          "bandcut"
                             <"minfreq" <"fix" value | 
                                         "list" [\"column\" col] | 
                                         "fixcolumn" <colname | colnum> | 
                                         "expr" expression>>
                             <"maxfreq" <"fix" value | 
                                         "list" [\"column\" col] | 
                                         "fixcolumn" <colname | colnum> | 
                                         "expr" expression>> >
                        ["filterexpr" expr] ["fullspec"] ["forcefft"]
                        ["ofourier" outdir ["nameformat" format]] */

      else if(!strcmp(argv[i],"-harmonicfilter"))
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].require_sort = 1;
	  c[cn].require_distinct = 1;
	  c[cn].cnum = CNUM_HARMONICFILTER;
	  if((c[cn].HarmonicFilter = (_HarmonicFilter *) malloc(sizeof(_HarmonicFilter))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(ParseHarmonicFilterCommand(&i, argc, argv, p, c[cn].HarmonicFilter, cn))
	    listcommands(argv[iterm],p);
	  cn++;
	}

      /* -resample <"nearest" | 
                    "linear"  |
                    "spline"  ["left" yp1] ["right" ypn] |
                    "splinemonotonic" |
                    "bspline" ["nbreaks" nbreaks] ["order" order] >
                   ["file" <"fix" times_file ["column" time_column] |
                            "list" ["listcolumn" col] ["tcolumn" time_column]
                           > |
                    ["tstart" <"fix" tstart | "fixcolumn" <colname | colnum> |
                               "list" ["column" col] | "expr" expression > ]
                    ["tstop" <"fix" tstop | "fixcolumn" <colname | colnum> |
                               "list" ["column" col] | "expr" expression > ]
                    [["delt" <"fix" delt | "fixcolumn" <colname | colnum> |
                               "list" ["column" col] | "expr" expression > ]
                     | ["Npoints" <"fix" Np | "fixcolumn" <colname | colnum> |
                               "list" ["column" col] | "expr" expression > ]]] 
                    ["gaps" <"fix" value | "fixcolumn" <colname | colnum> |
                               "list" ["column" col] | "expr" expression |
                             "frac_min_sep" value | "frac_med_sep" value |
                             "percentile_sep" value >
                     <"nearest" | "linear" | "spline" ["left" yp1] ["right" ypn]
                      | "splinemonotonic" | "bspline" ["nbreaks" nbreaks]
                      ["order" order]>]
                    ["extrap" 
                     <"nearest" | "linear" | "spline" ["left" yp1] ["right" ypn]
                      | "splinemonotonic" | "bspline" ["nbreaks" nbreaks]
                      ["order" order]>]

      */
      else if(!strcmp(argv[i],"-resample"))
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].require_sort = 1;
	  c[cn].require_distinct = 1;
	  c[cn].cnum = CNUM_RESAMPLE;
	  if((c[cn].Resample = (_Resample *) malloc(sizeof(_Resample))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(ParseResampleCommand(&i, argc, argv, p, c[cn].Resample, cn))
	    listcommands(argv[iterm],p);
	  cn++;
	}

      /* -restricttimes ["exclude"] <"JDrange" minJD maxJD | 
                         "JDrangebylc"
                            <"fix" minJD | "list" ["column" col] | 
                             "fixcolumn" <colname | colnum> |
                             "expr" expression>
                            <"fix" maxJD | "list" ["column" col] | 
                             "fixcolumn" <colname | colnum> |
                             "expr" expression> |
                         "JDlist" JDfilename |
                         "imagelist" imagefilename |
                         "expr" eval_expression>
      */
      else if(!strcmp(argv[i],"-restricttimes"))
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_RESTRICTTIMES;
	  if((c[cn].RestrictTimes = (_RestrictTimes *) malloc(sizeof(_RestrictTimes))) == NULL)
	    error(ERR_MEMALLOC);
	  c[cn].RestrictTimes->exclude = 0;
	  c[cn].RestrictTimes->saveexcludedpoints = 0;
	  i++;
	  if(i >= argc)
	    listcommands(argv[iterm],p);
	  if(!strcmp(argv[i],"exclude")) {
	    c[cn].RestrictTimes->exclude = 1;
	    i++;
	    if(i >= argc)
	      listcommands(argv[iterm],p);
	  }
	  if(!strcmp(argv[i],"JDrange")) {
	    c[cn].RestrictTimes->restricttype = VARTOOLS_RESTRICTTIMES_JDRANGE;
	    i++;
	    if(i >= argc)
	      listcommands(argv[iterm],p);
	    c[cn].RestrictTimes->minJDtype = PERTYPE_FIX;
	    c[cn].RestrictTimes->minJDfixval = atof(argv[i]);
	    i++;
	    if(i >= argc)
	      listcommands(argv[iterm],p);
	    c[cn].RestrictTimes->maxJDtype = PERTYPE_FIX;
	    c[cn].RestrictTimes->maxJDfixval = atof(argv[i]);
	  }
	  else if(!strcmp(argv[i],"JDrangebylc")) {
	    c[cn].RestrictTimes->restricttype = VARTOOLS_RESTRICTTIMES_JDRANGE;
	    i++;
	    if(i >= argc)
	      listcommands(argv[iterm],p);
	    if(!strcmp(argv[i],"fix"))
	      {
		c[cn].RestrictTimes->minJDtype = PERTYPE_FIX;
		i++;
		if(i < argc)
		  {
		    c[cn].RestrictTimes->minJDfixval = atof(argv[i]);
		  }
		else
		  listcommands(argv[iterm],p);
	      }
	    else if(!strcmp(argv[i],"fixcolumn"))
	      {
		c[cn].RestrictTimes->minJDtype = PERTYPE_FIXCOLUMN;
		i++;
		if(i < argc)
		  increaselinkedcols(p, &(c[cn].RestrictTimes->minJD_linkedcolumn), argv[i], cn);
		else
		  listcommands(argv[iterm],p);
	      }
	    else if(!strcmp(argv[i],"list"))
	      {
		c[cn].RestrictTimes->minJDtype = PERTYPE_SPECIFIED;
		k = 0;
		i++;
		if(i < argc) {
		  if(!strcmp(argv[i],"column")) {
		    i++;
		    if(i < argc)
		      k = atoi(argv[i]);
		    else
		      listcommands(argv[iterm],p);
		  } else i--;
		} else i--;
		RegisterDataFromInputList(p,
					  (void *) (&(c[cn].RestrictTimes->minJD)),
					  VARTOOLS_TYPE_DOUBLE,
					  0, cn, 0, 0, NULL, k,
					  "RESTRICTTIMES_MINJD");
	      }
	    else if(!strcmp(argv[i],"expr"))
	      {
		c[cn].RestrictTimes->minJDtype = PERTYPE_EXPR;
		i++;
		if(i >= argc)
		  listcommands(argv[iterm],p);
		RestrictTimes_ParseExpr(&i, argc, argv, p, c[cn].RestrictTimes, 0);
	      }
	    else
	      listcommands(argv[iterm],p);
	    i++;
	    if(i >= argc)
	      listcommands(argv[iterm],p);
	    if(!strcmp(argv[i],"fix"))
	      {
		c[cn].RestrictTimes->maxJDtype = PERTYPE_FIX;
		i++;
		if(i < argc)
		  {
		    c[cn].RestrictTimes->maxJDfixval = atof(argv[i]);
		  }
		else
		  listcommands(argv[iterm],p);
	      }
	    else if(!strcmp(argv[i],"fixcolumn"))
	      {
		c[cn].RestrictTimes->maxJDtype = PERTYPE_FIXCOLUMN;
		i++;
		if(i < argc)
		  increaselinkedcols(p, &(c[cn].RestrictTimes->maxJD_linkedcolumn), argv[i], cn);
		else
		  listcommands(argv[iterm],p);
	      }
	    else if(!strcmp(argv[i],"list"))
	      {
		c[cn].RestrictTimes->maxJDtype = PERTYPE_SPECIFIED;
		k = 0;
		i++;
		if(i < argc) {
		  if(!strcmp(argv[i],"column")) {
		    i++;
		    if(i < argc)
		      k = atoi(argv[i]);
		    else
		      listcommands(argv[iterm],p);
		  } else i--;
		} else i--;
		RegisterDataFromInputList(p,
					  (void *) (&(c[cn].RestrictTimes->maxJD)),
					  VARTOOLS_TYPE_DOUBLE,
					  0, cn, 0, 0, NULL, k,
					  "RESTRICTTIMES_MAXJD");
	      }
	    else if(!strcmp(argv[i],"expr"))
	      {
		c[cn].RestrictTimes->maxJDtype = PERTYPE_EXPR;
		i++;
		if(i >= argc)
		  listcommands(argv[iterm],p);
		RestrictTimes_ParseExpr(&i, argc, argv, p, c[cn].RestrictTimes, 1);
	      }
	    else
	      listcommands(argv[iterm],p);
	  }
	  else if(!strcmp(argv[i],"JDlist")) {
	    c[cn].RestrictTimes->restricttype = VARTOOLS_RESTRICTTIMES_JDLIST;
	    i++;
	    if(i >= argc)
	      listcommands(argv[iterm],p);
	    RestrictTimes_readJDlist(argv[i],&(c[cn].RestrictTimes->JD_restrictlist),&(c[cn].RestrictTimes->N_restrictlist));
	  }
	  else if(!strcmp(argv[i],"imagelist")) {
	    c[cn].RestrictTimes->restricttype = VARTOOLS_RESTRICTTIMES_IMAGELIST;
	    i++;
	    if(i >= argc)
	      listcommands(argv[iterm],p);
	    RestrictTimes_readimagelist(argv[i],&(c[cn].RestrictTimes->image_restrictlist),&(c[cn].RestrictTimes->image_restrictlist_indx),&(c[cn].RestrictTimes->N_restrictlist));
	    p->requirestringid = 1;
	  }
	  else if(!strcmp(argv[i],"expr")) {
	    c[cn].RestrictTimes->restricttype = VARTOOLS_RESTRICTTIMES_EXPR;
	    i++;
	    if(i >= argc)
	      listcommands(argv[iterm],p);
	    if((c[cn].RestrictTimes->restrictexprstring = (char *) malloc((strlen(argv[i]+1))*sizeof(char))) == NULL)
	      error(ERR_MEMALLOC);
	    sprintf(c[cn].RestrictTimes->restrictexprstring,"%s",argv[i]);
	  }
	  else
	    listcommands(argv[iterm],p);
	  cn++;
	}

      /* -restoretimes prior_restricttimes_command */
      else if(!strcmp(argv[i],"-restoretimes")) {
	iterm = i;
	increaseNcommands(p,&c);
	c[cn].cnum = CNUM_RESTORETIMES;
	  if((c[cn].RestoreTimes = (_RestoreTimes *) malloc(sizeof(_RestoreTimes))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    {
	      c[cn].RestoreTimes->restrictnum = atoi(argv[i]);
	      /* Get the -savelc command corresponding to this savenumber */
	      m = 0;
	      /* Get the index for the last -restricttimes command */
	      for(l = 0; l < cn; l++)
		{
		  if(c[l].cnum == CNUM_RESTRICTTIMES)
		    m++;
		  if(m == c[cn].RestoreTimes->restrictnum) {
		    c[cn].RestoreTimes->RestrictTimes = c[l].RestrictTimes;
		    c[l].RestrictTimes->saveexcludedpoints = 1;
		    break;
		  }
		}
	      if(l >= cn)
		error2(ERR_MISSINGRESTRICTTIMES,argv[i]);
	    }
	  else
	    listcommands(argv[iterm],p);	      
	  cn++;
      }
      
      /* -rms : Calculate un-binned rms */
      else if(!strncmp(argv[i],"-rms",4) && strlen(argv[i]) == 4)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_RMS_NOBIN;
	  if((c[cn].RMS_NoBin = (_RMS_NoBin *) malloc(sizeof(_RMS_NoBin))) == NULL)
	    error(ERR_MEMALLOC);
	  cn++;
	}

      /* -rmsbin Nbin bintime1 ... bintimen : Calculate rms after applying a moving mean filter */
      else if(!strncmp(argv[i],"-rmsbin",7) && strlen(argv[i]) == 7)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_RMS_BIN;
	  if((c[cn].RMS_Bin = (_RMS_Bin *) malloc(sizeof(_RMS_Bin))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    c[cn].RMS_Bin->Nbin = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if((c[cn].RMS_Bin->bintimes = (double *) malloc(c[cn].RMS_Bin->Nbin * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	  for(j=0;j<c[cn].RMS_Bin->Nbin;j++)
	    {
	      i++;
	      if(i < argc)
		c[cn].RMS_Bin->bintimes[j] = atof(argv[i]) / MINUTESPERDAY;
	      else
		listcommands(argv[iterm],p);
	    }
	  cn++;
	}

      /* -Jstet Jstet_time dates */
      else if(!strncmp(argv[i],"-Jstet",6) && strlen(argv[i]) == 6)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_JSTET;
	  c[cn].require_sort = 1;
	  c[cn].require_distinct = 1;
	  if((c[cn].Jstet = (_Jstet *) malloc(sizeof(_Jstet))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    c[cn].Jstet->Jstet_time = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    sprintf(c[cn].Jstet->datesname,"%s",argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  cn++;
	}

      /* -addnoise <"gamma" <"fix" val | "list">> <"sig_red" <"fix" val | "list">> <"sig_white <"fix" val | "list">> */
      else if(!strncmp(argv[i],"-addnoise",9) && strlen(argv[i]) == 9)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_ADDNOISE;
	  if((c[cn].AddNoise = (_AddNoise *) malloc(sizeof(_AddNoise))) == NULL)
	    error(ERR_MEMALLOC);

	  i++;
	  if(i >= argc) listcommands(argv[iterm],p);
	  if(!strcmp(argv[i],"white")) {
	    c[cn].AddNoise->noise_type = VARTOOLS_ADDNOISE_WHITE;
	    i++;
	    if(i < argc)
	      {
		if(!strncmp(argv[i],"sig_white",9) && strlen(argv[i]) == 9)
		  {
		    i++;
		    if(i < argc)
		      {
			if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			  {
			    c[cn].AddNoise->sig_w_type = PERTYPE_FIX;
			    i++;
			    if(i < argc)
			      c[cn].AddNoise->sig_w_fix = atof(argv[i]);
			    else
			      listcommands(argv[iterm],p);
			  }
			else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			  {
			    c[cn].AddNoise->sig_w_type = PERTYPE_SPECIFIED;
			    k = 0;
			    i++;
			    if(i < argc) {
			      if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
				{
				  i++;
				  if(i < argc)
				    k = atoi(argv[i]);
				  else
				    listcommands(argv[iterm],p);
				}
			      else
				i--;
			    }
			    else
			      i--;
			    RegisterDataFromInputList(p, 
						      (void *)(&(c[cn].AddNoise->sig_w)),
						      VARTOOLS_TYPE_DOUBLE,
						      1, cn, 0, 0, NULL, k,
						      "ADDNOISE_SIGMA_WHITE");
			  }
			else
			  listcommands(argv[iterm],p);
		      }
		    else
		      listcommands(argv[iterm],p);
		  }
		else
		  listcommands(argv[iterm],p);
	      }
	    else
	      listcommands(argv[iterm],p);
	  }
	  else if(!strcmp(argv[i],"squareexp") ||
		  !strcmp(argv[i],"exp")) {
	    c[cn].require_distinct = 1;
	    c[cn].require_sort = 1;
	    if(!strcmp(argv[i],"squareexp")) {
	      c[cn].AddNoise->noise_type = VARTOOLS_ADDNOISE_COVAR_SQUAREDEXPONENTIAL;
	    }
	    else if(!strcmp(argv[i],"exp")) {
	      c[cn].AddNoise->noise_type = VARTOOLS_ADDNOISE_COVAR_EXPONENTIAL;
	    }
	    i++;
	    if(i < argc)
	      {
		if(!strcmp(argv[i],"rho"))
		  {
		    i++;
		    if(i < argc)
		      {
			if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			  {
			    c[cn].AddNoise->rho_r_type = PERTYPE_FIX;
			    i++;
			    if(i < argc)
			      c[cn].AddNoise->rho_r_fix = atof(argv[i]);
			    else
			      listcommands(argv[iterm],p);
			  }
			else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			  {
			    c[cn].AddNoise->rho_r_type = PERTYPE_SPECIFIED;
			    k = 0;
			    i++;
			    if(i < argc)
			      {
				if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
				  {
				    i++;
				    if(i < argc)
				      k = atoi(argv[i]);
				    else
				      listcommands(argv[iterm],p);
				  }
				else
				  i--;
			      }
			    else
			      i--;
			    RegisterDataFromInputList(p, 
						      (void *)(&(c[cn].AddNoise->rho_r)),
						      VARTOOLS_TYPE_DOUBLE,
						      1, cn, 0, 0, NULL, k,
						      "ADDNOISE_RHO");
			  }
			else
			  listcommands(argv[iterm],p);
		      }
		    else
		      listcommands(argv[iterm],p);
		  }
		else
		  listcommands(argv[iterm],p);
	      }
	    else
	      listcommands(argv[iterm],p);
	    i++;
	    if(i < argc)
	      {
		if(!strncmp(argv[i],"sig_red",7) && strlen(argv[i]) == 7)
		  {
		    i++;
		    if(i < argc)
		      {
			if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			  {
			    c[cn].AddNoise->sig_r_type = PERTYPE_FIX;
			    i++;
			    if(i < argc)
			      c[cn].AddNoise->sig_r_fix = atof(argv[i]);
			    else
			      listcommands(argv[iterm],p);
			  }
			else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			  {
			    c[cn].AddNoise->sig_r_type = PERTYPE_SPECIFIED;
			    k = 0;
			    i++;
			    if(i < argc) {
			      if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
				{
				  i++;
				  if(i < argc)
				    k = atoi(argv[i]);
				  else
				    listcommands(argv[iterm],p);
				}
			      else
				i--;
			    }
			    else
			      i--;
			    RegisterDataFromInputList(p, 
						      (void *)(&(c[cn].AddNoise->sig_r)),
						      VARTOOLS_TYPE_DOUBLE,
						      1, cn, 0, 0, NULL, k,
						      "ADDNOISE_SIGMA_RED");
			  }
			else
			  listcommands(argv[iterm],p);
		      }
		    else
		      listcommands(argv[iterm],p);
		  }
		else
		  listcommands(argv[iterm],p);
	      }
	    else
	      listcommands(argv[iterm],p);
	    i++;
	    if(i < argc)
	      {
		if(!strncmp(argv[i],"sig_white",9) && strlen(argv[i]) == 9)
		  {
		    i++;
		    if(i < argc)
		      {
			if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			  {
			    c[cn].AddNoise->sig_w_type = PERTYPE_FIX;
			    i++;
			    if(i < argc)
			      c[cn].AddNoise->sig_w_fix = atof(argv[i]);
			    else
			      listcommands(argv[iterm],p);
			  }
			else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			  {
			    c[cn].AddNoise->sig_w_type = PERTYPE_SPECIFIED;
			    k = 0;
			    i++;
			    if(i < argc) {
			      if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
				{
				  i++;
				  if(i < argc)
				    k = atoi(argv[i]);
				  else
				    listcommands(argv[iterm],p);
				}
			      else
				i--;
			    }
			    else
			      i--;
			    RegisterDataFromInputList(p, 
						      (void *)(&(c[cn].AddNoise->sig_w)),
						      VARTOOLS_TYPE_DOUBLE,
						      1, cn, 0, 0, NULL, k,
						      "ADDNOISE_SIGMA_WHITE");
			  }
			else
			  listcommands(argv[iterm],p);
		      }
		    else
		      listcommands(argv[iterm],p);
		  }
		else
		  listcommands(argv[iterm],p);
	      }
	    else
	      listcommands(argv[iterm],p);

	    i++;
	    if(i < argc)
	      {
		if(!strcmp(argv[i],"bintime"))
		  {
		    i++;
		    if(i < argc)
		      {
			if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			  {
			    c[cn].AddNoise->bintime_type = PERTYPE_FIX;
			    i++;
			    if(i < argc)
			      c[cn].AddNoise->bintime_fix = atof(argv[i]);
			    else
			      listcommands(argv[iterm],p);
			  }
			else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			  {
			    c[cn].AddNoise->bintime_type = PERTYPE_SPECIFIED;
			    k = 0;
			    i++;
			    if(i < argc) {
			      if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
				{
				  i++;
				  if(i < argc)
				    k = atoi(argv[i]);
				  else
				    listcommands(argv[iterm],p);
				}
			      else
				i--;
			    }
			    else
			      i--;
			    RegisterDataFromInputList(p, 
						      (void *)(&(c[cn].AddNoise->bintime)),
						      VARTOOLS_TYPE_DOUBLE,
						      1, cn, 0, 0, NULL, k,
						      "ADDNOISE_BINTIME");
			  }
			else
			  listcommands(argv[iterm],p);
		      }
		    else
		      listcommands(argv[iterm],p);
		  }
		else {
		  i--;
		  c[cn].AddNoise->bintime_type = PERTYPE_FIX;
		  c[cn].AddNoise->bintime_fix = -1.0;
		}
	      }
	    else {
	      i--;
	      c[cn].AddNoise->bintime_type = PERTYPE_FIX;
	      c[cn].AddNoise->bintime_fix = -1.0;
	    }
	  }
	  else if(!strcmp(argv[i],"matern")) {
	    c[cn].require_distinct = 1;
	    c[cn].require_sort = 1;
	    c[cn].AddNoise->noise_type = VARTOOLS_ADDNOISE_COVAR_MATERN;
	    i++;
	    if(i < argc)
	      {
		if(!strcmp(argv[i],"nu"))
		  {
		    i++;
		    if(i < argc)
		      {
			if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			  {
			    c[cn].AddNoise->nu_r_type = PERTYPE_FIX;
			    i++;
			    if(i < argc)
			      c[cn].AddNoise->nu_r_fix = atof(argv[i]);
			    else
			      listcommands(argv[iterm],p);
			  }
			else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			  {
			    c[cn].AddNoise->nu_r_type = PERTYPE_SPECIFIED;
			    k = 0;
			    i++;
			    if(i < argc)
			      {
				if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
				  {
				    i++;
				    if(i < argc)
				      k = atoi(argv[i]);
				    else
				      listcommands(argv[iterm],p);
				  }
				else
				  i--;
			      }
			    else
			      i--;
			    RegisterDataFromInputList(p, 
						      (void *)(&(c[cn].AddNoise->nu_r)),
						      VARTOOLS_TYPE_DOUBLE,
						      1, cn, 0, 0, NULL, k,
						      "ADDNOISE_NU");
			  }
			else
			  listcommands(argv[iterm],p);
		      }
		    else
		      listcommands(argv[iterm],p);
		  }
		else
		  listcommands(argv[iterm],p);
	      }
	    else
	      listcommands(argv[iterm],p);
	    i++;
	    if(i < argc)
	      {
		if(!strcmp(argv[i],"rho"))
		  {
		    i++;
		    if(i < argc)
		      {
			if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			  {
			    c[cn].AddNoise->rho_r_type = PERTYPE_FIX;
			    i++;
			    if(i < argc)
			      c[cn].AddNoise->rho_r_fix = atof(argv[i]);
			    else
			      listcommands(argv[iterm],p);
			  }
			else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			  {
			    c[cn].AddNoise->rho_r_type = PERTYPE_SPECIFIED;
			    k = 0;
			    i++;
			    if(i < argc)
			      {
				if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
				  {
				    i++;
				    if(i < argc)
				      k = atoi(argv[i]);
				    else
				      listcommands(argv[iterm],p);
				  }
				else
				  i--;
			      }
			    else
			      i--;
			    RegisterDataFromInputList(p, 
						      (void *)(&(c[cn].AddNoise->rho_r)),
						      VARTOOLS_TYPE_DOUBLE,
						      1, cn, 0, 0, NULL, k,
						      "ADDNOISE_RHO");
			  }
			else
			  listcommands(argv[iterm],p);
		      }
		    else
		      listcommands(argv[iterm],p);
		  }
		else
		  listcommands(argv[iterm],p);
	      }
	    else
	      listcommands(argv[iterm],p);
	    i++;
	    if(i < argc)
	      {
		if(!strncmp(argv[i],"sig_red",7) && strlen(argv[i]) == 7)
		  {
		    i++;
		    if(i < argc)
		      {
			if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			  {
			    c[cn].AddNoise->sig_r_type = PERTYPE_FIX;
			    i++;
			    if(i < argc)
			      c[cn].AddNoise->sig_r_fix = atof(argv[i]);
			    else
			      listcommands(argv[iterm],p);
			  }
			else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			  {
			    c[cn].AddNoise->sig_r_type = PERTYPE_SPECIFIED;
			    k = 0;
			    i++;
			    if(i < argc) {
			      if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
				{
				  i++;
				  if(i < argc)
				    k = atoi(argv[i]);
				  else
				    listcommands(argv[iterm],p);
				}
			      else
				i--;
			    }
			    else
			      i--;
			    RegisterDataFromInputList(p, 
						      (void *)(&(c[cn].AddNoise->sig_r)),
						      VARTOOLS_TYPE_DOUBLE,
						      1, cn, 0, 0, NULL, k,
						      "ADDNOISE_SIGMA_RED");
			  }
			else
			  listcommands(argv[iterm],p);
		      }
		    else
		      listcommands(argv[iterm],p);
		  }
		else
		  listcommands(argv[iterm],p);
	      }
	    else
	      listcommands(argv[iterm],p);
	    i++;
	    if(i < argc)
	      {
		if(!strncmp(argv[i],"sig_white",9) && strlen(argv[i]) == 9)
		  {
		    i++;
		    if(i < argc)
		      {
			if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			  {
			    c[cn].AddNoise->sig_w_type = PERTYPE_FIX;
			    i++;
			    if(i < argc)
			      c[cn].AddNoise->sig_w_fix = atof(argv[i]);
			    else
			      listcommands(argv[iterm],p);
			  }
			else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			  {
			    c[cn].AddNoise->sig_w_type = PERTYPE_SPECIFIED;
			    k = 0;
			    i++;
			    if(i < argc) {
			      if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
				{
				  i++;
				  if(i < argc)
				    k = atoi(argv[i]);
				  else
				    listcommands(argv[iterm],p);
				}
			      else
				i--;
			    }
			    else
			      i--;
			    RegisterDataFromInputList(p, 
						      (void *)(&(c[cn].AddNoise->sig_w)),
						      VARTOOLS_TYPE_DOUBLE,
						      1, cn, 0, 0, NULL, k,
						      "ADDNOISE_SIGMA_WHITE");
			  }
			else
			  listcommands(argv[iterm],p);
		      }
		    else
		      listcommands(argv[iterm],p);
		  }
		else
		  listcommands(argv[iterm],p);
	      }
	    else
	      listcommands(argv[iterm],p);

	    i++;
	    if(i < argc)
	      {
		if(!strcmp(argv[i],"bintime"))
		  {
		    i++;
		    if(i < argc)
		      {
			if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			  {
			    c[cn].AddNoise->bintime_type = PERTYPE_FIX;
			    i++;
			    if(i < argc)
			      c[cn].AddNoise->bintime_fix = atof(argv[i]);
			    else
			      listcommands(argv[iterm],p);
			  }
			else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			  {
			    c[cn].AddNoise->bintime_type = PERTYPE_SPECIFIED;
			    k = 0;
			    i++;
			    if(i < argc) {
			      if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
				{
				  i++;
				  if(i < argc)
				    k = atoi(argv[i]);
				  else
				    listcommands(argv[iterm],p);
				}
			      else
				i--;
			    }
			    else
			      i--;
			    RegisterDataFromInputList(p, 
						      (void *)(&(c[cn].AddNoise->bintime)),
						      VARTOOLS_TYPE_DOUBLE,
						      1, cn, 0, 0, NULL, k,
						      "ADDNOISE_BINTIME");
			  }
			else
			  listcommands(argv[iterm],p);
		      }
		    else
		      listcommands(argv[iterm],p);
		  }
		else {
		  i--;
		  c[cn].AddNoise->bintime_type = PERTYPE_FIX;
		  c[cn].AddNoise->bintime_fix = -1.0;
		}
	      }
	    else {
	      i--;
	      c[cn].AddNoise->bintime_type = PERTYPE_FIX;
	      c[cn].AddNoise->bintime_fix = -1.0;
	    }
	  }
#ifdef _HAVE_GSL
	  else if(!strcmp(argv[i],"wavelet")) {
	    c[cn].require_distinct = 1;
	    c[cn].require_sort = 1;
	    c[cn].AddNoise->noise_type = VARTOOLS_ADDNOISE_WAVELET;
	    i++;
	    if(i < argc)
	      {
		if(!strncmp(argv[i],"gamma",5) && strlen(argv[i]) == 5)
		  {
		    i++;
		    if(i < argc)
		      {
			if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			  {
			    c[cn].AddNoise->gammaval_type = PERTYPE_FIX;
			    i++;
			    if(i < argc)
			      c[cn].AddNoise->gammaval_fix = atof(argv[i]);
			    else
			      listcommands(argv[iterm],p);
			  }
			else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			  {
			    c[cn].AddNoise->gammaval_type = PERTYPE_SPECIFIED;
			    k = 0;
			    i++;
			    if(i < argc)
			      {
				if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
				  {
				    i++;
				    if(i < argc)
				      k = atoi(argv[i]);
				    else
				      listcommands(argv[iterm],p);
				  }
				else
				  i--;
			      }
			    else
			      i--;
			    RegisterDataFromInputList(p, 
						      (void *)(&(c[cn].AddNoise->gammaval)),
						      VARTOOLS_TYPE_DOUBLE,
						      1, cn, 0, 0, NULL, k,
						      "ADDNOISE_GAMMA");
			  }
			else
			  listcommands(argv[iterm],p);
		      }
		    else
		      listcommands(argv[iterm],p);
		  }
		else
		  listcommands(argv[iterm],p);
	      }
	    else
	      listcommands(argv[iterm],p);
	    i++;
	    if(i < argc)
	      {
		if(!strncmp(argv[i],"sig_red",7) && strlen(argv[i]) == 7)
		  {
		    i++;
		    if(i < argc)
		      {
			if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			  {
			    c[cn].AddNoise->sig_r_type = PERTYPE_FIX;
			    i++;
			    if(i < argc)
			      c[cn].AddNoise->sig_r_fix = atof(argv[i]);
			    else
			      listcommands(argv[iterm],p);
			  }
			else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			  {
			    c[cn].AddNoise->sig_r_type = PERTYPE_SPECIFIED;
			    k = 0;
			    i++;
			    if(i < argc) {
			      if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
				{
				  i++;
				  if(i < argc)
				    k = atoi(argv[i]);
				  else
				    listcommands(argv[iterm],p);
				}
			      else
				i--;
			    }
			    else
			      i--;
			    RegisterDataFromInputList(p, 
						      (void *)(&(c[cn].AddNoise->sig_r)),
						      VARTOOLS_TYPE_DOUBLE,
						      1, cn, 0, 0, NULL, k,
						      "ADDNOISE_SIGMA_RED");
			  }
			else
			  listcommands(argv[iterm],p);
		      }
		    else
		      listcommands(argv[iterm],p);
		  }
		else
		  listcommands(argv[iterm],p);
	      }
	    else
	      listcommands(argv[iterm],p);
	    i++;
	    if(i < argc)
	      {
		if(!strncmp(argv[i],"sig_white",9) && strlen(argv[i]) == 9)
		  {
		    i++;
		    if(i < argc)
		      {
			if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			  {
			    c[cn].AddNoise->sig_w_type = PERTYPE_FIX;
			    i++;
			    if(i < argc)
			      c[cn].AddNoise->sig_w_fix = atof(argv[i]);
			    else
			      listcommands(argv[iterm],p);
			  }
			else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			  {
			    c[cn].AddNoise->sig_w_type = PERTYPE_SPECIFIED;
			    k = 0;
			    i++;
			    if(i < argc) {
			      if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
				{
				  i++;
				  if(i < argc)
				    k = atoi(argv[i]);
				  else
				    listcommands(argv[iterm],p);
				}
			      else
				i--;
			    }
			    else
			      i--;
			    RegisterDataFromInputList(p, 
						      (void *)(&(c[cn].AddNoise->sig_w)),
						      VARTOOLS_TYPE_DOUBLE,
						      1, cn, 0, 0, NULL, k,
						      "ADDNOISE_SIGMA_WHITE");
			  }
			else
			  listcommands(argv[iterm],p);
		      }
		    else
		      listcommands(argv[iterm],p);
		  }
		else
		  listcommands(argv[iterm],p);
	      }
	    else
	      listcommands(argv[iterm],p);
	  }
#endif
	  else
	    listcommands(argv[iterm],p);
	  cn++;
	}

      /* -alarm */
      else if(!strncmp(argv[i],"-alarm",6) && strlen(argv[i]) == 6)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].require_sort = 1;
	  c[cn].require_distinct = 1;
	  c[cn].cnum = CNUM_ALARM;
	  if((c[cn].Alarm = (_Alarm *) malloc(sizeof(_Alarm))) == NULL)
	    error(ERR_MEMALLOC);
	  cn++;
	}
      
      /* -aov ["Nbin" Nbin] minp maxp subsample finetune Npeaks operiodogram [outdir] [\"whiten\"] [\"clip\" clip clipiter] [\"uselog\"] [\"fixperiodSNR\" <\"aov\" | \"ls\" | \"injectharm\" | \"fix\" period | \"list\" [\"column\" col] | \"fixcolumn\" <colname | colnum>>] */
      else if(!strncmp(argv[i],"-aov",4) && strlen(argv[i]) == 4)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].require_sort = 1;
	  c[cn].require_distinct = 1;
	  c[cn].cnum = CNUM_AOV;
	  if((c[cn].Aov = (_Aov *) malloc(sizeof(_Aov))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"Nbin",4) && strlen(argv[i]) == 4)
		{
		  i++;
		  if(i < argc)
		    c[cn].Aov->Nbin = atoi(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Aov->minp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		{
		  c[cn].Aov->Nbin = 0;
		  c[cn].Aov->minp = atof(argv[i]);
		}
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Aov->maxp = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Aov->subsample = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Aov->finetune = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Aov->Npeaks = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Aov->operiodogram = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].Aov->operiodogram)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].Aov->outdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	    }
	  sprintf(c[cn].Aov->suffix,".aov");
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"whiten",6) && strlen(argv[i]) == 6)
		{
		  c[cn].Aov->whiten = 1;
		}
	      else
		{
		  i--;
		  c[cn].Aov->whiten = 0;
		}
	    }
	  else
	    {
	      i--;
	      c[cn].Aov->whiten = 0;
	    }

	  c[cn].Aov->clip = 5.0;
	  c[cn].Aov->clipiter = 1;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"clip",4) && strlen(argv[i]) == 4)
		{
		  i++;
		  if(i < argc)
		    c[cn].Aov->clip = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Aov->clipiter = atoi(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  c[cn].Aov->uselog = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"uselog",6) && strlen(argv[i]) == 6)
		c[cn].Aov->uselog = 1;
	      else
		i--;
	    }
	  else
	    i--;

	  c[cn].Aov->fixperiodSNR = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"fixperiodSNR",12) && strlen(argv[i]) == 12)
		{
		  c[cn].Aov->fixperiodSNR = 1;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"aov",3) && strlen(argv[i]) == 3)
			{
			  c[cn].Aov->fixperiodSNR_pertype = PERTYPE_AOV;
			  m = -1;
			  /* Get the index for the last aov command */
			  for(l = 0; l < cn; l++)
			    {
			      if(c[l].cnum == CNUM_AOV || c[l].cnum == CNUM_HARMAOV)
				m = l;
			    }
			  if(m < 0)
			    error(ERR_KILLHARM_NOAOV);
			  else
			    c[cn].Aov->fixperiodSNR_lastaovindex = m;
			}
		      else if(!strncmp(argv[i],"ls",2) && strlen(argv[i]) == 2)
			{
			  c[cn].Aov->fixperiodSNR_pertype = PERTYPE_LS;
			  m = -1;
			  /* Get the index for the last aov command */
			  for(l = 0; l < cn; l++)
			    {
			      if(c[l].cnum == CNUM_LS)
				m = l;
			    }
			  if(m < 0 || m == cn)
			    error(ERR_KILLHARM_NOLS);
			  else
			    c[cn].Aov->fixperiodSNR_lastaovindex = m;
			}
		      else if(!strncmp(argv[i],"injectharm",10) && strlen(argv[i]) == 10)
			{
			  c[cn].Aov->fixperiodSNR_pertype = PERTYPE_INJECTHARM;
			  m = -1;
			  /* Get the index for the last aov command */
			  for(l = 0; l < cn; l++)
			    {
			      if(c[l].cnum == CNUM_INJECTHARM)
				m = l;
			    }
			  if(m < 0)
			    error(ERR_KILLHARM_NOINJECTHARM);
			  else
			    c[cn].Aov->fixperiodSNR_lastaovindex = m;
			}
		      else if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			{
			  c[cn].Aov->fixperiodSNR_pertype = PERTYPE_FIX;
			  i++;
			  if(i < argc)
			    {
			      c[cn].Aov->fixperiodSNR_fixedperiod = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"fixcolumn",9) && strlen(argv[i]) == 9)
			{
			  c[cn].Aov->fixperiodSNR_pertype = PERTYPE_FIXCOLUMN;
			  i++;
			  if(i < argc)
			    increaselinkedcols(p, &(c[cn].Aov->fixperiodSNR_linkedcolumn), argv[i], cn);
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			{
			  c[cn].Aov->fixperiodSNR_pertype = PERTYPE_SPECIFIED;
			  k = 0;
			  i++;
			  if(i < argc) {
			    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    } else i--;
			  } else i--;
			  RegisterDataFromInputList(p,
						    (void *) (&(c[cn].Aov->fixperiodSNR_periods)),
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "AOV_FIXPERIODSNR_PERIOD");
			}
		      else
			listcommands(argv[iterm],p);
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  cn++;
	}

      /* -aov_harm Nharm minp maxp subsample finetune Npeaks operiodogram [outdir] [\"whiten\"] [\"clip\" clip clipiter] [\"fixperiodSNR\" <\"aov\" | \"ls\" | \"injectharm\" | \"fix\" period | \"list\" [\"column\" col] | \"fixcolumn\" <colname | colnum>>] */
      else if(!strncmp(argv[i],"-aov_harm",9) && strlen(argv[i]) == 9)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].require_sort = 1;
	  c[cn].require_distinct = 1;
	  c[cn].cnum = CNUM_HARMAOV;
	  if((c[cn].AovHarm = (_AovHarm *) malloc(sizeof(_AovHarm))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    c[cn].AovHarm->Nharm = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].AovHarm->minp = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].AovHarm->maxp = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].AovHarm->subsample = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].AovHarm->finetune = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].AovHarm->Npeaks = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].AovHarm->operiodogram = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].AovHarm->operiodogram)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].AovHarm->outdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	    }
	  sprintf(c[cn].AovHarm->suffix,".aov_harm");
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"whiten",6) && strlen(argv[i]) == 6)
		{
		  c[cn].AovHarm->whiten = 1;
		}
	      else
		{
		  i--;
		  c[cn].AovHarm->whiten = 0;
		}
	    }
	  else
	    {
	      i--;
	      c[cn].AovHarm->whiten = 0;
	    }
	  c[cn].AovHarm->clip = 5.0;
	  c[cn].AovHarm->clipiter = 1;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"clip",4) && strlen(argv[i]) == 4)
		{
		  i++;
		  if(i < argc)
		    c[cn].AovHarm->clip = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].AovHarm->clipiter = atoi(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;

	  c[cn].AovHarm->fixperiodSNR = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"fixperiodSNR",12) && strlen(argv[i]) == 12)
		{
		  c[cn].AovHarm->fixperiodSNR = 1;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"aov",3) && strlen(argv[i]) == 3)
			{
			  c[cn].AovHarm->fixperiodSNR_pertype = PERTYPE_AOV;
			  m = -1;
			  /* Get the index for the last aov command */
			  for(l = 0; l < cn; l++)
			    {
			      if(c[l].cnum == CNUM_AOV || c[l].cnum == CNUM_HARMAOV)
				m = l;
			    }
			  if(m < 0)
			    error(ERR_KILLHARM_NOAOV);
			  else
			    c[cn].AovHarm->fixperiodSNR_lastaovindex = m;
			}
		      else if(!strncmp(argv[i],"ls",2) && strlen(argv[i]) == 2)
			{
			  c[cn].AovHarm->fixperiodSNR_pertype = PERTYPE_LS;
			  m = -1;
			  /* Get the index for the last aov command */
			  for(l = 0; l < cn; l++)
			    {
			      if(c[l].cnum == CNUM_LS)
				m = l;
			    }
			  if(m < 0 || m == cn)
			    error(ERR_KILLHARM_NOLS);
			  else
			    c[cn].AovHarm->fixperiodSNR_lastaovindex = m;
			}
		      else if(!strncmp(argv[i],"injectharm",10) && strlen(argv[i]) == 10)
			{
			  c[cn].AovHarm->fixperiodSNR_pertype = PERTYPE_INJECTHARM;
			  m = -1;
			  /* Get the index for the last aov command */
			  for(l = 0; l < cn; l++)
			    {
			      if(c[l].cnum == CNUM_INJECTHARM)
				m = l;
			    }
			  if(m < 0)
			    error(ERR_KILLHARM_NOINJECTHARM);
			  else
			    c[cn].AovHarm->fixperiodSNR_lastaovindex = m;
			}
		      else if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			{
			  c[cn].AovHarm->fixperiodSNR_pertype = PERTYPE_FIX;
			  i++;
			  if(i < argc)
			    {
			      c[cn].AovHarm->fixperiodSNR_fixedperiod = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"fixcolumn",9) && strlen(argv[i]) == 9)
			{
			  c[cn].AovHarm->fixperiodSNR_pertype = PERTYPE_FIXCOLUMN;
			  i++;
			  if(i < argc)
			    increaselinkedcols(p, &(c[cn].AovHarm->fixperiodSNR_linkedcolumn), argv[i], cn);
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			{
			  c[cn].AovHarm->fixperiodSNR_pertype = PERTYPE_SPECIFIED;
			  k = 0;
			  i++;
			  if(i < argc) {
			    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    } else i--;
			  } else i--;
			  RegisterDataFromInputList(p,
						    (void *) (&(c[cn].AovHarm->fixperiodSNR_periods)),
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "AOVHARM_FIXPERIODSNR_PERIOD");
			}
		      else
			listcommands(argv[iterm],p);
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  cn++;
	}

      /* -autocorrelation start stop step match-tolerance outdir */
      else if(!strncmp(argv[i],"-autocorrelation",16) && strlen(argv[i]) == 16)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].require_sort = 1;
	  c[cn].cnum = CNUM_AUTOCORR;
	  if((c[cn].Autocorr = (_Autocorr *) malloc(sizeof(_Autocorr))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    c[cn].Autocorr->start = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Autocorr->stop = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Autocorr->step = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    sprintf(c[cn].Autocorr->outdir,"%s",argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  sprintf(c[cn].Autocorr->suffix,".autocorr");
	  cn++;
	}
      
      /*-converttime
    <"input" <"mjd" | "jd" | "hjd" | "bjd" >>
    ["inputsubtract" value] ["inputsys-tdb" | "inputsys-utc"]
    <"output" <"mjd" | "jd" | "hjd" | "bjd" >>
    ["outputsubtract" value] ["outputsys-tdb" | "outputsys-utc"]
    ["radec" <"list" ["column" col] | "fix" raval decval> ["epoch" epoch]]
    ["ppm" <"list" ["column" col] | "fix" mu_ra mu_dec>]
    ["input-radec" <"list" ["column" col] | "fix" raval decval> ["epoch" epoch]]
    ["input-ppm" <"list" ["column" col] | "fix" mu_ra mu_dec>] 
    ["ephemfile" file] ["leapsecfile" file] ["planetdatafile" file]
    ["observatory" < code | "show-codes">
        | "coords"
            <"fix" latitude[deg] longitude[deg_W] altitude[m]
            | "list" ["column" collat collong colalt]
            | "fromlc" collat collong colalt>] */
      else if(!strncmp(argv[i],"-converttime",12) && strlen(argv[i]) == 12)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_CONVERTTIME;
	  if((c[cn].ConvertTime = (_ConvertTime *) malloc(sizeof(_ConvertTime))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    {
	      if(strncmp(argv[i],"input",5) || strlen(argv[i]) != 5)
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"mjd",3) && strlen(argv[i]) == 3)
		{
		  c[cn].ConvertTime->inputtimetype = TIMETYPE_MJD;
		}
	      else if(!strncmp(argv[i],"jd",2) && strlen(argv[i]) == 2)
		{
		  c[cn].ConvertTime->inputtimetype = TIMETYPE_JD;
		}
	      else if(!strncmp(argv[i],"hjd",3) && strlen(argv[i]) == 3)
		{
		  c[cn].ConvertTime->inputtimetype = TIMETYPE_HJD;
		}
	      else if(!strncmp(argv[i],"bjd",3) && strlen(argv[i]) == 3)
		{
		  c[cn].ConvertTime->inputtimetype = TIMETYPE_BJD;
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  c[cn].ConvertTime->inputsubtractval = 0.;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"inputsubtract",13) && strlen(argv[i]) == 13)
		{
		  i++;
		  if(i < argc)
		    {
		      c[cn].ConvertTime->inputsubtractval = atof(argv[i]);
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  c[cn].ConvertTime->inputsys = TIMESYSTEM_UTC;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"inputsys-tdb",12) && strlen(argv[i]) == 12)
		c[cn].ConvertTime->inputsys = TIMESYSTEM_TDB;
	      else if(!strncmp(argv[i],"inputsys-utc",12) && strlen(argv[i]) == 12)
		c[cn].ConvertTime->inputsys = TIMESYSTEM_UTC;
	      else
		i--;
	    }
	  else
	    i--;
	  i++;
	  if(i >= argc || (strncmp(argv[i],"output",6) || strlen(argv[i]) != 6))
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"mjd",3) && strlen(argv[i]) == 3)
		{
		  c[cn].ConvertTime->outputtimetype = TIMETYPE_MJD;
		}
	      else if(!strncmp(argv[i],"jd",2) && strlen(argv[i]) == 2)
		{
		  c[cn].ConvertTime->outputtimetype = TIMETYPE_JD;
		}
	      else if(!strncmp(argv[i],"hjd",3) && strlen(argv[i]) == 3)
		{
		  c[cn].ConvertTime->outputtimetype = TIMETYPE_HJD;
		}
	      else if(!strncmp(argv[i],"bjd",3) && strlen(argv[i]) == 3)
		{
		  c[cn].ConvertTime->outputtimetype = TIMETYPE_BJD;
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  c[cn].ConvertTime->outputsubtractval = 0.;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"outputsubtract",14) && strlen(argv[i]) == 14)
		{
		  i++;
		  if(i < argc)
		    {
		      c[cn].ConvertTime->outputsubtractval = atof(argv[i]);
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  c[cn].ConvertTime->outputsys = c[cn].ConvertTime->inputsys;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"outputsys-tdb",13) && strlen(argv[i]) == 13)
		c[cn].ConvertTime->outputsys = TIMESYSTEM_TDB;
	      else if(!strncmp(argv[i],"outputsys-utc",13) && strlen(argv[i]) == 13)
		c[cn].ConvertTime->outputsys = TIMESYSTEM_UTC;
	      else
		i--;
	    }
	  else
	    i--;

	  c[cn].ConvertTime->useradec = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"radec",5) && strlen(argv[i]) == 5)
		{
		  c[cn].ConvertTime->useradec = 1;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			{
			  c[cn].ConvertTime->radec_source = VARTOOLS_SOURCE_INLIST;
			  k = 0;
			  i++;
			  if(i < argc) {
			    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    }
			    else
			      i--;
			  }
			  else i--;
			  RegisterDataFromInputList(p, 
						    (void *)(&(c[cn].ConvertTime->ravals)), 
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "CONVERTTIME_RA");
			  if(k > 0) k++;
			  RegisterDataFromInputList(p, 
						    (void *)(&(c[cn].ConvertTime->decvals)), 
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "CONVERTTIME_DEC");
			}
		      else if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			{
			  c[cn].ConvertTime->radec_source = VARTOOLS_SOURCE_FIXED;
			  i++;
			  if(i < argc)
			    {
			      c[cn].ConvertTime->raval_fix = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			  i++;
			  if(i < argc)
			    {
			      c[cn].ConvertTime->decval_fix = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			listcommands(argv[iterm],p);
		    }
		  else
		    listcommands(argv[iterm],p);
		  c[cn].ConvertTime->radecepoch = 2000.0;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"epoch",5) && strlen(argv[i]) == 5)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].ConvertTime->radecepoch = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			i--;
		    }
		  else
		    i--;
		}
	      else
		i--;
	    }
	  else
	    i--;
	  c[cn].ConvertTime->useppm = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"ppm",3) && strlen(argv[i]) == 3)
		{
		  c[cn].ConvertTime->useppm = 1;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			{
			  c[cn].ConvertTime->ppm_source = VARTOOLS_SOURCE_INLIST;
			  k = 0;
			  i++;
			  if(i < argc) {
			    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    }
			    else
			      i--;
			  }
			  else i--;
			  RegisterDataFromInputList(p, 
						    (void *)(&(c[cn].ConvertTime->ppm_mu_ra_vals)), 
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "CONVERTTIME_PPM_MU_RA");
			  if(k > 0) k++;
			  RegisterDataFromInputList(p, 
						    (void *)(&(c[cn].ConvertTime->ppm_mu_dec_vals)), 
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "CONVERTTIME_PPM_MU_DEC");
			}
		      else if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			{
			  c[cn].ConvertTime->ppm_source = VARTOOLS_SOURCE_FIXED;
			  i++;
			  if(i < argc)
			    {
			      c[cn].ConvertTime->ppm_mu_ra_fix = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			  i++;
			  if(i < argc)
			    c[cn].ConvertTime->ppm_mu_dec_fix = atof(argv[i]);
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			listcommands(argv[iterm],p);
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;

	  c[cn].ConvertTime->useinput_radec = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"input-radec",11) && strlen(argv[i]) == 11)
		{
		  c[cn].ConvertTime->useinput_radec = 1;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			{
			  c[cn].ConvertTime->inputradec_source = VARTOOLS_SOURCE_INLIST;
			  k = 0;
			  i++;
			  if(i < argc) {
			    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    }
			    else
			      i--;
			  }
			  else i--;
			  RegisterDataFromInputList(p, 
						    (void *)(&(c[cn].ConvertTime->inputravals)), 
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "CONVERTTIME_INPUT_RA");
			  if(k > 0) k++;
			  RegisterDataFromInputList(p, 
						    (void *)(&(c[cn].ConvertTime->inputdecvals)), 
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "CONVERTTIME_INPUT_DEC");

			}
		      else if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			{
			  c[cn].ConvertTime->inputradec_source = VARTOOLS_SOURCE_FIXED;
			  i++;
			  if(i < argc)
			    {
			      c[cn].ConvertTime->inputraval_fix = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			  i++;
			  if(i < argc)
			    {
			      c[cn].ConvertTime->inputdecval_fix = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			listcommands(argv[iterm],p);
		    }
		  else
		    listcommands(argv[iterm],p);
		  c[cn].ConvertTime->inputradecepoch = 2000.0;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"epoch",5) && strlen(argv[i]) == 5)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].ConvertTime->inputradecepoch = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			i--;
		    }
		  else
		    i--;
		}
	      else
		i--;
	    }
	  else
	    i--;
	  c[cn].ConvertTime->useinputppm = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"input-ppm",9) && strlen(argv[i]) == 9)
		{
		  c[cn].ConvertTime->useinputppm = 1;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			{
			  c[cn].ConvertTime->inputppm_source = VARTOOLS_SOURCE_INLIST;
			  k = 0;
			  i++;
			  if(i < argc) {
			    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    }
			    else
			      i--;
			  }
			  else i--;
			  RegisterDataFromInputList(p, 
						    (void *)(&(c[cn].ConvertTime->inputppm_mu_ra_vals)), 
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "CONVERTTIME_INPUT_PPM_MU_RA");
			  if(k > 0) k++;
			  RegisterDataFromInputList(p, 
						    (void *)(&(c[cn].ConvertTime->inputppm_mu_dec_vals)), 
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "CONVERTTIME_INPUT_PPM_MU_DEC");
			}
		      else if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			{
			  c[cn].ConvertTime->inputppm_source = VARTOOLS_SOURCE_FIXED;
			  i++;
			  if(i < argc)
			    {
			      c[cn].ConvertTime->inputppm_mu_ra_fix = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			  i++;
			  if(i < argc)
			    c[cn].ConvertTime->inputppm_mu_dec_fix = atof(argv[i]);
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			listcommands(argv[iterm],p);
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
#ifdef _HAVE_CSPICE	  
	  c[cn].ConvertTime->ephemfile[0] = '\0';
	  c[cn].ConvertTime->leapsecfile[0] = '\0';
	  c[cn].ConvertTime->planetdatafile[0] = '\0';
	  i++;
	  if(i < argc) {
	    if(!strncmp(argv[i],"ephemfile",9) && strlen(argv[i]) == 9) {
	      i++;
	      if(i < argc)
		sprintf(c[cn].ConvertTime->ephemfile,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	    } else i--;
	  } else i--;
	  i++;
	  if(i < argc) {
	    if(!strncmp(argv[i],"leapsecfile",11) && strlen(argv[i]) == 11) {
	      i++;
	      if(i < argc)
		sprintf(c[cn].ConvertTime->leapsecfile,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	    } else i--;
	  } else i--;
	  i++;
	  if(i < argc) {
	    if(!strncmp(argv[i],"planetdatafile",14) && strlen(argv[i]) == 14) {
	      i++;
	      if(i < argc)
		sprintf(c[cn].ConvertTime->planetdatafile,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	    } else i--;
	  } else i--;
	  c[cn].ConvertTime->source_obs_coords = VARTOOLS_SOURCE_NONE;
	  i++;
	  if(i < argc) {
	    if(!strncmp(argv[i],"observatory",11) && strlen(argv[i]) == 11) {
	      c[cn].ConvertTime->source_obs_coords = VARTOOLS_SOURCE_FIXED;
	      i++;
	      if(i < argc) {
		ParseObservatoryCode(argv[i],
				     &(c[cn].ConvertTime->obslat_fixval),
				     &(c[cn].ConvertTime->obslong_fixval),
				     &(c[cn].ConvertTime->obsalt_fixval));
	      } else
		listcommands(argv[iterm],p);
	    }
	    else if(!strncmp(argv[i],"coords",6) && strlen(argv[i]) == 6) {
	      i++;
	      if(i < argc) {
		if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3) {
		  c[cn].ConvertTime->source_obs_coords = VARTOOLS_SOURCE_FIXED;
		  i++;
		  if(i < argc)
		    c[cn].ConvertTime->obslat_fixval = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].ConvertTime->obslong_fixval = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].ConvertTime->obsalt_fixval = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
		else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4) {
		  c[cn].ConvertTime->source_obs_coords = VARTOOLS_SOURCE_INLIST;
		  j = 0; k = 0; l = 0;
		  i++;
		  if(i < argc) {
		    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
		      i++;
		      if(i < argc)
			j = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		      i++;
		      if(i < argc)
			k = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		      i++;
		      if(i < argc)
			l = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		    } else i--;
		  } else i--;
		  RegisterDataFromInputList(p,
					    (void *) &(c[cn].ConvertTime->obslat_listvals),
					    VARTOOLS_TYPE_DOUBLE,
					    0, cn, 0, 0, NULL, j, 
					    "CONVERTTIME_OBSLAT");
		  RegisterDataFromInputList(p,
					    (void *) &(c[cn].ConvertTime->obslong_listvals),
					    VARTOOLS_TYPE_DOUBLE,
					    0, cn, 0, 0, NULL, k, 
					    "CONVERTTIME_OBSLONG");
		  RegisterDataFromInputList(p,
					    (void *) &(c[cn].ConvertTime->obsalt_listvals),
					    VARTOOLS_TYPE_DOUBLE,
					    0, cn, 0, 0, NULL, l, 
					    "CONVERTTIME_OBSALT");
		}
		else if(!strncmp(argv[i],"fromlc",6) && strlen(argv[i]) == 6) {
		  c[cn].ConvertTime->source_obs_coords = VARTOOLS_SOURCE_LC;
		  i++;
		  if(i < argc) {
		    c[cn].ConvertTime->obslat_lc_col = atoi(argv[i]);
		    RegisterDataFromLightCurve(p, 
					       &(c[cn].ConvertTime->obslat_lcvals), 
					       VARTOOLS_TYPE_DOUBLE, 0, 0, cn, 
					       0, 0, NULL, NULL,
					       c[cn].ConvertTime->obslat_lc_col,
					       "ObsLatitude");
		  }
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc) {
		    c[cn].ConvertTime->obslong_lc_col = atoi(argv[i]);
		    RegisterDataFromLightCurve(p, 
					       &(c[cn].ConvertTime->obslong_lcvals), 
					       VARTOOLS_TYPE_DOUBLE, 0, 0, cn, 
					       0, 0, NULL, NULL,
					       c[cn].ConvertTime->obslong_lc_col,
					       "ObsLongitude");
		  }
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc) {
		    c[cn].ConvertTime->obsalt_lc_col = atoi(argv[i]);
		    RegisterDataFromLightCurve(p, 
					       &(c[cn].ConvertTime->obsalt_lcvals), 
					       VARTOOLS_TYPE_DOUBLE, 0, 0, cn, 
					       0, 0, NULL, NULL,
					       c[cn].ConvertTime->obsalt_lc_col,
					       "ObsAltitude");
		  }
		  else
		    listcommands(argv[iterm],p);
		}
		else
		  listcommands(argv[iterm],p);
	      }
	      else
		listcommands(argv[iterm],p);
	    }
	    else i--;
	  }
	  else i--;

	  load_cspice_kernels(c[cn].ConvertTime->ephemfile,
			      c[cn].ConvertTime->leapsecfile,
			      c[cn].ConvertTime->planetdatafile,
			      c[cn].ConvertTime->inputtimetype,
			      c[cn].ConvertTime->outputtimetype,
			      c[cn].ConvertTime->inputsys,
			      c[cn].ConvertTime->outputsys);
	  if(!c[cn].ConvertTime->useradec) {
	    if(c[cn].ConvertTime->inputtimetype == TIMETYPE_BJD ||
	       c[cn].ConvertTime->outputtimetype == TIMETYPE_BJD ||
	       c[cn].ConvertTime->inputtimetype == TIMETYPE_HJD ||
	       c[cn].ConvertTime->outputtimetype == TIMETYPE_HJD)
	      error(ERR_CONVERTTIME_NORADEC);
	  }
	  
#endif
	  cn++;
	}


      /* -linfit function paramlist [\"modelvar\" varname] [\"correctlc\"] [\"omodel\" model_outdir]" */
      else if(!strcmp(argv[i],"-linfit"))
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_LINFIT;
	  if((c[cn].Linfit = (_Linfit *) malloc(sizeof(_Linfit))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(ParseLinfitCommand(&i, argc, argv, p, c[cn].Linfit))
	    listcommands(argv[iterm],p);
	  cn++;
	}

      /* -nonlinfit function paramlist [\"priors\" var1=-2ln(P1),var2=-2ln(P2),...] [\"constraints\" constraint_expr1,constraint_expr2,...] <\"amoeba\" [\"tolerance\" tol] [\"maxsteps\" max] | \"mcmc\" [\"Naccept\" Naccept | \"Nlinkstotal\" Notal] [\"fracburnin\" frac_burnin] [\"eps\" eps] [\"skipamoeba\"] [\"chainstats\" expression_list  stats_list] [\"maxmemstore\" max_mem_store] [\"outchains\" outdir [\"printevery\" printevery]] > [\"modelvar\" varname] [\"correctlc\"] [\"omodel\" model_outdir]" */
      else if(!strcmp(argv[i],"-nonlinfit"))
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_NONLINFIT;
	  if((c[cn].Nonlinfit = (_Nonlinfit *) malloc(sizeof(_Nonlinfit))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(ParseNonlinfitCommand(&i, argc, argv, p, c[cn].Nonlinfit))
	    listcommands(argv[iterm],p);
	  if(c[cn].Nonlinfit->use_covar) {
	    c[cn].require_sort = 1;
	    c[cn].require_distinct = 1;
	  }
	  cn++;
	}

      /* -if command1.... [-elif command...] [-else] [-fi] */
      else if(!strcmp(argv[i],"-if") ||
	      !strcmp(argv[i],"-elif") ||
	      !strcmp(argv[i],"-else") ||
	      !strcmp(argv[i],"-fi"))
	{
	  iterm = i;
	  increaseNcommands(p, &c);
	  c[cn].cnum = CNUM_IF;
	  if(ParseIfCommand(&i, argc, argv, cn, p, c))
	    listcommands(argv[iterm],p);
	  cn++;
	}

      /* -LS minp maxp subsample Npeaks operiodogram [outdir] [\"noGLS\"] [\"whiten\"] [\"clip\" clip clipiter] [\"fixperiodSNR\" <\"aov\" | \"ls\" | \"injectharm\" | \"fix\" period | \"list\" ["column" col] | \"fixcolumn\" <colname | colnum>>] [\"bootstrap\" Nbootstrap]*/
      else if(!strncmp(argv[i],"-LS",3) && strlen(argv[i]) == 3)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].require_sort = 1;
	  c[cn].require_distinct = 1;
	  c[cn].cnum = CNUM_LS;
	  if((c[cn].Ls = (_Ls *) malloc(sizeof(_Ls))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    c[cn].Ls->minp = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Ls->maxp = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Ls->subsample = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Ls->Npeaks = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Ls->operiodogram = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].Ls->operiodogram)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].Ls->outdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	    }
	  sprintf(c[cn].Ls->suffix,".ls");
	  c[cn].Ls->use_orig_ls = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strcmp(argv[i],"noGLS"))
		{
		  c[cn].Ls->use_orig_ls = 1;
		}
	      else
		{
		  i--;
		}
	    }
	  else
	    {
	      i--;
	    }
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"whiten",6) && strlen(argv[i]) == 6)
		{
		  c[cn].Ls->whiten = 1;
		}
	      else
		{
		  i--;
		  c[cn].Ls->whiten = 0;
		}
	    }
	  else
	    {
	      i--;
	      c[cn].Ls->whiten = 0;
	    }
	  c[cn].Ls->clip = 5.0;
	  c[cn].Ls->clipiter = 1;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"clip",4) && strlen(argv[i]) == 4)
		{
		  i++;
		  if(i < argc)
		    c[cn].Ls->clip = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Ls->clipiter = atoi(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  c[cn].Ls->fixperiodSNR = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"fixperiodSNR",12) && strlen(argv[i]) == 12)
		{
		  c[cn].Ls->fixperiodSNR = 1;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"aov",3) && strlen(argv[i]) == 3)
			{
			  c[cn].Ls->fixperiodSNR_pertype = PERTYPE_AOV;
			  m = -1;
			  /* Get the index for the last aov command */
			  for(l = 0; l < cn; l++)
			    {
			      if(c[l].cnum == CNUM_AOV || c[l].cnum == CNUM_HARMAOV)
				m = l;
			    }
			  if(m < 0)
			    error(ERR_KILLHARM_NOAOV);
			  else
			    c[cn].Ls->fixperiodSNR_lastaovindex = m;
			}
		      else if(!strncmp(argv[i],"ls",2) && strlen(argv[i]) == 2)
			{
			  c[cn].Ls->fixperiodSNR_pertype = PERTYPE_LS;
			  m = -1;
			  /* Get the index for the last aov command */
			  for(l = 0; l < cn; l++)
			    {
			      if(c[l].cnum == CNUM_LS)
				m = l;
			    }
			  if(m < 0 || m == cn)
			    error(ERR_KILLHARM_NOLS);
			  else
			    c[cn].Ls->fixperiodSNR_lastaovindex = m;
			}
		      else if(!strncmp(argv[i],"injectharm",10) && strlen(argv[i]) == 10)
			{
			  c[cn].Ls->fixperiodSNR_pertype = PERTYPE_INJECTHARM;
			  m = -1;
			  /* Get the index for the last aov command */
			  for(l = 0; l < cn; l++)
			    {
			      if(c[l].cnum == CNUM_INJECTHARM)
				m = l;
			    }
			  if(m < 0)
			    error(ERR_KILLHARM_NOINJECTHARM);
			  else
			    c[cn].Ls->fixperiodSNR_lastaovindex = m;
			}
		      else if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			{
			  c[cn].Ls->fixperiodSNR_pertype = PERTYPE_FIX;
			  i++;
			  if(i < argc)
			    {
			      c[cn].Ls->fixperiodSNR_fixedperiod = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"fixcolumn",9) && strlen(argv[i]) == 9)
			{
			  c[cn].Ls->fixperiodSNR_pertype = PERTYPE_FIXCOLUMN;
			  i++;
			  if(i < argc)
			    increaselinkedcols(p, &(c[cn].Ls->fixperiodSNR_linkedcolumn), argv[i], cn);
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			{
			  c[cn].Ls->fixperiodSNR_pertype = PERTYPE_SPECIFIED;
			  k = 0;
			  i++;
			  if(i < argc) {
			    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    } else i--;
			  } else i--;
			  RegisterDataFromInputList(p,
						    (void *) (&(c[cn].Ls->fixperiodSNR_periods)),
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "LS_FIXPERIODSNR_PERIOD");
			}
		      else
			listcommands(argv[iterm],p);
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  c[cn].Ls->dobootstrapfap = 0;
	  c[cn].Ls->Nbootstrap = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strcmp(argv[i],"bootstrap"))
		{
		  c[cn].Ls->dobootstrapfap = 1;
		  i++;
		  if(i < argc)
		    c[cn].Ls->Nbootstrap = atoi(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		{
		  i--;
		}
	    }
	  else
	    {
	      i--;
	    }
	  cn++;
	}

      /* -GetLSAmpThresh <\"ls\" | \"list\" [\"column\" col]> minp thresh <\"harm\" Nharm Nsubharm | \"file\" listfile> [\"noGLS\"] */
      else if(!strncmp(argv[i],"-GetLSAmpThresh",15) && strlen(argv[i]) == 15)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].require_sort = 1;
	  c[cn].require_distinct = 1;
	  c[cn].cnum = CNUM_GETLSAMPTHRESH;
	  if((c[cn].GetLSAmpThresh = (_GetLSAmpThresh *) malloc(sizeof(_GetLSAmpThresh))) == NULL)
	    error(ERR_MEMALLOC);
	  c[cn].GetLSAmpThresh->line_size = MAXLEN;
	  c[cn].GetLSAmpThresh->line = malloc(c[cn].GetLSAmpThresh->line_size);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"ls",2) && strlen(argv[i]) == 2)
		{
		  c[cn].GetLSAmpThresh->pertype = PERTYPE_LS;
		  m = -1;
		  /* Get the index for the last ls command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_LS)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOLS);
		  else
		    c[cn].GetLSAmpThresh->lastlsindex = m;
		}
	      else
		{
		  c[cn].GetLSAmpThresh->pertype = PERTYPE_SPECIFIED;
		  k = 0;
		  i++;
		  if(i < argc) {
		    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
		      i++;
		      if(i < argc)
			k = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		    } else i--;
		  } else i--;
		  RegisterDataFromInputList(p,
					    (void *) (&(c[cn].GetLSAmpThresh->period)),
					    VARTOOLS_TYPE_DOUBLE,
					    1, cn, 0, 0, NULL, k,
					    "GETLSAMPTHRESH_PERIOD");
		}
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].GetLSAmpThresh->minPer = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].GetLSAmpThresh->thresh = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"harm",4) && strlen(argv[i]) == 4)
		{
		  c[cn].GetLSAmpThresh->harm_specsigflag = 0;
		  i++;
		  if(i < argc)
		    c[cn].GetLSAmpThresh->Nharm = atoi(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].GetLSAmpThresh->Nsubharm = atoi(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"file",4) && strlen(argv[i]) == 4)
		{
		  c[cn].GetLSAmpThresh->harm_specsigflag = 1;
		  i++;
		  if(i < argc)
		    sprintf(c[cn].GetLSAmpThresh->listfilename,"%s",argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  if((c[cn].GetLSAmpThresh->listfile = fopen(c[cn].GetLSAmpThresh->listfilename,"r")) == NULL)
		    {
		      error2(ERR_FILENOTFOUND,c[cn].GetLSAmpThresh->listfilename);
		    }
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  c[cn].GetLSAmpThresh->use_orig_ls = 0;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"noGLS")) {
	      c[cn].GetLSAmpThresh->use_orig_ls = 1;
	    } else {
	      i--;
	    }
	  } else
	    i--;
	  cn++;
	}
      
      /* -o outdir : Ouput light curves */
      else if(!strncmp(argv[i],"-o",2) && strlen(argv[i]) == 2)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_OUTPUTLCS;
	  if((c[cn].Outputlcs = (_Outputlcs *) malloc(sizeof(_Outputlcs))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    sprintf(c[cn].Outputlcs->outdir,"%s",argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  c[cn].Outputlcs->useformat = 0;
	  c[cn].Outputlcs->usecolumnformat = 0;
	  c[cn].Outputlcs->sepchar = ' ';
	  c[cn].Outputlcs->variables = NULL;
	  c[cn].Outputlcs->Nvar = 0;
	  c[cn].Outputlcs->printfformats = NULL;
	  c[cn].Outputlcs->varnames = NULL;
	  c[cn].Outputlcs->outfits = 0;
	  c[cn].Outputlcs->noclobber = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strcmp(argv[i],"nameformat"))
		{
		  c[cn].Outputlcs->useformat = 1;
		  i++;
		  if(i < argc)
		    sprintf(c[cn].Outputlcs->format,"%s",argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  i++;
	  if(i < argc)
	    {
	      if(!strcmp(argv[i],"columnformat"))
		{
		  c[cn].Outputlcs->usecolumnformat = 1;
		  i++;
		  if(i < argc) {
		    sprintf(c[cn].Outputlcs->columnformat,"%s",argv[i]);
		    ParseOutputColumnFormat(c[cn].Outputlcs);
		  }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  i++;
	  if(i < argc)
	    {
	      if(!strcmp(argv[i],"delimiter"))
		{
		  i++;
		  if(i < argc) {
		    c[cn].Outputlcs->sepchar =  argv[i][0];
		  }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
#ifdef USECFITSIO
	  i++;
	  if(i < argc)
	    {
	      if(!strcmp(argv[i],"fits"))
		{
		  c[cn].Outputlcs->outfits = 1;
		}
	      else
		i--;
	    }
	  else
	    i--;
#endif
	  i++;
	  if(i < argc)
	    {
	      if(!strcmp(argv[i],"noclobber"))
		{
		  c[cn].Outputlcs->noclobber = 1;
		}
	      else
		i--;
	    }
	  else
	    i--;
	  cn++;
	}

      /* -Injectharm <\"list\" [\"column\" col] | \"fix\" per | \"rand\" minp maxp | \"logrand\" minp maxp | \"randfreq\" minf maxf | \"lograndfreq\" minf maxf> Nharm (<\"amplist\" [\"column\" col] | \"ampfix\" amp | \"amprand\" minamp maxamp | \"amplogrand\" minamp maxamp> [\"amprel\"] <\"phaselist\" [\"column\" col] | \"phasefix\" phase | \"phaserand\"> [\"phaserel\"])0...Nharm Nsubharm (<\"amplist\" [\"column\" col] | \"ampfix\" amp | \"amprand\" minamp maxamp | \"amplogrand\" minamp maxamp> [\"amprel\"] <\"phaselist\" [\"column\" col] | \"phasefix\" phase | \"phaserand\"> [\"phaserel\"])1...Nsubharm omodel [modeloutdir] */
      else if(!strncmp(argv[i],"-Injectharm",11) && strlen(argv[i]) == 11)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_INJECTHARM;
	  if((c[cn].Injectharm = (_Injectharm *) malloc(sizeof(_Injectharm))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
		{
		  c[cn].Injectharm->pertype = PERTYPE_SPECIFIED;
		  k = 0;
		  i++;
		  if(i < argc) {
		    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
		      i++;
		      if(i < argc) {
			k = atoi(argv[i]);
		      }
		      else
			listcommands(argv[iterm],p);
		    }
		    else i--;
		  } else i--;
		  RegisterDataFromInputList(p,
					    (void *)(&(c[cn].Injectharm->periods)),
					    VARTOOLS_TYPE_DOUBLE,
					    1, cn, 0, 0, NULL, k,
					    "INJECTHARM_PERIOD");
		}
	      else if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
		{
		  c[cn].Injectharm->pertype = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->fixperiod = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"rand",4) && strlen(argv[i]) == 4)
		{
		  c[cn].Injectharm->pertype = PERTYPE_UNIFORMRAND;
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->minp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->maxp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"logrand",7) && strlen(argv[i]) == 7)
		{
		  c[cn].Injectharm->pertype = PERTYPE_LOGRAND;
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->minp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->maxp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"randfreq",8) && strlen(argv[i]) == 8)
		{
		  c[cn].Injectharm->pertype = PERTYPE_UNIFORMRANDFREQ;
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->minf = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->maxf = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"lograndfreq",11) && strlen(argv[i]) == 11)
		{
		  c[cn].Injectharm->pertype = PERTYPE_LOGRANDFREQ;
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->minf = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->maxf = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Injectharm->Nharm = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if((c[cn].Injectharm->harm_amptype = (int *) malloc((c[cn].Injectharm->Nharm + 1) * sizeof(int))) == NULL ||
	     (c[cn].Injectharm->harm_amprel = (int *) malloc((c[cn].Injectharm->Nharm + 1) * sizeof(int))) == NULL ||
	     (c[cn].Injectharm->harm_ampfix = (double *) malloc((c[cn].Injectharm->Nharm + 1) * sizeof(double))) == NULL ||
	     (c[cn].Injectharm->harm_ampspec = (double ***) malloc((c[cn].Injectharm->Nharm + 1) * sizeof(double **))) == NULL ||
	     (c[cn].Injectharm->harm_minamp = (double *) malloc((c[cn].Injectharm->Nharm + 1) * sizeof(double))) == NULL ||
	     (c[cn].Injectharm->harm_maxamp = (double *) malloc((c[cn].Injectharm->Nharm + 1) * sizeof(double))) == NULL ||
	     (c[cn].Injectharm->harm_phasetype = (int *) malloc((c[cn].Injectharm->Nharm + 1) * sizeof(int))) == NULL ||
	     (c[cn].Injectharm->harm_phaserel = (int *) malloc((c[cn].Injectharm->Nharm + 1) * sizeof(int))) == NULL ||
	     (c[cn].Injectharm->harm_phasefix = (double *) malloc((c[cn].Injectharm->Nharm + 1) * sizeof(double))) == NULL ||
	     (c[cn].Injectharm->harm_phasespec = (double ***) malloc((c[cn].Injectharm->Nharm + 1) * sizeof(double **))) == NULL)
	    error(ERR_MEMALLOC);
	  for(l=0; l <= c[cn].Injectharm->Nharm; l++)
	    {
	      i++;
	      if(i >= argc)
		listcommands(argv[iterm],p);
	      if(!strncmp(argv[i],"amplist",7) && strlen(argv[i]) == 7)
		{
		  c[cn].Injectharm->harm_amptype[l] = PERTYPE_SPECIFIED;
		  k = 0;
		  i++;
		  if(i < argc) {
		    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
		      i++;
		      if(i < argc) {
			k = atoi(argv[i]);
		      }
		      else
			listcommands(argv[iterm],p);
		    }
		    else i--;
		  } else i--;
		  RegisterDataFromInputList(p,
					    (void *)(&(c[cn].Injectharm->harm_ampspec[l])),
					    VARTOOLS_TYPE_DOUBLE,
					    1, cn, 0, 0, NULL, k,
					    "INJECTHARM_HARM_AMP");
		}
	      else if(!strncmp(argv[i],"ampfix",6) && strlen(argv[i]) == 6)
		{
		  c[cn].Injectharm->harm_amptype[l] = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->harm_ampfix[l] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"amprand",7) && strlen(argv[i]) == 7)
		{
		  c[cn].Injectharm->harm_amptype[l] = PERTYPE_UNIFORMRAND;
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->harm_minamp[l] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->harm_maxamp[l] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"amplogrand",10) && strlen(argv[i]) == 10)
		{
		  c[cn].Injectharm->harm_amptype[l] = PERTYPE_LOGRAND;
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->harm_minamp[l] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->harm_maxamp[l] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		listcommands(argv[iterm],p);

	      c[cn].Injectharm->harm_amprel[l] = 0;
	      i++;
	      if(i < argc)
		{
		  if(!strncmp(argv[i],"amprel",6) && strlen(argv[i]) == 6)
		    {
		      if(l)
			c[cn].Injectharm->harm_amprel[l] = 1;
		    }
		  else
		    i--;
		}
	      else
		i--;

	      i++;
	      if(i >= argc)
		listcommands(argv[iterm],p);
	      if(!strncmp(argv[i],"phaselist",9) && strlen(argv[i]) == 9)
		{
		  c[cn].Injectharm->harm_phasetype[l] = PERTYPE_SPECIFIED;
		  k = 0;
		  i++;
		  if(i < argc) {
		    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
		      i++;
		      if(i < argc) {
			k = atoi(argv[i]);
		      }
		      else
			listcommands(argv[iterm],p);
		    }
		    else i--;
		  } else i--;
		  RegisterDataFromInputList(p,
					    (void *)(&(c[cn].Injectharm->harm_phasespec[l])),
					    VARTOOLS_TYPE_DOUBLE,
					    1, cn, 0, 0, NULL, k,
					    "INJECTHARM_HARM_PHASE");
		}
	      else if(!strncmp(argv[i],"phasefix",8) && strlen(argv[i]) == 8)
		{
		  c[cn].Injectharm->harm_phasetype[l] = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->harm_phasefix[l] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"phaserand",9) && strlen(argv[i]) == 9)
		{
		  c[cn].Injectharm->harm_phasetype[l] = PERTYPE_UNIFORMRAND;
		}
	      else
		listcommands(argv[iterm],p);

	      c[cn].Injectharm->harm_phaserel[l] = 0;
	      i++;
	      if(i < argc)
		{
		  if(!strncmp(argv[i],"phaserel",8) && strlen(argv[i]) == 8)
		    {
		      if(l)
			c[cn].Injectharm->harm_phaserel[l] = 1;
		    }
		  else
		    i--;
		}
	      else
		i--;
	    }
	  i++;
	  if(i < argc)
	    c[cn].Injectharm->Nsubharm = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].Injectharm->Nsubharm > 0)
	    {
	      if((c[cn].Injectharm->subharm_amptype = (int *) malloc((c[cn].Injectharm->Nsubharm) * sizeof(int))) == NULL ||
		 (c[cn].Injectharm->subharm_amprel = (int *) malloc((c[cn].Injectharm->Nsubharm) * sizeof(int))) == NULL ||
		 (c[cn].Injectharm->subharm_ampfix = (double *) malloc((c[cn].Injectharm->Nsubharm))) == NULL ||
		 (c[cn].Injectharm->subharm_ampspec = (double ***) malloc((c[cn].Injectharm->Nharm + 1) * sizeof(double **))) == NULL ||
		 (c[cn].Injectharm->subharm_minamp = (double *) malloc((c[cn].Injectharm->Nsubharm) * sizeof(double))) == NULL ||
		 (c[cn].Injectharm->subharm_maxamp = (double *) malloc((c[cn].Injectharm->Nsubharm) * sizeof(double))) == NULL ||
		 (c[cn].Injectharm->subharm_phasetype = (int *) malloc((c[cn].Injectharm->Nsubharm) * sizeof(int))) == NULL ||
		 (c[cn].Injectharm->subharm_phaserel = (int *) malloc((c[cn].Injectharm->Nsubharm) * sizeof(int))) == NULL ||
		 (c[cn].Injectharm->subharm_phasefix = (double *) malloc((c[cn].Injectharm->Nsubharm) * sizeof(double))) == NULL ||
		 (c[cn].Injectharm->subharm_phasespec = (double ***) malloc((c[cn].Injectharm->Nharm + 1) * sizeof(double **))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  for(l=0; l < c[cn].Injectharm->Nsubharm; l++)
	    {
	      i++;
	      if(i >= argc)
		listcommands(argv[iterm],p);
	      if(!strncmp(argv[i],"amplist",7) && strlen(argv[i]) == 7)
		{
		  c[cn].Injectharm->subharm_amptype[l] = PERTYPE_SPECIFIED;
		  k = 0;
		  i++;
		  if(i < argc) {
		    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
		      i++;
		      if(i < argc) {
			k = atoi(argv[i]);
		      }
		      else
			listcommands(argv[iterm],p);
		    }
		    else i--;
		  } else i--;
		  RegisterDataFromInputList(p,
					    (void *)(&(c[cn].Injectharm->subharm_ampspec[l])),
					    VARTOOLS_TYPE_DOUBLE,
					    1, cn, 0, 0, NULL, k,
					    "INJECTHARM_SUBHARM_AMP");
		}
	      else if(!strncmp(argv[i],"ampfix",6) && strlen(argv[i]) == 6)
		{
		  c[cn].Injectharm->subharm_amptype[l] = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->subharm_ampfix[l] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"amprand",7) && strlen(argv[i]) == 7)
		{
		  c[cn].Injectharm->subharm_amptype[l] = PERTYPE_UNIFORMRAND;
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->subharm_minamp[l] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->subharm_maxamp[l] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"amplogrand",10) && strlen(argv[i]) == 10)
		{
		  c[cn].Injectharm->subharm_amptype[l] = PERTYPE_LOGRAND;
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->subharm_minamp[l] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->subharm_maxamp[l] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		listcommands(argv[iterm],p);

	      c[cn].Injectharm->subharm_amprel[l] = 0;
	      i++;
	      if(i < argc)
		{
		  if(!strncmp(argv[i],"amprel",6) && strlen(argv[i]) == 6)
		    c[cn].Injectharm->subharm_amprel[l] = 1;
		  else
		    i--;
		}
	      else
		i--;

	      i++;
	      if(i >= argc)
		listcommands(argv[iterm],p);
	      if(!strncmp(argv[i],"phaselist",9) && strlen(argv[i]) == 9)
		{
		  c[cn].Injectharm->subharm_phasetype[l] = PERTYPE_SPECIFIED;
		  k = 0;
		  i++;
		  if(i < argc) {
		    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
		      i++;
		      if(i < argc) {
			k = atoi(argv[i]);
		      }
		      else
			listcommands(argv[iterm],p);
		    }
		    else i--;
		  } else i--;
		  RegisterDataFromInputList(p,
					    (void *)(&(c[cn].Injectharm->subharm_phasespec[l])),
					    VARTOOLS_TYPE_DOUBLE,
					    1, cn, 0, 0, NULL, k,
					    "INJECTHARM_SUBHARM_PHASE");
		}
	      else if(!strncmp(argv[i],"phasefix",8) && strlen(argv[i]) == 8)
		{
		  c[cn].Injectharm->subharm_phasetype[l] = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    c[cn].Injectharm->subharm_phasefix[l] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"phaserand",9) && strlen(argv[i]) == 9)
		{
		  c[cn].Injectharm->subharm_phasetype[l] = PERTYPE_UNIFORMRAND;
		}
	      else
		listcommands(argv[iterm],p);

	      c[cn].Injectharm->subharm_phaserel[l] = 0;
	      i++;
	      if(i < argc)
		{
		  if(!strncmp(argv[i],"phaserel",8) && strlen(argv[i]) == 8)
		    c[cn].Injectharm->subharm_phaserel[l] = 1;
		  else
		    i--;
		}
	      else
		i--;
	    }
	  i++;
	  if(i < argc)
	    c[cn].Injectharm->omodel = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].Injectharm->omodel)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].Injectharm->modeloutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].Injectharm->modelsuffix,".injectharm.model");
	    }
	  cn++;
	}

      /* -Injecttransit <\"Plist\" [\"column\" col] | \"Pfix\" per | \"Prand\" minp maxp | \"Plogrand\" minp maxp | \"randfreq\" minf maxf | \"lograndfreq\" minf maxf> <\"Rplist\" [\"column\" col] | \"Rpfix\" Rp | \"Rprand\" minRp maxRp | \"Rplogrand\" minRp maxRp> <\"Mplist\" [\"column\" col] | \"Mpfix\" Mp | \"Mprand\" minMp maxMp | \"Mplogrand\" minMp maxMp> <\"phaselist\" [\"column\" col] | \"phasefix\" phase | \"phaserand> <\"sinilist\" [\"column\" col] | \"sinifix\" sin_i | \"sinirand\"> <\"eomega <\"elist\" [\"column\" col] | \"efix\" e | \"erand\"> <\"olist\" [\"column\" col] | \"ofix\" omega | \"orand\"> | \"hk\" <\"hlist\" [\"column\" col] | \"hfix\" h | \"hrand\"> <\"klist\" [\"column\" col] | \"kfix\" k | \"krand\">> <\"Mstarlist\" [\"column\" col] | \"Mstarfix\" Mstar> <\"Rstarlist\" [\"column\" col] | \"Rstarfix\" Rstar> <\"quad\" | \"nonlin\"> <\"ldlist\" [\"column\" col] | \"ldfix\" ld1 ... ldn> [\"dilute\" <\"list\" [\"column\" col] | \"fix\" dilute>] omodel [modeloutdir] */
      else if(!strncmp(argv[i],"-Injecttransit",14) && strlen(argv[i]) == 14)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_INJECTTRANSIT;
	  if((c[cn].Injecttransit = (_Injecttransit *) malloc(sizeof(_Injecttransit))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"Plist",5) && strlen(argv[i]) == 5) {
		c[cn].Injecttransit->paramtype[INJECTTR_IDX_PERIOD] = PERTYPE_SPECIFIED;
		k = 0;
		i++;
		if(i < argc) {
		  if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
		    {
		      i++;
		      if(i < argc)
			k = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		    } else i--;
		} else i--;
		RegisterDataFromInputList(p,
					  (void *) (&(c[cn].Injecttransit->paramspec[INJECTTR_IDX_PERIOD])),
					  VARTOOLS_TYPE_DOUBLE,
					  1, cn, 0, 0, NULL, k,
					  "INJECTTRANSIT_PERIOD");
	      }
	      else if(!strncmp(argv[i],"Pfix",4) && strlen(argv[i]) == 4)
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_PERIOD] = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->paramfix[INJECTTR_IDX_PERIOD] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"Prand",5) && strlen(argv[i]) == 5)
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_PERIOD] = PERTYPE_UNIFORMRAND;
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->minp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->maxp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"Plogrand",8) && strlen(argv[i]) == 8)
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_PERIOD] = PERTYPE_LOGRAND;
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->minp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->maxp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strcmp(argv[i],"Pexpr"))
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_PERIOD] = PERTYPE_EXPR;
		  i++;
		  if(i >= argc) 
		    listcommands(argv[iterm],p);
		  parse_setparam_expr(&(c[cn]), argv[i], &(c[cn].Injecttransit->paramexpr[INJECTTR_IDX_PERIOD]));
		}
	      else if(!strncmp(argv[i],"randfreq",8) && strlen(argv[i]) == 8)
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_PERIOD] = PERTYPE_UNIFORMRANDFREQ;
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->minf = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->maxf = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"lograndfreq",11) && strlen(argv[i]) == 11)
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_PERIOD] = PERTYPE_LOGRANDFREQ;
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->minf = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->maxf = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"Rplist",6) && strlen(argv[i]) == 6) {
		c[cn].Injecttransit->paramtype[INJECTTR_IDX_RP] = PERTYPE_SPECIFIED;
		k = 0;
		i++;
		if(i < argc) {
		  if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
		    {
		      i++;
		      if(i < argc)
			k = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		    } else i--;
		} else i--;
		RegisterDataFromInputList(p,
					  (void *) (&(c[cn].Injecttransit->paramspec[INJECTTR_IDX_RP])),
					  VARTOOLS_TYPE_DOUBLE,
					  1, cn, 0, 0, NULL, k,
					  "INJECTTRANSIT_RADIUS_PLANET");
	      }
	      else if(!strncmp(argv[i],"Rpfix",5) && strlen(argv[i]) == 5)
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_RP] = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->paramfix[INJECTTR_IDX_RP] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"Rprand",6) && strlen(argv[i]) == 6)
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_RP] = PERTYPE_UNIFORMRAND;
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->minRp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->maxRp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"Rplogrand",9) && strlen(argv[i]) == 9)
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_RP] = PERTYPE_LOGRAND;
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->minRp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->maxRp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strcmp(argv[i],"Rpexpr"))
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_RP] = PERTYPE_EXPR;
		  i++;
		  if(i >= argc) 
		    listcommands(argv[iterm],p);
		  parse_setparam_expr(&(c[cn]), argv[i], &(c[cn].Injecttransit->paramexpr[INJECTTR_IDX_RP]));
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"Mplist",6) && strlen(argv[i]) == 6) {
		c[cn].Injecttransit->paramtype[INJECTTR_IDX_MP] = PERTYPE_SPECIFIED;
		k = 0;
		i++;
		if(i < argc) {
		  if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
		    {
		      i++;
		      if(i < argc)
			k = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		    } else i--;
		} else i--;
		RegisterDataFromInputList(p,
					  (void *) (&(c[cn].Injecttransit->paramspec[INJECTTR_IDX_MP])),
					  VARTOOLS_TYPE_DOUBLE,
					  1, cn, 0, 0, NULL, k,
					  "INJECTTRANSIT_MASS_PLANET");
	      }
	      else if(!strncmp(argv[i],"Mpfix",5) && strlen(argv[i]) == 5)
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_MP] = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->paramfix[INJECTTR_IDX_MP] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"Mprand",6) && strlen(argv[i]) == 6)
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_MP] = PERTYPE_UNIFORMRAND;
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->minMp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->maxMp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"Mplogrand",9) && strlen(argv[i]) == 9)
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_MP] = PERTYPE_LOGRAND;
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->minMp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->maxMp = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strcmp(argv[i],"Mpexpr"))
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_MP] = PERTYPE_EXPR;
		  i++;
		  if(i >= argc) 
		    listcommands(argv[iterm],p);
		  parse_setparam_expr(&(c[cn]), argv[i], &(c[cn].Injecttransit->paramexpr[INJECTTR_IDX_MP]));
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"phaselist",9) && strlen(argv[i]) == 9) {
		c[cn].Injecttransit->paramtype[INJECTTR_IDX_PHASE] = PERTYPE_SPECIFIED;
		k = 0;
		i++;
		if(i < argc) {
		  if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
		    {
		      i++;
		      if(i < argc)
			k = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		    } else i--;
		} else i--;
		RegisterDataFromInputList(p,
					  (void *) (&(c[cn].Injecttransit->paramspec[INJECTTR_IDX_PHASE])),
					  VARTOOLS_TYPE_DOUBLE,
					  1, cn, 0, 0, NULL, k,
					  "INJECTTRANSIT_PHASE");
	      }
	      else if(!strncmp(argv[i],"phasefix",8) && strlen(argv[i]) == 8)
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_PHASE] = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->paramfix[INJECTTR_IDX_PHASE] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"phaserand",9) && strlen(argv[i]) == 9)
		c[cn].Injecttransit->paramtype[INJECTTR_IDX_PHASE] = PERTYPE_UNIFORMRAND;
	      else if(!strcmp(argv[i],"phaseexpr"))
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_PHASE] = PERTYPE_EXPR;
		  i++;
		  if(i >= argc) 
		    listcommands(argv[iterm],p);
		  parse_setparam_expr(&(c[cn]), argv[i], &(c[cn].Injecttransit->paramexpr[INJECTTR_IDX_PHASE]));
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"sinilist",8) && strlen(argv[i]) == 8) {
		c[cn].Injecttransit->paramtype[INJECTTR_IDX_SINI] = PERTYPE_SPECIFIED;
		k = 0;
		i++;
		if(i < argc) {
		  if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
		    {
		      i++;
		      if(i < argc)
			k = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		    } else i--;
		} else i--;
		RegisterDataFromInputList(p,
					  (void *) (&(c[cn].Injecttransit->paramspec[INJECTTR_IDX_SINI])),
					  VARTOOLS_TYPE_DOUBLE,
					  1, cn, 0, 0, NULL, k,
					  "INJECTTRANSIT_SINI");
	      }
	      else if(!strncmp(argv[i],"sinifix",7) && strlen(argv[i]) == 7)
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_SINI] = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->paramfix[INJECTTR_IDX_SINI] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"sinirand",8) && strlen(argv[i]) == 8)
		c[cn].Injecttransit->paramtype[INJECTTR_IDX_SINI] = PERTYPE_UNIFORMRAND;
	      else if(!strcmp(argv[i],"siniexpr"))
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_SINI] = PERTYPE_EXPR;
		  i++;
		  if(i >= argc) 
		    listcommands(argv[iterm],p);
		  parse_setparam_expr(&(c[cn]), argv[i], &(c[cn].Injecttransit->paramexpr[INJECTTR_IDX_SINI]));
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"eomega",6) && strlen(argv[i]) == 6)
		{
		  c[cn].Injecttransit->eomegatype = 0;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"elist",5) && strlen(argv[i]) == 5) {
			c[cn].Injecttransit->paramtype[INJECTTR_IDX_E] = PERTYPE_SPECIFIED;
			k = 0;
			i++;
			if(i < argc) {
			  if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
			    {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    } else i--;
			} else i--;
			RegisterDataFromInputList(p,
						  (void *) (&(c[cn].Injecttransit->paramspec[INJECTTR_IDX_E])),
						  VARTOOLS_TYPE_DOUBLE,
						  1, cn, 0, 0, NULL, k,
						  "INJECTTRANSIT_ECCEN");
		      }
		      else if(!strncmp(argv[i],"efix",4) && strlen(argv[i]) == 4)
			{
			  c[cn].Injecttransit->paramtype[INJECTTR_IDX_E] = PERTYPE_FIX;
			  i++;
			  if(i < argc)
			    c[cn].Injecttransit->paramfix[INJECTTR_IDX_E] = atof(argv[i]);
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"erand",5) && strlen(argv[i]) == 5)
			c[cn].Injecttransit->paramtype[INJECTTR_IDX_E] = PERTYPE_UNIFORMRAND;
		      else if(!strcmp(argv[i],"eexpr"))
			{
			  c[cn].Injecttransit->paramtype[INJECTTR_IDX_E] = PERTYPE_EXPR;
			  i++;
			  if(i >= argc) 
			    listcommands(argv[iterm],p);
			  parse_setparam_expr(&(c[cn]), argv[i], &(c[cn].Injecttransit->paramexpr[INJECTTR_IDX_E]));
			}
		      else
			listcommands(argv[iterm],p);
		    }
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"olist",5) && strlen(argv[i]) == 5) {
			c[cn].Injecttransit->paramtype[INJECTTR_IDX_OMEGA] = PERTYPE_SPECIFIED;
			k = 0;
			i++;
			if(i < argc) {
			  if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
			    {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    } else i--;
			} else i--;
			RegisterDataFromInputList(p,
						  (void *) (&(c[cn].Injecttransit->paramspec[INJECTTR_IDX_OMEGA])),
						  VARTOOLS_TYPE_DOUBLE,
						  1, cn, 0, 0, NULL, k,
						  "INJECTTRANSIT_OMEGA");
		      }
		      else if(!strncmp(argv[i],"ofix",4) && strlen(argv[i]) == 4)
			{
			  c[cn].Injecttransit->paramtype[INJECTTR_IDX_OMEGA] = PERTYPE_FIX;
			  i++;
			  if(i < argc)
			    c[cn].Injecttransit->paramfix[INJECTTR_IDX_OMEGA] = atof(argv[i]);
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"orand",9) && strlen(argv[i]) == 5)
			c[cn].Injecttransit->paramtype[INJECTTR_IDX_OMEGA] = PERTYPE_UNIFORMRAND;
		      else if(!strcmp(argv[i],"oexpr"))
			{
			  c[cn].Injecttransit->paramtype[INJECTTR_IDX_OMEGA] = PERTYPE_EXPR;
			  i++;
			  if(i >= argc) 
			    listcommands(argv[iterm],p);
			  parse_setparam_expr(&(c[cn]), argv[i], &(c[cn].Injecttransit->paramexpr[INJECTTR_IDX_OMEGA]));
			}
		      else
			listcommands(argv[iterm],p);
		    }
		  else
		    listcommands(argv[iterm],p);
		  
		}
	      else if(!strncmp(argv[i],"hk",2) && strlen(argv[i]) == 2)
		{
		  c[cn].Injecttransit->eomegatype = 1;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"hlist",5) && strlen(argv[i]) == 5) {
			c[cn].Injecttransit->paramtype[INJECTTR_IDX_H] = PERTYPE_SPECIFIED;
			k = 0;
			i++;
			if(i < argc) {
			  if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
			    {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    } else i--;
			} else i--;
			RegisterDataFromInputList(p,
						  (void *) (&(c[cn].Injecttransit->paramspec[INJECTTR_IDX_H])),
						  VARTOOLS_TYPE_DOUBLE,
						  1, cn, 0, 0, NULL, k,
						  "INJECTTRANSIT_H");
		      }
		      else if(!strncmp(argv[i],"hfix",4) && strlen(argv[i]) == 4)
			{
			  c[cn].Injecttransit->paramtype[INJECTTR_IDX_H] = PERTYPE_FIX;
			  i++;
			  if(i < argc)
			    c[cn].Injecttransit->paramfix[INJECTTR_IDX_H] = atof(argv[i]);
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"hrand",5) && strlen(argv[i]) == 5)
			c[cn].Injecttransit->paramtype[INJECTTR_IDX_H] = PERTYPE_UNIFORMRAND;
		      else if(!strcmp(argv[i],"hexpr"))
			{
			  c[cn].Injecttransit->paramtype[INJECTTR_IDX_H] = PERTYPE_EXPR;
			  i++;
			  if(i >= argc) 
			    listcommands(argv[iterm],p);
			  parse_setparam_expr(&(c[cn]), argv[i], &(c[cn].Injecttransit->paramexpr[INJECTTR_IDX_H]));
			}
		      else
			listcommands(argv[iterm],p);
		    }
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"klist",5) && strlen(argv[i]) == 5) {
			c[cn].Injecttransit->paramtype[INJECTTR_IDX_K] = PERTYPE_SPECIFIED;
			k = 0;
			i++;
			if(i < argc) {
			  if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
			    {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    } else i--;
			} else i--;
			RegisterDataFromInputList(p,
						  (void *) (&(c[cn].Injecttransit->paramspec[INJECTTR_IDX_K])),
						  VARTOOLS_TYPE_DOUBLE,
						  1, cn, 0, 0, NULL, k,
						  "INJECTTRANSIT_K");
		      }
		      else if(!strncmp(argv[i],"kfix",4) && strlen(argv[i]) == 4)
			{
			  c[cn].Injecttransit->paramtype[INJECTTR_IDX_K] = PERTYPE_FIX;
			  i++;
			  if(i < argc)
			    c[cn].Injecttransit->paramfix[INJECTTR_IDX_K] = atof(argv[i]);
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"krand",5) && strlen(argv[i]) == 5)
			c[cn].Injecttransit->paramtype[INJECTTR_IDX_K] = PERTYPE_UNIFORMRAND;
		      else if(!strcmp(argv[i],"kexpr"))
			{
			  c[cn].Injecttransit->paramtype[INJECTTR_IDX_K] = PERTYPE_EXPR;
			  i++;
			  if(i >= argc) 
			    listcommands(argv[iterm],p);
			  parse_setparam_expr(&(c[cn]), argv[i], &(c[cn].Injecttransit->paramexpr[INJECTTR_IDX_K]));
			}
		      else
			listcommands(argv[iterm],p);
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"Mstarlist",9) && strlen(argv[i]) == 9) {
		c[cn].Injecttransit->paramtype[INJECTTR_IDX_MSTAR] = PERTYPE_SPECIFIED;
		k = 0;
		i++;
		if(i < argc) {
		  if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
		    {
		      i++;
		      if(i < argc)
			k = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		    } else i--;
		} else i--;
		RegisterDataFromInputList(p,
					  (void *) (&(c[cn].Injecttransit->paramspec[INJECTTR_IDX_MSTAR])),
					  VARTOOLS_TYPE_DOUBLE,
					  1, cn, 0, 0, NULL, k,
					  "INJECTTRANSIT_MASS_STAR");
	      }
	      else if(!strncmp(argv[i],"Mstarfix",8) && strlen(argv[i]) == 8)
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_MSTAR] = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->paramfix[INJECTTR_IDX_MSTAR] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strcmp(argv[i],"Mstarexpr"))
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_MSTAR] = PERTYPE_EXPR;
		  i++;
		  if(i >= argc) 
		    listcommands(argv[iterm],p);
		  parse_setparam_expr(&(c[cn]), argv[i], &(c[cn].Injecttransit->paramexpr[INJECTTR_IDX_MSTAR]));
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"Rstarlist",9) && strlen(argv[i]) == 9) {
		c[cn].Injecttransit->paramtype[INJECTTR_IDX_RSTAR] = PERTYPE_SPECIFIED;
		k = 0;
		i++;
		if(i < argc) {
		  if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
		    {
		      i++;
		      if(i < argc)
			k = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		    } else i--;
		} else i--;
		RegisterDataFromInputList(p,
					  (void *) (&(c[cn].Injecttransit->paramspec[INJECTTR_IDX_RSTAR])),
					  VARTOOLS_TYPE_DOUBLE,
					  1, cn, 0, 0, NULL, k,
					  "INJECTTRANSIT_RADIUS_STAR");
	      }
	      else if(!strncmp(argv[i],"Rstarfix",8) && strlen(argv[i]) == 8)
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_RSTAR] = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    c[cn].Injecttransit->paramfix[INJECTTR_IDX_RSTAR] = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strcmp(argv[i],"Rstarexpr"))
		{
		  c[cn].Injecttransit->paramtype[INJECTTR_IDX_RSTAR] = PERTYPE_EXPR;
		  i++;
		  if(i >= argc) 
		    listcommands(argv[iterm],p);
		  parse_setparam_expr(&(c[cn]), argv[i], &(c[cn].Injecttransit->paramexpr[INJECTTR_IDX_RSTAR]));
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"quad",4) && strlen(argv[i]) == 4)
		{
		  c[cn].Injecttransit->ldtype = 0;
		  c[cn].Injecttransit->Nld = 2;
		}
	      else if(!strncmp(argv[i],"nonlin",6) && strlen(argv[i]) == 6)
		{
		  c[cn].Injecttransit->ldtype = 1;
		  c[cn].Injecttransit->Nld = 4;
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"ldlist",6) && strlen(argv[i]) == 6)
		{
		  k = 0;
		  i++;
		  if(i < argc) {
		    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
		      {
			i++;
			if(i < argc)
			  k = atoi(argv[i]);
			else
			  listcommands(argv[iterm],p);
		      } else i--;
		  } else i--;
		  for(j=0;j<c[cn].Injecttransit->Nld;j++) {
		    c[cn].Injecttransit->paramtype[INJECTTR_IDX_LD + j] = PERTYPE_SPECIFIED;
		    RegisterDataFromInputList(p,
					      (void *) (&(c[cn].Injecttransit->paramspec[INJECTTR_IDX_LD + j])),
					      VARTOOLS_TYPE_DOUBLE,
					      1, cn, 0, 0, NULL, k,
					      "INJECTTRANSIT_LD");
		    if(k > 0) k++;
		  }
		}
	      else if(!strncmp(argv[i],"ldfix",5) && strlen(argv[i]) == 5)
		{
		  for(j=0;j<c[cn].Injecttransit->Nld;j++)
		    {
		      c[cn].Injecttransit->paramtype[INJECTTR_IDX_LD + j] = PERTYPE_FIX;
		      i++;
		      if(i < argc)
			c[cn].Injecttransit->paramfix[INJECTTR_IDX_LD + j] = atof(argv[i]);
		      else
			listcommands(argv[iterm],p);
		    }
		}
	      else if(!strcmp(argv[i],"ldexpr"))
		{
		  for(j=0;j<c[cn].Injecttransit->Nld;j++)
		    {
		      c[cn].Injecttransit->paramtype[INJECTTR_IDX_LD + j] = PERTYPE_EXPR;
		      i++;
		      if(i >= argc) 
			listcommands(argv[iterm],p);
		      parse_setparam_expr(&(c[cn]), argv[i], &(c[cn].Injecttransit->paramexpr[INJECTTR_IDX_LD + j]));
		    }
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  c[cn].Injecttransit->Nparam = INJECTTR_IDX_LD + c[cn].Injecttransit->Nld;
	  i++;
	  if(i < argc) {
	    if(!strncmp(argv[i],"dilute", 6) && strlen(argv[i]) == 6)
	      {
		i++;
		if(i < argc) {
		  if(!strncmp(argv[i],"list", 4) && strlen(argv[i]) == 4) {
		    c[cn].Injecttransit->paramtype[INJECTTR_IDX_DILUTE] = PERTYPE_SPECIFIED;
		    k = 0;
		    i++;
		    if(i < argc) {
		      if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
			{
			  i++;
			  if(i < argc)
			    k = atoi(argv[i]);
			  else
			    listcommands(argv[iterm],p);
			} else i--;
		    } else i--;
		    RegisterDataFromInputList(p,
					      (void *) (&(c[cn].Injecttransit->paramspec[INJECTTR_IDX_DILUTE])),
					      VARTOOLS_TYPE_DOUBLE,
					      1, cn, 0, 0, NULL, k,
					      "INJECTTRANSIT_DILUTE");

		  }
		  else if(!strncmp(argv[i],"fix", 3) && strlen(argv[i]) == 3) {
		    c[cn].Injecttransit->paramtype[INJECTTR_IDX_DILUTE] = PERTYPE_FIX;
		    i++;
		    if(i < argc)
		      c[cn].Injecttransit->paramfix[INJECTTR_IDX_DILUTE] = atof(argv[i]);
		    else
		      listcommands(argv[iterm],p);
		  }
		  else if(!strcmp(argv[i],"expr"))
		    {
		      c[cn].Injecttransit->paramtype[INJECTTR_IDX_DILUTE] = PERTYPE_EXPR;
		      i++;
		      if(i >= argc) 
			listcommands(argv[iterm],p);
		      parse_setparam_expr(&(c[cn]), argv[i], &(c[cn].Injecttransit->paramexpr[INJECTTR_IDX_DILUTE]));
		    }
		  else
		    listcommands(argv[iterm],p);
		}
		else
		  listcommands(argv[iterm],p);
	      }
	    else 
	      {
		i--;
		c[cn].Injecttransit->paramtype[INJECTTR_IDX_DILUTE] = PERTYPE_FIX;
		c[cn].Injecttransit->paramfix[INJECTTR_IDX_DILUTE] = 1.0;
	      }
	  }
	  else {
	    i--;
		c[cn].Injecttransit->paramtype[INJECTTR_IDX_DILUTE] = PERTYPE_FIX;
		c[cn].Injecttransit->paramfix[INJECTTR_IDX_DILUTE] = 1.0;
	  }
	  i++;
	  if(i < argc)
	    c[cn].Injecttransit->omodel = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].Injecttransit->omodel)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].Injecttransit->modeloutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].Injecttransit->modelsuffix,".injecttransit.model");
	    }
	  cn++;
	}

      /* -Killharm <\"aov\" | \"ls\" | \"both\" | \"injectharm\" | \"fix\" Nper per1 ... perN | \"list\" Nper [\"column\" col1]> Nharm Nsubharm omodel [modeloutdir] [\"fitonly\"] [\"outampphase\" | \"outampradphase\" | \"outRphi\" | \"outRradphi\"] [\"clip\" val]*/
      else if(!strncmp(argv[i],"-Killharm",9) && strlen(argv[i]) == 9)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_KILLHARM;
	  if((c[cn].Killharm = (_Killharm *) malloc(sizeof(_Killharm))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"aov",3) && strlen(argv[i]) == 3)
		{
		  c[cn].Killharm->pertype = PERTYPE_AOV;
		  c[cn].Killharm->Nper = 1;
		  m = -1;
		  /* Get the index for the last aov command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_AOV || c[l].cnum == CNUM_HARMAOV)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOAOV);
		  else
		    c[cn].Killharm->lastaovindex = m;
		}
	      else if(!strncmp(argv[i],"ls",2) && strlen(argv[i]) == 2)
		{
		  c[cn].Killharm->pertype = PERTYPE_LS;
		  c[cn].Killharm->Nper = 1;
		  m = -1;
		  /* Get the index for the last ls command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_LS)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOLS);
		  else
		    c[cn].Killharm->lastlsindex = m;
		}
	      else if(!strncmp(argv[i],"both",4) && strlen(argv[i]) == 4)
		{
		  c[cn].Killharm->pertype = PERTYPE_BOTH;
		  c[cn].Killharm->Nper = 2;
		  m = -1;
		  /* Get the index for the last aov command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_AOV || c[l].cnum == CNUM_HARMAOV)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOBOTH);
		  else
		    c[cn].Killharm->lastaovindex = m;
		  m = -1;
		  /* Get the index for the last ls command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_LS)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOBOTH);
		  else
		    c[cn].Killharm->lastlsindex = m;
		}
	      else if(!strncmp(argv[i],"injectharm",10) && strlen(argv[i]) == 10)
		{
		  c[cn].Killharm->pertype = PERTYPE_INJECTHARM;
		  c[cn].Killharm->Nper = 1;
		  m = -1;
		  /* Get the index for the last injectharm command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_INJECTHARM)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOINJECTHARM);
		  else
		    c[cn].Killharm->lastaovindex = m;
		}
	      else if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
		{
		  c[cn].Killharm->pertype = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    c[cn].Killharm->Nper = atoi(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  if(c[cn].Killharm->Nper <= 0)
		    error(ERR_KILLHARM_WRONGNPER);
		  if((c[cn].Killharm->fixedperiods = (double *) malloc(c[cn].Killharm->Nper * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		  for(j=0;j<c[cn].Killharm->Nper;j++)
		    {
		      i++;
		      if(i < argc)
			c[cn].Killharm->fixedperiods[j] = atof(argv[i]);
		      else
			listcommands(argv[iterm],p);
		      if(c[cn].Killharm->fixedperiods[j] <= 0)
			error(ERR_KILLHARM_NEGATIVEPERIOD);
		    }
		}
	      else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
		{
		  c[cn].Killharm->pertype = PERTYPE_SPECIFIED;
		  i++;
		  if(i < argc)
		    c[cn].Killharm->Nper = atoi(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  if(c[cn].Killharm->Nper <= 0)
		    error(ERR_KILLHARM_WRONGNPER);
		  k = 0;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6)
			{
			  i++;
			  if(i < argc)
			    k = atoi(argv[i]);
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			i--;
		    }
		  else
		    i--;
		  RegisterDataFromInputList(p, 
					    (void *) &(c[cn].Killharm->periods),
					    VARTOOLS_TYPE_DOUBLE,
					    c[cn].Killharm->Nper,
					    cn, 0, 0, NULL, k, "Killharm_Input_Period");
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Killharm->Nharm = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Killharm->Nsubharm = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Killharm->omodel = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].Killharm->omodel)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].Killharm->modeloutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].Killharm->modelsuffix,".killharm.model");
	    }
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"fitonly",7) && strlen(argv[i]) == 7)
		{
		  c[cn].Killharm->fitonly = 1;
		}
	      else
		{
		  c[cn].Killharm->fitonly = 0;
		  i--;
		}
	    }
	  else
	    {
	      c[cn].Killharm->fitonly = 0;
	      i--;
	    }
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"outampphase",11) && strlen(argv[i]) == 11)
		{
		  c[cn].Killharm->outtype = KILLHARM_OUTTYPE_AMPPHASE;
		}
	      else if(!strncmp(argv[i],"outampradphase",14) && strlen(argv[i]) == 14)
		{
		  c[cn].Killharm->outtype = KILLHARM_OUTTYPE_AMPRADPHASE;
		}
	      else if(!strncmp(argv[i],"outRphi",7) && strlen(argv[i]) == 7)
		{
		  c[cn].Killharm->outtype = KILLHARM_OUTTYPE_RPHI;
		}
	      else if(!strncmp(argv[i],"outRradphi",10) && strlen(argv[i]) == 10)
		{
		  c[cn].Killharm->outtype = KILLHARM_OUTTYPE_RRADPHI;
		}
	      else
		{
		  c[cn].Killharm->outtype = KILLHARM_OUTTYPE_DEFAULT;
		  i--;
		}
	    }
	  else
	    {
	      c[cn].Killharm->outtype = KILLHARM_OUTTYPE_DEFAULT;
	      i--;
	    }

	  c[cn].Killharm->clip = -1.;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"clip")) {
	      i++;
	      if(i >= argc)
		listcommands(argv[iterm],p);
	      c[cn].Killharm->clip = atof(argv[i]);
	    } else {
	      i--;
	    }
	  } else {
	    i--;
	  }
	  cn++;
	}

      /* -Starspot <"aov" | "ls" | "list" ["column" col] | "fix" period | "fixcolumn" <colname | colnum>> a0 b0 alpha0 inc0 chi0 psi00 mconst0 correctlc */
      else if(!strncmp(argv[i],"-Starspot",9) && strlen(argv[i]) == 9)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_STARSPOT;
	  if((c[cn].Starspot = (_Starspot *) malloc(sizeof(_Starspot))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"aov",3) && strlen(argv[i]) == 3)
		{
		  c[cn].Starspot->pertype = PERTYPE_AOV;
		  m = -1;
		  /* Get the index for the last aov command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_AOV || c[l].cnum == CNUM_HARMAOV)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOAOV);
		  else
		    c[cn].Starspot->lastaovindex = m;
		}
	      else if(!strncmp(argv[i],"ls",2) && strlen(argv[i]) == 2)
		{
		  c[cn].Starspot->pertype = PERTYPE_LS;
		  m = -1;
		  /* Get the index for the last ls command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_LS)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOLS);
		  else
		    c[cn].Starspot->lastlsindex = m;
		}
	      else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
		{
		  c[cn].Starspot->pertype = PERTYPE_SPECIFIED;
		  k = 0;
		  i++;
		  if(i < argc) {
		    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
		      i++;
		      if(i < argc)
			k = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		    } else i--;
		  } else i--;
		  RegisterDataFromInputList(p,
					    (void *) (&(c[cn].Starspot->period)),
					    VARTOOLS_TYPE_DOUBLE,
					    1, cn, 0, 0, NULL, k,
					    "STARSPOT_PERIOD");
		}
	      else if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
		{
		  c[cn].Starspot->pertype = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    {
		      c[cn].Starspot->fixedperiod = atof(argv[i]);
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"fixcolumn",9) && strlen(argv[i]) == 9)
		{
		  c[cn].Starspot->pertype = PERTYPE_FIXCOLUMN;
		  i++;
		  if(i < argc)
		    increaselinkedcols(p, &(c[cn].Starspot->linkedcolumn), argv[i], cn);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->a0 = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->b0 = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->alpha0 = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->inclination0 = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->chi0 = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->psi00 = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->mconst0 = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->fitP = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->fita = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->fitb = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->fitalpha = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->fiti = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->fitchi = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->fitpsi0 = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->fitmconst = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->correctlc = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Starspot->omodel = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].Starspot->omodel)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].Starspot->modeloutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].Starspot->modelsuffix,".starspot.model");
	    }
	  cn++;
	} 

      /* -stats var1,var2,... stats1,stats2,... */
      else if(!strcmp(argv[i],"-stats"))
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_STATS;
	  if((c[cn].Stats = (_Stats *) malloc(sizeof(_Stats))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i >= argc)
	    listcommands(argv[iterm],p);
	  if(ParseStatsCommand(&i, argc, argv, p, c[cn].Stats))
	    listcommands(argv[iterm],p);
	  cn++;
	}

      /* -BLS < \"r\" rmin rmax | \"q\" qmin qmax | \"density\" rho min_expected_duration_frac max_expected_duration_frac > minper maxper nfreq nbins timezone Npeak outperiodogram [outdir] omodel [modeloutdir] correctlc [\"fittrap\"] [\"nobinnedrms\"] [\"ophcurve\" phmin phmax phstep] [\"ojdcurve\" jdstep] [\"stepP\" | \"steplogP\"] [\"adjust-qmin-by-mindt\" [\"reduce-nbins\"]] [\"reportharmonics\"]*/
      else if(!strncmp(argv[i],"-BLS",4) && strlen(argv[i]) == 4)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].require_sort = 1;
	  c[cn].require_distinct = 1;
	  c[cn].cnum = CNUM_BLS;
	  if((c[cn].Bls = (_Bls *) malloc(sizeof(_Bls))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"r",1) && strlen(argv[i]) == 1)
		{
		  c[cn].Bls->rflag = 1;
		  i++;
		  if(i < argc)
		    c[cn].Bls->rmin = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Bls->rmax = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
      	      else if(!strncmp(argv[i],"q",1) && strlen(argv[i]) == 1)
		{
		  c[cn].Bls->rflag = 0;
		  i++;
		  if(i < argc)
		    c[cn].Bls->qmin = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Bls->qmax = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strcmp(argv[i],"density"))
		{
		  c[cn].Bls->rflag = 1;
		  i++;
		  if(i < argc)
		    c[cn].Bls->rho = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Bls->rmin = pow(((0.0848203*atof(argv[i])*pow(c[cn].Bls->rho,(-1.0/3.0)))/0.076),1.5);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].Bls->rmax = pow(((0.0848203*atof(argv[i])*pow(c[cn].Bls->rho,(-1.0/3.0)))/0.076),1.5);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  i++;
	  if(i < argc)
	    c[cn].Bls->minper = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Bls->maxper = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Bls->nf = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
#ifndef PARALLEL
	  if((c[cn].Bls->p = (double *) malloc((c[cn].Bls->nf+1)*sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
#endif
	  c[cn].Bls->df = ((1./c[cn].Bls->minper) - (1./c[cn].Bls->maxper)) / (c[cn].Bls->nf - 1);
	  i++;
	  if(i < argc)
	    c[cn].Bls->nbins = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Bls->timezone = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Bls->Npeak = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Bls->operiodogram = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].Bls->operiodogram)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].Bls->outdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].Bls->suffix,".bls");
	    }
	  i++;
	  if(i < argc)
	    c[cn].Bls->omodel = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].Bls->omodel)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].Bls->modeloutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].Bls->modelsuffix,".bls.model");
	    }
	  i++;
	  if(i < argc)
	    c[cn].Bls->correctlc = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  c[cn].Bls->fittrap = 0;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"fittrap")) {
	      c[cn].Bls->fittrap = 1;
	    }
	    else
	      i--;
	  }
	  else
	    i--;
	  c[cn].Bls->nobinnedrms = 0;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"nobinnedrms")) {
	      c[cn].Bls->nobinnedrms = 1;
	    }
	    else
	      i--;
	  }
	  else
	    i--;
	  c[cn].Bls->ophcurve = 0;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"ophcurve")) {
	      c[cn].Bls->ophcurve = 1;
	      i++;
	      if(i < argc)
		sprintf(c[cn].Bls->ophcurveoutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].Bls->ophcurvesuffix,".bls.phcurve");
	      i++;
	      if(i < argc)
		c[cn].Bls->phmin = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		c[cn].Bls->phmax = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		c[cn].Bls->phstep = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	    }
	    else
	      i--;
	  }
	  else
	    i--;
	  c[cn].Bls->ojdcurve = 0;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"ojdcurve")) {
	      c[cn].Bls->ojdcurve = 1;
	      i++;
	      if(i < argc)
		sprintf(c[cn].Bls->ojdcurveoutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].Bls->ojdcurvesuffix,".bls.jdcurve");
	      i++;
	      if(i < argc)
		c[cn].Bls->jdstep = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	    }
	    else
	      i--;
	  }
	  else
	    i--;
	  c[cn].Bls->freqsteptype = VARTOOLS_FREQSTEPTYPE_FREQ;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"stepP")) {
	      c[cn].Bls->freqsteptype = VARTOOLS_FREQSTEPTYPE_PERIOD;
	      c[cn].Bls->df = (c[cn].Bls->maxper - c[cn].Bls->minper)/ (c[cn].Bls->nf - 1);
	    }
	    else if(!strcmp(argv[i],"steplogP")) {
	      c[cn].Bls->freqsteptype = VARTOOLS_FREQSTEPTYPE_LOGPERIOD;
	      c[cn].Bls->df = (log(c[cn].Bls->maxper) - log(c[cn].Bls->minper))/ (c[cn].Bls->nf - 1);
	    }
	    else
	      i--;
	  }
	  else
	    i--;
	  c[cn].Bls->adjust_qmin_mindt = 0;
	  c[cn].Bls->reduce_nb = 0;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"adjust-qmin-by-mindt")) {
	      c[cn].Bls->adjust_qmin_mindt = 1;
	      i++;
	      if(i < argc) {
		if(!strcmp(argv[i],"reduce-nbins")) {
		  c[cn].Bls->reduce_nb = 1;
		} 
		else
		  i--;
	      }
	      else
		i--;
	    }
	    else
	      i--;
	  }
	  else
	    i--;
	  c[cn].Bls->reportharmonics = 0;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"reportharmonics")) {
	      c[cn].Bls->reportharmonics = 1;
	      }
	    else
	      i--;
	  }
	  else
	    i--;
	  cn++;
	}

      /* -BLSFixPer <\"aov\" | \"ls\" | \"list\" [\"column\" col] | \"fix\" period | \"fixcolumn\" <colname | colnum> | \"expr\" expr> <\"r\" rmin rmax | \"q\" qmin qmax > nbins timezone omodel [model_outdir] correctlc [\"fittrap\"] */
      else if(!strncmp(argv[i],"-BLSFixPer",10) && strlen(argv[i]) == 10)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].require_sort = 1;
	  c[cn].require_distinct = 1;
	  c[cn].cnum = CNUM_FIXPERBLS;
	  if((c[cn].BlsFixPer = (_BlsFixPer *) malloc(sizeof(_BlsFixPer))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"aov",3) && strlen(argv[i]) == 3)
		{
		  c[cn].BlsFixPer->pertype = PERTYPE_AOV;
		  m = -1;
		  /* Get the index for the last aov command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_AOV || c[l].cnum == CNUM_HARMAOV)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOAOV);
		  else
		    c[cn].BlsFixPer->lastaovindex = m;
		}
	      else if(!strncmp(argv[i],"ls",2) && strlen(argv[i]) == 2)
		{
		  c[cn].BlsFixPer->pertype = PERTYPE_LS;
		  m = -1;
		  /* Get the index for the last ls command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_LS)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOLS);
		  else
		    c[cn].BlsFixPer->lastlsindex = m;
		}
	      else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
		{
		  c[cn].BlsFixPer->pertype = PERTYPE_SPECIFIED;
		  k = 0;
		  i++;
		  if(i < argc) {
		    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
		      i++;
		      if(i < argc)
			k = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		    } else i--;
		  } else i--;
		  RegisterDataFromInputList(p,
					    (void *) (&(c[cn].BlsFixPer->period)),
					    VARTOOLS_TYPE_DOUBLE,
					    1, cn, 0, 0, NULL, k,
					    "BLSFIXPER_PERIOD");
		}
	      else if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
		{
		  c[cn].BlsFixPer->pertype = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    {
		      c[cn].BlsFixPer->perfix = atof(argv[i]);
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strcmp(argv[i],"fixcolumn"))
		{
		  c[cn].BlsFixPer->pertype = PERTYPE_FIXCOLUMN;
		  i++;
		  if(i < argc)
		    increaselinkedcols(p, &(c[cn].BlsFixPer->linkedcolumn), argv[i], cn);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strcmp(argv[i],"expr"))
		{
		  c[cn].BlsFixPer->pertype = PERTYPE_EXPR;
		  i++;
		  if(i < argc)
		    parse_setparam_expr(&(c[cn]), argv[i], &(c[cn].BlsFixPer->perexpr));
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"r",1) && strlen(argv[i]) == 1)
		{
		  c[cn].BlsFixPer->rflag = 1;
		  i++;
		  if(i < argc)
		    c[cn].BlsFixPer->rmin = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].BlsFixPer->rmax = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"q",1) && strlen(argv[i]) == 1)
		{
		  c[cn].BlsFixPer->rflag = 0;
		  i++;
		  if(i < argc)
		    c[cn].BlsFixPer->qmin = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].BlsFixPer->qmax = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  i++;
	  if(i < argc)
	    c[cn].BlsFixPer->nbins = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].BlsFixPer->timezone = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].BlsFixPer->omodel = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].BlsFixPer->omodel)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].BlsFixPer->modeloutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].BlsFixPer->modelsuffix,".blsfixper.model");
	    }
	  i++;
	  if(i < argc)
	    c[cn].BlsFixPer->correctlc = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  c[cn].BlsFixPer->fittrap = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strcmp(argv[i],"fittrap")) {
		c[cn].BlsFixPer->fittrap = 1;
	      }
	      else
		i--;
	    }
	  else
	    i--;
	  cn++;
	}

      /* -BLSFixDurTc <\"duration\" < \"fix\" dur | \"fixcolumn\" < colname | colnum > | \"list\" [\"column\" col] >> <\"Tc\" < \"fix\" tcval | \"fixcolumn\" < colname | colnum > | \"list\" [\"column\" col] >> [ \"fixdepth\" <\"fix\" depth | \"fixcolumn\" < colname | colnum > | \"list\" [\"column\" col] [ \"qgress\" <\"fix\" depth | \"fixcolumn\" < colname | colnum > | \"list\" [\"column\" col]] ] minper maxper nfreq nbins timezone Npeak outperiodogram [outdir] omodel [modeloutdir] correctlc [\"fittrap\"] [\"nobinnedrms\"] [\"ophcurve\" phmin phmax phstep] [\"ojdcurve\" jdstep] */

      else if(!strcmp(argv[i],"-BLSFixDurTc"))
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].require_sort = 1;
	  c[cn].require_distinct = 1;
	  c[cn].cnum = CNUM_BLSFIXDURTC;
	  if((c[cn].BlsFixDurTc = (_BlsFixDurTc *) malloc(sizeof(_BlsFixDurTc))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i >= argc)
	    listcommands(argv[iterm],p);
	  if(strcmp(argv[i],"duration"))
	    listcommands(argv[iterm],p);
	  i++;
	  if(i >= argc)
	    listcommands(argv[iterm],p);
	  if(!strcmp(argv[i],"fix")) {
	    c[cn].BlsFixDurTc->durtype = PERTYPE_FIX;
	    i++;
	    if(i >= argc)
	      listcommands(argv[iterm],p);
	    c[cn].BlsFixDurTc->fixdur = atof(argv[i]);
	  }
	  else if(!strcmp(argv[i],"fixcolumn")) {
	    c[cn].BlsFixDurTc->durtype = PERTYPE_FIXCOLUMN;
	    i++;
	    if(i >= argc)
	      listcommands(argv[iterm],p);
	    increaselinkedcols(p, &(c[cn].BlsFixDurTc->fixdur_linkedcolumn), argv[i], cn);
	  }
	  else if(!strcmp(argv[i],"list")) {
	    c[cn].BlsFixDurTc->durtype = PERTYPE_SPECIFIED;
	    k = 0;
	    i++;
	    if(i < argc) {
	      if(!strcmp(argv[i],"column")) {
		i++;
		if(i < argc)
		  k = atoi(argv[i]);
		else
		  listcommands(argv[iterm],p);
	      } else i--;
	    } else i--;
	    RegisterDataFromInputList(p,
				      (void *) (&(c[cn].BlsFixDurTc->inputdur)),
				      VARTOOLS_TYPE_DOUBLE,
				      0, cn, 0, 0, NULL, k,
				      "BLSFIXDURPHASE_DURATION");
	  }
	  else
	    listcommands(argv[iterm],p);

	  i++;
	  if(i >= argc)
	    listcommands(argv[iterm],p);
	  if(strcmp(argv[i],"Tc"))
	    listcommands(argv[iterm],p);
	  i++;
	  if(i >= argc)
	    listcommands(argv[iterm],p);
	  if(!strcmp(argv[i],"fix")) {
	    c[cn].BlsFixDurTc->TCtype = PERTYPE_FIX;
	    i++;
	    if(i >= argc)
	      listcommands(argv[iterm],p);
	    c[cn].BlsFixDurTc->fixTC = atof(argv[i]);
	  }
	  else if(!strcmp(argv[i],"fixcolumn")) {
	    c[cn].BlsFixDurTc->TCtype = PERTYPE_FIXCOLUMN;
	    i++;
	    if(i >= argc)
	      listcommands(argv[iterm],p);
	    increaselinkedcols(p, &(c[cn].BlsFixDurTc->fixTC_linkedcolumn), argv[i], cn);
	  }
	  else if(!strcmp(argv[i],"list")) {
	    c[cn].BlsFixDurTc->TCtype = PERTYPE_SPECIFIED;
	    k = 0;
	    i++;
	    if(i < argc) {
	      if(!strcmp(argv[i],"column")) {
		i++;
		if(i < argc)
		  k = atoi(argv[i]);
		else
		  listcommands(argv[iterm],p);
	      } else i--;
	    } else i--;
	    RegisterDataFromInputList(p,
				      (void *) (&(c[cn].BlsFixDurTc->inputTC)),
				      VARTOOLS_TYPE_DOUBLE,
				      0, cn, 0, 0, NULL, k,
				      "BLSFIXDURPHASE_TC");
	  }
	  else
	    listcommands(argv[iterm],p);

	  c[cn].BlsFixDurTc->fixdepth = 0;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"fixdepth")) {
	      c[cn].BlsFixDurTc->fixdepth = 1;
	      i++;
	      if(i >= argc)
		listcommands(argv[iterm],p);
	      if(!strcmp(argv[i],"fix")) {
		c[cn].BlsFixDurTc->depthtype = PERTYPE_FIX;
		i++;
		if(i >= argc)
		  listcommands(argv[iterm],p);
		c[cn].BlsFixDurTc->fixdepthval = atof(argv[i]);
	      }
	      else if(!strcmp(argv[i],"fixcolumn")) {
		c[cn].BlsFixDurTc->depthtype = PERTYPE_FIXCOLUMN;
		i++;
		if(i >= argc)
		  listcommands(argv[iterm],p);
		increaselinkedcols(p, &(c[cn].BlsFixDurTc->fixdepth_linkedcolumn), argv[i], cn);
	      }
	      else if(!strcmp(argv[i],"list")) {
		c[cn].BlsFixDurTc->depthtype = PERTYPE_SPECIFIED;
		k = 0;
		i++;
		if(i < argc) {
		  if(!strcmp(argv[i],"column")) {
		    i++;
		    if(i < argc)
		      k = atoi(argv[i]);
		    else
		      listcommands(argv[iterm],p);
		  } else i--;
		} else i--;
		RegisterDataFromInputList(p,
					  (void *) (&(c[cn].BlsFixDurTc->inputdepth)),
				      VARTOOLS_TYPE_DOUBLE,
					  0, cn, 0, 0, NULL, k,
					  "BLSFIXDURPHASE_DEPTH");
	      }
	      else
		listcommands(argv[iterm],p);
	      
	      c[cn].BlsFixDurTc->qgresstype = PERTYPE_FIX;
	      c[cn].BlsFixDurTc->qgressval = 0.0;
	      i++;
	      if(i < argc) {
		if(!strcmp(argv[i],"qgress")) {
		  i++;
		  if(i >= argc)
		    listcommands(argv[iterm],p);
		  if(!strcmp(argv[i],"fix")) {
		    c[cn].BlsFixDurTc->qgresstype = PERTYPE_FIX;
		    i++;
		    if(i >= argc)
		      listcommands(argv[iterm],p);
		    c[cn].BlsFixDurTc->qgressval = atof(argv[i]);
		  }
		  else if(!strcmp(argv[i],"fixcolumn")) {
		    c[cn].BlsFixDurTc->qgresstype = PERTYPE_FIXCOLUMN;
		    i++;
		    if(i >= argc)
		      listcommands(argv[iterm],p);
		    increaselinkedcols(p, &(c[cn].BlsFixDurTc->fixqgress_linkedcolumn), argv[i], cn);
		  }
		  else if(!strcmp(argv[i],"list")) {
		    c[cn].BlsFixDurTc->qgresstype = PERTYPE_SPECIFIED;
		    k = 0;
		    i++;
		    if(i < argc) {
		      if(!strcmp(argv[i],"column")) {
			i++;
			if(i < argc)
			  k = atoi(argv[i]);
			else
			  listcommands(argv[iterm],p);
		      } else i--;
		    } else i--;
		    RegisterDataFromInputList(p,
					      (void *) (&(c[cn].BlsFixDurTc->inputqgress)),
					      VARTOOLS_TYPE_DOUBLE,
					      0, cn, 0, 0, NULL, k,
					      "BLSFIXDURPHASE_DEPTH");
		  }
		  else
		    listcommands(argv[iterm],p);
		} else 
		  i--;
	      } else
		i--;
	    } else 
	      i--;
	  } else
	    i--;
	  

	  i++;
	  if(i < argc)
	    c[cn].BlsFixDurTc->minper = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].BlsFixDurTc->maxper = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].BlsFixDurTc->nf = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
#ifndef PARALLEL
	  if((c[cn].BlsFixDurTc->p = (double *) malloc((c[cn].BlsFixDurTc->nf+1)*sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
#endif
	  c[cn].BlsFixDurTc->df = ((1./c[cn].BlsFixDurTc->minper) - (1./c[cn].BlsFixDurTc->maxper)) / (c[cn].BlsFixDurTc->nf - 1);
	  i++;
	  if(i < argc)
	    c[cn].BlsFixDurTc->timezone = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].BlsFixDurTc->Npeak = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].BlsFixDurTc->operiodogram = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].BlsFixDurTc->operiodogram)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].BlsFixDurTc->outdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].BlsFixDurTc->suffix,".blsfixdurtc");
	    }
	  i++;
	  if(i < argc)
	    c[cn].BlsFixDurTc->omodel = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].BlsFixDurTc->omodel)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].BlsFixDurTc->modeloutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].BlsFixDurTc->modelsuffix,".blsfixdurtc.model");
	    }
	  i++;
	  if(i < argc)
	    c[cn].BlsFixDurTc->correctlc = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  c[cn].BlsFixDurTc->fittrap = 0;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"fittrap")) {
	      c[cn].BlsFixDurTc->fittrap = 1;
	    }
	    else
	      i--;
	  }
	  else
	    i--;
	  c[cn].BlsFixDurTc->ophcurve = 0;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"ophcurve")) {
	      c[cn].BlsFixDurTc->ophcurve = 1;
	      i++;
	      if(i < argc)
		sprintf(c[cn].BlsFixDurTc->ophcurveoutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].BlsFixDurTc->ophcurvesuffix,".blsfixdurtc.phcurve");
	      i++;
	      if(i < argc)
		c[cn].BlsFixDurTc->phmin = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		c[cn].BlsFixDurTc->phmax = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		c[cn].BlsFixDurTc->phstep = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	    }
	    else
	      i--;
	  }
	  else
	    i--;
	  c[cn].BlsFixDurTc->ojdcurve = 0;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"ojdcurve")) {
	      c[cn].BlsFixDurTc->ojdcurve = 1;
	      i++;
	      if(i < argc)
		sprintf(c[cn].BlsFixDurTc->ojdcurveoutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].BlsFixDurTc->ojdcurvesuffix,".blsfixdurtc.jdcurve");
	      i++;
	      if(i < argc)
		c[cn].BlsFixDurTc->jdstep = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	    }
	    else
	      i--;
	  }
	  else
	    i--;
	  cn++;
	}
      else if(!strncmp(argv[i],"-BLSFixPer",10) && strlen(argv[i]) == 10)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].require_sort = 1;
	  c[cn].require_distinct = 1;
	  c[cn].cnum = CNUM_FIXPERBLS;
	  if((c[cn].BlsFixPer = (_BlsFixPer *) malloc(sizeof(_BlsFixPer))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"aov",3) && strlen(argv[i]) == 3)
		{
		  c[cn].BlsFixPer->pertype = PERTYPE_AOV;
		  m = -1;
		  /* Get the index for the last aov command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_AOV || c[l].cnum == CNUM_HARMAOV)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOAOV);
		  else
		    c[cn].BlsFixPer->lastaovindex = m;
		}
	      else if(!strncmp(argv[i],"ls",2) && strlen(argv[i]) == 2)
		{
		  c[cn].BlsFixPer->pertype = PERTYPE_LS;
		  m = -1;
		  /* Get the index for the last ls command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_LS)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOLS);
		  else
		    c[cn].BlsFixPer->lastlsindex = m;
		}
	      else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
		{
		  c[cn].BlsFixPer->pertype = PERTYPE_SPECIFIED;
		  k = 0;
		  i++;
		  if(i < argc) {
		    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
		      i++;
		      if(i < argc)
			k = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		    } else i--;
		  } else i--;
		  RegisterDataFromInputList(p,
					    (void *) (&(c[cn].BlsFixPer->period)),
					    VARTOOLS_TYPE_DOUBLE,
					    1, cn, 0, 0, NULL, k,
					    "BLSFIXPER_PERIOD");
		}
	      else if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
		{
		  c[cn].BlsFixPer->pertype = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    {
		      c[cn].BlsFixPer->perfix = atof(argv[i]);
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"r",1) && strlen(argv[i]) == 1)
		{
		  c[cn].BlsFixPer->rflag = 1;
		  i++;
		  if(i < argc)
		    c[cn].BlsFixPer->rmin = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].BlsFixPer->rmax = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"q",1) && strlen(argv[i]) == 1)
		{
		  c[cn].BlsFixPer->rflag = 0;
		  i++;
		  if(i < argc)
		    c[cn].BlsFixPer->qmin = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].BlsFixPer->qmax = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  i++;
	  if(i < argc)
	    c[cn].BlsFixPer->nbins = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].BlsFixPer->timezone = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].BlsFixPer->omodel = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].BlsFixPer->omodel)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].BlsFixPer->modeloutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].BlsFixPer->modelsuffix,".blsfixper.model");
	    }
	  i++;
	  if(i < argc)
	    c[cn].BlsFixPer->correctlc = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  c[cn].BlsFixPer->fittrap = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strcmp(argv[i],"fittrap")) {
		c[cn].BlsFixPer->fittrap = 1;
	      }
	      else
		i--;
	    }
	  else
	    i--;
	  cn++;
	}
	  

      /* -SoftenedTransit <"bls" | "blsfixper" | P0 T00 eta0 delta0 mconst0 cval0> fitephem fiteta fitcval fitdelta fitmconst correctlc omodel [modeloutdir] fit_harm [nharm nsubharm] */
      else if(!strncmp(argv[i],"-SoftenedTransit",16) && strlen(argv[i]) == 16)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_SOFTENEDTRANSIT;
	  if((c[cn].SoftenedTransit = (_SoftenedTransit *) malloc(sizeof(_SoftenedTransit))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)

	    {
	      if(!strncmp(argv[i],"bls",3) && strlen(argv[i]) == 3)
		{
		  c[cn].SoftenedTransit->frombls = 1;
		  c[cn].SoftenedTransit->fromblsfixper = 0;
		  c[cn].SoftenedTransit->cval0 = 100.;
		  m = -1;
		  /* Get the index for the last bls command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_BLS)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOBLS);
		  else
		    c[cn].SoftenedTransit->lastblsindex = m;
		}
	      else if(!strncmp(argv[i],"blsfixper",9) && strlen(argv[i]) == 9)
		{
		  c[cn].SoftenedTransit->frombls = 0;
		  c[cn].SoftenedTransit->fromblsfixper = 1;
		  c[cn].SoftenedTransit->cval0 = 100.;
		  m = -1;
		  /* Get the index for the last blsfixper command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_FIXPERBLS)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOBLSFIXPER);
		  else
		    c[cn].SoftenedTransit->lastblsfixperindex = m;
		}
	      else
		{
		  c[cn].SoftenedTransit->frombls = 0;
		  c[cn].SoftenedTransit->fromblsfixper = 0;
		  c[cn].SoftenedTransit->period0 = atof(argv[i]);
		  i++;
		  if(i < argc)
		    c[cn].SoftenedTransit->T00 = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].SoftenedTransit->eta0 = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].SoftenedTransit->delta0 = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].SoftenedTransit->mconst0 = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].SoftenedTransit->cval0 = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].SoftenedTransit->fitephem = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].SoftenedTransit->fiteta = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].SoftenedTransit->fitcval = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].SoftenedTransit->fitdelta = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].SoftenedTransit->fitmconst = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].SoftenedTransit->correctlc = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].SoftenedTransit->omodel = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].SoftenedTransit->omodel)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].SoftenedTransit->modeloutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].SoftenedTransit->modelsuffix,".softenedtransit.model");
	    }
	  i++;
	  if(i < argc)
	    c[cn].SoftenedTransit->dokillharm = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].SoftenedTransit->dokillharm)
	    {
	      i++;
	      if(i < argc)
		{
		  if(!strncmp(argv[i],"aov",3) && strlen(argv[i]) == 3)
		    {
		      c[cn].SoftenedTransit->pertype = PERTYPE_AOV;
		      m = -1;
		      /* Get the index for the last aov command */
		      for(l = 0; l < cn; l++)
			{
			  if(c[l].cnum == CNUM_AOV || c[l].cnum == CNUM_HARMAOV)
			    m = l;
			}
		      if(m < 0)
			error(ERR_KILLHARM_NOAOV);
		      else
			c[cn].SoftenedTransit->lastaovindex = m;
		    }
		  else if(!strncmp(argv[i],"ls",2) && strlen(argv[i]) == 2)
		    {
		      c[cn].SoftenedTransit->pertype = PERTYPE_LS;
		      m = -1;
		      /* Get the index for the last aov command */
		      for(l = 0; l < cn; l++)
			{
			  if(c[l].cnum == CNUM_LS)
			    m = l;
			}
		      if(m < 0)
			error(ERR_KILLHARM_NOLS);
		      else
			c[cn].SoftenedTransit->lastlsindex = m;
		    }
		  else if(!strncmp(argv[i],"bls",3) && strlen(argv[i]) == 3)
		    {
		      c[cn].SoftenedTransit->pertype = PERTYPE_BLS;
		      m = -1;
		      /* Get the index for the last aov command */
		      for(l = 0; l < cn; l++)
			{
			  if(c[l].cnum == CNUM_BLS)
			    m = l;
			}
		      if(m < 0)
			error(ERR_KILLHARM_NOBLS);
		      else
			c[cn].SoftenedTransit->lastlsindex = m;
		    }
		  else if(!strcmp(argv[i],"fix"))
		    {
		      c[cn].SoftenedTransit->pertype = PERTYPE_FIX;
		      i++;
		      if(i < argc)
			c[cn].SoftenedTransit->per_harm = atof(argv[i]);
		      else
			listcommands(argv[iterm],p);
		    }
		  else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
		    {
		      c[cn].SoftenedTransit->pertype = PERTYPE_SPECIFIED;
		      k = 0;
		      i++;
		      if(i < argc) {
			if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
			  i++;
			  if(i < argc)
			    k = atoi(argv[i]);
			  else
			    listcommands(argv[iterm],p);
			} else i--;
		      } else i--;
		      RegisterDataFromInputList(p,
						(void *) (&(c[cn].SoftenedTransit->per_harm_spec)),
						VARTOOLS_TYPE_DOUBLE,
						0, cn, 0, 0, NULL, k,
						"SOFTENEDTRANSIT_PERHARM");
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		c[cn].SoftenedTransit->nharm = atoi(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		c[cn].SoftenedTransit->nsubharm = atoi(argv[i]);
	      else
		listcommands(argv[iterm],p);
	    }
	  cn++;
	} 

      /* -MandelAgolTransit <bls | blsfixper | P0 T00 r0 a0 <\"i\" inclination | \"b\" bimpact> e0 omega0 mconst0> <\"quad\" | \"nonlin\"> ldcoeff1_0 ... ldcoeffn_0 fitephem fitr fita fitinclterm fite fitomega fitmconst fitldcoeff1 ... fitldcoeffn fitRV [RVinputfile RVmodeloutfile K0 gamma0 fitK fitgamma] correctlc omodel [model_outdir] [\"modelvar\" var] [\"ophcurve\" model_outdir phmin phmax phstep] [\"ojdcurve\" model_outdir jdstep]*/
      else if(!strncmp(argv[i],"-MandelAgolTransit",18) && strlen(argv[i]) == 18)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_MANDELAGOLTRANSIT;
	  if((c[cn].MandelAgolTransit = (_MandelAgolTransit *) malloc(sizeof(_MandelAgolTransit))) == NULL)
	    error(ERR_MEMALLOC);
	  c[cn].MandelAgolTransit->refititer = 1;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"bls",3) && strlen(argv[i]) == 3)
		{
		  c[cn].MandelAgolTransit->frombls = 1;
		  c[cn].MandelAgolTransit->fromblsfixper = 0;
		  m = -1;
		  /* Get the index for the last bls command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_BLS)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOBLS);
		  else
		    c[cn].MandelAgolTransit->lastblsindex = m;
		}
	      else if(!strncmp(argv[i],"blsfixper",9) && strlen(argv[i]) == 9)
		{
		  c[cn].MandelAgolTransit->frombls = 0;
		  c[cn].MandelAgolTransit->fromblsfixper = 1;
		  m = -1;
		  /* Get the index for the last blsfixper command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_FIXPERBLS)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOBLSFIXPER);
		  else
		    c[cn].MandelAgolTransit->lastblsfixperindex = m;
		}
	      else
		{
		  c[cn].MandelAgolTransit->frombls = 0;
		  c[cn].MandelAgolTransit->fromblsfixper = 0;
		  c[cn].MandelAgolTransit->P0 = atof(argv[i]);
		  i++;
		  if(i < argc)
		    c[cn].MandelAgolTransit->T00 = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].MandelAgolTransit->r0 = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].MandelAgolTransit->a0 = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    {
		      if(!strcmp(argv[i],"i"))
			{
			  i++;
			  if(i < argc) {
			    c[cn].MandelAgolTransit->inc0 = atof(argv[i]);
			    c[cn].MandelAgolTransit->inputinclterm = 0;
			  }
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strcmp(argv[i],"b"))
			{
			  i++;
			  if(i < argc) {
			    c[cn].MandelAgolTransit->bimpact0 = atof(argv[i]);
			    c[cn].MandelAgolTransit->inputinclterm = 1;
			  }
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			listcommands(argv[iterm],p);
		    }
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].MandelAgolTransit->e0 = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].MandelAgolTransit->omega0 = atof(argv[i])*M_PI/180.;
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].MandelAgolTransit->mconst0 = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	    }
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"quad",4) && strlen(argv[i]) == 4)
		{
		  c[cn].MandelAgolTransit->type = 0;
		  c[cn].MandelAgolTransit->nldcoeff = 2;
		}
	      else if(!strncmp(argv[i],"nonlin",6) && strlen(argv[i]) == 6)
		{
		  c[cn].MandelAgolTransit->type = 1;
		  c[cn].MandelAgolTransit->nldcoeff = 4;
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  i++;
	  if(i < argc)
	    c[cn].MandelAgolTransit->ldcoeffs0[0] = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].MandelAgolTransit->ldcoeffs0[1] = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].MandelAgolTransit->type == 1)
	    {
	      i++;
	      if(i < argc)
		c[cn].MandelAgolTransit->ldcoeffs0[2] = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		c[cn].MandelAgolTransit->ldcoeffs0[3] = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    {
	      c[cn].MandelAgolTransit->ldcoeffs0[2] = 0.;
	      c[cn].MandelAgolTransit->ldcoeffs0[3] = 0.;
	    }
	  i++;
	  if(i < argc)
	    c[cn].MandelAgolTransit->fitephem = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].MandelAgolTransit->fitr = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].MandelAgolTransit->fita = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].MandelAgolTransit->fitinclterm = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].MandelAgolTransit->fite = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].MandelAgolTransit->fitomega = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].MandelAgolTransit->fitmconst = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].MandelAgolTransit->fitldcoeffs[0] = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].MandelAgolTransit->fitldcoeffs[1] = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].MandelAgolTransit->type == 1)
	    {
	      i++;
	      if(i < argc)
		c[cn].MandelAgolTransit->fitldcoeffs[2] = atoi(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		c[cn].MandelAgolTransit->fitldcoeffs[3] = atoi(argv[i]);
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    {
	      c[cn].MandelAgolTransit->fitldcoeffs[2] = 0;
	      c[cn].MandelAgolTransit->fitldcoeffs[3] = 0;
	    }
	  i++;
	  if(i < argc)
	    c[cn].MandelAgolTransit->fitRV = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].MandelAgolTransit->fitRV)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].MandelAgolTransit->RVinputfile,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		sprintf(c[cn].MandelAgolTransit->RVmodeloutfile,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		c[cn].MandelAgolTransit->K0 = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		c[cn].MandelAgolTransit->gamma0 = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		c[cn].MandelAgolTransit->fitK = atoi(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		c[cn].MandelAgolTransit->fitgamma = atoi(argv[i]);
	      else
		listcommands(argv[iterm],p);
	    }
	  i++;
	  if(i < argc)
	    c[cn].MandelAgolTransit->correctlc = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].MandelAgolTransit->omodel = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].MandelAgolTransit->omodel)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].MandelAgolTransit->modeloutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].MandelAgolTransit->modelsuffix,".mandelagoltransit.model");
	    }
	  c[cn].MandelAgolTransit->modelvarname = NULL;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"modelvar")) {
	      i++;
	      if(i >= argc) listcommands(argv[iterm],p);
	      if((c[cn].MandelAgolTransit->modelvarname = (char *) malloc((strlen(argv[i])+1))) == NULL)
		error(ERR_MEMALLOC);
	      sprintf(c[cn].MandelAgolTransit->modelvarname,"%s",argv[i]);
	    } else
	      i--;
	  } else 
	    i--;
	  c[cn].MandelAgolTransit->ophcurve = 0;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"ophcurve")) {
	      c[cn].MandelAgolTransit->ophcurve = 1;
	      i++;
	      if(i < argc)
		sprintf(c[cn].MandelAgolTransit->ophcurveoutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].MandelAgolTransit->ophcurvesuffix,".mandelagoltransit.phcurve");
	      i++;
	      if(i < argc)
		c[cn].MandelAgolTransit->phmin = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		c[cn].MandelAgolTransit->phmax = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		c[cn].MandelAgolTransit->phstep = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	    }
	    else
	      i--;
	  }
	  else
	    i--;
	  c[cn].MandelAgolTransit->ojdcurve = 0;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"ojdcurve")) {
	      c[cn].MandelAgolTransit->ojdcurve = 1;
	      i++;
	      if(i < argc)
		sprintf(c[cn].MandelAgolTransit->ojdcurveoutdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].MandelAgolTransit->ojdcurvesuffix,".mandelagoltransit.jdcurve");
	      i++;
	      if(i < argc)
		c[cn].MandelAgolTransit->jdstep = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	    }
	    else
	      i--;
	  }
	  else
	    i--;
	  cn++;
	} 
      
      /*-microlens [\"f0\" [\"fix\" fixval | \"list\" [\"column\" col] | \"fixcolumn\" <colname | colnum> | \"auto\"] [\"step\" initialstepsize] [\"novary\"]] [\"f1\" [\"fix\" fixval | \"list\" [\"column\" col] | \"fixcolumn\" <colname | colnum> | \"auto\"] [\"step\" initialstepsize] [\"novary\"]] [\"u0\" [\"fix\" fixval | \"list\" [\"column\" col] | \"fixcolumn\" <colname | colnum> | \"auto\"] [\"step\" initialstepsize] [\"novary\"]] [\"t0\" [\"fix\" fixval | \"list\" [\"column\" col] | \"fixcolumn\" <colname | colnum> | \"auto\"] [\"step\" initialstepsize] [\"novary\"]] [\"tmax\" [\"fix\" fixval | \"list\" [\"column\" col] | \"fixcolumn\" <colname | colnum> | \"auto\"] [\"step\" initialstepsize] [\"novary\"]] [\"step\" initialstepsize] [\"novary\"]] */
      else if(!strncmp(argv[i],"-microlens",10) && strlen(argv[i]) == 10)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_MICROLENS;
	  if((c[cn].MicroLens = (_MicroLens *) malloc(sizeof(_MicroLens))) == NULL)
	    error(ERR_MEMALLOC);
	  /* Set the initial values */
	  c[cn].MicroLens->f0_source = PERTYPE_AUTOFIND;
	  c[cn].MicroLens->f1_source = PERTYPE_AUTOFIND;
	  c[cn].MicroLens->u0_source = PERTYPE_AUTOFIND;
	  c[cn].MicroLens->t0_source = PERTYPE_AUTOFIND;
	  c[cn].MicroLens->tmax_source = PERTYPE_AUTOFIND;
	  c[cn].MicroLens->fitf0 = 1;
	  c[cn].MicroLens->fitf1 = 1;
	  c[cn].MicroLens->fitu0 = 1;
	  c[cn].MicroLens->fitt0 = 1;
	  c[cn].MicroLens->fittmax = 1;
	  c[cn].MicroLens->correctlc = 0;
	  c[cn].MicroLens->omodel = 0;
	  c[cn].MicroLens->f0_initialstep = 0;
	  c[cn].MicroLens->f1_initialstep = 0;
	  c[cn].MicroLens->u0_initialstep = 0;
	  c[cn].MicroLens->t0_initialstep = 0;
	  c[cn].MicroLens->tmax_initialstep = 0;
	  /* Check for options */
	  /* do f0 */
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"f0",2) && strlen(argv[i]) == 2)
		{
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].MicroLens->f0_source = PERTYPE_FIX;
			      c[cn].MicroLens->f00_fix = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			{
			  c[cn].MicroLens->f0_source = PERTYPE_SPECIFIED;
			  k = 0;
			  i++;
			  if(i < argc) {
			    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    } else i--;
			  } else i--;
			  RegisterDataFromInputList(p,
						    (void *) (&(c[cn].MicroLens->f00)),
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "MICROLENS_F0");
			}
		      else if(!strncmp(argv[i],"fixcolumn",9) && strlen(argv[i]) == 9)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].MicroLens->f0_source = PERTYPE_FIXCOLUMN;
			      increaselinkedcols(p, &(c[cn].MicroLens->f0_linkedcolumn), argv[i], cn);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"auto",4) && strlen(argv[i]) == 4)
			{
			  c[cn].MicroLens->f0_source = PERTYPE_AUTOFIND;
			}
		      else
			i--;
		    }
		  else
		    i--;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"step",4) && strlen(argv[i]) == 4)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].MicroLens->f0_initialstep = 1;
			      c[cn].MicroLens->f0_initialstepval = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			i--;
		    }
		  else
		    i--;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"novary",6) && strlen(argv[i]) == 6)
			{
			  c[cn].MicroLens->fitf0 = 0;
			}
		      else
			i--;
		    }
		  else
		    i--;
		}
	      else
		i--;
	    }
	  else
	    i--;
	  /* do f1 */
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"f1",2) && strlen(argv[i]) == 2)
		{
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].MicroLens->f1_source = PERTYPE_FIX;
			      c[cn].MicroLens->f10_fix = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			{
			  c[cn].MicroLens->f1_source = PERTYPE_SPECIFIED;
			  k = 0;
			  i++;
			  if(i < argc) {
			    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    } else i--;
			  } else i--;
			  RegisterDataFromInputList(p,
						    (void *) (&(c[cn].MicroLens->f10)),
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "MICROLENS_F1");
			}
		      else if(!strncmp(argv[i],"fixcolumn",9) && strlen(argv[i]) == 9)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].MicroLens->f1_source = PERTYPE_FIXCOLUMN;
			      increaselinkedcols(p, &(c[cn].MicroLens->f1_linkedcolumn), argv[i], cn);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"auto",4) && strlen(argv[i]) == 4)
			{
			  c[cn].MicroLens->f1_source = PERTYPE_AUTOFIND;
			}
		      else
			i--;
		    }
		  else
		    i--;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"step",4) && strlen(argv[i]) == 4)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].MicroLens->f1_initialstep = 1;
			      c[cn].MicroLens->f1_initialstepval = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			i--;
		    }
		  else
		    i--;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"novary",6) && strlen(argv[i]) == 6)
			{
			  c[cn].MicroLens->fitf1 = 0;
			}
		      else
			i--;
		    }
		  else
		    i--;
		}
	      else
		i--;
	    }
	  else
	    i--;
	  /* do u0 */
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"u0",2) && strlen(argv[i]) == 2)
		{
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].MicroLens->u0_source = PERTYPE_FIX;
			      c[cn].MicroLens->u00_fix = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			{
			  c[cn].MicroLens->u0_source = PERTYPE_SPECIFIED;
			  k = 0;
			  i++;
			  if(i < argc) {
			    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    } else i--;
			  } else i--;
			  RegisterDataFromInputList(p,
						    (void *) (&(c[cn].MicroLens->u00)),
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "MICROLENS_U0");
			}
		      else if(!strncmp(argv[i],"fixcolumn",9) && strlen(argv[i]) == 9)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].MicroLens->u0_source = PERTYPE_FIXCOLUMN;
			      increaselinkedcols(p, &(c[cn].MicroLens->u0_linkedcolumn), argv[i], cn);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"auto",4) && strlen(argv[i]) == 4)
			{
			  c[cn].MicroLens->u0_source = PERTYPE_AUTOFIND;
			}
		      else
			i--;
		    }
		  else
		    i--;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"step",4) && strlen(argv[i]) == 4)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].MicroLens->u0_initialstep = 1;
			      c[cn].MicroLens->u0_initialstepval = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			i--;
		    }
		  else
		    i--;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"novary",6) && strlen(argv[i]) == 6)
			{
			  c[cn].MicroLens->fitu0 = 0;
			}
		      else
			i--;
		    }
		  else
		    i--;
		}
	      else
		i--;
	    }
	  else
	    i--;
	  /* do t0 */
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"t0",2) && strlen(argv[i]) == 2)
		{
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].MicroLens->t0_source = PERTYPE_FIX;
			      c[cn].MicroLens->t00_fix = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			{
			  c[cn].MicroLens->t0_source = PERTYPE_SPECIFIED;
			  k = 0;
			  i++;
			  if(i < argc) {
			    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    } else i--;
			  } else i--;
			  RegisterDataFromInputList(p,
						    (void *) (&(c[cn].MicroLens->t00)),
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "MICROLENS_T0");
			}
		      else if(!strncmp(argv[i],"fixcolumn",9) && strlen(argv[i]) == 9)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].MicroLens->t0_source = PERTYPE_FIXCOLUMN;
			      increaselinkedcols(p, &(c[cn].MicroLens->t0_linkedcolumn), argv[i], cn);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"auto",4) && strlen(argv[i]) == 4)
			{
			  c[cn].MicroLens->t0_source = PERTYPE_AUTOFIND;
			}
		      else
			i--;
		    }
		  else
		    i--;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"step",4) && strlen(argv[i]) == 4)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].MicroLens->t0_initialstep = 1;
			      c[cn].MicroLens->t0_initialstepval = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			i--;
		    }
		  else
		    i--;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"novary",6) && strlen(argv[i]) == 6)
			{
			  c[cn].MicroLens->fitt0 = 0;
			}
		      else
			i--;
		    }
		  else
		    i--;
		}
	      else
		i--;
	    }
	  else
	    i--;
	  /* do tmax */
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"tmax",4) && strlen(argv[i]) == 4)
		{
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].MicroLens->tmax_source = PERTYPE_FIX;
			      c[cn].MicroLens->tmax0_fix = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			{
			  c[cn].MicroLens->tmax_source = PERTYPE_SPECIFIED;
			  k = 0;
			  i++;
			  if(i < argc) {
			    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    } else i--;
			  } else i--;
			  RegisterDataFromInputList(p,
						    (void *) (&(c[cn].MicroLens->tmax0)),
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "MICROLENS_TMAX");
			}
		      else if(!strncmp(argv[i],"fixcolumn",9) && strlen(argv[i]) == 9)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].MicroLens->tmax_source = PERTYPE_FIXCOLUMN;
			      increaselinkedcols(p, &(c[cn].MicroLens->tmax_linkedcolumn), argv[i], cn);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"auto",4) && strlen(argv[i]) == 4)
			{
			  c[cn].MicroLens->tmax_source = PERTYPE_AUTOFIND;
			}
		      else
			i--;
		    }
		  else
		    i--;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"step",4) && strlen(argv[i]) == 4)
			{
			  i++;
			  if(i < argc)
			    {
			      c[cn].MicroLens->tmax_initialstep = 1;
			      c[cn].MicroLens->tmax_initialstepval = atof(argv[i]);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			i--;
		    }
		  else
		    i--;
		  i++;
		  if(i < argc)
		    {
		      if(!strncmp(argv[i],"novary",6) && strlen(argv[i]) == 6)
			{
			  c[cn].MicroLens->fittmax = 0;
			}
		      else
			i--;
		    }
		  else
		    i--;
		}
	      else
		i--;
	    }
	  else
	    i--;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"correctlc",9) && strlen(argv[i]) == 9)
		{
		  c[cn].MicroLens->correctlc = 1;
		}
	      else
		i--;
	    }
	  else
	    i--;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"omodel",6) && strlen(argv[i]) == 6)
		{
		  i++;
		  if(i < argc)
		    {
		      c[cn].MicroLens->omodel = 1;
		      sprintf(c[cn].MicroLens->modeloutdir,"%s",argv[i]);
		      sprintf(c[cn].MicroLens->modelsuffix,".microlens");
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		i--;
	    }
	  else
	    i--;
	  cn++;
	}

      /* -SYSREM Ninput_color [\"column\" col1] Ninput_airmass initial_airmass_file sigma_clip1 sigma_clip2 saturation correctlc omodel [model_outdir] otrends [trend_outfile] useweights */
      else if(!strncmp(argv[i],"-SYSREM",7) && strlen(argv[i]) == 7)
	{
	  iterm = i;
	  p->readallflag = 1;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_SYSREM;
	  if((c[cn].Sysrem = (_Sysrem *) malloc(sizeof(_Sysrem))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    c[cn].Sysrem->Nsysrem_color = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  k = 0;
	  i++;
	  if(i < argc) {
	    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
	      i++;
	      if(i < argc)
		k = atoi(argv[i]);
	      else
		listcommands(argv[iterm],p);
	    } else i--;
	  } else i--;
	  if(c[cn].Sysrem->Nsysrem_color > 0) {
	    RegisterDataFromInputList(p,
				      (void *)(&(c[cn].Sysrem->initial_colors_readin)),
				      VARTOOLS_TYPE_DOUBLE,
				      c[cn].Sysrem->Nsysrem_color,
				      cn, 0, 0, NULL, k,
				      "SYSREM_INITIAL_COLOR_TERM");
	  }
	  i++;
	  if(i < argc)
	    c[cn].Sysrem->Nsysrem_airmass = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  c[cn].Sysrem->Nsysrem_total = c[cn].Sysrem->Nsysrem_color + c[cn].Sysrem->Nsysrem_airmass;
	  if(c[cn].Sysrem->Nsysrem_color < 0 || c[cn].Sysrem->Nsysrem_airmass < 0)
	    listcommands(argv[iterm],p);
	  if(c[cn].Sysrem->Nsysrem_total <= 0)
	    listcommands(argv[iterm],p);
	  if((c[cn].Sysrem->final_colors = (double **) malloc(c[cn].Sysrem->Nsysrem_total * sizeof(double *))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    sprintf(c[cn].Sysrem->dates_name,"%s",argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Sysrem->sigma_clip1 = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Sysrem->sigma_clip2 = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Sysrem->saturation = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Sysrem->correctlc = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].Sysrem->omodel = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].Sysrem->omodel)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].Sysrem->model_outdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].Sysrem->model_suffix,".sysrem.model");
	    }
	  i++;
	  if(i < argc)
	    c[cn].Sysrem->otrend = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].Sysrem->otrend)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].Sysrem->trends_outname,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	    }
	  i++;
	  if(i < argc)
	    c[cn].Sysrem->useweights = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].Sysrem->useweights < 0 || c[cn].Sysrem->useweights > 2)
	    error(ERR_SYSREMUSEWEIGHTS);
	  cn++;
	}

      /* -TFA trendlist [\"readformat\" Nskip jdcol magcol] dates_file pixelsep [\"xycol\" xcol ycol] correctlc ocoeff [coeff_outdir] omodel [model_outdir] */
      else if(!strncmp(argv[i],"-TFA",4) && strlen(argv[i]) == 4)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_TFA;
	  if((c[cn].TFA = (_TFA *) malloc(sizeof(_TFA))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    sprintf(c[cn].TFA->trend_list_name,"%s",argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  c[cn].TFA->Nskip_trend = 0;
	  c[cn].TFA->JDcol_trend = 1;
	  c[cn].TFA->magcol_trend = 2;
	  c[cn].TFA->jdcol_isfromheader = 0;
	  c[cn].TFA->magcol_isfromheader = 0;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"readformat",10) && strlen(argv[i]) == 10)
		{
		  i++;
		  if(i < argc)
		    c[cn].TFA->Nskip_trend = atoi(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc) {
		    if(isstringint(argv[i])) {
		      c[cn].TFA->JDcol_trend = atoi(argv[i]);
		    } else {
		      sprintf(c[cn].TFA->jdcol_headername,"%s",argv[i]);
		      c[cn].TFA->jdcol_isfromheader = 1;
		    }
		  }
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc) {
		    if(isstringint(argv[i])) {
		      c[cn].TFA->magcol_trend = atoi(argv[i]);
		    } else {
		      sprintf(c[cn].TFA->magcol_headername,"%s",argv[i]);
		      c[cn].TFA->magcol_isfromheader = 1;
		    }
		  } else
		    listcommands(argv[iterm],p);
		  i++;
		}
	    }
	  if(i < argc)
	    sprintf(c[cn].TFA->dates_name,"%s",argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].TFA->pixelsep = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  k = 0;
	  l = 0;
	  i++;
	  if(i < argc) {
	    if(!strncmp(argv[i],"xycol",5) && strlen(argv[i]) == 5) {
	      i++;
	      if(i < argc)
		k = atoi(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		l = atoi(argv[i]);
	      else
		listcommands(argv[iterm],p);
	    } else i--;
	  } else i--;
	  RegisterDataFromInputList(p,
				    (void *) (&(c[cn].TFA->lcx)),
				    VARTOOLS_TYPE_DOUBLE,
				    1, cn, 0, 0, NULL, k, 
				    "TFA_STAR_X_POS");
	  RegisterDataFromInputList(p,
				    (void *) (&(c[cn].TFA->lcy)),
				    VARTOOLS_TYPE_DOUBLE,
				    1, cn, 0, 0, NULL, l, 
				    "TFA_STAR_Y_POS");
	  i++;
	  if(i < argc)
	    c[cn].TFA->correctlc = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].TFA->ocoeff = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].TFA->ocoeff)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].TFA->coeff_outdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].TFA->coeff_suffix,".tfa.coeff");
	    }
	  i++;
	  if(i < argc)
	    c[cn].TFA->omodel = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].TFA->omodel)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].TFA->model_outdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].TFA->model_suffix,".tfa.model");
	    }
	  cn++;
	}

      /* -TFA_SR trendlist [\"readformat\" Nskip jdcol magcol] dates_file [\"decorr\" iterativeflag Nlcterms lccolumn1 lcorder1 ...] pixelsep [\"xycol\" colx coly] correctlc ocoeff [coeff_outdir] omodel [model_outdir] dotfafirst tfathresh maxiter <\"bin\" nbins [\"period\" <\"aov\" | \"ls\" | \"bls\" | \"list\" | \"fix\" period>] | \"signal\" filename | \"harm\" Nharm Nsubharm [\"period\" <\"aov\" | \"ls\" | \"bls\" | \"list\" | \"fix\" period>]> */
      else if(!strncmp(argv[i],"-TFA_SR",7) && strlen(argv[i]) == 7)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_TFA_SR;
	  if((c[cn].TFA_SR = (_TFA_SR *) malloc(sizeof(_TFA_SR))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    sprintf(c[cn].TFA_SR->trend_list_name,"%s",argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  c[cn].TFA_SR->Nskip_trend = 0;
	  c[cn].TFA_SR->JDcol_trend = 1;
	  c[cn].TFA_SR->magcol_trend = 2;
	  c[cn].TFA_SR->jdcol_isfromheader = 0;
	  c[cn].TFA_SR->magcol_isfromheader = 0;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"readformat",10) && strlen(argv[i]) == 10)
		{
		  i++;
		  if(i < argc)
		    c[cn].TFA_SR->Nskip_trend = atoi(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc) {
		    if(isstringint(argv[i])) {
		      c[cn].TFA_SR->JDcol_trend = atoi(argv[i]);
		    } else {
		      sprintf(c[cn].TFA_SR->jdcol_headername,"%s",argv[i]);
		      c[cn].TFA_SR->jdcol_isfromheader = 1;
		    }
		  }
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc) {
		    if(isstringint(argv[i])) {
		      c[cn].TFA_SR->magcol_trend = atoi(argv[i]);
		    } else {
		      sprintf(c[cn].TFA_SR->magcol_headername,"%s",argv[i]);
		      c[cn].TFA_SR->magcol_isfromheader = 1;
		    }
		  }
		  else
		    listcommands(argv[iterm],p);
		  i++;
		}
	    }
	  if(i < argc)
	    sprintf(c[cn].TFA_SR->dates_name,"%s",argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  c[cn].TFA_SR->decorrflag = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"decorr",6) && strlen(argv[i]) == 6)
		{
		  c[cn].TFA_SR->decorrflag = 1;
		  i++;
		  if(i < argc)
		    c[cn].TFA_SR->decorr_iterate = atoi(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  i++;
		  if(i < argc)
		    c[cn].TFA_SR->decorr_Nlcterms = atoi(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		  if(c[cn].TFA_SR->decorr_Nlcterms > 0) {
		    if((c[cn].TFA_SR->lcdecorr_terms_in = (double ***) malloc(c[cn].TFA_SR->decorr_Nlcterms * sizeof(double **))) == NULL ||
		       (c[cn].TFA_SR->decorr_lc_order = (int *) malloc(c[cn].TFA_SR->decorr_Nlcterms * sizeof(int))) == NULL ||
		       (c[cn].TFA_SR->decorr_lc_columns = (int *) malloc(c[cn].TFA_SR->decorr_Nlcterms * sizeof(int))) == NULL)
		      error(ERR_MEMALLOC);
		  }
		  for(j=0;j<c[cn].TFA_SR->decorr_Nlcterms;j++)
		    {
		      i++;
		      if(i < argc)
			c[cn].TFA_SR->decorr_lc_columns[j] = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		      i++;
		      if(i < argc)
			c[cn].TFA_SR->decorr_lc_order[j] = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		      if(c[cn].TFA_SR->decorr_lc_order[j] <= 0)
			error(ERR_WRONGORDER);
		    }
		  i++;
		}
	    }
	  if(i < argc)
	    c[cn].TFA_SR->pixelsep = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  k = 0;
	  l = 0;
	  i++;
	  if(i < argc) {
	    if(!strncmp(argv[i],"xycol",5) && strlen(argv[i]) == 5) {
	      i++;
	      if(i < argc)
		k = atoi(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc)
		l = atoi(argv[i]);
	      else
		listcommands(argv[iterm],p);
	    } else i--;
	  } else i--;
	  RegisterDataFromInputList(p,
				    (void *) (&(c[cn].TFA_SR->lcx)),
				    VARTOOLS_TYPE_DOUBLE,
				    1, cn, 0, 0, NULL, k, 
				    "TFA_SR_STAR_X_POS");
	  RegisterDataFromInputList(p,
				    (void *) (&(c[cn].TFA_SR->lcy)),
				    VARTOOLS_TYPE_DOUBLE,
				    1, cn, 0, 0, NULL, l, 
				    "TFA_SR_STAR_Y_POS");
	  i++;
	  if(i < argc)
	    c[cn].TFA_SR->correctlc = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].TFA_SR->ocoeff = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].TFA_SR->ocoeff)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].TFA_SR->coeff_outdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].TFA_SR->coeff_suffix,".tfa.coeff");
	    }
	  i++;
	  if(i < argc)
	    c[cn].TFA_SR->omodel = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].TFA_SR->omodel)
	    {
	      i++;
	      if(i < argc)
		sprintf(c[cn].TFA_SR->model_outdir,"%s",argv[i]);
	      else
		listcommands(argv[iterm],p);
	      sprintf(c[cn].TFA_SR->model_suffix,".tfa.model");
	    }
	  i++;
	  if(i < argc)
	    c[cn].TFA_SR->dotfafirst = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].TFA_SR->tfathresh = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    c[cn].TFA_SR->maxiter = atoi(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  c[cn].TFA_SR->use_bin = 0;
	  c[cn].TFA_SR->use_period = 0;
	  c[cn].TFA_SR->use_harm = 0;
	  i++;
	  if(i < argc) {
	    if(!strncmp(argv[i],"bin",3) && strlen(argv[i]) == 3)
	      {
		c[cn].TFA_SR->use_bin = 1;
		i++;
		if(i < argc)
		  c[cn].TFA_SR->nbins = atoi(argv[i]);
		else
		  listcommands(argv[iterm],p);
		i++;
		if(i < argc) {
		  if(!strncmp(argv[i],"period",6) && strlen(argv[i]) == 6)
		    {
		      c[cn].TFA_SR->use_period = 1;
		      i++;
		      if(i < argc)
			{
			  if(!strncmp(argv[i],"aov",3) && strlen(argv[i]) == 3)
			    {
			      c[cn].TFA_SR->pertype = PERTYPE_AOV;
			      m = -1;
			      /* Get the index for the last aov command */
			      for(l = 0; l < cn; l++)
				{
				  if(c[l].cnum == CNUM_AOV || c[l].cnum == CNUM_HARMAOV)
				    m = l;
				}
			      if(m < 0)
				error(ERR_KILLHARM_NOAOV);
			      else
				c[cn].TFA_SR->lastindex = m;
			    }
			  else if(!strncmp(argv[i],"ls",2) && strlen(argv[i]) == 2)
			    {
			      c[cn].TFA_SR->pertype = PERTYPE_LS;
			      m = -1;
			      /* Get the index for the last aov command */
			      for(l = 0; l < cn; l++)
				{
				  if(c[l].cnum == CNUM_LS)
				    m = l;
				}
			      if(m < 0)
				error(ERR_KILLHARM_NOLS);
			      else
				c[cn].TFA_SR->lastindex = m;
			    }
			  else if(!strncmp(argv[i],"bls",3) && strlen(argv[i]) == 3)
			    {
			      c[cn].TFA_SR->pertype = PERTYPE_BLS;
			      m = -1;
			      /* Get the index for the last aov command */
			      for(l = 0; l < cn; l++)
				{
				  if(c[l].cnum == CNUM_BLS)
				    m = l;
				}
			      if(m < 0)
				error(ERR_KILLHARM_NOBLS);
			      else
				c[cn].TFA_SR->lastindex = m;
			    }
			  else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			    {
			      c[cn].TFA_SR->pertype = PERTYPE_SPECIFIED;
			      k = 0;
			      i++;
			      if(i < argc) {
				if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
				  i++;
				  if(i < argc)
				    k = atoi(argv[i]);
				  else
				    listcommands(argv[iterm],p);
				} else i--;
			      } else i--;
			      RegisterDataFromInputList(p,
							(void *) (&(c[cn].TFA_SR->periods)),
							VARTOOLS_TYPE_DOUBLE,
							1, cn, 0, 0, NULL, k,
							"TFA_SR_PERIOD");
			    }
			  else if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			    {
			      c[cn].TFA_SR->pertype = PERTYPE_FIX;
			      i++;
			      if(i < argc)
				c[cn].TFA_SR->fixperiod = atof(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			listcommands(argv[iterm],p);
		    }
		  else
		    i--;
		}
		else
		  i--;
	      }
	    else if(!strncmp(argv[i],"signal",6) && strlen(argv[i]) == 6)
	      {
		i++;
		if(i < argc)
		  {
		    sprintf(c[cn].TFA_SR->signal_listname,"%s",argv[i]);
		  }
		else
		  listcommands(argv[iterm],p);
	      }
	    else if(!strncmp(argv[i],"harm",4) && strlen(argv[i]) == 4)
	      {
		c[cn].TFA_SR->use_harm = 1;
		i++;
		if(i < argc)
		  c[cn].TFA_SR->Nharm = atoi(argv[i]);
		else
		  listcommands(argv[iterm],p);
		i++;
		if(i < argc)
		  c[cn].TFA_SR->Nsubharm = atoi(argv[i]);
		else
		  listcommands(argv[iterm],p);
		i++;
		if(i < argc) {
		  if(!strncmp(argv[i],"period",6) && strlen(argv[i]) == 6)
		    {
		      c[cn].TFA_SR->use_period = 1;
		      i++;
		      if(i < argc)
			{
			  if(!strncmp(argv[i],"aov",3) && strlen(argv[i]) == 3)
			    {
			      c[cn].TFA_SR->pertype = PERTYPE_AOV;
			      m = -1;
			      /* Get the index for the last aov command */
			      for(l = 0; l < cn; l++)
				{
				  if(c[l].cnum == CNUM_AOV || c[l].cnum == CNUM_HARMAOV)
				    m = l;
				}
			      if(m < 0)
				error(ERR_KILLHARM_NOAOV);
			      else
				c[cn].TFA_SR->lastindex = m;
			    }
			  else if(!strncmp(argv[i],"ls",2) && strlen(argv[i]) == 2)
			    {
			      c[cn].TFA_SR->pertype = PERTYPE_LS;
			      m = -1;
			      /* Get the index for the last aov command */
			      for(l = 0; l < cn; l++)
				{
				  if(c[l].cnum == CNUM_LS)
				    m = l;
				}
			      if(m < 0)
				error(ERR_KILLHARM_NOLS);
			      else
				c[cn].TFA_SR->lastindex = m;
			    }
			  else if(!strncmp(argv[i],"bls",3) && strlen(argv[i]) == 3)
			    {
			      c[cn].TFA_SR->pertype = PERTYPE_BLS;
			      m = -1;
			      /* Get the index for the last aov command */
			      for(l = 0; l < cn; l++)
				{
				  if(c[l].cnum == CNUM_BLS)
				    m = l;
				}
			      if(m < 0)
				error(ERR_KILLHARM_NOBLS);
			      else
				c[cn].TFA_SR->lastindex = m;
			    }
			  else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			    {
			      c[cn].TFA_SR->pertype = PERTYPE_SPECIFIED;
			      k = 0;
			      i++;
			      if(i < argc) {
				if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
				  i++;
				  if(i < argc)
				    k = atoi(argv[i]);
				  else
				    listcommands(argv[iterm],p);
				} else i--;
			      } else i--;
			      RegisterDataFromInputList(p,
							(void *) (&(c[cn].TFA_SR->periods)),
							VARTOOLS_TYPE_DOUBLE,
							1, cn, 0, 0, NULL, k,
							"TFA_SR_PERIOD");
			    }
			  else if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			    {
			      c[cn].TFA_SR->pertype = PERTYPE_FIX;
			      i++;
			      if(i < argc)
				c[cn].TFA_SR->fixperiod = atof(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    }
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			listcommands(argv[iterm],p);
		    }
		  else
		    i--;
		}
		else
		  i--;
	      }
	  }
	  else
	    listcommands(argv[iterm],p);
	    
	  cn++;
	}

      /* -binlc combinetype <\"binsize\" binsize | \"nbins\" nbins> [\"bincolumns\" var1[:stats1][,var2[:stats2],...]] [\"firstbinshift\" firstbinshift] <\"tcenter\" | \"taverage\" | \"tmedian\" | \"tnoshrink\" \"bincolumnsonly\">*/
      else if(!strncmp(argv[i],"-binlc",6) && strlen(argv[i]) == 6)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].require_sort = 1;
	  c[cn].cnum = CNUM_BINLC;
	  if((c[cn].Binlc = (_Binlc *) malloc(sizeof(_Binlc))) == NULL)
	    error(ERR_MEMALLOC);
	  c[cn].Binlc->only_bin_columns = 0;
	  i++;
	  if(i < argc)
	    {
	      if(!strcmp(argv[i],"average")) {
		c[cn].Binlc->medflag = VARTOOLS_BINLC_BINTYPE_AVERAGE;
	      }
	      else if(!strcmp(argv[i],"median")) {
		c[cn].Binlc->medflag = VARTOOLS_BINLC_BINTYPE_MEDIAN;
	      }
	      else if(!strcmp(argv[i],"weightedaverage")) {
		c[cn].Binlc->medflag = VARTOOLS_BINLC_BINTYPE_WEIGHTEDAVERAGE;
	      }
	      else
		c[cn].Binlc->medflag = atoi(argv[i]);
	    }
	  else
	    listcommands(argv[iterm],p);
	  if(c[cn].Binlc->medflag < VARTOOLS_BINLC_BINTYPE_AVERAGE || c[cn].Binlc->medflag > VARTOOLS_BINLC_BINTYPE_WEIGHTEDAVERAGE)
	    listcommands(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"binsize",7) && strlen(argv[i]) == 7)
		{
		  c[cn].Binlc->binsize_Nbins_flag = 0;
		  i++;
		  if(i < argc)
		    c[cn].Binlc->binsize = atof(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"nbins",5) && strlen(argv[i]) == 5)
		{
		  c[cn].Binlc->binsize_Nbins_flag = 1;
		  i++;
		  if(i < argc)
		    c[cn].Binlc->Nbins = atoi(argv[i]);
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		listcommands(argv[iterm],p);
	    }
	  else
	    listcommands(argv[iterm],p);
	  c[cn].Binlc->binvarstring = NULL;
	  c[cn].Binlc->Nvar = 0;
	  i++;
	  if(i >= argc)
	    listcommands(argv[iterm],p);
	  if(!strcmp(argv[i],"bincolumns")) {
	    i++;
	    if(i >= argc)
	      listcommands(argv[iterm],p);
	    if((c[cn].Binlc->binvarstring = (char *) malloc((strlen(argv[i])+1)*sizeof(char))) == NULL)
	      error(ERR_MEMALLOC);
	    sprintf(c[cn].Binlc->binvarstring,"%s",argv[i]);
	    if(binlc_parsevarstring(c[cn].Binlc))
	      listcommands(argv[iterm],p);
	  } else 
	    i--;
	  i++;
	  if(i >= argc)
	    listcommands(argv[iterm],p);
	  if(!strncmp(argv[i],"firstbinshift",13) && strlen(argv[i]) == 13)
	    {
	      c[cn].Binlc->firstbinflag = 1;
	      i++;
	      if(i < argc)
		c[cn].Binlc->firstbin = atof(argv[i]);
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i >= argc)
		listcommands(argv[iterm],p);
	    }
	  else
	    c[cn].Binlc->firstbinflag = 0;
	  if(!strcmp(argv[i],"tcenter"))
	    c[cn].Binlc->tflag = VARTOOLS_BINLC_TIMETYPE_CENTER;
	  else if(!strcmp(argv[i],"taverage"))
	    c[cn].Binlc->tflag = VARTOOLS_BINLC_TIMETYPE_AVERAGE;
	  else if(!strcmp(argv[i],"tmedian"))
	    c[cn].Binlc->tflag = VARTOOLS_BINLC_TIMETYPE_MEDIAN;
	  else if(!strcmp(argv[i],"tnoshrink")) {
	    c[cn].Binlc->tflag = VARTOOLS_BINLC_TIMETYPE_NOSHRINK;
	    i++;
	    if(i < argc) {
	      if(!strcmp(argv[i],"bincolumnsonly")) {
		c[cn].Binlc->only_bin_columns = 1;
	      } else
		i--;
	    } else
	      i--;
	  }
	  else {
	    c[cn].Binlc->tflag = atoi(argv[i]);
	    if(c[cn].Binlc->tflag == VARTOOLS_BINLC_TIMETYPE_NOSHRINK) {
	      i++;
	      if(i < argc) {
		if(!strcmp(argv[i],"bincolumnsonly")) {
		  c[cn].Binlc->only_bin_columns = 1;
		} else
		  i--;
	      } else
		i--;
	    }
	  }
	  if(c[cn].Binlc->tflag < VARTOOLS_BINLC_TIMETYPE_CENTER || 
	     c[cn].Binlc->tflag > VARTOOLS_BINLC_TIMETYPE_NOSHRINK)
	    listcommands(argv[iterm],p);
	  cn++;
	}
	    
      /* -medianfilter time [\"average\" | \"weightedaverage\"] */
      else if(!strncmp(argv[i],"-medianfilter",13) && strlen(argv[i]) == 13)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].require_sort = 1;
	  c[cn].cnum = CNUM_MEDIANFILTER;
	  if((c[cn].MedianFilter = (_MedianFilter *) malloc(sizeof(_MedianFilter))) == NULL)
	    error(ERR_MEMALLOC);
	  i++;
	  if(i < argc)
	    c[cn].MedianFilter->time = atof(argv[i]);
	  else
	    listcommands(argv[iterm],p);
	  c[cn].MedianFilter->usemean=0;
	  c[cn].MedianFilter->replace=0;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"average",7) && strlen(argv[i]) == 7)
		{
		  c[cn].MedianFilter->usemean=1;
		}
	      else if(!strncmp(argv[i],"weightedaverage",15) && strlen(argv[i]) == 15)
		{
		  c[cn].MedianFilter->usemean=2;
		}
	      else
		i--;
	    }
	  else
	    i--;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"replace",7) && strlen(argv[i]) == 7)
		{
		  c[cn].MedianFilter->replace=1;
		}
	      else
		i--;
	    }
	  else
	    i--;
	  cn++;
	}
	    
      /* -wwz <\"maxfreq\" <\"auto\" | maxfreq>> <\"freqsamp\" freqsamp> <\"tau0\" <\"auto\" | tau0>> <\"tau1\" <\"auto\" | tau1>> <\"dtau\" <\"auto\" | dtau>> [\"c\" cval] [\"outfulltransform\" outdir [\"fits\"] [\"format\" format]] [\"outmaxtransform\" outdir [\"format\" format]] */
      else if(!strcmp(argv[i],"-wwz"))
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_WWZ;
	  if((c[cn].WWZ = (_WWZ *) malloc(sizeof(_WWZ))) == NULL)
	    error(ERR_MEMALLOC);
	  if(ParseWWZCommand(&i, argc, argv, p, c[cn].WWZ))
	    listcommands(argv[iterm],p);
	  cn++;
	}

      /* -Phase <"aov" | "ls" | "bls" | "fixcolumn" <colname | colnum> | "list" ["column" col] | "fix" P> ["T0" <"bls" phaseTc | "fixcolumn" <colname | colnum> | "list" ["column" col] | "fix" T0>] ["phasevar" var] ["startphase" value]*/
      else if((!strncmp(argv[i],"-Phase",6) || !strncmp(argv[i],"-phase",6)) && strlen(argv[i]) == 6)
	{
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_PHASE;
	  if((c[cn].Phase = (_Phase *) malloc(sizeof(_Phase))) == NULL)
	    error(ERR_MEMALLOC);
	  c[cn].Phase->t0type = PERTYPE_AUTOFIND;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"aov",3) && strlen(argv[i]) == 3)
		{
		  c[cn].Phase->pertype = PERTYPE_AOV;
		  m = -1;
		  /* Get the index for the last aov command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_AOV || c[l].cnum == CNUM_HARMAOV)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOAOV);
		  else
		    c[cn].Phase->lastaovindex = m;
		}
	      else if(!strncmp(argv[i],"ls",2) && strlen(argv[i]) == 2)
		{
		  c[cn].Phase->pertype = PERTYPE_LS;
		  m = -1;
		  /* Get the index for the last ls command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_LS)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOLS);
		  else
		    c[cn].Phase->lastlsindex = m;
		}
	      else if(!strncmp(argv[i],"bls",3) && strlen(argv[i]) == 3)
		{
		  c[cn].Phase->pertype = PERTYPE_BLS;
		  m = -1;
		  /* Get the index for the last bls command */
		  for(l = 0; l < cn; l++)
		    {
		      if(c[l].cnum == CNUM_BLS)
			m = l;
		    }
		  if(m < 0)
		    error(ERR_KILLHARM_NOBLS);
		  else
		    c[cn].Phase->lastblsindex = m;
		}
	      else if(!strncmp(argv[i],"fixcolumn",9) && strlen(argv[i]) == 9)
		{
		  c[cn].Phase->pertype = PERTYPE_FIXCOLUMN;
		  i++;
		  if(i < argc)
		    increaselinkedcols(p, &(c[cn].Phase->period_linkedcolumn), argv[i], cn);
		  else
		    listcommands(argv[iterm],p);
		}
	      else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
		{
		  c[cn].Phase->pertype = PERTYPE_SPECIFIED;
		  k = 0;
		  i++;
		  if(i < argc) {
		    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
		      i++;
		      if(i < argc)
			k = atoi(argv[i]);
		      else
			listcommands(argv[iterm],p);
		    } else i--;
		  } else i--;
		  RegisterDataFromInputList(p,
					    (void *) (&(c[cn].Phase->period)),
					    VARTOOLS_TYPE_DOUBLE,
					    1, cn, 0, 0, NULL, k,
					    "PHASE_PERIOD");
		}
	      else if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
		{
		  c[cn].Phase->pertype = PERTYPE_FIX;
		  i++;
		  if(i < argc)
		    {
		      c[cn].Phase->fixperiod = atof(argv[i]);
		    }
		  else
		    listcommands(argv[iterm],p);
		}
	      else
		listcommands(argv[iterm],p);
	      i++;
	      if(i < argc) {
		if(!strncmp(argv[i],"T0",2) && strlen(argv[i]) == 2) 
		  {
		    i++;
		    if(i < argc) {
		      if(!strncmp(argv[i],"bls",3) && strlen(argv[i]) == 3)
			{
			  c[cn].Phase->t0type = PERTYPE_BLS;
			  if(c[cn].Phase->pertype != PERTYPE_BLS)
			    {
			      m = -1;
			      /* Get the index for the last bls command */
			      for(l = 0; l < cn; l++)
				{
				  if(c[l].cnum == CNUM_BLS)
				    m = l;
				}
			      if(m < 0)
				error(ERR_KILLHARM_NOBLS);
			      else
				c[cn].Phase->lastblsindex = m;
			    }
			  i++;
			  if(i < argc) {
			    c[cn].Phase->phaseTc = atof(argv[i]);
			  }
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"fixcolumn",9) && strlen(argv[i]) == 9)
			{
			  c[cn].Phase->t0type = PERTYPE_FIXCOLUMN;
			  i++;
			  if(i < argc)
			    increaselinkedcols(p, &(c[cn].Phase->T0_linkedcolumn), argv[i], cn);
			  else
			    listcommands(argv[iterm],p);
			}
		      else if(!strncmp(argv[i],"list",4) && strlen(argv[i]) == 4)
			{
			  c[cn].Phase->t0type = PERTYPE_SPECIFIED;
			  k = 0;
			  i++;
			  if(i < argc) {
			    if(!strncmp(argv[i],"column",6) && strlen(argv[i]) == 6) {
			      i++;
			      if(i < argc)
				k = atoi(argv[i]);
			      else
				listcommands(argv[iterm],p);
			    } else i--;
			  } else i--;
			  RegisterDataFromInputList(p,
						    (void *) (&(c[cn].Phase->T0)),
						    VARTOOLS_TYPE_DOUBLE,
						    1, cn, 0, 0, NULL, k,
						    "PHASE_T0");
			}
		      else if(!strncmp(argv[i],"fix",3) && strlen(argv[i]) == 3)
			{
			  c[cn].Phase->t0type = PERTYPE_FIX;
			  i++;
			  if(i < argc) {
			    c[cn].Phase->fixT0 = atof(argv[i]);
			  }
			  else
			    listcommands(argv[iterm],p);
			}
		      else
			listcommands(argv[iterm],p);
		    }
		    else
		      listcommands(argv[iterm],p);
		  }
		else
		  i--;
	      }
	      else 
		i--;
	    }
	  else
	    listcommands(argv[iterm],p);
	  c[cn].Phase->phasevarname = NULL;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"phasevar")) {
	      i++;
	      if(i >= argc) listcommands(argv[iterm],p);
	      if((c[cn].Phase->phasevarname = (char *) malloc((strlen(argv[i])+1))) == NULL)
		error(ERR_MEMALLOC);
	      sprintf(c[cn].Phase->phasevarname,"%s",argv[i]);
	    } else
	      i--;
	  }
	  else
	    i--;
	  c[cn].Phase->startphase = 0.0;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"startphase")) {
	      i++;
	      if(i >= argc) listcommands(argv[iterm],p);
	      c[cn].Phase->startphase = atof(argv[i]);
	    } else
	      i--;
	  }
	  else
	    i--;
	  cn++;
	}
      
#ifdef _HAVE_PYTHON
#ifdef DYNAMICLIB
      /*
        -python < "fromfile" commandfile | commandstring >
             [ "init" < "file" initializationfile | initializationstring > |
                "continueprocess" prior_python_command_number ]
             [ "vars" variablelist | 
               [ "invars" inputvariablelist ] [ "outvars" outputvariablelist ] ]
             [ "outputcolumns" variablelist ] */
      else if(!strcmp(argv[i],"-python"))
	{
	  //if(!p->pythonlibraryloaded)
	  //  LoadVartoolsRunPythonLibrary(p);
	  iterm = i;
	  increaseNcommands(p,&c);
	  c[cn].cnum = CNUM_PYTHON;
	  c[cn].PythonCommand = CreatePythonCommandStruct(p, argv[0]);
	  i++;
	  if(i >= argc)
	    listcommands(argv[iterm],p);
	  if(ParsePythonCommand(&i, argc, argv, p, c[cn].PythonCommand, c, cn))
	    listcommands(argv[iterm],p);
	  cn++;
	}
#endif
#endif

      /* -header */
      else if(!strncmp(argv[i],"-header",7) && strlen(argv[i]) == 7)
	{
	  p->header = 1;
	}
      
      /* -headeronly */
      else if(!strncmp(argv[i],"-headeronly",11) && strlen(argv[i]) == 11)
	{
	  p->headeronly = 1;
	}
      
      /* -showinputlistformat */
      else if(!strcmp(argv[i],"-showinputlistformat"))
	{
	  p->inputlistformat = 1;
	}

      /* -showinputlcformat */
      else if(!strcmp(argv[i],"-showinputlcformat"))
	{
	  p->showinputlcformat = 1;
	}

      /* -tab */
      else if(!strncmp(argv[i],"-tab",4) && strlen(argv[i]) == 4)
	{
	  p->tabflag = 1;
	}
      
      /* -basename */
      else if(!strncmp(argv[i],"-basename",9) && strlen(argv[i]) == 9)
	{
	  p->basename = 1;
	}

      /* -functionlist */
      else if(!strcmp(argv[i],"-functionlist"))
	{
	  PrintVartoolsFunctionList(p);
	}

      /* -inputlcformat var1:col1[:type1[:fmt1][,var2:col2[:vtype2[:fmt2]],...]]
	 ["skipnum" Nskip] ["skipchar" <skipchar1[,skipchar2,...]>]
         ["delimiter" delimiter] */
      else if(!strcmp(argv[i],"-inputlcformat"))
	{
	  if(p->readformatused == 1 || p->inputlcformatused == 1) {
	    error(ERR_TOOMANYREADFORMATINPUTLCFORMAT);
	  }
	  p->inputlcformatused = 1;
	  iterm = i;
	  i++;
	  if(i >= argc)
	    help(argv[iterm],p);
	  if(ParseInputLCFormatString(argv[i], p))
	    help(argv[iterm],p);
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"skipnum")) {
	      i++;
	      if(i >= argc)
		help(argv[iterm],p);
	      p->Nskip = atoi(argv[i]);
	    } else
	      i--;
	  } else
	    i--;
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"skipchar")) {
	      i++;
	      if(i >= argc)
		help(argv[iterm],p);
	      ParseSkipCharString(argv[i],&p->Nskipchar,&p->skipchars);
	    } else {
	      p->Nskipchar = 1;
	      p->skipchars = malloc(1);
	      p->skipchars[0] = '#';
	      i--;
	    }
	  } else {
	    p->Nskipchar = 1;
	    p->skipchars = malloc(1);
	    p->skipchars[0] = '#';
	    i--;
	  }
	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"delimiter")) {
	      i++;
	      if(i >= argc)
		help(argv[iterm],p);
	      if(strlen(argv[i]) == 1) {
		p->lcdelimtype = VARTOOLS_LC_DELIMTYPE_CHAR;
		p->delimchar = argv[i][0];
	      } else {
		p->lcdelimtype = VARTOOLS_LC_DELIMTYPE_STRING;
		p->delimstring = malloc(strlen(argv[i])+1);
		sprintf(p->delimstring,"%s",argv[i]);
	      }
	    }
	    else
	      i--;
	  }
	  else
	    i--;
	}
      
      /* -inlistvars var1:col1[:vtype1[:fmt1][,var2:col2[:vtyp2[:fmt2]],....]] */
     else if(!strcmp(argv[i],"-inlistvars"))
	{
	  if(p->inlistvars == 1) {
	    error(ERR_TOOMANYINLISTVARS);
	  }
	  p->inlistvars = 1;
	  iterm = i;
	  i++;
	  if(i >= argc)
	    help(argv[iterm],p);
	  if(ParseInListVarsString(argv[i], p))
	    help(argv[iterm],p);
	}

      /* -readformat Nskip ["stringid" colstringid] ["inpututc" format] coljd colmag colsig */
      else if(!strncmp(argv[i],"-readformat",11) && strlen(argv[i]) == 11)
	{
	  if(p->readformatused == 1 || p->inputlcformatused == 1) {
	    error(ERR_TOOMANYREADFORMATINPUTLCFORMAT);
	  }
	  p->readformatused = 1;
	  iterm = i;
	  i++;
	  if(i < argc)
	    p->Nskip = atoi(argv[i]);
	  else
	    help(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"stringid",8) && strlen(argv[i]) == 8)
		{
		  p->readimagestring = 1;
		  i++;
		  if(i < argc)
		    p->colstringid = atoi(argv[i]);
		  else
		    help(argv[iterm],p);
		  i++;
		}
	    }
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"inpututc",8) && strlen(argv[i]) == 8)
		{
		  p->inputUTC = 1;
		  i++;
		  if(i < argc)
		    sprintf(p->UTCformat,"%s",argv[i]);
		  else
		    help(argv[iterm],p);
		  i++;
		}
	    }
	  if(i < argc)
	    p->coljd = atoi(argv[i]);
	  else
	    help(argv[iterm],p);
	  i++;
	  if(i < argc)
	    p->colmag = atoi(argv[i]);
	  else
	    help(argv[iterm],p);
	  i++;
	  if(i < argc)
	    p->colsig = atoi(argv[i]);
	  else
	    help(argv[iterm],p);

	  if(p->coljd < 0 || p->colmag < 0 || p->colsig < 0)
	    error(ERR_READFORMAT);
	  if(p->readimagestring)
	    {
	      if(p->colstringid <= 0)
		error(ERR_READFORMAT);
	    }

	}

      /* -readall */
      else if(!strncmp(argv[i],"-readall",8) && strlen(argv[i]) == 8)
	{
	  p->readallflag = 1;
	}

      /* -ascii */
      else if(!strncmp(argv[i],"-ascii",6) && strlen(argv[i]) == 6)
	{
	  p->ascii = 1;
	}

      else if(!strcmp(argv[i],"-binaryperiodogram"))
	{
	  p->ascii = 0;
	}

      /* -jdtol jdtol */
      else if(!strncmp(argv[i],"-jdtol",6) && strlen(argv[i]) == 6)
	{
	  iterm = i;
	  i++;
	  if(i < argc)
	    JDTOL = atof(argv[i]);
	  else
	    help(argv[iterm],p);
	}
      
      /* -matchstringid */
      else if(!strncmp(argv[i],"-matchstringid",14) && strlen(argv[i]) == 14)
	{
	  p->matchstringid = 1;
	}

      /* -nobuffer */
      else if(!strncmp(argv[i],"-nobuffer",9) && strlen(argv[i]) == 9)
	{
	  setvbuf(stdout, NULL, _IOLBF, 512 );
	}

      /* -bufferlines */
      else if(!strcmp(argv[i],"-bufferlines"))
	{
	  iterm = i;
	  i++;
	  if(i < argc)
	    p->Nbuffs_free = atoi(argv[i]);
	  else
	    help(argv[iterm],p);
	}

      /* -listcommands */
      else if(!strncmp(argv[i],"-listcommands",13) && strlen(argv[i]) == 13)
	{
	  if(i + 1 == argc)
	    listcommands(NULL,p);
	  else
	    listcommands(argv[i+1],p);
	}

      /* -help */
      else if(!strncmp(argv[i],"-help",5) && strlen(argv[i]) == 5)
	{
	  if(i + 1 == argc)
	    help(NULL, p);
	  else
	    help(argv[i+1],p);
	}

      /* -example */
      else if(!strncmp(argv[i],"-example",8) && strlen(argv[i]) == 8)
	{
	  if(i + 1 < argc)
	    example(argv[i+1],p);
	  else
	    help(argv[i],p);
	}

      /* -quiet */
      else if(!strncmp(argv[i],"-quiet",6) && strlen(argv[i]) == 6)
	{
	  p->quiet_mode=1;
	}

      /* -randseed seed */
      else if(!strncmp(argv[i],"-randseed",9) && strlen(argv[i]) == 9)
	{
	  iterm = i;
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"time",4) && strlen(argv[i]) == 4)
		{
		  p->randseed = time(NULL);
		}
	      else
		p->randseed = atoi(argv[i]);
	    }
	  else
	    help(argv[iterm],p);
	}
      
      /* -numbercolumns */
      else if(!strncmp(argv[i],"-numbercolumns",14) && strlen(argv[i]) == 14)
	{
	  p->numbercolumns=1;
	}

      /* -redirectstats statsfile [\"append\"] */
      else if(!strncmp(argv[i],"-redirectstats",14) && strlen(argv[i]) == 14)
	{
	  p->redirectstats = 1;
	  iterm = i;
	  i++;
	  if(i < argc)
	    {
	      sprintf(p->redirectstatsname,"%s",argv[i]);
	    }
	  else
	    help(argv[iterm],p);
	  i++;
	  if(i < argc)
	    {
	      if(!strncmp(argv[i],"append",6) && strlen(argv[i]) == 6)
		p->redirectstatsappend = 1;
	      else
		i--;
	    }
	  else i--;
	}

      /* -oneline */
      else if(!strncmp(argv[i],"-oneline",8) && strlen(argv[i]) == 8)
	{
	  p->oneline=1;
	}

#ifdef PARALLEL
      /* -parallel Nproc */
      else if(!strncmp(argv[i],"-parallel",9) && strlen(argv[i]) == 9)
	{
	  iterm = i;
	  i++;
	  if(i < argc)
	    {
	      p->Nproc_allow = atoi(argv[i]);
	      if(p->Nproc_allow <= 0) {
		error(ERR_NPROC_TOO_SMALL);
	      }
	    }
	  else
	    help(argv[iterm],p);
	}
#endif

#ifdef VARTOOLS_VERSION
      /* -version */
      else if(!strcmp(argv[i],"-version"))
	{
	  p->showversion = 1;
	}
#endif

      /* -log-command-line */
      else if(!strcmp(argv[i],"-log-command-line"))
	{
	  p->logcmd = 1;
	}

#ifdef DYNAMICLIB
      /* -L libraryfile */
      else if(!strcmp(argv[i],"-L"))
	{
	  iterm = i;
	  i++;
	  if(i < argc)
	    {
	      if(load_user_library(argv[i],p,0))
		error2(ERR_OPEN_LIBRARY,argv[i]);
	    }
	  else
	    help(argv[iterm],p);
	}

      /* -F libraryfile */
      else if(!strcmp(argv[i],"-F"))
	{
	  iterm = i;
	  i++;
	  if(i < argc)
	    {
	      if(load_userfunction_library(argv[i],p))
		error2(ERR_OPEN_LIBRARY,argv[i]);
	    }
	  else
	    help(argv[iterm],p);
	}

      /* -f functiondefinition */
      else if(!strcmp(argv[i],"-f"))
	{
	  iterm = i;
	  i++;
	  if(i >= argc)
	    help(argv[iterm],p);
	  ParseDefineAnalyticUserFunction(p,argv[i]);
	}

      /* User Commands */
      else if(CheckIfUserCommandIsCalled(p,c,cn,argv[i]))
	{
	  iterm = i;

	  i++;
	  ParseCL_UserCommand(p,&(c[cn]),&i,&(argv[0]),argc);

	  /* Set up any linked columns */
	  SetLinkedColumns_UserCommand(p,&(c[cn]),cn);

	  cn++;
	}
#endif

      else if(!strcmp(argv[i],"-skipmissing"))
	{
	  p->skipmissing = 1;
	}
	  
      else if(!strcmp(argv[i],"-noskipempty"))
	{
	  p->skipempty = 0;
	}

      else
	listcommands(argv[i],p);

    }

  Ncommands = p->Ncommands;

  /* Quit if there isn't anything to read in, or if there are no commands */
  if(!p->fileflag && !p->listflag && !p->headeronly && !p->inputlistformat
     && !p->showinputlcformat)
    usage(argv[0]);

  if(Ncommands < 1)
    error(ERR_USAGE);

  /* Set the default skipchar if the inputlcformat option was not given */
  if(p->Nskipchar == 0 && p->skipchars == NULL) {
    p->Nskipchar = 1;
    p->skipchars = malloc(1);
    p->skipchars[0] = '#';
  }

  /* Quit if we're matching on string id but the string id column isn't specified */
  if((p->matchstringid || p->requirestringid) && !p->readimagestring)
    error(ERR_NOSTRINGIDCOLUMN);

  /* Make sure that -readall flag is not set if there are any -copylc commands */
  if(p->readallflag && p->Ncopycommands > 0)
    error(ERR_READALL_ANDCOPYLC);

  /* If we don't have a list, make sure that there aren't any commands that require a list */
  if(p->fileflag && !p->headeronly && !p->inputlistformat 
     && !p->showinputlcformat)
    {
      for(j=0;j<Ncommands;j++)
	{
	  if(c[j].cnum == CNUM_ENSEMBLERESCALESIG)
	    error(ERR_NEEDLIST);
#ifdef _HAVE_GSL
	  else if(c[j].cnum == CNUM_ADDNOISE)
	    {
	      if(c[j].AddNoise->gammaval_type == PERTYPE_SPECIFIED ||
		 c[j].AddNoise->sig_r_type == PERTYPE_SPECIFIED ||
		 c[j].AddNoise->sig_w_type == PERTYPE_SPECIFIED)
		{
		  error(ERR_NEEDLIST);
		}
	    }
#endif
#ifdef DYNAMICLIB
	  else if(c[j].cnum == CNUM_USERCOMMAND)
	    {
	      if(c[j].UserCommand->Ninlist > 0)
		error(ERR_NEEDLIST);
	    }
#endif
	  else if(c[j].cnum == CNUM_CONVERTTIME)
	    {
	      if((c[j].ConvertTime->useradec 
		  && c[j].ConvertTime->radec_source == VARTOOLS_SOURCE_INLIST) ||
		 (c[j].ConvertTime->useppm
		  && c[j].ConvertTime->ppm_source == VARTOOLS_SOURCE_INLIST) ||
		 (c[j].ConvertTime->useinput_radec
		  && c[j].ConvertTime->inputradec_source == VARTOOLS_SOURCE_INLIST) ||
		 (c[j].ConvertTime->useinputppm
		  && c[j].ConvertTime->inputppm_source == VARTOOLS_SOURCE_INLIST))
		error(ERR_NEEDLIST);
	    }
	  else if(c[j].cnum == CNUM_FINDBLENDS)
	    error(ERR_NEEDLIST);
	  else if(c[j].cnum == CNUM_KILLHARM)
	    { 
	      if(c[j].Killharm->pertype == PERTYPE_SPECIFIED)
		error(ERR_NEEDLIST);
	    }
	  else if(c[j].cnum == CNUM_INJECTHARM)
	    {
	      if(c[j].Injectharm->pertype == PERTYPE_SPECIFIED)
		error(ERR_NEEDLIST);
	      for(l=0;l<=c[j].Injectharm->Nharm;l++)
		{
		  if(c[j].Injectharm->harm_amptype[l] == PERTYPE_SPECIFIED)
		    error(ERR_NEEDLIST);
		  if(c[j].Injectharm->harm_phasetype[l] == PERTYPE_SPECIFIED)
		    error(ERR_NEEDLIST);
		}
	      for(l=0;l<c[j].Injectharm->Nsubharm;l++)
		{
		  if(c[j].Injectharm->subharm_amptype[l] == PERTYPE_SPECIFIED)
		    error(ERR_NEEDLIST);
		  if(c[j].Injectharm->subharm_phasetype[l] == PERTYPE_SPECIFIED)
		    error(ERR_NEEDLIST);
		}
	    }
	  else if(c[j].cnum == CNUM_INJECTTRANSIT)
	    {
	      for(l=0;l<c[j].Injecttransit->Nparam;l++)
		{
		  if(c[j].Injecttransit->paramtype[l] == PERTYPE_SPECIFIED)
		    error(ERR_NEEDLIST);
		}
	    }
	  else if(c[j].cnum == CNUM_STARSPOT)
	    {
	      if(c[j].Starspot->pertype == PERTYPE_SPECIFIED)
		error(ERR_NEEDLIST);
	    }
	  else if(c[j].cnum == CNUM_PHASE)
	    {
	      if(c[j].Phase->pertype == PERTYPE_SPECIFIED || 
		 c[j].Phase->t0type == PERTYPE_SPECIFIED)
		error(ERR_NEEDLIST);
	    }
	  else if(c[j].cnum == CNUM_FIXPERBLS)
	    {
	      if(c[j].BlsFixPer->pertype == PERTYPE_SPECIFIED)
		error(ERR_NEEDLIST);
	    }
	  else if(c[j].cnum == CNUM_GETLSAMPTHRESH)
	    {
	      if(c[j].GetLSAmpThresh->pertype == PERTYPE_SPECIFIED)
		error(ERR_NEEDLIST);
	    }
	  else if(c[j].cnum == CNUM_MICROLENS)
	    {
	      if(c[j].MicroLens->f0_source == PERTYPE_SPECIFIED ||
		 c[j].MicroLens->f1_source == PERTYPE_SPECIFIED ||
		 c[j].MicroLens->u0_source == PERTYPE_SPECIFIED ||
		 c[j].MicroLens->t0_source == PERTYPE_SPECIFIED ||
		 c[j].MicroLens->tmax_source == PERTYPE_SPECIFIED)
		error(ERR_NEEDLIST);
	    }
	  else if(c[j].cnum == CNUM_DIFFFLUXTOMAG)
	    error(ERR_NEEDLIST);
	  else if(c[j].cnum == CNUM_TFA)
	    error(ERR_NEEDLIST);
	  else if(c[j].cnum == CNUM_TFA_SR)
	    error(ERR_NEEDLIST);
	  else if(c[j].cnum == CNUM_SYSREM)
	    error(ERR_NEEDLIST);
	}
      if(p->readallflag)
	error(ERR_NEEDLIST);
    }

  *cptr = c;
}
