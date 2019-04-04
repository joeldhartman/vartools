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
/* This file holds the main function of the VARTOOLS program which is where
   execution of the code begins. Most of the work is carried out by functions
   found in other files. The parallel process handling is implemented in this
   file.
*/

#include "commands.h"
#include "programdata.h"
#include "functions.h"

#ifdef PARALLEL
#include <pthread.h>
#include <semaphore.h>

typedef struct {
  ProgramData *p;
  Command *c;
  FILE *outfile;
  int lcnumber;
  int threadnumber;
  int finished;
  sem_t newthreaddata;
} _threaddata;

void *ParallelProcessOneLC(void *arg)
{
  _threaddata *t;
  int i, k;
  int readlc_retval, cnum_start;
  _StringBuffer *buf = NULL;

  t = (_threaddata *) arg;

  while(!t->finished) {
    sem_wait(&(t->newthreaddata));
    while(pthread_mutex_trylock(&(t->p->Nproc_mutex)));
    if(t->finished) {
      pthread_mutex_unlock(&(t->p->Nproc_mutex));
      break;
    }
    pthread_mutex_unlock(&(t->p->Nproc_mutex));

    if(t->p->Ncopycommands > 0) {
      if(t->p->start_cnum[t->lcnumber] <= 0) {
	readlc_retval = ReadSingleLightCurve(t->p,t->c,t->lcnumber,t->threadnumber);
	cnum_start = 0;
      } else {
	getlccopy(t->p, t->c, t->threadnumber, t->lcnumber);
	cnum_start = t->p->start_cnum[t->lcnumber];
	readlc_retval = 0;
      }
    } else {
      readlc_retval = ReadSingleLightCurve(t->p,t->c,t->lcnumber,t->threadnumber);
      cnum_start = 0;
    }

    if(readlc_retval ||
       (t->p->NJD[t->threadnumber] == 0 &&
	t->p->skipempty)) {
      /* Skipping this light curve, indicate that this thread is free */
      if(t->p->Ncopycommands > 0) {
	/* Turn off any copies which will depend on this light curve */
	turnoffcopies(t->p, t->c, cnum_start, t->threadnumber, t->lcnumber);
      }
      while(pthread_mutex_trylock(&(t->p->Nproc_mutex)));
      t->p->threadsinuse[t->threadnumber] = 0;
      //t->p->Nproc = t->p->Nproc - 1;
      pthread_mutex_unlock(&(t->p->Nproc_mutex));
      sem_post(&(t->p->threadfree));
      continue;
    }
    if(t->p->basename)
      Switchtobasename(t->p,t->lcnumber);
    if(cnum_start == 0)
      SetTimeMagSigPointers(t->p,t->threadnumber);
    if(t->p->isifcommands) {
      if(cnum_start == 0) {
	while(t->p->IfStack[t->threadnumber]->curpos > 0)
	  popIfStack(t->p->IfStack[t->threadnumber]);
	t->p->IfStack[t->threadnumber]->curpos = 0;
      }
    }
    for(i=cnum_start;i<t->p->Ncommands;i++)
      ProcessCommandSingle(t->p,&(t->c[i]),t->lcnumber,i,t->threadnumber);
    if(!t->p->quiet_mode)
      {

	while(pthread_mutex_trylock(&(t->p->outfile_mutex)));
	if(t->p->Nbuffs_free > 0) {
	  buf = popFreeBuffer(t->p);
	  pthread_mutex_unlock(&(t->p->outfile_mutex));
	} else {
	  emptyresults_buffer(t->p,t->outfile);
	  buf = popFreeBuffer(t->p);
	  pthread_mutex_unlock(&(t->p->outfile_mutex));
	  if(buf == NULL) {
	    do {
	      while(pthread_mutex_trylock(&(t->p->outfile_mutex)));
	      emptyresults_buffer(t->p,t->outfile);
	      buf = popFreeBuffer(t->p);
	      pthread_mutex_unlock(&(t->p->outfile_mutex));
	    } while(buf == NULL);
	  }
	}
	printresults_buffer_new(t->p, t->threadnumber, t->lcnumber,&(buf->data),&(buf->len));
	while(pthread_mutex_trylock(&(t->p->outfile_mutex)));
	pushFullBuffer(t->p,buf);
	pthread_mutex_unlock(&(t->p->outfile_mutex));
	
	//printresults(p.tabflag,p.lcnames[j],c,p.Ncommands,0,j);
	/*while(pthread_mutex_trylock(&(t->p->outfile_mutex)));
	printresults_new(t->p, t->threadnumber, t->lcnumber,t->outfile);
	pthread_mutex_unlock(&(t->p->outfile_mutex));*/
      }
    while(pthread_mutex_trylock(&(t->p->Nproc_mutex)));
    t->p->threadsinuse[t->threadnumber] = 0;
    //t->p->Nproc = t->p->Nproc - 1;
    pthread_mutex_unlock(&(t->p->Nproc_mutex));
    sem_post(&(t->p->threadfree));
  }
  return NULL;
}
#endif

int main(int argc, char **argv)
{
  ProgramData p;
  Command *c;

  FILE *outfile;

  int i, j, k, kk;

  int cnum_start, readlc_retval;

  struct ListNode *testnode, *testnode2;

#ifdef PARALLEL
  _threaddata *threaddata;
  int checkthreadsout;
#endif

  /* Initialize the program variables */
  p.Nlcs = 0;
  p.Ncommands = 0;
  p.NJD = 0;
  p.Nthread = 0;
  p.readallflag = 0;
  p.listflag = 0;
  p.fileflag = 0;
  p.tabflag = 0;
  p.header = 0;
  p.headeronly = 0;
  p.inputlistformat = 0;
  p.showinputlcformat = 0;
  p.basename = 0;
  p.Nskip = 0;
  p.readimagestring = 0;
  p.colstringid = 0;
  p.matchstringid = 0;
  p.requirestringid = 0;
  p.inputUTC = 0;
  p.coljd = 1;
  p.colmag = 2;
  p.colsig = 3;
  p.decorrflag = 0;
  p.ascii = 1;
  p.quiet_mode = 0;
  p.sizecommandvector = INITCOMMANDSIZE;
  p.randseed = 1;
  p.numbercolumns = 0;
  p.oneline = 0;
  p.Ncolstolink = 0;
  p.redirectstats = 0;
  p.redirectstatsappend = 0;
  p.maxinputcolumn = 0;
  p.inputcolumn_iter_index = 1;
  p.DataFromInputList = NULL;
  p.NDataFromInputList = 0;
  p.DefinedVariables = NULL;
  p.NDefinedVariables = 0;
  p.NScalarData = 0;
  p.ScalarData = NULL;
  p.NInternalVars = 0;
  p.ismultifilt = 0;
  p.isifcommands = 0;
  p.max_colcommand = -1;
  p.size_colcommandvec = -1;
  p.showversion = 0;
  p.logcmd = 0;
  p.storecmd = 0;
  p.skipmissing = 0;
  p.skipempty = 1;

  p.NDataFromLightCurve = 0;
  p.maxinputlccolumn = 0;
  p.DataFromLightCurve = NULL;

  p.use_lc_open_exec_command = 0;

  p.lcdelimtype = VARTOOLS_LC_DELIMTYPE_WHITESPACE;

#ifdef _USEBINARY_LC
  p.binarylcinput = 0;
#endif

  p.lc_getcolumnsfromheader = 0;
  p.lc_getcolumnsfromheader_notyetset = 1;

  //  sizeHISTvector = 0;
  JDTOL = DEFAULT_JDTOL;

#ifdef PARALLEL
  p.Nproc_allow = 1;
#endif
  p.Nbuffs_free = VARTOOLS_DEFAULT_NOUTPUT_BUFFERS;

#ifdef DYNAMICLIB
  p.NUserLib = 0;
  p.NUserFunc = 0;
  p.UserFunc = NULL;
  p.NAnalyticUserFunc = 0;
  p.AnalyticUserFunc = NULL;
#endif

  InitLinkedList(&(p.lcs_to_proc));

  if((c = (Command *) malloc(p.sizecommandvector * sizeof(Command))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0; i < p.sizecommandvector; i++) {
    c[i].require_sort = 0;
    c[i].require_distinct = 0;
    c[i].N_setparam_expr = 0;
  }

#ifdef DYNAMICLIB

  if(lt_dlinit()) {
    error2(ERR_OPEN_LIBRARY,"Cannot initialize libtool for opening a library.\n");
  }
#ifdef VARTOOLSLIB_USERLIBSDIR
  if(lt_dladdsearchdir(VARTOOLSLIB_USERLIBSDIR)) {
    lt_dlerror();
    exit(1);
  }
#endif
#ifdef VARTOOLSLIB_USERFUNCSDIR
  if(lt_dladdsearchdir(VARTOOLSLIB_USERFUNCSDIR)) {
    lt_dlerror();
    exit(1);
  }
#endif

#endif

  /* Parse the command line and initialize the necessary variables */

  parsecommandline(argc, argv, &p, &c);

  /* Seed the random number generator */
  srand(p.randseed);

  /* Read in the input list, initialize variables, set up the output columns,
     compile expressions, and allocate memory */
  InitCommands(&p, c);

  /* Prepare the output columns, and link them to command inputs if needed */
  /*
  p.Ncolumns = 0;
  CreateOutputColumns(&p, c, p.Ncommands);
  linkcolumns(&p);
  */

  /* Get the file to output the statistics to */
  if(!p.redirectstats)
    {
      outfile = stdout;
    }
  else
    {
      if(!p.redirectstatsappend)
	{
	  if((outfile = fopen(p.redirectstatsname,"w")) == NULL)
	    error2(ERR_CANNOTWRITE,p.redirectstatsname);
	}
      else
	{
	  if((outfile = fopen(p.redirectstatsname,"a")) == NULL)
	    error2(ERR_CANNOTWRITE,p.redirectstatsname);
	}
    }

  /*Output the version and or the command-line if asked to */
#ifdef VARTOOLS_VERSION
  if(p.showversion) {
    fprintf(outfile,"# Output produced by VARTOOLS version %s\n",VARTOOLS_VERSION);
  }
#endif
  if(p.logcmd) {
    fprintf(outfile,"#");
    for(i=0; i < argc; i++) {
      k = 0;
      j = 0;
      while(argv[i][j] != '\0') {
	if(check_isspecialchar(argv[i][j])){
	  k = 1;
	  break;
	}
	j++;
      }
      if(k) {
	fprintf(outfile," \'%s\'", argv[i]);
      }
      else {
	fprintf(outfile," %s", argv[i]);
      }
    }
    fprintf(outfile,"\n");
  }
  if(p.storecmd) {
    p.sizecmdline = 0;
    for(i=0; i < argc; i++) {
      k = 0;
      j = 0;
      while(argv[i][j] != '\0') {
	if(check_isspecialchar(argv[i][j])){
	  k = 1;
	  break;
	}
	j++;
      }
      if(k) {
	p.sizecmdline += (3 + strlen(argv[i]));
      }
      else {
	p.sizecmdline += (1 + strlen(argv[i]));
      }
    }
    if((p.cmdline = (char *) malloc((p.sizecmdline + 1))) == NULL)
      error(ERR_MEMALLOC);
    kk = 0;
    for(i=0; i < argc; i++) {
      k = 0;
      j = 0;
      while(argv[i][j] != '\0') {
	if(check_isspecialchar(argv[i][j])){
	  k = 1;
	  break;
	}
	j++;
      }
      if(k) {
	sprintf(&(p.cmdline[kk])," \'%s\'", argv[i]);
	kk += (3 + strlen(argv[i]));
      }
      else {
	sprintf(&(p.cmdline[kk])," %s", argv[i]);
	kk += (1 + strlen(argv[i]));
      }
    }
  }    

  /* Print out the header if we're supposed to. */
  if(p.header || p.tabflag || p.headeronly)
    //printheader(p.tabflag,c,p.Ncommands,p.numbercolumns);
    printheader_new(&p,outfile);

  /* Print out the expected input list format */
  if(p.inputlistformat)
    printinputlistformat(&p,outfile);

  /* Print out the expected input lc format */
  if(p.inputlistformat)
    fprintf(outfile,"\n");
  if(p.showinputlcformat)
    printinputlcformat(&p,outfile);

  /* Now quit if all we wanted to do was print out the header or the
     input list format */
  if(p.headeronly || p.inputlistformat || p.showinputlcformat)
    exit(0);

  /* Determine the input columns */
  //DetermineColumns(&p, c);

  /* Compile any analytic expressions which will need to be evaluated */
  //CompileAllExpressions(&p, c);


  /* Read the global files for decorrelation */
  ReadGlobalDecorr(&p, c);

  /* Read in the dates files for Jstet */
  ReadDatesFiles(&p, c);

  /* The program will be run differently depending on whether or not we're reading in all the light curves at once */
  if(p.readallflag)
    {
      InitializeMemAllocDataFromLightCurve(&p, p.Nlcs);
      /* Read in the light curves */
      if(ReadAllLightCurves(&p, c)) {
	RemoveEmptyLightCurves(&p, c);
      }

      /* If the basename switch was set, then remove the directory names from all the light curves */
      if(p.basename)
	Switchtobasename(&p,-1);

      /* Process each command one at a time */
      for(j=0;j<p.Nlcs;j++)
	SetTimeMagSigPointers(&p,j);
      for(i=0;i<p.Ncommands;i++)
	ProcessCommandAll(&p,&c[i],i);

      /* Print out the results for all the stars */
      if(!p.quiet_mode)
	{
	  for(j=0;j<p.Nlcs;j++)
	    //printresults(p.tabflag,p.lcnames[j],c,p.Ncommands,j,j);
	    printresults_new(&p, j, j, outfile);
	}
    }
  else
    {
#ifdef PARALLEL
      if(p.Nproc_allow > 1) {
	/*Cycle through the light curves, processing each command one at a time,
	  Light curves are processed in parallel in this branch */
	if((threaddata = (_threaddata *) malloc(p.Nproc_allow*sizeof(_threaddata))) == NULL)
	  error(ERR_MEMALLOC);
	sem_init(&p.threadfree, 0, p.Nproc_allow);
	for(k=0;k<p.Nproc_allow;k++) {
	  threaddata[k].p = &p;
	  threaddata[k].c = c;
	  threaddata[k].outfile = outfile;
	  threaddata[k].finished = 0;
	}
	InitializeMemAllocDataFromLightCurve(&p, p.Nproc_allow);
	InitializeOutputBufferStacks(&p);
	p.Nproc = 0;
	for(j=0; j < p.Nproc_allow; j++)
	  p.pth_init[j] = 0;
	for(j=p.Nlcs-1;j>=0;j--)
	  PushNode(&(p.lcs_to_proc),j);
	while(p.lcs_to_proc.Nnodes > 0)
	  {
	    sem_wait(&(p.threadfree));
	    while(pthread_mutex_trylock(&(p.Nproc_mutex)));
	    for(k=0; k < p.Nproc_allow; k++) {
	      if(p.threadsinuse[k] <= 0)
		break;
	    }
	    /* Find an lc that can be processed */
	    testnode = p.lcs_to_proc.lastnode;
	    while(pthread_mutex_trylock(&(p.is_lc_ready_mutex)));
	    while(1) {
	      j = testnode->val;
	      if(p.Ncopycommands > 0 ? p.is_lc_ready[j] > 0 : 1) {
		RemoveNode(&(p.lcs_to_proc),testnode);
		break;
	      } else {
		testnode2 = testnode;
		testnode = testnode->priornode;
		if(p.is_lc_ready[j] < 0)
		  RemoveNode(&(p.lcs_to_proc),testnode2);
		if(testnode == NULL) {
		  j = -1;
		  break;
		}
	      }
	    }
	    pthread_mutex_unlock(&(p.is_lc_ready_mutex));
	    if(j >= 0) {
	      if(!p.pth_init[k]) {
		p.Nproc += 1;
		threaddata[k].lcnumber = j;
		p.threadsinuse[k] = 1;
		threaddata[k].threadnumber = k;
		sem_init(&(threaddata[k].newthreaddata), 0, 1);
		pthread_create(&(p.pth[k]), NULL, ParallelProcessOneLC, &(threaddata[k]));
		p.pth_init[k] = 1;
	      } else {
		threaddata[k].lcnumber = j;
		p.threadsinuse[k] = 1;
		threaddata[k].threadnumber = k;
		sem_post(&(threaddata[k].newthreaddata));
	      }
	    } else {
	      sem_post(&(p.threadfree));
	    }
	    pthread_mutex_unlock(&(p.Nproc_mutex));
	  }
	do {
	  while(pthread_mutex_trylock(&(p.Nproc_mutex)));
	  checkthreadsout = 0;
	  for(k=0; k < p.Nproc; k++) {
	    if(p.threadsinuse[k])
	      checkthreadsout = 1;
	    else if(!threaddata[k].finished) {
	      threaddata[k].finished = 1;
	      sem_post(&(threaddata[k].newthreaddata));
	    }
	  }
	  pthread_mutex_unlock(&(p.Nproc_mutex));
	} while(checkthreadsout);
	for(j=0; j < p.Nproc_allow; j++)
	  {
	    if(p.pth_init[j])
	      {
		pthread_join(p.pth[j],NULL);
	      }
	  }
	if(p.Nbuffs_full > 0) {
	  emptyresults_buffer(&p,outfile);
	}
      } else {
#endif
      /*Cycle through the light curves, processing each command one at a time */
      /* This branch is if we do not compile with parallel-processing enabled */
      /* or if running using only 1 process */
      InitializeMemAllocDataFromLightCurve(&p, 1);
      for(j=0;j<p.Nlcs;j++)
	{
	  if(p.Ncopycommands > 0) {
	    if(p.is_lc_ready[j] < 0)
	      continue;
	    if(p.start_cnum[j] <= 0) {
	      readlc_retval = ReadSingleLightCurve(&p,c,j,0);
	      cnum_start = 0;
	    } else {
	      getlccopy(&p, c, 0, j);
	      cnum_start = p.start_cnum[j];
	      readlc_retval = 0;
	    }
	  } else {
	    readlc_retval = ReadSingleLightCurve(&p,c,j,0);
	    cnum_start = 0;
	  }
	  if(readlc_retval ||
	     (p.NJD[0] == 0 && p.skipempty)) {
	    /* Skip this light curve */
	    if(p.Ncopycommands > 0) {
	      turnoffcopies(&p, c, cnum_start, 0, j);
	    }
	    continue;
	  }
	  if(p.basename)
	    Switchtobasename(&p,j);
	  if(cnum_start == 0)
	    SetTimeMagSigPointers(&p, 0);
	  if(p.isifcommands) {
	    if(cnum_start == 0) {
	      while(p.IfStack[0]->curpos > 0)
		popIfStack(p.IfStack[0]);
	      p.IfStack[0]->curpos = 0;
	    }
	  }
	  for(i=cnum_start;i<p.Ncommands;i++)
	    ProcessCommandSingle(&p,&c[i],j,i,0);
	  if(!p.quiet_mode)
	    //printresults(p.tabflag,p.lcnames[j],c,p.Ncommands,0,j);
	    printresults_new(&p, 0, j,outfile);
	}
#ifdef PARALLEL
      }
#endif
    }
  if(outfile != stdout)
    fclose(outfile);
#ifdef _HAVE_PYTHON
  KillAllPythonProcesses(&p, &c);
#endif
  exit(0);
}








