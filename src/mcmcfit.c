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
/*                                                                            */
/* This file includes routines to run a differential evolution markov chain   */
/* monte carlo fitting procedure.                                             */
#include "commands.h"
#include "programdata.h"
#include "functions.h"
#include "mcmcfit.h"

/* This file includes functions to run a downhill simplex fit to a transit system */


void MCMC_CleanUp(_MCMC_Chain **c)
{
  int i;
  for(i=0; i < (*c)->Nlinks_allocated; i++) {
    free((*c)->links[i].p);
    if((*c)->Nauxil > 0) {
      free((*c)->links[i].auxil);
    }
  }
  free((*c)->links);
  for(i=0; i < (*c)->Nparam; i++) {
    free((*c)->sample_covar_sum2[i]);
    free((*c)->covar[i]);
  }
  free((*c)->ptrial);
  free((*c)->pinit);
  free((*c)->sample_covar_sum1);
  free((*c)->sample_covar_sum2);
  free((*c)->covar);
  free((*c)->rvec);
  free((*c)->delp);
  if((*c)->outchainfilename != NULL) {
    free((*c)->outchainfilename);
  }
  if((*c)->outchainheader != NULL) {
    free((*c)->outchainheader);
  }
  free((*c));
  *c = NULL;
}

void MCMC_GrowChain(_MCMC_Chain *c)
{
  int Nnew = 0;
  int i;
  Nnew = c->Nlinks_allocated * 2;
  if(Nnew > c->maxNalloc) Nnew = c->maxNalloc;
  if(Nnew > c->Nlinks_allocated) {
    if((c->links = (_MCMC_Link *) realloc(c->links, Nnew*sizeof(_MCMC_Link))) == NULL)
      error(ERR_MEMALLOC);
    for(i = c->Nlinks_allocated; i < Nnew; i++) {
      if((c->links[i].p = (double *) malloc(c->Nparam * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      if(c->Nauxil > 0) {
	if((c->links[i].auxil = (double *) malloc(c->Nauxil * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
      }
    }
    c->Nlinks_allocated = Nnew;
  } else {
    /* Out of space, just write over the beginning of the chain */
    c->activeindex = 0;
  }
}

void MCMC_find_initial_stepsizes(_MCMC_Chain *c);

void MCMC_initialize_chain(_MCMC_Chain *c, double eps, double maxmem, int doprintchain, char *outchainfilename, char *outchainheader, int printevery, int Nauxil, double **auxil_params, double (*funk)(double *, int, int, double *, double *, double *, void *), int N, double *t, double *mag, double *err, void *userparam) 
/* This function sets up an MCMC chain, taking as input:
   c - pointer to an _MCMC_Chain struct which has already had all of
       its parameters set using the MCMC_increment_parameter function.
   eps - a quantity used to scale random perturbations added to the 
         trial MCMC links.
   maxmem - the maximum amount of memory in GB to allocate to the chain
   doprintchain - set to 1 to output the chain to a file, to 0 to not output it.
   outchainfilename - the name of the file to output the chain to ("-" is
       stdout). This can be NULL if doprintchain is 0.
   outchainheader - an optional header to be written at the top of the output
       chain file. If this is '\0' or NULL it will be skipped.
   printevery - how often the links will be printed (1 - every link is printed;
      2 - every other link is printed; 3 - every third link printed; etc.)
   Nauxil - By default all free parameters and chi2 will be output to the chain
      file. To include additional columns specify the number (Nauxil) and 
      then provide a pointer to a vector which will store the columns. This
      vector should be defined in the userparam struct, and then filled out
      by the chi2 function.
   auxil_params - A pointer to a vector which will store the auxiliary
      parameters, and be updated after each call to chi2.
   funk - a pointer to a function which provides chi2.
       funk has the expected format:
     double funk(double *param, int Nparameters, int N, double *t, double *mag,
                  double *sig, void *userparams)
           where param is a set of parameter values to calculate chi2 for,
           and the function should return chi2.
      The other input parameters to funk are as described below.
   The following terms are passed directly to funk and not actually used
   elsewhere within this function, they can, but do not have to, 
   take on the context given below:
    N - number of points in the time series that is being fit.
    t - times of observation
    mag - set of observed magnitude values being fit
    sig - set of magnitude uncertainties
    userparams - structure containing other user-defined parameters.
         must be cast to the appropriate type by funk.
*/
{
  int nchains, i, j;
  int Nparam;

  Nparam = c->Nparam;

  nchains = 2*Nparam;
  if(nchains < MCMC_MIN_NCHAINS) nchains = MCMC_MIN_NCHAINS;
  
  c->gamma = 1.7/sqrt(2.0*((double) Nparam));
  c->eps = eps;
  c->chi2funk = funk;
  c->NJD = N;
  c->t = t;
  c->mag = mag;
  c->err = err;
  c->userparam = userparam;

  c->doprintchain = doprintchain;
  if(!doprintchain) {
    c->outchainfilename = NULL;
    c->outchainheader = NULL;
    c->printevery = 0;
    c->Nauxil = 0;
    c->auxil_params = NULL;
  } else {
    if(outchainfilename == NULL ||
       outchainfilename[0] == '\0') {
      c->outchainfilename = NULL;
      c->doprintchain = 0;
      c->outchainheader = NULL;
      c->printevery = 0;
      c->Nauxil = 0;
      c->auxil_params = NULL;
    } else {
      c->outchainfilename = (char *) malloc((strlen(outchainfilename)+1));
      sprintf(c->outchainfilename,"%s",outchainfilename);
      if(outchainheader == NULL || outchainheader[0] == '\0') {
	c->outchainheader = NULL;
      } else {
	c->outchainheader = (char *) malloc((strlen(outchainheader)+1));
	sprintf(c->outchainheader,"%s",outchainheader);
      }
      c->printevery = printevery;
      c->Nauxil = Nauxil;
      c->auxil_params = auxil_params;
    }
  }

  c->Nchains_running = nchains;
  c->Ntrials_accepted = 0;
  c->Nlinks_run = 0;
  c->activeindex = nchains;
  c->chainindex = 0;
  c->Nlinks_covar = 0;

  c->maxNalloc = floor(maxmem*1.0e9/((Nparam + c->Nauxil + 1)*sizeof(double)));
  
  c->Nlinks_allocated = (MCMC_DEFAULT_INITIAL_CHAINALLOCATION < c->maxNalloc ?
                         MCMC_DEFAULT_INITIAL_CHAINALLOCATION :
                         c->maxNalloc);

  
  if((c->links = (_MCMC_Link *) malloc(c->Nlinks_allocated*(sizeof(_MCMC_Link)))) == NULL ||
     (c->ptrial = (double *) malloc(Nparam * sizeof(double))) == NULL ||
     (c->sample_covar_sum1 = (double *) malloc(Nparam * sizeof(double))) == NULL ||
     (c->sample_covar_sum2 = (double **) malloc(Nparam * sizeof(double *))) == NULL ||
     (c->covar = (double **) malloc(Nparam * sizeof(double *))) == NULL ||
     (c->rvec = (double *) malloc(Nparam * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  
  for(i=0; i < c->Nlinks_allocated; i++) {
    if((c->links[i].p = (double *) malloc(Nparam*sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
    if(c->Nauxil > 0) {
      if((c->links[i].auxil = (double *) malloc(c->Nauxil * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
  }
  for(i=0; i < Nparam; i++) {
    if((c->sample_covar_sum2[i] = (double *) malloc(Nparam * sizeof(double))) == NULL ||
       (c->covar[i] = (double *) malloc(Nparam * sizeof(double))) == NULL) {
      error(ERR_MEMALLOC);
    }
  }

  MCMC_find_initial_stepsizes(c);

  for(j=0; j < Nparam; j++) {
    c->links[0].p[j] = c->pinit[j];
  }
  c->links[0].L = (*funk)(c->links[0].p, Nparam, N, t, mag, err, userparam);
  for(j=0; j < c->Nauxil; j++) {
    c->links[0].auxil[j] = (*(c->auxil_params))[j];
  }

  for(i=1; i < nchains; i++) {
    for(j=0; j < Nparam; j++) {
      c->links[i].p[j] = c->pinit[j] + c->delp[j]*gasdev();
    }
    c->links[i].L = (*funk)(c->links[i].p, Nparam, N, t, mag, err, userparam);
    for(j=0; j < c->Nauxil; j++) {
      c->links[0].auxil[j] = (*(c->auxil_params))[j];
    }
  }
  for(i = 0; i < Nparam; i++) {
    c->sample_covar_sum1[i] = 0.;
    for(j = 0; j < Nparam; j++) {
      c->covar[i][j] = 0.;
      c->sample_covar_sum2[i][j] = 0.;
    }
    c->covar[i][i] = c->delp[i]*c->delp[i];
  }
  choldc(c->covar, Nparam);
}

void MCMC_DifferentialEvolution_GenerateProposal(_MCMC_Chain *c){
  int i, i0, itest;
  int j, k;
  
  itest = c->activeindex - c->Nchains_running;
  if(itest < 0) {
    itest = c->Nlinks_allocated + itest;
  }
  i0 = c->activeindex - c->chainindex - c->Nchains_running;
  if(i0 < 0) {
    i0 = c->Nlinks_allocated + i0;
  }

  j = rand() % (c->Nchains_running - 1);
  k = rand() % (c->Nchains_running - 2);
  if(j >= c->chainindex) j++;
  if(k >= (j < c->chainindex ? j : c->chainindex)) k++;
  if(k >= (j < c->chainindex ? c->chainindex : j)) k++;

  j = i0 + j;
  if(j >= c->Nlinks_allocated) { j -= c->Nlinks_allocated; }
  k = i0 + k;
  if(k >= c->Nlinks_allocated) { k -= c->Nlinks_allocated; }

  for(i=0; i < c->Nparam; i++) {
    c->rvec[i] = c->eps*gasdev();
  }

  for(i=0; i < c->Nparam; i++) {
    c->ptrial[i] = c->links[itest].p[i] + c->gamma*(c->links[j].p[i] - c->links[k].p[i]);
    for(j=0; j <= i; j++) {
      c->ptrial[i] += c->covar[i][j]*c->rvec[j];
    }
  }
}

void MCMC_DifferentialEvolution_NextLink(_MCMC_Chain *c){
  
  int i, j, k;
  double chi2new;
  double chi2old;
  int itest;
  double u, alpha, vx;
  int test = 0;


  if(c->activeindex >= c->Nlinks_allocated) {
    MCMC_GrowChain(c);
  }

  MCMC_DifferentialEvolution_GenerateProposal(c);

  itest = c->activeindex - c->Nchains_running;
  if(itest < 0) {
    itest = c->Nlinks_allocated + itest;
  }
  
  chi2old = c->links[itest].L;

  chi2new = (*c->chi2funk)(c->ptrial, c->Nparam, c->NJD, c->t, c->mag, c->err, c->userparam);

  j = c->activeindex;

  /* proposal is in an invalid region, accept the old proposal */
  if(isnan(chi2new)) {
    for(i=0; i < c->Nparam; i++) {
      c->links[j].p[i] = c->links[itest].p[i];
    }
    c->links[j].L = c->links[itest].L;
    for(i=0; i < c->Nauxil; i++) {
      c->links[j].auxil[i] = c->links[itest].auxil[i];
    }
  }
  else {
    vx = exp(-0.5*(chi2new - chi2old));
    alpha = (vx < 1.0 ? vx : 1.0);
    u = (rand() / (RAND_MAX + 0.0));
    if(alpha > 0.0 && u <= alpha) {
      /* the new proposal is accepted */
      for(i=0; i < c->Nparam; i++) {
	c->links[j].p[i] = c->ptrial[i];
      }
      c->links[j].L = chi2new;
      for(i=0; i < c->Nauxil; i++) {
	c->links[j].auxil[i] = (*(c->auxil_params))[i];
      }
      c->Ntrials_accepted += 1;
    }
    else {
      /* The new proposal is rejected */
      for(i=0; i < c->Nparam; i++) {
	c->links[j].p[i] = c->links[itest].p[i];
      }
      c->links[j].L = c->links[itest].L;
      for(i=0; i < c->Nauxil; i++) {
	c->links[j].auxil[i] = c->links[itest].auxil[i];
      }
    }
  }


  test = 0;
  for(i=0; i < c->Nparam; i++) {
    c->sample_covar_sum1[i] += (c->pinit[i] - c->links[j].p[i]);
    if(c->sample_covar_sum1[i]/(c->Nlinks_covar + 1) >= 0.1*DBL_MAX ||
       c->sample_covar_sum1[i]/(c->Nlinks_covar + 1) <= -0.1*DBL_MAX) {
      test = 1;
    }
    for(k=0; k < c->Nparam; k++) {
      c->sample_covar_sum2[i][k] += (c->pinit[i] - c->links[j].p[i]) * (c->pinit[k] - c->links[j].p[k]);
      if(c->sample_covar_sum2[i][k]/(c->Nlinks_covar + 1) >= 0.1*DBL_MAX ||
	 c->sample_covar_sum2[i][k]/(c->Nlinks_covar + 1) <= -0.1*DBL_MAX) {
	test = 1;
      }
    }
  }
  if(test) {
    j = c->activeindex - c->Nlinks_covar;
    if(j < 0) {
      j += c->Nlinks_allocated;
    }
  }
  else {
    if(c->Nlinks_run >= c->Nlinks_allocated) {
      if(c->activeindex < c->Nlinks_run) {
	j = c->activeindex + 1;
	if(j >= c->Nlinks_allocated) j = 0;
      }
      else {
	j = 0;
      }
      test = 1;
    }
  }
  if(test) {
    for(i=0; i < c->Nparam; i++) {
      c->sample_covar_sum1[i] -= (c->pinit[i] - c->links[j].p[i]);
      for(k=0; k < c->Nparam; k++) {
	c->sample_covar_sum2[i][k] -= (c->pinit[i] - c->links[j].p[i]) * (c->pinit[k] - c->links[j].p[k]);
      }
    }
  }

  if(!test)
    c->Nlinks_covar += 1;
  c->Nlinks_run += 1;
  c->activeindex += 1;
  c->chainindex += 1;
  if(c->chainindex >= c->Nchains_running) c->chainindex = 0;
}

void MCMC_get_covarmatrix(_MCMC_Chain *c)
{
  int j, k;

  for(j=0; j < c->Nparam; j++) {
    for(k = 0; k < c->Nparam; k++) {
      c->covar[j][k] = (c->sample_covar_sum2[j][k]/(c->Nlinks_covar+1)) - (c->sample_covar_sum1[j]/(c->Nlinks_covar+1))*(c->sample_covar_sum1[k]/(c->Nlinks_covar+1));
    }
  }
  choldc(c->covar, c->Nparam);
}

void MCMC_DifferentialEvolution_RunMCMC(int maxNtrials_accepted, int maxNlinks_run, int N_recalc_covar, _MCMC_Chain *c)
{
  FILE *outchainfile;
  int j, i;
  if(c->doprintchain) {
    if((outchainfile = fopen(c->outchainfilename,"w")) == NULL)
      error2(ERR_CANNOTWRITE,c->outchainfilename);
    if(c->outchainheader != NULL) {
      if(c->outchainheader[0] != '\0') {
	fprintf(outchainfile,"%s\n",c->outchainheader);
      }
    }
  }
  do {
    if(c->Nlinks_run % N_recalc_covar == N_recalc_covar - 1) {
      MCMC_get_covarmatrix(c);
    }
    MCMC_DifferentialEvolution_NextLink(c);
    if(c->doprintchain) {
      if(c->Nlinks_run % c->printevery == 0) {
	j = c->activeindex - 1;
	if(j < 0) {
	  j = c->Nlinks_allocated - 1;
	}
	for(i = 0; i < c->Nparam; i++) {
	  if(i)
	    fprintf(outchainfile," ");
	  fprintf(outchainfile,"%.17g",c->links[j].p[i]);
	}
	for(i = 0; i < c->Nauxil; i++) {
	  fprintf(outchainfile," %.17g",c->links[j].auxil[i]);
	}
	fprintf(outchainfile," %.17g\n",c->links[j].L);
      }
    }
  } while ((maxNtrials_accepted < 0 || 
	    c->Ntrials_accepted < maxNtrials_accepted)
	   &&
	   (maxNlinks_run < 0 ||
	    c->Nlinks_run < maxNlinks_run));
  if(c->doprintchain) {
    fclose(outchainfile);
  }
}

int MCMC_GetBestLink(_MCMC_Chain *c)
/* Returns the index of the lowest chi2 link in the Markov chain */
{
  int i, ibest, maxN;
  double Lbest;
  maxN = (c->Nlinks_run < c->Nlinks_allocated ? c->Nlinks_run 
	  : c->Nlinks_allocated);
  if(maxN <= 0) {
    return -1;
  }
  ibest = 0;
  Lbest = c->links[0].L;
  for(i=1; i < maxN; i++) {
    if(c->links[i].L < Lbest) {
      ibest = i;
      Lbest = c->links[i].L;
    }
  }
  return ibest;
}

void MCMC_find_initial_stepsizes(_MCMC_Chain *c)
/* Adjust the initial step sizes such that delta-chi2 = 1 for each one */
{
  int i, j, Ntest;
  double chi20, chi2test;
  for(i=0; i < c->Nparam; i++) {
    for(j=0; j < c->Nparam; j++) {
      c->ptrial[j] = c->pinit[j];
    }
    if(!i) {
      chi20 = (*c->chi2funk)(c->ptrial, c->Nparam, c->NJD, c->t, c->mag, c->err, c->userparam);
    }
    if(isnan(chi20))
      return;
    if(!c->delp[i]) {
      if(c->pinit[i] != 0) {
	c->delp[i] = 0.001*c->pinit[i];
      }
      else
	c->delp[i] = 0.001;
    }
    c->ptrial[i] = c->pinit[i] + c->delp[i];
    Ntest = 0;
    chi2test = (*c->chi2funk)(c->ptrial, c->Nparam, c->NJD, 
			      c->t, c->mag, c->err,
			      c->userparam);
    if(chi2test == chi20)
      continue;
    while(Ntest < 100 ? (isnan(chi2test) ? 1 : (fabs(chi2test - chi20) < 0.99 ||
						fabs(chi2test - chi20) > 1.01))
	  : 0) {
      while(Ntest < 100 && isnan(chi2test)) {
	c->delp[i] = -c->delp[i];
	c->ptrial[i] = c->pinit[i] + c->delp[i];
	chi2test = (*c->chi2funk)(c->ptrial, c->Nparam, c->NJD, 
				  c->t, c->mag, c->err,
				  c->userparam);
	Ntest++;
	if(!isnan(chi2test)) break;
	c->delp[i] = c->delp[i]/10.;
	c->ptrial[i] = c->pinit[i] + c->delp[i];
	chi2test = (*c->chi2funk)(c->ptrial, c->Nparam, c->NJD, 
				  c->t, c->mag, c->err,
				  c->userparam);
	Ntest++;
      }
      if(chi2test == chi20) break;
      c->delp[i] = c->delp[i]/(sqrt(fabs(chi2test - chi20)));
      c->ptrial[i] = c->pinit[i] + c->delp[i];
      chi2test = (*c->chi2funk)(c->ptrial, c->Nparam, c->NJD, 
				c->t, c->mag, c->err,
				c->userparam);
      Ntest++;
    }
  }
}

void MCMC_increment_parameter(_MCMC_Chain **cptr, double value, double delp)
/* This functions adds a new parameter to an MCMC chain.
   cptr - is a pointer to a pointer storing the chain. The first time this is
          called, *cptr should be NULL. MCMC_increment_parameter will
          allocate memory for the object.
   value - the initial value for the new parameter.
   delp - the initial estimated uncertainty for the parameter.

   For example, to use this function you might have:
   _MCMC_Chain *chain = NULL;
   MCMC_increment_parameter(&chain, 1.0, 1.0);
*/
{
  _MCMC_Chain *c;
  if(*cptr == NULL) {
    if(((*cptr) = (_MCMC_Chain *) malloc(sizeof(_MCMC_Chain))) == NULL)
      error(ERR_MEMALLOC);
    (*cptr)->delp = NULL;
    (*cptr)->ptrial = NULL;
    (*cptr)->pinit = NULL;
    (*cptr)->links = NULL;
    (*cptr)->rvec = NULL;
    (*cptr)->sample_covar_sum1 = NULL;
    (*cptr)->sample_covar_sum2 = NULL;
    (*cptr)->covar = NULL;
    (*cptr)->Nparam = 0;
    (*cptr)->Ntrials_accepted = 0;
    (*cptr)->Nlinks_run = 0;
    (*cptr)->Nlinks_allocated = 0;
    (*cptr)->activeindex = -1;
    (*cptr)->Nchains_running = 0;
  }
  c = *cptr;
  if(!c->Nparam) {
    if((c->delp = (double *) malloc(sizeof(double))) == NULL ||
       (c->pinit = (double *) malloc(sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
  } else {
    if((c->delp = (double *) realloc(c->delp, (c->Nparam + 1)*sizeof(double))) == NULL ||
       (c->pinit = (double *) realloc(c->pinit, (c->Nparam + 1)*sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
  }
  c->pinit[c->Nparam] = value;
  c->delp[c->Nparam] = delp;
  c->Nparam += 1;
}

