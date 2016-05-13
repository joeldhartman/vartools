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
#ifndef _MCMC_HEADER_DEFS
#define _MCMC_HEADER_DEFS

#include <float.h>

#define MCMC_DEFAULT_INITIAL_CHAINALLOCATION 16384
#define MCMC_DEFAULT_RECALC_COVAR 10000
#define MCMC_MIN_NCHAINS 3

typedef struct {
  double *p;
  double *auxil;
  double L;
} _MCMC_Link;
typedef struct {
  int Nparam;
  int Ntrials_accepted;
  int Nlinks_covar;
  int Nlinks_run;
  int Nlinks_allocated;
  int chainindex;
  int activeindex;
  int Nchains_running;
  double gamma;
  double eps;
  double *delp;
  double *ptrial;
  double *pinit;
  int maxNalloc;
  _MCMC_Link *links;
  double (*chi2funk)(double *, int, int, double *, double *, double *, void *);
  int NJD;
  double *t;
  double *mag;
  double *err;
  void *userparam;

  double *rvec;
  double *sample_covar_sum1;
  double **sample_covar_sum2;
  double **covar;

  int doprintchain;
  char *outchainfilename;
  char *outchainheader;
  int printevery;
  int Nauxil;
  double **auxil_params;
} _MCMC_Chain;
#endif
