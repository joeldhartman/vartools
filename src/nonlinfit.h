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
#include "mcmcfit.h"

typedef struct {
  int lcid, threadid;
  ProgramData *p;
  _Nonlinfit *c;
  int mcmcevalnum;
  FILE *mcmcoutfile;
} _NonLinFit_Chi2Struct;

#define VARTOOLS_NONLINFIT_DEFAULT_AMOEBA_TOLERANCE 1.0e-10
#define VARTOOLS_NONLINFIT_DEFAULT_AMOEBA_MAXSTEPS 5000
#define VARTOOLS_NONLINFIT_DEFAULT_MCMC_NLINKSTOTAL 100000
#define VARTOOLS_NONLINFIT_DEFAULT_MCMC_BURNIN 0.1
#define VARTOOLS_NONLINFIT_DEFAULT_MCMC_EPS 0.1
#define VARTOOLS_NONLINFIT_DEFAULT_MAX_MEM_STORE 4.0

#define VARTOOLS_COVARIANCE_KERNEL_SQUAREDEXPONENTIAL 0
#define VARTOOLS_COVARIANCE_KERNEL_EXPONENTIAL 1
#define VARTOOLS_COVARIANCE_KERNEL_MATERN 2

#define VARTOOLS_LOG_COVARIANCE_TINY -23
