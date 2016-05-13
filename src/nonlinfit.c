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
#include "nonlinfit.h"

double EvaluateChi2NonLinFit_Final(double *param, int Nvar, int NJD, double *t, double *mag, double *sig, void *fitstruct);

double EvaluateChi2NonLinFit(double *param, int Nvar, int NJD, double *t, double *mag, double *sig, void *fitstruct);

double Nonlinfit_GetAmoeba_Uncertainties(double *param, int paramindx, int Nparam, double chi20, double errval, int NJD, double *t, double *mag, double *err, void *userparam);

void Nonlinfit_Initialize_DataCovarMat(int NJD, double *t, double ***Cov_base, double ***Cov_calc, int **Nvec, int *sizemat, int *store_NJD, int CovKernelType) {
  int s, i, j;
  double **cov_tmp1;
  double **cov_tmp2;
  int *nvec_tmp;
  s = *sizemat;
  *store_NJD = NJD;
  if(NJD > s) {
    if(!s) {
      if((cov_tmp1 = (double **) malloc(NJD * sizeof(double *))) == NULL ||
	 (cov_tmp2 = (double **) malloc(NJD * sizeof(double *))) == NULL ||
	 (nvec_tmp = (int *) malloc(NJD * sizeof(int))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0; i < NJD; i++) {
	if((cov_tmp1[i] = (double *) malloc(NJD * sizeof(double))) == NULL ||
	   (cov_tmp2[i] = (double *) malloc(NJD * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
      }
      *sizemat = NJD;
      *Cov_base = cov_tmp1;
      *Cov_calc = cov_tmp2;
      *Nvec = nvec_tmp;
    } else {
      cov_tmp1 = *Cov_base;
      cov_tmp2 = *Cov_calc;
      nvec_tmp = *Nvec;
      if((cov_tmp1 = (double **) realloc(cov_tmp1, NJD * sizeof(double *))) == NULL ||
	 (cov_tmp2 = (double **) realloc(cov_tmp2, NJD * sizeof(double *))) == NULL ||
	 (nvec_tmp = (int *) realloc(nvec_tmp, NJD * sizeof(int))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0; i < s; i++) {
	if((cov_tmp1[i] = (double *) realloc(cov_tmp1[i], NJD * sizeof(double))) == NULL ||
	   (cov_tmp2[i] = (double *) realloc(cov_tmp2[i], NJD * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
      }
      for(i=s; i < NJD; i++) {
	if((cov_tmp1[i] = (double *) malloc(NJD * sizeof(double))) == NULL ||
	   (cov_tmp2[i] = (double *) malloc(NJD * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
      }
    }
    *sizemat = NJD;
    *Cov_base = cov_tmp1;
    *Cov_calc = cov_tmp2;
    *Nvec = nvec_tmp;
  } else {
    cov_tmp1 = *Cov_base;
    cov_tmp2 = *Cov_calc;
  }
  switch(CovKernelType){
  case VARTOOLS_COVARIANCE_KERNEL_SQUAREDEXPONENTIAL:
    for(i=0;i<NJD;i++) {
      for(j=i; j < NJD; j++) {
	cov_tmp1[i][j] = -(t[i]-t[j])*(t[i]-t[j]);
	cov_tmp2[i][j] = 0.0;
      }
    }
    break;
  case VARTOOLS_COVARIANCE_KERNEL_EXPONENTIAL:
    for(i=0;i<NJD;i++) {
      for(j=i; j < NJD; j++) {
	cov_tmp1[i][j] = -fabs(t[i]-t[j]);
	cov_tmp2[i][j] = 0.0;
      }
    }
    break;
  case VARTOOLS_COVARIANCE_KERNEL_MATERN:
    for(i=0;i<NJD;i++) {
      for(j=i; j < NJD; j++) {
	cov_tmp1[i][j] = fabs(t[i]-t[j]);
	cov_tmp2[i][j] = 0.0;
      }
    }
    break;
  default:
    error(ERR_CODEERROR);
  }
}

double Nonlinfit_Use_DataCovarMat(int NJD, double *t, double *obs, double *mod, double *err, double Cov_amp, double Cov_rho, double Cov_nu, double **Cov_base, double **Cov_calc, int *Nvec, int CovKernelType, double *lndet, int Output_Pred, double *GP_pred_mean, double *GP_pred_stddev) {
  double *x, *p, *b, retval;
  double ri, rk, rip, rkp, x1, x2, constterm;
  int i, j;
  double *K, *CovK, **Cov_calc_store;
  if((x = (double *) malloc(NJD*sizeof(double))) == NULL ||
     (p = (double *) malloc(NJD*sizeof(double))) == NULL ||
     (b = (double *) malloc(NJD*sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  if(Output_Pred) {
    if((CovK = (double *) malloc(NJD * sizeof(double))) == NULL ||
       (K = (double *) malloc(NJD * sizeof(double))) == NULL ||
       (Cov_calc_store = (double **) malloc(NJD * sizeof(double *))) == NULL)
      error(ERR_MEMALLOC);
    for(i=0; i < NJD; i++) {
      if((Cov_calc_store[i] = (double *) malloc(NJD * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
  }

  if(Cov_amp < 0.0) Cov_amp = -Cov_amp;
  if(Cov_rho < 0.0) Cov_rho = -Cov_rho;
  if(CovKernelType == VARTOOLS_COVARIANCE_KERNEL_MATERN) {
    if(Cov_nu < 0.0) Cov_nu = -Cov_nu;
  }

  switch(CovKernelType){
  case VARTOOLS_COVARIANCE_KERNEL_SQUAREDEXPONENTIAL:
    Cov_rho = Cov_rho*Cov_rho;
    for(i=0; i < NJD; i++) {
      b[i] = (obs[i]-mod[i]);
      Cov_calc[i][i] = err[i]*err[i] + Cov_amp;
      if(Cov_rho > 0.0) {
	for(j=i+1; j<NJD; j++) {
	  if(Cov_base[i][j]/2.0/Cov_rho < VARTOOLS_LOG_COVARIANCE_TINY) {
	    Cov_calc[i][j] = 0.0;
	    /* If this term is too small, all additional terms will be too
               small, so quit now */
	    break;
	  }
	  Cov_calc[i][j] = Cov_amp*exp(Cov_base[i][j]/2.0/Cov_rho);
	}
      }
    }
    break;
  case VARTOOLS_COVARIANCE_KERNEL_EXPONENTIAL:
    for(i=0; i < NJD; i++) {
      b[i] = (obs[i]-mod[i]);
      Cov_calc[i][i] = err[i]*err[i] + Cov_amp;
      if(Cov_rho > 0.0) {
	for(j=i+1; j<NJD; j++) {
	  if(Cov_base[i][j]/Cov_rho < VARTOOLS_LOG_COVARIANCE_TINY) {
	    Cov_calc[i][j] = 0.0;
	    /* If this term is too small, all additional terms will be too
               small, so quit now */
	    break;
	  }
	  Cov_calc[i][j] = Cov_amp*exp(Cov_base[i][j]/Cov_rho);
	}
      }
    }
    break;
  case VARTOOLS_COVARIANCE_KERNEL_MATERN:
    if(Cov_rho > 0.0) {
      x1 = sqrt(2.0*Cov_nu)/Cov_rho;
      constterm = Cov_amp/(exp(gammln(Cov_nu))*pow(2.0,(Cov_nu-1.0)));
    }
    for(i=0; i < NJD; i++) {
      b[i] = (obs[i]-mod[i]);
      Cov_calc[i][i] = err[i]*err[i] + Cov_amp;
      if(Cov_rho > 0.0) {
	for(j=i+1; j<NJD; j++) {
	  x2 = x1*Cov_base[i][j];
	  if(-x2 < VARTOOLS_LOG_COVARIANCE_TINY) {
	    /* If this term is too small, all additional terms will be too
               small, so quit now */
	    Cov_calc[i][j] = 0.0;
	    break;
	  }
	  bessik(x2,Cov_nu,&ri,&rk,&rip,&rkp);
	  Cov_calc[i][j] = Cov_amp*constterm*(pow(x2,Cov_nu))*rk;
	}
      }
    }
    break;
  default:
    error(ERR_CODEERROR);
    break;
  }

  /* Make a back-up copy of the covariance matrix if necessary */
  if(Output_Pred) {
    for(i=0; i < NJD; i++) {
      for(j=0; j < NJD; j++) {
	Cov_calc_store[i][j] = Cov_calc[i][j];
      }
      Cov_calc_store[i][i] -= err[i]*err[i];
    }
  }

  /* We need to compute x^T*(Cov^{-1})*x and det(Cov), first find the
     Cholesky decomposition of Cov taking advantage of the fact that
     the matrix is sparse and near diagonal, then do back-substitution
     to get (Cov^{-1})*x. Finally dot with x^T and compute the
     log determinant of Cov from the sum of the log of the diagonal
     elements of the decomposed matrix */
  choldc_sparse_neardiag(Cov_calc, NJD, Nvec, p);
  cholsl_sparse_neardiag(Cov_calc, NJD, Nvec, p, b, x);
  *lndet = log(p[0]);
  retval = x[0]*b[0];
  for(i=1; i < NJD; i += 1) {
    *lndet = (*lndet) + log(p[i]);
    retval += x[i]*b[i];
  }

  if(Output_Pred) {
    /* Compute the predicted mean and standard deviation of the Gaussian
       Process if requested */
    for(i=0; i < NJD; i++) {
      GP_pred_mean[i] = 0.;
      for(j=i; j < NJD; j++) {
	if(Cov_calc_store[i][j] > 0.) {
	  GP_pred_mean[i] += Cov_calc_store[i][j]*x[j];
	  K[j] = Cov_calc_store[i][j];
	} else {
	  K[j] = 0.;
	}
      }
      for(j=i-1; j >= 0; j--) {
	if(Cov_calc_store[j][i] > 0.) {
	  GP_pred_mean[i] += Cov_calc_store[j][i]*x[j];
	  K[j] = Cov_calc_store[j][i];
	} else {
	  K[j] = 0.;
	}
      }
      cholsl_sparse_neardiag(Cov_calc, NJD, Nvec, p, K, CovK);
      GP_pred_stddev[i] = K[i] + err[i]*err[i];
      for(j=i; j < NJD; j++) {
	if(K[j] > 0.) {
	  GP_pred_stddev[i] -= K[j]*CovK[j];
	} else {
	  break;
	}
      }
      for(j=i-1; j >= 0; j--) {
	if(K[j] > 0.) {
	  GP_pred_stddev[i] -= K[j]*CovK[j];
	} else {
	  break;
	}
      }
      if(GP_pred_stddev[i] > 0.)
	GP_pred_stddev[i] = sqrt(GP_pred_stddev[i]);
      else
	GP_pred_stddev[i] = 0.;
    }
    free(K);
    free(CovK);
    for(i=0; i < NJD; i++) free(Cov_calc_store[i]);
    free(Cov_calc_store);
  }
      

  free(x); free(p); free(b);
  return retval;
}

int Nonlinfit_RunAmoebaFit(ProgramData *p, _Nonlinfit *c, int threadid, int lcid) {
  /* This will carry-out the amoeba downhill-simplex fit for the nonlinfit
     command */

  _NonLinFit_Chi2Struct s;
  double **simplex = NULL;
  double *chi2vals = NULL;
  int Nparam = 0;
  int Ntovary = 0;
  int *ia = NULL;
  double initval, errval, chi2best;
  int i, ibest;
  int NJD;
  double *t, *mag, *err;

  int amoeba_val;
  int nfunkeval;

  s.lcid = lcid;
  s.threadid = threadid;
  s.p = p;
  s.c = c;
  
  NJD = p->NJD[threadid];
  t = p->t[threadid];
  mag = p->mag[threadid];
  err = p->sig[threadid];

  /* Initialize the correlation parameters */
  if(c->use_covar) {
    if(c->Corr_amp_expr != NULL) {
      initval = EvaluateExpression(lcid, threadid, 0, c->Corr_amp_expr);
      if(initval <= 0.0) {
	error(ERR_NEGATIVE_COVAR_PARAM);
      }
      SetVariable_Value_Double(lcid, threadid, 0, c->Corr_amp_var, initval);
    }
    if(c->Corr_rho_expr != NULL) {
      initval = EvaluateExpression(lcid, threadid, 0, c->Corr_rho_expr);
      if(initval <= 0.0) {
	error(ERR_NEGATIVE_COVAR_PARAM);
      }
      SetVariable_Value_Double(lcid, threadid, 0, c->Corr_rho_var, initval);
    }
    if(c->covar_type == VARTOOLS_COVARIANCE_KERNEL_MATERN) {
      if(c->Corr_nu_expr != NULL) {
	initval = EvaluateExpression(lcid, threadid, 0, c->Corr_nu_expr);
	if(initval <= 0.0) {
	  error(ERR_NEGATIVE_COVAR_PARAM);
	}
	SetVariable_Value_Double(lcid, threadid, 0, c->Corr_nu_var, initval);
      }
    }
  }

  /* Initialize the simplex using incrementparameters_foramoeba */
  for(i=0; i < c->Nparams; i++) {
    /* Evaluate the initialization expressions for the parameters */
    initval = EvaluateExpression(lcid, threadid, 0, c->paraminit_expressions[i]);
    errval = EvaluateExpression(lcid, threadid, 0, c->paramerr_expressions[i]);
    incrementparameters_foramoeba(&Nparam, &Ntovary, &simplex, &ia, 1, initval, errval);
  }
  /* Set the initial chi2 values */
  amoeba_initializesimplexchi2(Nparam, Ntovary, simplex, &chi2vals, &EvaluateChi2NonLinFit, NJD, t, mag, err, (void *) (&s));

  /* Run the fit */
  amoeba_val = amoeba(simplex, chi2vals, ia, Nparam, c->amoeba_tol,
		      &EvaluateChi2NonLinFit, &nfunkeval, c->amoeba_maxsteps,
		      NJD, t, mag, err, (void *) (&s));

  if(c->fittype == VARTOOLS_NONLINFIT_FITTYPE_AMOEBA) {
    if(amoeba_val) {
      c->amoeba_isconverged[threadid] = 0;
    } else {
      c->amoeba_isconverged[threadid] = 1;
    }
  }
  
  /* Find the optimum parameters and store the results */
  chi2best = chi2vals[0];
  ibest = 0;
  for(i=1; i < Ntovary+1; i++) {
    if(!isnan(chi2vals[i]) && (chi2vals[i] < chi2best || isnan(chi2best))) {
      chi2best = chi2vals[i];
      ibest = i;
    }
  }
  
  /* Run the model one last time with the optimum parameters to update the
     model variables etc. */
  c->chi2out[threadid] = EvaluateChi2NonLinFit_Final(simplex[ibest], Nparam, NJD, t, mag, err, (void *) (&s));
  
  for(i=0; i < c->Nparams; i++) {
    c->param_outvals[threadid][i] = EvaluateVariable_Double(lcid, threadid, 0, c->params[i]);
  }

  
  if(c->fittype == VARTOOLS_NONLINFIT_FITTYPE_AMOEBA) {
    if(!amoeba_val) {
      for(i=0; i < c->Nparams; i++) {
	errval = EvaluateExpression(lcid, threadid, 0, 
				    c->paramerr_expressions[i]);
	c->param_uncertainties[threadid][i] = 
	  Nonlinfit_GetAmoeba_Uncertainties(c->param_outvals[threadid], i, 
					    Nparam,
					    c->chi2out[threadid], errval,
					    NJD, t, mag,
					    err, (void *) (&s));
      }
    } else {
      for(i=0; i < c->Nparams; i++) {
	c->param_uncertainties[threadid][i] = 0.;
      }
    }
  }
	  

  amoeba_cleanup(&Nparam, &Ntovary, &simplex, &ia, &chi2vals);
   
}

double Nonlinfit_GetAmoeba_Uncertainties(double *param, int paramindx, int Nparam, double chi20, double errval, int NJD, double *t, double *mag, double *err, void *userparam) {
  int i, j, Ntest;
  double chi2test;
  double delp;
  double *ptrial;
  if((ptrial = (double *) malloc(Nparam * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0; i < Nparam; i++)
    ptrial[i] = param[i];
  if(!errval) {
    if(param[paramindx] != 0) {
      delp = 0.001*param[paramindx];
    } else
      delp = 0.001;
  }
  else
    delp = errval;
  ptrial[paramindx] = param[paramindx] + delp;
  Ntest = 0;
  chi2test = EvaluateChi2NonLinFit(ptrial, Nparam, NJD, t, mag, err, userparam);
  if(chi2test == chi20) {
    free(ptrial);
    return 0.0;
  }
  while(Ntest < 100 ? (isnan(chi2test) ? 1 : (fabs(chi2test - chi20) < 0.99 ||
					      fabs(chi2test - chi20) > 1.01))
	: 0) {
    while(Ntest < 100 && isnan(chi2test)) {
      delp = -delp;
      ptrial[paramindx] = param[paramindx] + delp;
      chi2test = EvaluateChi2NonLinFit(ptrial, Nparam, NJD, t, mag, err, userparam);
      Ntest++;
      if(!isnan(chi2test)) break;
      delp = delp/10.;
      ptrial[paramindx] = param[paramindx] + delp;
      chi2test = EvaluateChi2NonLinFit(ptrial, Nparam, NJD, t, mag, err, userparam);
      Ntest++;
    }
    if(chi2test == chi20) {
      free(ptrial);
      return 0.0;
    }
    delp = delp/(sqrt(fabs(chi2test - chi20)));
    ptrial[paramindx] = param[paramindx] + delp;
    chi2test = EvaluateChi2NonLinFit(ptrial, Nparam, NJD, t, mag, err, userparam);
    Ntest++;
  }
  free(ptrial);
  return delp;
}

void Nonlinfit_GetMCMC_ChainStats(ProgramData *p, _Nonlinfit *c, int threadid, int lcid, _MCMC_Chain *chain) 
{
  int maxNchain, istart, iend, Nlinks_stat, i, j, k, m, z, Npct;
  double *statsvec;
  
  /* Calculate the MCMC statistics and output the chain if requested */
  maxNchain = (chain->Nlinks_run < chain->Nlinks_allocated ? 
	       chain->Nlinks_run : chain->Nlinks_allocated);
  istart = (chain->Nlinks_run < chain->Nlinks_allocated ?
	    floor(c->mcmc_burninfrac*maxNchain) :
	    (chain->activeindex + floor(c->mcmc_burninfrac*maxNchain) 
	     < chain->Nlinks_allocated ?
	     chain->activeindex + floor(c->mcmc_burninfrac*maxNchain) :
	     chain->activeindex + floor(c->mcmc_burninfrac*maxNchain) -
	     chain->Nlinks_allocated));
  iend = istart + floor((1.0 - c->mcmc_burninfrac)*maxNchain);
  iend = iend < chain->Nlinks_allocated ? iend :
    iend - chain->Nlinks_allocated;

  Nlinks_stat = (istart <= iend ? iend - istart + 1 :
		 iend + 1 + (chain->Nlinks_allocated - istart));
  if((statsvec = (double *) malloc(Nlinks_stat * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0, z=0; i < c->N_mcmc_chain_expressions; i++) {
    m=0;
    if(istart <= iend) {
      for(j=istart; j <= iend; j++, m++) {
	for(k=0; k < chain->Nparam; k++) {
	  SetVariable_Value_Double(lcid, threadid, 0, c->params[k], chain->links[j].p[k]);
	}
	statsvec[m] = EvaluateExpression(lcid, threadid, 0, c->mcmc_chain_stats_expressions[i]);
      }
    } else {
      for(j=0; j <= iend; j++, m++) {
	for(k=0; k < chain->Nparam; k++) {
	  SetVariable_Value_Double(lcid, threadid, 0, c->params[k], chain->links[j].p[k]);
	}
	statsvec[m] = EvaluateExpression(lcid, threadid, 0, c->mcmc_chain_stats_expressions[i]);
      }
      for(j=istart; j < chain->Nlinks_allocated; j++, m++) {
	for(k=0; k < chain->Nparam; k++) {
	  SetVariable_Value_Double(lcid, threadid, 0, c->params[k], chain->links[j].p[k]);
	}
	statsvec[m] = EvaluateExpression(lcid, threadid, 0, c->mcmc_chain_stats_expressions[i]);
      }
    }
    for(k=0, Npct=0; k < c->N_mcmc_chain_stats; k++, z++) {
      switch(c->mcmc_statstocalc[k]) {
      case VARTOOLS_STATSTYPE_MEAN:
	c->mcmc_statsout[threadid][z] = getmean(m,statsvec);
	break;
      case VARTOOLS_STATSTYPE_WEIGHTEDMEAN:
	c->mcmc_statsout[threadid][z] = getmean(m,statsvec);
	break;
      case VARTOOLS_STATSTYPE_MEDIAN:
	c->mcmc_statsout[threadid][z] = median(m,statsvec);
	break;
      case VARTOOLS_STATSTYPE_MEDIAN_WEIGHT:
	c->mcmc_statsout[threadid][z] = median(m,statsvec);
	break;
      case VARTOOLS_STATSTYPE_STDDEV:
	c->mcmc_statsout[threadid][z] = stddev(m,statsvec);
	break;
      case VARTOOLS_STATSTYPE_MEDDEV:
	c->mcmc_statsout[threadid][z] = meddev(m,statsvec);
	break;
      case VARTOOLS_STATSTYPE_MEDMEDDEV:
	c->mcmc_statsout[threadid][z] = medmeddev(m,statsvec);
	break;
      case VARTOOLS_STATSTYPE_MAD:
	c->mcmc_statsout[threadid][z] = MAD(m,statsvec);
	break;
      case VARTOOLS_STATSTYPE_KURTOSIS:
	c->mcmc_statsout[threadid][z] = kurtosis(m,statsvec);
	break;
      case VARTOOLS_STATSTYPE_SKEWNESS:
	c->mcmc_statsout[threadid][z] = skewness(m,statsvec);
	break;
      case VARTOOLS_STATSTYPE_PERCENTILE:
	c->mcmc_statsout[threadid][z] = percentile(m,statsvec,
						   c->pctval[Npct]);
	Npct++;
	break;
      case VARTOOLS_STATSTYPE_PERCENTILE_WEIGHT:
	c->mcmc_statsout[threadid][z] = percentile(m,statsvec,
						   c->pctval[Npct]);
	Npct++;
	break;
      case VARTOOLS_STATSTYPE_MAXIMUM:
	c->mcmc_statsout[threadid][z] = getmaximum(m,statsvec);
	break;
      case VARTOOLS_STATSTYPE_MINIMUM:
	c->mcmc_statsout[threadid][z] = getminimum(m,statsvec);
	break;
      case VARTOOLS_STATSTYPE_SUM:
	c->mcmc_statsout[threadid][z] = getsum(m,statsvec);
	break;
      default:
	error(ERR_CODEERROR);
      }
    }
  }
  free(statsvec);
}

int Nonlinfit_RunMCMCFit(ProgramData *p, _Nonlinfit *c, int threadid, int lcid) {
  /* This will carry-out the Differential Evolution MCMC fit for the nonlinfit
     command */

  _NonLinFit_Chi2Struct s;
  _MCMC_Chain *chain = NULL;

  double chi2best, initval, errval;

  int i, ibest, istart, iend;
  int NJD, maxNchain, Nlinks_stat;
  double *t, *mag, *err, maxmemperchain;
  double *statsvec;
  char outfilename[MAXLEN];
  char *outheader;
  int Nauxil;
  int omodel_store, correctlc_store;
  double **auxil_params = NULL;
  OutText str;
  str.s = NULL;
  str.space = 0;
  str.len_s = 0;
  str.Nchar_cur_line = 0;

  s.lcid = lcid;
  s.threadid = threadid;
  s.p = p;
  s.c = c;
  s.mcmcevalnum = 0;
  
  if(c->mcmc_outchains) {
    GetOutputFilename(outfilename,p->lcnames[lcid],c->mcmc_outchains_dir,
		      "mcmc",
		      c->mcmc_outchains_format, lcid);
  } else {
    outfilename[0] = '\0';
  }

  NJD = p->NJD[threadid];
  t = p->t[threadid];
  mag = p->mag[threadid];
  err = p->sig[threadid];

  /* Initialize the correlation parameters */
  if(c->use_covar) {
    if(c->Corr_amp_expr != NULL) {
      initval = EvaluateExpression(lcid, threadid, 0, c->Corr_amp_expr);
      if(initval <= 0.0) {
	error(ERR_NEGATIVE_COVAR_PARAM);
      }
      SetVariable_Value_Double(lcid, threadid, 0, c->Corr_amp_var, initval);
    }
    if(c->Corr_rho_expr != NULL) {
      initval = EvaluateExpression(lcid, threadid, 0, c->Corr_rho_expr);
      if(initval <= 0.0) {
	error(ERR_NEGATIVE_COVAR_PARAM);
      }
      SetVariable_Value_Double(lcid, threadid, 0, c->Corr_rho_var, initval);
    }
    if(c->covar_type == VARTOOLS_COVARIANCE_KERNEL_MATERN) {
      if(c->Corr_nu_expr != NULL) {
	initval = EvaluateExpression(lcid, threadid, 0, c->Corr_nu_expr);
	if(initval <= 0.0) {
	  error(ERR_NEGATIVE_COVAR_PARAM);
	}
	SetVariable_Value_Double(lcid, threadid, 0, c->Corr_nu_var, initval);
      }
    }
  }

  /* Run downhill simplex first unless it is turned off */
  if(!c->mcmc_skipamoeba) {
    omodel_store = c->omodel;
    correctlc_store = c->correctlc;
    c->omodel = 0;
    c->correctlc = 0;
    Nonlinfit_RunAmoebaFit(p, c, threadid, lcid);
    c->omodel = omodel_store;
    c->correctlc = correctlc_store;
  }

  if(p->readallflag) {
    maxmemperchain = (c->mcmc_max_mem_store / ((double) p->Nlcs));
  } else {
#ifdef PARALLEL
    if(p->Nproc_allow > 1) {
      maxmemperchain = (c->mcmc_max_mem_store / ((double) p->Nproc_allow));
    } else {
#endif
      maxmemperchain = c->mcmc_max_mem_store;
#ifdef PARALLEL
    }
#endif
  }

  /* Setup the parameters in the chain using MCMC_increment_parameter */
  for(i=0; i < c->Nparams; i++) {
    /* Evaluate the initialization expressions for the parameters */
    if(c->mcmc_skipamoeba)
      initval = EvaluateExpression(lcid, threadid, 0, c->paraminit_expressions[i]);
    else
      initval = c->param_outvals[threadid][i];
    errval = EvaluateExpression(lcid, threadid, 0, c->paramerr_expressions[i]);

    MCMC_increment_parameter(&chain, initval, errval);
  }

  if(c->mcmc_outchains) {
    printtostring_nowrap(&str,
			 "#");
    for(i=0; i < c->Nparams; i++) {
      printtostring_nowrap(&str," ");
      printtostring_nowrap(&str,c->paramnames[i]);
    }
    if(c->uselinfit) {
      for(i=0; i < c->linfit->Nparams; i++) {
	printtostring_nowrap(&str," ");
	printtostring_nowrap(&str,c->linfit->paramnames[i]);
      }
      Nauxil = c->linfit->Nparams;
      auxil_params = &(c->linfit->param_outvals[threadid]);
    }
    else {
      Nauxil = 0;
      auxil_params = NULL;
    }
    printtostring_nowrap(&str," -2ln(L)");
    outheader = str.s;
  }
  else {
    outheader = NULL;
    Nauxil = 0;
    auxil_params = NULL;
  }

  /* Initialize the chain now that all parameters have been set */
  MCMC_initialize_chain(chain, c->mcmc_eps, maxmemperchain, c->mcmc_outchains,
			outfilename, outheader, c->mcmc_outchains_print_every,
			Nauxil, auxil_params, &EvaluateChi2NonLinFit, 
			NJD, t, mag, err, (void *) (&s));

  /* Run the MCMC procedure */
  MCMC_DifferentialEvolution_RunMCMC(c->mcmc_Naccept, c->mcmc_Nlinkstotal, 
				     MCMC_DEFAULT_RECALC_COVAR, chain);


  /* Calculate the statistics from the chain to report on the output
     table */
  Nonlinfit_GetMCMC_ChainStats(p, c, threadid, lcid, chain);

  /* Find the optimum parameters */
  ibest = MCMC_GetBestLink(chain);

  /* Run the model one last time with the optimum parameters to update the
     model variables etc. */
  c->chi2out[threadid] = EvaluateChi2NonLinFit_Final(chain->links[ibest].p, chain->Nparam, NJD, t, mag, err, (void *) (&s));
  
  for(i=0; i < c->Nparams; i++) {
    c->param_outvals[threadid][i] = EvaluateVariable_Double(lcid, threadid, 0, c->params[i]);
  }

  if(str.s != NULL) free(str.s);

  /* Clean up */
  MCMC_CleanUp(&chain);
   
}


void DoNonlinfit(ProgramData *p, _Nonlinfit *c, int threadid, int lcid) {
  /* This function executes the non-linear fit command. */

  if(c->use_covar) {
    Nonlinfit_Initialize_DataCovarMat(p->NJD[threadid], p->t[threadid], &(c->Corr_mat1[threadid]), &(c->Corr_mat2[threadid]), &(c->Corr_Nvec[threadid]), &(c->Corr_sizemat[threadid]), &(c->Corr_store_NJD[threadid]), c->covar_type);
  }
  if(c->fittype == VARTOOLS_NONLINFIT_FITTYPE_AMOEBA)
    Nonlinfit_RunAmoebaFit(p, c, threadid, lcid);
  else if(c->fittype == VARTOOLS_NONLINFIT_FITTYPE_DEMCMC)
    Nonlinfit_RunMCMCFit(p, c, threadid, lcid);
  return;
}

/* Return the chi2 value for the non-linear fit function, and also
   write out the model if asked to */
double EvaluateChi2NonLinFit_Final(double *param, int Nvar, int NJD, double *t, double *mag, double *sig, void *fitstruct)
{
  _NonLinFit_Chi2Struct *f;
  int i, j, Ngood;
  double model, obs, err, chi2out, term, *model_vec, *err_vec, *obs_vec, *t_vec;
  double *GP_pred_mean, *GP_pred_stddev;
  double Corr_rho_val, Corr_amp_val, Corr_nu_val, lndet;
  FILE *outfile;
  char lcoutname[MAXLEN];
  _Linfit linfitcpy;
  _Linfit *linfittouse;
  f = (_NonLinFit_Chi2Struct *) fitstruct;

  /* update the fit variables to have the trial values */
  for(i=0; i < Nvar; i++) {
    SetVariable_Value_Double(f->lcid, f->threadid, 0, f->c->params[i], param[i]);
  }

  /* Check the parameter constraints, return NaN if one of them is not met */
  if(f->c->Nconstraints) {
    for(i=0; i < f->c->Nconstraints; i++) {
      if(!EvaluateExpression(f->lcid, f->threadid, 0, f->c->constraint_expressions[i])) {
	term = sqrt(-1.0);
	return term;
      }
    }
  }

  /* Do linear optimization if requested, here we will make a dummy version
     if the Linfit structure with omodel and/or correctlc set to 1 if used by
     the nonlinfit command, so that these appropriate actions will be
     executed by the DoLinfit function. */
  if(f->c->uselinfit) {
    if(f->c->omodel || f->c->correctlc) {
      memcpy(&linfitcpy, f->c->linfit, sizeof(_Linfit));
      linfittouse = &linfitcpy;
      linfitcpy.omodel = f->c->omodel;
      linfitcpy.correctlc = f->c->correctlc;
      if(f->c->omodel) {
	free(linfitcpy.outfile_extension);
	free(linfitcpy.outfilename_format);
	linfitcpy.outdir = f->c->outdir;
	linfitcpy.outfile_extension = f->c->outfile_extension;
	linfitcpy.outfilename_format = f->c->outfilename_format;
      }
    }
    else {
      linfittouse = f->c->linfit;
    }
    DoLinfit(f->p, linfittouse, f->threadid, f->lcid);
    if(f->c->omodel || f->c->correctlc) {
      for(i=0; i < linfittouse->Nparams; i++) {
	f->c->linfit->param_outvals[f->threadid][i] = linfittouse->param_outvals[f->threadid][i];
	f->c->linfit->param_uncertainties[f->threadid][i] = linfittouse->param_uncertainties[f->threadid][i];
      }
      f->c->linfit->chi2out[f->threadid] = linfittouse->chi2out[f->threadid];
    }
    chi2out = linfittouse->chi2out[f->threadid];
  }
  else {
    /* Otherwise evaluate the model expression directly */
    chi2out = 0.;
    if(f->c->omodel) {
      GetOutputFilename(lcoutname,f->p->lcnames[f->lcid],f->c->outdir,
			f->c->outfile_extension,
			f->c->outfilename_format, f->lcid);
      if((outfile = fopen(lcoutname,"w")) == NULL) {
	fprintf(stderr,"Cannot write to %s\n", lcoutname);
	exit(ERR_CANNOTWRITE);
      }
      fprintf(outfile,"#Time Mag_lc Err Mag_model");
      if(f->c->use_covar) {
	fprintf(outfile," GP_Mean GP_Stddev");
      }
      fprintf(outfile,"\n");
    }

    if(!(f->c->use_covar)) {
      for(j=0; j < NJD; j++) {
	model = EvaluateExpression(f->lcid, f->threadid, j, f->c->functionexpression);
	if(f->c->modelvarname != NULL) {
	  SetVariable_Value_Double(f->lcid, f->threadid, j, f->c->modelvar, model);
	}
	obs = mag[j];
	if(f->c->errorstring == NULL)
	  err = sig[j];
	else
	  err = EvaluateExpression(f->lcid, f->threadid, j, f->c->errorexpression);
	if(!isnan(obs) && ! isnan(err) && err > 0.) {
	  if(f->c->omodel) {
	    fprintf(outfile,"%.17g %.17g %.17g %.17g\n", t[j], obs, err, model);
	  }
	  term = (model - obs)/err;
	  if(f->c->correctlc) {
	    mag[j] -= model;
	  }
	  chi2out += term*term;
	  if(f->c->errorstring != NULL) {
	    chi2out += (2.0*log(err));
	  }
	}
      }
    } else {
      if((model_vec = (double *) malloc(NJD * sizeof(double))) == NULL ||
	 (err_vec = (double *) malloc(NJD * sizeof(double))) == NULL ||
	 (obs_vec = (double *) malloc(NJD * sizeof(double))) == NULL ||
	 (t_vec = (double *) malloc(NJD * sizeof(double))) == NULL ||
	 (GP_pred_mean = (double *) malloc(NJD * sizeof(double))) == NULL ||
	 (GP_pred_stddev = (double *) malloc(NJD * sizeof(double))) == NULL) {
	error(ERR_MEMALLOC);
      }
      for(j=0, Ngood=0; j < NJD; j++) {
	model_vec[Ngood] = EvaluateExpression(f->lcid, f->threadid, j, f->c->functionexpression);
	if(f->c->modelvarname != NULL) {
	  SetVariable_Value_Double(f->lcid, f->threadid, j, f->c->modelvar, model_vec[Ngood]);
	}
	t_vec[Ngood] = t[j];
	obs_vec[Ngood] = mag[Ngood];
	if(f->c->errorstring == NULL)
	  err_vec[Ngood] = sig[j];
	else
	  err_vec[Ngood] = EvaluateExpression(f->lcid, f->threadid, j, f->c->errorexpression);
	if(!isnan(obs_vec[Ngood]) && ! isnan(err_vec[Ngood]) && err_vec[Ngood] > 0.) {
	  /*if(f->c->omodel) {
	    fprintf(outfile,"%.17g %.17g %.17g %.17g\n", t[j], obs_vec[Ngood], err_vec[Ngood], model_vec[Ngood]);
	    }*/
	  if(f->c->correctlc) {
	    mag[j] -= model_vec[Ngood];
	  }
	  Ngood++;
	}
      }
      Corr_amp_val = EvaluateVariable_Double(f->lcid, f->threadid, 0, f->c->Corr_amp_var);
      Corr_rho_val = EvaluateVariable_Double(f->lcid, f->threadid, 0, f->c->Corr_rho_var);
      if(f->c->covar_type == VARTOOLS_COVARIANCE_KERNEL_MATERN) {
	Corr_nu_val = EvaluateVariable_Double(f->lcid, f->threadid, 0, f->c->Corr_nu_var);
      } else {
	Corr_nu_val = 0.0;
      }
      if(f->c->Corr_store_NJD[f->threadid] != Ngood) {
	Nonlinfit_Initialize_DataCovarMat(Ngood, t_vec, &(f->c->Corr_mat1[f->threadid]), &(f->c->Corr_mat2[f->threadid]), &(f->c->Corr_Nvec[f->threadid]), &(f->c->Corr_sizemat[f->threadid]), &(f->c->Corr_store_NJD[f->threadid]), f->c->covar_type);
      }
      chi2out = Nonlinfit_Use_DataCovarMat(Ngood, t_vec, obs_vec, model_vec, err_vec, Corr_amp_val, Corr_rho_val, Corr_nu_val, f->c->Corr_mat1[f->threadid], f->c->Corr_mat2[f->threadid], f->c->Corr_Nvec[f->threadid], f->c->covar_type, &lndet, f->c->omodel, GP_pred_mean, GP_pred_stddev);

      fprintf(stderr,"#lndet =  %.17g\n", lndet);

      if(f->c->omodel) {
	for(j=0; j < Ngood; j++) {
	  fprintf(outfile,"%.17g %.17g %.17g %.17g %.17g %.17g\n", t_vec[j], obs_vec[j], err_vec[j], model_vec[j], GP_pred_mean[j], GP_pred_stddev[j]);
	}
      }
      chi2out += lndet;
      free(model_vec);
      free(err_vec);
      free(obs_vec);
      free(t_vec);
      free(GP_pred_mean);
      free(GP_pred_stddev);
    }
    if(f->c->omodel) {
      fclose(outfile);
    }
  }


  /* Evaluate the contribution of the priors to chi2 */
  for(i=0; i < f->c->Npriors; i++) {
    chi2out += EvaluateExpression(f->lcid, f->threadid, 0, f->c->prior_expressions[i]);
  }

  return chi2out;
}

double EvaluateChi2NonLinFit(double *param, int Nvar, int NJD, double *t, double *mag, double *sig, void *fitstruct)
{
  _NonLinFit_Chi2Struct *f;
  int i, j, Ngood;
  double model, obs, err, chi2out, term, *model_vec, *err_vec, *obs_vec, *t_vec;
  double Corr_rho_val, Corr_amp_val, Corr_nu_val, lndet;
  f = (_NonLinFit_Chi2Struct *) fitstruct;

  /* update the fit variables to have the trial values */
  for(i=0; i < Nvar; i++) {
    SetVariable_Value_Double(f->lcid, f->threadid, 0, f->c->params[i], param[i]);
  }

  /* Check the parameter constraints, return NaN if one of them is not met */
  if(f->c->Nconstraints) {
    for(i=0; i < f->c->Nconstraints; i++) {
      if(!EvaluateExpression(f->lcid, f->threadid, 0, f->c->constraint_expressions[i])) {
	term = sqrt(-1.0);
	return term;
      }
    }
  }

  /* Do linear optimization if requested */
  if(f->c->uselinfit) {
    DoLinfit(f->p, f->c->linfit, f->threadid, f->lcid);
    chi2out = f->c->linfit->chi2out[f->threadid];
  }
  else {
    /* Otherwise evaluate the model expression directly */
    if(!f->c->use_covar) {
      chi2out = 0.;
      for(j=0; j < NJD; j++) {
	model = EvaluateExpression(f->lcid, f->threadid, j, f->c->functionexpression);
	obs = mag[j];
	if(f->c->errorstring == NULL)
	  err = sig[j];
	else
	  err = EvaluateExpression(f->lcid, f->threadid, j, f->c->errorexpression);
	if(!isnan(obs) && ! isnan(err) && err > 0.) {
	  term = (model - obs)/err;
	  chi2out += term*term;
	  if(f->c->errorstring != NULL) {
	    chi2out += (2.0*log(err));
	  }
	}
      }
    } 
    else {
      if((model_vec = (double *) malloc(NJD * sizeof(double))) == NULL ||
	 (err_vec = (double *) malloc(NJD * sizeof(double))) == NULL ||
	 (obs_vec = (double *) malloc(NJD * sizeof(double))) == NULL ||
	 (t_vec = (double *) malloc(NJD * sizeof(double))) == NULL) {
	error(ERR_MEMALLOC);
      }
      for(j=0, Ngood=0; j < NJD; j++) {
	model_vec[Ngood] = EvaluateExpression(f->lcid, f->threadid, j, f->c->functionexpression);
	t_vec[Ngood] = t[j];
	obs_vec[Ngood] = mag[Ngood];
	if(f->c->errorstring == NULL)
	  err_vec[Ngood] = sig[j];
	else
	  err_vec[Ngood] = EvaluateExpression(f->lcid, f->threadid, j, f->c->errorexpression);
	if(!isnan(obs_vec[Ngood]) && ! isnan(err_vec[Ngood]) && err_vec[Ngood] > 0.) {
	  Ngood++;
	}
      }
      Corr_amp_val = EvaluateVariable_Double(f->lcid, f->threadid, 0, f->c->Corr_amp_var);
      Corr_rho_val = EvaluateVariable_Double(f->lcid, f->threadid, 0, f->c->Corr_rho_var);
      if(f->c->covar_type == VARTOOLS_COVARIANCE_KERNEL_MATERN) {
	Corr_nu_val = EvaluateVariable_Double(f->lcid, f->threadid, 0, f->c->Corr_nu_var);
      } else {
	Corr_nu_val = 0.0;
      }
      if(f->c->Corr_store_NJD[f->threadid] != Ngood) {
	Nonlinfit_Initialize_DataCovarMat(Ngood, t_vec, &(f->c->Corr_mat1[f->threadid]), &(f->c->Corr_mat2[f->threadid]), &(f->c->Corr_Nvec[f->threadid]), &(f->c->Corr_sizemat[f->threadid]), &(f->c->Corr_store_NJD[f->threadid]), f->c->covar_type);
      }
      chi2out = Nonlinfit_Use_DataCovarMat(Ngood, t_vec, obs_vec, model_vec, err_vec, Corr_amp_val, Corr_rho_val, Corr_nu_val, f->c->Corr_mat1[f->threadid], f->c->Corr_mat2[f->threadid], f->c->Corr_Nvec[f->threadid], f->c->covar_type, &lndet, 0., NULL, NULL);
      chi2out += lndet;
      free(model_vec);
      free(err_vec);
      free(obs_vec);
      free(t_vec);
    }
  }

  /* Evaluate the contribution of the priors to chi2 */
  for(i=0; i < f->c->Npriors; i++) {
    chi2out += EvaluateExpression(f->lcid, f->threadid, 0, f->c->prior_expressions[i]);
  }

  return chi2out;
}

void InitNonlinfit(ProgramData *p, _Nonlinfit *c, int cnum) {
  /* Initialize the Nonlinfit command */
  
  int nparam, s;
  char **paramnames;
  char oldval;
  int i, j;
  _Variable *v;

  /* Register the data to store the output parameter values and
     the final chi2 values, and the convergence flags */
  RegisterScalarData(p, (void *)(&(c->param_outvals)), VARTOOLS_TYPE_DOUBLE,
		     c->Nparams);
  RegisterScalarData(p, (void *)(&(c->chi2out)), VARTOOLS_TYPE_DOUBLE, 0);
  if(c->fittype == VARTOOLS_NONLINFIT_FITTYPE_AMOEBA) {
    RegisterScalarData(p, (void *)(&(c->amoeba_isconverged)), VARTOOLS_TYPE_INT,
		       0);
    RegisterScalarData(p, (void *)(&(c->param_uncertainties)), 
		       VARTOOLS_TYPE_DOUBLE,
		       c->Nparams);
  }

  if(c->uselinfit) {
    RegisterScalarData(p, (void *)(&(c->linfit->param_outvals)), VARTOOLS_TYPE_DOUBLE, c->linfit->Nparams);
    RegisterScalarData(p, (void *)(&(c->linfit->param_uncertainties)), VARTOOLS_TYPE_DOUBLE, c->linfit->Nparams);
  }

  if(c->fittype == VARTOOLS_NONLINFIT_FITTYPE_DEMCMC) {
    RegisterScalarData(p, (void *)(&(c->mcmc_statsout)), VARTOOLS_TYPE_DOUBLE,
		       c->N_mcmc_chain_statstot);
  }
  
  /* Create any new variables, and make sure that old variables are of
     the correct type */
  nparam = c->Nparams;

  for(i=0; i < nparam; i++) {
    for(j=0; j < p->NDefinedVariables; j++) {
      v = p->DefinedVariables[j];
      if(!strcmp(c->paramnames[i],v->varname)) {
	/* This is an existing variable, make sure it is of the correct type */
	if((v->vectortype != VARTOOLS_VECTORTYPE_INLIST &&
	    v->vectortype != VARTOOLS_VECTORTYPE_SCALAR) ||
	   v->datatype != VARTOOLS_TYPE_DOUBLE) {
	  error2(ERR_INVALIDVARIABLEFORNONLINFIT,v->varname);
	}
	c->params[i] = v;
	break;
      }
    }
    if(j == p->NDefinedVariables) {
      /* This is a new variable, create it */
      c->params[i] = SetupScalarVariable(p, c->paramnames[i], VARTOOLS_TYPE_DOUBLE);
    }
  }
  
  if(c->use_covar) {
    c->Corr_amp_var = NULL;
    for(i=0; i < nparam; i++) {
      if(!strcmp(c->Corr_amp_varname, c->paramnames[i])) {
	c->Corr_amp_var = c->params[i];
	if(c->Corr_amp_exprstring[0] != '\0') {
	  error2(ERR_INVALID_PARAMETERVALUE,"the amp_var covariance parameter provided in the -nonlinfit was given an initialization expression, but this variable also appears in the -nonlinfit parameter list. The initialization expression should not be provided in this case.");
	}
	break;
      }
    }
    if(c->Corr_amp_var == NULL) {
      /* Check if this is a new variable */
      for(j=0; j < p->NDefinedVariables; j++) {
	v = p->DefinedVariables[j];
	if(!strcmp(c->Corr_amp_varname,v->varname)) {
	  /* This is an existing variable, make sure it is of the correct type */
	  if((v->vectortype != VARTOOLS_VECTORTYPE_INLIST &&
	      v->vectortype != VARTOOLS_VECTORTYPE_SCALAR) ||
	     v->datatype != VARTOOLS_TYPE_DOUBLE) {
	    error2(ERR_INVALIDVARIABLEFORNONLINFIT,v->varname);
	  }
	  c->Corr_amp_var = v;
	  break;
	}
      }
      if(j == p->NDefinedVariables) {
	/* This is a new variable, create it */
	c->Corr_amp_var = SetupScalarVariable(p, c->Corr_amp_varname, VARTOOLS_TYPE_DOUBLE);
      }
      if(c->Corr_amp_exprstring[0] == '\0') {
	error2(ERR_INVALID_PARAMETERVALUE,"the amp_var covariance parameter provided in the -nonlinfit must be given an initialization expression if it does not appear in the -nonlinfit parameter list.");
      }
    }

    c->Corr_rho_var = NULL;
    for(i=0; i < nparam; i++) {
      if(!strcmp(c->Corr_rho_varname, c->paramnames[i])) {
	c->Corr_rho_var = c->params[i];
	if(c->Corr_rho_exprstring[0] != '\0') {
	  error2(ERR_INVALID_PARAMETERVALUE,"the rho_var covariance parameter provided in the -nonlinfit was given an initialization expression, but this variable also appears in the -nonlinfit parameter list. The initialization expression should not be provided in this case.");
	}
	break;
      }
    }
    if(c->Corr_rho_var == NULL) {
      /* Check if this is a new variable */
      for(j=0; j < p->NDefinedVariables; j++) {
	v = p->DefinedVariables[j];
	if(!strcmp(c->Corr_rho_varname,v->varname)) {
	  /* This is an existing variable, make sure it is of the correct type */
	  if((v->vectortype != VARTOOLS_VECTORTYPE_INLIST &&
	      v->vectortype != VARTOOLS_VECTORTYPE_SCALAR) ||
	     v->datatype != VARTOOLS_TYPE_DOUBLE) {
	    error2(ERR_INVALIDVARIABLEFORNONLINFIT,v->varname);
	  }
	  c->Corr_rho_var = v;
	  break;
	}
      }
      if(j == p->NDefinedVariables) {
	/* This is a new variable, create it */
	c->Corr_rho_var = SetupScalarVariable(p, c->Corr_rho_varname, VARTOOLS_TYPE_DOUBLE);
      }
      if(c->Corr_rho_exprstring[0] == '\0') {
	error2(ERR_INVALID_PARAMETERVALUE,"the rho_var covariance parameter provided in the -nonlinfit must be given an initialization expression if it does not appear in the -nonlinfit parameter list.");
      }
    }

    if(c->covar_type == VARTOOLS_COVARIANCE_KERNEL_MATERN) {
      c->Corr_nu_var = NULL;
      for(i=0; i < nparam; i++) {
	if(!strcmp(c->Corr_nu_varname, c->paramnames[i])) {
	  c->Corr_nu_var = c->params[i];
	  if(c->Corr_nu_exprstring[0] != '\0') {
	    error2(ERR_INVALID_PARAMETERVALUE,"the nu_var covariance parameter provided in the -nonlinfit was given an initialization expression, but this variable also appears in the -nonlinfit parameter list. The initialization expression should not be provided in this case.");
	  }
	  break;
	}
      }
      if(c->Corr_nu_var == NULL) {
	/* Check if this is a new variable */
	for(j=0; j < p->NDefinedVariables; j++) {
	  v = p->DefinedVariables[j];
	  if(!strcmp(c->Corr_nu_varname,v->varname)) {
	    /* This is an existing variable, make sure it is of the correct type */
	    if((v->vectortype != VARTOOLS_VECTORTYPE_INLIST &&
		v->vectortype != VARTOOLS_VECTORTYPE_SCALAR) ||
	       v->datatype != VARTOOLS_TYPE_DOUBLE) {
	      error2(ERR_INVALIDVARIABLEFORNONLINFIT,v->varname);
	    }
	    c->Corr_nu_var = v;
	    break;
	  }
	}
	if(j == p->NDefinedVariables) {
	  /* This is a new variable, create it */
	  c->Corr_nu_var = SetupScalarVariable(p, c->Corr_nu_varname, VARTOOLS_TYPE_DOUBLE);
	}
	if(c->Corr_nu_exprstring[0] == '\0') {
	  error2(ERR_INVALID_PARAMETERVALUE,"the nu_var covariance parameter provided in the -nonlinfit must be given an initialization expression if it does not appear in the -nonlinfit parameter list.");
	}
      }
    }
  }

  if(nparam > 0) {
    if((c->paraminit_expressions = (_Expression **) malloc(nparam * sizeof(_Expression *))) == NULL ||
       (c->paramerr_expressions = (_Expression **) malloc(nparam * sizeof(_Expression *))) == NULL) 
      error(ERR_MEMALLOC);
  }
  if(c->Npriors > 0) {
    if((c->prior_expressions = (_Expression **) malloc(c->Npriors * sizeof(_Expression *))) == NULL)
      error(ERR_MEMALLOC);
  }
  if(c->Nconstraints > 0) {
    if((c->constraint_expressions = (_Expression **) malloc(c->Nconstraints * sizeof(_Expression *))) == NULL)
      error(ERR_MEMALLOC);
  }

  if(c->fittype == VARTOOLS_NONLINFIT_FITTYPE_DEMCMC) {
    if(c->N_mcmc_chain_statstot > 0) {
      if((c->mcmc_chain_stats_expressions = (_Expression **) malloc(c->N_mcmc_chain_expressions * sizeof(_Expression *))) == NULL)
	error(ERR_MEMALLOC);
    }
  }

  /* Initialize the linear fit if needed */
  if(c->uselinfit) {
    InitLinfit(p, c->linfit, cnum);
  }

  /* Set-up the model variable if it is being using */
  if(c->modelvarname != NULL && ! c->uselinfit) {
    for(j=0; j < p->NDefinedVariables; j++) {
      v = p->DefinedVariables[j];
      if(!strcmp(c->modelvarname,v->varname)) {
	/* This is an existing variable, make sure it is the correct type */
	if(v->vectortype != VARTOOLS_VECTORTYPE_LC ||
	   v->datatype != VARTOOLS_TYPE_DOUBLE) {
	  error2(ERR_INVALIDVARIABLEFORNONLINFIT,c->modelvarname);
	}
	c->modelvar = v;
	break;
      }
    }
    if(j == p->NDefinedVariables) {
      /* This is a new variable, create it */
      c->modelvar = CreateVariable(p, c->modelvarname, VARTOOLS_TYPE_DOUBLE,
				   VARTOOLS_VECTORTYPE_LC, NULL);
      RegisterDataFromLightCurve(p,
				c->modelvar->dataptr,
				VARTOOLS_TYPE_DOUBLE,
				 0, 0, 0, 0, 0, NULL,
				 c->modelvar,
				-1, c->modelvarname);
    }
  }
  else if(c->modelvarname != NULL && c->uselinfit) {
    c->modelvar = c->linfit->modelvar;
  }

}

void SetupNonlinfitExpression(ProgramData *p, _Nonlinfit *c)
{
  int i, j;
  if(c->uselinfit) {
    SetupLinfitExpression(p, c->linfit);
  }
  else {
    c->functionexpression = ParseExpression(c->functionstring, p);
  }
  if(c->errorstring != NULL) {
    c->errorexpression = ParseExpression(c->errorstring, p);
  }
  for(i = 0; i < c->Nparams; i++) {
    c->paraminit_expressions[i] = ParseExpression(c->paraminitstrings[i], p);
    c->paramerr_expressions[i] = ParseExpression(c->paramerrstrings[i], p);
  }
  for(i = 0; i < c->Npriors; i++) {
    c->prior_expressions[i] = ParseExpression(c->priorstrings[i], p);
  }
  for(i = 0; i < c->Nconstraints; i++) {
    c->constraint_expressions[i] = ParseExpression(c->constraintstrings[i], p);
  }
  if(c->fittype == VARTOOLS_NONLINFIT_FITTYPE_DEMCMC) {
    for(i=0; i < c->N_mcmc_chain_expressions; i++) {
      c->mcmc_chain_stats_expressions[i] =
	ParseExpression(c->mcmc_chain_expr_strings[i],p);
    }
  }
  if(c->use_covar) {
    if(c->Corr_amp_exprstring[0] != '\0') {
      c->Corr_amp_expr = ParseExpression(c->Corr_amp_exprstring, p);
    } else {
      c->Corr_amp_expr = NULL;
    }
    if(c->Corr_rho_exprstring[0] != '\0') {
      c->Corr_rho_expr = ParseExpression(c->Corr_rho_exprstring, p);
    } else {
      c->Corr_rho_expr = NULL;
    }
    if(c->covar_type == VARTOOLS_COVARIANCE_KERNEL_MATERN) {
      if(c->Corr_nu_exprstring[0] != '\0') {
	c->Corr_nu_expr = ParseExpression(c->Corr_nu_exprstring, p);
      } else {
	c->Corr_nu_expr = NULL;
      }
    }
  }
}


int nonlinfit_parse_paraminit_string(char *ins, char *s1, char *s2, char *s3) {
  int i, j;
  char *cpystring;
  if((cpystring = (char *) malloc((strlen(ins)+1))) == NULL) {
    error(ERR_MEMALLOC);
  }
  sprintf(cpystring,"%s",ins);
  i=0; j=0;
  while(cpystring[j] != '\0' && cpystring[j] != '=') j++;
  if(!j) {
    free(cpystring);
    return 0;
  }
  cpystring[j] = '\0';
  sprintf(s1,"%s",&(cpystring[i]));
  if(ins[j] == '\0') {
    free(cpystring);
    return 1;
  }
  j++;
  i=j;
  while(cpystring[j] != '\0' && cpystring[j] != ':') j++;
  if(j == i) {
    free(cpystring);
    return 1;
  }
  cpystring[j] = '\0';;
  sprintf(s2,"%s",&(cpystring[i]));
  if(ins[j] == '\0') {
    free(cpystring);
    return 2;
  }
  j++;
  i=j;
  while(cpystring[j] != '\0') j++;
  if(j == i) {
    free(cpystring);
    return 1;
  }
  cpystring[j] = '\0';
  sprintf(s3,"%s",&(cpystring[i]));
  free(cpystring);
  return 3;
}

int ParseNonlinfitCommand(int *iret, int argc, char **argv, ProgramData *p,
			  _Nonlinfit *c)
{
  int i, k, j, ii, nparam, s, nscan, Nstat;
  char oldval;

  char *tmpstring;
  char **linfit_args;
  int lentmp = 256;
  int i1, i2;
  double pctval;
  int Npct = 0;

  i = *iret;

  if(i >= argc)
    return(1);

  if((tmpstring = (char *) malloc(lentmp)) == NULL)
    error(ERR_MEMALLOC);

  c->priorliststring = NULL;
  c->constraintliststring = NULL;
  c->modelvarname = NULL;
  c->correctlc = 0;
  c->omodel = 0;
  c->fittype = VARTOOLS_NONLINFIT_FITTYPE_AMOEBA;
  c->uselinfit = 0;
  c->Npriors = 0;
  c->Nconstraints = 0;

  c->mcmc_chain_exprliststring = NULL;
  c->mcmc_chain_statsliststring = NULL;

  c->N_mcmc_chain_expressions = 0;
  c->N_mcmc_chain_stats = 0;
  c->N_mcmc_chain_statstot = 0;
  c->mcmc_skipamoeba = 0;

  c->outfile_extension = (char *) malloc(20);
  sprintf(c->outfile_extension,"nonlinfit.model");
  c->outfilename_format = NULL;
  c->mcmc_outchains_format = NULL;


  if((c->functionstring = (char *) malloc((strlen(argv[i])+1))) == NULL)
    error(ERR_MEMALLOC);
  sprintf(c->functionstring,"%s",argv[i]);

  i++;
  if(i >= argc) {*iret = i; return 1;}
  if((c->paramliststring = (char *) malloc((strlen(argv[i])+1))) == NULL)
    error(ERR_MEMALLOC);
  sprintf(c->paramliststring,"%s",argv[i]);

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"linfit")) {
      c->uselinfit = 1;
      i++;
      if(i >= argc) {*iret = i; return 1;}
      /* Create a set of options to the linfit command */
      if((linfit_args = (char **) malloc(4*sizeof(char *))) == NULL)
	error(ERR_MEMALLOC);
      if((linfit_args[0] = (char *) malloc((strlen(c->functionstring)+1))) == NULL)
	error(ERR_MEMALLOC);
      if((linfit_args[1] = (char *) malloc((strlen(argv[i])+1))) == NULL)
	error(ERR_MEMALLOC);
      if((linfit_args[2] = (char *) malloc(9)) == NULL)
	error(ERR_MEMALLOC);
      sprintf(linfit_args[0],"%s",c->functionstring);
      sprintf(linfit_args[1],"%s",argv[i]);
    }
    else
      i--;
  } else
    i--;

  c->errorstring = NULL;
  c->use_covar = 0;
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"errors")) {
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if((c->errorstring = (char *) malloc((strlen(argv[i])+1))) == NULL)
	error(ERR_MEMALLOC);
      sprintf(c->errorstring,"%s",argv[i]);
    } else
      i--;
  } else
    i--;

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"covariance")) {
      c->use_covar = 1;
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if(!strcmp(argv[i],"squareexp")) {
	c->covar_type = VARTOOLS_COVARIANCE_KERNEL_SQUAREDEXPONENTIAL;
      }
      else if(!strcmp(argv[i],"exp")) {
	c->covar_type = VARTOOLS_COVARIANCE_KERNEL_EXPONENTIAL;
      }
      else if(!strcmp(argv[i],"matern")) {
	c->covar_type = VARTOOLS_COVARIANCE_KERNEL_MATERN;
      }
      else {
	*iret = i; return 1;
      }
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if((c->Corr_amp_varname = (char *) malloc((strlen(argv[i])+1))) == NULL)
	error(ERR_MEMALLOC);
      if((c->Corr_amp_exprstring = (char *) malloc((strlen(argv[i])+1))) == NULL)
	error(ERR_MEMALLOC);
      nscan = nonlinfit_parse_paraminit_string(argv[i],c->Corr_amp_varname,c->Corr_amp_exprstring,c->Corr_amp_exprstring);
      if(nscan != 1 && nscan != 2) {
	*iret = i; return 1;
      }
      if(nscan == 1) c->Corr_amp_exprstring[0]='\0';
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if((c->Corr_rho_varname = (char *) malloc((strlen(argv[i])+1))) == NULL)
	error(ERR_MEMALLOC);
      if((c->Corr_rho_exprstring = (char *) malloc((strlen(argv[i])+1))) == NULL)
	error(ERR_MEMALLOC);
      nscan = nonlinfit_parse_paraminit_string(argv[i],c->Corr_rho_varname,c->Corr_rho_exprstring,c->Corr_rho_exprstring);
      if(nscan != 1 && nscan != 2) {
	*iret = i; return 1;
      }
      if(nscan == 1) c->Corr_rho_exprstring[0]='\0';
      if(c->covar_type == VARTOOLS_COVARIANCE_KERNEL_MATERN) {
	i++;
	if(i >= argc) {*iret = i; return 1;}
	if((c->Corr_nu_varname = (char *) malloc((strlen(argv[i])+1))) == NULL)
	  error(ERR_MEMALLOC);
	if((c->Corr_nu_exprstring = (char *) malloc((strlen(argv[i])+1))) == NULL)
	  error(ERR_MEMALLOC);
	nscan = nonlinfit_parse_paraminit_string(argv[i],c->Corr_nu_varname,c->Corr_nu_exprstring,c->Corr_nu_exprstring);
	if(nscan != 1 && nscan != 2) {
	  *iret = i; return 1;
	}
	if(nscan == 1) c->Corr_nu_exprstring[0]='\0';
      }
    } else
      i--;
  } else
    i--;

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"priors")) {
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if((c->priorliststring = (char *) malloc((strlen(argv[i])+1))) == NULL)
	error(ERR_MEMALLOC);
      sprintf(c->priorliststring,"%s",argv[i]);
    } else
      i--;
  }
  else
    i--;

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"constraints")) {
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if((c->constraintliststring = (char *) malloc((strlen(argv[i])+1))) == NULL)
	error(ERR_MEMALLOC);
      sprintf(c->constraintliststring,"%s",argv[i]);
    } else
      i--;
  }
  else
    i--;

  i++;
  if(i >= argc) {*iret = i; return 1;}
  if(!strcmp(argv[i],"amoeba")) {
    c->fittype = VARTOOLS_NONLINFIT_FITTYPE_AMOEBA;
    c->amoeba_tol = VARTOOLS_NONLINFIT_DEFAULT_AMOEBA_TOLERANCE;
    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"tolerance")) {
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->amoeba_tol = atof(argv[i]);
      } else
	i--;
    } else
      i--;
    c->amoeba_maxsteps = VARTOOLS_NONLINFIT_DEFAULT_AMOEBA_MAXSTEPS;
    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"maxsteps")) {
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->amoeba_maxsteps = atol(argv[i]);
      } else
	i--;
    } else
      i--;
  }
  else if(!strcmp(argv[i],"mcmc")) {
    c->fittype = VARTOOLS_NONLINFIT_FITTYPE_DEMCMC;
    c->mcmc_Naccept = -1;
    c->mcmc_Nlinkstotal = VARTOOLS_NONLINFIT_DEFAULT_MCMC_NLINKSTOTAL;
    c->mcmc_burninfrac = VARTOOLS_NONLINFIT_DEFAULT_MCMC_BURNIN;
    c->mcmc_eps = VARTOOLS_NONLINFIT_DEFAULT_MCMC_EPS;
    c->mcmc_chain_exprliststring = NULL;
    c->mcmc_chain_statsliststring = NULL;
    c->mcmc_max_mem_store = VARTOOLS_NONLINFIT_DEFAULT_MAX_MEM_STORE;
    c->mcmc_outchains = 0;
    c->mcmc_outchains_print_every = 1;
    c->amoeba_maxsteps = VARTOOLS_NONLINFIT_DEFAULT_AMOEBA_MAXSTEPS;
    c->amoeba_tol = VARTOOLS_NONLINFIT_DEFAULT_AMOEBA_TOLERANCE;
    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"Naccept")) {
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->mcmc_Naccept = atol(argv[i]);
      } 
      else if(!strcmp(argv[i],"Nlinkstotal")) {
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->mcmc_Nlinkstotal = atol(argv[i]);
      }
      else
	i--;
    }
    else 
      i--;
    
    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"fracburnin")) {
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->mcmc_burninfrac = atof(argv[i]);
      } else
	i--;
    } else
      i--;
    
    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"eps")) {
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->mcmc_eps = atof(argv[i]);
      } else
	i--;
    } else
      i--;
    
    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"skipamoeba")) {
	c->mcmc_skipamoeba = 1;
      }
      else
	i--;
    }
    else
      i--;

    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"chainstats")) {
	i++;
	if(i >= argc) {*iret = i; return 1;}
	if((c->mcmc_chain_exprliststring = (char *) malloc((strlen(argv[i])+1))) == NULL)
	  error(ERR_MEMALLOC);
	sprintf(c->mcmc_chain_exprliststring,"%s",argv[i]);
	i++;
	if(i >= argc) {*iret = i; return 1;}
	if((c->mcmc_chain_statsliststring = (char *) malloc((strlen(argv[i])+1))) == NULL)
	  error(ERR_MEMALLOC);
	sprintf(c->mcmc_chain_statsliststring,"%s",argv[i]);
      } else
	i--;
    } else
      i--;
    
    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"maxmemstore")) {
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->mcmc_max_mem_store = atof(argv[i]);
	if(c->mcmc_max_mem_store <= 0)
	  error(ERR_MCMCINVALIDMAXMEM);
      } else
	i--;
    } else
      i--;

    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"outchains")) {
	c->mcmc_outchains = 1;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	if((c->mcmc_outchains_dir = (char *) malloc((strlen(argv[i])+1))) == NULL)
	  error(ERR_MEMALLOC);
	sprintf(c->mcmc_outchains_dir,"%s",argv[i]);

	if((c->mcmc_outchains_format = (char *) malloc(1)) == NULL)
	  error(ERR_MEMALLOC);
	c->mcmc_outchains_format[0] = '\0';
	  
	/* Check if the user gave the "format" keyword */
	i++;
	if(i < argc) {
	  if(!strcmp(argv[i],"format")) {
	    i++;
	    if(i < argc) {
	      c->mcmc_outchains_format = (char *) realloc(c->mcmc_outchains_format, (strlen(argv[i])+1)*sizeof(char));
	      sprintf(c->mcmc_outchains_format,"%s",argv[i]);
	    } else {
	      *iret = i; return 1;
	    }
	  } else
	    i--;
	} else
	  i--;

	i++;
	if(i < argc) {
	  if(!strcmp(argv[i],"printevery")) {
	    i++;
	    if(i >= argc) {*iret = i; return 1;}
	    c->mcmc_outchains_print_every = atoi(argv[i]);
	  } else
	    i--;
	} else
	  i--;
      } else
	i--;
    }
    else
      i--;
  }
  else {*iret = i; return 1;}
  
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"modelvar")) {
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if((c->modelvarname = (char *) malloc((strlen(argv[i])+1))) == NULL)
	error(ERR_MEMALLOC);
      sprintf(c->modelvarname,"%s",argv[i]);
    } else
      i--;
  }
  else
    i--;

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"correctlc")) {
      c->correctlc = 1;
    } else i--;
  } else i--;
  
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"omodel")) {
      c->omodel = 1;
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if((c->outdir = (char *) malloc((strlen(argv[i])+1))) == NULL)
	error(ERR_MEMALLOC);
      sprintf(c->outdir,"%s",argv[i]);
      /* Check if the user gave the "format" keyword */
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"format")) {
	  i++;
	  if(i < argc) {
	    c->outfilename_format = (char *) realloc(c->outfilename_format, (strlen(argv[i])+1)*sizeof(char));
	    sprintf(c->outfilename_format,"%s",argv[i]);
	  } else {
	    *iret = i; return 1;
	  }
	} else
	  i--;
      } else
	i--;
    } else
      i--;
  } else
    i--;

  if(c->uselinfit && c->use_covar) {
    error2(ERR_INVALID_PARAMETERVALUE,"-nonlinfit, you cannot use both linear optimization and a non-diagonal covariance matrix at the same time.\n");
  }

  /* Parse the parameter string to get the number of output parameters,
     this is needed for the CreateOutputColumns function. Later
     CompileAllExpressions will create the variables. */
  k = 0; j = 0; nparam = 0;
  do {
    if(c->paramliststring[k] == '\0' || c->paramliststring[k] == ',') {
      oldval = c->paramliststring[k];
      c->paramliststring[k] = '\0';
      if(!nparam) {
	if((c->paramnames = (char **) malloc(sizeof(char *))) == NULL ||
	   (c->paraminitstrings = (char **) malloc(sizeof(char *))) == NULL ||
	   (c->paramerrstrings = (char **) malloc(sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
      } else {
	if((c->paramnames = (char **) realloc(c->paramnames, (nparam + 1)*sizeof(char *))) == NULL ||
	   (c->paraminitstrings = (char **) realloc(c->paraminitstrings, (nparam + 1)*sizeof(char *))) == NULL ||
	   (c->paramerrstrings = (char **) realloc(c->paramerrstrings, (nparam + 1)*sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
      }
      if((s = strlen(&(c->paramliststring[j]))) == 0)
	error2(ERR_BADVARIABLENAME,"");
      if((c->paramnames[nparam] = (char *) malloc((s+1))) == NULL ||
	 (c->paraminitstrings[nparam] = (char *) malloc((s+1))) == NULL ||
	 (c->paramerrstrings[nparam] = (char *) malloc((s+1))) == NULL)
	error(ERR_MEMALLOC);
      nscan = nonlinfit_parse_paraminit_string(&(c->paramliststring[j]),c->paramnames[nparam],c->paraminitstrings[nparam],c->paramerrstrings[nparam]);
      if(nscan != 3) {
	error2(ERR_BADNONLINFITPARAMINIT,&(c->paramliststring[j]));
      }
      j = k+1;
      c->paramliststring[k] = oldval;
      nparam++;
    }
    k++;
  } while(c->paramliststring[k-1] != '\0');
  c->Nparams = nparam;
  if((c->params = (_Variable **) malloc(nparam * sizeof(_Variable *))) == NULL)
    error(ERR_MEMALLOC);

  /* Parse the prior list string */
  if(c->priorliststring != NULL) {
    k = 0; j = 0; nparam = 0;
    do {
      if(c->priorliststring[k] == '\0' || c->priorliststring[k] == ',') {
	oldval = c->priorliststring[k];
	c->priorliststring[k] = '\0';
	if(!nparam) {
	  if((c->priorvarnames = (char **) malloc(sizeof(char *))) == NULL ||
	     (c->priorstrings = (char **) malloc(sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	} else {
	  if((c->priorvarnames = (char **) realloc(c->priorvarnames, (nparam + 1)*sizeof(char *))) == NULL ||
	     (c->priorstrings = (char **) realloc(c->priorstrings, (nparam + 1)*sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	}
	if((s = strlen(&(c->priorliststring[j]))) == 0)
	  error2(ERR_BADVARIABLENAME,"");
	if((c->priorvarnames[nparam] = (char *) malloc((s+1))) == NULL ||
	   (c->priorstrings[nparam] = (char *) malloc((s+1))) == NULL)
	  error(ERR_MEMALLOC);
	nscan = nonlinfit_parse_paraminit_string(&(c->priorliststring[j]),c->priorvarnames[nparam],c->priorstrings[nparam],c->priorstrings[nparam]);
	if(nscan != 2) {
	  error2(ERR_BADNONLINFITPRIORINIT,&(c->priorliststring[j]));
	}
	/* Make sure that the variable given in the prior string is actually
	   one of the variables fromt he parameter list */
	for(ii = 0; ii < c->Nparams; ii++) {
	  if(!strcmp(c->paramnames[ii],c->priorvarnames[nparam]))
	    break;
	}
	if(ii >= c->Nparams) {
	  error2(ERR_BADNONLINFITPRIORINIT,&(c->priorliststring[j]));
	}	  
	j = k+1;
	c->priorliststring[k] = oldval;
	nparam++;
      }
      k++;
    } while(c->priorliststring[k-1] != '\0');
    c->Npriors = nparam;
  }

  /* Parse the constraint list string */
  if(c->constraintliststring != NULL) {
    k = 0; j = 0; nparam = 0;
    do {
      if(c->constraintliststring[k] == '\0' || c->constraintliststring[k] == ',') {
	oldval = c->constraintliststring[k];
	c->constraintliststring[k] = '\0';
	if(!nparam) {
	  if((c->constraintstrings = (char **) malloc(sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	} else {
	  if((c->constraintstrings = (char **) realloc(c->constraintstrings, (nparam + 1)*sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	}
	if((s = strlen(&(c->constraintliststring[j]))) == 0)
	  error2(ERR_BADVARIABLENAME,"");
	if((c->constraintstrings[nparam] = (char *) malloc((s+1))) == NULL)
	  error(ERR_MEMALLOC);
	sprintf(c->constraintstrings[nparam],"%s",&(c->constraintliststring[j]));
	j = k+1;
	c->constraintliststring[k] = oldval;
	nparam++;
      }
      k++;
    } while(c->constraintliststring[k-1] != '\0');
    c->Nconstraints = nparam;
  }
  if(c->use_covar) {
    if(c->covar_type == VARTOOLS_COVARIANCE_KERNEL_MATERN) {
      if(c->Corr_nu_exprstring[0] == '\0') {
	if(!c->Nconstraints) {
	  if((c->constraintstrings = (char **) malloc(sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	} else {
	  if((c->constraintstrings = (char **) realloc(c->constraintstrings,(c->Nconstraints+1)*sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	}
	s = strlen(c->Corr_nu_varname) + 3;
	if((c->constraintstrings[c->Nconstraints] = (char *) malloc(s*sizeof(char))) == NULL)
	  error(ERR_MEMALLOC);
	sprintf(c->constraintstrings[c->Nconstraints],"%s>0",c->Corr_nu_varname);
	c->Nconstraints += 1;
      }
    }
    if(c->Corr_amp_exprstring[0] == '\0') {
      if(!c->Nconstraints) {
	if((c->constraintstrings = (char **) malloc(sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
      } else {
	if((c->constraintstrings = (char **) realloc(c->constraintstrings,(c->Nconstraints+1)*sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
      }
      s = strlen(c->Corr_amp_varname) + 3;
      if((c->constraintstrings[c->Nconstraints] = (char *) malloc(s*sizeof(char))) == NULL)
	error(ERR_MEMALLOC);
      sprintf(c->constraintstrings[c->Nconstraints],"%s>0",c->Corr_amp_varname);
      c->Nconstraints += 1;
    }
    if(c->Corr_rho_exprstring[0] == '\0') {
      if(!c->Nconstraints) {
	if((c->constraintstrings = (char **) malloc(sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
      } else {
	if((c->constraintstrings = (char **) realloc(c->constraintstrings,(c->Nconstraints+1)*sizeof(char *))) == NULL)
	  error(ERR_MEMALLOC);
      }
      s = strlen(c->Corr_rho_varname) + 3;
      if((c->constraintstrings[c->Nconstraints] = (char *) malloc(s*sizeof(char))) == NULL)
	error(ERR_MEMALLOC);
      sprintf(c->constraintstrings[c->Nconstraints],"%s>0",c->Corr_rho_varname);
      c->Nconstraints += 1;
    }
  }

  /* Set the default list of statistics to be reported for MCMC fitting if
     this has not already been specified */
  if(c->fittype == VARTOOLS_NONLINFIT_FITTYPE_DEMCMC && c->mcmc_chain_exprliststring == NULL) {
    nparam = 0;
    for(j = 0; j < c->Nparams; j++) {
      nparam += (strlen(c->paramnames[j]) + 1);
    }
    nparam++;
    if((c->mcmc_chain_exprliststring = (char *) malloc(nparam)) == NULL) {
      error(ERR_MEMALLOC);
    }
    c->mcmc_chain_exprliststring[0] = '\0';
    for(j=0; j < c->Nparams; j++) {
      if(j > 0) {
	sprintf(c->mcmc_chain_exprliststring,"%s,",c->mcmc_chain_exprliststring);
      }
      sprintf(c->mcmc_chain_exprliststring,"%s%s",c->mcmc_chain_exprliststring,c->paramnames[j]);
    }
  }
  if(c->fittype == VARTOOLS_NONLINFIT_FITTYPE_DEMCMC && c->mcmc_chain_statsliststring == NULL) {
    if((c->mcmc_chain_statsliststring = (char *) malloc(14)) == NULL)
      error(ERR_MEMALLOC);
    sprintf(c->mcmc_chain_statsliststring,"median,stddev");
  }


  /* Parse the mcmc chain expression list */
  if(c->mcmc_chain_exprliststring != NULL) {
    k = 0; j = 0; nparam = 0;
    do {
      if(c->mcmc_chain_exprliststring[k] == '\0' || c->mcmc_chain_exprliststring[k] == ',') {
	oldval = c->mcmc_chain_exprliststring[k];
	c->mcmc_chain_exprliststring[k] = '\0';
	if(!nparam) {
	  if((c->mcmc_chain_expr_strings = (char **) malloc(sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	} else {
	  if((c->mcmc_chain_expr_strings = (char **) realloc(c->mcmc_chain_expr_strings, (nparam + 1)*sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	}
	if((s = strlen(&(c->mcmc_chain_exprliststring[j]))) == 0)
	  error2(ERR_BADVARIABLENAME,"");
	if((c->mcmc_chain_expr_strings[nparam] = (char *) malloc((s+1))) == NULL)
	  error(ERR_MEMALLOC);
	sprintf(c->mcmc_chain_expr_strings[nparam],"%s",&(c->mcmc_chain_exprliststring[j]));
	j = k+1;
	c->mcmc_chain_exprliststring[k] = oldval;
	nparam++;
      }
      k++;
    } while(c->mcmc_chain_exprliststring[k-1] != '\0');
    c->N_mcmc_chain_expressions = nparam;
  }

  if(c->mcmc_chain_statsliststring != NULL) {
    
    Nstat = 1;
    j = 0;
    while(c->mcmc_chain_statsliststring[j] != '\0') {
      if(c->mcmc_chain_statsliststring[j] == ',')
	Nstat++;
      j++;
    }
    if((c->mcmc_statstocalc = (int *) malloc(Nstat * sizeof(int))) == NULL)
      error(ERR_MEMALLOC);
    i1 = 0;
    i2 = 0;
    for(k = 0; k < Nstat; k++) {
      while(c->mcmc_chain_statsliststring[i2] != '\0' && c->mcmc_chain_statsliststring[i2] != ',') {
	i2++;
      }

      if((i2 - i1 + 1) > lentmp) {
	lentmp *= 2;
	if((tmpstring = (char *) realloc(tmpstring, lentmp*sizeof(char))) == NULL)
	  error(ERR_MEMALLOC);
      }
      for(j=i1; j < i2; j++) {
	tmpstring[j-i1] = c->mcmc_chain_statsliststring[j];
      }
      tmpstring[j-i1] = '\0';

      if(!strcmp(tmpstring,"mean")) {
	c->mcmc_statstocalc[k] = VARTOOLS_STATSTYPE_MEAN;
      }
      else if(!strcmp(tmpstring,"weightedmean")) {
	c->mcmc_statstocalc[k] = VARTOOLS_STATSTYPE_WEIGHTEDMEAN;
      }
      else if(!strcmp(tmpstring,"median")) {
	c->mcmc_statstocalc[k] = VARTOOLS_STATSTYPE_MEDIAN;
      }
      else if(!strcmp(tmpstring,"wmedian")) {
	c->mcmc_statstocalc[k] = VARTOOLS_STATSTYPE_MEDIAN_WEIGHT;
      }
      else if(!strcmp(tmpstring,"stddev")) {
	c->mcmc_statstocalc[k] = VARTOOLS_STATSTYPE_STDDEV;
      }
      else if(!strcmp(tmpstring,"meddev")) {
	c->mcmc_statstocalc[k] = VARTOOLS_STATSTYPE_MEDDEV;
      }
      else if(!strcmp(tmpstring,"medmeddev")) {
	c->mcmc_statstocalc[k] = VARTOOLS_STATSTYPE_MEDMEDDEV;
      }
      else if(!strcmp(tmpstring,"MAD")) {
	c->mcmc_statstocalc[k] = VARTOOLS_STATSTYPE_MAD;
      }
      else if(!strcmp(tmpstring,"kurtosis")) {
	c->mcmc_statstocalc[k] = VARTOOLS_STATSTYPE_KURTOSIS;
      }
      else if(!strcmp(tmpstring,"skewness")) {
	c->mcmc_statstocalc[k] = VARTOOLS_STATSTYPE_SKEWNESS;
      }
      else if(sscanf(tmpstring,"pct%lf",&pctval) == 1) {
	c->mcmc_statstocalc[k] = VARTOOLS_STATSTYPE_PERCENTILE;
	if(!Npct) {
	  c->pctval = (double *) malloc(sizeof(double));
	} else {
	  c->pctval = (double *) realloc(c->pctval, (Npct + 1)*sizeof(double));
	}
	c->pctval[Npct] = pctval;
	Npct++;
      }
      else if(sscanf(tmpstring,"wpct%lf",&pctval) == 1) {
	c->mcmc_statstocalc[k] = VARTOOLS_STATSTYPE_PERCENTILE_WEIGHT;
	if(!Npct) {
	  c->pctval = (double *) malloc(sizeof(double));
	} else {
	  c->pctval = (double *) realloc(c->pctval, (Npct + 1)*sizeof(double));
	}
	c->pctval[Npct] = pctval;
	Npct++;
      }
      else if(!strcmp(tmpstring,"max")) {
	c->mcmc_statstocalc[k] = VARTOOLS_STATSTYPE_MAXIMUM;
      }
      else if(!strcmp(tmpstring,"min")) {
	c->mcmc_statstocalc[k] = VARTOOLS_STATSTYPE_MINIMUM;
      }
      else if(!strcmp(tmpstring,"sum")) {
	c->mcmc_statstocalc[k] = VARTOOLS_STATSTYPE_SUM;
      }
      else {
	error2(ERR_INVALIDSTATISTIC,tmpstring);
      }

      i1 = i2+1;
      i2 = i2+1;
    }

    free(tmpstring);

    c->N_mcmc_chain_stats = Nstat;

    c->N_mcmc_chain_statstot = c->N_mcmc_chain_expressions * c->N_mcmc_chain_stats;
  }

  /* Setup the linfit structure if it is being used */
  if(c->uselinfit) {
    if((c->linfit = (_Linfit *) malloc(sizeof(_Linfit))) == NULL)
      error(ERR_MEMALLOC);
    if(c->modelvarname != NULL) {
      sprintf(linfit_args[2],"modelvar");
      linfit_args[3] = c->modelvarname;
      j = 0;
      if(ParseLinfitCommand(&j, 4, linfit_args, p, c->linfit)) {
	*iret = i; return 1;}
    }
    else {
      j = 0;
      if(ParseLinfitCommand(&j, 2, linfit_args, p, c->linfit)) {
	*iret = i; return 1;}
    }
    c->linfit->calcchi2out = 1;
    free(linfit_args[2]);
    free(linfit_args[1]);
    free(linfit_args[0]);
    free(linfit_args);
  }


  *iret = i;
  return 0;
}

