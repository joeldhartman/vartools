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
#include <stdio.h>
#include <stdlib.h>

int VARTOOOLS_ParseParameter(ProgramData *p,
			    Command *c,
			    int *iret,
			    char **argv,
			    int argc,
			    const char *keyword,
			    int Nvec,
			    ...);

int VARTOOLS_ParseOutNameKeyword(ProgramData *p,
				 Command *c,
				 int *iret,
				 char **argv,
				 int argc,
				 const char *keyword,
				 int *outputflag,
				 char **outdir,
				 int *formatflag,
				 char **format);

int VARTOOLS_ParseConstantParameter(ProgramData *p,
				    Command *c,
				    int *iret,
				    char **argv,
				    int argc,
				    const char *keyword,
				    int datatype,
				    void *dataptr,
				    int Ncolumns);

int VARTOOLS_ParseFixSpecFixcolumn(ProgramData *p,
				  Command *c,
				  int *iret,
				  char **argv,
				  int argc,
				  int Nvec,
				  ...);

int VARTOOLS_amoeba(double **, double *, int *, int, double, double (*)(double *, int, int, double *, double *, double *, void *), int *, int, int, double *, double *, double *, void *);

#ifdef PARALLEL
void VARTOOLS_mrqmin(double *x, double *y, double *sig, int ndata, double *a, int *ia, int ma, double **covar, double **alpha, double *chisq, double *alamda, void (*funcs)(double *, double *, double *, double **, int, int, int, double *, double **, double *, double *, int *, void *), int Nlin_coeff, double **Design_Matrix, double *lin_coeffs, int *varylin_coeffs, int mfit, double *ochisq, double *atry, double *beta, double *da, double **oneda, void *userparams);
#else
void VARTOOLS_mrqmin(double *x, double *y, double *sig, int ndata, double *a, int *ia, int ma, double **covar, double **alpha, double *chisq, double *alamda, void (*funcs)(double *, double *, double *, double **, int, int, int, double *, double **, double *, double *, int *, void *), int Nlin_coeff, double **Design_Matrix, double *lin_coeffs, int *varylin_coeffs, void *userparams);
#endif


void VARTOOLS_incrementparameters_foramoeba(int *Nparameters, int *Ntovary, double ***p, int **ia, int varyparam, double initialval, double stepval);

void VARTOOLS__amoeba_initializesimplexchi2(int Nparameters, int Ntovary, double **p, double **chi2vals, double (*funk)(double *, int, int, double *, double *, double *, void *), int N, double *t, double *mag, double *err, void * userparam);

void VARTOOLS_amoeba_cleanup(int *Nparameters, int *Ntovary, double ***p, int **ia, double **chi2vals);

void VARTOOLS_RegisterDataVector(ProgramData *p, Command *c, void *dataptr,
				 int datatype, int Ncolumns, int source,
				 int output, char *outname, ...);

void VARTOOLS_GetOutputFilename(char *lcoutname, char *lcname, char *outdir,
				char *suffix, char *format, int lc_name_num);

void VARTOOLS_integratemandelagoltransitmodel(double exptime_phase, int Npoints,
					      double *phase, double *outlc, 
					      int type, double *ldcoeffs, 
					      double sin_i, double a, double e,
					      double p, double omega, 
					      int Nresamp);

void VARTOOLS_mandelagoltransitmodel(int Npoints, double *phase, double *outlc,
				     int type, double *ldcoeffs, double sin_i,
				     double a, double e, double p, 
				     double omega);

void VARTOOLS_spline(double *x,double *y,int n,double yp1,double ypn,
		     double *y2,double *u);

void VARTOOLS_splint(double *xa,double *ya,double *y2a,int n,double x,
		     double *y);

void VARTOOLS_spline_monotonic(int N, double *x, double *y, double *yprime);

double VARTOOLS_splint_monotonic(int N, double *x, double *y, double *yprime, 
				 double xt);

void VARTOOLS_medianfilter(int N, double *t, double *mag, double *sig, 
			   double timesize, int meanflag, int replace);

double VARTOOLS_getweightedmean(int n, double *data, double *sig);

double VARTOOLS_getmean(int n, double *data);

double VARTOOLS_median(int n, double *data);

double VARTOOLS_MAD(int n, double *data);

double VARTOOLS_stddev(int n, double *data);

double VARTOOLS_kurtosis(int n, double *data);

double VARTOOLS_skewness(int n, double *data);

double VARTOOLS_percentile(int n, double *data, double pct);

void VARTOOLS_error(int errflag);

void VARTOOLS_error2(int errflag, char *s);

double VARTOOLS_fitpoly(int N, double *x, double *y, double *sig, int order,
			int subtractfit, double *fitparams, double *paramerrs);

double VARTOOLS_chi2(int N, double *t, double *mag, double *err, 
		     double *weighted_average, int *Ngood);

int VARTOOLS_isDifferentPeriods (double period1, double period2, double TimeSpan);


void VARTOOLS_sort_generic(int N, int isreverse, int *index, int Nms, ...);

void VARTOOLS_sortvec_double(int N, double *data1);

int VARTOOLS_RegisterUserFunction(ProgramData *, char *, int, double (*)(double *));

void VARTOOLS_occultquad(double *z0, double u1, double u2, double p, double *muo1, double *mu0, int nz);
