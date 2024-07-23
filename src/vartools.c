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
#include "vartools.h"
#include "vartools_functionpointers.h"

/*
typedef struct {
  int (*ParseFixSpecFixcolumn)(ProgramData *, Command *, int *,
			       char **, int, int, _ParseFixSpecFixcolumnStruct *);
  int (*amoeba)(double **, double *, int *, int, double, double (*)(double *, int, int, double *, double *, double *, void *), int *, int, double *, double *, double *, void *);
  void (*vRegisterDataVector)(ProgramData *, Command *, void *, int, int, int,
			     int, char *, va_list argp);
  void (*GetOutputFilename)(char *, char *, char *, char *, char *, int);
  void (*incrementparameters_foramoeba)(int *, int *, double ***, int **, int, double, double);
  void (*integratemandelagoltransitmodel)(double, int, double *, double *, int, double *, double, double, double, double, double, int);
  void (*mandelagoltransitmodel)(int, double *, double *, int, double *, double,
 double, double, double, double);
  void (*spline)(double *,double *,int,double,double,double *,double *);
  void (*splint)(double *,double *,double *,int,double,double *);
  void (*spline_monotonic)(int, double *, double *, double *);
  double (*splint_monotonic)(int, double *, double *, double *, double);
  void (*medianfilter)(int, double *, double *, double *, double, int, int);
  double (*getweightedmean)(int, double *, double *);
  double (*getmean)(int, double *);
  double (*median)(int, double *);
  double (*MAD)(int, double *);
  double (*stddev)(int, double *);
  double (*kurtosis)(int, double *);
  double (*skewness)(int, double *);  
} _VARTOOLS_FUNCTION_POINTER_STRUCT;
*/

_VARTOOLS_FUNCTION_POINTER_STRUCT VARTOOLS_FUNCTION_POINTER_STRUCT;

void VARTOOLS_Set_Function_Pointers(
				    void (*Set_Function_Pointers_Callback_)
				    (_VARTOOLS_FUNCTION_POINTER_STRUCT *)) {
  Set_Function_Pointers_Callback_(&VARTOOLS_FUNCTION_POINTER_STRUCT);
}

/*
void VARTOOLS_Set_Function_Pointers(int (*ParseFixSpecFixcolumn_)(ProgramData *,
								 Command *,
								 int *,
								 char **,
								 int,
								 int,
								  _ParseFixSpecFixcolumnStruct *),
				    int (*amoeba_)(double **, double *, int *, int, double, double (*)(double *, int, int, double *, double *, double *, void *), int *, int, double *, double *, double *, void *),
				    void (*vRegisterDataVector_)(ProgramData *,
								Command *,
								void *,
								int,
								int,
								int,
								int,
								char *,
								 va_list argp),
				    void (*GetOutputFilename_)(char *,
							       char *,
							       char *,
							       char *,
							       char *,
							       int),
				    void (*incrementparameters_foramoeba_)(int *,
									   int *,
									   double ***,
									   int **,
									   int,
									   double,
									   double),
				    void (*integratemandelagoltransitmodel_)(double, int, double *, double *, int, double *, double, double, double, double, double, int),
				    void (*mandelagoltransitmodel_)(int, double *, double *, int, double *, double,
 double, double, double, double),
				    void (*spline_)(double *,double *,int,double,double,double *,double *),
				    void (*splint_)(double *,double *,double *,int,double,double *),
				    void (*spline_monotonic_)(int, double *, double *, double *),
				    double (*splint_monotonic_)(int, double *, double *, double *, double),
				    void (*medianfilter_)(int, double *, double *, double *, double, int, int)){
  VARTOOLS_FUNCTION_POINTER_STRUCT.ParseFixSpecFixcolumn = ParseFixSpecFixcolumn_;
  VARTOOLS_FUNCTION_POINTER_STRUCT.amoeba = amoeba_;
  VARTOOLS_FUNCTION_POINTER_STRUCT.vRegisterDataVector = vRegisterDataVector_;
  VARTOOLS_FUNCTION_POINTER_STRUCT.GetOutputFilename = GetOutputFilename_;
  VARTOOLS_FUNCTION_POINTER_STRUCT.incrementparameters_foramoeba = incrementparameters_foramoeba_;
  VARTOOLS_FUNCTION_POINTER_STRUCT.integratemandelagoltransitmodel = integratemandelagoltransitmodel_;
  VARTOOLS_FUNCTION_POINTER_STRUCT.mandelagoltransitmodel = mandelagoltransitmodel_;
  VARTOOLS_FUNCTION_POINTER_STRUCT.spline = spline_;
  VARTOOLS_FUNCTION_POINTER_STRUCT.splint = splint_;
  VARTOOLS_FUNCTION_POINTER_STRUCT.spline_monotonic = spline_monotonic_;
  VARTOOLS_FUNCTION_POINTER_STRUCT.splint_monotonic = splint_monotonic_;
  VARTOOLS_FUNCTION_POINTER_STRUCT.medianfilter = medianfilter_;
}
*/
				   
int VARTOOLS_amoeba(double **p, double *y, int *ia, int ndim, double ftol, double (*funk)(double *, int, int, double *, double *, double *, void *), int *nfunk, int maxeval, int N, double *t, double *mag, double *sig, void *userparams)
{
  int retval;
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.amoeba(p, y, ia, ndim, ftol, funk, nfunk, maxeval, N, t, mag, sig, userparams);
  return(retval);
}

#ifdef PARALLEL
void VARTOOLS_mrqmin(double *x, double *y, double *sig, int ndata, double *a, int *ia, int ma, double **covar, double **alpha, double *chisq, double *alamda, void (*funcs)(double *, double *, double *, double **, int, int, int, double *, double **, double *, double *, int *, void *), int Nlin_coeff, double **Design_Matrix, double *lin_coeffs, int *varylin_coeffs, int mfit, double *ochisq, double *atry, double *beta, double *da, double **oneda, void *userparams)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, chisq, alamda, funcs, Nlin_coeff, Design_Matrix, lin_coeffs, varylin_coeffs, mfit, ochisq, atry, beta, da, oneda, userparams);
}
#else
void VARTOOLS_mrqmin(double *x, double *y, double *sig, int ndata, double *a, int *ia, int ma, double **covar, double **alpha, double *chisq, double *alamda, void (*funcs)(double *, double *, double *, double **, int, int, int, double *, double **, double *, double *, int *, void *), int Nlin_coeff, double **Design_Matrix, double *lin_coeffs, int *varylin_coeffs, void *userparams)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, chisq, alamda, funcs, Nlin_coeff, Design_Matrix, lin_coeffs, varylin_coeffs, userparams);
}
#endif

void VARTOOLS_amoeba_initializesimplexchi2(int Nparameters, int Ntovary, double **p, double **chi2vals, double (*funk)(double *, int, int, double *, double *, double *, void *), int N, double *t, double *mag, double *err, void * userparam)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.amoeba_initializesimplexchi2(Nparameters, Ntovary, p, chi2vals, funk, N, t, mag, err, userparam);
}

void VARTOOLS_incrementparameters_foramoeba(int *Nparameters, int *Ntovary, double ***p, int **ia, int varyparam, double initialval, double stepval)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.incrementparameters_foramoeba(Nparameters, Ntovary, p, ia, varyparam, initialval, stepval);
}

void VARTOOLS_amoeba_cleanup(int *Nparameters, int *Ntovary, double ***p, int **ia, double **chi2vals)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.amoeba_cleanup(Nparameters, Ntovary, p, ia, chi2vals);
}

int VARTOOLS_ParseFixSpecFixcolumn(ProgramData *p,
				   Command *c,
				   int *iret,
				   char **argv,
				   int argc,
				   int Nvec,
				   ...)
{
  va_list varlist;
  int datatype;
  void *dataptr;
  int Ncolumns;
  int output;
  char *name;
  int j, ret;

  _ParseFixSpecFixcolumnStruct *s;

  if((s = (_ParseFixSpecFixcolumnStruct *) malloc(Nvec*sizeof(_ParseFixSpecFixcolumnStruct))) == NULL)
    {
      fprintf(stderr,"Memory Allocation Error\n");
      exit(ERR_MEMALLOC);
    }

  va_start(varlist, Nvec);

  for(j=0; j < Nvec; j++) {
    s[j].datatype = va_arg(varlist,int);
    s[j].dataptr = va_arg(varlist,void *);
    s[j].Ncolumns = va_arg(varlist,int);
    s[j].output = va_arg(varlist,int);
    name = va_arg(varlist,char *);
    if(name != NULL) {
      sprintf(s[j].name,"%s",name);
    } else {
      s[j].name[0] = '\0';
    }
  }

  ret = VARTOOLS_FUNCTION_POINTER_STRUCT.ParseFixSpecFixcolumn(p, c, iret, argv, argc, Nvec, s);
  free(s);
  return(ret);
}

int VARTOOLS_ParseOutNameKeyword(ProgramData *p,
				 Command *c,
				 int *iret,
				 char **argv,
				 int argc,
				 const char *keyword,
				 int *outputflag,
				 char **outdir,
				 int *formatflag,
				 char **format)
{
  int retval;
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.ParseOutNameKeyword(p, c, iret, argv, argc, keyword, outputflag, outdir, formatflag, format);
  return retval;
}

int VARTOOLS_ParseConstantParameter(ProgramData *p,
				    Command *c,
				    int *iret,
				    char **argv,
				    int argc,
				    const char *keyword,
				    int datatype,
				    void *dataptr,
				    int Ncolumns)
{
  int retval;
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.ParseConstantParameter(p, c, iret, argv, argc, keyword, datatype, dataptr, Ncolumns);
  return retval;
}

int VARTOOLS_ParseParameter(ProgramData *p,
			    Command *c,
			    int *iret,
			    char **argv,
			    int argc,
			    const char *keyword,
			    int Nvec,
			    ...)
{
  va_list varlist;
  int datatype;
  void *dataptr;
  int Ncolumns;
  int output;
  char *name;
  int j, ret;

  int Ndblval = 0;
  int Nfltval = 0;
  int Nintval = 0;
  int Nlongval = 0;
  int Nshortval = 0;
  int Ncharval = 0;
  int Nstringval = 0;

  double *dblval = NULL;
  float *fltval = NULL;
  int *intval = NULL;
  long *longval = NULL;
  short *shortval = NULL;
  char *charval = NULL;
  char **stringval = NULL;

  _ParseFixSpecFixcolumnStruct *s;
  _ParseParameter_InitializeStruct *ppstruct;

  if((s = (_ParseFixSpecFixcolumnStruct *) malloc(Nvec*sizeof(_ParseFixSpecFixcolumnStruct))) == NULL)
    {
      fprintf(stderr,"Memory Allocation Error\n");
      exit(ERR_MEMALLOC);
    }
  if((ppstruct = (_ParseParameter_InitializeStruct *) malloc(Nvec*sizeof(_ParseParameter_InitializeStruct))) == NULL)
    {
      fprintf(stderr,"Memory Allocation Error\n");
      exit(ERR_MEMALLOC);
    }
  

  va_start(varlist, Nvec);

  for(j=0; j < Nvec; j++) {
    s[j].datatype = va_arg(varlist,int);
    s[j].dataptr = va_arg(varlist,void *);
    s[j].Ncolumns = va_arg(varlist,int);
    s[j].output = va_arg(varlist,int);
    name = va_arg(varlist,char *);
    if(name != NULL) {
      sprintf(s[j].name,"%s",name);
    } else {
      s[j].name[0] = '\0';
    }
    ppstruct[j].doinitialize = va_arg(varlist,int);
    if(ppstruct[j].doinitialize) {
      switch(s[j].datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	if(!Ndblval) {
	  if((dblval = (double *) malloc(sizeof(double))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	} else {
	  if((dblval = (double *) realloc(dblval, (Ndblval+1)*sizeof(double))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	}
	dblval[Ndblval] = va_arg(varlist,double);
	ppstruct[j].fixptr = (void *) &(dblval[Ndblval]);
	Ndblval++;
	break;
      case VARTOOLS_TYPE_FLOAT:
	if(!Nfltval) {
	  if((fltval = (float *) malloc(sizeof(float))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	} else {
	  if((fltval = (float *) realloc(fltval, (Nfltval+1)*sizeof(float))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	}
	fltval[Nfltval] = (float) va_arg(varlist,double);
	ppstruct[j].fixptr = (void *) &(fltval[Nfltval]);
	Nfltval++;
	break;
      case VARTOOLS_TYPE_INT:
	if(!Nintval) {
	  if((intval = (int *) malloc(sizeof(int))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	} else {
	  if((intval = (int *) realloc(intval, (Nintval+1)*sizeof(int))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	}
	intval[Nintval] = va_arg(varlist,int);
	ppstruct[j].fixptr = (void *) &(intval[Nintval]);
	Nintval++;
	break;
      case VARTOOLS_TYPE_LONG:
	if(!Nlongval) {
	  if((longval = (long *) malloc(sizeof(long))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	} else {
	  if((longval = (long *) realloc(longval, (Nlongval+1)*sizeof(long))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	}
	longval[Nintval] = va_arg(varlist,long);
	ppstruct[j].fixptr = (void *) &(longval[Nlongval]);
	Nlongval++;
	break;
      case VARTOOLS_TYPE_SHORT:
	if(!Nshortval) {
	  if((shortval = (short *) malloc(sizeof(short))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	} else {
	  if((shortval = (short *) realloc(shortval, (Nshortval+1)*sizeof(short))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	}
	shortval[Nintval] = (int) va_arg(varlist,int);
	ppstruct[j].fixptr = (void *) &(shortval[Nshortval]);
	Nshortval++;
	break;
      case VARTOOLS_TYPE_CHAR:
	if(!Ncharval) {
	  if((charval = (char *) malloc(sizeof(char))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	} else {
	  if((charval = (char *) realloc(charval, (Ncharval+1)*sizeof(char))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	}
	charval[Ncharval] = (char) va_arg(varlist,int);
	ppstruct[j].fixptr = (void *) &(charval[Ncharval]);
	Ncharval++;
	break;
      case VARTOOLS_TYPE_STRING:
	if(!Nstringval) {
	  if((stringval = (char **) malloc(sizeof(char *))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	} else {
	  if((stringval = (char **) realloc(stringval, (Nstringval+1)*sizeof(char *))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	}
	stringval[Nstringval] = va_arg(varlist,char *);
	ppstruct[j].fixptr = (void *) &(stringval[Nstringval]);
	Nstringval++;
	break;
      default:
	VARTOOLS_error(ERR_BADTYPE);
      }
    }
  }

  ret = VARTOOLS_FUNCTION_POINTER_STRUCT.ParseParameter(p, c, iret, argv, argc, keyword, Nvec, s, ppstruct);
  free(s);
  free(ppstruct);
  if(dblval != NULL)
    free(dblval);
  if(fltval != NULL)
    free(fltval);
  if(intval != NULL)
    free(intval);
  if(longval != NULL)
    free(longval);
  if(shortval != NULL)
    free(shortval);
  if(charval != NULL)
    free(charval);
  if(stringval != NULL)
    free(stringval);

  return(ret);
}


void VARTOOLS_RegisterDataVector(ProgramData *p, Command *c, void *dataptr,
				int datatype, int Ncolumns, int source,
				int output, char *outname, ...)
{
  va_list varlist;
  va_start(varlist, outname);
  VARTOOLS_FUNCTION_POINTER_STRUCT.vRegisterDataVector(p, c, dataptr,
						      datatype, Ncolumns,
						      source, output, outname,
						       varlist);
  va_end(varlist);
}

void VARTOOLS_GetOutputFilename(char *lcoutname, char *lcname, char *outdir,
				char *suffix, char *format, int lc_name_num)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.GetOutputFilename(lcoutname, lcname, outdir,
						     suffix, format, lc_name_num);
}

void VARTOOLS_integratemandelagoltransitmodel(double exptime_phase, int Npoints,
					      double *phase, double *outlc, 
					      int type, double *ldcoeffs, 
					      double sin_i, double a, double e,
					      double p, double omega, 
					      int Nresamp)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.integratemandelagoltransitmodel(exptime_phase, Npoints, phase, outlc, type, ldcoeffs, sin_i, a, e, p, omega, Nresamp);
}

void VARTOOLS_mandelagoltransitmodel(int Npoints, double *phase, double *outlc,
				     int type, double *ldcoeffs, double sin_i,
				     double a, double e, double p, 
				     double omega)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.mandelagoltransitmodel(Npoints, phase, outlc,
							  type, ldcoeffs, sin_i,
							  a, e, p, omega);
}

void VARTOOLS_spline(double *x,double *y,int n,double yp1,double ypn,
		     double *y2,double *u)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.spline(x, y, n, yp1, ypn, y2, u);
}

void VARTOOLS_splint(double *xa,double *ya,double *y2a,int n,double x,
		     double *y)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.splint(xa, ya, y2a, n, x, y);
}


void VARTOOLS_spline_monotonic(int N, double *x, double *y, double *yprime)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.spline_monotonic(N, x, y, yprime);
}

double VARTOOLS_splint_monotonic(int N, double *x, double *y, double *yprime, 
				 double xt)
{
  double retval;
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.splint_monotonic(N, x, y, yprime, xt);
  return retval;
}

void VARTOOLS_medianfilter(int N, double *t, double *mag, double *sig, 
			   double timesize, int meanflag, int replace)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.medianfilter(N, t, mag, sig, timesize, meanflag, replace);
}

double VARTOOLS_getweightedmean(int n, double *data, double *sig)
{
  double retval;
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.getweightedmean(n, data, sig);
  return retval;
}

double VARTOOLS_getmean(int n, double *data)
{
  double retval;
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.getmean(n, data);
  return retval;
}

double VARTOOLS_median(int n, double *data)
{
  double retval;
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.median(n, data);
  return retval;
}

double VARTOOLS_MAD(int n, double *data)
{
  double retval;
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.MAD(n, data);
  return retval;
}

double VARTOOLS_stddev(int n, double *data)
{
  double retval;
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.stddev(n, data);
  return retval;
}

double VARTOOLS_kurtosis(int n, double *data)
{
  double retval;
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.kurtosis(n, data);
  return retval;
}

double VARTOOLS_skewness(int n, double *data)
{
  double retval;
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.skewness(n, data);
  return retval;
}

double VARTOOLS_percentile(int n, double *data, double pct)
{
  double retval;
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.percentile(n, data, pct);
  return retval;
}

void VARTOOLS_error(int errflag)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.error(errflag);
}

void VARTOOLS_error2(int errflag, char *s)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.error2(errflag, s);
}

double VARTOOLS_fitpoly(int N, double *x, double *y, double *sig, int order,
			int subtractfit, double *fitparams, double *paramerrs)
{
  double retval;
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.fitpoly(N, x, y, sig, order, subtractfit, fitparams, paramerrs);
  return retval;
}

double VARTOOLS_chi2(int N, double *t, double *mag, double *err, 
		     double *weighted_average, int *Ngood)
{
  double retval;
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.chi2(N, t, mag, err, 
						 weighted_average, Ngood);
  return retval;
}

int VARTOOLS_isDifferentPeriods(double period1, double period2, double TimeSpan)
{
  int retval;
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.isDifferentPeriods(period1, period2, TimeSpan);
}

void VARTOOLS_sort_generic(int N, int isreverse, int *index, int Nms, ...)
/* This function can be used to sort one or more data vectors of various types.
   The syntax is as follows:
      N - number of datapoints in each of the vectors.
      isreverse - 0 = sort from low to high.
                  1 = sort from high to low.
      index - = A array of integer indices to be provided if we don't want to
                transpose the data for at least one of the vectors. Set this
                to NULL if all of the vectors will be transposed. If this is
                not null, it should initially store the values 0 to N-1 in
                order. The values of index will then be transposed such that
                dataN[index[i=0 ... N-1]] will be sorted based on the key.
      Nms - The number of data vectors. This must be >= 1 or else VARTOOLS
            will abort with an error.
      ...
          - For each data vector provide the additional input arguments
            given below. Note that the first vector described will be the
            key used for the sorting. Additional vectors will be re-ordered
            with the key.
            
            int datatype = the datatype for this vector. Allowed values are
                           VARTOOLS_TYPE_DOUBLE, VARTOOLS_TYPE_STRING, 
                           VARTOOLS_TYPE_INT, VARTOOLS_TYPE_FLOAT,
                           VARTOOLS_TYPE_LONG, VARTOOLS_TYPE_CHAR,
                           VARTOOLS_TYPE_SHORT, VARTOOLS_TYPE_USERDEF (for
                               an array of user-defined structures).
            void *dataptr = the data vector (e.g. if the vector is
                           s of type (double *), one would pass (void *)s
                           for this argument).
            int useindex = 0 - the data in the vector will be transposed.
                           1 - the data will not be transposed and instead
                               will be indexed. Note that if this is one,
                               then index above cannot be NULL.

          - If datatype is VARTOOLS_TYPE_STRING the following additional 
            argument is needed.

            size_t sizeobj = the size of single string element.

          - If datatype is VARTOOLS_TYPE_USERDEF the following additional
            two arguments are needed.

            size_t sizeobj = the size of a single element.

            int (*comp_func)((void *),(void *)) = a pointer to a function
                used to compare two elements. This function is of the form
                used by the qsort function (e.g. strcmp for strings).
                If this vector is not the key then this can be NULL.
*/                          
{
  mysort_generic_struct *s;
  int i;
  va_list varlist;
  int datatype;
  void *dataptr;
  comp_func_type comp_func;
  int useindex;
  size_t sizeobj;

  va_start(varlist, Nms);
  if(Nms <= 0) {
    VARTOOLS_error(ERR_MYSORT_GENERIC_BADCALL);
  }
  if((s = (mysort_generic_struct *) malloc(Nms * sizeof(mysort_generic_struct))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  for(i=0; i < Nms; i++) {
    datatype = va_arg(varlist,int);
    dataptr = va_arg(varlist,void *);
    useindex = va_arg(varlist,int);
    if(useindex && index == NULL)
      VARTOOLS_error(ERR_MYSORT_GENERIC_BADCALL);
    if(datatype == VARTOOLS_TYPE_STRING) {
      sizeobj = va_arg(varlist,size_t);
      comp_func = NULL;
    }
    else if(datatype == VARTOOLS_TYPE_USERDEF) {
      sizeobj = va_arg(varlist,size_t);
      comp_func = va_arg(varlist,comp_func_type);
    }
    else {
      sizeobj = 0;
      comp_func = NULL;
    }
    s[i].datatype = datatype;
    s[i].dataptr = dataptr;
    s[i].comp_func = comp_func;
    s[i].useindex = useindex;
    s[i].sizeobj = sizeobj;
  }
  va_end(varlist);
  
  VARTOOLS_FUNCTION_POINTER_STRUCT.vsort_generic(N, isreverse, index, Nms, s);

  free(s);
}

void VARTOOLS_sortvec_double(int N, double* data1)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.sortvec_double(N, data1);
}

int VARTOOLS_RegisterUserFunction(ProgramData *p, char *funcname, int Nexpr,
				  double (*func)(double *), int ishelp, ...)
{
  va_list varlist;
  int retval;
  va_start(varlist, ishelp);
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.vRegisterUserFunction(p,
								 funcname,
								  Nexpr, func,
								  ishelp,
								  varlist);
  va_end(varlist);
  return retval;
}

void VARTOOLS_occultquad(double *z0, double u1, double u2, double p, double *muo1, double *mu0, int nz)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.occultquad(z0, u1, u2, p, muo1, mu0, nz);
}

void VARTOOLS_occultnl(double rl, double c1, double c2, double c3, double c4, double *b0, double *mulimb0, double **mulimbf, int nb)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.occultnl(rl, c1, c2, c3, c4, b0, mulimb0, mulimbf, nb);
}

void VARTOOLS_MemAllocDataFromLightCurve(ProgramData *p, int threadid, int Nterm)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.memallocdatafromlightcurve(p, threadid, Nterm);
}

void VARTOOLS_MemAllocDataFromLightCurveMidProcess(ProgramData *p, int threadid, int Nterm)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.memallocdatafromlightcurvemidprocess(p, threadid, Nterm);
}

int VARTOOLS_gnu_getline(char **str, size_t *sizeval, FILE *infile)
{
  int retval;
  retval = VARTOOLS_FUNCTION_POINTER_STRUCT.gnu_getline(str, sizeval, infile);
  return retval;
}

void VARTOOLS_mysortstringint(int N, int sizestr, char **data1, int *data2)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.mysortstringint(N, sizestr, data1, data2);
}

void VARTOOLS_docorr(double *mag, double *err, int Npoints, int ndecorr, double **decorr, int *order, double *Avector, double *A_errvector, double mag_ave, int zeropoint)
{
  VARTOOLS_FUNCTION_POINTER_STRUCT.docorr(mag, err, Npoints, ndecorr, decorr, order, Avector, A_errvector, mag_ave, zeropoint);
}

