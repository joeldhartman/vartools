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
#ifndef _OUTTEXTSTRUCTDEFINE
#include "OutText.h"
#endif 
#include <stdio.h>
#include <stdlib.h>

int VARTOOLS_ParseParameter(ProgramData *p,
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
				    char datatype,
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

void VARTOOLS_amoeba_initializesimplexchi2(int Nparameters, int Ntovary, double **p, double **chi2vals, double (*funk)(double *, int, int, double *, double *, double *, void *), int N, double *t, double *mag, double *err, void * userparam);

void VARTOOLS_amoeba_cleanup(int *Nparameters, int *Ntovary, double ***p, int **ia, double **chi2vals);

void VARTOOLS_RegisterDataVector(ProgramData *p, Command *c, void *dataptr,
				 char datatype, int Ncolumns, int source,
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

int VARTOOLS_RegisterUserFunction(ProgramData *, char *, int, double (*)(double *), int, ...);

void VARTOOLS_occultquad(double *z0, double u1, double u2, double p, double *muo1, double *mu0, int nz);

void VARTOOLS_occultnl(double rl, double c1, double c2, double c3, double c4, double *b0, double *mulimb0, double **mulimbf, int nb);

void VARTOOLS_MemAllocDataFromLightCurve(ProgramData *p, int threadid, int Nterm);

void VARTOOLS_MemAllocDataFromLightCurveMidProcess(ProgramData *p, int threadid, int Nterm);

int VARTOOLS_gnu_getline(char **, size_t *, FILE *);

void VARTOOLS_mysortstringint(int, int, char **, int *);

void VARTOOLS_docorr(double *mag, double *err, int Npoints, int ndecorr, double **decorr, int *order, double *Avector, double *A_errvector, double mag_ave, int zeropoint, int usemask, _Variable *maskvar, int lcindex, int threadindex);

void VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(ProgramData *p, int lcnum,
						 char *keyname,
						 char *comment, int hdutouse,
						 int updateexisting,
						 int dtype, ...);
void VARTOOLS_Delete_Keyword_From_OutputLC_FitsHeader(ProgramData *p, int lcnum,
						      char *keyname, int hdutouse,
						      int prefixmatch);
int VARTOOLS_findX(double *, double, int, int);

int VARTOOLS_findX_string(char **, int *, char *, int, int);

void VARTOOLS_RegisterTrackedOpenFile(ProgramData *p, FILE *f);

int VARTOOLS_parseone(char *line, void *val, int vartype);

void VARTOOLS_printtostring(OutText *text, const char *stoadd);

/* ---------------------------------------------------------------------------
 * Pipeline library API (Step 2: library extraction)
 *
 * These three functions implement an init-once / process-many interface for
 * embedding vartools as an in-process library.
 *
 * Usage:
 *   ProgramData *p = vartools_init_pipeline(argc, argv);
 *   vartools_process_lc(p, t, mag, err, n, "lcname", outbuf, outbuf_size);
 *   vartools_free_pipeline(p);
 * ---------------------------------------------------------------------------
 */

/* Initialise a pipeline from a vartools command-line (without -i / -l).
 * Internally inserts "-i -" to configure for single-LC stdin mode.
 * Returns a heap-allocated opaque handle on success, NULL on error. */
ProgramData *vartools_init_pipeline(int argc, char **argv);

/* Inject one light curve and run the pipeline.
 * t, mag, err  – arrays of length n
 * lc_name      – name written into the output (may be NULL, defaults to "lc")
 * outbuf        – caller-provided buffer filled with -oneline output text
 * outbuf_size   – size of outbuf in bytes (e.g. 65536)
 * Returns 0 on success, non-zero error code on failure. */
int vartools_process_lc(ProgramData *p,
                        const double *t, const double *mag, const double *err,
                        int n, const char *lc_name,
                        char *outbuf, int outbuf_size);

/* Release all resources held by the pipeline handle. */
void vartools_free_pipeline(ProgramData *p);

/* Information about a single variable in the pipeline. */
typedef struct {
    const char *name;       /* variable name (e.g. "t", "mag", "err", "myvar") */
    int         datatype;   /* VARTOOLS_TYPE_DOUBLE, _INT, etc. */
    int         vectortype; /* VARTOOLS_VECTORTYPE_LC, _SCALAR, _PERSTARDATA */
    const void *dataptr;    /* pointer to raw data array for thread 0 */
    int         length;     /* NJD for LC vars, 1 for scalars */
} VartoolsVarInfo;

/* Enumerate all light-curve and per-star variables after processing.
 * Fills vars[0..min(count, max_vars)-1].  *n_vars is set to the number
 * of entries written.  Returns 0 if all fit, or the total count if
 * max_vars was too small.  Data pointers are valid until the next
 * vartools_process_lc() or vartools_free_pipeline(). */
int vartools_get_lc_variables(ProgramData *p,
                              int *n_vars, VartoolsVarInfo *vars,
                              int max_vars);

/* Return the number of LC points after processing (may differ from input). */
int vartools_get_njd(ProgramData *p);

/* Snapshot the current LC variables into the p->captured[] slot keyed
 * by *id*.  Used internally by DoOutputLightCurve when an
 * "-o <id> capture" command is reached.  Returns 0 on success, non-zero
 * if no slot matches *id* (which should not happen — slots are pre-
 * populated at vartools_init_pipeline() from the parsed argv). */
int vartools_capture_current_lc(ProgramData *p, const char *id);

/* Read a captured LC snapshot back out by *id*.  Same shape as
 * vartools_get_lc_variables: fills vars[0..min(count, max_vars)-1] and
 * sets *n_vars.  Returns 0 on success, the total count if max_vars was
 * too small, -1 if id not found, or -2 if the slot was never filled.
 * Data pointers are valid until the next vartools_process_lc() call. */
int vartools_get_captured_lc(ProgramData *p, const char *id,
                             VartoolsVarInfo *vars, int max_vars,
                             int *n_vars);

/* Return the NJD of a captured snapshot, or -1 (id not found) /
 * -2 (slot never filled). */
int vartools_get_captured_njd(ProgramData *p, const char *id);

/* Inject a full set of LC columns before processing.
 * col_names[0..2] must be "t", "mag", "err".  Additional columns are
 * matched by name against DefinedVariables (VECTORTYPE_LC, TYPE_DOUBLE).
 * Returns 0 on success, -1 if n_columns < 3. */
int vartools_set_lc_data(ProgramData *p,
                         int n_points, int n_columns,
                         const char **col_names,
                         const double **col_data,
                         const char *lc_name);

/* ------------------------------------------------------------------------ */
/* In-process -python callback hook (libvartoolspipeline only).             */
/*                                                                          */
/* When a callback is registered via vartools_register_python_callback(),   */
/* RunPythonCommand bypasses the per-thread Python sub-process / Unix-      */
/* socket path and instead calls the callback directly.  The callback       */
/* receives flat C arrays describing the named in-vars and out-vars.  This  */
/* lets pyvartools (running in library mode) execute -python user code in   */
/* the host Python interpreter, sharing a caller-supplied namespace dict.   */
/*                                                                          */
/* The interface is libpython-free at this boundary: namespace_ptr is an    */
/* opaque void* that vartools never dereferences; types and lengths are     */
/* primitives.  pyvartools' shim casts namespace_ptr back to PyObject* on   */
/* its side using its own libpython.                                        */
/*                                                                          */
/* Variable types use the existing VARTOOLS_TYPE_* and VARTOOLS_VECTORTYPE_*/
/* constants.  v1 supports DOUBLE/FLOAT/INT/LONG (scalar + LC vectors);     */
/* STRING and process_all_lcs=1 are NOT yet wired through this hook —       */
/* fall back to the subprocess path for those.                              */
/*                                                                          */
/* Returns 0 on success, 1 on user exception (error_buf populated).         */
typedef int (*vartools_python_callback_t)(
        void       *namespace_ptr,        /* opaque PyObject* */
        int         command_id,           /* unique per -python command */
        const char *code,                 /* code body to exec */
        int         is_init,              /* 1 = init code, 0 = per-LC body */
        int         n_invars,
        const char *const *invar_names,
        const int  *invar_types,          /* VARTOOLS_TYPE_* */
        const int  *invar_lengths,        /* 1 for scalar, NJD for vector */
        void *const *invar_data,
        int         n_outvars,
        const char *const *outvar_names,
        const int  *outvar_types,
        const int  *outvar_lengths,       /* expected length (1 for scalar) */
        void *const *outvar_data,         /* pre-allocated buffers */
        char       *error_buf,
        int         error_buf_size);

/* Register / clear the in-process callback.  Pass NULL to unregister and
 * return to the default subprocess behavior. */
void vartools_register_python_callback(vartools_python_callback_t cb);

/* Set the namespace pointer passed to the callback.  Stored as an opaque
 * void* and passed back unchanged in `namespace_ptr`. */
void vartools_set_python_namespace(void *namespace_ptr);
