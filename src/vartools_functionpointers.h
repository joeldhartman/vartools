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
typedef struct {
  int (*ParseFixSpecFixcolumn)(ProgramData *, Command *, int *,
			       char **, int, int, _ParseFixSpecFixcolumnStruct *);
  int (*ParseOutNameKeyword)(ProgramData *, Command *, int *, char **, int,
			     const char *, int *, char **, int *, char **);
  int (*ParseConstantParameter)(ProgramData *, Command *, int *,
				char **, int, const char *, char, void *,
				int);
  int (*ParseParameter)(ProgramData *, Command *, int *,
			char **, int, const char *,
			int, _ParseFixSpecFixcolumnStruct *,
			_ParseParameter_InitializeStruct *);
  int (*amoeba)(double **, double *, int *, int, double, double (*)(double *, int, int, double *, double *, double *, void *), int *, int, int, double *, double *, double *, void *);
#ifdef PARALLEL
  void (*mrqmin)(double *, double *, double *, int, double *, int *, int, double **, double **, double *, double *, void (*)(double *, double *, double *, double **, int, int, int, double *, double **, double *, double *, int *, void *), int, double **, double *, int *, int, double *, double *, double *, double *, double **, void *);
#else
  void (*mrqmin)(double *, double *, double *, int, double *, int *, int, double **, double **, double *, double *, void (*)(double *, double *, double *, double **, int, int, int, double *, double **, double *, double *, int *, void *), int, double **, double *, int *, void *);
#endif
  void (*vRegisterDataVector)(ProgramData *, Command *, void *, char, int, int,
			     int, char *, va_list argp);
  void (*GetOutputFilename)(char *, char *, char *, char *, char *, int);
  void (*incrementparameters_foramoeba)(int *, int *, double ***, int **, int, double, double);
  void (*amoeba_initializesimplexchi2)(int, int, double **, double **, double (*)(double *, int, int, double *, double *, double *, void *), int, double *, double *, double *, void *);
  void (*amoeba_cleanup)(int *, int *, double ***, int**, double **);
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
  double (*percentile)(int, double *, double);
  void (*error)(int);
  void (*error2)(int, char *);
  double (*fitpoly)(int, double *, double *, double *, int, int, double *, double *);
  double (*chi2)(int, double *, double *, double *, double *, int *);
  int (*isDifferentPeriods)(double, double, double);

  /* Sort functions */
  void (*vsort_generic)(int, int, int *, int, mysort_generic_struct *);
  void (*sortvec_double)(int, double *);
  int (*vRegisterUserFunction)(ProgramData *, char *, int, double (*)(double *), int, va_list argp);
  void (*occultquad)(double *, double, double, double, double *, double *, int);
  void (*occultnl)(double, double, double, double, double, double *, double *, double **, int);
  void (*memallocdatafromlightcurve)(ProgramData *, int, int);
  void (*memallocdatafromlightcurvemidprocess)(ProgramData *, int, int);
  int (*gnu_getline)(char **, size_t *, FILE *);
  void (*mysortstringint)(int, int, char **, int*);
  void (*docorr)(double *, double *, int, int, double **, int *, double *, double *, double, int, int, _Variable *, int, int);
  void (*vAdd_Keyword_To_OutputLC_FitsHeader)(ProgramData *, int, char *,
					     char *, int, int,
					     int, va_list argp);
  int (*findX)(double *, double, int, int);
  int (*findX_string)(char **, int *, char *, int, int);
  void (*RegisterTrackedOpenFile)(ProgramData *, FILE *);
  int (*parseone)(char *, void *, int);
  void (*printtostring)(OutText *, const char*);
} _VARTOOLS_FUNCTION_POINTER_STRUCT;
