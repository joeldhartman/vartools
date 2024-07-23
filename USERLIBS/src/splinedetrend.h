/* This header defines the _Splinedetrend data structure
   which is used to pass parameters from vartools
   to the splinedetrend library routines.
   It makes other constant definitions used by the library.
 */

#define SPLINEDETREND_DEFAULT_SIGMACLIP 5.0
#define SPLINEDETREND_DEFAULT_GROUPBYGAP -1.0
#define SPLINEDETREND_DEFAULT_SGRIDLOGMIN -4.0
#define SPLINEDETREND_DEFAULT_SGRIDLOGMAX -1.0
#define SPLINEDETREND_DEFAULT_SGRIDNUM 50
#define SPLINEDETREND_DEFAULT_SPLINEORDER 3

#define SPLINEDETREND_VARIABLE_FITTYPE_SPLINE 0
#define SPLINEDETREND_VARIABLE_FITTYPE_POLY 1
#define SPLINEDETREND_VARIABLE_FITTYPE_HARM 2

typedef struct {
  char *varname;
  int fittype;
  int order;
  int nharm;
  int groupbygap;
  double knotspacing;
  double gapsize;
  double **inputdata;
} _SplinedetrendVariable;

typedef struct {
  char *varname;
  _SplinedetrendVariable *invar;
  double **outputdata;
} _SplinedetrendOutputVariable;


typedef struct {
  int Nvar;
  _SplinedetrendVariable *detrendvars;
  double *sigmaclip;
  int outputmodel;
  char *modeloutdir;
  int outputmodel_useformat;
  char *outputmodel_format;
  int outputmodelcoeffs;
  char *modelcoeffsoutdir;
  int outputmodelcoeffs_useformat;
  char *outputmodelcoeffs_format;
  int Noutvar;
  _SplinedetrendOutputVariable *outvars;
  double *magmedian;
  int *Noutlier;
  int *Ngroups;
  int *Nparamtotal;
} _Splinedetrend;
