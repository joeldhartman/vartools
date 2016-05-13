/* This header defines the _Macula data structure
   which is used to pass parameters from vartools
   to the macula routines.
   It makes other constant definitions used by the library.
 */

/*
#ifndef MACULAFUNC
#define MACULAFUNC __maculamod_MOD_macula
#endif*/

void MACULAFUNC(double JD[], int*N, int *Nspot, int *Nsets, 
		int *derivatives, int *temporal, int *TdeltaV, 
		double Theta_star[], double Theta_spot_fortran[], 
		double Theta_inst_fortran[], double *tstart, 
		double *tend, double fluxout[], 
		double dFmod_star_fortran[], 
		double dFmod_spot_fortran[],
		double dFmod_inst_fortran[], 
		double dFmoddt[], double deltaratio[]);

#define MACULA_INJECT 0
#define MACULA_FIT 1

#define MACULA_FITTYPE_AMOEBA 0
#define MACULA_FITTYPE_LM 1

#define MACULA_OUTCURVE_STEP_DEFAULT 0.01

#define MACULA_NDIM_STAR_PARAM 12
#define MACULA_NDIM_SPOT_PARAM 8
#define MACULA_NDIM_INST_PARAM 2

typedef struct {
  int injectorfit; /* flag indicating whether to inject a spot model into
                      the light curve, or to fit a model to the light curve */

  int fittype; /* flag indicating what type of fitting to do */

  /* Variables storing parameters such as Period, T0, r1r2, etc are
     defined as vectors. VARTOOLS will allocate memory for these and
     initialize them based on the command-line issued by the user. For
     each parameter we also have an integer flag indicating whether or
     not the parameter is to be varied in a fit. */
  
  double *Istar;
  int Istar_vary;

  double *Prot;
  int Prot_vary;

  double *kappa2;
  int kappa2_vary;

  double *kappa4;
  int kappa4_vary;

  double *c1;
  int c1_vary;

  double *c2;
  int c2_vary;

  double *c3;
  int c3_vary;

  double *c4;
  int c4_vary;

  double *d1;
  int d1_vary;

  double *d2;
  int d2_vary;

  double *d3;
  int d3_vary;

  double *d4;
  int d4_vary;

  double *blend;
  int blend_vary;

  int Nspot;

  double **lambda0;
  int *lambda0_vary;

  double **phi0;
  int *phi0_vary;
  
  double **alphamax;
  int *alphamax_vary;

  double **fspot;
  int *fspot_vary;

  double **tmax;
  int *tmax_vary;

  double **life;
  int *life_vary;

  double **ingress;
  int *ingress_vary;

  double **egress;
  int *egress_vary;

  int outputmodel;
  char *modeloutdir;
  int outputmodel_useformat;
  char *outputmodel_format;
  int outputmodel_tdelv;

  int outputcurve;
  char *curveoutdir;
  int outputcurve_useformat;
  char *outputcurve_format;
  int outputcurve_tdelv;
  int outputcurve_stepgiven;
  double outputcurve_step;

  double *chi2val;
  int *Ndof;
  
  int correctlc;
  int fluxinput;
  int fluxoutput;

} _Macula;

/* Structure used to pass auxiliary parameters to macula_evalfunc_lm */
typedef struct {
  double *Theta_star;
  double **Theta_spot;
  double *dFmoddt;
  double *deltaratio;
  double **dFmod_star;
  double ***dFmod_spot;
  double ***dFmod_inst;
  double blend;
  double *Fmod;

  int Istar_vary;
  int Prot_vary;
  int kappa2_vary;
  int kappa4_vary;
  int c1_vary;
  int c2_vary;
  int c3_vary;
  int c4_vary;
  int d1_vary;
  int d2_vary;
  int d3_vary;
  int d4_vary;
  int blend_vary;

  int Nspot;

  int *lambda0_vary;
  int *phi0_vary;
  int *alphamax_vary;
  int *fspot_vary;
  int *tmax_vary;
  int *life_vary;
  int *ingress_vary;
  int *egress_vary;
  
} _MaculaLMFitStruct;

/* The structure below is used to pass auxiliary parameters to chisqmacula
   when fitting a macula model. */
typedef struct {

  double *magsim; /* Used to store the simulated lc at the model
		   output values */
} _Macula_fitstruct;
