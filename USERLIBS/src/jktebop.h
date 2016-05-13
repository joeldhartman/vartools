/* This header defines the _Jktebop data structure
   which is used to pass parameters from vartools
   to the jktebop library routines.
   It makes other constant definitions used by the library.
 */

#define JKTEBOP_INJECT 0
#define JKTEBOP_FIT 1

#define JKTEBOP_USEI 0
#define JKTEBOP_USEBIMPACT 1

#define JKTEBOP_LDLAW_LINEAR 1
#define JKTEBOP_LDLAW_QUAD 2
#define JKTEBOP_LDLAW_LOG 3
#define JKTEBOP_LDLAW_SQRT 4
#define JKTEBOP_LDLAW_LOCKLD1 5

#define JKTEBOP_DEFAULT_GRAVDARK 1.0 /* Default gravity darkening coefficient
                                        appropriate for a star with
                                        M > 1.0 M_sun */
#define JKTEBOP_DEFAULT_REFLECTION -1.0 /* A reflection coefficient <=
					   0 indicates that we want to
					   have the reflection
					   coefficients calculated */
#define JKTEBOP_DEFAULT_L3 0.0 /* By default assume no third light */
#define JKTEBOP_DEFAULT_TIDALLAG 0.0 /* By default assume no tidal lag angle */

#define JKTEBOP_OUTCURVE_TYPE_JD 0
#define JKTEBOP_OUTCURVE_TYPE_PHASE 1
#define JKTEBOP_OUTCURVE_STEP_DEFAULT 0.01

typedef struct {
  int injectorfit; /* flag indicating whether to inject an EB model into
                      the light curve, or to fit a model to the light curve */

  /* Variables storing parameters such as Period, T0, r1r2, etc are
     defined as vectors. VARTOOLS will allocate memory for these and
     initialize them based on the command-line issued by the user. For
     each parameter we also have an integer flag indicating whether or
     not the parameter is to be varied in a fit. */
  double *Period;
  int Period_vary;
  double *T0;
  int T0_vary;
  double *r1r2;
  int r1r2_vary;
  double *r2_r1;
  int r2_r1_vary;
  double *M2_M1;
  int M2_M1_vary;
  double *J2_J1;
  int J2_J1_vary;
  /* The user can either specify the inclination angle or the impact
     parameter of the primary eclipse */
  int use_i_or_b;
  double *incl;
  int incl_vary;
  double *bimpact;
  int bimpact_vary;
  double *esinomega;
  int esinomega_vary;
  double *ecosomega;
  int ecosomega_vary;

  /* These integers indicate the limb-darkening law to be used */
  int LD1law;
  int LD2law;
  
  int LD1_vary;
  int LD2_vary;
  /* We use double-pointers to store the LD coefficients.
     LD1_coeffs[lcindx][0] will be the first LD coefficient for star
     one for the light curve lcindx, and LD1_coeffs[lcindx][1] will be
     the second LD coefficient, etc . */
  double **LD1_coeffs;
  double **LD2_coeffs;
  
  double *gravdark1;
  int gravdark1_vary;
  double *gravdark2;
  int gravdark2_vary;

  double *reflection1;
  int reflection1_vary;
  double *reflection2;
  int reflection2_vary;

  double *L3;
  int L3_vary;

  double *tidallag;
  int tidallag_vary;

  /* Flag indicating whether or not to subtract the best-fit model
     from the light curve */
  int correctlc;

  /* Flag indicating whether or not to output the best-fit model
     light curve */
  int omodel;

  /* Variables storing the output directory name, and the format (if given) */
  char outdir[MAXLEN];
  char format[MAXLEN];

  /* Variables for outputting a model curve */
  int ocurve;
  int ocurvetype;
  double ocurvestep;
  char ocurve_outdir[MAXLEN];
  char ocurve_outdir_format[MAXLEN];

  /* Variables storing the output chi2 and Ndof */
  double *chi2val;
  int *Ndof;

} _Jktebop;

/* The structure below is used to pass auxiliary parameters to chisqjktebop
   when fitting a JKTEBOP model. */
typedef struct {
  int N1;  /* Initial eclipse number */
  int N2;  /* Final eclipse number */
  int use_i_or_b; /* Whether the inclination or impact param is fitted */
  double ival; /* Variables to store the calculated value of whichever param
		  is not fitted */
  double bval;
  int is_vary_P; /* Whether the period is allowed to vary */
  int is_vary_T0; /* Whether the initial eclipse time is allowed to vary */
  
  int LD1law; /* Limb darkening laws to use for each star */
  int LD2law;

  double *magsim; /* Used to store the simulated lc at the model
		   output values */
} _Jktebop_fitstruct;
