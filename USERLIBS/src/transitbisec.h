#define TRANSITBISEC_TSAMP_INTEGRATE 0.00034722222
#define TRANSITBISEC_DEFAULT_BISEC_HIGH_VAL 0.95
#define TRANSITBISEC_DEFAULT_BISEC_LOW_VAL 0.05
#define TRANSITBISEC_DEFAULT_BISEC_STEP 0.005
#define TRANSITBISEC_DEFAULT_SPAN_HIGH_VAL 0.9
#define TRANSITBISEC_DEFAULT_SPAN_LOW_VAL 0.1
#define TRANSITBISEC_DEFAULT_MAX_DELTA_T 0.03
#define TRANSITBISEC_DEFAULT_MIN_OOT_FRAC 0.25
#define TRANSITBISEC_DEFAULT_NERRORSIMS 100
#define TRANSITBISEC_DEFAULT_RESID_BIN_TIME 0.05
#define TRANSITBISEC_DEFAULT_RESID_FIT_POLY_ORDER 1
#define TRANSITBISEC_AMOEBA_CONVERGENCELIMIT 0.0001
#define TRANSITBISEC_DEFAULT_BISEC_HIGH_SAMP_FACTOR 10

typedef struct {
  int Nbisecs;
  double **tvals;
  int Ntvals;
  size_t sizebisec1;
  size_t sizebisec2;
  double bisec_high_val;
  double bisec_low_val;
  double bisec_step_val;
  double *spanvals;
  int *transitnums;
  double bisec_resid_bin_time;
  int bisec_resid_fit_poly_order;
  int bisec_usemoments;
  int bisec_high_samp_factor;
} Bisector_Data;

typedef struct {
  double Tc0;
  double P;
  double rprstar;
  double bimpact;
  double arstar;
  double eccen;
  double omega;
  double exptime;
  double *ldcoeffs;
  int ldtype;
} TransitParams;

typedef struct {
  TransitParams *transitparams;
  Bisector_Data *bisecdata;
  double BS_highval;
  double BS_lowval;
  char bisecoutdir[MAXLEN];
  char bisecspanoutdir[MAXLEN];
  int outputbisec;
  int outputbisecspan;
  int outputbisec_useformat;
  char outputbisec_format[MAXLEN];
  int outputspan_useformat;
  char outputspan_format[MAXLEN];
  int geterror;
  int Nerrorsims;
  double *out_BSrms;
  double *out_BSmed;
  int *out_Ntransits;
  double *out_BSrms_sim;
  double *out_BSmed_sim;
  double bisec_high_val;
  double bisec_low_val;
  double bisec_step_val;
  double bisec_resid_bin_time;
  int bisec_resid_fit_poly_order;
  double max_delta_t;
  double min_oot_frac;

  int bisec_high_samp_factor;

  double *Tc0;
  double *P;
  double *rprstar;
  double *bimpact;
  double *arstar;
  double *eccen;
  double *omega;
  double *exptime;
  double **ldcoeffs;
  int Nld;
  int ldtype;
  int refit_Tc;
  int bisec_usemoments;
} _TransitBisec;

typedef struct {
  TransitParams *t;
  Bisector_Data *b;
} _TransitBisec_passparams;
