#define FTUNEVEN_OUTPUTVECTORS 0
#define FTUNEVEN_OUTPUTFILE 1
#define FTUNEVEN_OUTPUTVECTORSANDFILE 2

#define FTUNEVEN_FREQTYPE_AUTO 0
#define FTUNEVEN_FREQTYPE_RANGE 1
#define FTUNEVEN_FREQTYPE_VARIABLE 2
#define FTUNEVEN_FREQTYPE_FILE 3

typedef struct {
  int outputtype;
  char *frequency_vecname;
  double **frequency_vec;
  char *FT_real_vecname;
  double **FT_real_vec;
  char *FT_imag_vecname;
  double **FT_imag_vec;
  char *Periodogram_vecname;
  double **Periodogram_vec;
  char *FT_outdir;
  int FT_nameformatflag;
  char *FT_nameformat;
  int freqtype;
  double *minfreq;
  double *maxfreq;
  double *freqstep;
  char *freqlist_varname;
  char *freqlist_filename;
  double **input_freqvar_vals;
  double ft_sign;
  double tt_zero;
  int Nfreq_file;
  int Nfreq_vecsize;
  double *freqvals_file;
  int changeinputvectors;
  char *input_tvar_varname;
  double **input_tvar_vals;
  char *input_real_varname;
  double **input_real_vals;
  char *input_imag_varname;
  double **input_imag_vals;
} _FTuneven;
