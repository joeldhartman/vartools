#define VARTOOLS_STITCH_METHOD_MEDIAN 0
#define VARTOOLS_STITCH_METHOD_MEAN 1
#define VARTOOLS_STITCH_METHOD_WEIGHTEDMEAN 2
#define VARTOOLS_STITCH_METHOD_POLY 3
#define VARTOOLS_STITCH_METHOD_HARM 4

typedef struct {
  int Nlcs;
  int *lcids;
  int lcnumval;
  int *timegroups;
  int Ntimegroups_uniqlist;
  int *timegroups_uniqlist;
  double shiftvalue;
} _StitchLightCurveGroup;

typedef struct {
  int Nlcs;
  int *lcids;
  int tbin;
  int *lcgids;
  int *lcglcindx;
  int Nlcgroups_uniqlist;
  int *lcgid_uniqlist;
  double mint;
  double maxt;
} _StitchTimeGroup;

typedef struct {
  int nstitchvar;
  int stitchmethod;
  int polyorder;
  int Nharm;
  double *harmperiodvals;
  _Variable *harmperiodvar;
  int groupbytime;
  int is_start_time;
  double start_time;
  double time_step;
  int fitonly;
  int save_fitted_parameters;
  char coeffoutdir[MAXLEN];
  char coeffoutformat[MAXLEN];
  int add_shifts_fitsheader;
  int add_shifts_fitsheader_hdutouse;
  int add_shifts_fitsheader_updateexistingkeyword;
  char keywordbase[MAXLEN];
  int add_stitchparams_fitsheader;
  int add_stitchparams_fitsheader_hdutouse;
  int add_stitchparams_fitsheader_updateexistingkeyword;  
  int **lcnumval;
  _Variable *lcnumvar;
  int userefnum;
  int **refnumval;
  _Variable *refnumvar;
  int *minrefnumindx;
  int *maxrefnumindx;

  char **stitchvarnames;
  double ***stitchvarvals;
  char **stitcherrnames;
  double ***stitcherrvals;
  char **stitchmasknames;
  double ***stitchmaskvals;

  _Variable **stitchvars;
  _Variable **stitcherrvars;
  _Variable **stitchmaskvars;

  int *Nlcgroups_used;
  int *Ntimegroups_used;
  int *Nparamtotal_used;
  double *final_time_step;

  int is_in_shifts_file;
  char **in_shifts_filename;

  int *N_in_shift_stars;
  int *size_in_shifts_file;

  int **N_shifts_per_star;
  char ***field_labels_vals;
  int **field_labels_vals_indx;
  int *field_labels_vals_issorted;
  _Variable *field_labels_var;
  int is_append_refnum_to_fieldlabel;
  char **starname_vals;
  _Variable *starname_var;
  char *starname_varname;
  char ***in_shift_starnames;
  int **in_shift_stars_found;
  int **in_shift_starnames_sortidx;
  char ****in_shift_labels;
  int ***in_shift_labels_sortidx;
  double ***in_shift_values;
  char ****in_shift_values_str;
  int ***Nobs_in_shifts;
  
  int is_nobs_refit;
  int nobs_refit;
  int is_inshifts_header_basename_only;

  int is_out_shifts_file;
  char **out_shifts_filename;
  FILE **out_shifts_file;
  int **N_shifts_per_star_out;
  char ****out_shift_labels;
  double ***out_shift_values;
  int ***Nobs_out_shifts;
  int include_missing_inputstars;
  int has_output_missing;

  
} _Stitch;
