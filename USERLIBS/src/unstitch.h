/* Data structure for the -unstitch user command.

   -unstitch is the inverse of -stitch: it adds previously-determined
   per-segment shifts back to a light curve to restore the original
   magnitudes.  Phase A implements the "in_shifts_file" source, which reads
   the shifts from a file written by -stitch's "out_shifts_file" option and
   keys them to points by a per-point field-label string plus a per-LC star
   name (the same machinery -stitch uses on input).

   Sign convention: -stitch apply does x_new = x_old - shiftvalue, and the
   shifts file stores +shiftvalue, so -unstitch does x_old = x_cur + value.
*/

#define UNSTITCH_SOURCE_INSHIFTSFILE 0
#define UNSTITCH_SOURCE_FITSHEADER   1

typedef struct {
  int nvar;                 /* number of variables to un-shift */

  int source;               /* UNSTITCH_SOURCE_INSHIFTSFILE or _FITSHEADER */

  /* The magnitude variable(s) to un-shift (modified in place) */
  char **varnames;
  double ***varvals;
  _Variable **vars;

  /* "fitsheader" source: read the shifts from keywords written by -stitch's
     add_shifts_fitsheader option into the input FITS header.  Points are
     keyed to shifts by lcnum (and refnum, if given), matched against the
     keyword comments ("Shift for variable V LCgroup L [REFID R]").  The
     header stores -shiftvalue, so unstitch subtracts the keyword value. */
  char keywordbase[MAXLEN];
  int **fits_lcnumval;
  _Variable *fits_lcnumvar;
  int fits_userefnum;
  int **fits_refnumval;
  _Variable *fits_refnumvar;
  int fits_hdu;             /* 1 = primary header, 2 = first extension */

  /* Optional mask: points with mask > VARTOOLS_MASK_TINY are "unmasked".
     Only unmasked points are required to have a shift (coverage check);
     masked points are exempt and left unchanged if they have no shift. */
  int usemask;
  char *maskname;
  double **maskvals;
  _Variable *maskvar;
  /* When set (and maskpoints is given), masked points are left completely
     unchanged -- never shifted, even if they match a shift.  This inverts a
     -stitch run that used its "noshiftmasked" option.  Without it, masked
     points that match a shift ARE shifted (inverting a default -stitch). */
  int noshiftmasked;

  /* Optional: strip the stitch keywords from the output FITS header.  Removes
     every keyword whose name starts with strip_keywordbase (and, if
     strip_stitchparams is set, the fixed STCH* stitch-parameter keywords) from
     the chosen HDU of any FITS light curve subsequently written. */
  int strip_fitsheader;
  char strip_keywordbase[MAXLEN];
  int strip_stitchparams;
  int strip_hdutouse;       /* 0 = primary, 1 = first extension */

  /* Per-point field-label string + optional refnum appended to it */
  char ***field_labels_vals;
  int **field_labels_vals_indx;
  int *field_labels_vals_issorted;
  _Variable *field_labels_var;
  int is_append_refnum_to_fieldlabel;
  int **refnumval;
  _Variable *refnumvar;

  /* Per-LC star name used to select the matching line of the shifts file */
  char **starname_vals;
  _Variable *starname_var;
  char *starname_varname;

  /* Shifts read from the input shifts file(s) (one file per variable).
     Mirrors the in_shift_* layout of -stitch. */
  char **in_shifts_filename;
  int *N_in_shift_stars;
  int *size_in_shifts_file;
  int **N_shifts_per_star;
  char ***in_shift_starnames;
  int **in_shift_stars_found;
  int **in_shift_starnames_sortidx;
  char ****in_shift_labels;
  int ***in_shift_labels_sortidx;
  double ***in_shift_values;
  char ****in_shift_values_str;
  int ***Nobs_in_shifts;

  /* Output column: number of points that received a shift (summed over the
     variables) for each light curve. */
  int *Npoints_shifted;

} _Unstitch;
