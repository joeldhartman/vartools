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

typedef struct {
  int nvar;                 /* number of variables to un-shift */

  /* The magnitude variable(s) to un-shift (modified in place) */
  char **varnames;
  double ***varvals;
  _Variable **vars;

  /* Optional mask: points with mask > VARTOOLS_MASK_TINY are "unmasked".
     Only unmasked points are required to have a shift (coverage check);
     masked points are exempt and left unchanged if they have no shift. */
  int usemask;
  char *maskname;
  double **maskvals;
  _Variable *maskvar;

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
