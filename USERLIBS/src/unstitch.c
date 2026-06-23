#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "../../src/vartools.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "unstitch.h"

/* The "fitsheader" source reads shift keywords directly from the input FITS
   header via cfitsio.  The cfitsio symbols are resolved at run time from
   libvartoolspipeline (which has cfitsio statically linked and exports its
   symbols); only the header is needed at compile time.  Everything cfitsio
   is guarded by HAVE_CFITSIO so the library still builds without it. */
#ifdef HAVE_CFITSIO
#include "fitsio.h"
#endif

/* User-defined command -unstitch: the inverse of -stitch.  It adds
   previously-determined per-segment shifts back to a light curve to restore
   the original (pre-stitch) magnitudes.

   Phase A implements the "in_shifts_file" source: shifts are read from a file
   written by -stitch's "out_shifts_file" option, keyed to points by a
   per-point field-label string + a per-LC star name.  Later phases will add a
   "fitsheader" source (read shifts from the input FITS header) and an option
   to strip the stitch keywords from the output FITS header.
*/

void DoUnstitch(ProgramData *p, _Unstitch *u, int lc_name_num, int lc_num);
void DoUnstitch_FitsHeader(ProgramData *p, _Unstitch *u, int lc_name_num, int lc_num);
void unstitch_ReadInshifts_File(ProgramData *p, Command *c, _Unstitch *u, int vv);

void unstitch_Initialize(char *commandname,
			 int *RequireReadAll,
			 int *RequireSortLC,
			 int *RequireDistinctTimes,
			 size_t *sizeuserdata)
{
  sprintf(commandname,"-unstitch");
  *RequireReadAll = 0;
  *RequireSortLC = 0;
  *RequireDistinctTimes = 0;
  *sizeuserdata = sizeof(_Unstitch);
}

/* Parse a comma-separated list in argv[i] into a freshly-allocated array of
   strings; returns the count via *nret and the array via *namesret. */
static int unstitch_parse_strlist(char *arg, char ***namesret, int *nret)
{
  int j, k, i1, i2, n;
  char **names;
  n = 1;
  j = 0;
  while(arg[j] != '\0') {
    if(arg[j] == ',') n++;
    j++;
  }
  if((names = (char **) malloc(n * sizeof(char *))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  i1 = 0; i2 = 0;
  for(k = 0; k < n; k++) {
    while(arg[i2] != '\0' && arg[i2] != ',') i2++;
    if((names[k] = (char *) malloc((i2 - i1 + 1)*sizeof(char))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    for(j = i1; j < i2; j++) names[k][j-i1] = arg[j];
    names[k][i2-i1] = '\0';
    i1 = i2+1; i2 = i2+1;
  }
  *namesret = names;
  *nret = n;
  return 0;
}

int unstitch_ParseCL(ProgramData *p, Command *c,
		     void *userdata, int *iret, char **argv, int argc)
/* Expected syntax:
   -unstitch unstitch_variable_list
       "in_shifts_file" fieldlabelsvar starnamevar file1[,file2,...]
           ["append_refnum_to_fieldlabel" refnum_var]
       ["maskpoints" maskvar]
*/
{
  int i, j, k, n;
  _Unstitch *u;
  u = (_Unstitch *) userdata;

  i = *iret;
  if(i >= argc) return 1;

  /* unstitch_variable_list */
  unstitch_parse_strlist(argv[i], &(u->varnames), &(u->nvar));
  if((u->varvals = (double ***) malloc(u->nvar * sizeof(double **))) == NULL ||
     (u->vars = (_Variable **) malloc(u->nvar * sizeof(_Variable *))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  for(k = 0; k < u->nvar; k++) {
    VARTOOLS_RegisterDataVector(p, c, (void *) (&(u->varvals[k])),
				VARTOOLS_TYPE_DOUBLE, 0,
				VARTOOLS_SOURCE_EXISTINGVARIABLE,
				0, NULL, u->varnames[k],
				(char) VARTOOLS_VECTORTYPE_LC,
				&(u->vars[k]));
  }

  /* Source dispatch: "in_shifts_file" (read shifts from a file) or
     "fitsheader" (read shifts from the input FITS header). */
  i++;
  if(i >= argc) return 1;
  if(!strcmp(argv[i],"in_shifts_file")) {
  u->source = UNSTITCH_SOURCE_INSHIFTSFILE;

  /* fieldlabelsvar */
  i++;
  if(i >= argc) return 1;
  VARTOOLS_RegisterDataVector(p, c, (void *) (&(u->field_labels_vals)),
			      VARTOOLS_TYPE_STRING, 0,
			      VARTOOLS_SOURCE_EXISTINGVARIABLE,
			      0, NULL, argv[i],
			      (char) VARTOOLS_VECTORTYPE_LC,
			      &(u->field_labels_var));
  VARTOOLS_RegisterDataVector(p, c, (void *) (&(u->field_labels_vals_indx)),
			      VARTOOLS_TYPE_INT, 0,
			      VARTOOLS_SOURCE_LC,
			      0, NULL, NULL, -1, NULL, NULL);
  VARTOOLS_RegisterDataVector(p, c, (void *) (&(u->field_labels_vals_issorted)),
			      VARTOOLS_TYPE_INT, 0,
			      VARTOOLS_SOURCE_EVALEXPRESSION,
			      0, NULL, "0");

  /* starnamevar */
  i++;
  if(i >= argc) return 1;
  VARTOOLS_RegisterDataVector(p, c, (void *) (&(u->starname_vals)),
			      VARTOOLS_TYPE_STRING, 0,
			      VARTOOLS_SOURCE_EXISTINGVARIABLE,
			      0, NULL, argv[i],
			      (char) VARTOOLS_VECTORTYPE_PERSTARDATA,
			      &(u->starname_var));
  if((u->starname_varname = (char *) malloc((strlen(argv[i])+1)*sizeof(char))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  sprintf(u->starname_varname,"%s",argv[i]);

  /* file list (one per variable) */
  i++;
  if(i >= argc) return 1;
  n = 1;
  j = 0;
  while(argv[i][j] != '\0') { if(argv[i][j] == ',') n++; j++; }
  if(n != u->nvar) {
    fprintf(stderr,"Error parsing the unstitch command - the number of files in the in_shifts_file list does not match the number in the unstitch_variable_list.\n");
    return 1;
  }
  if((u->in_shifts_filename = (char **) malloc(u->nvar*sizeof(char *))) == NULL ||
     (u->N_in_shift_stars = (int *) malloc(u->nvar*sizeof(int))) == NULL ||
     (u->N_shifts_per_star = (int **) malloc(u->nvar*sizeof(int *))) == NULL ||
     (u->size_in_shifts_file = (int *) malloc(u->nvar*sizeof(int))) == NULL ||
     (u->in_shift_starnames = (char ***) malloc(u->nvar*sizeof(char **))) == NULL ||
     (u->in_shift_stars_found = (int **) malloc(u->nvar*sizeof(int *))) == NULL ||
     (u->in_shift_starnames_sortidx = (int **) malloc(u->nvar*sizeof(int *))) == NULL ||
     (u->in_shift_labels = (char ****) malloc(u->nvar*sizeof(char ***))) == NULL ||
     (u->in_shift_labels_sortidx = (int ***) malloc(u->nvar*sizeof(int **))) == NULL ||
     (u->in_shift_values = (double ***) malloc(u->nvar*sizeof(double **))) == NULL ||
     (u->in_shift_values_str = (char ****) malloc(u->nvar*sizeof(char ***))) == NULL ||
     (u->Nobs_in_shifts = (int ***) malloc(u->nvar*sizeof(int **))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  {
    char **flist;
    int nf;
    unstitch_parse_strlist(argv[i], &flist, &nf);
    for(k = 0; k < u->nvar; k++) {
      u->in_shifts_filename[k] = flist[k];
      unstitch_ReadInshifts_File(p, c, u, k);
    }
    free(flist);
  }

  /* optional "append_refnum_to_fieldlabel" refnum_var */
  u->is_append_refnum_to_fieldlabel = 0;
  i++;
  if(i < argc && !strcmp(argv[i],"append_refnum_to_fieldlabel")) {
    u->is_append_refnum_to_fieldlabel = 1;
    i++;
    if(i >= argc) return 1;
    VARTOOLS_RegisterDataVector(p, c, (void *) (&(u->refnumval)),
				VARTOOLS_TYPE_INT, 0,
				VARTOOLS_SOURCE_EXISTINGVARIABLE,
				0, NULL, argv[i],
				(char) VARTOOLS_VECTORTYPE_LC,
				&(u->refnumvar));
  } else {
    i--;
  }

  } /* end of "in_shifts_file" source */
  else if(!strcmp(argv[i],"fitsheader")) {
#ifndef HAVE_CFITSIO
    fprintf(stderr,"Error parsing the unstitch command - the \"fitsheader\" source requires that vartools was compiled with cfitsio support, which is not the case for this build.\n");
    return 1;
#else
    u->source = UNSTITCH_SOURCE_FITSHEADER;

    /* keywordbase */
    i++;
    if(i >= argc) return 1;
    sprintf(u->keywordbase,"%s",argv[i]);

    /* lcnum_var */
    i++;
    if(i >= argc) return 1;
    VARTOOLS_RegisterDataVector(p, c, (void *) (&(u->fits_lcnumval)),
				VARTOOLS_TYPE_INT, 0,
				VARTOOLS_SOURCE_EXISTINGVARIABLE,
				0, NULL, argv[i],
				(char) VARTOOLS_VECTORTYPE_LC,
				&(u->fits_lcnumvar));

    /* optional "refnum_var" refnum_var */
    u->fits_userefnum = 0;
    i++;
    if(i < argc && !strcmp(argv[i],"refnum_var")) {
      u->fits_userefnum = 1;
      i++;
      if(i >= argc) return 1;
      VARTOOLS_RegisterDataVector(p, c, (void *) (&(u->fits_refnumval)),
				  VARTOOLS_TYPE_INT, 0,
				  VARTOOLS_SOURCE_EXISTINGVARIABLE,
				  0, NULL, argv[i],
				  (char) VARTOOLS_VECTORTYPE_LC,
				  &(u->fits_refnumvar));
    } else {
      i--;
    }

    /* optional "primary" | "extension" (default primary) */
    u->fits_hdu = 1;
    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"primary")) u->fits_hdu = 1;
      else if(!strcmp(argv[i],"extension")) u->fits_hdu = 2;
      else i--;
    } else {
      i--;
    }
#endif
  }
  else {
    return 1;
  }

  /* optional "maskpoints" maskvar */
  u->usemask = 0;
  i++;
  if(i < argc && !strcmp(argv[i],"maskpoints")) {
    u->usemask = 1;
    i++;
    if(i >= argc) return 1;
    VARTOOLS_RegisterDataVector(p, c, (void *) (&(u->maskvals)),
				VARTOOLS_TYPE_DOUBLE, 0,
				VARTOOLS_SOURCE_EXISTINGVARIABLE,
				0, NULL, argv[i],
				(char) VARTOOLS_VECTORTYPE_LC,
				&(u->maskvar));
  } else {
    i--;
  }

  /* optional "noshiftmasked" */
  u->noshiftmasked = 0;
  i++;
  if(i < argc && !strcmp(argv[i],"noshiftmasked")) {
    u->noshiftmasked = 1;
  } else {
    i--;
  }
  if(u->noshiftmasked && !u->usemask) {
    fprintf(stderr,"Error parsing the unstitch command - the \"noshiftmasked\" keyword requires that \"maskpoints\" is also given (there is otherwise no mask to identify masked points).\n");
    return 1;
  }

  /* optional "strip_fitsheader" keywordbase ["stitchparams"] ["primary"|"extension"] */
  u->strip_fitsheader = 0;
  u->strip_stitchparams = 0;
  u->strip_hdutouse = 0;
  i++;
  if(i < argc && !strcmp(argv[i],"strip_fitsheader")) {
    u->strip_fitsheader = 1;
    i++;
    if(i >= argc) return 1;
    sprintf(u->strip_keywordbase,"%s",argv[i]);
    i++;
    if(i < argc && !strcmp(argv[i],"stitchparams")) {
      u->strip_stitchparams = 1;
    } else {
      i--;
    }
    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"primary")) u->strip_hdutouse = 0;
      else if(!strcmp(argv[i],"extension")) u->strip_hdutouse = 1;
      else i--;
    } else {
      i--;
    }
  } else {
    i--;
  }

  /* Output column: number of points un-shifted per light curve */
  VARTOOLS_RegisterDataVector(p, c, (void *) (&(u->Npoints_shifted)),
			      VARTOOLS_TYPE_INT, 0, VARTOOLS_SOURCE_COMPUTED,
			      1, "Npoints_shifted");

  i++;
  *iret = i;
  return 0;
}

void unstitch_ShowSyntax(FILE *outfile)
{
  OutText s;
  s.s = NULL; s.space = 0; s.len_s = 0; s.Nchar_cur_line = 0;
  VARTOOLS_printtostring(&s,"-unstitch\n");
  VARTOOLS_printtostring(&s,"\tunstitch_variable_list\n");
  VARTOOLS_printtostring(&s,"\t<\"in_shifts_file\" fieldlabelsvar starnamevar file1[,file2,...]\n");
  VARTOOLS_printtostring(&s,"\t\t[\"append_refnum_to_fieldlabel\" refnum_var]\n");
  VARTOOLS_printtostring(&s,"\t | \"fitsheader\" keywordbase lcnum_var [\"refnum_var\" refnum_var]\n");
  VARTOOLS_printtostring(&s,"\t\t[\"primary\" | \"extension\"]>\n");
  VARTOOLS_printtostring(&s,"\t[\"maskpoints\" maskvar [\"noshiftmasked\"]]\n");
  VARTOOLS_printtostring(&s,"\t[\"strip_fitsheader\" keywordbase [\"stitchparams\"]\n");
  VARTOOLS_printtostring(&s,"\t\t[\"primary\" | \"extension\"]]\n");
  fprintf(outfile,"%s",s.s);
}

void unstitch_ShowHelp(FILE *outfile)
{
  OutText s;
  s.s = NULL; s.space = 0; s.len_s = 0; s.Nchar_cur_line = 0;
  VARTOOLS_printtostring(&s,"Undo a previous -stitch operation by adding the determined per-segment shifts back to the light curve, restoring the original magnitudes. The shifts are read either from a file produced by the \"out_shifts_file\" option of -stitch, or from the keywords that the \"add_shifts_fitsheader\" option of -stitch wrote into the FITS header of the input light curve.\n\n");
  VARTOOLS_printtostring(&s,"unstitch_variable_list - A comma-separated list of light curve variables to un-shift (typically mag). For the \"in_shifts_file\" source use one input shifts file per variable, in the same order.\n\n");
  VARTOOLS_printtostring(&s,"\"in_shifts_file\" - Read the shifts from a file. First provide the name of a light curve variable storing a string field identifier (or telescope ID) for each observation; this is used to match points to the shifts in the file. Then provide the name of an input list variable storing the string star name for each light curve, used to select the matching line of the file. Finally provide a comma-separated list of shifts files to read, one per variable being un-shifted. The file format is that written by the \"out_shifts_file\" option of -stitch. The stored shift is added back to each matching point (the inverse of how -stitch subtracts it).\n\n");
  VARTOOLS_printtostring(&s,"\"append_refnum_to_fieldlabel\" - If the shifts file was written with this option, give it here as well, followed by the refnum_var, so the field labels are reconstructed in the same way for matching.\n\n");
  VARTOOLS_printtostring(&s,"\"fitsheader\" - Read the shifts from the FITS header of the input light curve, as written by the \"add_shifts_fitsheader\" option of -stitch. Give the keyword basename that was used (e.g. SHFT), followed by the lcnum_var that identifies the light curve segment for each point. Optionally give \"refnum_var\" followed by the refnum variable if the shifts were determined with a refnum_var. By default the keywords are read from the primary header; give \"extension\" to read them from the first extension header instead. Each point is matched to a shift by its lcnum (and refnum) value against the keyword comments; the header stores the shift as the change that was applied, so the inverse (the original magnitude) is recovered by subtracting it. This source requires that vartools was compiled with cfitsio support and that the input light curve is in FITS format. Note that the shift is recovered to the precision stored in the FITS keyword, which may be slightly less than the full machine precision available through \"in_shifts_file\".\n\n");
  VARTOOLS_printtostring(&s,"\"maskpoints\" - Optionally give a mask variable. Points with a mask value greater than 0 are treated as in-use (unmasked); points with a mask value of 0 (or negative) are masked. Masked points are exempt from the coverage check below. By default a masked point that has a matching shift is still shifted (this correctly inverts a default -stitch, which shifts all points); a masked point with no matching shift is left unchanged.\n\n");
  VARTOOLS_printtostring(&s,"\"noshiftmasked\" - Only valid together with \"maskpoints\". When given, masked points are left completely unchanged -- they are never shifted, even if they match a shift. Use this to invert a -stitch run that used its own \"noshiftmasked\" option (so that masked points were not shifted by -stitch); otherwise -unstitch would over-correct those points.\n\n");
  VARTOOLS_printtostring(&s,"\"strip_fitsheader\" - Optionally remove the stitch keywords from the header of any light curve subsequently output in FITS format (for example, with -o ... fits copyheader, which carries the input FITS header through to the output). Give the keyword basename that -stitch used; every keyword whose name begins with it is removed. Optionally give \"stitchparams\" to also remove the fixed STCH* stitch-parameter keywords written by -stitch's add_stitchparams_fitsheader option. By default the keywords are removed from the primary header; give \"extension\" to remove them from the first extension header instead. This is a no-op if the output is not FITS or the keywords are not present, and is independent of which shift source is used.\n\n");
  VARTOOLS_printtostring(&s,"Coverage: -unstitch requires that every unmasked point has a matching shift. If any unmasked point has no matching shift (its star/field label is missing from the file, or no matching keyword is present in the header) the command quits with an error, since this indicates that the wrong set of shifts is being used to un-stitch the light curve.\n\n");
  VARTOOLS_printtostring(&s,"Output: the column Npoints_shifted reports the number of points that received a shift (summed over the variables) for each light curve.\n\n");
  fprintf(outfile,"%s",s.s);
}

void unstitch_ShowExample(FILE *outfile)
{
  fprintf(outfile,
	  "\nvartools -L USERLIBS/src/.libs/stitch.so \\\n"
	  "    -l EXAMPLES/lc_list_unstitch combinelcs lcnumvar lcnum \\\n"
	  "    -inlistvars 'field:2:combinelc:string,star:3:string' \\\n"
	  "    -expr 'mask=mag*0+1' \\\n"
	  "    -stitch mag err mask lcnum median shifts_file field star \\\n"
	  "        out_shifts_file EXAMPLES/OUTDIR1/shifts.txt -oneline\n\n"
	  "vartools -L USERLIBS/src/.libs/stitch.so -L USERLIBS/src/.libs/unstitch.so \\\n"
	  "    -l EXAMPLES/lc_list_unstitch combinelcs lcnumvar lcnum \\\n"
	  "    -inlistvars 'field:2:combinelc:string,star:3:string' \\\n"
	  "    -expr 'mask=mag*0+1' \\\n"
	  "    -rms \\\n"
	  "    -stitch mag err mask lcnum median \\\n"
	  "    -rms \\\n"
	  "    -unstitch mag in_shifts_file field star EXAMPLES/OUTDIR1/shifts.txt \\\n"
	  "    -rms -oneline\n\n"
	  "The list file EXAMPLES/lc_list_unstitch combines two light-curve segments (EXAMPLES/2 and EXAMPLES/2.shifted, the latter being EXAMPLES/2 with +0.3 mag added) into a single light curve, assigning per-segment field labels fA and fB (through the \"combinelc\" form of -inlistvars) and the star name star1. The first command stitches the segments together and writes the determined per-segment shifts to EXAMPLES/OUTDIR1/shifts.txt. The second command stitches again and then undoes it with -unstitch, reading those same shifts: the three -rms outputs show that the combined scatter is inflated by the inter-segment offset, drops after stitching, and returns to the inflated value after un-stitching, confirming that -unstitch restores the original magnitudes. In normal use the -stitch and -unstitch steps are run on separate occasions, with EXAMPLES/OUTDIR1/shifts.txt carried between them.\n");
}

void unstitch_RunCommand(ProgramData *p, void *userdata, int lc_name_num, int lc_num)
{
  _Unstitch *u;
  u = (_Unstitch *) userdata;
  if(u->source == UNSTITCH_SOURCE_FITSHEADER)
    DoUnstitch_FitsHeader(p, u, lc_name_num, lc_num);
  else
    DoUnstitch(p, u, lc_name_num, lc_num);

  /* Queue removal of the stitch keywords from the output FITS header (applied
     when the light curve is written; a no-op if the output is not FITS or the
     keywords are absent). */
  if(u->strip_fitsheader) {
    VARTOOLS_Delete_Keyword_From_OutputLC_FitsHeader(p, lc_num,
						     u->strip_keywordbase,
						     u->strip_hdutouse, 1);
    if(u->strip_stitchparams)
      VARTOOLS_Delete_Keyword_From_OutputLC_FitsHeader(p, lc_num, "STCH",
						       u->strip_hdutouse, 1);
  }
}

/* ------------------------------------------------------------------------
   Input shifts file reading and matching.  These mirror the corresponding
   -stitch routines (stitch_ParseInShifts_Line / stitch_ReadInshifts_File /
   stitch_Find_Star_in_InshiftFile / stitch_Find_FieldLabels_in_InshiftFile),
   adapted to the _Unstitch structure.
   ------------------------------------------------------------------------ */

int unstitch_ParseInShifts_Line(char *line, char **incols, int *size_incols)
{
  int j, i;
  j = 0; i = 0;
  while((line[j] == ' ' || line[j] == '\t') && line[j] != '\0' && line[j] != '\n')
    j++;
  if(line[j] == '#' || line[j] == '\0')
    return 1;
  while (i < 2 && line[j] != '\0' && line[j] != '\n') {
    if(!size_incols[i]) {
      if((incols[i] = (char *) malloc((strlen(&(line[j]))+1)*sizeof(char))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
      size_incols[i] = strlen(&(line[j]))+1;
    } else if(size_incols[i] < strlen(&line[j])+1) {
      if((incols[i] = (char *) realloc(incols[i], (strlen(&(line[j]))+1)*sizeof(char))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
      size_incols[i] = strlen(&(line[j]))+1;
    }
    j += VARTOOLS_parseone(&(line[j]),(void *) ((incols[i])), VARTOOLS_TYPE_STRING);
    i++;
  }
  if(i < 2)
    return 1;
  return 0;
}

void unstitch_ReadInshifts_File(ProgramData *p, Command *c, _Unstitch *u, int vv)
{
  FILE *shiftfile;
  char *inputline = NULL;
  size_t size_inputline = 0;
  char **incols;
  int *size_incols;
  int i, j, k, mm, ll, jj, testskip;

  if((shiftfile = fopen(u->in_shifts_filename[vv],"r")) == NULL) {
    fprintf(stderr,"Cannot open the input shifts file: %s\n", u->in_shifts_filename[vv]);
    exit(1);
  }
  u->N_in_shift_stars[vv] = 0;
  u->size_in_shifts_file[vv] = 256;
  if((u->N_shifts_per_star[vv] = (int *) malloc(u->size_in_shifts_file[vv]*sizeof(int))) == NULL ||
     (u->in_shift_starnames[vv] = (char **) malloc(u->size_in_shifts_file[vv]*sizeof(char *))) == NULL ||
     (u->in_shift_stars_found[vv] = (int *) malloc(u->size_in_shifts_file[vv]*sizeof(int))) == NULL ||
     (u->in_shift_starnames_sortidx[vv] = (int *) malloc(u->size_in_shifts_file[vv]*sizeof(int))) == NULL ||
     (u->in_shift_labels[vv] = (char ***) malloc(u->size_in_shifts_file[vv]*sizeof(char **))) == NULL ||
     (u->in_shift_labels_sortidx[vv] = (int **) malloc(u->size_in_shifts_file[vv]*sizeof(int *))) == NULL ||
     (u->in_shift_values[vv] = (double **) malloc(u->size_in_shifts_file[vv]*sizeof(double *))) == NULL ||
     (u->in_shift_values_str[vv] = (char ***) malloc(u->size_in_shifts_file[vv]*sizeof(char **))) == NULL ||
     (u->Nobs_in_shifts[vv] = (int **) malloc(u->size_in_shifts_file[vv]*sizeof(int *))) == NULL ||
     (incols = (char **) malloc(2*sizeof(char *))) == NULL ||
     (size_incols = (int *) malloc(2*sizeof(int))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  incols[0] = NULL; incols[1] = NULL;
  size_incols[0] = 0; size_incols[1] = 0;
  i = 0;
  while(VARTOOLS_gnu_getline(&(inputline),&(size_inputline), shiftfile) >= 0) {
    if(inputline[0] != '#') {
      u->N_in_shift_stars[vv] += 1;
      if(u->N_in_shift_stars[vv] >= u->size_in_shifts_file[vv]) {
	u->size_in_shifts_file[vv] *= 2;
	if((u->N_shifts_per_star[vv] = (int *) realloc(u->N_shifts_per_star[vv], u->size_in_shifts_file[vv]*sizeof(int))) == NULL ||
	   (u->in_shift_starnames[vv] = (char **) realloc(u->in_shift_starnames[vv], u->size_in_shifts_file[vv]*sizeof(char *))) == NULL ||
	   (u->in_shift_stars_found[vv] = (int *) realloc(u->in_shift_stars_found[vv], u->size_in_shifts_file[vv]*sizeof(int))) == NULL ||
	   (u->in_shift_starnames_sortidx[vv] = (int *) realloc(u->in_shift_starnames_sortidx[vv], u->size_in_shifts_file[vv]*sizeof(int))) == NULL ||
	   (u->in_shift_labels[vv] = (char ***) realloc(u->in_shift_labels[vv], u->size_in_shifts_file[vv]*sizeof(char **))) == NULL ||
	   (u->in_shift_labels_sortidx[vv] = (int **) realloc(u->in_shift_labels_sortidx[vv], u->size_in_shifts_file[vv]*sizeof(int *))) == NULL ||
	   (u->in_shift_values[vv] = (double **) realloc(u->in_shift_values[vv], u->size_in_shifts_file[vv]*sizeof(double *))) == NULL ||
	   (u->in_shift_values_str[vv] = (char ***) realloc(u->in_shift_values_str[vv], u->size_in_shifts_file[vv]*sizeof(char **))) == NULL ||
	   (u->Nobs_in_shifts[vv] = (int **) realloc(u->Nobs_in_shifts[vv], u->size_in_shifts_file[vv]*sizeof(int *))) == NULL)
	  VARTOOLS_error(ERR_MEMALLOC);
      }
      testskip = unstitch_ParseInShifts_Line(inputline, incols, size_incols);
      if(testskip) {
	fprintf(stderr,"Error parsing input shifts file: %s\n", u->in_shifts_filename[vv]);
	fclose(shiftfile);
	exit(1);
      }
      if((u->in_shift_starnames[vv][i] = (char *) malloc((strlen(incols[0])+1)*sizeof(char))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
      sprintf(u->in_shift_starnames[vv][i],"%s",incols[0]);
      u->in_shift_starnames_sortidx[vv][i] = i;
      u->in_shift_stars_found[vv][i] = 0;
      j = 0; k = 0; mm = 0;
      while(incols[1][j] != '\0' && incols[1][j] != '\n') {
	if(incols[1][j] == ',') mm++;
	if(incols[1][j] == ';') {
	  if(mm != 2) {
	    fprintf(stderr,"Error parsing input shifts file: %s\n", u->in_shifts_filename[vv]);
	    fclose(shiftfile);
	    exit(1);
	  }
	  k++; mm = 0;
	}
	j++;
      }
      if(mm != 2) {
	fprintf(stderr,"Error parsing input shifts file: %s\n", u->in_shifts_filename[vv]);
	fclose(shiftfile);
	exit(1);
      }
      u->N_shifts_per_star[vv][i] = k+1;
      if((u->in_shift_labels[vv][i] = (char **) malloc(u->N_shifts_per_star[vv][i]*sizeof(char *))) == NULL ||
	 (u->in_shift_labels_sortidx[vv][i] = (int *) malloc(u->N_shifts_per_star[vv][i]*sizeof(int))) == NULL ||
	 (u->in_shift_values[vv][i] = (double *) malloc(u->N_shifts_per_star[vv][i]*sizeof(double))) == NULL ||
	 (u->in_shift_values_str[vv][i] = (char **) malloc(u->N_shifts_per_star[vv][i]*sizeof(char *))) == NULL ||
	 (u->Nobs_in_shifts[vv][i] = (int *) malloc(u->N_shifts_per_star[vv][i]*sizeof(int))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
      j = 0; jj = 0; k = 0; mm = 0;
      while(incols[1][j] != '\0' && incols[1][j] != '\n') {
	if(incols[1][j] == ',') {
	  if(!mm) {
	    if((u->in_shift_labels[vv][i][k] = (char *) malloc((j - jj + 1)*sizeof(char))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(ll = jj; ll < j; ll++)
	      u->in_shift_labels[vv][i][k][ll-jj] = incols[1][ll];
	    u->in_shift_labels[vv][i][k][ll-jj] = '\0';
	    jj = j+1;
	  }
	  else if(mm == 1) {
	    incols[1][j] = '\0';
	    u->in_shift_values[vv][i][k] = atof(&(incols[1][jj]));
	    if((u->in_shift_values_str[vv][i][k] = (char *) malloc((j - jj + 1)*sizeof(char))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(ll = jj; ll < j; ll++)
	      u->in_shift_values_str[vv][i][k][ll-jj] = incols[1][ll];
	    u->in_shift_values_str[vv][i][k][ll-jj] = '\0';
	    jj = j+1;
	  }
	  mm++;
	}
	if(incols[1][j] == ';') {
	  incols[1][j] = '\0';
	  u->Nobs_in_shifts[vv][i][k] = atoi(&(incols[1][jj]));
	  k++; mm = 0; jj = j+1;
	}
	j++;
      }
      incols[1][j] = '\0';
      u->Nobs_in_shifts[vv][i][k] = atoi(&(incols[1][jj]));
      for(mm = 0; mm < u->N_shifts_per_star[vv][i]; mm++)
	u->in_shift_labels_sortidx[vv][i][mm] = mm;
      VARTOOLS_mysortstringint(u->N_shifts_per_star[vv][i], 0, u->in_shift_labels[vv][i], u->in_shift_labels_sortidx[vv][i]);
      i++;
    }
  }
  VARTOOLS_mysortstringint(u->N_in_shift_stars[vv], 0, u->in_shift_starnames[vv], u->in_shift_starnames_sortidx[vv]);
  fclose(shiftfile);
  if(incols[0] != NULL) free(incols[0]);
  if(incols[1] != NULL) free(incols[1]);
  free(incols);
  free(size_incols);
  if(inputline != NULL) free(inputline);
  return;
}

int unstitch_Find_Star_in_InshiftFile(_Unstitch *u, ProgramData *p, int lcnum, int lc_name_num, int vv)
{
  int i;
  if(!u->N_in_shift_stars[vv]) return -1;
  i = VARTOOLS_findX_string(u->in_shift_starnames[vv], u->in_shift_starnames_sortidx[vv], u->starname_vals[lcnum], 0, u->N_in_shift_stars[vv]);
  if(i >= 0) {
    i = u->in_shift_starnames_sortidx[vv][i];
    u->in_shift_stars_found[vv][i] = 1;
  }
  return i;
}

void unstitch_Find_FieldLabels_in_InshiftFile(_Unstitch *u, ProgramData *p, int lcnum, int lc_name_num, int vv, int *fieldlabels_indx, int shiftid)
{
  int i, j, jlastfound;
  if(!u->field_labels_vals_issorted[lcnum]) {
    for(i=0; i < p->NJD[lcnum]; i++) u->field_labels_vals_indx[lcnum][i] = i;
    VARTOOLS_mysortstringint(p->NJD[lcnum], 0, u->field_labels_vals[lcnum], u->field_labels_vals_indx[lcnum]);
    u->field_labels_vals_issorted[lcnum] = 1;
  }
  jlastfound = -1;
  for(i=0; i < p->NJD[lcnum]; i++) {
    if(jlastfound < 0) {
      j = VARTOOLS_findX_string(u->in_shift_labels[vv][shiftid], u->in_shift_labels_sortidx[vv][shiftid], u->field_labels_vals[lcnum][u->field_labels_vals_indx[lcnum][i]], 0, u->N_shifts_per_star[vv][shiftid]);
      if(j >= 0) {
	jlastfound = j;
	j = u->in_shift_labels_sortidx[vv][shiftid][j];
      }
    } else {
      j = VARTOOLS_findX_string(u->in_shift_labels[vv][shiftid], u->in_shift_labels_sortidx[vv][shiftid], u->field_labels_vals[lcnum][u->field_labels_vals_indx[lcnum][i]], jlastfound, u->N_shifts_per_star[vv][shiftid]);
      if(j >= 0) {
	jlastfound = j;
	j = u->in_shift_labels_sortidx[vv][shiftid][j];
      }
    }
    fieldlabels_indx[u->field_labels_vals_indx[lcnum][i]] = j;
  }
}

/* Add the matched shift back to each point (the inverse of -stitch).  A point
   is "covered" if its field label matched an entry in the file.  Quit with an
   error if any unmasked point is uncovered. */
void unstitch_Apply_InShifts(_Unstitch *u, ProgramData *p, int lcnum, int lc_name_num, int vv, int *fieldlabels_indx, int shiftid)
{
  int i;
  for(i=0; i < p->NJD[lcnum]; i++) {
    if(fieldlabels_indx[i] >= 0) {
      if(u->noshiftmasked && u->usemask && u->maskvals[lcnum][i] <= VARTOOLS_MASK_TINY)
	continue;   /* masked point left unchanged */
      u->varvals[vv][lcnum][i] += u->in_shift_values[vv][shiftid][fieldlabels_indx[i]];
      u->Npoints_shifted[lcnum] += 1;
    } else if(!u->usemask || u->maskvals[lcnum][i] > VARTOOLS_MASK_TINY) {
      fprintf(stderr,"Error in -unstitch: no shift found for field \"%s\" of star \"%s\" (variable %s) in the input shifts file %s. The supplied shifts do not fully cover the light curve; check that the correct shifts file is being used.\n",
	      u->field_labels_vals[lcnum][i], u->starname_vals[lcnum],
	      u->varnames[vv], u->in_shifts_filename[vv]);
      exit(1);
    }
  }
}

void DoUnstitch(ProgramData *p, _Unstitch *u, int lc_name_num, int lc_num)
{
  int vv, i;
  int inshift_starid;
  int *fieldlabels_indx = NULL;

  u->Npoints_shifted[lc_num] = 0;

  if(p->NJD[lc_num] <= 0)
    return;

  /* Reconstruct the field labels exactly as -stitch wrote them, if needed. */
  if(u->is_append_refnum_to_fieldlabel) {
    for(i=0; i < p->NJD[lc_num]; i++) {
      sprintf(&(u->field_labels_vals[lc_num][i][strlen(u->field_labels_vals[lc_num][i])]),"/%d",u->refnumval[lc_num][i]);
    }
  }

  for(vv = 0; vv < u->nvar; vv++) {
    inshift_starid = unstitch_Find_Star_in_InshiftFile(u, p, lc_num, lc_name_num, vv);
    if(inshift_starid < 0) {
      /* Star absent from the file: an error if any unmasked point exists. */
      for(i=0; i < p->NJD[lc_num]; i++) {
	if(!u->usemask || u->maskvals[lc_num][i] > VARTOOLS_MASK_TINY) {
	  fprintf(stderr,"Error in -unstitch: star \"%s\" is not present in the input shifts file %s (variable %s). The supplied shifts do not cover this light curve; check that the correct shifts file is being used.\n",
		  u->starname_vals[lc_num], u->in_shifts_filename[vv], u->varnames[vv]);
	  exit(1);
	}
      }
      continue;
    }
    if(fieldlabels_indx == NULL) {
      if((fieldlabels_indx = (int *) malloc(p->NJD[lc_num]*sizeof(int))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
    }
    unstitch_Find_FieldLabels_in_InshiftFile(u, p, lc_num, lc_name_num, vv, fieldlabels_indx, inshift_starid);
    unstitch_Apply_InShifts(u, p, lc_num, lc_name_num, vv, fieldlabels_indx, inshift_starid);
  }

  if(fieldlabels_indx != NULL) free(fieldlabels_indx);
}

/* "fitsheader" source: read the per-segment shifts from the input FITS header
   keywords that -stitch wrote with its add_shifts_fitsheader option, and add
   them back.  The keyword comments ("Shift for variable V LCgroup L [REFID R]")
   are the authoritative key -- the keyword character encoding is not uniquely
   invertible for the refnum case, so we match on the parsed comment instead.
   The header stores -shiftvalue, so we SUBTRACT the keyword value to restore
   the original magnitude. */
void DoUnstitch_FitsHeader(ProgramData *p, _Unstitch *u, int lc_name_num, int lc_num)
{
#ifdef HAVE_CFITSIO
  fitsfile *fptr = NULL;
  int status = 0;
  int nkeys, kk, hdutype;
  int vv, i, baselen;
  char keyname[FLEN_KEYWORD], value[FLEN_VALUE], comment[FLEN_COMMENT];
  /* table of shift records parsed from the header */
  int *map_V = NULL, *map_L = NULL, *map_R = NULL, *map_hasR = NULL;
  double *map_val = NULL;
  int nmap = 0, sizemap = 0;
  int V, L, R, Rq, hasR, np, matched;
  int *matchidx = NULL;
  int miss_found, miss_L0, miss_R0, has_second, second_L, second_R;

  u->Npoints_shifted[lc_num] = 0;
  if(p->NJD[lc_num] <= 0)
    return;

  baselen = strlen(u->keywordbase);

  fits_open_file(&fptr, p->lcnames[lc_name_num], READONLY, &status);
  if(status) {
    fprintf(stderr,"Error in -unstitch: cannot open the input light curve \"%s\" as a FITS file to read shift keywords (cfitsio status %d). The \"fitsheader\" source requires a FITS-format input light curve.\n", p->lcnames[lc_name_num], status);
    exit(1);
  }
  fits_movabs_hdu(fptr, u->fits_hdu, &hdutype, &status);
  if(status) {
    fprintf(stderr,"Error in -unstitch: cannot move to HDU %d of the input FITS file \"%s\" (cfitsio status %d).\n", u->fits_hdu, p->lcnames[lc_name_num], status);
    exit(1);
  }

  /* Collect every keyword whose name starts with keywordbase and whose comment
     is a -stitch shift record. */
  fits_get_hdrspace(fptr, &nkeys, NULL, &status);
  for(kk = 1; kk <= nkeys; kk++) {
    keyname[0] = '\0'; value[0] = '\0'; comment[0] = '\0';
    if(fits_read_keyn(fptr, kk, keyname, value, comment, &status)) { status = 0; continue; }
    if(baselen > 0 && strncmp(keyname, u->keywordbase, baselen) != 0) continue;
    np = sscanf(comment, "Shift for variable %d LCgroup %d REFID %d", &V, &L, &R);
    if(np == 3) {
      hasR = 1;
    } else {
      np = sscanf(comment, "Shift for variable %d LCgroup %d", &V, &L);
      if(np != 2) continue;   /* not a stitch shift keyword */
      hasR = 0;
      R = 0;
    }
    if(nmap >= sizemap) {
      sizemap = sizemap ? sizemap*2 : 16;
      if((map_V = (int *) realloc(map_V, sizemap*sizeof(int))) == NULL ||
	 (map_L = (int *) realloc(map_L, sizemap*sizeof(int))) == NULL ||
	 (map_R = (int *) realloc(map_R, sizemap*sizeof(int))) == NULL ||
	 (map_hasR = (int *) realloc(map_hasR, sizemap*sizeof(int))) == NULL ||
	 (map_val = (double *) realloc(map_val, sizemap*sizeof(double))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
    }
    map_V[nmap] = V; map_L[nmap] = L; map_R[nmap] = R;
    map_hasR[nmap] = hasR; map_val[nmap] = atof(value);
    nmap++;
  }
  status = 0;
  fits_close_file(fptr, &status);

  /* Apply, one variable at a time.  -stitch omits the reference segment's
     keyword (its shift is 0), so a single segment with no matching keyword is
     taken to be that reference and left unchanged.  If two or more distinct
     unmasked segments lack a keyword, the header does not cover the light
     curve (wrong keywordbase/HDU, or the LC was not stitched) and we quit. */
  if((matchidx = (int *) malloc(p->NJD[lc_num]*sizeof(int))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  for(vv = 0; vv < u->nvar; vv++) {
    miss_found = 0; has_second = 0;
    miss_L0 = 0; miss_R0 = 0; second_L = 0; second_R = 0;
    /* Pass 1: resolve each point and detect uncovered segments. */
    for(i = 0; i < p->NJD[lc_num]; i++) {
      L = u->fits_lcnumval[lc_num][i];
      Rq = u->fits_userefnum ? u->fits_refnumval[lc_num][i] : 0;
      matched = -1;
      for(kk = 0; kk < nmap; kk++) {
	if(map_V[kk] != vv || map_L[kk] != L) continue;
	if(u->fits_userefnum) {
	  if(!map_hasR[kk] || map_R[kk] != Rq) continue;
	} else {
	  if(map_hasR[kk]) continue;
	}
	matched = kk;
	break;
      }
      matchidx[i] = matched;
      if(matched < 0 && (!u->usemask || u->maskvals[lc_num][i] > VARTOOLS_MASK_TINY)) {
	if(!miss_found) {
	  miss_found = 1; miss_L0 = L; miss_R0 = Rq;
	} else if(!has_second && (L != miss_L0 || (u->fits_userefnum && Rq != miss_R0))) {
	  has_second = 1; second_L = L; second_R = Rq;
	}
      }
    }
    if(has_second) {
      fprintf(stderr,"Error in -unstitch: HDU %d of \"%s\" (keywordbase %s) is missing shift keywords for two or more distinct segments of variable %s (e.g. LCgroup %d%s and LCgroup %d%s); at most one segment (the reference, shift 0) may be absent. Check that the correct keywordbase, HDU, and lcnum/refnum variables are being used.\n",
	      u->fits_hdu, p->lcnames[lc_name_num], u->keywordbase, u->varnames[vv],
	      miss_L0, (u->fits_userefnum ? " (refnum given)" : ""),
	      second_L, (u->fits_userefnum ? " (refnum given)" : ""));
      exit(1);
    }
    /* Pass 2: apply the matched shifts (the single missing segment, if any,
       is the reference and is left unchanged). */
    for(i = 0; i < p->NJD[lc_num]; i++) {
      if(matchidx[i] >= 0) {
	if(u->noshiftmasked && u->usemask && u->maskvals[lc_num][i] <= VARTOOLS_MASK_TINY)
	  continue;   /* masked point left unchanged */
	u->varvals[vv][lc_num][i] -= map_val[matchidx[i]];
	u->Npoints_shifted[lc_num] += 1;
      }
    }
  }
  free(matchidx);

  if(map_V != NULL) free(map_V);
  if(map_L != NULL) free(map_L);
  if(map_R != NULL) free(map_R);
  if(map_hasR != NULL) free(map_hasR);
  if(map_val != NULL) free(map_val);
#endif
  return;
}
