#include "../../src/vartools.h"
#include <stdio.h>
#include <stdlib.h>
#include "stitch.h"

void DoStitch(ProgramData *p, _Stitch *stitch, int lc_name_num, 
			 int lc_num);

void stitch_ReadInshifts_File(ProgramData *p, Command *c, _Stitch *stitch, int vv);

/* This is the source code for a sample user-defined command
   to be used with vartools.

   This library defines the command -stitch which can be used
   to stitch light curves together
*/

void stitch_Initialize(char *commandname,
			  int *RequireReadAll,
			  int *RequireSortLC,
			  int *RequireDistinctTimes,
			  size_t *sizeuserdata)
/* This function defines five things which vartools
   will need to know about the command 

   Every library with the name $LIBNAME.so must contain a function with
   the name $LIBNAME_Initialize to set the five input variables.
*/
{
  /* Set the string used to call the command */
  sprintf(commandname,"-stitch");

  /* This is 1 if the command requires all light curves to be read
     at once, otherwise set it to 0. */
  *RequireReadAll = 0;

  /* This is 1 if the command requires input light curves to be sorted by
     time, otherwise set it to 0. */
  *RequireSortLC = 0;
  
  /* This is 1 if the command requires the times of observation in the
     input light curve to be unique */
  *RequireDistinctTimes = 0;

  /* You should define a structure to store the data needed for this
     command (in this example, a vector to hold the values to add to
     the magnitudes of each lc), below you would replace "_Stitch"
     with the type name of your structure.

     See magadd.h
*/
  *sizeuserdata = sizeof(_Stitch);
}

int stitch_ParseCL(ProgramData *p, Command *c,
			   void *userdata, int *iret, char **argv, int argc)
/* Every library with the name $LIBNAME.so must contain a function
   with the name $LIBNAME_ParseCL to check that the user has issued a
   valid syntax, and set-up the command-data accordingly.

   Note that any vectors/arrays used by the command should also be registered
   by this function so that memory will be allocated as needed.

Expected syntax for this command:
   -stitch stitch_variable_list uncertainty_variable_list mask_variable_list lcnum_var [\"refnum_var\" refnum_var] <\"median\" | \"mean\" | \"weightedmean\" | \"poly\" order | \"harmseries\" period_var Nharm> [\"groupbytime\" time_bin [\"start\" firstbintime]] [\"fitonly\"] [\"save_fitted_parameters\" <outdir [\"format\"] fmt>] [\"add_stitchparams_fitsheader\"  [\"primary\" | \"extension\"] [\"append\" | \"update\"]] [\"add_shifts_fitsheader\" keywordbase [\"primary\" | \"extension\"] [\"append\" | \"update\"]] [\"shifts_file\" fieldlabelsvar starnamevar [\"append_refnum_to_fieldlabel\"] [\"in_shifts_file\" inshiftsfile1[,inshiftsfile2,...] [\"nobs_refit\" nobs_refit] [\"header_basename_only\"]] [\"out_shifts_file\" outshiftsfile1[,outshiftsfile2,...] [\"include_missing\"]]]


  p = structure containing general program data. In most cases this 
      structure is only passed to other functions (like RegisterDataVector)
      and not used directly by the ParseCL function.

  c = pointer to a generic container for this command. This structure is
      needed by other functions that might be called here (RegisterDataVector)
      but should not be used directly by the ParseCL function.

  userdata = pointer to the structure storing data for this
      command. This should be cast to the appropriate type and the
      data there-in should be modified as needed by the ParseCL
      function. In this example, userdata corresponds to the _Magadd
      structure which stores the values to be added to the light
      curves.

  iret = command line argument index. This is initially zero, and on
      return should be set to the index of the last argument parsed by
      this command. Note that index 0 points to the name of the
      command, which does not need to be checked (this function is
      only called if the user issues the command). In this example, if
      the user issues the command "-magadd fix 5." the return value
      of iret should be 2.

  argv = array of command line arguments. argv[0] is the name of the
      command ("-magadd" in this case).

  argc = Number of command line arguments in the argv array (including
      0). It is the user's responsibility to check that i < argc
      before attempting to parse argv[i], failure to do so may lead to
      segmentation violations.
 */
{
  int i = 0;
  int j;
  int i1, i2, k;
  int check;

  /* Cast the input userdata pointer to the data structure type used
     for this command */
  _Stitch *stitch;
  stitch = (_Stitch *) userdata;

  /* We can use the function ParseFixSpecFixcolumn to parse a line of the form
     <"fix" value | "list" [\"column\" col] | "fixcolumn" <colname | colnum> | \"expr\" expr>

     After giving p, c, iret, argv, and argc, give the number of
     data-vectors to be set by this command. For each data-vector give
     the data-type for the vector and a pointer to the vector.

     The function will update iret.
  */
  i = *iret;
  if(i >= argc) {
    return 1;
  }
  
  /* Find the number of variables in the list */
  j = 0;
  stitch->stitchvarnames = NULL;
  stitch->nstitchvar = 1;
  while(argv[i][j] != '\0') {
    if(argv[i][j] == ',')
      stitch->nstitchvar += 1;
    j++;
  }

  /* Allocate memory for the variable names */
  if((stitch->stitchvarnames = (char **) malloc(stitch->nstitchvar * sizeof(char *))) == NULL ||
     (stitch->stitcherrnames = (char **) malloc(stitch->nstitchvar * sizeof(char *))) == NULL ||
     (stitch->stitchmasknames = (char **) malloc(stitch->nstitchvar * sizeof(char *))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);

  /* Parse the list of variables */
  i1 = 0;
  i2 = 0;
  for(k = 0; k < stitch->nstitchvar; k++) {
    while(argv[i][i2] != '\0' && argv[i][i2] != ',') {
      i2++;
    }
    if((stitch->stitchvarnames[k] = (char *) malloc((i2 - i1 + 1)*sizeof(char))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    for(j=i1; j < i2; j++) {
      stitch->stitchvarnames[k][j-i1] = argv[i][j];
    }
    stitch->stitchvarnames[k][j-i1] = '\0';
    i1 = i2+1;
    i2 = i2+1;
  }

  i++;
  if(i >= argc) {
    return 1;
  }

  /* Check whether the number of terms in the error list matches the
     number in the variable list */
  j = 0;
  k = 1;
  while(argv[i][j] != '\0') {
    if(argv[i][j] == ',')
      k += 1;
    j++;
  }
  
  if(k != stitch->nstitchvar) {
    fprintf(stderr,"Error parsing the stitch command - the number of variables in the uncertainty_variable_list doese not match the number in the stitch_variable_list.\n");
    return 1;
  }


  /* Parse the list of error variables */
  i1 = 0;
  i2 = 0;
  for(k = 0; k < stitch->nstitchvar; k++) {
    while(argv[i][i2] != '\0' && argv[i][i2] != ',') {
      i2++;
    }
    if((stitch->stitcherrnames[k] = (char *) malloc((i2 - i1 + 1)*sizeof(char))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    for(j=i1; j < i2; j++) {
      stitch->stitcherrnames[k][j-i1] = argv[i][j];
    }
    stitch->stitcherrnames[k][j-i1] = '\0';
    i1 = i2+1;
    i2 = i2+1;
  }

  i++;
  if(i >= argc) {
    return 1;
  }

  /* Check whether the number of terms in the mask list matches the
     number in the variable list */
  j = 0;
  k = 1;
  while(argv[i][j] != '\0') {
    if(argv[i][j] == ',')
      k += 1;
    j++;
  }
  
  if(k != stitch->nstitchvar) {
    fprintf(stderr,"Error parsing the stitch command - the number of variables in the mask_variable_list does not match the number in the stitch_variable_list.\n");
    return 1;
  }


  /* Parse the list of error variables */
  i1 = 0;
  i2 = 0;
  for(k = 0; k < stitch->nstitchvar; k++) {
    while(argv[i][i2] != '\0' && argv[i][i2] != ',') {
      i2++;
    }
    if((stitch->stitchmasknames[k] = (char *) malloc((i2 - i1 + 1)*sizeof(char))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    for(j=i1; j < i2; j++) {
      stitch->stitchmasknames[k][j-i1] = argv[i][j];
    }
    stitch->stitchmasknames[k][j-i1] = '\0';
    i1 = i2+1;
    i2 = i2+1;
  }

  /* Register all of the input variables */
  if((stitch->stitchvarvals = (double ***) malloc(stitch->nstitchvar * sizeof(double **))) == NULL ||
     (stitch->stitcherrvals = (double ***) malloc(stitch->nstitchvar * sizeof(double **))) == NULL ||
     (stitch->stitchmaskvals = (double ***) malloc(stitch->nstitchvar * sizeof(double **))) == NULL ||
     (stitch->stitchvars = (_Variable **) malloc(stitch->nstitchvar * sizeof(_Variable *))) == NULL ||
     (stitch->stitcherrvars = (_Variable **) malloc(stitch->nstitchvar * sizeof(_Variable *))) == NULL ||
     (stitch->stitchmaskvars = (_Variable **) malloc(stitch->nstitchvar * sizeof(_Variable *))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  
  for(k = 0; k < stitch->nstitchvar; k++) {
    VARTOOLS_RegisterDataVector(p, c, (void *) (&(stitch->stitchvarvals[k])),
				VARTOOLS_TYPE_DOUBLE, 0,
				VARTOOLS_SOURCE_EXISTINGVARIABLE,
				0, NULL, stitch->stitchvarnames[k],
				(char) VARTOOLS_VECTORTYPE_LC,
				&(stitch->stitchvars[k]));
    VARTOOLS_RegisterDataVector(p, c, (void *) (&(stitch->stitcherrvals[k])),
				VARTOOLS_TYPE_DOUBLE, 0,
				VARTOOLS_SOURCE_EXISTINGVARIABLE,
				0, NULL, stitch->stitcherrnames[k],
				(char) VARTOOLS_VECTORTYPE_LC,
				&(stitch->stitcherrvars[k]));
    VARTOOLS_RegisterDataVector(p, c, (void *) (&(stitch->stitchmaskvals[k])),
				VARTOOLS_TYPE_DOUBLE, 0,
				VARTOOLS_SOURCE_EXISTINGVARIABLE,
				0, NULL, stitch->stitchmasknames[k],
				(char) VARTOOLS_VECTORTYPE_LC,
				&(stitch->stitchmaskvars[k]));
  }

  i++;
  if(i >= argc) {
    return 1;
  }
  
  /* Parse the lcnum_var */
  VARTOOLS_RegisterDataVector(p, c, (void *) (&(stitch->lcnumval)),
			      VARTOOLS_TYPE_INT, 0,
			      VARTOOLS_SOURCE_EXISTINGVARIABLE,
			      0, NULL, argv[i],
			      (char) VARTOOLS_VECTORTYPE_LC,
			      &(stitch->lcnumvar));

  stitch->userefnum = 0;
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"refnum_var")) {
      stitch->userefnum = 1;
      i++;
      if(i >= argc) {
	return 1;
      }
      VARTOOLS_RegisterDataVector(p, c, (void *) (&(stitch->refnumval)),
				  VARTOOLS_TYPE_INT, 0,
				  VARTOOLS_SOURCE_EXISTINGVARIABLE,
				  0, NULL, argv[i],
				  (char) VARTOOLS_VECTORTYPE_LC,
				  &(stitch->refnumvar));

      VARTOOLS_RegisterDataVector(p, c, (void *) (&stitch->minrefnumindx),
				  VARTOOLS_TYPE_INT, 0, VARTOOLS_SOURCE_COMPUTED,
				  0, NULL);
      VARTOOLS_RegisterDataVector(p, c, (void *) (&stitch->maxrefnumindx),
				  VARTOOLS_TYPE_INT, 0, VARTOOLS_SOURCE_COMPUTED,
				  0, NULL);

    } else
      i--;
  } else
    i--;

  i++;
  if(i >= argc) {
    return 1;
  }
  if(!strcmp(argv[i],"median")) {
    stitch->stitchmethod = VARTOOLS_STITCH_METHOD_MEDIAN;
  }
  else if(!strcmp(argv[i],"mean")) {
    stitch->stitchmethod = VARTOOLS_STITCH_METHOD_MEAN;
  }
  else if(!strcmp(argv[i],"weightedmean")) {
    stitch->stitchmethod = VARTOOLS_STITCH_METHOD_WEIGHTEDMEAN;
  }
  else if(!strcmp(argv[i],"poly")) {
    stitch->stitchmethod = VARTOOLS_STITCH_METHOD_POLY;

    i++;
    if(i >= argc) {
      return 1;
    }
    
    stitch->polyorder = atoi(argv[i]);
  }
  else if(!strcmp(argv[i],"harmseries")) {
    stitch->stitchmethod = VARTOOLS_STITCH_METHOD_HARM;
    
    i++;
    if(i >= argc) {
      return 1;
    }

    VARTOOLS_RegisterDataVector(p, c, (void *) (&(stitch->harmperiodvals)),
				VARTOOLS_TYPE_DOUBLE, 0,
				VARTOOLS_SOURCE_EXISTINGVARIABLE,
				0, NULL, argv[i],
				(char) VARTOOLS_VECTORTYPE_PERSTARDATA,
				&(stitch->harmperiodvar));
    
    i++;
    if(i >= argc) {
      return 1;
    }

    stitch->Nharm = atoi(argv[i]);
  }
  else
    return 1;

  stitch->groupbytime = 0;
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"groupbytime")) {
      stitch->groupbytime = 1;

      i++;
      if(i >= argc) return 1;
      stitch->time_step = atof(argv[i]);
      
      stitch->is_start_time = 0;
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"start")) {
	  stitch->is_start_time = 1;
	  
	  i++;
	  if(i >= argc) return 1;
	  stitch->start_time = atof(argv[i]);
	} else
	  i--;
      } else
	i--;
    } else
      i--;
  } else
    i--;

  stitch->fitonly = 0;
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"fitonly")) {
      stitch->fitonly = 1;
    } else
      i--;
  } else
    i--;

  stitch->save_fitted_parameters = 0;
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"save_fitted_parameters")) {
      stitch->save_fitted_parameters = 1;
      
      i++;
      if(i >= argc) return 1;
      sprintf(stitch->coeffoutdir,"%s", argv[i]);
      
      stitch->coeffoutformat[0] = '\0';
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"format")) {
	  i++;
	  if(i >= argc) return 1;
	  sprintf(stitch->coeffoutformat,"%s",argv[i]);
	} else
	  i--;
      } else
	i--;
    } else
      i--;
  } else
    i--;
  
  stitch->add_stitchparams_fitsheader = 0;
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"add_stitchparams_fitsheader")) {
      stitch->add_stitchparams_fitsheader = 1;
            
      stitch->add_stitchparams_fitsheader_hdutouse = 0;
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"primary")) {
	  stitch->add_stitchparams_fitsheader_hdutouse = 0;
	}
	else if(!strcmp(argv[i],"extension")) {
	  stitch->add_stitchparams_fitsheader_hdutouse = 1;
	}
	else
	  i--;
      } else
	i--;
      
      stitch->add_stitchparams_fitsheader_updateexistingkeyword = 1;
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"append")) {
	  stitch->add_stitchparams_fitsheader_updateexistingkeyword = 0;
	}
	else if(!strcmp(argv[i],"update")) {
	  stitch->add_stitchparams_fitsheader_updateexistingkeyword = 1;
	}
	else
	  i--;
      } else
	i--;
    } else
      i--;
  } else
    i--;


  stitch->add_shifts_fitsheader = 0;
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"add_shifts_fitsheader")) {
      stitch->add_shifts_fitsheader = 1;
      
      i++;
      if(i >= argc) return 1;
      sprintf(stitch->keywordbase,"%s",argv[i]);
      
      stitch->add_shifts_fitsheader_hdutouse = 0;
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"primary")) {
	  stitch->add_shifts_fitsheader_hdutouse = 0;
	}
	else if(!strcmp(argv[i],"extension")) {
	  stitch->add_shifts_fitsheader_hdutouse = 1;
	}
	else
	  i--;
      } else
	i--;
      
      stitch->add_shifts_fitsheader_updateexistingkeyword = 1;
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"append")) {
	  stitch->add_shifts_fitsheader_updateexistingkeyword = 0;
	}
	else if(!strcmp(argv[i],"update")) {
	  stitch->add_shifts_fitsheader_updateexistingkeyword = 1;
	}
	else
	  i--;
      } else
	i--;
    } else
      i--;
  } else
    i--;

  stitch->is_in_shifts_file = 0;
  stitch->is_out_shifts_file = 0;
  stitch->is_nobs_refit = 0;
  stitch->is_append_refnum_to_fieldlabel = 0;
  stitch->include_missing_inputstars = 0;
  stitch->has_output_missing = 0;
  stitch->is_inshifts_header_basename_only = 0;

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"shifts_file")) {

      i++;
      if(i >= argc) return 1;
      VARTOOLS_RegisterDataVector(p, c, (void *) (&(stitch->field_labels_vals)),
				  VARTOOLS_TYPE_STRING, 0,
				  VARTOOLS_SOURCE_EXISTINGVARIABLE,
				  0, NULL, argv[i],
				  (char) VARTOOLS_VECTORTYPE_LC,
				  &(stitch->field_labels_var));

      VARTOOLS_RegisterDataVector(p, c, (void *) (&(stitch->field_labels_vals_indx)),
				  VARTOOLS_TYPE_INT, 0,
				  VARTOOLS_SOURCE_LC,
				  0, NULL, NULL, -1, NULL, NULL);

      VARTOOLS_RegisterDataVector(p, c, (void *) (&(stitch->field_labels_vals_issorted)),
				  VARTOOLS_TYPE_INT, 0,
				  VARTOOLS_SOURCE_EVALEXPRESSION,
				  0, NULL, "0");

      
      i++;
      if(i >= argc) return 1;
      VARTOOLS_RegisterDataVector(p, c, (void *) (&(stitch->starname_vals)),
				  VARTOOLS_TYPE_STRING, 0,
				  VARTOOLS_SOURCE_EXISTINGVARIABLE,
				  0, NULL, argv[i],
				  (char) VARTOOLS_VECTORTYPE_PERSTARDATA,
				  &(stitch->starname_var));

      if((stitch->starname_varname = (char *) malloc((strlen(argv[i])+1)*sizeof(char))) == NULL) 
	VARTOOLS_error(ERR_MEMALLOC);
      sprintf(stitch->starname_varname,"%s",argv[i]);
      

      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"append_refnum_to_fieldlabel")) {
	  stitch->is_append_refnum_to_fieldlabel = 1;
	}
	else
	  i--;
      } else
	i--;

      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"in_shifts_file")) {
	  stitch->is_in_shifts_file = 1;
	  
	  i++;
	  if(i >= argc) return 1;
	  
	  
	  if((stitch->in_shifts_filename = (char **) malloc(stitch->nstitchvar*sizeof(char *))) == NULL ||
	     (stitch->N_in_shift_stars = (int *) malloc(stitch->nstitchvar*sizeof(int))) == NULL ||
	     (stitch->N_shifts_per_star = (int **) malloc(stitch->nstitchvar*sizeof(int *))) == NULL ||
	     (stitch->size_in_shifts_file = (int *) malloc(stitch->nstitchvar*sizeof(int))) == NULL ||
	     (stitch->in_shift_starnames = (char ***) malloc(stitch->nstitchvar*sizeof(char **))) == NULL ||
	     (stitch->in_shift_stars_found = (int **) malloc(stitch->nstitchvar*sizeof(int *))) == NULL ||
	     (stitch->in_shift_starnames_sortidx = (int **) malloc(stitch->nstitchvar*sizeof(int *))) == NULL ||
	     (stitch->in_shift_labels = (char ****) malloc(stitch->nstitchvar*sizeof(char ***))) == NULL ||
	     (stitch->in_shift_labels_sortidx = (int ***) malloc(stitch->nstitchvar*sizeof(int **))) == NULL ||
	     (stitch->in_shift_values = (double ***) malloc(stitch->nstitchvar*sizeof(double **))) == NULL ||
	     (stitch->in_shift_values_str = (char ****) malloc(stitch->nstitchvar*sizeof(char ***))) == NULL ||
	     (stitch->Nobs_in_shifts = (int ***) malloc(stitch->nstitchvar*sizeof(int **))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);

	  /* Check whether the number of terms in the in_shifts_file list 
	     matches the
	     number in the variable list */
	  j = 0;
	  k = 1;
	  while(argv[i][j] != '\0') {
	    if(argv[i][j] == ',')
	      k += 1;
	    j++;
	  }
	  
	  if(k != stitch->nstitchvar) {
	    fprintf(stderr,"Error parsing the stitch command - the number of files in the in_shifts_file list does not match the number in the stitch_variable_list.\n");
	    return 1;
	  }
	  
	  
	  /* Parse the list of in shift files */
	  i1 = 0;
	  i2 = 0;
	  for(k = 0; k < stitch->nstitchvar; k++) {
	    while(argv[i][i2] != '\0' && argv[i][i2] != ',') {
	      i2++;
	    }
	    if((stitch->in_shifts_filename[k] = (char *) malloc((i2 - i1 + 1)*sizeof(char))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=i1; j < i2; j++) {
	      stitch->in_shifts_filename[k][j-i1] = argv[i][j];
	    }
	    stitch->in_shifts_filename[k][j-i1] = '\0';
	    i1 = i2+1;
	    i2 = i2+1;
	  
	    stitch_ReadInshifts_File(p, c, stitch, k);
	  }

	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"nobs_refit")) {
	      stitch->is_nobs_refit = 1;
	      i++;
	      if(i >= argc) return(1);
	      stitch->nobs_refit = atoi(argv[i]);
	    } else
	      i--;
	  } else
	    i--;

	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"header_basename_only")) {
	      stitch->is_inshifts_header_basename_only = 1;
	    } else
	      i--;
	  } else
	    i--;
	} else
	  i--;
      } else
	i--;

      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"out_shifts_file")) {
	  stitch->is_out_shifts_file = 1;

	  if((stitch->out_shifts_filename = (char **) malloc(stitch->nstitchvar*sizeof(char *))) == NULL ||
	     (stitch->out_shifts_file = (FILE **) malloc(stitch->nstitchvar*sizeof(FILE *))) == NULL ||
	     (stitch->N_shifts_per_star_out = (int **) malloc(stitch->nstitchvar*sizeof(int *))) == NULL ||
	     (stitch->out_shift_labels = (char ****) malloc(stitch->nstitchvar*sizeof(char ***))) == NULL ||
	     (stitch->out_shift_values = (double ***) malloc(stitch->nstitchvar*sizeof(double **))) == NULL ||
	     (stitch->Nobs_out_shifts = (int ***) malloc(stitch->nstitchvar*sizeof(int **))) == NULL) {
	    VARTOOLS_error(ERR_MEMALLOC);
	  }

	  i++;
	  if(i >= argc) return 1;


	  /* Check whether the number of terms in the out_shifts_file list 
	     matches the
	     number in the variable list */
	  j = 0;
	  k = 1;
	  while(argv[i][j] != '\0') {
	    if(argv[i][j] == ',')
	      k += 1;
	    j++;
	  }
	  
	  if(k != stitch->nstitchvar) {
	    fprintf(stderr,"Error parsing the stitch command - the number of files in the out_shifts_file list does not match the number in the stitch_variable_list.\n");
	    return 1;
	  }
	  
	  
	  /* Parse the list of out shift files */
	  i1 = 0;
	  i2 = 0;
	  for(k = 0; k < stitch->nstitchvar; k++) {
	    while(argv[i][i2] != '\0' && argv[i][i2] != ',') {
	      i2++;
	    }
	    if((stitch->out_shifts_filename[k] = (char *) malloc((i2 - i1 + 1)*sizeof(char))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=i1; j < i2; j++) {
	      stitch->out_shifts_filename[k][j-i1] = argv[i][j];
	    }
	    stitch->out_shifts_filename[k][j-i1] = '\0';
	    i1 = i2+1;
	    i2 = i2+1;
	  
	    if((stitch->out_shifts_file[k] = fopen(stitch->out_shifts_filename[k],"w")) == NULL) {
	      fprintf(stderr,"Cannot write the output stitch shifts to the file %s\n", stitch->out_shifts_filename[k]);
	      exit(1);
	    }
	    VARTOOLS_RegisterTrackedOpenFile(p,stitch->out_shifts_file[k]);
	  }

	  i++;
	  if(i < argc) {
	    if(!strcmp(argv[i],"include_missing")) {
	      stitch->include_missing_inputstars = 1;
	    } else
	      i--;
	  } else
	    i--;

	} else
	  i--;
      } else
	i--;
    } else
      i--;
  } else
    i--;

  /* Register the variables to store data that will appear in the output table */
  VARTOOLS_RegisterDataVector(p, c, (void *) (&stitch->Nlcgroups_used),
			      VARTOOLS_TYPE_INT, 0, VARTOOLS_SOURCE_COMPUTED,
			      1, "NLCGroups");
  VARTOOLS_RegisterDataVector(p, c, (void *) (&stitch->Ntimegroups_used),
			      VARTOOLS_TYPE_INT, 0, VARTOOLS_SOURCE_COMPUTED,
			      1, "NTimeGroups");
  VARTOOLS_RegisterDataVector(p, c, (void *) (&stitch->Nparamtotal_used),
			      VARTOOLS_TYPE_INT, 0, VARTOOLS_SOURCE_COMPUTED,
			      1, "NFitParamsTotal");

  if(stitch->groupbytime) {
    VARTOOLS_RegisterDataVector(p, c, (void *) (&stitch->final_time_step),
				VARTOOLS_TYPE_DOUBLE, 0, VARTOOLS_SOURCE_COMPUTED,
				1, "TimeStep");
  }

  i++;

  *iret = i;
  return 0;
}

void stitch_ShowSyntax(FILE *outfile)
/* Write out the expected syntax for this command to the file outfile */
{
  OutText s;
  s.s = NULL;
  s.space = 0;
  s.len_s = 0;
  s.Nchar_cur_line = 0;

  VARTOOLS_printtostring(&s,
		"-stitch\n");
  VARTOOLS_printtostring(&s,
		"\tstitch_variable_list uncertainty_variable_list\n");
  VARTOOLS_printtostring(&s,
		"\tmask_variable_list lcnum_var\n");
  VARTOOLS_printtostring(&s,
		"\t[\"refnum_var\" refnum_var]\n");
  VARTOOLS_printtostring(&s,
			 "\t<\"median\" | \"mean\" | \"weightedmean\"\n");
  VARTOOLS_printtostring(&s,
			 "\t\t| \"poly\" order | \"harmseries\" period_var Nharm>\n");
  VARTOOLS_printtostring(&s,
			 "\t[\"groupbytime\" time_bin [\"start\" firstbintime]]\n");
  VARTOOLS_printtostring(&s,
			 "\t[\"fitonly\"]\n");
  VARTOOLS_printtostring(&s,
			 "\t[\"save_fitted_parameters\" <outdir [\"format\"] fmt>]\n");
  VARTOOLS_printtostring(&s,
			 "\t[\"add_stitchparams_fitsheader\" [\"primary\" | \"extension\"]\n");
  VARTOOLS_printtostring(&s,
			 "\t\t[\"append\" | \"update\"]]\n");
  VARTOOLS_printtostring(&s,
			 "\t[\"add_shifts_fitsheader\" keywordbase [\"primary\" | \"extension\"]\n");
  VARTOOLS_printtostring(&s,
			 "\t\t[\"append\" | \"update\"]]\n");
  VARTOOLS_printtostring(&s,
			 "\t[\"shifts_file\" fieldlabelsvar starnamevar\n");
  VARTOOLS_printtostring(&s,
			 "\t\t[\"append_refnum_to_fieldlabel\"]\n");
  VARTOOLS_printtostring(&s,
			 "\t\t[\"in_shifts_file\" inshiftsfile1[,inshiftsfile2,...]\n");
  VARTOOLS_printtostring(&s,
			 "\t\t\t[\"nobs_refit\" nobs_refit] [\"header_basename_only\"]]\n");
  VARTOOLS_printtostring(&s,
			 "\t\t[\"out_shifts_file\" outshiftsfile1[,outshiftsfile2,...]\n");
  VARTOOLS_printtostring(&s,
			 "\t\t\t[\"include_missing\"]]]\n");
  fprintf(outfile,s.s);
}

void stitch_ShowHelp(FILE *outfile)
/* Output the help information for this command */
{
  /* Give the verbose help description */
  OutText s;
  s.s = NULL;
  s.space = 0;
  s.len_s = 0;
  s.Nchar_cur_line = 0;

  VARTOOLS_printtostring(&s,"Fit for, and remove, offsets between different light curve segments for a given source. This routine is meant to be used together with the \"combinelcs\" option to the -i or -l input flags.\n\n");
  VARTOOLS_printtostring(&s,"stitch_variable_list - A list of light curve variables to apply the stitch procedure to. This would typically be mag, but could be a comma-separated list of magnitude variables if using multiple apertures, and/or multiple detrending methods.\n\n");
  VARTOOLS_printtostring(&s,"uncertainty_variable_list - A comma-separated list of variables storing the uncertainties associated with the stitch_variable_list variables. This list must have the same number of elements as the prior list. In the typical use this would be err.\n\n");
  VARTOOLS_printtostring(&s,"mask_variable_list - A comma-separated list of variables to use for masking points from the fitting procedure. This list must have the same number of elements as the stitch_variable_list. Points with a mask value greater than 0 will be excluded from the fit.\n\n");
  VARTOOLS_printtostring(&s,"lcnum_var - A variable that identifies the light curve segment that each point is associated with. This variable might be one set using the \"lcnumvar\" keyword option with the \"combinelcs\" option to the -i command.\n\n");
  VARTOOLS_printtostring(&s,"\"refnum_var\" - Optionally include an additional \"reference\" number variable to further split the light curve into segments. Points with unique combinations of refnum_var and lcnum_var will be grouped into segments for the fitting.\n\n");
  VARTOOLS_printtostring(&s,"Following the variables, the method used for modelling the shifts between light curve segments must be specified. The options are as follows:\n\n");
  VARTOOLS_printtostring(&s,"\t\"median\" - take the median of each segment in each time group (see below for a discussion of optional time-grouping), and fit for shifts that minimize the differences between these medians.\n\n");
  VARTOOLS_printtostring(&s,"\t\"mean\" - as for the median, but use the mean of each segment.\n\n");
  VARTOOLS_printtostring(&s,"\t\"weightedmean\" - use the mean weighted with the inverse of the uncertainty squared.\n\n");
  VARTOOLS_printtostring(&s,"\t\"poly\" order - fit a polynomial in time, of specified order, to the light curve segment.\n\n");
  VARTOOLS_printtostring(&s,"\t\"harmseries\" period_var Nharm - fit a harmonic series in time, with a period given by the variable period_var and a specified number of harmonics (0 to fit a simple sinusoid).\n\n");
  VARTOOLS_printtostring(&s,"\"groupbytime\" - Give this keyword to group the light curve segments into time bins, the shift model will assume a shared intrinsic shape for the light curve (constant if using the median, mean, or weightedmean methods, a polynomial in time if using the poly method, or a harmonic series if using the harmseries method) within a given time-bin, and overall additive shifts will be fit to match the different light curve segments within their respective bins (the shifts are global across all bins, while the intrinsic shape is independent within each bin). Give the time_bin to use for binning the segments in the same units of time as used in the light curve. To set the start time for the first time bin, use the optional \"start\" keyword, followed by the start time to use. Note that the routine will make sure that all segments can be tied together, and the bin-size will be increased until that is true.  For example, if one is stitching two light curve segments together, and they are disjoint in time such that each time-bin contains points from only one light curve segment, then the time bin-size will be increased until at least one bin contains points from both segments.\n\n");
  VARTOOLS_printtostring(&s,"\"fitonly\" - fit for the shifts between light curve segments, but do not remove them.\n\n");
  VARTOOLS_printtostring(&s,"\"save_fitted_parameters\" - save the best fit parameters to a file. A separate file is used for each source in the input list. Specify the directory to save the files to. By default the coefficient files will have the same name as the input light curve (stripped of any leading directory names), with the \".stitch\" suffix appended. Use the \"format\" keyword to provide a rule for naming the files, where instances of \"%%s\" will be replaced with the input filename.\n\n");
  VARTOOLS_printtostring(&s,"\"add_stitchparams_fitsheader\" - log the vartools control parameters for stitching to the header of any light curve that is subsequently output in FITS format. Use the \"primary\" keyword to log the parameters to the primary header of the FITS file, or the \"extension\" keyword to log them in the first extension header. Use \"append\" to append the FITS header keywords to the header whether or not the keywords may already be present, or \"update\" to update any header keywords that may already be present.\n\n");
  VARTOOLS_printtostring(&s,"\"add_shifts_fitsheader\" - log the shifts determined by the -stitch command into the header of any light curve file that is subsequently output in FITS format. Here the basename of the keyword(s) used to store the shift values should be given. A good option is \"SHFT\". The values of lcnumvar associated with each shift will be indicated using two letters (AA for 0, AB for 1, AC for 2, etc). This will be followed by two letters indicating the value of refnumvar if relevant. Additional options are similar to those for the \"add_stitchparams_fitsheader\" option.\n\n");
  VARTOOLS_printtostring(&s,"\"shifts_file\" - optionally read-in from a file and apply any previously determined shifts for each star and light curve group, and/or write out the determined shifts for each star. First provide the name of a light curve variable storing a string field identifier for each observation. This field (or telescope ID) will be stored in the shifts file and used to match up points in a light curve with any previously measured shift for that field. Then provide the name of a input list variable storing the string star name for each light curve. This will also be stored in the shifts file and used to match up light curves with previously measured shifts for that star.  Optionally provide the \"append_refnum_to_fieldlabel\" to append the reference number to the end of the string field labels to provide a final string for matching.  To use an input shifts file give the \"in_shifts_file\" followed by a comma-separated list of shifts files to read-in previously measured shifts from. There should be one file per variable being shifted. Optionally give the \"nobs_refit\" followed by an integer to redetermine the shifts for any fields where the existing shift was based on fewer than nobs_refit observations, and where the new light curve contains at least that many observations of that field. The names of any input shift files used to determine shifts will be added to the header of subsequently output light curves if they are output in FITS format. Optionally give the \"header_basename_only\" to only list the basename of the shifts file, and not the full path.  Finally give the \"out_shifts_file\" to output a file with the measured shifts. Provide a list of files to use, one per variable being stitched. To include in the output shifts file any fields from the input list that had no observations in the current light curve, give the \"include_missing\" keyword.\n\n");
  fprintf(outfile,s.s);
}

void stitch_RunCommand(ProgramData *p, void *userdata, int lc_name_num, int lc_num)
/* This function runs the command on a light curve.  

   p = structure containing various general program data (in
       particular the light curves are contained in this structure).

   userdata = pointer to the structure containing the command specific
       data (including control parameters and vectors to store output
       results).

   lc_name_num = this is the index to use to access the light curve name.

   lc_num = this is the index to use to access the light curve and
   data from Registered vectors (e.g. the addval parameter for this
   command).

*/
{
  _Stitch *stitch;


  /* Translate pointers and structures input by vartools into
     easier to use forms */
  stitch = (_Stitch *) userdata;


  /* The number of points in the light curve, and vectors storing the
     time, magnitudes, and errors are accessed as illustrated
     below. Note that t, mag and err are vectors, t[0] is the time of
     the first observation in the light curve, t[NJD-1] is the time of
     the last observation, etc. */

  DoStitch(p, stitch, lc_name_num, lc_num);

}

int stitch_CloseCommand(ProgramData *p, void *userdata)
/* If a function named this way exists it will be called after all
   light curves have been processed, and can be used to perform
   any necessary actions at the end of the process.

   p = structure containing various general program data (in
       particular the light curves are contained in this structure).

   userdata = pointer to the structure containing the command specific
       data (including control parameters and vectors to store output
       results).

*/
{
  int i, j, vv;
  _Stitch *stitch;

  /* Translate pointers and structures input by vartools into
     easier to use forms */
  stitch = (_Stitch *) userdata;

#ifdef PARALLEL
  if(p->Nthread > 1) {
    while(pthread_mutex_trylock(&(p->outfile_mutex)));
  }
#endif

  if(stitch->is_out_shifts_file && stitch->include_missing_inputstars && 
     stitch->is_in_shifts_file)
    {
      for(vv=0; vv < stitch->nstitchvar; vv++) {
	for(i=0; i < stitch->N_in_shift_stars[vv]; i++) {
	  if(!stitch->in_shift_stars_found[vv][i]) {
	    fprintf(stitch->out_shifts_file[vv],"%s ",stitch->in_shift_starnames[vv][i]);
	    for(j=0; j < stitch->N_shifts_per_star[vv][i]; j++) {
	      if(j > 0) fprintf(stitch->out_shifts_file[vv],";");
	      fprintf(stitch->out_shifts_file[vv],"%s,%s,%d",stitch->in_shift_labels[vv][i][j],stitch->in_shift_values_str[vv][i][j],stitch->Nobs_in_shifts[vv][i][j]);
	    }
	    fprintf(stitch->out_shifts_file[vv],"\n");
	  }
	}
	fflush(stitch->out_shifts_file[vv]);
      }
    }

#ifdef PARALLEL
  if(p->Nthread > 1) {
    pthread_mutex_unlock(&(p->outfile_mutex));
  }
#endif

}

int stitch_ParseInShifts_Line(char *line, char **incols, int *size_incols)
{
  int j, i, k;
  j = 0;
  i = 0;
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

void stitch_ReadInshifts_File(ProgramData *p, Command *c, _Stitch *stitch, int vv)
{
  FILE *shiftfile;
  
  char *inputline = NULL;
  size_t size_inputline = 0;

  char **incols;
  int *size_incols;
  int i, j, k, mm, ll, jj, testskip;

  if((shiftfile = fopen(stitch->in_shifts_filename[vv],"r")) == NULL) {
    fprintf(stderr,"Cannot open the input shifts file: %s\n", stitch->in_shifts_filename[vv]);
    exit(1);
  }
  stitch->N_in_shift_stars[vv] = 0;
  stitch->size_in_shifts_file[vv] = 256;
  if((stitch->N_shifts_per_star[vv] = (int *) malloc(stitch->size_in_shifts_file[vv]*sizeof(int))) == NULL ||
     (stitch->in_shift_starnames[vv] = (char **) malloc(stitch->size_in_shifts_file[vv]*sizeof(char *))) == NULL ||
     (stitch->in_shift_stars_found[vv] = (int *) malloc(stitch->size_in_shifts_file[vv]*sizeof(int))) == NULL ||
     (stitch->in_shift_starnames_sortidx[vv] = (int *) malloc(stitch->size_in_shifts_file[vv]*sizeof(int))) == NULL ||
     (stitch->in_shift_labels[vv] = (char ***) malloc(stitch->size_in_shifts_file[vv]*sizeof(char **))) == NULL ||
     (stitch->in_shift_labels_sortidx[vv] = (int **) malloc(stitch->size_in_shifts_file[vv]*sizeof(int *))) == NULL ||
     (stitch->in_shift_values[vv] = (double **) malloc(stitch->size_in_shifts_file[vv]*sizeof(double *))) == NULL ||
     (stitch->in_shift_values_str[vv] = (char ***) malloc(stitch->size_in_shifts_file[vv]*sizeof(char **))) == NULL ||
     (stitch->Nobs_in_shifts[vv] = (int **) malloc(stitch->size_in_shifts_file[vv]*sizeof(int *))) == NULL ||
     (incols = (char **) malloc(2*sizeof(char *))) == NULL ||
     (size_incols = (int *) malloc(2*sizeof(int))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  incols[0] = NULL; incols[1] = NULL;
  size_incols[0] = 0; size_incols[1] = 0;
  i = 0;
  while(VARTOOLS_gnu_getline(&(inputline),&(size_inputline), shiftfile) >= 0) {
    if(inputline[0] != '#') {
      stitch->N_in_shift_stars[vv] += 1;
      if(stitch->N_in_shift_stars[vv] >= stitch->size_in_shifts_file[vv]) {
	stitch->size_in_shifts_file[vv] *= 2;
	if((stitch->N_shifts_per_star[vv] = (int *) realloc(stitch->N_shifts_per_star[vv], stitch->size_in_shifts_file[vv]*sizeof(int))) == NULL ||
	   (stitch->in_shift_starnames[vv] = (char **) realloc(stitch->in_shift_starnames[vv], stitch->size_in_shifts_file[vv]*sizeof(char *))) == NULL ||
	   (stitch->in_shift_stars_found[vv] = (int *) realloc(stitch->in_shift_stars_found[vv], stitch->size_in_shifts_file[vv]*sizeof(int))) == NULL ||
	   (stitch->in_shift_starnames_sortidx[vv] = (int *) realloc(stitch->in_shift_starnames_sortidx[vv], stitch->size_in_shifts_file[vv]*sizeof(int))) == NULL ||
	   (stitch->in_shift_labels[vv] = (char ***) realloc(stitch->in_shift_labels[vv], stitch->size_in_shifts_file[vv]*sizeof(char **))) == NULL ||
	   (stitch->in_shift_labels_sortidx[vv] = (int **) realloc(stitch->in_shift_labels_sortidx[vv], stitch->size_in_shifts_file[vv]*sizeof(int *))) == NULL ||
	   (stitch->in_shift_values[vv] = (double **) realloc(stitch->in_shift_values[vv], stitch->size_in_shifts_file[vv]*sizeof(double *))) == NULL ||
	   (stitch->in_shift_values_str[vv] = (char ***) realloc(stitch->in_shift_values_str[vv], stitch->size_in_shifts_file[vv]*sizeof(char **))) == NULL ||
	   (stitch->Nobs_in_shifts[vv] = (int **) realloc(stitch->Nobs_in_shifts[vv], stitch->size_in_shifts_file[vv]*sizeof(int *))) == NULL)
	  VARTOOLS_error(ERR_MEMALLOC);
      }
      testskip = stitch_ParseInShifts_Line(inputline, incols, size_incols);
      if(testskip) {
	fprintf(stderr,"Error parsing input shifts file: %s\n", stitch->in_shifts_filename[vv]);
	fclose(shiftfile);
	exit(1);
      }
      if((stitch->in_shift_starnames[vv][i] = (char *) malloc((strlen(incols[0])+1)*sizeof(char))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
      sprintf(stitch->in_shift_starnames[vv][i],"%s",incols[0]);
      stitch->in_shift_starnames_sortidx[vv][i] = i;
      stitch->in_shift_stars_found[vv][i] = 0;
      j = 0;
      k = 0;
      mm = 0;
      while(incols[1][j] != '\0' && incols[1][j] != '\n') {
	if(incols[1][j] == ',') mm++;
	if(incols[1][j] == ';') {
	  if(mm != 2) {
	    fprintf(stderr,"Error parsing input shifts file: %s\n", stitch->in_shifts_filename[vv]);
	    fclose(shiftfile);
	    exit(1);
	  }
	  k++;
	  mm = 0;
	}
	j++;
      }
      if(mm != 2) {
	fprintf(stderr,"Error parsing input shifts file: %s\n", stitch->in_shifts_filename[vv]);
	fclose(shiftfile);
	exit(1);
      }
      stitch->N_shifts_per_star[vv][i] = k+1;
      if((stitch->in_shift_labels[vv][i] = (char **) malloc(stitch->N_shifts_per_star[vv][i]*sizeof(char *))) == NULL ||
	 (stitch->in_shift_labels_sortidx[vv][i] = (int *) malloc(stitch->N_shifts_per_star[vv][i]*sizeof(int))) == NULL ||
	 (stitch->in_shift_values[vv][i] = (double *) malloc(stitch->N_shifts_per_star[vv][i]*sizeof(double))) == NULL ||
	 (stitch->in_shift_values_str[vv][i] = (char **) malloc(stitch->N_shifts_per_star[vv][i]*sizeof(char *))) == NULL ||
	 (stitch->Nobs_in_shifts[vv][i] = (int *) malloc(stitch->N_shifts_per_star[vv][i]*sizeof(int))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
      j = 0;
      jj = 0;
      k = 0;
      mm = 0;
      while(incols[1][j] != '\0' && incols[1][j] != '\n') {
	if(incols[1][j] == ',') {
	  if(!mm) {
	    if((stitch->in_shift_labels[vv][i][k] = (char *) malloc((j - jj + 1)*sizeof(char))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(ll = jj; ll < j; ll++)
	      stitch->in_shift_labels[vv][i][k][ll-jj] = incols[1][ll];
	    stitch->in_shift_labels[vv][i][k][ll-jj] = '\0';
	    jj = j+1;
	  }
	  else if(mm == 1) {
	    incols[1][j] = '\0';
	    stitch->in_shift_values[vv][i][k] = atof(&(incols[1][jj]));
	    if((stitch->in_shift_values_str[vv][i][k] = (char *) malloc((j - jj + 1)*sizeof(char))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(ll = jj; ll < j; ll++)
	      stitch->in_shift_values_str[vv][i][k][ll-jj] = incols[1][ll];
	    stitch->in_shift_values_str[vv][i][k][ll-jj] = '\0';
	    jj = j+1;
	  }
	  mm++;
	}
	if(incols[1][j] == ';') {
	  incols[1][j] = '\0';
	  stitch->Nobs_in_shifts[vv][i][k] = atoi(&(incols[1][jj]));
	  k++;
	  mm = 0;
	  jj = j+1;
	}
	j++;
      }
      incols[1][j] = '\0';
      stitch->Nobs_in_shifts[vv][i][k] = atoi(&(incols[1][jj]));
      for(mm = 0; mm < stitch->N_shifts_per_star[vv][i]; mm++) {
	stitch->in_shift_labels_sortidx[vv][i][mm] = mm;
      }
      VARTOOLS_mysortstringint(stitch->N_shifts_per_star[vv][i], 0, stitch->in_shift_labels[vv][i], stitch->in_shift_labels_sortidx[vv][i]);

      i++;
    }
  }
  VARTOOLS_mysortstringint(stitch->N_in_shift_stars[vv], 0, stitch->in_shift_starnames[vv], stitch->in_shift_starnames_sortidx[vv]);
  fclose(shiftfile);
  if(incols[0] != NULL) free(incols[0]);
  if(incols[1] != NULL) free(incols[1]);
  free(incols);
  free(size_incols);
  if(inputline != NULL)
    free(inputline);
  return;
}

int stitch_Find_Star_in_InshiftFile(_Stitch *stitch, ProgramData *p, int lcnum, int lc_name_num, int vv)
{
  int i;
  if(!stitch->N_in_shift_stars[vv]) return -1;
  i = VARTOOLS_findX_string(stitch->in_shift_starnames[vv], stitch->in_shift_starnames_sortidx[vv], stitch->starname_vals[lcnum], 0, stitch->N_in_shift_stars[vv]);
  if(i >= 0) {
    i = stitch->in_shift_starnames_sortidx[vv][i];
    stitch->in_shift_stars_found[vv][i] = 1;
  }
  return i;
}

void stitch_Find_FieldLabels_in_InshiftFile(_Stitch *stitch, ProgramData *p, int lcnum, int lc_name_num, int vv, int *fieldlabels_indx, int shiftid)
{
  int i, j, jlastfound;
  if(!stitch->field_labels_vals_issorted[lcnum]) {
    for(i=0; i < p->NJD[lcnum]; i++) stitch->field_labels_vals_indx[lcnum][i] = i;
    VARTOOLS_mysortstringint(p->NJD[lcnum], 0, stitch->field_labels_vals[lcnum], stitch->field_labels_vals_indx[lcnum]);
    stitch->field_labels_vals_issorted[lcnum] = 1;
  }
  jlastfound = -1;
  for(i=0; i < p->NJD[lcnum]; i++) {
    if(jlastfound < 0) {
      j = VARTOOLS_findX_string(stitch->in_shift_labels[vv][shiftid], stitch->in_shift_labels_sortidx[vv][shiftid], stitch->field_labels_vals[lcnum][stitch->field_labels_vals_indx[lcnum][i]], 0, stitch->N_shifts_per_star[vv][shiftid]);
      if(j >= 0) {
	jlastfound = j;
	j = stitch->in_shift_labels_sortidx[vv][shiftid][j];
      }
    } else {
      j = VARTOOLS_findX_string(stitch->in_shift_labels[vv][shiftid], stitch->in_shift_labels_sortidx[vv][shiftid], stitch->field_labels_vals[lcnum][stitch->field_labels_vals_indx[lcnum][i]], jlastfound, stitch->N_shifts_per_star[vv][shiftid]);
      if(j >= 0) {
	jlastfound = j;
	j = stitch->in_shift_labels_sortidx[vv][shiftid][j];
      }
    }
    fieldlabels_indx[stitch->field_labels_vals_indx[lcnum][i]] = j;
  }
}

int stitch_Apply_InShifts(_Stitch *stitch, ProgramData *p, int lcnum, int lc_name_num, int vv, int *fieldlabels_indx, int shiftid)
{
  int retval = 1;
  int i;
  for(i=0; i < p->NJD[lcnum]; i++) {
    if(fieldlabels_indx[i] >= 0) {
      stitch->stitchvarvals[vv][lcnum][i] -= stitch->in_shift_values[vv][shiftid][fieldlabels_indx[i]];
    } else if(stitch->stitchmaskvals[vv][lcnum][i] > VARTOOLS_MASK_TINY)
      retval = 0;
  }
  return retval;
}

void stitch_WriteOutShifts_Line(_Stitch *stitch, ProgramData *p, int lcnum, int lc_name_num, int vv, int N_shift_out, int size_out_shift, char **out_shift_labels, int *size_out_shift_labels, char **out_shift_values_str, int *size_out_shift_values_str, int *Nobs_out_shift)
{
  int i;

#ifdef PARALLEL
  if(p->Nthread > 1) {
    while(pthread_mutex_trylock(&(p->outfile_mutex)));
  }
#endif

  /****** Write out the line to the file *****/
  if(N_shift_out > 0) {
    fprintf(stitch->out_shifts_file[vv],"%s ",stitch->starname_vals[lcnum]);
    for(i=0; i < N_shift_out; i++) {
      if(i > 0) fprintf(stitch->out_shifts_file[vv],";");
      fprintf(stitch->out_shifts_file[vv],"%s,%s,%d",out_shift_labels[i],out_shift_values_str[i],Nobs_out_shift[i]);
    }
    fprintf(stitch->out_shifts_file[vv],"\n");
  }
  fflush(stitch->out_shifts_file[vv]);
  
#ifdef PARALLEL
  if(p->Nthread > 1) {
    pthread_mutex_unlock(&(p->outfile_mutex));
  }
#endif
}

void PrepareOutShifts(ProgramData *p, _Stitch *stitch, int lc_name_num, int lc_num, int NLCgroups,  _StitchLightCurveGroup *lcg, int Ntimegroups, int Nusedtimegroups, _StitchTimeGroup *tg, int vv, int *fieldlabels_indx, int inshift_starid, int all_points_inshift, int *fieldlabels_needing_update, int *fieldlabels_obscounts, int *ret_N_shift_out, int *ret_size_out_shift, char ***ret_out_shift_labels, int **ret_size_out_shift_labels, char ***ret_out_shift_values_str, int **ret_size_out_shift_values_str, int **ret_Nobs_out_shift)
{
  int N_shift_out = 0;
  int size_out_shift = 0;
  char **out_shift_labels = NULL;
  int *size_out_shift_labels = NULL;
  char **out_shift_values_str = NULL;
  int *size_out_shift_values_str = NULL;
  int *Nobs_out_shift = NULL;
  int i, j, ll;

  N_shift_out = *ret_N_shift_out;
  size_out_shift = *ret_size_out_shift;
  out_shift_labels = *ret_out_shift_labels;
  size_out_shift_labels = *ret_size_out_shift_labels;
  out_shift_values_str = *ret_out_shift_values_str;
  size_out_shift_values_str = *ret_size_out_shift_values_str;
  Nobs_out_shift = *ret_Nobs_out_shift;
  if(all_points_inshift) {
    /* Copy over the terms from the inshift file */
    N_shift_out = stitch->N_shifts_per_star[vv][inshift_starid];
    if(!size_out_shift) {
      size_out_shift = stitch->N_shifts_per_star[vv][inshift_starid];
      N_shift_out = stitch->N_shifts_per_star[vv][inshift_starid];
      if((out_shift_labels = (char **) malloc(size_out_shift * sizeof(char *))) == NULL ||
	 (size_out_shift_labels = (int *) malloc(size_out_shift * sizeof(int))) == NULL ||
	 (out_shift_values_str = (char **) malloc(size_out_shift * sizeof(char *))) == NULL ||
	 (size_out_shift_values_str = (int *) malloc(size_out_shift * sizeof(int))) == NULL ||
	 (Nobs_out_shift = (int * ) malloc(size_out_shift * sizeof(int))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
      for(i = 0; i < N_shift_out; i++) {
	size_out_shift_labels[i] = 0;
	out_shift_labels[i] = NULL;
	size_out_shift_values_str[i] = 0;
	out_shift_values_str[i] = NULL;
      }
    } else if(stitch->N_shifts_per_star[vv][inshift_starid] > size_out_shift) {
      N_shift_out = stitch->N_shifts_per_star[vv][inshift_starid];
      if((out_shift_labels = (char **) realloc(out_shift_labels, N_shift_out * sizeof(char *))) == NULL ||
	 (size_out_shift_labels = (int *) realloc(size_out_shift_labels, N_shift_out * sizeof(int))) == NULL ||
	 (out_shift_values_str = (char **) realloc(out_shift_values_str, N_shift_out * sizeof(char *))) == NULL ||
	 (size_out_shift_values_str = (int *) realloc(size_out_shift_values_str, N_shift_out * sizeof(int))) == NULL ||
	 (Nobs_out_shift = (int * ) realloc(Nobs_out_shift, N_shift_out * sizeof(int))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
      for(i = size_out_shift; i < N_shift_out; i++) {
	size_out_shift_labels[i] = 0;
	out_shift_labels[i] = NULL;
	size_out_shift_values_str[i] = 0;
	out_shift_values_str[i] = NULL;
      }
      size_out_shift = N_shift_out;
    }
    for(i=0; i < N_shift_out; i++) {
      if(strlen(stitch->in_shift_labels[vv][inshift_starid][i])+1 > size_out_shift_labels[i]) {
	if(!size_out_shift_labels[i]) {
	  size_out_shift_labels[i] = strlen(stitch->in_shift_labels[vv][inshift_starid][i])+1;
	  if((out_shift_labels[i] = (char *) malloc(size_out_shift_labels[i]*sizeof(char))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	} else {
	  size_out_shift_labels[i] = strlen(stitch->in_shift_labels[vv][inshift_starid][i])+1;
	  if((out_shift_labels[i] = (char *) realloc(out_shift_labels[i], size_out_shift_labels[i]*sizeof(char))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	}
      }
      sprintf(out_shift_labels[i], "%s", stitch->in_shift_labels[vv][inshift_starid][i]);
      if(strlen(stitch->in_shift_values_str[vv][inshift_starid][i])+1 > size_out_shift_values_str[i]) {
	if(!size_out_shift_values_str[i]) {
	  size_out_shift_values_str[i] = strlen(stitch->in_shift_values_str[vv][inshift_starid][i])+1;
	  if((out_shift_values_str[i] = (char *) malloc(size_out_shift_values_str[i]*sizeof(char))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	} else {
	  size_out_shift_values_str[i] = strlen(stitch->in_shift_values_str[vv][inshift_starid][i])+1;
	  if((out_shift_values_str[i] = (char *) realloc(out_shift_values_str[i], size_out_shift_values_str[i]*sizeof(char))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	}
      }
      sprintf(out_shift_values_str[i], "%s", stitch->in_shift_values_str[vv][inshift_starid][i]);
      Nobs_out_shift[i] = stitch->Nobs_in_shifts[vv][inshift_starid][i];
    }
  } else {
    N_shift_out = 0;
    if(inshift_starid >= 0) {
      /* We need to copy terms from the in-shift file that were used */
      for(i=0; i < stitch->N_shifts_per_star[vv][inshift_starid]; i++) {
	if(!stitch->is_nobs_refit ||
	   (stitch->is_nobs_refit && (!fieldlabels_needing_update[i] || (fieldlabels_obscounts[i] < stitch->nobs_refit)))) {
	  if(N_shift_out + 1 > size_out_shift) {
	    if(!size_out_shift) {
	      size_out_shift = 1;
	      if((out_shift_labels = (char **) malloc(sizeof(char *))) == NULL ||
		 (size_out_shift_labels = (int *) malloc(sizeof(int))) == NULL ||
		 (out_shift_values_str = (char **) malloc(sizeof(char *))) == NULL ||
		 (size_out_shift_values_str = (int *) malloc(sizeof(int))) == NULL ||
		 (Nobs_out_shift = (int * ) malloc(sizeof(int))) == NULL)
		VARTOOLS_error(ERR_MEMALLOC);
	      for(j = 0; j < 1; j++) {
		size_out_shift_labels[j] = 0;
		out_shift_labels[j] = NULL;
		size_out_shift_values_str[j] = 0;
		out_shift_values_str[j] = NULL;
	      }
	    } else {
	      if((out_shift_labels = (char **) realloc(out_shift_labels, (N_shift_out + 1)*sizeof(char *))) == NULL ||
		 (size_out_shift_labels = (int *) realloc(size_out_shift_labels, (N_shift_out + 1)*sizeof(int))) == NULL ||
		 (out_shift_values_str = (char **) realloc(out_shift_values_str, (N_shift_out + 1)*sizeof(char *))) == NULL ||
		 (size_out_shift_values_str = (int *) realloc(size_out_shift_values_str, (N_shift_out + 1)*sizeof(int))) == NULL ||
		 (Nobs_out_shift = (int * ) realloc(Nobs_out_shift, (N_shift_out + 1)*sizeof(int))) == NULL)
		VARTOOLS_error(ERR_MEMALLOC);
	      for(j = size_out_shift; j < (N_shift_out + 1); j++) {
		size_out_shift_labels[j] = 0;
		out_shift_labels[j] = NULL;
		size_out_shift_values_str[j] = 0;
		out_shift_values_str[j] = NULL;
	      }
	      size_out_shift = N_shift_out + 1;
	    }
	  }
	  if(strlen(stitch->in_shift_labels[vv][inshift_starid][i])+1 > size_out_shift_labels[N_shift_out]) {
	    if(!size_out_shift_labels[N_shift_out]) {
	      size_out_shift_labels[N_shift_out] = strlen(stitch->in_shift_labels[vv][inshift_starid][i])+1;
	      if((out_shift_labels[N_shift_out] = (char *) malloc(size_out_shift_labels[N_shift_out]*sizeof(char))) == NULL)
		VARTOOLS_error(ERR_MEMALLOC);
	    } else {
	      size_out_shift_labels[N_shift_out] = strlen(stitch->in_shift_labels[vv][inshift_starid][i])+1;
	      if((out_shift_labels[N_shift_out] = (char *) realloc(out_shift_labels[N_shift_out], size_out_shift_labels[N_shift_out]*sizeof(char))) == NULL)
		VARTOOLS_error(ERR_MEMALLOC);
	    }
	  }
	  sprintf(out_shift_labels[N_shift_out], "%s", stitch->in_shift_labels[vv][inshift_starid][i]);
	  if(strlen(stitch->in_shift_values_str[vv][inshift_starid][i])+1 > size_out_shift_values_str[N_shift_out]) {
	    if(!size_out_shift_values_str[N_shift_out]) {
	      size_out_shift_values_str[N_shift_out] = strlen(stitch->in_shift_values_str[vv][inshift_starid][i])+1;
	      if((out_shift_values_str[N_shift_out] = (char *) malloc(size_out_shift_values_str[N_shift_out]*sizeof(char))) == NULL)
		VARTOOLS_error(ERR_MEMALLOC);
	    } else {
	      size_out_shift_values_str[N_shift_out] = strlen(stitch->in_shift_values_str[vv][inshift_starid][i])+1;
	      if((out_shift_values_str[N_shift_out] = (char *) realloc(out_shift_values_str[N_shift_out], size_out_shift_values_str[N_shift_out]*sizeof(char))) == NULL)
		VARTOOLS_error(ERR_MEMALLOC);
	    }
	  }
	  sprintf(out_shift_values_str[N_shift_out], "%s", stitch->in_shift_values_str[vv][inshift_starid][i]);
	  Nobs_out_shift[N_shift_out] = stitch->Nobs_in_shifts[vv][inshift_starid][i];
	  N_shift_out += 1;
	}
      }
    }
    /* Add in any newly computed shifts */
    if(!N_shift_out) i = 0;
    else i = 1;
    for(; i < NLCgroups; i++) {
      if(N_shift_out + 1 > size_out_shift) {
	if(!size_out_shift) {
	  size_out_shift = 1;
	  if((out_shift_labels = (char **) malloc(sizeof(char *))) == NULL ||
	     (size_out_shift_labels = (int *) malloc(sizeof(int))) == NULL ||
	     (out_shift_values_str = (char **) malloc(sizeof(char *))) == NULL ||
	     (size_out_shift_values_str = (int *) malloc(sizeof(int))) == NULL ||
	     (Nobs_out_shift = (int * ) malloc(sizeof(int))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	  for(j = 0; j < 1; j++) {
	    size_out_shift_labels[j] = 0;
	    out_shift_labels[j] = NULL;
	    size_out_shift_values_str[j] = 0;
	    out_shift_values_str[j] = NULL;
	  }
	} else {
	  if((out_shift_labels = (char **) realloc(out_shift_labels, (N_shift_out + 1)*sizeof(char *))) == NULL ||
	     (size_out_shift_labels = (int *) realloc(size_out_shift_labels, (N_shift_out + 1)*sizeof(int))) == NULL ||
	     (out_shift_values_str = (char **) realloc(out_shift_values_str, (N_shift_out + 1)*sizeof(char *))) == NULL ||
	     (size_out_shift_values_str = (int *) realloc(size_out_shift_values_str, (N_shift_out + 1)*sizeof(int))) == NULL ||
	     (Nobs_out_shift = (int * ) realloc(Nobs_out_shift, (N_shift_out + 1)*sizeof(int))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	  for(j = size_out_shift; j < (N_shift_out + 1); j++) {
	    size_out_shift_labels[j] = 0;
	    out_shift_labels[j] = NULL;
	    size_out_shift_values_str[j] = 0;
	    out_shift_values_str[j] = NULL;
	  }
	  size_out_shift = N_shift_out + 1;
	}
      }
      if(strlen(stitch->field_labels_vals[lc_num][lcg[i].lcids[0]])+1 > size_out_shift_labels[N_shift_out])
	if(!size_out_shift_labels[N_shift_out]) {
	  size_out_shift_labels[N_shift_out] = strlen(stitch->field_labels_vals[lc_num][lcg[i].lcids[0]])+1;
	  if((out_shift_labels[N_shift_out] = (char *) malloc(size_out_shift_labels[N_shift_out]*sizeof(char))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	} else {
	  size_out_shift_labels[N_shift_out] = strlen(stitch->field_labels_vals[lc_num][lcg[i].lcids[0]])+1;
	  if((out_shift_labels[N_shift_out] = (char *) realloc(out_shift_labels[N_shift_out], size_out_shift_labels[N_shift_out]*sizeof(char))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	}
      sprintf(out_shift_labels[N_shift_out],"%s", stitch->field_labels_vals[lc_num][lcg[i].lcids[0]]);
      if(i > 0) {
	ll = snprintf(NULL, 0, "%.17g", lcg[i].shiftvalue) + 1;
      } else {
	ll = snprintf(NULL, 0, "%.17g", 0.0) + 1;
      }
      if(ll > size_out_shift_values_str[N_shift_out]) {
	if(!size_out_shift_values_str[N_shift_out]) {
	  size_out_shift_values_str[N_shift_out] = ll;
	  if((out_shift_values_str[N_shift_out] = (char *) malloc(size_out_shift_values_str[N_shift_out]*sizeof(char))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	} else {
	  size_out_shift_values_str[N_shift_out] = ll;
	  if((out_shift_values_str[N_shift_out] = (char *) realloc(out_shift_values_str[N_shift_out], size_out_shift_values_str[N_shift_out]*sizeof(char))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	}
      }
      if(i > 0) {
	sprintf(out_shift_values_str[N_shift_out],"%.17g",lcg[i].shiftvalue);
      } else {
	sprintf(out_shift_values_str[N_shift_out],"%.17g",0.0);
      }
      Nobs_out_shift[N_shift_out] = lcg[i].Nlcs;
      N_shift_out += 1;
    }
  }

  *ret_N_shift_out = N_shift_out;
  *ret_size_out_shift = size_out_shift;
  *ret_out_shift_labels = out_shift_labels;
  *ret_size_out_shift_labels = size_out_shift_labels;
  *ret_out_shift_values_str = out_shift_values_str;
  *ret_size_out_shift_values_str = size_out_shift_values_str;
  *ret_Nobs_out_shift = Nobs_out_shift;
  return;
}

void AddtoLCGroup(_StitchLightCurveGroup *lcg, int lcindx)
{
  lcg->Nlcs = lcg->Nlcs + 1;
  if((lcg->lcids = (int *) realloc(lcg->lcids, (lcg->Nlcs * sizeof(int)))) == NULL ||
     (lcg->timegroups = (int *) realloc(lcg->timegroups, (lcg->Nlcs * sizeof(int)))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  lcg->lcids[lcg->Nlcs-1] = lcindx;
}

void AddNewLCGroup(int *ret_NLCgroups, _StitchLightCurveGroup **ret_lcg, 
		   int lcnumval, int lcindx)
{
  _StitchLightCurveGroup *lcg;
  int NLCgroups;

  NLCgroups = *ret_NLCgroups;
  lcg = *ret_lcg;
  if(lcg == NULL) {
    if((lcg = (_StitchLightCurveGroup *) malloc(sizeof(_StitchLightCurveGroup))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    NLCgroups = 1;
  } else {
    if((lcg = (_StitchLightCurveGroup *) realloc(lcg, (NLCgroups + 1)*sizeof(_StitchLightCurveGroup))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    NLCgroups += 1;
  }
  lcg[NLCgroups-1].Nlcs = 1;
  lcg[NLCgroups-1].lcnumval = lcnumval;
  if((lcg[NLCgroups-1].lcids = (int *) malloc(sizeof(int))) == NULL ||
     (lcg[NLCgroups-1].timegroups = (int *) malloc(sizeof(int))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  lcg[NLCgroups-1].lcids[0] = lcindx;

  lcg[NLCgroups-1].Ntimegroups_uniqlist = 0;
  lcg[NLCgroups-1].timegroups_uniqlist = NULL;
  lcg[NLCgroups-1].shiftvalue = 0.;
  
  *ret_NLCgroups = NLCgroups;
  *ret_lcg = lcg;
}

void AddtoTimeGroup(_StitchTimeGroup *tg, int lcgindx, int lcglcindx, int tgindx, _StitchLightCurveGroup *lcg)
{
  int k;
  tg->Nlcs = tg->Nlcs + 1;
  if((tg->lcids = (int *) realloc(tg->lcids, (tg->Nlcs * sizeof(int)))) == NULL ||
     (tg->lcgids = (int *) realloc(tg->lcgids, (tg->Nlcs * sizeof(int)))) == NULL ||
     (tg->lcglcindx = (int *) realloc(tg->lcglcindx, (tg->Nlcs * sizeof(int)))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  tg->lcids[tg->Nlcs-1] = lcg->lcids[lcglcindx];
  tg->lcgids[tg->Nlcs-1] = lcgindx;
  tg->lcglcindx[tg->Nlcs-1] = lcglcindx;
  if(lcgindx > tg->lcgid_uniqlist[tg->Nlcgroups_uniqlist-1]) {
    tg->Nlcgroups_uniqlist += 1;
    if((tg->lcgid_uniqlist = (int *) realloc(tg->lcgid_uniqlist, (tg->Nlcgroups_uniqlist * sizeof(int)))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    tg->lcgid_uniqlist[tg->Nlcgroups_uniqlist - 1] = lcgindx;
  }
  lcg->timegroups[lcglcindx] = tgindx;
  if(lcg->Ntimegroups_uniqlist == 0) {
    lcg->Ntimegroups_uniqlist = 1;
    if((lcg->timegroups_uniqlist = (int *) malloc(sizeof(int))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    lcg->timegroups_uniqlist[lcg->Ntimegroups_uniqlist-1] = tgindx;
  }
  else {
    for(k = 0; k < lcg->Ntimegroups_uniqlist; k++) {
      if(tgindx == lcg->timegroups_uniqlist[k])
	break;
    }
    if(k >= lcg->Ntimegroups_uniqlist) {
      lcg->Ntimegroups_uniqlist += 1;
      if((lcg->timegroups_uniqlist = (int *) realloc(lcg->timegroups_uniqlist, (lcg->Ntimegroups_uniqlist)*sizeof(int))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
      lcg->timegroups_uniqlist[lcg->Ntimegroups_uniqlist - 1] = tgindx;
    }
  }
  return;
}
  
void AddNewTimeGroup(int *ret_Ntimegroups, _StitchTimeGroup **ret_tg,
		     int lcgindx, int lcglcindx, int tbin, 
		     _StitchLightCurveGroup *lcg, double mint, double maxt)
{
  _StitchTimeGroup *tg;
  int Ntimegroups, k;

  Ntimegroups = *ret_Ntimegroups;
  tg = *ret_tg;
  if(tg == NULL) {
    if((tg = (_StitchTimeGroup *) malloc(sizeof(_StitchTimeGroup))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    Ntimegroups = 1;
  } else {
    if((tg = (_StitchTimeGroup *) realloc(tg, (Ntimegroups + 1)*sizeof(_StitchTimeGroup))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    Ntimegroups += 1;
  }
  tg[Ntimegroups-1].Nlcs = 1;
  tg[Ntimegroups-1].tbin = tbin;
  tg[Ntimegroups-1].mint = mint;
  tg[Ntimegroups-1].maxt = maxt;
  if((tg[Ntimegroups-1].lcids = (int *) malloc(sizeof(int))) == NULL ||
     (tg[Ntimegroups-1].lcgids = (int *) malloc(sizeof(int))) == NULL ||
     (tg[Ntimegroups-1].lcglcindx = (int *) malloc(sizeof(int))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  tg[Ntimegroups-1].lcids[0] = lcg->lcids[lcglcindx];
  tg[Ntimegroups-1].lcgids[0] = lcgindx;
  tg[Ntimegroups-1].lcglcindx[0] = lcglcindx;
  tg[Ntimegroups-1].Nlcgroups_uniqlist = 1;
  if((tg[Ntimegroups-1].lcgid_uniqlist = (int *) malloc(sizeof(int))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  tg[Ntimegroups-1].lcgid_uniqlist[0] = lcgindx;
  lcg->timegroups[lcglcindx] = Ntimegroups - 1;
  if(lcg->Ntimegroups_uniqlist == 0) {
    if((lcg->timegroups_uniqlist = (int *) malloc(sizeof(int))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    lcg->timegroups_uniqlist[0] = Ntimegroups - 1;
    lcg->Ntimegroups_uniqlist = 1;
  } else {
    if(Ntimegroups - 1 > lcg->timegroups_uniqlist[lcg->Ntimegroups_uniqlist - 1]) {
      lcg->Ntimegroups_uniqlist += 1;
      if((lcg->timegroups_uniqlist = (int *) realloc(lcg->timegroups_uniqlist, lcg->Ntimegroups_uniqlist * sizeof(int))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
      lcg->timegroups_uniqlist[lcg->Ntimegroups_uniqlist - 1] = Ntimegroups - 1;
    }
  }

  *ret_Ntimegroups = Ntimegroups;
  *ret_tg = tg;
}

void StitchFillLinks(int NLCgroups, _StitchLightCurveGroup *lcg,
	       int Ntimegroups, _StitchTimeGroup *tg,
	       int *islinked, int i, int *Nlinked) {
  int j, k, ll, hh;
  if(!islinked[i]) {
    islinked[i] = 1;
    *Nlinked = (*Nlinked) + 1;
    if(*Nlinked == NLCgroups)
      return;
    for(j = 0; j < lcg[i].Ntimegroups_uniqlist; j++) {
      k = lcg[i].timegroups_uniqlist[j];
      for(ll = 0; ll < tg[k].Nlcgroups_uniqlist; ll++) {
	hh = tg[k].lcgid_uniqlist[ll];
	if(hh != i) {
	  if(!islinked[hh]) {
	    StitchFillLinks(NLCgroups, lcg, Ntimegroups, tg, islinked, hh, Nlinked);
	    if(*Nlinked == NLCgroups)
	      return;
	  }
	}
      }
    }
  }
}

int CheckAllLightCurveGroupsLinked(ProgramData *p, _Stitch *stitch, 
				   int lc_name_num,
				   int lc_num, int NLCgroups, 
				   _StitchLightCurveGroup *lcg,
				   int Ntimegroups, _StitchTimeGroup *tg)
{
  int retval;
  int i, j;
  int *islinked = NULL;
  int Nlinked = 0;
  if(NLCgroups < 1 || Ntimegroups < 1) return 1;
  if((islinked = (int *) malloc(NLCgroups * sizeof(int))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  for(i = 0; i < NLCgroups; i++) islinked[i] = 0;
  StitchFillLinks(NLCgroups, lcg, Ntimegroups, tg, islinked, 0, &Nlinked);
  retval = 1;
  for(i = 0; i < NLCgroups; i++) {
    if(!islinked[i]) {
      retval = 0;
      break;
    }
  }
  free(islinked);
  return(retval);
}

void MergeTwoTimeGroups(ProgramData *p, _Stitch *stitch, int lc_name_num,
			int lc_num, int NLCgroups, _StitchLightCurveGroup *lcg,
			int Ntimegroups, _StitchTimeGroup *tg, int tgindx1,
			int tgindx2)
{
  int i, j, k, kk;
  tg[tgindx1].mint = (tg[tgindx1].mint < tg[tgindx2].mint ? tg[tgindx1].mint : tg[tgindx2].mint);
  tg[tgindx1].maxt = (tg[tgindx1].maxt > tg[tgindx2].maxt ? tg[tgindx1].maxt : tg[tgindx2].maxt);
  if(tg[tgindx1].Nlcs > 0) {
    if(tg[tgindx2].Nlcs > 0) {
      if((tg[tgindx1].lcids = (int *) realloc(tg[tgindx1].lcids, (tg[tgindx1].Nlcs + tg[tgindx2].Nlcs)*sizeof(int))) == NULL ||
	 (tg[tgindx1].lcgids = (int *) realloc(tg[tgindx1].lcgids, (tg[tgindx1].Nlcs + tg[tgindx2].Nlcs)*sizeof(int))) == NULL ||
	 (tg[tgindx1].lcglcindx = (int *) realloc(tg[tgindx1].lcglcindx, (tg[tgindx1].Nlcs + tg[tgindx2].Nlcs)*sizeof(int))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
      
      for(i = tg[tgindx1].Nlcs, j = 0; j < tg[tgindx2].Nlcs; i++, j++) {
	tg[tgindx1].Nlcs++;
	tg[tgindx1].lcids[i] = tg[tgindx2].lcids[j];
	tg[tgindx1].lcgids[i] = tg[tgindx2].lcgids[j];
	tg[tgindx1].lcglcindx[i] = tg[tgindx2].lcglcindx[j];
	lcg[tg[tgindx1].lcgids[i]].timegroups[tg[tgindx1].lcglcindx[i]] = tgindx1;
	for(k = 0; k < tg[tgindx1].Nlcgroups_uniqlist; k++) {
	  if(tg[tgindx1].lcgids[i] == tg[tgindx1].lcgid_uniqlist[k])
	    break;
	}
	if(k >= tg[tgindx1].Nlcgroups_uniqlist) {
	  tg[tgindx1].Nlcgroups_uniqlist += 1;
	  if((tg[tgindx1].lcgid_uniqlist = (int *) realloc(tg[tgindx1].lcgid_uniqlist, tg[tgindx1].Nlcgroups_uniqlist*sizeof(int))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	  tg[tgindx1].lcgid_uniqlist[tg[tgindx1].Nlcgroups_uniqlist - 1] = tg[tgindx1].lcgids[i];
	}
	for(k = 0; k < lcg[tg[tgindx1].lcgids[i]].Ntimegroups_uniqlist; k++) {
	  if(lcg[tg[tgindx1].lcgids[i]].timegroups_uniqlist[k] == tgindx1)
	    break;
	}
	if(k >= lcg[tg[tgindx1].lcgids[i]].Ntimegroups_uniqlist) {
	  lcg[tg[tgindx1].lcgids[i]].Ntimegroups_uniqlist += 1;
	  if((lcg[tg[tgindx1].lcgids[i]].timegroups_uniqlist = (int *) realloc(lcg[tg[tgindx1].lcgids[i]].timegroups_uniqlist, lcg[tg[tgindx1].lcgids[i]].Ntimegroups_uniqlist*sizeof(int))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	  lcg[tg[tgindx1].lcgids[i]].timegroups_uniqlist[lcg[tg[tgindx1].lcgids[i]].Ntimegroups_uniqlist-1] = tgindx1;
	}
	for(k = 0; k < lcg[tg[tgindx1].lcgids[i]].Ntimegroups_uniqlist; k++) {
	  if(lcg[tg[tgindx1].lcgids[i]].timegroups_uniqlist[k] == tgindx2)
	    break;
	}
	if(k < lcg[tg[tgindx1].lcgids[i]].Ntimegroups_uniqlist) {
	  for(kk = k+1; kk < lcg[tg[tgindx1].lcgids[i]].Ntimegroups_uniqlist; kk++) {
	    lcg[tg[tgindx1].lcgids[i]].timegroups_uniqlist[k] = lcg[tg[tgindx1].lcgids[i]].timegroups_uniqlist[kk];
	    k += 1;
	  }
	  lcg[tg[tgindx1].lcgids[i]].Ntimegroups_uniqlist -= 1;
	  if((lcg[tg[tgindx1].lcgids[i]].timegroups_uniqlist = (int *) realloc(lcg[tg[tgindx1].lcgids[i]].timegroups_uniqlist, lcg[tg[tgindx1].lcgids[i]].Ntimegroups_uniqlist*sizeof(int))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	}
      }
      free(tg[tgindx2].lcids);
      free(tg[tgindx2].lcgids);
      free(tg[tgindx2].lcglcindx);
      free(tg[tgindx2].lcgid_uniqlist);
      tg[tgindx2].Nlcs = 0;
      tg[tgindx2].Nlcgroups_uniqlist = 0;
    }
  } else {
    if(tg[tgindx2].Nlcs > 0) {
      tg[tgindx1].Nlcs = tg[tgindx2].Nlcs;
      tg[tgindx1].lcids = tg[tgindx2].lcids;
      tg[tgindx1].lcgids = tg[tgindx2].lcgids;
      tg[tgindx1].lcglcindx = tg[tgindx2].lcglcindx;
      tg[tgindx1].Nlcgroups_uniqlist = tg[tgindx2].Nlcgroups_uniqlist;
      tg[tgindx1].lcgid_uniqlist = tg[tgindx2].lcgid_uniqlist;
      tg[tgindx2].Nlcs = 0;
      tg[tgindx2].lcids = NULL;
      tg[tgindx2].lcgids = NULL;
      tg[tgindx2].lcglcindx = NULL;
      tg[tgindx2].Nlcgroups_uniqlist = 0;
      tg[tgindx2].lcgid_uniqlist = NULL;
      for(i=0; i < tg[tgindx1].Nlcs; i++) {
	lcg[tg[tgindx1].lcgids[i]].timegroups[tg[tgindx1].lcglcindx[i]] = tgindx1;
      }
      for(i=0; i < tg[tgindx1].Nlcgroups_uniqlist; i++) {
	for(k=0; k < lcg[tg[tgindx1].lcgid_uniqlist[i]].Ntimegroups_uniqlist; k++) {
	  if(lcg[tg[tgindx1].lcgid_uniqlist[i]].timegroups_uniqlist[k] == tgindx2) {
	    lcg[tg[tgindx1].lcgid_uniqlist[i]].timegroups_uniqlist[k] = tgindx1;
	    break;
	  }
	}
      }
    }
  }
}     
									       
void MergeAllTimeGroups(ProgramData *p, _Stitch *stitch, int lc_name_num,
		     int lc_num, int NLCgroups, _StitchLightCurveGroup *lcg,
		     int Ntimegroups, _StitchTimeGroup *tg, int *Ntbinindx, 
		     int **tbinindx, int *tbin_factor)
{
  int i, j;
  int Ntbinindxv;
  int *tbinindxv;
  int tbin_factorv;
  int lcgindx;

  Ntbinindxv = *Ntbinindx;
  tbinindxv = *tbinindx;
  tbin_factorv = *tbin_factor;

  for(i = 0, j=0; i < Ntbinindxv; i += 2, j += 1) {
    if(i+1 < Ntbinindxv) {
      if(tbinindxv[i] > -1 && tbinindxv[i+1] > -1) {
	MergeTwoTimeGroups(p, stitch, lc_name_num, lc_num, NLCgroups, 
			   lcg, Ntimegroups, tg, tbinindxv[i], tbinindxv[i+1]);
      }
      if(tbinindxv[i] < 0 && tbinindxv[i+1] > -1) {
	tbinindxv[i] = tbinindxv[i+1];
	tbinindxv[i+1] = -1;
      }
    }
    tbinindxv[j] = tbinindxv[i];
  }
  Ntbinindxv = j;
  *Ntbinindx = Ntbinindxv;
  *tbin_factor = (*tbin_factor)*2;      
}

void FormTimeGroups(ProgramData *p, _Stitch *stitch, int lc_name_num,
		    int lc_num, int NLCgroups, _StitchLightCurveGroup *lcg,
		    int *ret_Ntimegroups, _StitchTimeGroup **ret_tg,
		    int *ret_Nusedtimegroups, int vv)
{
  _StitchTimeGroup *tg = NULL;
  int Ntimegroups = 0;
  int NJD, i, j, k, tmpint;
  double start_time;
  double time_step;
  double t;
  int tbin;
  int maxtbinindx;
  int mintbinindx;
  int maxtbinindx_neg;
  int mintbinindx_neg;
  int *tbinindx = NULL;
  int *tbinindx_neg = NULL;
  int Ntbinindx = -1;
  int tbin_factor = 1;
  int Nusedtimegroups = 0;

  maxtbinindx = -1;
  maxtbinindx_neg = -1;

  NJD = p->NJD[lc_num];

  if(NJD <= 0) return;

  if(stitch->is_start_time) {
    start_time = stitch->start_time;
  } else {
    i = 0;
    while(i < NJD && stitch->stitchmaskvals[vv][lc_num][i] <= VARTOOLS_MASK_TINY)
      i++;
    if(i >= NJD) return;
    start_time = p->t[lc_num][i];
    for(i = i+1; i < NJD; i++) {
      if(stitch->stitchmaskvals[vv][lc_num][i] <= VARTOOLS_MASK_TINY)
	continue;
      if(p->t[lc_num][i] < start_time)
	start_time = p->t[lc_num][i];
    }
  }
  time_step = stitch->time_step;
  if(time_step <= 0.)
    return;

  for(k = 0; k < NLCgroups; k++) {
    for(i = 0; i < lcg[k].Nlcs; i++) {
      t = p->t[lc_num][lcg[k].lcids[i]];
      tbin = (int) ((t - start_time)/time_step);
      if(tbin >= 0) {
	if(tbin > maxtbinindx) {
	  if(maxtbinindx < 0) {
	    maxtbinindx = tbin;
	    if((tbinindx = (int *) malloc((maxtbinindx+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j = 0; j <= maxtbinindx; j++)
	      tbinindx[j] = -1;
	  } else {
	    tmpint = maxtbinindx;
	    maxtbinindx = tbin;
	    if((tbinindx = (int *) realloc(tbinindx, (maxtbinindx+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=tmpint+1; j <= maxtbinindx; j++)
	      tbinindx[j] = -1;
	  }
	}
	if(tbinindx[tbin] < 0) {
	  AddNewTimeGroup(&Ntimegroups, &tg, k, i, tbin, &(lcg[k]), (double) (start_time + tbin*time_step), (double) (start_time + (tbin+1)*time_step));
	  tbinindx[tbin] = Ntimegroups - 1;
	} else {
	  AddtoTimeGroup(&(tg[tbinindx[tbin]]), k, i, tbinindx[tbin], &(lcg[k]));
	}
      } else {
	if(-tbin > maxtbinindx_neg) {
	  if(maxtbinindx_neg < 0) {
	    maxtbinindx_neg = -tbin;
	    if((tbinindx_neg = (int *) malloc((maxtbinindx_neg+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j = 0; j <= maxtbinindx_neg; j++)
	      tbinindx_neg[j] = -1;
	  } else {
	    tmpint = maxtbinindx_neg;
	    maxtbinindx_neg = -tbin;
	    if((tbinindx_neg = (int *) realloc(tbinindx_neg, (maxtbinindx_neg+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=tmpint+1; j <= maxtbinindx_neg; j++)
	      tbinindx_neg[j] = -1;
	  }
	}
	if(tbinindx_neg[-tbin] < 0) {
	  AddNewTimeGroup(&Ntimegroups, &tg, k, i, tbin, &(lcg[k]), (double) (start_time + tbin*time_step), (double) (start_time + (tbin+1)*time_step));
	  tbinindx_neg[-tbin] = Ntimegroups - 1;
	} else {
	  AddtoTimeGroup(&(tg[tbinindx_neg[-tbin]]), k, i, tbinindx_neg[-tbin], &(lcg[k]));
	}
      }
    }
  }
  
  if(maxtbinindx_neg > -1) {
    if(maxtbinindx > -1) {
      Ntbinindx = maxtbinindx + maxtbinindx_neg + 1;
      if((tbinindx = (int *) realloc(tbinindx, Ntbinindx*sizeof(int))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
      for(k = Ntbinindx-1; k >= maxtbinindx_neg; k--) {
	tbinindx[k] = tbinindx[k-maxtbinindx_neg];
      }
      for(k = maxtbinindx_neg-1; k >= 0; k--) {
	tbinindx[k] = tbinindx_neg[maxtbinindx_neg-k];
      }
    }
    else {
      mintbinindx_neg = maxtbinindx_neg;
      for(i = 1; i <= maxtbinindx_neg; i++) {
	if(tbinindx_neg[i] > -1) {
	  mintbinindx_neg = i;
	  break;
	}
      }
      Ntbinindx = maxtbinindx_neg-mintbinindx_neg+1;
      if((tbinindx = (int *) malloc(Ntbinindx*sizeof(int))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
      for(k = 0; k < Ntbinindx; k++) {
	tbinindx[k] = tbinindx_neg[maxtbinindx_neg-k];
      }
    }
  } else {
    if(maxtbinindx > -1) {
      mintbinindx = maxtbinindx;
      for(i = 0; i <= maxtbinindx; i++) {
	if(tbinindx[i] > -1) {
	  mintbinindx = i;
	  break;
	}
      }
      Ntbinindx = maxtbinindx - mintbinindx+1;
      if(i > 0) {
	for(k = 0; k < Ntbinindx; k++) {
	  tbinindx[k] = tbinindx[k+mintbinindx];
	}
      }
    }
  }

  while(Ntbinindx > 1 && 
	!CheckAllLightCurveGroupsLinked(p, stitch, lc_name_num,
					lc_num, NLCgroups, 
					lcg,
					Ntimegroups, tg)) {
    MergeAllTimeGroups(p, stitch, lc_name_num,
		       lc_num, NLCgroups, lcg,
		       Ntimegroups, tg, &Ntbinindx, 
		       &tbinindx, &tbin_factor);
  }

  Nusedtimegroups = 0;
  for(i = 0; i < Ntimegroups; i++) {
    if(tg[i].Nlcs > 0)
      Nusedtimegroups++;
  }

  stitch->final_time_step[lc_num] = ((double)tbin_factor)*stitch->time_step;

  if(tbinindx != NULL) free(tbinindx);
  if(tbinindx_neg != NULL) free(tbinindx_neg);
  *ret_Ntimegroups = Ntimegroups;
  *ret_Nusedtimegroups = Nusedtimegroups;
  *ret_tg = tg;
}

void FormOnlyOneTimeGroup(ProgramData *p, _Stitch *stitch, int lc_name_num,
			  int lc_num, int NLCgroups, _StitchLightCurveGroup *lcg,
			  int *ret_Ntimegroups, _StitchTimeGroup **ret_tg,
			  int *ret_Nusedtimegroups, int vv)
{
  _StitchTimeGroup *tg = NULL;
  int Ntimegroups = 0;
  int NJD, i, j, jj, k, NJDtouse;
  double start_time;
  double end_time;
  double time_step;
  double t;
  int tbin;
  int Nusedtimegroups = 0;

  NJD = p->NJD[lc_num];

  if(NJD <= 0) return;

  i = 0;
  while(i < NJD && stitch->stitchmaskvals[vv][lc_num][i] <= VARTOOLS_MASK_TINY)
    i++;
  if(i >= NJD) return;


  Ntimegroups = 1;
  NJDtouse = 1;
  start_time = p->t[lc_num][i];
  end_time = start_time;
  for(i = i+1; i < NJD; i++) {
    if(stitch->stitchmaskvals[vv][lc_num][i] <= VARTOOLS_MASK_TINY)
      continue;
    NJDtouse++;
    if(p->t[lc_num][i] < start_time)
      start_time = p->t[lc_num][i];
    if(p->t[lc_num][i] > end_time)
      end_time = p->t[lc_num][i];
  }
  
  if((tg = (_StitchTimeGroup *) malloc(sizeof(_StitchTimeGroup))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  
  time_step = end_time - start_time;

  tg->Nlcs = NJDtouse;
  if((tg->lcids = (int *) malloc(NJDtouse*sizeof(int))) == NULL ||
     (tg->lcgids = (int *) malloc(NJDtouse*sizeof(int))) == NULL ||
     (tg->lcglcindx = (int *) malloc(NJDtouse*sizeof(int))) == NULL ||
     (tg->lcgid_uniqlist = (int *) malloc(NLCgroups * sizeof(int))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  
  tg->mint = start_time;
  tg->maxt = end_time;
  tg->tbin = 0;

  j = 0;
  jj = 0;
  for(k = 0; k < NLCgroups; k++) {
    if(lcg[k].Nlcs > 0) {
      tg->lcgid_uniqlist[jj] = k;
      jj++;
    }
    for(i = 0; i < lcg[k].Nlcs; i++) {
      tg->lcids[j] = lcg[k].lcids[i];
      tg->lcgids[j] = k;
      tg->lcglcindx[j] = i;
      j++;
    }
  }

  tg->Nlcgroups_uniqlist = jj;

  Nusedtimegroups = 0;
  for(i = 0; i < Ntimegroups; i++) {
    if(tg[i].Nlcs > 0)
      Nusedtimegroups++;
  }

  *ret_Ntimegroups = Ntimegroups;
  *ret_Nusedtimegroups = Nusedtimegroups;
  *ret_tg = tg;
}
	

void FormLCGroups(ProgramData *p, _Stitch *stitch, int lc_name_num, 
			 int lc_num, int *ret_NLCgroups, _StitchLightCurveGroup **ret_lcg, int vv)
{
  int *groupidx;
  int NJD, i, j;
  int *lcindx;
  int Nlcindx;
  int tmpint;
  int maxlcindx;
  int maxlcindx_neg;
  int maxrefnumindx;
  int minrefnumindx;
  int *lcgrpindx = NULL;
  int *lcgrpindx_neg = NULL;
  int *Nlcindxvals;
  int lcnumvaltouse;

  _StitchLightCurveGroup *lcg = NULL;
  int NLCgroups = 0;

  NJD = p->NJD[lc_num];

  maxlcindx = -1;
  maxlcindx_neg = -1;

  if(NJD <= 0) return;

  if(!stitch->userefnum) {
    for(i = 0; i < NJD; i++) {
      if(stitch->stitchmaskvals[vv][lc_num][i] <= VARTOOLS_MASK_TINY) 
	continue;
      if(stitch->lcnumval[lc_num][i] < 0) {
	if(-(stitch->lcnumval[lc_num][i]) > maxlcindx_neg) {
	  if(maxlcindx_neg < 0) {
	    maxlcindx_neg = -(stitch->lcnumval[lc_num][i]);
	    if((lcgrpindx_neg = (int *) malloc((maxlcindx_neg+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=0; j <= maxlcindx_neg; j++)
	      lcgrpindx[j] = -1;
	  } else {
	    tmpint = maxlcindx_neg;
	    maxlcindx_neg = -(stitch->lcnumval[lc_num][i]);
	    if((lcgrpindx_neg = (int *) realloc(lcgrpindx_neg, (maxlcindx_neg+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=tmpint+1; j <= maxlcindx_neg; j++)
	      lcgrpindx_neg[j] = -1;
	  }
	}
	if(lcgrpindx_neg[-(stitch->lcnumval[lc_num][i])] < 0) {
	  AddNewLCGroup(&NLCgroups, &lcg, stitch->lcnumval[lc_num][i], i);
	  lcgrpindx_neg[-(stitch->lcnumval[lc_num][i])] = NLCgroups - 1;
	} else {
	  AddtoLCGroup(&(lcg[lcgrpindx_neg[-(stitch->lcnumval[lc_num][i])]]), i);
	}
      } else {
	if(stitch->lcnumval[lc_num][i] > maxlcindx) {
	  if(maxlcindx < 0) {
	    maxlcindx = stitch->lcnumval[lc_num][i];
	    if((lcgrpindx = (int *) malloc((maxlcindx+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=0; j <= maxlcindx; j++)
	      lcgrpindx[j] = -1;
	  } else {
	    tmpint = maxlcindx;
	    maxlcindx = stitch->lcnumval[lc_num][i];
	    if((lcgrpindx = (int *) realloc(lcgrpindx, (maxlcindx+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=tmpint+1; j <= maxlcindx; j++)
	      lcgrpindx[j] = -1;
	  }
	}
	if(lcgrpindx[stitch->lcnumval[lc_num][i]] < 0) {
	  AddNewLCGroup(&NLCgroups, &lcg, stitch->lcnumval[lc_num][i], i);
	  lcgrpindx[stitch->lcnumval[lc_num][i]] = NLCgroups - 1;
	} else {
	  AddtoLCGroup(&(lcg[lcgrpindx[stitch->lcnumval[lc_num][i]]]), i);
	}
      }
    }
  } else {
    for(i=0; i < NJD; i++) {
      if(stitch->stitchmaskvals[vv][lc_num][i] <= VARTOOLS_MASK_TINY) 
	continue;
      maxrefnumindx = stitch->refnumval[lc_num][i];
      minrefnumindx = stitch->refnumval[lc_num][i];
      break;
    }
    if(i == NJD) return;
    for(; i < NJD; i++) {
      if(stitch->stitchmaskvals[vv][lc_num][i] <= VARTOOLS_MASK_TINY) 
	continue;
      if(stitch->refnumval[lc_num][i] > maxrefnumindx)
	maxrefnumindx = stitch->refnumval[lc_num][i];
      if(stitch->refnumval[lc_num][i] < minrefnumindx)
	minrefnumindx = stitch->refnumval[lc_num][i];
    }
    stitch->minrefnumindx[lc_num] = minrefnumindx;
    stitch->maxrefnumindx[lc_num] = maxrefnumindx;
    for(i = 0; i < NJD; i++) {
      if(stitch->stitchmaskvals[vv][lc_num][i] <= VARTOOLS_MASK_TINY) 
	continue;
      lcnumvaltouse = stitch->lcnumval[lc_num][i]*(maxrefnumindx - minrefnumindx + 1) + (stitch->refnumval[lc_num][i] - minrefnumindx);
      if(lcnumvaltouse < 0) {
	if(-(lcnumvaltouse) > maxlcindx_neg) {
	  if(maxlcindx_neg < 0) {
	    maxlcindx_neg = -(lcnumvaltouse);
	    if((lcgrpindx_neg = (int *) malloc((maxlcindx_neg+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=0; j <= maxlcindx_neg; j++)
	      lcgrpindx[j] = -1;
	  } else {
	    tmpint = maxlcindx_neg;
	    maxlcindx_neg = -(lcnumvaltouse);
	    if((lcgrpindx_neg = (int *) realloc(lcgrpindx_neg, (maxlcindx_neg+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=tmpint+1; j <= maxlcindx_neg; j++)
	      lcgrpindx_neg[j] = -1;
	  }
	}
	if(lcgrpindx_neg[-(lcnumvaltouse)] < 0) {
	  AddNewLCGroup(&NLCgroups, &lcg, lcnumvaltouse, i);
	  lcgrpindx_neg[-(lcnumvaltouse)] = NLCgroups - 1;
	} else {
	  AddtoLCGroup(&(lcg[lcgrpindx_neg[-(lcnumvaltouse)]]), i);
	}
      } else {
	if(lcnumvaltouse > maxlcindx) {
	  if(maxlcindx < 0) {
	    maxlcindx = lcnumvaltouse;
	    if((lcgrpindx = (int *) malloc((maxlcindx+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=0; j <= maxlcindx; j++)
	      lcgrpindx[j] = -1;
	  } else {
	    tmpint = maxlcindx;
	    maxlcindx = lcnumvaltouse;
	    if((lcgrpindx = (int *) realloc(lcgrpindx, (maxlcindx+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=tmpint+1; j <= maxlcindx; j++)
	      lcgrpindx[j] = -1;
	  }
	}
	if(lcgrpindx[lcnumvaltouse] < 0) {
	  AddNewLCGroup(&NLCgroups, &lcg, lcnumvaltouse, i);
	  lcgrpindx[lcnumvaltouse] = NLCgroups - 1;
	} else {
	  AddtoLCGroup(&(lcg[lcgrpindx[lcnumvaltouse]]), i);
	}
      }
    }
  }
  *ret_NLCgroups = NLCgroups;
  *ret_lcg = lcg;
  if(lcgrpindx != NULL) free(lcgrpindx);
  if(lcgrpindx_neg != NULL) free(lcgrpindx_neg);
}

void FormLCGroups_IsInShift(ProgramData *p, _Stitch *stitch, int lc_name_num, 
			    int lc_num, int *ret_NLCgroups, 
			    _StitchLightCurveGroup **ret_lcg, int vv, 
			    int *fieldlabels_indx)
{
  int *groupidx;
  int NJD, i, j;
  int *lcindx;
  int Nlcindx;
  int tmpint;
  int maxlcindx;
  int maxlcindx_neg;
  int maxrefnumindx;
  int minrefnumindx;
  int *lcgrpindx = NULL;
  int *lcgrpindx_neg = NULL;
  int *Nlcindxvals;
  int lcnumvaltouse;

  _StitchLightCurveGroup *lcg = NULL;
  int NLCgroups = 0;
  int lcnumid_to_use_fordefault = 0;

  NJD = p->NJD[lc_num];

  maxlcindx = -1;
  maxlcindx_neg = -1;

  if(NJD <= 0) return;

  /* First set all of the points from a previous field with an accepted shift
     into the default group */
  for(i=0; i < NJD; i++) {
    if(stitch->stitchmaskvals[vv][lc_num][i] <= VARTOOLS_MASK_TINY) 
      continue;
    if(fieldlabels_indx[i] >= 0) {
      if(!NLCgroups) {
	AddNewLCGroup(&NLCgroups, &lcg, stitch->lcnumval[lc_num][i], i);
	lcnumid_to_use_fordefault = stitch->lcnumval[lc_num][i];
      } else {
	AddtoLCGroup(&(lcg[0]), i);
      }
    }
  }
  
  *ret_NLCgroups = NLCgroups;
  *ret_lcg = lcg;

  /* Next add any remaining points with un-determined shifts */

  if(!stitch->userefnum) {
    for(i = 0; i < NJD; i++) {
      if(stitch->stitchmaskvals[vv][lc_num][i] <= VARTOOLS_MASK_TINY) 
	continue;
      if(fieldlabels_indx[i] >= 0)
	continue;
      if(stitch->lcnumval[lc_num][i] < 0) {
	if(-(stitch->lcnumval[lc_num][i]) > maxlcindx_neg) {
	  if(maxlcindx_neg < 0) {
	    maxlcindx_neg = -(stitch->lcnumval[lc_num][i]);
	    if((lcgrpindx_neg = (int *) malloc((maxlcindx_neg+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=0; j <= maxlcindx_neg; j++)
	      lcgrpindx[j] = -1;
	  } else {
	    tmpint = maxlcindx_neg;
	    maxlcindx_neg = -(stitch->lcnumval[lc_num][i]);
	    if((lcgrpindx_neg = (int *) realloc(lcgrpindx_neg, (maxlcindx_neg+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=tmpint+1; j <= maxlcindx_neg; j++)
	      lcgrpindx_neg[j] = -1;
	  }
	}
	if(lcgrpindx_neg[-(stitch->lcnumval[lc_num][i])] < 0) {
	  AddNewLCGroup(&NLCgroups, &lcg, stitch->lcnumval[lc_num][i], i);
	  lcgrpindx_neg[-(stitch->lcnumval[lc_num][i])] = NLCgroups - 1;
	} else {
	  AddtoLCGroup(&(lcg[lcgrpindx_neg[-(stitch->lcnumval[lc_num][i])]]), i);
	}
      } else {
	if(stitch->lcnumval[lc_num][i] > maxlcindx) {
	  if(maxlcindx < 0) {
	    maxlcindx = stitch->lcnumval[lc_num][i];
	    if((lcgrpindx = (int *) malloc((maxlcindx+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=0; j <= maxlcindx; j++)
	      lcgrpindx[j] = -1;
	  } else {
	    tmpint = maxlcindx;
	    maxlcindx = stitch->lcnumval[lc_num][i];
	    if((lcgrpindx = (int *) realloc(lcgrpindx, (maxlcindx+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=tmpint+1; j <= maxlcindx; j++)
	      lcgrpindx[j] = -1;
	  }
	}
	if(lcgrpindx[stitch->lcnumval[lc_num][i]] < 0) {
	  AddNewLCGroup(&NLCgroups, &lcg, stitch->lcnumval[lc_num][i], i);
	  lcgrpindx[stitch->lcnumval[lc_num][i]] = NLCgroups - 1;
	} else {
	  AddtoLCGroup(&(lcg[lcgrpindx[stitch->lcnumval[lc_num][i]]]), i);
	}
      }
    }
  } else {
    for(i=0; i < NJD; i++) {
      if(stitch->stitchmaskvals[vv][lc_num][i] <= VARTOOLS_MASK_TINY) 
	continue;
      if(fieldlabels_indx[i] >= 0)
	continue;
      maxrefnumindx = stitch->refnumval[lc_num][i];
      minrefnumindx = stitch->refnumval[lc_num][i];
      break;
    }
    if(i == NJD) return;
    for(; i < NJD; i++) {
      if(stitch->stitchmaskvals[vv][lc_num][i] <= VARTOOLS_MASK_TINY) 
	continue;
      if(fieldlabels_indx[i] >= 0)
	continue;
      if(stitch->refnumval[lc_num][i] > maxrefnumindx)
	maxrefnumindx = stitch->refnumval[lc_num][i];
      if(stitch->refnumval[lc_num][i] < minrefnumindx)
	minrefnumindx = stitch->refnumval[lc_num][i];
    }
    stitch->minrefnumindx[lc_num] = minrefnumindx;
    stitch->maxrefnumindx[lc_num] = maxrefnumindx;
    for(i = 0; i < NJD; i++) {
      if(stitch->stitchmaskvals[vv][lc_num][i] <= VARTOOLS_MASK_TINY) 
	continue;
      if(fieldlabels_indx[i] >= 0)
	continue;
      lcnumvaltouse = stitch->lcnumval[lc_num][i]*(maxrefnumindx - minrefnumindx + 1) + (stitch->refnumval[lc_num][i] - minrefnumindx);
      if(lcnumvaltouse < 0) {
	if(-(lcnumvaltouse) > maxlcindx_neg) {
	  if(maxlcindx_neg < 0) {
	    maxlcindx_neg = -(lcnumvaltouse);
	    if((lcgrpindx_neg = (int *) malloc((maxlcindx_neg+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=0; j <= maxlcindx_neg; j++)
	      lcgrpindx[j] = -1;
	  } else {
	    tmpint = maxlcindx_neg;
	    maxlcindx_neg = -(lcnumvaltouse);
	    if((lcgrpindx_neg = (int *) realloc(lcgrpindx_neg, (maxlcindx_neg+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=tmpint+1; j <= maxlcindx_neg; j++)
	      lcgrpindx_neg[j] = -1;
	  }
	}
	if(lcgrpindx_neg[-(lcnumvaltouse)] < 0) {
	  AddNewLCGroup(&NLCgroups, &lcg, lcnumvaltouse, i);
	  lcgrpindx_neg[-(lcnumvaltouse)] = NLCgroups - 1;
	} else {
	  AddtoLCGroup(&(lcg[lcgrpindx_neg[-(lcnumvaltouse)]]), i);
	}
      } else {
	if(lcnumvaltouse > maxlcindx) {
	  if(maxlcindx < 0) {
	    maxlcindx = lcnumvaltouse;
	    if((lcgrpindx = (int *) malloc((maxlcindx+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=0; j <= maxlcindx; j++)
	      lcgrpindx[j] = -1;
	  } else {
	    tmpint = maxlcindx;
	    maxlcindx = lcnumvaltouse;
	    if((lcgrpindx = (int *) realloc(lcgrpindx, (maxlcindx+1)*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(j=tmpint+1; j <= maxlcindx; j++)
	      lcgrpindx[j] = -1;
	  }
	}
	if(lcgrpindx[lcnumvaltouse] < 0) {
	  AddNewLCGroup(&NLCgroups, &lcg, lcnumvaltouse, i);
	  lcgrpindx[lcnumvaltouse] = NLCgroups - 1;
	} else {
	  AddtoLCGroup(&(lcg[lcgrpindx[lcnumvaltouse]]), i);
	}
      }
    }
  }
  *ret_NLCgroups = NLCgroups;
  *ret_lcg = lcg;
  if(lcgrpindx != NULL) free(lcgrpindx);
  if(lcgrpindx_neg != NULL) free(lcgrpindx_neg);
}


void FreeLCgroups(int *ret_NLCgroups, _StitchLightCurveGroup **ret_lcg)
{
  int NLCgroups;
  int i;
  _StitchLightCurveGroup *lcg;
  lcg = *ret_lcg;
  NLCgroups = *ret_NLCgroups;
  for(i = 0; i < NLCgroups; i++) {
    if(lcg[i].Nlcs > 0) {
      free(lcg[i].lcids);
      free(lcg[i].timegroups);
    }
    if(lcg[i].Ntimegroups_uniqlist > 0)
      free(lcg[i].timegroups_uniqlist);
  }
  free(lcg);
  *ret_NLCgroups = 0;
  *ret_lcg = NULL;
}

void Freetimegroups(int *ret_Ntimegroups, _StitchTimeGroup **ret_tg)
{
  int Ntimegroups;
  int i;
  _StitchTimeGroup *tg;
  tg = *ret_tg;
  Ntimegroups = *ret_Ntimegroups;
  for(i = 0; i < Ntimegroups; i++) {
    if(tg[i].Nlcs > 0) {
      free(tg[i].lcids);
      free(tg[i].lcgids);
      free(tg[i].lcglcindx);
    }
    if(tg[i].Nlcgroups_uniqlist > 0)
      free(tg[i].lcgid_uniqlist);
  }
  free(tg);
  *ret_Ntimegroups = 0;
  *ret_tg = NULL;
}

double stitch_get_weightedmean_uncertainty(int N, double *errvals) {
  int j;
  double var2 = 0.;

  for(j = 0; j < N; j++)
    var2 += 1./errvals[j]/errvals[j];
  return sqrt(1./var2);
}

double stitch_get_mean_uncertainty(int N, double *errvals) {
  int j;
  double var1 = 0.;

  for(j = 0; j < N; j++)
    var1 += errvals[j]*errvals[j];
  return (sqrt(var1) / ((double) N));
}

void SetStitchResultsNothingToDo(ProgramData *p, _Stitch *stitch, 
				 int lc_name_num, int lc_num, int NLCgroups,
				 _StitchLightCurveGroup *lcg, int Ntimegroups, 
				 int Nusedtimegroups, _StitchTimeGroup *tg, 
				 int vv)
{
  stitch->Nlcgroups_used[lc_num] = NLCgroups;
  stitch->Ntimegroups_used[lc_num] = 0;
  stitch->Nparamtotal_used[lc_num] = 0;
  if(stitch->groupbytime) {
    stitch->final_time_step[lc_num] = 0.0;
  }
}

void StitchAddStitchParamsToHeader(ProgramData *p, _Stitch *stitch,
			       int lc_name_num, int lc_num)
{
  char fitshdrkeyword[MAXLEN];
  char fitshdrcomment[MAXLEN];
  char fitshdrstringval[MAXLEN];
  int fitshdrintval;
  double fitshdrdblval;
  if(!stitch->add_stitchparams_fitsheader) return;

  sprintf(fitshdrkeyword, "STCHMTHD");
  sprintf(fitshdrcomment, "Method used for lc stitching in vartools");
  switch(stitch->stitchmethod) {
  case VARTOOLS_STITCH_METHOD_MEDIAN: 
    sprintf(fitshdrstringval,"median");
    break;
  case VARTOOLS_STITCH_METHOD_MEAN:
    sprintf(fitshdrstringval,"mean");
    break;
  case VARTOOLS_STITCH_METHOD_WEIGHTEDMEAN:
    sprintf(fitshdrstringval,"weightedmean");
    break;
  case VARTOOLS_STITCH_METHOD_POLY:
    sprintf(fitshdrstringval,"poly");
    break;
  case VARTOOLS_STITCH_METHOD_HARM:
    sprintf(fitshdrstringval,"harmseries");
    break;
  default:
    VARTOOLS_error(ERR_CODEERROR);
    break;
  }
  VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_stitchparams_fitsheader_hdutouse, stitch->add_stitchparams_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_STRING, fitshdrstringval);

  if(stitch->stitchmethod == VARTOOLS_STITCH_METHOD_POLY) {
    sprintf(fitshdrkeyword, "STCHPLYO");
    sprintf(fitshdrcomment, "Poly order for lc stitching in vartools");
    fitshdrintval = stitch->polyorder;
    VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_stitchparams_fitsheader_hdutouse, stitch->add_stitchparams_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_INT, fitshdrintval);
  } else if(stitch->stitchmethod == VARTOOLS_STITCH_METHOD_HARM) {
    sprintf(fitshdrkeyword, "STCHHRMP");
    sprintf(fitshdrcomment, "Harmonic period for stitching in vartools");
    fitshdrdblval = stitch->harmperiodvals[lc_num];
    VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_stitchparams_fitsheader_hdutouse, stitch->add_stitchparams_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_DOUBLE, fitshdrdblval);
    
    sprintf(fitshdrkeyword, "STCHHRMN");
    sprintf(fitshdrcomment, "Nharm for harmonic stitching in vartools");
    fitshdrintval = stitch->Nharm;
    VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_stitchparams_fitsheader_hdutouse, stitch->add_stitchparams_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_INT, fitshdrintval);
  }

  sprintf(fitshdrkeyword, "STCHTGRP");
  sprintf(fitshdrcomment, "groupbytime used for stitching? (1=yes; 0=no)");
  fitshdrintval = stitch->groupbytime;
  VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_stitchparams_fitsheader_hdutouse, stitch->add_stitchparams_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_INT, fitshdrintval);

  if(stitch->groupbytime) {
    sprintf(fitshdrkeyword, "STCHTBIN");
    sprintf(fitshdrcomment, "Time bin used by vartools stitching");
    fitshdrdblval = stitch->time_step;
    VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_stitchparams_fitsheader_hdutouse, stitch->add_stitchparams_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_DOUBLE, fitshdrdblval);
    if(stitch->groupbytime && stitch->is_start_time) {
      sprintf(fitshdrkeyword, "STCHTSRT");
      sprintf(fitshdrcomment, "Time bin start used by vartools stitching");
      fitshdrdblval = stitch->start_time;
      VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_stitchparams_fitsheader_hdutouse, stitch->add_stitchparams_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_DOUBLE, fitshdrdblval);
    }
  }
}

void StitchAddInShiftsToHeader(ProgramData *p, _Stitch *stitch,
			       int lc_name_num, int lc_num, int NLCgroups,
			       _StitchLightCurveGroup *lcg, int Ntimegroups,
			       int Nusedtimegroups, _StitchTimeGroup *tg,
			       int vv, int *fieldlabels_indx, int inshift_starid)
{
  char fitshdrkeyword[MAXLEN];
  char fitshdrcomment[MAXLEN];
  char varnumchar;
  char lcgroupnumchar;
  char lcrefnumchar;
  char lcrefnumchar1;
  char lcrefnumchar2;
  int i, ll, kk, i1, i2;
  int *fieldlabeldoneyet = NULL;
  char *nextrefidchar1touse = NULL;
  char *nextrefidchar2touse = NULL;
  int minlcnum;
  int maxlcnum;
  int Nlcnum;

  if(inshift_starid < 0 || !stitch->add_shifts_fitsheader) return;

  if((fieldlabeldoneyet = (int *) malloc(stitch->N_shifts_per_star[vv][inshift_starid]*sizeof(int))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  for(i = 0; i < stitch->N_shifts_per_star[vv][inshift_starid]; i++) fieldlabeldoneyet[i] = 0;

  if(p->NJD[lc_num] > 0) {
    minlcnum = stitch->lcnumval[lc_num][0];
    maxlcnum = stitch->lcnumval[lc_num][0];
    for(i=1; i < p->NJD[lc_num]; i++) {
      if(stitch->lcnumval[lc_num][i] < minlcnum) minlcnum = stitch->lcnumval[lc_num][i];
      if(stitch->lcnumval[lc_num][i] > maxlcnum) maxlcnum = stitch->lcnumval[lc_num][i];
    }
    Nlcnum = (maxlcnum - minlcnum) + 1;
    if((nextrefidchar1touse = (char *) malloc(Nlcnum*sizeof(char))) == NULL ||
       (nextrefidchar2touse = (char *) malloc(Nlcnum*sizeof(char))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    for(i=0; i < Nlcnum; i++) {
      nextrefidchar1touse[i] = 'A';
      nextrefidchar2touse[i] = 'A';
    }
  }
  /* Add parameter values to the output light curve FITS header if requested */
  if(vv < 26)
    varnumchar = 'A' + (char) vv;
  else
    varnumchar = 'A' + (char) (vv % 26);
  sprintf(fitshdrkeyword, "INSHFT%c", varnumchar);

  sprintf(fitshdrcomment, "Input shift file for variable %d", vv+1);
  if(!stitch->is_inshifts_header_basename_only) {
    VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_shifts_fitsheader_hdutouse, stitch->add_shifts_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_STRING, stitch->in_shifts_filename[vv]);
  } else {
    i1 = 0;
    i2 = 0;
    while(stitch->in_shifts_filename[vv][i1] != '\n' && stitch->in_shifts_filename[vv][i1] != '\0') {
      if(stitch->in_shifts_filename[vv][i1] == '/')
	i2 = i1+1;
      i1++;
    }
    VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_shifts_fitsheader_hdutouse, stitch->add_shifts_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_STRING, &(stitch->in_shifts_filename[vv][i2]));
  }

  if(!stitch->userefnum) {
    for(i=0; i < p->NJD[lc_num]; i++) {
      if(fieldlabels_indx[i] < 0)
	continue;
      if(!fieldlabeldoneyet[fieldlabels_indx[i]]) {
	if(vv < 26)
	  varnumchar = 'A' + (char) vv;
	else
	  varnumchar = 'A' + (char) (vv % 26);
	if(stitch->lcnumval[lc_num][i] < 26)
	  lcgroupnumchar = 'A' + (char) (stitch->lcnumval[lc_num][i]);
	else
	  lcgroupnumchar = 'A' + (char) ((stitch->lcnumval[lc_num][i]) % 26);
	sprintf(fitshdrkeyword, "%s%c%c", stitch->keywordbase, varnumchar, lcgroupnumchar);
	sprintf(fitshdrcomment, "Shift for variable %d LCgroup %d", vv, stitch->lcnumval[lc_num][i]);
	VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_shifts_fitsheader_hdutouse, stitch->add_shifts_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_DOUBLE, stitch->in_shift_values[vv][inshift_starid][fieldlabels_indx[i]]);
	fieldlabeldoneyet[fieldlabels_indx[i]] = 1;
      }
    }
  } else {
    for(i=0; i < p->NJD[lc_num]; i++) {
      if(fieldlabels_indx[i] < 0)
	continue;
      if(!fieldlabeldoneyet[fieldlabels_indx[i]]) {
	if(vv < 26)
	  varnumchar = 'A' + (char) vv;
	else
	  varnumchar = 'A' + (char) (vv % 26);
	ll = stitch->lcnumval[lc_num][i];
	kk = stitch->refnumval[lc_num][i];
	if(ll < 26)
	  lcgroupnumchar = 'A' + (char) ll;
	else
	  lcgroupnumchar = 'A' + (char) (ll % 26);
	lcrefnumchar1 = nextrefidchar1touse[stitch->lcnumval[lc_num][i]-minlcnum];
	lcrefnumchar2 = nextrefidchar2touse[stitch->lcnumval[lc_num][i]-minlcnum];
	nextrefidchar2touse[stitch->lcnumval[lc_num][i]-minlcnum] += (char) 1;
	if(nextrefidchar2touse[stitch->lcnumval[lc_num][i]-minlcnum] > 'Z') {
	  nextrefidchar2touse[stitch->lcnumval[lc_num][i]-minlcnum] = 'A';
	  nextrefidchar1touse[stitch->lcnumval[lc_num][i]-minlcnum] += (char) 1;
	  if(nextrefidchar1touse[stitch->lcnumval[lc_num][i]-minlcnum] > 'Z') {
	    nextrefidchar1touse[stitch->lcnumval[lc_num][i]-minlcnum] = 'A';
	  }
	}
	sprintf(fitshdrkeyword, "%s%c%c%c%c", stitch->keywordbase, varnumchar, lcgroupnumchar, lcrefnumchar1, lcrefnumchar2);
	sprintf(fitshdrcomment, "Shift for variable %d LCgroup %d REFID %d", vv, ll, kk);
	VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_shifts_fitsheader_hdutouse, stitch->add_shifts_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_DOUBLE, stitch->in_shift_values[vv][inshift_starid][fieldlabels_indx[i]]);
	fieldlabeldoneyet[fieldlabels_indx[i]] = 1;
      }
    }
  }

  if(fieldlabeldoneyet != NULL)
    free(fieldlabeldoneyet);
  if(nextrefidchar1touse != NULL)
    free(nextrefidchar1touse);
  if(nextrefidchar2touse != NULL)
    free(nextrefidchar2touse);
  return;
}

void  StitchByStatistic(ProgramData *p, _Stitch *stitch, int lc_name_num, int lc_num, int NLCgroups, _StitchLightCurveGroup *lcg, int Ntimegroups, int Nusedtimegroups, _StitchTimeGroup *tg, int vv, FILE *coeffoutfile)
{
  double **stats_vecs = NULL, **err_vecs = NULL;
  double *stats_vals = NULL, *err_vals = NULL;
  int Nstats_vecs, j, i, k, ll, kk;
  int *length_stats_vecs = NULL;
  int *idx = NULL;
  int ndecorr;
  int npoints, ntgused;
  double *model_vec = NULL;
  double **decorr_terms = NULL;
  int *order_terms = NULL;
  double *Avector = NULL;
  double *A_errvector = NULL;
  double mint;
  double maxt;
  char fitshdrkeyword[MAXLEN];
  char fitshdrcomment[MAXLEN];
  char varnumchar;
  char lcgroupnumchar;
  char lcrefnumchar1;
  char lcrefnumchar2;
  int minlcnum;
  int maxlcnum;
  int Nlcnum;
  char *nextrefidchar1touse = NULL;
  char *nextrefidchar2touse = NULL;
  

  if(NLCgroups < 2) {
    /* Nothing to do */
    SetStitchResultsNothingToDo(p, stitch, lc_name_num, lc_num, NLCgroups,
				lcg, Ntimegroups, Nusedtimegroups, tg, vv);
    return;
  }
  Nstats_vecs = 0;
  for(i = 0; i < Ntimegroups; i++) {
    if(tg[i].Nlcgroups_uniqlist > 1) {
      Nstats_vecs += tg[i].Nlcgroups_uniqlist;
    }
  }
  if(Nstats_vecs < 1) {
    SetStitchResultsNothingToDo(p, stitch, lc_name_num, lc_num, NLCgroups,
				lcg, Ntimegroups, Nusedtimegroups, tg, vv);
    return;
  }


  /* Form the statistics vectors */
  if((stats_vecs = (double **) malloc(Nstats_vecs * sizeof(double *))) == NULL ||
     (err_vecs = (double **) malloc(Nstats_vecs * sizeof(double *))) == NULL ||
     (length_stats_vecs = (int *) malloc(Nstats_vecs * sizeof(int))) == NULL ||
     (idx = (int *) malloc(NLCgroups * sizeof(int))) == NULL ||
     (stats_vals = (double *) malloc(Nstats_vecs * sizeof(double))) == NULL ||
     (err_vals = (double *) malloc(Nstats_vecs * sizeof(double))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  ll = 0;
  ndecorr = 0;
  for(i = 0; i < Nstats_vecs; i++) {
    length_stats_vecs[i] = 0;
    stats_vecs[i] = NULL;
    err_vecs[i] = NULL;
  }
  for(i = 0; i < Ntimegroups; i++) {
    if(tg[i].Nlcgroups_uniqlist > 1) {
      ndecorr++;
      for(k = 0; k < NLCgroups; k++)
	idx[k] = -1;
      for(k = 0; k < tg[i].Nlcgroups_uniqlist; k++)
	idx[tg[i].lcgid_uniqlist[k]] = k;
      for(k = 0; k < tg[i].Nlcs; k++) {
	length_stats_vecs[ll + idx[tg[i].lcgids[k]]]++;
      }
      for(k = 0; k < tg[i].Nlcgroups_uniqlist; k++) {
	if(length_stats_vecs[ll + k] > 0) {
	  if((stats_vecs[ll + k] = (double *) malloc(length_stats_vecs[ll + k]*sizeof(double))) == NULL ||
	     (err_vecs[ll + k] = (double *) malloc(length_stats_vecs[ll + k]*sizeof(double))) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	}
	else {
	  stats_vecs[ll + k] = NULL;
	  err_vecs[ll + k] = NULL;
	}
	length_stats_vecs[ll+k] = 0;
      }
      for(k = 0; k < tg[i].Nlcs; k++) {
	stats_vecs[ll + idx[tg[i].lcgids[k]]][length_stats_vecs[ll + idx[tg[i].lcgids[k]]]] = stitch->stitchvarvals[vv][lc_num][tg[i].lcids[k]];
	err_vecs[ll + idx[tg[i].lcgids[k]]][length_stats_vecs[ll + idx[tg[i].lcgids[k]]]] = stitch->stitcherrvals[vv][lc_num][tg[i].lcids[k]];
	length_stats_vecs[ll + idx[tg[i].lcgids[k]]]++;
      }
      ll += tg[i].Nlcgroups_uniqlist;
    }
  }
  
  ntgused = ndecorr;
  ndecorr += (NLCgroups - 1);
  npoints = 0;
  for(i = 0; i < Nstats_vecs; i++) {
    if(length_stats_vecs[i] > 0) {
      switch(stitch->stitchmethod) {
      case VARTOOLS_STITCH_METHOD_MEDIAN:
	stats_vals[npoints] = VARTOOLS_median(length_stats_vecs[i], stats_vecs[i]);
	err_vals[npoints] = sqrt(M_PI / 2.0) * stitch_get_mean_uncertainty(length_stats_vecs[i], err_vecs[i]);
	break;
      case VARTOOLS_STITCH_METHOD_MEAN:
	stats_vals[npoints] = VARTOOLS_getmean(length_stats_vecs[i], stats_vecs[i]);
	err_vals[npoints] = stitch_get_mean_uncertainty(length_stats_vecs[i], err_vecs[i]);
	break;
      case VARTOOLS_STITCH_METHOD_WEIGHTEDMEAN:
	stats_vals[npoints] = VARTOOLS_getweightedmean(length_stats_vecs[i], stats_vecs[i], err_vecs[i]);
	err_vals[npoints] = stitch_get_weightedmean_uncertainty(length_stats_vecs[i], err_vecs[i]);
	break;
      default:
	VARTOOLS_error(ERR_CODEERROR);
	break;
      }
      npoints++;
    }
  }

  if(ndecorr < 1 || npoints < 1) {
    SetStitchResultsNothingToDo(p, stitch, lc_name_num, lc_num, NLCgroups,
				lcg, Ntimegroups, Nusedtimegroups, tg, vv);
    for(i=0; i < Nstats_vecs; i++) {
      if(stats_vecs[i] != NULL) free(stats_vecs[i]);
      if(err_vecs[i] != NULL) free(err_vecs[i]);
    }
    if(stats_vecs != NULL) free(stats_vecs);
    if(err_vecs != NULL) free(err_vecs);
    if(length_stats_vecs != NULL) free(length_stats_vecs);
    if(idx != NULL) free(idx);
    if(stats_vals != NULL) free(stats_vals);
    if(err_vals != NULL) free(err_vals);
    return;
  }
  
  if((decorr_terms = (double **) malloc(npoints*sizeof(double))) == NULL ||
     (order_terms = (int *) malloc(ndecorr*sizeof(int))) == NULL ||
     (Avector =  (double *) malloc(ndecorr*sizeof(double))) == NULL ||
     (A_errvector = (double *) malloc(ndecorr*sizeof(double))) == NULL ||
     (model_vec = (double *) malloc(npoints*sizeof(double))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  for(i = 0; i < npoints; i++) {
    if((decorr_terms[i] = (double *) malloc(ndecorr*sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    for(j=0; j < ndecorr; j++) {
      decorr_terms[i][j] = 0.;
    }
  }
  for(j=0; j < ndecorr; j++) order_terms[j] = 1;

  ll = 0;
  j = 0;
  for(i = 0; i < Ntimegroups; i++) {
    if(tg[i].Nlcgroups_uniqlist > 1) {
      for(k = 0; k < tg[i].Nlcgroups_uniqlist; k++) {
	decorr_terms[j][ll] = 1.;
	if(tg[i].lcgid_uniqlist[k] > 0) {
	  decorr_terms[j][ntgused + tg[i].lcgid_uniqlist[k] - 1] = 1.;
	}
	j++;
      }
      ll++;
    }
  }

  VARTOOLS_docorr(stats_vals, err_vals, npoints, ndecorr, decorr_terms, order_terms, Avector, A_errvector, 0., 0., 0, NULL, lc_name_num, lc_num);
  
  /* Fill in the results and correct the LC as needed */

  for(i = 1; i < NLCgroups; i++) {
    lcg[i].shiftvalue = Avector[ntgused + i - 1];
  }
  
  if(!stitch->fitonly) {
    for(i = 1; i < NLCgroups; i++) {
      for(k = 0; k < lcg[i].Nlcs; k++) {
	stitch->stitchvarvals[vv][lc_num][lcg[i].lcids[k]] -= Avector[ntgused + i - 1];
      }
    }
  }

  /* Save the parameters to an output file if requested */

  if(stitch->save_fitted_parameters && coeffoutfile != NULL) {
    fprintf(coeffoutfile,"# Parameters for stitch variable %d\n", vv+1);
    ll = 0;
    for(i = 0; i < Ntimegroups; i++) {
      if(tg[i].Nlcgroups_uniqlist > 1) {
	switch(stitch->stitchmethod) {
	case VARTOOLS_STITCH_METHOD_MEDIAN:
	  fprintf(coeffoutfile,"Median ");
	  break;
	case VARTOOLS_STITCH_METHOD_MEAN:
	  fprintf(coeffoutfile,"Mean ");
	  break;
	case VARTOOLS_STITCH_METHOD_WEIGHTEDMEAN:
	  fprintf(coeffoutfile,"WeightedMean ");
	  break;
	default:
	  VARTOOLS_error(ERR_CODEERROR);
	  break;
	}
	fprintf(coeffoutfile,"%.17g<t<%.17g: ", tg[i].mint, tg[i].maxt);
	fprintf(coeffoutfile,"%.17g\n", Avector[ll]);
	ll++;
      }
    }
    for(i=1; i < NLCgroups; i++) {
      fprintf(coeffoutfile,"LCgroup_%d shift: %.17g\n", i+1, -Avector[ntgused + i - 1]);
    }
    fprintf(coeffoutfile,"\n");
  }

  /* Add parameter values to the output light curve FITS header if requested */
  if(stitch->add_shifts_fitsheader) {
    if(!stitch->userefnum) {
      for(i=1; i < NLCgroups; i++) {
	if(vv < 26)
	  varnumchar = 'A' + (char) vv;
	else
	  varnumchar = 'A' + (char) (vv % 26);
	if(lcg[i].lcnumval < 26)
	  lcgroupnumchar = 'A' + (char) (lcg[i].lcnumval);
	else
	  lcgroupnumchar = 'A' + (char) ((lcg[i].lcnumval) % 26);
	sprintf(fitshdrkeyword, "%s%c%c", stitch->keywordbase, varnumchar, lcgroupnumchar);
	sprintf(fitshdrcomment, "Shift for variable %d LCgroup %d", vv, lcg[i].lcnumval);
	VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_shifts_fitsheader_hdutouse, stitch->add_shifts_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_DOUBLE, -Avector[ntgused + i - 1]);
      }
    }
    else {
      if(p->NJD[lc_num] > 0) {
	minlcnum = stitch->lcnumval[lc_num][0];
	maxlcnum = stitch->lcnumval[lc_num][0];
	for(i=1; i < p->NJD[lc_num]; i++) {
	  if(stitch->lcnumval[lc_num][i] < minlcnum) minlcnum = stitch->lcnumval[lc_num][i];
	  if(stitch->lcnumval[lc_num][i] > maxlcnum) maxlcnum = stitch->lcnumval[lc_num][i];
	}
	Nlcnum = (maxlcnum - minlcnum) + 1;
	if((nextrefidchar1touse = (char *) malloc(Nlcnum*sizeof(char))) == NULL ||
	   (nextrefidchar2touse = (char *) malloc(Nlcnum*sizeof(char))) == NULL)
	  VARTOOLS_error(ERR_MEMALLOC);
	for(i=0; i < Nlcnum; i++) {
	  nextrefidchar1touse[i] = 'A';
	  nextrefidchar2touse[i] = 'A';
	}
      }

      for(i=1; i < NLCgroups; i++) {
	if(vv < 26)
	  varnumchar = 'A' + (char) vv;
	else
	  varnumchar = 'A' + (char) (vv % 26);
	ll = (lcg[i].lcnumval / (stitch->maxrefnumindx[lc_num] - stitch->minrefnumindx[lc_num] + 1));
	kk = ((lcg[i].lcnumval % (stitch->maxrefnumindx[lc_num] - stitch->minrefnumindx[lc_num] + 1)));
	if(ll < 26)
	  lcgroupnumchar = 'A' + (char) ll;
	else
	  lcgroupnumchar = 'A' + (char) (ll % 26);
	lcrefnumchar1 = nextrefidchar1touse[ll];
	lcrefnumchar2 = nextrefidchar2touse[ll];
	nextrefidchar2touse[ll] += (char) 1;
	if(nextrefidchar2touse[ll] > 'Z') {
	  nextrefidchar2touse[ll] = 'A';
	  nextrefidchar1touse[ll] += (char) 1;
	  if(nextrefidchar1touse[ll] > 'Z') {
	    nextrefidchar1touse[ll] = 'A';
	  }
	}
	sprintf(fitshdrkeyword, "%s%c%c%c%c", stitch->keywordbase, varnumchar, lcgroupnumchar, lcrefnumchar1, lcrefnumchar2);
	sprintf(fitshdrcomment, "Shift for variable %d LCgroup %d REFID %d", vv, ll, kk + stitch->minrefnumindx[lc_num]);
	VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_shifts_fitsheader_hdutouse, stitch->add_shifts_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_DOUBLE, -Avector[ntgused + i - 1]);
      }
    }
  }

  stitch->Nlcgroups_used[lc_num] = NLCgroups;
  stitch->Ntimegroups_used[lc_num] = ntgused;
  stitch->Nparamtotal_used[lc_num] = ndecorr;

  /* Free allocated memory, and return  */
  for(i=0; i < Nstats_vecs; i++) {
    if(stats_vecs[i] != NULL) free(stats_vecs[i]);
    if(err_vecs[i] != NULL) free(err_vecs[i]);
  }
  if(stats_vecs != NULL) free(stats_vecs);
  if(err_vecs != NULL) free(err_vecs);
  if(length_stats_vecs != NULL) free(length_stats_vecs);
  if(idx != NULL) free(idx);
  if(stats_vals != NULL) free(stats_vals);
  if(err_vals != NULL) free(err_vals);

  for(i = 0; i < npoints; i++) {
    if(decorr_terms[i] != NULL) free(decorr_terms[i]);
  }

  if(decorr_terms != NULL) free(decorr_terms);
  if(order_terms != NULL) free(order_terms);
  if(Avector != NULL) free(Avector);
  if(A_errvector != NULL) free(A_errvector);
  if(model_vec != NULL) free(model_vec);
  if(nextrefidchar1touse != NULL) free(nextrefidchar1touse);
  if(nextrefidchar2touse != NULL) free(nextrefidchar2touse);

  return;

}

void  StitchByPoly(ProgramData *p, _Stitch *stitch, int lc_name_num, int lc_num, int NLCgroups, _StitchLightCurveGroup *lcg, int Ntimegroups, int Nusedtimegroups, _StitchTimeGroup *tg, int vv, FILE *coeffoutfile)
{
  double *data_vals = NULL, *err_vals = NULL;
  int j, i, k, ll, ii, kk;
  int dA;
  int ndecorr;
  int npoints, ntgused, ndecorrfrompoly;
  double **decorr_terms = NULL;
  int *order_terms = NULL;
  double *Avector = NULL;
  double *A_errvector = NULL;
  double mint;
  double maxt;
  double dt;
  double dtpow;
  char fitshdrkeyword[MAXLEN];
  char fitshdrcomment[MAXLEN];
  char varnumchar;
  char lcgroupnumchar;
  char lcrefnumchar;
  char lcrefnumchar1;
  char lcrefnumchar2;
  int minlcnum;
  int maxlcnum;
  int Nlcnum;
  char *nextrefidchar1touse = NULL;
  char *nextrefidchar2touse = NULL;

  if(NLCgroups < 2) {
    /* Nothing to do */
    SetStitchResultsNothingToDo(p, stitch, lc_name_num, lc_num, NLCgroups,
				lcg, Ntimegroups, Nusedtimegroups, tg, vv);
    return;
  }
  /* Determine the total number of points to be included in the fit */
  npoints = 0;
  for(i = 0; i < Ntimegroups; i++) {
    if(tg[i].Nlcgroups_uniqlist > 1) {
      npoints += tg[i].Nlcs;
    }
  }

  if(npoints < 2) {
    /* Nothing to do */
    SetStitchResultsNothingToDo(p, stitch, lc_name_num, lc_num, NLCgroups,
				lcg, Ntimegroups, Nusedtimegroups, tg, vv);
    return;
  }
  
  ll = 0;
  ndecorr = 0;

  for(i = 0; i < Ntimegroups; i++) {
    if(tg[i].Nlcgroups_uniqlist > 1) {
      ndecorr++;
    }
  }
  
  ntgused = ndecorr;

  dA = stitch->polyorder+1;
  ndecorr = ndecorr*dA;

  ndecorrfrompoly = ndecorr;

  ndecorr += (NLCgroups - 1);

  if(ndecorr < 1 || npoints < 1) {
    SetStitchResultsNothingToDo(p, stitch, lc_name_num, lc_num, NLCgroups,
				lcg, Ntimegroups, Nusedtimegroups, tg, vv);
    return;
  }
  
  if((data_vals = (double *) malloc(npoints*sizeof(double))) == NULL ||
     (err_vals = (double *) malloc(npoints*sizeof(double))) == NULL ||
     (decorr_terms = (double **) malloc(npoints*sizeof(double))) == NULL ||
     (order_terms = (int *) malloc(ndecorr*sizeof(int))) == NULL ||
     (Avector =  (double *) malloc(ndecorr*sizeof(double))) == NULL ||
     (A_errvector = (double *) malloc(ndecorr*sizeof(double))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  for(i = 0; i < npoints; i++) {
    if((decorr_terms[i] = (double *) malloc(ndecorr*sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    for(j=0; j < ndecorr; j++) {
      decorr_terms[i][j] = 0.;
    }
  }
  for(j=0; j < ndecorr; j++) order_terms[j] = 1;

  ll = 0;
  j = 0;
  ii = 0;
  for(i = 0; i < Ntimegroups; i++) {
    if(tg[i].Nlcgroups_uniqlist > 1) {
      for(k = 0; k < tg[i].Nlcs; k++) {
	data_vals[j] = stitch->stitchvarvals[vv][lc_num][tg[i].lcids[k]];
	err_vals[j] = stitch->stitcherrvals[vv][lc_num][tg[i].lcids[k]];

	dt = p->t[lc_num][tg[i].lcids[k]] - tg[i].mint;
	ll = ii*dA;
	decorr_terms[j][ll] = 1.;
	dtpow = dt;
	for(kk=1; kk <= stitch->polyorder; kk++) {
	  decorr_terms[j][ll+kk] = dtpow;
	  dtpow = dtpow*dt;
	}
	if(tg[i].lcgids[k] > 0) {
	  ll = ndecorrfrompoly + (tg[i].lcgids[k] - 1);
	  decorr_terms[j][ll] = 1.;
	}
	j++;
      }
      ii++;
    }
  }

  VARTOOLS_docorr(data_vals, err_vals, npoints, ndecorr, decorr_terms, order_terms, Avector, A_errvector, 0., 0., 0, NULL, lc_name_num, lc_num);
  
  /* Fill in the results and correct the LC as needed */

  for(i = 1; i < NLCgroups; i++) {
    lcg[i].shiftvalue = Avector[ndecorrfrompoly + i - 1];
  }
  
  if(!stitch->fitonly) {
    for(i = 1; i < NLCgroups; i++) {
      for(k = 0; k < lcg[i].Nlcs; k++) {
	stitch->stitchvarvals[vv][lc_num][lcg[i].lcids[k]] -= Avector[ndecorrfrompoly + i - 1];
      }
    }
  }

  /* Save the parameters to an output file if requested */

  if(stitch->save_fitted_parameters && coeffoutfile != NULL) {
    fprintf(coeffoutfile,"# Parameters for stitch variable %d\n", vv+1);
    ll = 0;
    for(i = 0; i < Ntimegroups; i++) {
      if(tg[i].Nlcgroups_uniqlist > 1) {
	for(k = 0; k <= stitch->polyorder; k++) {
	  fprintf(coeffoutfile,"Coeff for t^%d, ",k);
	  fprintf(coeffoutfile,"%.17g<t<%.17g: ", tg[i].mint, tg[i].maxt);
	  fprintf(coeffoutfile,"%.17g\n", Avector[ll]);
	  ll++;
	}
      }
    }
    for(i=1; i < NLCgroups; i++) {
      fprintf(coeffoutfile,"LCgroup_%d shift: %.17g\n", i+1, -Avector[ndecorrfrompoly + i - 1]);
    }
    fprintf(coeffoutfile,"\n");
  }

  /* Add parameter values to the output light curve FITS header if requested */
  if(stitch->add_shifts_fitsheader) {
    if(!stitch->userefnum) {
      for(i=1; i < NLCgroups; i++) {
	if(vv < 26)
	  varnumchar = 'A' + (char) vv;
	else
	  varnumchar = 'A' + (char) (vv % 26);
	if(lcg[i].lcnumval < 26)
	  lcgroupnumchar = 'A' + (char) (lcg[i].lcnumval);
	else
	  lcgroupnumchar = 'A' + (char) ((lcg[i].lcnumval) % 26);
	sprintf(fitshdrkeyword, "%s%c%c", stitch->keywordbase, varnumchar, lcgroupnumchar);
	sprintf(fitshdrcomment, "Shift for variable number %d LCgroup %d", vv, lcg[i].lcnumval);
	VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_shifts_fitsheader_hdutouse, stitch->add_shifts_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_DOUBLE, -Avector[ndecorrfrompoly + i - 1]);
      }
    } else {
      if(p->NJD[lc_num] > 0) {
	minlcnum = stitch->lcnumval[lc_num][0];
	maxlcnum = stitch->lcnumval[lc_num][0];
	for(i=1; i < p->NJD[lc_num]; i++) {
	  if(stitch->lcnumval[lc_num][i] < minlcnum) minlcnum = stitch->lcnumval[lc_num][i];
	  if(stitch->lcnumval[lc_num][i] > maxlcnum) maxlcnum = stitch->lcnumval[lc_num][i];
	}
	Nlcnum = (maxlcnum - minlcnum) + 1;
	if((nextrefidchar1touse = (char *) malloc(Nlcnum*sizeof(char))) == NULL ||
	   (nextrefidchar2touse = (char *) malloc(Nlcnum*sizeof(char))) == NULL)
	  VARTOOLS_error(ERR_MEMALLOC);
	for(i=0; i < Nlcnum; i++) {
	  nextrefidchar1touse[i] = 'A';
	  nextrefidchar2touse[i] = 'A';
	}
      }
      for(i=1; i < NLCgroups; i++) {
	if(vv < 26)
	  varnumchar = 'A' + (char) vv;
	else
	  varnumchar = 'A' + (char) (vv % 26);
	ll = (lcg[i].lcnumval / (stitch->maxrefnumindx[lc_num] - stitch->minrefnumindx[lc_num] + 1));
	kk = ((lcg[i].lcnumval % (stitch->maxrefnumindx[lc_num] - stitch->minrefnumindx[lc_num] + 1)));
	if(ll < 26)
	  lcgroupnumchar = 'A' + (char) ll;
	else
	  lcgroupnumchar = 'A' + (char) (ll % 26);
	lcrefnumchar1 = nextrefidchar1touse[ll];
	lcrefnumchar2 = nextrefidchar2touse[ll];
	nextrefidchar2touse[ll] += (char) 1;
	if(nextrefidchar2touse[ll] > 'Z') {
	  nextrefidchar2touse[ll] = 'A';
	  nextrefidchar1touse[ll] += (char) 1;
	  if(nextrefidchar1touse[ll] > 'Z') {
	    nextrefidchar1touse[ll] = 'A';
	  }
	}
	sprintf(fitshdrkeyword, "%s%c%c%c%c", stitch->keywordbase, varnumchar, lcgroupnumchar, lcrefnumchar1, lcrefnumchar2);
	sprintf(fitshdrcomment, "Shift for variable number %d LCgroup %d REFID %d", vv, ll, kk + stitch->minrefnumindx[lc_num]);
	VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_shifts_fitsheader_hdutouse, stitch->add_shifts_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_DOUBLE, -Avector[ndecorrfrompoly + i - 1]);
      }
    }
  }

  stitch->Nlcgroups_used[lc_num] = NLCgroups;
  stitch->Ntimegroups_used[lc_num] = ntgused;
  stitch->Nparamtotal_used[lc_num] = ndecorr;

  /* Free allocated memory, and return  */
  if(data_vals != NULL) free(data_vals);
  if(err_vals != NULL) free(err_vals);

  for(i = 0; i < npoints; i++) {
    if(decorr_terms[i] != NULL) free(decorr_terms[i]);
  }

  if(decorr_terms != NULL) free(decorr_terms);
  if(order_terms != NULL) free(order_terms);
  if(Avector != NULL) free(Avector);
  if(A_errvector != NULL) free(A_errvector);

  if(nextrefidchar1touse != NULL) free(nextrefidchar1touse);
  if(nextrefidchar2touse != NULL) free(nextrefidchar2touse);
  return;
}


void  StitchByHarm(ProgramData *p, _Stitch *stitch, int lc_name_num, int lc_num, int NLCgroups, _StitchLightCurveGroup *lcg, int Ntimegroups, int Nusedtimegroups, _StitchTimeGroup *tg, int vv, FILE *coeffoutfile)
{
  double *data_vals = NULL, *err_vals = NULL;
  int j, i, k, ll, ii, kk;
  int dA;
  int ndecorr;
  int npoints, ntgused, ndecorrfromharm;
  double **decorr_terms = NULL;
  int *order_terms = NULL;
  double *Avector = NULL;
  double *A_errvector = NULL;
  double mint;
  double maxt;
  double dt;
  double freq;
  double freq0;
  double periodtouse;
  char fitshdrkeyword[MAXLEN];
  char fitshdrcomment[MAXLEN];
  char varnumchar;
  char lcgroupnumchar;
  char lcrefnumchar;
  char lcrefnumchar1;
  char lcrefnumchar2;
  int minlcnum;
  int maxlcnum;
  int Nlcnum;
  char *nextrefidchar1touse = NULL;
  char *nextrefidchar2touse = NULL;
  if(NLCgroups < 2) {
    /* Nothing to do */
    SetStitchResultsNothingToDo(p, stitch, lc_name_num, lc_num, NLCgroups,
				lcg, Ntimegroups, Nusedtimegroups, tg, vv);
    return;
  }
  /* Determine the total number of points to be included in the fit */
  npoints = 0;
  for(i = 0; i < Ntimegroups; i++) {
    if(tg[i].Nlcgroups_uniqlist > 1) {
      npoints += tg[i].Nlcs;
    }
  }

  if(npoints < 2) {
    /* Nothing to do */
    SetStitchResultsNothingToDo(p, stitch, lc_name_num, lc_num, NLCgroups,
				lcg, Ntimegroups, Nusedtimegroups, tg, vv);
    return;
  }
  
  ll = 0;

  periodtouse = stitch->harmperiodvals[lc_num];

  ndecorr = 0;

  for(i = 0; i < Ntimegroups; i++) {
    if(tg[i].Nlcgroups_uniqlist > 1) {
      ndecorr++;
    }
  }
  
  ntgused = ndecorr;

  dA = (2*(stitch->Nharm+1))+1;
  ndecorr = ndecorr*dA;

  ndecorrfromharm = ndecorr;

  ndecorr += (NLCgroups - 1);

  if(ndecorr < 1 || npoints < 1) {
    SetStitchResultsNothingToDo(p, stitch, lc_name_num, lc_num, NLCgroups,
				lcg, Ntimegroups, Nusedtimegroups, tg, vv);
    return;
  }
  
  if((data_vals = (double *) malloc(npoints*sizeof(double))) == NULL ||
     (err_vals = (double *) malloc(npoints*sizeof(double))) == NULL ||
     (decorr_terms = (double **) malloc(npoints*sizeof(double))) == NULL ||
     (order_terms = (int *) malloc(ndecorr*sizeof(int))) == NULL ||
     (Avector =  (double *) malloc(ndecorr*sizeof(double))) == NULL ||
     (A_errvector = (double *) malloc(ndecorr*sizeof(double))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  for(i = 0; i < npoints; i++) {
    if((decorr_terms[i] = (double *) malloc(ndecorr*sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    for(j=0; j < ndecorr; j++) {
      decorr_terms[i][j] = 0.;
    }
  }
  for(j=0; j < ndecorr; j++) order_terms[j] = 1;

  ll = 0;
  j = 0;
  ii = 0;
  for(i = 0; i < Ntimegroups; i++) {
    if(tg[i].Nlcgroups_uniqlist > 1) {
      for(k = 0; k < tg[i].Nlcs; k++) {
	data_vals[j] = stitch->stitchvarvals[vv][lc_num][tg[i].lcids[k]];
	err_vals[j] = stitch->stitcherrvals[vv][lc_num][tg[i].lcids[k]];

	dt = p->t[lc_num][tg[i].lcids[k]] - tg[i].mint;
	ll = ii*dA;
	decorr_terms[j][ll] = 1.;
	ll++;
	freq0 = dt*2.0*M_PI/periodtouse;
	for(kk=1; kk <= stitch->Nharm + 1; kk++) {
	  freq = kk*freq0;
	  decorr_terms[j][ll+(2*(kk-1))] = cos(freq);
	  decorr_terms[j][ll+(2*(kk-1))+1] = sin(freq);
	}
	if(tg[i].lcgids[k] > 0) {
	  ll = ndecorrfromharm + (tg[i].lcgids[k] - 1);
	  decorr_terms[j][ll] = 1.;
	}
	j++;
      }
      ii++;
    }
  }

  VARTOOLS_docorr(data_vals, err_vals, npoints, ndecorr, decorr_terms, order_terms, Avector, A_errvector, 0., 0., 0, NULL, lc_name_num, lc_num);
  
  /* Fill in the results and correct the LC as needed */

  for(i = 1; i < NLCgroups; i++) {
    lcg[i].shiftvalue = Avector[ndecorrfromharm + i - 1];
  }
  
  if(!stitch->fitonly) {
    for(i = 1; i < NLCgroups; i++) {
      for(k = 0; k < lcg[i].Nlcs; k++) {
	stitch->stitchvarvals[vv][lc_num][lcg[i].lcids[k]] -= Avector[ndecorrfromharm + i - 1];
      }
    }
  }

  /* Save the parameters to an output file if requested */

  if(stitch->save_fitted_parameters && coeffoutfile != NULL) {
    fprintf(coeffoutfile,"# Parameters for stitch variable %d\n", vv+1);
    ll = 0;
    for(i = 0; i < Ntimegroups; i++) {
      if(tg[i].Nlcgroups_uniqlist > 1) {
	fprintf(coeffoutfile,"Coeff for 1, ");
	fprintf(coeffoutfile,"%.17g<t<%.17g: ", tg[i].mint, tg[i].maxt);
	fprintf(coeffoutfile,"%.17g\n", Avector[ll]);
	ll++;
	for(k = 1; k <= stitch->Nharm+1; k++) {
	  fprintf(coeffoutfile,"Coeff for cos(%d*t*f0), ",k);
	  fprintf(coeffoutfile,"%.17g<t<%.17g: ", tg[i].mint, tg[i].maxt);
	  fprintf(coeffoutfile,"%.17g\n", Avector[ll]);
	  ll++;
	  fprintf(coeffoutfile,"Coeff for sin(%d*t*f0), ",k);
	  fprintf(coeffoutfile,"%.17g<t<%.17g: ", tg[i].mint, tg[i].maxt);
	  fprintf(coeffoutfile,"%.17g\n", Avector[ll]);
	  ll++;
	}
      }
    }
    for(i=1; i < NLCgroups; i++) {
      fprintf(coeffoutfile,"LCgroup_%d shift: %.17g\n", i+1, -Avector[ndecorrfromharm + i - 1]);
    }
    fprintf(coeffoutfile,"\n");
  }

  /* Add parameter values to the output light curve FITS header if requested */
  if(stitch->add_shifts_fitsheader) {
    if(!stitch->userefnum) {
      for(i=1; i < NLCgroups; i++) {
	if(vv < 26)
	  varnumchar = 'A' + (char) vv;
	else
	  varnumchar = 'A' + (char) (vv % 26);
	if(lcg[i].lcnumval < 26)
	  lcgroupnumchar = 'A' + (char) (lcg[i].lcnumval);
	else
	  lcgroupnumchar = 'A' + (char) ((lcg[i].lcnumval) % 26);
	sprintf(fitshdrkeyword, "%s%c%c", stitch->keywordbase, varnumchar, lcgroupnumchar);
	sprintf(fitshdrcomment, "Shift for variable number %d LCgroup %d", vv, (lcg[i].lcnumval));
	VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_shifts_fitsheader_hdutouse, stitch->add_shifts_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_DOUBLE, -Avector[ndecorrfromharm + i - 1]);
      }
    } else {
      if(p->NJD[lc_num] > 0) {
	minlcnum = stitch->lcnumval[lc_num][0];
	maxlcnum = stitch->lcnumval[lc_num][0];
	for(i=1; i < p->NJD[lc_num]; i++) {
	  if(stitch->lcnumval[lc_num][i] < minlcnum) minlcnum = stitch->lcnumval[lc_num][i];
	  if(stitch->lcnumval[lc_num][i] > maxlcnum) maxlcnum = stitch->lcnumval[lc_num][i];
	}
	Nlcnum = (maxlcnum - minlcnum) + 1;
	if((nextrefidchar1touse = (char *) malloc(Nlcnum*sizeof(char))) == NULL ||
	   (nextrefidchar2touse = (char *) malloc(Nlcnum*sizeof(char))) == NULL)
	  VARTOOLS_error(ERR_MEMALLOC);
	for(i=0; i < Nlcnum; i++) {
	  nextrefidchar1touse[i] = 'A';
	  nextrefidchar2touse[i] = 'A';
	}
      }
      for(i=1; i < NLCgroups; i++) {
	if(vv < 26)
	  varnumchar = 'A' + (char) vv;
	else
	  varnumchar = 'A' + (char) (vv % 26);
	ll = (lcg[i].lcnumval / (stitch->maxrefnumindx[lc_num] - stitch->minrefnumindx[lc_num] + 1));
	kk = ((lcg[i].lcnumval % (stitch->maxrefnumindx[lc_num] - stitch->minrefnumindx[lc_num] + 1)));

	if(ll < 26)
	  lcgroupnumchar = 'A' + (char) ll;
	else
	  lcgroupnumchar = 'A' + (char) (ll % 26);
	lcrefnumchar1 = nextrefidchar1touse[ll];
	lcrefnumchar2 = nextrefidchar2touse[ll];
	nextrefidchar2touse[ll] += (char) 1;
	if(nextrefidchar2touse[ll] > 'Z') {
	  nextrefidchar2touse[ll] = 'A';
	  nextrefidchar1touse[ll] += (char) 1;
	  if(nextrefidchar1touse[ll] > 'Z') {
	    nextrefidchar1touse[ll] = 'A';
	  }
	}
	sprintf(fitshdrkeyword, "%s%c%c%c%c", stitch->keywordbase, varnumchar, lcgroupnumchar, lcrefnumchar1, lcrefnumchar2);
	
	/*ll = (lcg[i].lcnumval / (stitch->maxrefnumindx[lc_num] - stitch->minrefnumindx[lc_num] + 1));
	kk = ((lcg[i].lcnumval % (stitch->maxrefnumindx[lc_num] - stitch->minrefnumindx[lc_num] + 1)));*/

	sprintf(fitshdrcomment, "Shift for variable number %d LCgroup %d REFID %d", vv, ll, kk + stitch->minrefnumindx[lc_num]);
	VARTOOLS_Add_Keyword_To_OutputLC_FitsHeader(p, lc_num, fitshdrkeyword, fitshdrcomment, stitch->add_shifts_fitsheader_hdutouse, stitch->add_shifts_fitsheader_updateexistingkeyword, VARTOOLS_TYPE_DOUBLE, -Avector[ndecorrfromharm + i - 1]);
      }
    }
  }

  stitch->Nlcgroups_used[lc_num] = NLCgroups;
  stitch->Ntimegroups_used[lc_num] = ntgused;
  stitch->Nparamtotal_used[lc_num] = ndecorr;

  /* Free allocated memory, and return  */
  if(data_vals != NULL) free(data_vals);
  if(err_vals != NULL) free(err_vals);

  for(i = 0; i < npoints; i++) {
    if(decorr_terms[i] != NULL) free(decorr_terms[i]);
  }

  if(decorr_terms != NULL) free(decorr_terms);
  if(order_terms != NULL) free(order_terms);
  if(Avector != NULL) free(Avector);
  if(A_errvector != NULL) free(A_errvector);
  if(nextrefidchar2touse != NULL) free(nextrefidchar2touse);
  if(nextrefidchar1touse != NULL) free(nextrefidchar1touse);

  return;
}


void DoStitch(ProgramData *p, _Stitch *stitch, int lc_name_num, 
	      int lc_num)
{
  double *magvals;
  double *errvals;
  double *shiftvals;
  int *maskvals;
  int *groupidx;
  int Ngroup;
  int NLCgroups = 0;
  _StitchLightCurveGroup *lcg = NULL;
  int Ntimegroups = 0;
  int Nusedtimegroups = 0;
  _StitchTimeGroup *tg = NULL;
  int vv, testredogroups, isdiff, i, j;
  char coeffoutname[MAXLEN];
  FILE *coeffoutfile = NULL;
  int inshift_starid = -1;
  int *fieldlabels_indx = NULL;
  int *fieldlabels_needing_update = NULL;
  int *fieldlabels_obscounts = NULL;
  int sizefieldlabels_obscounts = 0;
  int sizefieldlabels_needing_update = 0;
  int Nfieldlabels_needing_update = 0;
  int Nfieldlabels_can_update = 0;
  int all_points_inshift = 0;

  int N_shift_out;
  char **out_shift_labels = NULL;
  int *size_out_shift_labels = NULL;
  char **out_shift_values_str = NULL;
  int *size_out_shift_values_str = NULL;
  int *Nobs_out_shift = NULL;
  int size_out_shift = 0;

  /**** setup the coefficient output file ************/
  if(stitch->save_fitted_parameters) {
    VARTOOLS_GetOutputFilename(coeffoutname, p->lcnames[lc_name_num],
			       stitch->coeffoutdir, "stitch",
			       stitch->coeffoutformat, lc_name_num);
    if((coeffoutfile = fopen(coeffoutname,"w")) == NULL) {
      fprintf(stderr,"Cannot write to %s\n", coeffoutname);
      exit(ERR_CANNOTWRITE);
    }
  }

  if(stitch->add_stitchparams_fitsheader)
    StitchAddStitchParamsToHeader(p, stitch, lc_name_num, lc_num);

  /***** If we need to append the refid to the field label, do it now ***/
  if(stitch->is_in_shifts_file || stitch->is_out_shifts_file) {
    if(stitch->userefnum && stitch->is_append_refnum_to_fieldlabel) {
      for(i=0; i < p->NJD[lc_num]; i++) {
	sprintf(&(stitch->field_labels_vals[lc_num][i][strlen(stitch->field_labels_vals[lc_num][i])]),"/%d",stitch->refnumval[lc_num][i]);
      }
    }
  }

  /***** Cycle through the variables to stitch *****/
  for(vv = 0; vv < stitch->nstitchvar; vv++) {

    all_points_inshift = 0;

    /***** If we are using an input shifts file, check if the star is in the file ****/
    if(stitch->is_in_shifts_file) {
      if(p->NJD[lc_num] > 0) {
	inshift_starid = stitch_Find_Star_in_InshiftFile(stitch, p, lc_num, lc_name_num, vv);
	if(inshift_starid >= 0) {
	  if(fieldlabels_indx == NULL) {
	    if((fieldlabels_indx = (int *) malloc(p->NJD[lc_num]*sizeof(int))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	  }
	  stitch_Find_FieldLabels_in_InshiftFile(stitch, p, lc_num, lc_name_num, vv, 
						   fieldlabels_indx, inshift_starid);
	  /* Check for any field labels in the input shift file that have too few observations,
             and see if the new light curve has enough observations to update them */
	  if(stitch->is_nobs_refit) {
	    if(!sizefieldlabels_needing_update) {
	      sizefieldlabels_needing_update = stitch->N_shifts_per_star[vv][inshift_starid];
	      if(sizefieldlabels_needing_update > 0) {
		if((fieldlabels_needing_update = (int *) malloc(sizefieldlabels_needing_update * sizeof(int))) == NULL || 
		   (fieldlabels_obscounts = (int *) malloc(sizefieldlabels_needing_update * sizeof(int))) == NULL)
		  VARTOOLS_error(ERR_MEMALLOC);
	      }
	    } else if(stitch->N_shifts_per_star[vv][inshift_starid] > sizefieldlabels_needing_update) {
	      sizefieldlabels_needing_update = stitch->N_shifts_per_star[vv][inshift_starid];
	      if((fieldlabels_needing_update = (int *) realloc(fieldlabels_needing_update, sizefieldlabels_needing_update * sizeof(int))) == NULL ||
		 (fieldlabels_obscounts = (int *) realloc(fieldlabels_obscounts, sizefieldlabels_needing_update * sizeof(int))) == NULL)
		VARTOOLS_error(ERR_MEMALLOC);
	    }
	    Nfieldlabels_needing_update = 0;
	    for(j=0; j < stitch->N_shifts_per_star[vv][inshift_starid]; j++) {
	      if(stitch->Nobs_in_shifts[vv][inshift_starid][j] < stitch->nobs_refit) {
		fieldlabels_needing_update[j] = 1;
		Nfieldlabels_needing_update++;
	      } else {
		fieldlabels_needing_update[j] = 0;
	      }
	    }
	    if(Nfieldlabels_needing_update > 0) {
	      /* Check if we have enough observations for any of the fields needing updates
		 to be able to perform the updates */
	      for(j=0; j < stitch->N_shifts_per_star[vv][inshift_starid]; j++) {
		fieldlabels_obscounts[j] = 0;
	      }
	      for(i=0; i < p->NJD[lc_num]; i++) {
		if(stitch->stitchmaskvals[vv][lc_num][i] <= VARTOOLS_MASK_TINY)
		  continue;
		if(fieldlabels_indx[i] >= 0) {
		  fieldlabels_obscounts[fieldlabels_indx[i]] += 1;
		}
	      }
	      for(i = 0; i < p->NJD[lc_num]; i++) {
		if(fieldlabels_indx[i] >= 0) {
		  if(fieldlabels_needing_update[fieldlabels_indx[i]] &&
		     fieldlabels_obscounts[fieldlabels_indx[i]] >= stitch->nobs_refit) {
		    /* We have enough observations of this field to update the fit; so
                       mark this observation as having no fields */
		    fieldlabels_indx[i] = -1;
		  }
		}
	      }
	    }
	  }
	
	  /* Apply the shifts; all_points_inshift will be 1 if all of the points have an
             input shift to apply */
	  all_points_inshift = stitch_Apply_InShifts(stitch, p, lc_num, lc_name_num, vv, fieldlabels_indx, inshift_starid);
	  StitchAddInShiftsToHeader(p, stitch,
				    lc_name_num, lc_num, NLCgroups,
				    lcg, Ntimegroups,
				    Nusedtimegroups, tg,
				    vv, fieldlabels_indx, inshift_starid);
	  if(all_points_inshift) {
	    SetStitchResultsNothingToDo(p, stitch, lc_name_num, lc_num, NLCgroups,
					  lcg, Ntimegroups, Nusedtimegroups, tg, vv);
	    if(stitch->is_out_shifts_file) {
	      /* Prepare the data for the output shifts file */
	      PrepareOutShifts(p, stitch, lc_name_num, lc_num, NLCgroups, lcg, Ntimegroups, Nusedtimegroups, tg, vv, fieldlabels_indx, inshift_starid, all_points_inshift, fieldlabels_needing_update, fieldlabels_obscounts, &N_shift_out, &size_out_shift, &out_shift_labels, &size_out_shift_labels, &out_shift_values_str, &size_out_shift_values_str, &Nobs_out_shift);
	      stitch_WriteOutShifts_Line(stitch, p, lc_num, lc_name_num, vv, N_shift_out, size_out_shift, out_shift_labels, size_out_shift_labels, out_shift_values_str, size_out_shift_values_str, Nobs_out_shift);
	    }
	    continue;
	  }
	}
      }
    }

    if(!vv || stitch->is_in_shifts_file) testredogroups = 1;
    else {
      /***** Check if we need to do the grouping - 
             this is necessary if we are working with the first variable,
             or if there are any different mask values from the last variable
             that we worked on ******/
      testredogroups = 0;
      for(i=0; i < p->NJD[lc_num]; i++) {
	if((stitch->stitchmaskvals[vv][lc_num][i] <= VARTOOLS_MASK_TINY && stitch->stitchmaskvals[vv-1][lc_num][i] > VARTOOLS_MASK_TINY) ||
	   (stitch->stitchmaskvals[vv][lc_num][i] > VARTOOLS_MASK_TINY && stitch->stitchmaskvals[vv-1][lc_num][i] <= VARTOOLS_MASK_TINY)) {
	  testredogroups = 1;
	  break;
	}
      }
    }
    if(testredogroups) {
      if(vv > 0) {
	FreeLCgroups(&NLCgroups, &lcg);
	Freetimegroups(&Ntimegroups, &tg);
      }
      /***** Form the light curve groups, then the time groups ****/
      if(stitch->is_in_shifts_file && p->NJD[lc_num] > 0 && inshift_starid >= 0) {
	FormLCGroups_IsInShift(p, stitch, lc_name_num, 
			       lc_num, &NLCgroups, &lcg, vv, fieldlabels_indx);
      } else {
	FormLCGroups(p, stitch, lc_name_num, lc_num, &NLCgroups, &lcg, vv);
      }
      if(NLCgroups > 0)
	if(stitch->groupbytime) {
	  FormTimeGroups(p, stitch, lc_name_num, lc_num, NLCgroups, lcg,
			 &Ntimegroups, &tg, &Nusedtimegroups, vv);
	} else {
	  FormOnlyOneTimeGroup(p, stitch, lc_name_num, lc_num, NLCgroups, lcg,
			       &Ntimegroups, &tg, &Nusedtimegroups, vv);
	}
    }    
  
    /* Now do the stitching */
    switch(stitch->stitchmethod) {
    case VARTOOLS_STITCH_METHOD_MEDIAN:
    case VARTOOLS_STITCH_METHOD_MEAN:
    case VARTOOLS_STITCH_METHOD_WEIGHTEDMEAN:
      StitchByStatistic(p, stitch, lc_name_num, lc_num, NLCgroups, lcg, Ntimegroups, Nusedtimegroups, tg, vv, coeffoutfile);
      break;
    case VARTOOLS_STITCH_METHOD_POLY:
      StitchByPoly(p, stitch, lc_name_num, lc_num, NLCgroups, lcg, Ntimegroups, Nusedtimegroups, tg, vv, coeffoutfile);
      break;
    case VARTOOLS_STITCH_METHOD_HARM:
      StitchByHarm(p, stitch, lc_name_num, lc_num, NLCgroups, lcg, Ntimegroups, Nusedtimegroups, tg, vv, coeffoutfile);
      break;
    default:
      VARTOOLS_error(ERR_CODEERROR);
      break;
    }

    if(stitch->is_out_shifts_file) {
      /* Prepare the data for the output shifts file */
      PrepareOutShifts(p, stitch, lc_name_num, lc_num, NLCgroups, lcg, Ntimegroups, Nusedtimegroups, tg, vv, fieldlabels_indx, inshift_starid, all_points_inshift, fieldlabels_needing_update, fieldlabels_obscounts, &N_shift_out, &size_out_shift, &out_shift_labels, &size_out_shift_labels, &out_shift_values_str, &size_out_shift_values_str, &Nobs_out_shift);
      stitch_WriteOutShifts_Line(stitch, p, lc_num, lc_name_num, vv, N_shift_out, size_out_shift, out_shift_labels, size_out_shift_labels, out_shift_values_str, size_out_shift_values_str, Nobs_out_shift);
    }
  }

  if(out_shift_labels != NULL) {
    for(i = 0; i < size_out_shift; i++) {
      if(out_shift_labels[i] != NULL)
	free(out_shift_labels[i]);
    }
    free(out_shift_labels);
  }
  if(size_out_shift_labels != NULL) free(size_out_shift_labels);
  if(out_shift_values_str != NULL) {
    for(i = 0; i < size_out_shift; i++) {
      if(out_shift_values_str[i] != NULL)
	free(out_shift_values_str[i]);
    }
    free(out_shift_values_str);
  }
  if(size_out_shift_values_str != NULL) free(size_out_shift_values_str);
  if(Nobs_out_shift != NULL) free(Nobs_out_shift);

  if(fieldlabels_indx != NULL)
    free(fieldlabels_indx);
  if(fieldlabels_needing_update != NULL)
    free(fieldlabels_needing_update);
  if(fieldlabels_obscounts != NULL)
    free(fieldlabels_obscounts);


  if(coeffoutfile != NULL)
    fclose(coeffoutfile);
  FreeLCgroups(&NLCgroups, &lcg);
  Freetimegroups(&Ntimegroups, &tg);
  
}


