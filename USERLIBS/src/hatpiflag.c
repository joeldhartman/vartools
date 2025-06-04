#include "../../src/vartools.h"
#include <stdio.h>
#include <stdlib.h>
#include "hatpiflag.h"

/* This is the source code for a sample user-defined command
   to be used with vartools.

   This library defines the command -magadd which can be used
   to add a constant to the magnitude values of a light curve 
*/

void hatpiflag_Initialize(char *commandname,
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
  sprintf(commandname,"-hatpiflag");

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
     the magnitudes of each lc), below you would replace "_Magadd"
     with the type name of your structure.

     See magadd.h
*/
  *sizeuserdata = sizeof(_Hatpiflag);
}

int hatpiflag_ParseCL(ProgramData *p, Command *c,
			   void *userdata, int *iret, char **argv, int argc)
/* Every library with the name $LIBNAME.so must contain a function
   with the name $LIBNAME_ParseCL to check that the user has issued a
   valid syntax, and set-up the command-data accordingly.

   Note that any vectors/arrays used by the command should also be registered
   by this function so that memory will be allocated as needed.

Expected syntax for this example:
   -hatpiflag <"fix" value | "list" | "fixcolumn" <colname | colnum> | "expr" expr>

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
  int check;

  /* Cast the input userdata pointer to the data structure type used
     for this command */
  _Hatpiflag *hf;
  hf = (_Hatpiflag *) userdata;

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
  
  VARTOOLS_RegisterDataVector(p, c, (void *) (&(hf->fiphotflagvals)),
			      VARTOOLS_TYPE_STRING, 0,
			      VARTOOLS_SOURCE_EXISTINGVARIABLE,
			      0, NULL, argv[i],
			      (char) VARTOOLS_VECTORTYPE_LC,
			      &(hf->fiphotflagvar));
  
  i++;
  if(i >= argc) {
    return 1;
  }
  
  VARTOOLS_RegisterDataVector(p, c, (void *) (&(hf->rejbadframemaskvals)),
			      VARTOOLS_TYPE_DOUBLE, 0,
			      VARTOOLS_SOURCE_EXISTINGVARIABLE,
			      0, NULL, argv[i],
			      (char) VARTOOLS_VECTORTYPE_LC,
			      &(hf->rejbadframemaskvar));

  i++;
  if(i >= argc) {
    return 1;
  }
  
  VARTOOLS_RegisterDataVector(p, c, (void *) (&(hf->tfaoutliermaskvals)),
			      VARTOOLS_TYPE_DOUBLE, 0,
			      VARTOOLS_SOURCE_EXISTINGVARIABLE,
			      0, NULL, argv[i],
			      (char) VARTOOLS_VECTORTYPE_LC,
			      &(hf->tfaoutliermaskvar));

  i++;
  if(i >= argc) {
    return 1;
  }

  VARTOOLS_RegisterDataVector(p, c, (void *) (&(hf->pointingoutlierflagvals)),
			      VARTOOLS_TYPE_DOUBLE, 0,
			      VARTOOLS_SOURCE_EXISTINGVARIABLE,
			      0, NULL, argv[i],
			      (char) VARTOOLS_VECTORTYPE_LC,
			      &(hf->pointingoutlierflagvar));

  i++;
  if(i >= argc) {
    return 1;
  }

  
  VARTOOLS_RegisterDataVector(p, c, (void *) (&(hf->outflagvals)),
			      VARTOOLS_TYPE_SHORT, 0,
			      VARTOOLS_SOURCE_LC,
			      0, NULL, NULL, -1, NULL,
			      argv[i]);

  i++;

  *iret = i;
  return 0;
}

void hatpiflag_ShowSyntax(FILE *outfile)
/* Write out the expected syntax for this command to the file outfile */
{
  fprintf(outfile,
	  "-hatpiflag fiphot_string_flag_var rejbadframe_mask_var TFA_outlier_mask_var pointing_outlier_flag_var output_flag_var\n");
}

void hatpiflag_ShowHelp(FILE *outfile)
/* Output the help information for this command */
{
  /* Give the verbose help description */
  fprintf(outfile,"Produce the output binary flag variable for the HATPI light curves. Takes as input a vector variable storing the string flags from fiphot, a variable used for masking rejected frames (0 = rejected, 1 = not rejected), a variable used for masking TFA outliers to exclude from the fit (0 = outlier, 1 = not outlier), and a variable used for flagging pointing outliers (1 = outlier, 0 = not an outlier). The resulting binary flag will be output to the output_flag_var variable.\n\n");
}

void hatpiflag_RunCommand(ProgramData *p, void *userdata, int lc_name_num, int lc_num)
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
  int NJD, i;
  _Hatpiflag *hf;


  /* Translate pointers and structures input by vartools into
     easier to use forms */
  hf = (_Hatpiflag *) userdata;


  /* The number of points in the light curve, and vectors storing the
     time, magnitudes, and errors are accessed as illustrated
     below. Note that t, mag and err are vectors, t[0] is the time of
     the first observation in the light curve, t[NJD-1] is the time of
     the last observation, etc. */
  NJD=p->NJD[lc_num];

  /* First set the bits from the fiphot string flag */
  for(i = 0; i < NJD; i++) {
    if(!strcmp(hf->fiphotflagvals[lc_num][i],"G")) {
      hf->outflagvals[lc_num][i] = 0; /* Good photometry - no bits set */
    }
    else if(!strcmp(hf->fiphotflagvals[lc_num][i],"X")) {
      hf->outflagvals[lc_num][i] = 1; /* Bad photometry - first bit set */
    }
    else if(!strcmp(hf->fiphotflagvals[lc_num][i],"C")) {
      hf->outflagvals[lc_num][i] = 2; /* Saturated or hot pixels in aperture - second bit */
    }
    else if(!strcmp(hf->fiphotflagvals[lc_num][i],"A")) {
      hf->outflagvals[lc_num][i] = 4; /* Asteroid in aperture - third bit */
    }
    else if(!strcmp(hf->fiphotflagvals[lc_num][i],"S")) {
      hf->outflagvals[lc_num][i] = 8; /* Satellite in aperture - fourth bit */
    }
    else if(!strcmp(hf->fiphotflagvals[lc_num][i],"H")) {
      hf->outflagvals[lc_num][i] = 4 | 8; /* Asteroid and Satellite - third and fourth bits */
    }
    else if(!strcmp(hf->fiphotflagvals[lc_num][i],"I")) {
      hf->outflagvals[lc_num][i] = 2 | 8; /* Saturated/hot pixel(s) and satellite - second and fourth bits */
    }
    else if(!strcmp(hf->fiphotflagvals[lc_num][i],"J")) {
      hf->outflagvals[lc_num][i] = 2 | 4; /* Saturated/hot pixel(s) and asteroid - second and third bits */
    }
    else if(!strcmp(hf->fiphotflagvals[lc_num][i],"K")) {
      hf->outflagvals[lc_num][i] = 2 | 4 | 8; /* Saturated/hot pixel(s) and asteroid and satellite - second and third and fourth bits */
    }
    else 
      hf->outflagvals[lc_num][i] = 0;
  
    if(hf->rejbadframemaskvals[lc_num][i] <= VARTOOLS_MASK_TINY) {
      hf->outflagvals[lc_num][i] = hf->outflagvals[lc_num][i] | 16; /* frames marked as having an
                                             excessive number of outlier points
                                             amongst all light curves have the
                                             fifth bit set */
    }
    if(hf->tfaoutliermaskvals[lc_num][i] <= VARTOOLS_MASK_TINY) {
      hf->outflagvals[lc_num][i] = hf->outflagvals[lc_num][i] | 32; /* frames excluded from the TFA fit
                                             due to the points being outliers
                                             have the sixth bit set */  
    }

    if(hf->pointingoutlierflagvals[lc_num][i] > VARTOOLS_MASK_TINY) {
      hf->outflagvals[lc_num][i] = hf->outflagvals[lc_num][i] | 64; /* frames flagged as having bad pointing have the seventh bit set */
    }
  }

}

