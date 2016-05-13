#include "../../src/vartools.h"
#include <stdio.h>
#include <stdlib.h>
#include "medlc.h"

/* This is the source code for a sample user-defined command
   to be used with vartools.

   This library defines the command -medlc which takes multiple
   light curves with the same time sampling, and computes the
   median light curve from them, and outputs the result to a
   file. This example illustrates a command which processes all
   of the light curves at once.
*/

void medlc_Initialize(char *commandname,
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
  sprintf(commandname,"-medlc");

  /* This is 1 if the command requires all light curves to be read
     at once, otherwise set it to 0. Here we set it to 1. */
  *RequireReadAll = 1;

  /* This is 1 if the command requires input light curves to be sorted by
     time, otherwise set it to 0. */
  *RequireSortLC = 1;
  
  /* This is 1 if the command requires the times of observation in the
     input light curve to be unique */
  *RequireDistinctTimes = 0;

  /* You should define a structure to store the data needed for this
     command (in this example, a vector to hold the values to add to
     the magnitudes of each lc), below you would replace "_Medlc"
     with the type name of your structure.

     See medlc.h
*/
  *sizeuserdata = sizeof(_Medlc);
}

int medlc_ParseCL(ProgramData *p, Command *c,
			   void *userdata, int *iret, char **argv, int argc)
/* Every library with the name $LIBNAME.so must contain a function
   with the name $LIBNAME_ParseCL to check that the user has issued a
   valid syntax, and set-up the command-data accordingly.

   Note that any vectors/arrays used by the command should also be registered
   by this function so that memory will be allocated as needed.

Expected syntax for this example:
   -medlc "output_file_name"

  p = structure containing general program data. In most cases this 
      structure is only passed to other functions (like RegisterDataVector)
      and not used directly by the ParseCL function.

  c = pointer to a generic container for this command. This structure is
      needed by other functions that might be called here (RegisterDataVector)
      but should not be used directly by the ParseCL function.

  userdata = pointer to the structure storing data for this
      command. This should be cast to the appropriate type and the
      data there-in should be modified as needed by the ParseCL
      function. In this example, userdata corresponds to the _Medlc
      structure which stores the values to be added to the light
      curves.

  iret = command line argument index. This is initially zero, and on
      return should be set to the index of the last argument parsed by
      this command. Note that index 0 points to the name of the
      command, which does not need to be checked (this function is
      only called if the user issues the command). In this example, if
      the user issues the command "-medlc fix 5." the return value
      of iret should be 2.

  argv = array of command line arguments. argv[0] is the name of the
      command ("-medlc" in this case).

  argc = Number of command line arguments in the argv array (including
      0). It is the user's responsibility to check that i < argc
      before attempting to parse argv[i], failure to do so may lead to
      segmentation violations.
 */
{
  int check;

  /* Cast the input userdata pointer to the data structure type used
     for this command */
  _Medlc *Medlc;
  Medlc = (_Medlc *) userdata;

  /* The expected command syntax is "-medlc outfilename", so we only need to
     scan a single constant parameter, the output file name. This can be done
     with the VARTOOLS_ParseConstantParameter function */

  check = VARTOOLS_ParseConstantParameter(p, c, iret, argv, argc, NULL,
					  VARTOOLS_TYPE_STRING,
					  (void *) (&Medlc->outfilename),
					  0);

  /* Check = 1 or 2 if there was an error in parsing the line, return
     a positive value if there is to tell vartools to print out the 
     syntax and quit */
  if(check)
    return 1;
  else
    return 0;
}

void medlc_ShowSyntax(FILE *outfile)
/* Write out the expected syntax for this command to the file outfile */
{
  fprintf(outfile,
	  "-medlc outfilename\n");
}

void medlc_ShowHelp(FILE *outfile)
/* Output the help information for this command */
{
  /* Give the verbose help description */
  fprintf(outfile,"Computes a median light curve from all of the input light curves, and outputs the result to the file \"outfilename\". This command requires all of the light curves to be read-in. Also, all of the light curves must be the same length.\n\n");
}

void medlc_ShowExample(FILE *outfile)
/* Output an example for this command */
{
  fprintf(outfile,"echo EXAMPLES/1 | \\\n"
	  "   gawk '{for(i=1; i <= 10; i += 1) print $1;}' | \\\n"
	  "   ./vartools -l - -L USERLIB/src/medlc.so \\\n"
          "      -expr 'mag=10.0+err*gauss()' \\\n"
          "      -o EXAMPLES/OUTDIR1/ nameformat \"sim%d.txt\" \\\n"
          "      -medlc EXAMPLES/OUTDIR1/testmedlc.txt \\\n"
          "      -quiet\n\n"
          "Example use of the -medlc command. We read in the light curve EXAMPLES/1 ten times (this is done with the initial gawk command). Then in each case the light curve magnitudes are replaced with gaussian random noise, simply for the purpose of simulating 10 noisy light curves with the same time sampling. The simulated light curves are output to EXAMPLES/OUTDIR1/sim1.txt ... EXAMPLES/OUTDIR1/sim10.txt. We then compute the median of these simulated light curves, with the result written to EXAMPLES/OUTDIR1/testmedlc.txt.\n\n");
}

void medlc_RunCommand(ProgramData *p, void *userdata, int lc_name_num, int lc_num)
/* This function runs the command on a light curve.  

   p = structure containing various general program data (in
       particular the light curves are contained in this structure).

   userdata = pointer to the structure containing the command specific
       data (including control parameters and vectors to store output
       results).

   lc_name_num = this is "-1" and can be ignored in this case because we have
   set RequireReadAll to 1.

   lc_num = this is "-1" and can be ignored in this case because we have set
   RequireReadAll to 1.

*/
{
  int i, j, NJD;
  int Nlcs;
  double **t, **mag, **err;

  double *magmed, *allmags;

  FILE *outfile;
  _Medlc *Medlc;

  /* Translate pointers and structures input by vartools into
     easier to use forms */
  Medlc = (_Medlc *)userdata;


  /* Get the number of light curves to process */
  Nlcs = p->Nlcs;

  /* Check that all of the light curves are the same length, quit with an 
     error if they are not */
  NJD = p->NJD[0];
  for(i=1; i < Nlcs; i++) {
    if(p->NJD[i] != NJD) {
      fprintf(stderr,"Error: all of the light curves must be of the same length to use the -medlc command.\n");
      exit(1);
    }
  }

  /* Open the file to output the median lc to */
  if((outfile = fopen(Medlc->outfilename,"w")) == NULL)
    {
      fprintf(stderr,"Error: cannot write to %s\n", Medlc->outfilename);
      exit(1);
    }

  /* The following are matrices holding the times, magnitudes and
     uncertainties */
  t = p->t; mag = p->mag; err = p->sig;

  /* Create temporary arrays to store the median light curve, as well as
     the set of points to calculate the median for, we will just do a
     direct call to malloc to allocate this memory. Remember to free these
     at the end of the function. */

  if((magmed = (double *) malloc(NJD * sizeof(double))) == NULL ||
     (allmags = (double *) malloc(Nlcs * sizeof(double))) == NULL)
    {
      fprintf(stderr,"Error: Memory allocation error in -medlc.\n");
      exit(1);
    }
  
  /* Compute the median magnitude at each time */
  for(j=0; j < NJD; j++) {
    for(i=0; i < Nlcs; i++)
      allmags[i] = mag[i][j];
    magmed[j] = VARTOOLS_median(Nlcs, allmags);
  }

  /* Output the results */

  for(j = 0; j < NJD; j++) {
    fprintf(outfile,"%.17g %.17g\n", t[0][j], magmed[j]);
  }

  /* clean up */
  fclose(outfile);
  free(magmed);
  free(allmags);  
}

