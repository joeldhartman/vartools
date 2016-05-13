#include "../../src/vartools.h"
#include <stdio.h>
#include <stdlib.h>
#include "magadd.h"

/* This is the source code for a sample user-defined command
   to be used with vartools.

   This library defines the command -magadd which can be used
   to add a constant to the magnitude values of a light curve 
*/

void magadd_Initialize(char *commandname,
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
  sprintf(commandname,"-magadd");

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
  *sizeuserdata = sizeof(_Magadd);
}

int magadd_ParseCL(ProgramData *p, Command *c,
			   void *userdata, int *iret, char **argv, int argc)
/* Every library with the name $LIBNAME.so must contain a function
   with the name $LIBNAME_ParseCL to check that the user has issued a
   valid syntax, and set-up the command-data accordingly.

   Note that any vectors/arrays used by the command should also be registered
   by this function so that memory will be allocated as needed.

Expected syntax for this example:
   -magadd <"fix" value | "list" | "fixcolumn" <colname | colnum> | "expr" expr>

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
  _Magadd *Magadd;
  Magadd = (_Magadd *) userdata;

  /* We can use the function ParseFixSpecFixcolumn to parse a line of the form
     <"fix" value | "list" [\"column\" col] | "fixcolumn" <colname | colnum> | \"expr\" expr>

     After giving p, c, iret, argv, and argc, give the number of
     data-vectors to be set by this command. For each data-vector give
     the data-type for the vector and a pointer to the vector.

     The function will update iret.
  */
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, iret, argv, argc, 1,
				VARTOOLS_TYPE_DOUBLE, 
				(void *) (&Magadd->addval),
				0,
				1,
				"addval");

  /* Check = 1 if there was an error in parsing the line, return
     a positive value if there is to tell vartools to print out the 
     syntax and quit */
  if(check)
    return 1;
  else
    return 0;
}

void magadd_ShowSyntax(FILE *outfile)
/* Write out the expected syntax for this command to the file outfile */
{
  fprintf(outfile,
	  "-magadd <\"fix\" value | \"list\" [\"column\" col] | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n");
}

void magadd_ShowHelp(FILE *outfile)
/* Output the help information for this command */
{
  /* Give the verbose help description */
  fprintf(outfile,"Add a constant to magnitudes of a light curve. If \"fix\" is given then the constant specified on the command line will be added for all light curves. If \"list\" is given, the constant will be read-in from the input-list file for each star, use the \"column\" keyword to optionally specify the column number to read from the input list. If \"fixcolumn\" is given then the constant will be taken from a previously computed output statistic. If \"expr\" is given then the constant will be determined by evaluating an analytic expression.\n\n");
}

void magadd_RunCommand(ProgramData *p, void *userdata, int lc_name_num, int lc_num)
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
  int NJD;
  double *t, *mag, *err;
  char *lcname;
  double addval;
  _Magadd *Magadd;

  void domagadd(int N, double *mag, double addval);

  /* Translate pointers and structures input by vartools into
     easier to use forms */
  Magadd = (_Magadd *)userdata;


  /* The number of points in the light curve, and vectors storing the
     time, magnitudes, and errors are accessed as illustrated
     below. Note that t, mag and err are vectors, t[0] is the time of
     the first observation in the light curve, t[NJD-1] is the time of
     the last observation, etc. */
  NJD=p->NJD[lc_num];
  t = p->t[lc_num];
  mag = p->mag[lc_num];
  err = p->sig[lc_num];

  /* This command does not actually use the name of the light curve, 
     but if you need it for some reason, you access it as follows: */
  lcname = p->lcnames[lc_name_num];

  /* Get the value to add to this light curve, the contents of the
     Magadd->addval vector will be set by VARTOOLS depending on how the
     user calls the command on the command-line. */
  addval = Magadd->addval[lc_num];

  /* Perform the actual routine */
  domagadd(NJD, mag, addval);
  
}

/* This is the function to do the actual work for this command.
   In this case the function is trivial, it just adds a constant
   to the magnitudes in a light curve. */
void domagadd(int N, double *mag, double addval)
{
  int i;
  for(i=0; i < N; i++)
    mag[i] += addval;
}
