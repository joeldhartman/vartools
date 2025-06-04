#include "../../src/vartools.h"
#include <stdio.h>
#include <stdlib.h>
#include "splinedetrend.h"
#include <gsl/gsl_bspline.h>

/* This is the source code for a sample user-defined command
   to be used with vartools.

   This library defines the command -splinedetrend which can be used
   to perform a multi-dimensional B-spline detrending of a light curve.

   This is a fairly involved example which, among other things,
   illustrates the inclusion of FORTRAN routines in a VARTOOLS
   library.

*/

int ParseSplineDetrendVarlist(char *varlist,  _Splinedetrend *Splinedetrend, ProgramData *p, Command *c);

int ParseSplineDetrendOutModelVar(char *varlist,  _Splinedetrend *Splinedetrend, ProgramData *p, Command *c);

void SplineDetrend_FillSplineTerms(int *ll, double **decorr_terms, int Np, double *datavec, double minval, double maxval, int nbreak, int splineorder, int *order_terms);

void SplineDetrend_FillPolyTerms(int *ll, double **decorr_terms, int Np, double *datavec, double minval, double maxval, int polyorder, int *order_terms);

void SplineDetrend_FillHarmTerms(int *ll, double **decorr_terms, int Np, double *datavec, double minval, double maxval, int nharm, int *order_terms);

void splinedetrend_Initialize(char *commandname,
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
  sprintf(commandname,"-splinedetrend");

  /* This is 1 if the command requires all light curves to be read
     at once, otherwise set it to 0. */
  *RequireReadAll = 0;

  /* This is 1 if the command requires input light curves to be sorted by
     time, otherwise set it to 0. */
  *RequireSortLC = 1;
  
  /* This is 1 if the command requires the times of observation in the
     input light curve to be unique */
  *RequireDistinctTimes = 1;

  /* You should define a structure to store the data needed for this
     command (in this example, a vector to hold the values to add to
     the magnitudes of each lc), below you would replace "_Jktebop"
     with the type name of your structure. This tells VARTOOLS how
     much memory to allocate for the structure storing the data for
     your command.

     See splinedetrend.h
*/
  *sizeuserdata = sizeof(_Splinedetrend);
}

int splinedetrend_ParseCL(ProgramData *p, Command *c,
			  void *userdata, int *iret, char **argv, int argc)
/* Every library with the name $LIBNAME.so must contain a function
   with the name $LIBNAME_ParseCL to check that the user has issued a
   valid syntax, and set-up the command-data accordingly.

   Note that any vectors/arrays used by the command should also be registered
   by this function so that memory will be allocated as needed.

Expected syntax for this example:
    -splinedetrend detrendvec1:<"spline":knotspacing:order|"poly":order|"harm":nharm>[:"groupbygap":gapsize][,detrendvec2:<"spline":knotspacing:order|"poly":order|"harm":nharm>[:"groupbygap":gapsize],...] \
        ["sigmaclip" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>] \
 	["omodel" outdir ["nameformat" format]]
 	["omodelcoeffs" outdir ["nameformat" format]]
	["omodelvariable" outvarname1[:inputvarsignal1[,outvarname2:inputvarsignal2...]]]

  p = structure containing general program data. In most cases this 
      structure is only passed to other functions (like RegisterDataVector)
      and not used directly by the ParseCL function.

  c = pointer to a generic container for this command. This structure is
      needed by other functions that might be called here (RegisterDataVector)
      but should not be used directly by the ParseCL function.

  userdata = pointer to the structure storing data for this
      command. This should be cast to the appropriate type and the
      data there-in should be modified as needed by the ParseCL
      function. In this example, userdata corresponds to the _Splinedetrend
      structure which stores the values to be added to the light
      curves.

  iret = command line argument index. Initially argv[*iret - 1] is the
      commandname ("-splinedetrend" in this case; this function will only be
      called if the user issues the command so this does not need to
      be verified by this function). On output *iret should be
      incremented by the number of terms parsed from the
      command-line.

  argv = array of command line arguments. argv[*iret - 1] is the name
      of the command ("-splinedetrend" in this case).

  argc = Number of command line arguments in the argv array (including
      0). It is the user's responsibility to check that i < argc
      before attempting to parse argv[i], failure to do so may lead to
      segmentation violations.
 */
{
  int i = 0, k=0, kstart = 0;
  char *varlist;
  int Nvar;
  int check;
  
  /* Cast the input userdata pointer to the data structure type used
     for this command */
  _Splinedetrend *Splinedetrend;
  Splinedetrend = (_Splinedetrend *) userdata;

  /* We'll use i rather than iret to index argv, this is just so we
     don't have to constantly dereference the pointer iret */
  i = *iret;

  /* Manually read the list of detrend variables and parse it */
  if(i >= argc) {
    /* There are no terms left on the command-line,
       return an error */
    return 1;
  }
  
  check = ParseSplineDetrendVarlist(argv[i],  Splinedetrend, p, c);
  if(check) {*iret = i; return 1;}

  i++;

  /* Now parse each of the optional parameters */

  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "sigmaclip", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Splinedetrend->sigmaclip), 0, 0,
				  "SigmaClip", 1, SPLINEDETREND_DEFAULT_SIGMACLIP);
  if(check == 2) {*iret = i; return 1;}


  check = VARTOOLS_ParseOutNameKeyword(p, c, &i, argv, argc, "omodel", 
				       &(Splinedetrend->outputmodel),
				       &(Splinedetrend->modeloutdir),
				       &(Splinedetrend->outputmodel_useformat),
				       &(Splinedetrend->outputmodel_format));
  if(check == 2) {*iret = i; return 1;}

  check = VARTOOLS_ParseOutNameKeyword(p, c, &i, argv, argc, "omodelcoeffs", 
				       &(Splinedetrend->outputmodelcoeffs),
				       &(Splinedetrend->modelcoeffsoutdir),
				       &(Splinedetrend->outputmodelcoeffs_useformat),
				       &(Splinedetrend->outputmodelcoeffs_format));
  if(check == 2) {*iret = i; return 1;}


  /* Check if the user is attempting to store the output model in a variable */
  Splinedetrend->Noutvar = 0;
  if(i < argc ? !strcmp(argv[i],"omodelvariable") : 0) {
    i++;
    if(i >= argc) {*iret = i; return 1;}
    check = ParseSplineDetrendOutModelVar(argv[i], Splinedetrend, p, c);
    if(check) {*iret = i; return 1;}
    i++;
  }

  /* Register the variables to store data that will appear in the output table */
  VARTOOLS_RegisterDataVector(p, c, (void *) (&Splinedetrend->magmedian),
			      VARTOOLS_TYPE_DOUBLE, 0, VARTOOLS_SOURCE_COMPUTED,
			      1, "MedianMagnitude");
  VARTOOLS_RegisterDataVector(p, c, (void *) (&Splinedetrend->Noutlier),
			      VARTOOLS_TYPE_INT, 0, VARTOOLS_SOURCE_COMPUTED,
			      1, "NOutliers");
  VARTOOLS_RegisterDataVector(p, c, (void *) (&Splinedetrend->Ngroups),
			      VARTOOLS_TYPE_INT, 0, VARTOOLS_SOURCE_COMPUTED,
			      1, "NDataGroups");
  VARTOOLS_RegisterDataVector(p, c, (void *) (&Splinedetrend->Nparamtotal),
			      VARTOOLS_TYPE_INT, 0, VARTOOLS_SOURCE_COMPUTED,
			      1, "NFitParamsTotal");

  *iret = i;
  return 0;
}

void splinedetrend_ShowSyntax(FILE *outfile)
/* Write out the expected syntax for this command to the file outfile */
{
  fprintf(outfile,
	  "-splinedetrend\n"
	  "\tdetrendvec1:<\"spline\":knotspacing:order|\"poly\":order|\"harm\":nharm>[:\"groupbygap\":gapsize][,detrendvec2:<\"spline\":knotspacing:order|\"poly\":order|\"harm\":nharm>[:\"groupbygap\":gapsize],...]\n"
	  "\t[\"sigmaclip\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>]\n"
	  "\t[\"omodel\" outdir [\"nameformat\" format]]\n"
	  "\t[\"omodelcoeffs\" outdir [\"nameformat\" format]]\n"
	  "\t[\"omodelvariable\" outvarname1[:inputvarsignal1[,outvarname2:inputvarsignal2...]]]\n");
}


void splinedetrend_ShowHelp(FILE *outfile)
/* Output the help information for this command */
{
  fprintf(outfile,"Perform a multivariable detrending of a light curve using linear basis functions. The light curve magnitudes will be fit as a linear combination of these basis functions, which are expressed in one or more auxiliary variables (e.g., t, x, y, etc.). The classes of functions supported include Basis splines, polynomials, and harmonics. Cross-terms between the detrending variables are not included in the model. The first argument to the command is a comma-separated list of variables to detrend the light curve magnitudes against, with additional keywords used to indicate what functional form will be used in the detrending, parameters controlling the function, and an option to split the fit at gaps in the variable. The arguments to supply are as follows:\n"
	  "\tdetrendvec - The name of the variable to detrend against\n"
	  "\t\t(e.g., \"t\" for detrending against time).\n"
	  "\n\tBasis type to use:\n"
	  "\n\t\t\"spline\" - Basis spline, computed using the GNU Scientific\n"
	  "\t\t           Library function gsl_bspline_eval.\n"
	  "\n\t\t\tknotspacing - the spacing to use between knots, in\n"
	  "\t\t\t              whatever units the variable is given in.\n"
	  "\n\t\t\torder - the order of the spline to use (e.g., 3 for\n"
	  "\t\t\t        a cubic spline).\n"
	  "\n\t\t\"poly\" - A polynomial in the detrendvec variable.\n"
	  "\n\t\t\torder - Order of the polynomial, e.g. 1 for a linear\n"
	  "\t\t\t        function, 2 for a quadratic, etc.\n"
	  "\n\t\t\"harm\" - A harmonic series in the detrendvec variable.\n"
	  "\n\t\t\tnharm - The number of harmonics to include. Use 0 to\n"
	  "\t\t\t        use only the fundamental (period equal to twice\n"
	  "\t\t\t        the range spanned by the variable), 1 to include\n"
	  "\t\t\t        also the first harmonic terms, 2 to include the\n"
	  "\t\t\t        second terms, etc.\n"
	  "\n\t\"groupbygap\" - An optional keyword to split the fit at gaps\n"
	  "\t               in the detrendvec variable.\n"
	  "\n\t\tgapsize - the maximimum space between consecutive points,\n"
	  "\t\t          above which a gap is defined, and the fit will be\n"
	  "\t\t          split into two groups.\n"
"\nAfter specifying the variables and functional form for the detrending, there are several optional keywords that can be supplied:\n"
	  "\n\t\"sigmaclip\" - Apply a sigma-clipping to the light curve\n"
	  "\t              to select points to exclude from the fit. The\n"
	  "\t              best-fit model will still be evaluated at these points\n"
	  "\t              and subtracted from the magnitude value.\n"
	  "\t              After giving the keyword, indicate the source of\n"
	  "\t              the sigma-clipping parameter to use. (use \"fix\"\n"
	  "\t              to fix it to a specific value for all light curves,\n"
	  "\t              \"list\" to specify the parameter in the input list,\n"
	  "\t              \"fixcolumn\" to take the value from the output of\n"
	  "\t              a previous command, or \"expr\" to use an analytic\n"
	  "\t              expression to initialize the parameter).\n"
	  "\n\t\"omodel\" - Output the best-fit model to a file. The file will\n"
	  "\t           be output to the directory outdir, and the name will be\n"
	  "\t           set to the input filename with the extension\n"
	  "\t           \"splinedetrend_model\" appended. Use the\n"
	  "\t           \"nameformat\" keyword to change the filename\n"
	  "\t           convention. See the help for the \"-o\" command for\n"
          "\t           the syntax.\n"
	  "\n\t\"omodelcoeffs\" - Output the best-fit coefficients to the\n"
	  "\t                 linear basis functions to a file. Options are\n"
          "\t                 as for the \"omodel\" keyword, but the default\n"
          "\t                 file extension is \"splinedetrend_modelcoeffs\".\n"
	  "\n\t\"omodelvariable\" - Store, in variables, the best-fit model\n"
	  "\t                   terms associated with one or more of the input\n"
	  "\t                   variables. Provide a comma-separated list with\n"
          "\t                   the name of the output variable followed by the\n"
          "\t                   name of the associated input variable used in\n"
          "\t                   the fit.\n");
}

void splinedetrend_ShowExample(FILE *outfile)
/* Output an example for this command */
{
  fprintf(outfile,
	  "\n> vartools -i EXAMPLES/6479535620075955328_llc.fits \\\n"
	  "\t-inputlcformat \\\n"
	  "\t\tt:TMID_BJD,mag:IRM1,err:IRE1,x:XIC,y:YIC,temp:CCDTEMP \\\n"
	  "\t-expr magorig=mag \\\n"
	  "\t-splinedetrend \\\n"
	  "\t\tt:spline:1.0:3:groupbygap:0.5,x:poly:1,y:poly:1,temp:poly:1 \\\n"
	  "\t\tsigmaclip fix 3.0 \\\n"
	  "\t\tomodel EXAMPLES/OUTDIR1/ \\\n"
	  "\t\tomodelcoeffs EXAMPLES/OUTDIR1/ \\\n"
	  "\t\tomodelvariable \\\n"
	  "\t\t\ttmod:t,xmod:x,ymod:y,tempmod:temp \\\n"
	  "\t-o EXAMPLES/OUTDIR1/6479535620075955328.splinedetrend.lc.txt \\\n"
	  "\t\tcolumnformat t,magorig,mag,err,x,y,temp,tmod,xmod,ymod,tempmod \\\n"
	  "\t-rms -oneline\n\n"
	  "Name                            = ./6479535620075955328_llc.fits\n"
          "Splinedetrend_MedianMagnitude_1 = 10.33609\n"
	  "Splinedetrend_NOutliers_1       = 0\n"
	  "Splinedetrend_NDataGroups_1     = 3\n"
	  "Splinedetrend_NFitParamsTotal_1 = 44\n"
	  "Mean_Mag_3                      =  10.33609\n"
	  "RMS_3                           =   0.00137\n"
	  "Expected_RMS_3                  =   0.00566\n"
	  "Npoints_3                       =  1149\n\n"
	  "Use the -splinedetrend command to detrend the a TESS sector one light curve for the star GAIA DR2 6479535620075955328. First we read the time, magnitude, uncertainty, x and y pixel position, and CCD temperature from the light curve. We then save the original magnitude values in the variable 'magorig.' We use the -splinedetrend command to perform the detrending. We detrend against the t, x, y, and temp variables. For t we use a third-order Basis spline, with a knot-spacing of 1.0 days. We split the fit into time groups at any gaps in time that are greater than 0.5 days. For x, y and temp we use a simple linear polynomial.  We apply 3 sigma clipping to the light curve before performing the fit. The best-fit model will be output to the file EXAMPLES/OUTDIR1/6479535620075955328_llc.fits.splinedetrend_model, and the optimized coefficients of the fit will be output to the file EXAMPLES/OUTDIR1/6479535620075955328_llc.fits.splinedetrend_modelcoeffs. The contribution of the t variable to the best-fit model will be stored in the variable tmod, the contribution from x will be stored in the xmod variable, from y in the ymod variable, and from temp in the tempmod variable. The detrended light curve will be output to the file EXAMPLES/OUTDIR1/6479535620075955328.splinedetrend.lc.txt, with the output including the time, the original magnitude, the detrended magnitude, the magnitude uncertainty, the x position, the y position, the CCD temperature, and the time, x, y and temperature contributions to the best-fit model. The four parameters given as output from the -splinedetrend command include: Splinedetrend_MedianMagnitude_1 - this is the median magnitude of the input light curve which is added back to the detrended magnitudes, after subtracted the best-fit detrending mode; Splinedetrend_NOutliers_1 - the number of outliers excluded from the fit; Splinedetrend_NDataGroups_1 - the number of separate fitting groups that the light curve is split to (as a result of the groupbygap option); Splinedetrend_NFitParamsTotal_1 - the total number of free parameters in the model that is fit to the light curve.\n");
}

void splinedetrend_RunCommand(ProgramData *p, void *userdata, int lc_name_num,
			      int lc_num)
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
  int Ngroup, *groupid = NULL, *groupidoutlier = NULL;
  int NJD, Ngood, i, k, j, ii, g, ll;
  int *goodflag = NULL;
  int *isoutlier = NULL;
  int *outlierindx = NULL;
  int **var_groupnum = NULL;
  int *var_Ngroup = NULL;
  int **var_sortindx = NULL;
  int *nbreak = NULL;
  double *t, *mag, *err;
  double magmedian, magstddev;

  double *grouperr = NULL;
  double *groupmag = NULL;
  double **groupvars = NULL;
  double **groupvars_minval = NULL;
  double **groupvars_maxval = NULL;
  int *Npoints_ingroup = NULL;
  int nsum = 0;
  int Nparam;
  int check;

  int Nparam_store = 0;
  double **decorr_terms = NULL;
  int *order_terms = NULL;
  double *Avector = NULL;
  double *A_errvector = NULL;

  int *origindx = NULL;
  double *model_mag = NULL;

  double outvarval;

  int *decorr_terms_startindx = NULL;
  int *decorr_terms_stopindx = NULL;

  _Splinedetrend *Splinedetrend;

  int isoutputmodelcoeffsopen = 0;
  FILE *outputmodelcoeffsfile = NULL;

  int Noutlier = 0;

  int Nextra = 0;

  char modelcoeffsoutname[MAXLEN];
  char modeloutname[MAXLEN];

  FILE *outmodel = NULL;


  /* Translate pointers and structures input by vartools into easier to use forms */
  Splinedetrend = (_Splinedetrend *)userdata;

  /* The number of points in the light curve, and vectors storing the
     time, magnitudes, and errors are accessed as illustrated
     below. Note that t, mag and err are vectors, t[0] is the time of
     the first observation in the light curve, t[NJD-1] is the time of
     the last observation, etc. */
  NJD=p->NJD[lc_num];
  t = p->t[lc_num];
  mag = p->mag[lc_num];
  err = p->sig[lc_num];
  
  if(NJD < 1) {
    /* Nothing to do, just exit */
    Splinedetrend->Noutlier[lc_num] = 0;
    Splinedetrend->Ngroups[lc_num] = 0;
    Splinedetrend->Nparamtotal[lc_num] = 0;
    Splinedetrend->magmedian[lc_num] = 0.;
    return;
  }

  if((model_mag = (double *) malloc(NJD * sizeof(double))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);

  if(Splinedetrend->Noutvar > 0) {
    if((decorr_terms_startindx = (int *) malloc(Splinedetrend->Noutvar*sizeof(int))) == NULL ||
       (decorr_terms_stopindx = (int *) malloc(Splinedetrend->Noutvar*sizeof(int))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    for(i=0; i < Splinedetrend->Noutvar; i++) {
      decorr_terms_startindx[i] = -1;
      decorr_terms_stopindx[i] = -1;
    }
  }


  if((goodflag = (int *) malloc(NJD*sizeof(int))) == NULL ||
     (groupid = (int *) malloc(NJD*sizeof(int))) == NULL ||
     (groupidoutlier = (int *) malloc(NJD*sizeof(int))) == NULL ||
     (origindx = (int *) malloc(NJD*sizeof(int))) == NULL || 
     (var_groupnum = (int **) malloc(Splinedetrend->Nvar*sizeof(int *))) == NULL ||
     (var_Ngroup = (int *) malloc(Splinedetrend->Nvar*sizeof(int))) == NULL ||
     (var_sortindx = (int **) malloc(Splinedetrend->Nvar*sizeof(int *))) == NULL ||
     (nbreak = (int *) malloc(Splinedetrend->Nvar*sizeof(int))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);

  for(k=0; k < Splinedetrend->Nvar; k++) {
    if((var_groupnum[k] = (int *) malloc(NJD*sizeof(int))) == NULL ||
       (var_sortindx[k] = (int *) malloc(NJD*sizeof(int))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
  }

  /* We will need the median later to add back to the subtracted signal */
  magmedian = VARTOOLS_median(NJD, mag);
  Splinedetrend->magmedian[lc_num] = magmedian;
  /* Figure out which points to exclude from the spline */
  if((isoutlier = (int *) malloc(NJD*sizeof(int))) == NULL ||
     (outlierindx = (int *) malloc(NJD * sizeof(int))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  Noutlier = 0;
  for(i=0; i < NJD; i++) isoutlier[i] = 0;
  if(Splinedetrend->sigmaclip[lc_num] > 0.) {
    /* Apply sigma-clipping to the light curves */
    magstddev = VARTOOLS_stddev(NJD, mag);
    for(i = 0; i < NJD; i++) {      
      if(!isnan(mag[i])) {
	if((fabs(mag[i] - magmedian) < Splinedetrend->sigmaclip[lc_num]*magstddev)) {
	  goodflag[i] = 1;
	  for(k=0; k < Splinedetrend->Nvar; k++) {
	    if(isnan(Splinedetrend->detrendvars[k].inputdata[lc_num][i])) {
	      goodflag[i] = 0;
	      break;
	    }
	  }
	} else {
	  goodflag[i] = 0;
	  isoutlier[i] = 1;
	  for(k=0; k < Splinedetrend->Nvar; k++) {
	    if(isnan(Splinedetrend->detrendvars[k].inputdata[lc_num][i])) {
	      isoutlier[i] = 0;
	      break;
	    }
	  }
	  if(isoutlier[i]) {
	    outlierindx[Noutlier] = i;
	    Noutlier++;
	  }
	}
      }
      else
	goodflag[i] = 0;
    }
  } else {
    for(i = 0; i < NJD; i++) {
      goodflag[i] = 1;
      for(k=0; k < Splinedetrend->Nvar; k++) {
	if(isnan(Splinedetrend->detrendvars[k].inputdata[lc_num][i])) {
	  goodflag[i] = 0;
	  break;
	}
      }
    }
  }
  
  Splinedetrend->Noutlier[lc_num] = Noutlier;

  /* Break the data into groups */
  Ngroup = 1;
  for(k = 0; k < Splinedetrend->Nvar; k++) {
    if(Splinedetrend->detrendvars[k].groupbygap) {
      Ngood = 0;
      for(i=0; i < NJD; i++) {
	if(goodflag[i]) {
	  var_sortindx[k][Ngood] = i;
	  Ngood++;
	}
      }
      VARTOOLS_sort_generic(Ngood, 0, var_sortindx[k], 1, VARTOOLS_TYPE_DOUBLE,
			    (void *) (Splinedetrend->detrendvars[k].inputdata[lc_num]),
			    1);
      i=0;
      while(i < Ngood ? goodflag[var_sortindx[k][i]] == 0 : 0) {
	var_groupnum[k][var_sortindx[k][i]] = -1;
	i++;
      }
      if(i < Ngood) {
	var_groupnum[k][var_sortindx[k][i]] = 1;
	var_Ngroup[k] = 1;
	j = i;
	i++;
	while(i < Ngood) {
	  if(goodflag[var_sortindx[k][i]]) {
	    if(Splinedetrend->detrendvars[k].inputdata[lc_num][var_sortindx[k][i]] - Splinedetrend->detrendvars[k].inputdata[lc_num][var_sortindx[k][j]] > Splinedetrend->detrendvars[k].gapsize) {
	      var_Ngroup[k]++;
	      var_groupnum[k][var_sortindx[k][i]] = var_Ngroup[k];
	      j = i;
	    } else {
	      var_groupnum[k][var_sortindx[k][i]] = var_Ngroup[k];
	      j = i;
	    }
	  }
	  i++;
	}
      } else {
	/* There are no good points - set any output model variables to 0, and
	   clean up the memory allocated up to this point */
	if(Splinedetrend->Noutvar > 0) {
	  for(i = 0; i < NJD; i++) {
	    for(k = 0; k < Splinedetrend->Noutvar; k++) {
	      Splinedetrend->outvars[k].outputdata[lc_num][i] = 0.;
	    }
	  }
	}
	if(model_mag != NULL) free(model_mag);
	
	if(decorr_terms_startindx != NULL) free(decorr_terms_startindx);
	
	if(decorr_terms_stopindx != NULL) free(decorr_terms_stopindx);
	
	if(goodflag != NULL) free(goodflag);
	
	if(groupid != NULL) free(groupid);
	
	if(origindx != NULL) free(origindx);
	
	if(var_groupnum != NULL) {
	  for(k=0; k < Splinedetrend->Nvar; k++) free(var_groupnum[k]);
	  free(var_groupnum);
	}
	
	if(var_Ngroup != NULL) free(var_Ngroup);
	
	if(var_sortindx != NULL) {
	  for(k=0; k < Splinedetrend->Nvar; k++) free(var_sortindx[k]);
	  free(var_sortindx);
	}
	
	if(nbreak != NULL) free(nbreak);
	
	if(groupmag != NULL) free(groupmag);
	
	if(grouperr != NULL) free(grouperr);
	
	if(groupvars != NULL) {
	  for(k = 0; k < Splinedetrend->Nvar; k++) free(groupvars[k]);
	  free(groupvars);
	}
	
	if(groupvars_minval != NULL) {
	  for(k=0; k < Splinedetrend->Nvar; k++) free(groupvars_minval[k]);
	  free(groupvars_minval);
	}
	
	if(groupvars_maxval != NULL) {
	  for(k=0; k < Splinedetrend->Nvar; k++) free(groupvars_maxval[k]);
	  free(groupvars_maxval);
	}
	
	if(Npoints_ingroup != NULL) free(Npoints_ingroup);
	
	if(decorr_terms != NULL) {
	  for(k=0; k < NJD; k++) free(decorr_terms[k]);
	  free(decorr_terms);
	}
	
	if(order_terms != NULL) free(order_terms);
	
	if(Avector != NULL) free(Avector);
	
	if(A_errvector != NULL) free(A_errvector);
	
	if(groupidoutlier != NULL) free(groupidoutlier);
	if(isoutlier != NULL) free(isoutlier);
	if(outlierindx != NULL) free(outlierindx);

	Splinedetrend->Ngroups[lc_num] = 0;
	Splinedetrend->Nparamtotal[lc_num] = 0;
	Splinedetrend->magmedian[lc_num] = 0.;
	return;
      }
    } else {
      var_Ngroup[k] = 1;
      for(i=0; i < NJD; i++) {
	if(goodflag[i]) {
	  var_groupnum[k][i] = 1;
	} else {
	  var_groupnum[k][i] = -1;
	}
      }
    }
  }
  
  Ngroup = 1;
  for(k = 0; k < Splinedetrend->Nvar; k++) {
    Ngroup = Ngroup*var_Ngroup[k];
  }


  if((groupmag = (double *) malloc(NJD * sizeof(double))) == NULL ||
     (grouperr = (double *) malloc(NJD * sizeof(double))) == NULL ||
     (groupvars = (double **) malloc(Splinedetrend->Nvar * sizeof(double *))) == NULL ||
     (groupvars_minval = (double **) malloc(Splinedetrend->Nvar * sizeof(double *))) == NULL ||
     (groupvars_maxval = (double **) malloc(Splinedetrend->Nvar * sizeof(double *))) == NULL ||
     (Npoints_ingroup = (int *) malloc(Ngroup * sizeof(int))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  
  for(g=0; g < Ngroup; g++) Npoints_ingroup[g] = 0;

  for(k = 0; k < Splinedetrend->Nvar; k++) {
    if((groupvars_minval[k] = (double *) malloc(Ngroup * sizeof(double))) == NULL ||
       (groupvars_maxval[k] = (double *) malloc(Ngroup * sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
  }


  for(i = 0; i < NJD; i++) {
    if(goodflag[i]) {
      groupid[i] = 0;
      nsum = 1;
      for(k=0; k < Splinedetrend->Nvar; k++) {
	if(var_groupnum[k][i] > 0) {
	  groupid[i] = groupid[i] + (var_groupnum[k][i] - 1)*nsum;
	} else {
	  groupid[i] = -1;
	  break;
	}
	nsum *= var_Ngroup[k];
      }
    } else {
      groupid[i] = -1;
    }
    if(groupid[i] >= 0) {
      if(!Npoints_ingroup[groupid[i]]) {
	for(k=0; k < Splinedetrend->Nvar; k++) {
	  groupvars_minval[k][groupid[i]] = Splinedetrend->detrendvars[k].inputdata[lc_num][i];
	  groupvars_maxval[k][groupid[i]] = Splinedetrend->detrendvars[k].inputdata[lc_num][i];
	}
      } else {
	for(k=0; k < Splinedetrend->Nvar; k++) {
	  if(Splinedetrend->detrendvars[k].inputdata[lc_num][i] < groupvars_minval[k][groupid[i]]) {
	    groupvars_minval[k][groupid[i]] = Splinedetrend->detrendvars[k].inputdata[lc_num][i];
	  }
	  if(Splinedetrend->detrendvars[k].inputdata[lc_num][i] > groupvars_maxval[k][groupid[i]]) {
	    groupvars_maxval[k][groupid[i]] = Splinedetrend->detrendvars[k].inputdata[lc_num][i];
	  }
	}
      }
      Npoints_ingroup[groupid[i]] += 1;
    }
  }

  /* Figure out which groups any outliers are in */
  if(Noutlier > 0) {
    for(j=0; j < Noutlier; j++) {
      i = outlierindx[j];
      groupidoutlier[i] = -1;
      for(g = 0; g < Ngroup; g++) {
	if(Npoints_ingroup[g] > 1) {
	  check = 1;
	  for(k=0; k < Splinedetrend->Nvar; k++) {
	    if(Splinedetrend->detrendvars[k].inputdata[lc_num][i] < groupvars_minval[k][g] ||
	       Splinedetrend->detrendvars[k].inputdata[lc_num][i] > groupvars_maxval[k][g]) {
	      check = 0;
	      break;
	    }
	  }
	  if(check) {
	    groupidoutlier[i] = g;
	    break;
	  }
	}
      }
    }
  }


  for(k=0; k < Splinedetrend->Nvar; k++) {
    if((groupvars[k] = (double *) malloc(NJD * sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
  }

  Splinedetrend->Ngroups[lc_num] = 0;

  Splinedetrend->Nparamtotal[lc_num] = 0;
  
  /* Now do the detrending group-by-group */
  for(g = 0; g < Ngroup; g++) {
    if(Npoints_ingroup[g] > 0) {
      Splinedetrend->Ngroups[lc_num] += 1;
      j = 0;
      Nparam = 0;
      Nextra = 0;
      for(i=0; i < NJD; i++) {
	if(groupid[i] == g) {
	  /* The point is in the group, and will be used in the fit */
	  origindx[j] = i;
	  groupmag[j] = mag[i];
	  grouperr[j] = err[i];
	  for(k=0; k < Splinedetrend->Nvar; k++) {
	    groupvars[k][j] = Splinedetrend->detrendvars[k].inputdata[lc_num][i];
	  }
	  j++;
	}
	else if(isoutlier[i] ? groupidoutlier[i] == g : 0) {
	  /* The point is in the group, but is an outlier that should not be
             used in the fit, but should have the model evaluated at the point
	  */
	  origindx[Npoints_ingroup[g] + Nextra] = i;
	  groupmag[Npoints_ingroup[g] + Nextra] = mag[i];
	  grouperr[Npoints_ingroup[g] + Nextra] = err[i];
	  for(k=0; k < Splinedetrend->Nvar; k++) {
	    groupvars[k][Npoints_ingroup[g] + Nextra] = Splinedetrend->detrendvars[k].inputdata[lc_num][i];
	  }
	  Nextra++;
	}
      }
      /* Build the decorrelation vectors */
      /* First determine the number of parameters */
      Nparam = 1;
      for(k = 0; k < Splinedetrend->Nvar; k++) {
	if(Splinedetrend->detrendvars[k].fittype == SPLINEDETREND_VARIABLE_FITTYPE_SPLINE) {
	  /* Determine the number of breaks needed */
	  nbreak[k] = floor((groupvars_maxval[k][g] - groupvars_minval[k][g])/Splinedetrend->detrendvars[k].knotspacing)+1;
	  if(nbreak[k] > 1)
	    Nparam += nbreak[k] + (Splinedetrend->detrendvars[k].order+1) - 2;
	}
	else if(Splinedetrend->detrendvars[k].fittype == SPLINEDETREND_VARIABLE_FITTYPE_POLY) {
	  Nparam += Splinedetrend->detrendvars[k].order;
	}
	else if(Splinedetrend->detrendvars[k].fittype == SPLINEDETREND_VARIABLE_FITTYPE_HARM) {
	  Nparam += 2*(Splinedetrend->detrendvars[k].nharm + 1);
	}
      }
      /* Check if we don't have enough data to decorrelate against */
      if(Npoints_ingroup[g] < Nparam) {
	Splinedetrend->Nparamtotal[lc_num] += 1;
	/*** Set the model equal to the weighted mean magnitude in this case **/
	outvarval = VARTOOLS_getweightedmean(Npoints_ingroup[g], groupmag, grouperr);
	
	if(Splinedetrend->outputmodelcoeffs) {
	  if(!isoutputmodelcoeffsopen) {
	    if(!Splinedetrend->outputmodelcoeffs_useformat) {
	      VARTOOLS_GetOutputFilename(modelcoeffsoutname, 
					 p->lcnames[lc_name_num],
					 Splinedetrend->modelcoeffsoutdir, 
					 "splinedetrend_modelcoeffs",
					 NULL, lc_name_num);
	    } else {
	      VARTOOLS_GetOutputFilename(modelcoeffsoutname, 
					 p->lcnames[lc_name_num],
					 Splinedetrend->modelcoeffsoutdir, 
					 "splinedetrend_modelcoeffs",
					 Splinedetrend->outputmodelcoeffs_format, 
					 lc_name_num);
	    }
	    if((outputmodelcoeffsfile = fopen(modelcoeffsoutname,"w")) == NULL) {
	      VARTOOLS_error2(ERR_CANNOTWRITE,modelcoeffsoutname);
	    }
	    isoutputmodelcoeffsopen = 1;
	  }
	  fprintf(outputmodelcoeffsfile,"#Model Coefficients for group %d\n",g);
	  fprintf(outputmodelcoeffsfile,"#Domain of group:\n");
	  for(k = 0; k < Splinedetrend->Nvar; k++) {
	    fprintf(outputmodelcoeffsfile,
		    "# %.17g <= %s <= %.17g\n", groupvars_minval[k][g],
		    Splinedetrend->detrendvars[k].varname,
		    groupvars_maxval[k][g]);
	  }
	  fprintf(outputmodelcoeffsfile,"#Coefficient values:\n");
	  fprintf(outputmodelcoeffsfile,"0.0 # constant term\n");
	  for(k = 0; k < Splinedetrend->Nvar; k++) {
	    if(Splinedetrend->detrendvars[k].fittype == SPLINEDETREND_VARIABLE_FITTYPE_SPLINE) {
	      if(nbreak[k] > 1) {
		for(ll = 0; ll < nbreak[k] + (Splinedetrend->detrendvars[k].order+1) - 2; ll++) {
		  fprintf(outputmodelcoeffsfile,
			  "0.0 # Basis spline %d for variable %s\n", 
			  ll, Splinedetrend->detrendvars[k].varname);
		}
	      }
	    }
	    else if(Splinedetrend->detrendvars[k].fittype == SPLINEDETREND_VARIABLE_FITTYPE_POLY) {
	      for(ll = 0; ll < Splinedetrend->detrendvars[k].order; ll++) {
		fprintf(outputmodelcoeffsfile,
			"0.0 # polynomial term (((%s)-%.17g)/%.17g)^%d\n", 
			Splinedetrend->detrendvars[k].varname,
			groupvars_minval[k][g], 
			(groupvars_maxval[k][g] - groupvars_minval[k][g]), 
			ll+1);
	      }
	    }
	    else if(Splinedetrend->detrendvars[k].fittype == SPLINEDETREND_VARIABLE_FITTYPE_HARM) {
	      for(ll = 0; ll < Splinedetrend->detrendvars[k].nharm; ll++) {
		fprintf(outputmodelcoeffsfile, 
			"0.0 # harmonic term sin(PI*%d*((%s)-%.17g)/%.17g)\n", 
			ll+1, Splinedetrend->detrendvars[k].varname, 
			groupvars_minval[k][g], 
			(groupvars_maxval[k][g] - groupvars_minval[k][g]));
		fprintf(outputmodelcoeffsfile, 
			"0.0 # harmonic term cos(PI*%d*((%s)-%.17g)/%.17g)\n", 
			ll+1, Splinedetrend->detrendvars[k].varname, 
			groupvars_minval[k][g], 
			(groupvars_maxval[k][g] - groupvars_minval[k][g]));
	      }
	    }
	  }
	  fprintf(outputmodelcoeffsfile,"\n");
	}

	/* Copy back the decorrelated magnitude values into the original magnitude
	   vector, store the model if needed, and evaluate any variables that
	   need to be stored */
	for(j = 0; j < Npoints_ingroup[g] + Nextra; j++) {
	  model_mag[origindx[j]] = outvarval;
	  mag[origindx[j]] = mag[origindx[j]] - outvarval + magmedian;
	  if(Splinedetrend->Noutvar > 0) {
	    for(k = 0; k < Splinedetrend->Noutvar; k++) {
	      Splinedetrend->outvars[k].outputdata[lc_num][origindx[j]] = 0.;
	    }
	  }
	}
      } else {

	Splinedetrend->Nparamtotal[lc_num] += Nparam;

	/* Allocate memory for performing the fit */
	if(Nparam > Nparam_store) {
	  if(!Nparam_store) {
	    if((decorr_terms = (double **) malloc(NJD*sizeof(double *))) == NULL ||
	       (order_terms = (int *) malloc(Nparam * sizeof(int))) == NULL ||
	       (Avector = (double *) malloc(Nparam * sizeof(double))) == NULL ||
	       (A_errvector = (double *) malloc(Nparam * sizeof(double))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	    for(i = 0; i < NJD; i++) {
	      if((decorr_terms[i] = (double *) malloc(Nparam*sizeof(double))) == NULL) {
		VARTOOLS_error(ERR_MEMALLOC);
	      }
	    }
	  } else {
	    if((order_terms = (int *) realloc(order_terms, Nparam * sizeof(int))) == NULL ||
	       (Avector = (double *) realloc(Avector, Nparam * sizeof(double))) == NULL ||
	       (A_errvector = (double *) realloc(A_errvector, Nparam * sizeof(double))) == NULL) {
	      VARTOOLS_error(ERR_MEMALLOC);
	    }
	    for(i = 0; i < NJD; i++) {
	      if((decorr_terms[i] = (double *) realloc(decorr_terms[i], Nparam * sizeof(double))) == NULL) {
		VARTOOLS_error(ERR_MEMALLOC);
	      }
	    }
	  }
	  Nparam_store = Nparam;
	}
	/* Now evaluate the decorrelation terms to fit against */
	ll = 0;
	for(i = 0; i < Npoints_ingroup[g] + Nextra; i++) {
	  decorr_terms[i][ll] = 1.;
	}
	order_terms[ll] = 1;
	ll++;
	for(k = 0; k < Splinedetrend->Nvar; k++) {
	  if(Splinedetrend->Noutvar > 0) {
	    for(i = 0; i < Splinedetrend->Noutvar; i++) {
	      if(Splinedetrend->outvars[i].invar == &(Splinedetrend->detrendvars[k])) {
		decorr_terms_startindx[i] = ll;
		break;
	      }
	    }
	  }
	  if(Splinedetrend->detrendvars[k].fittype == SPLINEDETREND_VARIABLE_FITTYPE_SPLINE) {
	    if(nbreak[k] > 1) {
	      SplineDetrend_FillSplineTerms(&ll, decorr_terms, Npoints_ingroup[g]+Nextra, groupvars[k], groupvars_minval[k][g], groupvars_maxval[k][g], nbreak[k], Splinedetrend->detrendvars[k].order, order_terms);
	    }
	  }
	  else if(Splinedetrend->detrendvars[k].fittype == SPLINEDETREND_VARIABLE_FITTYPE_POLY) {
	    SplineDetrend_FillPolyTerms(&ll, decorr_terms, Npoints_ingroup[g]+Nextra, groupvars[k], groupvars_minval[k][g], groupvars_maxval[k][g], Splinedetrend->detrendvars[k].order, order_terms);
	  }
	  else if(Splinedetrend->detrendvars[k].fittype == SPLINEDETREND_VARIABLE_FITTYPE_HARM) {
	    SplineDetrend_FillHarmTerms(&ll, decorr_terms, Npoints_ingroup[g]+Nextra, groupvars[k], groupvars_minval[k][g], groupvars_maxval[k][g], Splinedetrend->detrendvars[k].nharm, order_terms);
	  }
	  if(Splinedetrend->Noutvar > 0) {
	    for(i = 0; i < Splinedetrend->Noutvar; i++) {
	      if(Splinedetrend->outvars[i].invar == &(Splinedetrend->detrendvars[k])) {
		decorr_terms_stopindx[i] = ll;
		break;
	      }
	    }
	  }
	}
	/* Perform the fit */
	VARTOOLS_docorr(groupmag, grouperr, Npoints_ingroup[g], Nparam, decorr_terms, order_terms, Avector, A_errvector, 0., 0, 0, NULL, lc_name_num, lc_num);
      
	if(Splinedetrend->outputmodelcoeffs) {
	  if(!isoutputmodelcoeffsopen) {
	    if(!Splinedetrend->outputmodelcoeffs_useformat) {
	      VARTOOLS_GetOutputFilename(modelcoeffsoutname, 
					 p->lcnames[lc_name_num],
					 Splinedetrend->modelcoeffsoutdir, 
					 "splinedetrend_modelcoeffs",
					 NULL, lc_name_num);
	    } else {
	      VARTOOLS_GetOutputFilename(modelcoeffsoutname, 
					 p->lcnames[lc_name_num],
					 Splinedetrend->modelcoeffsoutdir, 
					 "splinedetrend_modelcoeffs",
					 Splinedetrend->outputmodelcoeffs_format, 
					 lc_name_num);
	    }
	    if((outputmodelcoeffsfile = fopen(modelcoeffsoutname,"w")) == NULL) {
	      VARTOOLS_error2(ERR_CANNOTWRITE,modelcoeffsoutname);
	    }
	    isoutputmodelcoeffsopen = 1;
	  }
	  fprintf(outputmodelcoeffsfile,"#Model Coefficients for group %d\n",g);
	  fprintf(outputmodelcoeffsfile,"#Domain of group:\n");
	  for(k = 0; k < Splinedetrend->Nvar; k++) {
	    fprintf(outputmodelcoeffsfile,
		    "# %.17g <= %s <= %.17g\n", groupvars_minval[k][g],
		    Splinedetrend->detrendvars[k].varname,
		    groupvars_maxval[k][g]);
	  }
	  fprintf(outputmodelcoeffsfile,"#Coefficient values:\n");
	  fprintf(outputmodelcoeffsfile,"%.17g # constant term\n", Avector[0]);
	  ii = 1;
	  for(k = 0; k < Splinedetrend->Nvar; k++) {
	    if(Splinedetrend->detrendvars[k].fittype == SPLINEDETREND_VARIABLE_FITTYPE_SPLINE) {
	      if(nbreak[k] > 1) {
		for(ll = 0; ll < nbreak[k] + (Splinedetrend->detrendvars[k].order+1) - 2; ll++) {
		  fprintf(outputmodelcoeffsfile,
			  "%.17g # Basis spline %d for variable %s\n", 
			  Avector[ii], ll, Splinedetrend->detrendvars[k].varname);
		  ii++;
		}
	      }
	    }
	    else if(Splinedetrend->detrendvars[k].fittype == SPLINEDETREND_VARIABLE_FITTYPE_POLY) {
	      for(ll = 0; ll < Splinedetrend->detrendvars[k].order; ll++) {
		fprintf(outputmodelcoeffsfile,
			"%.17g # polynomial term (((%s)-%.17g)/%.17g)^%d\n", 
			Avector[ii],
			Splinedetrend->detrendvars[k].varname,
			groupvars_minval[k][g], 
			(groupvars_maxval[k][g] - groupvars_minval[k][g]), 
			ll+1);
		ii++;
	      }
	    }
	    else if(Splinedetrend->detrendvars[k].fittype == SPLINEDETREND_VARIABLE_FITTYPE_HARM) {
	      for(ll = 0; ll < Splinedetrend->detrendvars[k].nharm; ll++) {
		fprintf(outputmodelcoeffsfile, 
			"%.17g # harmonic term sin(PI*%d*((%s)-%.17g)/%.17g)\n", 
			Avector[ii],
			ll+1, Splinedetrend->detrendvars[k].varname, 
			groupvars_minval[k][g], 
			(groupvars_maxval[k][g] - groupvars_minval[k][g]));
		ii++;
		fprintf(outputmodelcoeffsfile, 
			"%.17g # harmonic term cos(PI*%d*((%s)-%.17g)/%.17g)\n", 
			Avector[ii],
			ll+1, Splinedetrend->detrendvars[k].varname, 
			groupvars_minval[k][g], 
			(groupvars_maxval[k][g] - groupvars_minval[k][g]));
		ii++;
	      }
	    }
	  }
	  fprintf(outputmodelcoeffsfile,"\n");
	}
	
	/* Calculate the model magnitude values, and subtract these
           from the input magnitude values, store the model if needed,
           and evaluate any variables that need to be stored */
	for(j = 0; j < Npoints_ingroup[g]+Nextra; j++) {
	  model_mag[origindx[j]] = 0.;
	  for(k = 0; k < Nparam; k++) {
	    model_mag[origindx[j]] += decorr_terms[j][k]*Avector[k];
	  }
	  mag[origindx[j]] = mag[origindx[j]] - model_mag[origindx[j]] + magmedian;
	  if(Splinedetrend->Noutvar > 0) {
	    for(k = 0; k < Splinedetrend->Noutvar; k++) {
	      outvarval = 0.;
	      for(ii = decorr_terms_startindx[k]; ii < decorr_terms_stopindx[k]; ii++) {
		outvarval += decorr_terms[j][ii]*Avector[ii];
	      }
	      Splinedetrend->outvars[k].outputdata[lc_num][origindx[j]] = outvarval;
	    }
	  }
	}
      }
    }
  }

  /**** Close the output coefficients file if it is open *****/
  if(isoutputmodelcoeffsopen)
    fclose(outputmodelcoeffsfile);

  if(Splinedetrend->outputmodel) {
    if(!Splinedetrend->outputmodel_useformat) {
      VARTOOLS_GetOutputFilename(modeloutname, p->lcnames[lc_name_num],
				 Splinedetrend->modeloutdir, "splinedetrend_model",
				 NULL, lc_name_num);
    } else {
      VARTOOLS_GetOutputFilename(modeloutname, p->lcnames[lc_name_num],
				 Splinedetrend->modeloutdir, "splinedetrend_model",
				 Splinedetrend->outputmodel_format, lc_name_num);
    }
    if((outmodel = fopen(modeloutname, "w")) == NULL)
      VARTOOLS_error2(ERR_CANNOTWRITE,modeloutname);
    
    fprintf(outmodel,"#Time Mag_obs Mag_model Error\n");
    for(i=0; i < NJD; i++) {
      if(groupid[i] > -1) {
	fprintf(outmodel,"%.17g %.17g %.17g %.17g\n",
		t[i], mag[i]+model_mag[i]-magmedian, model_mag[i], err[i]);
      }
    }
    fclose(outmodel);
  }

  if(Splinedetrend->Noutvar > 0) {
    for(i=0; i < NJD; i++) {
      if(groupid[i] < 0) {
	if(isoutlier[i] ? groupidoutlier[i] < 0 : 1) {
	  for(k = 0; k < Splinedetrend->Noutvar; k++) {
	    outvarval = 0.;
	    Splinedetrend->outvars[k].outputdata[lc_num][i] = outvarval;
	  }
	}
      }
    }
  }
  
  /* Clean up the allocated memory and return */
  if(model_mag != NULL) free(model_mag);
  
  if(decorr_terms_startindx != NULL) free(decorr_terms_startindx);

  if(decorr_terms_stopindx != NULL) free(decorr_terms_stopindx);
  
  if(goodflag != NULL) free(goodflag);

  if(groupid != NULL) free(groupid);

  if(origindx != NULL) free(origindx);

  if(var_groupnum != NULL) {
    for(k=0; k < Splinedetrend->Nvar; k++) free(var_groupnum[k]);
    free(var_groupnum);
  }

  if(var_Ngroup != NULL) free(var_Ngroup);

  if(var_sortindx != NULL) {
    for(k=0; k < Splinedetrend->Nvar; k++) free(var_sortindx[k]);
    free(var_sortindx);
  }

  if(nbreak != NULL) free(nbreak);

  if(groupmag != NULL) free(groupmag);
  
  if(grouperr != NULL) free(grouperr);
  
  if(groupvars != NULL) {
    for(k = 0; k < Splinedetrend->Nvar; k++) free(groupvars[k]);
    free(groupvars);
  }

  if(groupvars_minval != NULL) {
    for(k=0; k < Splinedetrend->Nvar; k++) free(groupvars_minval[k]);
    free(groupvars_minval);
  }

  if(groupvars_maxval != NULL) {
    for(k=0; k < Splinedetrend->Nvar; k++) free(groupvars_maxval[k]);
    free(groupvars_maxval);
  }

  if(Npoints_ingroup != NULL) free(Npoints_ingroup);

  if(decorr_terms != NULL) {
    for(k=0; k < NJD; k++) free(decorr_terms[k]);
    free(decorr_terms);
  }

  if(order_terms != NULL) free(order_terms);

  if(Avector != NULL) free(Avector);

  if(A_errvector != NULL) free(A_errvector);

  if(groupidoutlier != NULL) free(groupidoutlier);

  if(isoutlier != NULL) free(isoutlier);
  
  if(outlierindx != NULL) free(outlierindx);
      
}

/* This function parses a string of the form:
   detrendvec1:<"spline":knotspacing:order|"poly":order|"harm":nharm>[:"groupbygap":gapsize][,detrendvec2:<"spline":knotspacing:order|"poly":order|"harm":nharm>[:"groupbygap":gapsize],...] */
int ParseSplineDetrendVarlist(char *varlist,  _Splinedetrend *Splinedetrend, ProgramData *p, Command *c)
{
    int i, j, k, termscan, startindex;
    char *parsecopy;
    double knotspacing;
    int order;
    int nharm;
    char *groupbygap;
    double gapsize;
    int sizestring;
    _SplinedetrendVariable *v;

    sizestring = strlen(varlist) + 1;

    if((parsecopy = (char *) malloc(sizestring*sizeof(char))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);

    Splinedetrend->Nvar = 0;
    termscan = 0;
    i = 0, j = 0;
    startindex = 0;
    do {
      if(varlist[i] == ',' || varlist[i] == ':' || varlist[i] == '\0') {
	if(startindex == i) {
	  free(parsecopy);
	  return 1;
	}
	parsecopy[j] = '\0';
	if(termscan == 0) {
	  Splinedetrend->Nvar += 1;
	  if(Splinedetrend->Nvar == 1) {
	    if((Splinedetrend->detrendvars = (_SplinedetrendVariable *) malloc(sizeof(_SplinedetrendVariable))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	  } else {
	    if((Splinedetrend->detrendvars = (_SplinedetrendVariable *) realloc(Splinedetrend->detrendvars, Splinedetrend->Nvar * sizeof(_SplinedetrendVariable))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	  }
	  v = &(Splinedetrend->detrendvars[Splinedetrend->Nvar-1]);
	  if((v->varname = (char *) malloc(strlen(parsecopy)+1)) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	  sprintf(v->varname, "%s", parsecopy);
	  v->groupbygap = 0;
	}
	else if(termscan == 1) {
	  if(!strcmp(parsecopy, "spline")) {
	    v->fittype = SPLINEDETREND_VARIABLE_FITTYPE_SPLINE;
	  } else if(!strcmp(parsecopy, "poly")) {
	    v->fittype = SPLINEDETREND_VARIABLE_FITTYPE_POLY;
	  } else if(!strcmp(parsecopy, "harm")) {
	    v->fittype = SPLINEDETREND_VARIABLE_FITTYPE_HARM;
	  }
	  else {
	    free(parsecopy);
	    return 1;
	  }
	}
	else if(termscan == 2) {
	  if(v->fittype == SPLINEDETREND_VARIABLE_FITTYPE_SPLINE) {
	    v->knotspacing = atof(parsecopy);
	  } else if(v->fittype == SPLINEDETREND_VARIABLE_FITTYPE_POLY) {
	    v->order = atoi(parsecopy);
	  } else if(v->fittype == SPLINEDETREND_VARIABLE_FITTYPE_HARM) {
	    v->nharm = atoi(parsecopy);
	  }
	}
	else if(termscan == 3) {
	  if(v->fittype == SPLINEDETREND_VARIABLE_FITTYPE_SPLINE) {
	    v->order = atoi(parsecopy);
	  } else {
	    if(!strcmp(parsecopy, "groupbygap")) {
	      v->groupbygap = 1;
	    } else {
	      free(parsecopy);
	      return 1;
	    }
	  }
	}
	else if(termscan == 4) {
	  if(v->fittype == SPLINEDETREND_VARIABLE_FITTYPE_SPLINE) {
	    if(!strcmp(parsecopy, "groupbygap")) {
	      v->groupbygap = 1;
	    } else {
	      free(parsecopy);
	      return 1;
	    }
	  } else {
	    if(v->groupbygap == 1) {
	      v->gapsize = atof(parsecopy);
	    } else {
	      free(parsecopy);
	      return 1;
	    }
	  }
	}
	else if(termscan == 5) {
	  if(v->fittype == SPLINEDETREND_VARIABLE_FITTYPE_SPLINE) {
	    if(v->groupbygap == 1) {
	      v->gapsize = atof(parsecopy);
	    } else {
	      free(parsecopy);
	      return 1;
	    }
	  } else {
	    free(parsecopy);
	    return 1;
	  }
	} else {
	  free(parsecopy);
	  return 1;
	}
	j = 0;
	startindex = i+1;
	if(varlist[i] == ':') {
	  termscan++;
	} else {
	  if((termscan < 2) || (v->fittype == SPLINEDETREND_VARIABLE_FITTYPE_SPLINE && termscan < 3) || (v->groupbygap == 1 && ((v->fittype == SPLINEDETREND_VARIABLE_FITTYPE_SPLINE && termscan < 5) || (v->fittype != SPLINEDETREND_VARIABLE_FITTYPE_SPLINE && termscan < 4)))) {
	    free(parsecopy);
	    return 1;
	  }
	  termscan = 0;
	  if(varlist[i] == '\0')
	    break;
	}
      } else {
	parsecopy[j] = varlist[i];
	j++;
      }
      i++;
    } while(1);

    /* Register the data vectors to store the input variable data for performing
       the detrending */
    for(k = 0; k < Splinedetrend->Nvar; k++) {
      VARTOOLS_RegisterDataVector(p, c, (void *) (&(Splinedetrend->detrendvars[k].inputdata)), VARTOOLS_TYPE_DOUBLE, 0, VARTOOLS_SOURCE_EVALEXPRESSION_LC, 0, NULL, Splinedetrend->detrendvars[k].varname, NULL);
    }
    
    free(parsecopy);
    return(0);
}

/* This function parses a string of the form:
	outvarname1[:inputvarsignal1[,outvarname2:inputvarsignal2...]]
 */
int ParseSplineDetrendOutModelVar(char *varlist,  _Splinedetrend *Splinedetrend, ProgramData *p, Command *c)
{
    int i, j, k, termscan, startindex;
    int ll;
    char *parsecopy;
    int sizestring;
    _SplinedetrendOutputVariable *v;

    sizestring = strlen(varlist) + 1;

    if((parsecopy = (char *) malloc(sizestring*sizeof(char))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);

    Splinedetrend->Noutvar = 0;
    termscan = 0;
    i = 0, j = 0;
    startindex = 0;
    do {
      if(varlist[i] == ',' || varlist[i] == ':' || varlist[i] == '\0') {
	if(startindex == i) {
	  free(parsecopy);
	  return 1;
	}
	parsecopy[j] = '\0';
	if(termscan == 0) {
	  Splinedetrend->Noutvar += 1;
	  if(Splinedetrend->Noutvar == 1) {
	    if((Splinedetrend->outvars = (_SplinedetrendOutputVariable *) malloc(sizeof(_SplinedetrendOutputVariable))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	  } else {
	    if((Splinedetrend->outvars = (_SplinedetrendOutputVariable *) realloc(Splinedetrend->outvars, Splinedetrend->Noutvar * sizeof(_SplinedetrendOutputVariable))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	  }
	  v = &(Splinedetrend->outvars[Splinedetrend->Noutvar-1]);
	  if((v->varname = (char *) malloc(strlen(parsecopy)+1)) == NULL)
	    VARTOOLS_error(ERR_MEMALLOC);
	  sprintf(v->varname, "%s", parsecopy);
	}
	else if(termscan == 1) {
	  for(ll=0; ll < Splinedetrend->Nvar; ll++) {
	    if(!strcmp(parsecopy,Splinedetrend->detrendvars[ll].varname)) {
	      v->invar = &(Splinedetrend->detrendvars[ll]);
	      break;
	    }
	  }
	  if(ll == Splinedetrend->Nvar) {
	    fprintf(stderr,"Error parsing the omodelvariable option to the -splinedetrend command. The variable %s is not one of the detrendvec variables used.\n", parsecopy);
	    exit(1);
	  }
	} else {
	  free(parsecopy);
	  return 1;
	}
	j = 0;
	startindex = i+1;
	if(varlist[i] == ':') {
	  termscan++;
	} else {
	  if(termscan < 1) {
	    free(parsecopy);
	    return 1;
	  }
	  termscan = 0;
	  if(varlist[i] == '\0')
	    break;
	}
      } else {
	parsecopy[j] = varlist[i];
	j++;
      }
      i++;
    } while(1);

    /* Register the data vectors to store the output variable data */
    for(k = 0; k < Splinedetrend->Noutvar; k++) {
      VARTOOLS_RegisterDataVector(p, c, (void *) (&Splinedetrend->outvars[k].outputdata),
				  VARTOOLS_TYPE_DOUBLE, 0, VARTOOLS_SOURCE_LC,
				  0, NULL, NULL, -1, NULL, 
				  Splinedetrend->outvars[k].varname);
    }
    
    free(parsecopy);
    return(0);
}


/* Fill in the decorrelation matrix with basis splines */
void SplineDetrend_FillSplineTerms(int *ll, double **decorr_terms, int Np, double *datavec, double minval, double maxval, int nbreak, int splineorder, int *order_terms) {
  gsl_bspline_workspace *bw;
  gsl_vector *B;
  int ncoeffs;
  int j, k, i;

  ncoeffs = nbreak + splineorder+1 - 2;

  /* allocate a Bspline workspace */
  bw = gsl_bspline_alloc(splineorder+1, nbreak);
  B = gsl_vector_alloc(ncoeffs);

  /* set up the Bspline knots */
  gsl_bspline_knots_uniform(minval, maxval, bw);
    
  /* Evaluate the basis splines and store the results in decorr_terms */
  for(j=0; j < Np; j++) {
    k = *ll;
    gsl_bspline_eval(datavec[j], B, bw);
    for(i = 0; i < ncoeffs; i++) {
      decorr_terms[j][k] = gsl_vector_get(B, i);
      k++;
    }
  }

  /* Set all of the order_terms components to 1; this is needed by the 
     VARTOOLS_docorr function */
  k = *ll;
  for(i = 0; i < ncoeffs; i++) {
    order_terms[k] = 1;
    k++;
  }

  k = *ll + ncoeffs;
  *ll = k;

  /* Clean up */
  gsl_bspline_free(bw);
  gsl_vector_free(B);

  return;
}

/* Fill in the decorrelation matrix with polynomial terms */
void SplineDetrend_FillPolyTerms(int *ll, double **decorr_terms, int Np, double *datavec, double minval, double maxval, int polyorder, int *order_terms) {
  double term;
  int j, k, i;

  /* Evaluate the polynomial and store the results in decorr_terms */
  for(j=0; j < Np; j++) {
    k = *ll;
    term = 1;
    for(i = 1; i <= polyorder; i++) {
      term = term*(datavec[j] - minval)/(maxval - minval);
      decorr_terms[j][k] = term;
      k++;
    }
  }

  /* Set all of the order_terms components to 1; this is needed by the 
     VARTOOLS_docorr function */
  k = *ll;
  for(i = 1; i <= polyorder; i++) {
    order_terms[k] = 1;
    k++;
  }

  k = *ll + polyorder;
  *ll = k;

  return;
}

/* Fill in the decorrelation matrix with harmonic terms */
void SplineDetrend_FillHarmTerms(int *ll, double **decorr_terms, int Np, double *datavec, double minval, double maxval, int nharm, int *order_terms) {
  double term;
  int j, k, i;

  /* Evaluate the harmonic series and store the results in decorr_terms */
  for(j=0; j < Np; j++) {
    k = *ll;
    for(i = 0; i <= nharm; i++) {
      term = sin(M_PI*(i+1)*(datavec[j] - minval)/(maxval - minval));
      decorr_terms[j][k] = term;
      k++;
      term = cos(M_PI*(i+1)*(datavec[j] - minval)/(maxval - minval));
      decorr_terms[j][k] = term;
      k++;
    }
  }

  /* Set all of the order_terms components to 1; this is needed by the 
     VARTOOLS_docorr function */
  k = *ll;
  for(i = 0; i <= nharm; i++) {
    order_terms[k] = 1;
    k++;
    order_terms[k] = 1;
    k++;
  }

  k = *ll + 2*(nharm+1);
  *ll = k;

  return;
}
