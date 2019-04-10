#include "../../src/vartools.h"
#include <stdio.h>
#include <stdlib.h>
#include "ftuneven.h"

/* This is the source code to implement the -ftuneven command which computes
   the complex Fourier Transform of data with arbitrary spacing.

   The code for the main function is based on matlab code provided by 
   Jeff Scargle, who also came up with the algorithm.
*/

void ftuneven_Initialize(char *commandname,
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
  sprintf(commandname,"-ftuneven");

  /* This is 1 if the command requires all light curves to be read
     at once, otherwise set it to 0. */
  *RequireReadAll = 0;

  /* This is 1 if the command requires input light curves to be sorted by
     time, otherwise set it to 0. */
  *RequireSortLC = 1;
  
  /* This is 1 if the command requires the times of observation in the
     input light curve to be unique */
  *RequireDistinctTimes = 0;

  /* You should define a structure to store the data needed for this
     command, below you would replace "_FTuneven"
     with the type name of your structure.

     See ftuneven.h
*/
  *sizeuserdata = sizeof(_FTuneven);
}

int ftuneven_ParseCL(ProgramData *p, Command *c,
		     void *userdata, int *iret, char **argv, int argc)
/* Every library with the name $LIBNAME.so must contain a function
   with the name $LIBNAME_ParseCL to check that the user has issued a
   valid syntax, and set-up the command-data accordingly.

   Note that any vectors/arrays used by the command should also be registered
   by this function so that memory will be allocated as needed.

Expected syntax for this example:
   -ftuneven <"outputvectors" frequency_vecname FT_real_vecname FT_imag_vecname Periodogram_vecname | "outputfile" <outdir ["format" fmt]> | "outputvectorsandfile" frequency_vecname FT_real_vecname FT_imag_vecname Periodogram_vecname <outdir ["format" fmt]>> <"freqauto" | "freqrange" minfreq maxfreq freqstep | "freqlist" <"variable" varname | "file" filename>> ["ft_sign" val] ["tt_zero" val]

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
  int j, k;

  char *line;
  size_t line_size = MAXLEN;
  
  FILE *freqfile;

  /* Cast the input userdata pointer to the data structure type used
     for this command */
  _FTuneven *FTuneven;
  FTuneven = (_FTuneven *) userdata;

  i = *iret;

  /* We will parse this command line manually */
  if(i >= argc) {
    /* There are no terms left on the command-line, return an error */
    return 1;
  }
  if(!strcmp(argv[i],"outputvectors")) {
    FTuneven->outputtype = FTUNEVEN_OUTPUTVECTORS;
    i++;
    if(i >= argc)
      return 1;
    FTuneven->frequency_vecname = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
    sprintf(FTuneven->frequency_vecname,"%s",argv[i]);
    VARTOOLS_RegisterDataVector(p, c, (void *) (&FTuneven->frequency_vec),
				VARTOOLS_TYPE_DOUBLE, 0, VARTOOLS_SOURCE_LC,
				0, NULL, FTuneven->frequency_vecname,
				-1, NULL, FTuneven->frequency_vecname);
    i++;
    if(i >= argc)
      return 1;
    FTuneven->FT_real_vecname = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
    sprintf(FTuneven->FT_real_vecname,"%s",argv[i]);
    VARTOOLS_RegisterDataVector(p, c, (void *) (&FTuneven->FT_real_vec),
				VARTOOLS_TYPE_DOUBLE, 0, VARTOOLS_SOURCE_LC,
				0, NULL, FTuneven->FT_real_vecname,
				-1, NULL, FTuneven->FT_real_vecname);
    i++;
    if(i >= argc)
      return 1;
    FTuneven->FT_imag_vecname = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
    sprintf(FTuneven->FT_imag_vecname,"%s",argv[i]);
    VARTOOLS_RegisterDataVector(p, c, (void *) (&FTuneven->FT_imag_vec),
				VARTOOLS_TYPE_DOUBLE, 0, VARTOOLS_SOURCE_LC,
				0, NULL, FTuneven->FT_imag_vecname,
				-1, NULL, FTuneven->FT_imag_vecname);
    i++;
    if(i >= argc)
      return 1;
    FTuneven->Periodogram_vecname = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
    sprintf(FTuneven->Periodogram_vecname,"%s",argv[i]);
    VARTOOLS_RegisterDataVector(p, c, (void *) (&FTuneven->Periodogram_vec),
				VARTOOLS_TYPE_DOUBLE, 0, VARTOOLS_SOURCE_LC,
				0, NULL, FTuneven->Periodogram_vecname,
				-1, NULL, FTuneven->Periodogram_vecname);
    i++;
  } 
  else if(!strcmp(argv[i],"outputfile")) {
    FTuneven->outputtype = FTUNEVEN_OUTPUTFILE;
    check = VARTOOLS_ParseOutNameKeyword(p, c, &i, argv, argc,
			       "outputfile", &j, &FTuneven->FT_outdir,
				&FTuneven->FT_nameformatflag,
				&FTuneven->FT_nameformat);
    if(FTuneven->FT_nameformatflag) i++;
    if(check) {
      *iret = i; return 1;
    }
  }
  else if(!strcmp(argv[i],"outputvectorsandfile")) {
    FTuneven->outputtype = FTUNEVEN_OUTPUTVECTORSANDFILE;
    check = VARTOOLS_ParseOutNameKeyword(p, c, &i, argv, argc,
			       "outputvectorsandfile", &j, &FTuneven->FT_outdir,
				&FTuneven->FT_nameformatflag,
				&FTuneven->FT_nameformat);
    if(FTuneven->FT_nameformatflag) i++;
    if(check) {
      *iret = i; return 1;
    }
    if(i >= argc)
      return 1;
    FTuneven->frequency_vecname = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
    sprintf(FTuneven->frequency_vecname,"%s",argv[i]);
    VARTOOLS_RegisterDataVector(p, c, (void *) (&FTuneven->frequency_vec),
				VARTOOLS_TYPE_DOUBLE, 0, VARTOOLS_SOURCE_LC,
				0, NULL, FTuneven->frequency_vecname,
				-1, NULL, FTuneven->frequency_vecname);
    i++;
    if(i >= argc)
      return 1;
    FTuneven->FT_real_vecname = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
    sprintf(FTuneven->FT_real_vecname,"%s",argv[i]);
    VARTOOLS_RegisterDataVector(p, c, (void *) (&FTuneven->FT_real_vec),
				VARTOOLS_TYPE_DOUBLE, 0, VARTOOLS_SOURCE_LC,
				0, NULL, FTuneven->FT_real_vecname,
				-1, NULL, FTuneven->FT_real_vecname);
    i++;
    if(i >= argc)
      return 1;
    FTuneven->FT_imag_vecname = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
    sprintf(FTuneven->FT_imag_vecname,"%s",argv[i]);
    VARTOOLS_RegisterDataVector(p, c, (void *) (&FTuneven->FT_imag_vec),
				VARTOOLS_TYPE_DOUBLE, 0, VARTOOLS_SOURCE_LC,
				0, NULL, FTuneven->FT_imag_vecname,
				-1, NULL, FTuneven->FT_imag_vecname);
    i++;
    if(i >= argc)
      return 1;
    FTuneven->Periodogram_vecname = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
    sprintf(FTuneven->Periodogram_vecname,"%s",argv[i]);
    VARTOOLS_RegisterDataVector(p, c, (void *) (&FTuneven->Periodogram_vec),
				VARTOOLS_TYPE_DOUBLE, 0, VARTOOLS_SOURCE_LC,
				0, NULL, FTuneven->Periodogram_vecname,
				-1, NULL, FTuneven->Periodogram_vecname);
    i++;
  }
  else {
    return 1;
  }

  if(i >= argc) {
    *iret = i-1; return 1;
  }
  if(!strcmp(argv[i],"freqauto")) {
    FTuneven->freqtype = FTUNEVEN_FREQTYPE_AUTO;
  }
  else if(!strcmp(argv[i],"freqrange")) {
    FTuneven->freqtype = FTUNEVEN_FREQTYPE_RANGE;
    i++;
    if(i >= argc) {
      *iret = i-1; return 1;
    }
    if(VARTOOLS_ParseParameter(p, c, &i, argv, argc, "minfreq", 1, VARTOOLS_TYPE_DOUBLE,
		       &(FTuneven->minfreq), 0, 0, "minfreq", 0, 0.)) {
      *iret = i; return 1;
    }
    if(VARTOOLS_ParseParameter(p, c, &i, argv, argc, "maxfreq", 1, VARTOOLS_TYPE_DOUBLE,
		       &(FTuneven->maxfreq), 0, 0, "maxfreq", 0, 0.)) {
      *iret = i; return 1;
    }
    if(VARTOOLS_ParseParameter(p, c, &i, argv, argc, "freqstep", 1, VARTOOLS_TYPE_DOUBLE,
		       &(FTuneven->freqstep), 0, 0, "freqstep", 0, 0.)) {
      *iret = i; return 1;
    }
    i--;
  } else if(!strcmp(argv[i],"freqvariable")) {
    FTuneven->freqtype = FTUNEVEN_FREQTYPE_VARIABLE;
    i++;
    if(i >= argc) {
      *iret = i-1; return 1;
    }
    FTuneven->freqlist_varname = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
    sprintf(FTuneven->freqlist_varname,"%s",argv[i]);
    VARTOOLS_RegisterDataVector(p, c, &(FTuneven->input_freqvar_vals), 
		       VARTOOLS_TYPE_DOUBLE,
		       0, VARTOOLS_SOURCE_EVALEXPRESSION_LC, 0, NULL, 
		       FTuneven->freqlist_varname, NULL);
  } else if(!strcmp(argv[i],"freqfile")) {
    FTuneven->freqtype = FTUNEVEN_FREQTYPE_FILE;
    i++;
    if(i >= argc) {
      *iret = i-1; return 1;
    }
    FTuneven->freqlist_filename = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
    sprintf(FTuneven->freqlist_filename,"%s",argv[i]);
    if((freqfile = fopen(FTuneven->freqlist_filename,"r")) == NULL) {
      fprintf(stderr,"Cannot open file %s\n",FTuneven->freqlist_filename);
      exit(ERR_CANNOTOPEN);
    }
    line = malloc(line_size);
    FTuneven->Nfreq_file = 0;
    FTuneven->Nfreq_vecsize = 256;
    if((FTuneven->freqvals_file = (double *) malloc(FTuneven->Nfreq_vecsize*sizeof(double))) == NULL) {
      fprintf(stderr,"Memory Allocation Error in reading the frequency file %s\n", FTuneven->freqlist_filename);
      exit(ERR_MEMALLOC);
    }
    while(VARTOOLS_gnu_getline(&line,&line_size,freqfile) >= 0) {
      k = 0;
      while(line[k] == ' ' || line[k] == '\t') k++;
      if(line[k] != '#') {
	if(FTuneven->Nfreq_file >= FTuneven->Nfreq_vecsize) {
	  FTuneven->Nfreq_vecsize *= 2;
	  if((FTuneven->freqvals_file = (double *) realloc(FTuneven->freqvals_file, FTuneven->Nfreq_vecsize*sizeof(double))) == NULL) {
	    fprintf(stderr,"Memory Allocation Error in reading the frequency file %s\n", FTuneven->freqlist_filename);
	    exit(ERR_MEMALLOC);
	  }
	}
	sscanf(&(line[k]),"%lf",&(FTuneven->freqvals_file[FTuneven->Nfreq_file]));
	FTuneven->Nfreq_file++;
      }
    }
    free(line);
    fclose(freqfile);
  } else {
    *iret = i; return 1;
  }

  i++;

  FTuneven->ft_sign = -1.0;
  check = VARTOOLS_ParseConstantParameter(p, c, &i, argv, argc, "ft_sign", VARTOOLS_TYPE_DOUBLE, &(FTuneven->ft_sign), 0);
  if(check == 2) {
    *iret = i; return 1;
  }

  if(FTuneven->ft_sign < 0.0) FTuneven->ft_sign = -1.0;
  else FTuneven->ft_sign = 1.0;

  FTuneven->tt_zero = 0.0;
  check = VARTOOLS_ParseConstantParameter(p, c, &i, argv, argc, "tt_zero", VARTOOLS_TYPE_DOUBLE, &(FTuneven->tt_zero), 0);
  if(check == 2) {
    *iret = i; return 1;
  }

  FTuneven->changeinputvectors = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"changeinputvectors")) {
      FTuneven->changeinputvectors = 1;
      i++;
      if(i >= argc) {
	*iret - i-1; return 1;
      }
      FTuneven->input_tvar_varname = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
      sprintf(FTuneven->input_tvar_varname,"%s",argv[i]);
      VARTOOLS_RegisterDataVector(p, c, &(FTuneven->input_tvar_vals), 
				  VARTOOLS_TYPE_DOUBLE,
				  0, VARTOOLS_SOURCE_EVALEXPRESSION_LC, 0, NULL, 
				  FTuneven->input_tvar_varname, NULL);
      i++;
      if(i >= argc) {
	*iret - i-1; return 1;
      }
      FTuneven->input_real_varname = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
      sprintf(FTuneven->input_real_varname,"%s",argv[i]);
      VARTOOLS_RegisterDataVector(p, c, &(FTuneven->input_real_vals), 
				  VARTOOLS_TYPE_DOUBLE,
				  0, VARTOOLS_SOURCE_EVALEXPRESSION_LC, 0, NULL, 
				  FTuneven->input_real_varname, NULL);
      i++;
      if(i >= argc) {
	*iret - i-1; return 1;
      }
      FTuneven->input_imag_varname = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
      sprintf(FTuneven->input_imag_varname,"%s",argv[i]);
      VARTOOLS_RegisterDataVector(p, c, &(FTuneven->input_imag_vals), 
				  VARTOOLS_TYPE_DOUBLE,
				  0, VARTOOLS_SOURCE_EVALEXPRESSION_LC, 0, NULL, 
				  FTuneven->input_imag_varname, NULL);
      i++;
    }
  }
  *iret = i;
  return 0;
}

void ftuneven_ShowSyntax(FILE *outfile)
/* Write out the expected syntax for this command to the file outfile */
{
  fprintf(outfile,
   "-ftuneven < \"outputvectors\" frequency_vecname FT_real_vecname\n"
   "			FT_imag_vecname Periodogram_vecname |\n"
   "            \"outputfile\" outdir [\"nameformat\" fmt] |\n"
   "            \"outputvectorsandfile\" outdir [\"nameformat\" fmt]\n"
   "                    frequency_vecname FT_real_vecname FT_imag_vecname\n"
   "                    Periodogram_vecname >\n"
   "          < \"freqauto\" |\n"
   "            \"freqrange\"\n"
   "                \"minfreq\" <\"fix\" value | \"list\" [\"column\" col] |\n"
   "                           \"fixcolumn\" <colname | colnum> |\n"
   "                           \"expr\" expression>\n"
   "                \"maxfreq\" <\"fix\" value | \"list\" [\"column\" col] |\n"
   "                           \"fixcolumn\" <colname | colnum> |\n"
   "                           \"expr\" expression>\n"
   "                \"freqstep\" <\"fix\" value | \"list\" [\"column\" col] |\n"
   "                            \"fixcolumn\" <colname | colnum> |\n"
   "                            \"expr\" expression>\n"
   "            \"freqvariable\" varname |\n"
   "            \"freqfile\" filename >\n"
   "          [\"ft_sign\" val] [\"tt_zero\" val]\n"
   "          [\"changeinputvectors\" tvec data_real_vec data_imag_vec]\n"
	  );
}

void ftuneven_ShowHelp(FILE *outfile)
/* Output the help information for this command */
{
  fprintf(outfile,
"Compute the complex Fourier transform of a time series with arbitrary spacing using a method developed by Jeff Scargle with support from NASA grant NNX16AL02G. The routine will return both the real and imaginary components of the Fourier transform, together with the absolute square of the result which is equivalent to the Lomb-Scargle periodogram. NOTE INPUT AND OUTPUT FREQUENCIES ARE IN UNITS OF RADIANS PER UNIT TIME.\n"
"The user must specify how these results will be returned. The following options are available:\n"
"\n\t\"outputvectors\" - return the results in a set of light curve vectors\n"
"\t\tthat can be used in subsequent vartools commands. The user\n"
"\t\tmust give the names for the vectors storing the frequency,\n"
"\t\tthe real component of the Fourier transform, the imaginary\n"
"\t\tcomponent, and the Lomb-Scargle periodogram. Note that\n"
"\t\tall light curve vectors will grow or shrink to match the\n"
"\t\tsize of the Fourier Transform vectors. Other vectors like\n"
"\t\tt, mag, or err may not be useful after if the Fourier\n"
"\t\ttransform is a different length.\n"
"\n\t\"outputfile\" - output the results to file. The user should specify\n"
"\t\tthe directory to write the file to. By default the output file\n"
"\t\twill have the name outdir/BASELC_NAME.ftuneven where\n"
"\t\tBASELC_NAME is the basename of the input light curve.\n"
"\t\tThe user may, however, give the \"nameformat\" keyword\n"
"\t\tfollowed by a format string to specify arbitrary filenames.\n"
"\t\tThe syntax is the same as for the \"-o\" VARTOOLS command.\n"
"\n\t\"outputvectorsandfile\" - output the results to a file, and to a set\n"
"\t\tof light curve vectors. First give the directory and possible\n"
"\t\tnameformat as in the \"outputfile\" option, and then give\n"
"\t\tthe list of vector names as in the \"outputvectors\" option.\n"
"\nThe user must then specify how the frequencies will be calculated. Options are:\n"
"\n\t\"freqauto\" - determine the frequencies automatically.\n"
"\n\t\"freqrange\" - use a uniformly spaced frequency grid with a specified\n"
"\t\trange and stepsize. Give the \"minfreq\", \"maxfreq\" and\n"
"\t\t\"freqstep\" parameters, in each case following standard\n"
"\t\tvartools syntax to indicate how the parameter should be set.\n"
"\n\t\"freqvariable\" - use an existing light curve vector for the\n"
"\t\tfrequencies. The user must give the name of the vector storing.\n"
"\t\tthe frequencies.\n"
"\n\t\"freqfile\" - read the frequencies from a file specified by the user.\n"
"\t\tThe file should be in white-space-delimited ascii format,\n"
"\t\twith the frequencies given in the first column of the file.\n"
"\t\tThe frequenceis read from this file will be used for all light\n"
"\t\tcurves processed.\n"
"\nThe user may set two optional parameters as well. These include \"ft_sign\", the sign of the Fourier transform. The default is -1 corresponding to a forward transform, but the user may change it +1 to carry out an inverse transform. The other parameter is \"tt_zero\", which is the origin of time for the transform. The default value is 0.\n"
"\nBy default the routine will apply the Fourier Transform to the data stored in the mag vector evaluated at time t. You can optionally use the \"changeinputvectors\" keyword to specify a different vector to use for t, and vectors to use for the real and imaginary components of the input data to be transformed. By using this option, together with the ft_sign option, you can perform an inverse Fourier transform to convert frequency-domain data back into time-domain data.\n"
"\nIf you use this routine please cite Scargle, 1989, 343, 874\n"
	  );
}

void ftuneven_RunCommand(ProgramData *p, void *userdata, int lc_name_num, int lc_num)
/* 
   This is the function called by VARTOOLS to run the command on a
   light curve. Typically this function is a wrapper which translates the
   data-types provided by VARTOOLS into what is needed for the command.

   p = structure containing various general program data (in
       particular the light curves are accessed through this structure).

   userdata = pointer to the structure containing the command specific
       data (including control parameters and vectors to store output
       results).

   lc_name_num = this is the index to use to access the light curve name.

   lc_num = this is the index to use to access the light curve and
   data from Registered vectors (e.g. the Period, T0, and other
   parameters for this command). In general this is not the same as
   lc_name_num. Its value will be 0 if VARTOOLS is reading a single
   light curve at a time and is not running in parallel, a number from
   0 to Nproc - 1 if VARTOOLS is running in parallel with Nproc processes,
   or it will equal lc_name_num if all light curves have been read-in.

*/
{
  int i;
  int NJD;
  double domega;
  double *t, *inputdata_real, *inputdata_imag;
  char *lcname;
  int Nfreq;
  double *freqvals;
  double *FT_real;
  double *FT_imag;
  double *FT_periodogram;

  FILE *outfile;

  /* Cast the input userdata pointer to the data structure type used
     for this command */
  _FTuneven *FTuneven;

  char lcoutname[MAXLEN];

  void do_ftuneven(int, double *, double *, double, double, int, double *, double *, double *, double *);

  void do_ftuneven_complexinput(int, double *, double *, double *, double, double, int, double *, double *, double *, double *);


  FTuneven = (_FTuneven *) userdata;

    /* The number of points in the light curve, and vectors storing the
     time, magnitudes, and errors are accessed as illustrated
     below. Note that t, mag, and err are vectors, t[0] is the time of
     the first observation in the light curve, t[NJD-1] is the time of
     the last observation, etc. */
  NJD=p->NJD[lc_num];

  if(!FTuneven->changeinputvectors) {
    t = p->t[lc_num];
    inputdata_real = p->mag[lc_num];
  } else {
    t = FTuneven->input_tvar_vals[lc_num];
    inputdata_real = FTuneven->input_real_vals[lc_num];
    inputdata_imag = FTuneven->input_imag_vals[lc_num];
  }

  /* If the light curve is too short, just quit */
  if(NJD <= 1) return;

  /* Access the light curve name as follows: */
  lcname = p->lcnames[lc_name_num];

  /* This command allows the user to optionally output model light curves.
     It uses a "standard" VARTOOLS syntax for specifying the output name,
     which is to provide a directory for the output, and then by default
     set the name to something like outdir/BASELC_NAME.suffix unless the
     user gives a format keyword, in which case the user can provide an
     arbitrary name. Any command using this convention can call the
     VARTOOLS_GetOutputFilename function to get the output filename. */
  if(FTuneven->outputtype == FTUNEVEN_OUTPUTFILE ||
     FTuneven->outputtype == FTUNEVEN_OUTPUTVECTORSANDFILE)
    VARTOOLS_GetOutputFilename(lcoutname,lcname,FTuneven->FT_outdir,"ftuneven",
			       FTuneven->FT_nameformat, lc_name_num);
  else 
    lcoutname[0] = '\0';

  /* Setup the frequencies to compute the Fourier Transform at */
  if(FTuneven->freqtype == FTUNEVEN_FREQTYPE_AUTO) {
    if((t[NJD-1] - t[0]) <= 0.0) return;
    domega = 2.0*M_PI/(t[NJD-1]-t[0]);
    Nfreq = 2*NJD + 1;
    if((freqvals = (double *) malloc(Nfreq*sizeof(double))) == NULL) {
      fprintf(stderr,"Memory allocation error for auto frequency initialization in ftuneven_RunCommand for lc_name_num=%d lc_num=%d\n", lc_name_num, lc_num);
      exit(ERR_MEMALLOC);
    }
    freqvals[0] = -NJD*domega;
    for(i=1; i < Nfreq; i++) {
      freqvals[i] = freqvals[i-1] + domega;
    }
  } else if(FTuneven->freqtype == FTUNEVEN_FREQTYPE_RANGE) {
    if((FTuneven->maxfreq[lc_num] <= FTuneven->minfreq[lc_num]) ||
       (FTuneven->freqstep[lc_num] <= 0.0)) return;
    Nfreq = floor((FTuneven->maxfreq[lc_num] - FTuneven->minfreq[lc_num])/FTuneven->freqstep[lc_num])+1;
    if((freqvals = (double *) malloc(Nfreq*sizeof(double))) == NULL) {
      fprintf(stderr,"Memory allocation error for range frequency initialization in ftuneven_RunCommand for lc_name_num=%d lc_num=%d\n", lc_name_num, lc_num);
      exit(ERR_MEMALLOC);
    }
    freqvals[0] = FTuneven->minfreq[lc_num];
    for(i=1; i < Nfreq; i++) {
      freqvals[i] = freqvals[i-1] + FTuneven->freqstep[lc_num];
    }
  } else if(FTuneven->freqtype == FTUNEVEN_FREQTYPE_VARIABLE) {
    Nfreq = NJD;
    if((freqvals = (double *) malloc(Nfreq*sizeof(double))) == NULL) {
      fprintf(stderr,"Memory allocation error for variable frequency initialization in ftuneven_RunCommand for lc_name_num=%d lc_num=%d\n", lc_name_num, lc_num);
      exit(ERR_MEMALLOC);
    }
    for(i=0; i < Nfreq; i++) {
      freqvals[i] = FTuneven->input_freqvar_vals[lc_num][i];
    }
  } else if(FTuneven->freqtype == FTUNEVEN_FREQTYPE_FILE) {
    Nfreq = FTuneven->Nfreq_file;
    if((freqvals = (double *) malloc(Nfreq*sizeof(double))) == NULL) {
      fprintf(stderr,"Memory allocation error for file frequency initialization in ftuneven_RunCommand for lc_name_num=%d lc_num=%d\n", lc_name_num, lc_num);
      exit(ERR_MEMALLOC);
    }
    for(i=0; i < Nfreq; i++) {
      freqvals[i] = FTuneven->freqvals_file[i];
    }
  } else
    return;

  if((FT_real = (double *) malloc(Nfreq*sizeof(double))) == NULL ||
     (FT_imag = (double *) malloc(Nfreq*sizeof(double))) == NULL ||
     (FT_periodogram = (double *) malloc(Nfreq*sizeof(double))) == NULL) {
      fprintf(stderr,"Memory allocation error ftuneven_RunCommand for lc_name_num=%d lc_num=%d\n", lc_name_num, lc_num);
      exit(ERR_MEMALLOC);
  }

  /* Compute the Fourier Transform */
  if(!FTuneven->changeinputvectors) {
    do_ftuneven(NJD, t, inputdata_real, FTuneven->ft_sign, FTuneven->tt_zero, Nfreq, freqvals, FT_real, FT_imag, FT_periodogram);
  } else {
    do_ftuneven_complexinput(NJD, t, inputdata_real, inputdata_imag, FTuneven->ft_sign, FTuneven->tt_zero, Nfreq, freqvals, FT_real, FT_imag, FT_periodogram);
  }

  /* Re-size the light curves if needed and store the output vectors */
  if(FTuneven->outputtype == FTUNEVEN_OUTPUTVECTORS ||
     FTuneven->outputtype == FTUNEVEN_OUTPUTVECTORSANDFILE) {
    if(Nfreq > NJD) {
      VARTOOLS_MemAllocDataFromLightCurveMidProcess(p, lc_num, Nfreq);
      p->NJD[lc_num] = Nfreq;
    } else if(Nfreq < NJD) {
      p->NJD[lc_num] = Nfreq;
    }
    if(Nfreq != NJD && p->readimagestring) {
      for(i=0;i<p->NJD[lc_num];i++)
	p->stringid_idx[lc_num][i] = i;
      VARTOOLS_mysortstringint(p->NJD[lc_num],MAXIDSTRINGLENGTH, p->stringid[lc_num], p->stringid_idx[lc_num]);
    }
    for(i=0; i < Nfreq; i++) {
      FTuneven->frequency_vec[lc_num][i] = freqvals[i];
      FTuneven->FT_real_vec[lc_num][i] = FT_real[i];
      FTuneven->FT_imag_vec[lc_num][i] = FT_imag[i];
      FTuneven->Periodogram_vec[lc_num][i] = FT_periodogram[i];
    }
  }
  
  /* Output the results to the file if requested */
  if(lcoutname[0] != '\0') {
    if((outfile = fopen(lcoutname,"w")) == NULL) {
      fprintf(stderr,"Error: cannot write to the file %s in the function ftuneven_RunCommand for lc_name_num=%d lc_num=%d\n", lcoutname, lc_name_num, lc_num);
      exit(ERR_CANNOTWRITE);
    }
    fprintf(outfile,"#Freq[Radians/time] FT_Real FT_Imag LS_Periodogram\n");
    for(i=0; i < Nfreq; i++) {
      fprintf(outfile,"%.17g %.17g %.17g %.17g\n", freqvals[i], FT_real[i], FT_imag[i], FT_periodogram[i]);
    }
    fclose(outfile);
  }
  
  free(freqvals); free(FT_real); free(FT_imag); free(FT_periodogram);

  return;
}

/* This is a translation to C of a matlab code written and provided by
   Jeff Scargle. The routine computes the Fourier Transform of
   arbitrarily spaced data using a modified version of the algorithm
   in ApJ 343, 1989, 874-887, Paper III. It also computes the
   Lomb-Scargle periodogram. */
void do_ftuneven(int NJD, double *t, double *mag, double ft_sign, double tt_zero, int Nfreq, double *freqvals, double *FT_real, double *FT_imag, double *FT_periodogram) {
  int ii_ww, ii_t;
  double magmean;
  double wrun, csum, ssum, wtau, sumr, sumi, scos2, ssin2, cval, sval;
  double ft_real, ft_imag, phi_this;

  magmean = 0.0;
  for(ii_t = 0; ii_t < NJD; ii_t++) {
    magmean += mag[ii_t];
  }
  magmean = magmean/(double) NJD;

  for(ii_ww = 0; ii_ww < Nfreq; ii_ww++) {
    wrun = freqvals[ii_ww];
    if(wrun == 0.0) {
      FT_real[ii_ww] = 0.0;
      for(ii_t = 0; ii_t < NJD; ii_t++) {
	FT_real[ii_ww] += (mag[ii_t]-magmean);
      }
      FT_real[ii_ww] = FT_real[ii_ww] / sqrt(NJD);
      FT_imag[ii_ww] = 0.0;
      FT_periodogram[ii_ww] = FT_real[ii_ww]*FT_real[ii_ww];
    } else {
      csum = 0.0; ssum = 0.0;
      for(ii_t = 0; ii_t < NJD; ii_t++) {
	csum += cos(2.0*wrun*t[ii_t]);
	ssum += sin(2.0*wrun*t[ii_t]);
      }
      wtau = 0.5*atan2(ssum, csum);
      
      /* Sum over the samples */
      sumr = 0.0; sumi = 0.0; scos2 = 0.0; ssin2 = 0.0;
      for(ii_t = 0; ii_t < NJD; ii_t++) {
	cval = cos(wrun*t[ii_t] - wtau);
	sval = sin(wrun*t[ii_t] - wtau);
	sumr += (mag[ii_t]-magmean)*cval;
	sumi += (mag[ii_t]-magmean)*sval;
	scos2 += cval*cval;
	ssin2 += sval*sval;
      }
      if(ft_sign > 0) {
	/* This is the inverse transform, apply an extra 1/sqrt(2) */
	sumr = sumr/sqrt(2.0);
	sumi = sumi/sqrt(2.0);
      }
      ft_real = sumr / (sqrt(2.0)*sqrt(scos2));
      ft_imag = ft_sign*sumi/(sqrt(2.0)*sqrt(ssin2));
      phi_this = wtau - wrun*tt_zero;
      FT_real[ii_ww] = ft_real*cos(phi_this) - ft_imag*sin(phi_this);
      FT_imag[ii_ww] = ft_imag*cos(phi_this) + ft_real*sin(phi_this);
      FT_periodogram[ii_ww] = (sumr*sumr/scos2) + (sumi*sumi/ssin2);
    }
  }
}

void do_ftuneven_complexinput(int NJD, double *t, double *input_real, double *input_imag, double ft_sign, double tt_zero, int Nfreq, double *freqvals, double *FT_real, double *FT_imag, double *FT_periodogram) {
  int ii_ww, ii_t;
  double magmean;
  double wrun, csum, ssum, wtau, sumr, sumi, scos2, ssin2, cval, sval;
  double ft_real, ft_imag, phi_this;

  for(ii_ww = 0; ii_ww < Nfreq; ii_ww++) {
    wrun = freqvals[ii_ww];
    if(wrun == 0.0) {
      FT_real[ii_ww] = 0.0;
      for(ii_t = 0; ii_t < NJD; ii_t++) {
	FT_real[ii_ww] += input_real[ii_t];
	FT_imag[ii_ww] += input_imag[ii_t];
      }
      FT_real[ii_ww] = FT_real[ii_ww] / sqrt(NJD);
      FT_imag[ii_ww] = FT_imag[ii_ww] / sqrt(NJD);
      FT_periodogram[ii_ww] = FT_real[ii_ww]*FT_real[ii_ww] + FT_imag[ii_ww]*FT_imag[ii_ww];
    } else {
      csum = 0.0; ssum = 0.0;
      for(ii_t = 0; ii_t < NJD; ii_t++) {
	csum += cos(2.0*wrun*t[ii_t]);
	ssum += sin(2.0*wrun*t[ii_t]);
      }
      wtau = 0.5*atan2(ssum, csum);
      
      /* Sum over the samples */
      sumr = 0.0; sumi = 0.0; scos2 = 0.0; ssin2 = 0.0;
      for(ii_t = 0; ii_t < NJD; ii_t++) {
	cval = cos(wrun*t[ii_t] - wtau);
	sval = sin(wrun*t[ii_t] - wtau);
	sumr += input_real[ii_t]*cval - ft_sign*input_imag[ii_t]*sval;
	sumi += ft_sign*input_real[ii_t]*sval + input_imag[ii_t]*cval;
	scos2 += cval*cval;
	ssin2 += sval*sval;
      }
      if(ft_sign > 0) {
	/* This is the inverse transform, apply an extra 1/sqrt(2) */
	sumr = sumr/sqrt(2.0);
	sumi = sumi/sqrt(2.0);
      }
      ft_real = sumr / (sqrt(2.0)*sqrt(scos2));
      ft_imag = sumi/(sqrt(2.0)*sqrt(ssin2));
      phi_this = wtau - wrun*tt_zero;
      FT_real[ii_ww] = ft_real*cos(phi_this) - ft_imag*sin(phi_this);
      FT_imag[ii_ww] = ft_imag*cos(phi_this) + ft_real*sin(phi_this);
      FT_periodogram[ii_ww] = (sumr*sumr/scos2) + (sumi*sumi/ssin2);
    }
  }
}
