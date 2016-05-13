/*
   This file is part of the VARTOOLS implementation of David Palmer's fastchi2 
   period search algorithm. The reference for this algorithm 
       Palmer 2009, ApJ, 695, 496.
   This file contains functions needed to interface the routine with VARTOOLS,
   as well as a driver function. The bulk of the code carrying out the algorithm
   is included in the file fastchi2_lib.c, which is a copy of the fastchi2.c
   file included in Palmer's own reference implementation of this algorithm.

   The driver function in this case is inspired by the chi2driver.c file
   included in Palmer's reference implementation of fastchi2.  This file is
   a modified/derivative version of that file. Below is the copyright
   notice from that file:

**************************************************************************
COPYRIGHT NOTICE FROM chi2driver.c included in FastChi2-1.03, which is the 
basis for this file. THIS FILE IS A MODIFIED, DERIVATIVE VERSION OF PALMER'S
ORIGINAL SOFTWARE.
**************************************************************************

* fastchi driver program
* main() for running fast-chi algorithm described in
* Palmer, D.M., 2007 (submitted)
* "A Fast Chi-squared Technique For Period Search of Irregularly Sampled Data."
*
* Usage information follows the copyright notice.
*/
/*
 Copyright 2007.  Los Alamos National Security, LLC. This material was 
 produced under U.S. Government contract DE-AC52-06NA25396 for 
 Los Alamos National Laboratory (LANL), which is operated by 
 Los Alamos National Security, LLC for the U.S. Department of Energy. 
 The U.S. Government has rights to use, reproduce, and distribute this software.  
 NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY 
 WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF 
 THIS SOFTWARE.  If software is modified to produce derivative works, 
 such modified software should be clearly marked, so as not to confuse 
 it with the version available from LANL.
 
 
 Additionally, this program is free software; you can redistribute it 
 and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation; version 2.0 of the License. 
 Accordingly, this program is distributed in the hope that it will be useful, 
 but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 for more details
*/

#include "../../src/vartools.h"
#include <stdio.h>
#include <stdlib.h>
#include "fastchi2.h"
#include "fastchi2_lib.h"

/* This library defines the command -fastchi2 which can be used to 
   apply Palmer's Fast Chi2 period search to a light curve.
*/

void fastchi2_Initialize(char *commandname,
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
  sprintf(commandname,"-fastchi2");

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
     the magnitudes of each lc), below you would replace "_Fastchi2"
     with the type name of your structure.

     See fastchi2.h
*/
  *sizeuserdata = sizeof(_Fastchi2);
}

int fastchi2_ParseCL(ProgramData *p, Command *c,
			   void *userdata, int *iret, char **argv, int argc)
/* Every library with the name $LIBNAME.so must contain a function
   with the name $LIBNAME_ParseCL to check that the user has issued a
   valid syntax, and set-up the command-data accordingly.

   Note that any vectors/arrays used by the command should also be registered
   by this function so that memory will be allocated as needed.

Expected syntax for this example:
   -fastchi2 
        <"Nharm" <"fix" val | "list" | "fixcolumn" <colname | colnum> | "expr" expr>
        <"freqmax" <"fix" val | "list" | "fixcolumn" <colname | colnum> | "expr" expr>
        ["freqmin" <"fix" val | "list" | "fixcolumn" <colname | colnum> | "expr" expr]
        ["detrendorder" <"fix" val | "list" | "fixcolumn" <colname | colnum> | "expr" expr]
        ["t0" <"fix" val | "list" | "fixcolumn" <colname | colnum> | "expr" expr]
        ["timespan" <"fix" val | "list" | "fixcolumn" <colname | colnum> | "expr" expr]
        ["oversample" <"fix" val | "list" | "fixcolumn" <colname | colnum> | "expr" expr]
        ["chimargin" <"fix" val | "list" | "fixcolumn" <colname | colnum> | "expr" expr]
	["Npeak" val]
        ["norefitpeak"]
        ["oper" outdir ["nameformat" format]]
        ["omodel" outdir ["nameformat" format]]

  p = structure containing general program data. In most cases this 
      structure is only passed to other functions (like RegisterDataVector)
      and not used directly by the ParseCL function.

  c = pointer to a generic container for this command. This structure is
      needed by other functions that might be called here (RegisterDataVector)
      but should not be used directly by the ParseCL function.

  userdata = pointer to the structure storing data for this
      command. This should be cast to the appropriate type and the
      data there-in should be modified as needed by the ParseCL
      function. In this example, userdata corresponds to the _Fastchi2
      structure which stores the values to be added to the light
      curves.

  iret = command line argument index. This is initially zero, and on
      return should be set to the index of the last argument parsed by
      this command. Note that index 0 points to the name of the
      command, which does not need to be checked (this function is
      only called if the user issues the command).

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
  _Fastchi2 *Fastchi2;
  Fastchi2 = (_Fastchi2 *) userdata;
  
  /* We can use the function ParseParameter to parse a line of the form
     "keyword" <"fix" value | "list" [\"column\" col] | "fixcolumn" <colname | colnum>>

     After giving p, c, iret, argv, argc, and the keyword to check
     for, give the number of data-vectors to be set by this
     command. For each data-vector give the data-type for the vector,
     a pointer to the vector, the number of columns in the array (if
     it is <= 0 then the pointer will be taken to be a vecotr,
     otherwise it is assumed to be an array), a flag indicating if
     this data should be included in the output ascii table, the root
     name of the vector for display in the output ascii table and for
     the input table (if the user gives the \"list\" keyword on the
     command line), a flag indicating if the vector should be
     initialized if the "keyword" is not given by the user (if the
     flag is set to 0, then memory will not be allocated for the
     vector/array if the keyword is not given), and if the flag is
     one, then the default value to fix the parameter to.

     The function will update iret.
     
     It will return 0 if the command was parsed successfully, 1 if the keyword
     was not given, and 2 if the keyword was given but the rest of the 
     command was not parsed successfully.
  */

  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "Nharm", 1,
				  VARTOOLS_TYPE_INT,
				  (void *) (&Fastchi2->Nharm), 0, 0,
				  "Nharm", 0);
  if(check) {*iret = i; return 1;}
  
  
  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "freqmax", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Fastchi2->freqmax),
				  0, 0, "Freqmax", 0);
  if(check) {*iret = i; return 1;}
 

  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "freqmin", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Fastchi2->freqmin),
				  0, 0, "Freqmin", 1, FASTCHI2_DEFAULT_FREQMIN);
  if(check == 2) {*iret = i; return 1;}
  else if(check == 1) Fastchi2->fixfreqmin=0;
  else Fastchi2->fixfreqmin=1;

 
  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "detrendorder", 1,
				  VARTOOLS_TYPE_INT,
				  (void *) (&Fastchi2->detrendorder),
				  0, 0, "DetrendOrder", 1, 
				  FASTCHI2_DEFAULT_DETRENDORDER);
  if(check == 2) {*iret = i; return 1;}
 
  
  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "t0", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Fastchi2->t0),
				  0, 0, "T0", 0);
  if(check == 2) {*iret = i; return 1;}
  else if(check == 1) Fastchi2->fixt0=0;
  else Fastchi2->fixt0=1;
      
  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "timespan", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Fastchi2->timespan),
				  0, 0, "Timespan", 0);
  if(check == 2) {*iret = i; return 1;}
  else if(check == 1) Fastchi2->fixtimespan=0;
  else Fastchi2->fixtimespan=1;
 
  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "oversample", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Fastchi2->oversample),
				  0, 0, "Oversample", 1,
				  FASTCHI2_DEFAULT_OVERSAMPLE);
  if(check == 2) {*iret = i; return 1;}
 
  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "chimargin", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Fastchi2->chimargin),
				  0, 0, "Chimargin", 1,
				  FASTCHI2_DEFAULT_CHIMARGIN);
  if(check == 2) {*iret = i; return 1;}
 
  check = VARTOOLS_ParseConstantParameter(p, c, &i, argv, argc, "Npeak",
					  VARTOOLS_TYPE_INT, 
					  &(Fastchi2->Npeak), 0);
  if(check == 2) {*iret = i; return 1;}
  else if(check == 1) {
    Fastchi2->Npeak = FASTCHI2_DEFAULT_NPEAK;
  }
 
  if(i < argc ? !strcmp(argv[i],"norefitpeak") : 0) {
    Fastchi2->dorefitpeak = 0;
    i++;
  } else {
    Fastchi2->dorefitpeak = 1;
  }
 
  check = VARTOOLS_ParseOutNameKeyword(p, c, &i, argv, argc, "oper",
				       &(Fastchi2->outputper),
				       &(Fastchi2->peroutdir),
				       &(Fastchi2->outputper_useformat),
				       &(Fastchi2->outputper_format));
  if(check == 2) {*iret = i; return 1;}
 
  check = VARTOOLS_ParseOutNameKeyword(p, c, &i, argv, argc, "omodel",
				       &(Fastchi2->outputmodel),
				       &(Fastchi2->modeloutdir),
				       &(Fastchi2->outputmodel_useformat),
				       &(Fastchi2->outputmodel_format));
  if(check == 2) {*iret = i; return 1;}

  /* Check if the user is attempting to store the output model in a variable */
  Fastchi2->saveoutputmodel = 0;
  if(i < argc ? !strcmp(argv[i],"omodelvariable") : 0) {
    i++;
    if(i >= argc) {*iret = i; return 1;}
    Fastchi2->saveoutputmodel = 1;
    VARTOOLS_RegisterDataVector(p, c, (void *) (&Fastchi2->outputmodeldata),
				VARTOOLS_TYPE_DOUBLE, 0, VARTOOLS_SOURCE_LC,
				0, NULL, NULL, -1, NULL, argv[i]);
    i++;
  }

  /* Register the variables that will store the output data */
  VARTOOLS_RegisterDataVector(p, c, (void *) (&Fastchi2->bestfreq),
			      VARTOOLS_TYPE_DOUBLE, Fastchi2->Npeak,
			      VARTOOLS_SOURCE_COMPUTED,
			      1, "Frequency");
  VARTOOLS_RegisterDataVector(p, c, (void *) (&Fastchi2->chireduction),
			      VARTOOLS_TYPE_FLOAT, Fastchi2->Npeak,
			      VARTOOLS_SOURCE_COMPUTED,
			      1, "Chi2Reduction");  
  *iret = i;
  return 0;
				       
}

void fastchi2_ShowSyntax(FILE *outfile)
/* Write out the expected syntax for this command to the file outfile */
{
  fprintf(outfile,
	  "-fastchi2\n"
	  "\t<\"Nharm\" <\"fix\" val | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
	  "\t<\"freqmax\" <\"fix\" val | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
	  "\t[\"freqmin\" <\"fix\" val | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr]\n"
	  "\t[\"detrendorder\" <\"fix\" val | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr]\n"
	  "\t[\"t0\" <\"fix\" val | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr]\n"
	  "\t[\"timespan\" <\"fix\" val | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr]\n"
	  "\t[\"oversample\" <\"fix\" val | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr]\n"
	  "\t[\"chimargin\" <\"fix\" val | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr]\n"
	  "\t[\"Npeak\" val]\n"
	  "\t[\"norefitpeak\"]\n"
	  "\t[\"oper\" outdir [\"nameformat\" format]]\n"
	  "\t[\"omodel\" outdir [\"nameformat\" format]]\n"
	  "\t[\"omodelvariable\" varname]\n");
}

void fastchi2_ShowHelp(FILE *outfile)
/* Output the help information for this command */
{
  fprintf(outfile,"Apply the Fastchi2 periodogram algorithm due to David Palmer (Palmer, D.M., 2009, ApJ, 695, 496). If you use this routine be sure to cite Palmer's paper. This implementation uses Palmer's code.  For each parameter given the keyword, and then the source of the parameter (use \"fix\" to fix it to a specific value for all light curves, \"list\" to specify the parameter in the input list, \"fixcolumn\" to take the value of the parameter from the output of a previous command, or \"expr\" to use an analytic expression to initialize the parameter). Nharm is the number of harmonics to use in the model (1 = just the fundamental, 2 = fundamental + first overtone, etc.). freqmax is the maximum frequency to search in cycles per day. These are the only two required parameters for this routine. freqmin is the minimum frequency to search (the default is 0.0). detrendorder is the order of a polynomial to use in detrending the light curve, before doing the period search (this is 0 by default). t0 can be used to fix the initial time for the detrending. timespan can be used to fix the timespan of the light curves (to use in calculating the Nyquist frequency). oversample is the factor by which the periodogram will be oversampled. chimargin is used in doing a finer search for the best frequencies around the peaks in the periodogram. Npeak is the number of peaks to find in the periodogram. If the \"norefitpeak\" keyword is given, then the fine search will not be carried out, and only the peaks in the calculated periodogram will be printed out. Use \"oper\" or \"omodel\" to output the periodogram and the harmonic function model, respectively.\n");
}

void fastchi2_ShowExample(FILE *outfile)
/* Output an example for this command */
{
  fprintf(outfile,
	  "\nvartools -i EXAMPLES/2 -oneline -ascii \\\n"
	  "\t-L USERLIBS/src/fastchi2.so \\\n"
	  "\t-fastchi2 Nharm fix 1 freqmax fix 10.0 \\\n"
	  "\tfreqmin fix 0.1 oper EXAMPLES/OUTDIR1/\n\n"
	  "Run Palmer's fast chi2 periodogram on the light curve EXAMPLES/2. Search up to a maximum frequency of 10 cycles per day and a minimum frequency of 0.1 cycles per day. Use one harmonic in the model, and output the periodogram to EXAMPLES/OUTDIR1.\n");
}

void fastchi2_RunCommand(ProgramData *p, void *userdata, int lc_name_num, int lc_num)
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
  int NJD, i, j, Ngood, nrealpoints;
  double *t, *mag, *err;
  double *tgood = NULL;
  double *magcorr_dbl = NULL;
  float *err_flt = NULL, *magcorr_flt = NULL;
  float *pbest_peakup = NULL;
  float *data_var = NULL, *one_var = NULL, *chireductionvalues = NULL;
  double *detrendcoeffs = NULL;
  char *lcname;
  double addval;
  _Fastchi2 *Fastchi2;
  char peroutname[MAXLEN];
  char modeloutname[MAXLEN];

  int Nharm, detrendorder, Npeak, dorefitpeak;
  double freqmax, freqmin;
  float **harmweights = NULL;
  double oversample, chimargin;

  FILE *outper, *outmodel;

  double obsrange, t0, tstart;
  double mintotalrange, deltat, totalrange, deltaf;
  double chi0, mean0, chid, v, v2;
  int nout, nchivalues;

  /* Translate pointers and structures input by vartools into
     easier to use forms */
  Fastchi2 = (_Fastchi2 *)userdata;


  /* The number of points in the light curve, and vectors storing the
     time, magnitudes, and errors are accessed as illustrated
     below. Note that t, mag and err are vectors, t[0] is the time of
     the first observation in the light curve, t[NJD-1] is the time of
     the last observation, etc. */
  NJD=p->NJD[lc_num];
  t = p->t[lc_num];
  mag = p->mag[lc_num];
  err = p->sig[lc_num];

  /* Get the parameters to use for this light curve */
  Nharm = Fastchi2->Nharm[lc_num];
  freqmax = Fastchi2->freqmax[lc_num];
  freqmin = Fastchi2->freqmin[lc_num];
  detrendorder = Fastchi2->detrendorder[lc_num];
  oversample = Fastchi2->oversample[lc_num];
  chimargin = Fastchi2->chimargin[lc_num];
  Npeak = Fastchi2->Npeak;
  dorefitpeak = Fastchi2->dorefitpeak;

  /* Allocate memory for the float versions of the data arrays */
  if(NJD > 0) {
    if((tgood = (double *) malloc(NJD * sizeof(double))) == NULL ||
       (magcorr_dbl = (double *) malloc(NJD * sizeof(double))) == NULL ||
       (err_flt = (float *) malloc(NJD * sizeof(float))) == NULL ||
       (magcorr_flt = (float *) malloc(NJD * sizeof(float))) == NULL ||
       (detrendcoeffs = (double *) malloc((1 + detrendorder)*sizeof(double))) == NULL ||
       (harmweights = (float **) malloc(Npeak * sizeof(float *))) == NULL ||
       (pbest_peakup = (float *) malloc(Npeak * sizeof(float))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
  }
  for(i=0; i < Npeak; i++) {
    if((harmweights[i] = (float *) malloc((2*Nharm + 1) * sizeof(float))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
  }

  /* Get the names of the output files if output will be done */
  if(Fastchi2->outputper) {
    if(!Fastchi2->outputper_useformat) {
      VARTOOLS_GetOutputFilename(peroutname, p->lcnames[lc_name_num],
				 Fastchi2->peroutdir, "fastchi2_per",
				 NULL, lc_name_num);
    } else {
      VARTOOLS_GetOutputFilename(peroutname, p->lcnames[lc_name_num],
				 Fastchi2->peroutdir, "fastchi2_per",
				 Fastchi2->outputper_format, lc_name_num);
    }
  }
  if(Fastchi2->outputmodel) {
    if(!Fastchi2->outputmodel_useformat) {
      VARTOOLS_GetOutputFilename(modeloutname, p->lcnames[lc_name_num],
				 Fastchi2->modeloutdir, "fastchi2_model",
				 NULL, lc_name_num);
    } else {
      VARTOOLS_GetOutputFilename(modeloutname, p->lcnames[lc_name_num],
				 Fastchi2->modeloutdir, "fastchi2_model",
				 Fastchi2->outputmodel_format, lc_name_num);
    }
  }

  /* Copy the magnitude errors from the input double vectors into the
     temporary float vectors */
  Ngood = 0;
  for(i=0; i < NJD; i++) {
    magcorr_dbl[i] = (double) mag[i];
  }
  
  /* Get the values of deltat, deltaf, data_var, one_var, and chireduction. */
  if(!Fastchi2->fixtimespan) {
    obsrange = t[NJD-1] - t[0];
  } else {
    obsrange = Fastchi2->timespan[lc_num];
  }
  if(!Fastchi2->fixt0) {
    t0 = 0.5*(t[NJD-1]+t[0]);
  } else {
    t0 = Fastchi2->t0[lc_num];
  }

  mintotalrange = obsrange * Nharm * oversample * 2;
  deltat = 1/(2*freqmax * 2*Nharm);

  nrealpoints = 1L << (int)(ceil(log2(mintotalrange/deltat)));
  if (nrealpoints <= 1) 
    {
      fprintf(stderr,"Error: too many time bins (%.0f) in FFT calculation in fastchi2\n",mintotalrange/deltat);
      exit(1);
    }

  totalrange = nrealpoints * deltat;
  deltaf = 1.0/totalrange;

  nchivalues = nrealpoints / (2 * 2 * Nharm);

  if((data_var = (float *) malloc(nrealpoints * sizeof(float))) == NULL ||
     (one_var = (float *) malloc(nrealpoints * sizeof(float))) == NULL ||
     (chireductionvalues = (float *) malloc(nchivalues * sizeof(float))) == NULL) {
    VARTOOLS_error(ERR_MEMALLOC);
  }

  if(NJD < detrendorder + 2*Nharm + 2) {
    /* Not enough points in this lc, skip it */
    if(tgood != NULL)
      free(tgood);
    if(magcorr_dbl != NULL)
      free(magcorr_dbl);
    if(err_flt != NULL)
      free(err_flt);
    if(magcorr_flt != NULL)
      free(magcorr_flt);
    if(data_var != NULL)
      free(data_var);
    if(one_var != NULL)
      free(one_var);
    if(chireductionvalues != NULL)
      free(chireductionvalues);
    if(detrendcoeffs != NULL)
      free(detrendcoeffs);
    if(harmweights != NULL) {
      for(i=0; i < Npeak; i++) {
	if(harmweights[i] != NULL)
	  free(harmweights[i]);
      }
      free(harmweights);
    }
    if(pbest_peakup != NULL)
      free(pbest_peakup);

    for(i=0; i < Npeak; i++) {
      Fastchi2->bestfreq[lc_num][i] = -1.0;
      Fastchi2->chireduction[lc_num][i] = 0.0;
    }
    return;
  }

  /* Get chi2 for a straaight line, and the weighted mean mag */

  chi0 = VARTOOLS_chi2(NJD, t, mag, err, &mean0, &nout);

  /* Subtract a polynomial fit to the LC */

  for(i=0; i < NJD; i++) {
    tgood[i] = t[i] - t0;
  }

  chid = VARTOOLS_fitpoly(NJD, tgood, magcorr_dbl, err, detrendorder, 
			  1, detrendcoeffs, NULL);

  Ngood = 0;
  for(i=0; i < NJD; i++) {
    if(!isnan(mag[i]) && err[i] > 0.0) {
      magcorr_flt[Ngood] = (float) magcorr_dbl[i];
      tgood[Ngood] = tgood[i];
      err_flt[Ngood] = (float) err[i];
      Ngood++;
    }
  }

  /* Run the fastchi2 algorithm */
  
  tstart = tgood[0] - deltat/2.0;  
  fastChi(tgood, magcorr_flt, err_flt, Ngood,
	  Nharm, deltat, nrealpoints, nchivalues,
	  tstart, data_var, one_var, chireductionvalues);

  /* Output the periodogram if asked to */
  if(Fastchi2->outputper) {
    if((outper = fopen(peroutname,"w")) == NULL) {
      VARTOOLS_error2(ERR_CANNOTWRITE,peroutname);
    }
    fprintf(outper,"#Frequency  Chireduction\n");
    for(i = freqmin/deltaf; i < nchivalues; i++) {
      fprintf(outper, "%.17g %.17g\n", i*deltaf, chireductionvalues[i]);
    }
    fclose(outper);
  }

  findBestFrequency(dorefitpeak, tgood, magcorr_flt, err_flt,
		    Ngood, Nharm, freqmin, deltaf,
		    tstart, chimargin, chireductionvalues,
		    nchivalues, pbest_peakup,
		    Fastchi2->bestfreq[lc_num], Npeak);

  for(i=0; i < Npeak; i++) {
    Fastchi2->chireduction[lc_num][i] = 
      singleChi(tgood, magcorr_flt, err_flt,
		Ngood, Nharm, tstart, 
		Fastchi2->bestfreq[lc_num][i],
		harmweights[i]);
  }

  /* Output the model if asked to */
  if(Fastchi2->outputmodel || Fastchi2->saveoutputmodel) {
    if(Fastchi2->outputmodel) {
      if((outmodel = fopen(modeloutname,"w")) == NULL) {
	VARTOOLS_error2(ERR_CANNOTWRITE,modeloutname);
      }
      /* First print the model function */
      fprintf(outmodel,"# Best-fit periodic function determined by the fastchi2 periodogram:\n");
      fprintf(outmodel,"# m(t) = ");
      
      v = harmweights[0][0];
      if(detrendorder >= 0) {
	v += detrendcoeffs[0];
      }
      fprintf(outmodel,"%.17g",v);
      for(i=1; i <= detrendorder; i++) {
	if(i == 1)
	  fprintf(outmodel," + %.17g*(t-t0)", detrendcoeffs[i]);
	else
	  fprintf(outmodel," + %.17g*(t-t0)^%d", detrendcoeffs[i], i);
      }
      for(i=1; i <= Nharm; i++) {
	if(i == 1)
	  fprintf(outmodel," + %.17g*sin(2*pi*(t-t00)*f) + %.17g*cos(2*pi*(t-t00)*f)", harmweights[0][2*i-1], harmweights[0][2*i]);
	else
	  fprintf(outmodel," + %.17g*sin(2*pi*(t-t00)*%d*f) + %.17g*cos(2*pi*(t-t00)*%d*f)", harmweights[0][2*i-1], i, harmweights[0][2*i], i);
      }
      fprintf(outmodel,"\n");
      fprintf(outmodel,"#f=%.17g  t0=%.17g  t00=%.17g\n", Fastchi2->bestfreq[lc_num][0], t0, t0+tstart);
      fprintf(outmodel,"#\n");
      fprintf(outmodel,"#Time  Mag_obs  Mag_model  Error\n");
    }
    /* Now print evaluate and print out the model */
    for(j=0; j < NJD; j++) {
      if((isnan(mag[j]) || err[j] <= 0.0) && !Fastchi2->saveoutputmodel)
	continue;
      v = harmweights[0][0];
      if(detrendorder >= 0) {
	v += detrendcoeffs[0];
      }
      v2 = 1;
      for(i=1; i <= detrendorder; i++) {
	v2 = v2 * (t[j] - t0);
	v += detrendcoeffs[i]*v2;
      }
      for(i=1; i <= Nharm; i++) {
	v += harmweights[0][2*i-1]*sin(2*M_PI*(t[j]-t0-tstart)*Fastchi2->bestfreq[lc_num][0]*i) + harmweights[0][2*i]*cos(2*M_PI*(t[j]-t0-tstart)*Fastchi2->bestfreq[lc_num][0]*i);
      }
      if(Fastchi2->outputmodel && !isnan(mag[j]) && err[j] > 0.0) {
	fprintf(outmodel,"%.17g %.17g %.17g %.17g\n", t[j], mag[j], v, err[j]);
      }
      if(Fastchi2->saveoutputmodel) {
	Fastchi2->outputmodeldata[lc_num][j] = v;
      }
    }
    if(Fastchi2->outputmodel)
      fclose(outmodel);
  }

  VARTOOLS_sort_generic(Npeak, 1, NULL, 2, 
			VARTOOLS_TYPE_FLOAT,
			(void *) Fastchi2->chireduction[lc_num],
			0,
			VARTOOLS_TYPE_DOUBLE,
			(void *) Fastchi2->bestfreq[lc_num],
			0);

  if(tgood != NULL)
    free(tgood);
  if(magcorr_dbl != NULL)
    free(magcorr_dbl);
  if(err_flt != NULL)
    free(err_flt);
  if(magcorr_flt != NULL)
    free(magcorr_flt);
  if(data_var != NULL)
    free(data_var);
  if(one_var != NULL)
    free(one_var);
  if(chireductionvalues != NULL)
    free(chireductionvalues);
  if(detrendcoeffs != NULL)
    free(detrendcoeffs);
  if(harmweights != NULL) {
    for(i=0; i < Npeak; i++) {
      if(harmweights[i] != NULL)
	free(harmweights[i]);
    }
    free(harmweights);
  }
  if(pbest_peakup != NULL)
    free(pbest_peakup);

}

