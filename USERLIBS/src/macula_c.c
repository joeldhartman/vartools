/*
  This file is part of the VARTOOLS implementation of David Kipping's
  macula spot modelling routine. The reference for macula is
      Kipping 2012, arXiv:1209.2985

  The website for macula is https://www.cfa.harvard.edu/~dkipping/macula.html

  The main body of the routine is contained in macula.f90, is taken
  directly from this website (accessed on 2012.1017), with minor
  modifications to allow it to interface with VARTOOLS.

  This file provides the interface between VARTOOLS and macula.f90.

  Please contact David Kipping for questions about the algorithm, also be
  sure to reference his paper if you use it.

*/

#include "../../src/vartools.h"
#include <stdio.h>
#include <stdlib.h>
#include "macula.h"

void macula_Initialize(char *commandname,
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
  sprintf(commandname,"-macula");

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
     the magnitudes of each lc), below you would replace "_Macula"
     with the type name of your structure.

     See macula.h
*/
  *sizeuserdata = sizeof(_Macula);
}

int macula_ParseCL(ProgramData *p, Command *c,
		   void *userdata, int *iret, char **argv, int argc)
/* Every library with the name $LIBNAME.so must contain a function
   with the name $LIBNAME_ParseCL to check that the user has issued a
   valid syntax, and set-up the command-data accordingly.

   Note that any vectors/arrays used by the command should also be registered
   by this function so that memory will be allocated as needed.

Expected syntax for this example:
   -macula < "inject" | "fit" <"amoeba" | "lm"> >
       <"Prot" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
           ["vary"]>
       <"istar" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
           ["vary"]>
       <"kappa2" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
           ["vary"]>
       <"kappa4" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
           ["vary"]>
       <"c1" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
           ["vary"]>
       <"c2" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
           ["vary"]>
       <"c3" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
           ["vary"]>
       <"c4" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
           ["vary"]>
       <"d1" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
           ["vary"]>
       <"d2" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
           ["vary"]>
       <"d3" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
           ["vary"]>
       <"d4" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
           ["vary"]>
       <"blend" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
           ["vary"]>
       <"Nspot" value>
          ... Repeat below Nspot times ....
          <"Lambda0" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
             ["vary"]> 
          <"Phi0" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
             ["vary"]> 
          <"alphamax" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
             ["vary"]> 
          <"fspot" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
             ["vary"]> 
          <"tmax" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
             ["vary"]> 
          <"life" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
             ["vary"]> 
          <"ingress" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
             ["vary"]> 
          <"egress" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
             ["vary"]> 
          .... Stop repeating here ....
      ["fluxinput"] ["fluxoutput"] ["correctlc"]
      ["omodel" <outdir ["nameformat" fmt]> ["tdelv"]]
      ["ocurve" <"outdir" outdir [ "nameformat" fmt]>
          ["tdelv"] ["step" stepsize]]

  p = structure containing general program data. In most cases this 
      structure is only passed to other functions (like RegisterDataVector)
      and not used directly by the ParseCL function.

  c = pointer to a generic container for this command. This structure is
      needed by other functions that might be called here (RegisterDataVector)
      but should not be used directly by the ParseCL function.

  userdata = pointer to the structure storing data for this
      command. This should be cast to the appropriate type and the
      data there-in should be modified as needed by the ParseCL
      function. In this example, userdata corresponds to the _Jktebop
      structure which stores the values to be added to the light
      curves.

  iret = command line argument index. Initially argv[*iret - 1] is the
      commandname ("-jktebop" in this case; this function will only be
      called if the user issues the command so this does not need to
      be verified by this function). On output *iret should be
      incremented by the number of terms parsed from the
      command-line. In this example, if the user issues a command like
      "-jktebop inject Period fix 1.0 T0 fix 0.0 r1+r2 fix 0.1 \
           r2/r1 fix 1.0 M2/M1 fix 1.0 J2/J1 fix 1.0 bimpact fix 0.0 \
           esinomega fix 0.0 ecosomega fix 0.0 LD1 linear fix 0.5 \
           LD2 lockLD1"
      The function will parse 34 terms and on output *iret should be
      incremented by 34.

  argv = array of command line arguments. argv[*iret - 1] is the name
      of the command ("-jktebop" in this case).

  argc = Number of command line arguments in the argv array (including
      0). It is the user's responsibility to check that i < argc
      before attempting to parse argv[i], failure to do so may lead to
      segmentation violations.
 */
{
  int i = 0, jspot;
  int check;

  char tmpstr[256];

  _Macula *Macula;
  Macula = (_Macula *) userdata;

  /* We'll use i rather than iret to index argv, this is just so we
     don't have to constantly dereference the pointer iret */
  i = *iret;

  /* Check if "inject" or "fit" is set */
  if(i >= argc) {
    /* There are no terms left on the command-line, return an error */
    return 1;
  }
  if(!strcmp(argv[i],"inject")) {
    Macula->injectorfit = MACULA_INJECT;
  } 
  else if(!strcmp(argv[i],"fit")) {
    Macula->injectorfit = MACULA_FIT;
    i++;
    if(i >= argc) {
      *iret = i;
      return 1;
    }
    if(!strcmp(argv[i],"amoeba")) {
      Macula->fittype = MACULA_FITTYPE_AMOEBA;
    }
    else if(!strcmp(argv[i],"lm")) {
      Macula->fittype = MACULA_FITTYPE_LM;
    }
    else {
      *iret = i;
      return 1;
    }
  }
  else {
    return 1;
  }

  i++;
  
  /* Now parse the options for each of the fitting parameters */
  
  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "Prot", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Macula->Prot), 0, 1,
				  "Prot", 0);
  if(check) {*iret = i; return 1;}
  Macula->Prot_vary = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Macula->Prot_vary = 1;
      i++;
    }
  }

  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "istar", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Macula->Istar), 0, 1,
				  "Istar", 0);
  if(check) {*iret = i; return 1;}
  Macula->Istar_vary = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Macula->Istar_vary = 1;
      i++;
    }
  }

  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "kappa2", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Macula->kappa2), 0, 1,
				  "kappa2", 0);
  if(check) {*iret = i; return 1;}
  Macula->kappa2_vary = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Macula->kappa2_vary = 1;
      i++;
    }
  }

  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "kappa4", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Macula->kappa4), 0, 1,
				  "kappa4", 0);
  if(check) {*iret = i; return 1;}
  Macula->kappa4_vary = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Macula->kappa4_vary = 1;
      i++;
    }
  }

  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "c1", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Macula->c1), 0, 1,
				  "c1", 0);
  if(check) {*iret = i; return 1;}
  Macula->c1_vary = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Macula->c1_vary = 1;
      i++;
    }
  }

  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "c2", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Macula->c2), 0, 1,
				  "c2", 0);
  if(check) {*iret = i; return 1;}
  Macula->c2_vary = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Macula->c2_vary = 1;
      i++;
    }
  }

  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "c3", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Macula->c3), 0, 1,
				  "c3", 0);
  if(check) {*iret = i; return 1;}
  Macula->c3_vary = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Macula->c3_vary = 1;
      i++;
    }
  }

  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "c4", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Macula->c4), 0, 1,
				  "c4", 0);
  if(check) {*iret = i; return 1;}
  Macula->c4_vary = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Macula->c4_vary = 1;
      i++;
    }
  }

  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "d1", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Macula->d1), 0, 1,
				  "d1", 0);
  if(check) {*iret = i; return 1;}
  Macula->d1_vary = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Macula->d1_vary = 1;
      i++;
    }
  }

  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "d2", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Macula->d2), 0, 1,
				  "d2", 0);
  if(check) {*iret = i; return 1;}
  Macula->d2_vary = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Macula->d2_vary = 1;
      i++;
    }
  }

  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "d3", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Macula->d3), 0, 1,
				  "d3", 0);
  if(check) {*iret = i; return 1;}
  Macula->d3_vary = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Macula->d3_vary = 1;
      i++;
    }
  }

  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "d4", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Macula->d4), 0, 1,
				  "d4", 0);
  if(check) {*iret = i; return 1;}
  Macula->d4_vary = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Macula->d4_vary = 1;
      i++;
    }
  }

  check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "blend", 1,
				  VARTOOLS_TYPE_DOUBLE,
				  (void *) (&Macula->blend), 0, 1,
				  "blend", 0);
  if(check) {*iret = i; return 1;}
  Macula->blend_vary = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Macula->blend_vary = 1;
      i++;
    }
  }

  check = VARTOOLS_ParseConstantParameter(p, c, &i, argv, argc, "Nspot", 
					  VARTOOLS_TYPE_INT, 
					  (void *) (&Macula->Nspot), 0);
  if(check) {*iret = i; return 1;}
  if(Macula->Nspot <= 0) {
    *iret = i; return 1;
  }

  /**** Allocate Memory for the arrays which will store the spot parameters ***/
  if((Macula->lambda0 = (double **) malloc(Macula->Nspot * sizeof(double *))) == NULL ||
     (Macula->phi0 = (double **) malloc(Macula->Nspot * sizeof(double *))) == NULL ||
     (Macula->alphamax = (double **) malloc(Macula->Nspot * sizeof(double *))) == NULL ||
     (Macula->fspot = (double **) malloc(Macula->Nspot * sizeof(double *))) == NULL ||
     (Macula->tmax = (double **) malloc(Macula->Nspot * sizeof(double *))) == NULL ||
     (Macula->life = (double **) malloc(Macula->Nspot * sizeof(double *))) == NULL ||
     (Macula->ingress = (double **) malloc(Macula->Nspot * sizeof(double *))) == NULL ||
     (Macula->egress = (double **) malloc(Macula->Nspot * sizeof(double *))) == NULL ||
     (Macula->lambda0_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
     (Macula->phi0_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
     (Macula->alphamax_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
     (Macula->fspot_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
     (Macula->tmax_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
     (Macula->life_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
     (Macula->ingress_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
     (Macula->egress_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL) {
    VARTOOLS_error(ERR_MEMALLOC);
  }

  /* Parse the spot parameters */
  for(jspot = 0; jspot < Macula->Nspot; jspot++) {

    sprintf(tmpstr,"lambda0_%d", jspot);
    check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "Lambda0", 1,
				    VARTOOLS_TYPE_DOUBLE,
				    (void *) (&(Macula->lambda0[jspot])), 0, 1,
				    tmpstr, 0);
    if(check) {*iret = i; return 1;}
    Macula->lambda0_vary[jspot] = 0;
    if(i < argc) {
      if(!strcmp(argv[i],"vary")) {
	Macula->lambda0_vary[jspot] = 1;
	i++;
      }
    }

    sprintf(tmpstr,"phi0_%d", jspot);
    check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "Phi0", 1,
				    VARTOOLS_TYPE_DOUBLE,
				    (void *) (&(Macula->phi0[jspot])), 0, 1,
				    tmpstr, 0);
    if(check) {*iret = i; return 1;}
    Macula->phi0_vary[jspot] = 0;
    if(i < argc) {
      if(!strcmp(argv[i],"vary")) {
	Macula->phi0_vary[jspot] = 1;
	i++;
      }
    }

    sprintf(tmpstr,"alphamax_%d", jspot);
    check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "alphamax", 1,
				    VARTOOLS_TYPE_DOUBLE,
				    (void *) (&(Macula->alphamax[jspot])), 0, 1,
				    tmpstr, 0);
    if(check) {*iret = i; return 1;}
    Macula->alphamax_vary[jspot] = 0;
    if(i < argc) {
      if(!strcmp(argv[i],"vary")) {
	Macula->alphamax_vary[jspot] = 1;
	i++;
      }
    }

    sprintf(tmpstr,"fspot_%d", jspot);
    check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "fspot", 1,
				    VARTOOLS_TYPE_DOUBLE,
				    (void *) (&(Macula->fspot[jspot])), 0, 1,
				    tmpstr, 0);
    if(check) {*iret = i; return 1;}
    Macula->fspot_vary[jspot] = 0;
    if(i < argc) {
      if(!strcmp(argv[i],"vary")) {
	Macula->fspot_vary[jspot] = 1;
	i++;
      }
    }

    sprintf(tmpstr,"tmax_%d", jspot);
    check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "tmax", 1,
				    VARTOOLS_TYPE_DOUBLE,
				    (void *) (&(Macula->tmax[jspot])), 0, 1,
				    tmpstr, 0);
    if(check) {*iret = i; return 1;}
    Macula->tmax_vary[jspot] = 0;
    if(i < argc) {
      if(!strcmp(argv[i],"vary")) {
	Macula->tmax_vary[jspot] = 1;
	i++;
      }
    }

    sprintf(tmpstr,"life_%d", jspot);
    check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "life", 1,
				    VARTOOLS_TYPE_DOUBLE,
				    (void *) (&(Macula->life[jspot])), 0, 1,
				    tmpstr, 0);
    if(check) {*iret = i; return 1;}
    Macula->life_vary[jspot] = 0;
    if(i < argc) {
      if(!strcmp(argv[i],"vary")) {
	Macula->life_vary[jspot] = 1;
	i++;
      }
    }

    sprintf(tmpstr,"ingress_%d", jspot);
    check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "ingress", 1,
				    VARTOOLS_TYPE_DOUBLE,
				    (void *) (&(Macula->ingress[jspot])), 0, 1,
				    tmpstr, 0);
    if(check) {*iret = i; return 1;}
    Macula->ingress_vary[jspot] = 0;
    if(i < argc) {
      if(!strcmp(argv[i],"vary")) {
	Macula->ingress_vary[jspot] = 1;
	i++;
      }
    }

    sprintf(tmpstr,"egress_%d", jspot);
    check = VARTOOLS_ParseParameter(p, c, &i, argv, argc, "egress", 1,
				    VARTOOLS_TYPE_DOUBLE,
				    (void *) (&(Macula->egress[jspot])), 0, 1,
				    tmpstr, 0);
    if(check) {*iret = i; return 1;}
    Macula->egress_vary[jspot] = 0;
    if(i < argc) {
      if(!strcmp(argv[i],"vary")) {
	Macula->egress_vary[jspot] = 1;
	i++;
      }
    }
    
  }

  Macula->fluxinput = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"fluxinput")) {
      Macula->fluxinput = 1;
      i++;
    }
  }

  Macula->fluxoutput = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"fluxoutput")) {
      Macula->fluxoutput = 1;
      i++;
    }
  }

  Macula->correctlc = 0; 
  if(i < argc) {
    if(!strcmp(argv[i],"correctlc")) {
      Macula->correctlc = 1;
      i++;
    }
  }

  check = VARTOOLS_ParseOutNameKeyword(p, c, &i, argv, argc, "omodel",
				       &(Macula->outputmodel),
				       &(Macula->modeloutdir),
				       &(Macula->outputmodel_useformat),
				       &(Macula->outputmodel_format));
  if(check == 2) {*iret = i; return 1;}
  Macula->outputmodel_tdelv = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"tdelv")) {
      Macula->outputmodel_tdelv = 1;
      i++;
    }
  }
  
  check = VARTOOLS_ParseOutNameKeyword(p, c, &i, argv, argc, "ocurve",
				       &(Macula->outputcurve),
				       &(Macula->curveoutdir),
				       &(Macula->outputcurve_useformat),
				       &(Macula->outputcurve_format));
  if(check == 2) {*iret = i; return 1;}
  Macula->outputcurve_tdelv = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"tdelv")) {
      Macula->outputcurve_tdelv = 1;
      i++;
    }
  }
  check = VARTOOLS_ParseConstantParameter(p, c, &i, argv, argc, "step", 
					  VARTOOLS_TYPE_DOUBLE, 
					  (void *) (&Macula->outputcurve_step), 0);
  if(!check) {
    Macula->outputcurve_stepgiven = 1;
  }
  else if(check == 1) {
    Macula->outputcurve_stepgiven = 0;
  }
  else {
    *iret = i; return 1;
  }
  
  /* If we are fitting the model, register chi2 and Ndof as computed 
     vectors */
  if(Macula->injectorfit == MACULA_FIT) {
    VARTOOLS_RegisterDataVector(p, c, (void *) (&Macula->chi2val),
				VARTOOLS_TYPE_DOUBLE,
				0, VARTOOLS_SOURCE_COMPUTED, 1, "CHI2");
    VARTOOLS_RegisterDataVector(p, c, (void *) (&Macula->Ndof),
				VARTOOLS_TYPE_INT,
				0, VARTOOLS_SOURCE_COMPUTED, 1, "NDOF");
  }
  
  *iret = i;
  
  return 0;
}

void macula_ShowSyntax(FILE *outfile)
/* Write out the expected syntax for this command to the file outfile */
{
  fprintf(outfile,
  "-macula < \"inject\" | \"fit\" <\"amoeba\" | \"lm\"> >\n"
  "     <\"Prot\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "         [\"vary\"]>\n"
  "     <\"istar\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "         [\"vary\"]>\n"
  "     <\"kappa2\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "         [\"vary\"]>\n"
  "     <\"kappa4\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "         [\"vary\"]>\n"
  "     <\"c1\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "         [\"vary\"]>\n"
  "     <\"c2\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "         [\"vary\"]>\n"
  "     <\"c3\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "         [\"vary\"]>\n"
  "     <\"c4\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "         [\"vary\"]>\n"
  "     <\"d1\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "         [\"vary\"]>\n"
  "     <\"d2\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "         [\"vary\"]>\n"
  "     <\"d3\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "         [\"vary\"]>\n"
  "     <\"d4\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "         [\"vary\"]>\n"
  "     <\"blend\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "         [\"vary\"]>\n"
  "     <\"Nspot\" value>\n"
  "        ... Repeat below Nspot times ....\n"
  "        <\"Lambda0\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "           [\"vary\"]>\n"
  "        <\"Phi0\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "           [\"vary\"]>\n"
  "        <\"alphamax\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "           [\"vary\"]>\n"
  "        <\"fspot\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "           [\"vary\"]>\n"
  "        <\"tmax\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "           [\"vary\"]>\n"
  "        <\"life\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "           [\"vary\"]>\n"
  "        <\"ingress\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "           [\"vary\"]>\n"
  "        <\"egress\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
  "           [\"vary\"]>\n"
  "        .... Stop repeating here ....\n"
  "    [\"fluxinput\"] [\"fluxoutput\"] [\"correctlc\"]\n"
  "    [\"omodel\" <outdir [\"nameformat\" fmt]> [\"tdelv\"]]\n"
  "    [\"ocurve\" <\"outdir\" outdir [ \"nameformat\" fmt]>\n"
  "        [\"tdelv\"] [\"step\" stepsize]]\n"
	  );
}

void macula_ShowHelp(FILE *outfile)
/* Output the help information for this command */
{
  fprintf(outfile,
"Fit or inject a Macula spot model [in]to the light curves. The reference for macula is \"Kipping 2012, arXiv:1209.2985\", please be sure to cite this paper if you make use of Macula. The user must give either the \"inject\" keyword or the \"fit\" keyword to indicate whether the model should be injected into the light curve, or fitted to it. If fitting, give the \"amoeba\" keyword to do downhill simplex fitting, or the \"lm\" to use the Levenberg-Marquardt fitting algorithm. After that the user specifies how the initialize the parameters used by this model. The parameters include:\n"
"\t\"Prot\" - the equatorial rotation period [input lc time units].\n"
"\t\"istar\" - the inclination of the star [radians].\n"
"\t\"kappa2\" - quadratic differential rotation coeff.\n"
"\t\"kappa4\" - quartic differential rotation coeff.\n"
"\t\"c1\" - 1st of four-coeff stellar LD terms.\n"
"\t\"c2\" - 2nd of four-coeff stellar LD terms.\n"
"\t\"c3\" - 3rd of four-coeff stellar LD terms.\n"
"\t\"c4\" - 4th of four-coeff stellar LD terms.\n"
"\t\"d1\" - 1st of four-coeff spot LD terms.\n"
"\t\"d2\" - 2nd of four-coeff spot LD terms.\n"
"\t\"d3\" - 3rd of four-coeff spot LD terms.\n"
"\t\"d4\" - 4th of four-coeff spot LD terms.\n"
"\t\"blend\" - blend parameter.\n"
"You must specify the number of spots to use, with the parameter NSpot.\n"
"For each spot, the following parameters are required:\n"
"\t\"Lambda0\" - Longitude of spot at its time of maximum size [radians].\n"
"\t\"Phi0\" - Latitude of spot at its time of maximum size [radians].\n"
"\t\"alphamax\" - Maximum angular size of spot [radians].\n"
"\t\"fspot\" - Spot-to-star flux contrast.\n"
"\t\"tmax\" - Reference time of maximum spot size [input lc time units].\n"
"\t\"life\" - Lifetime of spot (FWHM) [input lc time units].\n"
"\t\"ingress\" - ingress duration of spot [input lc time units].\n"
"\t\"egress\" - egress duration of spot [input lc time units].\n"
"For each parameter the user specifies the source for the initial value. This can either be the keyword \"fix\" followed by the value to use for all light curves, \"list\" to read the parameter from the light curve list (use the \"column\" keyword followed by a number to indicate the column number from the list to use, otherwise the next column in the list will be assumed), \"fixcolumn\" followed by the name or number of an output column from a previously executed command, or \"expr\" followed by an analytic expression. For each parameter the user can give the \"vary\" keyword which indicated that the parameter is to be varied in a fit. If not given, the parameter will be fixed in the fit. Following the parameters are number of optional keywords that control the behavior of the command. If the user gives the \"fluxinput\" command, then the light curve input into the command will be assumed to be in flux units, by default magnitudes are assumed. If the user gives the \"fluxoutput\" command, the output light curve and models will be in flux units, by default they are in magnitudes. The user may use the \"correctlc\" keyword to subtract the best-fit model from the light curve, the \"omodel\" keyword to output the model for each light curve, and the \"ocurve\" keyword to output a model curve file sampled at uniform spacing. If the \"omodel\" keyword is given, the user should specify the directory to output the model light curves to. By default the output light curve will have the name outdir/BASELC_NAME.macula where BASELC_Name is the basename of the input light curve. The user may, however, give the \"format\" keyword followed by a format string to specify arbitrary filenames. The syntax is the same as for the \"-o\" VARTOOLS command. If the \"tdelv\" keyword is given, then predicted transit depth variations will be included in the output model.  For \"ocurve\" the default output name is outdir/BASELC_NAME.maculacurve. In addition to the options available for \"omodel\", the user may also specify the time step size to use in generating the curve.\nAgain, if you use this routine, cite:\n\tKipping 2012, arXiv:1209.2985.\n\n"
	  );
}

void macula_ShowExample(FILE *outfile)
{
  fprintf(outfile,
          "Example 1:\n\n"
	  "./vartools -i EXAMPLES/3 -L USERLIB/src/macula.so \\\n"
	  "   -stats t min \\\n"
	  "   -expr 'mag=10.0+err*gauss()' \\\n"
	  "   -macula inject \\\n"
	  "      Prot fix 1.234567 \\\n"
	  "      istar fix 1.4567 \\\n"
	  "      kappa2 fix 0.0 \\\n"
	  "      kappa4 fix 0.0 \\\n"
	  "      c1 fix 0.2 c2 fix 0.1 c3 fix 0.0 c4 fix 0.0 \\\n"
	  "      d1 fix 0.2 d2 fix 0.1 d3 fix 0.0 d4 fix 0.0 \\\n"
	  "      blend fix 1.0 \\\n"
	  "      Nspot 1 \\\n"
	  "         Lambda0 fix 0.0 \\\n"
	  "         Phi0 fix 1.2345 \\\n"
          "         alphamax fix 0.2 \\\n"
          "         fspot fix 0.1 \\\n"
          "         tmax fixcolumn STATS_t_MIN_0 \\\n"
          "         life fix 1000.0 \\\n"
          "         ingress fix 0.1 \\\n"
          "         egress fix 0.1 \\\n"
          "   -o EXAMPLES/OUTDIR1/3.maculainject \\\n"
          "   -oneline\n\n"
	  );
  fprintf(outfile,
	  "Simulates the light curve of a spotted star using macula. We read-in the light curve EXAMPLES/3 for its timebase and magnitude uncertainties. The call to -stats is used to get the minimum time in the light curve, which we will use later in the call to -macula. The call to -expr replaces the magnitudes in the light curve with simple gaussian random noise. The call to -macula then injects the spot with all of the parameters fixed to values given on the command-line, save for \"tmax\" which we set to a previously determined quantity (STATS_t_MIN_0, which is the minimum time in the light curve). The simulated light curve is then output to the file EXAMPLES/OUTDIR1/3.maculainject.\n\n");
  fprintf(outfile,
          "Example 2:\n\n"
	  "./vartools -i EXAMPLES/OUTDIR1/3.maculainject -L USERLIB/src/macula.so \\\n"
	  "   -stats t min \\\n"
	  "   -LS 0.1 100 0.1 1 0 \\\n"
	  "   -macula fit amoeba \\\n"
	  "      Prot fixcolumn LS_Period_1_1 vary \\\n"
	  "      istar fix 1.6 vary \\\n"
	  "      kappa2 fix 0.0 \\\n"
	  "      kappa4 fix 0.0 \\\n"
	  "      c1 fix 0.2 c2 fix 0.1 c3 fix 0.0 c4 fix 0.0 \\\n"
	  "      d1 fix 0.2 d2 fix 0.1 d3 fix 0.0 d4 fix 0.0 \\\n"
	  "      blend fix 1.0 \\\n"
	  "      Nspot 1 \\\n"
	  "         Lambda0 fix 0.0 \\\n"
	  "         Phi0 fix 1.2345 \\\n"
          "         alphamax fix 0.2 \\\n"
          "         fspot fix 0.1 \\\n"
          "         tmax fixcolumn STATS_t_MIN_0 \\\n"
          "         life fix 1000.0 \\\n"
          "         ingress fix 0.1 \\\n"
          "         egress fix 0.1 \\\n"
          "      omodel EXAMPLES/OUTDIR1 \\\n"
          "      ocurve EXAMPLES/OUTDIR1 \\\n"
          "   -oneline\n\n"
	  );
  fprintf(outfile,
	  "Fit a macula starspot model to the light curve that was simulated in example 1. In this case we use the -LS command to first find the period of the light curve. When then run macula, using \"fit amoeba\" to run a fit with the downhill-simplex algorithm. We include the \"vary\" keyword for the \"Prot\" and \"istar\" parameters, all other parameters will be held constant in the fit. The initial guesses of the parameters are given, for \"Prot\" we use the period found by -LS. A model evaluated at the observed times is output using the \"omodel\" keyword, while a model evaluated at uniformly sampled times (for plotting) is output using the \"ocurve\" keyword.\n\n");

}

void macula_RunCommand(ProgramData *p, void *userdata, int lc_name_num, int lc_num)
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
  int NJD;
  double *t, *mag, *err;
  char *lcname, *format;
  _Macula *Macula;

  char lcoutname[MAXLEN];
  char curvename[MAXLEN];
  
  void domaculainject(int NJD, double *t, double *mag, double *err, 
		    _Macula *Macula, int lc_num, char *lcoutname, 
		      char *curvename);

  void domaculafit(int NJD, double *t, double *mag, double *err, 
		   _Macula *Macula, int lc_num, char *lcoutname, 
		   char *curvename);

  /* Cast the userdata pointer given by VARTOOLS into a _Macula type
     struct. Then we can access its contents. */
  Macula = (_Macula *)userdata;

  /* The number of points in the light curve, and vectors storing the
     time, magnitudes, and errors are accessed as illustrated
     below. Note that t, mag, and err are vectors, t[0] is the time of
     the first observation in the light curve, t[NJD-1] is the time of
     the last observation, etc. */
  NJD=p->NJD[lc_num];
  t = p->t[lc_num];
  mag = p->mag[lc_num];
  err = p->sig[lc_num];

  lcname = p->lcnames[lc_name_num];

    /* This command allows the user to optionally output model light curves.
     It uses a "standard" VARTOOLS syntax for specifying the output name,
     which is to provide a directory for the output, and then by default
     set the name to something like outdir/BASELC_NAME.suffix unless the
     user gives a format keyword, in which case the user can provide an
     arbitrary name. Any command using this convention can call the
     VARTOOLS_GetOutputFilename function to get the output filename. */
  if(Macula->outputmodel) {
    if(Macula->outputmodel_useformat)
      format = Macula->outputmodel_format;
    else
      format = NULL;
    VARTOOLS_GetOutputFilename(lcoutname,lcname,Macula->modeloutdir,"macula",
			       format, lc_name_num);
  }
  else 
    lcoutname[0] = '\0';

  if(Macula->outputcurve) {
    if(Macula->outputcurve_useformat)
      format = Macula->outputcurve_format;
    else
      format = NULL;
    VARTOOLS_GetOutputFilename(curvename,lcname,Macula->curveoutdir,"maculacurve",
			       format, lc_name_num);
  }
  else 
    curvename[0] = '\0';

  if(Macula->injectorfit == MACULA_INJECT) {
    domaculainject(NJD, t, mag, err, Macula, lc_num, lcoutname, curvename);
  } else {
    domaculafit(NJD, t, mag, err, Macula, lc_num, lcoutname, 
		curvename);
  }

}

/***************************************************************************
 Below are functions that are specific to this command. These are not
 required by a user library, but in this case are called through the
 macula_RunCommand function and perform the actual processing routine.
***************************************************************************/


void domaculainject(int NJD, double *t, double *mag, double *err, 
		    _Macula *Macula, int lc_num, char *lcoutname, 
		    char *curvename)
{
  /* Inject a Macula model into a light curve. Optionally output the
     model.  If the user requested that the model be subtracted with
     the correctlc keyword we will do it, but it is not obvious why
     they would want this. */
  
  double *fluxsim = NULL, *magsim = NULL,
    *deltaratio = NULL, *tsimcurve = NULL;

  double *Theta_star = NULL, **Theta_spot = NULL, *dFmoddt = NULL;

  FILE *outfile;
  
  double flux, magout, tstep, ttmp;

  double baseflux, blend;

  int i, j, k, testout, Nstep;

  int lengthsimvec = 0;

  void simulatemaculalc(int N, double *JD, double *fluxout, double *Theta_star, int Nspot, double **Theta_spot, double baseflux, double blend, int derivatives, int temporal, int TdeltaV, double *deltaratio, double **dFmod_star, double ***dFmod_spot, double ***dFmod_inst, double *dFmoddt);
  

  /* Create temporary arrays */
  if((Theta_star = (double *) malloc(MACULA_NDIM_STAR_PARAM * sizeof(double))) == NULL ||
     (Theta_spot = (double **) malloc(MACULA_NDIM_SPOT_PARAM * sizeof(double *))) == NULL ||
     (fluxsim = (double *) malloc(NJD * sizeof(double))) == NULL ||
     (dFmoddt = (double *) malloc(NJD * sizeof(double))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);

  if(Macula->Nspot > 0) {
    for(j=0; j < MACULA_NDIM_SPOT_PARAM; j++) {
      if((Theta_spot[j] = (double *) malloc(Macula->Nspot * sizeof(double))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
    }
  }
  
  /* Fill in parameters */
  Theta_star[0] = Macula->Istar[lc_num];
  Theta_star[1] = Macula->Prot[lc_num];
  Theta_star[2] = Macula->kappa2[lc_num];
  Theta_star[3] = Macula->kappa4[lc_num];
  Theta_star[4] = Macula->c1[lc_num];
  Theta_star[5] = Macula->c2[lc_num];
  Theta_star[6] = Macula->c3[lc_num];
  Theta_star[7] = Macula->c4[lc_num];
  Theta_star[8] = Macula->d1[lc_num];
  Theta_star[9] = Macula->d2[lc_num];
  Theta_star[10] = Macula->d3[lc_num];
  Theta_star[11] = Macula->d4[lc_num];

  for(k=0; k < Macula->Nspot; k++) {
    Theta_spot[0][k] = Macula->lambda0[k][lc_num];
    Theta_spot[1][k] = Macula->phi0[k][lc_num];
    Theta_spot[2][k] = Macula->alphamax[k][lc_num];
    Theta_spot[3][k] = Macula->fspot[k][lc_num];
    Theta_spot[4][k] = Macula->tmax[k][lc_num];
    Theta_spot[5][k] = Macula->life[k][lc_num];
    Theta_spot[6][k] = Macula->ingress[k][lc_num];
    Theta_spot[7][k] = Macula->egress[k][lc_num];
  }

  lengthsimvec = NJD;
  if(NJD > 0) {
    if((fluxsim = (double *) malloc(NJD * sizeof(double))) == NULL ||
       (magsim = (double *) malloc(NJD * sizeof(double))) == NULL ||
       (deltaratio = (double *) malloc(NJD * sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
  }

  baseflux = 1.0;
  blend = Macula->blend[lc_num];

  simulatemaculalc(NJD, t, fluxsim, Theta_star, Macula->Nspot, Theta_spot,
		   baseflux, blend, 0, 0, 1, deltaratio, NULL, NULL, NULL, 
		   dFmoddt);

  testout = 0;
  if(lcoutname != NULL ? (lcoutname[0] != '\0' ? 1 : 0) : 0) {
    if((outfile = fopen(lcoutname,"w")) == NULL) {
      VARTOOLS_error2(ERR_CANNOTWRITE,lcoutname);
    }
    testout = 1;
  }

  if(!Macula->fluxoutput && !Macula->fluxinput) {
    for(i=0; i < NJD; i++) {
      if(!isnan(mag[i])) {
	magsim[i] = -2.5*log10(fluxsim[i]);
	mag[i] += magsim[i];
	if(testout) {
	  fprintf(outfile,"%.17g %.17g %.17g %.17g", t[i], mag[i], err[i],
		  magsim[i]);
	  if(Macula->outputmodel_tdelv) {
	    fprintf(outfile," %.17g", deltaratio[i]);
	  }
	  fprintf(outfile,"\n");
	}
	if(Macula->correctlc) {
	  mag[i] -= magsim[i];
	}
      }
    }
  }
  else if(!Macula->fluxoutput && Macula->fluxinput) {
    for(i=0; i < NJD; i++) {
      if(!isnan(mag[i])) {
	magsim[i] = -2.5*log10(fluxsim[i]);
	err[i] = 1.0857*err[i]/mag[i];
	magout = 25.0 - 2.5*log10(mag[i]);
	magout += magsim[i];
	mag[i] = magout;
	if(testout) {
	  fprintf(outfile,"%.17g %.17g %.17g %.17g", t[i], mag[i], err[i], magsim[i]);
	  if(Macula->outputmodel_tdelv) {
	    fprintf(outfile," %.17g", deltaratio[i]);
	  }
	  fprintf(outfile,"\n");
	}
	if(Macula->correctlc) {
	  mag[i] -= magsim[i];
	}
      }
    }
  }
  else if(Macula->fluxoutput && !Macula->fluxinput) {
    for(i=0; i < NJD; i++) {
      if(!isnan(mag[i])) {
	flux = pow(10.0,(-0.4*(mag[i] - 25.0)));
	err[i] = err[i]*flux/1.0857;
	flux = flux*fluxsim[i];
	mag[i] = flux;
	if(testout) {
	  fprintf(outfile,"%.17g %.17g %.17g %.17g", t[i], mag[i], err[i], fluxsim[i]);
	  if(Macula->outputmodel_tdelv) {
	    fprintf(outfile," %.17g", deltaratio[i]);
	  }
	  fprintf(outfile,"\n");
	}
	if(Macula->correctlc) {
	  mag[i] = mag[i]/fluxsim[i];
	}
      }
    }
  }
  else if(Macula->fluxoutput && Macula->fluxinput) {
    for(i=0; i < NJD; i++) {
      if(!isnan(mag[i])) {
	mag[i] = mag[i]*fluxsim[i];
	if(testout) {
	  fprintf(outfile,"%.17g %.17g %.17g %.17g", t[i], mag[i], err[i], fluxsim[i]);
	  if(Macula->outputmodel_tdelv) {
	    fprintf(outfile," %.17g", deltaratio[i]);
	  }
	  fprintf(outfile,"\n");
	}
	if(Macula->correctlc) {
	  mag[i] = mag[i]/fluxsim[i];
	}
      }
    }
  }
  if(testout)
    fclose(outfile);

  /* write out the curve file if requested */
  if(Macula->outputcurve) {
    
    if((outfile = fopen(curvename,"w")) == NULL) {
      VARTOOLS_error2(ERR_CANNOTWRITE,curvename);
    }
    
    if(Macula->outputcurve_stepgiven) {
      tstep = Macula->outputcurve_step;
    }
    else {
      /* Set the step size equal to the minimum non-zero step between
	 points in the input t vector */
      tstep = t[NJD-1] - t[0];
      for(i=1; i < NJD; i++) {
	if(t[i] - t[i-1] > 0 && t[i] - t[i-1] < tstep) {
	  tstep = t[i] - t[i-1];
	}
      }
    }
    Nstep = ceil((t[NJD-1] - t[0])/tstep) + 1;
    if(Nstep > lengthsimvec) {
      if((fluxsim = (double *) realloc(fluxsim, Nstep*sizeof(double))) == NULL ||
	 (magsim = (double *) realloc(magsim, Nstep*sizeof(double))) == NULL ||
	 (deltaratio = (double *) realloc(deltaratio, Nstep*sizeof(double))) == NULL ||
	 (dFmoddt = (double *) realloc(dFmoddt, Nstep*sizeof(double))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
    }
    if((tsimcurve = (double *) malloc(Nstep*sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    Nstep = 0;
    ttmp = t[0];
    while(ttmp <= t[NJD-1]) {
      tsimcurve[Nstep] = ttmp;
      ttmp += tstep;
      Nstep++;
    }
    simulatemaculalc(Nstep, tsimcurve, fluxsim, Theta_star, Macula->Nspot, 
		     Theta_spot, baseflux, blend, 0, 0, 1, deltaratio, 
		     NULL, NULL, NULL, dFmoddt);
    
    if(!Macula->fluxoutput) {
      for(i=0; i < Nstep; i++) {
	magsim[i] = -2.5*log10(fluxsim[i]);
	fprintf(outfile,"%.17g %.17g", tsimcurve[i], magsim[i]);
	if(Macula->outputcurve_tdelv)
	  fprintf(outfile," %.17g", deltaratio[i]);
	fprintf(outfile,"\n");
      }
    } else {
      for(i=0; i < Nstep; i++) {
	fprintf(outfile,"%.17g %.17g", tsimcurve[i], fluxsim[i]);
	if(Macula->outputcurve_tdelv)
	  fprintf(outfile," %.17g", deltaratio[i]);
	fprintf(outfile,"\n");
      }
    }      
    fclose(outfile);
  }

  if(fluxsim != NULL)
    free(fluxsim);
  if(magsim != NULL)
    free(magsim);
  if(deltaratio != NULL)
    free(deltaratio);
  if(tsimcurve != NULL)
    free(tsimcurve);
  if(Theta_star != NULL)
    free(Theta_star);
  if(Theta_spot != NULL) {
    for(k=0; k < MACULA_NDIM_SPOT_PARAM; k++) {
      if(Theta_spot[k] != NULL)
	free(Theta_spot[k]);
    }
    free(Theta_spot);
  }
  if(dFmoddt != NULL)
    free(dFmoddt);
      
}

void domaculafit(int NJD, double *t, double *mag, double *err, 
		 _Macula *Macula, int lc_num, char *lcoutname, 
		 char *curvename)
{
  void domaculafit_lm(int NJD, double *t, double *mag, double *err, 
		      _Macula *Macula, int lc_num, char *lcoutname, 
		      char *curvename);
  void domaculafit_amoeba(int NJD, double *t, double *mag, double *err, 
		      _Macula *Macula, int lc_num, char *lcoutname, 
		      char *curvename);
  if(Macula->fittype == MACULA_FITTYPE_LM) {
    domaculafit_lm(NJD, t, mag, err, Macula, lc_num, lcoutname,
		   curvename);
  }
  else if(Macula->fittype == MACULA_FITTYPE_AMOEBA) {
    domaculafit_amoeba(NJD, t, mag, err, Macula, lc_num, lcoutname,
		       curvename);
  }
  return;
}

#define DELTACHI_CONVERGENCELIMIT 1e-5

void domaculafit_lm(int NJD, double *t, double *mag, double *err, 
		 _Macula *Macula, int lc_num, char *lcoutname, 
		    char *curvename)
{

  int *ia = NULL, ma, ma_max, i, j, k;
  double outmag;
  double *outflux = NULL;
  double *influx = NULL, *influx_err = NULL;
  double *a = NULL, **covar = NULL, **alpha = NULL, chisq, alamda = -1., ochisq,    oldchisq,
    *atry = NULL, *beta = NULL, *da = NULL, **oneda = NULL;

  _MaculaLMFitStruct MaculaFit;

  FILE *outfile;

  void macula_evalfunc_lm(double *t, double *a, double *yfit, double **dyda, int ma, int N, int Nlin_coeffs, double *lin_coeffs, double **Design_Matrix, double *y, double *sig, int *varylin_coeffs, void *userparams);

  void simulatemaculalc(int N, double *JD, double *fluxout, double *Theta_star, int Nspot, double **Theta_spot, double baseflux, double blend, int derivatives, int temporal, int TdeltaV, double *deltaratio, double **dFmod_star, double ***dFmod_spot, double ***dFmod_inst, double *dFmoddt);

  int iter = 0, testout, Nstep;
  double oalamda, deltachi, tstep, ttmp, magout;

  double *fluxsim = NULL, *deltaratio_sim = NULL, *dFmoddt_sim = NULL,
    *tsimcurve = NULL;

  ma_max = MACULA_NDIM_STAR_PARAM + MACULA_NDIM_SPOT_PARAM*Macula->Nspot 
    + MACULA_NDIM_INST_PARAM;

  if((a = (double *) malloc(ma_max * sizeof(double))) == NULL) {
    VARTOOLS_error(ERR_MEMALLOC);
  }

  /* Setup the MaculaFit structure */
  MaculaFit.Nspot = Macula->Nspot;
  if((MaculaFit.Theta_star = (double *) malloc(MACULA_NDIM_STAR_PARAM * sizeof(double))) == NULL ||
     (MaculaFit.Theta_spot = (double **) malloc(MACULA_NDIM_SPOT_PARAM * sizeof(double *))) == NULL ||
     (MaculaFit.dFmoddt = (double *) malloc(NJD * sizeof(double))) == NULL ||
     (MaculaFit.deltaratio = (double *) malloc(NJD * sizeof(double))) == NULL ||
     (MaculaFit.dFmod_star = (double **) malloc(NJD * sizeof(double *))) == NULL ||
     (MaculaFit.dFmod_spot = (double ***) malloc(NJD * sizeof(double **))) == NULL ||
     (MaculaFit.dFmod_inst = (double ***) malloc(NJD * sizeof(double **))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);

  if(Macula->Nspot > 0) {
    if((MaculaFit.lambda0_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
       (MaculaFit.phi0_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
       (MaculaFit.alphamax_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
       (MaculaFit.fspot_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
       (MaculaFit.tmax_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
       (MaculaFit.life_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
       (MaculaFit.ingress_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
       (MaculaFit.egress_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    for(j=0; j < MACULA_NDIM_SPOT_PARAM; j++) {
      if((MaculaFit.Theta_spot[j] = (double *) malloc(Macula->Nspot * sizeof(double))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
    }
  }
  for(i=0; i < NJD; i++) {
    if((MaculaFit.dFmod_star[i] = (double *) malloc(MACULA_NDIM_STAR_PARAM * sizeof(double))) == NULL ||
       (MaculaFit.dFmod_spot[i] = (double **) malloc(MACULA_NDIM_SPOT_PARAM * sizeof(double *))) == NULL ||
       (MaculaFit.dFmod_inst[i] = (double **) malloc(MACULA_NDIM_INST_PARAM * sizeof(double *))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    for(j=0; j < MACULA_NDIM_SPOT_PARAM; j++) {
      if((MaculaFit.dFmod_spot[i][j] = (double *) malloc(Macula->Nspot * sizeof(double))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
    }
    for(j=0; j < MACULA_NDIM_INST_PARAM; j++) {
      if((MaculaFit.dFmod_inst[i][j] = (double *) malloc(sizeof(double))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
    }
  }

  /* Fill in parameters */
  MaculaFit.Theta_star[0] = Macula->Istar[lc_num];
  MaculaFit.Theta_star[1] = Macula->Prot[lc_num];
  MaculaFit.Theta_star[2] = Macula->kappa2[lc_num];
  MaculaFit.Theta_star[3] = Macula->kappa4[lc_num];
  MaculaFit.Theta_star[4] = Macula->c1[lc_num];
  MaculaFit.Theta_star[5] = Macula->c2[lc_num];
  MaculaFit.Theta_star[6] = Macula->c3[lc_num];
  MaculaFit.Theta_star[7] = Macula->c4[lc_num];
  MaculaFit.Theta_star[8] = Macula->d1[lc_num];
  MaculaFit.Theta_star[9] = Macula->d2[lc_num];
  MaculaFit.Theta_star[10] = Macula->d3[lc_num];
  MaculaFit.Theta_star[11] = Macula->d4[lc_num];

  for(k=0; k < Macula->Nspot; k++) {
    MaculaFit.Theta_spot[0][k] = Macula->lambda0[k][lc_num];
    MaculaFit.Theta_spot[1][k] = Macula->phi0[k][lc_num];
    MaculaFit.Theta_spot[2][k] = Macula->alphamax[k][lc_num];
    MaculaFit.Theta_spot[3][k] = Macula->fspot[k][lc_num];
    MaculaFit.Theta_spot[4][k] = Macula->tmax[k][lc_num];
    MaculaFit.Theta_spot[5][k] = Macula->life[k][lc_num];
    MaculaFit.Theta_spot[6][k] = Macula->ingress[k][lc_num];
    MaculaFit.Theta_spot[7][k] = Macula->egress[k][lc_num];
  }

  MaculaFit.blend = Macula->blend[lc_num];

  /* Determine the number of free parameters */
  ma = 0;
  MaculaFit.Istar_vary = Macula->Istar_vary;
  if(Macula->Istar_vary) {
    a[ma] = Macula->Istar[lc_num];
    ma++;
  }
  MaculaFit.Prot_vary = Macula->Prot_vary;
  if(Macula->Prot_vary) {
    a[ma] = Macula->Prot[lc_num];
    ma++;
  }
  MaculaFit.kappa2_vary = Macula->kappa2_vary;
  if(Macula->kappa2_vary) {
    a[ma] = Macula->kappa2[lc_num];
    ma++;
  }
  MaculaFit.kappa4_vary = Macula->kappa4_vary;
  if(Macula->kappa4_vary) {
    a[ma] = Macula->kappa4[lc_num];
    ma++;
  }
  MaculaFit.c1_vary = Macula->c1_vary;
  if(Macula->c1_vary) {
    a[ma] = Macula->c1[lc_num];
    ma++;
  }
  MaculaFit.c2_vary = Macula->c2_vary;
  if(Macula->c2_vary) {
    a[ma] = Macula->c2[lc_num];
    ma++;
  }
  MaculaFit.c3_vary = Macula->c3_vary;
  if(Macula->c3_vary) {
    a[ma] = Macula->c3[lc_num];
    ma++;
  }
  MaculaFit.c4_vary = Macula->c4_vary;
  if(Macula->c4_vary) {
    a[ma] = Macula->c4[lc_num];
    ma++;
  }
  MaculaFit.d1_vary = Macula->d1_vary;
  if(Macula->d1_vary) {
    a[ma] = Macula->d1[lc_num];
    ma++;
  }
  MaculaFit.d2_vary = Macula->d2_vary;
  if(Macula->d2_vary) {
    a[ma] = Macula->d2[lc_num];
    ma++;
  }
  MaculaFit.d3_vary = Macula->d3_vary;
  if(Macula->d3_vary) {
    a[ma] = Macula->d3[lc_num];
    ma++;
  }
  MaculaFit.d4_vary = Macula->d4_vary;
  if(Macula->d4_vary) {
    a[ma] = Macula->d4[lc_num];
    ma++;
  }

  for(k=0; k < Macula->Nspot; k++) {
    MaculaFit.lambda0_vary[k] = Macula->lambda0_vary[k];
    if(Macula->lambda0_vary[k]) {
      a[ma] = Macula->lambda0[k][lc_num];
      ma++;
    }
    MaculaFit.phi0_vary[k] = Macula->phi0_vary[k];
    if(Macula->phi0_vary[k]) {
      a[ma] = Macula->phi0[k][lc_num];
      ma++;
    }
    MaculaFit.alphamax_vary[k] = Macula->alphamax_vary[k];
    if(Macula->alphamax_vary[k]) {
      a[ma] = Macula->alphamax[k][lc_num];
      ma++;
    }
    MaculaFit.fspot_vary[k] = Macula->fspot_vary[k];
    if(Macula->fspot_vary[k]) {
      a[ma] = Macula->fspot[k][lc_num];
      ma++;
    }
    MaculaFit.tmax_vary[k] = Macula->tmax_vary[k];
    if(Macula->tmax_vary[k]) {
      a[ma] = Macula->fspot[k][lc_num];
      ma++;
    }
    MaculaFit.life_vary[k] = Macula->life_vary[k];
    if(Macula->life_vary[k]) {
      a[ma] = Macula->life[k][lc_num];
      ma++;
    }
    MaculaFit.ingress_vary[k] = Macula->ingress_vary[k];
    if(Macula->ingress_vary[k]) {
      a[ma] = Macula->ingress[k][lc_num];
      ma++;
    }
    MaculaFit.egress_vary[k] = Macula->egress_vary[k];
    if(Macula->egress_vary[k]) {
      a[ma] = Macula->egress[k][lc_num];
      ma++;
    }
  }

  MaculaFit.blend_vary = Macula->blend_vary;
  if(Macula->blend_vary) {
    a[ma] = Macula->blend[lc_num];
    ma++;
  }

  /* Convert mag to flux if needed */
  if(!Macula->fluxinput) {
    if((influx = (double *) malloc(NJD * sizeof(double))) == NULL ||
       (influx_err = (double *) malloc(NJD * sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);

    for(i=0; i < NJD; i++) {
      influx[i] = pow(10.0,(-0.4*(mag[i] - 25.0)));
      influx_err[i] = err[i]*influx[i]/1.0857;
    }
  } else {
    influx = mag;
    influx_err = err;
  }

  /* Set the base flux level to the 95 percentile flux value */
  a[ma] = VARTOOLS_percentile(NJD, influx, 95.0);

  ma++;

  /* stopped working here... add call to mrqmin, see transitmodel.c
     for an example */

  /* Initialize other vectors used by mrqmin */
  if((ia = (int *) malloc(ma * sizeof(int))) == NULL ||
     (covar = (double **) malloc(ma * sizeof(double *))) == NULL ||
     (alpha = (double **) malloc(ma * sizeof(double *))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  for(i=0; i < ma; i++) {
    if((covar[i] = (double *) malloc(ma * sizeof(double))) == NULL ||
       (alpha[i] = (double *) malloc(ma * sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    ia[i] = 1;
  }
#ifdef PARALLEL
  if((atry = (double *) malloc(ma * sizeof(double))) == NULL ||
     (beta = (double *) malloc(ma * sizeof(double))) == NULL ||
     (da = (double *) malloc(ma * sizeof(double))) == NULL ||
     (oneda = (double **) malloc(ma * sizeof(double *))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  for(k=0; k < ma; k++)
    if((oneda[k] = (double *) malloc(sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);

  VARTOOLS_mrqmin(t, influx, influx_err, NJD, a, ia, ma, covar, alpha, &chisq, &alamda,
	 &macula_evalfunc_lm, 0, NULL, NULL, NULL, ma, &ochisq, atry, beta, 
	 da, oneda, (void *) &MaculaFit);
#else
  VARTOOLS_mrqmin(t, influx, influx_err, NJD, a, ia, ma, covar, alpha, &chisq, &alamda,
	 &macula_evalfunc_lm, 0, NULL, NULL, NULL, (void *) &MaculaFit);
#endif

  oldchisq = chisq;
  ochisq = chisq;
  iter = 0;

  do 
    {
      oalamda = alamda;
#ifdef PARALLEL
      VARTOOLS_mrqmin(t, influx, influx_err, NJD, a, ia, ma, covar, alpha, &chisq, 
	     &alamda, &macula_evalfunc_lm, 
	     0, NULL, NULL, NULL, ma, &ochisq, atry, beta, 
	     da, oneda, (void *) &MaculaFit);
#else
      VARTOOLS_mrqmin(t, influx, influx_err, NJD, a, ia, ma, covar, alpha, &chisq, 
	     &alamda, &macula_evalfunc_lm, 
	     0, NULL, NULL, NULL, (void *) &MaculaFit);
#endif
      deltachi = (oldchisq - chisq);
      oldchisq = chisq;
      ochisq = chisq;
      iter++;
    }
  while((deltachi < 0. || deltachi > DELTACHI_CONVERGENCELIMIT ||
	 iter < 10 || alamda > oalamda) && alamda != 0.0);
  
  Macula->chi2val[lc_num] = chisq;
  Macula->Ndof[lc_num] = NJD - ma;

  /* Update the parameters */

  i = 0;
  if(Macula->Istar_vary) {
    Macula->Istar[lc_num] = a[i];
    i++;
  }
  if(Macula->Prot_vary) {
    Macula->Prot[lc_num] = a[i];
    i++;
  }
  if(Macula->kappa2_vary) {
    Macula->kappa2[lc_num] = a[i];
    i++;
  }
  if(Macula->kappa4_vary) {
    Macula->kappa4[lc_num] = a[i];
    i++;
  }
  if(Macula->c1_vary) {
    Macula->c1[lc_num] = a[i];
    i++;
  }
  if(Macula->c2_vary) {
    Macula->c2[lc_num] = a[i];
    i++;
  }
  if(Macula->c3_vary) {
    Macula->c3[lc_num] = a[i];
    i++;
  }
  if(Macula->c4_vary) {
    Macula->c4[lc_num] = a[i];
    i++;
  }
  if(Macula->d1_vary) {
    Macula->d1[lc_num] = a[i];
    i++;
  }
  if(Macula->d2_vary) {
    Macula->d2[lc_num] = a[i];
    i++;
  }
  if(Macula->d3_vary) {
    Macula->d3[lc_num] = a[i];
    i++;
  }
  if(Macula->d4_vary) {
    Macula->d4[lc_num] = a[i];
    i++;
  }

  for(k=0; k < Macula->Nspot; k++) {
    if(Macula->lambda0_vary[k]) {
      Macula->lambda0[k][lc_num] = a[i];
      i++;
    }
    if(Macula->phi0_vary[k]) {
      Macula->phi0[k][lc_num] = a[i];
      i++;
    }
    if(Macula->alphamax_vary[k]) {
      Macula->alphamax[k][lc_num] = a[i];
      i++;
    }
    if(Macula->fspot_vary[k]) {
      Macula->fspot[k][lc_num] = a[i];
      i++;
    }
    if(Macula->tmax_vary[k]) {
      Macula->fspot[k][lc_num] = a[i];
      i++;
    }
    if(Macula->life_vary[k]) {
      Macula->life[k][lc_num] = a[i];
      i++;
    }
    if(Macula->ingress_vary[k]) {
      Macula->ingress[k][lc_num] = a[i];
      i++;
    }
    if(Macula->egress_vary[k]) {
      Macula->egress[k][lc_num] = a[i];
      i++;
    }
  }
  
  MaculaFit.blend_vary = Macula->blend_vary;
  if(Macula->blend_vary) {
    Macula->blend[lc_num] = a[i];
    i++;
  }

  /* Generate the model if needed for output or for correcting the input LC */
  testout = 0;
  if(lcoutname != NULL ? (lcoutname[0] != '\0' ? 1 : 0) : 0) {
    if((outfile = fopen(lcoutname,"w")) == NULL) {
      VARTOOLS_error2(ERR_CANNOTWRITE,lcoutname);
    }
    testout = 1;
  }
  if(testout || Macula->correctlc) {
    if((outflux = (double *) malloc(NJD * sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    simulatemaculalc(NJD, t, outflux, MaculaFit.Theta_star,
		     Macula->Nspot, MaculaFit.Theta_spot, a[ma-1],
		     Macula->blend[lc_num], 0, 0, Macula->outputmodel_tdelv,
		     MaculaFit.deltaratio, MaculaFit.dFmod_star,
		     MaculaFit.dFmod_spot, MaculaFit.dFmod_inst,
		     MaculaFit.dFmoddt);

    if(Macula->fluxinput && Macula->fluxoutput) {
      for(i=0; i < NJD; i++) {
	if(!isnan(mag[i])) {
	  if(testout) {
	    fprintf(outfile,"%.17g %.17g %.17g %.17g", t[i], mag[i], err[i], outflux[i]);
	    if(Macula->outputmodel_tdelv) {
	      fprintf(outfile," %.17g", MaculaFit.deltaratio[i]);
	    }
	    fprintf(outfile,"\n");
	  }
	  if(Macula->correctlc) {
	    mag[i] = mag[i]/outflux[i];
	  }
	}
      }
    }
    else if(Macula->fluxinput && !Macula->fluxoutput) {
      for(i=0; i < NJD; i++) {
	if(!isnan(mag[i])) {
	  outmag = 25.0-2.5*log10(outflux[i]);
	  err[i] = 1.0857*err[i]/mag[i];
	  mag[i] = 25.0-2.5*log10(mag[i]);
	  if(testout) {
	    fprintf(outfile,"%.17g %.17g %.17g %.17g", t[i], mag[i], err[i], outmag);
	    if(Macula->outputmodel_tdelv) {
	      fprintf(outfile," %.17g", MaculaFit.deltaratio[i]);
	    }
	    fprintf(outfile,"\n");
	  }
	  if(Macula->correctlc) {
	    mag[i] = mag[i] - outmag;
	  }
	}
      }
    }
    else if(!Macula->fluxinput && Macula->fluxoutput) {
      for(i=0; i < NJD; i++) {
	if(!isnan(mag[i])) {
	  mag[i] = influx[i];
	  err[i] = influx_err[i];
	  if(testout) {
	    fprintf(outfile,"%.17g %.17g %.17g %.17g", t[i], mag[i], err[i], outmag);
	    if(Macula->outputmodel_tdelv) {
	      fprintf(outfile," %.17g", MaculaFit.deltaratio[i]);
	    }
	    fprintf(outfile,"\n");
	  }
	  if(Macula->correctlc) {
	    mag[i] = mag[i] - outmag;
	  }
	}
      }
    }
    else if(!Macula->fluxinput && !Macula->fluxoutput) {
      for(i=0; i < NJD; i++) {
	if(!isnan(mag[i])) {
	  outmag = 25.0-2.5*log10(outflux[i]);
	  if(testout) {
	    fprintf(outfile,"%.17g %.17g %.17g %.17g", t[i], mag[i], err[i], outmag);
	    if(Macula->outputmodel_tdelv) {
	      fprintf(outfile," %.17g", MaculaFit.deltaratio[i]);
	    }
	    fprintf(outfile,"\n");
	  }
	  if(Macula->correctlc) {
	    mag[i] = mag[i] - outmag;
	  }
	}
      }
    }
  }
  if(testout)
    fclose(outfile);

  /* write out the curve file if requested */
  testout = 0;
  if(curvename != NULL ? (curvename[0] != '\0' ? 1 : 0) : 0) {
    if((outfile = fopen(curvename,"w")) == NULL) {
      VARTOOLS_error2(ERR_CANNOTWRITE,curvename);
    }
    testout = 1;
  }
  if(testout) {
    
    if(Macula->outputcurve_stepgiven) {
      tstep = Macula->outputcurve_step;
    }
    else {
      /* Set the step size equal to the minimum non-zero step between
	 points in the input t vector */
      tstep = t[NJD-1] - t[0];
      for(i=1; i < NJD; i++) {
	if(t[i] - t[i-1] > 0 && t[i] - t[i-1] < tstep) {
	  tstep = t[i] - t[i-1];
	}
      }
    }
    Nstep = ceil((t[NJD-1] - t[0])/tstep) + 1;
    if((fluxsim = (double *) malloc(Nstep*sizeof(double))) == NULL ||
       (deltaratio_sim = (double *) malloc(Nstep*sizeof(double))) == NULL ||
       (dFmoddt_sim = (double *) malloc(Nstep*sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    if((tsimcurve = (double *) malloc(Nstep*sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    Nstep = 0;
    ttmp = t[0];
    while(ttmp <= t[NJD-1]) {
      tsimcurve[Nstep] = ttmp;
      ttmp += tstep;
      Nstep++;
    }
    simulatemaculalc(Nstep, tsimcurve, fluxsim, MaculaFit.Theta_star, 
		     Macula->Nspot, 
		     MaculaFit.Theta_spot, a[ma-1], Macula->blend[lc_num], 0, 
		     0, 1, deltaratio_sim, 
		     NULL, NULL, NULL, dFmoddt_sim);
    
    if(!Macula->fluxoutput) {
      for(i=0; i < Nstep; i++) {
	magout = 25.0-2.5*log10(fluxsim[i]);
	fprintf(outfile,"%.17g %.17g", tsimcurve[i], magout);
	if(Macula->outputcurve_tdelv)
	  fprintf(outfile," %.17g", deltaratio_sim[i]);
	fprintf(outfile,"\n");
      }
    } else {
      for(i=0; i < Nstep; i++) {
	fprintf(outfile,"%.17g %.17g", tsimcurve[i], fluxsim[i]);
	if(Macula->outputcurve_tdelv)
	  fprintf(outfile," %.17g", deltaratio_sim[i]);
	fprintf(outfile,"\n");
      }
    }      

    free(fluxsim); free(deltaratio_sim);
    free(dFmoddt_sim); free(tsimcurve);
  }
  if(testout)
    fclose(outfile);

  /**** Free all malloc'd data *****/

  if(a != NULL) free(a);
  if(MaculaFit.Theta_star != NULL) free(MaculaFit.Theta_star);

  if(Macula->Nspot > 0) {
    free(MaculaFit.lambda0_vary);
    free(MaculaFit.phi0_vary);
    free(MaculaFit.alphamax_vary);
    free(MaculaFit.fspot_vary);
    free(MaculaFit.tmax_vary);
    free(MaculaFit.life_vary);
    free(MaculaFit.ingress_vary);
    free(MaculaFit.egress_vary);
    for(j=0; j < MACULA_NDIM_SPOT_PARAM; j++) {
      free(MaculaFit.Theta_spot[j]);
    }
  }
  free(MaculaFit.Theta_spot);
  free(MaculaFit.dFmoddt);
  free(MaculaFit.deltaratio);

  for(i=0; i < NJD; i++) {
    free(MaculaFit.dFmod_star[i]);
    if(Macula->Nspot > 0) {
      for(j=0; j < MACULA_NDIM_SPOT_PARAM; j++) {
	free(MaculaFit.dFmod_spot[i][j]);
      }
    }
    free(MaculaFit.dFmod_spot[i]);
    for(j=0; j < MACULA_NDIM_INST_PARAM; j++) {
      free(MaculaFit.dFmod_inst[i][j]);
    }
    free(MaculaFit.dFmod_inst[i]);
  }
  free(MaculaFit.dFmod_star);
  free(MaculaFit.dFmod_spot);
  free(MaculaFit.dFmod_inst);

  if(!Macula->fluxinput) {
    free(influx);
    free(influx_err);
  }

  for(i=0; i < ma; i++) {
    free(covar[i]);
    free(alpha[i]);
  }
  free(covar);
  free(alpha);
  free(ia);
#ifdef PARALLEL
  free(atry);
  free(beta);
  free(da);
  for(k=0; k < ma; k++)
    free(oneda[k]);
  free(oneda);
#endif

  return;
}

void domaculafit_amoeba(int NJD, double *t, double *mag, double *err, 
		 _Macula *Macula, int lc_num, char *lcoutname, 
		    char *curvename)
{

  double best_chi2, baseflux;
  int Ndof, Nparameters = 0, Ntovary = 0, nfunk, amoeba_val = 0;
  double ftol;

  double **p = NULL, *chi2vals = NULL, delta;

  int *ia = NULL, ma, ma_max, i, j, k;
  double outmag;
  double *outflux = NULL;
  double *influx = NULL, *influx_err = NULL;

  _MaculaLMFitStruct MaculaFit;

  FILE *outfile;

  double macula_evalfunc_amoeba(double *a, int ma, int N, double *t, double *y, double *sig, void *userparams);

  void simulatemaculalc(int N, double *JD, double *fluxout, double *Theta_star, int Nspot, double **Theta_spot, double baseflux, double blend, int derivatives, int temporal, int TdeltaV, double *deltaratio, double **dFmod_star, double ***dFmod_spot, double ***dFmod_inst, double *dFmoddt);

  int iter = 0, testout, Nstep;
  double oalamda, deltachi, tstep, ttmp, magout;

  double *fluxsim = NULL, *deltaratio_sim = NULL, *dFmoddt_sim = NULL,
    *tsimcurve = NULL;

  ma_max = MACULA_NDIM_STAR_PARAM + MACULA_NDIM_SPOT_PARAM*Macula->Nspot 
    + MACULA_NDIM_INST_PARAM;


  /* Setup the MaculaFit structure */
  MaculaFit.Nspot = Macula->Nspot;
  if((MaculaFit.Theta_star = (double *) malloc(MACULA_NDIM_STAR_PARAM * sizeof(double))) == NULL ||
     (MaculaFit.Theta_spot = (double **) malloc(MACULA_NDIM_SPOT_PARAM * sizeof(double *))) == NULL ||
     (MaculaFit.Fmod = (double *) malloc(NJD * sizeof(double))) == NULL ||
     (MaculaFit.dFmoddt = (double *) malloc(NJD * sizeof(double))) == NULL ||
     (MaculaFit.deltaratio = (double *) malloc(NJD * sizeof(double))) == NULL ||
     (MaculaFit.dFmod_star = (double **) malloc(NJD * sizeof(double *))) == NULL ||
     (MaculaFit.dFmod_spot = (double ***) malloc(NJD * sizeof(double **))) == NULL ||
     (MaculaFit.dFmod_inst = (double ***) malloc(NJD * sizeof(double **))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);

  if(Macula->Nspot > 0) {
    if((MaculaFit.lambda0_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
       (MaculaFit.phi0_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
       (MaculaFit.alphamax_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
       (MaculaFit.fspot_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
       (MaculaFit.tmax_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
       (MaculaFit.life_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
       (MaculaFit.ingress_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL ||
       (MaculaFit.egress_vary = (int *) malloc(Macula->Nspot * sizeof(int))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    for(j=0; j < MACULA_NDIM_SPOT_PARAM; j++) {
      if((MaculaFit.Theta_spot[j] = (double *) malloc(Macula->Nspot * sizeof(double))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
    }
  }
  for(i=0; i < NJD; i++) {
    if((MaculaFit.dFmod_star[i] = (double *) malloc(MACULA_NDIM_STAR_PARAM * sizeof(double))) == NULL ||
       (MaculaFit.dFmod_spot[i] = (double **) malloc(MACULA_NDIM_SPOT_PARAM * sizeof(double *))) == NULL ||
       (MaculaFit.dFmod_inst[i] = (double **) malloc(MACULA_NDIM_INST_PARAM * sizeof(double *))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    for(j=0; j < MACULA_NDIM_SPOT_PARAM; j++) {
      if((MaculaFit.dFmod_spot[i][j] = (double *) malloc(Macula->Nspot * sizeof(double))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
    }
    for(j=0; j < MACULA_NDIM_INST_PARAM; j++) {
      if((MaculaFit.dFmod_inst[i][j] = (double *) malloc(sizeof(double))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
    }
  }

  /* Fill in parameters */
  MaculaFit.Theta_star[0] = Macula->Istar[lc_num];
  MaculaFit.Theta_star[1] = Macula->Prot[lc_num];
  MaculaFit.Theta_star[2] = Macula->kappa2[lc_num];
  MaculaFit.Theta_star[3] = Macula->kappa4[lc_num];
  MaculaFit.Theta_star[4] = Macula->c1[lc_num];
  MaculaFit.Theta_star[5] = Macula->c2[lc_num];
  MaculaFit.Theta_star[6] = Macula->c3[lc_num];
  MaculaFit.Theta_star[7] = Macula->c4[lc_num];
  MaculaFit.Theta_star[8] = Macula->d1[lc_num];
  MaculaFit.Theta_star[9] = Macula->d2[lc_num];
  MaculaFit.Theta_star[10] = Macula->d3[lc_num];
  MaculaFit.Theta_star[11] = Macula->d4[lc_num];

  for(k=0; k < Macula->Nspot; k++) {
    MaculaFit.Theta_spot[0][k] = Macula->lambda0[k][lc_num];
    MaculaFit.Theta_spot[1][k] = Macula->phi0[k][lc_num];
    MaculaFit.Theta_spot[2][k] = Macula->alphamax[k][lc_num];
    MaculaFit.Theta_spot[3][k] = Macula->fspot[k][lc_num];
    MaculaFit.Theta_spot[4][k] = Macula->tmax[k][lc_num];
    MaculaFit.Theta_spot[5][k] = Macula->life[k][lc_num];
    MaculaFit.Theta_spot[6][k] = Macula->ingress[k][lc_num];
    MaculaFit.Theta_spot[7][k] = Macula->egress[k][lc_num];
  }

  MaculaFit.blend = Macula->blend[lc_num];

  /* Setup the initial simplex */
  MaculaFit.Istar_vary = Macula->Istar_vary;
  if(Macula->Istar_vary) {
    if(Macula->Istar[lc_num] < 0.5) delta = 0.1;
    else delta = -0.1;
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   1, Macula->Istar[lc_num], delta);
  }
  else {
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   0, Macula->Istar[lc_num], 0.);
  }
  MaculaFit.Prot_vary = Macula->Prot_vary;
  if(Macula->Prot_vary) {
    delta = 0.1*Macula->Prot[lc_num]*Macula->Prot[lc_num]/(t[NJD-1] - t[0]);
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   1, Macula->Prot[lc_num], delta);
  }
  else {
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   0, Macula->Prot[lc_num], 0.);
  }
  MaculaFit.kappa2_vary = Macula->kappa2_vary;
  if(Macula->kappa2_vary) {
    if(Macula->kappa2[lc_num] < 0.1) delta = 0.01;
    else delta = 0.1*Macula->kappa2[lc_num];
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   1, Macula->kappa2[lc_num], delta);
  }
  else {
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   0, Macula->kappa2[lc_num], 0.);
  }
  MaculaFit.kappa4_vary = Macula->kappa4_vary;
  if(Macula->kappa4_vary) {
    if(Macula->kappa4[lc_num] < 0.1) delta = 0.01;
    else delta = 0.1*Macula->kappa4[lc_num];
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   1, Macula->kappa4[lc_num], delta);
  }
  else {
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   0, Macula->kappa4[lc_num], 0.);
  }
  MaculaFit.c1_vary = Macula->c1_vary;
  if(Macula->c1_vary) {
    if(Macula->c1[lc_num] < 0.1) delta = 0.01;
    else delta = 0.1*Macula->c1[lc_num];
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   1, Macula->c1[lc_num], delta);
  }
  else {
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   0, Macula->c1[lc_num], 0.);
  }
  MaculaFit.c2_vary = Macula->c2_vary;
  if(Macula->c2_vary) {
    if(Macula->c2[lc_num] < 0.1) delta = 0.01;
    else delta = 0.1*Macula->c2[lc_num];
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   1, Macula->c2[lc_num], delta);
  }
  else {
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   0, Macula->c2[lc_num], 0.);
  }
  MaculaFit.c3_vary = Macula->c3_vary;
  if(Macula->c3_vary) {
    if(Macula->c3[lc_num] < 0.1) delta = 0.01;
    else delta = 0.1*Macula->c3[lc_num];
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   1, Macula->c3[lc_num], delta);
  }
  else {
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   0, Macula->c3[lc_num], 0.);
  }
  MaculaFit.c4_vary = Macula->c4_vary;
  if(Macula->c4_vary) {
    if(Macula->c4[lc_num] < 0.1) delta = 0.01;
    else delta = 0.1*Macula->c4[lc_num];
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   1, Macula->c4[lc_num], delta);
  }
  else {
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   0, Macula->c4[lc_num], 0.);
  }
  MaculaFit.d1_vary = Macula->d1_vary;
  if(Macula->d1_vary) {
    if(Macula->d1[lc_num] < 0.1) delta = 0.01;
    else delta = 0.1*Macula->d1[lc_num];
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   1, Macula->d1[lc_num], delta);
  }
  else {
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   0, Macula->d1[lc_num], 0.);
  }
  MaculaFit.d2_vary = Macula->d2_vary;
  if(Macula->d2_vary) {
    if(Macula->d2[lc_num] < 0.1) delta = 0.01;
    else delta = 0.1*Macula->d2[lc_num];
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   1, Macula->d2[lc_num], delta);
  }
  else {
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   0, Macula->d2[lc_num], 0.);
  }
  MaculaFit.d3_vary = Macula->d3_vary;
  if(Macula->d3_vary) {
    if(Macula->d3[lc_num] < 0.1) delta = 0.01;
    else delta = 0.1*Macula->d3[lc_num];
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   1, Macula->d3[lc_num], delta);
  }
  else {
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   0, Macula->d3[lc_num], 0.);
  }
  MaculaFit.d4_vary = Macula->d4_vary;
  if(Macula->d4_vary) {
    if(Macula->d4[lc_num] < 0.1) delta = 0.01;
    else delta = 0.1*Macula->d4[lc_num];
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   1, Macula->d4[lc_num], delta);
  }
  else {
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   0, Macula->d4[lc_num], 0.);
  }

  for(k=0; k < Macula->Nspot; k++) {
    MaculaFit.lambda0_vary[k] = Macula->lambda0_vary[k];
    if(Macula->lambda0_vary[k]) {
      if(Macula->lambda0[k][lc_num] < 0.5) delta = 0.1;
      else delta = -0.1;
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					     1, Macula->lambda0[k][lc_num], delta);
    }
    else {
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					     0, Macula->lambda0[k][lc_num], 0.);
    }
    MaculaFit.phi0_vary[k] = Macula->phi0_vary[k];
    if(Macula->phi0_vary[k]) {
      if(Macula->phi0[k][lc_num] < 0.5) delta = 0.1;
      else delta = -0.1;
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					     1, Macula->phi0[k][lc_num], delta);
    }
    else {
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					     0, Macula->phi0[k][lc_num], 0.);
    }
    MaculaFit.alphamax_vary[k] = Macula->alphamax_vary[k];
    if(Macula->alphamax_vary[k]) {
      if(Macula->alphamax[k][lc_num] > 0.0) delta = Macula->alphamax[k][lc_num]*0.1;
      else delta = 0.01;
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					     1, Macula->alphamax[k][lc_num], delta);
    }
    else {
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					     0, Macula->alphamax[k][lc_num], 0.);
    }
    MaculaFit.fspot_vary[k] = Macula->fspot_vary[k];
    if(Macula->fspot_vary[k]) {
      if(Macula->fspot[k][lc_num] > 0.9) delta = -Macula->fspot[k][lc_num]*0.1;
      else if(Macula->fspot[k][lc_num] > 0.1) delta = Macula->fspot[k][lc_num]*0.1;
      else delta = 0.01;
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					     1, Macula->fspot[k][lc_num], delta);
    }
    else {
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					     0, Macula->fspot[k][lc_num], 0.);
    }
    MaculaFit.tmax_vary[k] = Macula->tmax_vary[k];
    if(Macula->tmax_vary[k]) {
      if(Macula->life[k][lc_num] > 0) delta = Macula->life[k][lc_num]*0.05;
      else delta = 0.1;
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					     1, Macula->tmax[k][lc_num], delta);
    }
    else {
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					     0, Macula->tmax[k][lc_num], 0.);
    }
    MaculaFit.life_vary[k] = Macula->life_vary[k];
    if(Macula->life_vary[k]) {
      if(Macula->life[k][lc_num] > 0.1) delta = Macula->life[k][lc_num]*0.1;
      else delta = 0.01;
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					     1, Macula->life[k][lc_num], delta);
    }
    else {
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					     0, Macula->life[k][lc_num], 0.);
    }
    MaculaFit.ingress_vary[k] = Macula->ingress_vary[k];
    if(Macula->ingress_vary[k]) {
      if(Macula->ingress[k][lc_num] > 0) delta = Macula->ingress[k][lc_num]*0.1;
      else if(Macula->life[k][lc_num] > 0) delta = Macula->life[k][lc_num]*0.01;
      else delta = 0.01;
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					     1, Macula->ingress[k][lc_num], delta);
    }
    else {
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					     0, Macula->ingress[k][lc_num], 0.);
    }
    MaculaFit.egress_vary[k] = Macula->egress_vary[k];
    if(Macula->egress_vary[k]) {
      if(Macula->egress[k][lc_num] > 0) delta = Macula->egress[k][lc_num]*0.1;
      else if(Macula->life[k][lc_num] > 0) delta = Macula->life[k][lc_num]*0.01;
      else delta = 0.01;
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					     1, Macula->egress[k][lc_num], delta);
    }
    else {
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					     0, Macula->egress[k][lc_num], 0.);
    }
  }

  MaculaFit.blend_vary = Macula->blend_vary;
  if(Macula->blend_vary) {
    if(Macula->blend[lc_num] < 0.1) delta = 0.01;
    else if(Macula->blend[lc_num] >= 0.1 && Macula->blend[lc_num] <= 0.9)
      delta = 0.1*Macula->blend[lc_num];
    else
      delta = -0.1*Macula->blend[lc_num];
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   1, Macula->blend[lc_num], delta);
  }
  else {
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					   0, Macula->blend[lc_num], 0.);
  }

  /* Convert mag to flux if needed */
  if(!Macula->fluxinput) {
    if((influx = (double *) malloc(NJD * sizeof(double))) == NULL ||
       (influx_err = (double *) malloc(NJD * sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);

    for(i=0; i < NJD; i++) {
      influx[i] = pow(10.0,(-0.4*(mag[i] - 25.0)));
      influx_err[i] = err[i]*influx[i]/1.0857;
    }
  } else {
    influx = mag;
    influx_err = err;
  }

  /* Set the base flux level to the 95 percentile flux value */
  baseflux = VARTOOLS_percentile(NJD, influx, 95.0);
  delta = VARTOOLS_stddev(NJD, influx)/sqrt((double) NJD);
  VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					 1, baseflux, delta);

  Ndof = NJD - Ntovary;
  ftol = DELTACHI_CONVERGENCELIMIT;
  /* Get the chi2vals for the points on the simplex */
  VARTOOLS_amoeba_initializesimplexchi2(Nparameters, Ntovary, p, &chi2vals,
					&macula_evalfunc_amoeba, NJD, t,
					influx, influx_err,
					(void *) (&MaculaFit));
  /* Do the fit */
  amoeba_val = VARTOOLS_amoeba(p, chi2vals, ia, Nparameters, ftol, &macula_evalfunc_amoeba, &nfunk, 0, NJD, t, influx, influx_err, (void *) (&MaculaFit));
  if(1) {
    /* The solution converged, update the parameters */
    best_chi2 = chi2vals[0];
    j = 0;
    for(i=1; i < Ntovary+1; i++) {
      if(chi2vals[i] < best_chi2) {
	best_chi2 = chi2vals[i];
	j = i;
      }
    }
    
    /* The best parameters are stored in p[j] */
    
    Macula->chi2val[lc_num] = best_chi2;
    Macula->Ndof[lc_num] = Ndof;
    
    /* Update the parameters */
    
    i = 0;
    if(Macula->Istar_vary) {
      Macula->Istar[lc_num] = p[j][i];
    }
    i++;
    if(Macula->Prot_vary) {
      Macula->Prot[lc_num] = p[j][i];
    }
    i++;
    if(Macula->kappa2_vary) {
      Macula->kappa2[lc_num] = p[j][i];
    }
    i++;
    if(Macula->kappa4_vary) {
      Macula->kappa4[lc_num] = p[j][i];
    }
    i++;
    if(Macula->c1_vary) {
      Macula->c1[lc_num] = p[j][i];
    }
    i++;
    if(Macula->c2_vary) {
      Macula->c2[lc_num] = p[j][i];
    }
    i++;
    if(Macula->c3_vary) {
      Macula->c3[lc_num] = p[j][i];
    }
    i++;
    if(Macula->c4_vary) {
      Macula->c4[lc_num] = p[j][i];
    }
    i++;
    if(Macula->d1_vary) {
      Macula->d1[lc_num] = p[j][i];
    }
    i++;
    if(Macula->d2_vary) {
      Macula->d2[lc_num] = p[j][i];
    }
    i++;
    if(Macula->d3_vary) {
      Macula->d3[lc_num] = p[j][i];
    }
    i++;
    if(Macula->d4_vary) {
      Macula->d4[lc_num] = p[j][i];
    }
    i++;
    
    for(k=0; k < Macula->Nspot; k++) {
      if(Macula->lambda0_vary[k]) {
	Macula->lambda0[k][lc_num] = p[j][i];
      }
      i++;
      if(Macula->phi0_vary[k]) {
	Macula->phi0[k][lc_num] = p[j][i];
      }
      i++;
      if(Macula->alphamax_vary[k]) {
	Macula->alphamax[k][lc_num] = p[j][i];
      }
      i++;
      if(Macula->fspot_vary[k]) {
	Macula->fspot[k][lc_num] = p[j][i];
      }
      i++;
      if(Macula->tmax_vary[k]) {
	Macula->fspot[k][lc_num] = p[j][i];
      }
      i++;
      if(Macula->life_vary[k]) {
	Macula->life[k][lc_num] = p[j][i];
      }
      i++;
      if(Macula->ingress_vary[k]) {
	Macula->ingress[k][lc_num] = p[j][i];
      }
      i++;
      if(Macula->egress_vary[k]) {
	Macula->egress[k][lc_num] = p[j][i];
      }
      i++;
    }
    
    if(Macula->blend_vary) {
      Macula->blend[lc_num] = p[j][i];
    }
    i++;
  }
  else {
    Macula->chi2val[lc_num] = -1.;
  }
  
  /* update the best-fit parameters in the MaculaFit struct */
  Macula->chi2val[lc_num]  = macula_evalfunc_amoeba(p[j], Nparameters, NJD, t, influx, influx_err, (void *) (&MaculaFit));
  Macula->Ndof[lc_num ] = NJD - Ntovary;

  /* Generate the model if needed for output or for correcting the input LC */
  testout = 0;
  if(lcoutname != NULL ? (lcoutname[0] != '\0' ? 1 : 0) : 0) {
    if((outfile = fopen(lcoutname,"w")) == NULL) {
      VARTOOLS_error2(ERR_CANNOTWRITE,lcoutname);
    }
    testout = 1;
  }
  if(testout || Macula->correctlc) {
    if((outflux = (double *) malloc(NJD * sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    simulatemaculalc(NJD, t, outflux, MaculaFit.Theta_star,
		     Macula->Nspot, MaculaFit.Theta_spot, p[j][Nparameters-1],
		     Macula->blend[lc_num], 0, 0, Macula->outputmodel_tdelv,
		     MaculaFit.deltaratio, MaculaFit.dFmod_star,
		     MaculaFit.dFmod_spot, MaculaFit.dFmod_inst,
		     MaculaFit.dFmoddt);

    if(Macula->fluxinput && Macula->fluxoutput) {
      for(i=0; i < NJD; i++) {
	if(!isnan(mag[i])) {
	  if(testout) {
	    fprintf(outfile,"%.17g %.17g %.17g %.17g", t[i], mag[i], err[i], outflux[i]);
	    if(Macula->outputmodel_tdelv) {
	      fprintf(outfile," %.17g", MaculaFit.deltaratio[i]);
	    }
	    fprintf(outfile,"\n");
	  }
	  if(Macula->correctlc) {
	    mag[i] = mag[i]/outflux[i];
	  }
	}
      }
    }
    else if(Macula->fluxinput && !Macula->fluxoutput) {
      for(i=0; i < NJD; i++) {
	if(!isnan(mag[i])) {
	  outmag = 25.0-2.5*log10(outflux[i]);
	  err[i] = 1.0857*err[i]/mag[i];
	  mag[i] = 25.0-2.5*log10(mag[i]);
	  if(testout) {
	    fprintf(outfile,"%.17g %.17g %.17g %.17g", t[i], mag[i], err[i], outmag);
	    if(Macula->outputmodel_tdelv) {
	      fprintf(outfile," %.17g", MaculaFit.deltaratio[i]);
	    }
	    fprintf(outfile,"\n");
	  }
	  if(Macula->correctlc) {
	    mag[i] = mag[i] - outmag;
	  }
	}
      }
    }
    else if(!Macula->fluxinput && Macula->fluxoutput) {
      for(i=0; i < NJD; i++) {
	if(!isnan(mag[i])) {
	  mag[i] = influx[i];
	  err[i] = influx_err[i];
	  if(testout) {
	    fprintf(outfile,"%.17g %.17g %.17g %.17g", t[i], mag[i], err[i], outmag);
	    if(Macula->outputmodel_tdelv) {
	      fprintf(outfile," %.17g", MaculaFit.deltaratio[i]);
	    }
	    fprintf(outfile,"\n");
	  }
	  if(Macula->correctlc) {
	    mag[i] = mag[i] - outmag;
	  }
	}
      }
    }
    else if(!Macula->fluxinput && !Macula->fluxoutput) {
      for(i=0; i < NJD; i++) {
	if(!isnan(mag[i])) {
	  outmag = 25.0-2.5*log10(outflux[i]);
	  if(testout) {
	    fprintf(outfile,"%.17g %.17g %.17g %.17g", t[i], mag[i], err[i], outmag);
	    if(Macula->outputmodel_tdelv) {
	      fprintf(outfile," %.17g", MaculaFit.deltaratio[i]);
	    }
	    fprintf(outfile,"\n");
	  }
	  if(Macula->correctlc) {
	    mag[i] = mag[i] - outmag;
	  }
	}
      }
    }
  }
  if(testout)
    fclose(outfile);

  /* write out the curve file if requested */
  testout = 0;
  if(curvename != NULL ? (curvename[0] != '\0' ? 1 : 0) : 0) {
    if((outfile = fopen(curvename,"w")) == NULL) {
      VARTOOLS_error2(ERR_CANNOTWRITE,curvename);
    }
    testout = 1;
  }
  if(testout) {
    
    if(Macula->outputcurve_stepgiven) {
      tstep = Macula->outputcurve_step;
    }
    else {
      /* Set the step size equal to the minimum non-zero step between
	 points in the input t vector */
      tstep = t[NJD-1] - t[0];
      for(i=1; i < NJD; i++) {
	if(t[i] - t[i-1] > 0 && t[i] - t[i-1] < tstep) {
	  tstep = t[i] - t[i-1];
	}
      }
    }
    Nstep = ceil((t[NJD-1] - t[0])/tstep) + 1;
    if((fluxsim = (double *) malloc(Nstep*sizeof(double))) == NULL ||
       (deltaratio_sim = (double *) malloc(Nstep*sizeof(double))) == NULL ||
       (dFmoddt_sim = (double *) malloc(Nstep*sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    if((tsimcurve = (double *) malloc(Nstep*sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    Nstep = 0;
    ttmp = t[0];
    while(ttmp <= t[NJD-1]) {
      tsimcurve[Nstep] = ttmp;
      ttmp += tstep;
      Nstep++;
    }
    simulatemaculalc(Nstep, tsimcurve, fluxsim, MaculaFit.Theta_star, 
		     Macula->Nspot, 
		     MaculaFit.Theta_spot, p[j][Nparameters-1], Macula->blend[lc_num], 0, 
		     0, 1, deltaratio_sim, 
		     NULL, NULL, NULL, dFmoddt_sim);
    
    if(!Macula->fluxoutput) {
      for(i=0; i < Nstep; i++) {
	magout = 25.0-2.5*log10(fluxsim[i]);
	fprintf(outfile,"%.17g %.17g", tsimcurve[i], magout);
	if(Macula->outputcurve_tdelv)
	  fprintf(outfile," %.17g", deltaratio_sim[i]);
	fprintf(outfile,"\n");
      }
    } else {
      for(i=0; i < Nstep; i++) {
	fprintf(outfile,"%.17g %.17g", tsimcurve[i], fluxsim[i]);
	if(Macula->outputcurve_tdelv)
	  fprintf(outfile," %.17g", deltaratio_sim[i]);
	fprintf(outfile,"\n");
      }
    }      

    free(fluxsim); free(deltaratio_sim);
    free(dFmoddt_sim); free(tsimcurve);
  }
  if(testout)
    fclose(outfile);

  /**** Free all malloc'd data *****/

  VARTOOLS_amoeba_cleanup(&Nparameters, &Ntovary, &p, &ia, &chi2vals);
  if(MaculaFit.Theta_star != NULL) free(MaculaFit.Theta_star);

  if(Macula->Nspot > 0) {
    free(MaculaFit.lambda0_vary);
    free(MaculaFit.phi0_vary);
    free(MaculaFit.alphamax_vary);
    free(MaculaFit.fspot_vary);
    free(MaculaFit.tmax_vary);
    free(MaculaFit.life_vary);
    free(MaculaFit.ingress_vary);
    free(MaculaFit.egress_vary);
    for(j=0; j < MACULA_NDIM_SPOT_PARAM; j++) {
      free(MaculaFit.Theta_spot[j]);
    }
  }
  free(MaculaFit.Theta_spot);
  free(MaculaFit.dFmoddt);
  free(MaculaFit.deltaratio);
  free(MaculaFit.Fmod);

  for(i=0; i < NJD; i++) {
    free(MaculaFit.dFmod_star[i]);
    if(Macula->Nspot > 0) {
      for(j=0; j < MACULA_NDIM_SPOT_PARAM; j++) {
	free(MaculaFit.dFmod_spot[i][j]);
      }
    }
    free(MaculaFit.dFmod_spot[i]);
    for(j=0; j < MACULA_NDIM_INST_PARAM; j++) {
      free(MaculaFit.dFmod_inst[i][j]);
    }
    free(MaculaFit.dFmod_inst[i]);
  }
  free(MaculaFit.dFmod_star);
  free(MaculaFit.dFmod_spot);
  free(MaculaFit.dFmod_inst);

  if(!Macula->fluxinput) {
    free(influx);
    free(influx_err);
  }


  return;
}


  void macula_evalfunc_lm(double *t, double *a, double *yfit, double **dyda, int ma, int N, int Nlin_coeffs, double *lin_coeffs, double **Design_Matrix, double *y, double *sig, int *varylin_coeffs, void *userparams)
{
  /****** This function is in the form expected by mrqmin (the
levenberg-marquardt minimization routine. It is used to generate a
macula model for a given set of input parameters. *****/

  /* Parameters:

     t = time of observations.
     a = non-linear input parameters.
     yfit = the output model light curve.
     ma = the number of parameters.
     N = the number of data points.
     Nlin_coeffs through varylin_coeffs = optional coefficients that are 
        linear.  These are ignored here.
     userparams = the Macula struct

  */

  _MaculaLMFitStruct *MaculaFit;

  double *deltaratio = NULL;

  double *Theta_star = NULL, **Theta_spot = NULL, *dFmoddt = NULL,
    **dFmod_star = NULL, ***dFmod_spot = NULL, ***dFmod_inst = NULL;

  double baseflux, blend, derivatives;
  
  int temporal, TdeltaV;

  int i, j, k, acnt;

  void simulatemaculalc(int N, double *JD, double *fluxout, double *Theta_star, int Nspot, double **Theta_spot, double baseflux, double blend, int derivatives, int temporal, int TdeltaV, double *deltaratio, double **dFmod_star, double ***dFmod_spot, double ***dFmod_inst, double *dFmoddt);

  MaculaFit = (_MaculaLMFitStruct *) userparams;

  deltaratio = MaculaFit->deltaratio;
  Theta_star = MaculaFit->Theta_star;
  Theta_spot = MaculaFit->Theta_spot;
  dFmoddt = MaculaFit->dFmoddt;
  dFmod_star = MaculaFit->dFmod_star;
  dFmod_spot = MaculaFit->dFmod_spot;
  dFmod_inst = MaculaFit->dFmod_inst;
     
  /* Fill in varied parameters */
  acnt = 0;
  if(MaculaFit->Istar_vary) {
    Theta_star[0] = a[acnt];
    acnt++;
  }
  if(MaculaFit->Prot_vary) {
    Theta_star[1] = a[acnt];
    acnt++;
  }
  if(MaculaFit->kappa2_vary) {
    Theta_star[2] = a[acnt];
    acnt++;
  }
  if(MaculaFit->kappa4_vary) {
    Theta_star[3] = a[acnt];
    acnt++;
  }
  if(MaculaFit->c1_vary) {
    Theta_star[4] = a[acnt];
    acnt++;
  }
  if(MaculaFit->c2_vary) {
    Theta_star[5] = a[acnt];
    acnt++;
  }
  if(MaculaFit->c3_vary) {
    Theta_star[6] = a[acnt];
    acnt++;
  }
  if(MaculaFit->c4_vary) {
    Theta_star[7] = a[acnt];
    acnt++;
  }
  if(MaculaFit->d1_vary) {
    Theta_star[8] = a[acnt];
    acnt++;
  }
  if(MaculaFit->d2_vary) {
    Theta_star[9] = a[acnt];
    acnt++;
  }
  if(MaculaFit->d3_vary) {
    Theta_star[10] = a[acnt];
    acnt++;
  }
  if(MaculaFit->d4_vary) {
    Theta_star[11] = a[acnt];
    acnt++;
  }

  for(k=0; k < MaculaFit->Nspot; k++) {
    if(MaculaFit->lambda0_vary[k]) {
      Theta_spot[0][k] = a[acnt];
      acnt++;
    }
    if(MaculaFit->phi0_vary[k]) {
      Theta_spot[1][k] = a[acnt];
      acnt++;
    }
    if(MaculaFit->alphamax_vary[k]) {
      Theta_spot[2][k] = a[acnt];
      acnt++;
    }
    if(MaculaFit->fspot_vary[k]) {
      Theta_spot[3][k] = a[acnt];
      acnt++;
    }
    if(MaculaFit->tmax_vary[k]) {
      Theta_spot[4][k] = a[acnt];
      acnt++;
    }
    if(MaculaFit->life_vary[k]) {
      Theta_spot[5][k] = a[acnt];
      acnt++;
    }
    if(MaculaFit->ingress_vary[k]) {
      Theta_spot[6][k] = a[acnt];
      acnt++;
    }
    if(MaculaFit->egress_vary[k]) {
      Theta_spot[7][k] = a[acnt];
      acnt++;
    }
  }

  if(MaculaFit->blend_vary) {
    blend = a[acnt];
    acnt++;
  } else {
    blend = MaculaFit->blend;
  }
  
  baseflux = a[acnt];

  /* Simulate the lc */
  simulatemaculalc(N, t, yfit, Theta_star, MaculaFit->Nspot, Theta_spot,
		   baseflux, blend, 1, 0, 0, deltaratio, dFmod_star,
		   dFmod_spot, dFmod_inst, dFmoddt);

  /* fill out the dyda vector based on the parameters that are being varied */
  for(i=0; i < N; i++) {
    acnt = 0;
    if(MaculaFit->Istar_vary) {
      dyda[i][acnt] = dFmod_star[i][0];
      acnt++;
    }
    if(MaculaFit->Prot_vary) {
      dyda[i][acnt] = dFmod_star[i][1];
      acnt++;
    }
    if(MaculaFit->kappa2_vary) {
      dyda[i][acnt] = dFmod_star[i][2];
      acnt++;
    }
    if(MaculaFit->kappa4_vary) {
      dyda[i][acnt] = dFmod_star[i][3];
      acnt++;
    }
    if(MaculaFit->c1_vary) {
      dyda[i][acnt] = dFmod_star[i][4];
      acnt++;
    }
    if(MaculaFit->c2_vary) {
      dyda[i][acnt] = dFmod_star[i][5];
      acnt++;
    }
    if(MaculaFit->c3_vary) {
      dyda[i][acnt] = dFmod_star[i][6];
      acnt++;
    }
    if(MaculaFit->c4_vary) {
      dyda[i][acnt] = dFmod_star[i][7];
      acnt++;
    }
    if(MaculaFit->d1_vary) {
      dyda[i][acnt] = dFmod_star[i][4];
      acnt++;
    }
    if(MaculaFit->d2_vary) {
      dyda[i][acnt] = dFmod_star[i][5];
      acnt++;
    }
    if(MaculaFit->d3_vary) {
      dyda[i][acnt] = dFmod_star[i][6];
      acnt++;
    }
    if(MaculaFit->d4_vary) {
      dyda[i][acnt] = dFmod_star[i][7];
      acnt++;
    }
    
    for(k=0; k < MaculaFit->Nspot; k++) {
      if(MaculaFit->lambda0_vary[k]) {
	dyda[i][acnt] = dFmod_spot[i][0][k];
	acnt++;
      }
      if(MaculaFit->phi0_vary[k]) {
	dyda[i][acnt] = dFmod_spot[i][1][k];
	acnt++;
      }
      if(MaculaFit->alphamax_vary[k]) {
	dyda[i][acnt] = dFmod_spot[i][2][k];
	acnt++;
      }
      if(MaculaFit->fspot_vary[k]) {
	dyda[i][acnt] = dFmod_spot[i][3][k];
	acnt++;
      }
      if(MaculaFit->tmax_vary[k]) {
	dyda[i][acnt] = dFmod_spot[i][4][k];
	acnt++;
      }
      if(MaculaFit->life_vary[k]) {
	dyda[i][acnt] = dFmod_spot[i][5][k];
	acnt++;
      }
      if(MaculaFit->ingress_vary[k]) {
	dyda[i][acnt] = dFmod_spot[i][6][k];
	acnt++;
      }
      if(MaculaFit->egress_vary[k]) {
	dyda[i][acnt] = dFmod_spot[i][7][k];
	acnt++;
      }
    }
    
    if(MaculaFit->blend_vary) {
      dyda[i][acnt] = dFmod_inst[i][1][0];
      acnt++;
    }

    /* Base flux partial derivative */
    dyda[i][acnt] = dFmod_inst[i][0][0];

  }

  return;

}


double macula_evalfunc_amoeba(double *p, int Nparam, int N, double *t, double *y, double *sig, void *userparams)
{
  /****** This function is in the form expected by amoeba.. It is used
to generate a macula model for a given set of input parameters. The output
is chi2. *****/

  /* Parameters:

     p = non-linear input parameters.
     Nparam = the number of parameters.
     N = the number of data points.
     t = time of observations.
     y = the observations.
     sig = the uncertainties.
     userparams = the MaculaFit struct

  */

  _MaculaLMFitStruct *MaculaFit;

  double *yfit = NULL;
  double *deltaratio = NULL;

  double *Theta_star = NULL, **Theta_spot = NULL, *dFmoddt = NULL,
    **dFmod_star = NULL, ***dFmod_spot = NULL, ***dFmod_inst = NULL;

  double baseflux, blend, derivatives, chi2val;
  
  int temporal, TdeltaV;

  int i, j, k, acnt;

  void simulatemaculalc(int N, double *JD, double *fluxout, double *Theta_star, int Nspot, double **Theta_spot, double baseflux, double blend, int derivatives, int temporal, int TdeltaV, double *deltaratio, double **dFmod_star, double ***dFmod_spot, double ***dFmod_inst, double *dFmoddt);

  MaculaFit = (_MaculaLMFitStruct *) userparams;

  deltaratio = MaculaFit->deltaratio;
  Theta_star = MaculaFit->Theta_star;
  Theta_spot = MaculaFit->Theta_spot;
  dFmoddt = MaculaFit->dFmoddt;
  dFmod_star = MaculaFit->dFmod_star;
  dFmod_spot = MaculaFit->dFmod_spot;
  dFmod_inst = MaculaFit->dFmod_inst;
  yfit = MaculaFit->Fmod;

  /* Fill in varied parameters */
  for(acnt = 0; acnt <= 11; acnt++)
    Theta_star[acnt] = p[acnt];

  for(k=0; k < MaculaFit->Nspot; k++) {
    for(j=0; j <= 7; j++, acnt++) {
      Theta_spot[j][k] = p[acnt];
    }
  }
  blend = p[acnt];
  acnt++;
  baseflux = p[acnt];

  /* Simulate the lc */
  simulatemaculalc(N, t, yfit, Theta_star, MaculaFit->Nspot, Theta_spot,
		   baseflux, blend, 0, 0, 0, deltaratio, dFmod_star,
		   dFmod_spot, dFmod_inst, dFmoddt);

  chi2val = 0.;
  for(i=0; i < N; i++) {
    chi2val += (y[i] - yfit[i])*(y[i] - yfit[i])/sig[i]/sig[i];
  }

  return chi2val;

}


void simulatemaculalc(int N, double *JD, double *fluxout, double *Theta_star, int Nspot, double **Theta_spot, double baseflux, double blend, int derivatives, int temporal, int TdeltaV, double *deltaratio, double **dFmod_star, double ***dFmod_spot, double ***dFmod_inst, double *dFmoddt)
{
  /* These will store the matrices to interface with the fortran
     routine, note that multi-dimen arrays in fortran are declared as
     single-dimen arrays in c */

  double *Theta_spot_fortran = NULL, *Theta_inst_fortran = NULL,
    *dFmod_star_fortran = NULL, *dFmod_spot_fortran = NULL, 
    *dFmod_inst_fortran = NULL;

  double tstart, tend;

  int i, j, k, Nsets;

  size_t size_Theta_spot_fortran, size_Theta_inst_fortran,
    size_dFmod_star_fortran,
    size_dFmod_spot_fortran, size_dFmod_inst_fortran;

  size_Theta_spot_fortran = MACULA_NDIM_SPOT_PARAM * Nspot + 1;
  size_Theta_inst_fortran = MACULA_NDIM_INST_PARAM + 1;
  size_dFmod_star_fortran = N * MACULA_NDIM_STAR_PARAM + 1;
  size_dFmod_spot_fortran = N * MACULA_NDIM_SPOT_PARAM * Nspot + 1;
  size_dFmod_inst_fortran = N * MACULA_NDIM_INST_PARAM + 1;

  if(size_Theta_spot_fortran > 0) {
    if((Theta_spot_fortran = (double *) malloc(size_Theta_spot_fortran * sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
  }
  if(size_Theta_inst_fortran > 0) {
    if((Theta_inst_fortran = (double *) malloc(size_Theta_inst_fortran * sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
  }
  if(size_dFmod_star_fortran > 0) {
    if((dFmod_star_fortran = (double *) malloc(size_dFmod_star_fortran * sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
  }
  if(size_dFmod_spot_fortran > 0) {
    if((dFmod_spot_fortran = (double *) malloc(size_dFmod_spot_fortran * sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
  }
  if(size_dFmod_inst_fortran > 0) {
    if((dFmod_inst_fortran = (double *) malloc(size_dFmod_inst_fortran * sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
  }

  /* translate the Theta_spot parameters and inst parameters into
     fortran vector format */
  for(j=0; j < MACULA_NDIM_SPOT_PARAM; j++) {
    for(k=0; k < Nspot; k++) {
      Theta_spot_fortran[j + MACULA_NDIM_SPOT_PARAM*k] = Theta_spot[j][k];
    }
  }
  Theta_inst_fortran[0] = baseflux;
  Theta_inst_fortran[1] = blend;

  /* Other parameters used by macula */
  tstart = JD[0] - 0.001; tend = JD[N-1] + 0.001; Nsets = 1;

  /* call macula */
  MACULAFUNC (JD,&N,&Nspot,&Nsets,&derivatives,&temporal,&TdeltaV,
	     Theta_star, Theta_spot_fortran, Theta_inst_fortran, 
	     &tstart, &tend, fluxout, dFmod_star_fortran, 
	     dFmod_spot_fortran, dFmod_inst_fortran, dFmoddt, 
	     deltaratio);
  
  /* Convert the fortran matrices into c format */
  if(derivatives) {
    for(i=0; i < N; i++) {
      for(j=0; j < MACULA_NDIM_STAR_PARAM; j++)
	dFmod_star[i][j] = dFmod_star_fortran[i + N*j];
    
      for(j=0; j < MACULA_NDIM_SPOT_PARAM; j++) {
	for(k=0; k < Nspot; k++)
	  dFmod_spot[i][j][k] = 
	    dFmod_spot_fortran[i + N*(j + MACULA_NDIM_SPOT_PARAM*k)];
      }

      for(j=0; j < MACULA_NDIM_INST_PARAM; j++) {
	dFmod_inst[i][j][0] = dFmod_inst_fortran[i + N*j];
      }
    }
  }

  /* Free the temporary arrays */
  if(Theta_spot_fortran != NULL) free(Theta_spot_fortran);
  if(Theta_inst_fortran != NULL) free(Theta_inst_fortran);
  if(dFmod_star_fortran != NULL) free(dFmod_star_fortran);
  if(dFmod_spot_fortran != NULL) free(dFmod_spot_fortran);
  if(dFmod_inst_fortran != NULL) free(dFmod_inst_fortran);

  return;
}
