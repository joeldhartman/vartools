#include "../../src/vartools.h"
#include <stdio.h>
#include <stdlib.h>
#include "jktebop.h"
#include "jktebop_lib.h"

/* STOPPED WORKING HERE ---- Check for TBD below */

/* This is the source code for a sample user-defined command
   to be used with vartools.

   This library defines the command -jktebop which can be used
   to fit a detached eclipsing binary model to a light curve.

   The routine uses John Southworth's implementation of the EBOP
   program (http://www.astro.keele.ac.uk/jkt/codes/jktebop.html)
   originally due to Popper & Etzel (1981AJ.....86..102), Etzel
   (1981psbs.conf..111E), and Nelson & Davis
   (1972ApJ...174..617N). It has been modified to call it as a library
   rather than running it as a stand-alone program. We only use the
   model generation routines in the program, the model fitting is
   handled separately within VARTOOLS.

   This is a fairly involved example which, among other things,
   illustrates the inclusion of FORTRAN routines in a VARTOOLS
   library.

*/

void jktebop_Initialize(char *commandname,
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
  sprintf(commandname,"-jktebop");

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
     the magnitudes of each lc), below you would replace "_Jktebop"
     with the type name of your structure. This tells VARTOOLS how
     much memory to allocate for the structure storing the data for
     your command.

     See jktebop.h
*/
  *sizeuserdata = sizeof(_Jktebop);
}

int jktebop_ParseCL(ProgramData *p, Command *c,
			   void *userdata, int *iret, char **argv, int argc)
/* Every library with the name $LIBNAME.so must contain a function
   with the name $LIBNAME_ParseCL to check that the user has issued a
   valid syntax, and set-up the command-data accordingly.

   Note that any vectors/arrays used by the command should also be registered
   by this function so that memory will be allocated as needed.

Expected syntax for this example:
   -jktebop < "inject" | "fit" > 
      <"Period" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr> 
          ["vary"]> 
      <"T0" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
          ["vary"]>
      <"r1+r2" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
          ["vary"]>
      <"r2/r1" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
          ["vary"]>
      <"M2/M1" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
          ["vary"]>
      <"J2/J1" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
          ["vary"]>
      <<"i" | "bimpact"> <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
          ["vary"]>
      <"esinomega" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
          ["vary"]>
      <"ecosomega" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
          ["vary"]>
      <"LD1" <"linear" | "quad" | "log" | "sqrt"> 
          <"fix" value1 [value2] | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
          ["vary"]>
      <"LD2" <"lockLD1" | "linear" | "quad" | "log" | "sqrt">
          [<"fix" value [value2] | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
          ["vary"]>
      ["gravdark1" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
          ["vary"]]
      ["gravdark2" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
          ["vary"]]
      ["reflection1" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
          ["vary"]]
      ["reflection2" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
          ["vary"]]
      ["L3" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
          ["vary"]]
      ["tidallag" <"fix" value | "list" | "fixcolumn" <colname | colnum> | \"expr\" expr>
          ["vary"]]
      ["correctlc"]
      ["omodel" <outdir ["format" fmt]>]
      ["ocurve" <"jd" | "phase"> ["step" stepsize]
          <"outdir" outdir [ "format" fmt]>]

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
  int i = 0;
  int check;
  int Nldcoeff;

  double gravdark, reflection, L3, tidallag;

  /* Cast the input userdata pointer to the data structure type used
     for this command */
  _Jktebop *Jktebop;
  Jktebop = (_Jktebop *) userdata;

  /* We'll use i rather than iret to index argv, this is just so we
     don't have to constantly dereference the pointer iret */
  i = *iret;

  /* Our procedure here is to step through each term in the expected
     command-line syntax. When possible we use the 
     VARTOOLS_ParseFixSpecFixcolumn function to parse terms like

     <"fix" value | "list" | "fixcolumn" <colname | colnum> | "expr" expr>

     For other terms we first check to make sure that we haven't
     exceeded the number of command-line arguments (checks like if(i
     >= argc) return 1;), and then use strcmp to test whether the
     command-line term in argv[i] is equal to the expected string (note
     that strcmp returns 0 if the strings are equal). You can use the
     atof() or atoi() functions to parse a string into a double or
     integer respectively. If something on the command-line is not as
     expected return 1 to tell VARTOOLS that there was an error in
     parsing the command. */

  /* Manually check if "inject" or "fit" is set */
  if(i >= argc) {
    /* There are no terms left on the command-line, 
       return an error */
    return 1;
  }
  if(!strcmp(argv[i],"inject")) {
    /* The user gave inject on the command-line.
       We'll use the variable injectorfit in the _Jktebop struct
       to keep track of this */
    Jktebop->injectorfit = JKTEBOP_INJECT;
  } else if(!strcmp(argv[i],"fit")) {
    Jktebop->injectorfit = JKTEBOP_FIT;
  } else {
    /* The command requires the user to specify "inject" or "fit"
       if the next term is not one of these, return an error */
    return 1;
  }

  i++;

  /* Now parse the options for each of the fitting parameters */
  
  if(i >= argc) {
    return 1;
  }
  if(strcmp(argv[i],"Period")) {
    *iret = i; return 1;
  } else {
    i++;
  }
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
				VARTOOLS_TYPE_DOUBLE, 
				(void *) (&Jktebop->Period),
				0,
				1,
				"PERIOD");
  if(check) {*iret = i; return 1;}
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Jktebop->Period_vary = 1;
      i++;
    } else {
      Jktebop->Period_vary = 0;
    }
  }

  if(i >= argc)
    return 1;
  if(strcmp(argv[i],"T0")) {
    *iret = i; return 1;
  } else {
    i++;
  }
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
				VARTOOLS_TYPE_DOUBLE, 
				(void *) (&Jktebop->T0),
				0,
				1,
				"T0");
  if(check) {*iret = i; return 1;}
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Jktebop->T0_vary = 1;
      i++;
    } else {
      Jktebop->T0_vary = 0;
    }
  }


  if(i >= argc)
    return 1;
  if(strcmp(argv[i],"r1+r2")) {
    *iret = i; return 1;
  } else {
    i++;
  }
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
				VARTOOLS_TYPE_DOUBLE, 
				(void *) (&Jktebop->r1r2),
				0,
				1,
				"R1+R2");
  if(check) {*iret = i; return 1;}
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Jktebop->r1r2_vary = 1;
      i++;
    } else {
      Jktebop->r1r2_vary = 0;
    }
  }
  

  if(i >= argc)
    return 1;
  if(strcmp(argv[i],"r2/r1")) {
    *iret = i; return 1;
  } else {
    i++;
  }
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
				VARTOOLS_TYPE_DOUBLE, 
				(void *) (&Jktebop->r2_r1),
				0,
				1,
				"R2/R1");
  if(check) {*iret = i; return 1;}
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Jktebop->r2_r1_vary = 1;
      i++;
    } else {
      Jktebop->r2_r1_vary = 0;
    }
  }


  if(i >= argc)
    return 1;
  if(strcmp(argv[i],"M2/M1")) {
    *iret = i; return 1;
  } else {
    i++;
  }
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
				VARTOOLS_TYPE_DOUBLE, 
				(void *) (&Jktebop->M2_M1),
				0,
				1,
				"M2/M1");
  if(check) {*iret = i; return 1;}
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Jktebop->M2_M1_vary = 1;
      i++;
    } else {
      Jktebop->M2_M1_vary = 0;
    }
  }


  if(i >= argc)
    return 1;
  if(strcmp(argv[i],"J2/J1")) {
    *iret = i; return 1;
  } else {
    i++;
  }
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
				VARTOOLS_TYPE_DOUBLE, 
				(void *) (&Jktebop->J2_J1),
				0,
				1,
				"J2/J1");
  if(check) {*iret = i; return 1;}
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Jktebop->J2_J1_vary = 1;
      i++;
    } else {
      Jktebop->J2_J1_vary = 0;
    }
  }

  if(i >= argc) {*iret = i; return 1;}
  if(!strcmp(argv[i],"i")) {
    Jktebop->use_i_or_b = JKTEBOP_USEI;
    i++;
    check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					   VARTOOLS_TYPE_DOUBLE, 
					   (void *) (&Jktebop->incl),
					   0,
					   1,
					   "INCLINATION");
    if(check) {*iret = i; return 1;}
    if(i < argc) {
      if(!strcmp(argv[i],"vary")) {
	Jktebop->incl_vary = 1;
	i++;
      } else {
	Jktebop->incl_vary = 0;
      }
    }
    /* We did not scan the bimpact vector from the command-line, but we
       still need to allocate memory for it so that the calculated values
       of bimpact can be stored. To do this register it is a "computed"
       vector */
    VARTOOLS_RegisterDataVector(p, c, (void *) (&Jktebop->bimpact),
				VARTOOLS_TYPE_DOUBLE, 0,
				VARTOOLS_SOURCE_COMPUTED,
				1, "bimpact");

  } else if(!strcmp(argv[i],"bimpact")) {
    Jktebop->use_i_or_b = JKTEBOP_USEBIMPACT;
    i++;
    check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					   VARTOOLS_TYPE_DOUBLE, 
					   (void *) (&Jktebop->bimpact),
					   0,
					   1,
					   "BIMPACT");
    if(check) {*iret = i; return 1;}
    if(i < argc) {
      if(!strcmp(argv[i],"vary")) {
	Jktebop->bimpact_vary = 1;
	i++;
      } else {
	Jktebop->bimpact_vary = 0;
      }
    }
    /* We did not scan the incl vector from the command-line, but we
       still need to allocate memory for it so that the calculated values
       of incl can be stored. To do this register it is a "computed"
       vector */
    VARTOOLS_RegisterDataVector(p, c, (void *) (&Jktebop->incl),
				VARTOOLS_TYPE_DOUBLE, 0,
				VARTOOLS_SOURCE_COMPUTED,
				1, "INCLINATION");
 } else {*iret = i; return 1;}

  if(i >= argc)
    return 1;
  if(strcmp(argv[i],"esinomega")) {
    *iret = i; return 1;
  } else {
    i++;
  }
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					 VARTOOLS_TYPE_DOUBLE, 
					 (void *) (&Jktebop->esinomega),
					 0,
					 1,
					 "ESINOMEGA");
  if(check) {*iret = i; return 1;}
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Jktebop->esinomega_vary = 1;
      i++;
    } else {
      Jktebop->esinomega_vary = 0;
    }
  }


  if(i >= argc)
    return 1;
  if(strcmp(argv[i],"ecosomega")) {
    *iret = i; return 1;
  } else {
    i++;
  }
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					 VARTOOLS_TYPE_DOUBLE, 
					 (void *) (&Jktebop->ecosomega),
					 0,
					 1,
					 "ECOSOMEGA");
  if(check) {*iret = i; return 1;}
  if(i < argc) {
    if(!strcmp(argv[i],"vary")) {
      Jktebop->ecosomega_vary = 1;
      i++;
    } else {
      Jktebop->ecosomega_vary = 0;
    }
  }

  if(i < argc) {
    if(!strcmp(argv[i],"LD1")) {
      /* Manually parse the LD law to use for star 1 */
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"linear")) {
	  Jktebop->LD1law = JKTEBOP_LDLAW_LINEAR;
	  Nldcoeff = 1;
	}
	else if(!strcmp(argv[i],"quad")) {
	  Jktebop->LD1law = JKTEBOP_LDLAW_QUAD;
	  Nldcoeff = 2;
	}
	else if(!strcmp(argv[i],"log")) {
	  Jktebop->LD1law = JKTEBOP_LDLAW_LOG;
	  Nldcoeff = 2;
	}
	else if(!strcmp(argv[i],"sqrt")) {
	  Jktebop->LD1law = JKTEBOP_LDLAW_SQRT;
	  Nldcoeff = 2;
	}
	else {
	  *iret = i; return 1;
	}
      }
      else {
	*iret = i; return 1;
      }

      i++;
      /* Now parse in how the LD coefficients are to be set */
      check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					     VARTOOLS_TYPE_DOUBLE, 
					     (void *) (&Jktebop->LD1_coeffs),
					     Nldcoeff,
					     1,
					     "STAR1_LD_COEFF");
      if(check) {*iret = i; return 1;}
      if(i < argc) {
	if(!strcmp(argv[i],"vary")) {
	  Jktebop->LD1_vary = 1;
	  i++;
	} else {
	  Jktebop->LD1_vary = 0;
	}
      }
      else {
	Jktebop->LD1_vary = 0;
      }
    } else {
      *iret = i; return 1;
    }
  } else {
    *iret = i; return 1;
  }


  if(i < argc) {
    if(!strcmp(argv[i],"LD2")) {
      /* Manually parse the LD law to use for star 2 */
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"linear")) {
	  Jktebop->LD2law = JKTEBOP_LDLAW_LINEAR;
	  Nldcoeff = 1;
	}
	else if(!strcmp(argv[i],"quad")) {
	  Jktebop->LD2law = JKTEBOP_LDLAW_QUAD;
	  Nldcoeff = 2;
	}
	else if(!strcmp(argv[i],"log")) {
	  Jktebop->LD2law = JKTEBOP_LDLAW_LOG;
	  Nldcoeff = 2;
	}
	else if(!strcmp(argv[i],"sqrt")) {
	  Jktebop->LD2law = JKTEBOP_LDLAW_SQRT;
	  Nldcoeff = 2;
	}
	else if(!strcmp(argv[i],"lockLD1")) {
	  Jktebop->LD2law = JKTEBOP_LDLAW_LOCKLD1;
	  Nldcoeff = 0;
	}
	else {
	  *iret = i; return 1;
	}
      }
      else {
	*iret = i; return 1;
      }

      i++;

      if(Nldcoeff > 0) {
	/* Now parse in how the LD coefficients are to be set */
	check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					       VARTOOLS_TYPE_DOUBLE, 
					       (void *) (&Jktebop->LD2_coeffs),
					       Nldcoeff,
					       1,
					       "STAR2_LD_COEFF");
	if(check) {*iret = i; return 1;}
	if(i < argc) {
	  if(!strcmp(argv[i],"vary")) {
	    Jktebop->LD2_vary = 1;
	    i++;
	  } else {
	    Jktebop->LD2_vary = 0;
	  }
	}
	else {
	  Jktebop->LD2_vary = 0;
	}
      }
    } else {
      *iret = i; return 1;
    }
  } else {
    *iret = i; return 1;
  }
  
  if(i < argc) {
    if(!strcmp(argv[i],"gravdark1")) {
      i++;
      check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					     VARTOOLS_TYPE_DOUBLE, 
					     (void *) (&Jktebop->gravdark1),
					     0,
					     1,
					     "GRAVDARK1");
      if(check) {*iret = i; return 1;}
      if(i < argc) {
	if(!strcmp(argv[i],"vary")) {
	  Jktebop->gravdark1_vary = 1;
	  i++;
	} else {
	  Jktebop->gravdark1_vary = 0;
	}
      }
      else
	Jktebop->gravdark1_vary = 0;
    } else {
      /* The user did not give gravdark1 on the command line
	 so we will instead fix it to a default value. We use
	 the VARTOOLS_RegisterDataVector function to register
	 the gravdark1 vector, and we set it to the fixed default
	 gravity darkening coefficient value */
      gravdark = JKTEBOP_DEFAULT_GRAVDARK;

      VARTOOLS_RegisterDataVector(p, c, (void *) (&Jktebop->gravdark1),
				  VARTOOLS_TYPE_DOUBLE,
				  0, VARTOOLS_SOURCE_FIXED, 0, NULL,
				  (void *) (&gravdark));
      Jktebop->gravdark1_vary = 0;
    }
  } else {
      gravdark = JKTEBOP_DEFAULT_GRAVDARK;

      VARTOOLS_RegisterDataVector(p, c, (void *) (&Jktebop->gravdark1),
				  VARTOOLS_TYPE_DOUBLE,
				  0, VARTOOLS_SOURCE_FIXED, 0, NULL,
				  (void *) (&gravdark));
      Jktebop->gravdark1_vary = 0;
  }


  if(i < argc) {
    if(!strcmp(argv[i],"gravdark2")) {
      i++;
      check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					     VARTOOLS_TYPE_DOUBLE, 
					     (void *) (&Jktebop->gravdark2),
					     0,
					     1,
					     "GRAVDARK2");
      if(check) {*iret = i; return 1;}
      if(i < argc) {
	if(!strcmp(argv[i],"vary")) {
	  Jktebop->gravdark2_vary = 1;
	  i++;
	} else {
	  Jktebop->gravdark2_vary = 0;
	}
      }
      else
	Jktebop->gravdark2_vary = 0;
    } else {
      /* The user did not give gravdark2 on the command line
	 so we will instead fix it to a default value. We use
	 the VARTOOLS_RegisterDataVector function to register
	 the gravdark1 vector, and we set it to the fixed default
	 gravity darkening coefficient value */
      gravdark = JKTEBOP_DEFAULT_GRAVDARK;

      VARTOOLS_RegisterDataVector(p, c, (void *) (&Jktebop->gravdark2),
				  VARTOOLS_TYPE_DOUBLE,
				  0, VARTOOLS_SOURCE_FIXED, 0, NULL,
				  (void *) (&gravdark));
      Jktebop->gravdark2_vary = 0;
    }
  } else {
      gravdark = JKTEBOP_DEFAULT_GRAVDARK;

      VARTOOLS_RegisterDataVector(p, c, (void *) (&Jktebop->gravdark2),
				  VARTOOLS_TYPE_DOUBLE,
				  0, VARTOOLS_SOURCE_FIXED, 0, NULL,
				  (void *) (&gravdark));
      Jktebop->gravdark2_vary = 0;
  }


  if(i < argc) {
    if(!strcmp(argv[i],"reflection1")) {
      i++;
      check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					     VARTOOLS_TYPE_DOUBLE, 
					     (void *) (&Jktebop->reflection1),
					     0,
					     1,
					     "REFLECTION1");
      if(check) {*iret = i; return 1;}
      if(i < argc) {
	if(!strcmp(argv[i],"vary")) {
	  Jktebop->reflection1_vary = 1;
	  i++;
	} else {
	  Jktebop->reflection1_vary = 0;
	}
      }
      else
	Jktebop->reflection1_vary = 0;
    } else {
      /* The user did not give reflection1 on the command line
	 so we will instead fix it to a default value. We use
	 the VARTOOLS_RegisterDataVector function to register
	 the gravdark1 vector, and we set it to the fixed default
	 gravity darkening coefficient value */
      reflection = JKTEBOP_DEFAULT_REFLECTION;

      VARTOOLS_RegisterDataVector(p, c, (void *) (&Jktebop->reflection1),
				  VARTOOLS_TYPE_DOUBLE,
				  0, VARTOOLS_SOURCE_FIXED, 0, NULL,
				  (void *) (&reflection));
      Jktebop->reflection1_vary = 0;
    }
  } else {
      reflection = JKTEBOP_DEFAULT_REFLECTION;

      VARTOOLS_RegisterDataVector(p, c, (void *) (&Jktebop->reflection1),
				  VARTOOLS_TYPE_DOUBLE,
				  0, VARTOOLS_SOURCE_FIXED, 0, NULL,
				  (void *) (&reflection));
      Jktebop->reflection1_vary = 0;
  }


  if(i < argc) {
    if(!strcmp(argv[i],"reflection2")) {
      i++;
      check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					     VARTOOLS_TYPE_DOUBLE, 
					     (void *) (&Jktebop->reflection2),
					     0,
					     1,
					     "REFLECTION2");
      if(check) {*iret = i; return 1;}
      if(i < argc) {
	if(!strcmp(argv[i],"vary")) {
	  Jktebop->reflection2_vary = 1;
	  i++;
	} else {
	  Jktebop->reflection2_vary = 0;
	}
      }
      else
	Jktebop->reflection2_vary = 0;
    } else {
      /* The user did not give reflection2 on the command line
	 so we will instead fix it to a default value. We use
	 the VARTOOLS_RegisterDataVector function to register
	 the gravdark1 vector, and we set it to the fixed default
	 gravity darkening coefficient value */
      reflection = JKTEBOP_DEFAULT_REFLECTION;

      VARTOOLS_RegisterDataVector(p, c, (void *) (&Jktebop->reflection2),
				  VARTOOLS_TYPE_DOUBLE,
				  0, VARTOOLS_SOURCE_FIXED, 0, NULL,
				  (void *) (&reflection));
      Jktebop->reflection2_vary = 0;
    }
  } else {
      reflection = JKTEBOP_DEFAULT_REFLECTION;

      VARTOOLS_RegisterDataVector(p, c, (void *) (&Jktebop->reflection2),
				  VARTOOLS_TYPE_DOUBLE,
				  0, VARTOOLS_SOURCE_FIXED, 0, NULL,
				  (void *) (&reflection));
      Jktebop->reflection2_vary = 0;
  }


  if(i < argc) {
    if(!strcmp(argv[i],"L3")) {
      i++;
      check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					     VARTOOLS_TYPE_DOUBLE, 
					     (void *) (&Jktebop->L3),
					     0,
					     1,
					     "L3");
      if(check) {*iret = i; return 1;}
      if(i < argc) {
	if(!strcmp(argv[i],"vary")) {
	  Jktebop->L3_vary = 1;
	  i++;
	} else {
	  Jktebop->L3_vary = 0;
	}
      }
      else
	Jktebop->L3_vary = 0;
    } else {
      /* The user did not give L3 on the command line
	 so we will instead fix it to a default value. We use
	 the VARTOOLS_RegisterDataVector function to register
	 the gravdark1 vector, and we set it to the fixed default
	 gravity darkening coefficient value */
      L3 = JKTEBOP_DEFAULT_L3;

      VARTOOLS_RegisterDataVector(p, c, (void *) (&Jktebop->L3),
				  VARTOOLS_TYPE_DOUBLE,
				  0, VARTOOLS_SOURCE_FIXED, 0, NULL,
				  (void *) (&L3));
      Jktebop->L3_vary = 0;
    }
  } else {
      L3 = JKTEBOP_DEFAULT_L3;

      VARTOOLS_RegisterDataVector(p, c, (void *) (&Jktebop->L3),
				  VARTOOLS_TYPE_DOUBLE,
				  0, VARTOOLS_SOURCE_FIXED, 0, NULL,
				  (void *) (&L3));
      Jktebop->L3_vary = 0;
  }


  if(i < argc) {
    if(!strcmp(argv[i],"tidallag")) {
      i++;
      check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					     VARTOOLS_TYPE_DOUBLE, 
					     (void *) (&Jktebop->tidallag),
					     0,
					     1,
					     "TIDALLAG");
      if(check) {*iret = i; return 1;}
      if(i < argc) {
	if(!strcmp(argv[i],"vary")) {
	  Jktebop->tidallag_vary = 1;
	  i++;
	} else {
	  Jktebop->tidallag_vary = 0;
	}
      }
      else
	Jktebop->tidallag_vary = 0;
    } else {
      /* The user did not give tidallag on the command line
	 so we will instead fix it to a default value. We use
	 the VARTOOLS_RegisterDataVector function to register
	 the tidallag vector, and we set it to the fixed default
	 gravity darkening coefficient value */
      tidallag = JKTEBOP_DEFAULT_TIDALLAG;

      VARTOOLS_RegisterDataVector(p, c, (void *) (&Jktebop->tidallag),
				  VARTOOLS_TYPE_DOUBLE,
				  0, VARTOOLS_SOURCE_FIXED, 0, NULL,
				  (void *) (&tidallag));
      Jktebop->tidallag_vary = 0;
    }
  } else {
      L3 = JKTEBOP_DEFAULT_TIDALLAG;

      VARTOOLS_RegisterDataVector(p, c, (void *) (&Jktebop->tidallag),
				  VARTOOLS_TYPE_DOUBLE,
				  0, VARTOOLS_SOURCE_FIXED, 0, NULL,
				  (void *) (&tidallag));
      Jktebop->tidallag_vary = 0;
  }
  
  /* Check if the user has given the "correctlc" keyword */
  if(i < argc) {
    if(!strcmp(argv[i],"correctlc")) {
      Jktebop->correctlc = 1;
      i++;
    }
    else
      Jktebop->correctlc = 0;
  }
  else
    Jktebop->correctlc = 0;

  /* Check if the user has given the "omodel" keyword */
  if(i < argc) {
    if(!strcmp(argv[i],"omodel")) {
      Jktebop->omodel = 1;
      i++;
      if(i < argc) {
	/* Copy the outdir name, given at argv[i] to the outdir variable */
	sprintf(Jktebop->outdir,"%s",argv[i]);
	i++;
      }
      else {
	/* The user did not give the name of the output directory, return an
	   error */
	*iret = i; return 1;
      }
      /* Check if the user gave the "format" keyword */
      if(i < argc) {
	if(!strcmp(argv[i],"format")) {
	  i++;
	  if(i < argc) {
	    /* Copy the format string given on the command-line to the format
	       variable */
	    sprintf(Jktebop->format,"%s",argv[i]);
	    i++;
	  } else {
	    /* The user did not give the format string on the command-line,
	       return an error */
	    *iret = i; return 1;
	  }
	} else {
	  /* If the user didn't give the "format" keyword, make
	     Jktebop->format an empty string. */
	  Jktebop->format[0] = '\0';
	}
      } else {
	Jktebop->format[0] = '\0';
      }
    } else {
      Jktebop->omodel = 0;
    }
  } else {
    Jktebop->omodel = 0;
  }

  /* Check if the user has given the "ocurve" keyword */
  if(i < argc) {
    if(!strcmp(argv[i],"ocurve")) {
      Jktebop->ocurve = 1;
      i++;
      if(i >= argc) {
	*iret = i; return 1;
      }
      if(!strcmp(argv[i],"jd")) {
	Jktebop->ocurvetype = JKTEBOP_OUTCURVE_TYPE_JD;
      }
      else if(!strcmp(argv[i],"phase")) {
	Jktebop->ocurvetype = JKTEBOP_OUTCURVE_TYPE_PHASE;
      }
      else {
	*iret = i; return 1;
      }
      
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"step")) {
	  i++;
	  if(i >= argc) {
	    *iret = i; return 1;
	  }
	  Jktebop->ocurvestep = atof(argv[i]);
	}
	else {
	  Jktebop->ocurvestep = JKTEBOP_OUTCURVE_STEP_DEFAULT;
	  i--;
	}
      }
      else {
	Jktebop->ocurvestep = JKTEBOP_OUTCURVE_STEP_DEFAULT;
	i--;
      }
      
      i++;
      if(i >= argc) {
	*iret = i; return 1;
      }
      if(!strcmp(argv[i],"outdir")) {
	i++;
	if(i >= argc) {
	  *iret = i; return 1;
	}
	sprintf(Jktebop->ocurve_outdir,"%s",argv[i]);

	i++;
	if(i < argc) {
	  if(!strcmp(argv[i],"format")) {
	    i++;
	    if(i >= argc) {
	      *iret = i; return 1;
	    }
	    sprintf(Jktebop->ocurve_outdir_format,"%s",argv[i]);
	    i++;
	  }
	  else {
	    sprintf(Jktebop->ocurve_outdir_format,"");
	  }
	} else {
	  sprintf(Jktebop->ocurve_outdir_format,"");
	}
      }
      else {
	*iret = i; return 1;
      }
      
    } else {
      Jktebop->ocurve = 0;
    }
  } else {
    Jktebop->ocurve = 0;
  }

  /* If we are fitting the model, register chi2 and Ndof as computed vectors */
  if(Jktebop->injectorfit == JKTEBOP_FIT) {
    VARTOOLS_RegisterDataVector(p, c, (void *) (&Jktebop->chi2val),
				VARTOOLS_TYPE_DOUBLE,
				0, VARTOOLS_SOURCE_COMPUTED, 1, "CHI2");
    VARTOOLS_RegisterDataVector(p, c, (void *) (&Jktebop->Ndof),
				VARTOOLS_TYPE_INT,
				0, VARTOOLS_SOURCE_COMPUTED, 1, "NDOF");
  }
    
  /* Update the iret variable to point to the next command-line term */
  *iret = i;
  

  /* Check = 1 if there was an error in parsing the line, return
     a positive value if there is to tell vartools to print out the 
     syntax and quit */
  if(check)
    return 1;
  else
    return 0;
}

void jktebop_ShowSyntax(FILE *outfile)
/* Write out the expected syntax for this command to the file outfile */
{
  fprintf(outfile,
   "-jktebop < \"inject\" | \"fit\" >\n"
   "   <\"Period\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]>\n"
   "   <\"T0\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]>\n"
   "   <\"r1+r2\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]>\n"
   "   <\"r2/r1\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]>\n"
   "   <\"M2/M1\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]>\n"
   "   <\"J2/J1\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]>\n"
   "   <<\"i\" | \"bimpact\"> <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]>\n"
   "   <\"esinomega\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]>\n"
   "   <\"ecosomega\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]>\n"
   "   <\"LD1\" <\"linear\" | \"quad\" | \"log\" | \"sqrt\">\n"
   "       <\"fix\" value1 [value2] | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]>\n"
   "   <\"LD2\" <\"lockLD1\" | \"linear\" | \"quad\" | \"log\" | \"sqrt\">\n"
   "       [<\"fix\" value [value2] | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]>\n"
   "   [\"gravdark1\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]]\n"
   "   [\"gravdark2\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]]\n"
   "   [\"reflection1\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]]\n"
   "   [\"reflection2\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]]\n"
   "   [\"L3\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]]\n"
   "   [\"tidallag\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>\n"
   "       [\"vary\"]]\n"
   "   [\"correctlc\"]\n"
   "   [\"omodel\" <outdir [\"format\" fmt]>]\n"
   "   [\"ocurve\" <\"jd\" | \"phase\"> [\"step\" stepsize]\n"
   "       <\"outdir\" outdir [ \"format\" fmt]>]\n"
	  );
}

void jktebop_ShowHelp(FILE *outfile)
/* Output the help information for this command */
{
  fprintf(outfile,
"Fit or inject a JKTEBOP detached eclipsing binary model [in]to the light curves. The user must give either the \"inject\" keyword or the \"fit\" keyword to indicate whether the model should be injected into the light curve, or fitted to it. After that the user specifies how to initialize the parameters used by this model. The parameters include:\n"
"\t\"Period\" - the orbital period in days.\n"
"\t\"T0\" - central time of a primary eclipse.\n"
"\t\"r1+r2\" - the sum of the radii of the stars divided by the semi-major axis.\n"
"\t\"r2/r1\" - the ratio of the stellar radii.\n"
"\t\"M2/M1\" - the mass ratio.\n"
"\teither:\n"
"\t\t\"i\" - the inclination angle in degrees (90 degrees is an edge-on orbit)\n"
"\tor:\n"
"\t\t\"bimpact\" - the impact parameter of the primary eclipse. This is the\n"
"\t\t            projected distance between the centers of the two stars\n"
"\t\t            at the time T0, divided by the sum of the radii of the stars.\n"
"\t\t            It goes from 0 for a central eclipse, to 1 for a just\n"
"\t\t            grazing eclipse.\n"
"\t\"esinomega\" - the eccentricity times the sin of the argument of pericenter.\n"
"\t\"ecosomega\" - the eccentricity times the cos of the argument of pericenter.\n"
"\t\"LD[1-2]\" - the limb darkening coefficients for each star. After\n"
"\t\tgiving the keyword the user must specify the LD law to use\n"
"\t\twhich can be \"linear\", \"quad\", \"log\" or \"sqrt\". For\n"
"\t\t\"LD2\" the user can also specify \"lockLD1\" which forces\n"
"\t\tthe secondary star to have the same LD coefficients as the\n"
"\t\tprimary star. If a linear law is used, the user must specify\n"
"\t\tone coefficient, for the other laws the user specifies two\n"
"\t\tcoefficients. If \"lockLD1\" is used for \"LD2\", the user\n"
"\t\tdoes not provide any coefficients for \"LD2\".\n"
"\t\"gravdark[1-2]\" - Optional gravity darkening coefficients for each star.\n"
"\t\tIf not given on the command line, a value of 1.0 is assumed for each star.\n"
"\t\tSuggested values are 1.0 for stars with radiative envelopes, and 0.3 for\n"
"\t\tstars with convective envelopes.\n"
"\t\"reflection[1-2]\" - Optional reflection effect coefficients for each star.\n"
"\t\tIf not given on the command line they will be calculated. If a value <= 0\n"
"\t\tis given, they will also be calculated.\n"
"\t\"L3\" - Optional third light parameter. This is set to 0 by default.\n"
"\t\"tidallag\" - Optional tidal lag angle in degrees. This is set to 0\n"
"\t\tby default.\n"
"For each parameter the user specifies the source for the initial\n"
"parameter value. This can either be the keyword \"fix\" followed by\n"
"the value to use for all light curves, \"list\" to read the parameter\n"
"from the light curve list (use the \"column\" keyword followed by a\n"
"number to indicate the column number from the list to use, otherwise\n"
"the next column in the list will be assumed), \"fixcolumn\"\n"
"followed by the name or number of an output column from a previously\n"
"executed command, or \"expr\" followed by an analytic expression. For\n"
"each parameter the user can give the \"vary\"\n"
"keyword which indicates that the parameter is to be varied in a\n"
"fit. If not given the parameter will be fixed in the fit. Finally the\n"
"user may use the \"correctlc\" keyword to subtract the best-fit model\n"
"from the light curve, and the \"omodel\" keyword to output the model\n"
"for each light curve. If the \"omodel\" keyword is given, the user\n"
"should specify the directory to output the model light curves to. By\n"
"default the output light curve will have the name\n"
"outdir/BASELC_NAME.jktebop where BASELC_NAME is the basename of the\n"
"input light curve. The user may however give the \"format\" keyword\n"
"followed by a format string to specify arbitrary filenames. The\n"
"syntax is the same as for the \"-o\" VARTOOLS command. The user may\n"
"also output a model curve file sampled at uniform spacing by giving\n"
"the \"ocurve\" keyword. This must be followed by the keyword \"jd\"\n"
"or the keyword \"phase\" to indicate whether the curve is output with\n"
"JD as the independent variable, or if it is phase. One can optional\n"
"specify the step-size of the curve by giving the \"step\" keyword followed\n"
"by the value. The user must then give the \"outdir\" keyword and the\n"
"directory to output the curve files to. The \"format\" keyword may be\n"
"used to specify arbitrary filenames. The default output will be\n"
"outdir/BASELC_NAME.jktebopcurve. If you use this\n"
"routine the following references should be cited: Southworth et\n"
"al. (2004MNRAS.351.1277S), Popper & Etzel (1981AJ.....86..102), Etzel\n"
"(1981psbs.conf..111E), and/or Nelson & Davis (1972ApJ...174..617N)\n\n"
	  );
}

void jktebop_RunCommand(ProgramData *p, void *userdata, int lc_name_num, int lc_num)
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
  char *lcname;
  _Jktebop *Jktebop;

  char lcoutname[MAXLEN];

  void dojktebopinject(int, double *, double *, double *, _Jktebop*, int, char*);
  void dojktebopfit(int, double *, double *, double *, _Jktebop*, int, char*);


  /* Cast the userdata pointer given by VARTOOLS into a _Jktebop type
     struct. Then we can access its contents. */
  Jktebop = (_Jktebop *)userdata;

  /* The number of points in the light curve, and vectors storing the
     time, magnitudes, and errors are accessed as illustrated
     below. Note that t, mag, and err are vectors, t[0] is the time of
     the first observation in the light curve, t[NJD-1] is the time of
     the last observation, etc. */
  NJD=p->NJD[lc_num];
  t = p->t[lc_num];
  mag = p->mag[lc_num];
  err = p->sig[lc_num];

  /* Access the light curve name as follows: */
  lcname = p->lcnames[lc_name_num];

  /* This command allows the user to optionally output model light curves.
     It uses a "standard" VARTOOLS syntax for specifying the output name,
     which is to provide a directory for the output, and then by default
     set the name to something like outdir/BASELC_NAME.suffix unless the
     user gives a format keyword, in which case the user can provide an
     arbitrary name. Any command using this convention can call the
     VARTOOLS_GetOutputFilename function to get the output filename. */
  if(Jktebop->omodel)
    VARTOOLS_GetOutputFilename(lcoutname,lcname,Jktebop->outdir,"jktebop",
			       Jktebop->format, lc_name_num);
  else 
    lcoutname[0] = '\0';

  /* Perform the actual routine, we use separate functions for injecting
     and for fitting the light curve. Here we pass the Jktebop struct which
     contains all the parameter data, and the lc_num index which we need
     to reference the correct components of the parameter data arrays */
  if(Jktebop->injectorfit == JKTEBOP_INJECT) {
    dojktebopinject(NJD, t, mag, err, Jktebop, lc_num, lcoutname);
  } else {
    dojktebopfit(NJD, t, mag, err, Jktebop, lc_num, lcoutname);
  }
  
}

/***************************************************************************
 Below are functions that are specific to this command. These are not
 required by a user library, but in this case are called through the
 jktebop_RunCommand function and perform the actual processing routine.
***************************************************************************/

void dojktebopinject(int NJD, double *t, double *mag, double *err, _Jktebop *Jktebop, int lc_num, char *lcoutname)
{
  /* Inject a JKTEBOP model into a light curve. Optionally output the
     model.  If the user requested that the model be subtracted with
     the correctlc keyword we will do it, but it is not obvious why
     they would want this. */
  double *magsim;
  double Period, T0, r1r2, r2_r1, M2_M1, J2_J1, incl, bimpact, esinomega,
    ecosomega, LD1a, LD1b, LD2a, LD2b, gravdark1, gravdark2, reflection1,
    reflection2, L3, tidallag;
  int LD1law, LD2law, use_i_or_b;
  int i;

  FILE *outfile;

  void simulateEBOPEBlc(int N, double *JD, double *magout, double T0, double period, double esinomega, double ecosomega, int use_i_or_b, double *incl, double *bimpact, double r1r2, double r2_r1, double q, double J2J1, double L3, double grav1, double grav2, double reflection1, double reflection2, double quadA1, double quadB1, double quadA2, double quadB2, double tidallag, int ldtype1, int ldtype2);

  /* Get the parameters from the Jktebop data structure */
  Period = Jktebop->Period[lc_num];
  T0 = Jktebop->T0[lc_num];
  r1r2 = Jktebop->r1r2[lc_num];
  r2_r1 = Jktebop->r2_r1[lc_num];
  M2_M1 = Jktebop->M2_M1[lc_num];
  J2_J1 = Jktebop->J2_J1[lc_num];
  use_i_or_b = Jktebop->use_i_or_b;
  if(use_i_or_b == JKTEBOP_USEI) {
    incl = Jktebop->incl[lc_num];
    bimpact = 0.;
  } else {
    incl = 0.;
    bimpact = Jktebop->bimpact[lc_num];
  }
  esinomega = Jktebop->esinomega[lc_num];
  ecosomega = Jktebop->ecosomega[lc_num];
  LD1law = Jktebop->LD1law;
  switch(LD1law) {
  case JKTEBOP_LDLAW_LINEAR:
    LD1a = Jktebop->LD1_coeffs[lc_num][0];
    LD1b = 0.;
    break;
  default:
    LD1a = Jktebop->LD1_coeffs[lc_num][0];
    LD1b = Jktebop->LD1_coeffs[lc_num][1];
  }
  LD2law = Jktebop->LD2law;
  switch(LD2law) {
  case JKTEBOP_LDLAW_LOCKLD1:
    LD2a = LD1a;
    LD2b = LD1b;
    LD2law = LD1law;
    break;
  case JKTEBOP_LDLAW_LINEAR:
    LD2a = Jktebop->LD2_coeffs[lc_num][0];
    LD2b = 0.;
    break;
  default:
    LD2a = Jktebop->LD2_coeffs[lc_num][0];
    LD2b = Jktebop->LD2_coeffs[lc_num][1];
  }
  gravdark1 = Jktebop->gravdark1[lc_num];
  gravdark2 = Jktebop->gravdark2[lc_num];
  reflection1 = Jktebop->reflection1[lc_num];
  reflection2 = Jktebop->reflection2[lc_num];
  L3 = Jktebop->L3[lc_num];
  tidallag = Jktebop->tidallag[lc_num];
  
  /* Allocate memory to store the simulated light curve */
  if((magsim = (double *) malloc(NJD * sizeof(double))) == NULL)
    {
      fprintf(stderr,"Memory Allocation Error in function dojktebopinject\n");
      exit(ERR_MEMALLOC);
    }
  
  /* Simulate the light curve */
  simulateEBOPEBlc(NJD, t, magsim, T0, Period, esinomega, ecosomega, use_i_or_b, &incl, &bimpact, r1r2, r2_r1, M2_M1, J2_J1, L3, gravdark1, gravdark2, reflection1, reflection2, LD1a, LD1b, LD2a, LD2b, tidallag, LD1law, LD2law);

  /* Update the incl or bimpact value in the Jktebop structure, whichever was
     calculated */
  if(use_i_or_b == JKTEBOP_USEI) {
    Jktebop->bimpact[lc_num] = bimpact;
  } else {
    Jktebop->incl[lc_num] = incl;
  }

  /* Add the model to the light curve */
  for(i=0; i < NJD; i++) {
    mag[i] += magsim[i];
  }
  
  /* Output the model if asked to */
  if(Jktebop->omodel) {
    if((outfile = fopen(lcoutname,"w")) == NULL) {
      fprintf(stderr,"Cannot write to %s\n",lcoutname);
      exit(ERR_CANNOTWRITE);
    }
    fprintf(outfile,"#Time  Mag_lc   Err  Mag_model\n");
    for(i=0; i < NJD; i++) {
      if(!isnan(mag[i])) {
	fprintf(outfile,"%.17g %.17g %.17g %.17g\n", t[i], mag[i], err[i], magsim[i]);
      }
    }
    fclose(outfile);
  }
  
  /* Output the curve if asked to */
  /******** TBD --- IMPLEMENT THIS ******/
  
  /* Subtract the model if asked to */
  if(Jktebop->correctlc) {
    for(i=0; i < NJD; i++) {
      if(!isnan(mag[i]))
	mag[i] -= magsim[i];
    }
  }

  /* Free the memory allocated for the simulated lc */
  free(magsim);

}


void dojktebopfit(int NJD, double *t, double *mag, double *err, _Jktebop *Jktebop, int lc_num, char *lcoutname)
{
  int amoeba_val;
  double best_chi2;
  int Ndof, i;
  double *magsim;
  FILE *outfile;

  int fitjktebop_amoeba(int N, double *t, double *mag, double *err, _Jktebop *Jktebop, int lc_num, double *best_chi2, int *Ndof, double *magsim);

  /* Allocate memory to store the simulated light curve */
  if((magsim = (double *) malloc(NJD * sizeof(double))) == NULL)
    {
      fprintf(stderr,"Memory Allocation Error in function dojktebopinject\n");
      exit(ERR_MEMALLOC);
    }
  
  /* Run the amoeba fit */
  amoeba_val = fitjktebop_amoeba(NJD, t, mag, err, Jktebop, lc_num, &best_chi2, &Ndof, magsim);
  
  /* Store the chi2 value and the number of degrees of freedom */
  Jktebop->chi2val[lc_num] = best_chi2;
  Jktebop->Ndof[lc_num] = Ndof;

  if(!amoeba_val) {
    /* The fit worked, output the model if requested, and correct the lc if
       requested */
    if(Jktebop->omodel) {
      if((outfile = fopen(lcoutname,"w")) == NULL) {
	fprintf(stderr,"Cannot write to %s\n",lcoutname);
	exit(ERR_CANNOTWRITE);
      }
      fprintf(outfile,"#Time  Mag_lc   Err  Mag_model\n");
      for(i=0; i < NJD; i++) {
	if(!isnan(mag[i]))
	  fprintf(outfile,"%.17g %.17g %.17g %.17g\n", t[i], mag[i], err[i], magsim[i]);
      }
      fclose(outfile);
    }
    
    /* Output the curve if asked to */
    /******** TBD --- IMPLEMENT THIS ******/

    /* Subtract the model if asked to */
    if(Jktebop->correctlc) {
      for(i=0; i < NJD; i++) {
	if(!isnan(mag[i]))
	  mag[i] -= magsim[i];
      }
    }
  }

  free(magsim);
  return;
}

#define JKTEBOP_AMOEBA_CONVERGENCELIMIT 0.0001

int fitjktebop_amoeba(int N, double *t, double *mag, double *err, _Jktebop *Jktebop, int lc_num, double *best_chi2, int *Ndof, double *magsim)
{
  /* This will carry-out the amoeba downhill-simplex fit, on output the varied
     parameters in Jktebop will take on their best-fit values 
     Most of the code below is to initialize the simplex intelligently.
  */

  int i, *ia, Nparameters = 0, Ntovary = 0, N1, N2, nfunk, amoeba_val = 0;
  int j;
  double **p, *chi2vals = NULL, T1, T2, ph1, ph2, delta, sini, ftol;
  _Jktebop_fitstruct fitstruct;

  double chisqjktebop(double *p, int ma, int N, double *t, double *mag, double *sig, void *userparams);

  /* Set up the auxiliary parameter structure */
  fitstruct.use_i_or_b = Jktebop->use_i_or_b;
  fitstruct.is_vary_P = Jktebop->Period_vary;
  fitstruct.is_vary_T0 = Jktebop->T0_vary;
  fitstruct.LD1law = Jktebop->LD1law;
  fitstruct.LD2law = Jktebop->LD2law;
  fitstruct.magsim = magsim;

  /* Initialize the simplex */

  /* We will use the VARTOOLS_incrementparameters_foramoeba
     function, included in the vartools library, to add parameters
     to the simplex. 
     
     We go through each of the parameters, determine the step-size in that
     parameter for the initial simplex, and then add the parameter to the
     simplex
  */

  if(Jktebop->Period_vary && Jktebop->T0_vary) {
    /* If both the period and T0 are to be varied, then in practice
       vary T1 and T2 associated with eclipse numbers N1 and N2.
    */
    /* Find the first and last eclipse minima within the observational
     window*/
    ph1 = (t[0] - Jktebop->T0[lc_num])/Jktebop->Period[lc_num];
    ph1 = ph1 - floor(ph1);
    T1 = t[0] + (1. - ph1)*Jktebop->Period[lc_num];
    ph2 = (t[N-1] - Jktebop->T0[lc_num])/Jktebop->Period[lc_num];
    ph2 = ph2 - floor(ph2);
    T2 = t[N-1] - ph2*Jktebop->Period[lc_num];
    N1 = 0;
    /* N2 is the number of eclipses elapsed between T2 and T1 */
    N2 = (int) lround((T2 - T1)/Jktebop->Period[lc_num]);
    if(N2 > N1) {
      fitstruct.N1 = N1;
      fitstruct.N2 = N2;

      /* Allow the times to shift by a fraction of the eclipse duration */
      delta = 0.05*Jktebop->r1r2[lc_num]*(T2 - T1)/M_PI/((double) (N2 - N1));
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
						1, T1, delta);
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
						1, T2, delta);
    } else {
      /* We've only observed one eclipse, it doesn't make sense to vary the
	 period. We will overrule the wishes of the foolish user, and instead
	 only vary T0. */
      delta = 0.05*Jktebop->r1r2[lc_num]*Jktebop->Period[lc_num]/M_PI;
      N2 = N1 + 1;
      T2 = T1 + Jktebop->Period[lc_num];
      fitstruct.N1 = N1;
      fitstruct.N2 = N2;
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
						1, T1, delta);
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
						0, T2, 0.);
    }
  } else {
    /* We may be varying either the Period or T0 directly.
       Determine the appropriate time-steps to use for the Period and for
       T0 */
    if(Jktebop->T0_vary) {
      delta = 0.05*Jktebop->r1r2[lc_num]*Jktebop->Period[lc_num]/M_PI;
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
						1, Jktebop->T0[lc_num], delta);
    } else {
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
						0, Jktebop->T0[lc_num], 0.);
    }
    if(Jktebop->Period_vary) {
      ph1 = (t[0] - Jktebop->T0[lc_num])/Jktebop->Period[lc_num];
      ph1 = ph1 - floor(ph1);
      T1 = t[0] + (1. - ph1)*Jktebop->Period[lc_num];
      ph2 = (t[N-1] - Jktebop->T0[lc_num])/Jktebop->Period[lc_num];
      ph2 = ph2 - floor(ph2);
      T2 = t[N-1] - ph2*Jktebop->Period[lc_num];
      N1 = 0;
      N2 = (int) lround((T2 - T1)/Jktebop->Period[lc_num]);
      if(N2 > N1) {
	delta = 0.1*Jktebop->r1r2[lc_num]*Jktebop->Period[lc_num]/M_PI/((double)(N2 - N1));
      } else {
	T1 = Jktebop->T0[lc_num];
	N2 = (int) lround((T2 - T1)/Jktebop->Period[lc_num]);
	if(N2 == N1) {
	  delta = 0.1*Jktebop->Period[lc_num];
	} else {
	  delta = 0.1*Jktebop->r1r2[lc_num]*Jktebop->Period[lc_num]/M_PI/fabs((double)(N2 - N1));
	}
      }
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
						1, Jktebop->Period[lc_num], 
						delta);
    } else {
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
						0, Jktebop->Period[lc_num], 
						0.);
    }
  }
  
  /* delta for (r1+r2)/a */
  if(Jktebop->r1r2[lc_num] == 0.0) {
    delta = 0.1;
  } 
  else if(Jktebop->r1r2[lc_num] < 1.0/1.1) {
    delta = 0.1*Jktebop->r1r2[lc_num];
  } else {
    delta = -0.1*Jktebop->r1r2[lc_num];
  }
  VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					    Jktebop->r1r2_vary,
					    Jktebop->r1r2[lc_num], delta);
					    
  /* delta for R2/R1 */
  if(Jktebop->r2_r1[lc_num] == 0.0)
    delta = 0.1;
  else
    delta = 0.1*Jktebop->r2_r1[lc_num];
  VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					    Jktebop->r2_r1_vary,
					    Jktebop->r2_r1[lc_num], delta);
  
  /* delta for M2/M1 */
  if(Jktebop->M2_M1[lc_num] == 0.0)
    delta = 0.1;
  else
    delta = 0.1*Jktebop->M2_M1[lc_num];
  VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					    Jktebop->M2_M1_vary,
					    Jktebop->M2_M1[lc_num], delta);

  /* delta for J2/J1 */
  if(Jktebop->J2_J1[lc_num] == 0.0)
    delta = 0.1;
  else
    delta = 0.1*Jktebop->J2_J1[lc_num];
  VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					    Jktebop->J2_J1_vary,
					    Jktebop->J2_J1[lc_num], delta);

  /* Inclination, or bimpact */
  if(Jktebop->use_i_or_b == JKTEBOP_USEI) {
    sini = sin(M_PI*Jktebop->incl[lc_num]/180.);
    if(sini < 1.0) {
      delta = 0.1*(1.0 - sini);
    } else {
      delta = -0.01;
    }
    delta = 180.0*asin(sini + delta)/M_PI - Jktebop->incl[lc_num];
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					      Jktebop->incl_vary,
					      Jktebop->incl[lc_num], delta);
  } else {
    if(Jktebop->bimpact[lc_num] < 0.9)
      delta = 0.1;
    else
      delta = -0.1;
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					      Jktebop->bimpact_vary,
					      Jktebop->bimpact[lc_num], delta);
  }

  /* esinomega */
  if(Jktebop->esinomega[lc_num]*Jktebop->esinomega[lc_num] +
     Jktebop->ecosomega[lc_num]*Jktebop->ecosomega[lc_num] < 0.25) {
    delta = 0.1;
  } else {
    if(Jktebop->esinomega[lc_num] > 0.1) {
      delta = -0.1;
    } else {
      delta = 0.1;
    }
  }
  VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					    Jktebop->esinomega_vary,
					    Jktebop->esinomega[lc_num], delta);

  /* ecosomega */
  if(Jktebop->esinomega[lc_num]*Jktebop->esinomega[lc_num] +
     Jktebop->ecosomega[lc_num]*Jktebop->ecosomega[lc_num] < 0.25) {
    delta = 0.1;
  } else {
    if(Jktebop->ecosomega[lc_num] > 0.1) {
      delta = -0.1;
    } else {
      delta = 0.1;
    }
  }
  VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					    Jktebop->ecosomega_vary,
					    Jktebop->ecosomega[lc_num], delta);
  
  /* Limb darkening */
  delta = 0.05;
  VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					    Jktebop->LD1_vary,
					    Jktebop->LD1_coeffs[lc_num][0], 
					    delta);
  if(Jktebop->LD1law == JKTEBOP_LDLAW_LINEAR){
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					      0, 0., 0.);
  }
  else {
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					      Jktebop->LD1_vary,
					      Jktebop->LD1_coeffs[lc_num][1], 
					      delta);
  }
  if(Jktebop->LD2law == JKTEBOP_LDLAW_LOCKLD1){
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					      0, 0., 0.);
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					      0, 0., 0.);
  } else {
    VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					      Jktebop->LD2_vary,
					      Jktebop->LD2_coeffs[lc_num][0], 
					      delta);
    if(Jktebop->LD2law == JKTEBOP_LDLAW_LINEAR){
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					      0, 0., 0.);
    }
    else {
      VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
						Jktebop->LD2_vary,
						Jktebop->LD2_coeffs[lc_num][1], 
						delta);
    }
  }

  /* Gravity darkening */
  delta = 0.05;
  VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					    Jktebop->gravdark1_vary,
					    Jktebop->gravdark1[lc_num], 
					    delta);
  VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					    Jktebop->gravdark2_vary,
					    Jktebop->gravdark2[lc_num], 
					    delta);
  /* Reflection */
  delta = 0.05;
  VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					    Jktebop->reflection1_vary,
					    Jktebop->reflection1[lc_num], 
					    delta);
  VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					    Jktebop->reflection2_vary,
					    Jktebop->reflection2[lc_num], 
					    delta);

  /* L3 */
  delta = 0.05;
  VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					    Jktebop->L3_vary,
					    Jktebop->L3[lc_num], 
					    delta);
  
  /* tidal lag */
  delta = 10.0;
  VARTOOLS_incrementparameters_foramoeba(&Nparameters, &Ntovary, &p, &ia,
					    Jktebop->tidallag_vary,
					    Jktebop->tidallag[lc_num], 
					    delta);

  /* Count the number of degrees of freedom */
  *Ndof = -Ntovary;
  for(i=0; i < N; i++) {
    if(!isnan(mag[i])) {
      (*Ndof)++;
    }
  }

  if(Ntovary > 0)
    {
      ftol = JKTEBOP_AMOEBA_CONVERGENCELIMIT;
      /* Get the chi2vals for the points on the simplex */
      VARTOOLS_amoeba_initializesimplexchi2(Nparameters, Ntovary, p, &chi2vals, 
					    &chisqjktebop, N, t, mag, err, (void *)(&fitstruct));
      /* Do the fit */
      amoeba_val = VARTOOLS_amoeba(p, chi2vals, ia, Nparameters, ftol, &chisqjktebop, &nfunk, 0, N, t, mag, err, (void *)(&fitstruct));
      if(!amoeba_val) {
	/* The solution converged, update the parameters */
	*best_chi2 = chi2vals[0];
	j = 0;
	for(i=1; i < Ntovary+1; i++) {
	  if(chi2vals[i] < *best_chi2) {
	    *best_chi2 = chi2vals[i];
	    j = i;
	  }
	}
	/* The best parameters are stored in p[j], 
	   Run chi2 again to reset the model lc to the best-fit model
	   and then copy the parameters over to 
	   the Jktebop struct. */
	
	*best_chi2 = chisqjktebop(p[j], Nparameters, N, t, mag, err, (void *) (&fitstruct));

	if(fitstruct.is_vary_P && fitstruct.is_vary_T0) {
	  T1 = p[j][0]; T2 = p[j][1];
	  Jktebop->T0[lc_num] = T1;
	  Jktebop->Period[lc_num] = (T2 - T1)/((double)(fitstruct.N2 - fitstruct.N1));
	} else {
	  Jktebop->T0[lc_num] = p[j][0];
	  Jktebop->Period[lc_num] = p[j][1];
	}
	Jktebop->r1r2[lc_num] = p[j][2];
	Jktebop->r2_r1[lc_num] = p[j][3];
	Jktebop->M2_M1[lc_num] = p[j][4];
	Jktebop->J2_J1[lc_num] = p[j][5];

	if(fitstruct.use_i_or_b == JKTEBOP_USEI) {
	  Jktebop->incl[lc_num] = p[j][6];
	  Jktebop->bimpact[lc_num] = fitstruct.bval;
	}
	else {
	  Jktebop->bimpact[lc_num] = p[j][6];
	  Jktebop->incl[lc_num] = fitstruct.ival;
	}

	Jktebop->esinomega[lc_num] = p[j][7];
	Jktebop->ecosomega[lc_num] = p[j][8];
	
	if(Jktebop->LD1law == JKTEBOP_LDLAW_LINEAR) {
	  Jktebop->LD1_coeffs[lc_num][0] = p[j][9];
	} else {
	  Jktebop->LD1_coeffs[lc_num][0] = p[j][9];
	  Jktebop->LD1_coeffs[lc_num][1] = p[j][10];
	}
	if(!Jktebop->LD2law == JKTEBOP_LDLAW_LOCKLD1) {
	  if(Jktebop->LD2law == JKTEBOP_LDLAW_LINEAR) {
	    Jktebop->LD2_coeffs[lc_num][0] = p[j][11];
	  } else {
	    Jktebop->LD2_coeffs[lc_num][0] = p[j][11];
	    Jktebop->LD2_coeffs[lc_num][1] = p[j][12];
	  }
	}
	
	Jktebop->gravdark1[lc_num] = p[j][13];
	Jktebop->gravdark2[lc_num] = p[j][14];
	Jktebop->reflection1[lc_num] = p[j][15];
	Jktebop->reflection2[lc_num] = p[j][16];
	Jktebop->L3[lc_num] = p[j][17];
	Jktebop->tidallag[lc_num] = p[j][18];
	
      }
      else *best_chi2 = -1.;
    }
  else
    {
      /* Nothing is being varied, just do the calculation with the current
	 parameters to get its chi2 value */
      *best_chi2 = chisqjktebop(p[0], Nparameters, N, t, mag, err, (void *) (&fitstruct));
    }
  VARTOOLS_amoeba_cleanup(&Nparameters, &Ntovary, &p, &ia, &chi2vals);
  return amoeba_val;
}

#define MASSIVECHI2 1.0e99

double chisqjktebop(double *p, int ma, int N, double *t, double *mag, double *sig, void *userparams)
/* This function returns the chi-square value for a given JKTEBOP model.
   It has the form expected for use with the VARTOOLS_amoeba function.
   p -- vector storing parameters that might be varied.
   ma - Number of parameters.
   N - Number of points in the lc.
   t - light curve times.
   mag - light curve magnitudes.
   sig - light curve uncertainties.
   userparams - Any other optional user data that is passed along to this
                function. Here we define a _Jktebop_fitstruct type data
		structure which we use to pass information such as the
		limb darkening laws, and whether bimpact or the inclination
		is used as the fitting parameter.
 */
{
  _Jktebop_fitstruct *fitstruct;
  
  double P, T0, T1, T2, r1r2, r2_r1, M2_M1, J2_J1, bimpact, incl;
  double esinomega, ecosomega, LD1a, LD1b, LD2a, LD2b;
  double gravdark1, gravdark2, reflection1, reflection2, L3, tidallag;
  double chisqval;
  double magshift, val1, val2;
  int i, use_i_or_b, LD1law, LD2law;
  double *simlc;

  void simulateEBOPEBlc(int N, double *JD, double *magout, double T0, double period, double esinomega, double ecosomega, int use_i_or_bimpact, double *incl, double *bimpact, double r1r2, double r2_r1, double q, double J2J1, double L3, double grav1, double grav2, double reflection1, double reflection2, double quadA1, double quadB1, double quadA2, double quadB2, double tidallag, int ldtype1, int ldtype2);

  fitstruct = (_Jktebop_fitstruct *) userparams;

  /* Translate the parameters in the input p vector into the 
     physical parameters that are varied */

  /* If both P and T0 are varied, then the free parameters
     will actually be T1 and T2, which are times associated
     with fixed eclipse numbers N1 and N2. If they are not
     both varied then we assume the parameters are P and T0 */
  if(fitstruct->is_vary_P && fitstruct->is_vary_T0) {
    T1 = p[0]; T2 = p[1];
    T0 = T1; P = (T2 - T1)/((double)(fitstruct->N2 - fitstruct->N1));
  } else {
    T0 = p[0]; P = p[1];
  }
  r1r2 = p[2];
  r2_r1 = p[3];
  M2_M1 = p[4];
  J2_J1 = p[5];
  use_i_or_b = fitstruct->use_i_or_b;
  if(use_i_or_b == JKTEBOP_USEI) {
    incl = p[6];
    bimpact = 0.;
  } else {
    incl = 0.;
    bimpact = p[6];
  }
  esinomega = p[7];
  ecosomega = p[8];
  LD1law = fitstruct->LD1law;
  switch(LD1law) {
  case JKTEBOP_LDLAW_LINEAR:
    LD1a = p[9];
    LD1b = 0.;
    break;
  default:
    LD1a = p[9];
    LD1b = p[10];
  }
  LD2law = fitstruct->LD2law;
  switch(LD2law) {
  case JKTEBOP_LDLAW_LOCKLD1:
    LD2a = LD1a;
    LD2b = LD1b;
    LD2law = LD1law;
    break;
  case JKTEBOP_LDLAW_LINEAR:
    LD2a = p[11];
    LD2b = 0.;
    break;
  default:
    LD2a = p[11];
    LD2b = p[12];
  }
  gravdark1 = p[13];
  gravdark2 = p[14];
  reflection1 = p[15];
  reflection2 = p[16];
  L3 = p[17];
  tidallag = p[18];

  simlc = fitstruct->magsim;

  if(J2_J1 < 0) return MASSIVECHI2;
  if(use_i_or_b == JKTEBOP_USEI) {
    if(incl > 90.0 || incl < 0.) return MASSIVECHI2;
  } else {
    if(bimpact < 0) return MASSIVECHI2;
  }
  if(r2_r1 < 0) return MASSIVECHI2;
  if(M2_M1 < 0) return MASSIVECHI2;
  if(r1r2 < 0 || r1r2 > 1) return MASSIVECHI2;
  
  /* Generate the simulated lc */
  simulateEBOPEBlc(N, t, simlc, T0, P, esinomega, ecosomega, use_i_or_b, &incl, &bimpact, r1r2, r2_r1, M2_M1, J2_J1, L3, gravdark1, gravdark2, reflection1, reflection2, LD1a, LD1b, LD2a, LD2b, tidallag, LD1law, LD2law);

  /* Update the incl or bimpact, whichever was not varied */
  if(use_i_or_b == JKTEBOP_USEI) {
    fitstruct->bval = bimpact;
  } else {
    fitstruct->ival = incl;
  }

  /* Determine the magnitude shift to fit to the lc */
  val1 = 0.; val2 = 0.;
  for(i=0; i < N; i++) {
    if(!isnan(mag[i])) {
      val1 += (mag[i] - simlc[i])/sig[i]/sig[i];
      val2 += 1./sig[i]/sig[i];
    }
  }
  magshift = val1 / val2;

  /* Calculate chi^2 */
  chisqval = 0.;
  for(i=0; i < N; i++) {
    simlc[i] += magshift;
    if(!isnan(mag[i])) {
      val1 = (mag[i] - simlc[i]);
      val1 *= val1;
      chisqval += val1/sig[i]/sig[i];
    }
  }
  return(chisqval);

}

/* Use the F77_FUNC macro generated by the configure script and put in
   ../../config.h (which in turn in included here through vartools.h) to
   determine what name to use to call the fortran GETMODEL subroutine
   defined in jktebop_lib.f.  Note that the call to F77_FUNC takes two
   arguments, the first is the name of the subroutine all lower-case, the
   second is the name of the subroutine all upper-case. */
#define GETMODEL_FUNC F77_FUNC(getmodel,GETMODEL)

void simulateEBOPEBlc(int N, double *JD, double *magout, double T0, double period, double esinomega, double ecosomega, int use_i_or_bimpact, double *incl, double *bimpact, double r1r2, double r2_r1, double q, double J2J1, double L3, double grav1, double grav2, double reflection1, double reflection2, double quadA1, double quadB1, double quadA2, double quadB2, double tidallag, int ldtype1, int ldtype2)
{
  /* This function generates a simulated detached EB lc observed at
     the times given in JD by calling the FORTRAN GETMODEL subroutine
     included in libjktebop.f. The input parameters are:
   N - number of lc points.
   JD - observed JD times.
   magout - output simulated magnitudes.
   T0 - time of primary eclipse center in same units as JD.
   period - period of binary in same units as JD.
   esinomega, ecosomega - the usual eccentricity/arg of periastron parameters.
   use_i_or_bimpact - 0 to use the inclination, or 1 to use the impact parameter
                      the other quantity will be calculated and its value will
                      be returned. This is why incl and bimpact are passed as
                      pointers.
   incl - orbital inclination in degrees.
   bimpact -  the impact parameter of the primary eclipse. This is the
              projected distance between the centers of the two stars
              at the time T0, divided by the sum of the radii of the stars.
              It goes from 0 for a central eclipse, to 1 for a just
              grazing eclipse.
   r1r2 - R1 + R2, where R_j is the radius of star j as a fraction of the 
          semi-major axis.
   r2_r1 - R2 / R1.
   q = mass of star2 / mass of star1.
   J2J1 = surface brightness of star 2 / surface brightness of star 1.
   L3 - third light.
   grav1 - gravity darkening coefficient for primary star.
   grav2 - gravity darkening coefficient for secondary star.
   reflection1 - reflection coefficient for primary star, this will be calculated
                 for values <= 0.
   reflection2 - reflection coefficient for secondary star.
   quadA1, quadB1 - first and second limb darkening coefficients for the primary star.
   quadA2, quadB2 - first and second limb darkening coefficients for the secondary star.
   tidallag - tidal lead/lag angle (deg).
   ldtype1 - limb darkening law for star 1 (1 = linear
                                            2 = logarithmic
                                            3 = square-root
                                            4 = quadratic)
   ldtype2 - limb darkening law for star 2. (Previous options, plus 
                                             5 = coefficients fixed to
                                             that of star 1).
  */
  int i, type;
  double t0, mag;
  double cosi, eccen, omega;
  double V[22];
  int LDTYPE[2];
  double LP, LS;

  /* Prepare the parameter vector for passing to jktebop */
  V[0] = J2J1;
  V[1] = r1r2; /* sum of fractional radii */
  V[2] = r2_r1; /* ratio of stellar radii */
  V[3] = quadA1; /* first limb darkening coefficient for star A. */
  V[4] = quadA2; /* first limb darkening coefficient for star B. */

  /* Calculate the orbital inclination from the impact parameter, or
     vice versa */
  if(use_i_or_bimpact == JKTEBOP_USEI) {
    V[5] = *incl; /* orbital inclination */
    eccen = sqrt(fabs(ecosomega*ecosomega + esinomega*esinomega));
    if(eccen > 0.)
      omega = atan2(esinomega,ecosomega);
    else
      omega = 0.;
    cosi = cos((*incl)*M_PI/180.);
    *bimpact = cosi*(1. - eccen*eccen)/(1. + eccen*cos(M_PI - omega))/r1r2;
  } else {
    eccen = sqrt(fabs(ecosomega*ecosomega + esinomega*esinomega));
    if(eccen > 0.)
      omega = atan2(esinomega,ecosomega);
    else
      omega = 0.;
    cosi = (*bimpact)*(1. + eccen*cos(M_PI - omega))*r1r2/(1. - eccen*eccen);
    *incl = 180*acos(cosi)/M_PI;
    V[5] = *incl;
  }

  V[6] = ecosomega;
  V[7] = esinomega;
  /* Gravity darkening coefficients */
  V[8] = grav1;
  V[9] = grav2;

  /* Reflected light coefficients, these will be determined automatically. */
  if(reflection1 <= 0.)
    V[10] = 0.;
  else
    V[10] = reflection1;
  if(reflection2 <= 0.)
    V[11] = 0.;
  else
    V[11] = reflection2;

  V[12] = q; /* mass2 / mass1 */
  V[13] = tidallag; /* tidal lead/lag angle */
  V[14] = L3; /* Third light */
  V[15] = 0.; /* This is the phase correction factor, I don't see the
                 point of adjusting this rather than P and T0, so I
	         set it to 0.
              */
  V[16] = 0.0; /* Zero-point magnitude */

  V[17] = 5.; /* Integration ring size (deg) */

  V[18] = period; /* orbital period (days) */

  V[19] = T0; /* ephemeris timebase (days) */

  V[20] = quadB1; /* second limb darkening coeff. for star A. */
  V[21] = quadB2; /* second limb darkening coeff. for star B. */
  
  /* set limb darkening law types */
  LDTYPE[0] = ldtype1;
  LDTYPE[1] = ldtype2;

  /* Find reflection coefficients */
  t0 = T0;

  /* this will tell GETMODEL that we want the output in magnitudes */
  type = 1;

  /* This will call the GETMODEL function in the file jktebop_lib.f to
     set the reflection coefficients. Use GETMODEL_FUNC which is set
     using the #define before this function to whatever is the appropriate
     expression to use for calling this FORTRAN function from c */

  if(reflection1 <= 0. || reflection2 <= 0.) {
    mag = GETMODEL_FUNC(V,LDTYPE,&t0, &type, &LP, &LS);
    if(reflection1 <= 0.)
      V[10] = 0.4 * LS * pow((V[1]/(1.0 + V[2])),2.);
    if(reflection2 <= 0.)
      V[11] = 0.4 * LP * pow((V[1]/(1.0 + (1./V[2]))),2.);
  }
  
  /* Now get the model lc */
  for(i=0;i<N;i++)
    {
      magout[i] = GETMODEL_FUNC(V,LDTYPE,&(JD[i]),&type,&LP,&LS);
    }
}

