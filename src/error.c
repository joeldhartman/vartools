/*     This file is part of VARTOOLS version 1.31                      */
/*                                                                           */
/*     VARTOOLS is free software: you can redistribute it and/or modify      */
/*     it under the terms of the GNU General Public License as published by  */
/*     the Free Software Foundation, either version 3 of the License, or     */
/*     (at your option) any later version.                                   */
/*                                                                           */
/*     This program is distributed in the hope that it will be useful,       */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/*     GNU General Public License for more details.                          */
/*                                                                           */
/*     You should have received a copy of the GNU General Public License     */
/*     along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*                                                                           */
/*     Copyright 2007, 2008, 2009  Joel Hartman                              */
/*                                                                           */
#include "commands.h"
#include "programdata.h"
#include "functions.h"

void error(int errflag)
{
  switch(errflag)
    {
    case 0:
      exit(0);
    case ERR_USAGE:
      fprintf(stderr,"No commands issued\n");
      exit(ERR_USAGE);
    case ERR_FILENOTFOUND:
      fprintf(stderr,"Error: File not found\n");
      exit(ERR_FILENOTFOUND);
    case ERR_MEMALLOC:
      fprintf(stderr,"Memory Allocation Error\n");
      exit(ERR_MEMALLOC);
    case ERR_WRONGITER:
      fprintf(stderr,"Error: the clipping iteration flag must be 0 or 1\n");
      exit(ERR_WRONGITER);
    case ERR_WRONGORDER:
      fprintf(stderr,"Error: the order for all decorrelation terms must be greater than 0\n");
      exit(ERR_WRONGORDER);
    case ERR_KILLHARM_NOAOV:
      fprintf(stderr,"Error: to use aov for the killharm period you must run an aov or aov_harm command before running killharm\n");
      exit(ERR_KILLHARM_NOAOV);
    case ERR_KILLHARM_NOLS:
      fprintf(stderr,"Error: to use ls for the killharm period you must run an LS command before running killharm\n");
      exit(ERR_KILLHARM_NOLS);
    case ERR_KILLHARM_NOBOTH:
      fprintf(stderr,"Error: to use both aov and ls for the killharm periods you must run both an aov or aov_harm and LS commands before running killharm\n");
      exit(ERR_KILLHARM_NOBOTH);
    case ERR_KILLHARM_NOBLS:
      fprintf(stderr,"Error: to use bls for the phase period you must run a bls command before running phase\n");
      exit(ERR_KILLHARM_NOBLS);
    case ERR_KILLHARM_NOBLSFIXPER:
      fprintf(stderr,"Error: to use blsfixper for the phase period you must run a blsfixper command before running this command\n");
      exit(ERR_KILLHARM_NOBLSFIXPER);
    case ERR_KILLHARM_WRONGNPER:
      fprintf(stderr,"Error: the number of periods to use for killharm must be greater than 0\n");
      exit(ERR_KILLHARM_WRONGNPER);
    case ERR_READFORMAT:
      fprintf(stderr,"Error: coljd, colmag and colsig must be distinct and greater than zero\n");
      exit(ERR_READFORMAT);
    case ERR_NEEDLIST:
      fprintf(stderr,"Error: At least one of the commands requires an input list\n");
      exit(ERR_NEEDLIST);
    case ERR_INVALIDGLOBALDECORR:
      fprintf(stderr,"Error: The global_files for decorrelation should all be the same length!\n");
      exit(ERR_INVALIDGLOBALDECORR);
    case ERR_CODEERROR:
      fprintf(stderr,"If you receive this error message, there must be a bug in the code. Please notify J. Hartman giving details of the command that you issued so that the bug can be tracked down and eliminated.\n");
      exit(ERR_CODEERROR);
    case ERR_NOTENOUGHSTARS_ERS:
      fprintf(stderr,"Error: Not enough stars to perform ensemble sigma rescaling\n");
      exit(ERR_NOTENOUGHSTARS_ERS);
    case ERR_SINGULARMATRIX:
      fprintf(stderr,"Error: Singular matrix in routine gaussj. Aborting\n");
      exit(ERR_SINGULARMATRIX);
    case ERR_TOOMANYAMOEBAITERATIONS:
      fprintf(stderr,"Error: Too many iterations in the function amoeba. Aborting\n");
      exit(ERR_TOOMANYAMOEBAITERATIONS);
    case ERR_BLSNBMAX:
      fprintf(stderr,"Error: the number of bins for BLS must be less than 2000.\n");
      exit(ERR_BLSNBMAX);
    case ERR_BLSFMINTOOSMALL:
      fprintf(stderr,"Error: The minimum frequency is less than 1/T in BLS. Aborting\n");
      exit(ERR_BLSFMINTOOSMALL);
    case ERR_BLSNOFREQ:
      fprintf(stderr,"Error: Not enough frequencies survived the BLS clipping to find the specified number of peaks.\ny");
      exit(ERR_BLSNOFREQ);
    case ERR_SYSREMUSEWEIGHTS:
      fprintf(stderr,"Error: useweights for SYSREM should be 0, 1 or 2.\n");
      exit(ERR_SYSREMUSEWEIGHTS);
    case ERR_UNSORTEDLIGHTCURVE:
      fprintf(stderr,"Error: one of the light curves is not sorted in time.\n");
      exit(ERR_UNSORTEDLIGHTCURVE);
    case ERR_BINSIZEZERO:
      fprintf(stderr,"Error: The number of bins in the binlc routine is less than 1.\n");
      exit(ERR_BINSIZEZERO);
    case ERR_BADTYPE:
      fprintf(stderr,"Error: bad column type - if you get this error there is a bug in the program. Please report it to jhartman@cfa.harvard.edu, sorry for the inconvenience!.\n");
      exit(ERR_BADTYPE);
    case ERR_NOSTRINGIDCOLUMN:
      fprintf(stderr,"Error: -matchstringid specified but the stringid column has not been specified for one of the light curves.\n");
      exit(ERR_NOSTRINGIDCOLUMN);
    case ERR_KILLHARM_NOINJECTHARM:
      fprintf(stderr,"Error: to use injectharm for the killharm period you must run an Injectharm command before running killharm\n");
      exit(ERR_KILLHARM_NOINJECTHARM);
    case ERR_KILLHARM_NEGATIVEPERIOD:
      fprintf(stderr,"Error: fixed periods input to killharm must be greater than zero.\n");
      exit(ERR_KILLHARM_NEGATIVEPERIOD);
    case ERR_INVALIDOUTPUTFORMAT:
      fprintf(stderr,"Error: Invalid output format string.\n");
      exit(ERR_INVALIDOUTPUTFORMAT);
    case ERR_FITSERROR:
      fprintf(stderr,"Error: error reading fits file.\n");
      exit(ERR_FITSERROR);
#ifdef PARALLEL
    case ERR_NPROC_TOO_SMALL:
      fprintf(stderr,"Error: Nproc must be greater than 0.\n");
      exit(ERR_NPROC_TOO_SMALL);
#endif
    case ERR_EQUALTIMES:
      fprintf(stderr,"Error: Light Curves cannot contain points with identical time measurements.\n");
      exit(ERR_EQUALTIMES);
    case ERR_UTCINPUT_FITS:
      fprintf(stderr,"Error: UTC input is currently not supported for binary fits light curves.\n");
      exit(ERR_UTCINPUT_FITS);
    case ERR_NO_EPHEM_FILE:
      fprintf(stderr,"Error: Attempt to perform BJD and/or UTC/TDB time conversion without providing a solar system ephemeris file. See 'vartools -help -converttime' for details on where to obtain an ephemeris file.\n");
      exit(ERR_NO_EPHEM_FILE);
    case ERR_NO_LEAPSEC_FILE:
      fprintf(stderr,"Error: Attempt to perform UTC/TDB time conversion without providing a leap-second file. See 'vartools -help -converttime' for details on where to obtain a leap-second file.\n");
      exit(ERR_NO_LEAPSEC_FILE);
    case ERR_NO_PLANETDATA_FILE:
      fprintf(stderr,"Error: Attempt to perform BJD and/or UTC/TDB time conversion without providing a CSPICE planet-data file. See 'vartools -help -converttime' for details on where to obtain such a file.\n");
      exit(ERR_NO_PLANETDATA_FILE);
    case ERR_INVALID_UTC_FORMAT:
      fprintf(stderr,"Error: Invalid format for UTC input conversion.\n");
      exit(ERR_INVALID_UTC_FORMAT);
    case ERR_CONVERTTIME_NORADEC:
      fprintf(stderr,"Error: RA/DEC must be specified for conversions to/from BJD or HJD.\n");
      exit(ERR_CONVERTTIME_NORADEC);
#ifdef DYNAMICLIB
    case ERR_INVALIDUSEOFPARSEFIXSPECFIXCOLUMN:
      fprintf(stderr,"Error: Invalid use of the ParseFixSpecFixcolumn function, this error indicates a bug in a dynamically linked library.\n");
      exit(ERR_INVALIDUSEOFPARSEFIXSPECFIXCOLUMN);
#endif
    case ERR_TOOMANYREADFORMATINPUTLCFORMAT:
      fprintf(stderr,"Error: the -readformat and -inputlcformat options should be used only once, additionally they cannot both be used.\n");
      exit(ERR_TOOMANYREADFORMATINPUTLCFORMAT);
    case ERR_WRONGTYPEFORSPECIALVARIABLE:
      fprintf(stderr,"Error: attempt to define a special variable with the wrong type. The variable \"t\" must have type \"double\" or \"utc\", the variables \"mag\" and \"err\" must have type \"double\", and the variable \"id\" must have type \"string\".\n");
      exit(ERR_WRONGTYPEFORSPECIALVARIABLE);
    case ERR_INVALIDVARIABLEFORCHANGEVAR:
      fprintf(stderr,"Error: using a variable of the wrong type in the command -changevariable.\n");
      exit(ERR_INVALIDVARIABLEFORCHANGEVAR);
    case ERR_MYSORT_GENERIC_BADCALL:
      fprintf(stderr,"Error: bad call to the function mysort_generic. If you get this message there is a bug in the program.\n");
      exit(ERR_MYSORT_GENERIC_BADCALL);
    case ERR_CREATEVARIABLE_OUTCOLUMN_NEEDS_VPTRINPUT:
      fprintf(stderr,"Error: bad call to the function CreateVariable. If you get this message there is a bug in the program. When the outcolumn vectortype is passed to the data, a non-NULL input data pointer is required.\n");
      exit(ERR_CREATEVARIABLE_OUTCOLUMN_NEEDS_VPTRINPUT);
    case ERR_OUTCOLUMN_ON_LHS:
      fprintf(stderr,"Error: an output column variable cannot appear on the left-hand-side of an -expr command.\n");
      exit(ERR_OUTCOLUMN_ON_LHS);
    case ERR_BADVARIABLETYPE_STATSCOMMAND:
      fprintf(stderr,"Error: the stats command can only be used on light curve type vectors.\n");
      exit(ERR_BADVARIABLETYPE_STATSCOMMAND);
    case ERR_BADIFTHENELSE:
      fprintf(stderr,"Error: invalid -if -elif -else -fi conditional. This could be due to calling -elif, -else, or -fi without a preceding -if, or calling -elif after an -else.\n");
      exit(ERR_BADIFTHENELSE);
    case ERR_TOOMANYINLISTVARS:
      fprintf(stderr,"Error: only one call to -inlistvars can be used.\n");
      exit(ERR_TOOMANYINLISTVARS);
    case ERR_INVALID_CALL_GETCOLUMNNAMEFORRECENTCOMMAND:
      fprintf(stderr,"Error: invalid call to the function GetColumnNameForRecentCommand, attempt to pass it an unsupported command and/or unsupported parameter.\n");
      exit(ERR_INVALID_CALL_GETCOLUMNNAMEFORRECENTCOMMAND);
    case ERR_BINARYLC_KEYWORDTOOLONG:
      fprintf(stderr,"Error: the keyword in the header of an input binary-format light curve is longer than the space allowed.\n");
      exit(ERR_BINARYLC_KEYWORDTOOLONG);
    case ERR_BADTYPE_BINARYLC:
      fprintf(stderr,"Error: string and convertjd formats are not currently supported for binary light curves.\n");
      exit(ERR_BADTYPE_BINARYLC);
    case ERR_FITS_CREATETABLE:
      fprintf(stderr,"Error: cannot create a new fits table.\n");
      exit(ERR_FITS_CREATETABLE);
    case ERR_FITS_WRITECOLUMN:
      fprintf(stderr,"Error: cannot write a column to a fits table.\n");
      exit(ERR_FITS_WRITECOLUMN);
    case ERR_MCMCINVALIDMAXMEM:
      fprintf(stderr,"Error maxmemstore for mcmc must be > 0\n");
      exit(ERR_MCMCINVALIDMAXMEM);
    case ERR_READALL_ANDCOPYLC:
      fprintf(stderr,"Error the -copylc command cannot be used with the -readall option.\n");
      exit(ERR_READALL_ANDCOPYLC);
    case ERR_NEGATIVE_COVAR_PARAM:
      fprintf(stderr,"Error the initialization expressions for one of the covariance parameters used in a -nonlinfit command evaluates to a negative number, or to zero. These parameters must all be greater than zero.\n");
      exit(ERR_NEGATIVE_COVAR_PARAM);
    case ERR_NEGATIVE_COVAR_PARAM_ADDNOISE:
      fprintf(stderr,"Error: one of the covariance parameters in the -addnoise command is not greater than zero.\n");
      exit(ERR_NEGATIVE_COVAR_PARAM_ADDNOISE);
    case ERR_INTERP_NOTENOUGHPOINTS:
      fprintf(stderr,"Error: Not enough points in input light curve to perform interpolation. There must be at least two points in the light curve with unequal times.\n");
      exit(ERR_INTERP_NOTENOUGHPOINTS);
    case ERR_INVALIDEXECCOMMANDSTRFORMAT:
      fprintf(stderr,"Error: unable to parse the shell command given after the \"opencommand\" keyword.\n");
      exit(ERR_INVALIDEXECCOMMANDSTRFORMAT);
    case ERR_FUNCTIONCALL_LENINVALIDOPERAND:
      fprintf(stderr,"Error: the argument to the len() function must be a single variable or number, it does not accept expressions.\n");
      exit(ERR_FUNCTIONCALL_LENINVALIDOPERAND);
    case ERR_NEEDCSPICE:
      fprintf(stderr,"Error: To perform BJD and/or UTC/TDB time conversion you need to recompile VARTOOLS against the NASA CSPICE library.\n");
      exit(ERR_NEEDCSPICE);
    default:
      fprintf(stderr,"Error - Unspecified Error\n");
      exit(999);
    }
}

void error2_noexit(int errflag, char *s)
{
  switch(errflag)
    {
    case ERR_FILENOTFOUND:
      fprintf(stderr,"Cannot Open %s\n",s);
      break;
    case ERR_BINARYLIGHTCURVE_INVALIDFORMAT:
      fprintf(stderr,"Unable to read the binary light curve file %s. It appears to have an invalid format.\n", s);
      break;
    default:
      fprintf(stderr,"Undefined Error\n");
      break;
    }
}

void error2(int errflag, char *s)
{
  switch(errflag)
    {
    case ERR_FILENOTFOUND:
      fprintf(stderr,"Cannot Open %s\n",s);
      exit(ERR_FILENOTFOUND);
    case ERR_INPUTMISSINGCOLUMN:
      fprintf(stderr,"File %s does not have enough columns for the specified read-format and decorr commands.\n",s);
      exit(ERR_INPUTMISSINGCOLUMN);
    case ERR_CANNOTWRITE:
      fprintf(stderr,"Cannot write to %s\n",s);
      exit(ERR_CANNOTWRITE);
    case ERR_GETLSAMPTHRESH_FILETOSHORT:
      fprintf(stderr,"The signal list file %s does not have as many entries as there are light curves.\n",s);
      exit(ERR_GETLSAMPTHRESH_FILETOSHORT);
    case ERR_SIGFILEWRONGLENGTH:
      fprintf(stderr,"The signal file %s is not the same length as the corresponding light curve.\n",s);
      exit(ERR_SIGFILEWRONGLENGTH);
    case ERR_HELP_BADCOMMAND:
      fprintf(stderr,"Invalid command or option \"%s\" to a -listcommands or -help option\n",s);
      exit(ERR_HELP_BADCOMMAND);
    case ERR_EXAMPLE_BADCOMMAND:
      fprintf(stderr,"Error executing -example command. Either the command \"%s\" is invalid, or an example has not yet been included.\n",s);
      exit(ERR_EXAMPLE_BADCOMMAND);
    case ERR_MISSINGSAVELC:
      fprintf(stderr,"There are fewer than %s -savelc commands preceding a call to \"-restorelc %s\"\n",s,s);
      exit(ERR_MISSINGSAVELC);
    case ERR_NOCOLUMN:
      fprintf(stderr,"%s is not one of the output column numbers or names.\n",s);
      exit(ERR_NOCOLUMN);
    case ERR_BADCOLUMNLINK:
      fprintf(stderr,"Error: Attempt to link to column %s from a command that precedes the command that outputs this column.\n",s);
      exit(ERR_BADCOLUMNLINK);
    case ERR_CANNOTOPEN:
      fprintf(stderr,"Error: Cannot open the file %s\n",s);
      exit(ERR_CANNOTOPEN);
    case ERR_IMAGEHDU:
      fprintf(stderr,"Error: %s is a fits image, this routine only supports fits tables.\n",s);
      exit(ERR_IMAGEHDU);
    case ERR_UNKNOWNOBSERVATORY:
      fprintf(stderr,"Error: Unknown Observatory Code %s\n",s);
      exit(ERR_UNKNOWNOBSERVATORY);
#ifdef DYNAMICLIB
    case ERR_LIBRARY_MISSING_FUNCTION:
      fprintf(stderr,"Error: Dynamically Loaded Library does not have a function with the name %s\n",s);
      exit(ERR_LIBRARY_MISSING_FUNCTION);
    case ERR_INVALID_CALL_REGISTERDATAVECTOR:
      fprintf(stderr,"Error: An invalid data-source type has been sent to the VARTOOLS_RegisterDataVector function from a dynamically loaded library.\n");
      exit(ERR_INVALID_CALL_REGISTERDATAVECTOR);
#endif
    case ERR_FUNCTIONCALL_INVALIDNEXPR:
      fprintf(stderr,"Error: Attempt to call the function %s with an invalid number of expressions.\n",s);
      exit(ERR_FUNCTIONCALL_INVALIDNEXPR);
    case ERR_ANALYTICPARSE:
      fprintf(stderr,"Error: Cannot parse the analytic expression %s\n",s);
      exit(ERR_ANALYTICPARSE);
    case ERR_INVALIDARGUMENTTOEXPR:
      fprintf(stderr,"Error: Invalid argument \"%s\" supplied to the -expr command.\n",s);
      exit(ERR_INVALIDARGUMENTTOEXPR);
    case ERR_VARIABLEALREADYINUSE:
      fprintf(stderr,"Error: The variable \"%s\" is defined more than once in an -inputlcformat option.\n",s);
      exit(ERR_VARIABLEALREADYINUSE);
    case ERR_BADVARIABLENAME:
      fprintf(stderr,"Error: Invalid variable name \"%s\". Variables must contain only letters, numbers, and the character '_'. Variable names cannot begin with a number.\n",s);
      exit(ERR_BADVARIABLENAME);
    case ERR_BADCOLUMNFORMATSTRING:
      fprintf(stderr,"Error: Invalid format string \"%s\" given to the command \"-o\" after the keyword \"columnformat\". See \"vartools -help -o\" for the correct syntax.\n",s);
      exit(ERR_BADCOLUMNFORMATSTRING);
    case ERR_UNDEFINEDVARIABLE:
      fprintf(stderr,"Error: Attempt to use an undefined/initialized variable %s.\n",s);
      exit(ERR_UNDEFINEDVARIABLE);
    case ERR_INVALIDFUNCTIONLINFIT:
      fprintf(stderr,"Error: Invalid function \"%s\" given to linfit. Note that the function must be linear in its fitted parameters.\n", s);
      exit(ERR_INVALIDFUNCTIONLINFIT);
    case ERR_LINFITMISSINGPARAM:
      fprintf(stderr,"Error: \"%s\" is included in the list of parameters to a -linfit command, but it is not included in the function.\n", s);
      exit(ERR_LINFITMISSINGPARAM);
    case ERR_INVALIDSTATISTIC:
      fprintf(stderr,"Error: Invalid statistic name \"%s\" given to the -stats command or as an option to the -nonlinfit command.\n", s);
      exit(ERR_INVALIDSTATISTIC);
    case ERR_RESERVEDVARIABLENAME:
      fprintf(stderr,"Error: the variable name \"%s\" is reserved, and  cannot be used with the -inlistvars option.\n", s);
      exit(ERR_RESERVEDVARIABLENAME);
    case ERR_NO_PREVIOUS_COMMAND:
      fprintf(stderr,"Error: there is no previous %s command.\n", s);
      exit(ERR_NO_PREVIOUS_COMMAND);
    case ERR_BINARYLIGHTCURVE_INVALIDFORMAT:
      fprintf(stderr,"Unable to read the binary light curve file %s. It appears to have an invalid format.\n", s);
      exit(ERR_BINARYLIGHTCURVE_INVALIDFORMAT);
      break;
    case ERR_BINARYLIGHTCURVE_MISSINGDATA:
      fprintf(stderr,"The binary light curve file %s is shorter than what is indicated in its header\n", s);
      exit(ERR_BINARYLIGHTCURVE_MISSINGDATA);
      break;
    case ERR_MISSING_BINARYLC_HEADERNAME:
      fprintf(stderr,"The header keyword %s is not found in the first binary light curve read-in.\n", s);
      exit(ERR_MISSING_BINARYLC_HEADERNAME);
      break;
    case ERR_MISSING_FITSLC_HEADERNAME:
      fprintf(stderr,"No column with the name %s is found in the first fits light curve read-in.\n", s);
      exit(ERR_MISSING_FITSLC_HEADERNAME);
      break;
    case ERR_FILEEXISTS_NOCLOBBER:
      fprintf(stderr,"Attempt to write to existing file %s with noclobber flag enabled. Aborting.\n", s);
      exit(ERR_FILEEXISTS_NOCLOBBER);
      break;
    case ERR_INVALIDVARIABLEFORNONLINFIT:
      fprintf(stderr,"The variable \"%s\" has the wrong type to be used a parameter in the -nonlinfit command.\n", s);
      exit(ERR_INVALIDVARIABLEFORNONLINFIT);
      break;
    case ERR_BADNONLINFITPARAMINIT:
      fprintf(stderr,"Error parsing the term \"%s\" in reading the parameter list for a -nonlinfit command. The expected syntax is VARNAME=INITIALVALUE_EXPRESSION:INITIALUNCERTAINTY_EXPRESSION where VARNAME is the variable name of the parameter, INITIALVALUE_EXPRESSION is an analytic expression which will evaluate to the initial value to try for the parameter, and INITIALUNCERTAINTY_EXPRESSION is an analytic expression which will evaluate to the initial uncertainty to use for the parameter.\n", s);
      exit(ERR_BADNONLINFITPARAMINIT);
      break;
    case ERR_BADNONLINFITPRIORINIT:
      fprintf(stderr,"Error parsing the term \"%s\" in reading the list of priors for a -nonlinfit command. The expected syntax is VARNAME=PRIOR_EXPRESSION where VARNAME is the variable name of the parameter, and PRIOR_EXPRESSION is an analytic expression which will evaluate to -2*ln(P) where P is the prior probability for that parameter. Note that VARNAME must be one of the variables given on the parameter list. You will also get this same error message if you give a parameter in the prior list which is not in the parameter list for this command.\n",s);
      exit(ERR_BADNONLINFITPRIORINIT);
    case ERR_INVALID_PARAMETERVALUE:
      fprintf(stderr,"Error - Invalid Parameter Value Given on Command Line: %s.\n",s);
      exit(ERR_INVALID_PARAMETERVALUE);
    case ERR_NR:
      fprintf(stderr,"Error - The following error was encountered in a numerical recipes function: \"%s\"\n",s);
      exit(ERR_NR);
    case ERR_INVALIDVARIABLELCVARIABLE:
      fprintf(stderr,"The variable \"%s\" is being used with conflicting types on the command-line.\n", s);
      exit(ERR_INVALIDVARIABLEFORNONLINFIT);
      break;
    case ERR_INVALIDANALYTICFUNCTIONDEFINITION:
      fprintf(stderr,"Error - \"%s\" is not a valid format for defining an analytic function. The expected syntax has the form \"funcname(arg1,arg2,...,argN)=valid_analytic_expression\".\n",s);
      exit(ERR_INVALIDANALYTICFUNCTIONDEFINITION);
      break;
    case ERR_ANALYTICFUNCTIONDUPLICATEINPUTARG:
      fprintf(stderr,"Error - The variable name \"%s\" appears more than once in the input argument string to an analytic function defined through the \"-f\" option.\n",s);
      exit(ERR_ANALYTICFUNCTIONDUPLICATEINPUTARG);
      break;
    case ERR_FUNCNAMETOOLONG:
      fprintf(stderr,"Error - the function name \"%s\" is too long.\n",s);
      exit(ERR_FUNCNAMETOOLONG);
      break;
    case ERR_INDEXINGWRONGVARIABLETYPEINEXPRESSION:
      fprintf(stderr,"Error - attempting to index the variable \"%s\" which is not an array of other a type that supports indexing.\n",s);
      exit(ERR_INDEXINGWRONGVARIABLETYPEINEXPRESSION);
      break;
    case ERR_BADINDEXINGOFLHSVARIABLEINEXPRESSIONCOMMAND:
      fprintf(stderr,"Error - the term \"%s\", which includes a light curve vector, is given as part of an index range on the left hand side of an expression command. This is not allowed.\n",s);
      exit(ERR_BADINDEXINGOFLHSVARIABLEINEXPRESSIONCOMMAND);
      break;
    case ERR_BADVECTORTYPEFOROUTPUTCOLUMNVARIABLE:
      fprintf(stderr,"Error - the variable \"%s\" is not of the correct type to be written to the output ascii table. Light curve vectors, or variables storing ascii data (char or string) cannot be included in this list.\n",s);
      exit(ERR_BADVECTORTYPEFOROUTPUTCOLUMNVARIABLE);
      break;
    case ERR_PYTHONOUTPUTUNDEFINEDVARIABLE:
      fprintf(stderr,"Error - the variable \"%s\" included in the list of outputcolumns variables supplied to a -python command is not defined.\n",s);
      exit(ERR_PYTHONOUTPUTUNDEFINEDVARIABLE);
      break;
    case ERR_TFATEMPLATEFITSREADERROR:
      fprintf(stderr,"Error - an error was encountered in attempting to read the TFA template fits file \"%s\"\n",s);
      exit(ERR_TFATEMPLATEFITSREADERROR);
      break;
    default:
      fprintf(stderr,"Error - Unspecified Error\n");
      exit(999);
    }
}
