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
/*     This file is part of VARTOOLS version 1.152                      */
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
#include "outcolumn.h"

#ifdef PARALLEL
#include <pthread.h>
#include <semaphore.h>
#endif

#include "userlib.h"
#include "userfunc.h"

#include "ifelse.h"

#ifdef _USEBINARY_LC
#include "binarylcio.h"
#endif

#include "mcmcfit.h"

#include "doublelinklist.h"

#define INITCOMMANDSIZE 100

#define VARTOOLS_TYPE_DOUBLE 0
#define VARTOOLS_TYPE_STRING 1
#define VARTOOLS_TYPE_INT 2
#define VARTOOLS_TYPE_FLOAT 3
#define VARTOOLS_TYPE_LONG 4
#define VARTOOLS_TYPE_CHAR 5
#define VARTOOLS_TYPE_CONVERTJD 6
#define VARTOOLS_TYPE_SHORT 7
#define VARTOOLS_TYPE_USERDEF 8

#define VARTOOLS_SOURCE_INLIST 0
#define VARTOOLS_SOURCE_COMPUTED 1
#define VARTOOLS_SOURCE_FIXED 2
#define VARTOOLS_SOURCE_PRIORCOLUMN 3
#define VARTOOLS_SOURCE_LC 4
#define VARTOOLS_SOURCE_NONE 5
#define VARTOOLS_SOURCE_RECENTCOMMAND 6
#define VARTOOLS_SOURCE_EVALEXPRESSION 7

#define MAXIDSTRINGLENGTH 256

#define INCOLUMN_HEADER_INDICATOR 1

#define VARTOOLS_NONLINFIT_FITTYPE_AMOEBA 0
#define VARTOOLS_NONLINFIT_FITTYPE_DEMCMC 1
#define VARTOOLS_NONLINFIT_FITTYPE_GA 2

#define VARTOOLS_DEFAULT_NOUTPUT_BUFFERS 32

#define VARTOOLS_LC_DELIMTYPE_WHITESPACE 0
#define VARTOOLS_LC_DELIMTYPE_CHAR 1
#define VARTOOLS_LC_DELIMTYPE_STRING 2

//double sizeHISTvector;
double JDTOL;

typedef struct {
  char *data;
  int len;
} _StringBuffer;

typedef struct {
  int datatype;
  int Ncolumns;
  int *incolumns;
  char **incolumn_names;
  void *dataptr;
  char *scanformat;
  _Expression *expression;
  //_Variable *varptr;
} _DataFromInputList;

typedef struct {
  int maxstringlength;
  int datatype;
  int Ncolumns;
  int *incolumns;
  char **incolumn_names;
  char **incolumn_header_names;
  char *scanformat;
  int UTCindex[6];
  void *dataptr;
  _Variable *variable;
  _Expression *expression;
} _DataFromLightCurve;

typedef struct {
  int datatype;
  int Ncolumns;
  void *dataptr;
} _ScalarData;

#ifdef _USEBINARY_LC
typedef struct {
  int num_apertures;
  int num_pts;
  int hdr_size;
  int num_columns;
  int lc_record_size;
  char ***lc_header;
  int memsize_lc_header;
  char **lc_column_names;
  BinLC_OutputFormat *lc_column_format;
  int memsize_lc_columns;
} binarylightcurve;
#endif

typedef struct {
  struct LinkedList lcs_to_proc;


  int readformatused;
  int inputlcformatused;
  int inlistvars;
  int Nlcs;
  int Ncommands;
  int Nthread;
  double **mag, **sig, **t;
  char ***stringid;
  int **stringid_idx;
  int *NJD;
  _Variable *magvar, *sigvar, *tvar, *idvar;
  char lc_getcolumnsfromheader;
  char lc_getcolumnsfromheader_notyetset;
#ifdef PARALLEL
  pthread_mutex_t lc_getcolumnsfromheader_mutex;
#endif
#ifdef _USEBINARY_LC
  char binarylcinput;
  binarylightcurve *binlc;
#endif
  int readfromstdinflag;
  int readallflag;
  int listflag;
  int fileflag;
  int use_lc_open_exec_command;
  char *lc_open_exec_command_str;
  int tabflag;
  int header;
  int basename;
  int matchstringid;
  int requirestringid;
  int readimagestring;
  int colstringid;
  int Nskip;
  int Nskipchar;
  char *skipchars;
  int inputUTC;
  char UTCformat[MAXLEN];
  int UTCindex[6];
  int coljd;
  int colmag;
  int colsig;
  int *sizesinglelc;
  int ascii;
  int quiet_mode;

  int *Ncopies;
  int *copy_cnum;
  int Ncopiestotal;
  int Ncopycommands;
  int *copycommand_index;
  int *copy_origlc_index;
  int *is_lc_ready;
  int *start_cnum;

  char **lcnames;
  char lclistname[MAXLEN];
  int sizecommandvector;
  int Ncol;
  int *col;
  void ***colptr;
  int *coltype;
  int decorrflag;
  int headeronly;
  int inputlistformat;
  int showinputlcformat;
  int randseed;
  int numbercolumns;
  int oneline;
  int *col_commandstart;
  int *col_commandstop;
  int max_colcommand;
  int size_colcommandvec;
  OutColumn *outcolumns;
  int Ncolumns;
  int Ncolstolink;
  char **colnamestolink;
  OutColumn ***outcolumnstolink;
  int *columnstolink_cmdidx;
  int redirectstats;
  int redirectstatsappend;
  char redirectstatsname[MAXLEN];

  char isifcommands;
  _IfStack **IfStack;

  char showversion;
  char logcmd;

  /* Variable related to buffering the output table */
  _StringBuffer **free_buffer_stack;
  _StringBuffer **full_buffer_stack;
  int Nbuffs_free;
  int Nbuffs_free_stack;
  int Nbuffs_full;
  int free_buffer_stack_alloclen;
  int full_buffer_stack_alloclen;

  /* Variables related to handling multi-filter observations.
     Currently nothing is done with these */
  int ismultifilt;
  int *Nfilt;
  int **filt_id;
  char ***filtname;

#ifdef PARALLEL
  int *threadsinuse;
  int Nproc;
  int Nproc_allow;
  pthread_mutex_t Nproc_mutex;
  pthread_mutex_t outfile_mutex;
  pthread_t *pth;
  char *pth_init;
  sem_t threadfree;
  pthread_mutex_t cfitsio_mutex;
  pthread_mutex_t is_lc_ready_mutex;
#endif

#ifdef DYNAMICLIB
  _UserLib *UserLib;
  int NUserLib;

  _UserFunc *UserFunc;
  int NUserFunc;

  _AnalyticUserFunc *AnalyticUserFunc;
  int NAnalyticUserFunc;

#ifdef _HAVE_PYTHON
  _VartoolsPythonLibStruct VartoolsPythonLib;
#endif
#endif

  int maxinputcolumn;
  int inputcolumn_iter_index;
  _DataFromInputList *DataFromInputList;
  int NDataFromInputList;

  int NDefinedVariables;
  _Variable **DefinedVariables;
  int NDataFromLightCurve;
  _DataFromLightCurve *DataFromLightCurve;
  int inputlccolumn_iter_index;
  int maxinputlccolumn;

  _ScalarData *ScalarData;
  int NScalarData;

  int NInternalVars;

  char skipmissing;
  char skipempty;

  int pythonlibraryloaded;

  int lcdelimtype;
  char delimchar;
  char *delimstring;

} ProgramData;

#define ERR_USAGE 1
#define ERR_FILENOTFOUND 2
#define ERR_MEMALLOC 3
#define ERR_WRONGITER 4
#define ERR_WRONGORDER 5
#define ERR_KILLHARM_NOAOV 6
#define ERR_KILLHARM_NOLS 7
#define ERR_KILLHARM_NOBOTH 8
#define ERR_KILLHARM_WRONGNPER 9
#define ERR_READFORMAT 10
#define ERR_NEEDLIST 11
#define ERR_INVALIDGLOBALDECORR 12
#define ERR_INPUTMISSINGCOLUMN 13
#define ERR_CODEERROR 14
#define ERR_NOTENOUGHSTARS_ERS 15
#define ERR_CANNOTWRITE 16
#define ERR_SINGULARMATRIX 17
#define ERR_TOOMANYAMOEBAITERATIONS 18
#define ERR_BLSNBMAX 19
#define ERR_BLSFMINTOOSMALL 20
#define ERR_BLSNOFREQ 21
#define ERR_KILLHARM_NOBLS 22
#define ERR_KILLHARM_NOBLSFIXPER 23
#define ERR_SYSREMUSEWEIGHTS 24
#define ERR_UNSORTEDLIGHTCURVE 25
#define ERR_BINSIZEZERO 26
#define ERR_GETLSAMPTHRESH_FILETOSHORT 27
#define ERR_SIGFILEWRONGLENGTH 28
#define ERR_BADTYPE 29
#define ERR_NOSTRINGIDCOLUMN 30
#define ERR_KILLHARM_NOINJECTHARM 31
#define ERR_KILLHARM_NEGATIVEPERIOD 32
#define ERR_HELP_BADCOMMAND 33
#define ERR_MISSINGSAVELC 34
#define ERR_NOCOLUMN 35
#define ERR_BADCOLUMNLINK 36
#define ERR_INVALIDOUTPUTFORMAT 37
#define ERR_CANNOTOPEN 38
#define ERR_IMAGEHDU 39
#define ERR_FITSERROR 40
#define ERR_NPROC_TOO_SMALL 41
#define ERR_OPEN_LIBRARY 42
#define ERR_EQUALTIMES 43
#define ERR_EXAMPLE_BADCOMMAND 44
#define ERR_UTCINPUT_FITS 45
#define ERR_UNKNOWNOBSERVATORY 46
#define ERR_NO_EPHEM_FILE 47
#define ERR_NO_LEAPSEC_FILE 48
#define ERR_NO_PLANETDATA_FILE 49
#define ERR_INVALID_UTC_FORMAT 50
#define ERR_CONVERTTIME_NORADEC 51
#define ERR_LIBRARY_MISSING_FUNCTION 52
#define ERR_INVALIDUSEOFPARSEFIXSPECFIXCOLUMN 53
#define ERR_INVALID_CALL_REGISTERDATAVECTOR 54
#define ERR_FUNCTIONCALL_INVALIDNEXPR 55
#define ERR_ANALYTICPARSE 56
#define ERR_INVALIDARGUMENTTOEXPR 57
#define ERR_VARIABLEALREADYINUSE 58
#define ERR_BADVARIABLENAME 59
#define ERR_TOOMANYREADFORMATINPUTLCFORMAT 60
#define ERR_WRONGTYPEFORSPECIALVARIABLE 61
#define ERR_BADCOLUMNFORMATSTRING 62
#define ERR_UNDEFINEDVARIABLE 63
#define ERR_INVALIDVARIABLEFORCHANGEVAR 64
#define ERR_INVALIDVARIABLEFORLINFIT 65
#define ERR_INVALIDFUNCTIONLINFIT 66
#define ERR_LINFITMISSINGPARAM 67
#define ERR_MYSORT_GENERIC_BADCALL 68
#define ERR_CREATEVARIABLE_OUTCOLUMN_NEEDS_VPTRINPUT 69
#define ERR_OUTCOLUMN_ON_LHS 70
#define ERR_BADVARIABLETYPE_STATSCOMMAND 71
#define ERR_INVALIDSTATISTIC 72
#define ERR_BADIFTHENELSE 73
#define ERR_TOOMANYINLISTVARS 74
#define ERR_RESERVEDVARIABLENAME 75
#define ERR_NO_PREVIOUS_COMMAND 76
#define ERR_INVALID_CALL_GETCOLUMNNAMEFORRECENTCOMMAND 77
#define ERR_BINARYLC_KEYWORDTOOLONG 78
#define ERR_BINARYLIGHTCURVE_INVALIDFORMAT 79
#define ERR_BINARYLIGHTCURVE_MISSINGDATA 80
#define ERR_BADTYPE_BINARYLC 81
#define ERR_MISSING_BINARYLC_HEADERNAME 82
#define ERR_MISSING_FITSLC_HEADERNAME 83
#define ERR_FITS_CREATETABLE 84
#define ERR_FITS_WRITECOLUMN 85
#define ERR_FILEEXISTS_NOCLOBBER 86
#define ERR_INVALIDVARIABLEFORNONLINFIT 87
#define ERR_MCMCINVALIDMAXMEM 88
#define ERR_BADNONLINFITPARAMINIT 89
#define ERR_BADNONLINFITPRIORINIT 90
#define ERR_INVALID_PARAMETERVALUE 91
#define ERR_READALL_ANDCOPYLC 92
#define ERR_NR 93
#define ERR_NEGATIVE_COVAR_PARAM 94
#define ERR_NEGATIVE_COVAR_PARAM_ADDNOISE 95
#define ERR_INTERP_NOTENOUGHPOINTS 96
#define ERR_INVALIDVARIABLELCVARIABLE 97
#define ERR_INVALIDEXECCOMMANDSTRFORMAT 98
#define ERR_INVALIDANALYTICFUNCTIONDEFINITION 99
#define ERR_ANALYTICFUNCTIONDUPLICATEINPUTARG 100
#define ERR_FUNCNAMETOOLONG 101
#define ERR_INDEXINGWRONGVARIABLETYPEINEXPRESSION 102
#define ERR_BADINDEXINGOFLHSVARIABLEINEXPRESSIONCOMMAND 103
#define ERR_FUNCTIONCALL_LENINVALIDOPERAND 104
#define ERR_BADVECTORTYPEFOROUTPUTCOLUMNVARIABLE 105
#define ERR_PYTHONOUTPUTUNDEFINEDVARIABLE 106
#define ERR_MISSINGRESTRICTTIMES 107
#define ERR_NEEDCSPICE 108
#define ERR_TFATEMPLATEFITSREADERROR 109
