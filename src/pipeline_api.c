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
/* pipeline_api.c – init-once / process-many entry points for embedding
 * vartools as an in-process library.
 *
 * Public API (declared in vartools.h):
 *   ProgramData *vartools_init_pipeline(int argc, char **argv);
 *   int          vartools_process_lc(ProgramData *, const double *, ...);
 *   void         vartools_free_pipeline(ProgramData *);
 */

#include "commands.h"
#include "programdata.h"
#include "functions.h"

typedef struct {
    const char *name;
    int         datatype;
    int         vectortype;
    const void *dataptr;
    int         length;
} VartoolsVarInfo;

#include <stdlib.h>
#include <string.h>
#include <setjmp.h>

#ifdef DYNAMICLIB
#include <ltdl.h>
#endif

/* ---------------------------------------------------------------------------
 * Internal wrapper struct.
 *
 * ProgramData MUST be the first member so that a pointer to _PipelineState
 * can safely be cast to ProgramData * and vice versa (guaranteed by the C
 * standard).  The Command array cannot live inside ProgramData without
 * adding commands.h there, so we keep it here.
 * ---------------------------------------------------------------------------
 */
typedef struct {
  ProgramData p;    /* MUST be first */
  Command    *c;
} _PipelineState;


/* ---------------------------------------------------------------------------
 * Helper: zero-initialise a ProgramData and set the few non-zero defaults
 * that main() normally sets before calling parsecommandline().
 * ---------------------------------------------------------------------------
 */
static void _init_programdata(ProgramData *p)
{
  int i;
  /* calloc() already zeros everything – just set the non-zero defaults. */
  p->ascii                         = 1;
  p->skipempty                     = 1;
  p->sizecommandvector             = INITCOMMANDSIZE;
  p->randseed                      = 1;
  p->coljd                         = 1;
  p->colmag                        = 2;
  p->colsig                        = 3;
  p->inputcolumn_iter_index        = 1;
  p->lc_getcolumnsfromheader_notyetset = 1;
  p->lcdelimtype                   = VARTOOLS_LC_DELIMTYPE_WHITESPACE;
  p->JDTOL                         = DEFAULT_JDTOL;
  p->Nbuffs_free                   = VARTOOLS_DEFAULT_NOUTPUT_BUFFERS;
#ifdef PARALLEL
  p->Nproc_allow                   = 1;
#endif
  InitLinkedList(&(p->lcs_to_proc));
}


/* ---------------------------------------------------------------------------
 * vartools_init_pipeline
 *
 * Accepts a command-line of the form:
 *   argv[0]="vartools"  argv[1..]="-rms", "-oneline", ...
 * Inserts "-i -" at position 1 so that parsecommandline configures single-
 * LC stdin mode (Nlcs=1, lcnames[0]="stdin", readfromstdinflag=1).
 * ---------------------------------------------------------------------------
 */
ProgramData *vartools_init_pipeline(int argc, char **argv)
{
  _PipelineState *s;
  int             i, new_argc;
  char          **new_argv;

  /* Allocate and zero-initialise the state block. */
  s = (_PipelineState *) calloc(1, sizeof(_PipelineState));
  if (!s) return NULL;

  _init_programdata(&s->p);
  s->p.pipeline_mode = 1;

  /* Register as the active pipeline so error()/error2() longjmp here. */
  vartools_error_set_pipeline_context(&s->p);

  /* Set up the error-recovery landing pad.  If parsecommandline or
     InitCommands calls error(), we land here with exit_code set. */
  if (setjmp(s->p.exit_jmp)) {
    vartools_error_set_pipeline_context(NULL);
    /* Partial initialisation – free what we safely can. */
    if (s->c) free(s->c);
    free(s);
    return NULL;
  }

  /* Allocate the Command vector. */
  if ((s->c = (Command *) malloc(s->p.sizecommandvector * sizeof(Command)))
      == NULL)
    vt_error(ERR_MEMALLOC);
  for (i = 0; i < s->p.sizecommandvector; i++) {
    s->c[i].require_sort             = 0;
    s->c[i].require_distinct         = 0;
    s->c[i].N_setparam_expr          = 0;
    s->c[i].N_prior_vars             = 0;
    s->c[i].prior_var_datatypes      = NULL;
    s->c[i].prior_var_vectortypes    = NULL;
    s->c[i].prior_var_names          = NULL;
    s->c[i].prior_vars               = NULL;
    s->c[i].command_outcolumn_suffix = NULL;
  }

#ifdef DYNAMICLIB
  if (lt_dlinit()) {
    vt_error2(ERR_OPEN_LIBRARY,
           "Cannot initialize libtool for opening a library.\n");
  }
#ifdef VARTOOLSLIB_USERLIBSDIR
  if (lt_dladdsearchdir(VARTOOLSLIB_USERLIBSDIR)) {
    lt_dlerror();
    vt_error(ERR_OPEN_LIBRARY);
  }
#endif
#ifdef VARTOOLSLIB_USERFUNCSDIR
  if (lt_dladdsearchdir(VARTOOLSLIB_USERFUNCSDIR)) {
    lt_dlerror();
    vt_error(ERR_OPEN_LIBRARY);
  }
#endif
#endif /* DYNAMICLIB */

  /* Build a new argv with "-i -" inserted after argv[0] ("vartools").
   * This causes parsecommandline to set fileflag=1, Nlcs=1,
   * lcnames[0]="stdin", readfromstdinflag=1. */
  new_argc = argc + 2;
  new_argv = (char **) malloc((new_argc + 1) * sizeof(char *));
  if (!new_argv) {
    vt_error(ERR_MEMALLOC);
    return NULL;  /* unreachable, but keeps static analysis happy */
  }
  new_argv[0] = argc > 0 ? argv[0] : (char *)"vartools";
  new_argv[1] = (char *)"-i";
  new_argv[2] = (char *)"-";
  for (i = 1; i < argc; i++)
    new_argv[i + 2] = argv[i];
  new_argv[new_argc] = NULL;

  parsecommandline(new_argc, new_argv, &s->p, &s->c);
  free(new_argv);

  srand(s->p.randseed);
  InitCommands(&s->p, s->c);
  ReadGlobalDecorr(&s->p, s->c);
  ReadDatesFiles(&s->p, s->c);

  /* Set up the per-thread LC data structure (single thread). */
  InitializeMemAllocDataFromLightCurve(&s->p, s->c, 1);

  return (ProgramData *) s;
}


/* ---------------------------------------------------------------------------
 * vartools_process_lc
 *
 * Inject one light curve and run the command pipeline.  The output is
 * written in vartools -oneline format to outbuf (caller-provided,
 * NUL-terminated).
 *
 * Returns 0 on success, non-zero error code on failure.
 * ---------------------------------------------------------------------------
 */
int vartools_process_lc(ProgramData       *p,
                        const double      *t,
                        const double      *mag,
                        const double      *err,
                        int                n,
                        const char        *lc_name,
                        char              *outbuf,
                        int                outbuf_size)
{
  _PipelineState *s = (_PipelineState *) p;
  int    i;
  char  *printbuf      = NULL;
  int    printbuf_size = 0;

  /* Re-arm the error context so error() / error2() land back here. */
  vartools_error_set_pipeline_context(p);
  if (setjmp(p->exit_jmp)) {
    /* Re-arm for the next call. */
    vartools_error_set_pipeline_context(p);
    if (printbuf) free(printbuf);
    return p->exit_code;
  }

  /* ---- Reset per-call state ---- */
  p->skipfaillc[0] = 0;

  /* Grow the LC data arrays if this light curve is larger than any seen
   * so far.  MemAllocDataFromLightCurve is a no-op when already big enough. */
  MemAllocDataFromLightCurve(p, 0, n);

  /* Point p->t[0], p->mag[0], p->sig[0] at the allocated storage. */
  SetTimeMagSigPointers(p, 0);

  /* Inject the light curve data. */
  memcpy(p->t[0],   t,   (size_t)n * sizeof(double));
  memcpy(p->mag[0], mag, (size_t)n * sizeof(double));
  memcpy(p->sig[0], err, (size_t)n * sizeof(double));
  p->NJD[0] = n;

  /* Set the LC name (overwrite the "stdin" placeholder set by parsecommandline). */
  if (lc_name != NULL) {
    strncpy(p->lcnames[0], lc_name, MAXLEN - 1);
    p->lcnames[0][MAXLEN - 1] = '\0';
  } else {
    strncpy(p->lcnames[0], "lc", MAXLEN - 1);
  }

  /* Reset if-command stack if the pipeline uses -if/-fi. */
  if (p->isifcommands) {
    while (p->IfStack[0]->curpos > 0)
      popIfStack(p->IfStack[0]);
    p->IfStack[0]->curpos = 0;
  }

  Reset_outlc_fitsheader_additions(p, 0);

  /* ---- Run commands ---- */
  for (i = 0; i < p->Ncommands; i++) {
    if (!p->skipfaillc[0]) {
      ProcessCommandSingle(p, &(s->c[i]), 0, i, 0);
    } else {
      int ii;
      if (p->Ncopycommands > 0)
        turnoffcopies(p, s->c, i, 0, 0);
      for (ii = i; ii < p->Ncommands; ii++) {
        if (s->c[ii].cnum != CNUM_SAVELC && s->c[ii].cnum != CNUM_RESTORELC)
          SkipCommand(p, &(s->c[ii]), ii, 0, 0);
      }
      break;
    }
  }

  /* ---- Collect -oneline output into caller's buffer ---- */
  printresults_buffer_new(p, 0, 0, &printbuf, &printbuf_size);

  if (outbuf != NULL && outbuf_size > 0) {
    if (printbuf != NULL) {
      strncpy(outbuf, printbuf, (size_t)(outbuf_size - 1));
      outbuf[outbuf_size - 1] = '\0';
    } else {
      outbuf[0] = '\0';
    }
  }
  free(printbuf);

  return 0;
}


/* ---------------------------------------------------------------------------
 * vartools_free_pipeline
 *
 * Release all resources held by the pipeline handle.
 * ---------------------------------------------------------------------------
 */
void vartools_free_pipeline(ProgramData *p)
{
  _PipelineState *s;
  int i;

  if (!p) return;
  s = (_PipelineState *) p;

  /* Clear the global error context. */
  vartools_error_set_pipeline_context(NULL);

  /* Run any command-specific teardown (e.g. user-defined commands). */
  for (i = 0; i < p->Ncommands; i++) {
    if (s->c[i].cnum == CNUM_USERCOMMAND)
      CloseUserCommand(p, &(s->c[i]));
  }

  CloseTrackedOpenFiles(p);

#ifdef _HAVE_PYTHON
  KillAllPythonProcesses(p, s->c);
#endif
#ifdef _HAVE_R
  KillAllRProcesses(p, s->c);
#endif

  free(s->c);
  free(s);
}


/* ---------------------------------------------------------------------------
 * vartools_get_lc_variables
 *
 * After vartools_process_lc() has run, enumerate all light-curve variables
 * and per-star scalars.  The caller provides an array of VartoolsVarInfo
 * structs; this function fills them with metadata and a pointer to the
 * raw data array for thread 0.
 *
 * The data pointers are valid until the next vartools_process_lc() call
 * or vartools_free_pipeline().  The caller should copy any data it needs
 * before the next call.
 *
 * Returns 0 on success.  *n_vars is set to the number of entries written.
 * If max_vars is too small, only the first max_vars entries are written
 * and the function returns the total number of variables available.
 * ---------------------------------------------------------------------------
 */
int vartools_get_lc_variables(ProgramData  *p,
                              int          *n_vars,
                              VartoolsVarInfo *vars,
                              int           max_vars)
{
  int i, count = 0;
  _Variable *v;

  for (i = 0; i < p->NDefinedVariables; i++) {
    v = p->DefinedVariables[i];
    if (v == NULL) continue;

    if (v->vectortype == VARTOOLS_VECTORTYPE_LC) {
      if (count < max_vars) {
        vars[count].name = v->varname;
        vars[count].datatype = v->datatype;
        vars[count].vectortype = v->vectortype;
        vars[count].length = p->NJD[0];
        switch (v->datatype) {
        case VARTOOLS_TYPE_DOUBLE:
          vars[count].dataptr = (*((double ***) v->dataptr))[0];
          break;
        case VARTOOLS_TYPE_FLOAT:
          vars[count].dataptr = (*((float ***) v->dataptr))[0];
          break;
        case VARTOOLS_TYPE_INT:
          vars[count].dataptr = (*((int ***) v->dataptr))[0];
          break;
        case VARTOOLS_TYPE_LONG:
          vars[count].dataptr = (*((long ***) v->dataptr))[0];
          break;
        case VARTOOLS_TYPE_SHORT:
          vars[count].dataptr = (*((short ***) v->dataptr))[0];
          break;
        default:
          vars[count].dataptr = NULL;
          break;
        }
      }
      count++;
    }
    else if (v->vectortype == VARTOOLS_VECTORTYPE_SCALAR ||
             v->vectortype == VARTOOLS_VECTORTYPE_PERSTARDATA ||
             v->vectortype == VARTOOLS_VECTORTYPE_INLIST) {
      if (count < max_vars) {
        vars[count].name = v->varname;
        vars[count].datatype = v->datatype;
        vars[count].vectortype = v->vectortype;
        vars[count].length = 1;
        switch (v->datatype) {
        case VARTOOLS_TYPE_DOUBLE:
          vars[count].dataptr = &((*((double **) v->dataptr))[0]);
          break;
        case VARTOOLS_TYPE_INT:
          vars[count].dataptr = &((*((int **) v->dataptr))[0]);
          break;
        case VARTOOLS_TYPE_LONG:
          vars[count].dataptr = &((*((long **) v->dataptr))[0]);
          break;
        default:
          vars[count].dataptr = NULL;
          break;
        }
      }
      count++;
    }
  }

  *n_vars = count < max_vars ? count : count;
  return count > max_vars ? count : 0;
}


/* ---------------------------------------------------------------------------
 * vartools_get_njd
 *
 * Return the number of points in the light curve after processing
 * (may differ from the input if clipping, binning, etc. were applied).
 * ---------------------------------------------------------------------------
 */
int vartools_get_njd(ProgramData *p)
{
  return p->NJD[0];
}


/* ---------------------------------------------------------------------------
 * vartools_set_lc_data
 *
 * Inject a full set of light-curve columns into the pipeline before
 * calling vartools_process_lc().  This extends the basic t/mag/err
 * injection to include additional named columns.
 *
 * col_names[0..n_columns-1] are the variable names.  The first three
 * must be "t", "mag", "err" (or whatever the pipeline expects).
 * col_data[i] points to n_points doubles for column i.
 *
 * For columns beyond t/mag/err, the function looks up each name in
 * p->DefinedVariables.  If found and the variable is VECTORTYPE_LC,
 * it copies the data in.  If not found, it is silently skipped
 * (the variable does not exist in this pipeline's command set).
 *
 * Returns 0 on success, -1 if n_columns < 3.
 * ---------------------------------------------------------------------------
 */
int vartools_set_lc_data(ProgramData     *p,
                         int              n_points,
                         int              n_columns,
                         const char     **col_names,
                         const double   **col_data,
                         const char      *lc_name)
{
  int i, j;
  _Variable *v;
  double **arr;

  if (n_columns < 3) return -1;

  /* Grow internal arrays if needed */
  MemAllocDataFromLightCurve(p, 0, n_points);
  SetTimeMagSigPointers(p, 0);

  /* Inject t, mag, err (first three columns) */
  memcpy(p->t[0],   col_data[0], (size_t)n_points * sizeof(double));
  memcpy(p->mag[0], col_data[1], (size_t)n_points * sizeof(double));
  memcpy(p->sig[0], col_data[2], (size_t)n_points * sizeof(double));
  p->NJD[0] = n_points;

  /* Inject additional columns by variable name lookup */
  for (i = 3; i < n_columns; i++) {
    for (j = 0; j < p->NDefinedVariables; j++) {
      v = p->DefinedVariables[j];
      if (v == NULL) continue;
      if (v->vectortype != VARTOOLS_VECTORTYPE_LC) continue;
      if (v->datatype != VARTOOLS_TYPE_DOUBLE) continue;
      if (strcmp(v->varname, col_names[i]) != 0) continue;
      arr = *((double ***) v->dataptr);
      memcpy(arr[0], col_data[i], (size_t)n_points * sizeof(double));
      break;
    }
  }

  /* Set the LC name */
  if (lc_name != NULL) {
    strncpy(p->lcnames[0], lc_name, MAXLEN - 1);
    p->lcnames[0][MAXLEN - 1] = '\0';
  } else {
    strncpy(p->lcnames[0], "lc", MAXLEN - 1);
  }

  return 0;
}
