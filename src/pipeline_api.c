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
  p->Nbuffs_free_user_set          = 0;
  p->setlcname                     = NULL;
  p->Ncaptured                     = 0;
  p->Ncaptured_filled              = 0;
  p->captured                      = NULL;
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

  /* Pre-walk s->c[] to count cmd.o instances with capture_to_buffer set,
     and allocate the captured-LC slots once for the lifetime of the
     pipeline.  Slot indices match the order of CNUM_OUTPUTLCS commands
     in the command sequence.

     Thread-safety note: ``p->captured`` holds one set of buffers per
     capture id, not per thread.  vartools_process_lc is documented as
     single-threaded (one calling thread at a time per pipeline handle),
     and that's how libvartoolspipeline drives it today.  If someone
     ever wraps the API in a multi-threaded host that calls
     vartools_process_lc concurrently on the same handle, the slots
     would need to grow into a [thread][id] indexing scheme.  The
     standalone CLI path (no library init) leaves p->captured == NULL,
     so vartools_capture_current_lc early-returns and "-o ID capture"
     becomes a silent no-op -- safe even with -parallel N because each
     thread independently reads the NULL guard. */
  {
    int j, k;
    s->p.Ncaptured = 0;
    for (j = 0; j < s->p.Ncommands; j++) {
      if (s->c[j].cnum == CNUM_OUTPUTLCS &&
          s->c[j].Outputlcs->capture_to_buffer) {
        s->p.Ncaptured++;
      }
    }
    if (s->p.Ncaptured > 0) {
      s->p.captured = (_CapturedLC *) malloc(
          s->p.Ncaptured * sizeof(_CapturedLC));
      if (!s->p.captured) vt_error(ERR_MEMALLOC);
      for (j = 0; j < s->p.Ncaptured; j++) {
        s->p.captured[j].id[0]      = '\0';
        s->p.captured[j].filled     = 0;
        s->p.captured[j].njd        = 0;
        s->p.captured[j].n_vars     = 0;
        s->p.captured[j].varnames   = NULL;
        s->p.captured[j].datatypes  = NULL;
        s->p.captured[j].databufs   = NULL;
      }
      /* Pre-fill ids from the command sequence so a duplicate-id check
         can fire at init time (matching the Python wrapper, but also
         protecting direct CLI use of -o <id> capture). */
      k = 0;
      for (j = 0; j < s->p.Ncommands; j++) {
        if (s->c[j].cnum == CNUM_OUTPUTLCS &&
            s->c[j].Outputlcs->capture_to_buffer) {
          int kk;
          strncpy(s->p.captured[k].id,
                  s->c[j].Outputlcs->capture_id,
                  sizeof(s->p.captured[k].id) - 1);
          s->p.captured[k].id[sizeof(s->p.captured[k].id) - 1] = '\0';
          for (kk = 0; kk < k; kk++) {
            if (!strcmp(s->p.captured[kk].id, s->p.captured[k].id)) {
              fprintf(stderr,
                  "Error: -o capture id '%s' used by more than one -o "
                  "command in the pipeline.  Each capture id must be "
                  "unique.\n", s->p.captured[k].id);
              vt_error(ERR_USAGE);
            }
          }
          k++;
        }
      }
    }
  }

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

  /* Free any LC snapshots from the previous vartools_process_lc() call so
     captures from this run can be written into the same slots.  Memory
     usage stays bounded by one LC's worth of snapshots; nothing
     accumulates across batch iterations. */
  if (p->captured && p->Ncaptured > 0) {
    int j;
    for (j = 0; j < p->Ncaptured; j++) {
      _CapturedLC *cap = &p->captured[j];
      if (cap->databufs) {
        int v;
        for (v = 0; v < cap->n_vars; v++) {
          if (cap->databufs[v]) {
            if (cap->datatypes[v] == VARTOOLS_TYPE_STRING) {
              /* For strings the buffer is a malloc'd char ** of length njd
                 with one strdup'd entry per row. */
              char **strs = (char **) cap->databufs[v];
              int r;
              for (r = 0; r < cap->njd; r++) {
                if (strs[r]) free(strs[r]);
              }
            }
            free(cap->databufs[v]);
          }
        }
        free(cap->databufs);
        cap->databufs = NULL;
      }
      if (cap->varnames)  { free(cap->varnames);  cap->varnames  = NULL; }
      if (cap->datatypes) { free(cap->datatypes); cap->datatypes = NULL; }
      cap->n_vars = 0;
      cap->njd    = 0;
      cap->filled = 0;
    }
    p->Ncaptured_filled = 0;
  }

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

  /* Fire any -inputlcformat col=0 init expressions (e.g. phase:0:double:NR).
   * In file-reading mode these run inside the LC reader; the in-process
   * path skips that reader so we evaluate them here against the freshly
   * injected arrays. */
  EvaluateInputLCExpressions(p, 0, 0, 0);

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
 * vartools_capture_current_lc
 *
 * Internal helper invoked by DoOutputLightCurve when the
 * capture_to_buffer flag is set.  Walks the variable registry and
 * snapshots every VECTORTYPE_LC var into the matching p->captured[] slot
 * (looked up by id).  The snapshots are owned by the slot and freed at
 * the start of the next vartools_process_lc() call.
 *
 * Strings (VARTOOLS_TYPE_STRING) are deep-copied: the slot owns a
 * malloc'd char ** of length njd with one strdup'd entry per row.
 * Numeric types are flat-arrayed with one malloc'd buffer per var.
 *
 * Returns 0 on success, non-zero if the id is not found in the slot
 * table (which should never happen because slots are pre-populated at
 * pipeline init).
 * ---------------------------------------------------------------------------
 */
int vartools_capture_current_lc(ProgramData *p, const char *id)
{
  int slot, i, count, v, r;
  _Variable *var;
  _CapturedLC *cap;
  size_t bytes;
  int njd;

  if (!p->captured || p->Ncaptured == 0) return 1;

  /* Locate the slot by id (linear scan; Ncaptured is typically a
     small handful). */
  for (slot = 0; slot < p->Ncaptured; slot++) {
    if (!strcmp(p->captured[slot].id, id)) break;
  }
  if (slot == p->Ncaptured) return 1;

  cap = &p->captured[slot];
  njd = p->NJD[0];

  /* Count the LC variables we'll snapshot. */
  count = 0;
  for (i = 0; i < p->NDefinedVariables; i++) {
    var = p->DefinedVariables[i];
    if (var && var->vectortype == VARTOOLS_VECTORTYPE_LC) count++;
  }

  cap->njd       = njd;
  cap->n_vars    = count;
  cap->varnames  = (char **) malloc(count * sizeof(char *));
  cap->datatypes = (int *)   malloc(count * sizeof(int));
  cap->databufs  = (void **) malloc(count * sizeof(void *));
  if (!cap->varnames || !cap->datatypes || !cap->databufs)
    vt_error(ERR_MEMALLOC);
  for (v = 0; v < count; v++) cap->databufs[v] = NULL;

  v = 0;
  for (i = 0; i < p->NDefinedVariables; i++) {
    var = p->DefinedVariables[i];
    if (!var || var->vectortype != VARTOOLS_VECTORTYPE_LC) continue;
    cap->varnames[v]  = var->varname;     /* not owned */
    cap->datatypes[v] = var->datatype;
    switch (var->datatype) {
    case VARTOOLS_TYPE_DOUBLE:
      bytes = (size_t)njd * sizeof(double);
      cap->databufs[v] = malloc(bytes);
      if (!cap->databufs[v]) vt_error(ERR_MEMALLOC);
      memcpy(cap->databufs[v], (*((double ***) var->dataptr))[0], bytes);
      break;
    case VARTOOLS_TYPE_FLOAT:
      bytes = (size_t)njd * sizeof(float);
      cap->databufs[v] = malloc(bytes);
      if (!cap->databufs[v]) vt_error(ERR_MEMALLOC);
      memcpy(cap->databufs[v], (*((float ***) var->dataptr))[0], bytes);
      break;
    case VARTOOLS_TYPE_INT:
      bytes = (size_t)njd * sizeof(int);
      cap->databufs[v] = malloc(bytes);
      if (!cap->databufs[v]) vt_error(ERR_MEMALLOC);
      memcpy(cap->databufs[v], (*((int ***) var->dataptr))[0], bytes);
      break;
    case VARTOOLS_TYPE_LONG:
      bytes = (size_t)njd * sizeof(long);
      cap->databufs[v] = malloc(bytes);
      if (!cap->databufs[v]) vt_error(ERR_MEMALLOC);
      memcpy(cap->databufs[v], (*((long ***) var->dataptr))[0], bytes);
      break;
    case VARTOOLS_TYPE_SHORT:
      bytes = (size_t)njd * sizeof(short);
      cap->databufs[v] = malloc(bytes);
      if (!cap->databufs[v]) vt_error(ERR_MEMALLOC);
      memcpy(cap->databufs[v], (*((short ***) var->dataptr))[0], bytes);
      break;
    case VARTOOLS_TYPE_STRING: {
      /* Vartools stores strings as char ***[thread][row] -> char *.
         Snapshot is a malloc'd char ** of length njd, each entry
         strdup'd so subsequent in-place edits don't disturb it. */
      char **src = (*((char ****) var->dataptr))[0];
      char **dst = (char **) malloc(njd * sizeof(char *));
      if (!dst) vt_error(ERR_MEMALLOC);
      for (r = 0; r < njd; r++) {
        dst[r] = src[r] ? strdup(src[r]) : NULL;
      }
      cap->databufs[v] = (void *) dst;
      break;
    }
    default:
      cap->databufs[v] = NULL;   /* unsupported type — leave NULL */
      break;
    }
    v++;
  }

  cap->filled = 1;
  p->Ncaptured_filled++;
  return 0;
}


/* ---------------------------------------------------------------------------
 * vartools_get_captured_lc
 *
 * Read back a previously-captured LC snapshot by id, populating *vars*
 * with up to *max_vars* entries (one per LC variable).  *n_vars is set
 * to the number of entries written; the function returns 0 on success
 * or the total count if the caller's array was too small.  Returns -1
 * if the id is not found, or -2 if the slot was never filled.
 *
 * The dataptrs in *vars* point into the captured buffers and are valid
 * until the next vartools_process_lc() call (same lifetime contract as
 * vartools_get_lc_variables).  Callers must memcpy out anything they
 * want to keep across calls.
 * ---------------------------------------------------------------------------
 */
int vartools_get_captured_lc(ProgramData     *p,
                             const char      *id,
                             VartoolsVarInfo *vars,
                             int              max_vars,
                             int             *n_vars)
{
  int slot, v, count;
  _CapturedLC *cap;

  if (n_vars) *n_vars = 0;
  if (!p->captured || p->Ncaptured == 0) return -1;

  for (slot = 0; slot < p->Ncaptured; slot++) {
    if (!strcmp(p->captured[slot].id, id)) break;
  }
  if (slot == p->Ncaptured) return -1;

  cap = &p->captured[slot];
  if (!cap->filled) return -2;

  count = 0;
  for (v = 0; v < cap->n_vars; v++) {
    if (count < max_vars) {
      vars[count].name       = cap->varnames[v];
      vars[count].datatype   = cap->datatypes[v];
      vars[count].vectortype = VARTOOLS_VECTORTYPE_LC;
      vars[count].length     = cap->njd;
      vars[count].dataptr    = cap->databufs[v];
    }
    count++;
  }

  if (n_vars) *n_vars = count < max_vars ? count : max_vars;
  return count > max_vars ? count : 0;
}


/* ---------------------------------------------------------------------------
 * vartools_get_captured_njd
 *
 * Return the number of points in a captured LC snapshot, or -1 if id
 * is not found / -2 if the slot was never filled.
 * ---------------------------------------------------------------------------
 */
int vartools_get_captured_njd(ProgramData *p, const char *id)
{
  int slot;
  if (!p->captured || p->Ncaptured == 0) return -1;
  for (slot = 0; slot < p->Ncaptured; slot++) {
    if (!strcmp(p->captured[slot].id, id)) break;
  }
  if (slot == p->Ncaptured) return -1;
  return p->captured[slot].filled ? p->captured[slot].njd : -2;
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


/* ---------------------------------------------------------------------------
 * vartools_set_inlist_*
 *
 * Update the value of an INLIST variable (declared via -inlistvars) in the
 * pipeline's single-LC storage slot between successive vartools_process_lc()
 * calls.  Library mode allocates one inlist slot per variable at init time
 * (Nlcs=1 from the injected -i -); these setters overwrite that slot before
 * the next process_lc call so each LC sees its own per-LC value.
 *
 * Return values:
 *   0 on success
 *  -1 if the named variable doesn't exist or isn't an INLIST variable
 *  -2 if the variable exists but has a different datatype than the setter
 * --------------------------------------------------------------------------- */

/* Look up INLIST storage by name.  Searches two registration paths:
 *   1. -inlistvars declares CreateVariable-registered variables; the lookup
 *      matches v->varname directly.
 *   2. cmd.o(namefromlist) and similar auto-register a DataFromInputList
 *      entry WITHOUT going through CreateVariable; the lookup matches
 *      d->incolumn_names[0] (which has a "_<cnum>" suffix appended).
 * Returns 0 on success (populates *out_dataptr and *out_datatype); -1 if no
 * matching INLIST entry exists. */
static int _find_inlist_storage(ProgramData *p, const char *name,
                                 void **out_dataptr, int *out_datatype)
{
  int i;
  if (name == NULL) return -1;
  for (i = 0; i < p->NDefinedVariables; i++) {
    _Variable *v = p->DefinedVariables[i];
    if (v == NULL) continue;
    if (v->vectortype != VARTOOLS_VECTORTYPE_INLIST) continue;
    if (v->varname == NULL) continue;
    if (strcmp(v->varname, name) == 0) {
      *out_dataptr = v->dataptr;
      *out_datatype = v->datatype;
      return 0;
    }
  }
  for (i = 0; i < p->NDataFromInputList; i++) {
    _DataFromInputList *d = &p->DataFromInputList[i];
    if (d->Ncolumns != 0) continue;
    if (d->incolumn_names == NULL || d->incolumn_names[0] == NULL) continue;
    if (strcmp(d->incolumn_names[0], name) == 0) {
      *out_dataptr = d->dataptr;
      *out_datatype = d->datatype;
      return 0;
    }
  }
  return -1;
}

int vartools_set_inlist_double(ProgramData *p, const char *name, double value)
{
  void *dataptr; int datatype;
  if (_find_inlist_storage(p, name, &dataptr, &datatype) != 0) return -1;
  if (datatype != VARTOOLS_TYPE_DOUBLE) return -2;
  (*((double **) dataptr))[0] = value;
  return 0;
}

int vartools_set_inlist_int(ProgramData *p, const char *name, int value)
{
  void *dataptr; int datatype;
  if (_find_inlist_storage(p, name, &dataptr, &datatype) != 0) return -1;
  if (datatype != VARTOOLS_TYPE_INT) return -2;
  (*((int **) dataptr))[0] = value;
  return 0;
}

int vartools_set_inlist_long(ProgramData *p, const char *name, long value)
{
  void *dataptr; int datatype;
  if (_find_inlist_storage(p, name, &dataptr, &datatype) != 0) return -1;
  if (datatype != VARTOOLS_TYPE_LONG) return -2;
  (*((long **) dataptr))[0] = value;
  return 0;
}

int vartools_set_inlist_short(ProgramData *p, const char *name, short value)
{
  void *dataptr; int datatype;
  if (_find_inlist_storage(p, name, &dataptr, &datatype) != 0) return -1;
  if (datatype != VARTOOLS_TYPE_SHORT) return -2;
  (*((short **) dataptr))[0] = value;
  return 0;
}

int vartools_set_inlist_string(ProgramData *p, const char *name, const char *value)
{
  void *dataptr; int datatype;
  if (_find_inlist_storage(p, name, &dataptr, &datatype) != 0) return -1;
  if (datatype != VARTOOLS_TYPE_STRING) return -2;
  if (value == NULL) value = "";
  char *slot = (*((char ***) dataptr))[0];
  strncpy(slot, value, MAXLEN - 1);
  slot[MAXLEN - 1] = '\0';
  return 0;
}
