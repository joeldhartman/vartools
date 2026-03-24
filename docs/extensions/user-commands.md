# Writing User Extension Commands

Extension commands let you add entirely new operations to VARTOOLS without
modifying the core source.  Each command is compiled as a shared-object library
and loaded at runtime with the `-L` flag:

```bash
vartools -L /path/to/mycommand.so -mycommand arg1 arg2 ...
```

If the library is installed into the VARTOOLS `userlibs` data directory (which
`make install` does automatically when the library is listed in
`USERLIBS/src/Makefile.am`), the `-L` flag is optional:

```bash
vartools -mycommand arg1 arg2 ...
```

The simplest working example is `magadd`, which adds a constant to every
magnitude in a light curve.  Its source (`USERLIBS/src/magadd.c` and
`magadd.h`) is the canonical starting point.

---

## Prerequisites

- C (or any language callable from C).  Fortran is also supported — see
  `jktebop` and `macula` in `USERLIBS/src/` for worked examples.
- Familiarity with C pointers.
- The GNU autotools (`autoconf`, `automake`, `libtool`) for the recommended
  build path.  On Ubuntu: `sudo apt-get install libtool automake autoconf
  autotools-dev`.  On macOS: `brew install libtool automake autoconf`.

---

## Structure of a command library

A command library exposes exactly **five** functions to VARTOOLS.  Using
`mylib` as the library name throughout:

| Function | Purpose |
|----------|---------|
| `mylib_Initialize` | Tell VARTOOLS the command name and data-structure size |
| `mylib_ParseCL` | Validate syntax and register data vectors |
| `mylib_ShowSyntax` | One-line syntax string for `vartools -help` |
| `mylib_ShowHelp` | Verbose help text for `vartools -help mylib` |
| `mylib_RunCommand` | Execute the command on one light curve |

The library name must match the command name: library `mylib.so` provides
command `-mylib`.

### Recommended file layout

```
mylib.h     ← data structure definition
mylib.c     ← the five required functions + algorithm implementation
```

At the top of `mylib.c`:

```c
#include "../../src/vartools.h"   /* adjust path as needed */
#include <stdio.h>
#include <stdlib.h>
#include "mylib.h"
```

---

## Step 1 — Define the data structure (`mylib.h`)

All per-light-curve parameters and computed results are stored in a single
struct.  VARTOOLS allocates the memory; the struct is zeroed before use.

```c
/* mylib.h */
typedef struct {
    double *result;      /* one computed value per light curve */
    double *inputparam;  /* one input parameter per light curve */
} _Mylib;
```

---

## Step 2 — `mylib_Initialize`

```c
void mylib_Initialize(char *commandname,
                      int  *RequireReadAll,
                      int  *RequireSortLC,
                      int  *RequireDistinctTimes,
                      size_t *sizeuserdata)
{
    sprintf(commandname, "-mylib");   /* command the user types */
    *RequireReadAll        = 0;       /* 1 = need all LCs in memory at once */
    *RequireSortLC         = 0;       /* 1 = input LC must be time-sorted */
    *RequireDistinctTimes  = 0;       /* 1 = observations must have unique times */
    *sizeuserdata          = sizeof(_Mylib);
}
```

Set `RequireReadAll = 1` only if the algorithm genuinely needs simultaneous
access to all light curves (e.g., an ensemble normalisation).  Sequential
processing (`RequireReadAll = 0`) is faster and compatible with `-parallel`.

---

## Step 3 — `mylib_ParseCL`

This function is called when VARTOOLS encounters `-mylib` on the command line.
It must:

1. Parse and validate every argument.
2. Register all data vectors with VARTOOLS (so memory is allocated and the
   output table is populated automatically).
3. Return `0` on success, `1` on error (VARTOOLS will then print the syntax and
   exit).

```c
int mylib_ParseCL(ProgramData *p, Command *c, void *userdata,
                  int *iret, char **argv, int argc)
{
    _Mylib *mylib = (_Mylib *) userdata;

    /* Parse a "fix|list|fixcolumn|expr" parameter called "inputparam" */
    if (VARTOOLS_ParseFixSpecFixcolumn(p, c, iret, argv, argc, 1,
            VARTOOLS_TYPE_DOUBLE,
            (void *) (&mylib->inputparam),
            0,        /* Ncolumns: 0 = scalar per LC */
            0,        /* output: 0 = don't print in table */
            "inputparam"))
        return 1;

    /* Register a computed output value to appear in the output table */
    VARTOOLS_RegisterDataVector(p, c,
            (void *) (&mylib->result),
            VARTOOLS_TYPE_DOUBLE,
            0,                          /* Ncolumns */
            VARTOOLS_SOURCE_COMPUTED,   /* source */
            1,                          /* output: 1 = include in table */
            "MyResult");                /* column heading root */

    return 0;
}
```

`argv[*iret]` on entry is the command name (`"-mylib"`); the API functions
advance `*iret` past any tokens they consume.

---

## Step 4 — `mylib_ShowSyntax` and `mylib_ShowHelp`

```c
void mylib_ShowSyntax(FILE *outfile)
{
    fprintf(outfile,
        "-mylib <\"fix\" value | \"list\" [\"column\" col] "
        "| \"fixcolumn\" <colname|colnum> | \"expr\" expr>\n");
}

void mylib_ShowHelp(FILE *outfile)
{
    fprintf(outfile,
        "Compute something useful for each light curve.\n\n"
        "The single required argument provides the input parameter.  "
        "Use \"fix\" to supply a constant, \"list\" to read per-LC values "
        "from the input list file, \"fixcolumn\" to take the value from a "
        "prior command's output, or \"expr\" to evaluate an expression.\n\n");
}
```

---

## Step 5 — `mylib_RunCommand`

This function does the actual work.  It receives one light curve at a time
(when `RequireReadAll = 0`).

```c
void mylib_RunCommand(ProgramData *p, void *userdata,
                      int lc_name_num, int lc_num)
{
    _Mylib *mylib = (_Mylib *) userdata;

    int     N      = p->NJD[lc_num];     /* number of observations */
    double *t      = p->t[lc_num];       /* times of observation   */
    double *mag    = p->mag[lc_num];     /* magnitudes             */
    double *err    = p->sig[lc_num];     /* uncertainties          */
    char   *lcname = p->lcnames[lc_name_num];

    double param   = mylib->inputparam[lc_num];   /* registered input */

    /* --- run your algorithm here --- */
    double computed_result = 0.0;
    for (int i = 0; i < N; i++)
        computed_result += mag[i];
    computed_result /= N;   /* e.g. naive mean magnitude */

    mylib->result[lc_num] = computed_result;   /* stored → appears in table */
}
```

!!! warning "Thread safety"
    If you use `-parallel N`, `mylib_RunCommand` is called simultaneously from
    multiple threads.  **Never use static or global variables** inside this
    function or any function it calls.  All mutable state must live in the
    `userdata` struct, indexed by `lc_num`.

When `RequireReadAll = 1`, `lc_num` and `lc_name_num` are both `-1`.  Access
all light curves through `p->Nlcs`, `p->NJD[i]`, `p->t[i][j]`, etc.

---

## Registering data vectors

`VARTOOLS_RegisterDataVector` is the central API call.  Every variable that
holds per-LC data should be registered so VARTOOLS handles allocation,
input/output, and table formatting automatically.

```c
void VARTOOLS_RegisterDataVector(
    ProgramData *p,
    Command     *c,
    void        *dataptr,    /* e.g. (void *)(&mylib->result)  */
    int          datatype,   /* VARTOOLS_TYPE_DOUBLE, _INT, _STRING, … */
    int          Ncolumns,   /* 0 = scalar per LC; >0 = vector per LC */
    int          source,     /* see table below */
    int          output,     /* 1 = include in output ASCII table */
    char        *outname,    /* column heading root */
    ...                      /* source-specific extra args */
);
```

**Source constants:**

| Constant | Meaning | Extra arguments |
|----------|---------|-----------------|
| `VARTOOLS_SOURCE_COMPUTED` | Algorithm fills it in at run time | *(none)* |
| `VARTOOLS_SOURCE_INLIST` | Read one value per LC from the input list | `const char *name`, `int colnum` |
| `VARTOOLS_SOURCE_FIXED` | Same constant for every LC (set from the CLI) | `void *fixptr` |
| `VARTOOLS_SOURCE_PRIORCOLUMN` | Take value from a prior command's output column | `char *colname` |
| `VARTOOLS_SOURCE_LC` | Read a column from each light curve file | `const char *name`, `int colnum`, `char *scanformat`, `char *varnameout` |
| `VARTOOLS_SOURCE_EVALEXPRESSION` | Evaluate an analytic expression | `char *exprstring` |
| `VARTOOLS_SOURCE_RECENTCOMMAND` | Take value from most-recent instance of a command | `int cnum`, `char *paramname` |

---

## Parsing API helpers

These convenience functions parse common command-line patterns and automatically
call `VARTOOLS_RegisterDataVector` internally:

### `VARTOOLS_ParseFixSpecFixcolumn`

Parses: `"fix" value | "list" ["column" col] | "fixcolumn" <colname|colnum> | "expr" expression`

```c
int VARTOOLS_ParseFixSpecFixcolumn(
    ProgramData *p, Command *c, int *iret, char **argv, int argc,
    int Nvec,
    /* for each vector: */
    int datatype, void *dataptr, int Ncolumns, int output, char *name);
```

Returns `0` on success, non-zero on error.

### `VARTOOLS_ParseParameter`

Parses: `"keyword" <"fix" value | "list" [...] | "fixcolumn" ... | "expr" ...>`

```c
int VARTOOLS_ParseParameter(
    ProgramData *p, Command *c, int *iret, char **argv, int argc,
    const char *keyword, int Nvec,
    /* for each vector: */
    int datatype, void *dataptr, int Ncolumns,
    int output, char *name, int initialize, double defaultval);
```

Returns `0` = parsed, `1` = keyword absent, `2` = parse error.

### `VARTOOLS_ParseConstantParameter`

Parses: `["keyword"] value [value2 … valueN]`

```c
int VARTOOLS_ParseConstantParameter(
    ProgramData *p, Command *c, int *iret, char **argv, int argc,
    const char *keyword,   /* NULL to parse positionally */
    int datatype, void *dataptr, int Ncolumns);
```

### `VARTOOLS_ParseOutNameKeyword`

Parses: `"keyword" outdir ["nameformat" format]`

```c
int VARTOOLS_ParseOutNameKeyword(
    ProgramData *p, Command *c, int *iret, char **argv, int argc,
    const char *keyword,
    int *outputflag, char **outdir,
    int *formatflag, char **format);
```

### `VARTOOLS_GetOutputFilename`

Constructs an output filename following the standard vartools convention:

```c
void VARTOOLS_GetOutputFilename(
    char *lcoutname,   /* result stored here */
    char *lcname,      /* input light curve name */
    char *outdir,
    char *suffix,
    char *format,      /* printf-style format, or NULL */
    int   lc_name_num);
```

---

## libvartools utility functions

These functions from `libvartools` are available for use inside your command.

### Model fitting

```c
/* Downhill-simplex (Nelder-Mead) minimiser */
int VARTOOLS_amoeba(double **p, double *y, int *ia, int ndim,
    double ftol,
    double (*funk)(double *, int, int, double *, double *, double *, void *),
    int *nfunk, int maxeval,
    int N, double *t, double *mag, double *sig, void *userparams);

/* Set up the initial simplex, one call per parameter */
void VARTOOLS_incrementparameters_foramoeba(
    int *Nparameters, int *Ntovary,
    double ***p, int **ia,
    int varyparam, double initialval, double stepval);

void VARTOOLS_amoeba_cleanup(int *Nparameters, int *Ntovary,
    double ***p, int **ia);

/* Polynomial fit: y = sum_i a_i * x^i */
double VARTOOLS_fitpoly(int N, double *x, double *y, double *sig,
    int order, int subtractfit,
    double *fitparams, double *paramerrs);

/* Mandel & Agol transit model */
void VARTOOLS_mandelagoltransitmodel(int Npoints, double *phase, double *outlc,
    int type, double *ldcoeffs, double sin_i, double a, double e,
    double p, double omega);

/* As above with exposure-time integration */
void VARTOOLS_integratemandelagoltransitmodel(double exptime_phase,
    int Npoints, double *phase, double *outlc,
    int type, double *ldcoeffs, double sin_i, double a, double e,
    double p, double omega, int Nresamp);
```

### Statistics

```c
double VARTOOLS_getweightedmean(int n, double *data, double *sig);
double VARTOOLS_getmean(int n, double *data);
double VARTOOLS_median(int n, double *data);
double VARTOOLS_MAD(int n, double *data);
double VARTOOLS_stddev(int n, double *data);
double VARTOOLS_kurtosis(int n, double *data);
double VARTOOLS_skewness(int n, double *data);
double VARTOOLS_chi2(int N, double *t, double *mag, double *err,
    double *weighted_average, int *Ngood);
void   VARTOOLS_medianfilter(int N, double *t, double *mag, double *sig,
    double timesize, int meanflag, int replace);
```

### Splines

```c
void   VARTOOLS_spline(double *x, double *y, int n, double yp1, double ypn,
    double *y2, double *u);
void   VARTOOLS_splint(double *xa, double *ya, double *y2a, int n,
    double x, double *y);
void   VARTOOLS_spline_monotonic(int N, double *x, double *y, double *yprime);
double VARTOOLS_splint_monotonic(int N, double *x, double *y, double *yprime,
    double xt);
```

### Sorting

```c
/* Multi-key sort of one or more vectors */
void VARTOOLS_sort_generic(int N, int isreverse, int *index, int Nms, ...);

/* Sort a single double vector in place (low to high) */
void VARTOOLS_sortvec_double(int N, double *data);
```

### Error reporting

```c
void VARTOOLS_error(int errflag);
void VARTOOLS_error2(int errflag, char *s);
```

### Utilities

```c
/* Returns 1 if period1 and period2 are distinct given the time span */
int VARTOOLS_isDifferentPeriods(double period1, double period2,
    double TimeSpan);
```

---

## Building the library

### Recommended: automake

Edit `USERLIBS/src/Makefile.am` and add three lines for your library:

```makefile
userlibs_LTLIBRARIES += mylib.la

mylib_la_SOURCES = mylib.c mylib.h
mylib_la_LIBADD  = $(abs_top_srcdir)/src/libvartools.la
mylib_la_LDFLAGS = -module
```

Then regenerate and build from the VARTOOLS source root:

```bash
./autogen.sh
./configure
make
make install   # installs mylib.so so -mylib works without -L
```

For commands that require external libraries (e.g., GSL, FFTW, Fortran), follow
the pattern of `fastchi2` (GSL) or `jktebop` (Fortran) in `Makefile.am`.

### Manual build (Linux)

If you prefer not to use autotools:

```bash
gcc -fPIC -c mylib.c -I/path/to/vartools/src
gcc -shared -Wl,-soname,mylib.so.1 \
    -o mylib.so.1.0 mylib.o \
    -L/path/to/vartools -lvartools
ln -sf mylib.so.1.0 mylib.so
```

### Manual build (macOS)

```bash
gcc -fPIC -c mylib.c -I/path/to/vartools/src
gcc -dynamiclib -o mylib.dylib mylib.o \
    -L/path/to/vartools -lvartools
```

---

## Complete minimal example

The `magadd` command (add a scalar offset to every magnitude) is the simplest
working example.  The full source is in `USERLIBS/src/magadd.c` and
`USERLIBS/src/magadd.h`.

**magadd.h:**

```c
typedef struct {
    double *addval;   /* one value per light curve */
} _Magadd;
```

**magadd.c** (key excerpts):

```c
void magadd_Initialize(char *commandname, int *RequireReadAll,
                       int *RequireSortLC, int *RequireDistinctTimes,
                       size_t *sizeuserdata)
{
    sprintf(commandname, "-magadd");
    *RequireReadAll = *RequireSortLC = *RequireDistinctTimes = 0;
    *sizeuserdata = sizeof(_Magadd);
}

int magadd_ParseCL(ProgramData *p, Command *c, void *userdata,
                   int *iret, char **argv, int argc)
{
    _Magadd *Magadd = (_Magadd *) userdata;
    if (VARTOOLS_ParseFixSpecFixcolumn(p, c, iret, argv, argc, 1,
            VARTOOLS_TYPE_DOUBLE, (void *)(&Magadd->addval),
            0, 1, "addval"))
        return 1;
    return 0;
}

void magadd_RunCommand(ProgramData *p, void *userdata,
                       int lc_name_num, int lc_num)
{
    _Magadd *Magadd = (_Magadd *) userdata;
    int N       = p->NJD[lc_num];
    double *mag = p->mag[lc_num];
    double val  = Magadd->addval[lc_num];
    for (int i = 0; i < N; i++)
        mag[i] += val;
}
```

**Usage:**

```bash
# Fix: add 0.5 to every LC
vartools -l lc_list -L USERLIBS/src/magadd.so -magadd fix 0.5 -tab

# List: read per-LC value from column 2 of the input list
vartools -l lc_list -L USERLIBS/src/magadd.so -magadd list column 2 -tab

# expr: add the median magnitude
vartools -l lc_list -L USERLIBS/src/magadd.so -magadd expr "Mag_median" -tab
```

---

## Using from Python (pyvartools)

```python
import pyvartools as vt

# Quick one-off
pipe = vt.Pipeline([
    vt.UserCommand("USERLIBS/src/mylib.so", "mylib", "fix 1.0")
])
result = pipe.run(lc)

# Or build a proper class
class MyLib(vt.UserCommand):
    _lib_path = "USERLIBS/src/mylib.so"
    _cmd_name = "mylib"
    def __init__(self, value):
        super().__init__(self._lib_path, self._cmd_name, f"fix {value}")

pipe = vt.Pipeline([MyLib(1.0)])
```

See [User Extension Commands](../python/commands.md#user-extension-commands)
for the full pyvartools API.
