# Writing User Analytic Functions

VARTOOLS has a built-in analytic expression evaluator used throughout the
program — by `-expr`, the `"expr"` keyword in command parameters, analytic
signal descriptions in `-Injectharm`, and many other places.  You can extend
this evaluator with your own functions by compiling them into a shared-object
library and loading it with the `-F` flag:

```bash
vartools -F /path/to/myfuncs.so -expr "myfunc(t, mag)" ...
```

If the library is installed into the VARTOOLS `userfuncs` data directory
(`make install` does this automatically), the full path is not needed:

```bash
vartools -F myfuncs -expr "myfunc(t, mag)" ...
```

User function libraries are simpler than user command libraries.  They do not
interact with the pipeline or the output table — they just extend the vocabulary
of the expression evaluator.  A single library can register any number of
functions.

The simplest working example is `examplefunction` in
`USERFUNCS/src/examplefunction.c`.

---

## Prerequisites

- C (or any language callable from C).
- The vartools header: `#include "../../src/vartools.h"` (adjust path as
  needed).
- Link against `libvartools.a` when building.
- GNU autotools for the recommended build path (same as for user commands —
  see [Writing User Extension Commands](user-commands.md)).

---

## Structure of a function library

A function library requires exactly **one** specially-named function,
`$LIBNAME_Initialize`, plus one C function for each analytic function the
library provides.

### The analytic functions

Each function must have the signature:

```c
double myfunc(double *param);
```

- `param[0]` is the first argument, `param[1]` the second, and so on.
- The function returns a scalar `double`.

Arguments passed to the function in an expression are evaluated to `double`
before the call; inside the function they are plain numbers — no special
VARTOOLS types or structures are involved.

### `$LIBNAME_Initialize`

```c
void myfuncs_Initialize(ProgramData *p)
{
    /* declare each function first so the compiler knows its type */
    double myfunc(double *);
    double anotherfunc(double *);

    /* register each function with VARTOOLS */
    VARTOOLS_RegisterUserFunction(p, "myfunc",     2, &myfunc,     1,
        "brief description of what myfunc returns",
        "x", "first argument description",
        "y", "second argument description");

    VARTOOLS_RegisterUserFunction(p, "anotherfunc", 1, &anotherfunc, 0);
    /* last 0 = no help text provided */
}
```

`VARTOOLS_RegisterUserFunction` signature:

```c
void VARTOOLS_RegisterUserFunction(
    ProgramData *p,
    const char  *name,       /* name to use in expressions */
    int          nargs,      /* number of arguments the function takes */
    double      (*func)(double *),  /* pointer to the C function */
    int          provide_help,      /* 1 if help text follows, 0 if not */
    /* if provide_help == 1, supply: */
    const char  *summary,           /* one-line description of return value */
    /* then, for each argument: */
    const char  *arg_name,          /* symbolic name */
    const char  *arg_description,   /* brief description */
    ...
);
```

When `provide_help == 1` the function will appear with descriptions when the
user runs `vartools -F myfuncs -functionlist`.

---

## Complete minimal example

```c
/* myfuncs.c */
#include "../../src/vartools.h"   /* adjust path to vartools.h */
#include <math.h>

/* ------------------------------------------------------------------ */
/* Analytic functions                                                  */
/* ------------------------------------------------------------------ */

/* gaussian(x, mu, sigma) — unnormalized Gaussian evaluated at x */
double gaussian(double *p) {
    double x     = p[0];
    double mu    = p[1];
    double sigma = p[2];
    return exp(-0.5 * ((x - mu) / sigma) * ((x - mu) / sigma));
}

/* sinc(x) — sin(x)/x, with sinc(0) = 1 */
double sinc(double *p) {
    double x = p[0];
    if (x == 0.0) return 1.0;
    return sin(x) / x;
}

/* ------------------------------------------------------------------ */
/* Initialization                                                      */
/* ------------------------------------------------------------------ */

void myfuncs_Initialize(ProgramData *p)
{
    double gaussian(double *);
    double sinc(double *);

    VARTOOLS_RegisterUserFunction(p, "gaussian", 3, &gaussian, 1,
        "unnormalized Gaussian exp(-0.5*((x-mu)/sigma)^2)",
        "x",     "evaluation point",
        "mu",    "mean",
        "sigma", "standard deviation");

    VARTOOLS_RegisterUserFunction(p, "sinc", 1, &sinc, 1,
        "sinc(x) = sin(x)/x, with sinc(0) = 1",
        "x", "argument in radians");
}
```

**Build (see [Building](#building) below), then use:**

```bash
# Evaluate at fixed arguments
vartools -l lc_list -expr "mygaussian = gaussian(t, 2459000.0, 5.0)" -tab

# Load at runtime with -F
vartools -F ./myfuncs.so \
         -l lc_list \
         -expr "w = gaussian(mag, Mag_median, RMS)" \
         -tab
```

---

## Building

### Recommended: automake

Edit `USERFUNCS/src/Makefile.am`:

```makefile
userfuncs_LTLIBRARIES += myfuncs.la

myfuncs_la_SOURCES = myfuncs.c
myfuncs_la_LIBADD  = $(abs_top_srcdir)/src/libvartools.la
myfuncs_la_LDFLAGS = -module
```

Then from the VARTOOLS source root:

```bash
./autogen.sh
./configure
make
make install   # installs myfuncs.so so -F myfuncs works without a path
```

### Manual build (Linux)

```bash
gcc -fPIC -c myfuncs.c -I/path/to/vartools/src
gcc -shared -Wl,-soname,myfuncs.so.1 \
    -o myfuncs.so.1.0 myfuncs.o \
    -L/path/to/vartools -lvartools -lm
ln -sf myfuncs.so.1.0 myfuncs.so
```

### Manual build (macOS)

```bash
gcc -fPIC -c myfuncs.c -I/path/to/vartools/src
gcc -dynamiclib -o myfuncs.dylib myfuncs.o \
    -L/path/to/vartools -lvartools -lm
```

---

## Listing available functions

After loading a library, `vartools -F myfuncs -functionlist` prints every
registered function and its help text (when `provide_help == 1`).

---

## Tips

- **Multiple libraries in one call.** Multiple `-F` flags may be given; each
  file is loaded in order and all functions from all files become available.
- **Naming.** Function names share the same namespace as the built-in
  expression evaluator.  Avoid names that clash with common C math library
  names (`sin`, `cos`, `log`, etc.), which are already built in.
- **No state.** User functions receive only their arguments and must be
  stateless (pure functions).  They are evaluated repeatedly and may be
  called from parallel threads — never use static or global variables.
- **Integer arguments.** The evaluator passes all arguments as `double`.  If
  you need an integer-typed argument (e.g., a number of harmonics), just
  use `(int)param[0]` inside the function.
- **Calling libvartools utilities.** Math and statistics helpers from
  `libvartools` (splines, statistics, sorting — documented in
  [Writing User Extension Commands](user-commands.md#libvartools-utility-functions))
  are available in user function libraries too, since you link against
  `libvartools.la`.
