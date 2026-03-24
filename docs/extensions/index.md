# Extending VARTOOLS

VARTOOLS provides two extension mechanisms that let you add your own code
without modifying the core source:

| Extension type | CLI flag | What it adds |
|----------------|----------|--------------|
| [User Commands](user-commands.md) | `-L libname.so` | New pipeline commands (e.g., `-mycommand …`) |
| [User Analytic Functions](user-functions.md) | `-F libname.so` | New functions for the expression evaluator (e.g., `myfunc(t, mag)`) |

Both are compiled as shared-object libraries that VARTOOLS loads at runtime.
If the libraries are installed into the VARTOOLS data directory (`make install`
handles this automatically), the full path to the `.so` file is not required.

---

## User Commands (`-L`)

A **user command library** adds a new pipeline stage.  It can read input
parameters from the command line or the input list file, access the light curve
arrays, compute results, store computed values, and output them in the standard
VARTOOLS ASCII table.

Suitable for:

- New period-search or model-fitting algorithms
- Custom detrending or calibration steps
- Anything that reads or modifies the light curve and/or produces scalar output

→ [Writing User Extension Commands](user-commands.md)

---

## User Analytic Functions (`-F`)

A **user function library** adds new functions to the VARTOOLS analytic
expression evaluator.  Functions take a fixed number of scalar arguments and
return a scalar.  They can be used anywhere an expression is accepted: `-expr`,
`"expr"` keywords in command parameters, `-Injectharm` signal descriptions,
etc.

Suitable for:

- Custom window functions, apodisations, or signal shapes
- Domain-specific transformations (e.g., limb-darkening laws, extinction curves)
- Any pure mathematical function not already built into VARTOOLS

→ [Writing User Analytic Functions](user-functions.md)

---

## Bundled examples

| Library | Type | Location | What it does |
|---------|------|----------|--------------|
| `magadd` | Command | `USERLIBS/src/magadd.c` | Add a constant to magnitudes — simplest possible command example |
| `examplefunction` | Function | `USERFUNCS/src/examplefunction.c` | Adds `addvals(a,b)` — simplest possible function example |
| `fastchi2` | Command | `USERLIBS/src/fastchi2.c` | Palmer's Fast χ² periodogram (requires GSL) |
| `jktebop` | Command | `USERLIBS/src/jktebop.c` | Detached eclipsing binary model (Fortran core) |
| `macula` | Command | `USERLIBS/src/macula_c.c` | Kipping's analytic starspot model (Fortran 90 core) |
| `stitch` | Command | `USERLIBS/src/stitch.c` | Zero-point offset removal for combined multi-telescope LCs |
| `splinedetrend` | Command | `USERLIBS/src/splinedetrend.c` | Multivariate spline/polynomial detrending (requires GSL) |
| `astrofuncs` | Function | `USERFUNCS/src/astrofuncs.c` | Astronomical utility functions |
