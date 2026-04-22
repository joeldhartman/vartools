# Calling Python or R from VARTOOLS

VARTOOLS can invoke arbitrary user-supplied Python or R code on each
light curve (or on the full collection) via the `-python` and `-R`
commands. Both embed an interpreter at compile time, so the features are
only available if the corresponding development libraries were installed
when VARTOOLS was built (see the [Installation](../install.md) page).

## `-R`

```
-R
    < "fromfile" commandfile | commandstring >
    ["init" < "file" initializationfile | initializationstring >
        | "continueprocess" prior_R_command_number]
    ["vars" variablelist
        | ["invars" inputvariablelist] ["outvars" outputvariablelist]]
    ["outputcolumns" variablelist] ["process_all_lcs"] ["verbose"]
```

Execute arbitrary R code on each light curve. VARTOOLS embeds the user-supplied code in an R function, compiles it, and calls it once per light curve (or once for all light curves if `"process_all_lcs"` is given).

Variables are passed to R as **vectors**. VARTOOLS automatically imports `numpy` is **not** applicable here; R vectors are used natively.

!!! note "Environment requirement"
    The `R_HOME` environment variable must be set before calling VARTOOLS. Find the correct value by running `R RHOME` from the shell. Adding `export R_HOME=$(R RHOME)` to your `.bashrc` is recommended.

**Parameters**

- `"fromfile" commandfile` — Read R code from a file rather than the command line.
- `commandstring` — R code provided directly as a string on the command line.
- `"init"` — R code to execute **once** before processing begins (e.g., library imports or function definitions):
  - `"file" initializationfile` — Load initialization code from a file.
  - `initializationstring` — Provide initialization code as a string.
- `"continueprocess" prior_R_command_number` — Reuse the sub-process from the `prior_R_command_number`-th `-R` command (1-indexed). This allows subsequent `-R` calls to share state (global variables, loaded libraries) with the first call. When this is used, no initialization code should be provided.
- `"vars" variablelist` — Comma-separated list of variables to both pass **into** and **receive back from** the R function.
- `"invars" inputvariablelist` — Variables passed into the R function only.
- `"outvars" outputvariablelist` — Variables received back from the R function only.
- `"outputcolumns" variablelist` — Variables to include in the output ASCII statistics table. Light curve vectors and string variables cannot be listed here.
- `"process_all_lcs"` — Pass all light curves to R simultaneously. Vectors are supplied as lists of vectors; scalar variables are supplied as lists.
- `"verbose"` — Allow R to print messages to stdout (default: R runs in `--slave` mode). Note that verbose output may interfere with the VARTOOLS statistics table.

**Parallelism**

When VARTOOLS is run with `-parallel`, a separate R sub-process is launched per thread. Initialization code is executed separately for each thread; global variable state is not shared between threads.

---

## `-python`

```
-python
    < "fromfile" commandfile | commandstring >
    ["init" < "file" initializationfile | initializationstring >
        | "continueprocess" prior_python_command_number]
    ["vars" variablelist
        | ["invars" inputvariablelist] ["outvars" outputvariablelist]]
    ["outputcolumns" variablelist] ["process_all_lcs"] ["skipfail"]
```

Execute arbitrary Python code on each light curve. VARTOOLS embeds the user-supplied code in a Python function, compiles it via the Python C API, and calls it once per light curve. `numpy` is automatically imported and available without an explicit `import numpy` statement.

**Parameters**

- `"fromfile" commandfile` — Read Python code from a file rather than the command line.
- `commandstring` — Python code provided directly as a string on the command line.
- `"init"` — Python code to execute **once** before processing begins (e.g., function definitions or library imports):
  - `"file" initializationfile` — Load initialization code from a file.
  - `initializationstring` — Provide initialization code as a string.
- `"continueprocess" prior_python_command_number` — Reuse the sub-process from the `prior_python_command_number`-th `-python` command (1-indexed). Allows subsequent `-python` calls to share state with the first call. No initialization code should be provided when using this keyword.
- `"vars" variablelist` — Comma-separated list of variables to both pass **into** and **receive back from** the Python function.
- `"invars" inputvariablelist` — Variables passed into the Python function only.
- `"outvars" outputvariablelist` — Variables received back from the Python function only.
- `"outputcolumns" variablelist` — Variables to include in the output ASCII statistics table. Light curve vectors and string variables cannot be listed here.
- `"process_all_lcs"` — Pass all light curves to Python at once. Numeric vectors are supplied as lists of numpy arrays; scalar variables are supplied as numpy arrays.
- `"skipfail"` — If the Python code raises an error for a given light curve, skip subsequent processing of that light curve and continue to the next, rather than terminating VARTOOLS.

**Variable types in Python**

- Numeric vectors: `numpy` arrays.
- String data: Python lists.

**Parallelism**

When VARTOOLS is run with `-parallel`, a separate Python sub-process is spawned per thread, avoiding the Python Global Interpreter Lock. Initialization code is executed independently for each thread; global variable changes in one thread are not visible to others.

**Example**

```bash
# Compute the interquartile range using Python and add it to the output table
vartools -l EXAMPLES/lc_list \
    -python 'iqr = float(numpy.percentile(mag, 75) - numpy.percentile(mag, 25))' \
    invars mag \
    outvars iqr \
    outputcolumns iqr
```
