# Control Flow and Output

Commands and options for conditional execution, saving and restoring light curve state, writing output light curves, printing variables, modifying column names, and calling external Python or R code.

---

## `-if` / `-elif` / `-else` / `-fi`

```
-if <expression>

    [-command1 ... -commandN]

[-elif <expression>

    [-command1 ... -commandN]
]

[-else

    [-command ... -command]
]

-fi
```

Make execution of VARTOOLS commands conditional on the evaluation of an analytic expression. If `expression` evaluates to `0` (when cast to an integer), the commands in that branch are skipped; any non-zero integer value causes the branch to execute.

**Behavior**

- `-elif <expression>` and `-else` provide alternative branches, evaluated in order.
- The construct is terminated by `-fi`.
- Nested `-if`/`-elif`/`-else`/`-fi` constructs are allowed.
- Any conditional construct not explicitly terminated with `-fi` is assumed to be closed after the last command on the command line.

**Expressions**

The `expression` can reference any variable computed by a preceding command (use `-headeronly` to see variable names), vectors in the light curve, or scalar constants. Use `-functionlist` for the full list of supported operators, functions, and constants.

**Example**

```bash
vartools -l lclist.txt \
    -LS 0.5 20.0 4.0 1 0 \
    -if "LS_lnFAP_1_1 > 10.0" \
        -Killharm ls 2 0 0 /tmp/models \
    -fi
```

!!! caution
    Conditional constructs are **ignored** by commands that process all light curves simultaneously (e.g., `-SYSREM`, `-findblends`) and by the `-savelc` and `-restorelc` commands.

**Examples**

**Example 1.** Using nested `-if`, `-elif`, `-else`, `-fi` constructs to apply different statistics depending on the values returned by `-rms`.

```bash
vartools -l EXAMPLES/lc_list -rms \
    -if 'RMS_0>10*Expected_RMS_0' \
        -if 'RMS_0 > 0.1' \
            -stats mag stddev \
        -else \
            -stats mag pct30 \
        -fi \
    -elif 'Npoints_0>3900' \
        -stats mag kurtosis \
    -else \
        -rms \
    -fi \
    -header
```

Output:
```
table with columns Name, Mean_Mag_0, RMS_0, Expected_RMS_0, Npoints_0, STATS_mag_STDDEV_3, STATS_mag_PCT30.00_5, STATS_mag_KURTOSIS_8, Mean_Mag_10, RMS_10, Expected_RMS_10, Npoints_10 for 10 entries.
```

---

## `-savelc` / `-restorelc`

See the [Filtering and Detrending](filtering.md#-savelc) page for full documentation of these commands.

**Quick reference**

```
-savelc
```

Save the current light curve state.

```
-restorelc
    savenumber ["vars" var1,var2,...]
```

Restore the light curve to the state saved at the `savenumber`-th call to `-savelc`.

**Examples**

**Example 1.** Running a battery of variability selection algorithms on a list of light curves in parallel, using `-savelc` and `-restorelc` to branch from different saved states.

```bash
vartools -l EXAMPLES/lc_list -header -numbercolumns \
    -nobuffer -parallel 4 \
    -savelc \
    -clip 5.0 1 \
    -savelc \
    -LS 0.1 100. 0.1 3 0 clip 5. 1 \
    -aov 0.1 100. 0.1 0.01 1 0 clip 5. 1 \
    -restorelc 1 \
    -clip 10.0 1 \
    -BLS q 0.01 0.1 0.1 20. 10000 200 7 2 0 0 0 \
    -restorelc 2 \
    -changeerror \
    -autocorrelation 0. 30. 0.1 EXAMPLES/OUTDIR1/ \
  > EXAMPLES/OUTDIR1/vartools.out
```

---

## `-o`

```
-o
    <outdir | outname>
    ["nameformat" formatstring | "namecommand" command |
     "namefromlist" ["column" col]]
    ["columnformat" formatstring]
    ["delimiter" delimchar]
    ["fits"] ["copyheader"] ["logcommandline"] ["noclobber"]
```

Output the current light curve(s) to files or to stdout. When processing a list of light curves, provide an output directory (`outdir`); when processing a single light curve, provide a filename (`outname`). Use `"-"` for `outname` to write to stdout (combine with `-quiet` or `-redirectstats` to avoid mixing the statistics table with the light curve data).

**Naming conventions**

By default, output files are named `$outdir/$inname` where `$inname` is the base filename of the input light curve.

- `"nameformat" formatstring` — Override the naming convention using a format string:
  - `%s` → input basename
  - `%b` → input basename stripped of its last extension
  - `%d` → light curve number (starting from 1)
  - `%0nd` → zero-padded light curve number (`n` digits)
  - `%%` → literal `%`
- `"namecommand" command` — Execute a shell command to determine the output name. The string `echo $fulllcname $outdir $lcnum` is piped into `command`; its stdout is the output filename.
- `"namefromlist" ["column" col]` — Read the output filename from the input list file. `outdir` is prepended. By default the next unused column is used; use `"column"` to specify it.

**Column format**

By default the output has three columns: `t`, `mag`, `err`.

- `"columnformat" formatstring` — Comma-separated list of `varname[:printf_format]` entries. For example: `"columnformat t:%.17g,mag:%.5f,err:%.5f,xpos:%.3f"`. For FITS output, the text after the first `:` specifies the column units, a second `:` introduces a comment description, and an optional third `:` provides an alternative units string for the ASTROPY-SERIALIZED-COLUMNS section.

**Delimiter**

- `"delimiter" delimchar` — Character to use as column separator (default: single space).

**Output format**

- `"fits"` — Output in binary FITS table format. The `.fits` extension is appended if not already present.
- `"copyheader"` — Copy the primary FITS header from the input light curve (if it was FITS format) to the output.
- `"logcommandline"` — Log the VARTOOLS command line to the output file header.
- `"noclobber"` — Do not overwrite existing files; VARTOOLS terminates if a file already exists.

**Example**

```bash
# Write all light curves to /out with a custom name pattern
vartools -l lclist.txt -rms \
    -o /out nameformat %b.detrended.lc
```

---

## `-print`

```
-print
    var1[,var2,var3...]
    ["columnnames" col1[,col2,col3...]]
    ["format" fmt1[,fmt2,fmt3...]]
```

Print the value of one or more variables to the output statistics table. This is the primary way to include user-computed scalars or expressions in the output.

If a variable is a light-curve vector (length equal to the number of points in the light curve), only the value at the **first** point is output.

**Parameters**

- `var1[,var2,...]` — Comma-separated list of previously defined variable names to print.
- `"columnnames" col1[,col2,...]` — Override the default output column names. By default columns are named `Print_${varname}_${varnum}_${commandnum}`. The command-number suffix `_${commandnum}` is still appended unless `-columnsuffix` is used.
- `"format" fmt1[,fmt2,...]` — Printf-style format codes for each value (e.g., `%.5f`).

**Example**

```bash
vartools -l lclist.txt \
    -LS 0.5 20.0 4.0 1 0 \
    -print LS_Period_1_1,LS_lnFAP_1_1 \
           columnnames Period,lnFAP \
           format %.6f,%.3f
```

---

## `-addfitskeyword`

```
-addfitskeyword
    keyword
    ["combinelc" lcnumvar]
    < "TDOUBLE" | "TINT" | "TLONG" | "TSTRING" >
    < "fix" val | "var" variable >
    ["comment" commentstring]
    ["primary" | "extension"]
    ["append" | "update"]
```

Add a keyword to the header of any subsequently output FITS-format light curves.

**Parameters**

- `keyword` — The FITS header keyword to add. When `"combinelc"` is used, this should be a format string containing `%d` which will be replaced by the unique values of `lcnumvar`.
- `"combinelc" lcnumvar` — Add multiple keywords, one per unique value of `lcnumvar`. The value of the associated variable is taken at the first occurrence of each new unique value of `lcnumvar`.
- `"TDOUBLE" | "TINT" | "TLONG" | "TSTRING"` — Data type for the keyword value.
- `"fix" val` — Set the keyword value to a fixed constant for all light curves.
- `"var" variable` — Set the keyword value from a VARTOOLS variable.
- `"comment" commentstring` — Optional comment string for the FITS header entry.
- `"primary" | "extension"` — Place the keyword in the primary FITS header (`"primary"`, default) or in the binary table extension header (`"extension"`).
- `"append" | "update"` — If the keyword already exists: `"append"` adds a duplicate entry; `"update"` (default) overwrites the existing value.

---

## `-columnsuffix`

```
-columnsuffix
    suffix
```

Change the suffix appended to output column names for the **next** command called on the command line. By default VARTOOLS appends the command number (e.g., `_1`, `_2`) to every column name to keep them unique when the same command is called multiple times.

**Parameters**

- `suffix` — Replacement suffix string. Supply an empty string `""` to remove the suffix entirely.

**Notes**

- This option applies only to the immediately following command. Call it again before each subsequent command if needed.
- If the same command is called more than once with the same suffix, VARTOOLS will error because column names must be unique.

**Example**

```bash
vartools -l lclist.txt \
    -columnsuffix "" -rms \
    -columnsuffix "_detrend" -TFA trendlist ... 1 0 0 0 \
    -columnsuffix "" -rms
```

---

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
vartools -l lclist.txt \
    -python 'iqr = float(numpy.percentile(mag, 75) - numpy.percentile(mag, 25))' \
    invars mag \
    outputcolumns iqr
```
