# Control Flow and Output

Commands and options for conditional execution, saving and restoring light curve state, writing output light curves, printing variables, modifying column names, and calling external Python or R code.

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

## `-changevariable`

**Syntax**
```
-changevariable
    < "t" | "mag" | "err" | "id" > var
```

**Description**

Reassign the internal role of the time (`t`), magnitude (`mag`), magnitude uncertainty (`err`), or image-identifier (`id`) variable to a different named column. Subsequent commands will use the new assignment. The original variable still exists; to restore it, issue `-changevariable mag mag`.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `"t"` / `"mag"` / `"err"` / `"id"` | Which built-in role to reassign. |
| `var` | Name of the existing light-curve column to promote to that role. |

**Examples**

**Example 1.** Run an LS period search, then use `-changevariable` to store the current time in a `phase` variable, phase-fold on the LS period, then swap `t` back before writing output so the light curve is sorted by time.

```bash
vartools -l EXAMPLES/lc_list \
    -LS 0.1 100.0 0.1 1 0 \
    -expr 'phase=t' \
    -changevariable t phase \
    -Phase ls \
    -changevariable t t \
    -o EXAMPLES/OUTDIR1 nameformat "%s.phase.txt" \
        columnformat "t:%17.9f,mag:%9.5f,err:%9.5f,phase:%9.5f" \
    -header
```

Output: a table with columns `Name`, `LS_Period_1_0`, `Log10_LS_Prob_1_0`, `LS_SNR_1_0` for each light curve in the list, plus one phase-folded output file per light curve.

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

The `expression` can reference any variable computed by a preceding command (use `-headeronly` to see variable names), vectors in the light curve, or scalar constants. See the [Analytic Expressions](expressions.md) reference for the full list of supported operators, functions, and constants.

**Example**

```bash
vartools -l EXAMPLES/lc_list \
    -LS 0.5 20.0 4.0 1 0 \
    -if 'Log10_LS_Prob_1_0 < -10.0' \
        -Killharm ls 2 0 1 EXAMPLES/OUTDIR1 \
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

## `-o`

```
-o
    <outdir | outname>
    ["nameformat" formatstring | "namecommand" command |
     "namefromlist" ["column" col]]
    ["columnformat" formatstring | "allcols"]
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
- `"allcols"` — Write every light-curve-vector variable defined by commands *before* this `-o` in the sequence, in the order they were defined, using a type-appropriate default `printf` format. Mutually exclusive with `"columnformat"`. ASCII output includes a `# name1 name2 …` header line so downstream readers can recover the column names. Useful when you want to capture every vector a prior command created (e.g. `-Phase phasevar`, `-linfit modelvar`) without listing each one by hand.

**Delimiter**

- `"delimiter" delimchar` — Character to use as column separator (default: single space).

**Output format**

- `"fits"` — Output in binary FITS table format. The `.fits` extension is appended if not already present.
- `"copyheader"` — Copy the primary FITS header from the input light curve (if it was FITS format) to the output.
- `"logcommandline"` — Log the VARTOOLS command line to the output file header.
- `"noclobber"` — Do not overwrite existing files; VARTOOLS terminates if a file already exists.

**Example**

```bash
# Write all light curves to a directory with a custom name pattern
vartools -l EXAMPLES/lc_list -rms \
    -o EXAMPLES/OUTDIR1 nameformat "%s.detrended.lc"
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
vartools -l EXAMPLES/lc_list \
    -LS 0.5 20.0 4.0 1 0 \
    -print LS_Period_1_0,Log10_LS_Prob_1_0 \
           columnnames Period,Log10_FAP \
           format %.6f,%.3f
```

---

## `-savelc` / `-restorelc`

Checkpoint and restore the in-memory light-curve state. Useful for running several analysis branches from the same starting point without re-reading files, and for restoring the original light curve after a destructive transformation.

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

