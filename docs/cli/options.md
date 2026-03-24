# Input, Output and Processing Options

These flags are not analysis commands — they control how VARTOOLS reads light
curves, formats its output, and manages processing. They may be combined freely
with any sequence of commands.

---

## Input options

### `-i lcname`

Read a single light curve from file `lcname`. If `lcname` is `-` then input is
taken from stdin. Add the `"binary"` keyword after the filename to read a
light curve in Penev binary format. FITS binary tables are detected
automatically from the `.fits` file extension.

```bash
vartools -i EXAMPLES/2 -rms
```

### `-l lclist`

Read a list of light curves from the file `lclist`, one filename per line (by
default the first column). Use `-` to read the list from stdin.

Key sub-options:

| Sub-option | Meaning |
|---|---|
| `"column" col` | Read filenames from column `col` instead of column 1. Use `"all"` to treat every whitespace-separated token on each line as a filename. |
| `"combinelcs"` | Combine multiple files named on a single list line (comma-separated by default) into one light curve in memory. |
| `"delimiter" delim` | Change the delimiter used to separate filenames in a `"combinelcs"` list (default: comma). |
| `"lcnumvar" varname` | Store the origin-file index for each observation (0, 1, 2, …) in variable `varname` when `"combinelcs"` is active. |
| `"binary"` | Input files are in Penev binary format. |
| `"opencommand" command` | Pass each filename through `command` (run in a shell with `%s` replaced by the filename) and read the light curve from that command's stdout. |

```bash
vartools -l EXAMPLES/lc_list -rms -tab
```

### `-inputlcformat var1:col1[:type1[:fmt1]][,var2:col2[:type2[:fmt2]],...] ["skipnum" Nskip] ["skipchar" chars] ["delimiter" delim]`

Specify the column layout of the input light curve files. Each token
`varN:colN` maps a symbolic variable name to a column number (1-indexed).

The reserved names `t`, `mag`, `err`, and `id` refer to the time, magnitude,
magnitude uncertainty, and string image identifier. Any other name creates an
auxiliary variable that can be referenced by later commands.

Optional type codes: `double` (default), `float`, `int`, `long`, `short`,
`string`, `char`, `utc`. The `utc` type reads a UTC timestamp and converts it
to JD; you must supply a format string such as `"%Y-%M-%DT%h:%m:%s"`.

Setting `col` to `0` creates the variable without reading it from the file; the
format field is then an analytic initializer expression (special variable `NR`
gives the zero-based record number).

`"skipnum" Nskip` skips `Nskip` header lines. `"skipchar" chars` changes the
comment character(s) (default `#`). `"delimiter" delim` changes the column
separator (default: whitespace).

```bash
vartools -l EXAMPLES/lc_list \
    -inputlcformat t:1,mag:2,err:3,xpos:4,ypos:5 \
    -rms
```

For FITS or Penev binary files, column identifiers may be FITS column header
names rather than integers. Column positions are resolved once from the first
light curve and are not re-checked for subsequent files.

### `-inlistvars var1:col1[:"combinelc"][:vtype1[:fmt1]][,...]`

Read one or more scalar variables from the input list file rather than from the
light curves themselves. Each token `varN:colN` maps a name to a column in the
list. Variable names follow the same rules as `-inputlcformat`; the reserved
names `t`, `mag`, `err`, and `id` cannot be used.

Setting `col` to `0` creates the variable without reading it; in that case the
format field is an analytic expression that may reference previously defined
variables and the special variable `NF` (zero-based record number in the list).

When the `-l "combinelcs"` option is active, add the `"combinelc"` qualifier
after the column number to supply a separate value for each combined-file index.
The value in the list is then a comma-separated sequence matching the number of
input files for that line.

```bash
vartools -l EXAMPLES/lc_list \
    -inlistvars period:2:double \
    -LS fix 0.5 10.0 0.01 1 0
```

### `-readformat Nskip ["stringid" colstringid] ["inpututc" format] col_time col_mag col_sig`

!!! warning "Deprecated"
    `-readformat` was deprecated in VARTOOLS 1.3. Use `-inputlcformat` for all
    new scripts.

`Nskip` is the number of header lines to skip at the start of each file
(in addition to lines beginning with `#`). `col_time`, `col_mag`, and `col_sig`
are the columns for the time, magnitude, and magnitude uncertainty (defaults:
1, 2, 3). Setting a column to `0` suppresses reading it: times become
line numbers, magnitudes become 0, uncertainties become 1.

`"stringid" colstringid` reads a string image identifier from column
`colstringid` (required if `-matchstringid` is set).

`"inpututc" format` treats the time column as a UTC string and converts it to
JD. The format syntax is the same as for `-inputlcformat`.

---

## Output options

### `-oneline`

Print each output statistic on its own line in `Name = Value` format instead
of the default multi-column table. Most useful when processing a single light
curve for interactive inspection.

```
LS_Period_1_0                        =     1.23588006
Log10_LS_Prob_1_0                    = -707.23302
```

### `-tab`

Use a tab-delimited (starbase) format for the output table instead of the
default space-delimited layout. Column names appear in the header row.

### `-header` / `-headeronly`

`-header` prepends a single header line of column names to the output table.

`-headeronly` prints the header and then exits without processing any light
curves. Use this to inspect the output column schema for a given command set
before running a large batch job.

### `-numbercolumns`

Prefix each column name in the header with its sequential column number:

```
1:Name   2:Mean_Mag_0   3:RMS_0   ...
```

### `-redirectstats statsfile ["append"]`

Write the statistics table to `statsfile` rather than stdout. Add `"append"` to
append to an existing file. Useful when you also want to stream processed light
curves to stdout as part of a pipeline (e.g., combining with `-o -`).

### `-nobuffer`

Disable stdout buffering so that output lines appear immediately. By default
stdout is buffered.

### `-quiet`

Suppress the statistics table entirely. Light curve output files and
periodogram files requested with `-o`, `-oper`, etc. are still written.

### `-log-command-line`

Print the full VARTOOLS command line as a comment at the top of the output,
before the header. This embeds provenance information directly in the output
file.

### `-o <outdir | outname> ["nameformat" formatstring] ["columnformat" formatstring]`

Write the (processed) light curve to a file. If `outdir` is a directory,
output files are placed inside it using the base name of the input file. Use
`"nameformat"` to control the output filename pattern and `"columnformat"` to
select and format the columns written.

**Examples**

**Example 1.** Write phase-folded light curves with custom filename and column formats.

```bash
vartools -l EXAMPLES/lc_list -header \
    -LS 0.1 100.0 0.1 1 0 \
    -expr phase=t \
    -changevariable t phase \
    -Phase ls \
    -o EXAMPLES/OUTDIR1 \
        nameformat "file_%s_%05d_simout.txt" \
        columnformat "t:%11.5f,phase:%8.5f,mag:%7.4f,err:%7.4f"
```

Output table:
```
#Name LS_Period_1_0 Log10_LS_Prob_1_0 LS_SNR_1_0
EXAMPLES/1     0.97821072 -452.25157   41.33409
EXAMPLES/2     1.23440877 -704.49194   58.45119
EXAMPLES/3     1.14786351  -30.00548   15.74701
EXAMPLES/4    14.81290524  -59.52748   13.11947
EXAMPLES/5     7.40645262  -53.86771   10.01489
EXAMPLES/6     0.96306814  -42.42348   10.53479
EXAMPLES/7     0.32704113  -11.84669    4.77871
EXAMPLES/8     3.07991099  -88.30735   15.34709
EXAMPLES/9     7.23420953  -37.93155   14.15476
EXAMPLES/10     0.96906857  -40.55309   11.32727
```

Sample output (`head -3 EXAMPLES/OUTDIR1/file_1_00001_simout.txt`):
```
53725.17392  0.00000 10.0850  0.0012
53726.15280  0.00068 10.0886  0.0009
53726.15378  0.00169 10.0918  0.0009
```

This example demonstrates the `nameformat` and `columnformat` keywords for the `-o` command. Light curves are read in from the list and `-LS` locates periods. `-expr` defines a new vector `phase` initialized to the times in the light curves; `-changevariable` causes subsequent commands to use `phase` where time would normally be used. Together with `-Phase`, this stores the light curve phase for the LS period in `phase`. The light curves are then output to `EXAMPLES/OUTDIR1`. The `nameformat` keyword gives the rule for naming the output files — the first light curve (`EXAMPLES/1`) yields output to `EXAMPLES/OUTDIR1/file_1_00001_simout.txt`, the second to `EXAMPLES/OUTDIR1/file_2_00002_simout.txt`, and so on. The `columnformat` keyword specifies how data are formatted in the output light curve: four quantities (time, phase, magnitude, error) are included with printf-style format strings. Without `columnformat`, only `t`, `mag`, and `err` would be output, formatted as `%17.9f`, `%9.5f`, and `%9.5f` respectively.

---

## Help options

### `-help` / `-help all` / `-help commandname`

Print usage information. `-help` gives a brief summary; `-help all` lists all
commands and options; `-help commandname` prints the full documentation for a
specific command.

### `-example commandname`

Show an annotated example invocation for `commandname`.

### `-listcommands`

Print a terse list of all available commands and options, one per line.

### `-functionlist`

List all functions, operators, constants, and special variables understood by
the VARTOOLS analytic expression evaluator (used in `-expr`, `-linfit`, `-if`,
etc.).

### `-showinputlcformat`

Print the expected light curve format (columns and types) implied by the current
`-inputlcformat` specification and then exit.

### `-showinputlistformat`

Print the expected format of the input light curve list implied by the current
`-l` and `-inlistvars` settings and then exit. (Called `-inputlistformat` in
older versions.)

---

## Processing options

### `-parallel Nthreads`

Process up to `Nthreads` light curves simultaneously on a multi-core machine
using POSIX threads. Results are written as each thread finishes, so the output
order may differ from the input order.

```bash
vartools -l EXAMPLES/lc_list -parallel 8 -LS 0.5 10.0 0.01 1 0 -tab
```

### `-randseed seed`

Seed the random number generator used by commands such as `-addnoise` and
`-copylc`. `seed` must be a positive integer, or the string `"time"` to seed
from the system clock. If not specified, the seed defaults to 1 (fully
reproducible).

### `-jdtol jdtol`

Two time stamps are considered equal if they differ by less than `jdtol` days.
This tolerance is used by commands that match observations across light curves,
such as `-TFA`. The default value is `0.000010` days (about 0.86 seconds).

### `-matchstringid`

Match observations from different light curves by string image identifier rather
than by time. When this option is set you must declare the `id` column with
`-inputlcformat`.

### `-skipmissing`

Continue processing when a light curve file cannot be found or read, skipping
that entry, rather than aborting. Useful in large batch runs where a small
number of files may be missing.
