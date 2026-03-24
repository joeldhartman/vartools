# VARTOOLS CLI Overview

VARTOOLS is a command-line utility for processing and analyzing astronomical time-series data, especially light curves. It is designed primarily for batch processing large collections of light curves. A detailed scientific reference is provided in [Hartman and Bakos, 2016, Astronomy and Computing, 17, 1](http://adsabs.harvard.edu/abs/2016arXiv160506811H).

---

## How VARTOOLS Works

VARTOOLS operates as a **pipeline**. You supply a sequence of commands on a single invocation, and each command is applied in turn to every light curve. The output of one command is passed directly to the next. Statistics computed by each command are written to **stdout** as a space-delimited (or tab-delimited) ASCII table, one row per light curve.

```
vartools -l lc_list [global options] -cmd1 args -cmd2 args -cmd3 args ...
```

- Each light curve is read once, processed through the full command chain, and its result row is emitted.
- Commands that modify the light curve (e.g. `-clip`, `-TFA`, `-Killharm`) pass the modified version to the next command.
- Commands that compute statistics (e.g. `-rms`, `-LS`, `-BLS`) append columns to the output row without altering the light curve (unless instructed to with `correctlc`).
- Periodograms and model light curves are written to on-disk files; summary statistics go to stdout.

---

## Basic Invocation Patterns

### Single light curve

```
vartools -i path/to/lightcurve.lc [commands...]
```

Use `-i` to process a single light curve file.

### List of light curves

```
vartools -l path/to/lc_list [commands...]
```

Use `-l` to process many light curves at once. The list file contains one light curve filename per line (in the first column by default). Additional columns may carry per-light-curve parameters that commands can reference (e.g. periods, coordinates).

### Parallel processing

```
vartools -l lc_list -parallel 8 [commands...]
```

The `-parallel N` option processes up to `N` light curves simultaneously on a multi-core machine.

---

## Key Output Options

| Flag | Effect |
|------|--------|
| `-header` | Print a one-line column-name header before the data rows |
| `-tab` | Use tab-delimited (starbase) format instead of space-delimited |
| `-oneline` | Print each statistic on its own line (useful for single light curves) |
| `-numbercolumns` | Prefix each column name with its column number |
| `-quiet` | Suppress the output table entirely |
| `-redirectstats statsfile` | Write the statistics table to a file instead of stdout |
| `-headeronly` | Print the header and exit without processing any light curves |
| `-nobuffer` | Disable stdout buffering |
| `-noskipempty` | Include rows for empty light curves in the output table |
| `-basename` | Print only the basename of each light curve in the Name column |

---

## Output Column Names

Every column added to the output table by a command follows the naming convention:

```
CommandName_Stat_Peak_CmdIndex
```

where:

- **CommandName** is the name of the command (e.g. `LS`, `BLS`, `Killharm`)
- **Stat** describes the statistic (e.g. `Period`, `SNR`, `Depth`)
- **Peak** is the peak number when multiple peaks are reported (e.g. `1`, `2`)
- **CmdIndex** is the zero-based index of that command in the command chain

For example, if `-LS` is the first command (index 0) and reports two peaks:

```
LS_Period_1_0    LS_SNR_1_0    Log10_LS_Prob_1_0
LS_Period_2_0    LS_SNR_2_0    Log10_LS_Prob_2_0
```

If the same command appears more than once, the index distinguishes them (e.g. `RMS_0` from the first `-rms` call and `RMS_1` from the second). Use `-columnsuffix` to override the index with a custom string.

Column names can be used with the `-fixcolumn` parameter in later commands to pass computed values forward in the pipeline.

---

## Minimal End-to-End Example

Compute the RMS of a list of light curves, run a Lomb-Scargle period search, and write a header-labeled tab-delimited table:

```bash
vartools -l lc_list \
    -rms \
    -LS 0.1 30.0 0.1 1 0 \
    -header -tab
```

Sample output:

```
#Name  Mean_Mag_0  RMS_0  Expected_RMS_0  Npoints_0  LS_Period_1_1  Log10_LS_Prob_1_1  LS_SNR_1_1
EXAMPLES/1  10.247  0.0021  0.0010  3122  1.23588  -707.23  18.61
EXAMPLES/2  10.121  0.0500  0.0010  2800  1.23588  -490.11  16.42
...
```

To further process the light curve (e.g. whiten at the best-fit period and output the residuals):

```bash
vartools -l lc_list \
    -LS 0.1 30.0 0.1 1 0 \
    -Killharm ls 2 0 1 ./models \
    -rms \
    -header -tab
```

---

## Input Format Control

By default VARTOOLS expects light curve files with three whitespace-delimited columns: time, magnitude, and uncertainty. Use `-inputlcformat` to read non-standard formats:

```bash
vartools -i my.fits \
    -inputlcformat t:1,mag:2,err:3,xpos:4,ypos:5 \
    [commands...]
```

See `vartools -showinputlcformat` and `vartools -showinputlistformat` for format diagnostics.

---

## Further Reading

- [Light Curve Manipulation](manipulation.md) — binning, phase-folding, time conversion, filtering, flux-to-magnitude conversion, FFT, and more
- [Period Finding](period-finding.md) — Lomb-Scargle, AoV, BLS, WWZ, DFT/CLEAN, autocorrelation, and related commands
