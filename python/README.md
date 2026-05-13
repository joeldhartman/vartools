# pyvartools

Python API for the [VARTOOLS](https://www.astro.princeton.edu/~jhartman/vartools.html) light curve analysis program.

pyvartools wraps the `vartools` binary as a Python library, letting you build and run analysis pipelines entirely in Python without writing shell scripts.

## Requirements

- Python 3.8+
- `vartools` binary installed and on `PATH` (or configured via `vt.set_binary()`)
- numpy, pandas
- astropy (optional — required only for FITS I/O and `TimeSeries` interop)

## Installation

```bash
pip install -e /path/to/source/vartools/python
```

Or within HATpipe, `vartools` is already installed at `bin/vartools` and the package will find it automatically.

## Quick start

```python
import pyvartools as vt
from pyvartools import commands as cmd

# Load a light curve
lc = vt.LightCurve.from_file("EXAMPLES/2")

# Run a Lomb-Scargle periodogram
result = vt.Pipeline([cmd.LS(0.1, 10.0, 1e-3)]).run(lc)
best_period = float(result.vars["LS_Period_1_0"])
print(f"Best period: {best_period:.4f} d")

# Chain commands
pipe = vt.Pipeline([
    cmd.clip(sigclip=5.0),
    cmd.LS(0.1, 10.0, 1e-3),
])
result = pipe.run(lc)

# Retrieve the modified light curve
result = pipe.run(lc, capture_lc=True)
clipped_lc = result.lc

# Batch processing from memory
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 6)]
batch = vt.Pipeline([cmd.rms()]).run_batch(lcs, nthreads=4)
print(batch.var)  # one row per LC

# Batch processing directly from disk (no Python I/O)
batch = vt.Pipeline([cmd.rms()]).run_filelist("mylist.txt", nthreads=4)
```

## Binary discovery

pyvartools searches for the `vartools` binary in this order:

1. Path set by `vt.set_binary("/path/to/vartools")`
2. `VARTOOLS_BINARY` environment variable
3. HATpipe install: `<package>/../../../bin/vartools`
4. `vartools` on `PATH`

## LightCurve I/O

```python
# ASCII (vartools default format: time mag err [aux ...])
lc = vt.LightCurve.from_file("star.lc")

# FITS binary table
lc = vt.LightCurve.from_file("star.fits", t_col="BJD", mag_col="Mag", err_col="Err")

# From arrays
import numpy as np
lc = vt.LightCurve.from_arrays(t, mag, err, name="HAT-P-7")

# From pandas DataFrame (must have 't', 'mag', 'err' columns)
lc = vt.LightCurve.from_dataframe(df)
```

## Pipeline methods

| Method | Input | Returns |
|---|---|---|
| `pipe.run(lc)` | single `LightCurve` | `Result` |
| `pipe.run_file(path)` | path on disk | `Result` |
| `pipe.run_batch(lcs)` | list of `LightCurve` | `BatchResult` |
| `pipe.run_filelist(paths)` | list of paths or list-file | `BatchResult` |
| `pipe.run_combinelcs(groups)` | list of file-path groups | `BatchResult` |

All methods accept `capture_lc=True` to retrieve modified light curves.
Batch methods also accept `nthreads=N`, `timeout=seconds`, and `raise_on_error=False`.

All run methods accept these global vartools options:

| Parameter | Type | CLI flag | Notes |
|---|---|---|---|
| `randseed` | `int` | `-randseed N` | Reproducible random sequences |
| `skipmissing` | `bool` | `-skipmissing` | Skip missing input files silently; most useful in batch/list modes |
| `jdtol` | `float` | `-jdtol N` | Julian-date matching tolerance |
| `matchstringid` | `bool` | `-matchstringid` | Force string-based LC name matching |

### `run_combinelcs(groups, ...)`

Processes groups of files with vartools `-l … combinelcs`.  Each group is
combined into a single in-memory light curve and produces one row in
`result.vars`.  This is the natural entry point for multi-telescope stitching
workflows (e.g. with a `-stitch` user command).

```python
# Combine two per-telescope files into one in-memory LC and compute RMS
groups = [
    ["tel1/star001.lc", "tel2/star001.lc"],
    ["tel1/star002.lc", "tel2/star002.lc"],
]
result = vt.Pipeline([cmd.rms()]).run_combinelcs(groups)
print(result.vars)   # one row per group

# Track which file each observation came from
result = vt.Pipeline([cmd.rms()]).run_combinelcs(groups, lcnumvar="lcnum")
```

Key parameters:

| Parameter | Description |
|---|---|
| `groups` | Sequence of sequences of file paths; one group → one combined LC |
| `lcnumvar` | Variable name for per-observation file index (passed as `lcnumvar <name>` after `combinelcs`) |
| `delimiter` | Delimiter joining paths within a group in the list file (default `","`) |

## Available commands

```python
vt.list_commands()   # query the binary for the full list
```

Wrappers are available for all standard vartools commands, grouped by module:

- **periodicity**: `LS`, `aov`, `aov_harm`, `BLS`, `BLSFixPer`, `Phase`, `autocorrelation`, `dftclean`, `wwz`, ...
- **manipulation**: `clip`, `rms`, `stats`, `Killharm`, `Injectharm`, `sortlc`, `expr`, `FFT`, `decorr`, ...
- **fitting**: `TFA`, `TFA_SR`, `SYSREM`, `MandelAgolTransit`, `nonlinfit`, `addnoise`, ...
- **misc**: `converttime`, `binlc`, `match`, `R`, `columnsuffix`, `Raw`, ...

Each command mirrors its vartools `-help` documentation.  For example:

```python
help(cmd.LS)          # show Python docstring
cmd.LS.help()         # print vartools -LS -help output
cmd.LS.examples()     # print vartools -LS -example output
```

## Running tests

```bash
cd source/vartools/python
python3 -m pytest tests/ -p no:astropy_header -q
```
