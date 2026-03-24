# pyvartools — Python API Overview

**pyvartools** is the official Python package for VARTOOLS. It feeds light
curves to the vartools analysis engine and returns the results as structured
Python objects — pandas Series and DataFrames for statistics, `LightCurve`
objects for photometric data, and `pd.DataFrame` instances for auxiliary
output files such as periodograms.

When `libvartoolspipeline.so` is installed alongside the binary (which `make
install` does automatically), `Pipeline.run()` and `run_batch()` run
**in-process** — no subprocess is spawned, no I/O is performed — giving a
roughly 20× speedup per call. If the library is not present the package falls
back transparently to the original subprocess path.

Every command available from the command line has a corresponding Python class
in `pyvartools.commands`. Pipelines are built by composing these command
objects and then executed with one of several run methods.

---

## Installation

See the [Installation](../install.md) page for the full procedure. In brief:

```bash
pip install pyvartools      # from PyPI
# or, to install from the VARTOOLS source tree:
cd python && pip install -e .
```

pyvartools requires a working `vartools` binary on `PATH` (or configured via
`pyvartools.config.set_binary_path()`).

---

## Core Concepts

### LightCurve

A `LightCurve` wraps a pandas DataFrame. The three standard columns — time
(`t`), magnitude (`mag`), and magnitude uncertainty (`err`) — are treated
specially when present but **none are required**. Any column layout is
accepted; pyvartools constructs a `-inputlcformat` flag automatically so
vartools knows the mapping. Additional columns (e.g., airmass, *x*/*y*
position) can be supplied via the `aux` parameter of `from_arrays`, or are
loaded automatically from multi-column ASCII and FITS files.

See [LightCurve](lightcurve.md) for the full API.

### Pipeline

A `Pipeline` holds an ordered sequence of `VartoolsCommand` objects.
When a run method is called the pipeline translates the command objects into
a VARTOOLS command-line string, invokes the binary, and parses the output.

```python
import pyvartools as vt
from pyvartools import commands as cmd

pipe = vt.Pipeline([
    cmd.LS(minp=0.5, maxp=10.0, subsample=0.1),
    cmd.BLS(minper=0.5, maxper=10.0, qmin=0.01, qmax=0.1, nfreq=20000, nbins=200),
])
```

A pipeline is reusable: the same `Pipeline` instance can be called on many
different light curves.

### Result and BatchResult

`Result` is returned by single-LC runs. `BatchResult` is returned by batch
runs. Both are documented in detail on the [Results](results.md) page.

- **`Result.stats`** — `pd.Series` of output statistics, indexed by VARTOOLS
  column names (e.g., `'LS_Period_1_0'`).
- **`BatchResult.stats`** — `pd.DataFrame` with one row per light curve.
- Both carry an optional `.lc` / `.lcs` attribute (set when `capture_lc=True`)
  and a `.files` dict of auxiliary output DataFrames (e.g., periodogram tables).

### Commands

Each VARTOOLS command is a Python class in `pyvartools.commands`. Constructor
arguments correspond to the command's positional and keyword parameters. For
example:

```python
cmd.LS(minp=0.5, maxp=10.0, subsample=0.01)
cmd.BLS(minper=0.5, maxper=10.0, qmin=0.01, qmax=0.1, nfreq=20000, nbins=200)
cmd.Killharm(period="ls", nharm=3, nsubharm=0)
```

Commands that write output files (periodograms, model light curves, etc.)
accept an `outdir` argument; pyvartools captures these files and returns them
in `Result.files` / `BatchResult.files`.

User-developed extension commands (compiled `.so` libraries) are also
supported through `UserCommand`, `load_userlib()`, and `discover_userlibs()`.
See [User Extension Commands](commands.md#user-extension-commands) for the
full API.

---

## Quick Example

```python
import numpy as np
import pyvartools as vt
from pyvartools import commands as cmd

# --- Build a synthetic light curve ---
rng = np.random.default_rng(42)
t   = np.sort(rng.uniform(0, 30, 500))
mag = 10.0 + 0.05 * np.sin(2 * np.pi * t / 2.3) + rng.normal(0, 0.005, 500)
err = np.full(500, 0.005)
lc  = vt.LightCurve.from_arrays(t, mag, err, name="my_star")

# --- Define the pipeline ---
pipe = vt.Pipeline([
    cmd.LS(minp=0.5, maxp=10.0, subsample=0.1),
])

# --- Run on a single light curve ---
result = pipe.run(lc)
print(result.stats["LS_Period_1_0"])    # best-fit period
print(result.stats["Log10_LS_Prob_1_0"])  # log10 false-alarm probability

# --- Run on a batch ---
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 6)]
batch = pipe.run_batch(lcs)
print(batch.stats[["Name", "LS_Period_1_0"]].to_string())
```

---

## Run Methods

| Method | Input | Returns | Notes |
|---|---|---|---|
| `pipe.run(lc)` | `LightCurve` | `Result` | In-process (library mode) when `libvartoolspipeline.so` is available; otherwise passes LC via stdin. No temp files written on the Python side. |
| `pipe.run_file(path)` | file path (str or `Path`) | `Result` | Single LC read directly by the vartools binary; fastest option for large files on disk. |
| `pipe.run_batch(lcs)` | list of `LightCurve` | `BatchResult` | In-process (library mode) when available; otherwise writes LCs to a temporary list and feeds them to vartools in one subprocess call. |
| `pipe.run_filelist(paths)` | list of paths **or** path to a list-file | `BatchResult` | The binary reads all files directly. Pass a `str` naming a pre-built list file, or a Python list of path strings. |
| `pipe.run_combinelcs(groups)` | list of file-path groups | `BatchResult` | Uses vartools `-l … combinelcs` to merge multiple files into one in-memory LC per group. Natural entry point for multi-telescope stitching workflows. |

All run methods accept `capture_lc=True` to collect the final light curve
state in the returned `Result.lc` / `BatchResult.lcs`, and `nthreads=N` (batch
methods) to forward `-parallel N` to vartools. `run_file()` and
`run_filelist()` also accept a `columns` parameter to describe the on-disk
column layout; `run()` and `run_batch()` discover the column layout
automatically from the `LightCurve` object.

All run methods also accept four global vartools options: `randseed` (`int`),
`skipmissing` (`bool`), `jdtol` (`float`), and `matchstringid` (`bool`). These
emit the corresponding CLI flags (`-randseed`, `-skipmissing`, `-jdtol`,
`-matchstringid`) before the command chain. Any non-default value forces
subprocess mode even when library mode is active.

---

## Output Column Naming

pyvartools uses the same column-naming convention as the VARTOOLS CLI. Each
statistic name is formed from the command name, a descriptor, and a numeric
suffix that encodes the command's position in the pipeline (0-indexed):

```
LS_Period_1_0    → LS command, 2nd peak, command index 0
BLS_Period_1_1   → BLS command, 2nd peak, command index 1
```

When you call the same command more than once, the suffix changes and the names
remain unique. To give a predictable, human-readable suffix instead of a
positional number, insert a `columnsuffix()` call before the command:

```python
from pyvartools.commands import columnsuffix, LS

pipe = vt.Pipeline([
    columnsuffix("search1"),
    LS(minp=0.5, maxp=5.0, subsample=0.01),
    columnsuffix("search2"),
    LS(minp=5.0, maxp=50.0, subsample=0.1),
])
result = pipe.run(lc)
print(result.stats["LS_Period_1_search1"])
print(result.stats["LS_Period_1_search2"])
```

The suffix applies only to the immediately following command; `columnsuffix()`
must be repeated for each command whose name you want to control.
