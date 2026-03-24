# Pipeline

The `Pipeline` class is the central object in pyvartools. It holds an ordered list of command wrappers and executes them against one or more light curves. All output is returned as structured Python objects — no manual parsing required.

When `libvartoolspipeline.so` is installed (see [Installation](../install.md#library-mode-libvartoolspipelineso)), `run()` and `run_batch()` run vartools **in-process**: the command-line is parsed once at `Pipeline` construction time, and each subsequent call injects the light curve arrays directly into C memory. This eliminates subprocess startup and serialisation overhead, giving a roughly 20× speedup per call.

For cases that require output files (`save_periodogram=True`, `capture_lc=True`, `cmd.o(capture=True)`), or when a `timeout` is specified, the pipeline falls back automatically to the subprocess path.

---

## Construction

Pass a list of command objects to the constructor:

```python
import pyvartools as vt
from pyvartools import commands as cmd

pipe = vt.Pipeline([
    cmd.clip(5.0),
    cmd.LS(0.5, 10.0, 1e-3),
])
```

Commands are appended with `add()`, which returns the pipeline itself so calls can be chained:

```python
pipe = vt.Pipeline([cmd.clip(5.0)])
pipe.add(cmd.LS(0.5, 10.0, 1e-3)).add(cmd.rms())
```

---

## Run methods

### `run(lc, capture_lc=False, outdir=None, timeout=None, init_lc_vars=None, randseed=None, skipmissing=False, jdtol=None, matchstringid=False) → Result`

Run the pipeline on a single light curve held in memory.

| Parameter | Type | Description |
|-----------|------|-------------|
| `lc` | `LightCurve`, `pd.DataFrame`, or astropy `TimeSeries` | The light curve to process. DataFrames and TimeSeries objects are coerced to `LightCurve` automatically. |
| `capture_lc` | `bool` | If `True`, write the (possibly modified) output light curve to a temporary file and return it as `result.lc`. Default `False`. |
| `outdir` | `str` or `None` | Directory for command output files (e.g. periodogram files from `save_periodogram=True`). If `None` (default), a temporary directory is created and its contents are captured into `result.files` before it is deleted. |
| `timeout` | `int` or `None` | Maximum number of seconds to wait for the vartools process. Raises `RunError` if the process exceeds the limit. `None` means no limit. |
| `init_lc_vars` | `dict[str, LCVar]` or `None` | Per-observation variables to create and initialise via `-inputlcformat` col=0. See [Initialised LC variables](#initialised-lc-variables-init_lc_vars). |
| `randseed` | `int` or `None` | Pass `-randseed N` to vartools for reproducible random-number sequences. |
| `skipmissing` | `bool` | Pass `-skipmissing`; silently skip missing input files. Default `False`. |
| `jdtol` | `float` or `None` | Pass `-jdtol N` to set the Julian-date matching tolerance. |
| `matchstringid` | `bool` | Pass `-matchstringid` to force string-based LC name matching. Default `False`. |

If the `LightCurve` contains columns beyond the default three (`t`, `mag`, `err`), a `-inputlcformat` flag is constructed automatically and passed to vartools so the extra columns are accessible by name. See [Additional columns](#additional-columns-inputlcformat).

Returns a [`Result`](results.md) object.

---

### `run_file(path, capture_lc=False, outdir=None, timeout=None, columns=None, init_lc_vars=None, randseed=None, skipmissing=False, jdtol=None, matchstringid=False) → Result`

Run the pipeline on a light curve file already on disk. vartools reads the file directly — no Python I/O is performed.

| Parameter | Type | Description |
|-----------|------|-------------|
| `path` | `str` or `Path` | Path to the light curve file on disk. |
| `capture_lc` | `bool` | Capture the output light curve as `result.lc`. |
| `outdir` | `str` or `None` | Directory for command output files. |
| `timeout` | `int` or `None` | Timeout in seconds. |
| `columns` | `list[str]`, `dict`, or `None` | Column specification passed to vartools as `-inputlcformat`. See [Additional columns](#additional-columns-inputlcformat) below. |
| `init_lc_vars` | `dict[str, LCVar]` or `None` | Per-observation variables to create and initialise. See [Initialised LC variables](#initialised-lc-variables-init_lc_vars). |
| `randseed` | `int` or `None` | Pass `-randseed N` to vartools. |
| `skipmissing` | `bool` | Pass `-skipmissing`. Default `False`. |
| `jdtol` | `float` or `None` | Pass `-jdtol N`. |
| `matchstringid` | `bool` | Pass `-matchstringid`. Default `False`. |

The light curve name reported in `result.stats["Name"]` is taken from the file stem (i.e. the filename without directory or extension).

Returns a [`Result`](results.md) object.

---

### `run_batch(lcs, nthreads=1, capture_lc=False, outdir=None, timeout=None, raise_on_error=True, init_lc_vars=None, inlistvars=None, randseed=None, skipmissing=False, jdtol=None, matchstringid=False) → BatchResult`

Run the pipeline on a list of light curves in memory. All light curves are written to temporary files and processed in a single vartools invocation using `-l`.

| Parameter | Type | Description |
|-----------|------|-------------|
| `lcs` | list of `LightCurve`, `DataFrame`, or `TimeSeries` | The light curves to process. |
| `nthreads` | `int` | Number of parallel threads (vartools `-parallel` flag). Default `1`. |
| `capture_lc` | `bool` | If `True`, capture the modified output LC for each light curve and return them as `result.lcs`. |
| `outdir` | `str` or `None` | Directory for command output files. |
| `timeout` | `int` or `None` | Timeout in seconds for the whole batch. |
| `raise_on_error` | `bool` | If `True` (default), a vartools failure raises `RunError`. If `False`, the exception is stored in `result.error` and `result.stats` will be empty. |
| `init_lc_vars` | `dict[str, LCVar]` or `None` | Per-observation variables to create and initialise. See [Initialised LC variables](#initialised-lc-variables-init_lc_vars). |
| `inlistvars` | `dict[str, int \| ListVar]` or `None` | Per-star variables. For `run_batch()` the list file contains only LC paths (no extra columns), so only `col=0` (expression-initialised) `ListVar` entries are meaningful. For per-star values from file columns use `run_filelist()`. See [Per-star variables](#per-star-variables-inlistvars). |
| `randseed` | `int` or `None` | Pass `-randseed N` to vartools. |
| `skipmissing` | `bool` | Pass `-skipmissing`. Default `False`. |
| `jdtol` | `float` or `None` | Pass `-jdtol N`. |
| `matchstringid` | `bool` | Pass `-matchstringid`. Default `False`. |

Extra columns in the first `LightCurve` are used to construct a `-inputlcformat` flag automatically; all light curves in the batch are assumed to share the same column structure. See [Additional columns](#additional-columns-inputlcformat).

Returns a [`BatchResult`](results.md) object.

---

### `run_filelist(paths, nthreads=1, capture_lc=False, outdir=None, timeout=None, raise_on_error=True, columns=None, init_lc_vars=None, inlistvars=None, randseed=None, skipmissing=False, jdtol=None, matchstringid=False) → BatchResult`

Run the pipeline on a collection of light curve files on disk. No Python I/O is performed — vartools reads the files directly. This is the most efficient method for large surveys.

| Parameter | Type | Description |
|-----------|------|-------------|
| `paths` | `str`, `Path`, or list of `str`/`Path` | Either a path to an existing vartools list file (one LC file path per line), or a Python list of individual LC file paths. If a list is given, a temporary list file is written automatically. |
| `nthreads` | `int` | Number of parallel threads. |
| `capture_lc` | `bool` | Capture output LCs as `result.lcs`. |
| `outdir` | `str` or `None` | Directory for command output files. |
| `timeout` | `int` or `None` | Timeout in seconds. |
| `raise_on_error` | `bool` | If `False`, errors are stored in `result.error` rather than raised. |
| `columns` | `list[str]`, `dict`, or `None` | Column specification passed to vartools as `-inputlcformat`. See [Additional columns](#additional-columns-inputlcformat) below. |
| `init_lc_vars` | `dict[str, LCVar]` or `None` | Per-observation variables to create and initialise. See [Initialised LC variables](#initialised-lc-variables-init_lc_vars). |
| `inlistvars` | `dict[str, int \| ListVar]` or `None` | Per-star variables read from list file columns. See [Per-star variables](#per-star-variables-inlistvars). |
| `randseed` | `int` or `None` | Pass `-randseed N` to vartools. |
| `skipmissing` | `bool` | Pass `-skipmissing`. Default `False`. |
| `jdtol` | `float` or `None` | Pass `-jdtol N`. |
| `matchstringid` | `bool` | Pass `-matchstringid`. Default `False`. |

Returns a [`BatchResult`](results.md) object.

---

### `run_combinelcs(groups, nthreads=1, capture_lc=False, outdir=None, timeout=None, raise_on_error=True, columns=None, init_lc_vars=None, inlistvars=None, lcnumvar=None, delimiter=",", randseed=None, skipmissing=False, jdtol=None, matchstringid=False) → BatchResult`

Run the pipeline using vartools `-l … combinelcs` mode. Each entry in *groups* is a list of file paths that vartools combines into a single in-memory light curve before passing it to the command chain. The result contains one row in `batch.stats` per group.

This is the natural entry point for multi-telescope stitching workflows — for example, merging per-telescope files with a `-stitch` user command before running a period search.

| Parameter | Type | Description |
|-----------|------|-------------|
| `groups` | list of lists of `str`/`Path` | Each inner list is one group of files to combine. Files within a group are joined by `delimiter` to form one line in the vartools list file. |
| `nthreads` | `int` | Number of parallel threads. |
| `capture_lc` | `bool` | Capture the combined output LC for each group as `result.lcs`. |
| `outdir` | `str` or `None` | Directory for command output files. |
| `timeout` | `int` or `None` | Timeout in seconds. |
| `raise_on_error` | `bool` | If `False`, errors are stored in `result.error` rather than raised. |
| `columns` | `list[str]`, `dict`, or `None` | Column specification passed to vartools as `-inputlcformat`. |
| `init_lc_vars` | `dict[str, LCVar]` or `None` | Per-observation variables to create and initialise. |
| `inlistvars` | `dict[str, int \| ListVar]` or `None` | Per-star variables. See [Per-star variables](#per-star-variables-inlistvars). |
| `lcnumvar` | `str` or `None` | If given, pass `lcnumvar <name>` after `combinelcs` in the `-l` flag. vartools creates a per-observation integer variable recording which source file each point came from. |
| `delimiter` | `str` | Delimiter used to join paths within a group in the list file. Default `","` (the vartools `combinelcs` default). |
| `randseed` | `int` or `None` | Pass `-randseed N` to vartools. |
| `skipmissing` | `bool` | Pass `-skipmissing`. Default `False`. |
| `jdtol` | `float` or `None` | Pass `-jdtol N`. |
| `matchstringid` | `bool` | Pass `-matchstringid`. Default `False`. |

PerLC array parameters are not supported with `run_combinelcs()` — a `ValueError` is raised if any command in the pipeline carries a `PerLC` attribute.

Returns a [`BatchResult`](results.md) object.

#### Example — combine two per-telescope files and compute RMS

```python
import pyvartools as vt
from pyvartools import commands as cmd

groups = [
    ["tel1/star001.lc", "tel2/star001.lc"],
    ["tel1/star002.lc", "tel2/star002.lc"],
]

batch = vt.Pipeline([cmd.rms()]).run_combinelcs(groups)
print(batch.stats)   # one row per group
```

#### Example — track which file each point came from

```python
batch = vt.Pipeline([cmd.rms()]).run_combinelcs(
    groups,
    lcnumvar="lcnum",   # vartools creates an integer variable "lcnum" per observation
)
```

#### Example — multi-telescope stitch + period search

```python
from pyvartools import UserCommand, commands as cmd

Stitch = vt.load_userlib("/path/to/stitch.so", name="stitch")

batch = vt.Pipeline([
    Stitch("mag err mask lcnum median"),
    cmd.LS(0.5, 10.0, 1e-3, npeaks=1),
]).run_combinelcs(groups, lcnumvar="lcnum", nthreads=8)

print(batch.stats[["Name", "LS_Period_1_0"]])
```

---

## Examples

### Single light curve

```python
import pyvartools as vt
from pyvartools import commands as cmd

lc = vt.LightCurve.from_file("EXAMPLES/2")
pipe = vt.Pipeline([cmd.clip(5.0), cmd.LS(0.5, 10.0, 1e-3)])

result = pipe.run(lc)
print(result.stats["LS_Period_1_1"])
```

To get the modified light curve back (e.g. after sigma clipping):

```python
result = pipe.run(lc, capture_lc=True)
modified_lc = result.lc  # LightCurve with clipped points removed
```

### From disk

```python
result = pipe.run_file("EXAMPLES/2")
print(result.stats["LS_Period_1_1"])
```

### Batch from memory

```python
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 6)]
batch = pipe.run_batch(lcs, nthreads=4)
print(batch.stats)          # pd.DataFrame, one row per LC
print(batch.stats.columns)  # all output column names
```

### Batch from disk (most efficient for large surveys)

When you already have a list file or a directory of light curves, `run_filelist` avoids all Python-side I/O:

```python
# From an existing vartools list file
batch = pipe.run_filelist("EXAMPLES/lc_list", nthreads=4)

# Or supply an explicit Python list of paths
batch = pipe.run_filelist([f"EXAMPLES/{i}" for i in range(1, 6)], nthreads=4)

# Save the statistics table
batch.stats.to_csv("EXAMPLES/results.csv", index=False)
```

### Error handling in batch

```python
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 6)]
batch = pipe.run_batch(lcs, raise_on_error=False)
if batch.ok:
    print(batch.stats)
else:
    print(f"Pipeline failed: {batch.error}")
```

---

## Additional columns (`-inputlcformat`)

By default vartools treats a light curve file as having three columns: time (column 1), magnitude (column 2), and uncertainty (column 3). To make extra columns available to commands under a named variable, vartools accepts a `-inputlcformat` flag of the form:

```
-inputlcformat t:1,mag:2,err:3,airmass:4,xpos:5
```

pyvartools handles this automatically:

- **`run()` and `run_batch()`** — when the `LightCurve` object contains columns beyond `t`, `mag`, `err`, the format string is built from the DataFrame column list and passed to vartools automatically.
- **`run_file()` and `run_filelist()`** — the file is read directly by vartools, so the column layout must be supplied by the caller via the `columns` parameter.

### Specifying columns for disk-based runs

**List form** — variable names in column order:

```python
# Columns: t=1, mag=2, err=3
result = pipe.run_file(
    "EXAMPLES/2",
    columns=["t", "mag", "err"],
)
```

**Dict form with integer column numbers** — explicit mapping:

```python
result = pipe.run_file(
    "EXAMPLES/2",
    columns={"t": 1, "mag": 2, "err": 3},
)
```

**Dict form with FITS column names** — for FITS binary tables:

```python
result = pipe.run_file(
    "EXAMPLES/example.fits",
    columns={"t": "BJD", "mag": "Mag", "err": "Err"},
)
```

### Automatic discovery in memory-based runs

A `-inputlcformat` flag is emitted automatically whenever the `LightCurve`
column layout differs from the default `[t, mag, err]` — whether that means
extra columns are present, standard columns are absent, or the order is
non-standard.

```python
import numpy as np
import pyvartools as vt
from pyvartools import commands as cmd

t       = np.linspace(0, 30, 500)
mag     = 10.0 + 0.1 * np.sin(2 * np.pi * t / 2.3)
err     = np.full(500, 0.01)
airmass = 1.0 + 0.5 * np.abs(np.sin(np.pi * t / 15.0))

# Extra column
lc = vt.LightCurve.from_arrays(t, mag, err, aux={"airmass": airmass})
# → -inputlcformat t:1,mag:2,err:3,airmass:4  (auto-generated)
result = vt.Pipeline([cmd.rms()]).run(lc)

# Only t and mag — vartools receives -inputlcformat t:1,mag:2
lc2 = vt.LightCurve.from_arrays(t=t, mag=mag)
result2 = vt.Pipeline([cmd.rms()]).run(lc2)

# No standard columns at all — only custom vectors
phase_arr = (t % 2.3) / 2.3
flux_arr  = 1.0 - 0.01 * np.sin(2 * np.pi * phase_arr)
lc3 = vt.LightCurve.from_arrays(aux={"phase": phase_arr, "flux": flux_arr})
# → -inputlcformat phase:1,flux:2
```

The same applies to `run_batch()` — all light curves in the batch should share
the same column structure, as a single `-inputlcformat` is inferred from the
first light curve and applied to the whole batch.

---

## Initialised LC variables (`init_lc_vars`)

vartools supports a special column number `0` in `-inputlcformat`. Instead of reading a value from a file column, vartools **creates a per-observation variable** and initialises it from an analytic expression. This is useful for creating mask flags, synthetic indices, or any per-point derived quantity before the pipeline begins.

The expression is evaluated once per observation. The special variable `NR` holds the 0-based observation index.

### `LCVar`

```python
from pyvartools import LCVar

LCVar(type="double", init="0")
```

| Attribute | Type | Default | Description |
|-----------|------|---------|-------------|
| `type` | `str` | `"double"` | Variable type: `"double"`, `"float"`, `"int"`, `"long"`, `"short"`, `"string"`, `"char"`, or `"utc"`. |
| `init` | `str` | `"0"` | Initialisation expression. May reference `NR` (0-based obs index) or other variables already defined by earlier `-inputlcformat` entries. |

Pass a `dict[str, LCVar]` as `init_lc_vars` to any run method. The dictionary key is the variable name used inside vartools commands.

!!! note
    Supplying `init_lc_vars` **replaces** vartools' implicit default column mapping. pyvartools prepends `t:1,mag:2,err:3` automatically so the standard columns remain mapped. If your light curve has a non-standard column layout, also supply the `columns` parameter (`run_file`/`run_filelist`) or use a `LightCurve` with the correct column names (`run`/`run_batch`).

### Examples

**Create an integer mask variable, initially all zeros (unmasked):**

```python
import pyvartools as vt
from pyvartools import LCVar, commands as cmd

lc = vt.LightCurve.from_file("EXAMPLES/2")

result = vt.Pipeline([cmd.rms()]).run(
    lc,
    init_lc_vars={"mymask": LCVar(type="int", init="0")},
)
print(result.stats["RMS_0"])
```

**Index each observation (0-based) using `NR`:**

```python
result = vt.Pipeline([cmd.rms()]).run(
    lc,
    init_lc_vars={"obs_idx": LCVar(type="double", init="NR")},
)
```

**Use a per-observation flag to inflate errors for early observations (t < 10):**

```python
result = vt.Pipeline([
    cmd.expr("err = err * (1 + 9*early_mask)"),   # inflate errors 10× for t<10
    cmd.rms(),
]).run(
    lc,
    init_lc_vars={"early_mask": LCVar(type="int", init="t<10")},
)
print(result.stats["RMS_1"])
```

---

## Per-star variables (`inlistvars`)

vartools `-inlistvars` defines **per-star (scalar) variables** read from extra columns in the input list file, one value per light curve. These variables are accessible to commands like `LS`, `aov`, and `BLS` via the `var` keyword (e.g. `minp="myperiod"` in `cmd.LS`).

### `ListVar`

```python
from pyvartools import ListVar

ListVar(col=2, type="double", init=None, combinelc=False)
```

| Attribute | Type | Default | Description |
|-----------|------|---------|-------------|
| `col` | `int` | *(required)* | Column number in the list file (1-based). Use `0` to initialise from an expression instead of reading from a column. |
| `type` | `str` | `"double"` | Variable type: `"double"`, `"float"`, `"int"`, `"long"`, `"short"`, `"string"`, `"char"`, or `"utc"`. |
| `init` | `str` or `None` | `None` | Initialisation expression used when `col=0`. The special variable `NF` holds the 0-based line number. |
| `combinelc` | `bool` | `False` | If `True`, the same value is applied to all light curves that share the same LC name (used when combining multiple observations into a single LC). |

As a shorthand, `int` values are accepted in place of `ListVar` objects when only a column number is needed:

```python
inlistvars={"myperiod": 2}   # equivalent to ListVar(col=2)
```

### With `run_filelist()`

`run_filelist()` is the primary method for using `inlistvars`. The list file may contain extra columns beyond the LC file path. vartools reads those columns and maps them to the specified variable names.

**Example list file (`EXAMPLES/lc_list_periods`):**

```
EXAMPLES/1 0.573
EXAMPLES/2 0.412
EXAMPLES/3 1.891
```

**Use the period from column 2 as the minimum search period for LS:**

```python
import pyvartools as vt
from pyvartools import ListVar, commands as cmd

pipe = vt.Pipeline([
    cmd.LS("myperiod", 100.0, 0.1),
])
batch = pipe.run_filelist(
    "EXAMPLES/lc_list_periods",
    inlistvars={"myperiod": ListVar(col=2)},
)
print(batch.stats[["Name", "LS_Period_1_0"]])
```

Here `minp="myperiod"` is a bare identifier string, so pyvartools emits `var myperiod` — a per-star variable read from the list file.

**Expression-initialised per-star variable (`col=0`) using the line number `NF`:**

```python
batch = pipe.run_filelist(
    "EXAMPLES/lc_list_periods",
    inlistvars={
        "myperiod": ListVar(col=2),
        "lc_idx": ListVar(col=0, type="double", init="NF"),
    },
)
```

### With `run_batch()`

`run_batch()` builds a temporary list file containing only LC file paths — there are no extra columns available. Therefore `inlistvars` with `run_batch()` is only useful for `col=0` expression-initialised variables. For per-star values from a file, write the list file yourself and use `run_filelist()`.

### Reserved variable names

The names `t`, `mag`, `err`, and `id` are reserved by vartools and cannot be used as `inlistvars` variable names.

---

## Per-LC array parameters

Some numeric parameters on period-search commands accept a **per-LC array** — a different value for each light curve in a batch — as a direct alternative to a single fixed scalar. This lets you, for example, run LS on a hundred stars each with individually tailored search bounds, without writing a custom list file.

### Supported types

Any of the following may be passed for a supported parameter:

| Type | Example |
|------|---------|
| `numpy.ndarray` (1-D) | `np.array([0.1, 0.5, 0.2])` |
| `pyvartools.PerLC` | `PerLC([0.1, 0.5, 0.2])` |
| `pandas.Series` | `pd.Series([0.1, 0.5, 0.2])` |

`PerLC` is a thin wrapper that tags an array explicitly as a per-LC sequence, making intent clear in code. It behaves identically to a numpy array at runtime.

```python
from pyvartools import PerLC
```

### Supported (command, parameter) pairs

Per-LC arrays are supported wherever vartools itself accepts the `var varname` syntax. The complete list, confirmed against the installed binary:

| Command | Supported parameters |
|---------|---------------------|
| `LS` | `minp`, `maxp`, `subsample` |
| `aov` | `minp`, `maxp`, `subsample`, `finetune`, `nbin` |
| `aov_harm` | `nharm`, `minp`, `maxp`, `subsample`, `finetune` |
| `BLS` | `minper`, `maxper`, `rmin`, `rmax`, `qmin`, `qmax`, `stellar_density`, `min_exp_dur_frac`, `max_exp_dur_frac`, `nbins`, `subsample`, `nfreq`, `df` |

Parameters not listed here (e.g., `Phase.T0`, `restricttimes.minJD`, all `Injecttransit` parameters) do **not** support per-LC arrays — vartools has no `var`-keyword slot for those values.

### How it works

When `run_batch()` detects a per-LC array on any command parameter, it:

1. Adds the values as an extra column in the temporary list file it writes for vartools.
2. Passes `-inlistvars name:col` to vartools so the column is loaded as a named per-star variable.
3. Replaces the parameter value in the command with the variable name, emitting `var name` in the vartools CLI.

This is all automatic — no manual list file or `-inlistvars` specification is needed.

### Constraints

- **`run_batch()` and `run_filelist()` (list of paths) only.** `run()` and `run_file()` process a single light curve and raise `ValueError` if any per-LC array is present.
- **`run_filelist()` with a pre-built list file (a single path string) is not supported.** pyvartools cannot safely append extra columns to a list file it did not create. Pass a Python list of paths instead and let the pipeline write the list file.
- The array length must equal the number of light curves. A `ValueError` with a clear message is raised if it does not.

### Examples

**Different period search range for each LC:**

```python
import numpy as np
import pyvartools as vt
from pyvartools import commands as cmd

lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in [2, 3, 4, 5, 6]]

# Each star gets its own [minp, maxp] bounds
result = vt.Pipeline([
    cmd.LS(
        minp=np.array([0.1, 0.5, 0.1, 0.2, 0.3]),
        maxp=np.array([5.0, 10.0, 8.0, 6.0, 7.0]),
        subsample=0.001,
        npeaks=1,
    )
]).run_batch(lcs, nthreads=4)

print(result.stats[["Name", "LS_Period_1_0"]])
```

**Using the `PerLC` wrapper for readability:**

```python
from pyvartools import PerLC, commands as cmd

pipe = vt.Pipeline([
    cmd.aov(
        minp=PerLC([0.1, 0.5, 0.1, 0.2, 0.3]),
        maxp=10.0,               # scalar: same for all LCs
        subsample=0.001,
        finetune=2,
        npeaks=1,
    )
])
result = pipe.run_batch(lcs)
print(result.stats[["Name", "Period_1_0", "AOV_1_0"]])
```

**BLS with per-LC period bounds:**

```python
result = vt.Pipeline([
    cmd.BLS(
        minper=np.array([0.1, 0.5, 0.1, 0.2, 0.3]),
        maxper=np.array([5.0, 10.0, 8.0, 6.0, 7.0]),
        rmin=0.01,
        rmax=0.5,
        nbins=200,
        nfreq=5000,
        npeaks=1,
    )
]).run_batch(lcs)
print(result.stats[["Name", "BLS_Period_1_0", "BLS_SDE_1_0"]])
```

**Mixing per-LC and scalar parameters** — scalar values are broadcast to all light curves as usual:

```python
# minp varies per LC; maxp and subsample are fixed for all
pipe = vt.Pipeline([
    cmd.LS(
        minp=PerLC([0.1, 0.3, 0.05]),
        maxp=20.0,
        subsample=0.001,
    )
])
```

**Multiple per-LC parameters** — multiple parameters can vary simultaneously:

```python
# Both minp and maxp differ per LC
pipe = vt.Pipeline([
    cmd.LS(
        minp=np.array([0.1, 0.5, 0.1]),
        maxp=np.array([5.0, 20.0, 10.0]),
        subsample=np.array([0.001, 0.0005, 0.001]),
    )
])
```

---

## Output files

Commands that produce auxiliary output files (periodograms, fitted model curves, coefficient tables, etc.) expose a `save_*` parameter. Each such parameter accepts four modes controlled by a `bool`, a directory path string, or an [`Output`](commands.md#auxiliary-output-files) object:

| Value | Mode | Written to disk? | Captured in `result.files`? |
|-------|------|-----------------|------------------------------|
| `False` (default) | suppress | no | no |
| `True` | temp, capture | pipeline-managed temp dir (auto-deleted) | **yes** |
| `"/path/to/dir"` | disk only | that directory | no |
| `Output("/path/to/dir", capture=True)` | disk + capture | that directory | **yes** |

```python
import pyvartools as vt
from pyvartools import Output, commands as cmd

lc = vt.LightCurve.from_file("EXAMPLES/2")

pipe = vt.Pipeline([
    cmd.LS(0.5, 10.0, 1e-3, save_periodogram=True),   # Mode 1: temp dir, capture
])
result = pipe.run(lc)
pgram = result.files["LS_periodogram_0"]   # pd.DataFrame
pgram.to_csv("EXAMPLES/periodogram.csv", index=False)

# Mode 3: write to EXAMPLES/OUTDIR1/, do NOT capture into Python
pipe3 = vt.Pipeline([
    cmd.LS(0.5, 10.0, 1e-3, save_periodogram="EXAMPLES/OUTDIR1"),
])
result3 = pipe3.run(lc)
# EXAMPLES/OUTDIR1/stdin.ls is on disk; result3.files has no "LS_periodogram_0"

# Mode 2: write to EXAMPLES/OUTDIR1/ AND capture into Python
pipe2 = vt.Pipeline([
    cmd.LS(0.5, 10.0, 1e-3,
           save_periodogram=Output("EXAMPLES/OUTDIR1", capture=True)),
])
result2 = pipe2.run(lc)
pgram2 = result2.files["LS_periodogram_0"]   # read from EXAMPLES/OUTDIR1/stdin.ls
```

### Key naming

Captured files appear in `result.files` under keys of the form **`"{CommandName}_{logical}_{idx}"`**, where `idx` is the zero-based position of the command in the pipeline. This ensures that two period-search commands (or any two commands sharing a logical name) never overwrite each other's results.

All commands with `save_*` parameters and their corresponding key patterns (command at pipeline index `N`):

| Command | `save_*` param | Key in `result.files` | Description |
|---------|---------------|-----------------------|-------------|
| `LS` | `save_periodogram` | `"LS_periodogram_N"` | Frequency vs. LS power |
| `aov` | `save_periodogram` | `"aov_periodogram_N"` | Frequency vs. AOV statistic |
| `aov_harm` | `save_periodogram` | `"aov_harm_periodogram_N"` | Frequency vs. harmonic AOV statistic |
| `BLS` | `save_periodogram` | `"BLS_periodogram_N"` | BLS power spectrum |
| `BLS` | `save_model` | `"BLS_model_N"` | Best-fit transit model |
| `BLSFixPer` | `save_model` | `"BLSFixPer_model_N"` | Best-fit model at fixed period |
| `BLSFixDurTc` | `save_periodogram` | `"BLSFixDurTc_periodogram_N"` | BLS power spectrum |
| `BLSFixDurTc` | `save_model` | `"BLSFixDurTc_model_N"` | Best-fit model |
| `BLSFixPerDurTc` | `save_model` | `"BLSFixPerDurTc_model_N"` | Best-fit model at fixed parameters |
| `autocorrelation` | `save_result` | `"autocorrelation_result_N"` | Lag vs. autocorrelation value |
| `dftclean` | `save_dspec` | `"dftclean_dspec_N"` | Dirty spectrum |
| `dftclean` | `save_wfunc` | `"dftclean_wfunc_N"` | Window function |
| `dftclean` | `save_cspec` | `"dftclean_cspec_N"` | CLEAN spectrum |
| `wwz` | `save_transform` | `"wwz_transform_N"` | Full WWZ time-frequency map |
| `wwz` | `save_maxtransform` | `"wwz_maxtransform_N"` | WWZ maximum power vs. time |
| `Killharm` | `save_model` | `"Killharm_model_N"` | Fitted harmonic series |
| `Injectharm` | `save_model` | `"Injectharm_model_N"` | Injected signal model |
| `Injecttransit` | `save_model` | `"Injecttransit_model_N"` | Injected transit model |
| `linfit` | `save_model` | `"linfit_model_N"` | Fitted linear model curve |
| `nonlinfit` | `save_model` | `"nonlinfit_model_N"` | Fitted non-linear model curve |
| `decorr` | `save_model` | `"decorr_model_N"` | Fitted decorrelation model |
| `TFA` | `save_coeffs` | `"TFA_coeffs_N"` | TFA trend coefficients |
| `TFA` | `save_model` | `"TFA_model_N"` | Reconstructed TFA trend |
| `TFA_SR` | `save_coeffs` | `"TFA_SR_coeffs_N"` | TFA-SR trend coefficients |
| `TFA_SR` | `save_model` | `"TFA_SR_model_N"` | Reconstructed TFA-SR trend |
| `SYSREM` | `save_model` | `"SYSREM_model_N"` | SYSREM model |
| `SYSREM` | `save_trends` | `"SYSREM_trends_N"` | SYSREM trend vectors |
| `MandelAgolTransit` | `save_model` | `"MandelAgolTransit_model_N"` | Fitted transit model LC |
| `MandelAgolTransit` | `save_phcurve` | `"MandelAgolTransit_phcurve_N"` | Phase-folded transit model |
| `MandelAgolTransit` | `save_jdcurve` | `"MandelAgolTransit_jdcurve_N"` | Time-domain transit model |
| `SoftenedTransit` | `save_model` | `"SoftenedTransit_model_N"` | Fitted trapezoidal transit model |
| `Starspot` | `save_model` | `"Starspot_model_N"` | Fitted starspot model |
| `microlens` | `save_model` | `"microlens_model_N"` | Fitted microlensing model |
| `findblends` | `save_matches` | `"findblends_matches_N"` | Matched blend candidate list |

### `autocorrelation` — special case

`autocorrelation` always writes its output file (the vartools CLI provides no option to suppress it). `save_result=False` suppresses Python capture but the file is still written to a temp directory and discarded after the run.

### Batch runs

For batch runs, `batch.files[key]` is a **list** of DataFrames — one per light curve (or `None` where the file was not produced):

```python
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 6)]
batch = pipe.run_batch(lcs)
for i, pgram in enumerate(batch.files["LS_periodogram_0"]):
    pgram.to_csv(f"EXAMPLES/pgram_{i}.csv", index=False)
```

### The `Output` class

```python
from pyvartools import Output

Output(path=None, capture=True)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `path` | `str` or `None` | Directory to write the file to. `None` uses the pipeline-managed temp directory. |
| `capture` | `bool` | Whether to include the file in `result.files`. Default `True`. |

See [Auxiliary output files](commands.md#auxiliary-output-files) in the command reference for full examples.

---

## The `Result` and `BatchResult` objects

See the [Results](results.md) reference page for full attribute documentation. A brief summary:

**`Result`** (single LC):

- `result.stats` — `pd.Series` indexed by vartools output variable names
- `result.lc` — `LightCurve` or `None` (only if `capture_lc=True`)
- `result.files` — `dict[str, pd.DataFrame]`

**`BatchResult`** (multiple LCs):

- `batch.stats` — `pd.DataFrame`, one row per light curve
- `batch.lcs` — list of `LightCurve` or `None` (only if `capture_lc=True`)
- `batch.files` — `dict[str, list[pd.DataFrame]]`
- `batch.ok` — `True` if the run completed without error
- `batch.error` — `RunError` or `None`

---

## Performance: library mode

When `libvartoolspipeline.so` is on the library search path, pyvartools runs
vartools **in-process** rather than spawning a new subprocess for every call.
The pipeline is initialised once (command-line parsing, memory allocation) and
light curve arrays are injected directly into C memory on each subsequent call.

```python
import pyvartools as vt
from pyvartools import commands as cmd
import numpy as np

lc = vt.LightCurve.from_file("EXAMPLES/2")
pipe = vt.Pipeline([cmd.rms()])

# First call: initialises the in-process pipeline (~same cost as subprocess)
result = pipe.run(lc)

# Subsequent calls: inject arrays, skip subprocess overhead (~20× faster)
for lc in [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(3, 8)]:
    result = pipe.run(lc)
```

### Enabling and disabling library mode

Library mode is **on by default** when `libvartoolspipeline.so` can be found.
It falls back silently to subprocess mode when:

- the library is not installed
- `VARTOOLS_USE_LIBRARY=0` is set in the environment
- `timeout` is specified (library mode has no built-in timeout mechanism)
- `capture_lc=True` or output files are requested

To force subprocess mode globally:

```bash
export VARTOOLS_USE_LIBRARY=0
```

### Non-standard library path

If `libvartoolspipeline.so` is not on the default search path, point pyvartools
to it explicitly:

=== "Python"

    ```python
    import pyvartools as vt
    vt.set_library("/path/to/libvartoolspipeline.so")
    ```

=== "Environment variable"

    ```bash
    export VARTOOLS_LIBRARY=/path/to/libvartoolspipeline.so
    ```

### When library mode is not used

The pipeline falls back to the subprocess path transparently in any of these
cases:

| Situation | Reason |
|-----------|--------|
| `libvartoolspipeline.so` not found | Library not installed |
| `VARTOOLS_USE_LIBRARY=0` | Explicitly disabled |
| `timeout` is set | Library mode has no timeout support |
| `capture_lc=True` | Modified LC requires vartools to write a file |
| Any command has `save_*=True` | Output files require a working directory |
| `cmd.o(capture=True)` | Output LC capture requires filesystem |
| `init_lc_vars` is set | Library pipeline argv does not include `-inputlcformat` col=0 entries |
| Any global option is non-default (`randseed`, `jdtol`, `skipmissing`, `matchstringid`) | Global CLI flags are not threaded through the in-process library |

---

## Error handling

`RunError` is raised (or stored) when:

- vartools exits with a non-zero status
- vartools cannot be found on `PATH`
- The `timeout` is exceeded

```python
from pyvartools.results import RunError

try:
    result = pipe.run(lc, timeout=30)
except RunError as e:
    print(f"vartools failed: {e}")
```
