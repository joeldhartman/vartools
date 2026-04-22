# LightCurve

`pyvartools.LightCurve` is the primary data container for a single photometric
time series. Internally it stores data in a pandas DataFrame. The three
standard columns — `t` (time), `mag` (magnitude), and `err` (magnitude
uncertainty) — are treated specially when present, but **all three are
optional**. Any combination of columns is accepted; VARTOOLS reads whatever
is in the file and uses the column mapping supplied via `-inputlcformat`
(which pyvartools constructs automatically) to identify each vector. An
optional `name` attribute provides a string label that VARTOOLS uses as the
light curve identifier in its output table.

---

## Construction methods

### `LightCurve.from_file(path, format=None, t_col='BJD', mag_col='Mag', err_col='Err', hdu=1, name='')`

Load a light curve from disk. The file format is auto-detected from the file
extension:

- Files whose names end in `.fits`, `.fit`, or `.fts` (case-insensitive) are
  read as FITS binary tables using `astropy.io.fits`. The keyword arguments
  `t_col`, `mag_col`, `err_col`, and `hdu` select which FITS columns map to
  `t`, `mag`, and `err`, and which HDU to use (default: HDU 1, the first
  extension). All remaining table columns are preserved as auxiliary columns.
- All other files are treated as whitespace-delimited ASCII. If the file has
  three or more columns the first three are named `t`, `mag`, `err` and any
  further columns are named `col4`, `col5`, …  If the file has fewer than
  three columns the columns are named `col1`, `col2`, … (since the semantic
  meaning cannot be inferred). All columns are passed to vartools
  automatically via `-inputlcformat` when a `Pipeline` run method is called.

`name` defaults to the stem of `path` (filename without directory or
extension) if not supplied.

!!! note "Column layout for disk-based pipeline runs"
    `from_file` always reads columns 1–3 as `t`, `mag`, `err` for ASCII files.
    If your ASCII file has a different layout (e.g. time in column 3), load the
    file manually or use `from_arrays`, then pass the correct mapping to
    `Pipeline.run_file()` via its `columns` parameter rather than here.

```python
# ASCII, columns 1-2-3
lc = vt.LightCurve.from_file("EXAMPLES/2")

# FITS binary table with default column names (BJD, Mag, Err)
lc = vt.LightCurve.from_file("EXAMPLES/example.fits")

# FITS with non-default column names (here: time, mag, err instead of BJD, Mag, Err)
lc = vt.LightCurve.from_file(
    "EXAMPLES/2.fits", t_col="time", mag_col="mag", err_col="err"
)
```

---

### `LightCurve.from_arrays(t, mag, err, aux=None, name='')`

Construct a `LightCurve` directly from NumPy arrays (or anything that converts
to a 1-D NumPy array). All three arrays must have the same length.

`t`, `mag`, and `err` are all optional — pass `None` (or omit them) when
VARTOOLS will generate or ignore those vectors internally. `aux` is an
optional dict for any additional named columns. All columns present in the
DataFrame are written to the temporary input file; pyvartools constructs a
`-inputlcformat` flag automatically so VARTOOLS can identify each one by name.

```python
import numpy as np

t       = np.linspace(0, 30, 300)
mag     = 10.0 + 0.1 * np.sin(2 * np.pi * t / 2.3)
err     = np.full(300, 0.01)
airmass = 1.0 + 0.5 * np.abs(np.sin(np.pi * t / 15.0))

# Standard three columns plus an extra
lc = vt.LightCurve.from_arrays(t, mag, err, aux={"airmass": airmass}, name="my_star")
# → -inputlcformat t:1,mag:2,err:3,airmass:4 (built automatically)

# Only t and mag — no uncertainty column
lc2 = vt.LightCurve.from_arrays(t=t, mag=mag)
# → -inputlcformat t:1,mag:2

# Only auxiliary data — no standard columns at all
phase_arr = (t % 2.3) / 2.3          # phase-fold at 2.3 d period
flux_arr  = 1.0 - 0.01 * np.sin(2 * np.pi * phase_arr)
lc3 = vt.LightCurve.from_arrays(aux={"phase": phase_arr, "flux": flux_arr})
# → -inputlcformat phase:1,flux:2
```

---

### `LightCurve.from_dataframe(df, name='')`

Wrap an existing `pd.DataFrame`. Any DataFrame is accepted — columns named `t`,
`mag`, and `err` are treated as the standard VARTOOLS vectors when present;
all other columns are preserved as auxiliary variables. None of the three
standard columns are required.

```python
import pandas as pd

df = pd.read_csv("EXAMPLES/example.csv")   # columns: t, mag, err
lc = vt.LightCurve.from_dataframe(df, name="example")
```

---

### `LightCurve.from_timeseries(ts, mag_col='mag', err_col='err', name='')`

Build a `LightCurve` from an `astropy.timeseries.TimeSeries` object. The time
column is extracted from `ts.time` and converted to numeric days (BJD or the
native scale of the `Time` object). `mag_col` and `err_col` select the
magnitude and uncertainty columns by name.

```python
import numpy as np
from astropy.timeseries import TimeSeries
from astropy.time import Time
import astropy.units as u

# Build a TimeSeries from arrays (e.g. after reading a mission file with astropy)
times = Time(2459000 + np.linspace(0, 30, 300), format="jd", scale="tdb")
ts = TimeSeries(time=times)
ts["sap_flux"]     = (10.5 + 0.08 * np.sin(2 * np.pi * np.linspace(0, 30, 300) / 2.3)) * u.mag
ts["sap_flux_err"] = np.full(300, 0.005) * u.mag

lc = vt.LightCurve.from_timeseries(ts, mag_col="sap_flux", err_col="sap_flux_err")
```

---

## Properties

| Property | Type | Description |
|---|---|---|
| `.t` | `np.ndarray` or `None` | Time values (days), or `None` if the column is absent. |
| `.mag` | `np.ndarray` or `None` | Magnitude (or flux) values, or `None` if absent. |
| `.err` | `np.ndarray` or `None` | Per-point magnitude uncertainties, or `None` if absent. |
| `.name` | `str` | String label for this light curve. |
| `.scalars` | `dict[str, float]` | Per-star scalar variables (see below). |

All three array properties return `None` when the corresponding column is not
present in the underlying DataFrame. To modify data, create a new `LightCurve`
from the modified arrays.

### `.scalars`

The canonical per-star scalar store.  A dictionary of scalar variables
associated with this light curve, populated automatically by pyvartools
when a run captures its output LC (`capture_lc=True`) and used as input
to the next segment in a chain.

`Result.lcscalars` and `BatchResult.lcscalars` are shortcuts for reading from
this dict — the values live here, not in the `Result`.  Writing to
`result.lcscalars` does nothing persistent; mutate `result.lc.scalars` (or
build a new `LightCurve` with a different `scalars=` kwarg) to change
what flows into subsequent chain segments.

Its primary role is to carry state across chained pyvartools commands:
when a `Result.LS(...).expr(...)` chain runs, the prior segment's
output columns (and any user-defined scalars) are attached to the next
segment's input LightCurve via this dict, so downstream analytic
expressions can reference them by name.

Keys are the raw vartools variable names (for OUTCOLUMN values carried
forward this is the name *with* the `_N` suffix, e.g.
`"LS_Period_1_0"`; for user-defined scalars it is the bare name, e.g.
`"myvar"`). Values are Python `float` / `int` scalars.

```python
# Typical flow — users rarely touch .scalars directly; it's populated
# by the chain machinery and harvested from -printallscalars output.
r1 = lc.LS(0.5, 10.0, 0.1)              # r1.lc.scalars is empty on first run
r2 = r1.expr("doubled=2*LS_Period_1_0", vartype="scalar")
# Inside the .expr() call, pyvartools attaches r1's output values to the
# input LightCurve's .scalars and injects them into the new vartools run
# via `-expr const`.  See the Fluent API docs for details.

r2.lc.scalars["doubled"]       # == 2 * r1.vars["LS_Period_1_0"]
```

You can also construct a LightCurve with an explicit `scalars` dict — useful
for manual chaining or for seeding values in tests:

```python
lc = vt.LightCurve.from_arrays(
    t=t, mag=mag, err=err,
    scalars={"P0": 1.234, "offset": 0.05},
)
# Subsequent commands referencing "P0" or "offset" in expressions will
# resolve to these values.
```

---

## Conversion methods

### `.to_dataframe()` → `pd.DataFrame`

Return the full underlying DataFrame, including any auxiliary columns. The
index is a default integer index; the three core columns are always named `t`,
`mag`, and `err`.

```python
df = lc.to_dataframe()
print(df.head())
```

---

### `.to_arrays()` → `(t, mag, err)`

Return the three core arrays as a tuple of NumPy arrays. Equivalent to
`(lc.t, lc.mag, lc.err)` but convenient for unpacking.

```python
t, mag, err = lc.to_arrays()
```

---

### `.to_timeseries()` → `astropy.timeseries.TimeSeries`

Convert to an astropy `TimeSeries`. The `t` column is converted to an
`astropy.time.Time` object (scale `'tdb'`, format `'jd'`). Magnitude and
uncertainty are stored as columns named `mag` and `err`. Auxiliary columns are
included as additional columns.

```python
import tempfile

ts = lc.to_timeseries()
with tempfile.NamedTemporaryFile(suffix=".ecsv", delete=False) as f:
    ts.write(f.name, format="ascii.ecsv", overwrite=True)
```

---

## Special methods

`len(lc)` returns the number of data points. `repr(lc)` / `str(lc)` produces a
summary string:

```
LightCurve(name='my_star', n=300, cols=['t', 'mag', 'err'])
```

---

## Full example

```python
import numpy as np
import pyvartools as vt

# From file (auto-detects ASCII)
lc = vt.LightCurve.from_file("EXAMPLES/2")

# From FITS (default column names BJD, Mag, Err — no t_col/mag_col/err_col needed)
lc = vt.LightCurve.from_file("EXAMPLES/example.fits")

# From arrays
t = np.linspace(0, 30, 300)
mag = 10.0 + 0.1*np.sin(2*np.pi*t/2.3)
err = np.full(300, 0.01)
lc = vt.LightCurve.from_arrays(t, mag, err, name="my_star")

print(lc)       # LightCurve(name='my_star', n=300, cols=['t', 'mag', 'err'])
print(len(lc))  # 300

# Round-trip through a DataFrame
df = lc.to_dataframe()
lc2 = vt.LightCurve.from_dataframe(df, name="my_star")
```
