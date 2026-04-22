# Results

When a pipeline is executed, pyvartools returns either a `Result` (single
light curve) or a `BatchResult` (multiple light curves). Both objects give
structured access to the statistics table produced by VARTOOLS, to the
processed light curve data, and to any auxiliary output files that the commands
wrote to disk (periodograms, model files, etc.).

---

## Result

Returned by `pipe.run(lc)`, `pipe.run_file(path)`, `lc.CMD(...)`,
`vt.CMD(lc_input, ...)`, and similar single-LC execution methods.

### `.vars` — `pd.Series`

A pandas Series indexed by VARTOOLS column names. Each element corresponds to
one column in the VARTOOLS output table for the single light curve that was
processed.

```python
result = pipe.run(lc)

# Access by column name
period  = float(result.vars["LS_Period_1_0"])
log_fap = float(result.vars["Log10_LS_Prob_1_0"])
snr     = float(result.vars["LS_SNR_1_0"])

print(f"Best period: {period:.5f} d  (log FAP = {log_fap:.1f})")
```

Column names follow the convention `CommandStat_peak_commandindex`; see
[Output column naming](index.md#output-column-naming) for details.

### `.varobjs` — `VarsNamespace`

Structured, per-command access to output variables. Each command's results
are grouped under a namespace attribute named after the command.

```python
result = vt.Pipeline([cmd.LS(0.5, 10.0, 1e-3)]).run(lc)

# Access the top LS period directly
period  = result.varobjs.LS.Period_1        # e.g. 1.23534
log_fap = result.varobjs.LS.Log10_Prob_1   # e.g. -4222.3
```

When the same command appears more than once in the pipeline, index into
the list explicitly:

```python
result = vt.Pipeline([
    cmd.LS(0.5, 5.0, 1e-3),   # index 0
    cmd.LS(5.0, 50.0, 0.1),   # index 1
]).run(lc)

p_short = result.varobjs.LS[0].Period_1   # from first LS call
p_long  = result.varobjs.LS[1].Period_1   # from second LS call
```

For a single call, both `result.varobjs.LS.Period_1` and
`result.varobjs.LS[0].Period_1` work.

### Attribute shorthand

Any VARTOOLS output key can also be accessed directly as an attribute on
`result`:

```python
period  = result.LS_Period_1_0        # equivalent to result.vars["LS_Period_1_0"]  (attribute shorthand)
log_fap = result.Log10_LS_Prob_1_0
```

This is convenient for interactive use. A clean `AttributeError` is raised
for keys that do not exist.

### `.lc` — `LightCurve` or `None`

The light curve as it exists at the end of the pipeline — after all filtering,
detrending, and model corrections. This attribute is `None` unless
`capture_lc=True` was passed to the run method.

```python
result = pipe.run(lc, capture_lc=True)
if result.lc is not None:
    t, mag, err = result.lc.to_arrays()
```

### `.files` — `dict[str, pd.DataFrame]`

A dictionary mapping a logical file name to a `pd.DataFrame` containing the
contents of that output file. Keys follow the pattern `"{CommandName}_{logical}_{idx}"`
where `idx` is the zero-based position of the command in the pipeline
(e.g., `"LS_periodogram_0"`, `"BLS_model_1"`). This dict is empty when no commands write output files.

```python
from pyvartools import commands as cmd

pipe = vt.Pipeline([
    cmd.LS(0.5, 10.0, 0.01, save_periodogram=True),
])
result = pipe.run(lc)

pgram = result.files["LS_periodogram_0"]   # pd.DataFrame with frequency/power cols
pgram.plot(x="Frequency", y="Power")
```

### Per-star scalars — `result.lc.scalars`

Per-star scalars (vectortype `SCALAR`, `PERSTARDATA`, or `INLIST`) are
stored directly on the captured light curve at
[`LightCurve.scalars`](lightcurve.md#scalars).  These come from `-expr
scalar ...`, `-expr listvar ...`, or `-inlistvars`.  Keys are the raw
variable names with no `_N` suffix, in contrast to `result.vars`, which
holds OUTCOLUMN values whose names carry a `_N` suffix reflecting the
command's position (e.g. `"LS_Period_1_0"`).

pyvartools enables the underlying [`-printallscalars`](../cli/options.md#-printallscalars)
option automatically when running chained commands, so user-defined
scalars round-trip into `result.lc.scalars` with no extra configuration.

```python
# -expr scalar creates a SCALAR variable — it lives in lc.scalars, not .vars
r = lc.LS(0.5, 10.0, 0.1).expr("doubled=2*LS_Period_1_0", vartype="scalar")
r.vars["LS_Period_1_0"]     # 1.2344 — OUTCOLUMN, lives in .vars
r.lc.scalars["doubled"]     # 2.4688 — SCALAR on the captured LightCurve
```

If `result.lc is None` (no LC captured) there are no scalars to inspect —
check `result.lc` before reading `.scalars`.

### `.ok` and `.error`

`.ok` is `True` when the run completed without error. `.error` is `None` on
success or a `RunError` instance if the run failed.

---

## BatchResult

Returned by `pipe.run_batch(lcs)`, `pipe.run_filelist(paths)`, and
`LightCurveBatch.run()`.

### `.vars` — `pd.DataFrame`

A pandas DataFrame with one row per light curve. The `Name` column contains
the light curve identifier; remaining columns are the VARTOOLS statistics.
Column names are identical to those produced by the CLI with `-header`.

```python
batch = pipe.run_batch(lcs)

# Print the best period and log FAP for each light curve
print(batch.vars[["Name", "LS_Period_1_0", "Log10_LS_Prob_1_0"]])

# Find the light curve with the strongest detection
best = batch.vars.loc[batch.vars["Log10_LS_Prob_1_0"].idxmin()]
print(f"Most significant: {best['Name']}  P = {float(best['LS_Period_1_0']):.4f} d")
```

### Per-LC access

Individual `Result` objects can be retrieved by index or by iteration.
Slices return a new `BatchResult` containing only the selected light curves:

```python
batch = pipe.run_batch(lcs)

# Index access
r0 = batch[0]            # Result for the first LC
print(r0.vars["LS_Period_1_0"])
print(r0.varobjs.LS.Period_1)

# Slice access — returns a sub-BatchResult
sub = batch[1:4]         # LCs 1, 2, 3
every_other = batch[::2]

# Iteration
for result in batch:
    if result.ok:
        print(result.varobjs.LS.Period_1)
    else:
        print("failed:", result.error)

# Length
print(len(batch))        # number of light curves
```

For runs produced by `LightCurveBatch`, per-LC errors are also captured:

```python
batch_result = vt.LightCurveBatch(lcs).LS(0.5, 10.0, 1e-3).run()
for i, r in enumerate(batch_result):
    if not r.ok:
        print(f"LC {i} failed: {r.error}")
```

### `.lcs` — `list[LightCurve]`

List of processed light curves in the same order as the input.  Returns an
empty list when `capture_lc=False` was used (the default for
`Pipeline.run_batch()`), so ``for lc in batch.lcs:`` is always safe.
`LightCurveBatch.run()` defaults to `capture_lc=True`.

```python
batch = pipe.run_batch(lcs, capture_lc=True)
for source_lc in batch.lcs:
    print(source_lc.name, len(source_lc))
```

### `.files` — `dict[str, list[pd.DataFrame]]`

Maps each logical output file name to a list of DataFrames, one per input
light curve (same order as the input). Files that were not produced for a
particular light curve appear as `None` at that position.

```python
for name, df in zip([lc.name for lc in lcs], batch.files["LS_periodogram_0"]):
    if df is not None:
        df.to_csv(f"{name}_pgram.csv", index=False)
```

### `.lcscalars` — `pd.DataFrame`

A convenience view of every captured LC's
[`LightCurve.scalars`](lightcurve.md#scalars), packaged as a DataFrame
with one row per LC (input order) and one column per scalar variable
name.  Equivalent to building `pd.DataFrame([lc.scalars for lc in
batch.lcs])` — the scalars themselves live on the `LightCurve`
objects in `batch.lcs`, not in a separate store.

Column names are the raw variable names (no `_N` suffix); column values
may differ per LC (e.g. for a `listvar` that depends on the LC data).
The batched counterpart of accessing `result.lc.scalars` on a single-LC
`Result`.

```python
br = vt.LightCurveBatch(lcs).LS(0.5, 10.0, 0.1).run()
br2 = br.expr("doubled=2*LS_Period_1_0", vartype="scalar").run()
br2.lcscalars         # pd.DataFrame: one row per LC, column 'doubled' = 2 * LS_Period_1_0

# Same data via the canonical path:
br2.lcs[0].scalars["doubled"]   # same value as br2.lcscalars.iloc[0]["doubled"]
```

Returns an empty DataFrame when no LCs have scalars (e.g. `capture_lc=False`
so the output LCs weren't captured, or the commands produced no scalar
variables).

### `.ok` — `bool`

`True` when the overall batch run completed without error.

### `.error` — `RunError` or `None`

Set to a `RunError` instance when the entire batch run failed; `None` on
success.  For per-LC errors from `LightCurveBatch.run()`, check
`batch[i].error` instead.

---

## RunError

`pyvartools.RunError` is raised (and also stored in `BatchResult.error`) when
the `vartools` subprocess exits with a non-zero status. Its message includes
the reconstructed command string and the full stderr output from the binary,
making it straightforward to diagnose input format problems or invalid
parameter values.

```python
from pyvartools import RunError

try:
    result = pipe.run(lc)
except RunError as exc:
    print("vartools failed:")
    print(exc)          # prints command + stderr
```

For batch runs the exception is not raised by default; check `.ok` and
`.error` instead:

```python
batch = pipe.run_batch(lcs)
if not batch.ok:
    print(batch.error)
```

To raise on batch failure, pass `raise_on_error=True`:

```python
batch = pipe.run_batch(lcs, raise_on_error=True)
```

---

## Full example

```python
import pyvartools as vt
from pyvartools import commands as cmd

# Build pipeline
pipe = vt.Pipeline([
    cmd.LS(0.5, 10.0, 1e-3),
    cmd.Killharm(period="ls", nharm=2),
])

# --- Single LC ---
lc = vt.LightCurve.from_file("EXAMPLES/2")
result = pipe.run(lc, capture_lc=True)

# Three equivalent ways to read the LS period:
print(result.vars["LS_Period_1_0"])       # Series key access
print(result.varobjs.LS.Period_1)           # structured namespace
print(result.LS_Period_1_0)              # attribute shorthand

print(result.lc)                         # post-correction LightCurve
print(result.files.keys())               # auxiliary output files

# --- Batch ---
paths = [f"EXAMPLES/{i}" for i in range(1, 10)]
batch = pipe.run_filelist(paths)

if not batch.ok:
    raise RuntimeError(str(batch.error))

# Best-period summary
summary = batch.vars[["Name", "LS_Period_1_0", "Log10_LS_Prob_1_0"]]
summary.to_csv("ls_results.csv", index=False)
```
