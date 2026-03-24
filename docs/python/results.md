# Results

When a pipeline is executed, pyvartools returns either a `Result` (single
light curve) or a `BatchResult` (multiple light curves). Both objects give
structured access to the statistics table produced by VARTOOLS, to the
processed light curve data, and to any auxiliary output files that the commands
wrote to disk (periodograms, model files, etc.).

---

## Result

Returned by `pipe.run(lc)` and `pipe.run_file(path)`.

### `.stats` — `pd.Series`

A pandas Series indexed by VARTOOLS column names. Each element corresponds to
one column in the VARTOOLS output table for the single light curve that was
processed.

```python
result = pipe.run(lc)

# Access by column name
period  = float(result.stats["LS_Period_1_0"])
log_fap = float(result.stats["Log10_LS_Prob_1_0"])
snr     = float(result.stats["LS_SNR_1_0"])

print(f"Best period: {period:.5f} d  (log FAP = {log_fap:.1f})")
```

Column names follow the convention `CommandStat_peak_commandindex`; see
[Output column naming](index.md#output-column-naming) for details.

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
    cmd.LS(minp=0.5, maxp=10.0, stepsize=0.01, Nharm=1, Nexp=0,
           oper="periodogram"),
])
result = pipe.run(lc)

pgram = result.files["LS_periodogram_0"]   # pd.DataFrame with frequency/power cols
pgram.plot(x="Frequency", y="Power")
```

---

## BatchResult

Returned by `pipe.run_batch(lcs)` and `pipe.run_filelist(paths)`.

### `.stats` — `pd.DataFrame`

A pandas DataFrame with one row per light curve. The `Name` column contains
the light curve identifier; remaining columns are the VARTOOLS statistics.
Column names are identical to those produced by the CLI with `-header`.

```python
batch = pipe.run_batch(lcs)

# Print the best period and log FAP for each light curve
print(batch.stats[["Name", "LS_Period_1_0", "Log10_LS_Prob_1_0"]])

# Find the light curve with the strongest detection
best = batch.stats.loc[batch.stats["Log10_LS_Prob_1_0"].idxmin()]
print(f"Most significant: {best['Name']}  P = {float(best['LS_Period_1_0']):.4f} d")
```

### `.lcs` — `list[LightCurve]` or `None`

List of processed light curves in the same order as the input, or `None` if
`capture_lc=False` (the default).

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

### `.ok` — `bool`

`True` when the vartools subprocess exited with status 0. `False` if the run
failed.

### `.error` — `RunError` or `None`

Set to a `RunError` instance when the run failed; `None` on success.

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
    cmd.LS(minp=0.5, maxp=10.0, stepsize=0.01, Nharm=1, Nexp=0,
           oper="ls_periodogram"),
    cmd.Killharm(method="ls", nharmonics=0, nsubharmonics=0, npeaks=1,
                 omodel="killharm_model"),
])

# --- Single LC ---
lc = vt.LightCurve.from_file("EXAMPLES/2")
result = pipe.run(lc, capture_lc=True)

print(result.stats["LS_Period_1_0"])           # best LS period
print(result.stats["Killharm_Amplitude_1_1"])  # harmonic fit amplitude
print(result.lc)                               # post-correction LightCurve
print(result.files.keys())                     # {'ls_periodogram', 'killharm_model'}

# --- Batch ---
paths = [f"EXAMPLES/{i}" for i in range(1, 10)]
batch = pipe.run_filelist(paths)

if not batch.ok:
    raise RuntimeError(str(batch.error))

# Best-period summary
summary = batch.stats[["Name", "LS_Period_1_0", "Log10_LS_Prob_1_0"]]
summary.to_csv("ls_results.csv", index=False)
```
