# Command Wrappers

Every vartools command has a Python wrapper class in the `pyvartools.commands` module (imported as `cmd` by convention). All wrappers are subclasses of `VartoolsCommand` and implement `_to_cli_args()`, which translates the Python constructor arguments into the corresponding vartools command-line tokens.

```python
from pyvartools import commands as cmd
```

You can call `cmd.LS.help()` on any class to print the full vartools help text for that command, pulled live from the installed binary.

---

## Auxiliary output files

Many commands can write auxiliary files (periodograms, model curves, coefficient tables, etc.) to disk. Every such parameter is named `save_*` and accepts four modes, controlled by passing a `bool`, a path string, or an `Output` object.

### The four modes

| Value | Mode | Written to disk? | Captured in Python? |
|-------|------|-----------------|---------------------|
| `False` (default) | 4 — suppress | no | no |
| `True` | 1 — temp, capture | temp dir (auto-deleted) | **yes** — `result.files["key"]` |
| `"/path/to/dir"` | 3 — disk only | that directory | no |
| `Output("/path/to/dir", capture=True)` | 2 — disk + capture | that directory | **yes** — `result.files["key"]` |

### `Output` class

```python
from pyvartools import Output

Output(path=None, capture=True)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `path` | `str` or `None` | Directory to write the file to. `None` means use a pipeline-managed temp directory. |
| `capture` | `bool` | Whether to read the file into Python and include it in `result.files`. Default `True`. |

### Examples

```python
import pyvartools as vt
from pyvartools import commands as cmd, Output

lc = vt.LightCurve.from_file("EXAMPLES/2")

# Mode 1 (default True): temp dir, captured into result.files
result = vt.Pipeline([
    cmd.LS(0.1, 10.0, 1e-3, save_periodogram=True),
]).run(lc)
pgram = result.files["LS_periodogram_0"]   # pd.DataFrame

# Mode 3: written to EXAMPLES/OUTDIR1, NOT in result.files
result = vt.Pipeline([
    cmd.LS(0.1, 10.0, 1e-3, save_periodogram="EXAMPLES/OUTDIR1"),
]).run(lc)
# EXAMPLES/OUTDIR1/stdin.ls written; result.files has no "LS_periodogram_0"

# Mode 2: written to EXAMPLES/OUTDIR1 AND captured
result = vt.Pipeline([
    cmd.LS(0.1, 10.0, 1e-3,
           save_periodogram=Output("EXAMPLES/OUTDIR1", capture=True)),
]).run(lc)
pgram = result.files["LS_periodogram_0"]   # captured from EXAMPLES/OUTDIR1/stdin.ls

# Mode 4 (default False): nothing written, nothing captured
result = vt.Pipeline([
    cmd.LS(0.1, 10.0, 1e-3, save_periodogram=False),
]).run(lc)
# result.files has no "LS_periodogram_0"
```

### Output file keys

Captured files appear in `result.files` under a key of the form `"{CommandName}_{logical_name}_{idx}"`, where `idx` is the 0-based index of the command in the pipeline. For example, the first `LS` command's periodogram is `"LS_periodogram_0"`, and the second is `"LS_periodogram_1"`.

In a batch run (`run_batch` / `run_filelist`) each key maps to a list of DataFrames — one per light curve.

### Note on `autocorrelation`

`autocorrelation` is a special case: the vartools CLI always writes the output file regardless of the `save_result` setting. Passing `save_result=False` suppresses Python capture but the file is still written to a temp directory (and discarded after the run).

---

## Period finding

### `LS` — Generalized Lomb-Scargle

```python
cmd.LS(minp, maxp, subsample, npeaks=5, save_periodogram=False,
       noGLS=False, whiten=False, clip=None, clipiter=None,
       bootstrap=None, maskpoints=None, fixperiod_snr=None)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `minp`, `maxp` | `float`, `str`, numpy array, `PerLC`, or `pd.Series` | Period search range (same units as the time column, typically days). See [Variable and expression parameters](#variable-and-expression-parameters) below; for batch per-LC values see [Per-LC array parameters](pipeline.md#per-lc-array-parameters). |
| `subsample` | `float`, `str`, numpy array, `PerLC`, or `pd.Series` | Frequency step as a fraction of 1/T (time span). Smaller values = finer resolution. Typical: `1e-3`. Accepts variable names, expressions, and per-LC arrays. |
| `npeaks` | `int` | Number of peaks to report. Default `5`. |
| `save_periodogram` | `bool`, `str`, or `Output` | Controls auxiliary file output. `True` captures as `result.files["LS_periodogram_N"]`; a path string writes to that directory without capturing; `Output(path, capture=True)` does both. See [Auxiliary output files](#auxiliary-output-files). |
| `noGLS` | `bool` | Use the classical (non-generalised) Lomb-Scargle statistic instead of the generalised (Zechmeister & Kurster) version. |
| `whiten` | `bool` | Iteratively subtract the best-fit sinusoid and re-search. |
| `clip` | `float` | Sigma-clipping threshold applied during whitening iterations. |
| `bootstrap` | `int` | Number of bootstrap resamples for false-alarm probability estimation. |
| `maskpoints` | `str` or `None` | Name of a mask variable; points where the variable is non-zero are excluded from the periodogram. |
| `fixperiod_snr` | `float`, `int`, `str`, or `None` | Evaluate the periodogram at a known period and report its significance. See [fixperiod_snr — fixed-period significance](#fixperiod_snr--fixed-period-significance) below. |

Output statistics include `LS_Period_1_N`, `LS_Amplitude_1_N`, `Log10_LS_Prob_1_N` (log₁₀ of the false-alarm probability) for each of the `npeaks` peaks (N = 0-based command index, or the `columnsuffix` string).

CLI equivalent: `-LS minp maxp subsample npeaks [operiodogram dir] [noGLS] [whiten] ...`

#### fixperiod_snr — fixed-period significance

When `fixperiod_snr` is set, vartools evaluates the periodogram at a single known period and appends four extra output columns (N = 0-based pipeline index of the LS command):

| Column | Description |
|--------|-------------|
| `LS_PeriodFix_N` | The fixed period used (omitted when `fixperiod_snr` is a plain number, since it is already known). |
| `Log10_LS_Prob_PeriodFix_N` | Log₁₀ false-alarm probability at that period. |
| `LS_Periodogram_Value_PeriodFix_N` | LS power at that period. |
| `LS_SNR_PeriodFix_N` | SNR = (power − mean) / σ. |

The accepted forms and the CLI tokens they emit:

| Python value | Emitted CLI tokens | When to use |
|---|---|---|
| `1.234` (number) | `fixperiodSNR fix 1.234` | Period known at pipeline-construction time. |
| `"ls"` | `fixperiodSNR ls` | Use the best period found by the most recent prior `-LS` run. |
| `"aov"` | `fixperiodSNR aov` | Use the best period found by the most recent prior `-aov` run. |
| `"injectharm"` | `fixperiodSNR injectharm` | Use the injected-signal period from a prior injection run. |
| `"fixcolumn LS_Period_1_0"` | `fixperiodSNR fixcolumn LS_Period_1_0` | Read the period from a named per-star column. |
| `"list"` | `fixperiodSNR list` | Read the period from the current list-file column. |
| `"list column 2"` | `fixperiodSNR list column 2` | Read the period from column 2 of the list file. |

#### Variable and expression parameters

`minp`, `maxp`, and `subsample` each accept four forms:

| Value | Emitted CLI tokens | When to use |
|-------|--------------------|-------------|
| A number (`float` or `int`) | `0.5` | Fixed value known at pipeline-construction time. |
| A bare identifier string, e.g. `"minperiod"` | `var minperiod` | Value is read from a named per-star vartools variable — typically one loaded from a list file via `run_filelist`. |
| Any other string, e.g. `"tspan/200"` | `expr tspan/200` | Evaluated as a math expression using vartools' built-in expression engine, per light curve. |
| A numpy array, `PerLC`, or `pd.Series` | *(handled automatically)* | A different value for each light curve in a batch run. See [Per-LC array parameters](pipeline.md#per-lc-array-parameters). |

The identifier rule is: if the string matches `[A-Za-z_]\w*` it is treated as a variable name; otherwise it is treated as an expression.

!!! note "Defining variables for the `expr` form"
    The `expr` keyword evaluates an expression against vartools' internal variable registry at the time each light curve is processed. Variables such as `tspan` are *not* built-in; they must be defined by prior commands in the same pipeline. Use `cmd.stats` to compute per-star statistics and `cmd.expr` to derive new variables from them:

    ```python
    cmd.stats("t", "min,max")                         # → STATS_t_MIN_0, STATS_t_MAX_0
    cmd.expr("tspan=STATS_t_MAX_0-STATS_t_MIN_0")     # → tspan
    ```

    The `var` form similarly requires the named variable to exist in the per-star variable registry. This is most naturally supplied via `run_filelist` with a list file that includes per-star columns for `minp` and `maxp`.

**Examples**

```python
import pyvartools as vt
from pyvartools import commands as cmd

lc = vt.LightCurve.from_file("EXAMPLES/2")

# Fixed values — search periods 0.1–10 days, report top 5 peaks with whitening
result = vt.Pipeline([
    cmd.LS(0.1, 10.0, 0.1, npeaks=5, whiten=True, clip=5.0, clipiter=1,
           save_periodogram=True),
]).run(lc)
print(result.stats["LS_Period_1_0"])       # 1.23440877
print(result.stats["Log10_LS_Prob_1_0"])   # -4000.59209
pgram = result.files["LS_periodogram_0"]   # pd.DataFrame: frequency vs power

# Expression form — set period range relative to the time baseline of each LC.
# First compute min/max time with cmd.stats, then define tspan with cmd.expr,
# then pass expressions to LS.  LS is at pipeline index 2, so keys end in "_2".
result = vt.Pipeline([
    cmd.stats("t", "min,max"),                              # → STATS_t_MIN_0, STATS_t_MAX_0
    cmd.expr("tspan=STATS_t_MAX_0-STATS_t_MIN_0"),          # → tspan
    cmd.LS("tspan/200", "tspan/2", 1e-3, npeaks=1),         # expr tspan/200, expr tspan/2
]).run(lc)
print(result.stats["LS_Period_1_2"])       # 1.23534018

# Variable form — minp and maxp are per-star variables read from a list file.
# Each row in the list file supplies different search bounds for each LC.
# batch = pipe.run_filelist("lc_list.txt")   # list file has minp and maxp columns

# Batch: run on many light curves in parallel
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 11)]
batch = vt.Pipeline([cmd.LS(0.1, 10.0, 0.1, npeaks=1)]).run_batch(lcs, nthreads=4)
print(batch.stats[["Name", "LS_Period_1_0", "Log10_LS_Prob_1_0"]])

# fixperiod_snr — evaluate LS at a known period
lc = vt.LightCurve.from_file("EXAMPLES/2")

# Fixed number form: evaluate at period = 1.234
r = vt.Pipeline([cmd.LS(0.1, 10.0, 0.1, fixperiod_snr=1.234)]).run(lc)
print(r.stats["LS_SNR_PeriodFix_0"])            # SNR at period 1.234

# "ls" form: evaluate at the best peak from a prior LS search
r = vt.Pipeline([
    cmd.LS(0.1, 10.0, 0.1, npeaks=1),
    cmd.LS(0.1, 10.0, 0.1, fixperiod_snr="ls"),
]).run(lc)
print(r.stats["LS_SNR_PeriodFix_1"])            # SNR at period found by first LS

# "aov" form: evaluate LS at the best period from a prior AOV search
r = vt.Pipeline([
    cmd.aov(0.1, 10.0, 0.1, 0.01, npeaks=1),
    cmd.LS(0.1, 10.0, 0.1, fixperiod_snr="aov"),
]).run(lc)
print(r.stats["LS_SNR_PeriodFix_1"])

# "fixcolumn" form: read the period from a named per-star column
r = vt.Pipeline([
    cmd.LS(0.1, 10.0, 0.1, npeaks=1),
    cmd.LS(0.1, 10.0, 0.1, fixperiod_snr="fixcolumn LS_Period_1_0"),
]).run(lc)
print(r.stats["LS_PeriodFix_1"])
print(r.stats["Log10_LS_Prob_PeriodFix_1"])
```

---

### `aov` — Analysis of Variance

```python
cmd.aov(minp, maxp, subsample, finetune, npeaks=5, nbin=None,
        save_periodogram=False, whiten=False, clip=None,
        clipiter=None, uselog=False, maskpoints=None, fixperiod_snr=None)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `minp`, `maxp`, `subsample` | `float` or `str` | Period range and frequency step (same as LS). Accept variable names and expressions — see [Variable and expression parameters](#variable-and-expression-parameters).  Also accept per-LC arrays — see [Per-LC array parameters](pipeline.md#per-lc-array-parameters). |
| `finetune` | `float` or `str` | Oversampling factor applied near peak frequencies for fine-tuning. Accepts var/expr forms and per-LC arrays. |
| `npeaks` | `int` | Number of peaks to report. |
| `nbin` | `int`, `str`, or `None` | Number of phase bins. If `None`, vartools selects automatically. Accepts var/expr forms and per-LC arrays. |
| `save_periodogram` | `bool`, `str`, or `Output` | Auxiliary file output. `True` captures as `result.files["aov_periodogram_N"]`. See [Auxiliary output files](#auxiliary-output-files). |
| `uselog` | `bool` | Use the log of the AOV statistic. |
| `fixperiod_snr` | `float`, `int`, `str`, or `None` | Evaluate the AOV periodogram at a known period and report its significance. Same forms as `LS.fixperiod_snr`. When set, appends `aov_PeriodFix_N`, `Log10_aov_Prob_PeriodFix_N`, `aov_Periodogram_Value_PeriodFix_N`, `aov_SNR_PeriodFix_N`. |

Use AoV instead of LS when the signal is strictly periodic but non-sinusoidal (e.g. eclipsing binaries, pulsating stars). AoV is less sensitive to the shape of the variation.

CLI equivalent: `-aov [Nbin N] minp maxp subsample finetune npeaks ...`

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/2")

result = vt.Pipeline([
    cmd.aov(0.1, 10.0, 0.1, 0.01, npeaks=5, nbin=20,
            whiten=True, clip=5.0, clipiter=1, save_periodogram=True),
]).run(lc)
print(result.stats["AOV_Period_1_0"])   # 1.23583047
pgram = result.files["aov_periodogram_0"]   # pd.DataFrame: frequency vs AOV statistic

# fixperiod_snr — evaluate AOV at a known period
result = vt.Pipeline([
    cmd.aov(0.1, 10.0, 0.1, 0.01, npeaks=1),
    cmd.aov(0.1, 10.0, 0.1, 0.01, fixperiod_snr="aov"),
]).run(lc)
print(result.stats["aov_SNR_PeriodFix_1"])
```

---

### `aov_harm` — Multi-harmonic AoV

```python
cmd.aov_harm(nharm, minp, maxp, subsample, finetune, npeaks=5,
             save_periodogram=False, whiten=False, clip=None,
             clipiter=None, maskpoints=None, fixperiod_snr=None)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `nharm` | `int` or `str` | Number of harmonics to include in the model.  Also accepts variable names (`"nharmvar"`) and expressions (`"npeaks*2"`) — see [Variable and expression parameters](#variable-and-expression-parameters).  Also accepts per-LC arrays — see [Per-LC array parameters](pipeline.md#per-lc-array-parameters). |
| `minp`, `maxp`, `subsample`, `finetune` | `float` or `str` | Same as `aov`. Accept variable names, expressions, and per-LC arrays — see [Variable and expression parameters](#variable-and-expression-parameters) and [Per-LC array parameters](pipeline.md#per-lc-array-parameters). |
| `fixperiod_snr` | `float`, `int`, `str`, or `None` | Evaluate the AOV_HARM periodogram at a known period and report its significance. Same forms and output columns as `aov.fixperiod_snr`. |
| All others | — | Same as `aov`. |

Multi-harmonic AoV projects the phase-folded light curve onto a truncated Fourier series. It is better than plain AoV for highly non-sinusoidal signals such as RR Lyrae, Cepheids, and W UMa systems.

CLI equivalent: `-aov_harm nharm minp maxp subsample finetune npeaks ...`

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/2")

result = vt.Pipeline([
    cmd.aov_harm(1, 0.1, 10.0, 0.1, 0.01, npeaks=2,
                 whiten=True, clip=5.0, clipiter=1, save_periodogram=True),
]).run(lc)
print(result.stats["AOV_HARM_Period_1_0"])   # 1.23533969
pgram = result.files["aov_harm_periodogram_0"]

# fixperiod_snr — evaluate AOV_HARM at the period found by a prior LS search
result = vt.Pipeline([
    cmd.LS(0.1, 10.0, 0.1, npeaks=1),
    cmd.aov_harm(2, 0.1, 10.0, 0.1, 0.01, fixperiod_snr="ls"),
]).run(lc)
print(result.stats["aov_SNR_PeriodFix_1"])
```

---

### `BLS` — Box-Least-Squares transit search

```python
cmd.BLS(minper, maxper, rmin=0.01, rmax=0.1, nbins=200,
        timezone=0, npeaks=1, subsample=1.0, nfreq=None,
        qmin=None, qmax=None,
        density_mode=False, stellar_density=None,
        min_exp_dur_frac=0.5, max_exp_dur_frac=1.5,
        df=None, extraparams=False, nobinnedrms=False,
        freq_grid=None, adjust_qmin=False, reduce_nbins=False,
        reportharmonics=False,
        save_periodogram=False, save_model=False,
        save_phcurve=False, save_jdcurve=False,
        ophcurve_phmin=0, ophcurve_phmax=1, ophcurve_phstep=0.005,
        ojdcurve_jdstep=0.02,
        correct_lc=False, fittrap=False, maskpoints=None)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `minper`, `maxper` | `float` or `str` | Period search range.  Also accept variable names and expressions — see [Variable and expression parameters](#variable-and-expression-parameters).  Also accept per-LC arrays — see [Per-LC array parameters](pipeline.md#per-lc-array-parameters). |
| `rmin`, `rmax` | `float` or `str` | Minimum and maximum fractional transit duration (`r` mode, as a fraction of the orbital period).  Ignored when `qmin`/`qmax` are set.  Accept var/expr forms and per-LC arrays. |
| `qmin`, `qmax` | `float`, `str`, or `None` | Minimum and maximum fractional transit duration in `q` mode (`q` = ingress-to-egress fraction, not full duration).  When set, emits `"q" qmin qmax` instead of `"r" rmin rmax`.  Accept var/expr forms and per-LC arrays. |
| `nbins` | `int` or `str` | Number of phase bins.  Accepts var/expr forms and per-LC arrays. |
| `timezone` | `float` | Time-zone offset (0 for HJD/BJD). |
| `npeaks` | `int` | Number of transit candidates to report. |
| `subsample` | `float` or `str` | Frequency oversampling factor.  Accepts var/expr forms and per-LC arrays. |
| `nfreq` | `int`, `str`, or `None` | Fixed number of test frequencies (overrides `subsample`).  Accepts var/expr forms and per-LC arrays. |
| `density_mode` | `bool` | Use stellar density to set transit duration bounds instead of `rmin`/`rmax`. |
| `stellar_density` | `float`, `str`, or `None` | Stellar density (g/cm³) for density mode.  Accepts var/expr forms and per-LC arrays. |
| `min_exp_dur_frac`, `max_exp_dur_frac` | `float` or `str` | Expected-duration fractions for density mode (default `0.5` and `1.5`).  Accept var/expr forms and per-LC arrays. |
| `df` | `float`, `str`, or `None` | Fixed frequency step (alternative to `subsample`).  Accepts var/expr forms and per-LC arrays. |
| `extraparams` | `bool` | Compute extra BLS statistics (epoch, depth, etc.). |
| `nobinnedrms` | `bool` | Do not compute binned RMS statistics. |
| `freq_grid` | `str` or `None` | Frequency grid mode: `"stepP"` (uniform in period) or `"steplogP"` (log-uniform). |
| `adjust_qmin` | `bool` | Automatically adjust `qmin` based on the minimum cadence. |
| `reduce_nbins` | `bool` | Reduce `nbins` when `adjust_qmin` requires it (only active when `adjust_qmin=True`). |
| `reportharmonics` | `bool` | Report harmonic periods (½, ⅓, … of best period) as additional candidates. |
| `save_periodogram` | `bool`, `str`, or `Output` | Auxiliary file output. `True` captures as `result.files["BLS_periodogram_N"]`. See [Auxiliary output files](#auxiliary-output-files). |
| `save_model` | `bool`, `str`, or `Output` | Auxiliary file output. `True` captures as `result.files["BLS_model_N"]`. See [Auxiliary output files](#auxiliary-output-files). |
| `save_phcurve` | `bool`, `str`, or `Output` | Phase-folded model curve. `True` captures as `result.files["BLS_phcurve_N"]`. |
| `save_jdcurve` | `bool`, `str`, or `Output` | JD-sampled model curve. `True` captures as `result.files["BLS_jdcurve_N"]`. |
| `correct_lc` | `bool` | Subtract the best-fit box from the light curve. |
| `fittrap` | `bool` | Fit a trapezoidal rather than a box-shaped transit. |

Key output statistics: `BLS_Period_1_N`, `BLS_SDE_1_N` (signal detection efficiency), `BLS_SN_1_N` (signal-to-noise), `BLS_Depth_1_N`, `BLS_Qtran_1_N` (fractional duration), `BLS_Tc_1_N` (transit epoch).

CLI equivalent: `-BLS r rmin rmax minper maxper [nf N | optimal s] nbins timezone npeaks ...`

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/3.transit")

# Density mode (recommended): transit duration bounds set from stellar density.
# stellar_density in g/cm³; min/max_exp_dur_frac scale the expected duration.
result = vt.Pipeline([
    cmd.BLS(0.1, 20.0,
            density_mode=True, stellar_density=1.4,
            min_exp_dur_frac=0.5, max_exp_dur_frac=1.5,
            nbins=200, nfreq=100000, npeaks=1, fittrap=True,
            extraparams=True,
            save_periodogram=True, save_model=True),
]).run(lc)
print(result.stats["BLS_Period_1_0"])    # 2.12334706
print(result.stats["BLS_SN_1_0"])        # signal-to-noise
print(result.stats["BLS_SDE_1_0"])       # signal detection efficiency
pgram = result.files["BLS_periodogram_0"]   # pd.DataFrame: frequency vs BLS power

# stellar_density can also be a per-star variable name read from the list file
# (see Per-star variables in the Pipeline docs; use run_batch/run_filelist with inlistvars)
result2 = vt.Pipeline([
    cmd.BLS(0.1, 20.0,
            density_mode=True, stellar_density=1.4,
            min_exp_dur_frac=0.5, max_exp_dur_frac=1.5,
            nbins=200, nfreq=100000, npeaks=3),
]).run(lc)

# Fixed fractional duration (r mode) — use when stellar density is unknown
result3 = vt.Pipeline([
    cmd.BLS(0.1, 20.0, rmin=0.01, rmax=0.1, nbins=200, nfreq=100000, npeaks=1),
]).run(lc)

# Log-uniform frequency grid with auto qmin adjustment
result4 = vt.Pipeline([
    cmd.BLS(0.5, 10.0,
            density_mode=True, stellar_density=1.4,
            min_exp_dur_frac=0.5, max_exp_dur_frac=1.5,
            nbins=200, nfreq=50000,
            freq_grid="steplogP", adjust_qmin=True, reduce_nbins=True),
]).run(lc)
```

---

### `BLSFixPer` — BLS with fixed period

```python
cmd.BLSFixPer(period, rmin=0.01, rmax=0.1, nbins=200, timezone=0,
              qmin=None, qmax=None,
              save_model=False, correct_lc=False, fittrap=False,
              maskpoints=None)
```

Searches for the best transit epoch and duration at a known (fixed) period. Useful when a period has already been identified from a prior BLS or LS run.

| Parameter | Type | Description |
|-----------|------|-------------|
| `qmin`, `qmax` | `float` or `None` | Use `q` mode for duration bounds (ingress-to-egress fraction). When set, emits `"q" qmin qmax` instead of `"r" rmin rmax`. |
| All others | — | Same as `BLS`. |

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/3.transit")

# Compute RMS, fit transit at fixed period, then check RMS on residuals
pipe = vt.Pipeline([
    cmd.rms(),
    cmd.BLSFixPer("fix 2.12345", rmin=0.01, rmax=0.1, nbins=200, fittrap=True),
    cmd.rms(),
])
result = pipe.run(lc)
print(result.stats["BLSFixPer_Period_1"])    # 2.12345
print(result.stats["BLSFixPer_Depth_1"])     # transit depth
print(result.stats["BLSFixPer_Qtran_1"])     # fractional duration
```

---

### `BLSFixDurTc` — BLS with fixed duration and epoch, searching for period

```python
cmd.BLSFixDurTc(duration, Tc,
                minper=0.1, maxper=100.0, nfreq=10000,
                timezone=0, npeaks=1,
                fixdepth=None, qgress=None,
                save_periodogram=False, save_model=False,
                correct_lc=False, fittrap=False,
                save_phcurve=False, ophcurve_phmin=0.0,
                ophcurve_phmax=1.0, ophcurve_phstep=0.005,
                save_jdcurve=False, ojdcurve_jdstep=0.02,
                maskpoints=None)
```

Runs a BLS period search with the transit duration and epoch (Tc) held fixed.
The period that maximises the BLS statistic is found and reported.

`duration` and `Tc` each accept:

- A float → `fix <value>` (same value for all light curves)
- `"fixcolumn <colname>"` → read from a named per-star column
- `"list"` or `"list column <N>"` → read from an input-list column

| Parameter | Description |
|---|---|
| `duration` | Transit duration (days). |
| `Tc` | Epoch of transit center (JD/BJD). |
| `minper`, `maxper` | Period search range (days). |
| `nfreq` | Number of trial frequencies. |
| `timezone` | Add to JD to get local date (0 for UTC/BJD). |
| `npeaks` | Number of peaks to report. |
| `fixdepth` | Fix transit depth to this value (or column/list spec); `None` to optimise. |
| `qgress` | Fractional ingress/egress duration (requires `fixdepth`). |
| `save_periodogram` | Write BLS power spectrum to output dir. |
| `save_model` | Write best-fit transit model to output dir. |
| `correct_lc` | Subtract transit model from light curve. |
| `fittrap` | Fit a trapezoidal transit instead of a box. |

Output columns (suffix `N` is the pipeline command index):

| Column | Description |
|---|---|
| `BLSFixDurTc_Duration_N` | Fixed transit duration used. |
| `BLSFixDurTc_Tc_N` | Fixed epoch used. |
| `BLSFixDurTc_Period_1_N` | Best-fit period (first peak). |
| `BLSFixDurTc_SN_1_N` | Signal-to-noise of best peak. |
| `BLSFixDurTc_Depth_1_N` | Best-fit transit depth. |
| `BLSFixDurTc_Qtran_1_N` | Fractional transit duration. |
| `BLSFixDurTc_deltaChi2_1_N` | Δχ² of the best transit. |

```python
import pyvartools as vt
from pyvartools import commands as cmd

lc = vt.LightCurve.from_file("EXAMPLES/2")

# Search for the period with duration and Tc fixed
pipe = vt.Pipeline([
    cmd.BLSFixDurTc(duration=0.05, Tc=2450000.1,
                    minper=0.5, maxper=10.0, nfreq=5000,
                    timezone=0, npeaks=1),
])
result = pipe.run(lc)
print(result.stats["BLSFixDurTc_Period_1_0"])     # best-fit period
print(result.stats["BLSFixDurTc_SN_1_0"])         # signal-to-noise
print(result.stats["BLSFixDurTc_Depth_1_0"])      # transit depth
```

---

### `BLSFixPerDurTc` — BLS with fixed period, duration, and epoch

```python
cmd.BLSFixPerDurTc(period, duration, Tc,
                   timezone=0,
                   fixdepth=None, qgress=None,
                   save_model=False, correct_lc=False, fittrap=False,
                   save_phcurve=False, ophcurve_phmin=0.0,
                   ophcurve_phmax=1.0, ophcurve_phstep=0.005,
                   save_jdcurve=False, ojdcurve_jdstep=0.02,
                   maskpoints=None)
```

Computes BLS transit statistics for a fully specified signal — no period
search is performed. The period, duration, and Tc are all fixed; the depth
is optimised by default (or also fixed if `fixdepth` is given).

`period`, `duration`, and `Tc` each accept a float, `"fixcolumn <colname>"`,
or `"list"` / `"list column <N>"` (same forms as `BLSFixDurTc`).

| Parameter | Description |
|---|---|
| `period` | Transit period (days). |
| `duration` | Transit duration (days). |
| `Tc` | Epoch of transit center (JD/BJD). |
| `timezone` | Add to JD to get local date (0 for UTC/BJD). |
| `fixdepth` | Fix transit depth (or column/list spec); `None` to optimise. |
| `qgress` | Fractional ingress/egress duration (requires `fixdepth`). |
| `save_model` | Write best-fit transit model to output dir. |
| `correct_lc` | Subtract transit model from light curve. |
| `fittrap` | Fit a trapezoidal transit instead of a box. |

Output columns (suffix `N` is the pipeline command index):

| Column | Description |
|---|---|
| `BLSFixPerDurTc_Period_N` | Period used. |
| `BLSFixPerDurTc_Duration_N` | Duration used. |
| `BLSFixPerDurTc_Tc_N` | Epoch used. |
| `BLSFixPerDurTc_Depth_N` | Best-fit (or fixed) transit depth. |
| `BLSFixPerDurTc_Qtran_N` | Fractional transit duration. |
| `BLSFixPerDurTc_deltaChi2_N` | Δχ² of the transit signal. |
| `BLSFixPerDurTc_SN_N` | Signal-to-noise. |
| `BLSFixPerDurTc_fraconenight_N` | Fraction of Δχ² from one night. |
| `BLSFixPerDurTc_MeanMag_N` | Out-of-transit mean magnitude. |

```python
import pyvartools as vt
from pyvartools import commands as cmd

lc = vt.LightCurve.from_file("EXAMPLES/2")

# Evaluate BLS statistics at a fully fixed transit signal
pipe = vt.Pipeline([
    cmd.BLSFixPerDurTc(period=2.12345, duration=0.05, Tc=2450000.1,
                       timezone=0, correct_lc=False),
])
result = pipe.run(lc)
print(result.stats["BLSFixPerDurTc_Depth_0"])      # transit depth
print(result.stats["BLSFixPerDurTc_deltaChi2_0"])  # Δχ² of signal
print(result.stats["BLSFixPerDurTc_SN_0"])         # signal-to-noise
```

---

### `Phase` — Phase-fold the light curve

```python
cmd.Phase(period="ls", T0=None, phasevar=None, startphase=None)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `period` | `float` or `str` | Period to fold on. Can be a number, or a keyword like `"ls"`, `"aov"`, `"bls"` to inherit the period found by a prior pipeline command. |
| `T0` | `float` or `None` | Reference epoch. |
| `phasevar` | `str` or `None` | Name of the output phase variable. |
| `startphase` | `float` or `None` | Starting phase (default 0). |

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/2")

# Phase-fold at a known period
result = vt.Pipeline([
    cmd.Phase(period=1.2354),
]).run(lc, capture_lc=True)
phased_lc = result.lc   # time column replaced by orbital phase

# Use period found by BLS, set mid-transit at phase 0.5
lc_transit = vt.LightCurve.from_file("EXAMPLES/3.transit")
pipe = vt.Pipeline([
    cmd.BLS(0.5, 5.0, rmin=0.01, rmax=0.1, nbins=200, nfreq=20000, npeaks=1),
    cmd.Phase(period="bls", T0="bls"),
    cmd.binlc(method="median", nbins=200),
])
result = pipe.run(lc_transit, capture_lc=True)
phase_binned_lc = result.lc
```

---

### `autocorrelation` — Autocorrelation function

```python
cmd.autocorrelation(start, stop, step, save_result=True, maskpoints=None)
```

Computes the discrete autocorrelation function of the magnitude series.

| Parameter | Type | Description |
|-----------|------|-------------|
| `start`, `stop`, `step` | `float` | Lag range and step size. |
| `save_result` | `bool`, `str`, or `Output` | Controls Python capture of the output file. Default `True` (captured as `result.files["autocorrelation_result_N"]`). See [Auxiliary output files](#auxiliary-output-files) and the note below. |
| `maskpoints` | `str` or `None` | Name of a mask variable. |

!!! note "File is always written"
    The vartools CLI always writes the autocorrelation file to disk — there is no CLI option to suppress it. Setting `save_result=False` only suppresses Python capture; the file is still written to a temp directory and discarded after the run completes.

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/2")

# Default (save_result=True): ACF captured into result.files
result = vt.Pipeline([
    cmd.autocorrelation(0.0, 10.0, 0.05),
]).run(lc)
acf = result.files["autocorrelation_result_0"]   # pd.DataFrame: time-lag vs autocorrelation

# save_result=False: file written to temp dir but not captured
result = vt.Pipeline([
    cmd.autocorrelation(0.0, 10.0, 0.05, save_result=False),
]).run(lc)
# result.files has no "autocorrelation_result_0"

# Write to a specific directory and capture (Mode 2)
from pyvartools import Output
result = vt.Pipeline([
    cmd.autocorrelation(0.0, 10.0, 0.05,
                        save_result=Output("EXAMPLES/OUTDIR1", capture=True)),
]).run(lc)
acf = result.files["autocorrelation_result_0"]   # from EXAMPLES/OUTDIR1/

# Batch — result.files["autocorrelation_result_0"] is a list of DataFrames
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 4)]
batch = vt.Pipeline([
    cmd.autocorrelation(0.0, 10.0, 0.05),
]).run_batch(lcs)
acfs = batch.files["autocorrelation_result_0"]   # list of DataFrames, one per LC
```

---

### `dftclean` — CLEAN periodogram

```python
cmd.dftclean(nbeam, maxfreq=None, save_dspec=False, save_wfunc=False,
             save_cspec=False, gain=0.1, SNlimit=3.0, npeaks=None,
             finddirtypeaks=None, finddirtypeaks_clip=None,
             finddirtypeaks_clipiter=None,
             outcbeam=False, useampspec=False, verboseout=False,
             maskpoints=None)
```

Computes the CLEAN deconvolution periodogram (Roberts et al. 1987). Output files: `"dftclean_dspec_N"` (dirty spectrum), `"dftclean_wfunc_N"` (window function), `"dftclean_cspec_N"` (CLEAN spectrum), `"dftclean_cbeam_N"` (CLEAN beam, when `outcbeam` is set).

| Parameter | Type | Description |
|-----------|------|-------------|
| `finddirtypeaks` | `int` or `None` | Number of peaks to find in the dirty spectrum before CLEANing. |
| `finddirtypeaks_clip` | `float` or `None` | Sigma-clipping threshold for dirty-peak finding. |
| `finddirtypeaks_clipiter` | `int` or `None` | Number of sigma-clipping iterations. |
| `outcbeam` | `bool`, `str`, or `Output` | Write the CLEAN beam to a file. Captured as `result.files["dftclean_cbeam_N"]`. |
| `useampspec` | `bool` | Use the amplitude spectrum instead of power spectrum for CLEANing. |
| `verboseout` | `bool` | Write verbose diagnostics to the output file. |

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/2")

# Compute DFT CLEAN periodogram and find the top peak
result = vt.Pipeline([
    cmd.dftclean(4, maxfreq=10.0, npeaks=1, save_dspec=True),
]).run(lc)
print(result.stats["DFTCLEAN_DSPEC_PEAK_FREQ_0_0"])  # 0.81189711 cycles/day
dspec = result.files["dftclean_dspec_0"]   # pd.DataFrame: frequency vs power
```

---

### `wwz` — Weighted Wavelet Z-transform

```python
cmd.wwz(maxfreq="auto", freqsamp=None, tau0="auto", tau1="auto",
        dtau="auto", c=0.0125, save_transform=False,
        save_maxtransform=False,
        transform_format=None, transform_name=None,
        maxtransform_name=None, maskpoints=None)
```

Time-frequency analysis for non-stationary signals. Output files: `"wwz_transform_N"` (full time-frequency map), `"wwz_maxtransform_N"` (maximum power vs. time).

| Parameter | Type | Description |
|-----------|------|-------------|
| `transform_format` | `str` or `None` | Output format for the full transform: `"fits"` or `"pm3d"`. Only applies when `save_transform` is set. |
| `transform_name` | `str` or `None` | Naming format string for the full-transform output file (e.g. `"%s.wwz"`). |
| `maxtransform_name` | `str` or `None` | Naming format string for the max-transform output file. |

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/8")

result = vt.Pipeline([
    cmd.wwz(maxfreq=2.0, freqsamp=0.25, tau0="auto", tau1="auto",
            dtau=0.1, save_transform=True, save_maxtransform=True),
]).run(lc)
maxt = result.files["wwz_maxtransform_0"]   # pd.DataFrame: time vs peak frequency/power
```

---

### `GetLSAmpThresh` — LS amplitude threshold

```python
cmd.GetLSAmpThresh(period="ls", minp=0.1, thresh=10.0,
                   mode="harm", nharm=1, nsubharm=0,
                   listfile=None, noGLS=False)
```

Estimates the signal amplitude required to achieve a given detection threshold. Used in injection-recovery studies.

| Parameter | Type | Description |
|-----------|------|-------------|
| `mode` | `str` | Signal model: `"harm"` (Fourier series, default) or `"file"` (read template from `listfile`). |
| `nharm` | `int` | Number of harmonics (only when `mode="harm"`). |
| `nsubharm` | `int` | Number of sub-harmonics (only when `mode="harm"`). |
| `listfile` | `str` or `None` | Path to template signal file (required when `mode="file"`). |
| `noGLS` | `bool` | Use classical Lomb-Scargle instead of the generalized (GLS) form. |

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/2")

# Run LS, fit harmonic (fitonly), then compute minimum detectable amplitude
pipe = vt.Pipeline([
    cmd.LS(0.1, 10.0, 0.1, npeaks=1),
    cmd.Killharm("ls", nharm=0, nsubharm=0, fitonly=True),
    cmd.GetLSAmpThresh("ls", minp=0.1, thresh=-100.0, nharm=0, nsubharm=0),
])
result = pipe.run(lc)
print(result.stats["LS_Period_1_0"])           # 1.23440877
print(result.stats["LS_MinimumAmplitude_2"])   # 0.00248 mag
```

---

## Manipulation and statistics

### `clip` — Sigma clipping

```python
cmd.clip(sigclip, iterative=True, niter=None, median=False,
         markclip=None, noinitmark=False, maskpoints=None)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `sigclip` | `float` | Clipping threshold in units of standard deviation. |
| `iterative` | `bool` | Repeat clipping until no points are removed (default `True`). |
| `niter` | `int` or `None` | Clip at most this many times (overrides `iterative`). |
| `median` | `bool` | Clip relative to the median instead of the mean. |
| `markclip` | `str` or `None` | Variable name to record clipping mask (1 = kept, 0 = clipped). |

CLI equivalent: `-clip sigclip 1|0 [niter N] [median] ...`

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/5")

# Compute RMS, apply 3-sigma clipping, compute RMS of clipped LC
pipe = vt.Pipeline([
    cmd.rms(),
    cmd.clip(3.0),
    cmd.rms(),
])
result = pipe.run(lc)
print(result.stats["Nclip_1"])    # 51 points removed
print(result.stats["RMS_2"])      # RMS after clipping
```

---

### `rms` — Root mean square

```python
cmd.rms(maskpoints=None)
```

Computes the RMS and weighted RMS of the light curve. Output statistics: `Mean_Mag_N`, `RMS_N`, `Weighted_RMS_N`.

**Examples**

```python
# Single light curve
lc = vt.LightCurve.from_file("EXAMPLES/2")
result = vt.Pipeline([cmd.rms()]).run(lc)
print(result.stats["Mean_Mag_0"])
print(result.stats["RMS_0"])

# Batch: compute RMS for all 10 example light curves
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 11)]
batch = vt.Pipeline([cmd.rms()]).run_batch(lcs)
print(batch.stats[["Name", "Mean_Mag_0", "RMS_0", "Expected_RMS_0"]])
```

---

### `rmsbin` — Binned RMS

```python
cmd.rmsbin(nbin, bintimes, maskpoints=None)
```

Computes the binned RMS at a set of specified timescales. `bintimes` is a list of timescales (e.g. `[0.02, 0.05, 0.1]` days).

**Examples**

```python
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 11)]

# Compute binned RMS at 5, 10, 60, 1440, and 14400 minute windows
bintimes_days = [5/1440, 10/1440, 60/1440, 1.0, 10.0]
batch = vt.Pipeline([
    cmd.rmsbin(5, bintimes_days),
]).run_batch(lcs)
print(batch.stats)
```

---

### `chi2` — Chi-squared statistic

```python
cmd.chi2(maskpoints=None)
```

Computes the chi-squared statistic of the light curve. Output: `Chi2_N`, `Chi2_per_dof_N`.

**Examples**

```python
# Batch: chi-squared for all example light curves
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 11)]
batch = vt.Pipeline([cmd.chi2()]).run_batch(lcs)
print(batch.stats[["Name", "Chi2_0", "Weighted_Mean_Mag_0"]])
```

---

### `chi2bin` — Binned chi-squared

```python
cmd.chi2bin(nbin, bintimes, maskpoints=None)
```

Same as `rmsbin` but computes chi-squared instead of RMS at each timescale.

**Examples**

```python
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 11)]
bintimes_days = [5/1440, 10/1440, 60/1440, 1.0, 10.0]
batch = vt.Pipeline([
    cmd.chi2bin(5, bintimes_days),
]).run_batch(lcs)
print(batch.stats)
```

---

### `stats` — Generic statistics

```python
cmd.stats(variables, statistics, maskpoints=None)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `variables` | `str` or list of `str` | Variable name(s) to compute statistics for. |
| `statistics` | `str` or list of `str` | Statistics to compute: `"mean"`, `"median"`, `"stddev"`, `"min"`, `"max"`, etc. |

```python
cmd.stats("mag", "mean,median,stddev")
cmd.stats(["mag", "err"], ["mean", "stddev"])
```

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/3")

# Compute percentile and distribution statistics after adding Gaussian noise
pipe = vt.Pipeline([
    cmd.expr("mag2=mag+0.01*gauss()"),
    cmd.stats(
        ["mag", "mag2"],
        ["mean", "weightedmean", "median", "stddev", "MAD",
         "kurtosis", "skewness", "pct10", "pct90", "max", "min"],
    ),
])
result = pipe.run(lc)
print(result.stats["STATS_mag_MEAN_1"])
print(result.stats["STATS_mag_MEDIAN_1"])
print(result.stats["STATS_mag2_STDDEV_1"])
```

---

### `sortlc` — Sort observations

```python
cmd.sortlc(var=None, reverse=False)
```

Sort observations by time (default) or by any other variable. `reverse=True` sorts in descending order.

---

### `restricttimes` / `restoretimes` — Time windowing

```python
cmd.restricttimes(mode="JDrange", minJD=None, maxJD=None,
                  JDfilename=None, expression=None, exclude=False,
                  markrestrict=None, noinitmark=False)
cmd.restoretimes(prior_command=1)
```

`restricttimes` discards observations outside a time window or list. `restoretimes` undoes a prior restriction.

Modes: `"JDrange"` (min/max), `"JDlist"` (file of times), `"expr"` (boolean expression).

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/3")

# Restrict to a JD window, compute stats, then restore
pipe = vt.Pipeline([
    cmd.stats("t", ["min", "max"]),
    cmd.restricttimes(mode="JDrange", minJD=53740, maxJD=53750),
    cmd.stats("t", ["min", "max"]),
    cmd.restoretimes(prior_command=1),
    cmd.stats("t", ["min", "max"]),
])
result = pipe.run(lc)
print(result.stats["STATS_t_MIN_0"])   # original time range
print(result.stats["STATS_t_MIN_2"])   # restricted range
print(result.stats["STATS_t_MIN_4"])   # restored (original again)

# Restrict using a boolean expression on magnitude
pipe2 = vt.Pipeline([
    cmd.restricttimes(mode="expr",
                      expression="(mag>10.16311)&&(mag<10.17027)"),
])
result2 = pipe2.run(lc, capture_lc=True)
```

---

### `binlc` — Bin in time

```python
cmd.binlc(method="average", binsize=None, nbins=None,
          time_output="tcenter", bincolumns=None,
          bincolumnsonly=False, T0=None, firstbinshift=None,
          maskpoints=None)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `method` | `str` | `"average"`, `"median"`, or `"weightedaverage"`. |
| `binsize` | `float` | Bin width in time units. Either `binsize` or `nbins` must be given. |
| `nbins` | `int` | Number of bins (alternative to `binsize`). |
| `time_output` | `str` | Output time for each bin: `"tcenter"`, `"taverage"`, `"tmedian"`, or `"tnoshrink"`. |

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/2")

# Bin in time with 0.01-day bins (median combination)
result = vt.Pipeline([
    cmd.binlc(method="median", binsize=0.01),
]).run(lc, capture_lc=True)
binned_lc = result.lc

# Phase-fold then bin into 100 phase bins
pipe = vt.Pipeline([
    cmd.LS(0.1, 10.0, 0.1, npeaks=1),
    cmd.Phase(period="ls"),
    cmd.binlc(method="median", nbins=100),
])
result = pipe.run(lc, capture_lc=True)
phase_binned_lc = result.lc
```

---

### `expr` — Analytic expression

```python
cmd.expr(expression, outputcolumns=None)
```

Evaluate an analytic expression to create or update a variable. The expression string has the form `varname=formula`, e.g. `"residual=mag-model"`. Variables defined with `expr` can be passed to subsequent commands via `outputcolumns`.

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/1")

# Apply a mathematical transform in-place
result = vt.Pipeline([
    cmd.expr("mag=sqrt(mag+5)"),
]).run(lc, capture_lc=True)

# Convert to flux, normalise by median, then compute statistics
pipe = vt.Pipeline([
    cmd.expr("flux=10^(-0.4*(mag-25.0))"),
    cmd.stats("flux", ["median"]),
    cmd.expr("flux=flux/STATS_flux_MEDIAN_1"),
    cmd.stats(["flux", "mag"], ["median", "stddev"]),
])
result = pipe.run(lc)
print(result.stats["STATS_flux_MEDIAN_1"])   # original median flux
print(result.stats["STATS_flux_MEDIAN_3"])   # ≈ 1.0 after normalisation
```

---

### `decorr` — Decorrelation

```python
cmd.decorr(correct_lc=True, zeropointterm=1, subtractfirstterm=0,
           global_files=None, lc_columns=None, save_model=False,
           maskpoints=None)
```

Decorrelates the light curve against external trend vectors (polynomial fit). `global_files` is a list of `(filename, polynomial_order)` tuples for external trend files; `lc_columns` is a list of `(column_number, polynomial_order)` tuples for light-curve-internal trends. `save_model` accepts `bool`, `str`, or `Output` — see [Auxiliary output files](#auxiliary-output-files); the model file is captured as `result.files["decorr_model_N"]`.

**Examples**

```python
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 11)]

# Decorrelate against the time column using a quadratic polynomial
pipe = vt.Pipeline([
    cmd.rms(),
    cmd.decorr(correct_lc=True, zeropointterm=1, subtractfirstterm=1,
               lc_columns=[(1, 2)]),   # column 1 (JD), order 2
    cmd.rms(),
])
batch = pipe.run_batch(lcs)
print(batch.stats[["Name", "RMS_0", "RMS_2"]])
```

---

### `FFT` / `IFFT` — Fast Fourier Transform

```python
cmd.FFT(input_real, input_imag, output_real, output_imag)
cmd.IFFT(input_real, input_imag, output_real, output_imag)
```

Compute the FFT or inverse FFT of two variables (real and imaginary parts). Results are stored in the named output variables.

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/11")

# High-pass Fourier filter on a uniformly sampled light curve
pipe = vt.Pipeline([
    cmd.FFT("mag", "NULL", "fftreal", "fftimag"),
    cmd.rms(),
    # Zero low-frequency components (below 1/500 of full spectrum)
    cmd.expr("fftreal=(NR>(Npoints_1/500.0))*(NR<(Npoints_1*499.0/500.0))*fftreal"),
    cmd.expr("fftimag=(NR>(Npoints_1/500.0))*(NR<(Npoints_1*499.0/500.0))*fftimag"),
    cmd.IFFT("fftreal", "fftimag", "mag_filter", "NULL"),
])
result = pipe.run(lc, capture_lc=True)
```

---

### `resample` — Resample onto a new time grid

```python
cmd.resample(method="linear",
             left=None, right=None,
             nbreaks=None, order=None,
             file_times=None, file_column=None,
             gaps=None,
             tstart=None, tstop=None, delt=None, Npoints=None)
```

Methods: `"nearest"`, `"linear"`, `"spline"`, `"splinemonotonic"`, `"bspline"`. Specify the new grid with `delt` (step size), `Npoints` (number of points), or `file_times` (times from a file).

| Parameter | Type | Description |
|-----------|------|-------------|
| `left` | `float` or `None` | First-derivative boundary condition at the left edge of the spline. Only for `method="spline"` or `"splinemonotonic"`. |
| `right` | `float` or `None` | First-derivative boundary condition at the right edge of the spline. Only for `method="spline"` or `"splinemonotonic"`. |
| `nbreaks` | `int` or `None` | Number of interior break points for B-spline fitting. Only for `method="bspline"`. |
| `order` | `int` or `None` | Polynomial order of the B-spline. Only for `method="bspline"`. |
| `file_times` | `str` or `None` | Path to a file containing times (emits `"file" "fix" path`), or a string starting with `"list"` for list-column mode (e.g. `"list column 2"`). |
| `file_column` | `int` or `None` | Column number in the times file. Only used with the `"fix"` form (path given). |
| `gaps` | `str` or `None` | Gaps option string, e.g. `"fix"`, `"list"`, `"fixcolumn myvar"`, `"expr someexpr"`. |

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/2")

# Linear interpolation with default time grid
result = vt.Pipeline([
    cmd.resample(method="linear"),
]).run(lc, capture_lc=True)

# Monotonic spline onto a fixed time grid with 1000 points
result2 = vt.Pipeline([
    cmd.resample(method="splinemonotonic",
                 tstart=53726, tstop=53756, Npoints=1000),
]).run(lc, capture_lc=True)

# B-spline with 20 break points, order 3
result3 = vt.Pipeline([
    cmd.resample(method="bspline", nbreaks=20, order=3, Npoints=500),
]).run(lc, capture_lc=True)
```

---

### `difffluxtomag` / `fluxtomag` — Flux conversions

```python
cmd.difffluxtomag(mag_constant, offset=0.0, magcolumn=None)
cmd.fluxtomag(mag_constant, offset=0.0)
```

Convert differential or absolute flux to magnitude. `mag_constant` is the zero-point magnitude. `offset` is added to the flux before conversion.

**Examples**

```python
# Convert a FITS light curve in flux units to magnitudes
lc = vt.LightCurve.from_file(
    "EXAMPLES/example_flux.fits",
    t_col="TIME", mag_col="SAP_FLUX", err_col="SAP_FLUX_ERR",
)
result = vt.Pipeline([
    cmd.fluxtomag(25.0, offset=0.0),
]).run(lc, capture_lc=True)
```

---

### `medianfilter` — Median filtering

```python
cmd.medianfilter(time, method="median", replace=False)
```

Apply a sliding-window median (or mean) filter with window width `time`. `replace=True` writes the smoothed values back into the magnitude column.

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/1")

# High-pass and low-pass median filter: save LC, process both ways
pipe = vt.Pipeline([
    cmd.chi2(),
    cmd.savelc(),
    cmd.medianfilter(0.05),           # residuals (high-pass output)
    cmd.chi2(),
    cmd.restorelc(savenumber=1),
    cmd.medianfilter(0.05, replace=True),   # smoothed LC (low-pass)
    cmd.chi2(),
])
result = pipe.run(lc)
print(result.stats["Chi2_0"])   # original
print(result.stats["Chi2_3"])   # after high-pass filter
print(result.stats["Chi2_7"])   # after low-pass filter
```

---

### `changeerror` / `changevariable` — Column operations

```python
cmd.changeerror(maskpoints=None)
cmd.changevariable(column, var)
```

`changeerror` rescales measurement uncertainties. `changevariable` copies a named variable into one of the standard columns: `"t"`, `"mag"`, `"err"`, or `"id"`. This is useful when reading FITS files with non-standard column mappings.

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/4")

# Replace formal errors with RMS, then verify chi2 → ~1
pipe = vt.Pipeline([
    cmd.chi2(),
    cmd.changeerror(),
    cmd.chi2(),
])
result = pipe.run(lc)
print(result.stats["Chi2_0"])   # 5.19874 (original)
print(result.stats["Chi2_2"])   # ≈ 1.0 (after rescaling errors to RMS)

# Use changevariable to store LS phase then revert to time
pipe2 = vt.Pipeline([
    cmd.LS(0.1, 100.0, 0.1, npeaks=1),
    cmd.expr("phase=t"),
    cmd.changevariable("t", "phase"),
    cmd.Phase(period="ls"),
    cmd.changevariable("t", "t"),  # revert so output sorts by time
])
```

---

### `savelc` / `restorelc` / `copylc` — LC state management

```python
cmd.savelc()
cmd.restorelc(savenumber=1, vars=None)
cmd.copylc(ncopies)
```

`savelc` saves a snapshot of the current light curve state. `restorelc` restores a previous snapshot (useful for running multiple analyses on the same underlying data). `copylc` duplicates the light curve N times in the output.

**Examples**

```python
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 11)]

# Run LS and BLS with different clipping levels; restore between them
pipe = vt.Pipeline([
    cmd.savelc(),
    cmd.clip(5.0),
    cmd.savelc(),
    cmd.LS(0.1, 100.0, 0.1, npeaks=3, clip=5.0, clipiter=1),
    cmd.aov(0.1, 100.0, 0.1, 0.01, npeaks=1, clip=5.0, clipiter=1),
    cmd.restorelc(savenumber=1),   # back to pre-5σ-clip state
    cmd.clip(10.0),
    cmd.BLS(0.1, 20.0, rmin=0.01, rmax=0.1, nbins=200, nfreq=10000, npeaks=1),
    cmd.restorelc(savenumber=2),   # back to 5σ-clipped state
    cmd.changeerror(),
])
batch = pipe.run_batch(lcs, nthreads=4)
print(batch.stats)

# copylc: bootstrap false-alarm probability via noise copies
lc = vt.LightCurve.from_file("EXAMPLES/2")
pipe_bs = vt.Pipeline([
    cmd.LS(0.1, 10.0, 0.1, npeaks=1),
    cmd.copylc(100),
    cmd.expr("mag=err*gauss()"),
    cmd.LS(0.1, 10.0, 0.1, npeaks=1),
])
batch_bs = pipe_bs.run_batch([lc])
```

---

### `Jstet` — Stetson J-statistic

```python
cmd.Jstet(timescale, dates, maskpoints=None)
```

Computes the J variability index (Stetson 1996). Requires a dates file listing observation epochs.

**Examples**

```python
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 11)]
batch = vt.Pipeline([
    cmd.Jstet(0.5, "EXAMPLES/dates_tfa"),
]).run_batch(lcs)
print(batch.stats[["Name", "Jstet_0", "Kurtosis_0", "Lstet_0"]])
```

---

## Model fitting

### `Killharm` — Harmonic series subtraction

```python
cmd.Killharm(period="ls", nharm=3, nsubharm=0, save_model=False,
             fitonly=False, output_format=None, clip=None,
             maskpoints=None)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `period` | `float` or `str` | Period to fit. Can be a number or `"ls"`, `"aov"`, `"bls"`, `"both"`, or `"fix val1 val2..."` for multiple periods. |
| `nharm` | `int` | Number of harmonics. |
| `nsubharm` | `int` | Number of sub-harmonics. |
| `save_model` | `bool`, `str`, or `Output` | Auxiliary file output. `True` captures as `result.files["Killharm_model_N"]`. See [Auxiliary output files](#auxiliary-output-files). |
| `fitonly` | `bool` | Fit the model but do not subtract it (statistics are still computed). |
| `output_format` | `str` or `None` | Coefficient output format: `"outampphase"`, `"outampradphase"`, `"outRphi"`, or `"outRradphi"`. |
| `clip` | `float` or `None` | Sigma-clipping threshold: fit, clip outliers, then refit. |

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/2")

# Run LS, fit and subtract the best sinusoid, compare chi2 before/after
pipe = vt.Pipeline([
    cmd.LS(0.1, 10.0, 0.1, npeaks=1),          # index 0
    cmd.rms(),                                   # index 1
    cmd.chi2(),                                  # index 2
    cmd.Killharm("ls", nharm=1, nsubharm=0, save_model=True),  # index 3
    cmd.rms(),                                   # index 4
    cmd.chi2(),                                  # index 5
])
result = pipe.run(lc)
print(result.stats["Chi2_2"])         # before: 1709.50
print(result.stats["Chi2_5"])         # after:    6.51
model = result.files["Killharm_model_3"]   # pd.DataFrame of the fitted harmonic curve
```

---

### `linfit` — Linear combination fitting

```python
cmd.linfit(function, paramlist,
           modelvar=None,
           reject=None, reject_usemad=False, reject_iter=False,
           reject_fixednum=None,
           correct_lc=False, save_model=False,
           model_nameformat=None, fitmask=None)
```

Fit a linear combination of analytic basis functions. `function` is a vartools expression string; `paramlist` is a comma-separated list of free parameter names and initial values.  `fitmask` is the name of a mask variable; observations where the variable is non-zero are excluded from the fit (CLI token: `fitmask`). `save_model` accepts `bool`, `str`, or `Output` — see [Auxiliary output files](#auxiliary-output-files); the model file is captured as `result.files["linfit_model_N"]`.

| Parameter | Type | Description |
|-----------|------|-------------|
| `modelvar` | `str` or `None` | Variable name used to store the best-fit model values on the light curve. |
| `reject` | `float` or `None` | Sigma-clipping threshold: fit, clip outliers beyond this threshold (in σ), then refit. |
| `reject_usemad` | `bool` | Use MAD instead of standard deviation for the scatter estimate during rejection. |
| `reject_iter` | `bool` | Iteratively reject and refit until no more points are clipped. |
| `reject_fixednum` | `int` or `None` | Maximum number of rejection/refit iterations (requires `reject_iter=True`). |
| `model_nameformat` | `str` or `None` | Format string for the model output filename (e.g. `"%s.linfit.model"`). |

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/1")

# Fit a quadratic polynomial, using minimum time as reference epoch
pipe = vt.Pipeline([
    cmd.stats("t", ["min"]),
    cmd.expr("t0=STATS_t_MIN_0"),
    cmd.linfit("a*(t-t0)^2+b*(t-t0)+c", "a,b,c"),
])
result = pipe.run(lc)
print(result.stats["Linfit_a_2"])   # quadratic coefficient
print(result.stats["Linfit_b_2"])   # linear coefficient
print(result.stats["Linfit_c_2"])   # constant offset
```

---

### `nonlinfit` — Non-linear least-squares fitting

```python
cmd.nonlinfit(function, paramlist,
              optimizer="amoeba",
              linfit_params=None, errors=None,
              covariance=None, priors=None, constraints=None,
              amoeba_tolerance=None, amoeba_maxsteps=None,
              mcmc_naccept=None, mcmc_nlinkstotal=None,
              mcmc_fracburnin=None, mcmc_eps=None,
              mcmc_skipamoeba=False, mcmc_maxmemstore=None,
              mcmc_outchains=False, mcmc_chains_format=None,
              mcmc_chains_printevery=None,
              correct_lc=False, save_model=False,
              model_nameformat=None, modelvar=None,
              fitmask=None)
```

Fit an arbitrary analytic function using Nelder-Mead (`"amoeba"`) or MCMC (`"mcmc"`). `paramlist` has the form `name:initial[:step[:min:max]], ...`.  `fitmask` excludes non-zero-masked points from the fit. `save_model` accepts `bool`, `str`, or `Output` — see [Auxiliary output files](#auxiliary-output-files); the model file is captured as `result.files["nonlinfit_model_N"]`. MCMC chains are captured as `result.files["nonlinfit_chains_N"]`.

| Parameter | Type | Description |
|-----------|------|-------------|
| `linfit_params` | `str` or `None` | Space-separated list of parameter names to solve analytically (linear sub-problem solved exactly at each step). |
| `errors` | `str` or `None` | Expression or variable name providing per-point measurement errors (overrides the default error column). |
| `covariance` | `str` or `None` | Covariance model tokens, e.g. `"squareexp amp_v rho_v"` (passed verbatim after the `covariance` keyword). |
| `priors` | `str` or `None` | Prior expression string (passed verbatim after the `priors` keyword). |
| `constraints` | `str` or `None` | Constraint expression string (passed verbatim after `constraints`). |
| `amoeba_tolerance` | `float` or `None` | Convergence tolerance for the amoeba optimizer. |
| `amoeba_maxsteps` | `int` or `None` | Maximum number of steps for the amoeba optimizer. |
| `mcmc_naccept` | `int` or `None` | Number of accepted links to collect (mutually exclusive with `mcmc_nlinkstotal`). |
| `mcmc_nlinkstotal` | `int` or `None` | Total number of MCMC links (accepted + rejected). |
| `mcmc_fracburnin` | `float` or `None` | Fraction of links to discard as burn-in. |
| `mcmc_eps` | `float` or `None` | Initial step size for the MCMC proposal distribution. |
| `mcmc_skipamoeba` | `bool` | Skip the initial amoeba pre-optimisation in MCMC mode. |
| `mcmc_maxmemstore` | `int` or `None` | Maximum number of chain links to hold in memory simultaneously. |
| `mcmc_outchains` | `bool`, `str`, or `Output` | Write MCMC chain files. Captured as `result.files["nonlinfit_chains_N"]`. |
| `mcmc_chains_format` | `str` or `None` | Naming format string for chain output files. |
| `mcmc_chains_printevery` | `int` or `None` | Write every Nth accepted link to the chain file. |
| `model_nameformat` | `str` or `None` | Format string for the model output filename. |
| `modelvar` | `str` or `None` | Variable name used to store the best-fit model values on the light curve. |

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/1")

# Fit a sinusoid with free period using amoeba.
# paramlist format: "name=initial:step" (no spaces around commas)
result = vt.Pipeline([
    cmd.nonlinfit(
        "a*sin(2*pi*(t-t0)/P)+c",
        "a=0.01:0.001,P=1.23:0.01,t0=0.0:0.1,c=10.0:0.01",
        optimizer="amoeba",
        amoeba_tolerance=1e-8,
        amoeba_maxsteps=10000,
        save_model=True,
    ),
]).run(lc)
print(result.stats["Nonlinfit_P_BestFit_0"])   # best-fit period
model = result.files["nonlinfit_model_0"]

# MCMC posterior sampling (chains written to disk via mcmc_outchains)
result2 = vt.Pipeline([
    cmd.nonlinfit(
        "a*sin(2*pi*(t-t0)/P)+c",
        "a=0.01:0.001,P=1.23:0.01,t0=0.0:0.1,c=10.0:0.01",
        optimizer="mcmc",
        mcmc_naccept=10000,
        mcmc_fracburnin=0.2,
        mcmc_outchains="EXAMPLES/OUTDIR1",
    ),
]).run(lc)
print(result2.stats["Nonlinfit_P_MEDIAN_0"])   # posterior median period
```

---

### `MandelAgolTransit` — Mandel-Agol transit model

```python
cmd.MandelAgolTransit(P0, T00, r0=0.1, a0=10.0, inclination=90.0,
                      bimpact=None,
                      e0=0.0, omega0=0.0, mconst0=0.0,
                      ld_type="quad", ld_coeffs=None,
                      fitephem=1, fitr=1, fita=1,
                      fitinclterm=1, fite=0, fitomega=0,
                      fitmconst=1, fitldcoeffs=None,
                      rv_file=None, rv_model_file=None,
                      K0=None, gamma0=None, fitK=0, fitgamma=0,
                      correct_lc=False, save_model=False,
                      modelvar=None,
                      save_phcurve=False, save_jdcurve=False)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `P0`, `T00` | `float` | Initial period (days) and transit epoch (BJD). |
| `r0` | `float` | Initial planet-to-star radius ratio. |
| `a0` | `float` | Initial semi-major axis in stellar radii. |
| `inclination` | `float` | Initial orbital inclination in degrees. Used when `bimpact` is `None`. |
| `bimpact` | `float` or `None` | Initial impact parameter. When set, replaces `inclination` in the CLI (`"b" bimpact` instead of `"i" inclination`). |
| `ld_type` | `str` | Limb-darkening law: `"quad"` or `"nonlin"`. |
| `ld_coeffs` | list | Initial limb-darkening coefficients. Default `[0.3, 0.3]`. |
| `fitephem` | `int` | Fit the transit epoch and period together (0 = fixed, 1 = free).  The CLI `fitephem` flag controls both simultaneously — there is no separate period-only flag. |
| `fitr`, `fita`, `fitinclterm`, etc. | `int` | Toggle fitting of each parameter (0 = fixed, 1 = free). |
| `rv_file` | `str` or `None` | Path to an RV data file (columns: JD, RV, RV-error). When provided, simultaneous RV fitting is enabled. |
| `rv_model_file` | `str` or `None` | Output path for the best-fit RV model. |
| `K0` | `float` or `None` | Initial RV semi-amplitude (km/s). |
| `gamma0` | `float` or `None` | Initial systemic RV (km/s). |
| `fitK`, `fitgamma` | `int` | Fit K and γ respectively (0 = fixed, 1 = free). |
| `modelvar` | `str` or `None` | Variable name used to store the best-fit model on the light curve (requires `save_model` to be set). |
| `save_phcurve` | `bool`, `str`, or `Output` | Auxiliary file output. `True` captures as `result.files["MandelAgolTransit_phcurve_N"]`. See [Auxiliary output files](#auxiliary-output-files). |
| `save_jdcurve` | `bool`, `str`, or `Output` | Auxiliary file output. `True` captures as `result.files["MandelAgolTransit_jdcurve_N"]`. See [Auxiliary output files](#auxiliary-output-files). |

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/3.transit")

# Run BLS to find transit, then fit Mandel-Agol model
pipe = vt.Pipeline([
    cmd.BLS(0.5, 5.0, rmin=0.01, rmax=0.1, nbins=200, nfreq=20000, npeaks=1),
    cmd.MandelAgolTransit(
        P0="bls", T00="bls",     # initialise from BLS result
        ld_type="quad", ld_coeffs=[0.3471, 0.3180],
        fitephem=1, fitr=1, fita=1, fitinclterm=1,
        fite=0, fitomega=0, fitmconst=1,
        save_model=True, save_phcurve=True,
    ),
])
result = pipe.run(lc)
print(result.stats["MandelAgolTransit_Period_1"])   # 2.12328176
print(result.stats["MandelAgolTransit_r_1"])        # Rp/R* ≈ 0.098
model   = result.files["MandelAgolTransit_model_1"]    # fitted model LC
phcurve = result.files["MandelAgolTransit_phcurve_1"]  # phase-folded model
```

---

### `SoftenedTransit` — Softened (trapezoidal) transit

```python
cmd.SoftenedTransit(init_params="bls", fitephem=1, fiteta=1,
                    fitcval=1, fitdelta=1, fitmconst=1,
                    correct_lc=False, save_model=False,
                    fit_harm=0, fit_harm_method=None,
                    fit_harm_nharm=None, fit_harm_nsubharm=None)
```

Fits a softened trapezoid transit model. `init_params` can be `"bls"` or `"blsfixper"` to initialise from a prior BLS result, or a tuple `(P0, T00, eta0, delta0, mconst0, cval0)`.

| Parameter | Type | Description |
|-----------|------|-------------|
| `fit_harm` | `int` | Harmonic fitting flag. `0` = no harmonic component (default). When > 0, a harmonic series is fitted simultaneously. |
| `fit_harm_method` | `str` or `None` | Method for harmonic fitting, e.g. `"aov"`. Only used when `fit_harm > 0`. |
| `fit_harm_nharm` | `int` or `None` | Number of harmonics for the harmonic component. |
| `fit_harm_nsubharm` | `int` or `None` | Number of sub-harmonics for the harmonic component. |

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/3.transit")

pipe = vt.Pipeline([
    cmd.BLS(0.5, 5.0, rmin=0.01, rmax=0.1, nbins=200, nfreq=20000, npeaks=1),
    cmd.SoftenedTransit(init_params="bls",
                        fitephem=1, fiteta=1, fitcval=1,
                        fitdelta=1, fitmconst=1,
                        correct_lc=False, save_model=True),
])
result = pipe.run(lc)
print(result.stats["SoftenedTransit_Period_1"])     # 2.12322112
print(result.stats["SoftenedTransit_chi2perdof_1"])
model = result.files["SoftenedTransit_model_1"]   # SoftenedTransit is at index 1
```

---

### `Starspot` — Starspot model

```python
cmd.Starspot(period="ls",
             a0=0.1, b0=0.5, alpha0=20.0, i0=85.0,
             chi0=30.0, psi00=0.0, mconst0=0.0,
             fit_period=1, fit_a=1, fit_b=1, fit_alpha=1,
             fit_i=1, fit_chi=1, fit_psi=1, fit_mconst=1,
             correct_lc=False, save_model=False)
```

Fits a Dorren (1987) two-spot model to photometric variability.

| Parameter | Type | Description |
|-----------|------|-------------|
| `period` | `float` or `str` | Starting period. Can be a number or a keyword like `"ls"`, `"aov"`. |
| `a0` | `float` | Initial spot fractional radius (default `0.1`). |
| `b0` | `float` | Initial spot latitude in radians (default `0.5`). |
| `alpha0` | `float` | Initial spot longitude in degrees (default `20.0`). |
| `i0` | `float` | Initial stellar inclination in degrees (default `85.0`). |
| `chi0` | `float` | Initial spot contrast (default `30.0`). |
| `psi00` | `float` | Initial spot phase offset (default `0.0`). |
| `mconst0` | `float` | Initial magnitude offset (default `0.0`). |
| `fit_period` | `int` | Fit the period (0 = fixed, 1 = free; default `1`). |
| `fit_a`, `fit_b`, `fit_alpha`, `fit_i`, `fit_chi`, `fit_psi`, `fit_mconst` | `int` | Fit each model parameter (0 = fixed, 1 = free; all default `1`). |

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/3.starspot")

# Determine rotation period with AOV then fit Dorren starspot model
pipe = vt.Pipeline([
    cmd.aov(0.1, 10.0, 0.1, 0.01, npeaks=5, nbin=20),
    cmd.Starspot(
        period="aov",
        a0=0.0298, b0=0.08745, alpha0=20.0, i0=85.0,
        chi0=30.0, psi00=0.0, mconst0=-1.0,
        fit_period=1, fit_a=1, fit_b=1, fit_alpha=1,
        fit_i=1, fit_chi=1, fit_psi=1, fit_mconst=1,
        correct_lc=False,
        save_model=True,
    ),
])
result = pipe.run(lc)
print(result.stats["Starspot_Period_1"])
print(result.stats["Starspot_chi2perdof_1"])
model = result.files["Starspot_model_1"]   # Starspot is at index 1
```

---

### `microlens` — Microlensing model

```python
cmd.microlens(f0=None, f1=None, u0=None, t0=None, tmax=None,
              correct_lc=False, save_model=False,
              f0_step=None, f0_novary=False,
              f1_step=None, f1_novary=False,
              u0_step=None, u0_novary=False,
              t0_step=None, t0_novary=False,
              tmax_step=None, tmax_novary=False)
```

Fits a standard Paczynski microlensing light curve. Each parameter (`f0`, `f1`, `u0`, `t0`, `tmax`) can be a float (free-fit initial value), `"auto"` (vartools auto-estimate), a string passthrough (e.g. `"fixcolumn colname"`, `"list column 3"`), or `None` (omit). Use `{name}_step` to set the initial step size and `{name}_novary=True` to hold a parameter fixed during fitting.

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/4.microlensinject")

result = vt.Pipeline([
    cmd.microlens(f0="auto", f1="auto", u0="auto",
                  t0="auto", tmax="auto", save_model=True),
]).run(lc)
print(result.stats["Microlens_u0_0"])
print(result.stats["Microlens_tmax_0"])
print(result.stats["Microlens_chi2perdof_0"])
model = result.files["microlens_model_0"]
```

---

### `TFA` — Trend Filtering Algorithm

```python
cmd.TFA(trendlist, dates_file, pixelsep, correct_lc=True,
        save_coeffs=False, save_model=False, xycol=None,
        clip=None, usemedian=False, useMAD=False,
        readformat=None, trend_coeff_priors=None,
        weight_by_template_stddev=False, fitmask=None,
        outfitmask=None)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `trendlist` | `str` | Path to a file listing the trend (template) light curves, one per line. |
| `dates_file` | `str` | Path to the dates file (one epoch per line, matching the observation cadence). |
| `pixelsep` | `float` | Maximum pixel separation for selecting trend stars. Stars further than this threshold are excluded. |
| `correct_lc` | `bool` | Subtract the TFA model from the light curve. Default `True`. |
| `save_coeffs` | `bool`, `str`, or `Output` | Auxiliary file output. `True` captures as `result.files["TFA_coeffs_N"]`. See [Auxiliary output files](#auxiliary-output-files). |
| `xycol` | `(int, int)` or `None` | Column numbers `(xcol, ycol)` for pixel coordinates in the trend list. |
| `clip` | `float` or `None` | Sigma-clipping threshold during TFA fitting. |

---

### `TFA_SR` — TFA with signal reconstruction

```python
cmd.TFA_SR(trendlist, dates_file, pixelsep, dotfafirst=1,
           tfathresh=0.001, maxiter=10, signal_mode="bin",
           signal_params=None, signal_period=None,
           correct_lc=True, decorr_params=None, ...)
```

Simultaneous TFA detrending and signal reconstruction. `signal_mode` controls the signal model: `"bin"` (phase-binned), `"signal"` (from file), or `"harm"` (harmonic series with `signal_params=(nharm, nsubharm)`).

| Parameter | Type | Description |
|-----------|------|-------------|
| `signal_period` | `float`, `str`, or `None` | Period sub-option for `"bin"` or `"harm"` signal modes. Float emits `"period" val`; string passthrough (e.g. `"ls"` to use Lomb-Scargle period). |
| `decorr_params` | `str` or `None` | Raw token string for simultaneous EPD decorrelation, e.g. `"0 2 col1 1 col2 2"` (iterative_flag Nlcterms lccolumn1 lcorder1 ...). |

---

### `SYSREM` — Systematic noise removal

```python
cmd.SYSREM(ninput_color, ninput_airmass, initial_airmass_file,
           sigma_clip1=5.0, sigma_clip2=5.0, saturation=1e9,
           correct_lc=True, save_model=False, save_trends=False,
           useweights=1, col=None)
```

Tamuz et al. (2005) SYSREM algorithm. Removes systematic effects correlated with colour and airmass across an ensemble of light curves. Both `save_model` and `save_trends` accept `bool`, `str`, or `Output` — see [Auxiliary output files](#auxiliary-output-files). The model is captured as `result.files["SYSREM_model_N"]` and the trend vectors as `result.files["SYSREM_trends_N"]`.

---

### `addnoise` — Add synthetic noise

```python
cmd.addnoise(noise_type="white", sig_white=0.001, rho=None,
             sig_red=None, nu=None, gamma=None, bintime=None)
```

Adds synthetic noise drawn from a specified covariance model. `noise_type` controls the model; all amplitude/timescale parameters accept either a float or a vartools variable-name string.

| `noise_type` | Parameters | Description |
|---|---|---|
| `"white"` | `sig_white` | Independent Gaussian noise. |
| `"squareexp"` | `rho`, `sig_red`, `sig_white`, `bintime` | Squared-exponential (Gaussian) covariance. `bintime` enables integrated covariance. |
| `"exp"` | `rho`, `sig_red`, `sig_white`, `bintime` | Exponential covariance. |
| `"matern"` | `nu`, `rho`, `sig_red`, `sig_white` | Matérn covariance with smoothness `nu`. |
| `"wavelet"` | `gamma`, `sig_red`, `sig_white` | Wavelet (1/f-like) noise. |

| Parameter | Type | Description |
|-----------|------|-------------|
| `sig_white` | `float` or `str` | White noise amplitude (all models). |
| `rho` | `float` or `str` or `None` | Correlation timescale (`"squareexp"`, `"exp"`, `"matern"`). |
| `sig_red` | `float` or `str` or `None` | Red noise amplitude (all correlated models). |
| `nu` | `float` or `str` or `None` | Matérn smoothness parameter (`"matern"` only). |
| `gamma` | `float` or `str` or `None` | Wavelet decay exponent (`"wavelet"` only). |
| `bintime` | `float` or `str` or `None` | Bin integration time (`"squareexp"`, `"exp"`). |

**Examples**

```python
import numpy as np
import pyvartools as vt
from pyvartools import commands as cmd

# Build a zero-magnitude light curve with EXAMPLES/1 time sampling
lc_ref = vt.LightCurve.from_file("EXAMPLES/1")
t = lc_ref.t
lc_blank = vt.LightCurve.from_arrays(t, np.zeros_like(t), np.full_like(t, 0.005))

# Simulate wavelet (1/f-like) red noise + white noise
result = vt.Pipeline([
    cmd.addnoise(noise_type="wavelet", sig_red=0.005, sig_white=0.005),
]).run(lc_blank, capture_lc=True)
noisy_lc = result.lc

# Squared-exponential red noise with 0.01-day correlation timescale
result2 = vt.Pipeline([
    cmd.addnoise(noise_type="squareexp", rho=0.01,
                 sig_red=0.005, sig_white=0.001),
]).run(lc_blank, capture_lc=True)
```

---

### `Injectharm` — Inject a harmonic signal

```python
cmd.Injectharm(period, amplitude, nharm=1, phase=0.0,
               nsubharm=0, save_model=False)
```

Injects a sinusoidal (or multi-harmonic) signal with the specified period and amplitude. Used for injection-recovery tests.

> **Design note**: This class exposes only the most common injection modes.
> *Amplitude*: `"ampfix"` and `"amplogrand"` only. For `"amprand"` or `"amplist"` use `cmd.Raw()`.
> *Period*: `"fix"` and `"logrand"` (via the `period` parameter). For `"list"` or `"rand"` period modes use `cmd.Raw()`.

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/3")

# Inject a sine wave at a random log-uniform period, then try to recover it
pipe = vt.Pipeline([
    cmd.Injectharm(period="logrand 1.0 5.0",
                   amplitude="amplogrand 0.001 0.1",
                   nharm=0, phase="phaserand",
                   save_model=True),
    cmd.LS(0.5, 10.0, 0.1, npeaks=1),
])
result = pipe.run(lc)
print(result.stats["Injectharm_Period_0"])   # injected period
print(result.stats["LS_Period_1_1"])         # recovered period
```

---

### `Injecttransit` — Inject a transit signal

```python
cmd.Injecttransit(period, Rp, Mp, phase, sini, Mstar, Rstar,
                  e=0.0, omega=0.0,
                  hk=False, h=0.0, k=0.0,
                  dilute=None,
                  ld_type="quad", ld_coeffs=None, save_model=False)
```

All positional parameters accept floats (emits `{prefix}fix val`) or strings (passed through verbatim, e.g. `"Plogrand 0.2 2.0"`, `"phaserand"`, `"sinirand"`).

| Parameter | Type | Description |
|-----------|------|-------------|
| `period` | `float` or `str` | Orbital period. Float → `"Pfix val"`. String passthrough (e.g. `"Plogrand 0.2 2.0"`). |
| `Rp` | `float` or `str` | Planet-to-star radius ratio. Float → `"Rpfix val"`. String passthrough (e.g. `"Rplogrand 0.05 0.15"`). |
| `Mp` | `float` or `str` | Planet mass in solar masses. |
| `phase` | `float` or `str` | Orbital phase of transit centre (0–1). String e.g. `"phaserand"`. |
| `sini` | `float` or `str` | Sine of the orbital inclination. String e.g. `"sinirand"`. |
| `Mstar`, `Rstar` | `float` or `str` | Stellar mass (M☉) and radius (R☉). |
| `e` | `float` or `str` | Eccentricity (used in `eomega` mode, the default). |
| `omega` | `float` or `str` | Argument of periastron (used in `eomega` mode). |
| `hk` | `bool` | When `True`, use the `hk` eccentricity parameterisation (`h = e sin ω`, `k = e cos ω`) instead of `eomega`. |
| `h`, `k` | `float` or `str` | `h` and `k` eccentricity components. Used when `hk=True`. |
| `dilute` | `float`, `str`, or `None` | Dilution factor. Float → `["dilute", "fix", val]`. String passthrough (e.g. `"list"` or `"fix 0.5"`). |
| `ld_coeffs` | list | Limb-darkening coefficients. Default `[0.3, 0.3]`. |

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/4")

# Inject a Jupiter-sized transit at a random period, then search with BLS
pipe = vt.Pipeline([
    cmd.Injecttransit(
        period="Plogrand 0.2 2.0",
        Rp="Rpfix 0.1",     # Rp/R*
        Mp="Mpfix 0.001",   # M_sun
        phase="phaserand",
        sini="sinirand",
        Mstar="Mstarfix 1.0",
        Rstar="Rstarfix 1.0",
        ld_type="quad",
        ld_coeffs=[0.3471, 0.3180],
        save_model=True,
    ),
    cmd.BLS(0.1, 5.0, rmin=0.01, rmax=0.1, nbins=200, nfreq=20000, npeaks=1),
])
result = pipe.run(lc)
print(result.stats["Injecttransit_Period_0"])   # injected period
print(result.stats["BLS_Period_1_1"])           # recovered period
```

---

### `findblends` — Search for blended transits

```python
cmd.findblends(matchrad, period="list", radec=False, nharm=1,
               xycol=None, starlist=None, zeromag=None,
               nofluxconvert=False, save_matches=False)
```

Searches for nearby stars that could be the source of a blended transit signal. `matchrad` is the search radius in arcseconds. `save_matches` accepts `bool`, `str`, or `Output` — see [Auxiliary output files](#auxiliary-output-files); the matched-star list is captured as `result.files["findblends_matches_N"]`.

| Parameter | Type | Description |
|-----------|------|-------------|
| `period` | `float` or `str` | Transit period source. Valid values: `"list"` (from input list column), `"fix <period>"` (fixed value), `"fixcolumn <col>"` (column name/number). Default `"list"`. |
| `radec` | `bool` | Use RA/Dec coordinates instead of pixel coordinates for matching. |
| `xycol` | `tuple` or `None` | `(xcol, ycol)` — column names or numbers for pixel-coordinate columns. |
| `starlist` | `str` or `None` | Path to an additional star list file used for blending candidates. |
| `zeromag` | `float` or `None` | Zero-point magnitude used for flux conversion. |
| `nofluxconvert` | `bool` | Skip the flux-to-magnitude conversion step. |

---

## Miscellaneous

### `columnsuffix` — Set output column suffix

```python
cmd.columnsuffix(suffix)
```

By default, vartools appends a 0-based command-index suffix to output variable names (e.g. `LS_Period_1_0`, `LS_Period_1_1`). `columnsuffix` replaces the numeric suffix with a user-supplied string for all subsequent commands, making column names predictable regardless of pipeline position.

```python
lc = vt.LightCurve.from_file("EXAMPLES/2")
pipe = vt.Pipeline([
    cmd.columnsuffix("ls"),
    cmd.LS(0.5, 10.0, 1e-3),
])
result = pipe.run(lc)
period = float(result.stats["LS_Period_1_ls"])
```

If your pipeline contains multiple period-search commands with different suffixes, you can read each result unambiguously:

```python
lc = vt.LightCurve.from_file("EXAMPLES/2")
pipe = vt.Pipeline([
    cmd.columnsuffix("ls"),
    cmd.LS(0.5, 10.0, 1e-3),
    cmd.columnsuffix("aov"),
    cmd.aov(0.5, 10.0, 1e-3, 4.0),
])
result = pipe.run(lc)
ls_period  = float(result.stats["LS_Period_1_ls"])
aov_period = float(result.stats["Period_1_aov"])
```

---

### `Raw` — Pass arbitrary CLI arguments

```python
cmd.Raw(args)
```

An escape hatch for commands not yet wrapped in Python, or for experimental options. Pass either a string (split on whitespace) or a list of tokens:

```python
cmd.Raw("-LS 0.1 10.0 0.1 5 1 0")
cmd.Raw(["-LS", "0.1", "10.0", "0.1", "5", "1", "0"])
```

---

### `converttime` — Time system conversion

```python
cmd.converttime(input_format, output_format, ra=None, dec=None,
                input_subtract=None, output_subtract=None,
                input_sys=None, output_sys=None, ephemfile=None,
                leapsecfile=None)
```

Convert between time systems. `input_format` / `output_format` are `"mjd"`, `"jd"`, `"hjd"`, or `"bjd"`. For HJD/BJD conversions, provide `ra` and `dec` in degrees. `input_sys` / `output_sys` can be `"tdb"` or `"utc"`.

**Examples**

```python
lc = vt.LightCurve.from_file("EXAMPLES/1")

# Convert JD (minus 2400000) to HJD for a known sky position
result = vt.Pipeline([
    cmd.converttime(
        input_format="jd",
        output_format="hjd",
        ra=88.079166,
        dec=32.5533,
        input_subtract=2400000.0,
    ),
]).run(lc, capture_lc=True)
hjd_lc = result.lc
```

---

### `binlc` — Bin in time

See the [Manipulation section](#binlc-bin-in-time) above.

---

### `addfitskeyword` — Add a FITS keyword

```python
cmd.addfitskeyword(keyword, dtype, value, comment=None,
                   hdu=None, mode=None, combinelc=None)
```

Add a FITS header keyword to the output statistics table. `dtype` is one of `"TDOUBLE"`, `"TINT"`, `"TLONG"`, `"TSTRING"`. `value` can be a Python scalar (auto-wrapped with `"fix"`) or a full vartools token string like `"var myvar"`.

---

### `match` — Match against a catalog

```python
cmd.match(catalog, matchcolumn, addcolumns, missing="nanmissing",
          source="file", inlist_column=None, skipnum=None,
          skipchar=None, delimiter=None, opencommand=None)
```

Match each light curve against rows in a catalog file and add columns from the catalog. `addcolumns` is a comma-separated list like `"ra:2,dec:3"`.

| Parameter | Type | Description |
|-----------|------|-------------|
| `catalog` | `str` | Path to the catalog file. Ignored when `source="inlist"`. |
| `source` | `str` | `"file"` (default) — one catalog for all LCs; `"inlist"` — each LC specifies its own catalog in a list column. |
| `inlist_column` | `str` or `int` or `None` | Column number/name in the input list that holds per-LC catalog paths. Required when `source="inlist"`. |
| `matchcolumn` | `str` | Column specification for matching, e.g. `"t:1"` (variable:column) or `"1"` (column number). |
| `addcolumns` | `str` | Comma-separated `varname:colnum[:dtype]` specs for columns to import. |
| `missing` | `str` | How to handle unmatched rows: `"cullmissing"`, `"nanmissing"`, or `"missingval <value>"`. |

---

### `o` — Output light curve

```python
cmd.o(filename=None, nameformat=None, columnformat=None,
      fits=False, noclobber=False, copyheader=False,
      namecommand=None, namefromlist=None, delimiter=None,
      logcommandline=False, capture=False, key="o")
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `filename` | `str` or `None` | Output file path (or directory in batch mode). Required unless `capture=True`. |
| `nameformat` | `str` or `None` | Format string for output filenames in batch mode, e.g. `"file_%s_%05d.txt"` (`%s` = LC basename, `%d` = sequence number). |
| `columnformat` | `str` or `None` | Output column spec, e.g. `"t:%17.9f,mag:%9.5f,err:%9.5f"`. |
| `fits` | `bool` | Write output in FITS binary table format. |
| `noclobber` | `bool` | Do not overwrite an existing output file. |
| `copyheader` | `bool` | Copy the FITS header from the input file to the output file. |
| `namecommand` | `str` or `None` | Shell command used to generate the output filename dynamically. |
| `namefromlist` | `bool`, `str`, or `None` | Derive output filename from the input list. `True` uses the default column; a string specifies a column name/number (emits `namefromlist column <col>`). |
| `delimiter` | `str` or `None` | Column delimiter character for the output file (default: whitespace). |
| `logcommandline` | `bool` | Write the full vartools command line into the output file header. |
| `capture` | `bool` | If `True`, capture the written light curve into `result.files[key]`. For single-LC runs this is a `LightCurve`; for batch runs it is a list of `LightCurve` objects. When `filename=None`, the output goes to a temporary file that is cleaned up automatically. Default `False`. |
| `key` | `str` | Key under which the captured LC(s) appear in `result.files`. Default `"o"`. Use a unique key when the pipeline contains more than one `cmd.o(capture=True)`. |

`cmd.o` can be used in three modes:

- **Write to disk only** (`filename="path"`, `capture=False`): existing behaviour; the LC is saved to disk.
- **Capture only** (`capture=True`, `filename=None`): the LC is written to a temp file, captured into `result.files[key]`, and the temp file is cleaned up automatically.
- **Write and capture** (`capture=True`, `filename="path"`): the LC is both saved to disk and captured into `result.files[key]`.

**Examples**

```python
import pyvartools as vt
from pyvartools import commands as cmd

lc = vt.LightCurve.from_file("EXAMPLES/2")

# Capture the intermediate light curve state (no file left on disk)
pipe = vt.Pipeline([
    cmd.clip(5.0),
    cmd.o(capture=True, key="clipped"),   # captured, no disk file
    cmd.LS(0.1, 10.0, 0.1, npeaks=1),
])
result = pipe.run(lc)
clipped_lc = result.files["clipped"]     # LightCurve after sigma-clipping
print(result.stats["LS_Period_1_2"])     # clip=0, o=1, LS=2

# Write to disk AND capture
result2 = vt.Pipeline([
    cmd.clip(5.0),
    cmd.o(filename="EXAMPLES/OUTDIR1/2.clipped", capture=True, key="clipped"),
]).run(lc)
# File written to disk and also available as result2.files["clipped"]

# Multiple intermediate snapshots
pipe3 = vt.Pipeline([
    cmd.clip(5.0),
    cmd.o(capture=True, key="after_clip"),
    cmd.medianfilter(0.05),
    cmd.o(capture=True, key="after_filter"),
])
result3 = pipe3.run(lc)
after_clip   = result3.files["after_clip"]
after_filter = result3.files["after_filter"]

# Batch: result.files["o"] is a list of LightCurves, one per input LC
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 11)]
batch = vt.Pipeline([
    cmd.LS(0.1, 100.0, 0.1, npeaks=1),
    cmd.expr("phase=t"),
    cmd.changevariable("t", "phase"),
    cmd.Phase(period="ls"),
    cmd.o(capture=True, key="phased"),
]).run_batch(lcs)
phased_lcs = batch.files["phased"]   # list of 10 LightCurve objects

# Write to a named directory with custom column format (no capture)
vt.Pipeline([
    cmd.LS(0.1, 100.0, 0.1, npeaks=1),
    cmd.Phase(period="ls"),
    cmd.o("EXAMPLES/OUTDIR1",
          nameformat="file_%s_%05d_simout.txt",
          columnformat="t:%11.5f,phase:%8.5f,mag:%7.4f,err:%7.4f"),
]).run_batch(lcs)
```

---

### `ifcmd` — Conditional execution

```python
cmd.ifcmd(condition)
```

Conditionally execute the next command based on a vartools expression string. Corresponds to `-if condition` in the CLI.

**Examples**

```python
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 11)]

# Run different statistics depending on variability level
pipe = vt.Pipeline([
    cmd.rms(),
    cmd.ifcmd("RMS_0>10*Expected_RMS_0"),
        cmd.ifcmd("RMS_0 > 0.1"),
            cmd.stats("mag", "stddev"),
        cmd.ifcmd("else"),
            cmd.stats("mag", "pct30"),
        cmd.ifcmd("fi"),
    cmd.ifcmd("elif Npoints_0>3900"),
        cmd.stats("mag", "kurtosis"),
    cmd.ifcmd("else"),
        cmd.rms(),
    cmd.ifcmd("fi"),
])
batch = pipe.run_batch(lcs)
print(batch.stats)
```

!!! note
    Pass `"else"`, `"elif <condition>"`, and `"fi"` as the `condition` argument of `cmd.ifcmd` for those control-flow roles.

---

### `R` — Run R code

```python
cmd.R(command, fromfile=False, init=None, init_fromfile=False,
      vars=None, invars=None, outvars=None, outputcolumns=None,
      process_all_lcs=False, verbose=False, continueprocess=None)
```

Execute inline R code or an R script on each light curve. `vars` specifies variables to pass both into and out of R; `invars`/`outvars` allow separate control. `init` is R code run once before the batch loop begins.

**Examples**

```python
lcs = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 11)]

# Compute standard deviation via R for each light curve
batch = vt.Pipeline([
    cmd.R("b <- sd(mag)", invars="mag", outvars="b",
          outputcolumns="b"),
]).run_batch(lcs)
print(batch.stats[["Name", "R_b_0"]])

# Same computation but send all light curves to R at once
batch2 = vt.Pipeline([
    cmd.R(
        "b <- list(); for(i in 1:length(mag)) { b[[i]] <- sd(mag[[i]]) }",
        invars="mag", outvars="b", outputcolumns="b",
        process_all_lcs=True,
    ),
]).run_batch(lcs)
```

---

## User Extension Commands

vartools supports user-developed commands compiled as shared libraries (`.so` / `.la`). Examples bundled with the source tree include `-stitch` and `-macula` (see [Extension Commands](../cli/extensions.md) for the CLI reference). pyvartools provides a first-class integration layer for these extensions so they can be used from Python with the same convenience as built-in commands.

### How extensions are loaded

At the CLI level, an extension is loaded with `-L path/to/lib.so` placed immediately before the command flag:

```bash
vartools -l lc_list -L USERLIBS/src/stitch.so -stitch mag err mask lcnum median -tab
```

If the library is installed into the vartools userlibs data directory (e.g. via `make install` in `USERLIBS/src/`), vartools loads it automatically and the `-L` flag is not needed. pyvartools mirrors both behaviours.

### Usage pattern 1 — quick one-off (`UserCommand`)

`UserCommand` is the lowest-level entry point. Pass the library path, the command name, and the raw argument tokens:

```python
import pyvartools as vt

pipe = vt.Pipeline([
    vt.UserCommand(
        "USERLIBS/src/stitch.so",   # path to .so
        "stitch",                    # command name
        "mag err mask lcnum median", # raw args (str or list)
    )
])
result = pipe.run_batch(lcs)
```

When the library is installed and auto-loaded, omit the path:

```python
pipe = vt.Pipeline([
    vt.UserCommand(None, "stitch", "mag err mask lcnum median")
])
```

`UserCommand` constructor parameters:

| Parameter | Type | Description |
|-----------|------|-------------|
| `lib_path` | `str`, `Path`, or `None` | Path to the `.so` / `.la` file. `None` or `""` means the library is auto-loaded from the vartools userlibs directory. |
| `name` | `str`, optional | Command name (e.g. `"stitch"`). Inferred from the filename stem when `lib_path` is given. |
| `args` | `str` or `list[str]` | Raw CLI tokens passed after the command flag. A plain string is split on whitespace. Default: empty. |

Call `.help()` or `.examples()` on any instance to print the vartools help text for the command (requires the binary and library to be loadable):

```python
vt.UserCommand("USERLIBS/src/stitch.so", "stitch").help()
```

### Usage pattern 2 — named class (`load_userlib`)

`load_userlib()` creates a reusable `UserCommand` subclass with the library path and command name pre-bound. The returned class is functionally equivalent to a hand-written command wrapper and can be further subclassed.

```python
Stitch = vt.load_userlib("USERLIBS/src/stitch.so")

# Instantiate like any other command class
pipe = vt.Pipeline([Stitch("mag err mask lcnum median")])
result = pipe.run_batch(lcs)
```

`load_userlib()` parameters:

| Parameter | Type | Description |
|-----------|------|-------------|
| `lib_path` | `str` or `Path` | Path to the `.so` / `.la` file. Resolved to an absolute path. |
| `name` | `str`, optional | Command name. Defaults to the filename stem. |
| `cls_name` | `str`, optional | Python class name for the returned type. Defaults to the title-cased command name (e.g. `"Stitch"`). |

The class docstring is populated by running `vartools -L lib -help -name` at creation time, so `help(Stitch)` or `Stitch.help()` shows live vartools documentation.

### Usage pattern 3 — auto-discovery (`discover_userlibs`)

`discover_userlibs()` scans known directories and returns a `{name: cls}` dict of all installed extensions. The default search order is:

1. Paths in `$VARTOOLS_USERLIBS` (colon-separated).
2. `$prefix/share/vartools/userlibs/` derived from the binary location.
3. Any paths passed explicitly via `search_paths`.

```python
cmds = vt.discover_userlibs()          # auto-discover installed extensions
print(cmds)                            # {'stitch': <class Stitch>, ...}

pipe = vt.Pipeline([cmds["stitch"]("mag err mask lcnum median")])
result = pipe.run_batch(lcs)
```

Pass explicit paths to search additional directories:

```python
cmds = vt.discover_userlibs(search_paths=["USERLIBS/src"])
```

### Usage pattern 4 — full Python wrapper (subclass)

For a production-quality wrapper with named, typed arguments, subclass `UserCommand` directly:

```python
class Stitch(vt.UserCommand):
    """Fit and remove zero-point offsets between light curve segments."""

    def __init__(
        self,
        variables: str,
        errors: str,
        masks: str,
        lcnum_var: str,
        method: str = "median",
        lib_path: str = "/usr/local/share/vartools/userlibs/stitch.so",
    ):
        super().__init__(
            lib_path, "stitch",
            f"{variables} {errors} {masks} {lcnum_var} {method}",
        )

pipe = vt.Pipeline([
    Stitch("mag", "err", "mask", "lcnum", method="weightedmean"),
])
```

Alternatively, build the base class from the factory for a one-line definition:

```python
class Stitch(vt.load_userlib("USERLIBS/src/stitch.so", name="stitch")):
    def __init__(self, variables, errors, masks, lcnum, method="median"):
        super().__init__(f"{variables} {errors} {masks} {lcnum} {method}")
```

### Output statistics

Output statistics produced by user commands appear in `result.stats` automatically — vartools writes them to its standard output just like built-in commands, and pyvartools parses them in the same way. No special configuration is needed.

### Pipeline execution mode

Pipelines that contain `UserCommand` instances always run in **subprocess mode**, even when `libvartoolspipeline.so` is available. The in-process library does not support dynamically loaded extension libraries. This is handled transparently; no change to how you call the pipeline is required.

---

## Base class

All wrappers inherit from `VartoolsCommand`:

```python
from pyvartools._command import VartoolsCommand
```

Custom commands can be written by subclassing `VartoolsCommand` and implementing `_to_cli_args() → list[str]`. Override `_output_file_specs()` to declare output files that should be captured by the Pipeline.

```python
class MyCommand(VartoolsCommand):
    _vt_name = "mycommand"

    def __init__(self, param: float):
        self.param = param

    def _to_cli_args(self):
        return ["-mycommand", str(self.param)]
```
