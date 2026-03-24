# Period Finding

These recipes show how to search for periodic signals in a light curve using the Lomb-Scargle, AoV, and BLS algorithms.

---

## 1. Basic Lomb-Scargle search

Find the single strongest period between 0.5 and 10 days.

=== "CLI"

    ```bash
    vartools -i star.lc -LS 0.5 10.0 0.001 1 0 -oneline
    ```

    Flag breakdown:
    - `0.5 10.0` — period range (days)
    - `0.001` — frequency step as a fraction of 1/T
    - `1` — report the top 1 peak
    - `0` — do not save a periodogram file

    The output line contains `LS_Period_1_0`, `LS_Amplitude_1_0`, and `Log10_LS_Prob_1_0`.

=== "Python"

    ```python
    import pyvartools as vt
    from pyvartools import commands as cmd

    lc = vt.LightCurve.from_file("star.lc")
    pipe = vt.Pipeline([cmd.LS(0.5, 10.0, 1e-3, npeaks=1)])
    result = pipe.run(lc)

    print(f"Best period : {float(result.stats['LS_Period_1_0']):.5f} d")
    print(f"Amplitude   : {float(result.stats['LS_Amplitude_1_0']):.4f} mag")
    print(f"log10(FAP)  : {float(result.stats['Log10_LS_Prob_1_0']):.2f}")
    ```

---

## 2. LS with multiple peaks and SNR threshold

Report the top 5 peaks and assess their significance through the log₁₀ false-alarm probability (FAP).

=== "CLI"

    ```bash
    vartools -i star.lc -LS 0.1 20.0 0.001 5 0 -oneline
    ```

    The output contains `LS_Period_1_0` through `LS_Period_5_0` and the corresponding `Log10_LS_Prob` values.

=== "Python"

    ```python
    pipe = vt.Pipeline([cmd.LS(0.1, 20.0, 1e-3, npeaks=5)])
    result = pipe.run(lc)

    for i in range(1, 6):
        period = float(result.stats[f"LS_Period_{i}_0"])
        log_fap = float(result.stats[f"Log10_LS_Prob_{i}_0"])
        print(f"Peak {i}: P = {period:.5f} d,  log10(FAP) = {log_fap:.2f}")
        if log_fap < -5:
            print("  --> Significant detection")
    ```

!!! tip "Interpreting log10(FAP)"
    `Log10_LS_Prob` is the log₁₀ of the false-alarm probability under the null hypothesis (white noise). A value below −2 is noteworthy; below −5 is a strong candidate. Values near 0 indicate the peak is consistent with noise.

---

## 3. Save the periodogram to file

Write the full frequency–power spectrum for plotting or further analysis.

=== "CLI"

    ```bash
    vartools -i star.lc -LS 0.1 10.0 0.001 1 operiodogram /output 0 -oneline
    ```

    The file `/output/star.lc.ls` is written with columns: frequency, LS power.

=== "Python"

    ```python
    pipe = vt.Pipeline([
        cmd.LS(0.1, 10.0, 1e-3, npeaks=1, save_periodogram=True),
    ])
    result = pipe.run(lc)

    pgram = result.files["LS_periodogram_0"]   # pd.DataFrame
    pgram.columns = ["frequency", "power"]

    import matplotlib.pyplot as plt
    plt.plot(1.0 / pgram["frequency"], pgram["power"])
    plt.xlabel("Period (d)")
    plt.ylabel("LS power")
    plt.savefig("periodogram.png", dpi=150)
    ```

    The DataFrame is returned in memory — no output directory needs to be specified.

---

## 4. BLS transit search

Search for periodic box-shaped dips consistent with a transiting planet.

=== "CLI"

    ```bash
    vartools -i star.lc \
        -BLS r 0.01 0.1 0.5 20.0 optimal 1.0 200 0 1 operiodogram /output 0 1 0 \
        -oneline
    ```

    Key output columns: `BLS_Period_1_0`, `BLS_SDE_1_0` (signal detection efficiency), `BLS_SN_1_0` (signal-to-noise), `BLS_Depth_1_0`, `BLS_Qtran_1_0` (fractional duration), `BLS_Tc_1_0` (transit epoch).

=== "Python"

    ```python
    pipe = vt.Pipeline([
        cmd.BLS(
            minper=0.5,
            maxper=20.0,
            rmin=0.01,
            rmax=0.1,
            nbins=200,
            npeaks=1,
            save_periodogram=True,
        ),
    ])
    result = pipe.run(lc)

    period  = float(result.stats["BLS_Period_1_0"])
    sde     = float(result.stats["BLS_SDE_1_0"])
    depth   = float(result.stats["BLS_Depth_1_0"])
    qtran   = float(result.stats["BLS_Qtran_1_0"])
    tc      = float(result.stats["BLS_Tc_1_0"])

    print(f"Period:    {period:.5f} d")
    print(f"SDE:       {sde:.1f}")
    print(f"Depth:     {depth*1000:.2f} mmag")
    print(f"Duration:  {qtran * period * 24:.2f} h")
    print(f"Epoch:     {tc:.4f} BJD")
    ```

!!! tip "BLS threshold"
    A signal detection efficiency (SDE) above ~7 is typically considered a transit candidate worthy of follow-up. Always check the phase-folded light curve visually.

---

## 5. Analysis of Variance (AoV)

AoV projects the phase-folded light curve onto phase bins and measures how well the best period explains the variance. It makes fewer assumptions about signal shape than LS.

=== "CLI"

    ```bash
    vartools -i star.lc -aov 0.1 10.0 0.001 4.0 1 0 -oneline
    ```

    Flag order: `minp maxp subsample finetune npeaks save`. Key output: `AOV_Period_1_0`, `AOV_Statistic_1_0`.

=== "Python"

    ```python
    pipe = vt.Pipeline([
        cmd.aov(
            minp=0.1,
            maxp=10.0,
            subsample=1e-3,
            finetune=4.0,
            npeaks=1,
        ),
    ])
    result = pipe.run(lc)
    print(f"AoV period: {float(result.stats['AOV_Period_1_0']):.5f} d")
    ```

!!! note "AoV vs LS"
    LS is more powerful for sinusoidal signals; AoV is preferred when the waveform is non-sinusoidal (e.g. eclipsing binaries) or when you suspect strong aliasing.

---

## 6. Multi-harmonic AoV

For strongly non-sinusoidal signals — RR Lyrae, Cepheids, W UMa contact binaries — a multi-harmonic model provides a much better fit.

=== "CLI"

    ```bash
    vartools -i rr_lyrae.lc -aov_harm 6 0.1 1.0 0.001 4.0 1 0 -oneline
    ```

    The first argument after `-aov_harm` is the number of harmonics.

=== "Python"

    ```python
    pipe = vt.Pipeline([
        cmd.aov_harm(
            nharm=6,
            minp=0.1,
            maxp=1.0,
            subsample=1e-3,
            finetune=4.0,
            npeaks=1,
            save_periodogram=True,
        ),
    ])
    result = pipe.run(lc)
    print(f"AOV_harm period: {float(result.stats['AOV_HARM_Period_1_0']):.5f} d")

    pgram = result.files["aov_harm_periodogram_0"]
    ```

---

## 7. Using `columnsuffix` for readable output names

When a pipeline contains multiple period-search commands, the default numeric suffix can be confusing. Use `columnsuffix` to assign meaningful names.

=== "Python"

    ```python
    pipe = vt.Pipeline([
        cmd.columnsuffix("ls"),
        cmd.LS(0.5, 10.0, 1e-3),
        cmd.columnsuffix("aov"),
        cmd.aov(0.5, 10.0, 1e-3, 4.0),
    ])
    result = pipe.run(lc)

    ls_period  = float(result.stats["LS_Period_1_ls"])
    aov_period = float(result.stats["AOV_Period_1_aov"])
    print(f"LS period:  {ls_period:.5f} d")
    print(f"AoV period: {aov_period:.5f} d")
    ```

This is especially valuable in batch processing, where you iterate over `batch.stats` columns by name:

```python
batch = pipe.run_batch(lcs)
ls_periods = batch.stats["LS_Period_1_ls"]   # pd.Series, one value per LC
```

---

## 8. Chaining clip and period search

Sigma-clip outliers before running the period search to avoid false peaks driven by bad points.

=== "CLI"

    ```bash
    vartools -i star.lc \
        -clip 5.0 1 \
        -LS 0.5 10.0 0.001 1 0 \
        -oneline
    ```

    `-clip 5.0 1` clips iteratively at 5σ until convergence before the LS search.

=== "Python"

    ```python
    pipe = vt.Pipeline([
        cmd.clip(sigclip=5.0, iterative=True),
        cmd.LS(0.5, 10.0, 1e-3, npeaks=1),
    ])
    result = pipe.run(lc)
    print(f"Period (after clipping): {float(result.stats['LS_Period_1_0']):.5f} d")
    ```

To inspect the clipped light curve as well as the period:

```python
result = pipe.run(lc, capture_lc=True)
clipped_lc = result.lc
print(f"Points removed by clipping: {len(lc) - len(clipped_lc)}")
print(f"Period: {float(result.stats['LS_Period_1_0']):.5f} d")
```
