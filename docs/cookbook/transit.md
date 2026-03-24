# Transit Detection and Fitting

These recipes cover the full transit workflow: blind searching with BLS, detrending before searching, fitting a physical model, and injecting synthetic transits for recovery tests.

---

## 1. BLS transit search

A full BLS search with standard parameters suitable for typical ground-based photometry.

=== "CLI"

    ```bash
    vartools -i star.lc \
        -clip 5.0 1 \
        -BLS r 0.01 0.1 0.5 20.0 optimal 1.0 200 0 1 operiodogram /output 0 1 0 \
        -oneline
    ```

    Key flags: `r 0.01 0.1` sets the fractional duration range; `optimal 1.0` uses automatic frequency sampling; `200` phase bins; `1` peak reported; `operiodogram /output` saves the power spectrum.

=== "Python"

    ```python
    import pyvartools as vt
    from pyvartools import commands as cmd

    lc = vt.LightCurve.from_file("star.lc")

    pipe = vt.Pipeline([
        cmd.clip(5.0),
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

    period = float(result.stats["BLS_Period_1_0"])
    sde    = float(result.stats["BLS_SDE_1_0"])
    sn     = float(result.stats["BLS_SN_1_0"])
    depth  = float(result.stats["BLS_Depth_1_0"])
    qtran  = float(result.stats["BLS_Qtran_1_0"])
    tc     = float(result.stats["BLS_Tc_1_0"])

    print(f"Period:       {period:.5f} d")
    print(f"SDE:          {sde:.2f}")
    print(f"S/N:          {sn:.2f}")
    print(f"Depth:        {depth * 1000:.2f} mmag")
    print(f"Duration:     {qtran * period * 24:.2f} h")
    print(f"Transit epoch:{tc:.4f} BJD")

    # Plot the BLS power spectrum
    pgram = result.files["BLS_periodogram_1"]   # BLS is at index 1 (clip is at 0)
    pgram.columns = ["frequency", "power"]
    ```

---

## 2. Fitting a transit model (MandelAgolTransit)

After identifying a candidate period with BLS, fit the full Mandel-Agol transit model to measure the physical parameters.

=== "CLI"

    ```bash
    # Fix the period from BLS; fit epoch, radius ratio, impact parameter
    vartools -i star.lc \
        -BLS r 0.01 0.1 5.0 6.0 optimal 1.0 200 0 1 0 0 1 0 \
        -MandelAgolTransit bls 0.1 10.0 90.0 0.0 0.0 0.0 quad 0.3 0.3 \
            1 0 1 1 1 0 0 1 0 0 \
            1 omodel /output 0 0 \
        -oneline
    ```

    The `bls` keyword initialises the transit parameters from the prior BLS result.

=== "Python"

    ```python
    # Step 1: run BLS to find the period
    bls_pipe = vt.Pipeline([
        cmd.clip(5.0),
        cmd.BLS(minper=5.0, maxper=6.0, rmin=0.01, rmax=0.1),
    ])
    bls_result = bls_pipe.run(lc)

    period = float(bls_result.stats["BLS_Period_1_0"])
    tc     = float(bls_result.stats["BLS_Tc_1_0"])

    # Step 2: fit the Mandel-Agol model at the BLS period
    fit_pipe = vt.Pipeline([
        cmd.MandelAgolTransit(
            P0=period,
            T00=tc,
            r0=0.1,
            a0=10.0,
            inclination=88.0,
            e0=0.0,
            omega0=0.0,
            mconst0=0.0,
            ld_type="quad",
            ld_coeffs=[0.3, 0.3],
            fitephem=1,    # fit the epoch
            fitP=0,        # keep period fixed
            fitr=1,        # fit radius ratio
            fita=1,        # fit semi-major axis
            fitinclterm=1, # fit inclination
            fite=0,        # fix eccentricity = 0
            fitmconst=1,   # fit magnitude offset
            fitldcoeffs=[0, 0],  # keep LD coefficients fixed
            save_model=True,
            save_phcurve=True,
        ),
    ])
    fit_result = fit_pipe.run(lc)

    print(f"Fitted period:   {float(fit_result.stats['MandelAgolTransit_Period_1_0']):.6f} d")
    print(f"Radius ratio:    {float(fit_result.stats['MandelAgolTransit_Rp_1_0']):.4f}")
    print(f"Impact param:    {float(fit_result.stats['MandelAgolTransit_b_1_0']):.3f}")

    # Phase-folded model curve for plotting
    phcurve = fit_result.files["MandelAgolTransit_phcurve_0"]
    ```

---

## 3. Injecting a transit signal (Injecttransit)

Injection-recovery tests quantify pipeline sensitivity as a function of period, depth, and duration.

=== "CLI"

    ```bash
    vartools -i clean.lc \
        -Injecttransit fix 3.14159 Rpfix 0.1 Mpfix 0.001 phasefix 0.0 \
            sinifix 0.999 eomega efix 0 ofix 0 Mstarfix 1.0 Rstarfix 1.0 \
            quad 0.3 0.3 0 \
        -BLS r 0.01 0.1 1.0 10.0 optimal 1.0 200 0 1 0 0 1 0 \
        -oneline
    ```

=== "Python"

    ```python
    import numpy as np

    injection_pipe = vt.Pipeline([
        cmd.Injecttransit(
            period=3.14159,
            Rp=0.1,           # Rp/Rstar = 0.1 → ~1% depth
            Mp=0.001,         # negligible mass
            phase=0.0,
            sini=0.999,       # nearly edge-on
            Mstar=1.0,
            Rstar=1.0,
            ld_type="quad",
            ld_coeffs=[0.3, 0.3],
        ),
        cmd.BLS(minper=1.0, maxper=10.0, rmin=0.01, rmax=0.1),
    ])

    result = injection_pipe.run(lc)

    recovered_period = float(result.stats["BLS_Period_1_0"])
    sde = float(result.stats["BLS_SDE_1_0"])

    # Check period recovery
    injected_period = 3.14159
    recovered = abs(recovered_period - injected_period) / injected_period < 0.01
    print(f"Injected period: {injected_period:.5f} d")
    print(f"Recovered period:{recovered_period:.5f} d")
    print(f"SDE: {sde:.2f},  recovered: {recovered}")
    ```

For a full injection-recovery grid over period and radius ratio:

```python
import pandas as pd
import itertools

periods = np.linspace(1.0, 10.0, 20)
rp_values = [0.05, 0.08, 0.10, 0.12]
records = []

for period, rp in itertools.product(periods, rp_values):
    pipe = vt.Pipeline([
        cmd.Injecttransit(period=period, Rp=rp, Mp=0.001, phase=0.3,
                          sini=0.999, Mstar=1.0, Rstar=1.0),
        cmd.BLS(minper=0.8*period, maxper=1.2*period,
                rmin=0.005, rmax=0.2),
    ])
    r = pipe.run(lc)
    rec_period = float(r.stats["BLS_Period_1_0"])
    records.append({
        "injected_period": period,
        "rp": rp,
        "recovered_period": rec_period,
        "sde": float(r.stats["BLS_SDE_1_0"]),
        "recovered": abs(rec_period - period) / period < 0.02,
    })

recovery = pd.DataFrame(records)
recovery.to_csv("injection_recovery.csv", index=False)
```

---

## 4. TFA + BLS pipeline

Detrend with TFA before running BLS to reduce the number of systematic false positives.

=== "CLI"

    ```bash
    vartools -i target.lc \
        -savelc \
        -TFA trendlist.txt dates.txt 0 1 0 0 \
        -BLS r 0.01 0.1 0.5 20.0 optimal 1.0 200 0 1 operiodogram /output 0 1 0 \
        -oneline
    ```

=== "Python"

    ```python
    pipe = vt.Pipeline([
        cmd.savelc(),
        cmd.TFA(
            trendlist="trendlist.txt",
            dates_file="dates.txt",
            pixelsep=0,
            correct_lc=True,
        ),
        cmd.BLS(
            minper=0.5,
            maxper=20.0,
            rmin=0.01,
            rmax=0.1,
            nbins=200,
            save_periodogram=True,
        ),
    ])
    result = pipe.run(lc)

    print(f"BLS period: {float(result.stats['BLS_Period_1_0']):.5f} d")
    print(f"SDE:        {float(result.stats['BLS_SDE_1_0']):.2f}")
    ```

    The `savelc` step is not strictly required here (TFA subtraction does not destroy information), but is useful if you want to `restorelc` later to fit a transit model on the detrended-but-otherwise-unmodified light curve.

---

## 5. Interpreting BLS output

| Column | Meaning |
|--------|---------|
| `BLS_Period_1_N` | Best-fit period (days) |
| `BLS_SDE_1_N` | Signal Detection Efficiency — the BLS peak power normalised by the scatter in the power spectrum. Values above ~7 are worth inspecting. |
| `BLS_SN_1_N` | Signal-to-noise of the transit depth, accounting for the number of in-transit points. |
| `BLS_Depth_1_N` | Best-fit box depth (magnitudes; positive = dimming) |
| `BLS_Qtran_1_N` | Fractional transit duration (transit duration / period) |
| `BLS_Tc_1_N` | Transit epoch (same time system as the input light curve) |
| `BLS_i1_1_N`, `BLS_i2_1_N` | Phase-bin indices of ingress and egress |
| `BLS_PinkSN_1_N` | SDE computed accounting for correlated (pink) noise — more conservative than `BLS_SN`. |

!!! warning "False positives"
    High SDE alone does not confirm a transit. Common false positives include:

    - Eclipsing binaries at half the reported period
    - Systematic trends that survive sigma clipping
    - Aliased periods from the observing cadence

    Always inspect the phase-folded light curve, check for secondary eclipses at phase 0.5, and compare the `BLS_Depth` with the photometric noise level.
