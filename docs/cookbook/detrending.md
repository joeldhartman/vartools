# Detrending and Filtering

These recipes show how to remove trends, systematics, and outliers from light curves before or alongside period analysis.

---

## 1. Sigma clipping

Remove outlier points iteratively before any further analysis.

=== "CLI"

    ```bash
    vartools -i star.lc -clip 5.0 1 -rms -oneline
    ```

    `-clip 5.0 1` clips iteratively at 5σ until no more points are removed. Replace `1` with `0` for a single pass.

    To clip at most 3 times:

    ```bash
    vartools -i star.lc -clip 5.0 1 niter 3 -rms -oneline
    ```

=== "Python"

    ```python
    import pyvartools as vt
    from pyvartools import commands as cmd

    lc = vt.LightCurve.from_file("star.lc")

    # Iterative 5-sigma clip, then compute RMS
    pipe = vt.Pipeline([
        cmd.clip(sigclip=5.0, iterative=True),
        cmd.rms(),
    ])
    result = pipe.run(lc, capture_lc=True)

    print(f"RMS after clipping: {float(result.stats['RMS_0']):.5f} mag")
    print(f"Points retained:    {len(result.lc)}")
    ```

    Clip relative to the median (more robust for skewed distributions):

    ```python
    pipe = vt.Pipeline([cmd.clip(sigclip=4.0, median=True)])
    ```

---

## 2. TFA (Trend Filtering Algorithm)

TFA (Kovacs et al. 2005) detrends a target light curve by fitting a linear combination of template (trend) light curves from comparison stars observed in the same field.

**Required inputs:**

- `trendlist` — a text file listing the paths to the template light curves, one per line
- `dates_file` — a text file listing observation epochs (one BJD per line, matching the cadence of all light curves)

=== "CLI"

    ```bash
    vartools -i target.lc \
        -TFA trendlist.txt dates.txt 0 1 0 0 \
        -oneline
    ```

    Argument order after `-TFA trendlist dates_file pixelsep correct_lc save_coeffs save_model`.
    Here: pixel separation threshold = 0 (use all trend stars), correct LC = 1, no output files.

    With pixel-coordinate filtering (trend stars within 100 pixels):

    ```bash
    vartools -i target.lc \
        -TFA trendlist.txt dates.txt 100 xycol 3 4 1 0 0 \
        -oneline
    ```

=== "Python"

    ```python
    pipe = vt.Pipeline([
        cmd.TFA(
            trendlist="trendlist.txt",
            dates_file="dates.txt",
            pixelsep=0,          # 0 = use all trend stars
            correct_lc=True,
            save_model=True,     # capture the trend model
        ),
    ])
    result = pipe.run(lc)
    trend_model = result.files["TFA_model_0"]   # pd.DataFrame
    ```

    With pixel coordinates (columns 3 and 4 of the trend list hold x, y):

    ```python
    pipe = vt.Pipeline([
        cmd.TFA(
            trendlist="trendlist.txt",
            dates_file="dates.txt",
            pixelsep=100,
            xycol=(3, 4),
            correct_lc=True,
        ),
    ])
    ```

!!! note "Building the trend list and dates file"
    The trend list is simply a plain-text file with one path per line, pointing to comparison star light curves in the same format as the target. The dates file must have one line per observation epoch in the same order as the data columns. Both files must cover exactly the same set of observations as the target.

---

## 3. SYSREM

SYSREM (Tamuz et al. 2005) models and removes systematics that are correlated with colour and airmass across an ensemble of light curves.

=== "CLI"

    ```bash
    vartools -i target.lc \
        -SYSREM 1 1 airmass.dat 5.0 5.0 1e9 1 0 0 1 \
        -oneline
    ```

    Argument order: `ninput_color ninput_airmass airmass_file sigma_clip1 sigma_clip2 saturation correct_lc save_model save_trends useweights`.

=== "Python"

    ```python
    pipe = vt.Pipeline([
        cmd.SYSREM(
            ninput_color=1,
            ninput_airmass=1,
            initial_airmass_file="airmass.dat",
            sigma_clip1=5.0,
            sigma_clip2=5.0,
            correct_lc=True,
            save_model=True,
        ),
    ])
    result = pipe.run(lc)
    sysrem_model = result.files["SYSREM_model_0"]
    ```

---

## 4. Decorrelation

Polynomial decorrelation (`-decorr`) removes trends by fitting a polynomial function of one or more external vectors (e.g. airmass, CCD position, seeing) or light-curve-internal columns.

=== "CLI"

    ```bash
    # Decorrelate against an external airmass vector (1st-order polynomial)
    vartools -i star.lc \
        -decorr 1 1 0 1 airmass.dat 1 0 0 \
        -oneline
    ```

    Argument order after `-decorr`: `correct_lc zeropointterm subtractfirstterm nglobal [file order]... nlcterms [col order]... [omodel outdir | 0]`.

=== "Python"

    ```python
    # Decorrelate against an external airmass file (polynomial degree 1)
    pipe = vt.Pipeline([
        cmd.decorr(
            correct_lc=True,
            zeropointterm=1,
            subtractfirstterm=0,
            global_files=[("airmass.dat", 1)],
            save_model=True,
        ),
    ])
    result = pipe.run(lc)
    decorr_model = result.files["decorr_model_0"]
    ```

    Decorrelate against internal LC columns (e.g. x-position in column 4,
    y-position in column 5). If those extra columns were added via `aux` when
    constructing the `LightCurve`, pyvartools registers them with vartools
    automatically via `-inputlcformat`:

    ```python
    import numpy as np
    import pyvartools as vt
    from pyvartools import commands as cmd

    lc = vt.LightCurve.from_arrays(
        t, mag, err,
        aux={"xpos": xarr, "ypos": yarr},
    )

    pipe = vt.Pipeline([
        cmd.decorr(
            correct_lc=True,
            lc_columns=[(4, 1), (5, 1)],   # (column_number, polynomial_order)
        ),
    ])
    result = pipe.run(lc, capture_lc=True)
    ```

    For files read directly from disk, supply the layout via `columns`:

    ```python
    result = pipe.run_file(
        "star.lc",
        columns={"t": 1, "mag": 2, "err": 3, "xpos": 4, "ypos": 5},
    )
    ```

---

## 5. Median filtering

Smooth a light curve by replacing each point with the median of its neighbors in a sliding time window. Useful for estimating and removing slowly varying baselines.

=== "CLI"

    ```bash
    # Median filter with a 1-day window
    vartools -i star.lc -medianfilter 1.0 -rms -oneline
    ```

    To replace the magnitude values with the smoothed version (rather than storing the smoothed LC separately):

    ```bash
    vartools -i star.lc -medianfilter 1.0 replace -rms -oneline
    ```

=== "Python"

    ```python
    # Apply a 0.5-day median filter; capture the filtered LC
    pipe = vt.Pipeline([
        cmd.medianfilter(time=0.5, method="median", replace=True),
        cmd.rms(),
    ])
    result = pipe.run(lc, capture_lc=True)

    filtered_lc = result.lc
    print(f"RMS of filtered LC: {float(result.stats['RMS_0']):.5f} mag")
    ```

    Using a weighted mean filter instead:

    ```python
    pipe = vt.Pipeline([
        cmd.medianfilter(time=0.5, method="weightedaverage"),
    ])
    ```

---

## 6. Combining detrending with period search — `savelc`/`restorelc` pattern

A common pattern is to detrend the light curve, search for periods in the cleaned data, then restore the original to measure the signal amplitude or do a separate analysis.

=== "CLI"

    ```bash
    vartools -i target.lc \
        -savelc \
        -TFA trendlist.txt dates.txt 0 1 0 0 \
        -LS 0.5 10.0 0.001 1 0 \
        -restorelc 1 \
        -Killharm ls 3 0 0 \
        -rms \
        -oneline
    ```

    Flow:
    1. Save the raw LC
    2. TFA detrend
    3. LS search on the detrended LC
    4. Restore the raw LC
    5. Fit the LS period harmonically to the raw LC
    6. Compute RMS of the residuals

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
        cmd.LS(0.5, 10.0, 1e-3, npeaks=1),
        cmd.restorelc(savenumber=1),
        cmd.Killharm(period="ls", nharm=3),
        cmd.rms(),
    ])
    result = pipe.run(lc)

    period     = float(result.stats["LS_Period_1_0"])
    rms_resid  = float(result.stats["RMS_0"])
    print(f"Period:          {period:.5f} d")
    print(f"Residual RMS:    {rms_resid:.5f} mag")
    ```

For the TFA + BLS pattern used in transit searches, see [Transit Detection and Fitting](transit.md#4-tfa--bls-pipeline).
