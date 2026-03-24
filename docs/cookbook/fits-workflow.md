# Working with FITS Files

VARTOOLS can read and write FITS binary tables when compiled with cfitsio. These recipes cover reading FITS light curves from the CLI and Python, mapping non-standard column names, converting flux to magnitude, and working with Kepler and TESS data.

---

## 1. Reading a FITS light curve with the CLI

Use `-inputlcformat` to tell vartools which FITS columns to use as time, magnitude, and error.

=== "CLI"

    ```bash
    # Read a FITS LC with standard column names (BJD, Mag, Err) in HDU 1
    vartools -i star.fits \
        -inputlcformat fits bjd:BJD mag:Mag err:Err \
        -rms -oneline
    ```

    ```bash
    # Non-standard column names
    vartools -i kepler.fits \
        -inputlcformat fits bjd:TIME mag:PDCSAP_FLUX err:PDCSAP_FLUX_ERR \
        -rms -oneline
    ```

    If the light curve is in a different HDU:

    ```bash
    vartools -i star.fits \
        -inputlcformat fits hdu:2 bjd:BJD mag:Mag err:Err \
        -rms -oneline
    ```

    To rename FITS columns for use by subsequent commands, combine with `-changevariable`:

    ```bash
    vartools -i star.fits \
        -inputlcformat fits bjd:TIME mag:FLUX err:FLUX_ERR \
        -changevariable mag FLUX \
        -changevariable err FLUX_ERR \
        -rms -oneline
    ```

---

## 2. Reading a FITS light curve with Python

`LightCurve.from_file` detects FITS files automatically from the `.fits`, `.fit`, or `.fts` extension. Use `t_col`, `mag_col`, and `err_col` to specify which FITS columns map to time, magnitude, and uncertainty.

=== "Python"

    ```python
    import pyvartools as vt

    # Standard column names (defaults)
    lc = vt.LightCurve.from_file("star.fits")

    # Custom column names
    lc = vt.LightCurve.from_file(
        "star.fits",
        t_col="BJD",
        mag_col="Mag",
        err_col="Err",
    )

    # Non-default HDU
    lc = vt.LightCurve.from_file(
        "star.fits",
        t_col="BJD",
        mag_col="Mag",
        err_col="Err",
        hdu=2,
    )

    print(f"Name:   {lc.name}")
    print(f"Points: {len(lc)}")
    print(f"Columns:{list(lc.to_dataframe().columns)}")
    ```

    Any FITS columns not mapped to `t`, `mag`, or `err` are preserved as auxiliary columns in the `LightCurve` DataFrame and can be accessed with `lc.to_dataframe()`.

---

## 3. Writing output as FITS

=== "CLI"

    Use the `-of` flag (output format) together with `-o` to write the processed light curve as a FITS file:

    ```bash
    vartools -i star.fits \
        -inputlcformat fits bjd:BJD mag:Mag err:Err \
        -clip 5.0 1 \
        -LS 0.5 10.0 0.001 1 0 \
        -of fits \
        -o output.fits \
        -oneline
    ```

=== "Python"

    The Python wrapper handles output format through the pipeline's `outdir` parameter and the `cmd.o` command. For FITS output, use `Raw` to pass the `-of fits` flag directly:

    ```python
    from pyvartools import commands as cmd
    import pyvartools as vt

    pipe = vt.Pipeline([
        cmd.clip(5.0),
        cmd.rms(),
        cmd.Raw("-of fits"),
        cmd.o("output.fits"),
    ])
    result = pipe.run(lc)
    ```

---

## 4. Kepler/TESS FITS workflow

Kepler and TESS deliver light curves as FITS binary tables. The standard columns are `TIME` (BKJD or BTJD, i.e. BJD minus a constant), `PDCSAP_FLUX` (pre-search data conditioning flux), and `PDCSAP_FLUX_ERR`. These need to be converted to physical time and magnitude units before period searching.

=== "CLI"

    ```bash
    # Kepler LC: TIME is BKJD = BJD - 2454833.0
    # Convert differential flux to magnitude, then run LS

    vartools -i kplr012345678_llc.fits \
        -inputlcformat fits bjd:TIME mag:PDCSAP_FLUX err:PDCSAP_FLUX_ERR \
        -expr "t=t+2454833.0" \
        -difffluxtomag 12.0 0 \
        -clip 5.0 1 \
        -LS 0.5 30.0 0.001 1 0 \
        -oneline
    ```

    For TESS (`TIME` is BTJD = BJD - 2457000.0):

    ```bash
    vartools -i tess_lc.fits \
        -inputlcformat fits bjd:TIME mag:PDCSAP_FLUX err:PDCSAP_FLUX_ERR \
        -expr "t=t+2457000.0" \
        -difffluxtomag 10.0 0 \
        -clip 5.0 1 \
        -LS 0.1 13.5 0.001 1 0 \
        -oneline
    ```

=== "Python"

    ```python
    import pyvartools as vt
    from pyvartools import commands as cmd

    # Kepler long-cadence FITS
    lc = vt.LightCurve.from_file(
        "kplr012345678_llc.fits",
        t_col="TIME",
        mag_col="PDCSAP_FLUX",
        err_col="PDCSAP_FLUX_ERR",
    )

    # vartools expects BJD; Kepler TIME is BKJD = BJD - 2454833.0
    # Use expr to shift the time column, then convert flux to magnitude.
    pipe = vt.Pipeline([
        cmd.expr("t=t+2454833.0"),
        cmd.difffluxtomag(mag_constant=12.0, offset=0.0),
        cmd.clip(5.0),
        cmd.LS(0.5, 30.0, 1e-3, npeaks=1, save_periodogram=True),
    ])
    result = pipe.run(lc)

    print(f"LS period:  {float(result.stats['LS_Period_1_0']):.5f} d")
    print(f"log10(FAP): {float(result.stats['Log10_LS_Prob_1_0']):.2f}")

    pgram = result.files["LS_periodogram_3"]
    ```

    For TESS (BTJD = BJD − 2457000.0):

    ```python
    lc = vt.LightCurve.from_file(
        "tess_sector01_lc.fits",
        t_col="TIME",
        mag_col="PDCSAP_FLUX",
        err_col="PDCSAP_FLUX_ERR",
    )

    pipe = vt.Pipeline([
        cmd.expr("t=t+2457000.0"),
        cmd.difffluxtomag(mag_constant=10.0, offset=0.0),
        cmd.clip(5.0),
        cmd.BLS(
            minper=0.5,
            maxper=13.5,   # TESS sectors are ~27 d; reliable up to ~half that
            rmin=0.005,
            rmax=0.1,
            nbins=200,
        ),
    ])
    result = pipe.run(lc)

    print(f"BLS period: {float(result.stats['BLS_Period_1_0']):.5f} d")
    print(f"SDE:        {float(result.stats['BLS_SDE_1_0']):.2f}")
    ```

!!! tip "NaN handling"
    TESS and Kepler FITS files contain `NaN` values for cadences affected by momentum dumps, cosmic rays, or thruster fires. vartools ignores NaN values in the time and magnitude columns automatically — no pre-filtering is needed.

---

## 5. Specifying FITS columns

Summary of the `LightCurve.from_file` parameters for FITS input:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `t_col` | `"BJD"` | FITS column to use as the time axis. Case-insensitive match. |
| `mag_col` | `"Mag"` | FITS column to use as the magnitude/flux axis. |
| `err_col` | `"Err"` | FITS column to use as the uncertainty. |
| `hdu` | `1` | HDU index (integer, 0-based from the top of the file). HDU 0 is the primary header; data are usually in HDU 1. |
| `format` | auto | `"fits"` or `"ascii"`. Auto-detected from the file extension if omitted. |

Column matching is case-insensitive, so `t_col="bjd"` and `t_col="BJD"` both work.

If a requested column is not found, `LightCurve.from_file` raises a `ValueError` listing the available column names, which helps diagnose mismatches:

```python
# This will raise a helpful error showing available columns
try:
    lc = vt.LightCurve.from_file("star.fits", t_col="WRONG_NAME")
except ValueError as e:
    print(e)
# Column 'WRONG_NAME' not found in FITS HDU 1 of 'star.fits'.
# Available: ['TIME', 'TIMECORR', 'CADENCENO', 'SAP_FLUX', 'PDCSAP_FLUX', ...]
```

To inspect column names before loading:

```python
from astropy.io import fits

with fits.open("star.fits") as hdul:
    for hdu in hdul:
        if hasattr(hdu, "columns"):
            print(f"HDU {hdul.index_of(hdu.name)}: {[c.name for c in hdu.columns]}")
```
