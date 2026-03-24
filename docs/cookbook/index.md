# Cookbook

The cookbook collects practical, self-contained recipes for the most common VARTOOLS tasks. Every recipe shows both a CLI and a Python (`pyvartools`) approach side by side so that you can use whichever fits your workflow.

---

## Recipes

### [Period Finding](period-finding.md)

How to run Lomb-Scargle, BLS, and AoV period searches; interpret false-alarm probabilities and signal statistics; save full periodogram files; and chain sigma clipping with period finding.

- Basic Lomb-Scargle search
- LS with multiple peaks and SNR threshold
- Save the periodogram to file
- BLS transit search
- Analysis of Variance (AoV)
- Multi-harmonic AoV for non-sinusoidal signals
- Using `columnsuffix` for readable output names
- Chaining clip and period search

---

### [Transit Detection and Fitting](transit.md)

Full BLS transit search examples, interpreting BLS output statistics, injecting synthetic transits, and fitting the Mandel-Agol transit model.

- BLS transit search with parameter interpretation
- Fitting a Mandel-Agol transit model
- Injecting a transit signal for recovery tests
- TFA + BLS pipeline (detrend first, then search)
- Interpreting BLS output: SDE, S/N, depth, Qtran, pink noise

---

### [Detrending and Filtering](detrending.md)

Recipes for removing systematics and trends before period analysis, including sigma clipping, TFA, SYSREM, decorrelation, and median filtering.

- Sigma clipping
- TFA (Trend Filtering Algorithm)
- SYSREM systematic removal
- Polynomial decorrelation with `-decorr`
- Median filtering
- Combining detrending and period search with `savelc`/`restorelc`

---

### [Batch Processing](batch-processing.md)

How to efficiently process large numbers of light curves, collect results as DataFrames, handle errors, and use parallel threads.

- Process a list of files from disk
- Batch with parallel threads
- Get statistics as a DataFrame and save to CSV
- Capture modified light curves from a batch run
- Error handling with `raise_on_error=False`
- CLI parallel processing with `-parallel N`
- Survey-scale processing with a large list file

---

### [Working with FITS Files](fits-workflow.md)

How to read FITS light curves from the CLI and Python, map non-standard column names, convert flux to magnitude, and work with Kepler/TESS data.

- Reading a FITS light curve with the CLI
- Reading a FITS light curve with Python
- Writing output as FITS
- Kepler/TESS FITS workflow (PDCSAP flux to magnitude and LS search)
- Specifying FITS columns: `t_col`, `mag_col`, `err_col`, `hdu`
