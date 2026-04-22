# Command Index

Alphabetical index of every VARTOOLS command, with links to the
CLI reference and (where applicable) the pyvartools Python-API
reference.

The CLI column links to the command's section in the appropriate
`cli/*.md` page; the Python column links to its entry in
[`python/commands.md`](python/commands.md).  Commands flagged
*(USERLIB extension)* require the `-L` option to load the
corresponding library (see [Extension Commands](cli/extensions.md)).

!!! tip "Searching"
    The search box (top right) matches command names, descriptions,
    and doc content simultaneously — try it if you don't remember the
    exact spelling of a command.

| Command | Description | CLI | Python API |
|---|---|---|---|
| `-addfitskeyword` | Add a custom keyword to the output-light-curve FITS header | [CLI](cli/control-flow.md#-addfitskeyword) | [Python](python/commands.md#addfitskeyword-add-a-fits-keyword) |
| `-addnoise` | Add simulated white / red / correlated noise to light curves | [CLI](cli/simulation.md#-addnoise) | [Python](python/commands.md#addnoise-add-synthetic-noise) |
| `-alarm` | Alarm statistic (Kovacs, Bakos & Noyes 2005) | [CLI](cli/statistics.md#-alarm) | [Python](python/commands.md#alarm-alarm-statistic) |
| `-aov` | Phase-binned Analysis of Variance period search | [CLI](cli/period-finding.md#-aov-phase-binned-analysis-of-variance) | [Python](python/commands.md#aov-analysis-of-variance) |
| `-aov_harm` | Multi-harmonic Analysis of Variance period search | [CLI](cli/period-finding.md#-aov_harm-multi-harmonic-analysis-of-variance) | [Python](python/commands.md#aov_harm-multi-harmonic-aov) |
| `-autocorrelation` | Discrete autocorrelation function | [CLI](cli/statistics.md#-autocorrelation) | [Python](python/commands.md#autocorrelation-autocorrelation-function) |
| `-binlc` | Bin the light curve by time | [CLI](cli/manipulation.md#-binlc) | [Python](python/commands.md#binlc-bin-in-time_1) |
| `-BLS` | Box-fitting Least Squares transit search | [CLI](cli/period-finding.md#-bls-box-fitting-least-squares) | [Python](python/commands.md#bls-box-least-squares-transit-search) |
| `-BLSFixDurTc` | BLS with fixed transit duration and epoch | [CLI](cli/period-finding.md#-blsfixdurtc-bls-with-fixed-transit-duration-and-epoch) | [Python](python/commands.md#blsfixdurtc-bls-with-fixed-duration-and-epoch-searching-for-period) |
| `-BLSFixPer` | BLS at a fixed period | [CLI](cli/period-finding.md#-blsfixper-bls-at-a-fixed-period) | [Python](python/commands.md#blsfixper-bls-with-fixed-period) |
| `-BLSFixPerDurTc` | BLS with fixed period, duration, and epoch | [CLI](cli/period-finding.md#-blsfixperdurtc-bls-with-fixed-period-duration-and-epoch) | [Python](python/commands.md#blsfixperdurtc-bls-with-fixed-period-duration-and-epoch) |
| `-changeerror` | Replace the formal per-point uncertainty | [CLI](cli/manipulation.md#-changeerror) | [Python](python/commands.md#changeerror-changevariable-column-operations) |
| `-changevariable` | Reassign the internal t / mag / err / id roles | [CLI](cli/control-flow.md#-changevariable) | [Python](python/commands.md#changeerror-changevariable-column-operations) |
| `-chi2` | Weighted χ² relative to the mean | [CLI](cli/statistics.md#-chi2) | [Python](python/commands.md#chi2-chi-squared-statistic) |
| `-chi2bin` | Binned weighted χ² for several bin sizes | [CLI](cli/statistics.md#-chi2bin) | [Python](python/commands.md#chi2bin-binned-chi-squared) |
| `-clip` | Sigma-clip outliers from the light curve | [CLI](cli/filtering.md#-clip) | [Python](python/commands.md#clip-sigma-clipping) |
| `-converttime` | Convert between JD / MJD / HJD / BJD time systems | [CLI](cli/manipulation.md#-converttime) | [Python](python/commands.md#converttime-time-system-conversion) |
| `-copylc` | Replicate the in-memory light curve N times | [CLI](cli/simulation.md#-copylc) | [Python](python/commands.md#savelc-restorelc-copylc-lc-state-management) |
| `-decorr` | Polynomial decorrelation against external signals (deprecated; use -linfit) | [CLI](cli/model-fitting.md#-decorr) | [Python](python/commands.md#decorr-decorrelation) |
| `-dftclean` | Discrete Fourier Transform + CLEAN power spectrum | [CLI](cli/period-finding.md#-dftclean-dft-power-spectrum-clean) | [Python](python/commands.md#dftclean-clean-periodogram) |
| `-difffluxtomag` | Convert ISIS differential flux to magnitude | [CLI](cli/manipulation.md#-difffluxtomag) | [Python](python/commands.md#difffluxtomag-fluxtomag-flux-conversions) |
| `-ensemblerescalesig` | Rescale all light curves' errors so ⟨χ²/dof⟩ = 1 | [CLI](cli/manipulation.md#-ensemblerescalesig) | [Python](python/commands.md#rescalesig-ensemblerescalesig-rescale-per-point-uncertainties) |
| `-expr` | Evaluate an analytic expression; create or update a variable | [CLI](cli/manipulation.md#-expr) | [Python](python/commands.md#expr-analytic-expression) |
| `-fastchi2` | Palmer's fast χ² periodogram (USERLIB extension) | [CLI](cli/extensions.md#-fastchi2) | — |
| `-FFT` | Fast Fourier transform (evenly-sampled data) | [CLI](cli/manipulation.md#-fft) | [Python](python/commands.md#fft-ifft-fast-fourier-transform) |
| `-findblends` | Identify nearby variables that may blend into the target | [CLI](cli/misc.md#-findblends) | [Python](python/commands.md#findblends-search-for-blended-transits) |
| `-fluxtomag` | Convert flux to magnitude | [CLI](cli/manipulation.md#-fluxtomag) | [Python](python/commands.md#difffluxtomag-fluxtomag-flux-conversions) |
| `-ftuneven` | Complex Fourier transform of unevenly sampled data (USERLIB extension) | [CLI](cli/extensions.md#-ftuneven) | — |
| `-GetLSAmpThresh` | Minimum detectable amplitude at a given period (LS-based) | [CLI](cli/period-finding.md#-getlsampthresh-minimum-detectable-amplitude) | [Python](python/commands.md#getlsampthresh-ls-amplitude-threshold) |
| `-hatpiflag` | Combine HATPI quality flags into a single binary flag (USERLIB extension) | [CLI](cli/extensions.md#-hatpiflag) | — |
| `-if` | Conditional block (-if / -elif / -else / -fi) | [CLI](cli/control-flow.md#-if-elif-else-fi) | [Python](python/commands.md#ifcmd-conditional-execution) |
| `-IFFT` | Inverse fast Fourier transform | [CLI](cli/manipulation.md#-ifft) | [Python](python/commands.md#fft-ifft-fast-fourier-transform) |
| `-Injectharm` | Inject a harmonic (sinusoidal) signal | [CLI](cli/simulation.md#-injectharm) | [Python](python/commands.md#injectharm-inject-a-harmonic-signal) |
| `-Injecttransit` | Inject a Mandel-Agol transit model | [CLI](cli/simulation.md#-injecttransit) | [Python](python/commands.md#injecttransit-inject-a-transit-signal) |
| `-jktebop` | JKTEBOP detached eclipsing binary model (USERLIB extension) | [CLI](cli/extensions.md#-jktebop) | — |
| `-Jstet` | Stetson J variability statistic | [CLI](cli/statistics.md#-jstet) | [Python](python/commands.md#jstet-stetson-j-statistic) |
| `-Killharm` | Fit and subtract a harmonic (Fourier) model | [CLI](cli/model-fitting.md#-killharm) | [Python](python/commands.md#killharm-harmonic-series-subtraction) |
| `-linfit` | Linear regression against a user-specified basis | [CLI](cli/model-fitting.md#-linfit) | [Python](python/commands.md#linfit-linear-combination-fitting) |
| `-LS` | Generalized Lomb-Scargle periodogram | [CLI](cli/period-finding.md#-ls-generalized-lomb-scargle) | [Python](python/commands.md#ls-generalized-lomb-scargle) |
| `-macula` | Macula rotation + spot light-curve model (USERLIB extension) | [CLI](cli/extensions.md#-macula) | — |
| `-magadd` | Add a constant offset to the magnitudes (USERLIB extension) | [CLI](cli/extensions.md#-magadd) | — |
| `-MandelAgolTransit` | Mandel & Agol analytic transit model | [CLI](cli/model-fitting.md#-mandelagoltransit) | [Python](python/commands.md#mandelagoltransit-mandel-agol-transit-model) |
| `-match` | Match observations across light curves by time or string ID | [CLI](cli/manipulation.md#-match) | [Python](python/commands.md#match-match-against-a-catalog) |
| `-medianfilter` | Running-median high-pass or low-pass filter | [CLI](cli/filtering.md#-medianfilter) | [Python](python/commands.md#medianfilter-median-filtering) |
| `-microlens` | Point-source gravitational microlensing model | [CLI](cli/model-fitting.md#-microlens) | [Python](python/commands.md#microlens-microlensing-model) |
| `-nonlinfit` | Levenberg-Marquardt / MCMC nonlinear-model fitting | [CLI](cli/model-fitting.md#-nonlinfit) | [Python](python/commands.md#nonlinfit-non-linear-least-squares-fitting) |
| `-o` | Write the (possibly modified) light curve to a file | [CLI](cli/control-flow.md#-o) | [Python](python/commands.md#o-output-light-curve) |
| `-Phase` | Replace the time column with phase-folded coordinates | [CLI](cli/manipulation.md#-phase) | [Python](python/commands.md#phase-phase-fold-the-light-curve) |
| `-print` | Emit user-computed scalars as columns in the statistics table | [CLI](cli/control-flow.md#-print) | [Python](python/commands.md#print-emit-user-computed-variables-to-the-output-table) |
| `-python` | Run embedded Python code on each light curve | [CLI](cli/python-r.md#-python) | — |
| `-R` | Run embedded R code on each light curve | [CLI](cli/python-r.md#-r) | [Python](python/commands.md#r-run-r-code) |
| `-resample` | Interpolate the light curve onto a new time grid | [CLI](cli/manipulation.md#-resample) | [Python](python/commands.md#resample-resample-onto-a-new-time-grid) |
| `-rescalesig` | Rescale per-point errors to make χ²/dof = 1 | [CLI](cli/manipulation.md#-rescalesig) | [Python](python/commands.md#rescalesig-ensemblerescalesig-rescale-per-point-uncertainties) |
| `-restorelc` | Restore a light-curve state saved by -savelc | [CLI](cli/control-flow.md#-savelc-restorelc) | [Python](python/commands.md#savelc-restorelc-copylc-lc-state-management) |
| `-restoretimes` | Undo a prior -restricttimes filter | [CLI](cli/filtering.md#-restoretimes) | [Python](python/commands.md#restricttimes-restoretimes-time-windowing) |
| `-restricttimes` | Filter observations by JD range, JD list, image ID, or expression | [CLI](cli/filtering.md#-restricttimes) | [Python](python/commands.md#restricttimes-restoretimes-time-windowing) |
| `-rms` | Mean magnitude, RMS, and expected RMS | [CLI](cli/statistics.md#-rms) | [Python](python/commands.md#rms-root-mean-square) |
| `-rmsbin` | Binned RMS for several bin sizes | [CLI](cli/statistics.md#-rmsbin) | [Python](python/commands.md#rmsbin-binned-rms) |
| `-savelc` | Checkpoint the in-memory light curve for later -restorelc | [CLI](cli/control-flow.md#-savelc-restorelc) | [Python](python/commands.md#savelc-restorelc-copylc-lc-state-management) |
| `-SoftenedTransit` | Protopapas softened-transit empirical model | [CLI](cli/model-fitting.md#-softenedtransit) | [Python](python/commands.md#softenedtransit-softened-trapezoidal-transit) |
| `-sortlc` | Sort the light curve by time (or another variable) | [CLI](cli/manipulation.md#-sortlc) | [Python](python/commands.md#sortlc-sort-observations) |
| `-splinedetrend` | Spline-based detrending (USERLIB extension) | [CLI](cli/extensions.md#-splinedetrend) | — |
| `-Starspot` | Star-spot rotational modulation model | [CLI](cli/model-fitting.md#-starspot) | [Python](python/commands.md#starspot-starspot-model) |
| `-stats` | Compute arbitrary summary statistics on any variable | [CLI](cli/statistics.md#-stats) | [Python](python/commands.md#stats-generic-statistics) |
| `-stitch` | Stitch multi-segment light curves at offsets (USERLIB extension) | [CLI](cli/extensions.md#-stitch) | — |
| `-SYSREM` | Systematic-error removal (Tamuz, Mazeh & Zucker 2005) | [CLI](cli/filtering.md#-sysrem) | [Python](python/commands.md#sysrem-systematic-noise-removal) |
| `-TFA` | Trend Filtering Algorithm (Kovacs, Bakos & Noyes 2005) | [CLI](cli/filtering.md#-tfa) | [Python](python/commands.md#tfa-trend-filtering-algorithm) |
| `-TFA_SR` | TFA with signal-reconstruction | [CLI](cli/filtering.md#-tfa_sr) | [Python](python/commands.md#tfa_sr-tfa-with-signal-reconstruction) |
| `-wwz` | Weighted Wavelet Z-transform (Foster 1996) | [CLI](cli/period-finding.md#-wwz-weighted-wavelet-z-transform) | [Python](python/commands.md#wwz-weighted-wavelet-z-transform) |

## Python helper classes

These are not CLI commands but are pyvartools-side helpers that
appear in the Python API docs:

| Class / helper | Description |
|---|---|
| [`Output`](python/commands.md#output-class) | Auxiliary-output-file control passed to `save_*=...` parameters |

## Commands without a dedicated Python wrapper

The 9 commands in the first table that show `—` in
the Python column do not have a typed pyvartools class.  Both of
the groups they fall into are still reachable from Python — the
second group just takes a bit more manual setup than a one-line
`cmd.X(...)` call.

### USERLIB extension commands

`-fastchi2`, `-ftuneven`, `-hatpiflag`, `-jktebop`, `-macula`,
`-magadd`, `-splinedetrend`, and `-stitch` are all shipped as
dynamic-library extensions in `USERLIBS/src/`.  pyvartools provides
a generic [`UserCommand`](python/commands.md#usage-pattern-1-quick-one-off-usercommand)
that wraps any such extension without needing a dedicated subclass:

```python
from pyvartools import UserCommand, Pipeline
pipe = Pipeline([
    UserCommand("USERLIBS/src/.libs/fastchi2.so", "fastchi2",
                "Nharm fix 2 freqmax fix 24.0"),
])
result = pipe.run_file("EXAMPLES/2")
```

See the [User Extension Commands](python/commands.md#user-extension-commands)
section for three additional usage patterns: `load_userlib` (a named
class bound at import time), `discover_userlibs` (bulk
auto-registration), and the full subclass pattern for hand-written
typed wrappers.  Any of these would turn an entry in the table above
into a first-class `cmd.X(...)` call; a typed-wrapper subclass is
typically ~30 lines per extension.

### `-python`

`-python` evaluates arbitrary user Python code *inside* a VARTOOLS
run via an embedded interpreter.  A dedicated wrapper would be a bit
unusual, because the more idiomatic thing to do when you are already
in Python is to write a function that operates on a `LightCurve` and
call it directly (or wrap it in a custom `VartoolsCommand` subclass
that doesn't round-trip through a `-python` interpreter).  A
`cmd.python(code, invars=..., outvars=..., outputcolumns=...)` class
would be straightforward to add (~20 lines, matching the `-R`
counterpart that pyvartools already has) if there is demand.
