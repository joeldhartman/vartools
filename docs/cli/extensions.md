# Extension Commands

Extension commands are additional VARTOOLS commands that are compiled as
separate shared-object libraries (`.so` files) and loaded at runtime. They are
distributed in the `USERLIB/src` subdirectory of the VARTOOLS source tree.
After running `make` in that directory, each extension produces a `.so` file
that must be passed to VARTOOLS with the `-L` option before the command itself
is issued:

```bash
vartools -L USERLIB/src/fastchi2.so \
    -l EXAMPLES/lc_list \
    -fastchi2 "Nharm" "fix" 2 "freqmax" "fix" 24.0 -tab
```

The same `-L` flag can be given multiple times to load several extensions in a
single call.

---

## `-fastchi2`

**Palmer's Fast χ² periodogram** (Palmer 2009, ApJ, 695, 496).

```
-fastchi2
    <"Nharm" <"fix" val | "list" | "fixcolumn" <colname | colnum>>
    <"freqmax" <"fix" val | "list" | "fixcolumn" <colname | colnum>>
    ["freqmin" <"fix" val | "list" | "fixcolumn" <colname | colnum>>]
    ["detrendorder" <"fix" val | "list" | "fixcolumn" <colname | colnum>>]
    ["t0" <"fix" val | "list" | "fixcolumn" <colname | colnum>>]
    ["timespan" <"fix" val | "list" | "fixcolumn" <colname | colnum>>]
    ["oversample" <"fix" val | "list" | "fixcolumn" <colname | colnum>>]
    ["chimargin" <"fix" val | "list" | "fixcolumn" <colname | colnum>>]
    ["Npeak" val]
    ["norefitpeak"]
    ["oper" outdir ["nameformat" format]]
    ["omodel" outdir ["nameformat" format]]
    ["omodelvariable" varname]
```

Computes the Fast χ² periodogram using Palmer's algorithm, which searches for
the best-fitting multi-harmonic sinusoidal model at each trial frequency. This
implementation uses Palmer's original code.

For each parameter, specify the source with one of three keywords: `"fix"` to
set a constant for all light curves, `"list"` to read the value from the input
list, or `"fixcolumn"` to take it from the output of a previous command.

**Required parameters:**

| Parameter | Description |
|---|---|
| `Nharm` | Number of harmonics to include (1 = fundamental only, 2 = fundamental + first overtone, …). |
| `freqmax` | Maximum frequency to search, in cycles per day. |

**Optional parameters:**

| Parameter | Default | Description |
|---|---|---|
| `freqmin` | 0 | Minimum frequency to search. |
| `detrendorder` | 0 | Order of polynomial detrend applied before the search. |
| `t0` | — | Reference epoch for detrending. |
| `timespan` | — | Total time span, used for the Nyquist frequency. |
| `oversample` | — | Oversampling factor for the periodogram frequency grid. |
| `chimargin` | — | χ² margin for selecting periodogram peaks for fine search. |
| `Npeak` | — | Number of peaks to report. |

If `"norefitpeak"` is given, no refined search is performed around candidate
peaks; only the raw periodogram peaks are returned. Use `"oper"` to write the
periodogram to files, and `"omodel"` to write the best-fit harmonic model.

**Citation:** Palmer, D. M. 2009, ApJ, 695, 496.

---

## `-ftuneven`

**Complex Fourier transform of unevenly sampled data** (Scargle 1989, ApJ, 343, 874).

```
-ftuneven
    <"outputvectors" freq_vec FTreal_vec FTimag_vec Pgram_vec |
     "outputfile" outdir ["nameformat" fmt] |
     "outputvectorsandfile" outdir ["nameformat" fmt]
                     freq_vec FTreal_vec FTimag_vec Pgram_vec>
    <"freqauto" |
     "freqrange" "minfreq" <...> "maxfreq" <...> "freqstep" <...> |
     "freqvariable" varname |
     "freqfile" filename>
    ["ft_sign" val] ["tt_zero" val]
    ["changeinputvectors" tvec data_real_vec data_imag_vec]
```

Computes the complex Fourier transform of an unevenly sampled time series using
Scargle's method. Returns the real and imaginary components of the transform
and the absolute-square power spectrum (equivalent to the Lomb-Scargle
periodogram). **Frequencies are in radians per unit time.**

**Output mode** (choose one):

| Mode | Behavior |
|---|---|
| `"outputvectors"` | Store results in named light curve vectors for use by subsequent commands. All light curve vectors are resized to match the transform length. |
| `"outputfile"` | Write results to files in `outdir` (default name: `BASELC_NAME.ftuneven`). |
| `"outputvectorsandfile"` | Do both simultaneously. |

**Frequency specification** (choose one):

| Mode | Behavior |
|---|---|
| `"freqauto"` | Determine frequencies automatically from the time baseline and cadence. |
| `"freqrange"` | Uniform grid; supply `"minfreq"`, `"maxfreq"`, and `"freqstep"`. |
| `"freqvariable"` | Read frequencies from an existing light curve vector. |
| `"freqfile"` | Read frequencies from the first column of a whitespace-delimited ASCII file; used identically for every light curve. |

Optional parameters `"ft_sign"` (default −1, forward transform; use +1 for
inverse) and `"tt_zero"` (time origin, default 0) are also available. Use
`"changeinputvectors"` to apply the transform to vectors other than the default
`t` and `mag`.

**Citation:** Scargle, J. D. 1989, ApJ, 343, 874.

---

## `-hatpiflag`

**Binary flag generator for HATPI light curves.**

```
-hatpiflag
    fiphot_string_flag_var rejbadframe_mask_var TFA_outlier_mask_var
    pointing_outlier_flag_var output_flag_var
```

Combines several per-point quality indicators into a single binary flag
variable suitable for HATPI photometry pipelines.

| Argument | Description |
|---|---|
| `fiphot_string_flag_var` | Vector of string flags from fiphot. |
| `rejbadframe_mask_var` | Bad-frame mask (0 = rejected, 1 = keep). |
| `TFA_outlier_mask_var` | TFA outlier mask (0 = outlier, 1 = keep). |
| `pointing_outlier_flag_var` | Pointing outlier flag (1 = outlier, 0 = ok). |
| `output_flag_var` | Name of the variable to receive the combined binary flag. |

---

## `-jktebop`

**JKTEBOP detached eclipsing binary light curve model.**

```
-jktebop
    < "inject" | "fit" >
    <"Period"   <"fix" value | "list" | "fixcolumn" <colname|colnum>> ["vary"]>
    <"T0"       <"fix" value | "list" | "fixcolumn" <colname|colnum>> ["vary"]>
    <"r1+r2"    <"fix" value | "list" | "fixcolumn" <colname|colnum>> ["vary"]>
    <"r2/r1"    <"fix" value | "list" | "fixcolumn" <colname|colnum>> ["vary"]>
    <"M2/M1"    <"fix" value | "list" | "fixcolumn" <colname|colnum>> ["vary"]>
    <"J2/J1"    <"fix" value | "list" | "fixcolumn" <colname|colnum>> ["vary"]>
    <<"i" | "bimpact"> <"fix" value | ...> ["vary"]>
    <"esinomega" <...> ["vary"]>
    <"ecosomega" <...> ["vary"]>
    <"LD1" <"linear"|"quad"|"log"|"sqrt"> <"fix" value1 [value2] | ...> ["vary"]>
    <"LD2" <"lockLD1"|"linear"|"quad"|"log"|"sqrt"> [...] ["vary"]>
    ["gravdark1" <...> ["vary"]]  ["gravdark2" <...> ["vary"]]
    ["reflection1" <...> ["vary"]] ["reflection2" <...> ["vary"]]
    ["L3" <...> ["vary"]] ["tidallag" <...> ["vary"]]
    ["correctlc"]
    ["omodel" <outdir ["format" fmt]>]
    ["ocurve" <"jd" | "phase"> ["step" stepsize] <"outdir" outdir ["format" fmt]>]
```

Fits or injects a JKTEBOP detached eclipsing binary light curve model. Supply
`"inject"` to add the model to the light curve, or `"fit"` to optimize the
parameters to the data.

**Geometric parameters:**

| Parameter | Description |
|---|---|
| `Period` | Orbital period in days. |
| `T0` | Central time of a primary eclipse. |
| `r1+r2` | Sum of stellar radii divided by the semi-major axis. |
| `r2/r1` | Ratio of stellar radii. |
| `M2/M1` | Mass ratio. |
| `J2/J1` | Surface brightness ratio. |
| `i` | Orbital inclination in degrees (90° = edge-on). |
| `bimpact` | Impact parameter at primary eclipse (alternative to `i`; 0 = central, 1 = grazing). |
| `esinomega` | *e* sin ω. |
| `ecosomega` | *e* cos ω. |

**Limb darkening** (`LD1`, `LD2`): choose a law (`"linear"`, `"quad"`, `"log"`,
`"sqrt"`) and supply one or two coefficients as required. For `LD2` the keyword
`"lockLD1"` forces the secondary to share the primary's coefficients.

**Optional physical parameters:** `gravdark1/2` (gravity darkening, default 1.0),
`reflection1/2` (reflection coefficients, computed if absent or ≤ 0),
`L3` (third light, default 0), `tidallag` (tidal lag in degrees, default 0).

For each parameter, append `"vary"` to allow it to vary during fitting.

`"correctlc"` subtracts the best-fit model from the light curve. `"omodel"`
writes model files to `outdir`. `"ocurve"` writes a uniformly sampled model
curve in time (`"jd"`) or phase (`"phase"`).

**Citations:** Southworth et al. 2004, MNRAS, 351, 1277; Popper & Etzel 1981,
AJ, 86, 102; Etzel 1981, *Photometric and Spectroscopic Binary Systems*, 111;
Nelson & Davis 1972, ApJ, 174, 617.

---

## `-splinedetrend`

**Multivariate light curve detrending using linear basis functions.**

```
-splinedetrend
    detrendvec1:<"spline":knotspacing:order |
                 "poly":order |
                 "harm":nharm>[:"groupbygap":gapsize]
    [,detrendvec2:<...>[:"groupbygap":gapsize],...]
    ["sigmaclip" <"fix" val | "list" | "fixcolumn" <colname|colnum> | "expr" expr>]
    ["omodel" outdir ["nameformat" format]]
    ["omodelcoeffs" outdir ["nameformat" format]]
    ["omodelvariable" outvarname1[:inputvarsignal1][,...]]
```

Fits a multivariate linear model to the light curve magnitudes using one or
more auxiliary variables (e.g., time, CCD *x*/*y* position). Cross-terms
between variables are not included.

The first argument is a comma-separated list of `detrendvec:basistype` tokens.
Three basis types are supported:

| Basis | Parameters | Description |
|---|---|---|
| `"spline"` | `knotspacing`, `order` | B-spline basis (GSL `gsl_bspline_eval`). `knotspacing` in the same units as the variable; `order` = 3 for cubic. |
| `"poly"` | `order` | Polynomial; `order` = 1 for linear, 2 for quadratic, etc. |
| `"harm"` | `nharm` | Harmonic series; `nharm` = 0 for fundamental only (period = 2× variable range). |

Append `:"groupbygap":gapsize` to split the fit at gaps in the detrending
variable larger than `gapsize`.

**Optional keywords:**

| Keyword | Description |
|---|---|
| `"sigmaclip"` | Exclude outliers from the fit using the given sigma threshold; the model is still evaluated and subtracted at clipped points. |
| `"omodel"` | Write the best-fit model to files in `outdir` (default extension: `.splinedetrend_model`). |
| `"omodelcoeffs"` | Write the linear basis coefficients to files (default extension: `.splinedetrend_modelcoeffs`). |
| `"omodelvariable"` | Store the best-fit model contribution from a given input variable in a named output variable. |

---

## `-stitch`

**Fit and remove zero-point offsets between combined light curve segments.**

```
-stitch
    stitch_variable_list uncertainty_variable_list
    mask_variable_list lcnum_var
    ["refnum_var" refnum_var]
    <"median" | "mean" | "weightedmean" | "poly" order |
     "harmseries" period_var Nharm>
    ["groupbytime" time_bin ["start" firstbintime]]
    ["fitonly"]
    ["save_fitted_parameters" <outdir ["format" fmt]>]
    ...
```

Designed for use with the `-l "combinelcs"` input mode, `-stitch` fits for and
removes additive offsets between distinct light curve segments (e.g.,
observations from different telescopes or cameras).

**Required arguments:**

| Argument | Description |
|---|---|
| `stitch_variable_list` | Comma-separated list of magnitude variables to stitch (typically just `mag`). |
| `uncertainty_variable_list` | Comma-separated uncertainties for each magnitude variable (typically `err`). |
| `mask_variable_list` | Comma-separated mask vectors; points with value > 0 are excluded from the fit. |
| `lcnum_var` | Variable identifying which input segment each observation belongs to (e.g., set by `"lcnumvar"` in `combinelcs`). |

**Fitting methods** (choose one):

| Method | Description |
|---|---|
| `"median"` | Compute per-segment medians and solve for additive offsets. |
| `"mean"` | As above but using the mean. |
| `"weightedmean"` | Inverse-variance weighted mean. |
| `"poly" order` | Fit a polynomial in time to each segment. |
| `"harmseries" period_var Nharm` | Fit a harmonic series with period from `period_var` and `Nharm` harmonics. |

**Optional keywords:**

| Keyword | Description |
|---|---|
| `"refnum_var" varname` | Further subdivide segments by a second grouping variable. |
| `"groupbytime" time_bin` | Group segments into time bins; the bin size is automatically widened if necessary to ensure all segments can be inter-calibrated. |
| `"fitonly"` | Compute the shifts but do not subtract them from the light curve. |
| `"save_fitted_parameters" outdir` | Write per-source shift files to `outdir` (default suffix: `.stitch`). |
| `"shifts_file"` | Read previously determined shifts from files and/or write newly determined shifts out, for incremental re-processing of large datasets. |

---

## `-macula`

**Kipping's Macula analytic starspot model** (Kipping 2012, arXiv:1209.2985).

```
-macula
    <"inject" | "fit" <"amoeba" | "lm">>
    <"Prot"   <"fix" value | "list" | "fixcolumn" <colname|colnum>> ["vary"]>
    <"istar"  <...> ["vary"]>
    <"kappa2" <...> ["vary"]>  <"kappa4" <...> ["vary"]>
    <"c1" <...> ["vary"]>  <"c2" <...> ["vary"]>
    <"c3" <...> ["vary"]>  <"c4" <...> ["vary"]>
    <"d1" <...> ["vary"]>  <"d2" <...> ["vary"]>
    <"d3" <...> ["vary"]>  <"d4" <...> ["vary"]>
    <"blend" <...> ["vary"]>
    <"Nspot" value>
        <"Lambda0" <...> ["vary"]>  <"Phi0" <...> ["vary"]>
        <"alphamax" <...> ["vary"]> <"fspot" <...> ["vary"]>
        <"tmax" <...> ["vary"]>     <"life" <...> ["vary"]>
        <"ingress" <...> ["vary"]>  <"egress" <...> ["vary"]>
        ... (repeat Nspot times) ...
    ["fluxinput"] ["fluxoutput"] ["correctlc"]
    ["omodel" <outdir ["nameformat" fmt]> ["tdelv"]]
    ["ocurve" <"outdir" outdir ["nameformat" fmt]> ["tdelv"] ["step" stepsize]]
```

Fits or injects Kipping's Macula analytic model for starspot modulation. Specify
`"inject"` to add the model signal to the light curve, or `"fit"` to optimise
the parameters (use `"amoeba"` for Nelder-Mead simplex, `"lm"` for
Levenberg-Marquardt).

**Stellar parameters:**

| Parameter | Description |
|---|---|
| `Prot` | Equatorial rotation period (light curve time units). |
| `istar` | Stellar inclination (radians). |
| `kappa2`, `kappa4` | Quadratic and quartic differential-rotation coefficients. |
| `c1`–`c4` | Four-coefficient stellar limb-darkening terms. |
| `d1`–`d4` | Four-coefficient spot limb-darkening terms. |
| `blend` | Blend parameter. |

**Per-spot parameters** (repeated `Nspot` times):

| Parameter | Description |
|---|---|
| `Lambda0` | Longitude at maximum spot size (radians). |
| `Phi0` | Latitude at maximum spot size (radians). |
| `alphamax` | Maximum angular radius (radians). |
| `fspot` | Spot-to-star flux contrast. |
| `tmax` | Reference epoch of maximum spot size. |
| `life` | Spot lifetime, FWHM (light curve time units). |
| `ingress` | Spot growth duration. |
| `egress` | Spot decay duration. |

Append `"vary"` to any parameter to allow it to vary during fitting.

`"fluxinput"` and `"fluxoutput"` toggle flux vs. magnitude input/output
(default: magnitudes). `"correctlc"` subtracts the model. `"omodel"` writes
model files; adding `"tdelv"` includes predicted transit depth variations.
`"ocurve"` writes a uniformly sampled model curve.

**Citation:** Kipping, D. M. 2012, arXiv:1209.2985.

---

## `-magadd`

**Add a constant offset to light curve magnitudes.**

```
-magadd
    <"fix" value | "list" ["column" col] | "fixcolumn" <colname | colnum>>
```

Adds a scalar offset to every magnitude value in the light curve. This is a
simple template command included to demonstrate how to write user-defined
VARTOOLS extensions.

| Source keyword | Behavior |
|---|---|
| `"fix" value` | Add the same constant to all light curves. |
| `"list"` | Read the offset from the input list for each light curve; use `"column" col` to specify the column (default: next free column). |
| `"fixcolumn" colname/colnum` | Take the offset from a previously computed output statistic. |

---

## `-python`

**Embed Python expressions or scripts in a VARTOOLS pipeline.**

```
-python
    < "fromfile" commandfile | commandstring >
    ["init" <"file" initializationfile | initializationstring>
        | "continueprocess" prior_python_command_number]
    ["vars" variablelist
        | ["invars" inputvariablelist] ["outvars" outputvariablelist]]
    ["outputcolumns" variablelist] ["process_all_lcs"]
```

**Examples**

**Example 1.** Use Python to calculate the variance in magnitudes for each light curve.

```bash
vartools -l EXAMPLES/lc_list \
         -inputlcformat t:1,mag:2,err:3 \
         -header \
         -python 'b = numpy.var(mag)' \
                  invars mag outvars b outputcolumns b
```

Output: table showing variance values for each light curve file.

**Example 2.** Python script file `EXAMPLES/plotlc.py` containing matplotlib plotting function. Command executing Lomb-Scargle periodogram with conditional python execution for files meeting probability threshold. Uses the `-python` command in VARTOOLS to load initialization code and execute plotting on filtered light curves.

**Example 3.** Similar to Example 2, but uses the `process_all_lcs` keyword to send all of the light curves to python at once. Variables supplied as lists of numpy arrays with for-loop iteration.

---

## `-R`

**Embed R expressions or scripts in a VARTOOLS pipeline.**

```
-R
    < "fromfile" commandfile | commandstring >
    ["init" <"file" initializationfile | initializationstring>
        | "continueprocess" prior_R_command_number]
    ["vars" variablelist
        | ["invars" inputvariablelist] ["outvars" outputvariablelist]]
    ["outputcolumns" variablelist] ["process_all_lcs"]
```

**Examples**

**Example 1.** Use R to calculate the standard deviation in magnitudes for each light curve.

```bash
vartools -l EXAMPLES/lc_list \
         -inputlcformat t:1,mag:2,err:3 \
         -header \
         -R 'b <- sd(mag)' \
             invars mag outvars b outputcolumns b
```

Output:
```
#Name R_b_0
EXAMPLES/1 0.15946976931434592
EXAMPLES/2 0.036640196913116818
EXAMPLES/3 0.0048962905656505422
EXAMPLES/4 0.0020915710522042882
EXAMPLES/5 0.002880850234933455
EXAMPLES/6 0.0020898736803245783
EXAMPLES/7 0.003488095003079855
EXAMPLES/8 0.0022502571019889705
EXAMPLES/9 0.0018673694762206033
EXAMPLES/10 0.0023627959129451301
```

**Example 2.** Same calculation using `process_all_lcs` to send all light curves to R simultaneously.

```bash
vartools -l EXAMPLES/lc_list \
         -inputlcformat t:1,mag:2,err:3 \
         -header \
         -R 'b <- list();
             for(i in 1:length(mag)) {
                b[[i]] <- sd(mag[[i]]);
             }' \
             invars mag outvars b outputcolumns b process_all_lcs
```

Output: (same as Example 1)

**Example 3.** Apply ARIMA modeling using the `forecast` package. After reading light curves, the process saves the original magnitude vector, bins the data temporally with 0.05-day bins, and resamples onto a uniform grid. The R command creates a time series object, performs ARIMA fitting, and generates modeled values. Libraries are loaded via the `init` parameter. Subsequent resampling restores data to the original timebase, and the original magnitude vector is restored before output.

```bash
vartools -l EXAMPLES/lc_list \
    -inputlcformat t:1,mag:2,err:3 \
    -header \
    -savelc \
    -binlc average binsize 0.05 taverage \
    -resample linear delt fix 0.05 \
    -R \
       'mag_ts <- ts(mag, start=1, end=length(mag), frequency=1);
          arima_model <- auto.arima(mag_ts);
          mag_arima <- mag - as.vector(arima_model$residuals);' \
       init 'library(tseries); library(forecast);' \
       invars mag outvars mag_arima \
    -resample linear file list listcolumn 1 tcolumn 1 \
    -restorelc 1 vars mag \
    -o EXAMPLES/OUTDIR1 nameformat '%s.arimamodel' \
       columnformat t,mag,mag_arima \
  2> /dev/null
```

Output (first 3 lines of `EXAMPLES/OUTDIR1/1.arimamodel`):
```
53725.173920000001 10.085000000000001 10.069214881731755
53725.17654 10.0847 10.070167733349381
53725.17772 10.0825 10.070596880260995
```

**Example 4.** Sources initialization code from an external R file defining a custom function that performs ARIMA modeling and generates PNG plots. A Python call extracts basenames by removing directory paths, demonstrating language interoperability.

`EXAMPLES/Rexample4.R`:

```r
library(tseries)
library(forecast)

DoArimaFitPlot <- function(mag, plotoutdir, lcbasename) {
    mag_ts <- ts(mag, start=1, end=length(mag), frequency=1)
    arima_model <- auto.arima(mag_ts)
    mag_arima <- mag - as.vector(arima_model$residuals)
    mag_forecast <- forecast(arima_model)
    png(paste(plotoutdir, lcbasename, ".arimaforecast.png", sep=""), width=640, height=480)
    plot(mag_forecast)
    dev.off()
    png(paste(plotoutdir, lcbasename, ".arimaresiduals.png", sep=""), width=640, height=480)
    checkresiduals(arima_model)
    dev.off()
    return(mag_arima)
}
```

```bash
vartools -l EXAMPLES/lc_list \
    -inputlcformat t:1,mag:2,err:3 \
    -header \
    -savelc \
    -binlc average binsize 0.05 taverage \
    -resample linear delt fix 0.05 \
    -python 'lcbasename = Name.split("/")[-1]' \
        invars Name outvars lcbasename \
    -R 'mag_arima <- DoArimaFitPlot(mag, "EXAMPLES/OUTDIR1/", lcbasename)' \
        init file EXAMPLES/Rexample4.R \
        invars mag,t,lcbasename outvars mag_arima \
    -resample linear file list listcolumn 1 tcolumn 1 \
    -restorelc 1 vars mag \
    -o EXAMPLES/OUTDIR1 nameformat '%s.arimamodel' \
        columnformat t,mag,mag_arima \
  2> /dev/null
```

---

## `-generic`

**Template for writing user-defined VARTOOLS extension commands.**

The `-generic` command is an incomplete documentation template demonstrating the structure required to implement a custom VARTOOLS extension. Refer to the USERLIB source tree for working examples such as `-magadd`.
