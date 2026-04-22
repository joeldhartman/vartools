# Extension Commands

Extension commands are additional VARTOOLS commands that are compiled as
separate shared-object libraries (`.so` files) and loaded at runtime. They are
distributed in the `USERLIBS/src` subdirectory of the VARTOOLS source tree.
After running `make` in that directory, each extension produces a `.so` file
in `USERLIBS/src/.libs/` that must be passed to VARTOOLS with the `-L` option
before the command itself is issued:

```bash
# Load the magadd extension and use its -magadd command to add a constant
# offset to each light curve before measuring the RMS.
vartools -L USERLIBS/src/.libs/magadd.so \
    -l EXAMPLES/lc_list \
    -magadd fix 0.5 \
    -rms -tab
```

The same `-L` flag can be given multiple times to load several extensions in a
single call.

!!! note "Installed extensions"
    After `make install`, the compiled extension libraries are placed in
    `$(PREFIX)/share/vartools/USERLIBS/` and VARTOOLS finds them automatically
    — the explicit `-L` prefix is only needed when running from the source
    tree or when the library lives outside the installed search directory.

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
