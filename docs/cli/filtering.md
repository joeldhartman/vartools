# Filtering and Detrending

Commands for removing outliers, applying smoothing filters, restricting the time range of light curves, and saving and restoring light curve states.

---

## `-clip`

```
-clip <"var" sigclipvar | "expr" sigclipexpr | sigclip>
    <"var" itervar | "expr" iterexpr | iter>
    ["niter" <"var" nvar | "expr" nexpr | n>] ["median"]
    ["markclip" var ["noinitmark"]]
```

Sigma-clip the light curves. Points with errors ≤ 0 or NaN magnitude values are always removed. If `sigclip ≤ 0`, sigma-clipping is not performed but invalid points (err ≤ 0 or NaN magnitude) are still removed. The output table includes the number of points that were clipped.

**Parameters**

- `sigclip` — Clipping threshold in units of the standard deviation.
- `iter` — `1` for iterative clipping (continue until no further points are removed); `0` for a single pass.
- `"niter" n` — Clip for a fixed number of iterations `n`.
- `"median"` — Clip with respect to the median rather than the mean.
- `"markclip" var` — Instead of removing clipped points, mark them: points that survive clipping receive `var = 1`; clipped points receive `var = 0`. The light curve length is unchanged.
  - `"noinitmark"` — Use the existing values of `var` as an initial clipping mask. Only points with `var = 1` are considered for further clipping.

**Examples**

**Example 1.** Calculate the RMS before and after applying iterative 3-sigma clipping, then write the clipped light curve to a file.

```bash
vartools -oneline -i EXAMPLES/5 \
    -rms \
    -clip 3. 1 \
    -rms \
    -o EXAMPLES/OUTDIR1/5.clip.txt
```

Output:
```
Name           = EXAMPLES/5
Mean_Mag_0     =  10.43962
RMS_0          =   0.00288
Expected_RMS_0 =   0.00114
Npoints_0      =  3903
Nclip_1        =    51
Mean_Mag_2     =  10.43961
RMS_2          =   0.00267
Expected_RMS_2 =   0.00114
Npoints_2      =  3852
```

---

## `-medianfilter`

```
-medianfilter
    time ["average" | "weightedaverage"] ["replace"]
```

Apply a running high-pass or low-pass filter to the light curve.

**Parameters**

- `time` — Half-window width in the same units as the time coordinate. All points within `time` of a given observation contribute to the local statistic.
- Default (high-pass) behavior: the local median magnitude is **subtracted** from each observation.
- `"average"` — Use the running mean instead of the median.
- `"weightedaverage"` — Use the running weighted mean instead of the median.
- `"replace"` — **Replace** each point with the running statistic rather than subtracting it. This converts the command to a low-pass filter.

**Examples**

**Example 1.** Apply a high-pass and a low-pass median filter (0.05-day timescale) to a light curve, saving and restoring the light curve state between the two filter passes.

```bash
vartools -i EXAMPLES/1 -oneline -chi2 \
    -savelc \
    -medianfilter 0.05 \
    -chi2 -o EXAMPLES/OUTDIR1/1.medianhighpass \
    -restorelc 1 \
    -medianfilter 0.05 replace \
    -chi2 -o EXAMPLES/OUTDIR1/1.medianlowpass
```

Output:
```
Name                 = EXAMPLES/1
Chi2_0               = 34711.71793
Weighted_Mean_Mag_0  =  10.24430
Chi2_3               =   5.95454
Weighted_Mean_Mag_3  =  -0.00009
Chi2_7               = 34727.65120
Weighted_Mean_Mag_7  =  10.24440
```

---

## `-restricttimes`

```
-restricttimes
    ["exclude"]
    < "JDrange" minJD maxJD |
      "JDrangebylc"
        < "fix" minJD | "list" ["column" col] | "fixcolumn" <colname | colnum>
         "expr" expression >
        < "fix" maxJD | "list" ["column" col] | "fixcolumn" <colname | colnum>
         "expr" expression >
      "JDlist" JDfilename |
      "imagelist" imagefilename |
      "expr" expression >
    ["markrestrict" markvar ["noinitmark"]]
```

Filter observations from the light curve based on their time values or string IDs.

By default, only points **matching** the specified criterion are kept. Use `"exclude"` to instead **remove** matching points.

**Time selection methods**

| Keyword | Description |
|---------|-------------|
| `"JDrange" minJD maxJD` | Keep points with times in `[minJD, maxJD]` (same range for all light curves) |
| `"JDrangebylc"` | Like `"JDrange"` but allows a different range per light curve; values come from `"fix"`, `"list"`, `"fixcolumn"`, or `"expr"` |
| `"JDlist" JDfilename` | Keep (or exclude) points whose times appear in this file (column 1: JD) |
| `"imagelist" imagefilename` | Keep (or exclude) points whose string IDs appear in this file (column 1: image name) |
| `"expr" expression` | Keep (or exclude) points for which the analytic expression evaluates to a value > 0; e.g., `'-restricttimes expr "(mag>9.0)&&(mag<9.5)"'` keeps only points in the magnitude range 9.0–9.5 |

**Additional parameters**

- `"markrestrict" markvar` — Instead of removing points, mark them: points that would be kept receive `markvar = 1`; points that would be removed receive `markvar = 0`. Note: `-restoretimes` cannot be used with a `-restricttimes` that uses this keyword.
  - `"noinitmark"` — Use the existing `markvar` values as an initial mask.

!!! tip
    Use `-restricttimes` and `-restoretimes` together to apply modifications to isolated segments of a light curve.

---

## `-restoretimes`

```
-restoretimes
    prior_restricttimes_command
```

Restore observations that were filtered out by a prior `-restricttimes` command. The restored points are appended to the current light curve and the light curve is then sorted by time.

**Parameters**

- `prior_restricttimes_command` — Integer index of the `-restricttimes` command to restore from. `1` refers to the first `-restricttimes` call on the command line, `2` to the second, and so on.

!!! note
    Cannot be used with a `-restricttimes` command that used the `"markrestrict"` keyword.


---

## `-SYSREM`

```
-SYSREM
    Ninput_color ["column" col1]
    Ninput_airmass initial_airmass_file
    sigma_clip1 sigma_clip2 saturation correctlc
    omodel [model_outdir] otrends [trend_outfile]
    useweights
```

Run the SYSREM PCA-like algorithm to identify and remove ensemble trends from a set of light curves. This command requires a light curve list and automatically sets the `-readall` option.

**Parameters**

- `Ninput_color` — Number of initial color-term trends; their values are read from the input light curve list.
- `"column" col1` — Column in the input list for the first color term (subsequent terms follow in order).
- `Ninput_airmass` — Number of initial airmass-term trends.
- `initial_airmass_file` — File with initial airmass trends (column 1: JD; subsequent columns: trend values).
- `sigma_clip1` — σ-clipping for computing mean magnitudes.
- `sigma_clip2` — σ-clipping for determining which points contribute to the airmass/color terms.
- `saturation` — Points with magnitude below this value do not contribute to the fit.
- `correctlc` — `1` to subtract the model; `0` to compute without subtracting.
- `omodel` — `1` to output model light curves to `model_outdir`. Output format: `JD mag mag_model sig clip`. Suffix: `.sysrem.model`.
- `otrends` — `1` to output the final trend signals to `trend_outfile` (column 1: JD, subsequent columns: trend values).
- `useweights` — Include this flag to weight observations by their formal uncertainties.

**Citation:** Tamuz, Mazeh and Zucker, 2005, MNRAS, 356, 1466.

---

## `-TFA`

```
-TFA
    trendlist ["readformat" Nskip jdcol magcol]
    ["trend_coeff_priors" trend_coeff_prior_file
        ["use_lc_errors" | "weight_by_template_stddev"]]
    dates_file pixelsep ["xycol" xcol ycol]
    correctlc ocoeff [coeff_outdir] omodel [model_outdir]
    ["clip" sigclipfactor ["usemedian"] ["useMAD"]]
    ["fitmask" maskvar] ["outfitmask" outmaskvar]
```

Run the Trend Filtering Algorithm (TFA) on the light curves. TFA fits each light curve as a linear combination of a set of template (basis) light curves and subtracts the fit, yielding a filtered, detrended light curve.

A light curve list (`-l`) is required. The `x` and `y` pixel positions of each light curve must be given as columns in the list.

**Parameters**

- `trendlist` — File containing a list of basis vector files in the format: `trendname trendx trendy`. Files can be ASCII or binary FITS. Use `"readformat" Nskip jdcol magcol` to specify the format (defaults: `Nskip=0`, `jdcol=1`, `magcol=2`).
- `"trend_coeff_priors" trend_coeff_prior_file` — File containing Gaussian priors for the trend coefficients (columns: `trendname prior_mean prior_stddev`).
  - `"use_lc_errors"` — Weight light curve points by `1/err[i]` (more correct but slower).
  - `"weight_by_template_stddev"` — Weight points by `1/ave_template_stddev`.
- `dates_file` — File with the complete set of JDs for all light curves (column 1: filename/string-id, column 2: JD).
- `pixelsep` — Basis vectors with coordinates within `pixelsep` of the target are excluded from its detrending (to avoid self-filtering).
- `"xycol" xcol ycol` — Columns in the input list giving x and y positions (default: next two available columns).
- `correctlc` — `1` to apply the filter; `0` to compute but not subtract.
- `ocoeff` — `1` to output trend coefficients to `coeff_outdir`. Output suffix: `.tfa.coeff`.
- `omodel` — `1` to output the TFA model to `model_outdir`. Output suffix: `.tfa.model`.
- `"clip" sigclipfactor` — Clipping level for outlier rejection before fitting (default: 5σ). Add `"usemedian"` and/or `"useMAD"` to change the reference statistic.
- `"fitmask" maskvar` — Restrict points included in the trend fit (1 = include, 0 = exclude). Model is still evaluated and subtracted at excluded points.
- `"outfitmask" outmaskvar` — Store the actual fit mask (after clipping) in this variable.

**Citation:** Kovacs, Bakos and Noyes, 2005, MNRAS, 356, 557.

---

## `-TFA_SR`

```
-TFA_SR
    trendlist ["readformat" Nskip jdcol magcol] dates_file
    ["decorr" iterativeflag Nlcterms lccolumn1 lcorder1 ...]
    pixelsep ["xycol" colx coly]
    correctlc ocoeff [coeff_outdir] omodel [model_outdir]
    dotfafirst tfathresh maxiter
    < "bin" nbins
            ["period" < "aov" | "ls" | "bls" | "list" ["column" col] | "fix" period >]
        | "signal" filename
        | "harm" Nharm Nsubharm
            ["period" < "aov" | "ls" | "bls" | "list" ["column" col] | "fix" period >] >
    ["clip" sigclipfactor ["usemedian"] ["useMAD"]]
    ["fitmask" maskvar] ["outfitmask" outmaskvar]
```

Run TFA in Signal Reconstruction (SR) mode. TFA-SR iteratively applies TFA and fits a signal model to the light curve, allowing the algorithm to preserve astrophysical signal that would otherwise be partially filtered by TFA.

Most syntax is identical to [`-TFA`](#-tfa). Parameters specific to TFA-SR are described below.

**Parameters**

- `"decorr" iterativeflag Nlcterms lccolumn1 lcorder1 ...` — Simultaneously decorrelate against `Nlcterms` light-curve-specific signals. `iterativeflag=1` for iterative decorrelation and TFA (faster); `iterativeflag=0` for simultaneous (more correct but slower).
- `dotfafirst` — `1` to apply TFA first in each iteration, then fit the signal to the residual; `0` to subtract the signal first, then apply TFA to the residual.
- `tfathresh` — Stop iterating when the fractional change in RMS falls below this threshold.
- `maxiter` — Maximum number of iterations.
- Signal model (choose one):
  - `"bin" nbins` — Mean binned light curve with `nbins` bins. Use optional `"period"` keyword for phase-folding.
  - `"signal" filename` — Fixed signal form read from a file. The file contains a list of signal files (one per light curve), with the signal magnitude in the second column. Fits `a*signal + b`.
  - `"harm" Nharm Nsubharm` — Fourier series fit simultaneously with TFA (no iteration in this case). Use optional `"period"` to specify the period source.

**Citation:** Kovacs, Bakos and Noyes, 2005, MNRAS, 356, 557.

---

