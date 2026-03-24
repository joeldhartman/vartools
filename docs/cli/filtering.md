# Filtering and Detrending

Commands for removing outliers, applying smoothing filters, restricting the time range of light curves, and saving and restoring light curve states.

---

## `-clip`

```
-clip
    sigclip iter ["niter" n] ["median"]
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

## `-decorr`

```
-decorr
    correctlc zeropointterm subtractfirstterm
    Nglobalterms globalfile1 order1 ... globalfileN orderN
    Nlcterms lccolumn1 lcorder1 ... lccolumnN lcorderN
    omodel [modeloutdir] ["maskpoints" maskvar]
```

!!! warning "Deprecated"
    This command is deprecated as of VARTOOLS 1.3. Use [`-linfit`](model-fitting.md#-linfit) instead.

Decorrelate the light curves against specified external or light-curve-specific signals using polynomial regression.

**Parameters**

- `correctlc` — `1` to apply the decorrelation to the light curve; `0` to compute and output the coefficients and χ² without modifying the light curve.
- `zeropointterm` — `1` to include a zero-point offset term in the fit; `0` to omit it.
- `subtractfirstterm` — `1` to decorrelate against `(signal - signal[0])` rather than `signal` directly (useful for detrending against JD).
- `Nglobalterms` — Number of global signal files.
- `globalfile1 ... globalfileN` — Names of global signal files (format: `JD signal_value`).
- `order1 ... orderN` — Polynomial orders for each global signal (must be ≥ 1).
- `Nlcterms` — Number of light-curve-specific signals.
- `lccolumn1 ... lccolumnN` — Column indices in the light curve for each light-curve-specific signal.
- `lcorder1 ... lcorderN` — Polynomial orders for each light-curve-specific signal.
- `omodel` — `1` to output the decorrelation model to `modeloutdir`. Suffix: `.decorr.model`.
- `"maskpoints" maskvar` — Optional. Only points with `maskvar > 0` contribute to the fit.

**Examples**

**Example 1.** Fit quadratic polynomials to light curves using a JD-based light-curve term (column 1), including a zero-point offset, with the first term subtracted to reduce rounding errors. Report RMS before and after decorrelation.

```bash
vartools -l EXAMPLES/lc_list -header \
    -rms \
    -decorr 1 1 1 0 1 1 2 0 \
    -rms
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

## `-savelc`

```
-savelc
```

Save the current state of the light curve in memory. Use in conjunction with `-restorelc` to return the light curve to this state later. Multiple `-savelc` calls create a numbered sequence of checkpoints.

**Example**

Apply TFA with three different template lists without re-reading the original light curve:

```bash
vartools -l lclist ... -savelc \
    -TFA trendlist1 ... -LS ... \
    -restorelc 1 -TFA trendlist2 ... -LS ... \
    -restorelc 1 -TFA trendlist3 ... -LS ...
```

!!! caution
    Conditional constructs (`-if`/`-elif`/`-else`/`-fi`) are ignored by `-savelc`.

---

## `-restorelc`

```
-restorelc
    savenumber ["vars" var1,var2,...]
```

Restore the light curve to a state saved by a prior `-savelc` call.

**Parameters**

- `savenumber` — Which checkpoint to restore: `1` for the first `-savelc`, `2` for the second, etc.
- `"vars" var1,var2,...` — Optional. Restore only the listed light curve vectors instead of the entire light curve. Partial restoration only proceeds if the saved vectors are the same length as the current light curve; otherwise a warning is printed to stderr.

!!! caution
    Conditional constructs (`-if`/`-elif`/`-else`/`-fi`) are ignored by `-restorelc`.
