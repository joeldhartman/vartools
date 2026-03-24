# Light Curve Manipulation

This page documents the VARTOOLS commands that transform, filter, reformat, or inspect light curves. Commands are applied in the order given on the command line; each receives the (possibly modified) output of the previous command.

---

### `-binlc`

**Syntax**
```
-binlc
    < "average" | "median" | "weightedaverage" >
    < "binsize" binsize | "nbins" nbins >
    ["bincolumns" var1[:stats1][,var2[:stats2],...]]
    ["T0"
        < "fix" T0val | "list" ["column" col] | "fixcolumn" <colname | colnum> |
          "expr" expression >]
    ["firstbinshift" firstbinshift]
    < "tcenter" | "taverage" | "tmedian" | "tnoshrink" ["bincolumnsonly"] >
    ["maskpoints" maskvar]
```

**Description**

Bin the light curve in time (or in phase if a `-Phase` command has already been applied). All light curve vectors are binned together by default.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `average` / `median` / `weightedaverage` | Statistic used to combine points in a bin. Legacy integers `0`, `1`, `2` are also accepted. |
| `"binsize" binsize` | Bin width in units of the time (or phase) coordinate. |
| `"nbins" nbins` | Number of equal-width bins to divide the time span into. |
| `"bincolumns" var1[:stats1],...` | Override the binning statistic for specific named columns. Statistic names follow `-stats` conventions. |
| `"T0" ...` | Start time of the first bin. Sources: `"fix"` (command-line value), `"list"` (input list column), `"fixcolumn"` (prior output column), or `"expr"` (analytic expression). |
| `"firstbinshift" firstbinshift` | Shift the first bin start by `t0 - firstbinshift/binsize`. |
| `"tcenter"` | Output time for each bin is the bin center. |
| `"taverage"` | Output time is the average of times in the bin. |
| `"tmedian"` | Output time is the median of times in the bin. |
| `"tnoshrink"` | Replace every point with its binned value without reducing the length of the light curve. Append `"bincolumnsonly"` to restrict replacement to explicitly listed columns. |
| `"maskpoints" maskvar` | Only include points with `maskvar > 0` in the binning. Masked-out points still receive the binned value when `tnoshrink` is active. |

**Examples**

**Example 1.** Median-bin `EXAMPLES/2` in time with a bin width of 0.01 days.

```bash
vartools -i EXAMPLES/2 \
    -binlc median binsize 0.01 tcenter \
    -o EXAMPLES/OUTDIR1/2.bin.txt
```

**Example 2.** Run an LS period search, phase-fold on the best period, then bin the phased light curve into 100 equal phase bins.

```bash
vartools -i EXAMPLES/2 \
    -LS 0.1 10.0 0.1 1 0 \
    -Phase ls \
    -binlc median nbins 100 tcenter \
    -o EXAMPLES/OUTDIR1/2.phasebin.txt
```

---

### `-changeerror`

**Syntax**
```
-changeerror
    ["maskpoints" maskvar]
```

**Description**

Replace the formal per-point measurement uncertainties in the light curve with the RMS of the light curve. This is useful when formal errors are unavailable or unreliable.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `"maskpoints" maskvar` | Restrict the RMS calculation to points with `maskvar > 0`. |

**Examples**

**Example 1.** Replace formal errors with the light curve RMS and verify that χ²/dof becomes ≈ 1.

```bash
vartools -i EXAMPLES/4 -ascii -oneline \
    -chi2 \
    -changeerror \
    -chi2
```

Output:
```
Name                = EXAMPLES/4
Chi2_0              =      5.19874
Weighted_Mean_Mag_0 =  10.35137
Mean_Mag_1          =  10.35142
RMS_1               =   0.00209
Npoints_1           =  3227
Chi2_2              =      1.00031
Weighted_Mean_Mag_2 =  10.35142
```

---

### `-changevariable`

**Syntax**
```
-changevariable
    < "t" | "mag" | "err" | "id" > var
```

**Description**

Reassign the internal role of the time (`t`), magnitude (`mag`), magnitude uncertainty (`err`), or image-identifier (`id`) variable to a different named column. Subsequent commands will use the new assignment. The original variable still exists; to restore it, issue `-changevariable mag mag`.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `"t"` / `"mag"` / `"err"` / `"id"` | Which built-in role to reassign. |
| `var` | Name of the existing light-curve column to promote to that role. |

**Example**

```bash
vartools -i my.fits \
    -inputlcformat t:1,pdcsap_flux:8,pdcsap_flux_err:9 \
    -changevariable mag pdcsap_flux \
    -changevariable err pdcsap_flux_err \
    -fluxtomag 25.0 0 \
    -rms -oneline
```

**Examples**

**Example 1.** Run an LS period search, then use `-changevariable` to store the current time in a `phase` variable, phase-fold on the LS period, then swap `t` back before writing output so the light curve is sorted by time.

```bash
vartools -l EXAMPLES/lc_list \
    -LS 0.1 100.0 0.1 1 0 \
    -expr 'phase=t' \
    -changevariable t phase \
    -Phase ls \
    -changevariable t t \
    -o EXAMPLES/OUTDIR1 nameformat "%s.phase.txt" \
        columnformat "t:%17.9f,mag:%9.5f,err:%9.5f,phase:%9.5f" \
    -header
```

Output: a table with columns `Name`, `LS_Period_1_0`, `Log10_LS_Prob_1_0`, `LS_SNR_1_0` for each light curve in the list, plus one phase-folded output file per light curve.

---

### `-converttime`

**Syntax**
```
-converttime
    < "input" < "mjd" | "jd" | "hjd" | "bjd" > >
    ["inputsubtract" value] ["inputsys-tdb" | "inputsys-utc"]
    < "output" < "mjd" | "jd" | "hjd" | "bjd" > >
    ["outputsubtract" value] ["outputsys-tdb" | "outputsys-utc"]
    ["radec" < "list" ["column" col] | "fix" raval decval |
              "expr" raexpr decexpr | "perJD expr" raexpr decexpr |
              "perJD fromlc" ralccol declccol >
        ["epoch" epoch]]
    ["ppm" < "list" ["column" col] | "fix" mu_ra mu_dec >]
    ["input-radec" ...]
    ["input-ppm" ...]
    ["ephemfile" file] ["leapsecfile" file] ["planetdatafile" file]
    ["observatory" < code | "show-codes" >
      | "coords" ["xyz"]
          < "fix" lat lon alt | "list" [...] | "expr" [...] |
            "fromlc" col col col | "lcexpr" [...] >]
```

**Description**

Convert the time system of the light curve between Modified Julian Date (MJD), Julian Date (JD), Heliocentric Julian Date (HJD), and Barycentric Julian Date (BJD). BJD conversion requires VARTOOLS to be linked to the JPL NAIF cspice library. The internal precision near J2000.0 is approximately 0.1 milliseconds.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `"input" <system>` | Specify the input time system: `"mjd"`, `"jd"`, `"hjd"`, or `"bjd"`. |
| `"inputsubtract" value` | Constant subtracted from the stored input times (e.g. `2400000` if input is `HJD-2400000`). |
| `"inputsys-utc"` / `"inputsys-tdb"` | Whether input JDs were computed from UTC (default, typical for ground-based observatory headers) or from TDB (requires cspice). |
| `"output" <system>` | Desired output time system. |
| `"outputsubtract" value` | Subtract a constant from the output times. |
| `"radec" ...` | Source for the RA/Dec coordinates needed for HJD or BJD conversion. Options: `"list"`, `"fix" ra dec`, `"expr"`, `"perJD expr"` (per-observation expression), `"perJD fromlc"` (columns in the light curve). Coordinates in degrees; default epoch J2000.0. |
| `"ppm" ...` | Proper motion in mas/yr (RA, Dec). |
| `"input-radec"` / `"input-ppm"` | Separate RA/Dec and proper motion used for the input time system, if it differs from the target. |
| `"observatory" code` | Observatory code (use `"show-codes"` to list). |
| `"coords" ...` | Specify observer latitude (deg), longitude (deg E), altitude (m); or XYZ coordinates in the J2000 frame with `"xyz"`. |
| `"ephemfile"` / `"leapsecfile"` / `"planetdatafile"` | JPL NAIF kernel files for BJD/TDB conversions. Can also be set via environment variables `CSPICE_EPHEM_FILE`, `CSPICE_LEAPSEC_FILE`, `CSPICE_PLANETDATA_FILE`. |

**Examples**

**Example 1.** Convert a light curve from JD to HJD given RA/Dec coordinates, subtracting 2400000 from the input JD.

```bash
vartools -i EXAMPLES/1 -quiet \
    -converttime input jd inputsubtract 2400000. output hjd \
    radec 88.079166 32.5533 \
    -o EXAMPLES/OUTDIR1/1.hjd
```

**Example 2.** Convert from UTC timestamps to Barycentric Julian Date in the TDB reference frame, requiring CSPICE kernel files and an observatory specification.

```bash
vartools -i EXAMPLES/1.UTC -quiet \
    -readformat 0 inpututc '%Y-%M-%DT%h:%m:%s' 1 2 3 \
    -converttime input utc output bjd tdb \
    radec 88.079166 32.5533 \
    obs flwo \
    leapsecfile /path/to/naif0012.tls \
    spkfile /path/to/de430.bsp \
    -o EXAMPLES/OUTDIR1/1.bjd_tdb
```

---

### `-difffluxtomag`

**Syntax**
```
-difffluxtomag
    mag_constant offset ["magcolumn" col]
```

**Description**

Convert a light curve from ISIS image-subtraction differential flux to magnitudes. Requires a light curve list (`-l`) where the reference magnitude of each star (after aperture correction) is provided as an additional column. The conversion formula is:

```
mag = mag_constant - 2.5 * log10(ref_flux + diff_flux) + offset
```

This command produces no statistics output to stdout; call `-rms` or `-chi2` separately if statistics are needed.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `mag_constant` | Magnitude corresponding to a flux of 1 ADU. |
| `offset` | Additive constant applied to the output magnitudes. |
| `"magcolumn" col` | Column in the input list containing the reference magnitude. Default: next available column. |

**Examples**

**Example 1.** Convert ISIS differential fluxes to magnitudes. The reference magnitude (10.085) is piped via stdin as the second column of the list.

```bash
echo "EXAMPLES/1 10.085" | \
  vartools -l - -header \
           -difffluxtomag 25.0 0.0 -rms
```

Output:
```
#Name Mean_Mag_1 RMS_1 Expected_RMS_1 Npoints_1
EXAMPLES/1  10.08501   0.00000   0.00000  3122
```

**Example 2.** Same as Example 1, but explicitly specify that the reference magnitude is in column 2.

```bash
echo "EXAMPLES/1 10.085" | \
  vartools -l - -header \
           -difffluxtomag 25.0 0.0 magcolumn 2 -rms
```

---

### `-ensemblerescalesig`

**Syntax**
```
-ensemblerescalesig
    sigclip ["maskpoints" maskvar]
```

**Description**

Rescale magnitude uncertainties across an ensemble of light curves by fitting a linear relation between E(rms)² and χ²·E(rms)² for all light curves, where E(rms) is the expected RMS based on the photometric uncertainties. The result is that χ²/dof is distributed about unity across the ensemble. Requires a list input; all light curves are read into memory.

The output includes the average rescale factor for each light curve, defined as `sqrt(chi2_after / chi2_before)`.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `sigclip` | Sigma-clipping factor used to reject outliers in the ensemble χ² distribution. |
| `"maskpoints" maskvar` | Restrict the χ²/dof and expected RMS calculation to points with `maskvar > 0`. |

**Examples**

**Example 1.** Transform magnitude uncertainties across an ensemble of light curves, using `-chi2` before and after to demonstrate the rescaling effect.

```bash
vartools -l EXAMPLES/lc_list -header \
    -chi2 -ensemblerescalesig 3. -chi2
```

Output: table with columns Name, Chi2_0, Weighted_Mean_Mag_0, SigmaRescaleFactor_1, Chi2_2, Weighted_Mean_Mag_2 for EXAMPLES/1 through EXAMPLES/10.

---

### `-expr`

**Syntax**
```
-expr
    var"="expression
```

**Description**

Evaluate an analytic expression and assign the result to a named variable. If the variable does not yet exist it is created; it can be a scalar or a per-point vector depending on the right-hand side. The expression can reference any existing light curve vectors (`t`, `mag`, `err`, other named columns), scalars from prior commands, or output columns identified by their header names (with leading digits and underscores stripped). Run `vartools -functionlist` for the full list of supported functions, operators, and constants.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `var` | Name of the variable to create or update. Cannot be an output-column name. |
| `expression` | Any analytic expression supported by the VARTOOLS expression evaluator. |

**Example**

```bash
# Convert magnitudes to linear flux and store as a new column
vartools -i lc.txt -expr 'flux=10^(-0.4*(mag-25.0))' -rms -oneline
```

**Examples**

**Example 1.** Add a constant to all magnitude values and then take the square root.

```bash
vartools -i EXAMPLES/1 -expr 'mag=sqrt(mag+5)' -o EXAMPLES/1.add
```

Input (first 3 lines of EXAMPLES/1):
```
53725.173920 10.085 0.00119
53725.176540 10.0847 0.00144
53725.177720 10.0825 0.00123
```

Output (first 3 lines of EXAMPLES/1.add):
```
  53725.173920000   3.88394   0.00119
  53725.176540000   3.88390   0.00144
  53725.177720000   3.88362   0.00123
```

**Example 2.** Fit a sinusoid signal to each light curve and subtract it only when the delta chi² improvement is significant.

```bash
vartools -l EXAMPLES/lc_list -header \
    -LS 0.1 10. 0.1 1 0 \
    -rms -chi2 \
    -expr 'mag2=mag' \
    -Killharm ls 0 0 0 \
    -rms -chi2 \
    -expr \
        'mag=(Npoints_5*(Chi2_6-Chi2_2)<-10000)*mag+
            (Npoints_5*(Chi2_6-Chi2_2)>=-10000)*mag2' \
    -o EXAMPLES/OUTDIR1 nameformat '%s.cleanharm'
```

**Example 3.** Convert magnitudes to fluxes with `-expr`, compute the median flux, then normalize by the median.

```bash
vartools -i EXAMPLES/1 \
    -expr 'flux=10^(-0.4*(mag-25.0))' \
    -stats flux median \
    -expr 'flux=flux/STATS_flux_MEDIAN_1' \
    -stats flux,mag median,stddev \
    -oneline
```

Output:
```
Name                = EXAMPLES/1
STATS_flux_MEDIAN_1 = 842674.79516438092
STATS_flux_MEDIAN_3 = 1
STATS_flux_STDDEV_3 = 0.1290865428119792
STATS_mag_MEDIAN_3  = 10.18585
STATS_mag_STDDEV_3  = 0.15946976931434606
```

---

### `-FFT`

**Syntax**
```
-FFT
    input_real_var input_imag_var output_real_var output_imag_var
```

**Description**

Compute the Fast Fourier Transform of a light curve vector using the GSL function `gsl_fft_complex_forward()`. The input and output vectors follow the standard GSL complex layout where element `k` corresponds to frequency `k/(N·Δ)` (for `k < N/2`) or negative frequencies (for `k > N/2`), where `Δ` is the assumed uniform time step between points and `N` is the number of points.

Use `"NULL"` for either input component to substitute a zero vector. Use `"NULL"` for either output component to discard that component.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `input_real_var` | Name of the light curve vector holding the real part of the input signal, or `"NULL"`. |
| `input_imag_var` | Name of the light curve vector holding the imaginary part, or `"NULL"`. |
| `output_real_var` | Variable name to store the real part of the transform, or `"NULL"`. |
| `output_imag_var` | Variable name to store the imaginary part of the transform, or `"NULL"`. |

**Examples**

**Example 1.** High-pass Fourier filter a uniformly sampled time-series by computing the FFT, zeroing low-frequency components, and applying the inverse FFT.

```bash
vartools -i EXAMPLES/11 \
    -FFT mag NULL fftreal fftimag \
    -rms \
    -expr \
      'fftreal=(NR>(Npoints_1/500.0))*(NR<(Npoints_1*499.0/500.0))*fftreal' \
    -expr \
      'fftimag=(NR>(Npoints_1/500.0))*(NR<(Npoints_1*499.0/500.0))*fftimag' \
    -IFFT fftreal fftimag mag_filter NULL \
    -o EXAMPLES/11.highpassfftfilter columnformat t,mag_filter
```

Output:
```
EXAMPLES/11   0.00070   0.76113   0.00100 100000
```

**Example 2.** Apply FFT-based high-pass and low-pass filtering to a non-uniformly sampled light curve by resampling to a uniform grid before transforming, then resampling back to the original time base.

```bash
vartools -i EXAMPLES/2 \
    -resample splinemonotonic gaps percentile_sep 80 bspline \
    -FFT mag NULL fftreal fftimag \
    -rms \
    -expr 'fftreal1=(NR>(Npoints_2/10.0))*(NR<(Npoints_2*9.0/10.0))*fftreal' \
    -expr 'fftimag1=(NR>(Npoints_2/10.0))*(NR<(Npoints_2*9.0/10.0))*fftimag' \
    -IFFT fftreal1 fftimag1 mag_filter NULL \
    -expr 'fftreal2=fftreal-((NR>(Npoints_2/10.0))*(NR<(Npoints_2*9.0/10.0))*fftreal)' \
    -expr 'fftimag2=fftimag-((NR>(Npoints_2/10.0))*(NR<(Npoints_2*9.0/10.0))*fftimag)' \
    -IFFT fftreal2 fftimag2 mag_filter2 NULL \
    -resample linear file fix EXAMPLES/2 column 1 \
    -expr 'mag_filter=mag_filter+Mean_Mag_2' \
    -o EXAMPLES/2.fftfilter columnformat t,mag_filter,mag_filter2,mag
```

Output:
```
EXAMPLES/2  10.11155   0.02901   0.00028 66186
```

---

### `-IFFT`

**Syntax**
```
-IFFT
    input_real_var input_imag_var output_real_var output_imag_var
```

**Description**

Compute the Inverse Fast Fourier Transform of a light curve vector using the GSL function `gsl_fft_complex_backward()`. The parameter conventions are identical to those of `-FFT`, with input treated as frequency-domain data and output as time-domain data.

**Parameters**

Same as `-FFT`. See above.

**Examples**

**Example 1.** High-pass Fourier filter a uniformly sampled time-series using `-FFT` and `-IFFT` in combination.

```bash
vartools -i EXAMPLES/11 \
    -FFT mag NULL fftreal fftimag \
    -rms \
    -expr \
      'fftreal=(NR>(Npoints_1/500.0))*(NR<(Npoints_1*499.0/500.0))*fftreal' \
    -expr \
      'fftimag=(NR>(Npoints_1/500.0))*(NR<(Npoints_1*499.0/500.0))*fftimag' \
    -IFFT fftreal fftimag mag_filter NULL \
    -o EXAMPLES/11.highpassfftfilter columnformat t,mag_filter
```

Output:
```
EXAMPLES/11   0.00070   0.76113   0.00100 100000
```

**Example 2.** Apply high-pass and low-pass filtering to a non-uniformly sampled light curve by resampling to a uniform grid, transforming with `-FFT`, filtering frequency components, and inverting with `-IFFT`.

```bash
vartools -i EXAMPLES/2 \
    -resample splinemonotonic gaps percentile_sep 80 bspline \
    -FFT mag NULL fftreal fftimag \
    -rms \
    -expr 'fftreal1=(NR>(Npoints_2/10.0))*(NR<(Npoints_2*9.0/10.0))*fftreal' \
    -expr 'fftimag1=(NR>(Npoints_2/10.0))*(NR<(Npoints_2*9.0/10.0))*fftimag' \
    -IFFT fftreal1 fftimag1 mag_filter NULL \
    -expr 'fftreal2=fftreal-((NR>(Npoints_2/10.0))*(NR<(Npoints_2*9.0/10.0))*fftreal)' \
    -expr 'fftimag2=fftimag-((NR>(Npoints_2/10.0))*(NR<(Npoints_2*9.0/10.0))*fftimag)' \
    -IFFT fftreal2 fftimag2 mag_filter2 NULL \
    -resample linear file fix EXAMPLES/2 column 1 \
    -expr 'mag_filter=mag_filter+Mean_Mag_2' \
    -o EXAMPLES/2.fftfilter columnformat t,mag_filter,mag_filter2,mag
```

Output:
```
EXAMPLES/2  10.11155   0.02901   0.00028 66186
```

---

### `-fluxtomag`

**Syntax**
```
-fluxtomag
    mag_constant offset
```

**Description**

Convert light curve flux values to magnitudes. The conversion is:

```
mag = mag_constant - 2.5 * log10(flux) + offset
```

This command produces no output to stdout.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `mag_constant` | Magnitude of a source with a flux of 1 ADU. |
| `offset` | Additive constant applied to the output magnitudes. |

**Examples**

**Example 1.** Convert a Kepler public Q1 light curve from flux to magnitudes using a zero-point of 25.0 (1 ADU = magnitude 25).

```bash
vartools -i Q1_public/kplr000757076-2009166043257_llc.fits \
    -readformat 0 1 10 11 \
    -fluxtomag 25.0 0 \
    -o kplr000757076-2009166043257_llc.asc.txt
```

---

### `-match`

**Syntax**
```
-match
    < "file" filename | "inlist" inlistcolumn >
    ["opencommand" command] ["skipnum" Nskip]
    ["skipchar" <skipchar1[,skipchar2,...]>] ["delimiter" delimiter]
    < "matchcolumn" < varname:column | colnum > >
    < "addcolumns" varname1:colnum1[:coltype1[:colformat1]]
        [,varname2:colnum2[:coltype2[:colformat2]],...] >
    < "cullmissing" | "nanmissing" | "missingval" value >
```

**Description**

Perform a row-by-row match of an external data file to the light curve, merging columns from the file into the light curve. The match is performed on a specified variable (by default the time `t`). Float/double columns are matched within the tolerance set by `-jdtol`; all other types are matched exactly.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `"file" filename` | A single file matched to every light curve. Files ending in `.fits` are treated as binary FITS tables. |
| `"inlist" col` | Read the name of the match file for each light curve from column `col` of the input list. |
| `"opencommand" command` | Shell command to pre-process each match file; `%s` is replaced with the filename. Output is read from the command's stdout. |
| `"skipnum" Nskip` | Number of lines to skip at the top of each match file. |
| `"skipchar" chars` | Comma-separated list of comment characters (default: `#`). |
| `"delimiter" delim` | Column delimiter (default: whitespace). |
| `"matchcolumn" varname:colnum` | Light-curve variable and match-file column number (or FITS column name) to use as the match key. |
| `"addcolumns" varname:colnum[:type[:fmt]],...` | Columns to merge from the match file into the light curve. New variables are created; existing ones are overwritten. |
| `"cullmissing"` | Remove light curve rows that have no match. |
| `"nanmissing"` | Set unmatched rows to NaN (float/double) or 0 (integer/string). |
| `"missingval" value` | Set unmatched rows to the specified value. |

---

### `-Phase`

**Syntax**
```
-Phase
    < "aov" | "ls" | "bls" | "fixcolumn" <colname | colnum>
      | "list" ["column" col] | "fix" P >
    ["T0" < "bls" phaseTc | "fixcolumn" <colname | colnum>
           | "list" ["column" col] | "fix" T0 >]
    ["phasevar" var] ["startphase" startphase]
```

**Description**

Replace the time axis of the light curve with its phase and sort by phase. The phase is defined as `((t - T0) mod P) / P`, running from 0 to 1 (or from `startphase` to `startphase + 1`). After `-Phase`, time-based binning commands like `-binlc` operate in phase space.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `"aov"` | Use the best period from the most recent `-aov` or `-aov_harm` command. |
| `"ls"` | Use the best period from the most recent `-LS` command. |
| `"bls"` | Use the best period from the most recent `-BLS` command. |
| `"fixcolumn" col` | Use the value of a previously computed output column (by name or number). |
| `"list"` | Read the period from the input list (optionally specify the `"column"`). |
| `"fix" P` | Fix the period to the value `P` for all light curves. |
| `"T0" ...` | Reference epoch for phase zero. `"bls" phaseTc` assigns the BLS mid-transit time the phase `phaseTc`; other sources follow standard vartools conventions. Default: the earliest time in the light curve. |
| `"phasevar" var` | Store phases in the variable `var` instead of overwriting `t`. |
| `"startphase" startphase` | Start of the phase range (default 0; range becomes `[startphase, startphase+1)`). |

**Examples**

**Example 1.** Phase-fold the light curve on a fixed period of 1.2354 days and write the result to disk.

```bash
vartools -i EXAMPLES/2 -header \
    -Phase fix 1.2354 \
    -o EXAMPLES/OUTDIR1/2.phase.txt
```

**Example 2.** Detect a transit with BLS, phase-fold on the BLS period placing mid-transit at phase 0.5, write the phase-folded light curve, then median-bin it into 200 phase bins.

```bash
vartools -i EXAMPLES/3.transit -oneline \
    -BLS q 0.01 0.1 0.5 5.0 20000 200 7 1 0 0 0 \
    -Phase bls T0 bls 0.5 \
    -o EXAMPLES/OUTDIR1/3.phase.txt \
    -binlc median nbins 200 tcenter \
    -o EXAMPLES/OUTDIR1/3.phasebin.txt
```

---

### `-resample`

**Syntax**
```
-resample
    < "nearest" | "linear" | "spline" ["left" yp1] ["right" ypn]
      | "splinemonotonic" | "bspline" ["nbreaks" nbreaks] ["order" order] >
    ["file" < "fix" times_file ["column" col] |
              "list" ["listcolumn" col] ["tcolumn" col] > |
     ["tstart" < "fix" val | "fixcolumn" col | "list" | "expr" expr >]
     ["tstop"  < "fix" val | "fixcolumn" col | "list" | "expr" expr >]
     ["delt"   < "fix" val | "fixcolumn" col | "list" | "expr" expr > |
      "Npoints" < "fix" val | "fixcolumn" col | "list" | "expr" expr >]]
    ["gaps" <separation spec> <method>]
    ["extrap" <method>]
```

**Description**

Resample the light curve onto a new time base by interpolating all light curve vectors. The default output grid runs from the first to the last observed time with a step equal to the minimum observed time separation. String-type columns (e.g. image IDs) are always resampled with the `"nearest"` method.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `"nearest"` | Nearest-neighbor (step-function) interpolation. |
| `"linear"` | Linear interpolation. |
| `"spline"` | Cubic spline interpolation. Optional `"left"`/`"right"` set boundary first-derivative conditions. |
| `"splinemonotonic"` | Cubic spline constrained to be monotonic between observations. |
| `"bspline"` | Basis-spline interpolation. Optional `"nbreaks"` (default 15) and `"order"` (default 3). If `nbreaks < 2`, breaks are increased until χ²/dof ≤ 1 (can be slow). |
| `"file" ...` | Specify an arbitrary (non-uniform) time base from a file. |
| `"tstart"` / `"tstop"` / `"delt"` / `"Npoints"` | Override the start, stop, step size, or total number of resampled points. Each accepts `"fix"`, `"fixcolumn"`, `"list"`, or `"expr"`. |
| `"gaps" <sep> <method>` | Use a different interpolation method for points that are farther than `sep` from any observed time. Separation can be fixed, taken from a column, or derived from the data distribution (`"frac_min_sep"`, `"frac_med_sep"`, `"percentile_sep"`). |
| `"extrap" <method>` | Use a different method for extrapolated points (beyond the observed time range). |

**Examples**

**Example 1.** Resample using linear interpolation with default time sampling.

```bash
vartools -i EXAMPLES/2 -resample linear \
    -o EXAMPLES/2.resample.example1
```

**Example 2.** Resample using monotonic spline between fixed start/stop times with exactly 1000 output points.

```bash
vartools -i EXAMPLES/2 -resample splinemonotonic \
    tstart fix 53726 tstop fix 53756 Npoints fix 1000 \
    -o EXAMPLES/2.resample.example2
```

**Example 3.** Resample `EXAMPLES/4` onto the time grid of `EXAMPLES/8` using linear interpolation.

```bash
vartools -i EXAMPLES/4 -resample linear \
    file fix EXAMPLES/8 \
    -o EXAMPLES/4.resample.example3
```

**Example 4.** Resample onto a uniform 0.001-day grid, using a B-spline for normal gaps, nearest-neighbour for large gaps (> 80th-percentile separation), and nearest-neighbour for extrapolation.

```bash
vartools -i EXAMPLES/1 -resample splinemonotonic \
    tstart fix 53725 tstop fix 53757 delt fix 0.001 \
    gaps percentile_sep 80 bspline nbreaks 15 order 3 \
    extrap nearest \
    -o EXAMPLES/1.resample.example4
```

---

### `-rescalesig`

**Syntax**
```
-rescalesig
    ["maskpoints" maskvar]
```

**Description**

Rescale the magnitude uncertainties of each light curve independently so that χ²/dof = 1 for that light curve. The rescale factor applied to each light curve is included in the output table. Unlike `-ensemblerescalesig`, this operates per-light-curve without fitting an ensemble relation.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `"maskpoints" maskvar` | Restrict the χ²/dof calculation to points with `maskvar > 0`. |

**Examples**

**Example 1.** Rescale the formal errors in a light curve so that χ²/dof equals 1, demonstrating the effect with `-chi2` before and after.

```bash
vartools -i EXAMPLES/4 -oneline \
    -chi2 -rescalesig -chi2
```

Output:
```
Name                  = EXAMPLES/4
Chi2_0                = 5.19874
Weighted_Mean_Mag_0   = 10.35137
SigmaRescaleFactor_1  = 0.43858
Chi2_2                = 1.00000
Weighted_Mean_Mag_2   = 10.35137
```

---

### `-sortlc`

**Syntax**
```
-sortlc
    ["var" varname] ["reverse"]
```

**Description**

Sort the light curve. By default it is sorted by time. If any subsequent command requires time-sorted data and the light curve was re-sorted by another variable, it will be automatically re-sorted by time at the start of that command.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `"var" varname` | Sort by this named variable instead of time. |
| `"reverse"` | Sort in descending order. |

---

### `-restricttimes`

**Syntax**
```
-restricttimes
    ["exclude"]
    < "JDrange" minJD maxJD |
      "JDrangebylc"
          < "fix" minJD | "list" [...] | "fixcolumn" col | "expr" expr >
          < "fix" maxJD | "list" [...] | "fixcolumn" col | "expr" expr >
      "JDlist" JDfilename |
      "imagelist" imagefilename |
      "expr" expression >
    ["markrestrict" markvar ["noinitmark"]]
```

**Description**

Filter observations from the light curve based on their time values (or string IDs). By default, points **not** matching the specification are removed. If `"exclude"` is given, the **matching** points are removed instead.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `"exclude"` | Invert the filter: remove the specified times instead of keeping them. |
| `"JDrange" minJD maxJD` | Keep points with `minJD ≤ t ≤ maxJD`. The same range applies to all light curves. |
| `"JDrangebylc" ...` | Per-light-curve range. Each bound can be `"fix"`, `"list"`, `"fixcolumn"`, or `"expr"`. |
| `"JDlist" file` | File with JDs (first column) to keep (or exclude). |
| `"imagelist" file` | File with string image IDs (first column) to keep (or exclude). |
| `"expr" expression` | Keep points where the expression evaluates to a value > 0. For example, `'(mag>9.0)&&(mag<9.5)'`. |
| `"markrestrict" markvar` | Do not remove points; instead set `markvar = 1` for kept points and `markvar = 0` for filtered points. |
| `"noinitmark"` | Combine with `"markrestrict"`: treat the existing values of `markvar` as a prior mask and only update points that would be newly removed. |

**Examples**

**Example 1.** Keep only observations between JD 53740 and 53750 and verify with `-stats`.

```bash
vartools -i EXAMPLES/3 -stats t min,max \
    -restricttimes JDrange 53740 53750 \
    -stats t min,max -oneline
```

Output:
```
Name                  = EXAMPLES/3
STATS_t_MIN_0         = 53725.173920000001
STATS_t_MAX_0         = 53756.281021000003
RestrictTimes_MinJD_1 = 53740
RestrictTimes_MaxJD_1 = 53750
STATS_t_MIN_2         = 53740.336210000001
STATS_t_MAX_2         = 53745.478681000001
```

**Example 2.** Keep only points where `10.16311 < mag < 10.17027` using an expression filter.

```bash
vartools -i EXAMPLES/3 -stats mag min,max \
    -restricttimes expr '(mag>10.16311)&&(mag<10.17027)' \
    -stats mag min,max -oneline
```

**Example 3.** Compute the 20th and 80th percentile magnitudes, then restrict to that range dynamically using computed column names.

```bash
vartools -i EXAMPLES/3 -stats mag pct20.0,pct80.0 \
    -restricttimes expr \
      '(mag>STATS_mag_PCT20_00_0)&&(mag<STATS_mag_PCT80_00_00)' \
    -stats mag min,max -oneline
```

**Example 4.** Detect a transit with BLS, phase-fold, and *exclude* the in-transit points from the output light curve.

```bash
vartools -i EXAMPLES/3.transit -oneline \
    -BLS q 0.01 0.1 0.1 20.0 100000 200 0 1 \
         0 0 0 fittrap nobinnedrms \
    -expr 'minph=0.5-BLS_Qtran_1_0/2.0' \
    -expr 'maxph=0.5+BLS_Qtran_1_0/2.0' \
    -expr 'ph=t' \
    -changevariable t ph \
    -Phase bls T0 bls 0.5 \
    -restricttimes exclude JDrangebylc expr minph expr maxph \
    -o EXAMPLES/3.cliptransit columnformat t,mag,err,ph
```

---

### `-restoretimes`

**Syntax**
```
-restoretimes
    prior_restricttimes_command
```

**Description**

Restore observations that were removed by a prior `-restricttimes` command. The restored points are appended to the current light curve and then sorted by time. This allows commands to be applied to isolated time windows.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `prior_restricttimes_command` | Integer index (1-based) identifying which `-restricttimes` command to restore from. `1` = first, `2` = second, etc. Cannot be used if that `-restricttimes` used the `"markrestrict"` option. |

**Examples**

**Example 1.** Compute RMS on the full light curve, then restrict to a time window, compute RMS again, then restore and confirm the full dataset is recovered.

```bash
vartools -i EXAMPLES/3 \
    -rms \
    -restricttimes JDrange 53740 53750 \
    -rms \
    -restoretimes 1 \
    -rms -oneline
```

The full dataset contains 3417 points; the restricted window 738 points; after `-restoretimes 1` all 3417 points are recovered.

**Example 2.** Shift magnitudes only within the restricted window, then restore — points outside the window retain their original values.

```bash
vartools -i EXAMPLES/3 \
    -restricttimes JDrange 53740 53750 \
    -expr 'mag=mag+0.05' \
    -restoretimes 1 \
    -o EXAMPLES/OUTDIR1/3.restoretimes.txt -oneline
```

---

### `-copylc`

**Syntax**
```
-copylc
    Ncopies
```

**Description**

Replicate the current light curve `Ncopies` times. Each copy is processed independently through all subsequent commands. Output rows for commands preceding the `-copylc` call are also replicated. Copy names are appended with the suffix `_copy<cmdnum>.<copynum>`. Cannot be used with the `-readall` option.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `Ncopies` | Number of copies to create (integer ≥ 1). |

**Examples**

**Example 1.** Run the LS periodogram on a light curve, then make 100 copies with Gaussian noise substituted for the magnitudes and run LS on each simulation to assess false-alarm rates.

```bash
vartools -i EXAMPLES/2 -LS 0.1 10. 0.1 1 0 \
    -copylc 100 \
    -expr 'mag=err*gauss()' \
    -LS 0.1 10. 0.1 1 0 \
    -header
```

Output:
```
#Name LS_Period_1_0 Log10_LS_Prob_1_0 LS_SNR_1_0 LS_Period_1_3 Log10_LS_Prob_1_3 LS_SNR_1_3
EXAMPLES/2     1.23440877 -704.49194   58.45119     0.16051136   -0.65620    5.59463
EXAMPLES/2_copy1.1     1.23440877 -704.49194   58.45119     0.20478671   -0.03436    4.89785
EXAMPLES/2_copy1.2     1.23440877 -704.49194   58.45119     0.17224308   -0.00000    4.62340
EXAMPLES/2_copy1.3     1.23440877 -704.49194   58.45119     0.43264396   -0.13495    6.51119
EXAMPLES/2_copy1.4     1.23440877 -704.49194   58.45119     0.12629761   -1.21929    7.96670
EXAMPLES/2_copy1.5     1.23440877 -704.49194   58.45119     6.48064604   -0.44006    6.40004
...
```

---

### `-medianfilter`

**Syntax**
```
-medianfilter
    time ["average" | "weightedaverage"] ["replace"]
```

**Description**

Apply a running median (or mean) filter to the light curve. By default a **high-pass** filter is applied: the running median of all points within ±`time` days of each observation is subtracted from that observation. If `"replace"` is given, the command acts as a **low-pass** filter: each point is replaced by the running median (or mean) rather than having it subtracted.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `time` | Half-width of the running window in days. |
| `"average"` | Use the running mean instead of the median. |
| `"weightedaverage"` | Use the error-weighted running mean. |
| `"replace"` | Low-pass mode: replace each point with the running statistic. |

---

### `-decorr`

**Syntax**
```
-decorr
    correctlc zeropointterm subtractfirstterm
    Nglobalterms globalfile1 order1 ... globalfileN orderN
    Nlcterms lccolumn1 lcorder1 ... lccolumnN lcorderN
    omodel [modeloutdir] ["maskpoints" maskvar]
```

!!! warning "Deprecated"
    As of VARTOOLS 1.3, `-decorr` is deprecated. Use `-linfit` instead for general linear decorrelation.

**Description**

Decorrelate light curves against external signals. Signals may be global (provided in separate files with columns `JD signal_value`) or light-curve-specific (additional columns in the light curve files).

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `correctlc` | `1` to subtract the fitted model from the light curve; `0` to output statistics only. |
| `zeropointterm` | `1` to include a constant offset term in the fit. |
| `subtractfirstterm` | `1` to decorrelate against `signal - signal[0]` rather than `signal` (avoids round-off when using JD). |
| `Nglobalterms` | Number of global signal files. |
| `globalfile1 ... globalfileN` | Names of global signal files. |
| `order1 ... orderN` | Polynomial order for each global signal (must be ≥ 1). |
| `Nlcterms` | Number of light-curve-specific signals. |
| `lccolumn1 ... lccolumnN` | Column indices in the light curve for each signal. |
| `lcorder1 ... lcorderN` | Polynomial order for each light-curve signal. |
| `omodel` | `1` to output the model; `0` otherwise. |
| `modeloutdir` | Directory for model output files (suffix: `.decorr.model`). |
| `"maskpoints" maskvar` | Restrict the fit to points with `maskvar > 0`. |

---

### `-print`

**Syntax**
```
-print
    var1[,var2,var3...]
    ["columnnames" col1[,col2,col3...]]
    ["format" fmt1[,fmt2,fmt3...]]
```

**Description**

Print the value of one or more variables to the output statistics table. One value per light curve is output. If a variable is a light-curve vector (per-point), only the value at the first time point is reported. This is useful for passing computed values (e.g. parameters from `-linfit` or custom scalars from `-expr`) into the output table.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `var1,...` | Comma-separated list of variable names to print. |
| `"columnnames" col1,...` | Override the default column names (default pattern: `Print_<varname>_<varnum>_<cmdnum>`). The command index is still appended unless `-columnsuffix` is used. |
| `"format" fmt1,...` | Printf-style format strings for each variable (e.g. `%.6f`). |

**Example**

```bash
# Compute a custom statistic and include it in the output
vartools -l lc_list \
    -stats mag,mean \
    -expr 'delta=mag-Stats_mag_mean_0' \
    -stats delta,stddev \
    -print Stats_mag_mean_0,Stats_delta_stddev_1 \
    -header -tab
```
