# Period Finding

This page documents the VARTOOLS commands for detecting and characterizing periodic signals in light curves. Commands output period, significance, and signal-to-noise statistics to the output table; periodograms are optionally written to on-disk files.

---

### `-LS` — Generalized Lomb-Scargle

**Syntax**
```
-LS
    < "var" minpvar | "expr" minpexpr | minp >
    < "var" maxpvar | "expr" maxpexpr | maxp >
    < "var" subsamplevar | "expr" subsampleexpr | subsample >
    Npeaks operiodogram [outdir] ["noGLS"]
    ["whiten"] ["clip" clip clipiter]
    ["fixperiodSNR" < "aov" | "ls" | "injectharm" | "fix" period
                    | "list" ["column" col]
                    | "fixcolumn" <colname | colnum> >]
    ["bootstrap" Nbootstrap] ["maskpoints" maskvar]
```

**Description**

Perform a Generalized Lomb-Scargle (GLS) period search for sinusoidal signals. The search runs over frequencies from `fmin = 1/maxp` to `fmax = 1/minp` with a uniform frequency step `Δf = subsample/T`, where `T` is the time baseline. The GLS implementation of Zechmeister and Kürster (2009) allows a floating mean and heteroscedastic errors, unlike the traditional LS periodogram.

For each of the three search parameters (`minp`, `maxp`, `subsample`) you may either give a fixed value on the command line, use the `"var"` keyword followed by a variable name, or use the `"expr"` keyword to evaluate an analytic expression for each light curve.

The statistic reported is:

```
LS = (chi0^2 - chi(f)^2) / chi0^2
```

where `chi0^2` is χ² about the weighted mean and `chi(f)^2` is χ² about the best-fit sinusoid at frequency `f`. For the traditional (non-generalized) periodogram, `LS` is the standard un-normalized Lomb-Scargle power.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `minp` | Minimum period to search (days). |
| `maxp` | Maximum period to search (days). |
| `subsample` | Frequency oversampling factor. Typical value: 0.1 (one-tenth of the fundamental frequency resolution). |
| `Npeaks` | Number of highest peaks to find and report. |
| `operiodogram` | `1` to write the periodogram to `outdir/$basename.ls`; `0` to skip. |
| `outdir` | Directory for periodogram files (required when `operiodogram = 1`). |
| `"noGLS"` | Compute the traditional (non-generalized) Lomb-Scargle periodogram. |
| `"whiten"` | After each peak, whiten the light curve at that period before searching for the next peak. The SNR for each peak is computed on the whitened periodogram. |
| `"clip" clip clipiter` | Sigma-clipping parameters for the mean/RMS calculation used in the SNR estimate. `clip` is the clipping factor; `clipiter` is `1` for iterative clipping (default: iterative 5σ). |
| `"fixperiodSNR" ...` | Also output `log(FAP)` and SNR at a specified period. Sources: `"aov"` (last `-aov`), `"ls"` (last `-LS`), `"injectharm"` (last `-Injectharm`), `"fix" period`, `"list"`, or `"fixcolumn"`. |
| `"bootstrap" Nbootstrap` | Estimate the FAP via bootstrap resampling (`Nbootstrap` simulations). Each simulation uses the observed times with magnitudes drawn randomly with replacement. |
| `"maskpoints" maskvar` | Exclude points with `maskvar ≤ 0` from the periodogram calculation. |

**Output columns** (per peak `k`, command index `i`)

| Column | Description |
|--------|-------------|
| `LS_Period_k_i` | Best period of peak `k` (days). |
| `Log10_LS_Prob_k_i` | Log₁₀ of the formal false alarm probability. |
| `LS_SNR_k_i` | Spectroscopic SNR: `(LS - <LS>) / RMS(LS)`. |

**References**

Cite Zechmeister & Kürster 2009, A&A, 496, 577 and Press et al. 1992 (Numerical Recipes) for the GLS periodogram. For the traditional LS periodogram also cite Lomb 1976, Scargle 1982, and Press & Rybicki 1989.

**Examples**

**Example 1.** Run the Lomb-Scargle period search on a single light curve, reporting the top 5 peaks, writing periodograms, applying whitening, and using iterative 5σ clipping for SNR.

```bash
vartools -i EXAMPLES/2 -oneline -ascii \
    -LS 0.1 10. 0.1 5 1 EXAMPLES/OUTDIR1 whiten clip 5. 1
```

Output:
```
LS_Period_1_0 = 1.23440877 (Log10_LS_Prob = -704.49194, LS_SNR = 126.52865)
LS_Period_2_0 = 0.54573861 (Log10_LS_Prob = -60.19769, LS_SNR = 12.22747)
LS_Period_3_0 = 0.55647766 (Log10_LS_Prob = -52.64556, LS_SNR = 9.46714)
LS_Period_4_0 = 0.25922584 (Log10_LS_Prob = -28.26507, LS_SNR = 7.42006)
LS_Period_5_0 = 0.50172744 (Log10_LS_Prob = -25.89142, LS_SNR = 8.63338)
```

---

### `-aov` — Phase-Binned Analysis of Variance

**Syntax**
```
-aov
    ["Nbin" < "var" Nbinvar | "expr" Nbinexpr | Nbin >]
    < "var" minpvar | "expr" minpexpr | minp >
    < "var" maxpvar | "expr" maxpexpr | maxp >
    < "var" subsamplevar | "expr" subsampleexpr | subsample >
    < "var" finetunevar | "expr" finetuneexpr | finetune >
    Npeaks operiodogram [outdir]
    ["whiten"] ["clip" clip clipiter] ["uselog"]
    ["fixperiodSNR" < "aov" | "ls" | "injectharm" | "fix" period
                    | "list" ["column" col]
                    | "fixcolumn" <colname | colnum> >]
    ["maskpoints" maskvar]
```

**Description**

Perform an Analysis of Variance (AoV) period search using phase binning. For each trial frequency, the light curve is phase-folded and binned; the AoV statistic θ_aov measures how much variance is explained by the phase bins relative to the total variance. A high θ_aov indicates a good phase-coherent signal.

The initial search uses a frequency step of `subsample/T`. The top peaks are refined to a resolution of `finetune/T`.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `"Nbin"` | Number of phase bins (default: 8). |
| `minp` / `maxp` | Period search range (days). Accepts `"var"` or `"expr"` keywords. |
| `subsample` | Coarse frequency step factor. |
| `finetune` | Fine-tuning frequency step factor applied near peak periods. |
| `Npeaks` | Number of peaks to report. |
| `operiodogram` | `1` to write period vs. θ_aov to `outdir/$basename.aov`. |
| `"whiten"` | Whiten the light curve at each peak before searching for the next one. |
| `"clip" clip clipiter` | Clipping parameters for the SNR calculation. |
| `"uselog"` | Output `ln(θ_aov)` SNR: `(<-ln(θ_aov)> - ln(θ_aov)) / RMS(-ln(θ_aov))`. Also outputs the mean and RMS of `-ln(θ_aov)`. |
| `"fixperiodSNR" ...` | Output AoV statistic and SNR at a specified period. Syntax identical to `-LS`. |
| `"maskpoints" maskvar` | Exclude points with `maskvar ≤ 0`. |

**Output columns** (per peak `k`, command index `i`)

| Column | Description |
|--------|-------------|
| `Period_k_i` | Best period of peak `k` (days). |
| `AOV_k_i` | θ_aov statistic. |
| `AOV_SNR_k_i` | Signal-to-noise ratio in the periodogram. |
| `AOV_NEG_LOG_FAP_k_i` | `-ln(FAP)` (formal false alarm probability). |

**References**

Cite Schwarzenberg-Czerny 1989, MNRAS, 241, 153 and Devor 2005, ApJ, 628, 411.

**Examples**

**Example 1.** Run the phase-binning AoV period search on a light curve, using 20 phase bins, searching periods between 0.1 and 10.0 days, reporting the top 5 peaks with iterative 5σ clipping.

```bash
vartools -i EXAMPLES/2 -oneline -ascii \
    -aov Nbin 20 0.1 10. 0.1 0.01 5 1 EXAMPLES/OUTDIR1 whiten clip 5. 1
```

Output: Five detected periods with `Period_1_0 = 1.23583047` (`AOV_1_0 = 18330.55450`, AOV_SNR, AOV_NEG_LN_FAP) through `Period_5_0 = 0.12056969` (`AOV_5_0 = 16.69317`), with corresponding SNR and negative log FAP values for each.

---

### `-aov_harm` — Multi-Harmonic Analysis of Variance

**Syntax**
```
-aov_harm
    < "var" Nharmvar | "expr" Nharmexpr | Nharm >
    < "var" minpvar | "expr" minpexpr | minp >
    < "var" maxpvar | "expr" maxpexpr | maxp >
    < "var" subsamplevar | "expr" subsampleexpr | subsample >
    < "var" finetunevar | "expr" finetuneexpr | finetune >
    Npeaks operiodogram [outdir]
    ["whiten"] ["clip" clip clipiter]
    ["fixperiodSNR" < "aov" | "ls" | "injectharm" | "fix" period
                    | "list" ["column" col]
                    | "fixcolumn" <colname | colnum> >]
    ["maskpoints" maskvar]
```

**Description**

Perform an AoV period search fitting a multi-harmonic model instead of phase bins. The model signal has `Nharm` harmonics. If `Nharm < 1`, the number of harmonics is automatically varied to minimize the false alarm probability (accounting for the penalty for overfitting). All other parameters are identical to `-aov`.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `Nharm` | Number of harmonics in the model signal (≥ 1). Set to `0` or negative to allow automatic selection. Accepts `"var"` or `"expr"` keywords. |
| All others | Same as `-aov`. |

**Output columns**

Same structure as `-aov` with prefix `AOV_HARM`.

**References**

Cite Schwarzenberg-Czerny 1996, ApJ, 460, L107.

**Examples**

**Example 1.** Run the multi-harmonic AoV period search with 1 harmonic, searching periods between 0.1 and 10.0 days, reporting the top 2 peaks with whitening and iterative 5σ clipping.

```bash
vartools -i EXAMPLES/2 -oneline -ascii \
    -aov_harm 1 0.1 10. 0.1 0.01 2 1 EXAMPLES/OUTDIR1 whiten clip 5. 1
```

Output: Period values and AOV_HARM, SNR, and logarithmic FAP values for 2 identified peaks (1.23533969 days and 0.49981672 days).

---

### `-BLS` — Box-fitting Least Squares

**Syntax**
```
-BLS
    < "r" < "var" rminvar | "expr" rminexpr | rmin >
              < "var" rmaxvar | "expr" rmaxexpr | rmax > |
      "q" < "var" qminvar | "expr" qminexpr | qmin >
          < "var" qmaxvar | "expr" qmaxexpr | qmax > |
      "density" < "var" rhovar | ... | rho >
                < "var" minexpdurfracvar | ... | minexpdurfrac >
                < "var" maxexpdurfracvar | ... | maxexpdurfrac > >
    < "var" minpervar | "expr" minperexpr | minper >
    < "var" maxpervar | "expr" maxperexpr | maxper >
    < "nf" < ... | nfreq > |
      "df" < ... | dfreq > |
      "optimal" < ... | subsample > >
    < "var" nbinsvar | "expr" nbinsexpr | nbins >
    timezone Npeak outperiodogram [outdir] omodel [modeloutdir]
    correctlc ["extraparams"] ["fittrap"] ["nobinnedrms"]
    ["ophcurve" outdir phmin phmax phstep]
    ["ojdcurve" outdir jdstep]
    ["stepP" | "steplogP"]
    ["adjust-qmin-by-mindt" ["reduce-nbins"]]
    ["reportharmonics"] ["maskpoints" maskvar]
```

**Description**

Run the Box-Least Squares (BLS) transit search algorithm (Kovács, Zucker & Mazeh 2002). BLS searches for periodic box-shaped (trapezoidal) dips consistent with a transiting companion. The search is performed over a grid of trial periods and phase bins.

Three ways to specify the allowed range of transit durations:

- `"q" qmin qmax` — minimum and maximum fractional transit duration (fraction of the orbit spent in transit).
- `"r" rmin rmax` — minimum and maximum stellar radius in solar radii. The q range for each period is derived from `q = 0.076·R^(2/3)·P^(-2/3)`.
- `"density" rho minexpdurfrac maxexpdurfrac` — stellar density (g/cm³) with minimum and maximum fractions of the expected circular-orbit transit duration.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `minper` / `maxper` | Period search range in days. |
| `"nf" nfreq` | Total number of trial frequencies. Rule of thumb: `nfreq ≈ T · (fmax - fmin) / (0.25·qmin)`. |
| `"df" dfreq` | Explicit frequency step size. |
| `"optimal" subsample` | Use Ofir (2014) optimal frequency spacing (requires `"density"` mode). |
| `nbins` | Number of phase bins (≥ `2/qmin`). |
| `timezone` | Hours to add to UTC to obtain local time; used to compute the fraction of Δχ² from a single night. |
| `Npeak` | Number of peaks to find and report. |
| `outperiodogram` | `1` to output the BLS spectrum (SN vs. period) to `outdir/$basename.bls`. |
| `omodel` | `1` to output the best-fit box model to `modeloutdir/$basename.bls.model`. |
| `correctlc` | `1` to subtract the transit model before passing to the next command. |
| `"extraparams"` | Include additional false-positive diagnostic parameters in the output. |
| `"fittrap"` | Fit a trapezoidal transit at each BLS peak. Adds `qingress` (fraction of transit duration in ingress) and `OOTmag` to output. |
| `"nobinnedrms"` | Compute BLS_SN without binned RMS (faster, but SN is suppressed for high-significance peaks). |
| `"ophcurve" outdir phmin phmax phstep` | Output a model phase curve to `outdir/$basename.bls.phcurve`. |
| `"ojdcurve" outdir jdstep` | Output a model light curve (JD space) to `outdir/$basename.bls.jdcurve`. |
| `"stepP"` / `"steplogP"` | Sample the BLS spectrum at uniform P or log(P) steps instead of uniform frequency steps. |
| `"adjust-qmin-by-mindt"` | Adaptively increase `qmin` at each frequency to `max(qmin, mindt·f)`. |
| `"reduce-nbins"` | (With `adjust-qmin-by-mindt`) adaptively reduce `nbins` at each frequency. |
| `"reportharmonics"` | Report period harmonics even if a higher-power peak at a multiple of that frequency exists. |
| `"maskpoints" maskvar` | Exclude points with `maskvar ≤ 0` from the BLS spectrum. |

**Output columns** (per peak `k`, command index `i`)

| Column | Description |
|--------|-------------|
| `BLS_Period_k_i` | Best period of peak `k` (days). |
| `BLS_Tc_k_i` | Mid-transit epoch. |
| `BLS_SN_k_i` | Signal-to-noise ratio in the BLS spectrum. |
| `BLS_SR_k_i` | BLS spectral residual. |
| `BLS_SDE_k_i` | Signal detection efficiency. |
| `BLS_Depth_k_i` | Transit depth (magnitudes). |
| `BLS_Qtran_k_i` | Fractional transit duration `q`. |
| `BLS_deltaChi2_k_i` | Δχ² of the best-fit transit. |
| `BLS_SignaltoPinknoise_k_i` | Signal-to-pink-noise ratio. |
| `BLS_Ntransits_k_i` | Number of observed transits. |
| `BLS_Npointsintransit_k_i` | Points within the transit window. |
| `BLS_fraconenight_k_i` | Fraction of Δχ² from a single night. |
| `BLS_Rednoise_k_i` | Estimated red noise level. |
| `BLS_Whitenoise_k_i` | Estimated white noise level. |

When `"fittrap"` is given, `BLS_Qingress_k_i` and `BLS_OOTmag_k_i` are also included.

**References**

Cite Kovács, Zucker & Mazeh 2002, A&A, 391, 369. For the optimal frequency sampling cite Ofir 2014, A&A, 561, A138.

**Examples**

**Example 1.** Apply the BLS transit search to a light curve with an injected transit, scanning fractional durations 0.01–0.1, searching periods 0.1–20.0 days with 100,000 frequencies and 200 phase bins, fitting a trapezoid model and outputting a phase curve.

```bash
vartools -i EXAMPLES/3.transit -ascii -oneline \
    -BLS q 0.01 0.1 0.1 20.0 100000 200 0 1 \
        1 EXAMPLES/OUTDIR1/ 1 EXAMPLES/OUTDIR1/ 0 fittrap \
        nobinnedrms ophcurve EXAMPLES/OUTDIR1/ -0.1 1.1 0.001
```

Output:
```
Name                         = EXAMPLES/3.transit
BLS_Period_1_0               =     2.12334706
BLS_Tc_1_0                   = 53727.297293937358
BLS_SN_1_0                   =   7.26127
BLS_SR_1_0                   =   0.00238
BLS_SDE_1_0                  =   6.34195
BLS_Depth_1_0                =   0.01220
BLS_Qtran_1_0                =   0.03576
BLS_Qingress_1_0             =   0.19618
BLS_OOTmag_1_0               =  10.16686
BLS_i1_1_0                   =   0.98213
BLS_i2_1_0                   =   1.01790
BLS_deltaChi2_1_0            = -24217.21939
BLS_fraconenight_1_0         =   0.43155
BLS_Npointsintransit_1_0     =   165
BLS_Ntransits_1_0            =     4
BLS_Npointsbeforetransit_1_0 =   127
BLS_Npointsaftertransit_1_0  =   143
BLS_Rednoise_1_0             =   0.00151
BLS_Whitenoise_1_0           =   0.00489
BLS_SignaltoPinknoise_1_0    =  14.38935
BLS_Period_invtransit_0      =     1.14594782
BLS_deltaChi2_invtransit_0   = -3301.69183
BLS_MeanMag_0                =  10.16740
```

---

### `-BLSFixPer` — BLS at a Fixed Period

**Syntax**
```
-BLSFixPer
    < "aov" | "ls" | "list" ["column" col]
      | "fix" period | "fixcolumn" <colname | colnum>
      | "expr" expr >
    < "r" rmin rmax | "q" qmin qmax >
    nbins timezone omodel [model_outdir] correctlc ["fittrap"]
    ["maskpoints" maskvar]
```

**Description**

Run BLS at a single fixed period, searching only for the most transit-like signal at that period. Useful as a second pass after a full BLS or LS period search.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| Period source | `"aov"` (last `-aov`), `"ls"` (last `-LS`), `"list"`, `"fix" period`, `"fixcolumn"`, or `"expr"`. |
| `"r" rmin rmax` / `"q" qmin qmax` | Transit duration range (stellar radius bounds or fractional duration bounds). |
| `nbins` | Number of phase bins. |
| `timezone` | UTC offset for single-night fraction calculation. |
| `omodel` | `1` to output the model to `model_outdir` (suffix: `.blsfixper.model`). |
| `correctlc` | `1` to subtract the model from the light curve. |
| `"fittrap"` | Fit a trapezoid transit at the BLS peak. |
| `"maskpoints" maskvar` | Exclude points with `maskvar ≤ 0`. |

**References**

Cite Kovács, Zucker & Mazeh 2002, A&A, 391, 369.

**Examples**

**Example 1.** Fit a box-shaped (with trapezoid fitting) transit model at a fixed period of 2.12345 days, scanning fractional durations 0.01–0.1 with 200 phase bins, and compute RMS before and after subtracting the transit.

```bash
vartools -i EXAMPLES/3.transit -ascii -oneline \
    -rms \
    -BLSFixPer fix 2.12345 q 0.01 0.1 200 0 0 1 fittrap \
    -rms
```

Output:
```
Name                             = EXAMPLES/3.transit
Mean_Mag_0                       =  10.16727
RMS_0                            =   0.00542
Expected_RMS_0                   =   0.00104
Npoints_0                        =  3417
BLSFixPer_Period_1               =  2.12345000
BLSFixPer_Tc_1                   = 53727.29676321477
BLSFixPer_SR_1                   =   0.00238
BLSFixPer_Depth_1                =   0.01189
BLSFixPer_Qtran_1                =   0.03626
BLSFixPer_Qingress_1             =   0.20623
BLSFixPer_OOTmag_1               =  10.16687
BLSFixPer_i1_1                   =   0.98158
BLSFixPer_i2_1                   =   1.01785
BLSFixPer_deltaChi2_1            = -24228.56603
BLSFixPer_fraconenight_1         =   0.43366
BLSFixPer_Npointsintransit_1     =   166
BLSFixPer_Ntransits_1            =     4
BLSFixPer_Npointsbeforetransit_1 =   129
BLSFixPer_Npointsaftertransit_1  =   144
BLSFixPer_Rednoise_1             =   0.00151
BLSFixPer_Whitenoise_1           =   0.00489
BLSFixPer_SignaltoPinknoise_1    =  14.08946
BLSFixPer_deltaChi2_invtransit_1 = -2934.30109
BLSFixPer_MeanMag_1              =  10.16740
Mean_Mag_2                       =  10.16678
RMS_2                            =   0.00489
Expected_RMS_2                   =   0.00104
Npoints_2                        =  3417
```

---

### `-BLSFixDurTc` — BLS with Fixed Transit Duration and Epoch

**Syntax**
```
-BLSFixDurTc
    <"duration" <"fix" dur | "var" varname | "expr" expression
        | "fixcolumn" <colname | colnum>
        | "list" ["column" col]>>
    <"Tc" <"fix" Tc | "var" varname | "expr" expression
        | "fixcolumn" <colname | colnum>
        | "list" ["column" col]>>
    ["fixdepth" <"fix" depth | "var" varname | "expr" expression
        | "fixcolumn" <colname | colnum>
        | "list" ["column" col]>
        ["qgress" <"fix" qgress | "var" varname | "expr" expression
            | "fixcolumn" <colname | colnum>
            | "list" ["column" col]>]]
    <"var" minpvar | "expr" minpexpr | minper>
    <"var" maxpvar | "expr" maxpexpr | maxper>
    <"var" nfvar | "expr" nfexpr | nfreq> timezone
    Npeak outperiodogram [outdir] omodel [model_outdir]
    correctlc ["fittrap"]
    ["ophcurve" outdir phmin phmax phstep]
    ["ojdcurve" outdir jdstep] ["maskpoints" maskvar]
```

**Description**

Run BLS with the transit duration and a reference epoch fixed. The period is still searched over a grid from `minper` to `maxper`. Optionally the transit depth and ingress fraction (`qgress`) can be fixed as well. For `qgress`: `0` = box-shaped transit; `0.5` = V-shaped (grazing) transit.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `"duration"` | Transit duration. Source: `"fix"` (command line), `"fixcolumn"`, or `"list"`. |
| `"Tc"` | Reference transit epoch. Same source options as `"duration"`. |
| `"fixdepth"` | Optionally fix the transit depth. |
| `"qgress"` | Optionally fix the ingress fraction. |
| `minper` / `maxper` / `nfreq` | Period range and number of trial frequencies. |
| All others | Same as `-BLS`. |

**References**

Cite Kovács, Zucker & Mazeh 2002, A&A, 391, 369.

---

### `-BLSFixPerDurTc` — BLS with Fixed Period, Duration, and Epoch

**Syntax**
```
-BLSFixPerDurTc
    <"period" <"fix" per | "var" varname | "expr" expression
        | "fixcolumn" <colname | colnum>
        | "list" ["column" col]>>
    <"duration" <"fix" dur | "var" varname | "expr" expression
        | "fixcolumn" <colname | colnum>
        | "list" ["column" col]>>
    <"Tc" <"fix" Tc | "var" varname | "expr" expression
        | "fixcolumn" <colname | colnum>
        | "list" ["column" col]>>
    ["fixdepth" <"fix" depth | "var" varname | "expr" expression
        | "fixcolumn" <colname | colnum>
        | "list" ["column" col]>
        ["qgress" <"fix" qgress | "var" varname | "expr" expression
            | "fixcolumn" <colname | colnum>
            | "list" ["column" col]>]]
    timezone omodel [model_outdir]
    correctlc ["fittrap"]
    ["ophcurve" outdir phmin phmax phstep]
    ["ojdcurve" outdir jdstep] ["maskpoints" maskvar]
```

**Description**

Run BLS with the period, transit duration, and reference epoch all fixed. Only the phase of the transit (within a single period) is searched. All options are similar to `-BLSFixDurTc`, except the period is also specified.

**References**

Cite Kovács, Zucker & Mazeh 2002, A&A, 391, 369.

---

### `-dftclean` — DFT Power Spectrum + CLEAN

**Syntax**
```
-dftclean <"var" nbvar | "expr" nbexpr | nbeam>
    ["maxfreq" <"var" mfvar | "expr" mfexpr | maxf>]
    ["outdspec" dspec_outdir]
    ["finddirtypeaks" Npeaks ["clip" <"var" cvar | "expr" cexpr | clip>
    clipiter]]
    ["outwfunc" wfunc_outdir]
    ["clean" <"var" gvar | "expr" gexpr | gain> <"var" snvar | "expr" snexpr |
    SNlimit>
        ["outcbeam" cbeam_outdir]
    ["outcspec" cspec_outdir]
    ["findcleanpeaks" Npeaks ["clip" <"var" cvar | "expr" cexpr | clip>
    clipiter]]]
    ["useampspec"] ["verboseout"] ["maskpoints" maskvar]
```

**Description**

Compute the Discrete Fourier Transform (DFT) power spectrum of the light curve using the FDFT algorithm (Kurtz 1985) and optionally deconvolve it with the CLEAN algorithm (Roberts, Lehar & Dreher 1987) to remove aliasing due to the window function.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `nbeam` | Number of frequency samples per `1/T` frequency element (`T` = light curve baseline). Controls the spectral resolution. |
| `"maxfreq" maxf` | Maximum frequency (cycles/day). Default: `1 / (2 · min_time_separation)` (Nyquist). |
| `"outdspec" dspec_outdir` | Write the dirty (uncleaned) power spectrum to `dspec_outdir/$basename.dftclean.dspec`. |
| `"finddirtypeaks" Npeaks` | Find the top `Npeaks` peaks in the dirty spectrum. |
| `"outwfunc" wfunc_outdir` | Write the window function to `wfunc_outdir/$basename.dftclean.wfunc`. |
| `"clean" gain SNlimit` | Apply the CLEAN algorithm. `gain` ∈ [0.1, 1.0] (smaller = slower convergence, more thorough); iterations continue until the highest remaining peak is below `SNlimit` × noise. |
| `"outcbeam" cbeam_outdir` | Write the clean beam to `cbeam_outdir/$basename.dftclean.cbeam`. |
| `"outcspec" cspec_outdir` | Write the clean power spectrum to `cspec_outdir/$basename.dftclean.cspec`. |
| `"findcleanpeaks" Npeaks` | Find the top `Npeaks` peaks in the clean spectrum. |
| `"clip" clip clipiter` | Clipping parameters for peak SNR calculation (default: iterative 5σ). |
| `"useampspec"` | Compute SNR on the amplitude spectrum instead of the power spectrum. |
| `"verboseout"` | Include the mean and standard deviation of the spectrum before and after clipping in the output. |
| `"maskpoints" maskvar` | Exclude points with `maskvar ≤ 0`. |

**References**

Cite Kurtz 1985, MNRAS, 213, 773 for the FDFT algorithm. Cite Roberts, Lehar & Dreher 1987, AJ, 93, 4 for the CLEAN algorithm.

**Examples**

**Example 1.** Compute the DFT power spectrum of a light curve up to 10 cycles/day, find the top dirty-spectrum peak using iterative 5σ clipping.

```bash
vartools -i EXAMPLES/2 -oneline -ascii \
    -dftclean 4 maxfreq 10. outdspec EXAMPLES/OUTDIR1 \
        finddirtypeaks 1 clip 5. 1
```

Output:
```
Name                         = EXAMPLES/2
DFTCLEAN_DSPEC_PEAK_FREQ_0_0 =     0.81189711
DFTCLEAN_DSPEC_PEAK_POW_0_0  = 0.000687634
DFTCLEAN_DSPEC_PEAK_SNR_0_0  =   59.8532
```

**Example 2.** Inject three harmonic signals into a light curve, compute the DFT, find the top 3 dirty-spectrum peaks, then apply the CLEAN algorithm and find the top 3 clean-spectrum peaks.

```bash
vartools -i EXAMPLES/4 -oneline -ascii \
    -Injectharm fix 0.697516 0 ampfix 0.1 phaserand 0 0 \
    -Injectharm fix 2.123456 0 ampfix 0.05 phaserand 0 0 \
    -Injectharm fix 0.426515 0 ampfix 0.01 phaserand 0 0 \
    -dftclean 4 maxfreq 10. outdspec EXAMPLES/OUTDIR1 \
        finddirtypeaks 3 clip 5. 1 \
        outwfunc EXAMPLES/OUTDIR1 \
        clean 0.5 5.0 outcbeam EXAMPLES/OUTDIR1 \
        outcspec EXAMPLES/OUTDIR1 \
        findcleanpeaks 3 clip 5. 1 \
        verboseout
```

---

### `-wwz` — Weighted Wavelet Z-Transform

**Syntax**
```
-wwz <"maxfreq" <"auto" | "var" v | "expr" e | maxfreq>>
    <"freqsamp" <"var" v | "expr" e | freqsamp>>
    <"tau0" <"auto" | "var" v | "expr" e | tau0>>
    <"tau1" <"auto" | "var" v | "expr" e | tau1>>
    <"dtau" <"auto" | "var" v | "expr" e | dtau>>
    ["c" <"var" v | "expr" e | cval>]
    ["outfulltransform" outdir ["fits" | "pm3d"] ["format" format]]
    ["outmaxtransform" outdir ["format" format]]
    ["maskpoints" maskvar]
```

**Description**

Compute the Weighted Wavelet Z-Transform (WWZ) as defined by Foster (1996) using an abbreviated Morlet wavelet:

```
f(z) = exp(i·2π·f·(t-τ) - c·(2π·f)²·(t-τ)²)
```

The transform is computed for all combinations of trial frequency (up to `maxfreq`) and time shift (`tau0` to `tau1` in steps of `dtau`). This yields a time-frequency map of the signal power that is especially useful for non-stationary signals.

The decay constant `c` (default: `1/(8π²)`) controls the trade-off between time and frequency resolution.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `"maxfreq" auto/val` | Maximum frequency in cycles/day. `"auto"` = `1 / (2 · min_time_separation)`. |
| `"freqsamp" freqsamp` | Frequency sampling as a multiple of `1/T`. |
| `"tau0" auto/val` | Start time for the time-shift scan. `"auto"` = minimum time in the light curve. |
| `"tau1" auto/val` | End time for the time-shift scan. `"auto"` = maximum time in the light curve. |
| `"dtau" auto/val` | Step size in time shift. `"auto"` = minimum time separation. |
| `"c" cval` | Morlet wavelet decay constant (default: `1/(8π²)`). |
| `"outfulltransform" outdir` | Write the full 2D transform (time-shift × frequency) to `outdir/$basename.wwz`. Add `"fits"` for multi-extension FITS, or `"pm3d"` for gnuplot pm3d-compatible text format. |
| `"outmaxtransform" outdir` | Write the transform maximized over frequencies as a function of time-shift to `outdir/$basename.mwwz`. |
| `"format" format` | Override the default filename convention (same syntax as `-o "nameformat"`). |
| `"maskpoints" maskvar` | Exclude points with `maskvar ≤ 0`. |

**Output columns**

| Column | Description |
|--------|-------------|
| `WWZ_maxZ_i` | Maximum value of the Z-transform over all time-shifts and frequencies. |
| `WWZ_maxfreq_i` | Frequency corresponding to the maximum Z. |
| `WWZ_maxpower_i` | Power at the maximum Z. |
| `WWZ_maxamp_i` | Amplitude at the maximum Z. |
| `WWZ_maxNeff_i` | Effective number of data points at the maximum Z. |
| `WWZ_maxtau_i` | Time shift at the maximum Z. |
| `WWZ_maxmeanmag_i` | Local mean magnitude at the maximum Z. |
| `WWZ_medZ_i` | Median over time-shifts of the maximum-Z frequency's Z value. |
| (and corresponding `med*` columns for the other quantities) | |

**References**

Cite Foster 1996, AJ, 112, 1709.

**Examples**

**Example 1.** Compute the WWZ for a light curve, scanning frequencies up to 2.0 cycles/day with 0.25/T step, time-shifts spanning the full observation range in 0.1-day steps, and output the full 2D transform in gnuplot pm3d format as well as the maximum-Z projection over time-shifts.

```bash
vartools -i EXAMPLES/8 -oneline \
    -wwz maxfreq 2.0 freqsamp 0.25 tau0 auto tau1 auto dtau 0.1 \
        outfulltransform EXAMPLES/OUTDIR1/ pm3d \
        outmaxtransform EXAMPLES/OUTDIR1
```

Output: Name, MaxWWZ values, frequencies, time shifts, power measurements, amplitude, effective N values, and median statistics. A significant periodic signal is found around time 53735.17392 with a frequency of 0.30645 cycles per day (absent at series endpoints).

---

### `-GetLSAmpThresh` — Minimum Detectable Amplitude

**Syntax**
```
-GetLSAmpThresh
    < "ls" | "list" ["column" col] > minp thresh
    < "harm" Nharm Nsubharm | "file" listfile >
    ["noGLS"]
```

**Description**

Determine the minimum peak-to-peak amplitude that a signal at a given period must have to be detected by a Lomb-Scargle search with `-ln(FAP) > thresh`. The signal shape is either a Fourier series or read from a file. The threshold is computed by scaling the signal template until the LS statistic reaches the detection limit.

**Parameters**

| Parameter | Description |
|-----------|-------------|
| `"ls"` | Use the period from the most recent `-LS` command. |
| `"list"` | Read the period from the input list (optional `"column"` keyword). |
| `minp` | Minimum period that would be searched (sets the FAP scale). |
| `thresh` | Desired `-ln(FAP)` detection threshold. |
| `"harm" Nharm Nsubharm` | Describe the signal shape as a Fourier series with `Nharm` harmonics and `Nsubharm` sub-harmonics. The series is fit to the light curve at the specified period. |
| `"file" listfile` | Two-column file: `signal_file signal_amp`, one row per light curve. Each `signal_file` provides the signal magnitude in column 3; `signal_amp` is the peak-to-peak amplitude of that signal. |
| `"noGLS"` | Compute the threshold for the traditional (non-generalized) LS periodogram. |

**Output columns**

| Column | Description |
|--------|-------------|
| `GetLSAmpThresh_minfactor_i` | Minimum factor by which the current signal must be scaled to remain detectable. |
| `GetLSAmpThresh_minamp_i` | Corresponding minimum peak-to-peak amplitude (magnitudes). |

**Examples**

**Example 1.** Run LS to find the best period, extract the sinusoidal signal with `-Killharm`, then compute the minimum detectable amplitude at that period given a threshold of `-ln(FAP) > 100`.

```bash
vartools -i EXAMPLES/2 -oneline \
    -LS 0.1 10. 0.1 1 0 \
    -Killharm ls 0 0 0 fitonly \
    -GetLSAmpThresh ls 0.1 -100 harm 0 0
```

Output:
```
Name                                 = EXAMPLES/2
LS_Period_1_0                        =     1.23440877
Log10_LS_Prob_1_0                    = -704.49194
LS_SNR_1_0                           =   58.45119
Killharm_Mean_Mag_1                  =  10.12217
Killharm_Period_1_1                  =     1.23440877
Killharm_Per1_Fundamental_Sincoeff_1 =   0.05008
Killharm_Per1_Fundamental_Coscoeff_1 =  -0.00222
Killharm_Per1_Amplitude_1            =   0.10026
LS_AmplitudeScaleFactor_2            =   0.02473
LS_MinimumAmplitude_2                =   0.00248
```

