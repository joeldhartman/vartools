# Model Fitting

Commands for fitting analytic models to light curves, from simple harmonic series to full transit and microlensing models.

---

## `-linfit`

```
-linfit function paramlist ["modelvar" varname]
    ["reject" <"var" scvar | "expr" scexpr | sigclip> ["useMAD"]
    ["iter" ["fixednum" number]]]
    ["correctlc"]
    ["omodel" model_outdir ["format" nameformat]]
    ["fitmask" maskvar]
```

Fit a function that is linear in its free parameters to each light curve using least squares.

**Parameters**

- `function` — Analytic function to fit (e.g., `'a1*t^2+a2*t+a3'`).
- `paramlist` — Comma-separated list of free parameter names (e.g., `'a1,a2,a3'`). Parameter names must not conflict with any existing vector variable names (`t`, `mag`, `err`, etc.).
- `"modelvar" varname` — Store the best-fit model in the vector variable `varname` for use by later commands.
- `"reject" sigclip` — Reject outliers beyond `sigclip` sigma and refit.
  - `"useMAD"` — Use the MAD statistic rather than the standard deviation when measuring scatter.
  - `"iter"` — Iteratively reject until no additional outliers are found.
  - `"fixednum" number` — Perform at most `number` rejection iterations.
- `"correctlc"` — Subtract the best-fit model from the light curve.
- `"omodel" model_outdir` — Output the model to `model_outdir`. Default filename: `model_outdir/basefilename.linfit.model`. Use `"format"` keyword to override the naming convention.
- `"fitmask" maskvar` — Only include points with `maskvar > 0` in the fit.

**Examples**

**Example 1.** Fit a quadratic polynomial to a light curve. The minimum time value is first computed via `-stats` and stored as a reference point, then linear regression is performed using a normalized time variable to avoid numerical precision issues.

```bash
vartools -i EXAMPLES/1 \
    -stats t min \
    -expr t0=STATS_t_MIN_0 \
    -linfit 'a*(t-t0)^2+b*(t-t0)+c' 'a,b,c' \
    -oneline
```

Output:
```
Name          = EXAMPLES/1
STATS_t_MIN_0 = 53725.173920000001
Linfit_a_2    = 0.00025540627746042932
Linfit_erra_2 = 1.9561241332987699e-07
Linfit_b_2    = 0.0097933162509034055
Linfit_errb_2 = 5.9117874714733109e-06
Linfit_c_2    = 10.083037598482507
Linfit_errc_2 = 3.2584974556662493e-05
```

---

## `-nonlinfit`

```
-nonlinfit function paramlist ["linfit" linfitparams]
    ["errors" error_expr]
    ["covariance"
        <"squareexp" amp_var rho_var
         | "exp" amp_var rho_var
         | "matern" amp_var rho_var nu_var>]
    ["priors" priorlist] ["constraints" constraintlist]
    <"amoeba" ["tolerance" <"var" v | "expr" e | tol>]
              ["maxsteps" <"var" v | "expr" e | steps>]
     | "mcmc" ["Naccept" <"var" v | "expr" e | N>
              | "Nlinkstotal" <"var" v | "expr" e | N>]
            ["fracburnin" <"var" v | "expr" e | frac>]
            ["eps" <"var" v | "expr" e | eps>] ["skipamoeba"]
            ["chainstats" exprlist statslist]
            ["maxmemstore" maxmem]
            ["outchains" outdir ["format" format] ["printevery" N]] >
    ["modelvar" varname] ["correctlc"]
    ["omodel" model_outdir ["format" nameformat]]
```

Fit a function that is nonlinear in its free parameters to each light curve.

**Parameters**

- `function` — Analytic function to fit (e.g., `'a*exp(-(t-t0)^2/2/sigma^2)+b'`).
- `paramlist` — Comma-separated list of free parameters with initial guesses and step sizes in the form `var=init:step` (e.g., `'t0=5.0:2.0,sigma=10.0:8.0'`).
- `"linfit" linfitparams` — List any parameters that enter linearly; these will be optimized using linear least squares, speeding up the fit. These parameters are excluded from MCMC posterior distributions.
- `"errors" error_expr` — Analytic expression for the magnitude uncertainties used in the likelihood function.
- `"covariance"` — Allow for correlated errors using one of three Gaussian process covariance models:
  - `"squareexp"` — Squared-exponential: covariance ∝ `amp_var * exp(-(t_i-t_j)²/(2*rho_var)²)`.
  - `"exp"` — Exponential: covariance ∝ `amp_var * exp(-|t_i-t_j|/rho_var)`.
  - `"matern"` — Matérn: parameterized by `amp_var`, `rho_var` (correlation length), and `nu_var` (shape). Linear fitting of a subset of parameters is not permitted when `"covariance"` is used.
- `"priors" priorlist` — Comma-separated list of prior expressions, each evaluating to `-2*ln(P)`. Example: `'(t0-4.0)^2/3.0^2'` for a Gaussian prior on `t0` with mean 4.0 and σ 3.0.
- `"constraints" constraintlist` — Comma-separated list of constraint expressions (e.g., `'sigma>0'`).
- `"amoeba"` — Use the downhill simplex (Nelder-Mead) optimizer.
  - `"tolerance" tol` — Convergence tolerance (minimum fractional change in χ² between iterations).
  - `"maxsteps" steps` — Maximum number of iterations.
- `"mcmc"` — Use Differential Evolution Markov Chain Monte Carlo.
  - `"Naccept" N` — Number of accepted links to run (default: `"Nlinkstotal"` with N=100000).
  - `"Nlinkstotal" N` — Total number of links.
  - `"fracburnin" frac` — Fraction of chain to discard as burn-in (default: 0.1).
  - `"eps" eps` — Scale of differential evolution random variations (default: 0.001).
  - `"skipamoeba"` — Skip the initial downhill simplex optimization.
  - `"chainstats" exprlist statslist` — Change which statistics are reported from the MCMC posterior. `exprlist` is a comma-separated list of analytic expressions; `statslist` is a comma-separated list of statistics (see [`-stats`](statistics.md#-stats)).
  - `"maxmemstore" maxmem` — Maximum memory for chain storage in GB (default: 4.0).
  - `"outchains" outdir` — Output MCMC chains to `outdir`. Each file is named `outdir/lcname.mcmc`. Use `"format"` to override naming. Use `"printevery" N` to thin the chain output.
- `"modelvar" varname` — Store the best-fit model in vector variable `varname`.
- `"correctlc"` — Subtract the best-fit model from the light curve.
- `"omodel" model_outdir` — Output the model to a file. Default suffix: `.nonlinfit.model`.

---

## `-Killharm`

```
-Killharm <"aov" | "ls" | "both" | "injectharm"
    | "fix" Nper <"var" v | "expr" e | per1> ... <"var" v | "expr" e | perN>
    | "list" Nper ["column" col1]> Nharm Nsubharm
    omodel [model_outdir] ["fitonly"]
    ["outampphase" | "outampradphase" | "outRphi" | "outRradphi"]
    ["clip" <"var" cvar | "expr" cexpr | val>]
```

Fit and optionally subtract a harmonic series from the light curve (i.e., whiten the light curve against one or more periods). The model has the form:

```
sum_{i=1}^{Nper}(
    sum_{k=0}^{Nharm_i}(a_i_k*sin(2*pi*(k+1)*f_i*t) + b_i_k*cos(2*pi*(k+1)*f_i*t))
  + sum_{k=0}^{Nsubharm_i}(c_i_k*sin(2*pi*f_i*t/(k+1)) + d_i_k*cos(2*pi*f_i*t/(k+1)))
)
```

The mean magnitude, the period(s), and all cosine and sine coefficients are output. By default the whitened light curve is passed to the next command; use `"fitonly"` to suppress subtraction.

**Period source**

| Keyword | Description |
|---------|-------------|
| `"aov"` | Period from the most recent `-aov` or `-aov_harm` command |
| `"ls"` | Period from the most recent `-LS` command |
| `"both"` | Two periods: one from `-aov`, one from `-LS` |
| `"injectharm"` | Period from the most recent `-Injectharm` command |
| `"fix" Nper per1...perN` | `Nper` fixed periods given on the command line |
| `"list" Nper ["column" col1]` | `Nper` periods read from the input list file |

**Parameters**

- `Nharm` — Number of higher harmonics to include (frequencies `2f₀`, `3f₀`, ... `(Nharm+1)f₀`).
- `Nsubharm` — Number of sub-harmonics (frequencies `f₀/2`, `f₀/3`, ... `f₀/(Nsubharm+1)`).
- `omodel` — Flag (`1` or `0`) to output the model light curve to `modeloutdir`. Output suffix: `.killharm.model`.
- `"fitonly"` — If given, the model is computed but not subtracted from the light curve.
- `"outampphase"` — Output amplitudes `A_k = sqrt(a_k² + b_k²)` and phases `φ_k = atan2(-b_k, a_k)/(2π)` instead of raw coefficients.
- `"outampradphase"` — Same as `"outampphase"` but phases in radians.
- `"outRphi"` — Output relative amplitudes `R_k1 = A_k/A_1` and phases `φ_k1 = φ_k - k*φ_1` (in units of 0 to 1).
- `"outRradphi"` — Same as `"outRphi"` but phases in radians.
- `"clip" val` — Fit the model, clip residuals at `val` sigma, then refit to the surviving points.

---

## `-MandelAgolTransit`

```
-MandelAgolTransit <bls | blsfixper
        | <"var" v | "expr" e | P0> <"var" v | "expr" e | T00>
          <"var" v | "expr" e | r0> <"var" v | "expr" e | a0>
          <"i" <"var" v | "expr" e | inclination> | "b" <"var" v | "expr" e |
          bimpact>>
          <"var" v | "expr" e | e0> <"var" v | "expr" e | omega0>
          <"var" v | "expr" e | mconst0>>
    <"quad" | "nonlin"> <"var" v | "expr" e | ldcoeff1_0>
        ... <"var" v | "expr" e | ldcoeffn_0> fitephem
    fitr fita fitinclterm fite fitomega fitmconst fitldcoeff1 ... fitldcoeffn
    fitRV [RVinputfile RVmodeloutfile <"var" v | "expr" e | K0>
    <"var" v | "expr" e | gamma0> fitK fitgamma]
    correctlc omodel [model_outdir] ["modelvar" var]
    ["ophcurve" curve_outdir phmin phmax phstep]
    ["ojdcurve" curve_outdir jdstep]
```

Fit a Mandel and Agol (2002) transit model to the light curve. Initial parameters can come from a preceding `-BLS` run or be entered directly.

**Initial parameters** (when not using `"bls"` or `"blsfixper"`)

| Parameter | Description |
|-----------|-------------|
| `P0` | Orbital period |
| `T00` | Time of transit center |
| `r0` | Ratio of planet radius to stellar radius (Rp/R★) |
| `a0` | Ratio of semi-major axis to stellar radius (a/R★) |
| `"i" incl` | Inclination (degrees), or `"b" bimp` for impact parameter |
| `e0` | Eccentricity |
| `omega0` | Argument of periastron (degrees) |
| `mconst0` | Out-of-transit magnitude (if negative, estimated from the data) |

**Limb darkening**

- `"quad"` — Quadratic limb darkening; requires 2 coefficients.
- `"nonlin"` — Non-linear (Claret) limb darkening; requires 4 coefficients.

**Fit flags**

`fitephem fitr fita fitinclterm fite fitomega fitmconst fitldcoeff1 ... fitldcoeffn` — Each is `1` to allow that parameter to vary, `0` to hold it fixed.

**Optional radial velocity fitting**

Set `fitRV=1` to simultaneously fit an RV curve from `RVinputfile` (columns: JD, RV, RVerror). The model RV (evaluated at evenly spaced phase points) is written to `RVmodeloutfile`. Provide initial values `K0` and `gamma0` and flags `fitK`, `fitgamma`.

**Output control**

- `correctlc` — Set to `1` to subtract the transit model from the light curve.
- `omodel` — Set to `1` to output the model; provide `model_outdir`. Output suffix: `.mandelagoltransit.model`.
- `"modelvar" var` — Store the best-fit model in vector variable `var`.
- `"ophcurve" curve_outdir phmin phmax phstep` — Output a model phase curve to `curve_outdir` with phases from `phmin` to `phmax` in steps of `phstep`. Suffix: `.mandelagoltransit.phcurve`.
- `"ojdcurve" curve_outdir jdstep` — Output a model light curve evaluated at times spanning the observations with step size `jdstep`. Suffix: `.mandelagoltransit.jdcurve`.

**Citation:** Mandel, K. & Agol, E. 2002, ApJ, 580, L171.

**Examples**

**Example 1.** Use `-BLS` to identify a transit signal in the light curve `EXAMPLES/3.transit`, then fit a Mandel-Agol transit model to it. A quadratic limb-darkening law is used with coefficients 0.3471 and 0.3180. The ephemeris (P and T0), Rp/R★, a/R★, and the impact parameter are varied; eccentricity and argument of periastron are held fixed.

```bash
vartools -i EXAMPLES/3.transit -oneline \
    -BLS q 0.01 0.1 0.5 5.0 20000 200 7 1 0 0 0 \
    -MandelAgolTransit bls quad 0.3471 0.3180 \
        1 1 1 1 0 0 1 0 0 0 0 1 EXAMPLES/OUTDIR1
```

Output:
```
BLS_Period_1_0              = 2.12312625
BLS_Tc_1_0                  = 53727.297046247397
MandelAgolTransit_Period_1  = 2.12328176
MandelAgolTransit_r_1       = 0.09789
MandelAgolTransit_bimpact_1 = 0.33094
MandelAgolTransit_chi2_1    = 27.06054
```

---

## `-SoftenedTransit`

```
-SoftenedTransit
    < bls | blsfixper | P0 T00 eta0 delta0 mconst0 >
    cval0 fitP fitT0 fiteta fitcval fitdelta fitmconst
    correctlc omodel [model_outdir]
    fit_harm [< "aov" | "ls" | "bls"
            | "list" ["column" col] | "fix" Pharm >
        nharm nsubharm]
```

Fit a Protopapas, Jimenez and Alcock (2005) "softened" transit model to the light curve. Initial parameters may come from a preceding `-BLS` or `-BLSFixPer` run, or be entered directly.

**Parameters**

- `P0`, `T00`, `eta0`, `delta0`, `mconst0` — Initial period, time of transit, transit duration parameter, depth, and out-of-transit magnitude (use a negative value for `mconst0` to estimate it from the data).
- `cval0` — Initial value for the softening parameter.
- `fitP fitT0 fiteta fitcval fitdelta fitmconst` — Flags (`1` to vary, `0` to hold fixed) for each parameter.
- `correctlc` — Set to `1` to subtract the model from the light curve.
- `omodel` — Set to `1` to output the model to `model_outdir`. Suffix: `.softenedtransit.model`.
- `fit_harm` — Set to `1` to simultaneously fit a harmonic series; specify the period source and `nharm`, `nsubharm`.

**Citation:** Protopapas, P., Jimenez, R. & Alcock, C. 2005, MNRAS, 362, 460.

**Examples**

**Example 1.** Use `-BLS` to identify a transit signal in `EXAMPLES/3.transit`, then fit a Protopapas et al. 2005 softened transit model initialized from the BLS results. The ephemeris, eta, cval, delta, and mconst are varied; the model is not subtracted from the light curve.

```bash
vartools -i EXAMPLES/3.transit -oneline \
    -BLS q 0.01 0.1 0.5 5.0 20000 200 7 1 0 0 0 \
    -SoftenedTransit bls 1 1 1 1 1 0 1 EXAMPLES/OUTDIR1 0
```

Output:
```
Name = EXAMPLES/3.transit
BLS_Period_1_0               = 2.12312625
BLS_Tc_1_0                   = 53727.297046247397
BLS_SN_1_0                   = 38.39425
BLS_SR_1_0                   = 0.00237
BLS_SDE_1_0                  = 4.77204
BLS_Depth_1_0                = 0.01136
BLS_Qtran_1_0                = 0.03000
BLS_deltaChi2_1_0            = -24130.93833
BLS_Ntransits_1_0            = 4
BLS_Rednoise_1_0             = 0.00156
BLS_Whitenoise_1_0           = 0.00490
BLS_SignaltoPinknoise_1_0    = 12.89679
SoftenedTransit_Period_1     = 2.12322112
SoftenedTransit_T0_1         = 53727.29783160
SoftenedTransit_eta_1        = 0.06171206
SoftenedTransit_cval_1       = -10.87159958
SoftenedTransit_delta_1      = -0.01206461
SoftenedTransit_mconst_1     = 10.16686817
SoftenedTransit_chi2perdof_1 = 27.04335183
```

---

## `-microlens`

```
-microlens
    ["f0"
        ["fix" fixval | "var" varname | "expr" expression
            | "list" ["column" col]
            | "fixcolumn" <colname | colnum>
            | "auto"]
        ["step" initialstepsize] ["novary"]]
    ["f1" ... ] ["u0" ... ] ["t0" ... ] ["tmax" ... ]
    ["correctlc"] ["omodel" outdir]
```

Fit a simple (Wozniak 2001) microlensing model to the light curve using a downhill simplex optimizer.

**Parameters**

For each of the five model parameters (`f0`, `f1`, `u0`, `t0`, `tmax`), optionally specify:

- `"fix" fixval` — Fix the initial value for all light curves.
- `"list" ["column" col]` — Read the initial value from the input list.
- `"fixcolumn" <colname | colnum>` — Use a previously computed statistic as the initial value.
- `"auto"` — Automatically determine the initial value (default for all parameters).
- `"step" initialstepsize` — Initial step size for the downhill simplex.
- `"novary"` — Hold this parameter fixed at its initial value.

**Output control**

- `"correctlc"` — Subtract the best-fit model from the light curve.
- `"omodel" outdir` — Output the model to `outdir`. Output suffix: `.microlens`.

**Citation:** Wozniak, P.R. et al. 2001, AcA, 51, 175.

**Examples**

**Example 1.** Fit a simple microlensing model to the simulated light curve `EXAMPLES/4.microlensinject`. Initial values for all five parameters are set automatically, and the best-fit model is written to `EXAMPLES/OUTDIR1/`.

```bash
vartools -i EXAMPLES/4.microlensinject -oneline \
    -microlens f0 auto f1 auto u0 auto t0 auto tmax auto \
        omodel EXAMPLES/OUTDIR1
```

Output:
```
Name                   = EXAMPLES/4.microlensinject
Microlens_f0_0         = 7.242316197338e-05
Microlens_f1_0         = 7.5541525219661e-05
Microlens_u0_0         = 7.242316197338e-05
Microlens_t0_0         = 3.9109521538222
Microlens_tmax_0       = 53740.494617109
Microlens_chi2perdof_0 = 4.4674961258953
```


---

## `-Starspot`

```
-Starspot
    <aov | ls | list ["column" col] | "fix" period
        | "var" varname | "expr" expression
        | "fixcolumn" <colname | colnum>>
    <"var" v | "expr" e | a0> <"var" v | "expr" e | b0>
    <"var" v | "expr" e | alpha0> <"var" v | "expr" e | i0>
    <"var" v | "expr" e | chi0> <"var" v | "expr" e | psi00>
    <"var" v | "expr" e | mconst0> fitP fita fitb
    fitalpha fiti fitchi fitpsi fitmconst correctlc omodel [model_outdir]
```

!!! warning "Deprecated"
    This command is deprecated as of VARTOOLS 1.3. Use the `-macula` extension command instead.

Fit a single, circular, uniform-temperature starspot model to the light curve using the Dorren (1987) model. Parameters `a0`, `b0`, `alpha0`, `i0`, `chi0`, `psi00` are as defined in Dorren 1987. Set `mconst0` negative to estimate it automatically from the data. Fit flags (`fitP`, `fita`, etc.) are `1` to vary and `0` to fix the corresponding parameter.

**Citation:** Dorren 1987, ApJ, 320, 756.

**Examples**

**Example 1.** Determine the rotation period via AOV analysis, then apply Dorren 1987 single-starspot modeling. Initial parameters: a=0.0298, b=0.08745, spot radius=20°, inclination=85°, latitude=30°, longitude=0°. The unspotted magnitude is estimated automatically. Period, spot radius, inclination, latitude, longitude, and magnitude are varied; a and b remain fixed. The best-fit model is output to `EXAMPLES/OUTDIR1/`.

```bash
vartools -i EXAMPLES/3.starspot -oneline \
    -aov Nbin 20 0.1 10. 0.1 0.01 5 0 \
    -Starspot aov 0.0298 0.08745 20. 85. 30. 0. -1 \
        1 0 0 1 1 1 1 1 0 1 EXAMPLES/OUTDIR1/
```

Output:
```
Name                   = EXAMPLES/3.starspot
Period_1_0             = 3.07960303
AOV_1_0                = 2861.35783
AOV_SNR_1_0            = 605.83431
AOV_NEG_LN_FAP_1_0    = 4755.85353
Starspot_Period_1      = 3.12218969
Starspot_a_1           = 0.02980
Starspot_b_1           = 0.08745
Starspot_alpha_1       = 22.51312
Starspot_inclination_1 = 69.03963
Starspot_chi_1         = 30.00411
Starspot_psi0_1        = 0.00000
Starspot_mconst_1      = 10.16641
Starspot_chi2perdof_1  = 26.58796
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

