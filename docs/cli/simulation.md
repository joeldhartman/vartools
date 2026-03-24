# Simulation

Commands for injecting signals and noise into light curves, and for replicating light curves to support Monte Carlo experiments.

---

## `-addnoise`

```
-addnoise
    <   "white"
                <"sig_white" <"fix" val | "list" ["column" col]>>
      | "squareexp"
                <"rho" <"fix" val | "list" ["column" col]>>
                <"sig_red" <"fix" val | "list" ["column" col]>>
                <"sig_white" <"fix" val | "list" ["column" col]>>
                ["bintime" <"fix" val | "list" ["column" col]>]
      | "exp"
                <"rho" <"fix" val | "list" ["column" col]>>
                <"sig_red" <"fix" val | "list" ["column" col]>>
                <"sig_white" <"fix" val | "list" ["column" col]>>
                ["bintime" <"fix" val | "list" ["column" col]>]
      | "matern"
                <"nu" <"fix" val | "list" ["column" col]>>
                <"rho" <"fix" val | "list" ["column" col]>>
                <"sig_red" <"fix" val | "list" ["column" col]>>
                <"sig_white" <"fix" val | "list" ["column" col]>>
                ["bintime" <"fix" val | "list" ["column" col]>]
      | "wavelet"
                <"gamma" <"fix" val | "list" ["column" col]>>
                <"sig_red" <"fix" val | "list" ["column" col]>>
                <"sig_white" <"fix" val | "list" ["column" col]>>
    >
```

Add time-correlated Gaussian noise to the light curve. The user must choose one of five covariance models.

For every numerical parameter, supply the value using one of:

- `"fix" val` — Fixed value for all light curves.
- `"list" ["column" col]` — Read the value from the input light curve list. By default the next available column is used; use `"column" col` to specify explicitly.

### Noise models

#### `"white"` — Pure white (uncorrelated) noise

Adds independent Gaussian noise with standard deviation `sig_white` to each point.

| Parameter | Description |
|-----------|-------------|
| `sig_white` | Standard deviation of the white noise |

#### `"squareexp"` — Squared-exponential Gaussian process

Covariance between times *t_i* and *t_j*:

```
C(t_i, t_j) = sig_white² * δ_ij + sig_red² * exp(-(t_i-t_j)² / (2*rho²))
```

Both `rho` and `sig_red` must be greater than zero.

| Parameter | Description |
|-----------|-------------|
| `rho` | Correlation timescale |
| `sig_red` | Amplitude of the correlated (red noise) component |
| `sig_white` | Amplitude of the uncorrelated (white noise) component |
| `"bintime"` | Optional. Chunk the light curve into bins of this duration (same units as time) before simulating correlated noise in each bin independently. Substantially speeds up simulations when the light curve duration is much longer than `rho`. |

#### `"exp"` — Exponentially decaying Gaussian process

Covariance:

```
C(t_i, t_j) = sig_white² * δ_ij + sig_red² * exp(-|t_i-t_j| / rho)
```

Parameters and `"bintime"` option are the same as for `"squareexp"`.

#### `"matern"` — Matérn Gaussian process

Covariance:

```
C(t_i, t_j) = sig_white² * δ_ij + sig_red² * C(nu, x) * K_nu(x)
```

where `x = sqrt(2*nu) * |t_i-t_j| / rho`, `C(x,y) = (2^(1-x)/Gamma(x)) * y^x`, and `K_nu` is the modified Bessel function of the second kind. When `nu → ∞` the Matérn covariance converges to the squared-exponential; when `nu = 0.5` it equals the exponential covariance.

| Parameter | Description |
|-----------|-------------|
| `nu` | Shape parameter (must be > 0) |
| `rho` | Correlation timescale (must be > 0) |
| `sig_red` | Amplitude of the correlated component (must be > 0) |
| `sig_white` | Amplitude of the uncorrelated component |
| `"bintime"` | Optional binning acceleration (see `"squareexp"`) |

#### `"wavelet"` — 1/f^γ red noise + white noise

Generates noise as the sum of a red-noise component with power-spectral density proportional to `1/f^gamma` (γ must satisfy `-1 < gamma < 1`) with standard deviation `sig_red`, and an uncorrelated white-noise component with standard deviation `sig_white`. The red-noise is generated using the wavelet method of McCoy and Walden (1996).

| Parameter | Description |
|-----------|-------------|
| `gamma` | Power-law index for the red noise PSD; must satisfy `-1 < gamma < 1` |
| `sig_red` | Standard deviation of the red noise component |
| `sig_white` | Standard deviation of the white noise component |

**Examples**

**Example 1.** Simulate a light curve with time-correlated noise using the wavelet model. The red-noise component has power spectral density proportional to 1/f^0.99 and standard deviation 0.005; the white-noise component also has standard deviation 0.005.

```bash
gawk '{print $1, 0., 0.005}' EXAMPLES/1 | \
  vartools -i - -header -randseed 1 \
  -addnoise wavelet gamma fix 0.99 sig_red fix 0.005 sig_white fix 0.005 \
  -o EXAMPLES/OUTDIR1/noisesim.txt
```

**Example 2.** Same as above, using a squared-exponential model for the red-noise component with a correlation timescale of 0.01 days and standard deviation 0.005 mag. An additional white-noise component is included with standard deviation 0.001.

```bash
gawk '{print $1, 0., 0.005}' EXAMPLES/1 | \
  vartools -i - -header -randseed 1 \
  -addnoise squareexp rho fix 0.01 sig_red fix 0.005 sig_white fix 0.001 \
  -o EXAMPLES/OUTDIR1/noisesim.txt
```

---

## `-Injectharm`

```
-Injectharm
    < "list" ["column" col] | "fix" per
        | "rand" minp maxp | "logrand" minp maxp
        | "randfreq" minf maxf | "lograndfreq" minf maxf >
    Nharm
        (< "amplist" ["column" col] | "ampfix" amp
            | "amprand" minamp maxamp | "amplogrand" minamp maxamp >
         ["amprel"]
         < "phaselist" ["column" col] | "phasefix" phase | "phaserand" >
         ["phaserel"]) 0...Nharm
    Nsubharm
        (< "amplist" ["column" col] | "ampfix" amp
            | "amprand" minamp maxamp | "amplogrand" minamp maxamp >
         ["amprel"]
         < "phaselist" ["column" col] | "phasefix" phase | "phaserand" >
         ["phaserel"]) 1...Nsubharm
    omodel [modeloutdir]
```

Add a harmonic (Fourier) series signal to the light curve. The injected signal has the form:

```
A_1*cos(2*π*(t/P + φ_1))
    + sum_{k=2}^{Nharm+1} A_k*cos(2*π*(t*k/P + φ_k))
    + sum_{k=2}^{Nsubharm+1} A_k*cos(2*π*(t/k/P + φ_k))
```

**Period source**

| Keyword | Description |
|---------|-------------|
| `"list" ["column" col]` | Read from the input list |
| `"fix" per` | Fixed period for all light curves |
| `"rand" minp maxp` | Uniform random period in `[minp, maxp]` |
| `"logrand" minp maxp` | Uniform random period in log space |
| `"randfreq" minf maxf` | Uniform random frequency |
| `"lograndfreq" minf maxf` | Uniform random frequency in log space |

**Harmonic specification**

For each of the `Nharm+1` harmonics (fundamental = harmonic 1) and `Nsubharm` sub-harmonics, specify the amplitude and phase:

**Amplitude keywords**

| Keyword | Description |
|---------|-------------|
| `"amplist" ["column" col]` | Read from input list |
| `"ampfix" amp` | Fixed amplitude |
| `"amprand" minamp maxamp` | Uniform random amplitude |
| `"amplogrand" minamp maxamp` | Uniform log-random amplitude |
| `"amprel"` | Treat the specified amplitude as a ratio relative to the fundamental amplitude `A_k/A_1` |

**Phase keywords**

| Keyword | Description |
|---------|-------------|
| `"phaselist" ["column" col]` | Read from input list |
| `"phasefix" phase` | Fixed phase at `t=0` |
| `"phaserand"` | Uniform random phase in `[0, 1)` |
| `"phaserel"` | Treat the phase as relative to the fundamental: `φ_k1 = φ_k - k*φ_1` |

**Output**

- `omodel` — `1` to write the model light curve to `modeloutdir`. Output suffix: `.injectharm.model`.

**Examples**

**Example 1.** Inject a sinusoid into a light curve with a random period between 1.0–5.0 days, zero harmonic overtones, a uniform-log amplitude distribution, and a random phase.

**Example 2.** Inject an RR Lyrae signal into multiple light curves with a fixed 0.514333-day period and 10 harmonics. The detection success threshold is at fundamental amplitudes of 0.00195 or higher.

---

## `-Injecttransit`

```
-Injecttransit
    < "Plist" ["column" col] | "Pfix" per
        | "Pexpr" expr | "Prand" minp maxp | "Plogrand" minp maxp
        | "randfreq" minf maxf | "lograndfreq" minf maxf >
    < "Rplist" ["column" col] | "Rpfix" Rp |
        "Rpexpr" expr | "Rprand" minRp maxRp | "Rplogrand" minRp maxRp >
    < "Mplist" ["column" col] | "Mpfix" Mp |
        "Mpexpr" expr | "Mprand" minMp maxMp | "Mplogrand" minMp maxMp >
    < "phaselist" ["column" col] | "phasefix" phase |
        "phaseexpr" expr | "phaserand" >
    < "sinilist" ["column" col] | "sinifix" sin_i |
        "siniexpr" expr | "sinirand" >
    < "eomega"
            < "elist" ["column" col] | "efix" e |
                "eexpr" expr | "erand" >
            < "olist" ["column" col] | "ofix" omega |
                "oexpr" expr | "orand" >
        | "hk"
            < "hlist" ["column" col] | "hfix" h |
                "hexpr" expr | "hrand" >
            < "klist" ["column" col] | "kfix" k |
                "kexpr" expr | "krand" > >
    < "Mstarlist" ["column" col] | "Mstarfix" Mstar | "Mstarexpr" expr >
    < "Rstarlist" ["column" col] | "Rstarfix" Rstar | "Rstarexpr" expr >
    < "quad" | "nonlin" >
        < "ldlist" ["column" col] | "ldfix" ld1 ... ldn |
            "ldexpr" ld1 ... ldn >
    ["dilute" < "list" ["column" col] | "fix" dilute | "expr" diluteexpr >]
    omodel [modeloutdir]
```

Add a Mandel-Agol limb-darkened transit signal to the light curve.

**Parameters**

For each physical parameter, the source of the value must be specified. The keyword prefix determines the source:

| Prefix suffix | Description |
|---------------|-------------|
| `*list ["column" col]` | Read from the input light curve list |
| `*fix value` | Fixed value for all light curves |
| `*expr expr` | Analytic expression evaluated per light curve |
| `*rand min max` | Uniform random value in `[min, max]` |
| `*logrand min max` | Uniform log-random value |

**Physical parameters**

| Parameter | Units | Description |
|-----------|-------|-------------|
| `P` | days | Orbital period (or frequency in 1/day via `randfreq`/`lograndfreq`) |
| `Rp` | Jupiter radii | Planet radius |
| `Mp` | Jupiter masses | Planet mass (used to compute semi-major axis) |
| `phase` | — | Phase of transit center at `T=0`; `phase=0` corresponds to mid-transit |
| `sini` | — | Sine of the orbital inclination. Use `"sinirand"` to draw from a uniform orientation distribution constrained to produce a transit |
| `e`, `omega` | —, degrees | Eccentricity and argument of periastron (use `"eomega"` keyword group) |
| `h`, `k` | — | Alternatively specify `h = e*sin(omega)` and `k = e*cos(omega)` (use `"hk"` keyword group) |
| `Mstar` | solar masses | Stellar mass |
| `Rstar` | solar radii | Stellar radius |

**Limb darkening**

- `"quad"` — Quadratic limb darkening; 2 coefficients required.
- `"nonlin"` — Non-linear (Claret) limb darkening; 4 coefficients required.

**Dilution**

- `"dilute"` — Optional dilution factor (flux fraction from the target). Reduces the transit depth by this factor.

**Output**

- `omodel` — `1` to write the model light curve to `modeloutdir`. Suffix: `.injecttransit.model`.

**Examples**

**Example 1.** Inject a transit into a light curve and recover it using BLS analysis. The period is drawn from a uniform-log random distribution between 0.2 and 2.0 cycles per day. Planet radius and mass are fixed at 1.0 Jupiter values, phase and inclination are randomized, eccentricity and argument of periastron are set to zero, and stellar parameters are fixed at solar values with quadratic limb darkening coefficients 0.3471 and 0.3180. Output parameters include the injected transit properties and BLS detection results.

---

## `-copylc`

```
-copylc
    Ncopies
```

Replicate the current light curve `Ncopies` times in memory. Each copy is processed independently by all subsequent VARTOOLS commands. Data from commands preceding `-copylc` is replicated in the output table for each copy.

**Parameters**

- `Ncopies` — Number of copies to create.

Each copy has the suffix `_copy$copycommandnum.$copynum` appended to its name, where `$copycommandnum` is the index of the `-copylc` command that created it (useful when multiple `-copylc` commands are used) and `$copynum` runs from `0` to `Ncopies - 1`.

!!! note
    `-copylc` cannot be used together with the `-readall` option.

**Typical use**

Combine with `-addnoise` and a period-finding or signal-detection command to perform Monte Carlo false-alarm probability estimates:

```bash
vartools -l lclist.txt \
    -copylc 1000 \
    -addnoise white "sig_white" fix 0.005 \
    -LS 0.5 20.0 4.0 1 0 \
    -print LS_Period_1_1,LS_lnFAP_1_1
```
