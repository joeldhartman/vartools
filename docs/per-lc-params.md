# Per-LC variable/expression parameter catalogue

This file lists all numeric parameters across vartools commands that are candidates
for per-LC variable (`var varname`) or expression (`expr exprstring`) support.

Parameters are classified as:

| Tag | Meaning |
|-----|---------|
| **DONE** | Already upgraded — `_source`, `_var`, `_expr` companions exist in the struct |
| **CANDIDATE** | Safe to upgrade with `VT_PARAM_COMPANIONS` + `VT_PARSE_*` + `VT_EVAL_*` |
| **ALLOC_RISK** | Value feeds `malloc`/`calloc` sizing a per-command array — varying per-LC would require restructuring the allocation logic |
| **OUTPUT_FORMAT_RISK** | Value determines the number of output statistic columns — varying per-LC would break the fixed-width output table |
| **SKIP** | Internal/output/counter field — not appropriate for per-LC upgrade |

The macro header `src/vt_param_macros.h` provides `VT_PARAM_COMPANIONS`,
`VT_PARSE_DOUBLE/INT/FLOAT/LONG`, and `VT_EVAL_DOUBLE/INT/FLOAT/LONG`.

---

## Commands already fully or partially upgraded

### `-LS` — `_Ls`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `minp` | `double` | **DONE** | |
| `maxp` | `double` | **DONE** | |
| `subsample` | `double` | **DONE** | |

### `-aov` — `_Aov`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `minp` | `double` | **DONE** | |
| `maxp` | `double` | **DONE** | |
| `subsample` | `double` | **DONE** | |
| `finetune` | `double` | **DONE** | |
| `Nbin` | `int` | **DONE** | |

### `-aov_harm` — `_AovHarm`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `minp` | `double` | **DONE** | |
| `maxp` | `double` | **DONE** | |
| `subsample` | `double` | **DONE** | |
| `finetune` | `double` | **DONE** | |
| `Nharm` | `int` | **DONE** | Note: Nharm also controls `peakNharm` output array — verify no allocation issue |

### `-BLS` — `_Bls`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `minper` | `double` | **DONE** | |
| `maxper` | `double` | **DONE** | |
| `nbins` | `int` | **DONE** | |
| `rmin` | `double` | **DONE** | |
| `rmax` | `double` | **DONE** | |
| `qmin` | `double` | **DONE** | |
| `qmax` | `double` | **DONE** | |
| `rho` | `double` | **DONE** | |
| `df` | `double` | **DONE** | |
| `nf` | `int` | **DONE** | |
| `subsample` | `double` | **DONE** | |
| `minexpdurfrac` | `double` | **DONE** | |
| `maxexpdurfrac` | `double` | **DONE** | |

### `-harmonicfilter` — `_HarmonicFilter`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `maxfreq` | `double` | **DONE** | |
| `minfreq` | `double` | **DONE** | |

### `-resample` — `_Resample`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `tstart` | `double` | **DONE** | |
| `tstop` | `double` | **DONE** | |
| `Nresamp` | `int` | **DONE** | |
| `delt` | `double` | **DONE** | |
| `minsep` | `double` | **DONE** | |

---

## Candidate commands — proposed upgrades

### `-clip` — `_Clip`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `sigclip` | `double` | **CANDIDATE** | Sigma-clipping threshold |
| `iter` | `int` | **CANDIDATE** | Current iteration counter — check if this is an *input* iteration limit or internal counter |
| `niter` | `int` | **CANDIDATE** | Max iterations — safe input parameter |

### `-ensemblerescalesig` — `_Ensemblerescalesig`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `erssigclip` | `double` | **CANDIDATE** | Sigma-clipping level used internally |
| `a` | `double` | **SKIP** | Output fit coefficient |
| `b` | `double` | **SKIP** | Output fit coefficient |

### `-killharm` — `_Killharm`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `Nharm` | `int` | **ALLOC_RISK** | Sizes `harmA[Nharm]`, `harmB[Nharm]`, `amp[Nharm]` arrays and drives OUTPUT_FORMAT_RISK (harmonic amplitude/phase columns per peak) |
| `Nsubharm` | `int` | **ALLOC_RISK** | Sizes `subharmA[Nsubharm]`, `subharmB[Nsubharm]` |
| `clip` | `double` | **CANDIDATE** | Sigma-clipping threshold for harmonic fit |

### `-injectharm` — `_Injectharm`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `Nharm` | `int` | **ALLOC_RISK** | Sizes `harm_amp*`, `harm_phase*`, `harm_ampspec`, etc. arrays |
| `Nsubharm` | `int` | **ALLOC_RISK** | Sizes `subharm_*` arrays |
| `fixperiod` | `double` | **CANDIDATE** | Fixed period to inject (when `pertype == PERTYPE_FIX`) |
| `minp` | `double` | **CANDIDATE** | Min period for random injection |
| `maxp` | `double` | **CANDIDATE** | Max period for random injection |
| `minf` | `double` | **CANDIDATE** | Min frequency for random injection |
| `maxf` | `double` | **CANDIDATE** | Max frequency for random injection |

### `-injecttransit` — `_Injecttransit`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `Nld` | `int` | **ALLOC_RISK** | Sizes limb-darkening coefficient arrays |
| `minp` | `double` | **CANDIDATE** | Min period |
| `maxp` | `double` | **CANDIDATE** | Max period |
| `minf` | `double` | **CANDIDATE** | Min frequency |
| `maxf` | `double` | **CANDIDATE** | Max frequency |
| `minRp` | `double` | **CANDIDATE** | Min planet radius |
| `maxRp` | `double` | **CANDIDATE** | Max planet radius |
| `minMp` | `double` | **CANDIDATE** | Min planet mass |
| `maxMp` | `double` | **CANDIDATE** | Max planet mass |
| `paramfix[14]` | `double[14]` | **CANDIDATE** | Fixed initial values for each injection parameter; per-element upgrade possible but awkward — NEEDS_REVIEW for best approach |

### `-mandelagoltransit` — `_MandelAgolTransit`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `P0` | `double` | **CANDIDATE** | Initial/fixed period |
| `T00` | `double` | **CANDIDATE** | Initial/fixed transit epoch |
| `r0` | `double` | **CANDIDATE** | Initial/fixed planet-to-star radius ratio |
| `a0` | `double` | **CANDIDATE** | Initial/fixed scaled semi-major axis |
| `inc0` | `double` | **CANDIDATE** | Initial/fixed inclination |
| `bimpact0` | `double` | **CANDIDATE** | Initial/fixed impact parameter |
| `sin_i0` | `double` | **CANDIDATE** | Initial/fixed sin(i) — only one of inc0/bimpact0/sin_i0 is active depending on `inputinclterm` |
| `e0` | `double` | **CANDIDATE** | Initial/fixed eccentricity |
| `omega0` | `double` | **CANDIDATE** | Initial/fixed argument of periapsis |
| `mconst0` | `double` | **CANDIDATE** | Initial/fixed magnitude constant (out-of-transit level) |
| `K0` | `double` | **CANDIDATE** | Initial/fixed RV semi-amplitude |
| `gamma0` | `double` | **CANDIDATE** | Initial/fixed RV systematic velocity |
| `ldcoeffs0[4]` | `double[4]` | **CANDIDATE** | Limb-darkening coefficients; array is fixed size 4 — could expose as `ldcoeffs0_0` … `ldcoeffs0_3` or as a single vector parameter |
| `nldcoeff` | `int` | **ALLOC_RISK** | Sizes `ldcoeffs[nldcoeff]` output array and controls number of LD fit columns |

### `-starspot` — `_Starspot`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `a0` | `double` | **CANDIDATE** | Initial spot amplitude |
| `b0` | `double` | **CANDIDATE** | Initial spot size |
| `chi0` | `double` | **CANDIDATE** | Initial spot latitude |
| `inclination0` | `double` | **CANDIDATE** | Initial stellar inclination |
| `alpha0` | `double` | **CANDIDATE** | Initial spot longitude rate |
| `psi00` | `double` | **CANDIDATE** | Initial spot longitude |
| `mconst0` | `double` | **CANDIDATE** | Initial magnitude constant |
| `fixedperiod` | `double` | **CANDIDATE** | Fixed rotation period (when `pertype == PERTYPE_FIX`) |

### `-microlens` — `_MicroLens`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `f00_fix` | `double` | **CANDIDATE** | Initial baseline flux |
| `f10_fix` | `double` | **CANDIDATE** | Initial blend flux |
| `u00_fix` | `double` | **CANDIDATE** | Initial minimum impact parameter |
| `t00_fix` | `double` | **CANDIDATE** | Initial time of closest approach |
| `tmax0_fix` | `double` | **CANDIDATE** | Initial Einstein crossing time |

*Note: `_MicroLens` already uses `f0_source` / `f0_linkedcolumn` — check whether those use the var/expr mechanism or a different linked-column mechanism before implementing.*

### `-phase` — `_Phase`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `fixperiod` | `double` | **CANDIDATE** | Period used when `pertype == PERTYPE_FIX` |
| `fixT0` | `double` | **CANDIDATE** | Reference epoch |
| `phaseTc` | `double` | **CANDIDATE** | Phase of transit centre |
| `startphase` | `double` | **CANDIDATE** | Starting phase value for output |

### `-binlc` — `_Binlc`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `binsize` | `double` | **CANDIDATE** | Bin width (when `binsize_Nbins_flag == 0`) |
| `firstbin` | `double` | **CANDIDATE** | Starting edge of first bin |
| `Nbins` | `int` | **OUTPUT_FORMAT_RISK** | Number of bins (when `binsize_Nbins_flag == 1`) — determines output column count |
| `t0fixval` | `double` | **CANDIDATE** | Fixed reference epoch for phase-bin mode |

### `-medianfilter` — `_MedianFilter`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `time` | `double` | **CANDIDATE** | Filter half-width timescale |

### `-addnoise` — `_AddNoise`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `gammaval_fix` | `double` | **CANDIDATE** | White-noise amplitude (noise_type == WHITE) |
| `sig_r_fix` | `double` | **CANDIDATE** | Red-noise amplitude |
| `sig_w_fix` | `double` | **CANDIDATE** | White-noise component in correlated model |
| `rho_r_fix` | `double` | **CANDIDATE** | Correlation length for exponential/Matern |
| `nu_r_fix` | `double` | **CANDIDATE** | Matern smoothness parameter |
| `bintime_fix` | `double` | **CANDIDATE** | Binning timescale for wavelet model |

*Note: `_AddNoise` uses `gammaval_type`, `sig_r_type`, etc. as source flags but with a different mechanism (not the var/expr macro triple). Review whether to migrate to the standard triple or keep the existing scheme.*

### `-autocorr` — `_Autocorr`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `start` | `double` | **CANDIDATE** | Start lag |
| `stop` | `double` | **CANDIDATE** | Stop lag |
| `step` | `double` | **CANDIDATE** | Lag step — also drives output column count (NEEDS_REVIEW for OUTPUT_FORMAT_RISK) |
| `errsize` | `double` | **CANDIDATE** | Error bar size |

### `-jstet` — `_Jstet`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `Jstet_time` | `double` | **CANDIDATE** | Characteristic timescale |
| `wkmax` | `double` | **CANDIDATE** | Maximum kernel frequency |

### `-dftclean` — `_Dftclean`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `gain` | `double` | **CANDIDATE** | CLEAN loop gain |
| `SNlimit` | `double` | **CANDIDATE** | S/N stopping threshold |
| `maxfreq` | `double` | **CANDIDATE** | Maximum frequency searched |
| `clip_dirty` | `double` | **CANDIDATE** | Sigma-clipping for dirty spectrum peaks |
| `clip_clean` | `double` | **CANDIDATE** | Sigma-clipping for clean spectrum peaks |
| `nbeam` | `int` | **SKIP** | Number of beam harmonics — internal algorithm parameter |
| `Npeaks_dirty` | `int` | **OUTPUT_FORMAT_RISK** | Number of output dirty-spectrum peaks |
| `Npeaks_clean` | `int` | **OUTPUT_FORMAT_RISK** | Number of output clean-spectrum peaks |

### `-wwz` — `_WWZ`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `cterm` | `double` | **CANDIDATE** | WWZ decay constant |
| `maxfreq` | `double` | **CANDIDATE** | Maximum frequency |
| `freq_sample_factor` | `double` | **CANDIDATE** | Frequency oversampling factor |
| `tau0` | `double` | **CANDIDATE** | Start time (when not auto) |
| `tau1` | `double` | **CANDIDATE** | End time (when not auto) |
| `dtau` | `double` | **CANDIDATE** | Time step (when not auto) |

### `-nonlinfit` — `_Nonlinfit`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `amoeba_tol` | `double` | **CANDIDATE** | Amoeba convergence tolerance |
| `amoeba_maxsteps` | `long` | **CANDIDATE** | Max amoeba iterations |
| `mcmc_Naccept` | `long` | **CANDIDATE** | Target accepted MCMC steps |
| `mcmc_Nlinkstotal` | `long` | **CANDIDATE** | Total MCMC chain links |
| `mcmc_burninfrac` | `double` | **CANDIDATE** | Burn-in fraction |
| `mcmc_eps` | `double` | **CANDIDATE** | MCMC proposal step scale |
| `rejsigclip` | `double` | **CANDIDATE** | Outlier rejection sigma (when `rejectoutliers`) |

### `-BLSFixDurTc` — `_BlsFixDurTc`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `minper` | `double` | **CANDIDATE** | Minimum period searched |
| `maxper` | `double` | **CANDIDATE** | Maximum period searched |
| `fixdur` | `double` | **CANDIDATE** | Fixed transit duration (when `durtype` is fixed) |
| `fixTC` | `double` | **CANDIDATE** | Fixed transit centre (when `TCtype` is fixed) |
| `fixdepthval` | `double` | **CANDIDATE** | Fixed transit depth (when `fixdepth`) |
| `qgressval` | `double` | **CANDIDATE** | Fixed ingress/egress fraction (when `qgresstype` is fixed) |

### `-BLSFixPerDurTc` — `_BlsFixPerDurTc`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `fixper` | `double` | **CANDIDATE** | Fixed period |
| `fixdur` | `double` | **CANDIDATE** | Fixed transit duration |
| `fixTC` | `double` | **CANDIDATE** | Fixed transit centre |
| `fixdepthval` | `double` | **CANDIDATE** | Fixed transit depth |
| `qgressval` | `double` | **CANDIDATE** | Fixed ingress/egress fraction |

### `-difffluxtomag` — `_DiffFluxtomag`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `mag_constant1` | `double` | **CANDIDATE** | Reference magnitude constant |
| `offset` | `double` | **CANDIDATE** | Flux offset |

### `-fluxtomag` — `_Fluxtomag`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `mag_constant1` | `double` | **CANDIDATE** | Reference magnitude constant |
| `offset` | `double` | **CANDIDATE** | Flux offset |

### `-findblends` — `_FindBlends`

| Field | C type | Status | Notes |
|-------|--------|--------|-------|
| `matchrad` | `double` | **CANDIDATE** | Matching radius (arcsec) |
| `zeromag` | `double` | **CANDIDATE** | Reference zero-point magnitude |
| `fixperiod` | `double` | **CANDIDATE** | Fixed period for harmonic model (when `pertype == PERTYPE_FIX`) |
| `Nharm` | `int` | **ALLOC_RISK** | Sizes harmonic model arrays |

---

## Commands with no candidate scalar input parameters

| Command | Reason |
|---------|--------|
| `-rms` / `_RMS_NoBin` | No user-configurable scalar input parameters; `Nbin` in `_RMS_Bin` is OUTPUT_FORMAT_RISK |
| `-chi2` / `_Chi2_NoBin` | Same; `Nbin` in `_Chi2_Bin` is OUTPUT_FORMAT_RISK |
| `-rescalesig` | No scalar input parameters |
| `-decorr` | Externally-driven by column lists, not scalar thresholds |
| `-linfit` | Parameters defined by expression strings, not numeric scalars in struct |
| `-sysrem` | Internal ensemble method, parameters from file lists |
| `-tfa` | Parameters from template file, not struct scalars |
| `-converttime` | RA/Dec already handled via `_Expression`; time system by enum |
| `-savelc` / `-restorelc` | Pipeline state, no numeric parameters |
| `-if` / `-fi` | Expression-based, no numeric scalar parameters |
| `-restricttimes` | Already uses linked-column / expression mechanism for JD bounds |
| `-expression` | Pure expression command |
| `-print` | Output-only |
| `-o` | Output-only |
| `-sortlc` | Sort by variable, no numeric threshold |
| `-stats` | Statistics selection, no scalar thresholds |
| `-fft` | No scalar numeric parameters in struct |
| `-match` | String/column-based matching |

---

## Implementation notes

1. **Include the macro header** — add `#include "vt_param_macros.h"` near the top of
   `commands.h` (after `analytic.h` which defines `_Variable` and `_Expression`).

2. **Per-command checklist** for each CANDIDATE parameter:
   a. Add `VT_PARAM_COMPANIONS(field)` after the existing field in `commands.h`.
   b. Replace the parse-time `atof`/`atoi` call in `parsecommandline.c` with
      `VT_PARSE_DOUBLE(c[cn].Cmd, field, argv, i)` (or `_INT`, `_FLOAT`, `_LONG`).
   c. Replace the runtime use of `cmd->field` in the execution `.c` file with
      `VT_EVAL_DOUBLE(cmd, field, lcindex, threadindex)`.

3. **`ldcoeffs0[4]`** — fixed-size array in `_MandelAgolTransit`; the cleanest
   approach is to expose the individual elements as `ldcoeffs0_0` … `ldcoeffs0_3`
   using separate companion triples, or to skip per-LC support for those and keep
   only the scalar form.  Mark as NEEDS_REVIEW until the user decides.

4. **`_AddNoise`** — already has `_type` / `_fix` / `**vals` fields, but using a
   non-standard mechanism.  Migrating to `VT_PARAM_COMPANIONS` would require
   removing/replacing the existing `_type` fields to avoid naming conflicts.
   Treat as NEEDS_REVIEW.

5. **`_MicroLens`** — already has `f0_source`, `f0_linkedcolumn`, etc. via a
   linked-column mechanism.  Check whether those cover the var/expr use case before
   adding companions; if so, those parameters are already effectively DONE via a
   different path.
