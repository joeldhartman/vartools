# pyvartools Python API — Completeness Audit

**Purpose**: Systematic comparison of every Python command wrapper against the
actual CLI help text to find missing options, wrong keywords, and design gaps.

**How to read this document**

Each command that has gaps gets a section with a gap table.  Commands confirmed
complete are listed in the summary table only.  Use the severity codes below to
prioritise fix work.

| Code | Meaning |
|------|---------|
| **BUG** | Wrong CLI token emitted — current code is actively broken |
| **P1** | Missing option for a common use-case, easy to add |
| **P2** | Missing mode or feature, moderate effort |
| **P3** | Advanced / niche feature, complex effort |
| **DESIGN** | Intentional simplification; document limitation, use `Raw()` |

**Last updated**: 2026-03-22 (bugs fixed 2026-03-22; Output API implemented 2026-03-22; all Batch 1–5 gaps closed 2026-03-22; P2 deferred gaps closed 2026-03-22)

**Deployment status**: The pyvartools Python API has **not yet been deployed** and
will not be released until it is fully featured.  As a result:

- **Backward compatibility is not a concern** during development.  Parameters,
  method signatures, and CLI token patterns may be freely changed at any time
  without deprecation warnings or compatibility shims.
- **Work order priority**: favour completeness and internal consistency over
  minimal-touch fixes.  If fixing a gap naturally suggests cleaning up an adjacent
  inconsistency, do both in the same pass.
- **Design gaps** (DESIGN severity) should be resolved now rather than deferred —
  once the API ships, redesigning parameter signatures becomes a breaking change.

---

## Summary table

| Command | Status | BUG | P1 | P2 | P3 | DESIGN |
|---------|--------|-----|----|----|----|--------|
| `LS` | ✓ Complete | — | — | — | — | — |
| `aov` | ✓ Complete | — | — | — | — | — |
| `aov_harm` | ✓ Complete | — | ~~1~~ | — | — | — |
| `PDM` | ✓ Complete (new in 1.7) | — | — | — | — | — |
| `FTP` | ✓ Complete (new in 1.7) | — | — | — | — | — |
| `MatchedFilter` | ✓ Complete (new in 1.7) | — | — | — | — | — |
| `BLS` | ✓ Complete | — | ~~1~~ | ~~6~~ | ~~3~~ | — |
| `BLSFixPer` | ✓ Complete | — | — | ~~1~~ | — | — |
| `BLSFixPerDurTc` | ✗ Not implemented | — | — | ~~3~~ | — | — |
| `BLSFixDurTc` | ✗ Not implemented | — | — | — | — | — |
| `Phase` | ✓ Complete | — | — | — | — | — |
| `autocorrelation` | ✓ Complete | — | — | — | — | — |
| `dftclean` | ✓ Complete | — | — | ~~4~~ | — | — |
| `wwz` | ✓ Complete | — | — | ~~2~~ | — | — |
| `GetLSAmpThresh` | ✓ Complete | — | — | ~~1~~ | — | — |
| `Killharm` | ✓ Complete | — | — | ~~2~~ | — | — |
| `linfit` | ✓ Complete | — | ~~1~~ | ~~2~~ | — | — |
| `nonlinfit` | ✓ Complete | — | ~~2~~ | ~~4~~ | ~~3~~ | — |
| `Injectharm` | ✓ Complete (DESIGN) | — | — | — | — | ~~2~~ |
| `Injecttransit` | ✓ Complete | — | — | ~~3~~ | — | — |
| `resample` | ✓ Complete | — | ~~2~~ | ~~2~~ | — | — |
| `TFA_SR` | ✓ Complete | — | — | ~~2~~ | — | — |
| `SYSREM` | ✓ Complete | — | — | — | — | — |
| `MandelAgolTransit` | ✓ Complete | — | ~~1~~ | ~~1~~ | ~~1~~ | — |
| `SoftenedTransit` | ✓ Complete | — | — | ~~1~~ | — | — |
| `Starspot` | ✓ Complete | — | ~~1~~ | — | — | — |
| `microlens` | ✓ Complete | — | — | ~~2~~ | — | — |
| `addnoise` | ✓ Complete | — | — | ~~3~~ | — | — |
| `findblends` | ✓ Complete | — | ~~1~~ | ~~3~~ | — | — |
| `match` | ✓ Complete | — | ~~1~~ | — | — | — |
| `o` | ✓ Complete | — | ~~2~~ | ~~4~~ | ~~1~~ | — |
| `clip` | ✓ Complete | — | — | — | — | — |
| `rms` | ✓ Complete | — | — | — | — | — |
| `rmsbin` | ✓ Complete | — | — | — | — | — |
| `chi2` | ✓ Complete | — | — | — | — | — |
| `chi2bin` | ✓ Complete | — | — | — | — | — |
| `alarm` | ✓ Complete | — | — | — | — | — |
| `vonNeumann` | ✓ Complete (new in 1.7) | — | — | — | — | — |
| `percentileratios` | ✓ Complete (new in 1.7) | — | — | — | — | — |
| `rescalesig` | ✓ Complete | — | — | — | — | — |
| `ensemblerescalesig` | ✓ Complete | — | — | — | — | — |
| `stats` | ✓ Complete | — | — | — | — | — |
| `sortlc` | ✓ Complete | — | — | — | — | — |
| `restricttimes` | ✓ Complete | — | — | — | — | — |
| `restoretimes` | ✓ Complete | — | — | — | — | — |
| `savelc` | ✓ Complete | — | — | — | — | — |
| `restorelc` | ✓ Complete | — | — | — | — | — |
| `difffluxtomag` | ✓ Complete | — | — | — | — | — |
| `fluxtomag` | ✓ Complete | — | — | — | — | — |
| `magtoflux` | ✓ Complete | — | — | — | — | — |
| `changeerror` | ✓ Complete | — | — | — | — | — |
| `changevariable` | ✓ Complete | — | — | — | — | — |
| `copylc` | ✓ Complete | — | — | — | — | — |
| `medianfilter` | ✓ Complete | — | — | — | — | — |
| `expr` | ✓ Complete | — | — | — | — | — |
| `print_cols` | ✓ Complete | — | — | — | — | — |
| `FFT` / `IFFT` | ✓ Complete | — | — | — | — | — |
| `decorr` | ✓ Complete | — | — | — | — | — |
| `Jstet` | ✓ Complete | — | — | — | — | — |
| `TFA` | ✓ Complete | — | — | — | — | — |
| `addfitskeyword` | ✓ Complete | — | — | — | — | — |
| `R` | ✓ Complete | — | — | — | — | — |
| `ifcmd` | ✓ Complete | — | — | — | — | — |
| `columnsuffix` | ✓ Complete | — | — | — | — | — |
| `Raw` | ✓ Complete | — | — | — | — | — |

---

## Priority queue

### Bugs — **All fixed (2026-03-22)**

~~1. **`linfit`** — emitted `maskpoints` but CLI uses `fitmask`~~  ✓ Fixed: parameter renamed to `fitmask`; token corrected.

~~2. **`nonlinfit`** — emitted `maskpoints` but CLI uses `fitmask`~~  ✓ Fixed: same.

~~3. **`MandelAgolTransit`** — spurious `fitP` parameter that does not exist in CLI~~  ✓ Fixed: `fitP` removed; `fitephem` covers period+epoch jointly.

~~4. **`autocorrelation`** — `save_result=False` omitted required outdir~~  ✓ Fixed: `save_result` parameter removed; outdir always passed; result always captured.

Regression tests added for all four fixes (`test_linfit_fitmask_emits_fitmask_not_maskpoints`, `test_nonlinfit_fitmask_emits_fitmask_not_maskpoints`, `test_mandel_agol_no_fitp_param`, `test_autocorrelation_no_save_result_param`, `test_autocorrelation_always_includes_outdir`, `test_autocorrelation_file_capture_single`, `test_autocorrelation_file_capture_batch`).

### P1 — Missing common options (all fixed 2026-03-22)

| # | Command | Missing | Status |
|---|---------|---------|--------|
| ~~1~~ | ~~`aov_harm`~~ | ~~`nharm` `var`/`expr` forms~~ | ✓ Fixed |
| ~~2~~ | ~~`autocorrelation`~~ | ~~outdir always required; `save_result=False` is broken~~ | ✓ Fixed (Output API) |
| ~~3~~ | ~~`BLS`~~ | ~~`var`/`expr` for minper/maxper/rmin/rmax/nbins/subsample~~ | ✓ Fixed |
| ~~4~~ | ~~`linfit`~~ | ~~`modelvar` option~~ | ✓ Fixed |
| ~~5~~ | ~~`nonlinfit`~~ | ~~`modelvar` option~~ | ✓ Fixed |
| ~~6~~ | ~~`nonlinfit`~~ | ~~`reject`/`useMAD`/`iter`/`fixednum` options~~ | ✓ N/A — CLI has no reject for nonlinfit |
| ~~7~~ | ~~`SYSREM`~~ | ~~`save_trends` passes outdir; CLI expects a file path~~ | ✓ Fixed (Output API) |
| ~~8~~ | ~~`MandelAgolTransit`~~ | ~~`modelvar` option~~ | ✓ Fixed |
| ~~9~~ | ~~`resample`~~ | ~~`left`/`right` boundary values for `spline`~~ | ✓ Fixed |
| ~~10~~ | ~~`resample`~~ | ~~`nbreaks`/`order` for `bspline`~~ | ✓ Fixed |
| ~~11~~ | ~~`Starspot`~~ | ~~Named parameter API for `initial_params` / `fit_flags`~~ | ✓ Fixed |
| ~~12~~ | ~~`findblends`~~ | ~~Verify/fix whether `"bls"` is a valid period source~~ | ✓ Fixed — `"bls"` invalid; default changed to `"list"` |
| ~~13~~ | ~~`match`~~ | ~~Clarify/fix `"inlist"` API~~ | ✓ Fixed — `catalog` param + `inlist_column` kwarg |
| ~~14~~ | ~~`o`~~ | ~~Missing `fits` flag~~ | ✓ Fixed |
| ~~15~~ | ~~`o`~~ | ~~Missing `noclobber` flag~~ | ✓ Fixed |

### P2 — Missing modes / features

| # | Command | Missing | Status |
|---|---------|---------|--------|
| ~~1~~ | ~~`BLS`~~ | ~~`q` transit depth mode (qmin/qmax)~~ | ✓ Fixed |
| ~~2~~ | ~~`BLS`~~ | ~~`density` mode~~ | ✓ Fixed |
| ~~3~~ | ~~`BLS`~~ | ~~`df` frequency step mode~~ | ✓ Fixed |
| ~~4~~ | ~~`BLS`~~ | ~~`extraparams` flag~~ | ✓ Fixed |
| ~~5~~ | ~~`BLS`~~ | ~~`nobinnedrms` flag~~ | ✓ Fixed |
| ~~6~~ | ~~`BLS`~~ | ~~`ophcurve` / `ojdcurve` outputs~~ | ✓ Fixed |
| ~~7~~ | ~~`BLSFixPer`~~ | ~~`q` transit depth mode~~ | ✓ Fixed |
| ~~8~~ | ~~`BLSFixPerDurTc`~~ | ~~`fixdepth` option~~ | ✗ N/A — binary returns "Invalid command" |
| ~~9~~ | ~~`BLSFixPerDurTc`~~ | ~~`ophcurve` / `ojdcurve` outputs~~ | ✗ N/A — binary returns "Invalid command" |
| ~~10~~ | ~~`dftclean`~~ | ~~`finddirtypeaks`~~ | ✓ Fixed |
| ~~11~~ | ~~`dftclean`~~ | ~~`outcbeam`~~ | ✓ Fixed |
| ~~12~~ | ~~`dftclean`~~ | ~~`useampspec` / `verboseout` flags~~ | ✓ Fixed |
| ~~13~~ | ~~`wwz`~~ | ~~`outfulltransform` `fits`/`pm3d` format; `format` naming option~~ | ✓ Fixed |
| ~~14~~ | ~~`GetLSAmpThresh`~~ | ~~`"file" listfile` mode~~ | ✓ Fixed |
| ~~15~~ | ~~`Killharm`~~ | ~~Output format options (`outampphase` etc.)~~ | ✓ Fixed |
| ~~16~~ | ~~`Killharm`~~ | ~~`clip` option~~ | ✓ Fixed |
| ~~17~~ | ~~`linfit`~~ | ~~`reject` option~~ | ✓ Fixed |
| ~~18~~ | ~~`linfit`~~ | ~~`format` option for model naming~~ | ✓ Fixed |
| ~~19~~ | ~~`nonlinfit`~~ | ~~`linfit` companion parameter list~~ | ✓ Fixed |
| ~~20~~ | ~~`nonlinfit`~~ | ~~`errors` option~~ | ✓ Fixed |
| ~~21~~ | ~~`nonlinfit`~~ | ~~`amoeba` `tolerance` / `maxsteps`~~ | ✓ Fixed |
| ~~22~~ | ~~`nonlinfit`~~ | ~~`reject` / `useMAD` / `iter` / `fixednum`~~ | ✓ N/A — CLI has no reject for nonlinfit |
| ~~23~~ | ~~`Injecttransit`~~ | ~~Per-parameter `list`/`rand`/`logrand`/`expr` modes~~ | ✓ Fixed |
| ~~24~~ | ~~`Injecttransit`~~ | ~~`hk` eccentricity parameterization~~ | ✓ Fixed |
| ~~25~~ | ~~`Injecttransit`~~ | ~~`dilute` option~~ | ✓ Fixed |
| ~~26~~ | ~~`resample`~~ | ~~`file` mode (times from a file)~~ | ✓ Fixed |
| ~~27~~ | ~~`resample`~~ | ~~`gaps` option~~ | ✓ Fixed |
| ~~28~~ | ~~`TFA_SR`~~ | ~~`decorr` option (simultaneous EPD)~~ | ✓ Fixed |
| ~~29~~ | ~~`TFA_SR`~~ | ~~`bin`/`harm` signal modes missing `period` sub-option~~ | ✓ Fixed |
| ~~30~~ | ~~`MandelAgolTransit`~~ | ~~`"b" bimpact` alternative to `"i" inclination`~~ | ✓ Fixed |
| ~~31~~ | ~~`SoftenedTransit`~~ | ~~`fit_harm` currently hardcoded to `0`~~ | ✓ Fixed |
| ~~32~~ | ~~`microlens`~~ | ~~`list`/`fixcolumn` per-parameter modes~~ | ✓ Fixed |
| ~~33~~ | ~~`microlens`~~ | ~~`step` and `novary` per-parameter options~~ | ✓ Fixed |
| ~~34~~ | ~~`addnoise`~~ | ~~`matern` noise model~~ | ✓ Fixed |
| ~~35~~ | ~~`addnoise`~~ | ~~`wavelet` noise model~~ | ✓ Fixed |
| ~~36~~ | ~~`addnoise`~~ | ~~`bintime` for `squareexp`/`exp`~~ | ✓ Fixed |
| ~~37~~ | ~~`findblends`~~ | ~~`xycol`, `starlist`, `zeromag`, `nofluxconvert`~~ | ✓ Fixed |
| ~~38~~ | ~~`o`~~ | ~~`copyheader`, `namecommand`, `namefromlist`, `delimiter` flags~~ | ✓ Fixed |

### P3 — Advanced / niche features (all fixed 2026-03-22)

| # | Command | Missing | Status |
|---|---------|---------|--------|
| ~~1~~ | ~~`BLS`~~ | ~~`stepP` / `steplogP` frequency grid modes~~ | ✓ Fixed |
| ~~2~~ | ~~`BLS`~~ | ~~`adjust-qmin-by-mindt` / `reduce-nbins`~~ | ✓ Fixed |
| ~~3~~ | ~~`BLS`~~ | ~~`reportharmonics` flag~~ | ✓ Fixed |
| ~~4~~ | ~~`nonlinfit`~~ | ~~`mcmc` optimizer (with all sub-options)~~ | ✓ Fixed |
| ~~5~~ | ~~`nonlinfit`~~ | ~~`covariance` option (GP noise models)~~ | ✓ Fixed |
| ~~6~~ | ~~`nonlinfit`~~ | ~~`priors` / `constraints` options~~ | ✓ Fixed |
| ~~7~~ | ~~`MandelAgolTransit`~~ | ~~RV fitting options~~ | ✓ Fixed |
| ~~8~~ | ~~`o`~~ | ~~`logcommandline` flag~~ | ✓ Fixed |

### DESIGN — Intentional simplifications

| Command | Limitation | Workaround |
|---------|------------|------------|
| `Injectharm` | Only `ampfix`/`phasefix`; same for all harmonics | Use `Raw()` |
| `Injectharm` | Only `fix` period mode | Use `Raw()` |

---

## Detailed gap sections

---

### `aov_harm` — ✓ Fixed 2026-03-22

~~`nharm` accepts `int` only; CLI also accepts `var`/`expr` forms~~

**Fixed**: `nharm: Union[int, str]`; `_varexpr(self.nharm)` emitted.
Tests: `test_aov_harm_nharm_var`, `test_aov_harm_nharm_expr`.

---

### `BLS` — ✓ Fixed 2026-03-22

All gaps closed.  Added: `var`/`expr` for all numeric params; `qmin`/`qmax`; `density_mode`; `df`; `extraparams`; `nobinnedrms`; `save_phcurve`/`save_jdcurve`; `freq_grid`; `adjust_qmin`; `reduce_nbins`; `reportharmonics`.
Tests: `test_bls_minper_var`, `test_bls_rmin_var`, `test_bls_nbins_var`, `test_bls_qmin_qmax`, `test_bls_density_mode`, `test_bls_df`, `test_bls_extraparams`, `test_bls_nobinnedrms`, `test_bls_save_phcurve`, `test_bls_save_jdcurve`, `test_bls_freq_grid_steplogP`, `test_bls_adjust_qmin`, `test_bls_reduce_nbins`, `test_bls_reportharmonics`.

---

### `BLSFixPer` — ✓ Fixed 2026-03-22

~~Missing `q` transit depth mode (qmin/qmax)~~ — Fixed: `qmin`/`qmax` params added; emits `"q" qmin qmax` instead of `"r" rmin rmax` when set.
Test: `test_blsfixper_qmin_qmax`.

---

### `BLSFixPerDurTc` — ✗ Not implemented

`vartools -help -BLSFixPerDurTc` returns "Invalid command" in the installed binary.  `_to_cli_args()` raises `NotImplementedError`.  Use `Raw()` to pass tokens manually if needed.

---

### `BLSFixDurTc` — ✗ Not implemented

`vartools -help -BLSFixDurTc` returns "Invalid command" in the installed binary.  `_to_cli_args()` raises `NotImplementedError`.  Use `Raw()` to pass tokens manually if needed.

---

### `autocorrelation` — ✓ Updated 2026-03-22

`save_result` parameter **restored** with the four-mode Output API
(see [Output API](#output-api--four-mode-save_-control) below).

- `save_result=True` (default) — captures the result as
  `result.files["autocorrelation_result_N"]`.
- `save_result=False` — file still written (CLI always does so) but not
  captured into Python.
- `save_result="/path"` or `save_result=Output("/path", capture=True)` — write
  to the specified directory, optionally capturing.

The class-level `_mandatory_output = True` attribute marks this command so
`_has_output_reqs()` forces subprocess mode even when `save_result=False`.

End-to-end tests added: `test_autocorrelation_default_captured`,
`test_autocorrelation_no_capture`.

---

### `dftclean` — ✓ Fixed 2026-03-22

All P2 gaps closed: `finddirtypeaks`, `finddirtypeaks_clip`, `finddirtypeaks_clipiter`, `outcbeam`, `useampspec`, `verboseout`.
Tests: `test_dftclean_finddirtypeaks`, `test_dftclean_finddirtypeaks_clip`, `test_dftclean_finddirtypeaks_clipiter`, `test_dftclean_outcbeam`, `test_dftclean_useampspec`, `test_dftclean_verboseout`.

---

### `wwz` — ✓ Fixed 2026-03-22

All P2 gaps closed: `transform_format` (`"fits"`/`"pm3d"`), `transform_name` (format), `maxtransform_name` (format).
Tests: `test_wwz_transform_format_fits`, `test_wwz_transform_format_pm3d`, `test_wwz_transform_name_format`, `test_wwz_maxtransform_name_format`.

---

### `GetLSAmpThresh` — ✓ Fixed 2026-03-22

~~Missing `"file" listfile` mode~~ — Fixed: `mode="file"` param added; emits `"file" listfile` when set.
Test: `test_getlsampthresh_file_mode`.

---

### `Killharm` — ✓ Fixed 2026-03-22

All P2 gaps closed: `output_format` (`"outampphase"` etc.) and `clip`.
Tests: `test_killharm_outampphase`, `test_killharm_outampradphase`, `test_killharm_clip`.

---

### `linfit` — ✓ Fixed 2026-03-22

All gaps closed: `modelvar`, `reject`/`reject_usemad`/`reject_iter`/`reject_fixednum`, `model_nameformat`.
Tests: `test_linfit_modelvar`, `test_linfit_model_nameformat`, `test_linfit_reject`, `test_linfit_reject_usemad`, `test_linfit_reject_iter`, `test_linfit_reject_fixednum`.

---

### `nonlinfit` — ✓ Fixed 2026-03-22

All P1/P2/P3 gaps closed: `modelvar`, `errors`, `linfit_params`, `covariance`, `priors`, `constraints`; `amoeba` `tolerance`/`maxsteps`; `mcmc` with full sub-option set (`Naccept`, `Nlinkstotal`, `fracburnin`, `eps`, `skipamoeba`, `maxmemstore`, `outchains` with `format`/`printevery`); `model_nameformat`.

Note: CLI help for nonlinfit shows no `reject` option (unlike linfit) — P1 items 5/6 (reject) were not applicable.

Token order rewritten: `function paramlist [linfit] [errors] [covariance] [priors] [constraints] optimizer [opts] [correctlc] [omodel [format] [modelvar]] [fitmask]`.

Tests: `test_nonlinfit_modelvar`, `test_nonlinfit_errors`, `test_nonlinfit_linfit_params`, `test_nonlinfit_amoeba_tolerance`, `test_nonlinfit_amoeba_maxsteps`, `test_nonlinfit_mcmc_naccept`, `test_nonlinfit_mcmc_nlinkstotal`, `test_nonlinfit_mcmc_fracburnin`, `test_nonlinfit_mcmc_skipamoeba`, `test_nonlinfit_covariance`, `test_nonlinfit_priors`, `test_nonlinfit_constraints`, `test_nonlinfit_mcmc_model_nameformat`.

---

### `Injectharm` — DESIGN (✓ Documented 2026-03-22)

| Gap | Severity | Note |
|-----|----------|------|
| Only `ampfix`/`phasefix`; no per-harmonic amplitude modes | DESIGN | CLI: `amplist`, `amprand`, `amplogrand`, `amprel`, `phaselist`, `phaserand`, `phaserel` |
| Only `fix` period mode; no `rand`/`logrand`/`randfreq`/`lograndfreq` | DESIGN | CLI: `rand`, `logrand`, `randfreq`, `lograndfreq` |

**Resolution**: Documented as intentional simplifications in the class docstring.  Use `Raw()`
for full control.  The parameter space is large and the combinatorics get
unwieldy in Python — these restrictions keep the class manageable.

---

### `Injecttransit` — ✓ Fixed 2026-03-22

All P2 gaps closed.  Added `_injectparam` helper for per-parameter float/string
passthrough; `hk` mode (`hk=True`, `h`, `k` params); `dilute` float/string;
full string passthrough for all positional params (`period`, `Rp`, `Mp`,
`phase`, `sini`, `Mstar`, `Rstar`, `e`, `omega`).

~~| No per-parameter `list`/`rand`/`logrand`/`expr` modes | P2 | ... |~~
~~| Missing `hk` eccentricity parameterization | P2 | ... |~~
~~| Missing `dilute` option | P2 | ... |~~

Tests: `test_injecttransit_float_params`, `test_injecttransit_string_params`,
`test_injecttransit_hk_mode`, `test_injecttransit_dilute_float`,
`test_injecttransit_dilute_string`.

---

### `resample` — ✓ Fixed 2026-03-22

All gaps closed.  P1 gaps (`left`/`right`, `nbreaks`/`order`) fixed previously.
P2 gaps now fixed: `file_times` (path → `"file" "fix" path ["column" N]`; string
starting with `"list"` → `"file" "list" ...`); `gaps` string passthrough.

Tests: `test_resample_file_mode_fix`, `test_resample_file_mode_column`,
`test_resample_file_mode_list`, `test_resample_gaps`.

---

### `TFA_SR` — ✓ Fixed 2026-03-22

All P2 gaps closed.  Added `decorr_params: Optional[str]` (raw token string
passed after `"decorr"` keyword) and `signal_period: Optional[Union[float, str]]`
(emits `"period" val` after the signal mode tokens; supports float and string
passthrough via `_pval`).

Tests: `test_tfa_sr_decorr`, `test_tfa_sr_signal_period_bin`,
`test_tfa_sr_signal_period_harm`.

---

### `SYSREM` — ✓ Fixed 2026-03-22

~~`save_trends` passes outdir; CLI expects a file path~~ — **Fixed**: `save_trends`
now uses `_norm_save` / `_should_emit` and appends `["otrends", actual_outdir]`
only when emitting is requested.  `save_trends` is also registered in
`_output_file_specs()` as `".sysrem.trends"` so the file is properly captured
as `result.files["SYSREM_trends_N"]` when `save_trends=True`.

---

### `MandelAgolTransit` — ✓ Fixed 2026-03-22

All gaps closed: `modelvar`, `bimpact` (impact parameter mode), RV fitting (`rv_file`, `rv_model_file`, `K0`, `gamma0`, `fitK`, `fitgamma`).
Tests: `test_mandel_agol_modelvar`, `test_mandel_agol_default_inclination_token`, `test_mandel_agol_bimpact`, `test_mandel_agol_rv_fitting`.

---

### `SoftenedTransit` — ✓ Fixed 2026-03-22

~~`fit_harm` hardcoded to `0`~~ — Fixed: `fit_harm: int = 0`, `fit_harm_method`, `fit_harm_nharm`, `fit_harm_nsubharm` params added.
Tests: `test_softenedtransit_fit_harm_zero`, `test_softenedtransit_fit_harm`.

---

### `Starspot` — ✓ Fixed 2026-03-22

~~Positional `initial_params`/`fit_flags` lists~~ — Fixed: replaced with named kwargs `a0`, `b0`, `alpha0`, `i0`, `chi0`, `psi00`, `mconst0`, `fit_period`, `fit_a`, `fit_b`, `fit_alpha`, `fit_i`, `fit_chi`, `fit_psi`, `fit_mconst`.  Old `initial_params`/`fit_flags` kwargs raise `TypeError`.

CLI fit-flag order confirmed as `fitP fita fitb fitalpha fiti fitchi fitpsi fitmconst`.
Tests: `test_starspot_named_params`, `test_starspot_fit_flags_named`, `test_starspot_old_initial_params_raises`.

---

### `microlens` — ✓ Fixed 2026-03-22

All P2 gaps closed.  String passthrough for all per-parameter values (handles
`"fixcolumn colname"`, `"list column 3"`, etc.).  Added per-parameter
`{name}_step: Optional[float]` and `{name}_novary: bool` kwargs for all five
parameters (`f0`, `f1`, `u0`, `t0`, `tmax`).

Tests: `test_microlens_step`, `test_microlens_novary`, `test_microlens_fixcolumn`.

---

### `addnoise` — ✓ Fixed 2026-03-22

All P2 gaps closed: `matern` (`nu` param), `wavelet` (`gamma` param), `bintime` for `squareexp`/`exp`.
Tests: `test_addnoise_matern`, `test_addnoise_wavelet`, `test_addnoise_bintime`.

---

### `findblends` — ✓ Fixed 2026-03-22

All gaps closed: default period changed from invalid `"bls"` to `"list"`; added `xycol`, `starlist`, `zeromag`, `nofluxconvert`.
Tests: `test_findblends_period_fix`, `test_findblends_xycol`, `test_findblends_starlist`, `test_findblends_zeromag`, `test_findblends_nofluxconvert`.

---

### `match` — ✓ Fixed 2026-03-22

~~`"inlist"` API confusion~~ — Fixed: first positional param renamed `filename` → `catalog`; added `inlist_column: str | int | None`; when `source="inlist"`, emits `inlist_column` instead of `catalog`; raises `ValueError` when `inlist_column` is missing.
Tests: `test_match_inlist_column_emits_column_not_filename`, `test_match_inlist_missing_column_raises`.

---

### `o` (output command) — ✓ Fixed 2026-03-22

All gaps closed: `fits`, `noclobber`, `copyheader`, `namecommand`, `namefromlist` (bool → `["namefromlist"]` or str → `["namefromlist" "column" col]`), `delimiter`, `logcommandline`.
Tests: `test_o_fits`, `test_o_noclobber`, `test_o_copyheader`, `test_o_namecommand`, `test_o_namefromlist_bool`, `test_o_namefromlist_column`, `test_o_delimiter`, `test_o_logcommandline`.

---

---

## Output API — four-mode `save_*` control

**Implemented 2026-03-22.**

Every `save_*` parameter now accepts four modes via the `Output` class or shorthand forms:

| Value | Mode | Written to disk? | Captured in `result.files`? |
|-------|------|-----------------|------------------------------|
| `False` (default) | 4 — suppress | no | no |
| `True` | 1 — temp, capture | temp dir (auto-deleted) | **yes** |
| `"/path/to/dir"` | 3 — disk only | that directory | no |
| `Output("/path/to/dir", capture=True)` | 2 — disk + capture | that directory | **yes** |

### Files created/modified

| File | Change |
|------|--------|
| `pyvartools/_output.py` | **New** — `Output(path, capture)` class |
| `pyvartools/__init__.py` | Export `Output` |
| `pyvartools/commands/_helpers.py` | Added `_norm_save()`, `_should_emit()`; updated `_outtoken()` (backward-compatible) |
| `pyvartools/_command.py` | `_requested_outputs()` uses `_should_emit()` + `_mandatory_output` |
| `pyvartools/pipeline.py` | `_build_cmd()` builds per-output `_outdir_map`; `_collect_output_files()` gates on `spec.capture` |
| `pyvartools/commands/periodicity.py` | `autocorrelation` restored `save_result` + `_mandatory_output=True`; `dftclean`, `wwz` custom tokens |
| `pyvartools/commands/fitting.py` | `MandelAgolTransit`, `microlens`, `nonlinfit`, `findblends`, `SYSREM` updated; `findblends` + `SYSREM.trends` gained `_output_file_specs()` |
| `pyvartools/commands/manipulation.py` | `linfit`, `decorr` updated; both gained `_output_file_specs()` |

### Previously unregistered output specs — now fixed

The following commands had working `save_*` parameters (CLI tokens were emitted correctly) but the files were never captured into `result.files` because `_output_file_specs()` was not implemented.  This is now fixed:

| Command | Param | Key in `result.files` | Suffix |
|---------|-------|-----------------------|--------|
| `linfit` | `save_model` | `linfit_model_N` | `.linfit.model` |
| `decorr` | `save_model` | `decorr_model_N` | `.decorr.model` |
| `findblends` | `save_matches` | `findblends_matches_N` | `.findblends.matches` |
| `SYSREM` | `save_trends` | `SYSREM_trends_N` | `.sysrem.trends` |

!!! note "Suffix verification pending"
    The suffixes `.linfit.model`, `.decorr.model`, `.findblends.matches`, and
    `.sysrem.trends` are inferred from the CLI token patterns.  They should be
    verified against live vartools output before being relied upon in production.

---

## Notes on test coverage

Regression tests were added alongside the bug fixes (2026-03-22):

| Test | Verifies |
|------|---------|
| `test_linfit_fitmask_emits_fitmask_not_maskpoints` | `fitmask` token emitted, `maskpoints` absent |
| `test_nonlinfit_fitmask_emits_fitmask_not_maskpoints` | same for nonlinfit |
| `test_mandel_agol_no_fitp_param` | `fitP` kwarg now raises `TypeError` |
| ~~`test_autocorrelation_no_save_result_param`~~ | *Replaced* — `save_result` is now a valid parameter (Output API); this test was removed |
| `test_autocorrelation_always_includes_outdir` | outdir always present in CLI args |
| `test_autocorrelation_file_capture_single` | result file captured in single-LC run |
| `test_autocorrelation_file_capture_batch` | result files captured (one per LC) in batch run |

Output API tests added (2026-03-22) in `TestOutputAPICLI` and `TestOutputAPIEndToEnd`:

| Test | Verifies |
|------|---------|
| `test_ls_save_false_suppresses` | Mode 4: `"0"` token, no outdir |
| `test_ls_save_true_mode1` | Mode 1: `"1"` token + fallback outdir |
| `test_ls_save_path_mode3` | Mode 3: `"1"` token + user path |
| `test_ls_save_output_mode2` | Mode 2: `Output(path, capture=True)` |
| `test_norm_save_roundtrips` | `_norm_save()` handles all input forms |
| `test_should_emit` | `_should_emit()` returns correct bool for all modes |
| `test_dftclean_save_path` | custom `outdspec` token uses user path |
| `test_wwz_save_path` | custom `outfulltransform` token uses user path |
| `test_microlens_save_path` | `omodel` token uses user path |
| `test_nonlinfit_save_false_no_omodel` | Mode 4: no `omodel` token |
| `test_linfit_save_path` | `omodel` token uses user path |
| `test_decorr_save_false_emits_zero` | Mode 4: `"0"` sentinel for decorr |
| `test_sysrem_save_trends_path` | `otrends` token uses user path |
| `test_findblends_save_matches_path` | `omatches` token uses user path |
| `test_mandel_agol_save_phcurve_path` | `ophcurve` token uses user path |
| `test_autocorrelation_mandatory_output_class_attr` | `_mandatory_output = True` on `autocorrelation` |
| `test_linfit_output_file_specs` | `linfit._output_file_specs()` includes `"model"` |
| `test_decorr_output_file_specs` | `decorr._output_file_specs()` includes `"model"` |
| `test_findblends_output_file_specs` | `findblends._output_file_specs()` includes `"matches"` |
| `test_sysrem_output_file_specs_has_trends` | `SYSREM._output_file_specs()` includes `"trends"` |
| `test_output_repr` | `Output.__repr__` works |
| `test_ls_mode1_capture_only` (E2E) | Mode 1 unchanged from prior behaviour |
| `test_ls_mode3_writes_file_not_captured` (E2E) | Mode 3: file on disk, not in `result.files` |
| `test_ls_mode2_writes_and_captures` (E2E) | Mode 2: file on disk and in `result.files` |
| `test_ls_mode4_no_file_no_capture` (E2E) | Mode 4: no file, no capture |
| `test_autocorrelation_default_captured` (E2E) | `save_result=True` (default): captured |
| `test_autocorrelation_no_capture` (E2E) | `save_result=False`: not captured |
