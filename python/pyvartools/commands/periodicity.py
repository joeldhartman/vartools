"""Periodicity-search command wrappers."""

from __future__ import annotations
from typing import List, Optional, Union

from pyvartools._command import VartoolsCommand
from ._helpers import (_auto_or_varexpr, _bool, _flag, _fixperiodsnr_tokens,
                       _norm_save, _outtoken, _period_spec, _pval,
                       _should_emit, _varexpr)


class LS(VartoolsCommand):
    """Lomb-Scargle (GLS) periodogram.

    Parameters
    ----------
    minp, maxp : float or str
        Minimum and maximum period to search.  Can be:

        - A number — passed directly to vartools as a fixed value.
        - A bare identifier string (e.g. ``"minperiod"``) — vartools reads
          the value from a named variable (``var minperiod``).
        - Any other string (e.g. ``"tspan/100"``) — vartools evaluates it as
          a math expression (``expr tspan/100``).
    subsample : float or str
        Frequency step as a fraction of 1/T (time span).  Accepts the same
        number / variable-name / expression forms as *minp* and *maxp*.
    npeaks : int
        Number of peaks to report.
    save_periodogram : bool
        Write the periodogram spectrum to a file in the pipeline output dir.
    noGLS : bool
        Use the classical (non-generalised) Lomb-Scargle statistic.
    whiten : bool
        Iteratively whiten the light curve.
    clip : float, optional
        Sigma-clipping threshold applied during whitening.
    clipiter : int, optional
        Number of sigma-clipping iterations.
    bootstrap : int, optional
        Number of bootstrap resamples for false-alarm probability estimation.
    maskpoints : str, optional
        Name of a mask variable; masked points are excluded.
    fixperiod_snr : float, int, str, or None, optional
        Evaluate the periodogram at a known/fixed period and report its
        significance.  Forms:

        - A number (e.g. ``1.234``) — evaluate at that fixed period.
        - ``"ls"`` / ``"aov"`` / ``"injectharm"`` — use the best period found
          by a prior LS, AOV, or injection-recovery command.
        - ``"fixcolumn COLNAME"`` — read the period from a named per-star
          column.
        - ``"list"`` / ``"list column 2"`` — read the period from a list file
          column.

        When set, four extra output columns are appended:
        ``LS_PeriodFix_N``, ``Log10_LS_Prob_PeriodFix_N``,
        ``LS_Periodogram_Value_PeriodFix_N``, ``LS_SNR_PeriodFix_N``.
    """

    _vt_name = "LS"

    def __init__(
        self,
        minp: Union[float, str],
        maxp: Union[float, str],
        subsample: Union[float, str],
        npeaks: int = 5,
        save_periodogram: bool = False,
        noGLS: bool = False,
        whiten: bool = False,
        clip: Optional[float] = None,
        clipiter: Optional[int] = None,
        bootstrap: Optional[int] = None,
        maskpoints: Optional[str] = None,
        fixperiod_snr: Union[float, int, str, None] = None,
    ) -> None:
        # Numeric minp/maxp must be positive and ordered.  Variable
        # references / expressions / numpy arrays / PerLC are accepted
        # as-is because their values may only be known at run time.
        from pyvartools.perlc import PerLC
        import numpy as _np
        def _numeric(v) -> bool:
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                return True
            if isinstance(v, (_np.integer, _np.floating)):
                return True
            return False
        if _numeric(minp) and minp <= 0:
            raise ValueError(
                f"cmd.LS: minp must be > 0 (got {minp!r}); use a positive "
                f"period bound in days or pass a variable name / "
                f"expression / PerLC for per-LC search bounds."
            )
        if _numeric(maxp) and maxp <= 0:
            raise ValueError(
                f"cmd.LS: maxp must be > 0 (got {maxp!r})."
            )
        if _numeric(minp) and _numeric(maxp) and minp >= maxp:
            raise ValueError(
                f"cmd.LS: minp ({minp}) must be strictly less than "
                f"maxp ({maxp})."
            )
        self.minp = minp
        self.maxp = maxp
        self.subsample = subsample
        self.npeaks = npeaks
        self.save_periodogram = save_periodogram
        self.noGLS = noGLS
        self.whiten = whiten
        self.clip = clip
        self.clipiter = clipiter
        self.bootstrap = bootstrap
        self.maskpoints = maskpoints
        self.fixperiod_snr = fixperiod_snr

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = (["-LS"]
                + _varexpr(self.minp)
                + _varexpr(self.maxp)
                + _varexpr(self.subsample)
                + [str(self.npeaks)])
        args += _outtoken(self.save_periodogram, outdir)
        args += _bool("noGLS", self.noGLS)
        args += _bool("whiten", self.whiten)
        if self.clip is not None:
            args += ["clip", str(self.clip), str(self.clipiter or 3)]
        args += _fixperiodsnr_tokens(self.fixperiod_snr)
        if self.bootstrap is not None:
            args += ["bootstrap", str(self.bootstrap)]
        args += _flag("maskpoints", self.maskpoints)
        return args

    def _resolve_back_references(self, prev) -> None:
        from ._helpers import _resolve_period_backref
        self.fixperiod_snr = _resolve_period_backref(prev, self.fixperiod_snr)

    def _output_file_specs(self):
        return {"periodogram": (".ls", None)}


class aov(VartoolsCommand):
    """Analysis-of-Variance (AOV) periodogram.

    Parameters
    ----------
    minp, maxp : float or str
        Period search range.  Accepts the same number / variable-name /
        expression forms as :class:`LS` ``minp``/``maxp``.
    subsample : float or str
        Frequency step fraction.  Accepts var/expr forms.
    finetune : float or str
        Fine-tuning oversampling factor near peaks.  Accepts var/expr forms.
    npeaks : int
        Number of peaks to report.
    nbin : int or str, optional
        Number of phase bins.  If None, vartools chooses automatically.
        Accepts var/expr forms.
    save_periodogram : bool
        Write the periodogram spectrum to a file.
    whiten, clip, clipiter, uselog, maskpoints : see LS.
    """

    _vt_name = "aov"

    def __init__(
        self,
        minp: Union[float, str],
        maxp: Union[float, str],
        subsample: Union[float, str],
        finetune: Union[float, str],
        npeaks: int = 5,
        nbin: Optional[Union[int, str]] = None,
        save_periodogram: bool = False,
        whiten: bool = False,
        clip: Optional[float] = None,
        clipiter: Optional[int] = None,
        uselog: bool = False,
        maskpoints: Optional[str] = None,
        fixperiod_snr: Union[float, int, str, None] = None,
    ) -> None:
        self.minp = minp
        self.maxp = maxp
        self.subsample = subsample
        self.finetune = finetune
        self.npeaks = npeaks
        self.nbin = nbin
        self.save_periodogram = save_periodogram
        self.whiten = whiten
        self.clip = clip
        self.clipiter = clipiter
        self.uselog = uselog
        self.maskpoints = maskpoints
        self.fixperiod_snr = fixperiod_snr

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-aov"]
        if self.nbin is not None:
            args += ["Nbin"] + _varexpr(self.nbin)
        args += (_varexpr(self.minp) + _varexpr(self.maxp)
                 + _varexpr(self.subsample) + _varexpr(self.finetune)
                 + [str(self.npeaks)])
        args += _outtoken(self.save_periodogram, outdir)
        args += _bool("whiten", self.whiten)
        if self.clip is not None:
            args += ["clip", str(self.clip), str(self.clipiter or 3)]
        args += _fixperiodsnr_tokens(self.fixperiod_snr)
        args += _bool("uselog", self.uselog)
        args += _flag("maskpoints", self.maskpoints)
        return args

    def _resolve_back_references(self, prev) -> None:
        from ._helpers import _resolve_period_backref
        self.fixperiod_snr = _resolve_period_backref(prev, self.fixperiod_snr)

    def _output_file_specs(self):
        return {"periodogram": (".aov", None)}


class aov_harm(VartoolsCommand):
    """Harmonic AOV periodogram.

    Parameters
    ----------
    nharm : int
        Number of harmonics.
    minp, maxp : float or str
        Period search range.  Accepts var/expr forms (see :class:`LS`).
    subsample : float or str
        Frequency step fraction.  Accepts var/expr forms.
    finetune : float or str
        Fine-tuning oversampling factor.  Accepts var/expr forms.
    npeaks : int
        Number of peaks to report.
    save_periodogram : bool
    whiten, clip, clipiter, maskpoints : see LS.
    """

    _vt_name = "aov_harm"

    def __init__(
        self,
        nharm: Union[int, str],
        minp: Union[float, str],
        maxp: Union[float, str],
        subsample: Union[float, str],
        finetune: Union[float, str],
        npeaks: int = 5,
        save_periodogram: bool = False,
        whiten: bool = False,
        clip: Optional[float] = None,
        clipiter: Optional[int] = None,
        maskpoints: Optional[str] = None,
        fixperiod_snr: Union[float, int, str, None] = None,
    ) -> None:
        self.nharm = nharm
        self.minp = minp
        self.maxp = maxp
        self.subsample = subsample
        self.finetune = finetune
        self.npeaks = npeaks
        self.save_periodogram = save_periodogram
        self.whiten = whiten
        self.clip = clip
        self.clipiter = clipiter
        self.maskpoints = maskpoints
        self.fixperiod_snr = fixperiod_snr

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = (["-aov_harm"]
                + _varexpr(self.nharm)
                + _varexpr(self.minp) + _varexpr(self.maxp)
                + _varexpr(self.subsample) + _varexpr(self.finetune)
                + [str(self.npeaks)])
        args += _outtoken(self.save_periodogram, outdir)
        args += _bool("whiten", self.whiten)
        if self.clip is not None:
            args += ["clip", str(self.clip), str(self.clipiter or 3)]
        args += _fixperiodsnr_tokens(self.fixperiod_snr)
        args += _flag("maskpoints", self.maskpoints)
        return args

    def _resolve_back_references(self, prev) -> None:
        from ._helpers import _resolve_period_backref
        self.fixperiod_snr = _resolve_period_backref(prev, self.fixperiod_snr)

    def _output_file_specs(self):
        return {"periodogram": (".aov_harm", None)}


class BLS(VartoolsCommand):
    """Box-Least-Squares transit search.

    Parameters
    ----------
    minper, maxper : float or str
        Period search range.  Accepts var/expr forms (see :class:`LS`).
    rmin, rmax : float or str, optional
        Min/max fractional transit duration.  Used when ``qmin``/``qmax`` and
        ``density_mode`` are not set.  Accepts var/expr forms.
    qmin, qmax : float or str, optional
        Min/max fractional transit duration as fraction of orbit period.  When
        set, emits ``"q" qmin qmax`` instead of ``"r" rmin rmax``.
    density_mode : bool
        Use stellar density to set transit duration bounds.  When True,
        ``stellar_density``, ``min_exp_dur_frac``, ``max_exp_dur_frac`` must
        be set.
    stellar_density : float or str, optional
        Stellar density (g/cm³) for density mode.
    min_exp_dur_frac : float
        Minimum fraction of expected transit duration (density mode).
    max_exp_dur_frac : float
        Maximum fraction of expected transit duration (density mode).
    nbins : int or str
        Number of phase bins in the folded light curve.  Accepts var/expr forms.
    timezone : float
        Timezone offset (0 for HJD/BJD).
    npeaks : int
        Number of transit candidates to report.
    subsample : float or str
        Frequency oversampling factor (used with the ``"optimal"`` frequency
        mode).  Accepts var/expr forms.
    nfreq : int or str, optional
        Fixed number of frequencies (overrides subsample).  Accepts var/expr.
    df : float or str, optional
        Frequency step size (overrides both ``nfreq`` and ``subsample``).
    save_periodogram : bool
        Write the BLS power spectrum to a file.
    save_model : bool
        Write the best-fit transit model to a file.
    correct_lc : bool
        Subtract the best-fit transit from the light curve.
    extraparams : bool
        Compute and report additional transit parameters.
    fittrap : bool
        Fit a trapezoidal transit shape instead of a box.
    nobinnedrms : bool
        Do not report the binned-RMS diagnostic.
    save_phcurve : bool, str, or Output
        Write the phase-folded model light curve.
    ophcurve_phmin, ophcurve_phmax, ophcurve_phstep : float
        Phase range and step for the phase curve output.
    save_jdcurve : bool, str, or Output
        Write the model light curve in JD.
    ojdcurve_jdstep : float
        Time step for the JD curve output.
    freq_grid : str, optional
        ``"stepP"`` or ``"steplogP"`` frequency grid mode.
    adjust_qmin : bool
        Adjust qmin by the minimum time step.
    reduce_nbins : bool
        Reduce nbins by the minimum time step (requires ``adjust_qmin=True``).
    reportharmonics : bool
        Report harmonic periods.
    maskpoints : str, optional
    """

    _vt_name = "BLS"

    def __init__(
        self,
        minper: Union[float, str],
        maxper: Union[float, str],
        rmin: Union[float, str] = 0.01,
        rmax: Union[float, str] = 0.1,
        qmin: Union[float, str, None] = None,
        qmax: Union[float, str, None] = None,
        density_mode: bool = False,
        stellar_density: Union[float, str, None] = None,
        min_exp_dur_frac: float = 0.5,
        max_exp_dur_frac: float = 1.5,
        nbins: Union[int, str] = 200,
        timezone: float = 0,
        npeaks: int = 1,
        subsample: Union[float, str] = 1.0,
        nfreq: Union[int, str, None] = None,
        df: Union[float, str, None] = None,
        save_periodogram: bool = False,
        save_model: bool = False,
        correct_lc: bool = False,
        extraparams: bool = False,
        fittrap: bool = False,
        nobinnedrms: bool = False,
        save_phcurve=False,
        ophcurve_phmin: float = 0.0,
        ophcurve_phmax: float = 1.0,
        ophcurve_phstep: float = 0.005,
        save_jdcurve=False,
        ojdcurve_jdstep: float = 0.02,
        freq_grid: Optional[str] = None,
        adjust_qmin: bool = False,
        reduce_nbins: bool = False,
        reportharmonics: bool = False,
        maskpoints: Optional[str] = None,
    ) -> None:
        self.minper = minper
        self.maxper = maxper
        self.rmin = rmin
        self.rmax = rmax
        self.qmin = qmin
        self.qmax = qmax
        self.density_mode = density_mode
        self.stellar_density = stellar_density
        self.min_exp_dur_frac = min_exp_dur_frac
        self.max_exp_dur_frac = max_exp_dur_frac
        self.nbins = nbins
        self.timezone = timezone
        self.npeaks = npeaks
        self.subsample = subsample
        self.nfreq = nfreq
        self.df = df
        self.save_periodogram = save_periodogram
        self.save_model = save_model
        self.correct_lc = correct_lc
        self.extraparams = extraparams
        self.fittrap = fittrap
        self.nobinnedrms = nobinnedrms
        self.save_phcurve = save_phcurve
        self.ophcurve_phmin = ophcurve_phmin
        self.ophcurve_phmax = ophcurve_phmax
        self.ophcurve_phstep = ophcurve_phstep
        self.save_jdcurve = save_jdcurve
        self.ojdcurve_jdstep = ojdcurve_jdstep
        self.freq_grid = freq_grid
        self.adjust_qmin = adjust_qmin
        self.reduce_nbins = reduce_nbins
        self.reportharmonics = reportharmonics
        self.maskpoints = maskpoints

        # vartools' default "optimal" frequency grid only works when
        # density_mode=True (which uses stellar density to bound the
        # expected transit duration).  In r/q duration mode, the user
        # must specify nfreq= or df= explicitly.  Catch the invalid
        # combination at construction so the error fires before run
        # time, with a clear message pointing at both fix paths.
        if not self.density_mode and self.nfreq is None and self.df is None:
            raise ValueError(
                "cmd.BLS in r/q duration mode requires an explicit "
                "frequency grid: pass `nfreq=N` (number of "
                "frequencies) or `df=...` (frequency step).  The "
                "default `optimal` spacing is only valid when "
                "`density_mode=True` (which uses stellar density to "
                "set the expected transit duration)."
            )

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-BLS"]
        # Duration / transit-depth mode
        if self.density_mode:
            args += (["density"]
                     + _varexpr(self.stellar_density)
                     + _varexpr(self.min_exp_dur_frac)
                     + _varexpr(self.max_exp_dur_frac))
        elif self.qmin is not None or self.qmax is not None:
            args += ["q"] + _varexpr(self.qmin) + _varexpr(self.qmax)
        else:
            args += ["r"] + _varexpr(self.rmin) + _varexpr(self.rmax)
        args += _varexpr(self.minper) + _varexpr(self.maxper)
        # Frequency specification.  The valid-combination check fires
        # at construction time (in __init__), so any reachable state
        # here is one of: df set, nfreq set, or density_mode + optimal.
        if self.df is not None:
            args += ["df"] + _varexpr(self.df)
        elif self.nfreq is not None:
            args += ["nf"] + _varexpr(self.nfreq)
        else:
            # density_mode is True by the __init__ guard.
            args += ["optimal"] + _varexpr(self.subsample)
        args += _varexpr(self.nbins)
        args += [str(self.timezone), str(self.npeaks)]
        args += _outtoken(self.save_periodogram, outdir)
        args += _outtoken(self.save_model, outdir)
        args += ["1" if self.correct_lc else "0"]
        args += _bool("extraparams", self.extraparams)
        args += _bool("fittrap", self.fittrap)
        args += _bool("nobinnedrms", self.nobinnedrms)
        ph_spec = _norm_save(self.save_phcurve)
        if _should_emit(ph_spec):
            args += ["ophcurve", ph_spec.path or outdir,
                     str(self.ophcurve_phmin), str(self.ophcurve_phmax),
                     str(self.ophcurve_phstep)]
        jd_spec = _norm_save(self.save_jdcurve)
        if _should_emit(jd_spec):
            args += ["ojdcurve", jd_spec.path or outdir,
                     str(self.ojdcurve_jdstep)]
        if self.freq_grid is not None:
            args += [self.freq_grid]
        if self.adjust_qmin:
            args += ["adjust-qmin-by-mindt"]
            if self.reduce_nbins:
                args += ["reduce-nbins"]
        args += _bool("reportharmonics", self.reportharmonics)
        args += _flag("maskpoints", self.maskpoints)
        return args

    def _output_file_specs(self):
        return {
            "periodogram": (".bls", None),
            "model": (".bls.model", None),
            "phcurve": (".bls.phcurve", None),
            "jdcurve": (".bls.jdcurve", None),
        }


class BLSFixPer(VartoolsCommand):
    """BLS with a fixed (pre-known) period.

    Parameters
    ----------
    period : float or str
        Period source.
    rmin, rmax : float, optional
        Min/max fractional transit duration.  Used when ``qmin``/``qmax`` are
        not set.
    qmin, qmax : float, optional
        Min/max fractional transit duration as fraction of orbit.  When set,
        emits ``"q" qmin qmax`` instead of ``"r" rmin rmax``.
    """

    _vt_name = "BLSFixPer"

    def __init__(
        self,
        period,
        rmin: float = 0.01,
        rmax: float = 0.1,
        qmin: Optional[float] = None,
        qmax: Optional[float] = None,
        nbins: int = 200,
        timezone: float = 0,
        save_model: bool = False,
        correct_lc: bool = False,
        fittrap: bool = False,
        maskpoints: Optional[str] = None,
    ) -> None:
        self.period = period
        self.rmin = rmin
        self.rmax = rmax
        self.qmin = qmin
        self.qmax = qmax
        self.nbins = nbins
        self.timezone = timezone
        self.save_model = save_model
        self.correct_lc = correct_lc
        self.fittrap = fittrap
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-BLSFixPer"] + _period_spec(self.period)
        if self.qmin is not None or self.qmax is not None:
            args += ["q", str(self.qmin), str(self.qmax)]
        else:
            args += ["r", str(self.rmin), str(self.rmax)]
        args += [str(self.nbins), str(self.timezone)]
        args += _outtoken(self.save_model, outdir)
        args += ["1" if self.correct_lc else "0"]
        args += _bool("fittrap", self.fittrap)
        args += _flag("maskpoints", self.maskpoints)
        return args

    def _resolve_back_references(self, prev) -> None:
        from ._helpers import _resolve_period_backref
        self.period = _resolve_period_backref(prev, self.period)

    def _output_file_specs(self):
        return {"model": (".blsfixper.model", None)}


class BLSFixDurTc(VartoolsCommand):
    """BLS with fixed transit duration and epoch (Tc), searching for period.

    Parameters
    ----------
    duration : float or str
        Transit duration (same units as the period). Accepts a fixed float,
        ``"fixcolumn <colname>"``, or ``"list ["column" col]"`` forms.
    Tc : float or str
        Epoch of transit center. Accepts the same forms as *duration*.
    minper, maxper : float
        Period search range (days).
    nfreq : int
        Number of trial frequencies.
    timezone : float
        Timezone offset (add to JD to get local date; 0 for UTC/BJD).
    npeaks : int
        Number of peaks to report.
    fixdepth : float or str, optional
        Fix the transit depth to this value (or a ``"fixcolumn"``/``"list"``
        spec). When ``None`` the depth is optimised.
    qgress : float or str, optional
        Fractional ingress/egress duration (only used when *fixdepth* is set).
    save_periodogram : bool
        Write the BLS power spectrum to a file.
    save_model : bool
        Write the best-fit transit model to a file.
    correct_lc : bool
        Subtract the best-fit transit from the light curve.
    fittrap : bool
        Fit a trapezoidal transit shape instead of a box.
    save_phcurve : bool, str, or Output
        Write the phase-folded model light curve.
    ophcurve_phmin, ophcurve_phmax, ophcurve_phstep : float
        Phase range and step for phase curve output.
    save_jdcurve : bool, str, or Output
        Write the model light curve evaluated on a uniform JD grid.
    ojdcurve_jdstep : float
        Time step (days) for JD curve output.
    maskpoints : str, optional
        Name of a mask variable; points with mask > 0 are included.
    """

    _vt_name = "BLSFixDurTc"

    def __init__(
        self,
        duration,
        Tc,
        minper: float = 0.1,
        maxper: float = 100.0,
        nfreq: int = 10000,
        timezone: float = 0,
        npeaks: int = 1,
        fixdepth=None,
        qgress=None,
        save_periodogram: bool = False,
        save_model: bool = False,
        correct_lc: bool = False,
        fittrap: bool = False,
        save_phcurve=False,
        ophcurve_phmin: float = 0.0,
        ophcurve_phmax: float = 1.0,
        ophcurve_phstep: float = 0.005,
        save_jdcurve=False,
        ojdcurve_jdstep: float = 0.02,
        maskpoints: Optional[str] = None,
    ) -> None:
        self.duration = duration
        self.Tc = Tc
        self.minper = minper
        self.maxper = maxper
        self.nfreq = nfreq
        self.timezone = timezone
        self.npeaks = npeaks
        self.fixdepth = fixdepth
        self.qgress = qgress
        self.save_periodogram = save_periodogram
        self.save_model = save_model
        self.correct_lc = correct_lc
        self.fittrap = fittrap
        self.save_phcurve = save_phcurve
        self.ophcurve_phmin = ophcurve_phmin
        self.ophcurve_phmax = ophcurve_phmax
        self.ophcurve_phstep = ophcurve_phstep
        self.save_jdcurve = save_jdcurve
        self.ojdcurve_jdstep = ojdcurve_jdstep
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-BLSFixDurTc"]
        args += ["duration"] + _period_spec(self.duration)
        args += ["Tc"] + _period_spec(self.Tc)
        if self.fixdepth is not None:
            args += ["fixdepth"] + _period_spec(self.fixdepth)
            if self.qgress is not None:
                args += ["qgress"] + _period_spec(self.qgress)
        args += _varexpr(self.minper) + _varexpr(self.maxper) + _varexpr(self.nfreq)
        args += [str(self.timezone), str(self.npeaks)]
        args += _outtoken(self.save_periodogram, outdir)
        args += _outtoken(self.save_model, outdir)
        args += ["1" if self.correct_lc else "0"]
        args += _bool("fittrap", self.fittrap)
        ph_spec = _norm_save(self.save_phcurve)
        if _should_emit(ph_spec):
            args += ["ophcurve", ph_spec.path or outdir,
                     str(self.ophcurve_phmin), str(self.ophcurve_phmax),
                     str(self.ophcurve_phstep)]
        jd_spec = _norm_save(self.save_jdcurve)
        if _should_emit(jd_spec):
            args += ["ojdcurve", jd_spec.path or outdir,
                     str(self.ojdcurve_jdstep)]
        args += _flag("maskpoints", self.maskpoints)
        return args

    def _output_file_specs(self):
        return {
            "periodogram": (".blsfixdurtc", None),
            "model": (".blsfixdurtc.model", None),
            "phcurve": (".blsfixdurtc.phcurve", None),
            "jdcurve": (".blsfixdurtc.jdcurve", None),
        }


class BLSFixPerDurTc(VartoolsCommand):
    """BLS with fixed period, transit duration, and epoch (Tc).

    Computes BLS transit statistics for a fully specified box-transit signal
    without any period search.

    Parameters
    ----------
    period : float or str
        Transit period (days). Accepts a fixed float,
        ``"fixcolumn <colname>"``, or ``"list ["column" col]"`` forms.
    duration : float or str
        Transit duration (same units as the period). Accepts same forms.
    Tc : float or str
        Epoch of transit center. Accepts same forms.
    timezone : float
        Timezone offset (add to JD to get local date; 0 for UTC/BJD).
    fixdepth : float or str, optional
        Fix the transit depth to this value (or a ``"fixcolumn"``/``"list"``
        spec). When ``None`` the depth is optimised.
    qgress : float or str, optional
        Fractional ingress/egress duration (only used when *fixdepth* is set).
    save_model : bool
        Write the best-fit transit model to a file.
    correct_lc : bool
        Subtract the best-fit transit from the light curve.
    fittrap : bool
        Fit a trapezoidal transit shape instead of a box.
    save_phcurve : bool, str, or Output
        Write the phase-folded model light curve.
    ophcurve_phmin, ophcurve_phmax, ophcurve_phstep : float
        Phase range and step for phase curve output.
    save_jdcurve : bool, str, or Output
        Write the model light curve evaluated on a uniform JD grid.
    ojdcurve_jdstep : float
        Time step (days) for JD curve output.
    maskpoints : str, optional
        Name of a mask variable; points with mask > 0 are included.
    """

    _vt_name = "BLSFixPerDurTc"

    def __init__(
        self,
        period,
        duration,
        Tc,
        timezone: float = 0,
        fixdepth=None,
        qgress=None,
        save_model: bool = False,
        correct_lc: bool = False,
        fittrap: bool = False,
        save_phcurve=False,
        ophcurve_phmin: float = 0.0,
        ophcurve_phmax: float = 1.0,
        ophcurve_phstep: float = 0.005,
        save_jdcurve=False,
        ojdcurve_jdstep: float = 0.02,
        maskpoints: Optional[str] = None,
    ) -> None:
        self.period = period
        self.duration = duration
        self.Tc = Tc
        self.timezone = timezone
        self.fixdepth = fixdepth
        self.qgress = qgress
        self.save_model = save_model
        self.correct_lc = correct_lc
        self.fittrap = fittrap
        self.save_phcurve = save_phcurve
        self.ophcurve_phmin = ophcurve_phmin
        self.ophcurve_phmax = ophcurve_phmax
        self.ophcurve_phstep = ophcurve_phstep
        self.save_jdcurve = save_jdcurve
        self.ojdcurve_jdstep = ojdcurve_jdstep
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-BLSFixPerDurTc"]
        args += ["period"] + _period_spec(self.period)
        args += ["duration"] + _period_spec(self.duration)
        args += ["Tc"] + _period_spec(self.Tc)
        if self.fixdepth is not None:
            args += ["fixdepth"] + _period_spec(self.fixdepth)
            if self.qgress is not None:
                args += ["qgress"] + _period_spec(self.qgress)
        args += [str(self.timezone)]
        args += _outtoken(self.save_model, outdir)
        args += ["1" if self.correct_lc else "0"]
        args += _bool("fittrap", self.fittrap)
        ph_spec = _norm_save(self.save_phcurve)
        if _should_emit(ph_spec):
            args += ["ophcurve", ph_spec.path or outdir,
                     str(self.ophcurve_phmin), str(self.ophcurve_phmax),
                     str(self.ophcurve_phstep)]
        jd_spec = _norm_save(self.save_jdcurve)
        if _should_emit(jd_spec):
            args += ["ojdcurve", jd_spec.path or outdir,
                     str(self.ojdcurve_jdstep)]
        args += _flag("maskpoints", self.maskpoints)
        return args

    def _output_file_specs(self):
        return {
            "model": (".blsfixperdurtc.model", None),
            "phcurve": (".blsfixperdurtc.phcurve", None),
            "jdcurve": (".blsfixperdurtc.jdcurve", None),
        }


class autocorrelation(VartoolsCommand):
    """Autocorrelation function of the light curve.

    The CLI always writes the autocorrelation function to a file — there is no
    mode that suppresses file output entirely.

    Parameters
    ----------
    start : float
        Start of the lag range.
    stop : float
        End of the lag range.
    step : float
        Lag step size.
    save_result : bool, str, or Output
        Controls capture of the output file into Python.

        - ``True`` (default) — write to a temp dir and capture into
          ``result.files["autocorrelation_result_N"]``.
        - ``False`` — write to a temp dir but do **not** capture into Python
          (the file is still written because the CLI always does so).
        - A directory path string — write to that directory, no capture.
        - ``Output(path, capture=True)`` — write to that directory and capture.
    maskpoints : str, optional
    """

    _vt_name = "autocorrelation"
    # The CLI always writes the output file regardless of save_result.
    _mandatory_output = True

    def __init__(
        self,
        start: float,
        stop: float,
        step: float,
        save_result=True,
        maskpoints: Optional[str] = None,
    ) -> None:
        self.start = start
        self.stop = stop
        self.step = step
        self.save_result = save_result
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        spec = _norm_save(self.save_result)
        actual_outdir = spec.path if spec.path is not None else outdir
        args = ["-autocorrelation"] + _varexpr(self.start) + _varexpr(self.stop) + _varexpr(self.step) + [actual_outdir]
        args += _flag("maskpoints", self.maskpoints)
        return args

    def _output_file_specs(self):
        return {"result": (".autocorr", None)}


class dftclean(VartoolsCommand):
    """CLEAN deconvolution periodogram.

    Parameters
    ----------
    nbeam : int
        Number of frequencies in the beam (dirty spectrum).
    maxfreq : float, optional
        Maximum frequency.  If None, defaults to ``"auto"`` (vartools decides).
    save_dspec : bool
        Write the dirty spectrum.
    finddirtypeaks : int, optional
        Find the top N peaks in the dirty power spectrum.
    finddirtypeaks_clip : float, optional
        Sigma-clipping value for dirty-peak SNR (requires ``finddirtypeaks``).
    finddirtypeaks_clipiter : int, optional
        0 = no iterative clipping, 1 = iterative clipping (requires
        ``finddirtypeaks_clip``).
    save_wfunc : bool
        Write the window function.
    save_cspec : bool
        Write the CLEAN spectrum.
    gain : float, optional
        CLEAN gain factor.
    SNlimit : float, optional
        Stop CLEANing when the peak falls below this S/N.
    outcbeam : bool, str, or Output, optional
        Write the clean beam (requires the ``clean`` section to be active).
    npeaks : int, optional
        Number of peaks to find in the CLEAN spectrum.
    useampspec : bool
        Use the amplitude spectrum rather than the power spectrum for SNR.
    verboseout : bool
        Output average and standard deviation of spectrum before/after clipping.
    maskpoints : str, optional
    """

    _vt_name = "dftclean"

    def __init__(
        self,
        nbeam: int,
        maxfreq: Optional[float] = None,
        save_dspec: bool = False,
        finddirtypeaks: Optional[int] = None,
        finddirtypeaks_clip: Optional[float] = None,
        finddirtypeaks_clipiter: Optional[int] = None,
        save_wfunc: bool = False,
        save_cspec: bool = False,
        gain: float = 0.1,
        SNlimit: float = 3.0,
        outcbeam=False,
        npeaks: Optional[int] = None,
        useampspec: bool = False,
        verboseout: bool = False,
        maskpoints: Optional[str] = None,
    ) -> None:
        self.nbeam = nbeam
        self.maxfreq = maxfreq
        self.save_dspec = save_dspec
        self.finddirtypeaks = finddirtypeaks
        self.finddirtypeaks_clip = finddirtypeaks_clip
        self.finddirtypeaks_clipiter = finddirtypeaks_clipiter
        self.save_wfunc = save_wfunc
        self.save_cspec = save_cspec
        self.gain = gain
        self.SNlimit = SNlimit
        self.outcbeam = outcbeam
        # Alias so the capture pipeline (which keys on ``save_<logical_name>``)
        # picks up the request for the CLEAN beam file.
        self.save_cbeam = outcbeam
        self.npeaks = npeaks
        self.useampspec = useampspec
        self.verboseout = verboseout
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-dftclean"] + _varexpr(self.nbeam)
        if self.maxfreq is not None:
            args += ["maxfreq"] + _varexpr(self.maxfreq)
        dspec_spec = _norm_save(self.save_dspec)
        if _should_emit(dspec_spec):
            args += ["outdspec", dspec_spec.path or outdir]
        if self.finddirtypeaks is not None:
            args += ["finddirtypeaks", str(self.finddirtypeaks)]
            if self.finddirtypeaks_clip is not None:
                args += ["clip"] + _varexpr(self.finddirtypeaks_clip)
                if self.finddirtypeaks_clipiter is not None:
                    args += [str(self.finddirtypeaks_clipiter)]
        wfunc_spec = _norm_save(self.save_wfunc)
        if _should_emit(wfunc_spec):
            args += ["outwfunc", wfunc_spec.path or outdir]
        cspec_spec = _norm_save(self.save_cspec)
        cb_spec = _norm_save(self.outcbeam)
        if _should_emit(cspec_spec) or self.npeaks is not None or _should_emit(cb_spec):
            args += ["clean"] + _varexpr(self.gain) + _varexpr(self.SNlimit)
            if _should_emit(cb_spec):
                args += ["outcbeam", cb_spec.path or outdir]
            if _should_emit(cspec_spec):
                args += ["outcspec", cspec_spec.path or outdir]
            if self.npeaks is not None:
                args += ["findcleanpeaks", str(self.npeaks)]
        args += _bool("useampspec", self.useampspec)
        args += _bool("verboseout", self.verboseout)
        args += _flag("maskpoints", self.maskpoints)
        return args

    def _output_file_specs(self):
        return {
            "dspec": (".dftclean.dspec", None),
            "wfunc": (".dftclean.wfunc", None),
            "cspec": (".dftclean.cspec", None),
            "cbeam": (".dftclean.cbeam", None),
        }


class wwz(VartoolsCommand):
    """Weighted Wavelet Z-transform time-frequency analysis.

    Parameters
    ----------
    maxfreq : float or ``"auto"``
        Maximum frequency in cycles per day.  ``"auto"`` (default) sets
        it to ``1/(2*delmin)`` where ``delmin`` is the minimum spacing
        between consecutive light-curve points.
    freqsamp : float
        Frequency sampling step factor: the actual frequency step is
        ``freqsamp/T`` where ``T`` is the light-curve time baseline.
        Default ``0.25`` (Foster 1996 convention).  vartools does not
        accept ``"auto"`` for this parameter.
    tau0, tau1, dtau : float or ``"auto"``
        Start, end, and step of the time grid.  ``"auto"`` for
        ``tau0`` / ``tau1`` uses the LC's first/last time; ``"auto"``
        for ``dtau`` uses ``delmin``.
    c : float, optional
        Decay constant of the abbreviated Morlet wavelet (default
        0.0125, i.e. the Foster ``1/(8*pi^2)`` recommendation).
    save_transform : bool
        Write the full WWZ time-frequency map.
    transform_format : str, optional
        Output format for the full transform: ``"fits"`` or ``"pm3d"``.
    transform_name : str, optional
        Format string for the full transform output filename.
    save_maxtransform : bool
        Write the WWZ max-power transform.
    maxtransform_name : str, optional
        Format string for the max-transform output filename.
    maskpoints : str, optional
    """

    _vt_name = "wwz"

    def __init__(
        self,
        maxfreq="auto",
        freqsamp: float = 0.25,
        tau0="auto",
        tau1="auto",
        dtau="auto",
        c: float = 0.0125,
        save_transform: bool = False,
        transform_format: Optional[str] = None,
        transform_name: Optional[str] = None,
        save_maxtransform: bool = False,
        maxtransform_name: Optional[str] = None,
        maskpoints: Optional[str] = None,
    ) -> None:
        self.maxfreq = maxfreq
        self.freqsamp = freqsamp
        self.tau0 = tau0
        self.tau1 = tau1
        self.dtau = dtau
        self.c = c
        self.save_transform = save_transform
        self.transform_format = transform_format
        self.transform_name = transform_name
        self.save_maxtransform = save_maxtransform
        self.maxtransform_name = maxtransform_name
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-wwz"]
        # maxfreq / tau0 / tau1 / dtau accept the bare keyword ``auto`` as
        # a literal CLI option (route them through the auto-aware helper
        # rather than the plain _varexpr).  freqsamp does NOT — vartools
        # rejects ``freqsamp auto`` with "must be > 0".
        args += ["maxfreq"] + _auto_or_varexpr(self.maxfreq)
        args += ["freqsamp"] + _varexpr(self.freqsamp)
        args += ["tau0"] + _auto_or_varexpr(self.tau0)
        args += ["tau1"] + _auto_or_varexpr(self.tau1)
        args += ["dtau"] + _auto_or_varexpr(self.dtau)
        args += ["c"] + _varexpr(self.c)
        tr_spec = _norm_save(self.save_transform)
        if _should_emit(tr_spec):
            args += ["outfulltransform", tr_spec.path or outdir]
            if self.transform_format is not None:
                args += [self.transform_format]
            if self.transform_name is not None:
                args += ["format", self.transform_name]
        mt_spec = _norm_save(self.save_maxtransform)
        if _should_emit(mt_spec):
            args += ["outmaxtransform", mt_spec.path or outdir]
            if self.maxtransform_name is not None:
                args += ["format", self.maxtransform_name]
        args += _flag("maskpoints", self.maskpoints)
        return args

    def _output_file_specs(self):
        return {
            "transform": (".wwz", None),
            "maxtransform": (".mwwz", None),
        }


class GetLSAmpThresh(VartoolsCommand):
    """Get LS amplitude threshold for signal injection/recovery.

    Parameters
    ----------
    period : str
        Period source: ``"ls"`` or ``"list"`` column.
    minp : float
        Minimum period.
    thresh : float
        Detection threshold.
    mode : str
        ``"harm"`` (default) or ``"file"``.  When ``"harm"``, the signal is
        modelled as a Fourier series with ``nharm`` and ``nsubharm`` terms.
        When ``"file"``, the signal is read from ``listfile``.
    nharm : int
        Number of harmonics (used when ``mode="harm"``).
    nsubharm : int
        Number of sub-harmonics (used when ``mode="harm"``).
    listfile : str, optional
        Path to the signal list file (used when ``mode="file"``).
    noGLS : bool
        Use classical LS (not generalised).
    """

    _vt_name = "GetLSAmpThresh"

    def __init__(
        self,
        period: str = "ls",
        minp: float = 0.1,
        thresh: float = 10.0,
        mode: str = "harm",
        nharm: int = 1,
        nsubharm: int = 0,
        listfile: Optional[str] = None,
        noGLS: bool = False,
    ) -> None:
        self.period = period
        self.minp = minp
        self.thresh = thresh
        self.mode = mode
        self.nharm = nharm
        self.nsubharm = nsubharm
        self.listfile = listfile
        self.noGLS = noGLS

    def _to_cli_args(self) -> List[str]:
        args = ["-GetLSAmpThresh"] + _period_spec(self.period)
        args += [str(self.minp), str(self.thresh)]
        if self.mode == "file":
            if self.listfile is None:
                raise ValueError(
                    "GetLSAmpThresh(mode='file') requires the listfile "
                    "kwarg to point at a vartools signal-list file."
                )
            args += ["file", str(self.listfile)]
        elif self.mode == "harm":
            args += ["harm", str(self.nharm), str(self.nsubharm)]
        else:
            raise ValueError(
                f"GetLSAmpThresh: mode must be 'harm' or 'file'; "
                f"got {self.mode!r}"
            )
        args += _bool("noGLS", self.noGLS)
        return args

    def _resolve_back_references(self, prev) -> None:
        # The -GetLSAmpThresh CLI only accepts "ls" or "list" as its period
        # source.  Across a chain boundary the prior -LS command does not
        # survive, so the "ls" keyword must be resolved.  We cannot convert
        # to a bare number (the CLI rejects that), so we raise a clear error
        # and point the user at the direct form.
        if isinstance(self.period, str) and self.period.strip() == "ls":
            from ._helpers import _most_recent_lookup
            stats = _most_recent_lookup(prev, ["LS"])
            if stats is None:
                raise LookupError(
                    "GetLSAmpThresh back-reference 'ls' has no prior -LS "
                    "command in this chain"
                )
            raise NotImplementedError(
                "GetLSAmpThresh with period='ls' is not supported across "
                "chain boundaries (the CLI only accepts 'ls'/'list' period "
                "sources).  Call it within a single Pipeline alongside the "
                "-LS command."
            )


class Phase(VartoolsCommand):
    """Phase-fold the light curve on a period.

    Parameters
    ----------
    period : float or str
        Period to fold on.  Can be a number, or ``"ls"``, ``"aov"``, ``"bls"``
        to use the period found by a prior command, or a string like
        ``"fixcolumn colname"`` or ``"list"``.
    T0 : float or str, optional
        Reference epoch.  Number, or ``"bls phaseTc"``, or ``"fix T0"``, etc.
    phasevar : str, optional
        Name of the output phase variable.
    startphase : float, optional
        Starting phase (default 0).
    """

    _vt_name = "Phase"

    def __init__(
        self,
        period="ls",
        T0=None,
        phasevar: Optional[str] = None,
        startphase: Optional[float] = None,
    ) -> None:
        self.period = period
        self.T0 = T0
        self.phasevar = phasevar
        self.startphase = startphase

    def _to_cli_args(self) -> List[str]:
        args = ["-Phase"] + _period_spec(self.period)
        if self.T0 is not None:
            args += ["T0"] + _pval(self.T0, "fix")
        args += _flag("phasevar", self.phasevar)
        if self.startphase is not None:
            args += ["startphase", str(self.startphase)]
        return args

    def _resolve_back_references(self, prev) -> None:
        from ._helpers import _resolve_period_backref, _most_recent_lookup

        # period: resolve "ls"/"aov"/"bls"/"blsfixper"/"injectharm"/"fixcolumn"
        self.period = _resolve_period_backref(prev, self.period)

        # T0 accepts "bls <phaseTc>" — pull Tc from the prior BLS and add the
        # phase offset.  Example: "bls 0.5" → bls_Tc - 0.5*bls_Period.
        if isinstance(self.T0, str):
            s = self.T0.strip()
            if s.startswith("bls"):
                parts = s.split()
                phase_off = float(parts[1]) if len(parts) > 1 else 0.0
                stats = _most_recent_lookup(prev, ["BLS"])
                if stats is None:
                    raise LookupError(
                        "Back-reference T0='bls ...' has no prior -BLS "
                        "command in this chain"
                    )
                from ._helpers import _coerce_to_numeric
                period = _coerce_to_numeric(stats.Period_1)
                tc = _coerce_to_numeric(stats.Tc_1)
                # Scalar vs PerLC arithmetic
                from pyvartools.perlc import PerLC
                if isinstance(tc, PerLC) or isinstance(period, PerLC):
                    tc_vals = tc if isinstance(tc, PerLC) else PerLC(
                        [float(tc)] * len(period))
                    per_vals = period if isinstance(period, PerLC) else PerLC(
                        [float(period)] * len(tc_vals))
                    self.T0 = PerLC([
                        float(t) - phase_off * float(p)
                        for t, p in zip(tc_vals, per_vals)
                    ])
                else:
                    self.T0 = float(tc) - phase_off * float(period)
