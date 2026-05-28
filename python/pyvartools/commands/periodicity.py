"""Periodicity-search command wrappers."""

from __future__ import annotations
from typing import List, Optional, Union

from pyvartools._command import VartoolsCommand
from ._helpers import (_auto_or_varexpr, _bool, _flag, _fixperiodsnr_tokens,
                       _norm_save, _outtoken, _period_spec, _pval,
                       _should_emit, _varexpr)


class LS(VartoolsCommand):
    """Generalized Lomb-Scargle (GLS) periodogram.

    Search frequencies from ``1/maxp`` to ``1/minp`` with step
    ``Δf = subsample / T`` (``T`` = time baseline).  The reported
    statistic is ``LS = (χ0² − χ(f)²) / χ0²``; with ``noGLS=True`` the
    wrapper computes the standard un-normalised Lomb-Scargle power.

    Parameters
    ----------
    minp, maxp : float, str, numpy array, PerLC, or pd.Series
        Period search range (same units as the time column, typically
        days).  Numeric forms are validated at construction time:
        ``minp > 0``, ``maxp > 0``, and ``minp < maxp`` — a clear
        ``ValueError`` is raised otherwise.  Forms:

        - A number — passed directly as a fixed value.
        - A bare identifier string (e.g. ``"minperiod"``) — vartools
          reads from a named per-LC variable (``var minperiod``).
        - Any other string (e.g. ``"tspan/100"``) — evaluated as a math
          expression (``expr tspan/100``).
        - A numpy array, ``PerLC``, or ``pd.Series`` for per-LC batch
          values.
    subsample : float, str, numpy array, PerLC, or pd.Series
        Frequency step as a fraction of 1/T.  Smaller values give finer
        resolution; typical 1e-3.  Same value forms as *minp* / *maxp*.
    npeaks : int
        Number of highest peaks to report.  Default 5.
    save_periodogram : bool, str, or Output
        Auxiliary file output.  ``True`` captures as
        ``result.files["LS_periodogram_N"]``; a path string writes to that
        directory without capturing; ``Output(path, capture=True)`` does
        both.
    noGLS : bool
        Use the classical (non-generalised) Lomb-Scargle statistic
        instead of GLS.
    whiten : bool
        After each peak, whiten the light curve at that period before
        searching for the next.  The peak SNR is computed on the
        whitened periodogram.
    clip : float, optional
        Sigma-clipping factor for the mean / RMS used in the SNR noise
        estimate (default: iterative 5σ).
    clipiter : int, optional
        Number of clipping iterations.  ``1`` enables iterative clipping
        (the default when *clip* is set).
    bootstrap : int, optional
        Number of bootstrap resamples for empirical false-alarm
        probability estimation.
    maskpoints : str, optional
        Name of a mask variable; points where the variable is ``≤ 0``
        are excluded.
    fixperiod_snr : float, int, str, or None, optional
        Evaluate the periodogram at a known period and report its
        significance.  Forms:

        - A number (e.g. ``1.234``) — evaluate at that fixed period.
        - ``"ls"`` / ``"aov"`` / ``"injectharm"`` — back-reference to the
          best period from a prior LS, AOV, or injection-recovery
          command in the same chain.
        - ``"fixcolumn COLNAME"`` — read the period from a named per-star
          column.
        - ``"list"`` / ``"list column 2"`` — read the period from a list-
          file column.

        When set, four extra output columns are appended:
        ``LS_PeriodFix_N``, ``Log10_LS_Prob_PeriodFix_N``,
        ``LS_Periodogram_Value_PeriodFix_N``, ``LS_SNR_PeriodFix_N``.

    See Also
    --------
    CLI command: ``-LS``.
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
    """Phase-binned Analysis-of-Variance (AoV) periodogram.

    For each trial frequency the light curve is phase-folded and binned;
    the AoV statistic ``θ_aov`` measures the variance explained by the
    phase bins relative to the total variance.  Peaks at the correct
    period drive ``θ_aov`` high.  AoV is preferable to LS for strictly
    periodic but non-sinusoidal signals (eclipsing binaries, pulsators).

    Parameters
    ----------
    minp, maxp : float, str, numpy array, PerLC, or pd.Series
        Period search range (same units as the time column).  Forms:

        - A number — passed directly as a fixed value.
        - A bare identifier string (e.g. ``"minperiod"``) — vartools
          reads from a named per-LC variable (``var minperiod``).
        - Any other string (e.g. ``"tspan/100"``) — evaluated as a math
          expression (``expr tspan/100``).
        - A numpy array, ``PerLC``, or ``pd.Series`` for per-LC batch
          values.
    subsample : float, str, numpy array, PerLC, or pd.Series
        Coarse-grid frequency step as a fraction of 1/T (``T`` = time
        baseline).  Same value forms as *minp* / *maxp*.
    finetune : float, str, numpy array, PerLC, or pd.Series
        Fine-tune frequency-step fraction applied near peaks.  Same
        value forms as *minp* / *maxp*.
    npeaks : int
        Number of peaks to report.  Default 5.
    nbin : int, str, numpy array, PerLC, or pd.Series, optional
        Number of phase bins.  ``None`` uses the vartools default of 8.
        Accepts the same value forms as *minp* / *maxp*.
    save_periodogram : bool, str, or Output
        Auxiliary file output.  ``True`` captures as
        ``result.files["aov_periodogram_N"]``; a path string writes to
        that directory without capturing; ``Output(path, capture=True)``
        does both.
    whiten : bool
        After each peak, whiten the light curve at that period before
        searching for the next.
    clip : float, optional
        Sigma-clipping factor for the mean / RMS used in the SNR noise
        estimate (default: iterative 5σ).
    clipiter : int, optional
        Number of clipping iterations.  ``1`` enables iterative
        clipping (the default when *clip* is set).
    uselog : bool
        Use ``−ln(θ_aov)`` for the SNR statistic; also outputs the mean
        and RMS of ``−ln(θ_aov)``.
    maskpoints : str, optional
        Name of a mask variable; points where the variable is ``≤ 0``
        are excluded.
    fixperiod_snr : float, int, str, or None, optional
        Evaluate the AoV periodogram at a known period and report its
        significance.  Forms:

        - A number (e.g. ``1.234``) — evaluate at that fixed period.
        - ``"ls"`` / ``"aov"`` / ``"injectharm"`` — back-reference to the
          best period from a prior LS, AOV, or injection-recovery
          command in the same chain.
        - ``"fixcolumn COLNAME"`` — read the period from a named per-star
          column.
        - ``"list"`` / ``"list column 2"`` — read the period from a list-
          file column.

        When set, four extra output columns are appended:
        ``PeriodFix_N``, ``AOV_PeriodFix_N``, ``AOV_SNR_PeriodFix_N``,
        ``AOV_NEG_LN_FAP_PeriodFix_N``.  When ``uselog=True`` only
        ``PeriodFix_N`` and ``AOV_LOGSNR_PeriodFix_N`` are emitted.

    See Also
    --------
    CLI command: ``-aov``.
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
    """Multi-harmonic Analysis-of-Variance (AoV-Harm) periodogram.

    Replaces the phase-binned model of :class:`aov` with a Fourier
    series of ``nharm`` harmonics fit at each trial period.  Preferable
    for smoothly-varying non-sinusoidal signals (RR Lyrae, Cepheids,
    W UMa systems).

    Parameters
    ----------
    nharm : int, str, numpy array, PerLC, or pd.Series
        Number of harmonics in the model.  Set ``nharm <= 0`` to enable
        automatic selection — vartools then chooses ``nharm`` per peak to
        minimise the false-alarm probability (with an overfitting
        penalty).  Accepts variable names, expressions, and per-LC
        batch values.
    minp, maxp : float, str, numpy array, PerLC, or pd.Series
        Period search range (same units as the time column).  Forms:

        - A number — passed directly as a fixed value.
        - A bare identifier string (e.g. ``"minperiod"``) — vartools
          reads from a named per-LC variable (``var minperiod``).
        - Any other string (e.g. ``"tspan/100"``) — evaluated as a math
          expression (``expr tspan/100``).
        - A numpy array, ``PerLC``, or ``pd.Series`` for per-LC batch
          values.
    subsample : float, str, numpy array, PerLC, or pd.Series
        Coarse-grid frequency step as a fraction of 1/T.  Same value
        forms as *minp* / *maxp*.
    finetune : float, str, numpy array, PerLC, or pd.Series
        Fine-tune frequency-step fraction applied near peaks.  Same
        value forms as *minp* / *maxp*.
    npeaks : int
        Number of peaks to report.  Default 5.
    save_periodogram : bool, str, or Output
        Auxiliary file output.  ``True`` captures as
        ``result.files["aov_harm_periodogram_N"]``; a path string writes
        to that directory without capturing; ``Output(path, capture=True)``
        does both.
    whiten : bool
        After each peak, whiten the light curve at that period before
        searching for the next.
    clip : float, optional
        Sigma-clipping factor for the mean / RMS used in the SNR noise
        estimate (default: iterative 5σ).
    clipiter : int, optional
        Number of clipping iterations.  ``1`` enables iterative
        clipping (the default when *clip* is set).
    maskpoints : str, optional
        Name of a mask variable; points where the variable is ``≤ 0``
        are excluded.
    fixperiod_snr : float, int, str, or None, optional
        Evaluate the multi-harmonic AoV periodogram at a known period
        and report its significance.  Forms:

        - A number (e.g. ``1.234``) — evaluate at that fixed period.
        - ``"ls"`` / ``"aov"`` / ``"injectharm"`` — back-reference to the
          best period from a prior LS, AOV, or injection-recovery
          command in the same chain.
        - ``"fixcolumn COLNAME"`` — read the period from a named per-star
          column.
        - ``"list"`` / ``"list column 2"`` — read the period from a list-
          file column.

        When set, four extra columns are appended:
        ``PeriodFix_N``, ``AOV_HARM_PeriodFix_N``,
        ``AOV_HARM_SNR_PeriodFix_N``,
        ``AOV_HARM_NEG_LN_FAP_PeriodFix_N`` (the last only when
        ``nharm > 0``).

    See Also
    --------
    CLI command: ``-aov_harm``.
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


_PDM_VARIANTS = ("step", "linterp", "multicover", "tophat", "gauss")
_PDM_BINNED_VARIANTS = ("step", "linterp", "multicover")
_PDM_BINLESS_VARIANTS = ("tophat", "gauss")


class PDM(VartoolsCommand):
    """Phase Dispersion Minimization (PDM) periodogram.

    Stellingwerf-1978 / Schwarzenberg-Czerny-1997 phase-folding statistic.
    Five algorithm variants (selected via ``variant``): three binned and two
    binless.

    Parameters
    ----------
    variant : {"step", "linterp", "multicover", "tophat", "gauss"}
        ``step``: classic Stellingwerf fixed-bin theta.
        ``linterp``: cuvarbase-style linear interpolation between bin means.
        ``multicover``: Nc shifted bin sets averaged into a single theta.
        ``tophat``: binless, hard phase-window kernel of half-width ``dphi``.
        ``gauss``: binless, Gaussian phase kernel of sigma ``dphi``.
    minp, maxp : float, str, numpy array, PerLC, or pd.Series
        Period search range (same units as the time column).  Forms:

        - A number — passed directly as a fixed value.
        - A bare identifier string (e.g. ``"minperiod"``) — vartools
          reads from a named per-LC variable (``var minperiod``).
        - Any other string (e.g. ``"tspan/100"``) — evaluated as a math
          expression (``expr tspan/100``).
        - A numpy array, ``PerLC``, or ``pd.Series`` for per-LC batch
          values.
    subsample, finetune : float, str, numpy array, PerLC, or pd.Series
        Coarse-grid step and fine-tune resolution as fractions of 1/T
        (``T`` = time baseline).  Same value forms as *minp* / *maxp*.
    npeaks : int
        Number of peaks to report (default 5).
    nbin : int or str, optional
        Phase bin count for binned variants (step/linterp) or bins-per-cover
        for multicover.  Rejected with binless variants.  vartools defaults
        to 8 when not set.
    nc : int or str, optional
        Number of phase-shifted covers for the ``multicover`` variant.
        Rejected for other variants.  vartools defaults to 2.
    dphi : float or str, optional
        Half-width (``tophat``) or sigma (``gauss``) of the phase kernel for
        binless variants.  Rejected for binned variants.  vartools defaults
        to 0.05 (cuvarbase convention).
    save_periodogram : bool, str, or Output
        Auxiliary file output.  ``True`` captures as
        ``result.files["pdm_periodogram_N"]``; a path string writes to
        that directory without capturing; ``Output(path, capture=True)``
        does both.  ``False`` (default) suppresses.  When ``whiten=True``
        the file holds one column per whitening cycle.
    clip : float, optional
        Sigma-clip factor for the SNR noise estimate.  vartools defaults
        to 5.
    clipiter : int, optional
        0 or 1: disable / enable iterative clipping.  vartools defaults to 1.
    noerr : bool
        Force uniform weights instead of 1/sigma^2.
    whiten : bool
        Iterative pre-whitening between peaks (subtract step-bin phase model
        at each peak before searching for the next).
    fixperiod_snr : float, int, or str, optional
        Additional theta/SNR/FAP at a specified period.  Forms recognised:
        ``"aov"`` / ``"ls"`` / ``"pdm"`` / ``"injectharm"`` (back-reference to
        the most recent prior command of that type in the same pipeline),
        ``"fixcolumn <name>"``, ``"list [column N]"``, or a numeric value
        (treated as a literal period).
    bootstrap : int, optional
        Enable empirical-CDF FAP via ``Nboot`` shuffled-LC trials.  Replaces
        the analytic Schwarzenberg-Czerny Beta FAP when set.
    maskpoints : str, optional
        Name of an LC vector; points with maskvar > VARTOOLS_MASK_TINY are
        included, others excluded.

    See Also
    --------
    CLI command: ``-PDM``.
    Citations: Stellingwerf 1978 (ApJ 224, 953);
    Schwarzenberg-Czerny 1997 (ApJ 489, 941);
    Zalian, Chadid & Stellingwerf 2014 (MNRAS 440, 68).
    The ``linterp`` variant follows the cuvarbase reference implementation
    (https://github.com/johnh2o2/cuvarbase) authored by Attila Bodi.
    """

    _vt_name = "PDM"

    def __init__(
        self,
        variant: str,
        minp: Union[float, str],
        maxp: Union[float, str],
        subsample: Union[float, str],
        finetune: Union[float, str],
        *,
        npeaks: int = 5,
        nbin: Optional[Union[int, str]] = None,
        nc: Optional[Union[int, str]] = None,
        dphi: Optional[Union[float, str]] = None,
        save_periodogram=False,
        clip: Optional[float] = None,
        clipiter: Optional[int] = None,
        noerr: bool = False,
        whiten: bool = False,
        fixperiod_snr: Union[float, int, str, None] = None,
        bootstrap: Optional[int] = None,
        maskpoints: Optional[str] = None,
    ) -> None:
        # Constructor-time validation -- mirrors the strict-parser behaviour
        # in vartools so misuse is caught early at the pyvartools layer.
        if variant not in _PDM_VARIANTS:
            raise ValueError(
                f"PDM variant must be one of {_PDM_VARIANTS!r}, got {variant!r}"
            )
        if variant in _PDM_BINLESS_VARIANTS and nbin is not None:
            raise ValueError(
                "PDM: 'nbin' is not used by binless variants (tophat/gauss); "
                "use 'dphi' instead"
            )
        if variant != "multicover" and nc is not None:
            raise ValueError(
                "PDM: 'nc' is only valid with the multicover variant"
            )
        if variant in _PDM_BINNED_VARIANTS and dphi is not None:
            raise ValueError(
                "PDM: 'dphi' is only valid with the tophat or gauss variants"
            )
        if bootstrap is not None:
            try:
                _nboot = int(bootstrap)
            except (TypeError, ValueError):
                raise ValueError(
                    f"PDM: 'bootstrap' must be a positive integer, got {bootstrap!r}"
                )
            if _nboot < 1:
                raise ValueError(
                    f"PDM: 'bootstrap' must be >= 1, got {_nboot}"
                )
            bootstrap = _nboot

        self.variant = variant
        self.minp = minp
        self.maxp = maxp
        self.subsample = subsample
        self.finetune = finetune
        self.npeaks = npeaks
        self.nbin = nbin
        self.nc = nc
        self.dphi = dphi
        self.save_periodogram = save_periodogram
        self.clip = clip
        self.clipiter = clipiter
        self.noerr = noerr
        self.whiten = whiten
        self.fixperiod_snr = fixperiod_snr
        self.bootstrap = bootstrap
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-PDM", self.variant]
        # Variant-specific keyword params (order: Nbin, Nc, dphi).
        if self.nbin is not None:
            args += ["Nbin"] + _varexpr(self.nbin)
        if self.nc is not None:
            args += ["Nc"] + _varexpr(self.nc)
        if self.dphi is not None:
            args += ["dphi"] + _varexpr(self.dphi)
        # Required positionals.
        args += (_varexpr(self.minp) + _varexpr(self.maxp)
                 + _varexpr(self.subsample) + _varexpr(self.finetune)
                 + [str(self.npeaks)])
        args += _outtoken(self.save_periodogram, outdir)
        # Trailing keywords in the canonical (strict-parser) order:
        # clip / noerr / whiten / fixperiodSNR / bootstrap / maskpoints.
        if self.clip is not None:
            args += ["clip", str(self.clip), str(self.clipiter if self.clipiter is not None else 1)]
        args += _bool("noerr", self.noerr)
        args += _bool("whiten", self.whiten)
        args += _fixperiodsnr_tokens(self.fixperiod_snr)
        if self.bootstrap is not None:
            args += ["bootstrap", str(int(self.bootstrap))]
        args += _flag("maskpoints", self.maskpoints)
        return args

    def _resolve_back_references(self, prev) -> None:
        from ._helpers import _resolve_period_backref
        self.fixperiod_snr = _resolve_period_backref(prev, self.fixperiod_snr)

    def _output_file_specs(self):
        return {"periodogram": (".pdm", None)}


_FTP_TEMPLATE_SOURCES = ("file", "fitlc", "inline", "filelist")
_FTP_METHODS = ("auto", "brute", "poly", "verify")
_FTP_SUMS = ("auto", "direct", "nfft")


class FTP(VartoolsCommand):
    """Fast Template Periodogram (FTP).

    Hoffman et al. 2021 (arXiv:2101.12348) non-linear extension of the
    Generalized Lomb-Scargle periodogram that fits a known periodic template
    shape M(phi) at each trial period.  Four template-source modes select
    where the truncated Fourier series ``M(phi) = sum_{n=1..H} c_n cos(n phi)
    + s_n sin(n phi)`` comes from.

    Parameters
    ----------
    template_source : {"file", "fitlc", "inline", "filelist"}
        Selects how the template is sourced.  Each mode requires a different
        set of mode-specific keyword arguments; mixing is rejected at
        construction time:

        - ``"file"``:     requires ``template_file`` (path to a two-column
                          c_n s_n text file; H inferred from row count).
        - ``"fitlc"``:    requires ``lc_path``, ``lc_format``, ``t_col``,
                          ``mag_col``, ``err_col``, ``nharm``, ``period``.
                          Builds the template by fitting a Fourier series of
                          order ``nharm + 1`` to the LC at ``lc_path`` at
                          fixed ``period``.
        - ``"inline"``:   requires ``cn`` and ``sn`` lists of equal length
                          ``nharm + 1`` (``nharm`` inferred from ``len(cn)``).
                          Each entry can be a number or a bare identifier /
                          expression string (var/expr semantics).
        - ``"filelist"``: optionally ``filelist_column`` (1-indexed integer);
                          the template path is read from that column of the
                          ``-l`` input list.
    minp, maxp : float, str, numpy array, PerLC, or pd.Series
        Period search range (same units as the time column).  Forms:

        - A number — passed directly as a fixed value.
        - A bare identifier string (e.g. ``"minperiod"``) — vartools
          reads from a named per-LC variable (``var minperiod``).
        - Any other string (e.g. ``"tspan/100"``) — evaluated as a math
          expression (``expr tspan/100``).
        - A numpy array, ``PerLC``, or ``pd.Series`` for per-LC batch
          values.
    subsample, finetune : float, str, numpy array, PerLC, or pd.Series
        Coarse-grid step and fine-tune resolution as fractions of 1/T
        (``T`` = time baseline).  Same value forms as *minp* / *maxp*.
    npeaks : int
        Number of peaks to report (default 5).
    save_periodogram : bool, str, or Output
        Auxiliary file output.  ``True`` captures as
        ``result.files["ftp_periodogram_N"]``; a path string writes to
        that directory without capturing; ``Output(path, capture=True)``
        does both.  ``False`` (default) suppresses.  With ``whiten=True``
        the file holds one column per whitening cycle.

    template_file : str, optional (``file`` mode)
        Path to a two-column ``c_n s_n`` whitespace-separated text file.
    lc_path : str, optional (``fitlc`` mode)
        Path to a light curve from which the template will be built.
    lc_format : {"ascii", "fits"}, optional (``fitlc`` mode)
        Format of ``lc_path``.
    t_col, mag_col, err_col : int or str, optional (``fitlc`` mode)
        ASCII: 1-indexed integer column numbers; pass ``err_col=0`` for an
        unweighted fit.  FITS: column-name strings; pass
        ``err_col="none"`` or ``""`` for an unweighted fit.
    nharm : int, optional (``fitlc`` or ``inline`` mode)
        Harmonics ABOVE the fundamental (matches the ``-harmonicfilter`` /
        ``-Injectharm`` convention).  Total template harmonic count is
        ``nharm + 1``.
    period : float, optional (``fitlc`` mode)
        The fixed period at which the template Fourier series is fit.
        Numeric only (the vartools C parser uses ``atof``; var/expr are not
        accepted on this slot).
    cn, sn : list of (number or str), optional (``inline`` mode)
        Length-(``nharm + 1``) lists of the c_n and s_n coefficients.  Each
        entry is either a number or a string identifier / expression
        evaluated per LC.  ``nharm`` is inferred from ``len(cn)``.
    filelist_column : int, optional (``filelist`` mode)
        1-indexed column in the ``-l`` input list holding each LC's template
        path.  If omitted, the next available column is used.

    clip : float, optional
        Sigma-clip factor for the SNR noise estimate (default 5).
    clipiter : int, optional
        0 or 1: disable / enable iterative clipping (default 1).
    noerr : bool
        Force uniform weights instead of 1/sigma^2.
    posamponly : bool
        Skip negative-amplitude solutions during the search (a flipped
        template is generally not a valid match).
    whiten : bool
        Iterative pre-whitening between peaks: subtract the closed-form
        FTP-template fit from the LC at each peak before searching for
        the next.
    fixperiod_snr : float, int, or str, optional
        Additionally report the FTP statistic at a specified period.
        Recognised back-references: ``"aov"`` / ``"ls"`` / ``"pdm"`` /
        ``"ftp"`` / ``"injectharm"``, ``"fixcolumn <name>"``,
        ``"list [column N]"``, or a numeric value (treated as a literal).
    bootstrap : int, optional
        Empirical-CDF FAP via ``Nboot`` shuffled-LC trials.
    maskpoints : str, optional
        Name of an LC vector; points with maskvar > VARTOOLS_MASK_TINY are
        included, others excluded.
    method : {"auto", "brute", "poly", "verify"}, optional
        Per-frequency optimisation strategy.  ``auto`` (default) picks
        ``poly`` for H <= 2, else ``brute``.
    sums : {"auto", "direct", "nfft"}, optional
        Per-LC summation strategy.  ``auto`` (default) picks ``nfft`` if
        vartools was built with ``--with-nfft``, else ``direct``.

    See Also
    --------
    CLI command: ``-FTP``.
    Citations: Hoffman, VanderPlas, Hartman & Bakos 2021
    (arXiv:2101.12348).  Reference Python implementation:
    https://github.com/PrincetonUniversity/FastTemplatePeriodogram
    """

    _vt_name = "FTP"

    def __init__(
        self,
        template_source: str,
        minp: Union[float, str],
        maxp: Union[float, str],
        subsample: Union[float, str],
        finetune: Union[float, str],
        *,
        # file mode
        template_file: Optional[str] = None,
        # fitlc mode
        lc_path: Optional[str] = None,
        lc_format: Optional[str] = None,
        t_col: Optional[Union[int, str]] = None,
        mag_col: Optional[Union[int, str]] = None,
        err_col: Optional[Union[int, str]] = None,
        period: Optional[Union[float, int]] = None,
        # fitlc + inline modes
        nharm: Optional[int] = None,
        # inline mode
        cn: Optional[List] = None,
        sn: Optional[List] = None,
        # filelist mode
        filelist_column: Optional[int] = None,
        # required after positionals
        npeaks: int = 5,
        save_periodogram=False,
        # trailing
        clip: Optional[float] = None,
        clipiter: Optional[int] = None,
        noerr: bool = False,
        posamponly: bool = False,
        whiten: bool = False,
        fixperiod_snr: Union[float, int, str, None] = None,
        bootstrap: Optional[int] = None,
        maskpoints: Optional[str] = None,
        method: Optional[str] = None,
        sums: Optional[str] = None,
    ) -> None:
        if template_source not in _FTP_TEMPLATE_SOURCES:
            raise ValueError(
                f"FTP: template_source must be one of {_FTP_TEMPLATE_SOURCES!r}, "
                f"got {template_source!r}"
            )

        # Collect mode-specific kwargs to enforce mutual exclusion.
        _file_kw = (template_file,)
        _fitlc_kw = (lc_path, lc_format, t_col, mag_col, err_col, period)
        _inline_kw = (cn, sn)
        _filelist_kw = (filelist_column,)

        if template_source == "file":
            if template_file is None:
                raise ValueError("FTP file mode requires 'template_file'")
            if any(x is not None for x in (*_fitlc_kw, *_inline_kw,
                                            *_filelist_kw)):
                raise ValueError(
                    "FTP file mode rejects fitlc/inline/filelist kwargs"
                )
        elif template_source == "fitlc":
            if any(x is None for x in (lc_path, lc_format, t_col, mag_col,
                                        err_col, nharm, period)):
                raise ValueError(
                    "FTP fitlc mode requires lc_path, lc_format, t_col, "
                    "mag_col, err_col, nharm, and period"
                )
            if lc_format not in ("ascii", "fits"):
                raise ValueError(
                    f"FTP fitlc: lc_format must be 'ascii' or 'fits', got "
                    f"{lc_format!r}"
                )
            if not isinstance(nharm, int) or nharm < 0:
                raise ValueError(
                    f"FTP fitlc: nharm must be a non-negative int, got "
                    f"{nharm!r}"
                )
            if not isinstance(period, (int, float)):
                raise ValueError(
                    f"FTP fitlc: period must be numeric (vartools parses it "
                    f"as a literal float), got {period!r}"
                )
            if any(x is not None for x in (*_file_kw, *_inline_kw,
                                            *_filelist_kw)):
                raise ValueError(
                    "FTP fitlc mode rejects file/inline/filelist kwargs"
                )
        elif template_source == "inline":
            if cn is None or sn is None:
                raise ValueError("FTP inline mode requires cn and sn")
            if len(cn) != len(sn) or len(cn) < 1:
                raise ValueError(
                    f"FTP inline: cn and sn must have equal non-zero length, "
                    f"got len(cn)={len(cn)} len(sn)={len(sn)}"
                )
            if nharm is not None and nharm != len(cn) - 1:
                raise ValueError(
                    f"FTP inline: nharm={nharm} disagrees with len(cn)-1="
                    f"{len(cn) - 1}; either pass nharm or omit it (it is "
                    "inferred from len(cn))"
                )
            if any(x is not None for x in (*_file_kw, *_fitlc_kw,
                                            *_filelist_kw)):
                raise ValueError(
                    "FTP inline mode rejects file/fitlc/filelist kwargs"
                )
            nharm = len(cn) - 1
        elif template_source == "filelist":
            if filelist_column is not None and (
                not isinstance(filelist_column, int) or filelist_column < 1
            ):
                raise ValueError(
                    f"FTP filelist: filelist_column must be a 1-indexed int, "
                    f"got {filelist_column!r}"
                )
            if any(x is not None for x in (*_file_kw, *_fitlc_kw,
                                            *_inline_kw[:1], nharm)):
                raise ValueError(
                    "FTP filelist mode rejects file/fitlc/inline kwargs"
                )

        if method is not None and method not in _FTP_METHODS:
            raise ValueError(
                f"FTP: method must be one of {_FTP_METHODS!r}, got {method!r}"
            )
        if sums is not None and sums not in _FTP_SUMS:
            raise ValueError(
                f"FTP: sums must be one of {_FTP_SUMS!r}, got {sums!r}"
            )
        if bootstrap is not None:
            try:
                _nboot = int(bootstrap)
            except (TypeError, ValueError):
                raise ValueError(
                    f"FTP: bootstrap must be a positive int, got {bootstrap!r}"
                )
            if _nboot < 1:
                raise ValueError(
                    f"FTP: bootstrap must be >= 1, got {_nboot}"
                )
            bootstrap = _nboot

        self.template_source = template_source
        self.minp = minp
        self.maxp = maxp
        self.subsample = subsample
        self.finetune = finetune
        self.template_file = template_file
        self.lc_path = lc_path
        self.lc_format = lc_format
        self.t_col = t_col
        self.mag_col = mag_col
        self.err_col = err_col
        self.period = period
        self.nharm = nharm
        self.cn = cn
        self.sn = sn
        self.filelist_column = filelist_column
        self.npeaks = npeaks
        self.save_periodogram = save_periodogram
        self.clip = clip
        self.clipiter = clipiter
        self.noerr = noerr
        self.posamponly = posamponly
        self.whiten = whiten
        self.fixperiod_snr = fixperiod_snr
        self.bootstrap = bootstrap
        self.maskpoints = maskpoints
        self.method = method
        self.sums = sums

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-FTP"]
        if self.template_source == "file":
            args += ["file", str(self.template_file)]
        elif self.template_source == "fitlc":
            args += ["fitlc", str(self.lc_path), str(self.lc_format),
                     str(self.t_col), str(self.mag_col), str(self.err_col),
                     str(int(self.nharm)), str(float(self.period))]
        elif self.template_source == "inline":
            args += ["inline", str(int(self.nharm))]
            for c_val, s_val in zip(self.cn, self.sn):
                args += _varexpr(c_val) + _varexpr(s_val)
        elif self.template_source == "filelist":
            args += ["filelist"]
            if self.filelist_column is not None:
                args += ["column", str(int(self.filelist_column))]
        # Required positionals.
        args += (_varexpr(self.minp) + _varexpr(self.maxp)
                 + _varexpr(self.subsample) + _varexpr(self.finetune)
                 + [str(self.npeaks)])
        args += _outtoken(self.save_periodogram, outdir)
        # Trailing keywords in the canonical (strict-parser) order:
        # clip / noerr / posamponly / whiten / fixperiodSNR / bootstrap /
        # maskpoints / method / sums.
        if self.clip is not None:
            args += ["clip", str(self.clip),
                     str(self.clipiter if self.clipiter is not None else 1)]
        args += _bool("noerr", self.noerr)
        args += _bool("posamponly", self.posamponly)
        args += _bool("whiten", self.whiten)
        args += _fixperiodsnr_tokens(self.fixperiod_snr)
        if self.bootstrap is not None:
            args += ["bootstrap", str(int(self.bootstrap))]
        args += _flag("maskpoints", self.maskpoints)
        args += _flag("method", self.method)
        args += _flag("sums", self.sums)
        return args

    def _resolve_back_references(self, prev) -> None:
        from ._helpers import _resolve_period_backref
        self.fixperiod_snr = _resolve_period_backref(prev, self.fixperiod_snr)

    def _output_file_specs(self):
        return {"periodogram": (".ftp", None)}


class BLS(VartoolsCommand):
    """Box-fitting Least Squares (BLS) transit search.

    Search a grid of trial periods and phase bins for periodic
    box-shaped (or trapezoidal) dips.  Two grids must be specified: the
    transit-duration grid (one of ``r``, ``q``, ``density`` modes) and
    the trial-frequency grid (one of ``optimal``, ``nfreq``, ``df``
    modes).

    Transit-duration grid:

    - **r mode** (default): pass ``rmin`` / ``rmax`` as the
      minimum/maximum stellar radius in solar radii.  The fractional-
      duration range at each period is derived from
      ``q = 0.076 · R^(2/3) · P^(-2/3)``.
    - **q mode**: pass ``qmin`` / ``qmax`` directly as the fractional
      transit duration.
    - **density mode**: set ``density_mode=True`` and supply
      ``stellar_density`` (g/cm³) plus ``min_exp_dur_frac`` /
      ``max_exp_dur_frac`` to bracket the expected circular-orbit
      duration.

    Trial-frequency grid (mutually exclusive — pass exactly one of
    ``subsample``, ``nfreq``, ``df``):

    - **optimal mode**: Ofir 2014 frequency sampling optimal for transit
      search, controlled by ``subsample`` (oversampling factor).
      Available only with ``density_mode=True``.
    - **nfreq mode**: ``nfreq=N`` for a fixed number of trial
      frequencies on a uniform grid.
    - **df mode**: ``df=Δf`` for a fixed frequency step on a uniform
      grid.

    When ``density_mode=True``, ``optimal`` mode is the default unless
    ``nfreq`` or ``df`` is also set; in ``r`` / ``q`` duration mode
    ``nfreq`` or ``df`` is required (the wrapper raises ``ValueError``
    at construction time if both are omitted).

    Parameters
    ----------
    minper, maxper : float, str, numpy array, PerLC, or pd.Series
        Period search range (days).  Forms:

        - A number — passed directly as a fixed value.
        - A bare identifier string (e.g. ``"minperiod"``) — vartools
          reads from a named per-LC variable (``var minperiod``).
        - Any other string (e.g. ``"tspan/100"``) — evaluated as a math
          expression (``expr tspan/100``).
        - A numpy array, ``PerLC``, or ``pd.Series`` for per-LC batch
          values.
    rmin, rmax : float or str
        ``r``-mode duration bounds (default mode).  Ignored when
        ``qmin`` / ``qmax`` or ``density_mode`` is set.
    qmin, qmax : float, str, or None, optional
        ``q``-mode duration bounds (fractional transit duration).  When
        set, emits ``"q" qmin qmax`` instead of ``"r" rmin rmax``.
    density_mode : bool
        Use stellar density to set transit-duration bounds.  When
        ``True``, ``stellar_density``, ``min_exp_dur_frac``, and
        ``max_exp_dur_frac`` define the duration range.
    stellar_density : float, str, or None, optional
        Stellar density (g/cm³) for density mode.
    min_exp_dur_frac, max_exp_dur_frac : float or str
        Expected-duration fractions for density mode.  Defaults 0.5 and
        1.5.
    nbins : int or str
        Number of phase bins (≥ ``2/qmin``).  Default 200.  Accepts
        var/expr/PerLC forms.
    timezone : float
        Time-zone offset (0 for HJD/BJD); affects the single-night
        Δχ² fraction.
    npeaks : int
        Number of transit candidates to report.  Default 1.
    subsample : float or str
        Oversampling factor for the Ofir-2014 optimal frequency-
        sampling method.  Used only when ``density_mode=True`` and
        neither ``nfreq`` nor ``df`` is set.  Mutually exclusive with
        ``nfreq`` and ``df``.
    nfreq : int, str, or None, optional
        Fixed number of trial frequencies (uniform grid).  Mutually
        exclusive with ``subsample`` and ``df``.
    df : float, str, or None, optional
        Fixed frequency step (uniform grid).  Mutually exclusive with
        ``subsample`` and ``nfreq``.
    save_periodogram : bool, str, or Output
        BLS spectrum file.  ``True`` captures as
        ``result.files["BLS_periodogram_N"]``; a path string writes to
        that directory without capturing; ``Output(path, capture=True)``
        does both.
    save_model : bool, str, or Output
        Best-fit transit model.  Same value semantics as
        ``save_periodogram``; key ``BLS_model_N``.
    correct_lc : bool
        Subtract the best-fit transit from the light curve before
        passing to the next command.
    extraparams : bool
        Include additional false-positive diagnostic columns in the
        output (see ``BLS_SRSum_k_N`` etc.).
    fittrap : bool
        Fit a trapezoidal rather than box-shaped transit at each peak.
        Adds ``BLS_Qingress_k_N`` and ``BLS_OOTmag_k_N`` to the output.
    nobinnedrms : bool
        Adjust the way in which the BLS_SN statistic is calculated.
        The default mode of ``True`` yields a faster and more robust
        process.  Set to ``False`` to recover the historical VARTOOLS
        behavior.
    save_phcurve : bool, str, or Output
        Phase-folded model curve.  Same value semantics as
        ``save_periodogram``; key ``BLS_phcurve_N``.
    ophcurve_phmin, ophcurve_phmax, ophcurve_phstep : float
        Phase range and step for the phase-curve output.  Defaults
        0.0, 1.0, 0.005.
    save_jdcurve : bool, str, or Output
        JD-sampled model curve.  Same value semantics as
        ``save_periodogram``; key ``BLS_jdcurve_N``.
    ojdcurve_jdstep : float
        Time step (days) for the JD-curve output.  Default 0.02.
    freq_grid : str, optional
        ``"stepP"`` for uniform period sampling, ``"steplogP"`` for
        log-uniform.
    adjust_qmin : bool
        Adaptively increase ``qmin`` at each frequency to
        ``max(qmin, mindt·f)``.
    reduce_nbins : bool
        With ``adjust_qmin=True``, adaptively reduce ``nbins`` at each
        frequency.
    reportharmonics : bool
        Report period harmonics (½, ⅓, …) as additional candidates.
    maskpoints : str, optional
        Mask variable; points with ``maskvar ≤ 0`` are excluded from
        the BLS spectrum.

    See Also
    --------
    CLI command: ``-BLS``.
    Citations: Kovács, Zucker & Mazeh 2002 (A&A, 391, 369);
    Ofir 2014 (A&A, 561, A138) for the optimal frequency sampling.
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
        nobinnedrms: bool = True,
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
    """BLS at a fixed (pre-known) period.

    Search only the transit phase, depth, and duration at a single
    fixed period; no period grid is scanned.  Useful as a second pass
    after a full period search.

    Parameters
    ----------
    period : float or str
        Fixed period (days).  Forms:

        - A number — passed as ``"fix <value>"``.
        - ``"ls"`` / ``"aov"`` — back-reference to the best period from
          a prior LS or AOV command (works inside a single Pipeline or
          across chain steps).
        - ``"fixcolumn COLNAME"`` — read the period from a named
          per-star column.
        - ``"list"`` / ``"list column N"`` — read the period from a
          list-file column.

        A missing prior command for a back-reference raises
        ``LookupError``.
    rmin, rmax : float
        ``r``-mode duration bounds in solar radii.  Defaults 0.01 and
        0.1.  Ignored when ``qmin`` / ``qmax`` is set.
    qmin, qmax : float, optional
        ``q``-mode duration bounds (fractional transit duration).  When
        set, emits ``"q" qmin qmax`` instead of ``"r" rmin rmax``.
    nbins : int
        Number of phase bins.  Default 200.
    timezone : float
        Time-zone offset (0 for HJD/BJD); affects the single-night
        Δχ² fraction.
    save_model : bool, str, or Output
        Best-fit transit model.  ``True`` captures as
        ``result.files["BLSFixPer_model_N"]``; a path string writes to
        that directory without capturing; ``Output(path, capture=True)``
        does both.
    correct_lc : bool
        Subtract the best-fit transit from the light curve before
        passing to the next command.
    fittrap : bool
        Fit a trapezoidal transit shape instead of a box.
    maskpoints : str, optional
        Mask variable; points with ``maskvar ≤ 0`` are excluded.

    See Also
    --------
    CLI command: ``-BLSFixPer``.
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

    Search a period grid from ``minper`` to ``maxper`` with the
    transit duration and reference epoch held fixed; optionally also
    fix the transit depth and ingress fraction.

    Parameters
    ----------
    duration : float or str
        Transit duration (days).  Forms:

        - A number — passed as ``"fix <value>"``.
        - ``"fixcolumn <colname>"`` — read from a named per-star column.
        - ``"list"`` / ``"list column N"`` — read from a list-file column.
    Tc : float or str
        Mid-transit epoch (JD/BJD).  Same value forms as *duration*.
    minper, maxper : float
        Period search range (days).  Defaults 0.1 and 100.0.
    nfreq : int
        Number of trial frequencies.  Default 10000.
    timezone : float
        Time-zone offset (0 for UTC/BJD).
    npeaks : int
        Number of peaks to report.  Default 1.
    fixdepth : float, str, or None, optional
        Fix the transit depth to this value (or a ``"fixcolumn"`` /
        ``"list"`` spec).  When ``None`` the depth is optimised.
    qgress : float, str, or None, optional
        Fractional ingress/egress duration (``0`` = box, ``0.5`` =
        V-shaped/grazing).  Only meaningful when *fixdepth* is set.
    save_periodogram : bool, str, or Output
        BLS spectrum file.  ``True`` captures as
        ``result.files["BLSFixDurTc_periodogram_N"]``; a path string
        writes to that directory without capturing;
        ``Output(path, capture=True)`` does both.
    save_model : bool, str, or Output
        Best-fit transit model.  Same value semantics as
        ``save_periodogram``; key ``BLSFixDurTc_model_N``.
    correct_lc : bool
        Subtract the best-fit transit from the light curve before
        passing to the next command.
    fittrap : bool
        Fit a trapezoidal transit shape instead of a box.
    save_phcurve : bool, str, or Output
        Phase-folded model curve.  Same value semantics as
        ``save_periodogram``; key ``BLSFixDurTc_phcurve_N``.
    ophcurve_phmin, ophcurve_phmax, ophcurve_phstep : float
        Phase range and step for the phase-curve output.  Defaults
        0.0, 1.0, 0.005.
    save_jdcurve : bool, str, or Output
        JD-sampled model curve.  Same value semantics as
        ``save_periodogram``; key ``BLSFixDurTc_jdcurve_N``.
    ojdcurve_jdstep : float
        Time step (days) for the JD-curve output.  Default 0.02.
    maskpoints : str, optional
        Mask variable; points with ``maskvar ≤ 0`` are excluded.

    See Also
    --------
    CLI command: ``-BLSFixDurTc``.
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

    Compute BLS transit statistics for a fully specified signal — no
    period search is performed.  The depth is optimised by default
    (or fixed when *fixdepth* is given).

    Parameters
    ----------
    period : float or str
        Transit period (days).  Forms:

        - A number — passed as ``"fix <value>"``.
        - ``"fixcolumn <colname>"`` — read from a named per-star column.
        - ``"list"`` / ``"list column N"`` — read from a list-file column.
    duration : float or str
        Transit duration (days).  Same value forms as *period*.
    Tc : float or str
        Mid-transit epoch (JD/BJD).  Same value forms as *period*.
    timezone : float
        Time-zone offset (0 for UTC/BJD).
    fixdepth : float, str, or None, optional
        Fix the transit depth to this value (or a ``"fixcolumn"`` /
        ``"list"`` spec).  When ``None`` the depth is optimised.
    qgress : float, str, or None, optional
        Fractional ingress/egress duration.  Only meaningful when
        *fixdepth* is set.
    save_model : bool, str, or Output
        Best-fit transit model.  ``True`` captures as
        ``result.files["BLSFixPerDurTc_model_N"]``; a path string writes
        to that directory without capturing;
        ``Output(path, capture=True)`` does both.
    correct_lc : bool
        Subtract the best-fit transit from the light curve before
        passing to the next command.
    fittrap : bool
        Fit a trapezoidal transit shape instead of a box.
    save_phcurve : bool, str, or Output
        Phase-folded model curve.  Same value semantics as
        ``save_model``; key ``BLSFixPerDurTc_phcurve_N``.
    ophcurve_phmin, ophcurve_phmax, ophcurve_phstep : float
        Phase range and step for the phase-curve output.  Defaults
        0.0, 1.0, 0.005.
    save_jdcurve : bool, str, or Output
        JD-sampled model curve.  Same value semantics as ``save_model``;
        key ``BLSFixPerDurTc_jdcurve_N``.
    ojdcurve_jdstep : float
        Time step (days) for the JD-curve output.  Default 0.02.
    maskpoints : str, optional
        Mask variable; points with ``maskvar ≤ 0`` are excluded.

    See Also
    --------
    CLI command: ``-BLSFixPerDurTc``.
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

    The CLI always writes the autocorrelation function to a file —
    there is no mode that suppresses file output entirely.

    Parameters
    ----------
    start : float
        Start of the lag range (same units as the time column).
    stop : float
        End of the lag range.
    step : float
        Lag step size.
    save_result : bool, str, or Output
        Controls capture of the output file into Python.

        - ``True`` (default) — write to a temp dir and capture into
          ``result.files["autocorrelation_result_N"]``.
        - ``False`` — write to a temp dir but do **not** capture into
          Python (the file is still written because the CLI always does
          so).
        - A directory path string — write to that directory, no
          capture.
        - ``Output(path, capture=True)`` — write to that directory and
          capture.
    maskpoints : str, optional
        Mask variable; points with ``maskvar ≤ 0`` are excluded.

    See Also
    --------
    CLI command: ``-autocorrelation``.
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
    """Discrete Fourier Transform power spectrum + CLEAN deconvolution.

    Compute the DFT power spectrum of the light curve using the FDFT
    algorithm of Kurtz 1985, and optionally deconvolve it with the
    CLEAN algorithm of Roberts, Lehar & Dreher 1987 to remove aliasing
    due to the window function.  The CLEAN iteration subtracts a
    gain-scaled CLEAN beam at the strongest peak and repeats until the
    residual falls below ``SNlimit · noise``.

    Parameters
    ----------
    nbeam : int or str
        Number of frequency samples per ``1/T`` element (``T`` = time
        baseline).  Controls spectral resolution.
    maxfreq : float, str, or None, optional
        Maximum frequency (cycles/day).  ``None`` uses
        ``1 / (2 · min_time_separation)`` (Nyquist).
    save_dspec : bool, str, or Output
        Dirty (DFT) spectrum.  ``True`` captures as
        ``result.files["dftclean_dspec_N"]``; a path string writes to
        that directory without capturing; ``Output(path, capture=True)``
        does both.
    finddirtypeaks : int, optional
        Find the top N peaks in the dirty power spectrum.
    finddirtypeaks_clip : float, optional
        Sigma-clipping factor for the dirty-peak SNR noise estimate
        (requires *finddirtypeaks*; default 5σ).
    finddirtypeaks_clipiter : int, optional
        ``0`` = single-pass clipping, ``1`` = iterative clipping
        (requires *finddirtypeaks_clip*; default 1).
    save_wfunc : bool, str, or Output
        Window function.  Same value semantics as ``save_dspec``; key
        ``dftclean_wfunc_N``.
    save_cspec : bool, str, or Output
        CLEAN spectrum (activates the CLEAN iteration when set).  Same
        value semantics as ``save_dspec``; key ``dftclean_cspec_N``.
    gain : float
        CLEAN gain factor in ``[0.1, 1.0]`` controlling how much of each
        peak is subtracted per CLEAN iteration; smaller is slower but
        more thorough.  Default 0.1.
    SNlimit : float
        Stop CLEANing when the peak falls below this S/N threshold.
        Default 3.0.
    outcbeam : bool, str, or Output
        CLEAN beam.  Same value semantics as ``save_dspec``; key
        ``dftclean_cbeam_N``.  Requires the CLEAN section to be active
        (set ``save_cspec``, ``npeaks``, or ``finddirtypeaks``).
    npeaks : int, optional
        Number of peaks to find in the CLEAN spectrum.  Activates
        CLEAN when set.
    useampspec : bool
        Compute SNR on the amplitude spectrum instead of the power
        spectrum.
    verboseout : bool
        Include the mean and stddev of the spectrum (before and after
        clipping) in the output table.
    maskpoints : str, optional
        Mask variable; points with ``maskvar ≤ 0`` are excluded.

    See Also
    --------
    CLI command: ``-dftclean``.
    Citations: Kurtz 1985 (MNRAS 213, 773) for the FDFT algorithm;
    Roberts, Lehar & Dreher 1987 (AJ 93, 968) for CLEAN.
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

    Compute the WWZ of Foster 1996 with an abbreviated Morlet wavelet
    ``f(z) = exp(i·2π·f·(t − τ) − c·(2π·f)²·(t − τ)²)``.  The transform
    is computed for all combinations of trial frequency (up to
    ``maxfreq``) and time shift (``tau0`` to ``tau1`` in steps of
    ``dtau``).  Especially useful for non-stationary signals.

    Parameters
    ----------
    maxfreq : float, str, or ``"auto"``
        Maximum frequency in cycles per day.  ``"auto"`` (default) sets
        it to ``1 / (2 · delmin)`` where ``delmin`` is the minimum
        spacing between consecutive observations.
    freqsamp : float or str
        Frequency-sampling step factor; the actual frequency step is
        ``freqsamp / T`` where ``T`` is the time baseline.  Default
        ``0.25`` (Foster 1996 convention).  vartools does **not**
        accept ``"auto"`` for this parameter.
    tau0, tau1, dtau : float, str, or ``"auto"``
        Start, end, and step of the time-shift scan.  ``"auto"`` for
        *tau0* / *tau1* uses the LC's first/last time; ``"auto"`` for
        *dtau* uses ``delmin``.
    c : float or str
        Decay constant of the abbreviated Morlet wavelet (default
        ``1 / (8π²) ≈ 0.0125``, the Foster recommendation).  Controls
        the trade-off between time and frequency resolution.
    save_transform : bool, str, or Output
        Full WWZ time-frequency map.  ``True`` captures as
        ``result.files["wwz_transform_N"]``; a path string writes to
        that directory without capturing; ``Output(path, capture=True)``
        does both.
    transform_format : str, optional
        Output format for the full transform: ``"fits"`` or ``"pm3d"``.
        Only used when *save_transform* is set.
    transform_name : str, optional
        Filename format string for the full transform output (e.g.
        ``"%s.wwz"``).
    save_maxtransform : bool, str, or Output
        WWZ max-power projection over frequency.  Same value semantics
        as ``save_transform``; key ``wwz_maxtransform_N``.
    maxtransform_name : str, optional
        Filename format string for the max-transform output.
    maskpoints : str, optional
        Mask variable; points with ``maskvar ≤ 0`` are excluded.

    See Also
    --------
    CLI command: ``-wwz``.
    Citation: Foster 1996 (AJ 112, 1709).
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
    """Minimum detectable Lomb-Scargle amplitude for injection/recovery.

    Determine the minimum peak-to-peak amplitude that a signal at a
    given period must have to be detected by a Lomb-Scargle search with
    ``−ln(FAP) > thresh``.  The signal shape is either a Fourier
    series (``mode="harm"``) or read from a file (``mode="file"``); the
    threshold is found by scaling the template until the LS statistic
    reaches the detection limit.

    Parameters
    ----------
    period : str
        Reference period source.  The CLI accepts only ``"ls"`` (use
        the period from the most recent prior ``-LS``) or ``"list"`` —
        a bare numeric value is **not** accepted.  Default ``"ls"``.

        ``"ls"`` is meaningful only in single-Pipeline usage where the
        matching ``-LS`` is in the same vartools invocation.  Across a
        chain boundary the lookup is not supported and pyvartools raises
        ``NotImplementedError`` — fold the ``-LS`` step into the same
        Pipeline if you hit this.
    minp : float
        Minimum search period (days) that would be used in the LS
        search; sets the FAP scale.  Default 0.1.
    thresh : float
        Desired ``−ln(FAP)`` detection threshold.  Default 10.0.
    mode : str
        Signal model:

        - ``"harm"`` (default) — Fourier series with *nharm* and
          *nsubharm* terms.
        - ``"file"`` — signal template is read from *listfile*.
    nharm : int
        Number of harmonics (only when ``mode="harm"``).  Default 1.
    nsubharm : int
        Number of sub-harmonics (only when ``mode="harm"``).  Default 0.
    listfile : str, optional
        Path to the template-signal list file (required when
        ``mode="file"``).
    noGLS : bool
        Use the classical (non-generalised) Lomb-Scargle statistic.

    See Also
    --------
    CLI command: ``-GetLSAmpThresh``.
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

    Replaces the time variable with phase relative to a chosen period
    and (optionally) a reference epoch.  The original ``t`` column is
    overwritten with the folded phase.

    Parameters
    ----------
    period : float or str
        Period to fold on (days).  Forms:

        - A number — passed as ``"fix <value>"``.
        - ``"ls"`` / ``"aov"`` / ``"bls"`` / ``"blsfixper"`` /
          ``"injectharm"`` — back-reference to the best period from a
          prior command of that type (works inside a single Pipeline or
          across chain steps).
        - ``"fixcolumn COLNAME"`` — read the period from a named
          per-star column.
        - ``"list"`` / ``"list column N"`` — read from a list-file
          column.

        Default ``"ls"``.  A missing prior command for a back-reference
        raises ``LookupError``.
    T0 : float or str, optional
        Reference epoch.  Forms:

        - A number — passed as ``"fix <value>"``.
        - ``"bls"`` or ``"bls <phaseTc>"`` — use the prior BLS's
          ``Tc`` (optionally shifted by a phase offset:
          ``Tc - phaseTc · Period``).
        - ``"fix <value>"`` — same as a bare number.
        - ``"fixcolumn COLNAME"`` / ``"list"`` — read from external data.

        When ``None``, vartools does not subtract a reference time before
        folding.
    phasevar : str, optional
        Name of the output phase variable.  When set, the per-point
        phase is stored under this name in addition to overwriting the
        ``t`` column.
    startphase : float, optional
        Starting phase value (default 0; phase wraps into
        ``[startphase, startphase + 1)``).

    See Also
    --------
    CLI command: ``-Phase``.
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
