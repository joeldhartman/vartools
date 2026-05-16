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
    minp, maxp : float or str
        Period search range.  Accepts the same number / variable-name /
        expression forms as :class:`aov`.
    subsample, finetune : float or str
        Frequency-grid step and fine-tune resolution.  Accepts var/expr.
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
    save_periodogram : bool or str
        ``True`` captures the periodogram in memory; pass a path string to
        write it to disk; ``False`` (default) suppresses.
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
    minp, maxp : float or str
        Period search range.  Accepts var/expr forms.
    subsample, finetune : float or str
        Frequency-grid step and fine-tune resolution.  Accepts var/expr.
    npeaks : int
        Number of peaks to report (default 5).
    save_periodogram : bool or str
        ``True`` writes the periodogram to the pipeline outdir with suffix
        ``.ftp``; pass a path string to write to a specific directory;
        ``False`` (default) suppresses.

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
