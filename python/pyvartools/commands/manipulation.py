"""Light curve manipulation and statistics command wrappers."""

from __future__ import annotations
import re
from typing import List, Optional, Sequence, Union

from pyvartools._command import VartoolsCommand
from ._helpers import _bool, _flag, _injectparam, _norm_save, _outtoken, _period_spec, _pval, _should_emit, _varexpr


class clip(VartoolsCommand):
    """Sigma-clip outliers from the light curve.

    Parameters
    ----------
    sigclip : float
        Clipping threshold in sigma.
    iterative : bool
        If True (default), repeat clipping until no more points are removed
        (passes ``iter=1`` to vartools).  If False, clip only once.
    niter : int, optional
        If given, clip exactly this many times (passes ``"niter" n``).
        Ignored when ``iterative=False`` and ``niter`` is None.
    median : bool
        Clip relative to the median (default: clip relative to the mean).
    markclip : str, optional
        Variable name to record clipping mask (1=kept, 0=clipped).
    noinitmark : bool
        Do not initialise the markclip variable before clipping.
    maskpoints : str, optional
        Name of a mask variable; points with ``maskvar ≤ 0`` are
        excluded from the clipping statistics.

    See Also
    --------
    CLI command: ``-clip``.
    """

    _vt_name = "clip"

    def __init__(
        self,
        sigclip: float,
        iterative: bool = True,
        niter: Optional[int] = None,
        median: bool = False,
        markclip: Optional[str] = None,
        noinitmark: bool = False,
        maskpoints: Optional[str] = None,
    ) -> None:
        self.sigclip = sigclip
        self.iterative = iterative
        self.niter = niter
        self.median = median
        self.markclip = markclip
        self.noinitmark = noinitmark
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        args = ["-clip"] + _varexpr(self.sigclip) + _varexpr(1 if self.iterative else 0)
        if self.niter is not None:
            args += ["niter"] + _varexpr(self.niter)
        args += _bool("median", self.median)
        if self.markclip is not None:
            args += ["markclip", self.markclip]
            args += _bool("noinitmark", self.noinitmark)
        args += _flag("maskpoints", self.maskpoints)
        return args


class rms(VartoolsCommand):
    """Compute the RMS and expected RMS of the light curve.

    Emits the unweighted magnitude RMS, the error-weighted mean
    magnitude, and the expected RMS implied by the per-point
    uncertainties.

    Parameters
    ----------
    maskpoints : str, optional
        Name of a mask variable; points with ``maskvar ≤ 0`` are
        excluded.

    See Also
    --------
    CLI command: ``-rms``.
    """

    _vt_name = "rms"

    def __init__(self, maskpoints: Optional[str] = None) -> None:
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        return ["-rms"] + _flag("maskpoints", self.maskpoints)


class rmsbin(VartoolsCommand):
    """Compute the RMS of the light curve after moving-average binning.

    For each timescale the magnitudes are smoothed with a moving
    average of that width and the RMS of the smoothed series is
    reported — a measure of how much scatter survives binning.

    Parameters
    ----------
    nbin : int
        Number of timescales (must equal ``len(bintimes)``).
    bintimes : list of float
        Bin timescales (days) at which to compute the binned RMS.
    maskpoints : str, optional
        Name of a mask variable; points with ``maskvar ≤ 0`` are
        excluded.

    See Also
    --------
    CLI command: ``-rmsbin``.
    """

    _vt_name = "rmsbin"

    def __init__(
        self,
        nbin: int,
        bintimes: List[float],
        maskpoints: Optional[str] = None,
    ) -> None:
        self.nbin = nbin
        self.bintimes = bintimes
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        args = ["-rmsbin", str(self.nbin)] + [str(b) for b in self.bintimes]
        args += _flag("maskpoints", self.maskpoints)
        return args


class chi2(VartoolsCommand):
    """Compute the reduced chi-squared of the light curve about its
    weighted mean.

    Reports ``χ²/dof`` and the error-weighted mean magnitude.

    Parameters
    ----------
    maskpoints : str, optional
        Name of a mask variable; points with ``maskvar ≤ 0`` are
        excluded.

    See Also
    --------
    CLI command: ``-chi2``.
    """

    _vt_name = "chi2"

    def __init__(self, maskpoints: Optional[str] = None) -> None:
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        return ["-chi2"] + _flag("maskpoints", self.maskpoints)


class chi2bin(VartoolsCommand):
    """Compute reduced chi-squared after moving-average binning.

    For each timescale the magnitudes are smoothed with a moving
    average of that width and ``χ²/dof`` of the smoothed series is
    reported.

    Parameters
    ----------
    nbin : int
        Number of timescales (must equal ``len(bintimes)``).
    bintimes : list of float
        Bin timescales (days) at which to compute the binned χ².
    maskpoints : str, optional
        Name of a mask variable; points with ``maskvar ≤ 0`` are
        excluded.

    See Also
    --------
    CLI command: ``-chi2bin``.
    """

    _vt_name = "chi2bin"

    def __init__(
        self,
        nbin: int,
        bintimes: List[float],
        maskpoints: Optional[str] = None,
    ) -> None:
        self.nbin = nbin
        self.bintimes = bintimes
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        args = ["-chi2bin", str(self.nbin)] + [str(b) for b in self.bintimes]
        args += _flag("maskpoints", self.maskpoints)
        return args


class alarm(VartoolsCommand):
    """Tamuz, Mazeh & North (2006) alarm statistic.

    A time-correlated-variability indicator sensitive to runs of
    same-sign residuals about the weighted mean — useful for flagging
    systematics and coherent variability.

    Parameters
    ----------
    maskpoints : str, optional
        Name of a mask variable; points with ``maskvar ≤ 0`` are
        excluded.

    See Also
    --------
    CLI command: ``-alarm``.
    Citation: Tamuz, Mazeh & North 2006 (MNRAS 367, 1521).
    """

    _vt_name = "alarm"

    def __init__(self, maskpoints: Optional[str] = None) -> None:
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        return ["-alarm"] + _flag("maskpoints", self.maskpoints)


class vonNeumann(VartoolsCommand):
    """von Neumann (1941) ratio eta = delta^2 / s^2 as a variability indicator.

    For uncorrelated Gaussian noise E[eta] = 2 (variance ~ 4/N).  Smoothly
    varying (correlated) signals drive eta well below 2; anti-correlated
    signals push eta above 2.  Useful as a sparse-/uneven-sampling variability
    metric (Sokolovsky et al. 2017, MNRAS 464, 274).  The light curve is
    time-sorted automatically before the calculation.

    Parameters
    ----------
    weighted : bool
        If True, use inverse-variance weighting.  The weighted ratio is
        ``eta_w = (2N/(N-1)) * sum w_pair_i (y[i+1]-y[i])^2 / sum w_i (y_i - <y>_w)^2``
        with ``w_i = 1/sigma_i^2`` and ``w_pair_i = 1/(sigma_i^2 + sigma_{i+1}^2)``.
        The (2N/(N-1)) prefactor restores E[eta_w] ~ 2 for white noise under
        any sigma distribution (the raw ratio-of-weighted-averages instead
        converges to <w>/<w_pair>, which equals 2 only for homoscedastic
        errors).  Reduces exactly to the unweighted form when sigma is
        constant.  Points with NaN / non-positive uncertainty are dropped.
    maskpoints : str, optional
        Name of a vector; only points where ``maskvar > 0`` are included.

    See Also
    --------
    CLI command: ``-vonNeumann``.
    Cite von Neumann, J. 1941, Annals of Mathematical Statistics, 12, 367
    if you use this command; for astronomical applications see Sokolovsky,
    K. V., et al. 2017, MNRAS, 464, 274.
    """

    _vt_name = "vonNeumann"

    def __init__(self, weighted: bool = False,
                 maskpoints: Optional[str] = None) -> None:
        self.weighted = weighted
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        # Strict-parser order: weighted, then maskpoints.
        args = ["-vonNeumann"]
        args += _bool("weighted", self.weighted)
        args += _flag("maskpoints", self.maskpoints)
        return args


class percentileratios(VartoolsCommand):
    """Robust scatter statistics from the magnitude distribution.

    For each pair of percentiles ``(p, q)`` with ``0 < p < q < 100`` the
    command emits two statistics per light curve:

        amp_p_q  = pct(q) - pct(p)
        asym_p_q = (pct(q) - median) / (median - pct(p))

    plus one additional statistic that does not depend on the pair list:

        medmeddev_over_stddev = median(|x - median(x)|) / stddev(x)

    For any symmetric distribution ``asym -> 1.0``; positively-skewed
    distributions produce ``asym > 1`` and negatively-skewed distributions
    produce ``asym < 1``.  For independent Gaussian noise
    ``medmeddev / stddev -> 0.6745`` in the large-N limit.

    Percentile interpolation matches the ``-stats`` command, so values are
    directly comparable to the corresponding ``pct(p)`` columns from
    ``stats``.  NaN magnitudes are dropped before any statistic is
    computed; light curves with fewer than two finite magnitudes, and
    ratios with a zero denominator (e.g. ``median == pct(p)``, or
    ``stddev == 0``), produce NaN outputs.

    Parameters
    ----------
    percentilepairs : sequence of (float, float), optional
        Percentile pairs to use.  Defaults to ``[(5, 95), (1, 99)]``.
        Each pair must satisfy ``0 < p, q < 100`` and ``p != q``; pairs
        with ``p > q`` are silently canonicalized to ``p < q``; duplicate
        pairs (after canonicalization) are rejected.  Floating-point
        percentiles are accepted (e.g. ``(2.5, 97.5)``).
    maskpoints : str, optional
        Name of a light-curve vector; points with ``maskvar > 0`` are
        included in the calculation, others are excluded (applied alongside
        NaN rejection).

    See Also
    --------
    CLI command: ``-percentileratios``.
    """

    _vt_name = "percentileratios"

    _DEFAULT_PAIRS = ((5.0, 95.0), (1.0, 99.0))

    def __init__(
        self,
        percentilepairs: Optional[Sequence[Sequence[float]]] = None,
        maskpoints: Optional[str] = None,
    ) -> None:
        self.maskpoints = maskpoints
        if percentilepairs is None:
            self.percentilepairs = list(self._DEFAULT_PAIRS)
            return

        canonical = []
        seen = set()
        for idx, pair in enumerate(percentilepairs):
            try:
                p, q = pair
            except (TypeError, ValueError):
                raise ValueError(
                    f"percentilepairs[{idx}] must be a (p, q) pair, got {pair!r}"
                )
            p = float(p)
            q = float(q)
            if not (0.0 < p < 100.0) or not (0.0 < q < 100.0):
                raise ValueError(
                    f"percentilepairs[{idx}] = ({p}, {q}) is outside the open interval (0, 100)"
                )
            if p == q:
                raise ValueError(
                    f"percentilepairs[{idx}] has p == q ({p})"
                )
            if p > q:
                p, q = q, p
            key = (p, q)
            if key in seen:
                raise ValueError(
                    f"percentilepairs has duplicate pair {p}:{q} (after canonicalizing p < q)"
                )
            seen.add(key)
            canonical.append((p, q))
        self.percentilepairs = canonical

    def _to_cli_args(self) -> List[str]:
        args = ["-percentileratios"]
        if list(map(tuple, self.percentilepairs)) != list(self._DEFAULT_PAIRS):
            pair_str = ",".join(f"{p}:{q}" for p, q in self.percentilepairs)
            args += ["percentilepairs", pair_str]
        args += _flag("maskpoints", self.maskpoints)
        return args


class beyondNsigma(VartoolsCommand):
    """Fraction of magnitudes beyond N*sigma above/below the median.

    For each light curve and each N value in ``Nvalues``, two fractions are
    emitted::

        frac_above_N = #{ x : x > median + N*sigma } / N_rej
        frac_below_N = #{ x : x < median - N*sigma } / N_rej

    where ``N_rej`` is the number of finite magnitudes after NaN rejection.
    Comparisons are strict (``>`` and ``<``).

    By default ``sigma`` is the sample standard deviation.  When
    ``useMAD=True`` is given, ``sigma`` is taken to be
    ``1.483 * median(|x - median(x)|)`` instead -- the Gaussian-consistent
    calibration of the MAD.  The MAD-based scale is robust to heavy tails
    or outliers, which inflate the stddev and widen the threshold; using
    MAD recovers a tighter threshold and correctly flags outliers as such.

    NaN magnitudes are dropped before any statistic is computed; light
    curves with fewer than two finite magnitudes produce NaN outputs.  When
    sigma is exactly zero (e.g. every magnitude equals the median) the
    fractions are reported as zero, since no point strictly exceeds a zero
    threshold.

    The ``N=1`` instance of this statistic corresponds to the
    ``Beyond1Std`` feature of Nun et al. 2015, "FATS: Feature Analysis for
    Time Series" (arXiv:1506.00010), generalized here to an arbitrary list
    of N values and to a choice of stddev or MAD scale.

    Parameters
    ----------
    Nvalues : sequence of float, optional
        N values to evaluate.  Defaults to ``[1.0, 3.0, 5.0]``.  Each value
        must be strictly positive, and duplicates are rejected.
        Floating-point values are accepted (e.g. ``[1.5, 2.5, 4.0]``).
    useMAD : bool, optional
        If True, use ``1.483 * MAD`` as the scale instead of the sample
        standard deviation.  Defaults to False.
    maskpoints : str, optional
        Name of a light-curve vector; points with ``maskvar > 0`` are
        included in the calculation, others are excluded (applied alongside
        NaN rejection).

    See Also
    --------
    CLI command: ``-beyondNsigma``.
    Citation: Nun et al. 2015 (arXiv:1506.00010) for the ``Beyond1Std``
    feature this generalizes.
    """

    _vt_name = "beyondNsigma"

    _DEFAULT_NVALUES = (1.0, 3.0, 5.0)

    def __init__(
        self,
        Nvalues: Optional[Sequence[float]] = None,
        useMAD: bool = False,
        maskpoints: Optional[str] = None,
    ) -> None:
        if Nvalues is None:
            self.Nvalues = list(self._DEFAULT_NVALUES)
        else:
            canonical = []
            seen = set()
            for idx, val in enumerate(Nvalues):
                try:
                    v = float(val)
                except (TypeError, ValueError):
                    raise ValueError(
                        f"Nvalues[{idx}] must be a number, got {val!r}"
                    )
                if v <= 0.0:
                    raise ValueError(
                        f"Nvalues[{idx}] = {v} must be > 0"
                    )
                if v in seen:
                    raise ValueError(
                        f"Nvalues has duplicate value {v}"
                    )
                seen.add(v)
                canonical.append(v)
            self.Nvalues = canonical
        self.useMAD = bool(useMAD)
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        args = ["-beyondNsigma"]
        if [float(v) for v in self.Nvalues] != list(self._DEFAULT_NVALUES):
            args += ["Nvalues", ",".join(str(v) for v in self.Nvalues)]
        if self.useMAD:
            args += ["useMAD"]
        args += _flag("maskpoints", self.maskpoints)
        return args


class slopestats(VartoolsCommand):
    """Per-pair slope statistics over consecutive light-curve points.

    After time-sorting (and optionally binning) the filtered light curve,
    consecutive points are paired and slopes are formed::

        s_i = (m[i+1] - m[i]) / (t[i+1] - t[i])

    Five statistics are reported per (LC, binsize) combination::

        median_abs_dmdt = median(|s_i|)
        max_abs_dmdt    = max(|s_i|)
        mad_dmdt        = 1.483 * median(|s_i - median(s_i)|)
        frac_above_T    = #{ s_i : s_i - <s> > T*sigma } / N_pairs
        frac_below_T    = #{ s_i : s_i - <s> < -T*sigma } / N_pairs

    ``<s>`` is the mean of the slopes; deviations are mean-centered
    regardless of the choice of sigma flavor.  Comparisons against
    ``T*sigma`` are strict (``>`` and ``<``).

    By default ``sigma`` is the sample standard deviation of the slopes
    (normalized by ``N_pairs - 1``).  When ``useMAD=True``, ``sigma`` is
    ``1.483 * median(|s_i - median(s_i)|)`` instead (= ``mad_dmdt``),
    which is more robust to heavy tails or outliers.

    The optional ``bintime`` parameter takes a list of bin sizes in
    MINUTES (assuming the light curve time axis is in days; the kernel
    divides each bintime by 1440 before binning).  Each binsize produces
    its own column set.  The optional ``binshift`` parameter (only valid
    in combination with ``bintime``) shifts the first bin by a fraction
    of the binwidth: ``t0_bin = t[0] - binshift * binsize``.  Canonical
    use is ``0 <= binshift < 1``.

    The optional ``maxgap`` parameter (in days) drops consecutive pairs
    whose time separation exceeds ``maxgap``, suppressing spurious large
    slopes that span long observational gaps.

    The ``max_abs_dmdt`` statistic corresponds to the ``MaxSlope``
    feature of Richards et al. 2011, ApJ 733, 10.  The
    ``frac_above_T*sigma`` and ``frac_below_T*sigma`` statistics are
    slope-domain generalizations of the magnitude-domain ``Beyond1Std``
    feature of the same paper, extended to a user-supplied list of T
    values and split into signed above/below counts.  The ``mad_dmdt``
    statistic is the slope-domain analog of the magnitude-domain median
    absolute deviation (here with the Gaussian-calibrated 1.483 factor).

    Parameters
    ----------
    bintime : sequence of float, optional
        List of bin sizes in minutes.  Each value must be strictly
        positive; duplicates are rejected.  If ``None`` (the default),
        no binning is performed.
    binshift : float, optional
        Fraction of the binwidth by which to shift the first bin.
        Canonical use is ``0 <= binshift < 1``.  Only valid in
        combination with ``bintime``.
    threshold : sequence of float, optional
        List of T values for the threshold-fraction statistics.
        Defaults to ``[3.0]``.  Each value must be strictly positive;
        duplicates are rejected.
    maxgap : float, optional
        Drop consecutive pairs whose time separation exceeds
        ``maxgap``, in the same time units as the light curve (days).
        Must be strictly positive.
    useMAD : bool, optional
        If True, use ``1.483 * MAD`` of the slopes as ``sigma`` instead
        of the sample standard deviation.  Defaults to False.
    maskpoints : str, optional
        Name of a light-curve vector; points with ``maskvar > 0`` are
        included in the calculation, others are excluded (applied
        alongside NaN rejection on the magnitudes).

    See Also
    --------
    CLI command: ``-slopestats``.
    Citation: Richards et al. 2011 (ApJ 733, 10) for the ``MaxSlope``
    feature this generalizes.
    """

    _vt_name = "slopestats"

    _DEFAULT_THRESHOLD = (3.0,)

    def __init__(
        self,
        bintime: Optional[Sequence[float]] = None,
        binshift: Optional[float] = None,
        threshold: Optional[Sequence[float]] = None,
        maxgap: Optional[float] = None,
        useMAD: bool = False,
        maskpoints: Optional[str] = None,
    ) -> None:
        self.bintime = self._validate_positive_unique(bintime, "bintime")

        if binshift is None:
            self.binshift = None
        else:
            try:
                self.binshift = float(binshift)
            except (TypeError, ValueError):
                raise ValueError(
                    f"binshift must be a number, got {binshift!r}"
                )

        if self.binshift is not None and self.bintime is None:
            raise ValueError(
                "binshift requires bintime to be specified"
            )

        if threshold is None:
            self.threshold = list(self._DEFAULT_THRESHOLD)
        else:
            self.threshold = self._validate_positive_unique(
                threshold, "threshold"
            )

        if maxgap is None:
            self.maxgap = None
        else:
            try:
                m = float(maxgap)
            except (TypeError, ValueError):
                raise ValueError(
                    f"maxgap must be a number, got {maxgap!r}"
                )
            if m <= 0.0:
                raise ValueError(f"maxgap = {m} must be > 0")
            self.maxgap = m

        self.useMAD = bool(useMAD)
        self.maskpoints = maskpoints

    @staticmethod
    def _validate_positive_unique(values, name):
        if values is None:
            return None
        canonical = []
        seen = set()
        for idx, val in enumerate(values):
            try:
                v = float(val)
            except (TypeError, ValueError):
                raise ValueError(
                    f"{name}[{idx}] must be a number, got {val!r}"
                )
            if v <= 0.0:
                raise ValueError(f"{name}[{idx}] = {v} must be > 0")
            if v in seen:
                raise ValueError(f"{name} has duplicate value {v}")
            seen.add(v)
            canonical.append(v)
        return canonical

    def _to_cli_args(self) -> List[str]:
        args = ["-slopestats"]
        if self.bintime is not None:
            args += ["bintime", ",".join(str(v) for v in self.bintime)]
        if self.binshift is not None:
            args += ["binshift", str(self.binshift)]
        if [float(v) for v in self.threshold] != list(self._DEFAULT_THRESHOLD):
            args += ["threshold",
                     ",".join(str(v) for v in self.threshold)]
        if self.maxgap is not None:
            args += ["maxgap", str(self.maxgap)]
        if self.useMAD:
            args += ["useMAD"]
        args += _flag("maskpoints", self.maskpoints)
        return args


class CodyM(VartoolsCommand):
    """Flux-asymmetry statistic M of Cody et al. 2014, AJ, 147, 82 (Eq. 7).

    ``M`` quantifies whether a light curve is asymmetric with respect to
    reflection about the median magnitude::

        M = (<d_10%> - d_med) / sigma_d

    where ``d_med`` is the median of the long-term-detrended,
    outlier-filtered light curve, ``sigma_d`` is its standard deviation,
    and ``<d_10%>`` is the mean of the combined faintest-decile and
    brightest-decile values.  For magnitude-valued light curves ``M``
    moves in the positive direction for dipping signatures (asymmetric
    toward faint excursions) and in the negative direction for bursting
    signatures (asymmetric toward bright excursions); a value close to
    zero indicates a symmetric magnitude distribution.  The sign
    convention is reversed for flux input.

    The kernel filters NaN magnitudes and (when ``maskpoints`` is given)
    points with ``maskvar <= 0``, subtracts a boxcar-mean smooth of full
    width ``trendwindow``, identifies sigma-clip outliers, then evaluates
    ``M`` on the surviving values using deciles of ``round(0.1 * N')``
    points each.  Two outlier-identification schemes are supported:

    * **Two-stage (paper-faithful):** supply ``outlierwindow`` to build a
      short-timescale residual on which the sigma clip operates.  Real
      variability is removed in that residual, so the genuine dips and
      bursts that ``M`` quantifies are not themselves clipped.
    * **Single-stage:** omit ``outlierwindow``; the sigma clip operates
      directly on the trend-detrended curve.  Simpler, but a deep
      dip/burst can clip itself and bias ``M`` toward zero.

    Parameters
    ----------
    trendwindow : float or str
        Boxcar full width for the long-term detrend, in the same units
        as the light-curve time axis.  Required.  Accepts a number, a
        bare variable name (vartools ``var``), or an explicit
        ``"var NAME"`` / ``"expr EXPR"`` string for per-LC sourcing.
    outlierwindow : float or str, optional
        Short-timescale boxcar full width that enables the two-stage
        outlier-rejection scheme.  Same var/expr accepting form as
        ``trendwindow``.  Default: omitted (single-stage).
    sigclip : float or str, optional
        Outlier-rejection threshold in sigma.  Default ``5.0``; ``0``
        disables outlier rejection.  Accepts var/expr.
    maskpoints : str, optional
        Name of a light-curve vector; only points with ``maskvar > 0``
        contribute.

    See Also
    --------
    CLI command: ``-CodyM``.
    Citation: Cody et al. 2014 (AJ 147, 82).
    """

    _vt_name = "CodyM"

    _DEFAULT_SIGCLIP = 5.0

    def __init__(
        self,
        trendwindow,
        outlierwindow=None,
        sigclip=5.0,
        maskpoints: Optional[str] = None,
    ) -> None:
        self.trendwindow = trendwindow
        self.outlierwindow = outlierwindow
        self.sigclip = sigclip
        self.maskpoints = maskpoints

        # Parse-time validation only when a literal numeric value was
        # given.  var/expr-sourced values are validated per-LC at runtime
        # by the kernel.
        if isinstance(trendwindow, (int, float)) and trendwindow <= 0.0:
            raise ValueError(f"trendwindow = {trendwindow} must be > 0")
        if (outlierwindow is not None
                and isinstance(outlierwindow, (int, float))
                and outlierwindow <= 0.0):
            raise ValueError(f"outlierwindow = {outlierwindow} must be > 0")
        if isinstance(sigclip, (int, float)) and sigclip < 0.0:
            raise ValueError(
                f"sigclip = {sigclip} must be >= 0 (0 disables rejection)"
            )

    def _to_cli_args(self) -> List[str]:
        args = ["-CodyM", "trendwindow"] + _varexpr(self.trendwindow)
        if self.outlierwindow is not None:
            args += ["outlierwindow"] + _varexpr(self.outlierwindow)
        # Emit sigclip only when non-default, so the CLI surface stays
        # tidy for the common case.
        if not (isinstance(self.sigclip, (int, float))
                and float(self.sigclip) == self._DEFAULT_SIGCLIP):
            args += ["sigclip"] + _varexpr(self.sigclip)
        args += _flag("maskpoints", self.maskpoints)
        return args


class CodyQ(VartoolsCommand):
    """Quasi-periodicity statistic Q of Cody et al. 2014, AJ, 147, 82 (Eq. 6).

    ``Q`` measures how close the light curve is to the photometric noise
    floor before and after a phase-folded periodic model is subtracted::

        Q = (rms_resid^2 - sigma^2) / (rms_raw^2 - sigma^2)

    where ``rms_raw`` is the standard deviation of the long-term-detrended
    light curve, ``rms_resid`` is the standard deviation of the same
    curve after subtraction of a boxcar-smoothed phase model at the
    supplied period, and ``sigma^2`` is the mean of the per-point
    squared errors.  ``Q`` approaches 0 for a strictly periodic light
    curve (the phase model captures essentially all the variance) and
    approaches 1 for one with no detectable periodicity (the phase
    model removes nothing; the denominator can collapse to a
    non-positive value in this regime, yielding NaN); intermediate
    values indicate quasi-periodic variability.

    The kernel filters NaN magnitudes, ``sig <= 0`` points, and (when
    ``maskpoints`` is given) ``maskvar <= 0`` points; subtracts a
    boxcar-mean smooth of full width ``trendwindow``; phase-folds to
    ``period``; and builds a circular boxcar smooth in the phase domain
    with full width ``phasesmooth`` (a fraction of the period).  The
    phase smoother is invariant under a global phase shift, so the
    folding epoch is fixed at ``t[0]`` and is not exposed.

    Parameters
    ----------
    period : float or str
        Period source.  Accepts a number, a back-reference keyword
        (``"ls"``, ``"aov"``, ``"pdm"``, ``"ftp"``, ``"bls"``,
        ``"injectharm"``) to the most-recent corresponding command,
        ``"fix P"``, ``"fixcolumn NAME"``, ``"list ['column' N]"``, a
        bare variable name (vartools ``var``), or an explicit
        ``"var NAME"`` / ``"expr EXPR"`` string.
    trendwindow : float or str
        Boxcar full width for the long-term detrend, in the same units
        as the light-curve time axis.  Required.  Accepts var/expr.
    phasesmooth : float or str, optional
        Phase-domain boxcar full width as a fraction of the period.
        Default ``0.25`` (the paper value).  Must lie in ``(0, 1]``.
        Accepts var/expr.
    maskpoints : str, optional
        Name of a light-curve vector; only points with ``maskvar > 0``
        contribute.

    See Also
    --------
    CLI command: ``-CodyQ``.
    Citation: Cody et al. 2014 (AJ 147, 82).
    """

    _vt_name = "CodyQ"

    _DEFAULT_PHASESMOOTH = 0.25

    def __init__(
        self,
        period,
        trendwindow,
        phasesmooth=0.25,
        maskpoints: Optional[str] = None,
    ) -> None:
        self.period = period
        self.trendwindow = trendwindow
        self.phasesmooth = phasesmooth
        self.maskpoints = maskpoints

        # Parse-time validation only for literal numeric values; var/expr
        # values are validated per-LC at runtime.
        if isinstance(trendwindow, (int, float)) and trendwindow <= 0.0:
            raise ValueError(f"trendwindow = {trendwindow} must be > 0")
        if isinstance(phasesmooth, (int, float)) and (
                phasesmooth <= 0.0 or phasesmooth > 1.0):
            raise ValueError(
                f"phasesmooth = {phasesmooth} must be in (0, 1]"
            )
        if isinstance(period, (int, float)) and period <= 0.0:
            raise ValueError(f"period = {period} must be > 0")

    def _to_cli_args(self) -> List[str]:
        args = ["-CodyQ"] + _period_spec(self.period)
        args += ["trendwindow"] + _varexpr(self.trendwindow)
        if not (isinstance(self.phasesmooth, (int, float))
                and float(self.phasesmooth) == self._DEFAULT_PHASESMOOTH):
            args += ["phasesmooth"] + _varexpr(self.phasesmooth)
        args += _flag("maskpoints", self.maskpoints)
        return args

    def _resolve_back_references(self, prev) -> None:
        from ._helpers import _resolve_period_backref
        self.period = _resolve_period_backref(prev, self.period)


class structurefunction(VartoolsCommand):
    """Ensemble structure function with optional damped-random-walk fit.

    Compute the per-light-curve ensemble structure function (SF) on lag
    bins constructed from all pairs of observations.  Two estimators are
    supported (selected by ``estimator``):

    - ``"squared"`` (default) -- the squared-difference form with
      bin-averaged noise subtraction::

          SF^2(dt) = <(m_j - m_i)^2>_pairs in bin
                   - <sigma_i^2 + sigma_j^2>_pairs in bin

      This is the ``D^(1)(tau)`` of Simonetti, Cordes & Heeschen 1985,
      ApJ, 296, 46.  Bins where the noise-subtracted variance is
      non-positive are NaN.

    - ``"mad"`` -- the absolute-deviation form with per-pair noise
      subtraction (Hughes, Aller & Aller 1992; Schmidt et al. 2010
      Equation 2)::

          V(dt) = <(sqrt(pi)/2) * |m_j - m_i|
                  - sqrt(sigma_i^2 + sigma_j^2)>_pairs in bin

      More robust to outliers, biased low when the intrinsic variance
      approaches the photometric noise floor; bins whose averaged
      value goes non-positive are NaN.

    Per-bin error bars are reported as ``sigma_SF = SF / sqrt(2 N_eff)``
    (squared) or ``sigma_SF = SF / sqrt(N_eff / (pi/2 - 1))`` (mad),
    where ``N_eff = min(N_pairs, N_obs/2)`` bounds the number of
    effectively independent pairs.  This recovers the right
    ``~SF / sqrt(N_obs)`` precision floor at long lags rather than
    continuing to shrink as ``1 / sqrt(N_pairs)``.

    When ``fitDRW=True``, a damped-random-walk (Ornstein-Uhlenbeck /
    CAR(1)) model (Kelly, Bechtold & Siemiginowska 2009) is fit to the
    well-determined SF bins via a Nelder-Mead simplex on
    ``(log SF_inf, log tau)``, minimising
    ``chi^2 = sum [(SF_obs - SF_model) / sigma_SF]^2``.  The DRW SF is
    ``SF_DRW(dt) = SF_inf * sqrt(1 - exp(-dt / tau))`` with
    ``SF_inf = sqrt(2) * sigma_long``; the reported scalar amplitude is
    ``sigma_long`` (MacLeod et al. 2010 convention, mag).

    Parameters
    ----------
    bins : {"log", "linear", "edges"}
        Bin geometry selector.  ``"log"`` and ``"linear"`` require
        ``Nbins``; ``"edges"`` requires ``edges``.
    Nbins : int, optional
        Number of bins for ``"log"`` / ``"linear"``.  Must be >= 2.
    edges : sequence of float, optional
        Strictly increasing list of positive bin edges for ``"edges"``
        mode.  ``len(edges) - 1`` bins are produced.
    lagrange : tuple (lagmin, lagmax), optional
        Explicit lag-range edges for ``"log"`` / ``"linear"`` binning.
        Each element may be a number, a bare variable name (vartools
        ``var``), or an explicit ``"var NAME"`` / ``"expr EXPR"``
        string.  If omitted (default), ``lagmin`` is the smallest
        consecutive time spacing in the (filtered) light curve and
        ``lagmax`` is the full baseline.  Ignored for ``"edges"``.
    estimator : {"squared", "mad"}
        SF estimator family.  Default ``"squared"``.
    fitDRW : bool
        Fit a DRW model to the SF bins.  When ``True``, emits five
        scalar columns: ``STRUCTUREFUNCTION_SIGMA_N`` (= sigma_long,
        mag), ``_TAU_N``, ``_CHI2_N``, ``_DOF_N``, ``_CONVERGED_N``.
    sigma0 : float or str, optional
        Initial guess for the DRW ``sigma_long`` (only valid with
        ``fitDRW=True``).  Accepts var/expr.
    tau0 : float or str, optional
        Initial guess for the DRW damping time-scale.  Accepts var/expr.
    reportsfvalsintable : sequence of float, optional
        Strictly increasing list of positive lag values.  For each
        requested lag ``e_k`` the wrapper emits four scalar columns
        per light curve: ``STRUCTUREFUNCTION_DT_k_N`` (the actual bin
        centre containing ``e_k``), ``_SF_k_N``, ``_SIGMA_SF_k_N``,
        ``_NPAIRS_k_N``.  Out-of-range or noise-dominated bins yield
        NaN / 0.
    save_result : bool, str, or Output
        Control the ``.sf`` aux file (one row per bin: ``dt_center``,
        ``SF``, ``sigma_SF``, ``n_pairs``).

        - ``False`` (default) -- no aux file.
        - ``True`` -- write to a pipeline-managed temp dir and capture
          into ``result.files["structurefunction_result_N"]``.
        - path string -- write to that directory, no capture.
        - ``Output(path, capture=True)`` -- both.
    maskpoints : str, optional
        Mask variable; points with ``maskvar > 0`` contribute.

    See Also
    --------
    CLI command: ``-structurefunction``.
    Citations
        ``estimator="squared"``: Simonetti, Cordes & Heeschen 1985,
        ApJ, 296, 46.
        ``estimator="mad"``: Hughes, Aller & Aller 1992, ApJ, 396, 469;
        Schmidt et al. 2010, ApJ, 714, 1194 (Equation 2).
        ``fitDRW``: Kelly, Bechtold & Siemiginowska 2009, ApJ, 698, 895;
        MacLeod et al. 2010, ApJ, 721, 1014.
    """

    _vt_name = "structurefunction"
    _DEFAULT_ESTIMATOR = "squared"

    def __init__(
        self,
        bins: str,
        Nbins: Optional[int] = None,
        edges: Optional[Sequence[float]] = None,
        lagrange=None,
        estimator: str = "squared",
        fitDRW: bool = False,
        sigma0=None,
        tau0=None,
        reportsfvalsintable: Optional[Sequence[float]] = None,
        save_result=False,
        maskpoints: Optional[str] = None,
    ) -> None:
        # Bin mode + parameter consistency
        if bins not in ("log", "linear", "edges"):
            raise ValueError(
                f"bins must be 'log', 'linear', or 'edges'; got {bins!r}"
            )

        edges_list = None
        if bins in ("log", "linear"):
            if Nbins is None:
                raise ValueError(f"bins={bins!r} requires Nbins")
            if not isinstance(Nbins, int) or isinstance(Nbins, bool) \
                    or Nbins < 2:
                raise ValueError(
                    f"Nbins = {Nbins!r} must be an integer >= 2"
                )
            if edges is not None:
                raise ValueError(
                    "edges is only valid for bins='edges'"
                )
        else:  # edges
            if edges is None:
                raise ValueError("bins='edges' requires edges")
            edges_list = self._validate_edges(edges, "edges")
            if len(edges_list) < 3:
                raise ValueError(
                    "bins='edges' requires at least 3 edges (giving "
                    ">= 2 bins)"
                )
            if Nbins is not None:
                raise ValueError(
                    "Nbins is only valid for bins='log' or 'linear'"
                )

        if estimator not in ("squared", "mad"):
            raise ValueError(
                f"estimator must be 'squared' or 'mad'; got {estimator!r}"
            )

        if (sigma0 is not None or tau0 is not None) and not fitDRW:
            raise ValueError(
                "sigma0 / tau0 are only valid when fitDRW=True"
            )
        if isinstance(sigma0, (int, float)) and not isinstance(sigma0, bool) \
                and sigma0 <= 0.0:
            raise ValueError(f"sigma0 = {sigma0} must be > 0")
        if isinstance(tau0, (int, float)) and not isinstance(tau0, bool) \
                and tau0 <= 0.0:
            raise ValueError(f"tau0 = {tau0} must be > 0")

        if lagrange is not None:
            if len(lagrange) != 2:
                raise ValueError(
                    "lagrange must be a 2-tuple (lagmin, lagmax)"
                )
            lmin, lmax = lagrange
            if isinstance(lmin, (int, float)) and not isinstance(lmin, bool) \
                    and lmin <= 0:
                raise ValueError(f"lagrange lagmin = {lmin} must be > 0")
            if (isinstance(lmin, (int, float)) and not isinstance(lmin, bool)
                    and isinstance(lmax, (int, float))
                    and not isinstance(lmax, bool)
                    and lmax <= lmin):
                raise ValueError(
                    f"lagrange lagmax ({lmax}) must be > lagmin ({lmin})"
                )

        report_list = None
        if reportsfvalsintable is not None:
            report_list = self._validate_edges(
                reportsfvalsintable, "reportsfvalsintable"
            )

        self.bins = bins
        self.Nbins = Nbins
        self.edges = list(edges_list) if edges_list is not None else None
        self.lagrange = tuple(lagrange) if lagrange is not None else None
        self.estimator = estimator
        self.fitDRW = bool(fitDRW)
        self.sigma0 = sigma0
        self.tau0 = tau0
        self.reportsfvalsintable = report_list
        self.save_result = save_result
        self.maskpoints = maskpoints

    @staticmethod
    def _validate_edges(values, name):
        canonical = []
        prev = None
        for idx, val in enumerate(values):
            try:
                v = float(val)
            except (TypeError, ValueError):
                raise ValueError(
                    f"{name}[{idx}] must be a number, got {val!r}"
                )
            if v <= 0.0:
                raise ValueError(f"{name}[{idx}] = {v} must be > 0")
            if prev is not None and v <= prev:
                raise ValueError(
                    f"{name} must be strictly increasing "
                    f"(got {v} after {prev})"
                )
            canonical.append(v)
            prev = v
        return canonical

    def _to_cli_args(self) -> List[str]:
        args = ["-structurefunction", "bins"]
        if self.bins == "edges":
            args += ["edges", ",".join(repr(v) for v in self.edges)]
        else:
            args += [self.bins, str(self.Nbins)]

        if self.lagrange is not None:
            args += ["lagrange"]
            args += _varexpr(self.lagrange[0])
            args += _varexpr(self.lagrange[1])

        if self.estimator != self._DEFAULT_ESTIMATOR:
            args += ["estimator", self.estimator]

        if self.fitDRW:
            args += ["fitDRW"]
            if self.sigma0 is not None:
                args += ["sigma0"] + _varexpr(self.sigma0)
            if self.tau0 is not None:
                args += ["tau0"] + _varexpr(self.tau0)

        if self.reportsfvalsintable is not None:
            args += [
                "reportsfvalsintable",
                ",".join(repr(v) for v in self.reportsfvalsintable),
            ]

        if _should_emit(self.save_result):
            outdir = getattr(self, "_outdir", ".")
            spec = _norm_save(self.save_result)
            actual_outdir = spec.path if spec.path is not None else outdir
            args += ["save", actual_outdir]

        args += _flag("maskpoints", self.maskpoints)
        return args

    def _output_file_specs(self):
        if not _should_emit(self.save_result):
            return {}
        # Logical name must match the ``save_result`` attribute so the
        # pipeline collector (which reads ``save_<logical_name>``) actually
        # captures the file into ``result.files["structurefunction_result_N"]``.
        return {"result": (".sf", None)}


class drwfit(VartoolsCommand):
    """Direct maximum-likelihood damped-random-walk fit to a light curve.

    Fit a damped-random-walk (DRW) model -- the continuous-time
    first-order autoregressive CAR(1) / Ornstein-Uhlenbeck process,
    equivalently a Gaussian process with a Matern-1/2 (exponential)
    covariance kernel -- directly to a light curve by maximum likelihood.
    The likelihood is evaluated with the Kelly, Bechtold & Siemiginowska
    2009 (ApJ, 698, 895) state-space recursion (their Equations 6-13), so
    each evaluation is ``O(N)`` with no matrix inversion.  A downhill
    simplex (Nelder-Mead) minimises ``-2 ln L`` over
    ``(log sigma_long, log tau)`` (or jointly with ``mu`` when fitting the
    mean).

    This direct-likelihood method recovers ``tau`` substantially more
    accurately than fitting a DRW to the structure function
    (:class:`structurefunction` with ``fitDRW=True``), especially on short
    baselines (MacLeod et al. 2010 Section 4.2; Kelly 2009 Section 3.1).
    The reported amplitude ``DRWFIT_SIGMA_N`` is ``sigma_long`` (the
    long-term magnitude standard deviation, MacLeod 2010 convention, mag),
    matching what :class:`structurefunction` emits so the two commands are
    directly comparable on the same light curve.

    Scalar output columns (``N`` = command position): ``DRWFIT_SIGMA_N``
    (sigma_long, mag), ``_TAU_N`` (damping time-scale), ``_MU_N`` (fitted
    long-term mean; NaN when ``mean="subtract"``), ``_LNL_N`` (best-fit
    ln L), ``_DLNL_NOISE_N`` (ln L_best - ln L at sigma_long -> 0),
    ``_DLNL_INF_N`` (ln L_best - ln L at tau -> infinity), and
    ``_CONVERGED_N`` (1 if the amoeba converged).  The two Delta-ln-L
    likelihood-ratio indicators measure detection against the
    pure-noise and the tau -> infinity limits (see MacLeod 2010
    Section 3.1).

    Parameters
    ----------
    mean : {"fit", "subtract", "fix"}, optional
        Long-term-mean handling.  ``None`` (default) leaves it to the
        vartools default (``"fit"``).

        - ``"fit"`` -- fit ``mu`` jointly with ``sigma_long`` and ``tau``
          (3-parameter amoeba).
        - ``"subtract"`` -- subtract the weighted mean before fitting;
          the reported ``DRWFIT_MU_N`` is meaningless (NaN).
        - ``"fix"`` -- hold ``mu`` at ``mean_value`` (2-parameter amoeba).
    mean_value : float or str, optional
        The fixed ``mu`` for ``mean="fix"``.  Required for, and only valid
        with, ``mean="fix"``.  Accepts a number, a bare variable name
        (vartools ``var``), or an explicit ``"var NAME"`` / ``"expr EXPR"``
        string for per-light-curve sourcing.
    sigma0 : float or str, optional
        Initial guess for the amoeba's ``sigma_long``.  Accepts var/expr.
        Default: derived from the light-curve variance.
    tau0 : float or str, optional
        Initial guess for the amoeba's ``tau``.  Accepts var/expr.
        Default: ``0.1 * T_baseline``.
    mean0 : float or str, optional
        Initial guess for the amoeba's ``mu`` (only meaningful when the
        mean is fit).  Accepts var/expr.  Default: the weighted mean.
    save_result : bool, str, or Output
        Control the ``.drwfit`` aux file (one row per original light-curve
        point, NaN rows preserved for filtered-out points; eight columns
        ``t x sig_meas x_hat_fwd Omega_fwd chi_fwd x_smoothed
        Omega_smoothed`` from the forward Kalman filter and the RTS
        backward smoother).

        - ``False`` (default) -- no aux file.
        - ``True`` -- write to a pipeline-managed temp dir and capture
          into ``result.files["drwfit_result_N"]``.
        - path string -- write to that directory, no capture.
        - ``Output(path, capture=True)`` -- both.
    correctlc : {"smoothed", "forecast"}, optional
        Replace the in-memory light curve with the DRW residuals, in
        place, so downstream commands operate on the corrected curve.

        - ``"smoothed"`` -- subtract the RTS-smoothed model (uses past and
          future points; whitens the curve toward the noise floor).
        - ``"forecast"`` -- subtract the one-step-ahead Kalman forecast
          (uses only past points; retains unexplained variability).
    modelvar : sequence (mode, varname), optional
        Store the DRW model into a new light-curve variable ``varname``
        (usable by downstream commands such as ``-stats`` or ``-expr``)
        without modifying the light curve.  ``mode`` is ``"smoothed"`` or
        ``"forecast"`` with the same meaning as ``correctlc``.
    maskpoints : str, optional
        Name of a light-curve vector; only points with ``maskvar > 0``
        contribute to the fit.

    See Also
    --------
    structurefunction : structure-function-based DRW fit (less accurate
        ``tau`` recovery; useful as an independent cross-check).
    CLI command: ``-drwfit``.
    Citations
        Kelly, Bechtold & Siemiginowska 2009, ApJ, 698, 895;
        MacLeod et al. 2010, ApJ, 721, 1014.
    """

    _vt_name = "drwfit"

    _MEAN_MODES = ("fit", "subtract", "fix")
    _MODEL_MODES = ("smoothed", "forecast")

    def __init__(
        self,
        mean: Optional[str] = None,
        mean_value=None,
        sigma0=None,
        tau0=None,
        mean0=None,
        save_result=False,
        correctlc: Optional[str] = None,
        modelvar=None,
        maskpoints: Optional[str] = None,
    ) -> None:
        if mean is not None and mean not in self._MEAN_MODES:
            raise ValueError(
                f"mean must be one of {self._MEAN_MODES}; got {mean!r}"
            )
        if mean == "fix":
            if mean_value is None:
                raise ValueError("mean='fix' requires mean_value")
        elif mean_value is not None:
            raise ValueError(
                "mean_value is only valid when mean='fix'"
            )

        if isinstance(sigma0, (int, float)) and not isinstance(sigma0, bool) \
                and sigma0 <= 0.0:
            raise ValueError(f"sigma0 = {sigma0} must be > 0")
        if isinstance(tau0, (int, float)) and not isinstance(tau0, bool) \
                and tau0 <= 0.0:
            raise ValueError(f"tau0 = {tau0} must be > 0")

        if correctlc is not None and correctlc not in self._MODEL_MODES:
            raise ValueError(
                f"correctlc must be one of {self._MODEL_MODES}; "
                f"got {correctlc!r}"
            )

        modelvar_norm = None
        if modelvar is not None:
            if isinstance(modelvar, str) or len(modelvar) != 2:
                raise ValueError(
                    "modelvar must be a (mode, varname) pair, e.g. "
                    "('smoothed', 'drwmod')"
                )
            mv_mode, mv_name = modelvar
            if mv_mode not in self._MODEL_MODES:
                raise ValueError(
                    f"modelvar mode must be one of {self._MODEL_MODES}; "
                    f"got {mv_mode!r}"
                )
            if not isinstance(mv_name, str) or not mv_name:
                raise ValueError("modelvar varname must be a non-empty string")
            modelvar_norm = (mv_mode, mv_name)

        self.mean = mean
        self.mean_value = mean_value
        self.sigma0 = sigma0
        self.tau0 = tau0
        self.mean0 = mean0
        self.save_result = save_result
        self.correctlc = correctlc
        self.modelvar = modelvar_norm
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        args = ["-drwfit"]

        if self.mean is not None:
            args += ["mean", self.mean]
            if self.mean == "fix":
                args += _varexpr(self.mean_value)

        if self.sigma0 is not None:
            args += ["sigma0"] + _varexpr(self.sigma0)
        if self.tau0 is not None:
            args += ["tau0"] + _varexpr(self.tau0)
        if self.mean0 is not None:
            args += ["mean0"] + _varexpr(self.mean0)

        if _should_emit(self.save_result):
            outdir = getattr(self, "_outdir", ".")
            spec = _norm_save(self.save_result)
            actual_outdir = spec.path if spec.path is not None else outdir
            args += ["save", actual_outdir]

        if self.correctlc is not None:
            args += ["correctlc", self.correctlc]

        if self.modelvar is not None:
            args += ["modelvar", self.modelvar[0], self.modelvar[1]]

        args += _flag("maskpoints", self.maskpoints)
        return args

    def _output_file_specs(self):
        if not _should_emit(self.save_result):
            return {}
        return {"result": (".drwfit", None)}


class rescalesig(VartoolsCommand):
    """Rescale per-point uncertainties so that χ²/dof = 1.

    Multiply every uncertainty by a single factor chosen so the reduced
    chi-square about the weighted mean equals 1, and report that factor.

    Parameters
    ----------
    maskpoints : str, optional
        Name of a mask variable; points with ``maskvar ≤ 0`` are
        excluded.

    See Also
    --------
    CLI command: ``-rescalesig``.
    """

    _vt_name = "rescalesig"

    def __init__(self, maskpoints: Optional[str] = None) -> None:
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        return ["-rescalesig"] + _flag("maskpoints", self.maskpoints)


class ensemblerescalesig(VartoolsCommand):
    """Rescale uncertainties using the ensemble mag–χ² relation.

    Fit the magnitude-dependent scatter across the ensemble of light
    curves and rescale each light curve's uncertainties accordingly
    (sigma-clipped fit).

    Parameters
    ----------
    sigclip : float
        Sigma-clipping threshold for the ensemble fit.  Default 5.0.
    maskpoints : str, optional
        Name of a mask variable; points with ``maskvar ≤ 0`` are
        excluded.

    See Also
    --------
    CLI command: ``-ensemblerescalesig``.
    """

    _vt_name = "ensemblerescalesig"

    def __init__(
        self,
        sigclip: float = 5.0,
        maskpoints: Optional[str] = None,
    ) -> None:
        self.sigclip = sigclip
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        return ["-ensemblerescalesig", str(self.sigclip)] + _flag("maskpoints", self.maskpoints)


class stats(VartoolsCommand):
    """Compute statistics for one or more light curve variables.

    Parameters
    ----------
    variables : str or list of str
        Variable name(s) to compute statistics for (comma-separated in vartools).
    statistics : str or list of str
        Statistic name(s) to compute.  Each name may appear at most once
        per call (vartools rejects repeats; e.g.
        ``["mean","stddev","stddev"]`` errors out at parse time).

        Recognised names (full list, matching the vartools ``-stats``
        parser in ``src/statistics.c``):

        =================  ==================================================
        Name               Statistic
        =================  ==================================================
        ``mean``           Arithmetic mean
        ``weightedmean``   Inverse-variance-weighted mean (uses err)
        ``median``         Median (50th percentile)
        ``wmedian``        Weighted median (uses err)
        ``stddev``         Sample standard deviation (1/(N-1) normalisation)
        ``meddev``         Median absolute deviation from the **mean**
        ``medmeddev``      Median absolute deviation from the **median**
        ``MAD``            ``meddev``-style MAD scaled by 1.4826
                           (normal-consistent estimator of σ)
        ``kurtosis``       Sample kurtosis (Fisher normalisation, biased)
        ``skewness``       Sample skewness (biased)
        ``max``            Maximum value
        ``min``            Minimum value
        ``sum``            Sum
        ``pctXX``          XXth percentile (e.g. ``pct10``, ``pct90``,
                           ``pct50.5``)
        ``wpctXX``         Inverse-variance-weighted XXth percentile
        =================  ==================================================

        See the vartools ``-stats`` CLI docs for full definitions and
        the per-LC output-column names.
    maskpoints : str, optional
        Variable name; only points with ``maskvar > 0`` contribute.

    Examples
    --------
    ::

        stats("mag", "mean,median,stddev")
        stats(["mag", "err"], ["mean", "stddev"])
        stats("mag", ["pct10", "pct50", "pct90"])

    See Also
    --------
    CLI command: ``-stats``.
    """

    _vt_name = "stats"

    def __init__(
        self,
        variables: Union[str, List[str]],
        statistics: Union[str, List[str]],
        maskpoints: Optional[str] = None,
    ) -> None:
        self.variables = variables if isinstance(variables, str) else ",".join(variables)
        self.statistics = statistics if isinstance(statistics, str) else ",".join(statistics)
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        return (["-stats", self.variables, self.statistics]
                + _flag("maskpoints", self.maskpoints))


class harmonicfilter(VartoolsCommand):
    """Fit (and optionally subtract) a truncated Fourier series at one or
    more known periods — the ``-harmonicfilter`` vartools command.

    Output columns are emitted under the ``HarmonicFilter_*`` prefix when
    invoked via this class.  The legacy :class:`Killharm` subclass below
    invokes the same command under the ``-Killharm`` synonym and produces
    ``Killharm_*`` columns for backward compatibility.

    Parameters
    ----------
    period : float or str
        Period of the signal to fit (days).  Forms:

        - A number — fit at that fixed period.
        - ``"ls"`` / ``"aov"`` / ``"injectharm"`` — back-reference to the
          best period from a prior command of that type in the same
          chain.
        - ``"both"`` — fit at *both* the most-recent LS and AOV periods
          (a two-period fit).
        - ``"list"`` — read the period from a list-file column.
        - A bare identifier or expression string — per-LC variable or
          analytic expression.

        ``"bls"`` is **not** accepted (the vartools ``-harmonicfilter`` /
        ``-Killharm`` parser rejects it); use a numeric period or one of
        the keywords above.  Default ``"ls"``.
    nharm : int
        Number of harmonics to fit.  Default 3.
    nsubharm : int
        Number of sub-harmonics to fit.  Default 0.
    save_model : bool, str, or Output
        Fitted harmonic model.  ``True`` captures as
        ``result.files["harmonicfilter_model_N"]``; a path string writes
        to that directory without capturing; ``Output(path, capture=True)``
        does both.
    fitonly : bool
        Fit the model but do not subtract it from the light curve.
    output_format : str, optional
        Coefficient output format: ``"outampphase"``, ``"outampradphase"``,
        ``"outRphi"``, or ``"outRradphi"``.
    clip : float, optional
        Sigma-clipping threshold applied after fitting before refitting.
    maskpoints : str, optional
        Name of a mask variable; points with ``maskvar ≤ 0`` are
        excluded.

    See Also
    --------
    CLI command: ``-harmonicfilter`` (the :class:`Killharm` subclass
    emits the ``-Killharm`` synonym instead).
    """

    _vt_name = "harmonicfilter"
    # CLI token emitted by this class — subclasses override to swap to a
    # synonym (see Killharm).
    _cli_token = "-harmonicfilter"

    def __init__(
        self,
        period="ls",
        nharm: int = 3,
        nsubharm: int = 0,
        save_model: bool = False,
        fitonly: bool = False,
        output_format: Optional[str] = None,
        clip: Optional[float] = None,
        maskpoints: Optional[str] = None,
    ) -> None:
        self.period = period
        self.nharm = nharm
        self.nsubharm = nsubharm
        self.save_model = save_model
        self.fitonly = fitonly
        self.output_format = output_format
        self.clip = clip
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = [self._cli_token] + self._killharm_period_spec()
        args += [str(self.nharm), str(self.nsubharm)]
        args += _outtoken(self.save_model, outdir)
        args += _bool("fitonly", self.fitonly)
        if self.output_format is not None:
            args += [self.output_format]
        if self.clip is not None:
            args += ["clip"] + _varexpr(self.clip)
        args += _flag("maskpoints", self.maskpoints)
        return args

    def _killharm_period_spec(self) -> List[str]:
        """Build period spec tokens for Killharm.

        Killharm's "fix" spec has the form: fix Nper per1 ... perN, where
        each perN may be a number, ``var <name>``, or ``expr <expression>``.
        Keywords like "ls", "aov", "both", "injectharm", "list" are passed
        through verbatim.  A tuple/list of values emits the multi-period
        ``fix Nper <v1> ... <vN>`` form (used when a chained ``period="both"``
        back-reference resolves to a pair of periods).
        """
        # vartools' valid -harmonicfilter / -Killharm period keywords are
        # aov / ls / both / injectharm / list — *not* bls (the CLI
        # rejects "bls" with "undefined variable" because it routes the
        # token through the var-name path).  "both" is the vartools
        # CLI keyword that uses a paired (LS, AOV) period set.
        _KILLHARM_KEYWORDS = {"ls", "aov", "both", "injectharm", "list"}
        p = self.period
        if isinstance(p, (int, float)):
            return ["fix", "1", str(p)]
        if isinstance(p, (tuple, list)):
            return ["fix", str(len(p))] + [str(v) for v in p]
        s = str(p)
        tokens = s.split()
        first = tokens[0] if tokens else ""
        if first == "fix":
            # "fix <period> [period ...]" — insert the period count.  Each
            # period is assumed to be a bare value or pre-qualified token.
            periods = tokens[1:]
            return ["fix", str(len(periods))] + periods
        if first in _KILLHARM_KEYWORDS:
            return tokens
        # "var NAME" / "expr EXPR" pre-qualified forms (produced by the
        # per-LC substitution in Pipeline.run_batch) need the "fix 1"
        # prefix prepended, not another "expr"/"var" wrap.
        if first in ("var", "expr"):
            return ["fix", "1"] + tokens
        # Single period specified as a variable name (from -inlistvars or
        # -expr listvar) or an analytic expression.  Must be wrapped in
        # the "fix 1 var|expr" form since Killharm requires an explicit
        # period count.
        if re.match(r'^[A-Za-z_]\w*$', s):
            return ["fix", "1", "var", s]
        return ["fix", "1", "expr", s]

    def _resolve_back_references(self, prev) -> None:
        from ._helpers import (_resolve_period_backref, _most_recent_lookup,
                                _coerce_to_numeric)
        from pyvartools.perlc import PerLC
        # "both" is a Killharm-only back-ref that pulls the most-recent LS
        # period and the most-recent AOV period.  The result becomes a
        # (ls_period, aov_period) tuple that _killharm_period_spec renders
        # into `fix 2 <ls> <aov>`.
        if isinstance(self.period, str) and self.period.strip() == "both":
            ls_stats = _most_recent_lookup(prev, ["LS"])
            aov_stats = _most_recent_lookup(prev, ["aov", "aov_harm"])
            if ls_stats is None:
                raise LookupError(
                    "Back-reference 'both' has no prior -LS command in "
                    "this chain"
                )
            if aov_stats is None:
                raise LookupError(
                    "Back-reference 'both' has no prior -aov or -aov_harm "
                    "command in this chain"
                )
            ls_per = _coerce_to_numeric(ls_stats.Period_1)
            aov_per = _coerce_to_numeric(aov_stats.Period_1)
            # Both scalars → 2-tuple of floats.  Any PerLC input → 2-tuple of
            # PerLCs (multi-period batch emission handled elsewhere; for now
            # raise if either value is PerLC, since Killharm's per-LC pair
            # emission isn't plumbed through).
            if isinstance(ls_per, PerLC) or isinstance(aov_per, PerLC):
                raise NotImplementedError(
                    "Killharm(period='both') across a batch chain boundary "
                    "is not supported (per-LC pairs of periods would need "
                    "two -inlistvars columns).  Pass the two periods "
                    "explicitly via Raw() or a single-LC chain."
                )
            self.period = (float(ls_per), float(aov_per))
            return
        self.period = _resolve_period_backref(prev, self.period)

    def _output_file_specs(self):
        # Suffix follows the invoking CLI token — subclasses override.
        suffix = (".harmonicfilter.model"
                  if self._cli_token == "-harmonicfilter"
                  else ".killharm.model")
        return {"model": (suffix, None)}


class Killharm(harmonicfilter):
    """Legacy name for :class:`harmonicfilter`.  Accepted for backward
    compatibility; emits ``-Killharm`` on the command line and produces
    output columns under the ``Killharm_*`` prefix.  New code should use
    :class:`harmonicfilter`.

    See Also
    --------
    CLI command: ``-Killharm``.  Same parameters as
    :class:`harmonicfilter`.
    """

    _vt_name = "Killharm"
    _cli_token = "-Killharm"


class linfit(VartoolsCommand):
    """Fit a linear combination of analytic functions.

    Parameters
    ----------
    function : str
        Analytic function string (vartools expression syntax).
    paramlist : str
        Comma-separated list of free parameter names and initial values.
    modelvar : str, optional
        Variable name to store the best-fit model.
    reject : float, optional
        Sigma-clipping threshold for outlier rejection after fitting.
    reject_usemad : bool
        Use MAD instead of standard deviation for rejection scatter estimate.
    reject_iter : bool
        Iteratively reject outliers.
    reject_fixednum : int, optional
        Maximum number of rejection iterations (requires ``reject_iter=True``).
    correct_lc : bool
        Subtract the best-fit model from the light curve.
    save_model : bool, str, or Output
        Fitted model curve.  ``True`` captures as
        ``result.files["linfit_model_N"]``; a path string writes to that
        directory without capturing; ``Output(path, capture=True)`` does
        both.
    model_nameformat : str, optional
        Format string for the output model filename.
    fitmask : str, optional
        Name of a mask variable; points where the variable is non-zero are
        excluded from the fit.

    See Also
    --------
    CLI command: ``-linfit``.
    """

    _vt_name = "linfit"

    def __init__(
        self,
        function: str,
        paramlist: str,
        modelvar: Optional[str] = None,
        reject: Optional[float] = None,
        reject_usemad: bool = False,
        reject_iter: bool = False,
        reject_fixednum: Optional[int] = None,
        correct_lc: bool = False,
        save_model: bool = False,
        model_nameformat: Optional[str] = None,
        fitmask: Optional[str] = None,
    ) -> None:
        self.function = function
        self.paramlist = paramlist
        self.modelvar = modelvar
        self.reject = reject
        self.reject_usemad = reject_usemad
        self.reject_iter = reject_iter
        self.reject_fixednum = reject_fixednum
        self.correct_lc = correct_lc
        self.save_model = save_model
        self.model_nameformat = model_nameformat
        self.fitmask = fitmask

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-linfit", self.function, self.paramlist]
        if self.modelvar is not None:
            args += ["modelvar", self.modelvar]
        if self.reject is not None:
            args += ["reject"] + _varexpr(self.reject)
            args += _bool("useMAD", self.reject_usemad)
            if self.reject_iter:
                args += ["iter"]
                if self.reject_fixednum is not None:
                    args += ["fixednum", str(self.reject_fixednum)]
        args += _bool("correctlc", self.correct_lc)
        m_spec = _norm_save(self.save_model)
        if _should_emit(m_spec):
            args += ["omodel", m_spec.path or outdir]
            if self.model_nameformat is not None:
                args += ["format", self.model_nameformat]
        args += _flag("fitmask", self.fitmask)
        return args

    def _output_file_specs(self):
        return {"model": (".linfit.model", None)}


class Injectharm(VartoolsCommand):
    """Inject a harmonic (sinusoidal) signal into the light curve.

    Parameters
    ----------
    period : float or str
        Period of the signal to inject.  Float → emits ``"fix <val>"``.
        String passthrough (e.g. ``"logrand 0.2 2.0"``, ``"list"``).
    amplitude : float or str
        Fundamental-harmonic amplitude.  Float → ``"ampfix"``.  ``"rand"``
        → ``"amprand"``; ``"list"`` → ``"amplist"``; bare identifier →
        ``"ampvar <name>"``; any other string → ``"ampexpr <string>"``.
    nharm : int
        Total number of harmonics, including the fundamental.  ``nharm=1``
        is a pure sinusoid.  ``nharm=11`` is fundamental + 10 overtones.
    phase : float or str, optional
        Fundamental-harmonic phase.  Float → ``"phasefix"``.  ``"rand"``
        → ``"phaserand"``; bare identifier → ``"phasevar <name>"``; other
        string → ``"phaseexpr <string>"``.  Default 0.0.
    harmonic_amps_rel : sequence of float, optional
        Relative amplitudes for harmonics 2..``nharm`` (length ``nharm-1``).
        Each entry emits ``ampfix R amprel`` so the harmonic's amplitude is
        ``R`` times the fundamental amplitude.  Use this with the R values
        from a Fourier-series fit to inject a realistic templated signal.
    harmonic_phases_rel : sequence of float, optional
        Relative phases for harmonics 2..``nharm`` (length ``nharm-1``).
        Each entry emits ``phasefix phi phaserel`` so the harmonic phase
        is offset from the fundamental phase by ``phi``.
    nsubharm : int
        Number of sub-harmonics (default 0).  When non-zero, sub-harmonics
        share the fundamental's amp/phase mode (no per-sub-harmonic list).
    save_model : bool, str, or Output
        Injected-signal model.  ``True`` captures as
        ``result.files["Injectharm_model_N"]``; a path string writes to
        that directory without capturing; ``Output(path, capture=True)``
        does both.

    See Also
    --------
    CLI command: ``-Injectharm``.

    Examples
    --------
    Pure sinusoid::

        cmd.Injectharm(period=2.5, amplitude=0.05)

    10-harmonic RR Lyrae template, fundamental amplitude 0.1, random phase::

        cmd.Injectharm(period=0.514333, amplitude=0.1, phase="rand", nharm=11,
                       harmonic_amps_rel=[0.47, 0.36, 0.24, 0.16, 0.11, 0.06,
                                          0.04, 0.03, 0.02, 0.02],
                       harmonic_phases_rel=[0.61, 0.26, -0.07, 0.61, 0.29,
                                            0.22, 0.95, 0.59, 0.66, 0.94])
    """

    _vt_name = "Injectharm"

    # Special amplitude / phase keyword strings that map to bare CLI flags
    # rather than to ``ampexpr`` / ``phaseexpr``.  Both the short Pythonic
    # form (``"rand"``) and the verbatim CLI form (``"amprand"``) are
    # accepted to match the spelling in vartools' own help text.
    _AMP_KEYWORDS = {
        "rand": "amprand",        "amprand": "amprand",
        "list": "amplist",        "amplist": "amplist",
        "logrand": "amplogrand",  "amplogrand": "amplogrand",
    }
    _PHASE_KEYWORDS = {
        "rand": "phaserand",      "phaserand": "phaserand",
        "list": "phaselist",      "phaselist": "phaselist",
    }

    def __init__(
        self,
        period,
        amplitude,
        nharm: int = 1,
        phase=0.0,
        nsubharm: int = 0,
        save_model: bool = False,
        harmonic_amps_rel: Optional[Sequence[float]] = None,
        harmonic_phases_rel: Optional[Sequence[float]] = None,
    ) -> None:
        self.period = period
        self.amplitude = amplitude
        self.nharm = nharm
        self.phase = phase
        self.nsubharm = nsubharm
        self.save_model = save_model
        self.harmonic_amps_rel = (
            list(harmonic_amps_rel) if harmonic_amps_rel is not None else None
        )
        self.harmonic_phases_rel = (
            list(harmonic_phases_rel) if harmonic_phases_rel is not None else None
        )

        # Validate paired lengths up front so the error surfaces at
        # construction time rather than at run time.
        if self.harmonic_amps_rel is not None:
            if len(self.harmonic_amps_rel) != nharm - 1:
                raise ValueError(
                    f"harmonic_amps_rel must have length nharm-1 "
                    f"({nharm - 1}); got {len(self.harmonic_amps_rel)}.  "
                    f"Provide one relative amplitude per harmonic 2..{nharm}."
                )
        if self.harmonic_phases_rel is not None:
            if len(self.harmonic_phases_rel) != nharm - 1:
                raise ValueError(
                    f"harmonic_phases_rel must have length nharm-1 "
                    f"({nharm - 1}); got {len(self.harmonic_phases_rel)}."
                )
        # If only one of the two lists was supplied, the user almost
        # certainly wants both.  Reject early with a clear message.
        if (self.harmonic_amps_rel is None) != (self.harmonic_phases_rel is None):
            raise ValueError(
                "harmonic_amps_rel and harmonic_phases_rel must be supplied "
                "together (or both omitted)."
            )

    def _amp_tokens(self, amp) -> List[str]:
        """Build the ampfix/ampvar/ampexpr/amp<keyword> tokens for one amp value.

        The string form is whitespace-split so range-form keywords like
        ``"rand 0.01 0.1"`` map to ``["amprand", "0.01", "0.1"]`` rather
        than being treated as a single bareword.
        """
        if isinstance(amp, (int, float)):
            return ["ampfix", str(amp)]
        if isinstance(amp, str):
            tokens = amp.split()
            head = tokens[0] if tokens else ""
            if head in self._AMP_KEYWORDS:
                return [self._AMP_KEYWORDS[head]] + tokens[1:]
            if re.match(r'^[A-Za-z_]\w*$', amp):
                return ["ampvar", amp]
            return ["ampexpr", amp]
        return ["ampfix", str(amp)]

    def _phase_tokens(self, phase) -> List[str]:
        """Build the phasefix/phasevar/phaseexpr/phase<keyword> tokens for one phase value."""
        if isinstance(phase, (int, float)):
            return ["phasefix", str(phase)]
        if isinstance(phase, str):
            tokens = phase.split()
            head = tokens[0] if tokens else ""
            if head in self._PHASE_KEYWORDS:
                return [self._PHASE_KEYWORDS[head]] + tokens[1:]
            if re.match(r'^[A-Za-z_]\w*$', phase):
                return ["phasevar", phase]
            return ["phaseexpr", phase]
        return ["phasefix", str(phase)]

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-Injectharm"] + _period_spec(self.period)
        # vartools Nharm = nharm - 1: Nharm=0 means just the fundamental,
        # Nharm=1 means fundamental + 1st harmonic, etc.
        vt_nharm = self.nharm - 1
        args += [str(vt_nharm)]

        # Fundamental harmonic spec.
        args += self._amp_tokens(self.amplitude)
        args += self._phase_tokens(self.phase)

        # Higher harmonics.  When per-harmonic relative lists are supplied,
        # emit ``ampfix R amprel phasefix phi phaserel`` for each.
        # Otherwise, replicate the fundamental's amp/phase spec for every
        # harmonic — preserves prior behaviour where amplitude/phase are
        # shared across all harmonics.
        if self.harmonic_amps_rel is not None:
            for R, phi in zip(self.harmonic_amps_rel, self.harmonic_phases_rel):
                args += ["ampfix", str(R), "amprel",
                         "phasefix", str(phi), "phaserel"]
        else:
            for _ in range(self.nharm - 1):
                args += self._amp_tokens(self.amplitude)
                args += self._phase_tokens(self.phase)

        # Sub-harmonics — vartools expects Nsubharm followed by Nsubharm
        # copies of the (amp, phase) spec.
        args += [str(self.nsubharm)]
        for _ in range(self.nsubharm):
            args += self._amp_tokens(self.amplitude)
            args += self._phase_tokens(self.phase)
        args += _outtoken(self.save_model, outdir)
        return args

    def _output_file_specs(self):
        return {"model": (".injectharm.model", None)}


class Injecttransit(VartoolsCommand):
    """Inject a transit signal into the light curve.

    Parameters
    ----------
    period : float or str
        Orbital period.  Float → emits ``"Pfix <val>"``.  String is passed
        through as-is (e.g. ``"Plogrand 0.2 2.0"``, ``"Plist"``).
    Rp : float or str
        Planet-to-star radius ratio.  Float → ``"Rpfix <val>"``.  String
        passthrough (e.g. ``"Rplogrand 0.05 0.15"``).
    Mp : float or str
        Planet mass (solar masses).
    phase : float or str
        Orbital phase of transit centre (0–1).  String e.g. ``"phaserand"``.
    sini : float or str
        Sine of orbital inclination.  String e.g. ``"sinirand"``.
    Mstar : float or str
        Stellar mass (solar masses).
    Rstar : float or str
        Stellar radius (solar radii).
    e : float or str
        Eccentricity.  Used with ``eomega`` mode (default).
    omega : float or str
        Argument of periastron.  Used with ``eomega`` mode.
    hk : bool
        When ``True``, use the ``hk`` eccentricity parameterisation
        (emits ``"hk" h k`` instead of ``"eomega" e omega``).
    h : float or str
        ``h = e sin(omega)`` (used when ``hk=True``).
    k : float or str
        ``k = e cos(omega)`` (used when ``hk=True``).
    dilute : float or str, optional
        Dilution factor.  Float → ``["dilute", "fix", str(val)]``.
        String passthrough (e.g. ``"fix 0.5"`` or ``"list"``).
    ld_type : str
        Limb-darkening law: ``"quad"`` (default) or ``"nonlin"``.
    ld_coeffs : list of float, optional
        Limb-darkening coefficients.  Default ``[0.3, 0.3]`` (quadratic).
        Two coefficients for ``"quad"``, four for ``"nonlin"``.
    save_model : bool, str, or Output
        Injected-transit model.  ``True`` captures as
        ``result.files["Injecttransit_model_N"]``; a path string writes
        to that directory without capturing; ``Output(path, capture=True)``
        does both.

    See Also
    --------
    CLI command: ``-Injecttransit``.
    """

    _vt_name = "Injecttransit"

    def __init__(
        self,
        period,
        Rp,
        Mp,
        phase,
        sini,
        Mstar,
        Rstar,
        e: Union[float, str] = 0.0,
        omega: Union[float, str] = 0.0,
        hk: bool = False,
        h: Union[float, str] = 0.0,
        k: Union[float, str] = 0.0,
        dilute: Union[float, str, None] = None,
        ld_type: str = "quad",
        ld_coeffs: Optional[List[float]] = None,
        save_model: bool = False,
    ) -> None:
        self.period = period
        self.Rp = Rp
        self.Mp = Mp
        self.phase = phase
        self.sini = sini
        self.Mstar = Mstar
        self.Rstar = Rstar
        self.e = e
        self.omega = omega
        self.hk = hk
        self.h = h
        self.k = k
        self.dilute = dilute
        self.ld_type = ld_type
        self.ld_coeffs = ld_coeffs or [0.3, 0.3]
        self.save_model = save_model

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-Injecttransit"] + _injectparam("P", self.period)
        args += _injectparam("Rp", self.Rp)
        args += _injectparam("Mp", self.Mp)
        args += _injectparam("phase", self.phase)
        args += _injectparam("sini", self.sini)
        if self.hk:
            args += ["hk"] + _injectparam("h", self.h) + _injectparam("k", self.k)
        else:
            args += ["eomega"] + _injectparam("e", self.e) + _injectparam("o", self.omega)
        args += _injectparam("Mstar", self.Mstar)
        args += _injectparam("Rstar", self.Rstar)
        if self.dilute is not None:
            if isinstance(self.dilute, (int, float)):
                args += ["dilute", "fix", str(self.dilute)]
            else:
                args += ["dilute"] + str(self.dilute).split()
        args += [self.ld_type]
        ldc = list(self.ld_coeffs)
        if ldc and isinstance(ldc[0], str) and ldc[0] in (
                "ldlist", "ldfix", "ldvar", "ldexpr"):
            args += [str(c) for c in ldc]
        else:
            args += ["ldfix"] + [str(c) for c in ldc]
        args += _outtoken(self.save_model, outdir)
        return args

    def _output_file_specs(self):
        return {"model": (".injecttransit.model", None)}


class sortlc(VartoolsCommand):
    """Sort the light curve observations.

    Parameters
    ----------
    var : str, optional
        Variable to sort by.  Default is time (``"t"``).
    reverse : bool
        Sort in reverse (descending) order.

    See Also
    --------
    CLI command: ``-sortlc``.
    """

    _vt_name = "sortlc"

    def __init__(
        self,
        var: Optional[str] = None,
        reverse: bool = False,
    ) -> None:
        self.var = var
        self.reverse = reverse

    def _to_cli_args(self) -> List[str]:
        args = ["-sortlc"]
        args += _flag("var", self.var)
        args += _bool("reverse", self.reverse)
        return args


class restricttimes(VartoolsCommand):
    """Restrict the light curve to a time range or time list.

    Parameters
    ----------
    mode : str
        One of ``"JDrange"``, ``"JDrangebylc"``, ``"JDlist"``,
        ``"imagelist"``, ``"expr"``.
    minJD, maxJD : float or str, optional
        For ``"JDrange"`` / ``"JDrangebylc"`` modes.
    JDfilename : str, optional
        For ``"JDlist"`` mode.
    expression : str, optional
        For ``"expr"`` mode.
    exclude : bool
        Invert the restriction (exclude the specified range).
    markrestrict : str, optional
        Variable to record which points were restricted.
    noinitmark : bool
        Do not initialise the markrestrict variable.

    See Also
    --------
    CLI command: ``-restricttimes``.  Undo with :class:`restoretimes`.
    """

    _vt_name = "restricttimes"

    def __init__(
        self,
        mode: str = "JDrange",
        minJD=None,
        maxJD=None,
        JDfilename: Optional[str] = None,
        expression: Optional[str] = None,
        exclude: bool = False,
        markrestrict: Optional[str] = None,
        noinitmark: bool = False,
    ) -> None:
        self.mode = mode
        self.minJD = minJD
        self.maxJD = maxJD
        self.JDfilename = JDfilename
        self.expression = expression
        self.exclude = exclude
        self.markrestrict = markrestrict
        self.noinitmark = noinitmark

    def _to_cli_args(self) -> List[str]:
        args = ["-restricttimes"]
        args += _bool("exclude", self.exclude)
        if self.mode == "JDrange":
            args += ["JDrange", str(self.minJD), str(self.maxJD)]
        elif self.mode == "JDrangebylc":
            args += ["JDrangebylc"] + _pval(self.minJD, "fix") + _pval(self.maxJD, "fix")
        elif self.mode == "JDlist":
            args += ["JDlist", self.JDfilename]
        elif self.mode == "imagelist":
            args += ["imagelist", self.JDfilename]
        elif self.mode == "expr":
            args += ["expr", self.expression]
        if self.markrestrict is not None:
            args += ["markrestrict", self.markrestrict]
            args += _bool("noinitmark", self.noinitmark)
        return args


class restoretimes(VartoolsCommand):
    """Undo a prior -restricttimes command.

    Parameters
    ----------
    prior_command : int
        Index (1-based) of the -restricttimes command to undo.

    See Also
    --------
    CLI command: ``-restoretimes``.
    """

    _vt_name = "restoretimes"

    def __init__(self, prior_command: int = 1) -> None:
        self.prior_command = prior_command

    def _to_cli_args(self) -> List[str]:
        return ["-restoretimes", str(self.prior_command)]


class savelc(VartoolsCommand):
    """Save the current light curve state (for later restoration).

    Pipeline-stateful: meaningful only within a single ``Pipeline``
    invocation alongside a later :class:`restorelc`.

    See Also
    --------
    CLI command: ``-savelc``.  Restore with :class:`restorelc`.
    """

    _vt_name = "savelc"

    def _to_cli_args(self) -> List[str]:
        return ["-savelc"]


class restorelc(VartoolsCommand):
    """Restore a previously saved light curve state.

    Parameters
    ----------
    savenumber : int
        Index of the -savelc save point to restore (1-based).
    vars : str or list of str, optional
        Restore only specific variables instead of the full light curve.

    See Also
    --------
    CLI command: ``-restorelc``.  Save with :class:`savelc`.
    """

    _vt_name = "restorelc"

    def __init__(
        self,
        savenumber: int = 1,
        vars: Optional[Union[str, List[str]]] = None,
    ) -> None:
        self.savenumber = savenumber
        self.vars = vars

    def _to_cli_args(self) -> List[str]:
        args = ["-restorelc", str(self.savenumber)]
        if self.vars is not None:
            v = self.vars if isinstance(self.vars, str) else ",".join(self.vars)
            args += ["vars", v]
        return args


class difffluxtomag(VartoolsCommand):
    """Convert image-subtraction differential flux to magnitudes.

    Parameters
    ----------
    mag_constant : float or str
        Zero-point magnitude.  Accepts var/expr forms.
    offset : float or str
        Flux offset added before conversion.  Accepts var/expr forms.
        Default 0.0.
    magcolumn : int, optional
        Which magnitude column to convert (1-based).

    See Also
    --------
    CLI command: ``-difffluxtomag``.
    """

    _vt_name = "difffluxtomag"

    def __init__(
        self,
        mag_constant: float,
        offset: float = 0.0,
        magcolumn: Optional[int] = None,
    ) -> None:
        self.mag_constant = mag_constant
        self.offset = offset
        self.magcolumn = magcolumn

    def _to_cli_args(self) -> List[str]:
        args = ["-difffluxtomag"] + _varexpr(self.mag_constant) + _varexpr(self.offset)
        args += _flag("magcolumn", self.magcolumn)
        return args


class fluxtomag(VartoolsCommand):
    """Convert flux to magnitudes.

    Parameters
    ----------
    mag_constant : float or str
        Zero-point magnitude.  Accepts var/expr forms.
    offset : float or str
        Flux offset added before conversion.  Accepts var/expr forms.
        Default 0.0.

    See Also
    --------
    CLI command: ``-fluxtomag``.  Inverse: :class:`magtoflux`.
    """

    _vt_name = "fluxtomag"

    def __init__(self, mag_constant: float, offset: float = 0.0) -> None:
        self.mag_constant = mag_constant
        self.offset = offset

    def _to_cli_args(self) -> List[str]:
        return ["-fluxtomag"] + _varexpr(self.mag_constant) + _varexpr(self.offset)


class magtoflux(VartoolsCommand):
    """Convert magnitude to flux.  Inverse of :class:`fluxtomag`.

    Parameters
    ----------
    mag_constant : float or str, optional
        Zero-point magnitude.  Required unless ``normalize=True``.
        The conversion is ``flux = 10**((mag_constant - mag) / 2.5)``
        with ``sig_flux = flux * sig_mag / 1.0857``.
    normalize : bool, optional
        If True, compute fluxes with an arbitrary internal zero-point
        and then divide the flux and flux-uncertainty arrays by the
        median flux (NaNs rejected), so the output light curve has
        median flux 1.  Cannot be combined with ``mag_constant``.

    See Also
    --------
    CLI command: ``-magtoflux``.  Inverse: :class:`fluxtomag`.
    """

    _vt_name = "magtoflux"

    def __init__(
        self,
        mag_constant: Optional[float] = None,
        normalize: bool = False,
    ) -> None:
        if normalize and mag_constant is not None:
            raise ValueError(
                "magtoflux: cannot specify both mag_constant and normalize=True"
            )
        if not normalize and mag_constant is None:
            raise ValueError(
                "magtoflux: must specify mag_constant, or pass normalize=True"
            )
        self.mag_constant = mag_constant
        self.normalize = normalize

    def _to_cli_args(self) -> List[str]:
        if self.normalize:
            return ["-magtoflux", "normalize"]
        return ["-magtoflux"] + _varexpr(self.mag_constant)


class changeerror(VartoolsCommand):
    """Set the uncertainties equal to the RMS scatter of the magnitudes.

    Replace every per-point uncertainty with a single value: the RMS
    deviation of the (finite) magnitudes about their mean.  Useful when
    the input light curve has no usable error column.

    Parameters
    ----------
    maskpoints : str, optional
        Name of a mask variable; points with ``maskvar ≤ 0`` are
        excluded from the RMS calculation.

    See Also
    --------
    CLI command: ``-changeerror``.
    """

    _vt_name = "changeerror"

    def __init__(self, maskpoints: Optional[str] = None) -> None:
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        return ["-changeerror"] + _flag("maskpoints", self.maskpoints)


class changevariable(VartoolsCommand):
    """Copy a variable into one of the standard columns (t, mag, err, id).

    Parameters
    ----------
    column : str
        Which standard role to reassign: ``"t"``, ``"mag"``, ``"err"``,
        or ``"id"``.
    var : str
        Name of the existing variable that subsequent commands should
        use for that role.

    See Also
    --------
    CLI command: ``-changevariable``.
    """

    _vt_name = "changevariable"

    def __init__(self, column: str, var: str) -> None:
        self.column = column
        self.var = var

    def _to_cli_args(self) -> List[str]:
        return ["-changevariable", self.column, self.var]


class copylc(VartoolsCommand):
    """Replicate the light curve in memory for Monte-Carlo simulation.

    Make *ncopies* additional copies of the current light curve; each
    copy flows through the remaining pipeline independently and produces
    its own output row.  Commonly used ahead of ``-Injecttransit`` /
    ``-addnoise`` to build a per-LC simulation ensemble.

    Parameters
    ----------
    ncopies : int
        Number of replicated copies to create.

    See Also
    --------
    CLI command: ``-copylc``.  (Pipelines containing ``-copylc`` cannot
    use ``resume=True`` on streaming runs — one input row produces many
    output rows.)
    """

    _vt_name = "copylc"

    def __init__(self, ncopies: int) -> None:
        self.ncopies = ncopies

    def _to_cli_args(self) -> List[str]:
        return ["-copylc", str(self.ncopies)]


class medianfilter(VartoolsCommand):
    """Apply a moving median (or mean) filter to the light curve.

    By default the filtered series is subtracted from the magnitudes
    (a high-pass filter); with ``replace=True`` the magnitudes are
    replaced by the filtered series (a low-pass filter).

    Parameters
    ----------
    time : float or str
        Filter timescale (same units as the time column).  Accepts
        var/expr forms.
    method : str
        ``"median"`` (default), ``"average"``, or ``"weightedaverage"``.
    replace : bool
        Replace the magnitudes with the filtered series (low-pass)
        rather than subtracting it (high-pass).

    See Also
    --------
    CLI command: ``-medianfilter``.
    """

    _vt_name = "medianfilter"

    def __init__(
        self,
        time: float,
        method: str = "median",
        replace: bool = False,
    ) -> None:
        self.time = time
        self.method = method
        self.replace = replace

    def _to_cli_args(self) -> List[str]:
        args = ["-medianfilter"] + _varexpr(self.time)
        if self.method in ("average", "weightedaverage"):
            args.append(self.method)
        args += _bool("replace", self.replace)
        return args


class fourierfilter(VartoolsCommand):
    """Full-band Fourier-domain filter (``-fourierfilter``).

    Applies a band filter and/or an analytic ``filterexpr`` in
    frequency space via GSL's mixed-radix complex FFT (O(N log N))
    and requires a GSL-enabled vartools build.  Uniformly-sampled
    data is FFT-filtered directly; non-uniform data must be
    interpolated onto a uniform grid first using ``resample=<delta>``.
    This is distinct from :class:`harmonicfilter`, which fits
    harmonics of one or more *known* periods.

    Parameters
    ----------
    mode : str
        Filter type: ``"full"``, ``"highpass"``, ``"lowpass"``,
        ``"bandpass"``, or ``"bandcut"``.
    minfreq : float or str, optional
        Low-frequency cutoff.  Required for ``highpass``, ``bandpass``,
        ``bandcut``.  Accepts var/expr/fixcolumn forms.
    maxfreq : float or str, optional
        High-frequency cutoff.  Required for ``lowpass``, ``bandpass``,
        ``bandcut``.  Accepts var/expr/fixcolumn forms.
    filterexpr : str, optional
        Analytic filter applied to each Fourier coefficient as a
        function of frequency.  The frequency variable name defaults to
        ``"f"`` (e.g. ``"exp(-(f/0.5)**2)"``) and may be renamed via
        ``freqvar``.  The expression may only reference constants and
        per-star scalars, not light-curve vectors.
    freqvar : str, optional
        Override the default frequency-variable name (``"f"``) used in
        ``filterexpr``.
    fullspec : bool
        Compute Fourier coefficients across the full Nyquist range even
        when the selected band is narrower.  Useful with ``save_fouriercoeffs``.
    forcefft : bool
        Force the FFT fast path even when sampling is not detected as
        uniform (prints a warning if used on genuinely non-uniform
        data).  On uniformly-sampled data the FFT is used by default;
        ``forcefft=True`` is a no-op in that case.
    taper : str, optional
        Smooth-edge taper to apply at each cut edge instead of a
        brick-wall cutoff.  One of:

        - ``"linear"``: linear ramp (sharpest transition).
        - ``"cosine"`` (aliases ``"tukey"``, ``"hann"``):
          ``0.5*(1 - cos(pi*u))``.
        - ``"blackman"``: ``0.42 - 0.5*cos(pi*u) + 0.08*cos(2pi*u)``.
        - ``"kaiser"``: requires ``taper_beta``.

        The taper is centered on each cut edge and spans
        ``[edge - taper_deltafreq, edge + taper_deltafreq]``.  With
        ``mode="full"`` there are no edges and ``taper`` is a no-op
        (vartools prints a warning).
    taper_deltafreq : float, optional
        Half-width of the taper window, in frequency units.  Required
        when ``taper`` is given.
    taper_beta : float, optional
        Shape parameter for ``taper="kaiser"``; larger values give
        smoother but wider transitions (``beta ~= 5`` is Hann-like).
    resample : float, str, or None, optional
        Enable the resample → FFT → filter → IFFT → resample-back
        path (required for non-uniformly-sampled data).  Accepts:

        - ``"delmin"`` — use the minimum dt in the LC.
        - a positive ``float`` — fixed delta-t.
        - a string expression (e.g. ``"medt"``) — evaluated per LC.

        When absent and the input is non-uniform, vartools prints a
        warning to stderr and **skips the filter for that LC** (the
        mag column passes through unchanged; output columns are set
        to ``Mean_Mag=0``, ``RMS_Out=RMS_In``, ``Nfreqcalc=Nfreqfilt=0``).
        Subsequent LCs and subsequent pipeline commands are not
        affected.
    gapbreak_type : str, optional
        When ``resample`` is set, split the LC at any gap exceeding the
        specified threshold and filter each segment independently.
        Type is one of ``"fix"``, ``"expr"``, ``"frac_min_sep"``,
        ``"frac_med_sep"``, or ``"percentile_sep"``.
    gapbreak_value : float or str, optional
        Threshold value for the gap-break spec.  Units depend on
        ``gapbreak_type``: seconds (or whatever time unit the LC uses)
        for ``"fix"`` / ``"expr"``; a multiplier for
        ``"frac_min_sep"`` / ``"frac_med_sep"``; a percentile (0-100)
        for ``"percentile_sep"``.
    padmode : str, optional
        Edge-padding mode applied before the FFT (both the direct-FFT
        path and the resample path).  The FFT implicitly treats the
        signal as periodic, so if the first and last sample values
        disagree (common for astronomical LCs) the wrap-around injects
        spurious spectral power and the filtered output shows Gibbs-
        style ringing near the segment boundaries.  ``"wrap"``
        (default) uses the native FFT behaviour; ``"reflect"`` mirrors
        the signal at each edge; ``"zero"`` zero-extends (around the
        segment mean) at each edge.
    padfrac : float, optional
        Pad length per side, as a fraction of the segment length.
        Default ``0.5`` when ``padmode`` is ``"reflect"`` or
        ``"zero"`` (doubles the total FFT length); ignored for
        ``"wrap"``.  Clamped to at most one segment length for
        ``"reflect"``.
    save_fouriercoeffs : bool, str, or :class:`pyvartools.Output`, optional
        Write the Fourier cos/sin coefficients to a file.  Captures as
        ``result.files["fourierfilter_fouriercoeffs_N"]`` when truthy.
        See :doc:`Auxiliary output files <commands/index>`.
    nowarn : bool, optional
        When ``True``, suppress all per-LC runtime warnings from
        ``-fourierfilter`` (non-uniform advisory, within-segment gap
        vs. ``minfreq``, taper-edge overlap, forcefft on non-uniform,
        resample delta <= 0).  Useful in batch pipelines where the
        caller has vetted the data and doesn't need the repeated
        advisories.  Parse-time warnings about CLI misuse are not
        suppressed.

    See Also
    --------
    CLI command: ``-fourierfilter``.  For fitting harmonics of a
    *known* period instead, see :class:`harmonicfilter`.
    """

    _vt_name = "fourierfilter"

    _VALID_TAPERS = ("linear", "cosine", "tukey", "hann", "blackman", "kaiser")
    _VALID_GAPBREAK_TYPES = ("fix", "expr", "frac_min_sep",
                             "frac_med_sep", "percentile_sep")
    _VALID_PADMODES = ("wrap", "reflect", "zero")

    def __init__(
        self,
        mode: str = "full",
        minfreq: Union[float, str, None] = None,
        maxfreq: Union[float, str, None] = None,
        filterexpr: Optional[str] = None,
        freqvar: Optional[str] = None,
        fullspec: bool = False,
        forcefft: bool = False,
        taper: Optional[str] = None,
        taper_deltafreq: Optional[float] = None,
        taper_beta: Optional[float] = None,
        resample: Union[float, str, None] = None,
        gapbreak_type: Optional[str] = None,
        gapbreak_value: Union[float, str, None] = None,
        padmode: Optional[str] = None,
        padfrac: Optional[float] = None,
        nowarn: bool = False,
        save_fouriercoeffs=False,
    ) -> None:
        if mode not in ("full", "highpass", "lowpass", "bandpass", "bandcut"):
            raise ValueError(
                f"fourierfilter(): unknown mode {mode!r}; must be one of "
                f"full, highpass, lowpass, bandpass, bandcut"
            )
        if mode in ("highpass", "bandpass", "bandcut") and minfreq is None:
            raise ValueError(f"fourierfilter(mode={mode!r}): minfreq is required")
        if mode in ("lowpass", "bandpass", "bandcut") and maxfreq is None:
            raise ValueError(f"fourierfilter(mode={mode!r}): maxfreq is required")
        if freqvar is not None and filterexpr is None:
            raise ValueError(
                "fourierfilter(): freqvar has no effect without filterexpr"
            )
        if taper is not None:
            if taper not in self._VALID_TAPERS:
                raise ValueError(
                    f"fourierfilter(): unknown taper {taper!r}; must be one "
                    f"of {', '.join(self._VALID_TAPERS)}"
                )
            if taper_deltafreq is None or taper_deltafreq <= 0:
                raise ValueError(
                    "fourierfilter(): taper_deltafreq must be a positive "
                    "number when taper is given"
                )
            if taper == "kaiser" and (taper_beta is None or taper_beta <= 0):
                raise ValueError(
                    "fourierfilter(): taper='kaiser' requires a positive "
                    "taper_beta (try 5-8 for Hann-like behavior)"
                )
        elif taper_deltafreq is not None or taper_beta is not None:
            raise ValueError(
                "fourierfilter(): taper_deltafreq and taper_beta have no "
                "effect without taper"
            )
        # Gapbreak only makes sense alongside resample
        if gapbreak_type is not None:
            if resample is None:
                raise ValueError(
                    "fourierfilter(): gapbreak_type requires resample to "
                    "be set (gapbreak only applies on the resample path)"
                )
            if gapbreak_type not in self._VALID_GAPBREAK_TYPES:
                raise ValueError(
                    f"fourierfilter(): unknown gapbreak_type "
                    f"{gapbreak_type!r}; must be one of "
                    f"{', '.join(self._VALID_GAPBREAK_TYPES)}"
                )
            if gapbreak_value is None:
                raise ValueError(
                    "fourierfilter(): gapbreak_value is required when "
                    "gapbreak_type is given"
                )
        elif gapbreak_value is not None:
            raise ValueError(
                "fourierfilter(): gapbreak_value has no effect without "
                "gapbreak_type"
            )
        if padmode is not None:
            if padmode not in self._VALID_PADMODES:
                raise ValueError(
                    f"fourierfilter(): unknown padmode {padmode!r}; "
                    f"must be one of {', '.join(self._VALID_PADMODES)}"
                )
            if padfrac is not None and padfrac < 0:
                raise ValueError(
                    "fourierfilter(): padfrac must be >= 0"
                )
        elif padfrac is not None:
            raise ValueError(
                "fourierfilter(): padfrac has no effect without padmode"
            )
        # resample delmin|float|str validation
        if resample is not None:
            if isinstance(resample, str):
                if resample != "delmin" and not resample.replace(".", "", 1).lstrip("-").isdigit():
                    # treat it as an expression; pass through to CLI
                    pass
            elif isinstance(resample, (int, float)):
                if resample <= 0:
                    raise ValueError(
                        "fourierfilter(): resample delta must be positive"
                    )
            else:
                raise ValueError(
                    "fourierfilter(): resample must be a positive float, "
                    '"delmin", or a string expression'
                )
        self.mode = mode
        self.minfreq = minfreq
        self.maxfreq = maxfreq
        self.filterexpr = filterexpr
        self.freqvar = freqvar
        self.fullspec = fullspec
        self.forcefft = forcefft
        self.taper = taper
        self.taper_deltafreq = taper_deltafreq
        self.taper_beta = taper_beta
        self.resample = resample
        self.gapbreak_type = gapbreak_type
        self.gapbreak_value = gapbreak_value
        self.padmode = padmode
        self.padfrac = padfrac
        self.nowarn = nowarn
        self.save_fouriercoeffs = save_fouriercoeffs

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-fourierfilter", self.mode]
        # minfreq/maxfreq use the parser's keyword-prefixed form
        # (fix/var/expr/fixcolumn/list); _pval(..., keyword="fix")
        # prepends "fix" in front of bare numbers.
        if self.mode in ("highpass", "bandpass", "bandcut"):
            args += ["minfreq"] + _pval(self.minfreq, "fix")
        if self.mode in ("lowpass", "bandpass", "bandcut"):
            args += ["maxfreq"] + _pval(self.maxfreq, "fix")
        if self.filterexpr is not None:
            args += ["filterexpr", str(self.filterexpr)]
            if self.freqvar is not None:
                args += ["freqvar", str(self.freqvar)]
        if self.fullspec:
            args += ["fullspec"]
        if self.forcefft:
            args += ["forcefft"]
        if self.taper is not None:
            args += ["taper", self.taper, "deltafreq", str(self.taper_deltafreq)]
            if self.taper == "kaiser":
                args += ["beta", str(self.taper_beta)]
        if self.resample is not None:
            if isinstance(self.resample, str) and self.resample == "delmin":
                args += ["resample", "delmin"]
            elif isinstance(self.resample, (int, float)):
                args += ["resample", "fix", str(self.resample)]
            else:
                # arbitrary string → expression form
                args += ["resample", "expr", str(self.resample)]
            if self.gapbreak_type is not None:
                args += ["gapbreak", self.gapbreak_type, str(self.gapbreak_value)]
        if self.padmode is not None:
            args += ["padmode", self.padmode]
            if self.padfrac is not None:
                args += ["padfrac", str(self.padfrac)]
        # ofourier is a keyword-gated output: `ofourier <outdir>` (no 0/1
        # flag), so we emit it inline rather than using _outtoken.
        spec = _norm_save(self.save_fouriercoeffs)
        if _should_emit(spec):
            args += ["ofourier", spec.path if spec.path is not None else outdir]
        if self.nowarn:
            args += ["nowarn"]
        return args

    def _output_file_specs(self):
        return {"fouriercoeffs": (".fouriercoeffs", None)}


class expr(VartoolsCommand):
    """Evaluate an analytic expression to create or update a variable.

    Parameters
    ----------
    expression : str
        Expression of the form ``varname=expression``.
    vartype : str or None
        Variable type for the left-hand-side when creating a new variable:

        - ``None`` (default) — per-observation light-curve vector.
        - ``"listvar"`` — per-star variable (one value per light curve,
          persists across all LCs).  LC vectors on the RHS are evaluated
          at the first observation.
        - ``"scalar"`` — per-thread scalar.
        - ``"const"`` — global constant (single value, same for all LCs).

        If the variable already exists, its type is preserved regardless
        of this setting.
    outputcolumn : bool, optional
        If True, expose the LHS variable as an output column in the
        result table.  Only valid when ``vartype`` is one of
        ``"listvar"``, ``"scalar"``, or ``"const"``; raises ValueError
        otherwise (the value would otherwise be per-observation, not a
        single column).  The column name is ``Expr_<varname>_<idx>``
        (matching the vartools-C convention).

    Notes
    -----
    The expression engine supports aggregate functions that operate over
    all observations in the current light curve.  These are especially
    useful with ``vartype="listvar"`` to compute per-star summary
    statistics:

    - ``mean(x [, filter])``, ``median(x)``, ``stddev(x)``, ``MAD(x)``
    - ``vmin(x)``, ``vmax(x)``, ``sum(x)``
    - ``pct(x, pctval)``, ``wpct(x, w, pctval)``
    - ``weightedmean(x, w)``, ``wmedian(x, w)``
    - ``kurtosis(x)``, ``skewness(x)``, ``meddev(x)``, ``medmeddev(x)``

    All accept an optional filter argument, e.g. ``mean(mag, t>53730)``
    computes the mean of mag only for observations where t > 53730.

    Examples
    --------
    >>> cmd.expr("dmag=mag-10.0")                     # per-observation
    >>> cmd.expr("avg=mean(mag)", vartype="listvar")   # per-star mean
    >>> cmd.expr("pi=3.14159", vartype="const")        # global constant

    See Also
    --------
    CLI command: ``-expr``.
    """

    _vt_name = "expr"

    def __init__(
        self,
        expression: str,
        vartype: Optional[str] = None,
        outputcolumn: bool = False,
    ) -> None:
        if vartype is not None and vartype not in ("listvar", "scalar", "const"):
            raise ValueError(
                f"vartype must be None, 'listvar', 'scalar', or 'const', "
                f"got {vartype!r}"
            )
        if outputcolumn and vartype is None:
            raise ValueError(
                "expr: outputcolumn=True requires vartype to be one of "
                "'listvar', 'scalar', or 'const' (the default per-"
                "observation LC-vector type would produce one value per "
                "observation, not a single column)."
            )
        self.expression = expression
        self.vartype = vartype
        self.outputcolumn = outputcolumn

    def _to_cli_args(self) -> List[str]:
        args = ["-expr"]
        if self.vartype is not None:
            args.append(self.vartype)
        args.append(self.expression)
        if self.outputcolumn:
            args.append("outputcolumn")
        return args


class print_cols(VartoolsCommand):
    """Print selected variables to the output table as extra columns.

    Wraps the vartools ``-print`` command (the class is named
    ``print_cols`` to avoid shadowing the Python built-in; the
    pipeline-builder method is ``Pipeline.print``).

    Parameters
    ----------
    variables : str or list of str
        Variable name(s) to print.
    columnnames : str or list of str, optional
        Output column header names.
    fmt : str or list of str, optional
        printf-style format strings for each column.

    See Also
    --------
    CLI command: ``-print``.
    """

    _vt_name = "print"

    def __init__(
        self,
        variables: Union[str, List[str]],
        columnnames: Optional[Union[str, List[str]]] = None,
        fmt: Optional[Union[str, List[str]]] = None,
    ) -> None:
        self.variables = variables if isinstance(variables, str) else ",".join(variables)
        self.columnnames = columnnames
        self.fmt = fmt

    def _to_cli_args(self) -> List[str]:
        args = ["-print", self.variables]
        if self.columnnames is not None:
            cn = self.columnnames if isinstance(self.columnnames, str) else ",".join(self.columnnames)
            args += ["columnnames", cn]
        if self.fmt is not None:
            f = self.fmt if isinstance(self.fmt, str) else ",".join(self.fmt)
            args += ["format", f]
        return args


class FFT(VartoolsCommand):
    """Compute the Fast Fourier Transform of a pair of LC vectors.

    Transforms the (real, imaginary) input vector pair into a
    (real, imaginary) output vector pair.  For a purely real signal,
    pass a zero-filled vector as *input_imag*.

    Parameters
    ----------
    input_real, input_imag : str
        Input real and imaginary variable names.
    output_real, output_imag : str
        Output real and imaginary variable names (created if absent).

    See Also
    --------
    CLI command: ``-FFT``.  Inverse: :class:`IFFT`.
    """

    _vt_name = "FFT"

    def __init__(
        self,
        input_real: str,
        input_imag: str,
        output_real: str,
        output_imag: str,
    ) -> None:
        self.input_real = input_real
        self.input_imag = input_imag
        self.output_real = output_real
        self.output_imag = output_imag

    def _to_cli_args(self) -> List[str]:
        return ["-FFT", self.input_real, self.input_imag,
                self.output_real, self.output_imag]


class IFFT(VartoolsCommand):
    """Compute the inverse Fast Fourier Transform of a pair of LC vectors.

    Same parameter structure as :class:`FFT`: transforms the
    (real, imaginary) input vector pair into a (real, imaginary) output
    vector pair.

    Parameters
    ----------
    input_real, input_imag : str
        Input real and imaginary variable names.
    output_real, output_imag : str
        Output real and imaginary variable names (created if absent).

    See Also
    --------
    CLI command: ``-IFFT``.  Forward transform: :class:`FFT`.
    """

    _vt_name = "IFFT"

    def __init__(
        self,
        input_real: str,
        input_imag: str,
        output_real: str,
        output_imag: str,
    ) -> None:
        self.input_real = input_real
        self.input_imag = input_imag
        self.output_real = output_real
        self.output_imag = output_imag

    def _to_cli_args(self) -> List[str]:
        return ["-IFFT", self.input_real, self.input_imag,
                self.output_real, self.output_imag]


class resample(VartoolsCommand):
    """Resample the light curve onto a new time grid.

    Parameters
    ----------
    method : str
        Interpolation method: ``"nearest"``, ``"linear"``, ``"spline"``,
        ``"splinemonotonic"``, or ``"bspline"``.
    left : float, optional
        First-derivative boundary condition at the left end of the spline
        (only for ``method="spline"``).
    right : float, optional
        First-derivative boundary condition at the right end of the spline
        (only for ``method="spline"``).
    nbreaks : int, optional
        Number of break points for B-spline interpolation
        (only for ``method="bspline"``).
    order : int, optional
        Order of the B-spline (only for ``method="bspline"``).
    tstart, tstop : float or str, optional
        Start and stop of the new time grid.
    delt : float or str, optional
        Time step of the new grid.
    Npoints : int or str, optional
        Number of points in the new grid (alternative to delt).
    file_times : str, optional
        Source for the new time grid:

        * a path string → resample to the times in that file (CLI ``file fix <path>``);
        * the literal ``"list"`` → resample to the times in a per-LC file whose
          path is read from a column of the input list file (CLI ``file list``).
          Combine with ``list_column`` and ``t_column`` to control which
          columns are used.
    file_column : int, optional
        **Legacy.** Alias of ``t_column`` for the path-form (``file fix``)
        mode.  Prefer ``t_column``.
    list_column : int, optional
        Only meaningful with ``file_times="list"``.  Column number (1-based)
        in the *input list file* that holds the per-LC time-grid filename.
        Maps to the CLI ``listcolumn`` keyword.  When omitted, vartools
        consumes the next available list-file column.
    t_column : int, optional
        Column number (1-based) in the *time-grid file* that holds the time
        values.  Maps to the CLI ``column`` keyword (path mode) or
        ``tcolumn`` keyword (list mode).  Defaults to ``1``.
    gaps : str, optional
        Gap-handling spec passed through verbatim after the ``gaps``
        keyword (controls how interpolation behaves across large time
        gaps).

    See Also
    --------
    CLI command: ``-resample``.
    """

    _vt_name = "resample"

    def __init__(
        self,
        method: str = "linear",
        left: Optional[float] = None,
        right: Optional[float] = None,
        nbreaks: Optional[int] = None,
        order: Optional[int] = None,
        tstart=None,
        tstop=None,
        delt=None,
        Npoints=None,
        file_times: Optional[str] = None,
        file_column: Optional[int] = None,
        list_column: Optional[int] = None,
        t_column: Optional[int] = None,
        gaps: Optional[str] = None,
    ) -> None:
        self.method = method
        self.left = left
        self.right = right
        self.nbreaks = nbreaks
        self.order = order
        self.tstart = tstart
        self.tstop = tstop
        self.delt = delt
        self.Npoints = Npoints
        self.file_times = file_times
        self.file_column = file_column
        self.list_column = list_column
        self.t_column = t_column
        self.gaps = gaps

    def _to_cli_args(self) -> List[str]:
        args = ["-resample", self.method]
        if self.method == "spline":
            if self.left is not None:
                args += ["left", str(self.left)]
            if self.right is not None:
                args += ["right", str(self.right)]
        elif self.method == "bspline":
            if self.nbreaks is not None:
                args += ["nbreaks", str(self.nbreaks)]
            if self.order is not None:
                args += ["order", str(self.order)]
        # Resolve t_column from t_column or the legacy file_column alias.
        eff_t_column = self.t_column
        if eff_t_column is None and self.file_column is not None:
            eff_t_column = self.file_column
        if self.file_times is not None:
            toks = str(self.file_times).split()
            if toks[0] == "list":
                # New-grid times read from a per-LC file named in the input
                # list-file.  Accept either the kwarg form (preferred) or
                # extra tokens passed inside file_times itself (legacy).
                args += ["file", "list"]
                if len(toks) > 1:
                    # Legacy: user inlined "list listcolumn N tcolumn M" in
                    # the file_times string.  Pass the trailing tokens
                    # through verbatim and ignore the kwargs to avoid
                    # double-emission.
                    args += toks[1:]
                else:
                    if self.list_column is not None:
                        args += ["listcolumn", str(self.list_column)]
                    if eff_t_column is not None:
                        args += ["tcolumn", str(eff_t_column)]
            else:
                # Single file path → "fix" mode.
                args += ["file", "fix", str(self.file_times)]
                if eff_t_column is not None:
                    args += ["column", str(eff_t_column)]
        if self.gaps is not None:
            args += ["gaps"] + str(self.gaps).split()
        if self.tstart is not None:
            args += ["tstart"] + _pval(self.tstart, "fix")
        if self.tstop is not None:
            args += ["tstop"] + _pval(self.tstop, "fix")
        if self.delt is not None:
            args += ["delt"] + _pval(self.delt, "fix")
        elif self.Npoints is not None:
            args += ["Npoints"] + _pval(self.Npoints, "fix")
        return args


class decorr(VartoolsCommand):
    """Decorrelate the light curve against external trend vectors.

    Parameters
    ----------
    correct_lc : bool
        Subtract the decorrelation model from the light curve.
    zeropointterm : int
        Include (1) or exclude (0) a zero-point term per LC.
    subtractfirstterm : int
        Subtract (1) or keep (0) the first global term.
    global_files : list of (str, int)
        List of (filename, polynomial_order) for global trend vectors.
    lc_columns : list of (int, int)
        List of (column_number, polynomial_order) for LC-internal terms.
    save_model : bool, str, or Output, optional
        Directory for auxiliary file output.  A string is the
        **directory name** (not a filename); per-LC ``*.decorr.model``
        files are written inside it.  ``True`` captures the model files
        into ``result.files["decorr_model_N"]`` via a pyvartools-managed
        temp directory.  See the Auxiliary-output-files section of the
        pyvartools docs for the full ``Output`` semantics.  Default
        ``False`` (no model output).
    maskpoints : str, optional
        Name of a mask variable; points with ``maskvar ≤ 0`` are
        excluded from the fit.

    See Also
    --------
    CLI command: ``-decorr``.
    """

    _vt_name = "decorr"

    def __init__(
        self,
        correct_lc: bool = True,
        zeropointterm: int = 1,
        subtractfirstterm: int = 0,
        global_files: Optional[List] = None,
        lc_columns: Optional[List] = None,
        save_model: bool = False,
        maskpoints: Optional[str] = None,
    ) -> None:
        self.correct_lc = correct_lc
        self.zeropointterm = zeropointterm
        self.subtractfirstterm = subtractfirstterm
        self.global_files = global_files or []
        self.lc_columns = lc_columns or []
        self.save_model = save_model
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-decorr",
                "1" if self.correct_lc else "0",
                str(self.zeropointterm),
                str(self.subtractfirstterm),
                str(len(self.global_files))]
        for fname, order in self.global_files:
            args += [fname, str(order)]
        args += [str(len(self.lc_columns))]
        for col, order in self.lc_columns:
            args += [str(col), str(order)]
        # The -decorr CLI grammar treats `omodel` as a 0/1 flag in
        # this slot (parsecommandline.c:969 calls atoi() on it), not
        # as a literal keyword.  Help text writes "omodel" as the
        # placeholder name for the slot, which mis-led an earlier
        # version of this wrapper into emitting the literal word --
        # vartools then atoi()'d it to 0, skipped the path-reading
        # branch, and fell through with the path as a stray top-level
        # argument ("Invalid command or option <path>").
        m_spec = _norm_save(self.save_model)
        if _should_emit(m_spec):
            args += ["1", m_spec.path or outdir]
        else:
            args += ["0"]
        args += _flag("maskpoints", self.maskpoints)
        return args

    def _output_file_specs(self):
        return {"model": (".decorr.model", None)}


class Jstet(VartoolsCommand):
    """Stetson J / L variability statistics + kurtosis.

    Parameters
    ----------
    timescale : float
        Time gap (days) below which adjacent observations are treated as a
        single "close" pair (weight 1.0); otherwise as a singleton (weight 0.1).
    dates : str, optional
        Path to a file of JDs for *every possible* observation in the survey.
        Used to compute a survey-wide ``weight_max`` so the vartools J is
        ``J_stetson * (sum_w_actual / weight_max)`` — a multiplier that
        downweights LCs missing observations relative to the full schedule.
        Useful within a single survey; misleading when LCs come from surveys
        with different sampling. Mutually exclusive with ``skipnormalize``.
    skipnormalize : bool
        Skip the ``(sum_w / weight_max)`` rescaling and emit Stetson's
        original J and L. Use this when comparing across surveys / cadences,
        or when you want the textbook definition.
    maskpoints : str, optional
        Name of an LC vector; only points with ``maskvar > 0`` contribute.

    See Also
    --------
    CLI command: ``-Jstet``.
    Citation: Stetson 1996 (PASP 108, 851).
    """

    _vt_name = "Jstet"

    def __init__(
        self,
        timescale: float,
        dates: Optional[str] = None,
        skipnormalize: bool = False,
        maskpoints: Optional[str] = None,
    ) -> None:
        if (dates is None) == (not skipnormalize):
            raise ValueError(
                "Jstet: exactly one of 'dates' (path) or skipnormalize=True "
                "must be provided"
            )
        self.timescale = timescale
        self.dates = dates
        self.skipnormalize = skipnormalize
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        spec = "skipnormalize" if self.skipnormalize else self.dates
        return (["-Jstet", str(self.timescale), spec]
                + _flag("maskpoints", self.maskpoints))
