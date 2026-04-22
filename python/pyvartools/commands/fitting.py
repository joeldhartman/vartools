"""Fitting, modeling, and systematics command wrappers."""

from __future__ import annotations
from typing import List, Optional, Union

from pyvartools._command import VartoolsCommand
from ._helpers import _bool, _flag, _norm_save, _outtoken, _period_spec, _pval, _should_emit, _varexpr


class TFA(VartoolsCommand):
    """Trend Filtering Algorithm (TFA) detrending.

    Parameters
    ----------
    trendlist : str
        Path to the file listing trend light curves.
    dates_file : str
        Path to the dates file.
    pixelsep : float
        Pixel separation threshold for selecting trend stars.
    correct_lc : bool
        Subtract the TFA model from the light curve.
    save_coeffs : bool
        Write the TFA coefficients to a file.
    save_model : bool
        Write the TFA model to a file.
    xycol : tuple of int, optional
        ``(xcol, ycol)`` — column numbers for pixel coordinates.
    clip : float, optional
        Sigma-clipping threshold.
    usemedian : bool
        Use median instead of mean in clipping.
    useMAD : bool
        Use MAD in clipping.
    readformat : tuple, optional
        ``(Nskip, jdcol, magcol)`` — non-default light curve read format.
    trend_coeff_priors : str, optional
        Path to prior file for trend coefficients.
    weight_by_template_stddev : bool
        Weight by template standard deviation instead of LC errors.
    fitmask : str, optional
    outfitmask : str, optional
    """

    _vt_name = "TFA"

    def __init__(
        self,
        trendlist: str,
        dates_file: str,
        pixelsep: float,
        correct_lc: bool = True,
        save_coeffs: bool = False,
        save_model: bool = False,
        xycol: Optional[tuple] = None,
        clip: Optional[float] = None,
        usemedian: bool = False,
        useMAD: bool = False,
        readformat: Optional[tuple] = None,
        trend_coeff_priors: Optional[str] = None,
        weight_by_template_stddev: bool = False,
        fitmask: Optional[str] = None,
        outfitmask: Optional[str] = None,
    ) -> None:
        self.trendlist = trendlist
        self.dates_file = dates_file
        self.pixelsep = pixelsep
        self.correct_lc = correct_lc
        self.save_coeffs = save_coeffs
        self.save_model = save_model
        self.xycol = xycol
        self.clip = clip
        self.usemedian = usemedian
        self.useMAD = useMAD
        self.readformat = readformat
        self.trend_coeff_priors = trend_coeff_priors
        self.weight_by_template_stddev = weight_by_template_stddev
        self.fitmask = fitmask
        self.outfitmask = outfitmask

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-TFA", self.trendlist]
        if self.readformat is not None:
            args += ["readformat"] + [str(x) for x in self.readformat]
        if self.trend_coeff_priors is not None:
            args += ["trend_coeff_priors", self.trend_coeff_priors]
            if self.weight_by_template_stddev:
                args += ["weight_by_template_stddev"]
            else:
                args += ["use_lc_errors"]
        args += [self.dates_file, str(self.pixelsep)]
        if self.xycol is not None:
            args += ["xycol"] + [str(c) for c in self.xycol]
        args += ["1" if self.correct_lc else "0"]
        args += _outtoken(self.save_coeffs, outdir)
        args += _outtoken(self.save_model, outdir)
        if self.clip is not None:
            args += ["clip", str(self.clip)]
            args += _bool("usemedian", self.usemedian)
            args += _bool("useMAD", self.useMAD)
        args += _flag("fitmask", self.fitmask)
        args += _flag("outfitmask", self.outfitmask)
        return args

    def _output_file_specs(self):
        return {
            "coeffs": (".tfa.coeff", None),
            "model": (".tfa.model", None),
        }


class TFA_SR(VartoolsCommand):
    """TFA with signal reconstruction (TFA-SR).

    Parameters
    ----------
    trendlist : str
    dates_file : str
    pixelsep : float
    correct_lc : bool
    save_coeffs : bool
    save_model : bool
    dotfafirst : int
        Run TFA before SR (1) or SR before TFA (0).
    tfathresh : float
        Convergence threshold.
    maxiter : int
        Maximum number of TFA-SR iterations.
    signal_mode : str
        Signal model type: ``"bin"`` (nbins), ``"signal"`` (file),
        or ``"harm"`` (Nharm Nsubharm).
    signal_params : varies
        For ``"bin"``: nbins (int).
        For ``"signal"``: filename (str).
        For ``"harm"``: (Nharm, Nsubharm) tuple.
    clip, usemedian, useMAD, fitmask, outfitmask : see TFA.
    """

    _vt_name = "TFA_SR"

    def __init__(
        self,
        trendlist: str,
        dates_file: str,
        pixelsep: float,
        dotfafirst: int = 1,
        tfathresh: float = 0.001,
        maxiter: int = 10,
        signal_mode: str = "bin",
        signal_params=None,
        correct_lc: bool = True,
        save_coeffs: bool = False,
        save_model: bool = False,
        xycol: Optional[tuple] = None,
        clip: Optional[float] = None,
        usemedian: bool = False,
        useMAD: bool = False,
        readformat: Optional[tuple] = None,
        fitmask: Optional[str] = None,
        outfitmask: Optional[str] = None,
        decorr_params: Optional[str] = None,
        signal_period: Optional[Union[float, str]] = None,
    ) -> None:
        self.trendlist = trendlist
        self.dates_file = dates_file
        self.pixelsep = pixelsep
        self.dotfafirst = dotfafirst
        self.tfathresh = tfathresh
        self.maxiter = maxiter
        self.signal_mode = signal_mode
        self.signal_params = signal_params
        self.correct_lc = correct_lc
        self.save_coeffs = save_coeffs
        self.save_model = save_model
        self.xycol = xycol
        self.clip = clip
        self.usemedian = usemedian
        self.useMAD = useMAD
        self.readformat = readformat
        self.fitmask = fitmask
        self.outfitmask = outfitmask
        self.decorr_params = decorr_params
        self.signal_period = signal_period

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-TFA_SR", self.trendlist]
        if self.readformat is not None:
            args += ["readformat"] + [str(x) for x in self.readformat]
        args += [self.dates_file]
        if self.xycol is not None:
            args += ["xycol"] + [str(c) for c in self.xycol]
        args += [str(self.pixelsep)]
        args += ["1" if self.correct_lc else "0"]
        args += _outtoken(self.save_coeffs, outdir)
        args += _outtoken(self.save_model, outdir)
        args += [str(self.dotfafirst), str(self.tfathresh), str(self.maxiter)]
        # Signal model
        if self.signal_mode == "bin":
            nbins = self.signal_params if self.signal_params is not None else 100
            args += ["bin", str(nbins)]
        elif self.signal_mode == "signal":
            args += ["signal", str(self.signal_params)]
        elif self.signal_mode == "harm":
            nharm, nsubharm = self.signal_params or (3, 0)
            args += ["harm", str(nharm), str(nsubharm)]
        if self.signal_period is not None:
            args += ["period"] + _pval(self.signal_period, "fix")
        if self.decorr_params is not None:
            args += ["decorr"] + str(self.decorr_params).split()
        if self.clip is not None:
            args += ["clip", str(self.clip)]
            args += _bool("usemedian", self.usemedian)
            args += _bool("useMAD", self.useMAD)
        args += _flag("fitmask", self.fitmask)
        args += _flag("outfitmask", self.outfitmask)
        return args

    def _resolve_back_references(self, prev) -> None:
        from ._helpers import _resolve_period_backref
        self.signal_period = _resolve_period_backref(prev, self.signal_period)

    def _output_file_specs(self):
        return {
            "coeffs": (".tfa.coeff", None),
            "model": (".tfa.model", None),
        }


class SYSREM(VartoolsCommand):
    """SYSREM systematic noise removal.

    Parameters
    ----------
    ninput_color : int
        Number of colour terms.
    col : int, optional
        Column number for the colour.
    ninput_airmass : int
        Number of airmass terms.
    initial_airmass_file : str
        File with initial airmass values.
    sigma_clip1, sigma_clip2 : float
        Sigma-clipping thresholds.
    saturation : float
        Saturation level.
    correct_lc : bool
        Subtract the SYSREM model from the light curve.
    save_model : bool
        Write the SYSREM model to a file.
    save_trends : bool
        Write the trend vectors.
    useweights : int
        Use measurement weights (0 or 1).
    """

    _vt_name = "SYSREM"

    def __init__(
        self,
        ninput_color: int,
        ninput_airmass: int,
        initial_airmass_file: str,
        sigma_clip1: float = 5.0,
        sigma_clip2: float = 5.0,
        saturation: float = 1e9,
        correct_lc: bool = True,
        save_model: bool = False,
        save_trends: bool = False,
        useweights: int = 1,
        col: Optional[int] = None,
    ) -> None:
        self.ninput_color = ninput_color
        self.col = col
        self.ninput_airmass = ninput_airmass
        self.initial_airmass_file = initial_airmass_file
        self.sigma_clip1 = sigma_clip1
        self.sigma_clip2 = sigma_clip2
        self.saturation = saturation
        self.correct_lc = correct_lc
        self.save_model = save_model
        self.save_trends = save_trends
        self.useweights = useweights

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-SYSREM", str(self.ninput_color)]
        args += _flag("column", self.col)
        args += [str(self.ninput_airmass), self.initial_airmass_file,
                 str(self.sigma_clip1), str(self.sigma_clip2),
                 str(self.saturation),
                 "1" if self.correct_lc else "0"]
        args += _outtoken(self.save_model, outdir)
        tr_spec = _norm_save(self.save_trends)
        if _should_emit(tr_spec):
            args += ["otrends", tr_spec.path or outdir]
        else:
            args += ["0"]
        args += [str(self.useweights)]
        return args

    def _output_file_specs(self):
        return {
            "model": (".sysrem.model", None),
            "trends": (".sysrem.trends", None),
        }


class MandelAgolTransit(VartoolsCommand):
    """Mandel-Agol transit model fitting.

    Parameters
    ----------
    P0, T00, r0, a0 : float
        Initial period, epoch, radius ratio, semi-major axis (in stellar radii).
    inclination : float
        Initial orbital inclination (degrees).  Used when ``bimpact`` is None.
    bimpact : float, optional
        Initial impact parameter.  When set, replaces ``inclination`` and emits
        ``"b" bimpact`` instead of ``"i" inclination``.
    e0, omega0 : float
        Initial eccentricity and argument of periastron.
    mconst0 : float
        Out-of-transit magnitude constant.
    ld_type : str
        ``"quad"`` or ``"nonlin"`` limb darkening.
    ld_coeffs : list of float
        Initial limb-darkening coefficients.
    fitephem : int
        Fit the transit epoch (and period) (0 = fixed, 1 = free).
    fitr, fita, fitinclterm : int
        Fit radius ratio, semi-major axis, inclination.
    fite, fitomega, fitmconst : int
        Fit eccentricity, omega, magnitude constant.
    fitldcoeffs : list of int
        Fit each LD coefficient (0 or 1 per coefficient).
    rv_file : str, optional
        Path to RV input file (first col JD, second RV, third RV error).
        When set, fitRV=1 is emitted.
    rv_model_file : str, optional
        Path for RV model output.
    K0 : float, optional
        Initial RV semi-amplitude.
    gamma0 : float, optional
        Initial RV systemic velocity.
    fitK, fitgamma : int
        Whether to fit K and gamma (0 or 1).
    correct_lc : bool
    save_model : bool
    modelvar : str, optional
        Variable name to store the best-fit model.
    save_phcurve : bool
    ophcurve_phmin, ophcurve_phmax, ophcurve_phstep : float
        Phase range and step for the phcurve output (only used when
        ``save_phcurve`` is truthy).  Defaults are ``0.0``, ``1.0``, ``0.005``.
    save_jdcurve : bool
    ojdcurve_jdstep : float
        Time step for the jdcurve output (only used when ``save_jdcurve`` is
        truthy).  Default is ``0.02``.
    """

    _vt_name = "MandelAgolTransit"

    def __init__(
        self,
        P0: float,
        T00: float,
        r0: float = 0.1,
        a0: float = 10.0,
        inclination: float = 90.0,
        bimpact: Optional[float] = None,
        e0: float = 0.0,
        omega0: float = 0.0,
        mconst0: float = 0.0,
        ld_type: str = "quad",
        ld_coeffs: Optional[List[float]] = None,
        fitephem: int = 1,
        fitr: int = 1,
        fita: int = 1,
        fitinclterm: int = 1,
        fite: int = 0,
        fitomega: int = 0,
        fitmconst: int = 1,
        fitldcoeffs: Optional[List[int]] = None,
        rv_file: Optional[str] = None,
        rv_model_file: Optional[str] = None,
        K0: Optional[float] = None,
        gamma0: Optional[float] = None,
        fitK: int = 0,
        fitgamma: int = 0,
        correct_lc: bool = False,
        save_model: bool = False,
        modelvar: Optional[str] = None,
        save_phcurve: bool = False,
        ophcurve_phmin: float = 0.0,
        ophcurve_phmax: float = 1.0,
        ophcurve_phstep: float = 0.005,
        save_jdcurve: bool = False,
        ojdcurve_jdstep: float = 0.02,
    ) -> None:
        self.P0 = P0
        self.T00 = T00
        self.r0 = r0
        self.a0 = a0
        self.inclination = inclination
        self.bimpact = bimpact
        self.e0 = e0
        self.omega0 = omega0
        self.mconst0 = mconst0
        self.ld_type = ld_type
        self.ld_coeffs = ld_coeffs or [0.3, 0.3]
        self.fitephem = fitephem
        self.fitr = fitr
        self.fita = fita
        self.fitinclterm = fitinclterm
        self.fite = fite
        self.fitomega = fitomega
        self.fitmconst = fitmconst
        self.fitldcoeffs = fitldcoeffs or [0] * len(self.ld_coeffs)
        self.rv_file = rv_file
        self.rv_model_file = rv_model_file
        self.K0 = K0
        self.gamma0 = gamma0
        self.fitK = fitK
        self.fitgamma = fitgamma
        self.correct_lc = correct_lc
        self.save_model = save_model
        self.modelvar = modelvar
        self.save_phcurve = save_phcurve
        self.ophcurve_phmin = ophcurve_phmin
        self.ophcurve_phmax = ophcurve_phmax
        self.ophcurve_phstep = ophcurve_phstep
        self.save_jdcurve = save_jdcurve
        self.ojdcurve_jdstep = ojdcurve_jdstep

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-MandelAgolTransit"]
        # "bls"/"blsfixper" are single-token initial-param specs that pull
        # P0/T00/r0/a0/e0/omega0 from the prior BLS command.
        if isinstance(self.P0, str) and self.P0.strip() in ("bls", "blsfixper"):
            args += [self.P0.strip()]
        else:
            args += _varexpr(self.P0) + _varexpr(self.T00) + _varexpr(self.r0) + _varexpr(self.a0)
            if self.bimpact is not None:
                args += ["b"] + _varexpr(self.bimpact)
            else:
                args += ["i"] + _varexpr(self.inclination)
            args += _varexpr(self.e0) + _varexpr(self.omega0) + _varexpr(self.mconst0)
        args += [self.ld_type]
        for c in self.ld_coeffs:
            args += _varexpr(c)
        args += [str(self.fitephem), str(self.fitr), str(self.fita),
                 str(self.fitinclterm), str(self.fite), str(self.fitomega),
                 str(self.fitmconst)] + [str(f) for f in self.fitldcoeffs]
        if self.rv_file is not None:
            args += ["1", self.rv_file,
                     self.rv_model_file or "0"]
            args += _varexpr(self.K0 or 0.0) + _varexpr(self.gamma0 or 0.0)
            args += [str(self.fitK), str(self.fitgamma)]
        else:
            args += ["0"]  # fitRV=0
        args += ["1" if self.correct_lc else "0"]
        args += _outtoken(self.save_model, outdir)
        m_spec = _norm_save(self.save_model)
        if _should_emit(m_spec) and self.modelvar is not None:
            args += ["modelvar", self.modelvar]
        ph_spec = _norm_save(self.save_phcurve)
        if _should_emit(ph_spec):
            args += ["ophcurve", ph_spec.path or outdir,
                     str(self.ophcurve_phmin),
                     str(self.ophcurve_phmax),
                     str(self.ophcurve_phstep)]
        jd_spec = _norm_save(self.save_jdcurve)
        if _should_emit(jd_spec):
            args += ["ojdcurve", jd_spec.path or outdir,
                     str(self.ojdcurve_jdstep)]
        return args

    def _resolve_back_references(self, prev) -> None:
        # When P0 is "bls" / "blsfixper", pull P, Tc, depth, qtran from the
        # prior -BLS (or -BLSFixPer) and convert to (P0, T00, r0, a0).  The
        # C-side formulas are:  r0 = sqrt(depth), a0 = 1/(qtran*pi).  Other
        # initial params are set to the vartools defaults used on the C side.
        import math
        from ._helpers import _resolve_bls_transit_backref
        from pyvartools.perlc import PerLC

        if not (isinstance(self.P0, str)
                and self.P0.strip() in ("bls", "blsfixper")):
            return
        d = _resolve_bls_transit_backref(prev, self.P0.strip())

        def _sqrt(v):
            if isinstance(v, PerLC):
                return PerLC([math.sqrt(float(x)) for x in v])
            return math.sqrt(float(v))

        def _inv_pi(v):
            if isinstance(v, PerLC):
                return PerLC([1.0 / (float(x) * math.pi) for x in v])
            return 1.0 / (float(v) * math.pi)

        self.P0 = d["period"]
        self.T00 = d["T0"]
        self.r0 = _sqrt(d["depth"])
        self.a0 = _inv_pi(d["qtran"])
        self.e0 = 0.0
        self.omega0 = 0.0
        # Match the C defaults (inclination=90 unless bimpact is supplied).
        if self.bimpact is None:
            self.inclination = 90.0
        self.mconst0 = 0.0

    def _output_file_specs(self):
        return {
            "model": (".mandelagoltransit.model", None),
            "phcurve": (".mandelagoltransit.phcurve", None),
            "jdcurve": (".mandelagoltransit.jdcurve", None),
        }


class SoftenedTransit(VartoolsCommand):
    """Softened (trapezoidal) transit model fitting.

    Parameters
    ----------
    init_params : str or tuple
        Either ``"bls"``, ``"blsfixper"``, or a tuple
        ``(P0, T00, eta0, delta0, mconst0, cval0)``.
    fitephem, fiteta, fitcval, fitdelta, fitmconst : int
        Which parameters to fit (0 or 1).
    correct_lc : bool
    save_model : bool
    fit_harm : int
        Harmonic fit flag (0 = no harmonic fit, default).  When > 0, the
        following three parameters describe the harmonic fitting.
    fit_harm_method : str, optional
        Method for harmonic fitting, e.g. ``"aov"``.
    fit_harm_nharm : int, optional
        Number of harmonics for harmonic fitting.
    fit_harm_nsubharm : int, optional
        Number of sub-harmonics for harmonic fitting.
    """

    _vt_name = "SoftenedTransit"

    def __init__(
        self,
        init_params="bls",
        fitephem: int = 1,
        fiteta: int = 1,
        fitcval: int = 1,
        fitdelta: int = 1,
        fitmconst: int = 1,
        correct_lc: bool = False,
        save_model: bool = False,
        fit_harm: int = 0,
        fit_harm_method: Optional[str] = None,
        fit_harm_nharm: Optional[int] = None,
        fit_harm_nsubharm: Optional[int] = None,
    ) -> None:
        self.init_params = init_params
        self.fitephem = fitephem
        self.fiteta = fiteta
        self.fitcval = fitcval
        self.fitdelta = fitdelta
        self.fitmconst = fitmconst
        self.correct_lc = correct_lc
        self.save_model = save_model
        self.fit_harm = fit_harm
        self.fit_harm_method = fit_harm_method
        self.fit_harm_nharm = fit_harm_nharm
        self.fit_harm_nsubharm = fit_harm_nsubharm

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-SoftenedTransit"]
        if isinstance(self.init_params, str):
            args.append(self.init_params)
        else:
            args += [str(x) for x in self.init_params]
        args += [str(self.fitephem), str(self.fiteta), str(self.fitcval),
                 str(self.fitdelta), str(self.fitmconst),
                 "1" if self.correct_lc else "0"]
        args += _outtoken(self.save_model, outdir)
        if self.fit_harm > 0:
            args += [str(self.fit_harm)]
            if self.fit_harm_method is not None:
                # A resolved back-ref may be a numeric value — emit with the
                # "fix" keyword so the CLI parses it correctly.
                if isinstance(self.fit_harm_method, (int, float)):
                    args += ["fix", str(self.fit_harm_method)]
                else:
                    args += [str(self.fit_harm_method)]
            if self.fit_harm_nharm is not None:
                args += [str(self.fit_harm_nharm)]
            if self.fit_harm_nsubharm is not None:
                args += [str(self.fit_harm_nsubharm)]
        else:
            args += ["0"]
        return args

    def _resolve_back_references(self, prev) -> None:
        from ._helpers import (_resolve_bls_transit_backref,
                                _resolve_period_backref)
        from pyvartools.perlc import PerLC

        # init_params: "bls" / "blsfixper" → 6-tuple (P, T0, eta, delta,
        # mconst, cval) filled from the prior -BLS / -BLSFixPer.  The
        # vartools C defaults are mconst0=-1 (auto-estimate) and cval0=0.2.
        if isinstance(self.init_params, str) and self.init_params.strip() in (
                "bls", "blsfixper"):
            d = _resolve_bls_transit_backref(prev, self.init_params.strip())
            P, T0, eta, delta = d["period"], d["T0"], d["qtran"], d["depth"]
            if any(isinstance(v, PerLC) for v in (P, T0, eta, delta)):
                raise NotImplementedError(
                    "SoftenedTransit init_params='bls' across a batch chain "
                    "boundary is not supported (would need 4 per-LC columns "
                    "injected).  Use single-LC chaining or a Pipeline."
                )
            self.init_params = (float(P), float(T0), float(eta),
                                float(delta), -1.0, 0.2)

        # fit_harm_method: supports ls / aov / bls back-refs when fit_harm>0.
        if self.fit_harm and self.fit_harm_method:
            self.fit_harm_method = _resolve_period_backref(
                prev, self.fit_harm_method)

    def _output_file_specs(self):
        return {"model": (".softenedtransit.model", None)}


class Starspot(VartoolsCommand):
    """Starspot model fitting.

    Parameters
    ----------
    period : float or str
        Period for the starspot model.
    a0 : float
        Initial spot fractional radius.
    b0 : float
        Initial spot latitude.
    alpha0 : float
        Initial spot longitude.
    i0 : float
        Initial stellar inclination (degrees).
    chi0 : float
        Initial spot contrast.
    psi00 : float
        Initial spot phase offset.
    mconst0 : float
        Initial magnitude offset.
    fit_period : int
        Fit the period (0=fixed, 1=free).
    fit_a, fit_b, fit_alpha, fit_i, fit_chi, fit_psi, fit_mconst : int
        Fit each parameter (0=fixed, 1=free).
    correct_lc : bool
    save_model : bool
    """

    _vt_name = "Starspot"

    def __init__(
        self,
        period="ls",
        a0: float = 0.1,
        b0: float = 0.5,
        alpha0: float = 20.0,
        i0: float = 85.0,
        chi0: float = 30.0,
        psi00: float = 0.0,
        mconst0: float = 0.0,
        fit_period: int = 1,
        fit_a: int = 1,
        fit_b: int = 1,
        fit_alpha: int = 1,
        fit_i: int = 1,
        fit_chi: int = 1,
        fit_psi: int = 1,
        fit_mconst: int = 1,
        correct_lc: bool = False,
        save_model: bool = False,
    ) -> None:
        self.period = period
        self.a0 = a0
        self.b0 = b0
        self.alpha0 = alpha0
        self.i0 = i0
        self.chi0 = chi0
        self.psi00 = psi00
        self.mconst0 = mconst0
        self.fit_period = fit_period
        self.fit_a = fit_a
        self.fit_b = fit_b
        self.fit_alpha = fit_alpha
        self.fit_i = fit_i
        self.fit_chi = fit_chi
        self.fit_psi = fit_psi
        self.fit_mconst = fit_mconst
        self.correct_lc = correct_lc
        self.save_model = save_model

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-Starspot"] + _period_spec(self.period)
        # Initial parameter values: CLI order is a b alpha i chi psi mconst
        args += (_varexpr(self.a0) + _varexpr(self.b0) + _varexpr(self.alpha0)
                 + _varexpr(self.i0) + _varexpr(self.chi0) + _varexpr(self.psi00)
                 + _varexpr(self.mconst0))
        # Fit flags: CLI order is fitP fita fitb fitalpha fiti fitchi fitpsi fitmconst
        args += [str(self.fit_period), str(self.fit_a), str(self.fit_b),
                 str(self.fit_alpha), str(self.fit_i), str(self.fit_chi),
                 str(self.fit_psi), str(self.fit_mconst)]
        args += ["1" if self.correct_lc else "0"]
        args += _outtoken(self.save_model, outdir)
        return args

    def _resolve_back_references(self, prev) -> None:
        from ._helpers import _resolve_period_backref
        self.period = _resolve_period_backref(prev, self.period)

    def _output_file_specs(self):
        return {"model": (".starspot.model", None)}


class microlens(VartoolsCommand):
    """Microlensing event model fitting.

    Parameters
    ----------
    f0, f1, u0, t0, tmax : float or str or None
        Initial values for microlensing parameters.  Each can be a number
        (auto-fit), ``"auto"`` (vartools auto-estimates), or None (omit).
    correct_lc : bool
    save_model : bool
    """

    _vt_name = "microlens"

    def __init__(
        self,
        f0=None,
        f1=None,
        u0=None,
        t0=None,
        tmax=None,
        correct_lc: bool = False,
        save_model: bool = False,
        f0_step: Optional[float] = None,
        f0_novary: bool = False,
        f1_step: Optional[float] = None,
        f1_novary: bool = False,
        u0_step: Optional[float] = None,
        u0_novary: bool = False,
        t0_step: Optional[float] = None,
        t0_novary: bool = False,
        tmax_step: Optional[float] = None,
        tmax_novary: bool = False,
    ) -> None:
        self.f0 = f0
        self.f1 = f1
        self.u0 = u0
        self.t0 = t0
        self.tmax = tmax
        self.correct_lc = correct_lc
        self.save_model = save_model
        self.f0_step = f0_step
        self.f0_novary = f0_novary
        self.f1_step = f1_step
        self.f1_novary = f1_novary
        self.u0_step = u0_step
        self.u0_novary = u0_novary
        self.t0_step = t0_step
        self.t0_novary = t0_novary
        self.tmax_step = tmax_step
        self.tmax_novary = tmax_novary

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-microlens"]
        param_specs = [
            ("f0",   self.f0,   self.f0_step,   self.f0_novary),
            ("f1",   self.f1,   self.f1_step,   self.f1_novary),
            ("u0",   self.u0,   self.u0_step,   self.u0_novary),
            ("t0",   self.t0,   self.t0_step,   self.t0_novary),
            ("tmax", self.tmax, self.tmax_step, self.tmax_novary),
        ]
        for pname, val, step, novary in param_specs:
            if val is not None:
                args += [pname]
                if val == "auto":
                    args.append("auto")
                elif isinstance(val, (int, float)):
                    args += ["fix", str(val)]
                else:
                    args += str(val).split()
                if step is not None:
                    args += ["step", str(step)]
                if novary:
                    args += ["novary"]
        args += _bool("correctlc", self.correct_lc)
        m_spec = _norm_save(self.save_model)
        if _should_emit(m_spec):
            args += ["omodel", m_spec.path or outdir]
        return args

    def _resolve_back_references(self, prev) -> None:
        # Each of f0/f1/u0/t0/tmax may be given as "fixcolumn NAME" — resolve
        # each one independently against the prior step's output.
        from ._helpers import _resolve_period_backref
        self.f0 = _resolve_period_backref(prev, self.f0)
        self.f1 = _resolve_period_backref(prev, self.f1)
        self.u0 = _resolve_period_backref(prev, self.u0)
        self.t0 = _resolve_period_backref(prev, self.t0)
        self.tmax = _resolve_period_backref(prev, self.tmax)

    def _output_file_specs(self):
        return {"model": (".microlens", None)}


class nonlinfit(VartoolsCommand):
    """Non-linear least-squares fitting of an analytic function.

    Parameters
    ----------
    function : str
        Analytic function string.
    paramlist : str
        Parameter list string (name:initial[:step[:min:max]], ...).
    optimizer : str
        ``"amoeba"`` (default) or ``"mcmc"``.
    linfit_params : str, optional
        Linear parameters to solve analytically (space-separated names).
    errors : str, optional
        Expression for per-point errors (variable name or expression).
    covariance : str, optional
        Covariance model tokens, e.g. ``"squareexp amp_v rho_v"``
        (passed verbatim after the ``covariance`` keyword).
    priors : str, optional
        Prior expression (passed verbatim after the ``priors`` keyword).
    constraints : str, optional
        Constraint expression (passed verbatim after ``constraints``).
    amoeba_tolerance : float, optional
        Convergence tolerance for the amoeba optimizer.
    amoeba_maxsteps : int, optional
        Maximum number of steps for the amoeba optimizer.
    mcmc_naccept : int, optional
        Number of accepted links for MCMC.
    mcmc_nlinkstotal : int, optional
        Total number of MCMC links.
    mcmc_fracburnin : float, optional
        Fraction of links to discard as burn-in.
    mcmc_eps : float, optional
        Initial step size for MCMC.
    mcmc_skipamoeba : bool
        Skip the initial amoeba optimization in MCMC mode.
    mcmc_maxmemstore : int, optional
        Maximum number of links to store in memory.
    mcmc_outchains : bool or str or Output
        Write MCMC chains to a file.
    mcmc_chains_format : str, optional
        Naming format for chain output files.
    mcmc_chains_printevery : int, optional
        Print every Nth link to the chain file.
    correct_lc : bool
    save_model : bool or str or Output
    model_nameformat : str, optional
        Naming format for model output files.
    modelvar : str, optional
        Variable name to store the best-fit model.
    fitmask : str, optional
        Name of a mask variable; non-zero points are excluded from the fit.
    """

    _vt_name = "nonlinfit"

    def __init__(
        self,
        function: str,
        paramlist: str,
        optimizer: str = "amoeba",
        linfit_params: Optional[str] = None,
        errors: Optional[str] = None,
        covariance: Optional[str] = None,
        priors: Optional[str] = None,
        constraints: Optional[str] = None,
        amoeba_tolerance: Optional[float] = None,
        amoeba_maxsteps: Optional[int] = None,
        mcmc_naccept: Optional[int] = None,
        mcmc_nlinkstotal: Optional[int] = None,
        mcmc_fracburnin: Optional[float] = None,
        mcmc_eps: Optional[float] = None,
        mcmc_skipamoeba: bool = False,
        mcmc_maxmemstore: Optional[int] = None,
        mcmc_outchains=False,
        mcmc_chains_format: Optional[str] = None,
        mcmc_chains_printevery: Optional[int] = None,
        correct_lc: bool = False,
        save_model=False,
        model_nameformat: Optional[str] = None,
        modelvar: Optional[str] = None,
        fitmask: Optional[str] = None,
    ) -> None:
        self.function = function
        self.paramlist = paramlist
        self.optimizer = optimizer
        self.linfit_params = linfit_params
        self.errors = errors
        self.covariance = covariance
        self.priors = priors
        self.constraints = constraints
        self.amoeba_tolerance = amoeba_tolerance
        self.amoeba_maxsteps = amoeba_maxsteps
        self.mcmc_naccept = mcmc_naccept
        self.mcmc_nlinkstotal = mcmc_nlinkstotal
        self.mcmc_fracburnin = mcmc_fracburnin
        self.mcmc_eps = mcmc_eps
        self.mcmc_skipamoeba = mcmc_skipamoeba
        self.mcmc_maxmemstore = mcmc_maxmemstore
        self.mcmc_outchains = mcmc_outchains
        self.mcmc_chains_format = mcmc_chains_format
        self.mcmc_chains_printevery = mcmc_chains_printevery
        self.correct_lc = correct_lc
        self.save_model = save_model
        self.model_nameformat = model_nameformat
        self.modelvar = modelvar
        self.fitmask = fitmask

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-nonlinfit", self.function, self.paramlist]
        # Optional pre-optimizer args (CLI order: linfit errors covariance priors constraints)
        if self.linfit_params is not None:
            args += ["linfit", self.linfit_params]
        if self.errors is not None:
            args += ["errors", self.errors]
        if self.covariance is not None:
            args += ["covariance"] + str(self.covariance).split()
        if self.priors is not None:
            args += ["priors"] + str(self.priors).split()
        if self.constraints is not None:
            args += ["constraints"] + str(self.constraints).split()
        # Optimizer and optimizer-specific options
        args += [self.optimizer]
        if self.optimizer == "amoeba":
            if self.amoeba_tolerance is not None:
                args += ["tolerance"] + _varexpr(self.amoeba_tolerance)
            if self.amoeba_maxsteps is not None:
                args += ["maxsteps"] + _varexpr(self.amoeba_maxsteps)
        elif self.optimizer == "mcmc":
            if self.mcmc_naccept is not None:
                args += ["Naccept"] + _varexpr(self.mcmc_naccept)
            if self.mcmc_nlinkstotal is not None:
                args += ["Nlinkstotal"] + _varexpr(self.mcmc_nlinkstotal)
            if self.mcmc_fracburnin is not None:
                args += ["fracburnin"] + _varexpr(self.mcmc_fracburnin)
            if self.mcmc_eps is not None:
                args += ["eps"] + _varexpr(self.mcmc_eps)
            if self.mcmc_skipamoeba:
                args += ["skipamoeba"]
            if self.mcmc_maxmemstore is not None:
                args += ["maxmemstore", str(self.mcmc_maxmemstore)]
            chains_spec = _norm_save(self.mcmc_outchains)
            if _should_emit(chains_spec):
                args += ["outchains", chains_spec.path or outdir]
                if self.mcmc_chains_format is not None:
                    args += ["format", self.mcmc_chains_format]
                if self.mcmc_chains_printevery is not None:
                    args += ["printevery", str(self.mcmc_chains_printevery)]
        args += _bool("correctlc", self.correct_lc)
        m_spec = _norm_save(self.save_model)
        if _should_emit(m_spec):
            args += ["omodel", m_spec.path or outdir]
            if self.model_nameformat is not None:
                args += ["format", self.model_nameformat]
            if self.modelvar is not None:
                args += ["modelvar", self.modelvar]
        args += _flag("fitmask", self.fitmask)
        return args

    def _output_file_specs(self):
        return {
            "model": (".nonlinfit.model", None),
            "chains": (".nonlinfit.chains", None),
        }


class addnoise(VartoolsCommand):
    """Add synthetic noise to the light curve.

    Parameters
    ----------
    noise_type : str
        ``"white"``, ``"squareexp"``, ``"exp"``, ``"matern"``, or ``"wavelet"``.
    sig_white : float or str
        White noise amplitude.
    rho : float or str, optional
        Correlation length (for ``"squareexp"``, ``"exp"``, ``"matern"``).
    sig_red : float or str, optional
        Red noise amplitude (for correlated noise models).
    nu : float or str, optional
        Smoothness parameter for the Matern covariance (``"matern"`` only).
    gamma : float or str, optional
        Wavelet decay parameter (``"wavelet"`` only).
    bintime : float or str, optional
        Bin time for integrated covariance (``"squareexp"``, ``"exp"``).
    """

    _vt_name = "addnoise"

    def __init__(
        self,
        noise_type: str = "white",
        sig_white=0.001,
        rho=None,
        sig_red=None,
        nu=None,
        gamma=None,
        bintime=None,
    ) -> None:
        self.noise_type = noise_type
        self.sig_white = sig_white
        self.rho = rho
        self.sig_red = sig_red
        self.nu = nu
        self.gamma = gamma
        self.bintime = bintime

    def _to_cli_args(self) -> List[str]:
        args = ["-addnoise", self.noise_type]
        if self.noise_type == "white":
            args += ["sig_white"] + _pval(self.sig_white, "fix")
        elif self.noise_type == "matern":
            if self.nu is not None:
                args += ["nu"] + _pval(self.nu, "fix")
            if self.rho is not None:
                args += ["rho"] + _pval(self.rho, "fix")
            if self.sig_red is not None:
                args += ["sig_red"] + _pval(self.sig_red, "fix")
            args += ["sig_white"] + _pval(self.sig_white, "fix")
        elif self.noise_type == "wavelet":
            if self.gamma is not None:
                args += ["gamma"] + _pval(self.gamma, "fix")
            if self.sig_red is not None:
                args += ["sig_red"] + _pval(self.sig_red, "fix")
            args += ["sig_white"] + _pval(self.sig_white, "fix")
        else:
            # squareexp, exp
            if self.rho is not None:
                args += ["rho"] + _pval(self.rho, "fix")
            if self.sig_red is not None:
                args += ["sig_red"] + _pval(self.sig_red, "fix")
            if self.bintime is not None:
                args += ["bintime"] + _pval(self.bintime, "fix")
            args += ["sig_white"] + _pval(self.sig_white, "fix")
        return args


class findblends(VartoolsCommand):
    """Search for blended transit signals from nearby stars.

    Parameters
    ----------
    matchrad : float
        Search radius (arcsec).
    period : float or str
        Transit period source.  Valid values: ``"list"``,
        ``"fix <period>"``, ``"fixcolumn <colname|colnum>"``.
    radec : bool
        Use RA/Dec coordinates instead of pixel coordinates.
    nharm : int
        Number of harmonics.
    xycol : tuple of str, optional
        ``(xcol, ycol)`` column names/numbers for pixel coordinates.
    starlist : str, optional
        Path to a star list file.
    zeromag : float, optional
        Zero-point magnitude.
    nofluxconvert : bool
        Do not convert fluxes.
    save_matches : bool
        Write matched star list to a file.
    """

    _vt_name = "findblends"

    def __init__(
        self,
        matchrad: float,
        period="list",
        radec: bool = False,
        nharm: int = 1,
        xycol: Optional[tuple] = None,
        starlist: Optional[str] = None,
        zeromag: Optional[float] = None,
        nofluxconvert: bool = False,
        save_matches: bool = False,
    ) -> None:
        self.matchrad = matchrad
        self.period = period
        self.radec = radec
        self.nharm = nharm
        self.xycol = xycol
        self.starlist = starlist
        self.zeromag = zeromag
        self.nofluxconvert = nofluxconvert
        self.save_matches = save_matches

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-findblends", str(self.matchrad)]
        args += _bool("radec", self.radec)
        if self.xycol is not None:
            args += ["xycol"] + [str(c) for c in self.xycol]
        if self.starlist is not None:
            args += ["starlist", self.starlist]
        if self.zeromag is not None:
            args += ["zeromag", str(self.zeromag)]
        args += _bool("nofluxconvert", self.nofluxconvert)
        args += _period_spec(self.period)
        args += ["Nharm", str(self.nharm)]
        m_spec = _norm_save(self.save_matches)
        if _should_emit(m_spec):
            args += ["omatches", m_spec.path or outdir]
        return args

    def _resolve_back_references(self, prev) -> None:
        # Resolve "fixcolumn NAME" against the prior step's output.  Numeric
        # values and other specs pass through unchanged.
        from ._helpers import _resolve_period_backref
        self.period = _resolve_period_backref(prev, self.period)

    def _output_file_specs(self):
        return {"matches": (".findblends.matches", None)}
