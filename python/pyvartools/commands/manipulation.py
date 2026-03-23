"""Light curve manipulation and statistics command wrappers."""

from __future__ import annotations
from typing import List, Optional, Union

from pyvartools._command import VartoolsCommand
from ._helpers import _bool, _flag, _injectparam, _norm_save, _outtoken, _period_spec, _pval, _should_emit


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
        args = ["-clip", str(self.sigclip), "1" if self.iterative else "0"]
        if self.niter is not None:
            args += ["niter", str(self.niter)]
        args += _bool("median", self.median)
        if self.markclip is not None:
            args += ["markclip", self.markclip]
            args += _bool("noinitmark", self.noinitmark)
        args += _flag("maskpoints", self.maskpoints)
        return args


class rms(VartoolsCommand):
    """Compute RMS and weighted RMS of the light curve.

    Parameters
    ----------
    maskpoints : str, optional
    """

    _vt_name = "rms"

    def __init__(self, maskpoints: Optional[str] = None) -> None:
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        return ["-rms"] + _flag("maskpoints", self.maskpoints)


class rmsbin(VartoolsCommand):
    """Compute binned RMS at a series of timescales.

    Parameters
    ----------
    nbin : int
        Number of timescale bins.
    bintimes : list of float
        Timescales (e.g. in days) at which to compute the binned RMS.
    maskpoints : str, optional
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
    """Compute chi-squared statistic.

    Parameters
    ----------
    maskpoints : str, optional
    """

    _vt_name = "chi2"

    def __init__(self, maskpoints: Optional[str] = None) -> None:
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        return ["-chi2"] + _flag("maskpoints", self.maskpoints)


class chi2bin(VartoolsCommand):
    """Compute binned chi-squared at multiple timescales.

    Parameters
    ----------
    nbin : int
    bintimes : list of float
    maskpoints : str, optional
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
    """Alarm statistic for transit detection.

    Parameters
    ----------
    maskpoints : str, optional
    """

    _vt_name = "alarm"

    def __init__(self, maskpoints: Optional[str] = None) -> None:
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        return ["-alarm"] + _flag("maskpoints", self.maskpoints)


class rescalesig(VartoolsCommand):
    """Rescale measurement uncertainties to match the scatter.

    Parameters
    ----------
    maskpoints : str, optional
    """

    _vt_name = "rescalesig"

    def __init__(self, maskpoints: Optional[str] = None) -> None:
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        return ["-rescalesig"] + _flag("maskpoints", self.maskpoints)


class ensemblerescalesig(VartoolsCommand):
    """Rescale uncertainties using ensemble sigma-clipping.

    Parameters
    ----------
    sigclip : float
        Sigma-clipping threshold.
    maskpoints : str, optional
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
        Statistic name(s) to compute (e.g. ``"mean"``, ``"median"``,
        ``"stddev"``, ``"min"``, ``"max"``).
    maskpoints : str, optional

    Examples
    --------
    ::

        stats("mag", "mean,median,stddev")
        stats(["mag", "err"], ["mean", "stddev"])
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


class Killharm(VartoolsCommand):
    """Remove harmonic signals from the light curve.

    Parameters
    ----------
    period : float or str
        Period of the signal to remove.  Can be a number, ``"ls"``, ``"aov"``,
        ``"bls"``, etc.
    nharm : int
        Number of harmonics to fit.
    nsubharm : int
        Number of sub-harmonics to fit.
    save_model : bool
        Write the fitted harmonic model to a file.
    fitonly : bool
        Fit the model but do not subtract it from the light curve.
    output_format : str, optional
        Coefficient output format: ``"outampphase"``, ``"outampradphase"``,
        ``"outRphi"``, or ``"outRradphi"``.
    clip : float, optional
        Sigma-clipping threshold to apply after fitting before refitting.
    maskpoints : str, optional
    """

    _vt_name = "Killharm"

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
        args = ["-Killharm"] + self._killharm_period_spec()
        args += [str(self.nharm), str(self.nsubharm)]
        args += _outtoken(self.save_model, outdir)
        args += _bool("fitonly", self.fitonly)
        if self.output_format is not None:
            args += [self.output_format]
        if self.clip is not None:
            args += ["clip", str(self.clip)]
        args += _flag("maskpoints", self.maskpoints)
        return args

    def _killharm_period_spec(self) -> List[str]:
        """Build period spec tokens for Killharm.

        Killharm's "fix" spec has the form: fix Nper per1 ... perN.
        A plain float becomes: fix 1 <value>.
        A string like "fix 2.0 1.0" is expanded to: fix 2 2.0 1.0.
        Keywords like "ls", "aov", "both" are passed as-is.
        """
        p = self.period
        if isinstance(p, (int, float)):
            return ["fix", "1", str(p)]
        tokens = str(p).split()
        if tokens[0] == "fix":
            # bare "fix val" → insert count of periods
            periods = tokens[1:]
            return ["fix", str(len(periods))] + periods
        return tokens

    def _output_file_specs(self):
        return {"model": (".killharm.model", None)}


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
    save_model : bool
        Write the fitted model to a file.
    model_nameformat : str, optional
        Format string for the output model filename.
    fitmask : str, optional
        Name of a mask variable; points where the variable is non-zero are
        excluded from the fit.
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
            args += ["reject", str(self.reject)]
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
        Period of the signal to inject.
    nharm : int
        Number of harmonics.
    amplitude : float
        Signal amplitude (for a single harmonic, semi-amplitude in mag).
    phase : float, optional
        Initial phase (0–1).
    save_model : bool
        Write the injected signal model to a file.

    Notes
    -----
    By design, this class exposes only the most common injection modes:

    * **Amplitude**: ``"ampfix"`` and ``"amplogrand"`` only.  For ``"amprand"`` or
      ``"amplist"`` use :class:`~pyvartools.commands.Raw`.
    * **Period**: ``"fix"`` and ``"logrand"`` (via the ``period`` parameter).  For
      ``"list"`` or ``"rand"`` period modes use :class:`~pyvartools.commands.Raw`.

    These restrictions keep the parameter space manageable.  The full CLI can always
    be accessed via ``cmd.Raw("-Injectharm ...")``.
    """

    _vt_name = "Injectharm"

    def __init__(
        self,
        period,
        amplitude: float,
        nharm: int = 1,
        phase: float = 0.0,
        nsubharm: int = 0,
        save_model: bool = False,
    ) -> None:
        self.period = period
        self.amplitude = amplitude
        self.nharm = nharm
        self.phase = phase
        self.nsubharm = nsubharm
        self.save_model = save_model

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args = ["-Injectharm"] + _period_spec(self.period)
        # vartools Nharm = nharm - 1: Nharm=0 means just the fundamental,
        # Nharm=1 means fundamental + 1st harmonic, etc.
        vt_nharm = self.nharm - 1
        args += [str(vt_nharm)]
        # Repeat amp/phase spec for each of the nharm harmonics
        for _ in range(self.nharm):
            args += ["ampfix", str(self.amplitude), "phasefix", str(self.phase)]
        # sub-harmonics
        args += [str(self.nsubharm)]
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
        Limb-darkening law: ``"quad"`` or ``"nonlin"``.
    ld_coeffs : list of float
        Limb-darkening coefficients.
    save_model : bool
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
        args += [self.ld_type] + [str(c) for c in self.ld_coeffs]
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
    """

    _vt_name = "restoretimes"

    def __init__(self, prior_command: int = 1) -> None:
        self.prior_command = prior_command

    def _to_cli_args(self) -> List[str]:
        return ["-restoretimes", str(self.prior_command)]


class savelc(VartoolsCommand):
    """Save the current light curve state (for later restoration).

    Use :class:`restorelc` to restore.
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
    """Convert differential flux to magnitude.

    Parameters
    ----------
    mag_constant : float
        Reference magnitude (zero-point magnitude).
    offset : float
        Flux offset (added to the flux before conversion).
    magcolumn : int, optional
        Which magnitude column to convert (1-based).
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
        args = ["-difffluxtomag", str(self.mag_constant), str(self.offset)]
        args += _flag("magcolumn", self.magcolumn)
        return args


class fluxtomag(VartoolsCommand):
    """Convert flux to magnitude.

    Parameters
    ----------
    mag_constant : float
        Zero-point magnitude.
    offset : float
        Flux offset.
    """

    _vt_name = "fluxtomag"

    def __init__(self, mag_constant: float, offset: float = 0.0) -> None:
        self.mag_constant = mag_constant
        self.offset = offset

    def _to_cli_args(self) -> List[str]:
        return ["-fluxtomag", str(self.mag_constant), str(self.offset)]


class changeerror(VartoolsCommand):
    """Rescale measurement uncertainties by a constant factor.

    Parameters
    ----------
    maskpoints : str, optional
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
        Target column: ``"t"``, ``"mag"``, ``"err"``, or ``"id"``.
    var : str
        Source variable name.
    """

    _vt_name = "changevariable"

    def __init__(self, column: str, var: str) -> None:
        self.column = column
        self.var = var

    def _to_cli_args(self) -> List[str]:
        return ["-changevariable", self.column, self.var]


class copylc(VartoolsCommand):
    """Duplicate the light curve N times in the output.

    Parameters
    ----------
    ncopies : int
    """

    _vt_name = "copylc"

    def __init__(self, ncopies: int) -> None:
        self.ncopies = ncopies

    def _to_cli_args(self) -> List[str]:
        return ["-copylc", str(self.ncopies)]


class medianfilter(VartoolsCommand):
    """Median (or mean) filter the light curve.

    Parameters
    ----------
    time : float
        Filter timescale.
    method : str
        ``"median"`` (default), ``"average"``, or ``"weightedaverage"``.
    replace : bool
        Replace the magnitude values with the filtered version.
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
        args = ["-medianfilter", str(self.time)]
        if self.method in ("average", "weightedaverage"):
            args.append(self.method)
        args += _bool("replace", self.replace)
        return args


class expr(VartoolsCommand):
    """Evaluate an analytic expression to create or update a variable.

    Parameters
    ----------
    expression : str
        Expression of the form ``varname=expr``.
    outputcolumns : str, optional
        Comma-separated list of column names to output.
    """

    _vt_name = "expr"

    def __init__(
        self,
        expression: str,
        outputcolumns: Optional[str] = None,
    ) -> None:
        self.expression = expression
        self.outputcolumns = outputcolumns

    def _to_cli_args(self) -> List[str]:
        args = ["-expr", self.expression]
        args += _flag("outputcolumns", self.outputcolumns)
        return args


class print_cols(VartoolsCommand):
    """Print selected variables to stdout as additional columns.

    Parameters
    ----------
    variables : str or list of str
        Variable name(s) to print.
    columnnames : str or list of str, optional
        Output column header names.
    fmt : str or list of str, optional
        printf-style format strings for each column.
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
    """Compute the Fast Fourier Transform of two variables.

    Parameters
    ----------
    input_real, input_imag : str
        Input real and imaginary variable names.
    output_real, output_imag : str
        Output real and imaginary variable names.
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
    """Inverse FFT (same parameter structure as FFT)."""

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
        if self.file_times is not None:
            toks = str(self.file_times).split()
            if toks[0] == "list":
                args += ["file"] + toks
            else:
                # treat as a path → "fix" mode
                args += ["file", "fix", str(self.file_times)]
                if self.file_column is not None:
                    args += ["column", str(self.file_column)]
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
    save_model : bool
    maskpoints : str, optional
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
        m_spec = _norm_save(self.save_model)
        if _should_emit(m_spec):
            args += ["omodel", m_spec.path or outdir]
        else:
            args += ["0"]
        args += _flag("maskpoints", self.maskpoints)
        return args

    def _output_file_specs(self):
        return {"model": (".decorr.model", None)}


class Jstet(VartoolsCommand):
    """J-statistic (Stetson variability index).

    Parameters
    ----------
    timescale : float
        Timescale for variability index calculation.
    dates : str
        Path to the dates file.
    maskpoints : str, optional
    """

    _vt_name = "Jstet"

    def __init__(
        self,
        timescale: float,
        dates: str,
        maskpoints: Optional[str] = None,
    ) -> None:
        self.timescale = timescale
        self.dates = dates
        self.maskpoints = maskpoints

    def _to_cli_args(self) -> List[str]:
        return (["-Jstet", str(self.timescale), self.dates]
                + _flag("maskpoints", self.maskpoints))
