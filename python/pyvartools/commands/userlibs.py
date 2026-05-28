"""Typed Python wrappers for the USERLIB extension commands shipped with vartools.

These wrappers turn each extension (which is a dynamically-loaded ``.so`` /
``.la`` file under ``USERLIBS/src/``) into a first-class pyvartools command
with typed parameters, rather than requiring users to assemble raw argument
strings via :class:`~pyvartools.UserCommand`.

Each wrapper accepts an optional ``lib_path`` parameter; when given, an
``-L <lib_path>`` flag is prepended to the CLI invocation so the extension
can be located.  When omitted, vartools is expected to auto-load the
extension from its installed userlibs directory.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import List, Optional, Union

from pyvartools._command import VartoolsCommand
from ._helpers import _flag, _bool, _norm_save, _should_emit, _outtoken


def _extparam(val, vary: Optional[bool] = None) -> List[str]:
    """Build tokens for a ``<"fix" v | "list" | "fixcolumn" c | "expr" e>`` spec.

    - ``None`` → ``[]``
    - number → ``["fix", str(val)]``
    - string starting with ``fix/list/fixcolumn/expr`` → ``str.split()``
    - bare identifier or expression → ``["expr", val]``

    When *vary* is True the token ``"vary"`` is appended.
    """
    if val is None:
        return []
    if isinstance(val, (int, float)):
        out = ["fix", str(val)]
    else:
        s = str(val).strip()
        if re.match(r'^(fix|list|fixcolumn|expr)\b', s):
            out = s.split()
        else:
            out = ["expr", s]
    if vary:
        out = out + ["vary"]
    return out


class _UserLibCommand(VartoolsCommand):
    """Base class for typed USERLIB extension wrappers.

    Subclasses set ``_vt_name`` and implement ``_to_cli_args``.  If an instance
    attribute ``lib_path`` is truthy, the tokens ``["-L", lib_path]`` are
    prepended to ``_to_cli_args()``.
    """

    _vt_name = ""

    def _libprefix(self) -> List[str]:
        lp = getattr(self, "lib_path", None)
        return ["-L", str(lp)] if lp else []

    def _output_file_specs(self) -> dict:
        return {}


# -----------------------------------------------------------------------------
# magadd — simplest case: a single fix/list/fixcolumn/expr scalar.
# -----------------------------------------------------------------------------

class magadd(_UserLibCommand):
    """Add a constant offset to light-curve magnitudes (USERLIB ``-magadd``).

    Parameters
    ----------
    value : float or str
        Constant to add.  Number → ``fix value``; string is parsed as
        ``"list [column N]"``, ``"fixcolumn NAME"``, ``"expr EXPR"``, or a
        bare expression.
    lib_path : str, optional
        Path to ``magadd.so`` / ``magadd.la``.  Omit to rely on vartools
        auto-loading from the installed userlibs directory.

    See Also
    --------
    USERLIB extension command: ``-magadd``.  Shipped as the template /
    example extension; useful as a starting point for writing your own.
    """

    _vt_name = "magadd"

    def __init__(self, value, lib_path: Optional[str] = None) -> None:
        self.value = value
        self.lib_path = lib_path

    def _to_cli_args(self) -> List[str]:
        return self._libprefix() + ["-magadd"] + _extparam(self.value)


# -----------------------------------------------------------------------------
# hatpiflag — five variable names.
# -----------------------------------------------------------------------------

class hatpiflag(_UserLibCommand):
    """Combine HATPI quality flags into a single binary flag (USERLIB ``-hatpiflag``).

    Parameters
    ----------
    fiphot_string_flag_var : str
        Name of the vector variable holding fiphot string flags.
    rejbadframe_mask_var : str
        Mask variable: 0 = rejected frame, 1 = good.
    tfa_outlier_mask_var : str
        TFA outlier mask: 0 = outlier, 1 = not outlier.
    pointing_outlier_flag_var : str
        Pointing outlier flag: 1 = outlier, 0 = not an outlier.
    output_flag_var : str
        Name of the output binary flag variable to create.
    lib_path : str, optional
        Path to ``hatpiflag.so``.

    See Also
    --------
    USERLIB extension command: ``-hatpiflag``.
    """

    _vt_name = "hatpiflag"

    def __init__(
        self,
        fiphot_string_flag_var: str,
        rejbadframe_mask_var: str,
        tfa_outlier_mask_var: str,
        pointing_outlier_flag_var: str,
        output_flag_var: str,
        lib_path: Optional[str] = None,
    ) -> None:
        self.fiphot_string_flag_var = fiphot_string_flag_var
        self.rejbadframe_mask_var = rejbadframe_mask_var
        self.tfa_outlier_mask_var = tfa_outlier_mask_var
        self.pointing_outlier_flag_var = pointing_outlier_flag_var
        self.output_flag_var = output_flag_var
        self.lib_path = lib_path

    def _to_cli_args(self) -> List[str]:
        return self._libprefix() + [
            "-hatpiflag",
            self.fiphot_string_flag_var,
            self.rejbadframe_mask_var,
            self.tfa_outlier_mask_var,
            self.pointing_outlier_flag_var,
            self.output_flag_var,
        ]


# -----------------------------------------------------------------------------
# fastchi2 — Palmer (2009) fast chi2 periodogram.
# -----------------------------------------------------------------------------

class fastchi2(_UserLibCommand):
    """Palmer (2009) fast chi2 periodogram (USERLIB ``-fastchi2``).

    Parameters
    ----------
    Nharm : float or str
        Number of harmonics in the model.
    freqmax : float or str
        Maximum frequency to search (cycles/day).
    freqmin : float or str, optional
        Minimum frequency to search (default 0).
    detrendorder : int or str, optional
        Polynomial order for pre-detrending.
    t0, timespan, oversample, chimargin : float or str, optional
        Passed through to the ``fix/list/fixcolumn/expr`` parser.
    Npeak : int, optional
        Number of peaks to output (constant integer).
    norefitpeak : bool
        Skip the fine peak search.
    save_per : bool | str | Output, optional
        Output directory for periodograms (``oper`` keyword).
    save_model : bool | str | Output, optional
        Output directory for harmonic-function model light curves (``omodel``).
    omodelvariable : str, optional
        Name of a light-curve variable to write the model into.
    lib_path : str, optional
        Path to ``fastchi2.so``.

    See Also
    --------
    USERLIB extension command: ``-fastchi2``.
    Citation: Palmer 2009 (ApJ 695, 496).
    """

    _vt_name = "fastchi2"

    def __init__(
        self,
        Nharm,
        freqmax,
        freqmin=None,
        detrendorder=None,
        t0=None,
        timespan=None,
        oversample=None,
        chimargin=None,
        Npeak: Optional[int] = None,
        norefitpeak: bool = False,
        save_per=False,
        save_model=False,
        omodelvariable: Optional[str] = None,
        lib_path: Optional[str] = None,
    ) -> None:
        self.Nharm = Nharm
        self.freqmax = freqmax
        self.freqmin = freqmin
        self.detrendorder = detrendorder
        self.t0 = t0
        self.timespan = timespan
        self.oversample = oversample
        self.chimargin = chimargin
        self.Npeak = Npeak
        self.norefitpeak = norefitpeak
        self.save_per = save_per
        self.save_model = save_model
        self.omodelvariable = omodelvariable
        self.lib_path = lib_path

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args: List[str] = self._libprefix() + ["-fastchi2"]
        args += ["Nharm"] + _extparam(self.Nharm)
        args += ["freqmax"] + _extparam(self.freqmax)
        if self.freqmin is not None:
            args += ["freqmin"] + _extparam(self.freqmin)
        if self.detrendorder is not None:
            args += ["detrendorder"] + _extparam(self.detrendorder)
        if self.t0 is not None:
            args += ["t0"] + _extparam(self.t0)
        if self.timespan is not None:
            args += ["timespan"] + _extparam(self.timespan)
        if self.oversample is not None:
            args += ["oversample"] + _extparam(self.oversample)
        if self.chimargin is not None:
            args += ["chimargin"] + _extparam(self.chimargin)
        if self.Npeak is not None:
            args += ["Npeak", str(self.Npeak)]
        if self.norefitpeak:
            args += ["norefitpeak"]
        per_spec = _norm_save(self.save_per)
        if _should_emit(per_spec):
            args += ["oper", per_spec.path or outdir]
        model_spec = _norm_save(self.save_model)
        if _should_emit(model_spec):
            args += ["omodel", model_spec.path or outdir]
        if self.omodelvariable is not None:
            args += ["omodelvariable", self.omodelvariable]
        return args

    def _output_file_specs(self) -> dict:
        return {
            "per": (".fastchi2_per", None),
            "model": (".fastchi2_model", None),
        }


# -----------------------------------------------------------------------------
# splinedetrend — multivariate spline/poly/harm detrending.
# -----------------------------------------------------------------------------

class splinedetrend(_UserLibCommand):
    """Basis-spline / polynomial / harmonic detrending (USERLIB ``-splinedetrend``).

    Parameters
    ----------
    detrendvecs : str or list of str
        Detrending specs.  Either a single comma-joined string
        (e.g. ``"t:spline:0.1:3,x:poly:2"``) or a list of individual specs
        which will be joined with commas.  Each spec has the form
        ``VAR:<spline:knotspacing:order | poly:order | harm:nharm>[:groupbygap:gapsize]``.
    sigmaclip : float or str, optional
        Optional sigma-clipping threshold
        (``fix/list/fixcolumn/expr`` scalar).
    save_model : bool | str | Output, optional
        Output directory for best-fit models (``omodel``).
    save_coeffs : bool | str | Output, optional
        Output directory for model coefficients (``omodelcoeffs``).
    omodelvariable : str, optional
        Comma-separated ``outvar[:inputvarsignal]`` specs to store per-LC models.
    lib_path : str, optional
        Path to ``splinedetrend.so``.

    See Also
    --------
    USERLIB extension command: ``-splinedetrend``.
    """

    _vt_name = "splinedetrend"

    def __init__(
        self,
        detrendvecs: Union[str, List[str]],
        sigmaclip=None,
        save_model=False,
        save_coeffs=False,
        omodelvariable: Optional[str] = None,
        lib_path: Optional[str] = None,
    ) -> None:
        self.detrendvecs = detrendvecs
        self.sigmaclip = sigmaclip
        self.save_model = save_model
        self.save_coeffs = save_coeffs
        self.omodelvariable = omodelvariable
        self.lib_path = lib_path

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        if isinstance(self.detrendvecs, (list, tuple)):
            spec = ",".join(str(s) for s in self.detrendvecs)
        else:
            spec = str(self.detrendvecs)
        args: List[str] = self._libprefix() + ["-splinedetrend", spec]
        if self.sigmaclip is not None:
            args += ["sigmaclip"] + _extparam(self.sigmaclip)
        m_spec = _norm_save(self.save_model)
        if _should_emit(m_spec):
            args += ["omodel", m_spec.path or outdir]
        c_spec = _norm_save(self.save_coeffs)
        if _should_emit(c_spec):
            args += ["omodelcoeffs", c_spec.path or outdir]
        if self.omodelvariable is not None:
            args += ["omodelvariable", self.omodelvariable]
        return args

    def _output_file_specs(self) -> dict:
        return {
            "model": (".splinedetrend_model", None),
            "coeffs": (".splinedetrend_modelcoeffs", None),
        }


# -----------------------------------------------------------------------------
# ftuneven — complex Fourier transform of unevenly sampled data.
# -----------------------------------------------------------------------------

class ftuneven(_UserLibCommand):
    """Complex Fourier transform of unevenly-sampled data (USERLIB ``-ftuneven``).

    At least one **output mode** and one **frequency source** must be given.

    Output modes (supply exactly one):
      - ``output_vectors=(freq, ft_real, ft_imag, periodogram)`` — 4 tuple/list
        of vector variable names into which the FT is written.
      - ``output_file=True`` (or a path) — write the FT to a per-LC file in
        ``save_outdir`` (or the pipeline temp dir).
      - Both — pass ``output_vectors=(...)`` **and** ``output_file=True``.

    Frequency source (supply exactly one):
      - ``freqauto=True`` — choose the grid automatically.
      - ``freqrange=(minfreq, maxfreq, freqstep)`` — each element is passed
        through ``_extparam`` (so numbers or strings are both accepted).
      - ``freqvariable="name"`` — read the frequency grid from a variable.
      - ``freqfile="/path/to/file"`` — read the frequency grid from a file.

    Parameters
    ----------
    output_vectors : tuple of 4 str, optional
    output_file : bool | str, optional
    save_outdir : str, optional
        Used when ``output_file`` is True or ``output_vectors`` is combined
        with ``output_file`` — overrides the pipeline temp dir.
    nameformat : str, optional
        Format string for the per-LC output filename.
    freqauto : bool
    freqrange : tuple (min, max, step), optional
    freqvariable : str, optional
    freqfile : str, optional
    ft_sign : int, optional
    tt_zero : float, optional
    changeinputvectors : tuple of 3 str, optional
        ``(tvec, data_real_vec, data_imag_vec)``.
    lib_path : str, optional

    See Also
    --------
    USERLIB extension command: ``-ftuneven``.
    Citation: Scargle 1989 (ApJ 343, 874) for the unevenly-sampled
    Fourier-transform method.
    """

    _vt_name = "ftuneven"

    def __init__(
        self,
        output_vectors: Optional[tuple] = None,
        output_file: Union[bool, str] = False,
        save_outdir: Optional[str] = None,
        nameformat: Optional[str] = None,
        freqauto: bool = False,
        freqrange: Optional[tuple] = None,
        freqvariable: Optional[str] = None,
        freqfile: Optional[str] = None,
        ft_sign: Optional[int] = None,
        tt_zero: Optional[float] = None,
        changeinputvectors: Optional[tuple] = None,
        lib_path: Optional[str] = None,
    ) -> None:
        self.output_vectors = output_vectors
        self.output_file = output_file
        self.save_outdir = save_outdir
        self.nameformat = nameformat
        self.freqauto = freqauto
        self.freqrange = freqrange
        self.freqvariable = freqvariable
        self.freqfile = freqfile
        self.ft_sign = ft_sign
        self.tt_zero = tt_zero
        self.changeinputvectors = changeinputvectors
        self.lib_path = lib_path

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args: List[str] = self._libprefix() + ["-ftuneven"]

        # Output mode
        want_file = bool(self.output_file) or isinstance(self.output_file, str)
        want_vecs = self.output_vectors is not None
        if not (want_file or want_vecs):
            raise ValueError(
                "ftuneven: must specify output_vectors or output_file"
            )
        if want_vecs and len(self.output_vectors) != 4:
            raise ValueError(
                "ftuneven: output_vectors must be a 4-tuple"
                " (freq, ft_real, ft_imag, periodogram)"
            )

        filedir = (self.save_outdir
                   or (self.output_file if isinstance(self.output_file, str)
                       else outdir))

        if want_vecs and want_file:
            args += ["outputvectorsandfile", filedir]
            if self.nameformat is not None:
                args += ["nameformat", self.nameformat]
            args += list(self.output_vectors)
        elif want_vecs:
            args += ["outputvectors"] + list(self.output_vectors)
        else:
            args += ["outputfile", filedir]
            if self.nameformat is not None:
                args += ["nameformat", self.nameformat]

        # Frequency source (exactly one)
        freqspecs = [self.freqauto, self.freqrange is not None,
                     self.freqvariable is not None, self.freqfile is not None]
        if sum(bool(x) for x in freqspecs) != 1:
            raise ValueError(
                "ftuneven: must specify exactly one of freqauto, freqrange,"
                " freqvariable, freqfile"
            )
        if self.freqauto:
            args += ["freqauto"]
        elif self.freqrange is not None:
            mn, mx, st = self.freqrange
            args += ["freqrange"]
            args += ["minfreq"] + _extparam(mn)
            args += ["maxfreq"] + _extparam(mx)
            args += ["freqstep"] + _extparam(st)
        elif self.freqvariable is not None:
            args += ["freqvariable", self.freqvariable]
        else:
            args += ["freqfile", self.freqfile]

        if self.ft_sign is not None:
            args += ["ft_sign", str(self.ft_sign)]
        if self.tt_zero is not None:
            args += ["tt_zero", str(self.tt_zero)]
        if self.changeinputvectors is not None:
            if len(self.changeinputvectors) != 3:
                raise ValueError(
                    "ftuneven: changeinputvectors must be a 3-tuple"
                    " (tvec, data_real_vec, data_imag_vec)"
                )
            args += ["changeinputvectors"] + list(self.changeinputvectors)
        return args

    def _output_file_specs(self) -> dict:
        return {"ft": (".ftuneven", None)}


# -----------------------------------------------------------------------------
# stitch — fit for and remove offsets between light-curve segments.
# -----------------------------------------------------------------------------

_STITCH_HEADER_RE = re.compile(r"^#\s*Parameters for stitch variable\s+(\d+)\s*$")
_STITCH_SHIFT_RE  = re.compile(r"^LCgroup_(\d+)\s+shift:\s*(\S+)\s*$")
# Coefficient line, two flavours:
#   "Coeff for TERM, TMIN<t<TMAX: VALUE"   (poly, harmseries)
#   "Median TMIN<t<TMAX: VALUE"            (median / mean / weightedmean)
_STITCH_COEFF_RE = re.compile(
    r"^(?:Coeff for\s+(?P<term1>.+?),\s*"
    r"|(?P<term2>Median|Mean|WeightedMean)\s+)"
    r"(?P<tmin>\S+)<t<(?P<tmax>\S+):\s*(?P<value>\S+)\s*$"
)


def _parse_stitch_fitted_params(path):
    """Parse the structured text file written by ``-stitch save_fitted_parameters``.

    Returns a DataFrame with columns ``variable, kind, term, t_min, t_max,
    value``.  ``kind`` is ``"coeff"`` (one row per fitted basis-function
    coefficient and time bin) or ``"shift"`` (one row per LC group's
    additive offset).  ``t_min`` / ``t_max`` are NaN for shift rows.
    """
    import math
    import pandas as pd

    rows = []
    variable = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            m = _STITCH_HEADER_RE.match(line)
            if m:
                variable = int(m.group(1))
                continue
            m = _STITCH_SHIFT_RE.match(line)
            if m:
                rows.append(dict(variable=variable, kind="shift",
                                 term=f"LCgroup_{m.group(1)}",
                                 t_min=math.nan, t_max=math.nan,
                                 value=float(m.group(2))))
                continue
            m = _STITCH_COEFF_RE.match(line)
            if m:
                term = m.group("term1") or m.group("term2")
                rows.append(dict(variable=variable, kind="coeff",
                                 term=term,
                                 t_min=float(m.group("tmin")),
                                 t_max=float(m.group("tmax")),
                                 value=float(m.group("value"))))
    return pd.DataFrame(
        rows, columns=["variable", "kind", "term", "t_min", "t_max", "value"]
    )


class stitch(_UserLibCommand):
    """Stitch multi-segment light curves at offsets (USERLIB ``-stitch``).

    Parameters
    ----------
    stitch_variables : str or list of str
        Variable(s) to apply the stitch procedure to (typically ``"mag"``).
    uncertainty_variables : str or list of str
        Uncertainty variable(s) paired with *stitch_variables*.
    mask_variables : str or list of str
        Mask variable(s) — points with mask > 0 are excluded from fitting.
    lcnum_var : str
        Variable identifying the light-curve segment for each point.
    method : str
        One of ``"median"``, ``"mean"``, ``"weightedmean"``,
        ``"poly N"`` (where ``N`` is the polynomial order), or
        ``"harmseries PERIODVAR NHARM"``.
    refnum_var : str, optional
        Additional reference-number variable for finer segmentation.
    groupbytime : float, optional
        Group segments into bins of this width in the LC time units.
    groupbytime_start : float, optional
        Start time for the first time bin (only meaningful when
        *groupbytime* is set).
    fitonly : bool
        Fit shifts but do not subtract them.
    save_fitted_parameters : bool | str | Output, optional
        Output directory for per-source fitted-parameter files.
    fitted_parameters_nameformat : str, optional
        ``format`` string applied to the fitted-parameter filenames.
    add_stitchparams_fitsheader : bool or str, optional
        ``True`` or ``"primary"`` / ``"extension"``.
    add_stitchparams_mode : str, optional
        ``"append"`` or ``"update"``.
    add_shifts_fitsheader : str, optional
        Keyword base (e.g. ``"SHFT"``) to log shifts into FITS headers.
    add_shifts_hdu : str, optional
        ``"primary"`` or ``"extension"``.
    add_shifts_mode : str, optional
        ``"append"`` or ``"update"``.
    shifts_file : tuple of 2 str, optional
        ``(fieldlabelsvar, starnamevar)`` — enables ``shifts_file`` mode.
    append_refnum_to_fieldlabel : bool
    in_shifts_file : str or list of str, optional
        Pre-existing shifts file(s) to read.  One file per stitched
        variable: pass a single string when *stitch_variables* is a
        string, or a list of the same length as *stitch_variables* when
        it is a list.
    nobs_refit : int, optional
    header_basename_only : bool
    out_shifts_file : str or list of str, optional
        Output shifts file(s) to write.  One file per stitched variable
        (same shape rule as *in_shifts_file*).
    include_missing : bool
    lib_path : str, optional

    See Also
    --------
    USERLIB extension command: ``-stitch``.  Commonly used after
    ``run_combinelcs`` to merge multi-telescope segments.
    """

    _vt_name = "stitch"

    def __init__(
        self,
        stitch_variables: Union[str, List[str]],
        uncertainty_variables: Union[str, List[str]],
        mask_variables: Union[str, List[str]],
        lcnum_var: str,
        method: str,
        refnum_var: Optional[str] = None,
        groupbytime: Optional[float] = None,
        groupbytime_start: Optional[float] = None,
        fitonly: bool = False,
        save_fitted_parameters=False,
        fitted_parameters_nameformat: Optional[str] = None,
        add_stitchparams_fitsheader: Union[bool, str] = False,
        add_stitchparams_mode: Optional[str] = None,
        add_shifts_fitsheader: Optional[str] = None,
        add_shifts_hdu: Optional[str] = None,
        add_shifts_mode: Optional[str] = None,
        shifts_file: Optional[tuple] = None,
        append_refnum_to_fieldlabel: bool = False,
        in_shifts_file: Union[str, List[str], None] = None,
        nobs_refit: Optional[int] = None,
        header_basename_only: bool = False,
        out_shifts_file: Union[str, List[str], None] = None,
        include_missing: bool = False,
        lib_path: Optional[str] = None,
    ) -> None:
        self.stitch_variables = stitch_variables
        self.uncertainty_variables = uncertainty_variables
        self.mask_variables = mask_variables
        self.lcnum_var = lcnum_var
        self.method = method
        self.refnum_var = refnum_var
        self.groupbytime = groupbytime
        self.groupbytime_start = groupbytime_start
        self.fitonly = fitonly
        self.save_fitted_parameters = save_fitted_parameters
        self.fitted_parameters_nameformat = fitted_parameters_nameformat
        self.add_stitchparams_fitsheader = add_stitchparams_fitsheader
        self.add_stitchparams_mode = add_stitchparams_mode
        self.add_shifts_fitsheader = add_shifts_fitsheader
        self.add_shifts_hdu = add_shifts_hdu
        self.add_shifts_mode = add_shifts_mode
        self.shifts_file = shifts_file
        self.append_refnum_to_fieldlabel = append_refnum_to_fieldlabel
        self.in_shifts_file = in_shifts_file
        self.nobs_refit = nobs_refit
        self.header_basename_only = header_basename_only
        self.out_shifts_file = out_shifts_file
        self.include_missing = include_missing
        self.lib_path = lib_path

    @staticmethod
    def _joinlist(v) -> str:
        if isinstance(v, (list, tuple)):
            return ",".join(str(x) for x in v)
        return str(v)

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args: List[str] = self._libprefix() + ["-stitch"]
        args += [
            self._joinlist(self.stitch_variables),
            self._joinlist(self.uncertainty_variables),
            self._joinlist(self.mask_variables),
            self.lcnum_var,
        ]
        if self.refnum_var is not None:
            args += ["refnum_var", self.refnum_var]
        args += str(self.method).split()
        if self.groupbytime is not None:
            args += ["groupbytime", str(self.groupbytime)]
            if self.groupbytime_start is not None:
                args += ["start", str(self.groupbytime_start)]
        if self.fitonly:
            args += ["fitonly"]
        params_spec = _norm_save(self.save_fitted_parameters)
        if _should_emit(params_spec):
            args += ["save_fitted_parameters", params_spec.path or outdir]
            if self.fitted_parameters_nameformat is not None:
                args += ["format", self.fitted_parameters_nameformat]
        if self.add_stitchparams_fitsheader:
            args += ["add_stitchparams_fitsheader"]
            if isinstance(self.add_stitchparams_fitsheader, str):
                args += [self.add_stitchparams_fitsheader]
            if self.add_stitchparams_mode is not None:
                args += [self.add_stitchparams_mode]
        if self.add_shifts_fitsheader is not None:
            args += ["add_shifts_fitsheader", self.add_shifts_fitsheader]
            if self.add_shifts_hdu is not None:
                args += [self.add_shifts_hdu]
            if self.add_shifts_mode is not None:
                args += [self.add_shifts_mode]
        if self.shifts_file is not None:
            if len(self.shifts_file) != 2:
                raise ValueError(
                    "stitch: shifts_file must be (fieldlabelsvar, starnamevar)"
                )
            args += ["shifts_file"] + list(self.shifts_file)
            if self.append_refnum_to_fieldlabel:
                args += ["append_refnum_to_fieldlabel"]
            if self.in_shifts_file is not None:
                args += ["in_shifts_file", self._joinlist(self.in_shifts_file)]
                if self.nobs_refit is not None:
                    args += ["nobs_refit", str(self.nobs_refit)]
                if self.header_basename_only:
                    args += ["header_basename_only"]
            if self.out_shifts_file is not None:
                args += ["out_shifts_file",
                         self._joinlist(self.out_shifts_file)]
                if self.include_missing:
                    args += ["include_missing"]
        return args

    def _output_file_specs(self) -> dict:
        return {"fitted_parameters": (".stitch", _parse_stitch_fitted_params)}


# -----------------------------------------------------------------------------
# jktebop — detached-eclipsing-binary model.
# -----------------------------------------------------------------------------

class jktebop(_UserLibCommand):
    """JKTEBOP detached eclipsing-binary model (USERLIB ``-jktebop``).

    Each of the main parameters (``Period``, ``T0``, ``r1_r2``, ``r2_r1``,
    ``M2_M1``, ``J2_J1``, ``i`` or ``bimpact``, ``esinomega``, ``ecosomega``)
    is a value-spec accepted by :func:`_extparam` (number → ``fix``;
    ``"fix v"`` / ``"list"`` / ``"fixcolumn NAME"`` / ``"expr EXPR"``).
    Pass the corresponding ``vary_*`` boolean to free the parameter in the fit.

    LD1 and LD2 take a ``law`` (``linear/quad/log/sqrt``, or ``"lockLD1"``
    for LD2 to copy LD1) plus 1 or 2 coefficients.

    Parameters
    ----------
    mode : str
        ``"inject"`` or ``"fit"``.
    Period, T0, r1_r2, r2_r1, M2_M1, J2_J1 : value-spec
    i, bimpact : value-spec, exactly one required
    esinomega, ecosomega : value-spec
    vary_* : bool
        One flag per corresponding parameter.
    LD1_law, LD2_law : str
        ``"linear"``, ``"quad"``, ``"log"``, ``"sqrt"``, or (LD2 only)
        ``"lockLD1"``.
    LD1_coeffs, LD2_coeffs : tuple of value-spec
        Single value for ``linear``, two values for the other laws.
        (Empty for ``LD2_law="lockLD1"``.)
    vary_LD1, vary_LD2 : bool
    gravdark1, gravdark2 : value-spec, optional
    vary_gravdark1, vary_gravdark2 : bool
    reflection1, reflection2 : value-spec, optional
    vary_reflection1, vary_reflection2 : bool
    L3, tidallag : value-spec, optional
    vary_L3, vary_tidallag : bool
    correctlc : bool
    save_model : bool | str | Output, optional
        ``omodel`` output directory.
    model_nameformat : str, optional
    save_curve : bool | str | Output, optional
        ``ocurve`` output directory.
    curve_xaxis : str, optional
        ``"jd"`` or ``"phase"``.
    curve_step : float, optional
    curve_nameformat : str, optional
    lib_path : str, optional

    See Also
    --------
    USERLIB extension command: ``-jktebop``.
    Citations: Southworth, Maxted & Smalley 2004 (MNRAS 351, 1277) and
    Nelson & Davis 1972 (ApJ 174, 617) for the underlying EBOP model.
    """

    _vt_name = "jktebop"

    _MAIN_PARAMS = [
        ("Period",    "Period"),
        ("T0",        "T0"),
        ("r1_r2",     "r1+r2"),
        ("r2_r1",     "r2/r1"),
        ("M2_M1",     "M2/M1"),
        ("J2_J1",     "J2/J1"),
    ]

    def __init__(
        self,
        mode: str,
        Period=None, vary_Period: bool = False,
        T0=None, vary_T0: bool = False,
        r1_r2=None, vary_r1_r2: bool = False,
        r2_r1=None, vary_r2_r1: bool = False,
        M2_M1=None, vary_M2_M1: bool = False,
        J2_J1=None, vary_J2_J1: bool = False,
        i=None, vary_i: bool = False,
        bimpact=None, vary_bimpact: bool = False,
        esinomega=None, vary_esinomega: bool = False,
        ecosomega=None, vary_ecosomega: bool = False,
        LD1_law: str = "quad", LD1_coeffs=(0.3, 0.3), vary_LD1: bool = False,
        LD2_law: str = "lockLD1", LD2_coeffs=(), vary_LD2: bool = False,
        gravdark1=None, vary_gravdark1: bool = False,
        gravdark2=None, vary_gravdark2: bool = False,
        reflection1=None, vary_reflection1: bool = False,
        reflection2=None, vary_reflection2: bool = False,
        L3=None, vary_L3: bool = False,
        tidallag=None, vary_tidallag: bool = False,
        correctlc: bool = False,
        save_model=False, model_nameformat: Optional[str] = None,
        save_curve=False,
        curve_xaxis: Optional[str] = None,
        curve_step: Optional[float] = None,
        curve_nameformat: Optional[str] = None,
        lib_path: Optional[str] = None,
    ) -> None:
        self.mode = mode
        self.Period = Period; self.vary_Period = vary_Period
        self.T0 = T0; self.vary_T0 = vary_T0
        self.r1_r2 = r1_r2; self.vary_r1_r2 = vary_r1_r2
        self.r2_r1 = r2_r1; self.vary_r2_r1 = vary_r2_r1
        self.M2_M1 = M2_M1; self.vary_M2_M1 = vary_M2_M1
        self.J2_J1 = J2_J1; self.vary_J2_J1 = vary_J2_J1
        self.i = i; self.vary_i = vary_i
        self.bimpact = bimpact; self.vary_bimpact = vary_bimpact
        self.esinomega = esinomega; self.vary_esinomega = vary_esinomega
        self.ecosomega = ecosomega; self.vary_ecosomega = vary_ecosomega
        self.LD1_law = LD1_law
        self.LD1_coeffs = LD1_coeffs
        self.vary_LD1 = vary_LD1
        self.LD2_law = LD2_law
        self.LD2_coeffs = LD2_coeffs
        self.vary_LD2 = vary_LD2
        self.gravdark1 = gravdark1; self.vary_gravdark1 = vary_gravdark1
        self.gravdark2 = gravdark2; self.vary_gravdark2 = vary_gravdark2
        self.reflection1 = reflection1
        self.vary_reflection1 = vary_reflection1
        self.reflection2 = reflection2
        self.vary_reflection2 = vary_reflection2
        self.L3 = L3; self.vary_L3 = vary_L3
        self.tidallag = tidallag; self.vary_tidallag = vary_tidallag
        self.correctlc = correctlc
        self.save_model = save_model
        self.model_nameformat = model_nameformat
        self.save_curve = save_curve
        self.curve_xaxis = curve_xaxis
        self.curve_step = curve_step
        self.curve_nameformat = curve_nameformat
        self.lib_path = lib_path

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args: List[str] = self._libprefix() + ["-jktebop", self.mode]

        for pyname, cliname in self._MAIN_PARAMS:
            val = getattr(self, pyname)
            vary = getattr(self, f"vary_{pyname}")
            args += [cliname] + _extparam(val, vary=vary)

        # i or bimpact — exactly one
        if self.i is not None and self.bimpact is not None:
            raise ValueError("jktebop: specify exactly one of i or bimpact")
        if self.bimpact is not None:
            args += ["bimpact"] + _extparam(self.bimpact, vary=self.vary_bimpact)
        else:
            args += ["i"] + _extparam(self.i, vary=self.vary_i)

        args += ["esinomega"] + _extparam(self.esinomega, vary=self.vary_esinomega)
        args += ["ecosomega"] + _extparam(self.ecosomega, vary=self.vary_ecosomega)

        # LD1
        args += ["LD1", self.LD1_law]
        coeffs = self.LD1_coeffs if isinstance(self.LD1_coeffs, (list, tuple)) else [self.LD1_coeffs]
        if len(coeffs) == 0:
            raise ValueError("jktebop: LD1_coeffs must contain at least 1 coefficient")
        args += ["fix"] + [str(c) for c in coeffs]
        if self.vary_LD1:
            args += ["vary"]

        # LD2
        args += ["LD2", self.LD2_law]
        if self.LD2_law != "lockLD1":
            ld2c = self.LD2_coeffs if isinstance(self.LD2_coeffs, (list, tuple)) else [self.LD2_coeffs]
            if len(ld2c) == 0:
                raise ValueError(
                    "jktebop: LD2_coeffs required unless LD2_law='lockLD1'"
                )
            args += ["fix"] + [str(c) for c in ld2c]
            if self.vary_LD2:
                args += ["vary"]

        for pyname, cliname in [
            ("gravdark1", "gravdark1"), ("gravdark2", "gravdark2"),
            ("reflection1", "reflection1"), ("reflection2", "reflection2"),
            ("L3", "L3"), ("tidallag", "tidallag"),
        ]:
            val = getattr(self, pyname)
            if val is None:
                continue
            vary = getattr(self, f"vary_{pyname}")
            args += [cliname] + _extparam(val, vary=vary)

        if self.correctlc:
            args += ["correctlc"]
        model_spec = _norm_save(self.save_model)
        if _should_emit(model_spec):
            args += ["omodel", model_spec.path or outdir]
            if self.model_nameformat is not None:
                args += ["format", self.model_nameformat]
        curve_spec = _norm_save(self.save_curve)
        if _should_emit(curve_spec):
            args += ["ocurve"]
            if self.curve_xaxis is not None:
                args += [self.curve_xaxis]
            if self.curve_step is not None:
                args += ["step", str(self.curve_step)]
            args += ["outdir", curve_spec.path or outdir]
            if self.curve_nameformat is not None:
                args += ["format", self.curve_nameformat]
        return args

    def _output_file_specs(self) -> dict:
        # jktebop writes the model file with the suffix ".jktebop" and the
        # uniformly-sampled model curve (jd or phase) with ".jktebopcurve"
        # (see GetOutputFilename / output_jktebop_curve in USERLIBS/src/jktebop.c).
        return {
            "model": (".jktebop", None),
            "curve": (".jktebopcurve", None),
        }


# -----------------------------------------------------------------------------
# macula — Kipping 2012 spot model.
# -----------------------------------------------------------------------------

class macula(_UserLibCommand):
    """Kipping (2012) Macula rotation + spot light-curve model (USERLIB ``-macula``).

    A *spot* is a dict (or a positional 8-tuple) of value-specs for the
    eight per-spot parameters: ``Lambda0``, ``Phi0``, ``alphamax``,
    ``fspot``, ``tmax``, ``life``, ``ingress``, ``egress``.  Each may be a
    number (``fix``), a source string, or a ``(value, vary_bool)`` tuple to
    mark it as a free fit parameter.

    Parameters
    ----------
    mode : str
        ``"inject"`` or ``"fit amoeba"`` or ``"fit lm"``.
    Prot, istar, kappa2, kappa4, c1, c2, c3, c4, d1, d2, d3, d4, blend : value-spec
    vary_* : bool
    spots : list of dict
        One per spot.  Each dict's keys are the eight per-spot parameter
        names; values are value-specs.  To mark a parameter as free, pass a
        ``(value, True)`` tuple (or use the ``vary_<param>`` dict key with
        a boolean value).
    fluxinput, fluxoutput : bool
    correctlc : bool
    save_model : bool | str | Output
    model_tdelv : bool
    model_nameformat : str, optional
    save_curve : bool | str | Output
    curve_tdelv : bool
    curve_step : float, optional
    curve_nameformat : str, optional
    lib_path : str, optional

    See Also
    --------
    USERLIB extension command: ``-macula``.
    Citation: Kipping 2012 (MNRAS 427, 2487).
    """

    _vt_name = "macula"

    _MAIN_PARAMS = [
        "Prot", "istar", "kappa2", "kappa4",
        "c1", "c2", "c3", "c4",
        "d1", "d2", "d3", "d4",
        "blend",
    ]
    _SPOT_PARAMS = [
        "Lambda0", "Phi0", "alphamax", "fspot",
        "tmax", "life", "ingress", "egress",
    ]

    def __init__(
        self,
        mode: str,
        *,
        spots: Optional[list] = None,
        fluxinput: bool = False,
        fluxoutput: bool = False,
        correctlc: bool = False,
        save_model=False,
        model_tdelv: bool = False,
        model_nameformat: Optional[str] = None,
        save_curve=False,
        curve_tdelv: bool = False,
        curve_step: Optional[float] = None,
        curve_nameformat: Optional[str] = None,
        lib_path: Optional[str] = None,
        **kwargs,
    ) -> None:
        self.mode = mode
        self.spots = spots or []
        self.fluxinput = fluxinput
        self.fluxoutput = fluxoutput
        self.correctlc = correctlc
        self.save_model = save_model
        self.model_tdelv = model_tdelv
        self.model_nameformat = model_nameformat
        self.save_curve = save_curve
        self.curve_tdelv = curve_tdelv
        self.curve_step = curve_step
        self.curve_nameformat = curve_nameformat
        self.lib_path = lib_path
        for p in self._MAIN_PARAMS:
            setattr(self, p, kwargs.pop(p, None))
            setattr(self, f"vary_{p}", bool(kwargs.pop(f"vary_{p}", False)))
        if kwargs:
            raise TypeError(
                f"macula: unexpected keyword arguments: {sorted(kwargs)}"
            )

    @staticmethod
    def _spot_param(spec):
        """Return (value, vary).  Accepts value, (value,), (value, vary)."""
        if isinstance(spec, tuple):
            if len(spec) == 1:
                return spec[0], False
            if len(spec) == 2:
                return spec[0], bool(spec[1])
        return spec, False

    def _to_cli_args(self) -> List[str]:
        outdir = getattr(self, "_outdir", ".")
        args: List[str] = self._libprefix() + ["-macula"] + self.mode.split()
        for p in self._MAIN_PARAMS:
            val = getattr(self, p)
            vary = getattr(self, f"vary_{p}")
            args += [p] + _extparam(val, vary=vary)

        args += ["Nspot", str(len(self.spots))]
        for i, spot in enumerate(self.spots):
            if not isinstance(spot, dict):
                raise TypeError(
                    f"macula: spots[{i}] must be a dict of per-spot parameters"
                )
            for p in self._SPOT_PARAMS:
                if p not in spot:
                    raise ValueError(
                        f"macula: spots[{i}] missing required parameter {p!r}"
                    )
                vary_key = f"vary_{p}"
                val, vary = self._spot_param(spot[p])
                if vary_key in spot:
                    vary = bool(spot[vary_key])
                args += [p] + _extparam(val, vary=vary)

        if self.fluxinput:
            args += ["fluxinput"]
        if self.fluxoutput:
            args += ["fluxoutput"]
        if self.correctlc:
            args += ["correctlc"]
        model_spec = _norm_save(self.save_model)
        if _should_emit(model_spec):
            args += ["omodel", model_spec.path or outdir]
            if self.model_nameformat is not None:
                args += ["nameformat", self.model_nameformat]
            if self.model_tdelv:
                args += ["tdelv"]
        curve_spec = _norm_save(self.save_curve)
        if _should_emit(curve_spec):
            args += ["ocurve", curve_spec.path or outdir]
            if self.curve_nameformat is not None:
                args += ["nameformat", self.curve_nameformat]
            if self.curve_tdelv:
                args += ["tdelv"]
            if self.curve_step is not None:
                args += ["step", str(self.curve_step)]
        return args

    def _output_file_specs(self) -> dict:
        return {
            "model": (".macula", None),
            "curve": (".maculacurve", None),
        }
