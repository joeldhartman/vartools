"""
Pipeline: chains vartools commands and executes them via the binary.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import pandas as pd

from ._binary import get_binary
from ._command import VartoolsCommand
from .commands._helpers import _norm_save, _should_emit
from .lightcurve import LightCurve
from .perlc import PerLC
from .results import (BatchResult, Result, RunError, parse_oneline_output,
                      split_vars_and_scalars)


def _library_enabled() -> bool:
    """Return True if library (in-process) mode is active.

    Checks VARTOOLS_USE_LIBRARY env var and whether libvartoolspipeline is
    loadable.  Currently returns False until Step 3 (Python binding) is done.
    """
    use_lib = os.environ.get("VARTOOLS_USE_LIBRARY", "")
    if use_lib == "0":
        return False
    try:
        from pyvartools._libpipeline import LibPipeline  # noqa: F401
    except ImportError:
        return False
    return True

# Type accepted anywhere a LightCurve is expected
LightCurveInput = Union[LightCurve, pd.DataFrame]

# Type accepted as a file-path argument
FilePath = Union[str, Path]

# Unique per-star variable name injected when running with -parallel to
# track the original list-file row number.  After vartools finishes we sort
# by this column to restore input order, then strip it from the result.
# The name is deliberately obscure to avoid clashing with user variables.
_SEQ_VAR = "_vtpy_seq_"


def _reorder_stats_by_seq(
    stats: pd.DataFrame,
    lc_names: List[str],
) -> pd.DataFrame:
    """Sort *stats* by the injected sequence variable and fix the Name column.

    When ``_SEQ_VAR`` is present (injected for ``-parallel`` runs), rows are
    sorted ascending by it so output order matches the input list, the Name
    column is set from *lc_names* using the stored indices (correctly handling
    any rows dropped by ``-skipmissing``), and the seq column is stripped.

    When no seq column is found (single-threaded path), Name is assigned
    positionally as before.
    """
    seq_cols = [c for c in stats.columns if _SEQ_VAR in c]
    if seq_cols:
        seq_col = seq_cols[0]
        stats = stats.sort_values(seq_col).reset_index(drop=True)
        if "Name" in stats.columns:
            seq_indices = stats[seq_col].astype(int).tolist()
            stats["Name"] = [
                lc_names[i] if i < len(lc_names) else f"lc_{i}"
                for i in seq_indices
            ]
        stats = stats.drop(columns=seq_cols)
    elif not stats.empty and "Name" in stats.columns:
        stats["Name"] = lc_names[:len(stats)]
    return stats


def _apply_columnformat_names(lc: LightCurve, command) -> None:
    """Rename auto-generated ``col4``/``col5``/… columns on *lc* using the
    names declared in a ``cmd.o(..., columnformat=...)`` spec.

    When pyvartools captures a ``cmd.o``-written ASCII output file, any
    extra column beyond ``t``/``mag``/``err`` lands in the DataFrame under
    a placeholder name (``col4``, ``col5``, …).  The user already told us
    what those columns are called by passing ``columnformat``; this helper
    propagates those names so the captured DataFrame matches the
    user-declared column layout.

    No-op when the command has no columnformat or its declared names are
    already present in the DataFrame.
    """
    names = command._columnformat_names() if hasattr(
        command, "_columnformat_names") else None
    if not names:
        return
    have = list(lc._df.columns)
    # Build a rename map: for each user-declared name that isn't already a
    # column, but where a positional col-N placeholder exists, map col-N → name.
    rename: Dict[str, str] = {}
    for pos, desired in enumerate(names):
        if desired in have:
            continue
        placeholder = f"col{pos + 1}"
        if placeholder in have:
            rename[placeholder] = desired
    if rename:
        lc._df.rename(columns=rename, inplace=True)


def _to_lc(obj: LightCurveInput) -> LightCurve:
    """Coerce a DataFrame or LightCurve to a LightCurve."""
    if isinstance(obj, LightCurve):
        return obj
    if isinstance(obj, pd.DataFrame):
        return LightCurve.from_dataframe(obj)
    # Try astropy TimeSeries
    try:
        from astropy.timeseries import TimeSeries
        if isinstance(obj, TimeSeries):
            return LightCurve.from_timeseries(obj)
    except ImportError:
        pass
    raise TypeError(
        f"Expected LightCurve, DataFrame, or astropy TimeSeries; got {type(obj)}"
    )


# ---------------------------------------------------------------------------
# -inputlcformat helpers
# ---------------------------------------------------------------------------

def _inputlcformat_from_df(columns) -> Optional[str]:
    """Build an ``-inputlcformat`` string from a DataFrame's column list.

    Returns ``None`` only when the columns are exactly ``["t", "mag", "err"]``
    in that order — the vartools default, which needs no explicit flag.
    Any other layout (extra columns, missing standard columns, or a different
    order) generates an explicit format string so vartools knows the mapping.
    """
    col_list = list(columns)
    if col_list == ["t", "mag", "err"]:
        return None
    parts = [f"{name}:{i + 1}" for i, name in enumerate(col_list)]
    return ",".join(parts)


def _inputlcformat_from_spec(
    columns,
) -> str:
    """Build an ``-inputlcformat`` string from a user-supplied spec.

    Parameters
    ----------
    columns : list of str  **or**  dict
        * **list** — variable names in column order starting from 1, e.g.
          ``["t", "mag", "err", "airmass"]`` → ``"t:1,mag:2,err:3,airmass:4"``
        * **dict** — explicit mapping of vartools variable name to a
          column spec, where the value can be:

          * an ``int`` (1-based column number for ASCII files);
          * a ``str`` (FITS binary-table column name);
          * an :class:`LCColumn` instance for non-default type/format
            (e.g. ``LCColumn(col=4, type="string")``).
    """
    def _emit(name, val):
        if isinstance(val, LCColumn):
            tail = f"{name}:{val.col}:{val.type}"
            if val.format is not None:
                tail += f":{val.format}"
            return tail
        return f"{name}:{val}"

    if isinstance(columns, dict):
        parts = [_emit(name, col) for name, col in columns.items()]
    else:
        parts = [f"{name}:{i + 1}" for i, name in enumerate(columns)]
    return ",".join(parts)


# ---------------------------------------------------------------------------
# Per-observation and per-star variable descriptors
# ---------------------------------------------------------------------------

@dataclass
class LCVar:
    """Describe a new per-observation variable for ``-inputlcformat`` col=0.

    Pass instances in the ``init_lc_vars`` dict of any Pipeline run method
    to tell vartools to create a variable that is *not* read from the light
    curve file but instead initialised from an analytic expression.

    Parameters
    ----------
    type : str
        Variable type: ``"double"`` (default), ``"float"``, ``"int"``,
        ``"long"``, ``"short"``, ``"string"``, ``"char"``, or ``"utc"``.
    init : str
        Analytic expression evaluated once per observation to initialise
        the variable.  The special variable ``NR`` is the 0-based
        observation index within the light curve.  Defaults to ``"0"``.

    Examples
    --------
    Create a per-observation mask variable initialised to zero::

        pipe.run(lc, init_lc_vars={"mymask": vt.LCVar(type="double", init="0")})

    Create a phase variable based on record number::

        pipe.run(lc, init_lc_vars={"phase": vt.LCVar(init="NR/1000.0")})
    """

    type: str = "double"
    init: str = "0"


@dataclass
class LCColumn:
    """Describe a column read from a light-curve file by ``-inputlcformat``.

    Pass instances as values in the ``columns=`` dict of any Pipeline run
    method when a light-curve column has a non-default type or needs an
    explicit format string (e.g. UTC timestamps, string IDs).  Bare ``int``
    and ``str`` values still work as before for the common case of
    same-type ASCII or FITS columns.

    Parameters
    ----------
    col : int or str
        1-based column number for ASCII files, or column name for FITS
        binary tables.
    type : str
        Variable type — one of ``"double"`` (default), ``"float"``,
        ``"int"``, ``"long"``, ``"short"``, ``"string"``, ``"char"``,
        or ``"utc"``.  ``"utc"`` requires ``format``.
    format : str, optional
        Format string for ``"utc"`` columns (e.g. ``"%Y-%M-%DT%h:%m:%s"``).

    Examples
    --------
    Read a string-typed flag column::

        pipe.run_filelist(
            "lc_list.txt",
            columns={"t": 1, "mag": 2, "err": 3,
                     "fiphot_flag": vt.LCColumn(col=4, type="string")},
        )

    Read a UTC timestamp column with a format spec::

        pipe.run(lc, columns={"t": vt.LCColumn(col=1, type="utc",
                                               format="%Y-%M-%DT%h:%m:%s"),
                              "mag": 2, "err": 3})
    """

    col: Union[int, str]
    type: str = "double"
    format: Optional[str] = None


@dataclass
class ListVar:
    """Describe a per-star variable for ``-inlistvars``.

    Pass instances in the ``inlistvars`` dict of ``run_filelist()`` or
    ``run_batch()`` to tell vartools to read a per-star value from a column
    in the input list file, or to create and initialise the variable from
    an expression when ``col=0``.

    Parameters
    ----------
    col : int
        1-based column number in the input list file.  Use ``0`` to create
        the variable without reading it from the file (requires ``init``).
    type : str
        Variable type: ``"double"`` (default), ``"float"``, ``"int"``,
        ``"long"``, ``"short"``, ``"string"``, ``"char"``, or ``"utc"``.
    init : str, optional
        Analytic expression for initialisation when ``col=0``.  The special
        variable ``NF`` is the 0-based line number in the list file.
    combinelc : bool
        If ``True``, emit the ``combinelc`` keyword (for use with the
        ``-l combinelcs`` feature).  Default ``False``.

    Examples
    --------
    Read per-star period bounds from columns 2 and 3 of a list file::

        batch = pipe.run_filelist(
            "lc_list.txt",
            inlistvars={"minp": vt.ListVar(col=2), "maxp": vt.ListVar(col=3)},
        )

    Equivalently, using the shorthand int form::

        batch = pipe.run_filelist(
            "lc_list.txt",
            inlistvars={"minp": 2, "maxp": 3},
        )

    Initialise from an expression (no list column)::

        batch = pipe.run_filelist(
            "lc_list.txt",
            inlistvars={"minp": vt.ListVar(col=0, type="double", init="0.1")},
        )
    """

    col: int
    type: str = "double"
    init: Optional[str] = None
    combinelc: bool = False


def _inputlcformat_with_init(
    base_fmt: Optional[str],
    init_vars: Dict[str, LCVar],
) -> Optional[str]:
    """Append col=0 init-variable specs to an inputlcformat string.

    Returns ``None`` only when both *base_fmt* is ``None`` and *init_vars*
    is empty (i.e. no ``-inputlcformat`` flag is needed at all).

    When *init_vars* is non-empty but *base_fmt* is ``None`` (the standard
    three-column default), the function emits the explicit default mapping
    ``t:1,mag:2,err:3`` so that adding the col=0 variables does not
    accidentally suppress vartools' knowledge of which columns hold t/mag/err.
    """
    if not init_vars:
        return base_fmt
    # Supplying any -inputlcformat overrides the implicit default, so we must
    # write out the standard mapping explicitly when none was already given.
    base = base_fmt if base_fmt is not None else "t:1,mag:2,err:3"
    parts = [base]
    for varname, lcvar in init_vars.items():
        parts.append(f"{varname}:0:{lcvar.type}:{lcvar.init}")
    return ",".join(parts)


def _inlistvars_from_spec(
    inlistvars: Dict[str, Union[int, ListVar]]
) -> str:
    """Build the argument string for ``-inlistvars``.

    Parameters
    ----------
    inlistvars : dict mapping str to int or ListVar
        * ``int`` — shorthand for ``ListVar(col=N)`` with default type.
        * ``ListVar`` — full specification.
    """
    parts = []
    for varname, spec in inlistvars.items():
        if isinstance(spec, int):
            parts.append(f"{varname}:{spec}")
        else:
            tokens = [varname, str(spec.col)]
            if spec.combinelc:
                tokens.append("combinelc")
            tokens.append(spec.type)
            if spec.init is not None:
                tokens.append(spec.init)
            parts.append(":".join(tokens))
    return ",".join(parts)


class Pipeline:
    """A sequence of vartools commands to execute on one or more light curves.

    Parameters
    ----------
    commands : list of VartoolsCommand
        The commands to chain, in order.

    Examples
    --------
    Single light curve — builder form (preferred)::

        import pyvartools as vt
        lc = vt.LightCurve.from_file("EXAMPLES/2")
        pipe = vt.Pipeline().clip(5.0).LS(0.1, 10.0, 0.1, npeaks=5)
        result = pipe.run(lc)
        print(result.vars)

    Or equivalently, supply a pre-built list of commands::

        from pyvartools import commands as cmd
        pipe = vt.Pipeline([cmd.clip(5.0), cmd.LS(0.1, 10.0, 0.1, npeaks=5)])

    Batch::

        results = pipe.run_batch([lc1, lc2, lc3])
        print(results.vars)  # one row per LC
    """

    def __init__(
        self,
        commands: Optional[Sequence[VartoolsCommand]] = None,
    ) -> None:
        self.commands: List[VartoolsCommand] = list(commands or [])
        self._lib_pipeline = None  # lazily created LibPipeline when library mode is active

    # ------------------------------------------------------------------
    # Single-LC run
    # ------------------------------------------------------------------

    def run(
        self,
        lc: LightCurveInput,
        capture_lc: bool = False,
        outdir: Optional[str] = None,
        timeout: Optional[int] = None,
        init_lc_vars: Optional[Dict[str, LCVar]] = None,
        randseed: Optional[int] = None,
        skipmissing: bool = False,
        jdtol: Optional[float] = None,
        matchstringid: bool = False,
        _command_offset: int = 0,
    ) -> Result:
        """Run the pipeline on a single light curve.

        Parameters
        ----------
        lc : LightCurve | DataFrame | astropy TimeSeries
        capture_lc : bool
            If True, write the (possibly modified) output light curve to a
            temp file and return it as ``result.lc``.
        outdir : str, optional
            Directory for command output files.  Defaults to a fresh temp dir
            that is cleaned up after the run (output DataFrames are captured
            into ``result.files`` instead).
        timeout : int, optional
            Maximum seconds to wait for vartools.  Raises ``RunError`` if
            exceeded.  ``None`` (default) means no limit.
        init_lc_vars : dict mapping str to LCVar, optional
            Per-observation variables to create via ``-inputlcformat`` col=0.
            Each entry adds a ``varname:0:type:init`` token telling vartools
            to define the variable and initialise it from an analytic
            expression.  The special variable ``NR`` is the 0-based
            observation index.  Example::

                init_lc_vars={"mymask": vt.LCVar(type="double", init="0")}
        randseed : int, optional
            Pass ``-randseed N`` to vartools for reproducible random-number
            sequences.
        skipmissing : bool
            Pass ``-skipmissing`` to vartools to silently skip missing input
            files.  Most useful in batch/filelist runs; harmless but redundant
            for single-LC mode.
        jdtol : float, optional
            Pass ``-jdtol N`` to set the tolerance for Julian-date matching.
        matchstringid : bool
            Pass ``-matchstringid`` to force string-based LC name matching.

        Returns
        -------
        Result
        """
        perlc_attrs = self._collect_perlc_attrs()
        if perlc_attrs:
            params = [f"'{name}' in command {ci}" for (ci, name) in perlc_attrs]
            raise ValueError(
                f"Per-LC parameter values cannot be used with run() (single LC). "
                f"Affected parameters: {', '.join(params)}. Use run_batch([lc]) instead."
            )

        lc = _to_lc(lc)

        _has_global_opts = (randseed is not None or skipmissing
                            or jdtol is not None or matchstringid)

        # Fast path: in-process library mode when no output files are needed.
        # Falls back to subprocess when a timeout is requested (library mode
        # has no timeout support), or when init_lc_vars are supplied (library
        # mode does not pass -inputlcformat to vartools_init_pipeline), or
        # when the pipeline contains UserCommand instances (dynamically loaded
        # extensions are not supported by the in-process library).
        if (_library_enabled() and timeout is None and not init_lc_vars
                and not self._has_output_reqs()
                and not self._has_user_commands() and not _has_global_opts):
            if capture_lc:
                return self._run_library_capture(lc, command_offset=_command_offset)
            return self._run_library(lc, command_offset=_command_offset)

        with tempfile.TemporaryDirectory() as tmpdir:
            # Serialize LC to a string and pass via stdin (-i -).
            # vartools names the LC "stdin" internally, so output files
            # (periodograms etc.) will use "stdin" as their basename.
            lc_csv = lc._df.to_csv(sep=" ", header=False, index=False,
                                   float_format="%.10f")

            work_outdir = outdir or tmpdir
            out_lc_path = os.path.join(tmpdir, "output.lc") if capture_lc else None
            fmt = _inputlcformat_with_init(
                _inputlcformat_from_df(lc._df.columns),
                init_lc_vars or {},
            )

            self._assign_o_capture_paths(tmpdir, is_batch=False)

            cmd = self._build_cmd(
                input_flag=["-i", "-"],
                outdir=work_outdir,
                out_lc_path=out_lc_path,
                input_lc_format=fmt,
                randseed=randseed,
                skipmissing=skipmissing,
                jdtol=jdtol,
                matchstringid=matchstringid,
                scalars=lc.scalars,
                command_offset=_command_offset,
            )
            stdout, _ = self._execute(cmd, timeout=timeout, stdin_text=lc_csv)
            stats_full = parse_oneline_output(stdout)
            stats, scalars_df = split_vars_and_scalars(stats_full)

            # Replace the "stdin" name vartools writes as Name with the
            # LightCurve's own name.
            if not stats.empty and "Name" in stats.columns:
                stats["Name"] = lc.name

            # For a single LC, expose var as a Series so that
            # result.vars["LS_Period_1_0"] returns a scalar directly.
            if not stats.empty:
                stats = stats.iloc[0]

            out_lc = None
            if capture_lc and out_lc_path and os.path.isfile(out_lc_path):
                merged_scalars = dict(lc.scalars)
                if not scalars_df.empty:
                    merged_scalars.update(scalars_df.iloc[0].to_dict())
                out_lc = LightCurve.from_file(out_lc_path, name=lc.name)
                out_lc.scalars = merged_scalars

            # When reading from stdin, vartools uses "stdin" as the LC name
            # for output file naming (e.g. periodograms become stdin.ls).
            files = self._collect_output_files("stdin", work_outdir, tmpdir)
            files.update(self._collect_o_captures_single(lc.name))

        return Result(var=stats, lc=out_lc, files=files,
                      known_commands=[c._vt_name for c in self.commands])

    def run_file(
        self,
        lc_path: FilePath,
        capture_lc: bool = False,
        outdir: Optional[str] = None,
        timeout: Optional[int] = None,
        columns: Optional[Union[List[str], Dict[str, Union[int, str]]]] = None,
        init_lc_vars: Optional[Dict[str, LCVar]] = None,
        randseed: Optional[int] = None,
        skipmissing: bool = False,
        jdtol: Optional[float] = None,
        matchstringid: bool = False,
    ) -> Result:
        """Run the pipeline on a light curve file already on disk.

        vartools reads the file directly — no Python I/O is performed.

        Parameters
        ----------
        lc_path : str | Path
            Path to the light curve file.
        capture_lc : bool
            If True, capture the (possibly modified) output LC as ``result.lc``.
        outdir : str, optional
            Directory for command output files.  Defaults to a fresh temp dir.
        timeout : int, optional
            Maximum seconds to wait for vartools.
        columns : list of str  **or**  dict, optional
            Column specification passed to vartools as ``-inputlcformat``.
            * **list** — variable names in column order, e.g.
              ``["t", "mag", "err", "airmass"]``
            * **dict** — explicit mapping of variable name to column number or
              FITS column name, e.g. ``{"t": "BJD_TDB", "mag": "MAG", "err": "ERR"}``
            If ``None`` (default), vartools uses its built-in default mapping
            (columns 1, 2, 3 → t, mag, err).
        init_lc_vars : dict mapping str to LCVar, optional
            Per-observation variables to create via ``-inputlcformat`` col=0.
            Appended to the format string produced from *columns*.  See
            ``run()`` for details.

        Returns
        -------
        Result
        """
        lc_path = str(lc_path)
        lc_name = Path(lc_path).stem

        with tempfile.TemporaryDirectory() as tmpdir:
            work_outdir = outdir or tmpdir
            out_lc_path = os.path.join(tmpdir, "output.lc") if capture_lc else None
            base_fmt = _inputlcformat_from_spec(columns) if columns is not None else None
            fmt = _inputlcformat_with_init(base_fmt, init_lc_vars or {})

            self._assign_o_capture_paths(tmpdir, is_batch=False)

            cmd = self._build_cmd(
                input_flag=["-i", lc_path],
                outdir=work_outdir,
                out_lc_path=out_lc_path,
                input_lc_format=fmt,
                randseed=randseed,
                skipmissing=skipmissing,
                jdtol=jdtol,
                matchstringid=matchstringid,
            )
            stdout, _ = self._execute(cmd, timeout=timeout)
            stats_full = parse_oneline_output(stdout)
            stats, scalars_df = split_vars_and_scalars(stats_full)

            if not stats.empty and "Name" in stats.columns:
                stats["Name"] = lc_name
            if not stats.empty:
                stats = stats.iloc[0]

            out_lc = None
            if capture_lc and out_lc_path and os.path.isfile(out_lc_path):
                out_lc = LightCurve.from_file(out_lc_path, name=lc_name)
                if not scalars_df.empty:
                    out_lc.scalars = dict(scalars_df.iloc[0].to_dict())

            files = self._collect_output_files(lc_path, work_outdir, tmpdir)
            files.update(self._collect_o_captures_single(lc_name))

        return Result(var=stats, lc=out_lc, files=files,
                      known_commands=[c._vt_name for c in self.commands])

    def run_filelist(
        self,
        lc_paths: Union[FilePath, Sequence[FilePath]],
        nthreads: int = 1,
        capture_lc: bool = False,
        outdir: Optional[str] = None,
        timeout: Optional[int] = None,
        raise_on_error: bool = True,
        columns: Optional[Union[List[str], Dict[str, Union[int, str]]]] = None,
        init_lc_vars: Optional[Dict[str, LCVar]] = None,
        inlistvars: Optional[Dict[str, Union[int, ListVar]]] = None,
        combinelcs: bool = False,
        lcnumvar: Optional[str] = "lcnum",
        randseed: Optional[int] = None,
        skipmissing: bool = False,
        jdtol: Optional[float] = None,
        matchstringid: bool = False,
    ) -> BatchResult:
        """Run the pipeline on a collection of light curve files on disk.

        vartools reads the files directly — no Python I/O is performed.

        Parameters
        ----------
        lc_paths : str | Path | list of str | Path
            Either a path to an existing vartools list file (one LC path per
            line), or a sequence of individual LC file paths (a temporary list
            file is written automatically).
        nthreads : int
            Number of parallel threads (``-parallel``).
        capture_lc : bool
            If True, capture the (possibly modified) output LC for each light
            curve and return them as ``result.lcs``.
        outdir : str, optional
            Directory for command output files.  Defaults to a fresh temp dir.
        timeout : int, optional
            Maximum seconds to wait for vartools.
        raise_on_error : bool
            If False, a vartools failure is stored in ``result.error`` rather
            than raised.
        columns : list of str  **or**  dict, optional
            Column specification passed to vartools as ``-inputlcformat``.
            * **list** — variable names in column order, e.g.
              ``["t", "mag", "err", "airmass"]``
            * **dict** — explicit mapping of variable name to column number or
              FITS column name, e.g. ``{"t": "BJD_TDB", "mag": "MAG", "err": "ERR"}``
            If ``None`` (default), vartools uses its built-in default mapping.
        init_lc_vars : dict mapping str to LCVar, optional
            Per-observation variables to create via ``-inputlcformat`` col=0.
            Appended to the format string produced from *columns*.  See
            ``run()`` for details.
        inlistvars : dict mapping str to int or ListVar, optional
            Per-star variables passed to vartools via ``-inlistvars``.
            Each entry defines a variable read from a column of the input list
            file, or initialised from an expression when col=0.

            Use an ``int`` as shorthand for ``ListVar(col=N)``::

                inlistvars={"minp": 2, "maxp": 3}

            Use a ``ListVar`` for full control over type and initialisation::

                inlistvars={"minp": vt.ListVar(col=2, type="double")}

            Per-star variables defined here can then be referenced by name
            in LS (and other commands) via the ``var`` form, e.g.
            ``cmd.LS("minp", "maxp", 1e-3)``.
        combinelcs : bool
            If True, append ``combinelcs`` to the ``-l`` flag so vartools
            treats each line of the list file as a *group* of comma-separated
            paths combined into one in-memory light curve.  The list file (or
            list of strings passed in *lc_paths*) is responsible for the
            grouping; this flag does not split anything itself.
        lcnumvar : str, optional
            Only used when ``combinelcs=True``.  Name of the per-observation
            integer variable vartools creates to record which file each point
            came from.  Defaults to ``"lcnum"``; pass ``None`` to opt out.

        Returns
        -------
        BatchResult
        """
        perlc_attrs = self._collect_perlc_attrs()
        if combinelcs and perlc_attrs:
            params = [f"'{name}' in command {ci}" for (ci, name) in perlc_attrs]
            raise ValueError(
                f"PerLC parameter values cannot be used with combinelcs=True. "
                f"Use run_combinelcs() (which does not yet support PerLC), or "
                f"remove the PerLC values.  Affected parameters: "
                f"{', '.join(params)}."
            )
        if perlc_attrs and isinstance(lc_paths, (str, Path)):
            raise ValueError(
                "Per-LC parameter values (PerLC / numpy arrays) cannot be used with "
                "a pre-built list file. Either pass a list of file paths to "
                "run_filelist(), or use run_batch() with LightCurve objects."
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            # Resolve lc_paths to a list file and a list of individual paths.
            if isinstance(lc_paths, (str, Path)):
                list_path = str(lc_paths)
                with open(list_path) as f:
                    paths = [line.strip() for line in f if line.strip()]
            else:
                paths = [str(p) for p in lc_paths]

            if combinelcs:
                lc_names = [Path(p.split(",", 1)[0]).stem for p in paths]
            else:
                lc_names = [Path(p).stem for p in paths]
            work_outdir = outdir or tmpdir
            out_lc_dir = os.path.join(tmpdir, "lc_out") if capture_lc else None
            if out_lc_dir:
                os.makedirs(out_lc_dir, exist_ok=True)
            nth_args = ["-parallel", str(nthreads)] if nthreads > 1 else []
            base_fmt = _inputlcformat_from_spec(columns) if columns is not None else None
            fmt = _inputlcformat_with_init(base_fmt, init_lc_vars or {})

            col_assignments = {}
            perlc_subs = {}
            if perlc_attrs:
                batch_size = len(paths)
                for (ci, name), perlc in perlc_attrs.items():
                    if len(perlc) != batch_size:
                        raise ValueError(
                            f"PerLC parameter '{name}' in command {ci} has {len(perlc)} values "
                            f"but the batch has {batch_size} light curves."
                        )
                next_col = 2
                for key in sorted(perlc_attrs):
                    col_assignments[key] = next_col
                    next_col += 1
                perlc_subs = self._build_perlc_subs(col_assignments)

            list_path = os.path.join(tmpdir, "lclist.txt")
            if perlc_attrs:
                self._write_perlc_list_file(list_path, paths, perlc_attrs, col_assignments)
            else:
                with open(list_path, "w") as f:
                    for p in paths:
                        f.write(p + "\n")

            # Merge user-supplied inlistvars with auto-generated per-LC vars.
            # When running in parallel, also inject the sequence-index variable
            # so we can restore input order after vartools finishes.
            merged_inlistvars = dict(inlistvars) if inlistvars else {}
            if col_assignments:
                merged_inlistvars.update(self._build_perlc_inlistvars(col_assignments))
            use_seq = nthreads > 1
            if use_seq:
                merged_inlistvars[_SEQ_VAR] = ListVar(col=0, type="int", init="NF")
            inlistvars_str = _inlistvars_from_spec(merged_inlistvars) if merged_inlistvars else None

            self._assign_o_capture_paths(tmpdir, is_batch=True)

            input_flag = ["-l", list_path]
            if combinelcs:
                input_flag.append("combinelcs")
                # If the user opted out of lcnumvar but still requested
                # capture_lc=True, force the default back on so the captured
                # LC carries an lcnum column identifying source files.
                effective_lcnumvar = lcnumvar
                if capture_lc and not effective_lcnumvar:
                    effective_lcnumvar = "lcnum"
                if effective_lcnumvar:
                    input_flag += ["lcnumvar", effective_lcnumvar]
            cmd = self._build_cmd(
                input_flag=input_flag,
                outdir=work_outdir,
                out_lc_dir=out_lc_dir,
                nth_args=nth_args,
                input_lc_format=fmt,
                inlistvars_str=inlistvars_str,
                perlc_subs=perlc_subs,
                randseed=randseed,
                skipmissing=skipmissing,
                jdtol=jdtol,
                matchstringid=matchstringid,
                inject_print_var=_SEQ_VAR if use_seq else None,
            )
            try:
                stdout, _ = self._execute(cmd, timeout=timeout)
            except RunError as exc:
                if raise_on_error:
                    raise
                return BatchResult(var=pd.DataFrame(), error=exc)

            stats = parse_oneline_output(stdout)
            stats = _reorder_stats_by_seq(stats, lc_names)
            stats, scalars_df = split_vars_and_scalars(stats)

            # When combinelcs=True, vartools names per-LC output files after the
            # *first* path of each comma-joined group rather than the full line.
            if combinelcs:
                file_keys = [p.split(",", 1)[0] for p in paths]
            else:
                file_keys = list(paths)

            out_lcs = None
            if capture_lc:
                out_lcs = []
                for i, (lc_path, name) in enumerate(zip(file_keys, lc_names)):
                    opath = os.path.join(out_lc_dir, Path(lc_path).name)
                    if os.path.isfile(opath):
                        new_lc = LightCurve.from_file(opath, name=name)
                        if not scalars_df.empty and i < len(scalars_df):
                            new_lc.scalars = dict(scalars_df.iloc[i].to_dict())
                        out_lcs.append(new_lc)
                    else:
                        out_lcs.append(None)

            all_files: dict = {}
            for lc_path in file_keys:
                lc_files = self._collect_output_files(lc_path, work_outdir, tmpdir)
                for name, df in lc_files.items():
                    all_files.setdefault(name, []).append(df)
            for key, lc_list in self._collect_o_captures_batch(file_keys, lc_names).items():
                all_files[key] = lc_list

        return BatchResult(var=stats, lcs=out_lcs, files=all_files,
                               known_commands=[c._vt_name for c in self.commands])

    def run_combinelc(
        self,
        files: Sequence[FilePath],
        nthreads: int = 1,
        capture_lc: bool = False,
        outdir: Optional[str] = None,
        timeout: Optional[int] = None,
        raise_on_error: bool = True,
        columns: Optional[Union[List[str], Dict[str, Union[int, str]]]] = None,
        init_lc_vars: Optional[Dict[str, LCVar]] = None,
        inlistvars: Optional[Dict[str, Union[int, ListVar]]] = None,
        lcnumvar: Optional[str] = "lcnum",
        delimiter: str = ",",
        randseed: Optional[int] = None,
        skipmissing: bool = False,
        jdtol: Optional[float] = None,
        matchstringid: bool = False,
    ) -> Result:
        """Combine *files* into a single light curve and run the pipeline.

        Single-group convenience wrapper around :meth:`run_combinelcs`.  All
        keyword arguments forward to ``run_combinelcs``; the result is the
        first (and only) row, returned as a :class:`Result`.

        Parameters
        ----------
        files : sequence of str | Path
            Paths to combine into one in-memory light curve.

        Returns
        -------
        Result

        Examples
        --------
        ::

            result = (vt.Pipeline()
                      .stitch("mag err mask lcnum")
                      .run_combinelc(["seg1.txt", "seg2.txt", "seg3.txt"]))
        """
        files = list(files)
        if not files:
            raise ValueError("run_combinelc() requires at least one file path.")
        batch = self.run_combinelcs(
            groups=[files],
            nthreads=nthreads,
            capture_lc=capture_lc,
            outdir=outdir,
            timeout=timeout,
            raise_on_error=raise_on_error,
            columns=columns,
            init_lc_vars=init_lc_vars,
            inlistvars=inlistvars,
            lcnumvar=lcnumvar,
            delimiter=delimiter,
            randseed=randseed,
            skipmissing=skipmissing,
            jdtol=jdtol,
            matchstringid=matchstringid,
        )
        if batch.error is not None:
            return Result(
                var=pd.Series(dtype=object),
                error=batch.error,
                known_commands=[c._vt_name for c in self.commands],
            )
        return batch[0]

    def run_combinelcs(
        self,
        groups: Sequence[Sequence[FilePath]],
        nthreads: int = 1,
        capture_lc: bool = False,
        outdir: Optional[str] = None,
        timeout: Optional[int] = None,
        raise_on_error: bool = True,
        columns: Optional[Union[List[str], Dict[str, Union[int, str]]]] = None,
        init_lc_vars: Optional[Dict[str, LCVar]] = None,
        inlistvars: Optional[Dict[str, Union[int, ListVar]]] = None,
        lcnumvar: Optional[str] = "lcnum",
        delimiter: str = ",",
        randseed: Optional[int] = None,
        skipmissing: bool = False,
        jdtol: Optional[float] = None,
        matchstringid: bool = False,
    ) -> BatchResult:
        """Run the pipeline using vartools ``-l … combinelcs`` mode.

        Each entry in *groups* is a list of file paths that vartools combines
        into a single in-memory light curve.  The result contains one row in
        ``result.vars`` per group.

        Parameters
        ----------
        groups : sequence of sequences of str | Path
            Each inner sequence is one group of files to combine.  Within a
            group the paths are joined by *delimiter* to form one line in the
            vartools list file.
        nthreads : int
            Number of parallel threads (``-parallel``).
        capture_lc : bool
            If True, capture the (possibly modified) output LC for each group.
        outdir : str, optional
            Directory for command output files.  Defaults to a fresh temp dir.
        timeout : int, optional
            Maximum seconds to wait for vartools.
        raise_on_error : bool
            If False, a vartools failure is stored in ``result.error`` rather
            than raised.
        columns : list of str  **or**  dict, optional
            Column specification passed to vartools as ``-inputlcformat``.
        init_lc_vars : dict mapping str to LCVar, optional
            Per-observation variables to create via ``-inputlcformat`` col=0.
        inlistvars : dict mapping str to int or ListVar, optional
            Per-star variables passed to vartools via ``-inlistvars``.
        lcnumvar : str, optional
            Name of the per-observation integer variable vartools creates to
            record which file each point came from.  Defaults to ``"lcnum"``;
            pass ``None`` to opt out of emitting the ``lcnumvar`` qualifier.
        delimiter : str
            Delimiter used to join paths within each group in the list file.
            Default ``","`` (the vartools ``combinelcs`` default).
        randseed : int, optional
            Pass ``-randseed N`` to vartools.
        skipmissing : bool
            Pass ``-skipmissing`` to vartools.
        jdtol : float, optional
            Pass ``-jdtol N`` to vartools.
        matchstringid : bool
            Pass ``-matchstringid`` to vartools.

        Notes
        -----
        PerLC array parameter values are supported: each PerLC must have
        ``len(groups)`` entries (one value per combined group).  The values
        are appended as additional columns in the temporary list file and
        wired up via ``-inlistvars``.

        Returns
        -------
        BatchResult
        """
        perlc_attrs = self._collect_perlc_attrs()
        if perlc_attrs:
            batch_size = len(groups)
            for (ci, name), perlc in perlc_attrs.items():
                if len(perlc) != batch_size:
                    raise ValueError(
                        f"PerLC parameter '{name}' in command {ci} has "
                        f"{len(perlc)} values but the batch has {batch_size} "
                        f"groups."
                    )

        with tempfile.TemporaryDirectory() as tmpdir:
            # Assign extra columns for any PerLC params (col 2 onward — col 1
            # is the comma-joined paths field for each group).
            col_assignments = {}
            perlc_subs = {}
            if perlc_attrs:
                next_col = 2
                for key in sorted(perlc_attrs):
                    col_assignments[key] = next_col
                    next_col += 1
                perlc_subs = self._build_perlc_subs(col_assignments)

            # Build list file: one line per group, paths joined by delimiter,
            # followed by any PerLC value columns.
            list_path = os.path.join(tmpdir, "lclist.txt")
            group_strings = [delimiter.join(str(p) for p in group)
                             for group in groups]
            if perlc_attrs:
                self._write_perlc_list_file(
                    list_path, group_strings, perlc_attrs, col_assignments
                )
            else:
                with open(list_path, "w") as f:
                    for line in group_strings:
                        f.write(line + "\n")

            # Build the -l input flag with combinelcs.  If the user opted out
            # of lcnumvar but still requested capture_lc=True, force the
            # default name back on so the captured LC carries an lcnum column
            # identifying the source file of each point.
            effective_lcnumvar = lcnumvar
            if capture_lc and not effective_lcnumvar:
                effective_lcnumvar = "lcnum"
            input_flag = ["-l", list_path, "combinelcs"]
            if effective_lcnumvar:
                input_flag += ["lcnumvar", effective_lcnumvar]

            lc_names = [Path(group[0]).stem for group in groups]
            work_outdir = outdir or tmpdir
            out_lc_dir = os.path.join(tmpdir, "lc_out") if capture_lc else None
            if out_lc_dir:
                os.makedirs(out_lc_dir, exist_ok=True)
            nth_args = ["-parallel", str(nthreads)] if nthreads > 1 else []
            base_fmt = _inputlcformat_from_spec(columns) if columns is not None else None
            fmt = _inputlcformat_with_init(base_fmt, init_lc_vars or {})
            use_seq = nthreads > 1
            merged_inlistvars_comb = dict(inlistvars) if inlistvars else {}
            if col_assignments:
                merged_inlistvars_comb.update(
                    self._build_perlc_inlistvars(col_assignments)
                )
            if use_seq:
                merged_inlistvars_comb[_SEQ_VAR] = ListVar(col=0, type="int", init="NF")
            inlistvars_str = _inlistvars_from_spec(merged_inlistvars_comb) if merged_inlistvars_comb else None

            self._assign_o_capture_paths(tmpdir, is_batch=True)

            cmd = self._build_cmd(
                input_flag=input_flag,
                outdir=work_outdir,
                out_lc_dir=out_lc_dir,
                nth_args=nth_args,
                input_lc_format=fmt,
                inlistvars_str=inlistvars_str,
                perlc_subs=perlc_subs,
                randseed=randseed,
                skipmissing=skipmissing,
                jdtol=jdtol,
                matchstringid=matchstringid,
                inject_print_var=_SEQ_VAR if use_seq else None,
            )
            try:
                stdout, _ = self._execute(cmd, timeout=timeout)
            except RunError as exc:
                if raise_on_error:
                    raise
                return BatchResult(var=pd.DataFrame(), error=exc)

            stats = parse_oneline_output(stdout)
            stats = _reorder_stats_by_seq(stats, lc_names)
            stats, scalars_df = split_vars_and_scalars(stats)

            out_lcs = None
            if capture_lc:
                out_lcs = []
                for i, (group, name) in enumerate(zip(groups, lc_names)):
                    # vartools names the combined output after the first file
                    opath = os.path.join(out_lc_dir, Path(group[0]).name)
                    if os.path.isfile(opath):
                        new_lc = LightCurve.from_file(opath, name=name)
                        if not scalars_df.empty and i < len(scalars_df):
                            new_lc.scalars = dict(scalars_df.iloc[i].to_dict())
                        out_lcs.append(new_lc)
                    else:
                        out_lcs.append(None)

            all_files: dict = {}
            # Use first path of each group as the key for output file lookup
            first_paths = [str(group[0]) for group in groups]
            for lc_path in first_paths:
                lc_files = self._collect_output_files(lc_path, work_outdir, tmpdir)
                for name, df in lc_files.items():
                    all_files.setdefault(name, []).append(df)
            for key, lc_list in self._collect_o_captures_batch(first_paths, lc_names).items():
                all_files[key] = lc_list

        return BatchResult(var=stats, lcs=out_lcs, files=all_files,
                               known_commands=[c._vt_name for c in self.commands])

    # ------------------------------------------------------------------
    # Batch run
    # ------------------------------------------------------------------

    def run_batch(
        self,
        lcs: Sequence[LightCurveInput],
        nthreads: int = 1,
        capture_lc: bool = False,
        outdir: Optional[str] = None,
        timeout: Optional[int] = None,
        raise_on_error: bool = True,
        init_lc_vars: Optional[Dict[str, LCVar]] = None,
        inlistvars: Optional[Dict[str, Union[int, ListVar]]] = None,
        randseed: Optional[int] = None,
        skipmissing: bool = False,
        jdtol: Optional[float] = None,
        matchstringid: bool = False,
        _command_offset: int = 0,
    ) -> BatchResult:
        """Run the pipeline on a list of light curves.

        vartools processes them in a single invocation using ``-l``.

        Parameters
        ----------
        lcs : list of LightCurve | DataFrame | astropy TimeSeries
        nthreads : int
            Number of parallel threads (passed to vartools as ``-parallel``).
        capture_lc : bool
            If True, capture the (possibly modified) output LC for each light
            curve and return them as ``result.lcs``.
        outdir : str, optional
            Directory for output files.  Defaults to a temp dir.
        timeout : int, optional
            Maximum seconds to wait for vartools.  Raises ``RunError`` if
            exceeded (or stores it in ``result.error`` when
            ``raise_on_error=False``).
        raise_on_error : bool
            If False, a vartools failure is caught and stored in
            ``result.error`` rather than raised.  ``result.vars`` will be
            empty in that case.
        init_lc_vars : dict mapping str to LCVar, optional
            Per-observation variables to create via ``-inputlcformat`` col=0.
            Appended to the auto-generated format string.  See ``run()`` for
            details.
        inlistvars : dict mapping str to int or ListVar, optional
            Per-star variables passed to vartools via ``-inlistvars``.
            Note: for ``run_batch()`` the list file is a temporary file
            containing only LC paths (no extra columns).  This parameter is
            therefore only useful with ``col=0`` ``ListVar`` entries that
            initialise variables from expressions rather than reading from
            list columns.  To supply per-star values from a file use
            ``run_filelist()`` instead.

        Returns
        -------
        BatchResult
        """
        perlc_attrs = self._collect_perlc_attrs()
        lcs = [_to_lc(lc) for lc in lcs]

        _has_global_opts = (randseed is not None or skipmissing
                            or jdtol is not None or matchstringid)

        # Collect per-LC carried-forward scalars.  These are injected via the
        # -inlistvars mechanism (as INLIST variables) so each LC sees its own
        # value in downstream expressions — using -expr const would apply a
        # single value across the whole batch, which is wrong.
        batch_scalars = self._collect_batch_scalars(lcs)

        # Fast path: in-process library mode when no output files are needed
        # and parallel processing is not requested (library mode is single-threaded).
        # Also skip when UserCommand instances are present — dynamically loaded
        # extension libraries are not supported by the in-process library.
        # The library-mode fast path also does not support per-LC scalar
        # injection (no list-file machinery), so we route any batch with
        # carried-forward scalars or a non-zero chain offset to the subprocess
        # path unconditionally.
        if (_library_enabled() and nthreads == 1 and not init_lc_vars
                and not capture_lc and not self._has_output_reqs()
                and not perlc_attrs and not self._has_user_commands()
                and not _has_global_opts
                and not batch_scalars and _command_offset == 0):
            return self._run_batch_library(lcs, raise_on_error=raise_on_error)

        with tempfile.TemporaryDirectory() as tmpdir:
            # Write each LC to a numbered temp file; record names
            list_path = os.path.join(tmpdir, "lclist.txt")
            lc_paths = []
            for i, lc in enumerate(lcs):
                p = os.path.join(tmpdir, f"lc_{i:06d}.lc")
                lc._df.to_csv(p, sep=" ", header=False, index=False,
                              float_format="%.10f")
                lc_paths.append(p)

            col_assignments = {}
            perlc_subs = {}
            scalar_col_assignments: Dict[str, int] = {}
            next_col = 2
            if perlc_attrs:
                batch_size = len(lcs)
                for (ci, name), perlc in perlc_attrs.items():
                    if len(perlc) != batch_size:
                        raise ValueError(
                            f"PerLC parameter '{name}' in command {ci} has {len(perlc)} values "
                            f"but the batch has {batch_size} light curves."
                        )
                for key in sorted(perlc_attrs):
                    col_assignments[key] = next_col
                    next_col += 1
                perlc_subs = self._build_perlc_subs(col_assignments)

            # Assign list-file columns for carried-forward scalars, continuing
            # after any PerLC column allocations.
            for name in batch_scalars:
                scalar_col_assignments[name] = next_col
                next_col += 1

            # Build a single dict {col: per-LC values} that unifies PerLC and
            # scalar columns, then write the list file.
            col_to_values: Dict[int, list] = {}
            for (ci, name), col in col_assignments.items():
                col_to_values[col] = list(perlc_attrs[(ci, name)])
            for name, col in scalar_col_assignments.items():
                col_to_values[col] = batch_scalars[name]
            self._write_extra_cols_list_file(list_path, lc_paths, col_to_values)

            work_outdir = outdir or tmpdir
            out_lc_dir = os.path.join(tmpdir, "lc_out") if capture_lc else None
            if out_lc_dir:
                os.makedirs(out_lc_dir, exist_ok=True)
            nth_args = ["-parallel", str(nthreads)] if nthreads > 1 else []
            # Auto-discover extra columns from the first LC (all LCs in a batch
            # are expected to share the same column structure).
            base_fmt = _inputlcformat_from_df(lcs[0]._df.columns) if lcs else None
            fmt = _inputlcformat_with_init(base_fmt, init_lc_vars or {})
            # Merge user-supplied inlistvars with auto-generated per-LC vars
            # and carried-forward scalars.  Scalars are registered by their
            # actual variable names (e.g. "LS_Period_1_0") so downstream
            # expressions can reference them directly.
            merged_inlistvars = dict(inlistvars) if inlistvars else {}
            if col_assignments:
                merged_inlistvars.update(self._build_perlc_inlistvars(col_assignments))
            if scalar_col_assignments:
                merged_inlistvars.update(scalar_col_assignments)
            use_seq = nthreads > 1
            if use_seq:
                merged_inlistvars[_SEQ_VAR] = ListVar(col=0, type="int", init="NF")
            inlistvars_str = _inlistvars_from_spec(merged_inlistvars) if merged_inlistvars else None

            self._assign_o_capture_paths(tmpdir, is_batch=True)

            cmd = self._build_cmd(
                input_flag=["-l", list_path],
                outdir=work_outdir,
                out_lc_dir=out_lc_dir,
                nth_args=nth_args,
                input_lc_format=fmt,
                inlistvars_str=inlistvars_str,
                perlc_subs=perlc_subs,
                randseed=randseed,
                skipmissing=skipmissing,
                jdtol=jdtol,
                matchstringid=matchstringid,
                inject_print_var=_SEQ_VAR if use_seq else None,
                command_offset=_command_offset,
                harvest_scalars=bool(scalar_col_assignments),
            )
            try:
                stdout, _ = self._execute(cmd, timeout=timeout)
            except RunError as exc:
                if raise_on_error:
                    raise
                return BatchResult(var=pd.DataFrame(), error=exc)

            stats = parse_oneline_output(stdout)
            # Restore input order (may be scrambled by -parallel) and replace
            # temp-file paths in the Name column with the original LC names.
            stats = _reorder_stats_by_seq(stats, [lc.name for lc in lcs])
            stats, scalars_df = split_vars_and_scalars(stats)

            # Drop echoed INLIST values for scalars we injected — users
            # already have those on the input LCs, and echoing them as
            # "new" scalars would be noise.  Preserve genuinely new ones.
            if scalar_col_assignments and not scalars_df.empty:
                injected_names = set(scalar_col_assignments.keys())
                keep_cols = [c for c in scalars_df.columns if c not in injected_names]
                scalars_df = scalars_df[keep_cols] if keep_cols else pd.DataFrame(index=scalars_df.index)

            out_lcs = None
            if capture_lc:
                out_lcs = []
                for i, (lc, lc_path) in enumerate(zip(lcs, lc_paths)):
                    opath = os.path.join(out_lc_dir, Path(lc_path).name)
                    if os.path.isfile(opath):
                        new_lc = LightCurve.from_file(opath, name=lc.name)
                        merged = dict(lc.scalars)
                        if not scalars_df.empty and i < len(scalars_df):
                            merged.update(scalars_df.iloc[i].to_dict())
                        new_lc.scalars = merged
                        out_lcs.append(new_lc)
                    else:
                        out_lcs.append(None)

            # Collect per-LC output files if any commands requested them
            all_files: dict = {}
            for i, lc in enumerate(lcs):
                lc_files = self._collect_output_files(lc_paths[i], work_outdir, tmpdir)
                for name, df in lc_files.items():
                    all_files.setdefault(name, []).append(df)
            lc_names = [lc.name for lc in lcs]
            for key, lc_list in self._collect_o_captures_batch(lc_paths, lc_names).items():
                all_files[key] = lc_list

        return BatchResult(var=stats, lcs=out_lcs, files=all_files,
                               known_commands=[c._vt_name for c in self.commands])

    # ------------------------------------------------------------------
    # Library mode helpers
    # ------------------------------------------------------------------

    def _has_output_reqs(self) -> bool:
        """True if any command needs the subprocess path for file I/O.

        Any ``save_*`` directive that wants the file captured into
        ``result.files`` forces subprocess mode.  So does any ``cmd.o(...)``
        instance — the library-mode pipeline does not know how to construct
        per-LC output filenames from array-backed LCs, so both
        ``capture=True`` and plain ``filename=...`` forms of ``-o`` must run
        through the subprocess path.
        """
        from .commands.misc import o as OCommand
        for command in self.commands:
            if command._requested_outputs():
                return True
            if isinstance(command, OCommand):
                return True
        return False

    def _has_user_commands(self) -> bool:
        """True if any command is a UserCommand (forces subprocess mode)."""
        from pyvartools.userlib import UserCommand
        from pyvartools.commands.userlibs import _UserLibCommand
        return any(
            isinstance(c, (UserCommand, _UserLibCommand))
            for c in self.commands
        )

    def _scalar_injection_args(
        self,
        scalars: Optional[Dict[str, float]] = None,
    ) -> List[str]:
        """Return argv tokens that pre-register each scalar as -expr const.

        Each ``(name, value)`` pair becomes the pair of tokens
        ``["-expr", "const", "name=value"]`` (three tokens total per scalar
        since ``-expr`` with the ``const`` qualifier takes the qualifier as a
        separate argv word).  Emitted at the head of the command list so
        subsequent commands can resolve the names during parsecommandline.
        """
        if not scalars:
            return []
        tokens: List[str] = []
        for name, val in scalars.items():
            tokens += ["-expr", "const", f"{name}={val!r}"]
        return tokens

    def _commands_to_argv(
        self,
        scalars: Optional[Dict[str, float]] = None,
        command_offset: int = 0,
    ) -> List[str]:
        """Build a CLI arg list from pipeline commands (for LibPipeline init).

        Parameters
        ----------
        scalars : dict, optional
            Per-star scalars to pre-register (see :meth:`_scalar_injection_args`).
        command_offset : int
            When > 0, emit ``-columnsuffix <N>`` before each command so its
            output columns end in ``_<command_offset + i>`` instead of ``_<i>``.
            Used by chained calls to keep suffixes growing across segments and
            to avoid collision with injected scalar names of the form
            ``CMDNAME_descriptor_<i>``.
        """
        args: List[str] = []
        args += self._scalar_injection_args(scalars)
        for i, command in enumerate(self.commands):
            if command_offset > 0:
                args += ["-columnsuffix", str(command_offset + i)]
            args += command._to_cli_args()
        # -printallscalars is harmless when no scalars exist and enables round-
        # tripping of user-created scalars (from -expr scalar / listvar) into
        # result.lc.scalars.  Only emit when chained (command_offset > 0) or
        # when scalars were injected, to avoid changing output for existing
        # tests.
        if command_offset > 0 or scalars:
            args += ["-printallscalars"]
        args += ["-oneline"]
        return args

    def _lib_argv_with_format(
        self,
        lc: LightCurve,
        scalars: Optional[Dict[str, float]] = None,
        command_offset: int = 0,
    ) -> list:
        """Build LibPipeline argv including -inputlcformat for extra columns."""
        argv = self._commands_to_argv(scalars=scalars, command_offset=command_offset)
        fmt = _inputlcformat_from_df(lc._df.columns)
        if fmt is not None:
            argv = ["-inputlcformat", fmt] + argv
        return argv

    def _ensure_lib_pipeline(self, lc: LightCurve, command_offset: int = 0):
        """Create or recreate the LibPipeline if needed.

        The key covers both the LC column layout and the (scalars,
        command_offset) used to build the argv, so a change in injected-
        scalar names or in the chain position forces a rebuild.
        """
        from pyvartools._libpipeline import LibPipeline
        fmt = _inputlcformat_from_df(lc._df.columns)
        scalars_key = tuple(sorted(lc.scalars.keys())) if lc.scalars else ()
        key = (fmt, scalars_key, command_offset)
        current_key = getattr(self, "_lib_pipeline_key", None)
        if self._lib_pipeline is None or key != current_key:
            self._lib_pipeline = LibPipeline(
                self._lib_argv_with_format(
                    lc, scalars=lc.scalars, command_offset=command_offset))
            self._lib_pipeline_key = key
            self._lib_pipeline_fmt = fmt  # kept for backward compat

    def _run_library(self, lc: LightCurve, command_offset: int = 0) -> Result:
        """Execute one LC via LibPipeline (init-once, reused on subsequent calls)."""
        try:
            self._ensure_lib_pipeline(lc, command_offset=command_offset)
            extra_cols = {}
            for col in lc._df.columns:
                if col not in ("t", "mag", "err"):
                    extra_cols[col] = lc._df[col].values
            if extra_cols:
                stats, _, scalars = self._lib_pipeline.process_lc_capture(
                    lc.t, lc.mag, lc.err, name=lc.name,
                    extra_columns=extra_cols)
            else:
                stats = self._lib_pipeline.process_lc(lc.t, lc.mag, lc.err, name=lc.name)
                scalars = {}
        except RuntimeError as exc:
            self._lib_pipeline = None  # allow retry after failure
            raise RunError(str(exc)) from exc
        if isinstance(stats, pd.Series) and "Name" in stats.index:
            # Rebuild as object-dtype Series so the Name (string) is preserved
            # correctly — the raw output may have coerced it to a numeric value.
            data = stats.to_dict()
            data["Name"] = lc.name
            stats = pd.Series(data)
        # When capture_lc=False we still want prior scalars (carried on the
        # input lc) to be recoverable — but since no output LC is built here,
        # the scalars are dropped.  This matches capture_lc=False semantics
        # for LC data.  Users who need scalar round-tripping should capture.
        _ = scalars
        return Result(var=stats, lc=None, files={},
                      known_commands=[c._vt_name for c in self.commands])

    def _run_library_capture(self, lc: LightCurve, command_offset: int = 0) -> Result:
        """Execute one LC via LibPipeline and capture the modified LC."""
        try:
            self._ensure_lib_pipeline(lc, command_offset=command_offset)

            extra_cols = {}
            for col in lc._df.columns:
                if col not in ("t", "mag", "err"):
                    extra_cols[col] = lc._df[col].values

            stats, lc_columns, scalars = self._lib_pipeline.process_lc_capture(
                lc.t, lc.mag, lc.err, name=lc.name,
                extra_columns=extra_cols if extra_cols else None,
            )
        except RuntimeError as exc:
            self._lib_pipeline = None
            raise RunError(str(exc)) from exc

        if isinstance(stats, pd.Series) and "Name" in stats.index:
            data = stats.to_dict()
            data["Name"] = lc.name
            stats = pd.Series(data)

        # Merge input lc.scalars with harvested scalars so prior chain state
        # flows through even when the new run doesn't redefine those names.
        merged_scalars = dict(lc.scalars)
        merged_scalars.update(scalars or {})

        out_lc = None
        if lc_columns and "t" in lc_columns and "mag" in lc_columns:
            lc_df = pd.DataFrame(lc_columns)
            col_order = ["t", "mag"]
            if "err" in lc_columns:
                col_order.append("err")
            for c in sorted(lc_columns.keys()):
                if c not in col_order:
                    col_order.append(c)
            lc_df = lc_df[col_order]
            out_lc = LightCurve(lc_df, name=lc.name, scalars=merged_scalars)

        return Result(var=stats, lc=out_lc, files={},
                      known_commands=[c._vt_name for c in self.commands])

    def _run_batch_library(
        self, lcs: List[LightCurve], raise_on_error: bool = True
    ) -> BatchResult:
        """Execute a list of LCs via LibPipeline (init-once, loop per LC)."""
        from pyvartools._libpipeline import LibPipeline
        try:
            if self._lib_pipeline is None:
                self._lib_pipeline = LibPipeline(self._commands_to_argv())
        except RuntimeError as exc:
            self._lib_pipeline = None
            err = RunError(str(exc))
            if raise_on_error:
                raise err from exc
            return BatchResult(var=pd.DataFrame(), error=err)
        rows = []
        for lc in lcs:
            stats = self._lib_pipeline.process_lc(lc.t, lc.mag, lc.err, name=lc.name)
            if isinstance(stats, pd.Series) and "Name" in stats.index:
                data = stats.to_dict()
                data["Name"] = lc.name
                stats = pd.Series(data)
            rows.append(stats)
        df = pd.DataFrame(rows).reset_index(drop=True)
        return BatchResult(var=df, lcs=None, files={},
                           known_commands=[c._vt_name for c in self.commands])

    # ------------------------------------------------------------------
    # Per-LC parameter helpers
    # ------------------------------------------------------------------

    def _collect_perlc_attrs(self) -> dict:
        """Return {(cmd_idx, attr_name): PerLC} for all per-LC params in commands."""
        import numpy as np
        result = {}
        for i, command in enumerate(self.commands):
            for name, val in vars(command).items():
                if name.startswith('_'):
                    continue
                if isinstance(val, PerLC):
                    result[(i, name)] = val
                elif (isinstance(val, np.ndarray)
                      and val.ndim == 1
                      and np.issubdtype(val.dtype, np.number)):
                    result[(i, name)] = PerLC(val)
                else:
                    try:
                        import pandas as pd
                        if isinstance(val, pd.Series):
                            result[(i, name)] = PerLC(val)
                    except ImportError:
                        pass
        return result

    def _write_perlc_list_file(self, list_path, lc_paths, perlc_attrs, col_assignments):
        """Write the list file with per-LC value columns appended."""
        ordered_keys = sorted(col_assignments, key=lambda k: col_assignments[k])
        with open(list_path, "w") as f:
            for j, p in enumerate(lc_paths):
                parts = [p]
                for key in ordered_keys:
                    parts.append(f"{perlc_attrs[key][j]:.10g}")
                f.write(" ".join(parts) + "\n")

    def _write_extra_cols_list_file(self, list_path, lc_paths, col_to_perlc_values):
        """Write a list file with arbitrary extra columns appended per line.

        Parameters
        ----------
        col_to_perlc_values : dict
            Mapping of 1-based column number → list of per-LC values (length
            must equal ``len(lc_paths)``).  Columns are written in ascending
            column-number order.  Column numbers must start at 2 (column 1
            is the LC path itself) and be contiguous; missing numbers are
            not supported.
        """
        if not col_to_perlc_values:
            with open(list_path, "w") as f:
                for p in lc_paths:
                    f.write(p + "\n")
            return
        sorted_cols = sorted(col_to_perlc_values.keys())
        with open(list_path, "w") as f:
            for j, p in enumerate(lc_paths):
                parts = [p]
                for col in sorted_cols:
                    parts.append(f"{col_to_perlc_values[col][j]:.10g}")
                f.write(" ".join(parts) + "\n")

    def _collect_batch_scalars(self, lcs) -> dict:
        """Collect per-LC scalar values, keyed by scalar name.

        Every LightCurve in *lcs* must carry the same set of scalar names
        (they all originate from the same prior chain segment, so this is
        the expected invariant).  Returns ``{name: [per-LC values]}``.
        Raises ``ValueError`` if the name sets differ.
        """
        if not lcs:
            return {}
        name_sets = [set(lc.scalars.keys()) for lc in lcs]
        reference = name_sets[0]
        for i, ns in enumerate(name_sets[1:], 1):
            if ns != reference:
                only_ref = reference - ns
                only_other = ns - reference
                raise ValueError(
                    f"LC index {i} has a different set of carried-forward "
                    f"scalar names than LC 0.  LC 0 extras: {sorted(only_ref)}; "
                    f"LC {i} extras: {sorted(only_other)}.  All LCs in a batch "
                    f"must share the same .scalars key set."
                )
        return {name: [lc.scalars[name] for lc in lcs] for name in sorted(reference)}

    def _build_perlc_subs(self, col_assignments) -> dict:
        """Return {cmd_idx: {attr_name: "expr varname"}} substitution map.

        Each per-LC attribute is replaced by ``"expr <varname>"`` where
        ``<varname>`` is a unique name defined via ``-inlistvars``.  We use
        ``expr`` rather than ``var`` because a handful of vartools commands
        accept ``expr`` on their value-spec parameters but not ``var`` —
        e.g. ``-BLSFixPer`` on its period spec.  Passing a bare identifier
        through ``expr`` is semantically equivalent to ``var`` for every
        command that accepts both, so the change is safe across the board.
        """
        subs = {}
        for (ci, name), col in col_assignments.items():
            varname = f"_perlc_{ci}_{name}"
            subs.setdefault(ci, {})[name] = f"expr {varname}"
        return subs

    def _build_perlc_inlistvars(self, col_assignments) -> dict:
        """Return {varname: col} inlistvars dict for per-LC attributes."""
        result = {}
        for (ci, name), col in col_assignments.items():
            varname = f"_perlc_{ci}_{name}"
            result[varname] = col
        return result

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_cmd(
        self,
        input_flag: List[str],
        outdir: str,
        out_lc_path: Optional[str] = None,
        out_lc_dir: Optional[str] = None,
        nth_args: Optional[List[str]] = None,
        input_lc_format: Optional[str] = None,
        inlistvars_str: Optional[str] = None,
        perlc_subs: Optional[dict] = None,
        randseed: Optional[int] = None,
        skipmissing: bool = False,
        jdtol: Optional[float] = None,
        matchstringid: bool = False,
        inject_print_var: Optional[str] = None,
        scalars: Optional[Dict[str, float]] = None,
        command_offset: int = 0,
        harvest_scalars: bool = False,
    ) -> List[str]:
        """Assemble the full vartools command line."""
        binary = get_binary()
        cmd = [binary] + input_flag
        if input_lc_format:
            cmd += ["-inputlcformat", input_lc_format]
        if inlistvars_str:
            cmd += ["-inlistvars", inlistvars_str]
        if nth_args:
            cmd += nth_args
        if randseed is not None:
            cmd += ["-randseed", str(randseed)]
        if jdtol is not None:
            cmd += ["-jdtol", str(jdtol)]
        if matchstringid:
            cmd += ["-matchstringid"]
        if skipmissing:
            cmd += ["-skipmissing"]
        # Pre-register carried-forward scalars as -expr const before any user
        # command so subsequent commands can reference them by name.
        cmd += self._scalar_injection_args(scalars)
        for idx, command in enumerate(self.commands):
            # Give each command that writes output files its own subdirectory
            # so that two commands of the same type don't overwrite each other.
            # Commands that specify an explicit output path (Mode 2/3) get that
            # directory created instead; others fall back to a temp subdir.
            specs = command._output_file_specs()
            if specs:
                base_cmd_outdir = os.path.join(outdir, f"cmd_{idx}")
                outdir_map = {}
                needs_base = False
                for name, spec_tuple in specs.items():
                    # _output_file_specs entries are ``(suffix, parser)`` for
                    # per-LC directory-style outputs (the default) or
                    # ``(suffix, parser, "file")`` for commands like -SYSREM
                    # that write a single global file (otrends).  In the file
                    # case the user-supplied path *is* the output file and
                    # must not be makedirs-ed.
                    spec_mode = spec_tuple[2] if len(spec_tuple) >= 3 else "dir"
                    save_spec = _norm_save(getattr(command, f"save_{name}", False))
                    if save_spec.path is not None:
                        if spec_mode == "file":
                            parent = os.path.dirname(save_spec.path)
                            if parent:
                                os.makedirs(parent, exist_ok=True)
                        else:
                            os.makedirs(save_spec.path, exist_ok=True)
                        outdir_map[name] = save_spec.path
                    else:
                        needs_base = True
                        outdir_map[name] = base_cmd_outdir
                if needs_base:
                    os.makedirs(base_cmd_outdir, exist_ok=True)
                command._outdir = base_cmd_outdir
                command._outdir_map = outdir_map
            else:
                command._outdir = outdir
                command._outdir_map = {}
            # Emit an explicit -columnsuffix for each user command so its
            # output-column names end in "_<command_offset+idx>" rather than
            # "_<idx>".  Required in chain continuations to avoid collision
            # with injected scalar names like "LS_Period_1_0" that pyvartools
            # carried forward from a previous segment.
            if command_offset > 0:
                cmd += ["-columnsuffix", str(command_offset + idx)]
            subs = perlc_subs.get(idx, {}) if perlc_subs else {}
            if subs:
                cmd += command._to_cli_args_with_perlc(subs)
            else:
                cmd += command._to_cli_args()
        if inject_print_var:
            cmd += ["-print", inject_print_var]
        # -printallscalars rounds per-star scalar state (from -expr scalar /
        # listvar and injected -expr const / -inlistvars values) into the
        # stdout stream so pyvartools can re-capture it for the next chain
        # segment.
        if command_offset > 0 or scalars or harvest_scalars:
            cmd += ["-printallscalars"]
        cmd += ["-oneline"]
        if out_lc_path:
            # Single-LC mode: explicit output path.  "allcols" writes every
            # currently-registered LC-vector variable, mirroring the
            # library-mode capture path which also includes new vectors
            # (e.g. -Phase phasevar, -linfit modelvar).
            cmd += ["-o", out_lc_path, "allcols"]
        elif out_lc_dir:
            # Batch mode: vartools writes <out_lc_dir>/<input_basename> for each LC.
            cmd += ["-o", out_lc_dir, "allcols"]
        return cmd

    def _execute(
        self,
        cmd: List[str],
        timeout: Optional[int] = None,
        stdin_text: Optional[str] = None,
    ) -> tuple:
        """Run cmd, capture stdout/stderr, raise RunError on failure."""
        try:
            proc = subprocess.run(
                cmd,
                input=stdin_text,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except FileNotFoundError as e:
            raise RunError(f"Failed to launch vartools: {e}") from e
        except subprocess.TimeoutExpired as e:
            raise RunError(
                f"vartools timed out after {timeout}s.\n"
                f"Command: {' '.join(cmd)}"
            ) from e

        if proc.returncode != 0:
            raise RunError(
                f"vartools exited with status {proc.returncode}.\n"
                f"Command: {' '.join(cmd)}\n"
                f"stderr: {proc.stderr.strip()}"
            )
        return proc.stdout, proc.stderr

    def _collect_output_files(
        self, lc_name: str, outdir: str, tmpdir: str
    ) -> dict:
        """Read any output files requested via save_*=True on commands.

        Keys in the returned dict have the form
        ``"{command._vt_name}_{logical_name}_{idx}"`` where *idx* is the
        zero-based position of the command in the pipeline.  This ensures
        that two commands of the same type (e.g. two ``LS`` runs) or two
        different commands that share a logical name (e.g. ``LS`` and
        ``BLS`` both using ``"periodogram"``) never clobber each other.
        """
        files = {}
        base = Path(lc_name).name  # strip directory from lc name
        for idx, command in enumerate(self.commands):
            specs = command._output_file_specs()
            if not specs:
                continue
            mandatory = getattr(command, "_mandatory_output", False)
            for logical_name, spec_tuple in specs.items():
                # Optional 3rd element is the spec mode (e.g. ``"file"``
                # for SYSREM trends); not used for capture-side parsing.
                suffix, ncols = spec_tuple[0], spec_tuple[1]
                raw = getattr(command, f"save_{logical_name}", False)
                spec = _norm_save(raw)

                # Skip if not emitting at all (Mode 4), unless the command
                # mandates output regardless of save spec (e.g. autocorrelation).
                if not _should_emit(spec) and not mandatory:
                    continue
                # Skip if not capturing into Python.
                if not spec.capture:
                    continue

                # Locate the file using the per-output directory recorded in
                # _build_cmd, falling back to the old-style cmd_outdir.
                outdir_map = getattr(command, "_outdir_map", {})
                if logical_name in outdir_map:
                    actual_outdir = outdir_map[logical_name]
                else:
                    cmd_outdir = os.path.join(outdir, f"cmd_{idx}")
                    actual_outdir = cmd_outdir if os.path.isdir(cmd_outdir) else outdir

                candidate = os.path.join(actual_outdir, base + suffix)
                if os.path.isfile(candidate):
                    kwargs = {}
                    if ncols:
                        kwargs["usecols"] = list(range(ncols))
                    df = pd.read_csv(candidate, sep=r"\s+", comment="#",
                                     header=None, **kwargs)
                    key = f"{command._vt_name}_{logical_name}_{idx}"
                    files[key] = df
        return files

    # ------------------------------------------------------------------
    # cmd.o capture helpers
    # ------------------------------------------------------------------

    def _assign_o_capture_paths(self, tmpdir: str, is_batch: bool) -> None:
        """Set ``_capture_path`` on any ``cmd.o(capture=True)`` commands.

        Called before ``_build_cmd()`` so the paths are ready when
        ``_to_cli_args()`` is invoked.  Commands with an explicit
        ``filename`` already know their path and are left untouched.
        """
        from .commands.misc import o as OCommand
        for idx, command in enumerate(self.commands):
            if not (isinstance(command, OCommand) and command.capture
                    and command.filename is None):
                continue
            # Use a directory when in batch mode or when nameformat is set
            # (nameformat causes vartools to construct its own filename inside
            # the given directory, so we can't predict the exact path).
            if is_batch or command.nameformat is not None:
                cap_dir = os.path.join(tmpdir, f"_o_cap_{idx}")
                os.makedirs(cap_dir, exist_ok=True)
                command._capture_path = cap_dir
            else:
                command._capture_path = os.path.join(tmpdir, f"_o_cap_{idx}.lc")

    def _collect_o_captures_single(self, lc_name: str) -> dict:
        """Return captured ``cmd.o`` files for a single-LC run.

        Returns a dict mapping ``command.key`` → ``LightCurve``.
        """
        from .commands.misc import o as OCommand
        files = {}
        for command in self.commands:
            if not (isinstance(command, OCommand) and command.capture):
                continue
            path = (command.filename
                    if command.filename is not None
                    else command._capture_path)
            if path is None:
                continue
            lc = None
            if os.path.isdir(path):
                # nameformat case: single file somewhere inside the directory
                candidates = [
                    os.path.join(path, f)
                    for f in os.listdir(path)
                    if os.path.isfile(os.path.join(path, f))
                ]
                if len(candidates) == 1:
                    lc = LightCurve.from_file(candidates[0], name=lc_name)
            elif os.path.isfile(path):
                lc = LightCurve.from_file(path, name=lc_name)
            if lc is not None:
                _apply_columnformat_names(lc, command)
                files[command.key] = lc
        return files

    def _collect_o_captures_batch(
        self, lc_paths: List[str], lc_names: List[str]
    ) -> dict:
        """Return captured ``cmd.o`` files for a batch run.

        Returns a dict mapping ``command.key`` → list of ``LightCurve``
        (one per LC, ``None`` where the output file is missing).
        """
        from .commands.misc import o as OCommand
        files = {}
        for command in self.commands:
            if not (isinstance(command, OCommand) and command.capture):
                continue
            base_dir = (command.filename
                        if command.filename is not None
                        else command._capture_path)
            if base_dir is None:
                continue
            lc_list = []
            for lc_path, name in zip(lc_paths, lc_names):
                out_path = os.path.join(base_dir, Path(lc_path).name)
                if os.path.isfile(out_path):
                    one_lc = LightCurve.from_file(out_path, name=name)
                    _apply_columnformat_names(one_lc, command)
                    lc_list.append(one_lc)
                else:
                    lc_list.append(None)
            files[command.key] = lc_list
        return files

    # ------------------------------------------------------------------
    # Convenience: add commands and chain pipelines
    # ------------------------------------------------------------------

    def add(self, command: VartoolsCommand) -> "Pipeline":
        """Append a command and return self (for fluent chaining)."""
        self.commands.append(command)
        return self

    def __repr__(self) -> str:
        return f"Pipeline([{', '.join(repr(c) for c in self.commands)}])"
