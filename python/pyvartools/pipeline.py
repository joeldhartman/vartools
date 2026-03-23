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
from .results import BatchResult, Result, RunError, parse_oneline_output


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
    columns: Union[List[str], Dict[str, Union[int, str]]]
) -> str:
    """Build an ``-inputlcformat`` string from a user-supplied spec.

    Parameters
    ----------
    columns : list of str  **or**  dict mapping str → int | str
        * **list** — variable names in column order starting from 1, e.g.
          ``["t", "mag", "err", "airmass"]`` → ``"t:1,mag:2,err:3,airmass:4"``
        * **dict** — explicit mapping of vartools variable name to 1-based
          column number *or* FITS column name, e.g.
          ``{"t": "BJD_TDB", "mag": "MAG", "err": "ERR"}`` (FITS)
          or ``{"t": 1, "mag": 2, "err": 3, "airmass": 4}`` (ASCII).
    """
    if isinstance(columns, dict):
        parts = [f"{name}:{col}" for name, col in columns.items()]
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
    Single light curve::

        import pyvartools as vt
        lc = vt.LightCurve.from_file("EXAMPLES/2")
        pipe = vt.Pipeline([vt.commands.LS(0.1, 10.0, 0.1, 5, 1)])
        result = pipe.run(lc)
        print(result.stats)

    Batch::

        results = pipe.run_batch([lc1, lc2, lc3])
        print(results.stats)  # one row per LC
    """

    def __init__(self, commands: Sequence[VartoolsCommand]) -> None:
        self.commands: List[VartoolsCommand] = list(commands)
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
                and not capture_lc and not self._has_output_reqs()
                and not self._has_user_commands() and not _has_global_opts):
            return self._run_library(lc)

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
            )
            stdout, _ = self._execute(cmd, timeout=timeout, stdin_text=lc_csv)
            stats = parse_oneline_output(stdout)

            # Replace the "stdin" name vartools writes as Name with the
            # LightCurve's own name.
            if not stats.empty and "Name" in stats.columns:
                stats["Name"] = lc.name

            # For a single LC, expose stats as a Series so that
            # result.stats["LS_Period_1_0"] returns a scalar directly.
            if not stats.empty:
                stats = stats.iloc[0]

            out_lc = None
            if capture_lc and out_lc_path and os.path.isfile(out_lc_path):
                out_lc = LightCurve.from_file(out_lc_path, name=lc.name)

            # When reading from stdin, vartools uses "stdin" as the LC name
            # for output file naming (e.g. periodograms become stdin.ls).
            files = self._collect_output_files("stdin", work_outdir, tmpdir)
            files.update(self._collect_o_captures_single(lc.name))

        return Result(stats=stats, lc=out_lc, files=files)

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
            stats = parse_oneline_output(stdout)

            if not stats.empty and "Name" in stats.columns:
                stats["Name"] = lc_name
            if not stats.empty:
                stats = stats.iloc[0]

            out_lc = None
            if capture_lc and out_lc_path and os.path.isfile(out_lc_path):
                out_lc = LightCurve.from_file(out_lc_path, name=lc_name)

            files = self._collect_output_files(lc_path, work_outdir, tmpdir)
            files.update(self._collect_o_captures_single(lc_name))

        return Result(stats=stats, lc=out_lc, files=files)

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

        Returns
        -------
        BatchResult
        """
        perlc_attrs = self._collect_perlc_attrs()
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

            # Merge user-supplied inlistvars with auto-generated per-LC vars
            merged_inlistvars = dict(inlistvars) if inlistvars else {}
            if col_assignments:
                merged_inlistvars.update(self._build_perlc_inlistvars(col_assignments))
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
            )
            try:
                stdout, _ = self._execute(cmd, timeout=timeout)
            except RunError as exc:
                if raise_on_error:
                    raise
                return BatchResult(stats=pd.DataFrame(), error=exc)

            stats = parse_oneline_output(stdout)

            if not stats.empty and "Name" in stats.columns:
                stats["Name"] = lc_names[:len(stats)]

            out_lcs = None
            if capture_lc:
                out_lcs = []
                for lc_path, name in zip(paths, lc_names):
                    opath = os.path.join(out_lc_dir, Path(lc_path).name)
                    if os.path.isfile(opath):
                        out_lcs.append(LightCurve.from_file(opath, name=name))
                    else:
                        out_lcs.append(None)

            all_files: dict = {}
            for lc_path in paths:
                lc_files = self._collect_output_files(lc_path, work_outdir, tmpdir)
                for name, df in lc_files.items():
                    all_files.setdefault(name, []).append(df)
            for key, lc_list in self._collect_o_captures_batch(paths, lc_names).items():
                all_files[key] = lc_list

        return BatchResult(stats=stats, lcs=out_lcs, files=all_files)

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
        lcnumvar: Optional[str] = None,
        delimiter: str = ",",
        randseed: Optional[int] = None,
        skipmissing: bool = False,
        jdtol: Optional[float] = None,
        matchstringid: bool = False,
    ) -> BatchResult:
        """Run the pipeline using vartools ``-l … combinelcs`` mode.

        Each entry in *groups* is a list of file paths that vartools combines
        into a single in-memory light curve.  The result contains one row in
        ``result.stats`` per group.

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
            If given, pass ``lcnumvar <name>`` after ``combinelcs`` in the
            ``-l`` flag so vartools creates a per-observation integer variable
            recording which file each point came from.
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

        Returns
        -------
        BatchResult
        """
        perlc_attrs = self._collect_perlc_attrs()
        if perlc_attrs:
            params = [f"'{name}' in command {ci}" for (ci, name) in perlc_attrs]
            raise ValueError(
                f"PerLC parameter values cannot be used with run_combinelcs(). "
                f"Affected parameters: {', '.join(params)}."
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            # Build list file: one line per group, paths joined by delimiter
            list_path = os.path.join(tmpdir, "lclist.txt")
            with open(list_path, "w") as f:
                for group in groups:
                    f.write(delimiter.join(str(p) for p in group) + "\n")

            # Build the -l input flag with combinelcs
            input_flag = ["-l", list_path, "combinelcs"]
            if lcnumvar:
                input_flag += ["lcnumvar", lcnumvar]

            lc_names = [Path(group[0]).stem for group in groups]
            work_outdir = outdir or tmpdir
            out_lc_dir = os.path.join(tmpdir, "lc_out") if capture_lc else None
            if out_lc_dir:
                os.makedirs(out_lc_dir, exist_ok=True)
            nth_args = ["-parallel", str(nthreads)] if nthreads > 1 else []
            base_fmt = _inputlcformat_from_spec(columns) if columns is not None else None
            fmt = _inputlcformat_with_init(base_fmt, init_lc_vars or {})
            inlistvars_str = _inlistvars_from_spec(inlistvars) if inlistvars else None

            self._assign_o_capture_paths(tmpdir, is_batch=True)

            cmd = self._build_cmd(
                input_flag=input_flag,
                outdir=work_outdir,
                out_lc_dir=out_lc_dir,
                nth_args=nth_args,
                input_lc_format=fmt,
                inlistvars_str=inlistvars_str,
                randseed=randseed,
                skipmissing=skipmissing,
                jdtol=jdtol,
                matchstringid=matchstringid,
            )
            try:
                stdout, _ = self._execute(cmd, timeout=timeout)
            except RunError as exc:
                if raise_on_error:
                    raise
                return BatchResult(stats=pd.DataFrame(), error=exc)

            stats = parse_oneline_output(stdout)

            if not stats.empty and "Name" in stats.columns:
                stats["Name"] = lc_names[:len(stats)]

            out_lcs = None
            if capture_lc:
                out_lcs = []
                for group, name in zip(groups, lc_names):
                    # vartools names the combined output after the first file
                    opath = os.path.join(out_lc_dir, Path(group[0]).name)
                    if os.path.isfile(opath):
                        out_lcs.append(LightCurve.from_file(opath, name=name))
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

        return BatchResult(stats=stats, lcs=out_lcs, files=all_files)

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
            ``result.error`` rather than raised.  ``result.stats`` will be
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

        # Fast path: in-process library mode when no output files are needed
        # and parallel processing is not requested (library mode is single-threaded).
        # Also skip when UserCommand instances are present — dynamically loaded
        # extension libraries are not supported by the in-process library.
        if (_library_enabled() and nthreads == 1 and not init_lc_vars
                and not capture_lc and not self._has_output_reqs()
                and not perlc_attrs and not self._has_user_commands()
                and not _has_global_opts):
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
            if perlc_attrs:
                batch_size = len(lcs)
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
                self._write_perlc_list_file(list_path, lc_paths, perlc_attrs, col_assignments)
            else:
                with open(list_path, "w") as f:
                    for p in lc_paths:
                        f.write(p + "\n")

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
            merged_inlistvars = dict(inlistvars) if inlistvars else {}
            if col_assignments:
                merged_inlistvars.update(self._build_perlc_inlistvars(col_assignments))
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
            )
            try:
                stdout, _ = self._execute(cmd, timeout=timeout)
            except RunError as exc:
                if raise_on_error:
                    raise
                return BatchResult(stats=pd.DataFrame(), error=exc)

            stats = parse_oneline_output(stdout)

            # Replace temp-file paths in the Name column with the original LC names.
            if not stats.empty and "Name" in stats.columns:
                stats["Name"] = [lc.name for lc in lcs[:len(stats)]]

            out_lcs = None
            if capture_lc:
                out_lcs = []
                for lc, lc_path in zip(lcs, lc_paths):
                    opath = os.path.join(out_lc_dir, Path(lc_path).name)
                    if os.path.isfile(opath):
                        out_lcs.append(LightCurve.from_file(opath, name=lc.name))
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

        return BatchResult(stats=stats, lcs=out_lcs, files=all_files)

    # ------------------------------------------------------------------
    # Library mode helpers
    # ------------------------------------------------------------------

    def _has_output_reqs(self) -> bool:
        """True if any command requests output files (save_* or cmd.o capture)."""
        from .commands.misc import o as OCommand
        for command in self.commands:
            if command._requested_outputs():
                return True
            if isinstance(command, OCommand) and command.capture:
                return True
        return False

    def _has_user_commands(self) -> bool:
        """True if any command is a UserCommand (forces subprocess mode)."""
        from pyvartools.userlib import UserCommand
        return any(isinstance(c, UserCommand) for c in self.commands)

    def _commands_to_argv(self) -> List[str]:
        """Build a CLI arg list from pipeline commands (for LibPipeline init)."""
        args: List[str] = []
        for command in self.commands:
            args += command._to_cli_args()
        args += ["-oneline"]
        return args

    def _run_library(self, lc: LightCurve) -> Result:
        """Execute one LC via LibPipeline (init-once, reused on subsequent calls)."""
        from pyvartools._libpipeline import LibPipeline
        try:
            if self._lib_pipeline is None:
                self._lib_pipeline = LibPipeline(self._commands_to_argv())
            stats = self._lib_pipeline.process_lc(lc.t, lc.mag, lc.err, name=lc.name)
        except RuntimeError as exc:
            self._lib_pipeline = None  # allow retry after failure
            raise RunError(str(exc)) from exc
        if isinstance(stats, pd.Series) and "Name" in stats.index:
            # Rebuild as object-dtype Series so the Name (string) is preserved
            # correctly — the raw output may have coerced it to a numeric value.
            data = stats.to_dict()
            data["Name"] = lc.name
            stats = pd.Series(data)
        return Result(stats=stats, lc=None, files={})

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
            return BatchResult(stats=pd.DataFrame(), error=err)
        rows = []
        for lc in lcs:
            stats = self._lib_pipeline.process_lc(lc.t, lc.mag, lc.err, name=lc.name)
            if isinstance(stats, pd.Series) and "Name" in stats.index:
                data = stats.to_dict()
                data["Name"] = lc.name
                stats = pd.Series(data)
            rows.append(stats)
        df = pd.DataFrame(rows).reset_index(drop=True)
        return BatchResult(stats=df, lcs=None, files={})

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

    def _build_perlc_subs(self, col_assignments) -> dict:
        """Return {cmd_idx: {attr_name: varname}} substitution map.

        Each per-LC attribute is replaced by a unique variable name that will
        be defined via -inlistvars.  ``_varexpr(varname)`` then emits
        ``["var", varname]`` so vartools reads the value from the list column.
        """
        subs = {}
        for (ci, name), col in col_assignments.items():
            varname = f"_perlc_{ci}_{name}"
            subs.setdefault(ci, {})[name] = varname
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
                for name in specs:
                    save_spec = _norm_save(getattr(command, f"save_{name}", False))
                    if save_spec.path is not None:
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
            subs = perlc_subs.get(idx, {}) if perlc_subs else {}
            if subs:
                cmd += command._to_cli_args_with_perlc(subs)
            else:
                cmd += command._to_cli_args()
        cmd += ["-oneline"]
        if out_lc_path:
            # Single-LC mode: explicit output path.
            cmd += ["-o", out_lc_path]
        elif out_lc_dir:
            # Batch mode: vartools writes <out_lc_dir>/<input_basename> for each LC.
            cmd += ["-o", out_lc_dir]
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
            for logical_name, (suffix, ncols) in specs.items():
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
                    lc_list.append(LightCurve.from_file(out_path, name=name))
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
