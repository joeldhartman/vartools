"""
LightCurveBatch: fluent batch processing for a collection of LightCurves.

Runs one vartools invocation per light curve sequentially, collecting per-LC
results into a BatchResult that supports both DataFrame and per-LC access.

Usage::

    batch = LightCurveBatch([lc1, lc2, lc3])
    result = batch.LS(0.1, 10.0).rms().run()

    result.vars                     # DataFrame, one row per LC
    result[0].varobjs.LS.Period_1   # first LC's LS period
    result[1].error                 # RunError if second LC failed, else None

    for r in result:
        print(r.varobjs.LS.Period_1)

For power-user batch runs that need the efficiency of a single vartools
invocation, use ``Pipeline.run_batch()`` / ``Pipeline.run_filelist()``
directly.
"""

from __future__ import annotations

import copy
from typing import Any, Dict, Iterable, List, Optional, Union, TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from ._command import VartoolsCommand
    from .lightcurve import LightCurve
    from .results import BatchResult, Result, RunError


# ---------------------------------------------------------------------------
# Per-LC parameter helpers
# ---------------------------------------------------------------------------

def _is_perlc_value(val) -> bool:
    """Return True if *val* should be treated as a per-LC parameter array.

    Recognised types:

    * ``PerLC`` — always per-LC (explicit wrapper)
    * 1-D ``numpy.ndarray`` — any dtype, including object/string arrays
    * ``pandas.Series`` — any dtype; name-based index alignment is supported
      in :meth:`_extract_perlc_scalar`

    Plain Python :class:`list` values are intentionally **not** auto-detected
    to avoid ambiguity with fixed multi-valued parameters such as
    ``MandelAgolTransit(ld_coeffs=[0.236, 0.391])``.  Wrap them explicitly
    with ``PerLC([...])`` if per-LC behaviour is intended.
    """
    from .perlc import PerLC
    if isinstance(val, PerLC):
        return True
    if isinstance(val, np.ndarray) and val.ndim == 1:
        return True
    if isinstance(val, pd.Series):
        return True
    return False


def _extract_perlc_scalar(val, lc_idx: int, lc_name: str):
    """Return the scalar value for the LC at *lc_idx* from a per-LC array.

    For ``pandas.Series`` with a non-integer index the value is looked up by
    *lc_name* first; if the name is not in the index it falls back to
    positional access at *lc_idx*.  This lets callers pass a Series that is
    indexed by light-curve name (e.g. ``br1.vars.set_index("Name")["col"]``).

    ``PerLC`` values are always converted to Python ``float``.  numpy and
    pandas scalars are converted to the closest Python built-in type via
    ``.item()`` so that downstream CLI-rendering code sees plain ``int`` /
    ``float`` / ``str`` rather than numpy wrapper types.
    """
    from .perlc import PerLC

    if isinstance(val, PerLC):
        return float(val[lc_idx])

    if isinstance(val, pd.Series):
        if not pd.api.types.is_integer_dtype(val.index.dtype):
            # Non-integer (e.g. string) index — try name-based lookup first.
            if lc_name in val.index:
                item = val.loc[lc_name]
            else:
                item = val.iloc[lc_idx]
        else:
            item = val.iloc[lc_idx]
        # Unwrap numpy/pandas scalars to Python built-ins.
        if hasattr(item, "item"):
            return item.item()
        return item

    if isinstance(val, np.ndarray):
        item = val[lc_idx]
        if hasattr(item, "item"):
            return item.item()
        return item

    return val


class LightCurveList(list):
    """A ``list`` of captured ``LightCurve``\\s with positional + by-name access.

    Returned by ``BatchResult.files[key]`` for any ``cmd.o(capture=True)``
    output produced in batch mode.  Subclasses ``list`` so all of pandas
    / numpy / standard slicing works, and adds:

    * ``lcs[i]``        — same as ``list.__getitem__`` (int or slice).
    * ``lcs['name']``   — return the first ``LightCurve`` whose ``.name``
                          matches; ``KeyError`` if none.  ``None``
                          placeholders (missing-file slots) are skipped.
    * ``'name' in lcs`` — membership by name (string) or by identity
                          (``LightCurve`` instance).

    Positional alignment with the input batch is preserved — entries are
    ``None`` where the captured output file was missing.
    """

    def __getitem__(self, key):
        if isinstance(key, str):
            from .lightcurve import LightCurve
            for lc in self:
                if isinstance(lc, LightCurve) and lc.name == key:
                    return lc
            available = [lc.name for lc in self
                         if isinstance(lc, LightCurve)]
            raise KeyError(
                f"No LightCurve in list with name={key!r}.  "
                f"Available: {available}"
            )
        return super().__getitem__(key)

    def __contains__(self, item) -> bool:
        if isinstance(item, str):
            from .lightcurve import LightCurve
            return any(isinstance(lc, LightCurve) and lc.name == item
                       for lc in self)
        return super().__contains__(item)


class LightCurveBatch:
    """A collection of LightCurves with a fluent command-chaining interface.

    Parameters
    ----------
    lcs : iterable of LightCurve
        The light curves to process.  Also accepts ``*args`` form::

            LightCurveBatch(lc1, lc2, lc3)
            LightCurveBatch([lc1, lc2, lc3])

    Command methods (``LS``, ``BLS``, ``rms``, …) are attached at import time
    by ``_method_gen._attach_all_command_methods()``.  Each returns a new
    ``LightCurveBatch`` with the command appended; call ``.run()`` to execute.

    Immediate methods (``run_LS``, ``run_BLS``, …) run a single command on
    every LC in the collection and return a ``BatchResult``.
    """

    def __init__(
        self,
        lcs: Union[Iterable["LightCurve"], "LightCurve"],
        *more_lcs: "LightCurve",
        _commands: Optional[List["VartoolsCommand"]] = None,
        _global_opts: Optional[Dict[str, Any]] = None,
        _prior_batch: Optional["BatchResult"] = None,
    ) -> None:
        from .lightcurve import LightCurve

        # Accept both LightCurveBatch([lc1, lc2]) and LightCurveBatch(lc1, lc2)
        if isinstance(lcs, LightCurve):
            all_lcs = [lcs] + list(more_lcs)
        else:
            all_lcs = list(lcs) + list(more_lcs)

        self._lcs: List["LightCurve"] = all_lcs
        self._commands: List["VartoolsCommand"] = list(_commands or [])
        self._global_opts: Dict[str, Any] = dict(_global_opts or {})
        self._prior_batch: Optional["BatchResult"] = _prior_batch

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _append(self, cmd: "VartoolsCommand") -> "LightCurveBatch":
        """Return a new LightCurveBatch with *cmd* appended."""
        return LightCurveBatch(
            self._lcs,
            _commands=self._commands + [cmd],
            _global_opts=self._global_opts,
            _prior_batch=self._prior_batch,
        )

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_files(cls, paths, **lc_kwargs) -> "LightCurveBatch":
        """Create a LightCurveBatch from a list of file paths.

        Parameters
        ----------
        paths : iterable of str or Path
            File paths to load.  Each is passed to ``LightCurve.from_file``.
        **lc_kwargs
            Extra keyword arguments forwarded to ``LightCurve.from_file``
            (e.g. ``format``, ``t_col``, ``mag_col``, ``err_col``).

        Examples
        --------
        ::

            import glob
            batch = vt.LightCurveBatch.from_files(
                sorted(glob.glob("/data/lcs/*.lc"))
            )
            result = batch.LS(0.5, 10.0, 1e-3).run()
        """
        from .lightcurve import LightCurve
        lcs = [LightCurve.from_file(p, **lc_kwargs) for p in paths]
        return cls(lcs)

    # ------------------------------------------------------------------
    # Global options
    # ------------------------------------------------------------------

    def with_options(self, **global_opts) -> "LightCurveBatch":
        """Return a new LightCurveBatch carrying pipeline-level global options.

        Options are forwarded to each ``Pipeline.run()`` call.

        Parameters
        ----------
        randseed, skipmissing, jdtol, matchstringid, timeout, capture_lc
        """
        merged = dict(self._global_opts)
        merged.update(global_opts)
        return LightCurveBatch(
            self._lcs,
            _commands=self._commands,
            _global_opts=merged,
            _prior_batch=self._prior_batch,
        )

    # ------------------------------------------------------------------
    # Per-LC parameter support
    # ------------------------------------------------------------------

    def _collect_perlc(self) -> Dict[tuple, Any]:
        """Return ``{(cmd_idx, attr_name): val}`` for per-LC array params.

        Scans every command in ``self._commands`` for attributes whose value
        is a ``PerLC``, a 1-D numpy array, or a pandas Series.  These are
        treated as per-LC arrays: the i-th element is used when running the
        i-th light curve.
        """
        found: Dict[tuple, Any] = {}
        for ci, cmd in enumerate(self._commands):
            for name, val in vars(cmd).items():
                if name.startswith("_"):
                    continue
                if _is_perlc_value(val):
                    found[(ci, name)] = val
        return found

    def _validate_perlc(self, perlc_map: Dict[tuple, Any]) -> None:
        """Raise ``ValueError`` if any per-LC array has the wrong length."""
        n = len(self._lcs)
        for (ci, name), val in perlc_map.items():
            if len(val) != n:
                cmd_name = self._commands[ci].__class__.__name__
                raise ValueError(
                    f"Per-LC parameter '{name}' on {cmd_name} has {len(val)} "
                    f"values but the batch has {n} light curves.  "
                    f"Each element of the array is used for one LC in order.  "
                    f"If you intended a fixed multi-valued parameter, wrap it "
                    f"in a plain Python list instead of a numpy array or Series."
                )

    def _resolve_perlc_for_lc(
        self,
        lc_idx: int,
        lc_name: str,
        perlc_map: Dict[tuple, Any],
    ) -> List["VartoolsCommand"]:
        """Return a command list with per-LC arrays substituted by scalars.

        Commands that carry no per-LC parameters are returned as-is (shared
        reference, not copied).  Commands that do carry per-LC parameters are
        shallow-copied and their per-LC attributes are replaced by the scalar
        value appropriate for *lc_idx*.
        """
        if not perlc_map:
            return self._commands

        # Group by command index for efficient lookup.
        per_cmd: Dict[int, Dict[str, Any]] = {}
        for (ci, attr_name), val in perlc_map.items():
            per_cmd.setdefault(ci, {})[attr_name] = val

        resolved: List["VartoolsCommand"] = []
        for ci, cmd in enumerate(self._commands):
            if ci not in per_cmd:
                resolved.append(cmd)
            else:
                cmd_copy = copy.copy(cmd)
                for attr_name, val in per_cmd[ci].items():
                    scalar = _extract_perlc_scalar(val, lc_idx, lc_name)
                    setattr(cmd_copy, attr_name, scalar)
                resolved.append(cmd_copy)
        return resolved

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run(
        self,
        capture_lc: Optional[bool] = None,
        timeout: Optional[int] = None,
        perpoint_vars=None,
        perlc_vars=None,
        randseed: Optional[int] = None,
        skipmissing: bool = False,
        jdtol: Optional[float] = None,
        matchstringid: bool = False,
    ) -> "BatchResult":
        """Execute the command chain on each LC.

        Parameters
        ----------
        capture_lc : bool, optional
            Whether to capture output light curves into per-LC ``result.lc``.
            Defaults to ``True``, or to the value set by ``with_options(capture_lc=...)``.
            An explicit ``run(capture_lc=...)`` call overrides ``with_options``.
        timeout : int, optional
            Maximum seconds per LC.
        perpoint_vars : dict mapping str to PerPointVar, optional
            Per-observation init expressions, forwarded to the underlying run.
        perlc_vars : dict, optional
            Per-LC variables.  Schema entries (``int`` / ``PerLCColumn``)
            and values entries (sequence of length ``len(self)``, or a
            ``(values, type)`` tuple) are both accepted; see
            :meth:`Pipeline.run_batch` for the dispatch rules.  When set,
            pyvartools routes through ``Pipeline.run_batch()`` (single
            vartools invocation) rather than the per-LC loop, so each LC
            sees its own value through the ``-inlistvars`` mechanism.
        randseed : int, optional
        skipmissing : bool
        jdtol : float, optional
        matchstringid : bool

        Returns
        -------
        BatchResult
            Supports both ``batch.vars`` (DataFrame) and ``batch[i]`` (Result)
            access.  ``batch[i].error`` is set if that LC failed.
        """
        from ._chain import CommandChain
        from .results import BatchResult, Result, RunError

        # Merge stored global_opts, then override with kwargs.  Booleans are
        # always written so that run(skipmissing=False) can clear a
        # with_options(skipmissing=True) default.
        # capture_lc is kept separate so it is never accidentally passed twice.
        run_kwargs: Dict[str, Any] = dict(self._global_opts)
        run_kwargs.pop("capture_lc", None)   # handled explicitly below
        if timeout is not None:
            run_kwargs["timeout"] = timeout
        if perpoint_vars is not None:
            run_kwargs["perpoint_vars"] = perpoint_vars
        if perlc_vars is not None:
            run_kwargs["perlc_vars"] = perlc_vars
        if randseed is not None:
            run_kwargs["randseed"] = randseed
        run_kwargs["skipmissing"] = skipmissing
        if jdtol is not None:
            run_kwargs["jdtol"] = jdtol
        run_kwargs["matchstringid"] = matchstringid

        # Resolve effective capture_lc.  Priority (highest first):
        #   1. explicit run(capture_lc=...) argument
        #   2. with_options(capture_lc=...) default
        #   3. built-in default: True
        if capture_lc is not None:
            effective_capture_lc = capture_lc
        else:
            effective_capture_lc = self._global_opts.get("capture_lc", True)

        # Detect and validate per-LC array parameters before the loop so that
        # a length mismatch is reported once rather than partway through.
        perlc_map = self._collect_perlc()
        if perlc_map:
            self._validate_perlc(perlc_map)

        prior_cmds: List[str] = list(
            self._prior_batch._known_commands or []
        ) if self._prior_batch is not None else []
        new_cmds: List[str] = [c._vt_name for c in self._commands]
        known_commands = prior_cmds + new_cmds

        results: List[Result] = []

        # Compute the command-offset for this segment: if this batch is a
        # continuation (_prior_batch set), the new commands' output-column
        # suffixes must pick up where the prior batch left off so they don't
        # collide with injected scalar names.
        segment_offset = len(prior_cmds) if self._prior_batch is not None else 0

        # Route through Pipeline.run_batch (single vartools invocation) when
        # this is a chain continuation OR when the user supplied perlc_vars.
        # Both cases need per-LC variables injected via -inlistvars, which
        # the per-LC loop below cannot do (each loop iteration is its own
        # vartools invocation).
        use_run_batch = segment_offset > 0 or perlc_vars is not None

        if use_run_batch:
            from .lightcurve import LightCurve
            from .pipeline import Pipeline
            input_lcs: List["LightCurve"] = []
            for i, lc in enumerate(self._lcs):
                carry = dict(lc.scalars)
                if self._prior_batch is not None:
                    prior_r = self._prior_batch._result_at(i)
                    for k, v in prior_r.vars.items():
                        if k != "Name" and k not in carry:
                            carry[k] = v
                input_lcs.append(
                    LightCurve(lc._df, name=lc.name, scalars=carry))
            batched = Pipeline(self._commands).run_batch(
                input_lcs,
                capture_lc=effective_capture_lc,
                raise_on_error=False,
                _command_offset=segment_offset,
                **run_kwargs,
            )
            if batched.error is not None:
                # A batched RunError covers the whole invocation.  Surface it
                # as one-per-LC errors so the Result-list shape stays stable.
                for _ in self._lcs:
                    results.append(Result(
                        var=pd.Series(dtype=object),
                        lc=None,
                        files={},
                        known_commands=new_cmds,
                        error=batched.error,
                    ))
            else:
                for i in range(len(self._lcs)):
                    results.append(batched._result_at(i))
        else:
            for i, lc in enumerate(self._lcs):
                try:
                    resolved_cmds = self._resolve_perlc_for_lc(i, lc.name, perlc_map)
                    input_lc = lc
                    chain = CommandChain(input_lc, resolved_cmds, self._global_opts)
                    r = chain.run(capture_lc=effective_capture_lc, **run_kwargs)
                    results.append(r)
                except Exception as exc:
                    from .results import RunError as _RE
                    err = exc if isinstance(exc, _RE) else _RE(str(exc))
                    results.append(Result(
                        var=pd.Series(dtype=object),
                        lc=None,
                        files={},
                        known_commands=new_cmds,
                        error=err,
                    ))

        # Merge prior per-LC vars into each new result so that output variables
        # from all preceding BatchResult steps are preserved.  The new result's
        # vars already carry correctly-shifted suffixes (Pipeline.run was
        # invoked with _command_offset=segment_offset), so we concatenate
        # directly rather than invoking _merge_prior, which would double-shift.
        if self._prior_batch is not None:
            merged: List[Result] = []
            for i, r in enumerate(results):
                prior_r = self._prior_batch._result_at(i)
                prior_vars = prior_r.vars.drop("Name", errors="ignore") \
                    if not prior_r.vars.empty else prior_r.vars
                merged_vars = (pd.concat([prior_vars, r.vars])
                               if not prior_vars.empty else r.vars)
                merged_known = (list(prior_r._known_commands or [])
                                + list(r._known_commands or []))
                merged.append(Result(
                    var=merged_vars,
                    lc=r.lc,
                    files=r.files,
                    known_commands=merged_known,
                    error=r.error,
                ))
            results = merged

        # Build the flat DataFrame from successful results
        rows = []
        lcs_out: List[Optional["LightCurve"]] = []
        all_files: Dict[str, list] = {}

        for r in results:
            if r.ok and isinstance(r.vars, pd.Series) and not r.vars.empty:
                rows.append(r.vars)
            else:
                rows.append(pd.Series(dtype=object))
            lcs_out.append(r.lc)
            for key, val in r.files.items():
                all_files.setdefault(key, []).append(val)

        if rows and any(not s.empty for s in rows):
            df = pd.DataFrame(rows).reset_index(drop=True)
            # Keep Name as the first column regardless of merge order.
            if "Name" in df.columns:
                df = df[["Name"] + [c for c in df.columns if c != "Name"]]
        else:
            df = pd.DataFrame()

        return BatchResult(
            var=df,
            lcs=lcs_out if effective_capture_lc else None,
            files=all_files,
            error=None,
            known_commands=known_commands,
            _results=results,
        )

    # ------------------------------------------------------------------
    # Sequence interface
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return len(self._lcs)

    def __iter__(self):
        return iter(self._lcs)

    def __getitem__(self, key) -> "LightCurve":
        """Return a LightCurve by integer position or by ``name``.

        Examples
        --------
        >>> batch[0]            # first LC
        >>> batch['hat-123']    # LC whose .name == 'hat-123'

        Raises
        ------
        KeyError
            If a string *key* does not match any ``LightCurve.name`` in
            the batch.
        """
        if isinstance(key, str):
            for lc in self._lcs:
                if lc.name == key:
                    return lc
            raise KeyError(
                f"No LightCurve in batch with name={key!r}.  "
                f"Available: {[lc.name for lc in self._lcs]}"
            )
        return self._lcs[key]

    def __contains__(self, item) -> bool:
        if isinstance(item, str):
            return any(lc.name == item for lc in self._lcs)
        return item in self._lcs

    # ------------------------------------------------------------------
    # Repr
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        cmds = " → ".join(c.__class__.__name__ for c in self._commands)
        cmd_str = f"[{cmds}]" if cmds else "[]"
        return f"LightCurveBatch(n={len(self._lcs)}, commands={cmd_str})"
