"""
Result objects returned by Pipeline.run() and Pipeline.run_batch().
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd

from .lightcurve import LightCurve
from ._vars import VarsNamespace


class RunError(RuntimeError):
    """Raised when vartools exits with a non-zero status or writes to stderr."""
    pass


class Result:
    """Output from a single-light-curve pipeline run.

    Attributes
    ----------
    vars : pd.Series
        Series indexed by vartools output variable names
        (e.g. ``'LS_Period_1_0'``, ``'Log10_LS_Prob_1_0'``, ...) plus ``'Name'``.
        Always a ``pd.Series``; values are scalars, so
        ``result.vars["LS_Period_1_0"]`` returns a scalar directly.
    lc : LightCurve or None
        The (possibly modified) output light curve.  None if the pipeline
        does not modify the light curve (no -o or -os given).
    files : dict
        Mapping of logical name → DataFrame for any output files that were
        requested via ``save_*=True`` on individual commands.
        Empty dict if no output files were requested.
    """

    def __init__(
        self,
        var: Union[pd.Series, pd.DataFrame],
        lc: Optional[LightCurve] = None,
        files: Optional[Dict[str, pd.DataFrame]] = None,
        known_commands: Optional[List[str]] = None,
        error: Optional["RunError"] = None,
    ) -> None:
        # Normalize to always store a pd.Series.  The only time a DataFrame
        # arrives here is when parse_oneline_output returns an empty result
        # and the caller skips the .iloc[0] step.
        if isinstance(var, pd.DataFrame):
            self.vars: pd.Series = var.iloc[0] if not var.empty else pd.Series(dtype=object)
        else:
            self.vars = var
        self.lc = lc
        self.files = files or {}
        self._known_commands = known_commands
        self.error = error

    @property
    def ok(self) -> bool:
        """True if this result completed without error."""
        return self.error is None

    @property
    def varobjs(self) -> VarsNamespace:
        """Structured per-command access to output variables.

        Examples
        --------
        ::

            result.varobjs.LS.Period_1          # single LS call
            result.varobjs.LS[0].Period_1       # explicit index
            result.varobjs.BLS[0].SDE_0
        """
        return VarsNamespace.from_series(self.vars, self._known_commands)

    def __getattr__(self, name: str) -> object:
        """Shorthand access to individual var keys, e.g. result.LS_Period_1_0."""
        # Only reached when normal attribute lookup has failed.
        vars_ = object.__getattribute__(self, "vars")
        if isinstance(vars_, pd.Series) and name in vars_.index:
            return vars_[name]
        raise AttributeError(
            f"'Result' object has no attribute {name!r}. "
            f"Use result.vars[{name!r}] or result.varobjs.<command>.<key>."
        )

    def __repr__(self) -> str:
        keys = list(self.vars.index)
        return (f"Result(var_keys={keys}, "
                f"lc={'yes' if self.lc else 'no'}, "
                f"files={list(self.files.keys())})")

    def _repr_html_(self) -> str:
        """HTML repr for Jupyter notebooks."""
        name = self.vars.get("Name", "—") if hasattr(self.vars, "get") else "—"
        n_vars = max(0, len(self.vars) - (1 if "Name" in self.vars.index else 0))
        lc_note = (
            f"{len(self.lc)} points" if self.lc is not None else "not captured"
        )
        files_note = (
            ", ".join(self.files.keys()) if self.files else "none"
        )
        status = (
            '<span style="color:#c00">FAILED</span>'
            if self.error is not None else
            '<span style="color:#080">ok</span>'
        )

        # Vars table: drop "Name" since it's already in the header
        vars_for_table = (
            self.vars.drop("Name") if "Name" in self.vars.index else self.vars
        )
        vars_html = (
            vars_for_table.to_frame("value")._repr_html_()
            if len(vars_for_table) > 0
            else "<i>no variables</i>"
        )

        return (
            f'<div style="font-family:sans-serif;font-size:0.9em">'
            f'  <b>Result</b> &middot; '
            f'  Name: <code>{name}</code> &middot; '
            f'  vars: {n_vars} &middot; '
            f'  lc: {lc_note} &middot; '
            f'  files: {files_note} &middot; '
            f'  status: {status}'
            f'  <details style="margin-top:0.4em"><summary>vars</summary>'
            f'    {vars_html}'
            f'  </details>'
            f'</div>'
        )


class BatchResult:
    """Output from a multi-light-curve pipeline run.

    Attributes
    ----------
    vars : pd.DataFrame
        One row per light curve; columns are vartools output variable names.
    lcs : list of LightCurve
        The (possibly modified) output light curves, in input order.
        Empty list if light-curve capture was not requested.
    files : dict
        Mapping of logical name → list of DataFrames (one per light curve)
        for any output files that were requested.
    error : RunError or None
        If the entire run failed (e.g. for ``Pipeline.run_batch()``), the
        exception is stored here.  For per-LC errors from
        ``LightCurveBatch.run()``, check ``batch[i].error`` instead.

    Per-LC access
    -------------
    ``batch[i]``              — i-th ``Result`` object
    ``for result in batch:``  — iterate as individual ``Result`` objects
    ``batch[i].error``        — ``RunError`` if that LC failed, else ``None``
    """

    def __init__(
        self,
        var: pd.DataFrame,
        lcs: Optional[list] = None,
        files: Optional[Dict[str, list]] = None,
        error: Optional["RunError"] = None,
        known_commands: Optional[List[str]] = None,
        _results: Optional[List["Result"]] = None,
    ) -> None:
        self.vars = var
        self._lcs: Optional[List[Optional[LightCurve]]] = lcs  # None = not captured
        self.files = files or {}
        self.error = error
        self._known_commands = known_commands
        self._results = _results  # pre-built Results from LightCurveBatch.run()

    @property
    def ok(self) -> bool:
        """True if the overall run completed without error."""
        return self.error is None

    @property
    def lcs(self) -> List[Optional[LightCurve]]:
        """Processed output light curves, one per input LC.

        Returns an empty list when ``capture_lc=False`` was used (rather
        than ``None``), so ``for lc in batch.lcs:`` is always safe.
        Check ``batch._lcs is not None`` (or re-run with ``capture_lc=True``)
        if you need to distinguish "not captured" from "empty batch".
        """
        return self._lcs if self._lcs is not None else []

    @property
    def lcscalars(self) -> pd.DataFrame:
        """Per-star scalar variables as a DataFrame, one row per LC.

        A convenience view over the captured LCs' ``.scalars`` dicts —
        equivalent to ``pd.DataFrame([lc.scalars for lc in batch.lcs])``.
        The values live on the ``LightCurve`` objects in ``batch.lcs``;
        this property just aggregates them.

        Columns are scalar variable names (no ``_N`` suffix), rows follow
        the input LC order.  Returns an empty DataFrame when no LCs have
        any scalars (e.g. when ``capture_lc=False`` or before chained
        scalar injection is used).
        """
        lcs = self._lcs if self._lcs is not None else []
        rows = []
        for lc in lcs:
            if lc is None or not lc.scalars:
                rows.append({})
            else:
                rows.append(dict(lc.scalars))
        if not any(rows):
            return pd.DataFrame()
        return pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # Per-LC Result access
    # ------------------------------------------------------------------

    def _result_at(self, i: int) -> "Result":
        """Return the Result for the i-th light curve."""
        # If pre-built Results are stored (LightCurveBatch path), use them.
        if self._results is not None:
            return self._results[i]

        # Reconstruct a Result from the i-th DataFrame row.
        if self.vars.empty or i >= len(self.vars):
            return Result(
                var=pd.Series(dtype=object),
                lc=None,
                files={},
                known_commands=self._known_commands,
            )
        row = self.vars.iloc[i]
        lc = self.lcs[i] if self.lcs and i < len(self.lcs) else None
        per_lc_files = {k: v[i] for k, v in self.files.items()
                        if isinstance(v, list) and i < len(v)}
        return Result(
            var=row,
            lc=lc,
            files=per_lc_files,
            known_commands=self._known_commands,
        )

    def __getitem__(self, i):
        """Return the i-th Result, or a sub-BatchResult for a slice.

        Parameters
        ----------
        i : int or slice
            An integer index returns a single ``Result``.
            A slice returns a new ``BatchResult`` containing the selected LCs.

        Examples
        --------
        ::

            first = batch[0]           # Result
            last  = batch[-1]          # Result
            sub   = batch[1:4]         # BatchResult with LCs 1, 2, 3
            every_other = batch[::2]   # BatchResult
        """
        if isinstance(i, slice):
            indices = list(range(*i.indices(len(self))))
            new_vars = self.vars.iloc[indices].reset_index(drop=True) if not self.vars.empty else pd.DataFrame()
            new_lcs = [self._lcs[j] for j in indices] if self._lcs is not None else None
            new_results = ([self._results[j] for j in indices]
                           if self._results is not None else None)
            new_files = {
                k: [v[j] for j in indices]
                for k, v in self.files.items()
                if isinstance(v, list)
            }
            return BatchResult(
                var=new_vars,
                lcs=new_lcs,
                files=new_files,
                error=None,
                known_commands=self._known_commands,
                _results=new_results,
            )
        return self._result_at(i)

    def __iter__(self):
        for i in range(len(self)):
            yield self._result_at(i)

    def __len__(self) -> int:
        if self._results is not None:
            return len(self._results)
        return len(self.vars)

    # ------------------------------------------------------------------
    # Filtering
    # ------------------------------------------------------------------

    def filter(self, mask) -> "BatchResult":
        """Return a sub-BatchResult containing only the rows where *mask* is True.

        Parameters
        ----------
        mask : array-like of bool or pd.Series of bool
            One entry per light curve.  May be a boolean numpy array, a Python
            list of bools, or a boolean ``pd.Series`` (e.g. from a DataFrame
            comparison on ``batch.vars``).

        Returns
        -------
        BatchResult

        Examples
        --------
        ::

            # Keep only LCs with a significant LS detection
            significant = batch.filter(batch.vars["Log10_LS_Prob_1_0"] < -5)

            # Keep only the first three
            import numpy as np
            subset = batch.filter(np.array([True, True, True, False, False]))
        """
        mask_arr = np.asarray(mask, dtype=bool)
        indices = [i for i, m in enumerate(mask_arr) if m]

        new_vars = self.vars.iloc[indices].reset_index(drop=True) if not self.vars.empty else pd.DataFrame()
        new_lcs = [self._lcs[i] for i in indices] if self._lcs is not None else None
        new_results = ([self._results[i] for i in indices]
                       if self._results is not None else None)
        new_files = {
            k: [v[i] for i in indices]
            for k, v in self.files.items()
            if isinstance(v, list)
        }
        return BatchResult(
            var=new_vars,
            lcs=new_lcs,
            files=new_files,
            error=None,
            known_commands=self._known_commands,
            _results=new_results,
        )

    # ------------------------------------------------------------------
    # Repr
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        return (f"BatchResult(n={len(self)}, "
                f"cols={list(self.vars.columns) if not self.vars.empty else []}, "
                f"lcs={'yes' if self._lcs is not None else 'no'}, "
                f"ok={self.ok}, "
                f"files={list(self.files.keys())})")

    def _repr_html_(self) -> str:
        """HTML repr for Jupyter notebooks."""
        n = len(self)
        ncols = len(self.vars.columns) if not self.vars.empty else 0
        lc_note = (
            f"{sum(1 for x in self._lcs if x is not None)} captured"
            if self._lcs is not None else "not captured"
        )
        files_note = (
            ", ".join(self.files.keys()) if self.files else "none"
        )
        status = (
            '<span style="color:#c00">FAILED</span>'
            if self.error is not None else
            '<span style="color:#080">ok</span>'
        )

        # Show first/last rows of vars table; for short tables show all.
        if not self.vars.empty:
            if n <= 12:
                vars_html = self.vars._repr_html_()
            else:
                head = self.vars.head(5)
                tail = self.vars.tail(5)
                ellipsis = pd.DataFrame(
                    [["..."] * ncols],
                    columns=self.vars.columns,
                    index=["..."],
                )
                shown = pd.concat([head, ellipsis, tail])
                vars_html = shown._repr_html_()
        else:
            vars_html = "<i>no variables</i>"

        return (
            f'<div style="font-family:sans-serif;font-size:0.9em">'
            f'  <b>BatchResult</b> &middot; '
            f'  n: {n} &middot; '
            f'  cols: {ncols} &middot; '
            f'  lcs: {lc_note} &middot; '
            f'  files: {files_note} &middot; '
            f'  status: {status}'
            f'  <details style="margin-top:0.4em" open><summary>vars</summary>'
            f'    {vars_html}'
            f'  </details>'
            f'</div>'
        )


# ---------------------------------------------------------------------------
# Result-merging helpers (used when chaining from a Result)
# ---------------------------------------------------------------------------

_TRAILING_POS = re.compile(r'^(.*?)_(\d+)$')


def _shift_var_positions(series: pd.Series, offset: int) -> pd.Series:
    """Return a copy of *series* with every trailing ``_N`` position suffix
    incremented by *offset*.  The ``Name`` key is left unchanged.

    Used when merging a prior ``Result``'s vars with a newly produced one so
    that position indices remain unique and reflect the true command order.
    """
    if offset == 0:
        return series
    new_keys = []
    for key in series.index:
        if key == "Name":
            new_keys.append(key)
            continue
        m = _TRAILING_POS.match(key)
        if m:
            new_keys.append(f"{m.group(1)}_{int(m.group(2)) + offset}")
        else:
            new_keys.append(key)
    shifted = series.copy()
    shifted.index = new_keys
    return shifted


def _merge_prior(prior: "Result", new: "Result") -> "Result":
    """Return a new ``Result`` whose vars combine *prior* and *new*.

    The *new* result's position suffixes are shifted by
    ``len(prior._known_commands)`` so they do not collide with the prior
    result's positions.  The output ``lc``, ``files``, and ``error`` come
    from *new* (it holds the latest light curve state).

    If *prior* has no known_commands or an empty vars Series the merge is a
    no-op and *new* is returned unchanged.
    """
    prior_cmds = prior._known_commands or []
    offset = len(prior_cmds)

    if offset == 0 or prior.vars.empty:
        return new

    # Shift position suffixes in the newly produced result
    shifted_new_vars = _shift_var_positions(new.vars, offset)

    # Merge: prior vars (Name dropped — new's Name is authoritative) + shifted new vars
    prior_data = prior.vars.drop("Name", errors="ignore")
    merged_vars = pd.concat([prior_data, shifted_new_vars])

    merged_known_commands = list(prior_cmds) + list(new._known_commands or [])

    return Result(
        var=merged_vars,
        lc=new.lc,
        files=new.files,
        known_commands=merged_known_commands,
        error=new.error,
    )


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

# Lines starting with this prefix are per-star scalar variables emitted by
# vartools' ``-printallscalars`` option.  They are routed to
# ``LightCurve.scalars`` (and ``BatchResult.lcscalars`` for batch runs)
# rather than ``.vars``.
_SCALAR_PREFIX = "VARTOOLS_SCALAR:"


def split_vars_and_scalars(df: pd.DataFrame) -> tuple:
    """Split a parsed ``-oneline`` DataFrame into (vars_df, scalars_df).

    Columns whose name starts with :data:`_SCALAR_PREFIX` are moved to
    ``scalars_df`` with the prefix stripped.  All other columns (including
    ``Name``) stay in ``vars_df``.  Both DataFrames share the same row
    index and row count so per-LC correspondence is preserved.

    Returns two DataFrames even if one side is empty, so callers can treat
    the two uniformly.
    """
    if df.empty:
        return df, pd.DataFrame(index=df.index)
    scalar_cols = [c for c in df.columns if c.startswith(_SCALAR_PREFIX)]
    if not scalar_cols:
        return df, pd.DataFrame(index=df.index)
    vars_df = df.drop(columns=scalar_cols)
    scalars_df = df[scalar_cols].copy()
    scalars_df.columns = [c[len(_SCALAR_PREFIX):] for c in scalar_cols]
    return vars_df, scalars_df


def parse_oneline_output(stdout: str) -> pd.DataFrame:
    """Parse vartools ``-oneline`` stdout into a DataFrame.

    Each output line has the form::

        Name = EXAMPLES/2
        LS_Period_1_0 = 1.23440877
        ...
        Name = EXAMPLES/3
        ...

    Returns a DataFrame with one row per light curve.
    """
    rows = []
    current: Dict[str, object] = {}

    for line in stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        m = re.match(r'^([^\s=]+)\s*=\s*(.*)', line)
        if not m:
            continue
        key, val = m.group(1), m.group(2).strip()
        if key == "Name" and current:
            rows.append(current)
            current = {}
        current[key] = _coerce(val)

    if current:
        rows.append(current)

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(rows)


def _coerce(val: str):
    """Try to convert a string value to int, then float, then leave as str."""
    try:
        return int(val)
    except ValueError:
        pass
    try:
        return float(val)
    except ValueError:
        pass
    return val
