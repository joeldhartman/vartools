"""Shared helper functions for building vartools CLI argument lists."""

import re
from typing import List, Optional, Union


def _flag(name: str, val) -> List[str]:
    """Return [name, str(val)] if val is not None, else []."""
    return [name, str(val)] if val is not None else []


def _bool(name: str, val: bool) -> List[str]:
    """Return [name] if val is True, else []."""
    return [name] if val else []


def _norm_save(spec):
    """Normalise a save_* parameter value to an Output instance.

    Accepted forms
    --------------
    ``False`` or ``None``   → ``Output(path=None, capture=False)``  (Mode 4)
    ``True``                → ``Output(path=None, capture=True)``   (Mode 1)
    ``"/path/to/dir"``      → ``Output(path=str, capture=False)``   (Mode 3)
    ``Output(...)``         → returned as-is
    """
    from pyvartools._output import Output
    if spec is False or spec is None:
        return Output(path=None, capture=False)
    if spec is True:
        return Output(path=None, capture=True)
    if isinstance(spec, str):
        return Output(path=spec, capture=False)
    if isinstance(spec, Output):
        return spec
    raise TypeError(f"Invalid save spec: {spec!r}")


def _should_emit(spec) -> bool:
    """Return True if vartools should write this output file.

    True when capture=True (Mode 1/2) or path is set (Mode 2/3).
    False only for Mode 4 (capture=False, path=None).
    """
    s = _norm_save(spec)
    return s.capture or (s.path is not None)


def _outtoken(save_spec, fallback_outdir: str) -> List[str]:
    """Return ['1', actual_outdir] if the output should be emitted, else ['0'].

    actual_outdir is ``spec.path`` if the user specified a directory, else
    ``fallback_outdir`` (the pipeline-managed temp dir for this command).

    Backward compatible: ``_outtoken(True, outdir)`` and
    ``_outtoken(False, outdir)`` behave the same as before.
    """
    s = _norm_save(save_spec)
    if not (s.capture or s.path is not None):
        return ["0"]
    actual = s.path if s.path is not None else fallback_outdir
    return ["1", actual]


def _period_spec(period) -> List[str]:
    """Convert a period value to vartools period-spec tokens.

    Parameters
    ----------
    period : float, int, or str
        - A number → ``["fix", str(period)]``
        - A string whose first token is a known vartools keyword (``"ls"``,
          ``"aov"``, ``"bls"``, ``"fix 1.234"``, ``"list column 2"``,
          ``"fixcolumn COL"``, ``"var NAME"``, ``"expr EXPR"``, ...) →
          split and passed through verbatim.
        - A bare identifier (e.g. ``"myperiod"``) → ``["var", "myperiod"]``
        - Any other expression (e.g. ``"2*P"``) → ``["expr", "2*P"]``
    """
    if period is None:
        return []
    if isinstance(period, (int, float)):
        return ["fix", str(period)]
    s = str(period)
    if (s.split()[0] in ("ls", "aov", "bls", "both", "injectharm",
                          "rand", "logrand", "randfreq", "lograndfreq",
                          "fix", "list", "fixcolumn", "auto", "var", "expr")):
        return s.split()
    if re.match(r'^[A-Za-z_]\w*$', s):
        return ["var", s]
    return ["expr", s]


def _auto_or_varexpr(val) -> List[str]:
    """Like :func:`_varexpr`, but pass the bare keyword ``"auto"`` through
    verbatim (used by commands such as ``-wwz`` where ``auto`` is a literal
    CLI option, not a variable name)."""
    if isinstance(val, str) and val == "auto":
        return ["auto"]
    return _varexpr(val)


def _varexpr(val) -> List[str]:
    """Convert a numeric parameter that may be a variable name or expression.

    Parameters
    ----------
    val : int, float, or str
        - A number → ``[str(val)]`` (passed directly to vartools)
        - A bare identifier string (``[A-Za-z_]\\w*``) →
          ``["var", val]`` (vartools ``var`` keyword: reads a named variable)
        - A string starting with ``"list"`` or ``"fixcolumn"`` →
          split and passed through verbatim (e.g. ``"list column 2"`` →
          ``["list", "column", "2"]``)
        - Any other string → ``["expr", val]``
          (vartools ``expr`` keyword: evaluates a math expression)
    """
    if isinstance(val, (int, float)):
        return [str(val)]
    s = str(val)
    # Pass explicit keywords through verbatim (must be checked before the
    # bare-identifier test since "list" also matches that pattern).
    if re.match(r'^(list|fixcolumn|var|expr)\b', s):
        return s.split()
    if re.match(r'^[A-Za-z_]\w*$', s):
        return ["var", s]
    return ["expr", s]


def _fixperiodsnr_tokens(val) -> List[str]:
    """Convert a fixperiodSNR parameter to vartools CLI tokens.

    - None → [] (feature disabled)
    - float/int → ["fixperiodSNR", "fix", str(val)]
    - "aov" / "ls" / "injectharm" → ["fixperiodSNR", val]
    - "fixcolumn COLNAME" → ["fixperiodSNR", "fixcolumn", "COLNAME"]
    - "list" / "list column 2" → ["fixperiodSNR", "list", ...]
    """
    if val is None:
        return []
    if isinstance(val, (int, float)):
        return ["fixperiodSNR", "fix", str(val)]
    return ["fixperiodSNR"] + str(val).split()


def _pval(v, keyword: Optional[str] = None) -> List[str]:
    """Convert a Python value to a vartools fix/var/expr/bare parameter.

    Parameters
    ----------
    v : number or str or None
        - None → []
        - A number → [str(v)] (or [keyword, str(v)] if keyword given)
        - A bare identifier string → ["var", v]
        - An expression string → ["expr", v]
        - A string starting with fix/list/fixcolumn/var/expr → split as-is
    keyword : str, optional
        If given, prepend it when v is a bare number.
    """
    if v is None:
        return []
    if isinstance(v, (int, float)):
        return ([keyword] if keyword else []) + [str(v)]
    s = str(v)
    if re.match(r'^(fix|list|fixcolumn|var|expr|ls|aov|bls|both|auto|injectharm)\b', s):
        return s.split()
    if re.match(r'^[A-Za-z_]\w*$', s):
        return ["var", s]
    return ["expr", s]


def _injectparam(prefix: str, val) -> List[str]:
    """Build an Injecttransit per-parameter token sequence.

    - float/int → [f"{prefix}fix", str(val)]
    - bare identifier string → [f"{prefix}var", val]
    - expression string → [f"{prefix}expr", val]
    - str starting with a keyword → split as-is
      (e.g. "Rplogrand 0.05 0.15" or "phaserand")
    """
    if isinstance(val, (int, float)):
        return [f"{prefix}fix", str(val)]
    s = str(val)
    if re.match(r'^[A-Za-z_]\w*$', s) and not s.endswith("rand"):
        return [f"{prefix}var", s]
    if not re.match(r'^[A-Za-z]', s) or any(c in s for c in '+-*/()'):
        return [f"{prefix}expr", s]
    return s.split()


# ---------------------------------------------------------------------------
# Back-reference resolution helpers (for chained-command continuations)
#
# These take a prior `Result` or `BatchResult` and return a Python-side
# numeric value (or PerLC array, for batch) that replaces the CLI back-ref
# keyword.  Each resolver is tolerant of the case where *no* prior command
# of the requested type exists — it raises `LookupError` with a message that
# tells the user which chain step is missing.
# ---------------------------------------------------------------------------

_FIXCOLUMN_RE = re.compile(r'^\s*fixcolumn\s+(\S+)\s*$')


def _is_perstar_value(val) -> bool:
    """Detect per-star result values (arrays/Series) returned by a BatchResult."""
    try:
        import pandas as pd  # local import — keeps cold-start light
    except Exception:
        pd = None
    if pd is not None and isinstance(val, (pd.Series,)):
        return True
    try:
        import numpy as np
        if isinstance(val, np.ndarray):
            return val.ndim >= 1
    except Exception:
        pass
    return False


def _coerce_to_numeric(val):
    """Coerce a single-LC scalar or per-star vector to a form usable in _to_cli_args.

    - Single-LC: return ``float(val)``.
    - Per-star (batch): return a ``PerLC`` whose values are floats.
    """
    if _is_perstar_value(val):
        from pyvartools.perlc import PerLC
        return PerLC([float(x) for x in val])
    return float(val)


def _most_recent_lookup(prev, candidates):
    """Find the most-recent prior command of any of the given types.

    Parameters
    ----------
    prev : Result or BatchResult
        The continuation's prior step.  Must expose ``varobjs`` with
        ``CommandStats`` entries carrying a ``_position`` command-index.
    candidates : list[str]
        Command names to consider, e.g. ``["aov", "aov_harm"]``.  Any command
        name present in ``prev.varobjs`` whose ``_position`` is the greatest
        is returned.

    Returns
    -------
    CommandStats or None
        The matching stats namespace with the highest ``_position``, or None
        if none of *candidates* exists in *prev*.
    """
    best_pos = -1
    best_stats = None
    varobjs = getattr(prev, "varobjs", None)
    if varobjs is None:
        return None
    for name in candidates:
        ns = getattr(varobjs, name, None)
        if ns is None:
            continue
        # Single occurrence is a CommandStats; multiple becomes a list.
        entries = ns if isinstance(ns, list) else [ns]
        for stats in entries:
            pos = getattr(stats, "_position", -1)
            if pos > best_pos:
                best_pos = pos
                best_stats = stats
    return best_stats


def _resolve_period_backref(prev, spec):
    """Resolve a period-spec keyword against *prev*.

    Recognised specs: ``"ls"``, ``"aov"`` (matches either -aov or -aov_harm,
    most-recent wins), ``"pdm"``, ``"bls"``, ``"blsfixper"``, ``"injectharm"``,
    and ``"fixcolumn <name>"``.  All other strings pass through unchanged
    (the caller may still want to accept e.g. ``"fix 1.23"``).

    Returns a float / PerLC / the original value when the spec isn't a
    recognised back-ref.  Raises ``LookupError`` if a recognised back-ref
    has no prior counterpart in *prev*.
    """
    if spec is None or not isinstance(spec, str):
        return spec
    s = spec.strip()

    if s == "ls":
        stats = _most_recent_lookup(prev, ["LS"])
        if stats is None:
            raise LookupError(
                "Back-reference 'ls' has no prior -LS command in this chain"
            )
        return _coerce_to_numeric(stats.Period_1)

    if s == "aov":
        stats = _most_recent_lookup(prev, ["aov", "aov_harm"])
        if stats is None:
            raise LookupError(
                "Back-reference 'aov' has no prior -aov or -aov_harm command "
                "in this chain"
            )
        return _coerce_to_numeric(stats.Period_1)

    if s == "pdm":
        stats = _most_recent_lookup(prev, ["PDM"])
        if stats is None:
            raise LookupError(
                "Back-reference 'pdm' has no prior -PDM command in this chain"
            )
        return _coerce_to_numeric(stats.Period_1)

    if s == "bls":
        stats = _most_recent_lookup(prev, ["BLS"])
        if stats is None:
            raise LookupError(
                "Back-reference 'bls' has no prior -BLS command in this chain"
            )
        return _coerce_to_numeric(stats.Period_1)

    if s == "blsfixper":
        stats = _most_recent_lookup(prev, ["BLSFixPer"])
        if stats is None:
            raise LookupError(
                "Back-reference 'blsfixper' has no prior -BLSFixPer command "
                "in this chain"
            )
        return _coerce_to_numeric(stats.Period)

    if s == "injectharm":
        stats = _most_recent_lookup(prev, ["Injectharm"])
        if stats is None:
            raise LookupError(
                "Back-reference 'injectharm' has no prior -Injectharm command "
                "in this chain"
            )
        return _coerce_to_numeric(stats.Period)

    m = _FIXCOLUMN_RE.match(s)
    if m:
        return _resolve_fixcolumn(prev, m.group(1))

    # Not a back-reference we recognise — pass through so existing handling
    # (e.g. ``"fix 1.23"``, ``"2*tspan"``) continues to work.
    return spec


def _resolve_fixcolumn(prev, col):
    """Resolve ``fixcolumn NAME`` against *prev*.

    *col* is a bare column name as the user typed it.  Returns the column's
    value (float for single-LC, PerLC for batch).  Bare integers are not
    supported — they would require reconstructing vartools' column-number
    assignment.  Raises ``ValueError`` if *col* is numeric, or ``LookupError``
    if the named column is not present.
    """
    if col.lstrip("-").isdigit():
        raise ValueError(
            "fixcolumn with a numeric column index is not supported in chained "
            "calls; pass the column name instead, e.g. "
            "'fixcolumn LS_Period_1_0'."
        )
    vars_obj = getattr(prev, "vars", None)
    if vars_obj is None or col not in vars_obj:
        raise LookupError(
            f"Back-reference 'fixcolumn {col}' refers to a column that is not "
            f"present in the prior step's output."
        )
    return _coerce_to_numeric(vars_obj[col])


def _resolve_bls_transit_backref(prev, spec):
    """Resolve ``"bls"`` / ``"blsfixper"`` for full transit-parameter spawning.

    Returns a dict with keys ``period``, ``T0``, ``depth``, ``qtran`` pulled
    from the most-recent -BLS (or -BLSFixPer) in *prev*.  Each value is a
    float (single-LC) or a PerLC (batch).
    """
    if spec == "bls":
        stats = _most_recent_lookup(prev, ["BLS"])
        if stats is None:
            raise LookupError(
                "Back-reference 'bls' has no prior -BLS command in this chain"
            )
        return {
            "period": _coerce_to_numeric(stats.Period_1),
            "T0":     _coerce_to_numeric(stats.Tc_1),
            "depth":  _coerce_to_numeric(stats.Depth_1),
            "qtran":  _coerce_to_numeric(stats.Qtran_1),
        }
    if spec == "blsfixper":
        stats = _most_recent_lookup(prev, ["BLSFixPer"])
        if stats is None:
            raise LookupError(
                "Back-reference 'blsfixper' has no prior -BLSFixPer command "
                "in this chain"
            )
        return {
            "period": _coerce_to_numeric(stats.Period),
            "T0":     _coerce_to_numeric(stats.Tc),
            "depth":  _coerce_to_numeric(stats.Depth),
            "qtran":  _coerce_to_numeric(stats.Qtran),
        }
    raise ValueError(f"_resolve_bls_transit_backref: unknown spec {spec!r}")
