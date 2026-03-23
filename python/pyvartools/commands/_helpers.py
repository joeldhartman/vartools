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
        - A string like ``"ls"``, ``"aov"``, ``"bls"`` → split and passed as-is
        - A string like ``"fix 1.234"`` → split as-is
    """
    if period is None:
        return []
    if isinstance(period, (int, float)):
        return ["fix", str(period)]
    return str(period).split()


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
    # Pass "list [column N]" and "fixcolumn <name>" through verbatim
    # (must be checked before the bare-identifier test since "list" matches both)
    if re.match(r'^(list|fixcolumn)\b', s):
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
        - A string → split on whitespace (pass through var/expr/list specs)
    keyword : str, optional
        If given, prepend it when v is a bare number.
    """
    if v is None:
        return []
    if isinstance(v, (int, float)):
        return ([keyword] if keyword else []) + [str(v)]
    return str(v).split()


def _injectparam(prefix: str, val) -> List[str]:
    """Build an Injecttransit per-parameter token sequence.

    - float/int → [f"{prefix}fix", str(val)]
    - str       → str(val).split()  (already contains the prefix keyword,
                   e.g. "Rplogrand 0.05 0.15" or "phaserand")
    """
    if isinstance(val, (int, float)):
        return [f"{prefix}fix", str(val)]
    return str(val).split()
