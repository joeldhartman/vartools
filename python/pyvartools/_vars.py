"""
VarsNamespace: structured per-command access to vartools output variables.

Given a flat ``result.vars`` Series like::

    LS_Period_1_0           1.234
    Log10_LS_Prob_1_0      -4.56
    LS_Periodogram_Value_1_0  0.98
    Weighted_RMS_1          0.005

``result.varobjs`` exposes::

    result.varobjs.LS.Period_1            # LS at position 0
    result.varobjs.LS[0].Period_1         # explicit index
    result.varobjs.rms.Weighted_RMS       # rms at position 1

Key parsing
-----------
Vartools output column names follow the pattern::

    PREFIX_descriptor_N

where N is the 0-based command position in the pipeline.  The trailing ``_N``
position suffix is stripped to produce the attribute key.  The leading command
prefix is also stripped when the key starts with ``COMMANDNAME_`` (case-
insensitive).  Keys with a leading qualifier (e.g. ``Log10_LS_Prob_1_0``) keep
the qualifier in their attribute name.

The mapping from position N to command name is most reliable when the Pipeline
passes ``known_commands``.  Without it, command names are inferred from the
key prefixes, which is a best-effort heuristic.
"""

from __future__ import annotations

import re
from typing import Dict, Iterator, List, Optional, Union

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_TRAILING_INT = re.compile(r"^(.*?)_(\d+)$")


def _parse_position(key: str) -> Optional[int]:
    """Return the trailing integer position from a vartools column name, or None."""
    m = _TRAILING_INT.match(key)
    if m:
        try:
            return int(m.group(2))
        except ValueError:
            pass
    return None


def _add_peak_arrays(stripped: Dict[str, object]) -> Dict[str, object]:
    """Augment *stripped* with array forms of multi-peak descriptors.

    For every group of keys ``<base>_1, <base>_2, ..., <base>_N`` whose peak
    indices form a consecutive 1..N sequence and where ``<base>`` is not
    already a key, add ``<base>`` as a length-N numpy array of the values.

    The original ``<base>_i`` scalars are kept for backwards compatibility.
    """
    # Group keys by base prefix (everything before a trailing _<int>).
    groups: Dict[str, Dict[int, object]] = {}
    for key, val in stripped.items():
        m = _TRAILING_INT.match(key)
        if not m:
            continue
        base, peak_str = m.group(1), m.group(2)
        try:
            peak = int(peak_str)
        except ValueError:
            continue
        if peak < 1:
            continue
        groups.setdefault(base, {})[peak] = val

    out = dict(stripped)
    for base, peaks in groups.items():
        if base in stripped:
            # A bare <base> already exists — don't overwrite it.
            continue
        # Require a consecutive 1..N sequence (single-peak groups also qualify).
        n = max(peaks)
        if set(peaks.keys()) != set(range(1, n + 1)):
            continue
        # Ordered values for peaks 1..N.  Try to build a homogeneous numeric
        # array; fall back to object dtype if values are heterogeneous.
        ordered = [peaks[i] for i in range(1, n + 1)]
        try:
            arr = np.asarray(ordered, dtype=float)
        except (TypeError, ValueError):
            arr = np.asarray(ordered, dtype=object)
        out[base] = arr
    return out


def _strip_position_and_prefix(key: str, command: str, position: int) -> str:
    """Return the descriptor key with position suffix and command prefix stripped.

    Rules
    -----
    1. Strip trailing ``_<position>`` suffix.
    2. Strip leading ``COMMANDNAME_`` prefix (case-insensitive) only when the
       key starts with it — keys with a leading qualifier (e.g.
       ``Log10_LS_Prob_1``) keep their qualifier.

    Examples
    --------
    ::

        _strip("LS_Period_1_0",          "LS",  0) -> "Period_1"
        _strip("Log10_LS_Prob_1_0",      "LS",  0) -> "Log10_LS_Prob_1"
        _strip("LS_Periodogram_Value_1_0","LS", 0) -> "Periodogram_Value_1"
        _strip("RMS_1",                  "rms", 1) -> "RMS"
        _strip("Weighted_RMS_1",         "rms", 1) -> "Weighted_RMS"
    """
    # Step 1 — strip trailing _N
    sfx = f"_{position}"
    if key.endswith(sfx):
        key = key[: -len(sfx)]

    # Step 2 — strip leading COMMANDNAME_ (case-insensitive prefix match)
    pfx = command + "_"
    if key.upper().startswith(pfx.upper()):
        key = key[len(pfx):]

    return key


# ---------------------------------------------------------------------------
# CommandStats
# ---------------------------------------------------------------------------

class CommandStats:
    """Stats for a single occurrence of one vartools command.

    Attribute names are the vartools output variable descriptors with the
    position suffix (and command prefix where applicable) stripped.

    Access
    ------
    ::

        cs = result.varobjs.LS[0]
        cs.Period_1              # attribute access
        cs["Period_1"]           # item access
        cs.keys()                # list available keys
        cs.values()              # list available values
        cs.items()               # (key, value) pairs
    """

    def __init__(
        self,
        command: str,
        position: int,
        stripped: Dict[str, object],
    ) -> None:
        object.__setattr__(self, "_command", command)
        object.__setattr__(self, "_position", position)
        object.__setattr__(self, "_stripped", stripped)

    # ------------------------------------------------------------------
    # Access
    # ------------------------------------------------------------------

    def __getattr__(self, name: str) -> object:
        stripped = object.__getattribute__(self, "_stripped")
        if name in stripped:
            return stripped[name]
        cmd = object.__getattribute__(self, "_command")
        pos = object.__getattribute__(self, "_position")
        raise AttributeError(
            f"No output variable '{name}' for command '{cmd}' "
            f"(position {pos}). "
            f"Available: {list(stripped.keys())}"
        )

    def __getitem__(self, key: str) -> object:
        stripped = object.__getattribute__(self, "_stripped")
        if key in stripped:
            return stripped[key]
        raise KeyError(key)

    def __contains__(self, key: str) -> bool:
        return key in object.__getattribute__(self, "_stripped")

    def __len__(self) -> int:
        return len(object.__getattribute__(self, "_stripped"))

    def keys(self) -> list:
        return list(object.__getattribute__(self, "_stripped").keys())

    def values(self):
        return object.__getattribute__(self, "_stripped").values()

    def items(self):
        return object.__getattribute__(self, "_stripped").items()

    def __dir__(self):
        return list(object.__getattribute__(self, "_stripped").keys()) + object.__dir__(self)

    def to_series(self) -> pd.Series:
        """Return a pd.Series of this command's output variables."""
        return pd.Series(object.__getattribute__(self, "_stripped"))

    def __repr__(self) -> str:
        cmd = object.__getattribute__(self, "_command")
        pos = object.__getattribute__(self, "_position")
        keys = list(object.__getattribute__(self, "_stripped").keys())
        return f"CommandStats({cmd!r}, pos={pos}, keys={keys})"


# ---------------------------------------------------------------------------
# CommandStatsList
# ---------------------------------------------------------------------------

class CommandStatsList:
    """A list of CommandStats for one command name.

    When the list contains exactly one entry, attribute and item access are
    forwarded directly to that entry::

        result.varobjs.LS.Period_1     # shorthand for result.varobjs.LS[0].Period_1
        result.varobjs.LS[0].Period_1  # always valid

    When the list has more than one entry, bare attribute access raises a clear
    error pointing to the indexed form.
    """

    def __init__(self, command: str, entries: List[CommandStats]) -> None:
        object.__setattr__(self, "_command", command)
        object.__setattr__(self, "_entries", list(entries))

    # ------------------------------------------------------------------
    # List-like interface
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return len(object.__getattribute__(self, "_entries"))

    def __iter__(self) -> Iterator[CommandStats]:
        return iter(object.__getattribute__(self, "_entries"))

    def __getitem__(self, idx: int) -> CommandStats:
        return object.__getattribute__(self, "_entries")[idx]

    # ------------------------------------------------------------------
    # Shorthand access (only when len == 1)
    # ------------------------------------------------------------------

    def __getattr__(self, name: str) -> object:
        entries = object.__getattribute__(self, "_entries")
        command = object.__getattribute__(self, "_command")
        if len(entries) == 1:
            return getattr(entries[0], name)
        raise AttributeError(
            f"result.varobjs.{command} has {len(entries)} entries — "
            f"use result.varobjs.{command}[0].{name}, "
            f"result.varobjs.{command}[1].{name}, etc."
        )

    def __dir__(self):
        entries = object.__getattribute__(self, "_entries")
        if len(entries) == 1:
            return dir(entries[0])
        return object.__dir__(self)

    def __repr__(self) -> str:
        cmd = object.__getattribute__(self, "_command")
        entries = object.__getattribute__(self, "_entries")
        if len(entries) == 1:
            return repr(entries[0])
        return f"CommandStatsList({cmd!r}, {entries!r})"


# ---------------------------------------------------------------------------
# VarsNamespace
# ---------------------------------------------------------------------------

class VarsNamespace:
    """Per-command structured access to vartools output variables.

    Built from a flat ``pd.Series`` via ``VarsNamespace.from_series()``.

    Usage::

        result.varobjs.LS.Period_1          # single LS call
        result.varobjs.LS[0].Period_1       # explicit index
        result.varobjs.rms.Weighted_RMS     # rms command outputs
    """

    def __init__(self, by_command: Dict[str, CommandStatsList]) -> None:
        object.__setattr__(self, "_by_command", by_command)

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    @classmethod
    def from_series(
        cls,
        series: pd.Series,
        known_commands: Optional[List[str]] = None,
    ) -> "VarsNamespace":
        """Build from a flat pd.Series of vartools output variables.

        Parameters
        ----------
        series : pd.Series
            The ``result.vars`` Series.
        known_commands : list of str, optional
            The ordered list of command ``_vt_name`` values that produced this
            result (e.g. ``["LS", "rms"]``).  Position N in this list maps to
            keys whose trailing suffix is ``_N``.  When omitted, command names
            are inferred from key prefixes (best-effort).
        """
        # Step 1 — group keys by trailing position integer
        groups: Dict[int, List[tuple]] = {}
        for key in series.index:
            if key == "Name":
                continue
            pos = _parse_position(key)
            if pos is None:
                continue
            groups.setdefault(pos, []).append((key, series[key]))

        # Step 2 — for each position, determine the command name
        by_command: Dict[str, List[CommandStats]] = {}
        for pos, items in sorted(groups.items()):
            if known_commands is not None and pos < len(known_commands):
                cmd = known_commands[pos]
            else:
                cmd = _infer_command_from_items(items)

            stripped = {
                _strip_position_and_prefix(k, cmd, pos): v
                for k, v in items
            }
            stripped = _add_peak_arrays(stripped)
            cs = CommandStats(command=cmd, position=pos, stripped=stripped)
            by_command.setdefault(cmd, []).append(cs)

        return cls({
            cmd: CommandStatsList(cmd, entries)
            for cmd, entries in by_command.items()
        })

    # ------------------------------------------------------------------
    # Access
    # ------------------------------------------------------------------

    def __getattr__(self, name: str) -> CommandStatsList:
        by_command = object.__getattribute__(self, "_by_command")
        if name in by_command:
            return by_command[name]
        raise AttributeError(
            f"No output variables for command '{name}'. "
            f"Available: {list(by_command.keys())}"
        )

    def __contains__(self, name: str) -> bool:
        return name in object.__getattribute__(self, "_by_command")

    def __dir__(self):
        by_command = object.__getattribute__(self, "_by_command")
        return list(by_command.keys()) + object.__dir__(self)

    def commands(self) -> List[str]:
        """Return the command names that have output variables."""
        return list(object.__getattribute__(self, "_by_command").keys())

    def __repr__(self) -> str:
        by_command = object.__getattribute__(self, "_by_command")
        return f"VarsNamespace(commands={list(by_command.keys())})"


# ---------------------------------------------------------------------------
# Fallback command name inference (used when known_commands is not supplied)
# ---------------------------------------------------------------------------

def _infer_command_from_items(items: List[tuple]) -> str:
    """Infer a command name from a group of (key, value) pairs at one position.

    Heuristic: the command name is the longest common prefix of all keys in
    the group (up to the first ``_`` boundary) that appears to be a valid
    vartools command token.  Falls back to the prefix of the first key.
    """
    if not items:
        return "unknown"

    keys = [k for k, _ in items]

    # Collect all leading tokens (up to the first '_descriptor' boundary)
    # from each key by stripping the trailing _N and splitting.
    candidates: List[str] = []
    for key in keys:
        # The key ends with _N; strip that first.
        m = _TRAILING_INT.match(key)
        if not m:
            continue
        body = m.group(1)  # e.g. "LS_Period_1" or "Log10_LS_Prob_1"
        parts = body.split("_")
        # The command is the first part that looks like an identifier
        # (not a number, not a purely lowercase multi-word descriptor).
        # Take just the first part as the simplest heuristic.
        candidates.append(parts[0])

    if not candidates:
        return "unknown"

    # Return the most common candidate
    from collections import Counter
    return Counter(candidates).most_common(1)[0][0]
