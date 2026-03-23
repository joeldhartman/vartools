"""
Result objects returned by Pipeline.run() and Pipeline.run_batch().
"""

from __future__ import annotations

import re
from typing import Dict, Optional, Union

import pandas as pd

from .lightcurve import LightCurve


class RunError(RuntimeError):
    """Raised when vartools exits with a non-zero status or writes to stderr."""
    pass


class Result:
    """Output from a single-light-curve pipeline run.

    Attributes
    ----------
    stats : pd.Series
        Series indexed by vartools output variable names
        (e.g. 'LS_Period_1_0', 'Log10_LS_Prob_1_0', ...) plus 'Name'.
        Values are scalars, so ``result.stats["LS_Period_1_0"]`` works directly.
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
        stats: Union[pd.Series, pd.DataFrame],
        lc: Optional[LightCurve] = None,
        files: Optional[Dict[str, pd.DataFrame]] = None,
    ) -> None:
        self.stats = stats
        self.lc = lc
        self.files = files or {}

    def __repr__(self) -> str:
        keys = list(self.stats.index) if isinstance(self.stats, pd.Series) else list(self.stats.columns)
        return (f"Result(stats_keys={keys}, "
                f"lc={'yes' if self.lc else 'no'}, "
                f"files={list(self.files.keys())})")


class BatchResult:
    """Output from a multi-light-curve pipeline run.

    Attributes
    ----------
    stats : pd.DataFrame
        One row per light curve; columns are vartools output variable names.
    lcs : list of LightCurve or None
        The (possibly modified) output light curves, in input order.
        ``None`` if ``capture_lc=False`` was passed to the pipeline run method.
    files : dict
        Mapping of logical name → list of DataFrames (one per light curve)
        for any output files that were requested.
    error : RunError or None
        If the vartools process failed, the exception is stored here and
        ``stats`` will be empty.  For successful runs this is ``None``.
    """

    def __init__(
        self,
        stats: pd.DataFrame,
        lcs: Optional[list] = None,
        files: Optional[Dict[str, list]] = None,
        error: Optional["RunError"] = None,
    ) -> None:
        self.stats = stats
        self.lcs = lcs
        self.files = files or {}
        self.error = error

    @property
    def ok(self) -> bool:
        """True if the run completed without error."""
        return self.error is None

    def __repr__(self) -> str:
        return (f"BatchResult(n={len(self.stats)}, "
                f"cols={list(self.stats.columns)}, "
                f"lcs={'yes' if self.lcs else 'no'}, "
                f"ok={self.ok}, "
                f"files={list(self.files.keys())})")


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

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
