"""
Top-level command functions for pyvartools.

Exposes every vartools command as a module-level callable that accepts a light
curve in any of several convenient forms as its first argument::

    import pyvartools as vt

    result = vt.LS("EXAMPLES/2", 1.0, 2.0, 0.01)       # filename / Path
    result = vt.LS(df, 1.0, 2.0, 0.01)                 # pandas DataFrame
    result = vt.LS(lc, 1.0, 2.0, 0.01)                 # LightCurve object
    result = vt.LS(arr2d, 1.0, 2.0, 0.01)              # 2-D numpy array
    result = vt.LS((t, mag, err), 1.0, 2.0, 0.01)      # tuple of 1-D arrays

Each function coerces the first argument to a ``LightCurve``, then runs
``Pipeline([cmd]).run(lc, ...)`` and returns a ``Result``.

Run-time options (``capture_lc``, ``timeout``, ``randseed``, ``skipmissing``,
``jdtol``, ``matchstringid``) may be passed as keyword arguments alongside
the command parameters; they are forwarded to ``Pipeline.run()``.

``capture_lc`` defaults to ``True`` for all top-level functions.

Note: for FITS files that require non-default column names (``t_col``,
``mag_col``, ``err_col``, ``hdu``), construct a ``LightCurve`` object
explicitly with ``LightCurve.from_file(path, ...)`` and pass it in.
"""

from __future__ import annotations

import inspect
from pathlib import Path
from typing import List, TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from .lightcurve import LightCurve

# Standard column names used when interpreting a 2-D numpy array.
_STANDARD_COLS = ("t", "mag", "err")

# Command names that are pass-through / extension mechanisms, not standalone commands.
_SKIP = frozenset({"Raw", "UserCommand"})


# ---------------------------------------------------------------------------
# Light-curve input coercion
# ---------------------------------------------------------------------------

def _coerce_lc_input(obj) -> "LightCurve":
    """Convert any light-curve-like value to a ``LightCurve``.

    Accepted forms
    --------------
    LightCurve
        Returned as-is.
    str or Path
        Read via ``LightCurve.from_file()``.  Uses auto-detected format
        (ASCII or FITS by file extension) and default column names.
    pd.DataFrame
        Wrapped via ``LightCurve.from_dataframe()``.
    2-D np.ndarray
        Columns are mapped to ``t``, ``mag``, ``err``, ``col4``, … in order.
    tuple or list of 1-D arrays
        Unpacked as ``(t, mag, err)``.  Any trailing elements beyond the
        first three are ignored.
    astropy TimeSeries
        Converted via ``LightCurve.from_timeseries()``.
    """
    from .lightcurve import LightCurve

    if isinstance(obj, LightCurve):
        return obj

    if isinstance(obj, (str, Path)):
        return LightCurve.from_file(obj)

    if isinstance(obj, pd.DataFrame):
        return LightCurve.from_dataframe(obj)

    if isinstance(obj, np.ndarray):
        if obj.ndim != 2:
            raise TypeError(
                f"numpy array must be 2-D (got shape {obj.shape}). "
                f"To pass individual arrays, use a tuple: (t, mag, err)."
            )
        ncols = obj.shape[1]
        col_names = (
            list(_STANDARD_COLS[:min(ncols, 3)])
            + [f"col{i + 4}" for i in range(ncols - 3)]
        )
        return LightCurve.from_dataframe(pd.DataFrame(obj, columns=col_names))

    if isinstance(obj, (tuple, list)):
        arrs = list(obj)
        t   = np.asarray(arrs[0]) if len(arrs) > 0 else None
        mag = np.asarray(arrs[1]) if len(arrs) > 1 else None
        err = np.asarray(arrs[2]) if len(arrs) > 2 else None
        return LightCurve.from_arrays(t=t, mag=mag, err=err)

    # astropy TimeSeries (optional dependency)
    try:
        from astropy.timeseries import TimeSeries
        if isinstance(obj, TimeSeries):
            return LightCurve.from_timeseries(obj)
    except ImportError:
        pass

    raise TypeError(
        f"Cannot convert {type(obj).__name__!r} to LightCurve. "
        f"Pass a LightCurve, str/Path, pd.DataFrame, 2-D np.ndarray, "
        f"or a tuple of (t, mag, err) 1-D arrays."
    )


# ---------------------------------------------------------------------------
# Function factory
# ---------------------------------------------------------------------------

def _make_toplevel(cmd_cls):
    """Return a top-level function that runs *cmd_cls* on a flexible LC input."""
    vt_name = cmd_cls._vt_name

    def _fn(lc_input, *args, **kwargs):
        from ._method_gen import _split_kwargs
        run_opts, cmd_kwargs = _split_kwargs(kwargs)
        capture_lc = run_opts.pop("capture_lc", True)
        lc = _coerce_lc_input(lc_input)
        cmd = cmd_cls(*args, **cmd_kwargs)
        from .pipeline import Pipeline
        return Pipeline([cmd]).run(lc, capture_lc=capture_lc, **run_opts)

    _fn.__name__ = vt_name
    _fn.__qualname__ = f"pyvartools.{vt_name}"
    _fn.__doc__ = (
        f"Run -{vt_name} on a light curve and return a Result.\n\n"
        f"Parameters\n"
        f"----------\n"
        f"lc_input : LightCurve | str | Path | pd.DataFrame | 2-D np.ndarray | tuple\n"
        f"    The light curve to analyse.  Strings and Path objects are read via\n"
        f"    ``LightCurve.from_file()`` (auto-detects ASCII or FITS format).\n"
        f"    A 2-D numpy array is interpreted as columns (t, mag, err, ...).\n"
        f"    A tuple or list is unpacked as (t, mag, err) 1-D arrays.\n"
        f"    For FITS files with non-default column names, construct a\n"
        f"    ``LightCurve`` explicitly and pass it in.\n"
        f"*args, **kwargs\n"
        f"    Command parameters forwarded to ``commands.{vt_name}.__init__``.\n"
        f"    Run-time options ``capture_lc`` (default True), ``timeout``,\n"
        f"    ``randseed``, ``skipmissing``, ``jdtol``, ``matchstringid`` may\n"
        f"    also be passed as keyword arguments.\n\n"
        + (inspect.cleandoc(cmd_cls.__doc__) if cmd_cls.__doc__ else "")
    )
    return _fn


# ---------------------------------------------------------------------------
# Module attachment
# ---------------------------------------------------------------------------

def _attach_toplevel_commands(module) -> List[str]:
    """Generate and attach one top-level function per command to *module*.

    Called once at import time from ``pyvartools/__init__.py``.

    Returns
    -------
    list of str
        The command names that were added to the module (for ``__all__``).
    """
    from . import commands as _cmds

    added: List[str] = []
    for name in _cmds.__all__:
        if name in _SKIP or name == "VartoolsCommand":
            continue
        cls = getattr(_cmds, name, None)
        if cls is None:
            continue
        if not (isinstance(cls, type) and issubclass(cls, _cmds.VartoolsCommand)):
            continue
        vt_name = getattr(cls, "_vt_name", None)
        if not vt_name:
            continue
        fn = _make_toplevel(cls)
        setattr(module, vt_name, fn)
        added.append(vt_name)

    return added
