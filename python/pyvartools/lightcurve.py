"""
LightCurve class: in-memory representation of a vartools light curve.

Internally stored as a pandas DataFrame.  The first three columns are
always 't' (time), 'mag' (magnitude or flux), and 'err' (uncertainty).
Additional columns may be present for auxiliary variables.

Accepts construction from:
  - pandas DataFrame
  - numpy arrays (t, mag, err, and optional dict of aux arrays)
  - astropy TimeSeries
  - a file path (vartools three-column ASCII format)
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Dict, Optional, Union

import numpy as np
import pandas as pd

# Astropy is optional at import time so that the rest of the package works
# even without it installed (though pyproject.toml lists it as a dependency).
try:
    from astropy.timeseries import TimeSeries
    import astropy.units as u
    _HAVE_ASTROPY = True
except (ImportError, Exception):
    _HAVE_ASTROPY = False

_STANDARD_COLS = ("t", "mag", "err")


class LightCurve:
    """In-memory representation of a vartools light curve.

    Parameters
    ----------
    data : pd.DataFrame
        Any DataFrame whose columns form a valid vartools input.  The columns
        ``t``, ``mag``, and ``err`` are treated as the standard time,
        magnitude, and uncertainty vectors when present, but all three are
        optional.  If some or all are absent the caller is responsible for
        telling vartools how to interpret the columns (via the ``columns``
        parameter on the Pipeline run methods, or by using
        ``-inputlcformat`` on the CLI directly).  Additional columns beyond
        the standard three are preserved as auxiliary vartools variables.
    name : str, optional
        A label for this light curve (used as the vartools 'Name' field).
    """

    def __init__(self, data: pd.DataFrame, name: str = "") -> None:
        self._df = data.copy()
        self.name = name

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_arrays(
        cls,
        t: Optional[np.ndarray] = None,
        mag: Optional[np.ndarray] = None,
        err: Optional[np.ndarray] = None,
        aux: Optional[Dict[str, np.ndarray]] = None,
        name: str = "",
    ) -> "LightCurve":
        """Construct from numpy arrays.

        Parameters
        ----------
        t, mag, err : array-like or None
            Time, magnitude/flux, and uncertainty arrays.  Any or all may be
            ``None`` when vartools will generate them internally (e.g. when
            using ``-inputlcformat`` without these columns on the CLI).
        aux : dict, optional
            Mapping of column name → array for additional columns.
        name : str, optional
            Light curve label.
        """
        d = {}
        if t is not None:
            d["t"] = np.asarray(t)
        if mag is not None:
            d["mag"] = np.asarray(mag)
        if err is not None:
            d["err"] = np.asarray(err)
        if aux:
            for k, v in aux.items():
                d[k] = np.asarray(v)
        return cls(pd.DataFrame(d), name=name)

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame, name: str = "") -> "LightCurve":
        """Construct from a pandas DataFrame.

        Any DataFrame is accepted.  Columns named ``t``, ``mag``, and ``err``
        are treated as the standard vartools vectors when present; others are
        preserved as auxiliary columns.
        """
        return cls(df, name=name)

    @classmethod
    def from_timeseries(cls, ts: "TimeSeries", mag_col: str = "mag",
                        err_col: str = "err", name: str = "") -> "LightCurve":
        """Construct from an astropy TimeSeries.

        Parameters
        ----------
        ts : astropy.timeseries.TimeSeries
        mag_col : str
            Column name in *ts* to use as 'mag'.
        err_col : str
            Column name in *ts* to use as 'err'.
        name : str, optional
        """
        if not _HAVE_ASTROPY:
            raise ImportError("astropy is required to use from_timeseries().")
        t = ts.time.jd
        mag = np.asarray(ts[mag_col])
        err = np.asarray(ts[err_col])
        aux_cols = [c for c in ts.colnames if c not in ("time", mag_col, err_col)]
        aux = {c: np.asarray(ts[c]) for c in aux_cols}
        return cls.from_arrays(t, mag, err, aux=aux or None, name=name)

    @classmethod
    def from_file(
        cls,
        path: Union[str, Path],
        name: str = "",
        format: Optional[str] = None,
        t_col: str = "BJD",
        mag_col: str = "Mag",
        err_col: str = "Err",
        hdu: int = 1,
    ) -> "LightCurve":
        """Read a light curve from a file.

        Format is detected automatically from the file extension unless
        ``format`` is specified explicitly.

        Supported formats
        -----------------
        ``"ascii"`` (default for files without a recognised extension)
            Whitespace-separated columns: time  mag  err  [aux ...].
            Lines beginning with ``#`` are treated as comments.
            Columns are named ``t``, ``mag``, ``err``, ``col4``, ...

        ``"fits"`` (detected for ``.fits``, ``.fit``, ``.fts``)
            Reads a binary or ASCII table from the specified HDU.
            The columns named by *t_col*, *mag_col*, and *err_col* are mapped
            to ``t``, ``mag``, and ``err``.  All other table columns are kept
            as auxiliary columns.

        Parameters
        ----------
        path : str | Path
        name : str, optional
            Label for the light curve.  Defaults to the file stem.
        format : str, optional
            ``"ascii"`` or ``"fits"``.  Auto-detected from extension if omitted.
        t_col, mag_col, err_col : str
            Column names to use as time, magnitude, and error when reading FITS.
        hdu : int
            HDU index to read from when reading FITS (default 1).
        """
        path = Path(path)
        lc_name = name or path.stem

        if format is None:
            fmt = path.suffix.lower()
            if fmt in (".fits", ".fit", ".fts"):
                format = "fits"
            else:
                format = "ascii"

        if format == "fits":
            return cls._from_fits(path, lc_name, t_col, mag_col, err_col, hdu)
        else:
            return cls._from_ascii(path, lc_name)

    @classmethod
    def _from_ascii(cls, path: Path, name: str) -> "LightCurve":
        df = pd.read_csv(path, sep=r"\s+", comment="#", header=None)
        ncols = df.shape[1]
        if ncols >= 3:
            # Standard convention: first three columns are t, mag, err.
            col_names = list(_STANDARD_COLS) + [f"col{i+4}" for i in range(ncols - 3)]
        else:
            # Fewer than 3 columns — caller will need to supply column
            # semantics via the Pipeline `columns` parameter.
            col_names = [f"col{i+1}" for i in range(ncols)]
        df.columns = col_names
        return cls(df, name=name)

    @classmethod
    def _from_fits(
        cls,
        path: Path,
        name: str,
        t_col: str,
        mag_col: str,
        err_col: str,
        hdu: int,
    ) -> "LightCurve":
        try:
            from astropy.io import fits
        except ImportError as e:
            raise ImportError(
                "astropy is required to read FITS files."
            ) from e
        with fits.open(path) as hdul:
            table = hdul[hdu]
            cols = {c.name.upper(): c.name for c in table.columns}
            def _get(wanted: str) -> np.ndarray:
                key = cols.get(wanted.upper())
                if key is None:
                    raise ValueError(
                        f"Column '{wanted}' not found in FITS HDU {hdu} of "
                        f"'{path}'.  Available: {list(cols.values())}"
                    )
                return np.asarray(table.data[key], dtype=float)

            t   = _get(t_col)
            mag = _get(mag_col)
            err = _get(err_col)
            skip = {t_col.upper(), mag_col.upper(), err_col.upper()}
            aux = {
                c.name: np.asarray(table.data[c.name])
                for c in table.columns
                if c.name.upper() not in skip
            }
        return cls.from_arrays(t, mag, err, aux=aux or None, name=name)

    # ------------------------------------------------------------------
    # Conversion helpers
    # ------------------------------------------------------------------

    @property
    def t(self) -> Optional[np.ndarray]:
        return self._df["t"].to_numpy() if "t" in self._df.columns else None

    @property
    def mag(self) -> Optional[np.ndarray]:
        return self._df["mag"].to_numpy() if "mag" in self._df.columns else None

    @property
    def err(self) -> Optional[np.ndarray]:
        return self._df["err"].to_numpy() if "err" in self._df.columns else None

    def to_dataframe(self) -> pd.DataFrame:
        """Return a copy of the internal DataFrame."""
        return self._df.copy()

    def to_arrays(self) -> tuple:
        """Return (t, mag, err) as numpy arrays."""
        return self.t, self.mag, self.err

    def to_timeseries(self) -> "TimeSeries":
        """Convert to an astropy TimeSeries."""
        if not _HAVE_ASTROPY:
            raise ImportError("astropy is required to use to_timeseries().")
        from astropy.time import Time
        ts = TimeSeries(time=Time(self.t, format="jd"))
        for col in self._df.columns:
            if col != "t":
                ts[col] = self._df[col].to_numpy()
        return ts

    # ------------------------------------------------------------------
    # Serialisation (for passing to the vartools binary via temp files)
    # ------------------------------------------------------------------

    def to_tempfile(self, dir: Optional[str] = None) -> str:
        """Write the light curve to a named temp file and return its path.

        The caller is responsible for deleting the file when done.
        """
        fd, path = tempfile.mkstemp(suffix=".lc", dir=dir)
        try:
            with os.fdopen(fd, "w") as f:
                self._df.to_csv(f, sep=" ", header=False, index=False,
                                float_format="%.10f")
        except Exception:
            os.unlink(path)
            raise
        return path

    # ------------------------------------------------------------------
    # Dunder
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return len(self._df)

    def __repr__(self) -> str:
        return (f"LightCurve(name={self.name!r}, n={len(self)}, "
                f"cols={list(self._df.columns)})")
