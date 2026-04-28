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
import re
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Union, TYPE_CHECKING

import numpy as np
import pandas as pd

# Astropy is optional at import time so that the rest of the package works
# even without it installed (though pyproject.toml lists it as a dependency).
try:
    from astropy.timeseries import TimeSeries
    _HAVE_ASTROPY = True
except (ImportError, Exception):
    _HAVE_ASTROPY = False

if TYPE_CHECKING:
    from astropy.io.fits import Header as _FitsHeader

_STANDARD_COLS = ("t", "mag", "err")


def _read_ascii_header_names(path: Path):
    """Return column names from the last ``# name1 name2 ...`` comment line
    at the top of *path*, or ``None`` if no such line is present.

    Only the contiguous block of leading comment lines is scanned, so a
    ``#commandline`` log line (or any other `#...` line without space-
    separated tokens that look like identifiers) is ignored unless it
    happens to be the last header line and parses as a name list.
    """
    try:
        with open(path, "r") as fh:
            last_tokens = None
            for raw in fh:
                s = raw.lstrip()
                if not s:
                    continue
                if not s.startswith("#"):
                    break
                stripped = s[1:].strip()
                if not stripped:
                    continue
                tokens = stripped.split()
                if all(tok.replace("_", "").replace(".", "").isalnum()
                       and not tok[0].isdigit()
                       for tok in tokens):
                    last_tokens = tokens
            return last_tokens
    except OSError:
        return None

# -----------------------------------------------------------------------------
# FITS-header helpers — shared between read and write.
# -----------------------------------------------------------------------------

# Structural FITS keywords that must be (re)derived from the current
# DataFrame at write time.  These are filtered out on both read (when
# capturing a header onto LightCurve.fitsheader) and write (when emitting
# a preserved header back into an output file), so the user-visible
# ``fitsheader`` only carries observational / provenance metadata.
_STRUCTURAL_FIXED = frozenset({
    "SIMPLE", "BITPIX", "EXTEND", "XTENSION", "END",
    "NAXIS", "PCOUNT", "GCOUNT", "TFIELDS",
    "EXTNAME", "EXTVER", "EXTLEVEL",
})
# Per-column and per-axis indexed keywords: TTYPE1, TFORM2, NAXIS1, …
_STRUCTURAL_INDEXED_RE = re.compile(
    r"^(NAXIS|TTYPE|TFORM|TDISP|TSCAL|TZERO|TNULL|TBCOL|TUNIT|"
    r"TCTYP|TCRPX|TCRVL|TCDLT|TCUNI|TCROT|TDIM)\d+$"
)


def _is_structural_fits_key(keyword: str) -> bool:
    """True for FITS header keywords that must be redetermined from the data."""
    k = keyword.upper()
    return k in _STRUCTURAL_FIXED or bool(_STRUCTURAL_INDEXED_RE.match(k))


def _filter_structural(header: "_FitsHeader") -> "_FitsHeader":
    """Return a copy of *header* with structural keywords stripped."""
    from astropy.io.fits import Header
    out = Header()
    for card in header.cards:
        # Allow blank / COMMENT / HISTORY cards through unchanged; they carry
        # no structural meaning.
        kw = (card.keyword or "").strip()
        if not kw or kw in ("COMMENT", "HISTORY"):
            out.append(card, end=True)
            continue
        if not _is_structural_fits_key(kw):
            out.append(card, end=True)
    return out


def _coerce_fitsheader(value):
    """Convert *value* into an ``astropy.io.fits.Header`` (or return None).

    Accepts ``None`` (returns ``None``), an existing ``Header`` (returned as
    a copy), or any dict-like whose items can be used as FITS keyword/value
    pairs.  Raises ``TypeError`` for anything else.
    """
    if value is None:
        return None
    try:
        from astropy.io import fits
    except ImportError as exc:
        raise ImportError(
            "astropy is required to attach a FITS header to a LightCurve."
        ) from exc
    if isinstance(value, fits.Header):
        return value.copy()
    if isinstance(value, dict):
        return fits.Header(list(value.items()))
    raise TypeError(
        f"fitsheader must be None, an astropy.io.fits.Header, or a dict; "
        f"got {type(value).__name__}"
    )


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

    def __init__(
        self,
        data: pd.DataFrame,
        name: str = "",
        scalars: Optional[Dict[str, float]] = None,
        fitsheader: Optional[Any] = None,
    ) -> None:
        self._df = data.copy()
        self.name = name
        # Per-star scalar variables (VARTOOLS_VECTORTYPE_SCALAR /
        # PERSTARDATA / INLIST).  Canonical names, no ``_N`` suffix.
        # Populated by pyvartools during capture; also carried across chained
        # command invocations so subsequent commands can reference prior
        # results in analytic expressions.
        self.scalars: Dict[str, float] = dict(scalars or {})
        # Optional FITS header preserved from the input file (merged primary +
        # extension, structural keywords filtered).  Emitted back onto the
        # primary HDU of FITS output via ``to_file``.  ``None`` if the light
        # curve did not come from a FITS file.  Accepts an ``astropy.io.fits``
        # ``Header`` instance or ``None``; other truthy values are converted
        # to a ``Header`` when astropy is available.
        self.fitsheader = _coerce_fitsheader(fitsheader)

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
        scalars: Optional[Dict[str, float]] = None,
        fitsheader: Optional[Any] = None,
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
        fitsheader : astropy.io.fits.Header or dict, optional
            FITS header to round-trip through :meth:`to_file`.  Normally set
            automatically by :meth:`from_file` when reading a FITS file.
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
        return cls(pd.DataFrame(d), name=name, scalars=scalars,
                   fitsheader=fitsheader)

    @classmethod
    def from_dataframe(
        cls,
        df: pd.DataFrame,
        name: str = "",
        scalars: Optional[Dict[str, float]] = None,
        fitsheader: Optional[Any] = None,
    ) -> "LightCurve":
        """Construct from a pandas DataFrame.

        Any DataFrame is accepted.  Columns named ``t``, ``mag``, and ``err``
        are treated as the standard vartools vectors when present; others are
        preserved as auxiliary columns.
        """
        return cls(df, name=name, scalars=scalars, fitsheader=fitsheader)

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
    def from_files(
        cls,
        paths: Sequence[Union[str, Path]],
        name: str = "",
        lcnum_col: str = "lcnum",
        sort: bool = True,
        **read_kwargs: Any,
    ) -> "LightCurve":
        """Read several light-curve files and combine them into one ``LightCurve``.

        Each file is loaded via :meth:`from_file` (so ``read_kwargs`` such as
        ``format``, ``t_col``, ``mag_col``, ``err_col``, ``hdu`` may be
        forwarded to the FITS reader).  The resulting DataFrames are
        concatenated, an integer ``lcnum_col`` is filled in (0 for the first
        file, 1 for the second, …), and the combined frame is optionally
        time-sorted.

        This mirrors what vartools' ``-l ... combinelcs`` mode does internally,
        but does the merge entirely in Python — useful when you want to feed a
        single combined LC to commands that don't go through the file-list
        machinery (e.g. ``Pipeline.run(lc)`` after calling :func:`from_files`).

        Parameters
        ----------
        paths : sequence of str | Path
            Files to combine, in the order their points should be tagged.
        name : str, optional
            Label for the combined light curve.  Defaults to the stem of the
            first file.
        lcnum_col : str
            Name of the integer column added to identify the source file of
            each row.  Default ``"lcnum"``.  If the column already exists in
            any of the input frames, it is overwritten.
        sort : bool
            If True (default), sort the combined frame by ``t`` ascending.
        **read_kwargs
            Extra keyword arguments forwarded to :meth:`from_file` (e.g.
            ``t_col``, ``mag_col``, ``err_col``, ``hdu``).

        Returns
        -------
        LightCurve
        """
        paths = list(paths)
        if not paths:
            raise ValueError("from_files() requires at least one path.")
        frames = []
        for i, p in enumerate(paths):
            lc = cls.from_file(p, **read_kwargs)
            df = lc._df.copy()
            df[lcnum_col] = i
            frames.append(df)
        merged = pd.concat(frames, ignore_index=True)
        if sort and "t" in merged.columns:
            merged = merged.sort_values("t").reset_index(drop=True)
        if not name:
            name = Path(paths[0]).stem
        return cls(merged, name=name)

    @classmethod
    def _from_ascii(cls, path: Path, name: str) -> "LightCurve":
        df = pd.read_csv(path, sep=r"\s+", comment="#", header=None)
        ncols = df.shape[1]
        header_names = _read_ascii_header_names(path)
        if header_names is not None and len(header_names) == ncols:
            col_names = header_names
        elif ncols >= 3:
            col_names = list(_STANDARD_COLS) + [f"col{i+4}" for i in range(ncols - 3)]
        else:
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
            # Merge primary + data-HDU headers into a single observational
            # header, filtering out column/axis-structural keywords.  Primary
            # first so extension-HDU keys can override on conflict.
            merged = fits.Header(hdul[0].header.copy())
            if hdu != 0:
                merged.update(hdul[hdu].header)
            header = _filter_structural(merged)
        return cls.from_arrays(t, mag, err, aux=aux or None, name=name,
                               fitsheader=header)

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
        if self.t is None:
            raise ValueError(
                "LightCurve has no 't' column; cannot construct a TimeSeries. "
                "Supply a time array via from_arrays(t=...) or ensure the "
                "source file has a time column."
            )
        from astropy.time import Time
        ts = TimeSeries(time=Time(self.t, format="jd"))
        for col in self._df.columns:
            if col != "t":
                ts[col] = self._df[col].to_numpy()
        return ts

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    def to_file(
        self,
        path: Union[str, Path],
        format: Optional[str] = None,
    ) -> None:
        """Write the light curve to a file.

        Parameters
        ----------
        path : str | Path
        format : str, optional
            ``"ascii"`` (default for most extensions) or ``"fits"``.
            Auto-detected from the file extension when omitted.

        Notes
        -----
        ASCII output is whitespace-separated with 10 decimal places of
        precision and no header line, matching the vartools default format.
        FITS output requires astropy.
        """
        path = Path(path)
        if format is None:
            fmt = path.suffix.lower()
            if fmt in (".fits", ".fit", ".fts"):
                format = "fits"
            else:
                format = "ascii"

        if format == "fits":
            try:
                from astropy.io import fits
                from astropy.table import Table
            except ImportError as e:
                raise ImportError("astropy is required to write FITS files.") from e
            tbl = Table.from_pandas(self._df)
            bin_hdu = fits.BinTableHDU(data=tbl)
            primary = fits.PrimaryHDU()
            # When a preserved header is available, inject its non-structural
            # keywords onto the primary HDU.  Skip when fitsheader is None so
            # the default write path is byte-for-byte unchanged from before.
            if self.fitsheader is not None:
                for card in _filter_structural(self.fitsheader).cards:
                    primary.header.append(card, end=True)
            fits.HDUList([primary, bin_hdu]).writeto(str(path), overwrite=True)
        else:
            self._df.to_csv(path, sep=" ", header=False, index=False,
                            float_format="%.10f")

    def _to_tempfile(self, dir: Optional[str] = None) -> str:
        """Write the light curve to a named temp file and return its path.

        Internal helper — callers are responsible for deleting the file.
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
    # Visualisation
    # ------------------------------------------------------------------

    def plot(self, ax=None, **kwargs):
        """Quick-look plot of the light curve.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on.  A new figure and axes are created if omitted.
        **kwargs
            Passed to ``ax.errorbar`` (or ``ax.plot`` when there is no error
            column).  Override defaults such as ``fmt``, ``markersize``, etc.

        Returns
        -------
        matplotlib.axes.Axes

        Notes
        -----
        The y-axis is inverted automatically (standard astronomical magnitude
        convention).  Requires matplotlib.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError as e:
            raise ImportError("matplotlib is required for LightCurve.plot().") from e

        if ax is None:
            _, ax = plt.subplots()

        kw = dict(fmt=".", markersize=3, elinewidth=0.5, capsize=0)
        kw.update(kwargs)

        t = self.t
        mag = self.mag
        err = self.err

        if t is None or mag is None:
            raise ValueError(
                "LightCurve has no 't' or 'mag' column; cannot plot."
            )

        if err is not None:
            ax.errorbar(t, mag, err, **kw)
        else:
            # Strip errorbar-specific keys that aren't valid for plot()
            plot_kw = {k: v for k, v in kw.items()
                       if k not in ("elinewidth", "capsize", "ecolor", "barsabove")}
            fmt = plot_kw.pop("fmt", ".")
            ax.plot(t, mag, fmt, **plot_kw)

        ax.invert_yaxis()
        ax.set_xlabel("Time")
        ax.set_ylabel("Magnitude")
        if self.name:
            ax.set_title(self.name)

        return ax

    # ------------------------------------------------------------------
    # Dunder
    # ------------------------------------------------------------------

    @property
    def shape(self) -> tuple:
        """(n_points, n_columns) — mirrors ``DataFrame.shape``."""
        return self._df.shape

    def __len__(self) -> int:
        return len(self._df)

    def __repr__(self) -> str:
        extra = f", scalars={len(self.scalars)}" if self.scalars else ""
        return (f"LightCurve(name={self.name!r}, n={len(self)}, "
                f"cols={list(self._df.columns)}{extra})")
