"""ctypes binding to libvartoolspipeline — init-once / process-many API.

Usage (low-level)::

    lp = LibPipeline(["-rms", "-oneline"])
    var = lp.process_lc(t, mag, err, name="my_lc")  # pd.Series
    del lp  # calls vartools_free_pipeline
"""

from __future__ import annotations

import ctypes
from typing import List, Optional

import numpy as np
import pandas as pd

from .results import parse_oneline_output

# ---------------------------------------------------------------------------
# Library loading (lazy, module-level singleton)
# ---------------------------------------------------------------------------

_lib: Optional[ctypes.CDLL] = None


def _load_lib() -> ctypes.CDLL:
    """Load and configure libvartoolspipeline (cached after first call)."""
    global _lib
    if _lib is not None:
        return _lib

    from ._binary import find_library
    path = find_library()
    lib = ctypes.CDLL(path)

    lib.vartools_init_pipeline.restype = ctypes.c_void_p
    lib.vartools_init_pipeline.argtypes = [
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_char_p),
    ]

    lib.vartools_process_lc.restype = ctypes.c_int
    lib.vartools_process_lc.argtypes = [
        ctypes.c_void_p,                  # ProgramData *p
        ctypes.POINTER(ctypes.c_double),  # const double *t
        ctypes.POINTER(ctypes.c_double),  # const double *mag
        ctypes.POINTER(ctypes.c_double),  # const double *err
        ctypes.c_int,                     # int n
        ctypes.c_char_p,                  # const char *lc_name
        ctypes.c_void_p,                  # char *outbuf
        ctypes.c_int,                     # int outbuf_size
    ]

    lib.vartools_free_pipeline.restype = None
    lib.vartools_free_pipeline.argtypes = [ctypes.c_void_p]

    # --- Extended API for LC capture and injection ---

    class VartoolsVarInfo(ctypes.Structure):
        _fields_ = [
            ("name", ctypes.c_char_p),
            ("datatype", ctypes.c_int),
            ("vectortype", ctypes.c_int),
            ("dataptr", ctypes.c_void_p),
            ("length", ctypes.c_int),
        ]

    lib._VartoolsVarInfo = VartoolsVarInfo

    lib.vartools_get_lc_variables.restype = ctypes.c_int
    lib.vartools_get_lc_variables.argtypes = [
        ctypes.c_void_p,                              # ProgramData *p
        ctypes.POINTER(ctypes.c_int),                 # int *n_vars
        ctypes.POINTER(VartoolsVarInfo),              # VartoolsVarInfo *vars
        ctypes.c_int,                                 # int max_vars
    ]

    lib.vartools_get_njd.restype = ctypes.c_int
    lib.vartools_get_njd.argtypes = [ctypes.c_void_p]

    lib.vartools_set_lc_data.restype = ctypes.c_int
    lib.vartools_set_lc_data.argtypes = [
        ctypes.c_void_p,                              # ProgramData *p
        ctypes.c_int,                                 # int n_points
        ctypes.c_int,                                 # int n_columns
        ctypes.POINTER(ctypes.c_char_p),              # const char **col_names
        ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),  # const double **col_data
        ctypes.c_char_p,                              # const char *lc_name
    ]

    # --- Step D: in-memory LC capture for cmd.o(capture=True) ---

    lib.vartools_get_captured_lc.restype = ctypes.c_int
    lib.vartools_get_captured_lc.argtypes = [
        ctypes.c_void_p,                              # ProgramData *p
        ctypes.c_char_p,                              # const char *id
        ctypes.POINTER(VartoolsVarInfo),              # VartoolsVarInfo *vars
        ctypes.c_int,                                 # int max_vars
        ctypes.POINTER(ctypes.c_int),                 # int *n_vars
    ]

    lib.vartools_get_captured_njd.restype = ctypes.c_int
    lib.vartools_get_captured_njd.argtypes = [
        ctypes.c_void_p,                              # ProgramData *p
        ctypes.c_char_p,                              # const char *id
    ]

    _lib = lib
    return _lib


# VARTOOLS type/vectortype constants (from programdata.h / analytic.h)
VARTOOLS_TYPE_DOUBLE = 0
VARTOOLS_TYPE_STRING = 1
VARTOOLS_TYPE_INT = 2
VARTOOLS_TYPE_FLOAT = 3
VARTOOLS_TYPE_LONG = 4
VARTOOLS_TYPE_SHORT = 7

VARTOOLS_VECTORTYPE_CONSTANT    = 0
VARTOOLS_VECTORTYPE_SCALAR      = 1
VARTOOLS_VECTORTYPE_INLIST      = 2
VARTOOLS_VECTORTYPE_LC          = 3
VARTOOLS_VECTORTYPE_OUTCOLUMN   = 4
VARTOOLS_VECTORTYPE_PERSTARDATA = 5


# ---------------------------------------------------------------------------
# LibPipeline
# ---------------------------------------------------------------------------

class LibPipeline:
    """In-process vartools pipeline (init-once / process-many).

    Parameters
    ----------
    argv_list : list of str
        vartools command-line arguments *without* the program name, e.g.
        ``["-rms", "-oneline"]``.  "vartools" is prepended automatically.
    """

    def __init__(self, argv_list: List[str]) -> None:
        lib = _load_lib()

        # Prepend the program name as argv[0], as main() would have it.
        full_argv = ["vartools"] + list(argv_list)
        argc = len(full_argv)
        argv_enc = [a.encode() for a in full_argv]
        argv_c = (ctypes.c_char_p * argc)(*argv_enc)

        # Grab the optional diagnostic-message getter if the loaded library
        # is new enough to expose it.
        _getmsg = getattr(lib, "vartools_last_error_message", None)
        if _getmsg is not None:
            _getmsg.restype = ctypes.c_char_p
            _getmsg.argtypes = []

        self._p = lib.vartools_init_pipeline(argc, argv_c)
        if not self._p:
            msg = ""
            if _getmsg is not None:
                raw = _getmsg()
                if raw:
                    msg = raw.decode("utf-8", "replace").strip()
            detail = f": {msg}" if msg else ""
            raise RuntimeError(
                "vartools_init_pipeline failed" + detail
            )

    def process_lc(
        self,
        t: np.ndarray,
        mag: np.ndarray,
        err: np.ndarray,
        name: str = "lc",
    ) -> pd.Series:
        """Run the pipeline on one light curve and return stats as a Series.

        Parameters
        ----------
        t, mag, err : array-like
            Time, magnitude/flux, and uncertainty arrays (length n).
        name : str
            Label written into the output 'Name' field.

        Returns
        -------
        pd.Series indexed by vartools output variable names.
        """
        lib = _load_lib()

        t_arr   = np.ascontiguousarray(t,   dtype=np.float64)
        mag_arr = np.ascontiguousarray(mag, dtype=np.float64)
        err_arr = np.ascontiguousarray(err, dtype=np.float64)
        n = len(t_arr)

        outbuf = ctypes.create_string_buffer(65536)

        rc = lib.vartools_process_lc(
            self._p,
            t_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            mag_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            err_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_int(n),
            name.encode(),
            outbuf,
            ctypes.c_int(65536),
        )

        if rc != 0:
            raise RuntimeError(f"vartools_process_lc returned error code {rc}")

        text = outbuf.value.decode("utf-8", errors="replace")
        df = parse_oneline_output(text)
        if df.empty:
            return pd.Series(dtype=object)
        return df.iloc[0]

    def process_lc_capture(
        self,
        t: np.ndarray,
        mag: np.ndarray,
        err: np.ndarray,
        name: str = "lc",
        extra_columns: Optional[dict] = None,
    ) -> tuple:
        """Run the pipeline and return stats, LC columns, and scalars.

        Parameters
        ----------
        t, mag, err : array-like
            Input light curve arrays.
        name : str
            Label for the output 'Name' field.
        extra_columns : dict, optional
            Additional LC columns to inject, e.g. ``{"myvar": array}``.

        Returns
        -------
        (stats, lc_columns, scalars) where ``stats`` is a pd.Series of
        OUTCOLUMN values, ``lc_columns`` is a dict mapping LC-type
        variable name → numpy array, and ``scalars`` is a dict mapping
        SCALAR / PERSTARDATA / INLIST variable name → scalar value.
        """
        lib = _load_lib()

        if extra_columns:
            self._inject_lc_data(t, mag, err, name, extra_columns)
            # Run commands without re-injecting t/mag/err
            outbuf = ctypes.create_string_buffer(65536)
            n = len(np.ascontiguousarray(t, dtype=np.float64))
            rc = lib.vartools_process_lc(
                self._p,
                np.ascontiguousarray(t, dtype=np.float64).ctypes.data_as(
                    ctypes.POINTER(ctypes.c_double)),
                np.ascontiguousarray(mag, dtype=np.float64).ctypes.data_as(
                    ctypes.POINTER(ctypes.c_double)),
                np.ascontiguousarray(err, dtype=np.float64).ctypes.data_as(
                    ctypes.POINTER(ctypes.c_double)),
                ctypes.c_int(n), name.encode(), outbuf, ctypes.c_int(65536),
            )
        else:
            t_arr = np.ascontiguousarray(t, dtype=np.float64)
            mag_arr = np.ascontiguousarray(mag, dtype=np.float64)
            err_arr = np.ascontiguousarray(err, dtype=np.float64)
            n = len(t_arr)
            outbuf = ctypes.create_string_buffer(65536)
            rc = lib.vartools_process_lc(
                self._p,
                t_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                mag_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                err_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ctypes.c_int(n), name.encode(), outbuf, ctypes.c_int(65536),
            )

        if rc != 0:
            raise RuntimeError(f"vartools_process_lc returned error code {rc}")

        text = outbuf.value.decode("utf-8", errors="replace")
        df = parse_oneline_output(text)
        stats = df.iloc[0] if not df.empty else pd.Series(dtype=object)

        lc_columns, scalars = self._read_lc_variables()
        return stats, lc_columns, scalars

    def _read_lc_variables(self) -> tuple:
        """Read all per-observation and per-star variables from the pipeline.

        Returns a tuple ``(lc_columns, scalars)``:

        * ``lc_columns`` — dict mapping VECTORTYPE_LC variable name → numpy
          array of length NJD (one value per observation).
        * ``scalars`` — dict mapping VECTORTYPE_SCALAR / PERSTARDATA /
          INLIST variable name → scalar value (Python ``float`` / ``int``).

        Data is copied out of vartools' internal storage before returning,
        so values remain valid after subsequent pipeline calls.
        """
        lib = _load_lib()
        VarInfo = lib._VartoolsVarInfo
        max_vars = 256
        vars_array = (VarInfo * max_vars)()
        n_vars = ctypes.c_int(0)

        overflow = lib.vartools_get_lc_variables(
            self._p, ctypes.byref(n_vars), vars_array, max_vars)
        if overflow > 0:
            max_vars = overflow
            vars_array = (VarInfo * max_vars)()
            lib.vartools_get_lc_variables(
                self._p, ctypes.byref(n_vars), vars_array, max_vars)

        lc_columns: dict = {}
        scalars: dict = {}
        for i in range(n_vars.value):
            vi = vars_array[i]
            vname = vi.name.decode("utf-8") if vi.name else None
            if vname is None or vi.dataptr is None:
                continue

            if vi.vectortype == VARTOOLS_VECTORTYPE_LC:
                length = vi.length
                if vi.datatype == VARTOOLS_TYPE_DOUBLE:
                    arr = np.ctypeslib.as_array(
                        (ctypes.c_double * length).from_address(vi.dataptr))
                    lc_columns[vname] = arr.copy()
                elif vi.datatype == VARTOOLS_TYPE_INT:
                    arr = np.ctypeslib.as_array(
                        (ctypes.c_int * length).from_address(vi.dataptr))
                    lc_columns[vname] = arr.astype(np.float64).copy()
                elif vi.datatype == VARTOOLS_TYPE_FLOAT:
                    arr = np.ctypeslib.as_array(
                        (ctypes.c_float * length).from_address(vi.dataptr))
                    lc_columns[vname] = arr.astype(np.float64).copy()

            elif vi.vectortype in (
                    VARTOOLS_VECTORTYPE_SCALAR,
                    VARTOOLS_VECTORTYPE_PERSTARDATA,
                    VARTOOLS_VECTORTYPE_INLIST):
                if vi.datatype == VARTOOLS_TYPE_DOUBLE:
                    scalars[vname] = float(
                        ctypes.c_double.from_address(vi.dataptr).value)
                elif vi.datatype == VARTOOLS_TYPE_INT:
                    scalars[vname] = int(
                        ctypes.c_int.from_address(vi.dataptr).value)
                elif vi.datatype == VARTOOLS_TYPE_LONG:
                    scalars[vname] = int(
                        ctypes.c_long.from_address(vi.dataptr).value)
                elif vi.datatype == VARTOOLS_TYPE_SHORT:
                    scalars[vname] = int(
                        ctypes.c_short.from_address(vi.dataptr).value)

        return lc_columns, scalars

    def read_capture(self, key: str) -> Optional[dict]:
        """Read back an in-memory LC snapshot taken by a ``-o <key> capture``
        command during the previous ``vartools_process_lc()`` call.

        Returns a dict mapping LC-variable name to numpy array (one entry
        per VECTORTYPE_LC variable defined at the snapshot point), or
        ``None`` when the slot was never filled.  Strings come back as
        Python lists of ``str``; numeric types as numpy arrays.

        The vartools-side data buffers are valid until the next
        ``process_lc()``/``process_lc_capture()`` call; this method
        copies everything out before returning, so the result remains
        usable across subsequent calls.
        """
        lib = _load_lib()
        VarInfo = lib._VartoolsVarInfo

        njd = lib.vartools_get_captured_njd(self._p, key.encode())
        if njd < 0:
            return None  # -1 (id not found) or -2 (slot never filled)

        max_vars = 256
        vars_array = (VarInfo * max_vars)()
        n_vars = ctypes.c_int(0)
        overflow = lib.vartools_get_captured_lc(
            self._p, key.encode(), vars_array, max_vars,
            ctypes.byref(n_vars))
        if overflow > 0:
            max_vars = overflow
            vars_array = (VarInfo * max_vars)()
            lib.vartools_get_captured_lc(
                self._p, key.encode(), vars_array, max_vars,
                ctypes.byref(n_vars))

        out: dict = {}
        for i in range(n_vars.value):
            vi = vars_array[i]
            vname = vi.name.decode("utf-8") if vi.name else None
            if vname is None or vi.dataptr is None:
                continue
            length = vi.length
            if vi.datatype == VARTOOLS_TYPE_DOUBLE:
                arr = np.ctypeslib.as_array(
                    (ctypes.c_double * length).from_address(vi.dataptr))
                out[vname] = arr.copy()
            elif vi.datatype == VARTOOLS_TYPE_INT:
                arr = np.ctypeslib.as_array(
                    (ctypes.c_int * length).from_address(vi.dataptr))
                out[vname] = arr.astype(np.float64).copy()
            elif vi.datatype == VARTOOLS_TYPE_FLOAT:
                arr = np.ctypeslib.as_array(
                    (ctypes.c_float * length).from_address(vi.dataptr))
                out[vname] = arr.astype(np.float64).copy()
            elif vi.datatype == VARTOOLS_TYPE_LONG:
                arr = np.ctypeslib.as_array(
                    (ctypes.c_long * length).from_address(vi.dataptr))
                out[vname] = arr.astype(np.int64).copy()
            elif vi.datatype == VARTOOLS_TYPE_SHORT:
                arr = np.ctypeslib.as_array(
                    (ctypes.c_short * length).from_address(vi.dataptr))
                out[vname] = arr.astype(np.int64).copy()
            elif vi.datatype == VARTOOLS_TYPE_STRING:
                # The captured slot stores strings as a malloc'd char **
                # (length njd) with one strdup'd entry per row.  Walk it
                # and decode into Python strings.
                ptr_arr = (ctypes.c_char_p * length).from_address(vi.dataptr)
                out[vname] = [
                    ptr_arr[r].decode("utf-8", errors="replace")
                    if ptr_arr[r] is not None else ""
                    for r in range(length)
                ]
        return out

    def _inject_lc_data(
        self,
        t: np.ndarray,
        mag: np.ndarray,
        err: np.ndarray,
        name: str,
        extra_columns: dict,
    ) -> None:
        """Inject t/mag/err plus extra named columns into the pipeline."""
        lib = _load_lib()

        col_names_list = ["t", "mag", "err"] + list(extra_columns.keys())
        col_data_list = [
            np.ascontiguousarray(t, dtype=np.float64),
            np.ascontiguousarray(mag, dtype=np.float64),
            np.ascontiguousarray(err, dtype=np.float64),
        ] + [
            np.ascontiguousarray(v, dtype=np.float64)
            for v in extra_columns.values()
        ]

        n_points = len(col_data_list[0])
        n_columns = len(col_names_list)

        c_names = (ctypes.c_char_p * n_columns)(
            *[s.encode() for s in col_names_list])
        c_ptrs = (ctypes.POINTER(ctypes.c_double) * n_columns)(
            *[arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
              for arr in col_data_list])

        rc = lib.vartools_set_lc_data(
            self._p, ctypes.c_int(n_points), ctypes.c_int(n_columns),
            c_names, c_ptrs, name.encode())
        if rc != 0:
            raise RuntimeError(f"vartools_set_lc_data returned error code {rc}")

    def __del__(self) -> None:
        if getattr(self, "_p", None) and _lib is not None:
            _lib.vartools_free_pipeline(self._p)
            self._p = None
