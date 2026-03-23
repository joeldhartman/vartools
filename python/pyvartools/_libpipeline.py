"""ctypes binding to libvartoolspipeline — init-once / process-many API.

Usage (low-level)::

    lp = LibPipeline(["-rms", "-oneline"])
    stats = lp.process_lc(t, mag, err, name="my_lc")  # pd.Series
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

    _lib = lib
    return _lib


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

        self._p = lib.vartools_init_pipeline(argc, argv_c)
        if not self._p:
            raise RuntimeError(
                "vartools_init_pipeline failed — check command-line arguments."
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

    def __del__(self) -> None:
        if getattr(self, "_p", None) and _lib is not None:
            _lib.vartools_free_pipeline(self._p)
            self._p = None
