"""In-process callback for the ``-python`` command.

This module installs a C-callable callback into libvartoolspipeline.so so
that vartools' ``-python`` command can execute user code in *this* Python
interpreter rather than a per-thread sub-process spawned by vartools.

The callback runs user code via ``exec()`` against a caller-supplied
globals dict (``__main__.__dict__`` by default), so the user code can
see imports and variables defined in the host script.

This module is imported on demand by ``pyvartools.commands.python`` when
``inprocess=True`` is requested.  Import is *not* automatic at package
load time, so users who never opt in pay no init cost and don't risk a
spurious dlopen of libvartoolsrunpython.so.

Coverage (v1):
    * Numeric scalars and LC vectors of types DOUBLE, FLOAT, INT, LONG.
    * Init code (run once per command, before per-LC bodies).
    * ``process_all_lcs=False`` only.
    * ``continueprocess`` not supported (vartools falls back to
      subprocess; the wrapper checks at construction time).

The matching C side lives in ``src/runpython.c`` (see
``vt_python_run_inprocess``).  The callback signature is documented in
``src/vartools.h``.
"""
from __future__ import annotations

import ctypes
import sys
import threading
from typing import Optional

import numpy as np

from . import _binary

# vartools type-code constants (mirror src/programdata.h).
VARTOOLS_TYPE_DOUBLE     = 0
VARTOOLS_TYPE_STRING     = 1
VARTOOLS_TYPE_INT        = 2
VARTOOLS_TYPE_FLOAT      = 3
VARTOOLS_TYPE_LONG       = 4
VARTOOLS_TYPE_CHAR       = 5
VARTOOLS_TYPE_CONVERTJD  = 6
VARTOOLS_TYPE_SHORT      = 7

# Map vartools type → numpy dtype (for vector views).
_NUMPY_DTYPE = {
    VARTOOLS_TYPE_DOUBLE:    np.float64,
    VARTOOLS_TYPE_CONVERTJD: np.float64,
    VARTOOLS_TYPE_FLOAT:     np.float32,
    VARTOOLS_TYPE_INT:       np.int32,
    VARTOOLS_TYPE_LONG:      np.int64,
}

# Map vartools type → ctypes scalar type (for scalar pointer derefs).
_CTYPE_SCALAR = {
    VARTOOLS_TYPE_DOUBLE:    ctypes.c_double,
    VARTOOLS_TYPE_CONVERTJD: ctypes.c_double,
    VARTOOLS_TYPE_FLOAT:     ctypes.c_float,
    VARTOOLS_TYPE_INT:       ctypes.c_int,
    VARTOOLS_TYPE_LONG:      ctypes.c_long,
}

# Per-command cached globals dict (so init code state persists across
# per-LC calls within the same -python command).  Keyed on the C-side
# command_id parameter.
_command_globals: dict[int, dict] = {}
_command_globals_lock = threading.Lock()

# Holds the registered CFUNCTYPE callback so Python doesn't garbage-
# collect it while the C side still has the pointer.
_callback_keepalive: Optional[ctypes._CFuncPtr] = None
_lib_handle: Optional[ctypes.CDLL] = None


# C callback signature (must match the typedef in src/vartools.h).
_PYTHON_CALLBACK_T = ctypes.CFUNCTYPE(
    ctypes.c_int,                         # return: 0 OK, 1 user error
    ctypes.c_void_p,                      # namespace_ptr (PyObject* opaque)
    ctypes.c_int,                         # command_id
    ctypes.c_char_p,                      # code (NUL-terminated UTF-8)
    ctypes.c_int,                         # is_init flag
    ctypes.c_int,                         # n_invars
    ctypes.POINTER(ctypes.c_char_p),      # invar_names
    ctypes.POINTER(ctypes.c_int),         # invar_types
    ctypes.POINTER(ctypes.c_int),         # invar_lengths
    ctypes.POINTER(ctypes.c_void_p),      # invar_data
    ctypes.c_int,                         # n_outvars
    ctypes.POINTER(ctypes.c_char_p),      # outvar_names
    ctypes.POINTER(ctypes.c_int),         # outvar_types
    ctypes.POINTER(ctypes.c_int),         # outvar_lengths
    ctypes.POINTER(ctypes.c_void_p),      # outvar_data
    ctypes.c_void_p,                      # error_buf (writeable)
    ctypes.c_int,                         # error_buf_size
)


def _read_invar(name: str, dtype_code: int, length: int, data_ptr: int):
    """Materialise an invar into a Python value (scalar) or numpy view (vector)."""
    if length <= 1:
        # Scalar — read a single value through the typed scalar pointer.
        ctype = _CTYPE_SCALAR.get(dtype_code)
        if ctype is None:
            raise TypeError(
                f"in-process -python: scalar invar '{name}' has unsupported "
                f"vartools datatype {dtype_code}"
            )
        return ctype.from_address(data_ptr).value

    # Vector — return a writable numpy view over the C buffer so the user
    # can mutate it in place (and any writes flow back to vartools at no
    # extra cost).
    dtype = _NUMPY_DTYPE.get(dtype_code)
    if dtype is None:
        raise TypeError(
            f"in-process -python: vector invar '{name}' has unsupported "
            f"vartools datatype {dtype_code}"
        )
    arr_type = ctypes.c_byte * (length * np.dtype(dtype).itemsize)
    raw = arr_type.from_address(data_ptr)
    return np.frombuffer(raw, dtype=dtype, count=length)


def _write_outvar(name: str, value, dtype_code: int, length: int, data_ptr: int):
    """Write a Python value back to a C buffer for a named outvar."""
    if length <= 1:
        ctype = _CTYPE_SCALAR.get(dtype_code)
        if ctype is None:
            raise TypeError(
                f"in-process -python: scalar outvar '{name}' has unsupported "
                f"vartools datatype {dtype_code}"
            )
        # Python scalar (or 0-d numpy) → C scalar.
        ctype.from_address(data_ptr).value = ctype(value).value
        return

    dtype = _NUMPY_DTYPE.get(dtype_code)
    if dtype is None:
        raise TypeError(
            f"in-process -python: vector outvar '{name}' has unsupported "
            f"vartools datatype {dtype_code}"
        )
    target_type = ctypes.c_byte * (length * np.dtype(dtype).itemsize)
    raw = target_type.from_address(data_ptr)
    target = np.frombuffer(raw, dtype=dtype, count=length)
    src = np.asarray(value, dtype=dtype)
    if src.shape != (length,):
        raise ValueError(
            f"in-process -python: outvar '{name}' must be a vector of "
            f"length {length}, got shape {src.shape}"
        )
    target[:] = src


def _python_callback_impl(namespace_ptr, command_id, code_bytes, is_init,
                          n_invars, invar_names, invar_types,
                          invar_lengths, invar_data,
                          n_outvars, outvar_names, outvar_types,
                          outvar_lengths, outvar_data,
                          error_buf, error_buf_size) -> int:
    """The actual callback body — invoked from C.

    Returns 0 on success, 1 on user error (with the message written into
    the error buffer for vartools to print).
    """
    try:
        # Globals dict to use for this command.  In v1 each command shares
        # the same caller-supplied namespace (typically __main__), but we
        # still cache the dict per command_id so init writes don't leak
        # across distinct -python commands that ran concurrently.
        with _command_globals_lock:
            globals_dict = _command_globals.get(command_id)
            if globals_dict is None:
                globals_dict = sys.modules["__main__"].__dict__
                _command_globals[command_id] = globals_dict
        # numpy is the conventional helper exposed to vartools -python
        # bodies; ensure it's importable in the user's namespace.
        globals_dict.setdefault("numpy", np)

        code = (code_bytes or b"").decode("utf-8")

        if is_init:
            # Init code runs at module scope — pass globals as locals so
            # `def`/`import` statements land in the shared namespace and
            # are visible to subsequent per-LC bodies.  Init has no
            # invars/outvars (vartools always calls us with n=0).
            exec(code, globals_dict, globals_dict)
            return 0

        # Build a locals dict with the named invars.
        locals_dict: dict = {}
        for i in range(n_invars):
            name = invar_names[i].decode("utf-8")
            t    = invar_types[i]
            l    = invar_lengths[i]
            ptr  = invar_data[i]
            locals_dict[name] = _read_invar(name, t, l, ptr)

        # Run the user code.  exec() in (globals, locals) form: top-level
        # assignments land in locals_dict by default, which matches what
        # vartools' subprocess path does (each call's body runs in a
        # function with its own locals).  Functions defined by init
        # earlier are visible via globals_dict.
        exec(code, globals_dict, locals_dict)

        # Read outvars from the post-exec namespace.  Pull from locals
        # first, then fall back to globals so a `global foo; foo = …`
        # body still works.
        for i in range(n_outvars):
            name = outvar_names[i].decode("utf-8")
            t    = outvar_types[i]
            l    = outvar_lengths[i]
            ptr  = outvar_data[i]
            if name in locals_dict:
                value = locals_dict[name]
            elif name in globals_dict:
                value = globals_dict[name]
            else:
                raise NameError(
                    f"in-process -python: outvar '{name}' was not assigned "
                    f"in the user code"
                )
            _write_outvar(name, value, t, l, ptr)

        return 0

    except Exception as exc:
        # Format the error and copy into the C-side buffer.  error_buf
        # is a void* address; write the bytes + a NUL terminator so the
        # C side can fprintf() it as a plain string.
        msg = f"{type(exc).__name__}: {exc}".encode("utf-8")
        n = min(len(msg), error_buf_size - 1)
        if error_buf:
            ctypes.memmove(error_buf, msg, n)
            # NUL-terminate one byte past the message.
            ctypes.memmove(error_buf + n, b"\x00", 1)
        return 1


def register(library_path: Optional[str] = None) -> None:
    """Register the in-process callback with libvartoolspipeline.so.

    Idempotent: subsequent calls are no-ops.  Resolves the library via
    ``pyvartools._binary.get_library()`` unless an explicit path is
    supplied.
    """
    global _callback_keepalive, _lib_handle
    if _callback_keepalive is not None:
        return
    path = library_path or _binary.find_library()
    _lib_handle = ctypes.CDLL(path)
    _lib_handle.vartools_register_python_callback.argtypes = [_PYTHON_CALLBACK_T]
    _lib_handle.vartools_register_python_callback.restype = None
    _lib_handle.vartools_set_python_namespace.argtypes = [ctypes.c_void_p]
    _lib_handle.vartools_set_python_namespace.restype = None

    cb = _PYTHON_CALLBACK_T(_python_callback_impl)
    _lib_handle.vartools_register_python_callback(cb)
    _callback_keepalive = cb  # prevent GC while C holds the pointer


def set_namespace(namespace: dict) -> None:
    """Set the namespace dict the callback will use for ``exec()``.

    The caller passes any ``dict`` (typically ``__main__.__dict__`` or a
    custom sandbox); it's installed as the globals for subsequent user
    code.  Stored on the C side as an opaque ``void*`` and round-tripped
    back into the callback as ``namespace_ptr`` (currently the callback
    ignores this argument and uses its own per-command cache, but the
    pointer is held to keep the dict alive).
    """
    if _lib_handle is None:
        register()
    # Stash the dict in our own per-command map so the callback finds it.
    # We cast id() of the dict to a void*; the C side never dereferences
    # it, only stores it.
    _lib_handle.vartools_set_python_namespace(  # type: ignore[union-attr]
        ctypes.c_void_p(id(namespace))
    )
    # Also overwrite the cached default globals so the next callback call
    # picks up the new dict.  Use command_id=0 since there's only one
    # active -python command in the inprocess v1.
    with _command_globals_lock:
        _command_globals[0] = namespace
        # Make sure numpy is importable.
        namespace.setdefault("numpy", np)


def reset_command_globals() -> None:
    """Clear the per-command globals cache (mostly useful for tests)."""
    with _command_globals_lock:
        _command_globals.clear()
