"""
Binary discovery for the vartools executable.

Search order:
  1. Config file (path stored in pyvartools.config.set_binary())
  2. VARTOOLS_BINARY environment variable
  3. Install-time path recorded in _install_path.py by 'make install'
  4. PATH (shutil.which)
"""

import os
import re
import shutil
import subprocess
import sys
from typing import List

# Runtime override set by pyvartools.config.set_binary()
_configured_path: str = ""

# Runtime override set by pyvartools.config.set_library()
_configured_lib_path: str = ""

# Library filename: libtool produces .dylib on macOS, .so everywhere else.
_LIB_NAME = "libvartoolspipeline" + (".dylib" if sys.platform == "darwin" else ".so")


def set_binary(path: str) -> None:
    """Override the vartools binary path at runtime.

    This is the highest-priority entry in the discovery chain.  Call this
    at the top of your script or notebook if vartools is not on PATH and you
    have not set VARTOOLS_BINARY.

    Parameters
    ----------
    path : str
        Absolute path to the vartools executable.
    """
    global _configured_path
    _configured_path = path


def get_binary() -> str:
    """Return the path to the vartools binary.

    Raises
    ------
    FileNotFoundError
        If the binary cannot be found by any of the four mechanisms.
    """
    # 1. Runtime config (set_binary())
    if _configured_path:
        if os.path.isfile(_configured_path) and os.access(_configured_path, os.X_OK):
            return _configured_path
        raise FileNotFoundError(
            f"pyvartools.config.set_binary() was called with '{_configured_path}' "
            f"but that file does not exist or is not executable."
        )

    # 2. Environment variable
    env_path = os.environ.get("VARTOOLS_BINARY", "")
    if env_path:
        if os.path.isfile(env_path) and os.access(env_path, os.X_OK):
            return env_path
        raise FileNotFoundError(
            f"VARTOOLS_BINARY is set to '{env_path}' "
            f"but that file does not exist or is not executable."
        )

    # 3. Install-time path recorded by 'make install'
    try:
        from pyvartools._install_path import VARTOOLS_INSTALL_PATH
        if VARTOOLS_INSTALL_PATH:
            if os.path.isfile(VARTOOLS_INSTALL_PATH) and os.access(VARTOOLS_INSTALL_PATH, os.X_OK):
                return VARTOOLS_INSTALL_PATH
    except ImportError:
        pass

    # 4. PATH
    which = shutil.which("vartools")
    if which:
        return which

    raise FileNotFoundError(
        "Cannot find the vartools binary. Try one of:\n"
        "  - pyvartools.config.set_binary('/path/to/vartools')\n"
        "  - export VARTOOLS_BINARY=/path/to/vartools\n"
        "  - ensure vartools is on your PATH"
    )


def set_library(path: str) -> None:
    """Override the libvartoolspipeline path at runtime.

    Pass the absolute path to ``libvartoolspipeline.so`` (Linux) or
    ``libvartoolspipeline.dylib`` (macOS).
    """
    global _configured_lib_path
    _configured_lib_path = path


def _get_rpath_dirs(binary: str) -> List[str]:
    """Return RPATH / RUNPATH directories embedded in the vartools binary.

    On macOS the dynamic-loader equivalent is ``LC_RPATH`` load commands
    in the Mach-O header, surfaced by ``otool -l``.  Elsewhere we read
    ELF ``DT_RPATH`` / ``DT_RUNPATH`` tags via ``readelf -d``.  Both
    tools are best-effort: a missing binary helper just returns an
    empty list, leaving the rest of the discovery chain to handle the
    LC.
    """
    if sys.platform == "darwin":
        try:
            out = subprocess.check_output(
                ["otool", "-l", binary], stderr=subprocess.DEVNULL, text=True
            )
        except Exception:
            return []
        dirs: List[str] = []
        in_rpath = False
        for line in out.splitlines():
            if "LC_RPATH" in line:
                in_rpath = True
                continue
            if in_rpath:
                m = re.search(r'^\s*path\s+(.+?)(?:\s+\(offset|$)', line)
                if m:
                    dirs.append(m.group(1).strip())
                    in_rpath = False
        return dirs

    try:
        out = subprocess.check_output(
            ["readelf", "-d", binary], stderr=subprocess.DEVNULL, text=True
        )
    except Exception:
        return []
    dirs = []
    for line in out.splitlines():
        if "RPATH" in line or "RUNPATH" in line:
            m = re.search(r'\[(.+?)\]', line)
            if m:
                dirs.extend(m.group(1).split(":"))
    return dirs


def find_library() -> str:
    """Return the path to the libvartoolspipeline shared library.

    On macOS this resolves to ``libvartoolspipeline.dylib``; elsewhere
    to ``libvartoolspipeline.so``.

    Search order:
      1. Runtime config (set_library())
      2. VARTOOLS_LIBRARY environment variable
      3. RPATH dirs of the vartools binary
      4. Sibling lib/ directory relative to the binary
      5. System ldconfig / ctypes.util.find_library

    Raises
    ------
    FileNotFoundError
        If the library cannot be found.
    """
    import ctypes.util

    # 1. Runtime config
    if _configured_lib_path:
        if os.path.isfile(_configured_lib_path):
            return _configured_lib_path
        raise FileNotFoundError(
            f"pyvartools.config.set_library() was called with "
            f"'{_configured_lib_path}' but that file does not exist."
        )

    # 2. Environment variable
    env_path = os.environ.get("VARTOOLS_LIBRARY", "")
    if env_path:
        if os.path.isfile(env_path):
            return env_path
        raise FileNotFoundError(
            f"VARTOOLS_LIBRARY is set to '{env_path}' "
            f"but that file does not exist."
        )

    # 3. RPATH dirs from the vartools binary
    try:
        binary = get_binary()
        for d in _get_rpath_dirs(binary):
            candidate = os.path.join(d, _LIB_NAME)
            if os.path.isfile(candidate):
                return candidate
        # Also check common relative layouts: bin/../lib and bin/../data/*/
        bin_dir = os.path.dirname(os.path.realpath(binary))
        for rel in (".", "../lib", "../data/vartools/USERLIBS"):
            candidate = os.path.realpath(
                os.path.join(bin_dir, rel, _LIB_NAME)
            )
            if os.path.isfile(candidate):
                return candidate
    except Exception:
        pass

    # 4. ctypes.util.find_library (cross-platform; on macOS searches
    #    DYLD_LIBRARY_PATH and the standard dyld paths, on Linux searches
    #    LD_LIBRARY_PATH and the ldconfig cache).
    name = ctypes.util.find_library("vartoolspipeline")
    if name:
        return name

    ld_var = "DYLD_LIBRARY_PATH" if sys.platform == "darwin" else "LD_LIBRARY_PATH"
    raise FileNotFoundError(
        f"Cannot find {_LIB_NAME}.  Try one of:\n"
        f"  - pyvartools.config.set_library('/path/to/{_LIB_NAME}')\n"
        f"  - export VARTOOLS_LIBRARY=/path/to/{_LIB_NAME}\n"
        f"  - add the install dir to {ld_var}\n"
        "  - run 'make install' to install the library"
    )
