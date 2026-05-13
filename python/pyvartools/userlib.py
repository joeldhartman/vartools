"""Support for vartools user extension libraries (.so / .la).

User-developed vartools extensions are compiled shared libraries that
export five functions: ``$LIBNAME_Initialize``, ``$LIBNAME_ParseCL``,
``$LIBNAME_ShowSyntax``, ``$LIBNAME_ShowHelp``, ``$LIBNAME_RunCommand``.
They are loaded at runtime via ``vartools -L $path/$lib.so -$cmdname ...``,
or auto-loaded without ``-L`` when installed in the vartools userlibs dir.

Usage patterns
--------------
1. Quick one-off with raw args::

    pipe = vt.Pipeline([
        vt.UserCommand("USERLIBS/src/stitch.so", "stitch",
                       "mag err mask lcnum median")
    ])

2. Named class via factory::

    Stitch = vt.load_userlib("USERLIBS/src/stitch.so")
    pipe = vt.Pipeline([Stitch("mag err mask lcnum median")])

3. Auto-discover all installed extensions::

    cmds = vt.discover_userlibs()
    pipe = vt.Pipeline([cmds["stitch"]("mag err mask lcnum median")])

4. Full Python wrapper (subclass)::

    class Stitch(vt.UserCommand):
        _lib_path = "/usr/local/share/vartools/userlibs/stitch.so"
        _cmd_name = "stitch"
        def __init__(self, variables, errors, masks, lcnum, method="median"):
            super().__init__(self._lib_path, self._cmd_name,
                             f"{variables} {errors} {masks} {lcnum} {method}")
"""

from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Union

from pyvartools._command import VartoolsCommand

logger = logging.getLogger(__name__)


def _fetch_userlib_text(
    lib_path: Optional[str], cmd_name: str, flag: str
) -> str:
    """Run ``vartools [-L lib_path] <flag> -cmd_name`` and return stdout.

    Parameters
    ----------
    lib_path : str or None
        Path to the .so file, or ``None`` if the library is auto-loaded.
    cmd_name : str
        Command name (without leading dash).
    flag : str
        Either ``"-help"`` or ``"-example"``.

    Returns
    -------
    str
        Combined stdout from vartools, or empty string on failure.
    """
    if not cmd_name:
        return ""
    try:
        from pyvartools._binary import get_binary
        binary = get_binary()
        cmd: List[str] = [binary]
        if lib_path:
            cmd += ["-L", lib_path]
        cmd += [flag, f"-{cmd_name}"]
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=10
        )
        return result.stdout.strip()
    except Exception:
        return ""


class UserCommand(VartoolsCommand):
    """Wraps a user-developed vartools extension library command.

    Parameters
    ----------
    lib_path : str or Path or None
        Path to the compiled user library (.so or .la).  If the library is
        installed in the vartools userlibs directory it will be auto-loaded
        by vartools and *lib_path* may be omitted (pass ``None`` or ``""``).
    name : str, optional
        Command name (e.g. ``"stitch"``).  Defaults to the library filename
        stem when *lib_path* is given.
    args : list of str or str
        Raw CLI tokens for this command (everything after the command name).
        A plain string is split on whitespace.

    Notes
    -----
    When a :class:`UserCommand` is present in a :class:`~pyvartools.Pipeline`,
    the pipeline always uses subprocess mode (never in-process library mode),
    because the in-process library does not support dynamically loaded
    extension libraries.

    Duplicate ``-L`` flags in a single vartools invocation are safe; vartools
    detects and skips already-loaded libraries.
    """

    _vt_name = ""

    def __init__(
        self,
        lib_path: Union[str, Path, None],
        name: Optional[str] = None,
        args: Union[str, List[str]] = (),
    ) -> None:
        self._lib_path = str(lib_path) if lib_path else None
        if name is None and lib_path:
            name = Path(lib_path).name.split(".")[0]
        self._cmd_name = name or ""
        self._vt_name = self._cmd_name
        self._user_args: List[str] = (
            args.split() if isinstance(args, str) else list(args)
        )

    def _to_cli_args(self) -> List[str]:
        prefix = ["-L", self._lib_path] if self._lib_path else []
        return prefix + [f"-{self._cmd_name}"] + self._user_args

    def _output_file_specs(self) -> dict:
        return {}

    def help(self) -> None:  # type: ignore[override]
        """Fetch and print vartools help for this user command."""
        text = _fetch_userlib_text(self._lib_path, self._cmd_name, "-help")
        if text:
            print(text)
        else:
            print(f"No help available for -{self._cmd_name}.")

    def examples(self) -> None:  # type: ignore[override]
        """Fetch and print vartools examples for this user command."""
        text = _fetch_userlib_text(self._lib_path, self._cmd_name, "-example")
        if text:
            print(text)
        else:
            print(f"No examples available for -{self._cmd_name}.")


def load_userlib(
    lib_path: Union[str, Path],
    name: Optional[str] = None,
    cls_name: Optional[str] = None,
) -> type:
    """Create a :class:`UserCommand` subclass for a user extension library.

    The returned class constructor accepts a single positional argument *args*
    (raw CLI tokens for the command), with *lib_path* and *name* pre-bound.
    The class can itself be subclassed to write a full Python wrapper.

    Parameters
    ----------
    lib_path : str or Path
        Path to the compiled user library (.so or .la).
    name : str, optional
        Command name.  Defaults to the library filename stem.
    cls_name : str, optional
        Python class name for the returned type.  Defaults to the
        title-cased command name (e.g. ``"stitch"`` → ``"Stitch"``).

    Returns
    -------
    type
        A :class:`UserCommand` subclass.  Can be instantiated directly
        or further subclassed for a full Python wrapper.

    Examples
    --------
    ::

        Stitch = vt.load_userlib("USERLIBS/src/stitch.so")
        pipe = vt.Pipeline([Stitch("mag err mask lcnum median")])

        # Subclass for a full typed wrapper
        class Stitch(vt.load_userlib("stitch.so", name="stitch")):
            def __init__(self, variables, errors, masks, lcnum,
                         method="median", **kwargs):
                super().__init__(
                    f"{variables} {errors} {masks} {lcnum} {method}")
    """
    # Use absolute() rather than resolve() so symlink names like "stitch.so"
    # are preserved.  vartools' libtool dlopen needs the basename without a
    # version suffix (e.g. ".so.0.0.0") for symbol lookup to find
    # `<libbasename>_Initialize`, so following the symlink to the underlying
    # versioned filename breaks the load.
    resolved = Path(lib_path).absolute()
    if name is None:
        name = resolved.name.split(".")[0]
    if cls_name is None:
        cls_name = name.title()

    doc = _fetch_userlib_text(str(resolved), name, "-help") or (
        f"User extension command -{name} loaded from {resolved}."
    )

    _lib_path_str = str(resolved)
    _cmd_name = name

    def __init__(self, args: Union[str, List[str]] = ()) -> None:  # type: ignore[misc]
        UserCommand.__init__(self, _lib_path_str, _cmd_name, args)

    cls = type(
        cls_name,
        (UserCommand,),
        {
            "__init__": __init__,
            "__doc__": doc,
            "_vt_name": _cmd_name,
        },
    )
    return cls


def discover_userlibs(
    search_paths: Optional[List[Union[str, Path]]] = None,
) -> Dict[str, type]:
    """Return a dict mapping command name → :class:`UserCommand` subclass for
    all user libraries found in the given directories.

    Default search order:

    1. Paths listed in ``$VARTOOLS_USERLIBS`` (colon-separated).
    2. ``$prefix/share/vartools/userlibs/`` derived from the binary path.
    3. Any paths passed explicitly via *search_paths*.

    Parameters
    ----------
    search_paths : list of str or Path, optional
        Additional directories to search.  Pass ``[]`` to skip the default
        directories and search only these paths.

    Returns
    -------
    dict[str, type]
        e.g. ``{'stitch': <class Stitch>, 'macula': <class Macula>}``
    """
    dirs: List[Path] = []

    # 1. $VARTOOLS_USERLIBS
    env_val = os.environ.get("VARTOOLS_USERLIBS", "")
    for p in env_val.split(":"):
        p = p.strip()
        if p:
            dirs.append(Path(p))

    # 2. Derive from binary location: e.g. /usr/local/bin/vartools
    #    → /usr/local/share/vartools/userlibs/
    try:
        from pyvartools._binary import get_binary
        binary = get_binary()
        bin_dir = Path(binary).resolve().parent
        for rel in (
            "../share/vartools/userlibs",
            "../../share/vartools/userlibs",
        ):
            candidate = (bin_dir / rel).resolve()
            if candidate.is_dir():
                dirs.append(candidate)
    except Exception:
        pass

    # 3. Explicit search_paths argument
    if search_paths:
        for p in search_paths:
            dirs.append(Path(p))

    # Deduplicate while preserving order
    seen: set = set()
    unique_dirs: List[Path] = []
    for d in dirs:
        key = str(d.resolve()) if d.exists() else str(d)
        if key not in seen:
            seen.add(key)
            unique_dirs.append(d)

    result: Dict[str, type] = {}
    for d in unique_dirs:
        if not d.is_dir():
            continue
        candidates = (
            sorted(d.glob("*.so"))
            + sorted(d.glob("*.so.*"))
            + sorted(d.glob("*.la"))
        )
        for so in candidates:
            stem = so.name.split(".")[0]
            if not stem or stem in result:
                continue
            try:
                cls = load_userlib(so, name=stem)
                result[stem] = cls
            except Exception as exc:
                logger.debug("discover_userlibs: skipping %s: %s", so, exc)

    return result
