"""
Base class for all vartools command wrappers.

Auto-docstring mechanism
------------------------
Each subclass sets ``_vt_name`` to the vartools command name (e.g. ``"LS"``).
When the class is first *used* (lazy, not at import), ``_fetch_vt_help()``
calls ``vartools -help <name>`` and appends the output to the class
``__doc__``.  This means Python help() is always in sync with the installed
binary.  If the binary is not available the existing docstring is left as-is.

CLI rendering
-------------
Each subclass implements ``_to_cli_args()`` → list[str], returning the
command-line fragments for this command (e.g. ["-LS", "0.1", "10.0", ...]).

Output file declarations
------------------------
Commands that can produce output files declare them by implementing
``_output_file_specs()`` → dict[name, (suffix, ncols_hint)].  The Pipeline
uses this to locate and read output files when ``save_*=True`` is passed.
"""

from __future__ import annotations

import subprocess
import textwrap
from typing import ClassVar, Dict, List, Optional, Tuple


class VartoolsCommand:
    """Abstract base class for all pyvartools command wrappers."""

    # Subclasses must set this to the vartools command name, e.g. "LS"
    _vt_name: ClassVar[str] = ""

    # Set to True once the vartools help text has been appended to __doc__
    _help_fetched: ClassVar[bool] = False

    # ------------------------------------------------------------------
    # Auto-docstring (lazy, fetched on first call to help())
    # ------------------------------------------------------------------

    @classmethod
    def _fetch_vt_help(cls) -> None:
        """Append vartools -help output to this class's docstring (once)."""
        if cls._help_fetched or not cls._vt_name:
            return
        cls._help_fetched = True
        try:
            from pyvartools._binary import get_binary
            binary = get_binary()
            result = subprocess.run(
                [binary, "-help", f"-{cls._vt_name}"],
                capture_output=True, text=True, timeout=10,
            )
            vt_help = result.stdout.strip()
            if vt_help:
                separator = "\n\nvartools help\n" + "-" * 13 + "\n"
                existing = cls.__doc__ or ""
                cls.__doc__ = existing + separator + vt_help
            vt_examples = cls._fetch_vt_examples_text()
            if vt_examples:
                cls.__doc__ = (cls.__doc__ or "") + "\n\nExamples (vartools)\n" + "-" * 19 + "\n" + vt_examples
        except Exception:
            pass  # binary not found or other error; leave docstring unchanged

    @classmethod
    def _fetch_vt_examples_text(cls) -> str:
        """Return the text from vartools -example <name>, or empty string."""
        if not cls._vt_name:
            return ""
        try:
            from pyvartools._binary import get_binary
            binary = get_binary()
            result = subprocess.run(
                [binary, "-example", f"-{cls._vt_name}"],
                capture_output=True, text=True, timeout=10,
                # vartools informational commands exit non-zero; that is not an error
            )
            return result.stdout.strip()
        except Exception:
            return ""

    @classmethod
    def help(cls) -> None:
        """Print the full vartools help text for this command."""
        cls._fetch_vt_help()
        print(cls.__doc__ or f"No help available for -{cls._vt_name}.")

    @classmethod
    def examples(cls) -> None:
        """Print the vartools examples for this command."""
        text = cls._fetch_vt_examples_text()
        if text:
            print(text)
        else:
            print(f"No examples available for -{cls._vt_name}.")

    # ------------------------------------------------------------------
    # CLI rendering (subclasses implement this)
    # ------------------------------------------------------------------

    def _to_cli_args(self) -> List[str]:
        """Return command-line fragments for this command.

        Example: ["-LS", "0.1", "10.0", "0.1", "5", "1", "outdir", "whiten"]
        """
        raise NotImplementedError(
            f"{self.__class__.__name__} must implement _to_cli_args()"
        )

    def _to_cli_args_for_mode(self, mode: str) -> List[str]:
        """Return CLI fragments for this command, parameterised by run mode.

        ``mode`` is ``"single"`` when vartools will be invoked with ``-i``
        (run / run_file) and ``"list"`` for ``-l`` (run_filelist /
        run_batch / run_combinelcs).  The default implementation ignores
        ``mode`` and delegates to :meth:`_to_cli_args`; commands whose
        emission depends on the run mode (e.g. ``cmd.o`` whose output
        path is interpreted as a filename in single mode and as a
        directory in list mode) override this method instead.
        """
        return self._to_cli_args()

    def _to_cli_args_with_perlc(self, subs: dict,
                                mode: str = "single") -> List[str]:
        """Call _to_cli_args_for_mode(mode) with per-LC attrs temporarily replaced.

        Parameters
        ----------
        subs : dict mapping attr_name -> replacement_string
            e.g. {"minper": "list column 2", "maxper": "list column 3"}
        mode : str
            Forwarded to ``_to_cli_args_for_mode``.  Defaults to ``"single"``
            for callers that don't need mode-aware emission (e.g. legacy
            subprocess path); library-batch callers pass
            ``mode="library_batch"`` so that ``cmd.o`` emits ``forceoutdirmode``.
        """
        originals = {}
        for attr, replacement in subs.items():
            if hasattr(self, attr):
                originals[attr] = getattr(self, attr)
                object.__setattr__(self, attr, replacement)
        try:
            return self._to_cli_args_for_mode(mode)
        finally:
            for attr, orig in originals.items():
                object.__setattr__(self, attr, orig)

    # ------------------------------------------------------------------
    # Output file declarations (subclasses override as needed)
    # ------------------------------------------------------------------

    def _output_file_specs(self) -> Dict[str, Tuple[str, Optional[int]]]:
        """Return a dict describing output files this command may produce.

        Returns
        -------
        dict mapping logical_name → (file_suffix, ncols_hint)
            file_suffix : str
                Suffix appended to the light curve name to form the output
                filename, e.g. ".ls" for a Lomb-Scargle periodogram.
            ncols_hint : int or None
                Expected number of columns (for read_csv); None = auto-detect.
        """
        return {}

    # ------------------------------------------------------------------
    # Back-reference resolution for chained calls
    # ------------------------------------------------------------------

    def _resolve_back_references(self, prev) -> None:
        """Substitute CLI back-reference keywords with literal values from *prev*.

        When a command is issued as a continuation from a prior Result /
        BatchResult, string keywords like ``"ls"``, ``"aov"``, ``"bls"``,
        ``"both"``, ``"injectharm"``, and ``"fixcolumn <NAME>"`` that the
        vartools CLI would resolve against earlier commands in the same
        invocation cannot work out of the box — each chain step is a separate
        vartools invocation, so the "prior command" does not exist in the
        current one.

        Subclasses override this hook to rewrite their own parameters into
        numeric values pulled from *prev*.  The hook is called once, after
        ``__init__`` but before the command is placed in a Pipeline.  When
        *prev* is a ``BatchResult`` the substituted value may be a ``PerLC``
        array (one value per LC).

        The default implementation does nothing.
        """
        return None

    # ------------------------------------------------------------------
    # Requested output files tracking (set by Pipeline)
    # ------------------------------------------------------------------

    def _requested_outputs(self) -> List[str]:
        """Return the list of output file logical names requested by the user.

        Includes names where the file will be emitted (any non-suppressed mode)
        or where the command mandates output regardless of the save spec
        (e.g. ``autocorrelation``).
        """
        from .commands._helpers import _should_emit
        return [
            name for name in self._output_file_specs()
            if _should_emit(getattr(self, f"save_{name}", False))
               or getattr(self, "_mandatory_output", False)
        ]

    # ------------------------------------------------------------------
    # Repr
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        args = self._to_cli_args()
        return f"{self.__class__.__name__}({' '.join(args)})"
