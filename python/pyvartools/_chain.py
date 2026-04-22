"""
CommandChain: internal deferred sequence of vartools commands.

This class is used internally by ``LightCurveBatch`` to accumulate a list of
commands and execute them as a single vartools invocation per light curve.
It is not part of the public pyvartools API.

For single-LC work, call command methods directly on ``LightCurve`` or
``Result`` — they execute immediately and return a ``Result``.
For batch work, use ``LightCurveBatch``.
"""

from __future__ import annotations

import copy
import inspect
from typing import Any, Dict, List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from .lightcurve import LightCurve
    from .results import Result
    from ._command import VartoolsCommand


# ---------------------------------------------------------------------------
# Lambda helpers
# ---------------------------------------------------------------------------

def _has_lambda_params(cmd: "VartoolsCommand") -> bool:
    """Return True if any instance attribute of cmd is a callable (lambda/fn)."""
    return any(
        inspect.isfunction(v)
        for v in vars(cmd).values()
    )


def _resolve_lambdas(cmd: "VartoolsCommand", result: "Result") -> "VartoolsCommand":
    """Return a copy of cmd with any callable attribute resolved by calling it
    with *result*.

    The original cmd is not mutated.
    """
    resolved = copy.copy(cmd)
    for attr, val in vars(cmd).items():
        if inspect.isfunction(val):
            object.__setattr__(resolved, attr, val(result))
    return resolved


# ---------------------------------------------------------------------------
# CommandChain
# ---------------------------------------------------------------------------

class CommandChain:
    """A deferred sequence of vartools commands bound to a source LightCurve.

    Instances are created by calling a command method on a ``LightCurve`` or
    ``Result``::

        chain = lc.LS(0.1, 10.0)        # CommandChain with one command
        chain = chain.Killharm(nharm=2)  # CommandChain with two commands
        result = chain.run()             # execute

    Command methods are attached at import time by
    ``_method_gen._attach_command_methods()``.
    """

    def __init__(
        self,
        lc: "LightCurve",
        commands: Optional[List["VartoolsCommand"]] = None,
        global_opts: Optional[Dict[str, Any]] = None,
    ) -> None:
        self._lc = lc
        self._commands: List["VartoolsCommand"] = list(commands or [])
        self._global_opts: Dict[str, Any] = dict(global_opts or {})

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _append(self, cmd: "VartoolsCommand") -> "CommandChain":
        """Return a new CommandChain with *cmd* appended."""
        return CommandChain(self._lc, self._commands + [cmd], self._global_opts)

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run(
        self,
        capture_lc: bool = True,
        timeout: Optional[int] = None,
        init_lc_vars=None,
        randseed: Optional[int] = None,
        skipmissing: bool = False,
        jdtol: Optional[float] = None,
        matchstringid: bool = False,
    ) -> "Result":
        """Execute the accumulated command chain.

        Parameters
        ----------
        capture_lc : bool
            Whether to capture the output light curve into ``result.lc``.
            Defaults to ``True`` (unlike ``Pipeline.run()`` which defaults to
            ``False``) so the result can be used for further chaining.
        timeout : int, optional
            Maximum seconds to wait for each vartools invocation.
        init_lc_vars : dict, optional
            Per-observation variables (passed to Pipeline.run).
        randseed : int, optional
            ``-randseed N`` for reproducible random-number sequences.
        skipmissing : bool
            ``-skipmissing`` flag.
        jdtol : float, optional
            ``-jdtol N`` tolerance for Julian-date matching.
        matchstringid : bool
            ``-matchstringid`` flag.

        Returns
        -------
        Result
        """
        # Merge global_opts set via with_options(), then override with
        # any kwargs passed directly to run().  All non-None / explicitly
        # supplied values override the stored options, including False for
        # booleans (so skipmissing=False can clear a with_options default).
        run_kwargs: Dict[str, Any] = dict(self._global_opts)
        run_kwargs.pop("capture_lc", None)   # handled as explicit parameter
        if timeout is not None:
            run_kwargs["timeout"] = timeout
        if init_lc_vars is not None:
            run_kwargs["init_lc_vars"] = init_lc_vars
        if randseed is not None:
            run_kwargs["randseed"] = randseed
        run_kwargs["skipmissing"] = skipmissing
        if jdtol is not None:
            run_kwargs["jdtol"] = jdtol
        run_kwargs["matchstringid"] = matchstringid

        # Check for any lambda-bearing commands; split the chain if found.
        lambda_idx = next(
            (i for i, c in enumerate(self._commands) if _has_lambda_params(c)),
            None,
        )

        if lambda_idx is None:
            return self._run_segment(
                self._lc, self._commands, capture_lc=capture_lc, **run_kwargs
            )
        else:
            return self._run_with_lambda_split(
                capture_lc=capture_lc, **run_kwargs
            )

    def _run_segment(
        self,
        lc: "LightCurve",
        commands: List["VartoolsCommand"],
        capture_lc: bool = True,
        **run_kwargs,
    ) -> "Result":
        """Run a contiguous list of commands as a single Pipeline invocation."""
        from .pipeline import Pipeline
        return Pipeline(commands).run(lc, capture_lc=capture_lc, **run_kwargs)

    def _run_with_lambda_split(
        self,
        capture_lc: bool = True,
        **run_kwargs,
    ) -> "Result":
        """Execute the chain, splitting at each lambda-bearing command.

        Each split fires a separate vartools invocation.  The result from each
        segment is used to resolve the lambdas in the next lambda-bearing
        command.
        """
        current_lc = self._lc
        current_result = None
        remaining = list(self._commands)

        while remaining:
            # Find the next lambda boundary
            lambda_idx = next(
                (i for i, c in enumerate(remaining) if _has_lambda_params(c)),
                None,
            )

            if lambda_idx is None:
                # No more lambdas — run the rest in one shot.
                current_result = self._run_segment(
                    current_lc, remaining,
                    capture_lc=capture_lc, **run_kwargs
                )
                remaining = []
            elif lambda_idx == 0:
                # The very first command has a lambda — we need a preceding
                # result to resolve it.  If there is none, that means the
                # chain started with a lambda, which is a user error.
                if current_result is None:
                    raise ValueError(
                        "The first command in the chain has a lambda parameter "
                        "but there is no preceding result to resolve it from. "
                        "Run at least one command before using a lambda."
                    )
                resolved = _resolve_lambdas(remaining[0], current_result)
                remaining = [resolved] + remaining[1:]
            else:
                # Run everything up to (not including) the lambda command,
                # capturing the LC so it can be used in the next segment.
                current_result = self._run_segment(
                    current_lc, remaining[:lambda_idx],
                    capture_lc=True, **run_kwargs
                )
                if current_result.lc is None:
                    raise RuntimeError(
                        "Could not capture output light curve between lambda "
                        "split points.  Make sure the pipeline is not "
                        "suppressing LC output."
                    )
                current_lc = current_result.lc
                # Resolve the lambda command and continue
                resolved = _resolve_lambdas(remaining[lambda_idx], current_result)
                remaining = [resolved] + remaining[lambda_idx + 1:]

        return current_result

    # ------------------------------------------------------------------
    # Repr
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        cmds = " → ".join(c.__class__.__name__ for c in self._commands)
        return f"CommandChain([{cmds}], lc={self._lc!r})"
