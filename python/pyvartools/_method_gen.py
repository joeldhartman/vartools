"""
Programmatic generation of command methods for LightCurve, Result,
LightCurveBatch, and BatchResult.

Called once at import time from pyvartools/__init__.py via
``_attach_all_command_methods()``.

LightCurve and Result
---------------------
Every command is attached as a single **immediate** method that runs the
command and returns a ``Result``.  There is no deferred / chain-building
variant on these classes — run options are passed as keyword arguments:

    result = lc.LS(0.5, 10.0, 1e-3)
    result = lc.LS(0.5, 10.0, 1e-3, randseed=42, capture_lc=False)
    result = result.rms()            # prior vars are preserved

Pipeline-stateful commands (``savelc``, ``restorelc``, ``columnsuffix``,
``ifcmd``, ``o``) only work correctly within a single vartools invocation
and therefore raise ``NotImplementedError`` when called on ``LightCurve``
or ``Result``.  Use ``Pipeline([...]).run(lc)`` for these.

LightCurveBatch and BatchResult
---------------------------------
Both a deferred chain-building method and an immediate method are provided:

    batch.LS(...)      → ``LightCurveBatch`` (deferred; call ``.run()``)
    batch.run_LS(...)  → ``BatchResult`` (immediate)

Run-time options
----------------
The following keyword argument names are reserved for pipeline run options
and are separated from command parameters automatically:

    capture_lc, timeout, init_lc_vars, randseed,
    skipmissing, jdtol, matchstringid
"""

from __future__ import annotations

import inspect
from typing import FrozenSet, List, Type, TYPE_CHECKING

if TYPE_CHECKING:
    from ._command import VartoolsCommand

# ---------------------------------------------------------------------------
# Reserved run-option kwarg names
# ---------------------------------------------------------------------------

_RUN_OPT_NAMES: FrozenSet[str] = frozenset({
    "capture_lc",
    "timeout",
    "init_lc_vars",
    "randseed",
    "skipmissing",
    "jdtol",
    "matchstringid",
})

# Commands that are only meaningful within a single vartools invocation.
# Calling them directly on LightCurve or Result raises NotImplementedError.
_STATEFUL: FrozenSet[str] = frozenset({
    "savelc",
    "restorelc",
    "columnsuffix",
    "if",      # ifcmd class, _vt_name="if"
    "elif",    # elifcmd
    "else",    # elsecmd
    "fi",      # ficmd
    "o",
})


def _split_kwargs(kwargs: dict) -> tuple:
    """Split kwargs into (run_opts, cmd_kwargs)."""
    run_opts = {k: v for k, v in kwargs.items() if k in _RUN_OPT_NAMES}
    cmd_kwargs = {k: v for k, v in kwargs.items() if k not in _RUN_OPT_NAMES}
    return run_opts, cmd_kwargs


# ---------------------------------------------------------------------------
# Method factories — LightCurve and Result (immediate execution)
# ---------------------------------------------------------------------------

def _make_immediate_lc(cmd_cls: Type["VartoolsCommand"]):
    """Immediate method for LightCurve: lc.LS(...) → Result."""
    vt_name = cmd_cls._vt_name

    def _immediate(self, *args, **kwargs):
        run_opts, cmd_kwargs = _split_kwargs(kwargs)
        capture_lc = run_opts.pop("capture_lc", True)
        cmd = cmd_cls(*args, **cmd_kwargs)
        from .pipeline import Pipeline
        return Pipeline([cmd]).run(self, capture_lc=capture_lc, **run_opts)

    _immediate.__name__ = vt_name
    _immediate.__qualname__ = f"LightCurve.{vt_name}"
    _immediate.__doc__ = (
        f"Run -{vt_name} on this light curve and return a Result.\n\n"
        f"Run-time options (capture_lc, randseed, timeout, skipmissing, "
        f"jdtol, matchstringid) may be passed as keyword arguments alongside "
        f"the command parameters.\n\n"
        + (inspect.cleandoc(cmd_cls.__doc__) if cmd_cls.__doc__ else "")
    )
    return _immediate


def _make_immediate_result(cmd_cls: Type["VartoolsCommand"]):
    """Immediate method for Result: result.LS(...) → Result on result.lc,
    with prior vars preserved."""
    vt_name = cmd_cls._vt_name

    def _immediate(self, *args, **kwargs):
        if self.lc is None:
            raise AttributeError(
                f"Cannot call {vt_name}: result.lc is None.  "
                f"The preceding run must capture the output light curve "
                f"(capture_lc=True, which is the default)."
            )
        run_opts, cmd_kwargs = _split_kwargs(kwargs)
        capture_lc = run_opts.pop("capture_lc", True)
        cmd = cmd_cls(*args, **cmd_kwargs)
        # Back-reference keywords such as ``period="ls"`` / ``"aov"`` /
        # ``"bls"`` / ``"fixcolumn NAME"`` only work within a single vartools
        # invocation.  Chain boundaries would invalidate them, so resolve
        # them now against the prior Result before the command runs.
        cmd._resolve_back_references(self)
        # Carry forward prior OUTCOLUMN values and per-star scalars so that
        # the next segment can reference them via injected ``-expr const``.
        from .lightcurve import LightCurve
        from .pipeline import Pipeline
        import pandas as pd
        carry = dict(self.lc.scalars)
        for k, v in self.vars.items():
            if k != "Name" and k not in carry:
                carry[k] = v
        input_lc = LightCurve(self.lc._df, name=self.lc.name, scalars=carry)
        # Chain offset: the new segment's output-column suffixes should
        # continue the numbering from where the prior segment left off so
        # they don't collide with the injected scalars.  Pipeline.run uses
        # this to emit explicit -columnsuffix values, so the returned
        # Result already carries correctly-shifted suffixes.
        prior_known = list(self._known_commands or [])
        offset = len(prior_known)
        new_result = Pipeline([cmd]).run(
            input_lc, capture_lc=capture_lc,
            _command_offset=offset, **run_opts)
        # Concatenate prior vars with the new segment's (already-shifted)
        # vars.  No position shift here because the C side emitted the
        # suffixes we want.  Drop the prior "Name" key so the new one
        # (authoritative for this segment) wins.
        prior_vars = self.vars.drop("Name", errors="ignore")
        merged_vars = pd.concat([prior_vars, new_result.vars])
        merged_known = prior_known + list(new_result._known_commands or [])
        from .results import Result
        return Result(
            var=merged_vars,
            lc=new_result.lc,
            files=new_result.files,
            known_commands=merged_known,
            error=new_result.error,
        )

    _immediate.__name__ = vt_name
    _immediate.__qualname__ = f"Result.{vt_name}"
    _immediate.__doc__ = (
        f"Run -{vt_name} on this result's light curve and return a new Result.\n\n"
        f"Output variables from all prior commands are preserved in the returned "
        f"result's ``.vars`` and ``.varobjs``.\n"
        f"Run-time options (capture_lc, randseed, timeout, etc.) may be "
        f"passed as keyword arguments.\n\n"
        + (inspect.cleandoc(cmd_cls.__doc__) if cmd_cls.__doc__ else "")
    )
    return _immediate


def _make_stateful_error(cmd_cls: Type["VartoolsCommand"], owner: str):
    """Generate a method that raises NotImplementedError for pipeline-stateful
    commands when called on LightCurve or Result."""
    vt_name = cmd_cls._vt_name

    def _error(self, *args, **kwargs):
        raise NotImplementedError(
            f"-{vt_name} is a pipeline-stateful command that must be used "
            f"within a single vartools invocation.  Call it via "
            f"Pipeline([..., commands.{vt_name}(...), ...]).run(lc) instead."
        )

    _error.__name__ = vt_name
    _error.__qualname__ = f"{owner}.{vt_name}"
    _error.__doc__ = (
        f"-{vt_name} is a pipeline-stateful command.\n\n"
        f"It cannot be called directly on a {owner} because it only works "
        f"correctly within a single vartools invocation.  Use "
        f"``Pipeline().{vt_name}(...).run(lc)`` instead."
    )
    return _error


# ---------------------------------------------------------------------------
# Method factory — Pipeline builder
# ---------------------------------------------------------------------------

def _make_pipeline_builder(cmd_cls: Type["VartoolsCommand"]):
    """Builder method for Pipeline: pipe.LS(...) appends a command and returns self.

    Lets users construct a Pipeline with chained method calls::

        pipe = vt.Pipeline().clip(5.0).LS(0.5, 10.0, 1e-3)
        result = pipe.run(lc)

    This is equivalent to passing a list of command instances to the
    Pipeline constructor::

        pipe = vt.Pipeline([cmd.clip(5.0), cmd.LS(0.5, 10.0, 1e-3)])
    """
    vt_name = cmd_cls._vt_name

    def _builder(self, *args, **kwargs):
        self.commands.append(cmd_cls(*args, **kwargs))
        return self

    _builder.__name__ = vt_name
    _builder.__qualname__ = f"Pipeline.{vt_name}"
    _builder.__doc__ = (
        f"Append a -{vt_name} command to the pipeline and return ``self``.\n\n"
        f"Positional and keyword arguments are forwarded verbatim to "
        f"``commands.{vt_name}.__init__``.  The returned Pipeline is the "
        f"same instance the method was called on, so calls can be chained.\n\n"
        + (inspect.cleandoc(cmd_cls.__doc__) if cmd_cls.__doc__ else "")
    )
    return _builder


# ---------------------------------------------------------------------------
# Method factories — LightCurveBatch and BatchResult (deferred + immediate)
# ---------------------------------------------------------------------------

def _make_deferred_batch(cmd_cls: Type["VartoolsCommand"]):
    """Deferred method for LightCurveBatch: batch.LS(...) → LightCurveBatch."""
    vt_name = cmd_cls._vt_name

    def _deferred(self, *args, **kwargs):
        cmd = cmd_cls(*args, **kwargs)
        return self._append(cmd)

    _deferred.__name__ = vt_name
    _deferred.__qualname__ = f"LightCurveBatch.{vt_name}"
    _deferred.__doc__ = (
        f"Append -{vt_name} to the batch command chain.\n\n"
        f"Returns a new LightCurveBatch.  Call .run() to execute on all LCs.\n\n"
        + (inspect.cleandoc(cmd_cls.__doc__) if cmd_cls.__doc__ else "")
    )
    return _deferred


def _make_immediate_batch(cmd_cls: Type["VartoolsCommand"]):
    """Immediate method for LightCurveBatch: batch.run_LS(...) → BatchResult."""
    vt_name = cmd_cls._vt_name
    method_name = f"run_{vt_name}"

    def _immediate(self, *args, **kwargs):
        run_opts, cmd_kwargs = _split_kwargs(kwargs)
        capture_lc = run_opts.pop("capture_lc", True)
        cmd = cmd_cls(*args, **cmd_kwargs)
        from ._batch import LightCurveBatch
        one_cmd_batch = LightCurveBatch(
            self._lcs,
            _commands=[cmd],
            _global_opts=self._global_opts,
        )
        return one_cmd_batch.run(capture_lc=capture_lc, **run_opts)

    _immediate.__name__ = method_name
    _immediate.__qualname__ = f"LightCurveBatch.{method_name}"
    _immediate.__doc__ = (
        f"Run -{vt_name} immediately on all LCs and return a BatchResult.\n\n"
        f"Run-time options (capture_lc, randseed, timeout, etc.) may be "
        f"passed as keyword arguments alongside the command parameters.\n\n"
        + (inspect.cleandoc(cmd_cls.__doc__) if cmd_cls.__doc__ else "")
    )
    return _immediate


def _make_deferred_batch_result(cmd_cls: Type["VartoolsCommand"]):
    """Deferred method for BatchResult: batch.LS(...) → LightCurveBatch."""
    vt_name = cmd_cls._vt_name

    def _deferred(self, *args, **kwargs):
        if self._lcs is None:
            raise AttributeError(
                f"Cannot chain -{vt_name} from BatchResult: light curves were not "
                f"captured (capture_lc=False).  Re-run with capture_lc=True "
                f"(the default for LightCurveBatch) to enable further chaining."
            )
        from ._batch import LightCurveBatch
        cmd = cmd_cls(*args, **kwargs)
        # Resolve any back-reference keywords against the prior BatchResult.
        # Per-LC outputs come back as PerLC arrays so the batch-mode
        # ``-inlistvars`` plumbing can inject one value per light curve.
        cmd._resolve_back_references(self)
        return LightCurveBatch(self.lcs, _commands=[cmd], _prior_batch=self)

    _deferred.__name__ = vt_name
    _deferred.__qualname__ = f"BatchResult.{vt_name}"
    _deferred.__doc__ = (
        f"Continue chaining from this batch result's light curves with -{vt_name}.\n\n"
        f"Returns a LightCurveBatch.  Call .run() to execute on all LCs.\n"
        f"Requires that this BatchResult was produced with capture_lc=True.\n\n"
        + (inspect.cleandoc(cmd_cls.__doc__) if cmd_cls.__doc__ else "")
    )
    return _deferred


def _make_immediate_batch_result(cmd_cls: Type["VartoolsCommand"]):
    """Immediate method for BatchResult: batch.run_LS(...) → BatchResult."""
    vt_name = cmd_cls._vt_name
    method_name = f"run_{vt_name}"

    def _immediate(self, *args, **kwargs):
        if self._lcs is None:
            raise AttributeError(
                f"Cannot call {method_name} from BatchResult: light curves were "
                f"not captured (capture_lc=False).  Re-run with capture_lc=True "
                f"to enable further chaining."
            )
        run_opts, cmd_kwargs = _split_kwargs(kwargs)
        capture_lc = run_opts.pop("capture_lc", True)
        from ._batch import LightCurveBatch
        cmd = cmd_cls(*args, **cmd_kwargs)
        return LightCurveBatch(self.lcs, _commands=[cmd], _prior_batch=self).run(
            capture_lc=capture_lc, **run_opts
        )

    _immediate.__name__ = method_name
    _immediate.__qualname__ = f"BatchResult.{method_name}"
    _immediate.__doc__ = (
        f"Run -{vt_name} immediately on all LCs in this batch result.\n\n"
        f"Returns a new BatchResult.  Requires capture_lc=True on the "
        f"preceding run.\n\n"
        + (inspect.cleandoc(cmd_cls.__doc__) if cmd_cls.__doc__ else "")
    )
    return _immediate


# ---------------------------------------------------------------------------
# Public attachment function
# ---------------------------------------------------------------------------

def _attach_all_command_methods() -> None:
    """Generate and attach command methods to LightCurve, Result, Pipeline,
    LightCurveBatch, and BatchResult.

    Called once at the bottom of pyvartools/__init__.py.
    """
    from . import commands as _cmds
    from .lightcurve import LightCurve
    from ._batch import LightCurveBatch
    from .pipeline import Pipeline
    from .results import Result, BatchResult

    _SKIP = {"Raw", "UserCommand"}

    cmd_classes: List[Type["VartoolsCommand"]] = [
        cls for name in _cmds.__all__
        if name not in _SKIP and name != "VartoolsCommand"
        for cls in [getattr(_cmds, name)]
        if isinstance(cls, type) and issubclass(cls, _cmds.VartoolsCommand)
    ]

    for cmd_cls in cmd_classes:
        vt_name = cmd_cls._vt_name
        if not vt_name:
            continue

        run_name = f"run_{vt_name}"

        if vt_name in _STATEFUL:
            # Pipeline-stateful commands: raise a helpful error on LightCurve
            # and Result; work normally in LightCurveBatch (runs full pipeline).
            if not hasattr(LightCurve, vt_name):
                setattr(LightCurve, vt_name, _make_stateful_error(cmd_cls, "LightCurve"))
            if not hasattr(Result, vt_name):
                setattr(Result, vt_name, _make_stateful_error(cmd_cls, "Result"))
        else:
            # LightCurve: single immediate method
            if not hasattr(LightCurve, vt_name):
                setattr(LightCurve, vt_name, _make_immediate_lc(cmd_cls))

            # Result: single immediate method (merges prior vars)
            if not hasattr(Result, vt_name):
                setattr(Result, vt_name, _make_immediate_result(cmd_cls))

        # Pipeline: builder method that appends the command and returns self.
        # Use the class name (what callers write as ``cmd.X``) as the primary
        # method name; also alias the vt_name when it differs and is a valid
        # Python identifier.  Raise if a name would shadow an existing
        # Pipeline attribute (run, run_batch, commands, …) — silent skips
        # would hide a real API design mistake.
        cls_name = cmd_cls.__name__
        builder_names = {cls_name}
        if vt_name != cls_name and vt_name.isidentifier():
            builder_names.add(vt_name)
        for nm in builder_names:
            if hasattr(Pipeline, nm):
                existing = getattr(Pipeline, nm)
                if getattr(existing, "__qualname__", "") == f"Pipeline.{nm}":
                    continue  # already attached
                raise RuntimeError(
                    f"Cannot attach Pipeline.{nm} builder for {cls_name}: "
                    f"the name collides with an existing Pipeline attribute "
                    f"({existing!r}).  Rename the command or the attribute."
                )
            method = _make_pipeline_builder(cmd_cls)
            method.__name__ = nm
            method.__qualname__ = f"Pipeline.{nm}"
            setattr(Pipeline, nm, method)

        # LightCurveBatch: deferred chain-builder + immediate run_*
        if not hasattr(LightCurveBatch, vt_name):
            setattr(LightCurveBatch, vt_name, _make_deferred_batch(cmd_cls))
        if not hasattr(LightCurveBatch, run_name):
            setattr(LightCurveBatch, run_name, _make_immediate_batch(cmd_cls))

        # BatchResult: deferred chain-builder + immediate run_*
        if not hasattr(BatchResult, vt_name):
            setattr(BatchResult, vt_name, _make_deferred_batch_result(cmd_cls))
        if not hasattr(BatchResult, run_name):
            setattr(BatchResult, run_name, _make_immediate_batch_result(cmd_cls))
