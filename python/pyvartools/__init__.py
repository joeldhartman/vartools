"""
pyvartools — Python API for the VARTOOLS light curve analysis program.

Quick start
-----------
    import pyvartools as vt

    # Configure binary location if not on PATH
    # vt.config.set_binary("/path/to/vartools")

    # Simplest usage — pass a filename (or DataFrame / LightCurve / numpy array)
    result = vt.LS("EXAMPLES/2", 1.0, 2.0, 0.01)
    print(result.varobjs.LS.Period_1)

    # Fluent chaining on a LightCurve object
    lc = vt.LightCurve.from_file("EXAMPLES/2")
    result = lc.LS(1.0, 2.0, 0.01).rms().run()

    # Batch processing
    lc_list = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 11)]
    batch = vt.LightCurveBatch(lc_list).LS(1.0, 2.0, 0.01).run()
    print(batch.vars)

    # Pipeline mode for maximum throughput (single vartools invocation)
    from pyvartools import commands as cmd
    pipe = vt.Pipeline([cmd.clip(sigclip=5.0), cmd.LS(0.1, 10.0, 0.1)])
    batch = pipe.run_batch(lc_list, nthreads=4)
    print(batch.vars)
"""

from pyvartools import commands
from pyvartools._binary import get_binary, set_binary, set_library
from pyvartools._output import Output
from pyvartools._vars import CommandStats, VarsNamespace
from pyvartools._batch import LightCurveBatch
from pyvartools.lightcurve import LightCurve
from pyvartools.perlc import PerLC
from pyvartools.pipeline import LCColumn, LCVar, ListVar, Pipeline
from pyvartools.results import BatchResult, PipelineValidationError, Result, RunError
from pyvartools._utils import list_commands
from pyvartools.userlib import UserCommand, discover_userlibs, load_userlib

# Expose config.set_binary as vt.config.set_binary
from pyvartools import _binary as config

__all__ = [
    "LightCurve",
    "LCColumn",
    "LCVar",
    "ListVar",
    "Output",
    "PerLC",
    "Pipeline",
    "LightCurveBatch",
    "Result",
    "BatchResult",
    "PipelineValidationError",
    "RunError",
    "CommandStats",
    "VarsNamespace",
    "commands",
    "config",
    "get_binary",
    "set_binary",
    "set_library",
    "list_commands",
    "UserCommand",
    "load_userlib",
    "discover_userlibs",
]

# Attach fluent command methods to LightCurve, CommandChain, and Result.
# This runs once at import time and generates ~110 named methods.
from pyvartools._method_gen import _attach_all_command_methods
_attach_all_command_methods()

# Attach top-level command functions (vt.LS, vt.rms, ...) to this module.
# Each accepts a light curve as its first argument (LightCurve, filename,
# DataFrame, 2-D numpy array, or tuple of arrays).
import sys as _sys
from pyvartools._toplevel import _attach_toplevel_commands as _atc
_toplevel_names = _atc(_sys.modules[__name__])
__all__ = __all__ + _toplevel_names
del _sys, _atc, _toplevel_names
