"""
pyvartools — Python API for the VARTOOLS light curve analysis program.

Quick start
-----------
    import pyvartools as vt

    # Configure binary location if not on PATH
    # vt.config.set_binary("/path/to/vartools")

    # Load a light curve
    lc = vt.LightCurve.from_file("EXAMPLES/2")

    # Run a single command
    result = vt.Pipeline([vt.commands.LS(0.1, 10.0, 0.1, 5, 1)]).run(lc)
    print(result.stats)

    # Chain commands
    pipe = vt.Pipeline([
        vt.commands.clip(sigclip=5.0),
        vt.commands.LS(0.1, 10.0, 0.1, 5, 1),
    ])
    result = pipe.run(lc)

    # Batch processing
    lc_list = [vt.LightCurve.from_file(f"EXAMPLES/{i}") for i in range(1, 11)]
    batch = pipe.run_batch(lc_list, nthreads=4)
    print(batch.stats)
"""

from pyvartools import commands
from pyvartools._binary import get_binary, set_binary, set_library
from pyvartools._output import Output
from pyvartools.lightcurve import LightCurve
from pyvartools.perlc import PerLC
from pyvartools.pipeline import LCVar, ListVar, Pipeline
from pyvartools.results import BatchResult, Result, RunError
from pyvartools._utils import list_commands
from pyvartools.userlib import UserCommand, discover_userlibs, load_userlib

# Expose config.set_binary as vt.config.set_binary
from pyvartools import _binary as config

__all__ = [
    "LightCurve",
    "LCVar",
    "ListVar",
    "Output",
    "PerLC",
    "Pipeline",
    "Result",
    "BatchResult",
    "RunError",
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
