"""Utility functions for pyvartools."""

import subprocess
from typing import List

from ._binary import get_binary


def list_commands() -> List[str]:
    """Return the list of commands supported by the installed vartools binary.

    Calls ``vartools -listcommands``.  vartools exits with a non-zero status
    for this informational command, which is normal — the exit code is ignored.

    Returns
    -------
    list of str
        Command names, one per entry, in the order reported by vartools.
    """
    binary = get_binary()
    try:
        proc = subprocess.run(
            [binary, "-listcommands"],
            capture_output=True,
            text=True,
        )
        # vartools -listcommands exits non-zero but still prints to stdout
        output = proc.stdout or proc.stderr
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Cannot find vartools binary: {e}") from e

    commands = []
    for line in output.splitlines():
        # Command names appear as "-name" at the start of a line (no leading space)
        if line.startswith("-") and not line.startswith("- "):
            name = line.split()[0].lstrip("-")
            if name:
                commands.append(name)
    return commands
