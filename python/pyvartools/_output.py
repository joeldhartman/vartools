"""Output descriptor for auxiliary file handling in pyvartools commands."""

from __future__ import annotations

from typing import Optional


class Output:
    """Controls how a command's auxiliary output file is handled.

    Parameters
    ----------
    path : str or None
        Directory to write the file to on disk.  ``None`` (default) means use
        a pipeline-managed temporary directory; the file is deleted after the
        run completes (unless ``capture=True``, in which case it is read into
        Python first).
    capture : bool
        Whether to read the file into Python and include it in
        ``result.files``.  Default ``True``.

    Shorthand forms accepted by all ``save_*`` parameters
    -------------------------------------------------------
    ``False``
        Mode 4 — suppress (don't write at all).
    ``True``
        Mode 1 — temp dir, capture into Python (default when save enabled).
    ``"/path/to/dir"``
        Mode 3 — write to that directory on disk, no Python capture.
    ``Output("/path/to/dir", capture=True)``
        Mode 2 — write to disk AND capture into ``result.files``.

    Notes
    -----
    For ``autocorrelation``, Mode 4 (``save_result=False``) still writes the
    file to disk (the CLI always does so), but the file is not captured into
    Python.
    """

    def __init__(self, path: Optional[str] = None, capture: bool = True) -> None:
        self.path = path
        self.capture = capture

    def __repr__(self) -> str:
        return f"Output(path={self.path!r}, capture={self.capture!r})"
