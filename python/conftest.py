# pytest configuration for pyvartools tests.
# Placed here so pytest finds it when run from the python/ directory.

import os
import shutil

collect_ignore_glob = []

# Allow tests to be run from the vartools source root as well:
#   pytest python/tests/


def pytest_collection_finish(session):
    """Pin the freshest vartools binary for subprocess-mode tests.

    Several test modules call ``vt.set_binary(<in-tree src/vartools>)`` at
    import time.  In a working tree that is only *mirrored* and never built
    (e.g. the GIT mirror of the SVN source, where only the SVN tree is
    compiled and installed) that in-tree binary is stale and lacks the
    newest commands, which makes the subprocess-based ``validate()`` sweep
    reject them.  This hook runs once after collection -- i.e. after every
    test module's import-time ``set_binary`` -- and prefers the installed
    binary found on ``PATH`` (kept current per the install discipline).  In
    CI, where vartools is built in-tree but not installed on ``PATH``,
    ``shutil.which`` returns ``None`` and the in-tree binary (freshly built
    there) is left in place.  An explicit ``VARTOOLS_BINARY`` override is
    always respected.
    """
    if os.environ.get("VARTOOLS_BINARY"):
        return
    path_binary = shutil.which("vartools")
    if not path_binary:
        return
    try:
        import pyvartools as vt
        vt.set_binary(path_binary)
    except Exception:
        # Never let binary pinning abort the test session; modules that
        # need a binary already gate on their own _HAVE_BINARY checks.
        pass
