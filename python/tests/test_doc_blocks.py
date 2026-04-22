"""Regression tests for executable code blocks in the website docs.

Runs every ```python block in the Python-API pages through a helper harness
(see ``run_doc_blocks.py``).  Each page is one parametrized test; blocks
inside a page share a subprocess namespace so the cost of importing
pyvartools is paid once per page.

Set ``PYVARTOOLS_DOCS_ROOT`` to the docs directory (defaults to
``$HOME/src/vartools-site/docs/python``).  If that path doesn't exist — e.g.
in CI, or for a developer without the site checkout — every parametrized
test is skipped rather than failed.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

_HARNESS_PATH = Path(__file__).with_name("run_doc_blocks.py")
_spec = importlib.util.spec_from_file_location(
    "_pyvartools_docblocks_harness", _HARNESS_PATH
)
harness = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = harness
_spec.loader.exec_module(harness)


_DOCS_ROOT = harness._DEFAULT_DOCS_ROOT
_CWD = harness._DEFAULT_CWD
_SITE_AVAILABLE = _DOCS_ROOT.is_dir()


@pytest.mark.skipif(
    not _SITE_AVAILABLE,
    reason=f"Website docs not found at {_DOCS_ROOT}; "
           f"set PYVARTOOLS_DOCS_ROOT to enable.",
)
@pytest.mark.parametrize("page", harness.PAGES)
def test_doc_blocks(page):
    results = harness.run_page(page, _DOCS_ROOT, _CWD)
    failures = [r for r in results if r[1] != "PASS"]
    if failures:
        lines = [f"{len(failures)} failing block(s) in {page}:"]
        for label, status, detail in failures:
            lines.append(f"  {label}")
            for line in detail.splitlines():
                lines.append(f"    {line}")
        pytest.fail("\n".join(lines))
