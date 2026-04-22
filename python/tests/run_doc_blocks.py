#!/usr/bin/env python3
"""Run every ```python block in the pyvartools website docs and report pass/fail.

Paths (overridable):

    PYVARTOOLS_DOCS_ROOT  directory holding the Python-API markdown pages
                          (default: $HOME/src/vartools-site/docs/python)
    PYVARTOOLS_DOCS_CWD   working directory the code blocks execute in
                          (default: repository root — two levels up from this file)

Both can also be passed as CLI args (--docs-root, --cwd).  The script exits
non-zero if any block raises, 0 otherwise.  When DOCS_ROOT is missing the
script prints a skip message and exits 0 so CI can call it unconditionally.

- Blocks within a page share a single namespace (sequential execution).
- Each block is wrapped in try/except so one failure doesn't mask later blocks.
- Signature-only blocks (`cmd.X(args, ...)` with no assignment) are skipped.
"""
from __future__ import annotations

import argparse
import ast
import json
import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path
from typing import List, Tuple

_REPO_ROOT = Path(__file__).resolve().parents[2]
_DEFAULT_DOCS_ROOT = Path(
    os.environ.get(
        "PYVARTOOLS_DOCS_ROOT",
        str(Path.home() / "src" / "vartools-site" / "docs" / "python"),
    )
)
_DEFAULT_CWD = Path(os.environ.get("PYVARTOOLS_DOCS_CWD", str(_REPO_ROOT)))

PAGES = [
    "index.md",
    "lightcurve.md",
    "fluent.md",
    "pipeline.md",
    "results.md",
    "commands/index.md",
    "commands/manipulation.md",
    "commands/period-finding.md",
    "commands/statistics.md",
    "commands/model-fitting.md",
    "commands/filtering.md",
    "commands/simulation.md",
    "commands/control-flow.md",
    "commands/misc.md",
    "commands/python-r.md",
    "commands/extensions.md",
]

_PREAMBLE_TEMPLATE = textwrap.dedent("""
    import os, sys
    os.chdir(%r)
    import warnings, traceback, json
    warnings.filterwarnings("ignore")
    import numpy as np
    import pandas as pd
    import pyvartools as vt
    from pyvartools import commands as cmd
    from pyvartools import (
        Pipeline, PerLC, LightCurve, LightCurveBatch, Output, Result,
        BatchResult, LCVar, ListVar,
    )
    _RESULTS = []
    _NS = globals()
    def _run(label, code):
        try:
            exec(code, _NS)
        except Exception as exc:
            tb = traceback.format_exc().splitlines()[-3:]
            _RESULTS.append((label, "FAIL", "\\n".join(tb)))
            return
        _RESULTS.append((label, "PASS", ""))
""")

_FOOTER = "\nprint('__BLOCK_RESULTS__', json.dumps(_RESULTS))\n"
_FENCE_RE = re.compile(r"```(python|py)\s*\n(.*?)```", re.DOTALL)


def extract_blocks(md_path: Path) -> List[Tuple[int, str]]:
    text = md_path.read_text()
    out = []
    for m in _FENCE_RE.finditer(text):
        start_line = text[: m.start()].count("\n") + 1
        out.append((start_line, textwrap.dedent(m.group(2))))
    return out


def is_signature_only(src: str) -> bool:
    src = src.strip()
    if not src:
        return True
    if re.search(r'(^|[,\s(])\.\.\.([,\s)]|$)', src):
        return True
    try:
        tree = ast.parse(src)
    except SyntaxError:
        if re.match(r'^\s*(cmd|vt|lc|result|batch|pipe)\.\w+\(', src):
            return True
        return False

    def _is_cmd_call(node):
        if not isinstance(node, ast.Expr):
            return False
        call = node.value
        if not isinstance(call, ast.Call):
            return False
        return (isinstance(call.func, ast.Attribute)
                and isinstance(call.func.value, ast.Name)
                and call.func.value.id == "cmd")
    if tree.body and all(_is_cmd_call(n) for n in tree.body):
        return True
    return False


def _build_harness(rel: str, blocks, cwd: Path) -> str:
    parts = [_PREAMBLE_TEMPLATE % str(cwd)]
    for ln, src in blocks:
        label = f"{rel}:{ln}"
        parts.append(f"_run({label!r}, {src!r})")
    parts.append(_FOOTER)
    return "\n".join(parts)


def run_page(rel: str, docs_root: Path, cwd: Path) -> list:
    path = docs_root / rel
    if not path.exists():
        return []
    blocks = [(ln, src) for ln, src in extract_blocks(path)
              if not is_signature_only(src)]
    if not blocks:
        return []
    harness = _build_harness(rel, blocks, cwd)
    proc = subprocess.run(
        [sys.executable, "-c", harness],
        capture_output=True, text=True, timeout=900,
    )
    out = proc.stdout
    marker_idx = out.rfind("__BLOCK_RESULTS__")
    if marker_idx < 0:
        return [("HARNESS_ERROR", "FAIL", proc.stderr[-2000:])]
    line = out[marker_idx:].splitlines()[0]
    try:
        data = json.loads(line.split(" ", 1)[1])
    except Exception as exc:
        return [("HARNESS_JSON", "FAIL", str(exc) + "\n" + line[:2000])]
    return data


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--docs-root", type=Path, default=_DEFAULT_DOCS_ROOT,
                    help="Directory with the Python-API docs (.md pages).")
    ap.add_argument("--cwd", type=Path, default=_DEFAULT_CWD,
                    help="Working directory used when executing code blocks.")
    args = ap.parse_args()

    if not args.docs_root.is_dir():
        print(f"skip: docs root {args.docs_root} does not exist — "
              f"set PYVARTOOLS_DOCS_ROOT or pass --docs-root.")
        return 0

    total_fail = 0
    for rel in PAGES:
        rs = run_page(rel, args.docs_root, args.cwd)
        fails = [r for r in rs if r[1] != "PASS"]
        total_fail += len(fails)
        print(f"=== {rel}: {len(rs)} blocks, {len(fails)} failing ===")
        for label, status, detail in fails:
            print(f"  FAIL  {label}")
            for line in detail.splitlines():
                print(f"        {line}")
    print(f"\nTOTAL FAILING BLOCKS: {total_fail}")
    return 0 if total_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
