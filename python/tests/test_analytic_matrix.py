"""Comprehensive analytic-engine regression matrix.

Skipped by default.  To run::

    PYVARTOOLS_RUN_ANALYTIC_MATRIX=1 pytest tests/test_analytic_matrix.py -v

Run this whenever something in ``src/analytic.c`` changes -- it
exercises every built-in scalar / aggregate function, deeply nested
expressions, and constant-fold cases, in three invocation modes
(``-i`` single-LC, ``-l`` serial, ``-l -parallel 4``) on five
distinct LCs.  Each test verifies:

  1. **Mode consistency** -- the per-LC value matches across all
     three invocations (catches memo cache leaking across LCs in
     batch mode, parallel-thread cache contention, etc.).
  2. **Numerical agreement with numpy** (where the convention
     matches; MAD / kurtosis / skewness / pct / round use different
     bias / normalisation conventions, so for those we check
     mode-consistency only).

The slow part is invoking vartools as a subprocess three times per
test case * ~50 cases * 5 LCs ~= a few minutes, which is why this is
not part of the default ``pytest`` run.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile

import numpy as np
import pytest

# ----------------------------- skip gating -----------------------------

_RUN_FLAG = "PYVARTOOLS_RUN_ANALYTIC_MATRIX"
_SHOULD_RUN = os.environ.get(_RUN_FLAG) == "1"

try:
    _VARTOOLS = subprocess.check_output(
        ["which", "vartools"], stderr=subprocess.DEVNULL
    ).decode().strip()
    _HAVE_BIN = bool(_VARTOOLS)
except Exception:
    _VARTOOLS = ""
    _HAVE_BIN = False

pytestmark = [
    pytest.mark.skipif(
        not _SHOULD_RUN,
        reason=f"set {_RUN_FLAG}=1 to run the analytic-engine matrix.",
    ),
    pytest.mark.skipif(
        not _HAVE_BIN,
        reason="vartools binary not on PATH",
    ),
]

# ----------------------------- fixtures --------------------------------

_N_LC = 5
_NPTS = 200


def _make_lc(seed):
    """Deterministic LC.  mag in [10.0, 10.5] avoids log/sqrt traps."""
    rng = np.random.default_rng(seed)
    t = np.linspace(2459000.0, 2459000.5, _NPTS)
    mag = 10.0 + 0.5 * rng.random(_NPTS)
    err = 0.001 + 0.0001 * rng.random(_NPTS)
    return t, mag, err


@pytest.fixture(scope="module")
def lc_fixture():
    tmp = tempfile.mkdtemp(prefix="vt_amat_")
    lc_paths = []
    lc_data = []
    for i in range(_N_LC):
        t, mag, err = _make_lc(seed=42 + i)
        p = os.path.join(tmp, f"lc{i:03d}.dat")
        with open(p, "w") as f:
            for ti, mi, ei in zip(t, mag, err):
                f.write(f"{ti:.10f} {mi:.10f} {ei:.10f}\n")
        lc_paths.append(p)
        lc_data.append((t, mag, err))
    list_path = os.path.join(tmp, "list")
    with open(list_path, "w") as f:
        for p in lc_paths:
            f.write(p + "\n")
    yield {"tmp": tmp, "lc_paths": lc_paths, "lc_data": lc_data,
           "list_path": list_path}
    shutil.rmtree(tmp, ignore_errors=True)


# ----------------------------- helpers ---------------------------------

def _stats_for_expr(expr, mode, *, lc_path=None, list_path=None, parallel=0):
    """Run vartools -expr a=expr -stats a median.  Returns
    {lc_name: median_value}."""
    if mode == "i":
        cmd = [_VARTOOLS, "-i", lc_path,
               "-expr", f"a=({expr})", "-stats", "a", "median"]
    else:
        cmd = [_VARTOOLS, "-l", list_path]
        if parallel > 0:
            cmd += ["-parallel", str(parallel)]
        cmd += ["-expr", f"a=({expr})", "-stats", "a", "median"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, (
        f"vartools failed for `{expr}` in mode={mode}:\n{r.stderr}"
    )
    out = {}
    for line in r.stdout.strip().split("\n"):
        toks = line.split()
        if len(toks) >= 2:
            try:
                out[toks[0]] = float(toks[1])
            except ValueError:
                pass
    return out


# ----------------------------- test matrix -----------------------------

def _gen_cases():
    """Yield (id, label, expr, kind, expected_fn) -- expected_fn is None
    for `consistency_only` cases."""
    cases = [
        ("exp",      "exp(mag)",          "perpoint", lambda m, e: np.exp(m)),
        ("log",      "log(mag)",          "perpoint", lambda m, e: np.log(m)),
        ("log10",    "log10(mag)",        "perpoint", lambda m, e: np.log10(m)),
        ("sqrt",     "sqrt(mag)",         "perpoint", lambda m, e: np.sqrt(m)),
        ("abs",      "abs(mag-10.25)",    "perpoint", lambda m, e: np.abs(m - 10.25)),
        ("max2",     "max(mag,10.3)",     "perpoint", lambda m, e: np.maximum(m, 10.3)),
        ("min2",     "min(mag,10.3)",     "perpoint", lambda m, e: np.minimum(m, 10.3)),
        ("hypot",    "hypot(mag,err)",    "perpoint", lambda m, e: np.hypot(m, e)),
        ("sin",      "sin(mag)",          "perpoint", lambda m, e: np.sin(m)),
        ("cos",      "cos(mag)",          "perpoint", lambda m, e: np.cos(m)),
        ("tan",      "tan(mag*0.1)",      "perpoint", lambda m, e: np.tan(m * 0.1)),
        ("asin",     "asin(mag*0.05)",    "perpoint", lambda m, e: np.arcsin(m * 0.05)),
        ("acos",     "acos(mag*0.05)",    "perpoint", lambda m, e: np.arccos(m * 0.05)),
        ("atan2",    "atan2(mag,err)",    "perpoint", lambda m, e: np.arctan2(m, e)),
        ("ceil",     "ceil(mag)",         "perpoint", lambda m, e: np.ceil(m)),
        ("floor",    "floor(mag)",        "perpoint", lambda m, e: np.floor(m)),
        ("cosh",     "cosh(mag*0.1)",     "perpoint", lambda m, e: np.cosh(m * 0.1)),
        ("sinh",     "sinh(mag*0.1)",     "perpoint", lambda m, e: np.sinh(m * 0.1)),
        ("tanh",     "tanh(mag*0.1)",     "perpoint", lambda m, e: np.tanh(m * 0.1)),
        ("acosh",    "acosh(mag)",        "perpoint", lambda m, e: np.arccosh(m)),
        ("asinh",    "asinh(mag*0.1)",    "perpoint", lambda m, e: np.arcsinh(m * 0.1)),
        ("atanh",    "atanh(mag*0.05)",   "perpoint", lambda m, e: np.arctanh(m * 0.05)),
        ("theta",    "theta(mag-10.25)",  "perpoint", lambda m, e: (m >= 10.25).astype(float)),
        ("isnan",    "isnan(mag)",        "perpoint", lambda m, e: np.zeros_like(m)),
        # Aggregates (numpy match)
        ("mean",         "mean(mag)",                "scalar",
                          lambda m, e: np.mean(m)),
        ("median",       "median(mag)",              "scalar",
                          lambda m, e: np.median(m)),
        ("stddev",       "stddev(mag)",              "scalar",
                          lambda m, e: np.std(m, ddof=1)),
        ("sum",          "sum(mag)",                 "scalar",
                          lambda m, e: np.sum(m)),
        ("vmin",         "vmin(mag)",                "scalar",
                          lambda m, e: np.min(m)),
        ("vmax",         "vmax(mag)",                "scalar",
                          lambda m, e: np.max(m)),
        ("weightedmean", "weightedmean(mag,err)",    "scalar",
                          lambda m, e: np.sum(m / e**2) / np.sum(1.0 / e**2)),
        ("len",          "len(mag)",                 "scalar",
                          lambda m, e: float(len(m))),
        # Convention-dependent aggregates: mode-consistency only.
        ("MAD",       "MAD(mag)",       "consistency_only", None),
        ("kurtosis",  "kurtosis(mag)",  "consistency_only", None),
        ("skewness",  "skewness(mag)",  "consistency_only", None),
        ("pct",       "pct(mag,25)",    "consistency_only", None),
        ("round",     "round(mag*10)",  "consistency_only", None),
        # Aggregate-in-per-point (memo hot path)
        ("residual",     "mag - mean(mag)",                          "perpoint",
                          lambda m, e: m - np.mean(m)),
        ("z-score",      "(mag-mean(mag))/stddev(mag)",              "perpoint",
                          lambda m, e: (m - np.mean(m)) / np.std(m, ddof=1)),
        ("chi2-like",    "((mag-mean(mag))^2)/(err^2)",              "perpoint",
                          lambda m, e: ((m - np.mean(m)) ** 2) / (e ** 2)),
        ("3-deep nest",  "sqrt(abs(mag-median(mag-mean(mag))))",     "perpoint",
                          lambda m, e: np.sqrt(np.abs(m - np.median(m - np.mean(m))))),
        ("4-deep nest",  "log(1+abs(sqrt(abs(mag-mean(mag)))))",     "perpoint",
                          lambda m, e: np.log(1 + np.abs(np.sqrt(np.abs(m - np.mean(m)))))),
        # Constant folds
        ("fold-sqrt",    "sqrt(4.0)+0",                               "scalar",
                          lambda m, e: 2.0),
        ("fold-exp+sin", "exp(0)+sin(0)",                             "scalar",
                          lambda m, e: 1.0),
        ("fold-eq",      "(61==0) + (1==1)",                          "scalar",
                          lambda m, e: 1.0),
        ("fold-abs+log", "abs(-3.5)+log(1)",                          "scalar",
                          lambda m, e: 3.5),
        ("fold-pi",      "round(2.7)*pi/pi",                          "scalar",
                          lambda m, e: 3.0),
    ]
    # Optional scipy-backed cases
    try:
        from scipy.special import erf, erfc, gammaln, gamma as gamma_fn
        cases += [
            ("erf",    "erf(mag*0.1)",    "perpoint", lambda m, e: erf(m * 0.1)),
            ("erfc",   "erfc(mag*0.1)",   "perpoint", lambda m, e: erfc(m * 0.1)),
            ("lgamma", "lgamma(mag)",     "perpoint", lambda m, e: gammaln(m)),
            ("gamma",  "gamma(mag*0.5)",  "perpoint", lambda m, e: gamma_fn(m * 0.5)),
        ]
    except ImportError:
        pass
    return cases


@pytest.mark.parametrize("label,expr,kind,fn", _gen_cases(),
                         ids=[c[0] for c in _gen_cases()])
def test_analytic_engine(lc_fixture, label, expr, kind, fn):
    lc_paths = lc_fixture["lc_paths"]
    lc_data = lc_fixture["lc_data"]
    list_path = lc_fixture["list_path"]

    # 1. Single-LC reference: run -i per LC, build name->value map.
    ref_i = {}
    for p in lc_paths:
        ref_i.update(_stats_for_expr(expr, "i", lc_path=p))
    # 2. List serial.
    ref_l = _stats_for_expr(expr, "l", list_path=list_path, parallel=0)
    # 3. List parallel.
    ref_lp = _stats_for_expr(expr, "l", list_path=list_path, parallel=4)

    # Mode consistency: per-LC value identical across -i, -l serial,
    # -l parallel.  rtol 1e-12 because only float-summation order
    # could differ, and these all walk the same code path serially.
    for p in lc_paths:
        assert p in ref_i, f"{label}: -i mode missing entry for {p}"
        assert p in ref_l, f"{label}: -l mode missing entry for {p}"
        assert p in ref_lp, f"{label}: -l-parallel mode missing entry for {p}"
        assert np.isclose(ref_i[p], ref_l[p], rtol=1e-12, atol=1e-12), (
            f"{label}: -i vs -l mismatch on {p}: "
            f"-i={ref_i[p]:.15g} -l={ref_l[p]:.15g}"
        )
        assert np.isclose(ref_i[p], ref_lp[p], rtol=1e-12, atol=1e-12), (
            f"{label}: -i vs -l-parallel mismatch on {p}: "
            f"-i={ref_i[p]:.15g} -lp={ref_lp[p]:.15g}"
        )

    # Numpy ground-truth check (skip for convention-dependent cases).
    if kind == "consistency_only":
        return
    for p, (_, m, e) in zip(lc_paths, lc_data):
        if kind == "perpoint":
            exp = float(np.median(fn(m, e)))
        else:  # "scalar"
            exp = float(fn(m, e))
        # Loose tolerance: vartools' summation order can differ from
        # numpy's at ~1e-7 relative for per-point expressions.
        assert np.isclose(ref_i[p], exp, rtol=1e-7, atol=1e-9), (
            f"{label}: vartools vs numpy mismatch on {p}: "
            f"vartools={ref_i[p]:.15g} numpy={exp:.15g}"
        )
