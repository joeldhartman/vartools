"""First-time-user scenarios drawn from reading the webdocs.

Each test simulates something a user might naturally try after reading the
top-level API overview, the pipeline.md "what library mode supports"
matrix, the cmd.o(outname=PerLC([...])) section, and the LightCurveBatch
chaining docs.  Errors are evaluated for clarity, not just type.

Run with:
    conda run -n pyvartools python -m pytest python/tests/test_doc_user_scenarios.py -v
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

VARTOOLS_SRC = os.path.realpath(
    os.path.join(os.path.dirname(__file__), "..", "..")
)
PYTHON_DIR = os.path.join(VARTOOLS_SRC, "python")
if PYTHON_DIR not in sys.path:
    sys.path.insert(0, PYTHON_DIR)

import pyvartools as vt
from pyvartools import commands as cmd

# Point the binary at a real path if we can find one.
if not os.environ.get("VARTOOLS_BINARY"):
    for cand in (
        os.path.join(VARTOOLS_SRC, "..", "..", "bin", "vartools"),
        "/home/jhartman/HATpipebin/bin/vartools",
    ):
        cand = os.path.realpath(cand)
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            vt.set_binary(cand)
            break

try:
    _HAVE_BINARY = os.path.isfile(vt.get_binary())
except FileNotFoundError:
    _HAVE_BINARY = False

needs_binary = pytest.mark.skipif(not _HAVE_BINARY,
                                  reason="vartools binary required")


def _three_lcs(n=3, npts=150):
    """Three distinct sinusoidal LCs."""
    out = []
    for i in range(n):
        t = np.linspace(0, 20, npts)
        mag = 10.0 + 0.05 * np.sin(2 * np.pi * t / (1.0 + 0.3 * i))
        err = np.full(npts, 0.005)
        out.append(vt.LightCurve.from_arrays(t, mag, err, name=f"star{i}"))
    return out


# ---------------------------------------------------------------------------
# Happy paths drawn directly from the docs
# ---------------------------------------------------------------------------

class TestHappyPath:
    """Each test mirrors a documented example with the user's expected
    outcome.  A failure here means a doc example is broken."""

    @needs_binary
    def test_first_use_rms_top_level(self):
        """index.md: vt.rms(lc) returns a Result with RMS_0."""
        lc = _three_lcs(n=1)[0]
        result = vt.rms(lc)
        assert "RMS_0" in result.vars.index
        assert float(result.vars["RMS_0"]) > 0

    @needs_binary
    def test_pipeline_run_batch_basic(self):
        """pipeline.md: Pipeline().clip().rms().run_batch(lcs) returns
        one row per LC."""
        lcs = _three_lcs(n=3)
        batch = vt.Pipeline([cmd.clip(5.0), cmd.rms()]).run_batch(lcs)
        assert len(batch.vars) == 3
        assert "RMS_1" in batch.vars.columns  # rms is the 2nd command (_1)

    @needs_binary
    def test_capture_lc_returns_lightcurve(self):
        """pipeline.md: capture_lc=True populates result.lc."""
        lc = _three_lcs(n=1)[0]
        result = vt.Pipeline([cmd.clip(5.0)]).run(lc, capture_lc=True)
        assert result.lc is not None
        assert isinstance(result.lc, vt.LightCurve)
        assert len(result.lc.t) <= len(lc.t)  # clip may remove points

    @needs_binary
    def test_cmd_o_outname_perlc_documented_example(self, tmp_path):
        """pipeline.md / control-flow.md / chaining.md all show this
        pattern.  Files land at <outdir>/<perlc_value[i]>."""
        lcs = _three_lcs(n=3)
        outdir = tmp_path / "out"
        outdir.mkdir()
        names = [f"star{i}.lc" for i in range(3)]

        batch = vt.Pipeline([
            cmd.clip(5.0),
            cmd.o(outdir=str(outdir), outname=vt.PerLC(names),
                  allcols=True),
        ]).run_batch(lcs)

        assert sorted(os.listdir(outdir)) == sorted(names)
        assert len(batch.vars) == 3

    @needs_binary
    def test_chaining_on_lightcurve(self):
        """chaining.md: lc.LS(...) returns a Result."""
        lc = _three_lcs(n=1)[0]
        r = lc.LS(0.5, 10.0, 1e-3)
        assert "LS_Period_1_0" in r.vars.index

    @needs_binary
    def test_lightcurvebatch_chaining_with_perlc(self):
        """chaining.md: per-LC array parameters work in the chaining API."""
        lcs = _three_lcs(n=3)
        result = (vt.LightCurveBatch(lcs)
                  .LS(minp=np.array([0.3, 0.5, 0.7]),
                      maxp=10.0, subsample=0.1)
                  .run())
        assert len(result.vars) == 3
        # Three distinct LS periods (different per-LC minp shaped the search).
        periods = result.vars["LS_Period_1_0"].tolist()
        assert all(0.3 <= p <= 10.0 for p in periods)

    @needs_binary
    def test_batchresult_continuation_carries_scalars(self):
        """chaining.md: br1.expr(...).run() sees prior LS values."""
        lcs = _three_lcs(n=3)
        br1 = vt.LightCurveBatch(lcs).LS(0.5, 10.0, 1e-3).run()
        br2 = br1.expr("doubled=2*LS_Period_1_0",
                       vartype="scalar").run()
        # Each LC's continuation sees ITS OWN prior period.
        for i, lc_result in enumerate(br2.lcs):
            assert lc_result is not None
            orig = br1._result_at(i).vars["LS_Period_1_0"]
            assert abs(lc_result.scalars["doubled"] - 2 * orig) < 1e-9

    @needs_binary
    def test_cmd_o_capture_intermediate(self):
        """control-flow.md: cmd.o(capture=True) snapshots the LC mid-pipeline."""
        lc = _three_lcs(n=1)[0]
        result = vt.Pipeline([
            cmd.clip(5.0),
            cmd.o(capture=True, key="after_clip"),
            cmd.LS(0.5, 10.0, 1e-3, npeaks=1),
        ]).run(lc)
        assert "after_clip" in result.files
        clipped = result.files["after_clip"]
        assert isinstance(clipped, vt.LightCurve)
        assert "t" in clipped._df.columns


# ---------------------------------------------------------------------------
# Error message quality — things a user could plausibly get wrong
# ---------------------------------------------------------------------------

class TestErrorMessageQuality:
    """User does something the docs implicitly forbid.  Does the error tell
    them what they did wrong AND how to fix it?"""

    @needs_binary
    def test_outname_perlc_length_mismatch(self):
        """User gives 4 output names but only 3 LCs."""
        lcs = _three_lcs(n=3)
        pipe = vt.Pipeline([cmd.o(outdir="/tmp",
                                   outname=vt.PerLC(["a", "b", "c", "d"]))])
        with pytest.raises(ValueError) as excinfo:
            pipe.run_batch(lcs)
        msg = str(excinfo.value)
        # Useful error should mention: actual count, expected count, the
        # variable, and ideally the cmd.o context.
        assert "4 values" in msg or "4" in msg, f"missing length detail: {msg}"
        assert "3" in msg, f"missing batch-size detail: {msg}"

    @needs_binary
    def test_outname_perlc_without_outdir(self):
        """User forgets outdir=."""
        lcs = _three_lcs(n=2)
        pipe = vt.Pipeline([cmd.o(outname=vt.PerLC(["a", "b"]))])
        with pytest.raises(ValueError) as excinfo:
            pipe.run_batch(lcs)
        msg = str(excinfo.value)
        # User should see "outdir" in the suggestion.
        assert "outdir" in msg, f"error doesn't mention outdir: {msg}"

    @needs_binary
    def test_plain_list_not_autodetected_as_perlc(self):
        """chaining.md: 'plain Python lists are NOT auto-detected.'
        User passes minp=[0.3, 0.5, 0.7] (plain list, not numpy/PerLC).
        Should NOT be treated as per-LC; if vartools accepts it as a
        scalar, the user gets the wrong answer silently.  Verify the
        documented behaviour."""
        lcs = _three_lcs(n=3)
        # A plain list is treated as a fixed multi-valued parameter (e.g.
        # ld_coeffs).  For LS.minp which is scalar, it should error or be
        # passed verbatim and trip vartools' parser.
        with pytest.raises(Exception):
            vt.Pipeline([cmd.LS(minp=[0.3, 0.5, 0.7], maxp=10.0,
                                 subsample=0.1)]).run_batch(lcs)

    @needs_binary
    def test_dataframe_with_wrong_column_names_silently_fails(self):
        """**Known UX gap (pre-existing):** a DataFrame whose columns are
        'time' / 'magnitude' / 'error' (rather than 't'/'mag'/'err') is
        accepted without complaint, but produces garbage results
        (Npoints=0, RMS=-1).  No exception, no warning.

        The docstring for ``LightCurve.from_dataframe`` says "any
        DataFrame is accepted; columns named t/mag/err are treated as the
        standard vartools vectors when present" — but the failure mode
        when those columns are *missing* is a silent zero-row LC, which
        a casual user is unlikely to notice.

        This test pins the current behaviour so the gap is visible in
        the test suite; tighten it (validate column names, or
        auto-rename common aliases) as a separate UX fix."""
        df = pd.DataFrame({
            "time": np.linspace(0, 10, 50),
            "magnitude": np.full(50, 10.0),
            "error": np.full(50, 0.01),
        })
        result = vt.rms(df)  # silently succeeds
        # Sentinel values show the LC was empty for vartools.
        assert int(result.vars["Npoints_0"]) == 0
        assert float(result.vars["Mean_Mag_0"]) == -1.0
        assert float(result.vars["RMS_0"]) == -1.0

    @needs_binary
    def test_empty_batch(self):
        """User passes an empty list to run_batch."""
        # No assertion — just check what happens.  Worth knowing.
        try:
            batch = vt.Pipeline([cmd.rms()]).run_batch([])
            print(f"empty batch returned: vars={batch.vars.shape}, "
                  f"lcs={batch.lcs}")
        except Exception as e:
            print(f"empty batch raised: {type(e).__name__}: {e}")


# ---------------------------------------------------------------------------
# Surprising-but-reasonable: a user is exploring, not yet sure of the rules
# ---------------------------------------------------------------------------

class TestSurprisingButReasonable:

    @needs_binary
    def test_outname_perlc_in_single_lc_run(self, tmp_path):
        """User writes a Pipeline with outname=PerLC, then calls .run(lc)
        on a single LC instead of run_batch.  What happens?  The docs say
        outname=PerLC is for batch; we'd hope for a clear error pointing
        the user at run_batch."""
        lc = _three_lcs(n=1)[0]
        pipe = vt.Pipeline([cmd.o(outdir=str(tmp_path),
                                    outname=vt.PerLC(["only.lc"]))])
        # Single-LC run shouldn't apply the auto-rewrite (it's run_batch
        # only).  Expect an error.
        try:
            r = pipe.run(lc)
            print(f"single-LC + outname=PerLC unexpectedly succeeded: "
                  f"{r.vars}")
        except Exception as e:
            print(f"single-LC + outname=PerLC raised: {type(e).__name__}: {e}")

    @needs_binary
    def test_perlc_all_same_value(self, tmp_path):
        """Edge case: PerLC with all identical values.  Functionally
        equivalent to a scalar; should still work."""
        lcs = _three_lcs(n=3)
        outdir = tmp_path / "same"
        outdir.mkdir()
        names = ["same.lc", "same.lc", "same.lc"]
        # Three LCs writing to the same name — the third write overwrites
        # the first two.  We'd hope for either an explicit error OR a
        # documented "last wins" semantic.
        try:
            pipe = vt.Pipeline([
                cmd.o(outdir=str(outdir),
                      outname=vt.PerLC(names), allcols=True),
            ])
            batch = pipe.run_batch(lcs)
            written = sorted(os.listdir(outdir))
            print(f"PerLC with duplicate names — written: {written}, "
                  f"batch.vars rows: {len(batch.vars)}")
            # Worth knowing whether the user gets a single file with the
            # last LC's data, vs 3 distinct files, vs an error.
        except Exception as e:
            print(f"PerLC duplicate names raised: "
                  f"{type(e).__name__}: {e}")

    @needs_binary
    def test_perlc_on_npeaks_parameter(self):
        """User tries to vary LS.npeaks (an integer) per LC.  The CLI
        doesn't accept var/expr for npeaks — what does the user see?"""
        lcs = _three_lcs(n=3)
        try:
            batch = vt.Pipeline([
                cmd.LS(0.5, 10.0, 0.1, npeaks=vt.PerLC([1, 2, 3])),
            ]).run_batch(lcs)
            print(f"PerLC npeaks unexpectedly succeeded: "
                  f"columns={list(batch.vars.columns)[:5]}...")
        except Exception as e:
            print(f"PerLC npeaks raised: {type(e).__name__}: {e}")

    @needs_binary
    def test_chain_continuation_after_outname_perlc(self, tmp_path):
        """User has a batch with cmd.o(outname=PerLC), then continues
        the chain.  Does the PerLC affect carry-forward?"""
        lcs = _three_lcs(n=3)
        outdir = tmp_path / "step1"
        outdir.mkdir()
        names = [f"s{i}.lc" for i in range(3)]
        br1 = vt.Pipeline([
            cmd.clip(5.0),
            cmd.LS(0.5, 10.0, 1e-3),
            cmd.o(outdir=str(outdir), outname=vt.PerLC(names),
                  capture=True, key="captured"),
        ]).run_batch(lcs, capture_lc=True)

        # Continuation: should still see LS_Period_1_1 from prior segment.
        # (cmd.LS is at command index 1; cmd.o is at index 2.)
        # User expectation: this just works like any other chain.
        try:
            assert br1.lcs is not None
            from pyvartools._batch import LightCurveBatch
            br2 = LightCurveBatch(br1.lcs, _prior_batch=br1).rms().run()
            # rms is appended after the prior 3 commands → suffix _3
            assert "RMS_3" in br2.vars.columns, (
                f"expected RMS_3 in {list(br2.vars.columns)}")
            print(f"chain after outname=PerLC: PASS  "
                  f"({len(br2.vars)} rows)")
        except Exception as e:
            print(f"chain after outname=PerLC failed: "
                  f"{type(e).__name__}: {e}")
            raise
