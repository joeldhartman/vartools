"""Automation-script user — a different new-user persona.

This user is wiring pyvartools into a production-ish data pipeline.
Concerns:
  - Persist a Pipeline across processes (pickle/dill, save+load).
  - Persist a Result so they can refer back to it later.
  - Robust handling of weird inputs: NaN, inf, single-point LCs.
  - Per-LC error recovery — one bad LC shouldn't take down the batch.
  - Reproducibility: same Pipeline + same input = same output.
  - Concatenating / filtering BatchResult dataframes downstream.

Compared to v1 (generic exploration) and v2 (survey-scientist
ergonomics), this exercise probes:
  - Serialization (a fundamental Python expectation).
  - Result composition with the wider scientific Python stack.
  - Degenerate inputs that would crash naïve pipelines.

Run with:
    conda run -n pyvartools python -m pytest python/tests/test_doc_user_scenarios_v3.py -v -s
"""

from __future__ import annotations

import math
import os
import pickle
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


def _clean_lcs(n=3, npts=200):
    """Three sinusoidal LCs, distinct periods, mild noise."""
    rng = np.random.default_rng(0)
    out = []
    for i in range(n):
        period = 1.0 + 0.4 * i
        t = np.linspace(0, 30, npts)
        mag = 10.0 + 0.05 * np.sin(2 * np.pi * t / period)
        mag += rng.normal(0, 0.003, npts)
        err = np.full(npts, 0.005)
        out.append(vt.LightCurve.from_arrays(
            t, mag, err, name=f"obj{i}"))
    return out


# ---------------------------------------------------------------------------
# Persistence — pickle is a basic Python expectation
# ---------------------------------------------------------------------------

class TestPersistence:

    def test_pipeline_pickles_and_unpickles(self):
        """Round-trip a Pipeline through pickle.  After unpickle, the
        commands list should be intact and identical."""
        pipe = vt.Pipeline([
            cmd.clip(5.0),
            cmd.LS(0.5, 10.0, 0.01, npeaks=2),
            cmd.rms(),
        ])
        blob = pickle.dumps(pipe)
        loaded = pickle.loads(blob)
        assert len(loaded.commands) == len(pipe.commands)
        # Commands are equal if their reprs are equal — relies on the
        # Python-shaped repr from v2 cleanup.
        for a, b in zip(loaded.commands, pipe.commands):
            assert repr(a) == repr(b)

    @needs_binary
    def test_unpickled_pipeline_still_runs(self):
        """A Pipeline that's been pickled and reloaded produces the same
        results as the original — a real-world reproducibility check."""
        lcs = _clean_lcs(n=2)
        pipe = vt.Pipeline([cmd.LS(0.5, 10.0, 0.01, npeaks=1)])
        orig = pipe.run_batch(lcs, randseed=7)
        reloaded = pickle.loads(pickle.dumps(pipe))
        again = reloaded.run_batch(lcs, randseed=7)
        np.testing.assert_array_equal(
            orig.vars["LS_Period_1_0"].values,
            again.vars["LS_Period_1_0"].values,
        )

    @needs_binary
    def test_batchresult_pickles(self):
        """A BatchResult should survive pickle so users can stash it."""
        lcs = _clean_lcs(n=3)
        batch = vt.Pipeline([cmd.LS(0.5, 10.0, 0.01, npeaks=1)]).run_batch(lcs)
        try:
            blob = pickle.dumps(batch)
            loaded = pickle.loads(blob)
            pd.testing.assert_frame_equal(loaded.vars, batch.vars)
        except Exception as e:
            pytest.fail(f"BatchResult does not pickle cleanly: "
                        f"{type(e).__name__}: {e}")

    def test_lightcurve_pickles(self):
        """Plain LightCurve round-trips through pickle."""
        lc = _clean_lcs(n=1)[0]
        loaded = pickle.loads(pickle.dumps(lc))
        assert loaded.name == lc.name
        np.testing.assert_array_equal(loaded.t, lc.t)
        np.testing.assert_array_equal(loaded.mag, lc.mag)
        np.testing.assert_array_equal(loaded.err, lc.err)


# ---------------------------------------------------------------------------
# Edge-case inputs — what if the data is weird?
# ---------------------------------------------------------------------------

class TestWeirdInputs:

    @needs_binary
    def test_single_point_lc(self):
        """An LC with exactly 1 point.  Statistically nonsensical for most
        commands but should fail informatively, not crash."""
        lc = vt.LightCurve.from_arrays(
            t=np.array([0.0]), mag=np.array([10.0]), err=np.array([0.01]))
        try:
            r = vt.rms(lc)
            print(f"single-point rms: vars={dict(r.vars)}")
        except Exception as e:
            print(f"single-point rms raised: {type(e).__name__}: {e}")

    @needs_binary
    def test_lc_with_nan_in_mag(self):
        """Some LCs have NaN gaps.  Does vartools handle them?"""
        t = np.linspace(0, 10, 100)
        mag = 10.0 + 0.05 * np.sin(2 * np.pi * t / 1.7)
        mag[10:15] = np.nan
        err = np.full(100, 0.005)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="nan_lc")
        try:
            r = vt.rms(lc)
            print(f"NaN-in-mag rms: vars={dict(r.vars)}")
            # If vartools accepts NaN, RMS will likely be NaN too.
            # That's fine to document.
        except Exception as e:
            print(f"NaN-in-mag raised: {type(e).__name__}: {e}")

    @needs_binary
    def test_lc_with_inf_in_mag(self):
        """Pathological: ±inf in mag."""
        t = np.linspace(0, 10, 50)
        mag = np.full(50, 10.0)
        mag[25] = np.inf
        err = np.full(50, 0.005)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="inf_lc")
        try:
            r = vt.rms(lc)
            print(f"inf-in-mag rms: vars={dict(r.vars)}")
        except Exception as e:
            print(f"inf-in-mag raised: {type(e).__name__}: {e}")

    @needs_binary
    def test_lc_with_duplicate_timestamps(self):
        """Real-world LCs sometimes have duplicate timestamps from blends
        or stacked exposures.  Vartools' LS would normally still run; we
        want to confirm pyvartools doesn't reject the LC."""
        t = np.concatenate([np.linspace(0, 5, 50),
                            np.linspace(0, 5, 50)])  # duplicates
        mag = np.full(100, 10.0)
        err = np.full(100, 0.005)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="dup_t")
        r = vt.rms(lc)
        assert int(r.vars["Npoints_0"]) == 100

    @needs_binary
    def test_lc_with_unsorted_time(self):
        """Vartools generally expects time-sorted LCs.  What happens if
        the user supplies an unsorted one?"""
        t = np.array([1.0, 5.0, 2.0, 8.0, 3.0])
        mag = np.full(5, 10.0)
        err = np.full(5, 0.01)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="unsorted")
        try:
            r = vt.rms(lc)
            print(f"unsorted-time rms: vars={dict(r.vars)}")
        except Exception as e:
            print(f"unsorted-time raised: {type(e).__name__}: {e}")

    @needs_binary
    def test_empty_pipeline(self):
        """Pipeline with no commands.  Should it error, or pass-through?"""
        lc = _clean_lcs(n=1)[0]
        pipe = vt.Pipeline([])
        try:
            r = pipe.run(lc)
            print(f"empty pipeline: vars={dict(r.vars)}")
        except Exception as e:
            print(f"empty pipeline raised: {type(e).__name__}: {e}")


# ---------------------------------------------------------------------------
# Result composition
# ---------------------------------------------------------------------------

class TestResultComposition:

    @needs_binary
    def test_concat_two_batches(self):
        """User runs the pipeline on two disjoint LC sets and wants to
        combine the resulting BatchResult.vars DataFrames."""
        lcs_a = _clean_lcs(n=2)
        lcs_b = _clean_lcs(n=3)
        for i, lc in enumerate(lcs_b):
            lc.name = f"obj_b{i}"
        pipe = vt.Pipeline([cmd.LS(0.5, 10.0, 0.01, npeaks=1)])
        a = pipe.run_batch(lcs_a)
        b = pipe.run_batch(lcs_b)
        merged = pd.concat([a.vars, b.vars], ignore_index=True)
        assert len(merged) == 5
        assert set(merged["Name"]) == \
               set(["obj0", "obj1", "obj_b0", "obj_b1", "obj_b2"])

    @needs_binary
    def test_filter_batchresult_vars(self):
        """Pandas-style filtering on batch.vars."""
        lcs = _clean_lcs(n=5)
        batch = vt.Pipeline([cmd.LS(0.5, 10.0, 0.01, npeaks=1)]).run_batch(lcs)
        # Filter to LCs whose found period is short.
        short = batch.vars[batch.vars["LS_Period_1_0"] < 1.5]
        assert isinstance(short, pd.DataFrame)
        assert all(p < 1.5 for p in short["LS_Period_1_0"])

    @needs_binary
    def test_access_nonexistent_command_var_raises(self):
        """A user assumes there's an LS command but didn't add one.
        Accessing result.varobjs.LS should error helpfully, not return
        a None or empty object."""
        lc = _clean_lcs(n=1)[0]
        r = vt.Pipeline([cmd.rms()]).run(lc)
        try:
            v = r.varobjs.LS.Period_1
            print(f"varobjs.LS access without LS in pipeline: {v}")
        except AttributeError as e:
            print(f"varobjs.LS missing-cmd error: {e}")
        except Exception as e:
            print(f"varobjs.LS unexpected: {type(e).__name__}: {e}")


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

class TestReproducibility:

    @needs_binary
    def test_same_inputs_same_outputs(self):
        """Running the same Pipeline twice in the same process produces
        identical results."""
        lcs = _clean_lcs(n=3)
        pipe = vt.Pipeline([cmd.LS(0.5, 10.0, 0.01, npeaks=1), cmd.rms()])
        r1 = pipe.run_batch(lcs)
        r2 = pipe.run_batch(lcs)
        pd.testing.assert_frame_equal(r1.vars, r2.vars)

    @needs_binary
    def test_run_then_run_batch_consistent(self):
        """rms() on one LC via run() == that LC's row in run_batch([lc])."""
        lc = _clean_lcs(n=1)[0]
        pipe = vt.Pipeline([cmd.rms()])
        single = pipe.run(lc)
        batched = pipe.run_batch([lc])
        assert math.isclose(float(single.vars["RMS_0"]),
                            float(batched.vars["RMS_0"].iloc[0]),
                            rel_tol=1e-12)


# ---------------------------------------------------------------------------
# Per-LC error recovery in batch
# ---------------------------------------------------------------------------

class TestErrorRecovery:

    @needs_binary
    def test_one_bad_lc_in_batch_run_via_lightcurvebatch(self):
        """LightCurveBatch should isolate per-LC failures so other LCs
        complete."""
        good = _clean_lcs(n=2)
        empty = vt.LightCurve.from_arrays(
            np.array([]), np.array([]), np.array([]), name="empty")
        batch = vt.LightCurveBatch(good + [empty]).rms().run()
        oks = [r for r in batch if r.ok]
        bads = [r for r in batch if not r.ok]
        print(f"per-LC isolation: ok={len(oks)}, bad={len(bads)}, "
              f"bad_error={bads[0].error if bads else 'n/a'}")
        # Expectation: the empty LC is the failing one; the two good
        # LCs produced results.
        assert len(oks) >= 2, "good LCs should complete despite the bad one"

    @needs_binary
    def test_run_batch_bad_lc_raises_by_default(self):
        """run_batch with a bad LC raises when raise_on_error=True
        (the default).  The error message names the offending LC."""
        good = _clean_lcs(n=2)
        empty = vt.LightCurve.from_arrays(
            np.array([]), np.array([]), np.array([]), name="empty_lc")
        pipe = vt.Pipeline([cmd.rms()])
        with pytest.raises(Exception) as excinfo:
            pipe.run_batch(good + [empty])
        # Error message should help the user identify which LC failed.
        assert "empty_lc" in str(excinfo.value) or \
               "zero rows" in str(excinfo.value)

    @needs_binary
    def test_run_batch_bad_lc_with_raise_on_error_false(self):
        """raise_on_error=False suppresses the per-LC exception and
        surfaces it via BatchResult.error instead — even when the
        failure happens inside the library-mode per-LC loop."""
        good = _clean_lcs(n=2)
        empty = vt.LightCurve.from_arrays(
            np.array([]), np.array([]), np.array([]), name="empty_lc")
        pipe = vt.Pipeline([cmd.rms()])
        batch = pipe.run_batch(good + [empty], raise_on_error=False)
        assert batch.error is not None
        assert ("empty_lc" in str(batch.error)
                or "zero rows" in str(batch.error))
        # vars is empty because library batch is one-shot — partial
        # results aren't reported on failure.  Documenting that here.
        assert batch.vars.empty


# ---------------------------------------------------------------------------
# Numerical edge cases
# ---------------------------------------------------------------------------

class TestNumerical:

    @needs_binary
    def test_dtype_int_vs_float(self):
        """User builds an LC with int-typed arrays.  pyvartools should
        coerce or handle this without surprise."""
        t = np.arange(50, dtype=np.int32)
        mag = np.full(50, 10, dtype=np.int64)  # int mag!
        err = np.full(50, 1, dtype=np.int32)
        try:
            lc = vt.LightCurve.from_arrays(t, mag, err, name="intish")
            r = vt.rms(lc)
            print(f"int-typed arrays rms: vars={dict(r.vars)}")
            assert int(r.vars["Npoints_0"]) == 50
        except Exception as e:
            print(f"int-typed arrays raised: {type(e).__name__}: {e}")

    @needs_binary
    def test_very_short_lc(self):
        """LC of 3 points.  Vartools' LS / rms / etc. need a minimum number
        of points; what's the failure mode?"""
        t = np.array([0.0, 1.0, 2.0])
        mag = np.array([10.0, 10.1, 10.05])
        err = np.array([0.01, 0.01, 0.01])
        lc = vt.LightCurve.from_arrays(t, mag, err, name="tiny")
        try:
            r = vt.Pipeline([cmd.LS(0.5, 1.5, 0.01, npeaks=1)]).run(lc)
            print(f"3-point LS: vars={dict(r.vars)}")
        except Exception as e:
            print(f"3-point LS raised: {type(e).__name__}: {e}")
