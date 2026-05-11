"""
Feature C tests: library embedding (init-once / process-many).

These tests are written against the LibPipeline class and Pipeline's
library-mode dispatch.  They FAIL until:
  - libvartoolspipeline.so is built (Step 2 in the plan)
  - pyvartools/_libpipeline.py is written (Step 3)
  - Pipeline.run() dispatches to LibPipeline when the library is available

Run in library mode:    VARTOOLS_USE_LIBRARY=1 pytest tests/test_library_mode.py
Run in subprocess mode: VARTOOLS_USE_LIBRARY=0 pytest tests/test_library_mode.py
"""

from __future__ import annotations

import math
import os
import sys
from unittest.mock import patch, MagicMock

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------

VARTOOLS_SRC = os.path.realpath(
    os.path.join(os.path.dirname(__file__), "..", "..")
)
PYTHON_DIR = os.path.join(VARTOOLS_SRC, "python")
if PYTHON_DIR not in sys.path:
    sys.path.insert(0, PYTHON_DIR)

import pyvartools as vt
from pyvartools import commands as cmd
from pyvartools.results import RunError

# Point to the installed binary if not overridden via environment.
if not os.environ.get("VARTOOLS_BINARY"):
    _candidates = [
        os.path.join(VARTOOLS_SRC, "..", "..", "bin", "vartools"),
        os.path.join(VARTOOLS_SRC, "src", "vartools"),
    ]
    for _cand in _candidates:
        _cand = os.path.realpath(_cand)
        if os.path.isfile(_cand) and os.access(_cand, os.X_OK):
            vt.set_binary(_cand)
            break

try:
    _HAVE_BINARY = os.path.isfile(vt.get_binary())
except FileNotFoundError:
    _HAVE_BINARY = False

needs_binary = pytest.mark.skipif(
    not _HAVE_BINARY, reason="vartools binary not available"
)

# Check whether the C library is available
try:
    from pyvartools._libpipeline import LibPipeline
    _HAVE_LIBRARY = True
except ImportError:
    _HAVE_LIBRARY = False
    LibPipeline = None

needs_library = pytest.mark.skipif(
    not _HAVE_LIBRARY, reason="libvartoolspipeline not available (Step 2+3 not done)"
)

EXAMPLES_DIR = os.path.join(VARTOOLS_SRC, "EXAMPLES")
EXAMPLE_LC = os.path.join(EXAMPLES_DIR, "2")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_lc():
    t = np.linspace(0, 10, 100)
    mag = np.full(100, 10.0)
    err = np.full(100, 0.01)
    return vt.LightCurve.from_arrays(t, mag, err, name="testlc")


@pytest.fixture
def sinusoidal_lc():
    rng = np.random.default_rng(42)
    t = np.linspace(0, 30, 200)
    mag = 10.0 + 0.1 * np.sin(2 * np.pi * t / 1.5)
    err = np.full(200, 0.005)
    return vt.LightCurve.from_arrays(t, mag, err, name="sinlc")


def make_lcs(n=3, npts=150):
    """Return n distinct sinusoidal LightCurves with different RMS values."""
    lcs = []
    for i in range(n):
        t = np.linspace(0, 20, npts)
        mag = 10.0 + (i + 1) * 0.05 * np.sin(2 * np.pi * t / (1.0 + i * 0.3))
        err = np.full(npts, 0.005)
        lcs.append(vt.LightCurve.from_arrays(t, mag, err, name=f"lc{i}"))
    return lcs


# ---------------------------------------------------------------------------
# LibPipeline unit tests
# ---------------------------------------------------------------------------

class TestLibPipelineImport:

    def test_libpipeline_import(self):
        """from pyvartools._libpipeline import LibPipeline must not raise."""
        try:
            from pyvartools._libpipeline import LibPipeline as LP
        except ImportError as e:
            pytest.fail(
                f"ImportError: {e}\n"
                "libvartoolspipeline.so has not been built yet (Step 2)."
            )


@needs_library
class TestLibPipelineUnit:

    def test_libpipeline_init_rms(self):
        """LibPipeline initializes without error for a simple -rms pipeline."""
        lp = LibPipeline(["-rms", "-oneline"])
        assert lp._p is not None

    def test_libpipeline_process_lc_rms(self):
        """process_lc returns a Series with RMS_0 as a finite float."""
        lp = LibPipeline(["-rms", "-oneline"])
        t = np.linspace(0, 10, 100)
        mag = np.full(100, 10.0)
        err = np.full(100, 0.01)
        var = lp.process_lc(t, mag, err, name="test")
        assert "RMS_0" in var.index
        assert math.isfinite(float(var["RMS_0"]))

    def test_libpipeline_process_multiple_lcs_same_pipeline(self):
        """Init-once / process-many: three LCs give independent (different) RMS values."""
        lp = LibPipeline(["-rms", "-oneline"])
        results = []
        for i in range(1, 4):
            t = np.linspace(0, 10, 100)
            mag = 10.0 + i * 0.01 * np.sin(2 * np.pi * t)
            err = np.full(100, 0.005)
            var = lp.process_lc(t, mag, err, name=f"lc{i}")
            results.append(float(var["RMS_0"]))

        # All values should be finite
        for v in results:
            assert math.isfinite(v), f"RMS_0 is not finite: {v}"

        # Values should differ across LCs (different signals → different RMS)
        assert len(set(results)) > 1, (
            "All three RMS values are identical — "
            "parsecommandline may be re-running each call (init-once broken)"
        )

    def test_libpipeline_process_lc_ls(self, sinusoidal_lc):
        """LS period recovery: recovered period within 1% of injected 1.5 d."""
        lp = LibPipeline(["-LS", "0.5", "5.0", "0.001", "1", "0", "-oneline"])
        var = lp.process_lc(
            sinusoidal_lc.t, sinusoidal_lc.mag, sinusoidal_lc.err, name="sinlc"
        )
        assert "LS_Period_1_0" in var.index
        period = float(var["LS_Period_1_0"])
        assert math.isfinite(period)
        assert abs(period - 1.5) / 1.5 < 0.01, (
            f"Recovered period {period:.4f} d is not within 1% of 1.5 d"
        )

    def test_libpipeline_process_lc_clip_then_ls(self):
        """clip + LS pipeline: LS recovers correct period even with outliers."""
        rng = np.random.default_rng(7)
        t = np.linspace(0, 30, 300)
        mag = 10.0 + 0.1 * np.sin(2 * np.pi * t / 2.1)
        err = np.full(300, 0.005)
        # Inject outliers
        outlier_idx = rng.choice(300, 10, replace=False)
        mag[outlier_idx] += rng.uniform(0.5, 1.5, 10)

        lp = LibPipeline(["-clip", "5.0", "1", "-LS", "0.5", "5.0", "0.001", "1", "0",
                          "-oneline"])
        var = lp.process_lc(t, mag, err, name="cliplc")
        assert "LS_Period_1_1" in var.index
        period = float(var["LS_Period_1_1"])
        assert math.isfinite(period)
        assert abs(period - 2.1) / 2.1 < 0.01

    def test_libpipeline_large_lc(self):
        """10,000-point LC processes without error or segfault."""
        lp = LibPipeline(["-rms", "-oneline"])
        t = np.linspace(0, 1000, 10000)
        mag = 10.0 + 0.01 * np.sin(2 * np.pi * t / 2.0)
        err = np.full(10000, 0.005)
        var = lp.process_lc(t, mag, err, name="biglc")
        assert "RMS_0" in var.index
        assert math.isfinite(float(var["RMS_0"]))

    @needs_binary
    def test_libpipeline_stats_match_subprocess(self, sinusoidal_lc):
        """LibPipeline var agree with subprocess Pipeline to within 1e-6."""
        # Subprocess reference
        pipe = vt.Pipeline([cmd.rms()])
        result_sub = pipe.run(sinusoidal_lc)

        # Library
        lp = LibPipeline(["-rms", "-oneline"])
        stats_lib = lp.process_lc(
            sinusoidal_lc.t, sinusoidal_lc.mag, sinusoidal_lc.err,
            name=sinusoidal_lc.name
        )

        for key in ["RMS_0"]:
            v_lib = float(stats_lib[key])
            v_sub = float(result_sub.vars[key])
            assert abs(v_lib - v_sub) < 1e-6, (
                f"{key}: library={v_lib} vs subprocess={v_sub}"
            )

    def test_libpipeline_free(self):
        """Deleting a LibPipeline calls vartools_free_pipeline without segfault."""
        lp = LibPipeline(["-rms", "-oneline"])
        assert lp._p is not None
        del lp  # should not segfault

    def test_set_inlist_value_double_observable(self):
        """set_inlist_value('v', 3.14) is observable via -stats on the var."""
        t = np.linspace(0, 1, 10)
        mag = np.full(10, 10.0)
        err = np.full(10, 0.01)
        lp = LibPipeline(['-inlistvars', 'v:0:double:0.0',
                          '-expr', 'vv=v', '-stats', 'vv', 'mean', '-oneline'])
        for val in (3.14, 100.0, -1.5):
            lp.set_inlist_value('v', val)
            stats = lp.process_lc(t, mag, err, name='test')
            assert math.isclose(float(stats['STATS_vv_MEAN_1']), val,
                                rel_tol=1e-9, abs_tol=1e-9), (
                f"set v={val} but observed {stats['STATS_vv_MEAN_1']}")

    def test_set_inlist_value_int_dispatches_to_long(self):
        """Python int dispatches to vartools_set_inlist_long; bool→int."""
        t = np.linspace(0, 1, 10)
        mag = np.full(10, 10.0)
        err = np.full(10, 0.01)
        lp = LibPipeline(['-inlistvars', 'flag:0:int:0',
                          '-expr', 'ff=flag', '-stats', 'ff', 'mean', '-oneline'])
        # int → long path; the C int setter is selected by datatype, so this
        # only works because flag's declared type is int.  Bool gets routed
        # through the C int setter directly.
        lp.set_inlist_value('flag', True)
        assert int(lp.process_lc(t, mag, err)['STATS_ff_MEAN_1']) == 1
        lp.set_inlist_value('flag', False)
        assert int(lp.process_lc(t, mag, err)['STATS_ff_MEAN_1']) == 0

    def test_set_inlist_value_unwraps_numpy_scalars(self):
        """numpy.float64(42.0) is unwrapped via .item() to Python float."""
        t = np.linspace(0, 1, 10)
        mag = np.full(10, 10.0); err = np.full(10, 0.01)
        lp = LibPipeline(['-inlistvars', 'v:0:double:0.0',
                          '-expr', 'vv=v', '-stats', 'vv', 'mean', '-oneline'])
        lp.set_inlist_value('v', np.float64(42.0))
        assert math.isclose(float(lp.process_lc(t, mag, err)['STATS_vv_MEAN_1']),
                            42.0, rel_tol=1e-12)

    def test_set_inlist_value_missing_raises(self):
        """Setting a name that wasn't declared via -inlistvars raises."""
        lp = LibPipeline(['-rms', '-oneline'])
        with pytest.raises(RuntimeError, match="no INLIST variable named"):
            lp.set_inlist_value('does_not_exist', 1.0)

    def test_set_inlist_value_type_mismatch_raises(self):
        """Setting a wrong-typed value raises."""
        lp = LibPipeline(['-inlistvars', 'v:0:double:0.0', '-rms', '-oneline'])
        with pytest.raises(RuntimeError, match="datatype does not match"):
            lp.set_inlist_value('v', "hello")


# ---------------------------------------------------------------------------
# Pipeline library-mode dispatch tests
# ---------------------------------------------------------------------------

class TestPipelineLibraryDispatch:

    @needs_library
    def test_pipeline_uses_library_when_available(self, simple_lc):
        """With library available, run() should call LibPipeline, not subprocess."""
        from pyvartools import pipeline as _pipeline_mod

        pipe = vt.Pipeline([cmd.rms()])
        subprocess_called = []

        def mock_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return "", ""

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
                pipe.run(simple_lc)

        assert not subprocess_called, (
            "_execute (subprocess) was called even though library mode is active"
        )

    @needs_binary
    def test_pipeline_falls_back_to_subprocess_when_library_unavailable(
        self, simple_lc
    ):
        """When library is unavailable, run() falls back to subprocess."""
        from pyvartools import pipeline as _pipeline_mod

        pipe = vt.Pipeline([cmd.rms()])
        subprocess_called = []

        original_execute = pipe._execute

        def tracking_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return original_execute(command, timeout=timeout, stdin_text=stdin_text)

        with patch.object(pipe, "_execute", side_effect=tracking_execute):
            with patch.object(_pipeline_mod, "_library_enabled", return_value=False):
                pipe.run(simple_lc)

        assert subprocess_called, (
            "_execute (subprocess) was NOT called even though library is disabled"
        )

    @needs_binary
    def test_pipeline_library_mode_env_var_disable(self, simple_lc):
        """VARTOOLS_USE_LIBRARY=0 forces subprocess mode."""
        from pyvartools import pipeline as _pipeline_mod

        subprocess_called = []
        pipe = vt.Pipeline([cmd.rms()])
        original_execute = pipe._execute

        def tracking_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return original_execute(command, timeout=timeout, stdin_text=stdin_text)

        with patch.dict(os.environ, {"VARTOOLS_USE_LIBRARY": "0"}):
            with patch.object(pipe, "_execute", side_effect=tracking_execute):
                result = pipe.run(simple_lc)

        assert subprocess_called, (
            "VARTOOLS_USE_LIBRARY=0 should force subprocess mode"
        )
        assert math.isfinite(float(result.vars["RMS_0"]))

    @needs_library
    @needs_binary
    def test_pipeline_library_run_produces_same_result_as_subprocess(
        self, sinusoidal_lc
    ):
        """Library mode and subprocess mode yield identical var."""
        from pyvartools import pipeline as _pipeline_mod

        pipe = vt.Pipeline([cmd.rms()])

        # Library mode
        with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
            result_lib = pipe.run(sinusoidal_lc)

        # Subprocess mode
        with patch.object(_pipeline_mod, "_library_enabled", return_value=False):
            result_sub = pipe.run(sinusoidal_lc)

        for key in result_sub.vars.index:
            if key == "Name":
                continue
            v_lib = float(result_lib.vars[key])
            v_sub = float(result_sub.vars[key])
            assert abs(v_lib - v_sub) < 1e-6, (
                f"{key}: library={v_lib} vs subprocess={v_sub}"
            )

    @needs_library
    def test_pipeline_library_run_batch_uses_library_per_lc(self):
        """run_batch() in library mode calls LibPipeline, not subprocess."""
        from pyvartools import pipeline as _pipeline_mod

        lcs = make_lcs(n=3)
        pipe = vt.Pipeline([cmd.rms()])
        subprocess_called = []

        def mock_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return "", ""

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
                batch = pipe.run_batch(lcs)

        assert not subprocess_called, (
            "subprocess was called for run_batch() in library mode"
        )
        assert len(batch.vars) == len(lcs)

    @needs_library
    @needs_binary
    def test_pipeline_library_run_batch_nthreads_falls_back_to_subprocess(self):
        """run_batch(nthreads>1) must use subprocess (library mode is single-threaded)."""
        from pyvartools import pipeline as _pipeline_mod

        lcs = make_lcs(n=3)
        pipe = vt.Pipeline([cmd.rms()])
        subprocess_called = []

        original_execute = pipe._execute

        def tracking_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return original_execute(command, timeout=timeout, stdin_text=stdin_text)

        with patch.object(pipe, "_execute", side_effect=tracking_execute):
            with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
                batch = pipe.run_batch(lcs, nthreads=2)

        assert subprocess_called, (
            "nthreads>1 should fall back to subprocess even when library is available"
        )
        assert len(batch.vars) == len(lcs)

    @needs_library
    @needs_binary
    @pytest.mark.parametrize("opt_kwargs", [
        {"randseed": 7},
        {"jdtol": 0.001},
        {"skipmissing": True},      # no-op in library batch (no list file)
        {"matchstringid": True},    # no-op in library batch (no list file)
        {"randseed": 7, "jdtol": 0.001,
         "skipmissing": True, "matchstringid": True},
    ])
    def test_pipeline_library_run_batch_global_opts_use_library(self, opt_kwargs):
        """run_batch(randseed=, skipmissing=, jdtol=, matchstringid=) goes
        through library mode after commit 2 of the library-batch widening
        plan.  skipmissing / matchstringid are list-file-only and no-op in
        library mode but accepted for API uniformity.
        """
        from pyvartools import pipeline as _pipeline_mod

        lcs = make_lcs(n=3)
        pipe = vt.Pipeline([cmd.rms()])
        subprocess_called = []

        def mock_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return "", ""

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
                batch = pipe.run_batch(lcs, **opt_kwargs)

        assert not subprocess_called, (
            f"subprocess was called for run_batch({opt_kwargs!r})"
        )
        assert len(batch.vars) == len(lcs)
        # rms is deterministic so the values must match across calls regardless
        # of which global opt was set; sanity-check that the column landed.
        assert "RMS_0" in batch.vars.columns

    @needs_library
    @needs_binary
    def test_pipeline_library_run_batch_global_opts_value_parity(self):
        """Library batch with global opts produces the same vars values
        as subprocess.  Uses randseed (the only global opt that observably
        changes output for a deterministic command like rms is none, but
        we still verify they all run cleanly and agree)."""
        from pyvartools import pipeline as _pipeline_mod

        lcs = make_lcs(n=3)
        pipe = vt.Pipeline([cmd.rms()])

        with patch.object(_pipeline_mod, "_library_enabled", return_value=False):
            sub = pipe.run_batch(lcs, randseed=42, jdtol=0.0001)
        # New Pipeline so the cached _lib_pipeline doesn't carry over the argv.
        pipe = vt.Pipeline([cmd.rms()])
        with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
            lib = pipe.run_batch(lcs, randseed=42, jdtol=0.0001)

        for col in sub.vars.columns:
            if col == "Name":
                continue
            for i in range(len(lcs)):
                assert math.isclose(float(sub.vars[col].iloc[i]),
                                    float(lib.vars[col].iloc[i]),
                                    rel_tol=1e-12, abs_tol=1e-12), (
                    f"col={col} lc={i}: lib={lib.vars[col].iloc[i]} "
                    f"vs sub={sub.vars[col].iloc[i]}"
                )

    @needs_library
    @needs_binary
    def test_pipeline_library_run_batch_perpoint_vars_uses_library(self):
        """run_batch(perpoint_vars=...) goes through library mode after
        commit 3 of the library-batch widening plan.  The init-expression
        clause is appended to -inputlcformat at LibPipeline init time."""
        from pyvartools import pipeline as _pipeline_mod

        lcs = make_lcs(n=3)
        pipe = vt.Pipeline([cmd.stats("phase", "max")])
        subprocess_called = []

        def mock_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return "", ""

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
                batch = pipe.run_batch(
                    lcs,
                    perpoint_vars={"phase": vt.PerPointVar(
                        type="double", init="NR")})

        assert not subprocess_called, (
            "subprocess was called for run_batch(perpoint_vars=...)"
        )
        assert len(batch.vars) == len(lcs)

    @needs_library
    @needs_binary
    def test_pipeline_library_run_batch_perpoint_vars_value_parity(self):
        """Library batch with perpoint_vars matches subprocess output."""
        from pyvartools import pipeline as _pipeline_mod

        lcs = make_lcs(n=3)
        ppv = {"phase": vt.PerPointVar(type="double", init="NR")}

        def make_pipe():
            return vt.Pipeline([cmd.stats("phase", "max")])

        with patch.object(_pipeline_mod, "_library_enabled", return_value=False):
            sub = make_pipe().run_batch(lcs, perpoint_vars=ppv)
        with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
            lib = make_pipe().run_batch(lcs, perpoint_vars=ppv)

        for col in sub.vars.columns:
            if col == "Name":
                continue
            for i in range(len(lcs)):
                assert math.isclose(float(sub.vars[col].iloc[i]),
                                    float(lib.vars[col].iloc[i]),
                                    rel_tol=1e-12, abs_tol=1e-12), (
                    f"col={col} lc={i}: lib={lib.vars[col].iloc[i]} "
                    f"vs sub={sub.vars[col].iloc[i]}"
                )

    @needs_library
    @needs_binary
    def test_pipeline_library_run_batch_command_offset_uses_library(self):
        """run_batch(_command_offset=N) goes through library mode.

        Before commit 1 of the library-batch widening plan, _command_offset != 0
        forced subprocess.  Now it's just a -columnsuffix shift threaded through
        _commands_to_argv.
        """
        from pyvartools import pipeline as _pipeline_mod

        lcs = make_lcs(n=3)
        pipe = vt.Pipeline([cmd.rms()])
        subprocess_called = []

        def mock_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return "", ""

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
                batch = pipe.run_batch(lcs, _command_offset=2)

        assert not subprocess_called, (
            "subprocess was called for run_batch(_command_offset=2)"
        )
        assert len(batch.vars) == len(lcs)
        assert "RMS_2" in batch.vars.columns
        assert "RMS_0" not in batch.vars.columns

    @needs_library
    @needs_binary
    def test_pipeline_library_capture_lc(self, simple_lc):
        """run(lc, capture_lc=True) in library mode returns a LightCurve."""
        from pyvartools import pipeline as _pipeline_mod

        pipe = vt.Pipeline([cmd.clip(5.0)])
        with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
            result = pipe.run(simple_lc, capture_lc=True)

        assert result.lc is not None
        assert isinstance(result.lc, vt.LightCurve)
        assert len(result.lc) == len(simple_lc)

    @needs_library
    @needs_binary
    def test_pipeline_library_cmd_o_capture(self, simple_lc):
        """cmd.o(capture=True) in library mode returns result.files['o']."""
        from pyvartools import pipeline as _pipeline_mod

        pipe = vt.Pipeline([cmd.o(capture=True)])
        with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
            result = pipe.run(simple_lc)

        assert "o" in result.files
        assert isinstance(result.files["o"], vt.LightCurve)

    @needs_library
    @needs_binary
    def test_pipeline_library_save_periodogram(self, sinusoidal_lc):
        """save_periodogram=True in library mode returns a DataFrame."""
        from pyvartools import pipeline as _pipeline_mod

        pipe = vt.Pipeline([
            cmd.LS(0.5, 5.0, 0.01, npeaks=1, save_periodogram=True),
        ])
        with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
            result = pipe.run(sinusoidal_lc)

        key = "LS_periodogram_0"
        assert key in result.files, (
            f"Expected '{key}' in result.files; got {list(result.files.keys())}"
        )
        assert isinstance(result.files[key], pd.DataFrame)
        assert len(result.files[key]) > 0


# ---------------------------------------------------------------------------
# Step B: cmd.o(outname=...) in single-LC library mode
# ---------------------------------------------------------------------------

@needs_library
@needs_binary
class TestSingleLCOutputInLibraryMode:
    """``cmd.o(outname=PATH)`` runs in library mode for single-LC ``run(lc)``,
    not subprocess.  This is the first slice of the larger library-mode-with-
    file-output roadmap (Steps B / C3 / D)."""

    def test_outname_uses_library_mode(self, simple_lc, tmp_path):
        """A pipeline ending in ``cmd.o(outname=...)`` does NOT spawn a
        subprocess when library mode is available."""
        from pyvartools import pipeline as _pipeline_mod

        out = tmp_path / "out.lc"
        pipe = vt.Pipeline([cmd.clip(sigclip=5.0), cmd.o(outname=str(out))])
        subprocess_called = []

        def mock_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return "", ""

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
                pipe.run(simple_lc)

        assert not subprocess_called, (
            "_execute (subprocess) was called even though library mode "
            "should handle cmd.o(outname=...) in single-LC runs"
        )

    def test_outname_writes_file(self, simple_lc, tmp_path):
        """The library-mode write actually produces the file on disk."""
        from pyvartools import pipeline as _pipeline_mod
        out = tmp_path / "out.lc"
        with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
            vt.Pipeline([cmd.clip(sigclip=5.0),
                         cmd.o(outname=str(out))]).run(simple_lc)
        assert out.is_file()
        assert out.stat().st_size > 0
        # Re-read the written file and verify it has the expected number
        # of points (clip removes none for a clean sinusoid).
        df = pd.read_csv(str(out), sep=r"\s+", comment="#", header=None)
        assert len(df) == len(simple_lc)

    def test_outname_byte_parity_with_subprocess(self, simple_lc, tmp_path):
        """Library-mode and subprocess-mode writes are byte-equal."""
        from pyvartools import pipeline as _pipeline_mod

        lib_path = tmp_path / "lib.lc"
        sub_path = tmp_path / "sub.lc"

        with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
            vt.Pipeline([cmd.clip(sigclip=5.0),
                         cmd.o(outname=str(lib_path))]).run(simple_lc)

        with patch.object(_pipeline_mod, "_library_enabled", return_value=False):
            vt.Pipeline([cmd.clip(sigclip=5.0),
                         cmd.o(outname=str(sub_path))]).run(simple_lc)

        assert lib_path.read_bytes() == sub_path.read_bytes(), (
            "Library and subprocess paths produced different bytes for "
            "cmd.o(outname=...).  Investigate column ordering, headers, "
            "or numeric formatting drift."
        )

    def test_outname_with_capture_uses_library_mode(self, simple_lc, tmp_path):
        """``cmd.o(outname=path, capture=True)`` runs in library mode via
        the new ``capture_id`` keyword: vartools writes the file and
        also fills the in-memory slot in one call."""
        from pyvartools import pipeline as _pipeline_mod

        out = tmp_path / "out.lc"
        pipe = vt.Pipeline([cmd.clip(sigclip=5.0),
                            cmd.o(outname=str(out), capture=True, key="x")])
        subprocess_called = []

        def mock_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return "", ""

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch.object(_pipeline_mod, "_library_enabled",
                              return_value=True):
                result = pipe.run(simple_lc)
        assert not subprocess_called, (
            "capture+outname must run in library mode now (Step D follow-up)"
        )
        assert out.is_file()
        assert out.stat().st_size > 0
        assert "x" in result.files
        assert isinstance(result.files["x"], vt.LightCurve)

    def test_save_periodogram_with_outname_uses_library_mode(
            self, sinusoidal_lc, tmp_path):
        """A pipeline mixing cmd.o(outname=...) and save_periodogram=True
        runs in library mode now (the per-Pipeline tmpdir holds the
        save_* output files)."""
        from pyvartools import pipeline as _pipeline_mod
        out = tmp_path / "out.lc"
        pipe = vt.Pipeline([
            cmd.LS(0.5, 5.0, 0.01, npeaks=1, save_periodogram=True),
            cmd.o(outname=str(out)),
        ])
        subprocess_called = []

        def mock_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return "", ""

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch.object(_pipeline_mod, "_library_enabled",
                              return_value=True):
                result = pipe.run(sinusoidal_lc)
        assert not subprocess_called, (
            "save_*=True + cmd.o(outname=...) should now run in library mode"
        )
        assert out.is_file()
        assert "LS_periodogram_0" in result.files


# ---------------------------------------------------------------------------
# Step C3: cmd.o(outdir=...) in run_batch via the new forceoutdirmode keyword
# ---------------------------------------------------------------------------

def _make_lc_lib(n=64, period=2.0, name="lc"):
    t = np.linspace(0, 10, n)
    mag = 10.0 + 0.1 * np.sin(2 * np.pi * t / period)
    err = np.full(n, 0.01)
    lc = vt.LightCurve.from_arrays(t, mag, err)
    lc.name = name
    return lc


@needs_library
@needs_binary
class TestBatchOutputInLibraryMode:
    """``cmd.o(outdir=...)`` now runs through library mode in
    ``run_batch``, with vartools' new ``forceoutdirmode`` keyword
    flipping its single-file output writer into directory-naming mode."""

    def test_run_batch_uses_library_mode(self, tmp_path):
        from pyvartools import pipeline as _pipeline_mod

        lcs = [_make_lc_lib(name="alpha"),
               _make_lc_lib(name="beta"),
               _make_lc_lib(name="gamma")]
        outdir = tmp_path / "out"
        outdir.mkdir()
        pipe = vt.Pipeline([cmd.clip(sigclip=5.0),
                            cmd.o(outdir=str(outdir))])
        subprocess_called = []
        original_execute = pipe._execute

        def tracking_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return original_execute(command, timeout=timeout, stdin_text=stdin_text)

        with patch.object(pipe, "_execute", side_effect=tracking_execute):
            with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
                pipe.run_batch(lcs)
        assert not subprocess_called, (
            "_execute (subprocess) was called even though library mode "
            "should handle cmd.o(outdir=...) in run_batch"
        )
        files = sorted(p.name for p in outdir.iterdir())
        assert files == ["alpha", "beta", "gamma"]

    def test_run_batch_byte_parity_with_subprocess(self, tmp_path):
        """Library and subprocess paths produce byte-equal LC files."""
        from pyvartools import pipeline as _pipeline_mod

        lcs = [_make_lc_lib(name=n, period=p)
               for n, p in [("alpha", 1.5), ("beta", 2.5), ("gamma", 3.5)]]

        lib_dir = tmp_path / "lib"; lib_dir.mkdir()
        sub_dir = tmp_path / "sub"; sub_dir.mkdir()

        with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
            vt.Pipeline([cmd.clip(sigclip=5.0),
                         cmd.o(outdir=str(lib_dir))]).run_batch(lcs)
        with patch.object(_pipeline_mod, "_library_enabled", return_value=False):
            vt.Pipeline([cmd.clip(sigclip=5.0),
                         cmd.o(outdir=str(sub_dir))]).run_batch(lcs)

        for n in ("alpha", "beta", "gamma"):
            assert (lib_dir / n).read_bytes() == (sub_dir / n).read_bytes(), (
                f"{n}: library and subprocess paths produced different bytes."
            )

    def test_run_batch_disambiguates_duplicate_names(self, tmp_path):
        from pyvartools import pipeline as _pipeline_mod

        lc1 = _make_lc_lib(name="dup", period=1.5)
        lc2 = _make_lc_lib(name="dup", period=2.5)
        outdir = tmp_path / "out"; outdir.mkdir()
        with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
            vt.Pipeline([cmd.clip(sigclip=5.0),
                         cmd.o(outdir=str(outdir))]).run_batch([lc1, lc2])
        files = sorted(p.name for p in outdir.iterdir())
        assert files == ["dup", "dup_1"]
        for f in files:
            assert (outdir / f).stat().st_size > 0

    def test_run_batch_with_outname_still_subprocess(self, tmp_path):
        """``cmd.o(outname=...)`` is illegal in batch mode (the wrapper
        raises in the build step), so the gate stays in subprocess.
        This is here just to lock the boundary."""
        from pyvartools import pipeline as _pipeline_mod
        lcs = [_make_lc_lib(name="a"), _make_lc_lib(name="b")]
        out = tmp_path / "x"
        pipe = vt.Pipeline([cmd.clip(sigclip=5.0),
                            cmd.o(outname=str(out))])
        with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
            with pytest.raises(Exception):
                pipe.run_batch(lcs)

    def test_run_batch_capture_with_outdir_uses_library_mode(self, tmp_path):
        """``capture=True`` + ``outdir=DIR`` runs in library mode and
        produces both per-LC files on disk and per-LC captures in
        ``batch.files[key]`` (Step D follow-up)."""
        from pyvartools import pipeline as _pipeline_mod
        lcs = [_make_lc_lib(name="a"), _make_lc_lib(name="b")]
        outdir = tmp_path / "out"; outdir.mkdir()
        pipe = vt.Pipeline([cmd.clip(sigclip=5.0),
                            cmd.o(outdir=str(outdir),
                                  capture=True, key="cap")])
        subprocess_called = []

        def mock_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return "", ""

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch.object(_pipeline_mod, "_library_enabled",
                              return_value=True):
                batch = pipe.run_batch(lcs)
        assert not subprocess_called, (
            "capture+outdir must run in library mode now (Step D follow-up)"
        )
        # Both files written, both captures populated.
        files = sorted(p.name for p in outdir.iterdir())
        assert files == ["a", "b"]
        assert "cap" in batch.files
        assert len(batch.files["cap"]) == 2
        for c in batch.files["cap"]:
            assert isinstance(c, vt.LightCurve)


def test_forceoutdirmode_emitted_for_library_batch():
    """Pure-Python: ``cmd.o(outdir=...)._to_cli_args_for_mode('library_batch')``
    appends the new ``forceoutdirmode`` keyword and uses outdir as the path."""
    c = cmd.o(outdir="/tmp/dir")
    args = c._to_cli_args_for_mode("library_batch")
    assert args[0] == "-o"
    assert args[1] == "/tmp/dir"
    assert "forceoutdirmode" in args
    # Not emitted for "list" or "single" — those go through other branches.
    list_args = cmd.o(outdir="/tmp/dir")._to_cli_args_for_mode("list")
    assert "forceoutdirmode" not in list_args
    single_args = cmd.o(outname="/tmp/file.lc")._to_cli_args_for_mode("single")
    assert "forceoutdirmode" not in single_args


# ---------------------------------------------------------------------------
# Step D: cmd.o(capture=True) in library mode -- in-memory only
# ---------------------------------------------------------------------------

def test_capture_argv_emits_capture_keyword():
    """Pure-Python: ``cmd.o(capture=True, key="X")._to_cli_args_for_mode(
    'library_*')`` emits ``-o X capture`` and nothing else."""
    for mode in ("library_single", "library_batch"):
        c = cmd.o(capture=True, key="X")
        args = c._to_cli_args_for_mode(mode)
        assert args == ["-o", "X", "capture"], (
            f"mode={mode}: expected pure-capture argv, got {args}"
        )


def test_capture_argv_emits_capture_id_for_combined_mode():
    """Pure-Python: ``cmd.o(capture=True, key="X", outname=PATH)`` in
    library_single emits the path-based form *plus* ``capture_id X`` so
    vartools writes the file AND fills the slot.  Same shape for
    outdir+capture in library_batch."""
    c = cmd.o(capture=True, key="X", outname="/tmp/file.lc")
    args = c._to_cli_args_for_mode("library_single")
    assert args[0] == "-o"
    assert args[1] == "/tmp/file.lc"     # path, not the key
    assert "capture" not in args         # bare "capture" not emitted
    assert "capture_id" in args
    # capture_id is followed immediately by the key.
    idx = args.index("capture_id")
    assert args[idx + 1] == "X"

    c = cmd.o(capture=True, key="Y", outdir="/tmp/dir")
    args = c._to_cli_args_for_mode("library_batch")
    assert args[1] == "/tmp/dir"
    assert "forceoutdirmode" in args
    assert "capture_id" in args
    assert args[args.index("capture_id") + 1] == "Y"

    # In subprocess modes (single/list) capture_id is NOT emitted —
    # those go through the existing tmp-file capture path.
    c = cmd.o(capture=True, key="X", outname="/tmp/file.lc")
    args = c._to_cli_args_for_mode("single")
    assert "capture_id" not in args, (
        "subprocess single mode shouldn't emit capture_id"
    )


@needs_library
@needs_binary
class TestCaptureInLibraryMode:
    """``cmd.o(capture=True)`` with no disk path runs entirely in memory
    via the new C-side capture buffer (Step D).  No tmp directory is
    allocated; the captured arrays are pulled out via
    LibPipeline.read_capture()."""

    def test_single_lc_capture_uses_library_mode(self, simple_lc):
        from pyvartools import pipeline as _pipeline_mod
        pipe = vt.Pipeline([cmd.clip(sigclip=5.0),
                            cmd.o(capture=True, key="clipped")])
        subprocess_called = []

        def mock_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return "", ""

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch.object(_pipeline_mod, "_library_enabled",
                              return_value=True):
                result = pipe.run(simple_lc)
        assert not subprocess_called, (
            "cmd.o(capture=True) should run in library mode, not subprocess"
        )
        assert "clipped" in result.files
        captured = result.files["clipped"]
        # The clipped LC should have the same number of points as the
        # input (the simple_lc fixture is a clean sinusoid; clip 5σ
        # removes nothing) and the standard t/mag/err columns.
        assert isinstance(captured, vt.LightCurve)
        assert "t" in captured._df.columns
        assert "mag" in captured._df.columns
        assert "err" in captured._df.columns

    def test_capture_does_not_allocate_tmp_dirs(self, simple_lc, tmp_path):
        """No ``_o_cap_*`` tmp dirs are created in library-capture mode."""
        import tempfile
        # Snapshot the current set of /tmp entries so we can compare.
        before = set(os.listdir(tempfile.gettempdir()))
        from pyvartools import pipeline as _pipeline_mod
        with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
            for _ in range(3):
                vt.Pipeline([cmd.clip(sigclip=5.0),
                             cmd.o(capture=True, key="x")]).run(simple_lc)
        after = set(os.listdir(tempfile.gettempdir()))
        new = after - before
        # The temporary working directory pyvartools allocates *would*
        # show up here for subprocess paths, plus _o_cap_* subdirs.  In
        # library-capture mode neither should appear.
        leaked = [n for n in new if "_o_cap_" in n or "tmp" in n.lower()]
        assert not leaked, f"tmp leftovers: {leaked}"

    def test_single_lc_capture_value_parity_with_subprocess(self, simple_lc):
        """Library and subprocess paths produce equal captured LCs."""
        from pyvartools import pipeline as _pipeline_mod
        with patch.object(_pipeline_mod, "_library_enabled",
                          return_value=True):
            res_lib = vt.Pipeline([cmd.clip(sigclip=5.0),
                                   cmd.o(capture=True, key="x")]).run(simple_lc)
        with patch.object(_pipeline_mod, "_library_enabled",
                          return_value=False):
            res_sub = vt.Pipeline([cmd.clip(sigclip=5.0),
                                   cmd.o(capture=True, key="x")]).run(simple_lc)
        df_lib = res_lib.files["x"]._df
        df_sub = res_sub.files["x"]._df
        assert list(df_lib.columns) == list(df_sub.columns)
        assert len(df_lib) == len(df_sub)
        # Library mode preserves full %.17g precision; subprocess goes
        # through %.17g spill -> ASCII parse, so values should match
        # exactly.  Allow a tiny relative tolerance to be safe.
        for col in df_lib.columns:
            np.testing.assert_allclose(
                df_lib[col].values, df_sub[col].values,
                rtol=1e-12, atol=0,
                err_msg=f"column {col} drifted between library and subprocess"
            )

    def test_multi_capture_distinct_keys(self, sinusoidal_lc):
        """Two cmd.o(capture=True) at different points produce distinct
        snapshots, both captured in a single library call."""
        from pyvartools import pipeline as _pipeline_mod
        pipe = vt.Pipeline([
            cmd.o(capture=True, key="before"),
            cmd.clip(sigclip=2.0),     # aggressive clip — removes points
            cmd.o(capture=True, key="after"),
        ])
        subprocess_called = []

        def mock_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return "", ""

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch.object(_pipeline_mod, "_library_enabled",
                              return_value=True):
                result = pipe.run(sinusoidal_lc)
        assert not subprocess_called, "multi-capture must stay in library mode"
        assert "before" in result.files
        assert "after" in result.files
        assert len(result.files["before"]._df) >= len(result.files["after"]._df)

    @needs_binary
    def test_run_batch_capture_lc_uses_library_mode(self):
        """``run_batch(capture_lc=True)`` goes through library mode after
        commit 5 of the library-batch widening plan, instead of falling
        through to subprocess."""
        from pyvartools import pipeline as _pipeline_mod

        lcs = make_lcs(n=3)
        pipe = vt.Pipeline([cmd.clip(sigclip=5.0)])
        subprocess_called = []

        def mock_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return "", ""

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch.object(_pipeline_mod, "_library_enabled",
                              return_value=True):
                batch = pipe.run_batch(lcs, capture_lc=True)

        assert not subprocess_called, (
            "run_batch(capture_lc=True) should run in library mode"
        )
        assert len(batch.lcs) == 3
        for i, out_lc in enumerate(batch.lcs):
            assert isinstance(out_lc, vt.LightCurve)
            assert out_lc.name == lcs[i].name
            assert "t" in out_lc._df.columns
            assert "mag" in out_lc._df.columns
            assert "err" in out_lc._df.columns

    @needs_binary
    def test_run_batch_stats_file_uses_library_mode(self, tmp_path):
        """``run_batch(stats_file=PATH)`` goes through library mode after
        commit 6 of the library-batch widening plan."""
        from pyvartools import pipeline as _pipeline_mod

        lcs = make_lcs(n=3)
        pipe = vt.Pipeline([cmd.rms()])
        subprocess_called = []

        def mock_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return "", ""

        sf = tmp_path / "out.stats"
        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch.object(_pipeline_mod, "_library_enabled",
                              return_value=True):
                pipe.run_batch(lcs, stats_file=str(sf))

        assert not subprocess_called, (
            "run_batch(stats_file=...) should run in library mode"
        )
        # File has 1 header + 3 data rows.
        contents = sf.read_text()
        lines = contents.strip().split("\n")
        assert lines[0].startswith("#"), "header line missing '#' prefix"
        assert len(lines) == 1 + len(lcs)
        # Library batch writes lc.name in the Name column (not a temp path).
        for i, lc in enumerate(lcs):
            assert lc.name in lines[1 + i], (
                f"row {i} should contain lc.name {lc.name!r}, got {lines[1+i]!r}"
            )

    @needs_binary
    def test_run_batch_stats_file_resume_falls_back_to_subprocess(self, tmp_path):
        """``run_batch(stats_file=PATH, resume=True)`` routes through
        subprocess until the post-Tier-C byte-parity work lands.  The
        library-batch path doesn't emit a seq column so _resume_partition
        would have no key to match completed rows."""
        from pyvartools import pipeline as _pipeline_mod

        lcs = make_lcs(n=3)
        sf = tmp_path / "out.stats"
        # First run writes a partial subprocess stats_file (with seq column).
        with patch.object(_pipeline_mod, "_library_enabled",
                          return_value=False):
            vt.Pipeline([cmd.rms()]).run_batch(lcs[:1], stats_file=str(sf))

        # Resume run with library enabled.  The seq_col_remap gate keeps
        # routing to subprocess; verify by inspecting the final file —
        # subprocess output contains the internal _vtpy_seq_ column,
        # library output does not.
        with patch.object(_pipeline_mod, "_library_enabled",
                          return_value=True):
            vt.Pipeline([cmd.rms()]).run_batch(
                lcs, stats_file=str(sf), resume=True)

        header = sf.read_text().split("\n", 1)[0]
        assert "_vtpy_seq_" in header, (
            "resume=True should force subprocess; the persisted file should "
            "still contain the _vtpy_seq_ column.  Got header: " + header
        )

    @needs_binary
    def test_run_batch_capture_lc_parity_with_subprocess(self):
        """Library batch with capture_lc=True produces identical LC arrays
        to the subprocess path.  Catches any divergence in column ordering,
        numeric formatting, or array dtype between
        ``process_lc_capture`` (library) and the subprocess capture path."""
        from pyvartools import pipeline as _pipeline_mod

        lcs = make_lcs(n=3)

        def make_pipe():
            return vt.Pipeline([cmd.clip(sigclip=5.0)])

        with patch.object(_pipeline_mod, "_library_enabled",
                          return_value=False):
            sub = make_pipe().run_batch(lcs, capture_lc=True)
        with patch.object(_pipeline_mod, "_library_enabled",
                          return_value=True):
            lib = make_pipe().run_batch(lcs, capture_lc=True)

        for i in range(len(lcs)):
            lib_lc = lib.lcs[i]
            sub_lc = sub.lcs[i]
            assert lib_lc is not None and sub_lc is not None
            assert list(lib_lc._df.columns) == list(sub_lc._df.columns), (
                f"LC {i}: column order mismatch "
                f"(lib={list(lib_lc._df.columns)}, sub={list(sub_lc._df.columns)})"
            )
            for col in lib_lc._df.columns:
                np.testing.assert_allclose(
                    lib_lc._df[col].values, sub_lc._df[col].values,
                    rtol=1e-12, atol=1e-12,
                    err_msg=f"LC {i} col {col!r}: lib vs sub mismatch")

    def test_batch_capture_per_lc(self):
        """``run_batch`` with cmd.o(capture=True) returns ``files[key]`` as
        a list of LightCurves, one per input LC."""
        from pyvartools import pipeline as _pipeline_mod
        n = 100
        lcs = []
        for i in range(3):
            t = np.linspace(0, 10, n)
            mag = 10.0 + 0.05 * (i + 1) * np.sin(2 * np.pi * t / 2.0)
            err = np.full(n, 0.01)
            lc = vt.LightCurve.from_arrays(t, mag, err)
            lc.name = f"lc{i}"
            lcs.append(lc)
        with patch.object(_pipeline_mod, "_library_enabled",
                          return_value=True):
            batch = vt.Pipeline([cmd.clip(sigclip=5.0),
                                 cmd.o(capture=True, key="cap")]).run_batch(lcs)
        assert "cap" in batch.files
        captures = batch.files["cap"]
        assert len(captures) == 3
        for i, captured_lc in enumerate(captures):
            assert isinstance(captured_lc, vt.LightCurve)
            assert captured_lc.name == f"lc{i}"
            assert "t" in captured_lc._df.columns

    def test_capture_with_outname_writes_and_captures(self, simple_lc,
                                                      tmp_path):
        """``cmd.o(outname=path, capture=True)`` in library mode emits
        the new ``capture_id`` keyword: vartools writes the file *and*
        snapshots the post-write LC into ``result.files[key]`` in one
        call, byte-equal to the subprocess capture path."""
        from pyvartools import pipeline as _pipeline_mod
        out_lib = tmp_path / "lib.lc"
        out_sub = tmp_path / "sub.lc"
        with patch.object(_pipeline_mod, "_library_enabled",
                          return_value=True):
            res_lib = vt.Pipeline([cmd.clip(sigclip=5.0),
                                   cmd.o(outname=str(out_lib),
                                         capture=True,
                                         key="cap")]).run(simple_lc)
        with patch.object(_pipeline_mod, "_library_enabled",
                          return_value=False):
            res_sub = vt.Pipeline([cmd.clip(sigclip=5.0),
                                   cmd.o(outname=str(out_sub),
                                         capture=True,
                                         key="cap")]).run(simple_lc)
        # File written in both modes.
        assert out_lib.is_file()
        assert out_sub.is_file()
        # Library file matches subprocess file byte-for-byte.
        assert out_lib.read_bytes() == out_sub.read_bytes()
        # Captures populated in both, value-equal.
        df_lib = res_lib.files["cap"]._df
        df_sub = res_sub.files["cap"]._df
        for col in df_lib.columns:
            np.testing.assert_allclose(
                df_lib[col].values, df_sub[col].values,
                rtol=1e-12, atol=0,
                err_msg=f"library vs subprocess drift in column {col}"
            )


# ---------------------------------------------------------------------------
# CLI thread-safety: ``-o ID capture`` is a silent no-op outside library mode,
# safe under -parallel because each thread reads p->captured == NULL and
# early-returns.  Locks in the contract documented at the top of
# vartools_init_pipeline's pre-walk block.
# ---------------------------------------------------------------------------

@needs_binary
class TestCliCaptureKeywordSafety:
    def test_capture_keyword_silent_in_cli_single(self):
        """Standalone vartools binary with `-o ID capture` runs successfully
        and produces normal stats; the capture is silently discarded."""
        import subprocess
        proc = subprocess.run(
            [vt.get_binary(), "-i", "-",
             "-o", "discard", "capture",
             "-rms", "-oneline"],
            input="1.0 10.0 0.01\n2.0 10.1 0.01\n3.0 9.9 0.01\n",
            capture_output=True, text=True, timeout=15,
        )
        assert proc.returncode == 0, (
            f"vartools failed: stderr={proc.stderr!r}"
        )
        # Stats are emitted normally (capture didn't displace the -rms output).
        assert "RMS" in proc.stdout

    def test_capture_keyword_silent_in_cli_parallel(self, tmp_path):
        """Same CLI invocation but with -l listfile and -parallel 4: still
        a silent no-op, no race / crash from the multi-thread DoOutputLight
        Curve calls.  Pre-Step-D the new ProgramData fields were
        uninitialised on the CLI stack and could have been a non-NULL
        pointer; main.c now zero-inits them."""
        import subprocess
        # Build a tiny list of 6 quick LCs.
        lcs = []
        for i in range(6):
            p = tmp_path / f"lc{i}.txt"
            with open(p, "w") as f:
                for j in range(50):
                    f.write(f"{j*0.1} {10.0 + 0.01*i + 0.001*j} 0.01\n")
            lcs.append(str(p))
        listfile = tmp_path / "list.txt"
        listfile.write_text("\n".join(lcs) + "\n")

        proc = subprocess.run(
            [vt.get_binary(), "-l", str(listfile),
             "-o", "discard", "capture",
             "-rms", "-parallel", "4", "-oneline"],
            capture_output=True, text=True, timeout=30,
        )
        assert proc.returncode == 0, (
            f"vartools -parallel 4 with capture keyword failed: "
            f"stderr={proc.stderr!r}"
        )
        # Six rows of stats output, one per LC.
        rms_lines = [
            line for line in proc.stdout.splitlines() if "RMS" in line
        ]
        assert len(rms_lines) >= 6


# ---------------------------------------------------------------------------
# save_*=True outputs (save_periodogram et al.) in library mode
# ---------------------------------------------------------------------------

@needs_library
@needs_binary
class TestSaveOutputsInLibraryMode:
    """``save_*=True`` outputs (save_periodogram, save_model, ...) run in
    library mode by routing the C-side writers through a per-Pipeline
    tmpdir.  No subprocess fork; the existing subprocess parsers handle
    the file-readback path unchanged."""

    def test_save_periodogram_single_lc(self, sinusoidal_lc):
        from pyvartools import pipeline as _pipeline_mod

        pipe = vt.Pipeline([
            cmd.LS(0.5, 5.0, 0.01, npeaks=1, save_periodogram=True),
        ])
        subprocess_called = []

        def mock_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return "", ""

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch.object(_pipeline_mod, "_library_enabled",
                              return_value=True):
                result = pipe.run(sinusoidal_lc)
        assert not subprocess_called, (
            "save_periodogram should now run in library mode"
        )
        key = "LS_periodogram_0"
        assert key in result.files
        df = result.files[key]
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0

    def test_save_periodogram_value_parity_with_subprocess(self, sinusoidal_lc):
        """Library and subprocess produce equal periodogram DataFrames."""
        from pyvartools import pipeline as _pipeline_mod

        pipe_lib = vt.Pipeline([
            cmd.LS(0.5, 5.0, 0.01, npeaks=1, save_periodogram=True),
        ])
        pipe_sub = vt.Pipeline([
            cmd.LS(0.5, 5.0, 0.01, npeaks=1, save_periodogram=True),
        ])
        with patch.object(_pipeline_mod, "_library_enabled", return_value=True):
            res_lib = pipe_lib.run(sinusoidal_lc)
        with patch.object(_pipeline_mod, "_library_enabled", return_value=False):
            res_sub = pipe_sub.run(sinusoidal_lc)
        df_lib = res_lib.files["LS_periodogram_0"]
        df_sub = res_sub.files["LS_periodogram_0"]
        assert list(df_lib.columns) == list(df_sub.columns)
        assert len(df_lib) == len(df_sub)
        for col in df_lib.columns:
            np.testing.assert_allclose(
                df_lib[col].values, df_sub[col].values,
                rtol=1e-10, atol=0,
                err_msg=f"library/subprocess drift in periodogram col {col}",
            )

    def test_save_periodogram_batch(self):
        """``run_batch`` accumulates per-LC save_periodogram DataFrames."""
        from pyvartools import pipeline as _pipeline_mod

        n = 200
        lcs = []
        for i in range(3):
            t = np.linspace(0, 100, n)
            mag = 10.0 + 0.05 * np.sin(2*np.pi*t/(2.0+i*0.5))
            err = np.full(n, 0.01)
            lc = vt.LightCurve.from_arrays(t, mag, err)
            lc.name = f"lc{i}"
            lcs.append(lc)

        pipe = vt.Pipeline([
            cmd.LS(0.5, 5.0, 0.01, npeaks=1, save_periodogram=True),
        ])
        subprocess_called = []

        def mock_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return "", ""

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch.object(_pipeline_mod, "_library_enabled",
                              return_value=True):
                batch = pipe.run_batch(lcs)
        assert not subprocess_called
        per_lc = batch.files["LS_periodogram_0"]
        assert len(per_lc) == 3
        for df in per_lc:
            assert isinstance(df, pd.DataFrame)
            assert len(df) > 0


# ---------------------------------------------------------------------------
# UserCommand / userlib extensions in library mode
# ---------------------------------------------------------------------------
# vartools' user-extension .so files (fastchi2, splinedetrend, jktebop,
# macula, magadd, hatpiflag, stitch, ftuneven) load via lt_dlopen during
# parsecommandline.  In library mode pyvartools loads libvartoolspipeline
# with RTLD_GLOBAL so its statically-linked GSL/CFITSIO/CSPICE symbols
# are visible to those extensions when they're dlopen'd.

@needs_library
@needs_binary
class TestUserCommandsInLibraryMode:
    """``cmd.fastchi2`` and friends now run in library mode without a
    subprocess fallback (Step E)."""

    def test_fastchi2_library_mode(self):
        """fastchi2 (which calls into libgsl) runs in library mode and
        returns its own output columns alongside any subsequent ones."""
        from pyvartools import pipeline as _pipeline_mod

        t = np.linspace(0, 30, 500)
        mag = 10.0 + 0.05 * np.sin(2 * np.pi * t / 2.0)
        err = np.full(500, 0.005)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="testlc")

        pipe = vt.Pipeline().fastchi2(Nharm=2, freqmin=0.1,
                                      freqmax=10.0).rms()
        subprocess_called = []
        original_execute = pipe._execute

        def tracking_execute(command, timeout=None, stdin_text=None):
            subprocess_called.append(True)
            return original_execute(command, timeout=timeout,
                                    stdin_text=stdin_text)

        with patch.object(pipe, "_execute", side_effect=tracking_execute):
            with patch.object(_pipeline_mod, "_library_enabled",
                              return_value=True):
                result = pipe.run(lc)
        assert not subprocess_called, (
            "fastchi2 should run in library mode -- libvartoolspipeline "
            "is loaded with RTLD_GLOBAL so libgsl symbols resolve when "
            "fastchi2.so is dlopen'd"
        )
        # Pipeline should have produced fastchi2's output columns.
        fastchi2_cols = [c for c in result.vars.index if c.startswith("Fastchi2")]
        assert fastchi2_cols, (
            f"Expected Fastchi2_* columns in result.vars; got "
            f"{list(result.vars.index)}"
        )
        # And the rms call after it.
        assert "RMS_1" in result.vars.index
        assert math.isfinite(float(result.vars["RMS_1"]))

    def test_fastchi2_value_parity_with_subprocess(self):
        """Library and subprocess paths produce equal Fastchi2 outputs."""
        from pyvartools import pipeline as _pipeline_mod

        t = np.linspace(0, 30, 500)
        mag = 10.0 + 0.05 * np.sin(2 * np.pi * t / 2.0)
        err = np.full(500, 0.005)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="testlc")

        with patch.object(_pipeline_mod, "_library_enabled",
                          return_value=True):
            res_lib = vt.Pipeline().fastchi2(
                Nharm=2, freqmin=0.1, freqmax=10.0).rms().run(lc)
        with patch.object(_pipeline_mod, "_library_enabled",
                          return_value=False):
            res_sub = vt.Pipeline().fastchi2(
                Nharm=2, freqmin=0.1, freqmax=10.0).rms().run(lc)

        # Same column set, same values.
        assert set(res_lib.vars.index) == set(res_sub.vars.index)
        for col in res_lib.vars.index:
            if col == "Name":
                continue
            v_lib = res_lib.vars[col]
            v_sub = res_sub.vars[col]
            try:
                v_lib_f = float(v_lib)
                v_sub_f = float(v_sub)
            except (ValueError, TypeError):
                # Non-numeric -- skip
                continue
            assert math.isclose(v_lib_f, v_sub_f, rel_tol=1e-10), (
                f"{col}: library={v_lib_f}, subprocess={v_sub_f}"
            )
