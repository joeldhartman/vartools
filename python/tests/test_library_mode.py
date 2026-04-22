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
