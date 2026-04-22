"""
Feature A tests: stdin pipe for single-LC runs.

These tests verify that Pipeline.run() passes the light curve via stdin
(-i -) rather than writing a temp input.lc file, while producing results
identical to the old file-based path.
"""

from __future__ import annotations

import os
import sys
import math
import subprocess
from pathlib import Path
from unittest.mock import patch, MagicMock

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Path setup (allows running without installing the package)
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

EXAMPLES_DIR = os.path.join(VARTOOLS_SRC, "EXAMPLES")
EXAMPLE_LC = os.path.join(EXAMPLES_DIR, "2")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_lc():
    """A short, noise-free LightCurve for quick tests."""
    t = np.linspace(0, 10, 100)
    mag = np.full(100, 10.0)
    err = np.full(100, 0.01)
    return vt.LightCurve.from_arrays(t, mag, err, name="testlc")


@pytest.fixture
def sinusoidal_lc():
    """LightCurve with a sinusoidal signal of known period 1.5 d."""
    rng = np.random.default_rng(42)
    t = np.linspace(0, 30, 200)
    mag = 10.0 + 0.1 * np.sin(2 * np.pi * t / 1.5)
    err = np.full(200, 0.005)
    return vt.LightCurve.from_arrays(t, mag, err, name="sinlc")


# ---------------------------------------------------------------------------
# Structural tests (no binary required)
# ---------------------------------------------------------------------------

class TestStdinStructure:
    """Verify that run() uses -i - rather than a temp file."""

    def test_stdin_mode_uses_dash_i_dash(self, simple_lc):
        """_execute should be called with ['-i', '-'] not ['-i', '<path>']."""
        from pyvartools import pipeline as _pipeline_mod
        pipe = vt.Pipeline([cmd.rms()])
        captured_cmds = []

        original_execute = pipe._execute

        def mock_execute(command, timeout=None, stdin_text=None):
            captured_cmds.append(command)
            # Return a minimal valid -oneline response
            return "Name = testlc\nRMS_0 = 0.01\n", ""

        with patch.object(_pipeline_mod, "_library_enabled", return_value=False):
            with patch.object(pipe, "_execute", side_effect=mock_execute):
                pipe.run(simple_lc)

        assert len(captured_cmds) == 1
        cmd_list = captured_cmds[0]
        assert "-i" in cmd_list
        i_idx = cmd_list.index("-i")
        assert cmd_list[i_idx + 1] == "-", (
            f"Expected '-i -' for stdin, got '-i {cmd_list[i_idx + 1]}'"
        )

    def test_stdin_mode_passes_text_to_execute(self, simple_lc):
        """_execute should receive the LC data as stdin_text (non-None string)."""
        from pyvartools import pipeline as _pipeline_mod
        pipe = vt.Pipeline([cmd.rms()])
        received_stdin = []

        def mock_execute(command, timeout=None, stdin_text=None):
            received_stdin.append(stdin_text)
            return "Name = testlc\nRMS_0 = 0.01\n", ""

        with patch.object(_pipeline_mod, "_library_enabled", return_value=False):
            with patch.object(pipe, "_execute", side_effect=mock_execute):
                pipe.run(simple_lc)

        assert len(received_stdin) == 1
        assert received_stdin[0] is not None
        assert isinstance(received_stdin[0], str)
        # Should contain the numeric data
        assert "10.0" in received_stdin[0] or "10." in received_stdin[0]

    def test_no_input_lc_file_created(self, simple_lc, tmp_path):
        """run() must not write an input.lc temp file."""
        import tempfile

        created_files = []
        original_TemporaryDirectory = tempfile.TemporaryDirectory

        class CapturingTempDir:
            def __init__(self):
                self._td = original_TemporaryDirectory()
                self.name = self._td.name

            def __enter__(self):
                return self.name

            def __exit__(self, *args):
                # Capture files before cleanup
                try:
                    for fname in os.listdir(self.name):
                        created_files.append(fname)
                except Exception:
                    pass
                self._td.cleanup()

        pipe = vt.Pipeline([cmd.rms()])

        def mock_execute(command, timeout=None, stdin_text=None):
            return "Name = testlc\nRMS_0 = 0.01\n", ""

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch("pyvartools.pipeline.tempfile.TemporaryDirectory",
                       CapturingTempDir):
                pipe.run(simple_lc)

        assert "input.lc" not in created_files, (
            f"input.lc was written to tmpdir (files found: {created_files})"
        )

    def test_run_file_uses_path_not_stdin(self):
        """run_file() must pass the file path directly, not -i -."""
        pipe = vt.Pipeline([cmd.rms()])
        captured_cmds = []

        def mock_execute(command, timeout=None, stdin_text=None):
            captured_cmds.append((command, stdin_text))
            return "Name = testlc\nRMS_0 = 0.01\n", ""

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            pipe.run_file("/some/path/lc.dat")

        assert len(captured_cmds) == 1
        cmd_list, stdin_text = captured_cmds[0]
        assert "-i" in cmd_list
        i_idx = cmd_list.index("-i")
        assert cmd_list[i_idx + 1] != "-", (
            "run_file() should pass the file path, not stdin"
        )
        assert stdin_text is None

    def test_run_batch_does_not_use_stdin(self, simple_lc):
        """run_batch() uses -l listfile; stdin should not be used."""
        from pyvartools import pipeline as _pipeline_mod
        pipe = vt.Pipeline([cmd.rms()])
        captured_cmds = []

        def mock_execute(command, timeout=None, stdin_text=None):
            captured_cmds.append((command, stdin_text))
            return "Name = lc1\nRMS_0 = 0.01\nName = lc2\nRMS_0 = 0.01\n", ""

        with patch.object(_pipeline_mod, "_library_enabled", return_value=False):
            with patch.object(pipe, "_execute", side_effect=mock_execute):
                pipe.run_batch([simple_lc, simple_lc])

        assert len(captured_cmds) == 1
        cmd_list, stdin_text = captured_cmds[0]
        assert "-l" in cmd_list, "run_batch() should use -l listfile"
        assert stdin_text is None


# ---------------------------------------------------------------------------
# Functional tests (require vartools binary)
# ---------------------------------------------------------------------------

class TestStdinFunctional:

    @needs_binary
    def test_stdin_with_rms(self, simple_lc):
        """Smoke test: run with stdin, check RMS_0 is a finite float."""
        pipe = vt.Pipeline([cmd.rms()])
        result = pipe.run(simple_lc)
        assert "RMS_0" in result.vars.index
        assert math.isfinite(float(result.vars["RMS_0"]))

    @needs_binary
    def test_stdin_produces_same_stats_as_file_mode(self):
        """stdin and file-based runs on the same LC must yield identical var."""
        lc = vt.LightCurve.from_file(EXAMPLE_LC)
        pipe = vt.Pipeline([cmd.rms()])

        # stdin mode (current default)
        result_stdin = pipe.run(lc)

        # file mode: monkeypatch run() to write a temp file instead
        import tempfile as _tf
        with _tf.NamedTemporaryFile(suffix=".lc", mode="w", delete=False) as f:
            lc._df.to_csv(f, sep=" ", header=False, index=False,
                          float_format="%.10f")
            tmp_path = f.name
        try:
            result_file = pipe.run_file(tmp_path)
        finally:
            os.unlink(tmp_path)

        for key in result_stdin.vars.index:
            if key == "Name":
                continue
            v_stdin = float(result_stdin.vars[key])
            v_file = float(result_file.vars[key])
            assert abs(v_stdin - v_file) < 1e-9, (
                f"{key}: stdin={v_stdin} vs file={v_file}"
            )

    @needs_binary
    def test_stdin_with_clip_and_ls(self, sinusoidal_lc):
        """Chain clip + LS via stdin; LS_Period_1_1 should be a finite float."""
        pipe = vt.Pipeline([
            cmd.clip(5.0),
            cmd.LS(0.5, 5.0, 1e-3, npeaks=1),
        ])
        result = pipe.run(sinusoidal_lc)
        assert "LS_Period_1_1" in result.vars.index
        assert math.isfinite(float(result.vars["LS_Period_1_1"]))

    @needs_binary
    def test_stdin_with_extra_columns(self):
        """LightCurve with 5 columns; -inputlcformat is still emitted."""
        df = pd.DataFrame({
            "t": np.linspace(0, 10, 50),
            "mag": np.full(50, 10.0),
            "err": np.full(50, 0.01),
            "x": np.zeros(50),
            "y": np.ones(50),
        })
        lc = vt.LightCurve.from_dataframe(df)
        pipe = vt.Pipeline([cmd.rms()])

        captured_cmds = []
        original_execute = pipe._execute

        def capturing_execute(command, timeout=None, stdin_text=None):
            captured_cmds.append(command)
            return original_execute(command, timeout=timeout, stdin_text=stdin_text)

        from pyvartools import pipeline as _pipeline_mod
        with patch.object(_pipeline_mod, "_library_enabled", return_value=False):
            with patch.object(pipe, "_execute", side_effect=capturing_execute):
                result = pipe.run(lc)

        cmd_list = captured_cmds[0]
        assert "-inputlcformat" in cmd_list, (
            "-inputlcformat should be present for non-default column layout"
        )
        assert math.isfinite(float(result.vars["RMS_0"]))

    @needs_binary
    def test_stdin_with_capture_lc(self, simple_lc):
        """capture_lc=True still returns an output LightCurve when using stdin."""
        pipe = vt.Pipeline([cmd.clip(5.0)])
        result = pipe.run(simple_lc, capture_lc=True)
        assert result.lc is not None
        assert isinstance(result.lc, vt.LightCurve)
        assert len(result.lc) == len(simple_lc)

    @needs_binary
    def test_stdin_with_cmd_o_capture(self, simple_lc):
        """cmd.o(capture=True) returns result.files["o"] as a LightCurve."""
        pipe = vt.Pipeline([cmd.o(capture=True)])
        result = pipe.run(simple_lc)
        assert "o" in result.files
        assert isinstance(result.files["o"], vt.LightCurve)
        assert len(result.files["o"]) == len(simple_lc)

    @needs_binary
    def test_stdin_with_save_periodogram(self, sinusoidal_lc):
        """save_periodogram=True: periodogram DataFrame is returned from stdin run."""
        pipe = vt.Pipeline([
            cmd.LS(0.5, 5.0, 0.01, npeaks=1, save_periodogram=True),
        ])
        result = pipe.run(sinusoidal_lc)
        key = "LS_periodogram_0"
        assert key in result.files, (
            f"Expected '{key}' in result.files; got {list(result.files.keys())}"
        )
        df = result.files[key]
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0

    @needs_binary
    def test_stdin_name_is_replaced(self, simple_lc):
        """result.vars['Name'] should be the LightCurve's name, not 'stdin'."""
        pipe = vt.Pipeline([cmd.rms()])
        result = pipe.run(simple_lc)
        assert result.vars["Name"] == "testlc"
        assert result.vars["Name"] != "stdin"
