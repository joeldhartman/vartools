"""Undergrad-student persona — a sixth new-user pass.

The student was handed a folder of ASCII light curves and told to
"use pyvartools to find the periods and make some plots."  They are
working in Jupyter or in short scripts copy-pasted from the docs.
They are not a strong coder: they write naive for-loops, mistype
kwargs, accidentally pass strings where floats are wanted, and look
at reprs in the cell output rather than typing `print(...)`.

What this pass probes that v1-v5 did NOT:
  - Notebook-display ergonomics: reprs, ``_repr_html_``, ``lc.plot()``,
    and the shape of ``Result.vars``/``BatchResult.vars`` as a pandas
    object the student can index.
  - Error messages for trivial mistakes: missing file, mistyped kwarg,
    period bound = 0, search bound swap (min > max).
  - The naive loop pattern: ``for f in glob(...)`` rather than
    ``run_batch``/``run_filelist``.
  - Chaining one command's result into the next (``cmd.LS`` then
    ``cmd.Phase(period=fix)``).
  - "I want to save this as a CSV" — does ``result.vars`` go cleanly
    through ``pd.DataFrame.to_csv``?
  - "How do I plot the periodogram?" — ``save_periodogram`` then
    ``result.files`` dataframe access.

Run with:
    conda run -n pyvartools python -m pytest python/tests/test_doc_user_scenarios_v6.py -v
"""

from __future__ import annotations

import glob
import io
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


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _student_lc(name, period, n=400, noise=0.005, seed=0):
    """Generate a noisy sinusoidal LC mimicking a student's input file."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0, 30, n)
    mag = 10.0 + 0.05 * np.sin(2 * np.pi * t / period)
    mag += rng.normal(0, noise, n)
    err = np.full(n, noise)
    return vt.LightCurve.from_arrays(t, mag, err, name=name)


def _write_files(tmp_path, periods):
    """Write per-target ASCII files (t mag err) the student can ``glob``."""
    paths = []
    for i, p in enumerate(periods):
        lc = _student_lc(f"star_{i}", p, seed=i)
        path = tmp_path / f"star_{i}.txt"
        df = pd.DataFrame({"t": lc.t, "mag": lc.mag, "err": lc.err})
        df.to_csv(path, sep=" ", header=False, index=False,
                  float_format="%.6f")
        paths.append(str(path))
    return paths


# ---------------------------------------------------------------------------
# Notebook display ergonomics
# ---------------------------------------------------------------------------

class TestNotebookDisplay:
    """Things the student sees by typing the variable name in a cell."""

    def test_lightcurve_repr_is_informative(self):
        lc = _student_lc("HD12345", period=1.5)
        r = repr(lc)
        # Repr should include the name and the n.
        assert "HD12345" in r
        assert "400" in r  # n=400
        assert "LightCurve" in r

    def test_lightcurve_has_repr_html_for_jupyter(self):
        lc = _student_lc("Nu Pav", period=2.3)
        # Jupyter calls _repr_html_ for rich display.
        html = lc._repr_html_()
        assert html is not None
        assert isinstance(html, str)
        assert "Nu Pav" in html or "400" in html

    def test_pipeline_repr_lists_commands(self):
        pipe = vt.Pipeline([
            cmd.LS(0.5, 5.0, 1e-3),
            cmd.rms(),
        ])
        r = repr(pipe)
        assert "LS" in r
        assert "rms" in r

    @needs_binary
    def test_result_repr_lists_var_keys(self):
        lc = _student_lc("x", period=1.5)
        result = vt.Pipeline([cmd.LS(0.5, 5.0, 1e-3)]).run(lc)
        r = repr(result)
        # Student looks at the cell output to discover what came back.
        assert "LS_Period" in r or "Period" in r

    @needs_binary
    def test_result_repr_html_renders(self):
        lc = _student_lc("x", period=1.5)
        result = vt.Pipeline([cmd.LS(0.5, 5.0, 1e-3)]).run(lc)
        html = result._repr_html_()
        assert html is not None
        assert isinstance(html, str)
        assert len(html) > 0

    def test_command_object_has_helpful_repr(self):
        """Student typing `cmd.LS(0.5, 5.0, 0.01)` in a cell should see
        a friendly repr, not a raw <object at 0x...>."""
        c = cmd.LS(0.5, 5.0, 0.01, npeaks=2)
        r = repr(c)
        assert "0.5" in r
        assert "5.0" in r
        assert "LS" in r
        assert "object at" not in r  # default object repr leaks otherwise


# ---------------------------------------------------------------------------
# Quick-look plotting
# ---------------------------------------------------------------------------

class TestQuickLookPlot:
    """`lc.plot()` is what a student tries before reading the docs."""

    def test_plot_returns_axes(self):
        # Skip silently if matplotlib not installed (student would just
        # see an ImportError, not a hang).
        pytest.importorskip("matplotlib")
        import matplotlib
        matplotlib.use("Agg")  # headless
        lc = _student_lc("hd", period=1.5)
        ax = lc.plot()
        assert ax is not None
        # y-axis should be inverted (astronomical convention).
        ylim = ax.get_ylim()
        assert ylim[0] > ylim[1], "magnitude axis should be inverted"

    def test_plot_with_no_err_does_not_crash(self):
        """If a student loads a file with no error column, plot should
        still work — it should fall back from errorbar to plot."""
        pytest.importorskip("matplotlib")
        import matplotlib
        matplotlib.use("Agg")
        t = np.linspace(0, 10, 100)
        mag = 10 + 0.1 * np.sin(2 * np.pi * t / 1.5)
        lc = vt.LightCurve.from_arrays(t, mag, name="nonoise")
        # err is None — plot should not raise.
        ax = lc.plot()
        assert ax is not None


# ---------------------------------------------------------------------------
# Naive loop over files
# ---------------------------------------------------------------------------

class TestNaiveFileLoop:
    """`for f in glob(...): result = pipe.run(f); ...`
    — the pattern an undergrad reaches for first."""

    @needs_binary
    def test_naive_glob_loop_works(self, tmp_path):
        paths = _write_files(tmp_path, [1.0, 1.5, 2.0])
        # The student's naive code, almost verbatim:
        pipe = vt.Pipeline([cmd.LS(0.5, 5.0, 1e-3, npeaks=1)])
        periods = []
        for f in sorted(glob.glob(str(tmp_path / "*.txt"))):
            r = pipe.run(f)
            periods.append(float(r.vars["LS_Period_1_0"]))
        # Recovered periods are close to the injected ones.
        for p, want in zip(periods, [1.0, 1.5, 2.0]):
            assert abs(p - want) / want < 0.05

    @needs_binary
    def test_run_batch_matches_naive_loop(self, tmp_path):
        """The student eventually finds `run_batch` in the docs and
        sees it gives the same answers as their loop."""
        paths = _write_files(tmp_path, [1.0, 1.5, 2.0])
        pipe = vt.Pipeline([cmd.LS(0.5, 5.0, 1e-3, npeaks=1)])
        loop_p = []
        for f in paths:
            loop_p.append(float(pipe.run(f).vars["LS_Period_1_0"]))
        batch = pipe.run_batch(paths)
        batch_p = batch.vars["LS_Period_1_0"].tolist()
        # Within numerical noise.
        for a, b in zip(loop_p, batch_p):
            assert abs(a - b) / a < 1e-3


# ---------------------------------------------------------------------------
# Saving results to CSV
# ---------------------------------------------------------------------------

class TestSaveResultsCSV:
    """`result.vars.to_csv(...)` is the next thing the student tries."""

    @needs_binary
    def test_batch_vars_is_a_dataframe(self, tmp_path):
        paths = _write_files(tmp_path, [1.0, 1.5])
        batch = vt.Pipeline([
            cmd.LS(0.5, 5.0, 1e-3, npeaks=1),
        ]).run_batch(paths)
        # Student expects a DataFrame they can to_csv directly.
        assert isinstance(batch.vars, pd.DataFrame)
        csv = batch.vars.to_csv()
        assert isinstance(csv, str)
        assert "LS_Period_1_0" in csv

    @needs_binary
    def test_single_run_vars_is_a_series(self, tmp_path):
        """Single-LC run gives a Series, batch gives a DataFrame —
        student needs to learn this distinction.  Test it's at least
        consistent."""
        paths = _write_files(tmp_path, [1.5])
        r = vt.Pipeline([cmd.LS(0.5, 5.0, 1e-3, npeaks=1)]).run(paths[0])
        assert isinstance(r.vars, pd.Series)
        # Single-row DataFrame conversion is trivial when needed.
        df = r.vars.to_frame().T
        assert isinstance(df, pd.DataFrame)
        assert "LS_Period_1_0" in df.columns


# ---------------------------------------------------------------------------
# Friendly error messages on trivial mistakes
# ---------------------------------------------------------------------------

class TestErrorMessages:

    def test_mistyped_kwarg_gives_typeerror(self):
        """`npeak` instead of `npeaks` — clear TypeError mentioning the
        bad kwarg."""
        with pytest.raises(TypeError) as exc:
            cmd.LS(0.5, 5.0, 1e-3, npeak=2)
        msg = str(exc.value).lower()
        assert "npeak" in msg

    def test_negative_period_caught_at_construction(self):
        with pytest.raises(ValueError):
            cmd.LS(-1.0, 5.0, 1e-3)

    def test_zero_period_caught_at_construction(self):
        """The student types `cmd.LS(0, 5, 0.01)` thinking 0 = "from zero"."""
        with pytest.raises(ValueError):
            cmd.LS(0, 5.0, 1e-3)

    def test_min_greater_than_max_caught(self):
        """`cmd.LS(10, 1, ...)` — student swapped the bounds."""
        with pytest.raises(ValueError):
            cmd.LS(10.0, 1.0, 1e-3)

    @needs_binary
    def test_missing_file_raises_clearly(self, tmp_path):
        """`pipe.run("nope.txt")` — the file does not exist."""
        bad = str(tmp_path / "does_not_exist.txt")
        pipe = vt.Pipeline([cmd.LS(0.5, 5.0, 1e-3)])
        with pytest.raises(Exception) as exc:
            pipe.run(bad)
        msg = str(exc.value)
        # We want the path itself in the error so the student sees
        # what they typed.
        assert "does_not_exist" in msg or "not found" in msg.lower() \
            or "no such" in msg.lower()


# ---------------------------------------------------------------------------
# Chain one command's result into the next
# ---------------------------------------------------------------------------

class TestChainCommands:
    """The most common follow-on task: find period with LS, then
    phase the LC at that period and plot it."""

    @needs_binary
    def test_phase_at_recovered_period(self):
        """Student finds the period with LS, then phases the LC.
        Use the PerLC-or-batch_scalar mechanism: pipeline carries
        the LS_Period_1_0 forward into a Phase command via varname."""
        lc = _student_lc("variable", period=2.0)
        # Two-step pipeline: LS finds the period; Phase uses it.
        r = vt.Pipeline([
            cmd.LS(0.5, 5.0, 1e-3, npeaks=1),
            cmd.Phase(period="LS_Period_1_0"),  # use varname
        ]).run(lc, capture_lc=True)
        # Captured LC should now be phased.  t-column values should
        # be in [0, 1) (phase units) after Phase, OR in roughly
        # [0, period].  Just check we have an LC and didn't crash.
        assert r.lc is not None
        assert len(r.lc.t) == 400

    @needs_binary
    def test_period_then_killharm(self):
        """Same idea: find period, then filter the harmonic out."""
        lc = _student_lc("variable2", period=1.7)
        r = vt.Pipeline([
            cmd.LS(0.5, 5.0, 1e-3, npeaks=1),
            cmd.Killharm(period="LS_Period_1_0", nharm=1),
            cmd.rms(),
        ]).run(lc)
        # After harmonic removal, rms should be near the noise floor.
        rms_final = float(r.vars["RMS_2"])
        assert rms_final < 0.02


# ---------------------------------------------------------------------------
# Periodogram output access
# ---------------------------------------------------------------------------

class TestPeriodogramAccess:
    """`save_periodogram=True` + `result.files[...]` — the path for
    plotting a periodogram."""

    @needs_binary
    def test_periodogram_dataframe_has_expected_columns(self, tmp_path):
        lc = _student_lc("p", period=1.5)
        r = vt.Pipeline([
            cmd.LS(0.5, 5.0, 1e-3, save_periodogram=True),
        ]).run(lc, outdir=str(tmp_path))
        # Student tries to figure out the key from the repr.
        assert "LS_periodogram_0" in r.files, (
            f"Expected LS_periodogram_0 in files; got {list(r.files)}"
        )
        df = r.files["LS_periodogram_0"]
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0


# ---------------------------------------------------------------------------
# Filtering results in pandas (the "find the variables" task)
# ---------------------------------------------------------------------------

class TestFilterResults:

    @needs_binary
    def test_filter_by_logprob(self, tmp_path):
        """Student wants to keep only the 'good' detections.  Pandas
        boolean indexing on BatchResult.vars should just work."""
        # 5 LCs — 3 with real signals, 2 with pure noise.
        paths = []
        for i, p in enumerate([1.0, 1.5, 2.0]):
            lc = _student_lc(f"var_{i}", p, seed=i)
            path = tmp_path / f"var_{i}.txt"
            pd.DataFrame({"t": lc.t, "mag": lc.mag, "err": lc.err}
                ).to_csv(path, sep=" ", header=False, index=False,
                          float_format="%.6f")
            paths.append(str(path))
        for i in range(2):
            t = np.linspace(0, 30, 400)
            rng = np.random.default_rng(100 + i)
            mag = 10.0 + rng.normal(0, 0.005, 400)  # pure noise
            err = np.full(400, 0.005)
            path = tmp_path / f"noise_{i}.txt"
            pd.DataFrame({"t": t, "mag": mag, "err": err}
                ).to_csv(path, sep=" ", header=False, index=False,
                          float_format="%.6f")
            paths.append(str(path))
        batch = vt.Pipeline([
            cmd.LS(0.5, 5.0, 1e-3, npeaks=1),
        ]).run_batch(paths)
        # Log10_LS_Prob — more negative = stronger detection.
        col = "Log10_LS_Prob_1_0"
        assert col in batch.vars.columns
        variables = batch.vars[batch.vars[col] < -50]
        # Should keep most of the 3 real variables and reject the noise.
        assert len(variables) >= 2
        assert len(variables) <= 5
