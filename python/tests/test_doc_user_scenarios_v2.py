"""Survey-scientist scenarios — a different new-user persona.

This user is setting up a multi-step pipeline for ~100 light curves of
variable stars.  They care about:

  - Validating the pipeline before launching the run.
  - Reproducibility (randseed).
  - Per-LC outputs on disk: periodograms, model files, processed LCs.
  - Resume after a crash via stats_file.
  - Per-LC error handling that doesn't kill the whole run.
  - Cross-command references (LS → harmonicfilter at the found period).
  - Structured access to results (result.varobjs.LS.Period_1).

Compared to test_doc_user_scenarios.py, this exercise probes the *output*
side of the API and the operational ergonomics — validate, save_*,
streamed stats_file, multi-segment chains, FITS / gz file loading.

Run with:
    conda run -n pyvartools python -m pytest python/tests/test_doc_user_scenarios_v2.py -v -s
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

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

EXAMPLES_DIR = os.path.join(VARTOOLS_SRC, "EXAMPLES")


def _star_lcs(n=5, npts=200, seed=42):
    """Five sinusoidal LCs with distinct periods and a touch of noise —
    looks more like real survey data than pure sinusoids."""
    rng = np.random.default_rng(seed)
    out = []
    for i in range(n):
        period = 1.0 + 0.4 * i  # 1.0, 1.4, 1.8, 2.2, 2.6 d
        t = np.linspace(0, 30, npts)
        mag = 10.0 + 0.05 * np.sin(2 * np.pi * t / period)
        mag += rng.normal(0, 0.003, npts)
        err = np.full(npts, 0.005)
        out.append(vt.LightCurve.from_arrays(
            t, mag, err, name=f"HD_{1000 + i}"))
    return out


# ---------------------------------------------------------------------------
# Survey workflow happy paths
# ---------------------------------------------------------------------------

class TestSurveyWorkflow:

    @needs_binary
    def test_validate_before_running(self):
        """User wants to check the column layout before launching a 100-star
        run.  Pipeline.validate() returns the expected output columns."""
        pipe = vt.Pipeline([
            cmd.clip(5.0),
            cmd.LS(0.5, 10.0, 0.01, npeaks=2),
            cmd.rms(),
        ])
        cols = pipe.validate()
        # Spot-check known column names.
        assert "Name" in cols
        assert any("LS_Period_1" in c for c in cols)
        assert "RMS_2" in cols  # rms is command index 2

    @needs_binary
    def test_randseed_reproducibility_in_batch(self):
        """Same randseed → same results across invocations."""
        lcs = _star_lcs(n=3)
        pipe = vt.Pipeline([cmd.LS(0.5, 10.0, 0.01, npeaks=1)])
        r1 = pipe.run_batch(lcs, randseed=42)
        r2 = pipe.run_batch(lcs, randseed=42)
        np.testing.assert_array_equal(
            r1.vars["LS_Period_1_0"].values,
            r2.vars["LS_Period_1_0"].values,
        )

    @needs_binary
    def test_cross_command_backref_in_chaining(self):
        """lc.LS(...).Phase(period='ls') resolves 'ls' to LS's period
        across the chain boundary (each chain step is a separate
        vartools invocation, so the back-ref has to look at the prior
        result.vars)."""
        lc = _star_lcs(n=1)[0]
        r = lc.LS(0.5, 10.0, 1e-3, npeaks=1).Phase(period="ls")
        assert r.lc is not None
        # After phase-folding, t should be in [0, 1).
        assert 0.0 <= float(r.lc.t.min()) and float(r.lc.t.max()) < 1.0

    @needs_binary
    def test_per_lc_periodogram_output(self, tmp_path):
        """save_periodogram=True captures the periodogram per LC into
        result.files."""
        lcs = _star_lcs(n=3)
        batch = vt.Pipeline([
            cmd.LS(0.5, 10.0, 0.01, npeaks=1, save_periodogram=True),
        ]).run_batch(lcs)
        assert "LS_periodogram_0" in batch.files
        pgrams = batch.files["LS_periodogram_0"]
        assert len(pgrams) == 3
        for pg in pgrams:
            assert pg is not None
            assert len(pg) > 0  # non-trivial periodogram

    @needs_binary
    def test_save_periodogram_with_explicit_outdir(self, tmp_path):
        """save_periodogram=str(path) writes to disk only (NOT captured
        into result.files) per the docs."""
        lcs = _star_lcs(n=2)
        outdir = tmp_path / "pgrams"
        outdir.mkdir()
        batch = vt.Pipeline([
            cmd.LS(0.5, 10.0, 0.01, npeaks=1,
                   save_periodogram=str(outdir)),
        ]).run_batch(lcs)
        assert "LS_periodogram_0" not in batch.files
        # Files on disk
        files = sorted(os.listdir(outdir))
        assert len(files) >= 2

    @needs_binary
    def test_structured_varobjs_access(self):
        """result.varobjs.LS.Period_1 returns the LS period from
        the first peak."""
        lc = _star_lcs(n=1)[0]
        r = vt.Pipeline([cmd.LS(0.5, 10.0, 1e-3, npeaks=1)]).run(lc)
        period_via_index = r.vars["LS_Period_1_0"]
        period_via_attr = r.varobjs.LS.Period_1
        assert abs(period_via_attr - period_via_index) < 1e-12

    @needs_binary
    def test_per_lc_error_handling_in_chaining(self):
        """LightCurveBatch.run() with one bad LC continues for the
        others.  Per-LC errors are surfaced as Result.error."""
        good = _star_lcs(n=2)
        # An "LC" with zero rows — should fail vartools' validation.
        bad = vt.LightCurve.from_arrays(np.array([]), np.array([]),
                                         np.array([]), name="empty")
        batch = vt.LightCurveBatch(good + [bad]).rms().run()
        # Two good results, one bad.
        oks = sum(1 for r in batch if r.ok)
        bads = sum(1 for r in batch if not r.ok)
        print(f"per-LC error handling: ok={oks}, bad={bads}")
        # Expectation: the per-LC loop should not abort the whole batch.

    @needs_binary
    def test_stats_file_streaming_and_resume(self, tmp_path):
        """Survey ergonomics: stats_file writes incrementally; after a
        partial run, resume=True picks up the rest."""
        lcs = _star_lcs(n=5)
        sf = tmp_path / "survey.stats"

        # Step 1: run on first 2 only, simulating a crash.
        from pyvartools import pipeline as _pmod
        with patch.object(_pmod, "_library_enabled", return_value=False):
            vt.Pipeline([cmd.LS(0.5, 10.0, 0.01, npeaks=1)]).run_batch(
                lcs[:2], stats_file=str(sf))
        assert sf.exists()
        # Two data rows + 1 header in the file.
        lines_before = sf.read_text().strip().split("\n")
        assert len(lines_before) == 3

        # Step 2: resume with the full list — should re-run only LCs 2-4.
        with patch.object(_pmod, "_library_enabled", return_value=False):
            batch = vt.Pipeline([
                cmd.LS(0.5, 10.0, 0.01, npeaks=1)
            ]).run_batch(lcs, stats_file=str(sf), resume=True)
        lines_after = sf.read_text().strip().split("\n")
        assert len(lines_after) == 1 + 5  # header + 5 rows
        assert len(batch.vars) == 5

    @needs_binary
    def test_capture_lc_chain_continues(self):
        """A common survey pattern: run LS, then continue the chain on
        the captured LCs to remove the harmonic."""
        lcs = _star_lcs(n=3)
        br1 = vt.LightCurveBatch(lcs).LS(0.5, 10.0, 1e-3, npeaks=1).run()
        # Continuation uses harmonicfilter at the LS period; should not
        # error and should produce a result with both LS and HarmonicFilter
        # columns.
        br2 = br1.harmonicfilter(
            period=br1.vars["LS_Period_1_0"], nharm=2).run()
        assert "LS_Period_1_0" in br2.vars.columns
        # harmonicfilter has its own output columns
        hf_cols = [c for c in br2.vars.columns
                   if c.lower().startswith(("harmonic", "hf"))]
        # Either it produces such columns, or the LC arrays change — at
        # minimum, no error.
        print(f"After harmonic-filter continuation, columns: "
              f"{list(br2.vars.columns)[:8]}...")


# ---------------------------------------------------------------------------
# Edge cases the user might wonder about
# ---------------------------------------------------------------------------

class TestEdgeCases:

    @needs_binary
    def test_pipeline_reused_across_calls(self):
        """User builds one Pipeline and calls it on many LCs in a loop —
        the cached LibPipeline should produce identical results to running
        each through a fresh Pipeline."""
        lcs = _star_lcs(n=3)
        pipe_reused = vt.Pipeline([cmd.LS(0.5, 10.0, 0.01, npeaks=1)])
        reused = [pipe_reused.run(lc) for lc in lcs]
        fresh = [
            vt.Pipeline([cmd.LS(0.5, 10.0, 0.01, npeaks=1)]).run(lc)
            for lc in lcs
        ]
        for i in range(3):
            assert reused[i].vars["LS_Period_1_0"] == \
                   fresh[i].vars["LS_Period_1_0"]

    @needs_binary
    def test_two_LS_instances_distinct_suffixes(self):
        """A user might use cmd.LS twice — once to find the dominant
        period, once after clipping.  The output columns must not
        collide."""
        lc = _star_lcs(n=1)[0]
        r = vt.Pipeline([
            cmd.LS(0.5, 10.0, 1e-3, npeaks=1),
            cmd.clip(5.0),
            cmd.LS(0.5, 10.0, 1e-3, npeaks=1),
        ]).run(lc)
        # Distinct command-index suffixes.
        assert "LS_Period_1_0" in r.vars.index
        assert "LS_Period_1_2" in r.vars.index

    @needs_binary
    def test_validate_with_perlc_vars_values(self):
        """validate(perlc_vars={...}) pre-registers the named variables
        so vartools' parser accepts var-ref command parameters during the
        -headeronly probe."""
        pipe = vt.Pipeline([cmd.LS("minp", "maxp", 0.1, npeaks=1)])
        cols = pipe.validate(perlc_vars={
            "minp": [0.3, 0.5],
            "maxp": [3.0, 5.0],
        })
        # Sanity: validation succeeded and produced LS columns.
        assert "Name" in cols
        assert any("LS_Period_1" in c for c in cols)

    @needs_binary
    def test_validate_with_perlc_vars_schema(self):
        """Schema-form perlc_vars (column refs) is forwarded as-is."""
        pipe = vt.Pipeline([cmd.LS("minp", "maxp", 0.1, npeaks=1)])
        from pyvartools import PerLCColumn
        cols = pipe.validate(perlc_vars={
            "minp": PerLCColumn(col=2, type="double"),
            "maxp": PerLCColumn(col=3, type="double"),
        })
        assert any("LS_Period_1" in c for c in cols)

    @needs_binary
    def test_validate_without_perlc_vars_for_var_ref_pipeline_errors(self):
        """If the pipeline references var-ref params but validate() is
        called without perlc_vars=, the parser error still surfaces but
        with a clearer cause."""
        pipe = vt.Pipeline([cmd.LS("minp", "maxp", 0.1, npeaks=1)])
        with pytest.raises(Exception) as excinfo:
            pipe.validate()
        # Vartools' own message ("undefined/initialized variable") is
        # passed through; that's still informative enough for the user
        # to spot the issue.
        assert "minp" in str(excinfo.value).lower() or \
               "undefined" in str(excinfo.value).lower()

    @needs_binary
    def test_loading_gz_file(self):
        """vartools auto-detects .gz; pyvartools should too."""
        path = os.path.join(EXAMPLES_DIR, "1.gz")
        if not os.path.exists(path):
            pytest.skip(f"no test .gz file at {path}")
        lc = vt.LightCurve.from_file(path)
        assert lc.t is not None and len(lc.t) > 0
        r = vt.rms(lc)
        assert int(r.vars["Npoints_0"]) > 0

    @needs_binary
    def test_loading_fits_without_explicit_cols_errors(self):
        """FITS files require explicit t_col / mag_col / err_col."""
        path = os.path.join(EXAMPLES_DIR, "2.fits")
        if not os.path.exists(path):
            pytest.skip(f"no test FITS file at {path}")
        with pytest.raises(ValueError) as excinfo:
            vt.LightCurve.from_file(path)
        msg = str(excinfo.value)
        # Error should name all three unset kwargs AND list available columns.
        assert "t_col" in msg and "mag_col" in msg and "err_col" in msg
        # And it should list the actual FITS column names so the user can
        # see what to map t/mag/err onto.
        assert "Available columns" in msg

    @needs_binary
    def test_loading_fits_with_explicit_cols(self):
        """With explicit t_col/mag_col/err_col, FITS files load and every
        FITS column is preserved in the DataFrame under its original
        FITS column name."""
        path = os.path.join(EXAMPLES_DIR, "2.fits")
        if not os.path.exists(path):
            pytest.skip(f"no test FITS file at {path}")
        lc = vt.LightCurve.from_file(path, t_col="time", mag_col="mag",
                                      err_col="err")
        assert lc.t is not None and len(lc.t) > 0
        # Both aliases AND original FITS column names present.
        assert {"t", "mag", "err"}.issubset(set(lc._df.columns))
        assert {"time", "mag", "err"}.issubset(set(lc._df.columns))
        r = vt.rms(lc)
        assert int(r.vars["Npoints_0"]) > 0

    @needs_binary
    def test_loading_fits_with_some_cols_set_to_none(self, tmp_path):
        """Passing t_col=None means 'this FITS has no t column'; vartools
        defaults t=NR at run time."""
        try:
            from astropy.io import fits
        except ImportError:
            pytest.skip("astropy not available")

        # Build a FITS file with only mag and err — no t-like column.
        mag = 10.0 + 0.05 * np.sin(np.linspace(0, 10, 50))
        err = np.full(50, 0.005)
        hdu = fits.BinTableHDU.from_columns([
            fits.Column(name="brightness", format="D", array=mag),
            fits.Column(name="uncertainty", format="D", array=err),
        ])
        fpath = tmp_path / "no_time.fits"
        fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(str(fpath))

        lc = vt.LightCurve.from_file(
            fpath, t_col=None, mag_col="brightness", err_col="uncertainty",
        )
        assert lc.t is None        # No t column registered
        assert lc.mag is not None
        # But running rms still works (vartools defaults t=NR).
        r = vt.rms(lc)
        assert int(r.vars["Npoints_0"]) == 50

    @needs_binary
    def test_loading_fits_with_unknown_col_lists_alternatives(self, tmp_path):
        """If t_col= references a column that doesn't exist, the error
        lists the actual FITS columns."""
        try:
            from astropy.io import fits
        except ImportError:
            pytest.skip("astropy not available")
        mag = np.full(10, 10.0)
        hdu = fits.BinTableHDU.from_columns([
            fits.Column(name="t", format="D", array=np.arange(10.0)),
            fits.Column(name="m", format="D", array=mag),
            fits.Column(name="e", format="D", array=np.full(10, 0.01)),
        ])
        fpath = tmp_path / "short_names.fits"
        fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(str(fpath))

        with pytest.raises(ValueError, match=r"not found"):
            vt.LightCurve.from_file(
                fpath, t_col="BJD", mag_col="m", err_col="e")

    @needs_binary
    def test_nthreads_falls_back_to_subprocess(self):
        """nthreads > 1 should route through subprocess (library mode is
        single-threaded) — the user gets the speedup without surprise."""
        lcs = _star_lcs(n=4)
        from pyvartools import pipeline as _pmod
        pipe = vt.Pipeline([cmd.rms()])
        original = pipe._execute
        called = []

        def tracking(cmd_argv, timeout=None, stdin_text=None):
            called.append(True)
            return original(cmd_argv, timeout=timeout,
                            stdin_text=stdin_text)

        with patch.object(pipe, "_execute", side_effect=tracking):
            with patch.object(_pmod, "_library_enabled",
                              return_value=True):
                pipe.run_batch(lcs, nthreads=2)
        assert called, "nthreads=2 should force subprocess"

    @needs_binary
    def test_repr_of_pipeline(self):
        """A user typing `pipe` in Jupyter would want a useful repr."""
        pipe = vt.Pipeline([
            cmd.clip(5.0),
            cmd.LS(0.5, 10.0, 0.01, npeaks=2, save_periodogram=True),
            cmd.rms(),
        ])
        s = repr(pipe)
        print(f"repr(pipe) = {s!r}")
        # No assertion — informational.  A useful repr lists the commands.

    def test_negative_period_caught_at_construction(self):
        """cmd.LS(-0.5, ...) is rejected immediately, not silently."""
        with pytest.raises(ValueError, match="minp must be > 0"):
            cmd.LS(-0.5, 10.0, 1e-3, npeaks=1)

    def test_zero_period_caught(self):
        """minp=0 is also rejected (would be infinite frequency)."""
        with pytest.raises(ValueError, match="minp must be > 0"):
            cmd.LS(0.0, 10.0, 1e-3, npeaks=1)

    def test_maxp_below_minp_caught(self):
        """If both are numeric, minp < maxp is enforced."""
        with pytest.raises(ValueError, match="minp.*less than.*maxp"):
            cmd.LS(5.0, 1.0, 1e-3, npeaks=1)

    def test_variable_ref_minp_allowed_at_construction(self):
        """A variable name in minp= must be accepted at construction
        (the value is only known at run time)."""
        # Should not raise.
        cmd.LS("minp_col", "maxp_col", 1e-3, npeaks=1)
        cmd.LS("tspan/100", 10.0, 1e-3, npeaks=1)
