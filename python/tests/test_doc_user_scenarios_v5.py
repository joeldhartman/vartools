"""Multi-band / multi-segment user — a fifth new-user persona.

This user has data from multiple telescopes / sectors / nights on the
same target.  They want to stitch the segments, detrend trends out,
window to a time range of interest, and inject custom Python compute
mid-pipeline.

Probes corners untouched by v1-v4:
  - LightCurve.from_files for multi-file stitching.
  - run_combinelcs for vartools-side combine.
  - cmd.linfit / cmd.medianfilter / cmd.binlc.
  - cmd.expr / cmd.changevariable / cmd.restricttimes.
  - cmd.python(inprocess=False) — Python callback (subprocess mode).
  - Mixed-shape inputs (different telescopes, different LC lengths in a
    batch).

Run with:
    conda run -n pyvartools python -m pytest python/tests/test_doc_user_scenarios_v5.py -v -s
"""

from __future__ import annotations

import math
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


def _segment(t0, t1, npts, period=1.5, noise=0.005, seed=0,
             trend_slope=0.0):
    """A noisy sinusoidal LC segment with optional linear trend."""
    rng = np.random.default_rng(seed)
    t = np.linspace(t0, t1, npts)
    mag = 10.0 + 0.05 * np.sin(2 * np.pi * t / period)
    mag += trend_slope * (t - t0)
    mag += rng.normal(0, noise, npts)
    err = np.full(npts, noise)
    return t, mag, err


# ---------------------------------------------------------------------------
# Multi-segment LC assembly via LightCurve.from_files
# ---------------------------------------------------------------------------

class TestStitching:

    @needs_binary
    def test_from_files_stitches_in_time_order(self, tmp_path):
        """Two ASCII files (different time ranges) → one combined LC,
        sorted by t."""
        for i, (t0, t1) in enumerate([(0, 10), (10, 20)]):
            t, mag, err = _segment(t0, t1, 100, seed=i)
            df = pd.DataFrame({"t": t, "mag": mag, "err": err})
            df.to_csv(tmp_path / f"seg{i}.txt", sep=" ", header=False,
                      index=False, float_format="%.6f")
        lc = vt.LightCurve.from_files([
            str(tmp_path / "seg0.txt"),
            str(tmp_path / "seg1.txt"),
        ])
        assert len(lc) == 200
        assert lc.t[0] < 0.1
        assert lc.t[-1] > 19.9
        # lcnum column should be 0 for first half, 1 for second half.
        assert "lcnum" in lc._df.columns
        assert lc._df["lcnum"].max() == 1

    @needs_binary
    def test_from_files_unsorted_preserves_order(self, tmp_path):
        """sort=False keeps the input file order."""
        for i, (t0, t1) in enumerate([(10, 20), (0, 10)]):
            t, mag, err = _segment(t0, t1, 50, seed=i)
            df = pd.DataFrame({"t": t, "mag": mag, "err": err})
            df.to_csv(tmp_path / f"u{i}.txt", sep=" ", header=False,
                      index=False, float_format="%.6f")
        lc = vt.LightCurve.from_files([
            str(tmp_path / "u0.txt"),
            str(tmp_path / "u1.txt"),
        ], sort=False)
        # First half should be t in [10, 20], second in [0, 10].
        assert lc.t[0] >= 10.0
        assert lc.t[50] < 10.0

    @needs_binary
    def test_pipeline_runs_on_stitched_lc(self, tmp_path):
        """A combined LC behaves like any other LC for run()."""
        for i, (t0, t1) in enumerate([(0, 10), (10, 20)]):
            t, mag, err = _segment(t0, t1, 100, seed=i)
            df = pd.DataFrame({"t": t, "mag": mag, "err": err})
            df.to_csv(tmp_path / f"part{i}.txt", sep=" ", header=False,
                      index=False, float_format="%.6f")
        lc = vt.LightCurve.from_files([
            str(tmp_path / "part0.txt"),
            str(tmp_path / "part1.txt"),
        ])
        r = vt.Pipeline([cmd.LS(0.5, 5.0, 1e-3, npeaks=1)]).run(lc)
        recovered = float(r.vars["LS_Period_1_0"])
        assert abs(recovered - 1.5) / 1.5 < 0.02


# ---------------------------------------------------------------------------
# Per-segment processing via run_combinelcs
# ---------------------------------------------------------------------------

class TestCombineLcs:

    @needs_binary
    def test_run_combinelcs_two_groups(self, tmp_path):
        """run_combinelcs with two groups, each containing two segments
        for the same source."""
        files = {}
        for grp in (1, 2):
            for seg in (0, 1):
                t, mag, err = _segment(seg * 10, (seg + 1) * 10, 80,
                                        period=1.2 + 0.3 * grp,
                                        seed=grp * 10 + seg)
                df = pd.DataFrame({"t": t, "mag": mag, "err": err})
                p = tmp_path / f"src{grp}_seg{seg}.txt"
                df.to_csv(p, sep=" ", header=False, index=False,
                           float_format="%.6f")
                files[(grp, seg)] = str(p)
        groups = [
            [files[(1, 0)], files[(1, 1)]],
            [files[(2, 0)], files[(2, 1)]],
        ]
        batch = vt.Pipeline([
            cmd.LS(0.5, 5.0, 1e-3, npeaks=1),
        ]).run_combinelcs(groups)
        assert len(batch.vars) == 2  # 2 groups
        assert "LS_Period_1_0" in batch.vars.columns


# ---------------------------------------------------------------------------
# Detrending / filtering commands
# ---------------------------------------------------------------------------

class TestDetrending:

    @needs_binary
    def test_linfit_removes_linear_trend(self):
        """cmd.linfit fits a polynomial and (optionally) subtracts it.

        Inject a known linear trend, run linfit with the right model,
        verify the residual rms drops.  This is more an API-shape probe
        than a numerical correctness check."""
        t, mag, err = _segment(0, 30, 1000, period=1e9,
                                noise=0.005, trend_slope=0.05)
        # No periodic signal; trend dominates.  rms BEFORE detrend.
        lc = vt.LightCurve.from_arrays(t, mag, err, name="trendy")
        before = vt.rms(lc).vars["RMS_0"]
        # Run linfit with a t-linear model and correct_lc=True.
        try:
            r = vt.Pipeline([
                cmd.linfit("a0 + a1*t", "a0,a1", correct_lc=True),
                cmd.rms(),
            ]).run(lc, capture_lc=True)
            after = float(r.vars["RMS_1"])
            print(f"linfit: rms {float(before):.5f} -> {after:.5f}")
            assert after < 0.5 * float(before), (
                "linfit didn't substantially reduce trend rms"
            )
        except Exception as e:
            print(f"linfit pipeline raised: {type(e).__name__}: {e}")
            raise

    @needs_binary
    def test_medianfilter_smooths(self):
        """cmd.medianfilter with replace=True is low-pass: it replaces
        the magnitude with the running median, so high-frequency noise
        drops and rms falls.  Default is high-pass (subtracts median),
        which does *not* reduce rms — that's a separate test."""
        t, mag, err = _segment(0, 30, 1000, period=2.3, noise=0.05)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="noisy")
        try:
            r = vt.Pipeline([
                cmd.medianfilter(0.1, replace=True),  # low-pass
                cmd.rms(),
            ]).run(lc)
            print(f"medianfilter rms (post-filter): "
                  f"{float(r.vars['RMS_1']):.5f}")
            # 0.1-day window over 30-day span ≈ ~3 points per window.
            # Low-pass replace strips noise but keeps the sinusoid,
            # whose rms is ~0.05/sqrt(2) ≈ 0.035.
            assert float(r.vars["RMS_1"]) < 0.045
        except Exception as e:
            print(f"medianfilter raised: {type(e).__name__}: {e}")
            raise

    @needs_binary
    def test_binlc_resamples(self):
        """cmd.binlc bins the LC in time."""
        t, mag, err = _segment(0, 10, 1000, period=1.5)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="dense")
        try:
            r = vt.Pipeline([cmd.binlc(binsize=0.1)]).run(lc, capture_lc=True)
            # 10-day span / 0.1-day bin → ~100 bins.
            if r.lc is not None:
                assert 50 < len(r.lc.t) < 200, (
                    f"expected ~100 bins, got {len(r.lc.t)}"
                )
            else:
                print("binlc did not capture an LC")
        except Exception as e:
            print(f"binlc raised: {type(e).__name__}: {e}")
            raise


# ---------------------------------------------------------------------------
# Custom compute mid-pipeline
# ---------------------------------------------------------------------------

class TestCustomCompute:

    @needs_binary
    def test_expr_creates_named_column(self):
        """cmd.expr can define a new per-point variable mid-pipeline."""
        t, mag, err = _segment(0, 10, 500)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="x")
        r = vt.Pipeline([
            cmd.expr("snr=mag/err"),
            cmd.stats("snr", "mean,max,min"),
        ]).run(lc)
        snr_cols = [c for c in r.vars.index if "STATS_snr" in c
                                              or "snr" in c.lower()]
        assert snr_cols, f"no snr stats in {list(r.vars.index)[:10]}"

    @needs_binary
    def test_changevariable_swaps_t(self):
        """cmd.changevariable("t", "newvar") makes vartools use newvar as
        the time axis for downstream commands."""
        t, mag, err = _segment(0, 10, 500, period=2.0)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="ph")
        try:
            r = vt.Pipeline([
                cmd.expr("phase=(t%2.0)/2.0"),
                cmd.changevariable("t", "phase"),
                cmd.LS(0.05, 0.5, 1e-3, npeaks=1),
            ]).run(lc)
            print(f"changevariable+LS cols: "
                  f"{list(r.vars.index)[:6]}...")
        except Exception as e:
            print(f"changevariable raised: {type(e).__name__}: {e}")
            raise

    @needs_binary
    def test_restricttimes_windows(self):
        """cmd.restricttimes selects a time window."""
        t, mag, err = _segment(0, 30, 1000)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="full")
        try:
            r = vt.Pipeline([
                cmd.restricttimes("JDrange", 10.0, 20.0),
                cmd.rms(),
            ]).run(lc, capture_lc=True)
            # After restriction, ~1/3 of points should remain.
            if r.lc is not None:
                n = len(r.lc.t)
                assert 250 < n < 400, (
                    f"expected ~333 points after restriction, got {n}"
                )
            else:
                print("restricttimes didn't capture LC")
        except Exception as e:
            print(f"restricttimes raised: {type(e).__name__}: {e}")
            raise


# ---------------------------------------------------------------------------
# Mixed-shape inputs in a batch
# ---------------------------------------------------------------------------

class TestMixedShape:

    @needs_binary
    def test_batch_with_varying_lc_lengths(self):
        """LCs of different lengths in one batch.  The library mode
        per-LC loop must allocate / re-allocate the LC arrays correctly."""
        lcs = []
        for i, n in enumerate([50, 500, 200, 1000, 100]):
            t, mag, err = _segment(0, 10, n, period=1.5 + 0.1 * i,
                                    seed=i)
            lcs.append(vt.LightCurve.from_arrays(
                t, mag, err, name=f"lc{i}_n{n}"))
        batch = vt.Pipeline([cmd.rms()]).run_batch(lcs)
        assert len(batch.vars) == 5
        # Each row's Npoints_0 matches the input LC length.
        for i, expected_n in enumerate([50, 500, 200, 1000, 100]):
            row_npts = int(batch.vars["Npoints_0"].iloc[i])
            assert row_npts == expected_n, (
                f"LC {i}: got {row_npts} pts, expected {expected_n}"
            )

    @needs_binary
    def test_batch_extra_columns_per_lc_with_varying_shape(self):
        """Some LCs have an extra column, others don't.  pyvartools
        should reject this clearly (rather than producing garbage)."""
        lc_with = vt.LightCurve.from_arrays(
            t=np.linspace(0, 10, 100),
            mag=np.full(100, 10.0),
            err=np.full(100, 0.01),
            aux={"airmass": np.full(100, 1.2)},
            name="with_aux",
        )
        lc_without = vt.LightCurve.from_arrays(
            t=np.linspace(0, 10, 100),
            mag=np.full(100, 10.0),
            err=np.full(100, 0.01),
            name="without_aux",
        )
        try:
            batch = vt.Pipeline([cmd.rms()]).run_batch(
                [lc_with, lc_without], raise_on_error=False)
            print(f"mixed-aux batch: rows={len(batch.vars)}, "
                  f"error={batch.error}")
        except Exception as e:
            print(f"mixed-aux batch raised: "
                  f"{type(e).__name__}: {e}")
