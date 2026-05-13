"""Exoplanet-hunter scenarios — a fourth new-user persona.

This user is hunting for transiting exoplanets with BLS, injection-recovery
to characterize their sensitivity, transit fitting on recovered candidates,
and astropy/lightkurve interop for ingesting TESS-style data.

Probes corners untouched by the v1/v2/v3 passes:
  - BLS-specific output and back-references.
  - Injection-recovery flow (cmd.injecttransit → cmd.BLS).
  - Transit fitting (cmd.MandelAgolTransit).
  - cmd.Phase back-references ("bls", "ls").
  - astropy.timeseries / lightkurve interop the docs advertise.
  - Multi-stage chains where one command's output drives another.

Run with:
    conda run -n pyvartools python -m pytest python/tests/test_doc_user_scenarios_v4.py -v -s
"""

from __future__ import annotations

import math
import os
import sys
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


def _transit_lc(period=2.3, depth=0.01, duration_frac=0.04, npts=2000,
                noise=0.001, seed=0):
    """A noisy LC with a boxcar transit signal injected at *period*."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0, 30, npts)
    mag = np.full(npts, 10.0)
    # Boxcar transits: brighter source dims by `depth` mag for
    # `duration_frac * period` around each transit centre.
    phase = ((t - 0.3) % period) / period
    in_transit = (phase < duration_frac) | (phase > 1.0 - duration_frac)
    mag[in_transit] += depth
    mag += rng.normal(0, noise, npts)
    err = np.full(npts, noise)
    return vt.LightCurve.from_arrays(t, mag, err, name="planet_host")


# ---------------------------------------------------------------------------
# BLS-driven transit search
# ---------------------------------------------------------------------------

class TestBLSSearch:

    def test_bls_default_args_clear_error(self):
        """``cmd.BLS(minper, maxper)`` with all-default kwargs used to emit
        an argv vartools rejected with an opaque message (``r ... optimal``
        is invalid).  pyvartools now raises at construction time with a
        clear message pointing at both fix paths (density_mode= or
        nfreq=/df=)."""
        with pytest.raises(ValueError) as excinfo:
            cmd.BLS(0.5, 5.0)
        msg = str(excinfo.value)
        assert "density_mode" in msg
        assert "nfreq" in msg or "df" in msg

    @needs_binary
    def test_bls_recovers_injected_period(self):
        """BLS should recover an injected period to within ~0.5%."""
        lc = _transit_lc(period=2.3)
        r = vt.Pipeline([cmd.BLS(0.5, 5.0, npeaks=1, nfreq=20000)]).run(lc)
        # Find any BLS-period column.
        period_cols = [c for c in r.vars.index if "BLS_Period_1" in c]
        assert period_cols, f"no BLS_Period_1 in {list(r.vars.index)[:10]}"
        recovered = float(r.vars[period_cols[0]])
        rel_err = abs(recovered - 2.3) / 2.3
        assert rel_err < 0.005, (
            f"BLS missed the injected period: got {recovered:.4f}, "
            f"expected 2.3 (rel err {rel_err:.3%})"
        )

    @needs_binary
    def test_bls_save_periodogram_captured(self):
        """save_periodogram=True populates result.files."""
        lc = _transit_lc(period=2.3, npts=1000)
        r = vt.Pipeline([
            cmd.BLS(0.5, 5.0, npeaks=1, nfreq=10000,
                    save_periodogram=True),
        ]).run(lc)
        pgram_keys = [k for k in r.files
                      if "periodogram" in k.lower() and "BLS" in k]
        assert pgram_keys, f"no BLS periodogram in files: {list(r.files)}"

    @needs_binary
    def test_bls_then_phase_via_backref(self):
        """lc.BLS(...).Phase(period='bls') should resolve to BLS's
        best period across the chain boundary."""
        lc = _transit_lc(period=2.3, npts=1000)
        r = lc.BLS(0.5, 5.0, npeaks=1, nfreq=10000).Phase(period="bls")
        assert r.lc is not None
        # Phase-folded t in [0, 1).
        assert 0.0 <= float(r.lc.t.min()) < 1.0
        assert float(r.lc.t.max()) < 1.0

    @needs_binary
    def test_two_bls_instances_have_distinct_suffixes(self):
        """Two BLS commands in the same pipeline get distinct column
        suffixes — the user shouldn't have to worry about collisions."""
        lc = _transit_lc(period=2.3, npts=1000)
        r = vt.Pipeline([
            cmd.BLS(0.5, 5.0, npeaks=1, nfreq=20000),
            cmd.clip(5.0),
            cmd.BLS(0.5, 5.0, npeaks=1, nfreq=20000),
        ]).run(lc)
        idx = list(r.vars.index)
        # Suffix _0 from first BLS, _2 from second (after clip at index 1).
        period_0 = [c for c in idx if "BLS_Period_1" in c and c.endswith("_0")]
        period_2 = [c for c in idx if "BLS_Period_1" in c and c.endswith("_2")]
        assert period_0 and period_2, (
            f"missing BLS_Period_1_0 / BLS_Period_1_2: {idx[:10]}"
        )


# ---------------------------------------------------------------------------
# Injection-recovery (synthetic transit + BLS to detect)
# ---------------------------------------------------------------------------

# Injection-recovery and MandelAgolTransit speculative tests were
# trimmed: the cmd.injecttransit / cmd.MandelAgolTransit wrappers have
# many CLI-mirroring kwargs whose docstrings differ from the actual
# constructor signature in ways that need a focused investigation, not
# a drive-by test from a new-user pass.  The BLS-defaults bug fixed in
# this commit is what this exercise meant to surface.


# ---------------------------------------------------------------------------
# astropy / lightkurve interop (the docs advertise these — does it work?)
# ---------------------------------------------------------------------------

class TestEcosystemInterop:

    def test_from_timeseries_round_trip(self):
        """vt.LightCurve.from_timeseries / .to_timeseries should round-trip
        a LC's data faithfully (within float-conversion noise)."""
        try:
            from astropy.timeseries import TimeSeries
            from astropy.time import Time
            from astropy import units as u
        except ImportError:
            pytest.skip("astropy not available")

        t_vals = np.linspace(0, 10, 100)
        mag_vals = 10.0 + 0.05 * np.sin(2 * np.pi * t_vals / 1.7)
        err_vals = np.full(100, 0.005)
        ts = TimeSeries(time=Time(t_vals, format="jd"),
                        data={"mag": mag_vals, "err": err_vals})

        lc = vt.LightCurve.from_timeseries(ts)
        np.testing.assert_allclose(lc.t, t_vals, atol=1e-9)
        np.testing.assert_allclose(lc.mag, mag_vals, atol=1e-9)

        ts2 = lc.to_timeseries()
        # ts2.time may be either Time or a column depending on shape; the
        # mag column should be there and match.
        assert "mag" in ts2.colnames
        np.testing.assert_allclose(np.asarray(ts2["mag"]), mag_vals,
                                   atol=1e-9)

    @needs_binary
    def test_from_timeseries_and_run_pipeline(self):
        """A LC built via TimeSeries should run through a Pipeline just
        like any other LC."""
        try:
            from astropy.timeseries import TimeSeries
            from astropy.time import Time
        except ImportError:
            pytest.skip("astropy not available")

        t_vals = np.linspace(0, 10, 200)
        mag_vals = 10.0 + 0.05 * np.sin(2 * np.pi * t_vals / 1.5)
        err_vals = np.full(200, 0.005)
        ts = TimeSeries(time=Time(t_vals, format="jd"),
                        data={"mag": mag_vals, "err": err_vals})
        lc = vt.LightCurve.from_timeseries(ts)
        r = vt.Pipeline([cmd.LS(0.5, 5.0, 1e-3, npeaks=1)]).run(lc)
        recovered = float(r.vars["LS_Period_1_0"])
        assert abs(recovered - 1.5) / 1.5 < 0.01

    def test_from_lightkurve_with_minimal_synthetic_lk(self):
        """lightkurve.LightCurve interop — make a minimal stub and check
        the conversion runs without crashing.  Most users come from a
        real TESS / Kepler download, so we just probe the API surface."""
        try:
            import lightkurve as lk
        except ImportError:
            pytest.skip("lightkurve not available")
        t = np.linspace(0, 30, 500)
        flux = 1.0 + 0.001 * np.sin(2 * np.pi * t / 2.3)
        flux_err = np.full(500, 0.0005)
        lklc = lk.LightCurve(time=t, flux=flux, flux_err=flux_err)
        try:
            lc = vt.LightCurve.from_lightkurve(lklc)
            assert lc.t is not None
            assert lc.mag is not None
        except Exception as e:
            print(f"from_lightkurve raised: {type(e).__name__}: {e}")
            raise


# ---------------------------------------------------------------------------
# Cross-command back-references — chained "bls" / "ls" period passing
# ---------------------------------------------------------------------------

class TestBackReferences:

    @needs_binary
    def test_ls_then_phase_in_one_pipeline(self):
        """Single-Pipeline LS → Phase using period='ls' (no chain step)."""
        t = np.linspace(0, 30, 1000)
        mag = 10.0 + 0.1 * np.sin(2 * np.pi * t / 1.8)
        err = np.full(1000, 0.005)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="periodic")
        r = vt.Pipeline([
            cmd.LS(0.5, 5.0, 1e-3, npeaks=1),
            cmd.Phase(period="ls"),
        ]).run(lc, capture_lc=True)
        # Period should be recovered, LC should be phase-folded.
        period = float(r.vars["LS_Period_1_0"])
        assert abs(period - 1.8) / 1.8 < 0.01
        assert r.lc is not None
        assert 0.0 <= float(r.lc.t.min()) < 1.0
        assert float(r.lc.t.max()) < 1.0

    @needs_binary
    def test_aov_column_naming_quirk(self):
        """**Documented surprise**: vartools' AOV command produces output
        columns named ``Period_k_N`` / ``AOV_k_N`` / ``AOV_SNR_k_N`` —
        no ``AOV_`` prefix on the period column itself, unlike LS which
        emits ``LS_Period_k_N``.  A new user who reads the period-finder
        page and assumes a uniform ``<COMMAND>_Period_*`` convention
        will be tripped up.  Worth a doc note OR a wrapper that adds
        the prefix.  This test pins the current quirk."""
        lc = _transit_lc(period=2.1, npts=1000, noise=0.001)
        r = lc.aov(0.5, 5.0, 0.1, 0.01, npeaks=1)
        # No "AOV_Period_*" column — just "Period_*".
        assert "Period_1_0" in r.vars.index
        assert not any(c.startswith("AOV_Period") for c in r.vars.index), (
            f"AOV column naming changed: {list(r.vars.index)}"
        )

    @needs_binary
    def test_aov_then_killharm(self):
        """aov → Killharm chain.  Uses the documented period back-ref
        ``period='aov'`` so the killharm step doesn't need an explicit
        period value."""
        lc = _transit_lc(period=2.1, npts=1000, noise=0.001)
        r = lc.aov(0.5, 5.0, 0.1, 0.01, npeaks=1).Killharm(period="aov")
        # Killharm output columns appear (any column starting with
        # Killharm_ is good — the exact suffix layout depends on the
        # number of harmonics fit).
        kh_cols = [c for c in r.vars.index if "Killharm" in c]
        assert kh_cols, f"no Killharm cols in {list(r.vars.index)[:8]}"


# ---------------------------------------------------------------------------
# Realism: a 20K-point LC akin to TESS sector data
# ---------------------------------------------------------------------------

class TestRealistic:

    @needs_binary
    def test_tess_like_lc(self):
        """20K-point LC with a transit signal at TESS-sector size.  BLS
        on a boxcar transit may land on a 1x, 2x, or 0.5x harmonic
        depending on noise + frequency grid; accept any of those as a
        successful recovery.  The point is just a performance + non-
        crash check at realistic survey data size."""
        lc = _transit_lc(period=2.3, npts=20000, noise=0.0005)
        r = vt.Pipeline([
            cmd.BLS(0.5, 5.0, npeaks=1, nfreq=20000),
        ]).run(lc)
        period_cols = [c for c in r.vars.index if "BLS_Period_1" in c]
        recovered = float(r.vars[period_cols[0]])
        ratios = [recovered / 2.3, 2.3 / recovered]
        assert any(abs(r - 1) < 0.01 or abs(r - 2) < 0.05
                   or abs(r - 0.5) < 0.05 for r in ratios), (
            f"20K-point BLS recovered {recovered:.4f}d, neither close to "
            f"2.3 nor a 0.5x/2x harmonic"
        )
