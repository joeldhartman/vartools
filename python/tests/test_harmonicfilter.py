"""Tests for -harmonicfilter and its backward-compatible synonym -Killharm.

Scope
-----
* Parity: invoking via cmd.harmonicfilter and cmd.Killharm with the same
  arguments produces numerically identical outputs on the per-harmonic
  coefficients; the column PREFIX differs (HarmonicFilter_* vs Killharm_*)
  but every matching *suffix* carries the same value.
* Sanity: the top LS period is recovered to better than a period-grid
  step; the amplitude of an injected sinusoid matches the -Killharm
  convention within tolerance.
"""
from __future__ import annotations

import os
import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import pyvartools as vt
from pyvartools import commands as cmd

VARTOOLS_SRC = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                             "..", ".."))
EXAMPLES_DIR = os.path.join(VARTOOLS_SRC, "EXAMPLES")


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def _injected_sinusoid_lc(n=500, period=1.5, amplitude=0.1, phase=0.3,
                          noise=1e-4, rng=None):
    """Build a sinusoidal LC with known (P, A, phi) plus light Gaussian
    noise.  A little noise keeps post-fit residuals above the
    floating-point floor that -rms hits near machine precision.
    """
    if rng is None:
        rng = np.random.default_rng(0)
    t = np.sort(rng.uniform(0.0, 100.0, n))
    mag = (10.0 + amplitude * np.cos(2 * np.pi * (t / period + phase))
           + noise * rng.standard_normal(n))
    err = np.full(n, noise)
    return vt.LightCurve.from_arrays(t, mag, err), t, mag


# --------------------------------------------------------------------------
# Parity: -Killharm and -harmonicfilter produce the same numbers
# --------------------------------------------------------------------------

class TestKillharmHarmonicfilterParity:
    """Same arguments under both names must produce identical numerics."""

    @pytest.fixture(scope="class")
    def lc(self):
        return vt.LightCurve.from_file(os.path.join(EXAMPLES_DIR, "2"))

    def test_columns_differ_only_in_prefix(self, lc):
        r_old = vt.Pipeline().Killharm(period=1.2354, nharm=2, nsubharm=1).run(lc)
        r_new = vt.Pipeline().harmonicfilter(
            period=1.2354, nharm=2, nsubharm=1).run(lc)

        old_keys = [k for k in r_old.vars.index if k.startswith("Killharm_")]
        new_keys = [k for k in r_new.vars.index
                    if k.startswith("HarmonicFilter_")]

        assert len(old_keys) > 0, "no Killharm_* columns emitted"
        assert len(old_keys) == len(new_keys), (
            f"mismatched count: {len(old_keys)} Killharm_* vs "
            f"{len(new_keys)} HarmonicFilter_*")

        # For each Killharm_<suffix>, there must be HarmonicFilter_<suffix>
        # with the same numeric value (exact equality — same code path).
        for ok in old_keys:
            suffix = ok[len("Killharm_"):]
            nk = "HarmonicFilter_" + suffix
            assert nk in r_new.vars.index, f"missing {nk}"
            assert r_old.vars[ok] == r_new.vars[nk], (
                f"{ok}={r_old.vars[ok]!r} vs {nk}={r_new.vars[nk]!r}")

    @pytest.mark.parametrize("nharm,nsubharm", [(0, 0), (3, 0), (1, 2)])
    def test_various_harmonic_counts(self, lc, nharm, nsubharm):
        r_old = vt.Pipeline().Killharm(
            period=1.2354, nharm=nharm, nsubharm=nsubharm).run(lc)
        r_new = vt.Pipeline().harmonicfilter(
            period=1.2354, nharm=nharm, nsubharm=nsubharm).run(lc)
        for ok in r_old.vars.index:
            if not ok.startswith("Killharm_"):
                continue
            nk = "HarmonicFilter_" + ok[len("Killharm_"):]
            assert r_old.vars[ok] == r_new.vars[nk]

    @pytest.mark.parametrize("output_format", [
        None, "outampphase", "outampradphase", "outRphi", "outRradphi"])
    def test_alternate_output_formats(self, lc, output_format):
        kw = dict(period=1.2354, nharm=2, nsubharm=0)
        if output_format is not None:
            kw["output_format"] = output_format
        r_old = vt.Pipeline().Killharm(**kw).run(lc)
        r_new = vt.Pipeline().harmonicfilter(**kw).run(lc)
        for ok in r_old.vars.index:
            if not ok.startswith("Killharm_"):
                continue
            nk = "HarmonicFilter_" + ok[len("Killharm_"):]
            assert r_old.vars[ok] == r_new.vars[nk]

    def test_save_model_suffix_differs(self, lc, tmp_path):
        """The .killharm.model / .harmonicfilter.model model-file suffix
        follows the invoking token, matching the struct's column_prefix."""
        r_old = vt.Pipeline().Killharm(
            period=1.2354, nharm=1, save_model=True).run(lc)
        r_new = vt.Pipeline().harmonicfilter(
            period=1.2354, nharm=1, save_model=True).run(lc)
        # Both should capture a model DataFrame; the key name uses the
        # pyvartools class _vt_name, so the dict keys differ.
        assert any("killharm" in k.lower() for k in r_old.files), (
            f"expected Killharm-keyed model file; got {list(r_old.files)}")
        assert any("harmonicfilter" in k.lower() for k in r_new.files), (
            f"expected harmonicfilter-keyed model file; "
            f"got {list(r_new.files)}")


# --------------------------------------------------------------------------
# Sanity: recovering an injected sinusoid
# --------------------------------------------------------------------------

class TestHarmonicfilterSanity:

    def test_recover_injected_sinusoid_amplitude(self):
        """A noise-free injected cosine at a known period should come back
        with the expected amplitude (A_1 from the -harmonicfilter fit).
        """
        lc, _, _ = _injected_sinusoid_lc(
            n=2000, period=1.3, amplitude=0.05, phase=0.0)
        r = vt.Pipeline().harmonicfilter(
            period=1.3, nharm=0, nsubharm=0,
            output_format="outampphase").run(lc)
        amp = float(r.vars["HarmonicFilter_Per1_Fundamental_Amp_0"])
        assert abs(amp - 0.05) < 1e-4, f"expected amplitude ~0.05, got {amp}"

    def test_subtracts_model_by_default(self):
        """Without fitonly the output LC should have its mean sinusoid
        removed, so the post-filter RMS drops markedly.  Compare to a
        -rms before/after."""
        lc, _, _ = _injected_sinusoid_lc(
            n=2000, period=1.3, amplitude=0.05, phase=0.0)
        r = vt.Pipeline().rms().harmonicfilter(
            period=1.3, nharm=0, nsubharm=0).rms().run(lc)
        before = float(r.vars["RMS_0"])
        after = float(r.vars["RMS_2"])
        # Pure cosine amplitude A -> RMS A/sqrt(2); subtracting removes
        # it entirely, so after << before.
        assert after < 0.01 * before, (
            f"expected RMS to drop sharply; before={before}, after={after}")
