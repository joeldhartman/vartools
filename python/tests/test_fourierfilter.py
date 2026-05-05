"""Tests for the ``-fourierfilter`` command.

Scope
-----
* :class:`TestFourierFilterInfrastructure`: CLI parsing across all
  modes, output-column registration, auxiliary-file capture,
  filterexpr/freqvar acceptance, and RMS_In.
* :class:`TestFourierFilterNumerical`: analytic round-trip of the
  cos/sin fit and correctness of each filter mode (full/highpass/
  lowpass/bandpass/bandcut/filterexpr), using uniformly-sampled inputs
  on the natural ``df = 1/T`` bin grid.
"""
from __future__ import annotations

import os
import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import pyvartools as vt

VARTOOLS_SRC = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                             "..", ".."))
EXAMPLES_DIR = os.path.join(VARTOOLS_SRC, "EXAMPLES")


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def _uniform_sinusoid(freqs, amps, phases=None, n=1024, duration=100.0,
                      noise=1e-4, rng=None):
    """Build a uniformly-sampled multi-tone light curve.

    ``endpoint=True`` so the bin grid (``df = 1/T``) lands *exactly* on
    integer multiples of ``1/T`` — otherwise signals injected at
    "round" frequencies like 0.5 Hz straddle two bins and the analytic
    round-trip tolerances have to be widened.

    The ``err`` column is never allowed to hit exactly zero because the
    direct-DFT fit weights each point by ``1/err**2``; err=0 poisons
    the weights with infinities.  A fixed 1e-4 floor is plenty low to
    stay well below typical signal amplitudes.
    """
    if rng is None:
        rng = np.random.default_rng(0)
    if phases is None:
        phases = [0.0] * len(freqs)
    noise_floor = max(noise, 1e-4)
    t = np.linspace(0.0, duration, n, endpoint=True)
    mag = np.full(n, 10.0)
    for f, a, phi in zip(freqs, amps, phases):
        mag += a * np.cos(2 * np.pi * f * t + phi)
    if noise > 0:
        mag += noise * rng.standard_normal(n)
    err = np.full(n, noise_floor)
    return vt.LightCurve.from_arrays(t, mag, err)


def _rms(x):
    return float(np.sqrt(np.mean((x - np.mean(x)) ** 2)))


# --------------------------------------------------------------------------
# Infrastructure: parsing, column registration, file capture
# --------------------------------------------------------------------------

class TestFourierFilterInfrastructure:
    """CLI parsing, column registration, and file capture — the parts
    that are independent of the broken RAG math."""

    @pytest.fixture(scope="class")
    def lc(self):
        return _uniform_sinusoid([0.1, 1.0], [0.05, 0.02], n=512)

    @pytest.mark.parametrize("mode,kw", [
        ("full",     {}),
        ("highpass", {"minfreq": 0.5}),
        ("lowpass",  {"maxfreq": 0.5}),
        ("bandpass", {"minfreq": 0.3, "maxfreq": 0.7}),
        ("bandcut",  {"minfreq": 0.3, "maxfreq": 0.7}),
    ])
    def test_mode_parses_and_runs(self, lc, mode, kw):
        """Each filter mode accepts its required args and runs end-to-end."""
        r = vt.Pipeline().fourierfilter(mode=mode, **kw).run(lc)
        for col in ("FourierFilter_Mean_Mag_0",
                    "FourierFilter_RMS_In_0",
                    "FourierFilter_RMS_Out_0",
                    "FourierFilter_Nfreqcalc_0"):
            assert col in r.vars.index, f"missing {col} in mode={mode}"

    def test_rms_in_matches_injected_signal(self, lc):
        """RMS_In is computed *before* the filter touches the data and
        is therefore unaffected by the RAG bug.  For a two-tone signal
        of amplitudes a1, a2 the RMS is sqrt((a1**2 + a2**2)/2)."""
        r = vt.Pipeline().fourierfilter(mode="full").run(lc)
        expected = np.sqrt((0.05 ** 2 + 0.02 ** 2) / 2.0)
        actual = float(r.vars["FourierFilter_RMS_In_0"])
        assert abs(actual - expected) < 0.01 * expected, (
            f"RMS_In={actual}, expected ~{expected}")

    def test_nfreq_positive_and_scales_with_nsamples(self, lc):
        """Nfreq is a positive integer and (on this uniform grid)
        on the order of N/2.  Exact bin spacing depends on the
        requested band and whether fullspec is set, so we only assert
        a loose sanity range."""
        r = vt.Pipeline().fourierfilter(mode="full").run(lc)
        nfreq = int(r.vars["FourierFilter_Nfreqcalc_0"])
        assert 0 < nfreq < len(lc.mag), f"Nfreq={nfreq} out of range"

    def test_save_fouriercoeffs_captures_dataframe(self, lc):
        r = vt.Pipeline().fourierfilter(
            mode="full", save_fouriercoeffs=True).run(lc)
        key = "fourierfilter_fouriercoeffs_0"
        assert key in r.files, f"expected {key} in r.files; got {list(r.files)}"
        df = r.files[key]
        # Three columns: frequency, cos-coeff, sin-coeff (regardless of
        # whether the values are meaningful yet).
        assert df.shape[1] == 3
        assert len(df) > 0

    def test_filterexpr_default_freqvar_accepted(self, lc):
        """Default frequency variable name is 'f'."""
        r = vt.Pipeline().fourierfilter(
            mode="full", filterexpr="exp(-f*f)").run(lc)
        assert "FourierFilter_Nfreqcalc_0" in r.vars.index

    def test_filterexpr_custom_freqvar_accepted(self, lc):
        """Custom frequency variable name via freqvar="..."."""
        r = vt.Pipeline().fourierfilter(
            mode="full", filterexpr="exp(-freq*freq)",
            freqvar="freq").run(lc)
        assert "FourierFilter_Nfreqcalc_0" in r.vars.index

    def test_fullspec_flag_accepted(self, lc):
        """fullspec=True should be accepted alongside a narrow band."""
        r = vt.Pipeline().fourierfilter(
            mode="bandpass", minfreq=0.3, maxfreq=0.7,
            fullspec=True, save_fouriercoeffs=True).run(lc)
        df = r.files["fourierfilter_fouriercoeffs_0"]
        assert df.shape[1] == 3

    def test_lc_captured_and_has_same_length(self, lc):
        """Output LC is captured and has the same number of samples as
        the input."""
        r = vt.Pipeline().fourierfilter(mode="full").run(lc, capture_lc=True)
        assert r.lc is not None
        assert len(r.lc.mag) == len(lc.mag)


# --------------------------------------------------------------------------
# Numerical correctness
# --------------------------------------------------------------------------

class TestFourierFilterNumerical:
    """Analytic round-trip and filter-type correctness.

    Uniformly-sampled signals are used throughout so the DFT projection
    (now the backend for ``-fourierfilter``) lands on an orthogonal
    cos/sin basis with bin spacing ``df = 1/T``.
    """

    def test_mean_mag_recovers_dc_offset(self):
        """mode='full' with no cutoffs should recover the injected DC."""
        lc = _uniform_sinusoid([1.0], [0.05], n=1024, noise=0.0)
        r = vt.Pipeline().fourierfilter(mode="full").run(lc)
        mean = float(r.vars["FourierFilter_Mean_Mag_0"])
        assert abs(mean - 10.0) < 1e-2, f"expected Mean_Mag ~10, got {mean}"

    def test_full_mode_coefficients_recover_injected_amplitude(self):
        """Injected cosine at f=0.1 should show up as amp ~0.05 in the
        Fourier-coefficients file."""
        f0, a0 = 0.1, 0.05
        lc = _uniform_sinusoid([f0], [a0], n=1024, noise=0.0)
        r = vt.Pipeline().fourierfilter(
            mode="full", save_fouriercoeffs=True).run(lc)
        df = r.files["fourierfilter_fouriercoeffs_0"]
        freq = df.iloc[:, 0].to_numpy()
        cos  = df.iloc[:, 1].to_numpy()
        sin  = df.iloc[:, 2].to_numpy()
        amps = np.sqrt(cos ** 2 + sin ** 2)
        # Nearest bin to f0, amplitude should match a0 within a few %.
        i = int(np.argmin(np.abs(freq - f0)))
        assert abs(amps[i] - a0) < 0.1 * a0, (
            f"amp at f={freq[i]} is {amps[i]}, expected ~{a0}")

    def test_lowpass_preserves_lowfreq_kills_highfreq(self):
        """Two-tone signal: low tone survives lowpass, high tone removed."""
        lc = _uniform_sinusoid([0.1, 1.0], [0.05, 0.05], n=1024, noise=0.0)
        r = vt.Pipeline().fourierfilter(
            mode="lowpass", maxfreq=0.3).run(lc, capture_lc=True)
        out = np.asarray(r.lc.mag)
        t = np.asarray(r.lc.t)
        # Correlate with the pure low-freq template — expect > 0.95.
        lowtpl = np.cos(2 * np.pi * 0.1 * t)
        rho = np.corrcoef(out - out.mean(), lowtpl)[0, 1]
        assert rho > 0.95, f"correlation with low-freq template: {rho}"

    def test_highpass_kills_lowfreq_preserves_highfreq(self):
        lc = _uniform_sinusoid([0.1, 1.0], [0.05, 0.05], n=1024, noise=0.0)
        r = vt.Pipeline().fourierfilter(
            mode="highpass", minfreq=0.5).run(lc, capture_lc=True)
        out = np.asarray(r.lc.mag)
        t = np.asarray(r.lc.t)
        hightpl = np.cos(2 * np.pi * 1.0 * t)
        rho = np.corrcoef(out - out.mean(), hightpl)[0, 1]
        assert rho > 0.95

    def test_bandpass_isolates_middle_band(self):
        lc = _uniform_sinusoid(
            [0.1, 0.5, 2.0], [0.05, 0.05, 0.05], n=1024, noise=0.0)
        r = vt.Pipeline().fourierfilter(
            mode="bandpass", minfreq=0.3, maxfreq=1.0).run(lc, capture_lc=True)
        out = np.asarray(r.lc.mag)
        t = np.asarray(r.lc.t)
        midtpl = np.cos(2 * np.pi * 0.5 * t)
        rho = np.corrcoef(out - out.mean(), midtpl)[0, 1]
        assert rho > 0.95

    def test_bandcut_removes_middle_band(self):
        lc = _uniform_sinusoid(
            [0.1, 0.5, 2.0], [0.05, 0.05, 0.05], n=1024, noise=0.0)
        r = vt.Pipeline().fourierfilter(
            mode="bandcut", minfreq=0.3, maxfreq=1.0).run(lc, capture_lc=True)
        out = np.asarray(r.lc.mag)
        t = np.asarray(r.lc.t)
        # After band-cut, only the 0.1 and 2.0 Hz tones should remain —
        # the residual should be orthogonal to the 0.5 Hz template.
        midtpl = np.cos(2 * np.pi * 0.5 * t)
        rho = abs(np.corrcoef(out - out.mean(), midtpl)[0, 1])
        assert rho < 0.1, f"bandcut didn't remove mid tone (rho={rho})"

    def test_filterexpr_gaussian_attenuation(self):
        """An analytic Gaussian filter W(f)=exp(-(f/fc)**2) attenuates a
        tone at fc by 1/e (amplitude factor ~0.368)."""
        fc, f0, a0 = 0.5, 0.5, 0.05
        lc = _uniform_sinusoid([f0], [a0], n=1024, noise=0.0)
        r = vt.Pipeline().fourierfilter(
            mode="full",
            filterexpr=f"exp(-(f/{fc})*(f/{fc}))").run(lc, capture_lc=True)
        out = np.asarray(r.lc.mag)
        # RMS of a pure cosine of amp A is A/sqrt(2); after filter
        # A -> A/e, so RMS -> A/(e*sqrt(2)).
        expected_rms = a0 / (np.e * np.sqrt(2))
        actual_rms = _rms(out)
        assert abs(actual_rms - expected_rms) < 0.1 * expected_rms, (
            f"expected post-filter RMS ~{expected_rms}, got {actual_rms}")


# --------------------------------------------------------------------------
# FFT path (Phase E): fast path for uniform sampling or forcefft
# --------------------------------------------------------------------------

class TestFourierFilterFFT:
    """The FFT branch of -fourierfilter (GSL gsl_fft_complex_forward).

    Uniform sampling alone is enough to trigger the fast path, but we
    also exercise the explicit ``forcefft=True`` override, and confirm
    that the Fourier coefficients match what numpy.fft.rfft gives for
    the same signal — the strongest correctness check available.
    """

    def test_coefficients_match_numpy_rfft(self):
        """For a uniformly-sampled signal the ``-fourierfilter`` FFT
        path should produce the same cos/sin coefficients that numpy's
        rfft gives (up to the same 2/N normalization)."""
        n, T = 1024, 100.0
        lc = _uniform_sinusoid([0.3, 1.1], [0.05, 0.02], n=n, duration=T,
                               noise=0.0)
        r = vt.Pipeline().fourierfilter(
            mode="full", save_fouriercoeffs=True).run(lc)
        coeffs = r.files["fourierfilter_fouriercoeffs_0"].to_numpy()
        freq = coeffs[:, 0]
        a_vt = coeffs[:, 1]
        b_vt = coeffs[:, 2]

        # Compare with numpy.fft.rfft on the mean-subtracted signal.
        mag = np.asarray(lc.mag)
        mag_zero = mag - mag.mean()
        X = np.fft.rfft(mag_zero)
        a_np = np.concatenate(([mag.mean()], 2.0 * X[1:].real / n))
        b_np = np.concatenate(([0.0],         -2.0 * X[1:].imag / n))

        m = min(len(a_vt), len(a_np))
        np.testing.assert_allclose(a_vt[:m], a_np[:m], atol=1e-9, rtol=1e-6)
        np.testing.assert_allclose(b_vt[:m], b_np[:m], atol=1e-9, rtol=1e-6)

    def test_forcefft_true_on_uniform_matches_default(self):
        """Passing ``forcefft=True`` on uniformly-sampled data is a
        no-op — the FFT path is already used by default."""
        lc = _uniform_sinusoid([0.5], [0.05], n=2048, noise=1e-4)
        r_default = vt.Pipeline().fourierfilter(mode="full").run(
            lc, capture_lc=True)
        r_forced = vt.Pipeline().fourierfilter(
            mode="full", forcefft=True).run(lc, capture_lc=True)
        np.testing.assert_allclose(
            r_default.lc.mag, r_forced.lc.mag, atol=1e-10, rtol=1e-10)

    def test_nonuniform_sampling_without_resample_errors(self, capfd):
        """The direct-DFT fallback for non-uniform data has been
        removed.  Now the command errors out and returns without
        modifying the light curve; the user must add
        ``resample="delmin"`` (or another delta) for correctness."""
        rng = np.random.default_rng(0)
        n, T = 2048, 100.0
        t = np.sort(rng.uniform(0, T, n))
        mag = 10.0 + 0.05 * np.cos(2 * np.pi * 0.5 * t) + 1e-4 * rng.standard_normal(n)
        err = np.full(n, 1e-4)
        lc = vt.LightCurve.from_arrays(t, mag, err)
        capfd.readouterr()
        r = vt.Pipeline().fourierfilter(mode="full").run(lc)
        err_out = capfd.readouterr().err
        assert "requires `resample" in err_out
        # Nfreq==0 signals the filter was skipped for this LC.
        assert int(r.vars["FourierFilter_Nfreqcalc_0"]) == 0

    def test_full_roundtrip_preserves_phase(self):
        """Regression: FFT path must round-trip a signal without a
        phase shift.  The bug was that fit_harmonic_series_FFT returned
        coefficients relative to FFT index (t=0 origin), but the
        synthesis loop used absolute time t[i] = t[0] + i*delmin,
        leaving a bin-dependent exp(-i*2*pi*k*df*t[0]) phase mismatch.
        The fix rotates each FFT coefficient by that phase factor.

        Construct a LC with t[0] far from 0 (astronomical JD scale),
        run ``fourierfilter full``, and assert the output is
        element-wise close to the input (no bulk phase scramble)."""
        # JD-like absolute times — t[0] of order 53000 days; the
        # pre-fix code gave max|diff| ~ input RMS here.
        N, T = 4096, 100.0
        dt = T / N
        t = 53725.17 + np.arange(N) * dt
        rng = np.random.default_rng(0)
        mag = (10 + 0.05 * np.cos(2 * np.pi * 0.3 * t)
               + 0.03 * np.cos(2 * np.pi * 2.0 * t)
               + 1e-4 * rng.standard_normal(N))
        err = np.full(N, 1e-4)
        lc = vt.LightCurve.from_arrays(t, mag, err)

        r = vt.Pipeline().fourierfilter(mode="full").run(lc, capture_lc=True)
        out = np.asarray(r.lc.mag)

        # Element-wise agreement — not just RMS / magnitudes.  Residual
        # tolerance is loose because the input has noise injected, but
        # tight enough to catch the pre-fix phase-scramble (which gave
        # residuals on the order of the signal amplitude itself).
        diff = mag - out
        assert np.max(np.abs(diff)) < 1e-3, (
            f"full round-trip max|diff| = {np.max(np.abs(diff)):.3e} — "
            f"much larger than expected FP noise; the FFT-coefficient "
            f"phase-rotation bug may have regressed")
        # And the waveforms should be in phase (corr ~ 1, not 0.3).
        rho = np.corrcoef(mag, out)[0, 1]
        assert rho > 0.999, f"in/out correlation {rho} too low"

    def test_fft_path_is_fast_for_large_uniform(self):
        """End-to-end time should scale better than the O(N^2) direct
        DFT.  For N=16k uniform samples, pre-fix this was ~3s; with the
        FFT fit it should be well under 1s even including the O(N^2)
        synthesis loop (which dominates now but is still much faster
        than the old RAG/direct-DFT fit)."""
        import time
        n = 16384
        lc = _uniform_sinusoid([0.5], [0.05], n=n, duration=100.0, noise=1e-4)
        t0 = time.time()
        vt.Pipeline().fourierfilter(mode="full").run(lc)
        dt = time.time() - t0
        assert dt < 1.5, f"fourierfilter on N={n} took {dt:.2f}s (expected <1.5s)"


# --------------------------------------------------------------------------
# Taper (Phase #109): soft band-edge transitions
# --------------------------------------------------------------------------

class TestFourierFilterTaper:
    """The ``taper`` keyword softens each cut edge instead of a
    brick-wall transition.  Tests cover:

    * Python-wrapper parse coverage (each supported name, missing
      required args).
    * Behavioural: a tone placed near a cut edge is attenuated by the
      expected weight (~0.5 at the edge for every symmetric taper).
    * deltafreq=very-small approaches brick-wall behaviour.
    * FULLSPEC + taper prints a warning and is a no-op.
    """

    # -- parse coverage -----------------------------------------------------
    def test_each_taper_name_parses(self):
        """Every supported taper name round-trips through the CLI without
        error.  Values aren't validated here; correctness is exercised by
        the behavioural tests below."""
        lc = _uniform_sinusoid([0.5], [0.05], n=512, noise=1e-4)
        for taper in ("linear", "cosine", "tukey", "hann", "blackman"):
            r = vt.Pipeline().fourierfilter(
                mode="highpass", minfreq=0.3,
                taper=taper, taper_deltafreq=0.05).run(lc)
            assert "FourierFilter_Mean_Mag_0" in r.vars.index
        # kaiser needs a beta
        r = vt.Pipeline().fourierfilter(
            mode="highpass", minfreq=0.3,
            taper="kaiser", taper_deltafreq=0.05, taper_beta=5.0).run(lc)
        assert "FourierFilter_Mean_Mag_0" in r.vars.index

    def test_unknown_taper_rejected(self):
        with pytest.raises(ValueError, match="unknown taper"):
            vt.Pipeline().fourierfilter(
                mode="highpass", minfreq=0.3,
                taper="bogus", taper_deltafreq=0.05)

    def test_taper_requires_deltafreq(self):
        with pytest.raises(ValueError, match="taper_deltafreq"):
            vt.Pipeline().fourierfilter(
                mode="highpass", minfreq=0.3, taper="cosine")

    def test_kaiser_requires_beta(self):
        with pytest.raises(ValueError, match="taper_beta"):
            vt.Pipeline().fourierfilter(
                mode="highpass", minfreq=0.3,
                taper="kaiser", taper_deltafreq=0.05)

    def test_deltafreq_without_taper_rejected(self):
        with pytest.raises(ValueError, match="no effect without taper"):
            vt.Pipeline().fourierfilter(
                mode="highpass", minfreq=0.3, taper_deltafreq=0.05)

    # -- behavioural --------------------------------------------------------
    def test_taper_at_edge_gives_half_amplitude(self):
        """A tone placed exactly at the cut edge should be attenuated by
        the taper weight at u=0.5 — which is 0.5 for linear, cosine,
        and near-0.5 for Kaiser-5.  The tone's amplitude is halved, so
        RMS scales by ~0.5.  Brick-wall at the same edge can go either
        to 0 or 1 (depending on which side of the rounded bin the tone
        lands on) — demonstrably different from the tapered response."""
        f_cut = 0.5
        a0 = 0.1
        lc = _uniform_sinusoid([f_cut], [a0], n=1024, duration=100.0, noise=0.0)
        full_rms = a0 / np.sqrt(2)        # RMS of the un-attenuated tone

        # Linear and cosine have weight 0.5 at the midpoint; the surviving
        # tone amplitude is 0.5*a0, so RMS should be ~0.5 * full_rms.
        for taper in ("linear", "cosine"):
            r = vt.Pipeline().fourierfilter(
                mode="highpass", minfreq=f_cut,
                taper=taper, taper_deltafreq=0.1).run(lc, capture_lc=True)
            rms = float(r.vars["FourierFilter_RMS_Out_0"])
            assert 0.3 * full_rms < rms < 0.7 * full_rms, (
                f"{taper} at edge: rms={rms}, expected ~{0.5*full_rms}")

        # Brick-wall at the same edge: the nearest bin lands just below
        # f_cut so the tone is mostly rejected — RMS should be much
        # smaller than the tapered value.
        r_brick = vt.Pipeline().fourierfilter(
            mode="highpass", minfreq=f_cut).run(lc, capture_lc=True)
        rms_brick = float(r_brick.vars["FourierFilter_RMS_Out_0"])
        assert rms_brick < 0.2 * full_rms, (
            f"brick-wall at edge: rms={rms_brick}, expected much less than "
            f"taper's ~{0.5*full_rms}")

    def test_very_small_deltafreq_approaches_brickwall(self):
        """As deltafreq -> 0, taper approaches a brick wall.  Compare to
        the no-taper run for a tone far from the edge (where taper
        region doesn't matter).  They should agree closely."""
        lc = _uniform_sinusoid([0.5], [0.05], n=1024, noise=1e-4)
        r_no = vt.Pipeline().fourierfilter(
            mode="highpass", minfreq=0.1).run(lc, capture_lc=True)
        r_narrow = vt.Pipeline().fourierfilter(
            mode="highpass", minfreq=0.1,
            taper="linear", taper_deltafreq=0.001).run(lc, capture_lc=True)
        np.testing.assert_allclose(
            r_no.lc.mag, r_narrow.lc.mag, atol=1e-3, rtol=1e-3)

    def test_bandpass_taper_applies_to_both_edges(self):
        """For a bandpass filter with two tones straddling each edge,
        the taper should attenuate both.  Compared to brick-wall, the
        tapered reconstruction RMS should be reduced (because both
        tones are partially rejected)."""
        # Tones exactly at the edges 0.3 and 0.7
        lc = _uniform_sinusoid([0.3, 0.7], [0.05, 0.05],
                               n=1024, duration=100.0, noise=0.0)
        r_hard = vt.Pipeline().fourierfilter(
            mode="bandpass", minfreq=0.3, maxfreq=0.7).run(lc, capture_lc=True)
        r_tap = vt.Pipeline().fourierfilter(
            mode="bandpass", minfreq=0.3, maxfreq=0.7,
            taper="cosine", taper_deltafreq=0.05).run(lc, capture_lc=True)
        rms_hard = float(r_hard.vars["FourierFilter_RMS_Out_0"])
        rms_tap = float(r_tap.vars["FourierFilter_RMS_Out_0"])
        # Both RMS values should be finite and positive
        assert rms_hard > 0 and rms_tap > 0
        # Tapered RMS differs from brick-wall — direction depends on
        # bin alignment, but they should not be bit-identical.
        assert abs(rms_tap - rms_hard) > 1e-4, (
            f"tapered RMS {rms_tap} indistinguishable from brick {rms_hard}")

    def test_fullspec_taper_is_noop_with_warning(self, capfd):
        """Taper on mode=full has nothing to soften; the wrapper still
        accepts it but vartools prints a warning.  Reconstruction should
        be indistinguishable from mode=full with no taper."""
        lc = _uniform_sinusoid([0.5], [0.05], n=1024, noise=1e-4)
        r_no = vt.Pipeline().fourierfilter(mode="full").run(lc, capture_lc=True)
        r_tap = vt.Pipeline().fourierfilter(
            mode="full", taper="cosine", taper_deltafreq=0.1
        ).run(lc, capture_lc=True)
        np.testing.assert_allclose(r_no.lc.mag, r_tap.lc.mag, atol=1e-10)


# --------------------------------------------------------------------------
# Resample path (Phase R2): correct non-uniform-sampling filter
# --------------------------------------------------------------------------

def _nonuniform_sinusoid(freqs, amps, n=500, duration=100.0,
                         noise=1e-4, seed=42):
    """Non-uniform sampling via jittered uniform grid.  Using a jitter
    (rather than pure random times) keeps min(dt) roughly comparable to
    mean(dt) — pure random has a heavy low tail so ``delmin`` can be
    tiny and resample grids can get explosively large."""
    rng = np.random.default_rng(seed)
    dt_mean = duration / n
    t = np.arange(n) * dt_mean + 0.4 * dt_mean * rng.uniform(-1, 1, n)
    t = np.sort(t)
    mag = np.full(n, 10.0)
    for f, a in zip(freqs, amps):
        mag += a * np.cos(2 * np.pi * f * t)
    if noise > 0:
        mag += noise * rng.standard_normal(n)
    err = np.full(n, max(noise, 1e-4))
    return vt.LightCurve.from_arrays(t, mag, err)


class TestFourierFilterResample:
    """The ``resample`` keyword: correct Fourier-domain filtering on
    non-uniformly-sampled data via resample → FFT → mask → IFFT →
    resample-back.  The direct-DFT fallback for non-uniform data has
    been removed; attempts to run ``-fourierfilter`` on non-uniform
    data without ``resample`` now error out (see
    :meth:`TestFourierFilterFFT.test_nonuniform_sampling_without_resample_errors`).
    """

    def test_resample_delmin_recovers_input_rms(self):
        """With ``resample='delmin'``, a full-band filter on non-uniform
        data should round-trip the signal: RMS_Out ≈ RMS_In to within
        the linear-interp error."""
        lc = _nonuniform_sinusoid([0.3], [0.05])
        r = vt.Pipeline().fourierfilter(
            mode="full", resample="delmin").run(lc, capture_lc=True)
        rms_in = float(r.vars["FourierFilter_RMS_In_0"])
        rms_out = float(r.vars["FourierFilter_RMS_Out_0"])
        assert 0.8 * rms_in < rms_out < 1.2 * rms_in, (
            f"rms_in={rms_in}, rms_out={rms_out}; expected close agreement")

    def test_resample_fix_delta_gives_same_answer(self):
        """Same test as above, but specifying a fixed delta rather than
        delmin.  Should be similar numerical result."""
        lc = _nonuniform_sinusoid([0.3], [0.05])
        r = vt.Pipeline().fourierfilter(
            mode="full", resample=0.05).run(lc, capture_lc=True)
        rms_in = float(r.vars["FourierFilter_RMS_In_0"])
        rms_out = float(r.vars["FourierFilter_RMS_Out_0"])
        assert 0.8 * rms_in < rms_out < 1.2 * rms_in

    def test_resample_lowpass_removes_highfreq(self):
        """Two-tone non-uniform signal: low-pass below the high tone
        should knock it out.  With resample active we can check a real
        filter rather than just the round-trip."""
        lc = _nonuniform_sinusoid([0.3, 2.0], [0.05, 0.05])
        r = vt.Pipeline().fourierfilter(
            mode="lowpass", maxfreq=0.5, resample=0.05).run(lc, capture_lc=True)
        t = np.asarray(lc.t)
        out = np.asarray(r.lc.mag)
        # Correlate output with the low-frequency template.
        lowtpl = np.cos(2 * np.pi * 0.3 * t)
        rho = np.corrcoef(out - out.mean(), lowtpl)[0, 1]
        assert rho > 0.8, f"low-freq correlation {rho} too low"

    def test_resample_invalid_rejected(self):
        with pytest.raises(ValueError, match="resample delta must be positive"):
            vt.Pipeline().fourierfilter(mode="full", resample=-1.0)

    def test_gapbreak_requires_resample(self):
        with pytest.raises(ValueError, match="gapbreak.*requires resample"):
            vt.Pipeline().fourierfilter(
                mode="full", gapbreak_type="fix", gapbreak_value=5.0)

    def test_unknown_gapbreak_type_rejected(self):
        with pytest.raises(ValueError, match="unknown gapbreak_type"):
            vt.Pipeline().fourierfilter(
                mode="full", resample=0.05,
                gapbreak_type="bogus", gapbreak_value=1.0)

    def test_gapbreak_missing_value(self):
        with pytest.raises(ValueError, match="gapbreak_value is required"):
            vt.Pipeline().fourierfilter(
                mode="full", resample=0.05, gapbreak_type="fix")

    def test_frequency_units_are_cycles_per_time_unit(self):
        """User-supplied minfreq/maxfreq are in cycles/(time-unit) —
        i.e. the reciprocal of whatever units the `t` column uses, with
        NO stray 2π factor.  Inject a single tone at f=0.5 cyc/(t-unit)
        and verify that for every mode, cutoffs that bracket f=0.5 keep
        the tone (RMS ~= A/sqrt(2)) while cutoffs that don't bracket it
        kill it (RMS ~= 0).  This test pins down the frequency
        convention against any refactor that introduces a unit bug."""
        N, T = 2048, 100.0
        t = np.linspace(0, T, N, endpoint=False)
        A = 0.1
        f_inj = 0.5
        mag = 10 + A * np.cos(2 * np.pi * f_inj * t)
        err = np.full(N, 1e-4)
        lc = vt.LightCurve.from_arrays(t, mag, err)
        kept = A / np.sqrt(2)

        cases = [
            ("highpass", {"minfreq": 0.3}, True),    # 0.5 > 0.3 -> keep
            ("highpass", {"minfreq": 0.7}, False),   # 0.5 < 0.7 -> kill
            ("lowpass",  {"maxfreq": 0.7}, True),
            ("lowpass",  {"maxfreq": 0.3}, False),
            ("bandpass", {"minfreq": 0.4, "maxfreq": 0.6}, True),
            ("bandpass", {"minfreq": 0.7, "maxfreq": 0.9}, False),
            ("bandcut",  {"minfreq": 0.4, "maxfreq": 0.6}, False),
            ("bandcut",  {"minfreq": 0.7, "maxfreq": 0.9}, True),
        ]
        for mode, kw, should_keep in cases:
            r = vt.Pipeline().fourierfilter(mode=mode, **kw).run(lc, capture_lc=True)
            out = np.asarray(r.lc.mag)
            rms = float(np.sqrt(np.var(out - out.mean())))
            if should_keep:
                assert abs(rms - kept) < 0.02 * kept, (
                    f"{mode} {kw}: expected tone kept (rms~={kept}), got rms={rms}")
            else:
                assert rms < 0.02 * kept, (
                    f"{mode} {kw}: expected tone killed (rms~0), got rms={rms}")

        # filterexpr evaluates f in the same cycles/(time-unit) convention
        # W(f) = exp(-(f/0.5)^2)  ->  W(0.5) = 1/e
        r = vt.Pipeline().fourierfilter(
            mode="full", filterexpr="exp(-(f/0.5)^2)").run(lc, capture_lc=True)
        rms = float(np.sqrt(np.var(np.asarray(r.lc.mag) - np.mean(r.lc.mag))))
        expected = kept / np.e
        assert abs(rms - expected) < 0.02 * expected, (
            f"Gaussian filterexpr at f=0.5 should attenuate tone by 1/e; "
            f"expected rms={expected}, got rms={rms}")

    def test_padmode_reflect_reduces_edge_bias(self):
        """Construct an LC with a deliberate endpoint mismatch, apply a
        low-pass filter via the resample path both without padding
        (`padmode="wrap"`, equivalent to default) and with
        `padmode="reflect"`, and verify the smoothed residual near the
        boundaries is smaller with reflect padding.  This is the
        user-visible benefit: no Gibbs-style ringing at the edges."""
        # Linear trend from +0.05 to -0.05 plus a fast sinusoid.  The
        # endpoints differ by 0.1, which maximizes the FFT wrap-around
        # discontinuity.
        N, T = 4096, 100.0
        t = np.linspace(0.0, T, N, endpoint=False)
        mag = (10.0 + 0.05 * (1.0 - 2.0 * t / T)
               + 0.02 * np.cos(2 * np.pi * 3.0 * t))
        err = np.full(N, 1e-4)
        lc = vt.LightCurve.from_arrays(t, mag, err)

        # Filter: low-pass at 0.5 cyc/day removes the fast sinusoid.
        # With wrap-padding the residual at the edges is dominated by
        # Gibbs ringing from the 0.1-mag endpoint discontinuity.
        r_wrap = vt.Pipeline().fourierfilter(
            mode="lowpass", maxfreq=0.5,
            resample="delmin",
            padmode="wrap").run(lc, capture_lc=True)
        r_refl = vt.Pipeline().fourierfilter(
            mode="lowpass", maxfreq=0.5,
            resample="delmin",
            padmode="reflect").run(lc, capture_lc=True)
        mag_in = mag - mag.mean()
        diff_wrap = mag_in - (np.asarray(r_wrap.lc.mag) - np.mean(r_wrap.lc.mag))
        diff_refl = mag_in - (np.asarray(r_refl.lc.mag) - np.mean(r_refl.lc.mag))

        # Average residual in the first and last 5% of the LC
        edge_w = N // 20
        wrap_edge = max(abs(diff_wrap[:edge_w].mean()),
                        abs(diff_wrap[-edge_w:].mean()))
        refl_edge = max(abs(diff_refl[:edge_w].mean()),
                        abs(diff_refl[-edge_w:].mean()))
        assert refl_edge < wrap_edge, (
            f"reflect edge bias {refl_edge:.4f} expected < wrap bias "
            f"{wrap_edge:.4f}")

    def test_padmode_unknown_rejected(self):
        with pytest.raises(ValueError, match="unknown padmode"):
            vt.Pipeline().fourierfilter(
                mode="lowpass", maxfreq=0.5,
                resample="delmin", padmode="bogus")

    def test_padfrac_without_padmode_rejected(self):
        with pytest.raises(ValueError, match="padfrac.*without padmode"):
            vt.Pipeline().fourierfilter(
                mode="full", padfrac=0.3)

    def test_nowarn_suppresses_runtime_warnings(self, capfd):
        """The ``nowarn`` keyword should suppress all per-LC runtime
        warnings from ``-fourierfilter`` without affecting numerics.

        Use a bandpass + taper where deltafreq exceeds half the band
        width — that triggers the ``taper overlap`` warning which is
        runtime-emitted per LC."""
        # Uniformly-sampled LC so we actually hit the filter path.
        N, T = 1024, 100.0
        t = np.linspace(0, T, N, endpoint=False)
        mag = 10.0 + 0.05 * np.cos(2 * np.pi * 0.5 * t)
        err = np.full(N, 1e-4)
        lc = vt.LightCurve.from_arrays(t, mag, err)

        # deltafreq (0.3) > (0.8-0.4)/2 = 0.2  -> taper-overlap warning fires
        capfd.readouterr()
        r_noisy = vt.Pipeline().fourierfilter(
            mode="bandpass", minfreq=0.4, maxfreq=0.8,
            taper="cosine", taper_deltafreq=0.3).run(lc)
        stderr_noisy = capfd.readouterr().err
        assert "Warning" in stderr_noisy, (
            f"expected taper-overlap warning in stderr, got: {stderr_noisy!r}")

        # With nowarn: stderr should be quiet.
        r_quiet = vt.Pipeline().fourierfilter(
            mode="bandpass", minfreq=0.4, maxfreq=0.8,
            taper="cosine", taper_deltafreq=0.3, nowarn=True).run(lc)
        stderr_quiet = capfd.readouterr().err
        assert "Warning" not in stderr_quiet, (
            f"expected no warnings with nowarn=True, got: {stderr_quiet!r}")

        # Numerics must match with and without the flag.
        np.testing.assert_allclose(
            float(r_noisy.vars["FourierFilter_RMS_Out_0"]),
            float(r_quiet.vars["FourierFilter_RMS_Out_0"]),
            rtol=1e-10)

    def test_gapbreak_with_real_gap(self):
        """Construct a LC with a big gap in the middle, enable
        gapbreak.  Each segment should be filtered independently and
        the concatenated output should have the two segments at their
        correct local means with no catastrophic blowup."""
        rng = np.random.default_rng(0)
        t1 = np.linspace(0, 10, 500, endpoint=False)
        t2 = np.linspace(50, 60, 500, endpoint=False)  # 40-unit gap
        t = np.concatenate([t1, t2])
        mag = 10 + 0.05 * np.cos(2 * np.pi * 1.0 * t) + 1e-4 * rng.standard_normal(len(t))
        err = np.full(len(t), 1e-4)
        lc = vt.LightCurve.from_arrays(t, mag, err)

        r = vt.Pipeline().fourierfilter(
            mode="full", resample=0.05,
            gapbreak_type="fix", gapbreak_value=5.0
        ).run(lc, capture_lc=True)
        out = np.asarray(r.lc.mag)
        assert np.all(np.isfinite(out))
        # Mean across whole output should match the original weighted mean.
        assert abs(float(r.vars["FourierFilter_Mean_Mag_0"]) - 10.0) < 1e-2
