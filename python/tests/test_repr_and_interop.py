"""Tests for HTML reprs and lightkurve/astropy interop (Phase 2)."""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest

import pyvartools as vt
from pyvartools._binary import get_binary

try:
    get_binary()
    _HAVE_BINARY = True
except Exception:
    _HAVE_BINARY = False

try:
    import lightkurve as _lk  # noqa: F401
    _HAVE_LIGHTKURVE = True
except ImportError:
    _HAVE_LIGHTKURVE = False


_EXAMPLES = os.path.join(os.path.dirname(__file__), "..", "..", "EXAMPLES")


@pytest.fixture
def lc():
    return vt.LightCurve.from_file(os.path.join(_EXAMPLES, "2"))


# ---------------------------------------------------------------------------
# HTML reprs
# ---------------------------------------------------------------------------


class TestLightCurveHtmlRepr:
    def test_returns_html_string(self, lc):
        html = lc._repr_html_()
        assert isinstance(html, str)
        assert "<table" in html
        assert "<details" in html

    def test_includes_name_and_count(self, lc):
        html = lc._repr_html_()
        assert "LightCurve" in html
        assert str(len(lc)) in html
        assert lc.name in html

    def test_short_lc_shows_full_data(self):
        lc = vt.LightCurve.from_arrays(
            t=np.arange(5).astype(float), mag=np.arange(5).astype(float),
            err=np.full(5, 0.01), name="short")
        html = lc._repr_html_()
        # All 5 rows should be present (no ellipsis truncation)
        assert "..." not in html or html.count("...") < 4

    def test_scalars_section_when_present(self):
        lc = vt.LightCurve.from_arrays(
            t=np.arange(20).astype(float), mag=np.arange(20).astype(float),
            err=np.full(20, 0.01), name="x", scalars={"foo": 1.0, "bar": 2.0})
        html = lc._repr_html_()
        assert "scalars" in html
        assert "foo" in html and "bar" in html

    def test_no_scalars_section_when_empty(self, lc):
        html = lc._repr_html_()
        # Header section + data details — no scalars details.
        assert html.count("<details") == 1


@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools binary not available")
class TestArrayFormMultiPeak:
    """varobjs.<cmd>[i].<descriptor> returns a length-Npeaks array."""

    def test_ls_single_peak_yields_length1_array(self, lc):
        r = vt.Pipeline().LS(minp=0.1, maxp=10.0, subsample=0.1, npeaks=1).run(lc)
        ls = r.varobjs.LS[0]
        assert isinstance(ls.Period, np.ndarray)
        assert ls.Period.shape == (1,)
        # backwards-compat: scalar Period_1 still accessible
        assert np.isclose(ls.Period_1, ls.Period[0])

    def test_ls_multi_peak_array_matches_scalar_keys(self, lc):
        r = vt.Pipeline().LS(minp=0.1, maxp=10.0, subsample=0.1, npeaks=3).run(lc)
        ls = r.varobjs.LS[0]
        assert ls.Period.shape == (3,)
        # array elements equal the scalar keys
        for i in range(3):
            assert np.isclose(ls.Period[i], getattr(ls, f"Period_{i+1}"))

    def test_ls_qualified_descriptor_also_arrayified(self, lc):
        r = vt.Pipeline().LS(minp=0.1, maxp=10.0, subsample=0.1, npeaks=3).run(lc)
        ls = r.varobjs.LS[0]
        # Log10_LS_Prob_<peak>_<pos> -> stripped to Log10_LS_Prob_<peak>
        # so Log10_LS_Prob should also be a length-3 array
        assert isinstance(ls.Log10_LS_Prob, np.ndarray)
        assert ls.Log10_LS_Prob.shape == (3,)

    def test_bls_arrayifies_all_peak_descriptors(self, lc):
        r = vt.Pipeline().BLS(
            rmin=0.005, rmax=0.05, nfreq=200, nbins=200, npeaks=5,
            minper=0.5, maxper=10.0, qmin=0.01, qmax=0.1,
        ).run(lc)
        b = r.varobjs.BLS[0]
        for desc in ("Period", "Tc", "SN", "SR", "Depth", "Qtran", "SDE"):
            arr = getattr(b, desc)
            assert isinstance(arr, np.ndarray), f"{desc} is {type(arr)}"
            assert arr.shape == (5,), f"{desc} shape {arr.shape} != (5,)"
            # cross-check first scalar matches array[0]
            assert np.isclose(arr[0], getattr(b, f"{desc}_1"))


@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools binary not available")
class TestResultHtmlRepr:
    def test_result_repr_has_table_and_status(self, lc):
        r = vt.Pipeline().LS(minp=0.1, maxp=10.0, subsample=0.1, npeaks=1).run(lc)
        html = r._repr_html_()
        assert "<table" in html
        assert "Result" in html
        assert "ok" in html
        assert "LS_Period" in html

    def test_batchresult_repr_lists_n(self, lc):
        batch = vt.LightCurveBatch([lc, lc])
        br = batch.LS(minp=0.1, maxp=10.0, subsample=0.1, npeaks=1).run()
        html = br._repr_html_()
        assert "<table" in html
        assert "BatchResult" in html
        # the n=2 should be visible somewhere
        assert "n: 2" in html or "n=2" in html


# ---------------------------------------------------------------------------
# lightkurve interop
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not _HAVE_LIGHTKURVE, reason="lightkurve not installed")
class TestLightkurveInterop:
    def test_to_lightkurve_basic(self, lc):
        lklc = lc.to_lightkurve()
        # type from lightkurve
        import lightkurve as lk
        assert isinstance(lklc, lk.LightCurve)
        assert len(lklc) == len(lc)
        # data round-trips
        np.testing.assert_allclose(lklc.flux.value, lc.mag)
        np.testing.assert_allclose(lklc.flux_err.value, lc.err)

    def test_round_trip_preserves_values(self, lc):
        lklc = lc.to_lightkurve()
        lc2 = vt.LightCurve.from_lightkurve(lklc)
        assert len(lc2) == len(lc)
        np.testing.assert_allclose(lc.t, lc2.t)
        np.testing.assert_allclose(lc.mag, lc2.mag)
        np.testing.assert_allclose(lc.err, lc2.err)

    def test_from_lightkurve_carries_aux_columns(self):
        import lightkurve as lk
        from astropy.time import Time
        lklc = lk.LightCurve(
            time=Time(np.arange(10).astype(float), format="jd"),
            flux=np.linspace(1.0, 1.1, 10),
            flux_err=np.full(10, 0.001),
            quality=np.zeros(10, dtype=int),
        )
        lc = vt.LightCurve.from_lightkurve(lklc)
        # quality column should round-trip
        assert "quality" in lc._df.columns

    def test_from_lightkurve_picks_up_label_from_meta(self):
        import lightkurve as lk
        from astropy.time import Time
        lklc = lk.LightCurve(
            time=Time(np.arange(5).astype(float), format="jd"),
            flux=np.arange(5).astype(float),
            flux_err=np.full(5, 0.001),
            meta={"LABEL": "TIC123456"},
        )
        lc = vt.LightCurve.from_lightkurve(lklc)
        assert lc.name == "TIC123456"


@pytest.mark.skipif(_HAVE_LIGHTKURVE, reason="testing the missing-dep error path")
class TestLightkurveMissingError:
    def test_to_lightkurve_raises_clear_importerror(self):
        lc = vt.LightCurve.from_arrays(
            t=np.arange(3).astype(float),
            mag=np.arange(3).astype(float),
            err=np.full(3, 0.01),
        )
        with pytest.raises(ImportError, match="lightkurve is required"):
            lc.to_lightkurve()
