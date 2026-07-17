"""Tests for the ``unstitch`` USERLIB wrapper (the inverse of ``-stitch``).

Covers CLI-token assembly for both shift sources and all options,
constructor-time validation, Pipeline-method attachment, and end-to-end
round-trips against the real binary (skipped when it is unavailable).
"""
from __future__ import annotations

import os
import pathlib
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import pyvartools as vt
from pyvartools.commands import unstitch

VARTOOLS_SRC = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                             "..", ".."))
EXAMPLES_DIR = os.path.join(VARTOOLS_SRC, "EXAMPLES")

try:
    _HAVE_BINARY = os.path.isfile(vt.get_binary())
except FileNotFoundError:
    _HAVE_BINARY = False

needs_binary = pytest.mark.skipif(
    not _HAVE_BINARY, reason="vartools binary not available"
)


# ---------------------------------------------------------------------------
# CLI-token assembly
# ---------------------------------------------------------------------------

def test_tokens_in_shifts_file_minimal():
    u = unstitch("mag", "in_shifts_file", fieldlabelsvar="field",
                 starnamevar="star", in_shifts_file="shifts.txt")
    assert u._to_cli_args() == [
        "-unstitch", "mag", "in_shifts_file", "field", "star", "shifts.txt"]


def test_tokens_in_shifts_file_full():
    u = unstitch(["mag", "mag2"], "in_shifts_file", fieldlabelsvar="field",
                 starnamevar="star", in_shifts_file=["s1.txt", "s2.txt"],
                 append_refnum_to_fieldlabel="refn", maskpoints="m",
                 noshiftmasked=True, strip_fitsheader="SHFT",
                 strip_stitchparams=True, strip_hdu="extension")
    assert u._to_cli_args() == [
        "-unstitch", "mag,mag2", "in_shifts_file", "field", "star",
        "s1.txt,s2.txt", "append_refnum_to_fieldlabel", "refn",
        "maskpoints", "m", "noshiftmasked",
        "strip_fitsheader", "SHFT", "stitchparams", "extension"]


def test_tokens_fitsheader():
    u = unstitch("mag", "fitsheader", keywordbase="SHFT", lcnum_var="lcnum",
                 refnum_var="refn", hdu="extension")
    assert u._to_cli_args() == [
        "-unstitch", "mag", "fitsheader", "SHFT", "lcnum",
        "refnum_var", "refn", "extension"]


def test_tokens_fitsheader_with_strip():
    u = unstitch("mag", "fitsheader", keywordbase="SHFT", lcnum_var="lcnum",
                 strip_fitsheader="SHFT")
    assert u._to_cli_args() == [
        "-unstitch", "mag", "fitsheader", "SHFT", "lcnum",
        "strip_fitsheader", "SHFT"]


def test_tokens_strip_primary_default_no_hdu_token():
    # strip without an explicit hdu emits no primary/extension token
    u = unstitch("mag", "fitsheader", keywordbase="SHFT", lcnum_var="lcnum",
                 strip_fitsheader="SHFT", strip_stitchparams=True)
    args = u._to_cli_args()
    assert args[-3:] == ["strip_fitsheader", "SHFT", "stitchparams"]


def test_tokens_libprefix():
    u = unstitch("mag", "fitsheader", keywordbase="SHFT", lcnum_var="lcnum",
                 lib_path="/path/unstitch.so")
    assert u._to_cli_args()[:2] == ["-L", "/path/unstitch.so"]


# ---------------------------------------------------------------------------
# Constructor-time validation
# ---------------------------------------------------------------------------

def test_validate_bad_source():
    with pytest.raises(ValueError):
        unstitch("mag", "bogus")


@pytest.mark.parametrize("kw", [
    dict(fieldlabelsvar="f", starnamevar="s"),       # missing in_shifts_file
    dict(fieldlabelsvar="f", in_shifts_file="x"),    # missing starnamevar
    dict(starnamevar="s", in_shifts_file="x"),       # missing fieldlabelsvar
])
def test_validate_in_shifts_file_required(kw):
    with pytest.raises(ValueError):
        unstitch("mag", "in_shifts_file", **kw)


@pytest.mark.parametrize("kw", [
    dict(keywordbase="SHFT"),     # missing lcnum_var
    dict(lcnum_var="lcnum"),      # missing keywordbase
])
def test_validate_fitsheader_required(kw):
    with pytest.raises(ValueError):
        unstitch("mag", "fitsheader", **kw)


def test_validate_noshiftmasked_requires_maskpoints():
    with pytest.raises(ValueError):
        unstitch("mag", "fitsheader", keywordbase="SHFT", lcnum_var="lcnum",
                 noshiftmasked=True)


@pytest.mark.parametrize("field", ["hdu", "strip_hdu"])
def test_validate_bad_hdu(field):
    with pytest.raises(ValueError):
        unstitch("mag", "fitsheader", keywordbase="SHFT", lcnum_var="lcnum",
                 **{field: "middle"})


# ---------------------------------------------------------------------------
# Importability / Pipeline attachment
# ---------------------------------------------------------------------------

def test_attached_to_pipeline_and_toplevel():
    assert hasattr(vt.Pipeline, "unstitch")
    assert hasattr(vt, "unstitch")
    assert unstitch._vt_name == "unstitch"


# ---------------------------------------------------------------------------
# End-to-end against the real binary
# ---------------------------------------------------------------------------

def _examples_present():
    return (os.path.isfile(os.path.join(EXAMPLES_DIR, "2"))
            and os.path.isfile(os.path.join(EXAMPLES_DIR, "2.shifted")))


def _vkey(result, needle):
    """Return the value of the single result column whose name contains
    *needle*."""
    keys = [k for k in result.vars.index if needle in k]
    assert len(keys) == 1, (needle, list(result.vars.index))
    return result.vars[keys[0]]


@needs_binary
def test_end_to_end_in_shifts_file_roundtrip(tmp_path):
    """stitch (write the shifts) then stitch + unstitch reading them back:
    the magnitudes are restored (residual at machine precision) and every
    point is un-shifted."""
    if not _examples_present():
        pytest.skip("Example LCs not present")
    lc1 = os.path.join(EXAMPLES_DIR, "2")
    lc2 = os.path.join(EXAMPLES_DIR, "2.shifted")
    shifts_path = str(tmp_path / "shifts.txt")
    seg = {"field": ["fA", "fB"]}
    star = {"star": "s1"}

    # 1. determine and save the shifts
    r1 = (vt.Pipeline()
          .expr("mask=1")
          .stitch("mag", "err", "mask", "lcnum", method="median",
                  shifts_file=("field", "star"), out_shifts_file=shifts_path)
          ).run_combinelc([lc1, lc2], perlcsegment_vars=seg, perlc_vars=star)
    assert r1.error is None
    assert os.path.isfile(shifts_path)

    # 2. re-stitch and immediately un-stitch using the saved shifts
    r2 = (vt.Pipeline()
          .expr("mask=1")
          .expr("magorig=mag")
          .stitch("mag", "err", "mask", "lcnum", method="median")
          .unstitch("mag", source="in_shifts_file", fieldlabelsvar="field",
                    starnamevar="star", in_shifts_file=shifts_path)
          .expr("absdiff=abs(mag-magorig)")
          .stats("absdiff", "max")
          ).run_combinelc([lc1, lc2], perlcsegment_vars=seg, perlc_vars=star)
    assert r2.error is None
    assert int(_vkey(r2, "Unstitch_Npoints_shifted")) == 6626
    assert float(_vkey(r2, "STATS_absdiff_MAX")) < 1e-10


@needs_binary
def test_end_to_end_noshiftmasked_skips_masked(tmp_path):
    """With ``noshiftmasked`` the masked points are not un-shifted, so fewer
    points receive a shift than without it."""
    if not _examples_present():
        pytest.skip("Example LCs not present")
    lc1 = os.path.join(EXAMPLES_DIR, "2")
    lc2 = os.path.join(EXAMPLES_DIR, "2.shifted")
    shifts_path = str(tmp_path / "shifts.txt")
    seg = {"field": ["fA", "fB"]}
    star = {"star": "s1"}

    (vt.Pipeline()
     .expr("mask=1")
     .stitch("mag", "err", "mask", "lcnum", method="median",
             shifts_file=("field", "star"), out_shifts_file=shifts_path)
     ).run_combinelc([lc1, lc2], perlcsegment_vars=seg, perlc_vars=star)

    def n_shifted(noshiftmasked):
        r = (vt.Pipeline()
             .expr("mask=(t<53726.0)")
             .unstitch("mag", source="in_shifts_file", fieldlabelsvar="field",
                       starnamevar="star", in_shifts_file=shifts_path,
                       maskpoints="mask", noshiftmasked=noshiftmasked)
             ).run_combinelc([lc1, lc2], perlcsegment_vars=seg, perlc_vars=star)
        assert r.error is None
        return int(_vkey(r, "Unstitch_Npoints_shifted"))

    n_default = n_shifted(False)
    n_noshift = n_shifted(True)
    assert n_default == 6626                 # every matched point shifted
    assert 0 < n_noshift < n_default         # masked points skipped


@needs_binary
def test_end_to_end_coverage_error(tmp_path):
    """Un-stitching with a shifts file that does not cover the star fails."""
    if not _examples_present():
        pytest.skip("Example LCs not present")
    lc1 = os.path.join(EXAMPLES_DIR, "2")
    lc2 = os.path.join(EXAMPLES_DIR, "2.shifted")
    shifts_path = str(tmp_path / "shifts.txt")
    # a shifts file for a different star name
    with open(shifts_path, "w") as f:
        f.write("OTHERSTAR fA,0,10;fB,0.3,10\n")

    from pyvartools.results import RunError
    with pytest.raises(RunError):
        (vt.Pipeline()
         .expr("mask=1")
         .unstitch("mag", source="in_shifts_file", fieldlabelsvar="field",
                   starnamevar="star", in_shifts_file=shifts_path)
         ).run_combinelc([lc1, lc2],
                         perlcsegment_vars={"field": ["fA", "fB"]},
                         perlc_vars={"star": "s1"})
