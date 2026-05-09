"""Tests for the values-form of ``perlc_vars`` in run_batch and
LightCurveBatch.run, and the rejection path in run_filelist.

The schema-form of perlc_vars (int / PerLCColumn) is exercised across
test_infrastructure.py and test_all_commands.py via the renamed
inlistvars-style usage; this module covers the new sequence-form
(perlc_vars={"name": ["a", "b", "c"]}) added in commit 2 of the rename.
"""

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

VARTOOLS_SRC = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", ".."))
EXAMPLES_DIR = os.path.join(VARTOOLS_SRC, "EXAMPLES")


def _make_lcs(n=3, length=100):
    """Build a list of *n* synthetic light curves of length *length*."""
    out = []
    for i in range(n):
        t = np.linspace(0, 20, length)
        mag = 10.0 + 0.1 * np.sin(2 * np.pi * t / (1.5 + 0.3 * i))
        err = np.full(length, 0.01)
        lc = vt.LightCurve.from_arrays(t, mag, err)
        lc.name = f"lc{i}"
        out.append(lc)
    return out


# ---------------------------------------------------------------------------
# Schema-form perlc_vars on run_batch is unchanged (smoke test only — schema
# usage is broadly covered elsewhere)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools binary required")
def test_perlc_vars_schema_form_unchanged():
    """A col=0 schema entry remains supported in run_batch."""
    lcs = _make_lcs(2)
    pipe = vt.Pipeline().rms()
    batch = pipe.run_batch(
        lcs,
        perlc_vars={"row": vt.PerLCColumn(col=0, type="int", init="NF")},
    )
    assert len(batch.vars) == 2


# ---------------------------------------------------------------------------
# Values-form on run_batch
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools binary required")
def test_perlc_vars_values_string_namefromlist(tmp_path):
    """The original user repro: per-LC output names via cmd.o(namefromlist=)."""
    lcs = _make_lcs(3)
    outnames = [f"lc{i}_renamed" for i in range(3)]
    outdir = str(tmp_path / "out")
    os.makedirs(outdir, exist_ok=True)

    pipe = (vt.Pipeline()
            .clip(5.0)
            .o(outdir=outdir, allcols=True, namefromlist="outname"))
    batch = pipe.run_batch(lcs, perlc_vars={"outname": outnames})

    assert len(batch.vars) == 3
    # vartools writes one file per LC named after the perlc_vars value.
    written = sorted(os.listdir(outdir))
    expected = sorted(outnames)
    assert written == expected, (
        f"Expected files {expected} in outdir, got {written}"
    )


@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools binary required")
def test_perlc_vars_values_numeric_drives_LS_bounds():
    """Numeric per-LC values steer cmd.LS minp/maxp through -inlistvars."""
    lcs = _make_lcs(2)
    # Distinct minp / maxp per LC; the var-name forms reach into perlc_vars.
    pipe = vt.Pipeline().LS("minp", "maxp", 0.1)
    batch = pipe.run_batch(
        lcs,
        perlc_vars={"minp": [0.3, 0.6], "maxp": [3.0, 5.0]},
    )
    assert len(batch.vars) == 2
    # Sanity: LS produced a Period_1 column for both LCs.
    assert "LS_Period_1_0" in batch.vars.columns
    assert batch.vars["LS_Period_1_0"].notna().all()


@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools binary required")
def test_perlc_vars_values_with_explicit_type():
    """``(values, type)`` tuple form selects the vartools type tag."""
    lcs = _make_lcs(2)
    pipe = vt.Pipeline().rms()
    batch = pipe.run_batch(
        lcs,
        perlc_vars={"label": (["A", "B"], "string")},
    )
    assert len(batch.vars) == 2


@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools binary required")
def test_perlc_vars_mixed_with_perlc_cmd_attr():
    """Values-form perlc_vars co-exists with PerLC(...) on a command attr."""
    lcs = _make_lcs(2)
    pipe = vt.Pipeline().LS(
        vt.PerLC([0.3, 0.6]),  # per-LC minp via cmd-attr PerLC
        2.0,                    # fixed maxp
        0.1,
    ).rms()
    batch = pipe.run_batch(
        lcs,
        perlc_vars={"label": (["X", "Y"], "string")},
    )
    assert len(batch.vars) == 2


# ---------------------------------------------------------------------------
# Validation paths (no binary needed)
# ---------------------------------------------------------------------------


def test_perlc_vars_values_length_mismatch_run_batch():
    """A values list of the wrong length is rejected with a clear message."""
    lcs = _make_lcs(3)
    pipe = vt.Pipeline().rms()
    with pytest.raises(ValueError, match=r"perlc_vars\['x'\] has 2 values"):
        pipe.run_batch(lcs, perlc_vars={"x": [1.0, 2.0]})


def test_perlc_vars_values_rejected_in_run_filelist(tmp_path):
    """run_filelist refuses sequence-form perlc_vars and points at run_batch."""
    list_path = tmp_path / "list.txt"
    list_path.write_text(f"{os.path.join(EXAMPLES_DIR, '2')}\n")
    pipe = vt.Pipeline().rms()
    with pytest.raises(ValueError, match="run_batch"):
        pipe.run_filelist(
            str(list_path),
            perlc_vars={"x": [1.0]},
        )


def test_perlc_vars_schema_still_works_in_run_filelist(tmp_path):
    """Schema-form (int / PerLCColumn) entries still work in run_filelist."""
    # Build a list file with a numeric extra column the user wants to
    # reference via inlistvars (schema form).
    paths = [os.path.join(EXAMPLES_DIR, "2")]
    list_path = tmp_path / "list.txt"
    with open(list_path, "w") as f:
        for p in paths:
            f.write(f"{p} 0.5\n")
    pipe = vt.Pipeline().rms()
    if not _HAVE_BINARY:
        pytest.skip("vartools binary required")
    # int shorthand for column 2.
    batch = pipe.run_filelist(
        str(list_path),
        perlc_vars={"thr": 2},
    )
    assert len(batch.vars) == 1


# ---------------------------------------------------------------------------
# LightCurveBatch.run() mirror
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools binary required")
def test_lightcurvebatch_run_perlc_vars_routes_to_run_batch(tmp_path):
    """LightCurveBatch.run() with perlc_vars takes the run_batch path
    (single vartools invocation)."""
    lcs = _make_lcs(3)
    outnames = [f"lc{i}_via_lcb" for i in range(3)]
    outdir = str(tmp_path / "out_lcb")
    os.makedirs(outdir, exist_ok=True)

    batch = vt.LightCurveBatch(lcs).clip(5.0).o(
        outdir=outdir, allcols=True, namefromlist="outname"
    ).run(perlc_vars={"outname": outnames})

    assert len(batch.vars) == 3
    written = sorted(os.listdir(outdir))
    assert written == sorted(outnames)
