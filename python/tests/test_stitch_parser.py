"""Tests for the stitch-fitted-parameters file parser.

The ``-stitch`` user command writes a structured-text file when
``save_fitted_parameters`` is set; the format is not whitespace-tabular,
so pyvartools ships a dedicated parser
(:func:`pyvartools.commands.userlibs._parse_stitch_fitted_params`).
This file pins the parser against synthetic inputs that mirror every
``fprintf`` line currently emitted from ``USERLIBS/src/stitch.c``,
plus one end-to-end run of the real binary when available.

If the file format ever changes on either side (C writer or Python
parser) without a corresponding update to the other, this test will
fail loudly.
"""
from __future__ import annotations

import math
import os
import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import pyvartools as vt
from pyvartools.commands.userlibs import _parse_stitch_fitted_params

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
# Unit tests against synthetic input
# ---------------------------------------------------------------------------

def _write(tmp_path, text):
    p = tmp_path / "fitted.stitch"
    p.write_text(text)
    return str(p)


def test_poly_method(tmp_path):
    """Polynomial coefficients across two time bins, plus shifts."""
    text = (
        "# Parameters for stitch variable 1\n"
        "Coeff for t^0, 100.0<t<101.0: 10.0\n"
        "Coeff for t^1, 100.0<t<101.0: -0.5\n"
        "Coeff for t^2, 100.0<t<101.0: 1.5\n"
        "Coeff for t^0, 101.0<t<102.0: 11.0\n"
        "Coeff for t^1, 101.0<t<102.0: 0.25\n"
        "Coeff for t^2, 101.0<t<102.0: -2.0\n"
        "LCgroup_1 shift: 0.1\n"
        "LCgroup_2 shift: -0.3\n"
    )
    df = _parse_stitch_fitted_params(_write(tmp_path, text))

    assert list(df.columns) == ["variable", "kind", "term",
                                "t_min", "t_max", "value"]
    assert (df.variable == 1).all()

    coeffs = df[df.kind == "coeff"]
    assert len(coeffs) == 6
    assert list(coeffs.term) == ["t^0", "t^1", "t^2"] * 2
    assert list(coeffs.value) == [10.0, -0.5, 1.5, 11.0, 0.25, -2.0]
    assert list(coeffs.t_min) == [100.0] * 3 + [101.0] * 3
    assert list(coeffs.t_max) == [101.0] * 3 + [102.0] * 3

    shifts = df[df.kind == "shift"]
    assert len(shifts) == 2
    assert list(shifts.term) == ["LCgroup_1", "LCgroup_2"]
    assert list(shifts.value) == [0.1, -0.3]
    assert shifts.t_min.isna().all()
    assert shifts.t_max.isna().all()


@pytest.mark.parametrize("statname", ["Median", "Mean", "WeightedMean"])
def test_stat_methods(tmp_path, statname):
    """Median / Mean / WeightedMean each emit ``STATNAME TMIN<t<TMAX: VALUE``."""
    text = (
        "# Parameters for stitch variable 1\n"
        f"{statname} 50.0<t<60.0: 12.34\n"
        f"{statname} 60.0<t<70.0: 12.56\n"
        "LCgroup_1 shift: 0.0\n"
        "LCgroup_2 shift: -0.42\n"
    )
    df = _parse_stitch_fitted_params(_write(tmp_path, text))
    coeffs = df[df.kind == "coeff"]
    assert list(coeffs.term) == [statname, statname]
    assert list(coeffs.value) == [12.34, 12.56]
    assert list(coeffs.t_min) == [50.0, 60.0]
    assert list(coeffs.t_max) == [60.0, 70.0]


def test_harmseries_method(tmp_path):
    """Harmonic series: DC term plus cos/sin pairs at K=1, K=2."""
    text = (
        "# Parameters for stitch variable 1\n"
        "Coeff for 1, 0.0<t<1.0: 10.0\n"
        "Coeff for cos(1*t*f0), 0.0<t<1.0: 0.5\n"
        "Coeff for sin(1*t*f0), 0.0<t<1.0: -0.25\n"
        "Coeff for cos(2*t*f0), 0.0<t<1.0: 0.1\n"
        "Coeff for sin(2*t*f0), 0.0<t<1.0: -0.05\n"
        "LCgroup_1 shift: 0.0\n"
    )
    df = _parse_stitch_fitted_params(_write(tmp_path, text))
    coeffs = df[df.kind == "coeff"]
    assert list(coeffs.term) == [
        "1",
        "cos(1*t*f0)", "sin(1*t*f0)",
        "cos(2*t*f0)", "sin(2*t*f0)",
    ]
    assert list(coeffs.value) == [10.0, 0.5, -0.25, 0.1, -0.05]
    assert (coeffs.t_min == 0.0).all()
    assert (coeffs.t_max == 1.0).all()


def test_multiple_stitch_variables(tmp_path):
    """A stitch with two variables emits two ``# Parameters`` blocks separated
    by a blank line; the parser tags each row with its source variable."""
    text = (
        "# Parameters for stitch variable 1\n"
        "Coeff for t^0, 0.0<t<1.0: 1.0\n"
        "LCgroup_1 shift: 0.0\n"
        "LCgroup_2 shift: 0.1\n"
        "\n"
        "# Parameters for stitch variable 2\n"
        "Coeff for t^0, 0.0<t<1.0: 2.0\n"
        "LCgroup_1 shift: 0.0\n"
        "LCgroup_2 shift: 0.2\n"
    )
    df = _parse_stitch_fitted_params(_write(tmp_path, text))
    assert sorted(df.variable.unique().tolist()) == [1, 2]
    v1 = df[df.variable == 1]
    v2 = df[df.variable == 2]
    assert v1[v1.kind == "coeff"].iloc[0]["value"] == 1.0
    assert v2[v2.kind == "coeff"].iloc[0]["value"] == 2.0
    assert v1[v1.term == "LCgroup_2"]["value"].iloc[0] == 0.1
    assert v2[v2.term == "LCgroup_2"]["value"].iloc[0] == 0.2


def test_high_precision_values_round_trip(tmp_path):
    """Values are written by ``%.17g`` — parser must not round."""
    val = 1.2345678901234567e-08
    text = (
        "# Parameters for stitch variable 1\n"
        f"Coeff for t^0, 0.0<t<1.0: {val:.17g}\n"
        f"LCgroup_2 shift: {-val:.17g}\n"
    )
    df = _parse_stitch_fitted_params(_write(tmp_path, text))
    assert df[df.kind == "coeff"].value.iloc[0] == val
    assert df[df.kind == "shift"].value.iloc[0] == -val


def test_empty_file_returns_empty_frame(tmp_path):
    df = _parse_stitch_fitted_params(_write(tmp_path, ""))
    assert list(df.columns) == ["variable", "kind", "term",
                                "t_min", "t_max", "value"]
    assert len(df) == 0


# ---------------------------------------------------------------------------
# Integration test against the real binary
# ---------------------------------------------------------------------------

@needs_binary
def test_end_to_end_capture(tmp_path):
    """Run the binary, capture the file via ``save_fitted_parameters=True``,
    and check the structure matches what the parser produces."""
    lc1 = os.path.join(EXAMPLES_DIR, "2")
    lc2 = os.path.join(EXAMPLES_DIR, "2.shifted")
    if not (os.path.isfile(lc1) and os.path.isfile(lc2)):
        pytest.skip("Example LCs not present")

    pipe = (vt.Pipeline()
            .expr("mask=1")
            .stitch("mag", "err", "mask", "lcnum",
                    method="poly 3", groupbytime=0.5,
                    save_fitted_parameters=True))
    result = pipe.run_combinelcs([[lc1, lc2]])

    keys = [k for k in result.files if "fitted_parameters" in k]
    assert len(keys) == 1, result.files
    df = result.files[keys[0]][0]

    assert list(df.columns) == ["variable", "kind", "term",
                                "t_min", "t_max", "value"]
    assert {"coeff", "shift"} <= set(df.kind.unique())

    # poly 3 has 4 coefficients (t^0..t^3) per time bin
    coeff_terms = sorted(df[df.kind == "coeff"].term.unique())
    assert coeff_terms == ["t^0", "t^1", "t^2", "t^3"]

    # Two LCs combined → there should be exactly one LCgroup_2 shift row
    shifts = df[df.kind == "shift"]
    assert (shifts.term == "LCgroup_2").any()
    # The shifted LC was created by adding 0.3 mag, so the recovered
    # shift for LCgroup_2 should be approximately ±0.3.
    lcg2 = shifts[shifts.term == "LCgroup_2"].value.iloc[0]
    assert math.isclose(abs(lcg2), 0.3, abs_tol=0.05)
