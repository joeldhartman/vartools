"""Tests for ``segment_vars`` / ``lc_vars`` on ``run_combinelc(s)``.

These two kwargs let users attach Python data to each segment / group of
a combinelcs run without hand-rolling a vartools list file.  The tests
here cover:

* type inference and explicit ``(values, type)`` overrides;
* the rendered list-file content (numerically formatted floats, integer
  rendering, string rendering, comma-joined per-segment subcolumns);
* the inlistvars wiring (segment_vars get ``combinelc=True``, lc_vars do
  not);
* validation: shape mismatches, name collisions, embedded whitespace in
  string values;
* the singular-form ``run_combinelc`` auto-wrap;
* end-to-end against the real binary, exercising the stitch
  ``shifts_file`` path that motivated the feature.
"""
from __future__ import annotations

import os
import pathlib
import sys
import tempfile

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import pyvartools as vt
from pyvartools.pipeline import Pipeline, ListVar

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
# Helpers operating on the static Pipeline normalisation methods
# ---------------------------------------------------------------------------

class TestNormalize:
    def test_infer_type_string(self):
        vs, t = Pipeline._normalize_extravar_spec(["A", "B", "C"])
        assert t == "string"
        assert vs == ["A", "B", "C"]

    def test_infer_type_int(self):
        vs, t = Pipeline._normalize_extravar_spec([1, 2, 3])
        assert t == "int"

    def test_infer_type_double(self):
        vs, t = Pipeline._normalize_extravar_spec([1.0, 2.5, 3.0])
        assert t == "double"

    def test_infer_type_bool_maps_to_int(self):
        # vartools has no boolean type
        vs, t = Pipeline._normalize_extravar_spec([True, False])
        assert t == "int"

    def test_infer_skips_none(self):
        vs, t = Pipeline._normalize_extravar_spec([None, "X"])
        assert t == "string"

    def test_explicit_tuple_override(self):
        vs, t = Pipeline._normalize_extravar_spec((["1", "2"], "string"))
        assert t == "string"
        assert vs == ["1", "2"]

    def test_explicit_tuple_invalid_type_treated_as_values(self):
        # If the second element is not a recognised type tag, the spec is
        # interpreted as a bare 2-tuple of values, not a (values, type).
        vs, t = Pipeline._normalize_extravar_spec(("A", "B"))
        assert t == "string"
        assert vs == ["A", "B"]

    def test_infer_type_nested(self):
        # segment_vars values are nested list-of-lists
        vs, t = Pipeline._normalize_extravar_spec([["A", "B"], ["C", "D"]])
        assert t == "string"

    def test_infer_type_unsupported_raises(self):
        with pytest.raises(TypeError):
            Pipeline._normalize_extravar_spec([object()])

    def test_infer_type_all_none_raises(self):
        with pytest.raises(ValueError):
            Pipeline._normalize_extravar_spec([None, None])


class TestFormat:
    def test_string_value(self):
        assert Pipeline._format_extravar_value("ABC", "string") == "ABC"

    def test_string_with_whitespace_rejected(self):
        with pytest.raises(ValueError, match="whitespace"):
            Pipeline._format_extravar_value("a b", "string")

    def test_int_value(self):
        assert Pipeline._format_extravar_value(7, "int") == "7"
        assert Pipeline._format_extravar_value(7.9, "int") == "7"  # truncation

    def test_double_value(self):
        out = Pipeline._format_extravar_value(1.234567890123, "double")
        # %.10g formatting (10 significant digits, trailing zeros stripped)
        assert out == "1.23456789"

    def test_none_rejected(self):
        with pytest.raises(ValueError):
            Pipeline._format_extravar_value(None, "string")


# ---------------------------------------------------------------------------
# List-file content + inlistvars wiring (no binary needed)
# ---------------------------------------------------------------------------

class _ListFileSpy:
    """Run a Pipeline up to where the list file is written, then peek at it.

    We replace ``self._execute`` on a real Pipeline with a stub that throws
    a sentinel exception once the list file (and the assembled command) has
    been built.  The temporary directory and its contents are still alive
    inside the ``with tempfile.TemporaryDirectory()`` block at that point,
    so the spy reads them before re-raising.
    """

    def __init__(self):
        self.list_file_contents = None
        self.cmd = None

    def install(self, pipe: Pipeline):
        outer = self

        def fake_execute(self_pipe, cmd, timeout=None, stdin_text=None):
            list_idx = cmd.index("-l") + 1
            list_path = cmd[list_idx]
            with open(list_path) as f:
                outer.list_file_contents = f.read()
            outer.cmd = list(cmd)
            raise _Stop()

        pipe._execute = fake_execute.__get__(pipe, Pipeline)


class _Stop(Exception):
    pass


def _spy_combinelcs(pipe, **kwargs):
    spy = _ListFileSpy()
    spy.install(pipe)
    with pytest.raises(_Stop):
        pipe.run_combinelcs(**kwargs)
    return spy


class TestListFileWriting:
    def test_segment_vars_string(self):
        pipe = vt.Pipeline().rms()
        spy = _spy_combinelcs(pipe,
            groups=[["a.txt", "b.txt"], ["c.txt", "d.txt"]],
            segment_vars={"fld": [["A", "B"], ["C", "D"]]},
        )
        # Two rows, each: "<paths>  <comma-joined fields>"
        lines = spy.list_file_contents.strip().splitlines()
        assert len(lines) == 2
        assert lines[0] == "a.txt,b.txt A,B"
        assert lines[1] == "c.txt,d.txt C,D"

    def test_lc_vars_scalar(self):
        pipe = vt.Pipeline().rms()
        spy = _spy_combinelcs(pipe,
            groups=[["a.txt"], ["b.txt"]],
            lc_vars={"star": ["TIC1", "TIC2"]},
        )
        lines = spy.list_file_contents.strip().splitlines()
        assert lines == ["a.txt TIC1", "b.txt TIC2"]

    def test_segment_and_lc_vars_together(self):
        pipe = vt.Pipeline().rms()
        spy = _spy_combinelcs(pipe,
            groups=[["a.txt", "b.txt"]],
            segment_vars={"fld": [["A", "B"]]},
            lc_vars={"star": ["TIC9"]},
        )
        # Cols: 1=paths(comma), 2=fld(comma per seg), 3=star
        line = spy.list_file_contents.strip()
        assert line == "a.txt,b.txt A,B TIC9"

    def test_int_and_double_rendering(self):
        pipe = vt.Pipeline().rms()
        spy = _spy_combinelcs(pipe,
            groups=[["a.txt", "b.txt"]],
            segment_vars={"telnum": [[1, 2]]},
            lc_vars={"period": [1.234567890123]},
        )
        line = spy.list_file_contents.strip()
        # int rendered as bare digits, double as %.10g
        assert line == "a.txt,b.txt 1,2 1.23456789"

    def test_inlistvars_wiring(self):
        pipe = vt.Pipeline().rms()
        spy = _spy_combinelcs(pipe,
            groups=[["a.txt", "b.txt"]],
            segment_vars={"fld": [["A", "B"]]},
            lc_vars={"star": ["TIC9"]},
        )
        # Find the -inlistvars argument in the assembled cmd.
        i = spy.cmd.index("-inlistvars")
        ilv = spy.cmd[i + 1]
        # segment_vars get combinelc; lc_vars do not.
        # The comma-separated entry list looks like
        #   fld:2:combinelc:string,star:3:string
        entries = ilv.split(",")
        d = {}
        for ent in entries:
            d[ent.split(":")[0]] = ent
        assert d["fld"] == "fld:2:combinelc:string"
        assert d["star"] == "star:3:string"

    def test_custom_delimiter(self):
        pipe = vt.Pipeline().rms()
        spy = _spy_combinelcs(pipe,
            groups=[["a.txt", "b.txt"]],
            segment_vars={"fld": [["A", "B"]]},
            delimiter=";",
        )
        line = spy.list_file_contents.strip()
        # Both the path subcolumns and the segment_vars subcolumns use
        # the same delimiter as combinelcs itself.
        assert line == "a.txt;b.txt A;B"


# ---------------------------------------------------------------------------
# Validation errors
# ---------------------------------------------------------------------------

class TestValidation:
    def test_segment_vars_wrong_outer_length(self):
        pipe = vt.Pipeline().rms()
        with pytest.raises(ValueError, match="segment_vars"):
            pipe.run_combinelcs(
                groups=[["a.txt"], ["b.txt"]],
                segment_vars={"fld": [["A"]]},  # only 1 entry, need 2
            )

    def test_segment_vars_wrong_inner_length(self):
        pipe = vt.Pipeline().rms()
        with pytest.raises(ValueError, match="segment_vars"):
            pipe.run_combinelcs(
                groups=[["a.txt", "b.txt"]],
                segment_vars={"fld": [["A"]]},  # 1 value, need 2
            )

    def test_segment_vars_inner_not_a_list(self):
        pipe = vt.Pipeline().rms()
        with pytest.raises(TypeError, match="segment_vars"):
            pipe.run_combinelcs(
                groups=[["a.txt", "b.txt"]],
                segment_vars={"fld": ["A"]},  # bare string, not list-of-list
            )

    def test_lc_vars_wrong_length(self):
        pipe = vt.Pipeline().rms()
        with pytest.raises(ValueError, match="lc_vars"):
            pipe.run_combinelcs(
                groups=[["a.txt"], ["b.txt"]],
                lc_vars={"star": ["TIC1"]},
            )

    def test_name_collision_with_inlistvars(self):
        pipe = vt.Pipeline().rms()
        with pytest.raises(ValueError, match="both inlistvars"):
            pipe.run_combinelcs(
                groups=[["a.txt"]],
                inlistvars={"x": 7},
                lc_vars={"x": ["foo"]},
            )

    def test_name_collision_segment_vs_lc(self):
        pipe = vt.Pipeline().rms()
        with pytest.raises(ValueError, match="both segment_vars and lc_vars"):
            pipe.run_combinelcs(
                groups=[["a.txt", "b.txt"]],
                segment_vars={"x": [["A", "B"]]},
                lc_vars={"x": ["foo"]},
            )

    def test_string_with_whitespace_rejected(self):
        pipe = vt.Pipeline().rms()
        with pytest.raises(ValueError, match="whitespace"):
            pipe.run_combinelcs(
                groups=[["a.txt"]],
                lc_vars={"star": ["a b"]},
            )


# ---------------------------------------------------------------------------
# Singular run_combinelc auto-wrap
# ---------------------------------------------------------------------------

class TestSingularAutoWrap:
    def test_singular_segment_vars(self):
        """``run_combinelc(files, segment_vars={k: [...]})`` accepts a flat
        per-segment list and wraps it as a single-group entry for the
        plural form."""
        pipe = vt.Pipeline().rms()
        spy = _ListFileSpy()
        spy.install(pipe)
        with pytest.raises(_Stop):
            pipe.run_combinelc(
                ["a.txt", "b.txt"],
                segment_vars={"fld": ["A", "B"]},
                lc_vars={"star": "TIC9"},
            )
        line = spy.list_file_contents.strip()
        assert line == "a.txt,b.txt A,B TIC9"

    def test_singular_explicit_type(self):
        pipe = vt.Pipeline().rms()
        spy = _ListFileSpy()
        spy.install(pipe)
        with pytest.raises(_Stop):
            pipe.run_combinelc(
                ["a.txt", "b.txt"],
                segment_vars={"id": (["001", "002"], "string")},
                lc_vars={"period": (1.5, "double")},
            )
        line = spy.list_file_contents.strip()
        assert line == "a.txt,b.txt 001,002 1.5"

    def test_singular_empty_files_raises(self):
        with pytest.raises(ValueError):
            vt.Pipeline().rms().run_combinelc([])


# ---------------------------------------------------------------------------
# End-to-end against the real binary
# ---------------------------------------------------------------------------

@needs_binary
def test_end_to_end_stitch_shifts_file(tmp_path):
    """The motivating use case: per-segment field labels and per-LC
    starname feed ``-stitch shifts_file``, and the recovered shifts are
    written keyed by field label."""
    lc1 = os.path.join(EXAMPLES_DIR, "2")
    lc2 = os.path.join(EXAMPLES_DIR, "2.shifted")
    if not (os.path.isfile(lc1) and os.path.isfile(lc2)):
        pytest.skip("Example LCs not present")

    shifts_path = str(tmp_path / "shifts.txt")

    pipe = (vt.Pipeline()
            .expr("mask=1")
            .stitch("mag", "err", "mask", "lcnum",
                    method="poly 5", groupbytime=0.5,
                    groupbytime_start=54550.123, fitonly=True,
                    shifts_file=("fieldname", "starname"),
                    out_shifts_file=shifts_path))

    result = pipe.run_combinelc(
        [lc1, lc2],
        segment_vars={"fieldname": ["2_A", "2_B"]},
        lc_vars={"starname": "2"},
    )
    assert result.error is None

    with open(shifts_path) as f:
        contents = f.read().strip()
    # Format:  <starname> <field>,<shift>,<nobs>;<field>,<shift>,<nobs>
    assert contents.startswith("2 ")
    assert "2_A" in contents
    assert "2_B" in contents
    # The shifted segment was offset by 0.3 mag — the recovered shift
    # should be close to ±0.3 (sign convention is implementation detail).
    parts = contents.split()[1].split(";")
    seg_b = [p for p in parts if p.startswith("2_B")][0].split(",")
    assert abs(abs(float(seg_b[1])) - 0.3) < 0.05


@needs_binary
def test_end_to_end_plural_two_groups(tmp_path):
    """``run_combinelcs`` with two groups and per-group string vars."""
    lc1 = os.path.join(EXAMPLES_DIR, "2")
    lc2 = os.path.join(EXAMPLES_DIR, "2.shifted")
    if not (os.path.isfile(lc1) and os.path.isfile(lc2)):
        pytest.skip("Example LCs not present")

    shifts_path1 = str(tmp_path / "g1.txt")
    shifts_path2 = str(tmp_path / "g2.txt")
    # vartools writes one out_shifts_file per LC keyed by lc index in
    # the header; we just need a single shifts_file path here that the
    # stitch wrapper will reuse.  The actual content lands wherever
    # vartools writes per-LC.  For this regression it's sufficient that
    # both groups run without error and produce two stats rows.
    pipe = (vt.Pipeline()
            .expr("mask=1")
            .stitch("mag", "err", "mask", "lcnum",
                    method="poly 3", fitonly=True))
    batch = pipe.run_combinelcs(
        groups=[[lc1, lc2], [lc1, lc2]],
        segment_vars={"fld": [["A1", "B1"], ["A2", "B2"]]},
        lc_vars={"star": ["S1", "S2"]},
    )
    assert batch.error is None
    assert len(batch.vars) == 2
