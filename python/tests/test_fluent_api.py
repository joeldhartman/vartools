"""
Tests for the pyvartools command API.

Covers:
- lc.CMD(...) immediate execution → Result
- result.CMD(...) immediate execution → Result with merged vars
- Pipeline-stateful commands raise NotImplementedError on lc and Result
- LightCurveBatch deferred chaining and run()
- LightCurveBatch.run_CMD(...) immediate execution
- BatchResult.CMD(...) continuation chaining
- BatchResult.run_CMD(...) immediate execution
- VarsNamespace structured access
- result.vars dict-style access (pd.Series)
- result.<key> __getattr__ shorthand
- Top-level vt.CMD(lc_input, ...) functions
- BatchResult slicing
"""

import os
import pytest
import numpy as np
import pandas as pd

import pyvartools as vt
from pyvartools._batch import LightCurveBatch
from pyvartools.results import BatchResult, Result, RunError

VARTOOLS_SRC = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", ".."))
EXAMPLES_DIR = os.path.join(VARTOOLS_SRC, "EXAMPLES")


@pytest.fixture(scope="module")
def lc():
    return vt.LightCurve.from_file(os.path.join(EXAMPLES_DIR, "2"))


@pytest.fixture(scope="module")
def three_lcs():
    return [vt.LightCurve.from_file(os.path.join(EXAMPLES_DIR, str(i)))
            for i in range(1, 4)]


# ---------------------------------------------------------------------------
# lc.CMD() — immediate execution
# ---------------------------------------------------------------------------

class TestLcImmediate:

    def test_lc_cmd_returns_result(self, lc):
        result = lc.LS(0.1, 10.0, 1e-3)
        assert isinstance(result, Result)

    def test_lc_cmd_contains_output_key(self, lc):
        result = lc.LS(0.1, 10.0, 1e-3)
        assert "LS_Period_1_0" in result.vars.index

    def test_lc_cmd_result_ok(self, lc):
        result = lc.LS(0.1, 10.0, 1e-3)
        assert result.ok
        assert result.error is None

    def test_lc_cmd_capture_lc_default_true(self, lc):
        result = lc.LS(0.1, 10.0, 1e-3)
        assert result.lc is not None

    def test_lc_cmd_capture_lc_false(self, lc):
        result = lc.LS(0.1, 10.0, 1e-3, capture_lc=False)
        assert result.lc is None

    def test_lc_rms_returns_result(self, lc):
        result = lc.rms()
        assert isinstance(result, Result)
        assert "RMS_0" in result.vars.index

    def test_lc_varobjs(self, lc):
        result = lc.LS(0.1, 10.0, 1e-3)
        assert result.varobjs.LS.Period_1 > 0


# ---------------------------------------------------------------------------
# result.CMD() — immediate execution with var merging
# ---------------------------------------------------------------------------

class TestResultImmediate:

    def test_result_cmd_returns_result(self, lc):
        r1 = lc.LS(0.1, 10.0, 1e-3)
        r2 = r1.rms()
        assert isinstance(r2, Result)

    def test_result_cmd_merges_prior_vars(self, lc):
        """result.rms() should preserve LS vars from the prior result."""
        r1 = lc.LS(0.1, 10.0, 1e-3)
        r2 = r1.rms()
        # LS vars from r1 should appear in r2 (position 0)
        assert "LS_Period_1_0" in r2.vars.index
        # rms vars should appear at position 1
        assert any("RMS" in k for k in r2.vars.index)

    def test_result_cmd_varobjs_has_both_commands(self, lc):
        r1 = lc.LS(0.1, 10.0, 1e-3)
        r2 = r1.rms()
        cmds = r2.varobjs.commands()
        assert "LS" in cmds
        assert "rms" in cmds

    def test_result_cmd_matches_single_invocation(self, lc):
        """vt.clip(lc).LS() should give same varobjs as lc.clip().LS() (single invoc)."""
        # Single vartools invocation
        r_single = lc.clip(sigclip=5.0).LS(0.5, 10.0, 1e-3)
        # Chained via Result
        r_chained = vt.clip(lc, sigclip=5.0).LS(0.5, 10.0, 1e-3)
        assert sorted(r_single.vars.index.tolist()) == sorted(r_chained.vars.index.tolist())
        assert r_single.varobjs.LS.Period_1 == pytest.approx(r_chained.varobjs.LS.Period_1)

    def test_result_cmd_no_lc_raises(self, lc):
        result = lc.LS(0.1, 10.0, 1e-3, capture_lc=False)
        with pytest.raises(AttributeError, match="capture_lc"):
            result.rms()

    def test_three_step_chain(self, lc):
        """lc.clip().LS().rms() — three-step chain via Result merging."""
        r = lc.clip(sigclip=5.0).LS(0.1, 10.0, 1e-3).rms()
        assert isinstance(r, Result)
        cmds = r.varobjs.commands()
        assert "clip" in cmds
        assert "LS" in cmds
        assert "rms" in cmds


# ---------------------------------------------------------------------------
# Pipeline-stateful commands raise NotImplementedError
# ---------------------------------------------------------------------------

class TestStatefulCommandErrors:

    def test_savelc_on_lc_raises(self, lc):
        with pytest.raises(NotImplementedError, match="Pipeline"):
            lc.savelc()

    def test_restorelc_on_lc_raises(self, lc):
        with pytest.raises(NotImplementedError, match="Pipeline"):
            lc.restorelc()

    def test_columnsuffix_on_lc_raises(self, lc):
        with pytest.raises(NotImplementedError, match="Pipeline"):
            lc.columnsuffix("x")

    def test_ifcmd_on_lc_raises(self, lc):
        # ifcmd's _vt_name is "if", so it's attached as lc.if_ ... but "if"
        # is a Python keyword; the method is accessible via getattr.
        with pytest.raises(NotImplementedError, match="Pipeline"):
            getattr(lc, "if")("expr", "1")

    def test_savelc_on_result_raises(self, lc):
        r = lc.rms()
        with pytest.raises(NotImplementedError, match="Pipeline"):
            r.savelc()


# ---------------------------------------------------------------------------
# VarsNamespace and var access
# ---------------------------------------------------------------------------

class TestVarsAccess:

    def test_var_index_contains_ls_key(self, lc):
        result = lc.LS(0.1, 10.0, 1e-3)
        assert "LS_Period_1_0" in result.vars.index

    def test_getattr_shorthand(self, lc):
        result = lc.LS(0.1, 10.0, 1e-3)
        val = result.LS_Period_1_0
        assert isinstance(val, float)

    def test_getattr_missing_raises(self, lc):
        result = lc.LS(0.1, 10.0, 1e-3)
        with pytest.raises(AttributeError):
            _ = result.NonExistentKey_999

    def test_vars_namespace_ls(self, lc):
        result = lc.LS(0.1, 10.0, 1e-3)
        period = result.varobjs.LS.Period_1
        assert isinstance(period, float)
        assert period > 0

    def test_vars_namespace_single_and_indexed(self, lc):
        """Single LS call: varobjs.LS.Period_1 and varobjs.LS[0].Period_1 agree."""
        result = lc.LS(0.1, 10.0, 1e-3)
        assert result.varobjs.LS.Period_1 == result.varobjs.LS[0].Period_1


# ---------------------------------------------------------------------------
# LightCurveBatch
# ---------------------------------------------------------------------------

class TestLightCurveBatch:

    def test_batch_creation_list(self, three_lcs):
        batch = vt.LightCurveBatch(three_lcs)
        assert len(batch) == 3

    def test_batch_creation_varargs(self, three_lcs):
        batch = vt.LightCurveBatch(*three_lcs)
        assert len(batch) == 3

    def test_batch_cmd_returns_batch(self, three_lcs):
        result = vt.LightCurveBatch(three_lcs).LS(0.1, 10.0, 1e-3)
        assert isinstance(result, LightCurveBatch)

    def test_batch_chain(self, three_lcs):
        batch = vt.LightCurveBatch(three_lcs).LS(0.1, 10.0, 1e-3).rms()
        assert len(batch._commands) == 2

    def test_batch_run_returns_batch_result(self, three_lcs):
        result = vt.LightCurveBatch(three_lcs).LS(0.1, 10.0, 1e-3).run()
        assert isinstance(result, BatchResult)

    def test_batch_run_df_shape(self, three_lcs):
        result = vt.LightCurveBatch(three_lcs).LS(0.1, 10.0, 1e-3).run()
        assert len(result) == 3
        assert "LS_Period_1_0" in result.vars.columns

    def test_batch_run_capture_lc(self, three_lcs):
        result = vt.LightCurveBatch(three_lcs).LS(0.1, 10.0, 1e-3).run(capture_lc=True)
        assert len(result.lcs) == 3

    def test_batch_run_capture_lc_false(self, three_lcs):
        result = vt.LightCurveBatch(three_lcs).LS(0.1, 10.0, 1e-3).run(capture_lc=False)
        assert result.lcs == []

    def test_batch_per_lc_access(self, three_lcs):
        result = vt.LightCurveBatch(three_lcs).LS(0.1, 10.0, 1e-3).run()
        r0 = result[0]
        assert isinstance(r0, Result)
        assert "LS_Period_1_0" in r0.vars.index

    def test_batch_iter(self, three_lcs):
        result = vt.LightCurveBatch(three_lcs).LS(0.1, 10.0, 1e-3).run()
        items = list(result)
        assert len(items) == 3
        for r in items:
            assert isinstance(r, Result)

    def test_batch_run_ls_immediate(self, three_lcs):
        result = vt.LightCurveBatch(three_lcs).run_LS(0.1, 10.0, 1e-3)
        assert isinstance(result, BatchResult)
        assert len(result) == 3

    def test_batch_repr(self, three_lcs):
        batch = vt.LightCurveBatch(three_lcs).LS(0.1, 10.0, 1e-3)
        r = repr(batch)
        assert "LightCurveBatch" in r
        assert "n=3" in r

    def test_batch_per_lc_ok(self, three_lcs):
        result = vt.LightCurveBatch(three_lcs).LS(0.1, 10.0, 1e-3).run()
        for r in result:
            assert r.ok

    def test_batch_with_options(self, three_lcs):
        batch = vt.LightCurveBatch(three_lcs).with_options(capture_lc=False)
        result = batch.LS(0.1, 10.0, 1e-3).run()
        assert result.lcs == []

    def test_batch_slice(self, three_lcs):
        result = vt.LightCurveBatch(three_lcs).LS(0.1, 10.0, 1e-3).run()
        sub = result[0:2]
        assert isinstance(sub, BatchResult)
        assert len(sub) == 2


# ---------------------------------------------------------------------------
# BatchResult continuation chaining
# ---------------------------------------------------------------------------

class TestBatchResultContinuation:

    def test_batch_result_cmd_returns_batch(self, three_lcs):
        br = vt.LightCurveBatch(three_lcs).run_LS(0.1, 10.0, 1e-3)
        chain = br.LS(0.1, 10.0, 1e-3)
        assert isinstance(chain, LightCurveBatch)

    def test_batch_result_cmd_no_lcs_raises(self, three_lcs):
        br = vt.LightCurveBatch(three_lcs).run_LS(0.1, 10.0, 1e-3,
                                                   capture_lc=False)
        with pytest.raises(AttributeError, match="capture_lc"):
            br.LS(0.1, 10.0, 1e-3)

    def test_batch_result_run_cmd_returns_batch_result(self, three_lcs):
        br = vt.LightCurveBatch(three_lcs).run_LS(0.1, 10.0, 1e-3)
        br2 = br.run_LS(0.1, 10.0, 1e-3)
        assert isinstance(br2, BatchResult)
        assert len(br2) == 3

    def test_batch_result_run_cmd_no_lcs_raises(self, three_lcs):
        br = vt.LightCurveBatch(three_lcs).run_LS(0.1, 10.0, 1e-3,
                                                   capture_lc=False)
        with pytest.raises(AttributeError, match="capture_lc"):
            br.run_LS(0.1, 10.0, 1e-3)

    def test_batch_result_chain_and_run(self, three_lcs):
        br = vt.LightCurveBatch(three_lcs).run_LS(0.1, 10.0, 1e-3)
        br2 = br.LS(0.1, 10.0, 1e-3).run()
        assert isinstance(br2, BatchResult)
        assert len(br2) == 3
        assert "LS_Period_1_0" in br2.vars.columns


# ---------------------------------------------------------------------------
# Top-level vt.CMD(lc_input, ...) functions
# ---------------------------------------------------------------------------

class TestTopLevel:

    def test_lightcurve_input(self, lc):
        result = vt.LS(lc, 0.1, 10.0, 1e-3)
        assert isinstance(result, Result)
        assert result.ok
        assert result.varobjs.LS.Period_1 > 0

    def test_string_path_input(self):
        path = os.path.join(EXAMPLES_DIR, "2")
        result = vt.LS(path, 0.1, 10.0, 1e-3)
        assert isinstance(result, Result)
        assert result.ok
        assert result.varobjs.LS.Period_1 > 0

    def test_path_object_input(self):
        from pathlib import Path
        path = Path(EXAMPLES_DIR) / "2"
        result = vt.LS(path, 0.1, 10.0, 1e-3)
        assert isinstance(result, Result)
        assert result.ok

    def test_dataframe_input(self, lc):
        df = lc.to_dataframe()
        result = vt.LS(df, 0.1, 10.0, 1e-3)
        assert isinstance(result, Result)
        assert result.ok
        assert result.varobjs.LS.Period_1 > 0

    def test_2d_numpy_array_input(self, lc):
        arr = np.column_stack([lc.t, lc.mag, lc.err])
        result = vt.LS(arr, 0.1, 10.0, 1e-3)
        assert isinstance(result, Result)
        assert result.ok
        assert result.varobjs.LS.Period_1 > 0

    def test_tuple_of_arrays_input(self, lc):
        result = vt.LS((lc.t, lc.mag, lc.err), 0.1, 10.0, 1e-3)
        assert isinstance(result, Result)
        assert result.ok
        assert result.varobjs.LS.Period_1 > 0

    def test_second_command_rms(self, lc):
        result = vt.rms(lc)
        assert isinstance(result, Result)
        assert result.ok
        assert "RMS_0" in result.vars.index

    def test_capture_lc_false(self, lc):
        result = vt.LS(lc, 0.1, 10.0, 1e-3, capture_lc=False)
        assert result.lc is None

    def test_capture_lc_default_true(self, lc):
        result = vt.LS(lc, 0.1, 10.0, 1e-3)
        assert result.lc is not None

    def test_bad_input_type_raises(self, lc):
        with pytest.raises(TypeError, match="Cannot convert"):
            vt.LS(42, 0.1, 10.0, 1e-3)

    def test_1d_numpy_array_raises(self, lc):
        with pytest.raises(TypeError, match="2-D"):
            vt.LS(lc.t, 0.1, 10.0, 1e-3)

    def test_function_in_module_namespace(self):
        assert callable(vt.LS)
        assert callable(vt.rms)
        assert callable(vt.BLS)

    def test_function_has_docstring(self):
        assert vt.LS.__doc__ is not None
        assert "lc_input" in vt.LS.__doc__


# ---------------------------------------------------------------------------
# Per-LC array parameters on LightCurveBatch / BatchResult
# ---------------------------------------------------------------------------

class TestBatchPerLC:
    """Per-LC array parameters: different scalar value per LC in a batch.

    The same detection rules that Pipeline.run_batch() applies to PerLC /
    numpy arrays / pandas Series are applied by LightCurveBatch at run time,
    but without the -inlistvars restriction: any parameter (including those
    that do not accept "var varname" in vartools) can vary per LC.
    """

    # --- detection / validation (no binary needed) ---

    def test_detect_perlc_numpy(self, three_lcs):
        from pyvartools.commands import rms
        cmd = rms()
        # Attach a fake per-LC numpy array attribute (normally a real param)
        cmd.fake_param = np.array([1.0, 2.0, 3.0])
        batch = vt.LightCurveBatch(three_lcs, _commands=[cmd])
        perlc = batch._collect_perlc()
        assert (0, "fake_param") in perlc

    def test_detect_perlc_series(self, three_lcs):
        from pyvartools.commands import rms
        cmd = rms()
        cmd.fake_param = pd.Series([1.0, 2.0, 3.0])
        batch = vt.LightCurveBatch(three_lcs, _commands=[cmd])
        perlc = batch._collect_perlc()
        assert (0, "fake_param") in perlc

    def test_detect_perlc_explicit(self, three_lcs):
        from pyvartools.perlc import PerLC
        from pyvartools.commands import rms
        cmd = rms()
        cmd.fake_param = PerLC([1.0, 2.0, 3.0])
        batch = vt.LightCurveBatch(three_lcs, _commands=[cmd])
        perlc = batch._collect_perlc()
        assert (0, "fake_param") in perlc

    def test_plain_list_not_detected(self, three_lcs):
        """Plain Python lists are NOT auto-detected (avoids ld_coeffs ambiguity)."""
        from pyvartools.commands import rms
        cmd = rms()
        cmd.fake_param = [1.0, 2.0, 3.0]
        batch = vt.LightCurveBatch(three_lcs, _commands=[cmd])
        perlc = batch._collect_perlc()
        assert (0, "fake_param") not in perlc

    def test_length_mismatch_raises(self, three_lcs):
        from pyvartools.commands import LS
        cmd = LS(minp=np.array([0.1, 0.5]),  # 2 values for 3 LCs
                 maxp=10.0, subsample=0.1)
        batch = vt.LightCurveBatch(three_lcs, _commands=[cmd])
        perlc = batch._collect_perlc()
        with pytest.raises(ValueError, match="2 values but the batch has 3"):
            batch._validate_perlc(perlc)

    # --- extract scalar helpers ---

    def test_extract_perlc_positional(self, three_lcs):
        from pyvartools._batch import _extract_perlc_scalar
        s = pd.Series([10.0, 20.0, 30.0])
        assert _extract_perlc_scalar(s, 1, "star1") == 20.0

    def test_extract_perlc_name_based(self, three_lcs):
        from pyvartools._batch import _extract_perlc_scalar
        s = pd.Series([10.0, 20.0, 30.0], index=["1", "2", "3"])
        assert _extract_perlc_scalar(s, 0, "2") == 20.0  # name lookup wins

    def test_extract_perlc_name_fallback(self, three_lcs):
        from pyvartools._batch import _extract_perlc_scalar
        # Name not in index → fall back to positional
        s = pd.Series([10.0, 20.0, 30.0], index=["a", "b", "c"])
        assert _extract_perlc_scalar(s, 2, "not_there") == 30.0

    def test_extract_perlc_numpy(self, three_lcs):
        from pyvartools._batch import _extract_perlc_scalar
        arr = np.array([1.1, 2.2, 3.3])
        val = _extract_perlc_scalar(arr, 1, "star1")
        assert val == pytest.approx(2.2)
        assert isinstance(val, float)

    def test_extract_perlc_string_array(self, three_lcs):
        from pyvartools._batch import _extract_perlc_scalar
        arr = np.array(["ls", "aov", "fix 1.23"])
        assert _extract_perlc_scalar(arr, 2, "star3") == "fix 1.23"

    # --- end-to-end: per-LC parameters in LightCurveBatch.run() ---

    def test_batch_perlc_minp_numpy(self, three_lcs):
        """Different minp per LC via numpy array."""
        minp_arr = np.array([0.5, 0.5, 0.5])
        result = (
            vt.LightCurveBatch(three_lcs)
              .LS(minp=minp_arr, maxp=10.0, subsample=0.1)
              .run()
        )
        assert isinstance(result, BatchResult)
        assert len(result) == 3
        assert "LS_Period_1_0" in result.vars.columns

    def test_batch_perlc_minp_series(self, three_lcs):
        """Different minp per LC via pandas Series."""
        minp_series = pd.Series([0.5, 0.5, 0.5])
        result = (
            vt.LightCurveBatch(three_lcs)
              .LS(minp=minp_series, maxp=10.0, subsample=0.1)
              .run()
        )
        assert isinstance(result, BatchResult)
        assert len(result) == 3

    def test_batch_result_killharm_period_from_ls(self, three_lcs):
        """br.Killharm(period=br.vars[col]).run() — the primary use case."""
        br1 = vt.LightCurveBatch(three_lcs).run_LS(0.5, 10.0, 0.1)
        periods = br1.vars["LS_Period_1_0"]
        br2 = br1.Killharm(period=periods, nharm=2).run()
        assert isinstance(br2, BatchResult)
        assert len(br2) == 3
        # Killharm should produce output variables
        assert any("Killharm" in col for col in br2.vars.columns)

    def test_batch_result_killharm_period_named_series(self, three_lcs):
        """Name-indexed Series: values matched by LC name, not position."""
        br1 = vt.LightCurveBatch(three_lcs).run_LS(0.5, 10.0, 0.1)
        # Build Series indexed by LC name (reverse order to test that lookup
        # uses name, not position).
        names = list(br1.vars["Name"])
        periods = br1.vars["LS_Period_1_0"].values
        named = pd.Series(periods, index=names)
        named_reversed = named.iloc[::-1]  # reverse order
        br2 = br1.Killharm(period=named_reversed, nharm=2).run()
        assert isinstance(br2, BatchResult)
        assert len(br2) == 3

    def test_batch_perlc_length_mismatch_raises_at_run(self, three_lcs):
        """Length mismatch detected before any LC is processed."""
        wrong_arr = np.array([0.5, 0.5])  # 2 values, 3 LCs
        batch = vt.LightCurveBatch(three_lcs).LS(minp=wrong_arr, maxp=10.0,
                                                  subsample=0.1)
        with pytest.raises(ValueError, match="2 values but the batch has 3"):
            batch.run()


# ---------------------------------------------------------------------------
# Per-star scalar carry-forward (Step 1 of chained-state plumbing)
# ---------------------------------------------------------------------------

class TestScalars:
    """Tests for LightCurve.scalars, Result.lcscalars, BatchResult.lcscalars, and the
    parse_oneline_output split of ``VARTOOLS_SCALAR:`` rows.

    These cover the Python-only plumbing added in step 1 of the chained
    per-LC scalar feature.  Actual emission of ``VARTOOLS_SCALAR:`` lines
    requires the C-side ``-printallscalars`` option (step 2) and cross-
    chain referencing requires the C-side injection function (step 3);
    both are covered by separate end-to-end tests added in those steps.
    """

    def test_lightcurve_scalars_default_empty(self, lc):
        assert lc.scalars == {}

    def test_lightcurve_scalars_roundtrip_via_constructor(self):
        df = pd.DataFrame({"t": [1.0, 2.0], "mag": [10.0, 10.1], "err": [0.01, 0.01]})
        out = vt.LightCurve(df, name="x", scalars={"foo": 1.5, "bar": 42})
        assert out.scalars == {"foo": 1.5, "bar": 42}

    def test_lightcurve_scalars_roundtrip_via_from_arrays(self):
        out = vt.LightCurve.from_arrays(
            t=np.arange(5), mag=np.zeros(5), err=np.ones(5),
            scalars={"offset": 0.25},
        )
        assert out.scalars == {"offset": 0.25}

    def test_lightcurve_scalars_roundtrip_via_from_dataframe(self):
        df = pd.DataFrame({"t": [1.0], "mag": [10.0], "err": [0.01]})
        out = vt.LightCurve.from_dataframe(df, scalars={"z": 7.0})
        assert out.scalars == {"z": 7.0}

    def test_lightcurve_scalars_are_copied_not_aliased(self):
        s = {"a": 1.0}
        out = vt.LightCurve(pd.DataFrame({"t": [0.0]}), scalars=s)
        s["a"] = 999.0
        assert out.scalars == {"a": 1.0}

    def test_lightcurve_repr_includes_scalar_count(self):
        df = pd.DataFrame({"t": [0.0], "mag": [0.0], "err": [0.0]})
        out = vt.LightCurve(df, scalars={"x": 1.0, "y": 2.0})
        assert "scalars=2" in repr(out)

    def test_lightcurve_repr_omits_scalars_when_empty(self, lc):
        assert "scalars=" not in repr(lc)

    def test_result_lcscalars_empty_when_no_lc(self):
        r = Result(var=pd.Series({"Name": "foo"}), lc=None)
        assert r.lcscalars == {}

    def test_result_lcscalars_reads_from_captured_lc_scalars(self):
        df = pd.DataFrame({"t": [0.0], "mag": [0.0], "err": [0.0]})
        out_lc = vt.LightCurve(df, scalars={"myvar": 3.14})
        r = Result(var=pd.Series({"Name": "x"}), lc=out_lc)
        assert r.lcscalars == {"myvar": 3.14}

    def test_batchresult_lcscalars_empty_when_no_scalars(self):
        df = pd.DataFrame({"t": [0.0]})
        lcs = [vt.LightCurve(df, name=f"lc{i}") for i in range(3)]
        br = BatchResult(var=pd.DataFrame(), lcs=lcs)
        assert br.lcscalars.empty

    def test_batchresult_lcscalars_builds_dataframe(self):
        df = pd.DataFrame({"t": [0.0]})
        lcs = [
            vt.LightCurve(df, name="a", scalars={"s1": 1.0, "s2": 10.0}),
            vt.LightCurve(df, name="b", scalars={"s1": 2.0, "s2": 20.0}),
        ]
        br = BatchResult(var=pd.DataFrame(), lcs=lcs)
        out = br.lcscalars
        assert list(out.columns) == ["s1", "s2"]
        assert len(out) == 2
        assert out.iloc[0]["s1"] == 1.0
        assert out.iloc[1]["s2"] == 20.0

    # -- parse_oneline_output split helper ---------------------------------

    def test_parse_splits_vartools_scalar_rows(self):
        from pyvartools.results import parse_oneline_output, split_vars_and_scalars

        fake = (
            "Name = lc_a\n"
            "LS_Period_1_0 = 1.2345\n"
            "VARTOOLS_SCALAR:myvar = 42.0\n"
            "Name = lc_b\n"
            "LS_Period_1_0 = 6.789\n"
            "VARTOOLS_SCALAR:myvar = 99.0\n"
        )
        df = parse_oneline_output(fake)
        vars_df, scalars_df = split_vars_and_scalars(df)

        assert list(vars_df.columns) == ["Name", "LS_Period_1_0"]
        assert list(scalars_df.columns) == ["myvar"]
        assert len(vars_df) == 2
        assert scalars_df.iloc[0]["myvar"] == 42.0
        assert scalars_df.iloc[1]["myvar"] == 99.0

    def test_parse_no_scalar_rows_returns_empty_scalars_df(self):
        from pyvartools.results import parse_oneline_output, split_vars_and_scalars

        fake = "Name = lc\nLS_Period_1_0 = 1.0\n"
        df = parse_oneline_output(fake)
        _, scalars_df = split_vars_and_scalars(df)
        assert scalars_df.empty

    # -- chain carry-forward ----------------------------------------------

    def test_chain_carry_forward_populates_input_scalars(self, lc):
        """After `result.rms()`, the second segment's input LC carries the
        prior OUTCOLUMN values in its .scalars dict.

        Cannot observe this directly without hooking Pipeline, so we verify
        indirectly: the merged final Result's .lcscalars (from captured lc) at
        minimum contains the scalars carried in (library mode may also add
        harvested scalars once step 3 lands).
        """
        r1 = lc.LS(0.1, 10.0, 1e-3)
        r2 = r1.rms()
        # `r2.lc` carries the scalars that flowed through: prior LS_Period_1_0
        # should be among them, since _make_immediate_result injected it into
        # input_lc.scalars and Pipeline.run preserved it through to out_lc.
        assert r2.lc is not None
        assert "LS_Period_1_0" in r2.lcscalars
        assert r2.lcscalars["LS_Period_1_0"] == r1.vars["LS_Period_1_0"]


# ---------------------------------------------------------------------------
# Cross-chain expressions referencing injected scalars (Step 3)
# ---------------------------------------------------------------------------

class TestCrossChainScalars:
    """End-to-end tests verifying that -expr in a chained segment can
    reference OUTCOLUMN values produced by a prior segment.

    These exercise the full injection pipeline: pyvartools carries forward
    prior vars into the next run's lc.scalars, Pipeline.run prepends
    ``-expr const 'name=value'`` tokens and ``-columnsuffix <offset+i>``
    before each user command, vartools' -printallscalars dumps the new
    scalar state, and pyvartools parses it back into Result.lcscalars.
    """

    def test_chained_expr_references_prior_ls_period(self, lc):
        """`lc.LS().expr('doubled=2*LS_Period_1_0', vartype='scalar')` should
        produce doubled == 2 * r1.LS_Period_1_0."""
        r1 = lc.LS(0.5, 10.0, 0.1)
        period = r1.vars["LS_Period_1_0"]
        r2 = r1.expr("doubled=2*LS_Period_1_0", vartype="scalar")
        assert "doubled" in r2.lcscalars
        assert abs(r2.lcscalars["doubled"] - 2 * period) < 1e-9

    def test_chained_segment_suffix_shifts(self, lc):
        """Second segment's OUTCOLUMN suffixes continue the numbering."""
        r = lc.LS(0.1, 10.0, 1e-3).rms()
        # LS is segment 0 → suffix "_0", rms is segment 1 → suffix "_1".
        assert "LS_Period_1_0" in r.vars.index
        assert "RMS_1" in r.vars.index
        assert "Mean_Mag_1" in r.vars.index
        assert "Npoints_1" in r.vars.index
        # No collision: no "RMS_0" or "Mean_Mag_0" should appear.
        assert "RMS_0" not in r.vars.index
        assert "Mean_Mag_0" not in r.vars.index

    def test_chained_expr_then_further_chain(self, lc):
        """A value created by -expr scalar in segment 1 can be referenced in
        segment 2."""
        r1 = lc.LS(0.5, 10.0, 0.1)
        period = r1.vars["LS_Period_1_0"]
        r2 = r1.expr("halfper=LS_Period_1_0/2", vartype="scalar")
        # halfper is now in r2.lcscalars AND should be carried into the next segment
        r3 = r2.expr("quarterper=halfper/2", vartype="scalar")
        assert "quarterper" in r3.lcscalars
        assert abs(r3.lcscalars["quarterper"] - period / 4) < 1e-9

    def test_fresh_run_with_empty_scalars_unchanged(self, lc):
        """A non-chained run (no prior scalars, offset=0) produces output
        identical to what step-1 code produced — suffix still ``_0``, no
        injected -expr const tokens."""
        r = lc.LS(0.1, 10.0, 1e-3)
        assert "LS_Period_1_0" in r.vars.index
        # Nothing should be in .lcscalars for a fresh run with no prior scalars.
        assert r.lcscalars == {}

    def test_chained_via_result_method_merges_known_commands(self, lc):
        r = lc.LS(0.1, 10.0, 1e-3).rms()
        assert r._known_commands == ["LS", "rms"]

    def test_three_segment_chain(self, lc):
        """Three-segment chain keeps suffixes growing: _0, _1, _2."""
        r = lc.LS(0.1, 10.0, 1e-3).rms().stats("mag", "mean")
        assert "LS_Period_1_0" in r.vars.index
        assert "RMS_1" in r.vars.index
        # stats's output key ends in _2
        stats_keys = [k for k in r.vars.index if k.startswith("STATS_") or "mag" in k.lower()]
        assert any(k.endswith("_2") for k in stats_keys), (
            f"expected a stats key ending in _2, got: {stats_keys}"
        )
        assert r._known_commands == ["LS", "rms", "stats"]

    def test_batch_chain_continuation_suffix_shifts(self, three_lcs):
        """BatchResult → LightCurveBatch chain: each LC's rms output gets the
        shifted suffix."""
        br1 = vt.LightCurveBatch(three_lcs).LS(0.1, 10.0, 1e-3).run()
        br2 = br1.rms().run()
        # Every LC row should have both LS and rms columns, with shifted
        # suffixes on the rms side.
        assert "LS_Period_1_0" in br2.vars.columns
        assert "RMS_1" in br2.vars.columns
        assert "RMS_0" not in br2.vars.columns
        assert br2._known_commands == ["LS", "rms"]

    def test_batch_chain_inlistvars_per_lc_scalar(self, three_lcs):
        """Per-LC scalar values carried forward via -inlistvars: each LC
        in the continuation sees its own prior LS_Period_1_0, not a single
        shared constant.

        Verifies the user's concern directly — if the injection used
        -expr const, all three LCs' doubled values would equal 2 * LC[0]'s
        period, but with -inlistvars each gets its own.
        """
        br1 = vt.LightCurveBatch(three_lcs).LS(0.5, 10.0, 0.1).run()
        periods = list(br1.vars["LS_Period_1_0"])
        assert len(set(periods)) > 1, "test fixtures must have distinct LS periods"
        br2 = br1.expr("doubled=2*LS_Period_1_0", vartype="scalar").run()
        for i, lc in enumerate(br2.lcs):
            assert lc is not None
            assert "doubled" in lc.scalars
            assert abs(lc.scalars["doubled"] - 2 * periods[i]) < 1e-9, (
                f"LC {i}: got doubled={lc.scalars['doubled']}, expected "
                f"{2 * periods[i]} (2 * LS_Period_1_0={periods[i]})")


# ---------------------------------------------------------------------------
# Chained back-references — "ls" / "aov" / "bls" / "both" / "fixcolumn" when
# the prior chain step is in a separate vartools invocation.
# ---------------------------------------------------------------------------

class TestChainedBackReferences:

    def test_ls_backref_phase(self, lc):
        """lc.LS(...).Phase(period='ls') resolves 'ls' into the LS period."""
        r = lc.LS(0.5, 5.0, 0.001, npeaks=1).Phase(period="ls")
        assert r.lc is not None
        # After folding on the LS period the phase column starts near 0.
        assert 0.0 <= float(r.lc.t[0]) < 1.0

    def test_aov_backref_killharm(self, lc):
        """lc.aov(...).Killharm(period='aov') resolves 'aov' into the AOV period."""
        r_aov = lc.aov(0.5, 5.0, 0.1, 0.01, npeaks=1)
        r = r_aov.Killharm(period="aov", nharm=2)
        assert "Killharm_PeriodFix_1_1" in r.vars or "Killharm_Period_1_1" in r.vars \
            or any(k.startswith("Killharm_") for k in r.vars.index)

    def test_bls_backref_phase(self):
        """lc.BLS(...).Phase(period='bls', T0='bls 0.5') resolves both."""
        lc3 = vt.LightCurve.from_file(os.path.join(EXAMPLES_DIR, "3.transit"))
        r = (lc3.BLS(0.5, 5.0, rmin=0.01, rmax=0.1,
                     nbins=200, nfreq=20000, npeaks=1)
                .Phase(period="bls", T0="bls 0.5"))
        assert r.lc is not None

    def test_bls_backref_mandelagol_resolves_all_four_fields(self):
        """MandelAgolTransit(P0='bls') seeds period+T0+depth+qtran from BLS.

        Verifies the resolution alone — the subsequent MCMC fit can fail on
        numerical boundary conditions, but the important part is that the
        four transit initial-parameter fields get filled in from the prior
        BLS output.
        """
        import math
        lc3 = vt.LightCurve.from_file(os.path.join(EXAMPLES_DIR, "3.transit"))
        r_bls = lc3.BLS(0.5, 5.0, rmin=0.01, rmax=0.1,
                        nbins=200, nfreq=20000, npeaks=1)
        from pyvartools.commands import MandelAgolTransit
        cmd = MandelAgolTransit(P0="bls", T00=0)
        cmd._resolve_back_references(r_bls)
        # Period/T0 copied verbatim; Rp = sqrt(depth); a = 1/(qtran*pi).
        assert math.isclose(float(cmd.P0), float(r_bls.vars["BLS_Period_1_0"]))
        assert math.isclose(float(cmd.T00), float(r_bls.vars["BLS_Tc_1_0"]))
        assert math.isclose(float(cmd.r0),
                            math.sqrt(float(r_bls.vars["BLS_Depth_1_0"])))
        assert math.isclose(float(cmd.a0),
                            1.0 / (float(r_bls.vars["BLS_Qtran_1_0"]) * math.pi))

    def test_aov_ambiguity_picks_most_recent(self, lc):
        """aov().aov().aov_harm().Killharm(period='aov') picks aov_harm.

        This exercises the _position-based 'most-recent wins' rule — two
        -aov commands precede one -aov_harm, so the aov_harm period must
        be the one substituted into Killharm.
        """
        import math
        chain = (lc.aov(0.5, 5.0, 0.1, 0.01, npeaks=1)
                   .aov(0.3, 4.0, 0.1, 0.01, npeaks=1)
                   .aov_harm(2, 0.5, 5.0, 0.1, 0.001, npeaks=1))
        # Inspect resolution directly to avoid depending on Killharm's
        # output-column naming.
        from pyvartools.commands import Killharm
        cmd = Killharm(period="aov", nharm=2)
        cmd._resolve_back_references(chain)
        aov_harm_period = float(chain.vars["Period_1_2"])
        assert math.isclose(float(cmd.period), aov_harm_period, rel_tol=1e-6), (
            f"AOV back-reference should have picked aov_harm's period "
            f"({aov_harm_period}) but got {cmd.period}"
        )

    def test_fixcolumn_by_name(self, lc):
        """Phase(period='fixcolumn LS_Period_1_0') reads the named column."""
        r_ls = lc.LS(0.5, 5.0, 0.001, npeaks=1)
        r = r_ls.Phase(period="fixcolumn LS_Period_1_0")
        assert r.lc is not None

    def test_fixcolumn_missing_raises(self, lc):
        """An unknown column raises LookupError, not a segfault."""
        r_ls = lc.LS(0.5, 5.0, 0.001, npeaks=1)
        with pytest.raises(LookupError):
            r_ls.Phase(period="fixcolumn NOT_A_REAL_COLUMN")

    def test_missing_prior_cmd_raises(self, lc):
        """Phase(period='bls') without a prior BLS raises LookupError."""
        r_ls = lc.LS(0.5, 5.0, 0.001, npeaks=1)
        with pytest.raises(LookupError):
            r_ls.Phase(period="bls")

    def test_fixcolumn_numeric_index_rejected(self, lc):
        """A bare integer fixcolumn is not supported; user must pass a name."""
        r_ls = lc.LS(0.5, 5.0, 0.001, npeaks=1)
        with pytest.raises(ValueError, match="numeric column index"):
            r_ls.Phase(period="fixcolumn 3")

    def test_fixperiod_snr_ls_backref(self, lc):
        """lc.LS(...).LS(..., fixperiod_snr='ls') resolves into the prior LS period."""
        r = lc.LS(0.5, 5.0, 0.001, npeaks=1)
        r2 = r.LS(0.5, 5.0, 0.001, npeaks=1, fixperiod_snr="ls")
        assert any("PeriodFix" in k for k in r2.vars.index)

    def test_killharm_both_resolves_two_periods(self, lc):
        """Killharm(period='both') pulls the LS period AND the AOV period.

        When resolved on a prior Result, the substituted value becomes a
        tuple (ls_period, aov_period), and the CLI emission carries both.
        """
        import math
        r = (lc.LS(0.5, 5.0, 0.001, npeaks=1)
                .aov(0.5, 5.0, 0.1, 0.01, npeaks=1))
        from pyvartools.commands import Killharm
        cmd = Killharm(period="both", nharm=2)
        cmd._resolve_back_references(r)
        assert isinstance(cmd.period, tuple) and len(cmd.period) == 2
        ls_p, aov_p = (float(x) for x in cmd.period)
        assert math.isclose(ls_p, float(r.vars["LS_Period_1_0"]), rel_tol=1e-6)
        assert math.isclose(aov_p, float(r.vars["Period_1_1"]), rel_tol=1e-6)
