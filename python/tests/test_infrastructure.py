"""
Phase 1 infrastructure tests.

These tests require a working vartools binary on PATH or configured via
VARTOOLS_BINARY / pyvartools.config.set_binary().
"""

import os
import numpy as np
import pandas as pd
import pytest

import pyvartools as vt
from pyvartools.results import parse_oneline_output
from pyvartools.pipeline import _inputlcformat_from_df, _inputlcformat_from_spec

VARTOOLS_SRC = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", ".."))
EXAMPLES_DIR = os.path.join(VARTOOLS_SRC, "EXAMPLES")
EXAMPLE_LC = os.path.join(EXAMPLES_DIR, "2")


# ---------------------------------------------------------------------------
# Binary discovery
# ---------------------------------------------------------------------------

def test_get_binary_finds_something():
    binary = vt.get_binary()
    assert os.path.isfile(binary)
    assert os.access(binary, os.X_OK)


def test_set_binary_override(tmp_path):
    # Create a fake executable
    fake = tmp_path / "vartools"
    fake.write_text("#!/bin/sh\necho fake\n")
    fake.chmod(0o755)
    vt.set_binary(str(fake))
    assert vt.get_binary() == str(fake)
    # Reset
    vt.set_binary("")


def test_set_binary_missing_raises():
    vt.set_binary("/nonexistent/path/vartools")
    with pytest.raises(FileNotFoundError):
        vt.get_binary()
    vt.set_binary("")


# ---------------------------------------------------------------------------
# LightCurve I/O
# ---------------------------------------------------------------------------

def test_lightcurve_from_file():
    lc = vt.LightCurve.from_file(EXAMPLE_LC)
    assert len(lc) > 0
    assert lc.t.ndim == 1
    assert lc.mag.ndim == 1
    assert lc.err.ndim == 1


def test_lightcurve_from_arrays():
    t = np.linspace(0, 100, 200)
    mag = np.random.normal(10.0, 0.01, 200)
    err = np.full(200, 0.01)
    lc = vt.LightCurve.from_arrays(t, mag, err, name="test")
    assert len(lc) == 200
    assert lc.name == "test"


def test_lightcurve_from_dataframe():
    df = pd.DataFrame({"t": [1.0, 2.0], "mag": [10.0, 10.1], "err": [0.01, 0.01]})
    lc = vt.LightCurve.from_dataframe(df)
    assert len(lc) == 2


def test_lightcurve_missing_standard_cols_allowed():
    """t, mag, err are no longer required — partial columns are valid."""
    # Only t and mag (no err)
    lc = vt.LightCurve(pd.DataFrame({"t": [1.0, 2.0], "mag": [10.0, 10.1]}))
    assert lc.err is None
    assert lc.t is not None

    # Only custom columns (no standard cols at all)
    lc2 = vt.LightCurve(pd.DataFrame({"col1": [1.0], "col2": [2.0]}))
    assert lc2.t is None
    assert lc2.mag is None
    assert lc2.err is None


def test_lightcurve_properties_return_none_when_missing():
    lc = vt.LightCurve(pd.DataFrame({"col1": [1.0, 2.0]}))
    assert lc.t is None
    assert lc.mag is None
    assert lc.err is None


def test_lightcurve_from_arrays_partial():
    """from_arrays allows None for any of the three standard arrays."""
    t   = np.linspace(0, 10, 50)
    mag = np.full(50, 12.0)
    # No err
    lc = vt.LightCurve.from_arrays(t=t, mag=mag)
    assert lc.err is None
    assert len(lc) == 50

    # No t, mag, or err — only aux columns
    airmass = np.ones(30)
    lc2 = vt.LightCurve.from_arrays(aux={"airmass": airmass})
    assert lc2.t is None
    assert list(lc2._df.columns) == ["airmass"]


def test_lightcurve_roundtrip_tempfile(tmp_path):
    lc = vt.LightCurve.from_file(EXAMPLE_LC)
    path = lc._to_tempfile(dir=str(tmp_path))
    lc2 = vt.LightCurve.from_file(path)
    np.testing.assert_allclose(lc.t, lc2.t, rtol=1e-8)
    np.testing.assert_allclose(lc.mag, lc2.mag, rtol=1e-8)
    np.testing.assert_allclose(lc.err, lc2.err, rtol=1e-8)
    os.unlink(path)


def test_lightcurve_to_dataframe():
    lc = vt.LightCurve.from_file(EXAMPLE_LC)
    df = lc.to_dataframe()
    assert isinstance(df, pd.DataFrame)
    assert set(["t", "mag", "err"]).issubset(df.columns)


def test_lightcurve_to_arrays():
    lc = vt.LightCurve.from_file(EXAMPLE_LC)
    t, mag, err = lc.to_arrays()
    assert len(t) == len(lc)


# ---------------------------------------------------------------------------
# stdout parser
# ---------------------------------------------------------------------------

def test_parse_oneline_single():
    stdout = (
        "Name = EXAMPLES/2\n"
        "LS_Period_1_0            =     1.23440877\n"
        "Log10_LS_Prob_1_0        = -4000.59209\n"
    )
    df = parse_oneline_output(stdout)
    assert len(df) == 1
    assert df["Name"].iloc[0] == "EXAMPLES/2"
    assert abs(df["LS_Period_1_0"].iloc[0] - 1.23440877) < 1e-6


def test_parse_oneline_multiple():
    stdout = (
        "Name = lc1\n"
        "RMS = 0.005\n"
        "Name = lc2\n"
        "RMS = 0.007\n"
    )
    df = parse_oneline_output(stdout)
    assert len(df) == 2
    assert list(df["Name"]) == ["lc1", "lc2"]


def test_parse_oneline_empty():
    df = parse_oneline_output("")
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 0


# ---------------------------------------------------------------------------
# Pipeline — minimal end-to-end (requires binary)
# ---------------------------------------------------------------------------

class _RMSCommand(vt.commands.VartoolsCommand):
    """Minimal concrete command for testing: -rms"""
    _vt_name = "rms"

    def _to_cli_args(self):
        return ["-rms"]


def test_pipeline_run_single_lc():
    lc = vt.LightCurve.from_file(EXAMPLE_LC)
    pipe = vt.Pipeline([_RMSCommand()])
    result = pipe.run(lc)
    assert isinstance(result.vars, pd.Series)
    assert "Weighted_RMS" in result.vars.index or any(
        "rms" in c.lower() or "RMS" in c for c in result.vars.index
    )


def test_pipeline_run_batch():
    lcs = [vt.LightCurve.from_file(EXAMPLE_LC) for _ in range(3)]
    pipe = vt.Pipeline([_RMSCommand()])
    result = pipe.run_batch(lcs)
    assert isinstance(result.vars, pd.DataFrame)
    assert len(result.vars) == 3


def test_pipeline_bad_command_raises():
    class _BadCmd(vt.commands.VartoolsCommand):
        _vt_name = "notacommand"
        def _to_cli_args(self):
            return ["-notacommand"]

    lc = vt.LightCurve.from_file(EXAMPLE_LC)
    pipe = vt.Pipeline([_BadCmd()])
    with pytest.raises(vt.RunError):
        pipe.run(lc)


# ---------------------------------------------------------------------------
# capture_lc in batch methods
# ---------------------------------------------------------------------------

def _make_lc(n=100, period=2.0):
    import numpy as np
    t = np.linspace(0, 20, n)
    mag = 10.0 + 0.1 * np.sin(2 * np.pi * t / period)
    err = np.full(n, 0.01)
    return vt.LightCurve.from_arrays(t, mag, err)


def test_run_batch_capture_lc():
    lcs = [_make_lc(period=1.5 + i * 0.3) for i in range(3)]
    result = vt.Pipeline([vt.commands.clip(sigclip=5.0)]).run_batch(
        lcs, capture_lc=True
    )
    assert result.lcs is not None
    assert len(result.lcs) == 3
    for out_lc in result.lcs:
        assert out_lc is not None
        assert len(out_lc) > 0


def test_run_batch_capture_lc_names_preserved():
    lcs = [_make_lc() for _ in range(2)]
    lcs[0].name = "alpha"
    lcs[1].name = "beta"
    result = vt.Pipeline([vt.commands.clip(sigclip=5.0)]).run_batch(
        lcs, capture_lc=True
    )
    assert result.lcs[0].name == "alpha"
    assert result.lcs[1].name == "beta"


def test_run_filelist_capture_lc(tmp_path):
    lcs = [_make_lc(period=1.5 + i * 0.3) for i in range(3)]
    paths = []
    for i, lc in enumerate(lcs):
        p = tmp_path / f"lc{i}.lc"
        lc._df.to_csv(str(p), sep=" ", header=False, index=False,
                      float_format="%.10f")
        paths.append(p)

    result = vt.Pipeline([vt.commands.clip(sigclip=5.0)]).run_filelist(
        paths, capture_lc=True
    )
    assert result.lcs is not None
    assert len(result.lcs) == 3
    for out_lc in result.lcs:
        assert out_lc is not None


# ---------------------------------------------------------------------------
# raise_on_error / BatchResult.ok
# ---------------------------------------------------------------------------

def test_run_batch_raise_on_error_false():
    class _BadCmd(vt.commands.VartoolsCommand):
        _vt_name = "notacommand"
        def _to_cli_args(self):
            return ["-notacommand"]

    lcs = [_make_lc() for _ in range(2)]
    result = vt.Pipeline([_BadCmd()]).run_batch(lcs, raise_on_error=False)
    assert not result.ok
    assert result.error is not None
    assert isinstance(result.error, vt.RunError)
    assert result.vars.empty


def test_run_batch_ok_true_on_success():
    lcs = [_make_lc() for _ in range(2)]
    result = vt.Pipeline([vt.commands.rms()]).run_batch(lcs)
    assert result.ok
    assert result.error is None


# ---------------------------------------------------------------------------
# timeout
# ---------------------------------------------------------------------------

def test_run_timeout_raises():
    lc = vt.LightCurve.from_file(EXAMPLE_LC)
    # timeout=0 should always expire immediately
    with pytest.raises(vt.RunError, match="timed out"):
        vt.Pipeline([vt.commands.rms()]).run(lc, timeout=0)


# ---------------------------------------------------------------------------
# FITS I/O
# ---------------------------------------------------------------------------

def test_lightcurve_from_fits_roundtrip(tmp_path):
    try:
        from astropy.io import fits
        import numpy as np
    except ImportError:
        pytest.skip("astropy not available")

    t   = np.linspace(0, 10, 50)
    mag = 10.0 + 0.05 * np.sin(2 * np.pi * t / 2.3)
    err = np.full(50, 0.01)

    col_t   = fits.Column(name="BJD",  format="D", array=t)
    col_mag = fits.Column(name="Mag",  format="D", array=mag)
    col_err = fits.Column(name="Err",  format="D", array=err)
    hdu = fits.BinTableHDU.from_columns([col_t, col_mag, col_err])
    fpath = tmp_path / "star.fits"
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(str(fpath))

    lc = vt.LightCurve.from_file(fpath)   # auto-detect .fits
    assert len(lc) == 50
    np.testing.assert_allclose(lc.t,   t,   atol=1e-8)
    np.testing.assert_allclose(lc.mag, mag, atol=1e-8)
    np.testing.assert_allclose(lc.err, err, atol=1e-8)


def test_lightcurve_from_fits_custom_cols(tmp_path):
    try:
        from astropy.io import fits
        import numpy as np
    except ImportError:
        pytest.skip("astropy not available")

    t   = np.linspace(0, 5, 30)
    mag = np.full(30, 12.0)
    err = np.full(30, 0.005)

    col_t   = fits.Column(name="TIME",  format="D", array=t)
    col_mag = fits.Column(name="FLUX",  format="D", array=mag)
    col_err = fits.Column(name="EFLUX", format="D", array=err)
    hdu = fits.BinTableHDU.from_columns([col_t, col_mag, col_err])
    fpath = tmp_path / "lc.fits"
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(str(fpath))

    lc = vt.LightCurve.from_file(
        fpath, format="fits",
        t_col="TIME", mag_col="FLUX", err_col="EFLUX"
    )
    assert len(lc) == 30
    assert lc.name == "lc"


def test_lightcurve_from_fits_missing_col(tmp_path):
    try:
        from astropy.io import fits
        import numpy as np
    except ImportError:
        pytest.skip("astropy not available")

    col = fits.Column(name="X", format="D", array=np.zeros(10))
    hdu = fits.BinTableHDU.from_columns([col])
    fpath = tmp_path / "bad.fits"
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(str(fpath))

    with pytest.raises(ValueError, match="BJD"):
        vt.LightCurve.from_file(fpath)


# ---------------------------------------------------------------------------
# -inputlcformat helpers (unit tests)
# ---------------------------------------------------------------------------

def test_inputlcformat_from_df_no_extra():
    """Exactly t, mag, err in default order → None (no format string needed)."""
    df = pd.DataFrame({"t": [1.0], "mag": [10.0], "err": [0.01]})
    assert _inputlcformat_from_df(df.columns) is None


def test_inputlcformat_from_df_extra_cols():
    df = pd.DataFrame({
        "t": [1.0], "mag": [10.0], "err": [0.01], "airmass": [1.2], "fwhm": [2.5]
    })
    fmt = _inputlcformat_from_df(df.columns)
    assert fmt == "t:1,mag:2,err:3,airmass:4,fwhm:5"


def test_inputlcformat_from_df_missing_standard_cols():
    """Missing t/mag/err → format string is emitted."""
    df = pd.DataFrame({"t": [1.0], "mag": [10.0]})   # no err
    fmt = _inputlcformat_from_df(df.columns)
    assert fmt == "t:1,mag:2"

    df2 = pd.DataFrame({"airmass": [1.2], "xpos": [100.0]})
    fmt2 = _inputlcformat_from_df(df2.columns)
    assert fmt2 == "airmass:1,xpos:2"


def test_inputlcformat_from_spec_list():
    fmt = _inputlcformat_from_spec(["t", "mag", "err", "airmass"])
    assert fmt == "t:1,mag:2,err:3,airmass:4"


def test_inputlcformat_from_spec_dict_numeric():
    fmt = _inputlcformat_from_spec({"t": 1, "mag": 2, "err": 3, "airmass": 4})
    assert fmt == "t:1,mag:2,err:3,airmass:4"


def test_inputlcformat_from_spec_dict_fits_names():
    fmt = _inputlcformat_from_spec({"t": "BJD_TDB", "mag": "MAG", "err": "ERR"})
    assert fmt == "t:BJD_TDB,mag:MAG,err:ERR"


def test_inputlcformat_from_spec_lccolumn_string_type():
    """LCColumn lets the user attach a non-default type (e.g. ``"string"``)
    to a column, which is required for things like the fiphot string flag
    consumed by -hatpiflag.
    """
    fmt = _inputlcformat_from_spec({
        "t": 1, "mag": 2, "err": 3,
        "fiphot_flag": vt.LCColumn(col=4, type="string"),
    })
    assert fmt == "t:1,mag:2,err:3,fiphot_flag:4:string"


def test_inputlcformat_from_spec_lccolumn_utc_with_format():
    """LCColumn(type="utc", format=...) emits the format string as the
    fourth ``:``-separated field, matching vartools' own grammar.
    """
    fmt = _inputlcformat_from_spec({
        "t": vt.LCColumn(col=1, type="utc", format="%Y-%M-%DT%h:%m:%s"),
        "mag": 2, "err": 3,
    })
    assert fmt == "t:1:utc:%Y-%M-%DT%h:%m:%s,mag:2,err:3"


def test_inputlcformat_from_spec_lccolumn_default_type():
    """LCColumn with the default type emits ``name:col:double`` so the
    format string is unambiguous, not silently dropped.
    """
    fmt = _inputlcformat_from_spec({"airmass": vt.LCColumn(col=4)})
    # default type field is preserved so callers see what gets emitted.
    assert fmt == "airmass:4:double"


def test_inputlcformat_from_spec_mixed_lccolumn_and_int():
    """Bare ints/strings still work alongside LCColumn instances."""
    fmt = _inputlcformat_from_spec({
        "t": 1,
        "mag": 2,
        "err": 3,
        "x": "XIC",                                 # FITS column name
        "fiphot_flag": vt.LCColumn(col=4, type="string"),
    })
    assert fmt == ("t:1,mag:2,err:3,x:XIC,"
                   "fiphot_flag:4:string")


# ---------------------------------------------------------------------------
# -inlistvars helpers (unit tests)
# ---------------------------------------------------------------------------

from pyvartools.pipeline import _inlistvars_from_spec


def test_inlistvars_int_shorthand():
    spec = _inlistvars_from_spec({"minp": 2, "maxp": 3})
    assert spec == "minp:2,maxp:3"


def test_inlistvars_listvar_default_type():
    spec = _inlistvars_from_spec({"minp": vt.ListVar(col=2)})
    # The default-type field is preserved so the emitted token is unambiguous.
    assert spec == "minp:2:double"


def test_inlistvars_listvar_string_type():
    spec = _inlistvars_from_spec({"name": vt.ListVar(col=1, type="string")})
    assert spec == "name:1:string"


def test_inlistvars_listvar_init():
    """ListVar(col=0, init=...) creates a per-star variable from an
    expression (no list column) and emits the init expression as the
    fourth ``:``-separated field.
    """
    spec = _inlistvars_from_spec(
        {"seq": vt.ListVar(col=0, type="int", init="NF")}
    )
    assert spec == "seq:0:int:NF"


def test_inlistvars_combinelc_keyword_emitted():
    """ListVar(combinelc=True) inserts the literal ``combinelc`` token
    between the column number and the type, so vartools knows the
    column carries one comma-joined value per combined file (used by
    ``-l ... combinelcs`` mode).
    """
    spec = _inlistvars_from_spec(
        {"trendlist": vt.ListVar(col=2, type="string", combinelc=True)}
    )
    assert spec == "trendlist:2:combinelc:string"


def test_inlistvars_combinelc_default_type():
    """combinelc=True with the default ``"double"`` type still emits
    the type field (so the keyword position is unambiguous).
    """
    spec = _inlistvars_from_spec(
        {"x": vt.ListVar(col=2, combinelc=True)}
    )
    assert spec == "x:2:combinelc:double"


def test_inlistvars_combinelc_init():
    """combinelc=True composes with the init field, with combinelc
    appearing immediately after the column number — same order as the
    CLI grammar.
    """
    spec = _inlistvars_from_spec(
        {"v": vt.ListVar(col=0, type="int", init="NF", combinelc=True)}
    )
    assert spec == "v:0:combinelc:int:NF"


def test_inlistvars_mixed_combinelc_and_int_shorthand():
    """Bare ints continue to work alongside combinelc-qualified ListVars."""
    spec = _inlistvars_from_spec({
        "minp": 2,
        "trendlist": vt.ListVar(col=3, type="string", combinelc=True),
    })
    assert spec == "minp:2,trendlist:3:combinelc:string"


# ---------------------------------------------------------------------------
# -inputlcformat integration tests (require binary)
# ---------------------------------------------------------------------------

def _make_lc_with_aux(n=80, period=2.0):
    """LightCurve with an extra 'airmass' column."""
    t   = np.linspace(0, 20, n)
    mag = 10.0 + 0.1 * np.sin(2 * np.pi * t / period)
    err = np.full(n, 0.01)
    airmass = 1.0 + 0.5 * np.abs(np.sin(np.pi * t / 10.0))
    return vt.LightCurve.from_arrays(t, mag, err, aux={"airmass": airmass})


def test_run_with_extra_columns():
    """run() auto-discovers extra columns and injects -inputlcformat."""
    lc = _make_lc_with_aux()
    assert "airmass" in lc._df.columns
    result = vt.Pipeline([vt.commands.rms()]).run(lc)
    assert isinstance(result.vars, pd.Series)
    assert "Weighted_RMS" in result.vars.index or any(
        "rms" in c.lower() or "RMS" in c for c in result.vars.index
    )


def test_run_batch_with_extra_columns():
    """run_batch() auto-discovers extra columns from the first LC."""
    lcs = [_make_lc_with_aux(period=1.5 + i * 0.5) for i in range(3)]
    result = vt.Pipeline([vt.commands.rms()]).run_batch(lcs)
    assert isinstance(result.vars, pd.DataFrame)
    assert len(result.vars) == 3


def test_run_file_with_columns_list(tmp_path):
    """run_file() with a list-style columns spec."""
    lc = _make_lc_with_aux()
    p = tmp_path / "lc_aux.lc"
    lc._df.to_csv(str(p), sep=" ", header=False, index=False, float_format="%.10f")

    result = vt.Pipeline([vt.commands.rms()]).run_file(
        p, columns=["t", "mag", "err", "airmass"]
    )
    assert isinstance(result.vars, pd.Series)


def test_run_file_with_columns_dict(tmp_path):
    """run_file() with a dict-style columns spec."""
    lc = _make_lc_with_aux()
    p = tmp_path / "lc_aux2.lc"
    lc._df.to_csv(str(p), sep=" ", header=False, index=False, float_format="%.10f")

    result = vt.Pipeline([vt.commands.rms()]).run_file(
        p, columns={"t": 1, "mag": 2, "err": 3, "airmass": 4}
    )
    assert isinstance(result.vars, pd.Series)


def test_run_filelist_with_columns(tmp_path):
    """run_filelist() with a columns spec."""
    lcs = [_make_lc_with_aux(period=1.5 + i * 0.5) for i in range(3)]
    paths = []
    for i, lc in enumerate(lcs):
        p = tmp_path / f"lc_aux_{i}.lc"
        lc._df.to_csv(str(p), sep=" ", header=False, index=False, float_format="%.10f")
        paths.append(p)

    result = vt.Pipeline([vt.commands.rms()]).run_filelist(
        paths, columns=["t", "mag", "err", "airmass"]
    )
    assert isinstance(result.vars, pd.DataFrame)
    assert len(result.vars) == 3


# ---------------------------------------------------------------------------
# cmd.o capture
# ---------------------------------------------------------------------------

def test_o_requires_filename_or_capture():
    """cmd.o() without arguments should raise ValueError."""
    with pytest.raises(ValueError, match="requires either"):
        vt.commands.o()


def test_o_capture_single_no_filename():
    """capture=True with no filename: temp file cleaned up, LC in result.files."""
    lc = _make_lc()
    pipe = vt.Pipeline([
        vt.commands.clip(5.0),
        vt.commands.o(capture=True),
    ])
    result = pipe.run(lc)
    assert "o" in result.files
    captured = result.files["o"]
    assert isinstance(captured, vt.LightCurve)
    assert len(captured) > 0


def test_o_capture_single_custom_key():
    """capture=True with a custom key name."""
    lc = _make_lc()
    pipe = vt.Pipeline([
        vt.commands.o(capture=True, key="clipped_lc"),
    ])
    result = pipe.run(lc)
    assert "clipped_lc" in result.files
    assert isinstance(result.files["clipped_lc"], vt.LightCurve)


def test_o_capture_single_with_explicit_filename(tmp_path):
    """capture=True with an explicit filename writes to disk AND captures."""
    out_file = str(tmp_path / "output.lc")
    lc = _make_lc()
    pipe = vt.Pipeline([
        vt.commands.o(filename=out_file, capture=True),
    ])
    result = pipe.run(lc)
    # File written to disk
    assert os.path.isfile(out_file)
    # Also captured in result.files
    assert "o" in result.files
    assert isinstance(result.files["o"], vt.LightCurve)


def test_o_no_capture_no_disk_write_in_tmp():
    """capture=False with no filename should raise ValueError (not silently ignored)."""
    with pytest.raises(ValueError):
        vt.commands.o()


def test_o_capture_batch_no_filename():
    """capture=True in batch mode: result.files["o"] is a list of LightCurves."""
    lcs = [_make_lc(period=1.5 + i * 0.3) for i in range(3)]
    pipe = vt.Pipeline([
        vt.commands.clip(5.0),
        vt.commands.o(capture=True),
    ])
    batch = pipe.run_batch(lcs)
    assert "o" in batch.files
    captured = batch.files["o"]
    assert len(captured) == 3
    for lc in captured:
        assert isinstance(lc, vt.LightCurve)
        assert len(lc) > 0


def test_o_multiple_captures_different_keys():
    """Two cmd.o(capture=True) commands with distinct keys both appear in files."""
    lc = _make_lc()
    pipe = vt.Pipeline([
        vt.commands.clip(5.0),
        vt.commands.o(capture=True, key="after_clip"),
        vt.commands.rms(),
        vt.commands.o(capture=True, key="final"),
    ])
    result = pipe.run(lc)
    assert "after_clip" in result.files
    assert "final" in result.files
    assert isinstance(result.files["after_clip"], vt.LightCurve)
    assert isinstance(result.files["final"], vt.LightCurve)


# ---------------------------------------------------------------------------
# FITS-header round-trip via LightCurve.fitsheader (.from_file / .to_file).
# ---------------------------------------------------------------------------

FITS_LC = os.path.join(EXAMPLES_DIR, "10.fits")


def test_fitsheader_captured_on_read():
    """Reading a FITS light curve populates `fitsheader` with a Header."""
    from astropy.io import fits
    lc = vt.LightCurve.from_file(FITS_LC, t_col="time",
                                 mag_col="mag", err_col="err")
    assert isinstance(lc.fitsheader, fits.Header)


def test_fitsheader_filters_structural_on_read():
    """Structural keywords (TTYPE, TFORM, NAXIS, TFIELDS…) are removed."""
    lc = vt.LightCurve.from_file(FITS_LC, t_col="time",
                                 mag_col="mag", err_col="err")
    for structural in ("TFIELDS", "NAXIS1", "NAXIS2",
                       "TTYPE1", "TFORM1", "TTYPE2", "TUNIT1",
                       "XTENSION", "SIMPLE", "BITPIX"):
        assert structural not in lc.fitsheader, (
            f"{structural} leaked into preserved header"
        )


def test_fitsheader_ascii_source_is_none():
    """Reading an ASCII light curve leaves fitsheader = None."""
    lc = vt.LightCurve.from_file(EXAMPLE_LC)
    assert lc.fitsheader is None


def test_fitsheader_roundtrip(tmp_path):
    """Round-trip of observational keywords + COMMENT + HISTORY."""
    from astropy.io import fits
    lc = vt.LightCurve.from_file(FITS_LC, t_col="time",
                                 mag_col="mag", err_col="err")
    lc.fitsheader["TELESCOP"] = "HATPI"
    lc.fitsheader["OBJECT"]   = "V1234 Cyg"
    lc.fitsheader["EPOCH"]    = (2000.0, "J2000")
    lc.fitsheader.add_comment("round-trip COMMENT")
    lc.fitsheader.add_history("round-trip HISTORY")

    out = str(tmp_path / "out.fits")
    lc.to_file(out)

    lc2 = vt.LightCurve.from_file(out, t_col="t",
                                  mag_col="mag", err_col="err")
    assert lc2.fitsheader["TELESCOP"] == "HATPI"
    assert lc2.fitsheader["OBJECT"]   == "V1234 Cyg"
    assert float(lc2.fitsheader["EPOCH"]) == 2000.0
    comments = [str(c.value) for c in lc2.fitsheader.cards
                if c.keyword == "COMMENT"]
    histories = [str(c.value) for c in lc2.fitsheader.cards
                 if c.keyword == "HISTORY"]
    assert any("round-trip COMMENT" in c for c in comments)
    assert any("round-trip HISTORY" in h for h in histories)


def test_fitsheader_none_passthrough_to_fits(tmp_path):
    """A LightCurve with fitsheader=None still writes to FITS correctly."""
    from astropy.io import fits
    t = np.linspace(0, 10, 50)
    lc = vt.LightCurve.from_arrays(t, np.sin(t), np.full(50, 0.01))
    assert lc.fitsheader is None
    out = str(tmp_path / "blank.fits")
    lc.to_file(out)
    with fits.open(out) as hdul:
        # No observational keywords injected; only the astropy default
        # primary-HDU structure keywords.
        keys = list(hdul[0].header.keys())
        assert all(k in {"SIMPLE", "BITPIX", "NAXIS", "EXTEND", ""}
                   for k in keys), f"unexpected primary header keys: {keys}"
        # Data wrote correctly on the extension.
        assert hdul[1].header["TFIELDS"] == 3
        assert hdul[1].header["NAXIS2"] == 50


def test_fitsheader_structural_keys_in_user_header_dropped_on_write(tmp_path):
    """If the user attaches a Header that includes structural keys,
    those are silently dropped on write (they'd otherwise collide with
    astropy's auto-generated structure)."""
    from astropy.io import fits
    t = np.linspace(0, 10, 50)
    lc = vt.LightCurve.from_arrays(t, np.sin(t), np.full(50, 0.01))
    hdr = fits.Header()
    hdr["TELESCOP"] = "HATPI"
    hdr["TFIELDS"]  = 99           # structural — should be dropped
    hdr["TTYPE1"]   = "should_not_be_preserved"
    hdr["NAXIS2"]   = 999           # structural — should be dropped
    lc.fitsheader = hdr

    out = str(tmp_path / "out.fits")
    lc.to_file(out)
    with fits.open(out) as hdul:
        assert hdul[0].header.get("TELESCOP") == "HATPI"
        # The structural leaks must not appear in either HDU's header
        # with the user's bogus values:
        assert "TTYPE1" not in hdul[0].header
        assert hdul[1].header["TFIELDS"] == 3      # astropy's real value
        assert hdul[1].header["NAXIS2"] == 50
