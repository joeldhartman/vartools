"""Phase 2 tests: command CLI-arg generation and end-to-end execution."""

import os
import sys

import numpy as np
import pandas as pd
import pytest

# Locate the source tree so we can run without installing
VARTOOLS_SRC = os.path.realpath(
    os.path.join(os.path.dirname(__file__), "..", "..")
)
PYTHON_DIR = os.path.join(VARTOOLS_SRC, "python")
if PYTHON_DIR not in sys.path:
    sys.path.insert(0, PYTHON_DIR)

import pyvartools as vt
from pyvartools import commands as cmd

# Point to the installed binary if not overridden via environment.
# Search order: VARTOOLS_BINARY env → <src>/../../../bin/vartools (HATpipe install) → PATH.
if not os.environ.get("VARTOOLS_BINARY"):
    _candidates = [
        # HATpipe-style install: source/vartools → bin/vartools two levels up
        os.path.join(VARTOOLS_SRC, "..", "..", "bin", "vartools"),
        # In-tree wrapper after make (may need LD_LIBRARY_PATH)
        os.path.join(VARTOOLS_SRC, "src", "vartools"),
    ]
    for _cand in _candidates:
        _cand = os.path.realpath(_cand)
        if os.path.isfile(_cand) and os.access(_cand, os.X_OK):
            vt.set_binary(_cand)
            break

try:
    _HAVE_BINARY = os.path.isfile(vt.get_binary())
except FileNotFoundError:
    _HAVE_BINARY = False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_sinusoidal_lc(n=200, period=1.5, noise=0.01):
    """Return a LightCurve with a sinusoidal signal of known period."""
    t = np.linspace(0, 30, n)
    mag = 10.0 + 0.1 * np.sin(2 * np.pi * t / period)
    err = np.full(n, noise)
    return vt.LightCurve.from_arrays(t, mag, err)


# ---------------------------------------------------------------------------
# CLI-arg generation tests (no binary needed)
# ---------------------------------------------------------------------------

class TestCLIArgs:

    def test_ls_basic(self):
        c = cmd.LS(0.1, 10.0, 0.1)
        args = c._to_cli_args()
        assert args[0] == "-LS"
        assert "0.1" in args
        assert "10.0" in args

    def test_ls_save(self):
        c = cmd.LS(0.1, 10.0, 0.1, save_periodogram=True)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "1" in args
        assert "/tmp" in args

    def test_aov_basic(self):
        args = cmd.aov(0.1, 10.0, 0.1, finetune=2)._to_cli_args()
        assert args[0] == "-aov"

    def test_bls_basic(self):
        # nfreq= required when density_mode=False (the "optimal"
        # frequency grid is density-mode-only per vartools).
        args = cmd.BLS(0.5, 10.0, 1e-4, 0.01, 0.1, nfreq=1000)._to_cli_args()
        assert args[0] == "-BLS"

    def test_phase_basic(self):
        args = cmd.Phase(period=1.5)._to_cli_args()
        assert args[0] == "-Phase"
        assert "fix" in args
        assert "1.5" in args

    def test_clip_basic(self):
        args = cmd.clip(sigclip=3.0)._to_cli_args()
        assert args[0] == "-clip"
        assert "3.0" in args
        assert "1" in args  # iterative flag

    def test_rms_basic(self):
        args = cmd.rms()._to_cli_args()
        assert args[0] == "-rms"

    def test_stats_basic(self):
        args = cmd.stats("mag", "mean,stddev")._to_cli_args()
        assert args[0] == "-stats"

    def test_killharm_basic(self):
        args = cmd.Killharm(period="fix 1.5", nharm=3)._to_cli_args()
        assert args[0] == "-Killharm"
        assert "fix" in args
        assert "1.5" in args

    def test_expr_basic(self):
        args = cmd.expr("mag = mag - 0.01")._to_cli_args()
        assert args[0] == "-expr"

    def test_tfa_basic(self):
        # TFA requires trendlist, dates_file, pixelsep
        args = cmd.TFA("/tmp/trendlist.txt", "/tmp/dates.txt", 10.0)._to_cli_args()
        assert args[0] == "-TFA"
        assert "/tmp/trendlist.txt" in args

    def test_mandel_agol_basic(self):
        # MandelAgolTransit requires P0, T00
        args = cmd.MandelAgolTransit(P0=3.0, T00=1.0)._to_cli_args()
        assert args[0] == "-MandelAgolTransit"

    def test_addnoise_basic(self):
        args = cmd.addnoise(sig_white=0.005)._to_cli_args()
        assert args[0] == "-addnoise"

    def test_converttime_basic(self):
        args = cmd.converttime("hjd", "bjd")._to_cli_args()
        assert args[0] == "-converttime"
        assert "hjd" in args
        assert "bjd" in args

    def test_binlc_basic(self):
        args = cmd.binlc(binsize=0.01)._to_cli_args()
        assert args[0] == "-binlc"
        assert "binsize" in args
        assert "0.01" in args
        assert "average" in args
        assert "tcenter" in args

    def test_raw_string(self):
        args = cmd.Raw("-foo 1 2 3")._to_cli_args()
        assert args == ["-foo", "1", "2", "3"]

    def test_raw_list(self):
        args = cmd.Raw(["-bar", "x"])._to_cli_args()
        assert args == ["-bar", "x"]

    def test_sortlc(self):
        args = cmd.sortlc()._to_cli_args()
        assert args[0] == "-sortlc"

    def test_restricttimes(self):
        args = cmd.restricttimes(mode="JDrange", minJD=0.0, maxJD=10.0)._to_cli_args()
        assert args[0] == "-restricttimes"

    def test_decorr_basic(self):
        # decorr defaults should work
        args = cmd.decorr()._to_cli_args()
        assert args[0] == "-decorr"

    def test_fft_basic(self):
        # FFT requires input_real, input_imag, output_real, output_imag
        args = cmd.FFT("t", "mag", "freq", "amp")._to_cli_args()
        assert args[0] == "-FFT"


# ---------------------------------------------------------------------------
# End-to-end tests (require installed or in-tree binary)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools binary not found")
class TestEndToEnd:

    def test_rms_run(self):
        lc = make_sinusoidal_lc()
        result = vt.Pipeline([cmd.rms()]).run(lc)
        assert result.vars is not None
        # rms should report RMS or similar statistic
        assert len(result.vars) >= 1

    def test_clip_run(self):
        lc = make_sinusoidal_lc()
        result = vt.Pipeline([cmd.clip(sigclip=5.0)]).run(lc)
        assert result.vars is not None

    def test_ls_run(self):
        lc = make_sinusoidal_lc(n=300, period=1.5)
        result = vt.Pipeline([cmd.LS(0.5, 5.0, 1e-3)]).run(lc)
        assert result.vars is not None
        # The best LS period should be near 1.5 days
        period_keys = [k for k in result.vars.index if "Period" in str(k)]
        if period_keys:
            best = float(result.vars.loc[period_keys[0]])
            assert abs(best - 1.5) < 0.05, f"LS period {best} not near 1.5"

    def test_ls_save_output(self):
        lc = make_sinusoidal_lc(n=200, period=1.5)
        pipe = vt.Pipeline([cmd.LS(0.5, 5.0, 1e-3, save_periodogram=True)])
        result = pipe.run(lc)
        assert "LS_periodogram_0" in result.files, "LS periodogram output not captured"
        df = result.files["LS_periodogram_0"]
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0

    def test_stats_run(self):
        lc = make_sinusoidal_lc()
        result = vt.Pipeline([cmd.stats("mag", "mean,stddev")]).run(lc)
        assert result.vars is not None

    def test_clip_then_ls(self):
        lc = make_sinusoidal_lc(n=300, period=2.0)
        pipe = vt.Pipeline([
            cmd.clip(sigclip=5.0),
            cmd.LS(0.5, 5.0, 1e-3),
        ])
        result = pipe.run(lc)
        assert result.vars is not None

    def test_capture_lc(self):
        lc = make_sinusoidal_lc()
        result = vt.Pipeline([cmd.clip(sigclip=5.0)]).run(lc, capture_lc=True)
        assert result.lc is not None
        assert len(result.lc._df) > 0

    def test_batch_run(self):
        lcs = [make_sinusoidal_lc(period=1.0 + i * 0.3) for i in range(3)]
        result = vt.Pipeline([cmd.rms()]).run_batch(lcs)
        assert result.vars is not None
        assert len(result.vars) == 3

    def test_list_commands(self):
        cmds = vt.list_commands()
        assert isinstance(cmds, list)
        assert len(cmds) > 0
        assert any(c.upper() == "LS" for c in cmds)

    def test_addnoise_run(self):
        lc = make_sinusoidal_lc()
        result = vt.Pipeline([cmd.addnoise(sig_white=0.01)]).run(lc)
        assert result.vars is not None

    def test_sortlc_run(self):
        # Shuffle the light curve first, then sort
        lc = make_sinusoidal_lc()
        idx = np.random.permutation(len(lc._df))
        lc._df = lc._df.iloc[idx].reset_index(drop=True)
        result = vt.Pipeline([cmd.sortlc()]).run(lc, capture_lc=True)
        assert result.lc is not None
        times = result.lc._df.iloc[:, 0].values
        assert np.all(np.diff(times) >= 0), "sortlc did not sort times"
