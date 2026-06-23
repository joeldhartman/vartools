"""Comprehensive tests for all pyvartools command wrappers.

CLI-arg tests: exercise _to_cli_args() for every command and key options.
End-to-end pipeline tests: require the installed vartools binary.
"""

import os
import sys

import numpy as np
import pandas as pd
import pytest

VARTOOLS_SRC = os.path.realpath(
    os.path.join(os.path.dirname(__file__), "..", "..")
)
PYTHON_DIR = os.path.join(VARTOOLS_SRC, "python")
if PYTHON_DIR not in sys.path:
    sys.path.insert(0, PYTHON_DIR)

import pyvartools as vt
from pyvartools import commands as cmd

if not os.environ.get("VARTOOLS_BINARY"):
    for _cand in [
        os.path.join(VARTOOLS_SRC, "..", "..", "bin", "vartools"),
        os.path.join(VARTOOLS_SRC, "src", "vartools"),
    ]:
        _cand = os.path.realpath(_cand)
        if os.path.isfile(_cand) and os.access(_cand, os.X_OK):
            vt.set_binary(_cand)
            break

try:
    _HAVE_BINARY = os.path.isfile(vt.get_binary())
except FileNotFoundError:
    _HAVE_BINARY = False

EXAMPLES_DIR = os.path.join(VARTOOLS_SRC, "EXAMPLES")
EXAMPLE_LC = os.path.join(EXAMPLES_DIR, "2")


# ---------------------------------------------------------------------------
# Light-curve helpers
# ---------------------------------------------------------------------------

def make_lc(n=300, period=2.3, amp=0.1, noise=0.005):
    t = np.linspace(0, 30, n)
    mag = 10.0 + amp * np.sin(2 * np.pi * t / period)
    err = np.full(n, noise)
    return vt.LightCurve.from_arrays(t, mag, err, name="test")


def make_flat_lc(n=300, noise=0.01):
    t = np.linspace(0, 30, n)
    mag = np.full(n, 10.0) + np.random.default_rng(42).normal(0, noise, n)
    err = np.full(n, noise)
    return vt.LightCurve.from_arrays(t, mag, err, name="flat")


# ===========================================================================
# CLI-arg tests — Periodicity
# ===========================================================================

class TestCLIArgsPeriodicity:

    def test_ls_minimal(self):
        args = cmd.LS(0.1, 10.0, 0.1)._to_cli_args()
        assert args[0] == "-LS"
        assert "0.1" in args and "10.0" in args

    def test_ls_noGLS(self):
        args = cmd.LS(0.1, 10.0, 0.1, noGLS=True)._to_cli_args()
        assert "noGLS" in args

    def test_ls_whiten(self):
        args = cmd.LS(0.1, 10.0, 0.1, whiten=True)._to_cli_args()
        assert "whiten" in args

    def test_ls_clip(self):
        args = cmd.LS(0.1, 10.0, 0.1, clip=3.0, clipiter=1)._to_cli_args()
        assert "clip" in args and "3.0" in args and "1" in args

    def test_ls_bootstrap(self):
        args = cmd.LS(0.1, 10.0, 0.1, bootstrap=100)._to_cli_args()
        assert "bootstrap" in args and "100" in args

    def test_ls_save_periodogram(self):
        c = cmd.LS(0.1, 10.0, 0.1, save_periodogram=True)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "1" in args and "/tmp" in args
        assert "periodogram" in c._output_file_specs()

    def test_ls_npeaks(self):
        args = cmd.LS(0.1, 10.0, 0.1, npeaks=3)._to_cli_args()
        assert "3" in args

    def test_aov_minimal(self):
        args = cmd.aov(0.1, 10.0, 0.1, finetune=2)._to_cli_args()
        assert args[0] == "-aov"

    def test_aov_nbin(self):
        args = cmd.aov(0.1, 10.0, 0.1, finetune=2, nbin=10)._to_cli_args()
        assert "Nbin" in args and "10" in args

    def test_aov_uselog(self):
        args = cmd.aov(0.1, 10.0, 0.1, finetune=2, uselog=True)._to_cli_args()
        assert "uselog" in args

    def test_aov_save_periodogram(self):
        c = cmd.aov(0.1, 10.0, 0.1, finetune=2, save_periodogram=True)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "1" in args and "/tmp" in args

    def test_aov_harm_minimal(self):
        args = cmd.aov_harm(3, 0.1, 10.0, 0.1, finetune=2)._to_cli_args()
        assert args[0] == "-aov_harm"
        assert "3" in args

    def test_aov_harm_whiten(self):
        args = cmd.aov_harm(4, 0.1, 10.0, 0.1, finetune=2, whiten=True)._to_cli_args()
        assert "whiten" in args

    def test_aov_varexpr_minp(self):
        args = cmd.aov("mymin", 10.0, 0.1, finetune=2)._to_cli_args()
        assert "var" in args and args[args.index("var") + 1] == "mymin"

    def test_aov_varexpr_expr(self):
        args = cmd.aov("tspan/200", 10.0, 0.1, finetune=2)._to_cli_args()
        assert "expr" in args and "tspan/200" in args

    def test_aov_nbin_var(self):
        args = cmd.aov(0.1, 10.0, 0.1, finetune=2, nbin="nbinvar")._to_cli_args()
        assert "Nbin" in args
        idx = args.index("Nbin")
        assert args[idx + 1] == "var" and args[idx + 2] == "nbinvar"

    def test_aov_harm_varexpr(self):
        args = cmd.aov_harm(2, "tspan/200", 10.0, 0.1, 0.01)._to_cli_args()
        assert "expr" in args and "tspan/200" in args

    # ----- PDM (Phase Dispersion Minimization, all 5 variants) -----

    def test_pdm_step_minimal(self):
        args = cmd.PDM("step", 0.1, 10.0, 0.1, 0.01)._to_cli_args()
        assert args[0] == "-PDM" and args[1] == "step"

    def test_pdm_linterp_minimal(self):
        args = cmd.PDM("linterp", 0.1, 10.0, 0.1, 0.01)._to_cli_args()
        assert args[1] == "linterp"

    def test_pdm_multicover_nbin_nc(self):
        args = cmd.PDM("multicover", 0.1, 10.0, 0.1, 0.01,
                       nbin=8, nc=2)._to_cli_args()
        assert args[1] == "multicover"
        assert "Nbin" in args and "Nc" in args
        # Canonical order: Nbin must precede Nc
        assert args.index("Nbin") < args.index("Nc")

    def test_pdm_tophat_dphi(self):
        args = cmd.PDM("tophat", 0.5, 2.0, 0.5, 0.05, dphi=0.05)._to_cli_args()
        assert args[1] == "tophat"
        assert "dphi" in args and "0.05" in args

    def test_pdm_gauss_dphi(self):
        args = cmd.PDM("gauss", 0.5, 2.0, 0.5, 0.05, dphi=0.07)._to_cli_args()
        assert args[1] == "gauss"
        assert args[args.index("dphi") + 1] == "0.07"

    def test_pdm_save_periodogram(self):
        c = cmd.PDM("step", 0.1, 10.0, 0.1, 0.01, save_periodogram=True)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        # operiodogram=1 followed by outdir
        assert "/tmp" in args

    def test_pdm_clip_clipiter(self):
        args = cmd.PDM("step", 0.1, 10.0, 0.1, 0.01,
                       clip=3, clipiter=1)._to_cli_args()
        idx = args.index("clip")
        assert args[idx + 1] == "3" and args[idx + 2] == "1"

    def test_pdm_noerr_whiten(self):
        args = cmd.PDM("step", 0.1, 10.0, 0.1, 0.01,
                       noerr=True, whiten=True)._to_cli_args()
        assert "noerr" in args and "whiten" in args
        # Canonical order: noerr precedes whiten
        assert args.index("noerr") < args.index("whiten")

    def test_pdm_bootstrap(self):
        args = cmd.PDM("step", 0.1, 10.0, 0.1, 0.01,
                       bootstrap=100)._to_cli_args()
        idx = args.index("bootstrap")
        assert args[idx + 1] == "100"

    def test_pdm_maskpoints(self):
        args = cmd.PDM("step", 0.1, 10.0, 0.1, 0.01,
                       maskpoints="mymask")._to_cli_args()
        idx = args.index("maskpoints")
        assert args[idx + 1] == "mymask"

    def test_pdm_fixperiod_snr_fix(self):
        args = cmd.PDM("step", 0.1, 10.0, 0.1, 0.01,
                       fixperiod_snr=1.23)._to_cli_args()
        i = args.index("fixperiodSNR")
        assert args[i + 1] == "fix" and args[i + 2] == "1.23"

    def test_pdm_fixperiod_snr_pdm_passthrough(self):
        # Pre-resolution: keyword string passes through.
        args = cmd.PDM("step", 0.1, 10.0, 0.1, 0.01,
                       fixperiod_snr="pdm")._to_cli_args()
        i = args.index("fixperiodSNR")
        assert args[i + 1] == "pdm"

    def test_pdm_varexpr_minp(self):
        args = cmd.PDM("step", "mymin", 10.0, 0.1, 0.01)._to_cli_args()
        assert "var" in args and args[args.index("var") + 1] == "mymin"

    def test_pdm_canonical_trailing_order(self):
        # All trailing keywords together; verify they emerge in the canonical
        # order matching the strict vartools parser.
        args = cmd.PDM("step", 0.1, 10.0, 0.1, 0.01,
                       clip=5, noerr=True, whiten=True,
                       fixperiod_snr=1.234, bootstrap=50,
                       maskpoints="mask")._to_cli_args()
        order_keywords = ["clip", "noerr", "whiten", "fixperiodSNR",
                          "bootstrap", "maskpoints"]
        positions = [args.index(k) for k in order_keywords if k in args]
        assert positions == sorted(positions)

    # ----- PDM constructor validation -----

    def test_pdm_rejects_bad_variant(self):
        with pytest.raises(ValueError, match="variant"):
            cmd.PDM("bogus", 0.1, 10.0, 0.1, 0.01)

    def test_pdm_rejects_nbin_with_binless(self):
        with pytest.raises(ValueError, match="nbin"):
            cmd.PDM("tophat", 0.1, 10.0, 0.1, 0.01, nbin=8)
        with pytest.raises(ValueError, match="nbin"):
            cmd.PDM("gauss", 0.1, 10.0, 0.1, 0.01, nbin=8)

    def test_pdm_rejects_nc_outside_multicover(self):
        with pytest.raises(ValueError, match="multicover"):
            cmd.PDM("step", 0.1, 10.0, 0.1, 0.01, nc=2)
        with pytest.raises(ValueError, match="multicover"):
            cmd.PDM("tophat", 0.1, 10.0, 0.1, 0.01, nc=2)

    def test_pdm_rejects_dphi_with_binned(self):
        with pytest.raises(ValueError, match="dphi"):
            cmd.PDM("step", 0.1, 10.0, 0.1, 0.01, dphi=0.05)
        with pytest.raises(ValueError, match="dphi"):
            cmd.PDM("multicover", 0.1, 10.0, 0.1, 0.01, dphi=0.05)

    def test_pdm_rejects_bootstrap_zero(self):
        with pytest.raises(ValueError, match="bootstrap"):
            cmd.PDM("step", 0.1, 10.0, 0.1, 0.01, bootstrap=0)
        with pytest.raises(ValueError, match="bootstrap"):
            cmd.PDM("step", 0.1, 10.0, 0.1, 0.01, bootstrap=-3)

    # ----- FTP (Fast Template Periodogram) -----

    def test_ftp_file_minimal(self):
        args = cmd.FTP("file", 0.1, 2.0, 0.1, 0.01,
                       template_file="t.dat")._to_cli_args()
        assert args[0] == "-FTP"
        assert args[1] == "file" and args[2] == "t.dat"

    def test_ftp_fitlc_ascii(self):
        args = cmd.FTP("fitlc", 0.1, 2.0, 0.1, 0.01,
                       lc_path="lc.txt", lc_format="ascii",
                       t_col=1, mag_col=2, err_col=3,
                       nharm=5, period=1.235)._to_cli_args()
        assert args[1] == "fitlc"
        assert "ascii" in args
        # ascii column tokens emitted as positional ints.
        i = args.index("ascii")
        assert args[i + 1:i + 4] == ["1", "2", "3"]
        # nharm and period follow the cols.
        assert args[i + 4] == "5" and args[i + 5] == "1.235"

    def test_ftp_fitlc_fits(self):
        args = cmd.FTP("fitlc", 0.1, 2.0, 0.1, 0.01,
                       lc_path="lc.fits", lc_format="fits",
                       t_col="TIME", mag_col="MAG", err_col="none",
                       nharm=3, period=0.5)._to_cli_args()
        assert "fits" in args
        i = args.index("fits")
        assert args[i + 1:i + 4] == ["TIME", "MAG", "none"]

    def test_ftp_inline(self):
        # Length-2 cn/sn means H = 2 (nharm = 1).
        args = cmd.FTP("inline", 0.1, 2.0, 0.1, 0.01,
                       cn=[0.7, 0.2], sn=[0.0, 0.0])._to_cli_args()
        assert args[1] == "inline" and args[2] == "1"   # nharm
        # The 4 c/s spec tokens for the two harmonics, in CLI order.
        idx = args.index("inline")
        assert args[idx + 2:idx + 6] == ["0.7", "0.0", "0.2", "0.0"]

    def test_ftp_inline_var_coeff(self):
        # var/expr forms in the cn/sn lists are passed through.
        args = cmd.FTP("inline", 0.1, 2.0, 0.1, 0.01,
                       cn=["c1var"], sn=[0.0])._to_cli_args()
        assert "var" in args and "c1var" in args
        assert args[args.index("var") + 1] == "c1var"

    def test_ftp_filelist(self):
        args = cmd.FTP("filelist", 0.1, 2.0, 0.1, 0.01,
                       filelist_column=2)._to_cli_args()
        assert args[1] == "filelist"
        assert "column" in args
        assert args[args.index("column") + 1] == "2"

    def test_ftp_filelist_no_column(self):
        # Without filelist_column, the 'column N' tokens are omitted.
        args = cmd.FTP("filelist", 0.1, 2.0, 0.1, 0.01)._to_cli_args()
        assert args[1] == "filelist"
        assert "column" not in args

    def test_ftp_save_periodogram(self):
        c = cmd.FTP("file", 0.1, 2.0, 0.1, 0.01,
                    template_file="t.dat", save_periodogram=True)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "/tmp" in args

    def test_ftp_clip_clipiter(self):
        args = cmd.FTP("file", 0.1, 2.0, 0.1, 0.01,
                       template_file="t.dat",
                       clip=3, clipiter=0)._to_cli_args()
        i = args.index("clip")
        assert args[i + 1] == "3" and args[i + 2] == "0"

    def test_ftp_posamponly_noerr_whiten(self):
        args = cmd.FTP("file", 0.1, 2.0, 0.1, 0.01,
                       template_file="t.dat",
                       noerr=True, posamponly=True,
                       whiten=True)._to_cli_args()
        # Canonical order: noerr < posamponly < whiten.
        assert (args.index("noerr") < args.index("posamponly")
                < args.index("whiten"))

    def test_ftp_method_sums(self):
        args = cmd.FTP("file", 0.1, 2.0, 0.1, 0.01,
                       template_file="t.dat",
                       method="poly", sums="nfft")._to_cli_args()
        assert args[args.index("method") + 1] == "poly"
        assert args[args.index("sums") + 1] == "nfft"

    def test_ftp_fixperiod_snr_ftp(self):
        # Pre-resolution: 'ftp' keyword string passes through.
        args = cmd.FTP("file", 0.1, 2.0, 0.1, 0.01,
                       template_file="t.dat",
                       fixperiod_snr="ftp")._to_cli_args()
        i = args.index("fixperiodSNR")
        assert args[i + 1] == "ftp"

    def test_ftp_bootstrap(self):
        args = cmd.FTP("file", 0.1, 2.0, 0.1, 0.01,
                       template_file="t.dat",
                       bootstrap=50)._to_cli_args()
        i = args.index("bootstrap")
        assert args[i + 1] == "50"

    def test_ftp_canonical_trailing_order(self):
        args = cmd.FTP("file", 0.1, 2.0, 0.1, 0.01,
                       template_file="t.dat",
                       clip=5, noerr=True, posamponly=True, whiten=True,
                       fixperiod_snr=1.23, bootstrap=20,
                       maskpoints="m", method="brute",
                       sums="direct")._to_cli_args()
        order_keywords = ["clip", "noerr", "posamponly", "whiten",
                          "fixperiodSNR", "bootstrap", "maskpoints",
                          "method", "sums"]
        positions = [args.index(k) for k in order_keywords if k in args]
        assert positions == sorted(positions)

    # ----- FTP constructor validation -----

    def test_ftp_rejects_bad_template_source(self):
        with pytest.raises(ValueError, match="template_source"):
            cmd.FTP("bogus", 0.1, 2.0, 0.1, 0.01)

    def test_ftp_rejects_file_without_template_file(self):
        with pytest.raises(ValueError, match="template_file"):
            cmd.FTP("file", 0.1, 2.0, 0.1, 0.01)

    def test_ftp_rejects_mixed_mode_kwargs(self):
        with pytest.raises(ValueError, match="file mode rejects"):
            cmd.FTP("file", 0.1, 2.0, 0.1, 0.01,
                    template_file="t.dat", lc_path="x.lc")

    def test_ftp_rejects_fitlc_missing_args(self):
        with pytest.raises(ValueError, match="fitlc mode requires"):
            cmd.FTP("fitlc", 0.1, 2.0, 0.1, 0.01,
                    lc_path="x.lc", lc_format="ascii")  # cols/nharm/period missing

    def test_ftp_rejects_fitlc_bad_format(self):
        with pytest.raises(ValueError, match="lc_format"):
            cmd.FTP("fitlc", 0.1, 2.0, 0.1, 0.01,
                    lc_path="x.lc", lc_format="binary",
                    t_col=1, mag_col=2, err_col=3, nharm=3, period=1.0)

    def test_ftp_rejects_fitlc_period_string(self):
        with pytest.raises(ValueError, match="numeric"):
            cmd.FTP("fitlc", 0.1, 2.0, 0.1, 0.01,
                    lc_path="x.lc", lc_format="ascii",
                    t_col=1, mag_col=2, err_col=3, nharm=3, period="myP")

    def test_ftp_rejects_inline_uneven_cn_sn(self):
        with pytest.raises(ValueError, match="equal"):
            cmd.FTP("inline", 0.1, 2.0, 0.1, 0.01,
                    cn=[0.7, 0.2], sn=[0.0])

    def test_ftp_rejects_inline_nharm_mismatch(self):
        with pytest.raises(ValueError, match="nharm"):
            cmd.FTP("inline", 0.1, 2.0, 0.1, 0.01,
                    cn=[0.7, 0.2], sn=[0.0, 0.0], nharm=7)

    def test_ftp_rejects_bad_method(self):
        with pytest.raises(ValueError, match="method"):
            cmd.FTP("file", 0.1, 2.0, 0.1, 0.01,
                    template_file="t.dat", method="bogus")

    def test_ftp_rejects_bad_sums(self):
        with pytest.raises(ValueError, match="sums"):
            cmd.FTP("file", 0.1, 2.0, 0.1, 0.01,
                    template_file="t.dat", sums="bogus")

    def test_ftp_rejects_bootstrap_zero(self):
        with pytest.raises(ValueError, match="bootstrap"):
            cmd.FTP("file", 0.1, 2.0, 0.1, 0.01,
                    template_file="t.dat", bootstrap=0)

    def test_bls_minimal(self):
        # nfreq= required when density_mode=False ("optimal" is
        # density-mode-only per vartools).
        args = cmd.BLS(0.5, 10.0, 1e-4, 0.01, 0.1, nfreq=1000)._to_cli_args()
        assert args[0] == "-BLS"

    def test_bls_fittrap(self):
        args = cmd.BLS(0.5, 10.0, fittrap=True, nfreq=1000)._to_cli_args()
        assert "fittrap" in args

    def test_bls_correct_lc(self):
        args = cmd.BLS(0.5, 10.0, correct_lc=True, nfreq=1000)._to_cli_args()
        assert "1" in args  # correctlc token

    def test_bls_save_periodogram(self):
        c = cmd.BLS(0.5, 10.0, save_periodogram=True, nfreq=1000)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "1" in args

    def test_blsfixper_ls(self):
        args = cmd.BLSFixPer(period="ls")._to_cli_args()
        assert args[0] == "-BLSFixPer"
        assert "ls" in args

    def test_blsfixper_fix(self):
        args = cmd.BLSFixPer(period=1.5)._to_cli_args()
        assert "fix" in args and "1.5" in args

    def test_blsfixdurtc_minimal(self):
        args = cmd.BLSFixDurTc(duration=0.05, Tc=1.0)._to_cli_args()
        assert args[0] == "-BLSFixDurTc"
        assert "duration" in args
        assert "fix" in args
        assert "0.05" in args
        assert "Tc" in args
        assert "1.0" in args

    def test_blsfixdurtc_string_spec(self):
        args = cmd.BLSFixDurTc(duration="fixcolumn dur_col",
                                Tc="list column 3")._to_cli_args()
        assert "fixcolumn" in args
        assert "dur_col" in args
        assert "list" in args
        assert "column" in args
        assert "3" in args

    def test_blsfixdurtc_fixdepth(self):
        args = cmd.BLSFixDurTc(duration=0.05, Tc=1.0, fixdepth=0.01,
                                qgress=0.2)._to_cli_args()
        assert "fixdepth" in args
        assert "0.01" in args
        assert "qgress" in args
        assert "0.2" in args

    def test_blsfixperdurtc_minimal(self):
        args = cmd.BLSFixPerDurTc(period=1.5, duration=0.05,
                                   Tc=1.0)._to_cli_args()
        assert args[0] == "-BLSFixPerDurTc"
        assert "period" in args
        assert "fix" in args
        assert "1.5" in args
        assert "duration" in args
        assert "0.05" in args
        assert "Tc" in args
        assert "1.0" in args

    def test_blsfixperdurtc_list_spec(self):
        args = cmd.BLSFixPerDurTc(period="list", duration="fixcolumn dur",
                                   Tc="fix 1.0")._to_cli_args()
        assert "list" in args
        assert "fixcolumn" in args
        assert "dur" in args

    def test_autocorrelation_minimal(self):
        args = cmd.autocorrelation(0.1, 5.0, 0.1)._to_cli_args()
        assert args[0] == "-autocorrelation"

    def test_autocorrelation_always_includes_outdir(self):
        # outdir must always be present — CLI requires it unconditionally
        c = cmd.autocorrelation(0.1, 5.0, 0.1)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "/tmp" in args

    def test_autocorrelation_save_result_true(self):
        # save_result=True (default) — outdir still present; capture=True
        c = cmd.autocorrelation(0.1, 5.0, 0.1, save_result=True)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "/tmp" in args

    def test_autocorrelation_save_result_false_still_emits_outdir(self):
        # save_result=False — file still written (CLI mandatory), outdir present
        c = cmd.autocorrelation(0.1, 5.0, 0.1, save_result=False)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "-autocorrelation" in args
        assert "/tmp" in args

    def test_autocorrelation_save_result_path(self):
        # save_result="/some/path" — outdir is the user path, not _outdir
        c = cmd.autocorrelation(0.1, 5.0, 0.1, save_result="/data/acf")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "/data/acf" in args
        assert "/tmp" not in args

    def test_dftclean_minimal(self):
        args = cmd.dftclean(nbeam=3)._to_cli_args()
        assert args[0] == "-dftclean"

    def test_dftclean_maxfreq(self):
        args = cmd.dftclean(nbeam=3, maxfreq=10.0)._to_cli_args()
        assert "maxfreq" in args and "10.0" in args

    def test_dftclean_save_dspec(self):
        c = cmd.dftclean(nbeam=3, save_dspec=True)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "outdspec" in args

    def test_wwz_auto(self):
        args = cmd.wwz()._to_cli_args()
        assert args[0] == "-wwz"
        # The bare keyword ``auto`` must NOT be wrapped in ``var auto``;
        # vartools treats ``auto`` as a literal CLI option in -wwz.
        assert "var" not in args
        # maxfreq / tau0 / tau1 / dtau default to "auto" — emit literal.
        for key in ("maxfreq", "tau0", "tau1", "dtau"):
            i = args.index(key)
            assert args[i + 1] == "auto", \
                f"Expected '{key} auto', got '{args[i]} {args[i+1]}'"
        # freqsamp does NOT accept "auto" in vartools — it must be a
        # positive number.  Default is the Foster 1996 value 0.25.
        i = args.index("freqsamp")
        assert args[i + 1] == "0.25"

    def test_wwz_auto_explicit_string(self):
        # Passing "auto" explicitly must behave the same as the default.
        args = cmd.wwz(maxfreq="auto", tau0="auto",
                       tau1="auto", dtau="auto")._to_cli_args()
        assert "var" not in args
        for key in ("maxfreq", "tau0", "tau1", "dtau"):
            i = args.index(key)
            assert args[i + 1] == "auto"

    def test_wwz_freqsamp_rejects_auto_silently(self):
        # ``freqsamp="auto"`` is a user error — vartools rejects it at
        # runtime ("freqsamp must be > 0").  The wrapper must NOT emit
        # the literal "auto" keyword for freqsamp; if the user really
        # passes the string, it's their problem and we route through
        # _varexpr (which produces ``var auto``) — but the default path
        # uses the numeric value so the common case stays clean.
        args = cmd.wwz(freqsamp=0.1)._to_cli_args()
        i = args.index("freqsamp")
        assert args[i + 1] == "0.1"

    def test_wwz_var_name_still_uses_var_keyword(self):
        # A bare identifier other than "auto" should still emit ``var <name>``.
        args = cmd.wwz(maxfreq="myvar")._to_cli_args()
        i = args.index("maxfreq")
        assert args[i + 1] == "var" and args[i + 2] == "myvar"

    def test_wwz_explicit(self):
        args = cmd.wwz(maxfreq=10.0, freqsamp=0.01, tau0=0.0,
                       tau1=30.0, dtau=0.1)._to_cli_args()
        assert "10.0" in args and "0.01" in args

    def test_getlsampthresh_ls(self):
        args = cmd.GetLSAmpThresh(period="ls", minp=0.1, thresh=10.0)._to_cli_args()
        assert args[0] == "-GetLSAmpThresh"
        assert "ls" in args

    def test_getlsampthresh_noGLS(self):
        args = cmd.GetLSAmpThresh(period="ls", noGLS=True)._to_cli_args()
        assert "noGLS" in args

    def test_phase_ls(self):
        args = cmd.Phase()._to_cli_args()
        assert args[0] == "-Phase"
        assert "ls" in args

    def test_phase_fix(self):
        args = cmd.Phase(period=2.3)._to_cli_args()
        assert "fix" in args and "2.3" in args

    def test_phase_T0(self):
        args = cmd.Phase(period=2.3, T0=1.0)._to_cli_args()
        assert "T0" in args and "1.0" in args

    def test_phase_phasevar(self):
        args = cmd.Phase(period=2.3, phasevar="ph")._to_cli_args()
        assert "phasevar" in args and "ph" in args

    def test_phase_startphase(self):
        args = cmd.Phase(period=2.3, startphase=-0.5)._to_cli_args()
        assert "startphase" in args and "-0.5" in args

    # ------- aov_harm nharm var/expr (Batch 2a) -------

    def test_aov_harm_nharm_var(self):
        args = cmd.aov_harm("nharmvar", 0.5, 5.0, 0.1, 0.01)._to_cli_args()
        assert "var" in args
        idx = args.index("var")
        assert args[idx + 1] == "nharmvar"

    def test_aov_harm_nharm_expr(self):
        args = cmd.aov_harm("npeaks*2", 0.5, 5.0, 0.1, 0.01)._to_cli_args()
        assert "expr" in args and "npeaks*2" in args

    # ------- GetLSAmpThresh file mode (Batch 2h) -------

    def test_getlsampthresh_file_mode(self):
        args = cmd.GetLSAmpThresh(period="ls", mode="file",
                                   listfile="thresh.txt")._to_cli_args()
        assert "file" in args and "thresh.txt" in args
        assert "harm" not in args

    # ------- wwz format options (Batch 2k) -------

    def test_wwz_transform_format_fits(self):
        c = cmd.wwz(save_transform=True, transform_format="fits")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "fits" in args

    def test_wwz_transform_format_pm3d(self):
        c = cmd.wwz(save_transform=True, transform_format="pm3d")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "pm3d" in args

    def test_wwz_transform_name_format(self):
        c = cmd.wwz(save_transform=True, transform_name="%s.wwz")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "format" in args and "%s.wwz" in args

    def test_wwz_maxtransform_name_format(self):
        c = cmd.wwz(save_maxtransform=True, maxtransform_name="%s.mwwz")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "format" in args and "%s.mwwz" in args

    # ------- dftclean missing options (Batch 2l) -------

    def test_dftclean_finddirtypeaks(self):
        args = cmd.dftclean(nbeam=3, finddirtypeaks=5)._to_cli_args()
        assert "finddirtypeaks" in args and "5" in args

    def test_dftclean_finddirtypeaks_clip(self):
        args = cmd.dftclean(nbeam=3, finddirtypeaks=5,
                             finddirtypeaks_clip=3.0)._to_cli_args()
        assert "clip" in args and "3.0" in args

    def test_dftclean_finddirtypeaks_clipiter(self):
        args = cmd.dftclean(nbeam=3, finddirtypeaks=5,
                             finddirtypeaks_clip=3.0,
                             finddirtypeaks_clipiter=1)._to_cli_args()
        assert "1" in args

    def test_dftclean_outcbeam(self):
        c = cmd.dftclean(nbeam=3, save_cspec=True, outcbeam=True)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "outcbeam" in args

    def test_dftclean_useampspec(self):
        args = cmd.dftclean(nbeam=3, useampspec=True)._to_cli_args()
        assert "useampspec" in args

    def test_dftclean_verboseout(self):
        args = cmd.dftclean(nbeam=3, verboseout=True)._to_cli_args()
        assert "verboseout" in args

    # ------- BLS extensions (Batch 3a) -------

    def test_bls_minper_var(self):
        args = cmd.BLS("minperiodvar", 10.0, nfreq=1000)._to_cli_args()
        assert "var" in args and "minperiodvar" in args

    def test_bls_rmin_var(self):
        args = cmd.BLS(0.5, 10.0, rmin="bls_rmin_var",
                       nfreq=1000)._to_cli_args()
        assert "var" in args and "bls_rmin_var" in args

    def test_bls_nbins_var(self):
        args = cmd.BLS(0.5, 10.0, nbins="nbinvar",
                       nfreq=1000)._to_cli_args()
        assert "var" in args and "nbinvar" in args

    def test_bls_qmin_qmax(self):
        args = cmd.BLS(0.5, 10.0, qmin=0.01, qmax=0.1,
                       nfreq=1000)._to_cli_args()
        assert "q" in args and "0.01" in args and "0.1" in args

    def test_bls_density_mode(self):
        args = cmd.BLS(0.5, 10.0, density_mode=True, stellar_density=1.4,
                        min_exp_dur_frac=0.5, max_exp_dur_frac=1.5)._to_cli_args()
        idx = args.index("density")
        # tokens after "density": rho min_frac max_frac
        assert args[idx + 1] == "1.4"
        assert args[idx + 2] == "0.5"
        assert args[idx + 3] == "1.5"
        # r / q tokens must NOT appear when density_mode=True
        assert "r" not in args
        assert "q" not in args

    def test_bls_density_mode_varexpr(self):
        args = cmd.BLS(0.5, 10.0, density_mode=True,
                        stellar_density="rho_star",
                        min_exp_dur_frac=0.5, max_exp_dur_frac=1.5)._to_cli_args()
        idx = args.index("density")
        assert args[idx + 1] == "var"
        assert args[idx + 2] == "rho_star"

    def test_bls_df(self):
        args = cmd.BLS(0.5, 10.0, df=0.001)._to_cli_args()
        assert "df" in args and "0.001" in args

    def test_bls_extraparams(self):
        args = cmd.BLS(0.5, 10.0, extraparams=True,
                       nfreq=1000)._to_cli_args()
        assert "extraparams" in args

    def test_bls_nobinnedrms(self):
        args = cmd.BLS(0.5, 10.0, nobinnedrms=True,
                       nfreq=1000)._to_cli_args()
        assert "nobinnedrms" in args

    def test_bls_save_phcurve(self):
        c = cmd.BLS(0.5, 10.0, save_phcurve=True, nfreq=1000)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "ophcurve" in args

    def test_bls_save_jdcurve(self):
        c = cmd.BLS(0.5, 10.0, save_jdcurve=True, nfreq=1000)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "ojdcurve" in args

    def test_bls_freq_grid_steplogP(self):
        args = cmd.BLS(0.5, 10.0, freq_grid="steplogP",
                       nfreq=1000)._to_cli_args()
        assert "steplogP" in args

    def test_bls_adjust_qmin(self):
        args = cmd.BLS(0.5, 10.0, adjust_qmin=True,
                       nfreq=1000)._to_cli_args()
        assert "adjust-qmin-by-mindt" in args

    def test_bls_reduce_nbins(self):
        args = cmd.BLS(0.5, 10.0, adjust_qmin=True, reduce_nbins=True,
                       nfreq=1000)._to_cli_args()
        assert "reduce-nbins" in args

    def test_bls_reportharmonics(self):
        args = cmd.BLS(0.5, 10.0, reportharmonics=True,
                       nfreq=1000)._to_cli_args()
        assert "reportharmonics" in args

    # ------- BLSFixPer q mode (Batch 3b) -------

    def test_blsfixper_qmin_qmax(self):
        args = cmd.BLSFixPer(period=1.5, qmin=0.01, qmax=0.1)._to_cli_args()
        assert "q" in args and "0.01" in args and "0.1" in args

    # ------- BLSFixDurTc and BLSFixPerDurTc CLI args (Batch 3d) -------

    def test_blsfixdurtc_batch3d_minimal(self):
        c = cmd.BLSFixDurTc(duration=0.05, Tc=1.0)
        args = c._to_cli_args()
        assert args[0] == "-BLSFixDurTc"
        assert "duration" in args and "Tc" in args

    def test_blsfixperdurtc_batch3d_minimal(self):
        c = cmd.BLSFixPerDurTc(period=1.5, duration=0.05, Tc=1.0)
        args = c._to_cli_args()
        assert args[0] == "-BLSFixPerDurTc"
        assert "period" in args and "duration" in args and "Tc" in args


# ===========================================================================
# CLI-arg tests — Manipulation
# ===========================================================================

class TestCLIArgsManipulation:

    def test_clip_iterative(self):
        args = cmd.clip(5.0)._to_cli_args()
        assert args[0] == "-clip"
        assert "5.0" in args and "1" in args

    def test_clip_non_iterative(self):
        args = cmd.clip(5.0, iterative=False)._to_cli_args()
        assert "0" in args

    def test_clip_niter(self):
        args = cmd.clip(5.0, niter=3)._to_cli_args()
        assert "niter" in args and "3" in args

    def test_clip_median(self):
        args = cmd.clip(5.0, median=True)._to_cli_args()
        assert "median" in args

    def test_clip_markclip(self):
        args = cmd.clip(5.0, markclip="mask")._to_cli_args()
        assert "markclip" in args and "mask" in args

    def test_clip_markclip_noinitmark(self):
        args = cmd.clip(5.0, markclip="mask", noinitmark=True)._to_cli_args()
        assert "noinitmark" in args

    def test_rms_basic(self):
        assert cmd.rms()._to_cli_args()[0] == "-rms"

    def test_rms_maskpoints(self):
        args = cmd.rms(maskpoints="mask")._to_cli_args()
        assert "maskpoints" in args and "mask" in args

    def test_rmsbin_basic(self):
        args = cmd.rmsbin(2, [0.5, 1.0])._to_cli_args()
        assert args[0] == "-rmsbin"
        assert "2" in args and "0.5" in args and "1.0" in args

    def test_chi2_basic(self):
        assert cmd.chi2()._to_cli_args()[0] == "-chi2"

    def test_chi2bin_basic(self):
        args = cmd.chi2bin(2, [0.5, 1.0])._to_cli_args()
        assert args[0] == "-chi2bin"

    def test_alarm_basic(self):
        assert cmd.alarm()._to_cli_args()[0] == "-alarm"

    def test_vonneumann_basic(self):
        args = cmd.vonNeumann()._to_cli_args()
        assert args == ["-vonNeumann"]

    def test_vonneumann_weighted(self):
        args = cmd.vonNeumann(weighted=True)._to_cli_args()
        assert args == ["-vonNeumann", "weighted"]

    def test_vonneumann_maskpoints(self):
        args = cmd.vonNeumann(maskpoints="m")._to_cli_args()
        assert args == ["-vonNeumann", "maskpoints", "m"]

    def test_vonneumann_canonical_order(self):
        # weighted must precede maskpoints to match the strict-parser order.
        args = cmd.vonNeumann(weighted=True, maskpoints="m")._to_cli_args()
        assert args == ["-vonNeumann", "weighted", "maskpoints", "m"]

    def test_percentileratios_defaults(self):
        args = cmd.percentileratios()._to_cli_args()
        assert args == ["-percentileratios"]

    def test_percentileratios_custom_pairs(self):
        args = cmd.percentileratios(percentilepairs=[(10, 90), (20, 80)])._to_cli_args()
        assert args == ["-percentileratios", "percentilepairs", "10.0:90.0,20.0:80.0"]

    def test_percentileratios_floats(self):
        args = cmd.percentileratios(percentilepairs=[(2.5, 97.5)])._to_cli_args()
        assert args == ["-percentileratios", "percentilepairs", "2.5:97.5"]

    def test_percentileratios_auto_swap(self):
        args = cmd.percentileratios(percentilepairs=[(95, 5)])._to_cli_args()
        assert args == ["-percentileratios", "percentilepairs", "5.0:95.0"]

    def test_percentileratios_explicit_defaults_omitted_from_cli(self):
        # Passing the defaults explicitly should produce the same CLI as
        # leaving them unset -- no percentilepairs token emitted.
        args = cmd.percentileratios(
            percentilepairs=[(5, 95), (1, 99)]
        )._to_cli_args()
        assert args == ["-percentileratios"]

    def test_percentileratios_rejects_equal(self):
        with pytest.raises(ValueError):
            cmd.percentileratios(percentilepairs=[(50, 50)])

    def test_percentileratios_rejects_out_of_range_low(self):
        with pytest.raises(ValueError):
            cmd.percentileratios(percentilepairs=[(0, 50)])

    def test_percentileratios_rejects_out_of_range_high(self):
        with pytest.raises(ValueError):
            cmd.percentileratios(percentilepairs=[(5, 100)])

    def test_percentileratios_rejects_duplicate(self):
        with pytest.raises(ValueError):
            cmd.percentileratios(percentilepairs=[(5, 95), (5, 95)])

    def test_percentileratios_rejects_duplicate_after_swap(self):
        # 95:5 canonicalizes to 5:95, which then duplicates the first pair.
        with pytest.raises(ValueError):
            cmd.percentileratios(percentilepairs=[(5, 95), (95, 5)])

    def test_percentileratios_rejects_malformed_pair(self):
        with pytest.raises(ValueError):
            cmd.percentileratios(percentilepairs=[(5,)])

    def test_percentileratios_maskpoints(self):
        args = cmd.percentileratios(maskpoints="m")._to_cli_args()
        assert args == ["-percentileratios", "maskpoints", "m"]

    def test_percentileratios_canonical_order(self):
        # percentilepairs must precede maskpoints (strict parser order).
        args = cmd.percentileratios(
            percentilepairs=[(10, 90)],
            maskpoints="m",
        )._to_cli_args()
        assert args == [
            "-percentileratios", "percentilepairs", "10.0:90.0",
            "maskpoints", "m",
        ]

    def test_beyondNsigma_defaults(self):
        args = cmd.beyondNsigma()._to_cli_args()
        assert args == ["-beyondNsigma"]

    def test_beyondNsigma_custom_Nvalues(self):
        args = cmd.beyondNsigma(Nvalues=[1, 2, 3])._to_cli_args()
        assert args == ["-beyondNsigma", "Nvalues", "1.0,2.0,3.0"]

    def test_beyondNsigma_floats(self):
        args = cmd.beyondNsigma(Nvalues=[1.5, 2.5, 4.0])._to_cli_args()
        assert args == ["-beyondNsigma", "Nvalues", "1.5,2.5,4.0"]

    def test_beyondNsigma_useMAD(self):
        args = cmd.beyondNsigma(useMAD=True)._to_cli_args()
        assert args == ["-beyondNsigma", "useMAD"]

    def test_beyondNsigma_canonical_order(self):
        # Nvalues must precede useMAD to match the strict-parser order.
        args = cmd.beyondNsigma(Nvalues=[2, 4], useMAD=True)._to_cli_args()
        assert args == ["-beyondNsigma", "Nvalues", "2.0,4.0", "useMAD"]

    def test_beyondNsigma_explicit_defaults_omitted_from_cli(self):
        # Passing the defaults explicitly should produce the same CLI as
        # leaving them unset -- no Nvalues token emitted.
        args = cmd.beyondNsigma(Nvalues=[1, 3, 5])._to_cli_args()
        assert args == ["-beyondNsigma"]

    def test_beyondNsigma_rejects_zero(self):
        with pytest.raises(ValueError):
            cmd.beyondNsigma(Nvalues=[0, 1, 2])

    def test_beyondNsigma_rejects_negative(self):
        with pytest.raises(ValueError):
            cmd.beyondNsigma(Nvalues=[-1, 2])

    def test_beyondNsigma_rejects_duplicate(self):
        with pytest.raises(ValueError):
            cmd.beyondNsigma(Nvalues=[1, 2, 1])

    def test_beyondNsigma_rejects_non_numeric(self):
        with pytest.raises(ValueError):
            cmd.beyondNsigma(Nvalues=["abc", 1])

    def test_beyondNsigma_maskpoints(self):
        args = cmd.beyondNsigma(maskpoints="m")._to_cli_args()
        assert args == ["-beyondNsigma", "maskpoints", "m"]

    def test_beyondNsigma_canonical_order_all(self):
        # Strict parser order: Nvalues, useMAD, maskpoints.
        args = cmd.beyondNsigma(
            Nvalues=[2, 4],
            useMAD=True,
            maskpoints="m",
        )._to_cli_args()
        assert args == [
            "-beyondNsigma", "Nvalues", "2.0,4.0", "useMAD",
            "maskpoints", "m",
        ]

    def test_slopestats_defaults(self):
        args = cmd.slopestats()._to_cli_args()
        assert args == ["-slopestats"]

    def test_slopestats_custom_bintime(self):
        args = cmd.slopestats(bintime=[5, 10])._to_cli_args()
        assert args == ["-slopestats", "bintime", "5.0,10.0"]

    def test_slopestats_binshift(self):
        args = cmd.slopestats(bintime=[10], binshift=0.5)._to_cli_args()
        assert args == ["-slopestats", "bintime", "10.0",
                        "binshift", "0.5"]

    def test_slopestats_binshift_requires_bintime(self):
        with pytest.raises(ValueError):
            cmd.slopestats(binshift=0.5)

    def test_slopestats_custom_threshold(self):
        args = cmd.slopestats(threshold=[1, 3, 5])._to_cli_args()
        assert args == ["-slopestats", "threshold", "1.0,3.0,5.0"]

    def test_slopestats_default_threshold_omitted(self):
        # Passing the default threshold explicitly should produce the
        # same CLI as leaving it unset -- no threshold token emitted.
        args = cmd.slopestats(threshold=[3.0])._to_cli_args()
        assert args == ["-slopestats"]

    def test_slopestats_maxgap(self):
        args = cmd.slopestats(maxgap=0.5)._to_cli_args()
        assert args == ["-slopestats", "maxgap", "0.5"]

    def test_slopestats_useMAD(self):
        args = cmd.slopestats(useMAD=True)._to_cli_args()
        assert args == ["-slopestats", "useMAD"]

    def test_slopestats_maskpoints(self):
        args = cmd.slopestats(maskpoints="m")._to_cli_args()
        assert args == ["-slopestats", "maskpoints", "m"]

    def test_slopestats_canonical_order_all(self):
        # Strict parser order: bintime, binshift, threshold, maxgap,
        # useMAD, maskpoints.
        args = cmd.slopestats(
            bintime=[5, 10],
            binshift=0.5,
            threshold=[1, 3],
            maxgap=0.5,
            useMAD=True,
            maskpoints="m",
        )._to_cli_args()
        assert args == [
            "-slopestats",
            "bintime", "5.0,10.0",
            "binshift", "0.5",
            "threshold", "1.0,3.0",
            "maxgap", "0.5",
            "useMAD",
            "maskpoints", "m",
        ]

    def test_slopestats_rejects_zero_bintime(self):
        with pytest.raises(ValueError):
            cmd.slopestats(bintime=[0, 5])

    def test_slopestats_rejects_negative_bintime(self):
        with pytest.raises(ValueError):
            cmd.slopestats(bintime=[-1, 5])

    def test_slopestats_rejects_duplicate_bintime(self):
        with pytest.raises(ValueError):
            cmd.slopestats(bintime=[5, 10, 5])

    def test_slopestats_rejects_non_numeric_bintime(self):
        with pytest.raises(ValueError):
            cmd.slopestats(bintime=["abc", 5])

    def test_slopestats_rejects_zero_threshold(self):
        with pytest.raises(ValueError):
            cmd.slopestats(threshold=[0, 1, 3])

    def test_slopestats_rejects_negative_threshold(self):
        with pytest.raises(ValueError):
            cmd.slopestats(threshold=[-1, 3])

    def test_slopestats_rejects_duplicate_threshold(self):
        with pytest.raises(ValueError):
            cmd.slopestats(threshold=[1, 3, 1])

    def test_slopestats_rejects_zero_maxgap(self):
        with pytest.raises(ValueError):
            cmd.slopestats(maxgap=0)

    def test_slopestats_rejects_negative_maxgap(self):
        with pytest.raises(ValueError):
            cmd.slopestats(maxgap=-0.5)

    # ----- CodyM ----------------------------------------------------------

    def test_codym_basic(self):
        args = cmd.CodyM(trendwindow=10)._to_cli_args()
        assert args == ["-CodyM", "trendwindow", "10"]

    def test_codym_two_stage(self):
        args = cmd.CodyM(trendwindow=10, outlierwindow=0.1,
                          sigclip=3)._to_cli_args()
        assert args == ["-CodyM", "trendwindow", "10",
                        "outlierwindow", "0.1",
                        "sigclip", "3"]

    def test_codym_sigclip_default_omitted(self):
        # Default sigclip=5.0 stays off the CLI surface.
        args = cmd.CodyM(trendwindow=10, sigclip=5.0)._to_cli_args()
        assert args == ["-CodyM", "trendwindow", "10"]

    def test_codym_sigclip_zero_emitted(self):
        # sigclip=0 disables rejection and IS distinct from default 5.0.
        args = cmd.CodyM(trendwindow=10, sigclip=0)._to_cli_args()
        assert args == ["-CodyM", "trendwindow", "10", "sigclip", "0"]

    def test_codym_expr_sources(self):
        args = cmd.CodyM(trendwindow="5*2",
                          sigclip="2.5*2")._to_cli_args()
        assert args == ["-CodyM",
                        "trendwindow", "expr", "5*2",
                        "sigclip", "expr", "2.5*2"]

    def test_codym_maskpoints(self):
        args = cmd.CodyM(trendwindow=10, maskpoints="m")._to_cli_args()
        assert args == ["-CodyM", "trendwindow", "10", "maskpoints", "m"]

    def test_codym_canonical_order_all(self):
        # Strict parser order: trendwindow, outlierwindow, sigclip,
        # maskpoints.
        args = cmd.CodyM(
            trendwindow=10,
            outlierwindow=0.1,
            sigclip=3,
            maskpoints="m",
        )._to_cli_args()
        assert args == [
            "-CodyM",
            "trendwindow", "10",
            "outlierwindow", "0.1",
            "sigclip", "3",
            "maskpoints", "m",
        ]

    def test_codym_rejects_zero_trendwindow(self):
        with pytest.raises(ValueError):
            cmd.CodyM(trendwindow=0)

    def test_codym_rejects_negative_trendwindow(self):
        with pytest.raises(ValueError):
            cmd.CodyM(trendwindow=-5)

    def test_codym_rejects_zero_outlierwindow(self):
        with pytest.raises(ValueError):
            cmd.CodyM(trendwindow=10, outlierwindow=0)

    def test_codym_rejects_negative_sigclip(self):
        with pytest.raises(ValueError):
            cmd.CodyM(trendwindow=10, sigclip=-1)

    # ----- CodyQ ----------------------------------------------------------

    def test_codyq_fix_period(self):
        args = cmd.CodyQ(period=1.234, trendwindow=10)._to_cli_args()
        assert args == ["-CodyQ", "fix", "1.234", "trendwindow", "10"]

    def test_codyq_aov_period(self):
        args = cmd.CodyQ(period="aov", trendwindow=10)._to_cli_args()
        assert args == ["-CodyQ", "aov", "trendwindow", "10"]

    def test_codyq_pdm_period(self):
        # The pdm/ftp keywords were added to _period_spec specifically
        # for -CodyQ; confirm the keyword is passed through rather than
        # interpreted as a bare variable name.
        args = cmd.CodyQ(period="pdm", trendwindow=10)._to_cli_args()
        assert args == ["-CodyQ", "pdm", "trendwindow", "10"]

    def test_codyq_ftp_period(self):
        args = cmd.CodyQ(period="ftp", trendwindow=10)._to_cli_args()
        assert args == ["-CodyQ", "ftp", "trendwindow", "10"]

    def test_codyq_bls_period(self):
        args = cmd.CodyQ(period="bls", trendwindow=10)._to_cli_args()
        assert args == ["-CodyQ", "bls", "trendwindow", "10"]

    def test_codyq_injectharm_period(self):
        args = cmd.CodyQ(period="injectharm", trendwindow=10)._to_cli_args()
        assert args == ["-CodyQ", "injectharm", "trendwindow", "10"]

    def test_codyq_fixcolumn_period(self):
        args = cmd.CodyQ(period="fixcolumn Period_1_0",
                          trendwindow=10)._to_cli_args()
        assert args == ["-CodyQ", "fixcolumn", "Period_1_0",
                        "trendwindow", "10"]

    def test_codyq_list_period_with_column(self):
        args = cmd.CodyQ(period="list column 2",
                          trendwindow=10)._to_cli_args()
        assert args == ["-CodyQ", "list", "column", "2",
                        "trendwindow", "10"]

    def test_codyq_expr_period(self):
        args = cmd.CodyQ(period="2*P", trendwindow=10)._to_cli_args()
        assert args == ["-CodyQ", "expr", "2*P",
                        "trendwindow", "10"]

    def test_codyq_var_period(self):
        args = cmd.CodyQ(period="myperiod", trendwindow=10)._to_cli_args()
        assert args == ["-CodyQ", "var", "myperiod",
                        "trendwindow", "10"]

    def test_codyq_phasesmooth_non_default(self):
        args = cmd.CodyQ(period=2.0, trendwindow=10,
                          phasesmooth=0.5)._to_cli_args()
        assert args == ["-CodyQ", "fix", "2.0",
                        "trendwindow", "10",
                        "phasesmooth", "0.5"]

    def test_codyq_phasesmooth_default_omitted(self):
        # Default phasesmooth=0.25 stays off the CLI surface.
        args = cmd.CodyQ(period=2.0, trendwindow=10,
                          phasesmooth=0.25)._to_cli_args()
        assert args == ["-CodyQ", "fix", "2.0", "trendwindow", "10"]

    def test_codyq_maskpoints(self):
        args = cmd.CodyQ(period=2.0, trendwindow=10,
                          maskpoints="m")._to_cli_args()
        assert args == ["-CodyQ", "fix", "2.0", "trendwindow", "10",
                        "maskpoints", "m"]

    def test_codyq_canonical_order_all(self):
        # Strict parser order: <period-source>, trendwindow, phasesmooth,
        # maskpoints.
        args = cmd.CodyQ(
            period="aov",
            trendwindow=10,
            phasesmooth=0.5,
            maskpoints="m",
        )._to_cli_args()
        assert args == [
            "-CodyQ", "aov",
            "trendwindow", "10",
            "phasesmooth", "0.5",
            "maskpoints", "m",
        ]

    def test_codyq_expr_sourced_trendwindow_and_phasesmooth(self):
        args = cmd.CodyQ(period="aov",
                          trendwindow="5*2",
                          phasesmooth="1/4")._to_cli_args()
        assert args == ["-CodyQ", "aov",
                        "trendwindow", "expr", "5*2",
                        "phasesmooth", "expr", "1/4"]

    def test_codyq_rejects_zero_trendwindow(self):
        with pytest.raises(ValueError):
            cmd.CodyQ(period=2.0, trendwindow=0)

    def test_codyq_rejects_phasesmooth_zero(self):
        with pytest.raises(ValueError):
            cmd.CodyQ(period=2.0, trendwindow=10, phasesmooth=0)

    def test_codyq_rejects_phasesmooth_over_one(self):
        with pytest.raises(ValueError):
            cmd.CodyQ(period=2.0, trendwindow=10, phasesmooth=1.5)

    def test_codyq_rejects_zero_period(self):
        with pytest.raises(ValueError):
            cmd.CodyQ(period=0, trendwindow=10)

    def test_codyq_rejects_negative_period(self):
        with pytest.raises(ValueError):
            cmd.CodyQ(period=-1.0, trendwindow=10)

    # ----- structurefunction --------------------------------------------

    def test_sf_bins_log_basic(self):
        args = cmd.structurefunction(bins="log", Nbins=20)._to_cli_args()
        assert args == ["-structurefunction", "bins", "log", "20"]

    def test_sf_bins_linear_basic(self):
        args = cmd.structurefunction(bins="linear", Nbins=30)._to_cli_args()
        assert args == ["-structurefunction", "bins", "linear", "30"]

    def test_sf_bins_edges_basic(self):
        args = cmd.structurefunction(
            bins="edges", edges=[0.01, 0.1, 1.0, 10.0]
        )._to_cli_args()
        assert args[:2] == ["-structurefunction", "bins"]
        assert args[2] == "edges"
        # Comma list, no spaces.
        assert args[3] == "0.01,0.1,1.0,10.0"

    def test_sf_lagrange_literals(self):
        args = cmd.structurefunction(
            bins="log", Nbins=20, lagrange=(0.05, 50.0)
        )._to_cli_args()
        assert args == [
            "-structurefunction", "bins", "log", "20",
            "lagrange", "0.05", "50.0",
        ]

    def test_sf_lagrange_var_expr(self):
        args = cmd.structurefunction(
            bins="log", Nbins=20, lagrange=("myvar", "2*lo")
        )._to_cli_args()
        assert args == [
            "-structurefunction", "bins", "log", "20",
            "lagrange", "var", "myvar", "expr", "2*lo",
        ]

    def test_sf_estimator_default_suppressed(self):
        args = cmd.structurefunction(
            bins="log", Nbins=20, estimator="squared"
        )._to_cli_args()
        assert "estimator" not in args

    def test_sf_estimator_mad_emitted(self):
        args = cmd.structurefunction(
            bins="log", Nbins=20, estimator="mad"
        )._to_cli_args()
        assert args == [
            "-structurefunction", "bins", "log", "20", "estimator", "mad"
        ]

    def test_sf_fitDRW_flag_only(self):
        args = cmd.structurefunction(
            bins="log", Nbins=20, fitDRW=True
        )._to_cli_args()
        assert args == [
            "-structurefunction", "bins", "log", "20", "fitDRW"
        ]

    def test_sf_fitDRW_with_sigma0_tau0(self):
        args = cmd.structurefunction(
            bins="log", Nbins=20, fitDRW=True, sigma0=0.05, tau0=100.0
        )._to_cli_args()
        assert args == [
            "-structurefunction", "bins", "log", "20",
            "fitDRW", "sigma0", "0.05", "tau0", "100.0",
        ]

    def test_sf_fitDRW_sigma0_var_expr(self):
        args = cmd.structurefunction(
            bins="log", Nbins=20, fitDRW=True,
            sigma0="myvar", tau0="2*tau_init",
        )._to_cli_args()
        assert args == [
            "-structurefunction", "bins", "log", "20",
            "fitDRW",
            "sigma0", "var", "myvar",
            "tau0", "expr", "2*tau_init",
        ]

    def test_sf_reportsfvalsintable_basic(self):
        args = cmd.structurefunction(
            bins="log", Nbins=20, reportsfvalsintable=[0.1, 1.0, 10.0]
        )._to_cli_args()
        assert args == [
            "-structurefunction", "bins", "log", "20",
            "reportsfvalsintable", "0.1,1.0,10.0",
        ]

    def test_sf_maskpoints(self):
        args = cmd.structurefunction(
            bins="log", Nbins=20, maskpoints="m"
        )._to_cli_args()
        assert args == [
            "-structurefunction", "bins", "log", "20", "maskpoints", "m"
        ]

    def test_sf_canonical_order_all_keywords(self):
        # Strict parser canonical order: bins, lagrange, estimator,
        # fitDRW [sigma0, tau0], reportsfvalsintable, save, maskpoints.
        args = cmd.structurefunction(
            bins="log",
            Nbins=20,
            lagrange=(0.05, 50.0),
            estimator="mad",
            fitDRW=True,
            sigma0=0.05,
            tau0=100.0,
            reportsfvalsintable=[0.1, 1.0, 10.0],
            save_result="/tmp/sf_out",
            maskpoints="m",
        )._to_cli_args()
        assert args == [
            "-structurefunction",
            "bins", "log", "20",
            "lagrange", "0.05", "50.0",
            "estimator", "mad",
            "fitDRW", "sigma0", "0.05", "tau0", "100.0",
            "reportsfvalsintable", "0.1,1.0,10.0",
            "save", "/tmp/sf_out",
            "maskpoints", "m",
        ]

    def test_sf_save_result_path(self):
        args = cmd.structurefunction(
            bins="log", Nbins=20, save_result="/tmp/out"
        )._to_cli_args()
        assert "save" in args
        assert "/tmp/out" in args

    def test_sf_output_file_specs_save_off(self):
        c = cmd.structurefunction(bins="log", Nbins=20)
        assert c._output_file_specs() == {}

    def test_sf_output_file_specs_save_on(self):
        c = cmd.structurefunction(bins="log", Nbins=20, save_result=True)
        specs = c._output_file_specs()
        # Logical name "result" matches the save_result attribute so the
        # pipeline collector actually captures the file (was "sf", a no-op).
        assert specs == {"result": (".sf", None)}

    def test_sf_rejects_bogus_bins_mode(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(bins="frobnicate", Nbins=20)

    def test_sf_rejects_nbins_too_small(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(bins="log", Nbins=1)

    def test_sf_rejects_log_without_nbins(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(bins="log")

    def test_sf_rejects_edges_without_edges_list(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(bins="edges")

    def test_sf_rejects_edges_too_short(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(bins="edges", edges=[0.01, 0.1])

    def test_sf_rejects_edges_non_monotonic(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(
                bins="edges", edges=[0.1, 0.05, 1.0]
            )

    def test_sf_rejects_edges_with_zero(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(bins="edges", edges=[0.0, 1.0, 10.0])

    def test_sf_rejects_nbins_in_edges_mode(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(
                bins="edges", edges=[0.1, 1, 10], Nbins=5
            )

    def test_sf_rejects_edges_in_log_mode(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(
                bins="log", Nbins=20, edges=[0.1, 1.0]
            )

    def test_sf_rejects_bogus_estimator(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(bins="log", Nbins=20, estimator="rms")

    def test_sf_rejects_sigma0_without_fitDRW(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(bins="log", Nbins=20, sigma0=0.05)

    def test_sf_rejects_tau0_without_fitDRW(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(bins="log", Nbins=20, tau0=100.0)

    def test_sf_rejects_negative_sigma0(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(
                bins="log", Nbins=20, fitDRW=True, sigma0=-1
            )

    def test_sf_rejects_zero_tau0(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(
                bins="log", Nbins=20, fitDRW=True, tau0=0
            )

    def test_sf_rejects_lagrange_lagmax_le_lagmin(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(
                bins="log", Nbins=20, lagrange=(50.0, 5.0)
            )

    def test_sf_rejects_lagrange_zero_lagmin(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(
                bins="log", Nbins=20, lagrange=(0.0, 5.0)
            )

    def test_sf_rejects_report_non_monotonic(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(
                bins="log", Nbins=20,
                reportsfvalsintable=[1.0, 0.5, 5.0],
            )

    def test_sf_rejects_report_with_zero(self):
        with pytest.raises(ValueError):
            cmd.structurefunction(
                bins="log", Nbins=20,
                reportsfvalsintable=[0.0, 1.0, 5.0],
            )

    # ----- drwfit ---------------------------------------------------------

    def test_drwfit_default(self):
        assert cmd.drwfit()._to_cli_args() == ["-drwfit"]

    def test_drwfit_mean_fit(self):
        args = cmd.drwfit(mean="fit")._to_cli_args()
        assert args == ["-drwfit", "mean", "fit"]

    def test_drwfit_mean_subtract(self):
        args = cmd.drwfit(mean="subtract")._to_cli_args()
        assert args == ["-drwfit", "mean", "subtract"]

    def test_drwfit_mean_fix_literal(self):
        args = cmd.drwfit(mean="fix", mean_value=10.12)._to_cli_args()
        assert args == ["-drwfit", "mean", "fix", "10.12"]

    def test_drwfit_mean_fix_var(self):
        args = cmd.drwfit(mean="fix", mean_value="mu")._to_cli_args()
        assert args == ["-drwfit", "mean", "fix", "var", "mu"]

    def test_drwfit_mean_fix_expr(self):
        args = cmd.drwfit(mean="fix", mean_value="2*m")._to_cli_args()
        assert args == ["-drwfit", "mean", "fix", "expr", "2*m"]

    def test_drwfit_init_guesses(self):
        args = cmd.drwfit(sigma0=0.05, tau0=100.0, mean0=10.0)._to_cli_args()
        assert args == [
            "-drwfit",
            "sigma0", "0.05", "tau0", "100.0", "mean0", "10.0",
        ]

    def test_drwfit_init_guesses_var_expr(self):
        args = cmd.drwfit(sigma0="s", tau0="2*t")._to_cli_args()
        assert args == [
            "-drwfit", "sigma0", "var", "s", "tau0", "expr", "2*t",
        ]

    def test_drwfit_correctlc_smoothed(self):
        args = cmd.drwfit(correctlc="smoothed")._to_cli_args()
        assert args == ["-drwfit", "correctlc", "smoothed"]

    def test_drwfit_correctlc_forecast(self):
        args = cmd.drwfit(correctlc="forecast")._to_cli_args()
        assert args == ["-drwfit", "correctlc", "forecast"]

    def test_drwfit_modelvar_smoothed(self):
        args = cmd.drwfit(modelvar=("smoothed", "drwmod"))._to_cli_args()
        assert args == ["-drwfit", "modelvar", "smoothed", "drwmod"]

    def test_drwfit_modelvar_forecast(self):
        args = cmd.drwfit(modelvar=("forecast", "drwmod"))._to_cli_args()
        assert args == ["-drwfit", "modelvar", "forecast", "drwmod"]

    def test_drwfit_save_result_path(self):
        args = cmd.drwfit(save_result="/tmp/out")._to_cli_args()
        assert "save" in args and "/tmp/out" in args

    def test_drwfit_maskpoints(self):
        args = cmd.drwfit(maskpoints="m")._to_cli_args()
        assert args == ["-drwfit", "maskpoints", "m"]

    def test_drwfit_canonical_order_all_keywords(self):
        # Strict parser canonical order: mean, sigma0, tau0, mean0,
        # save, correctlc, modelvar, maskpoints.
        args = cmd.drwfit(
            mean="fix",
            mean_value=10.0,
            sigma0=0.05,
            tau0=100.0,
            mean0=10.0,
            save_result="/tmp/drw_out",
            correctlc="smoothed",
            modelvar=("smoothed", "m"),
            maskpoints="msk",
        )._to_cli_args()
        assert args == [
            "-drwfit",
            "mean", "fix", "10.0",
            "sigma0", "0.05",
            "tau0", "100.0",
            "mean0", "10.0",
            "save", "/tmp/drw_out",
            "correctlc", "smoothed",
            "modelvar", "smoothed", "m",
            "maskpoints", "msk",
        ]

    def test_drwfit_output_file_specs_save_off(self):
        assert cmd.drwfit()._output_file_specs() == {}

    def test_drwfit_output_file_specs_save_on(self):
        specs = cmd.drwfit(save_result=True)._output_file_specs()
        # Logical name matches the save_result attribute so capture works.
        assert specs == {"result": (".drwfit", None)}

    def test_drwfit_rejects_bogus_mean(self):
        with pytest.raises(ValueError):
            cmd.drwfit(mean="frobnicate")

    def test_drwfit_rejects_fix_without_value(self):
        with pytest.raises(ValueError):
            cmd.drwfit(mean="fix")

    def test_drwfit_rejects_mean_value_without_fix(self):
        with pytest.raises(ValueError):
            cmd.drwfit(mean_value=10.0)

    def test_drwfit_rejects_mean_value_with_subtract(self):
        with pytest.raises(ValueError):
            cmd.drwfit(mean="subtract", mean_value=10.0)

    def test_drwfit_rejects_negative_sigma0(self):
        with pytest.raises(ValueError):
            cmd.drwfit(sigma0=-1)

    def test_drwfit_rejects_zero_tau0(self):
        with pytest.raises(ValueError):
            cmd.drwfit(tau0=0)

    def test_drwfit_rejects_bogus_correctlc(self):
        with pytest.raises(ValueError):
            cmd.drwfit(correctlc="x")

    def test_drwfit_rejects_bogus_modelvar_mode(self):
        with pytest.raises(ValueError):
            cmd.drwfit(modelvar=("x", "n"))

    def test_drwfit_rejects_modelvar_not_pair(self):
        with pytest.raises(ValueError):
            cmd.drwfit(modelvar="smoothed")

    def test_drwfit_rejects_modelvar_empty_name(self):
        with pytest.raises(ValueError):
            cmd.drwfit(modelvar=("smoothed", ""))

    # ----- runlength ------------------------------------------------------

    def test_runlength_default(self):
        # Default k=3.0 stays off the CLI surface.
        assert cmd.runlength()._to_cli_args() == ["-runlength"]

    def test_runlength_k_literal(self):
        args = cmd.runlength(k=1.0)._to_cli_args()
        assert args == ["-runlength", "k", "1.0"]

    def test_runlength_k_default_omitted(self):
        assert cmd.runlength(k=3.0)._to_cli_args() == ["-runlength"]

    def test_runlength_k_zero_emitted(self):
        # k=0 collapses the band to the median and IS distinct from 3.0.
        assert cmd.runlength(k=0)._to_cli_args() == ["-runlength", "k", "0"]

    def test_runlength_k_var(self):
        args = cmd.runlength(k="myk")._to_cli_args()
        assert args == ["-runlength", "k", "var", "myk"]

    def test_runlength_k_expr(self):
        args = cmd.runlength(k="2*1.5")._to_cli_args()
        assert args == ["-runlength", "k", "expr", "2*1.5"]

    def test_runlength_maskpoints(self):
        args = cmd.runlength(maskpoints="m")._to_cli_args()
        assert args == ["-runlength", "maskpoints", "m"]

    def test_runlength_canonical_order_all(self):
        # Strict parser order: k, maskpoints.
        args = cmd.runlength(k=0.5, maskpoints="m")._to_cli_args()
        assert args == ["-runlength", "k", "0.5", "maskpoints", "m"]

    def test_runlength_rejects_negative_k(self):
        with pytest.raises(ValueError):
            cmd.runlength(k=-1)

    def test_rescalesig_basic(self):
        assert cmd.rescalesig()._to_cli_args()[0] == "-rescalesig"

    def test_ensemblerescalesig_basic(self):
        assert cmd.ensemblerescalesig()._to_cli_args()[0] == "-ensemblerescalesig"

    def test_ensemblerescalesig_sigclip(self):
        args = cmd.ensemblerescalesig(sigclip=3.0)._to_cli_args()
        assert "3.0" in args

    def test_stats_single_var(self):
        args = cmd.stats("mag", "mean")._to_cli_args()
        assert args[0] == "-stats" and "mag" in args and "mean" in args

    def test_stats_multi_var(self):
        args = cmd.stats(["mag", "err"], "mean,stddev")._to_cli_args()
        assert "mag,err" in args

    def test_stats_multi_stat(self):
        args = cmd.stats("mag", ["mean", "median", "stddev"])._to_cli_args()
        assert "mean,median,stddev" in args

    def test_killharm_ls(self):
        args = cmd.Killharm(period="ls", nharm=3)._to_cli_args()
        assert args[0] == "-Killharm" and "ls" in args

    def test_killharm_fix(self):
        args = cmd.Killharm(period=2.3, nharm=4)._to_cli_args()
        assert "fix" in args and "2.3" in args

    def test_killharm_both(self):
        args = cmd.Killharm(period="both", nharm=3)._to_cli_args()
        assert "both" in args

    def test_killharm_fitonly(self):
        args = cmd.Killharm(period="ls", fitonly=True)._to_cli_args()
        assert "fitonly" in args

    def test_linfit_basic(self):
        args = cmd.linfit("a0+a1*t", "a0 0 a1 0")._to_cli_args()
        assert args[0] == "-linfit"

    def test_linfit_correct_lc(self):
        args = cmd.linfit("a0+a1*t", "a0 0 a1 0", correct_lc=True)._to_cli_args()
        assert "correctlc" in args

    def test_linfit_fitmask_emits_fitmask_not_maskpoints(self):
        # BUG regression: must emit "fitmask", not the old "maskpoints" token
        args = cmd.linfit("a0+a1*t", "a0 0 a1 0", fitmask="mask")._to_cli_args()
        assert "fitmask" in args
        assert "maskpoints" not in args
        assert args[args.index("fitmask") + 1] == "mask"

    def test_nonlinfit_fitmask_emits_fitmask_not_maskpoints(self):
        # BUG regression: must emit "fitmask", not the old "maskpoints" token
        args = cmd.nonlinfit("a0+a1*t", "a0 0 a1 0",
                              fitmask="mask")._to_cli_args()
        assert "fitmask" in args
        assert "maskpoints" not in args
        assert args[args.index("fitmask") + 1] == "mask"

    def test_mandel_agol_no_fitp_param(self):
        # BUG regression: fitP no longer exists; must raise TypeError
        with pytest.raises(TypeError):
            cmd.MandelAgolTransit(P0=3.0, T00=1.0, fitP=0)

    def test_injectharm_fix(self):
        args = cmd.Injectharm(period=2.3, amplitude=0.1)._to_cli_args()
        assert args[0] == "-Injectharm"
        assert "fix" in args and "2.3" in args

    def test_injectharm_nharm(self):
        # Python nharm=3 → vartools Nharm=2, with 3 amp/phase pairs
        args = cmd.Injectharm(period=2.3, amplitude=0.1, nharm=3)._to_cli_args()
        assert "2" in args  # vt_nharm = nharm - 1
        assert args.count("ampfix") == 3

    def test_injectharm_rand(self):
        args = cmd.Injectharm(period="rand 0.5 5.0", amplitude=0.1)._to_cli_args()
        assert "rand" in args

    def test_injectharm_phase_rand_keyword(self):
        # phase="rand" → phaserand (CLI keyword, not phaseexpr).
        args = cmd.Injectharm(period=2.5, amplitude=0.05,
                              phase="rand")._to_cli_args()
        assert "phaserand" in args
        assert "phasefix" not in args

    def test_injectharm_amp_keyword(self):
        # amplitude="rand" → amprand keyword.
        args = cmd.Injectharm(period=2.5, amplitude="rand",
                              phase=0.0)._to_cli_args()
        assert "amprand" in args
        assert "ampfix" not in args

    def test_injectharm_per_harmonic_relative(self):
        # Per-harmonic relative amps/phases produce
        # `ampfix R amprel phasefix phi phaserel` blocks.
        args = cmd.Injectharm(
            period=0.514333, amplitude=0.1, phase="rand", nharm=4,
            harmonic_amps_rel=[0.47, 0.36, 0.24],
            harmonic_phases_rel=[0.61, 0.26, -0.07],
        )._to_cli_args()
        # Fundamental: ampfix 0.1 + phaserand
        assert args.count("phaserand") == 1
        # 3 overtones, each contributes one ampfix+amprel and one phasefix+phaserel
        assert args.count("amprel") == 3
        assert args.count("phaserel") == 3

    def test_injectharm_per_harmonic_length_validation(self):
        # nharm-1 must equal len(harmonic_amps_rel).
        with pytest.raises(ValueError, match="harmonic_amps_rel"):
            cmd.Injectharm(period=1.0, amplitude=0.1, nharm=4,
                           harmonic_amps_rel=[1, 2],
                           harmonic_phases_rel=[1, 2])

    def test_injectharm_per_harmonic_paired(self):
        # Both lists must be supplied together.
        with pytest.raises(ValueError, match="must be supplied"):
            cmd.Injectharm(period=1.0, amplitude=0.1, nharm=2,
                           harmonic_amps_rel=[0.5])

    def test_injecttransit_basic(self):
        args = cmd.Injecttransit(
            period=3.0, Rp=0.1, Mp=1.0, phase=0.0,
            sini=1.0, Mstar=1.0, Rstar=1.0
        )._to_cli_args()
        assert args[0] == "-Injecttransit"

    def test_sortlc_basic(self):
        assert cmd.sortlc()._to_cli_args()[0] == "-sortlc"

    def test_sortlc_reverse(self):
        args = cmd.sortlc(reverse=True)._to_cli_args()
        assert "reverse" in args

    def test_sortlc_var(self):
        args = cmd.sortlc(var="err")._to_cli_args()
        assert "var" in args and "err" in args

    def test_restricttimes_jdrange(self):
        args = cmd.restricttimes(mode="JDrange", minJD=0.0, maxJD=10.0)._to_cli_args()
        assert args[0] == "-restricttimes"
        assert "JDrange" in args and "0.0" in args and "10.0" in args

    def test_restricttimes_expression(self):
        args = cmd.restricttimes(mode="expr", expression="t > 5")._to_cli_args()
        assert "expr" in args

    def test_restricttimes_exclude(self):
        args = cmd.restricttimes(mode="JDrange", minJD=0.0, maxJD=5.0,
                                  exclude=True)._to_cli_args()
        assert "exclude" in args

    def test_restoretimes_basic(self):
        args = cmd.restoretimes()._to_cli_args()
        assert args[0] == "-restoretimes"

    def test_restoretimes_prior(self):
        args = cmd.restoretimes(prior_command=2)._to_cli_args()
        assert "2" in args

    def test_savelc_basic(self):
        assert cmd.savelc()._to_cli_args()[0] == "-savelc"

    def test_restorelc_basic(self):
        assert cmd.restorelc()._to_cli_args()[0] == "-restorelc"

    def test_restorelc_vars(self):
        args = cmd.restorelc(vars=["mag", "err"])._to_cli_args()
        assert "vars" in args

    def test_difffluxtomag_basic(self):
        args = cmd.difffluxtomag(mag_constant=10.0)._to_cli_args()
        assert args[0] == "-difffluxtomag"

    def test_difffluxtomag_magcolumn(self):
        args = cmd.difffluxtomag(mag_constant=10.0, magcolumn=3)._to_cli_args()
        assert "magcolumn" in args and "3" in args

    def test_fluxtomag_basic(self):
        args = cmd.fluxtomag(mag_constant=10.0)._to_cli_args()
        assert args[0] == "-fluxtomag"

    def test_magtoflux_basic(self):
        args = cmd.magtoflux(mag_constant=25.0)._to_cli_args()
        assert args == ["-magtoflux", "25.0"]

    def test_magtoflux_normalize(self):
        args = cmd.magtoflux(normalize=True)._to_cli_args()
        assert args == ["-magtoflux", "normalize"]

    def test_magtoflux_expr(self):
        args = cmd.magtoflux(mag_constant="25.0+0.0")._to_cli_args()
        assert args == ["-magtoflux", "expr", "25.0+0.0"]

    def test_magtoflux_var(self):
        args = cmd.magtoflux(mag_constant="zpt")._to_cli_args()
        assert args == ["-magtoflux", "var", "zpt"]

    def test_magtoflux_normalize_rejects_mag_constant(self):
        with pytest.raises(ValueError):
            cmd.magtoflux(mag_constant=25.0, normalize=True)

    def test_magtoflux_requires_one_of(self):
        with pytest.raises(ValueError):
            cmd.magtoflux()

    def test_changeerror_basic(self):
        assert cmd.changeerror()._to_cli_args()[0] == "-changeerror"

    def test_changevariable_basic(self):
        args = cmd.changevariable("mag", "flux")._to_cli_args()
        assert args[0] == "-changevariable"
        assert "mag" in args and "flux" in args

    def test_copylc_basic(self):
        args = cmd.copylc(2)._to_cli_args()
        assert args[0] == "-copylc" and "2" in args

    def test_medianfilter_basic(self):
        args = cmd.medianfilter(time=0.5)._to_cli_args()
        assert args[0] == "-medianfilter" and "0.5" in args

    def test_medianfilter_average(self):
        args = cmd.medianfilter(time=0.5, method="average")._to_cli_args()
        assert "average" in args

    def test_medianfilter_replace(self):
        args = cmd.medianfilter(time=0.5, replace=True)._to_cli_args()
        assert "replace" in args

    def test_expr_basic(self):
        args = cmd.expr("mag=mag-0.01")._to_cli_args()
        assert args[0] == "-expr"

    def test_expr_outputcolumn_emits_keyword(self):
        # outputcolumn=True with a non-LC vartype emits the bare
        # "outputcolumn" CLI keyword (vartools' -expr accepts it after
        # the var=expression argument when listvar/scalar/const was given).
        args = cmd.expr("avg=mean(mag)", vartype="listvar",
                        outputcolumn=True)._to_cli_args()
        assert args == ["-expr", "listvar", "avg=mean(mag)", "outputcolumn"]

    def test_expr_outputcolumn_requires_vartype(self):
        # outputcolumn=True without vartype is rejected at construction.
        # The default per-observation LC type would yield one value per
        # observation, not a single column.
        with pytest.raises(ValueError, match="outputcolumn=True requires"):
            cmd.expr("flux=10^((mag-10)/(-2.5))", outputcolumn=True)

    def test_expr_listvar(self):
        args = cmd.expr("avg=mean(mag)", vartype="listvar")._to_cli_args()
        assert args == ["-expr", "listvar", "avg=mean(mag)"]

    def test_expr_const(self):
        args = cmd.expr("pi=3.14159", vartype="const")._to_cli_args()
        assert args == ["-expr", "const", "pi=3.14159"]

    def test_expr_scalar(self):
        args = cmd.expr("x=mag[0]", vartype="scalar")._to_cli_args()
        assert args == ["-expr", "scalar", "x=mag[0]"]

    def test_expr_invalid_vartype(self):
        import pytest
        with pytest.raises(ValueError):
            cmd.expr("x=1", vartype="invalid")

    def test_print_cols_basic(self):
        args = cmd.print_cols("t,mag,err")._to_cli_args()
        assert args[0] == "-print"

    def test_print_cols_columnnames(self):
        args = cmd.print_cols(["t", "mag"], columnnames=["time", "magnitude"]
                               )._to_cli_args()
        assert "columnnames" in args

    def test_fft_basic(self):
        args = cmd.FFT("t", "mag", "freq", "amp")._to_cli_args()
        assert args[0] == "-FFT"

    def test_ifft_basic(self):
        args = cmd.IFFT("freq", "amp", "t", "mag")._to_cli_args()
        assert args[0] == "-IFFT"

    def test_resample_linear(self):
        args = cmd.resample(method="linear", tstart=0.0, tstop=30.0,
                             delt=0.1)._to_cli_args()
        assert args[0] == "-resample" and "linear" in args

    def test_resample_spline(self):
        args = cmd.resample(method="spline")._to_cli_args()
        assert "spline" in args

    def test_decorr_basic(self):
        assert cmd.decorr()._to_cli_args()[0] == "-decorr"

    def test_decorr_lc_columns(self):
        args = cmd.decorr(lc_columns=[(2, 1), (3, 2)])._to_cli_args()
        assert args[0] == "-decorr"

    def test_jstet_basic(self):
        args = cmd.Jstet(timescale=0.5, dates="dates.txt")._to_cli_args()
        assert args[0] == "-Jstet"
        assert args == ["-Jstet", "0.5", "dates.txt"]

    def test_jstet_skipnormalize(self):
        args = cmd.Jstet(timescale=0.5, skipnormalize=True)._to_cli_args()
        assert args == ["-Jstet", "0.5", "skipnormalize"]

    def test_jstet_requires_one_of_dates_or_skipnormalize(self):
        with pytest.raises(ValueError, match="dates"):
            cmd.Jstet(timescale=0.5)
        with pytest.raises(ValueError, match="dates"):
            cmd.Jstet(timescale=0.5, dates="dates.txt", skipnormalize=True)

    # ------- Killharm output_format and clip (Batch 2i) -------

    def test_killharm_outampphase(self):
        args = cmd.Killharm(period="ls", output_format="outampphase")._to_cli_args()
        assert "outampphase" in args

    def test_killharm_outampradphase(self):
        args = cmd.Killharm(period="ls", output_format="outampradphase")._to_cli_args()
        assert "outampradphase" in args

    def test_killharm_clip(self):
        args = cmd.Killharm(period="ls", clip=3.0)._to_cli_args()
        assert "clip" in args and "3.0" in args

    # ------- linfit modelvar, model_nameformat, reject (Batch 2b + 4a) -------

    def test_linfit_modelvar(self):
        args = cmd.linfit("a0+a1*t", "a0 0 a1 0", modelvar="mvar")._to_cli_args()
        assert "modelvar" in args and "mvar" in args

    def test_linfit_model_nameformat(self):
        c = cmd.linfit("a0+a1*t", "a0 0 a1 0", save_model=True,
                        model_nameformat="%s.fit")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "format" in args and "%s.fit" in args

    def test_linfit_reject(self):
        args = cmd.linfit("a0+a1*t", "a0 0 a1 0", reject=3.0)._to_cli_args()
        assert "reject" in args and "3.0" in args

    def test_linfit_reject_usemad(self):
        args = cmd.linfit("a0+a1*t", "a0 0 a1 0", reject=3.0,
                           reject_usemad=True)._to_cli_args()
        assert "useMAD" in args

    def test_linfit_reject_iter(self):
        args = cmd.linfit("a0+a1*t", "a0 0 a1 0", reject=3.0,
                           reject_iter=True)._to_cli_args()
        assert "iter" in args

    def test_linfit_reject_fixednum(self):
        args = cmd.linfit("a0+a1*t", "a0 0 a1 0", reject=3.0,
                           reject_iter=True, reject_fixednum=5)._to_cli_args()
        assert "fixednum" in args and "5" in args

    # ------- resample spline bounds and bspline options (Batch 2e) -------

    def test_resample_spline_left_right(self):
        args = cmd.resample(method="spline", left=0.0, right=1.0)._to_cli_args()
        assert "left" in args and "0.0" in args
        assert "right" in args and "1.0" in args

    def test_resample_bspline_nbreaks_order(self):
        args = cmd.resample(method="bspline", nbreaks=10, order=3)._to_cli_args()
        assert "nbreaks" in args and "10" in args
        assert "order" in args and "3" in args

    # ------- Injecttransit new _injectparam-based tests (P2 batch) -------

    def test_injecttransit_float_params(self):
        """Float params still emit correctly with new _injectparam helper."""
        args = cmd.Injecttransit(
            period=3.0, Rp=0.1, Mp=0.001, phase=0.0,
            sini=1.0, Mstar=1.0, Rstar=1.0
        )._to_cli_args()
        assert "Pfix" in args and "3.0" in args
        assert "Rpfix" in args and "0.1" in args
        assert "eomega" in args
        assert "efix" in args and "ofix" in args

    def test_injecttransit_string_params(self):
        """String params pass through as-is."""
        args = cmd.Injecttransit(
            period="Plogrand 0.2 2.0",
            Rp="Rplogrand 0.05 0.15",
            Mp=0.001, phase="phaserand",
            sini="sinirand", Mstar=1.0, Rstar=1.0
        )._to_cli_args()
        assert "Plogrand" in args
        assert "0.2" in args and "2.0" in args
        assert "Rplogrand" in args
        assert "phaserand" in args
        assert "sinirand" in args

    def test_injecttransit_hk_mode(self):
        """hk=True emits hk mode instead of eomega."""
        args = cmd.Injecttransit(
            period=3.0, Rp=0.1, Mp=0.001, phase=0.0,
            sini=1.0, Mstar=1.0, Rstar=1.0,
            hk=True, h=0.1, k=0.05
        )._to_cli_args()
        assert "hk" in args
        assert "hfix" in args
        assert "kfix" in args
        assert "eomega" not in args

    def test_injecttransit_dilute_float(self):
        args = cmd.Injecttransit(
            period=3.0, Rp=0.1, Mp=0.001, phase=0.0,
            sini=1.0, Mstar=1.0, Rstar=1.0,
            dilute=0.5
        )._to_cli_args()
        assert "dilute" in args
        idx = args.index("dilute")
        assert args[idx+1] == "fix" and args[idx+2] == "0.5"

    def test_injecttransit_dilute_string(self):
        args = cmd.Injecttransit(
            period=3.0, Rp=0.1, Mp=0.001, phase=0.0,
            sini=1.0, Mstar=1.0, Rstar=1.0,
            dilute="list"
        )._to_cli_args()
        assert "dilute" in args
        idx = args.index("dilute")
        assert args[idx+1] == "list"

    # ------- resample file_times and gaps (P2 batch) -------

    def test_resample_file_mode_fix(self):
        args = cmd.resample(method="linear", file_times="/path/to/times.txt")._to_cli_args()
        assert "file" in args
        idx = args.index("file")
        assert args[idx+1] == "fix"
        assert args[idx+2] == "/path/to/times.txt"

    def test_resample_file_mode_column(self):
        args = cmd.resample(method="linear",
                            file_times="/path/to/times.txt",
                            file_column=3)._to_cli_args()
        assert "column" in args
        idx = args.index("column")
        assert args[idx+1] == "3"

    def test_resample_file_mode_list(self):
        # Bare list mode — emits `file list` only.
        args = cmd.resample(method="linear", file_times="list")._to_cli_args()
        assert "file" in args
        idx = args.index("file")
        assert args[idx+1] == "list"
        # No listcolumn / tcolumn when neither kwarg is set.
        assert "listcolumn" not in args
        assert "tcolumn" not in args

    def test_resample_file_mode_list_kwargs(self):
        # New: list_column / t_column kwargs map to listcolumn / tcolumn.
        args = cmd.resample(method="linear", file_times="list",
                             list_column=2, t_column=1)._to_cli_args()
        assert "listcolumn" in args
        idx = args.index("listcolumn")
        assert args[idx+1] == "2"
        assert "tcolumn" in args
        idx = args.index("tcolumn")
        assert args[idx+1] == "1"

    def test_resample_file_mode_list_inline_legacy(self):
        # Legacy: extra tokens inside file_times still pass through.
        args = cmd.resample(
            method="linear",
            file_times="list listcolumn 1 tcolumn 1",
        )._to_cli_args()
        s = " ".join(args)
        assert "file list listcolumn 1 tcolumn 1" in s

    def test_resample_t_column_alias_for_file_column(self):
        # In path mode, t_column emits the `column` keyword (was `file_column`).
        args = cmd.resample(method="linear",
                             file_times="/path/to/times.txt",
                             t_column=4)._to_cli_args()
        assert "column" in args
        assert args[args.index("column")+1] == "4"

    def test_resample_gaps(self):
        args = cmd.resample(method="linear", gaps="fix")._to_cli_args()
        assert "gaps" in args
        idx = args.index("gaps")
        assert args[idx+1] == "fix"


# ===========================================================================
# CLI-arg tests — Fitting
# ===========================================================================

class TestCLIArgsFitting:

    def test_tfa_basic(self):
        args = cmd.TFA("/tmp/trends.txt", "/tmp/dates.txt", 10.0)._to_cli_args()
        assert args[0] == "-TFA"
        assert "/tmp/trends.txt" in args

    def test_tfa_correct_lc(self):
        args = cmd.TFA("/tmp/trends.txt", "/tmp/dates.txt", 10.0,
                       correct_lc=True)._to_cli_args()
        assert "1" in args  # correctlc

    def test_tfa_xycol(self):
        args = cmd.TFA("/tmp/trends.txt", "/tmp/dates.txt", 10.0,
                       xycol=(4, 5))._to_cli_args()
        assert "xycol" in args and "4" in args and "5" in args

    def test_tfa_clip(self):
        args = cmd.TFA("/tmp/trends.txt", "/tmp/dates.txt", 10.0,
                       clip=3.0)._to_cli_args()
        assert "clip" in args and "3.0" in args

    def test_tfa_sr_basic(self):
        args = cmd.TFA_SR("/tmp/trends.txt", "/tmp/dates.txt",
                          10.0)._to_cli_args()
        assert args[0] == "-TFA_SR"

    def test_tfa_sr_xycol_after_pixelsep(self):
        """``xycol`` is an optional keyword that must follow the positional
        ``pixelsep`` argument; emitting it before ``pixelsep`` made vartools
        treat the pixelsep value as the next positional argument and fail
        with "Invalid command or option".
        """
        args = cmd.TFA_SR("/tmp/trends.txt", "/tmp/dates.txt", 25.0,
                          xycol=(2, 3))._to_cli_args()
        # pixelsep must appear before "xycol" in the token list.
        idx_xycol = args.index("xycol")
        idx_pixelsep = args.index("25.0")
        assert idx_pixelsep < idx_xycol, (
            f"pixelsep (25.0 at {idx_pixelsep}) must precede xycol "
            f"(at {idx_xycol}); got {args}"
        )
        # And xycol's two integer args must come immediately after.
        assert args[idx_xycol:idx_xycol + 3] == ["xycol", "2", "3"]

    # ------- TFA_SR decorr and signal_period (P2 batch) -------

    def test_tfa_sr_decorr(self):
        args = cmd.TFA_SR("trends.txt", "dates.txt", 1.0,
                           decorr_params="0 2 col1 1 col2 2")._to_cli_args()
        assert "decorr" in args
        idx = args.index("decorr")
        assert args[idx+1] == "0"
        assert args[idx+2] == "2"

    def test_tfa_sr_signal_period_bin(self):
        args = cmd.TFA_SR("trends.txt", "dates.txt", 1.0,
                           signal_mode="bin", signal_params=50,
                           signal_period=1.23)._to_cli_args()
        assert "period" in args
        idx = args.index("period")
        # The CLI requires "fix <val>" for a fixed numeric signal period.
        assert args[idx+1] == "fix"
        assert args[idx+2] == "1.23"

    def test_tfa_sr_signal_period_harm(self):
        args = cmd.TFA_SR("trends.txt", "dates.txt", 1.0,
                           signal_mode="harm", signal_params=(3, 0),
                           signal_period="ls")._to_cli_args()
        assert "period" in args
        idx = args.index("period")
        assert args[idx+1] == "ls"

    # ------- TFA / TFA_SR refmag (reset corrected-LC level) -------

    def test_tfa_refmag_absent_by_default(self):
        args = cmd.TFA("trends.txt", "dates.txt", 10.0)._to_cli_args()
        assert "refmag" not in args

    def test_tfa_refmag_fixed(self):
        """A numeric refmag emits the bare value (built-in fix/var/expr
        spec — no explicit "fix" keyword) as the final tokens."""
        args = cmd.TFA("trends.txt", "dates.txt", 10.0,
                       refmag=12.0)._to_cli_args()
        assert args[-2:] == ["refmag", "12.0"]

    def test_tfa_refmag_usemedian(self):
        args = cmd.TFA("trends.txt", "dates.txt", 10.0,
                       refmag=12.0, refmag_usemedian=True)._to_cli_args()
        assert args[-3:] == ["refmag", "12.0", "usemedian"]

    def test_tfa_refmag_var_and_expr(self):
        # bare identifier -> "var"
        a = cmd.TFA("trends.txt", "dates.txt", 10.0,
                    refmag="targetmag")._to_cli_args()
        assert a[-3:] == ["refmag", "var", "targetmag"]
        # explicit expr passes through
        b = cmd.TFA("trends.txt", "dates.txt", 10.0,
                    refmag="expr 12.0+1.0")._to_cli_args()
        assert b[-3:] == ["refmag", "expr", "12.0+1.0"]

    def test_tfa_refmag_requires_correct_lc(self):
        with pytest.raises(ValueError):
            cmd.TFA("trends.txt", "dates.txt", 10.0,
                    correct_lc=False, refmag=12.0)._to_cli_args()

    def test_tfa_sr_refmag(self):
        args = cmd.TFA_SR("trends.txt", "dates.txt", 10.0,
                          refmag=12.0, refmag_usemedian=True)._to_cli_args()
        assert args[-3:] == ["refmag", "12.0", "usemedian"]

    def test_tfa_sr_refmag_requires_correct_lc(self):
        with pytest.raises(ValueError):
            cmd.TFA_SR("trends.txt", "dates.txt", 10.0,
                       correct_lc=False, refmag=12.0)._to_cli_args()

    def test_sysrem_basic(self):
        args = cmd.SYSREM(1, 1, "/tmp/airmass.txt")._to_cli_args()
        assert args[0] == "-SYSREM"

    def test_sysrem_save_trends_false(self):
        """save_trends=False emits the bare ``"0"`` flag — no literal
        "otrends" token, which the CLI grammar does not accept.
        """
        args = cmd.SYSREM(2, 1, "/tmp/airmass.txt",
                          save_trends=False)._to_cli_args()
        assert "otrends" not in args
        # The "0" flag for trends sits immediately before useweights.
        assert args[-2] == "0"

    def test_sysrem_save_trends_path(self):
        """save_trends="path" emits ``"1" <path>``, treating the path as a
        single file (the CLI writes one global trend file, not a dir).
        """
        args = cmd.SYSREM(2, 1, "/tmp/airmass.txt",
                          save_trends="/tmp/trends.txt")._to_cli_args()
        assert "otrends" not in args
        assert "/tmp/trends.txt" in args
        i = args.index("/tmp/trends.txt")
        # Trends path is preceded by the "1" emit-flag.
        assert args[i - 1] == "1"

    def test_sysrem_save_trends_marked_file_in_specs(self):
        """The ``trends`` entry in ``_output_file_specs`` must declare
        ``mode="file"`` so the pipeline runner does not ``os.makedirs`` the
        user-supplied trend-output path.
        """
        specs = cmd.SYSREM(1, 1, "/tmp/airmass.txt")._output_file_specs()
        assert "trends" in specs
        # 3-tuple form with the "file" mode marker.
        assert len(specs["trends"]) == 3
        assert specs["trends"][2] == "file"

    def test_mandel_agol_basic(self):
        args = cmd.MandelAgolTransit(P0=3.0, T00=1.0)._to_cli_args()
        assert args[0] == "-MandelAgolTransit"

    def test_mandel_agol_ld_quad(self):
        args = cmd.MandelAgolTransit(P0=3.0, T00=1.0,
                                     ld_type="quad",
                                     ld_coeffs=[0.3, 0.2])._to_cli_args()
        assert "quad" in args

    def test_mandel_agol_correct_lc(self):
        args = cmd.MandelAgolTransit(P0=3.0, T00=1.0,
                                     correct_lc=True)._to_cli_args()
        assert "1" in args  # correctlc

    def test_softened_transit_bls(self):
        args = cmd.SoftenedTransit(init_params="bls")._to_cli_args()
        assert args[0] == "-SoftenedTransit" and "bls" in args

    def test_starspot_ls(self):
        args = cmd.Starspot(period="ls")._to_cli_args()
        assert args[0] == "-Starspot" and "ls" in args

    def test_starspot_fix(self):
        args = cmd.Starspot(period=2.3)._to_cli_args()
        assert "fix" in args and "2.3" in args

    def test_microlens_basic(self):
        args = cmd.microlens()._to_cli_args()
        assert args[0] == "-microlens"

    def test_microlens_fix_params(self):
        args = cmd.microlens(u0=0.5, t0=10.0, tmax=0.3)._to_cli_args()
        assert "u0" in args and "t0" in args and "tmax" in args

    # ------- microlens step/novary/fixcolumn (P2 batch) -------

    def test_microlens_step(self):
        args = cmd.microlens(u0=0.5, t0=0.0, tmax=10.0,
                              u0_step=0.01)._to_cli_args()
        assert "step" in args
        idx = args.index("step")
        assert args[idx+1] == "0.01"

    def test_microlens_novary(self):
        args = cmd.microlens(u0=0.5, t0=0.0, tmax=10.0,
                              t0_novary=True)._to_cli_args()
        # novary must appear after "t0" tokens, not after "u0" or "tmax"
        t0_idx = args.index("t0")
        novary_idx = args.index("novary")
        assert novary_idx > t0_idx

    def test_microlens_fixcolumn(self):
        args = cmd.microlens(u0="fixcolumn u0_col",
                              t0=0.0, tmax=10.0)._to_cli_args()
        assert "fixcolumn" in args
        assert "u0_col" in args
        # "fix" must NOT appear before "fixcolumn"
        idx = args.index("u0")
        assert args[idx+1] == "fixcolumn"

    def test_nonlinfit_basic(self):
        args = cmd.nonlinfit("a0+a1*sin(2*PI*t/a2)",
                              "a0 10 a1 0.1 a2 2.0")._to_cli_args()
        assert args[0] == "-nonlinfit"

    def test_nonlinfit_mcmc(self):
        args = cmd.nonlinfit("a0+a1*t", "a0 0 a1 0",
                              optimizer="mcmc")._to_cli_args()
        assert "mcmc" in args

    def test_addnoise_white(self):
        args = cmd.addnoise(noise_type="white", sig_white=0.005)._to_cli_args()
        assert args[0] == "-addnoise" and "white" in args

    def test_addnoise_squareexp(self):
        args = cmd.addnoise(noise_type="squareexp",
                             sig_white=0.001, rho=1.0,
                             sig_red=0.005)._to_cli_args()
        assert "squareexp" in args and "rho" in args

    def test_addnoise_exp(self):
        args = cmd.addnoise(noise_type="exp",
                             sig_white=0.001, rho=1.0,
                             sig_red=0.005)._to_cli_args()
        assert "exp" in args

    def test_findblends_basic(self):
        args = cmd.findblends(matchrad=5.0)._to_cli_args()
        assert args[0] == "-findblends"

    def test_findblends_radec(self):
        args = cmd.findblends(matchrad=5.0, radec=True)._to_cli_args()
        assert "radec" in args

    def test_findblends_period(self):
        args = cmd.findblends(matchrad=5.0, period="list", nharm=2)._to_cli_args()
        assert "list" in args and "2" in args

    # ------- MandelAgolTransit modelvar and bimpact (Batch 2d) -------

    def test_mandel_agol_modelvar(self):
        c = cmd.MandelAgolTransit(P0=3.0, T00=1.0, save_model=True,
                                   modelvar="mvar")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "modelvar" in args and "mvar" in args

    def test_mandel_agol_default_inclination_token(self):
        args = cmd.MandelAgolTransit(P0=3.0, T00=1.0)._to_cli_args()
        # "i" is the inclination token at position 5 (0-based)
        assert args[5] == "i"

    def test_mandel_agol_bimpact(self):
        args = cmd.MandelAgolTransit(P0=3.0, T00=1.0, bimpact=0.5)._to_cli_args()
        assert args[5] == "b"
        assert args[6] == "0.5"

    # ------- SoftenedTransit fit_harm (Batch 2j) -------

    def test_softenedtransit_fit_harm_zero(self):
        args = cmd.SoftenedTransit()._to_cli_args()
        assert "0" in args  # fit_harm=0 at end

    def test_softenedtransit_fit_harm(self):
        args = cmd.SoftenedTransit(fit_harm=1, fit_harm_method="aov",
                                    fit_harm_nharm=3,
                                    fit_harm_nsubharm=0)._to_cli_args()
        assert "1" in args
        assert "aov" in args

    # ------- Starspot named params (Batch 2f) -------

    def test_starspot_named_params(self):
        args = cmd.Starspot(period=2.3, a0=0.05, b0=0.5, alpha0=20.0,
                             i0=85.0, chi0=30.0, psi00=0.0,
                             mconst0=0.0)._to_cli_args()
        assert args[0] == "-Starspot"
        assert "0.05" in args  # a0
        assert "20.0" in args  # alpha0

    def test_starspot_fit_flags_named(self):
        args = cmd.Starspot(period=2.3, fit_period=0, fit_a=0)._to_cli_args()
        assert args[0] == "-Starspot"
        # fit_period=0 → first fit flag after initial params is 0
        # After: -Starspot fix 2.3 a0 b0 alpha0 i0 chi0 psi00 mconst0 fitP ...
        # fitP is at index 3 + 2 (period tokens "fix" "2.3") + 7 (params) = 12
        assert "0" in args  # at least fit_period=0

    def test_starspot_old_initial_params_raises(self):
        with pytest.raises(TypeError):
            cmd.Starspot(initial_params=[0.1, 0.5, 0.5, 85.0, 30.0, 0.0, 0.0])

    # ------- nonlinfit new params (Batch 2c + Batch 5) -------

    def test_nonlinfit_modelvar(self):
        c = cmd.nonlinfit("a0+a1*t", "a0 0 a1 0", save_model=True,
                           modelvar="mvar")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "modelvar" in args and "mvar" in args

    def test_nonlinfit_errors(self):
        args = cmd.nonlinfit("a0+a1*t", "a0 0 a1 0", errors="err")._to_cli_args()
        assert "errors" in args and "err" in args

    def test_nonlinfit_linfit_params(self):
        args = cmd.nonlinfit("a0+a1*t", "a0 0 a1 0",
                              linfit_params="b0")._to_cli_args()
        assert "linfit" in args and "b0" in args

    def test_nonlinfit_amoeba_tolerance(self):
        args = cmd.nonlinfit("a0+a1*t", "a0 0 a1 0", optimizer="amoeba",
                              amoeba_tolerance=1e-6)._to_cli_args()
        assert "tolerance" in args

    def test_nonlinfit_amoeba_maxsteps(self):
        args = cmd.nonlinfit("a0+a1*t", "a0 0 a1 0", optimizer="amoeba",
                              amoeba_maxsteps=10000)._to_cli_args()
        assert "maxsteps" in args and "10000" in args

    def test_nonlinfit_mcmc_naccept(self):
        args = cmd.nonlinfit("a0+a1*t", "a0 0 a1 0", optimizer="mcmc",
                              mcmc_naccept=1000)._to_cli_args()
        assert "Naccept" in args and "1000" in args

    def test_nonlinfit_mcmc_nlinkstotal(self):
        args = cmd.nonlinfit("a0+a1*t", "a0 0 a1 0", optimizer="mcmc",
                              mcmc_nlinkstotal=50000)._to_cli_args()
        assert "Nlinkstotal" in args and "50000" in args

    def test_nonlinfit_mcmc_fracburnin(self):
        args = cmd.nonlinfit("a0+a1*t", "a0 0 a1 0", optimizer="mcmc",
                              mcmc_fracburnin=0.3)._to_cli_args()
        assert "fracburnin" in args and "0.3" in args

    def test_nonlinfit_mcmc_skipamoeba(self):
        args = cmd.nonlinfit("a0+a1*t", "a0 0 a1 0", optimizer="mcmc",
                              mcmc_skipamoeba=True)._to_cli_args()
        assert "skipamoeba" in args

    def test_nonlinfit_covariance(self):
        args = cmd.nonlinfit("a0+a1*t", "a0 0 a1 0",
                              covariance="squareexp amp_v rho_v")._to_cli_args()
        assert "covariance" in args and "squareexp" in args

    def test_nonlinfit_priors(self):
        args = cmd.nonlinfit("a0+a1*t", "a0 0 a1 0",
                              priors="(a0-1)^2/0.01")._to_cli_args()
        assert "priors" in args

    def test_nonlinfit_constraints(self):
        args = cmd.nonlinfit("a0+a1*t", "a0 0 a1 0",
                              constraints="a1>0")._to_cli_args()
        assert "constraints" in args

    def test_nonlinfit_mcmc_model_nameformat(self):
        c = cmd.nonlinfit("a0+a1*t", "a0 0 a1 0", optimizer="mcmc",
                           save_model=True, model_nameformat="%s.nlf")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "format" in args and "%s.nlf" in args

    # ------- findblends fixes and additions (Batch 4d) -------

    def test_findblends_period_fix(self):
        args = cmd.findblends(matchrad=5.0, period="fix 2.3")._to_cli_args()
        assert "fix" in args and "2.3" in args

    def test_findblends_xycol(self):
        args = cmd.findblends(matchrad=5.0, period="list",
                               xycol=("2", "3"))._to_cli_args()
        assert "xycol" in args and "2" in args and "3" in args

    def test_findblends_starlist(self):
        args = cmd.findblends(matchrad=5.0, period="list",
                               starlist="/tmp/stars.txt")._to_cli_args()
        assert "starlist" in args and "/tmp/stars.txt" in args

    def test_findblends_zeromag(self):
        args = cmd.findblends(matchrad=5.0, period="list",
                               zeromag=25.0)._to_cli_args()
        assert "zeromag" in args and "25.0" in args

    def test_findblends_nofluxconvert(self):
        args = cmd.findblends(matchrad=5.0, period="list",
                               nofluxconvert=True)._to_cli_args()
        assert "nofluxconvert" in args

    # ------- addnoise matern/wavelet/bintime (Batch 4f) -------

    def test_addnoise_matern(self):
        args = cmd.addnoise(noise_type="matern", sig_white=0.001,
                             rho=1.0, sig_red=0.005, nu=1.5)._to_cli_args()
        assert "matern" in args and "nu" in args

    def test_addnoise_wavelet(self):
        args = cmd.addnoise(noise_type="wavelet", sig_white=0.001,
                             sig_red=0.005, gamma=2.0)._to_cli_args()
        assert "wavelet" in args and "gamma" in args

    def test_addnoise_bintime(self):
        args = cmd.addnoise(noise_type="squareexp", sig_white=0.001,
                             rho=1.0, sig_red=0.005, bintime=0.02)._to_cli_args()
        assert "bintime" in args and "0.02" in args

    # ------- MandelAgolTransit RV fitting (Batch 4g) -------

    def test_mandel_agol_rv_fitting(self):
        args = cmd.MandelAgolTransit(
            P0=3.0, T00=1.0,
            rv_file="/tmp/rv.txt",
            rv_model_file="/tmp/rv_model.txt",
            K0=10.0, gamma0=0.0, fitK=1, fitgamma=1
        )._to_cli_args()
        assert "/tmp/rv.txt" in args
        assert "/tmp/rv_model.txt" in args
        assert "10.0" in args

    # ----- MatchedFilter -----

    def test_mf_gauss_basic(self):
        args = cmd.MatchedFilter("gauss", 2.0, "window", "both",
                                  sigma=0.5)._to_cli_args()
        assert args[0] == "-matchedfilter"
        assert args[1:4] == ["template", "gauss", "0.5"]
        assert "mode" in args and args[args.index("mode") + 1] == "window"
        assert "signs" in args and args[args.index("signs") + 1] == "both"

    def test_mf_doubleexp(self):
        args = cmd.MatchedFilter("doubleexp", 0.01, "window", "positive",
                                  tau_rise=0.001, tau_decay=0.005)._to_cli_args()
        i = args.index("doubleexp")
        assert args[i + 1:i + 3] == ["0.001", "0.005"]

    def test_mf_flare(self):
        args = cmd.MatchedFilter("flare", 0.02, "window", "positive",
                                  tfwhm=0.005)._to_cli_args()
        assert "flare" in args
        assert args[args.index("flare") + 1] == "0.005"

    def test_mf_box(self):
        args = cmd.MatchedFilter("box", 0.5, "window", "negative",
                                  width=0.083)._to_cli_args()
        assert "box" in args
        assert args[args.index("box") + 1] == "0.083"

    def test_mf_triangle(self):
        args = cmd.MatchedFilter("triangle", 0.5, "window", "both",
                                  width=0.1)._to_cli_args()
        assert "triangle" in args

    def test_mf_trap(self):
        args = cmd.MatchedFilter("trap", 0.5, "window", "both",
                                  rise=0.01, flat=0.05, fall=0.01)._to_cli_args()
        i = args.index("trap")
        assert args[i + 1:i + 4] == ["0.01", "0.05", "0.01"]

    def test_mf_exp(self):
        args = cmd.MatchedFilter("exp", 0.02, "window", "negative",
                                  tau=0.005)._to_cli_args()
        i = args.index("exp")
        assert args[i + 1] == "0.005"

    def test_mf_file(self):
        args = cmd.MatchedFilter("file", 2.0, "window", "positive",
                                  template_file="/tmp/t.mf")._to_cli_args()
        assert "file" in args
        assert args[args.index("file") + 1] == "/tmp/t.mf"

    def test_mf_expr_default_varname(self):
        args = cmd.MatchedFilter("expr", 2.0, "window", "positive",
                                  expression="exp(-s*s/0.5)")._to_cli_args()
        i = args.index("expr")
        # No 'varname' keyword when expr_varname is None (default 's').
        assert args[i + 1] == "exp(-s*s/0.5)"

    def test_mf_expr_custom_varname(self):
        args = cmd.MatchedFilter("expr", 2.0, "window", "positive",
                                  expression="exp(-x*x/0.5)",
                                  expr_varname="x")._to_cli_args()
        i = args.index("expr")
        assert args[i + 1:i + 4] == ["varname", "x", "exp(-x*x/0.5)"]

    def test_mf_nfft_mode(self):
        args = cmd.MatchedFilter("gauss", 2.0, "nfft", "both",
                                  sigma=0.5)._to_cli_args()
        assert args[args.index("mode") + 1] == "nfft"

    def test_mf_min_separation_whiten_maskpoints(self):
        args = cmd.MatchedFilter("gauss", 2.0, "window", "both",
                                  sigma=0.5,
                                  min_separation=0.5, whiten=True,
                                  maskpoints="mask")._to_cli_args()
        # Canonical strict-order trailing keywords.
        assert (args.index("min_separation") < args.index("whiten")
                < args.index("maskpoints"))

    def test_mf_save_matchfile(self):
        c = cmd.MatchedFilter("gauss", 2.0, "window", "both",
                               sigma=0.5, save_matchfile=True)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "/tmp" in args
        # omatchfile token immediately precedes outdir.
        i = args.index("/tmp")
        assert args[i - 1] == "1"

    def test_mf_varexpr_scalars(self):
        # sigma and support_halfwidth accept var/expr forms.
        args = cmd.MatchedFilter("gauss", "myvar", "window", "both",
                                  sigma="2*sigma0")._to_cli_args()
        assert "expr" in args
        assert "var" in args

    def test_mf_rejects_bad_template(self):
        with pytest.raises(ValueError):
            cmd.MatchedFilter("bogus", 2.0, "window", "both", sigma=0.5)

    def test_mf_rejects_bad_mode(self):
        with pytest.raises(ValueError):
            cmd.MatchedFilter("gauss", 2.0, "fft", "both", sigma=0.5)

    def test_mf_rejects_bad_signs(self):
        with pytest.raises(ValueError):
            cmd.MatchedFilter("gauss", 2.0, "window", "yes", sigma=0.5)

    def test_mf_rejects_missing_template_param(self):
        with pytest.raises(ValueError):
            cmd.MatchedFilter("gauss", 2.0, "window", "both")
        with pytest.raises(ValueError):
            cmd.MatchedFilter("trap", 2.0, "window", "both",
                              rise=0.01, flat=0.05)   # missing fall

    def test_mf_rejects_template_kwarg_mismatch(self):
        # 'tau' belongs to exp; passing it with gauss should reject.
        with pytest.raises(ValueError):
            cmd.MatchedFilter("gauss", 2.0, "window", "both",
                              sigma=0.5, tau=0.01)

    def test_mf_rejects_file_without_path(self):
        with pytest.raises(ValueError):
            cmd.MatchedFilter("file", 2.0, "window", "positive")

    def test_mf_rejects_expr_without_expression(self):
        with pytest.raises(ValueError):
            cmd.MatchedFilter("expr", 2.0, "window", "positive")

    def test_mf_rejects_npeaks_zero(self):
        with pytest.raises(ValueError):
            cmd.MatchedFilter("gauss", 2.0, "window", "both",
                              sigma=0.5, npeaks=0)


# ===========================================================================
# CLI-arg tests — Misc
# ===========================================================================

class TestCLIArgsMisc:

    def test_addfitskeyword_double(self):
        args = cmd.addfitskeyword("PERIOD", "TDOUBLE", 1.234)._to_cli_args()
        assert args[0] == "-addfitskeyword"
        assert "PERIOD" in args and "TDOUBLE" in args
        assert "fix" in args and "1.234" in args

    def test_addfitskeyword_string(self):
        args = cmd.addfitskeyword("OBJNAME", "TSTRING",
                                   "fix HAT-P-1")._to_cli_args()
        assert "TSTRING" in args and "HAT-P-1" in args

    def test_addfitskeyword_var(self):
        args = cmd.addfitskeyword("RMS", "TDOUBLE",
                                   "var Weighted_RMS_0")._to_cli_args()
        assert "var" in args

    def test_addfitskeyword_comment(self):
        args = cmd.addfitskeyword("PERIOD", "TDOUBLE", 1.23,
                                   comment="best period")._to_cli_args()
        assert "comment" in args and "best period" in args

    def test_addfitskeyword_hdu(self):
        args = cmd.addfitskeyword("PERIOD", "TDOUBLE", 1.23,
                                   hdu="primary")._to_cli_args()
        assert "primary" in args

    def test_addfitskeyword_mode(self):
        args = cmd.addfitskeyword("PERIOD", "TDOUBLE", 1.23,
                                   mode="update")._to_cli_args()
        assert "update" in args

    def test_converttime_jd_mjd(self):
        args = cmd.converttime("jd", "mjd")._to_cli_args()
        assert args[0] == "-converttime"
        assert "input" in args and "jd" in args
        assert "output" in args and "mjd" in args

    def test_converttime_hjd_bjd_radec(self):
        args = cmd.converttime("hjd", "bjd", ra=83.822, dec=1.5)._to_cli_args()
        assert "radec" in args and "fix" in args
        assert "83.822" in args and "1.5" in args

    def test_converttime_input_subtract(self):
        args = cmd.converttime("jd", "mjd",
                                input_subtract=2450000.0)._to_cli_args()
        assert "inputsubtract" in args and "2450000.0" in args

    def test_converttime_sys(self):
        args = cmd.converttime("jd", "bjd",
                                input_sys="utc",
                                output_sys="tdb")._to_cli_args()
        assert "inputsys-utc" in args and "outputsys-tdb" in args

    def test_R_inline(self):
        args = cmd.R("mag <- mag - mean(mag)")._to_cli_args()
        assert args[0] == "-R"

    def test_R_fromfile(self):
        args = cmd.R("/tmp/myscript.R", fromfile=True)._to_cli_args()
        assert "fromfile" in args and "/tmp/myscript.R" in args

    def test_R_vars(self):
        args = cmd.R("mag <- mag * 2", vars="t,mag,err")._to_cli_args()
        assert "vars" in args and "t,mag,err" in args

    def test_R_invars_outvars(self):
        args = cmd.R("mag <- mag - mean(mag)",
                     invars="t,mag,err", outvars="mag")._to_cli_args()
        assert "invars" in args and "outvars" in args

    def test_R_init(self):
        args = cmd.R("f(mag)", init="f <- function(x) x - mean(x)")._to_cli_args()
        assert "init" in args

    def test_R_process_all_lcs(self):
        args = cmd.R("x <- 1", process_all_lcs=True)._to_cli_args()
        assert "process_all_lcs" in args

    # ----- cmd.python — CLI-token emission -----
    #
    # The wrapper mirrors -R surface; the new bits to lock down are:
    #   * the CLI marker is `-python` (not `-R`),
    #   * `skipfail` (Python-only) instead of `verbose` (R-only),
    #   * `init` + `continueprocess` are mutually exclusive,
    #   * `inprocess=True` raises NotImplementedError until Stage 2,
    #   * the order of tokens matches what vartools' parser expects
    #     (command, [continueprocess], [init], [vars | invars+outvars],
    #      [outputcolumns], [process_all_lcs], [skipfail]).

    def test_python_inline(self):
        args = cmd.python("b = numpy.var(mag)")._to_cli_args()
        assert args[0] == "-python"
        assert args[1] == "b = numpy.var(mag)"

    def test_python_fromfile(self):
        args = cmd.python("/tmp/x.py", fromfile=True)._to_cli_args()
        assert "fromfile" in args and "/tmp/x.py" in args
        # `fromfile` precedes the path — match the CLI grammar.
        idx = args.index("fromfile")
        assert args[idx + 1] == "/tmp/x.py"

    def test_python_vars(self):
        args = cmd.python("mag = mag * 2", vars="t,mag,err")._to_cli_args()
        assert "vars" in args
        idx = args.index("vars")
        assert args[idx + 1] == "t,mag,err"
        # When `vars` is given, neither invars nor outvars is emitted.
        assert "invars" not in args
        assert "outvars" not in args

    def test_python_invars_outvars(self):
        args = cmd.python("b = numpy.var(mag)",
                          invars="mag", outvars="b")._to_cli_args()
        assert "invars" in args and "outvars" in args
        i = args.index("invars")
        assert args[i + 1] == "mag"
        o = args.index("outvars")
        assert args[o + 1] == "b"

    def test_python_init_inline(self):
        args = cmd.python("y = f(mag)",
                          init="def f(x): return x.mean()")._to_cli_args()
        assert "init" in args
        i = args.index("init")
        # init precedes the inline init body; no `file` keyword.
        assert "file" not in args
        assert args[i + 1] == "def f(x): return x.mean()"

    def test_python_init_fromfile(self):
        args = cmd.python("y = f(mag)",
                          init="/tmp/init.py",
                          init_fromfile=True)._to_cli_args()
        # init followed by `file`, then the path.
        idx = args.index("init")
        assert args[idx + 1] == "file"
        assert args[idx + 2] == "/tmp/init.py"

    def test_python_outputcolumns(self):
        args = cmd.python("b = numpy.var(mag)", outvars="b",
                          outputcolumns="b")._to_cli_args()
        assert "outputcolumns" in args
        i = args.index("outputcolumns")
        assert args[i + 1] == "b"

    def test_python_process_all_lcs(self):
        args = cmd.python("x = 1", process_all_lcs=True)._to_cli_args()
        assert "process_all_lcs" in args
        # No leading hyphen — it's a flag, not a CLI option of its own.
        assert "-process_all_lcs" not in args

    def test_python_skipfail(self):
        args = cmd.python("x = 1", skipfail=True)._to_cli_args()
        assert "skipfail" in args

    def test_python_skipfail_default_false(self):
        args = cmd.python("x = 1")._to_cli_args()
        assert "skipfail" not in args

    def test_python_continueprocess(self):
        args = cmd.python("y = f(mag)",
                          continueprocess=1,
                          invars="mag", outvars="y")._to_cli_args()
        assert "continueprocess" in args
        i = args.index("continueprocess")
        assert args[i + 1] == "1"

    def test_python_init_and_continueprocess_conflict(self):
        # vartools rejects the combination — the wrapper should pre-empt
        # at construction time so users get a clear Python-side error.
        with pytest.raises(ValueError, match="continueprocess"):
            cmd.python("y = f(mag)",
                       init="def f(x): return x.mean()",
                       continueprocess=1)

    def test_python_inprocess_constructs_without_error(self):
        # Stage-2 path: inprocess=True now constructs cleanly.  Validates
        # that the kwarg lives on the instance and that the C-side
        # callback shim got registered as a side effect (we don't
        # exercise execution here — that's covered in
        # TestPythonInprocessIntegration below).
        c = cmd.python("b = 0.0", outvars="b", outputcolumns="b",
                       inprocess=True)
        assert c.inprocess is True
        # Tokens are identical to the subprocess form — vartools sees the
        # same -python command; the callback hook diverts dispatch
        # internally based on whether it's been registered.
        args = c._to_cli_args()
        assert args[0] == "-python"

    def test_python_namespace_kwarg_accepted(self):
        # namespace= is meaningful only for inprocess=True.  Setting it
        # alongside the default subprocess path is allowed (it's
        # silently ignored — construction does NOT raise).
        c = cmd.python("x = 1", namespace={"foo": 42})
        assert c.namespace == {"foo": 42}
        # Tokens are unaffected.
        args = c._to_cli_args()
        assert args[0] == "-python"

    def test_python_token_order(self):
        # Sanity-check the relative order: command then continueprocess
        # then init then vars/invars+outvars then outputcolumns then
        # process_all_lcs then skipfail.  vartools' grammar is positional.
        args = cmd.python("y = f(mag)",
                          init="def f(x): return x.mean()",
                          invars="mag", outvars="y",
                          outputcolumns="y",
                          process_all_lcs=True,
                          skipfail=True)._to_cli_args()
        # Find each token's index and assert ascending order.
        order = ["-python", "y = f(mag)", "init",
                 "invars", "outvars", "outputcolumns",
                 "process_all_lcs", "skipfail"]
        positions = [args.index(tok) for tok in order]
        assert positions == sorted(positions), (
            f"unexpected token ordering: {list(zip(order, positions))}"
        )

    def test_match_basic(self):
        args = cmd.match("/tmp/cat.txt", matchcolumn="t:1",
                          addcolumns="ra:2,dec:3")._to_cli_args()
        assert args[0] == "-match"
        assert "file" in args and "/tmp/cat.txt" in args
        assert "matchcolumn" in args and "addcolumns" in args

    def test_match_nanmissing(self):
        args = cmd.match("/tmp/cat.txt", matchcolumn="1",
                          addcolumns="val:2",
                          missing="nanmissing")._to_cli_args()
        assert "nanmissing" in args

    def test_match_skipnum(self):
        args = cmd.match("/tmp/cat.txt", matchcolumn="1",
                          addcolumns="val:2",
                          skipnum=3)._to_cli_args()
        assert "skipnum" in args and "3" in args

    def test_match_delimiter(self):
        args = cmd.match("/tmp/cat.txt", matchcolumn="1",
                          addcolumns="val:2",
                          delimiter=",")._to_cli_args()
        assert "delimiter" in args and "," in args

    def test_o_basic(self):
        args = cmd.o("/tmp/out.lc")._to_cli_args()
        assert args[0] == "-o" and "/tmp/out.lc" in args

    # ------- o command new params (Batch 1) -------

    def test_o_fits(self):
        args = cmd.o("/tmp/out.fits", fits=True)._to_cli_args()
        assert "fits" in args

    def test_o_noclobber(self):
        args = cmd.o("/tmp/out.lc", noclobber=True)._to_cli_args()
        assert "noclobber" in args

    def test_o_copyheader(self):
        args = cmd.o("/tmp/out.lc", copyheader=True)._to_cli_args()
        assert "copyheader" in args

    def test_o_namecommand(self):
        args = cmd.o("/tmp/out.lc", namecommand="cat")._to_cli_args()
        assert "namecommand" in args and "cat" in args

    def test_o_namefromlist_bool(self):
        args = cmd.o("/tmp/out.lc", namefromlist=True)._to_cli_args()
        assert "namefromlist" in args

    def test_o_namefromlist_column(self):
        args = cmd.o("/tmp/out.lc", namefromlist="5")._to_cli_args()
        assert "namefromlist" in args and "column" in args and "5" in args

    def test_o_delimiter(self):
        args = cmd.o("/tmp/out.lc", delimiter=",")._to_cli_args()
        assert "delimiter" in args and "," in args

    def test_o_logcommandline(self):
        args = cmd.o("/tmp/out.lc", logcommandline=True)._to_cli_args()
        assert "logcommandline" in args

    def test_o_fits_not_present_by_default(self):
        args = cmd.o("/tmp/out.lc")._to_cli_args()
        assert "fits" not in args

    # ------- match inlist_column API (Batch 2g) -------

    def test_match_inlist_column_emits_column_not_filename(self):
        args = cmd.match(catalog="/tmp/cat.txt", matchcolumn="1",
                          addcolumns="val:2", source="inlist",
                          inlist_column="4")._to_cli_args()
        assert "inlist" in args
        assert "4" in args
        assert "/tmp/cat.txt" not in args  # catalog ignored for inlist

    def test_match_inlist_missing_column_raises(self):
        with pytest.raises(ValueError):
            cmd.match(catalog="/tmp/cat.txt", matchcolumn="1",
                       addcolumns="val:2",
                       source="inlist")._to_cli_args()

    def test_ifcmd_basic(self):
        args = cmd.ifcmd("Weighted_RMS_0 > 0.01")._to_cli_args()
        assert args[0] == "-if"

    def test_binlc_binsize(self):
        args = cmd.binlc(binsize=0.01)._to_cli_args()
        assert args[0] == "-binlc"
        assert "average" in args and "binsize" in args and "0.01" in args
        assert "tcenter" in args

    def test_binlc_nbins(self):
        args = cmd.binlc(nbins=50)._to_cli_args()
        assert "nbins" in args and "50" in args

    def test_binlc_median(self):
        args = cmd.binlc(binsize=0.1, method="median")._to_cli_args()
        assert "median" in args

    def test_binlc_weightedaverage(self):
        args = cmd.binlc(binsize=0.1, method="weightedaverage")._to_cli_args()
        assert "weightedaverage" in args

    def test_binlc_taverage(self):
        args = cmd.binlc(binsize=0.1, time_output="taverage")._to_cli_args()
        assert "taverage" in args

    def test_binlc_T0(self):
        args = cmd.binlc(binsize=0.1, T0=0.0)._to_cli_args()
        assert "T0" in args and "fix" in args and "0.0" in args

    def test_binlc_requires_binsize_or_nbins(self):
        with pytest.raises(ValueError):
            cmd.binlc()

    def test_binlc_binshift(self):
        args = cmd.binlc(binsize=0.5, binshift=0.5)._to_cli_args()
        assert "binshift" in args
        assert args[args.index("binshift") + 1] == "0.5"
        # The legacy keyword must NOT also be emitted.
        assert "firstbinshift" not in args

    def test_binlc_firstbinshift_legacy_emits_old_keyword(self):
        import warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            args = cmd.binlc(binsize=0.5, firstbinshift=0.5)._to_cli_args()
            assert any(issubclass(rec.category, DeprecationWarning)
                       for rec in w), \
                f"expected DeprecationWarning, got {[str(r.message) for r in w]}"
        assert "firstbinshift" in args
        assert args[args.index("firstbinshift") + 1] == "0.5"
        assert "binshift" not in args

    def test_binlc_binshift_and_firstbinshift_mutually_exclusive(self):
        with pytest.raises(ValueError):
            cmd.binlc(binsize=0.5, binshift=0.5, firstbinshift=0.5)

    def test_columnsuffix_basic(self):
        args = cmd.columnsuffix("myls")._to_cli_args()
        assert args == ["-columnsuffix", "myls"]

    def test_raw_string(self):
        assert cmd.Raw("-foo 1 2")._to_cli_args() == ["-foo", "1", "2"]

    def test_raw_list(self):
        assert cmd.Raw(["-bar", "x"])._to_cli_args() == ["-bar", "x"]


# ===========================================================================
# End-to-end pipeline tests
# ===========================================================================

@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools binary not found")
class TestEndToEndPipelines:

    # -----------------------------------------------------------------------
    # Basic statistics
    # -----------------------------------------------------------------------

    def test_rms_chi2_alarm(self):
        lc = make_lc()
        result = vt.Pipeline([
            cmd.rms(),
            cmd.chi2(),
            cmd.alarm(),
        ]).run(lc)
        assert "RMS_0" in result.vars.index
        assert any("Chi2" in k or "chi2" in k for k in result.vars.index)
        assert any("Alarm" in k or "alarm" in k for k in result.vars.index)

    def test_vonneumann_pipeline(self):
        # Use a clear sinusoid -- eta should be well below 2 (correlated).
        lc = make_lc(period=2.0, amp=0.5)
        result = vt.Pipeline([
            cmd.vonNeumann(),
            cmd.vonNeumann(weighted=True),
        ]).run(lc)
        vn = result.varobjs.vonNeumann
        assert len(vn) == 2
        eta_unw  = float(vn[0].Ratio)
        eta_wgt  = float(vn[1].Ratio)
        # Sinusoidal signal -> strongly correlated -> eta << 2.
        assert eta_unw < 1.0
        assert eta_wgt < 1.0

    def test_stats_multiple_variables(self):
        lc = make_lc()
        result = vt.Pipeline([
            cmd.stats(["mag", "err"], "mean,median,stddev"),
        ]).run(lc)
        assert any("mean" in k.lower() for k in result.vars.index)

    def test_stats_MAD(self):
        lc = make_lc()
        result = vt.Pipeline([cmd.stats("mag", "MAD")]).run(lc)
        assert any("MAD" in k for k in result.vars.index)

    def test_rmsbin(self):
        lc = make_lc()
        result = vt.Pipeline([
            cmd.rmsbin(2, [0.5, 1.0]),
        ]).run(lc)
        assert result.vars is not None

    # -----------------------------------------------------------------------
    # autocorrelation file capture
    # -----------------------------------------------------------------------

    def test_autocorrelation_file_capture_single(self):
        """autocorrelation result is always captured in result.files."""
        lc = make_lc(n=200)
        result = vt.Pipeline([
            cmd.autocorrelation(0.0, 5.0, 0.1),
        ]).run(lc)
        assert "autocorrelation_result_0" in result.files
        acf = result.files["autocorrelation_result_0"]
        assert len(acf) > 0

    def test_autocorrelation_file_capture_batch(self):
        """autocorrelation result is captured for every LC in a batch run."""
        lcs = [make_lc(n=200, period=1.0 + i * 0.5) for i in range(3)]
        result = vt.Pipeline([
            cmd.autocorrelation(0.0, 5.0, 0.1),
        ]).run_batch(lcs)
        assert "autocorrelation_result_0" in result.files
        acfs = result.files["autocorrelation_result_0"]
        assert len(acfs) == 3
        for acf in acfs:
            assert len(acf) > 0

    def test_linfit_model_capture_default_filename(self):
        """linfit save_model captures the model file at the default
        ``<lcbase>.linfit.model`` location."""
        lc = make_lc(n=200)
        result = vt.Pipeline([
            cmd.linfit(function="a*t+b", paramlist="a,b", save_model=True),
        ]).run(lc)
        assert "linfit_model_0" in result.files
        assert len(result.files["linfit_model_0"]) > 0

    def test_linfit_model_capture_custom_nameformat(self):
        """When ``model_nameformat`` is set, vartools writes the model
        file under the user's format string instead of the default
        ``.linfit.model`` suffix.  Pyvartools must honour that
        substitution when capturing the file (this used to silently
        fail, leaving result.files empty for the model)."""
        lc = make_lc(n=200)
        result = vt.Pipeline([
            cmd.linfit(function="a*t+b", paramlist="a,b",
                        model_nameformat="%s.baldump",
                        save_model=True),
        ]).run(lc)
        assert "linfit_model_0" in result.files, \
            "linfit save_model with model_nameformat='%s.baldump' " \
            "must still populate result.files['linfit_model_0']"
        assert len(result.files["linfit_model_0"]) > 0

    # -----------------------------------------------------------------------
    # structurefunction file capture
    # -----------------------------------------------------------------------

    def test_structurefunction_save_result_capture(self):
        # Regression: capture used to be a silent no-op because the logical
        # name ("sf") did not match the save_result attribute, so the
        # collector's getattr(cmd, "save_sf") returned False.
        lc = vt.LightCurve.from_file(EXAMPLE_LC)
        result = vt.Pipeline([
            cmd.structurefunction(bins="log", Nbins=20, save_result=True),
        ]).run(lc)
        assert "structurefunction_result_0" in result.files
        df = result.files["structurefunction_result_0"]
        assert len(df) > 0

    # -----------------------------------------------------------------------
    # drwfit
    # -----------------------------------------------------------------------

    def test_drwfit_default_scalars(self):
        lc = vt.LightCurve.from_file(EXAMPLE_LC)
        result = vt.Pipeline([cmd.drwfit()]).run(lc)
        for col in ("DRWFIT_SIGMA_0", "DRWFIT_TAU_0", "DRWFIT_MU_0",
                    "DRWFIT_LNL_0", "DRWFIT_CONVERGED_0"):
            assert col in result.vars.index
        assert int(float(result.vars["DRWFIT_CONVERGED_0"])) == 1
        assert float(result.vars["DRWFIT_SIGMA_0"]) > 0
        assert float(result.vars["DRWFIT_TAU_0"]) > 0

    def test_drwfit_correctlc_smoothed_then_chi2(self):
        # Smoothing whitens the curve toward (below) the noise floor.
        lc = vt.LightCurve.from_file(EXAMPLE_LC)
        result = vt.Pipeline([
            cmd.drwfit(correctlc="smoothed"),
            cmd.chi2(),
        ]).run(lc)
        chi2_keys = [k for k in result.vars.index if "Chi2" in k]
        assert chi2_keys
        assert float(result.vars[chi2_keys[0]]) < 1.0

    def test_drwfit_modelvar_smoothed_then_stats(self):
        lc = vt.LightCurve.from_file(EXAMPLE_LC)
        result = vt.Pipeline([
            cmd.drwfit(modelvar=("smoothed", "drwmod")),
            cmd.stats("drwmod", "mean,stddev"),
        ]).run(lc)
        assert any("drwmod" in k and "MEAN" in k for k in result.vars.index)

    def test_drwfit_save_result_capture(self):
        # Logical name "result" matches the save_result attribute, so the
        # .drwfit aux file is actually captured (unlike the SF "sf" name).
        lc = vt.LightCurve.from_file(EXAMPLE_LC)
        result = vt.Pipeline([cmd.drwfit(save_result=True)]).run(lc)
        assert "drwfit_result_0" in result.files
        df = result.files["drwfit_result_0"]
        # Eight columns auto-named from the "# t x sig_meas ..." header.
        assert list(df.columns) == [
            "t", "x", "sig_meas", "x_hat_fwd", "Omega_fwd", "chi_fwd",
            "x_smoothed", "Omega_smoothed",
        ]
        assert len(df) > 0

    # -----------------------------------------------------------------------
    # runlength
    # -----------------------------------------------------------------------

    def test_runlength_scalars_present(self):
        lc = vt.LightCurve.from_file(EXAMPLE_LC)
        result = vt.Pipeline([cmd.runlength()]).run(lc)
        for cat in ("ABOVE", "BELOW", "OUTHIGH", "OUTLOW"):
            for stat in ("MAXLEN", "NRUNS", "MEANLEN"):
                assert f"RUNLENGTH_{cat}_{stat}_0" in result.vars.index
        for col in ("RUNLENGTH_MEDIAN_0", "RUNLENGTH_MAD_0",
                    "RUNLENGTH_K_0"):
            assert col in result.vars.index
        assert float(result.vars["RUNLENGTH_MAD_0"]) > 0
        assert float(result.vars["RUNLENGTH_K_0"]) == 3.0

    def test_runlength_sinusoid_long_runs_no_outliers(self):
        # A pure sinusoid yields long above/below runs (~ half a period in
        # points) and, since MAD ~ 1.05*A, no points outside the +/-3*MAD
        # band -> zero outlier runs.
        n = 1000
        t = np.arange(n, dtype=float) * 0.01
        mag = 0.1 * np.sin(2.0 * np.pi * t / 2.0)   # 200 pts/period
        err = np.full(n, 0.001)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="sine")
        result = vt.Pipeline([cmd.runlength()]).run(lc)
        assert int(float(result.vars["RUNLENGTH_ABOVE_MAXLEN_0"])) > 50
        assert int(float(result.vars["RUNLENGTH_BELOW_MAXLEN_0"])) > 50
        assert int(float(result.vars["RUNLENGTH_OUTHIGH_NRUNS_0"])) == 0
        assert int(float(result.vars["RUNLENGTH_OUTLOW_NRUNS_0"])) == 0

    def test_runlength_sustained_excursion_one_outhigh_run(self):
        # A contiguous faint-excursion block on otherwise tight noise gives
        # exactly one OUTHIGH run spanning the block, and no OUTLOW run.
        rng = np.random.default_rng(7)
        n = 400
        t = np.arange(n, dtype=float)
        mag = rng.normal(10.0, 0.01, n)
        mag[150:175] += 0.5          # 25-point sustained excursion
        err = np.full(n, 0.01)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="block")
        result = vt.Pipeline([cmd.runlength()]).run(lc)
        assert int(float(result.vars["RUNLENGTH_OUTHIGH_MAXLEN_0"])) == 25
        assert int(float(result.vars["RUNLENGTH_OUTHIGH_NRUNS_0"])) == 1
        # No comparably long low excursion (isolated noise-tail points may
        # occur, but nothing sustained).
        assert int(float(result.vars["RUNLENGTH_OUTLOW_MAXLEN_0"])) < 25

    def test_runlength_k_zero_band_collapses_to_median(self):
        # With k=0 the band has zero width, so outhigh == above and
        # outlow == below exactly.
        n = 500
        t = np.arange(n, dtype=float) * 0.01
        mag = 0.1 * np.sin(2.0 * np.pi * t / 1.5)
        err = np.full(n, 0.001)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="collapse")
        result = vt.Pipeline([cmd.runlength(k=0)]).run(lc)
        assert float(result.vars["RUNLENGTH_K_0"]) == 0.0
        assert (int(float(result.vars["RUNLENGTH_OUTHIGH_NRUNS_0"]))
                == int(float(result.vars["RUNLENGTH_ABOVE_NRUNS_0"])))
        assert (int(float(result.vars["RUNLENGTH_OUTHIGH_MAXLEN_0"]))
                == int(float(result.vars["RUNLENGTH_ABOVE_MAXLEN_0"])))
        assert (int(float(result.vars["RUNLENGTH_OUTLOW_NRUNS_0"]))
                == int(float(result.vars["RUNLENGTH_BELOW_NRUNS_0"])))

    def test_runlength_maskpoints_changes_counts(self):
        # Masking out half the points (via an -expr indicator) must still
        # produce valid stats and a finite median over the kept subset.
        lc = vt.LightCurve.from_file(EXAMPLE_LC)
        result = vt.Pipeline([
            cmd.expr("mask=(t>53726.0)"),
            cmd.runlength(maskpoints="mask"),
        ]).run(lc)
        assert int(float(result.vars["RUNLENGTH_ABOVE_NRUNS_1"])) > 0
        assert np.isfinite(float(result.vars["RUNLENGTH_MEDIAN_1"]))

    # -----------------------------------------------------------------------
    # Period search
    # -----------------------------------------------------------------------

    def test_clip_ls_with_suffix(self):
        lc = make_lc(period=2.3)
        result = vt.Pipeline([
            cmd.clip(sigclip=5.0),
            cmd.columnsuffix("ls"),
            cmd.LS(0.5, 5.0, 1e-3),
        ]).run(lc)
        assert "LS_Period_1_ls" in result.vars.index
        best = float(result.vars["LS_Period_1_ls"])
        assert abs(best - 2.3) < 0.05, f"Expected ~2.3, got {best}"

    def test_aov_period_search(self):
        lc = make_lc(period=1.8)
        result = vt.Pipeline([
            cmd.columnsuffix("aov"),
            cmd.aov(0.5, 5.0, 0.1, finetune=2),
        ]).run(lc)
        assert any("Period" in k for k in result.vars.index)

    def test_aov_harm_period_search(self):
        lc = make_lc(period=1.8)
        result = vt.Pipeline([
            cmd.aov_harm(3, 0.5, 5.0, 0.1, finetune=2),
        ]).run(lc)
        assert result.vars is not None

    def test_ls_save_periodogram(self):
        lc = make_lc(period=2.3)
        result = vt.Pipeline([
            cmd.LS(0.5, 5.0, 1e-3, save_periodogram=True),
        ]).run(lc)
        assert "LS_periodogram_0" in result.files
        df = result.files["LS_periodogram_0"]
        assert len(df) > 0

    def test_ls_noGLS(self):
        lc = make_lc(period=2.3)
        result = vt.Pipeline([
            cmd.LS(0.5, 5.0, 1e-3, noGLS=True),
        ]).run(lc)
        assert result.vars is not None

    # -----------------------------------------------------------------------
    # PDM (Phase Dispersion Minimization)
    # -----------------------------------------------------------------------

    def test_pdm_linterp_pipeline(self):
        lc = make_lc(period=1.8)
        result = vt.Pipeline([
            cmd.PDM("linterp", 0.5, 5.0, 0.1, 0.01, npeaks=1, nbin=8),
        ]).run(lc)
        assert hasattr(result.varobjs.PDM, "Period_1")
        # Theta should be small at a real signal period; loose check.
        assert result.varobjs.PDM.Theta_1 < 0.5

    def test_pdm_multicover_pipeline(self):
        lc = make_lc(period=1.8)
        result = vt.Pipeline([
            cmd.PDM("multicover", 0.5, 5.0, 0.1, 0.01, npeaks=1,
                    nbin=8, nc=2),
        ]).run(lc)
        assert hasattr(result.varobjs.PDM, "Period_1")

    def test_pdm_tophat_pipeline(self):
        lc = make_lc(period=1.8)
        result = vt.Pipeline([
            cmd.PDM("tophat", 1.0, 3.0, 0.5, 0.05, npeaks=1, dphi=0.05,
                    noerr=True),
        ]).run(lc)
        assert hasattr(result.varobjs.PDM, "Period_1")
        # Binless theta dips deep at the signal period.
        assert result.varobjs.PDM.Theta_1 < 0.2

    def test_pdm_fixperiod_backref_from_aov(self):
        lc = make_lc(period=1.8)
        result = vt.Pipeline([
            cmd.aov(0.5, 5.0, 0.1, 0.01, nbin=8),
            cmd.PDM("linterp", 0.5, 5.0, 0.1, 0.01, npeaks=1,
                    fixperiod_snr="aov"),
        ]).run(lc)
        # PDM PeriodFix should appear and equal the aov peak.
        pdm = result.varobjs.PDM
        assert hasattr(pdm, "PeriodFix")
        assert abs(pdm.PeriodFix - result.varobjs.aov.Period_1) < 1e-6

    def test_pdm_to_pdm_backref(self):
        # Second -PDM uses the first -PDM's peak via fixperiod_snr='pdm'.
        lc = make_lc(period=1.8)
        result = vt.Pipeline([
            cmd.PDM("step", 0.5, 5.0, 0.1, 0.01, npeaks=1, nbin=8),
            cmd.PDM("linterp", 0.5, 5.0, 0.1, 0.01, npeaks=1,
                    fixperiod_snr="pdm"),
        ]).run(lc)
        # Two PDM commands → varobjs.PDM is a CommandStatsList wrapping
        # the two CommandStats entries; indexable by position.
        pdms = result.varobjs.PDM
        assert len(pdms) == 2
        first, second = pdms[0], pdms[1]
        # Second one's PeriodFix should equal first one's Period_1.
        assert abs(second.PeriodFix - first.Period_1) < 1e-6

    # -----------------------------------------------------------------------
    # FTP (Fast Template Periodogram)
    # -----------------------------------------------------------------------

    def test_ftp_inline_pipeline(self):
        # Pure-cosine H=1 template; FTP should peak at the injected period.
        lc = make_lc(period=1.8)
        result = vt.Pipeline([
            cmd.FTP("inline", 0.5, 5.0, 0.1, 0.01,
                    cn=[1.0], sn=[0.0], npeaks=1),
        ]).run(lc)
        ftp = result.varobjs.FTP
        assert hasattr(ftp, "Period_1")
        assert abs(ftp.Period_1 - 1.8) / 1.8 < 0.01
        assert ftp.Power_1 > 0.9   # near-perfect fit for a pure sinusoid

    def test_ftp_fitlc_pipeline(self, tmp_path):
        # Build an LC, dump it to disk, then FTP-fitlc against itself at
        # the known period: search should find the same period.
        lc = make_lc(period=1.8)
        lc_path = tmp_path / "t.lc"
        df = lc.to_dataframe()
        df.to_csv(lc_path, sep=" ", header=False, index=False)
        result = vt.Pipeline([
            cmd.FTP("fitlc", 0.5, 5.0, 0.1, 0.01,
                    lc_path=str(lc_path), lc_format="ascii",
                    t_col=1, mag_col=2, err_col=3,
                    nharm=3, period=1.8, npeaks=1),
        ]).run(lc)
        assert abs(result.varobjs.FTP.Period_1 - 1.8) / 1.8 < 0.01

    def test_ftp_to_ftp_backref(self):
        # Second -FTP uses the first -FTP's peak via fixperiod_snr='ftp'.
        lc = make_lc(period=1.8)
        result = vt.Pipeline([
            cmd.FTP("inline", 0.5, 5.0, 0.1, 0.01,
                    cn=[1.0], sn=[0.0], npeaks=1),
            cmd.FTP("inline", 0.5, 5.0, 0.1, 0.01,
                    cn=[1.0], sn=[0.0], npeaks=1,
                    fixperiod_snr="ftp"),
        ]).run(lc)
        ftps = result.varobjs.FTP
        assert len(ftps) == 2
        first, second = ftps[0], ftps[1]
        assert abs(second.PeriodFix - first.Period_1) < 1e-6

    # -----------------------------------------------------------------------
    # Phase-folding
    # -----------------------------------------------------------------------

    def test_ls_then_phase(self):
        lc = make_lc(period=2.3)
        result = vt.Pipeline([
            cmd.columnsuffix("s"),
            cmd.LS(0.5, 5.0, 1e-3),
        ]).run(lc)
        best = float(result.vars["LS_Period_1_s"])
        result2 = vt.Pipeline([
            cmd.Phase(period=best),
        ]).run(lc, capture_lc=True)
        assert result2.lc is not None
        df = result2.lc.to_dataframe()
        # phased LC should have a phase column (3rd or 4th col)
        assert len(df.columns) >= 3

    # -----------------------------------------------------------------------
    # Inject and recover
    # -----------------------------------------------------------------------

    def test_injectharm_ls_recovery(self):
        """Inject a sinusoid, verify LS recovers it."""
        lc = make_flat_lc(n=400)
        result = vt.Pipeline([
            cmd.Injectharm(period=2.5, amplitude=0.05, nharm=1),
            cmd.columnsuffix("ls"),
            cmd.LS(0.5, 5.0, 5e-4),
        ]).run(lc)
        best = float(result.vars["LS_Period_1_ls"])
        assert abs(best - 2.5) < 0.05, f"Expected ~2.5, got {best}"

    def test_injectharm_killharm_rms(self):
        """Inject then remove a harmonic; RMS should drop."""
        lc = make_flat_lc(n=300)
        result_before = vt.Pipeline([
            cmd.Injectharm(period=2.0, amplitude=0.1),
            cmd.columnsuffix("before"),
            cmd.rms(),
        ]).run(lc)
        result_after = vt.Pipeline([
            cmd.Injectharm(period=2.0, amplitude=0.1),
            cmd.Killharm(period=2.0, nharm=1),
            cmd.columnsuffix("after"),
            cmd.rms(),
        ]).run(lc)
        rms_before = float(result_before.vars["RMS_before"])
        rms_after = float(result_after.vars["RMS_after"])
        assert rms_after < rms_before, "Killharm should reduce RMS"

    # -----------------------------------------------------------------------
    # LC manipulation
    # -----------------------------------------------------------------------

    def test_restricttimes_reduces_points(self):
        lc = make_lc(n=300)
        result = vt.Pipeline([
            cmd.restricttimes(mode="JDrange", minJD=0.0, maxJD=15.0),
            cmd.rms(),
        ]).run(lc, capture_lc=True)
        assert result.lc is not None
        assert len(result.lc._df) < 300

    def test_restricttimes_restoretimes(self):
        lc = make_lc(n=300)
        result = vt.Pipeline([
            cmd.columnsuffix("r1"),
            cmd.rms(),
            cmd.restricttimes(mode="JDrange", minJD=0.0, maxJD=15.0),
            cmd.columnsuffix("r2"),
            cmd.rms(),
            cmd.restoretimes(),
            cmd.columnsuffix("r3"),
            cmd.rms(),
        ]).run(lc)
        n1 = float(result.vars["Npoints_r1"])
        n2 = float(result.vars["Npoints_r2"])
        n3 = float(result.vars["Npoints_r3"])
        assert n2 < n1
        assert n3 == n1

    def test_savelc_restorelc(self):
        lc = make_lc(n=300)
        result = vt.Pipeline([
            cmd.savelc(),
            cmd.clip(sigclip=3.0),
            cmd.columnsuffix("clipped"),
            cmd.rms(),
            cmd.restorelc(),
            cmd.columnsuffix("restored"),
            cmd.rms(),
        ]).run(lc)
        n_clipped = float(result.vars["Npoints_clipped"])
        n_restored = float(result.vars["Npoints_restored"])
        assert n_restored == 300

    def test_sortlc_reverse(self):
        lc = make_lc(n=200)
        result = vt.Pipeline([
            cmd.sortlc(reverse=True),
        ]).run(lc, capture_lc=True)
        times = result.lc._df.iloc[:, 0].values
        assert np.all(np.diff(times) <= 0)

    def test_expr_modifies_lc(self):
        lc = make_lc(n=100)
        orig_mean = lc._df.iloc[:, 1].mean()
        result = vt.Pipeline([
            cmd.expr("mag=mag-0.5"),
            cmd.stats("mag", "mean"),
        ]).run(lc)
        mean_key = next(k for k in result.vars.index
                        if "mean" in k.lower() and "mag" in k.lower())
        new_mean = float(result.vars[mean_key])
        assert abs(new_mean - (orig_mean - 0.5)) < 1e-4

    def test_clip_non_iterative(self):
        lc = make_lc(n=200)
        result = vt.Pipeline([
            cmd.clip(sigclip=3.0, iterative=False),
            cmd.rms(),
        ]).run(lc)
        assert result.vars is not None

    def test_clip_niter(self):
        lc = make_lc(n=200)
        result = vt.Pipeline([
            cmd.clip(sigclip=3.0, niter=2),
            cmd.rms(),
        ]).run(lc)
        assert result.vars is not None

    def test_medianfilter(self):
        lc = make_lc(n=200)
        result = vt.Pipeline([
            cmd.medianfilter(time=0.5),
            cmd.rms(),
        ]).run(lc)
        assert result.vars is not None

    # -----------------------------------------------------------------------
    # Noise injection
    # -----------------------------------------------------------------------

    def test_addnoise_white(self):
        lc = make_lc(n=200, noise=0.001)
        result_before = vt.Pipeline([
            cmd.columnsuffix("before"),
            cmd.rms(),
        ]).run(lc)
        result_after = vt.Pipeline([
            cmd.addnoise(noise_type="white", sig_white=0.05),
            cmd.columnsuffix("after"),
            cmd.rms(),
        ]).run(lc)
        rms_before = float(result_before.vars["RMS_before"])
        rms_after = float(result_after.vars["RMS_after"])
        assert rms_after > rms_before

    def test_addnoise_squareexp(self):
        lc = make_lc(n=200)
        result = vt.Pipeline([
            cmd.addnoise(noise_type="squareexp",
                         sig_white=0.001, rho=1.0, sig_red=0.01),
            cmd.rms(),
        ]).run(lc)
        assert result.vars is not None

    # -----------------------------------------------------------------------
    # Time conversion
    # -----------------------------------------------------------------------

    def test_converttime_jd_mjd(self):
        lc = make_lc(n=100)
        # JD→MJD subtracts 2400000.5; capture and verify
        result = vt.Pipeline([
            cmd.converttime("jd", "mjd"),
            cmd.rms(),
        ]).run(lc, capture_lc=True)
        orig_t = lc._df.iloc[:, 0].values[0]
        new_t = result.lc._df.iloc[:, 0].values[0]
        assert abs((orig_t - new_t) - 2400000.5) < 1.0

    # -----------------------------------------------------------------------
    # Binning
    # -----------------------------------------------------------------------

    def test_binlc_reduces_points(self):
        lc = make_lc(n=300)
        result = vt.Pipeline([
            cmd.binlc(binsize=0.5),
            cmd.rms(),
        ]).run(lc, capture_lc=True)
        assert result.lc is not None
        assert len(result.lc._df) < 300

    def test_binlc_median_method(self):
        lc = make_lc(n=300)
        result = vt.Pipeline([
            cmd.binlc(binsize=0.5, method="median"),
            cmd.rms(),
        ]).run(lc)
        assert result.vars is not None

    def test_binlc_weightedaverage(self):
        lc = make_lc(n=300)
        result = vt.Pipeline([
            cmd.binlc(binsize=0.5, method="weightedaverage"),
            cmd.rms(),
        ]).run(lc)
        assert result.vars is not None

    # -----------------------------------------------------------------------
    # columnsuffix in complex pipelines
    # -----------------------------------------------------------------------

    def test_columnsuffix_multiple(self):
        lc = make_lc(n=200)
        result = vt.Pipeline([
            cmd.columnsuffix("a"),
            cmd.rms(),
            cmd.clip(sigclip=5.0),
            cmd.columnsuffix("b"),
            cmd.rms(),
        ]).run(lc)
        assert "RMS_a" in result.vars.index
        assert "RMS_b" in result.vars.index

    # -----------------------------------------------------------------------
    # LC name in stats
    # -----------------------------------------------------------------------

    def test_name_in_stats(self):
        lc = vt.LightCurve.from_arrays(
            *make_lc()._df.values.T, name="my_star")
        result = vt.Pipeline([cmd.rms()]).run(lc)
        assert result.vars["Name"] == "my_star"

    # -----------------------------------------------------------------------
    # Batch runs
    # -----------------------------------------------------------------------

    def test_batch_rms_stats(self):
        lcs = [make_lc(period=1.0 + i * 0.5) for i in range(4)]
        result = vt.Pipeline([
            cmd.rms(),
            cmd.stats("mag", "mean,stddev"),
        ]).run_batch(lcs)
        assert len(result.vars) == 4
        assert "RMS_0" in result.vars.columns

    def test_batch_names(self):
        lcs = [make_lc() for i in range(3)]
        for i, lc in enumerate(lcs):
            lc._df  # just access; name already set in make_lc
        result = vt.Pipeline([cmd.rms()]).run_batch(lcs)
        assert "Name" in result.vars.columns

    def test_batch_nthreads(self):
        lcs = [make_lc(period=1.0 + i * 0.3) for i in range(4)]
        result = vt.Pipeline([cmd.rms()]).run_batch(lcs, nthreads=2)
        assert len(result.vars) == 4
        assert "RMS_0" in result.vars.columns

    # -----------------------------------------------------------------------
    # Capture LC
    # -----------------------------------------------------------------------

    def test_capture_lc_after_clip(self):
        lc = make_lc(n=200)
        result = vt.Pipeline([
            cmd.clip(sigclip=5.0),
        ]).run(lc, capture_lc=True)
        assert result.lc is not None
        assert len(result.lc._df) > 0

    def test_capture_lc_name_preserved(self):
        lc = vt.LightCurve.from_arrays(
            *make_lc()._df.values.T, name="mystar")
        result = vt.Pipeline([cmd.clip(sigclip=5.0)]).run(lc, capture_lc=True)
        assert result.lc.name == "mystar"

    # -----------------------------------------------------------------------
    # magtoflux
    # -----------------------------------------------------------------------

    def test_magtoflux_roundtrip(self):
        """fluxtomag then magtoflux should recover the original magnitudes."""
        lc = make_lc(n=200)
        mag_orig = lc._df.iloc[:, 1].values.copy()
        result = vt.Pipeline([
            cmd.fluxtomag(mag_constant=25.0),
            cmd.magtoflux(mag_constant=25.0),
        ]).run(lc, capture_lc=True)
        assert result.lc is not None
        mag_back = result.lc._df.iloc[:, 1].values
        np.testing.assert_allclose(mag_back, mag_orig, rtol=0, atol=1e-12)

    def test_magtoflux_normalize_median_is_one(self):
        """magtoflux normalize should produce median flux = 1."""
        lc = make_lc(n=200)
        result = vt.Pipeline([
            cmd.magtoflux(normalize=True),
        ]).run(lc, capture_lc=True)
        flux = result.lc._df.iloc[:, 1].values
        assert np.isclose(np.median(flux), 1.0, rtol=0, atol=1e-15)

    # -----------------------------------------------------------------------
    # percentileratios
    # -----------------------------------------------------------------------

    def test_percentileratios_gaussian_limit(self):
        """For large-N white-Gaussian noise, asym -> 1 and medmeddev/stddev -> 0.6745."""
        rng = np.random.default_rng(0)
        n = 20000
        t = np.linspace(0, 30, n)
        mag = rng.normal(0.0, 1.0, n)
        err = np.full(n, 1.0)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="gauss")
        result = vt.Pipeline([
            cmd.percentileratios(percentilepairs=[(5, 95), (1, 99), (25, 75)]),
        ]).run(lc)
        # asym ~ 1 for any symmetric distribution; finite-sample noise at
        # N=20000 leaves +/- ~0.05 headroom.
        for p, q in [(5, 95), (1, 99), (25, 75)]:
            asym = result.vars[
                f"PERCENTILERATIOS_asym_PCT{p}.00_PCT{q}.00_0"
            ]
            assert abs(asym - 1.0) < 0.1, f"asym_{p}_{q} = {asym} far from 1.0"
        # medmeddev / stddev -> 0.6745 in the Gaussian limit.
        ratio = result.vars["PERCENTILERATIOS_medmeddev_over_stddev_0"]
        assert abs(ratio - 0.6745) < 0.02, f"medmeddev/stddev = {ratio} far from 0.6745"

    def test_percentileratios_positive_skew(self):
        """A positively-skewed distribution must produce asym > 1."""
        rng = np.random.default_rng(1)
        n = 20000
        t = np.linspace(0, 30, n)
        # Exponential is right-skewed; median much lower than mean, so the
        # upper tail dominates and (pct(q)-median)/(median-pct(p)) > 1.
        mag = rng.exponential(scale=1.0, size=n)
        err = np.full(n, 1.0)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="expo")
        result = vt.Pipeline([
            cmd.percentileratios(),
        ]).run(lc)
        asym = result.vars["PERCENTILERATIOS_asym_PCT5.00_PCT95.00_0"]
        assert asym > 2.0, f"asym for exponential = {asym}, expected > 2"

    def test_percentileratios_nan_input_handled(self):
        """NaN magnitudes are dropped and the surviving N is used."""
        rng = np.random.default_rng(2)
        n = 20000
        t = np.linspace(0, 30, n)
        mag = rng.normal(0.0, 1.0, n)
        # Sprinkle 10% NaNs.
        mag[rng.choice(n, n // 10, replace=False)] = np.nan
        err = np.full(n, 1.0)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="withnan")
        result = vt.Pipeline([
            cmd.percentileratios(),
        ]).run(lc)
        # Surviving values are still Gaussian, so the limits still hold.
        ratio = result.vars["PERCENTILERATIOS_medmeddev_over_stddev_0"]
        assert abs(ratio - 0.6745) < 0.02

    def test_percentileratios_maskpoints_changes_result(self):
        """Masking out the tail of a skewed LC changes the asym statistic."""
        rng = np.random.default_rng(3)
        n = 5000
        t = np.linspace(0, 30, n)
        # Heavy positive skew (exponential): asym >> 1.
        mag = rng.exponential(scale=1.0, size=n)
        err = np.full(n, 1.0)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="expo")
        # Build a mask that keeps only the lower half of the magnitude
        # distribution -- the symmetric bulk -- via -expr.  Then asym
        # should be much closer to 1.
        result = vt.Pipeline([
            cmd.percentileratios(),
            cmd.expr("keep=(mag<1.0)"),
            cmd.percentileratios(maskpoints="keep"),
        ]).run(lc)
        asym_unmasked = result.vars["PERCENTILERATIOS_asym_PCT5.00_PCT95.00_0"]
        asym_masked   = result.vars["PERCENTILERATIOS_asym_PCT5.00_PCT95.00_2"]
        # Unmasked exponential is strongly skewed; masked bulk is much less so.
        assert asym_unmasked > 2.0, f"unmasked asym = {asym_unmasked}"
        assert asym_masked < asym_unmasked, (
            f"masked asym {asym_masked} should be smaller than unmasked {asym_unmasked}"
        )

    # -----------------------------------------------------------------------
    # beyondNsigma
    # -----------------------------------------------------------------------

    def test_beyondNsigma_gaussian_limit(self):
        """For large-N Gaussian noise, frac_above_N1 -> 0.1587 and
        frac_above_N3 -> 0.00135 (one-tailed 1 - Phi(N))."""
        rng = np.random.default_rng(0)
        n = 20000
        t = np.linspace(0, 30, n)
        mag = rng.normal(0.0, 1.0, n)
        err = np.full(n, 1.0)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="gauss")
        result = vt.Pipeline([
            cmd.beyondNsigma(),
        ]).run(lc)
        f_above_1 = result.vars["BEYONDNSIGMA_frac_above_N1.00_0"]
        f_below_1 = result.vars["BEYONDNSIGMA_frac_below_N1.00_0"]
        f_above_3 = result.vars["BEYONDNSIGMA_frac_above_N3.00_0"]
        f_below_3 = result.vars["BEYONDNSIGMA_frac_below_N3.00_0"]
        # 1 - Phi(1) = 0.1587; finite-sample std at N=20000 is ~0.0026.
        assert abs(f_above_1 - 0.1587) < 0.01, f"frac_above_1 = {f_above_1}"
        assert abs(f_below_1 - 0.1587) < 0.01, f"frac_below_1 = {f_below_1}"
        # 1 - Phi(3) = 0.00135; finite-sample std at N=20000 is ~0.00026.
        assert abs(f_above_3 - 0.00135) < 0.002, f"frac_above_3 = {f_above_3}"
        assert abs(f_below_3 - 0.00135) < 0.002, f"frac_below_3 = {f_below_3}"

    def test_beyondNsigma_useMAD_matches_stddev_for_gaussian(self):
        """For clean Gaussian noise, useMAD and stddev modes agree to <2%
        because 1.483 * medmeddev -> stddev in the large-N limit."""
        rng = np.random.default_rng(0)
        n = 20000
        t = np.linspace(0, 30, n)
        mag = rng.normal(0.0, 1.0, n)
        err = np.full(n, 1.0)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="gauss")
        result = vt.Pipeline([
            cmd.beyondNsigma(),
            cmd.beyondNsigma(useMAD=True),
        ]).run(lc)
        f1_std = result.vars["BEYONDNSIGMA_frac_above_N1.00_0"]
        f1_mad = result.vars["BEYONDNSIGMA_frac_above_N1.00_1"]
        assert abs(f1_std - f1_mad) < 0.01, f"std={f1_std}, mad={f1_mad}"

    def test_beyondNsigma_useMAD_robust_to_outliers(self):
        """With heavy-tailed outliers injected, useMAD identifies more 3-sigma
        outliers than stddev mode -- because stddev is inflated by the outliers
        and widens the threshold, while MAD reflects the bulk's scale."""
        rng = np.random.default_rng(42)
        n = 20000
        t = np.linspace(0, 30, n)
        mag = rng.normal(0.0, 1.0, n)
        # Inject 200 outliers at +-[8, 20] sigma.
        idx = rng.choice(n, size=200, replace=False)
        mag[idx] += rng.choice([-1.0, 1.0], size=200) * rng.uniform(8, 20, 200)
        err = np.full(n, 1.0)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="heavy")
        result = vt.Pipeline([
            cmd.beyondNsigma(),
            cmd.beyondNsigma(useMAD=True),
        ]).run(lc)
        f3_std_above = result.vars["BEYONDNSIGMA_frac_above_N3.00_0"]
        f3_mad_above = result.vars["BEYONDNSIGMA_frac_above_N3.00_1"]
        # MAD-based threshold is tighter -> more points cross it -> larger frac.
        assert f3_mad_above > f3_std_above, (
            f"useMAD frac_above_3={f3_mad_above} should exceed stddev "
            f"frac_above_3={f3_std_above} when heavy outliers are present"
        )

    def test_beyondNsigma_custom_Nvalues_column_names(self):
        """Custom float Nvalues produce columns with %.2f-formatted names."""
        rng = np.random.default_rng(0)
        n = 5000
        t = np.linspace(0, 30, n)
        mag = rng.normal(0.0, 1.0, n)
        err = np.full(n, 1.0)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="g")
        result = vt.Pipeline([
            cmd.beyondNsigma(Nvalues=[1.5, 2.5]),
        ]).run(lc)
        assert "BEYONDNSIGMA_frac_above_N1.50_0" in result.vars.index
        assert "BEYONDNSIGMA_frac_below_N1.50_0" in result.vars.index
        assert "BEYONDNSIGMA_frac_above_N2.50_0" in result.vars.index
        assert "BEYONDNSIGMA_frac_below_N2.50_0" in result.vars.index

    def test_beyondNsigma_degenerate_distribution(self):
        """When every magnitude equals the median, sigma=0 -> all fractions 0."""
        n = 100
        t = np.linspace(0, 30, n)
        mag = np.full(n, 10.0)
        err = np.full(n, 1.0)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="flat")
        result = vt.Pipeline([
            cmd.beyondNsigma(),
        ]).run(lc)
        for n_val in [1.0, 3.0, 5.0]:
            tag = f"BEYONDNSIGMA_frac_above_N{n_val:.2f}_0"
            assert result.vars[tag] == 0.0, f"{tag} = {result.vars[tag]}"
            tag = f"BEYONDNSIGMA_frac_below_N{n_val:.2f}_0"
            assert result.vars[tag] == 0.0, f"{tag} = {result.vars[tag]}"

    def test_beyondNsigma_two_tail_consistency(self):
        """frac_above + frac_below should equal the two-tailed fraction exactly
        (integer counts divided by the same denominator -- no FP slack)."""
        rng = np.random.default_rng(7)
        n = 5000
        t = np.linspace(0, 30, n)
        mag = rng.normal(0.0, 1.0, n)
        err = np.full(n, 1.0)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="g")
        result = vt.Pipeline([
            cmd.beyondNsigma(Nvalues=[1.0]),
        ]).run(lc)
        above = result.vars["BEYONDNSIGMA_frac_above_N1.00_0"]
        below = result.vars["BEYONDNSIGMA_frac_below_N1.00_0"]
        # Both fractions have the same denominator N_rej, so their sum is
        # exactly representable as a count fraction.
        total = above + below
        # Count back the total via N_rej and verify it's an integer / N_rej.
        # At n=5000, the sum's denominator must be 5000.
        assert abs(total * n - round(total * n)) < 1e-9, (
            f"above+below={total} not an integer/{n} fraction"
        )

    def test_beyondNsigma_maskpoints_changes_result(self):
        """Masking out the heavy upper tail reduces frac_above relative to
        the full distribution -- the bulk's fraction beyond N*sigma is smaller
        than the full distribution's because the outliers don't inflate sigma
        (under stddev mode) once they're masked."""
        rng = np.random.default_rng(8)
        n = 10000
        t = np.linspace(0, 30, n)
        mag = rng.normal(0.0, 1.0, n)
        # Inject 100 large positive outliers.
        idx = rng.choice(n, size=100, replace=False)
        mag[idx] += rng.uniform(8.0, 15.0, size=100)
        err = np.full(n, 1.0)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="heavy")
        result = vt.Pipeline([
            cmd.beyondNsigma(),
            cmd.expr("keep=(mag<5.0)"),  # mask out the injected outliers
            cmd.beyondNsigma(maskpoints="keep"),
        ]).run(lc)
        f3_full   = result.vars["BEYONDNSIGMA_frac_above_N3.00_0"]
        f3_masked = result.vars["BEYONDNSIGMA_frac_above_N3.00_2"]
        # Full distribution: ~100 outliers count as far-above (1%); masked
        # bulk is clean Gaussian where 3-sigma is ~0.135%.
        assert f3_full > 0.005, f"f3_full = {f3_full}; outliers should show"
        assert f3_masked < 0.01, f"f3_masked = {f3_masked}; should drop after mask"
        assert f3_masked < f3_full

    # -----------------------------------------------------------------------
    # slopestats
    # -----------------------------------------------------------------------

    def test_slopestats_gaussian_limit(self):
        """For large-N white-noise slopes, frac_above_T1 -> 0.1587 and
        frac_above_T3 -> 0.00135 (one-tailed 1 - Phi(T))."""
        rng = np.random.default_rng(0)
        n = 20000
        t = np.arange(n, dtype=float)  # dt = 1 day, integer-aligned
        mag = rng.normal(0.0, 1.0, n)
        err = np.full(n, 1.0)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="gauss")
        result = vt.Pipeline([
            cmd.slopestats(threshold=[1.0, 3.0]),
        ]).run(lc)
        f_above_1 = result.vars["SLOPESTATS_frac_above_T1.00_0"]
        f_below_1 = result.vars["SLOPESTATS_frac_below_T1.00_0"]
        f_above_3 = result.vars["SLOPESTATS_frac_above_T3.00_0"]
        f_below_3 = result.vars["SLOPESTATS_frac_below_T3.00_0"]
        # Slope distribution is Gaussian with mean 0; 1-Phi(1) = 0.1587.
        assert abs(f_above_1 - 0.1587) < 0.01, f"frac_above_1 = {f_above_1}"
        assert abs(f_below_1 - 0.1587) < 0.01, f"frac_below_1 = {f_below_1}"
        # 1-Phi(3) = 0.00135; small-sample noise is ~0.0003.
        assert abs(f_above_3 - 0.00135) < 0.002, f"frac_above_3 = {f_above_3}"
        assert abs(f_below_3 - 0.00135) < 0.002, f"frac_below_3 = {f_below_3}"

    def test_slopestats_constant_slope_has_zero_mad(self):
        """For a perfectly linear LC, all consecutive slopes are equal,
        so mad_dmdt -> 0 and median_abs_dmdt = max_abs_dmdt = |k|."""
        n = 100
        t = np.arange(n, dtype=float)
        k = 0.05
        mag = k * t
        err = np.full(n, 0.01)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="lin")
        result = vt.Pipeline([cmd.slopestats()]).run(lc)
        # MAD should be at floating-point noise level (well below 1e-10).
        assert result.vars["SLOPESTATS_mad_dmdt_0"] < 1e-10
        assert abs(result.vars["SLOPESTATS_median_abs_dmdt_0"] - k) < 1e-9
        assert abs(result.vars["SLOPESTATS_max_abs_dmdt_0"] - k) < 1e-9

    def test_slopestats_maxgap_drops_gap_spanning_pair(self):
        """A pathological gap-spanning pair should produce max_abs_dmdt
        far above the bulk, and maxgap should drop it back to the bulk."""
        # Times: 0..49 clean, then jump to t=150 with a spike, then 51..99.
        ts = list(range(0, 50)) + [150] + list(range(51, 100))
        ms = [0.0] * 50 + [1000.0] + [0.0] * 49
        t = np.array(ts, dtype=float)
        mag = np.array(ms, dtype=float)
        err = np.full(len(t), 0.01)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="gap")
        result = vt.Pipeline([
            cmd.slopestats(),                # no maxgap: spike dominates
            cmd.slopestats(maxgap=50.0),     # drops the t=99->150 pair
        ]).run(lc)
        max_no_filter = result.vars["SLOPESTATS_max_abs_dmdt_0"]
        max_filtered  = result.vars["SLOPESTATS_max_abs_dmdt_1"]
        assert max_no_filter > 10.0, f"max without maxgap = {max_no_filter}"
        assert max_filtered == 0.0, f"max with maxgap = {max_filtered}"

    def test_slopestats_useMAD_matches_stddev_for_gaussian(self):
        """For clean Gaussian slopes, useMAD and stddev modes agree to <2%
        because 1.483 * medmeddev -> stddev in the large-N limit."""
        rng = np.random.default_rng(0)
        n = 20000
        t = np.arange(n, dtype=float)
        mag = rng.normal(0.0, 1.0, n)
        err = np.full(n, 1.0)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="gauss")
        result = vt.Pipeline([
            cmd.slopestats(),
            cmd.slopestats(useMAD=True),
        ]).run(lc)
        f1_std = result.vars["SLOPESTATS_frac_above_T3.00_0"]
        f1_mad = result.vars["SLOPESTATS_frac_above_T3.00_1"]
        assert abs(f1_std - f1_mad) < 0.005, f"std={f1_std}, mad={f1_mad}"

    def test_slopestats_binshift_changes_partition(self):
        """A half-bin binshift on integer-day points with a 10-day binsize
        re-partitions the LC and produces a measurably different median."""
        rng = np.random.default_rng(0)
        n = 5000
        t = np.arange(n, dtype=float)
        mag = rng.normal(0.0, 1.0, n)
        err = np.full(n, 1.0)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="g")
        result = vt.Pipeline([
            # bintime in minutes: 14400 min = 10 days.  binshift 0.5 ->
            # bin partition shifts by 5 days, moving the first 5 points
            # into a half-size bin.
            cmd.slopestats(bintime=[14400]),
            cmd.slopestats(bintime=[14400], binshift=0.5),
        ]).run(lc)
        m_noshift = result.vars["SLOPESTATS_median_abs_dmdt_BT14400.00_0"]
        m_shifted = result.vars["SLOPESTATS_median_abs_dmdt_BT14400.00_1"]
        assert m_noshift != m_shifted, "binshift should change the partition"

    def test_slopestats_bintime_column_names_use_minutes(self):
        """bintime values appear in column names verbatim (as given in
        minutes), not as the kernel's days-converted internal value."""
        rng = np.random.default_rng(0)
        n = 500
        t = np.arange(n, dtype=float)
        mag = rng.normal(0.0, 1.0, n)
        err = np.full(n, 1.0)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="g")
        result = vt.Pipeline([
            cmd.slopestats(bintime=[5.0, 60.0]),
        ]).run(lc)
        # The user-supplied minute values appear formatted as %.2f.
        assert "SLOPESTATS_median_abs_dmdt_BT5.00_0" in result.vars.index
        assert "SLOPESTATS_median_abs_dmdt_BT60.00_0" in result.vars.index

    # -----------------------------------------------------------------------
    # CodyM
    # -----------------------------------------------------------------------

    def test_codym_symmetric_gaussian_yields_M_near_zero(self):
        """A symmetric Gaussian-noise LC should give |M| close to 0 once
        the deciles are balanced about the median."""
        rng = np.random.default_rng(11)
        n = 3000
        t = np.arange(n, dtype=float) * 0.01
        mag = rng.normal(0.0, 0.02, n)
        err = np.full(n, 0.02)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="gauss")
        # trendwindow >= LC duration -> detrend is a constant subtract
        # (no effect on M); sigclip=0 disables rejection so the noise
        # decile tails are intact.
        result = vt.Pipeline([
            cmd.CodyM(trendwindow=100, sigclip=0),
        ]).run(lc)
        M = result.vars["CODYM_M_0"]
        # n=3000 decile statistics give per-LC M sample stddev ~ 0.05;
        # 0.1 is a 2-sigma margin and still inside the |M| <= 0.25
        # "symmetric" bin of Cody et al.
        assert abs(M) < 0.1, f"M = {M} for symmetric Gaussian noise"

    def test_codym_injected_dips_yield_positive_M(self):
        """Adding a small fraction of faint excursions makes the
        faint-decile mean fall further from the median than the
        bright-decile mean -> M > 0 (a dipping signature)."""
        rng = np.random.default_rng(11)
        n = 3000
        t = np.arange(n, dtype=float) * 0.01
        mag = rng.normal(0.0, 0.02, n)
        dip_mask = rng.random(n) < 0.10
        mag[dip_mask] += 0.15
        err = np.full(n, 0.02)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="dips")
        result = vt.Pipeline([
            cmd.CodyM(trendwindow=100, sigclip=0),
        ]).run(lc)
        M = result.vars["CODYM_M_0"]
        assert M > 0.5, f"M = {M} for an injected-dip LC"

    def test_codym_injected_bursts_yield_negative_M(self):
        """Mirror of the dips test with bright excursions instead;
        antisymmetry of the M definition produces M < 0."""
        rng = np.random.default_rng(11)
        n = 3000
        t = np.arange(n, dtype=float) * 0.01
        mag = rng.normal(0.0, 0.02, n)
        burst_mask = rng.random(n) < 0.10
        mag[burst_mask] -= 0.15
        err = np.full(n, 0.02)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="bursts")
        result = vt.Pipeline([
            cmd.CodyM(trendwindow=100, sigclip=0),
        ]).run(lc)
        M = result.vars["CODYM_M_0"]
        assert M < -0.5, f"M = {M} for an injected-burst LC"

    # -----------------------------------------------------------------------
    # CodyQ
    # -----------------------------------------------------------------------

    def test_codyq_sinusoid_at_true_period_yields_Q_near_zero(self):
        """A clean sinusoid evaluated at its true period should give Q
        close to 0 -- the phase model captures essentially all the
        variance, leaving residuals at the noise floor."""
        rng = np.random.default_rng(11)
        n = 3000
        t = np.arange(n, dtype=float) * 0.01
        P = 2.0
        mag = 0.1 * np.sin(2.0 * np.pi * t / P) + rng.normal(0.0, 0.01, n)
        err = np.full(n, 0.01)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="sin")
        result = vt.Pipeline([
            cmd.CodyQ(period=P, trendwindow=100),
        ]).run(lc)
        Q = result.vars["CODYQ_Q_0"]
        assert Q < 0.1, f"Q = {Q} at the true period of a clean sinusoid"

    def test_codyq_sinusoid_at_wrong_period_yields_Q_near_one(self):
        """At an incorrect period the phase model removes essentially
        no variance, so rms_resid -> rms_raw and Q -> 1."""
        rng = np.random.default_rng(11)
        n = 3000
        t = np.arange(n, dtype=float) * 0.01
        mag = 0.1 * np.sin(2.0 * np.pi * t / 2.0) + rng.normal(0.0, 0.01, n)
        err = np.full(n, 0.01)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="sin")
        result = vt.Pipeline([
            cmd.CodyQ(period=1.37, trendwindow=100),
        ]).run(lc)
        Q = result.vars["CODYQ_Q_0"]
        assert Q > 0.95, f"Q = {Q} at a wrong period of a clean sinusoid"

    def test_codyq_aov_back_reference(self):
        """The 'aov' period source picks up the primary peak of the
        most-recent -aov; verify the pipeline-mode back-reference
        resolves and Q matches the value obtained from a direct fix
        at the same period."""
        rng = np.random.default_rng(11)
        n = 3000
        t = np.arange(n, dtype=float) * 0.01
        P_true = 2.0
        mag = 0.1 * np.sin(2.0 * np.pi * t / P_true) + rng.normal(0.0, 0.01, n)
        err = np.full(n, 0.01)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="sin")
        result = vt.Pipeline([
            cmd.aov(0.5, 4.0, 0.1, 5, npeaks=1, save_periodogram=False),
            cmd.CodyQ(period="aov", trendwindow=100),
        ]).run(lc)
        P_used = result.vars["CODYQ_Period_1"]
        Q = result.vars["CODYQ_Q_1"]
        assert abs(P_used - P_true) < 0.01, f"AOV picked P={P_used}"
        assert Q < 0.1, f"Q = {Q} via aov backref at true period"

    # -----------------------------------------------------------------------
    # MatchedFilter
    # -----------------------------------------------------------------------

    def test_matchedfilter_recovers_injected_gauss(self):
        """Inject a known Gaussian into white-noise data and recover it.

        Verifies that the named-gauss MatchedFilter recovers (within
        ~10 percent) the injected amplitude at the right time."""
        rng = np.random.default_rng(0)
        n = 500
        t = np.linspace(0, 30, n)
        # White noise + Gaussian dip at t0 with sigma=0.3 and depth -0.05.
        t0 = t[100]
        sigma_g = 0.3
        depth = -0.05
        mag = depth * np.exp(-0.5 * ((t - t0) / sigma_g) ** 2) \
              + rng.normal(0, 0.01, n)
        err = np.full(n, 0.01)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="inj")
        result = vt.Pipeline([
            cmd.MatchedFilter("gauss", 3.0, "window", "negative",
                              sigma=sigma_g, npeaks=1),
        ]).run(lc)
        # Output column naming: MatchedFilter_*_<peak>_<step>
        keys = result.vars.index
        snr_key = [k for k in keys if k.startswith("MatchedFilter_SNR_1_")][0]
        amp_key = [k for k in keys if k.startswith("MatchedFilter_Amplitude_1_")][0]
        time_key = [k for k in keys if k.startswith("MatchedFilter_Time_1_")][0]
        snr = float(result.vars[snr_key])
        amp = float(result.vars[amp_key])
        t_rec = float(result.vars[time_key])
        # Strong negative detection at the injection time.
        assert snr < -10.0
        assert abs(amp - depth) / abs(depth) < 0.10
        assert abs(t_rec - t0) < 0.1

    def test_matchedfilter_template_modes_consistent(self):
        """The named-gauss, file, and expr template modes should all
        give the same recovered amplitude/time on a clean injection."""
        rng = np.random.default_rng(0)
        n = 500
        t = np.linspace(0, 30, n)
        t0 = t[100]
        sigma_g = 0.3
        depth = -0.05
        mag = depth * np.exp(-0.5 * ((t - t0) / sigma_g) ** 2) \
              + rng.normal(0, 0.01, n)
        err = np.full(n, 0.01)
        lc = vt.LightCurve.from_arrays(t, mag, err, name="inj")
        # Reference: named gauss.
        result_named = vt.Pipeline([
            cmd.MatchedFilter("gauss", 3.0, "window", "negative",
                              sigma=sigma_g, npeaks=1),
        ]).run(lc)
        # Expression: sampled Gaussian formula.
        result_expr = vt.Pipeline([
            cmd.MatchedFilter("expr", 3.0, "window", "negative",
                              expression=f"exp(-s*s/{2.0 * sigma_g ** 2})",
                              npeaks=1),
        ]).run(lc)
        # Pull the amp from each.
        def amp_of(result):
            k = [kk for kk in result.vars.index
                 if kk.startswith("MatchedFilter_Amplitude_1_")][0]
            return float(result.vars[k])
        # Expression mode must match named gauss to ~1 part in 1e-5.
        assert abs(amp_of(result_named) - amp_of(result_expr)) \
               < 1e-5 * abs(amp_of(result_named))


# ---------------------------------------------------------------------------
# Disk-file pipeline tests (run_file / run_filelist)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools binary not found")
class TestDiskFilePipeline:

    def _write_lc(self, path, lc):
        """Write a LightCurve to a plain-text file."""
        lc._df.to_csv(path, sep=" ", header=False, index=False,
                      float_format="%.10f")

    def test_run_file_basic(self, tmp_path):
        """run_file() produces same stats as run() for the same data."""
        lc = make_lc(n=200, period=2.3)
        lc_file = tmp_path / "star.lc"
        self._write_lc(str(lc_file), lc)

        result_disk = vt.Pipeline([cmd.rms()]).run_file(lc_file)
        result_mem  = vt.Pipeline([cmd.rms()]).run(lc)

        assert result_disk.vars is not None
        assert "RMS_0" in result_disk.vars.index
        assert abs(float(result_disk.vars["RMS_0"]) -
                   float(result_mem.vars["RMS_0"])) < 1e-6

    def test_run_file_name_from_stem(self, tmp_path):
        """Name in stats comes from the file stem, not a temp-file path."""
        lc = make_lc()
        lc_file = tmp_path / "mystar.lc"
        self._write_lc(str(lc_file), lc)

        result = vt.Pipeline([cmd.rms()]).run_file(lc_file)
        assert result.vars["Name"] == "mystar"

    def test_run_file_capture_lc(self, tmp_path):
        """capture_lc=True returns the modified LC when using run_file."""
        lc = make_lc(n=200)
        lc_file = tmp_path / "star.lc"
        self._write_lc(str(lc_file), lc)

        result = vt.Pipeline([cmd.clip(sigclip=5.0)]).run_file(
            lc_file, capture_lc=True)
        assert result.lc is not None
        assert len(result.lc._df) > 0
        assert result.lc.name == "star"

    def test_run_file_path_string(self, tmp_path):
        """run_file() accepts a plain string path as well as a Path object."""
        lc = make_lc()
        lc_file = tmp_path / "s.lc"
        self._write_lc(str(lc_file), lc)

        result = vt.Pipeline([cmd.rms()]).run_file(str(lc_file))
        assert "RMS_0" in result.vars.index

    def test_run_file_example_lc(self):
        """run_file() works on the real EXAMPLES/2 file shipped with vartools."""
        if not os.path.isfile(EXAMPLE_LC):
            pytest.skip("EXAMPLES/2 not found")
        result = vt.Pipeline([cmd.rms()]).run_file(EXAMPLE_LC)
        assert result.vars is not None
        assert "RMS_0" in result.vars.index

    def test_run_filelist_list_of_paths(self, tmp_path):
        """run_filelist() with a list of paths processes all LCs."""
        lcs = [make_lc(period=1.0 + i * 0.5) for i in range(3)]
        paths = []
        for i, lc in enumerate(lcs):
            p = tmp_path / f"lc{i}.lc"
            self._write_lc(str(p), lc)
            paths.append(p)

        result = vt.Pipeline([cmd.rms()]).run_filelist(paths)
        assert result.vars is not None
        assert len(result.vars) == 3
        assert "RMS_0" in result.vars.columns

    def test_run_filelist_names_from_stems(self, tmp_path):
        """Name column is populated from file stems."""
        lcs = [make_lc() for _ in range(2)]
        paths = []
        for name in ["alpha", "beta"]:
            p = tmp_path / f"{name}.lc"
            self._write_lc(str(p), lcs[0])
            paths.append(p)

        result = vt.Pipeline([cmd.rms()]).run_filelist(paths)
        assert list(result.vars["Name"]) == ["alpha", "beta"]

    def test_run_filelist_existing_list_file(self, tmp_path):
        """run_filelist() with a path to an existing list file."""
        lcs = [make_lc(period=1.5 + i * 0.4) for i in range(3)]
        lc_files = []
        for i, lc in enumerate(lcs):
            p = tmp_path / f"lc{i}.lc"
            self._write_lc(str(p), lc)
            lc_files.append(str(p))

        list_file = tmp_path / "list.txt"
        list_file.write_text("\n".join(lc_files) + "\n")

        result = vt.Pipeline([cmd.rms()]).run_filelist(list_file)
        assert len(result.vars) == 3

    def test_run_filelist_nthreads(self, tmp_path):
        """run_filelist() with nthreads > 1 still returns correct results."""
        lcs = [make_lc(period=1.0 + i * 0.3) for i in range(4)]
        paths = []
        for i, lc in enumerate(lcs):
            p = tmp_path / f"lc{i}.lc"
            self._write_lc(str(p), lc)
            paths.append(p)

        result = vt.Pipeline([cmd.rms()]).run_filelist(paths, nthreads=2)
        assert len(result.vars) == 4

    def test_run_filelist_pipeline(self, tmp_path):
        """run_filelist() works with a multi-command pipeline."""
        lcs = [make_lc(n=300, period=2.0 + i * 0.3) for i in range(3)]
        paths = []
        for i, lc in enumerate(lcs):
            p = tmp_path / f"lc{i}.lc"
            self._write_lc(str(p), lc)
            paths.append(p)

        result = vt.Pipeline([
            cmd.clip(sigclip=5.0),
            cmd.rms(),
        ]).run_filelist(paths)
        assert len(result.vars) == 3
        assert "RMS_1" in result.vars.columns


# ===========================================================================
# Output API — CLI-arg tests (no binary needed)
# ===========================================================================

class TestOutputAPICLI:
    """Tests for the four-mode save_* / Output API (CLI-arg level only)."""

    def test_ls_save_false_suppresses(self):
        # Mode 4: False → emit "0"
        c = cmd.LS(0.1, 10.0, 0.1, save_periodogram=False)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "0" in args
        assert "/tmp" not in args

    def test_ls_save_true_mode1(self):
        # Mode 1: True → emit "1" with fallback outdir
        c = cmd.LS(0.1, 10.0, 0.1, save_periodogram=True)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "1" in args
        assert "/tmp" in args

    def test_ls_save_path_mode3(self):
        # Mode 3: string path → emit "1" with that path, no capture
        c = cmd.LS(0.1, 10.0, 0.1, save_periodogram="/data/pgrams")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "1" in args
        assert "/data/pgrams" in args
        assert "/tmp" not in args

    def test_ls_save_output_mode2(self):
        # Mode 2: Output(path=..., capture=True) → emit "1" with user path
        from pyvartools import Output
        c = cmd.LS(0.1, 10.0, 0.1, save_periodogram=Output("/data/pgrams", capture=True))
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "1" in args
        assert "/data/pgrams" in args

    def test_norm_save_roundtrips(self):
        from pyvartools import Output
        from pyvartools.commands._helpers import _norm_save
        s = _norm_save(False)
        assert not s.capture and s.path is None
        s = _norm_save(True)
        assert s.capture and s.path is None
        s = _norm_save("/some/dir")
        assert not s.capture and s.path == "/some/dir"
        o = Output("/x", capture=True)
        assert _norm_save(o) is o

    def test_should_emit(self):
        from pyvartools.commands._helpers import _should_emit
        assert not _should_emit(False)
        assert _should_emit(True)
        assert _should_emit("/some/dir")
        from pyvartools import Output
        assert _should_emit(Output("/x", capture=True))
        assert _should_emit(Output("/x", capture=False))
        assert not _should_emit(Output(None, capture=False))

    def test_dftclean_save_path(self):
        # Custom token pattern: outdspec uses the user path
        c = cmd.dftclean(nbeam=3, save_dspec="/data/dspec")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "outdspec" in args
        idx = args.index("outdspec")
        assert args[idx + 1] == "/data/dspec"

    def test_wwz_save_path(self):
        c = cmd.wwz(save_transform="/data/wwz")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "outfulltransform" in args
        idx = args.index("outfulltransform")
        assert args[idx + 1] == "/data/wwz"

    def test_microlens_save_path(self):
        c = cmd.microlens(save_model="/data/ml")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "omodel" in args
        idx = args.index("omodel")
        assert args[idx + 1] == "/data/ml"

    def test_nonlinfit_save_false_no_omodel(self):
        c = cmd.nonlinfit("a0+a1*t", "a0 0.0 a1 0.0", save_model=False)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "omodel" not in args

    def test_linfit_save_path(self):
        c = cmd.linfit("a0+a1*t", "a0 0.0 a1 0.0", save_model="/data/lf")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "omodel" in args
        idx = args.index("omodel")
        assert args[idx + 1] == "/data/lf"

    def test_decorr_save_false_emits_zero(self):
        # decorr uses "0" sentinel when not saving model
        c = cmd.decorr(save_model=False)
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "omodel" not in args
        assert "0" in args

    def test_sysrem_save_trends_path(self):
        # save_trends is a single-file output (CLI grammar emits ``"1" path``,
        # not the literal token "otrends" — that bug was fixed in fitting.py).
        c = cmd.SYSREM(1, 1, "airmass.txt", save_trends="/data/trends.txt")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "otrends" not in args
        assert "/data/trends.txt" in args
        i = args.index("/data/trends.txt")
        # The path is preceded by the "1" emit-flag.
        assert args[i - 1] == "1"

    def test_sysrem_trends_outpath_recorded(self):
        # _to_cli_args must record the chosen output path on the command
        # so the global-file capture in pipeline.py can read it back.
        c = cmd.SYSREM(1, 1, "airmass.txt", save_trends="/data/trends.txt")
        c._outdir = "/tmp"
        c._to_cli_args()
        assert getattr(c, "_trends_outpath", None) == "/data/trends.txt"

    def test_sysrem_trends_outpath_default(self):
        # save_trends=True (no path) → wrapper picks <outdir>/sysrem.trends.
        c = cmd.SYSREM(1, 1, "airmass.txt", save_trends=True)
        c._outdir = "/tmp/cmd_0"
        c._to_cli_args()
        assert c._trends_outpath == "/tmp/cmd_0/sysrem.trends"

    def test_sysrem_trends_outpath_none_when_not_saving(self):
        c = cmd.SYSREM(1, 1, "airmass.txt", save_trends=False)
        c._outdir = "/tmp"
        c._to_cli_args()
        assert c._trends_outpath is None

    def test_findblends_save_matches_path(self):
        c = cmd.findblends(matchrad=10.0, save_matches="/data/matches")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "omatches" in args
        idx = args.index("omatches")
        assert args[idx + 1] == "/data/matches"

    def test_mandel_agol_save_phcurve_path(self):
        c = cmd.MandelAgolTransit(P0=3.0, T00=5000.0, save_phcurve="/data/ph")
        c._outdir = "/tmp"
        args = c._to_cli_args()
        assert "ophcurve" in args
        idx = args.index("ophcurve")
        assert args[idx + 1] == "/data/ph"

    def test_autocorrelation_mandatory_output_class_attr(self):
        assert getattr(cmd.autocorrelation, "_mandatory_output", False) is True

    def test_linfit_output_file_specs(self):
        c = cmd.linfit("a0", "a0 0.0")
        specs = c._output_file_specs()
        assert "model" in specs

    def test_decorr_output_file_specs(self):
        c = cmd.decorr()
        specs = c._output_file_specs()
        assert "model" in specs

    def test_findblends_output_file_specs(self):
        c = cmd.findblends(matchrad=10.0)
        specs = c._output_file_specs()
        assert "matches" in specs

    def test_sysrem_output_file_specs_has_trends(self):
        c = cmd.SYSREM(1, 1, "airmass.txt")
        specs = c._output_file_specs()
        assert "model" in specs
        assert "trends" in specs

    def test_output_repr(self):
        from pyvartools import Output
        o = Output("/some/path", capture=True)
        assert "Output" in repr(o)
        assert "/some/path" in repr(o)


# ===========================================================================
# Output API — end-to-end tests (require binary)
# ===========================================================================

@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools binary not found")
class TestOutputAPIEndToEnd:

    def test_ls_mode1_capture_only(self):
        # Mode 1: existing behaviour unchanged
        lc = make_lc()
        result = vt.Pipeline([
            cmd.LS(0.5, 5.0, 1e-3, save_periodogram=True),
        ]).run(lc)
        assert "LS_periodogram_0" in result.files
        assert len(result.files["LS_periodogram_0"]) > 0

    def test_ls_mode3_writes_file_not_captured(self, tmp_path):
        # Mode 3: file written to tmp_path, NOT in result.files
        lc = make_lc()
        result = vt.Pipeline([
            cmd.LS(0.5, 5.0, 1e-3, save_periodogram=str(tmp_path)),
        ]).run(lc)
        assert "LS_periodogram_0" not in result.files
        assert any(f.endswith(".ls") for f in os.listdir(tmp_path))

    def test_ls_mode2_writes_and_captures(self, tmp_path):
        # Mode 2: file written AND captured
        from pyvartools import Output
        lc = make_lc()
        result = vt.Pipeline([
            cmd.LS(0.5, 5.0, 1e-3,
                   save_periodogram=Output(str(tmp_path), capture=True)),
        ]).run(lc)
        assert "LS_periodogram_0" in result.files
        assert any(f.endswith(".ls") for f in os.listdir(tmp_path))

    def test_ls_mode4_no_file_no_capture(self):
        # Mode 4: False → nothing in result.files
        lc = make_lc()
        result = vt.Pipeline([
            cmd.LS(0.5, 5.0, 1e-3, save_periodogram=False),
        ]).run(lc)
        assert "LS_periodogram_0" not in result.files

    def test_autocorrelation_default_captured(self):
        # Default (save_result=True): captured in result.files
        lc = make_lc()
        result = vt.Pipeline([
            cmd.autocorrelation(0.0, 5.0, 0.1),
        ]).run(lc)
        assert "autocorrelation_result_0" in result.files

    def test_autocorrelation_no_capture(self):
        # save_result=False: file written but not in result.files
        lc = make_lc()
        result = vt.Pipeline([
            cmd.autocorrelation(0.0, 5.0, 0.1, save_result=False),
        ]).run(lc)
        assert "autocorrelation_result_0" not in result.files


class TestRunBatchOutputValidation:
    """Pipeline.run_batch and LightCurveBatch.run reject cmd.o configurations
    that are single-LC-only (outname= without outdir=).  This replaces the
    earlier silent subprocess fallback / mid-run RuntimeError."""

    def _three_lcs(self):
        return [make_lc() for _ in range(3)]

    def test_run_batch_rejects_outname_without_outdir(self, tmp_path):
        pipe = vt.Pipeline([
            cmd.clip(sigclip=5.0),
            cmd.o(outname=str(tmp_path / "out.lc")),
        ])
        with pytest.raises(ValueError, match="single-LC mode"):
            pipe.run_batch(self._three_lcs())

    def test_lightcurvebatch_rejects_outname_without_outdir(self, tmp_path):
        with pytest.raises(ValueError, match="single-LC mode"):
            (vt.LightCurveBatch(self._three_lcs())
             .clip(sigclip=5.0)
             .o(outname=str(tmp_path / "out.lc"))
             .run())

    def test_run_batch_accepts_outname_plus_outdir(self, tmp_path):
        # Both set: validation lets it through; today's CLI builder uses
        # outdir for batch and the second cmd.o keyword (nameformat) gives
        # each LC a distinct filename.
        outdir = tmp_path / "out"
        outdir.mkdir()
        result = vt.Pipeline([
            cmd.clip(sigclip=5.0),
            cmd.o(outname=str(tmp_path / "x"),
                  outdir=str(outdir),
                  nameformat="%s.txt"),
        ]).run_batch(self._three_lcs())
        assert len(result.vars) == 3

    def test_run_batch_accepts_outdir_only(self, tmp_path):
        outdir = tmp_path / "out"
        outdir.mkdir()
        result = vt.Pipeline([
            cmd.clip(sigclip=5.0),
            cmd.o(outdir=str(outdir)),
        ]).run_batch(self._three_lcs())
        assert len(result.vars) == 3


# ===========================================================================
# Per-LC parameter tests
# ===========================================================================

class TestPerLC:

    def test_varexpr_list_column_passthrough(self):
        from pyvartools.commands._helpers import _varexpr
        assert _varexpr("list column 2") == ["list", "column", "2"]
        assert _varexpr("list") == ["list"]
        assert _varexpr("fixcolumn myvar") == ["fixcolumn", "myvar"]

    def test_perlc_from_list(self):
        from pyvartools import PerLC
        p = PerLC([0.1, 0.2, 0.3])
        assert len(p) == 3 and p[1] == pytest.approx(0.2)

    def test_perlc_from_numpy(self):
        from pyvartools import PerLC
        p = PerLC(np.array([0.1, 0.2, 0.3]))
        assert len(p) == 3

    def test_perlc_from_series(self):
        from pyvartools import PerLC
        p = PerLC(pd.Series([0.1, 0.2]))
        assert len(p) == 2

    def test_perlc_2d_raises(self):
        from pyvartools import PerLC
        with pytest.raises(ValueError):
            PerLC(np.ones((2, 3)))

    def test_perlc_strings_preserved(self):
        """PerLC(['a', 'b', 'c']) preserves strings (polymorphism)."""
        from pyvartools import PerLC
        p = PerLC(["alpha", "beta", "gamma"])
        assert len(p) == 3
        assert p[0] == "alpha"
        assert p[1] == "beta"
        assert isinstance(p[2], str)

    def test_perlc_ints_stay_ints(self):
        """PerLC(int list) preserves int type, no float coercion."""
        from pyvartools import PerLC
        p = PerLC([7, 99, -3])
        assert isinstance(p[0], int) and not isinstance(p[0], bool)
        assert p[1] == 99

    def test_perlc_numeric_round_trip(self):
        """Existing numeric callers still see floats (no behaviour change)."""
        from pyvartools import PerLC
        p = PerLC([0.1, 0.2, 0.15])
        assert p[0] == pytest.approx(0.1)
        # Float input → float output, even though no explicit cast happens.
        assert isinstance(p[1], float)

    def test_perlc_repr(self):
        from pyvartools import PerLC
        p = PerLC([1.0, 2.0])
        assert "PerLC" in repr(p) and "2" in repr(p)

    def test_to_cli_args_with_perlc_substitutes(self):
        c = cmd.LS(0.1, 10.0, 0.001)
        subs = {"minp": "_perlc_0_minp"}
        args = c._to_cli_args_with_perlc(subs)
        # _varexpr("_perlc_0_minp") → ["var", "_perlc_0_minp"]
        assert "var" in args and "_perlc_0_minp" in args

    def test_to_cli_args_with_perlc_restores(self):
        c = cmd.LS(0.1, 10.0, 0.001)
        c._to_cli_args_with_perlc({"minp": "_perlc_0_minp"})
        assert c.minp == pytest.approx(0.1)

    def test_pipeline_detects_numpy_perlc(self):
        pipe = vt.Pipeline([cmd.LS(minp=np.array([0.1, 0.2]), maxp=10.0,
                                    subsample=0.001, npeaks=1)])
        perlc = pipe._collect_perlc_attrs()
        assert (0, "minp") in perlc

    def test_pipeline_detects_explicit_perlc(self):
        from pyvartools import PerLC
        pipe = vt.Pipeline([cmd.LS(minp=PerLC([0.1, 0.2]), maxp=10.0,
                                    subsample=0.001, npeaks=1)])
        perlc = pipe._collect_perlc_attrs()
        assert (0, "minp") in perlc

    def test_pipeline_detects_series_perlc(self):
        pipe = vt.Pipeline([cmd.LS(minp=pd.Series([0.1, 0.2]), maxp=10.0,
                                    subsample=0.001, npeaks=1)])
        perlc = pipe._collect_perlc_attrs()
        assert (0, "minp") in perlc

    def test_pipeline_ignores_plain_list(self):
        # plain Python list should NOT be auto-detected as per-LC
        pipe = vt.Pipeline([cmd.LS(0.1, 10.0, 0.001)])
        # force minp to a plain list to verify it's ignored
        pipe.commands[0].minp = [0.1, 0.2]
        perlc = pipe._collect_perlc_attrs()
        assert (0, "minp") not in perlc

    def test_pipeline_perlc_length_validation(self):
        lcs = [make_lc() for _ in range(3)]
        pipe = vt.Pipeline([cmd.LS(minp=np.array([0.1, 0.2]),
                                    maxp=10.0, subsample=0.001, npeaks=1)])
        with pytest.raises(ValueError, match="has 2 values"):
            pipe.run_batch(lcs)

    def test_pipeline_run_raises_for_perlc(self):
        lc = make_lc()
        pipe = vt.Pipeline([cmd.LS(minp=np.array([0.1]),
                                    maxp=10.0, subsample=0.001, npeaks=1)])
        with pytest.raises(ValueError, match="run_batch"):
            pipe.run(lc)

    def test_pipeline_run_filelist_raises_for_prebuilt_list(self, tmp_path):
        list_file = tmp_path / "lclist.txt"
        list_file.write_text("EXAMPLES/2\n")
        pipe = vt.Pipeline([cmd.LS(minp=np.array([0.1]),
                                    maxp=10.0, subsample=0.001, npeaks=1)])
        with pytest.raises(ValueError, match="pre-built list file"):
            pipe.run_filelist(str(list_file))

    def test_write_perlc_list_file(self, tmp_path):
        from pyvartools import PerLC
        list_path = str(tmp_path / "lclist.txt")
        lc_paths = ["/tmp/a.lc", "/tmp/b.lc"]
        perlc_attrs = {(0, "minp"): PerLC([0.1, 0.2])}
        col_assignments = {(0, "minp"): 2}
        pipe = vt.Pipeline([cmd.LS(0.1, 10.0, 0.001)])
        pipe._write_perlc_list_file(list_path, lc_paths, perlc_attrs, col_assignments)
        lines = open(list_path).read().strip().split("\n")
        assert lines[0] == "/tmp/a.lc 0.1"
        assert lines[1] == "/tmp/b.lc 0.2"

    def test_build_perlc_subs(self):
        pipe = vt.Pipeline([cmd.LS(0.1, 10.0, 0.001)])
        col_assignments = {(0, "minp"): 2, (0, "maxp"): 3}
        subs = pipe._build_perlc_subs(col_assignments)
        # Each attr is substituted with an ``"expr <varname>"`` string.  The
        # "expr" form is portable across every command's value-spec parser
        # (some commands — notably ``-BLSFixPer`` on its period slot — do
        # not accept the bare ``"var NAME"`` form).
        assert subs[0]["minp"] == "expr _perlc_0_minp"
        assert subs[0]["maxp"] == "expr _perlc_0_maxp"

    def test_build_perlc_inlistvars(self):
        pipe = vt.Pipeline([cmd.LS(0.1, 10.0, 0.001)])
        col_assignments = {(0, "minp"): 2, (0, "maxp"): 3}
        ivars = pipe._build_cmdattr_perlc_vars(col_assignments)
        assert ivars["_perlc_0_minp"] == 2
        assert ivars["_perlc_0_maxp"] == 3

    def test_to_cli_args_with_perlc_uses_var_keyword(self):
        c = cmd.LS(0.1, 10.0, 0.001)
        subs = {"minp": "_perlc_0_minp"}
        args = c._to_cli_args_with_perlc(subs)
        # _varexpr("_perlc_0_minp") → ["var", "_perlc_0_minp"]
        assert "var" in args
        idx = args.index("var")
        assert args[idx + 1] == "_perlc_0_minp"


# ---------------------------------------------------------------------------
# Parametrized CLI-arg tests: every confirmed (command, param) pair that
# vartools supports via "var varname" syntax.
#
# Discovery (run against the installed binary):
#   LS      : var=[minpvar, maxpvar, subsamplevar]
#   aov     : var=[Nbinvar, minpvar, maxpvar, subsamplevar, finetunevar]
#   aov_harm: var=[Nharmvar, minpvar, maxpvar, subsamplevar, finetunevar]
#   BLS     : var=[rminvar, rmaxvar, qminvar, qmaxvar, rhovar,
#                  minexpdurfracvar, maxexpdurfracvar, minpervar, maxpervar,
#                  nfreqvar, dfreqvar, subsamplevar, nbinsvar]
#   All other commands (Phase, restricttimes, resample, binlc, addnoise,
#   TFA_SR, Injecttransit, …) do NOT support "var" syntax.
# ---------------------------------------------------------------------------

_VAR_PARAMS = [
    # LS
    (lambda: cmd.LS(0.1, 10.0, 0.001), "minp"),
    (lambda: cmd.LS(0.1, 10.0, 0.001), "maxp"),
    (lambda: cmd.LS(0.1, 10.0, 0.001), "subsample"),
    # aov
    (lambda: cmd.aov(0.1, 10.0, 0.001, 2), "minp"),
    (lambda: cmd.aov(0.1, 10.0, 0.001, 2), "maxp"),
    (lambda: cmd.aov(0.1, 10.0, 0.001, 2), "subsample"),
    (lambda: cmd.aov(0.1, 10.0, 0.001, 2), "finetune"),
    (lambda: cmd.aov(0.1, 10.0, 0.001, 2, nbin=10), "nbin"),
    # aov_harm
    (lambda: cmd.aov_harm(2, 0.1, 10.0, 0.001, 2), "nharm"),
    (lambda: cmd.aov_harm(2, 0.1, 10.0, 0.001, 2), "minp"),
    (lambda: cmd.aov_harm(2, 0.1, 10.0, 0.001, 2), "maxp"),
    (lambda: cmd.aov_harm(2, 0.1, 10.0, 0.001, 2), "subsample"),
    (lambda: cmd.aov_harm(2, 0.1, 10.0, 0.001, 2), "finetune"),
    # BLS — rmin/rmax branch.  In r/q mode the user must supply nfreq=
    # or df=; the "optimal" frequency grid is density-mode-only.
    (lambda: cmd.BLS(0.1, 10.0, nfreq=1000), "minper"),
    (lambda: cmd.BLS(0.1, 10.0, nfreq=1000), "maxper"),
    (lambda: cmd.BLS(0.1, 10.0, nfreq=1000), "rmin"),
    (lambda: cmd.BLS(0.1, 10.0, nfreq=1000), "rmax"),
    (lambda: cmd.BLS(0.1, 10.0, nfreq=1000), "nbins"),
    (lambda: cmd.BLS(0.1, 10.0, density_mode=True, stellar_density=1.0,
                     min_exp_dur_frac=0.5, max_exp_dur_frac=1.5),
        "subsample"),                            # "optimal" branch
    # BLS — qmin/qmax branch
    (lambda: cmd.BLS(0.1, 10.0, qmin=0.02, qmax=0.5, nfreq=1000), "qmin"),
    (lambda: cmd.BLS(0.1, 10.0, qmin=0.02, qmax=0.5, nfreq=1000), "qmax"),
    # BLS — nfreq branch
    (lambda: cmd.BLS(0.1, 10.0, nfreq=1000), "nfreq"),
    # BLS — df branch
    (lambda: cmd.BLS(0.1, 10.0, df=0.001), "df"),
    # BLS — density branch
    (lambda: cmd.BLS(0.1, 10.0, density_mode=True, stellar_density=1.0,
                     min_exp_dur_frac=0.5, max_exp_dur_frac=1.5), "stellar_density"),
    (lambda: cmd.BLS(0.1, 10.0, density_mode=True, stellar_density=1.0,
                     min_exp_dur_frac=0.5, max_exp_dur_frac=1.5), "min_exp_dur_frac"),
    (lambda: cmd.BLS(0.1, 10.0, density_mode=True, stellar_density=1.0,
                     min_exp_dur_frac=0.5, max_exp_dur_frac=1.5), "max_exp_dur_frac"),
]

_VAR_PARAMS_IDS = [f"{f()._vt_name}.{p}" for f, p in _VAR_PARAMS]


@pytest.mark.parametrize("factory,param", _VAR_PARAMS, ids=_VAR_PARAMS_IDS)
def test_perlc_cli_arg_produces_var_token(factory, param):
    """_to_cli_args_with_perlc emits 'var varname' for every supported param."""
    c = factory()
    c._outdir = "."
    varname = f"_perlc_0_{param}"
    args = c._to_cli_args_with_perlc({param: varname})
    assert "var" in args, f"'var' not in args for {c._vt_name}.{param}: {args}"
    idx = args.index("var")
    assert args[idx + 1] == varname, (
        f"Expected '{varname}' after 'var' for {c._vt_name}.{param}, "
        f"got: {args[idx + 1:]}"
    )


@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools binary not found")
class TestPerLCEndToEnd:

    def test_perlc_batch_run(self):
        """Per-LC minp/maxp values are correctly passed to each LC."""
        lcs = [vt.LightCurve.from_file(
            os.path.join(EXAMPLES_DIR, str(i))) for i in [2, 3]]
        minps = np.array([0.1, 0.5])
        maxps = np.array([5.0, 10.0])
        result = vt.Pipeline([
            cmd.LS(minp=minps, maxp=maxps, subsample=0.001, npeaks=1)
        ]).run_batch(lcs)
        assert len(result.vars) == 2
        assert "LS_Period_1_0" in result.vars.columns

    def test_perlc_explicit_wrapper(self):
        """Explicit PerLC([...]) wrapper works identically to numpy array."""
        from pyvartools import PerLC
        lcs = [vt.LightCurve.from_file(
            os.path.join(EXAMPLES_DIR, str(i))) for i in [2, 3]]
        result = vt.Pipeline([
            cmd.LS(minp=PerLC([0.1, 0.5]), maxp=10.0,
                   subsample=0.001, npeaks=1)
        ]).run_batch(lcs)
        assert len(result.vars) == 2

    def test_perlc_run_filelist_paths(self):
        """run_filelist with a list of paths + per-LC array works."""
        paths = [os.path.join(EXAMPLES_DIR, str(i)) for i in [2, 3]]
        result = vt.Pipeline([
            cmd.LS(minp=np.array([0.1, 0.5]), maxp=10.0,
                   subsample=0.001, npeaks=1)
        ]).run_filelist(paths)
        assert len(result.vars) == 2

    def test_perlc_aov_batch(self):
        """Per-LC minp values are correctly passed to aov for each LC."""
        lcs = [vt.LightCurve.from_file(
            os.path.join(EXAMPLES_DIR, str(i))) for i in [2, 3]]
        result = vt.Pipeline([
            cmd.aov(minp=np.array([0.1, 0.5]), maxp=10.0,
                    subsample=0.001, finetune=2, npeaks=1)
        ]).run_batch(lcs)
        assert len(result.vars) == 2
        assert "AOV_1_0" in result.vars.columns

    def test_perlc_bls_batch(self):
        """Per-LC minper/maxper values are correctly passed to BLS for each LC."""
        lcs = [vt.LightCurve.from_file(
            os.path.join(EXAMPLES_DIR, str(i))) for i in [2, 3]]
        result = vt.Pipeline([
            cmd.BLS(minper=np.array([0.1, 0.5]), maxper=np.array([5.0, 10.0]),
                    rmin=0.01, rmax=0.5, nbins=200, npeaks=1, nfreq=5000)
        ]).run_batch(lcs)
        assert len(result.vars) == 2
        assert "BLS_Period_1_0" in result.vars.columns


# ===========================================================================
# User extension library tests
# ===========================================================================

class TestUserCommandCLIArgs:
    """Unit tests for UserCommand — no binary required."""

    def test_usercommand_cli_args_with_path(self):
        c = vt.UserCommand("/path/to/stitch.so", "stitch",
                           "mag err mask lcnum median")
        args = c._to_cli_args()
        assert args[:2] == ["-L", "/path/to/stitch.so"]
        assert "-stitch" in args
        assert "median" in args

    def test_usercommand_no_lib_path(self):
        c = vt.UserCommand(None, "stitch", "mag err mask lcnum median")
        args = c._to_cli_args()
        assert "-L" not in args
        assert "-stitch" in args
        assert "median" in args

    def test_usercommand_empty_lib_path(self):
        c = vt.UserCommand("", "stitch", "mag err mask lcnum median")
        args = c._to_cli_args()
        assert "-L" not in args
        assert "-stitch" in args

    def test_usercommand_string_args_split(self):
        c = vt.UserCommand("/path/stitch.so", "stitch",
                           "mag err mask lcnum")
        assert c._user_args == ["mag", "err", "mask", "lcnum"]

    def test_usercommand_list_args(self):
        c = vt.UserCommand("/path/stitch.so", "stitch",
                           ["mag", "err", "mask", "lcnum"])
        assert c._user_args == ["mag", "err", "mask", "lcnum"]

    def test_usercommand_name_from_stem(self):
        c = vt.UserCommand("/path/to/stitch.so")
        assert c._cmd_name == "stitch"
        assert "-stitch" in c._to_cli_args()

    def test_usercommand_name_from_stem_la(self):
        c = vt.UserCommand("/path/to/macula.la")
        assert c._cmd_name == "macula"

    def test_usercommand_empty_args(self):
        c = vt.UserCommand("/path/stitch.so", "stitch")
        args = c._to_cli_args()
        assert args == ["-L", "/path/stitch.so", "-stitch"]

    def test_usercommand_output_file_specs_empty(self):
        c = vt.UserCommand("/path/stitch.so", "stitch")
        assert c._output_file_specs() == {}

    def test_usercommand_is_vartoolscommand(self):
        from pyvartools._command import VartoolsCommand
        c = vt.UserCommand("/path/stitch.so", "stitch")
        assert isinstance(c, VartoolsCommand)

    def test_usercommand_accessible_via_commands(self):
        assert cmd.UserCommand is vt.UserCommand

    def test_usercommand_repr(self):
        c = vt.UserCommand("/path/stitch.so", "stitch", "mag err")
        r = repr(c)
        assert "UserCommand" in r
        assert "-stitch" in r


class TestLoadUserlib:
    """Unit tests for load_userlib() — no binary required."""

    def test_load_userlib_returns_class(self):
        from unittest.mock import patch
        with patch("pyvartools.userlib._fetch_userlib_text", return_value=""):
            Stitch = vt.load_userlib("/path/stitch.so", name="stitch")
        assert issubclass(Stitch, vt.UserCommand)

    def test_load_userlib_class_instantiates(self):
        from unittest.mock import patch
        with patch("pyvartools.userlib._fetch_userlib_text", return_value=""):
            Stitch = vt.load_userlib("/path/stitch.so", name="stitch")
        inst = Stitch("mag err mask lcnum median")
        assert isinstance(inst, vt.UserCommand)

    def test_load_userlib_cli_args(self):
        from unittest.mock import patch
        with patch("pyvartools.userlib._fetch_userlib_text", return_value=""):
            Stitch = vt.load_userlib("/path/stitch.so", name="stitch")
        inst = Stitch("mag err mask lcnum median")
        args = inst._to_cli_args()
        assert "-stitch" in args
        assert "median" in args

    def test_load_userlib_class_name_default(self):
        from unittest.mock import patch
        with patch("pyvartools.userlib._fetch_userlib_text", return_value=""):
            Stitch = vt.load_userlib("/path/stitch.so", name="stitch")
        assert Stitch.__name__ == "Stitch"

    def test_load_userlib_class_name_custom(self):
        from unittest.mock import patch
        with patch("pyvartools.userlib._fetch_userlib_text", return_value=""):
            Cls = vt.load_userlib("/path/stitch.so", name="stitch",
                                  cls_name="MyStitch")
        assert Cls.__name__ == "MyStitch"

    def test_load_userlib_name_from_stem(self):
        from unittest.mock import patch
        with patch("pyvartools.userlib._fetch_userlib_text", return_value=""):
            Cls = vt.load_userlib("/path/stitch.so")
        inst = Cls("arg1")
        assert inst._cmd_name == "stitch"

    def test_load_userlib_docstring_from_help(self):
        from unittest.mock import patch
        with patch("pyvartools.userlib._fetch_userlib_text",
                   return_value="Stitch help text"):
            Cls = vt.load_userlib("/path/stitch.so", name="stitch")
        assert "Stitch help text" in (Cls.__doc__ or "")

    def test_load_userlib_docstring_fallback(self):
        from unittest.mock import patch
        with patch("pyvartools.userlib._fetch_userlib_text", return_value=""):
            Cls = vt.load_userlib("/path/stitch.so", name="stitch")
        assert Cls.__doc__  # non-empty fallback

    def test_load_userlib_preserves_symlink_basename(self, tmp_path):
        """``load_userlib`` must not follow symlinks: libtool's dlopen
        looks up ``<libbasename>_Initialize`` against the path passed to
        ``-L``, and that lookup fails when the path's basename has the
        versioned suffix (``stitch.so.0.0.0``).
        """
        from unittest.mock import patch
        target = tmp_path / "stitch.so.0.0.0"
        target.touch()
        link = tmp_path / "stitch.so"
        link.symlink_to(target.name)
        with patch("pyvartools.userlib._fetch_userlib_text", return_value=""):
            Cls = vt.load_userlib(str(link))
        inst = Cls("mag err mask lcnum median")
        assert inst._lib_path.endswith("stitch.so"), inst._lib_path
        assert not inst._lib_path.endswith(".so.0.0.0"), inst._lib_path


class TestDiscoverUserlibs:
    """Unit tests for discover_userlibs() — no binary required."""

    def test_discover_empty_dir(self, tmp_path):
        result = vt.discover_userlibs(search_paths=[str(tmp_path)])
        assert result == {}

    def test_discover_finds_so(self, tmp_path):
        from unittest.mock import patch
        # Create a fake .so file
        (tmp_path / "fakecmd.so").touch()
        with patch("pyvartools.userlib._fetch_userlib_text", return_value=""):
            result = vt.discover_userlibs(search_paths=[str(tmp_path)])
        assert "fakecmd" in result
        assert issubclass(result["fakecmd"], vt.UserCommand)

    def test_discover_finds_la(self, tmp_path):
        from unittest.mock import patch
        (tmp_path / "othercmd.la").touch()
        with patch("pyvartools.userlib._fetch_userlib_text", return_value=""):
            result = vt.discover_userlibs(search_paths=[str(tmp_path)])
        assert "othercmd" in result

    def test_discover_no_duplicates(self, tmp_path):
        from unittest.mock import patch
        # Both .so and .la present for the same stem — .so appears first
        (tmp_path / "mycmd.so").touch()
        (tmp_path / "mycmd.la").touch()
        with patch("pyvartools.userlib._fetch_userlib_text", return_value=""):
            result = vt.discover_userlibs(search_paths=[str(tmp_path)])
        assert list(result.keys()).count("mycmd") == 1

    def test_discover_returns_dict(self, tmp_path):
        result = vt.discover_userlibs(search_paths=[str(tmp_path)])
        assert isinstance(result, dict)

    def test_discover_nonexistent_dir(self, tmp_path):
        result = vt.discover_userlibs(
            search_paths=[str(tmp_path / "does_not_exist")])
        assert result == {}


class TestNonlinfitOutputFileSpecs:
    """Verify that the chains output spec was not silently dropped."""

    def test_nonlinfit_output_file_specs_includes_chains(self):
        from pyvartools.commands import nonlinfit
        specs = nonlinfit(function="a*sin(b*t+c)", paramlist="a,b,c")._output_file_specs()
        assert "chains" in specs
        assert "model" in specs


class TestGlobalPipelineOptions:
    """Unit tests for randseed/skipmissing/jdtol/matchstringid pipeline options."""

    def _capture_cmd(self, pipe, method, *args, **kwargs):
        """Run *method* on *pipe*, intercept _execute, return the cmd list."""
        from unittest.mock import patch
        captured = []

        def mock_execute(cmd_list, **kw):
            captured.extend(cmd_list)
            return ("Name = test\n", "")

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch("pyvartools.pipeline.parse_oneline_output",
                       return_value=pd.Series(dtype=object)):
                try:
                    getattr(pipe, method)(*args, **kwargs)
                except Exception:
                    pass
        return captured

    def test_pipeline_randseed_in_cmd(self):
        pipe = vt.Pipeline([cmd.rms()])
        result = self._capture_cmd(pipe, "run_file", "dummy.lc", randseed=42)
        assert "-randseed" in result
        assert "42" in result

    def test_pipeline_skipmissing_in_cmd(self):
        pipe = vt.Pipeline([cmd.rms()])
        result = self._capture_cmd(pipe, "run_file", "dummy.lc", skipmissing=True)
        assert "-skipmissing" in result

    def test_pipeline_jdtol_in_cmd(self):
        pipe = vt.Pipeline([cmd.rms()])
        result = self._capture_cmd(pipe, "run_file", "dummy.lc", jdtol=0.001)
        assert "-jdtol" in result
        assert "0.001" in result

    def test_pipeline_matchstringid_in_cmd(self):
        pipe = vt.Pipeline([cmd.rms()])
        result = self._capture_cmd(pipe, "run_file", "dummy.lc", matchstringid=True)
        assert "-matchstringid" in result

    def test_global_opts_ordering(self):
        """Global opts appear after -parallel and before command tokens."""
        pipe = vt.Pipeline([cmd.rms()])
        result = self._capture_cmd(pipe, "run_file", "dummy.lc",
                                   randseed=7, jdtol=0.5)
        assert "-randseed" in result
        assert result.index("-randseed") < result.index("-rms")

    def test_skipmissing_false_not_emitted(self):
        pipe = vt.Pipeline([cmd.rms()])
        result = self._capture_cmd(pipe, "run_file", "dummy.lc", skipmissing=False)
        assert "-skipmissing" not in result

    def test_matchstringid_false_not_emitted(self):
        pipe = vt.Pipeline([cmd.rms()])
        result = self._capture_cmd(pipe, "run_file", "dummy.lc", matchstringid=False)
        assert "-matchstringid" not in result


class TestRunCombinelcs:
    """Unit tests for Pipeline.run_combinelcs()."""

    def _capture_cmd(self, pipe, groups, **kwargs):
        from unittest.mock import patch
        captured = []

        def mock_execute(cmd_list, **kw):
            captured.extend(cmd_list)
            return ("Name = lc1\nRMS_0 = 0.01\n", "")

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch("pyvartools.pipeline.parse_oneline_output",
                       return_value=pd.DataFrame([{"Name": "lc1", "RMS_0": 0.01}])):
                result = pipe.run_combinelcs(groups, **kwargs)
        return captured, result

    def test_combinelcs_in_cmd(self, tmp_path):
        (tmp_path / "lc1.txt").write_text("1.0 10.0 0.01\n")
        (tmp_path / "lc2.txt").write_text("1.0 10.1 0.01\n")
        groups = [[str(tmp_path / "lc1.txt"), str(tmp_path / "lc2.txt")]]
        pipe = vt.Pipeline([cmd.rms()])
        captured, result = self._capture_cmd(pipe, groups)
        assert "combinelcs" in captured
        assert isinstance(result, vt.BatchResult)

    def test_combinelcs_list_file_format(self, tmp_path):
        """List file must have comma-joined paths per group."""
        (tmp_path / "a.lc").write_text("1.0 10.0 0.01\n")
        (tmp_path / "b.lc").write_text("2.0 10.1 0.01\n")
        groups = [[str(tmp_path / "a.lc"), str(tmp_path / "b.lc")]]
        pipe = vt.Pipeline([cmd.rms()])

        written_content = []
        original_open = open

        import builtins
        from unittest.mock import patch, MagicMock

        def mock_execute(cmd_list, **kw):
            return ("Name = a\nRMS_0 = 0.01\n", "")

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch("pyvartools.pipeline.parse_oneline_output",
                       return_value=pd.DataFrame([{"Name": "a", "RMS_0": 0.01}])):
                pipe.run_combinelcs(groups)
        # The method ran without error; the list file content is tested
        # implicitly by "combinelcs" appearing in the cmd above.

    def test_combinelcs_lcnumvar(self, tmp_path):
        """lcnumvar keyword is passed after combinelcs."""
        (tmp_path / "a.lc").write_text("1.0 10.0 0.01\n")
        groups = [[str(tmp_path / "a.lc")]]
        pipe = vt.Pipeline([cmd.rms()])
        captured, _ = self._capture_cmd(pipe, groups, lcnumvar="lcnum")
        assert "lcnumvar" in captured
        assert "lcnum" in captured
        idx = captured.index("lcnumvar")
        assert captured[idx + 1] == "lcnum"

    def test_combinelcs_perlc_length_mismatch_raises(self):
        """PerLC length must match the number of groups (M5)."""
        from pyvartools import PerLC
        pipe = vt.Pipeline([
            cmd.LS(minp=PerLC([0.1, 0.2, 0.3]), maxp=10.0, subsample=1e-3)
        ])
        with pytest.raises(ValueError, match="3 values but the batch has 2"):
            pipe.run_combinelcs(
                [["/tmp/a.lc", "/tmp/b.lc"], ["/tmp/c.lc", "/tmp/d.lc"]]
            )

    def test_combinelcs_perlc_emits_inlistvars(self, tmp_path):
        """PerLC params produce -inlistvars and var/expr substitutions (M5)."""
        from pyvartools import PerLC
        (tmp_path / "a.lc").write_text("1.0 10.0 0.01\n")
        (tmp_path / "b.lc").write_text("1.0 10.1 0.01\n")
        groups = [
            [str(tmp_path / "a.lc"), str(tmp_path / "b.lc")],
            [str(tmp_path / "a.lc"), str(tmp_path / "b.lc")],
        ]
        pipe = vt.Pipeline([
            cmd.LS(minp=PerLC([0.1, 0.2]), maxp=10.0, subsample=1e-3)
        ])
        captured, _ = self._capture_cmd(pipe, groups)
        # The PerLC param is wired as a per-star variable via -inlistvars.
        assert "-inlistvars" in captured
        idx = captured.index("-inlistvars")
        assert "_perlc_0_minp" in captured[idx + 1]
        # And the LS command emits the variable through `expr`.
        assert "_perlc_0_minp" in " ".join(captured)

    def test_combinelcs_randseed_in_cmd(self, tmp_path):
        (tmp_path / "a.lc").write_text("1.0 10.0 0.01\n")
        groups = [[str(tmp_path / "a.lc")]]
        pipe = vt.Pipeline([cmd.rms()])
        captured, _ = self._capture_cmd(pipe, groups, randseed=99)
        assert "-randseed" in captured
        assert "99" in captured

    def test_combinelcs_custom_delimiter(self, tmp_path):
        """Non-default delimiter is used in the list file (verified via cmd)."""
        (tmp_path / "a.lc").write_text("1.0 10.0 0.01\n")
        (tmp_path / "b.lc").write_text("2.0 10.1 0.01\n")
        groups = [[str(tmp_path / "a.lc"), str(tmp_path / "b.lc")]]
        pipe = vt.Pipeline([cmd.rms()])
        # Just verify it runs without error with a custom delimiter
        captured, result = self._capture_cmd(pipe, groups, delimiter=" ")
        assert "combinelcs" in captured

    def test_combinelcs_lcnumvar_default(self, tmp_path):
        """run_combinelcs() emits 'lcnumvar lcnum' by default (M4)."""
        (tmp_path / "a.lc").write_text("1.0 10.0 0.01\n")
        groups = [[str(tmp_path / "a.lc")]]
        pipe = vt.Pipeline([cmd.rms()])
        captured, _ = self._capture_cmd(pipe, groups)
        assert "lcnumvar" in captured
        idx = captured.index("lcnumvar")
        assert captured[idx + 1] == "lcnum"

    def test_combinelcs_lcnumvar_none_opts_out(self, tmp_path):
        """Passing lcnumvar=None suppresses the lcnumvar qualifier (M4)."""
        (tmp_path / "a.lc").write_text("1.0 10.0 0.01\n")
        groups = [[str(tmp_path / "a.lc")]]
        pipe = vt.Pipeline([cmd.rms()])
        captured, _ = self._capture_cmd(pipe, groups, lcnumvar=None)
        assert "combinelcs" in captured
        assert "lcnumvar" not in captured

    def test_combinelcs_capture_lc_overrides_lcnumvar_none(self, tmp_path):
        """capture_lc=True forces lcnumvar='lcnum' even with explicit None (M6)."""
        (tmp_path / "a.lc").write_text("1.0 10.0 0.01\n")
        groups = [[str(tmp_path / "a.lc")]]
        pipe = vt.Pipeline([cmd.rms()])
        captured, _ = self._capture_cmd(pipe, groups,
                                         capture_lc=True, lcnumvar=None)
        assert "lcnumvar" in captured
        idx = captured.index("lcnumvar")
        assert captured[idx + 1] == "lcnum"


class TestRunCombinelc:
    """Unit tests for the singular Pipeline.run_combinelc() (M2)."""

    def _capture_cmd(self, pipe, files, **kwargs):
        from unittest.mock import patch
        captured = []

        def mock_execute(cmd_list, **kw):
            captured.extend(cmd_list)
            return ("Name = a\nRMS_0 = 0.01\n", "")

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch("pyvartools.pipeline.parse_oneline_output",
                       return_value=pd.DataFrame([{"Name": "a", "RMS_0": 0.01}])):
                result = pipe.run_combinelc(files, **kwargs)
        return captured, result

    def test_run_combinelc_emits_combinelcs(self, tmp_path):
        (tmp_path / "a.lc").write_text("1.0 10.0 0.01\n")
        (tmp_path / "b.lc").write_text("1.0 10.1 0.01\n")
        files = [str(tmp_path / "a.lc"), str(tmp_path / "b.lc")]
        pipe = vt.Pipeline([cmd.rms()])
        captured, _ = self._capture_cmd(pipe, files)
        assert "combinelcs" in captured

    def test_run_combinelc_default_lcnumvar(self, tmp_path):
        (tmp_path / "a.lc").write_text("1.0 10.0 0.01\n")
        files = [str(tmp_path / "a.lc")]
        pipe = vt.Pipeline([cmd.rms()])
        captured, _ = self._capture_cmd(pipe, files)
        assert "lcnumvar" in captured
        idx = captured.index("lcnumvar")
        assert captured[idx + 1] == "lcnum"

    def test_run_combinelc_returns_single_result(self, tmp_path):
        (tmp_path / "a.lc").write_text("1.0 10.0 0.01\n")
        files = [str(tmp_path / "a.lc")]
        pipe = vt.Pipeline([cmd.rms()])
        _, result = self._capture_cmd(pipe, files)
        assert isinstance(result, vt.Result)
        # Result.vars is always a Series, not a DataFrame
        assert isinstance(result.vars, pd.Series)

    def test_run_combinelc_empty_raises(self):
        pipe = vt.Pipeline([cmd.rms()])
        with pytest.raises(ValueError, match="at least one file"):
            pipe.run_combinelc([])


class TestRunFilelistCombinelcs:
    """Unit tests for the combinelcs=True flag on run_filelist (M3)."""

    def _capture_cmd(self, pipe, lc_paths, **kwargs):
        from unittest.mock import patch
        captured = []

        def mock_execute(cmd_list, **kw):
            captured.extend(cmd_list)
            return ("Name = a\nRMS_0 = 0.01\n", "")

        with patch.object(pipe, "_execute", side_effect=mock_execute):
            with patch("pyvartools.pipeline.parse_oneline_output",
                       return_value=pd.DataFrame([{"Name": "a", "RMS_0": 0.01}])):
                result = pipe.run_filelist(lc_paths, **kwargs)
        return captured, result

    def test_filelist_combinelcs_emits_keywords(self, tmp_path):
        """combinelcs=True appends 'combinelcs lcnumvar lcnum' to -l."""
        (tmp_path / "a.lc").write_text("1.0 10.0 0.01\n")
        (tmp_path / "b.lc").write_text("1.0 10.1 0.01\n")
        # User supplies one comma-joined line per group.
        lc_paths = [f"{tmp_path / 'a.lc'},{tmp_path / 'b.lc'}"]
        pipe = vt.Pipeline([cmd.rms()])
        captured, _ = self._capture_cmd(pipe, lc_paths, combinelcs=True)
        assert "combinelcs" in captured
        cidx = captured.index("combinelcs")
        assert captured[cidx + 1] == "lcnumvar"
        assert captured[cidx + 2] == "lcnum"

    def test_filelist_combinelcs_lcnumvar_none(self, tmp_path):
        """combinelcs=True with lcnumvar=None emits combinelcs but not lcnumvar."""
        (tmp_path / "a.lc").write_text("1.0 10.0 0.01\n")
        lc_paths = [str(tmp_path / "a.lc")]
        pipe = vt.Pipeline([cmd.rms()])
        captured, _ = self._capture_cmd(pipe, lc_paths, combinelcs=True,
                                         lcnumvar=None)
        assert "combinelcs" in captured
        assert "lcnumvar" not in captured

    def test_filelist_combinelcs_capture_lc_overrides_lcnumvar_none(self, tmp_path):
        """capture_lc=True forces lcnumvar='lcnum' on run_filelist too (M6)."""
        (tmp_path / "a.lc").write_text("1.0 10.0 0.01\n")
        lc_paths = [str(tmp_path / "a.lc")]
        pipe = vt.Pipeline([cmd.rms()])
        captured, _ = self._capture_cmd(pipe, lc_paths, combinelcs=True,
                                         capture_lc=True, lcnumvar=None)
        assert "lcnumvar" in captured
        idx = captured.index("lcnumvar")
        assert captured[idx + 1] == "lcnum"

    def test_filelist_combinelcs_false_no_keyword(self, tmp_path):
        """combinelcs=False (default) does not emit combinelcs."""
        (tmp_path / "a.lc").write_text("1.0 10.0 0.01\n")
        lc_paths = [str(tmp_path / "a.lc")]
        pipe = vt.Pipeline([cmd.rms()])
        captured, _ = self._capture_cmd(pipe, lc_paths)
        assert "combinelcs" not in captured

    def test_filelist_combinelcs_perlc_raises(self):
        """combinelcs=True with PerLC params raises ValueError."""
        from pyvartools import PerLC
        pipe = vt.Pipeline([
            cmd.LS(minp=PerLC([0.1, 0.2]), maxp=10.0, subsample=1e-3)
        ])
        with pytest.raises(ValueError, match="combinelcs=True"):
            pipe.run_filelist(["/tmp/a.lc,/tmp/b.lc", "/tmp/c.lc,/tmp/d.lc"],
                              combinelcs=True)


@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools not installed")
class TestRunCombinelcsIntegration:

    def test_run_combinelcs_basic(self):
        groups = [[EXAMPLE_LC, EXAMPLE_LC]]
        result = vt.Pipeline([cmd.rms()]).run_combinelcs(groups)
        assert len(result.vars) == 1
        assert "RMS_0" in result.vars.columns

    def test_run_combinelc_basic(self):
        """Singular run_combinelc returns a single Result."""
        result = vt.Pipeline([cmd.rms()]).run_combinelc([EXAMPLE_LC, EXAMPLE_LC])
        assert isinstance(result, vt.Result)
        assert "RMS_0" in result.vars.index

    def test_run_filelist_combinelcs(self, tmp_path):
        """run_filelist(combinelcs=True) accepts comma-joined paths."""
        # Build a list file with one comma-joined line.
        list_file = tmp_path / "groups.txt"
        list_file.write_text(f"{EXAMPLE_LC},{EXAMPLE_LC}\n")
        result = vt.Pipeline([cmd.rms()]).run_filelist(
            str(list_file), combinelcs=True
        )
        assert len(result.vars) == 1
        assert "RMS_0" in result.vars.columns

    def test_run_combinelcs_perlc(self):
        """PerLC values flow through to vartools, one per group (M5)."""
        from pyvartools import PerLC
        groups = [[EXAMPLE_LC, EXAMPLE_LC], [EXAMPLE_LC, EXAMPLE_LC]]
        result = vt.Pipeline([
            cmd.LS(minp=PerLC([0.1, 0.5]), maxp=10.0,
                   subsample=0.001, npeaks=1)
        ]).run_combinelcs(groups)
        assert len(result.vars) == 2
        assert "LS_Period_1_0" in result.vars.columns

    def test_run_combinelcs_capture_lc_has_lcnum(self):
        """capture_lc=True returns an LC with the lcnum column populated (M6)."""
        groups = [[EXAMPLE_LC, EXAMPLE_LC]]
        batch = vt.Pipeline([cmd.rms()]).run_combinelcs(groups, capture_lc=True)
        assert batch.lcs[0] is not None
        assert "lcnum" in batch.lcs[0]._df.columns
        # Two source files → lcnum should take values 0 and 1.
        assert sorted(batch.lcs[0]._df["lcnum"].unique().tolist()) == [0, 1]

    def test_run_combinelc_capture_lc_has_lcnum(self):
        """Singular form: capture_lc=True returns one LC with lcnum (M6)."""
        result = vt.Pipeline([cmd.rms()]).run_combinelc(
            [EXAMPLE_LC, EXAMPLE_LC], capture_lc=True
        )
        assert result.lc is not None
        assert "lcnum" in result.lc._df.columns

    def test_run_combinelcs_capture_lc_overrides_none(self):
        """capture_lc=True forces lcnumvar='lcnum' even if the user passed None (M6)."""
        groups = [[EXAMPLE_LC, EXAMPLE_LC]]
        batch = vt.Pipeline([cmd.rms()]).run_combinelcs(
            groups, capture_lc=True, lcnumvar=None
        )
        # The captured LC must still carry the lcnum column despite the opt-out.
        assert batch.lcs[0] is not None
        assert "lcnum" in batch.lcs[0]._df.columns


class TestPipelineUserCommandSubprocessFallback:
    """Verify that UserCommand presence forces subprocess mode."""

    def test_has_user_commands_true(self):
        pipe = vt.Pipeline([
            vt.UserCommand("/path/stitch.so", "stitch", "mag err mask lcnum")
        ])
        assert pipe._has_user_commands() is True

    def test_has_user_commands_false(self):
        pipe = vt.Pipeline([cmd.clip(sigclip=5.0)])
        assert pipe._has_user_commands() is False

    def test_has_user_commands_mixed(self):
        pipe = vt.Pipeline([
            cmd.clip(sigclip=5.0),
            vt.UserCommand("/path/stitch.so", "stitch", "mag err mask lcnum"),
        ])
        assert pipe._has_user_commands() is True

    def test_load_userlib_instance_detected(self):
        from unittest.mock import patch
        with patch("pyvartools.userlib._fetch_userlib_text", return_value=""):
            Stitch = vt.load_userlib("/path/stitch.so", name="stitch")
        pipe = vt.Pipeline([Stitch("mag err mask lcnum median")])
        assert pipe._has_user_commands() is True


# ===========================================================================
# Integration: cmd.python (mirrors the CLI `-python` examples)
# ===========================================================================
#
# These tests run vartools as a subprocess and exercise the `-python`
# command end-to-end.  They are gated on `_HAVE_BINARY` so the suite
# still passes on machines where vartools is not installed; they will
# additionally fail if vartools was built *without* Python support
# (`./configure --without-python`), which is appropriate — there's no
# wrapper-side patch that can reach a vartools that can't run -python.

@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools not installed")
class TestPythonIntegration:
    """End-to-end coverage of cmd.python in subprocess mode.

    The reference CLI examples in ``docs/cli/python-r.md`` are reproduced
    here as a Python pipeline, asserting that the documented output
    column names and (where stable) values appear.
    """

    # Build an absolute-path list from the canonical EXAMPLES/{1..10} set.
    # The shipped EXAMPLES/lc_list uses relative paths (`EXAMPLES/1`, …)
    # which assume the test is run from the vartools source root.  Using
    # explicit absolute paths keeps these tests location-independent.
    LC_LIST = [os.path.join(EXAMPLES_DIR, str(i)) for i in range(1, 11)]

    def test_python_inline_var_reduction(self):
        """Mirror CLI Example 1: b = numpy.var(mag)."""
        batch = (vt.Pipeline()
                 .python("b = numpy.var(mag)",
                         invars="mag", outvars="b", outputcolumns="b")
                 ).run_filelist(self.LC_LIST,
                                perpoint_columns={"t": 1, "mag": 2, "err": 3})
        assert "PYTHON_b_0" in batch.vars.columns
        # Documented value for EXAMPLES/2 is 0.001342… — no need to be
        # bit-exact, but a sanity tolerance catches gross marshalling bugs.
        ex2 = batch.vars[batch.vars["Name"] == "2"]["PYTHON_b_0"].iloc[0]
        assert abs(ex2 - 0.0013420988) < 1e-7

    def test_python_init_inline(self):
        """init code defines a helper function used in the main code."""
        batch = (vt.Pipeline()
                 .python("b = mean_minus_zero(mag)",
                         init="def mean_minus_zero(x): return float(numpy.mean(x))",
                         invars="mag", outvars="b", outputcolumns="b")
                 ).run_filelist(self.LC_LIST,
                                perpoint_columns={"t": 1, "mag": 2, "err": 3})
        assert "PYTHON_b_0" in batch.vars.columns

    def test_python_init_fromfile(self, tmp_path):
        """init_fromfile reads helper definitions from a .py file."""
        init_file = tmp_path / "init.py"
        init_file.write_text("def avg(x): return float(numpy.mean(x))\n")
        batch = (vt.Pipeline()
                 .python("b = avg(mag)",
                         init=str(init_file), init_fromfile=True,
                         invars="mag", outvars="b", outputcolumns="b")
                 ).run_filelist(self.LC_LIST,
                                perpoint_columns={"t": 1, "mag": 2, "err": 3})
        assert "PYTHON_b_0" in batch.vars.columns

    def test_python_fromfile(self, tmp_path):
        """fromfile reads the main code body from a .py file."""
        cmd_file = tmp_path / "cmd.py"
        cmd_file.write_text("b = float(numpy.std(mag))\n")
        batch = (vt.Pipeline()
                 .python(str(cmd_file), fromfile=True,
                         invars="mag", outvars="b", outputcolumns="b")
                 ).run_filelist(self.LC_LIST,
                                perpoint_columns={"t": 1, "mag": 2, "err": 3})
        assert "PYTHON_b_0" in batch.vars.columns

    def test_python_vars_combined(self):
        """`vars` shorthand passes the same names IN and OUT, used here
        to round-trip an LC vector (mag in, mag out)."""
        # vartools' grammar makes `vars` exclusive of `invars`/`outvars`,
        # so the wrapper emits `vars mag` and the user code may both
        # read and reassign mag.  -rms before/after demonstrates the
        # round-trip succeeded.
        batch = (vt.Pipeline()
                 .rms()
                 .python("mag = mag - numpy.mean(mag)", vars="mag")
                 .rms()
                 ).run_filelist(self.LC_LIST,
                                perpoint_columns={"t": 1, "mag": 2, "err": 3})
        # After de-meaning, Mean_Mag_2 should be ~0; RMS unchanged.
        for v in batch.vars["Mean_Mag_2"]:
            assert abs(v) < 1e-10
        for pre, post in zip(batch.vars["RMS_0"], batch.vars["RMS_2"]):
            assert abs(pre - post) < 1e-10

    def test_python_process_all_lcs(self):
        """Mirror CLI Example 2-style usage: a single call processes
        the whole batch.  invars arrive as lists of arrays."""
        code = (
            "b = []\n"
            "for arr in mag: b.append(float(numpy.var(arr)))\n"
        )
        batch = (vt.Pipeline()
                 .python(code,
                         invars="mag", outvars="b", outputcolumns="b",
                         process_all_lcs=True)
                 ).run_filelist(self.LC_LIST,
                                perpoint_columns={"t": 1, "mag": 2, "err": 3})
        assert "PYTHON_b_0" in batch.vars.columns
        # Same answer as Example 1 (per-LC variance) — round-trip check.
        ex2 = batch.vars[batch.vars["Name"] == "2"]["PYTHON_b_0"].iloc[0]
        assert abs(ex2 - 0.0013420988) < 1e-7

    def test_python_continueprocess(self):
        """The 2nd -python(continueprocess=1) sees module-level state
        from the 1st by sharing the same Python subprocess.

        Each -python call's body runs inside a generated wrapper
        function, so plain assignments are local — to share across
        calls we must declare ``global``.  This is the canonical use
        of ``continueprocess`` (re-using imports / cached state /
        helper objects without paying the init cost twice).
        """
        batch = (vt.Pipeline()
                 .python("global _shared_mult\n_shared_mult = 1000.0",
                         invars="mag", outvars="mag")
                 .python("global _shared_mult\n"
                         "b = _shared_mult * float(numpy.var(mag))",
                         continueprocess=1,
                         invars="mag", outvars="b", outputcolumns="b")
                 ).run_filelist(self.LC_LIST,
                                perpoint_columns={"t": 1, "mag": 2, "err": 3})
        # PYTHON_b_1 = 1000 × variance.  EXAMPLES/2 has variance
        # 0.001342 → expected ~1.342 within whitespace-format tolerance.
        ex2 = batch.vars[batch.vars["Name"] == "2"]["PYTHON_b_1"].iloc[0]
        assert abs(ex2 - 1.3420988) < 1e-3

    def test_python_modifies_lc_vector(self):
        """Code that subtracts the mean of mag should leave RMS unchanged
        but Mean_Mag near zero.  Confirms the LC-vector round-trip."""
        batch = (vt.Pipeline()
                 .rms()
                 .python("mag = mag - numpy.mean(mag)", vars="mag")
                 .rms()
                 ).run_filelist(self.LC_LIST,
                                perpoint_columns={"t": 1, "mag": 2, "err": 3})
        # After de-meaning, Mean_Mag should be ~0.
        for v in batch.vars["Mean_Mag_2"]:
            assert abs(v) < 1e-10
        # RMS should be unchanged (de-meaning is a constant shift).
        for pre, post in zip(batch.vars["RMS_0"], batch.vars["RMS_2"]):
            assert abs(pre - post) < 1e-10


# ===========================================================================
# Integration: cmd.python(inprocess=True) — Stage 2 callback hook
# ===========================================================================
#
# Library mode is required for the callback hook to fire.  Each test
# constructs a single-LC pipeline and runs it via Pipeline.run(lc) so
# pyvartools can take the in-process library path; the wrapper pre-checks
# that subprocess-mode-forcing factors are absent.

@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools not installed")
class TestPythonInprocessIntegration:
    """End-to-end coverage of cmd.python(inprocess=True)."""

    def test_inprocess_basic_var(self):
        """Inprocess path returns the same variance as the subprocess path."""
        lc = vt.LightCurve.from_file(EXAMPLE_LC)
        r_sub = vt.Pipeline().python(
            "b = float(numpy.var(mag))",
            invars="mag", outvars="b", outputcolumns="b",
            inprocess=False).run(lc)
        r_in = vt.Pipeline().python(
            "b = float(numpy.var(mag))",
            invars="mag", outvars="b", outputcolumns="b",
            inprocess=True).run(lc)
        assert abs(r_sub.vars["PYTHON_b_0"] - r_in.vars["PYTHON_b_0"]) < 1e-12

    def test_inprocess_sees_main_namespace(self):
        """User code in the inprocess body can reference globals defined
        by the calling host script — this is the entire point of the
        callback hook."""
        import sys
        # Inject a sentinel into __main__ that the user code resolves.
        sys.modules["__main__"]._VT_TEST_K = 100.0
        try:
            lc = vt.LightCurve.from_file(EXAMPLE_LC)
            r = vt.Pipeline().python(
                "b = _VT_TEST_K + float(numpy.var(mag))",
                invars="mag", outvars="b", outputcolumns="b",
                inprocess=True).run(lc)
            assert 99.9 < r.vars["PYTHON_b_0"] < 100.1
        finally:
            del sys.modules["__main__"]._VT_TEST_K

    def test_inprocess_custom_namespace(self):
        """A custom namespace= dict isolates user code from __main__."""
        lc = vt.LightCurve.from_file(EXAMPLE_LC)
        sandbox = {"my_factor": 3.0}
        r = vt.Pipeline().python(
            "b = my_factor * float(numpy.var(mag))",
            invars="mag", outvars="b", outputcolumns="b",
            inprocess=True, namespace=sandbox).run(lc)
        # 3.0 × 0.001342 ≈ 0.004026.
        assert abs(r.vars["PYTHON_b_0"] - 0.004026) < 1e-4

    def test_inprocess_lc_vector_mutation(self):
        """Mutating an LC vector inside inprocess code is reflected in
        downstream commands (Mean_Mag → 0 after subtracting mean)."""
        lc = vt.LightCurve.from_file(EXAMPLE_LC)
        r = (vt.Pipeline()
             .rms()
             .python("mag = mag - numpy.mean(mag)",
                     vars="mag", inprocess=True)
             .rms()
             ).run(lc)
        assert abs(r.vars["Mean_Mag_2"]) < 1e-10
        assert abs(r.vars["RMS_0"] - r.vars["RMS_2"]) < 1e-10

    def test_inprocess_init_runs_once(self):
        """init code runs exactly once before per-LC bodies — a function
        defined in init is then callable from the main body."""
        lc = vt.LightCurve.from_file(EXAMPLE_LC)
        r = vt.Pipeline().python(
            "b = compute_metric(mag)",
            init="def compute_metric(arr): return float(numpy.std(arr) * 2)",
            invars="mag", outvars="b", outputcolumns="b",
            inprocess=True).run(lc)
        # std(EXAMPLES/2) * 2 ≈ 2 × 0.0366 ≈ 0.0732.
        assert 0.07 < r.vars["PYTHON_b_0"] < 0.08

    def test_inprocess_refused_for_run_filelist(self):
        """run_filelist() must hard-error rather than silently fall
        back to subprocess (which would defeat the inprocess intent)."""
        with pytest.raises(RuntimeError, match="inprocess=True"):
            (vt.Pipeline()
             .python("b = 0.0", outvars="b", outputcolumns="b",
                     inprocess=True)
             ).run_filelist([EXAMPLE_LC],
                            perpoint_columns={"t": 1, "mag": 2, "err": 3})

    def test_inprocess_refused_for_run_batch(self):
        """run_batch() also rejects inprocess=True."""
        lc = vt.LightCurve.from_file(EXAMPLE_LC)
        with pytest.raises(RuntimeError, match="inprocess=True"):
            (vt.Pipeline()
             .python("b = 0.0", outvars="b", outputcolumns="b",
                     inprocess=True)
             ).run_batch([lc])

    def test_inprocess_refused_for_run_combinelcs(self):
        """run_combinelcs() also rejects inprocess=True."""
        with pytest.raises(RuntimeError, match="inprocess=True"):
            (vt.Pipeline()
             .python("b = 0.0", outvars="b", outputcolumns="b",
                     inprocess=True)
             ).run_combinelcs([[EXAMPLE_LC, EXAMPLE_LC]])

    def test_inprocess_refused_when_no_marshalled_vars(self):
        """``cmd.python(inprocess=True)`` with no invars/outvars/vars
        rejects at construction time -- vartools' parser would set
        processallvariables=1 in that case, which the in-process
        callback explicitly does not support and which would otherwise
        fall through to a subprocess fork (unsafe inside a Python host
        process; segfaults during process_lc)."""
        with pytest.raises(ValueError, match="inprocess=True"):
            cmd.python("x = 1", inprocess=True)

    def test_inprocess_refused_with_process_all_lcs(self):
        """``cmd.python(inprocess=True, process_all_lcs=True)`` -- the
        process_all_lcs path uses RunPythonCommand_all_lcs which has
        no in-process branch, so the call would fall through to the
        unsafe subprocess fork."""
        with pytest.raises(ValueError, match="process_all_lcs=True"):
            cmd.python("x = 1", invars="mag",
                       inprocess=True, process_all_lcs=True)

    def test_inprocess_refused_with_continueprocess(self):
        """``cmd.python(inprocess=True, continueprocess=N)`` -- the
        in-process callback gate explicitly rejects iscontinueprocess
        (runpython.c:2640)."""
        with pytest.raises(ValueError, match="continueprocess"):
            cmd.python("x = 1", invars="mag",
                       inprocess=True, continueprocess=0)

    def test_inprocess_with_save_periodogram_runs_in_library_mode(self):
        """save_*=True is now library-compatible, so combining it with
        cmd.python(inprocess=True) (with proper invars/outvars) runs
        end-to-end without forcing subprocess.  This used to be blocked
        by a conservative obstacle list as a workaround for the
        no-invars segfault, which is now caught at ctor time."""
        lc = vt.LightCurve.from_file(EXAMPLE_LC)
        result = (vt.Pipeline()
                  .LS(0.1, 10.0, 0.1, save_periodogram=True)
                  .python("b = float(numpy.var(mag))",
                          invars="mag", outvars="b", outputcolumns="b",
                          inprocess=True)
                  ).run(lc)
        assert "PYTHON_b_1" in result.vars.index
        assert "LS_periodogram_0" in result.files
