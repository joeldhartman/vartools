"""Cross-implementation consistency: vartools -PDM vs cuvarbase CPU helpers.

cuvarbase (https://github.com/johnh2o2/cuvarbase) is a third-party Python
package whose CPU helpers implement the same Stellingwerf-1978 / Schwarzenberg-
Czerny PDM statistic that vartools' -PDM computes in C.  This test confirms
the two implementations agree at the same trial period on the same light
curve.

WHEN TO RUN
-----------

This suite is *opt-in* because each test invokes vartools through a
subprocess that runs a periodogram sweep on a multi-thousand-point light
curve, which makes the suite slow (minutes wall-clock).  It is intended
to be run when the PDM algorithm code (`src/pdm.c` or its dependencies)
changes -- not on every pytest invocation.

To run::

    RUN_CUVARBASE_TESTS=1 conda run -n pyvartools \\
        python -m pytest python/tests/test_cuvarbase_consistency.py -v

Without RUN_CUVARBASE_TESTS=1, all tests in this module are skipped.

Also skipped when cuvarbase isn't importable (e.g. no GPU / package not
installed in the active environment).

WHAT IT CHECKS
--------------

Confirms that vartools' theta values match cuvarbase's `pdm2_single_freq`
/ `binless_pdm_cpu` on EXAMPLES/2 at the canonical signal period.

Tolerances:

  binless tophat / gauss : ~1e-6 (bit-exact -- no bin-grid dependence)
  binned  step / linterp : ~1e-4 once vartools is given the same time
                           origin via `-expr 't=t-mean(t)'`.  Without the
                           t-shift the two disagree by ~10% because of
                           bin-grid alignment differences (cuvarbase
                           internally does `t -= np.mean(t)` before the
                           phase fold; vartools does not by default).
                           A separate `test_binned_disagrees_without_t_shift`
                           guards against silently regressing the t-shift
                           workaround.

Multicover has no cuvarbase reference -- it is exercised in unittest.sh.
"""

from __future__ import annotations

import os
import re
import subprocess

import numpy as np
import pytest

if os.environ.get("RUN_CUVARBASE_TESTS") != "1":
    pytest.skip(
        "cuvarbase consistency suite is slow (minutes wall-clock); set "
        "RUN_CUVARBASE_TESTS=1 to enable.  Run this when pdm.c or its "
        "dependencies change.",
        allow_module_level=True,
    )

cuvarbase_pdm = pytest.importorskip("cuvarbase.pdm")


_HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLES_LC = os.path.realpath(os.path.join(_HERE, "..", "..", "EXAMPLES", "2"))

if not os.path.exists(EXAMPLES_LC):
    pytest.skip(
        f"EXAMPLES/2 light curve not found at {EXAMPLES_LC}",
        allow_module_level=True,
    )

VARTOOLS = os.environ.get("VARTOOLS_BINARY", "vartools")

# Canonical dominant signal period in EXAMPLES/2; matches the unittest.sh
# regression tests and the value cuvarbase recovers as the periodogram peak.
P_ANCHOR = 1.23484930

# Narrow period range for the vartools subprocess -- we only need
# fixperiodSNR at P_ANCHOR; the surrounding sweep is overhead.  Binless
# variants are O(N^2) per trial period so the range matters a lot.
PMIN_BINNED, PMAX_BINNED = "1.0", "1.5"
PMIN_BINLESS, PMAX_BINLESS = "1.2", "1.3"


def _load_lc():
    d = np.loadtxt(EXAMPLES_LC)
    return d[:, 0].copy(), d[:, 1].copy(), d[:, 2].copy()


def _vartools_theta_at_fixperiod(variant, period, t_shift, *,
                                 Nbin=None, dphi=None, binless=False):
    """Run -PDM ... fixperiodSNR fix <period> and return theta at that
    period, parsed from the -oneline output.  The periodogram sweep range
    is narrowed around P_ANCHOR for speed; fixperiodSNR's theta does not
    depend on the sweep range, so this is a pure runtime optimisation."""
    cmd = [VARTOOLS, "-i", EXAMPLES_LC]
    if t_shift:
        cmd += ["-expr", "t=t-mean(t)"]
    cmd += ["-PDM", variant]
    if Nbin is not None:
        cmd += ["Nbin", str(Nbin)]
    if dphi is not None:
        cmd += ["dphi", str(dphi)]
    pmin = PMIN_BINLESS if binless else PMIN_BINNED
    pmax = PMAX_BINLESS if binless else PMAX_BINNED
    cmd += [
        pmin, pmax, "0.1", "0.01", "1", "0",
        "noerr",
        "fixperiodSNR", "fix", f"{period:.10f}",
        "-oneline",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, check=True)
    m = re.search(r"PDM_Theta_PeriodFix_\d+\s*=\s*(\S+)", r.stdout)
    if m is None:
        raise AssertionError(
            f"PDM_Theta_PeriodFix not found in vartools output.\n"
            f"Command: {' '.join(cmd)}\nstdout:\n{r.stdout}\nstderr:\n{r.stderr}"
        )
    return float(m.group(1))


# ---------------------------------------------------------------------------
# Binless variants: bit-exact agreement (no bin-grid dependence)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "variant,tophat,dphi",
    [("tophat", True, 0.05), ("gauss", False, 0.05)],
)
def test_binless_consistency(variant, tophat, dphi):
    t, mag, sig = _load_lc()
    N = len(t)
    w = np.ones(N) / N            # uniform weights to match -PDM noerr
    freq = 1.0 / P_ANCHOR

    t_c, mag_c = t.copy(), mag.copy()
    theta_cv = 1.0 - cuvarbase_pdm.binless_pdm_cpu(
        t_c, mag_c, w, [freq], dphi=dphi, tophat=tophat
    )[0]
    theta_vt = _vartools_theta_at_fixperiod(variant, P_ANCHOR, t_shift=False,
                                            dphi=dphi, binless=True)

    assert abs(theta_cv - theta_vt) < 1e-5, (
        f"{variant} dphi={dphi} disagreement: "
        f"cuvarbase={theta_cv:.6f}, vartools={theta_vt:.6f}, "
        f"|diff|={abs(theta_cv - theta_vt):.2e}"
    )


# ---------------------------------------------------------------------------
# Binned variants: agree once vartools is given the same time origin
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "variant,linterp,nbins",
    [("step", False, 8), ("linterp", True, 8)],
)
def test_binned_consistency_with_t_shift(variant, linterp, nbins):
    """vartools binned theta matches cuvarbase to ~1e-4 once the time
    origin is centred via -expr 't=t-mean(t)'.  cuvarbase implicitly does
    `t -= np.mean(t)` inside pdm2_single_freq; the t-shift makes vartools
    use the same phase grid."""
    t, mag, sig = _load_lc()
    N = len(t)
    w = np.ones(N) / N
    freq = 1.0 / P_ANCHOR

    t_c, mag_c = t.copy(), mag.copy()
    theta_cv = 1.0 - cuvarbase_pdm.pdm2_single_freq(
        t_c, mag_c, w, freq, nbins=nbins, linterp=linterp
    )
    theta_vt = _vartools_theta_at_fixperiod(variant, P_ANCHOR, t_shift=True,
                                            Nbin=nbins)

    assert abs(theta_cv - theta_vt) < 1e-4, (
        f"{variant} Nbin={nbins} disagreement (with -expr t=t-mean(t)): "
        f"cuvarbase={theta_cv:.6f}, vartools={theta_vt:.6f}, "
        f"|diff|={abs(theta_cv - theta_vt):.2e}"
    )


def test_binned_disagrees_without_t_shift():
    """Sanity check: without the -expr t=t-mean(t) workaround, vartools'
    binned theta differs from cuvarbase's by ~10% because the bin grids
    are offset by `mean(t)*freq mod 1`.  Both implementations are correct
    (they just use different phase origins) -- this is precisely the
    bin-edge sensitivity that the multicover variant is designed to
    average over.

    This test exists so a future refactor that internally centres t in
    vartools (eliminating the workaround) is flagged here for review; it
    should be updated to a matching assertion in that case rather than
    silently removed."""
    t, mag, sig = _load_lc()
    N = len(t)
    w = np.ones(N) / N
    freq = 1.0 / P_ANCHOR

    t_c, mag_c = t.copy(), mag.copy()
    theta_cv = 1.0 - cuvarbase_pdm.pdm2_single_freq(
        t_c, mag_c, w, freq, nbins=8, linterp=False
    )
    theta_vt_raw = _vartools_theta_at_fixperiod("step", P_ANCHOR,
                                                t_shift=False, Nbin=8)

    diff = abs(theta_cv - theta_vt_raw)
    assert diff > 1e-3, (
        f"Expected binned theta to disagree without t-shift; got "
        f"|diff|={diff:.2e}.  If vartools now centres t internally, "
        f"update this test."
    )
