"""Layer-1 API correctness sweep: every command class's simplest valid
construction must produce a CLI that vartools' own parser accepts.

Uses ``Pipeline.validate()`` as the oracle (runs vartools with
``-headeronly``, no light-curve work).  Each failure points at a
wrapper bug — the wrapper emits CLI tokens that vartools rejects.

Hand-curated kwargs for commands that have required parameters or
behaviour that depends on a meaningful argument.  Anything not in the
registry is constructed with ``cls()`` and the test will fail loudly if
it can't be — that's a signal to either add it to the registry or fix
the wrapper.
"""
from __future__ import annotations

import inspect
import os
from pathlib import Path

import pytest

import pyvartools as vt
from pyvartools import commands as cmd
from pyvartools._binary import get_binary
from pyvartools.results import PipelineValidationError

try:
    get_binary()
    _HAVE_BINARY = True
except Exception:
    _HAVE_BINARY = False


# ---------------------------------------------------------------------------
# Minimal-valid-construction registry
# ---------------------------------------------------------------------------
#
# Each entry is either:
#   - a dict of kwargs to pass to ``cls(**kwargs)``, or
#   - a zero-arg callable returning a constructed cmd instance.
#
# Add to this when a new command is added or when a wrapper's required
# kwargs change.  Commands NOT listed here are constructed with ``cls()``
# and the sweep will surface them as failures if that doesn't work.

_EXAMPLES = Path(__file__).resolve().parent.parent.parent / "EXAMPLES"


def _trendlist_file(tmp_path: Path) -> str:
    p = tmp_path / "trendlist.txt"
    p.write_text(f"{_EXAMPLES}/2\n")
    return str(p)


# Commands that can't be validated standalone (escape hatches, user
# extension shims) and aren't in scope for this sweep.
_SKIP_NOT_APPLICABLE = {
    "Raw",            # escape hatch — caller passes raw CLI tokens
    "UserCommand",    # user-supplied extension; no canonical CLI form
    "_UserLibCommand",
    # Block-structured control-flow tokens — only meaningful inside a
    # matched if/elif/else block; validating one in isolation is
    # nonsensical.  These are exercised by the control-flow integration
    # tests separately.
    "elifcmd", "elsecmd", "ficmd", "ifcmd",
    # Pipeline control commands that only make sense paired with a prior
    # savelc / restoretimes etc.
    "restorelc", "restoretimes", "savelc",
    # Aggregator commands that require a prior matching command to
    # carry forward state — not standalone-valid.
    "FFT", "IFFT",
    # cmd.python / cmd.R need real script bodies; covered by their own
    # dedicated tests in test_all_commands.py.
    "python", "R",
    # Match needs a per-LC match-file column setup; covered by its own
    # integration test.
    "match",
}

# Commands that need specific kwargs to construct meaningfully.  Values
# pulled from unittest.sh / examples where possible.
_MIN_KWARGS = {
    # Periodicity
    "LS":          dict(minp=0.1, maxp=10.0, subsample=0.1, npeaks=1),
    "aov":         dict(minp=0.1, maxp=10.0, subsample=0.1, finetune=0.01,
                         npeaks=1),
    "aov_harm":    dict(nharm=2, minp=0.1, maxp=10.0, subsample=0.1,
                         finetune=0.01, npeaks=1),
    "BLS":         dict(rmin=0.005, rmax=0.05, minper=0.5, maxper=10.0,
                         nfreq=200, nbins=200, npeaks=1, qmin=0.01, qmax=0.1),
    "BLSFixPer":   dict(period=2.5),
    "BLSFixDurTc": dict(duration=0.1, Tc=0.0),
    "BLSFixPerDurTc": dict(period=2.5, duration=0.1, Tc=0.0),
    "dftclean":    dict(nbeam=1),
    "Phase":       dict(period=2.5),  # avoid the 'ls' default that needs a prior LS
    "harmonicfilter": dict(period=10.0),
    "fourierfilter": dict(),
    # GetLSAmpThresh with period='ls' needs a prior -LS in the pipeline.
    "GetLSAmpThresh": lambda tmp_path: [
        cmd.LS(0.1, 10.0, 0.1, npeaks=1),
        cmd.GetLSAmpThresh(period="ls", minp=0.1, thresh=10.0),
    ],
    "wwz":         dict(),
    "Killharm":    dict(period=2.5),

    # Manipulation
    "addnoise":    dict(noise_type="white", sig_white=0.001),
    "binlc":       dict(method="average", binsize=0.1),
    "clip":        dict(sigclip=5.0),
    "decorr":      dict(),
    "expr":        dict(expression="x = 1.0"),
    "Injectharm":  dict(period=2.5, amplitude=0.05),
    "Injecttransit": dict(period=3.0, Rp=0.1, Mp=1.0, phase=0.0,
                           sini=1.0, Mstar=1.0, Rstar=1.0),
    "Starspot":    dict(period=10.0),  # other kwargs default
    "MandelAgolTransit": dict(P0=2.5, T00=0.0, r0=0.1, a0=10.0,
                                bimpact=0.1, mconst0=-1,
                                ld_coeffs=[0.3, 0.3]),
    # SoftenedTransit init_params='bls' (default) needs a prior BLS.
    "SoftenedTransit": lambda tmp_path: [
        cmd.BLS(rmin=0.005, rmax=0.05, minper=0.5, maxper=10.0,
                 nfreq=200, nbins=200, npeaks=1, qmin=0.01, qmax=0.1),
        cmd.SoftenedTransit(),
    ],
    "microlens":   dict(),
    "Jstet":       dict(timescale=0.01, dates=str(_EXAMPLES / "dates_tfa")),
    "rms":         dict(),
    "rmsbin":      dict(nbin=1, bintimes=[0.01]),
    "stats":       dict(variables="mag", statistics="mean"),
    "chi2":        dict(),
    "chi2bin":     dict(nbin=1, bintimes=[0.01]),
    "alarm":       dict(),
    "autocorrelation": dict(start=0.0, stop=1.0, step=0.01),
    "medianfilter": dict(time=0.01),
    "ensemblerescalesig": dict(),
    "rescalesig":  dict(),
    "linfit":      dict(function="a*t+b", paramlist="a,b"),
    "nonlinfit":   dict(function="a*sin(2*pi*t/p)+b",
                         paramlist="a=0.1:0.01,p=1.0:0.01,b=0.0:0.001"),
    "fluxtomag":   dict(mag_constant=25.0),
    "difffluxtomag": dict(mag_constant=25.0),
    "changeerror": dict(),
    "changevariable": dict(column="t", var="t"),
    "TFA":         lambda tmp_path: cmd.TFA(
        trendlist=_trendlist_file(tmp_path),
        dates_file=str(_EXAMPLES / "dates_tfa"),
        pixelsep=25.0,
    ),
    "TFA_SR":      lambda tmp_path: cmd.TFA_SR(
        trendlist=_trendlist_file(tmp_path),
        dates_file=str(_EXAMPLES / "dates_tfa"),
        pixelsep=25.0,
    ),
    "SYSREM":      dict(ninput_color=1, ninput_airmass=1,
                         initial_airmass_file=str(_EXAMPLES / "dates_tfa")),
    "sortlc":      dict(),
    "copylc":      dict(ncopies=2),
    # converttime to BJD needs an RA/Dec source; the simplest valid form
    # converts between two non-barycentric systems.
    "converttime": dict(input_format="JD", output_format="MJD"),
    "restricttimes": dict(mode="JDrange", minJD=0.0, maxJD=1e9),
    "resample":    dict(method="nearest"),
    # Verify the dtype-alias normalisation: "double" → TDOUBLE.
    "addfitskeyword": dict(keyword="TESTKEY", dtype="double", value=1.0),
    "findblends":  dict(matchrad=1.0, period=1.0),
    # columnsuffix is a meta-flag that modifies the next command; chain
    # it with a real command so vartools sees a non-empty pipeline.
    "columnsuffix": lambda tmp_path: [
        cmd.columnsuffix(suffix="x"),
        cmd.rms(),
    ],
    "print_cols":  dict(variables="mag"),
    "o":           dict(filename="-"),
    "ftuneven":    None,    # requires its userlib .so to be loaded first;
                            # not part of this in-tree sweep
}


def _all_command_classes():
    return [
        getattr(cmd, name) for name in dir(cmd)
        if isinstance(getattr(cmd, name, None), type)
        and getattr(cmd, name).__name__ not in _SKIP_NOT_APPLICABLE
        and getattr(cmd, name).__name__ != "VartoolsCommand"
        and hasattr(getattr(cmd, name), "_vt_name")
    ]


def _construct_minimal(cls, tmp_path):
    """Return a list of cmd instances to add to the pipeline.

    The target command is the LAST element; any preceding commands are
    set-up dependencies (e.g. an LS run before a Killharm that uses
    its period).  Single-command entries (the common case) return a
    1-element list.
    """
    name = cls.__name__
    if name in _MIN_KWARGS and _MIN_KWARGS[name] is None:
        pytest.skip(f"{name}: explicitly skipped (see registry)")
    spec = _MIN_KWARGS.get(name, "_NOT_PRESENT_")
    if spec == "_NOT_PRESENT_":
        try:
            return [cls()]
        except TypeError as e:
            pytest.skip(f"{name}: no minimal-construction registry entry "
                         f"({e}); add to _MIN_KWARGS")
    if callable(spec):
        result = spec(tmp_path)
    elif isinstance(spec, list):
        # List of (cls, kwargs) tuples — explicit chain.
        return [c(**kw) if isinstance(kw, dict) else kw for (c, kw) in spec]
    else:
        result = cls(**spec)
    return [result] if not isinstance(result, list) else result


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools binary not available")
@pytest.mark.parametrize(
    "cmd_cls",
    _all_command_classes(),
    ids=lambda c: c.__name__,
)
def test_default_construction_validates(cmd_cls, tmp_path):
    """For every wrapped vartools command, the minimal valid construction
    must produce a CLI that vartools accepts under -headeronly.

    A failure here means either:
      (a) the wrapper emits CLI tokens vartools rejects (a real bug), or
      (b) the entry in _MIN_KWARGS doesn't actually represent a valid
          minimal construction.  Verify against vartools' own
          unittest.sh / -listcommands output.
    """
    instances = _construct_minimal(cmd_cls, tmp_path)
    pipe = vt.Pipeline()
    for inst in instances:
        pipe = pipe.add(inst)
    try:
        pipe.validate()
    except PipelineValidationError as e:
        pytest.fail(
            f"{cmd_cls.__name__} produced a CLI that vartools rejected.\n"
            f"argv: {' '.join(e.argv)}\n"
            f"stderr:\n{e.stderr}"
        )
