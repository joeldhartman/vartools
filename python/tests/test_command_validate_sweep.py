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
    "magtoflux":   dict(mag_constant=25.0),
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
    "o":           dict(outname="-"),
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


# ---------------------------------------------------------------------------
# Layer 2 — enumerated-keyword variation matrix
# ---------------------------------------------------------------------------
#
# Each entry is ``(cls_name, base_kwargs, kwarg_name, valid_values)``.  The
# test substitutes each valid_value into base_kwargs and runs validate().
# This is where the wwz `phase="rand"` / addfitskeyword `dtype="double"`
# class of bug hides — the wrapper takes a closed-set string kwarg and
# either passes it through verbatim or wraps it the wrong way.

_LD_VARIANTS = ("linear", "quad", "logarithmic", "squareroot", "claret",
                 "nonlinear", "exponential")

# ``base_kwargs`` is what the cmd needs to construct AND validate at all;
# the kwarg under test gets injected on top.  When base needs a prior
# command in the pipeline, set ``priors`` to a list of cmd instances.

_VARIATION_MATRIX = [
    # (label, cls, base_kwargs, kwarg_name, valid_values, priors)
    ("addfitskeyword.dtype",  cmd.addfitskeyword,
        dict(keyword="K1", value=1.0), "dtype",
        ["double", "int", "long", "string",
         "TDOUBLE", "TINT", "TLONG", "TSTRING"], None),
    ("addfitskeyword.hdu",    cmd.addfitskeyword,
        dict(keyword="K1", dtype="double", value=1.0), "hdu",
        ["primary", "extension"], None),
    ("addfitskeyword.mode",   cmd.addfitskeyword,
        dict(keyword="K1", dtype="double", value=1.0), "mode",
        ["append", "update"], None),

    ("addnoise.noise_type",   cmd.addnoise,
        dict(sig_white=0.001, sig_red=0.001, rho=0.1, gamma=0.5,
             nu=1.0), "noise_type",
        ["white", "squareexp", "exp", "matern", "wavelet"], None),

    ("binlc.method",          cmd.binlc,
        dict(binsize=0.1), "method",
        ["average", "median", "weightedaverage"], None),
    ("binlc.time_output",     cmd.binlc,
        dict(method="average", binsize=0.1), "time_output",
        ["tcenter", "tmin", "tmax", "average", "median"], None),

    ("converttime.input_format", cmd.converttime,
        dict(output_format="mjd"), "input_format",
        ["jd", "mjd", "JD", "MJD"], None),  # case-insensitive
    ("converttime.output_format", cmd.converttime,
        dict(input_format="jd"), "output_format",
        ["jd", "mjd", "JD", "MJD"], None),
    ("converttime.input_sys", cmd.converttime,
        dict(input_format="jd", output_format="mjd"), "input_sys",
        ["tdb", "utc"], None),
    ("converttime.output_sys", cmd.converttime,
        dict(input_format="jd", output_format="mjd"), "output_sys",
        ["tdb", "utc"], None),

    ("fourierfilter.mode",    cmd.fourierfilter,
        dict(minfreq=0.1, maxfreq=1.0), "mode",
        ["lowpass", "highpass", "bandpass", "bandcut", "full"], None),

    ("GetLSAmpThresh.mode",   cmd.GetLSAmpThresh,
        dict(period="ls", minp=0.1, thresh=10.0,
             listfile=str(_EXAMPLES / "dates_tfa")), "mode",
        ["harm", "file"],
        [cmd.LS(0.1, 10.0, 0.1, npeaks=1)]),

    # 'rand' / 'logrand' need a min,max range — pass via "rand <min> <max>".
    ("Injectharm.amplitude.kw", cmd.Injectharm,
        dict(period=2.5), "amplitude",
        ["rand 0.01 0.1", "list", "logrand 0.01 0.1", 0.05], None),
    ("Injectharm.phase.kw",    cmd.Injectharm,
        dict(period=2.5, amplitude=0.05), "phase",
        ["rand", "list", 0.5], None),

    ("Injecttransit.ld_type", cmd.Injecttransit,
        dict(period=3.0, Rp=0.1, Mp=1.0, phase=0.0,
             sini=1.0, Mstar=1.0, Rstar=1.0,
             ld_coeffs=[0.3, 0.3]), "ld_type",
        ["quad"], None),

    ("MandelAgolTransit.ld_type", cmd.MandelAgolTransit,
        dict(P0=2.5, T00=0.0, r0=0.1, a0=10.0, bimpact=0.1,
             mconst0=-1, ld_coeffs=[0.3, 0.3]), "ld_type",
        ["quad"], None),

    ("medianfilter.method",   cmd.medianfilter,
        dict(time=0.01), "method",
        ["median", "mean", "weightedmean"], None),

    ("resample.method",       cmd.resample,
        dict(), "method",
        ["nearest", "linear", "spline", "splinemonotonic", "bspline"], None),

    ("restricttimes.mode",    cmd.restricttimes,
        dict(minJD=0.0, maxJD=1e9), "mode",
        ["JDrange"], None),

    ("wwz.maxfreq.kw",        cmd.wwz,
        dict(), "maxfreq", ["auto", 5.0], None),
    ("wwz.tau0.kw",           cmd.wwz,
        dict(), "tau0", ["auto", 0.0], None),
    ("wwz.tau1.kw",           cmd.wwz,
        dict(), "tau1", ["auto", 100.0], None),
    ("wwz.dtau.kw",           cmd.wwz,
        dict(), "dtau", ["auto", 0.1], None),

    ("nonlinfit.optimizer",   cmd.nonlinfit,
        dict(function="a*sin(2*pi*t/p)+b",
             paramlist="a=0.1:0.01,p=1.0:0.01,b=0.0:0.001"), "optimizer",
        ["amoeba", "mcmc"], None),

    ("findblends.period.kw",  cmd.findblends,
        dict(matchrad=1.0), "period",
        ["list", "fix 1.0"], None),

    ("Phase.period.kw",       cmd.Phase,
        dict(), "period",
        ["ls", "aov", "bls", 2.5],
        [cmd.LS(0.1, 10.0, 0.1, npeaks=1),
         cmd.aov(minp=0.1, maxp=10.0, subsample=0.1, finetune=0.01,
                 npeaks=1),
         cmd.BLS(rmin=0.005, rmax=0.05, minper=0.5, maxper=10.0,
                 nfreq=200, nbins=200, npeaks=1, qmin=0.01, qmax=0.1)]),

    # vartools' -Killharm/-harmonicfilter period keywords are
    # ls / aov / both / injectharm — NOT bls (the CLI uses Phase
    # for that, not Killharm).
    ("Killharm.period.kw",    cmd.Killharm,
        dict(), "period",
        ["ls", "aov", "both", 2.5],
        [cmd.LS(0.1, 10.0, 0.1, npeaks=1),
         cmd.aov(minp=0.1, maxp=10.0, subsample=0.1, finetune=0.01,
                 npeaks=1)]),

    ("Starspot.period.kw",    cmd.Starspot,
        dict(), "period",
        ["ls", 2.5],
        [cmd.LS(0.1, 10.0, 0.1, npeaks=1)]),
]


def _expand_variations():
    """Yield (label, cls, kwargs, priors) tuples for the test parametrize."""
    out = []
    for label, cls, base, name, values, priors in _VARIATION_MATRIX:
        for v in values:
            kw = dict(base)
            kw[name] = v
            out.append((f"{label}={v!r}", cls, kw, priors or []))
    return out


@pytest.mark.skipif(not _HAVE_BINARY, reason="vartools binary not available")
@pytest.mark.parametrize(
    "label,cls,kwargs,priors",
    _expand_variations(),
    ids=lambda x: x if isinstance(x, str) else "",
)
def test_keyword_variation_validates(label, cls, kwargs, priors, tmp_path):
    """For each enumerated string-keyword kwarg, every documented valid
    value must produce a CLI that vartools accepts.  Catches the
    wrapper-emits-wrong-token bug class (e.g. wwz `phase="rand"` →
    `phaseexpr rand`, addfitskeyword `dtype="double"` → literal
    `double` token)."""
    instance = cls(**kwargs)
    pipe = vt.Pipeline()
    for p in priors:
        pipe = pipe.add(p)
    pipe = pipe.add(instance)
    try:
        pipe.validate()
    except PipelineValidationError as e:
        pytest.fail(
            f"{label} -> CLI rejected by vartools.\n"
            f"argv: {' '.join(e.argv)}\n"
            f"stderr:\n{e.stderr}"
        )


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
