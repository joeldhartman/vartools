# VARTOOLS — Developer Guide for Claude

## Working agreement (read first)

These are **durable** preferences that override default tooling
behaviour.  They are also captured in the persistent memory store
(`feedback_no_push_without_approval.md`,
`feedback_confirm_before_big_project.md`); duplicated here so they
survive context compaction within a session.

* **Do not push to GitHub without explicit approval.**  Accumulate
  commits locally on the working branch (currently `develop`) and stop.
  Only run `git push` when the user says "push" (or equivalent).
  Local commits are fine — they're easy to amend or drop.  A push
  makes the change visible on GitHub Actions / collaborators.

* **Confirm before starting a multi-step development project.**
  Anything that introduces a new class / module, touches more than two
  files, or that you'd describe with bullet points ("Plan: 1… 2… 3…")
  should be paused on for explicit user approval before implementation
  begins.  A short proposal ("here's what I'd add, here's what I'd
  change, OK to proceed?") is preferred over diving in.  Quick fixes,
  bug repairs, single-file edits, and continuations of an in-progress
  task do not need this gate.

* **Prefer a sketch before code.**  When the path is non-obvious or
  the user's intent isn't yet 100% clear, write the design as a short
  prose paragraph or bullet list first.  Faster to course-correct on
  three sentences than on three files.

## Repository layout

```
configure.ac          autoconf input; regenerate configure with: autoreconf -f
Makefile.am           automake input
src/                  C source for the vartools binary and libvartools
python/               pyvartools Python package (see python/README.md)
  pyvartools/         package source
    commands/         one module per command family (periodicity, fitting, …)
    pipeline.py       Pipeline class and run methods
    lightcurve.py     LightCurve wrapper
    results.py        Result / BatchResult / parse_oneline_output
    _binary.py        binary discovery
    _libpipeline.py   in-process library fast-path (optional)
  tests/              pytest suite
  API_AUDIT.md        gap tracker — read before adding/changing commands
  pyproject.toml      package metadata; requires Python ≥ 3.8
EXAMPLES/             example light curves and expected outputs
USERFUNCS/            user-defined function extension mechanism
USERLIBS/             user-defined command extension mechanism
unittest.sh           vartools CLI regression tests (requires built binary)
.github/workflows/    GitHub Actions CI
```

## Authoritative source

The SVN working copy at `/home/jhartman/SVN/HATreduc/HATpipe/source/vartools`
is the primary source of truth.  The Git repo at `/home/jhartman/GIT/vartools`
is kept in sync manually (rsync from SVN → Git, then push).  Do **not** edit
only one side and forget to sync.

The website docs live separately at `/home/jhartman/src/vartools-site` (MkDocs).

## Build

```bash
# Standard build (from the vartools source root)
./configure
make

# Optional dependencies surfaced by configure:
#   --with-cfitsio    FITS file support
#   --with-RHOME      R integration (callrunR)
#   --with-python     Python scripting in vartools (callrunpython)
#   --with-cspice     NAIF SPICE toolkit

# After building, install the Python package:
make install
make install-python          # records install path into pyvartools
# or, for development:
cd python && pip install -e .
```

After editing `configure.ac`, regenerate with:
```bash
autoreconf -f
```

## Python package (pyvartools)

### Running tests

```bash
# From python/ directory, using the pyvartools conda environment:
conda run -n pyvartools python -m pytest tests/ -q

# Run a specific subset:
conda run -n pyvartools python -m pytest tests/test_all_commands.py -k "nonlinfit" -q
```

Tests that require a built vartools binary are skipped automatically when the
binary is not installed (`@pytest.mark.skipif(not _HAVE_BINARY, ...)`).

### Key conventions

- **`_build_cmd()`** in `pipeline.py` assembles the full CLI argument list.
  Global flags (`-randseed`, `-skipmissing`, `-jdtol`, `-matchstringid`) are
  emitted after `-parallel N` and before command tokens.
- **Library-mode fast-path**: `run()` and `run_batch()` can call
  `libvartoolspipeline.so` in-process.  Any non-default global option forces
  subprocess mode — guard with `_has_global_opts`.
- **`PerLC`**: wraps a list of per-light-curve values for batch runs.  Commands
  that receive a `PerLC` value emit the value inline per LC rather than a
  constant.  `run_combinelcs()` and `run_file()` must raise `ValueError` if any
  command has PerLC attributes.
- **Output file specs**: commands that write auxiliary files (periodograms,
  model LCs, MCMC chains, …) implement `_output_file_specs()` returning a dict
  of `key → (suffix, parser)`.  Check for accidental duplicate method
  definitions — Python silently uses the last one.
- **Column naming**: `COMMANDNAME_descriptor_suffix` where suffix is the
  0-indexed command position in the pipeline, or an explicit string set by
  `columnsuffix()`.

### Adding a new command

1. Add a class in the appropriate `commands/*.py` module, inheriting from
   `VartoolsCommand`.
2. Implement `_cli_tokens(self)` to return the list of CLI tokens.
3. If the command writes output files, implement `_output_file_specs(self)`.
4. Export the class from `commands/__init__.py`.
5. Add tests in `tests/test_all_commands.py`.
6. Update `API_AUDIT.md` to mark the gap closed.

## Sanitizer test pass (`make check-asan`)

`make check-asan` rebuilds vartools with `-fsanitize=address,undefined` and runs
`unittest.sh` against the instrumented binary.  Each ASan/UBSan hit gets written
to `asan-report.<pid>` files in the build root; the target exits non-zero if any
appear.

```bash
make check-asan                 # ~5 min full sanitizer pass
ls asan-report.*                # inspect any reports written
make clean && make              # restore the regular optimized build
```

When to run it:
- Before tagging a release.
- After editing memory-sensitive code (parsers, allocators, numerical loops,
  anything in `parselc.c`/`parseinputlist.c`/`runpython.c`/`runR.c`).
- When investigating a suspected memory bug.

When NOT to run it:
- Routine dev cycles — sanitizer-instrumented runs are 5–10x slower than the
  un-instrumented suite.  Use the regular `unittest.sh` path for that.
- In per-push CI — if/when CI gains an ASan job, wire it as a separate
  workflow gated on weekly schedule or release tags, not the default push.

The target uses `AM_CFLAGS`/`AM_LDFLAGS` overrides so the configured
`./configure` CFLAGS are preserved (anaconda includes, march flags, etc.) with
sanitizer flags layered on top.  `halt_on_error=0` and `abort_on_error=0` keep
the run going past the first hit so one pass surfaces all issues.

After a `check-asan` run, the build tree is in sanitizer mode.  Run `make clean
&& make` to return to the normal optimized binary before doing further dev
work.

## CI (GitHub Actions)

Workflow: `.github/workflows/test-compile-vartools-ubuntu.yml`

The workflow installs `gfortran`, `libcfitsio-dev`, and `r-base-dev`, then runs
`./configure && make`.  `r-base-dev` is required (not just `r-base`) because
`configure` probes for `libR.so`.

## Commit / sync workflow

1. Edit source in SVN (`/home/jhartman/SVN/HATreduc/HATpipe/source/vartools`).
2. `svn commit` from SVN working copy.
3. Rsync to Git repo, then `git add` / `git commit` / `git push` to GitHub.
   Generated files (`configure`, `Makefile.in`, `aclocal.m4`, …) **are**
   tracked in Git and should be synced.
4. Do not commit build artifacts (`config.h`, `Makefile`, `.deps/`, binaries,
   `*.o`, `*.la`), local-only files (`hatconf_jhlaptop.sh`), or tarballs.
