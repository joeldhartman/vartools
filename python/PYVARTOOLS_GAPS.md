# pyvartools API gaps surfaced while writing runnable Python examples

These are CLI features that have to fall through to `cmd.Raw` or `subprocess.run`
because there is no typed Python kwarg for them yet:


## `SYSREM` — `save_trends` produces invalid CLI arguments

The CLI form is `otrends [trend_outfile]` where `trend_outfile` is a single
**file path**, but `cmd.SYSREM(save_trends=...)` treats it like the
directory-style `Output` spec used by `save_model`.  Setting
`save_trends="path/to/file"` raises `FileExistsError` (pyvartools tries
`os.makedirs` on the file path); `save_trends=True` emits a directory path
that vartools then rejects with "Invalid command or option ...".

**Workaround:** call `-SYSREM` via `cmd.Raw(["-SYSREM", ...])` if the
trend-vector file is needed.

**Fix needed:** make `save_trends` accept a file path and emit it as
`otrends <path>` without `os.makedirs`.


## `TFA_SR` — `xycol` is emitted in the wrong position

`cmd.TFA_SR(..., xycol=(x, y))` emits the `xycol` keyword *before* the
positional `pixelsep` value (see `commands/fitting.py:192-194` in the
TFA_SR `_to_cli_args` method).  The CLI requires
`pixelsep ["xycol" colx coly]` in that order, so any pipeline that uses
the typed wrapper with `xycol` set fails with `Invalid command or option`.

**Workaround:** call `-TFA_SR` via `cmd.Raw([...])` or `subprocess.run([...])`.

**Fix needed:** swap the emission order so `pixelsep` is appended before
the `xycol` block, mirroring the layout already used by `TFA._to_cli_args`.



## `load_userlib` resolved symlinked .so paths to versioned filenames

`pyvartools/userlib.py:196` previously called `Path(lib_path).resolve()`,
which follows the `stitch.so → stitch.so.0.0.0` symlink and then passes the
versioned filename to vartools' `-L` flag.  vartools' libtool dlopen needs
the basename without the version suffix to find `<libbasename>_Initialize`,
so that flow failed with "Error - Unspecified Error".

**Status:** fixed in this session by changing `.resolve()` →
`.absolute()` so the symlink target stays as written.

## `inputlcformat` with non-default column types

The `columns=` kwarg of `Pipeline.run_filelist()` only emits `name:colN` —
it doesn't accept type tags like `:string` or `:utc`.  When a light curve
has a non-numeric column (e.g. the fiphot string flag in
`-hatpiflag`'s example), the wrapper currently has no way to express the
column type, so callers fall through to `cmd.Raw(["-inputlcformat", "..."])`.

**Fix needed:** extend `_inputlcformat_from_spec` to accept a per-column
spec object (similar to `ListVar`) that carries `type=` and `format=`
fields, then emit `name:col:type[:fmt]` tokens.

## `combinelcs` mode — closed

`pyvartools` now exposes the full `-l … combinelcs` workflow:

  - `Pipeline.run_combinelcs(groups, ...)` — multi-group form.
  - `Pipeline.run_combinelc(files, ...)` — single-group convenience wrapper
    that returns a `Result`.
  - `Pipeline.run_filelist(paths, combinelcs=True, ...)` — drives a
    user-supplied list file (one comma-joined line per group).
  - `LightCurve.from_files([f1, f2, ...])` — Python-side merge with an
    `lcnum` column, suitable for feeding `Pipeline.run(lc)`.

Defaults emit `lcnumvar lcnum` so commands like `-stitch` work without
extra wiring.  PerLC values are accepted on `run_combinelcs()` (one per
group).  `capture_lc=True` returns a `LightCurve` with the combined data
and the `lcnum` column populated.


## `cmd.resample` — closed

The wrapper now accepts ``file_times="list"`` plus ``list_column`` /
``t_column`` int kwargs that map to the CLI ``listcolumn`` / ``tcolumn``
keywords.  The legacy ``file_column`` kwarg still works as an alias for
``t_column`` in path-mode for backward compatibility.

The classic back-resampling pattern (bin to a coarse grid → uniform-grid
resample → ARIMA fit → resample back to the original LC times via the
list-form) can now be expressed without ``subprocess`` (see the
``resample`` example block in the docs).
