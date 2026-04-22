# vartools

VARTOOLS is an astronomical time-series analysis program providing a wide range of operations on light curves: period finding (Lomb-Scargle, BLS, AOV, WWZ, DFT-CLEAN, autocorrelation), signal injection and recovery, harmonic fitting, transit modeling (Mandel-Agol, trapezoidal), TFA/SYSREM detrending, microlensing, and more.

## Documentation

Full documentation, tutorials, and examples (including the `pyvartools`
Python API) are hosted at
<http://www.astro.princeton.edu/~jhartman/vartools_new/>.

The legacy single-page reference is at
<https://www.astro.princeton.edu/~jhartman/vartools.html>.

## Building from source

```bash
./configure --prefix=$HOME/.local
make
make install
```

See `README` for the full list of dependencies, `--with-*` configure flags, and build options; see `README.linux`, `README.mac`, or `README.windows` for platform-specific instructions (package-manager commands, CSPICE setup, compiler quirks, WSL recommendations for Windows).

Optional dependencies — each enables additional features and is autodetected when installed in standard locations:
- CFITSIO (FITS file I/O)
- GSL (numerical routines used by `-addnoise`, `-microlens`, `-resample`, `-FFT`, `-IFFT`)
- pthread (parallel processing via `-parallel`)
- CSPICE (BJD/UTC conversion via `-converttime`)
- Python 3 + NumPy (embedded `-python` command)
- R (embedded `-R` command)

## Python API

A Python package `pyvartools` is included in the `python/` subdirectory. It wraps the vartools binary (and the in-process `libvartoolspipeline` library installed alongside it) as a Python API, returning results as pandas objects and `LightCurve` wrappers with no manual parsing required.

```python
import pyvartools as vt
from pyvartools import commands as cmd

lc = vt.LightCurve.from_file("EXAMPLES/2")
result = vt.Pipeline([cmd.LS(0.1, 10.0, 1e-3)]).run(lc)
print(result.vars["LS_Period_1_0"])

# Batch processing from disk
batch = vt.Pipeline([cmd.rms()]).run_filelist("EXAMPLES/lc_list", nthreads=4)
print(batch.vars)
```

### Installing the Python package

```bash
cd python
pip install -e .
```

Requires Python 3.8+, numpy, and pandas. See `python/README.md` for full documentation.

## Examples

Example light curves and expected outputs are in the `EXAMPLES/` directory. See the website documentation for a walkthrough of each example.

## License

See `LICENSE`.
