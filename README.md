# vartools

VARTOOLS is an astronomical time series analysis program providing a wide range of operations on light curves: period finding (Lomb-Scargle, BLS, AOV, WWZ, DFT-CLEAN, autocorrelation), signal injection and recovery, harmonic fitting, transit modeling (Mandel-Agol, trapezoidal), TFA/SYSREM detrending, microlensing, and more.

## Documentation

Full documentation and examples are at https://www.astro.princeton.edu/~jhartman/vartools.html

## Building from source

```bash
./configure --bindir=/usr/local/bin
make
make install
```

See `INSTALL` for detailed build options including optional dependencies (CFITSIO, R, Python, CSPICE).

## Python API

A Python package `pyvartools` is included in the `python/` subdirectory. It wraps the vartools binary (and optionally the in-process `libvartoolspipeline` library) as a Python API, returning results as pandas DataFrames and `LightCurve` objects with no manual parsing required.

```python
import pyvartools as vt
from pyvartools import commands as cmd

lc = vt.LightCurve.from_file("EXAMPLES/2")
result = vt.Pipeline([cmd.LS(0.1, 10.0, 1e-3)]).run(lc)
print(result.stats["LS_Period_1_0"])

# Batch processing from disk
batch = vt.Pipeline([cmd.rms()]).run_filelist("EXAMPLES/lc_list", nthreads=4)
print(batch.stats)
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
