# Installation

VARTOOLS has two independent components that can be installed separately or together:

- **vartools** — the compiled C binary, required for all analyses.
- **pyvartools** — the Python package that wraps the binary; requires the binary to be installed first.

---

## C Program (vartools binary)

### Download the source

| Format | Link |
|--------|------|
| Source tarball | [vartools-1.52.tar.gz](http://www.astro.princeton.edu/~jhartman/vartools/vartools-1.52.tar.gz) |
| GitHub | [github.com/joeldhartman/vartools](https://github.com/joeldhartman/vartools) |

### Build and install

VARTOOLS uses the standard GNU Autotools build system.

```bash
tar xzf vartools-1.52.tar.gz
cd vartools-1.52
./configure --prefix=/usr/local
make
make install
```

After `make install`, the `vartools` executable is placed in `$PREFIX/bin`
(default `/usr/local/bin`) and the shared library `libvartoolspipeline.so`
is placed in `$PREFIX/lib`. Verify the installation:

```bash
vartools -help
```

!!! tip "Custom install prefix"
    If you do not have write access to `/usr/local`, install into your home directory:
    ```bash
    ./configure --prefix=$HOME/.local
    make
    make install
    ```
    Ensure `$HOME/.local/bin` is on your `PATH`.

### Platform notes

VARTOOLS has been tested on **Linux**, **macOS**, and **Windows (MinGW)**.

=== "macOS"

    Detailed step-by-step notes for building on macOS are available in the source distribution:

    [`MacInstallNotes_2019.0415.txt`](http://www.astro.princeton.edu/~jhartman/vartools/MacInstallNotes_2019.0415.txt)

    The main things to watch for are using a GNU-compatible `make` (e.g. from Homebrew: `brew install make`) and ensuring that optional library headers are on the compiler search path.

=== "Windows (MinGW)"

    Build instructions for Windows using MinGW are available in the source distribution:

    [`VartoolsWindowsInstallation-1.txt`](http://www.astro.princeton.edu/~jhartman/vartools/VartoolsWindowsInstallation-1.txt)

    Using the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/) and following the Linux instructions is often the easier path on modern Windows systems.

### Optional libraries

VARTOOLS has several optional dependencies. The `./configure` script will detect them automatically if they are installed in standard locations; otherwise pass the appropriate `--with-*` flag.

| Library | Purpose | Configure flag |
|---------|---------|----------------|
| [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/) | Read/write FITS files | `--with-cfitsio=/path` |
| [CSPICE](https://naif.jpl.nasa.gov/naif/toolkit.html) | Barycentric Julian Date (BJD) conversion | `--with-cspice=/path` |
| [GSL](https://www.gnu.org/software/gsl/) | Additional numerical routines | `--with-gsl=/path` |

!!! note "Python and R integration"
    The `-python` and `-R` inline-scripting commands require working `python3` and `Rscript` executables on your `PATH`. No special compile-time flags are needed — VARTOOLS invokes them at runtime.

---

## Python Package (pyvartools)

### Requirements

| Package | Version | Notes |
|---------|---------|-------|
| Python | 3.8+ | |
| numpy | any recent | |
| pandas | any recent | |
| astropy | optional | Required for FITS I/O via pyvartools |

### Install

The package is not yet on PyPI. Install directly from the source tree:

```bash
pip install -e /path/to/vartools/source/python
```

Alternatively, use the provided conda environment file (recommended if you want a fully reproducible environment):

```bash
conda env create -f environment.yml
conda activate pyvartools
```

### Library mode (`libvartoolspipeline.so`)

After `make install`, the shared library `libvartoolspipeline.so` is
installed alongside the binary. When pyvartools can find this library it
runs vartools **in-process** (no subprocess spawned per call), giving a
roughly 20× speedup for single-LC and batch runs. If the library is not
found, pyvartools falls back transparently to the subprocess path.

See the [Pipeline — Performance: library mode](python/pipeline.md#performance-library-mode)
section for details on controlling this behaviour.

### Configuring the vartools binary path

`pyvartools` locates the `vartools` binary by searching your `PATH` automatically. If you installed the binary in a non-standard location, you can override this in two ways:

=== "Python"

    ```python
    import pyvartools as vt
    vt.set_binary("/path/to/vartools")
    ```

    Call `set_binary` before creating any `Pipeline` or `LightCurve` objects.

=== "Environment variable"

    ```bash
    export VARTOOLS_BINARY=/path/to/vartools
    ```

    Set this in your shell profile (`.bashrc`, `.zshrc`, etc.) to make it permanent.

### Configuring the library path

`pyvartools` searches the vartools binary's RPATH and standard system paths
for `libvartoolspipeline.so`. If the library is installed in a non-standard
location, override the search:

=== "Python"

    ```python
    import pyvartools as vt
    vt.set_library("/path/to/libvartoolspipeline.so")
    ```

=== "Environment variable"

    ```bash
    export VARTOOLS_LIBRARY=/path/to/libvartoolspipeline.so
    ```

!!! note
    Even without `libvartoolspipeline.so`, all pyvartools functionality works
    via the subprocess fallback. The library only affects performance, not
    correctness.

---

## Verifying the full installation

Once both components are installed, run a quick end-to-end check:

```python
import pyvartools as vt
from pyvartools import commands as cmd

lc = vt.LightCurve.from_file("EXAMPLES/2")
result = vt.Pipeline([cmd.rms()]).run(lc)
print(result.stats)
```

If this prints a dictionary of statistics without errors, everything is working correctly.
