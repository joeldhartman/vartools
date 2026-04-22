# Installation

VARTOOLS has two components:

- **vartools** — the compiled C binary, required for all analyses.
- **pyvartools** — the Python package that wraps the binary; requires the binary to be installed first.

Full per-platform detail is in the [`README.linux`](#readme-linux), [`README.mac`](#readme-mac), and [`README.windows`](#readme-windows) files in the source distribution. The three recipes below give a complete end-to-end walkthrough for each platform (system prerequisites → CSPICE → vartools → pyvartools); the [general notes](#general-notes-and-configuration) further down expand on individual points.

---

## Download the source

| Format | Link |
|--------|------|
| Source tarball | [vartools-1.6.tar.gz](http://www.astro.princeton.edu/~jhartman/vartools/vartools-1.6.tar.gz) |
| GitHub | [github.com/joeldhartman/vartools](https://github.com/joeldhartman/vartools) |

---

## Quick-start recipes

Pick your platform. Each recipe installs the C binary, all optional
dependencies (CFITSIO, GSL, CSPICE, Python+NumPy, R), and the pyvartools
Python package.

=== "Linux"

    **Recommended — VARTOOLS is developed primarily on Linux.**

    **1. Install system prerequisites** (required + packaged optional
    dependencies).

    Debian / Ubuntu:

    ```bash
    sudo apt-get install -y \
        build-essential gfortran libtool csh wget \
        libcfitsio-dev libgsl-dev \
        python3-dev python3-numpy python3-pip \
        r-base-dev
    ```

    Fedora / RHEL / Rocky / AlmaLinux:

    ```bash
    sudo dnf install -y \
        gcc gcc-gfortran make libtool tcsh wget \
        cfitsio-devel gsl-devel \
        python3-devel python3-numpy python3-pip \
        R-core-devel
    ```

    **2. Build CSPICE from source** (not packaged; used by the
    `-converttime` command):

    ```bash
    mkdir -p ~/src && cd ~/src
    wget https://naif.jpl.nasa.gov/pub/naif/toolkit//C/PC_Linux_GCC_64bit/packages/cspice.tar.Z
    gunzip cspice.tar.Z
    tar -xf cspice.tar
    cd cspice
    csh makeall.csh
    ```

    This produces `~/src/cspice/lib/cspice.a` and headers in
    `~/src/cspice/include/`. You'll also want the ancillary kernel
    files — see [CSPICE kernel files](#cspice-kernel-files) below.

    **3. Unpack and build VARTOOLS:**

    ```bash
    cd ~/src
    wget http://www.astro.princeton.edu/~jhartman/vartools/vartools-1.6.tar.gz
    tar xzf vartools-1.6.tar.gz
    cd vartools-1.6

    ./configure \
        --prefix=$HOME/.local \
        --with-cspice=$HOME/src/cspice \
        --with-RHOME=$(R RHOME) \
        FC=gfortran F77=gfortran
    make
    make install
    ```

    Ensure `$HOME/.local/bin` is on your `PATH`:

    ```bash
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
    ```

    **4. Install pyvartools** (the Python wrapper):

    ```bash
    cd python
    pip install -e .
    ```

    **5. Verify:**

    ```bash
    vartools -help
    python -c "import pyvartools as vt; print(vt.LightCurve.from_file('EXAMPLES/2'))"
    ```

    Full details, including Arch Linux and troubleshooting, are in
    [`README.linux`](#readme-linux) in the source tarball.

=== "macOS"

    Supported on Apple Silicon and Intel Macs via [Homebrew](https://brew.sh).

    **1. Install Xcode Command Line Tools:**

    ```bash
    xcode-select --install
    ```

    **2. Install system prerequisites.** Apple's default `libtool` is
    **not** GNU libtool — Homebrew installs the correct one as
    `glibtool`; this only matters if you run `./autogen.sh`, not for
    normal tarball builds. `pthread` is built into macOS — no separate
    install needed.

    ```bash
    brew install \
        gcc               `# brings gfortran` \
        libtool wget \
        cfitsio gsl \
        python@3 numpy r
    ```

    **3. Build CSPICE from source** (not available from Homebrew).
    Pick the matching distribution for your CPU:

    For Apple Silicon (M1/M2/M3/M4):

    ```bash
    mkdir -p ~/src && cd ~/src
    curl -O https://naif.jpl.nasa.gov/pub/naif/toolkit//C/MacM1_OSX_clang_64bit/packages/cspice.tar.Z
    gunzip cspice.tar.Z
    tar -xf cspice.tar
    cd cspice
    sh makeall.csh
    ```

    For Intel Macs, download from
    `https://naif.jpl.nasa.gov/pub/naif/toolkit//C/PC_MacOSX_GCC_64bit/packages/cspice.tar.Z`
    and run the same commands. You'll also want the ancillary kernel
    files — see [CSPICE kernel files](#cspice-kernel-files) below.

    **4. Unpack and build VARTOOLS.** Homebrew installs into
    `/opt/homebrew` on Apple Silicon and `/usr/local` on Intel; point
    `./configure` at it either way:

    ```bash
    cd ~/src
    curl -O http://www.astro.princeton.edu/~jhartman/vartools/vartools-1.6.tar.gz
    tar xzf vartools-1.6.tar.gz
    cd vartools-1.6

    HOMEBREW_PREFIX=$(brew --prefix)
    ./configure \
        --prefix=$HOME/.local \
        --with-cspice=$HOME/src/cspice \
        --with-RHOME=$(R RHOME) \
        CPPFLAGS="-I${HOMEBREW_PREFIX}/include" \
        LDFLAGS="-L${HOMEBREW_PREFIX}/lib" \
        FC=gfortran F77=gfortran
    make
    make install
    ```

    Ensure `$HOME/.local/bin` is on your `PATH`:

    ```bash
    echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc
    source ~/.zshrc
    ```

    **5. Install pyvartools:**

    ```bash
    cd python
    pip install -e .
    ```

    **6. Verify:**

    ```bash
    vartools -help
    python -c "import pyvartools as vt; print(vt.LightCurve.from_file('EXAMPLES/2'))"
    ```

    A known macOS-specific issue with unnamed POSIX semaphores is worked
    around in the source (VARTOOLS 1.51+), so `-parallel` works
    correctly.

    Full details (Fortran-linker issues, CSPICE on Apple Silicon,
    troubleshooting) are in [`README.mac`](#readme-mac) in the source tarball.

=== "Windows"

    !!! caution "Instructions below not tested"
        The instructions below were generated by Claude and have not been
        tested on an actual Windows system.

    **Recommended path: WSL (Windows Subsystem for Linux).** WSL provides
    a full Linux environment in which every VARTOOLS feature works
    unchanged.

    **1. Install WSL** (Windows 10 2004+ or Windows 11). From an
    Administrator PowerShell:

    ```powershell
    wsl --install
    ```

    This installs Ubuntu by default; reboot when prompted.

    **2. Open the "Ubuntu" Start-menu entry and follow the full Linux
    recipe (first tab above) verbatim inside the WSL shell** — all the
    way through CSPICE, VARTOOLS, and pyvartools. Every feature works,
    including `-python`, `-R`, and `-parallel`.

    For best performance, keep your source and light-curve files on
    the WSL filesystem (`~/` etc.), not on Windows-side paths
    (`/mnt/c/...`).

    ---

    **Alternative: MSYS2 / MinGW-w64** (native Windows build). Core
    VARTOOLS compiles and runs, but the embedded `-python` and `-R`
    commands do **not** work — their implementation uses `fork()` and
    Unix-domain sockets, which are not available on native Windows.
    pyvartools itself (the Python wrapper) still works for everything
    except those two commands.

    **1. Install MSYS2** from <https://www.msys2.org/>, then open the
    **MINGW64** shell (not the generic MSYS2 shell):

    ```bash
    pacman -Syu      # run twice to complete updates
    pacman -Syu
    pacman -S --needed \
        mingw-w64-x86_64-gcc \
        mingw-w64-x86_64-gcc-fortran \
        mingw-w64-x86_64-make \
        mingw-w64-x86_64-libtool \
        mingw-w64-x86_64-cfitsio \
        mingw-w64-x86_64-gsl \
        mingw-w64-x86_64-python \
        mingw-w64-x86_64-python-numpy \
        mingw-w64-x86_64-python-pip \
        tcsh wget
    ```

    **2. Build CSPICE** using NAIF's `PC_Cygwin_GCC_64bit` distribution
    (works under MSYS2 / MinGW-w64):

    ```bash
    mkdir -p ~/src && cd ~/src
    curl -O https://naif.jpl.nasa.gov/pub/naif/toolkit//C/PC_Cygwin_GCC_64bit/packages/cspice.tar.Z
    gunzip cspice.tar.Z
    tar -xf cspice.tar
    cd cspice
    sh makeall.csh
    ```

    **3. Build VARTOOLS** — disable the POSIX-only features:

    ```bash
    cd ~/src
    curl -O http://www.astro.princeton.edu/~jhartman/vartools/vartools-1.6.tar.gz
    tar xzf vartools-1.6.tar.gz
    cd vartools-1.6

    ./configure \
        --prefix=$HOME/.local \
        --with-cspice=$HOME/src/cspice \
        --without-python \
        FC=gfortran F77=gfortran
    make
    make install
    ```

    **4. Install pyvartools:**

    ```bash
    cd python
    pip install -e .
    ```

    **5. Verify:**

    ```bash
    $HOME/.local/bin/vartools.exe -help
    python -c "import pyvartools as vt; print(vt.LightCurve.from_file('EXAMPLES/2'))"
    ```

    Visual Studio / MSVC builds are **not** supported.

    Full details in [`README.windows`](#readme-windows) in the source tarball.

---

## General notes and configuration

### Required prerequisites

A working C compiler, `make`, GNU `libtool`, and a Fortran compiler
(`gfortran`) — the last is needed to compile a small number of
extension commands written in Fortran.

### Optional dependencies

Every item below is autodetected by `./configure` and its corresponding
feature is silently disabled if the detection fails. Pass the listed
`--with-*` flag when a library lives outside standard system paths
(`/usr`, `/usr/local`, Homebrew prefix).

| Library | Purpose | Configure flag |
|---------|---------|----------------|
| [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/) | FITS light-curve I/O | `--with-cfitsio=/path` |
| [GSL](https://www.gnu.org/software/gsl/) | `-addnoise`, `-microlens`, `-resample`, `-FFT`, `-IFFT` commands | `--with-gsl` / `--with-gsl=no` |
| pthread | Parallel processing (`-parallel N`) | *(autodetect only)* |
| [CSPICE](https://naif.jpl.nasa.gov/naif/toolkit.html) | BJD/UTC conversion (`-converttime`) | `--with-cspice=/path` |
| Python 3 + NumPy | Embedded `-python` command | `--with-pythonconfig=/path`, `--with-pythonhome=/path` |
| R | Embedded `-R` command | `--with-RHOME=/path` |

!!! note "Python and R are compile-time dependencies"
    Unlike command-line tools that simply shell out, the `-python` and
    `-R` commands **embed** the Python and R interpreters directly into
    the VARTOOLS binary. Detection happens at `./configure` time —
    `configure` probes for `Python.h` + NumPy headers and for `libR.so`
    + the R headers, and the resulting features are compiled in or not.
    Installing Python or R after building VARTOOLS does not enable
    these commands; you must re-run `./configure` and rebuild.

### What `make install` places where

- The `vartools` executable is placed in `$PREFIX/bin`.
- The shared library `libvartoolspipeline.so` (used by `pyvartools` for
  in-process execution — see [below](#library-mode-libvartoolspipelineso))
  is placed in `$PREFIX/lib`.
- User-extension directories, example data, and documentation are
  installed under `$PREFIX/share/vartools/`.

### CSPICE kernel files

In addition to the CSPICE library, the `-converttime` command reads a
set of SPICE kernel files at runtime — a planetary ephemeris, a
leap-second table, and planet orientation data. Download the current
versions and place them in `$HOME/.vartools/` (the default location
VARTOOLS searches):

```bash
mkdir -p ~/.vartools && cd ~/.vartools
wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de432s.bsp
wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls
wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc
```

(On macOS and MSYS2, use `curl -O` in place of `wget`.) Browse the
parent directories periodically for newer versions — the leap-second
file in particular is updated whenever a new leap second is announced.

Set the environment variables below (e.g. in `.bashrc` / `.zshrc`) so
VARTOOLS can find the files without an explicit command-line path:

```bash
export CSPICE_EPHEM_FILE=$HOME/.vartools/de432s.bsp
export CSPICE_LEAPSEC_FILE=$HOME/.vartools/naif0012.tls
export CSPICE_PLANETDATA_FILE=$HOME/.vartools/pck00010.tpc
```

### Custom install prefix

If you do not have write access to `/usr/local`, install into your home
directory (or another location you own) with `--prefix`:

```bash
./configure --prefix=$HOME/.local
make
make install
```

Ensure `$PREFIX/bin` is on your `PATH`.

### Pointing at non-standard library locations

If an optional dependency lives somewhere `./configure` doesn't find,
two mechanisms are available:

1. **Dedicated `--with-*` flags** (see the table above). These are the
   cleanest option when applicable.
2. **Generic `CPPFLAGS` / `LDFLAGS`** for anything else. These are
   passed directly to the compiler and linker:

   ```bash
   ./configure \
       CPPFLAGS="-I$HOME/src/cfitsio -I/opt/local/include" \
       LDFLAGS="-L$HOME/src/cfitsio -L/opt/local/lib"
   ```

### Regenerating the build system

The tarball ships with a pre-generated `configure` script, so
`autoconf` / `automake` are **not** required on the build machine.

If `./configure` fails with errors about m4 macros or missing helper
scripts — which usually only happens when building from a fresh
checkout of the source repository rather than a release tarball —
regenerate with:

```bash
./autogen.sh
```

This wraps `autoreconf --install --force` and requires `autoconf`,
`automake`, and GNU `libtool` to be installed.

### Platform README files

Each platform's quick-start recipe above is a summary of the more
detailed `README.{linux,mac,windows}` files shipped in the source
distribution:

- <a id="readme-linux"></a>**`README.linux`** — Debian/Ubuntu, Fedora/RHEL, Arch; CSPICE from source; troubleshooting for missing `gfortran`, missing `r-base-dev`, `LD_LIBRARY_PATH` issues.
- <a id="readme-mac"></a>**`README.mac`** — Homebrew install flow; Apple Silicon vs. Intel prefix handling; `gfortran` via `brew install gcc`; `libtool` → `glibtool` for `./autogen.sh`; CSPICE for Apple Silicon; troubleshooting for linker errors and cross-architecture mixups.
- <a id="readme-windows"></a>**`README.windows`** — WSL instructions; MSYS2/MinGW-w64 native-build details and caveats (`--without-python` requirement; no MSVC support); interop with Windows-side files.

---

## Python Package (pyvartools) — reference

The platform recipes above already include the pyvartools install
step. This section covers configuration options that may be useful in
non-standard setups.

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
print(result.vars)
```

If this prints a pandas Series of statistics without errors, everything is working correctly.
