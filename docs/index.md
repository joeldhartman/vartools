# VARTOOLS

**Light curve analysis for astronomers — command-line C program and Python API.**

VARTOOLS is a fast, flexible program for analyzing astronomical photometric time series (light curves). It is written in C and runs from the command line, processing one or many light curves through a user-defined pipeline of statistical, filtering, period-finding, and model-fitting operations. A Python package, **pyvartools**, wraps the same binary so that the full command set is available from within Python scripts and notebooks, with results returned as structured objects instead of text output.

---

## Key Features

<div class="grid cards" markdown>

-   **Period Finding**

    ---

    Lomb-Scargle, BLS (box-fitting least squares), phase-folding, and more — all optimised for speed on large light curve sets.

-   **Model Fitting**

    ---

    Harmonic series, trapezoid transit, Mandel-Agol transit, eclipsing-binary, and user-supplied model fitting with robust statistics.

-   **Filtering & Detrending**

    ---

    Sigma clipping, TFA (Trend Filtering Algorithm), EPD, SysRem, polynomial detrending, and moving-average smoothing.

-   **Statistics**

    ---

    RMS, median absolute deviation, Stetson indices, percentile ratios, autocorrelation function, and many more variability metrics.

-   **Batch Processing**

    ---

    Process thousands of light curves in a single call using `-l`; output is one row per light curve for easy downstream analysis.

-   **Python API**

    ---

    `pyvartools` exposes every command as a Python object, returning results as pandas DataFrames and structured `Result` objects.

-   **FITS Support**

    ---

    Read and write FITS binary tables when compiled with cfitsio. Native support for Kepler, TESS, and HATNet formats.

-   **Extensible**

    ---

    Add custom commands in C or call external Python/R scripts inline with `-python` and `-R`.

</div>

---

## Quick Start

=== "CLI"

    ```bash
    vartools -i EXAMPLES/2 -LS 1.0 2.0 0.01 1 0 -oneline
    ```

    This command reads the light curve file `EXAMPLES/2`, searches for periods between 1.0 and 2.0 days with a step size of 0.01, and prints a single summary line to standard output.

=== "Python (pyvartools)"

    ```python
    import pyvartools as vt
    from pyvartools import commands as cmd

    lc = vt.LightCurve.from_file("EXAMPLES/2")
    result = vt.Pipeline([cmd.LS(1.0, 2.0, 0.01)]).run(lc)
    print(f"Best period: {float(result.stats['LS_Period_1_0']):.5f} d")
    ```

    The same Lomb-Scargle search is expressed as a `Pipeline` of `cmd.LS` objects. The `result.stats` dictionary maps output column names to their values.

---

## Download

!!! note "Current release: vartools 1.52"

    | Format | Link |
    |--------|------|
    | Source tarball | [vartools-1.52.tar.gz](http://www.astro.princeton.edu/~jhartman/vartools/vartools-1.52.tar.gz) |
    | GitHub | [github.com/joeldhartman/vartools](https://github.com/joeldhartman/vartools) |

See the [Installation](install.md) page for build instructions, optional dependencies, and platform-specific notes.

---

## Citation

If you use VARTOOLS in published research, please cite:

> Hartman, J. D. & Bakos, G. Á. 2016, *Astronomy and Computing*, **17**, 1

```bibtex
@article{Hartman2016,
  author  = {Hartman, Joel D. and Bakos, G{\'a}sp{\'a}r {\'A}.},
  title   = {{VARTOOLS}: A program for characterizing and searching for variable stars},
  journal = {Astronomy and Computing},
  year    = {2016},
  volume  = {17},
  pages   = {1},
  doi     = {10.1016/j.ascom.2016.05.006}
}
```

---

## Contact

Questions, bug reports, and feature requests are welcome.

**Joel D. Hartman** — jhartman AT astro DOT princeton DOT edu
