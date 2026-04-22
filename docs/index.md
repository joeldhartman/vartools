# VARTOOLS

**Light curve analysis tools — command-line C program and Python API.**

VARTOOLS is a collection of tools for analyzing astronomical photometric time series (light curves). It is written in C and runs from the command line, or via the Python package, **pyvartools**, processing one or many light curves through a user-defined pipeline of statistical, filtering, period-finding, and model-fitting operations. The tool is designed to enable efficient batch processing of a large collection of light curves.

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

    # Simplest form — pass a filename directly to the command
    result = vt.LS("EXAMPLES/2", 1.0, 2.0, 0.01)
    print(f"Best period: {result.varobjs.LS.Period_1:.5f} d")
    ```

    Every vartools command is available as `vt.CMD(lc_input, ...)`, where
    `lc_input` can be a filename, a `LightCurve` object, a pandas DataFrame,
    a 2-D numpy array, or a `(t, mag, err)` tuple.  For chaining multiple
    commands use `lc.LS(...).rms()` — see the [Python API docs](python/index.md).

---

## Download

!!! note "Current release: vartools 1.6"

    | Format | Link |
    |--------|------|
    | Source tarball | [vartools-1.6.tar.gz](http://www.astro.princeton.edu/~jhartman/vartools/vartools-1.6.tar.gz) |
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

## AI Usage

This website and recent versions of VARTOOLS (1.6+) have been developed with the assistance of [Claude Code](https://claude.com/claude-code) using the Claude Opus 4.7 (1M context) model. The file [`python/llms.txt`](https://github.com/joeldhartman/vartools/blob/master/python/llms.txt) (a plain-text API overview) included in the distribution can be used to provide context to an LLM to understand the usage of this package.

---

## Contact

Questions, bug reports, and feature requests are welcome.

**Joel D. Hartman** — jhartman AT astro DOT princeton DOT edu
