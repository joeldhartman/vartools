# pyvartools — examples notebooks

This directory holds runnable Jupyter notebooks that walk through the worked
examples from the [VARTOOLS docs site](https://hatpipe-doc.example/examples/).

## Files

| Notebook | Covers |
|----------|--------|
| `pyvartools_tour.ipynb` | All seven examples from the docs site (polynomial detrending, period finding, transit injection / search / fit, RR Lyrae recovery, batch variability, Kepler FITS), including the matplotlib plotting steps used to produce the docs-site figures. |

## Running the tour notebook

```bash
jupyter notebook python/examples/pyvartools_tour.ipynb
```

The notebook starts with a **setup cell** in which you must set the
`VARTOOLS_ROOT` variable to the absolute path of your vartools source root
(the directory that contains the `EXAMPLES/` subtree). The setup cell will
raise a clear error until you do, and the rest of the notebook will fail
because the example code uses relative paths like `EXAMPLES/2`.

After the setup cell, run cells top to bottom. Each example section below
the setup cell is independent — you can pick a single example and run only
its cells.

### Required packages

* `pyvartools` itself (install per `python/README.md`).
* `matplotlib` for the plotting cells.
* A built `vartools` binary on `PATH` or pointed to via `$VARTOOLS_BINARY`.

### What's not in the notebook

The transit-search-and-fit example on the docs site has two extra variants
that depend on [CETRA](https://github.com/leigh2/cetra), which requires a
CUDA-capable GPU and the CUDA toolkit. Those variants are intentionally
omitted from the notebook so it stays runnable on any machine with a CPU.
See the docs page if you want to try them.

## Updating the notebook

The notebook is generated from `docs/examples/*.md` on the docs-site repo.
If you change an example's Python code on the docs site, regenerate with the
build script (kept under `tools/` on the docs side, or rerun by hand). The
notebook is intended to track the docs — edit the docs first, then rebuild.
