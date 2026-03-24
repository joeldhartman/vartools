# Getting Started

This page walks through the most common VARTOOLS workflows, showing both the CLI and Python (`pyvartools`) approaches side by side. All examples assume you are working in the VARTOOLS source directory and that the `EXAMPLES/` directory is present.

!!! tip "Installation"
    If you have not yet installed VARTOOLS, see the [Installation](install.md) page before continuing.

---

## 1. Your First VARTOOLS Run

The simplest useful operation is computing the root-mean-square (RMS) scatter of a light curve.

=== "CLI"

    ```bash
    vartools -i EXAMPLES/2 -rms -oneline
    # Output:
    # Name        = EXAMPLES/2
    # Mean_Mag_0  = 10.12178
    # RMS_0       =  0.05012
    # ...
    ```

    - `-i EXAMPLES/2` — read a single light curve from the file `EXAMPLES/2`.
    - `-rms` — compute the RMS and mean magnitude.
    - `-oneline` — print results as a single space-separated line (suitable for piping into scripts or tables).

=== "Python"

    ```python
    import pyvartools as vt
    from pyvartools import commands as cmd

    lc = vt.LightCurve.from_file("EXAMPLES/2")
    result = vt.Pipeline([cmd.rms()]).run(lc)
    print(result.stats["RMS_0"])
    ```

    `result.stats` is a dictionary mapping output column names to their values. Individual statistics are accessible by key.

!!! note "Light curve file format"
    By default, VARTOOLS expects plain-text files with columns: `BJD  magnitude  magnitude_error`. Comment lines beginning with `#` are ignored. See the CLI reference for other supported formats and FITS input.

---

## 2. Chaining Commands

Commands are processed in sequence, so you can combine filtering with analysis in a single call. Here, sigma clipping is applied before the Lomb-Scargle period search.

=== "CLI"

    ```bash
    vartools -i EXAMPLES/2 -clip 5.0 1 -LS 0.5 5.0 0.001 1 0 -oneline
    ```

    - `-clip 5.0 1` — iteratively remove points more than 5σ from the mean; repeat until convergence (`1`).
    - `-LS 0.5 5.0 0.001 1 0` — Lomb-Scargle search from 0.5 to 5.0 days with a step of 0.001 days; return the top 1 peak; do not subtract the best-fit sinusoid.

=== "Python"

    ```python
    result = vt.Pipeline([
        cmd.clip(sigclip=5.0),
        cmd.LS(0.5, 5.0, 1e-3),
    ]).run(lc)
    best_period = float(result.stats["LS_Period_1_1"])
    ```

    Commands in the list are executed left to right. Each command receives the (possibly modified) light curve produced by the previous one.

!!! tip "Output column index"
    The trailing `_1` in `LS_Period_1_1` is the **command index** — this LS search is the second command (0-based index 1) in the pipeline. See [Understanding the Output](#4-understanding-the-output) below for full details.

---

## 3. Batch Processing

VARTOOLS is designed to process large sets of light curves efficiently. Pass a list file with `-l` on the CLI, or use `run_batch` / `run_filelist` in Python.

=== "CLI"

    ```bash
    vartools -l EXAMPLES/lc_list -rms
    ```

    `EXAMPLES/lc_list` is a plain-text file containing one light curve path per line. The output contains one line per light curve, with the `Name` column identifying each file.

=== "Python"

    ```python
    import glob

    lcs = [vt.LightCurve.from_file(f) for f in glob.glob("EXAMPLES/[0-9]*")]
    batch = vt.Pipeline([cmd.rms()]).run_batch(lcs)
    print(batch.stats)

    # Or process files directly from disk (no Python I/O):
    batch = vt.Pipeline([cmd.rms()]).run_filelist("EXAMPLES/lc_list")
    ```

    `run_filelist` passes the file list directly to the vartools binary, which is faster for large sets because Python does not load each light curve into memory.

!!! tip "Parallelism"
    For very large batches on the CLI, use `-header` once followed by multiple parallel `vartools` invocations (one per subset of the list), then concatenate the outputs. The Python `run_filelist` method accepts a `nproc` argument for multiprocessing.

---

## 4. Understanding the Output

### Column naming convention

Every output statistic produced by VARTOOLS follows the pattern:

```
CommandName_StatName_PeakNum_CmdIndex
```

| Field | Meaning |
|-------|---------|
| `CommandName` | The command that produced the statistic (e.g. `LS`, `RMS`, `BLS`) |
| `StatName` | The specific quantity (e.g. `Period`, `Power`, `SDE`) |
| `PeakNum` | For commands that return multiple peaks, the 1-based peak number |
| `CmdIndex` | The 0-based index of this command in the pipeline |

**Examples:**

| Column | Meaning |
|--------|---------|
| `RMS_0` | RMS from the first `-rms` command (index 0) |
| `LS_Period_1_0` | Period of the top peak from the first `-LS` command |
| `LS_Period_1_1` | Period of the top peak from the second command (index 1) |
| `LS_Period_2_0` | Period of the second-best peak from the first `-LS` command |

### Parsing `-oneline` output

With `-oneline`, VARTOOLS prints a single whitespace-separated line whose first field is always the light curve name (the `Name` column). A header line is printed when the `-header` flag is also supplied:

```bash
vartools -l EXAMPLES/lc_list -rms -header -oneline
```

This produces output suitable for reading directly into `numpy.loadtxt`, `pandas.read_csv` (with `sep=r'\s+'`), or similar tools.

---

## 5. Key Concepts

- **Pipeline model.** VARTOOLS reads each light curve, applies commands left to right, and emits statistics. Commands can modify the light curve in place (e.g. clipping, detrending) or simply compute statistics without changing it.

- **`-oneline` flag.** Without this flag, VARTOOLS prints multi-line, human-readable output for each light curve. With `-oneline`, output is one compact line per light curve — far easier to post-process programmatically.

- **The `Name` column.** The first output column is always the input file path (or the identifier field from a list file). It serves as the primary key when joining VARTOOLS output with other catalogs.

- **Output suffixes.** The `_CmdIndex` suffix makes it unambiguous which command produced each statistic when the same command type appears more than once in a pipeline (e.g. two separate `-LS` calls with different period ranges).

- **Stateless commands.** Each light curve is processed independently. There is no global state carried between light curves in a batch run.
