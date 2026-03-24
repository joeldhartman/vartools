# Batch Processing

VARTOOLS is designed for large-scale processing. A single invocation can process thousands of light curves, writing one output line per light curve. These recipes show the most efficient patterns for each use case.

---

## 1. Process a list of files from disk

The most efficient approach for existing light curve files is `run_filelist`, which avoids Python I/O entirely. Pass either a list file or a Python list of paths.

=== "CLI"

    ```bash
    # Process every LC listed in lc_list.txt
    vartools -l lc_list.txt \
        -clip 5.0 1 \
        -LS 0.5 10.0 0.001 1 0 \
        -oneline > results.txt
    ```

    Each line of `lc_list.txt` is a path to one light curve file. The output is written to `results.txt` with one row per light curve.

=== "Python"

    ```python
    import pyvartools as vt
    from pyvartools import commands as cmd

    pipe = vt.Pipeline([
        cmd.clip(5.0),
        cmd.LS(0.5, 10.0, 1e-3),
    ])

    # From an existing list file
    batch = pipe.run_filelist("/data/lc_list.txt")

    # Or from a Python list of paths
    import glob
    paths = sorted(glob.glob("/data/lcs/*.lc"))
    batch = pipe.run_filelist(paths)

    print(batch.stats)         # pd.DataFrame, one row per LC
    print(batch.stats.shape)   # (N_lcs, N_columns)
    ```

---

## 2. Batch with parallel threads

Pass `nthreads` to process multiple light curves simultaneously. The vartools `-parallel N` flag is used internally.

=== "CLI"

    ```bash
    vartools -l lc_list.txt -parallel 8 \
        -clip 5.0 1 \
        -LS 0.5 10.0 0.001 1 0 \
        -oneline > results.txt
    ```

=== "Python"

    ```python
    batch = pipe.run_filelist("/data/lc_list.txt", nthreads=8)
    ```

    For in-memory batches:

    ```python
    lcs = [vt.LightCurve.from_file(p) for p in paths]
    batch = pipe.run_batch(lcs, nthreads=4)
    ```

!!! tip "Choosing nthreads"
    A good starting point is the number of physical CPU cores. Disk I/O can become the bottleneck for very large files, so experiment with values between 4 and the machine's core count.

---

## 3. Get the statistics table and save to CSV

The `batch.stats` attribute is a pandas DataFrame with one row per light curve and one column per output statistic. Save it directly to CSV or feed it into your downstream analysis.

=== "CLI"

    ```bash
    vartools -l lc_list.txt -parallel 8 \
        -clip 5.0 1 \
        -LS 0.5 10.0 0.001 3 0 \
        -oneline > results.txt
    ```

    Parse the output with pandas:

    ```python
    import pandas as pd
    df = pd.read_csv("results.txt", sep=r"\s+", comment="#")
    df.to_csv("results.csv", index=False)
    ```

=== "Python"

    ```python
    pipe = vt.Pipeline([
        cmd.clip(5.0),
        cmd.columnsuffix("ls"),
        cmd.LS(0.5, 10.0, 1e-3, npeaks=3),
        cmd.columnsuffix("rms"),
        cmd.rms(),
    ])
    batch = pipe.run_filelist("/data/lc_list.txt", nthreads=8)

    # Save the full statistics table
    batch.stats.to_csv("results.csv", index=False)

    # Access specific columns
    print(batch.stats[["Name", "LS_Period_1_ls", "Log10_LS_Prob_1_ls", "RMS_rms"]])

    # Find top candidates
    candidates = batch.stats[batch.stats["Log10_LS_Prob_1_ls"] < -5]
    print(f"{len(candidates)} significant detections")
    ```

---

## 4. Capture modified light curves from a batch run

Set `capture_lc=True` to retrieve the (possibly detrended or clipped) light curves alongside the statistics. The result is available as `batch.lcs`, a list in the same order as the input.

=== "Python"

    ```python
    pipe = vt.Pipeline([
        cmd.clip(5.0),
        cmd.TFA("trendlist.txt", "dates.txt", 0),
    ])

    lcs = [vt.LightCurve.from_file(p) for p in paths]
    batch = pipe.run_batch(lcs, nthreads=4, capture_lc=True)

    for lc, detrended in zip(lcs, batch.lcs):
        if detrended is not None:
            print(f"{lc.name}: {len(detrended)} points after TFA")
    ```

    With `run_filelist`:

    ```python
    batch = pipe.run_filelist("/data/lc_list.txt", capture_lc=True)
    detrended_lcs = [lc for lc in batch.lcs if lc is not None]
    ```

    The detrended light curves can be written back to disk:

    ```python
    import os
    os.makedirs("detrended", exist_ok=True)
    for lc in detrended_lcs:
        lc.to_file(f"detrended/{lc.name}.lc")
    ```

---

## 5. Error handling in batch

By default, a vartools failure raises `RunError`. Set `raise_on_error=False` to instead store the exception in `batch.error` and continue your script.

=== "Python"

    ```python
    from pyvartools.results import RunError

    pipe = vt.Pipeline([cmd.LS(0.5, 10.0, 1e-3)])

    batch = pipe.run_filelist("/data/lc_list.txt", raise_on_error=False)

    if batch.ok:
        print(f"Processed {len(batch.stats)} light curves successfully")
        batch.stats.to_csv("results.csv", index=False)
    else:
        print(f"Batch failed: {batch.error}")
        # batch.stats is an empty DataFrame
    ```

    For more granular error handling, run each light curve individually:

    ```python
    results = []
    errors  = []

    for path in paths:
        try:
            r = pipe.run_file(path, timeout=60)
            results.append(r.stats)
        except RunError as e:
            errors.append((path, str(e)))

    import pandas as pd
    df = pd.DataFrame(results)
    df.to_csv("results.csv", index=False)

    if errors:
        print(f"\n{len(errors)} failures:")
        for path, msg in errors:
            print(f"  {path}: {msg}")
    ```

---

## 6. CLI parallel processing

vartools parallelises at the light-curve level: each thread processes one light curve at a time. Output lines are written in the order jobs complete (not necessarily input order), so always use the `Name` column to identify results.

=== "CLI"

    ```bash
    # 16 parallel threads; redirect both stdout and stderr
    vartools -l lc_list.txt -parallel 16 \
        -clip 5.0 1 \
        -LS 0.5 10.0 0.001 1 0 \
        -BLS r 0.01 0.1 1.0 20.0 optimal 1.0 200 0 1 0 0 1 0 \
        -oneline > results.txt 2>errors.log
    ```

    The `Name` column in `results.txt` gives the light curve identifier for each row.

---

## 7. Survey-scale processing

For a large survey with hundreds of thousands of light curves, use `run_filelist` with a pre-built list file. vartools processes the file in streaming fashion, so memory usage is bounded.

=== "CLI"

    ```bash
    # Generate the list file
    ls /survey/lcs/*.lc > survey_lcs.txt

    # Run the pipeline
    vartools -l survey_lcs.txt -parallel 32 \
        -clip 5.0 1 \
        -LS 0.5 20.0 0.001 1 0 \
        -BLS r 0.01 0.1 1.0 20.0 optimal 1.0 200 0 1 0 0 1 0 \
        -rms \
        -oneline > survey_results.txt
    ```

=== "Python"

    ```python
    pipe = vt.Pipeline([
        cmd.clip(5.0),
        cmd.columnsuffix("ls"),
        cmd.LS(0.5, 20.0, 1e-3, npeaks=1),
        cmd.columnsuffix("bls"),
        cmd.BLS(1.0, 20.0, rmin=0.01, rmax=0.1, nbins=200),
        cmd.columnsuffix("rms"),
        cmd.rms(),
    ])

    batch = pipe.run_filelist(
        "/survey/survey_lcs.txt",
        nthreads=32,
        timeout=3600,   # 1-hour timeout for the full batch
    )

    batch.stats.to_csv("survey_results.csv", index=False)
    print(f"Processed {len(batch.stats)} light curves")

    # Flag candidates for follow-up
    transit_cands = batch.stats[
        (batch.stats["BLS_SDE_1_bls"] > 7) &
        (batch.stats["BLS_Period_1_bls"] > 1.0)
    ]
    transit_cands.to_csv("transit_candidates.csv", index=False)
    print(f"{len(transit_cands)} transit candidates")
    ```

!!! tip "Splitting very large surveys"
    For lists with millions of light curves, split the list file into chunks and run each chunk separately, then concatenate the result CSV files. This also makes it easy to reprocess failed chunks.
