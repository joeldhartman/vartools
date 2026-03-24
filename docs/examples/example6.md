# Example 6: Parallel Variability Search

A script to run a battery of variability selection algorithms on a number of light curves in
parallel.

## Script

This example is included in the package as `EXAMPLES/example6.sh`.

```bash
#!/bin/sh
star="$1"

if [ -n "$star" ] ; then

    vartools -i $star \
        -savelc \
        -clip 5.0 1 \
        -savelc \
        -LS 0.1 100. 0.1 3 0 clip 5. 1 \
        -aov 0.1 100. 0.1 0.01 1 0 clip 5. 1 \
        -aov_harm 1 0.1 100. 0.1 0.01 1 0 clip 5. 1 \
        -restorelc 1 \
        -clip 10.0 1 \
        -BLS q 0.01 0.1 0.1 20. 10000 200 7 2 0 0 0 \
        -restorelc 2 \
        -changeerror \
        -autocorrelation 0. 30. 0.1 EXAMPLES/OUTDIR1/ \
        -nobuffer

else

    PEXEC=${HOME}/bin/pexec
    SELF=$0
    LCLIST=EXAMPLES/lc_list

    vartools \
        -savelc \
        -clip 5.0 1 \
        -savelc \
        -LS 0.1 100. 0.1 3 0 clip 5. 1 \
        -aov 0.1 100. 0.1 0.01 1 0 clip 5. 1 \
        -aov_harm 1 0.1 100. 0.1 0.01 1 0 clip 5. 1 \
        -restorelc 1 \
        -clip 10.0 1 \
        -BLS q 0.01 0.1 0.1 20. 10000 200 7 2 0 0 0 \
        -restorelc 2 \
        -changeerror \
        -autocorrelation 0. 30. 0.1 EXAMPLES/OUTDIR1/ \
        -headeronly \
        -nobuffer \
        -numbercolumns

    $PEXEC -f $LCLIST -e star -R -t -z 19 -n 4 -c eval $SELF \$star
fi
```

## Explanation

This script runs the `-LS`, `-aov`, `-aov_harm`, `-BLS`, and `-autocorrelation` commands on a
list of light curves in parallel. Parallelization uses the
[pexec program](http://www.gnu.org/software/pexec/) by Andras Pal.

The `if` statement splits the script into two branches:

**Per-star branch** (when a light curve name is passed as an argument): vartools processes that
single light curve and writes its results. The `-nobuffer` option ensures output is flushed
immediately, which is important when many processes write to the same output stream.

**Initialization branch** (when called with no arguments): vartools is invoked with no input
light curve and with `-headeronly`, which outputs just the column header. After printing the
header, `pexec` is called. The `-f` option tells `pexec` to read light curve names from
`EXAMPLES/lc_list`. Each name is assigned to the environment variable `star`, and `pexec` calls
this script passing `$star` as an argument. The call runs 4 processes in parallel at nice +19.

The processing pipeline for each light curve:

1. Save the initial light curve state (`-savelc`).
2. Apply iterative 5-sigma clipping (`-clip 5.0 1`).
3. Save the 5-sigma clipped state (`-savelc`).
4. Run `-LS`, `-aov`, and `-aov_harm` period searches (each also applies 5-sigma clipping
   internally).
5. Restore the light curve to its pre-5-sigma-clip state (`-restorelc 1`).
6. Apply 10-sigma clipping and run `-BLS` (a less aggressive clip avoids removing the
   eclipses we are searching for).
7. Restore the light curve to its 5-sigma-clipped state (`-restorelc 2`).
8. Replace the photometric errors with the light curve RMS (`-changeerror`).
9. Run `-autocorrelation`, writing the autocorrelation function to `EXAMPLES/OUTDIR1/`
   (e.g., `EXAMPLES/OUTDIR1/2.autocorr`, which is periodic with a first peak at 1.23 days).

## Python Equivalent

```python
import pyvartools as vt
from pyvartools import commands as cmd
from concurrent.futures import ProcessPoolExecutor

lc_names = open("EXAMPLES/lc_list").read().split()

pipeline = vt.Pipeline([
    cmd.savelc(),
    cmd.clip(5.0, niter=1),
    cmd.savelc(),
    cmd.LS(0.1, 100.0, 0.1, Npeaks=3, output_periodogram=False,
           clip=5.0, clip_niter=1),
    cmd.aov(0.1, 100.0, 0.1, 0.01, Npeaks=1, output_periodogram=False,
            clip=5.0, clip_niter=1),
    cmd.aov_harm(1, 0.1, 100.0, 0.1, 0.01, Npeaks=1, output_periodogram=False,
                 clip=5.0, clip_niter=1),
    cmd.restorelc(1),
    cmd.clip(10.0, niter=1),
    cmd.BLS(q_mode="q", qmin=0.01, qmax=0.1, pmin=0.1, pmax=20.0,
            nfreq=10000, nbins=200, utoffset=7, npeaks=2,
            output_periodogram=False, outdir=None),
    cmd.restorelc(2),
    cmd.changeerror(),
    cmd.autocorrelation(0.0, 30.0, 0.1, outdir="EXAMPLES/OUTDIR1/"),
])

def process_one(name):
    lc = vt.LightCurve.from_file(name)
    return pipeline.run(lc).stats

with ProcessPoolExecutor(max_workers=4) as pool:
    results = list(pool.map(process_one, lc_names))

for r in results:
    print(r)
```
