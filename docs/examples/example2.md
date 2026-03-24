# Example 2: Lomb-Scargle Period Search and Harmonic Fit

Performing a Lomb-Scargle period search on a light curve and fitting a harmonic series to the
light curve.

## Command

```bash
./vartools -i EXAMPLES/2 -LS 1.0 2.0 0.01 1 0 -Killharm ls 0 0 1 EXAMPLES/OUTDIR1 -oneline
```

## Output

```
Name                                 = EXAMPLES/2
LS_Period_1_0                        =     1.23588006
Log10_LS_Prob_1_0                    = -707.23302
LS_SNR_1_0                           =   18.61023
Killharm_Mean_Mag_1                  =  10.12178
Killharm_Period_1_1                  =     1.23588006
Killharm_Per1_Fundamental_Sincoeff_1 =   0.02004
Killharm_Per1_Fundamental_Coscoeff_1 =  -0.04580
Killharm_Per1_Amplitude_1            =   0.09998
```

## Explanation

Run the Lomb-Scargle period search on the light curve `EXAMPLES/2`. This light curve has had a
sine-curve with a period of 1.2354 days and an amplitude of 0.05 mag injected into it. After
using L-S to find the period, fit a harmonic series with no harmonics or sub-harmonics (i.e.,
just a sine-curve at the period found by L-S) to the light curve. Output the light curve with
the best-fit model in the third column to the file `EXAMPLES/OUTDIR1/2.killharm.model`.

In this case L-S finds the period of 1.23588006 days, and the `-Killharm` routine fits the
function:

```
0.02004*sin((JD - JD0)*2*pi/1.23588006)
  - 0.0458*cos((JD - JD0)*2*pi/1.23588006)
  + 10.12178
```

to the light curve. The killharm routine returns a peak-to-peak amplitude of 0.09998 mag.

This example also illustrates the use of the `-oneline` option, which causes the output to be
written in the format seen above (one line per statistic). This can be useful when processing
single light curves to get results in a more readable format.

## Python Equivalent

```python
import pyvartools as vt
from pyvartools import commands as cmd

lc = vt.LightCurve.from_file("EXAMPLES/2")

result = vt.Pipeline([
    cmd.LS(1.0, 2.0, 0.01, Npeaks=1, output_periodogram=False),
    cmd.Killharm(period_source="ls", nharmonics=0, nsubharmonics=0,
                 correct=True, outdir="EXAMPLES/OUTDIR1"),
]).run(lc)

print(result.stats)
```
