# Example 4: TFA Detrending and BLS Transit Search

De-trending a light curve with TFA and running BLS on the de-trended light curve.

## Command

```bash
./vartools -l EXAMPLES/lc_list_tfa -rms -TFA EXAMPLES/trendlist_tfa EXAMPLES/dates_tfa 25.0 1 0 0 \
  -BLS q 0.01 0.1 0.5 5.0 20000 200 7 2 0 1 EXAMPLES/OUTDIR1 1 fittrap -rms -oneline
```

## Output

```
Name                         = EXAMPLES/3.transit
Mean_Mag_0                   =  10.16727
RMS_0                        =   0.00542
Expected_RMS_0               =   0.00104
Npoints_0                    =  3417
TFA_MeanMag_1                =  10.16714
TFA_RMS_1                    =   0.00471
BLS_Period_1_2               =     2.12353204
BLS_Tc_1_2                   = 53727.296952650177
BLS_SN_1_2                   =  41.96283
BLS_SR_1_2                   =   0.00245
BLS_SDE_1_2                  =   4.91539
BLS_Depth_1_2                =   0.01154
BLS_Qtran_1_2                =   0.03365
BLS_Qingress_1_2             =   0.10122
BLS_OOTmag_1_2               =  10.16694
BLS_i1_1_2                   =   0.98294
BLS_i2_1_2                   =   1.01659
BLS_deltaChi2_1_2            = -25711.70352
BLS_fraconenight_1_2         =   0.43803
BLS_Npointsintransit_1_2     =   157
BLS_Ntransits_1_2            =     4
BLS_Npointsbeforetransit_1_2 =   121
BLS_Npointsaftertransit_1_2  =   138
BLS_Rednoise_1_2             =   0.00137
BLS_Whitenoise_1_2           =   0.00405
BLS_SignaltoPinknoise_1_2    =  15.23677
BLS_Period_2_2               =     1.35195976
BLS_Tc_2_2                   = 53725.941066649029
BLS_SN_2_2                   =  35.29095
BLS_SR_2_2                   =   0.00210
BLS_SDE_2_2                  =   3.76413
BLS_Depth_2_2                =   0.01094
BLS_Qtran_2_2                =   0.05499
BLS_Qingress_2_2             =   0.18241
BLS_OOTmag_2_2               =  10.16694
BLS_i1_2_2                   =   0.53994
BLS_i2_2_2                   =   0.59493
BLS_deltaChi2_2_2            = -18975.44628
BLS_fraconenight_2_2         =   0.53158
BLS_Npointsintransit_2_2     =   155
BLS_Ntransits_2_2            =     5
BLS_Npointsbeforetransit_2_2 =   164
BLS_Npointsaftertransit_2_2  =    89
BLS_Rednoise_2_2             =   0.00154
BLS_Whitenoise_2_2           =   0.00435
BLS_SignaltoPinknoise_2_2    =  14.18251
BLS_Period_invtransit_2      =     1.09096959
BLS_deltaChi2_invtransit_2   = -2374.38718
BLS_MeanMag_2                =  10.16739
Mean_Mag_3                   =  10.16664
RMS_3                        =   0.00405
Expected_RMS_3               =   0.00104
Npoints_3                    =  3417
```

## Explanation

Run TFA and BLS on the light curves listed in `EXAMPLES/lc_list_tfa` (the only light curve in
that list is `EXAMPLES/3.transit`, which was created in Example 3). For TFA we use the trendlist
given in `EXAMPLES/trendlist_tfa` (the second and third columns are the x and y positions for
each star). TFA only includes trend stars that are more than 25 pixels from `EXAMPLES/3.transit`
(the coordinates for this star are given in `EXAMPLES/lc_list_tfa`). We pass the corrected light
curve to the subsequent commands and do not output the TFA coefficients or the TFA model.

We then run BLS on the TFA-cleaned light curve. The search uses:

- A fixed q range of 0.01 to 0.1
- A period range of 0.5 to 5.0 days
- 20,000 frequency points with 200 bins per period
- +7 hours added to UT time to get local time (used for the one-night fraction statistic)
- The top 2 BLS peaks reported
- No periodogram output, but the BLS model light curve output to `EXAMPLES/OUTDIR1/3.transit.bls.model`
- The best-fit box-car transit subtracted before passing the light curve to the next command
- The `fittrap` keyword, which fits a trapezoid transit at each BLS peak

The `fittrap` fit adds the `Qingress` and `OOTmag` output parameters. `Qingress` is the fraction
of the transit duration taken up by ingress (0 for a perfect box-shaped transit, 0.5 for a
V-shaped transit). The RMS is computed before TFA and after BLS.

The best period (`BLS_Period_1_2`) is 2.12353204 days with a delta-chi2 of -25711.70352. The
second-best period is 1.35195976 days. The signal-to-pink-noise for the first peak is 15.23677.
There are 157 points in transit, 4 observed transits, and the estimated red and white noises
after subtracting the model transit are 0.00137 and 0.00405 mag, respectively. TFA reduces
the RMS from 0.00542 mag (`RMS_0`) to 0.00471 mag (`TFA_RMS_1`), and subtracting the transit
reduces it further to 0.00405 mag (`RMS_3`).

## Python Equivalent

```python
import pyvartools as vt
from pyvartools import commands as cmd

lc_list = vt.LightCurveList.from_file("EXAMPLES/lc_list_tfa")

result = vt.Pipeline([
    cmd.rms(),
    cmd.TFA(
        trendlist="EXAMPLES/trendlist_tfa",
        dates="EXAMPLES/dates_tfa",
        min_distance=25.0,
        correct=True,
        output_coefficients=False,
        output_model=False,
    ),
    cmd.BLS(
        q_mode="q", qmin=0.01, qmax=0.1,
        pmin=0.5, pmax=5.0,
        nfreq=20000, nbins=200,
        utoffset=7, npeaks=2,
        output_periodogram=False,
        outdir="EXAMPLES/OUTDIR1",
        correct=True,
        fittrap=True,
    ),
    cmd.rms(),
]).run(lc_list)

print(result.stats)
```
