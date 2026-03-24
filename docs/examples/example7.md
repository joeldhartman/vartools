# Example 7: Processing a Kepler FITS Light Curve

Processing a Kepler FITS light curve.

## Command

```bash
./vartools -i kplr001429092-2009166043257_llc.fits \
           -inputlcformat t:1,pdcsap_flux:8,pdcsap_flux_err:9 \
           -changevariable mag pdcsap_flux \
           -changevariable err pdcsap_flux_err \
           -fluxtomag 25.0 0 \
           -rms \
           -LS 0.1 30. 0.1 4 0 \
           -o tmp.lc \
           -oneline
```

## Output

```
Name              = kplr001429092-2009166043257_llc.fits
Mean_Mag_3        =  14.48415
RMS_3             =   0.00346
Expected_RMS_3    =   0.00037
Npoints_3         =  1624
LS_Period_1_4     =     8.16381396
Log10_LS_Prob_1_4 = -185.09080
LS_SNR_1_4        = 7939.08565
LS_Period_2_4     =     4.46288496
Log10_LS_Prob_2_4 =  -65.54306
LS_SNR_2_4        = 2887.81992
LS_Period_3_4     =     3.84731462
Log10_LS_Prob_3_4 =   -3.96317
LS_SNR_3_4        =  285.87713
LS_Period_4_4     =     2.88548597
Log10_LS_Prob_4_4 =   -0.01046
LS_SNR_4_4        =   94.27330
```

## Explanation

This example processes the Q1 Kepler public light curve for KIC 1429092, which can be obtained
from the Kepler MAST archive. The light curve is stored in a binary FITS table format.

The `-inputlcformat` command specifies which table columns to use: column 1 is BJD - 2454833,
column 8 is the PDC-corrected simple aperture photometry (PDCSAP) flux, and column 9 is the
corresponding uncertainty. The `-changevariable` commands associate the special variables `mag`
and `err` with `pdcsap_flux` and `pdcsap_flux_err`, respectively. In this particular example
we could have read columns 8 and 9 directly into `mag` and `err`, but explicitly naming the
columns makes it easier to track which quantities are being processed when there are multiple
photometric options in the file.

The `-fluxtomag` command converts the PDCSAP flux to magnitudes. The RMS of the resulting light
curve is then computed, and the Lomb-Scargle periodogram is applied over the range 0.1–30 days,
reporting the 4 highest peaks. The BJD, magnitude, and magnitude uncertainty are written to the
ASCII file `tmp.lc`.

## Python Equivalent

```python
import pyvartools as vt
from pyvartools import commands as cmd

lc = vt.LightCurve.from_file(
    "kplr001429092-2009166043257_llc.fits",
    t_col="TIME", mag_col="PDCSAP_FLUX", err_col="PDCSAP_FLUX_ERR"
)
result = vt.Pipeline([
    cmd.fluxtomag(25.0),
    cmd.rms(),
    cmd.LS(0.1, 30.0, 0.1, Npeaks=4),
]).run(lc)
print(result.stats)
```
