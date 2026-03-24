# Example 1: Fitting a Quadratic Polynomial in JD

Fitting a quadratic polynomial in JD to a list of light curves.

## Command

```bash
./vartools -l EXAMPLES/lc_list -rms -decorr 1 1 1 0 1 1 2 0 -rms -tab
```

## Output

```
Name Mean_Mag_0 RMS_0 Expected_RMS_0 Npoints_0 Decorr_constant_term_1 Decorr_constant_term_err_1 \
LCColumn_1_coeff_1_1 LCColumn_1_coeff_err_1_1 LCColumn_1_coeff_2_1 LCColumn_1_coeff_err_2_1 Decorr_chi2_1 \
Mean_Mag_2 RMS_2 Expected_RMS_2  Npoints_2
---- ---------- ----- -------------- --------- ---------------------- -------------------------- \
-------------------- ------------------------ -------------------- ------------------------ ------------- \
---------- ----- --------------  ---------
EXAMPLES/1 10.24745 0.15944 0.00101 3122 10.0830375984825 0.0000325849746 0.0097933162509 0.0000059117875 \
0.0002554062775 0.0000001956124 6.68601 10.24728 0.00211 0.00101 3122
...
```

## Explanation

Fit a quadratic polynomial to the light curves given in the file `EXAMPLES/lc_list`. To do this
we use no global terms and 1 lc term. For the lc term we use the first column in the light curve
(the JD) and fit a second order polynomial in this term to each light curve. We fit for the
zero-point term and correct the light curve (so that commands after `-decorr` will receive light
curves with the best-fit quadratic polynomial removed; note that when the light curve is corrected
the mean is kept constant). We also subtract the first term in the signal that we are decorrelating
against (use JD - JD0 rather than JD, since JD*JD runs into round-off problems whereas
(JD - JD0)*(JD - JD0) does not), but we do not output the corrected light curves to disk. The rms
is determined before and after the fit, and the `-tab` option causes the output table to be in
starbase format with headers for each column.

To interpret the output, note that for light curve 1 we find that the best-fit quadratic has the
form:

```
0.0002554062775*(JD - 53725.173920)*(JD - 53725.173920)
  + 0.0097933162509*(JD - 53725.173920)
  + 10.0830375984825
```

Fitting this equation reduces the RMS from 0.15944 mag to 0.00211 mag (a quadratic signal was
injected into this particular light curve).

## Python Equivalent

```python
import pyvartools as vt
from pyvartools import commands as cmd

lc = vt.LightCurve.from_file("EXAMPLES/1")

result = vt.Pipeline([
    cmd.rms(),
    cmd.decorr(nglobal=1, nlc=1, lcorder=2, correct=True, subtractfirst=True),
    cmd.rms(),
]).run(lc)

print(result.stats)
```
