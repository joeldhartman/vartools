# Example 3: Transit Injection

Injecting a transit signal into a light curve.

## Commands

```bash
./vartools -i EXAMPLES/3 -MandelAgolTransit 2.12345 53725.174 0.1 10.0 1.0 0. 0. 0 quad 0.236 0.391 \
  0 0 0 0 0 0 0 0 0 0 0 0 1 EXAMPLES/OUTDIR1 -tab

gawk '{print $1, $2 + $3, $4}' EXAMPLES/OUTDIR1/3.mandelagoltransit.model > EXAMPLES/3.transit
```

## Output

```
Name MandelAgolTransit_Period_0 MandelAgolTransit_T0_0 MandelAgolTransit_r_0 MandelAgolTransit_a_0 \
MandelAgolTransit_sin_i_0 MandelAgolTransit_e_0 MandelAgolTransit_omega_0 MandelAgolTransit_mconst_0 \
MandelAgolTransit_ldcoeff1_0 MandelAgolTransit_ldcoeff2_0 MandelAgolTransit_chi2_0
---- -------------------------- ---------------------- --------------------- --------------------- \
------------------------- --------------------- ------------------------- -------------------------- \
---------------------------- ---------------------------- ------------------------
EXAMPLES/3 2.12345000 53725.17400000 0.10000 10.00000 1.00000 0.00000 0.00000 0.00000 0.23600 \
0.39100 129730808.07162
```

## Explanation

These two commands effectively inject a Mandel & Agol transit into a light curve. The first
command calls the `-MandelAgolTransit` vartools command to create a simulated transit with:

- Period: 2.12345 days
- T0: 53725.174
- R_planet / R_star: 0.1
- a / R_star: 10.0
- sin(i): 1.0
- Circular orbit (e = 0, omega = 0)
- Quadratic limb darkening with coefficients 0.236 and 0.391

We do not fit for any of the parameters (10 zeros), nor do we fit an RV curve (the 11th zero),
nor do we correct the light curve (the 12th zero). The model is output in the third column of
`EXAMPLES/OUTDIR1/3.mandelagoltransit.model`.

In the second command we use `gawk` (or `awk`) to add the model transit signal to the light
curve and output the result to `EXAMPLES/3.transit`.

Note that this procedure may now also be done using the `-Injecttransit` command so that the
external call to `gawk` is not needed.

## Python Equivalent

```python
import pyvartools as vt
from pyvartools import commands as cmd

lc = vt.LightCurve.from_file("EXAMPLES/3")

result = vt.Pipeline([
    cmd.MandelAgolTransit(
        period=2.12345, T0=53725.174, r=0.1, a=10.0, sini=1.0,
        e=0.0, omega=0.0,
        limb_darkening="quad", ld1=0.236, ld2=0.391,
        fit=False, correct=False,
        outdir="EXAMPLES/OUTDIR1",
    ),
]).run(lc)

# Add model to the light curve to produce the injected transit
# (equivalent to the gawk step above, or use -Injecttransit instead)
print(result.stats)
```
