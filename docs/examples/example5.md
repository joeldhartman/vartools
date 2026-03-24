# Example 5: RR Lyrae Injection and Recovery

Injecting an RR Lyrae signal into a light curve and recovering it.

## Commands

### Step 1: Fit a Fourier series to a real RR Lyrae light curve

```bash
./vartools -i EXAMPLES/M3.V006.lc -Killharm fix 1 0.514333 10 0 1 \
   EXAMPLES/OUTDIR1/ fitonly outRphi -header
```

Output:

```
#Name Killharm_Mean_Mag_0 Killharm_Period_1_0 Killharm_Per1_Fundamental_Amp_2_0 \
Killharm_Per1_Fundamental_Phi_2_0 Killharm_Per1_Harm_R_2_1_0 Killharm_Per1_Harm_Phi_2_1_0 \
Killharm_Per1_Harm_R_3_1_0 Killharm_Per1_Harm_Phi_3_1_0 Killharm_Per1_Harm_R_4_1_0 \
Killharm_Per1_Harm_Phi_4_1_0 Killharm_Per1_Harm_R_5_1_0 Killharm_Per1_Harm_Phi_5_1_0 \
Killharm_Per1_Harm_R_6_1_0 Killharm_Per1_Harm_Phi_6_1_0 Killharm_Per1_Harm_R_7_1_0 \
Killharm_Per1_Harm_Phi_7_1_0 Killharm_Per1_Harm_R_8_1_0 Killharm_Per1_Harm_Phi_8_1_0 \
Killharm_Per1_Harm_R_9_1_0 Killharm_Per1_Harm_Phi_9_1_0 Killharm_Per1_Harm_R_10_1_0 \
Killharm_Per1_Harm_Phi_10_1_0 Killharm_Per1_Harm_R_11_1_0 Killharm_Per1_Harm_Phi_11_1_0 \
Killharm_Per1_Amplitude_0
EXAMPLES/M3.V006.lc  15.77123     0.51433300   0.38041  -0.07662   0.47077   0.60826   0.35917   \
0.26249   0.23631  -0.06843   0.16353   0.60682   0.10621   0.28738   0.06203   0.95751   0.03602   \
0.58867   0.02900   0.22322   0.01750   0.94258   0.00768   0.66560   1.11128
```

### Step 2: Inject the fitted signal at varying amplitudes and attempt recovery

```bash
echo EXAMPLES/4 | gawk '{amp = 0.25; for(i=1; i <= 10; i += 1) {print $1, amp; amp = amp*0.5;}}' \
     | ./vartools -l - -Injectharm fix 0.514333 10 \
        amplist phaserand \
        ampfix 0.47077 amprel phasefix 0.60826 phaserel \
        ampfix 0.35916 amprel phasefix 0.26249 phaserel \
        ampfix 0.23631 amprel phasefix -0.06843 phaserel \
        ampfix 0.16353 amprel phasefix 0.60682 phaserel \
        ampfix 0.10621 amprel phasefix 0.28738 phaserel \
        ampfix 0.06203 amprel phasefix 0.95751 phaserel \
        ampfix 0.03602 amprel phasefix 0.58867 phaserel \
        ampfix 0.02900 amprel phasefix 0.22322 phaserel \
        ampfix 0.01750 amprel phasefix 0.94258 phaserel \
        ampfix 0.00768 amprel phasefix 0.66560 phaserel \
        0 0 \
        -LS 0.1 10.0 0.01 2 0 \
        -aov_harm 2 0.1 10.0 0.1 0.01 2 0 -header -numbercolumns
```

Output (first two and last two lines shown):

```
#1_Name 2_Injectharm_Period_0 3_Injectharm_Fundamental_Amp_0 4_Injectharm_Fundamental_Phase_0 \
5_Injectharm_Harm_2_Amp_0 6_Injectharm_Harm_2_Phase_0 7_Injectharm_Harm_3_Amp_0 8_Injectharm_Harm_3_Phase_0 \
9_Injectharm_Harm_4_Amp_0 10_Injectharm_Harm_4_Phase_0 11_Injectharm_Harm_5_Amp_0 12_Injectharm_Harm_5_Phase_0 \
13_Injectharm_Harm_6_Amp_0 14_Injectharm_Harm_6_Phase_0 15_Injectharm_Harm_7_Amp_0 16_Injectharm_Harm_7_Phase_0 \
17_Injectharm_Harm_8_Amp_0 18_Injectharm_Harm_8_Phase_0 19_Injectharm_Harm_9_Amp_0 20_Injectharm_Harm_9_Phase_0 \
21_Injectharm_Harm_10_Amp_0 22_Injectharm_Harm_10_Phase_0 23_Injectharm_Harm_11_Amp_0 24_Injectharm_Harm_11_Phase_0 \
25_LS_Period_1_1 26_Log10_LS_Prob_1_1 27_LS_SNR_1_1 28_LS_Period_2_1 \
29_Log10_LS_Prob_2_1 30_LS_SNR_2_1 31_Period_1_2 32_AOV_HARM_1_2 \
33_AOV_HARM_SNR_1_2 34_AOV_HARM_NEG_LOG_FAP_1_2 35_Period_2_2 36_AOV_HARM_2_2 \
37_AOV_HARM_SNR_2_2 38_AOV_HARM_NEG_LOG_FAP_2_2 39_Mean_lnAOV_2 40_RMS_lnAOV_2 \
EXAMPLES/4 0.51433300 0.25000 0.84019 0.11769 2.28864 0.08979 2.78305 \
0.05908 3.29232 0.04088 4.80776 0.02655 5.32851 0.01551 6.83882 \
0.00901 7.31017 0.00725 7.78491 0.00438 9.34446 0.00192 9.90766 \
0.51493297 -492.06805 25.15612 1.06567664 -378.44559 19.24304 0.51432840 4322.49 \
149.531 2969.03 0.51500100 3826.26 132.25 2805.09 28.5848 28.7159
EXAMPLES/4 0.51433300 0.12500 0.39438 0.05885 1.39703 0.04489 1.44564 \
0.02954 1.50910 0.02044 2.57873 0.01328 2.65368 0.00775 3.71819 \
0.00450 3.74373 0.00363 3.77267 0.00219 4.88641 0.00096 5.00381 \
0.51408199 -487.91541 25.18381 0.33926383 -379.53859 19.49478 0.51432840 4609.61 \
137.471 3056.79 0.51414979 4585.38 136.744 3049.56 30.8468 33.3071
...
...
EXAMPLES/4 0.51433300 0.00098 0.27777 0.00046 1.16381 0.00035 1.09581 \
0.00023 1.04267 0.00016 1.99569 0.00010 1.95403 0.00006 2.90193 \
0.00004 2.81087 0.00003 2.72319 0.00002 3.72033 0.00001 3.72112 \
0.93527063 -69.42762 7.34581 0.33952304 -67.04912 7.08529 1.99136298 166.558 \
9.40744 291.44 1.98761005 157.378 8.84087 276.201 14.1343 16.2025
EXAMPLES/4 0.51433300 0.00049 0.55397 0.00023 1.71620 0.00018 1.92440 \
0.00012 2.14745 0.00008 3.37667 0.00005 3.61120 0.00003 4.83530 \
0.00002 5.02043 0.00001 5.20895 0.00001 6.48228 0.00000 6.75927 \
0.93667874 -57.95590 7.76038 1.04106764 -55.69804 7.45058 0.99348763 199.742 \
13.8977 345.37 0.99504450 188.754 13.0856 327.71 11.7027 13.5302
```

## Explanation

This example conducts a simple test of injecting an RR Lyrae signal into a light curve with a
fixed period and a range of amplitudes, then attempting to recover the signal.

**Step 1** fits a Fourier series to the light curve of an actual fundamental mode RR Lyrae star
(M3 Variable V006). The `-Killharm` command fits the fundamental mode plus 10 harmonics, with
the period fixed to the known value of 0.514333 days. The `outRphi` keyword outputs the Fourier
components as R_k1 = (A_k / A_1) and phi_k1 = (phi_k - k*phi_1) rather than the default sin
and cos coefficients. This format is convenient for use as a model signal to inject at varying
amplitudes while preserving the signal shape.

**Step 2** generates an input list via `echo`/`gawk` with the light curve name (`EXAMPLES/4`)
in the first column and the injection amplitude (halved at each step, from 0.25 down to ~0.00049)
in the second column. This is piped into vartools using `-l -` (reading the list from stdin).

Three commands are run on each iteration:

- `-Injectharm`: Injects the model Fourier series with a fixed period of 0.514333 days and
  10 harmonics. The fundamental amplitude is read from the input list (`amplist`) and the
  phase is chosen randomly (`phaserand`). The 10 harmonics use fixed R_k1 and phi_k1 values
  from the Fourier fit to the real RR Lyrae, specified with `amprel` and `phaserel`.
- `-LS`: Attempts period recovery with Lomb-Scargle over 0.1–10.0 days.
- `-aov_harm`: Attempts period recovery with multi-harmonic AoV over the same range.

Columns 3–24 of the output give the amplitude (A_k) and phase (phi_k) of each harmonic in the
injected signal. The most important is column 3, the amplitude of the fundamental mode.

The top L-S period (column 25) is recovered near 0.514–0.515 days for the first 8 simulations
(down to a fundamental amplitude of ~1.95 mmag, corresponding to a peak-to-peak amplitude of
~5.7 mmag). For the two lowest amplitude simulations (~0.98 and ~0.49 mmag fundamental, or
~2.86 and ~1.43 mmag peak-to-peak) the signal is not recovered with L-S. The same threshold
applies for multi-harmonic AoV (column 31), but AoV generally recovers a period closer to the
true value since the signal is not purely sinusoidal.

Note: the exact numbers you obtain may differ depending on your random number generator. To
obtain non-repeatable random phases use the `-randseed time` option.

## Python Equivalent

```python
import pyvartools as vt
from pyvartools import commands as cmd
import numpy as np

# Generate amplitude list (halving each step, 10 iterations)
amplitudes = [0.25 * 0.5**i for i in range(10)]

results = []
for amp in amplitudes:
    lc = vt.LightCurve.from_file("EXAMPLES/4")
    result = vt.Pipeline([
        cmd.Injectharm(
            period_mode="fix", period=0.514333, nharmonics=10,
            fundamental_amp_mode="ampfix", fundamental_amp=amp,
            fundamental_phase_mode="phaserand",
            harmonics=[
                {"ampfix": 0.47077, "amprel": True, "phasefix": 0.60826, "phaserel": True},
                # ... (remaining harmonics)
            ],
        ),
        cmd.LS(0.1, 10.0, 0.01, Npeaks=2, output_periodogram=False),
        cmd.aov_harm(2, 0.1, 10.0, 0.1, 0.01, Npeaks=2, output_periodogram=False),
    ]).run(lc)
    results.append(result.stats)

for r in results:
    print(r)
```
