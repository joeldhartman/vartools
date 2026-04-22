# Statistics

Commands for computing variability and scatter statistics on light curves.

---

## `-rms`

```
-rms
    ["maskpoints" maskvar]
```

Calculate the RMS of the light curves. The output includes the RMS, the mean magnitude, the expected RMS (derived from the formal photometric uncertainties), and the number of points in the light curve.

**Parameters**

- `"maskpoints" maskvar` — Optional. Only points with `maskvar > 0` are included in the calculation; all others are excluded.

**Examples**

**Example 1.** Calculate the mean magnitude, RMS, and expected RMS based on the formal magnitude uncertainties for all light curves in a list file.

```bash
vartools -l EXAMPLES/lc_list -header -rms
```

Output:
```
#Name Mean_Mag_0 RMS_0 Expected_RMS_0 Npoints_0
EXAMPLES/1   ...
EXAMPLES/2   ...
...EXAMPLES/10  ...
```

---

## `-rmsbin`

```
-rmsbin
    Nbin bintime1...bintimeN
    ["maskpoints" maskvar]
```

Calculate the RMS after applying a moving mean filter to the light curves. Similar to [`-chi2bin`](statistics.md#-chi2bin), this measures the correlated (red) noise component by binning on specified timescales. `Nbin` filters are applied, each producing a separate RMS estimate. The note that light curves passed to the next command are **unchanged** by this command.

**Parameters**

- `Nbin` — Number of time bins (filters) to apply.
- `bintime1...bintimeN` — Half-widths of each moving mean filter, in **minutes**. The full filter window for filter *i* is `2.0 * bintimei`.
- `"maskpoints" maskvar` — Optional. Only points with `maskvar > 0` are included.

**Examples**

**Example 1.** Apply moving-mean filters to the light curves and compute statistical measures for each filter. The filters operate by replacing each point in the light curve with the mean of all points that are within the specified number of minutes of that point. This example uses five time-window filters (5.0, 10.0, 60.0, 1440.0, and 14400.0 minutes). The output table displays computed RMS values alongside expected RMS values assuming white noise for each binning window. As the filtering window increases, RMS values generally decrease.

```bash
vartools -l EXAMPLES/lc_list -header -rmsbin 5 5.0 10.0 60.0 1440.0 14400.0
```

---

## `-chi2`

```
-chi2
    ["maskpoints" maskvar]
```

Calculate chi-squared per degree of freedom (χ²/dof) for the light curves. The output includes χ²/dof and the error-weighted mean magnitude.

**Parameters**

- `"maskpoints" maskvar` — Optional. Only points with `maskvar > 0` are included.

**Examples**

**Example 1.** Calculate chi-squared per degree of freedom and the weighted mean magnitude for all light curves in a list.

```bash
vartools -header -l EXAMPLES/lc_list -chi2
```

Output:
```
#Name Chi2_0 Weighted_Mean_Mag_0
EXAMPLES/1  34711.71793  10.24430
EXAMPLES/2   1709.50065  10.11178
EXAMPLES/3     27.06322  10.16684
EXAMPLES/4      5.19874  10.35137
EXAMPLES/5      8.26418  10.43932
EXAMPLES/6      3.94650  10.52748
EXAMPLES/7     10.39941  10.56951
EXAMPLES/8      4.19887  10.61132
EXAMPLES/9      2.67020  10.73129
EXAMPLES/10      3.72218  10.87763
```

---

## `-chi2bin`

```
-chi2bin
    Nbin bintime1...bintimeNbin
    ["maskpoints" maskvar]
```

Calculate χ²/dof after applying a moving mean filter to the light curves. As with `-rmsbin`, the light curves passed to the next command are unchanged. `Nbin` filters are used, producing `Nbin` separate estimates of χ²/dof and the error-weighted mean.

**Parameters**

- `Nbin` — Number of filters.
- `bintime1...bintimeNbin` — Half-widths of the moving mean filters, in **minutes**. The full window for filter *i* is `2.0 * bintimei`.
- `"maskpoints" maskvar` — Optional. Only points with `maskvar > 0` are included.

**Examples**

**Example 1.** Apply moving-mean filters to light curves and calculate chi-squared per degree of freedom along with the weighted mean magnitude for each filter. This example uses 5 filters with durations of 5.0, 10.0, 60.0, 1440.0, and 14400.0 minutes. As formal errors decrease according to white noise expectations, chi-squared values increase with larger filter sizes when red noise is present.

```bash
vartools -l EXAMPLES/lc_list -header -chi2bin 5 5.0 10.0 60.0 1440.0 14400.0
```

---

## `-stats`

```
-stats
    var1,var2,... stats1,stats2,...
    ["maskpoints" maskvar]
```

Compute one or more general statistics on one or more light-curve vectors (e.g., `t`, `mag`, `err`, or any user-defined variable). Every requested statistic is computed for every listed variable.

**Parameters**

- `var1,var2,...` — Comma-separated list of variable names to compute statistics on.
- `stats1,stats2,...` — Comma-separated list of statistics to compute. Every statistic is computed for every variable.
- `"maskpoints" maskvar` — Optional. Only points with `maskvar > 0` are used in the calculations.

**Available statistics strings**

| String | Description |
|--------|-------------|
| `mean` | Arithmetic mean |
| `weightedmean` | Mean weighted by 1/σ² |
| `median` | Median |
| `wmedian` | Median weighted by light curve uncertainties |
| `stddev` | Standard deviation with respect to the mean |
| `meddev` | Standard deviation with respect to the median |
| `medmeddev` | Median of the absolute deviations from the median |
| `MAD` | 1.483 × medmeddev. Equals stddev for a Gaussian distribution in the large-N limit |
| `kurtosis` | Kurtosis |
| `skewness` | Skewness |
| `pct%f` | The %f-th percentile, where %f is a floating point number between 0 and 100 (e.g., `pct25`) |
| `wpct%f` | Percentile including light curve uncertainties as weights |
| `max` | Maximum value (equivalent to `pct100`) |
| `min` | Minimum value (equivalent to `pct0`) |
| `sum` | Sum of all elements in the vector |

**Example**

```bash
vartools -l EXAMPLES/lc_list -stats mag,err mean,stddev,MAD
```

**Examples**

**Example 1.** Calculate various statistical measures for light curve magnitudes and magnitudes after adding Gaussian noise. The `-expr` parameter defines a new variable with noise added, while `-stats` specifies which variables and statistics to compute. Percentile statistics (`pct##`) represent the specified percentile values.

```bash
vartools -i EXAMPLES/3 \
    -oneline \
    -expr 'mag2=mag+0.01*gauss()' \
    -stats mag,mag2 \
        mean,weightedmean,median,stddev,meddev,medmeddev,MAD,kurtosis,skewness,pct10,pct20,pct80,pct90,max,min,sum
```

---

## `-alarm`

```
-alarm
    ["maskpoints" maskvar]
```

Calculate the alarm variability statistic for each light curve (Tamuz, Mazeh, and North 2006, MNRAS, 367, 1521). This statistic is designed to detect time-correlated variability.

**Parameters**

- `"maskpoints" maskvar` — Optional. Points with `maskvar > 0` are included; others are excluded.

**Citation:** Tamuz, Mazeh, and North 2006, MNRAS, 367, 1521.

---

## `-Jstet`

```
-Jstet
    timescale dates ["maskpoints" maskvar]
```

Calculate Stetson's J statistic, L statistic, and the kurtosis for each light curve. The J statistic measures time-correlated variability by comparing pairs of observations that are close in time.

**Parameters**

- `timescale` — Time in minutes that distinguishes between "near" (correlated) and "far" (uncorrelated) observation pairs.
- `dates` — File containing JDs for all possible observations in the first column. This is used to compute the maximum possible weight. Note: the J statistic here includes an extra factor of `(sum(weights)/weight_max)` compared to Stetson's original definition.
- `"maskpoints" maskvar` — Optional. Only points with `maskvar > 0` are included.

**Citation:** Stetson, P.B. 1996, PASP, 108, 851.

**Examples**

**Example 1.** Calculate Stetson's J statistic, L statistic, and kurtosis for all light curves in a list, using 0.5 days to distinguish between "near" and "far" observations.

```bash
vartools -l EXAMPLES/lc_list -header \
    -Jstet 0.5 EXAMPLES/dates_tfa
```

Output:
```
#Name Jstet_0 Kurtosis_0 Lstet_0
EXAMPLES/1  98.13279   0.96779  94.97154
EXAMPLES/2  30.19309   0.94719  28.59852
EXAMPLES/3   0.65597   0.92816   0.60885
EXAMPLES/4   0.34402   0.84500   0.29070
EXAMPLES/5   0.58730   0.92120   0.54102
EXAMPLES/6   0.34455   0.93794   0.32317
EXAMPLES/7   0.41754   0.92501   0.38623
EXAMPLES/8   0.46381   0.96124   0.44583
EXAMPLES/9   0.22075   0.80997   0.17880
EXAMPLES/10   0.25784   0.92806   0.23929
```

---

## `-autocorrelation`

```
-autocorrelation
    start stop step outdir ["maskpoints" maskvar]
```

Calculate the discrete auto-correlation function (Edelson and Krolik 1988, ApJ, 333, 646) for each light curve. The results are written to files in `outdir` with the suffix `.autocorr` (i.e., `outdir/$basename.autocorr`).

!!! note "Cross-reference"
    For period-finding using autocorrelation-based methods, see the [Period Finding](period-finding.md) page.

**Parameters**

- `start` — Start time for sampling the autocorrelation, in days.
- `stop` — Stop time for sampling the autocorrelation, in days.
- `step` — Step size for sampling, in days.
- `outdir` — Directory where output `.autocorr` files are written.
- `"maskpoints" maskvar` — Optional. Only points with `maskvar > 0` are included.

**Notes**

Rather than using the variance in the denominator (as in the Edelson and Krolik formula), the formal uncertainty is used. This avoids imaginary numbers when measurement errors are overestimated. To use the variance in the denominator instead, issue `-changeerror` before calling this command.

Due to binning, when the variance is used in the denominator the autocorrelation function may be smaller than 1 unless the time step is less than the minimum time difference between consecutive measurements.

**Citation:** Edelson, R.A. & Krolik, J.H. 1988, ApJ, 333, 646.

**Examples**

**Example 1.** Compute the discrete auto-correlation function (DACF) of a single light curve spanning time lags from 0 to 10.0 days with a step of 0.05 days. Output is written to `EXAMPLES/OUTDIR1/2.autocorr`.

```bash
vartools -i EXAMPLES/2 -header \
    -autocorrelation 0.0 10. 0.05 EXAMPLES/OUTDIR1
```

Output:
```
#Name
EXAMPLES/2
```
