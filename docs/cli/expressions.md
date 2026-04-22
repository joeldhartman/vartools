# Analytic Expressions

VARTOOLS includes a built-in expression evaluator used by several commands including [`-expr`](manipulation.md#-expr), [`-linfit`](model-fitting.md#-linfit), [`-nonlinfit`](model-fitting.md#-nonlinfit), [`-if`](control-flow.md), [`-restricttimes`](manipulation.md#-restricttimes), and the `"expr"` keyword for per-LC parameter values. This page documents all supported operators, functions, constants, and special variables.

Run `vartools -functionlist` to print this information from the command line.

---

## Variables

Expressions can reference any named variable in the VARTOOLS variable registry:

- **Light curve vectors** — `t`, `mag`, `err`, and any additional columns declared with `-inputlcformat`. These have one value per observation.
- **Output statistics** — Variables created by prior commands, named after their output column headers (with leading digits and underscores stripped, and non-alphanumeric characters replaced with `_`). Use `vartools -headeronly ...` to see the names. For example, the header `2_STATS_mag_PCT20.00_0` becomes the variable `STATS_mag_PCT20_00_0`.
- **User-defined variables** — Variables created by `-expr`, `-inlistvars`, or `-inputlcformat col=0`.

---

## Operators

| Operator | Description |
|----------|-------------|
| `a + b` | Addition |
| `a - b` | Subtraction |
| `a * b` | Multiplication |
| `a / b` | Division |
| `a % b` | Floating-point remainder (`fmod`) |
| `a ^ b` | Exponentiation |
| `a > b` | Greater than (returns 1 or 0) |
| `a >= b` | Greater than or equal |
| `a < b` | Less than |
| `a <= b` | Less than or equal |
| `a == b` | Logical equals |
| `a != b` | Logical not equal |
| `a && b` | Logical AND |
| `a \|\| b` | Logical OR |
| `!a` | Logical NOT |
| `~a` | Bitwise complement |
| `a & b` | Bitwise AND |
| `a \| b` | Bitwise inclusive OR |

Standard operator precedence applies. Use parentheses to override.

---

## Array indexing

| Syntax | Description |
|--------|-------------|
| `a[b]` | Returns the `b`-th element of vector `a` (0-indexed). `b` can be any expression and is rounded to the nearest integer. Out-of-bounds indices return 0.0 silently. |

When used on the left-hand side of `-expr` (e.g., `mag[0]=5.0`), indexing supports single indices, ranges (`mag[0:5]`), and boolean filters (`mag[t>5.0]`). See `vartools -help -expr` for details.

---

## Scalar functions

These operate element-wise on their arguments.

### Mathematical

| Function | Description |
|----------|-------------|
| `exp(x)` | Exponential |
| `log(x)` | Natural logarithm |
| `log10(x)` | Base-10 logarithm |
| `sqrt(x)` | Square root |
| `abs(x)` | Absolute value |
| `max(x, y)` | Larger of two values |
| `min(x, y)` | Smaller of two values |
| `hypot(x, y)` | `sqrt(x² + y²)` |

### Trigonometric (radians)

| Function | Description |
|----------|-------------|
| `sin(x)` | Sine |
| `cos(x)` | Cosine |
| `tan(x)` | Tangent |
| `asin(x)` | Inverse sine |
| `acos(x)` | Inverse cosine |
| `atan2(y, x)` | Four-quadrant inverse tangent |

### Trigonometric (degrees)

| Function | Description |
|----------|-------------|
| `sindegr(x)` | Sine (input in degrees) |
| `cosdegr(x)` | Cosine (input in degrees) |
| `tandegr(x)` | Tangent (input in degrees) |
| `asindegr(x)` | Inverse sine (output in degrees) |
| `acosdegr(x)` | Inverse cosine (output in degrees) |
| `atan2degr(y, x)` | Four-quadrant inverse tangent (output in degrees) |

### Hyperbolic

| Function | Description |
|----------|-------------|
| `sinh(x)` | Hyperbolic sine |
| `cosh(x)` | Hyperbolic cosine |
| `tanh(x)` | Hyperbolic tangent |
| `asinh(x)` | Inverse hyperbolic sine |
| `acosh(x)` | Inverse hyperbolic cosine |
| `atanh(x)` | Inverse hyperbolic tangent |

### Special functions and rounding

| Function | Description |
|----------|-------------|
| `erf(x)` | Error function |
| `erfc(x)` | Complementary error function |
| `lgamma(x)` | Log-gamma function |
| `gamma(x)` | Gamma function |
| `theta(x)` | Heaviside step: 1 for x ≥ 0, 0 for x < 0 |
| `round(x)` | Round to nearest integer |
| `ceil(x)` | Ceiling (smallest integer ≥ x) |
| `floor(x)` | Floor (largest integer ≤ x) |
| `isnan(x)` | Returns 1 if x is NaN, 0 otherwise |

### Random number generators

| Function | Description |
|----------|-------------|
| `rand()` | Uniform random number in [0, 1) |
| `gauss()` | Gaussian random number (mean 0, variance 1) |

---

## Vector functions

| Function | Description |
|----------|-------------|
| `len(x)` | Length of vector `x`. The argument must be a single variable name (not a compound expression). |

---

## Aggregate functions

Aggregate functions operate over **all observations** in the current light curve and return a single scalar value. They are especially useful with [`-expr listvar`](manipulation.md#-expr) to compute per-star summary statistics.

All aggregate functions accept **arbitrary expressions** as arguments — not just bare variable names. For example, `mean(mag*mag)`, `sum(sin(2*pi*t/P))`, and `stddev(mag-model)` are all valid.

Most aggregate functions accept an optional trailing **filter** argument. When provided, only observations where the filter expression evaluates to a value > 0 are included. For example, `mean(mag, t>53730)` computes the mean of `mag` using only observations with `t > 53730`.

### Basic statistics

| Function | Description |
|----------|-------------|
| `mean(x [, filter])` | Arithmetic mean |
| `median(x [, filter])` | Median |
| `stddev(x [, filter])` | Standard deviation |
| `meddev(x [, filter])` | Median deviation from the mean |
| `medmeddev(x [, filter])` | Median of absolute deviations from the median |
| `MAD(x [, filter])` | Median absolute deviation (1.483 × medmeddev) |
| `kurtosis(x [, filter])` | Excess kurtosis |
| `skewness(x [, filter])` | Skewness |

### Extrema and summation

| Function | Description |
|----------|-------------|
| `vmin(x [, filter])` | Minimum value. Named `vmin` to distinguish from the two-argument scalar `min(a, b)`. |
| `vmax(x [, filter])` | Maximum value. Named `vmax` to distinguish from the two-argument scalar `max(a, b)`. |
| `sum(x [, filter])` | Sum of all values |

### Weighted statistics

| Function | Description |
|----------|-------------|
| `weightedmean(x, w [, filter])` | Weighted mean of `x` with weights `w` (inverse-variance weighting) |
| `wmedian(x, w [, filter])` | Weighted median |

### Percentiles

| Function | Description |
|----------|-------------|
| `pct(x, pctval [, filter])` | Percentile of `x`. `pctval` is the percentile as a percentage (e.g. 50.0 for the median, 95.0 for the 95th percentile). |
| `wpct(x, w, pctval [, filter])` | Weighted percentile |

### Examples

```bash
# Per-star mean magnitude
vartools -l EXAMPLES/lc_list -expr listvar 'avg=mean(mag)' -oneline

# Mean of squared magnitudes
vartools -i EXAMPLES/2 -expr listvar 'meansq=mean(mag*mag)' -oneline

# Mean of mag for observations after JD 53730
vartools -i EXAMPLES/2 -expr listvar 'late_avg=mean(mag, t>53730)' -oneline

# 95th percentile of magnitude
vartools -i EXAMPLES/2 -expr listvar 'p95=pct(mag, 95.0)' -oneline

# Weighted mean using error column as weights
vartools -i EXAMPLES/2 -expr listvar 'wmag=weightedmean(mag, err)' -oneline
```

---

## Constants

| Name | Value |
|------|-------|
| `pi` | 3.14159265... |
| `e` | 2.71828182... |

---

## Special variables

| Variable | Description |
|----------|-------------|
| `NR` | Observation index within the current light curve (0-based) |
| `NF` | Light curve index in the input list (0-based) |

---

## Astronomical functions (`astrofuncs`)

The `astrofuncs` user-function library ships with the VARTOOLS distribution and provides orbital mechanics and transit model functions. Load it with `-F astrofuncs` on the command line. See [Extending VARTOOLS — User Functions](../extensions/user-functions.md) for details on the user-function mechanism.

```bash
# Compute the eccentric anomaly for each observation's mean-anomaly value.
# Any analytic-expression argument may reference light-curve vectors (t, mag,
# err, ...), constants, or variables defined by preceding commands.
vartools -F astrofuncs -i EXAMPLES/2 \
    -expr 'meananom=(t-53725.0)/1.234 * 2*pi' \
    -expr 'eccanom=EccentricAnomaly(meananom, 0.3)' \
    -rms -oneline
```

### Orbital mechanics

| Function | Description |
|----------|-------------|
| `EccentricAnomaly(M, e)` | Eccentric anomaly (radians) from mean anomaly `M` (radians) and eccentricity `e`. |
| `MeanAnomaly(dt, P)` | Mean anomaly (radians) from time since periastron `dt` and period `P`. |
| `MeanAnomalyConjunction(dt, P, e, omega)` | Mean anomaly (radians) from time since conjunction (transit). `omega` in degrees. |

### Transit models

| Function | Description |
|----------|-------------|
| `TransitQuadLD(dt, P, b, Rp/R*, a/R*, e, omega, u1, u2)` | Relative flux for a Mandel & Agol (2002) transit with quadratic limb darkening. Returns 1.0 out of transit, < 1.0 in transit. `dt` = time since transit center; `omega` in degrees. |
| `TransitNonlinLD(dt, P, b, Rp/R*, a/R*, e, omega, a1, a2, a3, a4)` | Same as above but with a 4-parameter non-linear limb darkening law. |

### Rossiter-McLaughlin effect

| Function | Description |
|----------|-------------|
| `BroadeningProfile(delv, dt, P, lambda, vsini, b, Rp/R*, a/R*, e, omega, u1, u2)` | Distorted stellar line broadening function during transit. `delv` = `(wl - wl0) * c / wl0` in km/s; `lambda` = projected obliquity (degrees); `vsini` in km/s. |
| `TransitProjectedX(dt, P, lambda, b, Rp/R*, a/R*, e, omega)` | Sky-projected X position of the planet center in units of the stellar radius (rotation axis along Y). |
| `TransitProjectedY(dt, P, lambda, b, Rp/R*, a/R*, e, omega)` | Sky-projected Y position of the planet center. |

### Radial velocity

| Function | Description |
|----------|-------------|
| `RV_E(E, e, omega, K)` | Radial velocity given eccentric anomaly `E` (radians). `omega` in degrees; result in same units as `K`. |
| `RV_M(M, e, omega, K)` | Radial velocity given mean anomaly `M` (radians). |
| `RV_dt(dt, P, e, omega, K)` | Radial velocity given time since conjunction `dt` and period `P`. |
| `RV_dtp(dtp, P, e, omega, K)` | Radial velocity given time since periastron `dtp` and period `P`. |
