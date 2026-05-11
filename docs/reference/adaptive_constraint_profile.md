# Build an Adaptive Constraint Profile from Observed Data

Inspects the response range, dynamic range, and scale to choose
appropriate bounds for nonlinear optimisation. The returned profile is
consumed by every `*_safe_constraint()` function.

## Usage

``` r
adaptive_constraint_profile(
  data,
  response_variable,
  is_log_response,
  antigen_settings
)
```

## Arguments

- data:

  Data frame with at least columns named by `response_variable` and
  `"concentration"`.

- response_variable:

  Character. Column name for the response.

- is_log_response:

  Logical. `TRUE` if data is already \\\log\_{10}\\-transformed.

- antigen_settings:

  List with constraint metadata, including `l_asy_min_constraint` and
  `l_asy_max_constraint`.

## Value

A named list:

- y_min, y_max:

  Observed response extremes (finite values only).

- dynamic_range:

  `y_max - y_min`.

- conc_range:

  Range of observed concentrations.

- scale_class:

  `"high"` (MFI-like, e.g. Luminex), `"medium"`, or `"low"`
  (OD/absorbance-like).

- slope_min, slope_max:

  Adapted bounds for scale parameter `b`.

- g_min, g_max:

  Adapted bounds for asymmetry parameter `g`.

- conc_pad_frac:

  Fraction beyond concentration range for inflection-point `c` bounds.

- d_margin_frac:

  Fraction of `y_max` for upper-asymptote `d` bounds.

## Details

Three scale classes are recognised:

|        |                |                 |                     |
|--------|----------------|-----------------|---------------------|
| Class  | Raw threshold  | Log10 threshold | Typical assay       |
| high   | `y_max > 1000` | `y_max > 2.5`   | Luminex MFI         |
| medium | `y_max > 10`   | `y_max > 0.5`   | moderate signals    |
| low    | otherwise      | otherwise       | ELISA OD/absorbance |

Narrower dynamic ranges (low scale class) receive wider slope and
asymmetry bounds to avoid near-singular Jacobians during
[`stats::nls()`](https://rdrr.io/r/stats/nls.html) fitting.

## See also

[`logistic5_safe_constraint()`](https://immunoplex.github.io/curveRfreq/reference/logistic5_safe_constraint.md),
[`logistic4_safe_constraint()`](https://immunoplex.github.io/curveRfreq/reference/logistic4_safe_constraint.md),
[`compute_robust_curves()`](https://immunoplex.github.io/curveRfreq/reference/compute_robust_curves.md)

Other safe-constraints:
[`gompertz4_safe_constraint()`](https://immunoplex.github.io/curveRfreq/reference/gompertz4_safe_constraint.md),
[`logistic4_safe_constraint()`](https://immunoplex.github.io/curveRfreq/reference/logistic4_safe_constraint.md),
[`logistic5_safe_constraint()`](https://immunoplex.github.io/curveRfreq/reference/logistic5_safe_constraint.md),
[`loglogistic4_safe_constraint()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic4_safe_constraint.md),
[`loglogistic5_safe_constraint()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic5_safe_constraint.md)

## Examples

``` r
fake_data <- data.frame(
  mfi = c(100, 500, 5000, 15000, 20000),
  concentration = c(0.01, 0.1, 1, 10, 100)
)
profile <- adaptive_constraint_profile(
  data = fake_data,
  response_variable = "mfi",
  is_log_response = FALSE,
  antigen_settings = list(l_asy_min_constraint = 0,
                          l_asy_max_constraint = 200)
)
str(profile)
#> List of 11
#>  $ y_min        : num 100
#>  $ y_max        : num 20000
#>  $ dynamic_range: num 19900
#>  $ conc_range   : num 100
#>  $ scale_class  : chr "high"
#>  $ slope_max    : num 2
#>  $ slope_min    : num 0.1
#>  $ g_min        : num 0.5
#>  $ g_max        : num 5
#>  $ conc_pad_frac: num 0.5
#>  $ d_margin_frac: num 0.5
```
