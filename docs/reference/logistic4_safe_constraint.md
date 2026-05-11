# Compute Safe Parameter Bounds for 4PL (logistic4) Fitting

Builds lower and upper bounds for all free parameters in the
[`logistic4()`](https://immunoplex.github.io/curveRfreq/reference/logistic4.md)
model using
[`adaptive_constraint_profile()`](https://immunoplex.github.io/curveRfreq/reference/adaptive_constraint_profile.md).

## Usage

``` r
logistic4_safe_constraint(
  data,
  y_min = 1,
  y_max,
  logistic4_formula,
  logistic4_free_vars,
  is_log_response,
  is_log_concentration,
  antigen_settings,
  constraint_profile = NULL
)
```

## Arguments

- data:

  Data frame. Must contain a `concentration` column and the response
  column.

- y_min:

  Numeric. Minimum observed response. Default `1`.

- y_max:

  Numeric. Maximum observed response.

- logistic4_formula:

  Formula. The logistic4 model formula.

- logistic4_free_vars:

  Character vector. Names of free parameters.

- is_log_response:

  Logical. Whether the response is \\\log\_{10}\\-transformed.

- is_log_concentration:

  Logical. Whether concentration is on log scale.

- antigen_settings:

  List. Must contain `l_asy_min_constraint` and `l_asy_max_constraint`.

- constraint_profile:

  List or `NULL`. Output of
  [`adaptive_constraint_profile()`](https://immunoplex.github.io/curveRfreq/reference/adaptive_constraint_profile.md).
  If `NULL`, one is built internally.

## Value

A list with `lower` and `upper` named numeric vectors.

## See also

Other safe-constraints:
[`adaptive_constraint_profile()`](https://immunoplex.github.io/curveRfreq/reference/adaptive_constraint_profile.md),
[`gompertz4_safe_constraint()`](https://immunoplex.github.io/curveRfreq/reference/gompertz4_safe_constraint.md),
[`logistic5_safe_constraint()`](https://immunoplex.github.io/curveRfreq/reference/logistic5_safe_constraint.md),
[`loglogistic4_safe_constraint()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic4_safe_constraint.md),
[`loglogistic5_safe_constraint()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic5_safe_constraint.md)
