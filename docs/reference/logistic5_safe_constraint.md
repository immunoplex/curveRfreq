# Compute Safe Parameter Bounds for 5PL (logistic5) Fitting

Builds lower and upper bounds for all free parameters in the
[`logistic5()`](https://immunoplex.github.io/curveRfreq/reference/logistic5.md)
model, adapted to the observed data scale via
[`adaptive_constraint_profile()`](https://immunoplex.github.io/curveRfreq/reference/adaptive_constraint_profile.md).
Used by
[`compute_robust_curves()`](https://immunoplex.github.io/curveRfreq/reference/compute_robust_curves.md)
to supply bounds to [`stats::nls()`](https://rdrr.io/r/stats/nls.html)
with `algorithm = "port"`.

## Usage

``` r
logistic5_safe_constraint(
  data,
  y_min = 1,
  y_max,
  logistic5_formula,
  logistic5_free_vars,
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

- logistic5_formula:

  Formula. The logistic5 model formula from
  [`select_model_formulas()`](https://immunoplex.github.io/curveRfreq/reference/select_model_formulas.md).

- logistic5_free_vars:

  Character vector. Names of free parameters in the formula (e.g.,
  `c("b", "c", "d", "g")` when `a` is fixed, or
  `c("a", "b", "c", "d", "g")` when free).

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

A list with elements `lower` and `upper`, each a named numeric vector
matching `logistic5_free_vars`.

## See also

[`adaptive_constraint_profile()`](https://immunoplex.github.io/curveRfreq/reference/adaptive_constraint_profile.md),
[`logistic5()`](https://immunoplex.github.io/curveRfreq/reference/logistic5.md),
[`compute_robust_curves()`](https://immunoplex.github.io/curveRfreq/reference/compute_robust_curves.md)

Other safe-constraints:
[`adaptive_constraint_profile()`](https://immunoplex.github.io/curveRfreq/reference/adaptive_constraint_profile.md),
[`gompertz4_safe_constraint()`](https://immunoplex.github.io/curveRfreq/reference/gompertz4_safe_constraint.md),
[`logistic4_safe_constraint()`](https://immunoplex.github.io/curveRfreq/reference/logistic4_safe_constraint.md),
[`loglogistic4_safe_constraint()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic4_safe_constraint.md),
[`loglogistic5_safe_constraint()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic5_safe_constraint.md)
