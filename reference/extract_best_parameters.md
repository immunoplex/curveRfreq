# Extract Parameters from the Best Fit

Pulls coefficient estimates, standard errors, and confidence intervals
from the selected best model.

## Usage

``` r
extract_best_parameters(
  ensemble,
  best_model_name,
  fixed_a = NULL,
  model_constraints = NULL,
  level = 0.95
)
```

## Arguments

- ensemble:

  Output of
  [`fit_ensemble_nls()`](https://immunoplex.github.io/curveRfreq/reference/fit_ensemble_nls.md).

- best_model_name:

  Character. Name of the selected model.

- fixed_a:

  Numeric or NULL. Fixed lower asymptote value (on the fitting scale).

- model_constraints:

  Constraints list (for bound info).

- level:

  Confidence level. Default 0.95.

## Value

Data frame with columns: term, estimate, std_error, ci_lower, ci_upper,
lower_bound, upper_bound.
