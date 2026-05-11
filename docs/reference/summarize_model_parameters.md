# Summarise Parameter Estimates and Confidence Intervals for All Models

Returns a tidy data frame with one row per parameter per model.
Confidence intervals are computed using
[`nlstools::confint2`](https://rdrr.io/pkg/nlstools/man/confint2.html);
if that fails, `NA` is reported.

## Usage

``` r
summarize_model_parameters(
  models_fit_list,
  level = 0.95,
  model_names = c("logistic5", "loglogistic5", "logistic4", "loglogistic4", "gompertz4"),
  verbose = TRUE
)
```

## Arguments

- models_fit_list:

  Named list from
  [`compute_robust_curves`](https://immunoplex.github.io/curveRfreq/reference/compute_robust_curves.md).

- level:

  Numeric. Confidence level (default `0.95`).

- model_names:

  Character vector of model names to include.

- verbose:

  Logical (default `TRUE`).

## Value

Data frame with columns: `model`, `parameter`, `estimate`, `conf.low`,
`conf.high`, `converged`.
