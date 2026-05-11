# Summarise Model Convergence and Information Criteria

Returns a data frame with one row per model, reporting convergence
status, RSS, degrees of freedom, number of parameters, AIC, and BIC.
Models that failed to converge appear with `converged = FALSE` and `NA`
statistics.

## Usage

``` r
summarize_model_fits(
  models_fit_list,
  model_names = c("logistic5", "loglogistic5", "logistic4", "loglogistic4", "gompertz4"),
  verbose = TRUE
)
```

## Arguments

- models_fit_list:

  Named list from
  [`compute_robust_curves`](https://immunoplex.github.io/curveRfreq/reference/compute_robust_curves.md).

- model_names:

  Character vector of model names to include (default:
  `c("logistic5","loglogistic5","logistic4","loglogistic4","gompertz4")`).

- verbose:

  Logical (default `TRUE`).

## Value

A data frame with columns: `model`, `converged`, `rss`, `df_resid`,
`n_params`, `AIC`, `BIC`.
