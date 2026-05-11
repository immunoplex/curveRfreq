# Fit All Candidate Sigmoid Models Using Robust Multi-Start Strategy

Attempts to fit each of the five candidate models (logistic5,
loglogistic5, logistic4, loglogistic4, gompertz4) using
[`nlsLM`](https://rdrr.io/pkg/minpack.lm/man/nlsLM.html). Multiple start
vectors are tried; the fit with the lowest AIC is retained. For
low-signal data two fallback strategies are attempted when all primary
starts fail:

1.  Relaxed bounds (50% wider).

2.  Base-R `nls` with the `"port"` algorithm.

## Usage

``` r
compute_robust_curves(
  prepped_data,
  response_variable,
  independent_variable,
  formulas,
  model_constraints,
  start_lists,
  verbose = TRUE
)
```

## Arguments

- prepped_data:

  Data frame of preprocessed standard curve data.

- response_variable:

  Character. Response column name.

- independent_variable:

  Character. Concentration column name.

- formulas:

  Named list of model formulae.

- model_constraints:

  Named list from
  [`obtain_model_constraints`](https://immunoplex.github.io/curveRfreq/reference/obtain_model_constraints.md).

- start_lists:

  Named list from
  [`make_start_lists`](https://immunoplex.github.io/curveRfreq/reference/make_start_lists.md).

- verbose:

  Logical (default `TRUE`).

## Value

Named list (one element per model). Each element is a list with:

- fit:

  `nlsLM` object, or `NULL` if fitting failed.

- data:

  The data used for fitting.
