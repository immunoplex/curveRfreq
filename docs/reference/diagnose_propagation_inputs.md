# Diagnose inputs for analytic error propagation

Prints diagnostic information to help debug issues in error propagation,
including parameter alignment, covariance structure, and gradient
evaluation.

Emits a structured diagnostic block summarising the model coefficients,
variance-covariance matrix, and optionally the inverse/gradient
functions evaluated at a test response value. Call this immediately
before
[`propagate_error_dataframe`](https://immunoplex.github.io/curveRfreq/reference/propagate_error_dataframe.md).

## Usage

``` r
diagnose_propagation_inputs(fit, model, fixed_a, y_test = NULL)

diagnose_propagation_inputs(fit, model, fixed_a, y_test = NULL)
```

## Arguments

- fit:

  Fitted model object with [`coef()`](https://rdrr.io/r/stats/coef.html)
  and [`vcov()`](https://rdrr.io/r/stats/vcov.html) methods.

- model:

  Character; model name passed to
  [`make_inv_and_grad_fixed()`](https://immunoplex.github.io/curveRfreq/reference/make_inv_and_grad_fixed.md).

- fixed_a:

  Numeric scalar or `NULL`.

- y_test:

  Optional numeric scalar test response value.

## Value

Invisibly returns `NULL`. Outputs diagnostic messages to console.

`NULL` invisibly.

## Details

This function is useful for identifying:

- Mismatches between gradient names and covariance matrix

- Incorrect handling of fixed vs free parameters

- Failures in inverse or gradient evaluation
