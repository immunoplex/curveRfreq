# Check Consistency Between fixed_a_result and coef(fit)

Validates that `fixed_a_result` and the free parameters in `coef(fit)`
are mutually consistent: exactly one of (a) *a* is fixed and absent from
`coef(fit)`, or (b) *a* is free and present in `coef(fit)`.

Call this just before
[`propagate_error_dataframe`](https://immunoplex.github.io/curveRfreq/reference/propagate_error_dataframe.md)
to catch mismatches early.

## Usage

``` r
check_fixed_a_fit_consistency(
  fit,
  fixed_a_result,
  context = "",
  verbose = TRUE
)

check_fixed_a_fit_consistency(
  fit,
  fixed_a_result,
  context = "",
  verbose = TRUE
)
```

## Arguments

- fit:

  Fitted model object with a
  [`coef()`](https://rdrr.io/r/stats/coef.html) method.

- fixed_a_result:

  Numeric scalar or `NULL`.

- context:

  Character string describing the call site.

- verbose:

  Logical; if `TRUE` emit informational messages (default `TRUE`).

## Value

Named list with:

- fixed_a_result:

  The (possibly corrected) fixed_a value.

- consistent:

  Logical.

- correctable:

  Logical. If `FALSE` propagation cannot proceed.

A named list with elements `fixed_a_result`, `consistent`, and
`correctable`.
