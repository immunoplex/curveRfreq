# Fit a Single NLS Model Using minpack.lm

Thin wrapper around
[`nlsLM`](https://rdrr.io/pkg/minpack.lm/man/nlsLM.html) with diagnostic
output and graceful error handling.

## Usage

``` r
nlsLM_fit(
  formula,
  data,
  start_values,
  lower = -Inf,
  upper = Inf,
  verbose = TRUE
)
```

## Arguments

- formula:

  A `formula` object.

- data:

  Data frame.

- start_values:

  Named list of starting parameter values.

- lower:

  Named numeric vector of lower bounds (default `-Inf`).

- upper:

  Named numeric vector of upper bounds (default `Inf`).

- verbose:

  Logical (default `TRUE`).

## Value

An `nlsLM` object, or `NULL` on failure.
