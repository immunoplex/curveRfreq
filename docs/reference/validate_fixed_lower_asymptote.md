# Validate a Raw Fixed Lower Asymptote Before Log Transformation

Checks that `fixed_a_result_raw` is a positive, finite scalar before it
is passed to [`log10()`](https://rdrr.io/r/base/Log.html) downstream.
Returns the original value if valid, or `NULL` to indicate the parameter
should be treated as free. Called by: predict_and_propagate_error(),
select_model_formulas() callers

Checks that `fixed_a_result_raw` is a positive, finite, scalar numeric
before it is `log10`-transformed downstream. Returns the original value
if valid, `NULL` if it should be treated as free.

## Usage

``` r
validate_fixed_lower_asymptote(fixed_a_result_raw, verbose = TRUE)

validate_fixed_lower_asymptote(fixed_a_result_raw, verbose = TRUE)
```

## Arguments

- fixed_a_result_raw:

  Raw numeric value to validate.

- verbose:

  Logical; if `TRUE` emit informational messages (default `TRUE`).

## Value

`fixed_a_result_raw` if valid; `NULL` otherwise.

The original `fixed_a_result_raw` value if valid, otherwise `NULL`.
