# Log10-Transform the Assay Response Variable

Applies [`log10()`](https://rdrr.io/r/base/Log.html) to the response
column when `is_log_response` is `TRUE`. Non-positive values should be
floored *before* calling this function (see
[`preprocess_robust_curves`](https://immunoplex.github.io/curveRfreq/reference/preprocess_robust_curves.md)).

## Usage

``` r
compute_log_response(data, response_variable, is_log_response = TRUE)
```

## Arguments

- data:

  Data frame.

- response_variable:

  Character. Name of the response column.

- is_log_response:

  Logical. If `TRUE` (default), apply `log10`.

## Value

`data` with the response column optionally log-transformed.
