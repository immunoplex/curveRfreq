# Compute lower and upper limits of detection (LOD)

Calculates the lower limit of detection (LLOD) and upper limit of
detection (ULOD) based on model confidence intervals and optional blank
standard error adjustment.

## Usage

``` r
generate_lods(best_fit, fixed_a_result, std_error_blank, verbose = TRUE)
```

## Arguments

- best_fit:

  A list containing model results with at least:

  best_ci

  :   A data.frame of parameter confidence intervals with columns
      `parameter`, `conf.low`, and `conf.high`.

  best_data

  :   A data.frame of the data used to fit the model.

- fixed_a_result:

  Optional numeric value for the fixed lower asymptote ("a"). If
  provided, the LLOD is computed using this value plus a margin of
  error.

- std_error_blank:

  Optional numeric standard error of the blank. If `NULL` or `NA`, it
  defaults to 0 when `fixed_a_result` is used.

- verbose:

  Logical; if `TRUE`, prints intermediate values used in the
  calculation.

## Value

A named list with:

- llod:

  Lower limit of detection.

- ulod:

  Upper limit of detection. Returns `NA` if invalid (e.g., negative or
  less than LLOD).

## Details

The ULOD is derived from the lower confidence bound of the upper
asymptote ("d"). The LLOD is either:

- The upper confidence bound of the lower asymptote ("a"), or

- A fixed asymptote value plus a margin of error based on a
  t-distribution and the blank standard error.

Invalid ULOD values (negative or less than LLOD) are set to `NA`.

## Examples

``` r
if (FALSE) { # \dontrun{
lods <- generate_lods(best_fit, fixed_a_result = NULL, std_error_blank = NULL)
} # }
```
