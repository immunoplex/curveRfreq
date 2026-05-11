# Compute Concentration Column from Dilution and Undiluted Standard

Converts the `dilution` column in `data` to an absolute concentration by
multiplying the reciprocal of each dilution factor by the undiluted
standard concentration. Optionally applies a log10 transform.

## Usage

``` r
compute_concentration(
  data,
  undiluted_sc_concentration,
  independent_variable,
  is_log_concentration = TRUE
)
```

## Arguments

- data:

  Data frame containing standards with a dilution column

- undiluted_sc_concentration:

  Numeric. Concentration of the undiluted standard (e.g. `10000`).

- independent_variable:

  Character. Name of the column to create or overwrite with
  concentration values.

- is_log_concentration:

  Logical. If `TRUE` (default), the computed concentration is
  log10-transformed.

## Value

`data` with the `independent_variable` column populated.
