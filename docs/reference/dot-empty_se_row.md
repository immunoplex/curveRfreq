# Create an empty summary row with NA values

Generates a single-row data.frame with grouping variables preserved and
summary statistics set to `NA` or zero. Used as a safe fallback when no
valid observations are available for a grouping.

## Usage

``` r
.empty_se_row(grouping, grouping_cols)
```

## Arguments

- grouping:

  A data.frame or named list representing grouping variables.

- grouping_cols:

  Character vector of grouping column names (unused but retained for
  consistency).

## Value

A data.frame with one row containing grouping values and NA summary
fields.
