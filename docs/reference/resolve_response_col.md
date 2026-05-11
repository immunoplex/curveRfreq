# Resolve the assay response column name

Determines the appropriate response column name in a data frame using
metadata when available. Falls back to commonly used response columns
for backward compatibility.

## Usage

``` r
resolve_response_col(df, default = "mfi")
```

## Arguments

- df:

  A data.frame containing assay data (e.g., standards, blanks, samples).

- default:

  Character. Fallback column name (default = `"mfi"`).

## Value

A character string indicating the column name to use for response
values.

## Details

Resolution priority:

1.  `assay_response_variable` column (if present and valid)

2.  Provided `default` column (if present in data)

3.  First available among common response columns: `"mfi"`,
    `"absorbance"`, `"fluorescence"`, `"od"`

4.  Fallback to `default`
