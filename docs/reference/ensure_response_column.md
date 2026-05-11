# Ensure the response variable column exists in a data frame. If the named column is missing, attempts to find it via assay_response_variable metadata or common response column names. Optionally coerces to numeric.

Ensure the response variable column exists in a data frame. If the named
column is missing, attempts to find it via assay_response_variable
metadata or common response column names. Optionally coerces to numeric.

## Usage

``` r
ensure_response_column(df, response_var, coerce_numeric = TRUE, context = "")
```

## Arguments

- df:

  Data frame to check

- response_var:

  Expected column name (e.g. "mfi", "absorbance")

- coerce_numeric:

  Logical; if TRUE, coerce the column to numeric

- context:

  Character label for diagnostic messages

## Value

A list with:

- df:

  The (possibly modified) data frame

- response_var:

  The resolved column name (may differ from input)

- ok:

  Logical: TRUE if a valid numeric response column was found
