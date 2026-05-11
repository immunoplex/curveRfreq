# Apply a Blank Operation to Standard Curve Data

Performs one of five blank-handling strategies on the standard curve
data:

- `"ignored"`:

  No blank adjustment (default).

- `"included"`:

  Append blank geometric mean as an extra curve point via
  [`include_blanks_conc`](https://immunoplex.github.io/curveRfreq/reference/include_blanks_conc.md).

- `"subtracted"`:

  Subtract the geometric mean of blanks from all responses.

- `"subtracted_3x"`:

  Subtract three times the geometric mean.

- `"subtracted_10x"`:

  Subtract ten times the geometric mean.

After subtraction, values that become zero or negative are floored at 0
(linear scale) or 1 (log scale) to prevent downstream errors.

## Usage

``` r
perform_blank_operation(
  blank_data,
  data,
  response_variable,
  independent_variable,
  is_log_response,
  blank_option = "ignored",
  verbose = TRUE
)
```

## Arguments

- blank_data:

  Data frame of blank measurements.

- data:

  Data frame of standard curve data.

- response_variable:

  Character. Name of the response column.

- independent_variable:

  Character. Name of the concentration column.

- is_log_response:

  Logical. Whether the response has been log10-transformed.

- blank_option:

  Character. One of `"ignored"`, `"included"`, `"subtracted"`,
  `"subtracted_3x"`, `"subtracted_10x"`.

- verbose:

  Logical. Emit a message describing which operation was applied
  (default `TRUE`).

## Value

`data` with the blank operation applied.
