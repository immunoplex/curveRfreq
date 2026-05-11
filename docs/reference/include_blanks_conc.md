# Include Blank Controls as an Extra Point on the Standard Curve

Appends a synthetic data row whose response is the geometric mean of the
blank controls and whose concentration is half of the minimum standard
concentration. This anchors the lower end of the fitted curve to the
background signal level.

## Usage

``` r
include_blanks_conc(
  blank_data,
  data,
  response_variable,
  independent_variable = "concentration"
)
```

## Arguments

- blank_data:

  Data frame of blank well measurements.

- data:

  Data frame of standard curve measurements.

- response_variable:

  Character. Name of the response column.

- independent_variable:

  Character. Name of the concentration column (default
  `"concentration"`).

## Value

`data` with one additional row representing the blank mean.
