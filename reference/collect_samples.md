# Extract All Sample Predictions from Multi-Curve Results

Concatenates sample prediction data frames from every curve, adding a
`curve_id` column for identification.

## Usage

``` r
collect_samples(mp)
```

## Arguments

- mp:

  A `calibration_result_multiplate` object.

## Value

Data frame of all sample predictions, or NULL if none.
