# Extract a Combined Summary Table from Multi-Curve Results

One row per curve_id with best model, parameters, and fit statistics.

## Usage

``` r
summary_table(mp)
```

## Arguments

- mp:

  A `calibration_result_multiplate` object.

## Value

Data frame with columns: `curve_id`, `best_model`, `converged`, `aic`,
`bic`, `rss`, `n_obs`, `n_samples`, `a`, `b`, `c`, `d`, `g`.
