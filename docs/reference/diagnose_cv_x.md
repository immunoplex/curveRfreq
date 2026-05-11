# Diagnose coefficient of variation (CV) behavior

Provides summary statistics and diagnostics for propagated CV values,
including detection of capped values and behavior within LOQ bounds.

## Usage

``` r
diagnose_cv_x(
  df,
  label = "pred_se",
  lloq = NULL,
  uloq = NULL,
  cv_x_max = 125,
  verbose = TRUE
)

diagnose_cv_x(
  df,
  label = "pred_se",
  lloq = NULL,
  uloq = NULL,
  cv_x_max = 125,
  verbose = TRUE
)
```

## Arguments

- df:

  data.frame containing at least `cv_x` and `predicted_concentration`
  columns.

- label:

  Character label used in messages (default `"pred_se"`).

- lloq:

  Numeric or `NULL`; lower limit of quantification.

- uloq:

  Numeric or `NULL`; upper limit of quantification.

- cv_x_max:

  Numeric cap applied to CV values (default `125`).

- verbose:

  Logical; if `FALSE` returns invisibly without printing (default
  `TRUE`).

## Value

Invisibly returns a list with summary statistics:

- min_cv:

  Minimum CV

- min_x:

  Concentration at minimum CV

- max_cv:

  Maximum CV

- mean_cv:

  Mean CV

- n_gt_20:

  Number of CV values \> 20

- n_at_cap:

  Number of capped CV values

A named list with summary statistics, returned invisibly.

## Details

Helps identify:

- Instability near asymptotes

- Excessive uncertainty within LOQ range

- Propagation failures leading to capped values
