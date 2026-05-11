# Propagate uncertainty for inverse predictions across a data frame

Applies analytic error propagation to a vector of observed responses,
returning inverse-predicted concentrations and associated uncertainty.

## Usage

``` r
propagate_error_dataframe(
  pred_df,
  fit,
  model = c("logistic4", "loglogistic4", "gompertz4", "logistic5", "loglogistic5"),
  y_col,
  se_col,
  fixed_a,
  cv_x_max = 125,
  is_log_x = TRUE,
  quiet = FALSE
)
```

## Arguments

- pred_df:

  Data frame containing observed responses and their errors.

- fit:

  A fitted `nlsLM` model object.

- model:

  Character string specifying the model form. One of `"logistic4"`,
  `"loglogistic4"`, `"gompertz4"`, `"logistic5"`, `"loglogistic5"`.

- y_col:

  Character. Column name containing observed responses.

- se_col:

  Character. Column name containing standard errors of responses.

- fixed_a:

  Optional numeric scalar. Fixed lower asymptote.

- cv_x_max:

  Numeric. Maximum allowed coefficient of variation (CV). Values above
  this are capped. Default is `125`.

- is_log_x:

  Logical. Whether `x_est` is on the log10 scale. Defaults to `TRUE`.

- quiet:

  Logical. If `FALSE`, prints progress and diagnostics.

## Value

The input data frame with added columns:

- predicted_concentration:

  Inverse-predicted concentration

- se_x:

  Standard error of predicted concentration

- cv_x:

  Coefficient of variation (capped at `cv_x_max`)

## Details

Uses the delta method to propagate uncertainty from:

- Model parameter covariance (`vcov(fit)`)

- Measurement error (`se_col`)

When `is_log_x = TRUE`, CV is computed on the linear scale: \$\$CV =
se_x \cdot \log(10) \cdot 100\$\$

This avoids instability when `x_est` is near zero on the log scale.
