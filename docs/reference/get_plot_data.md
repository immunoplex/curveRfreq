# Build prediction and diagnostic data for model comparison plots

For each model in `model_names`, generates predicted values on a fine
grid, residuals, first and second derivatives, and delta-method
confidence intervals. Appends AIC values to the parameter summary for
downstream plotting.

## Usage

``` r
get_plot_data(
  models_fit_list,
  prepped_data,
  fit_params,
  fixed_a_result,
  curve_id_lookup,
  model_names = c("logistic5", "loglogistic5", "logistic4", "loglogistic4", "gompertz4"),
  x_var = "concentration",
  y_var = "mfi",
  verbose = TRUE
)
```

## Arguments

- models_fit_list:

  Named list of model fit objects (or lists with a `$fit` element).
  Names should correspond to `model_names`.

- prepped_data:

  A `data.frame` containing at minimum columns named by `x_var` and
  `y_var`.

- fit_params:

  A `data.frame` of parameter estimates with columns `model`,
  `parameter`, `estimate`, `conf.low`, and `conf.high`, as returned by
  [`summarize_model_fits()`](https://immunoplex.github.io/curveRfreq/reference/summarize_model_fits.md).

- fixed_a_result:

  Numeric scalar or `NULL`. Fixed lower asymptote passed to
  [`compute_curve_ci()`](https://immunoplex.github.io/curveRfreq/reference/compute_curve_ci.md)
  and used to reconstruct full coefficient vectors.

- curve_id_lookup:

  the order of the elements of the curve_id.

- model_names:

  Character vector of model names to process. Must match keys in
  `models_fit_list`. Default is
  `c("logistic5","loglogistic5","logistic4","loglogistic4","gompertz4")`.

- x_var:

  Character string naming the independent variable column in
  `prepped_data`. Default is `"concentration"`.

- y_var:

  Character string naming the response variable column in
  `prepped_data`. Default is `"mfi"`.

- verbose:

  Logical. If `TRUE`, prints a completion message. Default is `TRUE`.

## Value

A named list with elements:

- dat:

  `data.frame` of the input data with added `x` and `y` columns.

- pred_df:

  `data.frame` of predicted values on the fine grid, with columns
  `model`, `x`, `yhat`.

- resid_df:

  `data.frame` of residuals, with columns `model`, `fitted`,
  `residuals`.

- fit_params:

  `data.frame` of parameter estimates as supplied in `fit_params`.

- d2xy_df:

  `data.frame` of second-derivative values, with columns `model`, `x`,
  `d2x_y`.

- dydx_df:

  `data.frame` of first-derivative values, with columns `model`, `x`,
  `dydx`.

- ci_df:

  `data.frame` of confidence intervals, with columns `model`, `x`,
  `ci_lo`, `ci_hi`.

- fit_params_aic:

  `data.frame` combining `fit_params` with AIC values for each converged
  model.

## See also

[`compute_curve_ci`](https://immunoplex.github.io/curveRfreq/reference/compute_curve_ci.md),
[`plot_model_comparisons`](https://immunoplex.github.io/curveRfreq/reference/plot_model_comparisons.md)

## Examples

``` r
if (FALSE) { # \dontrun{
plot_data <- get_plot_data(
  models_fit_list = my_fits,
  prepped_data    = my_data,
  fit_params      = my_params,
  fixed_a_result  = NULL
)
} # }
```
