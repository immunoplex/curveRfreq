# Predict Grid with Uncertainty (Frequentist)

For each point on the prediction grid: evaluates the forward model to
get predicted response and delta-method CI (via
[`curveRcore::compute_curve_ci()`](https://immunoplex.github.io/curveRcore/reference/compute_curve_ci.html)),
then back-calculates concentration and propagates uncertainty to get
`se_concentration` and `pcov`.

## Usage

``` r
predict_grid_freq(
  grid,
  fit,
  model_name,
  fixed_a = NULL,
  se_response = 0,
  cv_x_max = 150,
  pcov_threshold = 20,
  is_log_independent = TRUE,
  is_log_response = TRUE,
  independent_variable = "concentration",
  verbose = FALSE
)
```

## Arguments

- grid:

  Data frame from
  [`curveRcore::generate_prediction_grid()`](https://immunoplex.github.io/curveRcore/reference/generate_prediction_grid.html).
  Must contain an `x_fit` column (grid points on the fitting scale).

- fit:

  The best-fit `nlsLM` object.

- model_name:

  Character. Model name (must be in
  [`curveRcore::available_models()`](https://immunoplex.github.io/curveRcore/reference/available_models.html)).

- fixed_a:

  Numeric or NULL. Fixed lower asymptote on fitting scale. When non-NULL
  the `a` parameter is removed from `theta` and `Sigma` before gradient
  computation.

- se_response:

  Numeric. Residual SE of the response on the fitting scale. Typically
  `summary(fit)$sigma`. Default 0 (parameter uncertainty only).

- cv_x_max:

  Numeric. Hard cap for `pcov` and `pcov_rmse` (%). Values above this
  are clamped. Default 150.

- pcov_threshold:

  Numeric. Precision threshold (%) used to set `pcov_pass`. A grid point
  passes when `pcov < pcov_threshold`. Default 20. Must be
  `<= cv_x_max`.

- is_log_independent:

  Logical. Is the independent variable (x) on the log10 scale? When
  `TRUE`, a `log10_concentration` column is written as an alias for
  `x_fit`, and `pcov` is computed as `SE_x * log(10) * 100`. Default
  `TRUE`.

- is_log_response:

  Logical. Passed to
  [`curveRcore::enrich_grid_with_d2y()`](https://immunoplex.github.io/curveRcore/reference/enrich_grid_with_d2y.html).
  Default `TRUE`.

- independent_variable:

  Character. Name of the independent variable column in `grid`. Default
  `"concentration"`.

- verbose:

  Logical. Emit per-point diagnostic messages. Default `FALSE`.

## Value

The input `grid` data frame with the following columns appended:

- `log10_concentration`:

  Alias for `x_fit` when `is_log_independent = TRUE`; otherwise
  `NA_real_`. Present unconditionally so downstream code can always
  reference it.

- `predicted_response`:

  Forward model prediction \\\hat{y} = f(x, \hat{\theta})\\ on the
  fitting scale.

- `ci_lower`, `ci_upper`:

  Delta-method 95\\ scale.

- `predicted_concentration`:

  Back-calculated concentration from \\\hat{y}\\, on the fitting scale.

- `se_concentration`:

  Delta-method SE of `predicted_concentration`.

- `pcov`:

  Percent CV of back-calculated concentration (\\ at `cv_x_max`.

- `pcov_rmse`:

  Relative RMSE (\\ bias of `predicted_concentration` vs `x_fit`. Capped
  at `cv_x_max`.

- `pcov_pass`:

  Logical. `TRUE` when `pcov < pcov_threshold` and `pcov` is finite.

- `d2y_dx2`:

  Second derivative \\d^2(\log\_{10} y)/d(\log\_{10} x)^2\\ at each grid
  point, added by
  [`curveRcore::enrich_grid_with_d2y()`](https://immunoplex.github.io/curveRcore/reference/enrich_grid_with_d2y.html).

## Details

The `se_response` parameter controls observation-level noise in the
delta method. When non-zero, it adds the response variance term to the
concentration uncertainty via the O'Connell et al. (1993) formula:
`Var(x) = grad_theta' * Sigma * grad_theta + (dx/dy)^2 * sigma_y^2`.
