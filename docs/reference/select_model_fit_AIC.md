# Select the Best Model Fit by AIC

Identifies the model with the lowest AIC from a summary table and
returns the corresponding fit object and associated predictions.

## Usage

``` r
select_model_fit_AIC(
  fit_summary,
  fit_robust_lm,
  fit_params,
  plot_data,
  verbose = TRUE
)
```

## Arguments

- fit_summary:

  Data frame from
  [`summarize_model_fits`](https://immunoplex.github.io/curveRfreq/reference/summarize_model_fits.md).

- fit_robust_lm:

  Named list from
  [`compute_robust_curves`](https://immunoplex.github.io/curveRfreq/reference/compute_robust_curves.md).

- fit_params:

  Data frame from
  [`summarize_model_parameters`](https://immunoplex.github.io/curveRfreq/reference/summarize_model_parameters.md).

- plot_data:

  Named list containing `pred_df`, `d2xy_df`, `dydx_df`, and optionally
  `ci_df`.

- verbose:

  Logical (default `TRUE`).

## Value

Named list with:

- best_model_name:

  Character. Name of the selected model.

- best_fit:

  nlsLM fit object.

- best_data:

  Data used to fit the selected model.

- best_ci:

  Parameter confidence-interval data frame.

- best_pred:

  Prediction data frame.

- best_d2xy:

  Second-derivative data frame.

- best_dydx:

  First-derivative data frame.

- best_curve_ci:

  Curve confidence-interval data frame, or `NULL`.
