# Plot Standard Curve (Minimal + Uncertainty)

Generates an interactive Plotly visualization of a fitted standard
curve, including standards, samples, fitted curve, confidence intervals,
and measurement uncertainty (pCoV).

## Usage

``` r
plot_standard_curve(
  best_fit,
  is_display_log_response,
  is_display_log_independent,
  pcov_threshold,
  study_params,
  curve_id_lookup,
  response_variable = "mfi",
  independent_variable = "concentration"
)
```

## Arguments

- best_fit:

  A list containing model outputs including best_data, best_pred,
  best_glance, and optional best_curve_ci.

- is_display_log_response:

  Logical; whether to display response on log10 scale.

- is_display_log_independent:

  Logical; whether to display independent variable on log10 scale.

- pcov_threshold:

  Numeric threshold for percent coefficient of variation.

- study_params:

  study parameters including is_log_response and is_log_independent

- curve_id_lookup:

  lookup for the specific curve_id. Used to get the plate and antigen
  names

- response_variable:

  Character; name of response variable column (default "mfi").

- independent_variable:

  Character; name of independent variable column (default
  "concentration").

## Value

A plotly object.
