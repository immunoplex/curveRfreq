# Plot comparison of nonlinear model fits

Produces a three-panel patchwork figure (or a list of individual ggplot2
objects) showing (A) observed data with fitted curves, (B) residuals vs
fitted values, and (C) parameter estimates with confidence intervals
faceted by parameter name.

## Usage

``` r
plot_model_comparisons(
  plot_data,
  model_names = c("logistic5", "loglogistic5", "logistic4", "loglogistic4", "gompertz4"),
  x_var = "concentration",
  y_var = "mfi",
  is_display_log_response = TRUE,
  is_display_log_independent = TRUE,
  use_patchwork = TRUE
)
```

## Arguments

- plot_data:

  Named list as returned by
  [`get_plot_data`](https://immunoplex.github.io/curveRfreq/reference/get_plot_data.md).

- model_names:

  Character vector of model names to include in the plots. Default is
  `c("logistic5","loglogistic5","logistic4","loglogistic4","gompertz4")`.

- x_var:

  Character string naming the independent variable. Used for axis
  labelling after passing through
  [`format_assay_terms`](https://immunoplex.github.io/curveRfreq/reference/format_assay_terms.md).
  Default is `"concentration"`.

- y_var:

  Character string naming the response variable. Used for axis
  labelling. Default is `"mfi"`.

- is_display_log_response:

  Logical. If `TRUE`, the y-axis label is formatted as
  \\\log\_{10}(\text{y\\var})\\. Default is `TRUE`.

- is_display_log_independent:

  Logical. If `TRUE`, the x-axis label is formatted as
  \\\log\_{10}(\text{x\\var})\\. Default is `TRUE`.

- use_patchwork:

  Logical. If `TRUE` (default), returns a single combined patchwork
  figure. If `FALSE`, returns a named list of individual ggplot2 objects
  (`data_fit`, `resid`, `p_ci_params`).

## Value

When `use_patchwork = TRUE`, a patchwork object combining all three
panels. When `use_patchwork = FALSE`, a named list with elements
`data_fit`, `resid`, and `p_ci_params`.

## See also

[`get_plot_data`](https://immunoplex.github.io/curveRfreq/reference/get_plot_data.md),
[`format_assay_terms`](https://immunoplex.github.io/curveRfreq/reference/format_assay_terms.md)

## Examples

``` r
if (FALSE) { # \dontrun{
plot_model_comparisons(plot_data = my_plot_data,
                       is_display_log_response    = TRUE,
                       is_display_log_independent = TRUE)
} # }
```
