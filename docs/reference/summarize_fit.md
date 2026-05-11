# Summarize Fitted Model with QC Metrics, LOD, Detection Limits, and Curvature LOQs

Generates a one-row summary of a fitted model, including parameter
estimates, goodness-of-fit statistics, inflection point coordinates,
limits of detection (LOD), minimum/maximum detectable concentrations
(MDC), reliable detection limits (RDL), curvature-based limits of
quantification (LOQ), and metadata. The result is stored in
`best_fit$best_fit_summary`.

## Usage

``` r
summarize_fit(
  best_fit,
  response_variable,
  independent_variable,
  fixed_a_result,
  antigen_settings,
  curve_id_lookup,
  antigen_fit_options,
  verbose = TRUE
)
```

## Arguments

- best_fit:

  A list containing model results, including:

  - `best_fit$best_fit`: Fitted model object

  - `best_fit$best_data`: Data used for fitting

  - `best_fit$best_model_name`: the model name of the best fit (selected
    by minimizing the AIC score)

  - `best_fit$best_d2xy`: data.frame of second derivative values

- response_variable:

  Character string naming the response variable in `best_data` (e.g.
  mfi, absorbance).

- independent_variable:

  Character string naming the predictor variable used in the model.

- fixed_a_result:

  numeric value for parameter `a`. This a derived result from the
  [`select_antigen_plate`](https://immunoplex.github.io/curveRfreq/reference/select_antigen_plate.md)
  function.

- antigen_settings:

  List of antigen-specific settings. May include:

  - `std_error_blank`: Standard error of the blank used for LLOD
    calculation

- curve_id_lookup:

  lookup table for the selected curve_id associating with its
  components.

- antigen_fit_options:

  List of model fitting options, including:

  - `blank_option`

  - `is_log_response`

  - `is_log_concentration`

  - `apply_prozone`

- verbose:

  Logical; if `TRUE`, prints progress and diagnostic messages.

## Value

A list identical to `best_fit` with additional elements:

- `best_fit_summary`: A one-row `data.frame` containing model summary
  and QC metrics

- `lods`: A list with `llod` and `ulod`

- `mdc_rdl`: A list with `mindc`, `maxdc`, `minrdl`, and `maxrdl`

- `curv_loqs`: A list with `lloq`, `uloq`, `lloq_y`, and `uloq_y`

## Details

This function is designed to operate on the output of a model fitting
step and produce a consistent summary structure for downstream analysis,
QC filtering, or reporting.

The returned summary includes:

- Model parameters: `a, b, c, d, g`

- Inflection point coordinates: `inflect_x`, `inflect_y`

- Goodness-of-fit metrics: residual sum of squares (RSS), mean squared
  error (MSE), R-squared, Akaike information criterion (AIC), Bayesian
  information criterion (BIC), and log-likelihood

- Convergence diagnostics: iteration count and convergence status

- Limits of detection:

  - `llod`: Lower limit of detection

  - `ulod`: Upper limit of detection

- Detection limits derived from the fitted curve:

  - `mindc`, `maxdc`: Minimum and maximum detectable concentrations

  - `minrdl`, `maxrdl`: Reliable detection limits based on confidence
    intervals

- Curvature-based limits of quantification (LOQ):

  - `lloq`, `uloq`: LOQ values based on extrema of the second derivative

  - `lloq_y`, `uloq_y`: Corresponding predicted response values

- Metadata: source, transformation flags, and model formula

The inflection point is computed using
[`compute_inflection_point`](https://immunoplex.github.io/curveRfreq/reference/compute_inflection_point.md)
and corresponds to the location where the second derivative of the
fitted curve equals zero.

Limits of detection (LOD) are computed using
[`generate_lods`](https://immunoplex.github.io/curveRfreq/reference/generate_lods.md),
optionally incorporating blank standard error. MDC and RDL values are
computed using
[`generate_mdc_rdl`](https://immunoplex.github.io/curveRfreq/reference/generate_mdc_rdl.md)
via root-finding over the observed data range.

Curvature-based LOQs are computed using
[`compute_loqs`](https://immunoplex.github.io/curveRfreq/reference/compute_loqs.md),
which identifies local extrema in the second derivative of the fitted
curve and refines their locations using quadratic interpolation.

If no valid model, derivative data, or extrema are available,
corresponding values are returned as `NA`.

## See also

[`compute_inflection_point`](https://immunoplex.github.io/curveRfreq/reference/compute_inflection_point.md),
[`generate_lods`](https://immunoplex.github.io/curveRfreq/reference/generate_lods.md),
[`generate_mdc_rdl`](https://immunoplex.github.io/curveRfreq/reference/generate_mdc_rdl.md),
[`compute_loqs`](https://immunoplex.github.io/curveRfreq/reference/compute_loqs.md)

## Examples

``` r
if (FALSE) { # \dontrun{
best_fit <- fit_model(data)

best_fit <- summarize_fit(
  best_fit = best_fit,
  response_variable = "mfi",
  independent_variable = "concentration",
  fixed_a_result = NULL,
  antigen_settings = list(std_error_blank = 0.05),
  antigen_fit_options = list(
    blank_option = "none",
    is_log_response = FALSE,
    is_log_concentration = TRUE,
    apply_prozone = FALSE
  )
)

best_fit$best_fit_summary
best_fit$lods
best_fit$mdc_rdl
best_fit$curv_loqs
} # }
```
