# Predict Concentrations and Propagate Uncertainty for a Single Plate

Given a fitted standard-curve model (`best_fit`), this function performs
two related tasks:

1.  Back-calculates concentrations for every point on the standards
    prediction grid (`best_fit$best_pred`).

2.  Back-calculates concentrations for every sample well in
    `antigen_plate$plate_samples` and attaches propagated uncertainty
    estimates.

Uncertainty is propagated via the delta method, combining parameter
covariance from `vcov(fit)` with assay measurement error supplied
through `se_std_response`.

## Usage

``` r
predict_and_propagate_error(
  best_fit,
  response_var,
  antigen_plate,
  study_params,
  se_std_response,
  cv_x_max = 150,
  verbose = TRUE
)
```

## Arguments

- best_fit:

  Named list returned by the model-selection step. Required elements:

  `best_fit`

  :   The fitted model object (e.g. from
      [`minpack.lm::nlsLM()`](https://rdrr.io/pkg/minpack.lm/man/nlsLM.html))
      with [`coef()`](https://rdrr.io/r/stats/coef.html) and
      [`vcov()`](https://rdrr.io/r/stats/vcov.html) methods.

  `best_pred`

  :   data.frame of standards prediction grid; must contain a column
      named `"yhat"`.

  `best_data`

  :   data.frame of the standards data used to fit the model; must
      contain columns `study_accession`, `experiment_accession`,
      `nominal_sample_dilution`, `plateid`, `plate`, `antigen`, and
      `source`.

  `best_model_name`

  :   Character scalar; model identifier passed to
      [`propagate_error_dataframe()`](https://immunoplex.github.io/curveRfreq/reference/propagate_error_dataframe.md)
      and
      [`diagnose_propagation_inputs()`](https://immunoplex.github.io/curveRfreq/reference/diagnose_propagation_inputs.md).

  `best_glance`

  :   List or data.frame with optional elements `lloq` and `uloq` (used
      only for diagnostics).

- response_var:

  Character scalar; name of the response column in
  `antigen_plate$plate_samples` and `antigen_plate$plate_standard` (e.g.
  `"mfi"`).

- antigen_plate:

  Named list containing plate-level data. Required elements:

  `plate_samples`

  :   data.frame of sample wells; must contain columns `response_var`,
      `dilution`, and `well`.

  `plate_standard`

  :   data.frame of standard-curve wells; must contain column
      `response_var`. Used to estimate the reference MFI for SE scale
      conversion.

  `fixed_a_result`

  :   Raw (untransformed) numeric value for the fixed lower asymptote,
      or `NULL`. Validated and log10-transformed internally when
      `study_params$is_log_response` is `TRUE`.

  `antigen_settings`

  :   List with element `standard_curve_concentration`; the maximum
      standard concentration used to cap infinite predicted values.

- study_params:

  Named list of study-level modelling flags. Required elements:

  `is_log_response`

  :   Logical; `TRUE` when the model was fitted on `log10(response)`.

  `is_log_independent`

  :   Logical; `TRUE` when the concentration axis is on the `log10`
      scale. Passed as `is_log_x` to
      [`propagate_error_dataframe()`](https://immunoplex.github.io/curveRfreq/reference/propagate_error_dataframe.md)
      and used when computing `final_predicted_concentration`.

- se_std_response:

  Numeric scalar; pooled standard error of the assay response on the
  *raw* response scale (e.g. MFI units), typically obtained from
  [`compute_antigen_se_table`](https://immunoplex.github.io/curveRfreq/reference/compute_antigen_se_table.md)
  via
  [`lookup_antigen_se`](https://immunoplex.github.io/curveRfreq/reference/lookup_antigen_se.md).
  Supply `NA` or a non-positive value to trigger the fallback SE
  calculation.

- cv_x_max:

  Numeric scalar; upper cap applied to all \\ values (default `150`).
  Passed through to
  [`propagate_error_dataframe()`](https://immunoplex.github.io/curveRfreq/reference/propagate_error_dataframe.md).

- verbose:

  Logical; if `TRUE` (default) emit step-by-step diagnostic messages and
  print the head of `antigen_plate$plate_samples`.

## Value

The input `best_fit` list with two elements added or replaced:

- `best_pred`:

  The standards prediction data.frame augmented with
  `predicted_concentration`, `se_x`, `cv_x`, `pcov`, and plate/study
  metadata columns. Infinite predicted concentrations are replaced by
  `max_conc_standard` (positive infinity) or `0` (negative infinity).

- `sample_se`:

  data.frame with one row per sample well containing
  `final_predicted_concentration`, `raw_predicted_concentration`,
  `raw_assay_response`, `se_concentration`, `cv_x`, `pcov`, and all
  original columns from `antigen_plate$plate_samples` (except
  `response_var` and `dilution`, which are carried inside `sample_se`).
  Returns an empty data.frame if the `fixed_a`/`coef(fit)` consistency
  check fails and is not correctable.

## Details

When `study_params$is_log_response` is `TRUE` all response values are
assumed to be on the `log10` scale inside the model. `se_std_response`
must be supplied on the *raw* response scale (e.g. MFI units); the
function converts it to the log10 scale internally using the
delta-method approximation \\\mathrm{se}\_{\log\_{10}} =
\sqrt{\mathrm{se}\_{\mathrm{raw}} / (\bar{z}\_{\mathrm{ref}} \cdot
\ln(10) \cdot 100)}\\, where \\\bar{z}\_{\mathrm{ref}}\\ is the
geometric mean of the raw standard-curve responses on the plate.

## See also

[`compute_antigen_se_table`](https://immunoplex.github.io/curveRfreq/reference/compute_antigen_se_table.md),
[`lookup_antigen_se`](https://immunoplex.github.io/curveRfreq/reference/lookup_antigen_se.md),
[`propagate_error_dataframe`](https://immunoplex.github.io/curveRfreq/reference/propagate_error_dataframe.md),
[`validate_fixed_lower_asymptote`](https://immunoplex.github.io/curveRfreq/reference/validate_fixed_lower_asymptote.md),
[`check_fixed_a_fit_consistency`](https://immunoplex.github.io/curveRfreq/reference/check_fixed_a_fit_consistency.md),
[`diagnose_propagation_inputs`](https://immunoplex.github.io/curveRfreq/reference/diagnose_propagation_inputs.md),
[`diagnose_cv_x`](https://immunoplex.github.io/curveRfreq/reference/diagnose_cv_x.md)
