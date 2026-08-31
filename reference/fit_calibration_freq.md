# Fit a Frequentist Calibration Curve (Single Curve)

Fits an ensemble of nonlinear models to preprocessed standard curve data
for one curve_id, selects the best model by AIC subject to
quantification-eligibility gates, and produces grid and sample
predictions with propagated uncertainty.

## Usage

``` r
fit_calibration_freq(
  standards,
  blanks = NULL,
  samples = NULL,
  response_var,
  model_names = c("logistic4", "gompertz4"),
  is_log_response = TRUE,
  is_log_independent = TRUE,
  std_curve_conc,
  fixed_a = NULL,
  cv_x_max = 150,
  pcov_threshold = 20,
  min_dynamic_range_log10 = 0.5,
  max_rel_se = 5,
  bound_tol = 1e-04,
  n_grid = 200L,
  grid_min_conc = 1e-04,
  grid_max_conc = NULL,
  curve_id = "1",
  persist_draws = FALSE,
  verbose = FALSE
)
```

## Arguments

- standards:

  Data frame. Preprocessed standard curve data with a response column
  (named by `response_var`) and a `concentration` column, both already
  on the fitting scale.

- blanks:

  Data frame or NULL. Preprocessed blank data with a response column
  (named by `response_var`), on the fitting scale. Stored in
  `result$blanks` for QA and plotting only — not used in fitting. NULL
  (default) leaves `result$blanks` empty.

- samples:

  Data frame or NULL. Test sample data with the response column (on the
  **raw** measurement scale — log-transform is applied internally for
  back-calculation). Optionally a `dilution` column.

- response_var:

  Character. Name of the response column.

- model_names:

  Character vector. Models to fit. Must be a subset of
  [`curveRcore::available_models()`](https://immunoplex.github.io/curveRcore/reference/available_models.html).
  Default `c("logistic4", "gompertz4")`.

- is_log_response:

  Logical. Was the response log10-transformed during preprocessing?
  Needed to correctly back-transform sample predictions. Default TRUE.

- is_log_independent:

  Logical. Was concentration log10-transformed? Needed for pcov
  computation and grid generation. Default TRUE.

- std_curve_conc:

  Numeric. Undiluted standard concentration (raw scale). Used for grid
  generation.

- fixed_a:

  Numeric or NULL. Fixed lower asymptote on the **fitting** scale (i.e.
  already log-transformed if `is_log_response = TRUE`). NULL means `a`
  is estimated freely.

- cv_x_max:

  Numeric. Cap for percent CV of back-calculated concentration. Default
  150.

- pcov_threshold:

  Numeric. Percent CV threshold for LLOQ/ULOQ determination and the
  dynamic-range eligibility gate. Default 20.

- min_dynamic_range_log10:

  Numeric. Minimum dynamic range (log10 units) required for a model to
  be eligible for selection. Default 0.5.

- max_rel_se:

  Numeric. Maximum relative SE (SE/\|estimate\|) permitted for any
  parameter; models exceeding this are ineligible. Default 5.0.

- bound_tol:

  Numeric. Tolerance for the at-bound gate. Default 1e-4.

- n_grid:

  Integer. Number of prediction grid points. Default 200.

- grid_min_conc:

  Numeric. Minimum grid concentration (raw scale). Default 1e-4.

- grid_max_conc:

  Numeric or NULL. Maximum grid concentration. NULL uses
  `std_curve_conc`.

- curve_id:

  Character or numeric. Identifier for this curve. Default `"1"`.

- verbose:

  Logical. Default FALSE.

## Value

A `calibration_result` object (from curveRcore). The `$selection`
component now contains `$assessments` (per-model eligibility results),
`$eligible_models`, and `$fallback`.

## Details

**Important:** `standards` and 'blanks' must already be on the fitting
scale. Both the response column and the `concentration` column should be
transformed (e.g. log10) before calling this function. Use
[`curveRcore::preprocess_standards()`](https://immunoplex.github.io/curveRcore/reference/preprocess_standards.html)
upstream.
