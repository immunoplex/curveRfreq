# Fit Frequentist Calibration Curves Across Multiple Curves

Splits preprocessed `standards`, `blanks`, and optionally `samples` by
`curve_id`, calls
[`fit_calibration_freq()`](https://immunoplex.github.io/curveRfreq/reference/fit_calibration_freq.md)
for each curve, and collects the results.

## Usage

``` r
fit_calibration_freq_multiplate(
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
  n_grid = 200L,
  grid_min_conc = 1e-04,
  grid_max_conc = NULL,
  curve_ids = NULL,
  on_error = c("warn", "stop"),
  persist_draws = FALSE,
  verbose = FALSE
)
```

## Arguments

- standards:

  Data frame. Stacked preprocessed standard curve data. Must contain a
  `curve_id` column, a response column (named by `response_var`), and a
  `concentration` column — all already on the fitting scale.

- blanks:

  Data frame or NULL. Stacked preprocessed blank data. When non-NULL,
  must contain a `curve_id` column and a response column (named by
  `response_var`), both already on the fitting scale. Stored in each
  `result$blanks` for QA and plotting only — not used in fitting. NULL
  (default) leaves `result$blanks` empty for every curve.

- samples:

  Data frame or NULL. Stacked sample data with a `curve_id` column and
  the response column (raw measurement scale).

- response_var:

  Character. Name of the response column.

- model_names:

  Character vector. Models to fit. Default
  `c("logistic4", "gompertz4")`.

- is_log_response:

  Logical. Was the response log10-transformed? Default TRUE.

- is_log_independent:

  Logical. Was concentration log10-transformed? Default TRUE.

- std_curve_conc:

  Numeric. Undiluted standard concentration (raw scale). Used for grid
  generation.

- fixed_a:

  Numeric or NULL. Fixed lower asymptote on the **fitting** scale.

- cv_x_max:

  Numeric. Default 150.

- n_grid:

  Integer. Default 200.

- grid_min_conc:

  Numeric. Default 1e-4.

- grid_max_conc:

  Numeric or NULL.

- curve_ids:

  Optional vector. If supplied, only these curve_ids are fitted. Default
  NULL fits all found in `standards`.

- on_error:

  Character. `"warn"` (default) stores NULL and continues; `"stop"`
  raises on first failure.

- verbose:

  Logical.

## Value

A `calibration_result_multiplate` object (from curveRcore).

## Details

**Important:** `standards` and `blanks` must already be on the fitting
scale. Use
[`curveRcore::preprocess_standards()`](https://immunoplex.github.io/curveRcore/reference/preprocess_standards.html)
upstream on each curve's data before stacking, or preprocess the full
stacked frame if all curves share the same settings.
