# StandardCurve

Wraps the full `curveRfreq` standard-curve fitting pipeline into a
single stateful object. Data loading is external — supply a pre-loaded
data list via `$new()` or `$set_data()`, then call
`$set_curve_settings()`, `$fit()`, `$summarize()`, `$plot()`, and
`$propagate_error()` in sequence. All intermediate results are stored as
public fields and can be inspected at any stage.

## Public fields

- `loaded_data`:

  List of pre-loaded assay data (standards, blanks, samples,
  antigen_constraints, response_var, indep_var, ...). Produced
  externally, e.g. by `pull_data()`.

- `response_var`:

  Name of the response column (e.g. `"mfi"`).

- `indep_var`:

  Name of the independent variable column (e.g. `"concentration"`).

- `study_params`:

  Named list of study-level modelling parameters. Keys: `applyProzone`,
  `blank_option`, `is_log_response`, `is_log_independent`.

- `antigen_constraints`:

  Data frame with per-antigen constraint settings.

- `model_names`:

  Character vector of candidate model names to consider.

- `curve_col`:

  Name of the curve-ID column in the data. Default `"curve_id"`.

- `curve_id_element_order`:

  Character vector defining the positional order of fields inside a
  `curve_id` string.

- `is_display_log_response`:

  Logical. Display response on log scale in plots.

- `is_display_log_independent`:

  Logical. Display concentration on log scale in plots.

- `verbose`:

  Logical. Print diagnostic messages when `TRUE`.

- `antigen_plate`:

  List returned by
  [`select_antigen_plate()`](https://immunoplex.github.io/curveRfreq/reference/select_antigen_plate.md).

- `processed_data`:

  List returned by
  [`preprocess_robust_curves()`](https://immunoplex.github.io/curveRfreq/reference/preprocess_robust_curves.md).

- `formulas`:

  List of model formulas from
  [`select_model_formulas()`](https://immunoplex.github.io/curveRfreq/reference/select_model_formulas.md).

- `model_constraints`:

  Parameter bounds from
  [`obtain_model_constraints()`](https://immunoplex.github.io/curveRfreq/reference/obtain_model_constraints.md).

- `start_lists`:

  Multi-start parameter lists from
  [`make_start_lists()`](https://immunoplex.github.io/curveRfreq/reference/make_start_lists.md).

- `fit_robust_lm`:

  Raw list of candidate model fits.

- `fit_summary`:

  Data frame summarising each candidate fit (AIC, BIC, ...).

- `fit_params`:

  Data frame of per-model parameter estimates and CIs.

- `plot_data`:

  List used internally for plotting / model comparison.

- `best_fit`:

  List returned by
  [`select_model_fit_AIC()`](https://immunoplex.github.io/curveRfreq/reference/select_model_fit_AIC.md)
  and augmented by
  [`tidy.nlsLM()`](https://immunoplex.github.io/curveRfreq/reference/tidy.nlsLM.md)
  and
  [`summarize_fit()`](https://immunoplex.github.io/curveRfreq/reference/summarize_fit.md).

- `se_antigen_table`:

  Standard-error lookup table across all antigens.

- `current_se`:

  SE values for the currently selected antigen.

- `.selection_args`:

  Internal cache of arguments from the last `$set_curve_settings()`.

## Methods

### Public methods

- [`StandardCurve$new()`](#method-StandardCurve-new)

- [`StandardCurve$set_data()`](#method-StandardCurve-set_data)

- [`StandardCurve$set_curve_settings()`](#method-StandardCurve-set_curve_settings)

- [`StandardCurve$fit()`](#method-StandardCurve-fit)

- [`StandardCurve$propagate_error()`](#method-StandardCurve-propagate_error)

- [`StandardCurve$summarize()`](#method-StandardCurve-summarize)

- [`StandardCurve$plot()`](#method-StandardCurve-plot)

- [`StandardCurve$compare_models()`](#method-StandardCurve-compare_models)

- [`StandardCurve$get_results()`](#method-StandardCurve-get_results)

- [`StandardCurve$print()`](#method-StandardCurve-print)

- [`StandardCurve$clone()`](#method-StandardCurve-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new `StandardCurve` object.

#### Usage

    StandardCurve$new(
      loaded_data,
      study_params = list(applyProzone = TRUE, blank_option = "ignored", is_log_response =
        TRUE, is_log_independent = TRUE),
      antigen_constraints = NULL,
      model_names = NULL,
      curve_col = "curve_id",
      is_display_log_response = TRUE,
      is_display_log_independent = TRUE,
      verbose = TRUE
    )

#### Arguments

- `loaded_data`:

  List of pre-loaded assay data. Must contain at minimum: `$standards`,
  `$blanks`, `$samples`, `$response_var`, `$indep_var`. Optionally
  `$antigen_constraints` (used if `antigen_constraints` arg is `NULL`).

- `study_params`:

  Named list with keys: `applyProzone`, `blank_option`,
  `is_log_response`, `is_log_independent`.

- `antigen_constraints`:

  Data frame of per-antigen constraints. If `NULL`, falls back to
  `loaded_data$antigen_constraints`.

- `model_names`:

  Character vector of candidate model names. Defaults to all five
  supported models.

- `curve_col`:

  Name of the curve-ID column. Default `"curve_id"`.

- `is_display_log_response`:

  Logical. Default `TRUE`.

- `is_display_log_independent`:

  Logical. Default `TRUE`.

- `verbose`:

  Logical. Default `TRUE`.

- `curve_id_element_order`:

  Character vector defining field order inside `curve_id`. Defaults to
  the standard nine-element schema.

#### Returns

A `StandardCurve` object (invisibly).

#### Examples

    \dontrun{
    loaded_data <- pull_data(
      study_accession      = "MADI_01",
      experiment_accession = "IgG1",
      project_id           = 17,
      conn                 = conn
    )

    sc <- StandardCurve$new(
      loaded_data         = loaded_data,
      study_params        = list(
        applyProzone       = TRUE,
        blank_option       = "ignored",
        is_log_response    = TRUE,
        is_log_independent = TRUE
      ),
      antigen_constraints = my_constraints_df
    )
    }

------------------------------------------------------------------------

### Method `set_data()`

Replace the loaded data (and optionally constraints) without creating a
new object. Resets all downstream results automatically.

Useful when iterating over multiple studies with a single
`StandardCurve` instance.

#### Usage

    StandardCurve$set_data(loaded_data, antigen_constraints = NULL)

#### Arguments

- `loaded_data`:

  New pre-loaded data list.

- `antigen_constraints`:

  Optional replacement constraints data frame. Falls back to
  `loaded_data$antigen_constraints` or the existing value.

#### Returns

The `StandardCurve` object (invisibly).

#### Examples

    \dontrun{
    loaded_data_2 <- pull_data("MADI_02", "IgG1", 17, conn)
    sc$set_data(loaded_data_2)
    sc$select(...)$fit()
    }

------------------------------------------------------------------------

### Method `set_curve_settings()`

Resolve antigen-specific constraints and curve settings.

Wraps
[`resolve_curve_settings()`](https://immunoplex.github.io/curveRfreq/reference/resolve_curve_settings.md)
using `self$loaded_data`, which must already be filtered to the relevant
curve prior to calling this method. Filtering is performed externally
via `fetch_curve_id()` and `filter_by_curve_id()` before constructing
the object with `$new()` or updating it with `$set_data()`.

Results are stored in `$antigen_plate`. Any stale fit results from a
previous `$set_curve_settings()` call are cleared automatically.

#### Usage

    StandardCurve$set_curve_settings(wavelength = WL_NONE)

#### Arguments

- `wavelength`:

  Character. Defaults to `WL_NONE` (`"__none__"`) for bead arrays. Set
  explicitly for ELISA absorbance wavelengths, e.g. `"450"`.

#### Returns

The `StandardCurve` object (invisibly).

#### Examples

    \dontrun{
    curve_id <- fetch_curve_id(
      lookup                  = loaded_data$curve_id_lookup,
      project_id              = 17,
      study_accession         = "MADI_01",
      experiment_accession    = "IgG1",
      feature                 = "IgG1",
      source                  = "Standard",
      antigen                 = "victoria",
      plate                   = "plate_3",
      nominal_sample_dilution = "1000",
      wavelength              = "__none__",
      element_order           = curve_id_element_order
    )

    filtered <- filter_by_curve_id(loaded_data, curve_id = curve_id)

    sc$set_curve_settings()
    # or for ELISA:
    sc$set_curve_settings(wavelength = "450")
    }

------------------------------------------------------------------------

### Method `fit()`

Run the complete model-fitting pipeline on the selected antigen plate.

Sequentially calls:

1.  [`preprocess_robust_curves()`](https://immunoplex.github.io/curveRfreq/reference/preprocess_robust_curves.md)

2.  [`select_model_formulas()`](https://immunoplex.github.io/curveRfreq/reference/select_model_formulas.md)

3.  [`obtain_model_constraints()`](https://immunoplex.github.io/curveRfreq/reference/obtain_model_constraints.md)

4.  [`make_start_lists()`](https://immunoplex.github.io/curveRfreq/reference/make_start_lists.md)

5.  [`compute_robust_curves()`](https://immunoplex.github.io/curveRfreq/reference/compute_robust_curves.md)

6.  [`summarize_model_fits()`](https://immunoplex.github.io/curveRfreq/reference/summarize_model_fits.md)

7.  [`summarize_model_parameters()`](https://immunoplex.github.io/curveRfreq/reference/summarize_model_parameters.md)

8.  [`get_plot_data()`](https://immunoplex.github.io/curveRfreq/reference/get_plot_data.md)

9.  [`select_model_fit_AIC()`](https://immunoplex.github.io/curveRfreq/reference/select_model_fit_AIC.md)

10. [`tidy.nlsLM()`](https://immunoplex.github.io/curveRfreq/reference/tidy.nlsLM.md)

11. [`summarize_fit()`](https://immunoplex.github.io/curveRfreq/reference/summarize_fit.md)

#### Usage

    StandardCurve$fit(start_quantiles = c(low = 0.2, mid = 0.5, high = 0.8))

#### Arguments

- `start_quantiles`:

  Named numeric vector `c(low, mid, high)`. Default
  `c(low = 0.2, mid = 0.5, high = 0.8)`.

#### Returns

The `StandardCurve` object (invisibly).

#### Examples

    \dontrun{
    sc$fit()
    }

------------------------------------------------------------------------

### Method `propagate_error()`

Compute the SE lookup table and propagate measurement error from the
standard curve to the precision profile.

Calls
[`compute_antigen_se_table()`](https://immunoplex.github.io/curveRfreq/reference/compute_antigen_se_table.md),
[`lookup_antigen_se()`](https://immunoplex.github.io/curveRfreq/reference/lookup_antigen_se.md),
and
[`predict_and_propagate_error()`](https://immunoplex.github.io/curveRfreq/reference/predict_and_propagate_error.md).
Results stored in `$se_antigen_table`, `$current_se`, and
`$best_fit$sample_se`.

#### Usage

    StandardCurve$propagate_error()

#### Arguments

- `grouping_cols`:

  Character vector of columns for SE grouping.

#### Returns

The `StandardCurve` object (invisibly).

#### Examples

    \dontrun{
    sc$propagate_error()
    head(sc$best_fit$sample_se)
    }

------------------------------------------------------------------------

### Method `summarize()`

Print a human-readable summary of the fitted standard curve including QC
metrics, goodness-of-fit statistics, and parameter estimates.

#### Usage

    StandardCurve$summarize(digits = 4)

#### Arguments

- `digits`:

  Integer. Significant digits to display. Default `4`.

#### Returns

The `StandardCurve` object (invisibly).

#### Examples

    \dontrun{
    sc$summarize()
    }

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Plot the fitted standard curve with QC metric annotations.

Wraps
[`plot_standard_curve()`](https://immunoplex.github.io/curveRfreq/reference/plot_standard_curve.md).

#### Usage

    StandardCurve$plot(
      is_display_log_independent = self$is_display_log_independent,
      is_display_log_response = self$is_display_log_response
    )

#### Arguments

- `is_display_log_independent`:

  Logical. Log-scale x-axis. Defaults to the self of the display
  independent variable axis parameter.

- `is_display_log_response`:

  Logical. Log-scale y-axis. Defaults to the self of the display
  response parameter from constructor.

#### Returns

The ggplot object (invisibly).

#### Examples

    \dontrun{
    sc$plot()
    }

------------------------------------------------------------------------

### Method `compare_models()`

Multi-panel plot comparing all candidate model fits, residuals,
parameter estimates, and AIC scores.

Wraps
[`plot_model_comparisons()`](https://immunoplex.github.io/curveRfreq/reference/plot_model_comparisons.md).

#### Usage

    StandardCurve$compare_models(use_patchwork = TRUE)

#### Arguments

- `use_patchwork`:

  Logical. Use patchwork layout. Default `TRUE`.

#### Returns

The plot object (invisibly).

#### Examples

    \dontrun{
    sc$compare_models()
    }

------------------------------------------------------------------------

### Method `get_results()`

Return key result tables as a named list for downstream use or export.

#### Usage

    StandardCurve$get_results()

#### Returns

Named list: `fit_summary`, `best_parameters`, `best_fit_summary`,
`best_pred`,`best_standard`,`sample_se` (populated only after
`$propagate_error()`).

#### Examples

    \dontrun{
    res <- sc$get_results()
    write.csv(res$best_fit_summary, "summary.csv", row.names = FALSE)
    }

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print a concise pipeline-status overview.

#### Usage

    StandardCurve$print(...)

#### Arguments

- `...`:

  Ignored. Required for S3 print compatibility.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    StandardCurve$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
## ------------------------------------------------
## Method `StandardCurve$new`
## ------------------------------------------------

if (FALSE) { # \dontrun{
loaded_data <- pull_data(
  study_accession      = "MADI_01",
  experiment_accession = "IgG1",
  project_id           = 17,
  conn                 = conn
)

sc <- StandardCurve$new(
  loaded_data         = loaded_data,
  study_params        = list(
    applyProzone       = TRUE,
    blank_option       = "ignored",
    is_log_response    = TRUE,
    is_log_independent = TRUE
  ),
  antigen_constraints = my_constraints_df
)
} # }

## ------------------------------------------------
## Method `StandardCurve$set_data`
## ------------------------------------------------

if (FALSE) { # \dontrun{
loaded_data_2 <- pull_data("MADI_02", "IgG1", 17, conn)
sc$set_data(loaded_data_2)
sc$select(...)$fit()
} # }

## ------------------------------------------------
## Method `StandardCurve$set_curve_settings`
## ------------------------------------------------

if (FALSE) { # \dontrun{
curve_id <- fetch_curve_id(
  lookup                  = loaded_data$curve_id_lookup,
  project_id              = 17,
  study_accession         = "MADI_01",
  experiment_accession    = "IgG1",
  feature                 = "IgG1",
  source                  = "Standard",
  antigen                 = "victoria",
  plate                   = "plate_3",
  nominal_sample_dilution = "1000",
  wavelength              = "__none__",
  element_order           = curve_id_element_order
)

filtered <- filter_by_curve_id(loaded_data, curve_id = curve_id)

sc$set_curve_settings()
# or for ELISA:
sc$set_curve_settings(wavelength = "450")
} # }

## ------------------------------------------------
## Method `StandardCurve$fit`
## ------------------------------------------------

if (FALSE) { # \dontrun{
sc$fit()
} # }

## ------------------------------------------------
## Method `StandardCurve$propagate_error`
## ------------------------------------------------

if (FALSE) { # \dontrun{
sc$propagate_error()
head(sc$best_fit$sample_se)
} # }

## ------------------------------------------------
## Method `StandardCurve$summarize`
## ------------------------------------------------

if (FALSE) { # \dontrun{
sc$summarize()
} # }

## ------------------------------------------------
## Method `StandardCurve$plot`
## ------------------------------------------------

if (FALSE) { # \dontrun{
sc$plot()
} # }

## ------------------------------------------------
## Method `StandardCurve$compare_models`
## ------------------------------------------------

if (FALSE) { # \dontrun{
sc$compare_models()
} # }

## ------------------------------------------------
## Method `StandardCurve$get_results`
## ------------------------------------------------

if (FALSE) { # \dontrun{
res <- sc$get_results()
write.csv(res$best_fit_summary, "summary.csv", row.names = FALSE)
} # }
```
