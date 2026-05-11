# Tidy parameter estimates from a fitted nlsLM model

Extracts and tidies coefficient estimates from a fitted `nlsLM` model,
attaches parameter constraints, handles fixed lower asymptote values,
and appends study metadata. When the lower asymptote is fixed via the
`"range_of_blanks"` method and the response is log-transformed, the
fixed value and its constraint bounds are log10-transformed before
tidying.

## Usage

``` r
tidy.nlsLM(
  best_fit,
  fixed_a_result,
  model_constraints,
  antigen_settings,
  antigen_fit_options,
  verbose = TRUE
)
```

## Arguments

- best_fit:

  A named list as returned by the model selection workflow, containing
  at minimum:

  best_fit

  :   A converged `nlsLM` fit object.

  best_model_name

  :   Character string naming the selected model (e.g. `"logistic5"`).

  best_data

  :   A `data.frame` of the data used for fitting, containing columns
      `study_accession`, `experiment_accession`,
      `nominal_sample_dilution`, `antigen`, `plateid`, `plate`, and
      `source`.

- fixed_a_result:

  Numeric scalar or `NULL`. Fixed lower asymptote value. If `NULL` the
  lower asymptote is treated as a free parameter. When non-`NULL` and
  the response is log-transformed, the value is log10-transformed
  internally before being appended as a fixed row.

- model_constraints:

  Named list of parameter constraint matrices, one entry per model name.
  Each entry must be coercible to a `data.frame` with rownames
  corresponding to parameter terms and columns `lower` and `upper`.

- antigen_settings:

  A named list of antigen-level settings, requiring:

  l_asy_constraint_method

  :   Character. Method used to determine the lower asymptote constraint
      (e.g. `"range_of_blanks"`).

  l_asy_min_constraint

  :   Numeric. Lower bound for the asymptote constraint.

  l_asy_max_constraint

  :   Numeric. Upper bound for the asymptote constraint.

- antigen_fit_options:

  A named list of fitting options, requiring:

  is_log_response

  :   Logical. Whether the response variable is log10-transformed.

- verbose:

  Logical. If `TRUE` (default), prints a completion message via
  [`message()`](https://rdrr.io/r/base/message.html).

## Value

The input `best_fit` list with an additional element `best_parameters`:
a `tibble` with one row per parameter containing columns:

- curve_id:

  The curve identifier for the best fitted curve. It is associated with
  an identifier value in the lookup table.

- term:

  Parameter name.

- lower:

  Lower constraint bound from `model_constraints`.

- upper:

  Upper constraint bound from `model_constraints`.

- estimate:

  Parameter estimate.

- std_error:

  Standard error of the estimate (0 for fixed parameters).

- statistic:

  t-statistic (`NA` for fixed parameters).

- p_value:

  Two-sided p-value (`NA` for fixed parameters).

## Details

When `fixed_a_result` is non-`NULL`, a synthetic row for the fixed `"a"`
parameter is prepended to the tidy output with `std_error = 0` and
`statistic`/`p_value` set to `NA`. Constraint bounds for the fixed
asymptote are sourced from `antigen_settings$l_asy_min_constraint` and
`antigen_settings$l_asy_max_constraint`.

If `antigen_settings$l_asy_constraint_method == "range_of_blanks"` and
`antigen_fit_options$is_log_response` is `TRUE`, the fixed asymptote and
its bounds are log10-transformed before use. Constraint bounds that are
non-positive or non-finite are coerced to `NA_real_`.

## See also

[`attach_grouping_keys`](https://immunoplex.github.io/curveRfreq/reference/attach_grouping_keys.md),
[`validate_fixed_lower_asymptote`](https://immunoplex.github.io/curveRfreq/reference/validate_fixed_lower_asymptote.md)

## Examples

``` r
if (FALSE) { # \dontrun{
best_fit <- tidy.nlsLM(
  best_fit          = model_result,
  fixed_a_result    = 50,
  model_constraints = my_constraints,
  antigen_settings  = my_antigen_settings,
  antigen_fit_options = list(is_log_response = TRUE)
)
best_fit$best_parameters
} # }
```
