# Preprocess Standard Curve Data for Robust Curve Fitting

Orchestrates the full preprocessing chain:

1.  Compute concentrations from dilutions.

2.  Optionally apply prozone correction.

3.  Apply blank operation (ignore / include / subtract (1, 3, or 10
    times the geometric mean of blanks)).

4.  Floor non-positive values and optionally log10-transform the
    response.

## Usage

``` r
preprocess_robust_curves(
  data,
  antigen_settings,
  response_variable,
  independent_variable,
  is_log_response,
  blank_data = NULL,
  blank_option = "ignored",
  is_log_independent = TRUE,
  apply_prozone = TRUE,
  verbose = TRUE
)
```

## Arguments

- data:

  Data frame with columns `dilution` and the response variable.

- antigen_settings:

  Named list from
  [`obtain_lower_constraint`](https://immunoplex.github.io/curveRfreq/reference/obtain_lower_constraint.md).

- response_variable:

  Character. Response column name.

- independent_variable:

  Character. Concentration column name.

- is_log_response:

  Logical. Log10-transform the response?

- blank_data:

  Data frame of blank measurements, or `NULL`.

- blank_option:

  Character. Blank handling method.One of `"ignored"`, `"included"`,
  `"subtracted"`, `"subtracted_3x"`, `"subtracted_10x"`. (default
  `"ignored"`).

- is_log_independent:

  Logical. Log10-transform the concentration? (default `TRUE`).

- apply_prozone:

  Logical. Apply prozone correction? (default `TRUE`).

- verbose:

  Logical (default `TRUE`).

## Value

A list with:

- data:

  Preprocessed data frame.

- antigen_fit_options:

  Named list of the options used, for passing to downstream functions.
