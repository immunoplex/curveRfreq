# Select Best Multi-Start NLS Fit by AIC

Iterates over a named list of start-value sets, fits the model for each,
and returns the fit with the lowest AIC.

## Usage

``` r
select_nlsLM_aic(
  prepped_data,
  response_variable,
  independent_variable,
  formula,
  lower_model_constraints,
  upper_model_constraints,
  start_lists,
  verbose = TRUE
)
```

## Arguments

- prepped_data:

  Data frame.

- response_variable:

  Character. Response column name.

- independent_variable:

  Character. Concentration column name.

- formula:

  A `formula` object.

- lower_model_constraints:

  Named numeric vector of lower bounds.

- upper_model_constraints:

  Named numeric vector of upper bounds.

- start_lists:

  Named list of start-value lists.

- verbose:

  Logical (default `TRUE`).

## Value

Best `nlsLM` fit object, or `NULL` if all starts fail.
