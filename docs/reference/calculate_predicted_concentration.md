# Back-Calculate Concentration from Fitted Model for Each Sample

Applies the analytical inverse of the selected sigmoid model to each
observed response value in `plate_samples` to estimate the corresponding
concentration.

## Usage

``` r
calculate_predicted_concentration(
  model_name,
  fit,
  plate_samples,
  fixed_constraint,
  response_variable,
  is_log_response,
  verbose = TRUE
)
```

## Arguments

- model_name:

  Character. One of `"logistic5"`, `"loglogistic5"`, `"logistic4"`,
  `"loglogistic4"`, `"gompertz4"`.

- fit:

  Fitted `nlsLM` object.

- plate_samples:

  Data frame of sample wells with the response column.

- fixed_constraint:

  Numeric or `NULL`. Fixed lower asymptote on the model's fitting scale.

- response_variable:

  Character. Name of the response column.

- is_log_response:

  Logical. Was the response log10-transformed before fitting?

- verbose:

  Logical (default `TRUE`).

## Value

`plate_samples` with a new column `predicted_concentration`.
