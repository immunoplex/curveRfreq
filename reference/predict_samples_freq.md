# Predict Concentration for Test Samples

For each test sample, applies the inverse of the best-fit model to the
observed response to get predicted concentration, then propagates
uncertainty via the delta method to get `se_concentration` and `pcov`.

## Usage

``` r
predict_samples_freq(
  samples,
  fit,
  model_name,
  response_variable,
  fixed_a = NULL,
  is_log_response = TRUE,
  se_response = 0,
  cv_x_max = 150,
  is_log_independent = TRUE,
  verbose = FALSE
)
```

## Arguments

- samples:

  Data frame of test samples. Must contain the response column and
  `dilution`.

- fit:

  The best-fit `nlsLM` object.

- model_name:

  Character. Best model name.

- response_variable:

  Character. Response column name.

- fixed_a:

  Numeric or NULL. Fixed lower asymptote on fitting scale.

- is_log_response:

  Logical. Was the response log10-transformed?

- se_response:

  Numeric. Residual SE on the fitting scale. Typically
  `summary(fit)$sigma`. Default 0.

- cv_x_max:

  Numeric. Cap for pcov. Default 150.

- is_log_independent:

  Logical. Is x on log10 scale? Default TRUE.

- verbose:

  Logical.

## Value

Data frame with original sample columns plus: `observed_response_fit`,
`predicted_log10_concentration`, `predicted_concentration`,
`final_concentration`, `se_concentration`, `pcov`, `pcov_rmse`.
