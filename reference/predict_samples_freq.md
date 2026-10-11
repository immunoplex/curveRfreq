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
  pcov_threshold = 20,
  is_log_independent = TRUE,
  lloq_log10 = NA_real_,
  uloq_log10 = NA_real_,
  decision_cutoff = NULL,
  n_interval_draws = 1000L,
  interval_seed = NULL,
  param_draws = NULL,
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

- pcov_threshold:

  Numeric. Precision threshold (%) used to set `pcov_pass`. A grid point
  passes when `pcov < pcov_threshold`. Default 20. Must be
  `<= cv_x_max`.

- is_log_independent:

  Logical. Is x on log10 scale? Default TRUE.

- lloq_log10, uloq_log10:

  Numeric. This curve's LLOQ/ULOQ (log10, pre-dilution) for
  `p_below_lloq` / `p_above_uloq`. Default `NA`.

- decision_cutoff:

  Numeric or NULL. Dilution-corrected concentration for
  `p_below_cutoff`. Default NULL (column is `NA`).

- n_interval_draws:

  Integer. Monte Carlo draws per sample for the predictive interval and
  probabilities. Default 1000; 0 disables them (columns `NA`).

- interval_seed:

  Integer or NULL. Seed for those draws. The global RNG state is never
  changed either way.

- param_draws:

  Matrix or NULL. Pre-drawn parameter draws (rows = draws, columns named
  like `coef(fit)`), as made by
  [`fit_calibration_freq()`](https://immunoplex.github.io/curveRfreq/reference/fit_calibration_freq.md)
  so the intervals and the persisted `calib_draws` share one set. NULL
  draws them here from `coef(fit)` / `vcov(fit)`.

- verbose:

  Logical.

## Value

Data frame with original sample columns plus: `observed_response_fit`,
`predicted_log10_concentration`, `predicted_concentration`,
`final_concentration`, `se_concentration`, `pcov`, `pcov_rmse`, and the
[`curveRcore::summarize_conc_draws()`](https://immunoplex.github.io/curveRcore/reference/summarize_conc_draws.html)
columns.

## Details

The interval / probability columns come from a Monte Carlo predictive
distribution, not the delta method: parameters are drawn from their
asymptotic multivariate normal (`coef(fit)`, `vcov(fit)`), residual
noise `N(0, se_response)` is added to the observed response (the same
variance term the delta method uses), and each draw is inverted with
[`curveRcore::invert_with_bounds()`](https://immunoplex.github.io/curveRcore/reference/invert_with_bounds.html),
so draws beyond an asymptote count as censored.
[`curveRcore::summarize_conc_draws()`](https://immunoplex.github.io/curveRcore/reference/summarize_conc_draws.html)
then summarises them, exactly as for Bayesian fits. `se_concentration` /
`pcov` are unchanged (delta method).
