# Select the Best Model by AIC

Identifies the converged model with the lowest AIC and packages it into
the `selection` component of a `calibration_result`.

## Usage

``` r
select_best_aic(ensemble)
```

## Arguments

- ensemble:

  Output of
  [`fit_ensemble_nls()`](https://immunoplex.github.io/curveRfreq/reference/fit_ensemble_nls.md).

## Value

A named list with:

- best_model_name:

  Character, or NA if none converged.

- criterion:

  "AIC"

- weights:

  Data frame with model_name, aic, delta_aic, converged.
