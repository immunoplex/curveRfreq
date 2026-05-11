# Attach grouping keys to output data

Ensures that `wavelength` and `feature` columns exist in an output data
frame. Values are sourced from the corresponding `best_data` used to
generate the output.

## Usage

``` r
attach_grouping_keys(df, best_data, context = "")
```

## Arguments

- df:

  A data.frame containing model output or summary results.

- best_data:

  A data.frame used to derive grouping keys (typically input to model
  fitting).

- context:

  Optional character string used for logging/debugging context.

## Value

A data.frame with `wavelength` and `feature` columns guaranteed to
exist.

## Details

This function enforces consistent natural keys across all outputs
generated from model fitting workflows. If keys are missing:

- `wavelength` is filled using `best_data` or `WL_NONE`

- `feature` is filled using `best_data` or `FEAT_NONE`

The `wavelength` column is normalized using
[`normalize_wavelength()`](https://immunoplex.github.io/curveRfreq/reference/normalize_wavelength.md).

## See also

[`normalize_wavelength`](https://immunoplex.github.io/curveRfreq/reference/normalize_wavelength.md)
