# Compute Shape/curvature-based Lower and Upper Limits of Quantification (LOQ) based on the second derivative

Identifies the lower limit of quantification (LLOQ) and upper limit of
quantification (ULOQ) based on local extrema of the second derivative of
a fitted curve. Candidate extrema are refined using quadratic
interpolation, and corresponding response values are obtained from the
fitted model.

## Usage

``` r
compute_loqs(best_d2xy, fit, independent_variable, verbose = TRUE)
```

## Arguments

- best_d2xy:

  A data.frame containing second derivative (found in best fit object)
  information with:

  x

  :   Numeric vector of independent variable values.

  d2x_y

  :   Numeric vector of second derivative values evaluated at `x`.

- fit:

  A fitted model object supporting
  [`predict()`](https://rdrr.io/r/stats/predict.html).

- independent_variable:

  Character string naming the predictor variable used in the model.

- verbose:

  Logical; if `TRUE`, prints diagnostic messages.

## Value

A named list with:

- lloq:

  Estimated lower limit of quantification (x-value).

- uloq:

  Estimated upper limit of quantification (x-value).

- lloq_y:

  Predicted response at LLOQ.

- uloq_y:

  Predicted response at ULOQ.

## Details

The method proceeds as follows:

- Computes first differences of the second derivative to detect sign
  changes.

- Identifies candidate local maxima and minima of the second derivative.

- Applies quadratic interpolation using three neighboring points to
  refine the location of each extremum.

- Selects:

  - LLOQ as the x-value corresponding to the maximum interpolated second
    derivative

  - ULOQ as the x-value corresponding to the minimum interpolated second
    derivative

- Evaluates the fitted model at these x-values to obtain corresponding
  y-values.

If no valid extrema are found, the corresponding outputs are returned as
`NA`.

## See also

[`generate_lods`](https://immunoplex.github.io/curveRfreq/reference/generate_lods.md),
[`generate_mdc_rdl`](https://immunoplex.github.io/curveRfreq/reference/generate_mdc_rdl.md)

## Examples

``` r
if (FALSE) { # \dontrun{
loqs <- compute_loqs(
  best_d2xy = d2_data,
  fit = model_fit,
  independent_variable = "log_dilution"
)

loqs$lloq
loqs$uloq
} # }
```
