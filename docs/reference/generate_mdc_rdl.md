# Compute MDC and RDL values from fitted model and LODs

Computes minimum and maximum detectable concentrations (MDC) and
reliable detection limits (RDL) by solving for concentrations where the
fitted curve or its confidence bounds intersect the LODs.

## Usage

``` r
generate_mdc_rdl(best_fit, lods, independent_variable, verbose = TRUE)
```

## Arguments

- best_fit:

  A list containing:

  best_fit

  :   A fitted model object supporting
      [`predict()`](https://rdrr.io/r/stats/predict.html),
      [`coef()`](https://rdrr.io/r/stats/coef.html), and
      [`vcov()`](https://rdrr.io/r/stats/vcov.html).

  best_data

  :   A data.frame used for fitting the model.

- lods:

  A list containing `llod` and `ulod`, typically from
  [`generate_lods`](https://immunoplex.github.io/curveRfreq/reference/generate_lods.md).

- independent_variable:

  Character string naming the independent variable column in
  `best_data`.

- verbose:

  Logical; if `TRUE`, prints computed MDC and RDL values.

## Value

A named list with:

- mindc:

  Minimum detectable concentration (fitted curve = LLOD).

- maxdc:

  Maximum detectable concentration (fitted curve = ULOD).

- minrdl:

  Minimum reliable detection limit (lower CI = LLOD).

- maxrdl:

  Maximum reliable detection limit (upper CI = ULOD).

## Details

Root-finding is performed using
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html) over the observed
range of the independent variable.

The prediction standard error is computed using the delta method:
\$\$SE(\hat{y}) = \sqrt{\nabla f(\theta)^T V \nabla f(\theta)}\$\$ where
\\V\\ is the variance-covariance matrix of the fitted parameters.

Reliable detection limits (RDL) are defined using the confidence
interval:

- Lower bound: \\\hat{y}(x) - t\_{\alpha/2} \cdot SE(x)\\

- Upper bound: \\\hat{y}(x) + t\_{\alpha/2} \cdot SE(x)\\

If no valid root is found within the observed data range, the
corresponding value is returned as `NA`.

## See also

[`generate_lods`](https://immunoplex.github.io/curveRfreq/reference/generate_lods.md)

## Examples

``` r
if (FALSE) { # \dontrun{
lods <- generate_lods(best_fit)
mdc_rdl <- generate_mdc_rdl(best_fit, lods, independent_variable = "log_dilution")
} # }
```
