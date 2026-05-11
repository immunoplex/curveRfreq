# Compute delta-method 95% confidence interval bands for a fitted curve

Calculates confidence interval bands for a nonlinear least squares
fitted curve using the delta method. Performs exactly `p_free`
vectorised formula evaluations over the full `x_new` grid (one per free
parameter) rather than the naive `n × p_free` scalar
[`predict()`](https://rdrr.io/r/stats/predict.html) calls. For a
200-point grid and a 5-parameter model this is approximately 40x fewer
evaluations.

## Usage

``` r
compute_curve_ci(fit, x_new, x_var, fixed_a = NULL, level = 0.95)
```

## Arguments

- fit:

  A converged `nls` or `nlsLM` object.

- x_new:

  Numeric vector of x values at which to evaluate the curve.

- x_var:

  Character string naming the independent variable in the model formula.

- fixed_a:

  Numeric scalar or `NULL`. Fixed lower asymptote not included in
  [`vcov()`](https://rdrr.io/r/stats/vcov.html). Default is `NULL`.

- level:

  Numeric confidence level. Default is `0.95`.

## Value

A `data.frame` with columns:

- x:

  The x values supplied in `x_new`.

- ci_lo:

  Lower confidence band.

- ci_hi:

  Upper confidence band.

If `vcov(fit)` fails, `ci_lo` and `ci_hi` are `NA`.

## See also

[`nls`](https://rdrr.io/r/stats/nls.html),
[`vcov`](https://rdrr.io/r/stats/vcov.html)

## Examples

``` r
if (FALSE) { # \dontrun{
ci <- compute_curve_ci(fit = my_nls_fit, x_new = seq(0, 10, length.out = 200),
                       x_var = "concentration", fixed_a = NULL)
} # }
```
