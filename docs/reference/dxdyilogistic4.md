# \\dx/dy\\ for the 4PL Inverse

Scalar derivative of the
[`inv_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4.md)
function with respect to `y`, used by
[`propagate_error_analytic()`](https://immunoplex.github.io/curveRfreq/reference/propagate_error_analytic.md)
for error propagation through the measurement uncertainty component.

## Usage

``` r
dxdyilogistic4(y, a, b, c, d)
```

## Arguments

- y:

  Numeric scalar. Observed response.

- a, b, c, d:

  Numeric scalars. 4PL model parameters.

## Value

Numeric scalar \\dx/dy\\.

## See also

Other dxdy-helpers:
[`dxdyigompertz4()`](https://immunoplex.github.io/curveRfreq/reference/dxdyigompertz4.md),
[`dxdyilogistic5()`](https://immunoplex.github.io/curveRfreq/reference/dxdyilogistic5.md),
[`dxdyiloglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/dxdyiloglogistic4.md),
[`dxdyiloglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/dxdyiloglogistic5.md)
