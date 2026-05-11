# Analytical Gradient of the Inverse 4PL

Returns partial derivatives \\\partial x / \partial \theta\\ and
\\\partial x / \partial y\\ for the
[`inv_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4.md)
function when all four parameters are free. Consumed by
[`propagate_error_analytic()`](https://immunoplex.github.io/curveRfreq/reference/propagate_error_analytic.md)
for delta-method error propagation.

## Usage

``` r
grad_logistic4(y, a, b, c, d)
```

## Arguments

- y:

  Numeric scalar. Observed response value.

- a, b, c, d:

  Numeric scalars. Free 4PL model parameters.

## Value

A list with two elements:

- grad_theta:

  Named numeric vector `c(a=, b=, c=, d=)` — partial derivatives of `x`
  with respect to each model parameter.

- grad_y:

  Numeric scalar — \\\partial x / \partial y\\.

## See also

[`inv_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4.md),
[`grad_inv_logistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_logistic4_fixed.md)

Other gradient-functions:
[`grad_gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/grad_gompertz4.md),
[`grad_logistic5()`](https://immunoplex.github.io/curveRfreq/reference/grad_logistic5.md),
[`grad_loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/grad_loglogistic4.md),
[`grad_loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/grad_loglogistic5.md)
