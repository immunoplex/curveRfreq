# Analytical Gradient of the Inverse Gompertz

Returns partial derivatives for
[`inv_gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/inv_gompertz4.md)
with all parameters free.

## Usage

``` r
grad_gompertz4(y, a, b, c, d)
```

## Arguments

- y:

  Numeric scalar. Observed response.

- a, b, c, d:

  Numeric scalars. Free Gompertz model parameters.

## Value

A list with `grad_theta` (named vector `c(a=, b=, c=, d=)`) and scalar
`grad_y`.

## See also

[`inv_gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/inv_gompertz4.md),
[`grad_inv_gompertz4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_gompertz4_fixed.md)

Other gradient-functions:
[`grad_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/grad_logistic4.md),
[`grad_logistic5()`](https://immunoplex.github.io/curveRfreq/reference/grad_logistic5.md),
[`grad_loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/grad_loglogistic4.md),
[`grad_loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/grad_loglogistic5.md)
