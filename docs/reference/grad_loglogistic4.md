# Analytical Gradient of the Inverse loglogistic4

Returns partial derivatives for
[`inv_loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic4.md)
with all parameters free.

## Usage

``` r
grad_loglogistic4(y, a, b, c, d)
```

## Arguments

- y:

  Numeric scalar. Observed response.

- a, b, c, d:

  Numeric scalars. Free loglogistic4 model parameters.

## Value

A list with `grad_theta` (named vector `c(a=, b=, c=, d=)`) and scalar
`grad_y`.

## See also

[`inv_loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic4.md),
[`grad_inv_loglogistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_loglogistic4_fixed.md)

Other gradient-functions:
[`grad_gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/grad_gompertz4.md),
[`grad_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/grad_logistic4.md),
[`grad_logistic5()`](https://immunoplex.github.io/curveRfreq/reference/grad_logistic5.md),
[`grad_loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/grad_loglogistic5.md)
