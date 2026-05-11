# Analytical Gradient of the Inverse loglogistic5

Returns partial derivatives for
[`inv_loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic5.md)
with all five parameters free.

## Usage

``` r
grad_loglogistic5(y, a, b, c, d, g)
```

## Arguments

- y:

  Numeric scalar. Observed response.

- a, b, c, d:

  Numeric scalars. Free loglogistic5 model parameters.

- g:

  Numeric scalar. Asymmetry parameter.

## Value

A list with `grad_theta` (named vector `c(a=, b=, c=, d=, g=)`) and
scalar `grad_y`.

## See also

[`inv_loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic5.md),
[`grad_inv_loglogistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_loglogistic5_fixed.md)

Other gradient-functions:
[`grad_gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/grad_gompertz4.md),
[`grad_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/grad_logistic4.md),
[`grad_logistic5()`](https://immunoplex.github.io/curveRfreq/reference/grad_logistic5.md),
[`grad_loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/grad_loglogistic4.md)
