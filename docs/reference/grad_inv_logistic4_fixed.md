# Gradient of the Inverse 4PL with Fixed Lower Asymptote

Returns partial derivatives \\\partial x / \partial \theta\\ for
[`inv_logistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4_fixed.md)
where `a` is an externally fixed constant and only `b`, `c`, `d` are
free.

## Usage

``` r
grad_inv_logistic4_fixed(y, fixed_a, b, c, d)
```

## Arguments

- y:

  Numeric scalar. Observed response.

- fixed_a:

  Numeric scalar. Externally fixed lower asymptote.

- b, c, d:

  Numeric scalars. Free model parameters.

## Value

Named numeric vector `c(b=, c=, d=)` of partial derivatives.

## See also

[`grad_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/grad_logistic4.md),
[`inv_logistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4_fixed.md),
[`grad_y_logistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_y_logistic4_fixed.md)

Other fixed-a-gradients:
[`grad_inv_gompertz4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_gompertz4_fixed.md),
[`grad_inv_logistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_logistic5_fixed.md),
[`grad_inv_loglogistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_loglogistic4_fixed.md),
[`grad_inv_loglogistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_loglogistic5_fixed.md)
