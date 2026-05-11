# Gradient of the Inverse loglogistic4 with Fixed Lower Asymptote

Gradient of the Inverse loglogistic4 with Fixed Lower Asymptote

## Usage

``` r
grad_inv_loglogistic4_fixed(y, fixed_a, b, c, d)
```

## Arguments

- y:

  Numeric scalar. Observed response.

- fixed_a:

  Numeric scalar. Externally fixed lower asymptote.

- b, c, d:

  Numeric scalars. Free model parameters.

## Value

Named numeric vector `c(b=, c=, d=)`.

## See also

Other fixed-a-gradients:
[`grad_inv_gompertz4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_gompertz4_fixed.md),
[`grad_inv_logistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_logistic4_fixed.md),
[`grad_inv_logistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_logistic5_fixed.md),
[`grad_inv_loglogistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_loglogistic5_fixed.md)
