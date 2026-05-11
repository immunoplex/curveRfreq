# Gradient of the Inverse 5PL with Fixed Lower Asymptote

Gradient of the Inverse 5PL with Fixed Lower Asymptote

## Usage

``` r
grad_inv_logistic5_fixed(y, fixed_a, b, c, d, g)
```

## Arguments

- y:

  Numeric scalar. Observed response.

- fixed_a:

  Numeric scalar. Externally fixed lower asymptote.

- b, c, d:

  Numeric scalars. Free model parameters.

- g:

  Numeric scalar. Asymmetry parameter (free).

## Value

Named numeric vector `c(b=, c=, d=, g=)`.

## See also

Other fixed-a-gradients:
[`grad_inv_gompertz4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_gompertz4_fixed.md),
[`grad_inv_logistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_logistic4_fixed.md),
[`grad_inv_loglogistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_loglogistic4_fixed.md),
[`grad_inv_loglogistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_loglogistic5_fixed.md)
