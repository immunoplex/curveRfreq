# Inverse of the loglogistic5 Model

Solves for `x` given response `y` under
[`loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic5.md):
\$\$x = c - \frac{1}{b}\left(\log\\\left(\left(\frac{y - a}{d - a}
\right)^{-g} - 1\right) - \log g\right)\$\$

## Usage

``` r
inv_loglogistic5(y, a, b, c, d, g)
```

## Arguments

- y:

  Numeric vector. Observed response values.

- a, b, c, d:

  Numeric scalars. loglogistic5 model parameters.

- g:

  Numeric scalar. Asymmetry parameter.

## Value

Numeric vector of estimated `x` values.

## See also

[`loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic5.md),
[`inv_loglogistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic5_fixed.md)

Other inverse-functions:
[`inv_gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/inv_gompertz4.md),
[`inv_gompertz4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_gompertz4_fixed.md),
[`inv_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4.md),
[`inv_logistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4_fixed.md),
[`inv_logistic5()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic5.md),
[`inv_logistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic5_fixed.md),
[`inv_loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic4.md),
[`inv_loglogistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic4_fixed.md),
[`inv_loglogistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic5_fixed.md)
