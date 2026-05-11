# Inverse of the 5PL Model

Solves for `x` given response `y` under
[`logistic5()`](https://immunoplex.github.io/curveRfreq/reference/logistic5.md):
\$\$x = c + b\\\log\\\left(\left(\frac{a - d}{y - d}\right)^{1/g} -
1\right)\$\$

## Usage

``` r
inv_logistic5(y, a, b, c, d, g)
```

## Arguments

- y:

  Numeric vector. Observed response values.

- a, b, c, d:

  Numeric scalars. 5PL model parameters.

- g:

  Numeric scalar. Asymmetry parameter.

## Value

Numeric vector of estimated `x` values.

## See also

[`logistic5()`](https://immunoplex.github.io/curveRfreq/reference/logistic5.md),
[`inv_logistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic5_fixed.md)

Other inverse-functions:
[`inv_gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/inv_gompertz4.md),
[`inv_gompertz4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_gompertz4_fixed.md),
[`inv_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4.md),
[`inv_logistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4_fixed.md),
[`inv_logistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic5_fixed.md),
[`inv_loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic4.md),
[`inv_loglogistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic4_fixed.md),
[`inv_loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic5.md),
[`inv_loglogistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic5_fixed.md)
