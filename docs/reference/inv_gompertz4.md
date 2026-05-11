# Inverse of the Gompertz Model

Solves for `x` given response `y` under
[`gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/gompertz4.md):
\$\$x = c - \frac{1}{b}\\\log\\\left(-\log\frac{y - a}{d - a}\right)\$\$

## Usage

``` r
inv_gompertz4(y, a, b, c, d)
```

## Arguments

- y:

  Numeric vector. Observed response values.

- a, b, c, d:

  Numeric scalars. Gompertz model parameters.

## Value

Numeric vector of estimated `x` values.

## See also

[`gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/gompertz4.md),
[`inv_gompertz4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_gompertz4_fixed.md)

Other inverse-functions:
[`inv_gompertz4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_gompertz4_fixed.md),
[`inv_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4.md),
[`inv_logistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4_fixed.md),
[`inv_logistic5()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic5.md),
[`inv_logistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic5_fixed.md),
[`inv_loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic4.md),
[`inv_loglogistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic4_fixed.md),
[`inv_loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic5.md),
[`inv_loglogistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic5_fixed.md)
