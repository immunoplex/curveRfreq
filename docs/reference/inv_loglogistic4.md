# Inverse of the loglogistic4 Model

Solves for `x` given response `y` under the
[`loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic4.md)
parameterisation: \$\$x = c\\\left(\frac{d - a}{y - a} -
1\right)^{1/b}\$\$

## Usage

``` r
inv_loglogistic4(y, a, b, c, d)
```

## Arguments

- y:

  Numeric vector. Observed response values.

- a, b, c, d:

  Numeric scalars. loglogistic4 model parameters.

## Value

Numeric vector of estimated `x` values.

## See also

[`loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic4.md),
[`inv_loglogistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic4_fixed.md)

Other inverse-functions:
[`inv_gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/inv_gompertz4.md),
[`inv_gompertz4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_gompertz4_fixed.md),
[`inv_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4.md),
[`inv_logistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4_fixed.md),
[`inv_logistic5()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic5.md),
[`inv_logistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic5_fixed.md),
[`inv_loglogistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic4_fixed.md),
[`inv_loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic5.md),
[`inv_loglogistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic5_fixed.md)
