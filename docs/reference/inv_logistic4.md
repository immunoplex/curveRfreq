# Inverse of the 4PL Model

Given observed response `y`, solves analytically for `x`
(concentration): \$\$x = c + b\\\log\\\left(\frac{a - d}{y - d} -
1\right)\$\$

## Usage

``` r
inv_logistic4(y, a, b, c, d, tol = 1e-06)
```

## Arguments

- y:

  Numeric vector. Observed response values.

- a:

  Numeric scalar. Lower asymptote (free parameter from fit).

- b:

  Numeric scalar. Scale/slope parameter.

- c:

  Numeric scalar. Inflection-point location.

- d:

  Numeric scalar. Upper asymptote.

- tol:

  Numeric scalar. Buffer from asymptotes to prevent log-domain errors.
  Default `1e-6`.

## Value

Numeric vector of estimated `x` values, same length as `y`. `NA` for
out-of-range `y`.

## Details

Values of `y` outside the open interval \\(\min(a,d) + \mathrm{tol},\\
\max(a,d) - \mathrm{tol})\\ are undefined and return `NA`.

## See also

[`logistic4()`](https://immunoplex.github.io/curveRfreq/reference/logistic4.md),
[`inv_logistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4_fixed.md)
for the fixed-\\a\\ variant.

Other inverse-functions:
[`inv_gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/inv_gompertz4.md),
[`inv_gompertz4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_gompertz4_fixed.md),
[`inv_logistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4_fixed.md),
[`inv_logistic5()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic5.md),
[`inv_logistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic5_fixed.md),
[`inv_loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic4.md),
[`inv_loglogistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic4_fixed.md),
[`inv_loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic5.md),
[`inv_loglogistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic5_fixed.md)

## Examples

``` r
# Round-trip: forward then inverse recovers x
x <- seq(-1, 5, length.out = 50)
y <- logistic4(x, a = 100, b = 1.5, c = 2, d = 20000)
x_hat <- inv_logistic4(y, a = 100, b = 1.5, c = 2, d = 20000)
all.equal(x, x_hat)
#> [1] TRUE
```
