# Inverse of the 4PL Model with Fixed Lower Asymptote

Same algebra as
[`inv_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4.md)
but parameter `a` is supplied as a known constant (`fixed_a`) rather
than estimated from the fit. No domain checking is performed — caller
must ensure `y` is in range.

## Usage

``` r
inv_logistic4_fixed(y, fixed_a, b, c, d)
```

## Arguments

- y:

  Numeric vector. Observed response values.

- fixed_a:

  Numeric scalar. Externally fixed lower asymptote.

- b, c, d:

  Numeric scalars. Model parameters from
  [`stats::coef()`](https://rdrr.io/r/stats/coef.html).

## Value

Numeric vector of estimated `x` values.

## See also

[`inv_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4.md),
[`grad_inv_logistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_logistic4_fixed.md)

Other inverse-functions:
[`inv_gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/inv_gompertz4.md),
[`inv_gompertz4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_gompertz4_fixed.md),
[`inv_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4.md),
[`inv_logistic5()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic5.md),
[`inv_logistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic5_fixed.md),
[`inv_loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic4.md),
[`inv_loglogistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic4_fixed.md),
[`inv_loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic5.md),
[`inv_loglogistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic5_fixed.md)
