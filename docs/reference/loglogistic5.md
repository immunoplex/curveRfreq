# Five-Parameter Dose–Response (loglogistic5) Forward Function

A five-parameter generalised logistic (Richards) curve in dose–response
form with asymmetry parameter `g`: \$\$y = a + (d - a)\\\bigl(1 +
g\\\exp(-b\\(x - c))\bigr)^{-1/g}\$\$

## Usage

``` r
loglogistic5(x, a, b, c, d, g)
```

## Arguments

- x:

  Numeric vector. Independent variable (typically log-concentration).

- a:

  Numeric scalar. Lower (left) asymptote — baseline response at zero
  concentration or infinite dilution.

- b:

  Numeric scalar. Scale (slope) parameter controlling steepness of the
  transition region. Sign determines curve direction.

- c:

  Numeric scalar. Inflection-point location on the x-axis (midpoint).

- d:

  Numeric scalar. Upper (right) asymptote — maximum response at
  saturating concentration.

- g:

  Numeric scalar. Asymmetry (Richards) parameter.

## Value

Numeric vector of predicted response values.

## Details

When \\g = 1\\ this reduces to a standard 4PL dose–response.

## See also

[`loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic4.md),
[`logistic5()`](https://immunoplex.github.io/curveRfreq/reference/logistic5.md),
[`inv_loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic5.md),
[`dydxloglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/dydxloglogistic5.md)

Other forward-models:
[`gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/gompertz4.md),
[`logistic4()`](https://immunoplex.github.io/curveRfreq/reference/logistic4.md),
[`logistic5()`](https://immunoplex.github.io/curveRfreq/reference/logistic5.md),
[`loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic4.md)
