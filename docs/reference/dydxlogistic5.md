# First Derivative of the 5PL Model

Computes \\dy/dx\\ analytically for
[`logistic5()`](https://immunoplex.github.io/curveRfreq/reference/logistic5.md):
\$\$\frac{dy}{dx} = -\frac{g\\(a - d)\\u}{b\\(1 + u)^{g+1}}\$\$

## Usage

``` r
dydxlogistic5(x, a, b, c, d, g)
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

  Numeric scalar. Asymmetry parameter. Values \\\> 1\\ skew the curve
  toward the upper asymptote; values \\\< 1\\ skew it toward the lower
  asymptote.

## Value

Numeric vector of first-derivative values.

## See also

[`logistic5()`](https://immunoplex.github.io/curveRfreq/reference/logistic5.md),
[`d2xlogistic5()`](https://immunoplex.github.io/curveRfreq/reference/d2xlogistic5.md)

Other derivatives:
[`d2xgompertz4()`](https://immunoplex.github.io/curveRfreq/reference/d2xgompertz4.md),
[`d2xlogistic4()`](https://immunoplex.github.io/curveRfreq/reference/d2xlogistic4.md),
[`d2xlogistic5()`](https://immunoplex.github.io/curveRfreq/reference/d2xlogistic5.md),
[`d2xloglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/d2xloglogistic4.md),
[`d2xloglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/d2xloglogistic5.md),
[`dydxgompertz4()`](https://immunoplex.github.io/curveRfreq/reference/dydxgompertz4.md),
[`dydxlogistic4()`](https://immunoplex.github.io/curveRfreq/reference/dydxlogistic4.md),
[`dydxloglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/dydxloglogistic4.md),
[`dydxloglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/dydxloglogistic5.md)
