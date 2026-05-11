# First Derivative of the loglogistic4 Model

Computes \\dy/dx\\ analytically for
[`loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic4.md):
\$\$\frac{dy}{dx} = -\frac{b\\(d - a)\\r}{x\\(1 + r)^2} \quad\text{where
} r = (x/c)^b\$\$

## Usage

``` r
dydxloglogistic4(x, a, b, c, d)
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

## Value

Numeric vector of first-derivative values.

## See also

[`loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic4.md),
[`d2xloglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/d2xloglogistic4.md)

Other derivatives:
[`d2xgompertz4()`](https://immunoplex.github.io/curveRfreq/reference/d2xgompertz4.md),
[`d2xlogistic4()`](https://immunoplex.github.io/curveRfreq/reference/d2xlogistic4.md),
[`d2xlogistic5()`](https://immunoplex.github.io/curveRfreq/reference/d2xlogistic5.md),
[`d2xloglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/d2xloglogistic4.md),
[`d2xloglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/d2xloglogistic5.md),
[`dydxgompertz4()`](https://immunoplex.github.io/curveRfreq/reference/dydxgompertz4.md),
[`dydxlogistic4()`](https://immunoplex.github.io/curveRfreq/reference/dydxlogistic4.md),
[`dydxlogistic5()`](https://immunoplex.github.io/curveRfreq/reference/dydxlogistic5.md),
[`dydxloglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/dydxloglogistic5.md)
