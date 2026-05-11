# Second Derivative of the loglogistic5 Model (Numerical)

Approximates \\d^2y/dx^2\\ for
[`loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic5.md)
via central difference.

## Usage

``` r
d2xloglogistic5(x, a, b, c, d, g, h = 1e-05)
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

- h:

  Numeric scalar. Finite-difference step size. Default `1e-5`.

## Value

Numeric vector of second-derivative values.

## See also

Other derivatives:
[`d2xgompertz4()`](https://immunoplex.github.io/curveRfreq/reference/d2xgompertz4.md),
[`d2xlogistic4()`](https://immunoplex.github.io/curveRfreq/reference/d2xlogistic4.md),
[`d2xlogistic5()`](https://immunoplex.github.io/curveRfreq/reference/d2xlogistic5.md),
[`d2xloglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/d2xloglogistic4.md),
[`dydxgompertz4()`](https://immunoplex.github.io/curveRfreq/reference/dydxgompertz4.md),
[`dydxlogistic4()`](https://immunoplex.github.io/curveRfreq/reference/dydxlogistic4.md),
[`dydxlogistic5()`](https://immunoplex.github.io/curveRfreq/reference/dydxlogistic5.md),
[`dydxloglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/dydxloglogistic4.md),
[`dydxloglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/dydxloglogistic5.md)
