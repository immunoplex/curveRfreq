# First Derivative of the 4PL Model

Computes \\dy/dx\\ analytically for
[`logistic4()`](https://immunoplex.github.io/curveRfreq/reference/logistic4.md):
\$\$\frac{dy}{dx} = -\frac{(a - d)\\u}{b\\(1 + u)^2} \quad\text{where }
u = \exp\\\left(\frac{x - c}{b}\right)\$\$

## Usage

``` r
dydxlogistic4(x, a, b, c, d)
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

Numeric vector of first-derivative values, same length as `x`.

## See also

[`logistic4()`](https://immunoplex.github.io/curveRfreq/reference/logistic4.md),
[`d2xlogistic4()`](https://immunoplex.github.io/curveRfreq/reference/d2xlogistic4.md)

Other derivatives:
[`d2xgompertz4()`](https://immunoplex.github.io/curveRfreq/reference/d2xgompertz4.md),
[`d2xlogistic4()`](https://immunoplex.github.io/curveRfreq/reference/d2xlogistic4.md),
[`d2xlogistic5()`](https://immunoplex.github.io/curveRfreq/reference/d2xlogistic5.md),
[`d2xloglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/d2xloglogistic4.md),
[`d2xloglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/d2xloglogistic5.md),
[`dydxgompertz4()`](https://immunoplex.github.io/curveRfreq/reference/dydxgompertz4.md),
[`dydxlogistic5()`](https://immunoplex.github.io/curveRfreq/reference/dydxlogistic5.md),
[`dydxloglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/dydxloglogistic4.md),
[`dydxloglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/dydxloglogistic5.md)

## Examples

``` r
x <- seq(-2, 6, length.out = 200)
slope <- dydxlogistic4(x, a = 100, b = 1.5, c = 2, d = 20000)
plot(x, slope, type = "l", main = "4PL dy/dx")

```
