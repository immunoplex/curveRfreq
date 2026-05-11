# Four-Parameter Dose–Response (loglogistic4) Forward Function

An alternative 4-parameter logistic parameterisation where the
inflection point `c` is on the concentration scale (not log-scale) and
the slope `b` acts as a Hill coefficient: \$\$y = a + \frac{d - a}{1 +
(x / c)^b}\$\$

## Usage

``` r
loglogistic4(x, a, b, c, d)
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

Numeric vector of predicted response values.

## See also

[`logistic4()`](https://immunoplex.github.io/curveRfreq/reference/logistic4.md),
[`loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic5.md),
[`inv_loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_loglogistic4.md),
[`dydxloglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/dydxloglogistic4.md)

Other forward-models:
[`gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/gompertz4.md),
[`logistic4()`](https://immunoplex.github.io/curveRfreq/reference/logistic4.md),
[`logistic5()`](https://immunoplex.github.io/curveRfreq/reference/logistic5.md),
[`loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic5.md)

## Examples

``` r
x <- seq(0.01, 100, length.out = 200)
y <- loglogistic4(x, a = 50, b = -2, c = 10, d = 15000)
plot(x, y, type = "l", log = "x", main = "loglogistic4 Dose-Response")

```
