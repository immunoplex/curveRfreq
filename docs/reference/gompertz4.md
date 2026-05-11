# Four-Parameter Gompertz Forward Function

Computes the response for a Gompertz growth/saturation curve: \$\$y =
a + (d - a)\\\exp\\\bigl(-\exp(-b\\(x - c))\bigr)\$\$

## Usage

``` r
gompertz4(x, a, b, c, d)
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

## Details

The Gompertz is intrinsically asymmetric — unlike the 4PL it does not
require a fifth parameter for asymmetry.

## See also

[`logistic4()`](https://immunoplex.github.io/curveRfreq/reference/logistic4.md),
[`loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic4.md),
[`inv_gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/inv_gompertz4.md),
[`dydxgompertz4()`](https://immunoplex.github.io/curveRfreq/reference/dydxgompertz4.md)

Other forward-models:
[`logistic4()`](https://immunoplex.github.io/curveRfreq/reference/logistic4.md),
[`logistic5()`](https://immunoplex.github.io/curveRfreq/reference/logistic5.md),
[`loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic4.md),
[`loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic5.md)

## Examples

``` r
x <- seq(-2, 8, length.out = 200)
y <- gompertz4(x, a = 50, b = 1, c = 3, d = 15000)
plot(x, y, type = "l", main = "Gompertz Curve")

```
