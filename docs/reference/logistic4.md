# Four-Parameter Logistic (4PL) Function

Computes the response \\y\\ for a four-parameter logistic curve: \$\$y =
d + \frac{a - d}{1 + \exp\\\left(\frac{x - c}{b}\right)}\$\$

## Usage

``` r
logistic4(x, a, b, c, d)
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

Numeric vector of predicted response values, same length as `x`.

## Details

The classical 4PL model is the workhorse of immunoassay quantitation
(ELISA, Luminex bead arrays). It is symmetric about its inflection
point; for asymmetric data see
[`logistic5()`](https://immunoplex.github.io/curveRfreq/reference/logistic5.md).

## See also

[`logistic5()`](https://immunoplex.github.io/curveRfreq/reference/logistic5.md),
[`loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic4.md),
[`gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/gompertz4.md)
for alternative parameterisations;
[`inv_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4.md)
for the inverse;
[`dydxlogistic4()`](https://immunoplex.github.io/curveRfreq/reference/dydxlogistic4.md)
for the first derivative.

Other forward-models:
[`gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/gompertz4.md),
[`logistic5()`](https://immunoplex.github.io/curveRfreq/reference/logistic5.md),
[`loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic4.md),
[`loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic5.md)

## Examples

``` r
x <- seq(-2, 6, length.out = 200)
y <- logistic4(x, a = 100, b = 1.5, c = 2, d = 20000)
plot(x, y, type = "l", main = "4PL Standard Curve")

```
