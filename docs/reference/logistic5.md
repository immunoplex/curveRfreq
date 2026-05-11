# Five-Parameter Logistic (5PL) Function

Computes the response for a five-parameter logistic curve with asymmetry
parameter `g`: \$\$y = d + \frac{a - d}{\left(1 + \exp\\\left(\frac{x -
c}{b}\right)\right)^{\\g}}\$\$

## Usage

``` r
logistic5(x, a, b, c, d, g)
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

Numeric vector of predicted response values.

## Details

When \\g = 1\\ this reduces to
[`logistic4()`](https://immunoplex.github.io/curveRfreq/reference/logistic4.md).

## See also

[`logistic4()`](https://immunoplex.github.io/curveRfreq/reference/logistic4.md),
[`inv_logistic5()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic5.md),
[`dydxlogistic5()`](https://immunoplex.github.io/curveRfreq/reference/dydxlogistic5.md)

Other forward-models:
[`gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/gompertz4.md),
[`logistic4()`](https://immunoplex.github.io/curveRfreq/reference/logistic4.md),
[`loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic4.md),
[`loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic5.md)

## Examples

``` r
x <- seq(-2, 6, length.out = 200)
y <- logistic5(x, a = 80, b = 1.2, c = 2.5, d = 18000, g = 1.3)
plot(x, y, type = "l", main = "5PL Standard Curve")

```
