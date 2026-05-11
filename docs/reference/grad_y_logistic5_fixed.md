# \\\partial x / \partial y\\ for the 5PL Inverse with Fixed \\a\\

\\\partial x / \partial y\\ for the 5PL Inverse with Fixed \\a\\

## Usage

``` r
grad_y_logistic5_fixed(y, fixed_a, b, d, g)
```

## Arguments

- y:

  Numeric scalar. Observed response.

- fixed_a:

  Numeric scalar. Externally fixed lower asymptote.

- b:

  Numeric scalar. Slope parameter.

- d:

  Numeric scalar. Upper asymptote.

- g:

  Numeric scalar. Asymmetry parameter.

## Value

Numeric scalar \\\partial x / \partial y\\.

## See also

Other fixed-a-grad-y:
[`grad_y_gompertz4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_y_gompertz4_fixed.md),
[`grad_y_logistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_y_logistic4_fixed.md),
[`grad_y_loglogistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_y_loglogistic4_fixed.md),
[`grad_y_loglogistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_y_loglogistic5_fixed.md)
