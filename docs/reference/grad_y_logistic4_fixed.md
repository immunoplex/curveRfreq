# \\\partial x / \partial y\\ for the 4PL Inverse with Fixed \\a\\

Computes the derivative of the inverse 4PL with respect to the observed
response `y` when the lower asymptote is externally fixed. Used for the
measurement-uncertainty component of error propagation.

## Usage

``` r
grad_y_logistic4_fixed(y, fixed_a, b, d)
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

## Value

Numeric scalar \\\partial x / \partial y\\.

## See also

Other fixed-a-grad-y:
[`grad_y_gompertz4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_y_gompertz4_fixed.md),
[`grad_y_logistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_y_logistic5_fixed.md),
[`grad_y_loglogistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_y_loglogistic4_fixed.md),
[`grad_y_loglogistic5_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_y_loglogistic5_fixed.md)
