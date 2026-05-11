# Build Inverse, Gradient, and grad_y Closures for a Model

Constructs a list of three closures — `inv`, `grad`, and `grad_y` — that
evaluate the inverse function, its parameter gradient, and its
response-derivative for a given model at a specific observed response
`y`. The closures accept a single named parameter vector `p` (typically
from [`stats::coef()`](https://rdrr.io/r/stats/coef.html)) and are
consumed by
[`propagate_error_analytic()`](https://immunoplex.github.io/curveRfreq/reference/propagate_error_analytic.md).

## Usage

``` r
make_inv_and_grad_fixed(model, y, fixed_a)
```

## Arguments

- model:

  Character. One of `"logistic4"`, `"logistic5"`, `"loglogistic4"`,
  `"loglogistic5"`, `"gompertz4"`.

- y:

  Numeric scalar. The observed response value at which to evaluate the
  inverse and its derivatives.

- fixed_a:

  Numeric scalar or `NULL`. If non-`NULL`, the lower asymptote is
  treated as a fixed external constant.

## Value

A list with three closures:

- `inv(p)`:

  Returns the estimated `x` (concentration) given named parameter vector
  `p`.

- `grad(p)`:

  Returns a named numeric vector of \\\partial x / \partial \theta_i\\.

- `grad_y(p)`:

  Returns a numeric scalar \\\partial x / \partial y\\.

## Design

Two branches exist depending on whether `fixed_a` is supplied:

- `fixed_a` is non-`NULL`:

  Parameter `a` is treated as an externally known constant. The `grad`
  closure returns partials for `b`, `c`, `d` (and `g` for 5-parameter
  models) only. The `_fixed` family of functions is used.

- `fixed_a` is `NULL`:

  Parameter `a` is free and must appear in `p`. The full `grad_*()`
  functions are used, returning partials for `a`, `b`, `c`, `d` (and
  `g`).

## See also

[`propagate_error_analytic()`](https://immunoplex.github.io/curveRfreq/reference/propagate_error_analytic.md),
[`inv_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4.md),
[`grad_logistic4()`](https://immunoplex.github.io/curveRfreq/reference/grad_logistic4.md),
[`inv_logistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/inv_logistic4_fixed.md),
[`grad_inv_logistic4_fixed()`](https://immunoplex.github.io/curveRfreq/reference/grad_inv_logistic4_fixed.md)

## Examples

``` r
# Free-a example
fns <- make_inv_and_grad_fixed("logistic4", y = 5000, fixed_a = NULL)
p <- c(a = 100, b = 1.5, c = 2, d = 20000)
fns$inv(p)
#> [1] 0.3217775
fns$grad(p)
#>             a             b             c             d 
#> -0.0003061224 -1.1188149960  1.0000000000 -0.0001000000 
fns$grad_y(p)
#> [1] 0.0004061224

# Fixed-a example
fns2 <- make_inv_and_grad_fixed("logistic4", y = 5000, fixed_a = 100)
p2 <- c(b = 1.5, c = 2, d = 20000)
fns2$inv(p2)
#> [1] 0.3217775
fns2$grad(p2)
#>             b             c             d 
#> -1.118815e+00  1.000000e+00 -5.025126e-09 
```
