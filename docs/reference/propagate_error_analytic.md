# Analytic error propagation for inverse sigmoid models

Computes the inverse-predicted concentration (`x_est`) and its
propagated uncertainty (`se_x`) using the delta method for a fitted
nonlinear model.

## Usage

``` r
propagate_error_analytic(model, fit, y, se_y = 0, fixed_a, verbose = TRUE)
```

## Arguments

- model:

  Character string specifying the model form. One of `"logistic4"`,
  `"loglogistic4"`, `"gompertz4"`, `"logistic5"`, `"loglogistic5"`.

- fit:

  A fitted `nlsLM` model object (from `minpack.lm`).

- y:

  Numeric scalar. Observed response value.

- se_y:

  Numeric scalar. Standard error of the response `y`. Defaults to `0` if
  unknown.

- fixed_a:

  Optional numeric scalar. Fixed lower asymptote. If supplied, `a` is
  treated as fixed and excluded from uncertainty propagation.

- verbose:

  Logical. If `TRUE`, prints diagnostic information.

## Value

A list with components:

- x_est:

  Inverse-predicted concentration

- se_x:

  Standard error of `x_est`

- var_x:

  Variance of `x_est`

- grad_theta:

  Gradient w.r.t. model parameters

- grad_y:

  Gradient w.r.t. response `y`

## Details

Supports multiple sigmoid model forms (4PL, 5PL, Gompertz variants), and
optionally handles the case where the lower asymptote (`a`) is fixed.

Uses the delta method: \$\$Var(x) = \nabla\_\theta^T \Sigma
\nabla\_\theta + (\partial x / \partial y)^2 \cdot Var(y)\$\$

If `fixed_a` is supplied, gradients are computed using reduced parameter
space.
