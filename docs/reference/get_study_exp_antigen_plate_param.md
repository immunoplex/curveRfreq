# Extract Undiluted Standard Curve Concentration from Constraint List

Convenience accessor that retrieves the nominal concentration of the
undiluted standard from the constraint list produced by
[`obtain_lower_constraint`](https://immunoplex.github.io/curveRfreq/reference/obtain_lower_constraint.md).

## Usage

``` r
get_study_exp_antigen_plate_param(l_asy_constraints)
```

## Arguments

- l_asy_constraints:

  Named list returned by
  [`obtain_lower_constraint`](https://immunoplex.github.io/curveRfreq/reference/obtain_lower_constraint.md).

## Value

Numeric scalar: the undiluted standard curve concentration (e.g.
`10000`).
