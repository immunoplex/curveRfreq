# Determine Whether the Lower Asymptote Should Be Fixed

Tests whether the min and max constraints for *a* are identical. If they
are, the parameter is effectively fixed and its value is returned;
otherwise `NULL` is returned to indicate a free parameter.

## Usage

``` r
resolve_fixed_lower_asymptote(l_asy_constraints)
```

## Arguments

- l_asy_constraints:

  Named list (output of
  [`obtain_lower_constraint`](https://immunoplex.github.io/curveRfreq/reference/obtain_lower_constraint.md)).

## Value

Numeric scalar (the fixed value) if min == max; `NULL` if the parameter
should be estimated freely.
