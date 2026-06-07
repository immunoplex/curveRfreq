# Generate Multi-Start Lists for NLS Fitting

Creates candidate starting values evenly spread across the parameter
bounds. Low-signal data gets more starts with the slope biased small.

## Usage

``` r
generate_start_lists(model_constraints)
```

## Arguments

- model_constraints:

  Output of
  [`compute_model_constraints()`](https://immunoplex.github.io/curveRfreq/reference/compute_model_constraints.md).

## Value

Named list of lists: for each model, a list of named numeric vectors
suitable as `start` for `nlsLM()`.
