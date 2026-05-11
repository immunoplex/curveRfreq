# Build Multi-Start Lists for Each Model

Generates `n_starts` candidate start-value lists for each model by
distributing points evenly across the parameter bounds (Latin
hypercube-style). For low-signal data, extra starts are added and the
slope parameter is biased toward smaller values.

## Usage

``` r
make_start_lists(
  model_constraints,
  quants = c(low = 0.2, mid = 0.5, high = 0.8)
)
```

## Arguments

- model_constraints:

  Named list of model constraint objects as returned by
  [`obtain_model_constraints`](https://immunoplex.github.io/curveRfreq/reference/obtain_model_constraints.md).

- quants:

  Named numeric vector of quantiles (unused, retained for
  compatibility).

## Value

Named list of lists: for each model, a list of named numeric vectors
suitable for passing as `start` to `nlsLM`.
