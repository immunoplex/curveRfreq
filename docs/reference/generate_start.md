# Generate Start Value Lists for NLS Optimisation

Produces a list of lower and upper starting value vectors by spreading
candidate starts across the interior of the parameter bounds.

## Usage

``` r
generate_start(bounds, frac = 0.9)
```

## Arguments

- bounds:

  Named list with elements `lower` and `upper`, both named numeric
  vectors with the same parameter names.

- frac:

  Numeric in `(0, 1)`. The fraction of the total bound width kept inside
  the start interval (default `0.90`).

## Value

A list with `start_lower` and `start_upper`.
