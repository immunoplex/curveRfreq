# Compute Parameter Bounds for All Candidate Models

Inspects the response range to classify the assay scale
(high/medium/low), then builds adaptive bounds for a, b, c, d (and g for
5-param models).

## Usage

``` r
compute_model_constraints(
  data,
  formulas,
  response_variable,
  independent_variable,
  is_log_response,
  fixed_a = NULL,
  verbose = FALSE
)
```

## Arguments

- data:

  Data frame with response and concentration columns.

- formulas:

  Named list of formula objects.

- response_variable:

  Character. Response column name.

- independent_variable:

  Character. Concentration column name.

- is_log_response:

  Logical. Is the response log10-transformed?

- fixed_a:

  Numeric or NULL. If non-NULL, `a` is fixed and not bounded. If NULL,
  adaptive bounds for `a` are computed.

- verbose:

  Logical.

## Value

Named list of lists, each with `lower` and `upper` named vectors. Has
attribute `"constraint_profile"`.

## Details

When `fixed_a` is NULL, the lower asymptote `a` is free and its bounds
are set adaptively from the data. When `fixed_a` is non-NULL, `a` does
not appear in the formulas (it's baked in), so no bounds are needed.
