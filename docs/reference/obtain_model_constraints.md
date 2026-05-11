# Compute Model Constraints for All Candidate Models

Builds constraint lists (lower and upper bounds) for each of the five
candidate sigmoid models. A shared `constraint_profile` is computed once
from the data and reused across all models.

## Usage

``` r
obtain_model_constraints(
  data,
  formulas,
  response_variable,
  independent_variable,
  is_log_response,
  is_log_concentration,
  antigen_settings,
  max_response,
  min_response,
  verbose = TRUE
)
```

## Arguments

- data:

  Data frame of preprocessed standard curve data.

- formulas:

  Named list of model formulae
  (`logistic5, loglogistic5, logistic4, loglogistic4, gompertz4`).

- response_variable:

  Character. Response column name.

- independent_variable:

  Character. Concentration column name.

- is_log_response:

  Logical. Is the response log10-transformed?

- is_log_concentration:

  Logical. Is concentration log10-transformed?

- antigen_settings:

  Named list from
  [`obtain_lower_constraint`](https://immunoplex.github.io/curveRfreq/reference/obtain_lower_constraint.md).

- max_response:

  Numeric. Maximum observed response.

- min_response:

  Numeric. Minimum observed response.

- verbose:

  Logical (default `TRUE`).

## Value

Named list of constraint lists (one per model), with an
`"constraint_profile"` attribute attached for downstream use.
