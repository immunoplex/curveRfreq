# Fit Ensemble of NLS Models for One Plate

Fits each model in the ensemble using multi-start Levenberg-Marquardt
NLS, returning the best fit per model (by AIC). Handles convergence
failures gracefully — a model that fails to converge gets
`converged = FALSE` in the output.

## Usage

``` r
fit_ensemble_nls(
  data,
  formulas,
  model_constraints,
  start_lists,
  verbose = FALSE
)
```

## Arguments

- data:

  Data frame of preprocessed standards (already transformed by
  [`curveRcore::preprocess_standards()`](https://immunoplex.github.io/curveRcore/reference/preprocess_standards.html)).

- formulas:

  Named list of formula objects from
  [`curveRcore::build_nls_formulas()`](https://immunoplex.github.io/curveRcore/reference/build_nls_formulas.html).

- model_constraints:

  Named list from
  [`compute_model_constraints()`](https://immunoplex.github.io/curveRfreq/reference/compute_model_constraints.md).

- start_lists:

  Named list from
  [`generate_start_lists()`](https://immunoplex.github.io/curveRfreq/reference/generate_start_lists.md).

- verbose:

  Logical.

## Value

A named list (one element per model) of lists with:

- model_name:

  Character.

- converged:

  Logical.

- fit:

  The `nlsLM` object, or NULL.

- aic:

  Numeric (Inf if failed).

- data:

  The data used for fitting.
