# Obtain lower asymptote constraints for a given antigen

Returns a named list of lower-asymptote constraint parameters for one
antigen / plate combination. The method used is determined by
`antigen_constraints$l_asy_constraint_method`:

- `"default"`:

  Min = 0, max = observed data maximum.

- `"user_defined"`:

  Uses the values stored in `antigen_constraints`.

- `"range_of_blanks"`:

  Min/max set to the range of blank plate responses.

- `"geometric_mean_of_blanks"`:

  Both min and max set to the geometric mean of blank responses (i.e.
  the parameter is effectively fixed).

## Usage

``` r
obtain_lower_constraint(
  dat,
  antigen,
  study_accession,
  experiment_accession,
  plate,
  plate_blanks,
  antigen_constraints,
  response_col = NULL
)
```

## Arguments

- dat:

  Data frame of standards assay measurements for this antigen and plate
  (standards + samples).

- antigen:

  Character. Antigen identifier.

- study_accession:

  Character. Study accession identifier.

- experiment_accession:

  Character. Experiment accession identifier.

- plate:

  Character. Plate label.

- plate_blanks:

  Data frame of blank (buffer) wells for this plate.

- antigen_constraints:

  Data frame with exactly one row (or more — the first row is used)
  containing constraint columns: `l_asy_constraint_method`,
  `l_asy_min_constraint`, `l_asy_max_constraint`,
  `standard_curve_concentration`, `pcov_threshold`.

- response_col:

  Character or `NULL`. Name of the response column (`"mfi"` for bead
  arrays, `"absorbance"` for ELISA). If `NULL` it is resolved
  automatically via
  [`resolve_response_col()`](https://immunoplex.github.io/curveRfreq/reference/resolve_response_col.md).

## Value

A named list with elements:

- study_accession:

  Passed through.

- experiment_accession:

  Passed through.

- plate:

  Passed through.

- antigen:

  Passed through.

- l_asy_min_constraint:

  Numeric lower bound for the *a* parameter.

- l_asy_max_constraint:

  Numeric upper bound for the *a* parameter.

- l_asy_constraint_method:

  Character. The method used.

- std_error_blank:

  Numeric. Standard error of blank responses.

- standard_curve_concentration:

  Numeric. Concentration of the undiluted standard.

- pcov_threshold:

  Numeric. Percent-CV acceptance threshold.

Returns `NULL` if the constraint method is unrecognized.
