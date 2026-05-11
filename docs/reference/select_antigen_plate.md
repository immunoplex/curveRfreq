# Select and Prepare Antigen Plate Data

Filters and prepares standard curves, blanks, samples, and optional MCMC
outputs for a specific antigen plate, identified either by a full
`curve_id` or by its component fields.

## Usage

``` r
select_antigen_plate(
  loaded_data,
  project_id,
  study_accession,
  experiment_accession,
  feature,
  source,
  antigen,
  plate,
  nominal_sample_dilution,
  curve_id_element_order,
  wavelength = WL_NONE,
  antigen_constraints,
  verbose = TRUE
)
```

## Arguments

- loaded_data:

  A list containing datasets: `standards`, `blanks`, `samples`,
  `mcmc_samples`, `mcmc_pred`, and `antigen_constraints`.

- project_id:

  Character. Project identifier (used if `curve_id` is NULL).

- study_accession:

  Character. Study identifier (used if `curve_id` is NULL).

- experiment_accession:

  Character. Experiment identifier.

- feature:

  Character. Feature name (e.g., IgG1).

- source:

  Character. Sample source.

- antigen:

  Character. Antigen name.

- plate:

  Character. Plate identifier.

- nominal_sample_dilution:

  Numeric or character. Nominal dilution.

- curve_id_element_order:

  named vector of order of curve identifier.

- wavelength:

  Optional wavelength filter. Defaults to `WL_NONE`.

- antigen_constraints:

  Data frame or list containing antigen constraint rules.

- verbose:

  Logical. If `TRUE`, prints diagnostic messages.

## Value

A list containing:

- `plate_standard` Filtered standard curve data.

- `plate_blanks` Filtered blank data.

- `plate_samples` Filtered sample data.

- `plate_mcmc_samples` Filtered MCMC samples (if available).

- `plate_mcmc_pred` Filtered and sorted MCMC predictions (if available).

- `antigen_settings` Output of
  [`obtain_lower_constraint()`](https://immunoplex.github.io/curveRfreq/reference/obtain_lower_constraint.md).

- `fixed_a_result` Validated lower asymptote.

- `std_error_blank` Standard error of blank measurements.

- `curve_id` curve id object with attributes of the parts it uses to be
  constructed.

## Details

This function supports two input modes:

- **curve_id mode (preferred):** Provide a full `curve_id` string.

- **component mode:** Provide all individual fields used to construct
  the `curve_id`.

If both are supplied, `curve_id` takes precedence and missing component
arguments are automatically backfilled.

**Deprecated behavior:** Passing a `curve_id` string via
`study_accession` is deprecated and will trigger a warning. Use the
`curve_id` argument instead.

**Matching logic:** Filtering is performed using exact string matching
on the `curve_id` column. Attributes on `curve_id` are ignored, ensuring
compatibility with database inputs.

## Examples

``` r
# Preferred usage
if (FALSE) { # \dontrun{
select_antigen_plate(
  loaded_data = loaded_data,
  curve_id = "proj|study|exp|IgG1|serum|pt|plate1|100",
  antigen_constraints = loaded_data$antigen_constraints
)
} # }

# Component-based usage
if (FALSE) { # \dontrun{
select_antigen_plate(
  loaded_data = loaded_data,
  project_id = "proj",
  study_accession = "study",
  experiment_accession = "exp",
  feature = "IgG1",
  source = "serum",
  antigen = "pt",
  plate = "plate1",
  nominal_sample_dilution = 100,
  antigen_constraints = loaded_data$antigen_constraints
)
} # }
```
