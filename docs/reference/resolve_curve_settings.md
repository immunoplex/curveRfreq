# Resolve curve settings for an antigen on a plate

Computes modeling inputs and constraints for a selected antigen on a
plate, including lower asymptote handling and blank variability
estimates.

## Usage

``` r
resolve_curve_settings(
  loaded_data,
  antigen_constraints,
  wavelength = WL_NONE,
  verbose = TRUE
)
```

## Arguments

- loaded_data:

  A list containing plate-level datasets with the following elements:

  - `standards` Standard curve data.

  - `blanks` Blank/control data.

  - `samples` Sample data.

  - `curve_id_lookup` a single row containing the plate-level curve
    identifier including antigen, study accession and
    experiment_accession and the id number.

- antigen_constraints:

  Data frame or list containing antigen-specific constraint rules used
  in
  [`obtain_lower_constraint()`](https://immunoplex.github.io/curveRfreq/reference/obtain_lower_constraint.md).

- wavelength:

  Optional wavelength filter. If provided, all datasets are filtered to
  the matching normalized wavelength. Defaults to `WL_NONE`.

- verbose:

  Logical. If `TRUE`, prints diagnostic messages including row counts.

## Value

A list containing:

- `plate_standard` Filtered standard curve data.

- `plate_blanks` Filtered blank data.

- `plate_samples` Filtered sample data.

- `antigen_settings` Output of
  [`obtain_lower_constraint()`](https://immunoplex.github.io/curveRfreq/reference/obtain_lower_constraint.md).

- `fixed_a_result` Validated lower asymptote result.

- `std_error_blank` Estimated standard error of blank measurements.

## Details

This function assumes that plate-level datasets have already been
filtered and derives additional settings needed for curve fitting.

**Wavelength filtering:** If `wavelength` is provided (not `WL_NONE`),
all input datasets are filtered to rows matching the normalized
wavelength using
[`normalize_wavelength()`](https://immunoplex.github.io/curveRfreq/reference/normalize_wavelength.md).

**Assumptions:**

- Input data has already been subset to a single curve/plate context.

- The `standards` dataset must contain at least one row after filtering,
  otherwise an error is thrown.

## Examples

``` r
if (FALSE) { # \dontrun{
settings <- resolve_curve_settings(
  loaded_data = plate_data,
  antigen_constraints = constraints
)
} # }
```
