# Construct a Standardized curve_id String

Builds a colon-separated (`:`) `curve_id` string from named components,
enforcing a consistent field order defined by `order`. This ensures that
`curve_id` values are reproducible and schema-compliant regardless of
the order in which arguments are supplied.

## Usage

``` r
make_curve_id_string(..., sep = ":", order)
```

## Arguments

- ...:

  Named components used to construct the `curve_id`. Names must match
  the elements of `order`. All fields in `order` are required unless
  handled externally (e.g., optional defaults).

- sep:

  Character separator used to join components. Defaults to `":"`.

- order:

  Character vector defining the required fields and their order in the
  resulting `curve_id`.

## Value

A single character string representing the standardized `curve_id`.

## Details

All required fields defined in `order` must be provided as named
arguments. The function will reorder inputs internally to match the
specified schema.

This function does not attach attributes to the returned string.
Instead, parsing should be performed using
[`parse_curve_id`](https://immunoplex.github.io/curveRfreq/reference/parse_curve_id.md),
which relies on an explicit schema (`order`) for robustness and
reproducibility.

This design avoids issues with attribute loss during data manipulation,
storage, or serialization (e.g., database writes, CSV export).

## See also

[`parse_curve_id`](https://immunoplex.github.io/curveRfreq/reference/parse_curve_id.md)

## Examples

``` r
make_curve_id_string(
  project_id = "17",
  study_accession = "MADI_01",
  experiment_accession = "IgG1",
  feature = "IgG1",
  source = "Standard",
  antigen = "pt",
  plate = "plate1",
  nominal_sample_dilution = "1x",
  wavelength = "450",
  order =  c("project_id", "study_accession", "experiment_accession", "feature", "source",
           "antigen", "plate", "nominal_sample_dilution", "wavelength")
)
#> [1] "17:MADI_01:IgG1:IgG1:Standard:pt:plate1:1x:450"

# Order of arguments does not matter (must have order argument)
make_curve_id_string(
  antigen = "pt",
  plate = "plate1",
  feature = "IgG1",
  source = "Standard",
  project_id = "17",
  study_accession = "MADI_01",
  experiment_accession = "IgG1",
  nominal_sample_dilution = "1x",
  wavelength = "450",
  order =  c("project_id", "study_accession", "experiment_accession", "feature", "source",
           "antigen", "plate", "nominal_sample_dilution", "wavelength")
)
#> [1] "17:MADI_01:IgG1:IgG1:Standard:pt:plate1:1x:450"
```
