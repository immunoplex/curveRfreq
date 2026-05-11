# Parse Colon-Separated Curve Identifiers into Structured Columns

Splits a colon-delimited `curve_id` column into multiple fields based on
a user-supplied positional `order`, with optional validation and the
ability to retain only a subset of parsed fields.

## Usage

``` r
parse_curve_id(
  data,
  curve_col = "curve_id",
  order,
  schema = NULL,
  keep = NULL,
  sep = ":",
  validate = TRUE
)
```

## Arguments

- data:

  A data frame containing a column with encoded curve identifiers.

- curve_col:

  A character string specifying the name of the column containing the
  curve identifiers. Default is `"curve_id"`.

- order:

  A character vector defining the positional order of fields in the
  `curve_id`. The length of `order` must match the number of elements in
  each `curve_id` string.

- schema:

  Optional named list with elements:

  required

  :   Character vector of field names that must be present and
      non-missing.

  optional

  :   Character vector of optional field names.

  If provided and `validate = TRUE`, required fields are checked for
  missing or empty values.

- keep:

  Optional character vector specifying a subset of parsed fields to
  retain in the output. If `NULL` (default), all parsed fields are kept.

- sep:

  A character string used as the delimiter in `curve_id`. Default is
  `":"`.

- validate:

  Logical; if `TRUE`, performs validation checks:

  - Ensures each `curve_id` has the same number of fields as `order`.

  - If `schema` is provided, ensures required fields are present and
    non-missing.

  Default is `TRUE`.

## Value

A data frame containing the original `data` columns and additional
parsed columns defined by `order` (or a subset if `keep` is specified).

## Details

This function assumes that `curve_id` values are constructed using a
consistent positional encoding (e.g., via `paste(..., sep = ":")`),
where each position corresponds to a known field.

The mapping from `curve_id` to columns is purely positional: the i-th
element of the split string is assigned to the i-th name in `order`.
Therefore, correct parsing depends on consistent construction of
`curve_id` values.

The `keep` argument allows users to retain only specific parsed fields,
which is useful for downstream analysis, plotting, or modeling
workflows.

## Examples

``` r
df <- data.frame(
  curve_id = c(
    "proj1:study1:exp1:IgG1:sample:prn:plate1:1000:450",
    "proj2:study2:exp2:IgG2:control:pt:plate2:2000:__none__"
  )
)

order <- c(
  "project_id",
  "study_accession",
  "experiment_accession",
  "feature",
  "source",
  "antigen",
  "plate",
  "nominal_sample_dilution",
  "wavelength"
)

schema <- list(
  required = c("feature", "antigen", "plate", "wavelength"),
  optional = c("project_id")
)

# Parse all fields
parsed_all <- parse_curve_id(
  data = df,
  order = order,
  schema = schema
)

# Keep only selected fields
parsed_subset <- parse_curve_id(
  data = df,
  order = order,
  keep = c("feature", "antigen", "plate")
)
```
