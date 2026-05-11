# Compute Median Assay SE for Each Antigen/Feature Across All Plates

For each unique combination of study_accession, experiment_accession,
source, antigen, and feature, computes the standard error of assay
response at every dilution level across all plates, then returns the
median of those per-dilution SEs. This pooled median SE can be reused
for error propagation on each individual plate.

## Usage

``` r
compute_antigen_se_table(
  standards_data,
  curve_id_lookup,
  curve_col = "curve_id",
  response_col = "mfi",
  dilution_col = "dilution",
  plate_col = "plate_nom",
  grouping_cols = c("study_accession", "experiment_accession", "source_nom", "antigen",
    "feature"),
  min_reps = 2,
  verbose = FALSE
)
```

## Arguments

- standards_data:

  data.frame containing all standard curve data

- curve_id_lookup:

  the string of the names in order of the elements in the curve_id

- curve_col:

  name of the curve identifier colummn (default = "curve_id")

- response_col:

  name of the response column (e.g., "mfi")

- dilution_col:

  name of the dilution column (default = "dilution")

- plate_col:

  name of the plate identifier column (default = "plate_nom")

- grouping_cols:

  character vector of columns defining the grouping (default =
  c("study_accession", "experiment_accession", "source", "antigen",
  "feature"))

- min_reps:

  minimum number of non-missing plate replicates required at a dilution
  level for that dilution's SE to be included (default = 2)

- verbose:

  logical; if TRUE emit progress messages (default = FALSE)

## Value

A data.frame with one row per unique grouping containing:

- grouping_cols:

  the grouping columns

- median_se:

  median SE across all qualifying dilution levels

- n_dilutions_used:

  number of dilution levels with \>= min_reps non-missing observations
  that contributed to the median

- n_plates:

  number of distinct plates in the group

- total_obs:

  total number of non-missing response observations used
