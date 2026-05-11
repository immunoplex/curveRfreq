# Bead-based immunoassay example dataset

A named list containing simulated multi-plate bead-based immunoassay
data, including standard curves, blanks, and patient samples across
multiple time points. The dataset spans six assay plates across two
analytes (alpha and beta), each with three replicate plates.

## Usage

``` r
bead_assay_example
```

## Format

A named list with six elements:

- standards:

  A data frame with 60 rows and 8 columns containing standard curve
  data:

  curve_id

  :   Integer. Plate identifier (1–6).

  stype

  :   Character. Sample type; `"S"` for standard.

  sampleid

  :   Character. Standard well identifier (e.g. `"STD_01"`).

  well

  :   Character. Plate well position.

  dilution

  :   Numeric. Known analyte concentration of the standard.

  mfi

  :   Numeric. Median fluorescence intensity.

  assay_response_variable

  :   Character. Name of the response variable (`"mfi"`).

  assay_independent_variable

  :   Character. Name of the independent variable (`"concentration"`).

- blanks:

  A data frame with 24 rows and 7 columns containing blank well
  measurements (4 blanks per plate):

  curve_id

  :   Integer. Plate identifier (1–6).

  stype

  :   Character. Sample type; `"B"` for blank.

  well

  :   Character. Plate well position.

  dilution

  :   Numeric. Dilution factor (1 for blanks).

  mfi

  :   Numeric. Median fluorescence intensity.

  assay_response_variable

  :   Character. Name of the response variable (`"mfi"`).

  assay_independent_variable

  :   Character. Name of the independent variable (`"concentration"`).

- samples:

  A data frame with 120 rows and 13 columns containing patient sample
  measurements:

  curve_id

  :   Integer. Plate identifier (1–6).

  timeperiod

  :   Character. Study visit; one of `"baseline"`, `"month3"`, or
      `"month6"`.

  patientid

  :   Character. Patient identifier (e.g. `"PAT_001"`).

  well

  :   Character. Plate well position.

  stype

  :   Character. Sample type; `"X"` for unknown sample.

  sampleid

  :   Character. Sample identifier.

  agroup

  :   Character. Treatment group; `"GroupA"` or `"GroupB"`.

  dilution

  :   Numeric. Sample dilution factor.

  pctaggbeads

  :   Numeric. Percentage of aggregated beads.

  samplingerrors

  :   Logical. Flags for sampling errors; `NA` if none.

  mfi

  :   Numeric. Median fluorescence intensity.

  assay_response_variable

  :   Character. Name of the response variable (`"mfi"`).

  assay_independent_variable

  :   Character. Name of the independent variable (`"concentration"`).

- curve_id_lookup:

  A data frame with 6 rows and 2 columns mapping plate names to integer
  curve IDs:

  curve_id

  :   Character. Full plate name (e.g. `"alpha_STUDY_PLATE_01"`).

  curve_int

  :   Integer. Integer identifier corresponding to `curve_id` in other
      list elements.

  antigen

  :   Character. Antigen name cooresponding to the curve_id

  study_accession

  :   Character. Study accession that cooresponds to the curve_id

  experiment_accession

  :   Character. Experiment accession that cooresponds to the curve_id

  plate

  :   Plate identifier

- response_var:

  Character. Name of the assay response variable (`"mfi"`).

- indep_var:

  Character. Name of the assay independent variable (`"concentration"`).

## Source

Synthetic data generated for illustrative purposes. Fields of the
samples, standards, and blank controls are structured to align with the
ImmPort data model and can be used in the Interactive Serology Plate
Inspector. Curve parameters were informed by fitting real Luminex data
(MADI P3 GAPS study, ACT and PT antigens, renamed here as alpha and
beta). Standard curves were simulated using the fitted parameters, with
realistic plate-to-plate variability (~8–10\\ inflection point) and
heteroscedastic noise matching Luminex MFI behavior. Unknown sample
concentrations were drawn from across the calibration range. See
`data-raw/bead_assay_example.R` for the data generation script.

## Examples

``` r
data(bead_assay_example)

# Access standard curve data
head(bead_assay_example$standards)
#>   curve_id stype sampleid well    dilution     mfi assay_response_variable
#> 1        1     S   STD_01   A1 1000.000000   109.4                     mfi
#> 2        1     S   STD_02   B1  333.333333   316.9                     mfi
#> 3        1     S   STD_03   C1  100.000000  1133.0                     mfi
#> 4        1     S   STD_04   D1   33.333333  4156.1                     mfi
#> 5        1     S   STD_05   E1   10.000000 12458.1                     mfi
#> 6        1     S   STD_06   F1    3.333333 18933.4                     mfi
#>   assay_independent_variable
#> 1              concentration
#> 2              concentration
#> 3              concentration
#> 4              concentration
#> 5              concentration
#> 6              concentration

# Access patient samples
head(bead_assay_example$samples)
#>   curve_id timeperiod patientid well stype sampleid agroup dilution pctaggbeads
#> 1        1   baseline   PAT_001   A3     X     a001 GroupA     2000        2.49
#> 2        1   baseline   PAT_002   B3     X     a002 GroupB     2000        1.92
#> 3        1     month3   PAT_003   C3     X     a003 GroupA     2000        3.44
#> 4        1   baseline   PAT_004   D3     X     a004 GroupB     2000        3.70
#> 5        1   baseline   PAT_005   E3     X     a005 GroupB     2000        1.15
#> 6        1     month3   PAT_006   F3     X     a006 GroupA     2000        3.40
#>   samplingerrors     mfi assay_response_variable assay_independent_variable
#> 1             NA 18323.4                     mfi              concentration
#> 2             NA 19414.7                     mfi              concentration
#> 3             NA 20098.5                     mfi              concentration
#> 4             NA 19556.0                     mfi              concentration
#> 5             NA 20177.5                     mfi              concentration
#> 6             NA    70.1                     mfi              concentration
```
