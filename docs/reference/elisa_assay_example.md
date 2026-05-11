# ELISA example dataset

A named list containing simulated multi-plate ELISA (enzyme-linked
immunosorbent assay) data, including standard curves, blanks, and
patient samples across multiple time points. The dataset spans six assay
plates for a single analyte (alpha), where plates 1–3 use a
five-parameter logistic (5-PL) standard curve and plates 4–6 use a
Gompertz standard curve, reflecting realistic between-plate variability
in curve shape.

## Usage

``` r
elisa_assay_example
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

  :   Character. Plate well position (e.g. `"A1"`).

  dilution

  :   Numeric. Reciprocal of the known analyte concentration of the
      standard (i.e. `1 / concentration`).

  od

  :   Numeric. Measured optical density.

  assay_response_variable

  :   Character. Name of the response variable (`"od"`).

  assay_independent_variable

  :   Character. Name of the independent variable (`"concentration"`).

- blanks:

  A data frame with 24 rows and 7 columns containing blank well
  measurements (4 blanks per plate). Blanks contain buffer only (no
  antigen, no sample) and represent the lower-asymptote background
  signal of the assay:

  curve_id

  :   Integer. Plate identifier (1–6).

  stype

  :   Character. Sample type; `"B"` for blank.

  well

  :   Character. Plate well position.

  dilution

  :   Numeric. Dilution factor (`1` for blanks).

  od

  :   Numeric. Measured optical density.

  assay_response_variable

  :   Character. Name of the response variable (`"od"`).

  assay_independent_variable

  :   Character. Name of the independent variable (`"concentration"`).

- samples:

  A data frame with 120 rows and 12 columns containing patient sample
  measurements (20 samples per plate). Each sample is run at a fixed
  serum dilution of 1:400:

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

  :   Character. Treatment group; one of `"GroupA"` or `"GroupB"`.

  dilution

  :   Numeric. Fixed serum dilution factor (400).

  samplingerrors

  :   Logical. Flag for sampling errors; `NA` if none detected.

  od

  :   Numeric. Measured optical density.

  assay_response_variable

  :   Character. Name of the response variable (`"od"`).

  assay_independent_variable

  :   Character. Name of the independent variable (`"concentration"`).

- curve_id_lookup:

  A data frame with 6 rows and 4 columns mapping integer curve IDs to
  their plate metadata:

  curve_id

  :   Integer. Integer identifier corresponding to `curve_id` in all
      other list elements.

  antigen

  :   Character. Antigen name corresponding to the curve (`"alpha"` for
      all plates in this dataset).

  study_accession

  :   Character. Study accession corresponding to the curve (e.g.
      `"SDYexample"`).

  experiment_accession

  :   Character. Experiment accession corresponding to the curve (e.g.
      `"EXPexample"`).

  plate

  :   Character. Plate identifier (e.g. `"plate_1"`).

- response_var:

  Character scalar. Name of the assay response variable (`"od"`).

- indep_var:

  Character scalar. Name of the assay independent variable
  (`"concentration"`).

## Source

Synthetic data generated for illustrative purposes. Fields of the
samples, standards, and blank controls are structured to align with the
ImmPort data model and can be used in the Interactive Serology Plate
Inspector. Standard curves were simulated using biologically plausible
5-PL and Gompertz parameters with realistic plate-to-plate variability
(~3–5\\ noise matching typical ELISA plate-reader behaviour
(proportional CV 6\\ from a log-uniform distribution spanning the
calibration range. Samples were simulated at a fixed serum dilution of
1:400. See `data-raw/bead_assay_example.R` for the data generation
script.

## Details

The six plates measure a single antigen (`"alpha"`) but deliberately use
two different standard-curve families to test the robustness of
curve-fitting routines:

- **Plates 1–3** — five-parameter logistic (5-PL), reflecting the
  symmetric sigmoidal response typical of a well-optimised direct ELISA.

- **Plates 4–6** — Gompertz, reflecting the asymmetric response
  sometimes seen when a different substrate lot or incubation time is
  used.

Plate-to-plate variability is introduced by independently perturbing the
upper asymptote (\\d\\), EC50 (\\c\\), and slope (\\b\\) for each plate
using log-normal noise (approximately 3–5\\ follows a heteroscedastic
model (proportional CV of 6\\ additive floor of 0.005 OD), consistent
with typical ELISA plate-reader behaviour. OD values are bounded between
0.001 and 4.0.

## Examples

``` r
data(elisa_assay_example)

# Access standard curve data
head(elisa_assay_example$standards)
#>   curve_id stype sampleid well    dilution     od assay_response_variable
#> 1        1     S   STD_01   A1 1000.000000 0.0955                      od
#> 2        1     S   STD_02   B1  333.333333 0.0865                      od
#> 3        1     S   STD_03   C1  100.000000 0.0672                      od
#> 4        1     S   STD_04   D1   33.333333 0.0937                      od
#> 5        1     S   STD_05   E1   10.000000 0.1885                      od
#> 6        1     S   STD_06   F1    3.333333 0.7303                      od
#>   assay_independent_variable
#> 1              concentration
#> 2              concentration
#> 3              concentration
#> 4              concentration
#> 5              concentration
#> 6              concentration

# Inspect the curve-ID lookup table
elisa_assay_example$curve_id_lookup
#>   curve_id antigen study_accession experiment_accession   plate
#> 1        1   alpha SDYELISAexample           EXPexample plate_1
#> 2        2   alpha SDYELISAexample           EXPexample plate_2
#> 3        3   alpha SDYELISAexample           EXPexample plate_3
#> 4        4   alpha SDYELISAexample           EXPexample plate_4
#> 5        5   alpha SDYELISAexample           EXPexample plate_5
#> 6        6   alpha SDYELISAexample           EXPexample plate_6

# Access patient samples
head(elisa_assay_example$samples)
#>   curve_id timeperiod patientid well stype sampleid agroup dilution
#> 1        1   baseline   PAT_001   A3     X     a001 GroupA      400
#> 2        1   baseline   PAT_002   B3     X     a002 GroupB      400
#> 3        1     month3   PAT_003   C3     X     a003 GroupA      400
#> 4        1   baseline   PAT_004   D3     X     a004 GroupB      400
#> 5        1   baseline   PAT_005   E3     X     a005 GroupB      400
#> 6        1     month3   PAT_006   F3     X     a006 GroupA      400
#>   samplingerrors     od assay_response_variable assay_independent_variable
#> 1             NA 2.7663                      od              concentration
#> 2             NA 3.0137                      od              concentration
#> 3             NA 3.2028                      od              concentration
#> 4             NA 2.9290                      od              concentration
#> 5             NA 2.3813                      od              concentration
#> 6             NA 0.0598                      od              concentration

# Response and independent variable names
elisa_assay_example$response_var   # "od"
#> [1] "od"
elisa_assay_example$indep_var      # "concentration"
#> [1] "concentration"
```
