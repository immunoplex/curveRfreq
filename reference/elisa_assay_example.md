# Synthetic Multi-Plate ELISA Example Dataset

A synthetic dataset that mirrors the structure expected by
[`fit_calibration_freq_multiplate()`](https://immunoplex.github.io/curveRfreq/reference/fit_calibration_freq_multiplate.md)
for an ELISA (optical density) assay. Contains two antigens across three
replicate plates each (six curve_ids total).

A synthetic dataset that mirrors the structure expected by
[`fit_calibration_freq_multiplate()`](https://immunoplex.github.io/curveRfreq/reference/fit_calibration_freq_multiplate.md)
for an ELISA (optical density) assay. Contains two antigens across three
replicate plates each (six curve_ids total).

## Format

A list with elements analogous to
[bead_assay_example](https://immunoplex.github.io/curveRfreq/reference/bead_assay_example.md):

- standards:

  Data frame. One row per standard well with an `od` response column.

- blanks:

  Data frame. Blank wells per plate.

- samples:

  Data frame. Patient samples.

- curve_id_lookup:

  Data frame. Maps curve_id to antigen/plate.

- response_var:

  Character. Name of the response column ("od").

- indep_var:

  Character. Name of the independent variable ("concentration").

A list with elements analogous to
[bead_assay_example](https://immunoplex.github.io/curveRfreq/reference/bead_assay_example.md):

- standards:

  Data frame. One row per standard well with an `od` response column.

- blanks:

  Data frame. Blank wells per plate.

- samples:

  Data frame. Patient samples.

- curve_id_lookup:

  Data frame. Maps curve_id to antigen/plate.

- response_var:

  Character. Name of the response column ("od").

- indep_var:

  Character. Name of the independent variable ("concentration").

## Source

Synthetic data generated for package documentation.

Synthetic data generated for package documentation.

## Details

Simulated OD-based immunoassay data for testing and vignette examples.

Simulated OD-based immunoassay data for testing and vignette examples.
