# Synthetic Multi-Plate Bead Assay Example Dataset

A synthetic dataset that mirrors the structure expected by
[`fit_calibration_freq_multiplate()`](https://immunoplex.github.io/curveRfreq/reference/fit_calibration_freq_multiplate.md).
Contains two antigens (alpha and beta) across three replicate plates
each (six curve_ids total).

A synthetic dataset that mirrors the structure expected by
[`fit_calibration_freq_multiplate()`](https://immunoplex.github.io/curveRfreq/reference/fit_calibration_freq_multiplate.md).
Contains two antigens (alpha and beta) across three replicate plates
each (six curve_ids total).

## Format

A list with six elements:

- standards:

  Data frame, 60 rows. One row per standard well.

- blanks:

  Data frame, 24 rows. Four blank wells per plate.

- samples:

  Data frame, 120 rows. Patient samples at dilution 1:2000.

- curve_id_lookup:

  Data frame, 6 rows. Maps curve_id to antigen/plate.

- response_var:

  Character. Name of the response column ("mfi").

- indep_var:

  Character. Name of the independent variable ("concentration").

A list with six elements:

- standards:

  Data frame, 60 rows. One row per standard well.

- blanks:

  Data frame, 24 rows. Four blank wells per plate.

- samples:

  Data frame, 120 rows. Patient samples at dilution 1:2000.

- curve_id_lookup:

  Data frame, 6 rows. Maps curve_id to antigen/plate.

- response_var:

  Character. Name of the response column ("mfi").

- indep_var:

  Character. Name of the independent variable ("concentration").

## Source

Synthetic data generated for package documentation.

Synthetic data generated for package documentation.

## Details

The alpha curves were simulated with a Gompertz model; the beta curves
with a 5PL model, reflecting realistic between-antigen variation in
curve shape.

The alpha curves were simulated with a Gompertz model; the beta curves
with a 5PL model, reflecting realistic between-antigen variation in
curve shape.
