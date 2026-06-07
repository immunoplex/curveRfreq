# R/data.R
#' Synthetic Multi-Plate Bead Assay Example Dataset
#'
#' A synthetic dataset that mirrors the structure expected by
#' [fit_calibration_freq_multiplate()]. Contains two antigens (alpha and beta)
#' across three replicate plates each (six curve_ids total).
#'
#' @format A list with six elements:
#' \describe{
#'   \item{standards}{Data frame, 60 rows. One row per standard well.}
#'   \item{blanks}{Data frame, 24 rows. Four blank wells per plate.}
#'   \item{samples}{Data frame, 120 rows. Patient samples at dilution 1:2000.}
#'   \item{curve_id_lookup}{Data frame, 6 rows. Maps curve_id to antigen/plate.}
#'   \item{response_var}{Character. Name of the response column ("mfi").}
#'   \item{indep_var}{Character. Name of the independent variable ("concentration").}
#' }
#'
#' @details
#' The alpha curves were simulated with a Gompertz model; the beta curves
#' with a 5PL model, reflecting realistic between-antigen variation in
#' curve shape.
#'
#' @source Synthetic data generated for package documentation.
"bead_assay_example"
