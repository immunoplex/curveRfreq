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
#' @name bead_assay_example
#' @docType data
NULL

#' Synthetic Multi-Plate ELISA Example Dataset
#'
#' A synthetic dataset that mirrors the structure expected by
#' [fit_calibration_freq_multiplate()] for an ELISA (optical density)
#' assay. Contains two antigens across three replicate plates each
#' (six curve_ids total).
#'
#' @format A list with elements analogous to [bead_assay_example]:
#' \describe{
#'   \item{standards}{Data frame. One row per standard well with an `od`
#'     response column.}
#'   \item{blanks}{Data frame. Blank wells per plate.}
#'   \item{samples}{Data frame. Patient samples.}
#'   \item{curve_id_lookup}{Data frame. Maps curve_id to antigen/plate.}
#'   \item{response_var}{Character. Name of the response column ("od").}
#'   \item{indep_var}{Character. Name of the independent variable
#'     ("concentration").}
#' }
#'
#' @details
#' Simulated OD-based immunoassay data for testing and vignette examples.
#'
#' @source Synthetic data generated for package documentation.
#' @name elisa_assay_example
#' @docType data
NULL
