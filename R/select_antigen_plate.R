#' Select and Filter Data for a Specific Antigen-Plate Combination
#'
#' Filters standards, blanks, samples, MCMC samples, and MCMC predictions
#' from a loaded data object for a specific antigen, source, and plate
#' combination. Also resolves antigen settings, the fixed lower asymptote,
#' and blank standard error for downstream curve fitting.
#'
#' @param loaded_data A named list as returned by \code{pull_data()}, containing
#'   elements \code{standards}, \code{blanks}, \code{samples},
#'   \code{mcmc_samples}, and \code{mcmc_pred}.
#' @param study_accession Character. Study accession identifier.
#' @param experiment_accession Character. Experiment accession identifier.
#' @param source Character. Source label used to filter standards
#'   (matched against \code{source_nom} if present, otherwise \code{source}).
#' @param antigen Character. Antigen identifier.
#' @param plate Character. Plate label in \code{plate_nom} format
#'   (e.g. \code{"plate_13-1"}).
#' @param wavelength Character. Wavelength filter for ELISA multi-channel data.
#'   Defaults to \code{WL_NONE}, which skips wavelength filtering.
#' @param antigen_constraints Data frame of antigen-level constraints,
#'   pre-filtered to the relevant antigen. See \code{fetch_antigen_parameters()}.
#' @param verbose Logical. If \code{TRUE}, prints diagnostic information
#'   during filtering and constraint resolution. Default \code{TRUE}.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{plate_standard}{Data frame of filtered standard curve wells.}
#'     \item{plate_blanks}{Data frame of filtered blank (buffer) wells.}
#'     \item{plate_samples}{Data frame of filtered sample wells.}
#'     \item{plate_mcmc_samples}{Data frame of filtered MCMC sample
#'       concentration estimates, or an empty data frame if unavailable.}
#'     \item{plate_mcmc_pred}{Data frame of filtered MCMC prediction grid,
#'       sorted by \code{x}, or an empty data frame if unavailable.}
#'     \item{antigen_settings}{Named list of lower asymptote constraint
#'       parameters from \code{\link{obtain_lower_constraint}}.}
#'     \item{fixed_a_result}{Numeric scalar if the lower asymptote is fixed,
#'       or \code{NULL} if it is a free parameter.}
#'     \item{std_error_blank}{Numeric. Standard error of blank responses.}
#'   }
#'   Returns \code{NULL} if no standard curve data is found for the specified
#'   combination of \code{source}, \code{antigen}, and \code{plate}.
#'
#' @details
#' Filtering uses \code{source_nom} over \code{source} when available, which
#' is important for ELISA data where the source label encodes wavelength
#' information (e.g. \code{"Sando|450_nm"}).
#'
#' For ELISA data with multiple wavelength channels, pass the relevant
#' \code{wavelength} value to filter all data frames consistently.
#'
#' The \code{plate} argument should match the \code{plate_nom} column, which
#' is constructed as \code{paste(plate, nominal_sample_dilution, sep = "-")}.
#' The nominal dilution suffix is stripped internally when calling
#' \code{\link{obtain_lower_constraint}}.
#'
#' @seealso
#'   \code{\link{pull_data}},
#'   \code{\link{obtain_lower_constraint}},
#'   \code{\link{resolve_fixed_lower_asymptote}},
#'   \code{\link{validate_fixed_lower_asymptote}},
#'   \code{\link{get_blank_se}}
#' @export
select_antigen_plate <- function(loaded_data,
                                 study_accession,
                                 experiment_accession,
                                 source,
                                 antigen,
                                 plate,
                                 wavelength          = WL_NONE,
                                 antigen_constraints,
                                 verbose             = TRUE) {

  if (verbose) {
    message("[select_antigen_plate] plate: ", plate)
    message("[select_antigen_plate] wavelength: ", wavelength)
    message("[select_antigen_plate] antigens in standards: ",
            paste(unique(loaded_data$standards$antigen), collapse = ", "))
    message("[select_antigen_plate] source_nom values: ",
            paste(unique(loaded_data$standards$source_nom), collapse = ", "))
    message("[select_antigen_plate] plate_nom values: ",
            paste(unique(loaded_data$standards$plate_nom), collapse = ", "))
  }

  # ── Filter standards ───────────────────────────────────────────────
  if ("source_nom" %in% names(loaded_data$standards)) {
    plate_standard <- loaded_data$standards[
      loaded_data$standards$source_nom == source &
        loaded_data$standards$antigen    == antigen &
        loaded_data$standards$plate_nom  == plate, ]
  } else {
    plate_standard <- loaded_data$standards[
      loaded_data$standards$source    == source &
        loaded_data$standards$antigen   == antigen &
        loaded_data$standards$plate_nom == plate, ]
  }

  # ── Filter by wavelength for standards ────────────────────────────
  if ("wavelength" %in% names(plate_standard) &&
      !is.null(wavelength) && wavelength != WL_NONE) {
    plate_standard$wavelength <- normalize_wavelength(plate_standard$wavelength)
    wl_filter <- plate_standard$wavelength == normalize_wavelength(wavelength)
    if (any(wl_filter)) {
      plate_standard <- plate_standard[wl_filter, , drop = FALSE]
    } else {
      message(sprintf(
        "[select_antigen_plate] WARNING: wavelength '%s' matched 0 rows; keeping all %d rows. Wavelengths in data: %s",
        wavelength, nrow(plate_standard),
        paste(unique(plate_standard$wavelength), collapse = ", ")
      ))
    }
  }

  # ── Guard against empty plate_standard data ───────────────────────
  if (is.null(plate_standard) || nrow(plate_standard) == 0) {
    warning(paste("No standard curve data found for:",
                  "source =", source,
                  ", antigen =", antigen,
                  ", plate =", plate))
    return(NULL)
  }

  # ── Filter blanks ─────────────────────────────────────────────────
  plate_blanks <- loaded_data$blanks[
    loaded_data$blanks$antigen   == antigen &
      loaded_data$blanks$plate_nom == plate, ]

  # ── Filter samples ────────────────────────────────────────────────
  plate_samples <- loaded_data$samples[
    loaded_data$samples$antigen   == antigen &
      loaded_data$samples$plate_nom == plate, ]

  # ── Filter mcmc_samples ───────────────────────────────────────────
  plate_mcmc_samples <- if (!is.null(loaded_data$mcmc_samples) &&
                            nrow(loaded_data$mcmc_samples) > 0) {
    mcmc_df     <- loaded_data$mcmc_samples
    filter_mask <- mcmc_df$antigen == antigen & mcmc_df$plate_nom == plate
    if ("source_nom" %in% names(mcmc_df)) {
      filter_mask <- filter_mask & mcmc_df$source_nom == source
    }
    mcmc_df[filter_mask, , drop = FALSE]
  } else {
    data.frame()
  }

  # ── Filter mcmc_pred ──────────────────────────────────────────────
  plate_mcmc_pred <- if (!is.null(loaded_data$mcmc_pred) &&
                         nrow(loaded_data$mcmc_pred) > 0) {
    pred_df     <- loaded_data$mcmc_pred
    filter_mask <- pred_df$antigen == antigen & pred_df$plate_nom == plate
    if ("source_nom" %in% names(pred_df)) {
      filter_mask <- filter_mask & pred_df$source_nom == source
    }
    pred_df[filter_mask, , drop = FALSE]
  } else {
    data.frame()
  }

  # ── Filter all data frames by wavelength ──────────────────────────
  if (!is.null(wavelength) && wavelength != WL_NONE) {
    if ("wavelength" %in% names(plate_blanks) && nrow(plate_blanks) > 0) {
      wl_b <- plate_blanks$wavelength == normalize_wavelength(wavelength)
      if (any(wl_b)) plate_blanks <- plate_blanks[wl_b, , drop = FALSE]
    }
    if ("wavelength" %in% names(plate_samples) && nrow(plate_samples) > 0) {
      wl_s <- plate_samples$wavelength == normalize_wavelength(wavelength)
      if (any(wl_s)) plate_samples <- plate_samples[wl_s, , drop = FALSE]
    }
    if ("wavelength" %in% names(plate_mcmc_samples) && nrow(plate_mcmc_samples) > 0) {
      wl_m <- plate_mcmc_samples$wavelength == normalize_wavelength(wavelength)
      if (any(wl_m)) plate_mcmc_samples <- plate_mcmc_samples[wl_m, , drop = FALSE]
    }
    if ("wavelength" %in% names(plate_mcmc_pred) && nrow(plate_mcmc_pred) > 0) {
      wl_p <- plate_mcmc_pred$wavelength == normalize_wavelength(wavelength)
      if (any(wl_p)) plate_mcmc_pred <- plate_mcmc_pred[wl_p, , drop = FALSE]
    }
  }

  # ── Strip nominal dilution suffix for obtain_lower_constraint ─────
  plate_c <- sub("-.*$", "", plate)

  # ── Resolve response column ───────────────────────────────────────
  response_col <- resolve_response_col(plate_standard)

  # ── Antigen settings ──────────────────────────────────────────────
  antigen_settings <- obtain_lower_constraint(
    dat                  = plate_standard,
    antigen              = antigen,
    study_accession      = study_accession,
    experiment_accession = experiment_accession,
    plate                = plate_c,
    plateid              = unique(plate_standard$plateid),
    plate_blanks         = plate_blanks,
    antigen_constraints  = antigen_constraints,
    response_col         = response_col
  )

  # ── Fixed lower asymptote ─────────────────────────────────────────
  fixed_a_result <- resolve_fixed_lower_asymptote(antigen_settings)
  fixed_a_result <- validate_fixed_lower_asymptote(
    fixed_a_result_raw = fixed_a_result,
    verbose            = verbose
  )

  # ── Blank standard error ──────────────────────────────────────────
  std_error_blank <- get_blank_se(antigen_settings = antigen_settings)

  # ── Sort mcmc_pred by x for smooth line drawing ───────────────────
  if (nrow(plate_mcmc_pred) > 0 && "x" %in% names(plate_mcmc_pred)) {
    plate_mcmc_pred <- plate_mcmc_pred[order(plate_mcmc_pred$x), , drop = FALSE]
  }

  # ── Return ────────────────────────────────────────────────────────
  return(list(
    plate_standard     = plate_standard,
    plate_blanks       = plate_blanks,
    plate_samples      = plate_samples,
    plate_mcmc_samples = plate_mcmc_samples,
    plate_mcmc_pred    = plate_mcmc_pred,
    antigen_settings   = antigen_settings,
    fixed_a_result     = fixed_a_result,
    std_error_blank    = std_error_blank
  ))
}
