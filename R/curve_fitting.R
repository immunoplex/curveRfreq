# =============================================================================
# curve_fitting.R
#
# Standard curve fitting pipeline for immunoassay data (bead array & ELISA).
# Covers data preparation, blank handling, prozone correction, nonlinear
# least-squares model fitting (4PL / 5PL / Gompertz), model selection by
# AIC, and error propagation via the delta method.
#
# Dependencies: minpack.lm, nlstools, tibble, dplyr
# =============================================================================


# -----------------------------------------------------------------------------
# Section 1 -- Constraint helpers
# -----------------------------------------------------------------------------

#' Obtain lower asymptote constraints for a given antigen
#'
#' Returns a named list of lower-asymptote constraint parameters for one
#' antigen / plate combination.  The method used is determined by
#' \code{antigen_constraints$l_asy_constraint_method}:
#' \describe{
#'   \item{\code{"default"}}{Min = 0, max = observed data maximum.}
#'   \item{\code{"user_defined"}}{Uses the values stored in
#'         \code{antigen_constraints}.}
#'   \item{\code{"range_of_blanks"}}{Min/max set to the range of blank
#'         plate responses.}
#'   \item{\code{"geometric_mean_of_blanks"}}{Both min and max set to the
#'         geometric mean of blank responses (i.e. the parameter is effectively
#'         fixed).}
#' }
#'
#' @param dat               Data frame of standards assay measurements for this antigen
#'                          and plate (standards + samples).
#' @param antigen           Character. Antigen identifier.
#' @param study_accession   Character. Study accession identifier.
#' @param experiment_accession Character. Experiment accession identifier.
#' @param plate             Character. Plate label.
#' @param plate_blanks      Data frame of blank (buffer) wells for this plate.
#' @param antigen_constraints Data frame with exactly one row (or more --- the
#'                          first row is used) containing constraint columns:
#'                          \code{l_asy_constraint_method},
#'                          \code{l_asy_min_constraint},
#'                          \code{l_asy_max_constraint},
#'                          \code{standard_curve_concentration},
#'                          \code{pcov_threshold}.
#' @param response_col      Character or \code{NULL}.  Name of the response
#'                          column (\code{"mfi"} for bead arrays,
#'                          \code{"absorbance"} for ELISA).  If \code{NULL}
#'                          it is resolved automatically via
#'                          \code{resolve_response_col()}.
#'
#' @return A named list with elements:
#'   \item{study_accession}{Passed through.}
#'   \item{experiment_accession}{Passed through.}
#'   \item{plate}{Passed through.}
#'   \item{antigen}{Passed through.}
#'   \item{l_asy_min_constraint}{Numeric lower bound for the \emph{a}
#'         parameter.}
#'   \item{l_asy_max_constraint}{Numeric upper bound for the \emph{a}
#'         parameter.}
#'   \item{l_asy_constraint_method}{Character. The method used.}
#'   \item{std_error_blank}{Numeric. Standard error of blank responses.}
#'   \item{standard_curve_concentration}{Numeric. Concentration of the
#'         undiluted standard.}
#'   \item{pcov_threshold}{Numeric. Percent-CV acceptance threshold.}
#'   Returns \code{NULL} if the constraint method is unrecognized.
#'
#' @export
obtain_lower_constraint <- function(dat, antigen, study_accession, experiment_accession,
                                    plate, plate_blanks, antigen_constraints,
                                    response_col = NULL) {

  if (is.null(response_col)) response_col <- resolve_response_col(dat)

  # If multiple rows were passed, use only the first
  if (is.data.frame(antigen_constraints) && nrow(antigen_constraints) > 1) {
    # warning(paste("Multiple constraint rows found for antigen:", antigen,
    #               "- using first row. Consider deduplicating antigen_constraints."))
    antigen_constraints <- antigen_constraints[1, , drop = FALSE]
  }

  # Helper: safely extract first non-NA scalar value
  safe_extract <- function(x, default = NA) {
    if (is.null(x) || length(x) == 0) return(default)
    x <- x[!is.na(x)]
    if (length(x) == 0) return(default)
    return(x[1])
  }

  constraint_method <- safe_extract(trimws(antigen_constraints$l_asy_constraint_method), "default")
  l_asy_min   <- safe_extract(antigen_constraints$l_asy_min_constraint, 0)
  l_asy_max   <- safe_extract(antigen_constraints$l_asy_max_constraint, NA)
  std_curve_conc <- safe_extract(antigen_constraints$standard_curve_concentration, 10000)
  pcov_thresh <- safe_extract(antigen_constraints$pcov_threshold, 20)

  # Blank standard error (on the response scale)
  blank_response_col <- resolve_response_col(plate_blanks, default = response_col)
  if (nrow(plate_blanks) > 1) {
    se_blank_response <- sd(plate_blanks[[blank_response_col]], na.rm = TRUE) /
      sqrt(sum(!is.na(plate_blanks[[blank_response_col]])))
  } else {
    se_blank_response <- 0
  }

  if (constraint_method == "user_defined") {
    l_asy_constraints <- list(
      study_accession        = study_accession,
      experiment_accession   = experiment_accession,
      plate                  = plate,
      antigen                = antigen,
      l_asy_min_constraint   = l_asy_min,
      l_asy_max_constraint   = l_asy_max,
      l_asy_constraint_method = constraint_method,
      std_error_blank        = se_blank_response,
      standard_curve_concentration = std_curve_conc,
      pcov_threshold         = pcov_thresh
    )
  } else if (constraint_method == "default") {
    l_asy_constraints <- list(
      study_accession        = study_accession,
      experiment_accession   = experiment_accession,
      plate                  = plate,
      antigen                = antigen,
      l_asy_min_constraint   = 0,
      l_asy_max_constraint   = max(dat[[response_col]], na.rm = TRUE),
      l_asy_constraint_method = constraint_method,
      std_error_blank        = se_blank_response,
      standard_curve_concentration = std_curve_conc,
      pcov_threshold         = pcov_thresh
    )
  } else if (constraint_method == "range_of_blanks") {
    l_asy_constraints <- list(
      study_accession        = study_accession,
      experiment_accession   = experiment_accession,
      plate                  = plate,
      antigen                = antigen,
      l_asy_min_constraint   = min(plate_blanks[[blank_response_col]], na.rm = TRUE),
      l_asy_max_constraint   = max(plate_blanks[[blank_response_col]], na.rm = TRUE),
      l_asy_constraint_method = constraint_method,
      std_error_blank        = se_blank_response,
      standard_curve_concentration = std_curve_conc,
      pcov_threshold         = pcov_thresh
    )
  } else if (constraint_method == "geometric_mean_of_blanks") {
    geometric_mean <- exp(mean(log(plate_blanks[[blank_response_col]]), na.rm = TRUE))
    l_asy_constraints <- list(
      study_accession        = study_accession,
      experiment_accession   = experiment_accession,
      plate                  = plate,
      antigen                = antigen,
      l_asy_min_constraint   = geometric_mean,
      l_asy_max_constraint   = geometric_mean,
      l_asy_constraint_method = constraint_method,
      std_error_blank        = se_blank_response,
      standard_curve_concentration = std_curve_conc,
      pcov_threshold         = pcov_thresh
    )
  } else {
    return(NULL)
  }

  return(l_asy_constraints)
}


#' Extract Undiluted Standard Curve Concentration from Constraint List
#'
#' Convenience accessor that retrieves the nominal concentration of the
#' undiluted standard from the constraint list produced by
#' \code{\link{obtain_lower_constraint}}.
#'
#' @param l_asy_constraints Named list returned by
#'   \code{\link{obtain_lower_constraint}}.
#'
#' @return Numeric scalar: the undiluted standard curve concentration
#'   (e.g. \code{10000}).
#'
#' @export
get_study_exp_antigen_plate_param <- function(l_asy_constraints) {
  undiluted_sc_concentration <- l_asy_constraints$standard_curve_concentration
  return(undiluted_sc_concentration)
}


# -----------------------------------------------------------------------------
# Section 2 -- Concentration and response transforms
# -----------------------------------------------------------------------------

#' Compute Concentration Column from Dilution and Undiluted Standard
#'
#' Converts the \code{dilution} column in \code{data} to an absolute
#' concentration by multiplying the reciprocal of each dilution factor by the
#' undiluted standard concentration.  Optionally applies a log10 transform.
#'
#' @param data                     Data frame containing standards with a dilution column
#' @param undiluted_sc_concentration Numeric. Concentration of the undiluted
#'                                 standard (e.g. \code{10000}).
#' @param independent_variable     Character. Name of the column to create or
#'                                 overwrite with concentration values.
#' @param is_log_concentration     Logical. If \code{TRUE} (default), the
#'                                 computed concentration is log10-transformed.
#'
#' @return \code{data} with the \code{independent_variable} column populated.
#'
#' @export
compute_concentration <- function(data,
                                  undiluted_sc_concentration,
                                  independent_variable,
                                  is_log_concentration = TRUE) {
  independent_variable <- unique(independent_variable)
  data[[independent_variable]] <- (1 / data$dilution) * undiluted_sc_concentration

  if (is_log_concentration) {
    data[[independent_variable]] <- log10(data[[independent_variable]])
  }

  return(data)
}


#' Log10-Transform the Assay Response Variable
#'
#' Applies \code{log10()} to the response column when \code{is_log_response}
#' is \code{TRUE}.  Non-positive values should be floored \emph{before}
#' calling this function (see \code{\link{preprocess_robust_curves}}).
#'
#' @param data              Data frame.
#' @param response_variable Character. Name of the response column.
#' @param is_log_response   Logical. If \code{TRUE} (default), apply
#'                          \code{log10}.
#'
#' @return \code{data} with the response column optionally log-transformed.
#'
#' @export
compute_log_response <- function(data, response_variable, is_log_response = TRUE) {
  if (is_log_response) {
    data[[response_variable]] <- log10(data[[response_variable]])
  }
  return(data)
}


# -----------------------------------------------------------------------------
# Section 3 -- Prozone (hook effect) correction
# -----------------------------------------------------------------------------

#' Correct the Prozone (Hook) Effect in Standard Curve Data
#'
#' The prozone effect occurs when analyte concentration is so high that it
#' saturates antibody binding sites, causing the measured signal to decrease
#' at very high concentrations.  This function dampens the apparent decrease
#' beyond the signal peak by compressing the post-peak delta toward the peak
#' value.
#'
#' The method follows the computational approach described in:
#' \itemize{
#'   \item ACS Meas. Sci. Au 2024, 4, 4, 452--458.
#'   \item Sensors and Actuators B: Chemical, 324, 128756 (2020).
#'   \item Sensors and Actuators B: Chemical, 304, 127408 (2020).
#' }
#'
#' @param stdframe           Data frame of standard curve data for a single
#'                           dilution series.  Must contain
#'                           \code{response_variable} and
#'                           \code{independent_variable} columns with no NA
#'                           values.
#' @param prop_diff          Numeric. Dampening factor applied to the
#'                           post-peak signal delta (e.g. \code{0.1}).
#' @param dil_scale          Numeric. Dilution scale factor used in the
#'                           dampening formula (e.g. \code{2}).
#' @param response_variable  Character. Name of the response (y) column
#'                           (default \code{"mfi"}).
#' @param independent_variable Character. Name of the concentration (x)
#'                           column (default \code{"concentration"}).
#' @param verbose            Logical. If \code{TRUE} (default), prints
#'                           diagnostics to the console.
#'
#' @return \code{stdframe} with post-peak response values adjusted.
#'
#' @export
correct_prozone <- function(stdframe = NULL, prop_diff = NULL, dil_scale = NULL,
                            response_variable  = "mfi",
                            independent_variable = "concentration",
                            verbose = TRUE) {

  response_variable <- unique(response_variable)
  stdframe <- stdframe[!is.na(stdframe[[response_variable]]) &
                         !is.na(stdframe[[independent_variable]]), ]

  max_response       <- max(stdframe[[response_variable]], na.rm = TRUE)
  logc_at_max_response <- max(
    stdframe[stdframe[[response_variable]] == max_response, ][[independent_variable]]
  )

  if (verbose) cat("Peak MFI =", max_response, "at concentration =", logc_at_max_response, "\n")

  post_peak <- stdframe[[independent_variable]] > logc_at_max_response
  if (verbose) cat("Number of points beyond the peak:", sum(post_peak), "\n")

  stdframe[stdframe[[independent_variable]] > logc_at_max_response, ][[response_variable]] <-
    max_response + (
      (max_response - stdframe[stdframe[[independent_variable]] > logc_at_max_response, ][[response_variable]]) *
        prop_diff /
        ((stdframe[stdframe[[independent_variable]] > logc_at_max_response, ][[independent_variable]] -
            logc_at_max_response) * dil_scale)
    )

  return(stdframe)
}


# -----------------------------------------------------------------------------
# Section 4 -- Blank handling
# -----------------------------------------------------------------------------

#' Extract the Standard Error of the Blank from Antigen Settings
#'
#' @param antigen_settings Named list (output of
#'   \code{\link{obtain_lower_constraint}}).
#'
#' @return Numeric scalar: \code{std_error_blank}.
#'
#' @export
get_blank_se <- function(antigen_settings) {
  return(antigen_settings$std_error_blank)
}

#' Geometric Mean (Helper)
#'
#' Computes \code{exp(mean(log(x)))}, the geometric mean of a numeric vector.
#'
#' @param x     Numeric vector (all positive).
#' @param na.rm Logical. Remove \code{NA} values (default \code{TRUE}).
#'
#' @return Numeric scalar: the geometric mean of \code{x}.
#'
#' @keywords internal
geom_mean <- function(x, na.rm = TRUE) {
  exp(mean(log(x), na.rm = na.rm))
}


#' Include Blank Controls as an Extra Point on the Standard Curve
#'
#' Appends a synthetic data row whose response is the geometric mean of the
#' blank controls and whose concentration is half of the minimum standard
#' concentration.  This anchors the lower end of the fitted curve to the
#' background signal level.
#'
#' @param blank_data         Data frame of blank well measurements.
#' @param data               Data frame of standard curve measurements.
#' @param response_variable  Character. Name of the response column.
#' @param independent_variable Character. Name of the concentration column
#'                           (default \code{"concentration"}).
#'
#' @return \code{data} with one additional row representing the blank mean.
#'
#' @export
include_blanks_conc <- function(blank_data, data, response_variable,
                                independent_variable = "concentration") {

  data <- data[, !(names(data) %in% c("dilution_factor", "log_dilution"))]

  response_blank  <- geom_mean(blank_data[[response_variable]])
  min_concentration <- min(data[[independent_variable]], na.rm = TRUE)
  conc_blank      <- min_concentration - log10(2)

  new_point <- tibble::tibble(
    project_id              = unique(data$project_id),
    study_accession         = unique(data$study_accession),
    experiment_accession    = unique(data$experiment_accession),
    feature                 = unique(data$feature),
    source                  = unique(data$source),
    plateid                 = unique(data$plateid),
    plate                   = unique(data$plate),
    stype                   = "B",
    nominal_sample_dilution = unique(data$nominal_sample_dilution),
    sampleid                = "blank_mean",
    well                    = "geometric_mean_blank",
    dilution                = NA_real_,
    antigen                 = unique(data$antigen),
    !!response_variable    := response_blank,
    assay_response_variable = unique(data$assay_response_variable),
    assay_independent_variable = unique(data$assay_independent_variable),
    concentration           = conc_blank
  )

  # Carry optional columns forward
  if ("source_nom" %in% names(data))  new_point$source_nom  <- unique(data$source_nom)[1]
  if ("wavelength" %in% names(data))  new_point$wavelength  <- unique(data$wavelength)[1]

  if ("plate_nom" %in% names(data)) {
    new_point$plate_nom <- unique(data$plate_nom)[1]
    nm <- names(new_point)
    i  <- match("assay_independent_variable", nm)
    new_point <- new_point[, c(nm[1:i], "plate_nom", nm[(i + 1):(length(nm) - 1)])]
  }

  data_with_blank <- rbind(data, new_point)
  return(data_with_blank)
}


#' Apply a Blank Operation to Standard Curve Data
#'
#' Performs one of five blank-handling strategies on the standard curve data:
#' \describe{
#'   \item{\code{"ignored"}}{No blank adjustment (default).}
#'   \item{\code{"included"}}{Append blank geometric mean as an extra curve
#'         point via \code{\link{include_blanks_conc}}.}
#'   \item{\code{"subtracted"}}{Subtract the geometric mean of blanks from
#'         all responses.}
#'   \item{\code{"subtracted_3x"}}{Subtract three times the geometric mean.}
#'   \item{\code{"subtracted_10x"}}{Subtract ten times the geometric mean.}
#' }
#' After subtraction, values that become zero or negative are floored at 0
#' (linear scale) or 1 (log scale) to prevent downstream errors.
#'
#' @param blank_data         Data frame of blank measurements.
#' @param data               Data frame of standard curve data.
#' @param response_variable  Character. Name of the response column.
#' @param independent_variable Character. Name of the concentration column.
#' @param is_log_response    Logical. Whether the response has been
#'                           log10-transformed.
#' @param blank_option       Character. One of \code{"ignored"},
#'                           \code{"included"}, \code{"subtracted"},
#'                           \code{"subtracted_3x"}, \code{"subtracted_10x"}.
#' @param verbose            Logical. Emit a message describing which
#'                           operation was applied (default \code{TRUE}).
#'
#' @return \code{data} with the blank operation applied.
#'
#' @export
perform_blank_operation <- function(blank_data, data, response_variable,
                                    independent_variable, is_log_response,
                                    blank_option = "ignored", verbose = TRUE) {

  if (verbose) message("Blank Option Used: ", blank_option)

  valid_options <- c("ignored", "included", "subtracted",
                     "subtracted_3x", "subtracted_10x")

  if (!(blank_option %in% valid_options)) {
    message("Invalid value for blank_option. Must be one of: ",
            paste(valid_options, collapse = ", "), ".")
    return(data)
  }

  if (blank_option != "ignored" &&
      (is.null(blank_data) || nrow(blank_data) == 0)) {
    message("Blank data must be supplied when blank_option is not 'ignored'.")
    return(data)
  }

  if (blank_option == "included") {
    data <- include_blanks_conc(blank_data        = blank_data,
                                data              = data,
                                response_variable = response_variable,
                                independent_variable = independent_variable)
    if (verbose) message("Geometric mean of blanks included as an extra point.")
  }

  if (blank_option %in% c("subtracted", "subtracted_3x", "subtracted_10x")) {
    factor <- switch(blank_option,
                     "subtracted"    = 1,
                     "subtracted_3x" = 3,
                     "subtracted_10x" = 10)

    blank_mean <- geom_mean(blank_data[[response_variable]])
    data[[response_variable]] <- data[[response_variable]] - factor * blank_mean

    if (is_log_response) {
      data[data[[response_variable]] < 0, response_variable] <- 1
    } else {
      data[data[[response_variable]] < 0, response_variable] <- 0
    }

    if (verbose) message("Performed blank subtraction (x", factor, ") in linear space.")
  }

  return(data)
}




# -----------------------------------------------------------------------------
# Section 5 -- Lower asymptote constraint resolution
# -----------------------------------------------------------------------------

#' Determine Whether the Lower Asymptote Should Be Fixed
#'
#' Tests whether the min and max constraints for \emph{a} are identical.  If
#' they are, the parameter is effectively fixed and its value is returned;
#' otherwise \code{NULL} is returned to indicate a free parameter.
#'
#' @param l_asy_constraints Named list (output of
#'   \code{\link{obtain_lower_constraint}}).
#'
#' @return Numeric scalar (the fixed value) if min == max; \code{NULL} if the
#'   parameter should be estimated freely.
#'
#' @export
resolve_fixed_lower_asymptote <- function(l_asy_constraints) {
  if (l_asy_constraints$l_asy_min_constraint ==
      l_asy_constraints$l_asy_max_constraint) {
    return(l_asy_constraints$l_asy_min_constraint)
  } else {
    return(NULL)
  }
}


#' Validate a Raw Fixed Lower Asymptote Before Log Transformation
#'
#' Checks that \code{fixed_a_result_raw} is a positive, finite scalar before
#' it is passed to \code{log10()} downstream.  Returns the original value if
#' valid, or \code{NULL} to indicate the parameter should be treated as free.
#'Called by: predict_and_propagate_error(), select_model_formulas() callers
#'
#' @param fixed_a_result_raw Numeric scalar (or \code{NULL}).
#' @param verbose            Logical. Emit diagnostic messages (default
#'                           \code{TRUE}).
#'
#' @return \code{fixed_a_result_raw} if valid; \code{NULL} otherwise.
#'
#' @export
validate_fixed_lower_asymptote <- function(fixed_a_result_raw, verbose = TRUE) {

  if (is.null(fixed_a_result_raw)) return(NULL)

  if (!is.numeric(fixed_a_result_raw) || length(fixed_a_result_raw) != 1) {
    if (verbose) message(
      "[validate_fixed_lower_asymptote] Not a scalar numeric --- treating as NULL (free)."
    )
    return(NULL)
  }

  if (!is.finite(fixed_a_result_raw)) {
    if (verbose) message(sprintf(
      "[validate_fixed_lower_asymptote] = %s is not finite --- treating as NULL.",
      as.character(fixed_a_result_raw)
    ))
    return(NULL)
  }

  if (fixed_a_result_raw <= 0) {
    if (verbose) message(sprintf(
      "[validate_fixed_lower_asymptote] = %.6f is <= 0 --- log10() undefined. Treating as NULL.",
      fixed_a_result_raw
    ))
    return(NULL)
  }

  return(fixed_a_result_raw)
}


# -----------------------------------------------------------------------------
# Section 6 -- Model constraint helpers
# -----------------------------------------------------------------------------

#' Generate Start Value Lists for NLS Optimisation
#'
#' Produces a list of lower and upper starting value vectors by spreading
#' candidate starts across the interior of the parameter bounds.
#'
#' @param bounds Named list with elements \code{lower} and \code{upper},
#'               both named numeric vectors with the same parameter names.
#' @param frac   Numeric in \code{(0, 1)}.  The fraction of the total bound
#'               width kept inside the start interval (default \code{0.90}).
#'
#' @return A list with \code{start_lower} and \code{start_upper}.
#'
#' @keywords internal
generate_start <- function(bounds, frac = 0.90) {
  start_offset <- (1 - frac) / 2
  lower <- bounds$lower
  upper <- bounds$upper

  if (!all(names(lower) == names(upper))) {
    stop("Lower and upper bounds must have identical parameter names")
  }

  width       <- upper - lower
  start_lower <- lower + start_offset * width
  start_upper <- upper - start_offset * width

  list(start_lower = start_lower, start_upper = start_upper)
}


#' Construct Named Bounds Lists
#'
#' @param param_names Character vector of parameter names.
#' @param lower_vals  Numeric vector of lower bounds (same length).
#' @param upper_vals  Numeric vector of upper bounds (same length).
#'
#' @return List with named \code{lower} and \code{upper} vectors.
#'
#' @keywords internal
.make_bounds <- function(param_names, lower_vals, upper_vals) {
  if (length(param_names) != length(lower_vals) ||
      length(param_names) != length(upper_vals)) {
    stop("Lengths must match")
  }
  list(lower = lower_vals, upper = upper_vals)
}


#' Compute Middle-90% Bounds of the Response Range
#'
#' @param ymin Numeric. Minimum observed response.
#' @param ymax Numeric. Maximum observed response.
#'
#' @return Named numeric vector \code{c(low = ..., high = ...)}.
#'
#' @keywords internal
.y_mid_bounds <- function(ymin, ymax) {
  span <- ymax - ymin
  c(low  = ymin + 0.05 * span,
    high = ymin + 0.95 * span)
}


#' Identify Free Parameters in a List of Model Formulae
#'
#' @param formulas Named list of \code{formula} objects.
#' @param dep      Character. Dependent variable name (default \code{"mfi"}).
#' @param indep    Character. Independent variable name (default
#'                 \code{"concentration"}).
#'
#' @return Named list of character vectors, one per formula, containing the
#'   free parameter names in alphabetical order.
#'
#' @export
obtain_free_variables <- function(formulas, dep = "mfi", indep = "concentration") {
  lapply(formulas, function(f) {
    vars <- all.vars(f)
    sort(setdiff(vars, c(dep, indep)))
  })
}


#' Extract the Response Variable Name from a List of Formulae
#'
#' @param formulas Named list of \code{formula} objects.
#'
#' @return Character scalar: the response (LHS) variable name.  Throws an
#'   error if the formulas disagree on the response variable.
#'
#' @export
obtain_response_variable <- function(formulas) {
  response_vars <- sapply(formulas, function(f) as.character(f[[2]]),
                          USE.NAMES = TRUE)
  unique(response_vars)
}


#' Compute Model Constraints for All Candidate Models
#'
#' Builds constraint lists (lower and upper bounds) for each of the five
#' candidate sigmoid models.  A shared \code{constraint_profile} is computed
#' once from the data and reused across all models.
#'
#' @param data                 Data frame of preprocessed standard curve data.
#' @param formulas             Named list of model formulae
#'                             (\code{logistic5, loglogistic5, logistic4, loglogistic4, gompertz4}).
#' @param response_variable    Character. Response column name.
#' @param independent_variable Character. Concentration column name.
#' @param is_log_response      Logical. Is the response log10-transformed?
#' @param is_log_concentration Logical. Is concentration log10-transformed?
#' @param antigen_settings     Named list from
#'                             \code{\link{obtain_lower_constraint}}.
#' @param max_response         Numeric. Maximum observed response.
#' @param min_response         Numeric. Minimum observed response.
#' @param verbose              Logical (default \code{TRUE}).
#'
#' @return Named list of constraint lists (one per model), with an
#'   \code{"constraint_profile"} attribute attached for downstream use.
#'
#' @export
obtain_model_constraints <- function(data, formulas,
                                     response_variable,
                                     independent_variable,
                                     is_log_response,
                                     is_log_concentration,
                                     antigen_settings,
                                     max_response,
                                     min_response,
                                     verbose = TRUE) {

  constraint_profile <- adaptive_constraint_profile(
    data               = data,
    response_variable  = response_variable,
    is_log_response    = is_log_response,
    antigen_settings   = antigen_settings
  )

  if (verbose) {
    message(sprintf(
      "[obtain_model_constraints] scale_class=%s, dynamic_range=%.3f, slope=[%.3f, %.3f], g=[%.2f, %.2f]",
      constraint_profile$scale_class, constraint_profile$dynamic_range,
      constraint_profile$slope_min, constraint_profile$slope_max,
      constraint_profile$g_min, constraint_profile$g_max
    ))
  }

  free_variables <- obtain_free_variables(formulas, dep = response_variable,
                                          indep = independent_variable)

  active_models  <- names(formulas)
  constraint_models <- list()

  if ("logistic5" %in% active_models) {
    constraint_models$logistic5 <- logistic5_safe_constraint(
      data, y_min = min_response, y_max = max_response,
      logistic5_formula   = formulas$logistic5,
      logistic5_free_vars = free_variables$logistic5,
      is_log_response      = is_log_response,
      is_log_concentration = is_log_concentration,
      antigen_settings     = antigen_settings,
      constraint_profile   = constraint_profile
    )
  }

  if ("loglogistic5" %in% active_models) {
    constraint_models$loglogistic5 <- loglogistic5_safe_constraint(
      data, y_min = min_response, y_max = max_response,
      loglogistic5_formula   = formulas$loglogistic5,
      loglogistic5_free_vars = free_variables$loglogistic5,
      is_log_response      = is_log_response,
      is_log_concentration = is_log_concentration,
      antigen_settings     = antigen_settings,
      constraint_profile   = constraint_profile
    )
  }

  if ("logistic4" %in% active_models) {
    constraint_models$logistic4 <- logistic4_safe_constraint(
      data, y_min = min_response, y_max = max_response,
      logistic4_formula   = formulas$logistic4,
      logistic4_free_vars = free_variables$logistic4,
      is_log_response      = is_log_response,
      is_log_concentration = is_log_concentration,
      antigen_settings     = antigen_settings,
      constraint_profile   = constraint_profile
    )
  }

  if ("loglogistic4" %in% active_models) {
    constraint_models$loglogistic4 <- loglogistic4_safe_constraint(
      data, y_min = min_response, y_max = max_response,
      loglogistic4_formula   = formulas$loglogistic4,
      loglogistic4_free_vars = free_variables$loglogistic4,
      is_log_response      = is_log_response,
      is_log_concentration = is_log_concentration,
      antigen_settings     = antigen_settings,
      constraint_profile   = constraint_profile
    )
  }

  if ("gompertz4" %in% active_models) {
    constraint_models$gompertz4 <- gompertz4_safe_constraint(
      data, y_min = min_response, y_max = max_response,
      gompertz4_formula   = formulas$gompertz4,
      gompertz4_free_vars = free_variables$gompertz4,
      is_log_response      = is_log_response,
      is_log_concentration = is_log_concentration,
      antigen_settings     = antigen_settings,
      constraint_profile   = constraint_profile
    )
  }


  attr(constraint_models, "constraint_profile") <- constraint_profile

  if (verbose) print(constraint_models)

  return(constraint_models)
}


#' Build Multi-Start Lists for Each Model
#'
#' Generates \code{n_starts} candidate start-value lists for each model by
#' distributing points evenly across the parameter bounds (Latin
#' hypercube-style).  For low-signal data, extra starts are added and the
#' slope parameter is biased toward smaller values.
#'
#' @param model_constraints Named list of model constraint objects as
#'                          returned by \code{\link{obtain_model_constraints}}.
#' @param quants            Named numeric vector of quantiles (unused,
#'                          retained for compatibility).
#'
#' @return Named list of lists: for each model, a list of named numeric
#'   vectors suitable for passing as \code{start} to \code{nlsLM}.
#'
#' @export
make_start_lists <- function(model_constraints,
                             quants = c(low = 0.2, mid = 0.5, high = 0.8)) {

  profile     <- attr(model_constraints, "constraint_profile")
  start_lists <- list()

  for (model_name in names(model_constraints)) {
    mc     <- model_constraints[[model_name]]
    lower  <- mc$lower
    upper  <- mc$upper
    params <- names(lower)

    n_starts <- if (!is.null(profile) && profile$scale_class == "low") 5 else 3

    starts <- list()
    for (i in seq_len(n_starts)) {
      frac    <- (i - 0.5) / n_starts
      start_i <- setNames(lower + frac * (upper - lower), params)

      if ("b" %in% params && !is.null(profile) && profile$scale_class == "low") {
        b_range      <- upper["b"] - lower["b"]
        start_i["b"] <- lower["b"] + b_range * frac^2
      }

      starts[[i]] <- as.list(start_i)
    }

    start_lists[[model_name]] <- starts
  }

  return(start_lists)
}


# -----------------------------------------------------------------------------
# Section 7 -- Model fitting
# -----------------------------------------------------------------------------

#' Fit All Candidate Sigmoid Models Using Robust Multi-Start Strategy
#'
#' Attempts to fit each of the five candidate models (logistic5, loglogistic5, logistic4, loglogistic4,
#' gompertz4) using \code{\link[minpack.lm]{nlsLM}}.  Multiple start vectors
#' are tried; the fit with the lowest AIC is retained.  For low-signal data
#' two fallback strategies are attempted when all primary starts fail:
#' \enumerate{
#'   \item Relaxed bounds (50% wider).
#'   \item Base-R \code{nls} with the \code{"port"} algorithm.
#' }
#'
#' @param prepped_data         Data frame of preprocessed standard curve data.
#' @param response_variable    Character. Response column name.
#' @param independent_variable Character. Concentration column name.
#' @param formulas             Named list of model formulae.
#' @param model_constraints    Named list from
#'                             \code{\link{obtain_model_constraints}}.
#' @param start_lists          Named list from \code{\link{make_start_lists}}.
#' @param verbose              Logical (default \code{TRUE}).
#'
#' @return Named list (one element per model).  Each element is a list with:
#'   \item{fit}{\code{nlsLM} object, or \code{NULL} if fitting failed.}
#'   \item{data}{The data used for fitting.}
#'
#' @export
compute_robust_curves <- function(prepped_data,
                                  response_variable,
                                  independent_variable,
                                  formulas,
                                  model_constraints,
                                  start_lists,
                                  verbose = TRUE) {

  profile         <- attr(model_constraints, "constraint_profile")
  is_low_response <- !is.null(profile) && profile$scale_class == "low"

  models_fit_list <- list()

  for (formula_name in names(formulas)) {
    name   <- formula_name
    lower  <- model_constraints[[name]]$lower
    upper  <- model_constraints[[name]]$upper
    starts <- start_lists[[name]]

    if (verbose) message("\n Trying model: ", name)

    best_fit <- NULL
    best_aic <- Inf

    for (sl in starts) {
      fit <- tryCatch({
        minpack.lm::nlsLM(
          formula = formulas[[name]],
          data    = prepped_data,
          start   = sl,
          lower   = lower,
          upper   = upper,
          control = nls.lm.control(
            maxiter = if (is_low_response) 200 else 100,
            ftol    = if (is_low_response) 1e-8 else 1e-6,
            ptol    = if (is_low_response) 1e-8 else 1e-6
          )
        )
      }, error = function(e) NULL)

      if (!is.null(fit)) {
        current_aic <- tryCatch(AIC(fit), error = function(e) Inf)
        if (is.finite(current_aic) && current_aic < best_aic) {
          best_aic <- current_aic
          best_fit <- fit
        }
      }
    }

    # Fallback 1: relax bounds
    if (is.null(best_fit) && is_low_response) {
      if (verbose) message("  [fallback-1] Relaxing bounds for ", name)
      relaxed_lower <- lower - 0.5 * abs(lower)
      relaxed_upper <- upper + 0.5 * abs(upper)
      if ("b" %in% names(relaxed_lower)) relaxed_lower["b"] <- max(relaxed_lower["b"], 1e-6)
      mid_start <- as.list((relaxed_lower + relaxed_upper) / 2)

      best_fit <- tryCatch({
        minpack.lm::nlsLM(
          formula = formulas[[name]], data = prepped_data,
          start   = mid_start, lower = relaxed_lower, upper = relaxed_upper,
          control = nls.lm.control(maxiter = 300, ftol = 1e-10, ptol = 1e-10)
        )
      }, error = function(e) NULL)
    }

    # Fallback 2: base nls port algorithm
    if (is.null(best_fit) && is_low_response) {
      if (verbose) message("  [fallback-2] Trying port algorithm for ", name)
      mid_start <- as.list((lower + upper) / 2)
      best_fit <- tryCatch({
        nls(
          formula   = formulas[[name]], data = prepped_data,
          start     = mid_start, algorithm = "port",
          lower     = lower, upper = upper,
          control   = nls.control(maxiter = 200, tol = 1e-6)
        )
      }, error = function(e) NULL)
    }

    if (!is.null(best_fit)) {
      models_fit_list[[name]] <- list(fit = best_fit, data = prepped_data)
      if (verbose) message("  \u2713 ", name, " converged (AIC=", round(AIC(best_fit), 2), ")")
    } else {
      models_fit_list[[name]] <- list(fit = NULL, data = prepped_data)
      if (verbose) message("  \u2717 ", name, " failed to converge")
    }
  }

  return(models_fit_list)
}


#' Select Best Multi-Start NLS Fit by AIC
#'
#' Iterates over a named list of start-value sets, fits the model for each,
#' and returns the fit with the lowest AIC.
#'
#' @param prepped_data              Data frame.
#' @param response_variable         Character. Response column name.
#' @param independent_variable      Character. Concentration column name.
#' @param formula                   A \code{formula} object.
#' @param lower_model_constraints   Named numeric vector of lower bounds.
#' @param upper_model_constraints   Named numeric vector of upper bounds.
#' @param start_lists               Named list of start-value lists.
#' @param verbose                   Logical (default \code{TRUE}).
#'
#' @return Best \code{nlsLM} fit object, or \code{NULL} if all starts fail.
#'
#' @export
select_nlsLM_aic <- function(prepped_data,
                             response_variable,
                             independent_variable,
                             formula,
                             lower_model_constraints,
                             upper_model_constraints,
                             start_lists, verbose = TRUE) {

  fits <- lapply(names(start_lists), function(nm) {
    if (verbose) message("Fitting with start: ", nm)
    tryCatch(
      nlsLM_fit(formula      = formula,
                data         = prepped_data,
                start_values = start_lists[[nm]],
                lower        = lower_model_constraints,
                upper        = upper_model_constraints,
                verbose      = verbose),
      error = function(e) {
        if (verbose) message("  Start '", nm, "' failed: ", e$message)
        NULL
      }
    )
  })

  names(fits) <- names(start_lists)
  fits <- Filter(Negate(is.null), fits)

  if (length(fits) == 0) {
    if (verbose) message("  All starts failed for this formula.")
    return(NULL)
  }

  aic_vals  <- sapply(fits, AIC)
  best_name <- names(which.min(aic_vals))

  if (verbose) {
    message("Best fit: ", best_name)
    print(aic_vals)
  }

  return(fits[[best_name]])
}


#' Fit a Single NLS Model Using minpack.lm
#'
#' Thin wrapper around \code{\link[minpack.lm]{nlsLM}} with diagnostic
#' output and graceful error handling.
#'
#' @param formula      A \code{formula} object.
#' @param data         Data frame.
#' @param start_values Named list of starting parameter values.
#' @param lower        Named numeric vector of lower bounds (default
#'                     \code{-Inf}).
#' @param upper        Named numeric vector of upper bounds (default
#'                     \code{Inf}).
#' @param verbose      Logical (default \code{TRUE}).
#'
#' @return An \code{nlsLM} object, or \code{NULL} on failure.
#'
#' @export
nlsLM_fit <- function(formula, data, start_values,
                      lower = -Inf, upper = Inf, verbose = TRUE) {

  if (verbose) {
    message("nlsLM lower constraints"); print(lower)
    message("nlsLM start values");     print(start_values)
    message("nlsLM upper constraints"); print(upper)
  }

  fit <- tryCatch({
    minpack.lm::nlsLM(
      formula = formula, data = data, start = start_values,
      lower = lower, upper = upper,
      control = nls.lm.control(maxiter = 200)
    )
  }, error = function(e) {
    if (verbose) message("nlsLM failed: ", conditionMessage(e))
    NULL
  })

  if (!is.null(fit) && verbose) {
    message("Fit successful.")
    print(fit)
  }

  return(fit)
}


# -----------------------------------------------------------------------------
# Section 8 -- Preprocessing entry point
# -----------------------------------------------------------------------------

#' Preprocess Standard Curve Data for Robust Curve Fitting
#'
#' Orchestrates the full preprocessing chain:
#' \enumerate{
#'   \item Compute concentrations from dilutions.
#'   \item Optionally apply prozone correction.
#'   \item Apply blank operation (ignore / include / subtract (1, 3, or 10 times the geometric mean of blanks)).
#'   \item Floor non-positive values and optionally log10-transform the
#'         response.
#' }
#'
#' @param data                 Data frame with columns \code{dilution} and the
#'                             response variable.
#' @param antigen_settings     Named list from
#'                             \code{\link{obtain_lower_constraint}}.
#' @param response_variable    Character. Response column name.
#' @param independent_variable Character. Concentration column name.
#' @param is_log_response      Logical. Log10-transform the response?
#' @param blank_data           Data frame of blank measurements, or
#'                             \code{NULL}.
#' @param blank_option         Character. Blank handling method.One of \code{"ignored"},
#'                           \code{"included"}, \code{"subtracted"},
#'                           \code{"subtracted_3x"}, \code{"subtracted_10x"}.
#'                           (default
#'                             \code{"ignored"}).
#' @param is_log_independent   Logical. Log10-transform the concentration?
#'                             (default \code{TRUE}).
#' @param apply_prozone        Logical. Apply prozone correction?
#'                             (default \code{TRUE}).
#' @param verbose              Logical (default \code{TRUE}).
#'
#' @return A list with:
#'   \item{data}{Preprocessed data frame.}
#'   \item{antigen_fit_options}{Named list of the options used, for passing to
#'         downstream functions.}
#'
#' @export
preprocess_robust_curves <- function(data, antigen_settings, response_variable,
                                     independent_variable,
                                     is_log_response,
                                     blank_data           = NULL,
                                     blank_option         = "ignored",
                                     is_log_independent   = TRUE,
                                     apply_prozone        = TRUE,
                                     verbose              = TRUE) {

  undiluted_sc_concentration <- get_study_exp_antigen_plate_param(antigen_settings)

  data <- compute_concentration(data, undiluted_sc_concentration, independent_variable,
                                is_log_concentration = TRUE)

  if (apply_prozone) {
    if (verbose) message("Applying prozone correction")
    data <- correct_prozone(stdframe = data, prop_diff = 0.1, dil_scale = 2,
                            response_variable    = response_variable,
                            independent_variable = independent_variable,
                            verbose              = verbose)
  }

  data <- perform_blank_operation(blank_data = blank_data, data = data,
                                  response_variable    = response_variable,
                                  independent_variable = independent_variable,
                                  is_log_response      = is_log_response,
                                  blank_option         = blank_option,
                                  verbose              = verbose)

  if (is_log_response) {
    raw_vals      <- data[[response_variable]]
    positive_vals <- raw_vals[is.finite(raw_vals) & raw_vals > 0]
    adaptive_floor <- if (length(positive_vals) > 0) min(positive_vals) * 0.01 else 1e-6

    n_floored <- sum(raw_vals <= 0, na.rm = TRUE)
    if (n_floored > 0 && verbose) {
      message(sprintf("[preprocess] %d/%d values <= 0 floored to %.2e before log10",
                      n_floored, length(raw_vals), adaptive_floor))
    }

    data[[response_variable]][is.na(raw_vals) | raw_vals <= 0] <- adaptive_floor
    data[[response_variable]] <- log10(data[[response_variable]])
  }

  antigen_fit_options <- list(
    is_log_response      = is_log_response,
    blank_option         = blank_option,
    is_log_concentration = is_log_independent,
    apply_prozone        = apply_prozone
  )

  return(list(data = data, antigen_fit_options = antigen_fit_options))
}

#' Analytic error propagation for inverse sigmoid models
#'
#' Computes the inverse-predicted concentration (`x_est`) and its propagated
#' uncertainty (`se_x`) using the delta method for a fitted nonlinear model.
#'
#' Supports multiple sigmoid model forms (4PL, 5PL, Gompertz variants), and
#' optionally handles the case where the lower asymptote (`a`) is fixed.
#'
#' @param model Character string specifying the model form.
#'   One of `"logistic4"`, `"loglogistic4"`, `"gompertz4"`, `"logistic5"`, `"loglogistic5"`.
#' @param fit A fitted `nlsLM` model object (from `minpack.lm`).
#' @param y Numeric scalar. Observed response value.
#' @param se_y Numeric scalar. Standard error of the response `y`.
#'   Defaults to `0` if unknown.
#' @param fixed_a Optional numeric scalar. Fixed lower asymptote.
#'   If supplied, `a` is treated as fixed and excluded from uncertainty propagation.
#' @param verbose Logical. If `TRUE`, prints diagnostic information.
#'
#' @return A list with components:
#' \describe{
#'   \item{x_est}{Inverse-predicted concentration}
#'   \item{se_x}{Standard error of `x_est`}
#'   \item{var_x}{Variance of `x_est`}
#'   \item{grad_theta}{Gradient w.r.t. model parameters}
#'   \item{grad_y}{Gradient w.r.t. response `y`}
#' }
#'
#' @details
#' Uses the delta method:
#' \deqn{Var(x) = \nabla_\theta^T \Sigma \nabla_\theta + (\partial x / \partial y)^2 \cdot Var(y)}
#'
#' If `fixed_a` is supplied, gradients are computed using reduced parameter space.
#'
#' @export
propagate_error_analytic <- function(model,         # character: "logistic4","loglogistic4","gompertz4","logistic5","loglogistic5"
                                     fit,           # nlsLM object (already fitted)
                                     y,             # observed response
                                     se_y = 0,     # standard error of y (0 if unknown)
                                     fixed_a,  # the constrained lower asymptote
                                     verbose = TRUE
) {
  # ----- 1. Extract coefficients & covariance -------------------------------
  theta    <- coef(fit)          # named vector
  vcov_mat <- vcov(fit)          # covariance matrix

  if(!is.null(fixed_a)) {


    # ----- 2. Analytic gradient w.r.t. parameters & y -----------------------
    inv_and_grad <- make_inv_and_grad_fixed(model, y, fixed_a)
    # ----- 3. Evaluate inverse (point estimate) ------------------------------
    x_hat <- inv_and_grad$inv(theta)
    grad_theta <- inv_and_grad$grad(theta)   # named numeric vector
    grad_y <- inv_and_grad$grad_y(theta)    # scalar
  } else {
    # ----- 2. Evaluate inverse (point estimate) ------------------------------
    x_hat <- switch(model,
                    logistic4      = inv_logistic4(y, a = theta["a"], b = theta["b"], c = theta["c"], d = theta["d"]),
                    loglogistic4     = inv_loglogistic4(y, a = theta["a"], b = theta["b"], c = theta["c"], d = theta["d"]),
                    gompertz4  = inv_gompertz4(y, a = theta["a"], b = theta["b"], c = theta["c"], d = theta["d"]),
                    logistic5      = inv_logistic5(y, a = theta["a"], b = theta["b"], c = theta["c"], d = theta["d"], g = theta["g"]),
                    loglogistic5     = inv_loglogistic5(y, a = theta["a"], b = theta["b"], c = theta["c"], d = theta["d"], g = theta["g"]),
                    stop("Unsupported model name"))

    # ----- 3. Analytic gradient w.r.t. parameters & y -----------------------
    grads <- switch(model,
                    logistic4     = grad_logistic4(y, a = theta["a"], b = theta["b"], c = theta["c"], d = theta["d"]),
                    loglogistic4    = grad_loglogistic4(y, a = theta["a"], b = theta["b"], c = theta["c"], d = theta["d"]),
                    gompertz4 = grad_gompertz4(y, a = theta["a"], b = theta["b"], c = theta["c"], d = theta["d"]),
                    logistic5     = grad_logistic5(y, a = theta["a"], b = theta["b"], c = theta["c"], d = theta["d"], g = theta["g"]),
                    loglogistic5    = grad_loglogistic5(y, a = theta["a"], b = theta["b"], c = theta["c"], d = theta["d"], g = theta["g"]))

    grad_theta <- grads$grad_theta   # named vector (same order as theta)
    grad_y     <- grads$grad_y
  }


  # ----- 4. Delta‑method variance -----------------------------------------
  var_par <- as.numeric(t(grad_theta) %*% vcov_mat %*% grad_theta)
  var_y   <- (grad_y^2) * (se_y^2)
  var_x   <- var_par + var_y
  se_x    <- sqrt(var_x)

  # ----- 5. Return ---------------------------------------------------------
  list(x_est      = x_hat,
       se_x       = se_x,
       var_x      = var_x,
       grad_theta = grad_theta,
       grad_y     = grad_y)
}

#' Diagnose inputs for analytic error propagation
#'
#' Prints diagnostic information to help debug issues in error propagation,
#' including parameter alignment, covariance structure, and gradient evaluation.
#'
#' @param fit A fitted `nlsLM` model object.
#' @param model Character string specifying the model form.
#' @param fixed_a Optional numeric scalar. Fixed lower asymptote.
#' @param y_test Optional numeric value used to test inverse and gradient functions.
#'
#' @return Invisibly returns `NULL`. Outputs diagnostic messages to console.
#'
#' @details
#' This function is useful for identifying:
#' \itemize{
#'   \item Mismatches between gradient names and covariance matrix
#'   \item Incorrect handling of fixed vs free parameters
#'   \item Failures in inverse or gradient evaluation
#' }
#'
#' @export
diagnose_propagation_inputs <- function(fit, model, fixed_a, y_test = NULL) {

  params <- coef(fit)
  Sigma  <- vcov(fit)

  cat("\n=== Propagation Input Diagnosis ===\n")
  cat("Model         :", model, "\n")
  cat("coef(fit)     :", paste(names(params), "=", round(params, 5), collapse = ", "), "\n")
  cat("vcov dim      :", paste(dim(Sigma), collapse = " x "), "\n")
  cat("vcov rownames :", paste(rownames(Sigma), collapse = ", "), "\n")
  cat("fixed_a       :", if (is.null(fixed_a)) "NULL" else round(fixed_a, 6), "\n")

  # If fixed_a is supplied, 'a' should NOT be in coef(fit)
  if (!is.null(fixed_a) && "a" %in% names(params)) {
    cat("[!]  WARNING: fixed_a is supplied BUT 'a' is ALSO in coef(fit).\n")
    cat("   This causes the augmented-Sigma path to corrupt the gradient alignment.\n")
    cat("   Solution: fit the model WITHOUT 'a' as a free parameter when fixed_a is used.\n")
  }

  # Test inv and grad at a sample y value
  if (!is.null(y_test)) {
    cat("\nTesting inv/grad at y =", y_test, "\n")
    fns <- tryCatch(
      make_inv_and_grad_fixed(model = model, y = y_test,
                              fixed_a = if (!is.null(fixed_a)) fixed_a else params["a"]),
      error = function(e) { cat("make_inv_and_grad_fixed ERROR:", e$message, "\n"); NULL }
    )
    if (!is.null(fns)) {
      x_est <- tryCatch(fns$inv(params),    error = function(e) { cat("inv() ERROR:", e$message,"\n"); NA })
      g_t   <- tryCatch(fns$grad(params),   error = function(e) { cat("grad() ERROR:", e$message,"\n"); NULL })
      g_y   <- tryCatch(fns$grad_y(params), error = function(e) { cat("grad_y() ERROR:", e$message,"\n"); NA })

      cat("  x_est    :", x_est, "\n")
      cat("  grad_t   :", if (!is.null(g_t)) paste(names(g_t), "=", round(g_t,5), collapse=", ") else "NULL", "\n")
      cat("  grad_y   :", g_y, "\n")

      # Check alignment
      if (!is.null(g_t)) {
        common <- intersect(names(g_t), rownames(Sigma))
        cat("  grad_t names  :", paste(names(g_t),   collapse=", "), "\n")
        cat("  Sigma rownames:", paste(rownames(Sigma), collapse=", "), "\n")
        cat("  Common names  :", paste(common, collapse=", "), "\n")
        if (length(common) == length(g_t) && length(common) == nrow(Sigma)) {
          g_vec  <- g_t[common]
          S_sub  <- Sigma[common, common, drop=FALSE]
          var_x  <- as.numeric(t(g_vec) %*% S_sub %*% g_vec)
          cat("  var_x (param contribution) :", var_x, "\n")
          cat("  se_x  (param contribution) :", sqrt(max(var_x, 0)), "\n")
        } else {
          cat("  [!]  NAME MISMATCH --- this is the root cause of NA se_x!\n")
        }
      }
    }
  }
  cat("===================================\n\n")
}

#' Propagate uncertainty for inverse predictions across a data frame
#'
#' Applies analytic error propagation to a vector of observed responses,
#' returning inverse-predicted concentrations and associated uncertainty.
#'
#' @param pred_df Data frame containing observed responses and their errors.
#' @param fit A fitted `nlsLM` model object.
#' @param model Character string specifying the model form.
#'   One of `"logistic4"`, `"loglogistic4"`, `"gompertz4"`, `"logistic5"`, `"loglogistic5"`.
#' @param y_col Character. Column name containing observed responses.
#' @param se_col Character. Column name containing standard errors of responses.
#' @param fixed_a Optional numeric scalar. Fixed lower asymptote.
#' @param cv_x_max Numeric. Maximum allowed coefficient of variation (CV).
#'   Values above this are capped. Default is `125`.
#' @param is_log_x Logical. Whether `x_est` is on the log10 scale.
#'   Defaults to `TRUE`.
#' @param quiet Logical. If `FALSE`, prints progress and diagnostics.
#'
#' @return The input data frame with added columns:
#' \describe{
#'   \item{predicted_concentration}{Inverse-predicted concentration}
#'   \item{se_x}{Standard error of predicted concentration}
#'   \item{cv_x}{Coefficient of variation (capped at `cv_x_max`)}
#' }
#'
#' @details
#' Uses the delta method to propagate uncertainty from:
#' \itemize{
#'   \item Model parameter covariance (`vcov(fit)`)
#'   \item Measurement error (`se_col`)
#' }
#'
#' When `is_log_x = TRUE`, CV is computed on the linear scale:
#' \deqn{CV = se_x \cdot \log(10) \cdot 100}
#'
#' This avoids instability when `x_est` is near zero on the log scale.
#'
#' @export
propagate_error_dataframe <- function(pred_df,
                                      fit,
                                      model = c("logistic4","loglogistic4","gompertz4","logistic5","loglogistic5"),
                                      y_col,
                                      se_col,
                                      fixed_a,
                                      cv_x_max = 125,
                                      is_log_x  = TRUE,   # is x_est on log10 scale?
                                      quiet = FALSE) {
  model <- match.arg(model)

  # Validate is_log_x
  is_log_x <- isTRUE(is_log_x)

  if (!quiet) {
    message("[propagate] CV formula   : ",
            if (is_log_x) "LINEAR-scale (se_x * ln(10) * 100) --- avoids /0 at log10(conc)=0"
            else           "LOG-scale    (se_x / |x_est| * 100)")
  }

  # -- 1. Validate cv_x_max ------------------------------------
  cv_x_max <- if (isTRUE(is.finite(cv_x_max)) && cv_x_max > 0) {
    as.numeric(cv_x_max)[1]
  } else {
    message("[propagate] cv_x_max invalid; defaulting to 125.")
    125
  }

  # -- 2. Extract params and Sigma from fit ---------------------
  params <- coef(fit)
  Sigma  <- vcov(fit)

  if (!quiet) {
    message("[propagate] Model       : ", model)
    message("[propagate] Free params : ", paste(names(params), collapse = ", "))
    message("[propagate] Sigma rows  : ", paste(rownames(Sigma), collapse = ", "))
    message("[propagate] fixed_a     : ",
            if (is.null(fixed_a)) "NULL (a is free in coef)" else round(as.numeric(fixed_a), 6))
  }

  # -- 3. Decide which branch make_inv_and_grad_fixed will use --
  #
  #  fixed_a supplied as a real scalar -> truly fixed, pass as-is (as.numeric)
  #  fixed_a is NULL                   -> 'a' is free, pass NULL
  #
  #  We do NOT read 'a' from params and pass it as fixed_a here.
  #  When fixed_a = NULL, make_inv_and_grad_fixed Branch B reads
  #  p["a"] internally from coef(fit) and uses the full grad_* functions.

  use_fixed_a <- !is.null(fixed_a) && isTRUE(is.finite(as.numeric(fixed_a)))
  fna         <- if (use_fixed_a) as.numeric(fixed_a) else NULL

  # Sanity: when a is free it must be in coef(fit)
  if (!use_fixed_a && !"a" %in% names(params)) {
    stop("[propagate] fixed_a is NULL but 'a' not found in coef(fit).")
  }
  # Sanity: when a is fixed it must NOT be in coef(fit)
  if (use_fixed_a && "a" %in% names(params)) {
    warning("[propagate] fixed_a supplied but 'a' also in coef(fit). ",
            "The coef 'a' will be ignored; fixed_a will be used.")
    params <- params[names(params) != "a"]
    Sigma  <- Sigma[rownames(Sigma) != "a", colnames(Sigma) != "a", drop = FALSE]
  }

  # -- 4. Loop --------------------------------------------------
  n   <- nrow(pred_df)
  res <- vector("list", n)
  pb  <- if (!quiet) txtProgressBar(min = 0, max = n, style = 3) else NULL

  n_na_inv  <- 0L
  n_na_grad <- 0L
  n_na_vpar <- 0L
  n_capped  <- 0L
  n_ok      <- 0L

  for (i in seq_len(n)) {

    y_i    <- pred_df[[y_col]][i]
    se_y_i <- if (isTRUE(is.finite(pred_df[[se_col]][i]))) pred_df[[se_col]][i] else 0

    # Scale sanity check (log10 response only --- warn once)
    if (i == 1 && !quiet) {
      y_range <- diff(range(pred_df[[y_col]], na.rm = TRUE))
      se_median <- median(pred_df[[se_col]], na.rm = TRUE)
      if (isTRUE(is.finite(se_median)) && isTRUE(is.finite(y_range)) && y_range > 0) {
        se_to_range_ratio <- se_median / y_range
        if (se_to_range_ratio > 0.5) {
          message(sprintf(paste0(
            "\n[propagate] WARNING: median se_col (%.4f) is %.1f%% of the y_col range (%.4f).\n",
            "  This suggests se_col may be on a different scale than y_col.\n",
            "  If y_col is log10-transformed, se_col must also be in log10 units.\n",
            "  Convert: se_log10 = se_raw / (raw_value * log(10))"),
            se_median, se_to_range_ratio * 100, y_range
          ))
        }
      }
    }

    # Guard: skip non-finite y
    if (!isTRUE(is.finite(y_i))) {
      res[[i]] <- list(x_est = NA_real_, se_x = NA_real_, cv_x = cv_x_max)
      n_na_inv <- n_na_inv + 1L
      if (!quiet) setTxtProgressBar(pb, i)
      next
    }

    # Build closures --- pass fna (NULL or scalar)
    fns <- tryCatch(
      make_inv_and_grad_fixed(model = model, y = y_i, fixed_a = fna),
      error = function(e) {
        if (!quiet) message(sprintf("[propagate] closure failed row %d: %s", i, e$message))
        NULL
      }
    )
    if (is.null(fns)) {
      res[[i]] <- list(x_est = NA_real_, se_x = NA_real_, cv_x = cv_x_max)
      n_na_inv <- n_na_inv + 1L
      if (!quiet) setTxtProgressBar(pb, i)
      next
    }

    # Inverse prediction
    x_est <- tryCatch(fns$inv(params), error = function(e) NA_real_)
    if (!isTRUE(is.finite(x_est))) {
      res[[i]] <- list(x_est = x_est, se_x = NA_real_, cv_x = cv_x_max)
      n_na_inv <- n_na_inv + 1L
      if (!quiet) setTxtProgressBar(pb, i)
      next
    }

    # Gradient w.r.t. free parameters dx/dtheta
    grad_t <- tryCatch(
      fns$grad(params),
      error = function(e) { n_na_grad <<- n_na_grad + 1L; rep(NA_real_, length(params)) }
    )

    # Gradient w.r.t. response dx/dy
    grad_y_val <- tryCatch(fns$grad_y(params), error = function(e) NA_real_)

    # -- Delta-method variance ---------------------------------
    # Align grad_t names with Sigma --- they MUST match now that Branch B
    # returns a,b,c,d and Branch A returns b,c,d.
    var_par <- NA_real_

    if (all(is.finite(grad_t))) {
      common <- intersect(names(grad_t), rownames(Sigma))

      if (length(common) > 0) {
        if (length(common) < length(grad_t) && !quiet && i == 1) {
          message("[propagate] NOTE: grad_t has ", length(grad_t),
                  " names but only ", length(common),
                  " align with Sigma. Using: ", paste(common, collapse = ", "))
        }
        g_sub  <- grad_t[common]
        S_sub  <- Sigma[common, common, drop = FALSE]
        var_par <- tryCatch(
          as.numeric(t(g_sub) %*% S_sub %*% g_sub),
          error = function(e) NA_real_
        )
      } else {
        # Still a mismatch --- emit a clear one-time diagnostic
        if (!quiet && i == 1) {
          message("[propagate] WARNING: zero common names between grad_t and Sigma!")
          message("  grad_t names : ", paste(names(grad_t), collapse = ", "))
          message("  Sigma rows   : ", paste(rownames(Sigma), collapse = ", "))
          message("  This means 'a' is neither in Sigma (fixed) nor returned by grad_t (free).")
          message("  Check that fixed_a is correctly NULL or a scalar.")
        }
      }
    } else {
      n_na_grad <- n_na_grad + 1L
    }

    if (!isTRUE(is.finite(var_par))) n_na_vpar <- n_na_vpar + 1L

    # Measurement-error contribution  (dx/dy)^2 * se_y^2
    var_y <- if (isTRUE(is.finite(grad_y_val)) && se_y_i > 0)
      (grad_y_val^2) * (se_y_i^2) else 0

    var_x <- if (isTRUE(is.finite(var_par))) var_par + var_y else NA_real_
    se_x  <- if (isTRUE(is.finite(var_x)) && var_x >= 0) sqrt(var_x) else NA_real_

    # CV_x
    # -- CV_x computation -----------------------------------------------------
    #
    # STRATEGY: when x_est is on the log10 concentration scale, the standard
    # ratio cv = (se_x / |x_est|) * 100 diverges as x_est -> 0 (i.e. conc -> 1).
    # This is a mathematical artefact of the log scale passing through zero ---
    # NOT a real increase in uncertainty.
    #
    # Correct approach: propagate se_x to the LINEAR concentration scale first,
    # then compute the CV there.
    #
    #   x_linear     = 10^x_est
    #   se_x_linear  = se_x * 10^x_est * log(10)      [delta method]
    #   cv_x_linear  = (se_x_linear / x_linear) * 100
    #                = se_x * log(10) * 100
    #                = se_x * 230.259...
    #
    # This is independent of x_est, so it never diverges at x_est = 0.
    # The is_log_x flag controls which formula is used.
    #
    cv_x <- if (isTRUE(is.finite(se_x))) {

      if (is_log_x) {
        # Log10-scale x_est: use linear-scale CV (avoids /0 at x_est=0)
        raw_cv <- se_x * log(10) * 100   # = se_x * 230.26
      } else {
        # Linear-scale x_est: use standard ratio CV
        if (isTRUE(abs(x_est) > 1e-10)) {
          raw_cv <- (se_x / abs(x_est)) * 100
        } else {
          raw_cv <- Inf
        }
      }

      if (isTRUE(is.finite(raw_cv))) {
        if (raw_cv < cv_x_max) {
          n_ok <- n_ok + 1L
        } else {
          n_capped <- n_capped + 1L
        }
        min(raw_cv, cv_x_max)
      } else {
        n_capped <- n_capped + 1L
        cv_x_max
      }

    } else {
      n_capped <- n_capped + 1L
      cv_x_max
    }

    res[[i]] <- list(x_est = x_est, se_x = se_x, cv_x = cv_x)
    if (!quiet) setTxtProgressBar(pb, i)
  }

  if (!quiet) close(pb)

  # -- 5. Unpack -------------------------------------------------
  pred_df$predicted_concentration <- sapply(res, `[[`, "x_est")
  pred_df$se_x                    <- sapply(res, `[[`, "se_x")
  pred_df$cv_x                    <- sapply(res, `[[`, "cv_x")
  pred_df$cv_x[!is.finite(pred_df$cv_x)] <- cv_x_max

  # -- 6. Summary ------------------------------------------------
  if (!quiet) {
    message("\n-- propagate_error_dataframe summary ------------------")
    message(sprintf("  Rows processed     : %d", n))
    message(sprintf("  x_est finite       : %d", sum(is.finite(pred_df$predicted_concentration))))
    message(sprintf("  x_est NA           : %d", sum(!is.finite(pred_df$predicted_concentration))))
    message(sprintf("  se_x finite        : %d", sum(is.finite(pred_df$se_x))))
    message(sprintf("  se_x NA            : %d  (grad or var_par issue)", sum(!is.finite(pred_df$se_x))))
    message(sprintf("  cv_x < cap         : %d", n_ok))
    message(sprintf("  cv_x at cap (%3.0f) : %d", cv_x_max, n_capped))
    message(sprintf("  grad_t NA rows     : %d", n_na_grad))
    message(sprintf("  var_par NA rows    : %d", n_na_vpar))

    xv <- pred_df$predicted_concentration[is.finite(pred_df$predicted_concentration)]
    sv <- pred_df$se_x[is.finite(pred_df$se_x)]
    cv <- pred_df$cv_x[is.finite(pred_df$cv_x) & pred_df$cv_x < cv_x_max]

    if (length(xv)) message(sprintf("  x_est range        : [%.4f, %.4f]", min(xv), max(xv)))
    if (length(sv)) message(sprintf("  se_x  range        : [%.4f, %.4f]", min(sv), max(sv)))
    if (length(cv)) message(sprintf("  cv_x  range (excl cap): [%.2f, %.2f]", min(cv), max(cv)))
    message("-------------------------------------------------------")
  }

  pred_df
}



#' Diagnose coefficient of variation (CV) behavior
#'
#' Provides summary statistics and diagnostics for propagated CV values,
#' including detection of capped values and behavior within LOQ bounds.
#'
#' @param df Data frame containing `cv_x` and `predicted_concentration`.
#' @param label Character label for the dataset (used in output messages).
#' @param lloq Optional numeric. Lower limit of quantification.
#' @param uloq Optional numeric. Upper limit of quantification.
#' @param cv_x_max Numeric. Maximum CV cap used in propagation.
#' @param verbose Logical. If `TRUE`, prints diagnostic output.
#'
#' @return Invisibly returns a list with summary statistics:
#' \describe{
#'   \item{min_cv}{Minimum CV}
#'   \item{min_x}{Concentration at minimum CV}
#'   \item{max_cv}{Maximum CV}
#'   \item{mean_cv}{Mean CV}
#'   \item{n_gt_20}{Number of CV values > 20}
#'   \item{n_at_cap}{Number of capped CV values}
#' }
#'
#' @details
#' Helps identify:
#' \itemize{
#'   \item Instability near asymptotes
#'   \item Excessive uncertainty within LOQ range
#'   \item Propagation failures leading to capped values
#' }
#'
#' @export
diagnose_cv_x <- function(df, label = "pred_se",
                          lloq = NULL, uloq = NULL,
                          cv_x_max = 125,     # cap parameter
                          verbose = TRUE) {
  if (!verbose) return(invisible(NULL))
  if (!"cv_x" %in% names(df)) {
    message(sprintf("[cv_x diagnostic] '%s': cv_x column not found.", label))
    return(invisible(NULL))
  }

  cv <- df$cv_x
  xc <- df$predicted_concentration

  finite_mask <- is.finite(cv) & is.finite(xc)
  cv_f  <- cv[finite_mask]
  xc_f  <- xc[finite_mask]

  if (length(cv_f) == 0) {
    message(sprintf("[cv_x diagnostic] '%s': no finite cv_x values.", label))
    return(invisible(NULL))
  }

  min_idx  <- which.min(cv_f)
  min_cv   <- cv_f[min_idx]
  min_x    <- xc_f[min_idx]
  max_cv   <- max(cv_f, na.rm = TRUE)
  mean_cv  <- mean(cv_f, na.rm = TRUE)

  # Count rows that hit the cap exactly (were clamped) vs genuinely > 20
  n_at_cap <- sum(cv_f >= cv_x_max, na.rm = TRUE)
  n_gt_20  <- sum(cv_f >  20,       na.rm = TRUE)
  n_na_raw <- sum(!is.finite(df$cv_x))   # residual NAs before cap applied

  message(sprintf(
    "\n[cv_x diagnostic] --- %s ---
  cv_x_max (cap)   : %.1f
  N total          : %d
  N finite cv_x    : %d
  N non-finite raw : %d  (replaced with cap)
  Min  cv_x        : %.3f  at predicted_concentration = %.4f
  Max  cv_x        : %.3f
  Mean cv_x        : %.3f
  N cv_x > 20      : %d
  N cv_x at cap    : %d",
    label, cv_x_max,
    nrow(df), length(cv_f), n_na_raw,
    min_cv, min_x, max_cv, mean_cv,
    n_gt_20, n_at_cap
  ))

  if (!is.null(lloq) && !is.null(uloq) && isTRUE(is.finite(lloq)) && isTRUE(is.finite(uloq))) {
    in_loq    <- finite_mask &
      df$predicted_concentration >= lloq &
      df$predicted_concentration <= uloq
    cv_in_loq <- df$cv_x[in_loq]
    n_loq_cap <- sum(cv_in_loq >= cv_x_max, na.rm = TRUE)

    message(sprintf(
      "  Within [lloq=%.4f, uloq=%.4f]:
    N            = %d
    mean cv_x    = %.3f
    max  cv_x    = %.3f
    N at cap     = %d  %s",
      lloq, uloq,
      length(cv_in_loq),
      mean(cv_in_loq, na.rm = TRUE),
      max(cv_in_loq,  na.rm = TRUE),
      n_loq_cap,
      if (n_loq_cap > 0)
        "[WARNING] capped values inside LOQ window --- check curve fit near limits"
      else ""
    ))
  }

  if (n_at_cap > 0) {
    message(sprintf(
      "  [INFO] %d point(s) capped at cv_x_max=%.1f (asymptote proximity or failed propagation).",
      n_at_cap, cv_x_max
    ))
  }

  invisible(list(
    min_cv   = min_cv,
    min_x    = min_x,
    max_cv   = max_cv,
    mean_cv  = mean_cv,
    n_gt_20  = n_gt_20,
    n_at_cap = n_at_cap
  ))
}


# -----------------------------------------------------------------------------
# Section 9 -- Back-calculation
# -----------------------------------------------------------------------------

#' Back-Calculate Concentration from Fitted Model for Each Sample
#'
#' Applies the analytical inverse of the selected sigmoid model to each
#' observed response value in \code{plate_samples} to estimate the
#' corresponding concentration.
#'
#' @param model_name         Character. One of \code{"logistic5"}, \code{"loglogistic5"},
#'                           \code{"logistic4"}, \code{"loglogistic4"}, \code{"gompertz4"}.
#' @param fit                Fitted \code{nlsLM} object.
#' @param plate_samples      Data frame of sample wells with the response
#'                           column.
#' @param fixed_constraint   Numeric or \code{NULL}. Fixed lower asymptote
#'                           on the model's fitting scale.
#' @param response_variable  Character. Name of the response column.
#' @param is_log_response    Logical. Was the response log10-transformed
#'                           before fitting?
#' @param verbose            Logical (default \code{TRUE}).
#'
#' @return \code{plate_samples} with a new column
#'   \code{predicted_concentration}.
#'
#' @export
calculate_predicted_concentration <- function(model_name, fit,
                                              plate_samples,
                                              fixed_constraint,
                                              response_variable,
                                              is_log_response,
                                              verbose = TRUE) {

  if (is_log_response) {
    raw_vals <- plate_samples[[response_variable]]
    if (any(raw_vals < 0, na.rm = TRUE)) {
      if (verbose) message(
        "[calculate_predicted_concentration] WARNING: negative values with is_log_response=TRUE. ",
        "Skipping log10 transform."
      )
    } else {
      plate_samples[[response_variable]] <- log10(plate_samples[[response_variable]])
    }
  }

  params <- coef(fit)
  g <- if ("g" %in% names(params)) params["g"] else 1
  b <- params["b"]; c <- params["c"]; d <- params["d"]

  if (!is.null(fixed_constraint)) {
    message("Lower asymptote is fixed at", fixed_constraint)
    a <- fixed_constraint

    plate_samples$predicted_concentration <- tryCatch({
      switch(model_name,
             logistic5     = inv_logistic5_fixed(plate_samples[[response_variable]], fixed_a=a, b=b, c=c, d=d, g=g),
             loglogistic5    = inv_loglogistic5_fixed(plate_samples[[response_variable]], fixed_a=a, b=b, c=c, d=d, g=g),
             logistic4     = inv_logistic4_fixed(plate_samples[[response_variable]],  fixed_a=a, b=b, c=c, d=d),
             loglogistic4    = inv_loglogistic4_fixed(plate_samples[[response_variable]], fixed_a=a, b=b, c=c, d=d),
             gompertz4 = inv_gompertz4_fixed(plate_samples[[response_variable]], fixed_a=a, b=b, c=c, d=d)
      )
    }, error = function(e) { message("Error: ", e$message); rep(NA_real_, nrow(plate_samples)) })

  } else {
    a <- params["a"]

    plate_samples$predicted_concentration <- tryCatch({
      switch(model_name,
             logistic5     = inv_logistic5(plate_samples[[response_variable]],  a=a, b=b, c=c, d=d, g=g),
             loglogistic5    = inv_loglogistic5(plate_samples[[response_variable]], a=a, b=b, c=c, d=d, g=g),
             logistic4     = inv_logistic4(plate_samples[[response_variable]],  a=a, b=b, c=c, d=d),
             loglogistic4    = inv_loglogistic4(plate_samples[[response_variable]], a=a, b=b, c=c, d=d),
             gompertz4 = { message("gompertz4 predicted")
               inv_gompertz4(plate_samples[[response_variable]], a=a, b=b, c=c, d=d) }
      )
    }, error = function(e) { message("Error: ", e$message); rep(NA_real_, nrow(plate_samples)) })
  }

  return(plate_samples)
}


# -----------------------------------------------------------------------------
# Section 10 -- Model selection
# -----------------------------------------------------------------------------

#' Select the Best Model Fit by AIC
#'
#' Identifies the model with the lowest AIC from a summary table and returns
#' the corresponding fit object and associated predictions.
#'
#' @param fit_summary    Data frame from \code{\link{summarize_model_fits}}.
#' @param fit_robust_lm  Named list from \code{\link{compute_robust_curves}}.
#' @param fit_params     Data frame from
#'                       \code{\link{summarize_model_parameters}}.
#' @param plot_data      Named list containing \code{pred_df}, \code{d2xy_df},
#'                       \code{dydx_df}, and optionally \code{ci_df}.
#' @param verbose        Logical (default \code{TRUE}).
#'
#' @return Named list with:
#'   \item{best_model_name}{Character. Name of the selected model.}
#'   \item{best_fit}{nlsLM fit object.}
#'   \item{best_data}{Data used to fit the selected model.}
#'   \item{best_ci}{Parameter confidence-interval data frame.}
#'   \item{best_pred}{Prediction data frame.}
#'   \item{best_d2xy}{Second-derivative data frame.}
#'   \item{best_dydx}{First-derivative data frame.}
#'   \item{best_curve_ci}{Curve confidence-interval data frame, or \code{NULL}.}
#'
#' @export
select_model_fit_AIC <- function(fit_summary,
                                 fit_robust_lm,
                                 fit_params,
                                 plot_data,
                                 verbose = TRUE) {

  selected_model_name <- fit_summary[which.min(fit_summary$AIC), ]$model
  selected_fit        <- fit_robust_lm[[selected_model_name]]$fit
  selected_data       <- fit_robust_lm[[selected_model_name]]$data
  selected_params     <- fit_params[fit_params$model == selected_model_name, ]

  pred_df     <- plot_data$pred_df[plot_data$pred_df$model == selected_model_name, ]
  d2xy_df     <- plot_data$d2xy_df[plot_data$d2xy_df$model == selected_model_name, ]
  dydx_df     <- plot_data$dydx_df[plot_data$dydx_df$model == selected_model_name, ]
  curve_ci_df <- if (!is.null(plot_data$ci_df))
    plot_data$ci_df[plot_data$ci_df$model == selected_model_name, ]
  else NULL

  list(best_model_name = selected_model_name,
       best_fit   = selected_fit,
       best_data  = selected_data,
       best_ci    = selected_params,
       best_pred  = pred_df,
       best_d2xy  = d2xy_df,
       best_dydx  = dydx_df,
       best_curve_ci = curve_ci_df)
}


# -----------------------------------------------------------------------------
# Section 11 -- Model summaries
# -----------------------------------------------------------------------------

#' Summarise Model Convergence and Information Criteria
#'
#' Returns a data frame with one row per model, reporting convergence status,
#' RSS, degrees of freedom, number of parameters, AIC, and BIC.  Models that
#' failed to converge appear with \code{converged = FALSE} and \code{NA}
#' statistics.
#'
#' @param models_fit_list Named list from \code{\link{compute_robust_curves}}.
#' @param model_names     Character vector of model names to include (default:
#'                        \code{c("logistic5","loglogistic5","logistic4","loglogistic4","gompertz4")}).
#' @param verbose         Logical (default \code{TRUE}).
#'
#' @return A data frame with columns: \code{model}, \code{converged},
#'   \code{rss}, \code{df_resid}, \code{n_params}, \code{AIC}, \code{BIC}.
#'
#' @export
summarize_model_fits <- function(models_fit_list,
                                 model_names = c("logistic5","loglogistic5","logistic4","loglogistic4","gompertz4"),
                                 verbose = TRUE) {

  all_models <- unique(c(model_names, names(models_fit_list)))

  summary_list <- lapply(all_models, function(mname) {
    fit_obj <- models_fit_list[[mname]]$fit %||% models_fit_list[[mname]]

    if (is.null(fit_obj) || !inherits(fit_obj, "nls")) {
      return(data.frame(model = mname, converged = FALSE, rss = NA_real_,
                        df_resid = NA_integer_, n_params = NA_integer_,
                        AIC = NA_real_, BIC = NA_real_,
                        stringsAsFactors = FALSE))
    }

    data.frame(
      model     = mname,
      converged = TRUE,
      rss       = tryCatch(sum(residuals(fit_obj)^2),  error = function(e) NA_real_),
      df_resid  = tryCatch(df.residual(fit_obj),        error = function(e) NA_integer_),
      n_params  = tryCatch(length(coef(fit_obj)),        error = function(e) NA_integer_),
      AIC       = tryCatch(AIC(fit_obj),                error = function(e) NA_real_),
      BIC       = tryCatch(BIC(fit_obj),                error = function(e) NA_real_),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, summary_list)
}


#' Summarise Parameter Estimates and Confidence Intervals for All Models
#'
#' Returns a tidy data frame with one row per parameter per model.
#' Confidence intervals are computed using \code{nlstools::confint2};
#' if that fails, \code{NA} is reported.
#'
#' @param models_fit_list Named list from \code{\link{compute_robust_curves}}.
#' @param level           Numeric. Confidence level (default \code{0.95}).
#' @param model_names     Character vector of model names to include.
#' @param verbose         Logical (default \code{TRUE}).
#'
#' @return Data frame with columns: \code{model}, \code{parameter},
#'   \code{estimate}, \code{conf.low}, \code{conf.high}, \code{converged}.
#'
#' @export
summarize_model_parameters <- function(models_fit_list,
                                       level       = 0.95,
                                       model_names = c("logistic5","loglogistic5","logistic4","loglogistic4","gompertz4"),
                                       verbose     = TRUE) {

  all_models <- unique(names(models_fit_list))

  summary_list <- lapply(all_models, function(mname) {
    fit_obj <- models_fit_list[[mname]]$fit %||% models_fit_list[[mname]]

    if (is.null(fit_obj) || !inherits(fit_obj, "nls")) {
      return(data.frame(model = character(0), parameter = character(0),
                        estimate = numeric(0), conf.low = numeric(0),
                        conf.high = numeric(0), converged = logical(0),
                        stringsAsFactors = FALSE))
    }

    ci <- tryCatch(nlstools::confint2(fit_obj, level = level), error = function(e) NULL)

    if (verbose) { cat("confint2 output:"); print(mname); print(ci); cat("\n\n") }

    coefs     <- coef(fit_obj)
    par_names <- names(coefs)

    if (!is.null(ci)) {
      ci        <- ci[par_names, , drop = FALSE]
      conf.low  <- ci[, 1]
      conf.high <- ci[, 2]
    } else {
      conf.low  <- rep(NA_real_, length(coefs))
      conf.high <- rep(NA_real_, length(coefs))
    }

    data.frame(model = mname, parameter = par_names,
               estimate = as.numeric(coefs),
               conf.low = as.numeric(conf.low),
               conf.high = as.numeric(conf.high),
               converged = TRUE, stringsAsFactors = FALSE)
  })

  if (verbose) message("Summarized Parameters completed")
  do.call(rbind, summary_list)
}



# -----------------------------------------------------------------------------
# Section 12 -- Consistency check
# -----------------------------------------------------------------------------

#' Check Consistency Between fixed_a_result and coef(fit)
#'
#' Validates that \code{fixed_a_result} and the free parameters in
#' \code{coef(fit)} are mutually consistent: exactly one of (a) \emph{a} is
#' fixed and absent from \code{coef(fit)}, or (b) \emph{a} is free and
#' present in \code{coef(fit)}.
#'
#' @param fit            Fitted \code{nlsLM} object.
#' @param fixed_a_result Numeric scalar or \code{NULL}.
#' @param context        Character. Passed to diagnostic messages.
#' @param verbose        Logical (default \code{TRUE}).
#'
#' @return Named list with:
#'   \item{fixed_a_result}{The (possibly corrected) fixed_a value.}
#'   \item{consistent}{Logical.}
#'   \item{correctable}{Logical. If \code{FALSE} propagation cannot proceed.}
#'
#' @export
check_fixed_a_fit_consistency <- function(fit, fixed_a_result, context = "", verbose = TRUE) {

  has_a_in_coef  <- "a" %in% names(coef(fit))
  fixed_a_is_set <- !is.null(fixed_a_result) && is.finite(as.numeric(fixed_a_result))

  consistent <- (fixed_a_is_set && !has_a_in_coef) ||   # fixed: a baked in, not in coef
    (!fixed_a_is_set && has_a_in_coef)       # free:  a estimated, in coef

  if (!consistent) {
    msg <- sprintf(
      "[check_fixed_a_fit_consistency] MISMATCH in '%s':\n  fixed_a_result = %s\n  'a' in coef(fit) = %s\n  This will cause propagation to fail.\n  Likely cause: validate_fixed_lower_asymptote() nullified fixed_a AFTER fitting.\n  Fix: validate fixed_a_result BEFORE select_model_formulas() is called.",
      context,
      if (fixed_a_is_set) format(round(as.numeric(fixed_a_result), 5)) else "NULL",
      if (has_a_in_coef) "TRUE (free parameter)" else "FALSE (baked into formula)"
    )
    warning(msg)

    # Auto-correct: if 'a' is not in coef and fixed_a is NULL,
    # this means the formula was built with a fixed value but we lost it.
    # We cannot recover the original fixed value, so return a correction instruction.
    if (!fixed_a_is_set && !has_a_in_coef) {
      if (verbose) message(
        "[check_fixed_a_fit_consistency] Cannot propagate: 'a' is neither free nor provided as fixed_a.\n",
        "  Returning fixed_a = NULL with a warning. This plate will have NA propagation results."
      )
      return(list(fixed_a_result = NULL, consistent = FALSE, correctable = FALSE))
    }

    # If fixed_a is set but 'a' is also in coef --- use fixed_a, drop 'a' from coef
    # (propagate_error_dataframe already handles this case with a warning)
    return(list(fixed_a_result = fixed_a_result, consistent = FALSE, correctable = TRUE))
  }

  return(list(fixed_a_result = fixed_a_result, consistent = TRUE, correctable = TRUE))
}

#' Ensure the response variable column exists in a data frame.
#' If the named column is missing, attempts to find it via
#' assay_response_variable metadata or common response column names.
#' Optionally coerces to numeric.
#'
#' @param df            Data frame to check
#' @param response_var  Expected column name (e.g. "mfi", "absorbance")
#' @param coerce_numeric Logical; if TRUE, coerce the column to numeric
#' @param context       Character label for diagnostic messages
#' @return A list with:
#'   \item{df}{The (possibly modified) data frame}
#'   \item{response_var}{The resolved column name (may differ from input)}
#'   \item{ok}{Logical: TRUE if a valid numeric response column was found}
ensure_response_column <- function(df,
                                   response_var,
                                   coerce_numeric = TRUE,
                                   context = "") {

  prefix <- if (nzchar(context)) paste0("[", context, "] ") else ""

  # Guard: NULL or empty data frame
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    message(sprintf("%sData frame is NULL or empty.", prefix))
    return(list(df = df, response_var = response_var, ok = FALSE))
  }

  # Case 1: Column exists by name
  if (response_var %in% names(df)) {
    if (coerce_numeric && !is.numeric(df[[response_var]])) {
      message(sprintf(
        "%sCoercing '%s' from %s to numeric.",
        prefix, response_var, class(df[[response_var]])[1]
      ))
      df[[response_var]] <- suppressWarnings(as.numeric(df[[response_var]]))
    }
    n_finite <- sum(is.finite(df[[response_var]]))
    if (n_finite == 0) {
      message(sprintf("%s'%s' exists but has 0 finite values.", prefix, response_var))
      return(list(df = df, response_var = response_var, ok = FALSE))
    }
    return(list(df = df, response_var = response_var, ok = TRUE))
  }

  # Case 2: Try assay_response_variable metadata
  if ("assay_response_variable" %in% names(df)) {
    arv <- unique(df$assay_response_variable)
    arv <- arv[!is.na(arv) & arv != ""]
    for (candidate in arv) {
      if (candidate %in% names(df)) {
        message(sprintf(
          "%s'%s' not found; using '%s' from assay_response_variable.",
          prefix, response_var, candidate
        ))
        response_var <- candidate
        if (coerce_numeric && !is.numeric(df[[response_var]])) {
          df[[response_var]] <- suppressWarnings(as.numeric(df[[response_var]]))
        }
        return(list(df = df, response_var = response_var, ok = TRUE))
      }
    }
  }

  # Case 3: Try common response column names
  common_names <- c("mfi", "absorbance", "fluorescence", "od",
                    "MFI", "Absorbance", "Fluorescence", "OD")
  found <- intersect(common_names, names(df))
  if (length(found) > 0) {
    candidate <- found[1]
    message(sprintf(
      "%s'%s' not found; falling back to '%s'.",
      prefix, response_var, candidate
    ))
    response_var <- candidate
    if (coerce_numeric && !is.numeric(df[[response_var]])) {
      df[[response_var]] <- suppressWarnings(as.numeric(df[[response_var]]))
    }
    return(list(df = df, response_var = response_var, ok = TRUE))
  }

  # Case 4: Try to extract from the NLS formula LHS
  # (If there's a formula stored somewhere, we could parse it)

  message(sprintf(
    "%sCannot find response column '%s'. Available columns: %s",
    prefix, response_var, paste(names(df), collapse = ", ")
  ))
  return(list(df = df, response_var = response_var, ok = FALSE))
}



