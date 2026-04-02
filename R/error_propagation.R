
#' Compute Median Assay SE for Each Antigen/Feature Across All Plates
#'
#' For each unique combination of study_accession, experiment_accession,
#' source, antigen, and feature, computes the standard error of assay response
#' at every dilution level across all plates, then returns the median of those
#' per-dilution SEs. This pooled median SE can be reused for error propagation
#' on each individual plate.
#'
#' @param standards_data data.frame containing all standard curve data
#' @param response_col name of the response column (e.g., "mfi")
#' @param dilution_col name of the dilution column (default = "dilution")
#' @param plate_col name of the plate identifier column (default = "plate_nom")
#' @param grouping_cols character vector of columns defining the grouping
#'        (default = c("study_accession", "experiment_accession",
#'                     "source", "antigen", "feature"))
#' @param min_reps minimum number of non-missing plate replicates required at a
#'        dilution level for that dilution's SE to be included (default = 2)
#' @param verbose logical; if TRUE emit progress messages (default = FALSE)
#'
#' @return A data.frame with one row per unique grouping containing:
#'   \item{grouping_cols}{the grouping columns}
#'   \item{median_se}{median SE across all qualifying dilution levels}
#'   \item{n_dilutions_used}{number of dilution levels with >= min_reps
#'         non-missing observations that contributed to the median}
#'   \item{n_plates}{number of distinct plates in the group}
#'   \item{total_obs}{total number of non-missing response observations used}
#'
#' @export
compute_antigen_se_table <- function(
    standards_data,
    response_col  = "mfi",
    dilution_col  = "dilution",
    plate_col     = "plate_nom",
    grouping_cols = c("study_accession",
                      "experiment_accession",
                      "source_nom",
                      "antigen",
                      "feature"),
    min_reps = 2,
    verbose  = FALSE) {

  # ------------------------------------------------------------------
  # 1. Input validation
  # ------------------------------------------------------------------
  required_cols <- unique(c(grouping_cols, response_col, dilution_col, plate_col))
  missing_cols  <- setdiff(required_cols, colnames(standards_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.numeric(standards_data[[response_col]])) {
    stop("response_col '", response_col, "' must be numeric.")
  }

  # ------------------------------------------------------------------
  # 2. Identify unique groupings
  # ------------------------------------------------------------------
  unique_groupings <- unique(standards_data[, grouping_cols, drop = FALSE])
  unique_groupings <- unique_groupings[do.call(order, unique_groupings), , drop = FALSE]
  rownames(unique_groupings) <- NULL

  if (verbose) {
    message(sprintf("Computing median SE for %d unique groupings ...",
                    nrow(unique_groupings)))
  }

  # ------------------------------------------------------------------
  # 3. Compute median SE for each grouping
  # ------------------------------------------------------------------
  se_results <- lapply(seq_len(nrow(unique_groupings)), function(i) {

    grouping <- unique_groupings[i, , drop = FALSE]

    # --- 3a. Subset to this grouping across ALL plates ----------------
    mask <- rep(TRUE, nrow(standards_data))
    for (col in grouping_cols) {
      val  <- grouping[[col]]
      mask <- mask & (!is.na(standards_data[[col]]) &
                        standards_data[[col]] == val)
    }
    grp_data <- standards_data[mask, , drop = FALSE]

    # --- 3b. Early exit if no data ------------------------------------
    if (nrow(grp_data) == 0L) {
      return(.empty_se_row(grouping, grouping_cols))
    }

    # --- 3c. Pull relevant vectors (drop rows where key cols are NA) --
    keep <- !is.na(grp_data[[dilution_col]]) &
      !is.na(grp_data[[plate_col]])
    grp_data <- grp_data[keep, , drop = FALSE]

    if (nrow(grp_data) == 0L) {
      return(.empty_se_row(grouping, grouping_cols))
    }

    dilutions  <- grp_data[[dilution_col]]
    responses  <- grp_data[[response_col]]   # may contain NA
    plates     <- grp_data[[plate_col]]

    unique_dilutions <- sort(unique(dilutions))
    n_plates         <- length(unique(plates))

    # --- 3d. SE per dilution level ------------------------------------
    # SE = sd(response across plates) / sqrt(n_non_missing)
    per_dil_se <- vapply(unique_dilutions, function(dil) {
      idx    <- dilutions == dil          # rows for this dilution
      vals   <- responses[idx]            # may include NA
      vals   <- vals[!is.na(vals)]        # drop NA responses
      n      <- length(vals)
      if (n < min_reps) return(NA_real_)  # insufficient replicates
      if (n == 1L)      return(NA_real_)  # sd undefined
      sd(vals) / sqrt(n)
    }, numeric(1L))

    # --- 3e. Median SE across dilutions (ignore NA) -------------------
    valid_se       <- per_dil_se[!is.na(per_dil_se)]
    n_dil_used     <- length(valid_se)
    total_obs      <- sum(!is.na(responses))

    median_se <- if (n_dil_used == 0L) NA_real_ else median(valid_se)

    data.frame(
      grouping,
      median_se       = median_se,
      n_dilutions_used = n_dil_used,
      n_plates        = n_plates,
      total_obs       = total_obs,
      stringsAsFactors = FALSE,
      row.names        = NULL
    )
  })

  # ------------------------------------------------------------------
  # 4. Combine and return
  # ------------------------------------------------------------------
  se_table <- do.call(rbind, se_results)
  rownames(se_table) <- NULL

  if (verbose) {
    n_valid <- sum(!is.na(se_table$median_se))
    message(sprintf(
      "Done. %d / %d groupings have a valid median SE.",
      n_valid, nrow(se_table)
    ))
  }

  return(se_table)
}



#' Look Up SE for a Specific Antigen from the SE Table
#'
#' @param se_table data.frame from compute_antigen_se_table()
#' @param study_accession study identifier
#' @param experiment_accession experiment identifier
#' @param source source identifier
#' @param antigen antigen identifier
#' @param feature feature identifier
#'
#' @return numeric SE value, or NA if not found
#' @export
lookup_antigen_se <- function(se_table,
                              study_accession,
                              experiment_accession,
                              source,
                              antigen, feature) {

  if (is.null(se_table) || nrow(se_table) == 0) {
    return(NA_real_)
  }

  # Use source_nom column if available (for ELISA wavelength support),
  # fall back to source column for backward compatibility
  source_col <- if ("source_nom" %in% names(se_table)) "source_nom" else "source"

  idx <- which(
    se_table$study_accession == study_accession &
      se_table$experiment_accession == experiment_accession &
      se_table[[source_col]] == source &
      se_table$antigen == antigen &
      se_table$feature == feature
  )

  if (length(idx) == 0) {
    return(NA_real_)
  }

  # return(se_table$overall_se[idx[1]])
  return(se_table$median_se[idx[1]])
}

#' Predict Concentrations and Propagate Uncertainty for a Single Plate
#'
#' Given a fitted standard-curve model (\code{best_fit}), this function
#' performs two related tasks:
#' \enumerate{
#'   \item Back-calculates concentrations for every point on the
#'         standards prediction grid (\code{best_fit$best_pred}).
#'   \item Back-calculates concentrations for every sample well in
#'         \code{antigen_plate$plate_samples} and attaches propagated
#'         uncertainty estimates.
#' }
#' Uncertainty is propagated via the delta method, combining parameter
#' covariance from \code{vcov(fit)} with assay measurement error
#' supplied through \code{se_std_response}.
#'
#' When \code{study_params$is_log_response} is \code{TRUE} all response
#' values are assumed to be on the \code{log10} scale inside the model.
#' \code{se_std_response} must be supplied on the \emph{raw} response
#' scale (e.g. MFI units); the function converts it to the log10 scale
#' internally using the delta-method approximation
#' \eqn{\mathrm{se}_{\log_{10}} = \sqrt{\mathrm{se}_{\mathrm{raw}} /
#' (\bar{z}_{\mathrm{ref}} \cdot \ln(10) \cdot 100)}},
#' where \eqn{\bar{z}_{\mathrm{ref}}} is the geometric mean of the raw
#' standard-curve responses on the plate.
#'
#' @param best_fit Named list returned by the model-selection step.
#'   Required elements:
#'   \describe{
#'     \item{\code{best_fit}}{The fitted model object (e.g. from
#'       \code{minpack.lm::nlsLM()}) with \code{coef()} and
#'       \code{vcov()} methods.}
#'     \item{\code{best_pred}}{data.frame of standards prediction
#'       grid; must contain a column named \code{"yhat"}.}
#'     \item{\code{best_data}}{data.frame of the standards data used
#'       to fit the model; must contain columns
#'       \code{study_accession}, \code{experiment_accession},
#'       \code{nominal_sample_dilution}, \code{plateid},
#'       \code{plate}, \code{antigen}, and \code{source}.}
#'     \item{\code{best_model_name}}{Character scalar; model identifier
#'       passed to \code{propagate_error_dataframe()} and
#'       \code{diagnose_propagation_inputs()}.}
#'     \item{\code{best_glance}}{List or data.frame with optional
#'       elements \code{lloq} and \code{uloq} (used only for
#'       diagnostics).}
#'   }
#' @param response_var Character scalar; name of the response column in
#'   \code{antigen_plate$plate_samples} and
#'   \code{antigen_plate$plate_standard} (e.g. \code{"mfi"}).
#' @param antigen_plate Named list containing plate-level data.
#'   Required elements:
#'   \describe{
#'     \item{\code{plate_samples}}{data.frame of sample wells; must
#'       contain columns \code{response_var}, \code{dilution}, and
#'       \code{well}.}
#'     \item{\code{plate_standard}}{data.frame of standard-curve wells;
#'       must contain column \code{response_var}.  Used to estimate the
#'       reference MFI for SE scale conversion.}
#'     \item{\code{fixed_a_result}}{Raw (untransformed) numeric value
#'       for the fixed lower asymptote, or \code{NULL}.  Validated and
#'       log10-transformed internally when
#'       \code{study_params$is_log_response} is \code{TRUE}.}
#'     \item{\code{antigen_settings}}{List with element
#'       \code{standard_curve_concentration}; the maximum standard
#'       concentration used to cap infinite predicted values.}
#'   }
#' @param study_params Named list of study-level modelling flags.
#'   Required elements:
#'   \describe{
#'     \item{\code{is_log_response}}{Logical; \code{TRUE} when the
#'       model was fitted on \code{log10(response)}.}
#'     \item{\code{is_log_independent}}{Logical; \code{TRUE} when the
#'       concentration axis is on the \code{log10} scale.  Passed as
#'       \code{is_log_x} to \code{propagate_error_dataframe()} and
#'       used when computing
#'       \code{final_predicted_concentration}.}
#'   }
#' @param se_std_response Numeric scalar; pooled standard error of the
#'   assay response on the \emph{raw} response scale (e.g. MFI units),
#'   typically obtained from \code{\link{compute_antigen_se_table}} via
#'   \code{\link{lookup_antigen_se}}.  Supply \code{NA} or a
#'   non-positive value to trigger the fallback SE calculation.
#' @param cv_x_max Numeric scalar; upper cap applied to all \%CV
#'   values (default \code{150}).  Passed through to
#'   \code{propagate_error_dataframe()}.
#' @param verbose Logical; if \code{TRUE} (default) emit step-by-step
#'   diagnostic messages and print the head of
#'   \code{antigen_plate$plate_samples}.
#'
#' @return The input \code{best_fit} list with two elements added or
#'   replaced:
#'   \describe{
#'     \item{\code{best_pred}}{The standards prediction data.frame
#'       augmented with \code{predicted_concentration}, \code{se_x},
#'       \code{cv_x}, \code{pcov}, and plate/study metadata columns.
#'       Infinite predicted concentrations are replaced by
#'       \code{max_conc_standard} (positive infinity) or \code{0}
#'       (negative infinity).}
#'     \item{\code{sample_se}}{data.frame with one row per sample well
#'       containing \code{final_predicted_concentration},
#'       \code{raw_predicted_concentration},
#'       \code{raw_assay_response}, \code{se_concentration},
#'       \code{cv_x}, \code{pcov}, and all original columns from
#'       \code{antigen_plate$plate_samples} (except \code{response_var}
#'       and \code{dilution}, which are carried inside
#'       \code{sample_se}).  Returns an empty data.frame if the
#'       \code{fixed_a}/\code{coef(fit)} consistency check fails and
#'       is not correctable.}
#'   }
#'
#' @seealso
#'   \code{\link{compute_antigen_se_table}},
#'   \code{\link{lookup_antigen_se}},
#'   \code{\link{propagate_error_dataframe}},
#'   \code{\link{validate_fixed_lower_asymptote}},
#'   \code{\link{check_fixed_a_fit_consistency}},
#'   \code{\link{diagnose_propagation_inputs}},
#'   \code{\link{diagnose_cv_x}}
#'
#' @export
predict_and_propagate_error <- function(best_fit,
                                        response_var,
                                        antigen_plate,
                                        study_params,
                                        se_std_response,
                                        cv_x_max = 150,
                                        verbose = TRUE) {
  if (study_params$is_log_response) {

    # Step 1: validate the raw value before any log transform
    fixed_a_result_raw <- antigen_plate$fixed_a_result
    fixed_a_result_validated <- validate_fixed_lower_asymptote(
      fixed_a_result_raw, verbose = verbose
    )

    # Step 2: only apply log10 transform if validation passed
    if (is.null(fixed_a_result_validated)) {
      fixed_a_result <- NULL
    } else {
      .eps <- 0.000005
      fixed_a_result <- log10(fixed_a_result_validated + .eps)
    }

    # -- Consistency check: fixed_a_result must agree with what's in coef(fit) --
    consistency <- check_fixed_a_fit_consistency(
      fit            = best_fit$best_fit,
      fixed_a_result = fixed_a_result,
      context        = paste("predict_and_propagate_error",
                             unique(best_fit$best_data$antigen),
                             unique(best_fit$best_data$plate)),
      verbose        = verbose
    )

    if (!consistency$consistent) {
      if (!consistency$correctable) {
        # Cannot propagate -- skip gracefully rather than crashing
        warning(sprintf(
          "[predict_and_propagate_error] Skipping propagation for antigen '%s' plate '%s': fixed_a/coef mismatch is not correctable.",
          unique(best_fit$best_data$antigen),
          unique(best_fit$best_data$plate)
        ))
        # Return best_fit with empty sample_se to avoid downstream crash
        best_fit$sample_se <- data.frame()
        return(best_fit)
      }
      fixed_a_result <- consistency$fixed_a_result
    }

    log_plate_samples <- log10(antigen_plate$plate_samples[[response_var]])
  }

  # -- Compute overall_se_value --------------------------------------------
  # se_std_response is on the RAW response scale (e.g. MFI units).
  # When is_log_response is TRUE the model and all y values are on the
  # log10 scale.  The delta method needs se_y on the SAME scale as y.
  #
  # Conversion via the log10 derivative:
  #   if Y = log10(Z)  then  dY = dZ / (Z * ln(10))
  #   => se_log10 = se_mfi / (mean_mfi * ln(10))
  #
  # We use the geometric mean of the raw responses as the reference MFI
  # for the conversion, which is appropriate for a log-transformed variable.

  if (isTRUE(is.finite(se_std_response)) && se_std_response > 0) {

    if (study_params$is_log_response) {
      # Estimate mean raw MFI from standards (geometric mean on raw scale)
      raw_standards <- antigen_plate$plate_standard[[response_var]]
      raw_standards <- raw_standards[is.finite(raw_standards) & raw_standards > 0]

      if (length(raw_standards) > 0) {
        ref_mfi <- exp(mean(log(raw_standards)))   # geometric mean
      } else {
        # Fallback: back-transform from the log10 mean of plate samples
        raw_plate <- antigen_plate$plate_samples[[response_var]]
        raw_plate  <- raw_plate[is.finite(raw_plate) & raw_plate > 0]
        ref_mfi    <- if (length(raw_plate) > 0) exp(mean(log(raw_plate))) else 1
      }

      # Convert SE from raw MFI units to log10 units
      overall_se_value <- sqrt(se_std_response / (ref_mfi * log(10) * 100))

      if (verbose) {
        message(sprintf(
          "[predict_and_propagate] SE conversion: se_mfi=%.4f, ref_mfi=%.4f -> se_log10=%.6f",
          se_std_response, ref_mfi, overall_se_value
        ))
      }

    } else {
      # Response is NOT log-transformed -- use se_std_response directly
      overall_se_value <- sqrt(se_std_response / 100)
    }

  } else {
    # se_std_response is NA / non-finite / zero -- use a small fallback
    # based on the spread of the log-scale standards if available
    if (study_params$is_log_response) {
      log_stds <- log10(antigen_plate$plate_standard[[response_var]])
      log_stds <- log_stds[is.finite(log_stds)]
      overall_se_value <- if (length(log_stds) > 1) sd(log_stds) * 0.01 else 0.01
    } else {
      overall_se_value <- 0
    }

    if (verbose) {
      message(sprintf(
        "[predict_and_propagate] se_std_response not usable (%.4f); using fallback se=%.6f",
        if (is.finite(se_std_response)) se_std_response else NA_real_,
        overall_se_value
      ))
    }
  }

  # -- Validate that overall_se_value is now on the right scale ------------
  # On log10 scale, a realistic se_y is << 1 (typically 0.01 - 0.15).
  # Warn loudly if it still looks like raw MFI units.
  if (study_params$is_log_response && overall_se_value > 1) {
    warning(sprintf(
      "[predict_and_propagate] overall_se_value=%.4f is > 1 on the log10 scale. ",
      overall_se_value,
      "This will inflate se_x. Check that se_std_response is in raw response units."
    ))
  }

  lloq  <- if (!is.null(best_fit$best_glance$lloq)) as.numeric(best_fit$best_glance$lloq)[1] else NA_real_
  uloq  <- if (!is.null(best_fit$best_glance$uloq)) as.numeric(best_fit$best_glance$uloq)[1] else NA_real_

  # -- Standards prediction curve ------------------------------------------
  pred_se <- best_fit$best_pred

  if (verbose) {
    message(paste("pred_se has", nrow(pred_se), "row(s)"))
  }

  z <- qnorm(0.975)
  max_conc_standard <- ifelse(
    study_params$is_log_independent,
    log10(antigen_plate$antigen_settings$standard_curve_concentration),
    antigen_plate$antigen_settings$standard_curve_concentration
  )

  pred_se$overall_se <- overall_se_value

  diagnose_propagation_inputs(
    fit     = best_fit$best_fit,
    model   = best_fit$best_model_name,
    fixed_a = fixed_a_result,
    y_test  = NULL
  )

  pred_se <- propagate_error_dataframe(
    pred_df  = pred_se,
    fit      = best_fit$best_fit,
    model    = best_fit$best_model_name,
    y_col    = "yhat",
    se_col   = "overall_se",
    fixed_a  = fixed_a_result,
    cv_x_max = cv_x_max,
    is_log_x = study_params$is_log_independent   # TRUE when x is log10(conc)
  )

  diagnose_cv_x(
    df       = pred_se,
    label    = "standards pred_se",
    lloq     = lloq,
    uloq     = uloq,
    cv_x_max = cv_x_max,
    verbose  = verbose
  )

  pred_se$predicted_concentration <- ifelse(
    !is.infinite(pred_se$predicted_concentration),
    pred_se$predicted_concentration,
    ifelse(pred_se$predicted_concentration > 0, max_conc_standard, 0)
  )

  pred_se$pcov                     <- pred_se$cv_x
  pred_se$study_accession          <- unique(best_fit$best_data$study_accession)
  pred_se$experiment_accession     <- unique(best_fit$best_data$experiment_accession)
  pred_se$nominal_sample_dilution  <- unique(best_fit$best_data$nominal_sample_dilution)
  pred_se$plateid                  <- unique(best_fit$best_data$plateid)
  pred_se$plate                    <- unique(best_fit$best_data$plate)
  pred_se$antigen                  <- unique(best_fit$best_data$antigen)
  pred_se$source                   <- unique(best_fit$best_data$source)

  pred_se <- attach_grouping_keys(pred_se, best_fit$best_data, context = "predict_and_propagate_error/pred_se")

  # pred_se_v <<- pred_se
  best_fit$best_pred <- pred_se

  # -- Sample propagation --------------------------------------------------
  raw_assay_response <- antigen_plate$plate_samples[[response_var]]

  if (verbose) print(head(antigen_plate$plate_samples))

  sample_se <- data.frame(
    .row_id            = seq_len(length(log_plate_samples)),
    y_new              = log_plate_samples,
    raw_assay_response = raw_assay_response,
    dilution           = antigen_plate$plate_samples$dilution,
    well               = antigen_plate$plate_samples$well
  )

  if (verbose) {
    message(paste("sample_se has", nrow(sample_se), "row(s)"))
  }

  sample_se$overall_se      <- overall_se_value   # already on log10 scale
  sample_se[[response_var]] <- sample_se$y_new

  diagnose_propagation_inputs(
    fit     = best_fit$best_fit,
    model   = best_fit$best_model_name,
    fixed_a = fixed_a_result,
    y_test  = NULL
  )

  sample_se <- propagate_error_dataframe(
    pred_df  = sample_se,
    fit      = best_fit$best_fit,
    model    = best_fit$best_model_name,
    y_col    = response_var,
    se_col   = "overall_se",
    fixed_a  = fixed_a_result,
    cv_x_max = cv_x_max,
    is_log_x = study_params$is_log_independent
  )

  diagnose_cv_x(
    df       = sample_se,
    label    = "samples sample_se",
    lloq     = lloq,
    uloq     = uloq,
    cv_x_max = cv_x_max,
    verbose  = verbose
  )

  if (study_params$is_log_independent) {
    sample_se$final_predicted_concentration <-
      10^sample_se$predicted_concentration * sample_se$dilution
  } else {
    sample_se$final_predicted_concentration <-
      sample_se$predicted_concentration * sample_se$dilution
  }

  # Build the left-hand side for the join: plate_samples minus response_var and
  # dilution (already in sample_se), plus a .row_id to guarantee 1:1 matching.
  # Joining only on "well" causes many-to-many when multiple patients/timeperiods
  # share the same well label within the plate.
  plate_samples_for_join <- antigen_plate$plate_samples[
    , !(names(antigen_plate$plate_samples) %in% c(response_var, "dilution")),
    drop = FALSE
  ]
  plate_samples_for_join$.row_id <- seq_len(nrow(plate_samples_for_join))

  sample_se <- dplyr::inner_join(
    plate_samples_for_join,
    sample_se,
    by = c(".row_id", "well")
  )
  sample_se$.row_id <- NULL

  sample_se$pcov    <- sample_se$cv_x
  sample_se$source  <- unique(best_fit$best_data$source)

  # source_nom is an internal routing column (source + wavelength composite);
  # do NOT propagate it -- source and wavelength are stored as separate DB columns.
  # wavelength and feature are handled by attach_grouping_keys below.
  sample_se <- attach_grouping_keys(sample_se, best_fit$best_data, context = "predict_and_propagate_error/sample_se")

  # Remove intermediate columns, rename se_x
  sample_se <- sample_se[, !names(sample_se) %in% c("y_new")]
  names(sample_se)[names(sample_se) == "se_x"] <- "se_concentration"

  names(sample_se)[names(sample_se) == "predicted_concentration"] <-
    "raw_predicted_concentration"

  # sample_se$raw_robust_concentration   <- NA_real_
  # sample_se$final_robust_concentration <- NA_real_
  # sample_se$se_robust_concentration    <- NA_real_
  # sample_se$pcov_robust_concentration  <- NA_real_
  #
  # sample_se_v  <<- sample_se
  best_fit$sample_se <- sample_se

  if (verbose) message("Finished predict_and_propagate_error")

  return(best_fit)
}


#' Validate a Raw Fixed Lower-Asymptote Value
#'
#' Checks that \code{fixed_a_result_raw} is a positive, finite, scalar
#' numeric before it is \code{log10}-transformed downstream. Returns
#' the original value if valid, \code{NULL} if it should be treated as free.
#'
#' @param fixed_a_result_raw Raw numeric value to validate.
#' @param verbose Logical; if \code{TRUE} emit informational messages
#'   (default \code{TRUE}).
#'
#' @return The original \code{fixed_a_result_raw} value if valid,
#'   otherwise \code{NULL}.
#' @keywords internal
validate_fixed_lower_asymptote <- function(fixed_a_result_raw, verbose = TRUE) {
  if (is.null(fixed_a_result_raw)) {
    return(NULL)
  }

  if (!is.numeric(fixed_a_result_raw) || length(fixed_a_result_raw) != 1) {
    if (verbose) message(
      "[validate_fixed_lower_asymptote] fixed_a_result is not a scalar numeric -- treating as NULL (free)."
    )
    return(NULL)
  }

  if (!is.finite(fixed_a_result_raw)) {
    if (verbose) message(sprintf(
      "[validate_fixed_lower_asymptote] fixed_a_result = %s is not finite -- treating as NULL (free).",
      as.character(fixed_a_result_raw)
    ))
    return(NULL)
  }

  if (fixed_a_result_raw <= 0) {
    if (verbose) message(sprintf(
      "[validate_fixed_lower_asymptote] fixed_a_result = %.6f is <= 0; log10() would be undefined or extreme -- treating as NULL (free lower asymptote).",
      fixed_a_result_raw
    ))
    return(NULL)
  }

  # Value is safe: positive and finite
  return(fixed_a_result_raw)
}


#' Check Consistency Between fixed_a_result and coef(fit)
#'
#' Call this just before \code{\link{propagate_error_dataframe}} to
#' catch mismatches early.
#'
#' @param fit Fitted model object with a \code{coef()} method.
#' @param fixed_a_result Numeric scalar or \code{NULL}.
#' @param context Character string describing the call site.
#' @param verbose Logical; if \code{TRUE} emit informational messages
#'   (default \code{TRUE}).
#'
#' @return A named list with elements \code{fixed_a_result},
#'   \code{consistent}, and \code{correctable}.
#' @keywords internal
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

    # If fixed_a is set but 'a' is also in coef -- use fixed_a, drop 'a' from coef
    # (propagate_error_dataframe already handles this case with a warning)
    return(list(fixed_a_result = fixed_a_result, consistent = FALSE, correctable = TRUE))
  }

  return(list(fixed_a_result = fixed_a_result, consistent = TRUE, correctable = TRUE))
}

#' Diagnose Propagation Inputs Before Error Propagation
#'
#' Emits a structured diagnostic block summarising the model
#' coefficients, variance-covariance matrix, and optionally the
#' inverse/gradient functions evaluated at a test response value.
#' Call this immediately before \code{\link{propagate_error_dataframe}}.
#'
#' @param fit Fitted model object with \code{coef()} and \code{vcov()}
#'   methods.
#' @param model Character; model name passed to
#'   \code{make_inv_and_grad_fixed()}.
#' @param fixed_a Numeric scalar or \code{NULL}.
#' @param y_test Optional numeric scalar test response value.
#'
#' @return \code{NULL} invisibly.
#' @keywords internal
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
    cat("(!)  WARNING: fixed_a is supplied BUT 'a' is ALSO in coef(fit).\n")
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
          cat("  (!)  NAME MISMATCH -- this is the root cause of NA se_x!\n")
        }
      }
    }
  }
  cat("===================================\n\n")
}

#' Diagnose CV-x Distribution After Error Propagation
#'
#' @param df data.frame containing at least \code{cv_x} and
#'   \code{predicted_concentration} columns.
#' @param label Character label used in messages (default
#'   \code{"pred_se"}).
#' @param lloq Numeric or \code{NULL}; lower limit of quantification.
#' @param uloq Numeric or \code{NULL}; upper limit of quantification.
#' @param cv_x_max Numeric cap applied to CV values (default
#'   \code{125}).
#' @param verbose Logical; if \code{FALSE} returns invisibly without
#'   printing (default \code{TRUE}).
#'
#' @return A named list with summary statistics, returned invisibly.
#' @keywords internal
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
        "[WARNING] capped values inside LOQ window -- check curve fit near limits"
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


