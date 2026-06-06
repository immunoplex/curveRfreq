# ---- select_aic.R — AIC-based model selection from the fitted ensemble
# ----

#' Select the Best Model by AIC
#'
#' Identifies the converged model with the lowest AIC and packages
#' it into the `selection` component of a `calibration_result`.
#'
#' @param ensemble Output of [fit_ensemble_nls()].
#'
#' @return A named list with:
#'   \describe{
#'     \item{best_model_name}{Character, or NA if none converged.}
#'     \item{criterion}{"AIC"}
#'     \item{weights}{Data frame with model_name, aic, delta_aic, converged.}
#'   }
#'
#' @export
select_best_aic <- function(ensemble) {
  summary_df <- summarise_ensemble(ensemble)

  converged <- summary_df[summary_df$converged, , drop = FALSE]

  if (nrow(converged) == 0) {
    return(list(
      best_model_name = NA_character_,
      criterion       = "AIC",
      weights         = summary_df[, c("model", "aic", "converged"), drop = FALSE]
    ))
  }

  best_idx  <- which.min(converged$aic)
  best_name <- converged$model[best_idx]
  best_aic  <- converged$aic[best_idx]

  # Compute delta AIC and Akaike weights
  converged$delta_aic <- converged$aic - best_aic
  converged$weight    <- exp(-0.5 * converged$delta_aic)
  converged$weight    <- converged$weight / sum(converged$weight)

  weights <- merge(
    summary_df[, c("model", "converged"), drop = FALSE],
    converged[, c("model", "aic", "delta_aic", "weight"), drop = FALSE],
    by = "model", all.x = TRUE
  )
  names(weights)[names(weights) == "model"] <- "model_name"

  list(
    best_model_name = best_name,
    criterion       = "AIC",
    weights         = weights
  )
}


#' Extract Parameters from the Best Fit
#'
#' Pulls coefficient estimates, standard errors, and confidence intervals
#' from the selected best model.
#'
#' @param ensemble Output of [fit_ensemble_nls()].
#' @param best_model_name Character. Name of the selected model.
#' @param fixed_a Numeric or NULL. Fixed lower asymptote value (on the
#'   fitting scale).
#' @param model_constraints Constraints list (for bound info).
#' @param level Confidence level. Default 0.95.
#'
#' @return Data frame with columns: term, estimate, std_error, ci_lower,
#'   ci_upper, lower_bound, upper_bound.
#' @export
extract_best_parameters <- function(ensemble, best_model_name, fixed_a = NULL,
                                    model_constraints = NULL, level = 0.95) {

  if (is.na(best_model_name) || !(best_model_name %in% names(ensemble))) {
    return(data.frame(term = character(), estimate = numeric(),
                      std_error = numeric(), ci_lower = numeric(),
                      ci_upper = numeric(), stringsAsFactors = FALSE))
  }

  fit <- ensemble[[best_model_name]]$fit
  s   <- summary(fit)
  ct  <- as.data.frame(s$coefficients)

  params_df <- data.frame(
    term      = rownames(ct),
    estimate  = ct[, "Estimate"],
    std_error = ct[, "Std. Error"],
    stringsAsFactors = FALSE
  )

  # Confidence intervals via profile or Wald
  ci <- tryCatch({
    ci_mat <- stats::confint(fit, level = level)
    data.frame(ci_lower = ci_mat[, 1], ci_upper = ci_mat[, 2])
  }, error = function(e) {
    z <- stats::qnorm(1 - (1 - level) / 2)
    data.frame(ci_lower = params_df$estimate - z * params_df$std_error,
               ci_upper = params_df$estimate + z * params_df$std_error)
  })
  params_df <- cbind(params_df, ci)

  # Append fixed_a if applicable
  if (!is.null(fixed_a)) {
    a_row <- data.frame(
      term = "a", estimate = fixed_a, std_error = 0,
      ci_lower = fixed_a, ci_upper = fixed_a,
      stringsAsFactors = FALSE
    )
    params_df <- rbind(a_row, params_df)
  }

  # Add constraint bounds
  if (!is.null(model_constraints) && best_model_name %in% names(model_constraints)) {
    mc <- model_constraints[[best_model_name]]
    bounds <- data.frame(
      term = names(mc$lower),
      lower_bound = mc$lower,
      upper_bound = mc$upper,
      stringsAsFactors = FALSE
    )
    params_df <- merge(params_df, bounds, by = "term", all.x = TRUE)
  }

  rownames(params_df) <- NULL
  params_df
}
