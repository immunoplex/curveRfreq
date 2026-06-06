# ---- fit_one_plate.R — Fit an ensemble of NLS models for one plate
#
# The heart of curveRfreq. For each model in the ensemble, tries multiple
# starting values via minpack.lm::nlsLM, keeps the best by AIC, and
# falls back to relaxed bounds or port algorithm for low-signal data.
# -----

#' Fit Ensemble of NLS Models for One Plate
#'
#' Fits each model in the ensemble using multi-start Levenberg-Marquardt
#' NLS, returning the best fit per model (by AIC). Handles convergence
#' failures gracefully — a model that fails to converge gets
#' `converged = FALSE` in the output.
#'
#' @param data Data frame of preprocessed standards (already transformed
#'   by [curveRcore::preprocess_standards()]).
#' @param formulas Named list of formula objects from [build_nls_formulas()].
#' @param model_constraints Named list from [compute_model_constraints()].
#' @param start_lists Named list from [generate_start_lists()].
#' @param verbose Logical.
#'
#' @return A named list (one element per model) of lists with:
#'   \describe{
#'     \item{model_name}{Character.}
#'     \item{converged}{Logical.}
#'     \item{fit}{The `nlsLM` object, or NULL.}
#'     \item{aic}{Numeric (Inf if failed).}
#'     \item{data}{The data used for fitting.}
#'   }
#'
#' @export
fit_ensemble_nls <- function(data, formulas, model_constraints,
                             start_lists, verbose = FALSE) {

  profile <- attr(model_constraints, "constraint_profile")
  is_low  <- !is.null(profile) && profile$scale_class == "low"

  results <- list()

  for (nm in names(formulas)) {
    lo     <- model_constraints[[nm]]$lower
    hi     <- model_constraints[[nm]]$upper
    starts <- start_lists[[nm]]

    if (verbose) message("\n Trying model: ", nm)

    best_fit <- NULL
    best_aic <- Inf

    # Primary: multi-start nlsLM
    for (sl in starts) {
      fit <- tryCatch(
        minpack.lm::nlsLM(
          formula = formulas[[nm]],
          data    = data,
          start   = sl,
          lower   = lo,
          upper   = hi,
          control = minpack.lm::nls.lm.control(
            maxiter = if (is_low) 200 else 100,
            ftol    = if (is_low) 1e-8 else 1e-6,
            ptol    = if (is_low) 1e-8 else 1e-6
          )
        ),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        a <- tryCatch(stats::AIC(fit), error = function(e) Inf)
        if (is.finite(a) && a < best_aic) {
          best_aic <- a
          best_fit <- fit
        }
      }
    }

    # Fallback 1: relaxed bounds (50% wider)
    if (is.null(best_fit) && is_low) {
      if (verbose) message("  [fallback-1] Relaxing bounds for ", nm)
      r_lo <- lo - 0.5 * abs(lo)
      r_hi <- hi + 0.5 * abs(hi)
      if ("b" %in% names(r_lo)) r_lo["b"] <- max(r_lo["b"], 1e-6)
      mid <- as.list((r_lo + r_hi) / 2)
      best_fit <- tryCatch(
        minpack.lm::nlsLM(
          formula = formulas[[nm]], data = data,
          start = mid, lower = r_lo, upper = r_hi,
          control = minpack.lm::nls.lm.control(
            maxiter = 300, ftol = 1e-10, ptol = 1e-10)
        ),
        error = function(e) NULL
      )
    }

    # Fallback 2: base R nls port algorithm
    if (is.null(best_fit) && is_low) {
      if (verbose) message("  [fallback-2] Trying port algorithm for ", nm)
      mid <- as.list((lo + hi) / 2)
      best_fit <- tryCatch(
        stats::nls(
          formula = formulas[[nm]], data = data,
          start = mid, algorithm = "port",
          lower = lo, upper = hi,
          control = stats::nls.control(maxiter = 200, tol = 1e-6)
        ),
        error = function(e) NULL
      )
    }

    converged <- !is.null(best_fit)
    if (converged) {
      best_aic <- tryCatch(stats::AIC(best_fit), error = function(e) Inf)
      if (verbose) message("  \u2713 ", nm, " converged (AIC=", round(best_aic, 2), ")")
    } else {
      if (verbose) message("  \u2717 ", nm, " failed to converge")
    }

    results[[nm]] <- list(
      model_name = nm,
      converged  = converged,
      fit        = best_fit,
      aic        = best_aic,
      data       = data
    )
  }

  results
}


#' Summarise Ensemble Fit Statistics
#'
#' Extracts AIC, BIC, RSS, df, and convergence status for each model.
#'
#' @param ensemble Output of [fit_ensemble_nls()].
#'
#' @return Data frame with one row per model.
#' @export
summarise_ensemble <- function(ensemble) {
  rows <- lapply(ensemble, function(e) {
    if (!e$converged) {
      return(data.frame(model = e$model_name, converged = FALSE,
                        rss = NA_real_, df_resid = NA_integer_,
                        n_params = NA_integer_, aic = NA_real_,
                        bic = NA_real_, stringsAsFactors = FALSE))
    }
    fit <- e$fit
    rss <- sum(stats::residuals(fit)^2)
    data.frame(
      model    = e$model_name,
      converged = TRUE,
      rss       = rss,
      df_resid  = stats::df.residual(fit),
      n_params  = length(stats::coef(fit)),
      aic       = stats::AIC(fit),
      bic       = stats::BIC(fit),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
