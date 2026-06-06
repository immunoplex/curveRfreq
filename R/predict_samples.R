# ---- predict_samples.R — Back-calculate concentration for test samples
# ----

#' Predict Concentration for Test Samples
#'
#' For each test sample, applies the inverse of the best-fit model to
#' the observed response to get predicted concentration, then propagates
#' uncertainty via the delta method to get `se_concentration` and `pcov`.
#'
#' @param samples Data frame of test samples. Must contain the response
#'   column and `dilution`.
#' @param fit The best-fit `nlsLM` object.
#' @param model_name Character. Best model name.
#' @param response_variable Character. Response column name.
#' @param fixed_a Numeric or NULL. Fixed lower asymptote on fitting scale.
#' @param is_log_response Logical. Was the response log10-transformed?
#' @param se_response Numeric. Residual SE on the fitting scale.
#'   Typically `summary(fit)$sigma`. Default 0.
#' @param cv_x_max Numeric. Cap for pcov. Default 150.
#' @param is_log_independent Logical. Is x on log10 scale? Default TRUE.
#' @param verbose Logical.
#'
#' @return Data frame with original sample columns plus:
#'   `observed_response_fit`, `predicted_log10_concentration`,
#'   `predicted_concentration`, `final_concentration`,
#'   `se_concentration`, `pcov`, `pcov_rmse`.
#'
#' @export
predict_samples_freq <- function(samples, fit, model_name,
                                 response_variable,
                                 fixed_a = NULL,
                                 is_log_response = TRUE,
                                 se_response = 0,
                                 cv_x_max = 150,
                                 is_log_independent = TRUE,
                                 verbose = FALSE) {

  n <- nrow(samples)
  if (n == 0) return(samples)

  theta <- stats::coef(fit)
  Sigma <- stats::vcov(fit)

  use_fixed <- !is.null(fixed_a) && is.finite(fixed_a)
  if (use_fixed && "a" %in% names(theta)) {
    theta <- theta[names(theta) != "a"]
    Sigma <- Sigma[rownames(Sigma) != "a", colnames(Sigma) != "a", drop = FALSE]
  }

  # Get raw response and transform if needed
  raw_response <- samples[[response_variable]]
  if (is_log_response) {
    fit_response <- log10(pmax(raw_response, 1e-6))
  } else {
    fit_response <- raw_response
  }

  fna <- if (use_fixed) fixed_a else NULL

  # ── Build closures ONCE outside the loop ──
  fns <- tryCatch(
    curveRcore::make_inv_and_grad_fixed(model_name, fixed_a = fna),
    error = function(e) NULL
  )

  x_est <- se_x <- pcov_vec <- pcov_rmse_vec <- numeric(n)

  for (i in seq_len(n)) {
    y_i    <- fit_response[i]
    se_y_i <- se_response

    if (!is.finite(y_i)) {
      x_est[i] <- se_x[i] <- NA_real_
      pcov_vec[i] <- pcov_rmse_vec[i] <- cv_x_max
      next
    }

    if (is.null(fns)) {
      x_est[i] <- se_x[i] <- NA_real_
      pcov_vec[i] <- pcov_rmse_vec[i] <- cv_x_max
      next
    }

    x_est[i] <- tryCatch(fns$inv(y_i, theta), error = function(e) NA_real_)
    if (!is.finite(x_est[i])) {
      se_x[i] <- NA_real_
      pcov_vec[i] <- pcov_rmse_vec[i] <- cv_x_max
      next
    }

    grad_t <- tryCatch(fns$grad(y_i, theta),   error = function(e) NULL)
    grad_y <- tryCatch(fns$grad_y(y_i, theta), error = function(e) NA_real_)

    var_par <- if (!is.null(grad_t) && all(is.finite(grad_t))) {
      common <- intersect(names(grad_t), rownames(Sigma))
      if (length(common) > 0)
        as.numeric(t(grad_t[common]) %*%
                     Sigma[common, common, drop = FALSE] %*%
                     grad_t[common])
      else NA_real_
    } else NA_real_

    var_y <- if (is.finite(grad_y) && se_y_i > 0) (grad_y^2) * (se_y_i^2) else 0
    var_total <- if (is.finite(var_par)) var_par + var_y else NA_real_
    se_x[i] <- if (is.finite(var_total) && var_total >= 0) sqrt(var_total) else NA_real_

    # pcov: CV (%)
    pcov_vec[i] <- if (is.finite(se_x[i])) {
      raw_cv <- if (is_log_independent) se_x[i] * log(10) * 100
      else if (abs(x_est[i]) > 1e-10) (se_x[i] / abs(x_est[i])) * 100
      else Inf
      min(raw_cv, cv_x_max, na.rm = TRUE)
    } else cv_x_max

    # pcov_rmse: for samples, bias is not defined against a known truth,
    # so pcov_rmse equals pcov (SE = RMSE when bias = 0 by construction)
    pcov_rmse_vec[i] <- pcov_vec[i]
  }

  # Build output
  samples$raw_assay_response            <- raw_response
  samples$observed_response_fit         <- fit_response
  samples$predicted_log10_concentration <- if (is_log_independent) x_est else log10(pmax(x_est, 1e-20))
  samples$predicted_concentration       <- if (is_log_independent) x_est else x_est

  # final_concentration = back-transform × dilution
  dilution <- if ("dilution" %in% names(samples)) samples$dilution else 1
  if (is_log_independent) {
    samples$final_concentration <- 10^x_est * dilution
  } else {
    samples$final_concentration <- x_est * dilution
  }

  samples$se_concentration <- se_x
  samples$pcov             <- pcov_vec
  samples$pcov_rmse        <- pcov_rmse_vec
  samples$pcov_pass        <- !is.na(pcov_vec) & pcov_vec < cv_x_max

  samples
}
