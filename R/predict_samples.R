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
#' @param pcov_threshold Numeric. Precision threshold (%) used to set
#'   `pcov_pass`. A grid point passes when `pcov < pcov_threshold`. Default
#'   20. Must be `<= cv_x_max`.
#' @param is_log_independent Logical. Is x on log10 scale? Default TRUE.
#' @param lloq_log10,uloq_log10 Numeric. This curve's LLOQ/ULOQ (log10,
#'   pre-dilution) for `p_below_lloq` / `p_above_uloq`. Default `NA`.
#' @param decision_cutoff Numeric or NULL. Dilution-corrected concentration for
#'   `p_below_cutoff`. Default NULL (column is `NA`).
#' @param n_interval_draws Integer. Monte Carlo draws per sample for the
#'   predictive interval and probabilities. Default 1000; 0 disables them
#'   (columns `NA`).
#' @param interval_seed Integer or NULL. Seed for those draws. The global RNG
#'   state is never changed either way.
#' @param param_draws Matrix or NULL. Pre-drawn parameter draws (rows = draws,
#'   columns named like `coef(fit)`), as made by [fit_calibration_freq()] so the
#'   intervals and the persisted `calib_draws` share one set. NULL draws them
#'   here from `coef(fit)` / `vcov(fit)`.
#' @param verbose Logical.
#'
#' @details The interval / probability columns come from a Monte Carlo
#'   predictive distribution, not the delta method: parameters are drawn from
#'   their asymptotic multivariate normal (`coef(fit)`, `vcov(fit)`), residual
#'   noise `N(0, se_response)` is added to the observed response (the same
#'   variance term the delta method uses), and each draw is inverted with
#'   [curveRcore::invert_with_bounds()], so draws beyond an asymptote count as
#'   censored. [curveRcore::summarize_conc_draws()] then summarises them,
#'   exactly as for Bayesian fits. `se_concentration` / `pcov` are unchanged
#'   (delta method).
#'
#' @return Data frame with original sample columns plus:
#'   `observed_response_fit`, `predicted_log10_concentration`,
#'   `predicted_concentration`, `final_concentration`,
#'   `se_concentration`, `pcov`, `pcov_rmse`, and the
#'   [curveRcore::summarize_conc_draws()] columns.
#'
#' @export
predict_samples_freq <- function(samples, fit, model_name,
                                 response_variable,
                                 fixed_a = NULL,
                                 is_log_response = TRUE,
                                 se_response = 0,
                                 cv_x_max = 150,
                                 pcov_threshold = 20,
                                 is_log_independent = TRUE,
                                 lloq_log10 = NA_real_,
                                 uloq_log10 = NA_real_,
                                 decision_cutoff = NULL,
                                 n_interval_draws = 1000L,
                                 interval_seed = NULL,
                                 param_draws = NULL,
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
  samples$pcov_pass <- !is.na(pcov_vec) & pcov_vec < pcov_threshold
  # samples$pcov_pass        <- !is.na(pcov_vec) & pcov_vec < cv_x_max ### error found Sept 2026

  if (is.null(param_draws) && n_interval_draws > 0)
    param_draws <- .fit_param_draws(fit, n_interval_draws, interval_seed, fixed_a)
  lx <- .mc_conc_draws(fit_response, param_draws, fns, fna, se_response,
                       interval_seed, is_log_independent)
  if (is.null(lx)) lx <- matrix(NA_real_, nrow = n, ncol = 1L)
  summ <- curveRcore::summarize_conc_draws(lx, dilution, lloq_log10, uloq_log10,
                                           decision_cutoff)
  samples[names(summ)] <- summ
  samples
}


# Monte Carlo predictive draws of log10 concentration for each sample: the
# given MVN parameter draws P (see .fit_param_draws) plus N(0, se_response)
# response noise, inverted with out-of-range draws kept as -Inf/+Inf. Returns
# an n x S matrix, or NULL when draws are unavailable. The noise uses its own
# stream (seed + 1, so it never reuses the parameter-draw stream) and the
# caller's RNG is never moved.
.mc_conc_draws <- function(y, P, fns, fixed_a, se_response, seed, is_log_x) {
  if (is.null(fns) || is.null(P) || !nrow(P) || !length(y)) return(NULL)
  S <- nrow(P)
  noise_sd <- if (length(se_response) == 1L && is.finite(se_response) &&
                  se_response > 0) se_response else 0
  draw_e <- function() matrix(stats::rnorm(length(y) * S, 0, noise_sd), nrow = length(y))
  E <- if (is.null(seed)) withr::with_preserve_seed(draw_e())
       else withr::with_seed(seed + 1L, draw_e())

  X <- matrix(NA_real_, nrow = length(y), ncol = S)
  for (s in seq_len(S)) {
    p <- P[s, ]
    a <- if (!is.null(fixed_a)) fixed_a else p[["a"]]
    X[, s] <- curveRcore::invert_with_bounds(function(yy) fns$inv(yy, p),
                                             y + E[, s], a, p[["d"]])
  }
  X[!is.finite(y), ] <- NA_real_
  if (!is_log_x) {
    fin <- is.finite(X)
    X[fin] <- ifelse(X[fin] > 0, log10(X[fin]), -Inf)
  }
  X
}
