# ---- predict_grid.R — Predict response and back-calculate concentration + pcov
#                  across the standard grid for the best-fit model
# ----

#' Predict Grid with Uncertainty (Frequentist)
#'
#' For each point on the prediction grid: evaluates the forward model
#' to get predicted response and delta-method CI (via
#' [curveRcore::compute_curve_ci()]), then back-calculates concentration
#' and propagates uncertainty to get `se_concentration` and `pcov`.
#'
#' The `se_response` parameter controls observation-level noise in the
#' delta method. When non-zero, it adds the response variance term
#' to the concentration uncertainty via the O'Connell et al. (1993)
#' formula: `Var(x) = grad_theta' * Sigma * grad_theta + (dx/dy)^2 * sigma_y^2`.
#'
#' @param grid Data frame from [curveRcore::generate_prediction_grid()].
#'   Must contain an `x_fit` column (grid points on the fitting scale).
#' @param fit The best-fit `nlsLM` object.
#' @param model_name Character. Model name (must be in
#'   `curveRcore::available_models()`).
#' @param fixed_a Numeric or NULL. Fixed lower asymptote on fitting scale.
#'   When non-NULL the `a` parameter is removed from `theta` and `Sigma`
#'   before gradient computation.
#' @param se_response Numeric. Residual SE of the response on the fitting
#'   scale. Typically `summary(fit)$sigma`. Default 0 (parameter uncertainty
#'   only).
#' @param cv_x_max Numeric. Hard cap for `pcov` and `pcov_rmse` (%). Values
#'   above this are clamped. Default 150.
#' @param pcov_threshold Numeric. Precision threshold (%) used to set
#'   `pcov_pass`. A grid point passes when `pcov < pcov_threshold`. Default
#'   20. Must be `<= cv_x_max`.
#' @param is_log_independent Logical. Is the independent variable (x) on the
#'   log10 scale? When `TRUE`, a `log10_concentration` column is written as an
#'   alias for `x_fit`, and `pcov` is computed as
#'   `SE_x * log(10) * 100`. Default `TRUE`.
#' @param is_log_response Logical. Passed to [enrich_grid_with_d2y()].
#'   Default `TRUE`.
#' @param independent_variable Character. Name of the independent variable
#'   column in `grid`. Default `"concentration"`.
#' @param verbose Logical. Emit per-point diagnostic messages. Default
#'   `FALSE`.
#'
#' @return The input `grid` data frame with the following columns appended:
#'   \describe{
#'     \item{`log10_concentration`}{Alias for `x_fit` when
#'       `is_log_independent = TRUE`; otherwise `NA_real_`. Present
#'       unconditionally so downstream code can always reference it.}
#'     \item{`predicted_response`}{Forward model prediction
#'       \eqn{\hat{y} = f(x, \hat{\theta})} on the fitting scale.}
#'     \item{`ci_lower`, `ci_upper`}{Delta-method 95\% CI on the response
#'       scale.}
#'     \item{`predicted_concentration`}{Back-calculated concentration from
#'       \eqn{\hat{y}}, on the fitting scale.}
#'     \item{`se_concentration`}{Delta-method SE of `predicted_concentration`.}
#'     \item{`pcov`}{Percent CV of back-calculated concentration (\%), capped
#'       at `cv_x_max`.}
#'     \item{`pcov_rmse`}{Relative RMSE (\%) — incorporates both variance and
#'       bias of `predicted_concentration` vs `x_fit`. Capped at `cv_x_max`.}
#'     \item{`pcov_pass`}{Logical. `TRUE` when `pcov < pcov_threshold` and
#'       `pcov` is finite.}
#'     \item{`d2y_dx2`}{Second derivative
#'       \eqn{d^2(\log_{10} y)/d(\log_{10} x)^2} at each grid point,
#'       added by [enrich_grid_with_d2y()].}
#'   }
#'
#' @export
predict_grid_freq <- function(grid,
                              fit,
                              model_name,
                              fixed_a            = NULL,
                              se_response        = 0,
                              cv_x_max           = 150,
                              pcov_threshold     = 20,
                              is_log_independent = TRUE,
                              is_log_response    = TRUE,
                              independent_variable = "concentration",
                              verbose            = FALSE) {

  # ── Input validation ──────────────────────────────────────────────────────
  if (!is.data.frame(grid))
    stop("`grid` must be a data frame.", call. = FALSE)
  if (!"x_fit" %in% names(grid))
    stop("`grid` must contain an `x_fit` column.", call. = FALSE)
  if (!is.numeric(cv_x_max) || length(cv_x_max) != 1L || cv_x_max <= 0)
    stop("`cv_x_max` must be a single positive number.", call. = FALSE)
  if (!is.numeric(pcov_threshold) || length(pcov_threshold) != 1L ||
      pcov_threshold <= 0)
    stop("`pcov_threshold` must be a single positive number.", call. = FALSE)
  if (pcov_threshold > cv_x_max)
    warning(
      "`pcov_threshold` (", pcov_threshold, ") exceeds `cv_x_max` (",
      cv_x_max, "); every finite pcov will be marked pcov_pass = TRUE.",
      call. = FALSE
    )

  n <- nrow(grid)

  # ── log10_concentration alias ─────────────────────────────────────────────
  # Write unconditionally so downstream consumers (vignette, plots, eligibility
  # gates) can always reference the column by name. When is_log_independent is
  # FALSE x_fit is on the natural scale and the alias is set to NA_real_ with
  # a comment in the data — callers should use x_fit directly in that case.
  grid$log10_concentration <- if (is_log_independent) grid$x_fit else NA_real_

  # ── Forward prediction + CI via curveRcore ────────────────────────────────
  ci_df <- curveRcore::compute_curve_ci(
    grid                 = grid,
    model_name           = model_name,
    fit                  = fit,
    independent_variable = independent_variable
  )
  grid$predicted_response <- ci_df$yhat
  grid$ci_lower           <- ci_df$ci_lower
  grid$ci_upper           <- ci_df$ci_upper

  # ── Inverse + delta-method SE ─────────────────────────────────────────────
  theta <- stats::coef(fit)
  Sigma <- stats::vcov(fit)

  use_fixed <- !is.null(fixed_a) && is.finite(fixed_a)
  if (use_fixed && "a" %in% names(theta)) {
    theta <- theta[names(theta) != "a"]
    Sigma <- Sigma[rownames(Sigma) != "a", colnames(Sigma) != "a",
                   drop = FALSE]
  }
  fna <- if (use_fixed) fixed_a else NULL

  # Build closures ONCE outside the loop — make_inv_and_grad_fixed() returns
  # functions of (y, p) so they are reused across all grid points.
  fns <- tryCatch(
    curveRcore::make_inv_and_grad_fixed(model_name, fixed_a = fna),
    error = function(e) NULL
  )

  x_est <- se_x <- pcov_vec <- pcov_rmse_vec <- numeric(n)

  for (i in seq_len(n)) {

    y_i    <- ci_df$yhat[i]
    x_true <- grid$x_fit[i]

    # ── Propagate NULL fns (model has no analytical inverse) ──
    if (is.null(fns)) {
      x_est[i] <- se_x[i] <- pcov_vec[i] <- pcov_rmse_vec[i] <- NA_real_
      next
    }

    # ── Back-calculate concentration ──
    x_est[i] <- tryCatch(fns$inv(y_i, theta), error = function(e) NA_real_)
    if (!is.finite(x_est[i])) {
      se_x[i] <- pcov_vec[i] <- pcov_rmse_vec[i] <- NA_real_
      next
    }

    # ── Gradient w.r.t. parameters and w.r.t. response ──
    grad_t <- tryCatch(fns$grad(y_i, theta),   error = function(e) NULL)
    grad_y <- tryCatch(fns$grad_y(y_i, theta), error = function(e) NA_real_)

    # Parameter-uncertainty component: grad_theta' * Sigma * grad_theta
    var_par <- if (!is.null(grad_t) && all(is.finite(grad_t))) {
      common <- intersect(names(grad_t), rownames(Sigma))
      if (length(common) > 0L)
        as.numeric(
          t(grad_t[common]) %*%
            Sigma[common, common, drop = FALSE] %*%
            grad_t[common]
        )
      else NA_real_
    } else NA_real_

    # Observation-noise component: (dx/dy)^2 * sigma_y^2
    var_y <- if (is.finite(grad_y) && se_response > 0) {
      (grad_y^2) * (se_response^2)
    } else 0

    var_total <- if (is.finite(var_par)) var_par + var_y else NA_real_
    se_x[i]  <- if (is.finite(var_total) && var_total >= 0L) {
      sqrt(var_total)
    } else NA_real_

    # ── pcov: percent CV of back-calculated concentration ──
    pcov_vec[i] <- if (is.finite(se_x[i])) {
      raw_cv <- if (is_log_independent) {
        se_x[i] * log(10) * 100          # log10-scale: SE_x * ln10 * 100 %
      } else if (abs(x_est[i]) > 1e-10) {
        (se_x[i] / abs(x_est[i])) * 100  # natural scale: relative SE
      } else {
        Inf
      }
      min(raw_cv, cv_x_max, na.rm = TRUE)
    } else cv_x_max

    # ── pcov_rmse: relative RMSE — variance + bias of x_est vs x_true ──
    bias_sq <- (x_est[i] - x_true)^2
    mse     <- if (is.finite(var_total)) var_total + bias_sq else NA_real_
    pcov_rmse_vec[i] <- if (is.finite(mse) && mse >= 0) {
      rmse      <- sqrt(mse)
      raw_rrmse <- if (is_log_independent) {
        rmse * log(10) * 100
      } else if (abs(x_true) > 1e-10) {
        (rmse / abs(x_true)) * 100
      } else {
        Inf
      }
      min(raw_rrmse, cv_x_max, na.rm = TRUE)
    } else cv_x_max
  }

  # ── Attach results ────────────────────────────────────────────────────────
  grid$predicted_concentration <- x_est
  grid$se_concentration        <- se_x
  grid$pcov                    <- pcov_vec
  grid$pcov_rmse               <- pcov_rmse_vec

  # pcov_pass uses pcov_threshold (precision budget), NOT cv_x_max (hard cap).
  # cv_x_max is an overflow guard; marking everything below 150 % as "pass"
  # would be misleading — the usable range is defined by pcov_threshold.
  grid$pcov_pass <- !is.na(pcov_vec) & pcov_vec < pcov_threshold

  # ── Second derivative profile (shape-LOQ support) ────────────────────────
  grid <- curveRcore::enrich_grid_with_d2y(grid, is_log_response = is_log_response)

  grid
}
