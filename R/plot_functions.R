#' Compute delta-method 95% confidence interval bands for a fitted curve
#'
#' Calculates confidence interval bands for a nonlinear least squares fitted
#' curve using the delta method. Performs exactly \code{p_free} vectorised
#' formula evaluations over the full \code{x_new} grid (one per free
#' parameter) rather than the naive \code{n × p_free} scalar
#' \code{predict()} calls. For a 200-point grid and a 5-parameter model
#' this is approximately 40x fewer evaluations.
#'
#' @param fit A converged \code{nls} or \code{nlsLM} object.
#' @param x_new Numeric vector of x values at which to evaluate the curve.
#' @param x_var Character string naming the independent variable in the model
#'   formula.
#' @param fixed_a Numeric scalar or \code{NULL}. Fixed lower asymptote not
#'   included in \code{vcov()}. Default is \code{NULL}.
#' @param level Numeric confidence level. Default is \code{0.95}.
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{x}{The x values supplied in \code{x_new}.}
#'     \item{ci_lo}{Lower confidence band.}
#'     \item{ci_hi}{Upper confidence band.}
#'   }
#'   If \code{vcov(fit)} fails, \code{ci_lo} and \code{ci_hi} are \code{NA}.
#'
#' @examples
#' \dontrun{
#' ci <- compute_curve_ci(fit = my_nls_fit, x_new = seq(0, 10, length.out = 200),
#'                        x_var = "concentration", fixed_a = NULL)
#' }
#'
#' @seealso \code{\link[stats]{nls}}, \code{\link[stats]{vcov}}
#' @export
compute_curve_ci <- function(fit, x_new, x_var, fixed_a = NULL, level = 0.95) {
  n      <- length(x_new)
  nd     <- setNames(data.frame(x_new), x_var)
  theta  <- coef(fit)                      # free parameters only
  V      <- tryCatch(vcov(fit), error = function(e) NULL)
  if (is.null(V)) return(data.frame(x = x_new, ci_lo = NA_real_, ci_hi = NA_real_))

  p_free <- length(theta)
  rhs    <- as.list(formula(fit))[[3]]     # model formula RHS

  # Baseline predictions via predict() for reliability
  yhat <- tryCatch(
    as.numeric(predict(fit, newdata = nd)),
    error = function(e) rep(NA_real_, n)
  )

  # Reference predictions via formula eval (used as base for finite differences)
  base_env <- c(as.list(theta),
                if (!is.null(fixed_a)) list(a = fixed_a) else list(),
                as.list(nd))
  yhat_ref <- tryCatch(
    as.numeric(eval(rhs, envir = base_env)),
    error = function(e) yhat
  )

  # Jacobian (n x p_free): one vectorised formula eval per free parameter
  # sqrt(.Machine$double.eps) is the near-optimal step for forward differences
  eps <- sqrt(.Machine$double.eps)
  J <- matrix(NA_real_, nrow = n, ncol = p_free)
  for (j in seq_len(p_free)) {
    theta_j    <- theta
    theta_j[j] <- theta[j] + eps
    env_j <- c(as.list(theta_j),
               if (!is.null(fixed_a)) list(a = fixed_a) else list(),
               as.list(nd))
    yhat_j <- tryCatch(
      as.numeric(eval(rhs, envir = env_j)),
      error = function(e) rep(NA_real_, n)
    )
    J[, j] <- (yhat_j - yhat_ref) / eps
  }

  # SE via efficient diag(J V J') = rowSums((J %*% V) * J)
  JV      <- J %*% V
  se_pred <- sqrt(pmax(rowSums(JV * J), 0))

  df_resid <- max(length(residuals(fit)) - p_free, 1)
  t_crit   <- qt((1 + level) / 2, df = df_resid)

  data.frame(
    x     = x_new,
    ci_lo = yhat - t_crit * se_pred,
    ci_hi = yhat + t_crit * se_pred
  )
}

#' Build prediction and diagnostic data for model comparison plots
#'
#' For each model in \code{model_names}, generates predicted values on a
#' fine grid, residuals, first and second derivatives, and delta-method
#' confidence intervals. Appends AIC values to the parameter summary for
#' downstream plotting.
#'
#' @param models_fit_list Named list of model fit objects (or lists with a
#'   \code{$fit} element). Names should correspond to \code{model_names}.
#' @param prepped_data A \code{data.frame} containing at minimum columns
#'   named by \code{x_var} and \code{y_var}.
#' @param fit_params A \code{data.frame} of parameter estimates with columns
#'   \code{model}, \code{parameter}, \code{estimate}, \code{conf.low}, and
#'   \code{conf.high}, as returned by \code{summarize_model_fits()}.
#' @param fixed_a_result Numeric scalar or \code{NULL}. Fixed lower asymptote
#'   passed to \code{compute_curve_ci()} and used to reconstruct full
#'   coefficient vectors.
#' @param model_names Character vector of model names to process. Must match
#'   keys in \code{models_fit_list}. Default is
#'   \code{c("Y5","Yd5","Y4","Yd4","Ygomp4")}.
#' @param x_var Character string naming the independent variable column in
#'   \code{prepped_data}. Default is \code{"concentration"}.
#' @param y_var Character string naming the response variable column in
#'   \code{prepped_data}. Default is \code{"mfi"}.
#' @param verbose Logical. If \code{TRUE}, prints a completion message.
#'   Default is \code{TRUE}.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{dat}{\code{data.frame} of the input data with added \code{x} and
#'       \code{y} columns.}
#'     \item{pred_df}{\code{data.frame} of predicted values on the fine grid,
#'       with columns \code{model}, \code{x}, \code{yhat}.}
#'     \item{resid_df}{\code{data.frame} of residuals, with columns
#'       \code{model}, \code{fitted}, \code{residuals}.}
#'     \item{fit_params}{\code{data.frame} of parameter estimates as supplied
#'       in \code{fit_params}.}
#'     \item{d2xy_df}{\code{data.frame} of second-derivative values, with
#'       columns \code{model}, \code{x}, \code{d2x_y}.}
#'     \item{dydx_df}{\code{data.frame} of first-derivative values, with
#'       columns \code{model}, \code{x}, \code{dydx}.}
#'     \item{ci_df}{\code{data.frame} of confidence intervals, with columns
#'       \code{model}, \code{x}, \code{ci_lo}, \code{ci_hi}.}
#'     \item{fit_params_aic}{\code{data.frame} combining \code{fit_params}
#'       with AIC values for each converged model.}
#'   }
#'
#' @examples
#' \dontrun{
#' plot_data <- get_plot_data(
#'   models_fit_list = my_fits,
#'   prepped_data    = my_data,
#'   fit_params      = my_params,
#'   fixed_a_result  = NULL
#' )
#' }
#'
#' @seealso \code{\link{compute_curve_ci}}, \code{\link{plot_model_comparisons}}
#' @export
get_plot_data <- function(models_fit_list,
                          prepped_data,
                          fit_params,
                          fixed_a_result,
                          model_names = c("Y5","Yd5","Y4","Yd4","Ygomp4"),
                          x_var = "concentration",
                          y_var = "mfi",
                          verbose = TRUE) {
  # Extract data
  dat <- prepped_data
  dat$x <- dat[[x_var]]
  dat$y <- dat[[y_var]]
  x <- dat[[x_var]]
  y <- dat[[y_var]]

  ## 1. Build prediction data for all models
  # A fine grid over the range of x for smooth curves
  x_new <- seq(min(x, na.rm = TRUE),
               max(x, na.rm = TRUE),
               length.out = 200)

  pred_list <- list()
  resid_list <- list()
  d2xy_list <- list()
  dydx_list <- list()
  ci_list   <- list()

  for (mname in model_names) {
    entry <- models_fit_list[[mname]]
    fit_obj <- if (is.list(entry)) entry$fit else entry

    if (is.null(fit_obj) || !inherits(fit_obj, "nls")) {
      next
    }

    # prediction data.frame must have same column name as x_var
    newdata <- data.frame(x_new)
    names(newdata) <- x_var

    # Predictions + residuals
    y_pred <- tryCatch({
      predict(fit_obj, newdata = newdata)
    }, error = function(e) rep(NA_real_, length(x_new)))

    if (is.null(fixed_a_result)) {
      fit_obj_coef <- coef(fit_obj)
    } else {
      fit_obj_coef <- c(a = fixed_a_result, coef(fit_obj))
    }

    d2x_y <- tryCatch({
      if (mname == "Y5") {
        do.call(d2xY5, c(list(x = x_new), fit_obj_coef))
      } else if (mname == "Y4") {
        do.call(d2xY4, c(list(x = x_new), fit_obj_coef))
      } else if (mname == "Yd5") {
        do.call(d2xYd5, c(list(x = x_new), fit_obj_coef))
      } else if (mname == "Yd4") {
        do.call(d2xYd4, c(list(x = x_new), fit_obj_coef))
      } else if (mname == "Ygomp4") {
        do.call(d2xYgomp4, c(list(x = x_new), fit_obj_coef))
      }
    }, error = function(e) rep (NA_real_, length(x_new)))

    dydx <- tryCatch({
      if (mname == "Y5") {
        do.call(dydxY5, c(list(x = x_new), fit_obj_coef))
      } else if (mname == "Y4") {
        do.call(dydxY4, c(list(x = x_new), fit_obj_coef))
      } else if (mname == "Yd5") {
        do.call(dydxYd5, c(list(x = x_new), fit_obj_coef))
      } else if (mname == "Yd4") {
        do.call(dydxYd4, c(list(x = x_new), fit_obj_coef))
      } else if (mname == "Ygomp4") {
        do.call(dydxYgomp4, c(list(x = x_new), fit_obj_coef))
      }
    }, error = function(e) rep (NA_real_, length(x_new)))

    y_fit <- tryCatch({
      fitted(fit_obj)
    }, error = function(e) rep(NA_real_, length(x)))

    resid <- tryCatch({
      residuals(fit_obj)
    }, error = function(e) rep(NA_real_, length(x)))

    pred_list[[mname]] <- data.frame(
      model = mname,
      x     = x_new,
      yhat  = y_pred
    )

    resid_list[[mname]] <- data.frame(
      model     = mname,
      fitted    = y_fit,
      residuals = resid
    )

    d2xy_list[[mname]] <- data.frame(
      model = mname,
      x = x_new,
      d2x_y = d2x_y
    )

    dydx_list[[mname]] <- data.frame(
      model = mname,
      x = x_new,
      dydx = dydx
    )

    ci_raw <- tryCatch(
      compute_curve_ci(fit_obj, x_new, x_var, fixed_a = fixed_a_result),
      error = function(e) data.frame(x = x_new, ci_lo = NA_real_, ci_hi = NA_real_)
    )
    ci_list[[mname]] <- data.frame(model = mname, ci_raw)
  }


  pred_df  <- do.call(rbind, pred_list)
  pred_df$yhat <- as.numeric(pred_df$yhat)
  resid_df <- do.call(rbind, resid_list)
  d2xy_df  <- do.call(rbind, d2xy_list)
  dydx_df  <- do.call(rbind, dydx_list)
  ci_df    <- do.call(rbind, ci_list)


  fit_summary <- summarize_model_fits(models_fit_list, model_names)
  fit_summary_long <- reshape2::melt(
    fit_summary,
    id.vars = c("model", "converged"),
    measure.vars = c("AIC"),
    variable.name = "criterion",
    value.name = "value"
  )
  fit_summary_long <- subset(fit_summary_long, converged & is.finite(value))
  names(fit_summary_long)[names(fit_summary_long) == "criterion"] <- "parameter"
  names(fit_summary_long)[names(fit_summary_long) == "value"] <- "estimate"
  fit_summary_long$conf.low <- fit_summary_long$estimate
  fit_summary_long$conf.high <- fit_summary_long$estimate

  fit_params_aic <- rbind(fit_summary_long, fit_params)
  # ── Drop rows where parameter is NA (from non-converged models) ──
  fit_params_aic <- fit_params_aic[!is.na(fit_params_aic$parameter), , drop = FALSE]

  if (verbose) {
    message("Plot Data Completed")
  }
  return(list(dat = dat, pred_df = pred_df,
              resid_df = resid_df,
              fit_params = fit_params,
              d2xy_df = d2xy_df,
              dydx_df = dydx_df,
              ci_df = ci_df,
              fit_params_aic = fit_params_aic))
}



#' Plot comparison of nonlinear model fits
#'
#' Produces a three-panel \pkg{patchwork} figure (or a list of individual
#' \pkg{ggplot2} objects) showing (A) observed data with fitted curves,
#' (B) residuals vs fitted values, and (C) parameter estimates with
#' confidence intervals faceted by parameter name.
#'
#' @param plot_data Named list as returned by \code{\link{get_plot_data}}.
#' @param model_names Character vector of model names to include in the plots.
#'   Default is \code{c("Y5","Yd5","Y4","Yd4","Ygomp4")}.
#' @param x_var Character string naming the independent variable. Used for
#'   axis labelling after passing through \code{\link{format_assay_terms}}.
#'   Default is \code{"concentration"}.
#' @param y_var Character string naming the response variable. Used for axis
#'   labelling. Default is \code{"mfi"}.
#' @param is_display_log_response Logical. If \code{TRUE}, the y-axis label
#'   is formatted as \eqn{\log_{10}(\text{y\_var})}. Default is \code{TRUE}.
#' @param is_display_log_independent Logical. If \code{TRUE}, the x-axis
#'   label is formatted as \eqn{\log_{10}(\text{x\_var})}. Default is
#'   \code{TRUE}.
#' @param use_patchwork Logical. If \code{TRUE} (default), returns a single
#'   combined \pkg{patchwork} figure. If \code{FALSE}, returns a named list
#'   of individual \pkg{ggplot2} objects (\code{data_fit}, \code{resid},
#'   \code{p_ci_params}).
#'
#' @return When \code{use_patchwork = TRUE}, a \pkg{patchwork} object
#'   combining all three panels. When \code{use_patchwork = FALSE}, a named
#'   list with elements \code{data_fit}, \code{resid}, and \code{p_ci_params}.
#'
#' @examples
#' \dontrun{
#' plot_model_comparisons(plot_data = my_plot_data,
#'                        is_display_log_response    = TRUE,
#'                        is_display_log_independent = TRUE)
#' }
#'
#' @seealso \code{\link{get_plot_data}}, \code{\link{format_assay_terms}}
#' @export
plot_model_comparisons <- function(plot_data,
                                   model_names = c("Y5","Yd5","Y4","Yd4","Ygomp4"),
                                   x_var = "concentration",
                                   y_var = "mfi",
                                   is_display_log_response = TRUE,
                                   is_display_log_independent = TRUE,
                                   use_patchwork = TRUE) {
  ## 1.Extract data
  pred_df <- plot_data$pred_df
  resid_df <- plot_data$resid_df
  fit_params_df <- plot_data$fit_params
  fit_params_aic <- plot_data$fit_params_aic
  dat <- plot_data$dat
  x <- dat[[x_var]]
  y <- dat[[y_var]]

  ## 2. format names of cols
  x_var <- format_assay_terms(x_var)
  y_var <- format_assay_terms(y_var)
  if (is_display_log_independent) {
    x_var <- bquote(log[10]~.(x_var))

  } else {
    x_var <- x_var
  }
  if(is_display_log_response) {
    y_var <- bquote(log[10]~.(y_var))

  } else {
    y_var <- y_var
  }


  ## 2. Plot: data + fitted curves
  p_data_fit <- ggplot(dat[ , c("x","y")], aes(x = x, y = y)) +
    geom_point(alpha = 0.7) +
    geom_line(data = pred_df, aes(x = x, y = yhat, color = model, group = model)) +
    labs(title = "Observed data with fitted curves",
         x = x_var,
         y = y_var,
         color = "Model") +
    theme_bw()

  ## 3. Residual vs fitted
  p_resid <- ggplot(resid_df,
                    aes(x = fitted, y = residuals, color = model)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.6) +
    labs(title = "Residuals vs Fitted",
         x = "Fitted values",
         y = "Residuals",
         color = "Model") +
    theme_bw()


  # ── Drop any rows with NA parameter to prevent phantom facet panels ──
  fit_params_aic <- fit_params_aic[!is.na(fit_params_aic$parameter), , drop = FALSE]

  p_ci_params <- ggplot(fit_params_aic) +
    geom_point(aes(as.factor(model), estimate), color = "black") +
    facet_wrap(~ parameter, scales = 'free_x', ncol = 6) +
    geom_linerange(aes(as.factor(model), ymin = conf.low, ymax = conf.high)) +
    coord_flip() +
    scale_y_continuous(
      breaks = function(x) pretty(x, n = 3),
      expand = expansion(mult = c(0.1, 0.2))
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'top') +
    xlab('Model') +
    ylab('Parameter Estimate')


  if (!use_patchwork) {
    # Return a list of plots so caller can arrange as desired
    return(list(
      data_fit = p_data_fit,
      resid    = p_resid,
      p_ci_params = p_ci_params
    ))
  }

  # Use patchwork for a single panel display
  # Arrange: top row = data+fits, info; bottom row = resid, qq
  combined <- (p_data_fit) /
    (p_resid)/
    (p_ci_params) +
    plot_annotation(title = paste("Comparision of Model Fits for", unique(dat$antigen), "on", unique(dat$plate_nom)), tag_levels = "A",
                    theme = theme(
                      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                      plot.tag   = element_text(size = 12, face = "bold")
                    ))

  combined
}



#' Format common assay variable names for display
#'
#' Normalises and maps assay-related variable name strings (e.g. column names
#' such as \code{"mfi"} or \code{"concentration"}) to their conventional
#' display forms. Matching is case-insensitive. Unrecognized values are
#' returned unchanged.
#'
#' @param x Character vector of variable names to format.
#'
#' @return Character vector the same length as \code{x} with display-ready
#'   names. Known mappings are: \code{"mfi"} → \code{"MFI"},
#'   \code{"absorbance"} → \code{"Absorbance"},
#'   \code{"fluorescence"} → \code{"Fluorescence"},
#'   \code{"od"} → \code{"OD"},
#'   \code{"concentration"} → \code{"Concentration"}.
#' @export
format_assay_terms <- function(x) {
  lookup <- c(
    MFI = "MFI",
    Absorbance = "Absorbance",
    Fluorescence = "Fluorescence",
    OD = "OD",
    Concentration = "Concentration"
  )

  # normalize lookup keys for matching
  names(lookup) <- tolower(names(lookup))

  sapply(x, function(v) {
    key <- trimws(tolower(v))
    if (key %in% names(lookup)) lookup[[key]] else v
  }, USE.NAMES = FALSE)
}
