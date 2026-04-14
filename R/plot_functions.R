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
#' @param curve_id_lookup the order of the elements of the curve_id.
#' @param model_names Character vector of model names to process. Must match
#'   keys in \code{models_fit_list}. Default is
#'   \code{c("logistic5","loglogistic5","logistic4","loglogistic4","gompertz4")}.
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
                          curve_id_lookup,
                          model_names = c("logistic5","loglogistic5","logistic4","loglogistic4","gompertz4"),
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
      if (mname == "logistic5") {
        do.call(d2xlogistic5, c(list(x = x_new), fit_obj_coef))
      } else if (mname == "logistic4") {
        do.call(d2xlogistic4, c(list(x = x_new), fit_obj_coef))
      } else if (mname == "loglogistic5") {
        do.call(d2xloglogistic5, c(list(x = x_new), fit_obj_coef))
      } else if (mname == "loglogistic4") {
        do.call(d2xloglogistic4, c(list(x = x_new), fit_obj_coef))
      } else if (mname == "gompertz4") {
        do.call(d2xgompertz4, c(list(x = x_new), fit_obj_coef))
      }
    }, error = function(e) rep (NA_real_, length(x_new)))

    dydx <- tryCatch({
      if (mname == "logistic5") {
        do.call(dydxlogistic5, c(list(x = x_new), fit_obj_coef))
      } else if (mname == "logistic4") {
        do.call(dydxlogistic4, c(list(x = x_new), fit_obj_coef))
      } else if (mname == "loglogistic5") {
        do.call(dydxloglogistic5, c(list(x = x_new), fit_obj_coef))
      } else if (mname == "loglogistic4") {
        do.call(dydxloglogistic4, c(list(x = x_new), fit_obj_coef))
      } else if (mname == "gompertz4") {
        do.call(dydxgompertz4, c(list(x = x_new), fit_obj_coef))
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


  # parse out
  dat <- merge(dat, curve_id_lookup, by = "curve_id", all.x = TRUE)

  # dat <- parse_curve_id(dat, order = curve_id_element_order, keep = c("antigen", "feature", "plate", "nominal_sample_dilution"))

  if (all(c("plate", "nominal_sample_dilution") %in% names(dat))) {

    dat$plate_nom <- ifelse(
      is.na(dat$nominal_sample_dilution) | dat$nominal_sample_dilution == "",
      dat$plate,
      paste(dat$plate, dat$nominal_sample_dilution, sep = "-")
    )

  }

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
#'   Default is \code{c("logistic5","loglogistic5","logistic4","loglogistic4","gompertz4")}.
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
                                   model_names = c("logistic5","loglogistic5","logistic4","loglogistic4","gompertz4"),
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

  if (all(c("plate", "nominal_sample_dilution") %in% names(df))) {

    df$plate_nom <- ifelse(
      is.na(df$nominal_sample_dilution) | df$nominal_sample_dilution == "",
      df$plate,
      paste(df$plate, df$nominal_sample_dilution, sep = "-")
    )

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


#' Plot Standard Curve (Minimal + Uncertainty)
#'
#' Generates an interactive Plotly visualization of a fitted standard curve,
#' including standards, samples, fitted curve, confidence intervals, and
#' measurement uncertainty (pCoV).
#'
#' @param best_fit A list containing model outputs including best_data,
#' best_pred, best_glance, and optional best_curve_ci.
#' @param is_display_log_response Logical; whether to display response on log10 scale.
#' @param is_display_log_independent Logical; whether to display independent variable on log10 scale.
#' @param pcov_threshold Numeric threshold for percent coefficient of variation.
#' @param study_params study parameters including is_log_response and is_log_independent
#' @param response_variable Character; name of response variable column (default "mfi").
#' @param independent_variable Character; name of independent variable column (default "concentration").
#'
#' @return A plotly object.
#' @export
plot_standard_curve <- function(best_fit,
                                is_display_log_response,
                                is_display_log_independent,
                                pcov_threshold,
                                study_params,
                                # curve_id_element_order,
                                # curve_col = "curve_id",
                                response_variable = "mfi",
                                independent_variable = "concentration"
                                ) {
  p <- plotly::plot_ly()
  # best_fit_v <<- best_fit
  # mcmc_samples_in <<- mcmc_samples
  # mcmc_pred_in <<- mcmc_pred

  ## expand antigen and plate columns in the dataset
  #best_fit$best_data <-  parse_curve_id(best_fit$best_data, curve_col = curve_col, order = curve_id_element_order, keep = c("plate", "antigen"))

  # ── Resolve response column ────────────────────────────────────────
  resolved <- ensure_response_column(
    df           = best_fit$best_data,
    response_var = response_variable,
    coerce_numeric = TRUE,
    context      = "plot_standard_curve/best_data"
  )
  best_fit$best_data <- resolved$df
  response_variable  <- resolved$response_var

  if (!resolved$ok) {
    return(
      plotly::plot_ly() %>%
        plotly::layout(
          title = "Cannot plot: response variable not found",
          annotations = list(
            text = paste0(
              "Column '", response_variable,
              "' not found or has no finite values in standard data.<br>",
              "Available columns: ",
              paste(names(best_fit$best_data), collapse = ", ")
            ),
            xref = "paper", yref = "paper",
            x = 0.5, y = 0.5, showarrow = FALSE
          )
        )
    )
  }

  # ── Resolve independent variable ───────────────────────────────────
  if (!independent_variable %in% names(best_fit$best_data)) {
    if ("concentration" %in% names(best_fit$best_data)) {
      independent_variable <- "concentration"
    } else {
      return(plotly::plot_ly() %>%
               plotly::layout(title = "Missing independent variable column"))
    }
  }

  # ── Ensure stype exists ────────────────────────────────────────────
  if (!"stype" %in% names(best_fit$best_data)) {
    best_fit$best_data$stype <- "S"
  }

  # ── Resolve response column in sample_se too ───────────────────────
  samples_predicted_conc <- best_fit$sample_se
  if (!is.null(samples_predicted_conc) && nrow(samples_predicted_conc) > 0) {
    samp_resolved <- ensure_response_column(
      df           = samples_predicted_conc,
      response_var = response_variable,
      coerce_numeric = TRUE,
      context      = "plot_standard_curve/sample_se"
    )
    samples_predicted_conc <- samp_resolved$df
    if (samp_resolved$ok && samp_resolved$response_var != response_variable) {
      samples_predicted_conc[[response_variable]] <-
        samples_predicted_conc[[samp_resolved$response_var]]
    }
    samples_predicted_conc <- samples_predicted_conc[
      !is.nan(samples_predicted_conc$raw_predicted_concentration) &
        is.finite(samples_predicted_conc$raw_predicted_concentration), ,
      drop = FALSE
    ]
  } else {
    samples_predicted_conc <- data.frame(
      raw_predicted_concentration = numeric(0),
      pcov = numeric(0),
      stringsAsFactors = FALSE
    )
    samples_predicted_conc[[response_variable]] <- numeric(0)
  }

  best_fit$best_pred$pcov_threshold <- pcov_threshold

  safe_glance <- function(field, default = NA_real_) {
    val <- best_fit$best_glance[[field]]
    if (is.null(val) || length(val) == 0) return(default)
    val <- unlist(val)
    if (all(is.na(val))) return(default)
    val[1]
  }

  ### 1. RESPONSE VARIABLE (Y) log transform
  log_response_status <- isTRUE(as.logical(study_params$is_log_response))
  if (log_response_status && !isTRUE(is_display_log_response)) {
    best_fit$best_data[[response_variable]] <- 10^best_fit$best_data[[response_variable]]
    best_fit$best_pred$yhat               <- 10^best_fit$best_pred$yhat
    best_fit$best_glance$llod             <- 10^safe_glance("llod")
    best_fit$best_glance$ulod             <- 10^safe_glance("ulod")
    best_fit$best_glance$lloq_y           <- 10^safe_glance("lloq_y")
    best_fit$best_glance$uloq_y           <- 10^safe_glance("uloq_y")
    best_fit$best_glance$inflect_y        <- 10^safe_glance("inflect_y")
    best_fit$best_d2xy$d2x_y              <- 10^best_fit$best_d2xy$d2x_y
    if (!is.null(best_fit$best_curve_ci)) {
      best_fit$best_curve_ci$ci_lo        <- 10^best_fit$best_curve_ci$ci_lo
      best_fit$best_curve_ci$ci_hi        <- 10^best_fit$best_curve_ci$ci_hi
    }
    if (nrow(samples_predicted_conc) > 0 &&
        response_variable %in% names(samples_predicted_conc)) {
      samples_predicted_conc[[response_variable]] <-
        10^samples_predicted_conc[[response_variable]]
    }
  }

  ### 2. INDEPENDENT VARIABLE (X) log transform
  log_x_status <- isTRUE(as.logical(study_params$is_log_independent))
  if (log_x_status && !isTRUE(is_display_log_independent)) {
    best_fit$best_data$concentration       <- 10^best_fit$best_data$concentration
    best_fit$best_pred$x                   <- 10^best_fit$best_pred$x
    best_fit$best_glance$lloq              <- 10^safe_glance("lloq")
    best_fit$best_glance$uloq              <- 10^safe_glance("uloq")
    best_fit$best_glance$inflect_x         <- 10^safe_glance("inflect_x")
    best_fit$best_d2xy$x                   <- 10^best_fit$best_d2xy$x
    if (nrow(samples_predicted_conc) > 0) {
      samples_predicted_conc$raw_predicted_concentration <-
        10^samples_predicted_conc$raw_predicted_concentration
    }
    best_fit$best_glance$mindc             <- 10^safe_glance("mindc")
    best_fit$best_glance$maxdc             <- 10^safe_glance("maxdc")
    best_fit$best_glance$minrdl            <- 10^safe_glance("minrdl")
    best_fit$best_glance$maxrdl            <- 10^safe_glance("maxrdl")
    if (!is.null(best_fit$best_curve_ci)) {
      best_fit$best_curve_ci$x             <- 10^best_fit$best_curve_ci$x
    }
  }

  y3_label <- "Precision Coefficient of Variation (pCoV %)"

  if (is_display_log_response) {
    response_formatted <- format_assay_terms(response_variable)
    cat("FORMATTED:", response_formatted, "\n")
    y_label <- paste0("log<sub>10</sub> ", response_formatted)
  } else {
    y_label <- format_assay_terms(response_variable)
  }
  # y_label_v <<- y_label
  #
  if (is_display_log_independent) {
    x_label <- paste0("log<sub>10</sub> ", format_assay_terms(independent_variable))
  } else {
    x_label <- format_assay_terms(independent_variable)
  }

  ### 3. MODEL NAME
  model_name <- best_fit$best_model_name
  title_model_name <- switch(
    model_name,
    "logistic4" = "4-parameter Logistic",
    "loglogistic4" = "4-parameter Log-Logistic",
    "gompertz4" = "4-parameter Gompertz type",
    "logistic5" = "5-parameter Logistic",
    "loglogistic5" = "5-parameter Log-Logistic",
    model_name
  )

  ## 3b. PREPARE SAMPLE-UNCERTAINTY (single scaling)
  print(names(best_fit$best_pred))
  se_model   <- best_fit$best_pred$pcov
  se_samples <- samples_predicted_conc$pcov
  se_all     <- c(best_fit$best_pred$pcov, samples_predicted_conc$pcov)
  se_range   <- range(se_all, na.rm = TRUE)

  se_max <- 125
  se_min <- -2
  se_axis_limits <- c(se_min, se_max * 1.1)
  dtick <- ifelse(se_max > 19, ifelse(se_max > 35, 10, 5), 1)

  ### 4. RAW POINTS
  plot_std <- best_fit$best_data

  # if (is_display_log_independent) {
  #   glance_fda_lloq_conc <- log10(best_fit$best_glance$lloq_fda2018_concentration)
  #   glance_fda_2018_uloq_conc <- log10(best_fit$best_glance$uloq_fda2018_concentration)
  # } else {
  #   glance_fda_lloq_conc <- best_fit$best_glance$lloq_fda2018_concentration
  #   glance_fda_2018_uloq_conc <- best_fit$best_glance$uloq_fda2018_concentration
  # }
  #
  # plot_std$fda2018_class <- ifelse(
  #   plot_std[[independent_variable]] >= glance_fda_lloq_conc &
  #     plot_std[[independent_variable]] <= glance_fda_2018_uloq_conc,
  #   "Standards (+ FDA 2018)",
  #   "Standards (- FDA 2018)"
  # )
  std_in <- plot_std[plot_std$stype == "S",]
  #std_in <- plot_std[
  #   plot_std$stype == "S" &
  #     plot_std$fda2018_class == "Standards (+ FDA 2018)", ]
  # std_out <- plot_std[
  #   plot_std$stype == "S" &
  #     plot_std$fda2018_class == "Standards (- FDA 2018)", ]
   blanks <- plot_std[plot_std$stype == "B", ]
  #
  ### Standards inside FDA range (circle)
  p <- p %>% plotly::add_trace(
    data = std_in,
    x = std_in[[independent_variable]],
    y = std_in[[response_variable]],
    type = "scatter",
    mode = "markers",
    name = "Standards",
    legendgroup = "standards",
    marker = list(color = "#2b3d26", symbol = "circle"),
    text = ~paste0(
      "<br>", format_assay_terms(independent_variable), ": ",
      std_in[[independent_variable]],
      "<br>Dilution Factor: ", dilution,
      "<br>", format_assay_terms(response_variable), ": ",
      std_in[[response_variable]]
     #, "<br> FDA 2018 Status: ", gsub("Standards ", "", std_in$fda2018_class)
    ),
    hoverinfo = "text"
  )

  # ### Standards outside FDA range (triangle)
  # p <- p %>% plotly::add_trace(
  #   data = std_out,
  #   x = std_out[[independent_variable]],
  #   y = std_out[[response_variable]],
  #   type = "scatter",
  #   mode = "markers",
  #   name = "Standards (- FDA 2018)",
  #   legendgroup = "standards",
  #   marker = list(color = "#2b3d26", symbol = "triangle-up", size = 8),
  #   text = ~paste0(
  #     "<br>", format_assay_terms(independent_variable), ": ",
  #     std_out[[independent_variable]],
  #     "<br>Dilution Factor: ", dilution,
  #     "<br>", format_assay_terms(response_variable), ": ",
  #     std_out[[response_variable]],
  #     "<br>FDA 2018 Status: ", gsub("Standards ", "", std_out$fda2018_class)
  #   ),
  #   hoverinfo = "text"
  # )
  #
  ### Blanks
  p <- p %>% plotly::add_trace(
    data = blanks,
    x = blanks[[independent_variable]],
    y = blanks[[response_variable]],
    type = "scatter",
    mode = "markers",
    name = "Geometric Mean of Blanks",
    marker = list(color = "#c2b280", symbol = "circle"),
    text = ~paste0(
      "<br>", format_assay_terms(independent_variable), ": ",
      blanks[[independent_variable]],
      "<br>Dilution Factor: ", dilution,
      "<br>", format_assay_terms(response_variable), ": ",
      blanks[[response_variable]]
    ),
    hoverinfo = "text"
  )

  ### 5. FITTED CURVE
  p <- p %>% add_lines(
    x = best_fit$best_pred$x,
    y = best_fit$best_pred$yhat,
    name = "Fitted Curve",
    legendgroup = "fitted_curve",
    showlegend = TRUE,
    line = list(color = "#2b3d26")
  )

  ### 5b. 95% CI BANDS (delta method)
  if (!is.null(best_fit$best_curve_ci)) {
    p <- p %>% add_lines(
      x           = best_fit$best_curve_ci$x,
      y           = best_fit$best_curve_ci$ci_lo,
      name        = "95% CI",
      line        = list(color = "#2b3d26", dash = "dash"),
      legendgroup = "fitted_curve"
    ) %>% add_lines(
      x           = best_fit$best_curve_ci$x,
      y           = best_fit$best_curve_ci$ci_hi,
      name        = "",
      line        = list(color = "#2b3d26", dash = "dash"),
      legendgroup = "fitted_curve",
      showlegend  = FALSE
    )
  }


  ### 6. LOD lines (horizontal)
  p <- p %>% add_lines(
    x = best_fit$best_pred$x,
    y = best_fit$best_fit_summary$ulod,
    name = paste("Upper LOD: (",
                 round(best_fit$best_fit_summary$maxdc, 3), ",",
                 round(best_fit$best_fit_summary$ulod, 3), ")"),
    line = list(color = "#e25822", dash = "dash"),
    legendgroup = "linked_ulod",
    visible = "legendonly"
  )

  p <- p %>% add_lines(
    x = best_fit$best_pred$x,
    y = best_fit$best_fit_summary$llod,
    name = paste("Lower LOD: (",
                 round(best_fit$best_fit_summary$mindc, 3), ",",
                 round(best_fit$best_fit_summary$llod, 3), ")"),
    line = list(color = "#e25822", dash = "dash"),
    legendgroup = "linked_llod",
    visible = "legendonly"
  )

  ###
  y_min <- min(best_fit$best_data[[response_variable]], na.rm = TRUE)
  y_max <- max(best_fit$best_data[[response_variable]], na.rm = TRUE)

  ### 6b. MDC / RDL vertical lines
  if (!is.na(best_fit$best_fit_summary$mindc)) {
    p <- p %>% add_lines(
      x = c(best_fit$best_fit_summary$mindc, best_fit$best_fit_summary$mindc),
      y = c(y_min, y_max),
      name = paste("Lower DC:", round(best_fit$best_fit_summary$mindc, 3)),
      line = list(color = "#e25822", dash = "dash"),
      legendgroup = "linked_llod",
      showlegend = FALSE, hoverinfo = "text", visible = "legendonly"
    )
  }

  if (!is.na(best_fit$best_fit_summary$minrdl)) {
    p <- p %>% add_lines(
      x = c(best_fit$best_fit_summary$minrdl, best_fit$best_fit_summary$minrdl),
      y = c(y_min, y_max),
      name = paste("Lower RDL:", round(best_fit$best_fit_summary$minrdl, 3)),
      line = list(color = "#e25822"),
      legendgroup = "linked_llod",
      showlegend = TRUE, hoverinfo = "text", visible = "legendonly"
    )
  }

  if (!is.na(best_fit$best_fit_summary$maxdc)) {
    p <- p %>% add_lines(
      x = c(best_fit$best_fit_summary$maxdc, best_fit$best_fit_summary$maxdc),
      y = c(y_min, y_max),
      name = paste("Upper DC:", round(best_fit$best_fit_summary$maxdc, 3)),
      line = list(color = "#e25822", dash = "dash"),
      legendgroup = "linked_ulod",
      showlegend = FALSE, hoverinfo = "text", visible = "legendonly"
    )
  }

  if (!is.na(best_fit$best_fit_summary$maxrdl)) {
    p <- p %>% add_lines(
      x = c(best_fit$best_fit_summary$maxrdl, best_fit$best_fit_summary$maxrdl),
      y = c(y_min, y_max),
      name = paste("Upper RDL:", round(best_fit$best_fit_summary$maxrdl, 3)),
      line = list(color = "#e25822"),
      legendgroup = "linked_ulod",
      showlegend = TRUE, hoverinfo = "text", visible = "legendonly"
    )
  }

  ## LLOQ vertical line
  p <- p %>% add_lines(
    x = c(best_fit$best_fit_summary$lloq),
    y = c(y_min, y_max),
    name = paste("Lower LOQ: (",
                 round(best_fit$best_fit_summary$lloq, 3), ",",
                 round(best_fit$best_fit_summary$lloq_y, 3), ")"),
    line = list(color = "#875692"),
    legendgroup = "linked_lloq",
    hoverinfo = "text", visible = "legendonly"
  )

  ### ULOQ vertical line
  p <- p %>% add_lines(
    x = c(best_fit$best_fit_summary$uloq),
    y = c(y_min, y_max),
    name = paste("Upper LOQ: (",
                 round(best_fit$best_fit_summary$uloq, 3), ",",
                 round(best_fit$best_fit_summary$uloq_y, 3), ")"),
    line = list(color = "#875692"),
    legendgroup = "linked_uloq",
    hoverinfo = "text", visible = "legendonly"
  )

  ### Horizontal LOQ lines
  p <- p %>% add_lines(
    x = best_fit$best_pred$x,
    y = best_fit$best_fit_summary$uloq_y,
    name = "",
    legendgroup = "linked_uloq", showlegend = FALSE,
    line = list(color = "#875692"), visible = "legendonly"
  )

  p <- p %>% add_lines(
    x = best_fit$best_pred$x,
    y = best_fit$best_fit_summary$lloq_y,
    name = "",
    legendgroup = "linked_lloq", showlegend = FALSE,
    line = list(color = "#875692"), visible = "legendonly"
  )

  # ### 8a. SECOND DERIVATIVE (y2 axis)
  if ("best_d2xy" %in% names(best_fit)) {
      p <- p %>% add_lines(
        x = best_fit$best_d2xy$x,
        y = best_fit$best_d2xy$d2x_y,
        name = "2nd Derivative of x given y",
        yaxis = "y2",
        line = list(color = "#604e97"),
        visible = "legendonly"
      )
  }

  ## 9. Samples - interpolated
  p <- p %>% add_trace(
    data = samples_predicted_conc,
    x = ~raw_predicted_concentration,
    y = samples_predicted_conc[[response_variable]],
    type = "scatter",
    mode = "markers",
    name = "Samples",
    marker = list(color = "#d1992a", symbol = "circle"),
    text = ~paste("Predicted", x_label, ":", raw_predicted_concentration,
                  "<br>", y_label, ":", samples_predicted_conc[[response_variable]],
                  "<br>Patient ID:", patientid,
                  "<br> Timepoint:", timeperiod,
                  "<br>Well:", well
                  # ,"<br>LOQ Gate Class:", samples_predicted_conc$gate_class_loq,
                  # "<br>LOD Gate Class:", samples_predicted_conc$gate_class_lod,
                  # "<br> PCOV Gate Class:", samples_predicted_conc$gate_class_pcov
                  ),
    hovertemplate = "%{text}<extra></extra>"
  )
  ### 8b. Sample uncertainty (y3 axis) — interpolated
  unc_col <- list(color = "#e68fac")
  p <- p %>% add_lines(
    x = best_fit$best_pred$x,
    y = best_fit$best_pred$pcov,
    name = "Measurement Uncertainty",
    yaxis = "y3",
    line = unc_col,
    legendgroup = "linked_interp_uncertainty",
    visible = "legendonly"
  ) %>% add_trace(
    data = samples_predicted_conc,
    x = ~raw_predicted_concentration,
    y = ~pcov,
    type = "scatter",
    mode = "markers",
    name = "",
    marker = list(color = "#800032", symbol = "circle"),
    text = ~paste("Predicted", x_label, ":", raw_predicted_concentration,
                  "<br>Coefficient of Variation (pCoV):", round(pcov, 2), "%"),
    yaxis = "y3",
    legendgroup = "linked_interp_uncertainty",
    showlegend = FALSE,
    hovertemplate = "%{text}<extra></extra>",
    visible = "legendonly"
  ) %>% add_lines(
    x = best_fit$best_pred$x,
    y = best_fit$best_pred$pcov_threshold,
    name = paste0("pCoV Threshold: ", best_fit$best_pred$pcov_threshold, "%"),
    yaxis = "y3",
    line = list(color = "#e68fac", dash = "dash"),
    legendgroup = "linked_interp_uncertainty",
    visible = "legendonly"
  )

  # MCMC ROBUST SAMPLES


  #   ### 9c. MCMC pCoV scatter points at sample locations (y3 axis)
  #   if ("pcov_robust_concentration" %in% names(mcmc_samples)) {
  #     pcov_valid <- is.finite(mcmc_samples$pcov_robust_concentration) &
  #       is.finite(mcmc_x)
  #
  #     if (any(pcov_valid)) {
  #       p <- p %>% plotly::add_trace(
  #         x = mcmc_x[pcov_valid],
  #         y = mcmc_samples$pcov_robust_concentration[pcov_valid],
  #         type = "scatter",
  #         mode = "markers",
  #         name = "",
  #         marker = list(color = "#800032", symbol = "diamond", size = 5),
  #         text = paste0(
  #           "MCMC ", x_label, ": ", round(mcmc_x[pcov_valid], 4),
  #           "<br>MCMC pCoV: ",
  #           formatC(mcmc_samples$pcov_robust_concentration[pcov_valid], format = "g", digits = 4), "%"
  #         ),
  #         yaxis = "y3",
  #         legendgroup = "linked_mcmc_uncertainty",
  #         showlegend = FALSE,
  #         hovertemplate = "%{text}<extra></extra>",
  #         visible = "legendonly"
  #       )
  #     }
  #   }
  # }



  # ### 10. INFLECTION POINT
  if ("best_fit_summary" %in% names(best_fit)) {

    p <- p %>% add_trace(
      x = best_fit$best_fit_summary$inflect_x,
      y = best_fit$best_fit_summary$inflect_y,
      type = "scatter",
      mode = "markers",
      name = paste("Inflection Point: (",
                   round(best_fit$best_fit_summary$inflect_x, 3), ",",
                   round(best_fit$best_fit_summary$inflect_y, 3), ")"),
      legendgroup = "fitted_curve",
      showlegend = TRUE,
      marker = list(color = "#2724F0", size = 8)
    )

  }

  ### 11. LAYOUT
  p <- p %>% layout(
    title = paste(
      "Fitted", title_model_name, "Model (",
      unique(best_fit$best_fit_summary$plate), ",",
      unique(best_fit$best_fit_summary$antigen), ")"
    ),
    xaxis = list(
      title    = x_label,
      showgrid = TRUE,
      zeroline = FALSE
    ),
    yaxis = list(
      title    = y_label,
      showgrid = TRUE,
      zeroline = TRUE
    ),
    legend = list(
      x       = 1.1,
      y       = 1,
      xanchor = "left"
    ),
    font = list(size = 12),
    yaxis2 = list(
      showticklabels = FALSE,
      title          = "",
      tickmode       = "linear",
      dtick          = 10,
      overlaying     = "y",
      side           = "right",
      showgrid       = FALSE,
      zeroline       = FALSE
    ),
    yaxis3 = list(
      overlaying     = "y",
      side           = "right",
      title          = y3_label,
      range          = se_axis_limits,
      tickmode       = "linear",
      type           = "linear",
      dtick          = dtick,
      showgrid       = FALSE,
      zeroline       = FALSE,
      showticklabels = TRUE
    )
  )

  return(p)
}
