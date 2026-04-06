#' Compute Inflection Point of Fitted Curve
#'
#' Calculates the inflection point (x and y coordinates) for a fitted nonlinear model.
#' The x-coordinate is computed analytically based on the model type, and the y-coordinate
#' is obtained by evaluating the fitted model at that x-value (via \code{predict()} when possible).
#'
#' Supported model types include 4- and 5-parameter logistic and log-logistic models,
#' as well as the 4-parameter Gompertz model.
#'
#' @param model_name Character string specifying the model type.
#'   Supported values include:
#'   \itemize{
#'     \item \code{"logistic5"}
#'     \item \code{"logistic4"}
#'     \item \code{"loglogistic5"}
#'     \item \code{"loglogistic4"}
#'     \item \code{"gompertz4"}
#'   }
#' @param fit A fitted model object with coefficients accessible via \code{coef()}
#'   and compatible with \code{predict()}.
#' @param fixed_a_result Optional numeric value for parameter \code{a}. If provided,
#'   this overrides the estimated \code{a} from the model coefficients.
#' @param independent_variable Character string giving the name of the predictor
#'   variable used in the model. This is required to construct \code{newdata} for prediction.
#' @param verbose Logical; if \code{TRUE}, prints a message including the computed
#'   inflection point coordinates.
#'
#' @details
#' The inflection point is defined as the point on the curve where the second derivative equals zero.
#' For supported models, the x-coordinate is computed analytically:
#' \itemize{
#'   \item Logistic 5-parameter: \eqn{x = c - b \log(g)}
#'   \item Logistic 4-parameter: \eqn{x = c}
#'   \item Log-logistic 5-parameter: \eqn{x = c + \log(g)/b}
#'   \item Log-logistic 4-parameter: \eqn{x = c}
#'   \item Gompertz 4-parameter: \eqn{x = c}
#' }
#'
#' The y-coordinate is computed by evaluating the fitted model at the inflection x-value.
#' If \code{predict()} fails, an analytical expression is used as a fallback.
#'
#' When \code{verbose = TRUE}, the function prints the inflection point as a coordinate:
#' \preformatted{
#' [compute_inflection_point] Inflection point: (x = ..., y = ...)
#' }
#'
#' @return A named list with:
#' \itemize{
#'   \item \code{inflect_x} Numeric x-coordinate of the inflection point
#'   \item \code{inflect_y} Numeric y-coordinate of the inflection point
#' }
#'
#' @examples
#' \dontrun{
#' result <- compute_inflection_point(
#'   model_name = "logistic5",
#'   fit = fitted_model,
#'   fixed_a_result = NULL,
#'   independent_variable = "log_dilution"
#' )
#'
#' result$inflect_x
#' result$inflect_y
#' }
#'
#' @export
compute_inflection_point <- function(model_name, fit, fixed_a_result, independent_variable,  verbose = TRUE) {
  params <- coef(fit)
  g <- as.numeric(if ("g" %in% names(params)) params["g"] else 1)# auto default
  a <-  as.numeric(ifelse(!is.null(fixed_a_result), fixed_a_result, params["a"]))
  b <- as.numeric(params["b"])
  c <- as.numeric(params["c"])
  d <- as.numeric(params["d"])

  # Calculate x-coordinate of inflection point analytically
  # The inflection point is where the second derivative equals zero
  inflect_x  <- tryCatch({
    if (model_name == "logistic5") {
      # For 5 parameter model: x_inflect = c - b*ln(g)
      c - b * log(g)
    } else if (model_name == "logistic4") {
      # For 4 parameter logistic: x_inflect = c
      c
    } else if (model_name == "loglogistic5") {
      # For 5 parameter decreasing logistic: x_inflect = c + ln(g)/b
      c + (log(g) / b)
    } else if (model_name == "loglogistic4") {
      # For 4 parameter log-logistic: x_inflect = c
      c
    } else if (model_name == "gompertz4") {
      # For Gompertz: x_inflect = c (where second derivative = 0)
      c
    }
  }, error = function(e) NA)

  inflect_x <- as.numeric(inflect_x)

  # Calculate y-coordinate by evaluating the fitted model at inflect_x
  # This ensures the inflection point lies exactly on the fitted curve
  inflect_y <- tryCatch({
    # Use predict() to evaluate the fitted model at the inflection point
    # This is more robust than analytical formulas as it uses the actual fitted model
    newdata <- setNames(data.frame(x = inflect_x), independent_variable)
    predicted_y <- predict(fit, newdata = newdata)
    as.numeric(predicted_y)
  }, error = function(e) {
    # Fallback to analytical calculation if predict fails
    if (verbose) message("predict() failed, using analytical formula for inflect_y")
    tryCatch({
      switch(model_name,
             "logistic5" = d + (a - d) / (1 + exp((inflect_x - c) / b))^g,
             "logistic4" = d + (a - d) / (1 + exp((inflect_x - c) / b)),
             "loglogistic5" = a + (d - a) * (1 + g * exp(-b * (inflect_x - c)))^(-1/g),
             "loglogistic4" = a + (d - a) / (1 + exp(-b * (inflect_x - c))),
             "gompertz4" = a + (d - a) * exp(-exp(-b * (inflect_x - c))),
             NA
      )
    }, error = function(e2) NA)
  })

  inflect_y <- as.numeric(inflect_y)

  if (verbose) {
    message(sprintf(
      "[compute_inflection_point] Inflection point: (x = %.6f, y = %.6f)",
      inflect_x, inflect_y
    ))
  }

  return(list(inflect_x = inflect_x, inflect_y = inflect_y))
}




#' Summarize Fitted Model with QC Metrics
#'
#' Generates a one-row summary of a fitted model, including parameter estimates,
#' goodness-of-fit statistics, inflection point coordinates, and metadata.
#' The result is stored in \code{best_fit$fit_summary}.
#'
#' This function is designed to operate on the output of a model fitting step
#' and produce a consistent summary structure for downstream analysis,
#' QC filtering, or reporting.
#'
#' @param best_fit A list containing model results, including:
#'   \itemize{
#'     \item \code{best_fit$best_fit}: Fitted model object
#'     \item \code{best_fit$best_data}: Data used for fitting
#'     \item \code{best_fit$best_model_name}: Character model identifier
#'   }
#' @param response_variable Character string naming the response variable in \code{best_data}.
#' @param independent_variable Character string naming the predictor variable used in the model.
#' @param fixed_a_result Optional numeric value for parameter \code{a}. Overrides the fitted value if provided.
#' @param antigen_settings List of antigen-specific settings (reserved for future use).
#' @param antigen_fit_options List of model fitting options, including:
#'   \itemize{
#'     \item \code{blank_option}
#'     \item \code{is_log_response}
#'     \item \code{is_log_concentration}
#'     \item \code{apply_prozone}
#'   }
#' @param verbose Logical; if \code{TRUE}, prints progress and diagnostic messages.
#'
#' @details
#' The returned summary includes:
#' \itemize{
#'   \item Model parameters: \code{a, b, c, d, g}
#'   \item Inflection point coordinates: \code{inflect_x}, \code{inflect_y}
#'   \item Goodness-of-fit metrics: residual sum of squares (RSS), mean squared error (MSE),
#'         R-squared, Akaike information criterion (AIC), Bayesian information criterion (BIC),
#'         and log-likelihood
#'   \item Convergence diagnostics: iteration count and convergence status
#'   \item Metadata: source, transformation flags, and model formula
#' }
#'
#' The inflection point is computed using \code{compute_inflection_point()} and
#' corresponds to the location where the second derivative of the fitted curve equals zero.
#'
#' If no valid model or data is available, a placeholder summary with \code{NA} values
#' is returned via \code{.make_na_glance()}.
#'
#' @return A list identical to \code{best_fit} with an additional element:
#' \itemize{
#'   \item \code{fit_summary}: A one-row \code{data.frame} containing model summary and QC metrics
#' }
#'
#' @examples
#' \dontrun{
#' best_fit <- fit_model(data)
#'
#' best_fit <- summarize_fit(
#'   best_fit = best_fit,
#'   response_variable = "response",
#'   independent_variable = "log_dilution",
#'   fixed_a_result = NULL,
#'   antigen_settings = list(),
#'   antigen_fit_options = list(
#'     blank_option = "none",
#'     is_log_response = FALSE,
#'     is_log_concentration = TRUE,
#'     apply_prozone = FALSE
#'   )
#' )
#'
#' best_fit$best_fit_summary
#' }
#'
#' @seealso \code{\link{compute_inflection_point}}
#' @export
summarize_fit <- function(best_fit,
                          response_variable,
                          independent_variable,
                          fixed_a_result,
                          antigen_settings,
                          antigen_fit_options,
                          verbose = TRUE) {

  fit        <- best_fit$best_fit
  best_data  <- best_fit$best_data
  model_name <- best_fit$best_model_name

  # # ── Early exit if no valid fit ──
  # if (is.null(fit) || is.null(best_data) || nrow(best_data) == 0) {
  #   if (verbose) message("[summarize_fit] No valid fit or data -> NA fit summary")
  #   best_fit$best_glance <- .make_na_glance(best_data, model_name,
  #                                           fixed_a_result, antigen_fit_options)
  #   return(best_fit)
  # }

  # Extract model coefficients
  coefs <- tryCatch(coef(fit), error = function(e) c(a=NA, b=NA, c=NA, d=NA, g=NA))

  inflect <- tryCatch({
    compute_inflection_point(
      model_name = model_name,
      fit = fit,
      fixed_a_result = fixed_a_result,
      independent_variable = independent_variable,
      verbose = verbose
    )
  }, error = function(e) list(inflect_x = NA_real_, inflect_y = NA_real_))

  # Goodness-of-fit
  rss       <- tryCatch(sum(residuals(fit)^2), error = function(e) NA_real_)
  df_resid  <- tryCatch(df.residual(fit), error = function(e) NA_real_)
  nobs_fit  <- tryCatch(nobs(fit), error = function(e) NA_integer_)
  aic_val   <- tryCatch(AIC(fit), error = function(e) NA_real_)
  bic_val   <- tryCatch(BIC(fit), error = function(e) NA_real_)
  loglik_val <- tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_)
  converged   <- fit$convInfo$isConv
  iter        <- fit$convInfo$finIter
  mse_val   <- if (!is.na(rss) && !is.na(df_resid) && df_resid > 0) rss / df_resid else NA_real_

  # R-squared
  y_obs     <- tryCatch(best_data[[response_variable]], error = function(e) numeric(0))
  y_pred    <- tryCatch(predict(fit, newdata = best_data), error = function(e) numeric(0))
  rsq       <- if (length(y_obs) > 1 && length(y_pred) == length(y_obs)) {
    ss_res <- sum((y_obs - y_pred)^2, na.rm = TRUE)
    ss_tot <- sum((y_obs - mean(y_obs, na.rm = TRUE))^2, na.rm = TRUE)
    if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  } else NA_real_

  cv_val <- if (!is.na(mse_val) && length(y_obs) > 0) {
    sqrt(mse_val) / mean(y_obs, na.rm = TRUE) * 100
  } else NA_real_

  glance_df <- data.frame(
    curve_id                = safe_unique(best_data$curve_id),
    iter                    = iter,
    status                  = converged,
    crit                    = as.character(model_name),
    a  = as.numeric(coefs["a"] %||% fixed_a_result %||% NA_real_),
    b  = as.numeric(coefs["b"] %||% NA_real_),
    c  = as.numeric(coefs["c"] %||% NA_real_),
    d  = as.numeric(coefs["d"] %||% NA_real_),
    g  = as.numeric(coefs["g"] %||% NA_real_),
    inflect_x = inflect$inflect_x,
    inflect_y = inflect$inflect_y,
    stringsAsFactors = FALSE
  )

  # Add remaining fit stats
  glance_df$dfresidual  <- df_resid
  glance_df$nobs        <- nobs_fit
  glance_df$rsquare_fit <- rsq
  glance_df$aic         <- aic_val
  glance_df$bic         <- bic_val
  glance_df$loglik      <- loglik_val
  glance_df$mse         <- mse_val
  glance_df$cv          <- cv_val
  glance_df$source      <- safe_unique(best_data$source)
  glance_df$bkg_method  <- antigen_fit_options$blank_option %||% NA_character_
  glance_df$is_log_response <- antigen_fit_options$is_log_response %||% NA
  glance_df$is_log_x    <- antigen_fit_options$is_log_concentration %||% NA
  glance_df$apply_prozone <- antigen_fit_options$apply_prozone %||% NA
  glance_df$formula     <- tryCatch(sub("I\\((.*)\\)", "\\1", paste(deparse(formula(fit)), collapse = " ")), error = function(e) NA_character_)
  glance_df$last_concentration_calc_method <- "interpolated"


  best_fit$best_fit_summary <- glance_df

  return(best_fit)
}
