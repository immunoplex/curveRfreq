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

#' Compute lower and upper limits of detection (LOD)
#'
#' Calculates the lower limit of detection (LLOD) and upper limit of detection (ULOD)
#' based on model confidence intervals and optional blank standard error adjustment.
#'
#' @param best_fit A list containing model results with at least:
#' \describe{
#'   \item{best_ci}{A data.frame of parameter confidence intervals with columns
#'   \code{parameter}, \code{conf.low}, and \code{conf.high}.}
#'   \item{best_data}{A data.frame of the data used to fit the model.}
#' }
#' @param fixed_a_result Optional numeric value for the fixed lower asymptote ("a").
#' If provided, the LLOD is computed using this value plus a margin of error.
#' @param std_error_blank Optional numeric standard error of the blank. If \code{NULL}
#' or \code{NA}, it defaults to 0 when \code{fixed_a_result} is used.
#' @param verbose Logical; if \code{TRUE}, prints intermediate values used in the calculation.
#'
#' @return A named list with:
#' \describe{
#'   \item{llod}{Lower limit of detection.}
#'   \item{ulod}{Upper limit of detection. Returns \code{NA} if invalid (e.g., negative or less than LLOD).}
#' }
#'
#' @details
#' The ULOD is derived from the lower confidence bound of the upper asymptote ("d").
#' The LLOD is either:
#' \itemize{
#'   \item The upper confidence bound of the lower asymptote ("a"), or
#'   \item A fixed asymptote value plus a margin of error based on a t-distribution
#'   and the blank standard error.
#' }
#'
#' Invalid ULOD values (negative or less than LLOD) are set to \code{NA}.
#'
#' @examples
#' \dontrun{
#' lods <- generate_lods(best_fit, fixed_a_result = NULL, std_error_blank = NULL)
#' }
#'
#' @export
generate_lods <- function(best_fit, fixed_a_result, std_error_blank,  verbose = TRUE) {

  best_ci <- best_fit$best_ci
  best_data <- best_fit$best_data

  ulod <- best_ci[best_ci$parameter == "d",]$conf.low

  if (!is.null(fixed_a_result)) {
    if (is.null(std_error_blank) || is.na(std_error_blank)) {
      std_error_blank <- 0
    }
    critical_value <- qt(0.975, df = nrow(best_data) - length(best_ci$parameter))
    cat("critical value:\n")
    print(critical_value)
    cat("Blank SE:\n")
    print(std_error_blank)
    margin_of_error <- critical_value * std_error_blank
    llod <- fixed_a_result + margin_of_error
  } else {
    llod <- best_ci[best_ci$parameter == "a",]$conf.high
  }

  if (ulod < 0 || ulod < llod) {
    ulod  <- NA_real_
  }
  return(list(llod = llod, ulod = ulod))

}

#' Compute MDC and RDL values from fitted model and LODs
#'
#' Computes minimum and maximum detectable concentrations (MDC) and
#' reliable detection limits (RDL) by solving for concentrations where
#' the fitted curve or its confidence bounds intersect the LODs.
#'
#' @param best_fit A list containing:
#' \describe{
#'   \item{best_fit}{A fitted model object supporting \code{predict()}, \code{coef()}, and \code{vcov()}.}
#'   \item{best_data}{A data.frame used for fitting the model.}
#' }
#' @param lods A list containing \code{llod} and \code{ulod}, typically from
#' \code{\link{generate_lods}}.
#' @param independent_variable Character string naming the independent variable column
#' in \code{best_data}.
#' @param verbose Logical; if \code{TRUE}, prints computed MDC and RDL values.
#'
#' @return A named list with:
#' \describe{
#'   \item{mindc}{Minimum detectable concentration (fitted curve = LLOD).}
#'   \item{maxdc}{Maximum detectable concentration (fitted curve = ULOD).}
#'   \item{minrdl}{Minimum reliable detection limit (lower CI = LLOD).}
#'   \item{maxrdl}{Maximum reliable detection limit (upper CI = ULOD).}
#' }
#'
#' @details
#' Root-finding is performed using \code{uniroot()} over the observed range of the
#' independent variable.
#'
#' The prediction standard error is computed using the delta method:
#' \deqn{SE(\hat{y}) = \sqrt{\nabla f(\theta)^T V \nabla f(\theta)}}
#' where \eqn{V} is the variance-covariance matrix of the fitted parameters.
#'
#' Reliable detection limits (RDL) are defined using the confidence interval:
#' \itemize{
#'   \item Lower bound: \eqn{\hat{y}(x) - t_{\alpha/2} \cdot SE(x)}
#'   \item Upper bound: \eqn{\hat{y}(x) + t_{\alpha/2} \cdot SE(x)}
#' }
#'
#' If no valid root is found within the observed data range, the corresponding
#' value is returned as \code{NA}.
#'
#' @examples
#' \dontrun{
#' lods <- generate_lods(best_fit)
#' mdc_rdl <- generate_mdc_rdl(best_fit, lods, independent_variable = "log_dilution")
#' }
#'
#' @seealso \code{\link{generate_lods}}
#'
#' @export
generate_mdc_rdl <- function(best_fit, lods,
                             independent_variable, verbose = TRUE) {

  # Refactor 1: accept pre-computed lods instead of recomputing
  llod <- as.numeric(lods$llod)
  ulod <- as.numeric(lods$ulod)

  fit        <- best_fit$best_fit
  best_data  <- best_fit$best_data

  # x-range of the standards (search bounds for uniroot)
  x_data <- best_data[[independent_variable]]
  x_lo   <- min(x_data, na.rm = TRUE)
  x_hi   <- max(x_data, na.rm = TRUE)

  # Degrees of freedom for t-quantile (used by CI calculations)
  n_params <- length(coef(fit))
  n_obs    <- nrow(best_data)
  t_crit   <- qt(0.975, df = n_obs - n_params)

  # Variance-covariance matrix of fitted parameters
  V <- vcov(fit)

  # Refactor 2: hoist invariants out of pred_se inner loop
  theta <- coef(fit)
  p     <- length(theta)
  rhs   <- as.list(formula(fit))[[3]]

  # --- helper: predicted y at a single x ----------------------------------
  pred_y <- function(x_val) {
    nd <- setNames(data.frame(x_val), independent_variable)
    as.numeric(predict(fit, newdata = nd))
  }

  # --- helper: SE of predicted y via delta method -------------------------
  #     grad(theta) evaluated numerically; se = sqrt(grad' V grad)
  #     theta, p, rhs, V are captured from enclosing scope (hoisted)
  pred_se <- function(x_val, eps = 1e-6) {
    y0   <- pred_y(x_val)
    nd   <- setNames(data.frame(x_val), independent_variable)
    grad <- vapply(seq_len(p), function(j) {
      theta_j    <- theta
      theta_j[j] <- theta[j] + eps
      env <- c(as.list(theta_j), as.list(nd))
      (as.numeric(eval(rhs, envir = env)) - y0) / eps
    }, numeric(1))
    sqrt(as.numeric(crossprod(grad, V %*% grad)))
  }

  # --- helper: safe uniroot wrapper --------------------------------------
  safe_uniroot <- function(f, lower, upper) {
    tryCatch({
      f_lo <- f(lower)
      f_hi <- f(upper)
      if (is.na(f_lo) || is.na(f_hi)) return(NA_real_)
      if (sign(f_lo) == sign(f_hi)) return(NA_real_)
      uniroot(f, lower = lower, upper = upper, tol = .Machine$double.eps^0.5)$root
    }, error = function(e) NA_real_)
  }

  # --- 1. mindc: fitted curve == llod -------------------------------------
  mindc <- NA_real_
  if (!is.na(llod)) {
    mindc <- safe_uniroot(function(x_val) pred_y(x_val) - llod,
                          lower = x_lo, upper = x_hi)
  }

  # --- 2. maxdc: fitted curve == ulod -------------------------------------
  maxdc <- NA_real_
  if (!is.na(ulod)) {
    maxdc <- safe_uniroot(function(x_val) pred_y(x_val) - ulod,
                          lower = x_lo, upper = x_hi)
  }

  # --- 3. minrdl: 2.5% CI of fitted curve == llod ------------------------
  #        lower CI = pred_y(x) - t * se(x)
  minrdl <- NA_real_
  if (!is.na(llod)) {
    minrdl <- safe_uniroot(
      function(x_val) (pred_y(x_val) - t_crit * pred_se(x_val)) - llod,
      lower = x_lo, upper = x_hi
    )
  }

  # --- 4. maxrdl: 97.5% CI of fitted curve == ulod -----------------------
  #        upper CI = pred_y(x) + t * se(x)
  maxrdl <- NA_real_
  if (!is.na(ulod)) {
    maxrdl <- safe_uniroot(
      function(x_val) (pred_y(x_val) + t_crit * pred_se(x_val)) - ulod,
      lower = x_lo, upper = x_hi
    )
  }

  if (verbose) {
    message(sprintf("MDC/RDL - mindc: %s, maxdc: %s, minrdl: %s, maxrdl: %s",
                    format(mindc, digits = 4), format(maxdc, digits = 4),
                    format(minrdl, digits = 4), format(maxrdl, digits = 4)))
  }

  list(
    mindc  = as.numeric(mindc),
    maxdc  = as.numeric(maxdc),
    minrdl = as.numeric(minrdl),
    maxrdl = as.numeric(maxrdl)
  )
}


#' Compute Shape/curvature-based Lower and Upper Limits of Quantification (LOQ)
#' based on the second derivative
#'
#' Identifies the lower limit of quantification (LLOQ) and upper limit of
#' quantification (ULOQ) based on local extrema of the second derivative
#' of a fitted curve. Candidate extrema are refined using quadratic
#' interpolation, and corresponding response values are obtained from
#' the fitted model.
#'
#' @param best_d2xy A data.frame containing second derivative (found in best fit object) information with:
#' \describe{
#'   \item{x}{Numeric vector of independent variable values.}
#'   \item{d2x_y}{Numeric vector of second derivative values evaluated at \code{x}.}
#' }
#' @param fit A fitted model object supporting \code{predict()}.
#' @param independent_variable Character string naming the predictor variable
#' used in the model.
#' @param verbose Logical; if \code{TRUE}, prints diagnostic messages.
#'
#' @return A named list with:
#' \describe{
#'   \item{lloq}{Estimated lower limit of quantification (x-value).}
#'   \item{uloq}{Estimated upper limit of quantification (x-value).}
#'   \item{lloq_y}{Predicted response at LLOQ.}
#'   \item{uloq_y}{Predicted response at ULOQ.}
#' }
#'
#' @details
#' The method proceeds as follows:
#' \itemize{
#'   \item Computes first differences of the second derivative to detect sign changes.
#'   \item Identifies candidate local maxima and minima of the second derivative.
#'   \item Applies quadratic interpolation using three neighboring points to refine
#'         the location of each extremum.
#'   \item Selects:
#'     \itemize{
#'       \item LLOQ as the x-value corresponding to the maximum interpolated second derivative
#'       \item ULOQ as the x-value corresponding to the minimum interpolated second derivative
#'     }
#'   \item Evaluates the fitted model at these x-values to obtain corresponding y-values.
#' }
#'
#' If no valid extrema are found, the corresponding outputs are returned as \code{NA}.
#'
#' @examples
#' \dontrun{
#' loqs <- compute_loqs(
#'   best_d2xy = d2_data,
#'   fit = model_fit,
#'   independent_variable = "log_dilution"
#' )
#'
#' loqs$lloq
#' loqs$uloq
#' }
#'
#' @seealso \code{\link{generate_lods}}, \code{\link{generate_mdc_rdl}}
#'
#' @export
compute_loqs <- function(best_d2xy, fit, independent_variable, verbose = TRUE) {
  y <- as.numeric(best_d2xy$d2x_y)
  x <- as.numeric(best_d2xy$x)

  n <- length(y)
  if (n < 3) stop("Need at least 3 points to detect local extrema.")

  # first differences
  dy <- diff(y)

  # candidate interior indices where slope changes sign
  idx_max <- which(dy[-1] < 0 & dy[-length(dy)] > 0) + 1  # local max neighborhood
  idx_min <- which(dy[-1] > 0 & dy[-length(dy)] < 0) + 1  # local min neighborhood

  interpolate_vertex <- function(i) {
    # use points (i-1, i, i+1)
    xi <- x[(i-1):(i+1)]
    yi <- y[(i-1):(i+1)]

    # fit quadratic: y = a*x^2 + b*x + c
    X <- cbind(xi^2, xi, 1)
    coef <- solve(t(X) %*% X, t(X) %*% yi)  # least squares for robustness
    a <- coef[1]; b <- coef[2]; c <- coef[3]

    if (a == 0) {
      # fallback: no curvature, just return middle point
      return(list(x = xi[2], y = yi[2]))
    }

    # vertex of parabola: x* = -b / (2a)
    xv <- -b / (2 * a)

    # clamp to local interval [x_{i-1}, x_{i+1}] to avoid nonsense extrapolation
    xv <- max(min(xv, max(xi)), min(xi))

    yv <- a * xv^2 + b * xv + c
    list(x = xv, y = yv)
  }

  # interpolate for each candidate
  if (length(idx_max) > 0) {
    max_list <- lapply(idx_max, interpolate_vertex)
    max_df <- data.frame(
      x = vapply(max_list, `[[`, numeric(1), "x"),
      y = vapply(max_list, `[[`, numeric(1), "y"),
      i_center = idx_max
    )
  } else {
    max_df <- data.frame(x = numeric(0), y = numeric(0), i_center = integer(0))
  }

  if (length(idx_min) > 0) {
    min_list <- lapply(idx_min, interpolate_vertex)
    min_df <- data.frame(
      x = vapply(min_list, `[[`, numeric(1), "x"),
      y = vapply(min_list, `[[`, numeric(1), "y"),
      i_center = idx_min
    )
  } else {
    min_df <- data.frame(x = numeric(0), y = numeric(0), i_center = integer(0))
  }

  # Handle empty dataframes and ensure scalar values are returned
  if (nrow(max_df) > 0) {
    lloq_x <- as.numeric(max_df[which.max(max_df$y), "x"][1])
  } else {
    lloq_x <- NA_real_
  }

  if (nrow(min_df) > 0) {
    uloq_x <- as.numeric(min_df[which.min(min_df$y), "x"][1])
  } else {
    uloq_x <- NA_real_
  }

  y_loq <- tryCatch({
    predict(fit, newdata = setNames(data.frame(x = c(lloq_x, uloq_x)), independent_variable))
  }, error = function(e) rep(NA_real_, 2))

  # Ensure lloq_y and uloq_y are scalar numeric values
  lloq_y <- if (length(y_loq) >= 1 && !all(is.na(y_loq))) as.numeric(min(y_loq, na.rm = TRUE)) else NA_real_
  uloq_y <- if (length(y_loq) >= 1 && !all(is.na(y_loq))) as.numeric(max(y_loq, na.rm = TRUE)) else NA_real_

  return(list(
    lloq = lloq_x,
    uloq = uloq_x,
    lloq_y = lloq_y,
    uloq_y = uloq_y
  ))
}

#' Summarize Fitted Model with QC Metrics, LOD, Detection Limits, and Curvature LOQs
#'
#' Generates a one-row summary of a fitted model, including parameter estimates,
#' goodness-of-fit statistics, inflection point coordinates, limits of detection (LOD),
#' minimum/maximum detectable concentrations (MDC), reliable detection limits (RDL),
#' curvature-based limits of quantification (LOQ), and metadata. The result is stored
#' in \code{best_fit$best_fit_summary}.
#'
#' This function is designed to operate on the output of a model fitting step
#' and produce a consistent summary structure for downstream analysis,
#' QC filtering, or reporting.
#'
#' @param best_fit A list containing model results, including:
#'   \itemize{
#'     \item \code{best_fit$best_fit}: Fitted model object
#'     \item \code{best_fit$best_data}: Data used for fitting
#'     \item \code{best_fit$best_model_name}: the model name of the best fit (selected by minimizing the AIC score)
#'     \item \code{best_fit$best_d2xy}: data.frame of second derivative values
#'   }
#' @param response_variable Character string naming the response variable in \code{best_data} (e.g. mfi, absorbance).
#' @param independent_variable Character string naming the predictor variable used in the model.
#' @param fixed_a_result numeric value for parameter \code{a}. This a derived result from the \code{\link{select_antigen_plate}} function.
#' @param antigen_settings List of antigen-specific settings. May include:
#'   \itemize{
#'     \item \code{std_error_blank}: Standard error of the blank used for LLOD calculation
#'   }
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
#'   \item Limits of detection:
#'     \itemize{
#'       \item \code{llod}: Lower limit of detection
#'       \item \code{ulod}: Upper limit of detection
#'     }
#'   \item Detection limits derived from the fitted curve:
#'     \itemize{
#'       \item \code{mindc}, \code{maxdc}: Minimum and maximum detectable concentrations
#'       \item \code{minrdl}, \code{maxrdl}: Reliable detection limits based on confidence intervals
#'     }
#'   \item Curvature-based limits of quantification (LOQ):
#'     \itemize{
#'       \item \code{lloq}, \code{uloq}: LOQ values based on extrema of the second derivative
#'       \item \code{lloq_y}, \code{uloq_y}: Corresponding predicted response values
#'     }
#'   \item Metadata: source, transformation flags, and model formula
#' }
#'
#' The inflection point is computed using \code{\link{compute_inflection_point}}
#' and corresponds to the location where the second derivative of the fitted curve equals zero.
#'
#' Limits of detection (LOD) are computed using \code{\link{generate_lods}}, optionally
#' incorporating blank standard error. MDC and RDL values are computed using
#' \code{\link{generate_mdc_rdl}} via root-finding over the observed data range.
#'
#' Curvature-based LOQs are computed using \code{\link{compute_loqs}}, which identifies
#' local extrema in the second derivative of the fitted curve and refines their locations
#' using quadratic interpolation.
#'
#' If no valid model, derivative data, or extrema are available, corresponding values
#' are returned as \code{NA}.
#'
#' @return A list identical to \code{best_fit} with additional elements:
#' \itemize{
#'   \item \code{best_fit_summary}: A one-row \code{data.frame} containing model summary and QC metrics
#'   \item \code{lods}: A list with \code{llod} and \code{ulod}
#'   \item \code{mdc_rdl}: A list with \code{mindc}, \code{maxdc}, \code{minrdl}, and \code{maxrdl}
#'   \item \code{curv_loqs}: A list with \code{lloq}, \code{uloq}, \code{lloq_y}, and \code{uloq_y}
#' }
#'
#' @examples
#' \dontrun{
#' best_fit <- fit_model(data)
#'
#' best_fit <- summarize_fit(
#'   best_fit = best_fit,
#'   response_variable = "mfi",
#'   independent_variable = "concentration",
#'   fixed_a_result = NULL,
#'   antigen_settings = list(std_error_blank = 0.05),
#'   antigen_fit_options = list(
#'     blank_option = "none",
#'     is_log_response = FALSE,
#'     is_log_concentration = TRUE,
#'     apply_prozone = FALSE
#'   )
#' )
#'
#' best_fit$best_fit_summary
#' best_fit$lods
#' best_fit$mdc_rdl
#' best_fit$curv_loqs
#' }
#'
#' @seealso \code{\link{compute_inflection_point}},
#'   \code{\link{generate_lods}},
#'   \code{\link{generate_mdc_rdl}},
#'   \code{\link{compute_loqs}}
#'
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

  # --- LODs -------------------
  lods <- tryCatch({
    generate_lods(
      best_fit = best_fit,
      fixed_a_result = fixed_a_result,
      std_error_blank = antigen_settings$std_error_blank %||% NA_real_,
      verbose = verbose
    )
  }, error = function(e) {
    if (verbose) message("[summarize_fit] LOD calculation failed: ", e$message)
    list(llod = NA_real_, ulod = NA_real_)
  })

  # ---- MDC / RDL -------------------------------
  mdc_rdl <- tryCatch({
    generate_mdc_rdl(
      best_fit = best_fit,
      lods = lods,
      independent_variable = independent_variable,
      verbose = verbose
    )
  }, error = function(e) {
    if (verbose) message("[summarize_fit] MDC/RDL calculation failed: ", e$message)
    list(mindc = NA_real_, maxdc = NA_real_, minrdl = NA_real_, maxrdl = NA_real_)
  })

  # --- Curvature LOQs -------
  curv_loqs <- tryCatch({
    if (!is.null(best_fit$best_d2xy) && nrow(best_fit$best_d2xy) >= 3) {
      compute_loqs(
        best_d2xy            = best_fit$best_d2xy,
        fit                  = fit,
        independent_variable = independent_variable,
        verbose              = verbose
      )
    } else {
      if (verbose) message("[summarize_fit] No best_d2xy available for curvature LOQs")
      list(lloq = NA_real_, uloq = NA_real_, lloq_y = NA_real_, uloq_y = NA_real_)
    }
  }, error = function(e) {
    if (verbose) message("[summarize_fit] compute_loqs error: ", conditionMessage(e))
    list(lloq = NA_real_, uloq = NA_real_, lloq_y = NA_real_, uloq_y = NA_real_)
  })


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
    llod = lods$llod,
    ulod = lods$ulod,
    mindc = mdc_rdl$mindc,
    maxdc  = mdc_rdl$maxdc,
    minrdl = mdc_rdl$minrdl,
    maxrdl = mdc_rdl$maxrdl,
    lloq  = curv_loqs$lloq,
    uloq  = curv_loqs$uloq,
    lloq_y = curv_loqs$lloq_y,
    uloq_y = curv_loqs$uloq_y,
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
  # glance_df$source      <- safe_unique(best_data$source)
  glance_df$bkg_method  <- antigen_fit_options$blank_option %||% NA_character_
  glance_df$is_log_response <- antigen_fit_options$is_log_response %||% NA
  glance_df$is_log_x    <- antigen_fit_options$is_log_concentration %||% NA
  glance_df$apply_prozone <- antigen_fit_options$apply_prozone %||% NA
  glance_df$formula     <- tryCatch(sub("I\\((.*)\\)", "\\1", paste(deparse(formula(fit)), collapse = " ")), error = function(e) NA_character_)
  glance_df$last_concentration_calc_method <- "interpolated"


  best_fit$best_fit_summary <- glance_df

  return(best_fit)
}
