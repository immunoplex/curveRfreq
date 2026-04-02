# ---- model_functions.R ----


#' Four-Parameter Logistic (4PL) Function
#'
#' Computes the response \eqn{y} for a four-parameter logistic curve:
#' \deqn{y = d + \frac{a - d}{1 + \exp\!\left(\frac{x - c}{b}\right)}}
#'
#' @param x Numeric vector. Independent variable (typically log-concentration).
#' @param a Numeric scalar. Lower (left) asymptote — baseline response at
#'   zero concentration or infinite dilution.
#' @param b Numeric scalar. Scale (slope) parameter controlling steepness
#'   of the transition region. Sign determines curve direction.
#' @param c Numeric scalar. Inflection-point location on the x-axis (midpoint).
#' @param d Numeric scalar. Upper (right) asymptote — maximum response at
#'   saturating concentration.
#'
#' @return Numeric vector of predicted response values, same length as `x`.
#'
#' @details
#' The classical 4PL model is the workhorse of immunoassay quantitation
#' (ELISA, Luminex bead arrays). It is symmetric about its inflection
#' point; for asymmetric data see [logistic5()].
#'
#' @family forward-models
#' @seealso [logistic5()], [loglogistic4()], [gompertz4()] for alternative parameterisations;
#'   [inv_logistic4()] for the inverse; [dydxlogistic4()] for the first derivative.
#'
#' @examples
#' x <- seq(-2, 6, length.out = 200)
#' y <- logistic4(x, a = 100, b = 1.5, c = 2, d = 20000)
#' plot(x, y, type = "l", main = "4PL Standard Curve")
#'
#' @export
logistic4 <- function(x, a, b, c, d) {
  d + (a - d) / (1 + exp((x - c) / b))
}


#' First Derivative of the 4PL Model
#'
#' Computes \eqn{dy/dx} analytically for [logistic4()]:
#' \deqn{\frac{dy}{dx} = -\frac{(a - d)\,u}{b\,(1 + u)^2}
#'   \quad\text{where } u = \exp\!\left(\frac{x - c}{b}\right)}
#'
#' @inheritParams logistic4
#'
#' @return Numeric vector of first-derivative values, same length as `x`.
#'
#' @family derivatives
#' @seealso [logistic4()], [d2xlogistic4()]
#'
#' @examples
#' x <- seq(-2, 6, length.out = 200)
#' slope <- dydxlogistic4(x, a = 100, b = 1.5, c = 2, d = 20000)
#' plot(x, slope, type = "l", main = "4PL dy/dx")
#'
#' @export
dydxlogistic4 <- function(x, a, b, c, d) {
  u <- exp((x - c) / b)
  -(a - d) * u / (b * (1 + u)^2)
}


#' Second Derivative of the 4PL Model
#'
#' Computes \eqn{d^2y/dx^2} analytically for [logistic4()]:
#' \deqn{\frac{d^2y}{dx^2} = \frac{(a - d)\,u\,(u - 1)}{b^2\,(1 + u)^3}}
#'
#' @inheritParams logistic4
#'
#' @return Numeric vector of second-derivative values, same length as `x`.
#'
#' @family derivatives
#' @seealso [logistic4()], [dydxlogistic4()]
#'
#' @keywords internal
d2xlogistic4 <- function(x, a, b, c, d) {
  u <- exp((x - c) / b)
  (a - d) * u * (u - 1) / (b^2 * (1 + u)^3)
}


#' Five-Parameter Logistic (5PL) Function
#'
#' Computes the response for a five-parameter logistic curve with
#' asymmetry parameter `g`:
#' \deqn{y = d + \frac{a - d}{\left(1 +
#'   \exp\!\left(\frac{x - c}{b}\right)\right)^{\!g}}}
#'
#' When \eqn{g = 1} this reduces to [logistic4()].
#'
#' @inheritParams logistic4
#' @param g Numeric scalar. Asymmetry parameter. Values \eqn{> 1} skew the
#'   curve toward the upper asymptote; values \eqn{< 1} skew it toward the
#'   lower asymptote.
#'
#' @return Numeric vector of predicted response values.
#'
#' @family forward-models
#' @seealso [logistic4()], [inv_logistic5()], [dydxlogistic5()]
#'
#' @examples
#' x <- seq(-2, 6, length.out = 200)
#' y <- logistic5(x, a = 80, b = 1.2, c = 2.5, d = 18000, g = 1.3)
#' plot(x, y, type = "l", main = "5PL Standard Curve")
#'
#' @export
logistic5 <- function(x, a, b, c, d, g) {
  d + (a - d) / (1 + exp((x - c) / b))^g
}


#' First Derivative of the 5PL Model
#'
#' Computes \eqn{dy/dx} analytically for [logistic5()]:
#' \deqn{\frac{dy}{dx} = -\frac{g\,(a - d)\,u}{b\,(1 + u)^{g+1}}}
#'
#' @inheritParams logistic5
#'
#' @return Numeric vector of first-derivative values.
#'
#' @family derivatives
#' @seealso [logistic5()], [d2xlogistic5()]
#'
#' @export
dydxlogistic5 <- function(x, a, b, c, d, g) {
  u <- exp((x - c) / b)
  -g * (a - d) * u / (b * (1 + u)^(g + 1))
}


#' Second Derivative of the 5PL Model (Numerical)
#'
#' Approximates \eqn{d^2y/dx^2} for [logistic5()] using a central-difference
#' finite-difference scheme.
#'
#' @inheritParams logistic5
#' @param h Numeric scalar. Step size for the finite difference.
#'   Default `1e-5`.
#'
#' @return Numeric vector of second-derivative values.
#'
#' @family derivatives
#' @seealso [logistic5()], [dydxlogistic5()]
#'
#' @keywords internal
d2xlogistic5 <- function(x, a, b, c, d, g, h = 1e-5) {
  (logistic5(x + h, a, b, c, d, g) - 2 * logistic5(x, a, b, c, d, g) +
     logistic5(x - h, a, b, c, d, g)) / h^2
}


#' Four-Parameter Dose–Response (loglogistic4) Forward Function
#'
#' An alternative 4-parameter logistic parameterisation where the
#' inflection point `c` is on the concentration scale (not log-scale)
#' and the slope `b` acts as a Hill coefficient:
#' \deqn{y = a + \frac{d - a}{1 + (x / c)^b}}
#'
#' @inheritParams logistic4
#'
#' @return Numeric vector of predicted response values.
#'
#' @family forward-models
#' @seealso [logistic4()], [loglogistic5()], [inv_loglogistic4()], [dydxloglogistic4()]
#'
#' @examples
#' x <- seq(0.01, 100, length.out = 200)
#' y <- loglogistic4(x, a = 50, b = -2, c = 10, d = 15000)
#' plot(x, y, type = "l", log = "x", main = "loglogistic4 Dose-Response")
#'
#' @export
loglogistic4 <- function(x, a, b, c, d) {
  a + (d - a) / (1 + (x / c)^b)
}


#' First Derivative of the loglogistic4 Model
#'
#' Computes \eqn{dy/dx} analytically for [loglogistic4()]:
#' \deqn{\frac{dy}{dx} = -\frac{b\,(d - a)\,r}{x\,(1 + r)^2}
#'   \quad\text{where } r = (x/c)^b}
#'
#' @inheritParams loglogistic4
#'
#' @return Numeric vector of first-derivative values.
#'
#' @family derivatives
#' @seealso [loglogistic4()], [d2xloglogistic4()]
#'
#' @export
dydxloglogistic4 <- function(x, a, b, c, d) {
  r <- (x / c)^b
  -b * (d - a) * r / (x * (1 + r)^2)
}


#' Second Derivative of the loglogistic4 Model (Numerical)
#'
#' Approximates \eqn{d^2y/dx^2} for [loglogistic4()] via central difference.
#'
#' @inheritParams loglogistic4
#' @param h Numeric scalar. Finite-difference step size. Default `1e-5`.
#'
#' @return Numeric vector of second-derivative values.
#'
#' @family derivatives
#' @keywords internal
d2xloglogistic4 <- function(x, a, b, c, d, h = 1e-5) {
  (loglogistic4(x + h, a, b, c, d) - 2 * loglogistic4(x, a, b, c, d) +
     loglogistic4(x - h, a, b, c, d)) / h^2
}


#' Five-Parameter Dose–Response (loglogistic5) Forward Function
#'
#' A five-parameter generalised logistic (Richards) curve in dose–response
#' form with asymmetry parameter `g`:
#' \deqn{y = a + (d - a)\,\bigl(1 + g\,\exp(-b\,(x - c))\bigr)^{-1/g}}
#'
#' When \eqn{g = 1} this reduces to a standard 4PL dose–response.
#'
#' @inheritParams logistic4
#' @param g Numeric scalar. Asymmetry (Richards) parameter.
#'
#' @return Numeric vector of predicted response values.
#'
#' @family forward-models
#' @seealso [loglogistic4()], [logistic5()], [inv_loglogistic5()], [dydxloglogistic5()]
#'
#' @export
loglogistic5 <- function(x, a, b, c, d, g) {
  a + (d - a) * (1 + g * exp(-b * (x - c)))^(-1 / g)
}


#' First Derivative of the loglogistic5 Model
#'
#' Computes \eqn{dy/dx} analytically for [loglogistic5()]:
#' \deqn{\frac{dy}{dx} = b\,(d - a)\,\exp(-b(x-c))\,
#'   u^{-1/g - 1}
#'   \quad\text{where } u = 1 + g\,\exp(-b(x-c))}
#'
#' @inheritParams loglogistic5
#'
#' @return Numeric vector of first-derivative values.
#'
#' @family derivatives
#' @seealso [loglogistic5()], [d2xloglogistic5()]
#'
#' @export
dydxloglogistic5 <- function(x, a, b, c, d, g) {
  u <- 1 + g * exp(-b * (x - c))
  b * (d - a) * exp(-b * (x - c)) * u^(-1 / g - 1)
}


#' Second Derivative of the loglogistic5 Model (Numerical)
#'
#' Approximates \eqn{d^2y/dx^2} for [loglogistic5()] via central difference.
#'
#' @inheritParams loglogistic5
#' @param h Numeric scalar. Finite-difference step size. Default `1e-5`.
#'
#' @return Numeric vector of second-derivative values.
#'
#' @family derivatives
#' @keywords internal
d2xloglogistic5 <- function(x, a, b, c, d, g, h = 1e-5) {
  (loglogistic5(x + h, a, b, c, d, g) - 2 * loglogistic5(x, a, b, c, d, g) +
     loglogistic5(x - h, a, b, c, d, g)) / h^2
}


#' Four-Parameter Gompertz Forward Function
#'
#' Computes the response for a Gompertz growth/saturation curve:
#' \deqn{y = a + (d - a)\,\exp\!\bigl(-\exp(-b\,(x - c))\bigr)}
#'
#' The Gompertz is intrinsically asymmetric — unlike the 4PL it does not
#' require a fifth parameter for asymmetry.
#'
#' @inheritParams logistic4
#'
#' @return Numeric vector of predicted response values.
#'
#' @family forward-models
#' @seealso [logistic4()], [loglogistic4()], [inv_gompertz4()], [dydxgompertz4()]
#'
#' @examples
#' x <- seq(-2, 8, length.out = 200)
#' y <- gompertz4(x, a = 50, b = 1, c = 3, d = 15000)
#' plot(x, y, type = "l", main = "Gompertz Curve")
#'
#' @export
gompertz4 <- function(x, a, b, c, d) {
  a + (d - a) * exp(-exp(-b * (x - c)))
}


#' First Derivative of the Gompertz Model
#'
#' Computes \eqn{dy/dx} analytically for [gompertz4()]:
#' \deqn{\frac{dy}{dx} = b\,(d - a)\,u\,\exp(-u)
#'   \quad\text{where } u = \exp(-b(x - c))}
#'
#' @inheritParams gompertz4
#'
#' @return Numeric vector of first-derivative values.
#'
#' @family derivatives
#' @seealso [gompertz4()], [d2xgompertz4()]
#'
#' @export
dydxgompertz4 <- function(x, a, b, c, d) {
  u <- exp(-b * (x - c))
  b * (d - a) * u * exp(-u)
}


#' Second Derivative of the Gompertz Model
#'
#' Computes \eqn{d^2y/dx^2} analytically for [gompertz4()]:
#' \deqn{\frac{d^2y}{dx^2} = b^2\,(d - a)\,e_2\,(e_2 - 1)\,\exp(-e_2)
#'   \quad\text{where } e_2 = \exp(-b(x - c))}
#'
#' @inheritParams gompertz4
#'
#' @return Numeric vector of second-derivative values.
#'
#' @family derivatives
#' @keywords internal
d2xgompertz4 <- function(x, a, b, c, d) {
  e2 <- exp(-(b * (x - c)))
  b^2 * (d - a) * e2 * (e2 - 1) * exp(-e2)
}


# ============================================================================
# INVERSE FUNCTIONS
# ============================================================================

#' Inverse of the 4PL Model
#'
#' Given observed response `y`, solves analytically for `x` (concentration):
#' \deqn{x = c + b\,\log\!\left(\frac{a - d}{y - d} - 1\right)}
#'
#' Values of `y` outside the open interval
#' \eqn{(\min(a,d) + \mathrm{tol},\; \max(a,d) - \mathrm{tol})} are
#' undefined and return `NA`.
#'
#' @param y Numeric vector. Observed response values.
#' @param a Numeric scalar. Lower asymptote (free parameter from fit).
#' @param b Numeric scalar. Scale/slope parameter.
#' @param c Numeric scalar. Inflection-point location.
#' @param d Numeric scalar. Upper asymptote.
#' @param tol Numeric scalar. Buffer from asymptotes to prevent
#'   log-domain errors. Default `1e-6`.
#'
#' @return Numeric vector of estimated `x` values, same length as `y`.
#'   `NA` for out-of-range `y`.
#'
#' @family inverse-functions
#' @seealso [logistic4()], [inv_logistic4_fixed()] for the fixed-\eqn{a} variant.
#'
#' @examples
#' # Round-trip: forward then inverse recovers x
#' x <- seq(-1, 5, length.out = 50)
#' y <- logistic4(x, a = 100, b = 1.5, c = 2, d = 20000)
#' x_hat <- inv_logistic4(y, a = 100, b = 1.5, c = 2, d = 20000)
#' all.equal(x, x_hat)
#'
#' @export
inv_logistic4 <- function(y, a, b, c, d, tol = 1e-6) {
  y_min <- min(a, d) + tol
  y_max <- max(a, d) - tol
  result <- rep(NA_real_, length(y))
  valid  <- !is.na(y) & y > y_min & y < y_max
  if (any(valid))
    result[valid] <- c + b * log((a - d) / (y[valid] - d) - 1)
  result
}


#' Inverse of the 4PL Model with Fixed Lower Asymptote
#'
#' Same algebra as [inv_logistic4()] but parameter `a` is supplied as a known
#' constant (`fixed_a`) rather than estimated from the fit. No domain
#' checking is performed — caller must ensure `y` is in range.
#'
#' @param y Numeric vector. Observed response values.
#' @param fixed_a Numeric scalar. Externally fixed lower asymptote.
#' @param b,c,d Numeric scalars. Model parameters from [stats::coef()].
#'
#' @return Numeric vector of estimated `x` values.
#'
#' @family inverse-functions
#' @seealso [inv_logistic4()], [grad_inv_logistic4_fixed()]
#'
#' @keywords internal
inv_logistic4_fixed <- function(y, fixed_a, b, c, d) {
  c + b * log((fixed_a - d) / (y - d) - 1)
}


#' Inverse of the loglogistic4 Model
#'
#' Solves for `x` given response `y` under the [loglogistic4()] parameterisation:
#' \deqn{x = c\,\left(\frac{d - a}{y - a} - 1\right)^{1/b}}
#'
#' @param y Numeric vector. Observed response values.
#' @param a,b,c,d Numeric scalars. loglogistic4 model parameters.
#'
#' @return Numeric vector of estimated `x` values.
#'
#' @family inverse-functions
#' @seealso [loglogistic4()], [inv_loglogistic4_fixed()]
#'
#' @export
inv_loglogistic4 <- function(y, a, b, c, d) {
  c * (((d - a) / (y - a)) - 1)^(1 / b)
}


#' Inverse of the loglogistic4 Model with Fixed Lower Asymptote
#'
#' @inheritParams inv_logistic4_fixed
#' @param b,c,d Numeric scalars. Model parameters.
#'
#' @return Numeric vector of estimated `x` values.
#'
#' @family inverse-functions
#' @keywords internal
inv_loglogistic4_fixed <- function(y, fixed_a, b, c, d) {
  c * (((d - fixed_a) / (y - fixed_a)) - 1)^(1 / b)
}


#' Inverse of the Gompertz Model
#'
#' Solves for `x` given response `y` under [gompertz4()]:
#' \deqn{x = c - \frac{1}{b}\,\log\!\left(-\log\frac{y - a}{d - a}\right)}
#'
#' @param y Numeric vector. Observed response values.
#' @param a,b,c,d Numeric scalars. Gompertz model parameters.
#'
#' @return Numeric vector of estimated `x` values.
#'
#' @family inverse-functions
#' @seealso [gompertz4()], [inv_gompertz4_fixed()]
#'
#' @export
inv_gompertz4 <- function(y, a, b, c, d) {
  c - (1 / b) * log(-log((y - a) / (d - a)))
}


#' Inverse of the Gompertz Model with Fixed Lower Asymptote
#'
#' @inheritParams inv_logistic4_fixed
#' @param b,c,d Numeric scalars. Model parameters.
#'
#' @return Numeric vector of estimated `x` values.
#'
#' @family inverse-functions
#' @keywords internal
inv_gompertz4_fixed <- function(y, fixed_a, b, c, d) {
  c - (1 / b) * log(-log((y - fixed_a) / (d - fixed_a)))
}


#' Inverse of the 5PL Model
#'
#' Solves for `x` given response `y` under [logistic5()]:
#' \deqn{x = c + b\,\log\!\left(\left(\frac{a - d}{y - d}\right)^{1/g}
#'   - 1\right)}
#'
#' @param y Numeric vector. Observed response values.
#' @param a,b,c,d Numeric scalars. 5PL model parameters.
#' @param g Numeric scalar. Asymmetry parameter.
#'
#' @return Numeric vector of estimated `x` values.
#'
#' @family inverse-functions
#' @seealso [logistic5()], [inv_logistic5_fixed()]
#'
#' @export
inv_logistic5 <- function(y, a, b, c, d, g) {
  c + b * log(((a - d) / (y - d))^(1 / g) - 1)
}


#' Inverse of the 5PL Model with Fixed Lower Asymptote
#'
#' @inheritParams inv_logistic4_fixed
#' @param b,c,d Numeric scalars. Model parameters.
#' @param g Numeric scalar. Asymmetry parameter.
#'
#' @return Numeric vector of estimated `x` values.
#'
#' @family inverse-functions
#' @keywords internal
inv_logistic5_fixed <- function(y, fixed_a, b, c, d, g) {
  c + b * log(((fixed_a - d) / (y - d))^(1 / g) - 1)
}


#' Inverse of the loglogistic5 Model
#'
#' Solves for `x` given response `y` under [loglogistic5()]:
#' \deqn{x = c - \frac{1}{b}\left(\log\!\left(\left(\frac{y - a}{d - a}
#'   \right)^{-g} - 1\right) - \log g\right)}
#'
#' @param y Numeric vector. Observed response values.
#' @param a,b,c,d Numeric scalars. loglogistic5 model parameters.
#' @param g Numeric scalar. Asymmetry parameter.
#'
#' @return Numeric vector of estimated `x` values.
#'
#' @family inverse-functions
#' @seealso [loglogistic5()], [inv_loglogistic5_fixed()]
#'
#' @export
inv_loglogistic5 <- function(y, a, b, c, d, g) {
  c - (1 / b) * (log(((y - a) / (d - a))^(-g) - 1) - log(g))
}


#' Inverse of the loglogistic5 Model with Fixed Lower Asymptote
#'
#' @inheritParams inv_logistic4_fixed
#' @param b,c,d Numeric scalars. Model parameters.
#' @param g Numeric scalar. Asymmetry parameter.
#'
#' @return Numeric vector of estimated `x` values.
#'
#' @family inverse-functions
#' @keywords internal
inv_loglogistic5_fixed <- function(y, fixed_a, b, c, d, g) {
  c - (1 / b) * (log(((y - fixed_a) / (d - fixed_a))^(-g) - 1) - log(g))
}


# ============================================================================
# dx/dy HELPERS (used in propagate_error_analytic)
# ============================================================================

#' \eqn{dx/dy} for the 4PL Inverse
#'
#' Scalar derivative of the [inv_logistic4()] function with respect to `y`,
#' used by [propagate_error_analytic()] for error propagation through
#' the measurement uncertainty component.
#'
#' @param y Numeric scalar. Observed response.
#' @param a,b,c,d Numeric scalars. 4PL model parameters.
#'
#' @return Numeric scalar \eqn{dx/dy}.
#'
#' @family dxdy-helpers
#' @keywords internal
dxdyilogistic4 <- function(y, a, b, c, d) {
  ratio <- (a - d) / (y - d)
  b * ratio / ((ratio - 1) * (y - d))
}


#' \eqn{dx/dy} for the loglogistic4 Inverse
#'
#' @inheritParams dxdyilogistic4
#'
#' @return Numeric scalar \eqn{dx/dy}.
#'
#' @family dxdy-helpers
#' @keywords internal
dxdyiloglogistic4 <- function(y, a, b, c, d) {
  inner <- ((d - a) / (y - a)) - 1
  -c * (1 / b) * inner^(1 / b - 1) * (d - a) / (y - a)^2
}


#' \eqn{dx/dy} for the Gompertz Inverse
#'
#' @param y Numeric scalar. Observed response.
#' @param a,b,c,d Numeric scalars. Gompertz model parameters.
#'
#' @return Numeric scalar \eqn{dx/dy}.
#'
#' @family dxdy-helpers
#' @keywords internal
dxdyigompertz4 <- function(y, a, b, c, d) {
  ratio <- (y - a) / (d - a)
  (1 / b) / ((y - a) * log(ratio))
}


#' \eqn{dx/dy} for the 5PL Inverse
#'
#' @param y Numeric scalar. Observed response.
#' @param a,b,c,d Numeric scalars. 5PL model parameters.
#' @param g Numeric scalar. Asymmetry parameter.
#'
#' @return Numeric scalar \eqn{dx/dy}.
#'
#' @family dxdy-helpers
#' @keywords internal
dxdyilogistic5 <- function(y, a, b, c, d, g) {
  ratio <- (a - d) / (y - d)
  b * (1 / g) * ratio^(1 / g) / ((ratio^(1 / g) - 1) * (y - d))
}


#' \eqn{dx/dy} for the loglogistic5 Inverse
#'
#' @param y Numeric scalar. Observed response.
#' @param a,b,c,d Numeric scalars. loglogistic5 model parameters.
#' @param g Numeric scalar. Asymmetry parameter.
#'
#' @return Numeric scalar \eqn{dx/dy}.
#'
#' @family dxdy-helpers
#' @keywords internal
dxdyiloglogistic5 <- function(y, a, b, c, d, g) {
  ratio <- (y - a) / (d - a)
  (1 / b) * g * ratio^(-g) / ((ratio^(-g) - 1) * (y - a))
}


# ============================================================================
# FULL GRADIENT FUNCTIONS (free a) — return ALL partials
# ============================================================================

#' Analytical Gradient of the Inverse 4PL
#'
#' Returns partial derivatives \eqn{\partial x / \partial \theta} and
#' \eqn{\partial x / \partial y} for the [inv_logistic4()] function when all
#' four parameters are free. Consumed by [propagate_error_analytic()]
#' for delta-method error propagation.
#'
#' @param y Numeric scalar. Observed response value.
#' @param a,b,c,d Numeric scalars. Free 4PL model parameters.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{grad_theta}{Named numeric vector `c(a=, b=, c=, d=)` —
#'       partial derivatives of `x` with respect to each model parameter.}
#'     \item{grad_y}{Numeric scalar — \eqn{\partial x / \partial y}.}
#'   }
#'
#' @family gradient-functions
#' @seealso [inv_logistic4()], [grad_inv_logistic4_fixed()]
#'
#' @export
grad_logistic4 <- function(y, a, b, c, d) {
  A  <- (a - d) / (y - d) - 1
  da <-  b / (A * (y - d))
  db <-  log(A)
  dc <-  1
  dd <-  b * (a - y) / (A * (y - d)^2)
  dy <- -b * (a - d) / (A * (y - d)^2)
  list(grad_theta = c(a = da, b = db, c = dc, d = dd), grad_y = dy)
}


#' Analytical Gradient of the Inverse loglogistic4
#'
#' Returns partial derivatives for [inv_loglogistic4()] with all parameters free.
#'
#' @param y Numeric scalar. Observed response.
#' @param a,b,c,d Numeric scalars. Free loglogistic4 model parameters.
#'
#' @return A list with `grad_theta` (named vector `c(a=, b=, c=, d=)`)
#'   and scalar `grad_y`.
#'
#' @family gradient-functions
#' @seealso [inv_loglogistic4()], [grad_inv_loglogistic4_fixed()]
#'
#' @export
grad_loglogistic4 <- function(y, a, b, c, d) {
  Q  <- ((d - a) / (y - a)) - 1
  p  <- 1 / b
  x  <- c * Q^p
  da <-  (c / b) * Q^(p - 1) * (d - y) / (y - a)^2
  db <- -x * log(Q) / b^2
  dc <-  Q^p
  dd <-  (c / b) * Q^(p - 1) / (y - a)
  dy <- -(c / b) * Q^(p - 1) * (d - a) / (y - a)^2
  list(grad_theta = c(a = da, b = db, c = dc, d = dd), grad_y = dy)
}


#' Analytical Gradient of the Inverse Gompertz
#'
#' Returns partial derivatives for [inv_gompertz4()] with all parameters free.
#'
#' @param y Numeric scalar. Observed response.
#' @param a,b,c,d Numeric scalars. Free Gompertz model parameters.
#'
#' @return A list with `grad_theta` (named vector `c(a=, b=, c=, d=)`)
#'   and scalar `grad_y`.
#'
#' @family gradient-functions
#' @seealso [inv_gompertz4()], [grad_inv_gompertz4_fixed()]
#'
#' @export
grad_gompertz4 <- function(y, a, b, c, d) {
  R  <- (y - a) / (d - a)
  B  <- -log(R)
  da <-  1 / (b * B * (y - a))
  db <-  log(B) / b^2
  dc <-  1
  dd <- -1 / (b * B * (d - a))
  dy <-  1 / (b * B * (y - a))
  list(grad_theta = c(a = da, b = db, c = dc, d = dd), grad_y = dy)
}


#' Analytical Gradient of the Inverse 5PL
#'
#' Returns partial derivatives for [inv_logistic5()] with all five parameters free.
#'
#' @param y Numeric scalar. Observed response.
#' @param a,b,c,d Numeric scalars. Free 5PL model parameters.
#' @param g Numeric scalar. Asymmetry parameter.
#'
#' @return A list with `grad_theta` (named vector `c(a=, b=, c=, d=, g=)`)
#'   and scalar `grad_y`.
#'
#' @family gradient-functions
#' @seealso [inv_logistic5()], [grad_inv_logistic5_fixed()]
#'
#' @export
grad_logistic5 <- function(y, a, b, c, d, g) {
  H  <- ((a - d) / (y - d))^(1 / g) - 1
  da <-  b / (g * H * (a - d))
  db <-  log(H)
  dc <-  1
  dd <-  b * (-1 / (g * H * (y - d)) - (a - d) / ((y - d)^2 * H))
  dg <- -b * log((a - d) / (y - d)) / (g^2 * H)
  dy <-  b * (a - d) / (g * H * (y - d)^2)
  list(grad_theta = c(a = da, b = db, c = dc, d = dd, g = dg), grad_y = dy)
}


#' Analytical Gradient of the Inverse loglogistic5
#'
#' Returns partial derivatives for [inv_loglogistic5()] with all five parameters free.
#'
#' @param y Numeric scalar. Observed response.
#' @param a,b,c,d Numeric scalars. Free loglogistic5 model parameters.
#' @param g Numeric scalar. Asymmetry parameter.
#'
#' @return A list with `grad_theta` (named vector `c(a=, b=, c=, d=, g=)`)
#'   and scalar `grad_y`.
#'
#' @family gradient-functions
#' @seealso [inv_loglogistic5()], [grad_inv_loglogistic5_fixed()]
#'
#' @export
grad_loglogistic5 <- function(y, a, b, c, d, g) {
  Q  <- ((y - a) / (d - a))^(-g) - 1
  da <-  (1 / (b * Q * g)) * ((y - a) / (d - a))^(-g - 1) * (1 / (d - a))
  db <-  (log(g) - log(Q)) / b^2
  dc <-  1
  dd <-  (1 / (b * Q * g)) * ((y - a) / (d - a))^(-g - 1) * (y - a) / (d - a)^2
  dg <-  (1 / (b * g^2)) * (1 - log(g * Q))
  dy <-  ((y - a) / (d - a))^(-g - 1) / (b * g * Q * (d - a))
  list(grad_theta = c(a = da, b = db, c = dc, d = dd, g = dg), grad_y = dy)
}


# ============================================================================
# FIXED-a GRADIENT FUNCTIONS (only b, c, d [, g] partials)
# ============================================================================

#' Gradient of the Inverse 4PL with Fixed Lower Asymptote
#'
#' Returns partial derivatives \eqn{\partial x / \partial \theta} for
#' [inv_logistic4_fixed()] where `a` is an externally fixed constant and only
#' `b`, `c`, `d` are free.
#'
#' @param y Numeric scalar. Observed response.
#' @param fixed_a Numeric scalar. Externally fixed lower asymptote.
#' @param b,c,d Numeric scalars. Free model parameters.
#'
#' @return Named numeric vector `c(b=, c=, d=)` of partial derivatives.
#'
#' @family fixed-a-gradients
#' @seealso [grad_logistic4()], [inv_logistic4_fixed()], [grad_y_logistic4_fixed()]
#'
#' @keywords internal
grad_inv_logistic4_fixed <- function(y, fixed_a, b, c, d) {
  y <- as.numeric(y); fixed_a <- as.numeric(fixed_a)
  b <- as.numeric(b); c <- as.numeric(c); d <- as.numeric(d)
  inside <- (fixed_a - d) / (y - d) - 1
  db <- log(inside)
  dc <- 1.0
  dd <- (b / (fixed_a - y)) * (1 / (fixed_a - d) - 1 / (y - d))
  c(b = db, c = dc, d = dd)
}


#' Gradient of the Inverse loglogistic4 with Fixed Lower Asymptote
#'
#' @inheritParams grad_inv_logistic4_fixed
#'
#' @return Named numeric vector `c(b=, c=, d=)`.
#'
#' @family fixed-a-gradients
#' @keywords internal
grad_inv_loglogistic4_fixed <- function(y, fixed_a, b, c, d) {
  y <- as.numeric(y); fixed_a <- as.numeric(fixed_a)
  b <- as.numeric(b); c <- as.numeric(c); d <- as.numeric(d)
  inside <- (d - fixed_a) / (y - fixed_a) - 1
  Z  <- inside^(1 / b)
  db <- -c * Z * log(inside) / (b^2)
  dc <- Z
  dd <- (c / b) * Z / ((d - fixed_a) * inside)
  c(b = db, c = dc, d = dd)
}


#' Gradient of the Inverse Gompertz with Fixed Lower Asymptote
#'
#' @inheritParams grad_inv_logistic4_fixed
#'
#' @return Named numeric vector `c(b=, c=, d=)`.
#'
#' @family fixed-a-gradients
#' @keywords internal
grad_inv_gompertz4_fixed <- function(y, fixed_a, b, c, d) {
  y <- as.numeric(y); fixed_a <- as.numeric(fixed_a)
  b <- as.numeric(b); c <- as.numeric(c); d <- as.numeric(d)
  L  <- -log((y - fixed_a) / (d - fixed_a))
  db <-  log(L) / (b^2)
  dc <-  1.0
  dd <- -1.0 / (b * L * (d - fixed_a))
  c(b = db, c = dc, d = dd)
}


#' Gradient of the Inverse 5PL with Fixed Lower Asymptote
#'
#' @param y Numeric scalar. Observed response.
#' @param fixed_a Numeric scalar. Externally fixed lower asymptote.
#' @param b,c,d Numeric scalars. Free model parameters.
#' @param g Numeric scalar. Asymmetry parameter (free).
#'
#' @return Named numeric vector `c(b=, c=, d=, g=)`.
#'
#' @family fixed-a-gradients
#' @keywords internal
grad_inv_logistic5_fixed <- function(y, fixed_a, b, c, d, g) {
  y <- as.numeric(y); fixed_a <- as.numeric(fixed_a)
  b <- as.numeric(b); c <- as.numeric(c)
  d <- as.numeric(d); g <- as.numeric(g)
  ratio   <- (fixed_a - d) / (y - d)
  ratio_g <- ratio^(1 / g)
  U  <- ratio_g - 1
  db <- log(U)
  dc <- 1.0
  dd <- (b / (g * U)) * ratio_g * (1 / (y - d) - 1 / (fixed_a - d))
  dg <- -(b / g^2) * log(ratio) * ratio_g / U
  c(b = db, c = dc, d = dd, g = dg)
}


#' Gradient of the Inverse loglogistic5 with Fixed Lower Asymptote
#'
#' @inheritParams grad_inv_logistic5_fixed
#'
#' @return Named numeric vector `c(b=, c=, d=, g=)`.
#'
#' @family fixed-a-gradients
#' @keywords internal
grad_inv_loglogistic5_fixed <- function(y, fixed_a, b, c, d, g) {
  y <- as.numeric(y); fixed_a <- as.numeric(fixed_a)
  b <- as.numeric(b); c <- as.numeric(c)
  d <- as.numeric(d); g <- as.numeric(g)
  ratio   <- (y - fixed_a) / (d - fixed_a)
  ratio_g <- ratio^(-g)
  V  <- ratio_g - 1
  db <- (log(V) - log(g)) / (b^2)
  dc <- 1.0
  dd <- -(1 / b) * (ratio_g / V) * (g / (d - fixed_a))
  dg <- (1 / b) * ((ratio_g / V) * log(ratio) - 1 / g)
  c(b = db, c = dc, d = dd, g = dg)
}


# ============================================================================
# FIXED-a grad_y FUNCTIONS (∂x/∂y when a is externally fixed)
# ============================================================================

#' \eqn{\partial x / \partial y} for the 4PL Inverse with Fixed \eqn{a}
#'
#' Computes the derivative of the inverse 4PL with respect to the
#' observed response `y` when the lower asymptote is externally fixed.
#' Used for the measurement-uncertainty component of error propagation.
#'
#' @param y Numeric scalar. Observed response.
#' @param fixed_a Numeric scalar. Externally fixed lower asymptote.
#' @param b Numeric scalar. Slope parameter.
#' @param d Numeric scalar. Upper asymptote.
#'
#' @return Numeric scalar \eqn{\partial x / \partial y}.
#'
#' @family fixed-a-grad-y
#' @keywords internal
grad_y_logistic4_fixed <- function(y, fixed_a, b, d) {
  y <- as.numeric(y); fixed_a <- as.numeric(fixed_a)
  b <- as.numeric(b); d <- as.numeric(d)
  inside <- (fixed_a - d) / (y - d) - 1
  -b * (fixed_a - d) / ((y - d)^2 * inside)
}


#' \eqn{\partial x / \partial y} for the loglogistic4 Inverse with Fixed \eqn{a}
#'
#' @param y Numeric scalar. Observed response.
#' @param fixed_a Numeric scalar. Externally fixed lower asymptote.
#' @param b,c,d Numeric scalars. Free model parameters.
#'
#' @return Numeric scalar \eqn{\partial x / \partial y}.
#'
#' @family fixed-a-grad-y
#' @keywords internal
grad_y_loglogistic4_fixed <- function(y, fixed_a, b, c, d) {
  y <- as.numeric(y); fixed_a <- as.numeric(fixed_a)
  b <- as.numeric(b); c <- as.numeric(c); d <- as.numeric(d)
  inside <- (d - fixed_a) / (y - fixed_a) - 1
  -c / b * inside^(1 / b - 1) * (d - fixed_a) / (y - fixed_a)^2
}


#' \eqn{\partial x / \partial y} for the Gompertz Inverse with Fixed \eqn{a}
#'
#' @param y Numeric scalar. Observed response.
#' @param fixed_a Numeric scalar. Externally fixed lower asymptote.
#' @param b Numeric scalar. Slope parameter.
#' @param d Numeric scalar. Upper asymptote.
#'
#' @return Numeric scalar \eqn{\partial x / \partial y}.
#'
#' @family fixed-a-grad-y
#' @keywords internal
grad_y_gompertz4_fixed <- function(y, fixed_a, b, d) {
  y <- as.numeric(y); fixed_a <- as.numeric(fixed_a)
  b <- as.numeric(b); d <- as.numeric(d)
  ratio <- (y - fixed_a) / (d - fixed_a)
  L <- -log(ratio)
  1.0 / (b * L * (y - fixed_a))
}


#' \eqn{\partial x / \partial y} for the 5PL Inverse with Fixed \eqn{a}
#'
#' @param y Numeric scalar. Observed response.
#' @param fixed_a Numeric scalar. Externally fixed lower asymptote.
#' @param b Numeric scalar. Slope parameter.
#' @param d Numeric scalar. Upper asymptote.
#' @param g Numeric scalar. Asymmetry parameter.
#'
#' @return Numeric scalar \eqn{\partial x / \partial y}.
#'
#' @family fixed-a-grad-y
#' @keywords internal
grad_y_logistic5_fixed <- function(y, fixed_a, b, d, g) {
  y <- as.numeric(y); fixed_a <- as.numeric(fixed_a)
  b <- as.numeric(b); d <- as.numeric(d); g <- as.numeric(g)
  ratio   <- (fixed_a - d) / (y - d)
  ratio_g <- ratio^(1 / g)
  U <- ratio_g - 1
  b * (fixed_a - d) / (g * U * (y - d)^2)
}


#' \eqn{\partial x / \partial y} for the loglogistic5 Inverse with Fixed \eqn{a}
#'
#' @param y Numeric scalar. Observed response.
#' @param fixed_a Numeric scalar. Externally fixed lower asymptote.
#' @param b Numeric scalar. Slope parameter.
#' @param d Numeric scalar. Upper asymptote.
#' @param g Numeric scalar. Asymmetry parameter.
#'
#' @return Numeric scalar \eqn{\partial x / \partial y}.
#'
#' @family fixed-a-grad-y
#' @keywords internal
grad_y_loglogistic5_fixed <- function(y, fixed_a, b, d, g) {
  y <- as.numeric(y); fixed_a <- as.numeric(fixed_a)
  b <- as.numeric(b); d <- as.numeric(d); g <- as.numeric(g)
  ratio   <- (y - fixed_a) / (d - fixed_a)
  ratio_g <- ratio^(-g)
  V <- ratio_g - 1
  g / (b * (d - fixed_a)) * ratio^(-g - 1) / V
}


# ============================================================================
# DISPATCH: make_inv_and_grad_fixed
# ============================================================================

#' Build Inverse, Gradient, and grad_y Closures for a Model
#'
#' Constructs a list of three closures — `inv`, `grad`, and `grad_y` —
#' that evaluate the inverse function, its parameter gradient, and its
#' response-derivative for a given model at a specific observed response
#' `y`. The closures accept a single named parameter vector `p`
#' (typically from [stats::coef()]) and are consumed by
#' [propagate_error_analytic()].
#'
#' @section Design:
#' Two branches exist depending on whether `fixed_a` is supplied:
#' \describe{
#'   \item{`fixed_a` is non-`NULL`}{Parameter `a` is treated as an
#'     externally known constant. The `grad` closure returns partials
#'     for `b`, `c`, `d` (and `g` for 5-parameter models) only.
#'     The `_fixed` family of functions is used.}
#'   \item{`fixed_a` is `NULL`}{Parameter `a` is free and must appear
#'     in `p`. The full `grad_*()` functions are used, returning
#'     partials for `a`, `b`, `c`, `d` (and `g`).}
#' }
#'
#' @param model Character. One of `"logistic4"`, `"logistic5"`, `"loglogistic4"`, `"loglogistic5"`,
#'   `"gompertz4"`.
#' @param y Numeric scalar. The observed response value at which to
#'   evaluate the inverse and its derivatives.
#' @param fixed_a Numeric scalar or `NULL`. If non-`NULL`, the lower
#'   asymptote is treated as a fixed external constant.
#'
#' @return A list with three closures:
#'   \describe{
#'     \item{`inv(p)`}{Returns the estimated `x` (concentration) given
#'       named parameter vector `p`.}
#'     \item{`grad(p)`}{Returns a named numeric vector of
#'       \eqn{\partial x / \partial \theta_i}.}
#'     \item{`grad_y(p)`}{Returns a numeric scalar
#'       \eqn{\partial x / \partial y}.}
#'   }
#'
#' @seealso [propagate_error_analytic()], [inv_logistic4()], [grad_logistic4()],
#'   [inv_logistic4_fixed()], [grad_inv_logistic4_fixed()]
#'
#' @examples
#' # Free-a example
#' fns <- make_inv_and_grad_fixed("logistic4", y = 5000, fixed_a = NULL)
#' p <- c(a = 100, b = 1.5, c = 2, d = 20000)
#' fns$inv(p)
#' fns$grad(p)
#' fns$grad_y(p)
#'
#' # Fixed-a example
#' fns2 <- make_inv_and_grad_fixed("logistic4", y = 5000, fixed_a = 100)
#' p2 <- c(b = 1.5, c = 2, d = 20000)
#' fns2$inv(p2)
#' fns2$grad(p2)
#'
#' @export
make_inv_and_grad_fixed <- function(model, y, fixed_a) {
  # Strip any accidental names from scalar inputs
  y <- as.numeric(y)

  # ── Branch A: 'a' is a TRUE external constant ──────────────
  if (!is.null(fixed_a)) {
    fixed_a <- as.numeric(fixed_a)
    return(switch(model,
                  logistic4 = list(
                    inv    = function(p) inv_logistic4_fixed(y, fixed_a,
                                                      as.numeric(p["b"]), as.numeric(p["c"]), as.numeric(p["d"])),
                    grad   = function(p) grad_inv_logistic4_fixed(y, fixed_a,
                                                           as.numeric(p["b"]), as.numeric(p["c"]), as.numeric(p["d"])),
                    grad_y = function(p) grad_y_logistic4_fixed(y, fixed_a,
                                                         as.numeric(p["b"]), as.numeric(p["d"]))
                  ),
                  loglogistic4 = list(
                    inv    = function(p) inv_loglogistic4_fixed(y, fixed_a,
                                                       as.numeric(p["b"]), as.numeric(p["c"]), as.numeric(p["d"])),
                    grad   = function(p) grad_inv_loglogistic4_fixed(y, fixed_a,
                                                            as.numeric(p["b"]), as.numeric(p["c"]), as.numeric(p["d"])),
                    grad_y = function(p) grad_y_loglogistic4_fixed(y, fixed_a,
                                                          as.numeric(p["b"]), as.numeric(p["c"]), as.numeric(p["d"]))
                  ),
                  gompertz4 = list(
                    inv    = function(p) inv_gompertz4_fixed(y, fixed_a,
                                                          as.numeric(p["b"]), as.numeric(p["c"]), as.numeric(p["d"])),
                    grad   = function(p) grad_inv_gompertz4_fixed(y, fixed_a,
                                                               as.numeric(p["b"]), as.numeric(p["c"]), as.numeric(p["d"])),
                    grad_y = function(p) grad_y_gompertz4_fixed(y, fixed_a,
                                                             as.numeric(p["b"]), as.numeric(p["d"]))
                  ),
                  logistic5 = list(
                    inv    = function(p) inv_logistic5_fixed(y, fixed_a,
                                                      as.numeric(p["b"]), as.numeric(p["c"]),
                                                      as.numeric(p["d"]), as.numeric(p["g"])),
                    grad   = function(p) grad_inv_logistic5_fixed(y, fixed_a,
                                                           as.numeric(p["b"]), as.numeric(p["c"]),
                                                           as.numeric(p["d"]), as.numeric(p["g"])),
                    grad_y = function(p) grad_y_logistic5_fixed(y, fixed_a,
                                                         as.numeric(p["b"]), as.numeric(p["d"]), as.numeric(p["g"]))
                  ),
                  loglogistic5 = list(
                    inv    = function(p) inv_loglogistic5_fixed(y, fixed_a,
                                                       as.numeric(p["b"]), as.numeric(p["c"]),
                                                       as.numeric(p["d"]), as.numeric(p["g"])),
                    grad   = function(p) grad_inv_loglogistic5_fixed(y, fixed_a,
                                                            as.numeric(p["b"]), as.numeric(p["c"]),
                                                            as.numeric(p["d"]), as.numeric(p["g"])),
                    grad_y = function(p) grad_y_loglogistic5_fixed(y, fixed_a,
                                                          as.numeric(p["b"]), as.numeric(p["d"]), as.numeric(p["g"]))
                  ),
                  stop("Unsupported model: ", model)
    ))
  }

  # ── Branch B: 'a' is FREE — must be in p (coef(fit)) ───────
  switch(model,
         logistic4 = list(
           inv    = function(p) {
             a <- as.numeric(p["a"])
             inv_logistic4_fixed(y, a, as.numeric(p["b"]),
                          as.numeric(p["c"]), as.numeric(p["d"]))
           },
           grad   = function(p) {
             grads <- grad_logistic4(y,
                              a = as.numeric(p["a"]), b = as.numeric(p["b"]),
                              c = as.numeric(p["c"]), d = as.numeric(p["d"]))
             grads$grad_theta
           },
           grad_y = function(p) {
             grad_logistic4(y,
                     a = as.numeric(p["a"]), b = as.numeric(p["b"]),
                     c = as.numeric(p["c"]), d = as.numeric(p["d"]))$grad_y
           }
         ),
         loglogistic4 = list(
           inv    = function(p) {
             a <- as.numeric(p["a"])
             inv_loglogistic4_fixed(y, a, as.numeric(p["b"]),
                           as.numeric(p["c"]), as.numeric(p["d"]))
           },
           grad   = function(p) {
             grads <- grad_loglogistic4(y,
                               a = as.numeric(p["a"]), b = as.numeric(p["b"]),
                               c = as.numeric(p["c"]), d = as.numeric(p["d"]))
             grads$grad_theta
           },
           grad_y = function(p) {
             grad_loglogistic4(y,
                      a = as.numeric(p["a"]), b = as.numeric(p["b"]),
                      c = as.numeric(p["c"]), d = as.numeric(p["d"]))$grad_y
           }
         ),
         gompertz4 = list(
           inv    = function(p) {
             a <- as.numeric(p["a"])
             inv_gompertz4_fixed(y, a, as.numeric(p["b"]),
                              as.numeric(p["c"]), as.numeric(p["d"]))
           },
           grad   = function(p) {
             grads <- grad_gompertz4(y,
                                  a = as.numeric(p["a"]), b = as.numeric(p["b"]),
                                  c = as.numeric(p["c"]), d = as.numeric(p["d"]))
             grads$grad_theta
           },
           grad_y = function(p) {
             grad_gompertz4(y,
                         a = as.numeric(p["a"]), b = as.numeric(p["b"]),
                         c = as.numeric(p["c"]), d = as.numeric(p["d"]))$grad_y
           }
         ),
         logistic5 = list(
           inv    = function(p) {
             a <- as.numeric(p["a"])
             inv_logistic5_fixed(y, a, as.numeric(p["b"]), as.numeric(p["c"]),
                          as.numeric(p["d"]), as.numeric(p["g"]))
           },
           grad   = function(p) {
             grads <- grad_logistic5(y,
                              a = as.numeric(p["a"]), b = as.numeric(p["b"]),
                              c = as.numeric(p["c"]), d = as.numeric(p["d"]),
                              g = as.numeric(p["g"]))
             grads$grad_theta
           },
           grad_y = function(p) {
             grad_logistic5(y,
                     a = as.numeric(p["a"]), b = as.numeric(p["b"]),
                     c = as.numeric(p["c"]), d = as.numeric(p["d"]),
                     g = as.numeric(p["g"]))$grad_y
           }
         ),
         loglogistic5 = list(
           inv    = function(p) {
             a <- as.numeric(p["a"])
             inv_loglogistic5_fixed(y, a, as.numeric(p["b"]), as.numeric(p["c"]),
                           as.numeric(p["d"]), as.numeric(p["g"]))
           },
           grad   = function(p) {
             grads <- grad_loglogistic5(y,
                               a = as.numeric(p["a"]), b = as.numeric(p["b"]),
                               c = as.numeric(p["c"]), d = as.numeric(p["d"]),
                               g = as.numeric(p["g"]))
             grads$grad_theta
           },
           grad_y = function(p) {
             grad_loglogistic5(y,
                      a = as.numeric(p["a"]), b = as.numeric(p["b"]),
                      c = as.numeric(p["c"]), d = as.numeric(p["d"]),
                      g = as.numeric(p["g"]))$grad_y
           }
         ),
         stop("Unsupported model: ", model)
  )
}


# ============================================================================
# MODEL FORMULA SELECTION
# ============================================================================

#' Select Model Formulas for Standard Curve Fitting
#'
#' Returns a named list of [stats::nls()]-compatible formulas for all
#' candidate models (logistic5, loglogistic5, logistic4, loglogistic4, gompertz4). When `fixed_constraint`
#' is supplied, parameter `a` (lower asymptote) is substituted as a
#' numeric constant in every formula; otherwise `a` is a free parameter.
#'
#' @param fixed_constraint Numeric or `NULL`. If non-`NULL`, the lower
#'   asymptote `a` is fixed to this value in every formula. Typically
#'   derived from blank-well measurements. Values that are non-positive
#'   or non-finite are rejected with a message, falling back to free-`a`
#'   formulas.
#' @param response_variable Character. Column name for the assay response
#'   (e.g., `"mfi"`, `"od"`). Used as the left-hand side of each formula.
#' @param is_log_response Logical. If `TRUE`, `fixed_constraint` is
#'   \eqn{\log_{10}}-transformed before substitution (with a small
#'   epsilon to avoid \eqn{\log(0)}).
#'
#' @param model_names a vector of model formulas to be considered.
#' @return Named list of formulas keyed by model name:
#'   `"logistic5"`, `"loglogistic5"`, `"logistic4"`, `"loglogistic4"`, `"gompertz4"`.
#'
#' @details
#' The concentration column is hard-coded as `concentration` in all
#' formulas. The `I()` wrappers protect complex sub-expressions from
#' formula parsing by [stats::nls()].
#'
#' @seealso [logistic4()], [logistic5()], [loglogistic4()], [loglogistic5()], [gompertz4()] for the
#'   underlying model functions; [compute_robust_curves()] which
#'   consumes these formulas.
#'
#' @examples
#' # Fixed lower asymptote
#' forms <- select_model_formulas(fixed_constraint = 50,
#'                                response_variable = "mfi",
#'                                is_log_response = FALSE,
#'                                model_names = c("logistic4", "logistic5",
#'                                "loglogistic4", "loglogistic5",
#'                                "gompertz4"))
#' names(forms)
#' forms$logistic4
#'
#' # Free lower asymptote
#' forms_free <- select_model_formulas(fixed_constraint = NULL,
#'                                     response_variable = "mfi",
#'                                     is_log_response = FALSE,
#'                                     model_names = c("logistic4", "logistic5",
#'                                      "loglogistic4", "loglogistic5",
#'                                       "gompertz4"))
#' forms_free$logistic4
#'
#' @export
select_model_formulas <- function(fixed_constraint, response_variable, is_log_response, model_names) {

  if (!is.null(fixed_constraint)) {
    if (!is.numeric(fixed_constraint) ||
        !is.finite(fixed_constraint) ||
        fixed_constraint <= 0) {
      message(sprintf(
        "[select_model_formulas] fixed_constraint = %s is invalid (<=0 or non-finite). Using free-a formulas.",
        if (is.numeric(fixed_constraint)) format(fixed_constraint) else "non-numeric"
      ))
      fixed_constraint <- NULL
    }
  }

  if (!is.null(fixed_constraint)) {
    if (is_log_response) {
      .eps <- 0.00005
      fixed_constraint <- log10(fixed_constraint + .eps)
    }
    message("Lower asymptote is fixed at", fixed_constraint)
    fixed_value <- fixed_constraint

    standard_curve_formulas <- list(
      logistic5 = as.formula(substitute(
        lhs ~ d + (((fixed_a) - d) / (1 + exp((concentration - c) / b))^g),
        list(lhs = as.name(response_variable), fixed_a = fixed_value)
      ), env = parent.frame()),

      loglogistic5 = as.formula(substitute(
        lhs ~ (fixed_a) + (d - (fixed_a)) * I((1 + g * exp(-b * (concentration - c)))^(-1 / g)),
        list(lhs = as.name(response_variable), fixed_a = fixed_value)
      ), env = parent.frame()),

      logistic4 = as.formula(substitute(
        lhs ~ d + (((fixed_a) - d) / I((1 + exp((concentration - c) / b)))),
        list(lhs = as.name(response_variable), fixed_a = fixed_value)
      ), env = parent.frame()),

      loglogistic4 = as.formula(substitute(
        lhs ~ (fixed_a) + (d - (fixed_a)) / I((1 + (concentration / c)^b)),
        list(lhs = as.name(response_variable), fixed_a = fixed_value)
      ), env = parent.frame()),

      gompertz4 = as.formula(substitute(
        lhs ~ (fixed_a) + (d - (fixed_a)) * I(exp(-exp(-b * (concentration - c)))),
        list(lhs = as.name(response_variable), fixed_a = fixed_value)
      ), env = parent.frame())
    )

  } else {
    standard_curve_formulas <- list(
      logistic5 = as.formula(substitute(
        lhs ~ d + (a - d) / (I((1 + exp((concentration - c) / b))^g)),
        list(lhs = as.name(response_variable))),
        env = parent.frame()),

      loglogistic5 = as.formula(substitute(
        lhs ~ a + (d - a) * I((1 + g * exp(-b * (concentration - c)))^(-1 / g)),
        list(lhs = as.name(response_variable))),
        env = parent.frame()),

      logistic4 = as.formula(substitute(
        lhs ~ d + (a - d) / I(1 + exp((concentration - c) / b)),
        list(lhs = as.name(response_variable))),
        env = parent.frame()),

      loglogistic4 = as.formula(substitute(
        lhs ~ a + (d - a) / I(1 + (concentration / c)^b),
        list(lhs = as.name(response_variable))),
        env = parent.frame()),

      gompertz4 = as.formula(substitute(
        lhs ~ a + (d - a) * I(exp(-exp(-b * (concentration - c)))),
        list(lhs = as.name(response_variable))),
        env = parent.frame())
    )
  }

  model_name_map <- c(
    "logistic5" = "logistic5",
    "logistic4" = "logistic4",
    "gompertz4" = "gompertz4",
    "loglogistic5" = "loglogistic5",
    "loglogistic4" = "loglogistic4"
  )

  # Replace user-supplied names with internal ones
  model_names_internal <- unname(model_name_map[model_names])

  # Subset formulas to only the models specified
  standard_curve_formulas <- standard_curve_formulas[model_names_internal]

  return(standard_curve_formulas)
}


# ============================================================================
# SAFE CONSTRAINT FUNCTIONS
# ============================================================================

#' Compute Safe Parameter Bounds for 5PL (logistic5) Fitting
#'
#' Builds lower and upper bounds for all free parameters in the [logistic5()]
#' model, adapted to the observed data scale via
#' [adaptive_constraint_profile()]. Used by [compute_robust_curves()]
#' to supply bounds to [stats::nls()] with `algorithm = "port"`.
#'
#' @param data Data frame. Must contain a `concentration` column and
#'   the response column.
#' @param y_min Numeric. Minimum observed response. Default `1`.
#' @param y_max Numeric. Maximum observed response.
#' @param logistic5_formula Formula. The logistic5 model formula from
#'   [select_model_formulas()].
#' @param logistic5_free_vars Character vector. Names of free parameters in
#'   the formula (e.g., `c("b", "c", "d", "g")` when `a` is fixed,
#'   or `c("a", "b", "c", "d", "g")` when free).
#' @param is_log_response Logical. Whether the response is
#'   \eqn{\log_{10}}-transformed.
#' @param is_log_concentration Logical. Whether concentration is on
#'   log scale.
#' @param antigen_settings List. Must contain `l_asy_min_constraint`
#'   and `l_asy_max_constraint`.
#' @param constraint_profile List or `NULL`. Output of
#'   [adaptive_constraint_profile()]. If `NULL`, one is built internally.
#'
#' @return A list with elements `lower` and `upper`, each a named
#'   numeric vector matching `logistic5_free_vars`.
#'
#' @family safe-constraints
#' @seealso [adaptive_constraint_profile()], [logistic5()],
#'   [compute_robust_curves()]
#'
#' @keywords internal
logistic5_safe_constraint <- function(data,
                               y_min = 1,
                               y_max,
                               logistic5_formula,
                               logistic5_free_vars,
                               is_log_response,
                               is_log_concentration,
                               antigen_settings,
                               constraint_profile = NULL) {
  .eps <- 1e-5

  if (is.null(constraint_profile)) {
    constraint_profile <- adaptive_constraint_profile(
      data, response_variable = names(data)[1],
      is_log_response = is_log_response,
      antigen_settings = antigen_settings
    )
  }

  .slope_max <- constraint_profile$slope_max
  .slope_min <- constraint_profile$slope_min
  .g_min     <- constraint_profile$g_min
  .g_max     <- constraint_profile$g_max

  formula_vars <- all.vars(logistic5_formula)

  # ── Lower asymptote a ──
  a_lower <- antigen_settings$l_asy_min_constraint
  a_upper <- antigen_settings$l_asy_max_constraint
  if (is_log_response) {
    a_lower <- log10(a_lower + .eps)
    a_upper <- log10(a_upper + .eps)
  }

  # ── Upper asymptote d — adaptive ──
  dr     <- constraint_profile$dynamic_range
  margin <- constraint_profile$d_margin_frac
  d_lower <- y_min + dr * (0.5 - margin)
  d_lower <- max(d_lower, y_min + .eps)
  d_upper <- y_max + dr * margin
  if (d_upper <= d_lower) d_upper <- d_lower + abs(dr) * 0.5 + .eps

  # ── Inflection point c ──
  midpoint_concentration <- mean(range(data$concentration, na.rm = TRUE))
  conc_range <- diff(range(data$concentration, na.rm = TRUE))
  pad     <- constraint_profile$conc_pad_frac
  c_lower <- midpoint_concentration - pad * conc_range
  c_upper <- midpoint_concentration + pad * conc_range

  # ── Slope b ──
  b_lower <- .slope_min
  b_upper <- .slope_max

  # ── Asymmetry g ──
  g_lower <- .g_min
  g_upper <- .g_max

  if ("a" %in% formula_vars) {
    lower <- c(a = a_lower, b = b_lower, c = c_lower, d = d_lower, g = g_lower)
    upper <- c(a = a_upper, b = b_upper, c = c_upper, d = d_upper, g = g_upper)
  } else {
    lower <- c(b = b_lower, c = c_lower, d = d_lower, g = g_lower)
    upper <- c(b = b_upper, c = c_upper, d = d_lower, g = g_upper)
  }

  return(.make_bounds(logistic5_free_vars, lower_vals = lower, upper_vals = upper))
}


#' Compute Safe Parameter Bounds for loglogistic5 Fitting
#'
#' Builds lower and upper bounds for all free parameters in the [loglogistic5()]
#' model using [adaptive_constraint_profile()].
#'
#' @inheritParams logistic5_safe_constraint
#' @param loglogistic5_formula Formula. The loglogistic5 model formula.
#' @param loglogistic5_free_vars Character vector. Names of free parameters.
#'
#' @return A list with `lower` and `upper` named numeric vectors.
#'
#' @family safe-constraints
#' @keywords internal
loglogistic5_safe_constraint <- function(data, y_min = 1, y_max, loglogistic5_formula,
                                loglogistic5_free_vars, is_log_response,
                                is_log_concentration, antigen_settings,
                                constraint_profile = NULL) {
  .eps <- 1e-5

  if (is.null(constraint_profile)) {
    constraint_profile <- adaptive_constraint_profile(
      data, response_variable = names(data)[1],
      is_log_response = is_log_response,
      antigen_settings = antigen_settings
    )
  }

  .slope_max <- constraint_profile$slope_max
  .slope_min <- constraint_profile$slope_min
  .g_min     <- constraint_profile$g_min
  .g_max     <- constraint_profile$g_max

  formula_vars <- all.vars(loglogistic5_formula)

  # ── Lower asymptote a ──
  a_lower <- antigen_settings$l_asy_min_constraint
  a_upper <- antigen_settings$l_asy_max_constraint
  if (is_log_response) {
    a_lower <- log10(max(a_lower, .eps) + .eps)
    a_upper <- log10(max(a_upper, .eps) + .eps)
  }

  # ── Upper asymptote d — adaptive ──
  dr     <- constraint_profile$dynamic_range
  margin <- constraint_profile$d_margin_frac
  d_lower <- y_min + dr * (0.5 - margin)
  d_lower <- max(d_lower, y_min + .eps)
  d_upper <- y_max + dr * margin
  if (d_upper <= d_lower) d_upper <- d_lower + abs(dr) * 0.5 + .eps

  # ── Inflection point c ──
  midpoint_concentration <- mean(range(data$concentration, na.rm = TRUE))
  conc_range <- diff(range(data$concentration, na.rm = TRUE))
  pad     <- constraint_profile$conc_pad_frac
  c_lower <- midpoint_concentration - pad * conc_range
  c_upper <- midpoint_concentration + pad * conc_range

  # ── Slope b ──
  b_lower <- .slope_min
  b_upper <- .slope_max

  # ── Asymmetry g ──
  g_lower <- .g_min
  g_upper <- .g_max

  if ("a" %in% formula_vars) {
    lower <- c(a = a_lower, b = b_lower, c = c_lower, d = d_lower, g = g_lower)
    upper <- c(a = a_upper, b = b_upper, c = c_upper, d = d_upper, g = g_upper)
  } else {
    lower <- c(b = b_lower, c = c_lower, d = d_lower, g = g_lower)
    upper <- c(b = b_upper, c = c_upper, d = d_upper, g = g_upper)
  }

  return(.make_bounds(loglogistic5_free_vars, lower_vals = lower, upper_vals = upper))
}


#' Compute Safe Parameter Bounds for 4PL (logistic4) Fitting
#'
#' Builds lower and upper bounds for all free parameters in the [logistic4()]
#' model using [adaptive_constraint_profile()].
#'
#' @inheritParams logistic5_safe_constraint
#' @param logistic4_formula Formula. The logistic4 model formula.
#' @param logistic4_free_vars Character vector. Names of free parameters.
#'
#' @return A list with `lower` and `upper` named numeric vectors.
#'
#' @family safe-constraints
#' @keywords internal
logistic4_safe_constraint <- function(data, y_min = 1, y_max, logistic4_formula,
                               logistic4_free_vars, is_log_response,
                               is_log_concentration, antigen_settings,
                               constraint_profile = NULL) {
  .eps <- 1e-5

  if (is.null(constraint_profile)) {
    constraint_profile <- adaptive_constraint_profile(
      data, response_variable = names(data)[1],
      is_log_response = is_log_response,
      antigen_settings = antigen_settings
    )
  }

  .slope_max <- constraint_profile$slope_max
  .slope_min <- constraint_profile$slope_min

  formula_vars <- all.vars(logistic4_formula)

  # ── Lower asymptote a ──
  a_lower <- antigen_settings$l_asy_min_constraint
  a_upper <- antigen_settings$l_asy_max_constraint
  if (is_log_response) {
    a_lower <- log10(max(a_lower, .eps) + .eps)
    a_upper <- log10(max(a_upper, .eps) + .eps)
  }

  # ── Upper asymptote d — adaptive ──
  dr     <- constraint_profile$dynamic_range
  margin <- constraint_profile$d_margin_frac
  d_lower <- y_min + dr * (0.5 - margin)
  d_lower <- max(d_lower, y_min + .eps)
  d_upper <- y_max + dr * margin
  if (d_upper <= d_lower) d_upper <- d_lower + abs(dr) * 0.5 + .eps

  # ── Inflection point c ──
  midpoint_concentration <- mean(range(data$concentration, na.rm = TRUE))
  conc_range <- diff(range(data$concentration, na.rm = TRUE))
  pad     <- constraint_profile$conc_pad_frac
  c_lower <- midpoint_concentration - pad * conc_range
  c_upper <- midpoint_concentration + pad * conc_range

  # ── Slope b ──
  b_lower <- .slope_min
  b_upper <- .slope_max

  if ("a" %in% formula_vars) {
    lower <- c(a = a_lower, b = b_lower, c = c_lower, d = d_lower)
    upper <- c(a = a_upper, b = b_upper, c = c_upper, d = d_upper)
  } else {
    lower <- c(b = b_lower, c = c_lower, d = d_lower)
    upper <- c(b = b_upper, c = c_upper, d = d_upper)
  }

  return(.make_bounds(logistic4_free_vars, lower_vals = lower, upper_vals = upper))
}


#' Compute Safe Parameter Bounds for loglogistic4 Fitting
#'
#' Builds lower and upper bounds for all free parameters in the [loglogistic4()]
#' model using [adaptive_constraint_profile()].
#'
#' @inheritParams logistic5_safe_constraint
#' @param loglogistic4_formula Formula. The loglogistic4 model formula.
#' @param loglogistic4_free_vars Character vector. Names of free parameters.
#'
#' @return A list with `lower` and `upper` named numeric vectors.
#'
#' @family safe-constraints
#' @keywords internal
loglogistic4_safe_constraint <- function(data, y_min = 1, y_max, loglogistic4_formula,
                                loglogistic4_free_vars, is_log_response,
                                is_log_concentration, antigen_settings,
                                constraint_profile = NULL) {
  .eps <- 1e-5

  if (is.null(constraint_profile)) {
    constraint_profile <- adaptive_constraint_profile(
      data, response_variable = names(data)[1],
      is_log_response = is_log_response,
      antigen_settings = antigen_settings
    )
  }

  .slope_max <- constraint_profile$slope_max
  .slope_min <- constraint_profile$slope_min

  formula_vars <- all.vars(loglogistic4_formula)

  # ── Lower asymptote a ──
  a_lower <- antigen_settings$l_asy_min_constraint
  a_upper <- antigen_settings$l_asy_max_constraint
  if (is_log_response) {
    a_lower <- log10(max(a_lower, .eps) + .eps)
    a_upper <- log10(max(a_upper, .eps) + .eps)
  }

  # ── Upper asymptote d — adaptive ──
  dr     <- constraint_profile$dynamic_range
  margin <- constraint_profile$d_margin_frac
  d_lower <- y_min + dr * (0.5 - margin)
  d_lower <- max(d_lower, y_min + .eps)
  d_upper <- y_max + dr * margin
  if (d_upper <= d_lower) d_upper <- d_lower + abs(dr) * 0.5 + .eps

  # ── Inflection point c ──
  midpoint_concentration <- mean(range(data$concentration, na.rm = TRUE))
  conc_range <- diff(range(data$concentration, na.rm = TRUE))
  pad     <- constraint_profile$conc_pad_frac
  c_lower <- midpoint_concentration - pad * conc_range
  c_upper <- midpoint_concentration + pad * conc_range

  # ── Slope b ──
  b_lower <- .slope_min
  b_upper <- .slope_max

  if ("a" %in% formula_vars) {
    lower <- c(a = a_lower, b = b_lower, c = c_lower, d = d_lower)
    upper <- c(a = a_upper, b = b_upper, c = c_upper, d = d_upper)
  } else {
    lower <- c(b = b_lower, c = c_lower, d = d_lower)
    upper <- c(b = b_upper, c = c_upper, d = d_upper)
  }

  return(.make_bounds(loglogistic4_free_vars, lower_vals = lower, upper_vals = upper))
}


#' Compute Safe Parameter Bounds for Gompertz Fitting
#'
#' Builds lower and upper bounds for all free parameters in the
#' [gompertz4()] model using [adaptive_constraint_profile()].
#'
#' @inheritParams logistic5_safe_constraint
#' @param gompertz4_formula Formula. The Gompertz model formula.
#' @param gompertz4_free_vars Character vector. Names of free parameters.
#'
#' @return A list with `lower` and `upper` named numeric vectors.
#'
#' @family safe-constraints
#' @keywords internal
gompertz4_safe_constraint <- function(data, y_min = 1, y_max, gompertz4_formula,
                                   gompertz4_free_vars, is_log_response,
                                   is_log_concentration, antigen_settings,
                                   constraint_profile = NULL) {
  .eps <- 1e-5

  if (is.null(constraint_profile)) {
    constraint_profile <- adaptive_constraint_profile(
      data, response_variable = names(data)[1],
      is_log_response = is_log_response,
      antigen_settings = antigen_settings
    )
  }

  .slope_max <- constraint_profile$slope_max
  .slope_min <- constraint_profile$slope_min

  formula_vars <- all.vars(gompertz4_formula)

  # ── Lower asymptote a ──
  a_lower <- antigen_settings$l_asy_min_constraint
  a_upper <- antigen_settings$l_asy_max_constraint
  if (is_log_response) {
    a_lower <- log10(max(a_lower, .eps) + .eps)
    a_upper <- log10(max(a_upper, .eps) + .eps)
  }

  # ── Upper asymptote d — adaptive ──
  dr     <- constraint_profile$dynamic_range
  margin <- constraint_profile$d_margin_frac
  d_lower <- y_min + dr * (0.5 - margin)
  d_lower <- max(d_lower, y_min + .eps)
  d_upper <- y_max + dr * margin
  if (d_upper <= d_lower) d_upper <- d_lower + abs(dr) * 0.5 + .eps

  # ── Inflection point c ──
  midpoint_concentration <- mean(range(data$concentration, na.rm = TRUE))
  conc_range <- diff(range(data$concentration, na.rm = TRUE))
  pad     <- constraint_profile$conc_pad_frac
  c_lower <- midpoint_concentration - pad * conc_range
  c_upper <- midpoint_concentration + pad * conc_range

  # ── Slope b ──
  b_lower <- .slope_min
  b_upper <- .slope_max

  if ("a" %in% formula_vars) {
    lower <- c(a = a_lower, b = b_lower, c = c_lower, d = d_lower)
    upper <- c(a = a_upper, b = b_upper, c = c_upper, d = d_upper)
  } else {
    lower <- c(b = b_lower, c = c_lower, d = d_lower)
    upper <- c(b = b_upper, c = c_upper, d = d_upper)
  }

  return(.make_bounds(gompertz4_free_vars, lower_vals = lower, upper_vals = upper))
}


# ============================================================================
# ADAPTIVE CONSTRAINT PROFILE
# ============================================================================

#' Build an Adaptive Constraint Profile from Observed Data
#'
#' Inspects the response range, dynamic range, and scale to choose
#' appropriate bounds for nonlinear optimisation. The returned profile
#' is consumed by every `*_safe_constraint()` function.
#'
#' @param data Data frame with at least columns named by
#'   `response_variable` and `"concentration"`.
#' @param response_variable Character. Column name for the response.
#' @param is_log_response Logical. `TRUE` if data is already
#'   \eqn{\log_{10}}-transformed.
#' @param antigen_settings List with constraint metadata, including
#'   `l_asy_min_constraint` and `l_asy_max_constraint`.
#'
#' @return A named list:
#'   \describe{
#'     \item{y_min, y_max}{Observed response extremes (finite values only).}
#'     \item{dynamic_range}{\code{y_max - y_min}.}
#'     \item{conc_range}{Range of observed concentrations.}
#'     \item{scale_class}{`"high"` (MFI-like, e.g. Luminex),
#'       `"medium"`, or `"low"` (OD/absorbance-like).}
#'     \item{slope_min, slope_max}{Adapted bounds for scale parameter `b`.}
#'     \item{g_min, g_max}{Adapted bounds for asymmetry parameter `g`.}
#'     \item{conc_pad_frac}{Fraction beyond concentration range for
#'       inflection-point `c` bounds.}
#'     \item{d_margin_frac}{Fraction of `y_max` for upper-asymptote
#'       `d` bounds.}
#'   }
#'
#' @details
#' Three scale classes are recognised:
#'
#' | Class | Raw threshold | Log10 threshold | Typical assay |
#' |-------|--------------|-----------------|---------------|
#' | high | `y_max > 1000` | `y_max > 2.5` | Luminex MFI |
#' | medium | `y_max > 10` | `y_max > 0.5` | moderate signals |
#' | low | otherwise | otherwise | ELISA OD/absorbance |
#'
#' Narrower dynamic ranges (low scale class) receive wider slope and
#' asymmetry bounds to avoid near-singular Jacobians during [stats::nls()]
#' fitting.
#'
#' @family safe-constraints
#' @seealso [logistic5_safe_constraint()], [logistic4_safe_constraint()],
#'   [compute_robust_curves()]
#'
#' @examples
#' fake_data <- data.frame(
#'   mfi = c(100, 500, 5000, 15000, 20000),
#'   concentration = c(0.01, 0.1, 1, 10, 100)
#' )
#' profile <- adaptive_constraint_profile(
#'   data = fake_data,
#'   response_variable = "mfi",
#'   is_log_response = FALSE,
#'   antigen_settings = list(l_asy_min_constraint = 0,
#'                           l_asy_max_constraint = 200)
#' )
#' str(profile)
#'
#' @export
adaptive_constraint_profile <- function(data,
                                        response_variable,
                                        is_log_response,
                                        antigen_settings) {
  y_vals <- data[[response_variable]]
  y_vals <- y_vals[is.finite(y_vals)]
  y_min  <- min(y_vals)
  y_max  <- max(y_vals)
  dynamic_range <- y_max - y_min

  conc_vals  <- data$concentration[is.finite(data$concentration)]
  conc_range <- diff(range(conc_vals))

  # Classify the response scale
  if (is_log_response) {
    scale_class <- if (y_max > 2.5) {
      "high"
    } else if (y_max > 0.5) {
      "medium"
    } else {
      "low"
    }
  } else {
    scale_class <- if (y_max > 1000) "high" else if (y_max > 10) "medium" else "low"
  }

  # Adapt slope bounds to dynamic range
  slope_max <- switch(scale_class,
                      high   = 2.0,
                      medium = 3.0,
                      low    = 5.0
  )
  slope_min <- switch(scale_class,
                      high   = 0.1,
                      medium = 0.05,
                      low    = 0.01
  )

  # Asymmetry parameter g
  g_min <- switch(scale_class, high = 0.5, medium = 0.3, low = 0.1)
  g_max <- switch(scale_class, high = 5.0, medium = 7.0, low = 10.0)

  # How far beyond data range to allow the inflection point c
  conc_pad_frac <- switch(scale_class, high = 0.5, medium = 0.7, low = 1.0)

  # Upper asymptote d: margin above/below y_max
  d_margin_frac <- switch(scale_class, high = 0.5, medium = 0.3, low = 0.15)

  list(
    y_min          = y_min,
    y_max          = y_max,
    dynamic_range  = dynamic_range,
    conc_range     = conc_range,
    scale_class    = scale_class,
    slope_max      = slope_max,
    slope_min      = slope_min,
    g_min          = g_min,
    g_max          = g_max,
    conc_pad_frac  = conc_pad_frac,
    d_margin_frac  = d_margin_frac
  )
}
