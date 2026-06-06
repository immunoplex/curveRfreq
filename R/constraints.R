# ---- constraints.R — Adaptive parameter bounds and multi-start generation
#
# Builds lower/upper bounds for each model's free parameters based on
# observed data range and scale class, then generates multi-start lists.
# ----

#' Compute Parameter Bounds for All Candidate Models
#'
#' Inspects the response range to classify the assay scale (high/medium/low),
#' then builds adaptive bounds for a, b, c, d (and g for 5-param models).
#'
#' When `fixed_a` is NULL, the lower asymptote `a` is free and its bounds
#' are set adaptively from the data. When `fixed_a` is non-NULL, `a` does
#' not appear in the formulas (it's baked in), so no bounds are needed.
#'
#' @param data Data frame with response and concentration columns.
#' @param formulas Named list of formula objects.
#' @param response_variable Character. Response column name.
#' @param independent_variable Character. Concentration column name.
#' @param is_log_response Logical. Is the response log10-transformed?
#' @param fixed_a Numeric or NULL. If non-NULL, `a` is fixed and not
#'   bounded. If NULL, adaptive bounds for `a` are computed.
#' @param verbose Logical.
#'
#' @return Named list of lists, each with `lower` and `upper` named vectors.
#'   Has attribute `"constraint_profile"`.
#' @export
compute_model_constraints <- function(data, formulas,
                                      response_variable,
                                      independent_variable,
                                      is_log_response,
                                      fixed_a = NULL,
                                      verbose = FALSE) {

  # Adaptive profile from curveRcore
  # Pass dummy antigen_settings — only y_min/y_max/scale_class matter
  antigen_settings_for_profile <- list(
    l_asy_min_constraint = 0,
    l_asy_max_constraint = 0
  )

  profile <- curveRcore::adaptive_constraint_profile(
    data, response_variable, is_log_response, antigen_settings_for_profile
  )

  if (verbose) {
    message(sprintf("[constraints] scale=%s, dr=%.3f, slope=[%.3f,%.3f], g=[%.2f,%.2f]",
                    profile$scale_class, profile$dynamic_range,
                    profile$slope_min, profile$slope_max,
                    profile$g_min, profile$g_max))
  }

  y_min <- profile$y_min
  y_max <- profile$y_max
  dr    <- profile$dynamic_range
  eps   <- 1e-5

  # ── Lower asymptote a bounds ──
  # When fixed_a is non-NULL, 'a' is baked into the formula and won't

  # appear as a free parameter — so these bounds are irrelevant.
  # When fixed_a is NULL, 'a' is free and needs data-adaptive bounds.
  if (is.null(fixed_a)) {
    # a should be below the data minimum with some margin
    a_lo <- y_min - dr * 0.5
    a_hi <- y_min + dr * 0.1
    # Ensure a_hi > a_lo
    if (a_hi <= a_lo) a_hi <- a_lo + eps
  } else {
    # Won't be used (a not in formula), but set for safety
    a_lo <- fixed_a
    a_hi <- fixed_a
  }

  # ── Upper asymptote d bounds (adaptive) ──
  margin <- profile$d_margin_frac
  d_lo   <- max(y_min + dr * (0.5 - margin), y_min + eps)
  d_hi   <- y_max + dr * margin
  if (d_hi <= d_lo) d_hi <- d_lo + abs(dr) * 0.5 + eps

  # ── Inflection point c bounds ──
  conc_mid   <- mean(range(data[[independent_variable]], na.rm = TRUE))
  conc_range <- diff(range(data[[independent_variable]], na.rm = TRUE))
  pad        <- profile$conc_pad_frac
  c_lo       <- conc_mid - pad * conc_range
  c_hi       <- conc_mid + pad * conc_range

  # ── Slope b bounds ──
  b_lo <- profile$slope_min
  b_hi <- profile$slope_max

  # ── Asymmetry g bounds ──
  g_lo <- profile$g_min
  g_hi <- profile$g_max

  # Build per-model bounds based on which params are free in the formula
  constraints <- list()
  for (nm in names(formulas)) {
    f_vars <- setdiff(all.vars(formulas[[nm]]),
                      c(response_variable, independent_variable))
    lo <- c(a = a_lo, b = b_lo, c = c_lo, d = d_lo, g = g_lo)
    hi <- c(a = a_hi, b = b_hi, c = c_hi, d = d_hi, g = g_hi)
    constraints[[nm]] <- list(lower = lo[f_vars], upper = hi[f_vars])
  }

  attr(constraints, "constraint_profile") <- profile
  constraints
}


#' Generate Multi-Start Lists for NLS Fitting
#'
#' Creates candidate starting values evenly spread across the parameter
#' bounds. Low-signal data gets more starts with the slope biased small.
#'
#' @param model_constraints Output of [compute_model_constraints()].
#'
#' @return Named list of lists: for each model, a list of named numeric
#'   vectors suitable as `start` for `nlsLM()`.
#' @export
generate_start_lists <- function(model_constraints) {
  profile <- attr(model_constraints, "constraint_profile")
  is_low  <- !is.null(profile) && profile$scale_class == "low"
  n_starts <- if (is_low) 5L else 3L

  start_lists <- list()
  for (nm in names(model_constraints)) {
    lo <- model_constraints[[nm]]$lower
    hi <- model_constraints[[nm]]$upper
    params <- names(lo)

    starts <- lapply(seq_len(n_starts), function(i) {
      frac <- (i - 0.5) / n_starts
      s <- setNames(lo + frac * (hi - lo), params)
      if ("b" %in% params && is_low)
        s["b"] <- lo["b"] + (hi["b"] - lo["b"]) * frac^2
      as.list(s)
    })
    start_lists[[nm]] <- starts
  }
  start_lists
}
