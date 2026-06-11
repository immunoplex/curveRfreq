# =============================================================================
# fit_calibration.R — Single-curve frequentist calibration
#
# Receives ALREADY-PREPROCESSED standards. The response column and the
# concentration column are both on their final fitting scale (log or raw).
# All preprocessing (concentration computation, prozone, blanks,
# log-transform) is done upstream via curveRcore::preprocess_standards()
# or the individual curveRcore transform functions.
#
# curveRfreq does one thing: fit nonlinear calibration models.
#
# Formulas come from curveRcore::build_nls_formulas().
# Identified by curve_id only — no antigen, plate, or feature metadata.
# =============================================================================


#' Fit a Frequentist Calibration Curve (Single Curve)
#'
#' Fits an ensemble of nonlinear models to preprocessed standard curve
#' data for one curve_id, selects the best model by AIC subject to
#' quantification-eligibility gates, and produces grid and sample
#' predictions with propagated uncertainty.
#'
#' **Important:** `standards` and 'blanks' must already be on the fitting scale.
#' Both the response column and the `concentration` column should be
#' transformed (e.g. log10) before calling this function. Use
#' [curveRcore::preprocess_standards()] upstream.
#'
#' @param standards Data frame. Preprocessed standard curve data with a
#'   response column (named by `response_var`) and a `concentration`
#'   column, both already on the fitting scale.
#' @param blanks Data frame or NULL. Preprocessed blank data with a
#'   response column (named by `response_var`), on the fitting scale.
#'   Stored in `result$blanks` for QA and plotting only — not used in
#'   fitting. NULL (default) leaves `result$blanks` empty.
#' @param samples Data frame or NULL. Test sample data with the response
#'   column (on the **raw** measurement scale — log-transform is applied
#'   internally for back-calculation). Optionally a `dilution` column.
#' @param response_var Character. Name of the response column.
#' @param model_names Character vector. Models to fit. Must be a subset
#'   of [curveRcore::available_models()].
#'   Default `c("logistic4", "gompertz4")`.
#' @param is_log_response Logical. Was the response log10-transformed
#'   during preprocessing? Needed to correctly back-transform sample
#'   predictions. Default TRUE.
#' @param is_log_independent Logical. Was concentration log10-transformed?
#'   Needed for pcov computation and grid generation. Default TRUE.
#' @param std_curve_conc Numeric. Undiluted standard concentration (raw
#'   scale). Used for grid generation.
#' @param fixed_a Numeric or NULL. Fixed lower asymptote on the **fitting**
#'   scale (i.e. already log-transformed if `is_log_response = TRUE`).
#'   NULL means `a` is estimated freely.
#' @param cv_x_max Numeric. Cap for percent CV of back-calculated
#'   concentration. Default 150.
#' @param pcov_threshold Numeric. Percent CV threshold for LLOQ/ULOQ
#'   determination and the dynamic-range eligibility gate. Default 20.
#' @param min_dynamic_range_log10 Numeric. Minimum dynamic range (log10
#'   units) required for a model to be eligible for selection. Default 0.5.
#' @param max_rel_se Numeric. Maximum relative SE (SE/|estimate|) permitted
#'   for any parameter; models exceeding this are ineligible. Default 5.0.
#' @param bound_tol Numeric. Tolerance for the at-bound gate. Default 1e-4.
#' @param n_grid Integer. Number of prediction grid points. Default 200.
#' @param grid_min_conc Numeric. Minimum grid concentration (raw scale).
#'   Default 1e-4.
#' @param grid_max_conc Numeric or NULL. Maximum grid concentration.
#'   NULL uses `std_curve_conc`.
#' @param curve_id Character or numeric. Identifier for this curve.
#'   Default `"1"`.
#' @param verbose Logical. Default FALSE.
#'
#' @return A `calibration_result` object (from curveRcore).  The
#'   `$selection` component now contains `$assessments` (per-model
#'   eligibility results), `$eligible_models`, and `$fallback`.
#'
#' @export
fit_calibration_freq <- function(standards,
                                 blanks = NULL,
                                 samples = NULL,
                                 response_var,
                                 model_names = c("logistic4", "gompertz4"),
                                 is_log_response = TRUE,
                                 is_log_independent = TRUE,
                                 std_curve_conc,
                                 fixed_a = NULL,
                                 cv_x_max = 150,
                                 pcov_threshold = 20,
                                 min_dynamic_range_log10 = 0.5,
                                 max_rel_se = 5.0,
                                 bound_tol = 1e-4,
                                 n_grid = 200L,
                                 grid_min_conc = 1e-4,
                                 grid_max_conc = NULL,
                                 curve_id = "1",
                                 verbose = FALSE) {

  indep_var <- "concentration"

  # ── 1. Validate inputs ── (unchanged)
  if (!(response_var %in% names(standards)))
    stop("response_var '", response_var, "' not found in standards")
  if (!(indep_var %in% names(standards)))
    stop("standards must contain a '", indep_var, "' column (preprocessed)")

  # ── 2. Resolve effective models ── (unchanged)
  if (is_log_independent) {
    model_names <- setdiff(model_names, "loglogistic4")
  }
  if (length(model_names) == 0)
    stop("No models to fit after resolving redundant parameterisations")

  # ── 3. Build formulas from curveRcore ── (unchanged)
  formulas <- curveRcore::build_nls_formulas(
    model_names          = model_names,
    response_variable    = response_var,
    independent_variable = indep_var,
    fixed_a              = NULL,
    is_log_response      = is_log_response
  )

  if (!is.null(fixed_a)) {
    formulas <- curveRcore::build_nls_formulas(
      model_names          = model_names,
      response_variable    = response_var,
      independent_variable = indep_var,
      fixed_a              = fixed_a,
      is_log_response      = FALSE
    )
  }

  # ── 4. Compute constraints and start lists ── (unchanged)
  constraints <- compute_model_constraints(
    data                 = standards,
    formulas             = formulas,
    response_variable    = response_var,
    independent_variable = indep_var,
    is_log_response      = is_log_response,
    fixed_a              = fixed_a,
    verbose              = verbose
  )

  starts <- generate_start_lists(constraints)

  # ── 5. Fit ensemble ── (unchanged)
  ensemble_raw <- fit_ensemble_nls(
    data              = standards,
    formulas          = formulas,
    model_constraints = constraints,
    start_lists       = starts,
    verbose           = verbose
  )

  # ── 6. AIC ranking — unchanged, always computed ──
  # select_best_aic() is preserved as a pure AIC ranking step.
  # Its output is retained in selection$weights for full auditability.
  aic_selection <- select_best_aic(ensemble_raw)

  # ── 6a. Per-model precision grids ──────────────────────────────────────────
  # Compute a pcov grid for every converged model. This enables the
  # dynamic-range eligibility gate and per-model precision reporting.
  # Uses the same generate_prediction_grid + predict_grid_freq pipeline
  # that the best-model grid (Section 9) uses.

  base_grid_for_gates <- curveRcore::generate_prediction_grid(
    std_curve_conc     = std_curve_conc,
    n_grid             = n_grid,
    grid_min_conc      = grid_min_conc,
    grid_max_conc      = grid_max_conc,
    is_log_independent = is_log_independent
  )

  # Named list: model_name -> grid data frame (or NULL if not converged)
  ensemble_grids <- lapply(names(ensemble_raw), function(nm) {
    e <- ensemble_raw[[nm]]
    if (!e$converged) return(NULL)
    sigma_fit <- tryCatch(summary(e$fit)$sigma, error = function(e) 0)
    tryCatch(
      predict_grid_freq(
        grid                 = base_grid_for_gates,
        fit                  = e$fit,
        model_name           = nm,
        fixed_a              = fixed_a,
        se_response          = sigma_fit,
        cv_x_max             = cv_x_max,
        is_log_independent = is_log_independent,
        is_log_response = is_log_response,
        independent_variable = indep_var
      ),
      error = function(e) NULL
    )
  })
  names(ensemble_grids) <- names(ensemble_raw)

  # ── 6b. Per-model eligibility assessment ───────────────────────────────────
  # Applies the four gates: at_bound, vcov_condition, rel_se, dynamic_range.
  # The frequentist path supplies constraints and vcov_matrix, so all four
  # gates are active. Results are stored for downstream reporting.

  assessments <- lapply(names(ensemble_raw), function(nm) {
    e <- ensemble_raw[[nm]]

    # Parameters data frame: use summary coefficients (has std_error)
    params_df <- if (e$converged) {
      s  <- summary(e$fit)
      ct <- as.data.frame(s$coefficients)
      data.frame(term      = rownames(ct),
                 estimate  = ct[, "Estimate"],
                 std_error = ct[, "Std. Error"],
                 stringsAsFactors = FALSE)
    } else NULL

    # vcov matrix for condition-number gate
    vcov_mat <- if (e$converged) {
      tryCatch(stats::vcov(e$fit), error = function(e) NULL)
    } else NULL

    # Constraints for this model (lower/upper bounds on free params)
    model_con <- if (nm %in% names(constraints)) constraints[[nm]] else NULL

    # pcov profile from the per-model grid
    g <- ensemble_grids[[nm]]
    pcov_vec <- if (!is.null(g)) g$pcov else NULL
    grid_x   <- if (!is.null(g)) g$log10_concentration else NULL

    if (!e$converged || is.null(params_df)) {
      # Non-converged models: ineligible by definition, no gates run
      return(list(
        model_name          = nm,
        eligible            = FALSE,
        gates               = data.frame(
          gate   = "convergence",
          passed = FALSE,
          detail = "Model did not converge",
          stringsAsFactors = FALSE
        ),
        dynamic_range_log10 = NA_real_,
        lloq                = NA_real_,
        uloq                = NA_real_
      ))
    }

    curveRcore::assess_model_eligibility(
      model_name              = nm,
      parameters              = params_df,
      constraints             = model_con,
      pcov_profile            = pcov_vec,
      grid_x                  = grid_x,
      pcov_threshold          = pcov_threshold,
      bound_tol               = bound_tol,
      max_rel_se              = max_rel_se,
      min_dynamic_range_log10 = min_dynamic_range_log10,
      vcov_matrix             = vcov_mat
    )
  })
  names(assessments) <- names(ensemble_raw)

  # ── 6b+. Attach shape-LOQ to each assessment from the d2y_dx2-enriched grid ─
  for (nm in names(assessments)) {
    g <- ensemble_grids[[nm]]
    if (!is.null(g) && "d2y_dx2" %in% names(g)) {
      shape_loq <- curveRcore::compute_shape_loq_from_grid(g)
      assessments[[nm]] <- c(assessments[[nm]], shape_loq)
    }
  }

  if (verbose) {
    for (nm in names(assessments)) {
      a <- assessments[[nm]]
      status <- if (a$eligible) "\u2713 eligible" else "\u2717 ineligible"
      message(sprintf("  [eligibility] %-14s %s", nm, status))
      if (!a$eligible && nrow(a$gates) > 0) {
        failed <- a$gates[!a$gates$passed, , drop = FALSE]
        for (j in seq_len(nrow(failed))) {
          message(sprintf("    gate %-16s %s", failed$gate[j], failed$detail[j]))
        }
      }
    }
  }

  # ── 6c. Eligible model selection ───────────────────────────────────────────
  # Build the AIC ranking scores (lower is better) for converged models.
  aic_scores <- vapply(names(ensemble_raw), function(nm) {
    e <- ensemble_raw[[nm]]
    if (!e$converged) return(NA_real_)
    tryCatch(stats::AIC(e$fit), error = function(e) NA_real_)
  }, numeric(1))

  eligible_selection <- curveRcore::select_best_eligible(
    assessments      = assessments,
    ranking_scores   = aic_scores,
    criterion        = "AIC+eligibility",
    higher_is_better = FALSE
  )

  # Merge the AIC weights table from select_best_aic() into the eligibility
  # selection so that callers have the full picture in one place.
  eligible_selection$weights   <- aic_selection$weights
  eligible_selection$aic_best  <- aic_selection$best_model_name

  best_name <- eligible_selection$best_model_name

  if (verbose && !is.na(best_name)) {
    fb_tag <- if (isTRUE(eligible_selection$fallback)) " [fallback]" else ""
    message(sprintf("  [selection] best = %s%s", best_name, fb_tag))
  }

  # ── 7. Extract parameters ── (unchanged, now uses eligibility-gated best_name)
  params_df <- extract_best_parameters(
    ensemble_raw, best_name, fixed_a = fixed_a,
    model_constraints = constraints
  )

  # ── 8. Build ensemble output ────────────────────────────────────────────────
  # Attach per-model grids and eligibility assessments to each ensemble entry.
  # The ensemble_grids computed in 6a are reused here — no recomputation.

  ensemble_out <- lapply(names(ensemble_raw), function(nm) {
    e <- ensemble_raw[[nm]]

    if (!e$converged) {
      return(list(
        model_name  = e$model_name,
        converged   = FALSE,
        parameters  = NULL,
        fit_stats   = NULL,
        raw_fit     = NULL,
        grid        = NULL,
        eligibility = assessments[[nm]]
      ))
    }

    fit <- e$fit
    s   <- summary(fit)
    ct  <- as.data.frame(s$coefficients)

    list(
      model_name = e$model_name,
      converged  = TRUE,
      parameters = data.frame(
        term      = rownames(ct),
        estimate  = ct[, 1],
        std_error = ct[, 2],
        stringsAsFactors = FALSE
      ),
      fit_stats = list(
        n_obs    = stats::nobs(fit),
        n_params = length(stats::coef(fit)),
        df_resid = stats::df.residual(fit),
        rss      = sum(stats::residuals(fit)^2),
        aic      = stats::AIC(fit),
        bic      = stats::BIC(fit),
        log_lik  = as.numeric(stats::logLik(fit))
      ),
      raw_fit     = fit,
      grid        = ensemble_grids[[nm]],    # per-model precision grid
      eligibility = assessments[[nm]]        # per-model gate results
    )
  })
  names(ensemble_out) <- names(ensemble_raw)

  # ── 9. Best-model grid (full resolution) ────────────────────────────────────
  # The per-model grids (Section 6a) used the same base_grid_for_gates, so the
  # plate-level $grid is simply taken from the best-model ensemble entry —
  # no recomputation needed.

  grid <- if (!is.na(best_name) && !is.null(ensemble_out[[best_name]]$grid)) {
    ensemble_out[[best_name]]$grid
  } else {
    # Fallback: generate and predict (handles edge cases where best has no grid)
    g <- curveRcore::generate_prediction_grid(
      std_curve_conc     = std_curve_conc,
      n_grid             = n_grid,
      grid_min_conc      = grid_min_conc,
      grid_max_conc      = grid_max_conc,
      is_log_independent = is_log_independent
    )
    if (!is.na(best_name)) {
      best_fit  <- ensemble_raw[[best_name]]$fit
      sigma_fit <- tryCatch(summary(best_fit)$sigma, error = function(e) 0)
      g <- predict_grid_freq(
        grid                 = g,
        fit                  = best_fit,
        model_name           = best_name,
        fixed_a              = fixed_a,
        se_response          = sigma_fit,
        cv_x_max             = cv_x_max,
        is_log_independent            = is_log_independent,
        independent_variable = indep_var,
        verbose              = verbose
      )
    }
    g
  }

  # ── 10. Sample prediction ── (unchanged, uses eligibility-gated best_name)
  samples_out <- NULL
  if (!is.null(samples) && nrow(samples) > 0 && !is.na(best_name)) {
    best_fit  <- ensemble_raw[[best_name]]$fit
    sigma_fit <- tryCatch(summary(best_fit)$sigma, error = function(e) 0)

    samples_out <- predict_samples_freq(
      samples           = samples,
      fit               = best_fit,
      model_name        = best_name,
      response_variable = response_var,
      fixed_a           = fixed_a,
      is_log_response   = is_log_response,
      se_response       = sigma_fit,
      cv_x_max          = cv_x_max,
      is_log_independent = is_log_independent,
      verbose           = verbose
    )
  }

  # ── 11. Build calibration_result ─────────────────────────────────────────
  meta <- list(
    method             = "frequentist",
    package            = "curveRfreq",
    version            = as.character(utils::packageVersion("curveRfreq")),
    timestamp          = Sys.time(),
    curve_id           = as.character(curve_id),
    response_var       = response_var,
    independent_var    = indep_var,
    is_log_response    = is_log_response,
    is_log_independent = is_log_independent,
    n_standards        = nrow(standards),
    n_samples          = if (!is.null(samples)) nrow(samples) else 0L,
    pcov_threshold     = pcov_threshold
  )

  cr <- curveRcore::new_calibration_result(
    meta      = meta,
    ensemble  = ensemble_out,
    selection = eligible_selection,
    grid      = grid,
    samples   = samples_out,
    standards = standards,
    blanks = blanks
  )

  # ── 12. Detection limits (LODs, MDC, RDL) for the best eligible model ────
  cr <- curveRcore::compute_detection_limits(cr, verbose = verbose)

  cr
}
