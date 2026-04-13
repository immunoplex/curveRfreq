# =============================================================================
# StandardCurve.R
#
# R6 class providing a simplified, chainable interface to the curveRfreq
# immunoassay standard-curve fitting workflow.
#
# Data loading is intentionally outside this class. Pull your data first,
# then hand it in via $new() or $set_data().
#
# Dependencies: R6, curveRfreq (devtools::load_all() or installed)
# =============================================================================

StandardCurve <- R6Class(
  classname = "StandardCurve",

  # ---------------------------------------------------------------------------
  # Public fields
  # ---------------------------------------------------------------------------
  public = list(

    # ── Input data ────────────────────────────────────────────────────────────

    #' @field loaded_data List of pre-loaded assay data (standards, blanks,
    #'   samples, antigen_constraints, response_var, indep_var, ...).
    #'   Produced externally, e.g. by `pull_data()`.
    loaded_data = NULL,

    #' @field response_var Name of the response column (e.g. `"mfi"`).
    response_var = NULL,

    #' @field indep_var Name of the independent variable column (e.g. `"concentration"`).
    indep_var = NULL,

    # ── Config ────────────────────────────────────────────────────────────────

    #' @field study_params Named list of study-level modelling parameters.
    #'   Keys: `applyProzone`, `blank_option`, `is_log_response`,
    #'   `is_log_independent`.
    study_params = NULL,

    #' @field antigen_constraints Data frame with per-antigen constraint settings.
    antigen_constraints = NULL,

    #' @field model_names Character vector of candidate model names to consider.
    model_names = NULL,

    #' @field curve_col Name of the curve-ID column in the data. Default `"curve_id"`.
    curve_col = NULL,

    #' @field curve_id_element_order Character vector defining the positional
    #'   order of fields inside a `curve_id` string.
    curve_id_element_order = NULL,

    #' @field is_display_log_response Logical. Display response on log scale in plots.
    is_display_log_response = NULL,

    #' @field is_display_log_independent Logical. Display concentration on log scale in plots.
    is_display_log_independent = NULL,

    #' @field verbose Logical. Print diagnostic messages when `TRUE`.
    verbose = NULL,

    # ── Intermediate results (populated step-by-step) ─────────────────────────

    #' @field antigen_plate List returned by `select_antigen_plate()`.
    antigen_plate = NULL,

    #' @field processed_data List returned by `preprocess_robust_curves()`.
    processed_data = NULL,

    #' @field formulas List of model formulas from `select_model_formulas()`.
    formulas = NULL,

    #' @field model_constraints Parameter bounds from `obtain_model_constraints()`.
    model_constraints = NULL,

    #' @field start_lists Multi-start parameter lists from `make_start_lists()`.
    start_lists = NULL,

    #' @field fit_robust_lm Raw list of candidate model fits.
    fit_robust_lm = NULL,

    #' @field fit_summary Data frame summarising each candidate fit (AIC, BIC, ...).
    fit_summary = NULL,

    #' @field fit_params Data frame of per-model parameter estimates and CIs.
    fit_params = NULL,

    #' @field plot_data List used internally for plotting / model comparison.
    plot_data = NULL,

    #' @field best_fit List returned by `select_model_fit_AIC()` and augmented
    #'   by `tidy.nlsLM()` and `summarize_fit()`.
    best_fit = NULL,

    #' @field se_antigen_table Standard-error lookup table across all antigens.
    se_antigen_table = NULL,

    #' @field current_se SE values for the currently selected antigen.
    current_se = NULL,

    #' @field .selection_args Internal cache of arguments from the last `$select()`.
    .selection_args = NULL,

    # =========================================================================
    # initialize()
    # =========================================================================

    #' @description
    #' Create a new `StandardCurve` object.
    #'
    #' @param loaded_data List of pre-loaded assay data. Must contain at minimum:
    #'   `$standards`, `$blanks`, `$samples`, `$response_var`, `$indep_var`.
    #'   Optionally `$antigen_constraints` (used if `antigen_constraints` arg
    #'   is `NULL`).
    #' @param study_params Named list with keys: `applyProzone`, `blank_option`,
    #'   `is_log_response`, `is_log_independent`.
    #' @param antigen_constraints Data frame of per-antigen constraints. If
    #'   `NULL`, falls back to `loaded_data$antigen_constraints`.
    #' @param model_names Character vector of candidate model names. Defaults
    #'   to all five supported models.
    #' @param curve_col Name of the curve-ID column. Default `"curve_id"`.
    #' @param curve_id_element_order Character vector defining field order inside
    #'   `curve_id`. Defaults to the standard nine-element schema.
    #' @param is_display_log_response Logical. Default `TRUE`.
    #' @param is_display_log_independent Logical. Default `TRUE`.
    #' @param verbose Logical. Default `TRUE`.
    #'
    #' @return A `StandardCurve` object (invisibly).
    #'
    #' @examples
    #' \dontrun{
    #' loaded_data <- pull_data(
    #'   study_accession      = "MADI_01",
    #'   experiment_accession = "IgG1",
    #'   project_id           = 17,
    #'   conn                 = conn
    #' )
    #'
    #' sc <- StandardCurve$new(
    #'   loaded_data         = loaded_data,
    #'   study_params        = list(
    #'     applyProzone       = TRUE,
    #'     blank_option       = "ignored",
    #'     is_log_response    = TRUE,
    #'     is_log_independent = TRUE
    #'   ),
    #'   antigen_constraints = my_constraints_df
    #' )
    #' }
    initialize = function(
    loaded_data,
    study_params = list(
      applyProzone       = TRUE,
      blank_option       = "ignored",
      is_log_response    = TRUE,
      is_log_independent = TRUE
    ),
    antigen_constraints      = NULL,
    model_names              = NULL,
    curve_col                = "curve_id",
    curve_id_element_order   = c("project_id", "study_accession",
                                 "experiment_accession", "feature",
                                 "source", "antigen", "plate",
                                 "nominal_sample_dilution", "wavelength"),
    is_display_log_response    = TRUE,
    is_display_log_independent = TRUE,
    verbose                    = TRUE
    ) {
      self$loaded_data                <- loaded_data
      self$response_var               <- loaded_data$response_var
      self$indep_var                  <- loaded_data$indep_var
      self$study_params               <- study_params
      self$model_names                <- model_names
      self$curve_col                  <- curve_col
      self$curve_id_element_order     <- curve_id_element_order
      self$is_display_log_response    <- is_display_log_response
      self$is_display_log_independent <- is_display_log_independent
      self$verbose                    <- verbose

      # Resolve constraints: explicit arg > loaded_data > error
      self$antigen_constraints <- antigen_constraints %||%
        loaded_data$antigen_constraints %||%
        stop("[StandardCurve] antigen_constraints not found. ",
             "Supply via the antigen_constraints argument or include in loaded_data.")

      private$.check_study_params()

      if (self$verbose) {
        message("[StandardCurve] Initialized.")
        message("  response variable   : ", self$response_var)
        message("  independent variable: ", self$indep_var)
      }

      invisible(self)
    },

    # =========================================================================
    # set_data()
    # =========================================================================

    #' @description
    #' Replace the loaded data (and optionally constraints) without creating a
    #' new object. Resets all downstream results automatically.
    #'
    #' Useful when iterating over multiple studies with a single `StandardCurve`
    #' instance.
    #'
    #' @param loaded_data New pre-loaded data list.
    #' @param antigen_constraints Optional replacement constraints data frame.
    #'   Falls back to `loaded_data$antigen_constraints` or the existing value.
    #'
    #' @return The `StandardCurve` object (invisibly).
    #'
    #' @examples
    #' \dontrun{
    #' loaded_data_2 <- pull_data("MADI_02", "IgG1", 17, conn)
    #' sc$set_data(loaded_data_2)
    #' sc$select(...)$fit()
    #' }
    set_data = function(loaded_data, antigen_constraints = NULL) {
      self$loaded_data         <- loaded_data
      self$response_var        <- loaded_data$response_var
      self$indep_var           <- loaded_data$indep_var
      self$antigen_constraints <- antigen_constraints %||%
        loaded_data$antigen_constraints %||%
        self$antigen_constraints

      private$.reset_downstream()

      if (self$verbose)
        message("[StandardCurve] Data replaced. All downstream results cleared.")

      invisible(self)
    },

    # =========================================================================
    # select()
    # =========================================================================

    #' @description
    #' Filter the loaded data to a single antigen / plate combination and
    #' apply antigen-specific constraints.
    #'
    #' Wraps `select_antigen_plate()`. Results are stored in `$antigen_plate`.
    #'
    #' @param project_id Character or integer.
    #' @param study_accession Character.
    #' @param experiment_accession Character.
    #' @param feature Character. E.g. `"IgG1"`.
    #' @param source Character. E.g. `"Standard"`.
    #' @param antigen Character. Antigen label.
    #' @param plate Character. Plate label, e.g. `"plate_3"`.
    #' @param nominal_sample_dilution Character or numeric.
    #' @param wavelength Character. Defaults to `"__none__"` (bead arrays).
    #'
    #' @return The `StandardCurve` object (invisibly).
    #'
    #' @examples
    #' \dontrun{
    #' sc$select(
    #'   project_id              = 17,
    #'   study_accession         = "MADI_01",
    #'   experiment_accession    = "IgG1",
    #'   feature                 = "IgG1",
    #'   source                  = "Standard",
    #'   antigen                 = "victoria",
    #'   plate                   = "plate_3",
    #'   nominal_sample_dilution = "1000"
    #' )
    #' }
    select = function(
    project_id,
    study_accession,
    experiment_accession,
    feature,
    source,
    antigen,
    plate,
    nominal_sample_dilution,
    wavelength = WL_NONE
    ) {
      private$.step_banner("SELECT ANTIGEN PLATE")

      # Cache for propagate_error() and print()
      self$.selection_args <- list(
        project_id              = project_id,
        study_accession         = study_accession,
        experiment_accession    = experiment_accession,
        feature                 = feature,
        source                  = source,
        antigen                 = antigen,
        plate                   = plate,
        nominal_sample_dilution = nominal_sample_dilution,
        wavelength              = wavelength
      )

      self$antigen_plate <- select_antigen_plate(
        loaded_data             = self$loaded_data,
        project_id              = project_id,
        study_accession         = study_accession,
        experiment_accession    = experiment_accession,
        feature                 = feature,
        source                  = source,
        antigen                 = antigen,
        plate                   = plate,
        nominal_sample_dilution = nominal_sample_dilution,
        wavelength              = wavelength,
        antigen_constraints     = self$antigen_constraints,
        curve_id_element_order  = self$curve_id_element_order,
        verbose                 = self$verbose
      )

      # Clear stale fit results so they are never silently returned
      private$.reset_fit()

      invisible(self)
    },

    # =========================================================================
    # fit()
    # =========================================================================

    #' @description
    #' Run the complete model-fitting pipeline on the selected antigen plate.
    #'
    #' Sequentially calls:
    #' 1. `preprocess_robust_curves()`
    #' 2. `select_model_formulas()`
    #' 3. `obtain_model_constraints()`
    #' 4. `make_start_lists()`
    #' 5. `compute_robust_curves()`
    #' 6. `summarize_model_fits()`
    #' 7. `summarize_model_parameters()`
    #' 8. `get_plot_data()`
    #' 9. `select_model_fit_AIC()`
    #' 10. `tidy.nlsLM()`
    #' 11. `summarize_fit()`
    #'
    #' @param start_quantiles Named numeric vector `c(low, mid, high)`.
    #'   Default `c(low = 0.2, mid = 0.5, high = 0.8)`.
    #'
    #' @return The `StandardCurve` object (invisibly).
    #'
    #' @examples
    #' \dontrun{
    #' sc$fit()
    #' }
    fit = function(start_quantiles = c(low = 0.2, mid = 0.5, high = 0.8)) {
      private$.require_selected()
      private$.step_banner("FIT MODELS")

      ap <- self$antigen_plate
      sp <- self$study_params
      rv <- self$response_var
      iv <- self$indep_var

      # 1 ── Preprocess ───────────────────────────────────────────────────────
      if (self$verbose) message("[1/8] Preprocessing data ...")
      self$processed_data <- preprocess_robust_curves(
        data                 = ap$plate_standard,
        antigen_settings     = ap$antigen_settings,
        response_variable    = rv,
        independent_variable = iv,
        is_log_response      = sp$is_log_response,
        blank_data           = ap$plate_blanks,
        blank_option         = sp$blank_option,
        is_log_independent   = sp$is_log_independent,
        apply_prozone        = sp$applyProzone,
        verbose              = self$verbose
      )

      pd <- self$processed_data$data

      # 2 ── Model formulas ───────────────────────────────────────────────────
      if (self$verbose) message("[2/8] Building model formulas ...")
      self$formulas <- select_model_formulas(
        fixed_constraint  = ap$fixed_a_result,
        response_variable = rv,
        is_log_response   = sp$is_log_response,
        model_names       = self$model_names
      )

      # 3 ── Parameter constraints ────────────────────────────────────────────
      if (self$verbose) message("[3/8] Computing parameter constraints ...")
      self$model_constraints <- obtain_model_constraints(
        data                 = pd,
        formulas             = self$formulas,
        independent_variable = iv,
        response_variable    = rv,
        is_log_response      = sp$is_log_response,
        is_log_concentration = sp$is_log_independent,
        antigen_settings     = ap$antigen_settings,
        max_response         = max(pd[[rv]], na.rm = TRUE),
        min_response         = min(pd[[rv]], na.rm = TRUE),
        verbose              = self$verbose
      )

      # 4 ── Starting values ──────────────────────────────────────────────────
      if (self$verbose) message("[4/8] Generating start lists ...")
      self$start_lists <- make_start_lists(
        model_constraints = self$model_constraints,
        quants            = start_quantiles
      )

      # 5 ── Robust curve fitting ─────────────────────────────────────────────
      if (self$verbose) message("[5/8] Fitting candidate models ...")
      self$fit_robust_lm <- compute_robust_curves(
        prepped_data         = pd,
        response_variable    = rv,
        independent_variable = iv,
        formulas             = self$formulas,
        model_constraints    = self$model_constraints,
        start_lists          = self$start_lists,
        verbose              = self$verbose
      )

      # 6 ── Fit summary ──────────────────────────────────────────────────────
      if (self$verbose) message("[6/8] Summarising model fits ...")
      self$fit_summary <- summarize_model_fits(
        self$fit_robust_lm,
        verbose = self$verbose
      )

      # 7 ── Parameter estimates ──────────────────────────────────────────────
      if (self$verbose) message("[7/8] Extracting parameter estimates ...")
      self$fit_params <- summarize_model_parameters(
        models_fit_list = self$fit_robust_lm,
        level           = 0.95,
        model_names     = self$model_names
      )

      # 8 ── Plot data + AIC selection ────────────────────────────────────────
      if (self$verbose) message("[8/8] Building plot data & selecting best model ...")
      self$plot_data <- get_plot_data(
        models_fit_list        = self$fit_robust_lm,
        prepped_data           = pd,
        fit_params             = self$fit_params,
        fixed_a_result         = ap$fixed_a_result,
        curve_id_element_order = self$curve_id_element_order,
        model_names            = self$model_names,
        x_var                  = iv,
        y_var                  = rv
      )

      self$best_fit <- select_model_fit_AIC(
        fit_summary   = self$fit_summary,
        fit_robust_lm = self$fit_robust_lm,
        fit_params    = self$fit_params,
        plot_data     = self$plot_data,
        verbose       = self$verbose
      )

      # # 9 ── Tidy parameter table ─────────────────────────────────────────────
      # self$best_fit <- tidy.nlsLM(
      #   best_fit            = self$best_fit,
      #   fixed_a_result      = ap$fixed_a_result,
      #   model_constraints   = self$model_constraints,
      #   antigen_settings    = ap$antigen_settings,
      #   antigen_fit_options = self$processed_data$antigen_fit_options,
      #   verbose             = self$verbose
      # )
      #
      # # 10 ── QC summary ──────────────────────────────────────────────────────
      # self$best_fit <- summarize_fit(
      #   best_fit             = self$best_fit,
      #   response_variable    = rv,
      #   independent_variable = iv,
      #   fixed_a_result       = ap$fixed_a_result,
      #   antigen_settings     = ap$antigen_settings,
      #   antigen_fit_options  = self$processed_data$antigen_fit_options
      # )

      if (self$verbose) {
        message("\nBest model: ", self$best_fit$best_model_name %||% "unknown")
        message("Done. Call $summarize(), $plot(), or $propagate_error().")
      }

      invisible(self)
    },

    # =========================================================================
    # propagate_error()
    # =========================================================================

    #' @description
    #' Compute the SE lookup table and propagate measurement error from the
    #' standard curve to the precision profile.
    #'
    #' Calls `compute_antigen_se_table()`, `lookup_antigen_se()`, and
    #' `predict_and_propagate_error()`. Results stored in `$se_antigen_table`,
    #' `$current_se`, and `$best_fit$sample_se`.
    #'
    #' @param grouping_cols Character vector of columns for SE grouping.
    #'
    #' @return The `StandardCurve` object (invisibly).
    #'
    #' @examples
    #' \dontrun{
    #' sc$propagate_error()
    #' head(sc$best_fit$sample_se)
    #' }
    propagate_error = function(
    grouping_cols = c("project_id", "study_accession",
                      "experiment_accession", "source",
                      "antigen", "feature")
    ) {
      private$.require_fitted()
      private$.step_banner("PROPAGATE ERROR")

      args <- self$.selection_args

      self$se_antigen_table <- compute_antigen_se_table(
        standards_data         = self$loaded_data$standards,
        curve_id_element_order = self$curve_id_element_order,
        curve_col              = self$curve_col,
        response_col           = self$response_var,
        dilution_col           = "dilution",
        plate_col              = "plate",
        grouping_cols          = grouping_cols,
        verbose                = self$verbose
      )

      self$current_se <- lookup_antigen_se(
        se_table             = self$se_antigen_table,
        study_accession      = args$study_accession,
        experiment_accession = args$experiment_accession,
        source               = args$source,
        antigen              = args$antigen,
        feature              = args$feature
      )

      self$best_fit <- predict_and_propagate_error(
        best_fit        = self$best_fit,
        response_var    = self$response_var,
        antigen_plate   = self$antigen_plate,
        study_params    = self$study_params,
        se_std_response = self$current_se,
        verbose         = self$verbose
      )

      invisible(self)
    },

    # =========================================================================
    # summarize()
    # =========================================================================

    #' @description
    #' Print a human-readable summary of the fitted standard curve including
    #' QC metrics, goodness-of-fit statistics, and parameter estimates.
    #'
    #' @param digits Integer. Significant digits to display. Default `4`.
    #'
    #' @return The `StandardCurve` object (invisibly).
    #'
    #' @examples
    #' \dontrun{
    #' sc$summarize()
    #' }
    summarize = function(digits = 4) {
      private$.require_fitted()

      bf  <- self$best_fit
      bfs <- bf$best_fit_summary


      ap <- self$antigen_plate
      sp <- self$study_params
      rv <- self$response_var
      iv <- self$indep_var

      cat("\n================================================================\n")
      cat("  Standard Curve Summary\n")
      cat("================================================================\n")
      cat(sprintf("  Best model   : %s\n", bf$best_model_name %||% "-"))
      cat(sprintf("  curve_id     : %s\n",
                  tryCatch(unique(bf$best_data[[self$curve_col]])[1],
                           error = function(e) "-")))
      cat(sprintf("  Response var : %s  |  Independent var: %s\n",
                  self$response_var, self$indep_var))


      # 9 ── Tidy parameter table ─────────────────────────────────────────────
      self$best_fit <- tidy.nlsLM(
        best_fit            = self$best_fit,
        fixed_a_result      = ap$fixed_a_result,
        model_constraints   = self$model_constraints,
        antigen_settings    = ap$antigen_settings,
        antigen_fit_options = self$processed_data$antigen_fit_options,
        verbose             = self$verbose
      )

      # 10 ── QC summary ──────────────────────────────────────────────────────
      self$best_fit <- summarize_fit(
        best_fit             = self$best_fit,
        response_variable    = rv,
        independent_variable = iv,
        fixed_a_result       = ap$fixed_a_result,
        antigen_settings     = ap$antigen_settings,
        antigen_fit_options  = self$processed_data$antigen_fit_options
      )
      # if (!is.null(bf$best_tidy)) {
      #   message("\n─── Parameter Estimates ──────────────────────────────────────")
      #   print(bf$best_tidy)
      # }
      #
      # if (!is.null(bf$best_fit_summary)) {
      #   message("\n─── Fit Statistics & Curve Characteristics ───────────────────")
      #   print(bf$best_fit_summary)
      # }
      # cat("\n  -- Fit Statistics ------------------------------------------\n")
      # if (!is.null(bfs)) {
      #   for (col in intersect(c("rsquare_fit", "aic", "bic", "mse", "cv",
      #                           "nobs", "dfresidual"), names(bfs))) {
      #     val <- bfs[[col]][1]
      #     cat(sprintf("  %-20s: %s\n", col,
      #                 if (is.numeric(val)) signif(val, digits) else val))
      #   }
      # }
      #
      # cat("\n  -- Curve Characteristics ------------------------------------\n")
      # if (!is.null(bfs)) {
      #   for (col in intersect(
      #     c("inflect_x", "inflect_y",
      #       "llod", "ulod", "mindc", "maxdc",
      #       "minrdl", "maxrdl", "lloq", "uloq"),
      #     names(bfs)
      #   )) {
      #     val <- bfs[[col]][1]
      #     cat(sprintf("  %-20s: %s\n", col,
      #                 if (is.numeric(val)) signif(val, digits) else val))
      #   }
      # }
      #
      # cat("\n  -- Model Parameters -----------------------------------------\n")
      # if (!is.null(bf$best_tidy)) {
      #   print(
      #     bf$best_tidy[, intersect(
      #       c("term", "estimate", "std_error", "p_value", "lower", "upper"),
      #       names(bf$best_tidy)
      #     ), drop = FALSE],
      #     digits    = digits,
      #     row.names = FALSE
      #   )
      # }

      cat("================================================================\n\n")
      invisible(self)
    },

    # =========================================================================
    # plot()
    # =========================================================================

    #' @description
    #' Plot the fitted standard curve with QC metric annotations.
    #'
    #' Wraps `plot_standard_curve()`.
    #'
    #' @param is_display_log_independent Logical. Log-scale x-axis. Defaults to the self of the display independent variable axis parameter.
    #' @param is_display_log_response Logical. Log-scale y-axis. Defaults to the self of the display response parameter from constructor.
    #'
    #' @return The ggplot object (invisibly).
    #'
    #' @examples
    #' \dontrun{
    #' sc$plot()
    #' }
    plot = function(is_display_log_independent = self$is_display_log_independent,
                    is_display_log_response = self$is_display_log_response) {
      private$.require_fitted()

      pct <- tryCatch(self$antigen_constraints$pcov_threshold[1],
                      error = function(e) 15)

      p <- plot_standard_curve(
        best_fit                   = self$best_fit,
        is_display_log_independent = is_display_log_independent,
        is_display_log_response    = is_display_log_response,
        pcov_threshold             = pct,
        study_params               = self$study_params,
        curve_id_element_order     = self$curve_id_element_order,
        curve_col                  = self$curve_col,
        response_variable          = self$response_var,
        independent_variable       = self$indep_var
      )

      p
    },

    # =========================================================================
    # compare_models()
    # =========================================================================

    #' @description
    #' Multi-panel plot comparing all candidate model fits, residuals,
    #' parameter estimates, and AIC scores.
    #'
    #' Wraps `plot_model_comparisons()`.
    #'
    #' @param use_patchwork Logical. Use patchwork layout. Default `TRUE`.
    #'
    #' @return The plot object (invisibly).
    #'
    #' @examples
    #' \dontrun{
    #' sc$compare_models()
    #' }
    compare_models = function(use_patchwork = TRUE) {
      private$.require_fitted()

      p <- plot_model_comparisons(
        plot_data     = self$plot_data,
        model_names   = self$model_names,
        x_var         = self$indep_var,
        y_var         = self$response_var,
        use_patchwork = use_patchwork
      )

      print(p)
      invisible(p)
    },

    # =========================================================================
    # get_results()
    # =========================================================================

    #' @description
    #' Return key result tables as a named list for downstream use or export.
    #'
    #' @return Named list: `fit_summary`, `best_tidy`, `best_fit_summary`,
    #'   `sample_se` (populated only after `$propagate_error()`).
    #'
    #' @examples
    #' \dontrun{
    #' res <- sc$get_results()
    #' write.csv(res$best_fit_summary, "summary.csv", row.names = FALSE)
    #' }
    get_results = function() {
      private$.require_fitted()

      list(
        fit_summary      = self$fit_summary,
        best_tidy        = self$best_fit$best_tidy,
        best_fit_summary = self$best_fit$best_fit_summary,
        sample_se        = self$best_fit$sample_se %||% NULL
      )
    },

    # =========================================================================
    # print()
    # =========================================================================

    #' @description Print a concise pipeline-status overview.
    #' @param ... Ignored. Required for S3 print compatibility.
    print = function(...) {
      status <- private$.pipeline_status()

      cat("\n<StandardCurve>\n")
      cat(sprintf("  Stage    : %s\n", status$stage))
      if (!is.null(status$model))
        cat(sprintf("  Model    : %s\n", status$model))
      if (!is.null(status$curve_id))
        cat(sprintf("  curve_id : %s\n", status$curve_id))

      cat("\n  Steps:\n")
      for (nm in names(status$steps)) {
        tick <- if (status$steps[[nm]]) "[v]" else "[ ]"
        cat(sprintf("    %s %s\n", tick, nm))
      }
      cat("\n")
      invisible(self)
    }
  ),

  # ---------------------------------------------------------------------------
  # Private helpers
  # ---------------------------------------------------------------------------
  private = list(

    .step_banner = function(title) {
      if (self$verbose) {
        bar <- paste(rep("-", 55), collapse = "")
        pad <- strrep(" ", max(0L, floor((55L - nchar(title) - 2L) / 2L)))
        message("\n", bar, "\n", pad, " ", title, "\n", bar)
      }
    },

    .check_study_params = function() {
      required <- c("applyProzone", "blank_option",
                    "is_log_response", "is_log_independent")
      missing  <- setdiff(required, names(self$study_params))
      if (length(missing) > 0)
        stop("[StandardCurve] study_params missing keys: ",
             paste(missing, collapse = ", "))
    },

    .require_selected = function() {
      if (is.null(self$antigen_plate))
        stop("[StandardCurve] Call $select() before $fit().")
    },

    .require_fitted = function() {
      if (is.null(self$best_fit))
        stop("[StandardCurve] Call $fit() before this method.")
    },

    # Clears everything downstream of select()
    .reset_fit = function() {
      self$processed_data    <- NULL
      self$formulas          <- NULL
      self$model_constraints <- NULL
      self$start_lists       <- NULL
      self$fit_robust_lm     <- NULL
      self$fit_summary       <- NULL
      self$fit_params        <- NULL
      self$plot_data         <- NULL
      self$best_fit          <- NULL
      self$se_antigen_table  <- NULL
      self$current_se        <- NULL
    },

    # Clears everything downstream of $new() / set_data()
    .reset_downstream = function() {
      self$antigen_plate   <- NULL
      self$.selection_args <- NULL
      private$.reset_fit()
    },

    .pipeline_status = function() {
      steps <- list(
        "select()"          = !is.null(self$antigen_plate),
        "fit()"             = !is.null(self$best_fit),
        "propagate_error()" = !is.null(self$current_se)
      )
      last_done <- rev(names(which(unlist(steps))))[1]
      stage     <- last_done %||% "initialized"

      list(
        stage    = stage,
        steps    = steps,
        model    = self$best_fit$best_model_name %||% NULL,
        curve_id = tryCatch(
          unique(self$best_fit$best_data[[self$curve_col]])[1],
          error = function(e) NULL
        )
      )
    }
  )
)

# -----------------------------------------------------------------------------
# Export — placed here so roxygen2 sees the object, not R6Class() arguments
# -----------------------------------------------------------------------------

#' @title StandardCurve
#' @name StandardCurve
#' @description
#' Wraps the full \code{curveRfreq} standard-curve fitting pipeline into a
#' single stateful object. Data loading is external — supply a pre-loaded
#' data list via \code{$new()} or \code{$set_data()}, then call
#' \code{$select()}, \code{$fit()}, \code{$summarize()}, \code{$plot()}, and
#' \code{$propagate_error()} in sequence. All intermediate results are stored
#' as public fields and can be inspected at any stage.
#' @usage NULL
#' @format NULL
#' @export
StandardCurve
