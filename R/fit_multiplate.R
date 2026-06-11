# =============================================================================
# fit_multiplate.R — Multi-curve loop wrapper
#
# Accepts ALREADY-PREPROCESSED stacked data frames (standards, blanks, and
# optionally samples) with a curve_id column. Splits by curve_id, calls
# fit_calibration_freq() for each, collects results.
#
# No preprocessing. No lookup table. No antigen/plate metadata.
# Just curve_ids.
# =============================================================================


#' Fit Frequentist Calibration Curves Across Multiple Curves
#'
#' Splits preprocessed `standards`, `blanks`, and optionally `samples` by
#' `curve_id`, calls [fit_calibration_freq()] for each curve, and
#' collects the results.
#'
#' **Important:** `standards` and `blanks` must already be on the fitting
#' scale. Use [curveRcore::preprocess_standards()] upstream on each curve's
#' data before stacking, or preprocess the full stacked frame if all
#' curves share the same settings.
#'
#' @param standards Data frame. Stacked preprocessed standard curve data.
#'   Must contain a `curve_id` column, a response column (named by
#'   `response_var`), and a `concentration` column — all already on the
#'   fitting scale.
#' @param blanks Data frame or NULL. Stacked preprocessed blank data. When
#'   non-NULL, must contain a `curve_id` column and a response column (named
#'   by `response_var`), both already on the fitting scale. Stored in each
#'   `result$blanks` for QA and plotting only — not used in fitting.
#'   NULL (default) leaves `result$blanks` empty for every curve.
#' @param samples Data frame or NULL. Stacked sample data with a
#'   `curve_id` column and the response column (raw measurement scale).
#' @param response_var Character. Name of the response column.
#' @param model_names Character vector. Models to fit.
#'   Default `c("logistic4", "gompertz4")`.
#' @param is_log_response Logical. Was the response log10-transformed?
#'   Default TRUE.
#' @param is_log_independent Logical. Was concentration log10-transformed?
#'   Default TRUE.
#' @param std_curve_conc Numeric. Undiluted standard concentration (raw
#'   scale). Used for grid generation.
#' @param fixed_a Numeric or NULL. Fixed lower asymptote on the
#'   **fitting** scale.
#' @param cv_x_max Numeric. Default 150.
#' @param n_grid Integer. Default 200.
#' @param grid_min_conc Numeric. Default 1e-4.
#' @param grid_max_conc Numeric or NULL.
#' @param curve_ids Optional vector. If supplied, only these curve_ids
#'   are fitted. Default NULL fits all found in `standards`.
#' @param on_error Character. `"warn"` (default) stores NULL and continues;
#'   `"stop"` raises on first failure.
#' @param verbose Logical.
#'
#' @return A `calibration_result_multiplate` object (from curveRcore).
#'
#' @export
fit_calibration_freq_multiplate <- function(standards,
                                            blanks = NULL,
                                            samples = NULL,
                                            response_var,
                                            model_names = c("logistic4", "gompertz4"),
                                            is_log_response = TRUE,
                                            is_log_independent = TRUE,
                                            std_curve_conc,
                                            fixed_a = NULL,
                                            cv_x_max = 150,
                                            n_grid = 200L,
                                            grid_min_conc = 1e-4,
                                            grid_max_conc = NULL,
                                            curve_ids = NULL,
                                            on_error = c("warn", "stop"),
                                            verbose = FALSE) {

  on_error <- match.arg(on_error)

  # ── Validate inputs ──
  if (!("curve_id" %in% names(standards)))
    stop("standards must contain a 'curve_id' column")
  if (!(response_var %in% names(standards)))
    stop("response_var '", response_var, "' not found in standards")
  if (!("concentration" %in% names(standards)))
    stop("standards must contain a 'concentration' column (preprocessed)")

  if (!is.null(blanks)) {
    if (!is.data.frame(blanks))
      stop("blanks must be a data frame")
    if (!("curve_id" %in% names(blanks)))
      stop("blanks must contain a 'curve_id' column")
    if (!(response_var %in% names(blanks)))
      stop("response_var '", response_var, "' not found in blanks")
  }

  if (!is.null(samples) && !("curve_id" %in% names(samples)))
    stop("samples must contain a 'curve_id' column when provided")

  # ── Discover curve_ids ──
  all_cids <- sort(unique(standards$curve_id))

  if (!is.null(curve_ids)) {
    missing <- setdiff(curve_ids, all_cids)
    if (length(missing) > 0)
      warning("curve_ids not found in standards: ",
              paste(missing, collapse = ", "))
    all_cids <- intersect(all_cids, curve_ids)
  }

  n <- length(all_cids)
  if (n == 0) stop("No curve_ids to fit")

  # ── Split data ──
  std_splits   <- split(standards, standards$curve_id)
  blank_splits <- if (!is.null(blanks)) {
    split(blanks, blanks$curve_id)
  } else {
    list()
  }

  # Warn about any curve_ids that have standards but no blanks
  if (!is.null(blanks)) {
    missing_blank_cids <- setdiff(as.character(all_cids),
                                  names(blank_splits))
    if (length(missing_blank_cids) > 0)
      warning("curve_ids present in standards but not in blanks: ",
              paste(missing_blank_cids, collapse = ", "),
              ". An empty data frame will be used for those curves.")
  }

  samp_splits <- if (!is.null(samples)) {
    split(samples, samples$curve_id)
  } else {
    NULL
  }

  # ── Fit each curve ──
  plates      <- vector("list", n)
  plate_names <- as.character(all_cids)

  for (i in seq_len(n)) {
    cid <- as.character(all_cids[i])

    if (verbose) message(sprintf("[%d/%d] curve_id=%s", i, n, cid))

    this_std   <- std_splits[[cid]]
    this_blank <- if (!is.null(blanks)) {
      if (!is.null(blank_splits[[cid]])) {
        blank_splits[[cid]]
      } else {
        # blanks supplied but this curve_id is absent — use zero-row frame
        blanks[0L, , drop = FALSE]
      }
    } else {
      NULL
    }
    this_samp  <- if (!is.null(samp_splits)) samp_splits[[cid]] else NULL

    if (is.null(this_samp) ||
        (is.data.frame(this_samp) && nrow(this_samp) == 0))
      this_samp <- NULL

    result <- tryCatch(
      withCallingHandlers(
        fit_calibration_freq(
          standards          = this_std,
          blanks             = this_blank,
          samples            = this_samp,
          response_var       = response_var,
          model_names        = model_names,
          is_log_response    = is_log_response,
          is_log_independent = is_log_independent,
          std_curve_conc     = std_curve_conc,
          fixed_a            = fixed_a,
          cv_x_max           = cv_x_max,
          n_grid             = n_grid,
          grid_min_conc      = grid_min_conc,
          grid_max_conc      = grid_max_conc,
          curve_id           = cid,
          verbose            = verbose
        ),
        warning = function(w) {
          if (verbose) message("  [warning] ", conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        msg <- sprintf("curve_id %s failed: %s", cid, conditionMessage(e))
        if (on_error == "stop") stop(msg, call. = FALSE)
        warning(msg, call. = FALSE)
        NULL
      }
    )

    plates[i] <- list(result)
  }

  names(plates) <- plate_names

  # ── Build multiplate meta ──
  multi_meta <- list(
    method             = "frequentist",
    package            = "curveRfreq",
    curve_ids          = all_cids,
    n_curves           = n,
    response_var       = response_var,
    is_log_response    = is_log_response,
    is_log_independent = is_log_independent,
    timestamp          = Sys.time()
  )

  curveRcore::new_calibration_result_multiplate(
    meta   = multi_meta,
    plates = plates
  )
}


# =============================================================================
# Summary extraction
# =============================================================================

#' Extract a Combined Summary Table from Multi-Curve Results
#'
#' One row per curve_id with best model, parameters, and fit statistics.
#'
#' @param mp A `calibration_result_multiplate` object.
#'
#' @return Data frame with columns: `curve_id`, `best_model`, `converged`,
#'   `aic`, `bic`, `rss`, `n_obs`, `n_samples`, `a`, `b`, `c`, `d`, `g`.
#'
#' @export
summary_table <- function(mp) {
  stopifnot(inherits(mp, "calibration_result_multiplate"))

  rows <- lapply(seq_along(mp$plates), function(i) {
    cr  <- mp$plates[[i]]
    cid <- names(mp$plates)[i]

    if (is.null(cr)) {
      return(data.frame(
        curve_id = cid,
        best_model = NA_character_, converged = FALSE,
        aic = NA_real_, bic = NA_real_, rss = NA_real_,
        n_obs = NA_integer_, n_samples = NA_integer_,
        a = NA_real_, b = NA_real_, c = NA_real_,
        d = NA_real_, g = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    best <- cr$selection$best_model_name
    ens  <- if (!is.na(best)) cr$ensemble[[best]] else NULL

    pv <- if (!is.null(ens) && !is.null(ens$parameters)) {
      stats::setNames(ens$parameters$estimate, ens$parameters$term)
    } else NULL

    fs <- if (!is.null(ens)) ens$fit_stats else NULL

    data.frame(
      curve_id   = cid,
      best_model = if (is.na(best)) NA_character_ else best,
      converged  = !is.na(best),
      aic        = if (!is.null(fs)) fs$aic else NA_real_,
      bic        = if (!is.null(fs)) fs$bic else NA_real_,
      rss        = if (!is.null(fs)) fs$rss else NA_real_,
      n_obs      = cr$meta$n_standards,
      n_samples  = cr$meta$n_samples,
      a = if (!is.null(pv) && "a" %in% names(pv)) pv[["a"]] else NA_real_,
      b = if (!is.null(pv) && "b" %in% names(pv)) pv[["b"]] else NA_real_,
      c = if (!is.null(pv) && "c" %in% names(pv)) pv[["c"]] else NA_real_,
      d = if (!is.null(pv) && "d" %in% names(pv)) pv[["d"]] else NA_real_,
      g = if (!is.null(pv) && "g" %in% names(pv)) pv[["g"]] else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}


#' Extract All Sample Predictions from Multi-Curve Results
#'
#' Concatenates sample prediction data frames from every curve, adding
#' a `curve_id` column for identification.
#'
#' @param mp A `calibration_result_multiplate` object.
#'
#' @return Data frame of all sample predictions, or NULL if none.
#'
#' @export
collect_samples <- function(mp) {
  stopifnot(inherits(mp, "calibration_result_multiplate"))

  dfs <- lapply(seq_along(mp$plates), function(i) {
    cr <- mp$plates[[i]]
    if (is.null(cr) || is.null(cr$samples) || nrow(cr$samples) == 0)
      return(NULL)
    s <- cr$samples
    s$curve_id <- names(mp$plates)[i]
    s
  })

  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs) == 0) return(NULL)
  do.call(rbind, dfs)
}
