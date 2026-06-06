# test-fitting.R — Tests for curveRfreq fitting pipeline
#
# Uses bead_assay_example data from curveRcore, preprocessed upstream.
# Tests all five models (where applicable) and the full pipeline.

library(curveRcore)

# =============================================================================
# Test data setup — preprocess once, reuse everywhere
# =============================================================================
data("bead_assay_example", package = "curveRcore")

std_raw <- bead_assay_example$standards[
  bead_assay_example$standards$curve_id %in% c(1, 2, 3), ]

prepped <- curveRcore::preprocess_standards(
  data                 = std_raw,
  antigen_settings     = list(standard_curve_concentration = 10000),
  response_variable    = "mfi",
  independent_variable = "concentration",
  is_log_response      = TRUE,
  is_log_independent   = TRUE,
  apply_prozone        = TRUE
)

standards <- prepped$data

samples <- bead_assay_example$samples[
  bead_assay_example$samples$curve_id %in% c(1, 2, 3), ]

# Single curve for single-plate tests
std_single <- standards[standards$curve_id == 1, ]
samp_single <- samples[samples$curve_id == 1, ]


# =============================================================================
# Single-curve fitting: logistic4 + gompertz4 (default)
# =============================================================================
test_that("fit_calibration_freq works for single curve with defaults", {
  result <- fit_calibration_freq(
    standards      = std_single,
    samples        = samp_single,
    response_var   = "mfi",
    std_curve_conc = 10000,
    curve_id       = "1",
    verbose        = FALSE
  )
  expect_s3_class(result, "calibration_result")
  expect_equal(result$meta$method, "frequentist")
  expect_equal(result$meta$curve_id, "1")
  expect_true(!is.na(result$selection$best_model_name))
  expect_true(nrow(result$grid) > 0)
  expect_true(nrow(result$samples) > 0)
})


# =============================================================================
# Single-curve fitting: all five models
# =============================================================================
test_that("fit_calibration_freq works with all five model names", {
  # On log scale, loglogistic4 is dropped (redundant with logistic4)
  result <- fit_calibration_freq(
    standards      = std_single,
    response_var   = "mfi",
    model_names    = curveRcore::available_models(),
    std_curve_conc = 10000,
    curve_id       = "1"
  )
  # Should have 4 models (loglogistic4 removed on log scale)
  expect_equal(length(result$ensemble), 4)
  converged <- vapply(result$ensemble, function(e) e$converged, logical(1))
  expect_true(sum(converged) >= 1, "At least one model should converge")
})


# =============================================================================
# Parameter sanity: a is estimated freely when fixed_a = NULL
# =============================================================================
test_that("a parameter is estimated freely (not pinned at ~-5)", {
  result <- fit_calibration_freq(
    standards      = std_single,
    response_var   = "mfi",
    std_curve_conc = 10000,
    fixed_a        = NULL,
    curve_id       = "1"
  )
  best <- result$selection$best_model_name
  params <- result$ensemble[[best]]$parameters
  a_est <- params$estimate[params$term == "a"]

  # a should be in a reasonable range on log10(MFI) scale, not near -5
  expect_true(a_est > 0, info = paste("a =", a_est, "should be > 0"))
  expect_true(a_est < 3, info = paste("a =", a_est, "should be < 3"))
})


# =============================================================================
# Parameter sanity: fixed_a constrains a
# =============================================================================
test_that("fixed_a bakes a into formula (not estimated)", {
  fixed_val <- 1.5
  result <- fit_calibration_freq(
    standards      = std_single,
    response_var   = "mfi",
    std_curve_conc = 10000,
    fixed_a        = fixed_val,
    curve_id       = "1"
  )
  best <- result$selection$best_model_name
  params <- result$ensemble[[best]]$parameters
  a_row <- params[params$term == "a", ]
  if (nrow(a_row) > 0) {
    # fixed_a appears in the parameter table with exact value and zero SE
    expect_equal(a_row$estimate, fixed_val, tolerance = 1e-6)
    expect_equal(a_row$std_error, 0)
  } else {
    # fixed_a was baked into the formula — a does not appear as a free param
    # Verify it's not among the estimated terms
    expect_false("a" %in% params$term)
  }
})


# =============================================================================
# Multi-curve fitting
# =============================================================================
test_that("fit_calibration_freq_multiplate works with stacked data", {
  result <- fit_calibration_freq_multiplate(
    standards      = standards,
    samples        = samples,
    response_var   = "mfi",
    std_curve_conc = 10000,
    verbose        = FALSE
  )
  expect_s3_class(result, "calibration_result_multiplate")
  expect_equal(length(result$plates), 3)
  expect_equal(sort(names(result$plates)), c("1", "2", "3"))
})


# =============================================================================
# Summary table
# =============================================================================
test_that("summary_table returns one row per curve_id", {
  result <- fit_calibration_freq_multiplate(
    standards      = standards,
    response_var   = "mfi",
    std_curve_conc = 10000
  )
  tbl <- summary_table(result)
  expect_equal(nrow(tbl), 3)
  expect_true("curve_id" %in% names(tbl))
  expect_true("best_model" %in% names(tbl))
  expect_true("a" %in% names(tbl))
  expect_true(all(tbl$converged))
})


# =============================================================================
# collect_samples
# =============================================================================
test_that("collect_samples concatenates all curve predictions", {
  result <- fit_calibration_freq_multiplate(
    standards      = standards,
    samples        = samples,
    response_var   = "mfi",
    std_curve_conc = 10000
  )
  all_samp <- collect_samples(result)
  expect_true(!is.null(all_samp))
  expect_true("curve_id" %in% names(all_samp))
  expect_true("predicted_log10_concentration" %in% names(all_samp))
})


# =============================================================================
# Grid prediction has expected columns
# =============================================================================
test_that("grid output has all expected columns", {
  result <- fit_calibration_freq(
    standards      = std_single,
    response_var   = "mfi",
    std_curve_conc = 10000,
    curve_id       = "1"
  )
  grid_cols <- c("log10_concentration", "concentration", "x_fit",
                 "predicted_response", "ci_lower", "ci_upper",
                 "predicted_concentration", "se_concentration",
                 "pcov", "pcov_pass")
  for (col in grid_cols) {
    expect_true(col %in% names(result$grid), info = paste("Missing:", col))
  }
})


# =============================================================================
# Curve_ids subset
# =============================================================================
test_that("curve_ids argument filters correctly", {
  result <- fit_calibration_freq_multiplate(
    standards      = standards,
    response_var   = "mfi",
    std_curve_conc = 10000,
    curve_ids      = c(1, 3)
  )
  expect_equal(length(result$plates), 2)
})


# =============================================================================
# Error handling
# =============================================================================
test_that("missing response_var raises error", {
  expect_error(
    fit_calibration_freq(
      standards      = std_single,
      response_var   = "nonexistent",
      std_curve_conc = 10000
    ),
    "not found"
  )
})

test_that("missing concentration column raises error", {
  bad_df <- std_single
  bad_df$concentration <- NULL
  expect_error(
    fit_calibration_freq(
      standards      = bad_df,
      response_var   = "mfi",
      std_curve_conc = 10000
    ),
    "concentration"
  )
})

test_that("missing curve_id in multiplate raises error", {
  bad_df <- standards
  bad_df$curve_id <- NULL
  expect_error(
    fit_calibration_freq_multiplate(
      standards      = bad_df,
      response_var   = "mfi",
      std_curve_conc = 10000
    ),
    "curve_id"
  )
})
