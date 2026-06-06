# =============================================================================
# test-elisa-fitting.R — Tests using the ELISA example dataset
#
# The ELISA data has OD response (0–4 range), 6 curve_ids, single analyte.
# Plates 1-3 were generated from a 5PL model, plates 4-6 from Gompertz.
# This tests a different response scale (low dynamic range) and exercises
# the adaptive constraint profile's medium/low scale handling.
# =============================================================================

library(curveRcore)


# =============================================================================
# Setup: preprocess ELISA data
# =============================================================================

data("elisa_assay_example", package = "curveRcore")

elisa_std_raw <- elisa_assay_example$standards
elisa_samples <- elisa_assay_example$samples

# Preprocess all 6 plates together
elisa_prepped <- curveRcore::preprocess_standards(
  data                 = elisa_std_raw,
  antigen_settings     = list(standard_curve_concentration = 10000),
  response_variable    = "od",
  independent_variable = "concentration",
  is_log_response      = TRUE,
  is_log_independent   = TRUE,
  apply_prozone        = TRUE
)

elisa_standards <- elisa_prepped$data


# =============================================================================
# Single curve: fit curve_id=1 (5PL-generated data)
# =============================================================================

test_that("ELISA curve_id=1 fits with all applicable models", {
  std1 <- elisa_standards[elisa_standards$curve_id == 1, ]

  result <- fit_calibration_freq(
    standards      = std1,
    response_var   = "od",
    model_names    = curveRcore::available_models(),
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    curve_id       = "1",
    verbose        = FALSE
  )

  expect_s3_class(result, "calibration_result")
  expect_equal(result$meta$curve_id, "1")

  # At least one model should converge
  converged <- vapply(result$ensemble, function(e) e$converged, logical(1))
  expect_true(sum(converged) >= 1,
              info = paste("Converged:", paste(names(converged[converged]),
                                               collapse = ", ")))
})


# =============================================================================
# Single curve: fit curve_id=4 (Gompertz-generated data)
# =============================================================================

test_that("ELISA curve_id=4 (Gompertz-generated) fits correctly", {
  std4 <- elisa_standards[elisa_standards$curve_id == 4, ]

  result <- fit_calibration_freq(
    standards      = std4,
    response_var   = "od",
    model_names    = c("logistic4", "gompertz4"),
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    curve_id       = "4",
    verbose        = FALSE
  )

  expect_s3_class(result, "calibration_result")
  converged <- vapply(result$ensemble, function(e) e$converged, logical(1))
  expect_true(sum(converged) >= 1)

  # Parameter a should be reasonable on log10(OD) scale
  best <- result$selection$best_model_name
  if (!is.na(best)) {
    params <- result$ensemble[[best]]$parameters
    a_val <- params$estimate[params$term == "a"]
    if (length(a_val) > 0) {
      # log10(OD) lower asymptote: OD ~ 0.01–0.1 → log10 ~ -2 to -1
      expect_true(a_val > -3, info = paste("a =", a_val))
      expect_true(a_val < 1, info = paste("a =", a_val))
    }
  }
})


# =============================================================================
# Multi-curve: all 6 ELISA plates
# =============================================================================

test_that("ELISA multiplate fits all 6 curve_ids", {
  result <- fit_calibration_freq_multiplate(
    standards      = elisa_standards,
    samples        = elisa_samples,
    response_var   = "od",
    model_names    = c("logistic4", "gompertz4"),
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    verbose        = FALSE
  )

  expect_s3_class(result, "calibration_result_multiplate")
  expect_equal(length(result$plates), 6)

  tbl <- summary_table(result)
  expect_equal(nrow(tbl), 6)

  # Most plates should converge
  expect_true(sum(tbl$converged) >= 4,
              info = paste("Converged:", sum(tbl$converged), "of 6"))
})


# =============================================================================
# Multi-curve: subset of plates
# =============================================================================

test_that("ELISA multiplate works with curve_ids subset", {
  result <- fit_calibration_freq_multiplate(
    standards      = elisa_standards,
    response_var   = "od",
    model_names    = c("logistic4", "gompertz4"),
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    curve_ids      = c(1, 4),
    verbose        = FALSE
  )

  expect_equal(length(result$plates), 2)
  tbl <- summary_table(result)
  expect_equal(nrow(tbl), 2)
})


# =============================================================================
# Sample predictions: ELISA
# =============================================================================

test_that("ELISA sample predictions are produced", {
  result <- fit_calibration_freq_multiplate(
    standards      = elisa_standards,
    samples        = elisa_samples,
    response_var   = "od",
    model_names    = c("logistic4", "gompertz4"),
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    curve_ids      = c(1, 2),
    verbose        = FALSE
  )

  all_samp <- collect_samples(result)
  expect_true(!is.null(all_samp))
  expect_true("curve_id" %in% names(all_samp))
  expect_true("predicted_log10_concentration" %in% names(all_samp))
  expect_true(nrow(all_samp) > 0)
})


# =============================================================================
# Grid: ELISA data produces monotonic predicted response
# =============================================================================

test_that("ELISA grid predicted_response is monotonic", {
  std1 <- elisa_standards[elisa_standards$curve_id == 1, ]

  result <- fit_calibration_freq(
    standards      = std1,
    response_var   = "od",
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    curve_id       = "1"
  )

  y <- result$grid$predicted_response
  if (all(is.finite(y))) {
    expect_true(all(diff(y) >= -1e-8),
                info = "Predicted response not monotonically increasing")
  }
})


# =============================================================================
# 5PL vs Gompertz: plates 1-3 should prefer 5PL or logistic5
# =============================================================================

test_that("ELISA plates 1-3 (5PL-generated) can fit logistic5", {
  std_5pl <- elisa_standards[elisa_standards$curve_id %in% c(1, 2, 3), ]

  # Fit each separately with logistic4 + logistic5 + gompertz4
  for (cid in c(1, 2, 3)) {
    this_std <- std_5pl[std_5pl$curve_id == cid, ]

    result <- fit_calibration_freq(
      standards      = this_std,
      response_var   = "od",
      model_names    = c("logistic4", "logistic5", "gompertz4"),
      is_log_response    = TRUE,
      is_log_independent = TRUE,
      std_curve_conc = 10000,
      curve_id       = as.character(cid),
      verbose        = FALSE
    )

    # logistic5 should at least converge on 5PL-generated data
    l5_converged <- isTRUE(result$ensemble$logistic5$converged)
    expect_true(l5_converged,
                info = paste("logistic5 failed to converge on curve_id", cid))
  }
})


# =============================================================================
# Error handling: wrong response variable for ELISA
# =============================================================================

test_that("ELISA fitting fails with wrong response_var", {
  std1 <- elisa_standards[elisa_standards$curve_id == 1, ]

  expect_error(
    fit_calibration_freq(
      standards      = std1,
      response_var   = "mfi",
      std_curve_conc = 10000
    ),
    "not found"
  )
})
