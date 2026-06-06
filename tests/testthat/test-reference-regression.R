# =============================================================================
# test-reference-regression.R
#
# Regression test: fit_calibration_freq() on bead_assay_example curve_id=1
# versus a frozen reference RDS snapshot.
#
# The reference must be regenerated after any intentional change to the
# fitting pipeline:
#   Rscript data-raw/generate_reference_freq.R
# =============================================================================

library(curveRcore)

# ---- Tolerance constants ----
TOL_PARAM      <- 1e-4
TOL_RESPONSE   <- 1e-4
TOL_CONC       <- 1e-3
TOL_PCOV       <- 0.5
TOL_SAMPLE_REL <- 1e-3

pcov_threshold <- 20


# ============================================================================
# Helper: preprocess and fit exactly as the reference was generated
# ============================================================================

run_freq_curve1 <- function() {
  data("bead_assay_example", package = "curveRcore", envir = environment())

  std_raw <- bead_assay_example$standards[
    bead_assay_example$standards$curve_id == 1, ]

  prepped <- curveRcore::preprocess_standards(
    data                 = std_raw,
    antigen_settings     = list(standard_curve_concentration = 10000),
    response_variable    = "mfi",
    independent_variable = "concentration",
    is_log_response      = TRUE,
    is_log_independent   = TRUE,
    apply_prozone        = TRUE
  )

  samples <- bead_assay_example$samples[
    bead_assay_example$samples$curve_id == 1, ]

  fit_calibration_freq(
    standards          = prepped$data,
    samples            = samples,
    response_var       = "mfi",
    model_names        = c("logistic4", "logistic5", "gompertz4"),
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc     = 10000,
    n_grid             = 200,
    cv_x_max           = 150,
    grid_min_conc      = 1e-4,
    curve_id           = "1",
    verbose            = FALSE
  )
}

ref_path <- function() {
  testthat::test_path("fixtures", "reference_freq_curve1.rds")
}


# ============================================================================
# 0. Smoke test
# ============================================================================

test_that("fit_calibration_freq runs on curve_id=1 without error", {
  result <- run_freq_curve1()
  expect_s3_class(result, "calibration_result")
  expect_equal(result$meta$method, "frequentist")
  expect_equal(result$meta$curve_id, "1")
})


# ============================================================================
# 1. Reference file
# ============================================================================

test_that("reference RDS file exists and loads", {
  skip_if_not(file.exists(ref_path()),
              "Reference not found. Run: Rscript data-raw/generate_reference_freq.R")
  ref <- readRDS(ref_path())
  expect_s3_class(ref, "calibration_result")
})


# ============================================================================
# 2. Model selection
# ============================================================================

test_that("best model name matches reference", {
  skip_if_not(file.exists(ref_path()))
  ref    <- readRDS(ref_path())
  result <- run_freq_curve1()

  expect_equal(result$selection$best_model_name,
               ref$selection$best_model_name,
               info = "Best model selection changed")
})


# ============================================================================
# 3. Parameter estimates
# ============================================================================

test_that("best-model parameter estimates match reference", {
  skip_if_not(file.exists(ref_path()))
  ref    <- readRDS(ref_path())
  result <- run_freq_curve1()

  best_name <- result$selection$best_model_name
  ref_params <- ref$ensemble[[best_name]]$parameters
  cur_params <- result$ensemble[[best_name]]$parameters

  expect_false(is.null(ref_params))
  expect_false(is.null(cur_params))
  expect_equal(sort(ref_params$term), sort(cur_params$term))

  ref_params <- ref_params[order(ref_params$term), ]
  cur_params <- cur_params[order(cur_params$term), ]

  expect_equal(cur_params$estimate, ref_params$estimate,
               tolerance = TOL_PARAM, info = "Parameter estimates shifted")
  expect_equal(cur_params$std_error, ref_params$std_error,
               tolerance = TOL_PARAM, info = "Parameter SEs shifted")
})


# ============================================================================
# 4. Ensemble convergence pattern
# ============================================================================

test_that("ensemble convergence pattern matches reference", {
  skip_if_not(file.exists(ref_path()))
  ref    <- readRDS(ref_path())
  result <- run_freq_curve1()

  ref_conv <- vapply(ref$ensemble, function(e) isTRUE(e$converged), logical(1))
  cur_conv <- vapply(result$ensemble, function(e) isTRUE(e$converged), logical(1))

  expect_equal(sort(names(ref_conv)), sort(names(cur_conv)))
  common <- intersect(names(ref_conv), names(cur_conv))
  for (nm in common) {
    expect_equal(cur_conv[[nm]], ref_conv[[nm]],
                 info = paste("Convergence changed for:", nm))
  }
})


# ============================================================================
# 5. Fit statistics
# ============================================================================

test_that("fit statistics match reference for converged models", {
  skip_if_not(file.exists(ref_path()))
  ref    <- readRDS(ref_path())
  result <- run_freq_curve1()

  for (nm in names(ref$ensemble)) {
    ref_e <- ref$ensemble[[nm]]
    cur_e <- result$ensemble[[nm]]
    if (!isTRUE(ref_e$converged)) next

    expect_true(isTRUE(cur_e$converged),
                info = paste(nm, "was converged in reference but not now"))
    expect_equal(cur_e$fit_stats$aic, ref_e$fit_stats$aic,
                 tolerance = 1e-2, info = paste("AIC changed for", nm))
    expect_equal(cur_e$fit_stats$rss, ref_e$fit_stats$rss,
                 tolerance = 1e-4, info = paste("RSS changed for", nm))
  }
})


# ============================================================================
# 6. Grid predictions
# ============================================================================

test_that("grid predicted_response matches reference", {
  skip_if_not(file.exists(ref_path()))
  ref    <- readRDS(ref_path())
  result <- run_freq_curve1()

  expect_equal(nrow(result$grid), nrow(ref$grid))
  expect_equal(result$grid$predicted_response,
               ref$grid$predicted_response,
               tolerance = TOL_RESPONSE, info = "Grid predicted_response shifted")
})

test_that("grid predicted_concentration matches reference", {
  skip_if_not(file.exists(ref_path()))
  ref    <- readRDS(ref_path())
  result <- run_freq_curve1()

  ref_pc <- ref$grid$predicted_concentration
  cur_pc <- result$grid$predicted_concentration
  both_ok <- is.finite(ref_pc) & is.finite(cur_pc)

  if (any(both_ok)) {
    expect_equal(cur_pc[both_ok], ref_pc[both_ok],
                 tolerance = TOL_CONC)
  }
  expect_equal(is.na(cur_pc), is.na(ref_pc),
               info = "NA pattern changed")
})


# ============================================================================
# 7. Sample predictions
# ============================================================================

test_that("sample predictions match reference", {
  skip_if_not(file.exists(ref_path()))
  ref    <- readRDS(ref_path())
  result <- run_freq_curve1()

  if (is.null(ref$samples)) {
    expect_null(result$samples)
    skip("No samples in reference")
  }

  expect_equal(nrow(result$samples), nrow(ref$samples))

  if ("predicted_log10_concentration" %in% names(ref$samples)) {
    ref_plc <- ref$samples$predicted_log10_concentration
    cur_plc <- result$samples$predicted_log10_concentration
    both_ok <- is.finite(ref_plc) & is.finite(cur_plc)
    if (any(both_ok)) {
      expect_equal(cur_plc[both_ok], ref_plc[both_ok],
                   tolerance = TOL_CONC)
    }
  }
})


# ============================================================================
# 8. Metadata
# ============================================================================

test_that("metadata matches reference (excluding volatile fields)", {
  skip_if_not(file.exists(ref_path()))
  ref    <- readRDS(ref_path())
  result <- run_freq_curve1()

  expect_equal(result$meta$method,             ref$meta$method)
  expect_equal(result$meta$package,            ref$meta$package)
  expect_equal(result$meta$curve_id,           ref$meta$curve_id)
  expect_equal(result$meta$response_var,       ref$meta$response_var)
  expect_equal(result$meta$is_log_response,    ref$meta$is_log_response)
  expect_equal(result$meta$is_log_independent, ref$meta$is_log_independent)
  expect_equal(result$meta$n_standards,        ref$meta$n_standards)
  expect_equal(result$meta$n_samples,          ref$meta$n_samples)
})


# ============================================================================
# 9. Structural invariants
# ============================================================================

test_that("grid pcov_pass is consistent with pcov < cv_x_max", {
  result <- run_freq_curve1()
  finite_pcov <- !is.na(result$grid$pcov)
  if (any(finite_pcov)) {
    expected_pass <- result$grid$pcov < pcov_threshold
    expect_equal(result$grid$pcov_pass[finite_pcov],
                 expected_pass[finite_pcov])
  }
})

test_that("grid predicted_response is monotonically non-decreasing", {
  result <- run_freq_curve1()
  y <- result$grid$predicted_response
  if (all(is.finite(y))) {
    expect_true(all(diff(y) >= -1e-8))
  }
})


# ============================================================================
# 10. Eligibility gating — new behaviour
# ============================================================================

test_that("logistic5 is gated out and logistic4 is selected", {
  result <- run_freq_curve1()

  # The eligibility gate should reject logistic5 (boundary-pinned g, poor DR)
  l5_elig <- result$ensemble$logistic5$eligibility
  expect_false(isTRUE(l5_elig$eligible),
               info = "logistic5 should be ineligible after gating")

  # logistic4 should be eligible and selected
  l4_elig <- result$ensemble$logistic4$eligibility
  expect_true(isTRUE(l4_elig$eligible),
              info = "logistic4 should be eligible")
  expect_equal(result$selection$best_model_name, "logistic4")
  expect_false(isTRUE(result$selection$fallback),
               info = "Selection should not be a fallback")

  # selection should carry the full AIC weights table for auditability
  expect_true("weights" %in% names(result$selection))
  expect_true("model_name" %in% names(result$selection$weights))
})

test_that("selection criterion is AIC+eligibility", {
  result <- run_freq_curve1()
  expect_equal(result$selection$criterion, "AIC+eligibility")
})

test_that("all ensemble entries have eligibility assessments", {
  result <- run_freq_curve1()
  for (nm in names(result$ensemble)) {
    ens <- result$ensemble[[nm]]
    if (!isTRUE(ens$converged)) next
    expect_true(!is.null(ens$eligibility),
                info = paste("Missing eligibility for", nm))
    expect_true("eligible" %in% names(ens$eligibility),
                info = paste("Missing eligible field for", nm))
    expect_true("gates" %in% names(ens$eligibility),
                info = paste("Missing gates for", nm))
  }
})

test_that("all converged ensemble entries have per-model grids", {
  result <- run_freq_curve1()
  grid_cols <- c("log10_concentration", "pcov", "pcov_pass")
  for (nm in names(result$ensemble)) {
    ens <- result$ensemble[[nm]]
    if (!isTRUE(ens$converged)) next
    expect_true(!is.null(ens$grid),
                info = paste("Missing grid for", nm))
    for (col in grid_cols) {
      expect_true(col %in% names(ens$grid),
                  info = paste("Missing", col, "in grid for", nm))
    }
  }
})
