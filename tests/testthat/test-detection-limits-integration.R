# =============================================================================
# test-detection-limits-integration.R
#
# Integration tests for detection limits, shape-LOQ, and d2y_dx2 grid
# enrichment as computed through the curveRfreq fitting pipeline.
#
# These test that the Phase 2 and Phase 3 insertions in predict_grid_freq()
# and fit_calibration_freq() produce the expected output structure with
# real data from both bead_assay_example and elisa_assay_example.
# =============================================================================

library(curveRcore)

# ── Shared test data (bead assay, reused from test-fitting.R) ────────────────

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
samples   <- bead_assay_example$samples[
  bead_assay_example$samples$curve_id %in% c(1, 2, 3), ]

std_single  <- standards[standards$curve_id == 1, ]
samp_single <- samples[samples$curve_id == 1, ]


# =============================================================================
# Phase 2: d2y_dx2 column is present on all per-model grids
# =============================================================================

test_that("d2y_dx2 column is present on all converged model grids", {
  result <- fit_calibration_freq(
    standards      = std_single,
    samples        = samp_single,
    response_var   = "mfi",
    std_curve_conc = 10000,
    curve_id       = "1",
    verbose        = FALSE
  )

  for (nm in names(result$ensemble)) {
    ens <- result$ensemble[[nm]]
    if (!isTRUE(ens$converged)) next
    expect_true("d2y_dx2" %in% names(ens$grid),
                info = paste("d2y_dx2 missing from grid of", nm))
  }
})

test_that("d2y_dx2 column is present on the plate-level grid", {
  result <- fit_calibration_freq(
    standards      = std_single,
    response_var   = "mfi",
    std_curve_conc = 10000,
    curve_id       = "1",
    verbose        = FALSE
  )
  expect_true("d2y_dx2" %in% names(result$grid))
})

test_that("d2y_dx2 has the expected shape: peak, zero-crossing, valley", {
  result <- fit_calibration_freq(
    standards      = std_single,
    response_var   = "mfi",
    std_curve_conc = 10000,
    curve_id       = "1",
    verbose        = FALSE
  )

  best_nm <- result$selection$best_model_name
  g       <- result$ensemble[[best_nm]]$grid
  ok      <- is.finite(g$d2y_dx2)

  # Should have mostly finite values (all but 2 boundary NAs)
  expect_true(sum(ok) >= nrow(g) - 2,
              info = sprintf("%d of %d finite", sum(ok), nrow(g)))

  d2  <- g$d2y_dx2[ok]
  x   <- g$log10_concentration[ok]

  # Peak should be positive, valley should be negative
  expect_gt(max(d2), 0)
  expect_lt(min(d2), 0)

  # Peak before valley (monotone increasing sigmoid)
  expect_lt(x[which.max(d2)], x[which.min(d2)])
})

test_that("d2y_dx2 is present on multiplate grids", {
  result <- fit_calibration_freq_multiplate(
    standards      = standards,
    response_var   = "mfi",
    std_curve_conc = 10000,
    verbose        = FALSE
  )

  for (cid in names(result$plates)) {
    cr <- result$plates[[cid]]
    if (is.null(cr)) next
    expect_true("d2y_dx2" %in% names(cr$grid),
                info = paste("d2y_dx2 missing from plate-level grid, curve_id", cid))
  }
})


# =============================================================================
# Phase 3a: shape-LOQ fields in eligibility
# =============================================================================

test_that("shape-LOQ fields are present in eligibility for converged models", {
  result <- fit_calibration_freq(
    standards      = std_single,
    response_var   = "mfi",
    std_curve_conc = 10000,
    curve_id       = "1",
    verbose        = FALSE
  )

  shape_fields <- c("shape_lloq_log10", "shape_uloq_log10",
                    "shape_lloq_conc",  "shape_uloq_conc",
                    "shape_lloq_response", "shape_uloq_response")

  for (nm in names(result$ensemble)) {
    ens <- result$ensemble[[nm]]
    if (!isTRUE(ens$converged)) next
    for (fld in shape_fields) {
      expect_true(fld %in% names(ens$eligibility),
                  info = paste(fld, "missing from eligibility of", nm))
    }
  }
})

test_that("shape-LLOQ < inflection < shape-ULOQ for the best model", {
  result <- fit_calibration_freq(
    standards      = std_single,
    response_var   = "mfi",
    std_curve_conc = 10000,
    curve_id       = "1",
    verbose        = FALSE
  )

  best_nm <- result$selection$best_model_name
  elig    <- result$ensemble[[best_nm]]$eligibility

  if (is.finite(elig$shape_lloq_log10) && is.finite(elig$shape_uloq_log10)) {
    # Shape-LLOQ should be less than shape-ULOQ
    expect_lt(elig$shape_lloq_log10, elig$shape_uloq_log10)

    # Both should be within the grid range
    grid_range <- range(result$grid$log10_concentration)
    expect_gte(elig$shape_lloq_log10, grid_range[1])
    expect_lte(elig$shape_uloq_log10, grid_range[2])
  }
})

test_that("shape-LOQ natural concentrations are consistent with log10 values", {
  result <- fit_calibration_freq(
    standards      = std_single,
    response_var   = "mfi",
    std_curve_conc = 10000,
    curve_id       = "1",
    verbose        = FALSE
  )

  best_nm <- result$selection$best_model_name
  elig    <- result$ensemble[[best_nm]]$eligibility

  if (is.finite(elig$shape_lloq_log10)) {
    expect_equal(elig$shape_lloq_conc, 10^elig$shape_lloq_log10,
                 tolerance = 1e-10)
  }
  if (is.finite(elig$shape_uloq_log10)) {
    expect_equal(elig$shape_uloq_conc, 10^elig$shape_uloq_log10,
                 tolerance = 1e-10)
  }
})


# =============================================================================
# Phase 3b: $detection_limits on the calibration_result
# =============================================================================

test_that("$detection_limits is present after fitting", {
  result <- fit_calibration_freq(
    standards      = std_single,
    samples        = samp_single,
    response_var   = "mfi",
    std_curve_conc = 10000,
    curve_id       = "1",
    verbose        = FALSE
  )

  expect_true("detection_limits" %in% names(result))
  dl <- result$detection_limits

  expect_equal(dl$model_name, result$selection$best_model_name)
  expect_equal(dl$method, "frequentist")
  expect_equal(dl$alpha, 0.05)
})

test_that("detection_limits contains lods and mdc_rdl sublists", {
  result <- fit_calibration_freq(
    standards      = std_single,
    response_var   = "mfi",
    std_curve_conc = 10000,
    curve_id       = "1",
    verbose        = FALSE
  )

  dl <- result$detection_limits
  expect_true("lods" %in% names(dl))
  expect_true("mdc_rdl" %in% names(dl))

  # LOD fields
  lod_fields <- c("lower_lod_response", "upper_lod_response",
                  "lower_lod_log10_conc", "upper_lod_log10_conc",
                  "lower_lod_conc", "upper_lod_conc")
  for (fld in lod_fields) {
    expect_true(fld %in% names(dl$lods), info = paste("Missing:", fld))
  }

  # MDC/RDL fields
  mdc_fields <- c("mdc_lower_log10", "mdc_upper_log10",
                  "mdc_lower_conc", "mdc_upper_conc",
                  "rdl_lower_log10", "rdl_upper_log10",
                  "rdl_lower_conc", "rdl_upper_conc")
  for (fld in mdc_fields) {
    expect_true(fld %in% names(dl$mdc_rdl), info = paste("Missing:", fld))
  }
})

test_that("LODs are derived from asymptote CIs correctly", {
  result <- fit_calibration_freq(
    standards      = std_single,
    response_var   = "mfi",
    std_curve_conc = 10000,
    curve_id       = "1",
    verbose        = FALSE
  )

  dl      <- result$detection_limits
  best_nm <- result$selection$best_model_name
  params  <- result$ensemble[[best_nm]]$parameters

  z    <- qnorm(0.975)
  a    <- params$estimate[params$term == "a"]
  se_a <- params$std_error[params$term == "a"]
  d    <- params$estimate[params$term == "d"]
  se_d <- params$std_error[params$term == "d"]

  expect_equal(dl$lods$lower_lod_response, a + z * se_a, tolerance = 1e-10)
  expect_equal(dl$lods$upper_lod_response, d - z * se_d, tolerance = 1e-10)
})

test_that("LOD concentrations are valid inversions", {
  result <- fit_calibration_freq(
    standards      = std_single,
    response_var   = "mfi",
    std_curve_conc = 10000,
    curve_id       = "1",
    verbose        = FALSE
  )

  dl      <- result$detection_limits
  best_nm <- result$selection$best_model_name
  params  <- result$ensemble[[best_nm]]$parameters
  est     <- stats::setNames(params$estimate, params$term)

  if (is.finite(dl$lods$lower_lod_log10_conc)) {
    # Manually invert and compare
    manual <- curveRcore::inv_logistic4(
      dl$lods$lower_lod_response,
      est["a"], est["b"], est["c"], est["d"])
    expect_equal(dl$lods$lower_lod_log10_conc, as.numeric(manual),
                 tolerance = 1e-8)

    # Natural-scale consistency
    expect_equal(dl$lods$lower_lod_conc, 10^dl$lods$lower_lod_log10_conc,
                 tolerance = 1e-10)
  }
})

test_that("RDL lower >= MDC lower (compressed curve shifts LLOQ right)", {
  result <- fit_calibration_freq(
    standards      = std_single,
    response_var   = "mfi",
    std_curve_conc = 10000,
    curve_id       = "1",
    verbose        = FALSE
  )

  dl <- result$detection_limits
  if (is.finite(dl$mdc_rdl$rdl_lower_log10) &&
      is.finite(dl$mdc_rdl$mdc_lower_log10)) {
    expect_gte(dl$mdc_rdl$rdl_lower_log10, dl$mdc_rdl$mdc_lower_log10)
  }
})

test_that("$detection_limits is present on all multiplate plates", {
  result <- fit_calibration_freq_multiplate(
    standards      = standards,
    samples        = samples,
    response_var   = "mfi",
    std_curve_conc = 10000,
    verbose        = FALSE
  )

  for (cid in names(result$plates)) {
    cr <- result$plates[[cid]]
    if (is.null(cr)) next
    expect_true("detection_limits" %in% names(cr),
                info = paste("detection_limits missing from curve_id", cid))
    expect_equal(cr$detection_limits$method, "frequentist")
  }
})


# =============================================================================
# ELISA: detection limits with a different response scale
# =============================================================================

test_that("detection limits work with ELISA (OD) data", {
  data("elisa_assay_example", package = "curveRcore")

  elisa_prepped <- curveRcore::preprocess_standards(
    data                 = elisa_assay_example$standards,
    antigen_settings     = list(standard_curve_concentration = 10000),
    response_variable    = "od",
    independent_variable = "concentration",
    is_log_response      = TRUE,
    is_log_independent   = TRUE,
    apply_prozone        = TRUE
  )

  std1 <- elisa_prepped$data[elisa_prepped$data$curve_id == 1, ]

  result <- fit_calibration_freq(
    standards      = std1,
    response_var   = "od",
    model_names    = c("logistic4", "gompertz4"),
    is_log_response    = TRUE,
    is_log_independent = TRUE,
    std_curve_conc = 10000,
    curve_id       = "1",
    verbose        = FALSE
  )

  # d2y_dx2 should be on the grid
  expect_true("d2y_dx2" %in% names(result$grid))

  # detection_limits should be present
  expect_true("detection_limits" %in% names(result))
  expect_equal(result$detection_limits$method, "frequentist")

  # shape-LOQ should be in eligibility for converged models
  best_nm <- result$selection$best_model_name
  if (!is.na(best_nm)) {
    expect_true("shape_lloq_log10" %in% names(result$ensemble[[best_nm]]$eligibility))
  }
})


# =============================================================================
# Grid column completeness (updated to include d2y_dx2)
# =============================================================================

test_that("grid output has all expected columns including d2y_dx2", {
  result <- fit_calibration_freq(
    standards      = std_single,
    response_var   = "mfi",
    std_curve_conc = 10000,
    curve_id       = "1"
  )
  grid_cols <- c("log10_concentration", "concentration", "x_fit",
                 "predicted_response", "ci_lower", "ci_upper",
                 "predicted_concentration", "se_concentration",
                 "pcov", "pcov_pass", "d2y_dx2")
  for (col in grid_cols) {
    expect_true(col %in% names(result$grid), info = paste("Missing:", col))
  }
})
