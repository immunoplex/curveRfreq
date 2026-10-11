# Censoring, predictive intervals and threshold probabilities for test
# samples, including samples whose response is beyond the curve's asymptotes
# (previously untested).

library(curveRcore)
data("bead_assay_example", package = "curveRcore")

oor_std <- curveRcore::preprocess_standards(
  data = bead_assay_example$standards[bead_assay_example$standards$curve_id == 1, ],
  antigen_settings = list(standard_curve_concentration = 10000),
  response_variable = "mfi", independent_variable = "concentration",
  is_log_response = TRUE, is_log_independent = TRUE, apply_prozone = TRUE)$data

oor_samp <- bead_assay_example$samples[bead_assay_example$samples$curve_id == 1, ]
oor_samp <- oor_samp[seq_len(min(8, nrow(oor_samp))), ]
# one sample far below the lower asymptote, one far above the upper
extremes <- oor_samp[1:2, ]
extremes$sampleid <- c("very_low", "very_high")
extremes$mfi <- c(1e-3, 1e9)
oor_samp <- rbind(oor_samp, extremes)

new_cols <- c("conc_q_lo", "conc_q_med", "conc_q_hi", "p_below_lloq", "p_above_uloq",
              "p_below_cutoff", "frac_draws_below_a", "frac_draws_above_d")

fit_oor <- function(...) {
  fit_calibration_freq_multiplate(standards = oor_std, samples = oor_samp,
                                  response_var = "mfi", std_curve_conc = 10000, ...)
}

test_that("every sample gets the interval / probability columns", {
  mp <- fit_oor(decision_cutoff = 50, interval_seed = 11)
  s  <- mp$plates[["1"]]$samples
  expect_true(all(new_cols %in% names(s)))
  ok <- is.finite(s$conc_q_lo) & is.finite(s$conc_q_hi)
  expect_true(any(ok))
  expect_true(all(s$conc_q_lo[ok] <= s$conc_q_med[ok] & s$conc_q_med[ok] <= s$conc_q_hi[ok]))
  expect_true(all(s$p_below_cutoff >= 0 & s$p_below_cutoff <= 1, na.rm = TRUE))
})

test_that("out-of-asymptote samples are censored on the correct side", {
  mp <- fit_oor(interval_seed = 11)
  mp <- curveRcore::compute_detection_limits_multiplate(mp)
  mp <- curveRcore::classify_pcov_gates_multiplate(mp)
  mp <- curveRcore::classify_censoring_multiplate(mp)
  s  <- mp$plates[["1"]]$samples
  lo <- s[s$sampleid == "very_low", ]; hi <- s[s$sampleid == "very_high", ]

  expect_true(lo$censor_class %in% c("below_lod", "below_lloq"))
  expect_equal(lo$censor_type, "left")
  expect_true(is.finite(lo$censor_upper))
  expect_equal(lo$frac_draws_below_a, 1)
  expect_equal(lo$conc_q_hi, 0)
  expect_equal(lo$p_below_lloq, 1)

  expect_true(hi$censor_class %in% c("saturated", "above_uloq"))
  expect_equal(hi$censor_type, "right")
  expect_equal(hi$frac_draws_above_d, 1)
  expect_true(is.na(hi$conc_q_lo))
  expect_equal(hi$p_above_uloq, 1)
})

test_that("se_concentration / pcov are unchanged and the global RNG is untouched", {
  set.seed(42); base <- fit_oor(n_interval_draws = 0); after_base <- runif(1)
  set.seed(42); full <- fit_oor(interval_seed = 3);    after_full <- runif(1)
  sb <- base$plates[["1"]]$samples; sf <- full$plates[["1"]]$samples
  expect_identical(sb$se_concentration, sf$se_concentration)
  expect_identical(sb$pcov, sf$pcov)
  expect_identical(sb$pcov_pass, sf$pcov_pass)
  expect_identical(after_base, after_full)
  expect_true(all(is.na(sb$conc_q_med)))            # disabled -> NA columns
})

test_that("parameter draws are persisted by default, named like calib_param terms", {
  mp  <- fit_oor(interval_seed = 9, n_interval_draws = 400)
  cr  <- mp$plates[["1"]]
  best <- cr$selection$best_model_name
  dr  <- cr$ensemble[[best]]$draws
  expect_false(is.null(dr))                               # persist_draws not needed
  expect_true(all(names(dr) %in% cr$ensemble[[best]]$parameters$term))
  expect_false(any(grepl("^p[0-9]+$", names(dr))))        # the old .mvn_draws() names
  expect_true(all(lengths(dr) == 400))
})

test_that("the stored draws reproduce the samples' intervals exactly", {
  mp  <- fit_oor(interval_seed = 9, n_interval_draws = 400)
  cr  <- mp$plates[["1"]]; best <- cr$selection$best_model_name
  P   <- do.call(cbind, cr$ensemble[[best]]$draws)
  expect_equal(cr$meta$interval_seed, 10)               # multiplate: seed + curve index
  s   <- cr$samples
  fns <- curveRcore::make_inv_and_grad_fixed(best, fixed_a = NULL)
  se_resp <- cr$ensemble[[best]]$parameters$estimate[cr$ensemble[[best]]$parameters$term == "sigma_resid"]
  lx  <- curveRfreq:::.mc_conc_draws(s$observed_response_fit, P, fns, NULL,
                                     se_resp, cr$meta$interval_seed, TRUE)
  re  <- curveRcore::summarize_conc_draws(lx, s$dilution,
           cr$ensemble[[best]]$eligibility$lloq, cr$ensemble[[best]]$eligibility$uloq)
  expect_equal(re$conc_q_med, s$conc_q_med)
  expect_equal(re$frac_draws_above_d, s$frac_draws_above_d)
})

test_that("a fixed lower asymptote is recorded, not drawn", {
  mp <- fit_oor(interval_seed = 9, n_interval_draws = 200, fixed_a = 1)
  cr <- mp$plates[["1"]]; best <- cr$selection$best_model_name
  skip_if(is.na(best))
  expect_false("a" %in% names(cr$ensemble[[best]]$draws))
  expect_equal(cr$ensemble[[best]]$draws_fixed$a, 1)
})

test_that("interval_seed makes the draws reproducible", {
  a <- fit_oor(interval_seed = 7)$plates[["1"]]$samples
  b <- fit_oor(interval_seed = 7)$plates[["1"]]$samples
  expect_identical(a[new_cols], b[new_cols])
})
