library(curveRcore)
library(curveRfreq)

data("bead_assay_example", package = "curveRcore")
std_raw <- bead_assay_example$standards[bead_assay_example$standards$curve_id == 1, ]
prepped <- preprocess_standards(
  data = std_raw,
  antigen_settings = list(standard_curve_concentration = 10000),
  response_variable = "mfi", independent_variable = "concentration",
  is_log_response = TRUE, is_log_independent = TRUE, apply_prozone = TRUE
)
samples <- bead_assay_example$samples[bead_assay_example$samples$curve_id == 1, ]

ref <- fit_calibration_freq(
  standards = prepped$data, samples = samples, response_var = "mfi",
  model_names = c("logistic4", "logistic5", "gompertz4"),
  is_log_response = TRUE, is_log_independent = TRUE,
  std_curve_conc = 10000, n_grid = 200, cv_x_max = 150,
  grid_min_conc = 1e-4, curve_id = "1"
)

saveRDS(ref, "tests/testthat/fixtures/reference_freq_curve1.rds")
cat("Saved reference_freq_curve1.rds\n")
