## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>"
)

## ----setup, message = FALSE---------------------------------------------------
devtools::load_all()

## ----load-data----------------------------------------------------------------
script_path <- system.file("vignette_helpers", "db_functions.R", package = "curveRfreq")
source(script_path)
library(glue)
conn <- get_db_connection()

curve_id_elements = c("project_id", "study_accession", "experiment_accession", "feature", "source",
                                     "antigen", "plate", "nominal_sample_dilution", "wavelength")

loaded_data <- pull_data(study_accession = "MADI_01",
                         experiment_accession = "IgG1", 
                         project_id = 17, 
                         curve_id_elements = curve_id_elements,
                         conn = conn)

head(loaded_data$curve_id_lookup)

## -----------------------------------------------------------------------------
curve_id <- 1118
curve_data <- filter_by_curve_id(loaded_data, curve_id = curve_id)

## -----------------------------------------------------------------------------
model_names <- c("logistic5", "loglogistic5",
                                 "logistic4", "loglogistic4",
                                 "gompertz4")

is_display_log_response <- TRUE
is_display_log_independent <- TRUE
verbose <- TRUE

## ----antigen-constraints------------------------------------------------------
antigen_constraints <- data.frame(
  antigen                      = "victoria",
  l_asy_min_constraint         = 0,
  l_asy_max_constraint         = 0,
  l_asy_constraint_method      = "default",
  standard_curve_concentration = 10000,
  pcov_threshold               = 15,
  stringsAsFactors             = FALSE
)

print(antigen_constraints)

## ----study-params-------------------------------------------------------------
study_params <- list(
  applyProzone       = TRUE,
  blank_option       = "ignored",
  is_log_response    = TRUE,
  is_log_independent = TRUE
)

## ----new----------------------------------------------------------------------
sc <- StandardCurve$new(
  loaded_data         = curve_data,
  study_params        = study_params,
  antigen_constraints = antigen_constraints,
  model_names = model_names,
  is_display_log_response    = is_display_log_response,
  is_display_log_independent = is_display_log_independent,
  verbose = verbose
)

# Printing the object shows pipeline status at any time
sc

## ----select-------------------------------------------------------------------
sc$set_curve_settings()

## ----select-inspect-----------------------------------------------------------
head(sc$antigen_plate$plate_standard)
sc$antigen_plate$antigen_settings

## ----fit----------------------------------------------------------------------
sc$fit()

## ----fit-internals------------------------------------------------------------
# Raw candidate fits
sc$fit_summary          # AIC, BIC, RSS for each model
sc$fit_params           # parameter estimates and CIs per model

# Best model outputs
# sc$best_fit$best_model_name
# sc$best_fit$best_tidy        # tidy parameter table
# sc$best_fit$best_fit_summary # QC metrics + fit stats (one row)

## ----summarize----------------------------------------------------------------
sc$summarize()

sc$best_fit$best_model_name
sc$best_fit$best_parameters        # tidy parameter table
sc$best_fit$best_fit_summary # QC metrics + fit stats (one row)

## ----compare-models, fig.width = 10, fig.height = 12, out.width = "100%"------
sc$compare_models()


## ----propagate-error----------------------------------------------------------
sc$propagate_error()

head(sc$best_fit$sample_se)

## ----plot, fig.width = 14, fig.height = 8, out.width = "100%"-----------------
sc$plot()


## ----plot-linear, fig.width = 14, fig.height = 8, out.width = "100%"----------
sc$plot(is_display_log_independent = FALSE, is_display_log_response = FALSE)

## ----get-results--------------------------------------------------------------
results <- sc$get_results()

# Tables available
names(results)

# One-row QC + fit statistics summary
results$best_fit_summary

# parameter estimates
results$best_parameters

# Precision profile (populated after $propagate_error())
head(results$sample_se)

## ----multi-antigen, eval = FALSE----------------------------------------------
# # antigens <- c("victoria", "perth", "brisbane")
# # plates   <- c("plate_1", "plate_2", "plate_3")
# #
# # all_summaries <- vector("list", length(antigens) * length(plates))
# # i <- 1
# #
# # for (ant in antigens) {
# #   for (plt in plates) {
# #     sc$select(
# #       project_id              = 17,
# #       study_accession         = "MADI_01",
# #       experiment_accession    = "IgG1",
# #       feature                 = "IgG1",
# #       source                  = "Standard",
# #       antigen                 = ant,
# #       plate                   = plt,
# #       nominal_sample_dilution = "1000"
# #     )$fit()
# #
# #     all_summaries[[i]] <- sc$get_results()$best_fit_summary
# #     i <- i + 1
# #   }
# # }
# #
# # combined_summary <- do.call(rbind, all_summaries)

## ----set-data, eval = FALSE---------------------------------------------------
# # loaded_data_2 <- pull_data("MADI_02", "IgG1", 17, conn)
# #
# # sc$set_data(loaded_data_2)
# #
# # sc$select(...)$fit()$summarize()

## ----quick-ref, eval = FALSE--------------------------------------------------
# # # Full single-antigen workflow
# # loaded_data <- pull_data(study_accession, experiment_accession, project_id, conn)
# #
# # sc <- StandardCurve$new(loaded_data, study_params, antigen_constraints)
# #
# # sc$select(project_id, study_accession, experiment_accession,
# #           feature, source, antigen, plate, nominal_sample_dilution)
# #
# # sc$fit()
# # sc$summarize()
# # sc$compare_models()
# # sc$plot()
# # sc$propagate_error()
# #
# # results <- sc$get_results()
# #
# # # Check status at any time
# # sc   # or print(sc)

