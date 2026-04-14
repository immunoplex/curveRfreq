## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
devtools::load_all()
# other packages for loading
library(glue)
script_path <- system.file("vignette_helpers", "db_functions.R", package = "curveRfreq")
if (file.exists(script_path)) {
  source(script_path)
} else {
  warning("Script file not found at expected path: ", script_path)
}
conn <- get_db_connection()
param_user = "seamus.owen.stein@dartmouth.edu"


## ----configuration------------------------------------------------------------
verbose = TRUE
model_names = c("logistic5", "loglogistic5", "logistic4", "loglogistic4", "gompertz4") 


curve_col = "curve_id"
# define the order of all elements in the curve_id
curve_id_element_order = c("project_id", "study_accession", "experiment_accession", "feature", "source",
                                     "antigen", "plate", "nominal_sample_dilution", "wavelength")

is_display_log_response = TRUE
is_display_log_independent = TRUE


project_id = 17
study_accession = "MADI_01"
experiment_accession ="IgG1"
feature = "IgG1"
plate = "plate_3"
nominal_sample_dilution <- "1000"
wavelength <- "__none__"
source = "Standard"
antigen = "victoria"


## ----load data----------------------------------------------------------------
loaded_data <- pull_data(study_accession = study_accession, experiment_accession = experiment_accession, project_id = project_id, curve_id_elements = curve_id_element_order, conn = conn)

## -----------------------------------------------------------------------------
curve_id <- fetch_curve_id(
  lookup = loaded_data$curve_id_lookup,
  project_id = project_id,
  study_accession = study_accession,
  experiment_accession = experiment_accession,
  feature = feature,
  source = source,
  antigen = antigen,
  plate = plate,
  nominal_sample_dilution = nominal_sample_dilution,
  wavelength = wavelength,
  element_order = curve_id_element_order
)
curve_id

## -----------------------------------------------------------------------------
head(loaded_data$curve_id_lookup)

## -----------------------------------------------------------------------------
filtered <- filter_by_curve_id(loaded_data, curve_id = curve_id)

## ----response and independent variables---------------------------------------
response_var <- filtered$response_var
indep_var <- filtered$indep_var
cat("response variable:", response_var, "\n")
cat("independent variable:", indep_var)

## ----antigen constraints------------------------------------------------------
# important to go into the package. 
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

## ----study parameters---------------------------------------------------------
# importamt that gets fed into the package
study_params <- list(
  applyProzone       = TRUE,
  blank_option       = "ignored",
  is_log_response    = TRUE,
  is_log_independent = TRUE
)

print(study_params)

## ----resolve curve settings---------------------------------------------------
antigen_plate <- resolve_curve_settings(
  loaded_data = filtered,
  antigen_constraints = antigen_constraints,
  verbose = TRUE
)

## ----select antigen on the plate----------------------------------------------
# ── Select antigen data on a plate plate ──────────────────────────────────────────
# Builds curve_id for the selected values
# # triggers autodetect and all scalar arguments are unpacked automatically.
# antigen_plate <- select_antigen_plate(
#   loaded_data         = loaded_data,
#   project_id              = project_id,
#   study_accession         = study_accession,
#   experiment_accession    = experiment_accession,
#   feature                 = feature,
#   source                  = source,
#   antigen                 = antigen,
#   plate                   = plate,
#   nominal_sample_dilution = nominal_sample_dilution,
#   wavelength              = wavelength,
#   antigen_constraints = loaded_data$antigen_constraints,
#   curve_id_element_order = curve_id_element_order
# )


## ----preproccess data---------------------------------------------------------
processed_data <- preprocess_robust_curves(data = antigen_plate$plate_standard, antigen_settings = antigen_plate$antigen_settings,
                                     response_variable = response_var,
                                     independent_variable = indep_var,
                                     is_log_response = study_params$is_log_response,
                                     blank_data = antigen_plate$plate_blanks,
                                     blank_option = study_params$blank_option,
                                     is_log_independent = study_params$is_log_independent,
                                     apply_prozone = study_params$applyProzone,
                                     verbose = verbose)
processed_data$data
processed_data$antigen_fit_options

## ----model formulas-----------------------------------------------------------
formulas <- select_model_formulas(fixed_constraint = antigen_plate$fixed_a_result, response_variable = response_var,
                                    is_log_response = study_params$is_log_response, model_names = model_names)
print(formulas)

## ----model constraints--------------------------------------------------------
model_constraints <- obtain_model_constraints(data = processed_data$data,
                                                formulas = formulas,
                                                independent_variable = indep_var,
                                                response_variable = response_var,
                                                is_log_response = study_params$is_log_response,
                                                is_log_concentration = study_params$is_log_independent,
                                                antigen_settings = antigen_plate$antigen_settings,
                                                max_response = max(processed_data$data[[response_var]], na.rm = TRUE),
                                                min_response = min(processed_data$data[[response_var]], na.rm = TRUE),
                                                verbose = verbose
                                                )
model_constraints

## ----starting lists-----------------------------------------------------------
start_lists <- make_start_lists(model_constraints = model_constraints,
                                  quants = c(low = 0.2, mid = 0.5, high = 0.8))

## ----compute robust curves, collapse= TRUE------------------------------------
fit_robust_lm <- compute_robust_curves(prepped_data = processed_data$data,
                                         response_variable = response_var,
                                         independent_variable = indep_var,
                                         formulas = formulas,
                                         model_constraints = model_constraints,
                                         start_lists = start_lists,
                                         verbose = verbose)

## ----summary of fits----------------------------------------------------------
fit_summary <- summarize_model_fits(fit_robust_lm, verbose = verbose)
fit_summary

## -----------------------------------------------------------------------------
fit_params <- summarize_model_parameters(models_fit_list = fit_robust_lm, level = 0.95, model_names = model_names)

## -----------------------------------------------------------------------------
plot_data <- get_plot_data(models_fit_list = fit_robust_lm,
                             prepped_data = processed_data$data,
                             fit_params = fit_params,
                             fixed_a_result = antigen_plate$fixed_a_result,
                             curve_id_lookup = antigen_plate$curve_id_lookup,
                            # curve_id_element_order = curve_id_element_order,
                             model_names = model_names,
                             x_var = indep_var,
                             y_var = response_var)

## ----model comparisions, fig.width = 10, fig.height = 12----------------------
plot_model_comparisons(plot_data = plot_data,
                         model_names = model_names,
                         x_var = indep_var,
                         y_var = response_var,
                         use_patchwork = TRUE)

## -----------------------------------------------------------------------------
best_fit <- select_model_fit_AIC(fit_summary = fit_summary,
                                   fit_robust_lm = fit_robust_lm,
                                   fit_params = fit_params,
                                   plot_data = plot_data,
                                   verbose = verbose)

# only view first 9 rows of data frames for display purposes.
lapply(best_fit, function(x) {
  if (is.data.frame(x)) head(x, 9) else x
})

## -----------------------------------------------------------------------------
## add the tidy to the best fit object
best_fit <- tidy.nlsLM(best_fit = best_fit, fixed_a_result = antigen_plate$fixed_a_result, model_constraints =
               model_constraints, antigen_settings = antigen_plate$antigen_settings,
               antigen_fit_options = prepped_data$antigen_fit_options, 
               verbose = verbose)

best_fit$best_tidy

## ----best_fit_summary, collapse= TRUE-----------------------------------------
best_fit <- summarize_fit(best_fit = best_fit, response_variable = response_var,
                          independent_variable = indep_var,
                          fixed_a_result = antigen_plate$fixed_a_result, 
                          antigen_settings = antigen_plate$antigen_settings,
                          curve_id_lookup  = antigen_plate$curve_id_lookup,
                          antigen_fit_options = processed_data$antigen_fit_options)

best_fit$best_fit_summary

## -----------------------------------------------------------------------------
se_antigen_table <- compute_antigen_se_table(
      standards_data = loaded_data$standards,
      curve_id_lookup = loaded_data$curve_id_lookup,
      #curve_id_element_order = curve_id_element_order,
      curve_col = curve_col,
      response_col   = response_var,
      dilution_col   = "dilution",
      plate_col      = "plate",
      grouping_cols  = c("project_id", "study_accession", "experiment_accession",
                         "source", "antigen", "feature"),
      verbose        = TRUE
    )

se_antigen_table

## -----------------------------------------------------------------------------
current_se <-  lookup_antigen_se(
        se_table             = se_antigen_table,
        study_accession      = study_accession,
        experiment_accession = experiment_accession,
        source = source,
        antigen = antigen,
        feature = feature
      )

current_se


## ----collapse= TRUE-----------------------------------------------------------
best_fit <- predict_and_propagate_error(
        best_fit        = best_fit,
        response_var    = response_var,
        antigen_plate   = antigen_plate,
        study_params    = study_params,
        se_std_response = current_se,
        verbose         = verbose
      )


print(head(best_fit$sample_se))

## ----fig.width=14, fig.height= 8, out.width="100%"----------------------------
plot_standard_curve(best_fit = best_fit, is_display_log_independent = is_display_log_independent,
                    is_display_log_response = is_display_log_response,
                    pcov_threshold = antigen_constraints$pcov_threshold,
                    study_params = study_params,
                    # curve_id_element_order = curve_id_element_order,
                    # curve_col = curve_col,
                    response_variable = response_var,
                    independent_variable = indep_var)

