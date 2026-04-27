## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>"
)

## ----setup, message = FALSE---------------------------------------------------
devtools::load_all()

## ----load-data----------------------------------------------------------------
data(bead_assay_example)

# Inspect the curve/plate lookup table
head(bead_assay_example$curve_id_lookup)

## ----filter-function----------------------------------------------------------
filter_by_curve_id <- function(loaded_data,
                               curve_id,
                               target_names = c("standards", "blanks",
                                                "samples", "curve_id_lookup"),
                               verbose = FALSE) {
  
  filtered <- loaded_data
  
  filtered$curve_id_whole_lookup <- filtered$curve_id_lookup
  filtered$whole_standards <- filtered$standards
  filtered[target_names] <- lapply(filtered[target_names], function(df) {
    
    # Skip non-data frames
    if (!is.data.frame(df)) {
      if (verbose) message("[filter_by_curve_id] Skipping non-data.frame")
      return(df)
    }
    
    # Skip empty data
    if (nrow(df) == 0) {
      if (verbose) message("[filter_by_curve_id] Empty data.frame")
      return(df)
    }
    
    # Skip if no curve_id column
    if (!"curve_id" %in% names(df)) {
      if (verbose) message("[filter_by_curve_id] No curve_id column")
      return(df)
    }
    
    # Filter
    df[as.character(df$curve_id) == as.character(curve_id), , drop = FALSE]
  })
  
  return(filtered)
}

## ----filter-data--------------------------------------------------------------
curve_id   <- 6
antigen <- bead_assay_example$curve_id_lookup[bead_assay_example$curve_id_lookup$curve_id == curve_id,]$antigen
curve_data <- filter_by_curve_id(bead_assay_example, curve_id = curve_id)

## ----model-settings-----------------------------------------------------------
model_names <- c("logistic5", "loglogistic5",
                 "logistic4", "loglogistic4",
                 "gompertz4")

is_display_log_response    <- TRUE
is_display_log_independent <- TRUE
verbose                    <- TRUE

## ----antigen-constraints------------------------------------------------------
antigen_constraints <- data.frame(
  antigen                      = antigen,
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
  loaded_data                = curve_data,
  study_params               = study_params,
  antigen_constraints        = antigen_constraints,
  model_names                = model_names,
  is_display_log_response    = is_display_log_response,
  is_display_log_independent = is_display_log_independent,
  verbose                    = verbose
)

# Printing the object shows pipeline status at any time
sc

## ----set-curve-settings-------------------------------------------------------
sc$set_curve_settings()

## ----set-curve-settings-inspect-----------------------------------------------
head(sc$antigen_plate$plate_standard)
sc$antigen_plate$antigen_settings

## ----fit----------------------------------------------------------------------
sc$fit()

## ----fit-internals------------------------------------------------------------
# Raw candidate fits
sc$fit_summary          # AIC, BIC, RSS for each model
sc$fit_params           # parameter estimates and CIs per model

## ----summarize----------------------------------------------------------------
sc$summarize()

sc$best_fit$best_model_name
sc$best_fit$best_parameters        # parameter table
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

# Parameter estimates
results$best_parameters

# prediction grid from the best fitting model.
head(results$best_pred)

# Standards used for fitting
head(results$best_standard)
# Precision profile (populated after $propagate_error())
head(results$sample_se)

## ----multi-antigen, eval = FALSE----------------------------------------------
# # Define per-antigen constraints as a named list
# antigen_constraints_list <- list(
#   alpha = data.frame(
#     antigen                      = "alpha",
#     l_asy_min_constraint         = 0,
#     l_asy_max_constraint         = 0,
#     l_asy_constraint_method      = "default",
#     standard_curve_concentration = 10000,
#     pcov_threshold               = 15,
#     stringsAsFactors             = FALSE
#   ),
#   beta = data.frame(
#     antigen                      = "beta",
#     l_asy_min_constraint         = 0,
#     l_asy_max_constraint         = 0,
#     l_asy_constraint_method      = "default",
#     standard_curve_concentration = 10000,
#     pcov_threshold               = 15,
#     stringsAsFactors             = FALSE
#   )
# )
# 
# all_summaries <- vector("list", nrow(bead_assay_example$curve_id_lookup))
# 
# for (i in seq_len(nrow(bead_assay_example$curve_id_lookup))) {
# 
#   cid        <- bead_assay_example$curve_id_lookup$curve_id[i]
#   antigen    <- bead_assay_example$curve_id_lookup$antigen[i]   # pull antigen name
#   curve_data <- filter_by_curve_id(bead_assay_example, curve_id = cid)
# 
#   # Update the StandardCurve object with new data AND new constraints
#   sc$set_data(curve_data)
#   sc$antigen_constraints <- antigen_constraints_list[[antigen]]  # swap constraints
# 
#   sc$set_curve_settings()$fit()$summarize()
#   all_summaries[[i]] <- sc$get_results()$best_fit_summary
# }
# 
# combined_summary <- do.call(rbind, all_summaries)
# combined_summary

## ----quick-ref, eval = FALSE--------------------------------------------------
# # # Full single-plate workflow using the built-in example data
# # data(bead_assay_example)
# #
# # curve_data <- filter_by_curve_id(bead_assay_example, curve_id = 1)
# #
# # sc <- StandardCurve$new(curve_data, study_params, antigen_constraints,
# #                         model_names = model_names,
# #                         is_display_log_response    = TRUE,
# #                         is_display_log_independent = TRUE)
# #
# # sc$set_curve_settings()
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

