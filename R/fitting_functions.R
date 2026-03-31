  #' #' Obtain the lower constraint
  #' #' This function returns the lower and upper constraints for an antigen given its method
  #' #' methods are  ['default','user_defined','range_of_blanks', 'geometric_mean_of_blanks']
  #' #' @param dat standard curve data
  #' #' @param antigen name of antigen
  #' #' @param study_accession name of study accession
  #' #' @param experiment_accession name of experiment accession
  #' #' @param plate name of plate  from standard curve data
  #' #' @param plate_blanks the blanks on the std curvw dara
  #' obtain_lower_constraint <- function(dat, antigen, study_accession, experiment_accession,
  #'                                     plate, plateid, plate_blanks, antigen_constraints,
  #'                                     response_col = NULL) {
  #'
  #'   # Resolve the response column name dynamically (mfi for bead array, absorbance for ELISA)
  #'   if (is.null(response_col)) response_col <- resolve_response_col(dat)
  #'
  #'   # Handle case where antigen_constraints is a dataframe with multiple rows
  #'   # Take the first row to ensure scalar values for all constraint parameters
  #'   if (is.data.frame(antigen_constraints) && nrow(antigen_constraints) > 1) {
  #'     warning(paste("Multiple constraint rows found for antigen:", antigen,
  #'                   "- using first row. Consider deduplicating antigen_constraints."))
  #'     antigen_constraints <- antigen_constraints[1, , drop = FALSE]
  #'   }
  #'
  #'   # Extract scalar values from antigen_constraints to avoid "condition has length > 1" errors
  #'   # Use helper function to safely extract first non-NA value
  #'   safe_extract <- function(x, default = NA) {
  #'     if (is.null(x) || length(x) == 0) return(default)
  #'     x <- x[!is.na(x)]
  #'     if (length(x) == 0) return(default)
  #'     return(x[1])
  #'   }
  #'
  #'   constraint_method <- safe_extract(trimws(antigen_constraints$l_asy_constraint_method), "default")
  #'   l_asy_min <- safe_extract(antigen_constraints$l_asy_min_constraint, 0)
  #'   l_asy_max <- safe_extract(antigen_constraints$l_asy_max_constraint, NA)
  #'   std_curve_conc <- safe_extract(antigen_constraints$standard_curve_concentration, 10000)
  #'   pcov_thresh <- safe_extract(antigen_constraints$pcov_threshold, 20)
  #'
  #'   # Blank SE calculation using dynamic response column
  #'   blank_response_col <- resolve_response_col(plate_blanks, default = response_col)
  #'   if (nrow(plate_blanks) > 1) {
  #'     se_blank_response <- sd(plate_blanks[[blank_response_col]], na.rm = TRUE) /
  #'       sqrt(sum(!is.na(plate_blanks[[blank_response_col]])))
  #'   } else {
  #'     se_blank_response <- 0
  #'   }
  #'
  #'   if (constraint_method == "user_defined") {
  #'     l_asy_constraints <- list(
  #'       study_accession = study_accession,
  #'       experiment_accession = experiment_accession,
  #'       plate = plate,
  #'       antigen = antigen,
  #'       l_asy_min_constraint = l_asy_min,
  #'       l_asy_max_constraint = l_asy_max,
  #'       l_asy_constraint_method = constraint_method,
  #'       std_error_blank = se_blank_response,
  #'       standard_curve_concentration = std_curve_conc,
  #'       pcov_threshold = pcov_thresh
  #'     )
  #'   } else if (constraint_method == "default") {
  #'     l_asy_constraints <- list(
  #'       study_accession = study_accession,
  #'       experiment_accession = experiment_accession,
  #'       plate = plate,
  #'       antigen = antigen,
  #'       l_asy_min_constraint = 0, # lower bound is set to 0
  #'       l_asy_max_constraint = max(dat[[response_col]], na.rm = TRUE),
  #'       l_asy_constraint_method = constraint_method,
  #'       std_error_blank = se_blank_response,
  #'       standard_curve_concentration = std_curve_conc,
  #'       pcov_threshold = pcov_thresh
  #'     )
  #'   } else if (constraint_method == "range_of_blanks") {
  #'     l_asy_constraints <- list(
  #'       study_accession = study_accession,
  #'       experiment_accession = experiment_accession,
  #'       plate = plate,
  #'       antigen = antigen,
  #'       l_asy_min_constraint = min(plate_blanks[[blank_response_col]], na.rm = TRUE),
  #'       l_asy_max_constraint = max(plate_blanks[[blank_response_col]], na.rm = TRUE),
  #'       l_asy_constraint_method = constraint_method,
  #'       std_error_blank = se_blank_response,
  #'       standard_curve_concentration = std_curve_conc,
  #'       pcov_threshold = pcov_thresh
  #'     )
  #'   } else if (constraint_method == 'geometric_mean_of_blanks') {
  #'     geometric_mean <- exp(mean(log(plate_blanks[[blank_response_col]]), na.rm = TRUE))
  #'     l_asy_constraints <- list(
  #'       study_accession = study_accession,
  #'       experiment_accession = experiment_accession,
  #'       #plateid = plateid,
  #'       plate = plate,
  #'       antigen = antigen,
  #'       l_asy_min_constraint = geometric_mean,
  #'       l_asy_max_constraint = geometric_mean,
  #'       l_asy_constraint_method = constraint_method,
  #'       std_error_blank = se_blank_response,
  #'       standard_curve_concentration = std_curve_conc,
  #'       pcov_threshold = pcov_thresh
  #'     )
  #'   } else {
  #'     return(NULL)
  #'   }
  #'   return(l_asy_constraints)
  #' }
