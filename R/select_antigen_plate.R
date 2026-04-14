#' Resolve curve settings for an antigen on a plate
#'
#' Computes modeling inputs and constraints for a selected antigen on a plate,
#' including lower asymptote handling and blank variability estimates.
#'
#' This function assumes that plate-level datasets have already been filtered
#' and derives additional settings needed
#' for curve fitting.
#'
#' @param loaded_data A list containing plate-level datasets with the following
#' elements:
#' \itemize{
#'   \item `standards` Standard curve data.
#'   \item `blanks` Blank/control data.
#'   \item `samples` Sample data.
#'   \item `curve_id_lookup` a single row containing the plate-level curve identifier
#'   including antigen, study accession and experiment_accession and the id number.
#' }
#' @param antigen_constraints Data frame or list containing antigen-specific
#' constraint rules used in `obtain_lower_constraint()`.
#' @param wavelength Optional wavelength filter. If provided, all datasets are
#' filtered to the matching normalized wavelength. Defaults to `WL_NONE`.
#' @param verbose Logical. If `TRUE`, prints diagnostic messages including row counts.
#'
#' @return A list containing:
#' \itemize{
#'   \item `plate_standard` Filtered standard curve data.
#'   \item `plate_blanks` Filtered blank data.
#'   \item `plate_samples` Filtered sample data.
#'   \item `antigen_settings` Output of `obtain_lower_constraint()`.
#'   \item `fixed_a_result` Validated lower asymptote result.
#'   \item `std_error_blank` Estimated standard error of blank measurements.
#' }
#'
#' @details
#' \strong{Wavelength filtering:}
#' If `wavelength` is provided (not `WL_NONE`), all input datasets are filtered
#' to rows matching the normalized wavelength using `normalize_wavelength()`.
#'
#' \strong{Assumptions:}
#' \itemize{
#'   \item Input data has already been subset to a single curve/plate context.
#'   \item The `standards` dataset must contain at least one row after filtering,
#'         otherwise an error is thrown.
#' }
#'
#' @examples
#' \dontrun{
#' settings <- resolve_curve_settings(
#'   loaded_data = plate_data,
#'   antigen_constraints = constraints
#' )
#' }
#'
#' @export
#'
#' @export
resolve_curve_settings <- function(
    loaded_data,
    antigen_constraints,
    wavelength = WL_NONE,
    verbose = TRUE
) {

  # ── Extract datasets ───────────────────────────────────────────────
  plate_standard <- loaded_data$standards
  plate_blanks   <- loaded_data$blanks
  plate_samples  <- loaded_data$samples
  curve_id_lookup <- loaded_data$curve_id_lookup

  # ── Wavelength filtering ───────────────────────────────────────────
  if (!is.null(wavelength) && wavelength != WL_NONE) {

    .filter_wl <- function(df) {
      if (!is.null(df) && nrow(df) > 0 && "wavelength" %in% names(df)) {
        df$wavelength <- normalize_wavelength(df$wavelength)
        mask <- df$wavelength == normalize_wavelength(wavelength)
        if (any(mask)) return(df[mask, , drop = FALSE])
      }
      df
    }

    plate_standard <- .filter_wl(plate_standard)
    plate_blanks   <- .filter_wl(plate_blanks)
    plate_samples  <- .filter_wl(plate_samples)
  }

  # ── Guard: ensure standard data exists AFTER filtering ─────────────
  if (is.null(plate_standard) || nrow(plate_standard) == 0) {
    stop("[resolve_curve_settings] No standard data after filtering")
  }

  # ── Strip dilution suffix from plate ───────────────────────────────
  # plate_c <- sub("-.*$", "", plate)

  # ── Resolve response column ────────────────────────────────────────
  response_col <- resolve_response_col(plate_standard)

  # ── Antigen settings ───────────────────────────────────────────────
  antigen_settings <- obtain_lower_constraint(
    dat                  = plate_standard,
    antigen              = curve_id_lookup$antigen,
    study_accession      = curve_id_lookup$study_accession,
    experiment_accession = curve_id_lookup$experiment_accession,
    plate                = curve_id_lookup$plate, #plate_c,
    plate_blanks         = plate_blanks,
    antigen_constraints  = antigen_constraints,
    response_col         = response_col
  )

  # ── Fixed lower asymptote ──────────────────────────────────────────
  fixed_a_result <- resolve_fixed_lower_asymptote(antigen_settings)
  fixed_a_result <- validate_fixed_lower_asymptote(
    fixed_a_result_raw = fixed_a_result,
    verbose            = verbose
  )

  # ── Blank standard error ───────────────────────────────────────────
  std_error_blank <- get_blank_se(antigen_settings = antigen_settings)

  # ── Logging ────────────────────────────────────────────────────────
  if (verbose) {
    counts <- c(
      standard = nrow(plate_standard),
      blanks   = nrow(plate_blanks),
      samples  = nrow(plate_samples)
    )

    message(
      "[resolve_curve_settings] counts: ",
      paste(names(counts), counts, sep = "=", collapse = ", ")
    )
  }

  # ── Return ─────────────────────────────────────────────────────────
  return(list(
    plate_standard   = plate_standard,
    plate_blanks     = plate_blanks,
    plate_samples    = plate_samples,
    antigen_settings = antigen_settings,
    fixed_a_result   = fixed_a_result,
    std_error_blank  = std_error_blank,
    curve_id_lookup  = curve_id_lookup
  ))
}


#' Select and Prepare Antigen Plate Data
#'
#' Filters and prepares standard curves, blanks, samples, and optional MCMC outputs
#' for a specific antigen plate, identified either by a full `curve_id` or by its
#' component fields.
#'
#' This function supports two input modes:
#' \itemize{
#'   \item \strong{curve_id mode (preferred):} Provide a full `curve_id` string.
#'   \item \strong{component mode:} Provide all individual fields used to construct
#'         the `curve_id`.
#' }
#'
#' If both are supplied, `curve_id` takes precedence and missing component arguments
#' are automatically backfilled.
#'
#' @param loaded_data A list containing datasets: `standards`, `blanks`, `samples`,
#' `mcmc_samples`, `mcmc_pred`, and `antigen_constraints`.
#' @param project_id Character. Project identifier (used if `curve_id` is NULL).
#' @param study_accession Character. Study identifier (used if `curve_id` is NULL).
#' @param experiment_accession Character. Experiment identifier.
#' @param feature Character. Feature name (e.g., IgG1).
#' @param source Character. Sample source.
#' @param antigen Character. Antigen name.
#' @param plate Character. Plate identifier.
#' @param nominal_sample_dilution Numeric or character. Nominal dilution.
#' @param curve_id_element_order named vector of order of curve identifier.
#' @param wavelength Optional wavelength filter. Defaults to `WL_NONE`.
#' @param antigen_constraints Data frame or list containing antigen constraint rules.
#' @param verbose Logical. If `TRUE`, prints diagnostic messages.
#'
#' @return A list containing:
#' \itemize{
#'   \item `plate_standard` Filtered standard curve data.
#'   \item `plate_blanks` Filtered blank data.
#'   \item `plate_samples` Filtered sample data.
#'   \item `plate_mcmc_samples` Filtered MCMC samples (if available).
#'   \item `plate_mcmc_pred` Filtered and sorted MCMC predictions (if available).
#'   \item `antigen_settings` Output of `obtain_lower_constraint()`.
#'   \item `fixed_a_result` Validated lower asymptote.
#'   \item `std_error_blank` Standard error of blank measurements.
#'   \item `curve_id` curve id object with attributes of the parts it uses to be constructed.
#' }
#'
#' @details
#' \strong{Deprecated behavior:}
#' Passing a `curve_id` string via `study_accession` is deprecated and will trigger
#' a warning. Use the `curve_id` argument instead.
#'
#' \strong{Matching logic:}
#' Filtering is performed using exact string matching on the `curve_id` column.
#' Attributes on `curve_id` are ignored, ensuring compatibility with database inputs.
#'
#' @examples
#' # Preferred usage
#' \dontrun{
#' select_antigen_plate(
#'   loaded_data = loaded_data,
#'   curve_id = "proj|study|exp|IgG1|serum|pt|plate1|100",
#'   antigen_constraints = loaded_data$antigen_constraints
#' )
#' }
#'
#' # Component-based usage
#' \dontrun{
#' select_antigen_plate(
#'   loaded_data = loaded_data,
#'   project_id = "proj",
#'   study_accession = "study",
#'   experiment_accession = "exp",
#'   feature = "IgG1",
#'   source = "serum",
#'   antigen = "pt",
#'   plate = "plate1",
#'   nominal_sample_dilution = 100,
#'   antigen_constraints = loaded_data$antigen_constraints
#' )
#' }
#'
#' @export
select_antigen_plate <- function(
    loaded_data,
    # curve_id = NULL,

    # Component fields (used if curve_id is NULL)
    project_id,
    study_accession,
    experiment_accession,
    feature,
    source,
    antigen,
    plate,
    nominal_sample_dilution,
    curve_id_element_order,

    wavelength          = WL_NONE,
    antigen_constraints,
    verbose             = TRUE
    ) {

      # ── Backward compatibility: deprecated usage ─────────────────────
  curve_id <-  make_curve_id_string(
    project_id              = project_id,
    study_accession         = study_accession,
    experiment_accession    = experiment_accession,
    feature                 = feature,
    source                  = source,
    antigen                 = antigen,
    plate                   = plate,
    nominal_sample_dilution = nominal_sample_dilution,
    wavelength              = wavelength,
    order                   = curve_id_element_order
  )
      # if (is.null(curve_id) &&
      #     !missing(study_accession) &&
      #     is.character(study_accession) &&
      #     length(study_accession) == 1 &&
      #     grepl("\\:", study_accession)) {
      #
      #   warning("[select_antigen_plate] Passing curve_id via `study_accession` is deprecated. Use `curve_id=` instead.")
      #   curve_id <- study_accession
      # }

      target_curve_id <- NULL
      parsed <- NULL

      # ── Mode 1: curve_id provided ────────────────────────────────────
      if (!is.null(curve_id)) {

        parsed_df <- tryCatch(
          parse_curve_id(
            data = data.frame(curve_id = curve_id, stringsAsFactors = FALSE),
            order = curve_id_element_order,
            validate = TRUE
          ),
          error = function(e) {
            stop("[select_antigen_plate] Invalid curve_id: ", conditionMessage(e))
          }
        )

        parsed <- as.list(parsed_df[1, , drop = TRUE])
        target_curve_id <- as.character(curve_id)

        if (verbose) {
          message("[select_antigen_plate] Using curve_id; parsed fields: ",
                  paste(names(parsed), parsed, sep = "=", collapse = ", "))
        }

        # Backfill missing args from parsed values
        if (missing(project_id))              project_id              <- parsed[["project_id"]]
        if (missing(study_accession))         study_accession         <- parsed[["study_accession"]]
        if (missing(experiment_accession))    experiment_accession    <- parsed[["experiment_accession"]]
        if (missing(feature))                 feature                 <- parsed[["feature"]]
        if (missing(source))                  source                  <- parsed[["source"]]
        if (missing(antigen))                 antigen                 <- parsed[["antigen"]]
        if (missing(plate))                   plate                   <- parsed[["plate"]]
        if (missing(nominal_sample_dilution)) nominal_sample_dilution <- parsed[["nominal_sample_dilution"]]
        if (missing(wavelength))              wavelength              <- parsed[["wavelength"]]

      } else {

        # ── Mode 2: build curve_id from components ─────────────────────
        required_args <- c(
          "project_id", "study_accession", "experiment_accession",
          "feature", "source", "antigen", "plate", "nominal_sample_dilution"
        )

        missing_args <- required_args[!vapply(required_args, exists, logical(1), inherits = FALSE)]

        if (length(missing_args) > 0) {
          stop("[select_antigen_plate] Missing required arguments: ",
               paste(missing_args, collapse = ", "))
        }

        # ensure wavelength is defined
        if (missing(wavelength) || is.null(wavelength)) {
          wavelength <- WL_NONE
        }

        target_curve_id <- make_curve_id_string(
          project_id              = project_id,
          study_accession         = study_accession,
          experiment_accession    = experiment_accession,
          feature                 = feature,
          source                  = source,
          antigen                 = antigen,
          plate                   = plate,
          nominal_sample_dilution = nominal_sample_dilution,
          wavelength              = wavelength,
          order                   = curve_id_element_order
        )

        # also construct parsed list for consistency
        parsed <- as.list(setNames(
          strsplit(target_curve_id, ":", fixed = TRUE)[[1]],
          curve_id_element_order
        ))
      }

      # ── Logging ──────────────────────────────────────────────────────
      if (verbose) {
        message("[select_antigen_plate] target curve_id: ", target_curve_id)
        message("[select_antigen_plate] wavelength: ", wavelength)
      }

      # ── Helper: filter by curve_id ───────────────────────────────────
      .filter_by_curve_id <- function(df, cid) {
        if (is.null(df) || nrow(df) == 0 || !"curve_id" %in% names(df)) {
          return(df)
        }
        df[df$curve_id == cid, , drop = FALSE]
      }

      # ── Filter datasets ──────────────────────────────────────────────
      plate_standard     <- .filter_by_curve_id(loaded_data$standards, target_curve_id)
      plate_blanks       <- .filter_by_curve_id(loaded_data$blanks, target_curve_id)
      plate_samples      <- .filter_by_curve_id(loaded_data$samples, target_curve_id)

      plate_mcmc_samples <- if (!is.null(loaded_data$mcmc_samples) &&
                                nrow(loaded_data$mcmc_samples) > 0) {
        .filter_by_curve_id(loaded_data$mcmc_samples, target_curve_id)
      } else {
        data.frame()
      }

      plate_mcmc_pred <- if (!is.null(loaded_data$mcmc_pred) &&
                             nrow(loaded_data$mcmc_pred) > 0) {
        .filter_by_curve_id(loaded_data$mcmc_pred, target_curve_id)
      } else {
        data.frame()
      }

      # ── Wavelength filtering ─────────────────────────────────────────
      if (!is.null(wavelength) && wavelength != WL_NONE) {

        .filter_wl <- function(df) {
          if ("wavelength" %in% names(df) && nrow(df) > 0) {
            df$wavelength <- normalize_wavelength(df$wavelength)
            mask <- df$wavelength == normalize_wavelength(wavelength)
            if (any(mask)) return(df[mask, , drop = FALSE])
          }
          df
        }

        plate_standard     <- .filter_wl(plate_standard)
        plate_blanks       <- .filter_wl(plate_blanks)
        plate_samples      <- .filter_wl(plate_samples)
        plate_mcmc_samples <- .filter_wl(plate_mcmc_samples)
        plate_mcmc_pred    <- .filter_wl(plate_mcmc_pred)
      }

      # ── Guard: no standard data ──────────────────────────────────────
      if (is.null(plate_standard) || nrow(plate_standard) == 0) {
        warning("[select_antigen_plate] No standard curve data found for curve_id: ",
                target_curve_id)
        return(NULL)
      }

      # ── Strip dilution suffix ────────────────────────────────────────
      plate_c <- sub("-.*$", "", plate)

      # ── Resolve response column ──────────────────────────────────────
      response_col <- resolve_response_col(plate_standard)

      # ── Antigen settings ─────────────────────────────────────────────
      antigen_settings <- obtain_lower_constraint(
        dat                  = plate_standard,
        antigen              = antigen,
        study_accession      = study_accession,
        experiment_accession = experiment_accession,
        plate                = plate_c,
        plate_blanks         = plate_blanks,
        antigen_constraints  = antigen_constraints,
        response_col         = response_col
      )

      # ── Fixed lower asymptote ────────────────────────────────────────
      fixed_a_result <- resolve_fixed_lower_asymptote(antigen_settings)
      fixed_a_result <- validate_fixed_lower_asymptote(
        fixed_a_result_raw = fixed_a_result,
        verbose            = verbose
      )

      # ── Blank SE ─────────────────────────────────────────────────────
      std_error_blank <- get_blank_se(antigen_settings = antigen_settings)

      # ── Sort predictions ─────────────────────────────────────────────
      if (nrow(plate_mcmc_pred) > 0 && "x" %in% names(plate_mcmc_pred)) {
        plate_mcmc_pred <- plate_mcmc_pred[order(plate_mcmc_pred$x), , drop = FALSE]
      }

      # ── Logging counts ───────────────────────────────────────────────
      if (verbose) {
        counts <- c(
          standard     = nrow(plate_standard),
          blanks       = nrow(plate_blanks),
          samples      = nrow(plate_samples),
          mcmc_samples = nrow(plate_mcmc_samples),
          mcmc_pred    = nrow(plate_mcmc_pred)
        )

        message(
          " | counts: ",
          paste(names(counts), counts, sep = "=", collapse = ", "),
          " \nCompleted"
        )
      }

      # ── Return ───────────────────────────────────────────────────────
      return(list(
        plate_standard     = plate_standard,
        plate_blanks       = plate_blanks,
        plate_samples      = plate_samples,
        plate_mcmc_samples = plate_mcmc_samples,
        plate_mcmc_pred    = plate_mcmc_pred,
        antigen_settings   = antigen_settings,
        fixed_a_result     = fixed_a_result,
        std_error_blank    = std_error_blank,

        curve_id           = target_curve_id,
        curve_id_fields    = parsed
      ))
    }
# select_antigen_plate <- function(
#     loaded_data,
#     curve_id = NULL,
#
#     # Component fields (used if curve_id is NULL)
#     project_id,
#     study_accession,
#     experiment_accession,
#     feature,
#     source,
#     antigen,
#     plate,
#     nominal_sample_dilution,
#
#     wavelength          = WL_NONE,
#     antigen_constraints,
#     verbose             = TRUE
# ) {
#
#   # ── Backward compatibility: detect deprecated usage ──────────────
#   if (is.null(curve_id) &&
#       !missing(study_accession) &&
#       is.character(study_accession) &&
#       length(study_accession) == 1 &&
#       grepl("\\:", study_accession)) {
#
#     warning("[select_antigen_plate] Passing curve_id via `study_accession` is deprecated. Use `curve_id=` instead.")
#     curve_id <- study_accession
#   }
#
#   target_curve_id <- NULL
#
#   # ── Mode 1: curve_id provided ────────────────────────────────────
#   if (!is.null(curve_id)) {
#
#     parsed_df <- parse_curve_id(
#       data = data.frame(curve_id = curve_id, stringsAsFactors = FALSE),
#       order = curve_id_element_order,
#       validate = TRUE
#     )
#
#     parsed <- as.list(parsed_df[1, curve_id_element_order])
#
#     # parsed <- tryCatch(
#     #   # parse_curve_id_string(curve_id),
#     #   error = function(e) {
#     #     stop("[select_antigen_plate] Invalid curve_id: ", conditionMessage(e))
#     #   }
#     # )
#
#     if (verbose) {
#       message("[select_antigen_plate] Using curve_id; parsed fields: ",
#               paste(names(parsed), parsed, sep = "=", collapse = ", "))
#     }
#
#     # preserve its attributes
#     curve_id_object <- curve_id
#     target_curve_id <- as.character(curve_id)
#
#     # Backfill only missing args
#     if (missing(project_id))              project_id              <- parsed[["project_id"]]
#     if (missing(study_accession))         study_accession         <- parsed[["study_accession"]]
#     if (missing(experiment_accession))    experiment_accession    <- parsed[["experiment_accession"]]
#     if (missing(feature))                 feature                 <- parsed[["feature"]]
#     if (missing(source))                  source                  <- parsed[["source"]]
#     if (missing(antigen))                 antigen                 <- parsed[["antigen"]]
#     if (missing(plate))                   plate                   <- parsed[["plate"]]
#     if (missing(nominal_sample_dilution)) nominal_sample_dilution <- parsed[["nominal_sample_dilution"]]
#
#   } else {
#
#     # ── Mode 2: build curve_id from components ─────────────────────
#     required_args <- c(
#       "project_id", "study_accession", "experiment_accession",
#       "feature", "source", "antigen", "plate", "nominal_sample_dilution"
#     )
#
#     missing_args <- required_args[!vapply(required_args, exists, logical(1), inherits = FALSE)]
#
#     if (length(missing_args) > 0) {
#       stop("[select_antigen_plate] Missing required arguments: ",
#            paste(missing_args, collapse = ", "))
#     }
#
#     target_curve_id <- as.character(make_curve_id_string(
#       project_id              = project_id,
#       study_accession         = study_accession,
#       experiment_accession    = experiment_accession,
#       feature                 = feature,
#       source                  = source,
#       antigen                 = antigen,
#       plate                   = plate,
#       nominal_sample_dilution = nominal_sample_dilution
#     ))
#   }
#
#   # ── Logging ──────────────────────────────────────────────────────
#   if (verbose) {
#     message("[select_antigen_plate] target curve_id: ", target_curve_id)
#     message("[select_antigen_plate] wavelength: ", wavelength)
#   }
#
#   # ── Helper: filter by curve_id ───────────────────────────────────
#   .filter_by_curve_id <- function(df, cid) {
#     if (is.null(df) || nrow(df) == 0 || !"curve_id" %in% names(df)) {
#       return(df)
#     }
#     df[df$curve_id == cid, , drop = FALSE]
#   }
#
#   # ── Filter datasets ──────────────────────────────────────────────
#   plate_standard     <- .filter_by_curve_id(loaded_data$standards, target_curve_id)
#   plate_blanks       <- .filter_by_curve_id(loaded_data$blanks, target_curve_id)
#   plate_samples      <- .filter_by_curve_id(loaded_data$samples, target_curve_id)
#
#   plate_mcmc_samples <- if (!is.null(loaded_data$mcmc_samples) &&
#                             nrow(loaded_data$mcmc_samples) > 0) {
#     .filter_by_curve_id(loaded_data$mcmc_samples, target_curve_id)
#   } else {
#     data.frame()
#   }
#
#   plate_mcmc_pred <- if (!is.null(loaded_data$mcmc_pred) &&
#                          nrow(loaded_data$mcmc_pred) > 0) {
#     .filter_by_curve_id(loaded_data$mcmc_pred, target_curve_id)
#   } else {
#     data.frame()
#   }
#
#   # ── Wavelength filtering ─────────────────────────────────────────
#   if (!is.null(wavelength) && wavelength != WL_NONE) {
#
#     .filter_wl <- function(df) {
#       if ("wavelength" %in% names(df) && nrow(df) > 0) {
#         df$wavelength <- normalize_wavelength(df$wavelength)
#         mask <- df$wavelength == normalize_wavelength(wavelength)
#         if (any(mask)) return(df[mask, , drop = FALSE])
#       }
#       df
#     }
#
#     plate_standard     <- .filter_wl(plate_standard)
#     plate_blanks       <- .filter_wl(plate_blanks)
#     plate_samples      <- .filter_wl(plate_samples)
#     plate_mcmc_samples <- .filter_wl(plate_mcmc_samples)
#     plate_mcmc_pred    <- .filter_wl(plate_mcmc_pred)
#   }
#
#   # ── Guard: no standard data ──────────────────────────────────────
#   if (is.null(plate_standard) || nrow(plate_standard) == 0) {
#     warning("[select_antigen_plate] No standard curve data found for curve_id: ",
#             target_curve_id)
#     return(NULL)
#   }
#
#   # ── Strip dilution suffix ────────────────────────────────────────
#   plate_c <- sub("-.*$", "", plate)
#
#   # ── Resolve response column ──────────────────────────────────────
#   response_col <- resolve_response_col(plate_standard)
#
#   # ── Antigen settings ─────────────────────────────────────────────
#   antigen_settings <- obtain_lower_constraint(
#     dat                  = plate_standard,
#     antigen              = antigen,
#     study_accession      = study_accession,
#     experiment_accession = experiment_accession,
#     plate                = plate_c,
#     plate_blanks         = plate_blanks,
#     antigen_constraints  = antigen_constraints,
#     response_col         = response_col
#   )
#
#   # ── Fixed lower asymptote ────────────────────────────────────────
#   fixed_a_result <- resolve_fixed_lower_asymptote(antigen_settings)
#   fixed_a_result <- validate_fixed_lower_asymptote(
#     fixed_a_result_raw = fixed_a_result,
#     verbose            = verbose
#   )
#
#   # ── Blank SE ─────────────────────────────────────────────────────
#   std_error_blank <- get_blank_se(antigen_settings = antigen_settings)
#
#   # ── Sort predictions ─────────────────────────────────────────────
#   if (nrow(plate_mcmc_pred) > 0 && "x" %in% names(plate_mcmc_pred)) {
#     plate_mcmc_pred <- plate_mcmc_pred[order(plate_mcmc_pred$x), , drop = FALSE]
#   }
#
#   if (verbose) {
#     counts <- c(
#       standard     = nrow(plate_standard),
#       blanks       = nrow(plate_blanks),
#       samples      = nrow(plate_samples),
#       mcmc_samples = nrow(plate_mcmc_samples),
#       mcmc_pred    = nrow(plate_mcmc_pred)
#     )
#
#     message(
#       " | counts: ",
#       paste(names(counts), counts, sep = "=", collapse = ", "),
#       " \n Completed"
#     )
#   }
#
#   # ── Return ───────────────────────────────────────────────────────
#   return(list(
#     plate_standard     = plate_standard,
#     plate_blanks       = plate_blanks,
#     plate_samples      = plate_samples,
#     plate_mcmc_samples = plate_mcmc_samples,
#     plate_mcmc_pred    = plate_mcmc_pred,
#     antigen_settings   = antigen_settings,
#     fixed_a_result     = fixed_a_result,
#     std_error_blank    = std_error_blank,
#
#     curve_id           = curve_id_object,
#     curve_id_fields    = parsed
#     # curve_id_fields    = if (!is.null(attr(curve_id_object, "fields"))) {
#     #   tryCatch(parse_curve_id_string(curve_id_object), error = function(e) NULL)
#     # } else {
#     #   NULL
#    # }
#   ))
# }



#' Construct a Standardized curve_id String
#'
#' Builds a colon-separated (`:`) `curve_id` string from named components,
#' enforcing a consistent field order defined by `order`. This ensures that
#' `curve_id` values are reproducible and schema-compliant regardless of the
#' order in which arguments are supplied.
#'
#' All required fields defined in `order` must be provided as named arguments.
#' The function will reorder inputs internally to match the specified schema.
#'
#' @param ... Named components used to construct the `curve_id`. Names must match
#'   the elements of `order`. All fields in `order` are required unless handled
#'   externally (e.g., optional defaults).
#' @param sep Character separator used to join components. Defaults to `":"`.
#' @param order Character vector defining the required fields and their order
#'   in the resulting `curve_id`.
#'
#' @return A single character string representing the standardized `curve_id`.
#'
#' @details
#' This function does not attach attributes to the returned string. Instead,
#' parsing should be performed using \code{\link{parse_curve_id}}, which relies
#' on an explicit schema (`order`) for robustness and reproducibility.
#'
#' This design avoids issues with attribute loss during data manipulation,
#' storage, or serialization (e.g., database writes, CSV export).
#'
#' @examples
#' make_curve_id_string(
#'   project_id = "17",
#'   study_accession = "MADI_01",
#'   experiment_accession = "IgG1",
#'   feature = "IgG1",
#'   source = "Standard",
#'   antigen = "pt",
#'   plate = "plate1",
#'   nominal_sample_dilution = "1x",
#'   wavelength = "450",
#'   order =  c("project_id", "study_accession", "experiment_accession", "feature", "source",
#'            "antigen", "plate", "nominal_sample_dilution", "wavelength")
#' )
#'
#' # Order of arguments does not matter (must have order argument)
#' make_curve_id_string(
#'   antigen = "pt",
#'   plate = "plate1",
#'   feature = "IgG1",
#'   source = "Standard",
#'   project_id = "17",
#'   study_accession = "MADI_01",
#'   experiment_accession = "IgG1",
#'   nominal_sample_dilution = "1x",
#'   wavelength = "450",
#'   order =  c("project_id", "study_accession", "experiment_accession", "feature", "source",
#'            "antigen", "plate", "nominal_sample_dilution", "wavelength")
#' )
#'
#' @seealso \code{\link{parse_curve_id}}
#'
#' @export
make_curve_id_string <- function(..., sep = ":", order) {
  args <- list(...)

  if (missing(order) || length(order) == 0) {
    stop("[make_curve_id_string] 'order' must be supplied and non-empty")
  }

  if (!all(order %in% names(args))) {
    stop("[make_curve_id_string] Missing required fields: ",
         paste(setdiff(order, names(args)), collapse = ", "))
  }

  values <- unlist(args[order])
  paste(values, collapse = sep)
}

# make_curve_id_string <- function(..., sep = ":") {
# args <- list(...)
#
# # Ensure all args are named
# if (is.null(names(args)) || any(names(args) == "")) {
#   stop("[make_curve_id_string] All arguments must be named.")
# }
#
# required_names <- c("feature", "antigen", "plate", "wavelength")
#
# # Check missing required fields
# missing <- setdiff(required_names, names(args))
# if (length(missing) > 0) {
#   stop("[make_curve_id_string] Missing required arguments: ",
#        paste(missing, collapse = ", "))
# }
#
# # # Optional: enforce consistent order (required first, then extras)
# # ordered_names <- c(required_names, setdiff(names(args), required_names))
# # args <- args[ordered_names]
#
# values   <- unlist(args)
# curve_id <- paste(values, collapse = sep)
# attr(curve_id, "fields") <- names(args)
#
# curve_id
# }





#' #' Parse a curve ID string into its component fields
#' #'
#' #' Splits a \code{curve_id} string (created by \code{make_curve_id_string()})
#' #' on a separator and returns a one-row data frame containing the original
#' #' \code{curve_id} plus one column per field.
#' #'
#' #' @param curve_id A character scalar with a \code{"fields"} attribute,
#' #'   as produced by \code{make_curve_id_string()}.  If \code{curve_id} was
#' #'   not created by that function, use \code{curve_id_map} instead.
#' #' @param df A data frame. Currently unused; reserved for future use or
#' #'   downstream joining.
#' #' @param sep A single character used as the field separator in
#' #'   \code{curve_id}.  Defaults to \code{":"}.
#' #'
#' #' @return A one-row \code{data.frame} with columns:
#' #'   \describe{
#' #'     \item{curve_id}{The original \code{curve_id} string.}
#' #'     \item{...}{One column per field named in the \code{"fields"} attribute,
#' #'       containing the corresponding parsed value.}
#' #'   }
#' #'
#' #' @seealso \code{\link{make_curve_id_string}}, \code{\link{curve_id_map}}
#' #'
#' #' @examples
#' #' cid <- "A:B:C"
#' #' attr(cid, "fields") <- c("study", "dose", "rep")
#' #' parse_curve_id_string(cid)
#' #' #   curve_id study dose rep
#' #' # 1    A:B:C     A    B   C
#' #'
#' #' @export
#' parse_curve_id_string <- function(curve_id, df, sep = ":") {
#'   fields <- attr(curve_id, "fields")
#'   if (is.null(fields)) {
#'     stop("[parse_curve_id_string] curve_id has no 'fields' attribute. ",
#'          "Was it created with make_curve_id_string()?  ",
#'          "If curve_id came from elsewhere (numeric, external system), ",
#'          "use curve_id_map instead.")
#'   }
#'   parts <- strsplit(curve_id, paste0("\\", sep))[[1]]
#'   if (length(parts) != length(fields)) {
#'     stop("[parse_curve_id_string] Number of fields (", length(fields),
#'          ") does not match number of parts in curve_id (", length(parts), ")")
#'   }
#'
#'   parts_df <- setNames(as.data.frame(t(parts)), fields)
#'
#'   return(cbind(data.frame(curve_id = curve_id), parts_df))
#' }

