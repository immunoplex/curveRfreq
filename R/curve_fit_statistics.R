#' Tidy parameter estimates from a fitted nlsLM model
#'
#' Extracts and tidies coefficient estimates from a fitted \code{nlsLM} model,
#' attaches parameter constraints, handles fixed lower asymptote values, and
#' appends study metadata. When the lower asymptote is fixed via the
#' \code{"range_of_blanks"} method and the response is log-transformed, the
#' fixed value and its constraint bounds are log10-transformed before tidying.
#'
#' @param best_fit A named list as returned by the model selection workflow,
#'   containing at minimum:
#'   \describe{
#'     \item{best_fit}{A converged \code{nlsLM} fit object.}
#'     \item{best_model_name}{Character string naming the selected model
#'       (e.g. \code{"Y5"}).}
#'     \item{best_data}{A \code{data.frame} of the data used for fitting,
#'       containing columns \code{study_accession}, \code{experiment_accession},
#'       \code{nominal_sample_dilution}, \code{antigen}, \code{plateid},
#'       \code{plate}, and \code{source}.}
#'   }
#' @param fixed_a_result Numeric scalar or \code{NULL}. Fixed lower asymptote
#'   value. If \code{NULL} the lower asymptote is treated as a free parameter.
#'   When non-\code{NULL} and the response is log-transformed, the value is
#'   log10-transformed internally before being appended as a fixed row.
#' @param model_constraints Named list of parameter constraint matrices, one
#'   entry per model name. Each entry must be coercible to a \code{data.frame}
#'   with rownames corresponding to parameter terms and columns \code{lower}
#'   and \code{upper}.
#' @param antigen_settings A named list of antigen-level settings, requiring:
#'   \describe{
#'     \item{l_asy_constraint_method}{Character. Method used to determine the
#'       lower asymptote constraint (e.g. \code{"range_of_blanks"}).}
#'     \item{l_asy_min_constraint}{Numeric. Lower bound for the asymptote
#'       constraint.}
#'     \item{l_asy_max_constraint}{Numeric. Upper bound for the asymptote
#'       constraint.}
#'   }
#' @param antigen_fit_options A named list of fitting options, requiring:
#'   \describe{
#'     \item{is_log_response}{Logical. Whether the response variable is
#'       log10-transformed.}
#'   }
#' @param verbose Logical. If \code{TRUE} (default), prints a completion
#'   message via \code{message()}.
#'
#' @return The input \code{best_fit} list with an additional element
#'   \code{best_tidy}: a \code{tibble} with one row per parameter containing
#'   columns:
#'   \describe{
#'     \item{term}{Parameter name.}
#'     \item{lower}{Lower constraint bound from \code{model_constraints}.}
#'     \item{upper}{Upper constraint bound from \code{model_constraints}.}
#'     \item{estimate}{Parameter estimate.}
#'     \item{std_error}{Standard error of the estimate (0 for fixed
#'       parameters).}
#'     \item{statistic}{t-statistic (\code{NA} for fixed parameters).}
#'     \item{p_value}{Two-sided p-value (\code{NA} for fixed parameters).}
#'     \item{study_accession}{Study identifier.}
#'     \item{experiment_accession}{Experiment identifier.}
#'     \item{nominal_sample_dilution}{Nominal sample dilution.}
#'     \item{antigen}{Antigen name.}
#'     \item{plateid}{Plate identifier.}
#'     \item{plate}{Plate name.}
#'     \item{source}{Data source.}
#'     \item{wavelength}{Wavelength grouping key (via
#'       \code{\link{attach_grouping_keys}}).}
#'     \item{feature}{Feature grouping key (via
#'       \code{\link{attach_grouping_keys}}).}
#'   }
#'
#' @details
#' When \code{fixed_a_result} is non-\code{NULL}, a synthetic row for the
#' fixed \code{"a"} parameter is prepended to the tidy output with
#' \code{std_error = 0} and \code{statistic}/\code{p_value} set to \code{NA}.
#' Constraint bounds for the fixed asymptote are sourced from
#' \code{antigen_settings$l_asy_min_constraint} and
#' \code{antigen_settings$l_asy_max_constraint}.
#'
#' If \code{antigen_settings$l_asy_constraint_method == "range_of_blanks"} and
#' \code{antigen_fit_options$is_log_response} is \code{TRUE}, the fixed
#' asymptote and its bounds are log10-transformed before use. Constraint bounds
#' that are non-positive or non-finite are coerced to \code{NA_real_}.
#'
#' @seealso \code{\link{attach_grouping_keys}},
#'   \code{\link{validate_fixed_lower_asymptote}}
#'
#' @examples
#' \dontrun{
#' best_fit <- tidy.nlsLM(
#'   best_fit          = model_result,
#'   fixed_a_result    = 50,
#'   model_constraints = my_constraints,
#'   antigen_settings  = my_antigen_settings,
#'   antigen_fit_options = list(is_log_response = TRUE)
#' )
#' best_fit$best_tidy
#' }
#'
#' @export

tidy.nlsLM <- function(best_fit, fixed_a_result, model_constraints, antigen_settings, antigen_fit_options,  verbose = TRUE) {

  if (antigen_settings$l_asy_constraint_method == "range_of_blanks" &&
      antigen_fit_options$is_log_response) {

    # Validate before log-transforming
    fixed_a_result_validated <- validate_fixed_lower_asymptote(
      fixed_a_result, verbose = verbose
    )

    if (!is.null(fixed_a_result_validated)) {
      .eps <- 0.000005
      fixed_a_result <- log10(fixed_a_result_validated + .eps)
    } else {
      fixed_a_result <- NULL  # treat as free if invalid
    }

    # Guard constraint bounds too — only log-transform if positive
    min_c <- antigen_settings$l_asy_min_constraint
    max_c <- antigen_settings$l_asy_max_constraint

    antigen_settings$l_asy_min_constraint <- if (
      is.numeric(min_c) && is.finite(min_c) && min_c > 0
    ) log10(min_c) else NA_real_

    antigen_settings$l_asy_max_constraint <- if (
      is.numeric(max_c) && is.finite(max_c) && max_c > 0
    ) log10(max_c) else NA_real_
  }

  m_constraints <- model_constraints[[best_fit$best_model_name]]
  m_constraints_df <- as.data.frame(m_constraints)
  m_constraints_df$term <- rownames(m_constraints_df)
  rownames(m_constraints_df) <- NULL
  m_constraints_df <- m_constraints_df[, c("term", "lower", "upper")]
  if (!is.null(fixed_a_result)) {
    a_fixed_constraint <- tibble::tibble(
      term = "a",
      lower = antigen_settings$l_asy_min_constraint,
      upper = antigen_settings$l_asy_max_constraint,
    )
    m_constraints_df <-  rbind(a_fixed_constraint, m_constraints_df)
  }

  s <- summary(best_fit$best_fit)
  out <- as.data.frame(s$coefficients)
  tidy_df <- tibble::tibble(
    term = rownames(out),
    estimate = out[, "Estimate"],
    std.error = out[, "Std. Error"],
    statistic = out[, "t value"],
    p.value = out[, "Pr(>|t|)"]
  )

  tidy_df$study_accession <- unique(best_fit$best_data$study_accession)
  tidy_df$experiment_accession <- unique(best_fit$best_data$experiment_accession)
  tidy_df$nominal_sample_dilution <- unique(best_fit$best_data$nominal_sample_dilution)
  tidy_df$antigen <- unique(best_fit$best_data$antigen)
  tidy_df$plateid <- unique(best_fit$best_data$plateid)
  tidy_df$plate <- unique(best_fit$best_data$plate)
  tidy_df$source <- unique(best_fit$best_data$source)

  if (!is.null(fixed_a_result)) {
    a_fixed <- tibble::tibble(
      term = "a",
      estimate = fixed_a_result,
      std.error = 0,
      statistic = NA_real_,
      p.value = NA_real_,
      study_accession = unique(best_fit$best_data$study_accession),
      experiment_accession = unique(best_fit$best_data$experiment_accession),
      nominal_sample_dilution = unique(best_fit$best_data$nominal_sample_dilution),
      antigen = unique(best_fit$best_data$antigen),
      plateid = unique(best_fit$best_data$plateid),
      plate = unique(best_fit$best_data$plate),
      source = unique(best_fit$best_data$source)
    )
    tidy_df <- rbind(a_fixed, tidy_df)
  }

  # rename standard error column And p-value column
  names(tidy_df)[names(tidy_df) == "std.error"] <- "std_error"
  names(tidy_df)[names(tidy_df) == "p.value"] <- "p_value"


  tidy_df <- merge(tidy_df, m_constraints_df, by = "term", all.x = TRUE)
  other_cols <- setdiff(colnames(tidy_df), c("term", c("lower", "upper")))

  # New order: term → lower, upper → rest
  tidy_df <- tidy_df[, c("term", "lower", "upper", other_cols)]

  tidy_df <- attach_grouping_keys(tidy_df, best_fit$best_data, context = "tidy.nlsLM")

  best_fit$best_tidy <- tidy_df
  if (verbose) {
    message("Finished tidy.nlsLM")
  }
  return(best_fit)
}


