## ────────────────────────────────────────────────────────────────────
## Feature sentinel — same pattern as wavelength.
## Use AFTER running migration_phase2_not_null_nk.sql.
## ────────────────────────────────────────────────────────────────────

#' Sentinel value for missing feature
#'
#' Constant used to represent a missing or undefined feature in grouping keys.
#' Mirrors the behavior of wavelength sentinel values to maintain consistency
#' in natural key construction.
#'
#' @format A character scalar.
#' @examples
#' FEAT_NONE
#' @export
FEAT_NONE <- "__none__"

#' Sentinel value for missing wavelength
#'
#' Constant used to represent a missing or undefined wavelength in grouping keys.
#' This avoids the use of `NA` or `NULL`, ensuring compatibility with SQL
#' `UNIQUE` constraints and stable joins in R workflows, particularly for
#' bead-array data where wavelength may not be defined.
#'
#' @details
#' This value must match the sentinel used for assays that do
#' not have a wavelength (bead arrays)
#'
#' @format A character scalar.
#'
#' @examples
#' WL_NONE
#'
#' @export
WL_NONE <- "__none__"


### ---- HELPERS ----



#' Resolve the assay response column name
#'
#' Determines the appropriate response column name in a data frame using
#' metadata when available. Falls back to commonly used response columns
#' for backward compatibility.
#'
#' @param df A data.frame containing assay data (e.g., standards, blanks, samples).
#' @param default Character. Fallback column name (default = `"mfi"`).
#'
#' @return A character string indicating the column name to use for response values.
#'
#' @details
#' Resolution priority:
#' \enumerate{
#'   \item `assay_response_variable` column (if present and valid)
#'   \item Provided `default` column (if present in data)
#'   \item First available among common response columns:
#'   `"mfi"`, `"absorbance"`, `"fluorescence"`, `"od"`
#'   \item Fallback to `default`
#' }
#'
#' @export
resolve_response_col <- function(df, default = "mfi") {
  if (!is.null(df) && "assay_response_variable" %in% names(df)) {
    rv <- unique(df$assay_response_variable)
    rv <- rv[!is.na(rv) & rv != ""]
    if (length(rv) == 1 && rv %in% names(df)) return(rv)
  }

  # Fall back to default if the column exists
  if (!is.null(df) && default %in% names(df)) return(default)

  # Last resort: find any response-like column
  if (!is.null(df)) {
    candidates <- intersect(c("mfi", "absorbance", "fluorescence", "od"), names(df))
    if (length(candidates) > 0) return(candidates[1])
  }

  return(default)
}


#' Create an empty summary row with NA values
#'
#' Generates a single-row data.frame with grouping variables preserved and
#' summary statistics set to `NA` or zero. Used as a safe fallback when
#' no valid observations are available for a grouping.
#'
#' @param grouping A data.frame or named list representing grouping variables.
#' @param grouping_cols Character vector of grouping column names (unused but retained for consistency).
#'
#' @return A data.frame with one row containing grouping values and NA summary fields.
#'
#' @keywords internal
.empty_se_row <- function(grouping, grouping_cols) {
  data.frame(
    grouping,
    median_se        = NA_real_,
    n_dilutions_used = 0L,
    n_plates         = 0L,
    total_obs        = 0L,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
}

#' Normalize wavelength value for bead-arrays.
#'
#' Replace NA or empty wavelength values with the sentinel.
#'
#' @param x Character or NA. Wavelength value to normalize.
#' @return Character vector.
#' @export
normalize_wavelength <- function(x) {
  ifelse(is.na(x) | trimws(x) == "", WL_NONE, as.character(x))
}


#' Attach grouping keys to output data
#'
#' Ensures that `wavelength` and `feature` columns exist in an output data frame.
#' Values are sourced from the corresponding `best_data` used to generate the output.
#'
#' @param df A data.frame containing model output or summary results.
#' @param best_data A data.frame used to derive grouping keys (typically input to model fitting).
#' @param context Optional character string used for logging/debugging context.
#'
#' @return A data.frame with `wavelength` and `feature` columns guaranteed to exist.
#'
#' @details
#' This function enforces consistent natural keys across all outputs generated
#' from model fitting workflows. If keys are missing:
#' \itemize{
#'   \item `wavelength` is filled using `best_data` or `WL_NONE`
#'   \item `feature` is filled using `best_data` or `FEAT_NONE`
#' }
#'
#' The `wavelength` column is normalized using `normalize_wavelength()`.
#'
#' @seealso \code{\link{normalize_wavelength}}
#'
#' @export
attach_grouping_keys <- function(df, best_data, context = "") {
  if (!"wavelength" %in% names(df)) {
    wl <- if ("wavelength" %in% names(best_data)) unique(best_data$wavelength)[1] else WL_NONE
    df$wavelength <- wl
  }
  df$wavelength <- normalize_wavelength(df$wavelength)

  if (!"feature" %in% names(df)) {
    feat <- if ("feature" %in% names(best_data)) unique(best_data$feature)[1] else FEAT_NONE
    df$feature <- feat
  }

  if (context != "") {
    message(sprintf(
      "[attach_grouping_keys] %s: wavelength=%s, feature=%s",
      context,
      paste(unique(df$wavelength), collapse = ","),
      paste(unique(df$feature), collapse = ",")
    ))
  }

  df
}


#' Null-coalescing operator
#'
#' Returns the left-hand side if it is not `NULL`, otherwise returns the right-hand side.
#'
#' @param a First value.
#' @param b Fallback value if `a` is `NULL`.
#'
#' @return The first non-NULL value.
#'
#' @examples
#' NULL %||% 5
#' 10 %||% 5
#'
#' @export
`%||%` <- function(a, b) if (!is.null(a)) a else b
