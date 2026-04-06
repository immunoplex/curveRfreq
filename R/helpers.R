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

#' Parse Colon-Separated Curve Identifiers into Structured Columns
#'
#' Splits a colon-delimited `curve_id` column into multiple fields based on a
#' user-supplied positional `order`, with optional validation and the ability
#' to retain only a subset of parsed fields.
#'
#' This function assumes that `curve_id` values are constructed using a
#' consistent positional encoding (e.g., via \code{paste(..., sep = ":")}),
#' where each position corresponds to a known field.
#'
#' @param data A data frame containing a column with encoded curve identifiers.
#' @param curve_col A character string specifying the name of the column
#'   containing the curve identifiers. Default is `"curve_id"`.
#' @param order A character vector defining the positional order of fields in
#'   the `curve_id`. The length of `order` must match the number of elements in
#'   each `curve_id` string.
#' @param schema Optional named list with elements:
#'   \describe{
#'     \item{required}{Character vector of field names that must be present and non-missing.}
#'     \item{optional}{Character vector of optional field names.}
#'   }
#'   If provided and \code{validate = TRUE}, required fields are checked for
#'   missing or empty values.
#' @param keep Optional character vector specifying a subset of parsed fields
#'   to retain in the output. If \code{NULL} (default), all parsed fields are kept.
#' @param sep A character string used as the delimiter in `curve_id`.
#'   Default is `":"`.
#' @param validate Logical; if \code{TRUE}, performs validation checks:
#'   \itemize{
#'     \item Ensures each `curve_id` has the same number of fields as `order`.
#'     \item If `schema` is provided, ensures required fields are present and non-missing.
#'   }
#'   Default is \code{TRUE}.
#'
#' @return A data frame containing the original `data` columns and additional
#'   parsed columns defined by `order` (or a subset if `keep` is specified).
#'
#' @details
#' The mapping from `curve_id` to columns is purely positional:
#' the i-th element of the split string is assigned to the i-th name in `order`.
#' Therefore, correct parsing depends on consistent construction of `curve_id`
#' values.
#'
#' The `keep` argument allows users to retain only specific parsed fields,
#' which is useful for downstream analysis, plotting, or modeling workflows.
#'
#' @examples
#' df <- data.frame(
#'   curve_id = c(
#'     "proj1:study1:exp1:IgG1:sample:prn:plate1:1000:450",
#'     "proj2:study2:exp2:IgG2:control:pt:plate2:2000:__none__"
#'   )
#' )
#'
#' order <- c(
#'   "project_id",
#'   "study_accession",
#'   "experiment_accession",
#'   "feature",
#'   "source",
#'   "antigen",
#'   "plate",
#'   "nominal_sample_dilution",
#'   "wavelength"
#' )
#'
#' schema <- list(
#'   required = c("feature", "antigen", "plate", "wavelength"),
#'   optional = c("project_id")
#' )
#'
#' # Parse all fields
#' parsed_all <- parse_curve_id(
#'   data = df,
#'   order = order,
#'   schema = schema
#' )
#'
#' # Keep only selected fields
#' parsed_subset <- parse_curve_id(
#'   data = df,
#'   order = order,
#'   keep = c("feature", "antigen", "plate")
#' )
#'
#' @export
parse_curve_id <- function(data,
                           curve_col = "curve_id",
                           order,
                           schema = NULL,
                           keep = NULL,
                           sep = ":",
                           validate = TRUE) {

  # checks
  if (!curve_col %in% names(data)) {
    stop(paste0("Column '", curve_col, "' not found in data"))
  }

  parts <- strsplit(data[[curve_col]], sep, fixed = TRUE)

  # validate length
  if (validate) {
    actual_lengths <- lengths(parts)

    if (any(actual_lengths != length(order))) {
      stop(
        paste0(
          "Field mismatch: expected ", length(order),
          " but found values: ",
          paste(unique(actual_lengths), collapse = ", ")
        )
      )
    }
  }

  # build df
  parsed_df <- as.data.frame(do.call(rbind, parts),
                             stringsAsFactors = FALSE)

  names(parsed_df) <- order

  # optional schema validation
  if (!is.null(schema) && validate) {
    if (!all(schema$required %in% names(parsed_df))) {
      stop("Required fields missing from parsed data")
    }

    bad_required <- apply(parsed_df[, schema$required, drop = FALSE],
                          1,
                          function(x) any(is.na(x) | x == ""))

    if (any(bad_required)) {
      stop(paste0("Missing required values in rows: ",
                  paste(which(bad_required), collapse = ", ")))
    }
  }

  # 🔥 keep subset if requested
  if (!is.null(keep)) {
    if (!all(keep %in% names(parsed_df))) {
      stop("Some 'keep' fields are not in parsed columns")
    }
    parsed_df <- parsed_df[, keep, drop = FALSE]
  }

  cbind(data, parsed_df)
}


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


#' Safely Collapse Unique Values
#'
#' Returns the unique non-missing values of a vector as a single
#' semicolon-separated string. If no non-missing values are present,
#' returns \code{NA_character_}.
#'
#' This is useful for summarizing identifiers or metadata fields
#' that may contain repeated or missing values within grouped data.
#'
#' @param x A vector of values (character, numeric, or factor).
#'
#' @details
#' The function removes \code{NA} values before computing uniqueness.
#' If multiple unique values remain, they are concatenated using
#' a semicolon (\code{";"}) separator.
#'
#' @return A character string containing unique non-missing values
#' separated by semicolons, or \code{NA_character_} if no such values exist.
#'
#' @examples
#' safe_unique(c("A", "A", "B", NA))
#' # "A;B"
#'
#' safe_unique(c(NA, NA))
#' # NA
#'
#' safe_unique(1:3)
#' # "1;2;3"
#'
#' @export
safe_unique <- function(x) {
  u <- unique(x); u <- u[!is.na(u)]
  if (length(u) == 0) NA_character_ else paste(u, collapse = ";")
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
#' @name op-null-default
#' @rdname op-null-default
#'
#' @examples
#' NULL %||% 5
#' 10 %||% 5
#'
#' @export
`%||%` <- function(a, b) if (!is.null(a)) a else b
