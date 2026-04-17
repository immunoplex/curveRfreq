#' Bead-based immunoassay example dataset
#'
#' A named list containing simulated multi-plate bead-based immunoassay data,
#' including standard curves, blanks, and patient samples across multiple
#' time points. The dataset spans six assay plates across two analytes
#' (alpha and beta), each with three replicate plates.
#'
#' @format A named list with six elements:
#' \describe{
#'   \item{standards}{A data frame with 60 rows and 8 columns containing
#'     standard curve data:
#'     \describe{
#'       \item{curve_id}{Integer. Plate identifier (1--6).}
#'       \item{stype}{Character. Sample type; \code{"S"} for standard.}
#'       \item{sampleid}{Character. Standard well identifier (e.g. \code{"STD_01"}).}
#'       \item{well}{Character. Plate well position.}
#'       \item{dilution}{Numeric. Known analyte concentration of the standard.}
#'       \item{mfi}{Numeric. Median fluorescence intensity.}
#'       \item{assay_response_variable}{Character. Name of the response variable (\code{"mfi"}).}
#'       \item{assay_independent_variable}{Character. Name of the independent variable (\code{"concentration"}).}
#'     }
#'   }
#'   \item{blanks}{A data frame with 24 rows and 7 columns containing blank
#'     well measurements (4 blanks per plate):
#'     \describe{
#'       \item{curve_id}{Integer. Plate identifier (1--6).}
#'       \item{stype}{Character. Sample type; \code{"B"} for blank.}
#'       \item{well}{Character. Plate well position.}
#'       \item{dilution}{Numeric. Dilution factor (1 for blanks).}
#'       \item{mfi}{Numeric. Median fluorescence intensity.}
#'       \item{assay_response_variable}{Character. Name of the response variable (\code{"mfi"}).}
#'       \item{assay_independent_variable}{Character. Name of the independent variable (\code{"concentration"}).}
#'     }
#'   }
#'   \item{samples}{A data frame with 120 rows and 13 columns containing
#'     patient sample measurements:
#'     \describe{
#'       \item{curve_id}{Integer. Plate identifier (1--6).}
#'       \item{timeperiod}{Character. Study visit; one of \code{"baseline"},
#'         \code{"month3"}, or \code{"month6"}.}
#'       \item{patientid}{Character. Patient identifier (e.g. \code{"PAT_001"}).}
#'       \item{well}{Character. Plate well position.}
#'       \item{stype}{Character. Sample type; \code{"X"} for unknown sample.}
#'       \item{sampleid}{Character. Sample identifier.}
#'       \item{agroup}{Character. Treatment group; \code{"GroupA"} or \code{"GroupB"}.}
#'       \item{dilution}{Numeric. Sample dilution factor.}
#'       \item{pctaggbeads}{Numeric. Percentage of aggregated beads.}
#'       \item{samplingerrors}{Logical. Flags for sampling errors; \code{NA} if none.}
#'       \item{mfi}{Numeric. Median fluorescence intensity.}
#'       \item{assay_response_variable}{Character. Name of the response variable (\code{"mfi"}).}
#'       \item{assay_independent_variable}{Character. Name of the independent variable (\code{"concentration"}).}
#'     }
#'   }
#'   \item{curve_id_lookup}{A data frame with 6 rows and 2 columns mapping
#'     plate names to integer curve IDs:
#'     \describe{
#'       \item{curve_id}{Character. Full plate name (e.g. \code{"alpha_STUDY_PLATE_01"}).}
#'       \item{curve_int}{Integer. Integer identifier corresponding to \code{curve_id}
#'         in other list elements.}
#'       \item{antigen}{Character. Antigen name cooresponding to the curve_id}
#'       \item{study_accession}{Character. Study accession that cooresponds to the curve_id}
#'       \item{plate}{Plate identifier}
#'     }
#'   }
#'   \item{response_var}{Character. Name of the assay response variable (\code{"mfi"}).}
#'   \item{indep_var}{Character. Name of the assay independent variable (\code{"concentration"}).}
#' }
#'
#' @examples
#' data(bead_assay_example)
#'
#' # Access standard curve data
#' head(bead_assay_example$standards)
#'
#' # Access patient samples
#' head(bead_assay_example$samples)
#'
#'@source Synthetic data generated for illustrative purposes.  Fields of the samples, standards, and blank controls
#'   are structured to align with the ImmPort data model and can be used in the Interactive Serology Plate Inspector. 
#'   Curve parameters were informed by fitting real Luminex data (MADI P3 GAPS
#'   study, ACT and PT antigens, renamed here as alpha and beta). Standard
#'   curves were simulated using the fitted parameters, with realistic
#'   plate-to-plate variability (~8--10\% random shifts in saturation and
#'   inflection point) and heteroscedastic noise matching Luminex MFI behavior.
#'   Unknown sample concentrations were drawn from across the calibration range. See
#'   \code{data-raw/bead_assay_example.R} for the data generation script.
"bead_assay_example"