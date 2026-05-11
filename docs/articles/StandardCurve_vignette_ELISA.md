# Frequentist Modeling of Immunoassay Standard Curves with the StandardCurve Class: ELISA Example

## Overview

This vignette demonstrates a complete workflow for fitting standard
curves for immunoassay data, including both ELISA data, and deriving
detection metrics using the `curveRfreq` package. The workflow is
encapsulated in the `StandardCurve` R6 class, which wraps the full
`curveRfreq` pipeline into a single stateful object with a clean,
verb-driven API: the class handles everything from antigen selection
through model fitting, QC metrics, error propagation, and plotting.

We walk through the complete workflow:

1.  **Load and prepare data** — load the built-in `elisa_assay_example`
    dataset and construct the `StandardCurve` object with study settings
    and antigen constraints.
2.  **Set Curve Settings** (`$set_curve_settings()`) — Resolves
    antigen-specific curve settings for the data loaded into the
    StandardCurve object, including lower asymptote handling (fix,
    constrained, etc.) and blank variability estimates.
3.  **Fit** (`$fit()`) — fit candidate nonlinear models with the
    constants set in the object and select the best model grounded in
    information theory.
4.  **Compare** (`$compare_models()`) — Provides a visual comparison of
    the converging models fit including an overlay, residuals vs. fitted
    values, the parameter estimates, and the AIC.
5.  **Summarize** (`$summarize()`) — inspect fitted parameters and key
    QC metrics for the standard curve:
    - Inflection Point
    - Limits of Detection (LOD)
    - Minimum/Maximum Detectable Concentrations (MDC)
    - Reliable Detection Limits (RDL)
    - Curvature-based Limits of Quantification (LOQ)
6.  **Plot** (`$plot()`) — visualize the final standard curve
7.  **Propagate error** (`$propagate_error()`) — compute sample
    concentration predictions, measurement error, and the Precision
    Profile
8.  **Extract** (`$get_results()`) — pull tables for downstream use

## Setup

We begin by loading the development version of the package.

``` r
devtools::load_all()
```

## Load Example Data

The package ships with `elisa_assay_example`, a synthetic multi-plate
ELISA immunoassay dataset containing standard curves, blanks, and
patient samples across one antigen (`alpha`) run on six plates.

``` r
data(elisa_assay_example)

# Inspect the curve/plate lookup table
head(elisa_assay_example$curve_id_lookup)
#>   curve_id antigen study_accession experiment_accession   plate
#> 1        1   alpha SDYELISAexample           EXPexample plate_1
#> 2        2   alpha SDYELISAexample           EXPexample plate_2
#> 3        3   alpha SDYELISAexample           EXPexample plate_3
#> 4        4   alpha SDYELISAexample           EXPexample plate_4
#> 5        5   alpha SDYELISAexample           EXPexample plate_5
#> 6        6   alpha SDYELISAexample           EXPexample plate_6
```

`elisa_assay_example` is a plain list. The `StandardCurve` class expects
it to contain `$standards`, `$blanks`, `$samples`, `$response_var`, and
`$indep_var`.

We choose a `curve_id` from the lookup table and use a custom filtering
function (`filter_by_curve_id)` to subset the dataset of interest
(`elisa_assay_example` in this vignette). This function takes the full
dataset and the chosen `curve_id`, filters all relevant components, and
returns a reduced dataset containing only records associated with that
specific curve. The filtered dataset is then passed to the constructor.

``` r
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
```

``` r
curve_id   <- 6
antigen <- elisa_assay_example$curve_id_lookup[elisa_assay_example$curve_id_lookup$curve_id == curve_id,]$antigen
curve_data <- filter_by_curve_id(elisa_assay_example, curve_id = curve_id)
```

## Configuration

- `model_names`: Candidate nonlinear models considered during the
  fitting algorithm. Model names are specified as a vector of strings
  limited to:
  - logistic5
  - loglogistic5
  - logistic4
  - loglogistic4
  - gompertz4

**To manually exclude a model form from consideration do not include it
in the list.**

- `is_display_log_response` and `is_display_log_independent` flags
  control whether plotting results are displayed on log scales or not.

``` r
model_names <- c("logistic5", "loglogistic5",
                 "logistic4", "loglogistic4",
                 "gompertz4")

is_display_log_response    <- TRUE
is_display_log_independent <- TRUE
verbose                    <- TRUE
```

## Define Antigen Constraints and Study Parameters

These are specified once and passed to `$new()`. They apply to every
`$set_curve_settings()` / `$fit()` call on this object.

### Antigen Constraints

For the particular antigen in a study and experiment, constraints on the
lower asymptote of the standard curve and constraint method as well as
the standard curve concentration is set. The lower asymptote
`l_asy_constraint_method` can be set to the following:

- `default`: the default constraints built into the algorithm
- `user_defined`: uses `l_asy_min_constraint` and `l_asy_max_constraint`
  to constrain the lower asymptote. This is helpful if one wishes to fix
  the lower asymptote to a constant such as an MFI value of 0.
- `range_of_blanks`: uses the range of the blanks for the associated
  antigen and plate identifier as the lower asymptote.
- `geometric_mean_of_blanks`: constrains the lower asymptote for the
  associated antigen and plate identifier to the geometric mean of the
  blanks.

``` r
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
#>   antigen l_asy_min_constraint l_asy_max_constraint l_asy_constraint_method
#> 1   alpha                    0                    0                 default
#>   standard_curve_concentration pcov_threshold
#> 1                        10000             15
```

### Study Parameters

Parameters for the entire study are set in a named list before fitting
the standard curve and influence its fit. The following are the
parameters and their accepted values and definitions. The example below
uses the suggested defaults.

- **`applyProzone`**  
  Logical (`TRUE` or `FALSE`).  
  Determines whether prozone (hook effect) detection and correction is
  applied.

  - `TRUE`: Detects and adjusts for potential signal suppression at low
    dilutions.  
  - `FALSE`: No prozone correction is performed.

- **`blank_option`**  
  Character string. Controls how blank wells are handled. Options
  include:

  - `"ignored"`: Blank wells are excluded from analysis.  
  - `"subtracted"`: The geometric mean of blank responses is subtracted
    from observed responses.
  - `"subtracted_3x"`: 3 times the geometric mean of blank responses is
    subtracted from observed responses.
  - `"subtracted_10x"`: 10 times the geometric mean of blank responses
    is subtracted from observed responses.

  MFI values from test samples below zero are not possible (fluorescence
  intensity is physically constrained to non-negative values) and are
  set to zero.

- **`is_log_response`**  
  Logical (`TRUE` or `FALSE`).  
  Indicates whether the response variable (e.g., MFI or absorbance)
  should be log-transformed.

  - `TRUE`: Fit models on the log-transformed response.  
  - `FALSE`: Use the raw response values.

- **`is_log_independent`**  
  Logical (`TRUE` or `FALSE`).  
  Specifies whether the independent variable (concentration) is
  log-transformed.

  - `TRUE`: Use log-scale for the independent variable.  
  - `FALSE`: Use the original scale.

``` r
study_params <- list(
  applyProzone       = TRUE,
  blank_option       = "ignored",
  is_log_response    = TRUE,
  is_log_independent = TRUE
)
```

## Construct the StandardCurve Object

Pass the filtered data, study parameters, and constraints to `$new()`.
The class reads `response_var` and `indep_var` directly from
`curve_data`.

``` r
sc <- StandardCurve$new(
  loaded_data                = curve_data,
  study_params               = study_params,
  antigen_constraints        = antigen_constraints,
  model_names                = model_names,
  is_display_log_response    = is_display_log_response,
  is_display_log_independent = is_display_log_independent,
  verbose                    = verbose
)
#> [StandardCurve] Initialized.
#>   response variable   : od
#>   independent variable: concentration

# Printing the object shows pipeline status at any time
sc
#> 
#> <StandardCurve>
#>   Stage    : NA
#> 
#>   Steps:
#>     [ ] set_curve_settings()
#>     [ ] fit()
#>     [ ] propagate_error()
```

## Set Curve Settings

`$set_curve_settings()` receives the data in the `StandardCurve` object
containing a single antigen / plate combination and resolves
antigen-specific constraints. It resets any stale fit results
automatically so re-selecting is always safe.

``` r
sc$set_curve_settings()
#> 
#> -------------------------------------------------------
#>                   SET CURVE SETTINGS
#> -------------------------------------------------------
#> [resolve_curve_settings] counts: standard=10, blanks=4, samples=20
```

The filtered standards and antigen settings are now available for
inspection if desired.

``` r
head(sc$antigen_plate$plate_standard)
#>    curve_id stype sampleid well    dilution     od assay_response_variable
#> 51        6     S   STD_01   A1 1000.000000 0.0374                      od
#> 52        6     S   STD_02   B1  333.333333 0.0462                      od
#> 53        6     S   STD_03   C1  100.000000 0.0639                      od
#> 54        6     S   STD_04   D1   33.333333 0.0833                      od
#> 55        6     S   STD_05   E1   10.000000 0.1311                      od
#> 56        6     S   STD_06   F1    3.333333 0.5268                      od
#>    assay_independent_variable
#> 51              concentration
#> 52              concentration
#> 53              concentration
#> 54              concentration
#> 55              concentration
#> 56              concentration
sc$antigen_plate$antigen_settings
#> $study_accession
#> [1] "SDYELISAexample"
#> 
#> $experiment_accession
#> [1] "EXPexample"
#> 
#> $plate
#> [1] "plate_6"
#> 
#> $antigen
#> [1] "alpha"
#> 
#> $l_asy_min_constraint
#> [1] 0
#> 
#> $l_asy_max_constraint
#> [1] 3.6596
#> 
#> $l_asy_constraint_method
#> [1] "default"
#> 
#> $std_error_blank
#> [1] 0.003397671
#> 
#> $standard_curve_concentration
#> [1] 10000
#> 
#> $pcov_threshold
#> [1] 15
```

## Model Fitting and Selection

`$fit()` runs the entire fitting pipeline in one call across eight
internal stages: preprocessing, formula construction, constraint
computation, starting value generation, multi-start nonlinear least
squares, model comparison, AIC-based selection, and QC metric
calculation.

``` r
sc$fit()
#> 
#> -------------------------------------------------------
#>                       FIT MODELS
#> -------------------------------------------------------
#> [1/8] Preprocessing data ...
#> Applying prozone correction
#> Peak MFI = 3.6596 at concentration = 4.477121 
#> Number of points beyond the peak: 2
#> Blank Option Used: ignored
#> [2/8] Building model formulas ...
#> [3/8] Computing parameter constraints ...
#> [obtain_model_constraints] scale_class=medium, dynamic_range=1.996, slope=[0.050, 3.000], g=[0.30, 7.00]
#> $logistic5
#> $logistic5$lower
#>          a          b          c          d          g 
#> -5.0000000  0.0500000  0.1045757 -1.0279847  0.3000000 
#> 
#> $logistic5$upper
#>         a         b         c         d         g 
#> 0.5634348 3.0000000 6.3725455 1.1673056 7.0000000 
#> 
#> 
#> $loglogistic5
#> $loglogistic5$lower
#>          a          b          c          d          g 
#> -4.6989700  0.0500000  0.1045757 -1.0279847  0.3000000 
#> 
#> $loglogistic5$upper
#>         a         b         c         d         g 
#> 0.5634348 3.0000000 6.3725455 1.1673056 7.0000000 
#> 
#> 
#> $logistic4
#> $logistic4$lower
#>          a          b          c          d 
#> -4.6989700  0.0500000  0.1045757 -1.0279847 
#> 
#> $logistic4$upper
#>         a         b         c         d 
#> 0.5634348 3.0000000 6.3725455 1.1673056 
#> 
#> 
#> $loglogistic4
#> $loglogistic4$lower
#>          a          b          c          d 
#> -4.6989700  0.0500000  0.1045757 -1.0279847 
#> 
#> $loglogistic4$upper
#>         a         b         c         d 
#> 0.5634348 3.0000000 6.3725455 1.1673056 
#> 
#> 
#> $gompertz4
#> $gompertz4$lower
#>          a          b          c          d 
#> -4.6989700  0.0500000  0.1045757 -1.0279847 
#> 
#> $gompertz4$upper
#>         a         b         c         d 
#> 0.5634348 3.0000000 6.3725455 1.1673056 
#> 
#> 
#> attr(,"constraint_profile")
#> attr(,"constraint_profile")$y_min
#> [1] -1.427128
#> 
#> attr(,"constraint_profile")$y_max
#> [1] 0.5685901
#> 
#> attr(,"constraint_profile")$dynamic_range
#> [1] 1.995718
#> 
#> attr(,"constraint_profile")$conc_range
#> [1] 4.477121
#> 
#> attr(,"constraint_profile")$scale_class
#> [1] "medium"
#> 
#> attr(,"constraint_profile")$slope_max
#> [1] 3
#> 
#> attr(,"constraint_profile")$slope_min
#> [1] 0.05
#> 
#> attr(,"constraint_profile")$g_min
#> [1] 0.3
#> 
#> attr(,"constraint_profile")$g_max
#> [1] 7
#> 
#> attr(,"constraint_profile")$conc_pad_frac
#> [1] 0.7
#> 
#> attr(,"constraint_profile")$d_margin_frac
#> [1] 0.3
#> [4/8] Generating start lists ...
#> [5/8] Fitting candidate models ...
#> 
#>  Trying model: logistic5
#>   ✓ logistic5 converged (AIC=-0.89)
#> 
#>  Trying model: loglogistic5
#>   ✓ loglogistic5 converged (AIC=-18.74)
#> 
#>  Trying model: logistic4
#>   ✓ logistic4 converged (AIC=-17.43)
#> 
#>  Trying model: loglogistic4
#>   ✓ loglogistic4 converged (AIC=20.7)
#> 
#>  Trying model: gompertz4
#>   ✓ gompertz4 converged (AIC=-10.18)
#> [6/8] Summarising model fits ...
#> [7/8] Extracting parameter estimates ...
#> confint2 output:[1] "logistic5"
#>        2.5 %     97.5 %
#> a -1.8411760 -0.9296467
#> b -0.8313549  1.6583909
#> c -0.1668203  5.6957863
#> d -1.4540962  3.7887074
#> g -1.6675737  2.2675737
#> 
#> 
#> confint2 output:[1] "loglogistic5"
#>        2.5 %     97.5 %
#> a -1.5807761 -1.1746800
#> b  0.1056239  5.8943761
#> c  3.0796969  3.8819648
#> d  0.4573644  0.7700219
#> g -1.7515208  5.1574361
#> 
#> 
#> confint2 output:[1] "logistic4"
#>        2.5 %     97.5 %
#> a -1.4424143 -1.1997214
#> b  0.2384511  0.4963688
#> c  3.2597344  3.5410593
#> d  0.4743834  0.7532152
#> 
#> 
#> confint2 output:[1] "loglogistic4"
#>        2.5 %     97.5 %
#> a  -4.890275  6.0171451
#> b -14.910455 20.9104550
#> c  -4.347089 10.9755285
#> d  -2.871808  0.8158381
#> 
#> 
#> confint2 output:[1] "gompertz4"
#>       2.5 %     97.5 %
#> a -1.412669 -1.1222497
#> b  0.965319  3.1428054
#> c  3.045551  3.4094717
#> d  0.408897  0.8675226
#> Summarized Parameters completed
#> [8/8] Building plot data & selecting best model ...
#> Plot Data Completed
#> 
#> Best model: loglogistic5
#> Done. Call $summarize(), $plot(), or $propagate_error().
```

The object now holds all intermediate results as accessible fields —
nothing is hidden.

``` r
# Raw candidate fits
sc$fit_summary          # AIC, BIC, RSS for each model
#>          model converged        rss df_resid n_params         AIC         BIC
#> 1    logistic5      TRUE 0.16140677        5        5  -0.8853559   0.9301547
#> 2 loglogistic5      TRUE 0.02705860        5        5 -18.7445845 -16.9290739
#> 3    logistic4      TRUE 0.03770427        6        4 -17.4268983 -15.9139728
#> 4 loglogistic4      TRUE 1.70667624        6        4  20.6983973  22.2113228
#> 5    gompertz4      TRUE 0.07779406        6        4 -10.1839824  -8.6710569
sc$fit_params           # parameter estimates and CIs per model
#>           model parameter   estimate    conf.low  conf.high converged
#> 1     logistic5         a -1.3854114  -1.8411760 -0.9296467      TRUE
#> 2     logistic5         b  0.4135180  -0.8313549  1.6583909      TRUE
#> 3     logistic5         c  2.7644830  -0.1668203  5.6957863      TRUE
#> 4     logistic5         d  1.1673056  -1.4540962  3.7887074      TRUE
#> 5     logistic5         g  0.3000000  -1.6675737  2.2675737      TRUE
#> 6  loglogistic5         a -1.3777281  -1.5807761 -1.1746800      TRUE
#> 7  loglogistic5         b  3.0000000   0.1056239  5.8943761      TRUE
#> 8  loglogistic5         c  3.4808309   3.0796969  3.8819648      TRUE
#> 9  loglogistic5         d  0.6136932   0.4573644  0.7700219      TRUE
#> 10 loglogistic5         g  1.7029577  -1.7515208  5.1574361      TRUE
#> 11    logistic4         a -1.3210679  -1.4424143 -1.1997214      TRUE
#> 12    logistic4         b  0.3674099   0.2384511  0.4963688      TRUE
#> 13    logistic4         c  3.4003969   3.2597344  3.5410593      TRUE
#> 14    logistic4         d  0.6137993   0.4743834  0.7532152      TRUE
#> 15 loglogistic4         a  0.5634348  -4.8902754  6.0171451      TRUE
#> 16 loglogistic4         b  3.0000000 -14.9104550 20.9104550      TRUE
#> 17 loglogistic4         c  3.3142196  -4.3470893 10.9755285      TRUE
#> 18 loglogistic4         d -1.0279847  -2.8718075  0.8158381      TRUE
#> 19    gompertz4         a -1.2674595  -1.4126693 -1.1222497      TRUE
#> 20    gompertz4         b  2.0540622   0.9653190  3.1428054      TRUE
#> 21    gompertz4         c  3.2275114   3.0455511  3.4094717      TRUE
#> 22    gompertz4         d  0.6382098   0.4088970  0.8675226      TRUE
```

## Summarize

`$summarize()` computes and prints fit statistics, curve characteristics
(inflection point, LOD, MDC, RDL, LOQ), and the tidy parameter table to
the console — including the best model name, candidate model comparison
(AIC, BIC, RSS), and QC metrics. A summary of the fitted model
parameters alongside their lower and upper constraints is stored in
`$best_fit$best_parameters` as a data frame, making it available for
reporting and downstream analysis.

``` r
sc$summarize()
#> 
#> ================================================================
#>   Standard Curve Summary
#> ================================================================
#>   Best model   : loglogistic5
#>   curve_id     : 6
#>   Response var : od  |  Independent var: concentration
#> Finished tidy.nlsLM
#> [compute_inflection_point] Inflection point: (x = 3.658286, y = -0.052184)
#> MDC/RDL - mindc: 2.369, maxdc: 4.292, minrdl: 2.675, maxrdl: 4.02
#> ================================================================

sc$best_fit$best_model_name
#> [1] "loglogistic5"
sc$best_fit$best_parameters  # parameter table
#>   curve_id term      lower     upper   estimate  std_error  statistic
#> 1        6    a -4.6989700 0.5634348 -1.3777281 0.07898914 -17.441995
#> 2        6    b  0.0500000 3.0000000  3.0000000 1.12596146   2.664390
#> 3        6    c  0.1045757 6.3725455  3.4808309 0.15604791  22.306168
#> 4        6    d -1.0279847 1.1673056  0.6136932 0.06081453  10.091225
#> 5        6    g  0.3000000 7.0000000  1.7029577 1.34385081   1.267222
#>        p_value
#> 1 1.135398e-05
#> 2 4.464802e-02
#> 3 3.364127e-06
#> 4 1.636600e-04
#> 5 2.608899e-01
sc$best_fit$best_fit_summary # QC metrics + fit stats (one row)
#>   curve_id iter status   model_name         a b        c         d        g
#> 1        6   10   TRUE loglogistic5 -1.377728 3 3.480831 0.6136932 1.702958
#>   inflect_x  inflect_y     llod      ulod    mindc    maxdc   minrdl   maxrdl
#> 1  3.658286 -0.0521839 -1.17468 0.4573644 2.369143 4.292039 2.675491 4.019506
#>       lloq     uloq    lloq_y    uloq_y dfresidual nobs rsquare_fit       aic
#> 1 2.980659 3.980714 -0.816233 0.2704457          5   10   0.9959069 -18.74458
#>         bic   loglik        mse        cv bkg_method is_log_response is_log_x
#> 1 -16.92907 15.37229 0.00541172 -17.67497    ignored            TRUE     TRUE
#>   apply_prozone
#> 1          TRUE
#>                                                             formula
#> 1 od ~ a + (d - a) * (1 + g * exp(-b * (concentration - c)))^(-1/g)
```

Each row corresponds to a model parameter whose definitions are as
follows.

### Model Parameters

| Parameter | Description                                   |
|-----------|-----------------------------------------------|
| `a`       | Lower asymptote                               |
| `b`       | Slope parameter                               |
| `c`       | Inflection point (midpoint)                   |
| `d`       | Upper asymptote                               |
| `g`       | Asymmetry parameter (5-parameter models only) |

The columns include:

- `term`: Name of the model parameter
- `lower`, `upper`: Bounds for the parameter based on constraints
- `estimate`: Estimated value of the parameter
- `std_error`: Standard error of the estimate
- `statistic`: Test statistic (estimate divided by standard error)
- `p_value`: p-value of the parameter
- `curve_id`: Unique identifier for the fitted curve

The resulting `best_fit$best_fit_summary` contains one row per fitted
curve and includes:

**Curve Characteristics**

| Metric                                          | Notation                 | Description                                                                                                                                                                                                                                                                                                                                                                               |
|-------------------------------------------------|--------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Inflection Point                                | `(inflect_x, inflect_y)` | The point on the standard curve where the concavity transitions from concave up to concave down. It is the point where the assay is most sensitive to measurement errors in the measured response of the assay.                                                                                                                                                                           |
| Limits of Detection (LOD)                       | `(llod, ulod)`           | Lower and upper LODs are defined as the upper 97.5% confidence bound of the lower asymptote and the lower 2.5% confidence bound of the upper asymptote, respectively (Rajam et al.). Limits of Detection correspond to the y-coordinate in the legend, as they are defined on the response axis.                                                                                          |
| Minimum/Maximum Detectable Concentrations (MDC) | `(mindc, maxdc)`         | The smallest antibody concentration that produces a signal the assay can detect above background (Rajam et al.). This corresponds to the x-coordinate of the Lower Limit of Detection in the legend, as it is on the concentration axis.                                                                                                                                                  |
| Reliable Detection Limits (RDL)                 | `(minrdl, maxrdl)`       | Lower RDL: The lowest concentration at which the assay consistently produces a signal above background with 95% confidence based on the fit of the standard curve (Rajam et al.). Upper RDL: Analogously, the highest concentration at which the assay consistently produces a signal below the upper asymptote (saturation) with 95% confidence, based on the fit of the standard curve. |
| Limits of Quantification (LOQ)                  | `(lloq, uloq)`           | Defines a region of assay response (MFI) and concentration where sample estimates have less measurement error. Limits of Quantification are derived from the local minimum and maximum of the second derivative of x given y of the standard curve (Daly et al.), (Jeanne L Sebaugh and P. D. McCray), (Sanz et al.).                                                                     |

------------------------------------------------------------------------

**Model Fit Statistics**

| Metric                      | Notation      | Description                                                      |
|-----------------------------|---------------|------------------------------------------------------------------|
| Residual Degrees of Freedom | `dfresidual`  | Degrees of freedom remaining after model fitting                 |
| Number of Observations      | `nobs`        | Total number of data points used in the fit                      |
| R-squared                   | `rsquare_fit` | Goodness-of-fit measure                                          |
| AIC                         | `aic`         | Akaike Information Criterion for model selection                 |
| BIC                         | `bic`         | Bayesian Information Criterion for model selection               |
| Log-likelihood              | `loglik`      | Likelihood of the fitted model                                   |
| Mean Squared Error          | `mse`         | Average squared difference between observed and predicted values |
| Coefficient of Variation    | `cv`          | Relative variability of the residuals                            |

------------------------------------------------------------------------

**Model Parameters**  
Model parameter values (e.g., `a`, `b`, `c`, `d`, `g`) are included.

**Metadata**

| Field                         | Notation          | Description                                                     |
|-------------------------------|-------------------|-----------------------------------------------------------------|
| Model Type                    | `model_name`      | Selected model identifier                                       |
| Convergence Status            | `status`          | Indicates whether the model successfully converged              |
| Standard curve formula        | `formula`         | Model formula used for fitting                                  |
| Log Response Flag             | `is_log_response` | Whether response variable is log-transformed                    |
| Log independent variable Flag | `is_log_x`        | Whether independent variable (concentration) is log-transformed |
| Prozone Correction            | `apply_prozone`   | Indicates if prozone correction was applied                     |
| Background Method             | `bkg_method`      | Method used for blank/background handling                       |

## Compare Candidate Models

`$compare_models()` produces a multi-panel plot showing fitted curves,
residuals vs. fitted values, parameter estimates with confidence
intervals, and AIC scores for all converged candidate models. Use this
to validate that the AIC-selected model is sensible.

``` r
sc$compare_models()
```

![](StandardCurve_vignette_ELISA_files/figure-html/compare-models-1.png)

## Propagate Measurement Error

`$propagate_error()` constructs an SE lookup table from the standards
data, retrieves the values relevant to the selected antigen, and
propagates measurement error from the standard curve to produce the
precision profile (stored in `$best_fit$sample_se`).

``` r
sc$propagate_error()
#> 
#> -------------------------------------------------------
#>                     PROPAGATE ERROR
#> -------------------------------------------------------
#> Computing median SE for 1 unique groupings ...
#> Done. 1 / 1 groupings have a valid median SE.
#> [predict_and_propagate] se_std_response not usable (NA); using fallback se=0.008381
#> pred_se has 200 row(s)
#> 
#> === Propagation Input Diagnosis ===
#> Model         : loglogistic5 
#> coef(fit)     : a = -1.37773, b = 3, c = 3.48083, d = 0.61369, g = 1.70296 
#> vcov dim      : 5 x 5 
#> vcov rownames : a, b, c, d, g 
#> fixed_a       : NULL 
#> ===================================
#> [propagate] CV formula   : LINEAR-scale (se_x * ln(10) * 100) --- avoids /0 at log10(conc)=0
#> [propagate] Model       : loglogistic5
#> [propagate] Free params : a, b, c, d, g
#> [propagate] Sigma rows  : a, b, c, d, g
#> [propagate] fixed_a     : NULL (a is free in coef)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |==============================                                        |  44%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |=================================================================     |  94%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
#> 
#> -- propagate_error_dataframe summary ------------------
#>   Rows processed     : 200
#>   x_est finite       : 200
#>   x_est NA           : 0
#>   se_x finite        : 200
#>   se_x NA            : 0  (grad or var_par issue)
#>   cv_x < cap         : 70
#>   cv_x at cap (150) : 130
#>   grad_t NA rows     : 0
#>   var_par NA rows    : 0
#>   x_est range        : [1.0000, 5.4771]
#>   se_x  range        : [0.0697, 2.5278]
#>   cv_x  range (excl cap): [16.04, 147.93]
#> -------------------------------------------------------
#> 
#> [cv_x diagnostic] --- standards pred_se ---
#>   cv_x_max (cap)   : 150.0
#>   N total          : 200
#>   N finite cv_x    : 200
#>   N non-finite raw : 0  (replaced with cap)
#>   Min  cv_x        : 16.043  at predicted_concentration = 3.3398
#>   Max  cv_x        : 150.000
#>   Mean cv_x        : 124.692
#>   N cv_x > 20      : 195
#>   N cv_x at cap    : 130
#>   [INFO] 130 point(s) capped at cv_x_max=150.0 (asymptote proximity or failed propagation).
#>     curve_id timeperiod patientid well stype sampleid agroup dilution
#> 101        6     month6   PAT_101   A3     X     a101 GroupB      400
#> 102        6     month3   PAT_102   B3     X     a102 GroupB      400
#> 103        6     month3   PAT_103   C3     X     a103 GroupB      400
#> 104        6   baseline   PAT_104   D3     X     a104 GroupA      400
#> 105        6     month6   PAT_105   E3     X     a105 GroupB      400
#> 106        6     month3   PAT_106   F3     X     a106 GroupA      400
#>     samplingerrors     od assay_response_variable assay_independent_variable
#> 101             NA 3.3980                      od              concentration
#> 102             NA 0.0708                      od              concentration
#> 103             NA 2.9100                      od              concentration
#> 104             NA 3.3102                      od              concentration
#> 105             NA 0.0515                      od              concentration
#> 106             NA 2.8245                      od              concentration
#> sample_se has 20 row(s)
#> 
#> === Propagation Input Diagnosis ===
#> Model         : loglogistic5 
#> coef(fit)     : a = -1.37773, b = 3, c = 3.48083, d = 0.61369, g = 1.70296 
#> vcov dim      : 5 x 5 
#> vcov rownames : a, b, c, d, g 
#> fixed_a       : NULL 
#> ===================================
#> [propagate] CV formula   : LINEAR-scale (se_x * ln(10) * 100) --- avoids /0 at log10(conc)=0
#> [propagate] Model       : loglogistic5
#> [propagate] Free params : a, b, c, d, g
#> [propagate] Sigma rows  : a, b, c, d, g
#> [propagate] fixed_a     : NULL (a is free in coef)
#>   |                                                                              |                                                                      |   0%  |                                                                              |====                                                                  |   5%  |                                                                              |=======                                                               |  10%  |                                                                              |==========                                                            |  15%  |                                                                              |==============                                                        |  20%  |                                                                              |==================                                                    |  25%  |                                                                              |=====================                                                 |  30%  |                                                                              |========================                                              |  35%  |                                                                              |============================                                          |  40%  |                                                                              |================================                                      |  45%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  55%  |                                                                              |==========================================                            |  60%  |                                                                              |==============================================                        |  65%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  75%  |                                                                              |========================================================              |  80%  |                                                                              |============================================================          |  85%  |                                                                              |===============================================================       |  90%  |                                                                              |==================================================================    |  95%  |                                                                              |======================================================================| 100%
#> 
#> -- propagate_error_dataframe summary ------------------
#>   Rows processed     : 20
#>   x_est finite       : 20
#>   x_est NA           : 0
#>   se_x finite        : 20
#>   se_x NA            : 0  (grad or var_par issue)
#>   cv_x < cap         : 5
#>   cv_x at cap (150) : 15
#>   grad_t NA rows     : 0
#>   var_par NA rows    : 0
#>   x_est range        : [0.8916, 4.5231]
#>   se_x  range        : [0.0701, 2.7577]
#>   cv_x  range (excl cap): [16.14, 77.49]
#> -------------------------------------------------------
#> 
#> [cv_x diagnostic] --- samples sample_se ---
#>   cv_x_max (cap)   : 150.0
#>   N total          : 20
#>   N finite cv_x    : 20
#>   N non-finite raw : 0  (replaced with cap)
#>   Min  cv_x        : 16.139  at predicted_concentration = 3.3332
#>   Max  cv_x        : 150.000
#>   Mean cv_x        : 121.936
#>   N cv_x > 20      : 18
#>   N cv_x at cap    : 15
#>   [INFO] 15 point(s) capped at cv_x_max=150.0 (asymptote proximity or failed propagation).
#> Finished predict_and_propagate_error

head(sc$best_fit$sample_se)
#>   curve_id timeperiod patientid well stype sampleid agroup samplingerrors
#> 1        6     month6   PAT_101   A3     X     a101 GroupB             NA
#> 2        6     month3   PAT_102   B3     X     a102 GroupB             NA
#> 3        6     month3   PAT_103   C3     X     a103 GroupB             NA
#> 4        6   baseline   PAT_104   D3     X     a104 GroupA             NA
#> 5        6     month6   PAT_105   E3     X     a105 GroupB             NA
#> 6        6     month3   PAT_106   F3     X     a106 GroupA             NA
#>   assay_response_variable assay_independent_variable raw_assay_response
#> 1                      od              concentration             3.3980
#> 2                      od              concentration             0.0708
#> 3                      od              concentration             2.9100
#> 4                      od              concentration             3.3102
#> 5                      od              concentration             0.0515
#> 6                      od              concentration             2.8245
#>   dilution         od raw_predicted_concentration se_concentration
#> 1      400  0.5312234                    4.523121        0.8813822
#> 2      400 -1.1499667                    2.435850        0.7808687
#> 3      400  0.4638930                    4.307874        0.7521044
#> 4      400  0.5198542                    4.477366        0.8564575
#> 5      400 -1.2881928                    1.899142        1.2922224
#> 6      400  0.4509416                    4.277023        0.7314766
#>   final_predicted_concentration pcov
#> 1                   13340770.41  150
#> 2                     109121.50  150
#> 3                    8127065.46  150
#> 4                   12006774.66  150
#> 5                      31710.44  150
#> 6                    7569776.63  150
```

## Plot the Final Standard Curve

`$plot()` renders the selected best-fit curve with QC annotations
including the inflection point, LOD, MDC, RDL, and LOQ bands.

``` r
sc$plot()
#> FORMATTED: OD 
#> [1] "model_name"              "yhat_response"          
#> [3] "predicted_concentration" "se_concentration"       
#> [5] "pcov"                    "pcov_threshold"
```

Log-scale axes can be toggled without refitting:

``` r
sc$plot(is_display_log_independent = FALSE, is_display_log_response = FALSE)
#> [1] "model_name"              "yhat_response"          
#> [3] "predicted_concentration" "se_concentration"       
#> [5] "pcov"                    "pcov_threshold"
```

## Extract Results

`$get_results()` returns a named list of the key output tables for
downstream analysis or export.

``` r
results <- sc$get_results()

# Tables available
names(results)
#> [1] "fit_summary"          "best_parameters"      "best_fit_summary"    
#> [4] "sample_se"            "best_pred"            "best_standard"       
#> [7] "candidate_parameters" "candidate_residuals"  "second_derivative"

# One-row QC + fit statistics summary
results$best_fit_summary
#>   curve_id iter status   model_name         a b        c         d        g
#> 1        6   10   TRUE loglogistic5 -1.377728 3 3.480831 0.6136932 1.702958
#>   inflect_x  inflect_y     llod      ulod    mindc    maxdc   minrdl   maxrdl
#> 1  3.658286 -0.0521839 -1.17468 0.4573644 2.369143 4.292039 2.675491 4.019506
#>       lloq     uloq    lloq_y    uloq_y dfresidual nobs rsquare_fit       aic
#> 1 2.980659 3.980714 -0.816233 0.2704457          5   10   0.9959069 -18.74458
#>         bic   loglik        mse        cv bkg_method is_log_response is_log_x
#> 1 -16.92907 15.37229 0.00541172 -17.67497    ignored            TRUE     TRUE
#>   apply_prozone
#> 1          TRUE
#>                                                             formula
#> 1 od ~ a + (d - a) * (1 + g * exp(-b * (concentration - c)))^(-1/g)

# Parameter estimates
results$best_parameters
#>   curve_id term      lower     upper   estimate  std_error  statistic
#> 1        6    a -4.6989700 0.5634348 -1.3777281 0.07898914 -17.441995
#> 2        6    b  0.0500000 3.0000000  3.0000000 1.12596146   2.664390
#> 3        6    c  0.1045757 6.3725455  3.4808309 0.15604791  22.306168
#> 4        6    d -1.0279847 1.1673056  0.6136932 0.06081453  10.091225
#> 5        6    g  0.3000000 7.0000000  1.7029577 1.34385081   1.267222
#>        p_value
#> 1 1.135398e-05
#> 2 4.464802e-02
#> 3 3.364127e-06
#> 4 1.636600e-04
#> 5 2.608899e-01

# prediction grid from the best fitting model.
head(results$best_pred)
#>                curve_id   model_name yhat_response predicted_concentration
#> loglogistic5.1        6 loglogistic5     -1.359308                1.000000
#> loglogistic5.2        6 loglogistic5     -1.358563                1.022498
#> loglogistic5.3        6 loglogistic5     -1.357789                1.044996
#> loglogistic5.4        6 loglogistic5     -1.356983                1.067494
#> loglogistic5.5        6 loglogistic5     -1.356145                1.089992
#> loglogistic5.6        6 loglogistic5     -1.355272                1.112490
#>                se_concentration pcov
#> loglogistic5.1         2.527817  150
#> loglogistic5.2         2.483678  150
#> loglogistic5.3         2.440643  150
#> loglogistic5.4         2.398663  150
#> loglogistic5.5         2.357692  150
#> loglogistic5.6         2.317686  150

# Standards used for fitting
head(results$best_standard)
#>    curve_id stype sampleid well    dilution         od assay_response_variable
#> 51        6     S   STD_01   A1 1000.000000 -1.4271284                      od
#> 52        6     S   STD_02   B1  333.333333 -1.3353580                      od
#> 53        6     S   STD_03   C1  100.000000 -1.1944991                      od
#> 54        6     S   STD_04   D1   33.333333 -1.0793550                      od
#> 55        6     S   STD_05   E1   10.000000 -0.8823973                      od
#> 56        6     S   STD_06   F1    3.333333 -0.2783542                      od
#>    assay_independent_variable concentration
#> 51              concentration      1.000000
#> 52              concentration      1.477121
#> 53              concentration      2.000000
#> 54              concentration      2.477121
#> 55              concentration      3.000000
#> 56              concentration      3.477121
# Precision profile (populated after $propagate_error())
head(results$sample_se)
#>   curve_id timeperiod patientid well stype sampleid agroup samplingerrors
#> 1        6     month6   PAT_101   A3     X     a101 GroupB             NA
#> 2        6     month3   PAT_102   B3     X     a102 GroupB             NA
#> 3        6     month3   PAT_103   C3     X     a103 GroupB             NA
#> 4        6   baseline   PAT_104   D3     X     a104 GroupA             NA
#> 5        6     month6   PAT_105   E3     X     a105 GroupB             NA
#> 6        6     month3   PAT_106   F3     X     a106 GroupA             NA
#>   assay_response_variable assay_independent_variable raw_assay_response
#> 1                      od              concentration             3.3980
#> 2                      od              concentration             0.0708
#> 3                      od              concentration             2.9100
#> 4                      od              concentration             3.3102
#> 5                      od              concentration             0.0515
#> 6                      od              concentration             2.8245
#>   dilution         od raw_predicted_concentration se_concentration
#> 1      400  0.5312234                    4.523121        0.8813822
#> 2      400 -1.1499667                    2.435850        0.7808687
#> 3      400  0.4638930                    4.307874        0.7521044
#> 4      400  0.5198542                    4.477366        0.8564575
#> 5      400 -1.2881928                    1.899142        1.2922224
#> 6      400  0.4509416                    4.277023        0.7314766
#>   final_predicted_concentration pcov
#> 1                   13340770.41  150
#> 2                     109121.50  150
#> 3                    8127065.46  150
#> 4                   12006774.66  150
#> 5                      31710.44  150
#> 6                    7569776.63  150

# Model Comparisons can be derived from the following tables. 
# This is all converged candidate fits and not solely the best fit
head(results$candidate_parameters)
#>   curve_id        model converged parameter    estimate    conf.low   conf.high
#> 1        6    logistic5      TRUE       AIC  -0.8853559  -0.8853559  -0.8853559
#> 2        6 loglogistic5      TRUE       AIC -18.7445845 -18.7445845 -18.7445845
#> 3        6    logistic4      TRUE       AIC -17.4268983 -17.4268983 -17.4268983
#> 4        6 loglogistic4      TRUE       AIC  20.6983973  20.6983973  20.6983973
#> 5        6    gompertz4      TRUE       AIC -10.1839824 -10.1839824 -10.1839824
#> 6        6    logistic5      TRUE         a  -1.3854114  -1.8411760  -0.9296467
#>   is_best_model
#> 1         FALSE
#> 2          TRUE
#> 3         FALSE
#> 4         FALSE
#> 5         FALSE
#> 6         FALSE

# residuals vs fitted values
head(results$candidate_residuals)
#>             curve_id     model     fitted   residuals is_best_model
#> logistic5.1        6 logistic5 -1.3747686 -0.05235977         FALSE
#> logistic5.2        6 logistic5 -1.3523157  0.01695772         FALSE
#> logistic5.3        6 logistic5 -1.2758631  0.08136392         FALSE
#> logistic5.4        6 logistic5 -1.0934427  0.01408769         FALSE
#> logistic5.5        6 logistic5 -0.7136483 -0.16874898         FALSE
#> logistic5.6        6 logistic5 -0.2817147  0.00336045         FALSE

## The second derivative for the curve is also available 
head(results$second_derivative)
#>                curve_id        model        x      d2x_y
#> loglogistic5.1        6 loglogistic5 1.000000 0.05709433
#> loglogistic5.2        6 loglogistic5 1.022498 0.05939471
#> loglogistic5.3        6 loglogistic5 1.044996 0.06178835
#> loglogistic5.4        6 loglogistic5 1.067494 0.06427969
#> loglogistic5.5        6 loglogistic5 1.089992 0.06686873
#> loglogistic5.6        6 loglogistic5 1.112490 0.06956435
```

## A more advanced example: Fitting across Multiple Plates

Since `elisa_assay_data` contains one antigen each across six plates, we
can iterate over all combinations. Call `filter_by_curve_id()` with a
new `curve_id` from the lookup table and pass the result to
`$set_data()` — it replaces the data and clears all downstream state
automatically.

``` r
# Define per-antigen constraints as a named list
antigen_constraints_list <- list(
  alpha = data.frame(
    antigen                      = "alpha",
    l_asy_min_constraint         = 0,
    l_asy_max_constraint         = 0,
    l_asy_constraint_method      = "default",
    standard_curve_concentration = 10000,
    pcov_threshold               = 15,
    stringsAsFactors             = FALSE
  )
)

all_summaries <- vector("list", nrow(elisa_assay_example$curve_id_lookup))

for (i in seq_len(nrow(elisa_assay_example$curve_id_lookup))) {

  cid        <- elisa_assay_example$curve_id_lookup$curve_id[i]
  antigen    <- elisa_assay_example$curve_id_lookup$antigen[i]   # pull antigen name
  curve_data <- filter_by_curve_id(elisa_assay_example, curve_id = cid)

  # Update the StandardCurve object with new data AND new constraints
  sc$set_data(curve_data)
  sc$antigen_constraints <- antigen_constraints_list[[antigen]]  # swap constraints

  sc$set_curve_settings()$fit()$summarize()
  all_summaries[[i]] <- sc$get_results()$best_fit_summary
}

combined_summary <- do.call(rbind, all_summaries)
combined_summary
```

Key Notes:

- Pull the antigen name from the lookup table on each iteration so you
  know which constraint set to apply. This assumes `curve_id_lookup` has
  an antigen column.

- Swap `sc$antigen_constraints` before `$set_curve_settings()` because
  `set_curve_settings()` reads from antigen_constraints to configure the
  model bounds. If you swap after, the new constraints will not take
  effect until the next call.

- `$set_data()` clears downstream state (fitted models, summaries, etc.)
  but it does not reset antigen_constraints, so the swap step is still
  necessary and will persist safely until overwritten on the next
  iteration.

- Alternatively, if the constraint differences are substantial or you
  want a cleaner separation, instantiate a fresh `StandardCurve` object
  per antigen inside the loop using `StandardCurve$new()` with the
  appropriate constraint data frame. This avoids any risk of stale state
  but is slightly more verbose.

## Quick Reference

``` r
# # Full single-plate workflow using the built-in example data
# data(bead_assay_example)
# 
# curve_data <- filter_by_curve_id(bead_assay_example, curve_id = 1)
# 
# sc <- StandardCurve$new(curve_data, study_params, antigen_constraints,
#                         model_names = model_names,
#                         is_display_log_response    = TRUE,
#                         is_display_log_independent = TRUE)
# 
# sc$set_curve_settings()
# sc$fit()
# sc$summarize()
# sc$compare_models()
# sc$plot()
# sc$propagate_error()
# 
# results <- sc$get_results()
# 
# # Check status at any time
# sc   # or print(sc)
```
