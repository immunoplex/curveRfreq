


# ============================================
# DB HELPERS FOR VIGNETTES (MINIMAL VERSION)
# ============================================

library(DBI)
# -------------------------------
# Connection
# -------------------------------
get_db_connection <- function() {
  dbConnect(
    RPostgres::Postgres(),
    dbname   = Sys.getenv("db"),
    host     = Sys.getenv("db_host"),
    port     = Sys.getenv("db_port"),
    user     = Sys.getenv("db_userid_x"),
    password = Sys.getenv("db_pwd_x"),
    sslmode  = "require",
    options  = "-c search_path=madi_results"
  )
}

# -------------------------------
# Minimal data pull
# -------------------------------
pull_data <- function(study_accession, experiment_accession, project_id, conn = conn) {
  plates <- fetch_db_header(study_accession = study_accession,
                            experiment_accession = experiment_accession,
                            project_id = project_id,
                            conn = conn)
  plates$plate_id <- trimws(plates$plate_id)

  plates$plate_nom <- paste(plates$plate, plates$nominal_sample_dilution, sep = "-")


  antigen_constraints <- fetch_antigen_parameters(study_accession = study_accession,
                                                  experiment_accession = experiment_accession,
                                                  project_id = project_id,
                                                  conn = conn)

  standard_curve_data <- fetch_db_standards(study_accession = study_accession,
                                            experiment_accession = experiment_accession,
                                            project_id = project_id,
                                            conn = conn)
  standard_curve_data$plate_id <- trimws(standard_curve_data$plate_id)

  # cat("plates\n")
  # print(names(plates))
  # cat("standard_curve data\n")
  # print(names(standard_curve_data))

  standards <- inner_join(standard_curve_data, plates[, c("study_accession", "experiment_accession" ,"plateid", "plate", "plate_id" , "assay_response_variable" ,"assay_independent_variable"
                                                          ,"project_id")], by = c("study_accession", "experiment_accession","plate_id"))[ ,c("project_id", "study_accession","experiment_accession","feature", "source", "wavelength","plateid",
                                                                                                                                             "plate", "stype", "nominal_sample_dilution",
                                                                                                                                             "sampleid","well","dilution","antigen","mfi",
                                                                                                                                             "assay_response_variable", "assay_independent_variable")]
  if (nrow(standards) == 0) {
    warning(paste("[pull_data] No standards found after join for study:",
                  study_accession, "experiment:", experiment_accession,
                  "- returning NULL"))
    return(NULL)
  }
  standards$plate_nom <- paste(standards$plate, standards$nominal_sample_dilution, sep = "-")

  blanks <- inner_join(fetch_db_buffer(study_accession = study_accession,
                                       experiment_accession = experiment_accession,
                                       project_id = project_id,
                                       conn = conn) %>% dplyr::mutate(plate_id = trimws(as.character(plate_id))),
                       plates[, c("study_accession", "experiment_accession" ,"plateid", "plate", "plate_id" , "assay_response_variable" ,"assay_independent_variable"
                                  ,"project_id")], by = c("study_accession", "experiment_accession","plate_id"))

  blanks$plate_nom <- paste(blanks$plate, blanks$nominal_sample_dilution, sep = "-")

  samples <- inner_join(fetch_db_samples(study_accession = study_accession,
                                         experiment_accession = experiment_accession,
                                         project_id = project_id,
                                         conn = conn) %>% dplyr::mutate(plate_id = trimws(as.character(plate_id))),
                        plates[, c("study_accession", "experiment_accession" ,"plateid", "plate", "plate_id" , "assay_response_variable" ,"assay_independent_variable"
                                   ,"project_id")], by = c("study_accession", "experiment_accession","plate_id"))

  samples$plate_nom <- paste(samples$plate, samples$nominal_sample_dilution, sep = "-")

  response_var = unique(plates$assay_response_variable)
  indep_var = unique(plates$assay_independent_variable)

  # ====================================================================
  # RESPONSE VARIABLE COLUMN ALIGNMENT
  # ====================================================================
  # The DB always stores the measurement in 'antibody_mfi' (aliased as 'mfi').
  # For bead arrays, assay_response_variable = "mfi" → column name matches.
  # For ELISA, assay_response_variable = "absorbance" → column name mismatch.
  # Rename the data column so downstream code can use response_var directly.
  if (length(response_var) == 1 && response_var != "mfi") {
    if ("mfi" %in% names(standards)) names(standards)[names(standards) == "mfi"] <- response_var
    if ("mfi" %in% names(blanks))    names(blanks)[names(blanks) == "mfi"] <- response_var
    if ("mfi" %in% names(samples))   names(samples)[names(samples) == "mfi"] <- response_var
    cat("  ✓ Renamed 'mfi' column to '", response_var, "' for ELISA compatibility\n", sep = "")
  }

  # ====================================================================
  # WAVELENGTH DETECTION AND source_nom CONSTRUCTION
  # ====================================================================
  # For ELISA, wavelength identifies distinct measurement channels.
  # For bead array, wavelength is not applicable.
  # Try to get wavelength from the data; if not available, derive from context.

  # Check if wavelength column exists in any of the data frames
  has_wavelength_std <- "wavelength" %in% names(standards)
  has_wavelength_blk <- "wavelength" %in% names(blanks)
  has_wavelength_smp <- "wavelength" %in% names(samples)

  # Add wavelength column if not present
  if (!has_wavelength_std) standards$wavelength <- WL_NONE
  if (!has_wavelength_blk) blanks$wavelength    <- WL_NONE
  if (!has_wavelength_smp) samples$wavelength   <- WL_NONE

  # Normalize any NA/empty wavelength values from DB to sentinel
  standards$wavelength <- normalize_wavelength(standards$wavelength)
  blanks$wavelength    <- normalize_wavelength(blanks$wavelength)
  samples$wavelength   <- normalize_wavelength(samples$wavelength)

  # Construct source_nom: combines source with wavelength for ELISA
  # For bead array (wavelength == WL_NONE), source_nom = source
  build_source_nom <- function(df) {
    src <- if ("source" %in% names(df)) as.character(df$source) else rep(NA_character_, nrow(df))
    src[is.na(src) | trimws(src) == ""] <- "unknown"
    wl <- as.character(df$wavelength)
    ifelse(
      is.na(wl) | trimws(wl) == "" | wl == WL_NONE,
      src,
      paste0(src, "|", wl, "_nm")
    )
  }

  standards$source_nom <- build_source_nom(standards)
  blanks$source_nom    <- build_source_nom(blanks)
  samples$source_nom   <- build_source_nom(samples)

  cat("  source_nom values (standards):", paste(unique(standards$source_nom), collapse = ", "), "\n")

  # ====================================================================
  # CHECK IF SOURCE PREFIXES DIFFER FROM STANDARDS
  # ====================================================================

  std_prefix <- unique(sub("\\|.*$", "", standards$source_nom))

  if (length(std_prefix) == 1) {

    sample_prefix <- unique(sub("\\|.*$", "", samples$source_nom))
    blank_prefix  <- unique(sub("\\|.*$", "", blanks$source_nom))

    needs_fix <- !(all(sample_prefix == std_prefix) & all(blank_prefix == std_prefix))

    if (needs_fix) {


      samples <- fix_source_nom(samples, std_prefix)
      blanks  <- fix_source_nom(blanks, std_prefix)

      cat("✓ source_nom prefixes updated to match standards:", std_prefix, "\n")

    } else {
      cat("✓ source_nom prefixes already match standards\n")
    }

  } else {
    cat("⚠ Multiple standard source prefixes detected:",
        paste(std_prefix, collapse = ", "), "\n")
  }

  mcmc_samples <- fetch_best_sample_robust_concentrations(study_accession = study_accession,
                                                          experiment_accession = experiment_accession,
                                                          project_id = project_id,
                                                          conn = conn)

  #mcmc_samples_read <<- mcmc_samples

  mcmc_pred <- fetch_best_pred_robust_concentrations(study_accession = study_accession,
                                                     experiment_accession = experiment_accession,
                                                     project_id = project_id,
                                                     conn = conn)

  # ── Normalize MCMC data for ELISA wavelength compatibility ──────────
  # Apply the same wavelength normalization and source_nom construction
  # that standards/blanks/samples receive above, so that
  # select_antigen_plate can filter MCMC data consistently.
  if (nrow(mcmc_samples) > 0) {
    if (!"wavelength" %in% names(mcmc_samples)) {
      mcmc_samples$wavelength <- WL_NONE
    }
    mcmc_samples$wavelength <- normalize_wavelength(mcmc_samples$wavelength)
    mcmc_samples$source_nom <- build_source_nom(mcmc_samples)

    # Ensure plate_nom exists and is consistent
    if (all(c("plate", "nominal_sample_dilution") %in% names(mcmc_samples))) {
      mcmc_samples$plate_nom <- paste(
        mcmc_samples$plate, mcmc_samples$nominal_sample_dilution, sep = "-"
      )
    }

    # Fix source_nom prefix to match standards (same logic as above)
    if (length(std_prefix) == 1) {
      mcmc_samp_prefix <- unique(sub("\\|.*$", "", mcmc_samples$source_nom))
      if (!all(mcmc_samp_prefix == std_prefix)) {
        mcmc_samples <- fix_source_nom(mcmc_samples, std_prefix)
        cat("  ✓ mcmc_samples source_nom prefixes updated to match standards\n")
      }
    }

    # Rename response column for ELISA compatibility
    if (length(response_var) == 1 && response_var != "mfi" &&
        "mfi" %in% names(mcmc_samples)) {
      names(mcmc_samples)[names(mcmc_samples) == "mfi"] <- response_var
    }
  }

  if (nrow(mcmc_pred) > 0) {
    if (!"wavelength" %in% names(mcmc_pred)) {
      mcmc_pred$wavelength <- WL_NONE
    }
    mcmc_pred$wavelength <- normalize_wavelength(mcmc_pred$wavelength)
    mcmc_pred$source_nom <- build_source_nom(mcmc_pred)

    if (all(c("plate", "nominal_sample_dilution") %in% names(mcmc_pred))) {
      mcmc_pred$plate_nom <- paste(
        mcmc_pred$plate, mcmc_pred$nominal_sample_dilution, sep = "-"
      )
    }

    if (length(std_prefix) == 1) {
      mcmc_pred_prefix <- unique(sub("\\|.*$", "", mcmc_pred$source_nom))
      if (!all(mcmc_pred_prefix == std_prefix)) {
        mcmc_pred <- fix_source_nom(mcmc_pred, std_prefix)
        cat("  ✓ mcmc_pred source_nom prefixes updated to match standards\n")
      }
    }
  }

  #mcmc_samples_2 <<- mcmc_samples

  return(list(plates=plates, standards=standards,
              blanks=blanks, samples=samples,
              mcmc_samples = mcmc_samples,
              mcmc_pred = mcmc_pred,
              antigen_constraints=antigen_constraints,
              response_var = response_var,
              indep_var = indep_var)
  )
}
# -------------------------------
# Safe wrapper for vignette use
# -------------------------------
load_example_from_db <- function(study_accession,
                                 experiment_accession,
                                 project_id) {

  # if (Sys.getenv("RUN_DB_EXAMPLES") != "true") {
  #   message("DB examples disabled. Set RUN_DB_EXAMPLES=true to enable.")
  #   return(NULL)
  # }

  conn <- get_db_connection()
  on.exit(DBI::dbDisconnect(conn), add = TRUE)

  pull_curve_data(
    study_accession,
    experiment_accession,
    project_id,
    conn
  )
}


fetch_db_standards <- function(study_accession, experiment_accession, project_id, conn) {
  query <- glue("SELECT study_accession, experiment_accession, feature, plate_id, stype, source, wavelength, sampleid, well, dilution, antigen, antibody_mfi AS mfi, nominal_sample_dilution
  FROM madi_results.xmap_standard
WHERE project_id = {project_id}
AND study_accession = '{study_accession}'
AND experiment_accession = '{experiment_accession}'
")
  standard_df  <- dbGetQuery(conn, query)
  standard_df <- distinct(standard_df)
  return(standard_df)
}

fetch_db_buffer <- function(study_accession, experiment_accession, project_id, conn) {
  query <- glue("SELECT study_accession, experiment_accession, plate_id, stype, source, wavelength, well, antigen, dilution,
  feature, antibody_mfi AS mfi, nominal_sample_dilution
  FROM madi_results.xmap_buffer
WHERE project_id = {project_id}
AND study_accession = '{study_accession}'
AND experiment_accession = '{experiment_accession}'
")
  blank_data <- dbGetQuery(conn, query)
  blank_data <- distinct(blank_data)
  return(blank_data)
}

fetch_db_controls <- function(study_accession, experiment_accession, project_id, conn) {
  query <- glue("SELECT study_accession, experiment_accession, plate_id, well, stype, source, wavelength, dilution, pctaggbeads, samplingerrors, antigen, antibody_mfi as MFI, antibody_n
                    feature, project_id, plateid, nominal_sample_dilution, plate
                  	FROM madi_results.xmap_control
              WHERE project_id = {project_id}
              AND study_accession = '{study_accession}'
              AND experiment_accession = '{experiment_accession}';")

  control_data <- dbGetQuery(conn, query)
  control_data <- distinct(control_data)

  return(control_data)
}

fetch_db_samples <- function(study_accession, experiment_accession, project_id, conn) {
  query <- glue("SELECT study_accession,
experiment_accession, plate_id, timeperiod, patientid,
well, stype, source, wavelength, sampleid,  agroup, dilution, pctaggbeads, samplingerrors, antigen, antibody_mfi AS mfi,
antibody_n, nominal_sample_dilution, feature FROM madi_results.xmap_sample
WHERE project_id = {project_id}
AND study_accession = '{study_accession}'
AND experiment_accession = '{experiment_accession}'
")
  sample_data <- dbGetQuery(conn, query)
  sample_data <- distinct(sample_data)
  return(sample_data)
}

fix_source_nom <- function(df, std_prefix) {

  suffix <- sub("^[^|]*", "", df$source_nom)
  df$source_nom <- paste0(std_prefix, suffix)

  df
}


fetch_study_parameters <- function(study_accession, param_user, param_group = "standard_curve_options", project_id, conn) {
  query <- glue("
  SELECT study_accession, param_name, param_boolean_value, param_character_value
	FROM madi_results.xmap_study_config
  WHERE project_id = {project_id}
  AND study_accession = '{study_accession}'
  AND param_user = '{param_user}'
  AND param_group = '{param_group}';
")
  study_parameters <- dbGetQuery(conn, query)
  return(list(
    applyProzone = study_parameters[study_parameters$param_name=="applyProzone", "param_boolean_value"],
    blank_option = study_parameters[study_parameters$param_name=="blank_option", "param_character_value"],
    standard_source = study_parameters[study_parameters$param_name=="default_source", "param_character_value"],
    is_log_response = study_parameters[study_parameters$param_name=="is_log_mfi_axis", "param_boolean_value"],
    is_log_independent = TRUE,
    mean_mfi = study_parameters[study_parameters$param_name=="mean_mfi", "param_boolean_value"]
  ))
}

fetch_antigen_parameters <- function(study_accession, experiment_accession, project_id, conn) {
  query <- glue("
  SELECT
    xmap_antigen_family_id,
    study_accession,
    experiment_accession,
    antigen,
    l_asy_min_constraint,
    l_asy_max_constraint,
    l_asy_constraint_method,
    standard_curve_concentration,
    pcov_threshold
  FROM madi_results.xmap_antigen_family
  WHERE project_id = {project_id}
  AND study_accession = '{study_accession}'
  AND experiment_accession = '{experiment_accession}'
  AND l_asy_constraint_method IS NOT NULL;
")
  antigen_constraints <- dbGetQuery(conn, query)
  return(antigen_constraints=antigen_constraints)
}


fetch_db_header <- function(study_accession, experiment_accession, project_id, conn) {
  query <- glue("SELECT study_accession, experiment_accession, plateid, plate, nominal_sample_dilution,plate_id,
  assay_response_variable, assay_independent_variable, nominal_sample_dilution, project_id
  FROM madi_results.xmap_header
WHERE project_id = {project_id}
AND study_accession = '{study_accession}'
AND experiment_accession = '{experiment_accession}'
")
  header_data <- dbGetQuery(conn, query)
  header_data <- distinct(header_data)
  return(header_data)
}

fetch_db_header_experiments <- function(study_accession, conn, verbose = TRUE) {
  query <- glue("SELECT study_accession, experiment_accession, plateid, plate, nominal_sample_dilution, plate_id,
  assay_response_variable, assay_independent_variable
  FROM madi_results.xmap_header
WHERE study_accession = '{study_accession}'
")
  header_data <- dbGetQuery(conn, query)
  header_data <- distinct(header_data)
  return(header_data)
}

fetch_best_sample_robust_concentrations <- function(study_accession, experiment_accession, project_id, conn) {
  query <- glue("SELECT
    -- Sample identifiers
    bss.project_id,
    bss.study_accession,
    bss.experiment_accession,
    bss.patientid,
    bss.sampleid,
    bss.well,
    bss.stype,
    bss.agroup,
    bss.timeperiod,
    bss.antigen,
    bss.plateid,
    bss.plate,
    bss.nominal_sample_dilution,
    CONCAT(bss.plate, '-', bss.nominal_sample_dilution) AS plate_nom,
    -- Feature info
    bss.feature,
    bss.wavelength,
    bss.source,
    -- Robust concentration results
	bss.assay_response,
    bss.raw_robust_concentration,
    bss.final_robust_concentration,
    bss.se_robust_concentration,
    bss.pcov_robust_concentration,
    -- Key LOQ/LOD gating only
    bss.gate_class_loq,
    bss.gate_class_lod,
    -- Curve quality from glance
    bga.last_concentration_calc_method

FROM madi_results.best_sample_se_all bss
INNER JOIN madi_results.best_glance_all bga
    ON  bss.best_glance_all_id = bga.best_glance_all_id
    AND bss.study_accession    = bga.study_accession
    AND bss.plateid             = bga.plateid
    AND bss.antigen             = bga.antigen
    AND bss.feature             = bga.feature
	AND bss.source              = bga.source
    AND bss.wavelength          = bga.wavelength

WHERE bss.study_accession                = '{study_accession}'
  AND bss.experiment_accession               = '{experiment_accession}'
  AND bss.project_id                     = {project_id}
  AND bga.last_concentration_calc_method = 'mcmc_robust';")

  robust_sample_concentrations <- dbGetQuery(conn, query)
  return(robust_sample_concentrations)
}

fetch_best_pred_robust_concentrations <- function(study_accession, experiment_accession, project_id, conn) {
  query <- glue("SELECT
    bpa.project_id,
    bpa.study_accession,
    bpa.experiment_accession,
    bpa.x,
    bpa.yhat,
    bpa.pcov,
    bpa.antigen,
    bpa.plateid,
    bpa.plate,
    bpa.nominal_sample_dilution,
    CONCAT(bpa.plate, '-', bpa.nominal_sample_dilution) AS plate_nom,
    bpa.feature,
    bpa.wavelength,
    bpa.source,
    -- MCMC robust results on the dense grid
    bpa.raw_robust_concentration,
    bpa.se_robust_concentration,
    bpa.pcov_robust_concentration,
    -- Curve quality from glance
    bga.last_concentration_calc_method
FROM madi_results.best_pred_all bpa
INNER JOIN madi_results.best_glance_all bga
    ON  bpa.best_glance_all_id = bga.best_glance_all_id
WHERE bpa.study_accession        = '{study_accession}'
  AND bpa.experiment_accession   = '{experiment_accession}'
  AND bpa.project_id             = {project_id}
  AND bga.last_concentration_calc_method = 'mcmc_robust';")
  robust_pred_concentrations <- dbGetQuery(conn, query)
  return(robust_pred_concentrations)
}
