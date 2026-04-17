## code to prepare `bead_assay_data` dataset 
set.seed(42)

# ── True (population-level) parameters ───────────────────────────────────────
ALPHA <- list(family = "gompertz", a = 18,  d = 20000, b = 1.20, c = 0.10)
BETA  <- list(family = "5pl",      a = 25,  d = 28000, b = 1.80, c = 1.00, g = 0.80)
CONC   <- c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30)  # AU/mL
PLATES <- paste0("SDYexample:plate_", 1:3)

# ── Plate layout (96-well, column-major fill) ─────────────────────────────────
PLATE_GRID <- paste0(
  rep(LETTERS[1:8], times = 12),
  rep(1:12,         each  = 8)
)
WELLS_STANDARDS <- PLATE_GRID[1:10]
WELLS_SAMPLES   <- PLATE_GRID[17:36]
WELLS_BLANKS    <- c("G11", "H11", "G12", "H12")

stopifnot(length(intersect(WELLS_STANDARDS, WELLS_SAMPLES)) == 0)
stopifnot(length(intersect(WELLS_STANDARDS, WELLS_BLANKS))  == 0)
stopifnot(length(intersect(WELLS_SAMPLES,   WELLS_BLANKS))  == 0)
cat("Well zones verified: no overlaps\n")
cat("Standards:", WELLS_STANDARDS, "\n")
cat("Samples:  ", WELLS_SAMPLES,   "\n")
cat("Blanks:   ", WELLS_BLANKS,    "\n\n")

# ── Forward models ────────────────────────────────────────────────────────────
f_5pl <- function(x, a, d, b, c, g) {
  z <- exp(b * (log(x) - log(c)))
  d + (a - d) / (1 + z)^g
}

f_gompertz <- function(x, a, d, b, c) {
  u <- b * (log(x) - log(c))
  d + (a - d) * exp(-exp(pmin(u, 30)))
}

forward <- function(x, p) {
  if (p$family == "gompertz") f_gompertz(x, p$a, p$d, p$b, p$c)
  else                         f_5pl(x, p$a, p$d, p$b, p$c, p$g)
}

# ── Noise model ───────────────────────────────────────────────────────────────
add_noise <- function(mu, seed_offset = 0) {
  set.seed(42 + seed_offset)
  sigma <- 5 + 0.04 * abs(mu)
  pmax(1, rnorm(length(mu), mu, sigma))
}

# ── Standards ─────────────────────────────────────────────────────────────────
make_standards <- function(params, antigen_label, seed_base) {
  n <- length(CONC)
  stopifnot(n <= length(WELLS_STANDARDS))
  
  do.call(rbind, lapply(seq_along(PLATES), function(i) {
    set.seed(seed_base + i)
    d_plt <- params$d * exp(rnorm(1, 0, 0.04))
    c_plt <- params$c * exp(rnorm(1, 0, 0.05))
    p_plt <- params; p_plt$d <- d_plt; p_plt$c <- c_plt
    true_mfi <- forward(CONC, p_plt)
    obs_mfi  <- add_noise(true_mfi, seed_offset = seed_base * 10 + i)
    
    data.frame(
      curve_id                   = paste(antigen_label, PLATES[i], sep = ":"),
      stype                      = "S",
      sampleid                   = paste0("STD_", formatC(seq_len(n), width = 2, flag = "0")),
      well                       = WELLS_STANDARDS[seq_len(n)],
      dilution                   = 1 / CONC,
      mfi                        = round(obs_mfi, 1),
      assay_response_variable    = "mfi",
      assay_independent_variable = "concentration",
      stringsAsFactors = FALSE
    )
  }))
}

# ── Blanks ────────────────────────────────────────────────────────────────────
make_blanks <- function(params, antigen_label, n_per_plate = 4, seed_base = 200) {
  stopifnot(n_per_plate <= length(WELLS_BLANKS))
  
  do.call(rbind, lapply(seq_along(PLATES), function(i) {
    set.seed(seed_base + i)
    mfi_vals <- round(rnorm(n_per_plate, mean = params$a, sd = 2), 1)
    
    data.frame(
      curve_id                   = paste(antigen_label, PLATES[i], sep = ":"),
      stype                      = "B",
      well                       = WELLS_BLANKS[seq_len(n_per_plate)],
      dilution                   = 1L,
      mfi                        = mfi_vals,
      assay_response_variable    = "mfi",
      assay_independent_variable = "concentration",
      stringsAsFactors = FALSE
    )
  }))
}

# ── Samples ───────────────────────────────────────────────────────────────────
make_samples <- function(params, antigen_label,
                         n_per_plate    = 20,
                         seed_base      = 300,
                         timeperiods    = c("baseline", "month3", "month6"),
                         agroups        = c("GroupA", "GroupB")) {
  stopifnot(n_per_plate <= length(WELLS_SAMPLES))
  
  set.seed(seed_base)
  N           <- n_per_plate * length(PLATES)
  eff_conc    <- exp(runif(N, log(min(CONC) * 0.5), log(max(CONC) * 1.5)))
  patient_ids <- paste0("PAT_", formatC(seq_len(N), width = 3, flag = "0"))
  
  do.call(rbind, lapply(seq_along(PLATES), function(i) {
    idx <- ((i - 1) * n_per_plate + 1):(i * n_per_plate)
    set.seed(seed_base + i * 7)
    d_plt <- params$d * exp(rnorm(1, 0, 0.04))
    c_plt <- params$c * exp(rnorm(1, 0, 0.05))
    p_plt <- params; p_plt$d <- d_plt; p_plt$c <- c_plt
    true_mfi <- forward(eff_conc[idx], p_plt)
    obs_mfi  <- add_noise(true_mfi, seed_offset = seed_base * 100 + i)
    n        <- length(idx)
    dil      <- 2000L
    
    data.frame(
      curve_id                   = paste(antigen_label, PLATES[i], sep = ":"),
      timeperiod                 = sample(timeperiods, n, replace = TRUE),
      patientid                  = patient_ids[idx],
      well                       = WELLS_SAMPLES[seq_len(n)],
      stype                      = "X",
      sampleid                   = paste0(substr(antigen_label, 1, 1),
                                          formatC(idx, width = 3, flag = "0")),
      agroup                     = sample(agroups, n, replace = TRUE),
      dilution                   = dil,
      pctaggbeads                = round(runif(n, 0, 5), 2),
      samplingerrors             = NA_real_,
      mfi                        = round(obs_mfi, 1),
      effective_conc             = round(eff_conc[idx], 5),          # conc at detection
      true_serum_conc            = round(eff_conc[idx] * dil, 3),    # before dilution
      assay_response_variable    = "mfi",
      assay_independent_variable = "concentration",
      stringsAsFactors = FALSE
    )
  }))
}

# ── Generate all data ─────────────────────────────────────────────────────────
std_alpha   <- make_standards(ALPHA, "alpha", seed_base = 10)
std_beta    <- make_standards(BETA,  "beta",  seed_base = 20)
blank_alpha <- make_blanks(ALPHA, "alpha", seed_base = 200)
blank_beta  <- make_blanks(BETA,  "beta",  seed_base = 210)
samp_alpha  <- make_samples(ALPHA, "alpha", n_per_plate = 20, seed_base = 300)
samp_beta   <- make_samples(BETA,  "beta",  n_per_plate = 20, seed_base = 400)

# ── Assemble curve_data list ──────────────────────────────────────────────────
all_samples <- rbind(samp_alpha, samp_beta)

# samples: drop true conc cols
samples_clean <- all_samples
samples_clean$effective_conc  <- NULL
samples_clean$true_serum_conc <- NULL

# validation: keep true conc cols, rename to match generate_synthetic_data.R
validation <- all_samples
names(validation)[names(validation) == "effective_conc"]  <- "true_conc_at_detection_AU_mL"
names(validation)[names(validation) == "true_serum_conc"] <- "true_serum_conc_AU_mL"
validation <- validation[order(validation$curve_id,
                               validation$sampleid), ]

curve_data <- list(
  standards  = rbind(std_alpha, std_beta),
  blanks     = rbind(blank_alpha, blank_beta),
  samples    = samples_clean,
  validation = validation
)

# ── Replace curve_id string with integer ─────────────────────────────────────
curve_levels   <- unique(curve_data$standards$curve_id)
curve_id_lookup <- data.frame(
  curve_id  = seq_along(curve_levels),
  curve_str = curve_levels,
  stringsAsFactors = FALSE
)

curve_data <- lapply(curve_data, function(df) {
  df$curve_id <- as.integer(factor(df$curve_id, levels = curve_levels))
  df
})

curve_data <- c(curve_data, list(curve_id_lookup = curve_id_lookup))

# ── Verify column names ───────────────────────────────────────────────────────
cat("names(curve_data$standards)\n");  print(names(curve_data$standards))
cat("names(curve_data$blanks)\n");     print(names(curve_data$blanks))
cat("names(curve_data$samples)\n");    print(names(curve_data$samples))
cat("names(curve_data$validation)\n"); print(names(curve_data$validation))

# ── Verify well separation (within a single plate) ───────────────────────────
cat("\n── Well separation check ──\n")
one_plate_id <- curve_id_lookup$curve_int[curve_id_lookup$curve_id == 1,]
w_std <- curve_data$standards$well[curve_data$standards$curve_id == one_plate_id]
w_blk <- curve_data$blanks$well[curve_data$blanks$curve_id        == one_plate_id]
w_smp <- curve_data$samples$well[curve_data$samples$curve_id      == one_plate_id]
cat("Standard wells :", w_std, "\n")
cat("Blank wells    :", w_blk, "\n")
cat("Sample wells   :", w_smp, "\n")
cat("STD ∩ BLK:", intersect(w_std, w_blk), "→ should be empty\n")
cat("STD ∩ SMP:", intersect(w_std, w_smp), "→ should be empty\n")
cat("BLK ∩ SMP:", intersect(w_blk, w_smp), "→ should be empty\n")


curve_data
curve_id_elements <- c("antigen", "study_accession", "plate")
parts <- strsplit(curve_data$curve_id_lookup$curve_str, ":", fixed = TRUE)
curve_data$curve_id_lookup[, curve_id_elements] <- do.call(rbind, parts)

curve_data$curve_id_lookup$curve_str <- NULL
# do not save validation data for R package. 
curve_data_for_package <- curve_data

curve_data_for_package["validation"] <- NULL

response_var <- unique(curve_data_for_package$standards$assay_response_variable)
indep_var <- unique(curve_data_for_package$standards$assay_independent_variable)

bead_assay_example <- c(curve_data_for_package, list(response_var = response_var, indep_var = indep_var))

# usethis::use_data(bead_assay_data, overwrite = TRUE)
