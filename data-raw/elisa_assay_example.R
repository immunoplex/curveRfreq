## --------------------------------------------------------------
## ELISA synthetic data – 6 plates (3 × 5‑PL, 3 × Gompertz)
## --------------------------------------------------------------
set.seed(42)

## -----------------------------------------------------------------
## True (population‑level) parameters
## -----------------------------------------------------------------

## Gompertz (used on plates 4‑6)
alpha_params <- list(
  family = "gompertz",
  a      =  0.05,   # background OD (lower asymptote)
  d      =  3.20,   # saturation OD (upper asymptote)
  b      =  1.50,   # growth rate
  c      =  1.00    # EC50 in AU/mL
)

## 5‑PL  (plates 1‑3)
beta_params <- list(
  family = "5pl",
  a      =  0.08,   # background OD
  d      =  2.80,   # saturation OD
  b      =  2.00,   # slope
  c      =  0.50,   # EC50
  g      =  0.90    # asymmetry factor
)
## -----------------------------------------------------------------
##   Concentrations & plate layout (96‑well, column‑major)
## -----------------------------------------------------------------
CONC   <- c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0)   # AU/mL
PLATES <- paste0("SDYELISAexample:plate_", 1:6)   # six plates

PLATE_GRID <- paste0(
  rep(LETTERS[1:8], times = 12),
  rep(1:12, each = 8)
)

WELLS_STANDARDS <- PLATE_GRID[1:10]      # A1:B5  (10 wells)
WELLS_SAMPLES   <- PLATE_GRID[17:36]    # block of 20 wells per plate
WELLS_BLANKS    <- c("G11", "H11", "G12", "H12")

stopifnot(length(intersect(WELLS_STANDARDS, WELLS_SAMPLES)) == 0)
stopifnot(length(intersect(WELLS_STANDARDS, WELLS_BLANKS )) == 0)  # typo fixed
stopifnot(length(intersect(WELLS_SAMPLES,   WELLS_BLANKS )) == 0)

cat("Plate layout OK – no overlapping wells\n")

## -----------------------------------------------------------------
## 2  Forward models (5‑PL & Gompertz)
## -----------------------------------------------------------------
f_5pl <- function(x, a, d, b, c, g) {
  z <- exp(b * (log(x) - log(c)))
  d + (a - d) / (1 + z)^g
}

f_gompertz <- function(x, a, d, b, c) {
  u <- b * (log(x) - log(c))
  d + (a - d) * exp(-exp(pmin(u, 30)))   # cap u to avoid overflow
}

forward <- function(x, p) {
  if (p$family == "gompertz") {
    f_gompertz(x, p$a, p$d, p$b, p$c)
  } else {                               # default = 5‑PL
    f_5pl(x, p$a, p$d, p$b, p$c, p$g)
  }
}

## -----------------------------------------------------------------
## 3  Noise model (optical density)
## -----------------------------------------------------------------
add_noise <- function(mu, cv = 0.06, seed_offset = 0) {
  set.seed(42 + seed_offset)
  sigma <- cv * abs(mu) + 0.005
  pmin(4.0, pmax(0.001, rnorm(length(mu), mu, sigma)))
}

## -----------------------------------------------------------------
## 4  Plate → parameter mapping (first 3 = 5‑PL, last 3 = Gompertz)
## -----------------------------------------------------------------
params_by_plate <- list(
  beta_params, beta_params, beta_params,   # plates 1‑3  → 5‑PL
  alpha_params, alpha_params, alpha_params# plates 4‑6  → Gompertz
)

## -----------------------------------------------------------------
## 5  Standards (different curve family per plate)
## -----------------------------------------------------------------
make_standards <- function(params_list, antigen_label, seed_base = 10) {
  n <- length(CONC)
  stopifnot(n <= length(WELLS_STANDARDS))
  
  do.call(rbind,
          lapply(seq_along(PLATES), function(i) {
            set.seed(seed_base + i)
            p <- params_list[[i]]                     # <-- pick the right set
            
            ## plate‑level perturbations (upper asymptote & EC50)
            d_plt <- p$d * exp(rnorm(1, 0, 0.05))
            c_plt <- p$c * exp(rnorm(1, 0, 0.04))
            b_plt <- p$b * exp(rnorm(1, 0, 0.03))
            
            p_plt <- p
            p_plt$d <- d_plt
            p_plt$c <- c_plt
            p_plt$b <- b_plt
            
            true_od <- forward(CONC, p_plt)
            obs_od  <- add_noise(true_od,
                                 seed_offset = seed_base * 10 + i)
            
            data.frame(
              curve_id   = paste(antigen_label, PLATES[i], sep = ":"),
              stype      = "S",
              sampleid   = paste0("STD_",
                                  formatC(seq_len(n), width = 2, flag = "0")),
              well       = WELLS_STANDARDS[seq_len(n)],
              dilution   = 1 / CONC,
              od         = round(obs_od, 4),
              assay_response_variable    = "od",
              assay_independent_variable = "concentration",
              stringsAsFactors = FALSE
            )
          })
  )
}

## -----------------------------------------------------------------
## 6️  Blanks (same parameters for all plates – use α lower asymptote)
## -----------------------------------------------------------------
make_blanks <- function(params, antigen_label,
                        n_per_plate = 4, seed_base = 200) {
  stopifnot(n_per_plate <= length(WELLS_BLANKS))
  
  do.call(rbind,
          lapply(seq_along(PLATES), function(i) {
            set.seed(seed_base + i)
            od_vals <- round(
              pmax(0.001,
                   rnorm(n_per_plate,
                         mean = params$a,
                         sd = 0.005)),
              4)
            
            data.frame(
              curve_id   = paste(antigen_label, PLATES[i], sep = ":"),
              stype      = "B",
              well       = WELLS_BLANKS[seq_len(n_per_plate)],
              dilution   = 1L,
              od         = od_vals,
              assay_response_variable    = "od",
              assay_independent_variable = "concentration",
              stringsAsFactors = FALSE
            )
          })
  )
}

## -----------------------------------------------------------------
## 7️  Samples (all belong to the α antigen → use α parameters)
## -----------------------------------------------------------------
make_samples <- function(params, antigen_label,
                         n_per_plate = 20,
                         seed_base = 300,
                         timeperiods = c("baseline", "month3", "month6"),
                         agroups = c("GroupA", "GroupB")) {
  stopifnot(n_per_plate <= length(WELLS_SAMPLES))
  
  N <- n_per_plate * length(PLATES)
  
  ## Effective concentration AFTER the fixed serum dilution
  set.seed(seed_base)
  eff_conc <- exp(runif(N,
                        log(min(CONC) * 0.5),   # 0.5 × smallest standard
                        log(max(CONC) * 1.5)))  # 1.5 × largest standard
  
  patient_ids <- paste0("PAT_",
                        formatC(seq_len(N), width = 3, flag = "0"))
  
  do.call(rbind,
          lapply(seq_along(PLATES), function(i) {
            idx <- ((i - 1) * n_per_plate + 1):(i * n_per_plate)
            
            set.seed(seed_base + i * 7)
            
            ## plate‑level perturbations (same sources as standards)
            d_plt <- params$d * exp(rnorm(1, 0, 0.05))
            c_plt <- params$c * exp(rnorm(1, 0, 0.04))
            b_plt <- params$b * exp(rnorm(1, 0, 0.03))
            
            p_plt <- params
            p_plt$d <- d_plt
            p_plt$c <- c_plt
            p_plt$b <- b_plt
            
            true_od <- forward(eff_conc[idx], p_plt)
            obs_od  <- add_noise(true_od,
                                 seed_offset = seed_base * 100 + i)
            
            data.frame(
              curve_id   = paste(antigen_label, PLATES[i], sep = ":"),
              timeperiod = sample(timeperiods, n_per_plate, replace = TRUE),
              patientid  = patient_ids[idx],
              well       = WELLS_SAMPLES[seq_len(n_per_plate)],
              stype      = "X",
              sampleid   = paste0(substr(antigen_label, 1, 1),
                                  formatC(idx, width = 3, flag = "0")),
              agroup     = sample(agroups, n_per_plate, replace = TRUE),
              dilution   = 400L,                     # typical serum dilution
              samplingerrors = NA_real_,
              od         = round(obs_od, 4),
              effective_conc = round(eff_conc[idx], 5),
              true_serum_conc = round(eff_conc[idx] * 400, 3),
              assay_response_variable    = "od",
              assay_independent_variable = "concentration",
              stringsAsFactors = FALSE
            )
          })
  )
}

## -----------------------------------------------------------------
## 8️  Build the full data set
## -----------------------------------------------------------------
std_all   <- make_standards(params_by_plate, "alpha", seed_base = 10)
blank_all <- make_blanks(alpha_params, "alpha", seed_base = 200)   # blanks use α lower asymptote
samp_all  <- make_samples(alpha_params, "alpha",                 # samples use α (Gompertz) parameters
                          n_per_plate = 20, seed_base = 300)

## --------------------------------------------------------------
## 9️  Assemble the list that the package expects
## --------------------------------------------------------------
all_samples <- samp_all

samples_clean <- all_samples
samples_clean$effective_conc  <- NULL
samples_clean$true_serum_conc <- NULL

validation <- all_samples
names(validation)[names(validation) == "effective_conc"]  <- "true_conc_at_detection_AU_mL"
names(validation)[names(validation) == "true_serum_conc"] <- "true_serum_conc_AU_mL"
validation <- validation[order(validation$curve_id,
                               validation$sampleid), ]

curve_data <- list(
  standards  = std_all,
  blanks     = blank_all,
  samples    = samples_clean,
  validation = validation
)

## --------------------------------------------------------------
## 10️  Integer `curve_id` + lookup table
## --------------------------------------------------------------
curve_levels    <- unique(curve_data$standards$curve_id)
curve_id_lookup <- data.frame(
  curve_id  = seq_along(curve_levels),   # integer version
  curve_str = curve_levels,
  stringsAsFactors = FALSE
)

curve_data <- lapply(curve_data, function(df) {
  df$curve_id <- as.integer(factor(df$curve_id,
                                   levels = curve_levels))
  df
})
curve_data$curve_id_lookup <- curve_id_lookup

## --------------------------------------------------------------
## 11️  Split the compound ID into its three components
## --------------------------------------------------------------
parts <- strsplit(curve_data$curve_id_lookup$curve_str,
                  split = ":", fixed = TRUE)

curve_data$curve_id_lookup[, c("antigen", "study_accession", "plate")] <-
  do.call(rbind, parts)

curve_data$curve_id_lookup$curve_str <- NULL

## --------------------------------------------------------------
## 12️  Final package‑ready object
## --------------------------------------------------------------
curve_data_for_package <- curve_data
curve_data_for_package["validation"] <- NULL   # keep validation only for testing

response_var <- unique(curve_data_for_package$standards$assay_response_variable)
indep_var    <- unique(curve_data_for_package$standards$assay_independent_variable)

elisa_assay_example <- c(
  curve_data_for_package,
  list(response_var = response_var,
       indep_var    = indep_var)
)

# usethis::use_data(elisa_assay_example, overwrite = TRUE)