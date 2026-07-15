# ==========================================================================
# SENSITIVITY ANALYSIS: Gc & Leaf Water Potential (ΨL)
# Evaluates kr_max, ks_max, kl_max, isohydricity, volumetric_capacitance,
# sapwood_theta_sat for 4 species × 2 conditions
#
# Usage (after all sensitivity runs complete and outputs are copied):
#   Rscript sensitivity_analysis_Gc_PSI_hydraulic.R
#
# Output:
#   Figures/lpj_guess_stem_storage/sensitivity_hydraulic/
#     best_parameters_Gc.csv          — best parameters per species for Gc
#     best_parameters_PSI.csv         — best parameters per species for ΨL
#     all_metrics_Gc.csv              — all metrics, all runs (Gc)
#     all_metrics_PSI.csv             — all metrics, all runs (ΨL)
#     parameter_sensitivity_Gc_*.png  — heatmaps & ranking/response plots (Gc)
#     parameter_sensitivity_PSI_*.png — heatmaps & ranking/response plots (ΨL)
#     combined_best_parameters.csv    — combined Gc+ΨL score ranking
# ==========================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(purrr)

# ==========================================================================
# 0. PATHS & SETUP
# ==========================================================================
BASE_DIR       <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD"
RESULTS_BASE   <- file.path(BASE_DIR, "results_lpj/results_sensitivity")
OUTPUT_DIR     <- file.path(BASE_DIR, "Figures/lpj_guess_stem_storage/sensitivity_hydraulic")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

setwd(BASE_DIR)

# ==========================================================================
# 1. DEFINE PARAMETER RANGES — COMMON OAT VALUES FOR ALL SPECIES
# ==========================================================================
# All four species share the same OAT test values (matching submit_all_sensitivity_hydraulic.sh).
# Baseline values are species-specific (from .ins files).
# volumetric_capacitance and sapwood_theta_sat baselines differ between control and drought.

parameter_bases <- list(
  Oak = list(
    control = list(kr = 1.8, ks = 2.1, kl = 10,   iso = 0.4, vc = 200, theta = 0.55),
    drought = list(kr = 1.8, ks = 2.1, kl = 10,   iso = 0.4, vc = 150, theta = 0.65)
  ),
  Beech = list(
    control = list(kr = 0.8, ks = 1.0, kl = 25,   iso = 0.6, vc = 200, theta = 0.55),
    drought = list(kr = 0.8, ks = 1.0, kl = 25,   iso = 0.6, vc = 150, theta = 0.70)
  ),
  Spruce = list(
    control = list(kr = 0.5, ks = 0.7, kl = 33.1, iso = 0.7, vc = 300, theta = 0.55),
    drought = list(kr = 0.5, ks = 0.7, kl = 33.1, iso = 0.7, vc = 300, theta = 0.55)
  ),
  Pine = list(
    control = list(kr = 0.6, ks = 0.8, kl = 12.5, iso = 0.7, vc = 300, theta = 0.60),
    drought = list(kr = 0.6, ks = 0.8, kl = 12.5, iso = 0.7, vc = 300, theta = 0.60)
  )
)

# Common OAT test values (matching shell script exactly)
oat_kr    <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5)
oat_ks    <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9)
oat_kl    <- c(5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35)
oat_iso   <- c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
oat_vc    <- c(150, 200, 250, 300, 350, 400)
oat_theta <- c(0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70)

# ==========================================================================
# 2. METRIC COMPUTATION FUNCTIONS
# ==========================================================================

safe_cor <- function(x, y) {
  if (length(x) < 3 || sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0)
    return(NA_real_)
  cor(x, y, use = "complete.obs")
}

compute_kge <- function(sim, obs) {
  # Kling-Gupta Efficiency
  # KGE = 1 - sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2)
  valid <- complete.cases(sim, obs)
  if (sum(valid) < 3) return(NA_real_)
  sim_v <- sim[valid]
  obs_v <- obs[valid]
  if (sd(obs_v) == 0 || mean(obs_v) == 0) return(NA_real_)

  r     <- cor(sim_v, obs_v)
  alpha <- sd(sim_v) / sd(obs_v)
  beta  <- mean(sim_v) / mean(obs_v)
  1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
}

compute_metrics <- function(sim, obs) {
  valid <- complete.cases(sim, obs)
  if (sum(valid) < 3) {
    return(data.frame(
      n = sum(valid), pearson_r = NA_real_, r_squared = NA_real_,
      rmse = NA_real_, nrmse_pct = NA_real_, bias = NA_real_,
      slope = NA_real_, kge = NA_real_
    ))
  }
  s <- sim[valid]
  o <- obs[valid]

  r_val   <- cor(s, o)
  fit     <- lm(s ~ o)
  sl      <- coef(fit)[2]
  rmse_v  <- sqrt(mean((s - o)^2))
  nrmse_v <- (rmse_v / mean(o)) * 100
  bias_v  <- mean(s - o)
  kge_v   <- compute_kge(s, o)

  data.frame(
    n          = sum(valid),
    pearson_r  = r_val,
    r_squared  = r_val^2,
    rmse       = rmse_v,
    nrmse_pct  = nrmse_v,
    bias       = bias_v,
    slope      = sl,
    kge        = kge_v
  )
}

# ==========================================================================
# 3. READ LPJ OUTPUT FOR A SINGLE RUN
# ==========================================================================
# Reads dgc.out (Gc) and dpsileaf.out (psi_leaf) directly from a run directory

read_lpj_run <- function(run_dir, species_pft) {
  # species_pft: e.g. "Que_rob", "Fag_syl", "Pic_abi", "Pin_syl"
  dgc_file      <- file.path(run_dir, "dgc.out")
  dpsileaf_file <- file.path(run_dir, "dpsileaf.out")

  if (!file.exists(dgc_file) || !file.exists(dpsileaf_file)) {
    warning("Missing output files in: ", run_dir)
    return(NULL)
  }

  gc_data <- tryCatch(read.table(dgc_file, header = TRUE, check.names = FALSE),
                      error = function(e) NULL)
  pl_data <- tryCatch(read.table(dpsileaf_file, header = TRUE, check.names = FALSE),
                      error = function(e) NULL)

  if (is.null(gc_data) || is.null(pl_data)) return(NULL)

  # Trim whitespace from column names (LPJ-GUESS output may have padding)
  colnames(gc_data) <- trimws(colnames(gc_data))
  colnames(pl_data) <- trimws(colnames(pl_data))

  # Verify the PFT column exists
  if (!species_pft %in% colnames(gc_data)) {
    warning("Column '", species_pft, "' not found in dgc.out. Available: ",
            paste(colnames(gc_data), collapse = ", "))
    return(NULL)
  }
  if (!species_pft %in% colnames(pl_data)) {
    warning("Column '", species_pft, "' not found in dpsileaf.out. Available: ",
            paste(colnames(pl_data), collapse = ", "))
    return(NULL)
  }

  # Convert LPJ date format (Year, Day) to Date and rename PFT column
  # unname() needed because species_pft_map[sp] returns a named vector
  pft_name <- unname(species_pft)
  colnames(gc_data)[colnames(gc_data) == pft_name] <- "Gc"
  gc_data <- gc_data %>%
    mutate(date = as.Date(Day, origin = paste0(Year, "-01-01"))) %>%
    select(date, Gc)

  colnames(pl_data)[colnames(pl_data) == pft_name] <- "psi_leaf"
  pl_data <- pl_data %>%
    mutate(date = as.Date(Day, origin = paste0(Year, "-01-01"))) %>%
    select(date, psi_leaf)

  inner_join(gc_data, pl_data, by = "date") %>%
    filter(date >= as.Date("2023-01-01"))
}

# ==========================================================================
# 4. READ OBSERVED DATA
# ==========================================================================

# Gc observations
obs_gc <- tryCatch({
  read.csv(file.path(BASE_DIR, "SCCII/sap_flux_gc_daytime_Climate_filter.csv")) %>%
    mutate(
      date      = as.Date(date),
      species   = tolower(species),
      treatment = case_when(treatment == "control" ~ "control",
                            treatment == "treatment" ~ "drought")
    ) %>%
    filter(date >= as.Date("2023-01-01"), !is.na(G_ms), G_ms <= 12)
}, error = function(e) NULL)

# Leaf water potential observations (midday)
obs_psi <- tryCatch({
  bind_rows(
    read.csv(file.path(BASE_DIR, "SCCII/psiL_hoelstein_drought.csv")) %>%
      mutate(treatment = "drought"),
    read.csv(file.path(BASE_DIR, "SCCII/psiL_hoelstein_control.csv")) %>%
      mutate(treatment = "control")
  ) %>%
    mutate(date = as.Date(date), month = month(date)) %>%
    filter(month >= 6, month <= 9) %>%
    rename(species = species_name) %>%
    group_by(date, species, treatment) %>%
    summarise(psi_leaf_md_obs = mean(md_wp_av, na.rm = TRUE), .groups = "drop")
}, error = function(e) NULL)

# ==========================================================================
# 5. METEO DATA (for climate filter if needed)
# ==========================================================================
read_meteo <- function(treatment) {
  meteo_path <- ifelse(
    treatment == "control",
    file.path(BASE_DIR, "MeteoSwiss/MeteoSwiss_station/all_stations_RUE_replaced_daytime_control.csv"),
    file.path(BASE_DIR, "MeteoSwiss/MeteoSwiss_station/all_stations_RUE_replaced_daytime_drought.csv")
  )
  if (!file.exists(meteo_path)) return(NULL)
  read.csv(meteo_path) %>%
    filter(station_abbr == "RUE") %>%
    mutate(date = as.Date(date))
}

# ==========================================================================
# 6. SPECIES NAME MAPPING
# ==========================================================================
species_pft_map <- c(
  Oak    = "Que_rob",
  Beech  = "Fag_syl",
  Spruce = "Pic_abi",
  Pine   = "Pin_syl"
)

species_lower_map <- c(
  Oak    = "oak",
  Beech  = "beech",
  Spruce = "spruce",
  Pine   = "pine"
)

# ==========================================================================
# 7. MAIN EVALUATION LOOP
# ==========================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  SENSITIVITY ANALYSIS: Gc & Psi_Leaf Parameter Evaluation\n")
cat("  Parameters: kr_max, ks_max, kl_max, isohydricity, volumetric_capacitance, sapwood_theta_sat\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

all_metrics_gc  <- list()
all_metrics_psi <- list()

for (sp in names(parameter_bases)) {

  cat("\n--- Processing species:", sp, "---\n")

  pft    <- species_pft_map[sp]
  sp_low <- species_lower_map[sp]

  for (cond in c("control", "drought")) {

    cat("  Condition:", cond, "\n")

    # ---- Build parameter combinations ----
    # Helper: format tag values matching shell script
    fmt_tag_1dp <- function(x) {
      # For kr, ks, iso: 1 decimal place (e.g. "kr0.1", "iso-0.3")
      x <- round(x, digits = 8)
      ifelse(x == round(x), sprintf("%.1f", x), as.character(x))
    }
    fmt_tag_2dp <- function(x) {
      # For theta: 2 decimal places (e.g. "theta0.30")
      sprintf("%.2f", round(x, digits = 8))
    }

    bases <- parameter_bases[[sp]][[cond]]

    combos <- list()

    # 1) Baseline run (species×condition-specific values for all 6 params)
    combos[[1]] <- list(
      tag = "baseline",
      kr_max = bases$kr,
      ks_max = bases$ks,
      kl_max = bases$kl,
      isohydricity = bases$iso,
      volumetric_capacitance = bases$vc,
      sapwood_theta_sat = bases$theta
    )

    combo_idx <- 2

    # 2) kr_max OAT (common values, other params at species×condition base)
    for (kr in oat_kr) {
      combos[[combo_idx]] <- list(
        tag = paste0("kr", fmt_tag_1dp(kr)),
        kr_max = kr, ks_max = bases$ks,
        kl_max = bases$kl, isohydricity = bases$iso,
        volumetric_capacitance = bases$vc, sapwood_theta_sat = bases$theta
      )
      combo_idx <- combo_idx + 1
    }

    # 3) ks_max OAT
    for (ks in oat_ks) {
      combos[[combo_idx]] <- list(
        tag = paste0("ks", fmt_tag_1dp(ks)),
        kr_max = bases$kr, ks_max = ks,
        kl_max = bases$kl, isohydricity = bases$iso,
        volumetric_capacitance = bases$vc, sapwood_theta_sat = bases$theta
      )
      combo_idx <- combo_idx + 1
    }

    # 4) kl_max OAT
    for (kl in oat_kl) {
      combos[[combo_idx]] <- list(
        tag = paste0("kl", as.character(kl)),
        kr_max = bases$kr, ks_max = bases$ks,
        kl_max = kl, isohydricity = bases$iso,
        volumetric_capacitance = bases$vc, sapwood_theta_sat = bases$theta
      )
      combo_idx <- combo_idx + 1
    }

    # 5) isohydricity OAT
    for (iso in oat_iso) {
      combos[[combo_idx]] <- list(
        tag = paste0("iso", fmt_tag_1dp(iso)),
        kr_max = bases$kr, ks_max = bases$ks,
        kl_max = bases$kl, isohydricity = iso,
        volumetric_capacitance = bases$vc, sapwood_theta_sat = bases$theta
      )
      combo_idx <- combo_idx + 1
    }

    # 6) volumetric_capacitance OAT
    for (vc in oat_vc) {
      combos[[combo_idx]] <- list(
        tag = paste0("vc", as.character(vc)),
        kr_max = bases$kr, ks_max = bases$ks,
        kl_max = bases$kl, isohydricity = bases$iso,
        volumetric_capacitance = vc, sapwood_theta_sat = bases$theta
      )
      combo_idx <- combo_idx + 1
    }

    # 7) sapwood_theta_sat OAT
    for (theta in oat_theta) {
      combos[[combo_idx]] <- list(
        tag = paste0("theta", fmt_tag_2dp(theta)),
        kr_max = bases$kr, ks_max = bases$ks,
        kl_max = bases$kl, isohydricity = bases$iso,
        volumetric_capacitance = bases$vc, sapwood_theta_sat = theta
      )
      combo_idx <- combo_idx + 1
    }

    cat("    Evaluating", length(combos), "parameter combinations\n")

    for (cmb in combos) {
      # Build the run directory path (after copy_results_to_shared.sh)
      dir_name <- sprintf("%s_%s_%s", sp, cond, cmb$tag)
      run_dir  <- file.path(RESULTS_BASE, cond, sp, dir_name)

      # Check if run has outputs
      if (!dir.exists(run_dir)) {
        cat("      [SKIP] Directory not found:", dir_name, "\n")
        next
      }

      # Read LPJ output
      lpj_data <- read_lpj_run(run_dir, pft)
      if (is.null(lpj_data) || nrow(lpj_data) == 0) {
        cat("      [SKIP] No valid data in:", dir_name, "\n")
        next
      }

      # ---- Gc Evaluation ----
      if (!is.null(obs_gc)) {
        obs_sub <- obs_gc %>%
          filter(species == sp_low, treatment == cond) %>%
          group_by(date) %>%
          summarise(G_obs = mean(G_ms, na.rm = TRUE), .groups = "drop")

        joined_gc <- lpj_data %>%
          inner_join(obs_sub, by = "date") %>%
          filter(Gc <= 12)

        if (nrow(joined_gc) > 3) {
          m_gc <- compute_metrics(joined_gc$Gc, joined_gc$G_obs)
          m_gc$species      <- sp
          m_gc$condition    <- cond
          m_gc$run_tag      <- cmb$tag
          m_gc$kr_max       <- cmb$kr_max
          m_gc$ks_max       <- cmb$ks_max
          m_gc$kl_max       <- cmb$kl_max
          m_gc$isohydricity <- cmb$isohydricity
          m_gc$volumetric_capacitance <- cmb$volumetric_capacitance
          m_gc$sapwood_theta_sat      <- cmb$sapwood_theta_sat
          m_gc$variable     <- "Gc"
          all_metrics_gc[[length(all_metrics_gc) + 1]] <- m_gc

          cat(sprintf("      [OK] %s | Gc R²=%.3f KGE=%.3f\n",
                      cmb$tag, m_gc$r_squared, m_gc$kge))
        } else {
          cat("      [SKIP] Insufficient Gc data for:", cmb$tag, "\n")
        }
      }

      # ---- ΨL Evaluation (midday) ----
      if (!is.null(obs_psi)) {
        obs_psi_sub <- obs_psi %>%
          filter(species == sp, treatment == cond)

        joined_psi <- lpj_data %>%
          inner_join(obs_psi_sub, by = "date") %>%
          filter(!is.na(psi_leaf), !is.na(psi_leaf_md_obs))

        if (nrow(joined_psi) > 3) {
          m_psi <- compute_metrics(joined_psi$psi_leaf, joined_psi$psi_leaf_md_obs)
          m_psi$species      <- sp
          m_psi$condition    <- cond
          m_psi$run_tag      <- cmb$tag
          m_psi$kr_max       <- cmb$kr_max
          m_psi$ks_max       <- cmb$ks_max
          m_psi$kl_max       <- cmb$kl_max
          m_psi$isohydricity <- cmb$isohydricity
          m_psi$volumetric_capacitance <- cmb$volumetric_capacitance
          m_psi$sapwood_theta_sat      <- cmb$sapwood_theta_sat
          m_psi$variable     <- "psi_leaf"
          all_metrics_psi[[length(all_metrics_psi) + 1]] <- m_psi

          cat(sprintf("      [OK] %s | ΨL R²=%.3f KGE=%.3f\n",
                      cmb$tag, m_psi$r_squared, m_psi$kge))
        } else {
          cat("      [SKIP] Insufficient ΨL data for:", cmb$tag, "\n")
        }
      }
    }
  }
}

# ==========================================================================
# 8. COMPILE RESULTS
# ==========================================================================
metrics_gc  <- bind_rows(all_metrics_gc)
metrics_psi <- bind_rows(all_metrics_psi)

if (nrow(metrics_gc) == 0 && nrow(metrics_psi) == 0) {
  cat("\n*** No sensitivity results found. Have the runs completed and outputs been copied? ***\n")
  cat("Expected directory structure:\n")
  cat("  ", RESULTS_BASE, "/{condition}/{Species}/{Species}_{condition}_{tag}/\n")
  quit(save = "no", status = 1)
}

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  RESULTS SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# ==========================================================================
# 9. RANK PARAMETERS - ONE-AT-A-TIME ANALYSIS
# ==========================================================================

analyze_oat <- function(metrics_df, var_name) {
  # For each species×condition, identify best single-parameter value
  # by comparing OAT runs to baseline

  if (nrow(metrics_df) == 0) return(data.frame())

  all_params <- c("kr_max", "ks_max", "kl_max", "isohydricity",
                  "volumetric_capacitance", "sapwood_theta_sat")

  results <- list()

  for (sp in unique(metrics_df$species)) {
    for (cond in unique(metrics_df$condition)) {
      sub <- metrics_df %>% filter(species == sp, condition == cond)

      # Get baseline
      baseline <- sub %>% filter(run_tag == "baseline")
      if (nrow(baseline) == 0) next

      # For each parameter, evaluate which value gives best KGE/R²
      for (param in all_params) {
        # Get OAT runs for this parameter (tag starts with param abbreviation)
        param_abbrev <- case_when(
          param == "kr_max"                 ~ "kr",
          param == "ks_max"                 ~ "ks",
          param == "kl_max"                 ~ "kl",
          param == "isohydricity"           ~ "iso",
          param == "volumetric_capacitance" ~ "vc",
          param == "sapwood_theta_sat"      ~ "theta"
        )

        oat_runs <- sub %>%
          filter(grepl(paste0("^", param_abbrev, "[0-9.-]"), run_tag))

        if (nrow(oat_runs) == 0) next

        # Find best value based on composite score
        # Composite = 0.4*R² + 0.3*KGE + 0.3*(1 - |bias|/range)
        best <- oat_runs %>%
          mutate(
            composite = 0.4 * r_squared + 0.3 * kge +
              0.3 * (1 - pmin(abs(bias) / max(abs(bias), 0.001), 1))
          ) %>%
          slice_max(composite, n = 1, with_ties = FALSE)

        if (nrow(best) > 0) {
          results[[length(results) + 1]] <- data.frame(
            species             = sp,
            condition           = cond,
            parameter           = param,
            best_value          = best[[param]][1],
            baseline_value      = baseline[[param]][1],
            best_r_squared      = best$r_squared[1],
            best_kge            = best$kge[1],
            best_rmse           = best$rmse[1],
            best_slope          = best$slope[1],
            baseline_r_squared  = baseline$r_squared[1],
            baseline_kge        = baseline$kge[1],
            improvement_kge     = best$kge[1] - baseline$kge[1]
          )
        }
      }
    }
  }

  bind_rows(results)
}

oat_gc  <- analyze_oat(metrics_gc,  "Gc")
oat_psi <- analyze_oat(metrics_psi, "ΨL")

cat("\n--- Best Parameters by OAT Analysis (Gc) ---\n")
if (nrow(oat_gc) > 0) print(as.data.frame(oat_gc))

cat("\n--- Best Parameters by OAT Analysis (ΨL) ---\n")
if (nrow(oat_psi) > 0) print(as.data.frame(oat_psi))

# ==========================================================================
# 10. BEST OVERALL PARAMETER SET PER SPECIES
# ==========================================================================

select_best_params <- function(metrics_df) {
  # Select the single best run for each species×condition
  if (nrow(metrics_df) == 0) return(data.frame())

  metrics_df %>%
    group_by(species, condition) %>%
    mutate(
      # Composite score favoring high R², high KGE, slope near 1
      composite = 0.35 * r_squared +
        0.30 * kge +
        0.15 * (1 - pmin(abs(1 - slope), 1)) +
        0.10 * (1 - pmin(abs(bias), 1)) +
        0.10 * (1 - pmin(nrmse_pct / 200, 1))
    ) %>%
    slice_max(composite, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(species, condition, run_tag,
           kr_max, ks_max, kl_max, isohydricity,
           volumetric_capacitance, sapwood_theta_sat,
           r_squared, kge, rmse, nrmse_pct, bias, slope, composite)
}

best_gc  <- select_best_params(metrics_gc)
best_psi <- select_best_params(metrics_psi)

cat("\n--- Best Parameter Sets (Gc) ---\n")
if (nrow(best_gc) > 0) print(as.data.frame(best_gc))

cat("\n--- Best Parameter Sets (ΨL) ---\n")
if (nrow(best_psi) > 0) print(as.data.frame(best_psi))

# Combined recommendation (average of Gc and ΨL rankings)
if (nrow(metrics_gc) > 0 && nrow(metrics_psi) > 0) {
  # Score each run by both Gc and ΨL metrics
  combined_scores <- bind_rows(
    metrics_gc %>% mutate(source = "Gc"),
    metrics_psi %>% mutate(source = "PSI")
  ) %>%
    group_by(species, condition, run_tag, kr_max, ks_max, kl_max, isohydricity,
             volumetric_capacitance, sapwood_theta_sat) %>%
    summarise(
      n_sources         = n(),
      mean_r_squared    = mean(r_squared, na.rm = TRUE),
      mean_kge          = mean(kge, na.rm = TRUE),
      mean_slope        = mean(slope, na.rm = TRUE),
      mean_rmse         = mean(rmse, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      combined_score = 0.4 * mean_r_squared +
        0.4 * mean_kge +
        0.2 * (1 - pmin(abs(1 - mean_slope), 1))
    ) %>%
    group_by(species, condition) %>%
    slice_max(combined_score, n = 1, with_ties = FALSE) %>%
    ungroup()

  cat("\n--- COMBINED Best Parameters (Gc + ΨL) ---\n")
  print(as.data.frame(combined_scores))
}

# ==========================================================================
# 11. EXPORT TABLES
# ==========================================================================

write.csv(metrics_gc,  file.path(OUTPUT_DIR, "all_metrics_Gc.csv"),  row.names = FALSE)
write.csv(metrics_psi, file.path(OUTPUT_DIR, "all_metrics_PSI.csv"), row.names = FALSE)
if (nrow(oat_gc)  > 0) write.csv(oat_gc,  file.path(OUTPUT_DIR, "oat_best_Gc.csv"),  row.names = FALSE)
if (nrow(oat_psi) > 0) write.csv(oat_psi, file.path(OUTPUT_DIR, "oat_best_PSI.csv"), row.names = FALSE)
if (nrow(best_gc)  > 0) write.csv(best_gc,  file.path(OUTPUT_DIR, "best_parameters_Gc.csv"),  row.names = FALSE)
if (nrow(best_psi) > 0) write.csv(best_psi, file.path(OUTPUT_DIR, "best_parameters_PSI.csv"), row.names = FALSE)

if (exists("combined_scores") && nrow(combined_scores) > 0) {
  write.csv(combined_scores, file.path(OUTPUT_DIR, "combined_best_parameters.csv"), row.names = FALSE)
}

# ==========================================================================
# 12. DIAGNOSTIC PLOTS
# ==========================================================================

# --- 12a. Parameter Sensitivity Heatmap (R²) ---
create_heatmap <- function(metrics_df, var_label) {
  if (nrow(metrics_df) == 0) return(NULL)

  plot_data <- metrics_df %>%
    filter(!grepl("^combo", run_tag)) %>%  # OAT runs only
    mutate(
      param_changed = case_when(
        grepl("^kr", run_tag)    ~ "kr_max",
        grepl("^ks", run_tag)    ~ "ks_max",
        grepl("^kl", run_tag)    ~ "kl_max",
        grepl("^iso", run_tag)   ~ "isohydricity",
        grepl("^vc", run_tag)    ~ "volumetric_capacitance",
        grepl("^theta", run_tag) ~ "sapwood_theta_sat",
        run_tag == "baseline"    ~ "baseline"
      ),
      param_value = case_when(
        param_changed == "kr_max"                 ~ kr_max,
        param_changed == "ks_max"                 ~ ks_max,
        param_changed == "kl_max"                 ~ kl_max,
        param_changed == "isohydricity"           ~ isohydricity,
        param_changed == "volumetric_capacitance" ~ volumetric_capacitance,
        param_changed == "sapwood_theta_sat"      ~ sapwood_theta_sat,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(param_changed))

  p <- ggplot(plot_data, aes(x = factor(param_value), y = factor(species),
                             fill = r_squared)) +
    geom_tile(color = "white", linewidth = 0.5) +
    facet_grid(condition ~ param_changed, scales = "free_x") +
    scale_fill_viridis_c(name = "R²", option = "D", limits = c(0, 1)) +
    labs(
      title = paste("Parameter Sensitivity: R² for", var_label),
      subtitle = "One-at-a-time analysis. Baseline marked with *.",
      x = "Parameter value", y = "Species"
    ) +
    theme_minimal() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
      strip.text      = element_text(size = 10, face = "bold"),
      panel.grid      = element_blank(),
      plot.title      = element_text(hjust = 0.5),
      plot.subtitle   = element_text(hjust = 0.5, color = "grey40")
    )

  # Mark baseline
  baseline_data <- plot_data %>% filter(run_tag == "baseline")
  if (nrow(baseline_data) > 0) {
    p <- p + geom_point(
      data = baseline_data,
      aes(x = factor(param_value), y = species),
      shape = 8, size = 3, color = "red"
    )
  }

  p
}

# --- 12b. KGE Bar Chart per Species ---
create_kge_ranking <- function(metrics_df, var_label) {
  if (nrow(metrics_df) == 0) return(NULL)

  metrics_df %>%
    filter(!grepl("^combo", run_tag)) %>%
    group_by(species, condition, run_tag) %>%
    summarise(kge = mean(kge, na.rm = TRUE), r_squared = mean(r_squared, na.rm = TRUE),
              .groups = "drop") %>%
    ggplot(aes(x = reorder(run_tag, kge), y = kge, fill = species)) +
    geom_col() +
    facet_grid(species ~ condition, scales = "free_y") +
    coord_flip() +
    scale_fill_manual(values = c("Oak"="#E69F00", "Beech"="#0072B2",
                                 "Spruce"="#009E73", "Pine"="#F0E442")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = paste("KGE Ranking:", var_label),
      x = "Parameter combination", y = "KGE"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}

# --- 12c. Parameter Response Curve ---
create_response_curve <- function(metrics_df, var_label) {
  if (nrow(metrics_df) == 0) return(NULL)

  oat_data <- metrics_df %>%
    filter(!grepl("^combo", run_tag)) %>%
    mutate(
      param = case_when(
        grepl("^kr", run_tag)    ~ "kr_max",
        grepl("^ks", run_tag)    ~ "ks_max",
        grepl("^kl", run_tag)    ~ "kl_max",
        grepl("^iso", run_tag)   ~ "isohydricity",
        grepl("^vc", run_tag)    ~ "volumetric_capacitance",
        grepl("^theta", run_tag) ~ "sapwood_theta_sat",
        run_tag == "baseline"    ~ "baseline"
      ),
      pval = case_when(
        param == "kr_max"                 ~ kr_max,
        param == "ks_max"                 ~ ks_max,
        param == "kl_max"                 ~ kl_max,
        param == "isohydricity"           ~ isohydricity,
        param == "volumetric_capacitance" ~ volumetric_capacitance,
        param == "sapwood_theta_sat"      ~ sapwood_theta_sat,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(param), param != "baseline")

  ggplot(oat_data, aes(x = pval, y = r_squared, color = species)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_grid(condition ~ param, scales = "free_x") +
    scale_color_manual(values = c("Oak"="#E69F00", "Beech"="#0072B2",
                                  "Spruce"="#009E73", "Pine"="#F0E442")) +
    labs(
      title = paste("Parameter Response Curves:", var_label),
      x = "Parameter value", y = "R²",
      color = "Species"
    ) +
    theme_minimal() +
    theme(
      strip.text   = element_text(size = 9, face = "bold")
    )
}

# Generate all plots
if (nrow(metrics_gc) > 0) {
  p_gc_heat  <- create_heatmap(metrics_gc, "Canopy Conductance (Gc)")
  p_gc_rank  <- create_kge_ranking(metrics_gc, "Canopy Conductance (Gc)")
  p_gc_resp  <- create_response_curve(metrics_gc, "Canopy Conductance (Gc)")

  ggsave(file.path(OUTPUT_DIR, "parameter_sensitivity_Gc_heatmap.png"),
         p_gc_heat, width = 16, height = 8, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, "parameter_sensitivity_Gc_kge_ranking.png"),
         p_gc_rank, width = 12, height = 10, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, "parameter_response_curve_Gc.png"),
         p_gc_resp, width = 16, height = 9, dpi = 300)
}

if (nrow(metrics_psi) > 0) {
  p_psi_heat <- create_heatmap(metrics_psi, "Leaf Water Potential (ΨL)")
  p_psi_rank <- create_kge_ranking(metrics_psi, "Leaf Water Potential (ΨL)")
  p_psi_resp <- create_response_curve(metrics_psi, "Leaf Water Potential (ΨL)")

  ggsave(file.path(OUTPUT_DIR, "parameter_sensitivity_PSI_heatmap.png"),
         p_psi_heat, width = 16, height = 8, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, "parameter_sensitivity_PSI_kge_ranking.png"),
         p_psi_rank, width = 12, height = 10, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, "parameter_response_curve_PSI.png"),
         p_psi_resp, width = 16, height = 9, dpi = 300)
}

# ==========================================================================
# 13. FINAL RECOMMENDATIONS
# ==========================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  FINAL PARAMETER RECOMMENDATIONS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("To apply best parameters, update the .ins files with:\n\n")

if (exists("combined_scores") && nrow(combined_scores) > 0) {
  for (i in seq_len(nrow(combined_scores))) {
    r <- combined_scores[i, ]
    cat(sprintf("  %-8s %-10s: kr_max=%-5.1f ks_max=%-5.1f kl_max=%-5.1f isohydricity=%-5.1f vc=%-4.0f theta=%-5.2f (score=%.3f)\n",
                r$species, r$condition,
                r$kr_max, r$ks_max, r$kl_max, r$isohydricity,
                r$volumetric_capacitance, r$sapwood_theta_sat, r$combined_score))
  }
} else {
  # Fallback: show best Gc and PSI separately
  if (nrow(best_gc) > 0) {
    cat("--- Gc-based best parameters ---\n")
    for (i in seq_len(nrow(best_gc))) {
      r <- best_gc[i, ]
      cat(sprintf("  %-8s %-10s: kr_max=%-5.1f ks_max=%-5.1f kl_max=%-5.1f isohydricity=%-5.1f vc=%-4.0f theta=%-5.2f\n",
                  r$species, r$condition,
                  r$kr_max, r$ks_max, r$kl_max, r$isohydricity,
                  r$volumetric_capacitance, r$sapwood_theta_sat))
    }
  }
  if (nrow(best_psi) > 0) {
    cat("\n--- ΨL-based best parameters ---\n")
    for (i in seq_len(nrow(best_psi))) {
      r <- best_psi[i, ]
      cat(sprintf("  %-8s %-10s: kr_max=%-5.1f ks_max=%-5.1f kl_max=%-5.1f isohydricity=%-5.1f vc=%-4.0f theta=%-5.2f\n",
                  r$species, r$condition,
                  r$kr_max, r$ks_max, r$kl_max, r$isohydricity,
                  r$volumetric_capacitance, r$sapwood_theta_sat))
    }
  }
}

cat("\n*** Sensitivity analysis complete. Results saved to:", OUTPUT_DIR, "***\n")
