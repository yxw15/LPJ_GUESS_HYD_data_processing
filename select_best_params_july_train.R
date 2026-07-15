# ==========================================================================
# PARAMETER SELECTION ON JULY + SEPTEMBER (TRAINING DATA)
# Trains on July AND September (all years) from the sensitivity analysis runs,
# outputs best parameters per species × condition.
#
# These parameters are then used to run LPJ-GUESS and validated
# on August data (see run_pipeline_validation_august.R).
#
# Usage:
#   Rscript select_best_params_july_train.R
#
# Output:
#   Figures/lpj_guess_stem_storage/cross_validation/
#     best_params_trained_on_july_sept.csv  — best parameters per species × treatment
# ==========================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)

# ==========================================================================
# 0. PATHS & SETUP
# ==========================================================================
BASE_DIR       <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD"
RESULTS_BASE   <- file.path(BASE_DIR, "results_lpj/results_sensitivity")
OUTPUT_DIR     <- file.path(BASE_DIR, "Figures/lpj_guess_stem_storage/cross_validation")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(BASE_DIR)

TRAIN_MONTHS <- c(7, 9)  # July and September

# ==========================================================================
# 1. METRIC FUNCTIONS
# ==========================================================================
compute_kge <- function(sim, obs) {
  valid <- complete.cases(sim, obs)
  if (sum(valid) < 3) return(NA_real_)
  sim_v <- sim[valid]; obs_v <- obs[valid]
  if (sd(obs_v) == 0 || mean(obs_v) == 0) return(NA_real_)
  r <- cor(sim_v, obs_v)
  alpha <- sd(sim_v) / sd(obs_v)
  beta  <- mean(sim_v) / mean(obs_v)
  1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
}

compute_metrics <- function(sim, obs) {
  valid <- complete.cases(sim, obs)
  if (sum(valid) < 3) {
    return(data.frame(n = sum(valid), pearson_r = NA_real_, r_squared = NA_real_,
                      rmse = NA_real_, nrmse_pct = NA_real_, bias = NA_real_,
                      slope = NA_real_, kge = NA_real_))
  }
  s <- sim[valid]; o <- obs[valid]
  r_val <- cor(s, o)
  fit   <- lm(s ~ o)
  sl    <- coef(fit)[2]
  rmse_v <- sqrt(mean((s - o)^2))
  nrmse_v <- (rmse_v / mean(o)) * 100
  bias_v <- mean(s - o)
  data.frame(n = sum(valid), pearson_r = r_val, r_squared = r_val^2,
             rmse = rmse_v, nrmse_pct = nrmse_v, bias = bias_v,
             slope = sl, kge = compute_kge(s, o))
}

# ==========================================================================
# 2. READ OBSERVATIONS (ALL YEARS, JULY + SEPTEMBER ONLY)
# ==========================================================================

# Gc observations — July + September only
obs_gc_train <- tryCatch({
  read.csv(file.path(BASE_DIR, "SCCII/sap_flux_gc_daytime_Climate_filter.csv")) %>%
    mutate(
      date      = as.Date(date),
      species   = tolower(species),
      treatment = case_when(treatment == "control" ~ "control",
                            treatment == "treatment" ~ "drought")
    ) %>%
    filter(month(date) %in% TRAIN_MONTHS, !is.na(G_ms), G_ms <= 12)
}, error = function(e) NULL)

# PSI observations — July + September only
obs_psi_train <- tryCatch({
  bind_rows(
    read.csv(file.path(BASE_DIR, "SCCII/psiL_hoelstein_drought.csv")) %>%
      mutate(treatment = "drought"),
    read.csv(file.path(BASE_DIR, "SCCII/psiL_hoelstein_control.csv")) %>%
      mutate(treatment = "control")
  ) %>%
    mutate(date = as.Date(date), month = month(date)) %>%
    filter(month %in% TRAIN_MONTHS) %>%
    rename(species = species_name) %>%
    group_by(date, species, treatment) %>%
    summarise(psi_leaf_md_obs = mean(md_wp_av, na.rm = TRUE), .groups = "drop")
}, error = function(e) NULL)

# ==========================================================================
# 3. MAPPINGS & PARAMETER DEFINITIONS
# ==========================================================================
species_pft_map <- c(Oak = "Que_rob", Beech = "Fag_syl", Spruce = "Pic_abi", Pine = "Pin_syl")
species_lower_map <- c(Oak = "oak", Beech = "beech", Spruce = "spruce", Pine = "pine")

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

oat_kr    <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5)
oat_ks    <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9)
oat_kl    <- c(5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35)
oat_iso   <- c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
oat_vc    <- c(150, 200, 250, 300, 350, 400)
oat_theta <- c(0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70)

build_combos <- function(bases) {
  combos <- list()
  fmt_tag_1dp <- function(x) {
    x <- round(x, digits = 8)
    ifelse(x == round(x), sprintf("%.1f", x), as.character(x))
  }
  fmt_tag_2dp <- function(x) sprintf("%.2f", round(x, digits = 8))

  combos[[1]] <- list(tag = "baseline", kr_max = bases$kr, ks_max = bases$ks,
    kl_max = bases$kl, isohydricity = bases$iso,
    volumetric_capacitance = bases$vc, sapwood_theta_sat = bases$theta)
  idx <- 2
  for (v in oat_kr) {
    combos[[idx]] <- list(tag = paste0("kr", fmt_tag_1dp(v)),
      kr_max = v, ks_max = bases$ks, kl_max = bases$kl,
      isohydricity = bases$iso, volumetric_capacitance = bases$vc,
      sapwood_theta_sat = bases$theta); idx <- idx + 1
  }
  for (v in oat_ks) {
    combos[[idx]] <- list(tag = paste0("ks", fmt_tag_1dp(v)),
      kr_max = bases$kr, ks_max = v, kl_max = bases$kl,
      isohydricity = bases$iso, volumetric_capacitance = bases$vc,
      sapwood_theta_sat = bases$theta); idx <- idx + 1
  }
  for (v in oat_kl) {
    combos[[idx]] <- list(tag = paste0("kl", as.character(v)),
      kr_max = bases$kr, ks_max = bases$ks, kl_max = v,
      isohydricity = bases$iso, volumetric_capacitance = bases$vc,
      sapwood_theta_sat = bases$theta); idx <- idx + 1
  }
  for (v in oat_iso) {
    combos[[idx]] <- list(tag = paste0("iso", fmt_tag_1dp(v)),
      kr_max = bases$kr, ks_max = bases$ks, kl_max = bases$kl,
      isohydricity = v, volumetric_capacitance = bases$vc,
      sapwood_theta_sat = bases$theta); idx <- idx + 1
  }
  for (v in oat_vc) {
    combos[[idx]] <- list(tag = paste0("vc", as.character(v)),
      kr_max = bases$kr, ks_max = bases$ks, kl_max = bases$kl,
      isohydricity = bases$iso, volumetric_capacitance = v,
      sapwood_theta_sat = bases$theta); idx <- idx + 1
  }
  for (v in oat_theta) {
    combos[[idx]] <- list(tag = paste0("theta", fmt_tag_2dp(v)),
      kr_max = bases$kr, ks_max = bases$ks, kl_max = bases$kl,
      isohydricity = bases$iso, volumetric_capacitance = bases$vc,
      sapwood_theta_sat = v); idx <- idx + 1
  }
  combos
}

# ==========================================================================
# 4. READ LPJ OUTPUT (reuse from sensitivity analysis)
# ==========================================================================
read_lpj_run <- function(run_dir, species_pft) {
  dgc_file      <- file.path(run_dir, "dgc.out")
  dpsileaf_file <- file.path(run_dir, "dpsileaf.out")

  if (!file.exists(dgc_file) || !file.exists(dpsileaf_file)) return(NULL)

  gc_data <- tryCatch(read.table(dgc_file, header = TRUE, check.names = FALSE),
                      error = function(e) NULL)
  pl_data <- tryCatch(read.table(dpsileaf_file, header = TRUE, check.names = FALSE),
                      error = function(e) NULL)

  if (is.null(gc_data) || is.null(pl_data)) return(NULL)
  colnames(gc_data) <- trimws(colnames(gc_data))
  colnames(pl_data) <- trimws(colnames(pl_data))

  if (!species_pft %in% colnames(gc_data) || !species_pft %in% colnames(pl_data))
    return(NULL)

  pft_name <- unname(species_pft)
  colnames(gc_data)[colnames(gc_data) == pft_name] <- "Gc"
  gc_data <- gc_data %>%
    mutate(date = as.Date(Day, origin = paste0(Year, "-01-01"))) %>%
    select(date, Gc)

  colnames(pl_data)[colnames(pl_data) == pft_name] <- "psi_leaf"
  pl_data <- pl_data %>%
    mutate(date = as.Date(Day, origin = paste0(Year, "-01-01"))) %>%
    select(date, psi_leaf)

  inner_join(gc_data, pl_data, by = "date")
}

# ==========================================================================
# 5. EVALUATE ALL PARAMETER COMBINATIONS ON JULY + SEPTEMBER DATA
# ==========================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  PARAMETER SELECTION ON JULY + SEPTEMBER (TRAINING MONTHS)\n")
cat("  Months:", paste(TRAIN_MONTHS, collapse = ", "), "| All available years\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Summary of available data
cat("--- Training Observation Counts ---\n")
if (!is.null(obs_gc_train)) {
  gc_counts <- obs_gc_train %>%
    group_by(species, treatment) %>%
    summarise(n_days = n_distinct(date), n_obs = n(), .groups = "drop")
  cat("Gc (daily):\n")
  print(as.data.frame(gc_counts))
}
if (!is.null(obs_psi_train)) {
  psi_counts <- obs_psi_train %>%
    group_by(species, treatment) %>%
    summarise(n_days = n_distinct(date), .groups = "drop")
  cat("\nPSI (midday):\n")
  print(as.data.frame(psi_counts))
}
cat("\n")

all_results <- list()

for (sp in names(parameter_bases)) {
  cat("--- Processing:", sp, "---\n")
  pft      <- species_pft_map[sp]
  sp_low   <- species_lower_map[sp]

  for (cond in c("control", "drought")) {
    cat(sprintf("  Condition: %s\n", cond))

    bases  <- parameter_bases[[sp]][[cond]]
    combos <- build_combos(bases)

    best_score <- -Inf
    best_combo <- NULL

    for (cmb in combos) {
      dir_name <- sprintf("%s_%s_%s", sp, cond, cmb$tag)
      run_dir  <- file.path(RESULTS_BASE, cond, sp, dir_name)

      if (!dir.exists(run_dir)) next

      lpj_data <- read_lpj_run(run_dir, pft)
      if (is.null(lpj_data) || nrow(lpj_data) == 0) next

      # --- Gc evaluation on TRAINING MONTHS ---
      gc_score <- NA_real_
      gc_metrics <- NULL
      if (!is.null(obs_gc_train)) {
        obs_sub <- obs_gc_train %>%
          filter(species == sp_low, treatment == cond) %>%
          group_by(date) %>%
          summarise(G_obs = mean(G_ms, na.rm = TRUE), .groups = "drop")

        joined_gc <- lpj_data %>%
          filter(month(date) %in% TRAIN_MONTHS) %>%
          inner_join(obs_sub, by = "date") %>%
          filter(Gc <= 12)

        if (nrow(joined_gc) >= 5) {
          m <- compute_metrics(joined_gc$Gc, joined_gc$G_obs)
          gc_score <- 0.4 * m$r_squared + 0.3 * m$kge +
            0.15 * (1 - pmin(abs(1 - m$slope), 1)) +
            0.15 * (1 - pmin(abs(m$bias), 1))
          gc_metrics <- m
        }
      }

      # --- PSI evaluation on TRAINING MONTHS ---
      psi_score <- NA_real_
      psi_metrics <- NULL
      if (!is.null(obs_psi_train)) {
        obs_psi_sub <- obs_psi_train %>%
          filter(species == sp, treatment == cond)

        joined_psi <- lpj_data %>%
          filter(month(date) %in% TRAIN_MONTHS) %>%
          inner_join(obs_psi_sub, by = "date") %>%
          filter(!is.na(psi_leaf), !is.na(psi_leaf_md_obs))

        if (nrow(joined_psi) >= 3) {
          m <- compute_metrics(joined_psi$psi_leaf, joined_psi$psi_leaf_md_obs)
          psi_score <- 0.4 * m$r_squared + 0.3 * m$kge +
            0.15 * (1 - pmin(abs(1 - m$slope), 1)) +
            0.15 * (1 - pmin(abs(m$bias), 1))
          psi_metrics <- m
        }
      }

      # --- Combined score ---
      # Weight: Gc gets 0.7 (dense daily data), PSI gets 0.3 (sparse)
      if (!is.na(gc_score) && !is.na(psi_score)) {
        combined_score <- 0.7 * gc_score + 0.3 * psi_score
      } else if (!is.na(gc_score)) {
        combined_score <- gc_score
      } else if (!is.na(psi_score)) {
        combined_score <- psi_score
      } else {
        next
      }

      if (combined_score > best_score) {
        best_score <- combined_score
        best_combo <- list(
          tag = cmb$tag,
          kr_max = cmb$kr_max, ks_max = cmb$ks_max,
          kl_max = cmb$kl_max, isohydricity = cmb$isohydricity,
          volumetric_capacitance = cmb$volumetric_capacitance,
          sapwood_theta_sat = cmb$sapwood_theta_sat,
          combined_score = combined_score,
          gc_score = gc_score, psi_score = psi_score,
          gc_r2 = if(!is.null(gc_metrics)) gc_metrics$r_squared else NA_real_,
          gc_kge = if(!is.null(gc_metrics)) gc_metrics$kge else NA_real_,
          gc_rmse = if(!is.null(gc_metrics)) gc_metrics$rmse else NA_real_,
          gc_n = if(!is.null(gc_metrics)) gc_metrics$n else NA_integer_,
          psi_r2 = if(!is.null(psi_metrics)) psi_metrics$r_squared else NA_real_,
          psi_kge = if(!is.null(psi_metrics)) psi_metrics$kge else NA_real_,
          psi_rmse = if(!is.null(psi_metrics)) psi_metrics$rmse else NA_real_,
          psi_n = if(!is.null(psi_metrics)) psi_metrics$n else NA_integer_
        )
      }
    }

    if (!is.null(best_combo)) {
      cat(sprintf("    BEST: %-15s  kr=%-5.1f ks=%-5.1f kl=%-5.1f iso=%-5.1f vc=%-4.0f theta=%-5.2f\n",
                  best_combo$tag, best_combo$kr_max, best_combo$ks_max, best_combo$kl_max,
                  best_combo$isohydricity, best_combo$volumetric_capacitance,
                  best_combo$sapwood_theta_sat))
      cat(sprintf("    Train Gc:  R²=%.3f  KGE=%.3f  n=%d\n",
                  best_combo$gc_r2, best_combo$gc_kge, best_combo$gc_n))
      cat(sprintf("    Train PSI: R²=%.3f  KGE=%.3f  n=%d\n",
                  best_combo$psi_r2, best_combo$psi_kge, best_combo$psi_n))

      all_results[[length(all_results) + 1]] <- data.frame(
        species       = sp,
        condition     = cond,
        train_months  = paste(TRAIN_MONTHS, collapse = "+"),
        best_run_tag  = best_combo$tag,
        kr_max        = best_combo$kr_max,
        ks_max        = best_combo$ks_max,
        kl_max        = best_combo$kl_max,
        isohydricity  = best_combo$isohydricity,
        volumetric_capacitance = best_combo$volumetric_capacitance,
        sapwood_theta_sat      = best_combo$sapwood_theta_sat,
        combined_score = best_combo$combined_score,
        train_gc_r2    = best_combo$gc_r2,
        train_gc_kge   = best_combo$gc_kge,
        train_gc_rmse  = best_combo$gc_rmse,
        train_gc_n     = best_combo$gc_n,
        train_psi_r2   = best_combo$psi_r2,
        train_psi_kge  = best_combo$psi_kge,
        train_psi_rmse = best_combo$psi_rmse,
        train_psi_n    = best_combo$psi_n
      )
    } else {
      cat(sprintf("    WARNING: No valid runs found for %s %s\n", sp, cond))
    }
  }
}

best_params_df <- bind_rows(all_results)

# ==========================================================================
# 6. BASELINE COMPARISON
# ==========================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  COMPARISON: BASELINE vs BEST (JULY+SEPT TRAINING)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Show what changed from baseline
if (nrow(best_params_df) > 0) {
  for (i in seq_len(nrow(best_params_df))) {
    r <- best_params_df[i, ]
    bases <- parameter_bases[[r$species]][[r$condition]]
    cat(sprintf("  %-8s %-10s:\n", r$species, r$condition))
    cat(sprintf("    kr_max:    baseline=%-5.1f → best=%-5.1f\n", bases$kr, r$kr_max))
    cat(sprintf("    ks_max:    baseline=%-5.1f → best=%-5.1f\n", bases$ks, r$ks_max))
    cat(sprintf("    kl_max:    baseline=%-5.1f → best=%-5.1f\n", bases$kl, r$kl_max))
    cat(sprintf("    isohydricity: baseline=%-5.1f → best=%-5.1f\n", bases$iso, r$isohydricity))
    cat(sprintf("    vol_capac: baseline=%-4.0f → best=%-4.0f\n", bases$vc, r$volumetric_capacitance))
    cat(sprintf("    theta_sat: baseline=%-5.2f → best=%-5.2f\n", bases$theta, r$sapwood_theta_sat))
    cat(sprintf("    Combined score: %.3f  (Train Gc R²=%.3f, PSI R²=%.3f)\n\n",
                r$combined_score, r$train_gc_r2, r$train_psi_r2))
  }
}

# ==========================================================================
# 7. EXPORT
# ==========================================================================

# Main output: best parameters
write.csv(best_params_df, file.path(OUTPUT_DIR, "best_params_trained_on_july_sept.csv"),
          row.names = FALSE)

# Also generate a .ins-ready parameter table
# This matches the format needed for sed substitution in the shell script
ins_params <- best_params_df %>%
  select(species, condition, kr_max, ks_max, kl_max, isohydricity,
         volumetric_capacitance, sapwood_theta_sat)

# Also save the parameters needed for the shell script's generate_run function
# (cav_slope, psi50_xylem, delta_psi_max, max_stem_shrinkage are not varied
#  in this analysis — they stay at species-specific values)
write.csv(ins_params, file.path(OUTPUT_DIR, "best_params_for_ins_files.csv"),
          row.names = FALSE)

# ==========================================================================
# 8. .INS FILE UPDATE INSTRUCTIONS
# ==========================================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  TO APPLY THESE PARAMETERS:\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("Option A: Manual update\n")
cat("  Update each species .ins file in LPJ-GUESS/ins_lrz/ins_control/ and ins_drought/\n")
cat("  with the values above.\n\n")

cat("Option B: Use sed commands (from LPJ_GUESS_HYD_TWD directory):\n\n")

# Generate sed commands
if (nrow(best_params_df) > 0) {
  for (i in seq_len(nrow(best_params_df))) {
    r <- best_params_df[i, ]
    species <- r$species
    cond <- r$condition
    ins_dir <- ifelse(cond == "control",
                      "LPJ-GUESS/ins_lrz/ins_control",
                      "LPJ-GUESS/ins_lrz/ins_drought")
    ins_file <- file.path(ins_dir, paste0(species, ".ins"))

    cat(sprintf("# %-8s %-10s\n", species, cond))
    cat(sprintf("sed -i \\\n"))
    cat(sprintf("  -e 's/^[[:space:]]*kr_max .*/    kr_max %.1f/' \\\n", r$kr_max))
    cat(sprintf("  -e 's/^[[:space:]]*ks_max .*/    ks_max %.1f/' \\\n", r$ks_max))
    cat(sprintf("  -e 's/^[[:space:]]*kl_max .*/    kl_max %.1f/' \\\n", r$kl_max))
    cat(sprintf("  -e 's/^[[:space:]]*isohydricity .*/    isohydricity %.1f/' \\\n", r$isohydricity))
    cat(sprintf("  -e 's/^[[:space:]]*volumetric_capacitance .*/    volumetric_capacitance    %.0f/' \\\n", r$volumetric_capacitance))
    cat(sprintf("  -e 's/^[[:space:]]*sapwood_theta_sat .*/    sapwood_theta_sat         %.2f/' \\\n", r$sapwood_theta_sat))
    cat(sprintf("  %s\n\n", ins_file))
  }
}

cat("\nThen run LPJ-GUESS with the updated .ins files.\n")
cat("After runs complete, validate on August with:\n")
cat("  Rscript R_scripts/run_pipeline_validation_august.R\n")

cat("\n*** Parameter selection complete. Results saved to:", OUTPUT_DIR, "***\n")
