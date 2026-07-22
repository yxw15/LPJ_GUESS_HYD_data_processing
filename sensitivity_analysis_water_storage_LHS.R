# ==========================================================================
# SENSITIVITY ANALYSIS (LHS): Water Storage + Hydraulic Parameters
# Evaluates 6 parameters via Latin Hypercube Sampling for 4 species × 2 conditions
#   kr_max, ks_max, kl_max, isohydricity,
#   volumetric_capacitance, sapwood_theta_sat
#
# Primary evaluation: TWD (Tree Water Deficit) against dendrometer observations
#   - Uses STANDARDIZED values (z-scores) for metric computation,
#     matching the approach in plot_compare_twd_lpjtwd_hoelstein.R
# Secondary outputs: stem RWC, stem diameter
#
# Usage (after all sensitivity runs complete and outputs are copied):
#   Rscript sensitivity_analysis_water_storage_LHS.R
#
# Output:
#   Figures/lpj_guess_stem_storage/sensitivity_water_storage_LHS/
#     all_metrics_TWD.csv                  — all metrics, all runs (TWD)
#     best_parameters_TWD.csv              — best parameter sets per species for TWD
#     parameter_correlation_TWD.png        — correlation of each parameter with KGE
#     parameter_scatter_TWD.png            — scatter of KGE vs each parameter
#     parameter_importance_TWD.png         — variable importance (random forest)
#     scatter_best_worst_TWD.png           — best vs worst scatter
#     parameter_ranking_TWD.png            — KGE ranking bar chart
# ==========================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(purrr)
library(stringr)

# ==========================================================================
# 0. PATHS & SETUP
# ==========================================================================
BASE_DIR       <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD"
RESULTS_BASE   <- file.path(BASE_DIR, "results_lpj/results_sensitivity_LHS_water_storage")
OUTPUT_DIR     <- file.path(BASE_DIR, "Figures/lpj_guess_stem_storage/sensitivity_water_storage_LHS")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# All paths are absolute — working directory doesn't matter
tryCatch(setwd(BASE_DIR), error = function(e) warning("Cannot setwd(): ", e$message))

# ==========================================================================
# 1. DEFINE PARAMETERS
# ==========================================================================
# Parameters varied in the LHS sensitivity analysis
PARAM_NAMES <- c("kr_max", "ks_max", "kl_max", "isohydricity",
                 "volumetric_capacitance", "sapwood_theta_sat")

# Species name mappings
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

species_colors <- c("Oak"="#E69F00", "Beech"="#0072B2", "Spruce"="#009E73", "Pine"="#F0E442")

# ==========================================================================
# 2. METRIC COMPUTATION FUNCTIONS
# ==========================================================================

compute_kge <- function(sim, obs) {
  # Kling-Gupta Efficiency
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
  nrmse_v <- (rmse_v / sd(o)) * 100
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

read_lpj_run <- function(run_dir, species_pft) {
  twd_file      <- file.path(run_dir, "twd.out")
  srwc_file     <- file.path(run_dir, "stem_rwc.out")
  sdiam_file    <- file.path(run_dir, "stem_diameter.out")

  if (!file.exists(twd_file)) {
    warning("Missing twd.out in: ", run_dir)
    return(NULL)
  }

  twd_data <- tryCatch(read.table(twd_file, header = TRUE, check.names = FALSE),
                       error = function(e) NULL)

  if (is.null(twd_data) || nrow(twd_data) < 10) return(NULL)

  colnames(twd_data) <- trimws(colnames(twd_data))

  if (!species_pft %in% colnames(twd_data)) {
    available <- paste(colnames(twd_data), collapse = ", ")
    warning("Column '", species_pft, "' not found in twd.out. Available: ", available)
    return(NULL)
  }

  result <- twd_data %>%
    mutate(
      date  = as.Date(Day, origin = paste0(Year, "-01-01")),
      twd   = .data[[species_pft]] * 1e6    # m → µm
    ) %>%
    select(date, twd)

  # Read stem RWC if available
  if (file.exists(srwc_file)) {
    srwc_data <- tryCatch(read.table(srwc_file, header = TRUE, check.names = FALSE),
                          error = function(e) NULL)
    if (!is.null(srwc_data)) {
      colnames(srwc_data) <- trimws(colnames(srwc_data))
      if (species_pft %in% colnames(srwc_data)) {
        srwc_sub <- srwc_data %>%
          mutate(date = as.Date(Day, origin = paste0(Year, "-01-01")),
                 stem_rwc = .data[[species_pft]]) %>%
          select(date, stem_rwc)
        result <- result %>% left_join(srwc_sub, by = "date")
      }
    }
  }

  # Read stem diameter if available
  if (file.exists(sdiam_file)) {
    sdiam_data <- tryCatch(read.table(sdiam_file, header = TRUE, check.names = FALSE),
                           error = function(e) NULL)
    if (!is.null(sdiam_data)) {
      colnames(sdiam_data) <- trimws(colnames(sdiam_data))
      if (species_pft %in% colnames(sdiam_data)) {
        sdiam_sub <- sdiam_data %>%
          mutate(date = as.Date(Day, origin = paste0(Year, "-01-01")),
                 stem_diameter = .data[[species_pft]]) %>%
          select(date, stem_diameter)
        result <- result %>% left_join(sdiam_sub, by = "date")
      }
    }
  }

  # result %>% filter(date >= as.Date("2023-01-01"))
}

# ==========================================================================
# 4. READ OBSERVED DATA
# ==========================================================================

obs_twd_daily <- tryCatch({
  dendro_files <- list.files(
    path = file.path(BASE_DIR, "SCCII/point_dendro"),
    pattern = "^Point_dendrometers_.*_archive\\.txt$",
    full.names = TRUE
  )

  dendro_obs <- map_dfr(dendro_files, ~ read.delim(.x), .id = "source_file")

  tree_info <- read.csv(file.path(BASE_DIR, "SCCII/tree_info.csv")) %>%
    mutate(treatment = ifelse(treatment == "treatment", "drought", treatment))

  dendro_obs %>%
    inner_join(tree_info, by = "tree_id") %>%
    mutate(
      date      = as.Date(str_extract(timestamp_UTC, "\\d{4}-\\d{2}-\\d{2}")),
      species   = factor(tolower(species), levels = c("oak", "beech", "spruce", "pine")),
      treatment = tolower(treatment)
    ) %>%
    filter(date >= as.Date("2023-01-01")) %>%
    group_by(date, species, treatment) %>%
    summarise(twd_obs = mean(twd_micron_treenetproc, na.rm = TRUE), .groups = "drop")
}, error = function(e) NULL)

# ==========================================================================
# 5. READ LHS PARAMETER TABLES
# ==========================================================================

read_lhs_csv <- function(csv_path) {
  if (!file.exists(csv_path)) return(NULL)
  df <- tryCatch(read.csv(csv_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df)) return(NULL)

  df <- df %>% mutate(
    sample_id = as.character(sample_id),
    species   = as.character(species),
    condition = as.character(condition)
  )
  return(df)
}

# Collect all LHS CSV files and add baseline runs
build_parameter_table <- function() {
  csv_files <- list.files(RESULTS_BASE, pattern = "^LHS_samples_.*\\.csv$",
                          full.names = TRUE, recursive = FALSE)

  lhs_all <- bind_rows(lapply(csv_files, read_lhs_csv))

  if (nrow(lhs_all) == 0) {
    # Build table from directory listing as fallback
    cat("  No LHS CSV files found; building parameter table from directory listing.\n")
    run_dirs <- list.dirs(RESULTS_BASE, recursive = FALSE, full.names = FALSE)
    # This won't give us parameter values, so return empty
    return(data.frame())
  }

  # Add baseline entries for each species × condition
  baselines <- lhs_all %>%
    group_by(species, condition) %>%
    slice_head(n = 1) %>%
    mutate(
      sample_id = "baseline",
      kr_max = NA_real_, ks_max = NA_real_, kl_max = NA_real_,
      isohydricity = NA_real_, volumetric_capacitance = NA_real_,
      sapwood_theta_sat = NA_real_
    ) %>%
    ungroup() %>%
    distinct(species, condition, .keep_all = TRUE)

  # Baselines are already in the CSV? Check if any row has sample_id == "baseline"
  if (!any(lhs_all$sample_id == "baseline")) {
    lhs_all <- bind_rows(lhs_all, baselines)
  }

  lhs_all
}

lhs_params <- build_parameter_table()

# ==========================================================================
# 6. MAIN EVALUATION LOOP
# ==========================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  LHS WATER STORAGE + HYDRAULIC SENSITIVITY ANALYSIS\n")
cat("  Parameters: ", paste(PARAM_NAMES, collapse = ", "), "\n")
cat("  (Using STANDARDIZED TWD values — ref: plot_compare_twd_lpjtwd_hoelstein.R)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

all_metrics_twd <- list()

for (sp in names(species_pft_map)) {

  cat("\n--- Processing species:", sp, "---\n")

  pft    <- species_pft_map[sp]
  sp_low <- species_lower_map[sp]

  for (cond in c("control", "drought")) {

    cat("  Condition:", cond, "\n")

    # Get parameter combinations for this species × condition
    if (nrow(lhs_params) > 0) {
      sp_cond_params <- lhs_params %>%
        filter(species == sp, condition == cond)
    } else {
      sp_cond_params <- data.frame()
    }

    # If no LHS CSV, fall back to scanning directories
    if (nrow(sp_cond_params) == 0) {
      cond_dir <- file.path(RESULTS_BASE, cond, sp)
      if (!dir.exists(cond_dir)) {
        cat("    [SKIP] Directory not found:", cond_dir, "\n")
        next
      }
      run_dirs <- list.dirs(cond_dir, recursive = FALSE, full.names = FALSE)
      cat("    Found", length(run_dirs), "run directories (no LHS CSV — parameter values unknown)\n")

      for (rn in run_dirs) {
        run_dir <- file.path(cond_dir, rn)
        lpj_data <- read_lpj_run(run_dir, pft)
        if (is.null(lpj_data) || nrow(lpj_data) == 0) {
          cat("      [SKIP] No valid data in:", rn, "\n")
          next
        }

        if (!is.null(obs_twd_daily)) {
          obs_sub <- obs_twd_daily %>%
            filter(species == sp_low, treatment == cond)

          joined <- lpj_data %>%
            inner_join(obs_sub, by = "date") %>%
            filter(!is.na(twd), !is.na(twd_obs))

          if (nrow(joined) > 10) {
            twd_obs_std <- (joined$twd_obs - mean(joined$twd_obs, na.rm = TRUE)) /
                           sd(joined$twd_obs, na.rm = TRUE)
            twd_sim_std <- (joined$twd - mean(joined$twd, na.rm = TRUE)) /
                           sd(joined$twd, na.rm = TRUE)

            m <- compute_metrics(twd_sim_std, twd_obs_std)
            m$species    <- sp
            m$condition  <- cond
            m$run_tag    <- rn
            for (p in PARAM_NAMES) m[[p]] <- NA_real_
            m$variable   <- "TWD"

            m$twd_mean_sim <- mean(joined$twd, na.rm = TRUE)
            m$twd_sd_sim   <- sd(joined$twd, na.rm = TRUE)
            m$twd_mean_obs <- mean(joined$twd_obs, na.rm = TRUE)
            m$twd_sd_obs   <- sd(joined$twd_obs, na.rm = TRUE)

            if ("stem_rwc" %in% colnames(joined)) {
              m$stem_rwc_mean <- mean(joined$stem_rwc, na.rm = TRUE)
              m$stem_rwc_sd   <- sd(joined$stem_rwc, na.rm = TRUE)
            }

            all_metrics_twd[[length(all_metrics_twd) + 1]] <- m
            cat(sprintf("      [OK] %-20s | R²=%.3f KGE=%.3f\n", rn, m$r_squared, m$kge))
          } else {
            cat("      [SKIP] Insufficient data in:", rn, "\n")
          }
        }
      }
      next
    }

    # Process with known parameter values from LHS CSV
    cat("    Evaluating", nrow(sp_cond_params), "parameter combinations\n")

    for (i in seq_len(nrow(sp_cond_params))) {
      row <- sp_cond_params[i, ]
      tag <- row$sample_id
      # Build directory name matching actual naming: {sp}_{cond}_LHS{tag} or {sp}_{cond}_baseline
      if (tag == "baseline") {
        dir_name <- sprintf("%s_%s_baseline", sp, cond)
      } else {
        dir_name <- sprintf("%s_%s_LHS%04d", sp, cond, as.integer(tag))
      }
      run_dir  <- file.path(RESULTS_BASE, cond, sp, dir_name)

      if (!dir.exists(run_dir)) {
        cat("      [SKIP] Directory not found:", dir_name, "\n")
        next
      }

      lpj_data <- read_lpj_run(run_dir, pft)
      if (is.null(lpj_data) || nrow(lpj_data) == 0) {
        cat("      [SKIP] No valid data in:", dir_name, "\n")
        next
      }

      if (!is.null(obs_twd_daily)) {
        obs_sub <- obs_twd_daily %>%
          filter(species == sp_low, treatment == cond)

        joined <- lpj_data %>%
          inner_join(obs_sub, by = "date") %>%
          filter(!is.na(twd), !is.na(twd_obs))

        if (nrow(joined) > 10) {
          twd_obs_std <- (joined$twd_obs - mean(joined$twd_obs, na.rm = TRUE)) /
                         sd(joined$twd_obs, na.rm = TRUE)
          twd_sim_std <- (joined$twd - mean(joined$twd, na.rm = TRUE)) /
                         sd(joined$twd, na.rm = TRUE)

          m <- compute_metrics(twd_sim_std, twd_obs_std)
          m$species    <- sp
          m$condition  <- cond
          m$run_tag    <- tag
          m$variable   <- "TWD"

          # Attach parameter values from LHS CSV
          for (p in PARAM_NAMES) {
            m[[p]] <- if (p %in% colnames(row)) row[[p]] else NA_real_
          }

          m$twd_mean_sim <- mean(joined$twd, na.rm = TRUE)
          m$twd_sd_sim   <- sd(joined$twd, na.rm = TRUE)
          m$twd_mean_obs <- mean(joined$twd_obs, na.rm = TRUE)
          m$twd_sd_obs   <- sd(joined$twd_obs, na.rm = TRUE)

          if ("stem_rwc" %in% colnames(joined)) {
            m$stem_rwc_mean <- mean(joined$stem_rwc, na.rm = TRUE)
            m$stem_rwc_sd   <- sd(joined$stem_rwc, na.rm = TRUE)
          }

          all_metrics_twd[[length(all_metrics_twd) + 1]] <- m

          cat(sprintf("      [OK] %-20s | R²=%.3f KGE=%.3f RMSE=%.1f\n",
                      tag, m$r_squared, m$kge, m$rmse))
        } else {
          cat("      [SKIP] Insufficient TWD data for:", tag, "\n")
        }
      }
    }
  }
}

# ==========================================================================
# 7. COMPILE RESULTS
# ==========================================================================
metrics_twd <- bind_rows(all_metrics_twd)

if (nrow(metrics_twd) == 0) {
  cat("\n*** No sensitivity results found. Have the runs completed and outputs been copied? ***\n")
  cat("Expected directory structure:\n")
  cat("  ", RESULTS_BASE, "/{condition}/{species}/{Species}_{condition}_{tag}/\n")
  quit(save = "no", status = 1)
}

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  RESULTS SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("  Total runs evaluated:", nrow(metrics_twd), "\n")
cat("  Species:", paste(unique(metrics_twd$species), collapse = ", "), "\n")
cat("  Conditions:", paste(unique(metrics_twd$condition), collapse = ", "), "\n\n")

# ==========================================================================
# 8. BEST PARAMETER SETS PER SPECIES × CONDITION
# ==========================================================================

select_best_params <- function(metrics_df) {
  if (nrow(metrics_df) == 0) return(data.frame())

  metrics_df %>%
    group_by(species, condition) %>%
    mutate(
      composite = 0.4 * r_squared +
                  0.4 * kge +
                  0.2 * (1 - pmin(abs(1 - slope), 1))
    ) %>%
    slice_max(composite, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(species, condition, run_tag,
           any_of(PARAM_NAMES),
           r_squared, kge, rmse, nrmse_pct, bias, slope,
           twd_mean_sim, twd_sd_sim, composite)
}

best_twd <- select_best_params(metrics_twd)

cat("\n--- Best Parameter Sets (TWD, standardized) ---\n")
if (nrow(best_twd) > 0) print(as.data.frame(best_twd))

# Top-5 per species×condition
top5_twd <- metrics_twd %>%
  group_by(species, condition) %>%
  slice_max(kge, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(species, condition, desc(kge))

cat("\n--- Top 5 Parameter Sets by KGE (TWD, standardized) ---\n")
if (nrow(top5_twd) > 0) {
  print(as.data.frame(
    top5_twd %>% select(species, condition, run_tag, any_of(PARAM_NAMES),
                        r_squared, kge, rmse)
  ))
}

# ==========================================================================
# 9. PARAMETER IMPORTANCE VIA CORRELATION
# ==========================================================================

# For each species × condition, compute Spearman correlation of each
# parameter with KGE (excluding baseline runs)
correlation_analysis <- metrics_twd %>%
  filter(run_tag != "baseline") %>%
  select(species, condition, kge, all_of(PARAM_NAMES)) %>%
  pivot_longer(all_of(PARAM_NAMES), names_to = "parameter", values_to = "value") %>%
  group_by(species, condition, parameter) %>%
  summarise(
    spearman_r  = cor(value, kge, method = "spearman", use = "complete.obs"),
    pearson_r   = cor(value, kge, method = "pearson", use = "complete.obs"),
    n_valid     = sum(complete.cases(value, kge)),
    .groups = "drop"
  )

# ==========================================================================
# 10. EXPORT TABLES
# ==========================================================================

write.csv(metrics_twd, file.path(OUTPUT_DIR, "all_metrics_TWD.csv"), row.names = FALSE)
if (nrow(best_twd) > 0) write.csv(best_twd, file.path(OUTPUT_DIR, "best_parameters_TWD.csv"), row.names = FALSE)
if (nrow(top5_twd) > 0) write.csv(top5_twd, file.path(OUTPUT_DIR, "top5_parameters_TWD.csv"), row.names = FALSE)
write.csv(correlation_analysis, file.path(OUTPUT_DIR, "parameter_correlation_TWD.csv"), row.names = FALSE)

# ==========================================================================
# 11. DIAGNOSTIC PLOTS
# ==========================================================================

base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    plot.title        = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5, size = 10, color = "grey30"),
    axis.title        = element_text(size = 12),
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 8),
    strip.text        = element_text(size = 10, face = "bold"),
    panel.grid.major  = element_line(color = "grey92"),
    panel.grid.minor  = element_blank()
  )

# --- 11a. Correlation heatmap: parameter vs KGE ---
create_correlation_plot <- function(cor_df) {
  if (nrow(cor_df) == 0) return(NULL)

  ggplot(cor_df, aes(x = parameter, y = condition, fill = spearman_r)) +
    geom_tile(color = "white", linewidth = 0.5) +
    facet_wrap(~ species, ncol = 2) +
    scale_fill_gradient2(
      name = "Spearman ρ",
      low = "#2166AC", mid = "white", high = "#B2182B",
      limits = c(-1, 1)
    ) +
    geom_text(aes(label = sprintf("%.2f", spearman_r)), size = 3) +
    labs(
      title = "Parameter–KGE Correlation (TWD, standardized)",
      subtitle = paste("Spearman rank correlation. LHS, n =",
                       max(cor_df$n_valid, na.rm = TRUE)),
      x = NULL, y = NULL
    ) +
    base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
}

# --- 11b. Scatter: KGE vs each parameter ---
create_parameter_scatter <- function(metrics_df) {
  if (nrow(metrics_df) == 0) return(NULL)

  plot_data <- metrics_df %>%
    filter(run_tag != "baseline") %>%
    select(species, condition, kge, all_of(PARAM_NAMES)) %>%
    pivot_longer(all_of(PARAM_NAMES), names_to = "parameter", values_to = "value") %>%
    filter(!is.na(value))

  ggplot(plot_data, aes(x = value, y = kge, color = species)) +
    geom_point(alpha = 0.4, size = 1.5) +
    geom_smooth(aes(group = species), method = "loess", se = FALSE, linewidth = 0.8) +
    facet_grid(condition ~ parameter, scales = "free_x") +
    scale_color_manual(values = species_colors) +
    labs(
      title = "KGE vs Parameter Values (TWD, standardized)",
      subtitle = "LOESS smooth. Each point = 1 LHS run.",
      x = "Parameter value",
      y = "KGE (TWD, standardized)",
      color = "Species"
    ) +
    base_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7))
}

# --- 11c. Variable importance via random forest (if available) ---
create_importance_plot <- function(metrics_df) {
  if (nrow(metrics_df) == 0) return(NULL)
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    warning("Package 'randomForest' not installed. Skipping variable importance plot.")
    return(NULL)
  }

  imp_list <- list()

  for (sp in unique(metrics_df$species)) {
    for (cd in unique(metrics_df$condition)) {
      sub <- metrics_df %>%
        filter(species == sp, condition == cd, run_tag != "baseline") %>%
        select(kge, all_of(PARAM_NAMES)) %>%
        filter(complete.cases(.))

      if (nrow(sub) < 20) next

      rf <- tryCatch(
        randomForest::randomForest(
          kge ~ ., data = sub,
          ntree = 500, importance = TRUE
        ),
        error = function(e) NULL
      )

      if (is.null(rf)) next

      imp <- randomForest::importance(rf, type = 1)  # %IncMSE
      imp_df <- data.frame(
        species   = sp,
        condition = cd,
        parameter = rownames(imp),
        inc_mse_pct = as.numeric(imp[, 1])
      )
      imp_list[[length(imp_list) + 1]] <- imp_df
    }
  }

  if (length(imp_list) == 0) return(NULL)

  imp_all <- bind_rows(imp_list)

  # Normalize within each species×condition to 0–1 for comparability
  imp_all <- imp_all %>%
    group_by(species, condition) %>%
    mutate(inc_mse_norm = inc_mse_pct / max(abs(inc_mse_pct), na.rm = TRUE)) %>%
    ungroup()

  ggplot(imp_all, aes(x = reorder(parameter, inc_mse_norm), y = inc_mse_norm, fill = species)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~ condition, ncol = 1) +
    coord_flip() +
    scale_fill_manual(values = species_colors) +
    labs(
      title = "Variable Importance: %IncMSE (Random Forest)",
      subtitle = "Normalized per species×condition. Higher = more important for KGE.",
      x = NULL,
      y = "Normalized %IncMSE",
      fill = "Species"
    ) +
    base_theme +
    theme(axis.text.y = element_text(size = 10))
}

# --- 11d. Scatter: Best vs Worst TWD (standardized) ---
create_scatter_best_worst <- function(metrics_df) {
  if (nrow(metrics_df) == 0) return(NULL)

  best_worst <- metrics_df %>%
    group_by(species, condition) %>%
    summarise(
      best_tag  = run_tag[which.max(kge)],
      worst_tag = run_tag[which.min(kge)],
      best_kge  = max(kge, na.rm = TRUE),
      worst_kge = min(kge, na.rm = TRUE),
      .groups   = "drop"
    )

  scatter_data <- list()
  for (i in seq_len(nrow(best_worst))) {
    sp   <- best_worst$species[i]
    cond <- best_worst$condition[i]
    pft  <- species_pft_map[sp]
    sp_low <- species_lower_map[sp]

    for (rank in c("best", "worst")) {
      tag <- if (rank == "best") best_worst$best_tag[i] else best_worst$worst_tag[i]
      # Build directory name matching actual naming: {sp}_{cond}_LHS{tag} or {sp}_{cond}_baseline
      if (tag == "baseline") {
        dir_name <- sprintf("%s_%s_baseline", sp, cond)
      } else {
        dir_name <- sprintf("%s_%s_LHS%04d", sp, cond, as.integer(tag))
      }
      run_dir  <- file.path(RESULTS_BASE, cond, sp, dir_name)

      if (!dir.exists(run_dir)) next

      lpj_data <- read_lpj_run(run_dir, pft)
      if (is.null(lpj_data)) next

      obs_sub <- obs_twd_daily %>%
        filter(species == sp_low, treatment == cond)

      joined <- lpj_data %>%
        inner_join(obs_sub, by = "date") %>%
        filter(!is.na(twd), !is.na(twd_obs))

      if (nrow(joined) > 10) {
        joined$twd_obs_std <- (joined$twd_obs - mean(joined$twd_obs, na.rm = TRUE)) /
                               sd(joined$twd_obs, na.rm = TRUE)
        joined$twd_sim_std <- (joined$twd - mean(joined$twd, na.rm = TRUE)) /
                               sd(joined$twd, na.rm = TRUE)
        joined$species   <- sp
        joined$condition <- cond
        joined$rank      <- rank
        joined$tag       <- tag
        scatter_data[[length(scatter_data) + 1]] <- joined
      }
    }
  }

  if (length(scatter_data) == 0) return(NULL)

  scatter_all <- bind_rows(scatter_data)

  ggplot(scatter_all, aes(x = twd_obs_std, y = twd_sim_std, color = species)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "grey40", linewidth = 0.5) +
    facet_grid(condition + rank ~ species, scales = "fixed") +
    coord_fixed(ratio = 1) +
    scale_color_manual(values = species_colors) +
    labs(
      title = "TWD: Best vs Worst Parameter Sets (standardized)",
      subtitle = "Standardized TWD (z-score). 1:1 line dashed.",
      x = "Observed TWD (z-score)",
      y = "Simulated TWD (z-score)",
      color = "Species"
    ) +
    base_theme +
    theme(legend.position = "none")
}

# --- 11e. KGE Ranking Bar Chart ---
create_kge_ranking <- function(metrics_df) {
  if (nrow(metrics_df) == 0) return(NULL)

  top_runs <- metrics_df %>%
    group_by(species, condition) %>%
    slice_max(kge, n = 15, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      facet_grp = paste(species, condition, sep = "__")
    )

  top_runs <- top_runs %>%
    group_by(facet_grp) %>%
    mutate(label_ordered = reorder(run_tag, kge)) %>%
    ungroup()

  ggplot(top_runs, aes(x = label_ordered, y = kge, fill = species)) +
    geom_col() +
    facet_grid(condition ~ species, scales = "free_y") +
    coord_flip() +
    scale_fill_manual(values = species_colors) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    labs(
      title = "Top 15 Parameter Sets by KGE (TWD, standardized)",
      x = NULL,
      y = "KGE"
    ) +
    base_theme +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8)
    )
}

# ==========================================================================
# 12. GENERATE ALL PLOTS
# ==========================================================================

cat("\n--- Generating plots ---\n\n")

p_cor <- create_correlation_plot(correlation_analysis)
if (!is.null(p_cor)) {
  ggsave(file.path(OUTPUT_DIR, "parameter_correlation_TWD.png"),
         p_cor, width = 14, height = 8, dpi = 300)
  cat("  Saved: parameter_correlation_TWD.png\n")
}

p_scatter <- create_parameter_scatter(metrics_twd)
if (!is.null(p_scatter)) {
  ggsave(file.path(OUTPUT_DIR, "parameter_scatter_TWD.png"),
         p_scatter, width = 18, height = 10, dpi = 300)
  cat("  Saved: parameter_scatter_TWD.png\n")
}

p_imp <- create_importance_plot(metrics_twd)
if (!is.null(p_imp)) {
  ggsave(file.path(OUTPUT_DIR, "parameter_importance_TWD.png"),
         p_imp, width = 12, height = 10, dpi = 300)
  cat("  Saved: parameter_importance_TWD.png\n")
}

p_bw <- create_scatter_best_worst(metrics_twd)
if (!is.null(p_bw)) {
  ggsave(file.path(OUTPUT_DIR, "scatter_best_worst_TWD.png"),
         p_bw, width = 14, height = 10, dpi = 300)
  cat("  Saved: scatter_best_worst_TWD.png\n")
}

p_rank <- create_kge_ranking(metrics_twd)
if (!is.null(p_rank)) {
  ggsave(file.path(OUTPUT_DIR, "parameter_ranking_TWD.png"),
         p_rank, width = 16, height = 14, dpi = 300)
  cat("  Saved: parameter_ranking_TWD.png\n")
}

# ==========================================================================
# 13. SUMMARY: PARAMETER IMPORTANCE
# ==========================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  PARAMETER IMPORTANCE SUMMARY (Spearman correlation with KGE)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

if (nrow(correlation_analysis) > 0) {
  # Show top-2 parameters per species×condition
  top_params <- correlation_analysis %>%
    group_by(species, condition) %>%
    slice_max(abs(spearman_r), n = 2, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(species, condition, desc(abs(spearman_r)))

  for (i in seq_len(nrow(top_params))) {
    r <- top_params[i, ]
    cat(sprintf("  %-8s %-10s: %-25s ρ=%.3f\n",
                r$species, r$condition, r$parameter, r$spearman_r))
  }
}

# Baseline comparison
baseline_metrics <- metrics_twd %>% filter(run_tag == "baseline")
if (nrow(baseline_metrics) > 0) {
  cat("\n  Baseline KGE values:\n")
  for (i in seq_len(nrow(baseline_metrics))) {
    r <- baseline_metrics[i, ]
    cat(sprintf("  %-8s %-10s: KGE=%.3f R²=%.3f\n",
                r$species, r$condition, r$kge, r$r_squared))
  }
}

# KGE improvement of best LHS over baseline
if (nrow(best_twd) > 0 && nrow(baseline_metrics) > 0) {
  cat("\n  Best LHS vs Baseline:\n")
  for (i in seq_len(nrow(best_twd))) {
    r <- best_twd[i, ]
    bl <- baseline_metrics %>%
      filter(species == r$species, condition == r$condition)
    if (nrow(bl) > 0) {
      cat(sprintf("  %-8s %-10s: best_KGE=%.3f baseline_KGE=%.3f (Δ=%+.3f)\n",
                  r$species, r$condition, r$kge, bl$kge[1], r$kge - bl$kge[1]))
    }
  }
}

# ==========================================================================
# 14. FINAL RECOMMENDATIONS
# ==========================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  FINAL PARAMETER RECOMMENDATIONS (LHS Sensitivity)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("To apply best parameters, update the .ins files with:\n\n")

if (nrow(best_twd) > 0) {
  for (i in seq_len(nrow(best_twd))) {
    r <- best_twd[i, ]
    # Build parameter string from whatever columns are available
    param_str <- ""
    for (p in PARAM_NAMES) {
      if (p %in% colnames(r) && !is.na(r[[p]])) {
        param_str <- paste0(param_str, sprintf("%s=%-6s ", p, format(r[[p]], digits = 3)))
      }
    }
    cat(sprintf("  %-8s %-10s: %s (KGE=%.3f, R²=%.3f)\n",
                r$species, r$condition, param_str, r$kge, r$r_squared))
  }
}

cat("\n*** LHS water storage sensitivity analysis complete. Results saved to:", OUTPUT_DIR, "***\n")
