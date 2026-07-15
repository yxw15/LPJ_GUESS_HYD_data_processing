# ==========================================================================
# SENSITIVITY ANALYSIS: Water Storage Parameters
# Evaluates volumetric_capacitance × sapwood_theta_sat for 4 species × 2 conditions
#
# Primary evaluation: TWD (Tree Water Deficit) against dendrometer observations
#   - Uses STANDARDIZED values (z-scores) for metric computation,
#     matching the approach in plot_compare_twd_lpjtwd_hoelstein.R
# Secondary outputs: stem RWC, stem diameter
#
# Usage (after all sensitivity runs complete and outputs are copied):
#   Rscript sensitivity_analysis_water_storage.R
#
# Output:
#   Figures/lpj_guess_stem_storage/sensitivity_water_storage/
#     all_metrics_TWD.csv                — all metrics, all runs (TWD)
#     best_parameters_TWD.csv            — best parameter sets per species for TWD
#     parameter_sensitivity_TWD_heatmap_R2.png  — 2D heatmap (R²)
#     parameter_sensitivity_TWD_heatmap_KGE.png — 2D heatmap (KGE)
#     parameter_response_curve_TWD_vc.png       — response curves vs vc
#     parameter_response_curve_TWD_theta.png    — response curves vs theta
#     parameter_main_effects_TWD_vc.png         — main effects (vc)
#     parameter_main_effects_TWD_theta.png      — main effects (theta)
#     scatter_best_worst_TWD.png               — best vs worst scatter
#     parameter_ranking_TWD.png                — KGE ranking bar chart
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
RESULTS_BASE   <- file.path(BASE_DIR, "results_lpj/results_sensitivity_water_storage")
OUTPUT_DIR     <- file.path(BASE_DIR, "Figures/lpj_guess_stem_storage/sensitivity_water_storage")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

setwd(BASE_DIR)

# ==========================================================================
# 1. DEFINE PARAMETER GRID
# ==========================================================================
# Full factorial: volumetric_capacitance × sapwood_theta_sat
VC_VALS    <- c(50, 100, 150, 200)
THETA_VALS <- c(0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70)

# Species-specific baseline values for the two water storage parameters VARIED in
# the sensitivity analysis. All other hydraulic parameters (kr_max, ks_max, kl_max,
# isohydricity, cav_slope, psi50_xylem, delta_psi_max, max_stem_shrinkage) are
# taken as-is from the template .ins files and are NOT listed or overridden here.
parameter_bases <- list(
  Oak    = list(vc = 200, theta = 0.55),
  Beech  = list(vc = 200, theta = 0.55),
  Spruce = list(vc = 300, theta = 0.55),
  Pine   = list(vc = 300, theta = 0.55)
)

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
  nrmse_v <- (rmse_v / sd(o)) * 100  # NRMSE relative to obs SD (for standardized data)
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

  # Extract TWD and convert to microns (LPJ outputs in m, ×1e6 → µm)
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

  result %>% filter(date >= as.Date("2023-01-01"))
}

# ==========================================================================
# 4. READ OBSERVED DATA
# ==========================================================================

# TWD observations from point dendrometers
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
# 5. BUILD PARAMETER COMBINATIONS
# ==========================================================================

build_combos <- function(sp) {
  bases <- parameter_bases[[sp]]
  combos <- list()

  # 1) Baseline (species-specific vc and theta = 0.55)
  combos[[1]] <- list(
    tag   = "baseline",
    vc    = bases$vc,
    theta = bases$theta
  )

  # 2) Full factorial grid
  combo_idx <- 2
  for (vc in VC_VALS) {
    for (theta in THETA_VALS) {
      tag <- sprintf("vc%d_th%.2f", vc, theta)
      combos[[combo_idx]] <- list(
        tag   = tag,
        vc    = vc,
        theta = theta
      )
      combo_idx <- combo_idx + 1
    }
  }

  combos
}

# ==========================================================================
# 6. MAIN EVALUATION LOOP
# ==========================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  WATER STORAGE SENSITIVITY ANALYSIS: TWD Parameter Evaluation\n")
cat("  (Using STANDARDIZED TWD values — ref: plot_compare_twd_lpjtwd_hoelstein.R)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

all_metrics_twd <- list()

for (sp in names(parameter_bases)) {

  cat("\n--- Processing species:", sp, "---\n")

  bases  <- parameter_bases[[sp]]
  pft    <- species_pft_map[sp]
  sp_low <- species_lower_map[sp]

  combos <- build_combos(sp)

  for (cond in c("control", "drought")) {

    cat("  Condition:", cond, "\n")
    cat("    Evaluating", length(combos), "parameter combinations\n")

    for (cmb in combos) {
      dir_name <- sprintf("%s_%s_%s", sp, cond, cmb$tag)
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

      # ---- TWD Evaluation with standardization ----
      # Standardization follows plot_compare_twd_lpjtwd_hoelstein.R:
      # each variable is z-scored independently so the evaluation
      # focuses on temporal pattern (correlation structure) rather
      # than absolute magnitude.
      if (!is.null(obs_twd_daily)) {
        obs_sub <- obs_twd_daily %>%
          filter(species == sp_low, treatment == cond)

        joined <- lpj_data %>%
          inner_join(obs_sub, by = "date") %>%
          filter(!is.na(twd), !is.na(twd_obs))

        if (nrow(joined) > 10) {
          # Z-score standardization (independently for sim and obs)
          twd_obs_std <- (joined$twd_obs - mean(joined$twd_obs, na.rm = TRUE)) /
                         sd(joined$twd_obs, na.rm = TRUE)
          twd_sim_std <- (joined$twd - mean(joined$twd, na.rm = TRUE)) /
                         sd(joined$twd, na.rm = TRUE)

          m <- compute_metrics(twd_sim_std, twd_obs_std)
          m$species                <- sp
          m$condition              <- cond
          m$run_tag                <- cmb$tag
          m$volumetric_capacitance <- cmb$vc
          m$sapwood_theta_sat      <- cmb$theta
          m$variable               <- "TWD"

          # Store raw (non-standardized) stats for supplementary analysis
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
                      cmb$tag, m$r_squared, m$kge, m$rmse))
        } else {
          cat("      [SKIP] Insufficient TWD data for:", cmb$tag, "\n")
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
           volumetric_capacitance, sapwood_theta_sat,
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
    top5_twd %>% select(species, condition, run_tag, volumetric_capacitance,
                        sapwood_theta_sat, r_squared, kge, rmse)
  ))
}

# ==========================================================================
# 9. MAIN EFFECTS ANALYSIS
# ==========================================================================

main_effects <- metrics_twd %>%
  filter(run_tag != "baseline") %>%
  group_by(species, condition, volumetric_capacitance, sapwood_theta_sat) %>%
  summarise(
    mean_r2  = mean(r_squared, na.rm = TRUE),
    mean_kge = mean(kge, na.rm = TRUE),
    sd_r2    = sd(r_squared, na.rm = TRUE),
    .groups  = "drop"
  )

# Marginal effect of vc (averaged over theta)
vc_effect <- metrics_twd %>%
  filter(run_tag != "baseline") %>%
  group_by(species, condition, volumetric_capacitance) %>%
  summarise(
    mean_r2  = mean(r_squared, na.rm = TRUE),
    mean_kge = mean(kge, na.rm = TRUE),
    sd_r2    = sd(r_squared, na.rm = TRUE),
    .groups  = "drop"
  )

# Marginal effect of theta (averaged over vc)
theta_effect <- metrics_twd %>%
  filter(run_tag != "baseline") %>%
  group_by(species, condition, sapwood_theta_sat) %>%
  summarise(
    mean_r2  = mean(r_squared, na.rm = TRUE),
    mean_kge = mean(kge, na.rm = TRUE),
    sd_r2    = sd(r_squared, na.rm = TRUE),
    .groups  = "drop"
  )

# ==========================================================================
# 10. EXPORT TABLES
# ==========================================================================

write.csv(metrics_twd, file.path(OUTPUT_DIR, "all_metrics_TWD.csv"), row.names = FALSE)
if (nrow(best_twd) > 0) write.csv(best_twd, file.path(OUTPUT_DIR, "best_parameters_TWD.csv"), row.names = FALSE)
if (nrow(top5_twd) > 0) write.csv(top5_twd, file.path(OUTPUT_DIR, "top5_parameters_TWD.csv"), row.names = FALSE)
write.csv(main_effects, file.path(OUTPUT_DIR, "main_effects_TWD.csv"), row.names = FALSE)

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

# --- 11a. 2D Heatmap: R² / KGE across vc × theta grid ---
create_2d_heatmap <- function(metrics_df, metric_var = "r_squared", metric_label = "R²") {
  if (nrow(metrics_df) == 0) return(NULL)

  plot_data <- metrics_df %>%
    filter(run_tag != "baseline") %>%
    mutate(
      vc_label    = factor(paste0("vc=", volumetric_capacitance),
                           levels = paste0("vc=", sort(unique(volumetric_capacitance)))),
      theta_label = factor(paste0("θ=", sprintf("%.2f", sapwood_theta_sat)),
                           levels = paste0("θ=", sprintf("%.2f", sort(unique(sapwood_theta_sat)))))
    )

  metric_val <- plot_data[[metric_var]]

  p <- ggplot(plot_data, aes(x = vc_label, y = theta_label)) +
    geom_tile(aes(fill = metric_val), color = "white", linewidth = 0.5) +
    facet_grid(condition ~ species) +
    scale_fill_viridis_c(
      name = metric_label,
      option = "D",
      limits = if (metric_var == "r_squared") c(0, 1) else c(min(metric_val, na.rm = TRUE), 1)
    ) +
    labs(
      title = paste("Water Storage Sensitivity:", metric_label, "for TWD (standardized)"),
      subtitle = "Full factorial: volumetric_capacitance × sapwood_theta_sat. Baseline: ★",
      x = "Volumetric Capacitance (kg m⁻³ MPa⁻¹)",
      y = "Sapwood Theta Sat (m³ m⁻³)"
    ) +
    base_theme

  # Mark baseline if it falls within the grid
  baseline_data <- metrics_df %>%
    filter(run_tag == "baseline") %>%
    mutate(
      vc_label    = factor(paste0("vc=", volumetric_capacitance),
                           levels = levels(plot_data$vc_label)),
      theta_label = factor(paste0("θ=", sprintf("%.2f", sapwood_theta_sat)),
                           levels = levels(plot_data$theta_label))
    ) %>%
    filter(vc_label %in% levels(plot_data$vc_label),
           theta_label %in% levels(plot_data$theta_label))

  if (nrow(baseline_data) > 0) {
    p <- p + geom_point(
      data = baseline_data,
      aes(x = vc_label, y = theta_label),
      shape = 8, size = 4, color = "red", stroke = 1.2,
      inherit.aes = FALSE
    )
  }

  p
}

# --- 11b. Response Curves: KGE vs vc, colored by theta ---
create_response_vc <- function(metrics_df) {
  if (nrow(metrics_df) == 0) return(NULL)

  plot_data <- metrics_df %>%
    filter(run_tag != "baseline") %>%
    mutate(theta_group = factor(paste0("θ=", sprintf("%.2f", sapwood_theta_sat))))

  p <- ggplot(plot_data, aes(x = volumetric_capacitance, y = kge,
                              color = theta_group, group = theta_group)) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    geom_point(size = 1.5) +
    facet_grid(condition ~ species) +
    scale_color_viridis_d(name = "Sapwood θ_sat", option = "D") +
    labs(
      title = "TWD KGE vs Volumetric Capacitance (standardized)",
      subtitle = "Lines colored by sapwood_theta_sat. Baseline: ★",
      x = "Volumetric Capacitance (kg m⁻³ MPa⁻¹)",
      y = "KGE (TWD, standardized)"
    ) +
    base_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9))

  baseline_data <- metrics_df %>% filter(run_tag == "baseline")
  if (nrow(baseline_data) > 0) {
    p <- p + geom_point(
      data = baseline_data %>% filter(volumetric_capacitance %in% VC_VALS),
      aes(x = volumetric_capacitance, y = kge),
      shape = 8, size = 4, color = "red", inherit.aes = FALSE
    )
  }
  p
}

# --- 11c. Response Curves: KGE vs theta, colored by vc ---
create_response_theta <- function(metrics_df) {
  if (nrow(metrics_df) == 0) return(NULL)

  plot_data <- metrics_df %>%
    filter(run_tag != "baseline") %>%
    mutate(vc_group = factor(paste0("vc=", volumetric_capacitance)))

  p <- ggplot(plot_data, aes(x = sapwood_theta_sat, y = kge,
                              color = vc_group, group = vc_group)) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    geom_point(size = 1.5) +
    facet_grid(condition ~ species) +
    scale_color_viridis_d(name = "Vol. Capacitance", option = "D") +
    labs(
      title = "TWD KGE vs Sapwood Theta Sat (standardized)",
      subtitle = "Lines colored by volumetric_capacitance. Baseline: ★",
      x = "Sapwood Theta Sat (m³ m⁻³)",
      y = "KGE (TWD, standardized)"
    ) +
    base_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9))

  baseline_data <- metrics_df %>% filter(run_tag == "baseline")
  if (nrow(baseline_data) > 0) {
    p <- p + geom_point(
      data = baseline_data %>% filter(sapwood_theta_sat %in% THETA_VALS),
      aes(x = sapwood_theta_sat, y = kge),
      shape = 8, size = 4, color = "red", inherit.aes = FALSE
    )
  }
  p
}

# --- 11d. Main Effects Plot ---
create_main_effects <- function(vc_eff, theta_eff) {
  if (nrow(vc_eff) == 0 && nrow(theta_eff) == 0) return(NULL)

  p_vc <- ggplot(vc_eff, aes(x = factor(volumetric_capacitance), y = mean_kge, fill = species)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(
      aes(ymin = mean_kge - sd_r2, ymax = mean_kge + sd_r2),
      position = position_dodge(width = 0.8), width = 0.2
    ) +
    facet_wrap(~ condition, ncol = 1) +
    scale_fill_manual(values = species_colors) +
    labs(
      title = "Main Effect: Volumetric Capacitance on TWD KGE",
      subtitle = "KGE averaged over all sapwood_theta_sat values (±1 SD)",
      x = "Volumetric Capacitance (kg m⁻³ MPa⁻¹)",
      y = "Mean KGE (TWD, standardized)",
      fill = "Species"
    ) +
    base_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10))

  p_theta <- ggplot(theta_eff, aes(x = factor(sapwood_theta_sat), y = mean_kge, fill = species)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(
      aes(ymin = mean_kge - sd_r2, ymax = mean_kge + sd_r2),
      position = position_dodge(width = 0.8), width = 0.2
    ) +
    facet_wrap(~ condition, ncol = 1) +
    scale_fill_manual(values = species_colors) +
    labs(
      title = "Main Effect: Sapwood Theta Sat on TWD KGE",
      subtitle = "KGE averaged over all volumetric_capacitance values (±1 SD)",
      x = "Sapwood Theta Sat (m³ m⁻³)",
      y = "Mean KGE (TWD, standardized)",
      fill = "Species"
    ) +
    base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

  list(vc = p_vc, theta = p_theta)
}

# --- 11e. Scatter: Best vs Worst TWD (standardized) ---
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
      dir_name <- sprintf("%s_%s_%s", sp, cond, tag)
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
        # Standardize for display
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

# --- 11f. KGE Ranking Bar Chart ---
create_kge_ranking <- function(metrics_df) {
  if (nrow(metrics_df) == 0) return(NULL)

  top_runs <- metrics_df %>%
    group_by(species, condition) %>%
    slice_max(kge, n = 15, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      label = paste0(run_tag, " (vc=", volumetric_capacitance,
                     ", θ=", sprintf("%.2f", sapwood_theta_sat), ")"),
      # Create ordering labels that respect species×condition groups
      facet_grp = paste(species, condition, sep = "__")
    )

  # Build ordered factor manually (avoids tidytext dependency)
  top_runs <- top_runs %>%
    group_by(facet_grp) %>%
    mutate(label_ordered = reorder(label, kge)) %>%
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

p_heat_r2   <- create_2d_heatmap(metrics_twd, "r_squared", "R²")
p_heat_kge  <- create_2d_heatmap(metrics_twd, "kge", "KGE")

if (!is.null(p_heat_r2)) {
  ggsave(file.path(OUTPUT_DIR, "parameter_sensitivity_TWD_heatmap_R2.png"),
         p_heat_r2, width = 16, height = 9, dpi = 300)
  cat("  Saved: parameter_sensitivity_TWD_heatmap_R2.png\n")
}

if (!is.null(p_heat_kge)) {
  ggsave(file.path(OUTPUT_DIR, "parameter_sensitivity_TWD_heatmap_KGE.png"),
         p_heat_kge, width = 16, height = 9, dpi = 300)
  cat("  Saved: parameter_sensitivity_TWD_heatmap_KGE.png\n")
}

p_resp_vc    <- create_response_vc(metrics_twd)
p_resp_theta <- create_response_theta(metrics_twd)

if (!is.null(p_resp_vc)) {
  ggsave(file.path(OUTPUT_DIR, "parameter_response_curve_TWD_vc.png"),
         p_resp_vc, width = 16, height = 9, dpi = 300)
  cat("  Saved: parameter_response_curve_TWD_vc.png\n")
}

if (!is.null(p_resp_theta)) {
  ggsave(file.path(OUTPUT_DIR, "parameter_response_curve_TWD_theta.png"),
         p_resp_theta, width = 16, height = 9, dpi = 300)
  cat("  Saved: parameter_response_curve_TWD_theta.png\n")
}

me_plots <- create_main_effects(vc_effect, theta_effect)
if (!is.null(me_plots)) {
  ggsave(file.path(OUTPUT_DIR, "parameter_main_effects_TWD_vc.png"),
         me_plots$vc, width = 12, height = 10, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, "parameter_main_effects_TWD_theta.png"),
         me_plots$theta, width = 14, height = 10, dpi = 300)
  cat("  Saved: parameter_main_effects_TWD_vc.png\n")
  cat("  Saved: parameter_main_effects_TWD_theta.png\n")
}

p_scatter <- create_scatter_best_worst(metrics_twd)
if (!is.null(p_scatter)) {
  ggsave(file.path(OUTPUT_DIR, "scatter_best_worst_TWD.png"),
         p_scatter, width = 14, height = 10, dpi = 300)
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
cat("  PARAMETER IMPORTANCE SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# For each species×condition: how much does KGE vary when changing
# vc (at median theta) vs changing theta (at median vc)?
importance_summary <- metrics_twd %>%
  filter(run_tag != "baseline") %>%
  group_by(species, condition) %>%
  summarise(
    kge_range_vc = {
      med_theta <- THETA_VALS[round(length(THETA_VALS) / 2)]
      rows <- .[abs(.$sapwood_theta_sat - med_theta) < 0.001, ]
      if (nrow(rows) > 1) diff(range(rows$kge, na.rm = TRUE)) else NA_real_
    },
    kge_range_theta = {
      med_vc <- VC_VALS[round(length(VC_VALS) / 2)]
      rows <- .[abs(.$volumetric_capacitance - med_vc) < 1, ]
      if (nrow(rows) > 1) diff(range(rows$kge, na.rm = TRUE)) else NA_real_
    },
    best_kge   = max(kge, na.rm = TRUE),
    best_vc    = volumetric_capacitance[which.max(kge)],
    best_theta = sapwood_theta_sat[which.max(kge)],
    .groups = "drop"
  ) %>%
  mutate(
    dominant_param = ifelse(kge_range_vc > kge_range_theta,
                            "volumetric_capacitance", "sapwood_theta_sat")
  )

# Add baseline KGE separately
baseline_kge <- metrics_twd %>%
  filter(run_tag == "baseline") %>%
  select(species, condition, baseline_kge = kge)

importance_summary <- importance_summary %>%
  left_join(baseline_kge, by = c("species", "condition")) %>%
  mutate(kge_improvement = best_kge - baseline_kge)

if (nrow(importance_summary) > 0) {
  print(as.data.frame(importance_summary))
  write.csv(importance_summary, file.path(OUTPUT_DIR, "parameter_importance_TWD.csv"), row.names = FALSE)
}

# ==========================================================================
# 14. FINAL RECOMMENDATIONS
# ==========================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  FINAL PARAMETER RECOMMENDATIONS (Water Storage)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("To apply best water storage parameters, update the .ins files with:\n\n")

if (nrow(best_twd) > 0) {
  for (i in seq_len(nrow(best_twd))) {
    r <- best_twd[i, ]
    cat(sprintf("  %-8s %-10s: volumetric_capacitance=%-5.0f sapwood_theta_sat=%-5.2f (KGE=%.3f, R²=%.3f)\n",
                r$species, r$condition,
                r$volumetric_capacitance, r$sapwood_theta_sat,
                r$kge, r$r_squared))
  }
}

cat("\n")
cat("  Dominant parameter per species×condition:\n")
if (nrow(importance_summary) > 0) {
  for (i in seq_len(nrow(importance_summary))) {
    r <- importance_summary[i, ]
    cat(sprintf("  %-8s %-10s: %s (KGE improvement over baseline: %+.3f)\n",
                r$species, r$condition, r$dominant_param, r$kge_improvement))
  }
}

cat("\n*** Water storage sensitivity analysis complete. Results saved to:", OUTPUT_DIR, "***\n")

