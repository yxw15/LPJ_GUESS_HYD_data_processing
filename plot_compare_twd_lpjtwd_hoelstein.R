# ==============================================================================
# 1. SETUP, INITIALIZATION & GLOBAL CONFIGURATIONS
# ==============================================================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(stringr)
library(purrr)
library(scales)

# Set working directory to the project root
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

# Establish structural directory targets
dir.create("Figures/lpj_guess_hyd_twd", recursive = TRUE, showWarnings = FALSE)

# Global Configuration Constants
species_levels <- c("oak", "beech", "spruce", "pine")
species_colors <- c("oak" = "#E69F00", "beech" = "#0072B2", "spruce" = "#009E73", "pine" = "#F0E442")

base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text       = element_text(color = "black", size = 11),
    legend.position   = "bottom",
    legend.box        = "vertical",
    plot.title        = element_text(hjust = 0.5, size = 14, face = "bold", color = "black"),
    plot.subtitle     = element_text(hjust = 0.5, size = 10, color = "grey30"),
    axis.title        = element_text(size = 12),
    axis.text.x       = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y       = element_text(angle = 0, hjust = 0.5, size = 10),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(size = 11, face = "bold")
  )

# ==============================================================================
# 2. DATA INGESTION & PIPELINE MERGING (PRESERVE FULL YEAR RAW RECORDS)
# ==============================================================================

# Load model output reference - PRESERVE ALL MONTHS
lpj_raw_twd <- read.csv("lpj_guess/lpj_control_drought_Gc_psiL_psiX_psiS_stem_diameter_twd_stem_rwc.csv") %>%
  mutate(date = as.Date(date), month = month(date)) %>%
  filter(!is.na(twd), !is.na(species), !is.na(treatment)) %>%
  mutate(
    species = tolower(species),
    twd_um_model = twd * 1e6
  ) %>%
  filter(species %in% species_levels) %>%
  mutate(species = factor(species, levels = species_levels)) %>%
  select(date, species, treatment, twd_um_model)

# Identify all field data archives
dendro_obs_files <- list.files(
  path = "SCCII/point_dendro", 
  pattern = "^Point_dendrometers_.*_archive\\.txt$",
  full.names = TRUE
)

dendro_obs <- dendro_obs_files %>%
  set_names() %>%
  map_dfr(~ read.delim(.x), .id = "source_file")

# Load and update treatment designations
tree_info <- read.csv("SCCII/tree_info.csv")

tree_info_updated <- tree_info %>%
  mutate(treatment = case_when(
    treatment == "treatment" ~ "drought",
    treatment == "control"   ~ "control",
    TRUE ~ treatment
  ))

# Raw observed aggregation step over ALL recorded seasonal timestamps
obs_twd_raw_aggregated <- dendro_obs %>%
  inner_join(tree_info_updated, by = "tree_id") %>%
  mutate(
    timestamp_clean = str_extract(timestamp_UTC, "\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}"),
    datetime = as.POSIXct(timestamp_clean, format = "%Y-%m-%d %H:%M:%S"),
    date = as.Date(datetime),
    month = month(date),
    species = tolower(species),
    treatment = tolower(str_extract(treatment, "^[a-zA-Z]+"))
  ) %>%
  filter(!is.na(datetime), !is.na(twd_micron_treenetproc), species %in% species_levels) %>%
  group_by(date, species, treatment) %>%
  summarise(
    twd_mean_obs  = mean(twd_micron_treenetproc, na.rm = TRUE),
    twd_top10_obs = mean(twd_micron_treenetproc[twd_micron_treenetproc >= quantile(twd_micron_treenetproc, 0.90, na.rm = TRUE)], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(species = factor(species, levels = species_levels))

# ==============================================================================
# 3. GENERATE DUAL CHRONOLOGICAL SETS (FULL YEAR VS SUMMER ONLY 6-9)
# ==============================================================================

# --- COHORT A: FULL YEAR TIME-SERIES ARRAYS ---
common_dates_fullyear <- inner_join(
  obs_twd_raw_aggregated %>% distinct(date, species, treatment),
  lpj_raw_twd %>% distinct(date, species, treatment),
  by = c("date", "species", "treatment")
)

daily_combined_fullyear <- common_dates_fullyear %>%
  inner_join(obs_twd_raw_aggregated, by = c("date", "species", "treatment")) %>%
  inner_join(lpj_raw_twd %>% select(date, species, treatment, twd_um_model), by = c("date", "species", "treatment")) %>%
  rename(sim_twd = twd_um_model)

monthly_combined_fullyear <- daily_combined_fullyear %>%
  mutate(month_date = floor_date(date, "month")) %>%
  group_by(month_date, species, treatment) %>%
  summarise(
    twd_mean_obs  = mean(twd_mean_obs, na.rm = TRUE),
    twd_top10_obs = mean(twd_top10_obs, na.rm = TRUE),
    sim_twd       = mean(sim_twd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(date = month_date)

# --- COHORT B: SUMMER SEASON ONLY TIME-SERIES ARRAYS (MONTHS 6 TO 9) ---
daily_combined_summer <- daily_combined_fullyear %>%
  filter(month(date) >= 6 & month(date) <= 9)

monthly_combined_summer <- monthly_combined_fullyear %>%
  filter(month(date) >= 6 & month(date) <= 9)

# ==============================================================================
# 4. TIME-SERIES PLOTTING FUNCTIONS (UNIFIED Y-AXIS RANGE & CHRONOLOGICAL AXIS)
# ==============================================================================

# Unified Single Y-Axis Timeline Plot Function Engine
plot_unified_axis_twd <- function(plot_data, time_scale_name, seasonal_name, date_breaks_config, date_labels_config, is_summer = FALSE) {
  processed_axis_data <- plot_data %>%
    rename(
      target_mean  = any_of(c("twd_mean_obs", "obs_mean")),
      target_top10 = any_of(c("twd_top10_obs", "obs_top10")),
      target_sim   = any_of(c("sim_twd"))
    ) %>%
    mutate(year = factor(year(date))) # Extract year for faceting
  
  # Dynamically calculate the global maximum across all lines to set unified limits
  global_max_y <- max(c(processed_axis_data$target_mean, processed_axis_data$target_top10, processed_axis_data$target_sim), na.rm = TRUE) * 1.05
  
  p <- ggplot(processed_axis_data) +
    geom_line(aes(x = date, y = target_mean, color = "observed daily mean"), linewidth = 0.8) +
    geom_line(aes(x = date, y = target_top10, color = "observed top 10% mean"), linewidth = 0.8) +
    geom_line(aes(x = date, y = target_sim, color = species), linewidth = 1.1)
  
  # Conditional faceting: Split by year if processing the summer subsets
  if (is_summer) {
    p <- p + facet_grid(species ~ treatment + year, scales = "free_x")
  } else {
    p <- p + facet_grid(species ~ treatment)
  }
  
  p <- p +
    scale_y_continuous(
      name = "tree water deficit (\u00b5m)",
      limits = c(0, global_max_y),
      expand = c(0, 0)
    ) +
    scale_color_manual(
      name = "data tracking lines",
      values = c("observed daily mean" = "black", "observed top 10% mean" = "grey60", species_colors),
      breaks = c("observed daily mean", "observed top 10% mean", species_levels)
    ) +
    scale_x_date(date_breaks = date_breaks_config, date_labels = date_labels_config) +
    labs(
      title = paste0("tree water deficit: ", tolower(time_scale_name), " (", tolower(seasonal_name), ")"), 
      subtitle = "unified axis calibration: observed data tracking vs lpj-guess ecosystem simulation", 
      x = "timeline"
    ) +
    base_theme
  
  # Clean up duplicate/cluttered x-axis labels if split into year panels
  if (is_summer) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  return(p)
}

# B. Relative Scale Trend Synchrony Function Engine (0-100% Normalized Alone per line/panel)
plot_trend_synchrony <- function(plot_data, time_scale_name, seasonal_name, date_breaks_config, date_labels_config, is_summer = FALSE) {
  processed_sync_data <- plot_data %>%
    rename(
      target_mean  = any_of(c("twd_mean_obs", "obs_mean")),
      target_top10 = any_of(c("twd_top10_obs", "obs_top10")),
      target_sim   = any_of(c("sim_twd"))
    ) %>%
    mutate(year = factor(year(date)))
  
  panel_scaled_data <- processed_sync_data %>%
    group_by(species, treatment) %>%
    mutate(
      min_mean  = if(all(is.na(target_mean)))  NA_real_ else min(target_mean, na.rm = TRUE),  
      max_mean  = if(all(is.na(target_mean)))  NA_real_ else max(target_mean, na.rm = TRUE),
      min_top10 = if(all(is.na(target_top10))) NA_real_ else min(target_top10, na.rm = TRUE), 
      max_top10 = if(all(is.na(target_top10))) NA_real_ else max(target_top10, na.rm = TRUE),
      min_sim   = if(all(is.na(target_sim)))   NA_real_ else min(target_sim, na.rm = TRUE),        
      max_sim   = if(all(is.na(target_sim)))   NA_real_ else max(target_sim, na.rm = TRUE),
      
      obs_norm   = if_else(!is.na(max_mean)  & max_mean > min_mean,   (target_mean - min_mean) / (max_mean - min_mean), 0),
      top10_norm = if_else(!is.na(max_top10) & max_top10 > min_top10, (target_top10 - min_top10) / (max_top10 - min_top10), 0),
      sim_norm   = if_else(!is.na(max_sim)   & max_sim > min_sim,     (target_sim - min_sim) / (max_sim - min_sim), 0)
    ) %>%
    ungroup() %>%
    filter(!is.na(obs_norm), !is.na(top10_norm), !is.na(sim_norm)) %>%
    mutate(species = factor(tolower(species), levels = tolower(species_levels)))
  
  lowercase_colors <- species_colors
  names(lowercase_colors) <- tolower(names(lowercase_colors))
  
  p <- ggplot(panel_scaled_data) +
    geom_line(aes(x = date, y = obs_norm, color = "observed daily mean"), linewidth = 0.8) +
    geom_line(aes(x = date, y = top10_norm, color = "observed top 10% mean"), linewidth = 0.8) +
    geom_line(aes(x = date, y = sim_norm, color = species), linewidth = 1.1) +
    scale_y_continuous(
      name = "relative line metric range (0% to 100% normalized alone)",
      labels = scales::percent_format(accuracy = 1)
    ) +
    scale_color_manual(
      name = "trend comparison lines",
      values = c("observed daily mean" = "black", "observed top 10% mean" = "grey60", lowercase_colors),
      breaks = c("observed daily mean", "observed top 10% mean", tolower(species_levels))
    ) +
    scale_x_date(date_breaks = date_breaks_config, date_labels = date_labels_config) +
    labs(
      title = paste0("tree water deficit trend synchrony: ", tolower(time_scale_name), " (", tolower(seasonal_name), ")"), 
      subtitle = "each line is independently standardized per species and treatment panel (0% = panel min, 100% = panel max)", 
      x = "timeline"
    ) +
    base_theme
  
  if (is_summer) {
    p <- p + facet_grid(species ~ treatment + year, scales = "free_x") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    p <- p + facet_grid(species ~ treatment)
  }
  
  return(p)
}

# --- SAVE PLOTS FOR FULL YEAR ARCHIVES ---
ggsave("Figures/lpj_guess_hyd_twd/twd_unified_axis_daily_fullyear.png", plot_unified_axis_twd(daily_combined_fullyear, "Daily", "Full Year", "1 year", "%Y", is_summer = FALSE), width = 16, height = 11, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/twd_unified_axis_monthly_fullyear.png", plot_unified_axis_twd(monthly_combined_fullyear, "Monthly", "Full Year", "1 year", "%Y", is_summer = FALSE), width = 16, height = 11, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/twd_trend_synchrony_daily_fullyear.png", plot_trend_synchrony(daily_combined_fullyear, "Daily", "Full Year", "1 year", "%Y", is_summer = FALSE), width = 16, height = 11, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/twd_trend_synchrony_monthly_fullyear.png", plot_trend_synchrony(monthly_combined_fullyear, "Monthly", "Full Year", "1 year", "%Y", is_summer = FALSE), width = 16, height = 11, dpi = 300)

# --- SAVE PLOTS FOR SUMMER ONLY WINDOWS (6-9) ---
ggsave("Figures/lpj_guess_hyd_twd/twd_unified_axis_daily_summer.png", plot_unified_axis_twd(daily_combined_summer, "Daily", "Summer Months 6-9", "1 year", "%Y", is_summer = TRUE), width = 18, height = 11, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/twd_unified_axis_monthly_summer.png", plot_unified_axis_twd(monthly_combined_summer, "Monthly", "Summer Months 6-9", "1 year", "%Y", is_summer = TRUE), width = 18, height = 11, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/twd_trend_synchrony_daily_summer.png", plot_trend_synchrony(daily_combined_summer, "Daily", "Summer Months 6-9", "1 year", "%Y", is_summer = TRUE), width = 18, height = 11, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/twd_trend_synchrony_monthly_summer.png", plot_trend_synchrony(monthly_combined_summer, "Monthly", "Summer Months 6-9", "1 year", "%Y", is_summer = TRUE), width = 18, height = 11, dpi = 300)
# ==============================================================================
# 5. DATA PAIR COHORT RESTRUCTURING ENGINE
# ==============================================================================
pair_levels_twd <- c("twd_vs_mean", "twd_vs_top10")

generate_flat_pairs <- function(source_daily_df) {
  bind_rows(
    source_daily_df %>% transmute(date, species, treatment, pair_id = "twd_vs_mean",  model_val = sim_twd, obs_val = twd_mean_obs),
    source_daily_df %>% transmute(date, species, treatment, pair_id = "twd_vs_top10", model_val = sim_twd, obs_val = twd_top10_obs)
  ) %>%
    filter(!is.na(model_val) & !is.na(obs_val)) %>%
    mutate(
      species = factor(species, levels = species_levels),
      pair_id = factor(pair_id, levels = pair_levels_twd),
      diff = model_val - obs_val
    )
}

flat_pairs_fullyear <- generate_flat_pairs(daily_combined_fullyear)
flat_pairs_summer   <- generate_flat_pairs(daily_combined_summer)

# ==============================================================================
# 6. CALCULATE SCATTER STATISTICS TRACKS PER SCENARIO
# ==============================================================================
calculate_scatter_statistics <- function(flat_paired_df) {
  flat_paired_df %>%
    group_by(species, treatment, pair_id) %>%
    summarise(
      n = n(),
      pearson_r = cor(obs_val, model_val, use = "complete.obs"),
      pearson_r2 = cor(obs_val, model_val, use = "complete.obs")^2,
      rmse = sqrt(mean((model_val - obs_val)^2, na.rm = TRUE)),
      nrmse = (sqrt(mean((model_val - obs_val)^2, na.rm = TRUE)) / mean(obs_val, na.rm = TRUE)) * 100,
      bias = mean(model_val - obs_val, na.rm = TRUE),
      slope = if(n() > 1 && sd(obs_val, na.rm=TRUE) > 0) coef(lm(model_val ~ obs_val))[2] else NA,
      .groups = "drop"
    )
}

stats_fullyear <- calculate_scatter_statistics(flat_pairs_fullyear)
stats_summer   <- calculate_scatter_statistics(flat_pairs_summer)

write.csv(stats_fullyear, "Figures/lpj_guess_hyd_twd/scatter_statistics_twd_2pairs_fullyear.csv", row.names = FALSE)
write.csv(stats_summer, "Figures/lpj_guess_hyd_twd/scatter_statistics_twd_2pairs_summer.csv", row.names = FALSE)

# ==============================================================================
# 7. GENERATE ANNOTATED ISOMETRIC 1:1 SCATTER PLOTS
# ==============================================================================
create_scatter_twd_plot <- function(flat_paired_df, stats_df, target_treatment, plot_title, subtitle_extension) {
  data_subset <- flat_paired_df %>% filter(treatment == target_treatment)
  stats_subset <- stats_df %>% filter(treatment == target_treatment)
  
  pair_labels <- c(
    "twd_vs_mean"  = "model twd vs obs mean",
    "twd_vs_top10" = "model twd vs obs top 10% mean"
  )
  
  data_subset <- data_subset %>% mutate(pair_label = pair_labels[pair_id])
  
  annotation_data <- stats_subset %>%
    mutate(
      pair_label = pair_labels[pair_id],
      text_summary = paste0(
        "n = ", n, "\n",
        "r = ", round(pearson_r, 2), "\n",
        "r\u00b2 = ", round(pearson_r2, 2), "\n",
        "rmse = ", round(rmse, 2), "\n",
        "bias = ", round(bias, 2), "\n",
        "slope = ", round(slope, 1)
      )
    )
  
  axis_max <- max(c(data_subset$obs_val, data_subset$model_val), na.rm = TRUE) * 1.05
  
  ggplot(data_subset, aes(x = obs_val, y = model_val, color = species)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.1, linewidth = 0.6) +
    facet_grid(pair_label ~ species) +
    scale_color_manual(values = species_colors) +
    coord_fixed(ratio = 1, xlim = c(0, axis_max), ylim = c(0, axis_max)) +
    geom_text(
      data = annotation_data,
      aes(x = Inf, y = -Inf, label = text_summary),
      hjust = 1.05, vjust = -0.1, size = 2.8, color = "black", inherit.aes = FALSE
    ) +
    labs(
      title = tolower(plot_title),
      subtitle = paste0("dashed line = 1:1 identity framework | ", tolower(subtitle_extension)),
      x = "observed tree water deficit (\u00b5m)",
      y = "lpj-guess simulated tree water deficit (\u00b5m)"
    ) +
    base_theme +
    theme(aspect.ratio = 1, strip.text.y = element_text(angle = -90, size = 9))
}

# --- GENERATE AND EXPORT COHORT SCATTERS ---
scatter_control_fullyear <- create_scatter_twd_plot(flat_pairs_fullyear, stats_fullyear, "control", "tree water deficit 1:1 scatter: control (full year)", "full year matched records")
scatter_drought_fullyear <- create_scatter_twd_plot(flat_pairs_fullyear, stats_fullyear, "drought", "tree water deficit 1:1 scatter: drought (full year)", "full year matched records")
scatter_control_summer   <- create_scatter_twd_plot(flat_pairs_summer, stats_summer, "control", "tree water deficit 1:1 scatter: control (summer 6-9)", "june-september records only")
scatter_drought_summer   <- create_scatter_twd_plot(flat_pairs_summer, stats_summer, "drought", "tree water deficit 1:1 scatter: drought (summer 6-9)", "june-september records only")

ggsave("Figures/lpj_guess_hyd_twd/scatter_control_twd_2pairs_fullyear.png", scatter_control_fullyear, width = 14, height = 8, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/scatter_drought_twd_2pairs_fullyear.png", scatter_drought_fullyear, width = 14, height = 8, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/scatter_control_twd_2pairs_summer.png",   scatter_control_summer,   width = 14, height = 8, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/scatter_drought_twd_2pairs_summer.png",   scatter_drought_summer,   width = 14, height = 8, dpi = 300)

# ==============================================================================
# 8. METADATA STATISTICAL DICTIONARY FILE GENERATION
# ==============================================================================
cat("\n=== logging metadata explanations ===\n")
sink("Figures/lpj_guess_hyd_twd/statistical_variables_explanation_twd.txt")
cat("================================================================================\n")
cat("statistical variables explanation dictionary\n")
cat("tree water deficit evaluation matrix: field observations vs lpj-guess simulated\n")
cat("================================================================================\n\n")
cat("n (sample cohort size): paired strict chronological date match records.\n")
cat("pearson_r (pearson correlation): direction and linear strength (-1 to +1).\n")
cat("pearson_r2 (coefficient of determination): explained variation proportion (0 to 1).\n")
cat("rmse (root mean square error): absolute quadratic fit penalty metric (micrometers).\n")
cat("bias (mean directional bias): mean delta calculation (simulated - observed) (micrometers).\n")
cat("slope (fitted linear slope): regression slope. ideal matching target = 1.0.\n")
sink()
cat("\u2713 statistical dictionary written successfully.\n")