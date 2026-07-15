# ==========================================================================
# 1. SETUP: GLOBAL CONSTANTS & THEMES (AUGUST VALIDATION ONLY)
# ==========================================================================
library(dplyr)
library(ggplot2)
library(lubridate)
library(scales)
library(tidyr)

VALIDATION_MONTH <- 8  # August only

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

dir.create("Figures/lpj_guess_stem_storage/validation_august/Gc", recursive = TRUE, showWarnings = FALSE)

species_levels <- c("oak", "beech", "spruce", "pine")
species_colors <- c("oak" = "#E69F00", "beech" = "#0072B2", "spruce" = "#009E73", "pine" = "#F0E442")

base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text       = element_text(color = "black", size = 12),
    legend.position   = "bottom",
    plot.title        = element_text(hjust = 0.5, size = 16, color = "black"),
    axis.title        = element_text(size = 14),
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y       = element_text(angle = 0, hjust = 0.5, size = 10),
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_line(color = "grey92", linewidth = 0.25),
    panel.border      = element_blank(),
    strip.text        = element_text(size = 12, face = "bold")
  )

# ==========================================================================
# 2. DATA PREPARATION (AUGUST ONLY)
# ==========================================================================
lpj_output_filter <- read.csv("lpj_guess/lpj_guess_stem_storage/lpj_control_drought_ET_Gc_psi_leaf_psi_soil_psi_xylem_hydraulic_lag_kappy_s_min_mort_mort_cav_mort_greff_mort_min_stem_diameter_stem_rwc_twd_august_climate_filter.csv")
sap_flux_gc_filter <- read.csv("SCCII/sap_flux_gc_daytime_Climate_filter.csv")

sap_flux_gc_filter <- sap_flux_gc_filter %>%
  mutate(treatment = case_when(
    treatment == "control" ~ "control",
    treatment == "treatment" ~ "drought"
  ))

# ==========================================================================
# 3. DATA PROCESSING: AUGUST ONLY
# ==========================================================================

observed_data <- sap_flux_gc_filter %>%
  mutate(date = as.Date(date),
         species = tolower(species)) %>%
  filter(month(date) == VALIDATION_MONTH) %>%   # AUGUST ONLY
  filter(!is.na(G_ms)) %>%
  filter(!is.na(species)) %>%
  filter(!is.na(treatment)) %>%
  filter(G_ms <= 12) %>%
  mutate(species = factor(species, levels = species_levels))

lpj_processed <- lpj_output_filter %>%
  mutate(date = as.Date(date),
         species = tolower(species)) %>%
  filter(month(date) == VALIDATION_MONTH) %>%   # AUGUST ONLY
  filter(!is.na(Gc)) %>%
  filter(!is.na(species)) %>%
  filter(!is.na(treatment)) %>%
  filter(Gc <= 12) %>%
  mutate(species = factor(species, levels = species_levels)) %>%
  rename(G_cond_mean = Gc)

# Find common dates
observed_dates <- observed_data %>%
  distinct(treatment, species, date)
lpj_dates <- lpj_processed %>%
  distinct(treatment, species, date)
common_dates <- inner_join(observed_dates, lpj_dates,
                           by = c("treatment", "species", "date"))

cat("\n=== August validation: common dates ===\n")
cat("total common date-treatment-species combinations:", nrow(common_dates), "\n")

observed_data_filtered <- observed_data %>%
  inner_join(common_dates, by = c("treatment", "species", "date"))
lpj_processed_filtered <- lpj_processed %>%
  inner_join(common_dates, by = c("treatment", "species", "date"))

start_date <- min(common_dates$date)
end_date <- max(common_dates$date)
cat("date range:", as.character(start_date), "to", as.character(end_date), "\n")

# ==========================================================================
# 4. OBSERVED DATA AGGREGATION
# ==========================================================================
observed_daytime_mean <- observed_data_filtered %>%
  group_by(date, species, treatment) %>%
  summarise(G_obs_mean = mean(G_ms, na.rm = TRUE), .groups = "drop")

observed_top10 <- observed_data_filtered %>%
  group_by(date, species, treatment) %>%
  mutate(q90 = quantile(G_ms, 0.90, na.rm = TRUE)) %>%
  filter(G_ms >= q90) %>%
  summarise(G_obs_top10 = mean(G_ms, na.rm = TRUE), .groups = "drop")

# ==========================================================================
# 5. TIME SERIES PLOTS
# ==========================================================================

create_timeseries_plot <- function(treatment_name, plot_title) {
  lpj_subset <- lpj_processed_filtered %>%
    filter(treatment == treatment_name) %>%
    mutate(year = year(date))
  obs_daytime_subset <- observed_daytime_mean %>%
    filter(treatment == treatment_name) %>%
    mutate(year = year(date))
  obs_top10_subset <- observed_top10 %>%
    filter(treatment == treatment_name) %>%
    mutate(year = year(date))
  treatment_dates <- common_dates %>% filter(treatment == treatment_name)
  t_start <- min(treatment_dates$date)
  t_end <- max(treatment_dates$date)

  ggplot() +
    geom_point(data = lpj_subset,
               aes(x = date, y = G_cond_mean, color = species),
               alpha = 1, size = 1.2) +
    geom_point(data = obs_daytime_subset,
               aes(x = date, y = G_obs_mean, shape = "observed daytime mean"),
               color = "grey70", alpha = 0.7, size = 1.8, stroke = 1) +
    geom_point(data = obs_top10_subset,
               aes(x = date, y = G_obs_top10, shape = "observed top 10%"),
               color = "grey20", alpha = 0.7, size = 1.8, stroke = 1) +
    facet_grid(species ~ year, scales = "free") +
    scale_color_manual(name = "lpj-guess simulation", values = species_colors) +
    scale_shape_manual(name = "field measurements (sccii)",
                       values = c("observed daytime mean" = 2, "observed top 10%" = 1)) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b") +
    coord_cartesian(ylim = c(0, 12)) +
    labs(title = tolower(plot_title),
         subtitle = paste0("AUGUST VALIDATION | dates: ", t_start, " to ", t_end),
         x = "timeline", y = "conductance (m s-1)") +
    base_theme +
    theme(legend.position = "bottom", strip.text = element_text(size = 12, face = "bold"),
          legend.box = "vertical", axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10))
}

create_timeseries_plot_line <- function(treatment_name, plot_title) {
  lpj_subset <- lpj_processed_filtered %>%
    filter(treatment == treatment_name) %>%
    mutate(year = year(date))
  obs_daytime_subset <- observed_daytime_mean %>%
    filter(treatment == treatment_name) %>%
    mutate(year = year(date))
  treatment_dates <- common_dates %>% filter(treatment == treatment_name)

  ggplot() +
    geom_line(data = lpj_subset,
              aes(x = date, y = G_cond_mean, color = species),
              linetype = "solid", linewidth = 1, alpha = 1) +
    geom_line(data = obs_daytime_subset,
              aes(x = date, y = G_obs_mean),
              color = "grey20", linetype = "dashed", linewidth = 0.8, alpha = 1) +
    facet_grid(species ~ year, scales = "free") +
    scale_color_manual(name = "lpj-guess simulation", values = species_colors) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b") +
    coord_cartesian(ylim = c(0, 12)) +
    labs(title = tolower(plot_title),
         subtitle = paste0("AUGUST VALIDATION | dates: ", min(treatment_dates$date), " to ", max(treatment_dates$date)),
         x = "timeline", y = "conductance (m s-1)") +
    base_theme +
    theme(legend.position = "bottom", strip.text = element_text(size = 12, face = "bold"),
          legend.box = "vertical", axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10))
}

plot_control_ts <- create_timeseries_plot("control", "canopy conductance: control [AUGUST]")
plot_drought_ts <- create_timeseries_plot("drought", "canopy conductance: drought [AUGUST]")
plot_control_ts_line <- create_timeseries_plot_line("control", "canopy conductance: control [AUGUST]")
plot_drought_ts_line <- create_timeseries_plot_line("drought", "canopy conductance: drought [AUGUST]")

# ==========================================================================
# 6. SCATTER PLOTS
# ==========================================================================

scatter_data_daytime <- observed_daytime_mean %>%
  inner_join(lpj_processed_filtered, by = c("date", "species", "treatment")) %>%
  rename(g_sim = G_cond_mean, g_obs_mean = G_obs_mean)

scatter_data_top10 <- observed_top10 %>%
  inner_join(lpj_processed_filtered, by = c("date", "species", "treatment")) %>%
  rename(g_sim = G_cond_mean, g_obs_mean = G_obs_top10)

safe_cor <- function(x, y) {
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) NA_real_ else cor(x, y, use = "complete.obs")
}
safe_slope <- function(x, y) {
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) NA_real_ else coef(lm(y ~ x))[2]
}

scatter_stats_daytime <- scatter_data_daytime %>%
  group_by(species, treatment) %>%
  summarise(
    n = n(), pearson_r = safe_cor(g_obs_mean, g_sim),
    pearson_r2 = safe_cor(g_obs_mean, g_sim)^2,
    rmse = sqrt(mean((g_sim - g_obs_mean)^2, na.rm = TRUE)),
    nrmse = (sqrt(mean((g_sim - g_obs_mean)^2, na.rm = TRUE)) / mean(g_obs_mean, na.rm = TRUE)) * 100,
    bias = mean(g_sim - g_obs_mean, na.rm = TRUE),
    slope = safe_slope(g_obs_mean, g_sim),
    .groups = "drop"
  )

scatter_stats_top10 <- scatter_data_top10 %>%
  group_by(species, treatment) %>%
  summarise(
    n = n(), pearson_r = safe_cor(g_obs_mean, g_sim),
    pearson_r2 = safe_cor(g_obs_mean, g_sim)^2,
    rmse = sqrt(mean((g_sim - g_obs_mean)^2, na.rm = TRUE)),
    nrmse = (sqrt(mean((g_sim - g_obs_mean)^2, na.rm = TRUE)) / mean(g_obs_mean, na.rm = TRUE)) * 100,
    bias = mean(g_sim - g_obs_mean, na.rm = TRUE),
    slope = safe_slope(g_obs_mean, g_sim),
    .groups = "drop"
  )

create_scatter_plot <- function(data, stats, treatment_name, plot_title, max_limit = 12) {
  data_subset <- data %>% filter(treatment == treatment_name)
  stats_subset <- stats %>% filter(treatment == treatment_name)

  annotation_data <- stats_subset %>%
    mutate(label = paste0(
      "n = ", n, "\n", "r = ", round(pearson_r, 3), "\n",
      "r² = ", round(pearson_r2, 3), "\n",
      "rmse = ", round(rmse, 4), "\n", "nrmse = ", round(nrmse, 1), "%\n",
      "bias = ", round(bias, 4), "\n", "slope = ", round(slope, 2)
    ))

  ggplot(data_subset, aes(x = g_obs_mean, y = g_sim, color = species)) +
    geom_point(alpha = 0.6, size = 2.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 0.8) +
    facet_wrap(~ species, ncol = 2) +
    geom_text(data = annotation_data,
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.05, vjust = 1.1, size = 3.2, color = "black", inherit.aes = FALSE) +
    scale_color_manual(values = species_colors, name = "species") +
    coord_fixed(ratio = 1, xlim = c(0, max_limit), ylim = c(0, max_limit)) +
    labs(title = tolower(plot_title),
         subtitle = "AUGUST VALIDATION | 1:1 line (dashed), linear regression (solid)",
         x = "observed conductance (m s-1)", y = "lpj simulated conductance (m s-1)") +
    base_theme +
    theme(legend.position = "bottom", strip.text = element_text(size = 11, face = "bold"))
}

scatter_control_daytime <- create_scatter_plot(scatter_data_daytime, scatter_stats_daytime,
                                               "control", "Gc scatter (daytime mean): control [AUGUST]")
scatter_drought_daytime <- create_scatter_plot(scatter_data_daytime, scatter_stats_daytime,
                                               "drought",  "Gc scatter (daytime mean): drought [AUGUST]")
scatter_control_top10 <- create_scatter_plot(scatter_data_top10, scatter_stats_top10,
                                             "control", "Gc scatter (top 10%): control [AUGUST]")
scatter_drought_top10 <- create_scatter_plot(scatter_data_top10, scatter_stats_top10,
                                             "drought",  "Gc scatter (top 10%): drought [AUGUST]")

# ==========================================================================
# 7. EXPORT
# ==========================================================================
out_dir <- "Figures/lpj_guess_stem_storage/validation_august/Gc"

# Time series
ggsave(file.path(out_dir, "timeseries_control.png"), plot_control_ts, width = 16, height = 11, dpi = 300)
ggsave(file.path(out_dir, "timeseries_drought.png"), plot_drought_ts, width = 16, height = 11, dpi = 300)
ggsave(file.path(out_dir, "timeseries_control_lines.png"), plot_control_ts_line, width = 16, height = 11, dpi = 300)
ggsave(file.path(out_dir, "timeseries_drought_lines.png"), plot_drought_ts_line, width = 16, height = 11, dpi = 300)

# Scatter
ggsave(file.path(out_dir, "scatter_control_daytime_mean.png"), scatter_control_daytime, width = 11, height = 9, dpi = 300)
ggsave(file.path(out_dir, "scatter_drought_daytime_mean.png"), scatter_drought_daytime, width = 11, height = 9, dpi = 300)
ggsave(file.path(out_dir, "scatter_control_top10.png"), scatter_control_top10, width = 11, height = 9, dpi = 300)
ggsave(file.path(out_dir, "scatter_drought_top10.png"), scatter_drought_top10, width = 11, height = 9, dpi = 300)

# Stats
write.csv(scatter_stats_daytime, file.path(out_dir, "scatter_statistics_daytime_mean_august.csv"), row.names = FALSE)
write.csv(scatter_stats_top10, file.path(out_dir, "scatter_statistics_top10_august.csv"), row.names = FALSE)

cat("\n*** August Gc validation complete ***\n")
cat("Plots saved to:", out_dir, "\n")
