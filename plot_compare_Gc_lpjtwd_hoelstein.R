# ==========================================================================
# 1. setup: global constants & themes
# ==========================================================================
library(dplyr)
library(ggplot2)
library(lubridate)
library(scales)
library(tidyr)

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

# create directory if it doesn't exist
dir.create("Figures/compare_Gc_lpjtwd_hoelstein", recursive = TRUE, showWarnings = FALSE)

# aesthetics
species_levels <- c("Oak", "Beech", "Spruce", "Pine")
species_colors <- c("Oak" = "#E69F00", "Beech" = "#0072B2", "Spruce" = "#009E73", "Pine" = "#F0E442")

base_theme <- theme_minimal() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(color = "black", size = 12),
    legend.position = "bottom",
    plot.title  = element_text(hjust = 0.5, size = 16, color = "black"),
    axis.title  = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 10),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    panel.border = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )

# ==========================================================================
# 2. data preparation & filter definition
# ==========================================================================
lpj_output_filter <- read.csv("lpj_guess/lpj_control_drought_Gc_psiL_psiX_psiS_stem_diameter_twd_stem_rwc_climate_filter.csv")
sap_flux_gc_filter <- read.csv("SCCII/sap_daily_filter.csv")

sap_flux_gc_filter <- sap_flux_gc_filter %>%
  mutate(treatment = case_when(
    treatment == "control" ~ "control",
    treatment == "treatment" ~ "drought"
  ))

# ==========================================================================
# 3. data processing: find common dates (exact date matching) and exclude >12
# ==========================================================================

# process observed data and exclude G_ms > 12
observed_data <- sap_flux_gc_filter %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= as.Date("2023-01-01")) %>%
  filter(!is.na(G_ms)) %>%
  filter(!is.na(species)) %>%
  filter(!is.na(treatment)) %>%
  filter(G_ms <= 12) %>%
  mutate(species = factor(species, levels = species_levels))

# process lpj data and exclude Gc > 12
lpj_processed <- lpj_output_filter %>%
  mutate(date = as.Date(date)) %>%
  filter(date >= as.Date("2023-01-01")) %>%
  filter(!is.na(Gc)) %>%
  filter(!is.na(species)) %>%
  filter(!is.na(treatment)) %>%
  filter(Gc <= 12) %>%
  mutate(species = factor(species, levels = species_levels)) %>%
  rename(G_cond_mean = Gc)

# find common dates
observed_dates <- observed_data %>%
  distinct(treatment, species, date)

lpj_dates <- lpj_processed %>%
  distinct(treatment, species, date)

common_dates <- inner_join(observed_dates, lpj_dates, 
                           by = c("treatment", "species", "date"))

cat("\n=== common dates (exact matching, values <=12) ===\n")
cat("total common date-treatment-species combinations:", nrow(common_dates), "\n")

# filter both datasets to only include common dates
observed_data_filtered <- observed_data %>%
  inner_join(common_dates, by = c("treatment", "species", "date"))

lpj_processed_filtered <- lpj_processed %>%
  inner_join(common_dates, by = c("treatment", "species", "date"))

# get date range
start_date <- min(common_dates$date)
end_date <- max(common_dates$date)
cat("\ndate range:", as.character(start_date), "to", as.character(end_date), "\n")

# ==========================================================================
# 4. data processing: observed data (both methods)
# ==========================================================================

# method 1: simple daily mean
observed_daily_mean <- observed_data_filtered %>%
  group_by(date, species, treatment) %>%
  summarise(G_obs_mean = mean(G_ms, na.rm = TRUE), 
            .groups = "drop")

# method 2: daily top 10% quantile mean
observed_top10 <- observed_data_filtered %>%
  group_by(date, species, treatment) %>%
  mutate(q90 = quantile(G_ms, 0.90, na.rm = TRUE)) %>%
  filter(G_ms >= q90) %>%
  summarise(G_obs_top10 = mean(G_ms, na.rm = TRUE), 
            .groups = "drop")

# ==========================================================================
# 5. TIME SERIES POINT SCATTER PLOTS
# ==========================================================================

# function to create time series plot
create_timeseries_plot <- function(treatment_name, plot_title) {
  
  lpj_subset <- lpj_processed_filtered %>% filter(treatment == treatment_name)
  obs_daily_subset <- observed_daily_mean %>% filter(treatment == treatment_name)
  obs_top10_subset <- observed_top10 %>% filter(treatment == treatment_name)
  
  observed_combined <- bind_rows(
    obs_daily_subset %>% mutate(obs_type = "daily mean", G_value = G_obs_mean),
    obs_top10_subset %>% mutate(obs_type = "top 10%", G_value = G_obs_top10)
  )
  
  treatment_dates <- common_dates %>% filter(treatment == treatment_name)
  t_start <- min(treatment_dates$date)
  t_end <- max(treatment_dates$date)
  
  p <- ggplot() +
    geom_point(data = lpj_subset,
               aes(x = date, y = G_cond_mean, color = species),
               alpha = 0.7, size = 1.5) +
    geom_point(data = observed_combined,
               aes(x = date, y = G_value, shape = obs_type),
               color = "grey40", alpha = 0.6, size = 1.5, stroke = 1) +
    facet_wrap(~ species, ncol = 2, scales = "free_y") +
    scale_color_manual(name = "lpj-guess (colored)", values = species_colors) +
    scale_shape_manual(name = "observed (grey)",
                       values = c("daily mean" = 1, "top 10%" = 2),
                       labels = c("daily mean" = "daily mean (open circle)", 
                                  "top 10%" = "top 10% quantile (open triangle)")) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
    coord_cartesian(ylim = c(0, 12)) +
    labs(title = plot_title,
         subtitle = paste0("common dates: ", as.character(t_start), " to ", as.character(t_end)),
         x = "timeline",
         y = expression(conductance ~ (m ~ s^{-1}))) +
    base_theme +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 12, face = "bold"),
          legend.box = "vertical")
  
  return(p)
}

# create time series plots
plot_control_ts <- create_timeseries_plot("control", "canopy conductance: control treatment")
plot_drought_ts <- create_timeseries_plot("drought", "canopy conductance: drought treatment")

# ==========================================================================
# 6. prepare data for scatter plots
# ==========================================================================

# daily mean scatter data
scatter_data_daily <- observed_daily_mean %>%
  inner_join(lpj_processed_filtered, by = c("date", "species", "treatment")) %>%
  rename(g_sim = G_cond_mean, g_obs_mean = G_obs_mean)

# top 10% scatter data
scatter_data_top10 <- observed_top10 %>%
  inner_join(lpj_processed_filtered, by = c("date", "species", "treatment")) %>%
  rename(g_sim = G_cond_mean, g_obs_mean = G_obs_top10)

# ==========================================================================
# 7. calculate statistics PER SPECIES AND PER TREATMENT
# ==========================================================================

scatter_stats_daily <- scatter_data_daily %>%
  group_by(species, treatment) %>%
  summarise(
    n = n(),
    pearson_r = cor(g_obs_mean, g_sim, use = "complete.obs"),
    pearson_r2 = cor(g_obs_mean, g_sim, use = "complete.obs")^2,
    rmse = sqrt(mean((g_sim - g_obs_mean)^2, na.rm = TRUE)),
    nrmse = (sqrt(mean((g_sim - g_obs_mean)^2, na.rm = TRUE)) / mean(g_obs_mean, na.rm = TRUE)) * 100,
    bias = mean(g_sim - g_obs_mean, na.rm = TRUE),
    slope = coef(lm(g_sim ~ g_obs_mean))[2],
    .groups = "drop"
  )

scatter_stats_top10 <- scatter_data_top10 %>%
  group_by(species, treatment) %>%
  summarise(
    n = n(),
    pearson_r = cor(g_obs_mean, g_sim, use = "complete.obs"),
    pearson_r2 = cor(g_obs_mean, g_sim, use = "complete.obs")^2,
    rmse = sqrt(mean((g_sim - g_obs_mean)^2, na.rm = TRUE)),
    nrmse = (sqrt(mean((g_sim - g_obs_mean)^2, na.rm = TRUE)) / mean(g_obs_mean, na.rm = TRUE)) * 100,
    bias = mean(g_sim - g_obs_mean, na.rm = TRUE),
    slope = coef(lm(g_sim ~ g_obs_mean))[2],
    .groups = "drop"
  )

# ==========================================================================
# 8. SCATTER PLOTS with CORRECT per-species statistics
# ==========================================================================

create_scatter_plot <- function(data, stats, treatment_name, plot_title, max_limit = 12) {
  
  data_subset <- data %>% filter(treatment == treatment_name)
  stats_subset <- stats %>% filter(treatment == treatment_name)
  
  annotation_data <- stats_subset %>%
    mutate(
      label = paste0(
        "n = ", n, "\n",
        "r = ", round(pearson_r, 3), "\n",
        "r² = ", round(pearson_r2, 3), "\n",
        "rmse = ", round(rmse, 4), "\n",
        "nrmse = ", round(nrmse, 1), "%\n",
        "bias = ", round(bias, 4), "\n",
        "slope = ", round(slope, 2)
      )
    )
  
  p <- ggplot(data_subset, aes(x = g_obs_mean, y = g_sim, color = species)) +
    geom_point(alpha = 0.6, size = 2.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 0.8) +
    facet_wrap(~ species, ncol = 2) +
    scale_color_manual(values = species_colors, name = "species") +
    coord_fixed(ratio = 1, xlim = c(0, max_limit), ylim = c(0, max_limit)) +
    labs(title = plot_title,
         subtitle = "1:1 line (dashed), linear regression (solid), values ≤12",
         x = expression(observed ~ conductance ~ (m ~ s^{-1})),
         y = expression(lpj ~ simulated ~ conductance ~ (m ~ s^{-1}))) +
    base_theme +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 11, face = "bold"))
  
  # add statistics for each species in its own facet
  for(spp in species_levels) {
    ann_data <- annotation_data %>% filter(species == spp)
    if(nrow(ann_data) > 0) {
      p <- p + geom_text(
        data = ann_data,
        aes(x = -Inf, y = Inf, label = label),
        hjust = -0.05, vjust = 1.1,
        size = 3.2, color = "black",
        inherit.aes = FALSE
      )
    }
  }
  
  return(p)
}

# create all scatter plots
scatter_control_daily <- create_scatter_plot(scatter_data_daily, scatter_stats_daily, 
                                             "control", "scatter plot (daily mean): control treatment")
scatter_drought_daily <- create_scatter_plot(scatter_data_daily, scatter_stats_daily, 
                                             "drought", "scatter plot (daily mean): drought treatment")
scatter_control_top10 <- create_scatter_plot(scatter_data_top10, scatter_stats_top10, 
                                             "control", "scatter plot (top 10% quantile): control treatment")
scatter_drought_top10 <- create_scatter_plot(scatter_data_top10, scatter_stats_top10, 
                                             "drought", "scatter plot (top 10% quantile): drought treatment")

# ==========================================================================
# 9. prepare data for difference calculations
# ==========================================================================

daily_combined <- observed_daily_mean %>%
  inner_join(lpj_processed_filtered, by = c("date", "species", "treatment")) %>%
  rename(G_sim = G_cond_mean, G_obs = G_obs_mean) %>%
  mutate(diff = G_sim - G_obs, method = "daily mean")

top10_combined <- observed_top10 %>%
  inner_join(lpj_processed_filtered, by = c("date", "species", "treatment")) %>%
  rename(G_sim = G_cond_mean, G_obs = G_obs_top10) %>%
  mutate(diff = G_sim - G_obs, method = "top 10%")

# ==========================================================================
# 10. MEAN DIFFERENCE BAR PLOTS
# ==========================================================================

diff_daily_summary <- daily_combined %>%
  group_by(species, treatment) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff = sd(diff, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(method = "daily mean")

diff_top10_summary <- top10_combined %>%
  group_by(species, treatment) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff = sd(diff, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(method = "top 10%")

diff_summary <- bind_rows(diff_daily_summary, diff_top10_summary)

y_max <- max(abs(diff_summary$mean_diff) + diff_summary$se_diff) * 1.1
y_min <- -y_max

plot_diff_summary <- ggplot(diff_summary, aes(x = species, y = mean_diff, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = mean_diff - se_diff, ymax = mean_diff + se_diff),
                width = 0.25, position = position_dodge(0.9)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  facet_wrap(~ method, ncol = 2) +
  scale_fill_manual(name = "treatment", values = c("control" = "#1f77b4", "drought" = "#d62728")) +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(title = "mean difference: simulated - observed conductance",
       subtitle = "positive = overestimation, negative = underestimation",
       x = "species", y = expression(mean ~ difference ~ (m ~ s^{-1}))) +
  base_theme + theme(axis.text.x = element_text(angle = 0, size = 11))

plot_diff_by_species <- ggplot(diff_summary, aes(x = method, y = mean_diff, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = mean_diff - se_diff, ymax = mean_diff + se_diff),
                width = 0.25, position = position_dodge(0.9)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  facet_wrap(~ species, ncol = 2) +
  scale_fill_manual(name = "treatment", values = c("control" = "#1f77b4", "drought" = "#d62728")) +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(title = "mean difference: simulated - observed conductance",
       subtitle = "positive = overestimation, negative = underestimation",
       x = "method", y = expression(mean ~ difference ~ (m ~ s^{-1}))) +
  base_theme + theme(axis.text.x = element_text(angle = 0, size = 11))

# ==========================================================================
# 11. DAILY DIFFERENCE BAR PLOTS
# ==========================================================================

plot_daily_diff_daily <- ggplot(daily_combined, aes(x = date, y = diff, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  facet_wrap(~ species, ncol = 2, scales = "free_y") +
  scale_fill_manual(name = "treatment", values = c("control" = "#1f77b4", "drought" = "#d62728")) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(title = "daily difference (simulated - observed): daily mean method",
       subtitle = "positive = overestimation, negative = underestimation",
       x = "timeline", y = expression(difference ~ (m ~ s^{-1}))) +
  base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_daily_diff_top10 <- ggplot(top10_combined, aes(x = date, y = diff, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  facet_wrap(~ species, ncol = 2, scales = "free_y") +
  scale_fill_manual(name = "treatment", values = c("control" = "#1f77b4", "drought" = "#d62728")) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(title = "daily difference (simulated - observed): top 10% quantile method",
       subtitle = "positive = overestimation, negative = underestimation",
       x = "timeline", y = expression(difference ~ (m ~ s^{-1}))) +
  base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ==========================================================================
# 12. BOXPLOTS OF DIFFERENCES
# ==========================================================================

plot_boxplot_daily <- ggplot(daily_combined, aes(x = species, y = diff, fill = treatment)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(name = "treatment", values = c("control" = "#1f77b4", "drought" = "#d62728")) +
  labs(title = "distribution of daily differences: daily mean method",
       subtitle = "boxplot shows median, quartiles, and outliers",
       x = "species", y = expression(difference ~ (m ~ s^{-1}))) +
  base_theme + theme(axis.text.x = element_text(angle = 0, size = 11))

plot_boxplot_top10 <- ggplot(top10_combined, aes(x = species, y = diff, fill = treatment)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(name = "treatment", values = c("control" = "#1f77b4", "drought" = "#d62728")) +
  labs(title = "distribution of daily differences: top 10% quantile method",
       subtitle = "boxplot shows median, quartiles, and outliers",
       x = "species", y = expression(difference ~ (m ~ s^{-1}))) +
  base_theme + theme(axis.text.x = element_text(angle = 0, size = 11))

# ==========================================================================
# 13. DISPLAY ALL PLOTS
# ==========================================================================

# time series
print(plot_control_ts)
print(plot_drought_ts)

# scatter plots
print(scatter_control_daily)
print(scatter_drought_daily)
print(scatter_control_top10)
print(scatter_drought_top10)

# difference plots
print(plot_diff_summary)
print(plot_diff_by_species)
print(plot_daily_diff_daily)
print(plot_daily_diff_top10)
print(plot_boxplot_daily)
print(plot_boxplot_top10)

# ==========================================================================
# 14. SAVE ALL PLOTS
# ==========================================================================

# time series
ggsave("Figures/compare_Gc_lpjtwd_hoelstein/timeseries_control.png", plot_control_ts, width = 12, height = 10, dpi = 300)
ggsave("Figures/compare_Gc_lpjtwd_hoelstein/timeseries_drought.png", plot_drought_ts, width = 12, height = 10, dpi = 300)

# scatter plots
ggsave("Figures/compare_Gc_lpjtwd_hoelstein/scatter_control_daily_mean.png", scatter_control_daily, width = 11, height = 9, dpi = 300)
ggsave("Figures/compare_Gc_lpjtwd_hoelstein/scatter_drought_daily_mean.png", scatter_drought_daily, width = 11, height = 9, dpi = 300)
ggsave("Figures/compare_Gc_lpjtwd_hoelstein/scatter_control_top10.png", scatter_control_top10, width = 11, height = 9, dpi = 300)
ggsave("Figures/compare_Gc_lpjtwd_hoelstein/scatter_drought_top10.png", scatter_drought_top10, width = 11, height = 9, dpi = 300)

# difference plots
ggsave("Figures/compare_Gc_lpjtwd_hoelstein/mean_difference_by_method.png", plot_diff_summary, width = 10, height = 8, dpi = 300)
ggsave("Figures/compare_Gc_lpjtwd_hoelstein/mean_difference_by_species.png", plot_diff_by_species, width = 10, height = 8, dpi = 300)
ggsave("Figures/compare_Gc_lpjtwd_hoelstein/daily_difference_daily_mean.png", plot_daily_diff_daily, width = 14, height = 10, dpi = 300)
ggsave("Figures/compare_Gc_lpjtwd_hoelstein/daily_difference_top10.png", plot_daily_diff_top10, width = 14, height = 10, dpi = 300)
ggsave("Figures/compare_Gc_lpjtwd_hoelstein/boxplot_differences_daily_mean.png", plot_boxplot_daily, width = 10, height = 8, dpi = 300)
ggsave("Figures/compare_Gc_lpjtwd_hoelstein/boxplot_differences_top10.png", plot_boxplot_top10, width = 10, height = 8, dpi = 300)

# ==========================================================================
# 15. SAVE STATISTICS
# ==========================================================================

write.csv(scatter_stats_daily, "Figures/compare_Gc_lpjtwd_hoelstein/scatter_statistics_daily_mean.csv", row.names = FALSE)
write.csv(scatter_stats_top10, "Figures/compare_Gc_lpjtwd_hoelstein/scatter_statistics_top10.csv", row.names = FALSE)
write.csv(diff_summary, "Figures/compare_Gc_lpjtwd_hoelstein/difference_statistics.csv", row.names = FALSE)
write.csv(daily_combined, "Figures/compare_Gc_lpjtwd_hoelstein/daily_differences_daily_mean.csv", row.names = FALSE)
write.csv(top10_combined, "Figures/compare_Gc_lpjtwd_hoelstein/daily_differences_top10.csv", row.names = FALSE)

# ==========================================================================
# 16. STATISTICAL VARIABLES EXPLANATION FILE
# ==========================================================================

create_stats_explanation <- function(output_file = "Figures/compare_Gc_lpjtwd_hoelstein/statistical_variables_explanation.txt") {
  sink(output_file)
  cat("================================================================================\n")
  cat("STATISTICAL VARIABLES EXPLANATION\n")
  cat("Canopy Conductance Comparison: Observed vs LPJ-GUESS Simulated\n")
  cat("================================================================================\n\n")
  cat("Date created:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  cat("n (Number of observations)\n")
  cat("  - Definition: Total number of paired data points\n")
  cat("  - Ideal: As large as possible\n\n")
  
  cat("pearson_r (Pearson correlation coefficient)\n")
  cat("  - Definition: Measures linear relationship between observed and simulated\n")
  cat("  - Range: -1 to +1, Ideal: +1\n\n")
  
  cat("pearson_r2 (Coefficient of determination)\n")
  cat("  - Definition: Proportion of variance explained by the model\n")
  cat("  - Range: 0 to 1, Ideal: 1\n\n")
  
  cat("rmse (Root Mean Square Error)\n")
  cat("  - Definition: Standard deviation of prediction errors\n")
  cat("  - Units: m s⁻¹, Ideal: 0\n\n")
  
  cat("nrmse (Normalized RMSE)\n")
  cat("  - Definition: RMSE as percentage of observed mean\n")
  cat("  - Units: %, Ideal: 0%, Acceptable: <30%\n\n")
  
  cat("bias (Mean Bias)\n")
  cat("  - Definition: Average difference (simulated - observed)\n")
  cat("  - Units: m s⁻¹, Ideal: 0\n\n")
  
  cat("slope (Regression Slope)\n")
  cat("  - Definition: Slope of regression line\n")
  cat("  - Ideal: 1 (perfect 1:1 relationship)\n\n")
  
  cat("mean_diff (Mean Difference)\n")
  cat("  - Definition: Average (simulated - observed) across all time points\n")
  cat("  - Units: m s⁻¹, Ideal: 0\n\n")
  
  cat("diff (Daily Difference)\n")
  cat("  - Definition: Daily (simulated - observed)\n")
  cat("  - Units: m s⁻¹, Positive = overestimation\n\n")
  
  cat("================================================================================\n")
  cat("END OF DOCUMENT\n")
  cat("================================================================================\n")
  sink()
  cat("✓ Statistical variables explanation saved to:", output_file, "\n")
}

create_stats_explanation()

# ==========================================================================
# 17. SUMMARY TABLE
# ==========================================================================

create_stats_table <- function(output_csv = "Figures/compare_Gc_lpjtwd_hoelstein/statistical_variables_summary.csv") {
  stats_summary <- data.frame(
    Variable = c("n", "pearson_r", "pearson_r2", "rmse", "nrmse", "bias", "slope", "mean_diff", "diff"),
    Full_Name = c("Number of observations", "Pearson correlation", "R-squared", 
                  "Root mean square error", "Normalized RMSE", "Mean bias", 
                  "Regression slope", "Mean difference", "Daily difference"),
    Units = c("count", "dimensionless", "dimensionless", "m s⁻¹", "%", "m s⁻¹", "dimensionless", "m s⁻¹", "m s⁻¹"),
    Ideal_Value = c("As high as possible", "1.0", "1.0", "0", "0%", "0", "1.0", "0", "0")
  )
  write.csv(stats_summary, output_csv, row.names = FALSE)
  cat("✓ Statistical variables summary saved to:", output_csv, "\n")
}

create_stats_table()

# ==========================================================================
# 18. FINAL SUMMARY
# ==========================================================================

cat("\n================================================================================\n")
cat("ALL FIGURES COMPLETED\n")
cat("================================================================================\n")
cat("\n✓ Figures saved to: Figures/compare_Gc_lpjtwd_hoelstein/\n")
cat("\nFiles created:\n")
cat("  TIME SERIES PLOTS:\n")
cat("    - timeseries_control.png\n")
cat("    - timeseries_drought.png\n")
cat("  SCATTER PLOTS:\n")
cat("    - scatter_control_daily_mean.png\n")
cat("    - scatter_drought_daily_mean.png\n")
cat("    - scatter_control_top10.png\n")
cat("    - scatter_drought_top10.png\n")
cat("  DIFFERENCE PLOTS:\n")
cat("    - mean_difference_by_method.png\n")
cat("    - mean_difference_by_species.png\n")
cat("    - daily_difference_daily_mean.png\n")
cat("    - daily_difference_top10.png\n")
cat("    - boxplot_differences_daily_mean.png\n")
cat("    - boxplot_differences_top10.png\n")
cat("  STATISTICS FILES:\n")
cat("    - scatter_statistics_daily_mean.csv\n")
cat("    - scatter_statistics_top10.csv\n")
cat("    - difference_statistics.csv\n")
cat("    - daily_differences_daily_mean.csv\n")
cat("    - daily_differences_top10.csv\n")
cat("    - statistical_variables_explanation.txt\n")
cat("    - statistical_variables_summary.csv\n")
cat("\n================================================================================\n")

cat("\n=== DAILY MEAN STATISTICS (per species per treatment) ===\n")
print(scatter_stats_daily)

cat("\n=== TOP 10% STATISTICS (per species per treatment) ===\n")
print(scatter_stats_top10)
