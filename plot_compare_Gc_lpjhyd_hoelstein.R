# ==========================================================================
# 1. SETUP: GLOBAL CONSTANTS & THEMES
# ==========================================================================
library(dplyr)
library(ggplot2)
library(lubridate)
library(scales)
library(tidyr)
library(patchwork)

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

# aesthetics (lowercased factor levels to ensure lowercase strip texts)
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
# 2. DATA PREPARATION & FILTER DEFINITION
# ==========================================================================
lpj_output_filter <- read.csv("lpj_guess/lpj_guess_stem_storage/lpj_control_drought_ET_plant_ET_total_Gc_psi_leaf_psi_soil_psi_xylem_hydraulic_lag_kappy_s_min_mort_mort_cav_mort_greff_mort_min_stem_diameter_stem_rwc_twd_climate_filter.csv")
sap_flux_gc_filter <- read.csv("SCCII/sap_daily_filter.csv")

sap_flux_gc_filter <- sap_flux_gc_filter %>%
  mutate(treatment = case_when(
    treatment == "control" ~ "control",
    treatment == "treatment" ~ "drought"
  ))

# ==========================================================================
# 3. DATA PROCESSING: FIND COMMON DATES (EXACT DATE MATCHING) AND EXCLUDE >12
# ==========================================================================

# process observed data, lowercase species names, and exclude G_ms > 12
observed_data <- sap_flux_gc_filter %>%
  mutate(date = as.Date(date),
         species = tolower(species)) %>%
  filter(!is.na(G_ms)) %>%
  filter(!is.na(species)) %>%
  filter(!is.na(treatment)) %>%
  filter(G_ms <= 12) %>%
  mutate(species = factor(species, levels = species_levels))

# process lpj data, lowercase species names, and exclude Gc > 12
lpj_processed <- lpj_output_filter %>%
  mutate(date = as.Date(date),
         species = tolower(species)) %>%
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
# 4. DATA PROCESSING: OBSERVED DATA (BOTH METHODS)
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
# 5. TIME SERIES PLOT FUNCTIONS (SCATTER AND LINE METHODS)
# ==========================================================================

# Method A: Scatter Point Function
create_timeseries_plot <- function(treatment_name, plot_title) {
  
  lpj_subset <- lpj_processed_filtered %>% 
    filter(treatment == treatment_name) %>%
    mutate(year = year(date))
  
  obs_daily_subset <- observed_daily_mean %>% 
    filter(treatment == treatment_name) %>%
    mutate(year = year(date))
  
  obs_top10_subset <- observed_top10 %>% 
    filter(treatment == treatment_name) %>%
    mutate(year = year(date))
  
  treatment_dates <- common_dates %>% filter(treatment == treatment_name)
  t_start <- min(treatment_dates$date)
  t_end <- max(treatment_dates$date)
  
  p <- ggplot() +
    # 1. LPJ-GUESS Simulated Data (Points colored natively by species)
    geom_point(data = lpj_subset,
               aes(x = date, y = G_cond_mean, color = species),
               alpha = 1, size = 1.2) +
    
    # 2. Observed Daily Mean (Blue open triangles, shape = 2)
    geom_point(data = obs_daily_subset,
               aes(x = date, y = G_obs_mean, shape = "observed daily mean"),
               color = "grey70", alpha = 0.7, size = 1.8, stroke = 1) +
    
    # 3. Observed Top 10% Mean (Red open circles, shape = 1)
    geom_point(data = obs_top10_subset,
               aes(x = date, y = G_obs_top10, shape = "observed top 10%"),
               color = "grey20", alpha = 0.7, size = 1.8, stroke = 1) +
    
    facet_grid(species ~ year, scales = "free") +
    
    scale_color_manual(name = "lpj-guess simulation", values = species_colors) +
    
    scale_shape_manual(name = "field measurements (sccii)",
                       values = c("observed daily mean" = 2, "observed top 10%" = 1),
                       labels = c("observed daily mean" = "daily mean (blue open triangle)", 
                                  "observed top 10%" = "top 10% quantile (red open circle)")) +
    
    scale_x_date(date_breaks = "3 months", date_labels = "%b") +
    coord_cartesian(ylim = c(0, 12)) +
    labs(title = tolower(plot_title),
         subtitle = paste0("common dates: ", as.character(t_start), " to ", as.character(t_end)),
         x = "timeline",
         y = "conductance (m s-1)") +
    base_theme +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 12, face = "bold"),
          legend.box = "vertical",
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10))
  
  return(p)
}

# Method B: Trend Line Function
create_timeseries_plot_line <- function(treatment_name, plot_title) {
  
  lpj_subset <- lpj_processed_filtered %>% 
    filter(treatment == treatment_name) %>%
    mutate(year = year(date))
  
  obs_daily_subset <- observed_daily_mean %>% 
    filter(treatment == treatment_name) %>%
    mutate(year = year(date))
  
  obs_top10_subset <- observed_top10 %>% 
    filter(treatment == treatment_name) %>%
    mutate(year = year(date))
  
  treatment_dates <- common_dates %>% filter(treatment == treatment_name)
  t_start <- min(treatment_dates$date)
  t_end <- max(treatment_dates$date)
  
  p <- ggplot() +
    # 1. Model Simulation line (solid line, colored by species)
    geom_line(data = lpj_subset,
              aes(x = date, y = G_cond_mean, color = species),
              linetype = "solid", linewidth = 1, alpha = 1) +
    
    # 2. Observed Daily Mean line (dashed line, matches grey70 color)
    geom_line(data = obs_daily_subset,
              aes(x = date, y = G_obs_mean),
              color = "grey70", linetype = "dashed", linewidth = 0.8, alpha = 1) +
    
    # 3. Observed Top 10% line (dotted line, matches dark grey20 color)
    # geom_line(data = obs_top10_subset,
    #           aes(x = date, y = G_obs_top10),
    #           color = "grey20", linetype = "dotted", linewidth = 0.8, alpha = 1) +
    
    facet_grid(species ~ year, scales = "free") +
    
    scale_color_manual(name = "lpj-guess simulation", values = species_colors) +
    
    scale_x_date(date_breaks = "3 months", date_labels = "%b") +
    coord_cartesian(ylim = c(0, 12)) +
    labs(title = tolower(plot_title),
         subtitle = paste0("common dates: ", as.character(t_start), " to ", as.character(t_end)),
         x = "timeline",
         y = "conductance (m s-1)") +
    base_theme +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 12, face = "bold"),
          legend.box = "vertical",
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10))
  
  return(p)
}

# generate all time series instances
plot_control_ts <- create_timeseries_plot("control", "canopy conductance: control treatment")
plot_drought_ts <- create_timeseries_plot("drought", "canopy conductance: drought treatment")

plot_control_ts_line <- create_timeseries_plot_line("control", "canopy conductance: control treatment")
plot_drought_ts_line <- create_timeseries_plot_line("drought", "canopy conductance: drought treatment")

# ==========================================================================
# 6. PREPARE DATA FOR SCATTER PLOTS
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
# 7. CALCULATE STATISTICS PER SPECIES AND PER TREATMENT
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
# 8. SCATTER PLOTS WITH CORRECT PER-SPECIES STATISTICS
# ==========================================================================

create_scatter_plot <- function(data, stats, treatment_name, plot_title, max_limit = 12) {
  
  data_subset <- data %>% filter(treatment == treatment_name)
  stats_subset <- stats %>% filter(treatment == treatment_name)
  
  annotation_data <- stats_subset %>%
    mutate(
      label = paste0(
        "n = ", n, "\n",
        "r = ", round(pearson_r, 3), "\n",
        "r\u00b2 = ", round(pearson_r2, 3), "\n",
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
    geom_text(
      data = annotation_data,
      aes(x = -Inf, y = Inf, label = label),
      hjust = -0.05, vjust = 1.1,
      size = 3.2, color = "black",
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = species_colors, name = "species") +
    coord_fixed(ratio = 1, xlim = c(0, max_limit), ylim = c(0, max_limit)) +
    labs(title = tolower(plot_title),
         subtitle = "1:1 line (dashed), linear regression (solid), values \u226412",
         x = "observed conductance (m s-1)",
         y = "lpj simulated conductance (m s-1)") +
    base_theme +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 11, face = "bold"))
  
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
# 9. PREPARE DATA FOR DIFFERENCE CALCULATIONS
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
       x = "species", y = "mean difference (m s-1)") +
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
       x = "method", y = "mean difference (m s-1)") +
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
       x = "timeline", y = "difference (m s-1)") +
  base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_daily_diff_top10 <- ggplot(top10_combined, aes(x = date, y = diff, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  facet_wrap(~ species, ncol = 2, scales = "free_y") +
  scale_fill_manual(name = "treatment", values = c("control" = "#1f77b4", "drought" = "#d62728")) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(title = "daily difference (simulated - observed): top 10% quantile method",
       subtitle = "positive = overestimation, negative = underestimation",
       x = "timeline", y = "difference (m s-1)") +
  base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ==========================================================================
# 12. BOXPLOTS OF DIFFERENCES
# ==========================================================================

plot_boxplot_daily <- ggplot(daily_combined, aes(x = species, y = diff, fill = treatment)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(name = "treatment", values = c("control" = "#1f77b4", "drought" = "#d62728")) +
  labs(title = "distribution of daily differences: daily mean method",
       subtitle = "",
       x = "species", y = "difference Gc (m s-1)") +
  base_theme + theme(axis.text.x = element_text(angle = 0, size = 11))

plot_boxplot_top10 <- ggplot(top10_combined, aes(x = species, y = diff, fill = treatment)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(name = "treatment", values = c("control" = "#1f77b4", "drought" = "#d62728")) +
  labs(title = "distribution of daily differences: top 10% quantile method",
       subtitle = "",
       x = "species", y = "difference Gc (m s-1)") +
  base_theme + theme(axis.text.x = element_text(angle = 0, size = 11))


# 1. Determine common y-axis limits to ensure consistency
# You can find the min/max across both datasets
all_diffs <- c(daily_combined$diff, top10_combined$diff)
y_limits <- c(min(all_diffs, na.rm = TRUE), max(all_diffs, na.rm = TRUE))

# 2. Update your plots to use these limits
# Applying limits to both ensures they are visually comparable
p1 <- plot_boxplot_daily + ylim(y_limits)
p2 <- plot_boxplot_top10 + ylim(y_limits)

# 3. Combine them side-by-side
final_plot <- p1 + p2 + 
  plot_layout(guides = 'collect', ncol = 2) & 
  theme(legend.position = 'bottom')

# 5. Save the result
ggsave("Figures/lpj_guess_hyd_twd/combined_boxplots_diff_Gc.png", 
       final_plot, width = 14, height = 7, dpi = 300)

# ==========================================================================
# 13. DISPLAY ALL PLOTS
# ==========================================================================

# time series
print(plot_control_ts)
print(plot_drought_ts)
print(plot_control_ts_line)
print(plot_drought_ts_line)

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

# time series (points scatter)
ggsave("Figures/lpj_guess_hyd_twd/timeseries_control.png", plot_control_ts, width = 16, height = 11, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/timeseries_drought.png", plot_drought_ts, width = 16, height = 11, dpi = 300)

# time series (lines tracking)
ggsave("Figures/lpj_guess_hyd_twd/timeseries_control_lines.png", plot_control_ts_line, width = 16, height = 11, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/timeseries_drought_lines.png", plot_drought_ts_line, width = 16, height = 11, dpi = 300)

# scatter plots
ggsave("Figures/lpj_guess_hyd_twd/scatter_control_daily_mean.png", scatter_control_daily, width = 11, height = 9, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/scatter_drought_daily_mean.png", scatter_drought_daily, width = 11, height = 9, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/scatter_control_top10.png", scatter_control_top10, width = 11, height = 9, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/scatter_drought_top10.png", scatter_drought_top10, width = 11, height = 9, dpi = 300)

# difference plots
ggsave("Figures/lpj_guess_hyd_twd/mean_difference_by_method.png", plot_diff_summary, width = 10, height = 8, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/mean_difference_by_species.png", plot_diff_by_species, width = 10, height = 8, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/daily_difference_daily_mean.png", plot_daily_diff_daily, width = 14, height = 10, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/daily_difference_top10.png", plot_daily_diff_top10, width = 14, height = 10, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/boxplot_differences_daily_mean.png", plot_boxplot_daily, width = 10, height = 8, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/boxplot_differences_top10.png", plot_boxplot_top10, width = 10, height = 8, dpi = 300)

# ==========================================================================
# 15. SAVE STATISTICS
# ==========================================================================

write.csv(scatter_stats_daily, "Figures/lpj_guess_hyd_twd/scatter_statistics_daily_mean.csv", row.names = FALSE)
write.csv(scatter_stats_top10, "Figures/lpj_guess_hyd_twd/scatter_statistics_top10.csv", row.names = FALSE)
write.csv(diff_summary, "Figures/lpj_guess_hyd_twd/difference_statistics.csv", row.names = FALSE)
write.csv(daily_combined, "Figures/lpj_guess_hyd_twd/daily_differences_daily_mean.csv", row.names = FALSE)
write.csv(top10_combined, "Figures/lpj_guess_hyd_twd/daily_differences_top10.csv", row.names = FALSE)