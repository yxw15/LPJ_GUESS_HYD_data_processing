# ==========================================================================
# 1. SETUP, THEME, & PATHS
# ==========================================================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(scales)

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

species_levels <- c("Oak", "Beech", "Spruce", "Pine")
species_colors <- c(Oak = "#E69F00", Beech = "#0072B2", Spruce = "#009E73", Pine = "#F0E442")

base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text       = element_text(color = "black", size = 11),
    legend.position   = "bottom",
    legend.box        = "vertical",
    plot.title        = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5, size = 10, color = "grey30"),
    axis.title        = element_text(size = 12),
    axis.text.x       = element_text(angle = 0, hjust = 0.5, size = 9),
    strip.text        = element_text(size = 11, face = "bold"),
    panel.grid.major  = element_line(color = "grey92"),
    panel.grid.minor  = element_blank()
  )

dir.create("Figures/compare_Psi_lpjtwd_hoelstein", recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# 2. RAW DATA INGESTION & STANDARDIZATION (All Unified to MPa)
# ==========================================================================

# LPJ Model Output (Already in MPa)
lpj_raw <- read.csv("lpj_guess/lpj_control_drought_Gc_psiL_psiX_psiS_stem_diameter_twd_stem_rwc.csv") %>%
  mutate(date = as.Date(date), month = month(date)) %>%
  filter(month >= 6 & month <= 9) %>%
  filter(!is.na(psi_leaf), !is.na(psi_xylem), !is.na(psi_soil), !is.na(species), !is.na(treatment)) %>%
  mutate(species = factor(species, levels = species_levels)) %>%
  select(date, species, treatment, psi_soil_model = psi_soil, psi_xylem_model = psi_xylem, psi_leaf_model = psi_leaf)

# Observed Soil Data (Converted from kPa to MPa by dividing by 1000)
obs_soil <- bind_rows(
  read.csv("SCCII/psiS_hoelstein_drought.csv") %>% mutate(treatment = "drought"),
  read.csv("SCCII/psiS_hoelstein_control.csv") %>% mutate(treatment = "control")
) %>%
  mutate(date = as.Date(date), month = month(date)) %>%
  filter(month >= 6 & month <= 9 & !is.na(psiS_mean)) %>%
  mutate(psi_soil_obs = psiS_mean / 1000) %>% 
  select(date, treatment, psi_soil_obs)

# Observed Leaf Data (Already in MPa)
obs_leaf_raw <- bind_rows(
  read.csv("SCCII/psiL_hoelstein_drought.csv") %>% mutate(treatment = "drought"),
  read.csv("SCCII/psiL_hoelstein_control.csv") %>% mutate(treatment = "control")
) %>%
  mutate(date = as.Date(date), month = month(date)) %>%
  filter(month >= 6 & month <= 9) %>%
  rename(species = species_name) %>%
  mutate(species = factor(species, levels = species_levels))

# Extract clean daily aggregates for Midday and Predawn
obs_leaf_processed <- obs_leaf_raw %>%
  group_by(date, species, treatment) %>%
  summarise(
    psi_leaf_md_obs = mean(md_wp_av, na.rm = TRUE),
    psi_leaf_pd_obs = mean(pd_wp_av, na.rm = TRUE),
    .groups = "drop"
  )

# ==========================================================================
# 3. INTERSECTION MATRIX (Dates where ALL water potentials exist)
# ==========================================================================

common_dates_full <- lpj_raw %>%
  inner_join(obs_leaf_processed %>% filter(!is.na(psi_leaf_md_obs) & !is.na(psi_leaf_pd_obs)), 
             by = c("date", "species", "treatment")) %>%
  inner_join(obs_soil, by = c("date", "treatment")) %>%
  select(date, species, treatment) %>%
  distinct()

lpj_intersect <- lpj_raw %>% inner_join(common_dates_full, by = c("date", "species", "treatment"))
obs_intersect <- obs_leaf_processed %>% 
  inner_join(common_dates_full, by = c("date", "species", "treatment")) %>%
  left_join(obs_soil, by = c("date", "treatment")) 

# ==========================================================================
# 4. UNIFIED TIME-SERIES PLOT
# ==========================================================================
plot_common_intersection <- ggplot() +
  geom_line(data = lpj_intersect, aes(x = date, y = psi_soil_model, color = species), linewidth = 0.4, alpha = 0.6) +
  geom_line(data = lpj_intersect, aes(x = date, y = psi_xylem_model, color = species), linetype = "dotdash", linewidth = 0.6) +
  geom_line(data = lpj_intersect, aes(x = date, y = psi_leaf_model, color = species), linewidth = 1.1) +
  
  geom_line(data = obs_intersect, aes(x = date, y = psi_leaf_md_obs), color = "red", linewidth = 0.5) +
  geom_point(data = obs_intersect, aes(x = date, y = psi_leaf_md_obs), color = "red", shape = 17, size = 1.8) +
  
  geom_line(data = obs_intersect, aes(x = date, y = psi_leaf_pd_obs), color = "blue", linewidth = 0.5, linetype = "dashed") +
  geom_point(data = obs_intersect, aes(x = date, y = psi_leaf_pd_obs), color = "blue", shape = 15, size = 1.6) +
  
  geom_line(data = obs_intersect, aes(x = date, y = psi_soil_obs), color = "black", linetype = "dotted", linewidth = 0.8) +
  geom_point(data = obs_intersect, aes(x = date, y = psi_soil_obs), color = "black", shape = 18, size = 2.2) +
  
  facet_grid(treatment ~ species, scales = "free_y") +
  scale_color_manual(values = species_colors) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  ylim(-3.5, 0) +
  labs(
    title = "Water Potential Intersection Comparison (June-September)",
    subtitle = "Model: Thick Solid=ψL, Dotdash=ψX, Thin=ψS | Obs: ▲=ψL Midday (red), ■=ψL Predawn (blue), ◆=ψS Soil (black dotted)",
    x = "Year", y = expression(Psi ~ (MPa)), color = "Species Framework"
  ) +
  base_theme

ggsave("Figures/compare_Psi_lpjtwd_hoelstein/water_potential_common.png", 
       plot_common_intersection, width = 16, height = 10, dpi = 300)

# ==========================================================================
# 5. DATA PAIR COHORT RESTRUCTURING (Flat Matrix for All 5 Evaluation Pairs)
# ==========================================================================
combined_matrix <- lpj_intersect %>%
  inner_join(obs_intersect, by = c("date", "species", "treatment"))

pair_levels <- c("psiL_vs_md", "psiL_vs_pd", "psiS_vs_soil", "psiX_vs_md", "psiX_vs_pd")

flat_evaluation_pairs <- bind_rows(
  combined_matrix %>% transmute(date, species, treatment, pair_id = "psiL_vs_md",   model_val = psi_leaf_model,  obs_val = psi_leaf_md_obs),
  combined_matrix %>% transmute(date, species, treatment, pair_id = "psiL_vs_pd",   model_val = psi_leaf_model,  obs_val = psi_leaf_pd_obs),
  combined_matrix %>% transmute(date, species, treatment, pair_id = "psiS_vs_soil", model_val = psi_soil_model,  obs_val = psi_soil_obs),
  combined_matrix %>% transmute(date, species, treatment, pair_id = "psiX_vs_md",   model_val = psi_xylem_model, obs_val = psi_leaf_md_obs),
  combined_matrix %>% transmute(date, species, treatment, pair_id = "psiX_vs_pd",   model_val = psi_xylem_model, obs_val = psi_leaf_pd_obs)
) %>%
  filter(!is.na(model_val) & !is.na(obs_val)) %>%
  mutate(
    species = factor(species, levels = species_levels),
    pair_id = factor(pair_id, levels = pair_levels),
    diff = model_val - obs_val
  )

# ==========================================================================
# 6. CALCULATE SCATTER STATISTICS PER SPECIES, TREATMENT, AND PAIR
# ==========================================================================
scatter_stats_psi <- flat_evaluation_pairs %>%
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

write.csv(scatter_stats_psi, "Figures/compare_Psi_lpjtwd_hoelstein/scatter_statistics_psi_5pairs.csv", row.names = FALSE)

# ==========================================================================
# 7. GENERATE ANNOTATED 1:1 SCATTER PLOTS 
# ==========================================================================
create_scatter_psi_plot <- function(target_treatment, plot_title) {
  data_subset <- flat_evaluation_pairs %>% filter(treatment == target_treatment)
  stats_subset <- scatter_stats_psi %>% filter(treatment == target_treatment)
  
  # Map human-readable text labels for the facets
  pair_labels <- c(
    "psiL_vs_md"   = "Model ψL vs Obs Midday",
    "psiL_vs_pd"   = "Model ψL vs Obs Predawn",
    "psiS_vs_soil" = "Model ψS vs Obs Soil",
    "psiX_vs_md"   = "Model ψX vs Obs Midday",
    "psiX_vs_pd"   = "Model ψX vs Obs Predawn"
  )
  
  data_subset <- data_subset %>% mutate(pair_label = pair_labels[pair_id])
  
  annotation_data <- stats_subset %>%
    mutate(
      pair_label = pair_labels[pair_id],
      text_summary = paste0(
        "n = ", n, "\n",
        "r = ", round(pearson_r, 2), "\n",
        "r² = ", round(pearson_r2, 2), "\n",
        "rmse = ", round(rmse, 2), "\n",
        "bias = ", round(bias, 2), "\n",
        "slope = ", round(slope, 1)
      )
    )
  
  # Establish clean isometric square dimensions
  axis_min <- min(c(data_subset$obs_val, data_subset$model_val), na.rm = TRUE) * 1.05
  axis_max <- max(c(data_subset$obs_val, data_subset$model_val), na.rm = TRUE) * 0.95
  
  ggplot(data_subset, aes(x = obs_val, y = model_val, color = species)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.1, linewidth = 0.6) +
    facet_grid(pair_label ~ species) +
    scale_color_manual(values = species_colors) +
    coord_fixed(ratio = 1, xlim = c(axis_min, 0), ylim = c(axis_min, 0)) +
    geom_text(
      data = annotation_data,
      aes(x = -Inf, y = Inf, label = text_summary),
      hjust = -0.05, vjust = 1.1, size = 2.8, color = "black", inherit.aes = FALSE
    ) +
    labs(
      title = plot_title,
      subtitle = "Dashed line = 1:1 Identity framework | June-September common records",
      x = expression(Observed ~ Water ~ Potential ~ (MPa)),
      y = expression(LPJ-GUESS ~ Simulated ~ Water ~ Potential ~ (MPa))
    ) +
    base_theme +
    theme(aspect.ratio = 1, strip.text.y = element_text(angle = -90, size = 9))
}

scatter_control <- create_scatter_psi_plot("control", "Water Potential 1:1 Scatter: Control Cohort")
scatter_drought <- create_scatter_psi_plot("drought", "Water Potential 1:1 Scatter: Drought Cohort")

ggsave("Figures/compare_Psi_lpjtwd_hoelstein/scatter_control_5pairs.png", scatter_control, width = 14, height = 15, dpi = 300)
ggsave("Figures/compare_Psi_lpjtwd_hoelstein/scatter_drought_5pairs.png", scatter_drought, width = 14, height = 15, dpi = 300)

# ==========================================================================
# 8. DIFFERENCE PLOTS: BAR SUMMARY (MEAN DIFFERENCE BY PAIR)
# ==========================================================================
diff_summary_psi <- flat_evaluation_pairs %>%
  group_by(species, treatment, pair_id) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff = sd(diff, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

write.csv(diff_summary_psi, "Figures/compare_Psi_lpjtwd_hoelstein/difference_statistics_psi.csv", row.names = FALSE)

# Color scale to easily distinguish the 5 pairs (Reds for Midday, Blues for Predawn, Orange for Soil)
pair_colors <- c(
  "psiL_vs_md"   = "red2", 
  "psiL_vs_pd"   = "blue2", 
  "psiS_vs_soil" = "darkorange1", 
  "psiX_vs_md"   = "coral1", 
  "psiX_vs_pd"   = "deepskyblue"
)

plot_diff_summary_psi <- ggplot(diff_summary_psi, aes(x = species, y = mean_diff, fill = pair_id)) +
  geom_bar(stat = "identity", position = position_dodge(0.85), width = 0.8) +
  geom_errorbar(aes(ymin = mean_diff - se_diff, ymax = mean_diff + se_diff),
                width = 0.25, position = position_dodge(0.85)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  facet_wrap(~ treatment, ncol = 2) +
  scale_fill_manual(
    name = "Evaluation Pairs", 
    values = pair_colors,
    labels = c("Model ψL vs Obs MD", "Model ψL vs Obs PD", "Model ψS vs Obs Soil", "Model ψX vs Obs MD", "Model ψX vs Obs PD")
  ) +
  labs(
    title = "Mean Difference Summary: Simulated - Observed Water Potentials",
    subtitle = "Positive = Overestimation, Negative = Underestimation | Error bars represent ±1 SE",
    x = "Species Architecture", y = expression(Mean ~ Difference ~ (MPa))
  ) +
  base_theme + 
  theme(axis.text.x = element_text(angle = 0, size = 11))

ggsave("Figures/compare_Psi_lpjtwd_hoelstein/mean_difference_by_pair.png", plot_diff_summary_psi, width = 12, height = 7, dpi = 300)

# ==========================================================================
# 9. DIFFERENCE PLOTS: TIMELINE SERIATION (DAILY LINE FLUCTUATIONS)
# ==========================================================================
create_daily_diff_timeline <- function(target_treatment, plot_title) {
  data_subset <- flat_evaluation_pairs %>% filter(treatment == target_treatment)
  
  ggplot(data_subset, aes(x = date, y = diff, color = pair_id, group = pair_id)) +
    # Swapped geom_bar for lines and points to track daily trajectories clearly
    geom_line(linewidth = 0.5, alpha = 0.75) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    facet_wrap(~ species, ncol = 2, scales = "free_y") +
    scale_color_manual(
      name = "Evaluation Pairs", 
      values = pair_colors,
      labels = c("Model ψL vs Obs MD", "Model ψL vs Obs PD", "Model ψS vs Obs Soil", "Model ψX vs Obs MD", "Model ψX vs Obs PD")
    ) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
    labs(
      title = plot_title,
      subtitle = "Daily Difference Line Metric (Simulated - Observed) | Positive = Overestimation, Negative = Underestimation",
      x = "Timeline", y = expression(Daily ~ Delta ~ Psi ~ (MPa))
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}

timeline_control_diff <- create_daily_diff_timeline("control", "Daily Evaluation Differences: Control Experiment (Line)")
timeline_drought_diff <- create_daily_diff_timeline("drought", "Daily Evaluation Differences: Drought Experiment (Line)")

ggsave("Figures/compare_Psi_lpjtwd_hoelstein/daily_difference_control.png", timeline_control_diff, width = 14, height = 10, dpi = 300)
ggsave("Figures/compare_Psi_lpjtwd_hoelstein/daily_difference_drought.png", timeline_drought_diff, width = 14, height = 10, dpi = 300)


# ==========================================================================
# 10. DIFFERENCE PLOTS: BOXPLOTS (SEPARATED BY TREATMENT)
# ==========================================================================
create_boxplot_psi_plot <- function(target_treatment, plot_title) {
  data_subset <- flat_evaluation_pairs %>% filter(treatment == target_treatment)
  
  ggplot(data_subset, aes(x = species, y = diff, fill = pair_id)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8, position = position_dodge(0.85)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    scale_fill_manual(
      name = "Evaluation Pairs", 
      values = pair_colors,
      labels = c("Model ψL vs Obs MD", "Model ψL vs Obs PD", "Model ψS vs Obs Soil", "Model ψX vs Obs MD", "Model ψX vs Obs PD")
    ) +
    labs(
      title = plot_title,
      subtitle = "Boxplots display medians, interquartile ranges, and associated outliers",
      x = "Species Architecture", y = expression(Delta ~ Psi ~ (MPa))
    ) +
    base_theme + 
    theme(axis.text.x = element_text(angle = 0, size = 11))
}

boxplot_control <- create_boxplot_psi_plot("control", "Distribution of Daily Differences: Control Experiment")
boxplot_drought <- create_boxplot_psi_plot("drought", "Distribution of Daily Differences: Drought Experiment")

ggsave("Figures/compare_Psi_lpjtwd_hoelstein/boxplot_differences_psi_control.png", boxplot_control, width = 10, height = 7, dpi = 300)
ggsave("Figures/compare_Psi_lpjtwd_hoelstein/boxplot_differences_psi_drought.png", boxplot_drought, width = 10, height = 7, dpi = 300)

# ==========================================================================
# 11. METADATA STATISTICAL DICTIONARY FILE GENERATION
# ==========================================================================
cat("\n=== Logging Metadata Explanations ===\n")
sink("Figures/compare_Psi_lpjtwd_hoelstein/statistical_variables_explanation.txt")
cat("================================================================================\n")
cat("STATISTICAL VARIABLES EXPLANATION DICTIONARY\n")
cat("Water Potential Evaluation Matrix: Field Observations vs LPJ-GUESS Simulated\n")
cat("================================================================================\n\n")
cat("n (Sample Cohort Size): Paired strict chronological date match records.\n")
cat("pearson_r (Pearson Correlation): Direction and linear strength (-1 to +1).\n")
cat("pearson_r2 (Coefficient of Determination): Explained variation proportion (0 to 1).\n")
cat("rmse (Root Mean Square Error): Absolute quadratic fit penalty metric (MPa).\n")
cat("bias (Mean Directional Bias): Mean delta calculation (Simulated - Observed) (MPa).\n")
cat("slope (Fitted Linear Slope): Regression slope. Ideal matching target = 1.0.\n")
sink()
cat("✓ Statistical dictionary written successfully.\n")