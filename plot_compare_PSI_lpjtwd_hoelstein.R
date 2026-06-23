# ==========================================================================
# 1. setup, theme, & paths
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

dir.create("Figures/lpj_guess_hyd_twd", recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# 2. raw data ingestion & standardization (all unified to mpa)
# ==========================================================================

# LPJ Model Output (Already in MPa)
lpj_raw <- read.csv("lpj_guess/lpj_guess_twd/lpj_control_drought_Gc_psiL_psiX_psiS_stem_diameter_twd_stem_rwc.csv") %>%
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
# 3. intersection matrix (dates where ALL water potentials exist)
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

# ==========================================
# 4. unified time-series plot
# ==========================================
# Define the treatments you want to plot
treatments <- unique(lpj_intersect$treatment)

for (t in treatments) {
  
  # Filter both data frames for the current treatment
  lpj_sub <- lpj_intersect[lpj_intersect$treatment == t, ]
  obs_sub <- obs_intersect[obs_intersect$treatment == t, ]
  
  # Generate the plot
  plot_t <- ggplot() +
    geom_line(data = lpj_sub, aes(x = date, y = psi_soil_model, color = species), linewidth = 0.4, alpha = 0.6) +
    geom_line(data = lpj_sub, aes(x = date, y = psi_xylem_model, color = species), linetype = "dotdash", linewidth = 0.6) +
    geom_line(data = lpj_sub, aes(x = date, y = psi_leaf_model, color = species), linewidth = 1.1) +
    
    geom_line(data = obs_sub, aes(x = date, y = psi_leaf_md_obs), color = "grey", linewidth = 0.5) +
    geom_point(data = obs_sub, aes(x = date, y = psi_leaf_md_obs), color = "grey", shape = 17, size = 1.8) +
    
    geom_line(data = obs_sub, aes(x = date, y = psi_leaf_pd_obs), color = "grey40", linewidth = 0.5, linetype = "dashed") +
    geom_point(data = obs_sub, aes(x = date, y = psi_leaf_pd_obs), color = "grey40", shape = 15, size = 1.6) +
    
    geom_line(data = obs_sub, aes(x = date, y = psi_soil_obs), color = "black", linetype = "dotted", linewidth = 0.8) +
    geom_point(data = obs_sub, aes(x = date, y = psi_soil_obs), color = "black", shape = 18, size = 2.2) +
    
    # Update facet: remove 'treatment' from formula
    facet_grid(. ~ species, scales = "free_y") + 
    scale_color_manual(values = species_colors) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    ylim(-3.5, 0) +
    labs(
      title = paste("Water potential intersection:", t),
      subtitle = "Model: thick solid=\u03c8l, dotdash=\u03c8x, thin=\u03c8s | Obs: \u25b2=\u03c8l md, \u25a0=\u03c8l pd, \u25c6=\u03c8s",
      x = "year", y = expression(psi ~ (mpa)), color = "Species"
    ) +
    base_theme
  
  # Save with a dynamic filename
  ggsave(paste0("Figures/lpj_guess_hyd_twd/water_potential_", t, ".png"), 
         plot_t, width = 16, height = 5, dpi = 300) # Adjusted height since we removed a facet row
}

# ==========================================
# 5. data pair cohort restructuring (flat matrix for all 5 evaluation pairs)
# ==========================================
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

# ==========================================
# 6. calculate scatter statistics per species, treatment, and pair
# ==========================================
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

write.csv(scatter_stats_psi, "Figures/lpj_guess_hyd_twd/scatter_statistics_psi_5pairs.csv", row.names = FALSE)

# ==========================================
# 7. generate annotated 1:1 scatter plots 
# ==========================================
create_scatter_psi_plot <- function(target_treatment, plot_title) {
  data_subset <- flat_evaluation_pairs %>% filter(treatment == target_treatment)
  stats_subset <- scatter_stats_psi %>% filter(treatment == target_treatment)
  
  # Map human-readable text labels for the facets
  pair_labels <- c(
    "psiL_vs_md"   = "model \u03c8l vs obs midday",
    "psiL_vs_pd"   = "model \u03c8l vs obs predawn",
    "psiS_vs_soil" = "model \u03c8s vs obs soil",
    "psiX_vs_md"   = "model \u03c8x vs obs midday",
    "psiX_vs_pd"   = "model \u03c8x vs obs predawn"
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
  
  # Establish clean isometric square dimensions
  axis_min <- min(c(data_subset$obs_val, data_subset$model_val), na.rm = TRUE) * 1.05
  
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
      subtitle = "dashed line = 1:1 identity framework | june-september common records",
      x = expression(observed ~ water ~ potential ~ (mpa)),
      y = expression(lpj-guess ~ simulated ~ water ~ potential ~ (mpa))
    ) +
    base_theme +
    theme(aspect.ratio = 1, strip.text.y = element_text(angle = -90, size = 9))
}

scatter_control <- create_scatter_psi_plot("control", "water potential 1:1 scatter: control")
scatter_drought <- create_scatter_psi_plot("drought", "water potential 1:1 scatter: drought")

ggsave("Figures/lpj_guess_hyd_twd/scatter_control_5pairs.png", scatter_control, width = 14, height = 15, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/scatter_drought_5pairs.png", scatter_drought, width = 14, height = 15, dpi = 300)

# ==========================================
# 8. difference plots: bar summary (mean difference by pair)
# ==========================================
diff_summary_psi <- flat_evaluation_pairs %>%
  group_by(species, treatment, pair_id) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    se_diff = sd(diff, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

write.csv(diff_summary_psi, "Figures/lpj_guess_hyd_twd/difference_statistics_psi.csv", row.names = FALSE)

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
    name = "evaluation pairs", 
    values = pair_colors,
    labels = c("model \u03c8l vs obs md", "model \u03c8l vs obs pd", "model \u03c8s vs obs soil", "model \u03c8x vs obs md", "model \u03c8x vs obs pd")
  ) +
  labs(
    title = "mean difference summary: simulated - observed water potentials",
    subtitle = "positive = overestimation, negative = underestimation | error bars represent \u00b11 se",
    x = "species architecture", y = expression(mean ~ difference ~ (mpa))
  ) +
  base_theme + 
  theme(axis.text.x = element_text(angle = 0, size = 11))

ggsave("Figures/lpj_guess_hyd_twd/mean_difference_by_pair.png", plot_diff_summary_psi, width = 12, height = 7, dpi = 300)

# ==========================================
# 9. difference plots: timeline seriation (daily line fluctuations)
# ==========================================
create_daily_diff_timeline <- function(target_treatment, plot_title) {
  data_subset <- flat_evaluation_pairs %>% filter(treatment == target_treatment)
  
  ggplot(data_subset, aes(x = date, y = diff, color = pair_id, group = pair_id)) +
    geom_line(linewidth = 0.5, alpha = 0.75) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    facet_wrap(~ species, ncol = 2, scales = "free_y") +
    scale_color_manual(
      name = "evaluation pairs", 
      values = pair_colors,
      labels = c("model \u03c8l vs obs md", "model \u03c8l vs obs pd", "model \u03c8s vs obs soil", "model \u03c8x vs obs md", "model \u03c8x vs obs pd")
    ) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
    labs(
      title = plot_title,
      subtitle = "daily difference line metric (simulated - observed) | positive = overestimation, negative = underestimation",
      x = "timeline", y = expression(daily ~ delta ~ psi ~ (mpa))
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}

timeline_control_diff <- create_daily_diff_timeline("control", "daily evaluation differences: control experiment (line)")
timeline_drought_diff <- create_daily_diff_timeline("drought", "daily evaluation differences: drought experiment (line)")

ggsave("Figures/lpj_guess_hyd_twd/daily_difference_control.png", timeline_control_diff, width = 14, height = 10, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/daily_difference_drought.png", timeline_drought_diff, width = 14, height = 10, dpi = 300)

# ==========================================
# 10. difference plots: boxplots (separated by treatment)
# ==========================================
create_boxplot_psi_plot <- function(target_treatment, plot_title) {
  data_subset <- flat_evaluation_pairs %>% filter(treatment == target_treatment)
  
  ggplot(data_subset, aes(x = species, y = diff, fill = pair_id)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8, position = position_dodge(0.85)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    scale_fill_manual(
      name = "evaluation pairs", 
      values = pair_colors,
      labels = c("model \u03c8l vs obs md", "model \u03c8l vs obs pd", "model \u03c8s vs obs soil", "model \u03c8x vs obs md", "model \u03c8x vs obs pd")
    ) +
    labs(
      title = plot_title,
      subtitle = "boxplots display medians, interquartile ranges, and associated outliers",
      x = "species architecture", y = expression(delta ~ psi ~ (mpa))
    ) +
    base_theme + 
    theme(axis.text.x = element_text(angle = 0, size = 11))
}

boxplot_control <- create_boxplot_psi_plot("control", "distribution of daily differences: control experiment")
boxplot_drought <- create_boxplot_psi_plot("drought", "distribution of daily differences: drought experiment")

ggsave("Figures/lpj_guess_hyd_twd/boxplot_differences_psi_control.png", boxplot_control, width = 10, height = 7, dpi = 300)
ggsave("Figures/lpj_guess_hyd_twd/boxplot_differences_psi_drought.png", boxplot_drought, width = 10, height = 7, dpi = 300)

# ==========================================
# 11. metadata statistical dictionary file generation
# ==========================================
cat("\n=== logging metadata explanations ===\n")
sink("Figures/lpj_guess_hyd_twd/statistical_variables_explanation.txt")
cat("================================================================================\n")
cat("statistical variables explanation dictionary\n")
cat("water potential evaluation matrix: field observations vs lpj-guess simulated\n")
cat("================================================================================\n\n")
cat("n (sample cohort size): paired strict chronological date match records.\n")
cat("pearson_r (pearson correlation): direction and linear strength (-1 to +1).\n")
cat("pearson_r2 (coefficient of determination): explained variation proportion (0 to 1).\n")
cat("rmse (root mean square error): absolute quadratic fit penalty metric (mpa).\n")
cat("bias (mean directional bias): mean delta calculation (simulated - observed) (mpa).\n")
cat("slope (fitted linear slope): regression slope. ideal matching target = 1.0.\n")
sink()
cat("\u2713 statistical dictionary written successfully.\n")