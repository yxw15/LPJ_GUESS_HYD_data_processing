# ==========================================================================
# 1. SETUP, THEME, & PATHS (AUGUST VALIDATION ONLY)
# ==========================================================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(scales)

VALIDATION_MONTH <- 8  # August only

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

dir.create("Figures/lpj_guess_stem_storage/validation_august/PSI", recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# 2. RAW DATA INGESTION (AUGUST ONLY)
# ==========================================================================

lpj_raw <- read.csv("lpj_guess/lpj_guess_stem_storage/lpj_control_drought_ET_Gc_psi_leaf_psi_soil_psi_xylem_hydraulic_lag_kappy_s_min_mort_mort_cav_mort_greff_mort_min_stem_diameter_stem_rwc_twd.csv") %>%
  mutate(date = as.Date(date), month = month(date)) %>%
  filter(month == VALIDATION_MONTH) %>%   # AUGUST ONLY
  filter(!is.na(psi_leaf), !is.na(psi_xylem), !is.na(psi_soil), !is.na(species), !is.na(treatment)) %>%
  mutate(species = factor(species, levels = species_levels)) %>%
  select(date, species, treatment, psi_soil_model = psi_soil, psi_xylem_model = psi_xylem, psi_leaf_model = psi_leaf)

# Observed Soil Data
obs_soil <- bind_rows(
  read.csv("SCCII/psiS_hoelstein_drought.csv") %>% mutate(treatment = "drought"),
  read.csv("SCCII/psiS_hoelstein_control.csv") %>% mutate(treatment = "control")
) %>%
  mutate(date = as.Date(date), month = month(date)) %>%
  filter(month == VALIDATION_MONTH & !is.na(psiS_mean)) %>%   # AUGUST ONLY
  mutate(psi_soil_obs = psiS_mean / 1000) %>%
  select(date, treatment, psi_soil_obs)

# Observed Leaf Data
obs_leaf_raw <- bind_rows(
  read.csv("SCCII/psiL_hoelstein_drought.csv") %>% mutate(treatment = "drought"),
  read.csv("SCCII/psiL_hoelstein_control.csv") %>% mutate(treatment = "control")
) %>%
  mutate(date = as.Date(date), month = month(date)) %>%
  filter(month == VALIDATION_MONTH) %>%   # AUGUST ONLY
  rename(species = species_name) %>%
  mutate(species = factor(species, levels = species_levels))

obs_leaf_processed <- obs_leaf_raw %>%
  group_by(date, species, treatment) %>%
  summarise(
    psi_leaf_md_obs = mean(md_wp_av, na.rm = TRUE),
    psi_leaf_pd_obs = mean(pd_wp_av, na.rm = TRUE),
    .groups = "drop"
  )

# ==========================================================================
# 3. INTERSECTION MATRIX
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
# 4. TIME-SERIES PLOTS
# ==========================================================================
treatments <- unique(lpj_intersect$treatment)

for (t in treatments) {
  lpj_sub <- lpj_intersect[lpj_intersect$treatment == t, ]
  obs_sub <- obs_intersect[obs_intersect$treatment == t, ]

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
    facet_grid(. ~ species, scales = "free_y") +
    scale_color_manual(values = species_colors) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    ylim(-3.5, 0) +
    labs(
      title = paste("Water potential [AUGUST]:", t),
      subtitle = "AUGUST ONLY | Model: thick solid=ψl, dotdash=ψx, thin=ψs | Obs: ▲=ψl md, ■=ψl pd, ◆=ψs",
      x = "year", y = expression(psi ~ (mpa)), color = "Species"
    ) +
    base_theme

  ggsave(paste0("Figures/lpj_guess_stem_storage/validation_august/PSI/water_potential_", t, ".png"),
         plot_t, width = 16, height = 5, dpi = 300)
}

# ==========================================================================
# 5. 1:1 SCATTER PLOTS
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

safe_cor <- function(x, y) {
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) NA_real_ else cor(x, y, use = "complete.obs")
}

scatter_stats_psi <- flat_evaluation_pairs %>%
  group_by(species, treatment, pair_id) %>%
  summarise(
    n = n(),
    pearson_r = safe_cor(obs_val, model_val),
    pearson_r2 = safe_cor(obs_val, model_val)^2,
    rmse = sqrt(mean((model_val - obs_val)^2, na.rm = TRUE)),
    nrmse = (sqrt(mean((model_val - obs_val)^2, na.rm = TRUE)) / mean(obs_val, na.rm = TRUE)) * 100,
    bias = mean(model_val - obs_val, na.rm = TRUE),
    slope = if(n() > 1 && sd(obs_val, na.rm=TRUE) > 0) coef(lm(model_val ~ obs_val))[2] else NA,
    .groups = "drop"
  )

write.csv(scatter_stats_psi, "Figures/lpj_guess_stem_storage/validation_august/PSI/scatter_statistics_psi_5pairs.csv", row.names = FALSE)

create_scatter_psi_plot <- function(target_treatment, plot_title) {
  data_subset <- flat_evaluation_pairs %>% filter(treatment == target_treatment)
  stats_subset <- scatter_stats_psi %>% filter(treatment == target_treatment)

  pair_labels <- c(
    "psiL_vs_md" = "model ψl vs obs midday",
    "psiL_vs_pd" = "model ψl vs obs predawn",
    "psiS_vs_soil" = "model ψs vs obs soil",
    "psiX_vs_md" = "model ψx vs obs midday",
    "psiX_vs_pd" = "model ψx vs obs predawn"
  )

  data_subset <- data_subset %>% mutate(pair_label = pair_labels[pair_id])

  annotation_data <- stats_subset %>%
    mutate(
      pair_label = pair_labels[pair_id],
      text_summary = paste0(
        "n = ", n, "\n", "r = ", round(pearson_r, 2), "\n",
        "r² = ", round(pearson_r2, 2), "\n", "rmse = ", round(rmse, 2), "\n",
        "bias = ", round(bias, 2), "\n", "slope = ", round(slope, 1)
      )
    )

  axis_min <- min(c(data_subset$obs_val, data_subset$model_val), na.rm = TRUE) * 1.05

  ggplot(data_subset, aes(x = obs_val, y = model_val, color = species)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.1, linewidth = 0.6) +
    facet_grid(pair_label ~ species) +
    scale_color_manual(values = species_colors) +
    coord_fixed(ratio = 1, xlim = c(axis_min, 0), ylim = c(axis_min, 0)) +
    geom_text(data = annotation_data,
              aes(x = -Inf, y = Inf, label = text_summary),
              hjust = -0.05, vjust = 1.1, size = 2.8, color = "black", inherit.aes = FALSE) +
    labs(title = plot_title,
         subtitle = "AUGUST ONLY | dashed = 1:1, solid = lm",
         x = expression(observed ~ psi ~ (mpa)),
         y = expression(lpj ~ simulated ~ psi ~ (mpa))) +
    base_theme +
    theme(aspect.ratio = 1, strip.text.y = element_text(angle = -90, size = 9))
}

scatter_control <- create_scatter_psi_plot("control", "ψ 1:1 scatter: control [AUGUST]")
scatter_drought <- create_scatter_psi_plot("drought", "ψ 1:1 scatter: drought [AUGUST]")

ggsave("Figures/lpj_guess_stem_storage/validation_august/PSI/scatter_control_5pairs.png", scatter_control, width = 14, height = 15, dpi = 300)
ggsave("Figures/lpj_guess_stem_storage/validation_august/PSI/scatter_drought_5pairs.png", scatter_drought, width = 14, height = 15, dpi = 300)

# Dedicated psiL vs midday scatter
psiL_md_data <- flat_evaluation_pairs %>% filter(pair_id == "psiL_vs_md")
psiL_md_stats <- scatter_stats_psi %>% filter(pair_id == "psiL_vs_md") %>%
  mutate(text_label = paste0("n = ", n, "\n", "r = ", round(pearson_r, 2), "\n",
                              "R² = ", round(pearson_r2, 2), "\n",
                              "RMSE = ", round(rmse, 2), "\n", "slope = ", round(slope, 2)))

psiL_axis_min <- min(c(psiL_md_data$obs_val, psiL_md_data$model_val), na.rm = TRUE) * 1.05

p_psiL_md_scatter <- ggplot(psiL_md_data, aes(x = obs_val, y = model_val, color = species)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.7) +
  geom_text(data = psiL_md_stats,
            aes(x = -Inf, y = Inf, label = text_label),
            hjust = -0.05, vjust = 1.1, size = 3, color = "black", inherit.aes = FALSE) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = species_colors) +
  coord_fixed(ratio = 1, xlim = c(psiL_axis_min, 0), ylim = c(psiL_axis_min, 0)) +
  labs(title = "sim vs obs midday ψL [AUGUST]",
       subtitle = "AUGUST ONLY | dashed = 1:1 | solid = lm",
       x = expression(observed~midday~Psi[leaf]~(MPa)),
       y = expression(simulated~Psi[leaf]~(MPa))) +
  base_theme + theme(aspect.ratio = 1)

ggsave("Figures/lpj_guess_stem_storage/validation_august/PSI/scatter_psiL_vs_midday.png", p_psiL_md_scatter, width = 12, height = 8, dpi = 300)

cat("\n*** August PSI validation complete ***\n")
