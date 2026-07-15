# ==========================================================================
# 1. SETUP, THEME, & PATHS
# ==========================================================================
library(dplyr)
library(ggplot2)
library(tidyr)

# Standardize TWD: Z-score within species and treatment
lpj_plot_data <- read.csv("lpj_guess/lpj_guess_stem_storage/lpj_control_drought_ET_Gc_psi_leaf_psi_soil_psi_xylem_hydraulic_lag_kappy_s_min_mort_mort_cav_mort_greff_mort_min_stem_diameter_stem_rwc_twd.csv") %>%
  mutate(
    species = factor(tolower(species), levels = c("oak", "beech", "spruce", "pine")),
    treatment = factor(treatment, levels = c("control", "drought"))
  ) %>%
  group_by(species, treatment) %>%
  mutate(std_twd = (twd - mean(twd, na.rm = TRUE)) / sd(twd, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(species))

dir.create("Figures/lpj_guess_stem_storage/ET_VPD", recursive = TRUE, showWarnings = FALSE)

# Species color palette (consistent with other plot scripts)
species_colors <- c(
  "oak"    = "#E69F00",
  "beech"  = "#0072B2",
  "spruce" = "#009E73",
  "pine"   = "#F0E442"
)

# Base ggplot2 theme
base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.position   = "bottom",
    legend.text       = element_text(size = 11),
    plot.title        = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title        = element_text(size = 12),
    axis.text         = element_text(size = 10),
    axis.ticks        = element_line(color = "grey50"),
    strip.text        = element_text(size = 11, face = "bold"),
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_line(color = "grey92", linewidth = 0.25),
    panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )

# Relationships to plot: list(x, y, is_twd_std)
relationships <- list(
  list(x = "vpd", y = "ET", std = FALSE),
  list(x = "vpd", y = "Gc", std = FALSE),
  list(x = "psi_soil", y = "std_twd", std = TRUE),
  list(x = "psi_leaf", y = "std_twd", std = TRUE),
  list(x = "vpd", y = "std_twd", std = TRUE)
)

# ==========================================================================
# 2. GENERATE PLOTS
# ==========================================================================

for (rel in relationships) {
  # A. RAW SCATTER PLOT
  p_raw <- ggplot(lpj_plot_data, aes(x = .data[[rel$x]], y = .data[[rel$y]], color = species)) +
    geom_point(alpha = 0.1, size = 0.5) +
    # geom_smooth(method = "gam", formula = y ~ s(x, k = 3), color = "black", size = 0.6) +
    facet_grid(treatment ~ species) +
    scale_color_manual(values = species_colors) +
    labs(title = paste(rel$y, "vs", rel$x, "(Raw)"), x = rel$x, y = rel$y) +
    base_theme
  
  ggsave(paste0("Figures/lpj_guess_stem_storage/ET_VPD/raw_", rel$y, "_vs_", rel$x, ".png"), 
         p_raw, width = 10, height = 6, dpi = 300)
  
  # B. BINNED PLOT
  # Bin the X-axis into 20 equal intervals if not already categorical
  binned_data <- lpj_plot_data %>%
    filter(!is.na(.data[[rel$x]])) %>%
    mutate(x_bin = cut(.data[[rel$x]], breaks = 20)) %>%
    group_by(species, treatment, x_bin) %>%
    summarise(
      mean_y = mean(.data[[rel$y]], na.rm = TRUE),
      se_y = sd(.data[[rel$y]], na.rm = TRUE) / sqrt(n()),
      x_mid = mean(.data[[rel$x]], na.rm = TRUE),
      .groups = "drop"
    )
  
  p_binned <- ggplot(binned_data, aes(x = x_mid, y = mean_y, color = species, group = species)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_y - se_y, ymax = mean_y + se_y), width = 0.05) +
    facet_grid(treatment ~ species) +
    scale_color_manual(values = species_colors) +
    labs(title = paste(rel$y, "vs", rel$x, "(Binned Means)"), x = rel$x, y = paste("Mean", rel$y)) +
    base_theme
  
  ggsave(paste0("Figures/lpj_guess_stem_storage/ET_VPD/binned_", rel$y, "_vs_", rel$x, ".png"), 
         p_binned, width = 10, height = 6, dpi = 300)
}


# ==========================================================================
# 2. GENERATE PLOTS (UPDATED)
# ==========================================================================

for (rel in relationships) {
  # A. RAW SCATTER PLOT
  p_raw <- ggplot(lpj_plot_data, aes(x = .data[[rel$x]], y = .data[[rel$y]], color = species)) +
    geom_point(alpha = 0.1, size = 0.5) +
    facet_grid(treatment ~ species) +
    scale_color_manual(values = species_colors) +
    labs(title = paste(rel$y, "vs", rel$x, "(Raw)"), x = rel$x, y = rel$y) +
    base_theme
  
  ggsave(paste0("Figures/lpj_guess_stem_storage/ET_VPD/raw_", rel$y, "_vs_", rel$x, ".png"), 
         p_raw, width = 10, height = 6, dpi = 300)
  
  # B. BINNED PLOT
  binned_data <- lpj_plot_data %>%
    filter(!is.na(.data[[rel$x]])) %>%
    mutate(x_bin = cut(.data[[rel$x]], breaks = 20)) %>%
    group_by(species, treatment, x_bin) %>%
    summarise(
      mean_y = mean(.data[[rel$y]], na.rm = TRUE),
      se_y = sd(.data[[rel$y]], na.rm = TRUE) / sqrt(n()),
      x_mid = mean(.data[[rel$x]], na.rm = TRUE),
      .groups = "drop"
    )
  
  p_binned <- ggplot(binned_data, aes(x = x_mid, y = mean_y, color = species, group = species)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_y - se_y, ymax = mean_y + se_y), width = 0.05) +
    facet_grid(treatment ~ species) +
    scale_color_manual(values = species_colors) +
    labs(title = paste(rel$y, "vs", rel$x, "(Binned Means)"), x = rel$x, y = paste("Mean", rel$y)) +
    base_theme
  
  ggsave(paste0("Figures/lpj_guess_stem_storage/ET_VPD/binned_", rel$y, "_vs_", rel$x, ".png"), 
         p_binned, width = 10, height = 6, dpi = 300)
  
  # C. ALL SPECIES IN ONE PANEL (Combined)
  p_combined <- ggplot(binned_data, aes(x = x_mid, y = mean_y, color = species, group = species)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    # Facet only by treatment so species are overlaid in each treatment panel
    facet_wrap(~ treatment) + 
    scale_color_manual(values = species_colors) +
    labs(title = paste(rel$y, "vs", rel$x, "(All Species Combined)"), 
         x = rel$x, y = paste("Mean", rel$y)) +
    base_theme
  
  ggsave(paste0("Figures/lpj_guess_stem_storage/ET_VPD/combined_", rel$y, "_vs_", rel$x, ".png"), 
         p_combined, width = 8, height = 5, dpi = 300)
}