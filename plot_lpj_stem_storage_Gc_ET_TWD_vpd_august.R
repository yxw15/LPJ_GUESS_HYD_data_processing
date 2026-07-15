# ==========================================================================
# STEM STORAGE: Gc, ET, TWD vs VPD — AUGUST VALIDATION ONLY
# ==========================================================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)

VALIDATION_MONTH <- 8  # August only

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")
dir.create("Figures/lpj_guess_stem_storage/validation_august/ET_VPD", recursive = TRUE, showWarnings = FALSE)

lpj_plot_data <- read.csv("lpj_guess/lpj_guess_stem_storage/lpj_control_drought_ET_Gc_psi_leaf_psi_soil_psi_xylem_hydraulic_lag_kappy_s_min_mort_mort_cav_mort_greff_mort_min_stem_diameter_stem_rwc_twd.csv") %>%
  mutate(
    date = as.Date(date),
    species = factor(tolower(species), levels = c("oak", "beech", "spruce", "pine")),
    treatment = factor(treatment, levels = c("control", "drought"))
  ) %>%
  filter(month(date) == VALIDATION_MONTH) %>%
  group_by(species, treatment) %>%
  mutate(std_twd = (twd - mean(twd, na.rm = TRUE)) / sd(twd, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(species))

species_colors <- c("oak" = "#E69F00", "beech" = "#0072B2", "spruce" = "#009E73", "pine" = "#F0E442")

base_theme <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom", legend.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12), axis.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25)
  )

relationships <- list(
  list(x = "vpd", y = "ET", std = FALSE),
  list(x = "vpd", y = "Gc", std = FALSE),
  list(x = "psi_soil", y = "std_twd", std = TRUE),
  list(x = "psi_leaf", y = "std_twd", std = TRUE),
  list(x = "vpd", y = "std_twd", std = TRUE)
)

for (rel in relationships) {
  p <- ggplot(lpj_plot_data, aes_string(x = rel$x, y = rel$y, color = "species")) +
    geom_point(alpha = 0.4, size = 1.5) +
    stat_summary_bin(fun = "mean", geom = "point", size = 2.5) +
    stat_summary_bin(fun = "mean", geom = "line", linewidth = 0.9) +
    facet_wrap(~ treatment, ncol = 2) +
    scale_color_manual(values = species_colors) +
    labs(
      title = paste0("[AUGUST] ", rel$y, " vs ", rel$x),
      subtitle = "AUGUST VALIDATION",
      x = rel$x, y = ifelse(rel$std, paste0(rel$y, " (std)"), rel$y)
    ) +
    base_theme

  fname <- paste0("Figures/lpj_guess_stem_storage/validation_august/ET_VPD/",
                  rel$y, "_vs_", rel$x, ".png")
  ggsave(fname, p, width = 10, height = 7, dpi = 300)
}

cat("\n*** August Gc/ET/TWD vs VPD validation complete ***\n")
