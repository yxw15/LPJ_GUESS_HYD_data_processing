# ==============================================================================
# 1. SETUP & DATA LOADING
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

# ==============================================================================
# PHYSICAL CONSTANTS
# ==============================================================================

eta <- 44.6
T0  <- 273
h   <- 500

# ==============================================================================
# AESTHETIC CONFIGURATION
# ==============================================================================

species_order <- c("Oak", "Beech", "Spruce", "Pine")

cb_palette <- c(
  Oak   = "#E69F00",
  Beech = "#0072B2",
  Spruce = "#009E73",
  Pine  = "#F0E442"
)

pt_size <- 1.3

base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.position   = "bottom",
    legend.text       = element_text(size = 13),
    plot.title        = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title        = element_text(size = 16),
    axis.text         = element_text(size = 12),
    strip.text        = element_text(size = 14, face = "bold"),
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_line(color = "grey92", linewidth = 0.25),
    panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )

# ==============================================================================
# DATA INGESTION
# ==============================================================================

lpj_base <- read.csv("lpj_guess/lpj_Gc_psiL_psiX_psiS.csv") %>% mutate(date = as.Date(date))
psiL_control <- read.csv("SCCII/psiL_hoelstein_control.csv") %>% mutate(date = as.Date(date))
sap_control <- read.csv("SCCII/sap_hoelstein_control.csv") %>% mutate(date = as.Date(date))

# ==============================================================================
# CLIMATE FILTER
# ==============================================================================

climate_filter_dates <- lpj_base %>%
  filter(
    (temp_C + 273.15) > 287.15,
    precipitation < 1,
    global_radiation > 150,
    vpd > 0.3
  ) %>%
  pull(date) %>%
  unique()

climate_txt <- paste(
  "Temp > 14°C,", "Precip < 1 mm,", "Rad > 150 W/m²,", "VPD > 0.3 kPa"
)

# ==============================================================================
# MODEL DATA PROCESSING
# ==============================================================================

mod_proc <- lpj_base %>%
  mutate(
    gc_mmod = Gc * (eta * (T0 / (T0 + temp_C)) * exp(-0.00012 * h)) *1000,
  ) %>%
  select(date, species, gc_mmod, psiL = psiL) %>%
  mutate(species = factor(species, levels = species_order))


# ==============================================================================
# OBSERVATION PROCESSING
# ==============================================================================

# Midday leaf water potential
obs_md <- psiL_control %>%
  group_by(date, species_name) %>%
  summarise(psiL_md = mean(md_wp_av, na.rm = TRUE), .groups = "drop") %>%
  rename(species = species_name)

# Predawn leaf water potential
obs_pd <- psiL_control %>%
  group_by(date, species_name) %>%
  summarise(psiL_pd = mean(pd_wp_av, na.rm = TRUE), .groups = "drop") %>%
  rename(species = species_name)

# Sap flux conductance
obs_gc <- sap_control %>%
  group_by(date, species) %>%
  summarise(gc_obs = mean(G_asw, na.rm = TRUE), .groups = "drop")

# ==============================================================================
# COMBINE OBSERVATIONS
# ==============================================================================

obs_combined <- obs_gc %>%
  left_join(obs_md, by = c("date", "species")) %>%
  left_join(obs_pd, by = c("date", "species")) %>%
  mutate(species = factor(species, levels = species_order))

# ==============================================================================
# FILTERED OBSERVATIONS
# ==============================================================================

obs_filtered_climate <- obs_combined %>%
  filter(date %in% climate_filter_dates)

# ==============================================================================
# COMMON-TIME DATASET (both model and observations available)
# ==============================================================================

combined_data <- mod_proc %>%
  inner_join(obs_combined, by = c("date", "species")) %>%
  filter(
    date %in% climate_filter_dates,
    month(date) %in% c(6, 7, 8, 9)
  ) %>%
  na.omit()

# ==============================================================================
# FULL MODEL DATASET (climate-filtered, summer months only)
# ==============================================================================

data_full_model <- mod_proc %>%
  filter(
    date %in% climate_filter_dates,
    month(date) %in% c(6, 7, 8, 9)
  )

# ==============================================================================
# STANDARDIZATION (for common-time dataset)
# ==============================================================================

combined_std <- combined_data %>%
  group_by(species) %>%
  mutate(
    gc_max_obs = mean(gc_obs[gc_obs >= quantile(gc_obs, 0.90, na.rm = TRUE)], na.rm = TRUE),
    gc_max_mmod = mean(gc_mmod[gc_mmod >= quantile(gc_mmod, 0.90, na.rm = TRUE)], na.rm = TRUE),
    gc_rel_obs = gc_obs / gc_max_obs,
    gc_rel_mod = gc_mmod / gc_max_mmod
  ) %>%
  ungroup() %>%
  filter(gc_rel_obs >= 0, gc_rel_obs <= 1, gc_rel_mod >= 0, gc_rel_mod <= 1)

# ==============================================================================
# FULL-SERIES STANDARDIZATION
# ==============================================================================

mod_max_full <- data_full_model %>%
  group_by(species) %>%
  summarise(
    gc_max_mmod = mean(gc_mmod[gc_mmod >= quantile(gc_mmod, 0.90, na.rm = TRUE)], na.rm = TRUE)
  )

obs_max_full <- obs_filtered_climate %>%
  group_by(species) %>%
  summarise(
    gc_max_obs = mean(gc_obs[gc_obs >= quantile(gc_obs, 0.90, na.rm = TRUE)], na.rm = TRUE)
  )

data_full_model_std <- data_full_model %>%
  left_join(mod_max_full, by = "species") %>%
  mutate(gc_rel_mod = gc_mmod / gc_max_mmod) %>%
  filter(gc_rel_mod >= 0, gc_rel_mod <= 1)

data_obs_full_std <- obs_filtered_climate %>%
  left_join(obs_max_full, by = "species") %>%
  mutate(gc_rel_obs = gc_obs / gc_max_obs) %>%
  filter(gc_rel_obs >= 0, gc_rel_obs <= 1)

# ==============================================================================
# FIGURE 1 — ABSOLUTE Gc VS PSI_L (COMMON TIME)
# ==============================================================================

p_gc_psi_common <- ggplot(combined_data) +
  geom_point(aes(x = psiL_md, y = gc_obs), color = "black", alpha = 0.6, size = pt_size) +
  geom_point(aes(x = psiL_pd, y = gc_obs), color = "black", shape = 3, alpha = 0.8, size = pt_size) +
  geom_point(aes(x = psiL, y = gc_mmod, color = species), alpha = 0.7, size = pt_size) +
  facet_wrap(~species, ncol = 2) +
  scale_color_manual(values = cb_palette) +
  labs(
    title = "Conductance vs Leaf Water Potential (Common Time)",
    subtitle = paste0(climate_txt, " | June–Sept\nMidday = black dots | Predawn = black + | LPJ = colored"),
    x = expression(paste(Psi[L], " (MPa)")),
    y = expression(G[c]~(mol~m^{-2}~s^{-1})),
    color = "LPJ Species"
  ) +
  base_theme

print(p_gc_psi_common)

# ==============================================================================
# FIGURE 2 — ABSOLUTE Gc VS PSI_L (FULL SERIES)
# ==============================================================================

p_gc_psi_full <- ggplot() +
  geom_point(data = obs_filtered_climate, aes(x = psiL_md, y = gc_obs), 
             color = "black", alpha = 0.35, size = pt_size) +
  geom_point(data = obs_filtered_climate, aes(x = psiL_pd, y = gc_obs), 
             color = "black", shape = 3, alpha = 0.5, size = pt_size) +
  geom_point(data = data_full_model, aes(x = psiL, y = gc_mmod, color = species), 
             alpha = 0.45, size = pt_size) +
  facet_wrap(~species, ncol = 2) +
  scale_color_manual(values = cb_palette) +
  labs(
    title = "Conductance vs Leaf Water Potential (Full Time Series)",
    subtitle = "Midday = black dots | Predawn = black + | LPJ = colored",
    x = expression(paste(Psi[L], " (MPa)")),
    y = expression(G[c]~(mol~m^{-2}~s^{-1})),
    color = "LPJ Species"
  ) +
  base_theme

print(p_gc_psi_full)

# ==============================================================================
# FIGURE 3 — RELATIVE Gc/Gcmax (COMMON TIME)
# ==============================================================================

p_gc_rel_common <- ggplot(combined_std) +
  geom_point(aes(x = psiL_md, y = gc_rel_obs), color = "black", alpha = 0.6, size = pt_size) +
  geom_point(aes(x = psiL_pd, y = gc_rel_obs), color = "black", shape = 3, alpha = 0.8, size = pt_size) +
  geom_point(aes(x = psiL, y = gc_rel_mod, color = species), alpha = 0.7, size = pt_size) +
  facet_wrap(~species, ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Relative Conductance vs Leaf Water Potential (Common Time)",
    subtitle = "Midday = black dots | Predawn = black + | LPJ = colored",
    x = expression(paste(Psi[L], " (MPa)")),
    y = expression(G[c]/G[cmax]),
    color = "LPJ Species"
  ) +
  base_theme

print(p_gc_rel_common)

# ==============================================================================
# FIGURE 4 — RELATIVE Gc/Gcmax (FULL SERIES)
# ==============================================================================

p_gc_rel_full <- ggplot() +
  geom_point(data = data_obs_full_std, aes(x = psiL_md, y = gc_rel_obs), 
             color = "black", alpha = 0.35, size = pt_size) +
  geom_point(data = data_obs_full_std, aes(x = psiL_pd, y = gc_rel_obs), 
             color = "black", shape = 3, alpha = 0.5, size = pt_size) +
  geom_point(data = data_full_model_std, aes(x = psiL, y = gc_rel_mod, color = species), 
             alpha = 0.45, size = pt_size) +
  facet_wrap(~species, ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Relative Conductance vs Leaf Water Potential (Full Time Series)",
    subtitle = "Midday = black dots | Predawn = black + | LPJ = colored",
    x = expression(paste(Psi[L], " (MPa)")),
    y = expression(G[c]/G[cmax]),
    color = "LPJ Species"
  ) +
  base_theme

print(p_gc_rel_full)

# ==============================================================================
# FIGURE 5 — LPJ psiL vs OBSERVED MIDDAY psiL (COMMON TIME)
# ==============================================================================

psi_limits <- c(-4, 0)

p_psi_md_common <- ggplot(combined_data, aes(x = psiL_md, y = psiL, color = species)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_point(alpha = 0.7, size = pt_size) +
  facet_wrap(~species, ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_x_continuous(limits = psi_limits) +
  scale_y_continuous(limits = psi_limits) +
  coord_fixed() +
  labs(
    title = "LPJ vs Observed Midday Leaf Water Potential (Common Time)",
    subtitle = paste0("Dashed line = 1:1 relationship\n", climate_txt, " | June–Sept"),
    x = expression(paste("Observed ", Psi[L], " Midday (MPa)")),
    y = expression(paste("LPJ ", Psi[L], " (MPa)")),
    color = "Species"
  ) +
  base_theme

print(p_psi_md_common)

# ==============================================================================
# FIGURE 6 — LPJ psiL vs OBSERVED PREDAWN psiL (COMMON TIME)
# ==============================================================================

p_psi_pd_common <- ggplot(combined_data, aes(x = psiL_pd, y = psiL, color = species)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_point(alpha = 0.7, size = pt_size) +
  facet_wrap(~species, ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_x_continuous(limits = psi_limits) +
  scale_y_continuous(limits = psi_limits) +
  coord_fixed() +
  labs(
    title = "LPJ vs Observed Predawn Leaf Water Potential (Common Time)",
    subtitle = paste0("Dashed line = 1:1 relationship\n", climate_txt, " | June–Sept"),
    x = expression(paste("Observed ", Psi[L], " Predawn (MPa)")),
    y = expression(paste("LPJ ", Psi[L], " (MPa)")),
    color = "Species"
  ) +
  base_theme

print(p_psi_pd_common)

# ==============================================================================
# FIGURE 7 — Gc COMPARISON (ALL CONDITIONS)
# ==============================================================================

combined_data_withNA <- mod_proc %>%
  inner_join(obs_combined, by = c("date", "species"))

p_gc_all <- ggplot(combined_data_withNA, aes(x = gc_obs, y = gc_mmod, color = species)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_point(alpha = 0.6, size = pt_size) +
  facet_wrap(~species, ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_x_continuous(limits = c(0, 2000)) +
  scale_y_continuous(limits = c(0, 2000)) +
  coord_fixed() +
  labs(
    title = "LPJ vs Observed Conductance",
    subtitle = "X = G_asw | Y = LPJ gc_mmod (no climate filtering)\nDashed line = 1:1",
    x = expression(G[c]~Observed~(mol~m^{-2}~s^{-1})),
    y = expression(G[c]~LPJ~(mol~m^{-2}~s^{-1})),
    color = "Species"
  ) +
  base_theme

print(p_gc_all)

# ==============================================================================
# FIGURE 8 — Gc COMPARISON (CLIMATE-FILTERED)
# ==============================================================================

gc_compare_climate <- combined_data_withNA %>%
  filter(
    date %in% climate_filter_dates,
    month(date) %in% c(6, 7, 8, 9)
  ) %>%
  select(date, species, gc_mmod, gc_obs) %>%
  mutate(species = factor(species, levels = species_order))

p_gc_climate <- ggplot(gc_compare_climate, aes(x = gc_obs, y = gc_mmod, color = species)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.7, color = "black") +
  geom_point(alpha = 0.7, size = pt_size) +
  facet_wrap(~species, ncol = 2) +
  scale_color_manual(values = cb_palette) +
  coord_fixed() +
  labs(
    title = "LPJ vs Observed Conductance (Climate-Filtered Conditions)",
    subtitle = paste0("Filtered: ", climate_txt, "\nJune–September only\nDashed line = 1:1"),
    x = expression(G[c]~Observed~(mol~m^{-2}~s^{-1})),
    y = expression(G[c]~LPJ~(mol~m^{-2}~s^{-1})),
    color = "Species"
  ) +
  base_theme

print(p_gc_climate)

# ==============================================================================
# PRINT ALL FIGURES
# ==============================================================================

print(p_gc_psi_common)
print(p_gc_psi_full)
print(p_gc_rel_common)
print(p_gc_rel_full)
print(p_psi_md_common)
print(p_psi_pd_common)
print(p_gc_all)
print(p_gc_climate)

# ==============================================================================
# SAVE FIGURES
# ==============================================================================

ggsave("Figures/Hoelstein/compare_Gc_vs_PsiL_common_time.png", p_gc_psi_common, width = 13, height = 9, dpi = 300)
ggsave("Figures/Hoelstein/compare_Gc_vs_PsiL_full.png", p_gc_psi_full, width = 13, height = 9, dpi = 300)
ggsave("Figures/Hoelstein/compare_Gc_rel_vs_PsiL_common_time.png", p_gc_rel_common, width = 13, height = 9, dpi = 300)
ggsave("Figures/Hoelstein/compare_Gc_rel_vs_PsiL_full.png", p_gc_rel_full, width = 13, height = 9, dpi = 300)
ggsave("Figures/Hoelstein/compare_LPJ_vs_obs_md_common.png", p_psi_md_common, width = 13, height = 9, dpi = 300)
ggsave("Figures/Hoelstein/compare_LPJ_vs_obs_pd_common.png", p_psi_pd_common, width = 13, height = 9, dpi = 300)
ggsave("Figures/Hoelstein/compare_LPJ_vs_obs_Gc_all_summer.png", p_gc_all, width = 13, height = 9, dpi = 300)
ggsave("Figures/Hoelstein/compare_LPJ_vs_obs_Gc_climate_filtered.png", p_gc_climate, width = 13, height = 9, dpi = 300)

# ==============================================================================
# SINGLE PANEL VERSIONS (All species colored, unique shapes)
# ==============================================================================

# 1. Absolute Gc vs Psi_L (Common Time)
p_gc_psi_common_single <- ggplot(combined_data) +
  geom_point(aes(x = psiL_md, y = gc_obs, color = species, shape = "Obs Midday"), alpha = 0.8, size = pt_size) +
  geom_point(aes(x = psiL_pd, y = gc_obs, color = species, shape = "Obs Predawn"), alpha = 0.8, size = pt_size) +
  geom_point(aes(x = psiL, y = gc_mmod, color = species, shape = "LPJ Model"), alpha = 1, size = pt_size + 0.5) +
  scale_color_manual(values = cb_palette) +
  # Shape 4 = x, Shape 3 = +, Shape 16 = Solid Point
  scale_shape_manual(name = "Data Source", 
                     values = c("Obs Midday" = 4, "Obs Predawn" = 3, "LPJ Model" = 16)) +
  labs(
    title = "Conductance vs Leaf Water Potential (Common Time)",
    subtitle = "Color = Species | Midday = x | Predawn = + | LPJ = ●",
    x = expression(paste(Psi[L], " (MPa)")),
    y = expression(G[c]~(mol~m^{-2}~s^{-1})),
    color = "Species"
  ) +
  base_theme

print(p_gc_psi_common_single)

# 2. Absolute Gc vs Psi_L (Full Series)
p_gc_psi_full_single <- ggplot() +
  geom_point(data = obs_filtered_climate, aes(x = psiL_md, y = gc_obs, color = species, shape = "Obs Midday"), 
             alpha = 0.8, size = pt_size) +
  geom_point(data = obs_filtered_climate, aes(x = psiL_pd, y = gc_obs, color = species, shape = "Obs Predawn"), 
             alpha = 0.8, size = pt_size) +
  geom_point(data = data_full_model, aes(x = psiL, y = gc_mmod, color = species, shape = "LPJ Model"), 
             alpha = 0.6, size = pt_size) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "Data Source", 
                     values = c("Obs Midday" = 4, "Obs Predawn" = 3, "LPJ Model" = 16)) +
  labs(
    title = "Conductance vs Leaf Water Potential (Full Series)",
    subtitle = "Color = Species | Midday = x | Predawn = + | LPJ = ●",
    x = expression(paste(Psi[L], " (MPa)")),
    y = expression(G[c]~(mol~m^{-2}~s^{-1})),
    color = "Species"
  ) +
  base_theme

print(p_gc_psi_full_single)

# 3. Relative Gc (Common Time)
p_gc_rel_common_single <- ggplot(combined_std) +
  geom_point(aes(x = psiL_md, y = gc_rel_obs, color = species, shape = "Obs Midday"), alpha = 0.8, size = pt_size) +
  geom_point(aes(x = psiL_pd, y = gc_rel_obs, color = species, shape = "Obs Predawn"), alpha = 0.8, size = pt_size) +
  geom_point(aes(x = psiL, y = gc_rel_mod, color = species, shape = "LPJ Model"), alpha = 1, size = pt_size + 0.5) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "Data Source", 
                     values = c("Obs Midday" = 4, "Obs Predawn" = 3, "LPJ Model" = 16)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Relative Conductance vs Leaf Water Potential (Common Time)",
    subtitle = "Color = Species | Midday = x | Predawn = + | LPJ = ●",
    x = expression(paste(Psi[L], " (MPa)")),
    y = expression(G[c]/G[cmax]),
    color = "Species"
  ) +
  base_theme
print(p_gc_rel_common_single)

# 4. Relative Gc (Full Series)
p_gc_rel_full_single <- ggplot() +
  geom_point(data = data_obs_full_std, aes(x = psiL_md, y = gc_rel_obs, color = species, shape = "Obs Midday"), 
             alpha = 0.8, size = pt_size) +
  geom_point(data = data_obs_full_std, aes(x = psiL_pd, y = gc_rel_obs, color = species, shape = "Obs Predawn"), 
             alpha = 0.8, size = pt_size) +
  geom_point(data = data_full_model_std, aes(x = psiL, y = gc_rel_mod, color = species, shape = "LPJ Model"), 
             alpha = 0.6, size = pt_size) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "Data Source", 
                     values = c("Obs Midpoint" = 4, "Obs Predawn" = 3, "LPJ Model" = 16)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Relative Conductance vs Leaf Water Potential (Full Series)",
    subtitle = "Color = Species | Midday = x | Predawn = + | LPJ = ●",
    x = expression(paste(Psi[L], " (MPa)")),
    y = expression(G[c]/G[cmax]),
    color = "Species"
  ) +
  base_theme

print(p_gc_rel_full_single)

# ==============================================================================
# SAVE SINGLE PANEL FIGURES
# ==============================================================================

print(p_gc_psi_common_single)
print(p_gc_rel_common_single)
print(p_gc_psi_full_single)
print(p_gc_rel_full_single)

ggsave("Figures/Hoelstein/single_panel_Gc_vs_PsiL_common.png", p_gc_psi_common_single, width = 11, height = 8, dpi = 300)
ggsave("Figures/Hoelstein/single_panel_Gc_vs_PsiL_full.png", p_gc_psi_full_single, width = 11, height = 8, dpi = 300)
ggsave("Figures/Hoelstein/single_panel_Gc_rel_vs_PsiL_common.png", p_gc_rel_common_single, width = 11, height = 8, dpi = 300)
ggsave("Figures/Hoelstein/single_panel_Gc_rel_vs_PsiL_full.png", p_gc_rel_full_single, width = 11, height = 8, dpi = 300)
