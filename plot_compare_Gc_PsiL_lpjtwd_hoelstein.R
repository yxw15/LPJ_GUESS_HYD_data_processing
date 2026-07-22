# ==============================================================================
# 1. SETUP & DATA LOADING
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)

# Set working directory to the project root
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

# ==============================================================================
# 2. CONFIGURATIONS & AESTHETICS
# ==============================================================================

climate_txt <- paste(
  "temperature > 14°C,", "precipitation < 1 mm,", "global radiation > 150 W/m²,", "VPD > 0.3 kPa"
)

species_order = c("Oak", "Beech", "Spruce", "Pine")

cb_palette <- c(
  Oak    = "#E69F00",
  Beech  = "#0072B2",
  Spruce = "#009E73",
  Pine   = "#F0E442"
)

pt_size <- 1.3

base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.position   = "bottom",
    legend.text       = element_text(size = 13),
    plot.title        = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5, size = 11),
    axis.title        = element_text(size = 14),
    axis.text         = element_text(size = 12),
    strip.text        = element_text(size = 14, face = "bold"),
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_line(color = "grey92", linewidth = 0.25),
    panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )

# ==============================================================================
# 3. DATA INGESTION & FORMATTING
# ==============================================================================

lpj_output_filter <- read.csv("lpj_guess/lpj_guess_stem_storage/lpj_control_drought_ET_plant_ET_total_Gc_psi_leaf_psi_soil_psi_xylem_hydraulic_lag_kappy_s_min_mort_mort_cav_mort_greff_mort_min_stem_diameter_stem_rwc_twd_climate_filter.csv") %>%
  mutate(date = as.Date(date))

sap_flux_gc_filter <- read.csv("SCCII/sap_daily_filter.csv") %>%
  mutate(
    date = as.Date(date),
    treatment = case_when(
      treatment == "control" ~ "control",
      treatment == "treatment" ~ "drought",
      TRUE ~ treatment
    )
  )

obs_leaf_raw <- bind_rows(
  read.csv("SCCII/psiL_hoelstein_drought.csv") %>% mutate(treatment = "drought"),
  read.csv("SCCII/psiL_hoelstein_control.csv") %>% mutate(treatment = "control")
) %>%
  mutate(date = as.Date(date), month = month(date)) %>%
  filter(month >= 6 & month <= 9) %>%
  rename(species = if_else("species_name" %in% names(.), "species_name", "species")) %>%
  mutate(species = factor(species, levels = species_order))

# Safely extract climate context dates
climate_filter_dates <- unique(lpj_output_filter$date)

# ==============================================================================
# 4. DATA PROCESSING & ALIGNMENT
# ==============================================================================

# Clean and align model data
mod_proc <- lpj_output_filter %>%
  select(treatment, species, date, gc = Gc, psiL = psi_leaf) %>%
  mutate(species = factor(species, levels = species_order))

# Helper to catch edge cases if an unexpected column setup exists
mod_proc_monthly <- mod_proc %>%
  mutate(year = year(date), month = month(date)) %>%
  group_by(species, treatment, year, month) %>%
  summarise(gc = mean(gc, na.rm = TRUE), psiL = mean(psiL, na.rm = TRUE), .groups = "drop")

# Aggregate Leaf Water Potential (md_wp_av -> midday, pd_wp_av -> predawn)
obs_psi <- obs_leaf_raw %>%
  group_by(date, species, treatment) %>%
  summarise(
    psiL_md = mean(md_wp_av, na.rm = TRUE),
    psiL_pd = mean(pd_wp_av, na.rm = TRUE),
    .groups = "drop"
  )

# Aggregate Sap Flux Conductance (G_ms is used since LPJ runs in m/s units)
obs_gc <- sap_flux_gc_filter %>%
  group_by(date, species, treatment) %>%
  summarise(
    gc_obs = mean(G_ms, na.rm = TRUE), 
    .groups = "drop"
  )

# Combine daily observations
obs_combined <- obs_gc %>%
  left_join(obs_psi, by = c("date", "species", "treatment")) %>%
  mutate(species = factor(species, levels = species_order))

# Filter observations to target seasonal limits (Full Series Context)
obs_filtered_climate <- obs_combined %>%
  filter(date %in% climate_filter_dates, month(date) %in% c(6, 7, 8, 9))

# Combine datasets for Common-Time bounds with source column
combined_data_lpj_obs <- mod_proc %>%
  inner_join(obs_filtered_climate, by = c("date", "species", "treatment"))

# ==============================================================================
# 5. MINIMUM VALUE STANDARDIZATION & AGGREGATION
# ==============================================================================

# 1. Common Time Dataset: long format with source column, 90% quantile baseline
# Stack observed (midday + predawn) and simulated into long format
combined_long <- bind_rows(
  # Observed midday
  combined_data_lpj_obs %>%
    select(date, species, treatment, gc = gc_obs, psiL = psiL_md) %>%
    mutate(source = "obs", psiL_label = "midday"),
  # Observed predawn
  combined_data_lpj_obs %>%
    select(date, species, treatment, gc = gc_obs, psiL = psiL_pd) %>%
    mutate(source = "obs", psiL_label = "predawn"),
  # Simulated
  combined_data_lpj_obs %>%
    select(date, species, treatment, gc = gc, psiL = psiL) %>%
    mutate(source = "sim", psiL_label = "simulated")
)

# Drop rows where psiL is NA (no x-axis position, would skew min/max)
combined_long <- combined_long %>% filter(!is.na(psiL))

# Standardize: true min-max per species × treatment × source × psiL_label
combined_std_lpj_obs <- combined_long %>%
  group_by(species, treatment, source, psiL_label) %>%
  mutate(
    gc_min = min(gc, na.rm = TRUE),
    gc_max = max(gc, na.rm = TRUE),
    gc_rel = (gc - gc_min) / (gc_max - gc_min)
  ) %>%
  ungroup()

# --- Diagnostic: verify each group spans 0-1 ---
cat("\n=== Standardization range check (should all be min=0, max=1) ===\n")
combined_std_lpj_obs %>%
  group_by(species, treatment, source, psiL_label) %>%
  summarise(
    gc_abs_min = min(gc, na.rm = TRUE),
    gc_abs_max = max(gc, na.rm = TRUE),
    gc_rel_min = min(gc_rel, na.rm = TRUE),
    gc_rel_max = max(gc_rel, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  print(n = 100)

# Print out explicit counts across the structured grid matrix
point_counts <- combined_data_lpj_obs %>%
  group_by(species, treatment) %>%
  summarise(
    obs_midday   = sum(!is.na(psiL_md) & !is.na(gc_obs)),
    obs_predawn  = sum(!is.na(psiL_pd) & !is.na(gc_obs)),
    lpj_model    = sum(!is.na(psiL)    & !is.na(gc)),
    .groups = "drop"
  )

print("--- Data Point Counts per Facet Panel ---")
print(point_counts)

# 2. Full Simulated Series Standardization (true min-max per species × treatment)
data_full_model_std <- mod_proc %>%
  filter(date %in% climate_filter_dates, month(date) %in% c(6, 7, 8, 9)) %>%
  group_by(species, treatment) %>%
  mutate(
    gc_min = min(gc, na.rm = TRUE),
    gc_max = max(gc, na.rm = TRUE),
    gc_rel_mod = (gc - gc_min) / (gc_max - gc_min)
  ) %>%
  ungroup()

# 3. Full Observed Series: all summer obs (not restricted to LPJ dates), long format
data_obs_full_std <- obs_combined %>%
  filter(month(date) %in% c(6, 7, 8, 9))

data_obs_full_std <- bind_rows(
  data_obs_full_std %>%
    select(date, species, treatment, gc = gc_obs, psiL = psiL_md) %>%
    mutate(psiL_label = "midday"),
  data_obs_full_std %>%
    select(date, species, treatment, gc = gc_obs, psiL = psiL_pd) %>%
    mutate(psiL_label = "predawn")
) %>%
  mutate(source = "obs") %>%
  filter(!is.na(psiL)) %>%
  group_by(species, treatment, psiL_label) %>%
  mutate(
    gc_min = min(gc, na.rm = TRUE),
    gc_max = max(gc, na.rm = TRUE),
    gc_rel = (gc - gc_min) / (gc_max - gc_min)
  ) %>%
  ungroup()

# 4. Full Monthly Simulated Series Standardization
data_full_model_monthly_raw <- mod_proc %>%
  filter(date %in% climate_filter_dates, month(date) %in% c(6, 7, 8, 9)) %>%
  mutate(
    year  = year(date),
    month = month(date)
  ) %>%
  group_by(species, treatment, year, month) %>%
  summarise(
    gc   = mean(gc, na.rm = TRUE),
    psiL = mean(psiL, na.rm = TRUE),
    .groups = "drop"
  )

data_full_model_monthly_std <- data_full_model_monthly_raw %>%
  group_by(species, treatment) %>%
  mutate(
    gc_min     = min(gc, na.rm = TRUE),
    gc_max     = max(gc, na.rm = TRUE),
    gc_rel_mod = (gc - gc_min) / (gc_max - gc_min)
  ) %>%
  ungroup()

# ==============================================================================
# 6. VISUALIZATIONS: FACETED GRID DESIGNS (Rows = Tx, Cols = Sp)
# ==============================================================================

# Figure 1: Absolute Gc vs Psi_L (Common Time)
p_gc_psi_common <- ggplot(combined_data_lpj_obs) +
  geom_point(aes(x = psiL_md, y = gc_obs), shape = 2, color = "black", alpha = 0.8, size = pt_size, na.rm = TRUE) +
  geom_point(aes(x = psiL_pd, y = gc_obs), color = "black", shape = 3, alpha = 0.8, size = pt_size, na.rm = TRUE) +
  geom_point(aes(x = psiL, y = gc, color = species), alpha = 1, size = pt_size, na.rm = TRUE) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = cb_palette) +
  ylim(0, 12) +
  labs(
    title = "canopy conductance vs leaf water potential (common time)",
    subtitle = tolower(paste0(climate_txt, " | june–september\nmidday = open triangle | predawn = black + | lpj = colored")),
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]~(m~s^{-1})),
    color = ""
  ) + base_theme

print(p_gc_psi_common)

# Figure 2: Absolute Gc vs Psi_L (Full Series using monthly averages)
p_gc_psi_full <- ggplot() +
  geom_point(data = obs_filtered_climate, aes(x = psiL_md, y = gc_obs), shape = 2, color = "black", alpha = 0.8, size = pt_size, na.rm = TRUE) +
  geom_point(data = obs_filtered_climate, aes(x = psiL_pd, y = gc_obs), color = "black", shape = 3, alpha = 0.8, size = pt_size, na.rm = TRUE) +
  geom_point(data = mod_proc_monthly, aes(x = psiL, y = gc, color = species), alpha = 0.6, size = pt_size, na.rm = TRUE) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = cb_palette) +
  ylim(0, 12) +
  labs(
    title = "canopy conductance vs leaf water potential (monthly modelled means)",
    subtitle = "midday = open triangle | predawn = black + | lpj = colored monthly mean",
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]~(m~s^{-1})),
    color = ""
  ) + base_theme

print(p_gc_psi_full)

# Figure 3: Relative Gc/Gcmax (Common Time, 90% quantile per species × treatment × source)
p_gc_rel_common <- ggplot() +
  geom_point(
    data = combined_std_lpj_obs %>% filter(source == "obs", psiL_label == "midday"),
    aes(x = psiL, y = gc_rel), shape = 2, color = "black", alpha = 0.8, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = combined_std_lpj_obs %>% filter(source == "obs", psiL_label == "predawn"),
    aes(x = psiL, y = gc_rel), color = "black", shape = 3, alpha = 0.8, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = combined_std_lpj_obs %>% filter(source == "sim"),
    aes(x = psiL, y = gc_rel, color = species), alpha = 1, size = pt_size + 0.5, na.rm = TRUE
  ) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = cb_palette) +
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "standardized canopy conductance vs leaf water potential (common time)",
    subtitle = "true min-max per species × treatment × source × psiL_label\nmidday = open triangle | predawn = black + | lpj = colored",
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]/G[cmax]),
    color = ""
  ) + base_theme

print(p_gc_rel_common)

# Figure 4: Relative Gc/Gcmax (Full Series, true min-max per species × treatment × psiL_label)
p_gc_rel_full <- ggplot() +
  geom_point(
    data = data_obs_full_std %>% filter(psiL_label == "midday"),
    aes(x = psiL, y = gc_rel), shape = 2, color = "black", alpha = 0.8, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = data_obs_full_std %>% filter(psiL_label == "predawn"),
    aes(x = psiL, y = gc_rel), color = "black", shape = 3, alpha = 0.8, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = data_full_model_monthly_std,
    aes(x = psiL, y = gc_rel_mod, color = species),
    alpha = 0.65, size = pt_size + 0.5, na.rm = TRUE
  ) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = cb_palette) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "standardized canopy conductance vs leaf water potential (full monthly mean)",
    subtitle = "true min-max per species × treatment × psiL_label | lpj = colored monthly mean",
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]/G[cmax]),
    color = ""
  ) + base_theme

print(p_gc_rel_full)

# ==============================================================================
# 7. VISUALIZATIONS: SINGLE PANEL WRAP DESIGNS (control vs drought)
# ==============================================================================

# Single Panel 1: Absolute Values (Common Time)
p_gc_psi_common_single <- ggplot(combined_data_lpj_obs) +
  geom_point(aes(x = psiL_md, y = gc_obs, color = species, shape = "obs midday"), alpha = 0.6, size = pt_size, na.rm = TRUE) +
  geom_point(aes(x = psiL_pd, y = gc_obs, color = species, shape = "obs predawn"), alpha = 0.6, size = pt_size, na.rm = TRUE) +
  geom_point(aes(x = psiL, y = gc, color = species, shape = "lpj simulated"), alpha = 1, size = pt_size + 0.5, na.rm = TRUE) +
  facet_wrap(vars(treatment), ncol = 2) +
  ylim(0, 12) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "", values = c("obs midday" = 2, "obs predawn" = 3, "lpj simulated" = 16)) +
  labs(
    title = "canopy conductance vs leaf water potential (common time)",
    subtitle = tolower(paste0(climate_txt, "\nmidday = open triangle | predawn = + | lpj-guess-hyd = ●")),
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]~(m~s^{-1})),
    color = "species"
  ) + base_theme

print(p_gc_psi_common_single)

## Single Four panels: : Absolute Values (Common Time)
# Observations - midday
obs_md <- combined_data_lpj_obs %>%
  select(treatment, species,
         gc = gc_obs,
         psi = psiL_md) %>%
  mutate(source = "Observation",
         type = "Midday")

# Observations - predawn
obs_pd <- combined_data_lpj_obs %>%
  select(treatment, species,
         gc = gc_obs,
         psi = psiL_pd) %>%
  mutate(source = "Observation",
         type = "Predawn")

# Simulation
sim_dat <- combined_data_lpj_obs %>%
  select(treatment, species,
         gc,
         psi = psiL) %>%
  mutate(source = "Simulation",
         type = "Simulation")

plot_dat <- bind_rows(sim_dat, obs_md, obs_pd) %>%
  mutate(
    source = factor(source,
                    levels = c("Simulation", "Observation"))
  )

p_gc_psi_common_single_four <- ggplot(plot_dat,
       aes(x = psi, y = gc,
           color = species,
           shape = type)) +
  geom_point(size = pt_size, alpha = 0.7, na.rm = TRUE) +
  facet_grid(source ~ treatment) +
  scale_shape_manual(
    values = c(
      "Midday" = 15,       # solid sqaure
      "Predawn" = 2,      # open triangle
      "Simulation" = 16   # solid circle
    )
  ) +
  scale_color_manual(values = cb_palette) +
  ylim(0, 12) +
  labs(
    shape = "",
    color = "Species",
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]~(m~s^{-1}))
  ) +
  base_theme
print(p_gc_psi_common_single_four)

# Single Panel 2: Standardized Values (Common Time)
p_gc_rel_common_single <- ggplot() +
  geom_point(
    data = combined_std_lpj_obs %>% filter(source == "obs", psiL_label == "midday"),
    aes(x = psiL, y = gc_rel, color = species, shape = "obs midday"), alpha = 0.6, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = combined_std_lpj_obs %>% filter(source == "obs", psiL_label == "predawn"),
    aes(x = psiL, y = gc_rel, color = species, shape = "obs predawn"), alpha = 0.6, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = combined_std_lpj_obs %>% filter(source == "sim"),
    aes(x = psiL, y = gc_rel, color = species, shape = "lpj simulated"), alpha = 1, size = pt_size + 0.5, na.rm = TRUE
  ) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "", values = c("obs midday" = 2, "obs predawn" = 3, "lpj simulated" = 16)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "standardized canopy conductance vs leaf water potential (common time)",
    subtitle = tolower(paste0(climate_txt, " (true min-max per species × treatment × source × psiL_label)\nmidday = open triangle | predawn = + | lpj-guess-hyd = ●")),
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]/G[cmax]),
    color = "species"
  ) + base_theme

print(p_gc_rel_common_single)

# Single four panels: Standardized Values (Common Time)
plot_dat <- combined_std_lpj_obs %>%
  mutate(
    source = case_when(
      source == "sim" ~ "Simulation",
      source == "obs" ~ "Observation"
    ),
    shape_type = case_when(
      source == "Simulation"  ~ "Simulation",
      psiL_label == "midday"  ~ "Midday",
      psiL_label == "predawn" ~ "Predawn"
    ),
    source = factor(
      source,
      levels = c("Simulation", "Observation")
    )
  )

p_gc_rel_common_single_four <- ggplot(
  plot_dat,
  aes(
    x = psiL,
    y = gc_rel,
    color = species,
    shape = shape_type
  )
) +
  geom_point(
    alpha = 0.7,
    size = pt_size,
    na.rm = TRUE
  ) +
  facet_grid(source ~ treatment) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(
    name = "",
    values = c(
      "Midday" = 15,      # solid square
      "Predawn" = 2,      # open triangle
      "Simulation" = 16   # solid circle
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    title = "standardized canopy conductance vs leaf water potential (common time)",
    subtitle = tolower(
      paste0(
        climate_txt,
        " (true min-max per species × treatment × source × psiL_label)\n",
        "midday = ■ | predawn = △ | simulation = ●"
      )
    ),
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] / G[cmax]),
    color = "species"
  ) +
  base_theme

print(p_gc_rel_common_single_four)

# ==============================================================================
# 7b. SINGLE PANEL DESIGNS — FULL TIME (Monthly Mean & Daily)
# ==============================================================================

# Single Panel 3: Absolute Values (Full Monthly Mean)
p_gc_psi_full_single_monthly_mean <- ggplot() +
  geom_point(data = obs_filtered_climate, aes(x = psiL_md, y = gc_obs, color = species, shape = "obs midday"), alpha = 0.5, size = pt_size, na.rm = TRUE) +
  geom_point(data = obs_filtered_climate, aes(x = psiL_pd, y = gc_obs, color = species, shape = "obs predawn"), alpha = 0.5, size = pt_size, na.rm = TRUE) +
  geom_point(data = mod_proc_monthly, aes(x = psiL, y = gc, color = species, shape = "lpj monthly mean"), alpha = 0.8, size = pt_size + 0.5, na.rm = TRUE) +
  facet_wrap(vars(treatment), ncol = 2) +
  ylim(0, 12) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "", values = c("obs midday" = 2, "obs predawn" = 3, "lpj monthly mean" = 16)) +
  labs(
    title = "canopy conductance vs leaf water potential (full monthly mean)",
    subtitle = "midday = open triangle | predawn = + | lpj-guess-hyd monthly mean = ●",
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]~(m~s^{-1})),
    color = "species"
  ) + base_theme
print(p_gc_psi_full_single_monthly_mean)

# Single Panel 4: Standardized Values (Full Monthly Mean)
p_gc_rel_full_single_monthly_mean <- ggplot() +
  geom_point(
    data = data_obs_full_std %>% filter(psiL_label == "midday"),
    aes(x = psiL, y = gc_rel, color = species, shape = "obs midday"), alpha = 0.5, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = data_obs_full_std %>% filter(psiL_label == "predawn"),
    aes(x = psiL, y = gc_rel, color = species, shape = "obs predawn"), alpha = 0.5, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = data_full_model_monthly_std,
    aes(x = psiL, y = gc_rel_mod, color = species, shape = "lpj monthly mean"), alpha = 0.8, size = pt_size + 0.5, na.rm = TRUE
  ) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "", values = c("obs midday" = 2, "obs predawn" = 3, "lpj monthly mean" = 16)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "standardized canopy conductance vs leaf water potential (full monthly mean)",
    subtitle = "true min-max per species × treatment × psiL_label | lpj = colored monthly mean",
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]/G[cmax]),
    color = "species"
  ) + base_theme
print(p_gc_rel_full_single_monthly_mean)

# Single four panels: Standardized Values (Full Monthly Mean)
# Observation data
obs_md <- obs_filtered_climate %>%
  mutate(
    source = "Observation",
    shape_type = "Midday",
    psi = psiL_md,
    gc_plot = gc_obs
  ) %>%
  select(treatment, species, source, shape_type, psi, gc_plot)

obs_pd <- obs_filtered_climate %>%
  mutate(
    source = "Observation",
    shape_type = "Predawn",
    psi = psiL_pd,
    gc_plot = gc_obs
  ) %>%
  select(treatment, species, source, shape_type, psi, gc_plot)

# Simulation data
sim_dat <- mod_proc_monthly %>%
  mutate(
    source = "Simulation",
    shape_type = "Simulation",
    psi = psiL,
    gc_plot = gc
  ) %>%
  select(treatment, species, source, shape_type, psi, gc_plot)

# Combined dataframe
plot_dat <- bind_rows(sim_dat, obs_md, obs_pd) %>%
  mutate(
    source = factor(source,
                    levels = c("Simulation", "Observation")),
    treatment = factor(treatment,
                       levels = c("control", "drought"))
  )

p_gc_psi_full_single_monthly_mean_four <- ggplot(
  plot_dat,
  aes(x = psi,
      y = gc_plot,
      color = species,
      shape = shape_type)
) +
  geom_point(
    alpha = 0.7,
    size = pt_size,
    na.rm = TRUE
  ) +
  facet_grid(source ~ treatment) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(
    name = "",
    values = c(
      "Midday" = 15,      # solid square
      "Predawn" = 2,      # open triangle
      "Simulation" = 16   # solid circle
    )
  ) +
  coord_cartesian(ylim = c(0, 12)) +
  labs(
    title = "canopy conductance vs leaf water potential (full monthly mean)",
    subtitle = "midday = ■ | predawn = △ | lpj monthly mean = ●",
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] ~ (m~s^{-1})),
    color = "species"
  ) +
  base_theme

print(p_gc_psi_full_single_monthly_mean_four)


# Single Panel 5: Absolute Values (Full Daily Range)
p_gc_psi_full_single_daily <- ggplot() +
  geom_point(data = obs_filtered_climate, aes(x = psiL_md, y = gc_obs, color = species, shape = "obs midday"), alpha = 0.5, size = pt_size, na.rm = TRUE) +
  geom_point(data = obs_filtered_climate, aes(x = psiL_pd, y = gc_obs, color = species, shape = "obs predawn"), alpha = 0.5, size = pt_size, na.rm = TRUE) +
  geom_point(data = data_full_model_std, aes(x = psiL, y = gc, color = species, shape = "lpj daily simulated"), alpha = 0.3, size = pt_size, na.rm = TRUE) +
  ylim(0, 12) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "", values = c("obs midday" = 2, "obs predawn" = 3, "lpj daily simulated" = 16)) +
  labs(
    title = "canopy conductance vs leaf water potential (full daily data)",
    subtitle = "midday = open triangle | predawn = + | lpj daily simulated points = ●",
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]~(m~s^{-1})),
    color = "species"
  ) + base_theme
print(p_gc_psi_full_single_daily)

# Single four panels: Absolute Values (Full Daily Range)
# Observation data
obs_md <- obs_filtered_climate %>%
  mutate(
    source = "Observation",
    shape_type = "Midday",
    psi = psiL_md,
    gc_plot = gc_obs
  ) %>%
  select(treatment, species, source, shape_type, psi, gc_plot)

obs_pd <- obs_filtered_climate %>%
  mutate(
    source = "Observation",
    shape_type = "Predawn",
    psi = psiL_pd,
    gc_plot = gc_obs
  ) %>%
  select(treatment, species, source, shape_type, psi, gc_plot)

# Daily simulation data
sim_dat <- data_full_model_std %>%
  mutate(
    source = "Simulation",
    shape_type = "Simulation",
    psi = psiL,
    gc_plot = gc
  ) %>%
  select(treatment, species, source, shape_type, psi, gc_plot)

# Combined data
plot_dat <- bind_rows(sim_dat, obs_md, obs_pd) %>%
  mutate(
    source = factor(
      source,
      levels = c("Simulation", "Observation")
    ),
    treatment = factor(
      treatment,
      levels = c("control", "drought")
    )
  )

p_gc_psi_full_single_daily_four <- ggplot(
  plot_dat,
  aes(
    x = psi,
    y = gc_plot,
    color = species,
    shape = shape_type
  )
) +
  geom_point(
    alpha = 0.5,
    size = pt_size,
    na.rm = TRUE
  ) +
  facet_grid(source ~ treatment) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(
    name = "",
    values = c(
      "Midday" = 15,      # solid square
      "Predawn" = 2,      # open triangle
      "Simulation" = 16   # solid circle
    )
  ) +
  coord_cartesian(ylim = c(0, 12)) +
  labs(
    title = "canopy conductance vs leaf water potential (full daily data)",
    subtitle = "midday = ■ | predawn = △ | lpj daily simulated = ●",
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] ~ (m~s^{-1})),
    color = "species"
  ) +
  base_theme

print(p_gc_psi_full_single_daily_four)

# Single Panel 6: Standardized Values (Full Daily Range)
p_gc_rel_full_single_daily <- ggplot() +
  geom_point(
    data = data_obs_full_std %>% filter(psiL_label == "midday"),
    aes(x = psiL, y = gc_rel, color = species, shape = "obs midday"), alpha = 0.5, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = data_obs_full_std %>% filter(psiL_label == "predawn"),
    aes(x = psiL, y = gc_rel, color = species, shape = "obs predawn"), alpha = 0.5, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = data_full_model_std,
    aes(x = psiL, y = gc_rel_mod, color = species, shape = "lpj daily simulated"), alpha = 0.3, size = pt_size, na.rm = TRUE
  ) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "", values = c("obs midday" = 2, "obs predawn" = 3, "lpj daily simulated" = 16)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "standardized canopy conductance vs leaf water potential (full daily data)",
    subtitle = "true min-max per species × treatment × psiL_label",
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]/G[cmax]),
    color = "species"
  ) + base_theme
print(p_gc_rel_full_single_daily)

# Single four panels: Standardized Values (Full Daily Range)
# Observation data
obs_md <- data_obs_full_std %>%
  filter(psiL_label == "midday") %>%
  mutate(
    source = "Observation",
    shape_type = "Midday",
    gc_plot = gc_rel
  ) %>%
  select(treatment, species, source, shape_type, psiL, gc_plot)

obs_pd <- data_obs_full_std %>%
  filter(psiL_label == "predawn") %>%
  mutate(
    source = "Observation",
    shape_type = "Predawn",
    gc_plot = gc_rel
  ) %>%
  select(treatment, species, source, shape_type, psiL, gc_plot)

# Simulation data
sim_dat <- data_full_model_std %>%
  mutate(
    source = "Simulation",
    shape_type = "Simulation",
    gc_plot = gc_rel_mod
  ) %>%
  select(treatment, species, source, shape_type, psiL, gc_plot)

# Combined data
plot_dat <- bind_rows(sim_dat, obs_md, obs_pd) %>%
  mutate(
    source = factor(
      source,
      levels = c("Simulation", "Observation")
    ),
    treatment = factor(
      treatment,
      levels = c("control", "drought")
    )
  )

p_gc_rel_full_single_daily_four <- ggplot(
  plot_dat,
  aes(
    x = psiL,
    y = gc_plot,
    color = species,
    shape = shape_type
  )
) +
  geom_point(
    alpha = 0.5,
    size = pt_size,
    na.rm = TRUE
  ) +
  facet_grid(source ~ treatment) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(
    name = "",
    values = c(
      "Midday" = 15,      # solid square
      "Predawn" = 2,      # open triangle
      "Simulation" = 16   # solid circle
    )
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2)
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "standardized canopy conductance vs leaf water potential (full daily data)",
    subtitle = paste0(
      "true min-max per species × treatment × psiL_label\n",
      "midday = ■ | predawn = △ | lpj daily simulated = ●"
    ),
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] / G[cmax]),
    color = "species"
  ) +
  base_theme

print(p_gc_rel_full_single_daily_four)

# ==============================================================================
# 7c. THEIL-SEN SLOPE — SIMULATED DATA (absolute & relative)
# ==============================================================================

# Theil-Sen slope estimator (median of pairwise slopes)
theil_sen <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(NA_real_)
  n <- length(x)
  slopes <- c()
  for (i in 1:(n - 1)) {
    dx <- x[(i + 1):n] - x[i]
    slopes_i <- (y[(i + 1):n] - y[i]) / dx
    slopes <- c(slopes, slopes_i[is.finite(slopes_i)])
  }
  median(slopes)
}

# --- Absolute Gc Theil-Sen slopes (simulated, common-time) ---
ts_abs <- combined_data_lpj_obs %>%
  group_by(species, treatment) %>%
  summarise(
    ts_slope  = theil_sen(psiL, gc),
    n         = sum(!is.na(psiL) & !is.na(gc)),
    .groups   = "drop"
  )

print("=== Theil-Sen slopes — absolute Gc (simulated) ===")
print(as.data.frame(ts_abs))

p_ts_abs <- ggplot(ts_abs, aes(x = species, y = ts_slope, fill = species)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_fill_manual(values = cb_palette, guide = "none") +
  geom_text(aes(label = round(ts_slope, 3)), vjust = -0.4, size = 3.5) +
  labs(
    title = "Theil-Sen slope — absolute Gc vs PsiL (LPJ simulated, common time)",
    subtitle = "median of pairwise slopes",
    x = "Species",
    y = expression("Theil-Sen slope (m·s"^{-1}~"·MPa"^{-1}~")")
  ) + base_theme

print(p_ts_abs)

# --- Relative Gc/Gcmax Theil-Sen slopes (simulated, common-time) ---
ts_rel <- combined_std_lpj_obs %>%
  filter(source == "sim") %>%
  group_by(species, treatment) %>%
  summarise(
    ts_slope  = theil_sen(psiL, gc_rel),
    n         = n(),
    .groups   = "drop"
  )

print("=== Theil-Sen slopes — relative Gc/Gcmax (simulated) ===")
print(as.data.frame(ts_rel))

p_ts_rel <- ggplot(ts_rel, aes(x = species, y = ts_slope, fill = species)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_fill_manual(values = cb_palette, guide = "none") +
  geom_text(aes(label = round(ts_slope, 3)), vjust = -0.4, size = 3.5) +
  labs(
    title = "Theil-Sen slope — relative Gc/Gcmax vs PsiL (LPJ simulated, common time)",
    subtitle = "median of pairwise slopes | true min-max standardized",
    x = "Species",
    y = expression("Theil-Sen slope (MPa"^{-1}~")")
  ) + base_theme

print(p_ts_rel)

# ==============================================================================
# 8. BINNED MODELLED Gc vs PsiL (LPJ ONLY, NO OBSERVED DATA)
# ==============================================================================

# Bin width for psiL (MPa)
psiL_bin_width <- 0.2

# Bin the modelled data by psiL, compute mean Gc per bin
mod_binned <- data_full_model_std %>%
  mutate(psiL_bin = round(psiL / psiL_bin_width) * psiL_bin_width) %>%
  group_by(treatment, species, psiL_bin) %>%
  summarise(
    gc_mean  = mean(gc, na.rm = TRUE),
    gc_sd    = sd(gc, na.rm = TRUE),
    gc_se    = sd(gc, na.rm = TRUE) / sqrt(n()),
    psiL_bin_mean = mean(psiL, na.rm = TRUE),
    n        = n(),
    .groups  = "drop"
  ) %>%
  filter(n >= 3)  # exclude bins with too few points

# Print bin summary for diagnostics
cat("\n=== binned modelled Gc vs PsiL: bin width =", psiL_bin_width, "MPa ===\n")
cat("total bins:", nrow(mod_binned), "\n")
print(mod_binned %>% group_by(treatment, species) %>% summarise(bins = n(), .groups = "drop"))

# Plot: binned modelled Gc vs PsiL
p_binned_mod_gc_psiL <- ggplot(mod_binned, aes(x = psiL_bin, y = gc_mean, color = species)) +
  geom_errorbar(aes(ymin = gc_mean - gc_se, ymax = gc_mean + gc_se),
                width = 0.03, linewidth = 0.4, alpha = 0.5) +
  geom_point(size = 2.2, na.rm = TRUE) +
  geom_line(linewidth = 0.7, alpha = 0.7, na.rm = TRUE) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = cb_palette) +
  labs(
    title = "binned modelled canopy conductance vs leaf water potential (LPJ only)",
    subtitle = paste0("psiL binned at ", psiL_bin_width, " MPa intervals | points = bin-mean Gc | error bars = ±1 SE | line connects bins"),
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]~(m~s^{-1})),
    color = ""
  ) +
  base_theme +
  theme(legend.position = "none")

print(p_binned_mod_gc_psiL)

# Also create a single-panel version (control vs drought side-by-side)
p_binned_mod_gc_psiL_single <- ggplot(mod_binned, aes(x = psiL_bin, y = gc_mean, color = species)) +
  geom_errorbar(aes(ymin = gc_mean - gc_se, ymax = gc_mean + gc_se),
                width = 0.03, linewidth = 0.4, alpha = 0.5) +
  geom_point(size = 2.2, na.rm = TRUE) +
  geom_line(linewidth = 0.7, alpha = 0.7, na.rm = TRUE) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_color_manual(values = cb_palette) +
  labs(
    title = "binned modelled canopy conductance vs leaf water potential (LPJ only)",
    subtitle = paste0("psiL binned at ", psiL_bin_width, " MPa intervals | points = bin-mean Gc | error bars = ±1 SE | line connects bins"),
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]~(m~s^{-1})),
    color = "species"
  ) +
  base_theme

print(p_binned_mod_gc_psiL_single)

# ------------------------------------------------------------------------------
# 8b. BINNED & STANDARDIZED per species (y-axis 0-1 within each species)
# ------------------------------------------------------------------------------

# Standardize binned Gc to 0-1 range separately for each species
# (pooling both treatments so control vs drought curves are comparable per species)
mod_binned_std <- mod_binned %>%
  group_by(species) %>%
  mutate(
    gc_bin_min = min(gc_mean, na.rm = TRUE),
    gc_bin_max = max(gc_mean, na.rm = TRUE),
    gc_std     = (gc_mean - gc_bin_min) / (gc_bin_max - gc_bin_min),
    gc_std     = pmin(pmax(gc_std, 0), 1)
  ) %>%
  ungroup()

# Print standardization ranges per species
cat("\n=== per-species standardization ranges (binned) ===\n")
mod_binned_std %>%
  group_by(species) %>%
  summarise(gc_min = min(gc_mean), gc_max = max(gc_mean), .groups = "drop") %>%
  print()

# Faceted grid: rows = treatment, cols = species
p_binned_std <- ggplot(mod_binned_std, aes(x = psiL_bin, y = gc_std, color = species)) +
  geom_errorbar(aes(
    ymin = pmax(gc_std - gc_se / (gc_bin_max - gc_bin_min), 0),
    ymax = pmin(gc_std + gc_se / (gc_bin_max - gc_bin_min), 1)
  ), width = 0.03, linewidth = 0.4, alpha = 0.5) +
  geom_point(size = 2.2, na.rm = TRUE) +
  geom_line(linewidth = 0.7, alpha = 0.7, na.rm = TRUE) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = cb_palette) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "standardized binned Gc vs PsiL (LPJ only, per-species standardization)",
    subtitle = paste0("psiL binned at ", psiL_bin_width, " MPa | Gc standardized 0-1 within each species | error bars = ±1 SE (scaled)"),
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]~standardized~(0-1~per~species)),
    color = ""
  ) +
  base_theme +
  theme(legend.position = "none")

print(p_binned_std)

# Single-panel version: control vs drought side-by-side
p_binned_std_single <- ggplot(mod_binned_std, aes(x = psiL_bin, y = gc_std, color = species)) +
  geom_errorbar(aes(
    ymin = pmax(gc_std - gc_se / (gc_bin_max - gc_bin_min), 0),
    ymax = pmin(gc_std + gc_se / (gc_bin_max - gc_bin_min), 1)
  ), width = 0.03, linewidth = 0.4, alpha = 0.5) +
  geom_point(size = 2.2, na.rm = TRUE) +
  geom_line(linewidth = 0.7, alpha = 0.7, na.rm = TRUE) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "standardized binned Gc vs PsiL (LPJ only, per-species standardization)",
    subtitle = paste0("psiL binned at ", psiL_bin_width, " MPa | Gc standardized 0-1 within each species | error bars = ±1 SE (scaled)"),
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]~standardized~(0-1~per~species)),
    color = "species"
  ) +
  base_theme

print(p_binned_std_single)

# ------------------------------------------------------------------------------
# 8c. BINNED LPJ + RAW OBSERVED: Gc/Gcmax vs PsiL (full daily, standardized)
# ------------------------------------------------------------------------------

# Bin the modelled standardized data by psiL per species & treatment
lpj_binned_rel <- data_full_model_std %>%
  filter(gc_rel_mod > 0) %>%
  mutate(psiL_bin = round(psiL / psiL_bin_width) * psiL_bin_width) %>%
  group_by(treatment, species, psiL_bin) %>%
  summarise(
    gc_rel_mean = mean(gc_rel_mod, na.rm = TRUE),
    gc_rel_se   = sd(gc_rel_mod, na.rm = TRUE) / sqrt(n()),
    psiL_mean   = mean(psiL, na.rm = TRUE),
    n           = n(),
    .groups     = "drop"
  ) %>%
  filter(n >= 3)

cat("\n=== binned LPJ standardized Gc/Gcmax: bin width =", psiL_bin_width, "MPa ===\n")
cat("total bins:", nrow(lpj_binned_rel), "\n")

# --- Version A: Single-panel (matching p_gc_rel_full_single_daily style) ---
p_binned_rel_single <- ggplot() +
  # Observed midday (raw scatter)
  geom_point(
    data = data_obs_full_std %>% filter(psiL_label == "midday"),
    aes(x = psiL, y = gc_rel, color = species, shape = "obs midday"),
    alpha = 0.35, size = pt_size, na.rm = TRUE
  ) +
  # Observed predawn (raw scatter)
  geom_point(
    data = data_obs_full_std %>% filter(psiL_label == "predawn"),
    aes(x = psiL, y = gc_rel, color = species, shape = "obs predawn"),
    alpha = 0.35, size = pt_size, na.rm = TRUE
  ) +
  # LPJ binned: error bars + points + connecting line
  geom_errorbar(
    data = lpj_binned_rel,
    aes(x = psiL_bin, y = gc_rel_mean, color = species,
        ymin = pmax(gc_rel_mean - gc_rel_se, 0),
        ymax = pmin(gc_rel_mean + gc_rel_se, 1)),
    width = 0.04, linewidth = 0.4, alpha = 0.5
  ) +
  geom_point(
    data = lpj_binned_rel,
    aes(x = psiL_bin, y = gc_rel_mean, color = species, shape = "lpj binned"),
    size = 2.5, na.rm = TRUE
  ) +
  geom_line(
    data = lpj_binned_rel,
    aes(x = psiL_bin, y = gc_rel_mean, color = species),
    linewidth = 0.8, alpha = 0.8, na.rm = TRUE
  ) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(
    name = "",
    values = c("obs midday" = 2, "obs predawn" = 3, "lpj binned" = 16)
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "standardized Gc vs PsiL — binned LPJ + raw observed (full daily)",
    subtitle = paste0("LPJ psiL binned at ", psiL_bin_width, " MPa | binned points = mean Gc/Gcmax | error bars = ±1 SE"),
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]/G[cmax]),
    color = "species"
  ) +
  base_theme

print(p_binned_rel_single)

# --- Version B: Faceted grid (rows = treatment, cols = species) ---
p_binned_rel_grid <- ggplot() +
  geom_point(
    data = data_obs_full_std %>% filter(psiL_label == "midday"),
    aes(x = psiL, y = gc_rel, color = species, shape = "obs midday"),
    alpha = 0.3, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = data_obs_full_std %>% filter(psiL_label == "predawn"),
    aes(x = psiL, y = gc_rel, color = species, shape = "obs predawn"),
    alpha = 0.3, size = pt_size, na.rm = TRUE
  ) +
  geom_errorbar(
    data = lpj_binned_rel,
    aes(x = psiL_bin, y = gc_rel_mean, color = species,
        ymin = pmax(gc_rel_mean - gc_rel_se, 0),
        ymax = pmin(gc_rel_mean + gc_rel_se, 1)),
    width = 0.04, linewidth = 0.4, alpha = 0.5
  ) +
  geom_point(
    data = lpj_binned_rel,
    aes(x = psiL_bin, y = gc_rel_mean, color = species, shape = "lpj binned"),
    size = 2.5, na.rm = TRUE
  ) +
  geom_line(
    data = lpj_binned_rel,
    aes(x = psiL_bin, y = gc_rel_mean, color = species),
    linewidth = 0.8, alpha = 0.8, na.rm = TRUE
  ) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(
    name = "",
    values = c("obs midday" = 2, "obs predawn" = 3, "lpj binned" = 16)
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "standardized Gc vs PsiL — binned LPJ + raw observed (full daily)",
    subtitle = paste0("LPJ psiL binned at ", psiL_bin_width, " MPa | binned points = mean Gc/Gcmax | error bars = ±1 SE"),
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]/G[cmax]),
    color = ""
  ) +
  base_theme +
  theme(legend.position = "none")

print(p_binned_rel_grid)

# ==============================================================================
# 9. CURVE FITTING & SLOPE CALCULATION AT 0.5 Gc/Gcmax
# ==============================================================================

fit_exp_slope <- function(df) {
  if (nrow(df) < 5) return(data.frame(slope_at_05 = NA))
  mod <- lm(log(gc_rel) ~ psiL, data = df)
  b_parameter <- coef(mod)[2]
  data.frame(slope_at_05 = 0.5 * b_parameter)
}

stream_common_md <- combined_std_lpj_obs %>% filter(source == "obs", psiL_label == "midday", gc_rel > 0) %>% select(treatment, species, psiL, gc_rel) %>%
  mutate(data_type = "obs midday", timeline = "common time")
stream_common_pd <- combined_std_lpj_obs %>% filter(source == "obs", psiL_label == "predawn", gc_rel > 0) %>% select(treatment, species, psiL, gc_rel) %>%
  mutate(data_type = "obs predawn", timeline = "common time")
stream_common_mod<- combined_std_lpj_obs %>% filter(source == "sim", gc_rel > 0) %>% select(treatment, species, psiL, gc_rel) %>%
  mutate(data_type = "lpj simulated", timeline = "common time")

stream_full_md   <- data_obs_full_std %>% filter(psiL_label == "midday", gc_rel > 0) %>% select(treatment, species, psiL, gc_rel) %>%
  mutate(data_type = "obs midday", timeline = "full time range")
stream_full_pd   <- data_obs_full_std %>% filter(psiL_label == "predawn", gc_rel > 0) %>% select(treatment, species, psiL, gc_rel) %>%
  mutate(data_type = "obs predawn", timeline = "full time range")
stream_full_mth  <- data_full_model_monthly_std %>% filter(gc_rel_mod > 0) %>% select(treatment, species, psiL, gc_rel = gc_rel_mod)  %>%
  mutate(data_type = "lpj monthly mean", timeline = "full time range")
stream_full_dly  <- data_full_model_std %>% filter(gc_rel_mod > 0) %>% select(treatment, species, psiL, gc_rel = gc_rel_mod)          %>%
  mutate(data_type = "lpj daily simulated", timeline = "full time range")

master_long_dataset <- bind_rows(
  stream_common_md, stream_common_pd, stream_common_mod,
  stream_full_md, stream_full_pd, stream_full_mth, stream_full_dly
) %>% drop_na(psiL, gc_rel)

calculated_slopes <- master_long_dataset %>%
  group_by(timeline, treatment, species, data_type) %>%
  do(fit_exp_slope(.)) %>%
  ungroup() %>%
  mutate(
    species = factor(species, levels = species_order),
    data_type = factor(data_type, levels = c("obs midday", "obs predawn", "lpj simulated", "lpj monthly mean", "lpj daily simulated"))
  )

print("--- calculated exponential curve slopes at gc/gcmax = 0.5 ---")
print(as.data.frame(calculated_slopes))

p_slopes_bar <- ggplot(calculated_slopes, aes(x = data_type, y = slope_at_05, fill = species)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2, alpha = 0.9) +
  facet_grid(timeline ~ treatment) +
  scale_fill_manual(values = cb_palette) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
    strip.text  = element_text(size = 13, face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "calculated exponential curve slope at gc/gcmax = 0.5",
    subtitle = "derived using dy/dx = 0.5 * b parameter | sorted by treatment and configuration",
    x = "",
    y = "instantaneous slope (conductance sensitivity)",
    fill = "species"
  )
print(p_slopes_bar)

# ==============================================================================
# 10. PRINT & EXPORT TARGET DIRECTORY
# ==============================================================================

out_dir <- "Figures/lpj_guess_stem_storage/Gc_PsiL"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "compare_Gc_vs_PsiL_common_time.png"), p_gc_psi_common, width = 13, height = 9, dpi = 300)
ggsave(file.path(out_dir, "compare_Gc_vs_PsiL_full.png"), p_gc_psi_full, width = 13, height = 9, dpi = 300)
ggsave(file.path(out_dir, "compare_Gc_rel_vs_PsiL_common_time.png"), p_gc_rel_common, width = 13, height = 9, dpi = 300)
ggsave(file.path(out_dir, "compare_Gc_rel_vs_PsiL_full.png"), p_gc_rel_full, width = 13, height = 9, dpi = 300)
ggsave(file.path(out_dir, "single_panel_Gc_vs_PsiL_common.png"), p_gc_psi_common_single, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "single_panel_Gc_rel_vs_PsiL_common.png"), p_gc_rel_common_single, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "single_panel_Gc_vs_PsiL_full_monthly.png"), p_gc_psi_full_single_monthly_mean, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "single_panel_Gc_rel_vs_PsiL_full_monthly.png"), p_gc_rel_full_single_monthly_mean, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "single_panel_Gc_vs_PsiL_full_daily.png"), p_gc_psi_full_single_daily, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "single_panel_Gc_rel_vs_PsiL_full_daily.png"), p_gc_rel_full_single_daily, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "theil_sen_absolute.png"), p_ts_abs, width = 10, height = 8, dpi = 300)
ggsave(file.path(out_dir, "theil_sen_relative.png"), p_ts_rel, width = 10, height = 8, dpi = 300)
ggsave(file.path(out_dir, "calculated_slopes_comparison_bar.png"), p_slopes_bar, width = 12, height = 9, dpi = 300)
ggsave(file.path(out_dir, "binned_modelled_Gc_vs_PsiL.png"), p_binned_mod_gc_psiL, width = 13, height = 9, dpi = 300)
ggsave(file.path(out_dir, "binned_modelled_Gc_vs_PsiL_single.png"), p_binned_mod_gc_psiL_single, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "binned_std_Gc_vs_PsiL.png"), p_binned_std, width = 13, height = 9, dpi = 300)
ggsave(file.path(out_dir, "binned_std_Gc_vs_PsiL_single.png"), p_binned_std_single, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "binned_rel_Gc_vs_PsiL_single.png"), p_binned_rel_single, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "binned_rel_Gc_vs_PsiL_grid.png"), p_binned_rel_grid, width = 13, height = 9, dpi = 300)


# ==============================================================================
# 11. ADD 50% Gc/Gcmax LOSS AND FITTED CRUVE FOR OBSERVED DATA
# ==============================================================================

library(minpack.lm)  # robust non-linear least squares
library(patchwork)    # combining panel matrices

# ------------------------------------------------------------------------------
# A. PREPARE THE OBSERVED LONG DATASTREAM FOR FITTING ONLY
# ------------------------------------------------------------------------------
fit_stream_md <- combined_std_lpj_obs %>%
  filter(source == "obs", psiL_label == "midday") %>%
  select(treatment, species, psiL, gc_rel) %>%
  mutate(data_type = "obs midday")

fit_stream_pd <- combined_std_lpj_obs %>%
  filter(source == "obs", psiL_label == "predawn") %>%
  select(treatment, species, psiL, gc_rel) %>%
  mutate(data_type = "obs predawn")

obs_relative_data <- bind_rows(fit_stream_md, fit_stream_pd) %>%
  drop_na(psiL, gc_rel) %>%
  filter(gc_rel > 0 & gc_rel <= 1.1) %>% 
  mutate(X = -psiL) 

# ------------------------------------------------------------------------------
# D. ROBUST CURVE FITTING FUNCTION FOR RELATIVE DATA 
#    (Standard exponential OR Sigmoid for drought predawn)
# ------------------------------------------------------------------------------
fit_obs_rel_metrics <- function(sub_df, treatment_val, data_type_val) {
  fail_output <- data.frame(a = NA, b = NA, c = NA, psi50 = NA, slope50 = NA)
  if (nrow(sub_df) < 6) return(fail_output)
  
  # Special case: drought treatment + predawn -> use sigmoid function
  if (treatment_val == "drought" && data_type_val == "obs predawn") {
    tryCatch({
      # Sigmoid function: gc_rel = a + b / (1 + exp(c * (psiL - d)))
      # This allows initial plateau/rise before decline
      fit <- nlsLM(
        gc_rel ~ a + b / (1 + exp(c * (psiL - d))),
        data = sub_df,
        start = list(a = 0.05, b = 0.9, c = -2, d = -1),
        lower = c(a = -0.2, b = 0.3, c = -10, d = -3),
        upper = c(a = 0.5, b = 1.5, c = 0, d = 0),
        control = nls.lm.control(maxiter = 200)
      )
      
      cc <- coef(fit)
      a_val <- cc[["a"]]
      b_val <- cc[["b"]]
      c_val <- cc[["c"]]
      d_val <- cc[["d"]]
      
      # Calculate Psi50 (where gc_rel = 0.5)
      if (0.5 > a_val && 0.5 < (a_val + b_val)) {
        psi50 <- d_val + log((b_val / (0.5 - a_val)) - 1) / c_val
        # Slope at Psi50: derivative of sigmoid at midpoint
        slope50 <- -b_val * c_val * exp(c_val * (psi50 - d_val)) / 
          (1 + exp(c_val * (psi50 - d_val)))^2
        slope50 <- abs(slope50)
      } else {
        psi50 <- NA
        slope50 <- NA
      }
      
      return(data.frame(a = a_val, b = b_val, c = c_val, psi50 = psi50, slope50 = slope50))
      
    }, error = function(e) {
      # Fall back to exponential if sigmoid fails
      return(fit_obs_rel_exp(sub_df))
    })
  } else {
    # All other cases: standard exponential decay
    return(fit_obs_rel_exp(sub_df))
  }
}

# Standard exponential fitting function
fit_obs_rel_exp <- function(sub_df) {
  tryCatch({
    fit <- nlsLM(
      gc_rel ~ a + b * exp(-c * X), 
      data = sub_df %>% mutate(X = -psiL),
      start = list(a = 0.01, b = 0.95, c = 0.5),
      lower = c(a = -0.2,  b = 0.3,  c = 0.01),
      upper = c(a = 0.4,   b = 1.5,  c = 15.0),
      control = nls.lm.control(maxiter = 120)
    )
    
    cc <- coef(fit)
    a_val <- cc[["a"]]
    b_val <- cc[["b"]]
    c_val <- cc[["c"]]
    
    if ((0.5 - a_val) > 0 && b_val > 0) {
      x50 <- -log((0.5 - a_val) / b_val) / c_val
      psi50 <- -x50 
      slope50 <- c_val * (0.5 - a_val)
    } else {
      psi50 <- NA
      slope50 <- NA
    }
    
    return(data.frame(a = a_val, b = b_val, c = c_val, psi50 = psi50, slope50 = slope50))
    
  }, error = function(e) {
    return(data.frame(a = NA, b = NA, c = NA, psi50 = NA, slope50 = NA))
  })
}

# Apply fitting with treatment and data_type info
obs_rel_parameters <- obs_relative_data %>%
  group_by(treatment, species, data_type) %>%
  do(fit_obs_rel_metrics(., .$treatment[1], .$data_type[1])) %>%
  ungroup()

# ------------------------------------------------------------------------------
# E. GENERATE SMOOTH PREDICTION LINES FOR THE OBSERVED CURVES
# ------------------------------------------------------------------------------
obs_rel_predictions <- obs_rel_parameters %>%
  drop_na(psi50) %>%
  group_by(treatment, species, data_type) %>%
  do({
    sub_df <- obs_relative_data %>% 
      filter(treatment == .$treatment & species == .$species & data_type == .$data_type)
    
    psi_range <- seq(min(sub_df$psiL, na.rm = TRUE), max(sub_df$psiL, na.rm = TRUE), length.out = 100)
    
    # Check if this is the drought+predawn case (sigmoid)
    if (.$treatment[1] == "drought" && .$data_type[1] == "obs predawn" && !is.na(.$c[1]) && .$c[1] < 0) {
      # Sigmoid prediction
      gc_rel_pred <- .$a[1] + .$b[1] / (1 + exp(.$c[1] * (psi_range - .$psi50[1])))
    } else {
      # Exponential decay prediction
      gc_rel_pred <- .$a[1] + .$b[1] * exp(-.$c[1] * (-psi_range))
    }
    
    data.frame(
      psiL = psi_range,
      gc_rel_pred = gc_rel_pred
    )
  }) %>%
  ungroup()

# ------------------------------------------------------------------------------
# F. PANEL A: 4-PANEL SCATTER & CURVES MATRIX (Left Side)
# ------------------------------------------------------------------------------
obs_points_separated <- combined_std_lpj_obs %>%
  filter(source == "obs") %>%
  select(treatment, species, psiL, gc_rel, psiL_label) %>%
  mutate(
    data_type = paste0("obs ", psiL_label),
    point_type = paste0("obs ", psiL_label)
  )

# Use RAW LPJ values (not binned) for background points
raw_lpj_points <- combined_std_lpj_obs %>%
  filter(source == "sim") %>%
  select(treatment, species, psiL, gc_rel) %>%
  mutate(data_type = "obs midday", point_type = "lpj raw") %>%
  bind_rows(
    combined_std_lpj_obs %>%
      filter(source == "sim") %>%
      select(treatment, species, psiL, gc_rel) %>%
      mutate(data_type = "obs predawn", point_type = "lpj raw")
  )

p_curves <- ggplot() +
  # 1. Plot separated observation points
  geom_point(
    data = obs_points_separated,
    aes(x = psiL, y = gc_rel, color = species, shape = point_type),
    alpha = 0.3, size = pt_size, na.rm = TRUE
  ) +
  # 2. Contextual background RAW LPJ model points (not binned)
  geom_point(
    data = raw_lpj_points,
    aes(x = psiL, y = gc_rel, color = species, shape = point_type),
    alpha = 0.4, size = pt_size - 0.5, na.rm = TRUE
  ) +
  # 3. Fitted curves matching each specific matrix node
  geom_line(
    data = obs_rel_predictions,
    aes(x = psiL, y = gc_rel_pred, color = species),
    linewidth = 1.0, na.rm = TRUE
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  
  facet_grid(rows = vars(treatment), cols = vars(data_type)) + 
  
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "Data Points", 
                     values = c("obs midday" = 2, "obs predawn" = 3, "lpj raw" = 1)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]/G[cmax]),
    tag = "A"
  ) + 
  base_theme +
  theme(legend.position = "bottom")

# ------------------------------------------------------------------------------
# G. PANEL B: 4-PANEL SLOPE MATRIX (Right Side)
# ------------------------------------------------------------------------------
p_slopes <- ggplot(
  data = obs_rel_parameters %>% drop_na(slope50), 
  aes(x = species, y = slope50, fill = species)
) +
  geom_bar(stat = "identity", color = "black", width = 0.6, show.legend = FALSE) +
  
  facet_grid(rows = vars(treatment), cols = vars(data_type)) + 
  
  scale_fill_manual(values = cb_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Species",
    y = expression("Absolute Slope at " ~ Psi[50] ~ " (MPa"^{-1} ~ ")"),
    tag = "B"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    panel.grid.major.x = element_blank(),
    strip.text.y = element_blank()
  )

# ------------------------------------------------------------------------------
# H. COMBINE MATRICES AND EXPORT
# ------------------------------------------------------------------------------
composite_8panel_plot <- p_curves + p_slopes + 
  plot_layout(widths = c(1, 1), guides = "collect") & 
  theme(legend.position = "bottom")

composite_8panel_plot <- composite_8panel_plot +
  plot_annotation(
    title = "stomatal sensitivity observed with simulated from LPJ-GUESS-HYD (common)",
    subtitle = "drought predawn uses sigmoid function for spruce and pine; all others use exponential curve",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, color = "grey30", hjust = 0.5)
    )
  )

print(composite_8panel_plot)

ggsave(
  file.path(out_dir, "compare_stomamtal_sensitivity_common_time.png"),
  composite_8panel_plot,
  width = 16, height = 10, dpi = 300
)

# ==============================================================================
# 11b. FULL TIME REFERENCE - FIT CURVES + GENERATE 8-PANEL COMPOSITE FIGURE
# ==============================================================================

# ------------------------------------------------------------------------------
# A. USE FULL TIME STANDARDIZED DATA (NOT COMMON TIME)
# ------------------------------------------------------------------------------
# Use data_obs_full_std (long format) and data_full_model_std (wide format)

# ------------------------------------------------------------------------------
# B. PREPARE RAW LPJ DATA (FULL TIME) FOR BACKGROUND CONTEXT
# ------------------------------------------------------------------------------
raw_lpj_points_full <- data_full_model_std %>%
  select(treatment, species, psiL, gc_rel = gc_rel_mod) %>%
  mutate(data_type = "obs midday") %>%
  bind_rows(
    data_full_model_std %>%
      select(treatment, species, psiL, gc_rel = gc_rel_mod) %>%
      mutate(data_type = "obs predawn")
  )

# ------------------------------------------------------------------------------
# C. PREPARE THE OBSERVED LONG DATASTREAM FOR FITTING (FULL TIME)
# ------------------------------------------------------------------------------
fit_stream_md_full <- data_obs_full_std %>%
  filter(psiL_label == "midday") %>%
  select(treatment, species, psiL, gc_rel) %>%
  mutate(data_type = "obs midday")

fit_stream_pd_full <- data_obs_full_std %>%
  filter(psiL_label == "predawn") %>%
  select(treatment, species, psiL, gc_rel) %>%
  mutate(data_type = "obs predawn")

obs_relative_data_full <- bind_rows(fit_stream_md_full, fit_stream_pd_full) %>%
  drop_na(psiL, gc_rel) %>%
  filter(gc_rel > 0 & gc_rel <= 1.1) %>% 
  mutate(X = -psiL)

# ------------------------------------------------------------------------------
# D. ROBUST CURVE FITTING FUNCTION (SIGMOID FOR DROUGHT PREDAWN, EXPONENTIAL FOR OTHERS)
# ------------------------------------------------------------------------------
fit_obs_rel_metrics_full <- function(sub_df, treatment_val, data_type_val) {
  fail_output <- data.frame(a = NA, b = NA, c = NA, psi50 = NA, slope50 = NA)
  if (nrow(sub_df) < 6) return(fail_output)
  
  # Special case: drought treatment + predawn -> use sigmoid function
  if (treatment_val == "drought" && data_type_val == "obs predawn") {
    tryCatch({
      fit <- nlsLM(
        gc_rel ~ a + b / (1 + exp(c * (psiL - d))),
        data = sub_df,
        start = list(a = 0.05, b = 0.9, c = -2, d = -1),
        lower = c(a = -0.2, b = 0.3, c = -10, d = -3),
        upper = c(a = 0.5, b = 1.5, c = 0, d = 0),
        control = nls.lm.control(maxiter = 200)
      )
      
      cc <- coef(fit)
      a_val <- cc[["a"]]
      b_val <- cc[["b"]]
      c_val <- cc[["c"]]
      d_val <- cc[["d"]]
      
      if (0.5 > a_val && 0.5 < (a_val + b_val)) {
        psi50 <- d_val + log((b_val / (0.5 - a_val)) - 1) / c_val
        slope50 <- -b_val * c_val * exp(c_val * (psi50 - d_val)) / 
          (1 + exp(c_val * (psi50 - d_val)))^2
        slope50 <- abs(slope50)
      } else {
        psi50 <- NA
        slope50 <- NA
      }
      
      return(data.frame(a = a_val, b = b_val, c = c_val, psi50 = psi50, slope50 = slope50))
      
    }, error = function(e) {
      return(fit_obs_rel_exp_full(sub_df))
    })
  } else {
    return(fit_obs_rel_exp_full(sub_df))
  }
}

# Exponential fitting function (full time)
fit_obs_rel_exp_full <- function(sub_df) {
  tryCatch({
    fit <- nlsLM(
      gc_rel ~ a + b * exp(-c * X), 
      data = sub_df %>% mutate(X = -psiL),
      start = list(a = 0.01, b = 0.95, c = 0.5),
      lower = c(a = -0.2, b = 0.3, c = 0.01),
      upper = c(a = 0.4, b = 1.5, c = 15.0),
      control = nls.lm.control(maxiter = 120)
    )
    
    cc <- coef(fit)
    a_val <- cc[["a"]]
    b_val <- cc[["b"]]
    c_val <- cc[["c"]]
    
    if ((0.5 - a_val) > 0 && b_val > 0) {
      x50 <- -log((0.5 - a_val) / b_val) / c_val
      psi50 <- -x50 
      slope50 <- c_val * (0.5 - a_val)
    } else {
      psi50 <- NA
      slope50 <- NA
    }
    
    return(data.frame(a = a_val, b = b_val, c = c_val, psi50 = psi50, slope50 = slope50))
    
  }, error = function(e) {
    return(data.frame(a = NA, b = NA, c = NA, psi50 = NA, slope50 = NA))
  })
}

# Apply fitting
obs_rel_parameters_full <- obs_relative_data_full %>%
  group_by(treatment, species, data_type) %>%
  do(fit_obs_rel_metrics_full(., .$treatment[1], .$data_type[1])) %>%
  ungroup()

# ------------------------------------------------------------------------------
# E. GENERATE SMOOTH PREDICTION LINES (FULL TIME)
# ------------------------------------------------------------------------------
obs_rel_predictions_full <- obs_rel_parameters_full %>%
  drop_na(psi50) %>%
  group_by(treatment, species, data_type) %>%
  do({
    sub_df <- obs_relative_data_full %>% 
      filter(treatment == .$treatment & species == .$species & data_type == .$data_type)
    
    psi_range <- seq(min(sub_df$psiL, na.rm = TRUE), 
                     max(sub_df$psiL, na.rm = TRUE), 
                     length.out = 100)
    
    # Check for drought+predawn sigmoid case
    if (.$treatment[1] == "drought" && .$data_type[1] == "obs predawn" && !is.na(.$c[1]) && .$c[1] < 0) {
      gc_rel_pred <- .$a[1] + .$b[1] / (1 + exp(.$c[1] * (psi_range - .$psi50[1])))
    } else {
      gc_rel_pred <- .$a[1] + .$b[1] * exp(-.$c[1] * (-psi_range))
    }
    
    data.frame(psiL = psi_range, gc_rel_pred = gc_rel_pred)
  }) %>%
  ungroup()

# ------------------------------------------------------------------------------
# F. PREPARE OBSERVATION POINTS (FULL TIME)
# ------------------------------------------------------------------------------
obs_points_separated_full <- data_obs_full_std %>%
  select(treatment, species, psiL, gc_rel, psiL_label) %>%
  mutate(
    data_type = paste0("obs ", psiL_label),
    point_type = paste0("obs ", psiL_label)
  )

# ------------------------------------------------------------------------------
# G. PANEL A: 4-PANEL SCATTER & CURVES MATRIX (FULL TIME)
# ------------------------------------------------------------------------------
p_curves_full <- ggplot() +
  geom_point(
    data = obs_points_separated_full,
    aes(x = psiL, y = gc_rel, color = species, shape = point_type),
    alpha = 0.3, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = raw_lpj_points_full,
    aes(x = psiL, y = gc_rel, color = species, shape = "lpj raw"),
    alpha = 0.4, size = pt_size - 0.5, na.rm = TRUE
  ) +
  geom_line(
    data = obs_rel_predictions_full,
    aes(x = psiL, y = gc_rel_pred, color = species),
    linewidth = 1.0, na.rm = TRUE
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  
  facet_grid(rows = vars(treatment), cols = vars(data_type)) + 
  
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "Data Points", 
                     values = c("obs midday" = 2, "obs predawn" = 3, "lpj raw" = 1)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]/G[cmax]),
    tag = "A",
    title = "FULL TIME REFERENCE (true min-max per species × treatment)"
  ) + 
  base_theme +
  theme(legend.position = "bottom")

# ------------------------------------------------------------------------------
# H. PANEL B: 4-PANEL SLOPE MATRIX (FULL TIME)
# ------------------------------------------------------------------------------
p_slopes_full <- ggplot(
  data = obs_rel_parameters_full %>% drop_na(slope50), 
  aes(x = species, y = slope50, fill = species)
) +
  geom_bar(stat = "identity", color = "black", width = 0.6, show.legend = FALSE) +
  
  facet_grid(rows = vars(treatment), cols = vars(data_type)) + 
  
  scale_fill_manual(values = cb_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Species",
    y = expression("Absolute Slope at " ~ Psi[50] ~ " (MPa"^{-1} ~ ")"),
    tag = "B"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    panel.grid.major.x = element_blank(),
    strip.text.y = element_blank()
  )

# ------------------------------------------------------------------------------
# I. COMBINE AND EXPORT FULL TIME COMPOSITE
# ------------------------------------------------------------------------------
composite_8panel_plot_full <- p_curves_full + p_slopes_full + 
  plot_layout(widths = c(1, 1), guides = "collect") & 
  theme(legend.position = "bottom")

composite_8panel_plot_full <- composite_8panel_plot_full +
  plot_annotation(
    title = "stomatal sensitivity observed with simulated from LPJ-GUESS-HYD (full time series)",
    subtitle = "drought predawn uses sigmoid function for spruce and pine; all others use exponential curve | Normalized using true min-max per species × treatment (full time series)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, color = "grey30", hjust = 0.5)
    )
  )

print(composite_8panel_plot_full)

# Save
ggsave(
  file.path(out_dir, "compare_stomamtal_sensitivity_full_time.png"),
  composite_8panel_plot_full,
  width = 16, height = 10, dpi = 300
)

# ------------------------------------------------------------------------------
# J. OPTIONAL: Print comparison of Psi50 values between common time and full time
# ------------------------------------------------------------------------------
print("=== FULL TIME REFERENCE - PSI50 VALUES ===")
print(obs_rel_parameters_full %>% 
        select(treatment, species, data_type, psi50, slope50) %>% 
        drop_na(psi50) %>%
        arrange(treatment, data_type, species))

print("=== COMMON TIME REFERENCE - PSI50 VALUES ===")
print(obs_rel_parameters %>% 
        select(treatment, species, data_type, psi50, slope50) %>% 
        drop_na(psi50) %>%
        arrange(treatment, data_type, species))

