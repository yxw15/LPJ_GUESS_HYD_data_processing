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

lpj_output_filter <- read.csv("lpj_guess/lpj_control_drought_Gc_psiL_psiX_psiS_stem_diameter_twd_stem_rwc_climate_filter.csv") %>%
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

# Combine datasets for Common-Time bounds
combined_data_lpj_obs <- mod_proc %>%
  inner_join(obs_filtered_climate, by = c("date", "species", "treatment")) 

# ==============================================================================
# 5. MINIMUM VALUE STANDARDIZATION & AGGREGATION
# ==============================================================================

# 1. Common Time Dataset Standardization (True minimum value baseline)
combined_std_lpj_obs <- combined_data_lpj_obs %>%
  group_by(species, treatment) %>%
  mutate(
    # Changed from 10% quantile to absolute minimum value
    gc_min_obs = min(gc_obs, na.rm = TRUE),
    gc_min     = min(gc, na.rm = TRUE),
    
    gc_max_obs = mean(gc_obs[gc_obs >= quantile(gc_obs, 0.90, na.rm = TRUE)], na.rm = TRUE),
    gc_max     = mean(gc[gc >= quantile(gc, 0.90, na.rm = TRUE)], na.rm = TRUE),
    
    gc_rel_obs = (gc_obs - gc_min_obs) / (gc_max_obs - gc_min_obs),
    gc_rel_mod = (gc - gc_min) / (gc_max - gc_min),
    
    gc_rel_obs = pmin(pmax(gc_rel_obs, 0), 1),
    gc_rel_mod = pmin(pmax(gc_rel_mod, 0), 1)
  ) %>% 
  ungroup()

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

# 2. Full Simulated Series Standardization (Daily true minimum value baseline)
data_full_model_std <- mod_proc %>%
  filter(date %in% climate_filter_dates, month(date) %in% c(6, 7, 8, 9)) %>%
  group_by(species, treatment) %>%
  mutate(
    # Changed from 10% quantile to absolute minimum value
    gc_min = min(gc, na.rm = TRUE),
    gc_max = mean(gc[gc >= quantile(gc, 0.90, na.rm = TRUE)], na.rm = TRUE),
    gc_rel_mod = (gc - gc_min) / (gc_max - gc_min),
    gc_rel_mod = pmin(pmax(gc_rel_mod, 0), 1)
  ) %>%
  ungroup() 

# 3. Full Observed Series Standardization (Daily true minimum value baseline)
data_obs_full_std <- obs_filtered_climate %>%
  group_by(species, treatment) %>%
  mutate(
    # Changed from 10% quantile to absolute minimum value
    gc_min_obs = min(gc_obs, na.rm = TRUE),
    gc_max_obs = mean(gc_obs[gc_obs >= quantile(gc_obs, 0.90, na.rm = TRUE)], na.rm = TRUE),
    gc_rel_obs = (gc_obs - gc_min_obs) / (gc_max_obs - gc_min_obs),
    gc_rel_obs = pmin(pmax(gc_rel_obs, 0), 1)
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
    # Changed from 10% quantile to absolute minimum value of monthly means
    gc_min     = min(gc, na.rm = TRUE),
    gc_max     = mean(gc[gc >= quantile(gc, 0.90, na.rm = TRUE)], na.rm = TRUE),
    gc_rel_mod = (gc - gc_min) / (gc_max - gc_min),
    gc_rel_mod = pmin(pmax(gc_rel_mod, 0), 1) 
  ) %>%
  ungroup()


# ==============================================================================
# 6. VISUALIZATIONS: FACETED GRID DESIGNS (Rows = Tx, Cols = Sp)
# ==============================================================================

# Figure 1: Absolute Gc vs Psi_L (Common Time)
p_gc_psi_common <- ggplot(combined_data_lpj_obs) +
  geom_point(aes(x = psiL_md, y = gc_obs), shape = 2, color = "black", alpha = 0.8, size = pt_size, na.rm = TRUE) +
  geom_point(aes(x = psiL_pd, y = gc_obs), color = "black", shape = 3, alpha = 0.8, size = pt_size, na.rm = TRUE) +
  geom_point(aes(x = psiL, y = gc, color = species), alpha = 1, size = pt_size + 0.5, na.rm = TRUE) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = cb_palette) +
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

# Figure 3: Relative Gc/Gcmax (Common Time Bounds strictly 0-1)
p_gc_rel_common <- ggplot(combined_std_lpj_obs) +
  geom_point(aes(x = psiL_md, y = gc_rel_obs), shape = 2, color = "black", alpha = 0.8, size = pt_size, na.rm = TRUE) +
  geom_point(aes(x = psiL_pd, y = gc_rel_obs), color = "black", shape = 3, alpha = 0.8, size = pt_size, na.rm = TRUE) +
  geom_point(aes(x = psiL, y = gc_rel_mod, color = species), alpha = 1, size = pt_size + 0.5, na.rm = TRUE) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = cb_palette) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "standardized canopy conductance vs leaf water potential (common time)",
    subtitle = "min-max (true min to 90% quantile) normalized: showing internal values\nmidday = open triangle | predawn = black + | lpj = colored",
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]/G[cmax]),
    color = ""
  ) + base_theme

print(p_gc_rel_common)

# Figure 4: Relative Gc/Gcmax (Full Series - Independently Standardized Monthly Means)
p_gc_rel_full <- ggplot() +
  geom_point(
    data = data_obs_full_std %>% filter(gc_rel_obs > 0, gc_rel_obs < 1), 
    aes(x = psiL_md, y = gc_rel_obs), 
    shape = 2, color = "black", alpha = 0.8, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = data_obs_full_std %>% filter(gc_rel_obs > 0, gc_rel_obs < 1), 
    aes(x = psiL_pd, y = gc_rel_obs), 
    color = "black", shape = 3, alpha = 0.8, size = pt_size, na.rm = TRUE
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
    subtitle = "true min to 90% quantile normalization structure | lpj = colored monthly mean",
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

# Single Panel 2: Standardized Values (Common Time)
p_gc_rel_common_single <- ggplot() +
  geom_point(
    data = combined_std_lpj_obs,
    aes(x = psiL_md, y = gc_rel_obs, color = species, shape = "obs midday"), alpha = 0.6, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = combined_std_lpj_obs,
    aes(x = psiL_pd, y = gc_rel_obs, color = species, shape = "obs predawn"), alpha = 0.6, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = combined_std_lpj_obs,
    aes(x = psiL, y = gc_rel_mod, color = species, shape = "lpj simulated"), alpha = 1, size = pt_size + 0.5, na.rm = TRUE
  ) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "", values = c("obs midday" = 2, "obs predawn" = 3, "lpj simulated" = 16)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "standardized canopy conductance vs leaf water potential (common time)",
    subtitle = tolower(paste0(climate_txt, " (true min to 90% standardized bounds)\nmidday = open triangle | predawn = + | lpj-guess-hyd = ●")),
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]/G[cmax]),
    color = "species"
  ) + base_theme

print(p_gc_rel_common_single)

# ==============================================================================
# 7b. ADDITIONAL SINGLE PANEL DESIGNS (Full Monthly Mean vs Full Daily)
# ==============================================================================

# Single Panel 3: Absolute Values (Full Monthly Mean)
p_gc_psi_full_single_monthly_mean <- ggplot() +
  geom_point(data = obs_filtered_climate, aes(x = psiL_md, y = gc_obs, color = species, shape = "obs midday"), alpha = 0.5, size = pt_size, na.rm = TRUE) +
  geom_point(data = obs_filtered_climate, aes(x = psiL_pd, y = gc_obs, color = species, shape = "obs predawn"), alpha = 0.5, size = pt_size, na.rm = TRUE) +
  geom_point(data = mod_proc_monthly, aes(x = psiL, y = gc, color = species, shape = "lpj monthly mean"), alpha = 0.8, size = pt_size + 0.5, na.rm = TRUE) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_color_manual(values = cb_palette) +
  ylim(0, 12) +
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
    data = data_obs_full_std,
    aes(x = psiL_md, y = gc_rel_obs, color = species, shape = "obs midday"), alpha = 0.5, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = data_obs_full_std,
    aes(x = psiL_pd, y = gc_rel_obs, color = species, shape = "obs predawn"), alpha = 0.5, size = pt_size, na.rm = TRUE
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
    subtitle = "showing internal values using robust true min to 90% baseline thresholds",
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]/G[cmax]),
    color = "species"
  ) + base_theme

print(p_gc_rel_full_single_monthly_mean)

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

# Single Panel 6: Standardized Values (Full Daily Range)
p_gc_rel_full_single_daily <- ggplot() +
  geom_point(
    data = data_obs_full_std %>% filter(gc_rel_obs > 0, gc_rel_obs < 1),
    aes(x = psiL_md, y = gc_rel_obs, color = species, shape = "obs midday"), alpha = 0.5, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = data_obs_full_std %>% filter(gc_rel_obs > 0, gc_rel_obs < 1),
    aes(x = psiL_pd, y = gc_rel_obs, color = species, shape = "obs predawn"), alpha = 0.5, size = pt_size, na.rm = TRUE
  ) +
  geom_point(
    data = data_full_model_std %>% filter(gc_rel_mod > 0, gc_rel_mod < 1),
    aes(x = psiL, y = gc_rel_mod, color = species, shape = "lpj daily simulated"), alpha = 0.3, size = pt_size, na.rm = TRUE
  ) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "", values = c("obs midday" = 2, "obs predawn" = 3, "lpj daily simulated" = 16)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "standardized canopy conductance vs leaf water potential (full daily data)",
    subtitle = "showing internal values evaluated against robust true min to 90% baselines",
    x = expression(Psi["leaf"]~"(MPa)"),
    y = expression(G[c]/G[cmax]),
    color = "species"
  ) + base_theme
print(p_gc_rel_full_single_daily)

# ==============================================================================
# 9. CURVE FITTING & SLOPE CALCULATION AT 0.5 Gc/Gcmax
# ==============================================================================

fit_exp_slope <- function(df) {
  if (nrow(df) < 5) return(data.frame(slope_at_05 = NA))
  mod <- lm(log(gc_rel) ~ psiL, data = df)
  b_parameter <- coef(mod)[2]
  data.frame(slope_at_05 = 0.5 * b_parameter)
}

stream_common_md <- combined_std_lpj_obs %>% filter(gc_rel_obs > 0, gc_rel_obs < 1) %>% select(treatment, species, psiL = psiL_md, gc_rel = gc_rel_obs) %>% 
  mutate(data_type = "obs midday", timeline = "common time")
stream_common_pd <- combined_std_lpj_obs %>% filter(gc_rel_obs > 0, gc_rel_obs < 1) %>% select(treatment, species, psiL = psiL_pd, gc_rel = gc_rel_obs) %>% 
  mutate(data_type = "obs predawn", timeline = "common time")
stream_common_mod<- combined_std_lpj_obs %>% filter(gc_rel_mod > 0, gc_rel_mod < 1) %>% select(treatment, species, psiL, gc_rel = gc_rel_mod)          %>% 
  mutate(data_type = "lpj simulated", timeline = "common time")

stream_full_md   <- data_obs_full_std %>% filter(gc_rel_obs > 0, gc_rel_obs < 1) %>% select(treatment, species, psiL = psiL_md, gc_rel = gc_rel_obs)      %>% 
  mutate(data_type = "obs midday", timeline = "full time range")
stream_full_pd   <- data_obs_full_std %>% filter(gc_rel_obs > 0, gc_rel_obs < 1) %>% select(treatment, species, psiL = psiL_pd, gc_rel = gc_rel_obs)      %>% 
  mutate(data_type = "obs predawn", timeline = "full time range")
stream_full_mth  <- data_full_model_monthly_std %>% filter(gc_rel_mod > 0, gc_rel_mod < 1) %>% select(treatment, species, psiL, gc_rel = gc_rel_mod)  %>% 
  mutate(data_type = "lpj monthly mean", timeline = "full time range")
stream_full_dly  <- data_full_model_std %>% filter(gc_rel_mod > 0, gc_rel_mod < 1) %>% select(treatment, species, psiL, gc_rel = gc_rel_mod)          %>% 
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

out_dir <- "Figures/compare_Gc_PsiL_lpjtwd_hoelstein"
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
ggsave(file.path(out_dir, "calculated_slopes_comparison_bar.png"), p_slopes_bar, width = 12, height = 9, dpi = 300)

print(p_gc_rel_common_single)
print(p_gc_rel_full_single_daily)
