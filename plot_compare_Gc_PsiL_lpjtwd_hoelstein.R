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

lpj_output_filter <- read.csv("lpj_guess/lpj_guess_twd/lpj_control_drought_Gc_psiL_psiX_psiS_stem_diameter_twd_stem_rwc_climate_filter.csv") %>%
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

out_dir <- "Figures/lpj_guess_hyd_twd"
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
plot(p_gc_rel_full_single_monthly_mean)

# ==============================================================================
# 11. ADD 50% Gc/Gcmax LOSS AND FITTED CRUVE FOR OBSERVED DATA
# ==============================================================================

library(minpack.lm) # For robust non-linear least squares fitting
library(dplyr)
library(ggplot2)
library(patchwork) # For combining your panel matrices side-by-side

# ------------------------------------------------------------------------------
# A. PREPARE RAW LPJ DATA (NOT BINNED) FOR BACKGROUND CONTEXT
# ------------------------------------------------------------------------------
# Use raw LPJ modelled values - no binning
raw_lpj_points <- combined_std_lpj_obs %>%
  select(treatment, species, psiL, gc_rel_mod) %>%
  mutate(data_type = "obs midday")  # Will duplicate for both panels

# Duplicate for predawn panels
raw_lpj_points <- bind_rows(
  raw_lpj_points %>% mutate(data_type = "obs midday"),
  raw_lpj_points %>% mutate(data_type = "obs predawn")
)

# ------------------------------------------------------------------------------
# B. BINNING MODELLED DATA (0.05 MPa Bin Width for Optional Use)
# ------------------------------------------------------------------------------
combined_std_lpj_obs_binned <- combined_std_lpj_obs %>%
  mutate(psiL_bin = round(psiL / 0.05) * 0.05) %>%
  group_by(treatment, species, psiL_bin) %>%
  summarise(
    psiL       = median(psiL, na.rm = TRUE),
    gc_rel_mod = mean(gc_rel_mod, na.rm = TRUE),
    .groups    = "drop"
  )

# ------------------------------------------------------------------------------
# C. PREPARE THE OBSERVED LONG DATASTREAM FOR FITTING ONLY
# ------------------------------------------------------------------------------
fit_stream_md <- combined_std_lpj_obs %>%
  select(treatment, species, psiL = psiL_md, gc_rel = gc_rel_obs) %>%
  mutate(data_type = "obs midday")

fit_stream_pd <- combined_std_lpj_obs %>%
  select(treatment, species, psiL = psiL_pd, gc_rel = gc_rel_obs) %>%
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
obs_points_separated <- bind_rows(
  combined_std_lpj_obs %>% select(treatment, species, psiL = psiL_md, gc_rel_obs) %>% 
    mutate(data_type = "obs midday", point_type = "obs midday"),
  combined_std_lpj_obs %>% select(treatment, species, psiL = psiL_pd, gc_rel_obs) %>% 
    mutate(data_type = "obs predawn", point_type = "obs predawn")
)

# Use RAW LPJ values (not binned) for background points
raw_lpj_points <- bind_rows(
  combined_std_lpj_obs %>% 
    select(treatment, species, psiL, gc_rel_mod) %>% 
    mutate(data_type = "obs midday", point_type = "lpj raw"),
  combined_std_lpj_obs %>% 
    select(treatment, species, psiL, gc_rel_mod) %>% 
    mutate(data_type = "obs predawn", point_type = "lpj raw")
)

p_curves <- ggplot() +
  # 1. Plot separated observation points
  geom_point(
    data = obs_points_separated,
    aes(x = psiL, y = gc_rel_obs, color = species, shape = point_type), 
    alpha = 0.3, size = pt_size, na.rm = TRUE
  ) +
  # 2. Contextual background RAW LPJ model points (not binned)
  geom_point(
    data = raw_lpj_points,
    aes(x = psiL, y = gc_rel_mod, color = species, shape = point_type), 
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

library(minpack.lm)
library(dplyr)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
# A. USE FULL TIME STANDARDIZED DATA (NOT COMMON TIME)
# ------------------------------------------------------------------------------
# Use data_obs_full_std and data_full_model_std instead of combined_std_lpj_obs

# Prepare full time standardized data for fitting
full_time_data <- data_obs_full_std %>%
  left_join(
    data_full_model_std %>% select(date, species, treatment, psiL, gc_rel_mod),
    by = c("date", "species", "treatment")
  ) %>%
  rename(gc_rel_obs = gc_rel_obs, psiL_md = psiL_md, psiL_pd = psiL_pd)

# ------------------------------------------------------------------------------
# B. PREPARE RAW LPJ DATA (FULL TIME) FOR BACKGROUND CONTEXT
# ------------------------------------------------------------------------------
raw_lpj_points_full <- data_full_model_std %>%
  select(treatment, species, psiL, gc_rel_mod) %>%
  mutate(data_type = "obs midday") %>%
  bind_rows(
    data_full_model_std %>%
      select(treatment, species, psiL, gc_rel_mod) %>%
      mutate(data_type = "obs predawn")
  )

# ------------------------------------------------------------------------------
# C. PREPARE THE OBSERVED LONG DATASTREAM FOR FITTING (FULL TIME)
# ------------------------------------------------------------------------------
fit_stream_md_full <- data_obs_full_std %>%
  select(treatment, species, psiL = psiL_md, gc_rel = gc_rel_obs) %>%
  mutate(data_type = "obs midday")

fit_stream_pd_full <- data_obs_full_std %>%
  select(treatment, species, psiL = psiL_pd, gc_rel = gc_rel_obs) %>%
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
obs_points_separated_full <- bind_rows(
  data_obs_full_std %>% 
    select(treatment, species, psiL = psiL_md, gc_rel = gc_rel_obs) %>% 
    mutate(data_type = "obs midday", point_type = "obs midday"),
  data_obs_full_std %>% 
    select(treatment, species, psiL = psiL_pd, gc_rel = gc_rel_obs) %>% 
    mutate(data_type = "obs predawn", point_type = "obs predawn")
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
    aes(x = psiL, y = gc_rel_mod, color = species, shape = "lpj raw"), 
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
    title = "FULL TIME REFERENCE (true min to 90% quantile)"
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
    subtitle = "drought predawn uses sigmoid function for spruce and pine; all others use exponential curve | Normalized using true min to 90% quantile (full time series)",
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