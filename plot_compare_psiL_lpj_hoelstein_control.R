# ==========================================================================
# 1. SETUP & DATA LOADING
# ==========================================================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

# Aesthetic Configuration
species_levels <- c("Oak", "Beech", "Spruce", "Pine")
cb_palette     <- c(Oak = "#E69F00", Beech = "#0072B2", Spruce = "#009E73", Pine = "#F0E442")
custom_colors  <- c(cb_palette, "Observed" = "#999999") 

base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text       = element_text(color = "black", size = 14),
    legend.position   = "bottom",
    plot.title        = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title        = element_text(size = 16),
    axis.text.x       = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y       = element_text(angle = 0, hjust = 0.5, size = 12),
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_line(color = "grey92", linewidth = 0.25),
    strip.text        = element_text(size = 14, face = "bold")
  )

# Data Ingestion
lpj_base     <- read.csv("lpj_guess/lpj_Gc_psiL_psiX_psiS.csv") %>% mutate(date = as.Date(date))
psiL_control <- read.csv("SCCII/psiL_hoelstein_control.csv") %>% mutate(date = as.Date(date))
sap_control  <- read.csv("SCCII/sap_hoelstein_control.csv") %>% mutate(date = as.Date(date))
psiS_control <- read.csv("SCCII/psiS_hoelstein_control.csv") %>% mutate(date = as.Date(date), psiS_mean = psiS_mean / 1000,)

# ==========================================================================
# 2. SYNCHRONIZATION (COMMON TIME WINDOW)
# ==========================================================================

# 1. Calculate daily mean observed Leaf values
obs_psiL_daily <- psiL_control %>%
  group_by(date, species_name) %>%
  summarise(
    psiLmd_obs = mean(md_wp_av, na.rm = TRUE),
    psiLpd_obs = mean(pd_wp_av, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(species = species_name)

# 2. Prepare Observed Soil values (psiS_control)
# Note: Since psiS_control is site-level, we join it to all species later
obs_psiS_daily <- psiS_control %>%
  rename(psiS_obs = psiS_mean) %>%
  mutate(date = as.Date(date))

# 3. Combine everything using inner_join to find the intersection of all data
# This finds the "Common Time" across LPJ, Leaf Observations, and Soil Observations
combined_data <- lpj_base %>%
  inner_join(obs_psiL_daily, by = c("date", "species")) %>%
  inner_join(obs_psiS_daily, by = "date") %>%
  mutate(species = factor(species, levels = species_levels))

# ==========================================================================
# 3. VISUALIZATION
# ==========================================================================

p_psi_common <- ggplot(combined_data, aes(x = date)) +
  # --- LPJ Components (Colored by Species) ---
  geom_line(aes(y = psiL, color = species), linewidth = 1) + 
  geom_line(aes(y = psiX, color = species), linewidth = 0.7, linetype = "dotdash", alpha = 1) + 
  geom_line(aes(y = psiS, color = species), linewidth = 0.7, alpha = 1) + 
  
  # --- Observed Leaf Data (Grey Lines/Points) ---
  geom_line(aes(y = psiLmd_obs, color = "observed"), linewidth = 0.6, linetype = "solid") +
  geom_point(aes(y = psiLmd_obs, color = "observed"), size = 1.5, alpha = 0.8) +
  
  geom_line(aes(y = psiLpd_obs, color = "observed"), linewidth = 0.8, linetype = "dashed") +
  geom_point(aes(y = psiLpd_obs, color = "observed"), size = 1.5, shape = 1, alpha = 1) +
  
  # --- Observed Soil Data ---
  geom_line(aes(y = psiS_obs), color = "black", linewidth = 1, linetype = "dotted") +
  geom_point(aes(y = psiS_obs), color = "black", size = 1, shape = 4) + 
  
  # Faceting
  facet_wrap(~species, ncol = 2, scales = "free_y") +
  
  # Styling
  scale_color_manual(values = custom_colors) +
  labs(
    title = "water potential time series comparison (common time)",
    subtitle = paste0(climate_txt, "\ncolor: lpj (L solid, X dotdash, S thin) | grey: obs leaf (MD solid, PD dashed) | black: obs soil"),
    x = "date",
    y = expression(Psi~"(MPa)"),
    color = ""
  ) +
  base_theme

print(p_psi_common)

# Save and Print
if(!dir.exists("Figures/Hoelstein")) dir.create("Figures/Hoelstein", recursive = TRUE)
ggsave("Figures/Hoelstein/compare_psi_common_time.png", p_psi_common, width = 14, height = 10, bg = "white")

print(p_psi_common)

# ==========================================================================
# 4. CLIMATE TIME SERIES VISUALIZATION
# ==========================================================================

# Prepare Climate Data for the common time window
# We take only unique dates (climate is site-level, not species-level)
climate_common <- combined_data %>%
  distinct(date, temp_C, vpd, global_radiation, precipitation) %>%
  pivot_longer(
    cols = c(temp_C, vpd, global_radiation, precipitation),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(variable = factor(variable, levels = c("temp_C", "vpd", "global_radiation", "precipitation")))

# Define labels for facets
climate_labels <- c(
  temp_C = "Temp (°C)",
  vpd = "VPD (kPa)",
  global_radiation = "Radiation (W/m²)",
  precipitation = "Precip (mm)"
)

# ==========================================================================
# 4. CLIMATE TIME SERIES VISUALIZATION (ALL LINES)
# ==========================================================================

# Prepare labels for facets
climate_labels <- c(
  temp_C = "temperature (°C)",
  vpd = "VPD (kPa)",
  global_radiation = "global radiation (W/m²)",
  precipitation = "precipitation (mm)"
)

p_climate <- ggplot(climate_common, aes(x = date, y = value)) +
  # Use geom_line for all variables
  geom_line(color = "dodgerblue", linewidth = 0.7) +
  
  # Optional: Add a light area fill under precipitation to make it more visible
  # geom_area(data = subset(climate_common, variable == "precipitation"), 
  #           fill = "steelblue", alpha = 0.2) +
  
  # Faceting by climate variable
  facet_wrap(~variable, ncol = 1, scales = "free_y", 
             labeller = as_labeller(climate_labels)) +
  
  labs(
    title = "environmental drivers (common time window)",
    subtitle = "daily time series of meteorological forcing",
    x = "Date",
    y = NULL
  ) +
  base_theme +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

# Display and Save
print(p_climate)

ggsave("Figures/Hoelstein/climate_drivers_all_lines.png", p_climate, width = 12, height = 10, bg = "white")

# ==========================================================================
# 5. OPTIONAL: COMBINED PLOT (using patchwork)
# ==========================================================================
# If you have the patchwork library, you can stack them:
library(patchwork)
p_climate / p_psi_common + plot_layout(heights = c(1, 2))
