# ==========================================================================
# 1. SETUP: GLOBAL CONSTANTS & THEMES
# ==========================================================================
library(dplyr)
library(ggplot2)
library(lubridate)
library(scales)
library(tidyr)
library(tidyverse)
library(patchwork)
library(glue)

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

# Physical constants for unit conversion (Source: Capture.JPG)
eta   <- 44.6   # mol m-3 (molar density of air at STP)
T0    <- 273    # K
h     <- 500    # Altitude in meters

# Aesthetics
species_levels <- c("Oak", "Beech", "Spruce", "Pine")
cb_palette     <- c(Oak = "#E69F00", Beech = "#0072B2", Spruce = "#009E73", Pine = "#F0E442")
custom_colors  <- c(cb_palette, "Observed" = "#999999") # Grey for observed

base_theme <- theme_minimal() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(color = "black", size = 14),
    legend.position = "bottom",
    plot.title  = element_text(hjust = 0.5, size = 18, color = "black"),
    axis.title  = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    panel.border = element_blank(),
    strip.text = element_text(size = 14)
  )

# ==========================================================================
# 2. DATA PREPARATION & FILTER DEFINITION
# ==========================================================================
lpj_output <- read.csv("lpj_guess/lpj_control_drought_Gc_psiL_psiX_psiS_stem_diameter_twd_stem_rwc.csv")

# Convert date column to Date format
lpj_output$date <- as.Date(lpj_output$date)

# ==========================================================================
# 3. VARIABLE DEFINITIONS
# ==========================================================================

# LPJ hydraulic variables
var_lpj <- c(
  "Gc" = "Gc (m/s)",
  "psi_leaf" = "psi_leaf (MPa)",
  "psi_soil" = "psi_soil (MPa)",
  "psi_xylem" = "psi_xylem (MPa)",
  "twd" = "TWD (m/m)",
  "stem_rwc" = "stem RWC"
)

# Climate variables
var_climate <- c(
  "temp_C" = "T (°C)",
  "vpd" = "VPD (kPa)",
  "global_radiation" = "radiation (W/m²)",
  "precipitation" = "precipitation (mm)"
)

# Combine all variables for the full dataset
all_variables <- c(names(var_lpj), names(var_climate))

# ==========================================================================
# PLOT 1: LPJ HYDRAULIC VARIABLES (Variables as Rows, Species as Columns)
# ==========================================================================

# Ensure factors are strictly ordered for clean plotting layout
lpj_output <- lpj_output %>%
  mutate(
    species = factor(species, levels = species_levels),
    treatment = factor(treatment, levels = c("control", "drought"))
  ) %>% 
  filter(year(date) > 2020)

plot_lpj <- lpj_output %>%
  pivot_longer(cols = all_of(names(var_lpj)), 
               names_to = "variable", 
               values_to = "value") %>%
  ggplot(aes(x = date, y = value, color = treatment)) +
  geom_line(linewidth = 0.4, alpha = 0.7) +
  scale_color_manual(values = c("control" = "dodgerblue", "drought" = "orange")) +
  
  # Crucial Change: Rows = Variables, Columns = Species
  # free_y allows each variable row to have its own unique scale 
  facet_grid(variable ~ species, scales = "free_y", 
             labeller = labeller(variable = var_lpj)) +
  
  labs(x = "Date", y = NULL, color = "Treatment", 
       title = "LPJ-GUESS-HYD-TWD Simulated Variables across Species") +
  base_theme 
# Display LPJ plot
print(plot_lpj)

# Save LPJ plot 
# Adjusted width to 16 to comfortably fit the 4 species columns side-by-side
ggsave("Figures/lpj_guess_hyd_twd/time_series_lpj_variables.png", 
       plot_lpj, width = 16, height = 18, dpi = 300)

# ==========================================================================
# 5. PLOT 2: CLIMATE VARIABLES (One Variable Per Row)
# ==========================================================================

plot_climate <- lpj_output %>%
  pivot_longer(cols = all_of(names(var_climate)), 
               names_to = "variable", 
               values_to = "value") %>%
  mutate(treatment = factor(treatment, levels = c("control", "drought"))) %>%
  ggplot(aes(x = date, y = value, color = treatment)) +
  geom_line(linewidth = 0.3, alpha = 0.7) +
  scale_color_manual(values = c("control" = "dodgerblue", "drought" = "orange")) +
  # Using facet_grid to stack variables exactly one per row
  facet_grid(variable ~ ., scales = "free_y", 
             labeller = labeller(variable = var_climate)) +
  labs(x = "Date", y = NULL, color = "Treatment", 
       title = "Climate Variables") +
  base_theme 

# Display climate plot
print(plot_climate)

# Save climate plot (Increased height to accommodate 4 single-stacked rows)
ggsave("Figures/lpj_guess_hyd_twd/time_series_climate_variables.png", 
       plot_climate, width = 14, height = 12, dpi = 300)

