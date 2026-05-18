# ==========================================================================
# 1. SETUP: GLOBAL CONSTANTS & THEMES
# ==========================================================================
library(dplyr)
library(ggplot2)
library(lubridate)
library(scales)
library(tidyr)

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
lpj_Gc_psiL <- read.csv("lpj_guess/lpj_Gc_psiL_psiX_psiS_drought.csv")
lpj_base <- lpj_Gc_psiL %>% mutate(date = as.Date(date))

# Define the strict climate filter dates
climate_filter_dates <- lpj_base %>%
  filter(
    (temp_C + 273.15) > 287.15, # > 14°C
    precipitation < 1,           #
    global_radiation > 150,      #
    vpd > 0.3,                   #
    month(date) %in% c(6, 7, 8, 9) # Summer months
  ) %>%
  pull(date) %>% unique()

# ==========================================================================
# 3. PROCESSING FUNCTION
# ==========================================================================
process_gc_comparison <- function(apply_climate_filter = FALSE) {
  
  # A. Process Observations
  obs <- read.csv("SCCII/sap_daily.csv") %>%
    mutate(date = as.Date(date)) %>%
    filter(treatment == "control", vpd_crane > 0.05)
  
  if(apply_climate_filter) {
    obs <- obs %>% filter(date %in% climate_filter_dates)
  }
  
  obs_clean <- obs %>%
    mutate(
      sfd_kg = sfd * (10 / 3600), 
      G_c = ((115.8 + 0.4236 * temp_crane) * sfd_kg / vpd_crane) * 
        eta * (T0 / (T0 + temp_crane)) * exp(-0.00012 * h) #
    ) %>%
    group_by(date, species) %>%
    filter(G_c >= quantile(G_c, 0.90, na.rm = TRUE)) %>%
    summarise(gc_value = mean(G_c, na.rm = TRUE), .groups = "drop") %>%
    mutate(Source = "Observed (Sap Flux)")
  
  # B. Process Model (Synchronized to Observed dates)
  mod_clean <- lpj_base %>%
    filter(date %in% obs_clean$date) %>% 
    mutate(
      gc_mmol = Gc * (eta * (T0 / (T0 + temp_C)) * exp(-0.00012 * h)) * 1000 #
    ) %>%
    group_by(date, species) %>%
    filter(gc_mmol >= quantile(gc_mmol, 0.90, na.rm = TRUE)) %>%
    summarise(gc_value = mean(gc_mmol, na.rm = TRUE), .groups = "drop") %>%
    mutate(Source = "Model (LPJ)")
  
  # C. Combine
  bind_rows(obs_clean, mod_clean) %>%
    mutate(
      species = factor(species, levels = species_levels),
      color_group = ifelse(Source == "Observed (Sap Flux)", "Observed", as.character(species))
    )
}

# ==========================================================================
# 4. COORDINATE SYNCHRONIZATION
# ==========================================================================
# Generate datasets first to determine global ranges
data_full     <- process_gc_comparison(apply_climate_filter = FALSE) %>% filter (gc_value < 2000)
data_filtered <- process_gc_comparison(apply_climate_filter = TRUE) %>% filter (gc_value < 2000)

# Calculate global limits to ensure identical axis ranges
# y_limit_max <- max(data_full$gc_value, na.rm = TRUE) 
y_limit_max <- 2000
x_limit_range <- range(data_full$date, na.rm = TRUE)

# ==========================================================================
# 5. VISUALIZATION
# ==========================================================================

# --- Plot 1: Standard Synchronized Comparison ---
p1 <- ggplot(data_full, aes(x = date, y = gc_value, group = interaction(Source, species, year(date)))) +
  geom_point(aes(color = color_group, shape = Source), size = 1, alpha = 0.5) +
  geom_line(aes(color = color_group, linetype = Source), linewidth = 0.5, alpha = 0.7) + 
  facet_wrap(~species, ncol = 4) + # Fixed scales for better comparison
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(limits = c(0, y_limit_max)) +
  scale_x_date(limits = x_limit_range, date_labels = "%Y", date_breaks = "1 year") +
  labs(title = "canopy conductance: sap flux density vs. LPJ-GUESS-HYD (drought experiment since 2023)", 
       y = expression(G[c]~(mmol~m^{-2}~s^{-1}))) +
  base_theme

print(p1)

ggsave("Figures/Hoelstein/compare_obs_vs_lpj_Gc_full_drought.png", p1, width = 15, height = 8)

# --- Plot 2: Climate Filtered Comparison ---
p2 <- ggplot(data_filtered, aes(x = date, y = gc_value, group = interaction(Source, species, year(date)))) +
  geom_point(aes(color = color_group, shape = Source), size = 1.2, alpha = 0.6) +
  geom_line(aes(color = color_group, linetype = Source), linewidth = 0.6, alpha = 0.7) + 
  facet_wrap(~species, ncol = 4) + # Fixed scales for better comparison
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(limits = c(0, y_limit_max)) +
  scale_x_date(limits = x_limit_range, date_labels = "%Y", date_breaks = "1 year") +
  labs(title = "canopy conductance: sap flux density vs. LPJ-GUESS-HYD (drought experiment since 2023)", 
       y = expression(G[c]~(mmol~m^{-2}~s^{-1}))) +
  base_theme

ggsave("Figures/Hoelstein/compare_obs_vs_lpj_Gc_climate_filtered_drought.png", p2, width = 15, height = 8)

