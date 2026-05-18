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

base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text       = element_text(color = "black", size = 12),
    legend.position   = "bottom",
    plot.title        = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5, size = 12, color = "grey30"),
    axis.title        = element_text(size = 14),
    strip.text        = element_text(size = 12, face = "bold"),
    panel.grid.major  = element_line(color = "grey90")
  )

# Data Ingestion
lpj_base     <- read.csv("lpj_guess/lpj_Gc_psiL_psiX_psiS_drought.csv") %>% mutate(date = as.Date(date))
psiL_drought <- read.csv("SCCII/psiL_hoelstein_drought.csv") %>% mutate(date = as.Date(date))
psiS_drought <- read.csv("SCCII/psiS_hoelstein_drought.csv") %>% 
  mutate(date = as.Date(date), psiS_mean = psiS_mean / 1000)

# ==========================================================================
# 2. DATA PREPARATION (DAILY MEANS)
# ==========================================================================

# Aggregate Observed Leaf Data
obs_psiL_daily <- psiL_drought %>%
  group_by(date, species_name) %>%
  summarise(
    psiLmd_obs = mean(md_wp_av, na.rm = TRUE),
    psiLpd_obs = mean(pd_wp_av, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(species = species_name)

# Aggregate Observed Soil Data
obs_psiS_daily <- psiS_drought %>%
  group_by(date) %>%
  summarise(psiS_obs = mean(psiS_mean, na.rm = TRUE), .groups = "drop")

# ==========================================================================
# 3. FIGURE 1: COMMON TIME (INTERSECTION)
# ==========================================================================

data_common <- lpj_base %>%
  inner_join(obs_psiL_daily, by = c("date", "species")) %>%
  inner_join(obs_psiS_daily, by = "date") %>%
  mutate(species = factor(species, levels = species_levels))

p_common <- ggplot(data_common, aes(x = date)) +
  # --- 1. MODEL COMPONENTS (Colored by Species) ---
  geom_line(aes(y = psiS, color = species), linewidth = 0.5, alpha = 0.7) +
  geom_line(aes(y = psiX, color = species), linetype = "dotdash", linewidth = 0.7) +
  geom_line(aes(y = psiL, color = species), linewidth = 1) +
  
  # --- 2. OBSERVED COMPONENTS (Grey & Black Style) ---
  # Observed Midday Leaf (Grey Solid Line + Filled Triangle)
  geom_line(aes(y = psiLmd_obs), color = "grey40", linewidth = 0.6, linetype = "solid") +
  geom_point(aes(y = psiLmd_obs), color = "grey40", shape = 17, size = 1.8) + 
  
  # Observed Predawn Leaf (Grey Dashed Line + Open Circle)
  geom_line(aes(y = psiLpd_obs), color = "grey40", linewidth = 0.6, linetype = "dashed") +
  geom_point(aes(y = psiLpd_obs), color = "grey40", shape = 1, size = 1.8) + 
  
  # Observed Soil Data (Black Dotted Line + X shape)
  geom_line(aes(y = psiS_obs), color = "black", linetype = "dotted", linewidth = 0.8) +
  geom_point(aes(y = psiS_obs), color = "black", shape = 4, size = 1.5) +
  
  # Faceting & Styling
  facet_wrap(~species, ncol = 4, scales = "free_y") +
  scale_color_manual(values = cb_palette) +
  ylim(-3.2, 0) +
  labs(
    title = "water potential time series comparison (common time, drought experiment since 2023)",
    subtitle = "Model (Color): Solid=L, Dotdash=X, Thin=S | Obs (Grey): ▲=MD Solid, ○=PD Dashed | Obs (Black): ✖=Soil Dotted",
    x = NULL, 
    y = expression(Psi~"(MPa)"), 
    color = "Species Model"
  ) +
  base_theme

print(p_common)
ggsave("Figures/Hoelstein/psi_common_time_drought.png", p_common, width = 15, height = 8, bg = "white")

data_common_withleaf <- data_common %>% filter(psiL < -0.5)

p_common_withleaf <- ggplot(data_common_withleaf, aes(x = date)) +
  # --- 1. MODEL COMPONENTS (Colored by Species) ---
  geom_line(aes(y = psiS, color = species), linewidth = 0.5, alpha = 0.7) +
  geom_line(aes(y = psiX, color = species), linetype = "dotdash", linewidth = 0.7) +
  geom_line(aes(y = psiL, color = species), linewidth = 1) +
  
  # --- 2. OBSERVED COMPONENTS (Grey & Black Style) ---
  # Observed Midday Leaf (Grey Solid Line + Filled Triangle)
  geom_line(aes(y = psiLmd_obs), color = "grey40", linewidth = 0.6, linetype = "solid") +
  geom_point(aes(y = psiLmd_obs), color = "grey40", shape = 17, size = 1.8) + 
  
  # Observed Predawn Leaf (Grey Dashed Line + Open Circle)
  geom_line(aes(y = psiLpd_obs), color = "grey40", linewidth = 0.6, linetype = "dashed") +
  geom_point(aes(y = psiLpd_obs), color = "grey40", shape = 1, size = 1.8) + 
  
  # Observed Soil Data (Black Dotted Line + X shape)
  geom_line(aes(y = psiS_obs), color = "black", linetype = "dotted", linewidth = 0.8) +
  geom_point(aes(y = psiS_obs), color = "black", shape = 4, size = 1.5) +
  
  # Faceting & Styling
  facet_wrap(~species, ncol = 4, scales = "free_y") +
  scale_color_manual(values = cb_palette) +
  ylim(-3.2, 0) +
  labs(
    title = "water potential time series comparison (common time, drought experiment since 2023)",
    subtitle = "Model (Color): Solid=L, Dotdash=X, Thin=S | Obs (Grey): ▲=MD Solid, ○=PD Dashed | Obs (Black): ✖=Soil Dotted",
    x = NULL, 
    y = expression(Psi~"(MPa)"), 
    color = "Species Model"
  ) +
  base_theme

print(p_common_withleaf)
ggsave("Figures/Hoelstein/psi_common_time_drought_withleaf.png", p_common_withleaf, width = 15, height = 8, bg = "white")

# ==========================================================================
# 4. FIGURE 2: LONG-TERM (2018–2025)
# ==========================================================================

data_longterm <- lpj_base %>%
  filter(date >= "2018-01-01" & date <= "2025-12-31") %>%
  left_join(obs_psiL_daily, by = c("date", "species")) %>%
  left_join(obs_psiS_daily, by = "date") %>%
  mutate(species = factor(species, levels = species_levels))

p_longterm <- ggplot(data_longterm, aes(x = date)) +
  # --- 1. MODEL COMPONENTS (Colored by Species) ---
  geom_line(aes(y = psiS, color = species), linewidth = 0.4, alpha = 0.6) +
  geom_line(aes(y = psiX, color = species), linetype = "dotdash", linewidth = 0.5, alpha = 0.8) +
  geom_line(aes(y = psiL, color = species), linewidth = 0.8) +
  # --- 2. OBSERVED COMPONENTS (Specific Colors and Shapes) ---
  # Observed Soil (Red Triangles)
  geom_point(aes(y = psiS_obs), color = "grey", shape = 17, size = 0.8, alpha = 0.7, na.rm = TRUE) +
  # Observed Midday Leaf (Blue Triangles)
  geom_point(aes(y = psiLmd_obs), color = "black", shape = 17, size = 1.5, alpha = 0.8, na.rm = TRUE) + 
  # Observed Predawn Leaf (Black Open Circles)
  geom_point(aes(y = psiLpd_obs), color = "black", shape = 1, size = 1.5, na.rm = TRUE) + 
  ylim(-3, 0) +
  # Faceting & Styling
  facet_wrap(~species, ncol = 4, scales = "free_y") +
  scale_color_manual(values = cb_palette) +
  # Scale x-axis for better readability over long term
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title = "water potential time series (2018–2025, drought experiment since 2023)",
    subtitle = "LPJ-GUESS-HYD (color): solid=L, dotdash=X, thin=S | Obs: black ▲=MD, black ○=PD, grey ▲=soil",
    x = "", 
    y = expression(Psi~"(MPa)"), 
    color = ""
  ) +
  base_theme

print(p_longterm)
ggsave("Figures/Hoelstein/psi_longterm_2018_2025_drought.png", p_longterm, width = 15, height = 8, bg = "white")


# ==========================================================================
# 1. PREPARE DATA FOR SCATTER PLOT (5 PAIRS)
# ==========================================================================

# Define the pairs for comparison
plot_comparison <- data_longterm %>%
  filter(!is.na(psiLmd_obs) | !is.na(psiLpd_obs) | !is.na(psiS_obs)) %>%
  # 1. Create the primary pairs using pivot_longer
  pivot_longer(
    cols = c(psiLmd_obs, psiLpd_obs, psiS_obs),
    names_to = "obs_type",
    values_to = "obs_value"
  ) %>%
  mutate(model_value = case_when(
    obs_type == "psiLmd_obs" ~ psiL,      # psiL vs Midday Obs
    obs_type == "psiLpd_obs" ~ psiL,      # psiL vs Predawn Obs
    obs_type == "psiS_obs"   ~ psiS,      # psiS vs Soil Obs
    TRUE ~ NA_real_
  )) %>%
  # 2. Add the specific Xylem comparisons
  bind_rows(
    # Pair: psiX vs Predawn Obs
    data_longterm %>%
      filter(!is.na(psiLpd_obs)) %>%
      transmute(species, date, obs_type = "psiX_vs_pd_obs", 
                obs_value = psiLpd_obs, model_value = psiX),
    # Pair: psiX vs Midday Obs
    data_longterm %>%
      filter(!is.na(psiLmd_obs)) %>%
      transmute(species, date, obs_type = "psiX_vs_md_obs", 
                obs_value = psiLmd_obs, model_value = psiX)
  ) %>%
  filter(!is.na(obs_value) & !is.na(model_value)) %>%
  # Set species order and clean labels for the legend
  mutate(
    species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")),
    obs_type = factor(obs_type, levels = c("psiLmd_obs", "psiLpd_obs", "psiS_obs", 
                                           "psiX_vs_md_obs", "psiX_vs_pd_obs"))
  )

# ==========================================================================
# 2. VISUALIZATION (1:1 Comparison)
# ==========================================================================

# ==========================================================================
# 1. CALCULATE UNIFIED LIMITS
# ==========================================================================
# Find the absolute min and max across all plotted values to ensure square scales
all_values <- c(plot_comparison$obs_value, plot_comparison$model_value)
plot_range <- range(all_values, na.rm = TRUE)

# Optional: Add a small buffer (e.g., 5%) so points don't touch the edges
buffer <- diff(plot_range) * 0.05
final_limits <- c(plot_range[1] - buffer, plot_range[2] + buffer)

# ==========================================================================
# 2. VISUALIZATION (SQUARE 1:1)
# ==========================================================================

p_scatter_final <- ggplot(plot_comparison, aes(x = obs_value, y = model_value, color = obs_type, shape = obs_type)) +
  # 1:1 Reference Line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  
  # Data Points
  geom_point(size = 2.5, alpha = 0.7) +
  
  # Faceting by Species
  facet_wrap(~species, ncol = 4) +
  
  # FORCE EQUAL AXES: This makes the 1:1 line a true 45-degree angle
  coord_equal(xlim = final_limits, ylim = final_limits) +
  
  # Scale Configuration
  scale_color_manual(
    values = c(
      "psiLmd_obs"     = "green4", 
      "psiLpd_obs"     = "purple", 
      "psiS_obs"       = "orange", 
      "psiX_vs_md_obs" = "dodgerblue", 
      "psiX_vs_pd_obs" = "pink"
    ),
    labels = c(
      "psiLmd_obs"     = expression(Psi["  L"]~vs~Obs~MD),
      "psiLpd_obs"     = expression(Psi["  L"]~vs~Obs~PD),
      "psiS_obs"       = expression(Psi["  S"]~vs~Obs~Soil),
      "psiX_vs_md_obs" = expression(Psi["  X"]~vs~Obs~MD),
      "psiX_vs_pd_obs" = expression(Psi["  X"]~vs~Obs~PD)
    )
  ) +
  scale_shape_manual(
    values = c(17, 16, 15, 2, 1), 
    labels = c(
      "psiLmd_obs"     = expression(Psi["  L"]~vs~Obs~MD),
      "psiLpd_obs"     = expression(Psi["  L"]~vs~Obs~PD),
      "psiS_obs"       = expression(Psi["  S"]~vs~Obs~Soil),
      "psiX_vs_md_obs" = expression(Psi["  X"]~vs~Obs~MD),
      "psiX_vs_pd_obs" = expression(Psi["  X"]~vs~Obs~PD)
    )
  ) +
  
  labs(
    title = "water potential model vs observation: 1:1 scatter plot",
    subtitle = "drought experiment since 2023 | dashed = 1:1 reference line",
    x = expression(Observed~Psi[" "]~"(MPa)"),
    y = expression(Simulated~LPJ-GUESS~Psi[" "]~"(MPa)"),
    color = "comparison Pair",
    shape = "comparison Pair"
  ) +
  ylim(-4, 0) +
  base_theme +
  theme(
    legend.position = "bottom", 
    legend.text.align = 0,
    panel.grid.minor = element_blank(),
    aspect.ratio = 1, plot.margin = margin(20, 20, 20, 20)
  )

# Print and Save
print(p_scatter_final)
ggsave("Figures/Hoelstein/psi_1to1_drought.png", p_scatter_final, width = 16, height = 6, bg = "white")

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
