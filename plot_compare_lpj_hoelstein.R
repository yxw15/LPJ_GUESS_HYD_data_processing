setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

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


species_levels <- c("Oak", "Beech", "Spruce", "Pine")

cb_palette <- c(
  Oak      = "#E69F00",
  Beech    = "#0072B2",
  Spruce   = "#009E73",
  Pine     = "#F0E442"
)

library(dplyr)

psiL_hoelstein <- read.csv("SCCII/psiL_hoelstein.csv")
sap_daily <- read.csv("SCCII/sap_daily.csv")
phenology <- read.csv("SCCII/phenology.csv")


# Constants and parameters from Capture.JPG
eta <- 44.6       # mol m-3 (molar density of air at STP)
T0 <- 273         # K
P_std <- 101.3    # kPa (standard pressure)
h <- 500          # Altitude in meters

# Constants and parameters
eta <- 44.6       # mol m-3 (Molar density of air at STP)
T0  <- 273        # K
h   <- 500        # Altitude in meters

sap_daily <- sap_daily %>%
  filter(vpd_crane > 0.05) %>%
  mutate(
    # 1. Transform SFD from [cm3 cm-2 h-1] to [kg m-2 s-1]
    # Factor (10/3600) converts cm3/cm2/h to kg/m2/s assuming water density of 1g/cm3
    sfd_kg = sfd * (10 / 3600),
    
    # 2. Calculate G_asw_mmol (Eqn 1) - Result is in mmol m-2 s-1
    # Using Ta (temp_crane), Fd (sfd_kg), D (vpd_crane), and h (altitude)
    G_asw_mmol = ((115.8 + 0.4236 * temp_crane) * sfd_kg / vpd_crane) * 
      eta * (T0 / (T0 + temp_crane)) * 
      exp(-0.00012 * h),
    
    # 3. Convert mmol m-2 s-1 to m s-1 (G_ms)
    # G_ms (m/s) = G_asw (mol m-2 s-1) / [eta * (T_factor) * (P_factor)]
    # We divide G_asw_mmol by 1000 to convert to mol to match eta (mol m-3)
    G_ms = (G_asw_mmol / 1000) / (eta * (T0 / (T0 + temp_crane)) * exp(-0.00012 * h))
  )

# Preview results
head(sap_daily %>% select(date, species, sfd, G_asw, G_ms))

head(sap_daily)

summary(sap_daily$G_ms)

library(dplyr)
library(ggplot2)

plot_df <- sap_daily %>%
  filter(treatment == "control") %>%
  mutate(
    species = factor(species, levels = species_levels)
  )

# Constants from Capture.JPG and image_79942f.png
eta <- 44.6       # mol m-3
T0  <- 273        # K
h   <- 500        # Altitude in meters

sap_daily <- sap_daily %>%
  filter(vpd_crane > 0.05) %>%
  mutate(
    # 1. Transform SFD to kg m-2 s-1
    sfd_kg = sfd * (10 / 3600),
    
    # 2. Calculate G_asw in mmol m-2 s-1 (per image_79942f.png)
    G_asw = ((115.8 + 0.4236 * temp_crane) * sfd_kg / vpd_crane) * 
      eta * (T0 / (T0 + temp_crane)) * 
      exp(-0.00012 * h),
    
    # 3. Convert to m/s (per Capture.JPG logic)
    # Divide by 1000 to convert mmol to mol
    G_ms = (G_asw / 1000) / (eta * (T0 / (T0 + temp_crane)) * exp(-0.00012 * h))
  )


library(ggplot2)
library(scales)
library(dplyr)

sap_daily <- sap_daily %>%
  filter(vpd_crane > 0.05) %>%
  mutate(
    # FIX: Convert the date column from character to Date class
    date = as.Date(date),
    
    # 1. Transform SFD to kg m-2 s-1
    sfd_kg = sfd * (10 / 3600),
    
    # 2. Calculate G_asw in mmol m-2 s-1 (per image_79942f.png)
    G_asw = ((115.8 + 0.4236 * temp_crane) * sfd_kg / vpd_crane) * 
      eta * (T0 / (T0 + temp_crane)) * 
      exp(-0.00012 * h)
  )


library(dplyr)
library(ggplot2)

# 1. Prepare the data: Calculate the mean of the top 10% for each day
plot_df_top10_mean <- sap_daily %>%
  filter(treatment == "control") %>% # Fixed typo: 'treatment'
  group_by(date, species) %>%
  # Keep only the top 10% of values for that specific day and species
  filter(G_asw >= quantile(G_asw, 0.90, na.rm = TRUE)) %>%
  # Calculate the mean of those top values
  summarise(G_asw_mean = mean(G_asw, na.rm = TRUE), .groups = "drop")

# 2. Create the plot
p <- ggplot(plot_df_top10_mean, aes(x = date, y = G_asw_mean, color = species)) +
  geom_point(size = 0.8, alpha = 0.6) +
  # Optional: add a light line to see the seasonal trend of the means
  geom_line(aes(group = year(date)), linewidth = 0.4, alpha = 0.3) +
  facet_wrap(~species, ncol = 2, scales = "free_y") +
  scale_color_manual(values = cb_palette) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  labs(
    x = "Year",
    y = expression(G[c]~(mmol~m^{-2}~s^{-1})),
    title = "Daily Mean of Top 10% Canopy Conductance (Hoelstein Control)"
  ) +
  base_theme +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.line = element_line(colour = "black")
  )

ggsave(
  filename = "Figures/Hoelstein/plot_Gc_top10mean_ts_control_hoelstein.png",
  plot = p,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

# 3. Save the figure
if (!dir.exists("Figures/Hoelstein")) dir.create("Figures/Hoelstein", recursive = TRUE)

ggsave(
  filename = "Figures/Hoelstein/plot_Gc_top10_ts_hoelstein.png",
  plot = p,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

library(dplyr)
library(ggplot2)
library(lubridate)

# Standardize column names and add source tags
df_obs <- plot_df_top10_mean %>%
  rename(gc_value = G_asw_mean) %>%
  mutate(Source = "Observed (Sap Flux)")

df_mod <- lpj_filtered %>%
  rename(gc_value = gc_mmol_mean) %>%
  mutate(Source = "Model (LPJ)")

# Bind them together
combined_df <- bind_rows(df_obs, df_mod)

combined_df <- combined_df %>%
  mutate(
    color_group = ifelse(Source == "Observed (Sap Flux)", "Observed", species)
  )

# Define a custom color scale: Species get cb_palette, Observed gets grey
custom_colors <- c(cb_palette, "Observed" = "grey70")

p_comp <- ggplot(combined_df, aes(x = date, y = gc_value, group = interaction(Source, species, year(date)))) +
  # Points: Observed in grey, Model in species colors
  geom_point(aes(color = color_group, shape = Source), size = 1, alpha = 0.5) +
  
  # Lines: Observed in grey, Model in species colors
  geom_line(aes(color = color_group, linetype = Source), linewidth = 0.5, alpha = 0.7) + 
  
  facet_wrap(~species, ncol = 2, scales = "free_y") +
  
  # Apply the custom color scale
  scale_color_manual(values = custom_colors) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  
  labs(
    x = "Year",
    y = expression(G[c]~(mmol~m^{-2}~s^{-1})),
    title = "Comparison: Observed vs. LPJ Model Canopy Conductance",
    subtitle = "Daily Mean of Top 10% Values (Observed = Grey, LPJ = Colored)"
  ) +
  
  base_theme +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.line = element_line(colour = "black")
  )

# Save the figure
ggsave(
  filename = "Figures/Hoelstein/compare_obs_vs_lpj_Gc.png",
  plot = p_comp,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)



lpj_gcwater_psiL <- read.csv("lpj_guess/lpj_gcwater_psiL.csv")

# Constants from Capture.JPG
eta <- 44.6       # mol m-3
T0  <- 273        # K
h   <- 500        # Altitude in meters

lpj_gcwater_psiL <- lpj_gcwater_psiL %>%
  mutate(
    # Ensure date is class <Date> for plotting later
    date = as.Date(date),
    
    # Transfer m/s to mmol m-2 s-1
    # gcwater (m/s) * molar density (mol/m3) * 1000 (mmol/mol)
    gc_mmol = gcwater * 
      (eta * (T0 / (T0 + temp_C)) * exp(-0.00012 * h)) * 
      1000
  )

lpj_filtered <- lpj_gcwater_psiL %>%
  filter(date %in% sap_daily$date) %>%
  group_by(date, species) %>%
  # Calculate mean of top 10% for the model data
  filter(gc_mmol >= quantile(gc_mmol, 0.90, na.rm = TRUE)) %>%
  summarise(gc_mmol_mean = mean(gc_mmol, na.rm = TRUE), .groups = "drop")

library(lubridate)

p_lpj <- ggplot(lpj_filtered, aes(x = date, y = gc_mmol_mean, color = species)) +
  geom_point(size = 0.8, alpha = 0.6) +
  # Connect lines only within each year
  geom_line(aes(group = year(date)), linewidth = 0.4, alpha = 0.3) + 
  facet_wrap(~species, ncol = 2, scales = "free_y") +
  scale_color_manual(values = cb_palette) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  labs(
    x = "Year",
    y = expression(G[c]~(mmol~m^{-2}~s^{-1})),
    title = "LPJ Model: Daily Mean of Top 10% Canopy Conductance"
  ) +
  base_theme +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.line = element_line(colour = "black")
  )

# Save the figure
ggsave(
  filename = "Figures/Hoelstein/plot_Gc_top10mean_ts_control_lpj.png",
  plot = p_lpj,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)



sap_psiL_daily <- sap_daily %>%
  # 1. Join with the water potential data
  # We use left_join to keep all sap flow records
  left_join(psiL_hoelstein, by = c("date", "tree_id", "tree_nr")) %>%
  
  # 2. Rename columns to your specific requirements
  rename(
    psiLmd = md_wp_av,
    psiLpd = pd_wp_av
  ) %>%
  
  # 3. Select and reorder columns as requested
  # Note: 'species' from sap_daily is used for the species name
  select(
    date, 
    species, 
    tree_id, 
    tree_nr, 
    psiLmd, 
    psiLpd, 
    treatment, 
    relhum_crane, 
    temp_crane, 
    radiation_crane, 
    vpd_crane, 
    precip_crane, 
    temp_meteo, 
    relhum_meteo, 
    wind_meteo, 
    radiation_meteo, 
    vpd_meteo, 
    precip_meteo, 
    sfd, 
    G_asw, 
    G_ms
  )

# Preview the result
head(sap_psiL_daily)

sap_psiL_daily_clean <- na.omit(sap_psiL_daily)

write.csv(sap_psiL_daily_clean, file = "SCCII/sap_psiL_daily_clean.csv", row.names = FALSE)


library(dplyr)

# 1. Prepare the phenology ranges
# We create a clean lookup table for the specific stages you want
pheno_limits <- phenology %>%
  filter(
    (species_name == "Beech" & phenological_stage %in% c("lo", "lf")) |
      (species_name %in% c("Pine", "Spruce") & phenological_stage == "ns")
  ) %>%
  # Reshape so we have start and end on one row per species/year
  # For Pine/Spruce, start/end comes from the 'ns' row
  # For Beech/Oak, we need the start of 'lo' and start of 'lf'
  group_by(species_name, year) %>%
  summarize(
    season_start = ifelse("lo" %in% phenological_stage, 
                          start_date[phenological_stage == "lo"], 
                          start_date[phenological_stage == "ns"]),
    season_end = ifelse("lf" %in% phenological_stage, 
                        start_date[phenological_stage == "lf"], 
                        end_date[phenological_stage == "ns"]),
    .groups = "drop"
  )

# 2. Join and Filter the main data
sap_psiL_daily_phenology_filtered <- sap_psiL_daily_clean %>%
  mutate(
    year = as.numeric(format(as.Date(date), "%Y")),
    # Create a helper column to match Oak to Beech phenology
    pheno_match_species = ifelse(species == "Oak", "Beech", species)
  ) %>%
  left_join(pheno_limits, by = c("pheno_match_species" = "species_name", "year" = "year")) %>%
  filter(
    date >= as.Date(season_start) & 
      date <= as.Date(season_end)
  ) %>%
  select(-pheno_match_species, -season_start, -season_end, -year) # Clean up helper columns

# Check results
nrow(sap_psiL_daily_phenology_filtered)
head(sap_psiL_daily_phenology_filtered)

