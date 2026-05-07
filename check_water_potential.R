# ============================================================
# LIBRARIES
# ============================================================

library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)

# ============================================================
# 1) WORKING DIRECTORY
# ============================================================

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD/SCCII")

# ============================================================
# 2) LOAD DATA
# ============================================================

psi.leaf <- read.csv("psiL_hoelstein.csv")
psi.soil <- readRDS("soilwaterpotential/Hoelstein_SWP_2023_2025_ALL_rescaled.rds")
sap_flux_gc <- read.csv("sap_flux_gc.csv")

# ============================================================
# 3) TREE METADATA
# ============================================================

tree_info <- sap_flux_gc %>%
  distinct(tree_id, species, tree_nr, treatment)

write.csv(tree_info, "tree_info.csv", row.names = FALSE)

# ============================================================
# 4) STANDARDISE KEYS
# ============================================================

psi.leaf <- psi.leaf %>%
  mutate(tree_nr = as.character(tree_nr),
         date = as.Date(date))

tree_info <- tree_info %>%
  mutate(tree_nr = as.character(tree_nr))

# ============================================================
# 5) ADD TREATMENT TO LEAF DATA
# ============================================================

psi.leaf <- psi.leaf %>%
  left_join(
    tree_info %>% select(tree_nr, treatment),
    by = "tree_nr"
  )

# ============================================================
# 6) SOIL DATA PROCESSING (kPa → MPa + DAILY MEAN)
# ============================================================

psi.soil <- psi.soil %>%
  mutate(
    SWP_site_rescaled = SWP_site_rescaled / 1000,
    date = as.Date(timestamp_UTC)
  )

psi.soil_daily <- psi.soil %>%
  group_by(date, treatment) %>%
  summarise(
    psi_soil = mean(SWP_site_rescaled, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 7) MERGE LEAF + SOIL
# ============================================================

psi.leaf_soil <- psi.leaf %>%
  left_join(
    psi.soil_daily,
    by = c("date", "treatment")
  )

# ============================================================
# 8) CLEAN DATA + DERIVED VARIABLES
# ============================================================

psi.leaf_soil_clean <- psi.leaf_soil %>%
  filter(
    !is.na(psi_soil),
    !is.na(md_wp_av)
  ) %>%
  mutate(
    psi_Lmd = md_wp_av,
    psi_Lpd = pd_wp_av,
    delta_psi_Lmd_soil = md_wp_av - psi_soil,
    delta_psi_Lmd_Lpd = md_wp_av - pd_wp_av
  ) %>% 
  filter(psi_Lpd > -3.5)

# ============================================================
# 9) SPECIES ORDER
# ============================================================

psi.leaf_soil_clean$species_name <- factor(
  psi.leaf_soil_clean$species_name,
  levels = c("Oak", "Beech", "Spruce", "Pine")
)

# ============================================================
# ============================================================
# PLOT 1 — ABSOLUTE WATER POTENTIALS
# ============================================================

df_abs <- psi.leaf_soil_clean %>%
  select(date, species_name, treatment,
         psi_soil, psi_Lmd, psi_Lpd) %>%
  pivot_longer(
    cols = c(psi_soil, psi_Lmd, psi_Lpd),
    names_to = "variable",
    values_to = "value"
  )

ggplot(df_abs,
       aes(x = date, y = value,
           color = variable,
           linetype = treatment,
           group = interaction(variable, treatment))) +
  geom_line(linewidth = 1) +
  facet_wrap(~species_name, ncol = 2) +
  scale_color_manual(values = c(
    psi_soil = "orange",
    psi_Lmd = "green",
    psi_Lpd = "dodgerblue"
  )) +
  scale_linetype_manual(values = c(
    control = "solid",
    treatment = "dashed"
  )) +
  labs(
    x = "Date",
    y = "Water potential (MPa)",
    color = "Variable",
    linetype = "Treatment"
  ) +
  theme_classic()

# ============================================================
# PLOT 2 — Δψ (GRADIENTS)
# ============================================================

df_delta <- psi.leaf_soil_clean %>%
  select(date, species_name, treatment,
         delta_psi_Lmd_soil,
         delta_psi_Lmd_Lpd) %>%
  pivot_longer(
    cols = c(delta_psi_Lmd_soil, delta_psi_Lmd_Lpd),
    names_to = "variable",
    values_to = "value"
  )

ggplot(df_delta,
       aes(x = date, y = value,
           color = variable,
           linetype = treatment,
           group = interaction(variable, treatment))) +
  geom_line(linewidth = 1) +
  facet_wrap(~species_name, ncol = 2) +
  scale_color_manual(values = c(
    delta_psi_Lmd_soil = "orange",
    delta_psi_Lmd_Lpd = "dodgerblue"
  )) +
  scale_linetype_manual(values = c(
    control = "solid",
    treatment = "dashed"
  )) +
  labs(
    x = "Date",
    y = "Δ Water potential (MPa)",
    color = "Variable",
    linetype = "Treatment"
  ) +
  theme_classic()

delta_summary <- psi.leaf_soil_clean %>%
  group_by(species_name, treatment) %>%
  summarise(
    mean_delta_Lmd_soil = mean(delta_psi_Lmd_soil, na.rm = TRUE),
    mean_delta_Lmd_Lpd  = mean(delta_psi_Lmd_Lpd, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

print(delta_summary)
