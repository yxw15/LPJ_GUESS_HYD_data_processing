# --------------------------------------------------
# 1. load libraries
# --------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

# --------------------------------------------------
# 2. set working directory
# --------------------------------------------------
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD/")

# define and create the output directory for figures if it doesn't exist
output_dir <- "Figures/lpj_guess_hyd_twd"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --------------------------------------------------
# 3. read and clean raw data
# --------------------------------------------------
sm_rew_sel_2018_2025 <- readRDS("SCCII/Soil_water_content_and_plant_water_uptake_depth/SM_REW_SEL_2018_2025.rds")
hoelstein_rwud_sel   <- readRDS("SCCII/Soil_water_content_and_plant_water_uptake_depth/Hoelstein_RWUD_sel.rds")

sm_rew_sel_2018_2025_clean <- sm_rew_sel_2018_2025 %>% drop_na()
hoelstein_rwud_clean       <- hoelstein_rwud_sel %>% drop_na()

# --------------------------------------------------
# 4. filter, transform, and order species/treatments
# --------------------------------------------------
plot_data <- hoelstein_rwud_clean %>%
  # keep only targeted species
  filter(species %in% c("Fagus sylvatica", "Quercus sp", "Pinus sylvestris", "Picea abies")) %>%
  # rename columns and map names/treatments to clean lowercase values
  rename(rwud = RWUD) %>%
  mutate(
    species = case_when(
      species == "Quercus sp"       ~ "oak",
      species == "Fagus sylvatica"  ~ "beech",
      species == "Picea abies"      ~ "spruce",
      species == "Pinus sylvestris" ~ "pine"
    ),
    # explicitly force the specified panel order: oak, beech, spruce, pine
    species = factor(species, levels = c("oak", "beech", "spruce", "pine")),
    # normalize treatment string data
    treatment = tolower(treatment),
    treatment = case_when(
      treatment == "treatment" ~ "drought",
      TRUE                     ~ treatment
    )
  )

# --------------------------------------------------
# 5. execute lowercase boxplot
# --------------------------------------------------
p_boxplot <- ggplot(plot_data, aes(x = treatment, y = rwud, fill = treatment)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  # arrange panels side-by-side using the custom ordered species factor
  facet_wrap(~ species, nrow = 1) +
  # apply custom color palette
  scale_fill_manual(values = c("control" = "#2b8cbe", "drought" = "#e66101")) +
  labs(
    title = "distribution of relative water uptake depth (rwud) by species and treatment",
    x = "treatment type",
    y = "relative water uptake depth (rwud)",
    fill = "treatment"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "plain", size = 12),
    legend.title = element_text(face = "plain"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

print(p_boxplot)

# --------------------------------------------------
# 6. save figure
# --------------------------------------------------
ggsave(
  filename = file.path(output_dir, "boxplot_hoelstein_rwud.png"), 
  plot = p_boxplot, 
  width = 9, 
  height = 5, 
  dpi = 300
)