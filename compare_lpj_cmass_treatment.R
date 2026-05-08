# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# 1. Setup Environment and Metadata
setwd("~/Documents/LPJ_GUESS_hydraulic/results")

species_map <- tribble(
  ~species, ~colname, ~folder,
  "Beech",  "Fag_syl", "Beech",
  "Oak",    "Que_rob", "Oak",
  "Pine",   "Pin_syl", "Pine",
  "Spruce", "Pic_abi", "Spruce"
)

species_levels <- c("Oak", "Beech", "Spruce", "Pine")
cb_palette     <- c(Oak = "#E69F00", Beech = "#0072B2", Spruce = "#009E73", Pine = "#F0E442")

# 2. Define Theme
base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text       = element_text(color = "black", size = 14),
    legend.position   = "bottom",
    plot.title        = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5, size = 14),
    axis.title        = element_text(size = 16),
    axis.text.x       = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y       = element_text(angle = 0, hjust = 0.5, size = 12),
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_line(color = "grey92", linewidth = 0.25),
    strip.text        = element_text(size = 14, face = "bold")
  )

# 3. Helper Function to Read and Clean LPJ-GUESS Output
load_cmass <- function(path, scenario_name) {
  if (!file.exists(path)) {
    warning(paste("File not found:", path))
    return(NULL)
  }
  
  df <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
  
  # Ensure numeric conversion (LPJ-GUESS sometimes has messy white space)
  df <- df %>% mutate(across(everything(), as.numeric))
  
  # Pivot to long format to match our species_map
  df_long <- df %>%
    pivot_longer(cols = -c(Lon, Lat, Year), names_to = "colname", values_to = "cmass") %>%
    mutate(Scenario = scenario_name)
  
  return(df_long)
}

# ... (Keep previous libraries and theme definitions)

# 4. Loop through species and scenarios to load all data
all_data_list <- list()

for (i in 1:nrow(species_map)) {
  s_name   <- species_map$species[i]
  s_folder <- species_map$folder[i]
  
  path_ctrl <- file.path("control", s_folder, "cmass.out")
  path_drgt <- file.path("drought", s_folder, "cmass.out")
  
  all_data_list[[paste0(s_name, "_ctrl")]] <- load_cmass(path_ctrl, "Control")
  all_data_list[[paste0(s_name, "_drgt")]] <- load_cmass(path_drgt, "Drought")
}

# Combine and FILTER for years > 2000
full_data <- bind_rows(all_data_list) %>%
  left_join(species_map, by = "colname") %>%
  filter(!is.na(species)) %>%
  filter(Year > 2000) %>%  # <--- Filter applied here
  mutate(species = factor(species, levels = species_levels))

# 5. Create Time Series Plot
p <- ggplot(full_data, aes(x = Year, y = cmass, color = species, linetype = Scenario)) +
  geom_line(linewidth = 1.1, alpha = 0.8) +
  facet_wrap(~species, scales = "free_y") +
  scale_color_manual(values = cb_palette) +
  scale_linetype_manual(values = c("Control" = "solid", "Drought" = "dashed")) +
  scale_x_continuous(breaks = seq(2000, 2026, by = 5)) + # Adjust 'by' based on your end year
  labs(
    title = "LPJ-GUESS: Post-2000 Biomass Response",
    subtitle = "Impact of Drought Scenario on Carbon Mass (C-mass)",
    x = "Year",
    y = expression(C-mass~(kg~C~m^{-2})),
    color = "Species",
    linetype = "Scenario"
  ) +
  base_theme

# Print and Save
print(p)