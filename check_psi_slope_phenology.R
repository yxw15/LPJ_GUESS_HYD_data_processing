setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# -----------------------------
# 1. Load and prepare data
# -----------------------------
data <- read.delim("SCCII/psi_values_phenology/estimates_p12_p50_p88_archive.txt")

selected <- data %>%
  filter(species_id %in% c("Fs", "Pa", "Ps", "Qs")) %>%
  mutate(
    species_name = recode(species_id,
                          "Fs" = "Beech",
                          "Pa" = "Spruce",
                          "Ps" = "Pine",
                          "Qs" = "Oak"
    )
  )

# -----------------------------
# 2. Reshape to long format
# -----------------------------
selected_long <- selected %>%
  pivot_longer(
    cols = c(p12, p50, p88),
    names_to = "metric",
    values_to = "value"
  )

# -----------------------------
# 3. Plot all metrics together
# -----------------------------
ggplot(selected_long, aes(x = tree_id, y = value)) +
  geom_point() +
  facet_grid(metric ~ species_name, scales = "free_x") +
  labs(
    x = "Tree ID",
    y = "Value",
    title = "p12, p50, p88 per Tree by Species"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# -----------------------------
# 4. Mean per species (raw)
# -----------------------------
species_means <- selected %>%
  group_by(species_name) %>%
  summarise(
    mean_p12 = mean(p12, na.rm = TRUE),
    mean_p50 = mean(p50, na.rm = TRUE),
    mean_p88 = mean(p88, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

print(species_means)

# -----------------------------
# 5. Mean per species (tree-level first)
# -----------------------------
species_means_tree <- selected %>%
  group_by(species_name, tree_id) %>%
  summarise(
    p12 = mean(p12, na.rm = TRUE),
    p50 = mean(p50, na.rm = TRUE),
    p88 = mean(p88, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(species_name) %>%
  summarise(
    mean_p12 = mean(p12, na.rm = TRUE),
    mean_p50 = mean(p50, na.rm = TRUE),
    mean_p88 = mean(p88, na.rm = TRUE),
    n_trees = n(),
    .groups = "drop"
  )

print(species_means_tree)

# -----------------------------
# 6. Plot species means
# -----------------------------
species_means_long <- species_means_tree %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "metric",
    values_to = "value"
  )

ggplot(species_means_long, aes(x = species_name, y = value, fill = metric)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(
    x = "Species",
    y = "Mean value",
    title = "Mean p12, p50, p88 per Species"
  )

selected <- selected %>%
  mutate(
    m_cav = 2 / log10(p50 / p88)
  )

slope_species <- selected %>%
  group_by(species_name) %>%
  summarise(
    mean_m_cav = mean(m_cav, na.rm = TRUE),
    .groups = "drop"
  )

phenology_2022 <- read.delim("SCCII/psi_values_phenology/Phenology_2022_archive.txt")
phenology_2023 <- read.delim("SCCII/psi_values_phenology/Phenology_2023_archive.txt")

phenology_2022_2023 <- bind_rows(phenology_2022, phenology_2023) %>%
  filter(species_id %in% c("Fs", "Pa", "Ps", "Qs")) %>%
  mutate(species_name = case_when(
    species_id == "Fs" ~ "Beech",
    species_id == "Pa" ~ "Spruce",
    species_id == "Ps" ~ "Pine",
    species_id == "Qs" ~ "Oak",
    TRUE ~ NA_character_
  ))

phenology_2022_2023 <- phenology_2022_2023 %>%
  mutate(stage_desc = case_when(
    phenological_stage == "bb" ~ "Bud break",
    phenological_stage == "lo" ~ "Leaf opening",
    phenological_stage == "ld" ~ "Leaf discoloration",
    phenological_stage == "lf" ~ "Leaf fall",
    phenological_stage == "ns" ~ "Needle spreading",
    TRUE ~ NA_character_
  ))

unique_stages <- phenology_2022_2023 %>%
  select(species_name, phenological_stage, stage_desc) %>%
  distinct() %>%
  arrange(species_name, phenological_stage)

print(unique_stages)

stage_periods <- phenology_2022_2023 %>%
  group_by(species_name, year, phenological_stage, stage_desc) %>%
  filter(class >= 4) %>% 
  summarise(
    start_date = min(as.Date(date, format = "%d/%m/%Y")),
    end_date   = max(as.Date(date, format = "%d/%m/%Y")),
    .groups = "drop"
  ) %>%
  arrange(species_name, year, phenological_stage)

print(stage_periods)

write.csv(stage_periods, file = "SCCII/phenology.csv", row.names = FALSE)
