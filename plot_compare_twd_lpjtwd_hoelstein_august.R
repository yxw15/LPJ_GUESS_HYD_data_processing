# ==========================================================================
# TWD COMPARISON — AUGUST VALIDATION ONLY
# ==========================================================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(stringr)
library(purrr)
library(scales)
library(hydroGOF)

VALIDATION_MONTH <- 8  # August only

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")
dir.create("Figures/lpj_guess_stem_storage/validation_august/TWD", recursive = TRUE, showWarnings = FALSE)

species_levels <- c("oak", "beech", "spruce", "pine")
species_colors <- c("oak" = "#E69F00", "beech" = "#0072B2", "spruce" = "#009E73", "pine" = "#F0E442")

base_theme <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom", legend.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 11, face = "bold")
  )

# Load LPJ and obs TWD — filter to August
lpj_raw_twd <- read.csv("lpj_guess/lpj_guess_stem_storage/lpj_control_drought_ET_Gc_psi_leaf_psi_soil_psi_xylem_hydraulic_lag_kappy_s_min_mort_mort_cav_mort_greff_mort_min_stem_diameter_stem_rwc_twd.csv") %>%
  mutate(date = as.Date(date),
         species = factor(tolower(species), levels = species_levels),
         twd_model = twd * 1e6) %>%
  filter(month(date) == VALIDATION_MONTH, !is.na(species)) %>%
  select(date, species, treatment, twd_model)

dendro_obs <- list.files(path = "SCCII/point_dendro", pattern = "^Point_dendrometers_.*_archive\\.txt$", full.names = TRUE) %>%
  map_dfr(~ read.delim(.x), .id = "source_file")
tree_info <- read.csv("SCCII/tree_info.csv") %>% mutate(treatment = ifelse(treatment == "treatment", "drought", treatment))

obs_twd_daily <- dendro_obs %>% inner_join(tree_info, by = "tree_id") %>%
  mutate(date = as.Date(str_extract(timestamp_UTC, "\\d{4}-\\d{2}-\\d{2}")),
         species = factor(tolower(species), levels = species_levels),
         treatment = tolower(treatment)) %>%
  filter(month(date) == VALIDATION_MONTH) %>%
  group_by(date, species, treatment) %>%
  summarise(twd_mean_obs = mean(twd_micron_treenetproc, na.rm = TRUE), .groups = "drop")

combined_twd <- obs_twd_daily %>% inner_join(lpj_raw_twd, by = c("date", "species", "treatment"))

combined_std <- combined_twd %>% group_by(species, treatment) %>%
  mutate(std_twd_mean = (twd_mean_obs - mean(twd_mean_obs, na.rm=T))/sd(twd_mean_obs, na.rm=T),
         std_twd_model = (twd_model - mean(twd_model, na.rm=T))/sd(twd_model, na.rm=T)) %>%
  ungroup()

# KGE stats
stats_summary <- combined_std %>% filter(!is.na(std_twd_mean), !is.na(std_twd_model)) %>%
  group_by(species, treatment) %>%
  summarise(KGE = KGE(sim = std_twd_model, obs = std_twd_mean),
            n = n(), .groups = "drop")

# Time series plot
p_ts <- ggplot(combined_std, aes(x = date)) +
  geom_line(aes(y = std_twd_model, color = "LPJ"), linewidth = 0.7) +
  geom_line(aes(y = std_twd_mean, color = "Observed"), linewidth = 0.7, linetype = "dashed") +
  facet_grid(species ~ treatment) +
  scale_color_manual(values = c("LPJ" = "steelblue", "Observed" = "black")) +
  labs(title = "Standardized TWD: Model vs Observed [AUGUST]",
       subtitle = "AUGUST VALIDATION", y = "Standardized TWD", color = "") +
  base_theme

ggsave("Figures/lpj_guess_stem_storage/validation_august/TWD/timeseries_std.png", p_ts, width = 12, height = 10, dpi = 300)

# Scatter
p_scatter <- ggplot(combined_std, aes(x = std_twd_mean, y = std_twd_model, color = species)) +
  geom_point(alpha = 0.5) + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ treatment) +
  scale_color_manual(values = species_colors) +
  labs(title = "Standardized TWD Scatter [AUGUST]", x = "Observed (std)", y = "LPJ (std)") +
  base_theme + coord_fixed()

ggsave("Figures/lpj_guess_stem_storage/validation_august/TWD/scatter_std.png", p_scatter, width = 10, height = 7, dpi = 300)

write.csv(stats_summary, "Figures/lpj_guess_stem_storage/validation_august/TWD/kge_stats.csv", row.names = FALSE)

cat("\n*** August TWD validation complete ***\n")
