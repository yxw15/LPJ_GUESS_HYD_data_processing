# =========================================================
# 0. SETUP & PATHS
# =========================================================
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

library(tidyverse)
library(lubridate)
library(minpack.lm)
library(plantecophys)

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

# Create output directory
dir.create("Figures/lpj_gcwater_psiL_filter", recursive = TRUE, showWarnings = FALSE)

# Color Palette for Species
cb_palette <- c(
  Oak      = "#E69F00",
  Beech    = "#0072B2",
  Spruce   = "#009E73",
  Pine     = "#F0E442"
)


# =========================================================
# 1. METEO DATA (Hoelstein / RUE)
# =========================================================
data_RUE <- read.csv(
  "MeteoSwiss/MeteoSwiss_station_to_netcdf_daytime/all_stations_RUE_replaced_daytime.csv"
) %>%
  mutate(
    date = as.Date(date),
    temp_C = temperature - 273.15,
    rh_100 = relative_humidity *100,
    # Use plantecophys for VPD
    vpd = RHtoVPD(RH = rh_100, TdegC = temp_C)
  ) %>%
  filter(station_abbr == "RUE") %>% 
  filter(
    station_abbr == "RUE",
    temperature > 287.15, # >14°C
    precipitation < 1,
    global_radiation > 150,
    vpd > 0.3,
    month(date) %in% c(6, 7, 8, 9)
  )

# =========================================================
# 2. HELPERS
# =========================================================
# Convert LPJ Day/Year to R Date
convert_lpj_date <- function(df) {
  df %>% mutate(date = as.Date(Day + 1, origin = paste0(Year, "-01-01")))
}

# Generic LPJ File Reader
read_lpj_file <- function(path, colname, newname) {
  read.table(path, header = TRUE) %>%
    convert_lpj_date() %>%
    semi_join(data_RUE, by = "date") %>%
    transmute(date, value = .data[[colname]]) %>%
    rename(!!newname := value)
}

# =========================================================
# 3. DATA LOADING & MERGING
# =========================================================
species_map <- tribble(
  ~species,   ~colname,
  "Beech",    "Fag_syl",
  "Oak",      "Que_rob",
  "Pine",     "Pin_syl",
  "Spruce",   "Pic_abi"
)

final_daily <- map_dfr(1:nrow(species_map), function(i) {
  sp  <- species_map$species[i]
  col <- species_map$colname[i]
  base_path <- file.path("results/hoelstein", paste0(sp, "_hoelstein"))
  
  gc  <- read_lpj_file(file.path(base_path, "dgcwater.out"), col, "gcwater")
  psi <- read_lpj_file(file.path(base_path, "dpsileaf.out"), col, "psiL")
  
  full_join(gc, psi, by = "date") %>%
    mutate(species = sp)
}) %>%
  left_join(data_RUE %>% select(date, vpd, temp_C, global_radiation, precipitation), by = "date") %>%
  drop_na(gcwater, psiL)

# write.csv(final_daily, file = "lpj_guess/lpj_gcwater_psiL.csv", row.names = FALSE)
write.csv(final_daily, file = "lpj_guess/lpj_gcwater_psiL_climate_filter.csv", row.names = FALSE)

# =========================================================
# 4. NORMALIZATION & AGGREGATION
# =========================================================

# Daily Normalization
final_daily <- final_daily %>%
  group_by(species) %>%
  mutate(
    gcmax = quantile(gcwater, 0.95, na.rm = TRUE),
    gc_rel = pmin(gcwater / gcmax, 1) # Relative gc (0-1)
  ) %>%
  ungroup() %>% 
  mutate(species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")))

# Monthly Aggregation
final_monthly <- final_daily %>%
  mutate(month_date = floor_date(date, unit = "month")) %>%
  group_by(species, month_date) %>%
  summarise(
    psiL = mean(psiL, na.rm = TRUE),
    gcwater = mean(gcwater, na.rm = TRUE),
    gc_rel = mean(gc_rel, na.rm = TRUE),
    .groups = "drop"
  )

# =========================================================
# 5. VISUALIZATION - DAILY
# =========================================================

# 5a. Absolute Daily Gc
p_daily_abs <- ggplot(final_daily, aes(x = psiL, y = gcwater, color = species)) +
  geom_point(size = 1, alpha = 0.5) +
  scale_color_manual(values = cb_palette, name = NULL) +
  facet_wrap(~species) +
  base_theme +
  theme(
    panel.spacing.x = unit(1.5, "cm")
  ) +
  labs(
    title = expression(paste("daily absolute canopy conductance vs ", psi[paste("  ", italic(leaf))])),
    x = expression(psi[paste("  ", italic(leaf))]),
    y = expression("canopy conductance (m s"^{-1}*")")
  )

print(p_daily_abs)

# 5b. Relative Daily Gc
p_daily_rel <- ggplot(
  final_daily %>% filter(gc_rel < 1),
  aes(x = psiL, y = gc_rel, color = species)
) +
  geom_point(size = 1, alpha = 0.5) +
  scale_color_manual(values = cb_palette, name = NULL) +
  facet_wrap(~species, ncol = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  base_theme +
  theme(
    panel.spacing.x = unit(1.5, "cm")
  ) +
  labs(
    title = expression(paste("daily standardized canopy conductance vs ", psi[paste("  ", italic(leaf))])),
    x = expression(psi[paste("  ", italic(leaf))]),
    y = expression(G[c] / G[c[max]])
  )

print(p_daily_rel)

# 5c. Combined Daily (All lines in one)
p_daily_combined <- ggplot(
  final_daily %>% filter(gc_rel < 1),
  aes(x = psiL, y = gc_rel, color = species)
) +
  geom_point(size = 1, alpha = 0.5) +
  scale_color_manual(values = cb_palette, name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  base_theme +
  labs(
    title = expression(paste("species comparison: daily standardized canopy conductance vs ", psi[paste("  ", italic(leaf))])),
    x = expression(psi[paste("  ", italic(leaf))]),
    y = expression(G[c] / G[c[max]])
  )

print(p_daily_combined)

# =========================================================
# 6. VISUALIZATION - MONTHLY
# =========================================================

# 6a. Absolute Monthly Gc
p_monthly_abs <- ggplot(final_monthly, aes(x = psiL, y = gcwater, color = species)) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = cb_palette, name = NULL) +
  facet_wrap(~species) +
  base_theme +
  theme(
    panel.spacing.x = unit(1.5, "cm")
  ) +
  labs(
    title = expression(paste("monthly mean canopy conductance vs ", psi[paste("  ", italic(leaf))])),
    x = expression(psi[paste("  ", italic(leaf))]),
    y = expression("canopy conductance (m s"^{-1}*")")
  )

print(p_monthly_abs)

# 6b. Relative Monthly Gc
p_monthly_rel <- ggplot(
  final_monthly %>% filter(gc_rel < 1),
  aes(x = psiL, y = gc_rel, color = species)
) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = cb_palette, name = NULL) +
  facet_wrap(~species, ncol = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  base_theme +
  theme(
    panel.spacing.x = unit(1.5, "cm")
  ) +
  labs(
    title = expression(paste("monthly standardized canopy conductance vs ", psi[paste("  ", italic(leaf))])),
    x = expression(psi[paste("  ", italic(leaf))]),
    y = expression(G[c] / G[c[max]])
  )

print(p_monthly_rel)

# 6c. Combined Monthly (All lines in one)
p_monthly_combined <- ggplot(
  final_monthly %>% filter(gc_rel < 1),
  aes(x = psiL, y = gc_rel, color = species)
) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = cb_palette, name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  base_theme +
  labs(
    title = expression(paste("species comparison: monthly standardized canopy conductance vs ", psi[paste("  ", italic(leaf))])),
    x = expression(psi[paste("  ", italic(leaf))]),
    y = expression(G[c] / G[c[max]])
  )

print(p_monthly_combined)
# =========================================================
# 7. SAVE OUTPUTS
# =========================================================
ggsave("Figures/lpj_gcwater_psiL_filter/daily_absolute_facet.png", p_daily_abs, width = 10, height = 7)
ggsave("Figures/lpj_gcwater_psiL_filter/daily_relative_facet.png", p_daily_rel, width = 10, height = 7)
ggsave("Figures/lpj_gcwater_psiL_filter/daily_relative_comparison.png", p_daily_combined, width = 10, height = 7)

ggsave("Figures/lpj_gcwater_psiL_filter/monthly_absolute_facet.png", p_monthly_abs, width = 10, height = 7)
ggsave("Figures/lpj_gcwater_psiL_filter/monthly_relative_facet.png", p_monthly_rel, width = 10, height = 7)
ggsave("Figures/lpj_gcwater_psiL_filter/monthly_relative_comparison.png", p_monthly_combined, width = 10, height = 7)

