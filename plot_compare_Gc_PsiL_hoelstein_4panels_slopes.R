# ==============================================================================
# Gc vs PsiL plots with side-by-side Theil-Sen panels
# Split Observation Theil-Sen bars: Midday = solid, Predawn = hatched
# ==============================================================================

# ==============================================================================
# 1. SETUP & DATA LOADING
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(patchwork)
library(ggpattern)

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

# ==============================================================================
# 2. CONFIGURATIONS & AESTHETICS
# ==============================================================================

climate_txt <- paste(
  "temperature > 14°C,", "precipitation < 1 mm,",
  "global radiation > 150 W/m²,", "VPD > 0.3 kPa"
)

species_order <- c("Oak", "Beech", "Spruce", "Pine")

cb_palette <- c(
  Oak    = "#E69F00",
  Beech  = "#0072B2",
  Spruce = "#009E73",
  Pine   = "#F0E442"
)

pt_size <- 1.3

base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.position   = "bottom",
    legend.text       = element_text(size = 13),
    plot.title        = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5, size = 11),
    axis.title        = element_text(size = 14),
    axis.text         = element_text(size = 12),
    strip.text        = element_text(size = 14, face = "bold"),
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_line(color = "grey92", linewidth = 0.25),
    panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )

shape_values <- c(
  "Midday"     = 15,
  "Predawn"    = 2,
  "Simulation" = 16
)

type_patterns <- c(
  "Simulation" = "none",
  "Midday"     = "none",
  "Predawn"    = "stripe"
)

# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

minmax_standardise <- function(x) {
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  denom <- xmax - xmin
  if (!is.finite(denom) || denom == 0) return(rep(NA_real_, length(x)))
  (x - xmin) / denom
}

theil_sen <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3) return(NA_real_)
  
  slopes <- numeric(0)
  n <- length(x)
  for (i in seq_len(n - 1)) {
    dx <- x[(i + 1):n] - x[i]
    dy <- y[(i + 1):n] - y[i]
    slopes <- c(slopes, dy[dx != 0] / dx[dx != 0])
  }
  median(slopes[is.finite(slopes)], na.rm = TRUE)
}

make_four_panel_data <- function(obs_df,
                                 sim_df,
                                 obs_gc_col,
                                 obs_psi_md_col,
                                 obs_psi_pd_col,
                                 sim_gc_col,
                                 sim_psi_col) {
  obs_md <- obs_df %>%
    transmute(
      treatment,
      species,
      source  = "Observation",
      type    = "Midday",
      psi     = .data[[obs_psi_md_col]],
      gc_plot = .data[[obs_gc_col]]
    )
  
  obs_pd <- obs_df %>%
    transmute(
      treatment,
      species,
      source  = "Observation",
      type    = "Predawn",
      psi     = .data[[obs_psi_pd_col]],
      gc_plot = .data[[obs_gc_col]]
    )
  
  sim_dat <- sim_df %>%
    transmute(
      treatment,
      species,
      source  = "Simulation",
      type    = "Simulation",
      psi     = .data[[sim_psi_col]],
      gc_plot = .data[[sim_gc_col]]
    )
  
  bind_rows(sim_dat, obs_md, obs_pd) %>%
    mutate(
      species   = factor(species, levels = species_order),
      treatment = factor(treatment, levels = c("control", "drought")),
      source    = factor(source, levels = c("Simulation", "Observation")),
      type      = factor(type, levels = c("Simulation", "Midday", "Predawn"))
    ) %>%
    filter(!is.na(psi), !is.na(gc_plot))
}

make_four_panel_plot <- function(plot_dat,
                                 title,
                                 subtitle,
                                 y_lab,
                                 y_limits = NULL,
                                 y_breaks = waiver(),
                                 alpha_value = 0.7) {
  p <- ggplot(plot_dat, aes(x = psi, y = gc_plot, color = species, shape = type)) +
    geom_point(alpha = alpha_value, size = pt_size, na.rm = TRUE) +
    facet_grid(source ~ treatment) +
    scale_color_manual(values = cb_palette) +
    scale_shape_manual(name = "", values = shape_values) +
    scale_y_continuous(breaks = y_breaks) +
    labs(
      title = title,
      subtitle = subtitle,
      x = expression(Psi["leaf"] ~ "(MPa)"),
      y = y_lab,
      color = "species"
    ) +
    base_theme
  
  if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
  p
}

# Pooled observation Theil-Sen: Midday and Predawn together
calc_ts_by_panel <- function(plot_dat) {
  plot_dat %>%
    group_by(source, treatment, species) %>%
    summarise(
      ts_slope = theil_sen(psi, gc_plot),
      n = sum(complete.cases(psi, gc_plot)),
      .groups = "drop"
    ) %>%
    mutate(
      species   = factor(species, levels = species_order),
      treatment = factor(treatment, levels = c("control", "drought")),
      source    = factor(source, levels = c("Simulation", "Observation"))
    )
}

# Split observation Theil-Sen: Midday and Predawn separately
calc_ts_by_panel_split <- function(plot_dat) {
  plot_dat %>%
    group_by(source, treatment, species, type) %>%
    summarise(
      ts_slope = theil_sen(psi, gc_plot),
      n = sum(complete.cases(psi, gc_plot)),
      .groups = "drop"
    ) %>%
    mutate(
      species   = factor(species, levels = species_order),
      treatment = factor(treatment, levels = c("control", "drought")),
      source    = factor(source, levels = c("Simulation", "Observation")),
      type      = factor(type, levels = c("Simulation", "Midday", "Predawn"))
    )
}

get_ts_limits <- function(ts_dat, source_val = NULL, scale_mode = c("source", "global")) {
  scale_mode <- match.arg(scale_mode)
  
  vals <- if (scale_mode == "source") {
    ts_dat %>% filter(as.character(source) == source_val) %>% pull(ts_slope)
  } else {
    ts_dat$ts_slope
  }
  
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(c(-1, 1))
  
  rng <- range(c(vals, 0), na.rm = TRUE)
  pad <- if (diff(rng) == 0) {
    ifelse(abs(rng[2]) > 0, abs(rng[2]) * 0.15, 1)
  } else {
    diff(rng) * 0.15
  }
  
  c(rng[1] - pad, rng[2] + pad)
}

make_ts_bar_plot <- function(plot_dat, title, y_lab) {
  ts_dat <- calc_ts_by_panel(plot_dat)
  
  ggplot(ts_dat, aes(x = species, y = ts_slope, fill = species)) +
    geom_col(color = "black", width = 0.65, show.legend = FALSE) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey35") +
    facet_grid(source ~ treatment) +
    scale_fill_manual(values = cb_palette) +
    labs(
      title = title,
      subtitle = "Theil-Sen slope; observations pool Midday + Predawn",
      x = "Species",
      y = y_lab
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "none"
    )
}

make_ts_bar_plot_split <- function(plot_dat,
                                   title,
                                   y_lab,
                                   ts_scale_mode = c("source", "global")) {
  ts_scale_mode <- match.arg(ts_scale_mode)
  ts_dat <- calc_ts_by_panel_split(plot_dat)
  
  ggplot(
    ts_dat,
    aes(
      x = species,
      y = ts_slope,
      fill = species,
      pattern = type,
      group = type
    )
  ) +
    ggpattern::geom_col_pattern(
      position = position_dodge(width = 0.72),
      color = "black",
      width = 0.62,
      pattern_fill = "black",
      pattern_colour = "black",
      pattern_angle = 45,
      pattern_density = 0.35,
      pattern_spacing = 0.035,
      pattern_key_scale_factor = 0.65,
      show.legend = TRUE
    ) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey35") +
    facet_grid(
      source ~ treatment,
      scales = ifelse(ts_scale_mode == "global", "fixed", "free_y")
    ) +
    scale_fill_manual(values = cb_palette, name = "species") +
    scale_pattern_manual(values = type_patterns, name = "water potential") +
    labs(
      title = title,
      subtitle = "Observation bars split by type: Midday = solid, Predawn = hatched",
      x = "Species",
      y = y_lab
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "bottom"
    )
}

make_four_panel_with_ts_side <- function(plot_dat,
                                         title,
                                         subtitle,
                                         y_lab,
                                         ts_y_lab,
                                         y_limits = NULL,
                                         y_breaks = waiver(),
                                         alpha_value = 0.7,
                                         ts_scale_mode = c("source", "global")) {
  ts_scale_mode <- match.arg(ts_scale_mode)
  ts_dat <- calc_ts_by_panel(plot_dat)
  
  make_scatter_cell <- function(source_val, treatment_val, show_y = TRUE, show_x = TRUE) {
    p <- ggplot(
      plot_dat %>% filter(as.character(source) == source_val, as.character(treatment) == treatment_val),
      aes(x = psi, y = gc_plot, color = species, shape = type)
    ) +
      geom_point(alpha = alpha_value, size = pt_size, na.rm = TRUE) +
      scale_color_manual(values = cb_palette) +
      scale_shape_manual(name = "", values = shape_values) +
      scale_y_continuous(breaks = y_breaks) +
      labs(
        title = paste(source_val, treatment_val, sep = " - "),
        x = if (show_x) expression(Psi["leaf"] ~ "(MPa)") else NULL,
        y = if (show_y) y_lab else NULL,
        color = "species"
      ) +
      base_theme +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "bottom"
      )
    
    if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
    p
  }
  
  make_ts_cell <- function(source_val, treatment_val, show_y = TRUE, show_x = TRUE) {
    ts_lim <- get_ts_limits(ts_dat, source_val = source_val, scale_mode = ts_scale_mode)
    
    ggplot(
      ts_dat %>% filter(as.character(source) == source_val, as.character(treatment) == treatment_val),
      aes(x = species, y = ts_slope, fill = species)
    ) +
      geom_col(color = "black", width = 0.65, show.legend = FALSE) +
      geom_hline(yintercept = 0, linewidth = 0.4, color = "grey35") +
      scale_fill_manual(values = cb_palette) +
      coord_cartesian(ylim = ts_lim) +
      labs(
        title = "Theil-Sen",
        x = if (show_x) "Species" else NULL,
        y = if (show_y) ts_y_lab else NULL
      ) +
      base_theme +
      theme(
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.position = "none"
      )
  }
  
  row_sim <-
    make_scatter_cell("Simulation", "control", show_y = TRUE, show_x = FALSE) +
    make_ts_cell("Simulation", "control", show_y = TRUE, show_x = FALSE) +
    make_scatter_cell("Simulation", "drought", show_y = FALSE, show_x = FALSE) +
    make_ts_cell("Simulation", "drought", show_y = FALSE, show_x = FALSE) +
    plot_layout(widths = c(4, 1.25, 4, 1.25), guides = "collect")
  
  row_obs <-
    make_scatter_cell("Observation", "control", show_y = TRUE, show_x = TRUE) +
    make_ts_cell("Observation", "control", show_y = TRUE, show_x = TRUE) +
    make_scatter_cell("Observation", "drought", show_y = FALSE, show_x = TRUE) +
    make_ts_cell("Observation", "drought", show_y = FALSE, show_x = TRUE) +
    plot_layout(widths = c(4, 1.25, 4, 1.25), guides = "collect")
  
  row_sim / row_obs +
    plot_layout(heights = c(1, 1), guides = "collect") +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5)
      )
    ) &
    theme(legend.position = "bottom")
}

make_four_panel_with_ts_side_split <- function(plot_dat,
                                               title,
                                               subtitle,
                                               y_lab,
                                               ts_y_lab,
                                               y_limits = NULL,
                                               x_limits = NULL,
                                               y_breaks = waiver(),
                                               alpha_value = 0.7,
                                               ts_scale_mode = c("source", "global")) {
  ts_scale_mode <- match.arg(ts_scale_mode)
  ts_dat <- calc_ts_by_panel_split(plot_dat)
  
  make_scatter_cell <- function(source_val, treatment_val, show_y = TRUE, show_x = TRUE) {
    p <- ggplot(
      plot_dat %>% filter(as.character(source) == source_val, as.character(treatment) == treatment_val),
      aes(x = psi, y = gc_plot, color = species, shape = type)
    ) +
      geom_point(alpha = alpha_value, size = pt_size, na.rm = TRUE) +
      scale_color_manual(values = cb_palette) +
      scale_shape_manual(name = "", values = shape_values) +
      scale_y_continuous(breaks = y_breaks) +
      labs(
        title = paste(source_val, treatment_val, sep = " - "),
        x = if (show_x) expression(Psi["leaf"] ~ "(MPa)") else NULL,
        y = if (show_y) y_lab else NULL,
        color = "species"
      ) +
      base_theme +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "bottom"
      )
    
    if (!is.null(y_limits) & !is.null(x_limits)) p <- p + coord_cartesian(xlim = x_limits,
                                                                          ylim = y_limits)
    p
  }
  
  make_ts_cell <- function(source_val, treatment_val, show_y = TRUE, show_x = TRUE) {
    ts_lim <- get_ts_limits(ts_dat, source_val = source_val, scale_mode = ts_scale_mode)
    
    ggplot(
      ts_dat %>% filter(as.character(source) == source_val, as.character(treatment) == treatment_val),
      aes(
        x = species,
        y = ts_slope,
        fill = species,
        pattern = type,
        group = type
      )
    ) +
      ggpattern::geom_col_pattern(
        position = position_dodge(width = 0.72),
        color = "black",
        width = 0.62,
        pattern_fill = "black",
        pattern_colour = "black",
        pattern_angle = 45,
        pattern_density = 0.35,
        pattern_spacing = 0.035,
        pattern_key_scale_factor = 0.65,
        show.legend = FALSE
      ) +
      geom_hline(yintercept = 0, linewidth = 0.4, color = "grey35") +
      scale_fill_manual(values = cb_palette) +
      scale_pattern_manual(values = type_patterns) +
      coord_cartesian(ylim = ts_lim) +
      labs(
        title = "Theil-Sen",
        x = if (show_x) "Species" else NULL,
        y = if (show_y) ts_y_lab else NULL
      ) +
      base_theme +
      theme(
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.position = "none"
      )
  }
  
  row_sim <-
    make_scatter_cell("Simulation", "control", show_y = TRUE, show_x = FALSE) +
    make_ts_cell("Simulation", "control", show_y = TRUE, show_x = FALSE) +
    make_scatter_cell("Simulation", "drought", show_y = FALSE, show_x = FALSE) +
    make_ts_cell("Simulation", "drought", show_y = FALSE, show_x = FALSE) +
    plot_layout(widths = c(4, 1.4, 4, 1.4), guides = "collect")
  
  row_obs <-
    make_scatter_cell("Observation", "control", show_y = TRUE, show_x = TRUE) +
    make_ts_cell("Observation", "control", show_y = TRUE, show_x = TRUE) +
    make_scatter_cell("Observation", "drought", show_y = FALSE, show_x = TRUE) +
    make_ts_cell("Observation", "drought", show_y = FALSE, show_x = TRUE) +
    plot_layout(widths = c(4, 1.4, 4, 1.4), guides = "collect")
  
  row_sim / row_obs +
    plot_layout(heights = c(1, 1), guides = "collect") +
    plot_annotation(
      title = title,
      subtitle = paste0(subtitle, "\nObservation Theil-Sen bars: Midday = solid, Predawn = hatched"),
      theme = theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5)
      )
    ) &
    theme(legend.position = "bottom")
}

# ==============================================================================
# 4. DATA INGESTION & FORMATTING
# ==============================================================================

lpj_output_filter <- read.csv(
  "lpj_guess/lpj_guess_stem_storage/lpj_control_drought_ET_plant_ET_total_Gc_psi_leaf_psi_soil_psi_xylem_hydraulic_lag_kappy_s_min_mort_mort_cav_mort_greff_mort_min_stem_diameter_stem_rwc_twd_climate_filter.csv"
) %>%
  mutate(date = as.Date(date))

sap_flux_gc_filter <- read.csv("SCCII/sap_daily_filter.csv") %>%
  mutate(
    date = as.Date(date),
    treatment = case_when(
      treatment == "control" ~ "control",
      treatment == "treatment" ~ "drought",
      TRUE ~ treatment
    )
  )

obs_leaf_raw <- bind_rows(
  read.csv("SCCII/psiL_hoelstein_drought.csv") %>% mutate(treatment = "drought"),
  read.csv("SCCII/psiL_hoelstein_control.csv") %>% mutate(treatment = "control")
) %>%
  mutate(date = as.Date(date), month = month(date)) %>%
  filter(month >= 6, month <= 9) %>%
  rename(species = if_else("species_name" %in% names(.), "species_name", "species")) %>%
  mutate(species = factor(species, levels = species_order))

climate_filter_dates <- unique(lpj_output_filter$date)

# ==============================================================================
# 5. DATA PROCESSING & ALIGNMENT
# ==============================================================================

mod_proc <- lpj_output_filter %>%
  select(treatment, species, date, gc = Gc, psiL = psi_leaf) %>%
  mutate(
    treatment = factor(treatment, levels = c("control", "drought")),
    species = factor(species, levels = species_order)
  )

mod_proc_monthly <- mod_proc %>%
  mutate(year = year(date), month = month(date)) %>%
  group_by(species, treatment, year, month) %>%
  summarise(gc = mean(gc, na.rm = TRUE), psiL = mean(psiL, na.rm = TRUE), .groups = "drop")

obs_psi <- obs_leaf_raw %>%
  group_by(date, species, treatment) %>%
  summarise(
    psiL_md = mean(md_wp_av, na.rm = TRUE),
    psiL_pd = mean(pd_wp_av, na.rm = TRUE),
    .groups = "drop"
  )

obs_gc <- sap_flux_gc_filter %>%
  group_by(date, species, treatment) %>%
  summarise(gc_obs = mean(G_ms, na.rm = TRUE), .groups = "drop")

obs_combined <- obs_gc %>%
  left_join(obs_psi, by = c("date", "species", "treatment")) %>%
  mutate(species = factor(species, levels = species_order))

obs_filtered_climate <- obs_combined %>%
  filter(date %in% climate_filter_dates, month(date) %in% c(6, 7, 8, 9))

combined_data_lpj_obs <- mod_proc %>%
  inner_join(obs_filtered_climate, by = c("date", "species", "treatment"))

# ==============================================================================
# 6. STANDARDISATION
# ==============================================================================

combined_long <- bind_rows(
  combined_data_lpj_obs %>%
    transmute(date, species, treatment, gc = gc_obs, psiL = psiL_md,
              source = "obs", psiL_label = "midday"),
  combined_data_lpj_obs %>%
    transmute(date, species, treatment, gc = gc_obs, psiL = psiL_pd,
              source = "obs", psiL_label = "predawn"),
  combined_data_lpj_obs %>%
    transmute(date, species, treatment, gc = gc, psiL = psiL,
              source = "sim", psiL_label = "simulated")
) %>%
  filter(!is.na(psiL))

combined_std_lpj_obs <- combined_long %>%
  group_by(species, treatment, source, psiL_label) %>%
  mutate(
    gc_min = min(gc, na.rm = TRUE),
    gc_max = max(gc, na.rm = TRUE),
    gc_rel = minmax_standardise(gc)
  ) %>%
  ungroup()

cat("\n=== Standardisation range check; gc_rel should span 0-1 where possible ===\n")
combined_std_lpj_obs %>%
  group_by(species, treatment, source, psiL_label) %>%
  summarise(
    gc_abs_min = min(gc, na.rm = TRUE),
    gc_abs_max = max(gc, na.rm = TRUE),
    gc_rel_min = min(gc_rel, na.rm = TRUE),
    gc_rel_max = max(gc_rel, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  print(n = 100)

point_counts <- combined_data_lpj_obs %>%
  group_by(species, treatment) %>%
  summarise(
    obs_midday  = sum(!is.na(psiL_md) & !is.na(gc_obs)),
    obs_predawn = sum(!is.na(psiL_pd) & !is.na(gc_obs)),
    lpj_model   = sum(!is.na(psiL) & !is.na(gc)),
    .groups = "drop"
  )
print("--- Data point counts per species × treatment ---")
print(point_counts)

data_full_model_std <- mod_proc %>%
  filter(date %in% climate_filter_dates, month(date) %in% c(6, 7, 8, 9)) %>%
  group_by(species, treatment) %>%
  mutate(
    gc_min = min(gc, na.rm = TRUE),
    gc_max = max(gc, na.rm = TRUE),
    gc_rel_mod = minmax_standardise(gc)
  ) %>%
  ungroup()

data_obs_full_std_raw <- obs_combined %>%
  filter(month(date) %in% c(6, 7, 8, 9))

data_obs_full_std <- bind_rows(
  data_obs_full_std_raw %>%
    transmute(date, species, treatment, gc = gc_obs, psiL = psiL_md,
              psiL_label = "midday"),
  data_obs_full_std_raw %>%
    transmute(date, species, treatment, gc = gc_obs, psiL = psiL_pd,
              psiL_label = "predawn")
) %>%
  mutate(source = "obs") %>%
  filter(!is.na(psiL)) %>%
  group_by(species, treatment, psiL_label) %>%
  mutate(
    gc_min = min(gc, na.rm = TRUE),
    gc_max = max(gc, na.rm = TRUE),
    gc_rel = minmax_standardise(gc)
  ) %>%
  ungroup()

data_full_model_monthly_raw <- mod_proc %>%
  filter(date %in% climate_filter_dates, month(date) %in% c(6, 7, 8, 9)) %>%
  mutate(year = year(date), month = month(date)) %>%
  group_by(species, treatment, year, month) %>%
  summarise(gc = mean(gc, na.rm = TRUE), psiL = mean(psiL, na.rm = TRUE), .groups = "drop")

data_full_model_monthly_std <- data_full_model_monthly_raw %>%
  group_by(species, treatment) %>%
  mutate(
    gc_min = min(gc, na.rm = TRUE),
    gc_max = max(gc, na.rm = TRUE),
    gc_rel_mod = minmax_standardise(gc)
  ) %>%
  ungroup()

# ==============================================================================
# 7. FACETED GRID DESIGNS: ROWS = TREATMENT, COLUMNS = SPECIES
# ==============================================================================

p_gc_psi_common <- ggplot(combined_data_lpj_obs) +
  geom_point(aes(x = psiL_md, y = gc_obs), shape = 2, color = "black", alpha = 0.8, size = pt_size, na.rm = TRUE) +
  geom_point(aes(x = psiL_pd, y = gc_obs), shape = 3, color = "black", alpha = 0.8, size = pt_size, na.rm = TRUE) +
  geom_point(aes(x = psiL, y = gc, color = species), alpha = 1, size = pt_size, na.rm = TRUE) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = cb_palette) +
  coord_cartesian(ylim = c(0, 12)) +
  labs(
    title = "canopy conductance vs leaf water potential (common time)",
    subtitle = tolower(paste0(climate_txt, " | june-september\nmidday = open triangle | predawn = black + | lpj = colored")),
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] ~ (m ~ s^{-1})),
    color = ""
  ) + base_theme

p_gc_psi_full <- ggplot() +
  geom_point(data = obs_filtered_climate, aes(x = psiL_md, y = gc_obs), shape = 2, color = "black", alpha = 0.8, size = pt_size, na.rm = TRUE) +
  geom_point(data = obs_filtered_climate, aes(x = psiL_pd, y = gc_obs), shape = 3, color = "black", alpha = 0.8, size = pt_size, na.rm = TRUE) +
  geom_point(data = mod_proc_monthly, aes(x = psiL, y = gc, color = species), alpha = 0.6, size = pt_size, na.rm = TRUE) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = cb_palette) +
  coord_cartesian(ylim = c(0, 12)) +
  labs(
    title = "canopy conductance vs leaf water potential (monthly modelled means)",
    subtitle = "midday = open triangle | predawn = black + | lpj = colored monthly mean",
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] ~ (m ~ s^{-1})),
    color = ""
  ) + base_theme

p_gc_rel_common <- ggplot() +
  geom_point(data = combined_std_lpj_obs %>% filter(source == "obs", psiL_label == "midday"), aes(x = psiL, y = gc_rel), shape = 2, color = "black", alpha = 0.8, size = pt_size, na.rm = TRUE) +
  geom_point(data = combined_std_lpj_obs %>% filter(source == "obs", psiL_label == "predawn"), aes(x = psiL, y = gc_rel), shape = 3, color = "black", alpha = 0.8, size = pt_size, na.rm = TRUE) +
  geom_point(data = combined_std_lpj_obs %>% filter(source == "sim"), aes(x = psiL, y = gc_rel, color = species), alpha = 1, size = pt_size + 0.5, na.rm = TRUE) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = cb_palette) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "standardized canopy conductance vs leaf water potential (common time)",
    subtitle = "true min-max per species × treatment × source × psiL_label\nmidday = open triangle | predawn = black + | lpj = colored",
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] / G[cmax]),
    color = ""
  ) + base_theme

p_gc_rel_full <- ggplot() +
  geom_point(data = data_obs_full_std %>% filter(psiL_label == "midday"), aes(x = psiL, y = gc_rel), shape = 2, color = "black", alpha = 0.8, size = pt_size, na.rm = TRUE) +
  geom_point(data = data_obs_full_std %>% filter(psiL_label == "predawn"), aes(x = psiL, y = gc_rel), shape = 3, color = "black", alpha = 0.8, size = pt_size, na.rm = TRUE) +
  geom_point(data = data_full_model_monthly_std, aes(x = psiL, y = gc_rel_mod, color = species), alpha = 0.65, size = pt_size + 0.5, na.rm = TRUE) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = cb_palette) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "standardized canopy conductance vs leaf water potential (full monthly mean)",
    subtitle = "true min-max per species × treatment × psiL_label | lpj = colored monthly mean",
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] / G[cmax]),
    color = ""
  ) + base_theme

print(p_gc_psi_common)
print(p_gc_psi_full)
print(p_gc_rel_common)
print(p_gc_rel_full)

# ==============================================================================
# 8. SINGLE PANEL WRAP DESIGNS AND FOUR-PANEL DESIGNS
# ==============================================================================

# ------------------------------------------------------------------------------
# 8.1 Absolute values - common time
# ------------------------------------------------------------------------------
p_gc_psi_common_single <- ggplot(combined_data_lpj_obs) +
  geom_point(aes(x = psiL_md, y = gc_obs, color = species, shape = "obs midday"), alpha = 0.6, size = pt_size, na.rm = TRUE) +
  geom_point(aes(x = psiL_pd, y = gc_obs, color = species, shape = "obs predawn"), alpha = 0.6, size = pt_size, na.rm = TRUE) +
  geom_point(aes(x = psiL, y = gc, color = species, shape = "lpj simulated"), alpha = 1, size = pt_size + 0.5, na.rm = TRUE) +
  facet_wrap(vars(treatment), ncol = 2) +
  coord_cartesian(ylim = c(0, 12)) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "", values = c("obs midday" = 2, "obs predawn" = 3, "lpj simulated" = 16)) +
  labs(
    title = "canopy conductance vs leaf water potential (common time)",
    subtitle = tolower(paste0(climate_txt, "\nmidday = open triangle | predawn = + | lpj-guess-hyd = ●")),
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] ~ (m ~ s^{-1})),
    color = "species"
  ) + base_theme

plot_gc_psi_common_four <- make_four_panel_data(
  obs_df = combined_data_lpj_obs,
  sim_df = combined_data_lpj_obs,
  obs_gc_col = "gc_obs",
  obs_psi_md_col = "psiL_md",
  obs_psi_pd_col = "psiL_pd",
  sim_gc_col = "gc",
  sim_psi_col = "psiL"
)

p_gc_psi_common_single_four <- make_four_panel_plot(
  plot_gc_psi_common_four,
  title = "canopy conductance vs leaf water potential (common time)",
  subtitle = tolower(paste0(climate_txt, "\nmidday = ■ | predawn = △ | simulation = ●")),
  y_lab = expression(G[c] ~ (m ~ s^{-1})),
  y_limits = c(0, 12),
  alpha_value = 0.7
)

p_gc_psi_common_single_four_ts <- make_ts_bar_plot(
  plot_gc_psi_common_four,
  title = "Theil-Sen slopes",
  y_lab = expression("Theil-Sen slope (m s"^{-1} ~ "MPa"^{-1} ~ ")")
)

p_gc_psi_common_single_four_ts_split <- make_ts_bar_plot_split(
  plot_gc_psi_common_four,
  title = "Theil-Sen slopes split by observation type",
  y_lab = expression("Theil-Sen slope (m s"^{-1} ~ "MPa"^{-1} ~ ")"),
  ts_scale_mode = "source"
)

p_gc_psi_common_single_four_with_ts <- make_four_panel_with_ts_side(
  plot_dat = plot_gc_psi_common_four,
  title = "canopy conductance vs leaf water potential (common time)",
  subtitle = tolower(paste0(climate_txt, "\nscatter = Gc vs PsiL | right bars = pooled Theil-Sen slope per species")),
  y_lab = expression(G[c] ~ (m ~ s^{-1})),
  ts_y_lab = expression("Theil-Sen slope (m s"^{-1} ~ "MPa"^{-1} ~ ")"),
  y_limits = c(0, 12),
  alpha_value = 0.7,
  ts_scale_mode = "source"
)

p_gc_psi_common_single_four_with_ts_split <- make_four_panel_with_ts_side_split(
  plot_dat = plot_gc_psi_common_four,
  title = "canopy conductance vs leaf water potential (common time)",
  subtitle = tolower(paste0(climate_txt, "\nscatter = Gc vs PsiL | right bars = Theil-Sen slope per species")),
  y_lab = expression(G[c] ~ (m ~ s^{-1})),
  ts_y_lab = expression("Theil-Sen slope (m s"^{-1} ~ "MPa"^{-1} ~ ")"),
  y_limits = c(0, 12),
  alpha_value = 0.7,
  ts_scale_mode = "source"
)

print(p_gc_psi_common_single)
print(p_gc_psi_common_single_four)
print(p_gc_psi_common_single_four_ts)
print(p_gc_psi_common_single_four_ts_split)
print(p_gc_psi_common_single_four_with_ts)
print(p_gc_psi_common_single_four_with_ts_split)

# ------------------------------------------------------------------------------
# 8.2 Standardized values - common time
# ------------------------------------------------------------------------------
p_gc_rel_common_single <- ggplot() +
  geom_point(data = combined_std_lpj_obs %>% filter(source == "obs", psiL_label == "midday"), aes(x = psiL, y = gc_rel, color = species, shape = "obs midday"), alpha = 0.6, size = pt_size, na.rm = TRUE) +
  geom_point(data = combined_std_lpj_obs %>% filter(source == "obs", psiL_label == "predawn"), aes(x = psiL, y = gc_rel, color = species, shape = "obs predawn"), alpha = 0.6, size = pt_size, na.rm = TRUE) +
  geom_point(data = combined_std_lpj_obs %>% filter(source == "sim"), aes(x = psiL, y = gc_rel, color = species, shape = "lpj simulated"), alpha = 1, size = pt_size + 0.5, na.rm = TRUE) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "", values = c("obs midday" = 2, "obs predawn" = 3, "lpj simulated" = 16)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "standardized canopy conductance vs leaf water potential (common time)",
    subtitle = tolower(paste0(climate_txt, " (true min-max per species × treatment × source × psiL_label)\nmidday = open triangle | predawn = + | lpj-guess-hyd = ●")),
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] / G[cmax]),
    color = "species"
  ) + base_theme

plot_gc_rel_common_four <- combined_std_lpj_obs %>%
  transmute(
    treatment,
    species,
    source_raw = source,
    psiL_label,
    psi = psiL,
    gc_plot = gc_rel
  ) %>%
  mutate(
    source = case_when(
      source_raw == "sim" ~ "Simulation",
      source_raw == "obs" ~ "Observation"
    ),
    type = case_when(
      source_raw == "sim" ~ "Simulation",
      psiL_label == "midday" ~ "Midday",
      psiL_label == "predawn" ~ "Predawn"
    )
  ) %>%
  mutate(
    species = factor(species, levels = species_order),
    treatment = factor(treatment, levels = c("control", "drought")),
    source = factor(source, levels = c("Simulation", "Observation")),
    type = factor(type, levels = c("Simulation", "Midday", "Predawn"))
  ) %>%
  filter(!is.na(psi), !is.na(gc_plot), !is.na(type))

p_gc_rel_common_single_four <- make_four_panel_plot(
  plot_gc_rel_common_four,
  title = "standardized canopy conductance vs leaf water potential (common time)",
  subtitle = tolower(paste0(climate_txt, " (true min-max per species × treatment × source × psiL_label)\nmidday = ■ | predawn = △ | simulation = ●")),
  y_lab = expression(G[c] / G[cmax]),
  y_limits = c(0, 1),
  y_breaks = seq(0, 1, 0.2),
  alpha_value = 0.7
)

p_gc_rel_common_single_four_ts <- make_ts_bar_plot(
  plot_gc_rel_common_four,
  title = "Theil-Sen slopes",
  y_lab = expression("Theil-Sen slope (MPa"^{-1} ~ ")")
)

p_gc_rel_common_single_four_ts_split <- make_ts_bar_plot_split(
  plot_gc_rel_common_four,
  title = "Theil-Sen slopes split by observation type",
  y_lab = expression("Theil-Sen slope (MPa"^{-1} ~ ")"),
  ts_scale_mode = "source"
)

p_gc_rel_common_single_four_with_ts <- make_four_panel_with_ts_side(
  plot_dat = plot_gc_rel_common_four,
  title = "standardized canopy conductance vs leaf water potential (common time)",
  subtitle = tolower(paste0(climate_txt, " (true min-max per species × treatment × source × psiL_label)\nscatter = Gc/Gcmax vs PsiL | right bars = pooled Theil-Sen slope per species")),
  y_lab = expression(G[c] / G[cmax]),
  ts_y_lab = expression("Theil-Sen slope (MPa"^{-1} ~ ")"),
  y_limits = c(0, 1),
  y_breaks = seq(0, 1, 0.2),
  alpha_value = 0.7,
  ts_scale_mode = "global"
)

p_gc_rel_common_single_four_with_ts_split <- make_four_panel_with_ts_side_split(
  plot_dat = plot_gc_rel_common_four,
  title = "standardized canopy conductance vs leaf water potential (common time)",
  subtitle = tolower(paste0(climate_txt, " (true min-max per species × treatment × source × psiL_label)\nscatter = Gc/Gcmax vs PsiL | right bars = Theil-Sen slope per species")),
  y_lab = expression(G[c] / G[cmax]),
  ts_y_lab = expression("Theil-Sen slope (MPa"^{-1} ~ ")"),
  y_limits = c(0, 1),
  x_limits = c(-4, 0),
  y_breaks = seq(0, 1, 0.2),
  alpha_value = 0.7,
  ts_scale_mode = "source"
)

print(p_gc_rel_common_single)
print(p_gc_rel_common_single_four)
print(p_gc_rel_common_single_four_ts)
print(p_gc_rel_common_single_four_ts_split)
print(p_gc_rel_common_single_four_with_ts)
print(p_gc_rel_common_single_four_with_ts_split)

# ------------------------------------------------------------------------------
# 8.3 Absolute values - full monthly mean
# ------------------------------------------------------------------------------
p_gc_psi_full_single_monthly_mean <- ggplot() +
  geom_point(data = obs_filtered_climate, aes(x = psiL_md, y = gc_obs, color = species, shape = "obs midday"), alpha = 0.5, size = pt_size, na.rm = TRUE) +
  geom_point(data = obs_filtered_climate, aes(x = psiL_pd, y = gc_obs, color = species, shape = "obs predawn"), alpha = 0.5, size = pt_size, na.rm = TRUE) +
  geom_point(data = mod_proc_monthly, aes(x = psiL, y = gc, color = species, shape = "lpj monthly mean"), alpha = 0.8, size = pt_size + 0.5, na.rm = TRUE) +
  facet_wrap(vars(treatment), ncol = 2) +
  coord_cartesian(ylim = c(0, 12)) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "", values = c("obs midday" = 2, "obs predawn" = 3, "lpj monthly mean" = 16)) +
  labs(
    title = "canopy conductance vs leaf water potential (full monthly mean)",
    subtitle = "midday = open triangle | predawn = + | lpj-guess-hyd monthly mean = ●",
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] ~ (m ~ s^{-1})),
    color = "species"
  ) + base_theme

plot_gc_psi_full_monthly_four <- make_four_panel_data(
  obs_df = obs_filtered_climate,
  sim_df = mod_proc_monthly,
  obs_gc_col = "gc_obs",
  obs_psi_md_col = "psiL_md",
  obs_psi_pd_col = "psiL_pd",
  sim_gc_col = "gc",
  sim_psi_col = "psiL"
)

p_gc_psi_full_single_monthly_mean_four <- make_four_panel_plot(
  plot_gc_psi_full_monthly_four,
  title = "canopy conductance vs leaf water potential (full monthly mean)",
  subtitle = "midday = ■ | predawn = △ | lpj monthly mean = ●",
  y_lab = expression(G[c] ~ (m ~ s^{-1})),
  y_limits = c(0, 12),
  alpha_value = 0.7
)

p_gc_psi_full_single_monthly_mean_four_ts <- make_ts_bar_plot(
  plot_gc_psi_full_monthly_four,
  title = "Theil-Sen slopes",
  y_lab = expression("Theil-Sen slope (m s"^{-1} ~ "MPa"^{-1} ~ ")")
)

p_gc_psi_full_single_monthly_mean_four_ts_split <- make_ts_bar_plot_split(
  plot_gc_psi_full_monthly_four,
  title = "Theil-Sen slopes split by observation type",
  y_lab = expression("Theil-Sen slope (m s"^{-1} ~ "MPa"^{-1} ~ ")"),
  ts_scale_mode = "source"
)

p_gc_psi_full_single_monthly_mean_four_with_ts <- make_four_panel_with_ts_side(
  plot_dat = plot_gc_psi_full_monthly_four,
  title = "canopy conductance vs leaf water potential (full monthly mean)",
  subtitle = "scatter = Gc vs PsiL | right bars = pooled Theil-Sen slope per species",
  y_lab = expression(G[c] ~ (m ~ s^{-1})),
  ts_y_lab = expression("Theil-Sen slope (m s"^{-1} ~ "MPa"^{-1} ~ ")"),
  y_limits = c(0, 12),
  alpha_value = 0.7,
  ts_scale_mode = "source"
)

p_gc_psi_full_single_monthly_mean_four_with_ts_split <- make_four_panel_with_ts_side_split(
  plot_dat = plot_gc_psi_full_monthly_four,
  title = "canopy conductance vs leaf water potential (full monthly mean)",
  subtitle = "scatter = Gc vs PsiL | right bars = Theil-Sen slope per species",
  y_lab = expression(G[c] ~ (m ~ s^{-1})),
  ts_y_lab = expression("Theil-Sen slope (m s"^{-1} ~ "MPa"^{-1} ~ ")"),
  y_limits = c(0, 12),
  alpha_value = 0.7,
  ts_scale_mode = "source"
)

print(p_gc_psi_full_single_monthly_mean)
print(p_gc_psi_full_single_monthly_mean_four)
print(p_gc_psi_full_single_monthly_mean_four_ts)
print(p_gc_psi_full_single_monthly_mean_four_ts_split)
print(p_gc_psi_full_single_monthly_mean_four_with_ts)
print(p_gc_psi_full_single_monthly_mean_four_with_ts_split)

# ------------------------------------------------------------------------------
# 8.4 Standardized values - full monthly mean
# ------------------------------------------------------------------------------
p_gc_rel_full_single_monthly_mean <- ggplot() +
  geom_point(data = data_obs_full_std %>% filter(psiL_label == "midday"), aes(x = psiL, y = gc_rel, color = species, shape = "obs midday"), alpha = 0.5, size = pt_size, na.rm = TRUE) +
  geom_point(data = data_obs_full_std %>% filter(psiL_label == "predawn"), aes(x = psiL, y = gc_rel, color = species, shape = "obs predawn"), alpha = 0.5, size = pt_size, na.rm = TRUE) +
  geom_point(data = data_full_model_monthly_std, aes(x = psiL, y = gc_rel_mod, color = species, shape = "lpj monthly mean"), alpha = 0.8, size = pt_size + 0.5, na.rm = TRUE) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "", values = c("obs midday" = 2, "obs predawn" = 3, "lpj monthly mean" = 16)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "standardized canopy conductance vs leaf water potential (full monthly mean)",
    subtitle = "true min-max per species × treatment × psiL_label | lpj = colored monthly mean",
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] / G[cmax]),
    color = "species"
  ) + base_theme
print(p_gc_rel_full_single_monthly_mean)

# ------------------------------------------------------------------------------
# 8.5 Absolute values - full daily range
# ------------------------------------------------------------------------------
p_gc_psi_full_single_daily <- ggplot() +
  geom_point(data = obs_filtered_climate, aes(x = psiL_md, y = gc_obs, color = species, shape = "obs midday"), alpha = 0.5, size = pt_size, na.rm = TRUE) +
  geom_point(data = obs_filtered_climate, aes(x = psiL_pd, y = gc_obs, color = species, shape = "obs predawn"), alpha = 0.5, size = pt_size, na.rm = TRUE) +
  geom_point(data = data_full_model_std, aes(x = psiL, y = gc, color = species, shape = "lpj daily simulated"), alpha = 0.3, size = pt_size, na.rm = TRUE) +
  facet_wrap(vars(treatment), ncol = 2) +
  coord_cartesian(ylim = c(0, 12)) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "", values = c("obs midday" = 2, "obs predawn" = 3, "lpj daily simulated" = 16)) +
  labs(
    title = "canopy conductance vs leaf water potential (full daily data)",
    subtitle = "midday = open triangle | predawn = + | lpj daily simulated points = ●",
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] ~ (m ~ s^{-1})),
    color = "species"
  ) + base_theme

plot_gc_psi_full_daily_four <- make_four_panel_data(
  obs_df = obs_filtered_climate,
  sim_df = data_full_model_std,
  obs_gc_col = "gc_obs",
  obs_psi_md_col = "psiL_md",
  obs_psi_pd_col = "psiL_pd",
  sim_gc_col = "gc",
  sim_psi_col = "psiL"
)

p_gc_psi_full_single_daily_four <- make_four_panel_plot(
  plot_gc_psi_full_daily_four,
  title = "canopy conductance vs leaf water potential (full daily data)",
  subtitle = "midday = ■ | predawn = △ | lpj daily simulated = ●",
  y_lab = expression(G[c] ~ (m ~ s^{-1})),
  y_limits = c(0, 12),
  alpha_value = 0.5
)

p_gc_psi_full_single_daily_four_ts <- make_ts_bar_plot(
  plot_gc_psi_full_daily_four,
  title = "Theil-Sen slopes",
  y_lab = expression("Theil-Sen slope (m s"^{-1} ~ "MPa"^{-1} ~ ")")
)

p_gc_psi_full_single_daily_four_ts_split <- make_ts_bar_plot_split(
  plot_gc_psi_full_daily_four,
  title = "Theil-Sen slopes split by observation type",
  y_lab = expression("Theil-Sen slope (m s"^{-1} ~ "MPa"^{-1} ~ ")"),
  ts_scale_mode = "source"
)

p_gc_psi_full_single_daily_four_with_ts <- make_four_panel_with_ts_side(
  plot_dat = plot_gc_psi_full_daily_four,
  title = "canopy conductance vs leaf water potential (full daily data)",
  subtitle = "scatter = Gc vs PsiL | right bars = pooled Theil-Sen slope per species",
  y_lab = expression(G[c] ~ (m ~ s^{-1})),
  ts_y_lab = expression("Theil-Sen slope (m s"^{-1} ~ "MPa"^{-1} ~ ")"),
  y_limits = c(0, 12),
  alpha_value = 0.5,
  ts_scale_mode = "source"
)

p_gc_psi_full_single_daily_four_with_ts_split <- make_four_panel_with_ts_side_split(
  plot_dat = plot_gc_psi_full_daily_four,
  title = "canopy conductance vs leaf water potential (full daily data)",
  subtitle = "scatter = Gc vs PsiL | right bars = Theil-Sen slope per species",
  y_lab = expression(G[c] ~ (m ~ s^{-1})),
  ts_y_lab = expression("Theil-Sen slope (m s"^{-1} ~ "MPa"^{-1} ~ ")"),
  y_limits = c(0, 12),
  alpha_value = 0.5,
  ts_scale_mode = "source"
)

print(p_gc_psi_full_single_daily)
print(p_gc_psi_full_single_daily_four)
print(p_gc_psi_full_single_daily_four_ts)
print(p_gc_psi_full_single_daily_four_ts_split)
print(p_gc_psi_full_single_daily_four_with_ts)
print(p_gc_psi_full_single_daily_four_with_ts_split)

# ------------------------------------------------------------------------------
# 8.6 Standardized values - full daily range
# ------------------------------------------------------------------------------
p_gc_rel_full_single_daily <- ggplot() +
  geom_point(data = data_obs_full_std %>% filter(psiL_label == "midday"), aes(x = psiL, y = gc_rel, color = species, shape = "obs midday"), alpha = 0.5, size = pt_size, na.rm = TRUE) +
  geom_point(data = data_obs_full_std %>% filter(psiL_label == "predawn"), aes(x = psiL, y = gc_rel, color = species, shape = "obs predawn"), alpha = 0.5, size = pt_size, na.rm = TRUE) +
  geom_point(data = data_full_model_std, aes(x = psiL, y = gc_rel_mod, color = species, shape = "lpj daily simulated"), alpha = 0.3, size = pt_size, na.rm = TRUE) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "", values = c("obs midday" = 2, "obs predawn" = 3, "lpj daily simulated" = 16)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "standardized canopy conductance vs leaf water potential (full daily data)",
    subtitle = "true min-max per species × treatment × psiL_label",
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] / G[cmax]),
    color = "species"
  ) + base_theme

plot_gc_rel_full_daily_four <- bind_rows(
  data_obs_full_std %>%
    filter(psiL_label == "midday") %>%
    transmute(treatment, species, source = "Observation", type = "Midday", psi = psiL, gc_plot = gc_rel),
  data_obs_full_std %>%
    filter(psiL_label == "predawn") %>%
    transmute(treatment, species, source = "Observation", type = "Predawn", psi = psiL, gc_plot = gc_rel),
  data_full_model_std %>%
    transmute(treatment, species, source = "Simulation", type = "Simulation", psi = psiL, gc_plot = gc_rel_mod)
) %>%
  mutate(
    species   = factor(species, levels = species_order),
    treatment = factor(treatment, levels = c("control", "drought")),
    source    = factor(source, levels = c("Simulation", "Observation")),
    type      = factor(type, levels = c("Simulation", "Midday", "Predawn"))
  ) %>%
  filter(!is.na(psi), !is.na(gc_plot))

p_gc_rel_full_single_daily_four <- make_four_panel_plot(
  plot_gc_rel_full_daily_four,
  title = "standardized canopy conductance vs leaf water potential (full daily data)",
  subtitle = "true min-max per species × treatment × psiL_label\nmidday = ■ | predawn = △ | lpj daily simulated = ●",
  y_lab = expression(G[c] / G[cmax]),
  y_limits = c(0, 1),
  y_breaks = seq(0, 1, 0.2),
  alpha_value = 0.5
)

p_gc_rel_full_single_daily_four_ts <- make_ts_bar_plot(
  plot_gc_rel_full_daily_four,
  title = "Theil-Sen slopes for standardized Gc/Gcmax vs PsiL (full daily data)",
  y_lab = expression("Theil-Sen slope (MPa"^{-1} ~ ")")
)

p_gc_rel_full_single_daily_four_ts_split <- make_ts_bar_plot_split(
  plot_gc_rel_full_daily_four,
  title = "Theil-Sen slopes split by observation type",
  y_lab = expression("Theil-Sen slope (MPa"^{-1} ~ ")"),
  ts_scale_mode = "global"
)

p_gc_rel_full_single_daily_four_with_ts <- make_four_panel_with_ts_side(
  plot_dat = plot_gc_rel_full_daily_four,
  title = "standardized canopy conductance vs leaf water potential (full daily data)",
  subtitle = "scatter = Gc/Gcmax vs PsiL | right bars = pooled Theil-Sen slope per species",
  y_lab = expression(G[c] / G[cmax]),
  ts_y_lab = expression("Theil-Sen slope (MPa"^{-1} ~ ")"),
  y_limits = c(0, 1),
  y_breaks = seq(0, 1, 0.2),
  alpha_value = 0.5,
  ts_scale_mode = "global"
)

p_gc_rel_full_single_daily_four_with_ts_split <- make_four_panel_with_ts_side_split(
  plot_dat = plot_gc_rel_full_daily_four,
  title = "standardized canopy conductance vs leaf water potential (full daily data)",
  subtitle = "scatter = Gc/Gcmax vs PsiL | right bars = Theil-Sen slope per species",
  y_lab = expression(G[c] / G[cmax]),
  ts_y_lab = expression("Theil-Sen slope (MPa"^{-1} ~ ")"),
  y_limits = c(0, 1),
  y_breaks = seq(0, 1, 0.2),
  alpha_value = 0.5,
  ts_scale_mode = "global"
)

print(p_gc_rel_full_single_daily)
print(p_gc_rel_full_single_daily_four)
print(p_gc_rel_full_single_daily_four_ts)
print(p_gc_rel_full_single_daily_four_ts_split)
print(p_gc_rel_full_single_daily_four_with_ts)
print(p_gc_rel_full_single_daily_four_with_ts_split)

# ==============================================================================
# 9. ORIGINAL SIMULATION-ONLY THEIL-SEN SLOPE FIGURES
# ==============================================================================

ts_abs <- combined_data_lpj_obs %>%
  group_by(species, treatment) %>%
  summarise(ts_slope = theil_sen(psiL, gc), n = sum(!is.na(psiL) & !is.na(gc)), .groups = "drop")

print("=== Theil-Sen slopes - absolute Gc (simulated, common time) ===")
print(as.data.frame(ts_abs))

p_ts_abs <- ggplot(ts_abs, aes(x = species, y = ts_slope, fill = species)) +
  geom_col(color = "black", width = 0.6) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_fill_manual(values = cb_palette, guide = "none") +
  geom_text(aes(label = round(ts_slope, 3)), vjust = -0.4, size = 3.5) +
  labs(
    title = "Theil-Sen slope - absolute Gc vs PsiL (LPJ simulated, common time)",
    subtitle = "median of pairwise slopes",
    x = "Species",
    y = expression("Theil-Sen slope (m s"^{-1} ~ "MPa"^{-1} ~ ")")
  ) + base_theme

print(p_ts_abs)

ts_rel <- combined_std_lpj_obs %>%
  filter(source == "sim") %>%
  group_by(species, treatment) %>%
  summarise(ts_slope = theil_sen(psiL, gc_rel), n = n(), .groups = "drop")

print("=== Theil-Sen slopes - relative Gc/Gcmax (simulated, common time) ===")
print(as.data.frame(ts_rel))

p_ts_rel <- ggplot(ts_rel, aes(x = species, y = ts_slope, fill = species)) +
  geom_col(color = "black", width = 0.6) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.4) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_fill_manual(values = cb_palette, guide = "none") +
  geom_text(aes(label = round(ts_slope, 3)), vjust = -0.4, size = 3.5) +
  labs(
    title = "Theil-Sen slope - relative Gc/Gcmax vs PsiL (LPJ simulated, common time)",
    subtitle = "median of pairwise slopes | true min-max standardized",
    x = "Species",
    y = expression("Theil-Sen slope (MPa"^{-1} ~ ")")
  ) + base_theme

print(p_ts_rel)

# ==============================================================================
# 10. BINNED MODELLED Gc vs PsiL
# ==============================================================================

psiL_bin_width <- 0.2

mod_binned <- data_full_model_std %>%
  mutate(psiL_bin = round(psiL / psiL_bin_width) * psiL_bin_width) %>%
  group_by(treatment, species, psiL_bin) %>%
  summarise(
    gc_mean = mean(gc, na.rm = TRUE),
    gc_sd = sd(gc, na.rm = TRUE),
    gc_se = sd(gc, na.rm = TRUE) / sqrt(n()),
    psiL_bin_mean = mean(psiL, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 3)

cat("\n=== binned modelled Gc vs PsiL: bin width =", psiL_bin_width, "MPa ===\n")
cat("total bins:", nrow(mod_binned), "\n")
print(mod_binned %>% group_by(treatment, species) %>% summarise(bins = n(), .groups = "drop"))

p_binned_mod_gc_psiL <- ggplot(mod_binned, aes(x = psiL_bin, y = gc_mean, color = species)) +
  geom_errorbar(aes(ymin = gc_mean - gc_se, ymax = gc_mean + gc_se), width = 0.03, linewidth = 0.4, alpha = 0.5) +
  geom_point(size = 2.2, na.rm = TRUE) +
  geom_line(linewidth = 0.7, alpha = 0.7, na.rm = TRUE) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = cb_palette) +
  labs(
    title = "binned modelled canopy conductance vs leaf water potential (LPJ only)",
    subtitle = paste0("psiL binned at ", psiL_bin_width, " MPa intervals | points = bin-mean Gc | error bars = ±1 SE | line connects bins"),
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] ~ (m ~ s^{-1})),
    color = ""
  ) + base_theme +
  theme(legend.position = "none")

p_binned_mod_gc_psiL_single <- ggplot(mod_binned, aes(x = psiL_bin, y = gc_mean, color = species)) +
  geom_errorbar(aes(ymin = gc_mean - gc_se, ymax = gc_mean + gc_se), width = 0.03, linewidth = 0.4, alpha = 0.5) +
  geom_point(size = 2.2, na.rm = TRUE) +
  geom_line(linewidth = 0.7, alpha = 0.7, na.rm = TRUE) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_color_manual(values = cb_palette) +
  labs(
    title = "binned modelled canopy conductance vs leaf water potential (LPJ only)",
    subtitle = paste0("psiL binned at ", psiL_bin_width, " MPa intervals | points = bin-mean Gc | error bars = ±1 SE | line connects bins"),
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] ~ (m ~ s^{-1})),
    color = "species"
  ) + base_theme

mod_binned_std <- mod_binned %>%
  group_by(species) %>%
  mutate(
    gc_bin_min = min(gc_mean, na.rm = TRUE),
    gc_bin_max = max(gc_mean, na.rm = TRUE),
    gc_std = (gc_mean - gc_bin_min) / (gc_bin_max - gc_bin_min),
    gc_std = pmin(pmax(gc_std, 0), 1)
  ) %>%
  ungroup()

p_binned_std <- ggplot(mod_binned_std, aes(x = psiL_bin, y = gc_std, color = species)) +
  geom_errorbar(aes(
    ymin = pmax(gc_std - gc_se / (gc_bin_max - gc_bin_min), 0),
    ymax = pmin(gc_std + gc_se / (gc_bin_max - gc_bin_min), 1)
  ), width = 0.03, linewidth = 0.4, alpha = 0.5) +
  geom_point(size = 2.2, na.rm = TRUE) +
  geom_line(linewidth = 0.7, alpha = 0.7, na.rm = TRUE) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = cb_palette) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "standardized binned Gc vs PsiL (LPJ only, per-species standardization)",
    subtitle = paste0("psiL binned at ", psiL_bin_width, " MPa | Gc standardized 0-1 within each species | error bars = ±1 SE (scaled)"),
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] ~ standardized ~ (0 - 1 ~ per ~ species)),
    color = ""
  ) + base_theme +
  theme(legend.position = "none")

p_binned_std_single <- ggplot(mod_binned_std, aes(x = psiL_bin, y = gc_std, color = species)) +
  geom_errorbar(aes(
    ymin = pmax(gc_std - gc_se / (gc_bin_max - gc_bin_min), 0),
    ymax = pmin(gc_std + gc_se / (gc_bin_max - gc_bin_min), 1)
  ), width = 0.03, linewidth = 0.4, alpha = 0.5) +
  geom_point(size = 2.2, na.rm = TRUE) +
  geom_line(linewidth = 0.7, alpha = 0.7, na.rm = TRUE) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "standardized binned Gc vs PsiL (LPJ only, per-species standardization)",
    subtitle = paste0("psiL binned at ", psiL_bin_width, " MPa | Gc standardized 0-1 within each species | error bars = ±1 SE (scaled)"),
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] ~ standardized ~ (0 - 1 ~ per ~ species)),
    color = "species"
  ) + base_theme

lpj_binned_rel <- data_full_model_std %>%
  filter(gc_rel_mod > 0) %>%
  mutate(psiL_bin = round(psiL / psiL_bin_width) * psiL_bin_width) %>%
  group_by(treatment, species, psiL_bin) %>%
  summarise(
    gc_rel_mean = mean(gc_rel_mod, na.rm = TRUE),
    gc_rel_se = sd(gc_rel_mod, na.rm = TRUE) / sqrt(n()),
    psiL_mean = mean(psiL, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 3)

p_binned_rel_single <- ggplot() +
  geom_point(data = data_obs_full_std %>% filter(psiL_label == "midday"), aes(x = psiL, y = gc_rel, color = species, shape = "obs midday"), alpha = 0.35, size = pt_size, na.rm = TRUE) +
  geom_point(data = data_obs_full_std %>% filter(psiL_label == "predawn"), aes(x = psiL, y = gc_rel, color = species, shape = "obs predawn"), alpha = 0.35, size = pt_size, na.rm = TRUE) +
  geom_errorbar(data = lpj_binned_rel, aes(x = psiL_bin, y = gc_rel_mean, color = species, ymin = pmax(gc_rel_mean - gc_rel_se, 0), ymax = pmin(gc_rel_mean + gc_rel_se, 1)), width = 0.04, linewidth = 0.4, alpha = 0.5) +
  geom_point(data = lpj_binned_rel, aes(x = psiL_bin, y = gc_rel_mean, color = species, shape = "lpj binned"), size = 2.5, na.rm = TRUE) +
  geom_line(data = lpj_binned_rel, aes(x = psiL_bin, y = gc_rel_mean, color = species), linewidth = 0.8, alpha = 0.8, na.rm = TRUE) +
  facet_wrap(vars(treatment), ncol = 2) +
  scale_color_manual(values = cb_palette) +
  scale_shape_manual(name = "", values = c("obs midday" = 2, "obs predawn" = 3, "lpj binned" = 16)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "standardized Gc vs PsiL - binned LPJ + raw observed (full daily)",
    subtitle = paste0("LPJ psiL binned at ", psiL_bin_width, " MPa | binned points = mean Gc/Gcmax | error bars = ±1 SE"),
    x = expression(Psi["leaf"] ~ "(MPa)"),
    y = expression(G[c] / G[cmax]),
    color = "species"
  ) + base_theme

p_binned_rel_grid <- p_binned_rel_single +
  facet_grid(treatment ~ species) +
  labs(color = "") +
  theme(legend.position = "none")

print(p_binned_mod_gc_psiL)
print(p_binned_mod_gc_psiL_single)
print(p_binned_std)
print(p_binned_std_single)
print(p_binned_rel_single)
print(p_binned_rel_grid)

# ==============================================================================
# 11. CURVE-FITTING / PSI50 SECTION
# ==============================================================================
# Keep your nonlinear curve-fitting section here if needed.

# ==============================================================================
# 12. EXPORT
# ==============================================================================

out_dir <- "Figures/lpj_guess_stem_storage/Gc_PsiL"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "compare_Gc_vs_PsiL_common_time.png"), p_gc_psi_common, width = 13, height = 9, dpi = 300)
ggsave(file.path(out_dir, "compare_Gc_vs_PsiL_full.png"), p_gc_psi_full, width = 13, height = 9, dpi = 300)
ggsave(file.path(out_dir, "compare_Gc_rel_vs_PsiL_common_time.png"), p_gc_rel_common, width = 13, height = 9, dpi = 300)
ggsave(file.path(out_dir, "compare_Gc_rel_vs_PsiL_full.png"), p_gc_rel_full, width = 13, height = 9, dpi = 300)

ggsave(file.path(out_dir, "single_panel_Gc_vs_PsiL_common.png"), p_gc_psi_common_single, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "single_panel_Gc_rel_vs_PsiL_common.png"), p_gc_rel_common_single, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "single_panel_Gc_vs_PsiL_full_monthly.png"), p_gc_psi_full_single_monthly_mean, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "single_panel_Gc_rel_vs_PsiL_full_monthly.png"), p_gc_rel_full_single_monthly_mean, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "single_panel_Gc_vs_PsiL_full_daily.png"), p_gc_psi_full_single_daily, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "single_panel_Gc_rel_vs_PsiL_full_daily.png"), p_gc_rel_full_single_daily, width = 14, height = 8, dpi = 300)

# Four-panel scatter-only exports
ggsave(file.path(out_dir, "four_panel_Gc_vs_PsiL_common.png"), p_gc_psi_common_single_four, width = 14, height = 9, dpi = 300)
ggsave(file.path(out_dir, "four_panel_Gc_rel_vs_PsiL_common.png"), p_gc_rel_common_single_four, width = 14, height = 9, dpi = 300)
ggsave(file.path(out_dir, "four_panel_Gc_vs_PsiL_full_monthly.png"), p_gc_psi_full_single_monthly_mean_four, width = 14, height = 9, dpi = 300)
ggsave(file.path(out_dir, "four_panel_Gc_vs_PsiL_full_daily.png"), p_gc_psi_full_single_daily_four, width = 14, height = 9, dpi = 300)
ggsave(file.path(out_dir, "four_panel_Gc_rel_vs_PsiL_full_daily.png"), p_gc_rel_full_single_daily_four, width = 14, height = 9, dpi = 300)

# Theil-Sen bar-only exports, pooled and split
ggsave(file.path(out_dir, "ts_four_panel_Gc_vs_PsiL_common.png"), p_gc_psi_common_single_four_ts, width = 14, height = 7, dpi = 300)
ggsave(file.path(out_dir, "ts_four_panel_Gc_rel_vs_PsiL_common.png"), p_gc_rel_common_single_four_ts, width = 14, height = 7, dpi = 300)
ggsave(file.path(out_dir, "ts_four_panel_Gc_vs_PsiL_full_monthly.png"), p_gc_psi_full_single_monthly_mean_four_ts, width = 14, height = 7, dpi = 300)
ggsave(file.path(out_dir, "ts_four_panel_Gc_vs_PsiL_full_daily.png"), p_gc_psi_full_single_daily_four_ts, width = 14, height = 7, dpi = 300)
ggsave(file.path(out_dir, "ts_four_panel_Gc_rel_vs_PsiL_full_daily.png"), p_gc_rel_full_single_daily_four_ts, width = 14, height = 7, dpi = 300)

ggsave(file.path(out_dir, "ts_split_four_panel_Gc_vs_PsiL_common.png"), p_gc_psi_common_single_four_ts_split, width = 14, height = 7, dpi = 300)
ggsave(file.path(out_dir, "ts_split_four_panel_Gc_rel_vs_PsiL_common.png"), p_gc_rel_common_single_four_ts_split, width = 14, height = 7, dpi = 300)
ggsave(file.path(out_dir, "ts_split_four_panel_Gc_vs_PsiL_full_monthly.png"), p_gc_psi_full_single_monthly_mean_four_ts_split, width = 14, height = 7, dpi = 300)
ggsave(file.path(out_dir, "ts_split_four_panel_Gc_vs_PsiL_full_daily.png"), p_gc_psi_full_single_daily_four_ts_split, width = 14, height = 7, dpi = 300)
ggsave(file.path(out_dir, "ts_split_four_panel_Gc_rel_vs_PsiL_full_daily.png"), p_gc_rel_full_single_daily_four_ts_split, width = 14, height = 7, dpi = 300)

# Combined scatter + Theil-Sen side-by-side exports, pooled and split
ggsave(file.path(out_dir, "four_panel_Gc_vs_PsiL_common_with_TS.png"), p_gc_psi_common_single_four_with_ts, width = 18, height = 10, dpi = 300)
ggsave(file.path(out_dir, "four_panel_Gc_rel_vs_PsiL_common_with_TS.png"), p_gc_rel_common_single_four_with_ts, width = 18, height = 10, dpi = 300)
ggsave(file.path(out_dir, "four_panel_Gc_vs_PsiL_full_monthly_with_TS.png"), p_gc_psi_full_single_monthly_mean_four_with_ts, width = 18, height = 10, dpi = 300)
ggsave(file.path(out_dir, "four_panel_Gc_vs_PsiL_full_daily_with_TS.png"), p_gc_psi_full_single_daily_four_with_ts, width = 18, height = 10, dpi = 300)
ggsave(file.path(out_dir, "four_panel_Gc_rel_vs_PsiL_full_daily_with_TS.png"), p_gc_rel_full_single_daily_four_with_ts, width = 18, height = 10, dpi = 300)

ggsave(file.path(out_dir, "four_panel_Gc_vs_PsiL_common_with_TS_split.png"), p_gc_psi_common_single_four_with_ts_split, width = 18, height = 10, dpi = 300)
ggsave(file.path(out_dir, "four_panel_Gc_rel_vs_PsiL_common_with_TS_split.png"), p_gc_rel_common_single_four_with_ts_split, width = 18, height = 10, dpi = 300)
ggsave(file.path(out_dir, "four_panel_Gc_vs_PsiL_full_monthly_with_TS_split.png"), p_gc_psi_full_single_monthly_mean_four_with_ts_split, width = 18, height = 10, dpi = 300)
ggsave(file.path(out_dir, "four_panel_Gc_vs_PsiL_full_daily_with_TS_split.png"), p_gc_psi_full_single_daily_four_with_ts_split, width = 18, height = 10, dpi = 300)
ggsave(file.path(out_dir, "four_panel_Gc_rel_vs_PsiL_full_daily_with_TS_split.png"), p_gc_rel_full_single_daily_four_with_ts_split, width = 18, height = 10, dpi = 300)

# Other original exports
ggsave(file.path(out_dir, "theil_sen_absolute.png"), p_ts_abs, width = 10, height = 8, dpi = 300)
ggsave(file.path(out_dir, "theil_sen_relative.png"), p_ts_rel, width = 10, height = 8, dpi = 300)
ggsave(file.path(out_dir, "binned_modelled_Gc_vs_PsiL.png"), p_binned_mod_gc_psiL, width = 13, height = 9, dpi = 300)
ggsave(file.path(out_dir, "binned_modelled_Gc_vs_PsiL_single.png"), p_binned_mod_gc_psiL_single, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "binned_std_Gc_vs_PsiL.png"), p_binned_std, width = 13, height = 9, dpi = 300)
ggsave(file.path(out_dir, "binned_std_Gc_vs_PsiL_single.png"), p_binned_std_single, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "binned_rel_Gc_vs_PsiL_single.png"), p_binned_rel_single, width = 14, height = 8, dpi = 300)
ggsave(file.path(out_dir, "binned_rel_Gc_vs_PsiL_grid.png"), p_binned_rel_grid, width = 13, height = 9, dpi = 300)

cat("\nAll plots exported to:", out_dir, "\n")
