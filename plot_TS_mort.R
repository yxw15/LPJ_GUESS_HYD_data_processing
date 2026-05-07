# ============================================================
# Annual time series from the LAST DAY of each year
# Mortality-related variables + kappa variables
#
# Data location:
#   results/species_folder/*.out
#
# Figure output location:
#   Figures/annual_last_day_mort_kappa/
#
# Outputs:
#   1) all species together:
#      annual_last_day_plots$mort
#      annual_last_day_plots$kappa_s_min
#
#   2) one plot per species:
#      annual_last_day_plots_by_species$mort$Beech
#      annual_last_day_plots_by_species$kappa_s_today$Spruce
# ============================================================

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

suppressPackageStartupMessages({
  library(tidyverse)
})

# ---------------------------
# Settings
# ---------------------------
results_dir <- file.path("results")
out_dir     <- file.path("Figures", "annual_last_day_mort_kappa")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

year_min <- 1991
year_max <- 2022

plot_year_min <- 1991
plot_year_max <- 2022

species_vec <- c("Beech", "Oak_pub", "Oak_rob", "Pine", "Spruce")

cb_palette <- c(
  Oak_pub = "darkorange",
  Oak_rob = "#F0E442",
  Beech   = "dodgerblue",
  Spruce  = "green4",
  Pine    = "purple1"
)

vars_to_plot <- tribble(
  ~file,               ~var_label,                          ~ylab,                               ~stem,
  "mort.out",          "Total mortality",                   "Total mortality",                   "mort",
  "mort_cav_day.out",  "Hydraulic-failure mortality",       "Hydraulic mortality",               "mort_cav_day",
  "mort_min.out",      "Background mortality",              "Background mortality",              "mort_min",
  "mort_greff.out",    "Growth-efficiency mortality",       "Growth-efficiency mortality",       "mort_greff",
  "kappa_s_min.out",   "Maximum value of cavitated xylems", "Maximum value of cavitated xylems", "kappa_s_min",
  "kappa_s_today.out", "Fraction of cavitated xylem",       "Fraction of cavitated xylem",       "kappa_s_today"
)

# ============================================================
# Read one file for one species
# Assumes columns 3 and 4 are: year, day
# Uses the last numeric column after excluding year/day as value
# ============================================================
read_species_file <- function(filepath, species_name, year_min, year_max) {
  
  df <- read.table(filepath, header = TRUE, check.names = FALSE)
  
  names(df)[3:4] <- c("year", "day")
  
  value_cols <- names(df)[sapply(df, is.numeric)]
  value_cols <- setdiff(value_cols, c("year", "day"))
  
  if (length(value_cols) == 0) {
    stop(paste("No numeric value column found in:", filepath))
  }
  
  value_col <- tail(value_cols, 1)
  message("Using value column for ", species_name, " in ", basename(filepath), ": ", value_col)
  
  df %>%
    select(year, day, value = all_of(value_col)) %>%
    mutate(species = species_name) %>%
    filter(year >= year_min, year <= year_max)
}

# ============================================================
# Read the same variable across all species folders
# ============================================================
read_all_species <- function(filename, species_vec, year_min, year_max, results_dir) {
  
  out_list <- list()
  
  for (sp in species_vec) {
    filepath <- file.path(results_dir, sp, filename)
    
    if (!file.exists(filepath)) {
      warning("File not found: ", filepath)
      next
    }
    
    out_list[[sp]] <- read_species_file(
      filepath     = filepath,
      species_name = sp,
      year_min     = year_min,
      year_max     = year_max
    )
  }
  
  bind_rows(out_list)
}

# ============================================================
# Keep ONLY the last day of each year for each species
# ============================================================
keep_last_day_each_year <- function(df) {
  df %>%
    group_by(species, year) %>%
    arrange(day, .by_group = TRUE) %>%
    slice_tail(n = 1) %>%
    ungroup()
}

# ============================================================
# Plot annual time series: all species together
# ============================================================
plot_last_day_ts <- function(df, title_text, ylab_text) {
  ggplot(df, aes(x = year, y = value, color = species)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_color_manual(values = cb_palette) +
    scale_x_continuous(
      breaks = seq(plot_year_min, plot_year_max, by = 2),
      limits = c(plot_year_min, plot_year_max)
    ) +
    labs(
      title = title_text,
      x = "Year",
      y = ylab_text,
      color = "Species"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

# ============================================================
# Plot annual time series: one species only
# ============================================================
plot_last_day_ts_one_species <- function(df, species_name, title_text, ylab_text) {
  
  df_sp <- df %>% filter(species == species_name)
  
  ggplot(df_sp, aes(x = year, y = value)) +
    geom_line(linewidth = 0.9, color = cb_palette[[species_name]]) +
    geom_point(size = 2, color = cb_palette[[species_name]]) +
    scale_x_continuous(
      breaks = seq(plot_year_min, plot_year_max, by = 2),
      limits = c(plot_year_min, plot_year_max)
    ) +
    labs(
      title = paste0(title_text, " - ", species_name),
      x = "Year",
      y = ylab_text
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

# ============================================================
# Run
# ============================================================
annual_last_day_data             <- list()
annual_last_day_plots            <- list()
annual_last_day_plots_by_species <- list()

for (i in seq_len(nrow(vars_to_plot))) {
  
  this_file  <- vars_to_plot$file[i]
  this_label <- vars_to_plot$var_label[i]
  this_ylab  <- vars_to_plot$ylab[i]
  this_stem  <- vars_to_plot$stem[i]
  
  message("Processing: ", this_file)
  
  df_raw <- read_all_species(
    filename    = this_file,
    species_vec = species_vec,
    year_min    = year_min,
    year_max    = year_max,
    results_dir = results_dir
  )
  
  if (nrow(df_raw) == 0) {
    warning("No data read for file: ", this_file)
    next
  }
  
  df_last <- df_raw %>%
    keep_last_day_each_year() %>%
    filter(year >= plot_year_min, year <= plot_year_max)
  
  annual_last_day_data[[this_stem]] <- df_last
  
  # ---------------------------
  # All species in one figure
  # ---------------------------
  p_all <- plot_last_day_ts(
    df         = df_last,
    title_text = paste0(this_label, " (last day of each year)"),
    ylab_text  = this_ylab
  )
  
  annual_last_day_plots[[this_stem]] <- p_all
  
  ggsave(
    filename = file.path(out_dir, paste0(this_stem, "_annual_last_day_all_species.png")),
    plot     = p_all,
    width    = 9,
    height   = 5.5,
    dpi      = 300
  )
  
  # ---------------------------
  # One figure per species
  # ---------------------------
  annual_last_day_plots_by_species[[this_stem]] <- list()
  
  for (sp in species_vec) {
    
    if (!sp %in% unique(df_last$species)) {
      warning("Species not found in data for ", this_file, ": ", sp)
      next
    }
    
    p_sp <- plot_last_day_ts_one_species(
      df           = df_last,
      species_name = sp,
      title_text   = paste0(this_label, " (last day of each year)"),
      ylab_text    = this_ylab
    )
    
    annual_last_day_plots_by_species[[this_stem]][[sp]] <- p_sp
    
    ggsave(
      filename = file.path(out_dir, paste0(this_stem, "_annual_last_day_", sp, ".png")),
      plot     = p_sp,
      width    = 8,
      height   = 5,
      dpi      = 300
    )
  }
}

# ============================================================
# Show plots
# ============================================================

# all species together
annual_last_day_plots$mort
annual_last_day_plots$mort_cav_day
annual_last_day_plots$mort_min
annual_last_day_plots$mort_greff
annual_last_day_plots$kappa_s_min
annual_last_day_plots$kappa_s_today

# one species only
annual_last_day_plots_by_species$mort$Beech
annual_last_day_plots_by_species$mort$Oak_pub
annual_last_day_plots_by_species$mort$Oak_rob
annual_last_day_plots_by_species$mort$Pine
annual_last_day_plots_by_species$mort$Spruce

annual_last_day_plots_by_species$kappa_s_min$Beech
annual_last_day_plots_by_species$kappa_s_min$Spruce
annual_last_day_plots_by_species$kappa_s_today$Beech
annual_last_day_plots_by_species$kappa_s_today$Spruce