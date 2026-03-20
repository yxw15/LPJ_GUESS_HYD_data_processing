# ============================================================
# LPJ-GUESS post-processing
# Relationship plots ONLY:
#   1) gc      (leaf conductance)   vs psi_leaf
#   2) gcwater (canopy conductance) vs psi_leaf
#
# Outputs for each timescale and season selection:
#   - combined ALL species
#   - facet by species
#
# Timescales:
#   - daily raw
#   - monthly mean
#   - yearly mean
#
# Season selections:
#   - whole year
#   - July-August only
#
# Notes:
#   - facet plots use free x and free y scales
#   - smooth lines keep species color
# ============================================================
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")


suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
})

# ---------------------------
# Settings
# ---------------------------

out_dir_whole  <- file.path(base_dir, "Figures", "gc_gcwater_vs_psileaf_WholeYear")
out_dir_julaug <- file.path(base_dir, "Figures", "gc_gcwater_vs_psileaf_JulAug")

dir.create(out_dir_whole,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_julaug, recursive = TRUE, showWarnings = FALSE)

day_starts_at_zero <- TRUE

# Optional: restrict years
year_min <- NULL
year_max <- NULL

# Binning settings
bin_width_mpa <- 0.05
min_n_bin_daily   <- 0
min_n_bin_monthly <- 0
min_n_bin_yearly  <- 0

# ---------------------------
# Species mapping + colors
# ---------------------------
species_map <- tribble(
  ~species,   ~colname,
  "Beech",    "Fag_syl",
  "Oak_pub",  "Que_pub",
  "Oak_rob",  "Que_rob",
  "Pine",     "Pin_syl",
  "Spruce",   "Pic_abi"
)

cb_palette <- c(
  Oak_pub  = "darkorange",
  Oak_rob  = "#F0E442",
  Beech    = "dodgerblue",
  Spruce   = "green4",
  Pine     = "purple1"
)

# ---------------------------
# Helpers
# ---------------------------
as_guess_date <- function(Year, Day, day_starts_at_zero = TRUE) {
  as.Date(if (day_starts_at_zero) Day else Day - 1, origin = paste0(Year, "-01-01"))
}

theme_clean <- function() {
  theme_minimal() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      strip.text       = element_text(face = "bold", size = 12)
    )
}

pick_numeric_col <- function(df, wanted, file_for_msg = "") {
  if (!is.null(wanted) && wanted %in% names(df)) return(wanted)
  candidates <- setdiff(names(df), c("Year", "Day", "Lon", "Lat"))
  numeric_candidates <- candidates[sapply(df[candidates], is.numeric)]
  if (length(numeric_candidates) == 0) {
    stop("No numeric column found in ", file_for_msg, ". Available: ", paste(names(df), collapse = ", "))
  }
  message("Note: using numeric column '", numeric_candidates[1], "' from ", basename(file_for_msg))
  numeric_candidates[1]
}

read_daily_value <- function(species, colname, filename, value_name) {
  f <- file.path(base_dir, "results", species, filename)
  if (!file.exists(f)) {
    warning("Missing file: ", f)
    return(NULL)
  }
  
  df <- read.table(f, header = TRUE)
  
  if (!all(c("Year", "Day") %in% names(df))) {
    stop("File missing Year/Day columns: ", f)
  }
  
  col_use <- pick_numeric_col(df, colname, f)
  
  df %>%
    mutate(Date = as_guess_date(Year, Day, day_starts_at_zero)) %>%
    transmute(Date, Species = species, !!value_name := .data[[col_use]])
}

read_pair_psileaf_y <- function(species, colname, y_file, y_name) {
  psi <- read_daily_value(species, colname, "dpsileaf.out", "psi_leaf")
  yy  <- read_daily_value(species, colname, y_file, y_name)
  
  if (is.null(psi) || is.null(yy)) return(NULL)
  
  inner_join(psi, yy, by = c("Date", "Species"))
}

apply_time_filter <- function(df, months_keep = NULL) {
  if (is.null(df)) return(NULL)
  
  if (!is.null(year_min)) {
    df <- df %>% filter(lubridate::year(Date) >= year_min)
  }
  
  if (!is.null(year_max)) {
    df <- df %>% filter(lubridate::year(Date) <= year_max)
  }
  
  if (!is.null(months_keep)) {
    df <- df %>% filter(lubridate::month(Date) %in% months_keep)
  }
  
  df
}

# ---------------------------
# Aggregation
# ---------------------------
aggregate_monthly_mean <- function(df, y_col) {
  df %>%
    mutate(
      Year  = lubridate::year(Date),
      Month = lubridate::month(Date)
    ) %>%
    group_by(Species, Year, Month) %>%
    summarise(
      psi_leaf = mean(psi_leaf, na.rm = TRUE),
      value    = mean(.data[[y_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(!!y_col := value) %>%
    mutate(Species = factor(Species, levels = species_map$species))
}

aggregate_yearly_mean <- function(df, y_col) {
  df %>%
    mutate(Year = lubridate::year(Date)) %>%
    group_by(Species, Year) %>%
    summarise(
      psi_leaf = mean(psi_leaf, na.rm = TRUE),
      value    = mean(.data[[y_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(!!y_col := value) %>%
    mutate(Species = factor(Species, levels = species_map$species))
}

# ---------------------------
# Binning
# ---------------------------
bin_by_x <- function(df, x_col, y_col, bin_width, min_n) {
  df %>%
    filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]])) %>%
    mutate(xbin = floor(.data[[x_col]] / bin_width) * bin_width) %>%
    group_by(Species, xbin) %>%
    summarise(
      x_mean = mean(.data[[x_col]], na.rm = TRUE),
      y_mean = mean(.data[[y_col]], na.rm = TRUE),
      n      = n(),
      .groups = "drop"
    ) %>%
    filter(n >= min_n)
}

# ---------------------------
# Plot functions
# ---------------------------
plot_combined_binned <- function(raw_df, binned_df, y_col, title, ylab, filename, out_dir) {
  p <- ggplot() +
    geom_point(
      data = raw_df,
      aes(x = psi_leaf, y = .data[[y_col]], color = Species),
      alpha = 0.05, size = 0.2
    ) +
    geom_point(
      data = binned_df,
      aes(x = x_mean, y = y_mean, color = Species),
      size = 2
    ) +
    geom_smooth(
      data = binned_df,
      aes(x = x_mean, y = y_mean, color = Species, fill = Species),
      se = TRUE, linewidth = 1, alpha = 0.20
    ) +
    scale_color_manual(values = cb_palette) +
    scale_fill_manual(values = cb_palette) +
    labs(
      title = title,
      x = expression(paste("Leaf water potential ", psi[leaf], " (MPa)")),
      y = ylab,
      color = NULL,
      fill = NULL
    ) +
    theme_clean() +
    theme(legend.position = "top")
  
  ggsave(file.path(out_dir, filename), p, width = 11, height = 6, dpi = 300)
  p
}

plot_facet_binned <- function(raw_df, binned_df, y_col, title, ylab, filename, out_dir) {
  p <- ggplot() +
    geom_point(
      data = raw_df,
      aes(x = psi_leaf, y = .data[[y_col]], color = Species),
      alpha = 0.08, size = 0.5
    ) +
    geom_point(
      data = binned_df,
      aes(x = x_mean, y = y_mean, color = Species),
      size = 2
    ) +
    geom_smooth(
      data = binned_df,
      aes(x = x_mean, y = y_mean, color = Species, fill = Species),
      se = TRUE, linewidth = 1, alpha = 0.20
    ) +
    facet_wrap(~Species, ncol = 2, scales = "free") +
    scale_color_manual(values = cb_palette) +
    scale_fill_manual(values = cb_palette) +
    labs(
      title = title,
      x = expression(paste("Leaf water potential ", psi[leaf], " (MPa)")),
      y = ylab
    ) +
    theme_clean() +
    theme(legend.position = "none")
  
  ggsave(file.path(out_dir, filename), p, width = 11, height = 6, dpi = 300)
  p
}

# ---------------------------
# Generic wrapper
# ---------------------------
make_relationship_plots <- function(df, y_col, label_name, ylab, prefix, min_n_bin, out_dir, file_tag) {
  binned <- bin_by_x(df, "psi_leaf", y_col, bin_width_mpa, min_n_bin)
  
  p_all <- plot_combined_binned(
    raw_df    = df,
    binned_df = binned,
    y_col     = y_col,
    title     = paste0(label_name, " vs leaf water potential (", prefix, "; binned 50 kPa)"),
    ylab      = ylab,
    filename  = paste0(y_col, "_vs_psileaf_", file_tag, "_binned50kPa_ALL.png"),
    out_dir   = out_dir
  )
  
  p_facet <- plot_facet_binned(
    raw_df    = df,
    binned_df = binned,
    y_col     = y_col,
    title     = paste0(label_name, " vs leaf water potential (", prefix, "; species panels)"),
    ylab      = ylab,
    filename  = paste0(y_col, "_vs_psileaf_", file_tag, "_binned50kPa_facet.png"),
    out_dir   = out_dir
  )
  
  list(all = p_all, facet = p_facet, binned = binned)
}

# ============================================================
# MAIN: Build full daily datasets once
# ============================================================
df_gcwater_daily_full <- pmap_dfr(
  list(species_map$species, species_map$colname),
  ~ read_pair_psileaf_y(..1, ..2, y_file = "dgcwater.out", y_name = "gcwater")
) %>%
  filter(is.finite(psi_leaf), is.finite(gcwater)) %>%
  mutate(Species = factor(Species, levels = species_map$species))

df_gc_daily_full <- pmap_dfr(
  list(species_map$species, species_map$colname),
  ~ read_pair_psileaf_y(..1, ..2, y_file = "dgc.out", y_name = "gc")
) %>%
  filter(is.finite(psi_leaf), is.finite(gc)) %>%
  mutate(Species = factor(Species, levels = species_map$species))

# ============================================================
# Function to run one season selection
# ============================================================
run_all_plots <- function(df_gcwater_full, df_gc_full, months_keep, season_label, out_dir, file_tag) {
  
  # ---------------------------
  # Filter daily data
  # ---------------------------
  df_gcwater_daily <- apply_time_filter(df_gcwater_full, months_keep = months_keep)
  df_gc_daily      <- apply_time_filter(df_gc_full,      months_keep = months_keep)
  
  # ---------------------------
  # Monthly mean datasets
  # ---------------------------
  df_gcwater_monthly <- aggregate_monthly_mean(df_gcwater_daily, "gcwater")
  df_gc_monthly      <- aggregate_monthly_mean(df_gc_daily, "gc")
  
  # ---------------------------
  # Yearly mean datasets
  # ---------------------------
  df_gcwater_yearly <- aggregate_yearly_mean(df_gcwater_daily, "gcwater")
  df_gc_yearly      <- aggregate_yearly_mean(df_gc_daily, "gc")
  
  # ---------------------------
  # PLOTS: DAILY
  # ---------------------------
  res_gcwater_daily <- make_relationship_plots(
    df         = df_gcwater_daily,
    y_col      = "gcwater",
    label_name = "Canopy conductance (gcwater)",
    ylab       = "Canopy conductance (gcwater, m/s)",
    prefix     = paste0("Daily ", season_label),
    min_n_bin  = min_n_bin_daily,
    out_dir    = out_dir,
    file_tag   = paste0(file_tag, "_Daily")
  )
  
  res_gc_daily <- make_relationship_plots(
    df         = df_gc_daily,
    y_col      = "gc",
    label_name = "Leaf conductance (gc)",
    ylab       = "Leaf conductance (gc, mm/s)",
    prefix     = paste0("Daily ", season_label),
    min_n_bin  = min_n_bin_daily,
    out_dir    = out_dir,
    file_tag   = paste0(file_tag, "_Daily")
  )
  
  # ---------------------------
  # PLOTS: MONTHLY mean
  # ---------------------------
  res_gcwater_monthly <- make_relationship_plots(
    df         = df_gcwater_monthly,
    y_col      = "gcwater",
    label_name = "Canopy conductance (gcwater)",
    ylab       = "Canopy conductance (gcwater, m/s)",
    prefix     = paste0("Monthly mean ", season_label),
    min_n_bin  = min_n_bin_monthly,
    out_dir    = out_dir,
    file_tag   = paste0(file_tag, "_MonthlyMean")
  )
  
  res_gc_monthly <- make_relationship_plots(
    df         = df_gc_monthly,
    y_col      = "gc",
    label_name = "Leaf conductance (gc)",
    ylab       = "Leaf conductance (gc, mm/s)",
    prefix     = paste0("Monthly mean ", season_label),
    min_n_bin  = min_n_bin_monthly,
    out_dir    = out_dir,
    file_tag   = paste0(file_tag, "_MonthlyMean")
  )
  
  # ---------------------------
  # PLOTS: YEARLY mean
  # ---------------------------
  res_gcwater_yearly <- make_relationship_plots(
    df         = df_gcwater_yearly,
    y_col      = "gcwater",
    label_name = "Canopy conductance (gcwater)",
    ylab       = "Canopy conductance (gcwater, m/s)",
    prefix     = paste0("Yearly mean ", season_label),
    min_n_bin  = min_n_bin_yearly,
    out_dir    = out_dir,
    file_tag   = paste0(file_tag, "_YearlyMean")
  )
  
  res_gc_yearly <- make_relationship_plots(
    df         = df_gc_yearly,
    y_col      = "gc",
    label_name = "Leaf conductance (gc)",
    ylab       = "Leaf conductance (gc, mm/s)",
    prefix     = paste0("Yearly mean ", season_label),
    min_n_bin  = min_n_bin_yearly,
    out_dir    = out_dir,
    file_tag   = paste0(file_tag, "_YearlyMean")
  )
  
  list(
    gcwater_daily   = res_gcwater_daily,
    gc_daily        = res_gc_daily,
    gcwater_monthly = res_gcwater_monthly,
    gc_monthly      = res_gc_monthly,
    gcwater_yearly  = res_gcwater_yearly,
    gc_yearly       = res_gc_yearly
  )
}

# ============================================================
# Run 1: Whole year
# ============================================================
res_whole <- run_all_plots(
  df_gcwater_full = df_gcwater_daily_full,
  df_gc_full      = df_gc_daily_full,
  months_keep     = NULL,
  season_label    = "whole year",
  out_dir         = out_dir_whole,
  file_tag        = "WholeYear"
)

# ============================================================
# Run 2: July-August only
# ============================================================
res_julaug <- run_all_plots(
  df_gcwater_full = df_gcwater_daily_full,
  df_gc_full      = df_gc_daily_full,
  months_keep     = c(7, 8),
  season_label    = "July-August",
  out_dir         = out_dir_julaug,
  file_tag        = "JulAug"
)

# ============================================================
# show in RStudio
# ============================================================
# Whole year
res_whole$gcwater_daily$all
res_whole$gcwater_daily$facet
res_whole$gc_daily$all
res_whole$gc_daily$facet

res_whole$gcwater_monthly$all
res_whole$gcwater_monthly$facet
res_whole$gc_monthly$all
res_whole$gc_monthly$facet

res_whole$gcwater_yearly$all
res_whole$gcwater_yearly$facet
res_whole$gc_yearly$all
res_whole$gc_yearly$facet

# July-August
res_julaug$gcwater_daily$all
res_julaug$gcwater_daily$facet
res_julaug$gc_daily$all
res_julaug$gc_daily$facet

res_julaug$gcwater_monthly$all
res_julaug$gcwater_monthly$facet
res_julaug$gc_monthly$all
res_julaug$gc_monthly$facet

res_julaug$gcwater_yearly$all
res_julaug$gcwater_yearly$facet
res_julaug$gc_yearly$all
res_julaug$gc_yearly$facet

# ============================================================
# check the dataframe 
# ============================================================
res_whole$gcwater_daily$binned
res_whole$gcwater_monthly$binned
res_whole$gcwater_yearly$binned

