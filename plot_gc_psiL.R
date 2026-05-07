# ============================================================
# LPJ-GUESS post-processing
# Relationship plots ONLY:
#   1) gc      (leaf conductance)   vs psi_leaf
#   2) gcwater (canopy conductance) vs psi_leaf
#
# Outputs for each time selection:
#   - combined ALL species
#   - facet by species
#
# Timescales:
#   - daily raw
#   - monthly mean
#   - yearly mean
#
# Time selections:
#   - whole year
#   - July, and August of 2003, 2011, 2015, 2018, 2022
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
base_dir    <- getwd()
results_dir <- file.path(base_dir, "results")

out_dir_whole   <- file.path(base_dir, "Figures", "gc_gcwater_vs_psileaf_WholeYear")
out_dir_2015gs  <- file.path(base_dir, "Figures", "gc_gcwater_vs_psileaf_2015_Jul01_to_Sep15")

dir.create(out_dir_whole,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_2015gs, recursive = TRUE, showWarnings = FALSE)

day_starts_at_zero <- TRUE

# Optional: restrict years globally
year_min <- NULL
year_max <- NULL

# Binning settings
bin_width_mpa      <- 0.05
min_n_bin_daily    <- 1
min_n_bin_monthly  <- 1
min_n_bin_yearly   <- 1

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
  offset <- if (day_starts_at_zero) Day else Day - 1
  as.Date(offset, origin = paste0(Year, "-01-01"))
}

theme_clean <- function() {
  theme_minimal() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      strip.text       = element_text(face = "bold", size = 12)
    )
}

pick_numeric_col <- function(df, wanted = NULL, file_for_msg = "") {
  if (!is.null(wanted) && wanted %in% names(df)) return(wanted)
  
  candidates <- setdiff(names(df), c("Year", "Day", "Lon", "Lat"))
  numeric_candidates <- candidates[sapply(df[candidates], is.numeric)]
  
  if (length(numeric_candidates) == 0) {
    stop(
      "No numeric column found in ", file_for_msg,
      ". Available: ", paste(names(df), collapse = ", ")
    )
  }
  
  message(
    "Note: wanted column '", wanted, "' not found in ", basename(file_for_msg),
    ". Using first numeric column: '", numeric_candidates[1], "'"
  )
  
  numeric_candidates[1]
}

read_daily_value <- function(species, colname, filename, value_name) {
  f <- file.path(results_dir, species, filename)
  
  if (!file.exists(f)) {
    warning("Missing file: ", f)
    return(NULL)
  }
  
  df <- tryCatch(
    read.table(f, header = TRUE, check.names = FALSE),
    error = function(e) {
      warning("Could not read file: ", f, " | ", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (is.null(df)) return(NULL)
  
  if (ncol(df) < 4) {
    warning("File has fewer than 4 columns: ", f)
    return(NULL)
  }
  
  # LPJ-GUESS files often store year/day in columns 3 and 4
  names(df)[3:4] <- c("Year", "Day")
  
  col_use <- tryCatch(
    pick_numeric_col(df, colname, f),
    error = function(e) {
      warning(conditionMessage(e))
      return(NULL)
    }
  )
  
  if (is.null(col_use)) return(NULL)
  
  df %>%
    mutate(
      Year = as.numeric(Year),
      Day  = as.numeric(Day),
      Date = as_guess_date(Year, Day, day_starts_at_zero)
    ) %>%
    transmute(
      Date,
      Species = species,
      !!value_name := .data[[col_use]]
    )
}

read_pair_psileaf_y <- function(species, colname, y_file, y_name) {
  psi <- read_daily_value(species, colname, "dpsileaf.out", "psi_leaf")
  yy  <- read_daily_value(species, colname, y_file, y_name)
  
  if (is.null(psi) || is.null(yy)) return(NULL)
  
  joined <- inner_join(psi, yy, by = c("Date", "Species"))
  
  if (nrow(joined) == 0) {
    message("Join returned 0 rows for species: ", species, " and file: ", y_file)
    return(NULL)
  }
  
  joined
}

apply_time_filter <- function(df, start_date = NULL, end_date = NULL) {
  if (is.null(df) || nrow(df) == 0) return(df)
  
  if (!is.null(year_min)) {
    df <- df %>% filter(lubridate::year(Date) >= year_min)
  }
  
  if (!is.null(year_max)) {
    df <- df %>% filter(lubridate::year(Date) <= year_max)
  }
  
  if (!is.null(start_date)) {
    df <- df %>% filter(Date >= as.Date(start_date))
  }
  
  if (!is.null(end_date)) {
    df <- df %>% filter(Date <= as.Date(end_date))
  }
  
  df
}

# ---------------------------
# Aggregation
# ---------------------------
aggregate_monthly_mean <- function(df, y_col) {
  if (is.null(df) || nrow(df) == 0) return(df[0, ])
  
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
  if (is.null(df) || nrow(df) == 0) return(df[0, ])
  
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
      method = "gam",
      formula = y ~ s(x, k = 5),
      se = TRUE, linewidth = 1, alpha = 0.20
    ) +
    scale_color_manual(values = cb_palette, drop = FALSE) +
    scale_fill_manual(values = cb_palette, drop = FALSE) +
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
      method = "gam",
      formula = y ~ s(x, k = 5),
      se = TRUE, linewidth = 1, alpha = 0.20
    ) +
    facet_wrap(~Species, ncol = 2, scales = "free") +
    scale_color_manual(values = cb_palette, drop = FALSE) +
    scale_fill_manual(values = cb_palette, drop = FALSE) +
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
  if (is.null(df) || nrow(df) == 0) {
    warning("No data available for: ", file_tag)
    return(list(all = NULL, facet = NULL, binned = NULL))
  }
  
  binned <- bin_by_x(df, "psi_leaf", y_col, bin_width_mpa, min_n_bin)
  
  if (nrow(binned) == 0) {
    warning("No binned data available for: ", file_tag)
    return(list(all = NULL, facet = NULL, binned = binned))
  }
  
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
species_gcwater <- pmap(
  list(species_map$species, species_map$colname),
  ~ read_pair_psileaf_y(..1, ..2, y_file = "dgcwater.out", y_name = "gcwater")
)

species_gc <- pmap(
  list(species_map$species, species_map$colname),
  ~ read_pair_psileaf_y(..1, ..2, y_file = "dgc.out", y_name = "gc")
)

species_gcwater <- compact(species_gcwater)
species_gc      <- compact(species_gc)

if (length(species_gcwater) == 0) {
  stop("No valid gcwater dataset could be created.")
}

if (length(species_gc) == 0) {
  stop("No valid gc dataset could be created.")
}

df_gcwater_daily_full <- bind_rows(species_gcwater) %>%
  filter(is.finite(psi_leaf), is.finite(gcwater)) %>%
  mutate(Species = factor(Species, levels = species_map$species))

df_gc_daily_full <- bind_rows(species_gc) %>%
  filter(is.finite(psi_leaf), is.finite(gc)) %>%
  mutate(Species = factor(Species, levels = species_map$species))

# ============================================================
# Function to run one time selection
# ============================================================
run_all_plots <- function(df_gcwater_full, df_gc_full, start_date, end_date, season_label, out_dir, file_tag) {
  
  # ---------------------------
  # Filter daily data
  # ---------------------------
  df_gcwater_daily <- apply_time_filter(df_gcwater_full, start_date = start_date, end_date = end_date)
  df_gc_daily      <- apply_time_filter(df_gc_full,      start_date = start_date, end_date = end_date)
  
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
# Run 1: Whole period
# ============================================================
res_whole <- run_all_plots(
  df_gcwater_full = df_gcwater_daily_full,
  df_gc_full      = df_gc_daily_full,
  start_date      = NULL,
  end_date        = NULL,
  season_label    = "whole period",
  out_dir         = out_dir_whole,
  file_tag        = "WholeYear"
)

# ============================================================
# Run 2: 2015-07-01 to 2015-09-15
# ============================================================
res_2015_gs <- run_all_plots(
  df_gcwater_full = df_gcwater_daily_full,
  df_gc_full      = df_gc_daily_full,
  start_date      = "2015-07-01",
  end_date        = "2015-09-15",
  season_label    = "2015 Jul 1 to Sep 15",
  out_dir         = out_dir_2015gs,
  file_tag        = "2015Jul01toSep15"
)

# ============================================================
# show in RStudio
# ============================================================
# Whole period
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

# 2015-07-01 to 2015-09-15
res_2015_gs$gcwater_daily$all
res_2015_gs$gcwater_daily$facet
res_2015_gs$gc_daily$all
res_2015_gs$gc_daily$facet

res_2015_gs$gcwater_monthly$all
res_2015_gs$gcwater_monthly$facet
res_2015_gs$gc_monthly$all
res_2015_gs$gc_monthly$facet

res_2015_gs$gcwater_yearly$all
res_2015_gs$gcwater_yearly$facet
res_2015_gs$gc_yearly$all
res_2015_gs$gc_yearly$facet

# ============================================================
# check the dataframe
# ============================================================
res_whole$gcwater_daily$binned
res_whole$gcwater_monthly$binned
res_whole$gcwater_yearly$binned

res_2015_gs$gcwater_daily$binned
res_2015_gs$gcwater_monthly$binned
res_2015_gs$gcwater_yearly$binned