# ============================================================
# LPJ-GUESS post-processing
# Relationship plots: gc and gcwater vs psi_leaf
# Timescales: Daily, Monthly Mean, Yearly Mean
# Selections: Whole Period, Drought Years, and Growing Season
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

# Output Directories
out_dir_whole   <- file.path(base_dir, "Figures", "gc_gcwater_vs_psileaf_WholeYear")
out_dir_drought <- file.path(base_dir, "Figures", "gc_gcwater_vs_psileaf_Drought_JulAug")
out_dir_gs      <- file.path(base_dir, "Figures", "gc_gcwater_vs_psileaf_GrowingSeason_14_16")

dir.create(out_dir_whole,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_drought, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_gs,      recursive = TRUE, showWarnings = FALSE)

day_starts_at_zero <- TRUE

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
  if (length(numeric_candidates) == 0) stop("No numeric column found in ", file_for_msg)
  numeric_candidates[1]
}

read_daily_value <- function(species, colname, filename, value_name) {
  f <- file.path(results_dir, species, filename)
  if (!file.exists(f)) return(NULL)
  
  df <- read.table(f, header = TRUE, check.names = FALSE)
  names(df)[3:4] <- c("Year", "Day")
  col_use <- pick_numeric_col(df, colname, f)
  
  df %>%
    mutate(
      Year = as.numeric(Year),
      Day  = as.numeric(Day),
      Date = as_guess_date(Year, Day, day_starts_at_zero)
    ) %>%
    transmute(Date, Species = species, !!value_name := .data[[col_use]])
}

read_pair_psileaf_y <- function(species, colname, y_file, y_name) {
  psi <- read_daily_value(species, colname, "dpsileaf.out", "psi_leaf")
  yy  <- read_daily_value(species, colname, y_file, y_name)
  if (is.null(psi) || is.null(yy)) return(NULL)
  inner_join(psi, yy, by = c("Date", "Species"))
}

# ---------------------------
# Filtering Logic
# ---------------------------
apply_time_filter <- function(df, mode = "all") {
  if (is.null(df) || nrow(df) == 0) return(df)
  
  if (mode == "drought") {
    drought_years <- c(2003, 2011, 2015, 2018, 2022)
    df <- df %>% 
      filter(lubridate::year(Date) %in% drought_years,
             lubridate::month(Date) %in% c(7, 8))
    
  } else if (mode == "growing_season") {
    # April to October for years 2014, 2015, 2016
    df <- df %>% 
      filter(lubridate::year(Date) %in% 2014:2016,
             lubridate::month(Date) >= 4, 
             lubridate::month(Date) <= 10)
  }
  df
}

# ---------------------------
# Aggregation & Binning
# ---------------------------
aggregate_monthly_mean <- function(df, y_col) {
  df %>%
    group_by(Species, Year = year(Date), Month = month(Date)) %>%
    summarise(psi_leaf = mean(psi_leaf, na.rm = TRUE),
              value = mean(.data[[y_col]], na.rm = TRUE), .groups = "drop") %>%
    rename(!!y_col := value) %>%
    mutate(Species = factor(Species, levels = species_map$species))
}

aggregate_yearly_mean <- function(df, y_col) {
  df %>%
    group_by(Species, Year = year(Date)) %>%
    summarise(psi_leaf = mean(psi_leaf, na.rm = TRUE),
              value = mean(.data[[y_col]], na.rm = TRUE), .groups = "drop") %>%
    rename(!!y_col := value) %>%
    mutate(Species = factor(Species, levels = species_map$species))
}

bin_by_x <- function(df, x_col, y_col, bin_width, min_n) {
  df %>%
    filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]])) %>%
    mutate(xbin = floor(.data[[x_col]] / bin_width) * bin_width) %>%
    group_by(Species, xbin) %>%
    summarise(x_mean = mean(.data[[x_col]], na.rm = TRUE),
              y_mean = mean(.data[[y_col]], na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    filter(n >= min_n)
}

# ---------------------------
# Plotting
# ---------------------------
plot_relationship <- function(raw_df, binned_df, y_col, title, ylab, filename, out_dir, facet = FALSE) {
  p <- ggplot() +
    geom_point(data = raw_df, aes(x = psi_leaf, y = .data[[y_col]], color = Species), alpha = 0.05, size = 0.2) +
    geom_point(data = binned_df, aes(x = x_mean, y = y_mean, color = Species), size = 2) +
    geom_smooth(data = binned_df, aes(x = x_mean, y = y_mean, color = Species, fill = Species),
                method = "gam", formula = y ~ s(x, k = 5), se = TRUE, linewidth = 1, alpha = 0.20) +
    scale_color_manual(values = cb_palette, drop = FALSE) +
    scale_fill_manual(values = cb_palette, drop = FALSE) +
    labs(title = title, x = expression(paste(psi[leaf], " (MPa)")), y = ylab) +
    theme_clean()
  
  if (facet) {
    p <- p + facet_wrap(~Species, ncol = 2, scales = "free") + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "top")
  }
  
  ggsave(file.path(out_dir, filename), p, width = 11, height = 7, dpi = 300)
  p
}

make_relationship_plots <- function(df, y_col, label_name, ylab, prefix, min_n_bin, out_dir, file_tag) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  binned <- bin_by_x(df, "psi_leaf", y_col, bin_width_mpa, min_n_bin)
  if (nrow(binned) == 0) return(NULL)
  
  list(
    all = plot_relationship(df, binned, y_col, paste0(label_name, " (", prefix, ")"), ylab, 
                            paste0(y_col, "_vs_psi_", file_tag, "_ALL.png"), out_dir, facet = FALSE),
    facet = plot_relationship(df, binned, y_col, paste0(label_name, " (", prefix, ")"), ylab, 
                              paste0(y_col, "_vs_psi_", file_tag, "_facet.png"), out_dir, facet = TRUE),
    binned = binned
  )
}

# ---------------------------
# Analysis Wrapper
# ---------------------------
run_analysis <- function(df_gcw, df_gc, out_dir, tag, label) {
  # Daily
  res_gcw_d <- make_relationship_plots(df_gcw, "gcwater", "Canopy Conductance", "gcwater (m/s)", paste0("Daily ", label), min_n_bin_daily, out_dir, paste0(tag, "_Daily"))
  res_gc_d  <- make_relationship_plots(df_gc,  "gc",      "Leaf Conductance",   "gc (mm/s)",    paste0("Daily ", label), min_n_bin_daily, out_dir, paste0(tag, "_Daily"))
  
  # Monthly
  df_gcw_m  <- aggregate_monthly_mean(df_gcw, "gcwater")
  df_gc_m   <- aggregate_monthly_mean(df_gc, "gc")
  res_gcw_m <- make_relationship_plots(df_gcw_m, "gcwater", "Canopy Conductance", "gcwater (m/s)", paste0("Monthly ", label), min_n_bin_monthly, out_dir, paste0(tag, "_Monthly"))
  res_gc_m  <- make_relationship_plots(df_gc_m,  "gc",      "Leaf Conductance",   "gc (mm/s)",    paste0("Monthly ", label), min_n_bin_monthly, out_dir, paste0(tag, "_Monthly"))
  
  # Yearly
  df_gcw_y  <- aggregate_yearly_mean(df_gcw, "gcwater")
  df_gc_y   <- aggregate_yearly_mean(df_gc, "gc")
  res_gcw_y <- make_relationship_plots(df_gcw_y, "gcwater", "Canopy Conductance", "gcwater (m/s)", paste0("Yearly ", label), min_n_bin_yearly, out_dir, paste0(tag, "_Yearly"))
  res_gc_y  <- make_relationship_plots(df_gc_y,  "gc",      "Leaf Conductance",   "gc (mm/s)",    paste0("Yearly ", label), min_n_bin_yearly, out_dir, paste0(tag, "_Yearly"))
  
  list(gcw_d = res_gcw_d, gc_d = res_gc_d, gcw_m = res_gcw_m, gc_m = res_gc_m, gcw_y = res_gcw_y, gc_y = res_gc_y)
}


# ============================================================
# MAIN EXECUTION
# ============================================================

# 1. Load Data
species_gcwater <- pmap(list(species_map$species, species_map$colname), ~ read_pair_psileaf_y(..1, ..2, "dgcwater.out", "gcwater")) %>% compact()
species_gc      <- pmap(list(species_map$species, species_map$colname), ~ read_pair_psileaf_y(..1, ..2, "dgc.out", "gc")) %>% compact()

df_gcwater_all <- bind_rows(species_gcwater) %>% filter(is.finite(psi_leaf), is.finite(gcwater)) %>% mutate(Species = factor(Species, levels = species_map$species))
df_gc_all      <- bind_rows(species_gc)      %>% filter(is.finite(psi_leaf), is.finite(gc))      %>% mutate(Species = factor(Species, levels = species_map$species))

# 2. Run for Whole Period
results_whole <- run_analysis(df_gcwater_all, df_gc_all, out_dir_whole, "WholeYear", "Whole Period")

# 3. Run for Drought (Jul/Aug of 03, 11, 15, 18, 22)
df_gcw_drought <- apply_time_filter(df_gcwater_all, mode = "drought")
df_gc_drought  <- apply_time_filter(df_gc_all,      mode = "drought")
results_drought <- run_analysis(df_gcw_drought, df_gc_drought, out_dir_drought, "Drought", "Jul-Aug 03,11,15,18,22")

# 4. Run for Growing Season (Apr-Oct 2014, 2015, 2016)
df_gcw_gs <- apply_time_filter(df_gcwater_all, mode = "growing_season")
df_gc_gs  <- apply_time_filter(df_gc_all,      mode = "growing_season")
results_gs <- run_analysis(df_gcw_gs, df_gc_gs, out_dir_gs, "GrowingSeason", "Apr-Oct 2014-2016")

# ============================================================
# END OF SCRIPT
# ============================================================