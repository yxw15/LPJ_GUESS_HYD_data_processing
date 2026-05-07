# ============================================================
# LPJ-GUESS post-processing
# kappa_s_today.out vs dpsixylem.out
#
# Data location:
#   results/species_folder/*.out
#
# Outputs:
#   1) RAW scatter:
#      - combined ALL species
#      - facet by species
#
#   2) BINNED relationship:
#      - x bin width = 0.01
#      - mean x and mean y within each bin
#      - geom_smooth() based on binned points
#      - combined ALL species
#      - facet by species
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(mgcv)
})

# ---------------------------
# Working directory / paths
# ---------------------------
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

base_dir    <- getwd()
results_dir <- file.path(base_dir, "results")
out_dir     <- file.path(base_dir, "Figures", "kappa_s_today_vs_dpsixylem")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# If LPJ-GUESS Day starts at 0 (often true), keep TRUE
day_starts_at_zero <- TRUE

# Optional: restrict plotted years (set to NULL to disable)
year_min <- 1900
year_max <- NULL

# Binning settings
bin_width_psixylem <- 0.01
min_n_bin <- 1

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
  Oak_pub = "darkorange",
  Oak_rob = "#F0E442",
  Beech   = "dodgerblue",
  Spruce  = "green4",
  Pine    = "purple1"
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

apply_year_filter <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (!is.null(year_min)) df <- df %>% filter(lubridate::year(Date) >= year_min)
  if (!is.null(year_max)) df <- df %>% filter(lubridate::year(Date) <= year_max)
  df
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

# ---- Bin by x ----
bin_by_x <- function(df, x_col, y_col, bin_width, min_n = 1) {
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
# Build dataset safely
# ---------------------------
species_results <- pmap(
  list(species_map$species, species_map$colname),
  function(species, colname) {
    message("Processing species: ", species)
    
    kappa <- read_daily_value(species, colname, "kappa_s_today.out", "kappa_s_today")
    dpsi  <- read_daily_value(species, colname, "dpsixylem.out", "dpsixylem")
    
    if (is.null(kappa) || is.null(dpsi)) {
      message("  -> skipped (missing or unreadable input)")
      return(NULL)
    }
    
    joined <- inner_join(kappa, dpsi, by = c("Date", "Species"))
    
    if (nrow(joined) == 0) {
      message("  -> skipped (join returned 0 rows)")
      return(NULL)
    }
    
    message("  -> joined rows: ", nrow(joined))
    joined
  }
)

species_results <- compact(species_results)

if (length(species_results) == 0) {
  stop(
    "No valid merged data could be created.\n",
    "Please check:\n",
    "1) files exist under results/<species>/kappa_s_today.out and results/<species>/dpsixylem.out\n",
    "2) files contain at least 4 columns\n",
    "3) columns 3 and 4 correspond to Year and Day\n",
    "4) the target species columns exist or at least one numeric column is available\n",
    "5) Date values match between the two files"
  )
}

df_kappa_psixylem <- bind_rows(species_results) %>%
  filter(is.finite(kappa_s_today), is.finite(dpsixylem)) %>%
  mutate(Species = factor(Species, levels = species_map$species)) %>%
  apply_year_filter()

if (nrow(df_kappa_psixylem) == 0) {
  stop(
    "Merged dataset exists but has 0 rows after filtering.\n",
    "Check year_min/year_max and whether kappa_s_today / dpsixylem contain finite values."
  )
}

# Save merged raw daily data
write.csv(
  df_kappa_psixylem,
  file.path(out_dir, "kappa_s_today_vs_dpsixylem_merged_daily.csv"),
  row.names = FALSE
)

# ---------------------------
# RAW scatter plots
# ---------------------------
p_kappa_psixylem_raw_all <- ggplot(
  df_kappa_psixylem,
  aes(x = dpsixylem, y = kappa_s_today, color = Species)
) +
  geom_point(alpha = 0.25, size = 1) +
  scale_color_manual(values = cb_palette, drop = FALSE) +
  labs(
    title = "Cavitation fraction vs dpsixylem (RAW)",
    x = "dpsixylem",
    y = "Cavitation fraction (kappa_s_today)",
    color = NULL
  ) +
  theme_clean() +
  theme(legend.position = "top")

ggsave(
  file.path(out_dir, "kappa_s_today_vs_dpsixylem_scatter_RAW_ALL.png"),
  p_kappa_psixylem_raw_all,
  width = 10,
  height = 6,
  dpi = 300
)

p_kappa_psixylem_raw_facet <- ggplot(
  df_kappa_psixylem,
  aes(x = dpsixylem, y = kappa_s_today, color = Species)
) +
  geom_point(alpha = 0.25, size = 1) +
  facet_wrap(~Species, ncol = 2, scales = "free") +
  scale_color_manual(values = cb_palette, drop = FALSE) +
  labs(
    title = "Cavitation fraction vs dpsixylem (RAW) — species panels",
    x = "dpsixylem",
    y = "Cavitation fraction (kappa_s_today)"
  ) +
  theme_clean() +
  theme(legend.position = "none")

ggsave(
  file.path(out_dir, "kappa_s_today_vs_dpsixylem_scatter_RAW_facet.png"),
  p_kappa_psixylem_raw_facet,
  width = 11,
  height = 6,
  dpi = 300
)

# ---------------------------
# Bin data
# ---------------------------
b_kappa_psixylem <- bin_by_x(
  df = df_kappa_psixylem,
  x_col = "dpsixylem",
  y_col = "kappa_s_today",
  bin_width = bin_width_psixylem,
  min_n = min_n_bin
)

if (nrow(b_kappa_psixylem) == 0) {
  stop("Binned dataset has 0 rows. Try reducing min_n_bin or check the data range.")
}

# Save binned data
write.csv(
  b_kappa_psixylem,
  file.path(out_dir, "kappa_s_today_vs_dpsixylem_binned_0p01.csv"),
  row.names = FALSE
)

# ---------------------------
# BINNED + smooth plots
# raw points in grey
# binned means in species colors
# smooth fitted to binned means using GAM
# ---------------------------
p_kappa_psixylem_binned_all <- ggplot() +
  geom_point(
    data = df_kappa_psixylem,
    aes(x = dpsixylem, y = kappa_s_today),
    color = "grey60",
    alpha = 0.08,
    size = 0.5
  ) +
  geom_point(
    data = b_kappa_psixylem,
    aes(x = x_mean, y = y_mean, color = Species),
    size = 2
  ) +
  geom_smooth(
    data = b_kappa_psixylem,
    aes(x = x_mean, y = y_mean, color = Species, fill = Species),
    method = "gam",
    formula = y ~ s(x, k = 6),
    se = TRUE,
    linewidth = 1,
    alpha = 0.20
  ) +
  scale_color_manual(values = cb_palette, drop = FALSE) +
  scale_fill_manual(values = cb_palette, drop = FALSE) +
  labs(
    title = "Cavitation fraction vs dpsixylem (binned x = 0.01) — ALL species",
    x = "dpsixylem",
    y = "Cavitation fraction (kappa_s_today)",
    color = NULL,
    fill = NULL
  ) +
  theme_clean() +
  theme(legend.position = "top")

ggsave(
  file.path(out_dir, "kappa_s_today_vs_dpsixylem_binned0.01_smooth_ALL.png"),
  p_kappa_psixylem_binned_all,
  width = 11,
  height = 6,
  dpi = 300
)

p_kappa_psixylem_binned_facet <- ggplot() +
  geom_point(
    data = df_kappa_psixylem,
    aes(x = dpsixylem, y = kappa_s_today),
    color = "grey60",
    alpha = 0.08,
    size = 0.5
  ) +
  geom_point(
    data = b_kappa_psixylem,
    aes(x = x_mean, y = y_mean, color = Species),
    size = 2
  ) +
  geom_smooth(
    data = b_kappa_psixylem,
    aes(x = x_mean, y = y_mean, color = Species, fill = Species),
    method = "gam",
    formula = y ~ s(x, k = 6),
    se = TRUE,
    linewidth = 1,
    alpha = 0.20
  ) +
  facet_wrap(~Species, ncol = 2, scales = "free") +
  scale_color_manual(values = cb_palette, drop = FALSE) +
  scale_fill_manual(values = cb_palette, drop = FALSE) +
  labs(
    title = "Cavitation fraction vs dpsixylem (binned x = 0.01) — species panels",
    x = "dpsixylem",
    y = "Cavitation fraction (kappa_s_today)"
  ) +
  theme_clean() +
  theme(legend.position = "none")

ggsave(
  file.path(out_dir, "kappa_s_today_vs_dpsixylem_binned0.01_smooth_facet.png"),
  p_kappa_psixylem_binned_facet,
  width = 11,
  height = 6,
  dpi = 300
)

# ---------------------------
# Optional: print summary
# ---------------------------
cat("\nSaved files to:\n", out_dir, "\n")
cat("\nNumber of raw rows:", nrow(df_kappa_psixylem), "\n")
cat("Number of binned rows:", nrow(b_kappa_psixylem), "\n")

# ---------------------------
# Show in RStudio
# ---------------------------
print(p_kappa_psixylem_raw_all)
print(p_kappa_psixylem_raw_facet)
print(p_kappa_psixylem_binned_all)
print(p_kappa_psixylem_binned_facet)