# ============================================================
# LPJ-GUESS post-processing
# kappa_s_today.out vs dpsixylem.out
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

base_dir <- getwd()
out_dir  <- file.path(base_dir, "Figures", "kappa_s_today_vs_psixylem")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# If LPJ-GUESS Day starts at 0 (often true), keep TRUE
day_starts_at_zero <- TRUE

# Optional: restrict plotted years (set to NULL to disable)
year_min <- 1900
year_max <- NULL

# Binning settings
bin_width_psixylem <- 0.01
min_n_bin <- 0

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
  as.Date(
    if (day_starts_at_zero) Day else Day - 1,
    origin = paste0(Year, "-01-01")
  )
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
    stop(
      "No numeric column found in ", file_for_msg,
      ". Available: ", paste(names(df), collapse = ", ")
    )
  }
  
  message(
    "Note: using numeric column '", numeric_candidates[1],
    "' from ", basename(file_for_msg)
  )
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

apply_year_filter <- function(df) {
  if (is.null(df)) return(NULL)
  if (!is.null(year_min)) df <- df %>% filter(lubridate::year(Date) >= year_min)
  if (!is.null(year_max)) df <- df %>% filter(lubridate::year(Date) <= year_max)
  df
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
# Build dataset
# ---------------------------
df_kappa_psixylem <- pmap_dfr(
  list(species_map$species, species_map$colname),
  ~{
    kappa     <- read_daily_value(..1, ..2, "kappa_s_today.out", "kappa_s_today")
    psixylem  <- read_daily_value(..1, ..2, "dpsixylem.out", "psixylem")
    
    if (is.null(kappa) || is.null(psixylem)) return(NULL)
    
    inner_join(kappa, psixylem, by = c("Date", "Species"))
  }
) %>%
  filter(is.finite(kappa_s_today), is.finite(psixylem)) %>%
  mutate(Species = factor(Species, levels = species_map$species)) %>%
  apply_year_filter()

# Save merged raw daily data
write.csv(
  df_kappa_psixylem,
  file.path(out_dir, "kappa_s_today_vs_psixylem_merged_daily.csv"),
  row.names = FALSE
)

# ---------------------------
# RAW scatter plots
# ---------------------------
p_kappa_psixylem_raw_all <- ggplot(
  df_kappa_psixylem,
  aes(x = psixylem, y = kappa_s_today, color = Species)
) +
  geom_point(alpha = 0.25, size = 1) +
  scale_color_manual(values = cb_palette) +
  labs(
    title = "Cavitation fraction vs xylem water potential (RAW)",
    x = "Xylem water potential (psixylem)",
    y = "Cavitation fraction (kappa_s_today)",
    color = NULL
  ) +
  theme_clean() +
  theme(legend.position = "top")

ggsave(
  file.path(out_dir, "kappa_s_today_vs_psixylem_scatter_RAW_ALL.png"),
  p_kappa_psixylem_raw_all,
  width = 10,
  height = 6,
  dpi = 300
)

p_kappa_psixylem_raw_facet <- ggplot(
  df_kappa_psixylem,
  aes(x = psixylem, y = kappa_s_today, color = Species)
) +
  geom_point(alpha = 0.25, size = 1) +
  facet_wrap(~Species, ncol = 2, scales = "free") +
  scale_color_manual(values = cb_palette) +
  labs(
    title = "Cavitation fraction vs xylem water potential (RAW) — species panels",
    x = "Xylem water potential (psixylem)",
    y = "Cavitation fraction (kappa_s_today)"
  ) +
  theme_clean() +
  theme(legend.position = "none")

ggsave(
  file.path(out_dir, "kappa_s_today_vs_psixylem_scatter_RAW_facet.png"),
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
  x_col = "psixylem",
  y_col = "kappa_s_today",
  bin_width = bin_width_psixylem,
  min_n = min_n_bin
)

# Save binned data
write.csv(
  b_kappa_psixylem,
  file.path(out_dir, "kappa_s_today_vs_psixylem_binned_0p01.csv"),
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
    aes(x = psixylem, y = kappa_s_today),
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
    se = TRUE,
    linewidth = 1,
    alpha = 0.20
  ) +
  scale_color_manual(values = cb_palette) +
  scale_fill_manual(values = cb_palette) +
  labs(
    title = "Cavitation fraction vs xylem water potential (binned x = 0.01) — ALL species",
    x = "Xylem water potential (psixylem)",
    y = "Cavitation fraction (kappa_s_today)",
    color = NULL,
    fill = NULL
  ) +
  theme_clean() +
  theme(legend.position = "top")

ggsave(
  file.path(out_dir, "kappa_s_today_vs_psixylem_binned0.01_smooth_ALL.png"),
  p_kappa_psixylem_binned_all,
  width = 11,
  height = 6,
  dpi = 300
)

p_kappa_psixylem_binned_facet <- ggplot() +
  geom_point(
    data = df_kappa_psixylem,
    aes(x = psixylem, y = kappa_s_today),
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
    se = TRUE,
    linewidth = 1,
    alpha = 0.20
  ) +
  facet_wrap(~Species, ncol = 2, scales = "free") +
  scale_color_manual(values = cb_palette) +
  scale_fill_manual(values = cb_palette) +
  labs(
    title = "Cavitation fraction vs xylem water potential (binned x = 0.01) — species panels",
    x = "Xylem water potential (psixylem)",
    y = "Cavitation fraction (kappa_s_today)"
  ) +
  theme_clean() +
  theme(legend.position = "none")

ggsave(
  file.path(out_dir, "kappa_s_today_vs_psixylem_binned0.01_smooth_facet.png"),
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
p_kappa_psixylem_raw_all
p_kappa_psixylem_raw_facet
p_kappa_psixylem_binned_all
p_kappa_psixylem_binned_facet