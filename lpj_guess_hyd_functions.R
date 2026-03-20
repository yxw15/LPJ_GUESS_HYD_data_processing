# ============================================================
# LPJ-GUESS FUNCTIONS
# Logic:
# - native scale plots = raw values
# - coarser scale plots = aggregated mean/min/max
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(purrr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
})

# ============================================================
# Helpers
# ============================================================

# If LPJ-GUESS years are simulation years like 901, 902, ...
# and you want them interpreted as 1901, 1902, ...
# then set year_offset = 1000 in your function calls.
# If your files already contain real calendar years, use year_offset = 0.
to_calendar_year <- function(year_vec, year_offset = 0) {
  year_vec + year_offset
}

as_guess_date <- function(Year, Day, day_starts_at_zero = TRUE, year_offset = 0) {
  year_use <- to_calendar_year(Year, year_offset)
  as.Date(
    if (day_starts_at_zero) Day else Day - 1,
    origin = paste0(year_use, "-01-01")
  )
}

as_date_min <- function(y) as.Date(paste0(y, "-01-01"))
as_date_max <- function(y) as.Date(paste0(y, "-12-31"))

theme_clean <- function() {
  theme_minimal() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      strip.text = element_text(face = "bold", size = 12)
    )
}

summarise_fun <- function(x, stat = c("mean", "min", "max")) {
  stat <- match.arg(stat)
  x <- x[is.finite(x)]
  
  if (length(x) == 0) return(NA_real_)
  
  switch(
    stat,
    mean = mean(x, na.rm = TRUE),
    min  = min(x, na.rm = TRUE),
    max  = max(x, na.rm = TRUE)
  )
}

scale_label <- function(scale) {
  switch(
    scale,
    daily   = "Daily",
    monthly = "Monthly",
    yearly  = "Yearly",
    stop("Unknown scale: ", scale)
  )
}

species_tag <- function(species) {
  if (is.null(species)) "ALL" else paste(species, collapse = "_")
}

filter_plot_window <- function(df, plot_year_min = NULL, plot_year_max = NULL) {
  if (!is.null(plot_year_min)) {
    df <- df %>% filter(Date >= as.Date(paste0(plot_year_min, "-01-01")))
  }
  
  if (!is.null(plot_year_max)) {
    df <- df %>% filter(Date <= as.Date(paste0(plot_year_max, "-12-31")))
  }
  
  df
}

# ============================================================
# Column selectors
# ============================================================

guess_data_column <- function(df, requested_col, file_path = NULL, species = NULL) {
  meta_cols <- c(
    "Lon", "Lat", "Year", "Day",
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
  )
  
  if (!is.null(requested_col) && requested_col %in% names(df)) {
    return(requested_col)
  }
  
  candidates <- setdiff(names(df), meta_cols)
  numeric_candidates <- candidates[sapply(df[candidates], is.numeric)]
  
  if (length(numeric_candidates) == 0) {
    stop(
      "Requested column '", requested_col,
      "' not found and no numeric data column available",
      if (!is.null(file_path)) paste0(" in ", file_path) else "",
      ". Available columns: ", paste(names(df), collapse = ", ")
    )
  }
  
  if (!is.null(file_path)) {
    message(
      "Note: column '", requested_col, "' not found",
      if (!is.null(species)) paste0(" for species '", species, "'") else "",
      ". Using '", numeric_candidates[1], "' from ", basename(file_path), "."
    )
  }
  
  numeric_candidates[1]
}

guess_yearly_column <- function(df, requested_col, file_path = NULL, species = NULL) {
  meta_cols <- c("Lon", "Lat", "Year")
  
  if (!is.null(requested_col) && requested_col %in% names(df)) {
    return(requested_col)
  }
  
  if ("Total" %in% names(df) && is.numeric(df$Total)) {
    if (!is.null(file_path)) {
      message(
        "Note: column '", requested_col, "' not found",
        if (!is.null(species)) paste0(" for species '", species, "'") else "",
        ". Using 'Total' from ", basename(file_path), "."
      )
    }
    return("Total")
  }
  
  candidates <- setdiff(names(df), meta_cols)
  numeric_candidates <- candidates[sapply(df[candidates], is.numeric)]
  
  if (length(numeric_candidates) == 0) {
    stop(
      "Requested column '", requested_col,
      "' not found and no numeric data column available",
      if (!is.null(file_path)) paste0(" in ", file_path) else "",
      ". Available columns: ", paste(names(df), collapse = ", ")
    )
  }
  
  if (!is.null(file_path)) {
    message(
      "Note: column '", requested_col, "' not found",
      if (!is.null(species)) paste0(" for species '", species, "'") else "",
      ". Using '", numeric_candidates[1], "' from ", basename(file_path), "."
    )
  }
  
  numeric_candidates[1]
}

# ============================================================
# Generic aggregation
# Input df must contain: Date, Species, Variable, Value
# ============================================================

aggregate_timeseries <- function(df,
                                 species_levels,
                                 scale = c("monthly", "yearly"),
                                 stat = c("mean", "min", "max")) {
  scale <- match.arg(scale)
  stat  <- match.arg(stat)
  
  if (!all(c("Date", "Species", "Variable", "Value") %in% names(df))) {
    stop("Input df must contain: Date, Species, Variable, Value")
  }
  
  if (scale == "monthly") {
    out <- df %>%
      mutate(Date = floor_date(Date, "month")) %>%
      group_by(Variable, Species, Date) %>%
      summarise(Value = summarise_fun(Value, stat), .groups = "drop")
  }
  
  if (scale == "yearly") {
    out <- df %>%
      mutate(Year = year(Date)) %>%
      group_by(Variable, Species, Year) %>%
      summarise(Value = summarise_fun(Value, stat), .groups = "drop") %>%
      mutate(Date = as.Date(paste0(Year, "-01-01"))) %>%
      select(Variable, Species, Date, Value)
  }
  
  out %>%
    mutate(
      Scale   = factor(scale_label(scale), levels = c("Yearly", "Monthly", "Daily")),
      Stat    = stat,
      Species = factor(Species, levels = species_levels)
    )
}

# ============================================================
# Readers
# ============================================================

read_daily_one <- function(base_dir, species, colname, filename, var_label,
                           year_min, year_max,
                           day_starts_at_zero = TRUE,
                           year_offset = 0) {
  f <- file.path(base_dir, species, filename)
  
  if (!file.exists(f)) {
    warning("Missing file: ", f)
    return(NULL)
  }
  
  df <- read.table(f, header = TRUE)
  
  if (!all(c("Year", "Day") %in% names(df))) {
    stop("File missing Year/Day columns: ", f)
  }
  
  col_use <- guess_data_column(df, colname, f, species)
  
  df %>%
    mutate(
      Date = as_guess_date(
        Year = Year,
        Day = Day,
        day_starts_at_zero = day_starts_at_zero,
        year_offset = year_offset
      )
    ) %>%
    filter(Date >= as_date_min(year_min), Date <= as_date_max(year_max)) %>%
    transmute(
      Date,
      Species  = species,
      Variable = var_label,
      Value    = .data[[col_use]]
    ) %>%
    filter(is.finite(Value))
}

read_monthlywide_one <- function(base_dir, species, filename, var_label,
                                 year_min, year_max,
                                 year_offset = 0) {
  f <- file.path(base_dir, species, filename)
  
  if (!file.exists(f)) {
    warning("Missing file: ", f)
    return(NULL)
  }
  
  df <- read.table(f, header = TRUE)
  
  if (!("Year" %in% names(df))) {
    stop("Monthly file missing Year: ", f)
  }
  
  month_levels <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  
  if (!all(month_levels %in% names(df))) {
    stop(
      "Monthly file does not look like wide Jan..Dec format: ", f,
      "\nAvailable columns: ", paste(names(df), collapse = ", ")
    )
  }
  
  df %>%
    pivot_longer(
      cols = all_of(month_levels),
      names_to = "Month",
      values_to = "Value"
    ) %>%
    mutate(
      YearCal = to_calendar_year(Year, year_offset),
      Month = factor(Month, levels = month_levels),
      MonthNum = match(as.character(Month), month_levels),
      Date = as.Date(sprintf("%d-%02d-01", YearCal, MonthNum)),
      Species = species,
      Variable = var_label
    ) %>%
    filter(Date >= as_date_min(year_min), Date <= as_date_max(year_max)) %>%
    select(Date, Species, Variable, Value) %>%
    filter(is.finite(Value))
}

read_yearly_one <- function(base_dir, species, colname, filename, var_label,
                            year_min, year_max,
                            year_offset = 0) {
  f <- file.path(base_dir, species, filename)
  
  if (!file.exists(f)) {
    warning("Missing file: ", f)
    return(NULL)
  }
  
  df <- read.table(f, header = TRUE)
  
  if (!("Year" %in% names(df))) {
    stop("Yearly file missing Year: ", f)
  }
  
  col_use <- guess_yearly_column(df, colname, f, species)
  
  df %>%
    mutate(YearCal = to_calendar_year(Year, year_offset)) %>%
    filter(YearCal >= year_min, YearCal <= year_max) %>%
    transmute(
      Date     = as.Date(paste0(YearCal, "-01-01")),
      Species  = species,
      Variable = var_label,
      Value    = .data[[col_use]]
    ) %>%
    filter(is.finite(Value))
}

# ============================================================
# Builders
# ============================================================

build_daily_all <- function(base_dir, daily_vars, species_map,
                            year_min, year_max,
                            day_starts_at_zero = TRUE,
                            year_offset = 0) {
  pmap_dfr(
    list(daily_vars$file, daily_vars$var_label),
    function(f, vlab) {
      pmap_dfr(
        list(species_map$species, species_map$colname),
        ~ read_daily_one(
          base_dir = base_dir,
          species = ..1,
          colname = ..2,
          filename = f,
          var_label = vlab,
          year_min = year_min,
          year_max = year_max,
          day_starts_at_zero = day_starts_at_zero,
          year_offset = year_offset
        )
      )
    }
  )
}

build_monthly_all <- function(base_dir, filename, var_label, species_map,
                              year_min, year_max,
                              year_offset = 0) {
  map_dfr(
    species_map$species,
    ~ read_monthlywide_one(
      base_dir = base_dir,
      species = .x,
      filename = filename,
      var_label = var_label,
      year_min = year_min,
      year_max = year_max,
      year_offset = year_offset
    )
  ) %>%
    filter(is.finite(Value))
}

build_yearly_all <- function(base_dir, filename, var_label, species_map,
                             year_min, year_max,
                             year_offset = 0) {
  pmap_dfr(
    list(species_map$species, species_map$colname),
    ~ read_yearly_one(
      base_dir = base_dir,
      species = ..1,
      colname = ..2,
      filename = filename,
      var_label = var_label,
      year_min = year_min,
      year_max = year_max,
      year_offset = year_offset
    )
  ) %>%
    filter(is.finite(Value))
}

make_daily_aggregates <- function(daily_df, species_map, stats = c("mean", "min", "max")) {
  out <- list()
  
  for (st in stats) {
    out[[paste0("monthly_", st)]] <- aggregate_timeseries(
      daily_df,
      species_levels = species_map$species,
      scale = "monthly",
      stat = st
    )
    
    out[[paste0("yearly_", st)]] <- aggregate_timeseries(
      daily_df,
      species_levels = species_map$species,
      scale = "yearly",
      stat = st
    )
  }
  
  bind_rows(out)
}

make_monthly_aggregates <- function(monthly_df, species_map, stats = c("mean", "min", "max")) {
  out <- list()
  
  for (st in stats) {
    out[[paste0("yearly_", st)]] <- aggregate_timeseries(
      monthly_df,
      species_levels = species_map$species,
      scale = "yearly",
      stat = st
    )
  }
  
  bind_rows(out)
}

# ============================================================
# Plotters
# ============================================================

plot_raw_ts <- function(df,
                        cb_palette,
                        out_dir,
                        var_label,
                        ylab,
                        filename_stub,
                        scale_name = c("daily", "monthly", "yearly"),
                        species = NULL,
                        plot_year_min = NULL,
                        plot_year_max = NULL,
                        width = 11,
                        height = 6) {
  scale_name <- match.arg(scale_name)
  
  df_plot <- df %>%
    filter(Variable == var_label)
  
  if (!is.null(species)) {
    df_plot <- df_plot %>% filter(Species %in% species)
  }
  
  df_plot <- filter_plot_window(df_plot, plot_year_min, plot_year_max)
  
  if (nrow(df_plot) == 0) {
    warning("No data available for plot: ", var_label, " (", scale_name, ")")
    return(NULL)
  }
  
  p <- ggplot(df_plot, aes(Date, Value, color = Species)) +
    geom_line(linewidth = ifelse(scale_name == "yearly", 1, 0.8), alpha = 0.9) +
    scale_color_manual(values = cb_palette, drop = FALSE) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(
      title = paste(tools::toTitleCase(scale_name), "raw", var_label),
      x = "Time",
      y = ylab,
      color = NULL
    ) +
    theme_clean() +
    theme(legend.position = "top")
  
  if (scale_name == "yearly") {
    p <- p + geom_point(size = 1.6)
  }
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(
        filename_stub, "_", scale_name, "_raw_", species_tag(species), "_",
        ifelse(is.null(plot_year_min), "start", plot_year_min), "_",
        ifelse(is.null(plot_year_max), "end", plot_year_max),
        ".png"
      )
    ),
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )
  
  p
}

plot_aggregated_ts <- function(df,
                               cb_palette,
                               out_dir,
                               var_label,
                               ylab,
                               filename_stub,
                               scale = c("monthly", "yearly"),
                               stat = c("mean", "min", "max"),
                               species = NULL,
                               plot_year_min = NULL,
                               plot_year_max = NULL,
                               width = 11,
                               height = 6) {
  scale <- match.arg(scale)
  stat  <- match.arg(stat)
  
  lab_scale <- scale_label(scale)
  
  df_plot <- df %>%
    filter(
      Variable == var_label,
      Scale == lab_scale,
      Stat == stat
    )
  
  if (!is.null(species)) {
    df_plot <- df_plot %>% filter(Species %in% species)
  }
  
  df_plot <- filter_plot_window(df_plot, plot_year_min, plot_year_max)
  
  if (nrow(df_plot) == 0) {
    warning("No data available for aggregated plot: ", var_label, " / ", scale, " / ", stat)
    return(NULL)
  }
  
  p <- ggplot(df_plot, aes(Date, Value, color = Species)) +
    geom_line(linewidth = 0.8, alpha = 0.9) +
    scale_color_manual(values = cb_palette, drop = FALSE) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(
      title = paste(lab_scale, stat, var_label),
      x = "Time",
      y = ylab,
      color = NULL
    ) +
    theme_clean() +
    theme(legend.position = "top")
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(
        filename_stub, "_", scale, "_", stat, "_", species_tag(species), "_",
        ifelse(is.null(plot_year_min), "start", plot_year_min), "_",
        ifelse(is.null(plot_year_max), "end", plot_year_max),
        ".png"
      )
    ),
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )
  
  p
}