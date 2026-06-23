library(tidyverse)
library(lubridate)

# ============================================================
# 1. Global Settings & Directories
# ============================================================
base_lpj_dir <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD/results_lpj"
figures_output_dir <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD/Figures/Beech_mortality"
csv_dir <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD/results_lpj/results_csv"

dir.create(figures_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

experimental_runs <- tribble(
  ~input_folder,
  "Beech_mortality"
)

plot_year_min <- 2000 
plot_year_max <- 2022
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    plot.title        = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5, size = 12),
    axis.title        = element_text(size = 13),
    axis.text         = element_text(size = 11),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    strip.text        = element_text(size = 12, face = "bold"),
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_line(color = "grey92", linewidth = 0.25),
    panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )

# ============================================================
# 2. Variable Registry
# ============================================================
vars_to_plot <- tribble(
  ~file, ~type, ~processing_mode, ~var_label, ~ylab, ~stem,
  "anpp.out", "annual", "annual_raw", "Net Primary Production", "NPP (kg C m-2 yr-1)", "anpp",
  "agpp.out", "annual", "annual_raw", "Gross Primary Production", "GPP (kg C m-2 yr-1)", "agpp",
  "cmass.out", "annual", "annual_raw", "Carbon Biomass", "Biomass (kg C m-2)", "cmass",
  "cmass_mort.out", "annual", "annual_raw", "Mortality Biomass Loss", "Biomass (kg C m-2)", "cmass_mort",
  "diam.out", "annual", "annual_raw", "Mean Diameter", "Diameter (cm)", "diam",
  "dbh.out", "annual", "annual_raw", "DBH", "Diameter (cm)", "dbh",
  "total_basal_area.out", "annual", "annual_raw", "Total Basal Area", "Basal Area (m²)", "total_basal_area",
  "total_basal_area_fraction.out", "annual", "annual_raw", "Basal Area (m² ha⁻¹)", "Basal Area (m² ha⁻¹) ", "total_basal_area_fraction",
  "landcover_frac.out", "annual", "annual_raw", "Forest Fraction", "Fraction", "forest_sum",
  "mort.out", "daily", "last_day_of_year", "Total Mortality", "Mortality", "mort",
  "mort_cav.out", "daily", "last_day_of_year", "Hydraulic Mortality", "Mortality", "mort_cav_day",
  "mort_min.out", "daily", "last_day_of_year", "Background Mortality", "Mortality", "mort_min",
  "mort_greff.out", "daily", "last_day_of_year", "Growth Efficiency Mort.", "Mortality", "mort_greff",
  "kappa_s_min.out", "daily", "last_day_of_year", "Fraction Cavitated Xylem", "Fraction", "kappa_s_min"
)

# ============================================================
# 3. Spatial Reference Loading
# ============================================================
Beech_sites <- read.csv(file.path(base_lpj_dir, "Beech_mortality/Beech_sites_coordinates.csv"))
sites_ready <- Beech_sites %>%
  mutate(Lon_round = round(lon, 2), Lat_round = round(lat, 2)) %>%
  select(short_name, site_name, Lon_round, Lat_round)

# ============================================================
# 4. Universal Pipeline Data Processor
# ============================================================
process_lpj_file <- function(filepath, sites_df, stem_name, file_type, processing_mode) {
  if (!file.exists(filepath)) return(NULL)
  
  df <- read.table(filepath, header = TRUE, check.names = FALSE)
  df_matched <- df %>%
    left_join(sites_df, by = c("Lon" = "Lon_round", "Lat" = "Lat_round")) %>%
    filter(!is.na(short_name))
  
  if (nrow(df_matched) == 0) return(NULL)
  if ("Year" %in% names(df_matched)) df_matched <- rename(df_matched, year = Year)
  if ("Day" %in% names(df_matched))  df_matched <- rename(df_matched, day = Day)
  
  if (stem_name == "forest_sum" && "Forest_sum" %in% names(df_matched)) {
    value_col <- "Forest_sum"
  } else if ("Fag_syl" %in% names(df_matched)) {
    value_col <- "Fag_syl"
  } else {
    num_cols <- names(df_matched)[sapply(df_matched, is.numeric)]
    value_col <- tail(setdiff(num_cols, c("year", "day", "Lon", "Lat", "Total", "Natural_sum")), 1)
  }
  
  if (length(value_col) == 0 || is.na(value_col)) return(NULL)
  
  if (file_type == "daily") {
    daily_df <- df_matched %>%
      mutate(year_padded = sprintf("%04d", year), 
             Date = as.Date(paste0(year_padded, "-01-01")) + day) %>%
      select(short_name, site_name, year, day, Date, value = all_of(value_col)) %>%
      mutate(variable = stem_name, type = "daily")
    
    if (processing_mode == "last_day_of_year") {
      return(daily_df %>% group_by(short_name, site_name, year, variable) %>% slice_tail(n = 1) %>% ungroup() %>% mutate(type = "annual", day = NA, Date = as.Date(NA)))
    } else {
      return(daily_df)
    }
  } else {
    return(df_matched %>% select(short_name, site_name, year, value = all_of(value_col)) %>% mutate(variable = stem_name, type = "annual", day = NA, Date = as.Date(NA)))
  }
}

# ============================================================
# 5. Execution Loop
# ============================================================
master_data_accumulator <- list()

for (run_idx in seq_len(nrow(experimental_runs))) {
  run_info <- experimental_runs[run_idx, ]
  current_in_dir <- file.path(base_lpj_dir, run_info$input_folder)
  
  for (i in seq_len(nrow(vars_to_plot))) {
    file_info <- vars_to_plot[i, ]
    df_processed <- process_lpj_file(file.path(current_in_dir, file_info$file), sites_ready, file_info$stem, file_info$type, file_info$processing_mode)
    
    if (is.null(df_processed)) next
    
    df_filtered <- df_processed %>% 
      filter(year >= plot_year_min, year <= plot_year_max) %>% 
      mutate(experiment = run_info$input_folder)
    
    if (nrow(df_filtered) == 0) next
    
    master_data_accumulator[[paste0(run_info$input_folder, ".", file_info$stem)]] <- df_filtered
    
    global_ymin <- max(0, min(df_filtered$value, na.rm = TRUE) - (diff(range(df_filtered$value, na.rm = TRUE)) * 0.05))
    global_ymax <- max(df_filtered$value, na.rm = TRUE) + (diff(range(df_filtered$value, na.rm = TRUE)) * 0.05)
    
    p <- ggplot(df_filtered, aes(x = year, y = value, color = short_name)) +
      geom_line(linewidth = 0.8) + geom_point(size = 1.5, alpha = 0.6) +
      scale_color_manual(values = cb_palette) +
      facet_wrap(~ site_name, ncol = 3) +
      coord_cartesian(ylim = c(global_ymin, global_ymax)) +
      labs(title = paste(file_info$var_label, "- by Site"), y = file_info$ylab, x = "Year") +
      base_theme + theme(legend.position = "none")
    
    ggsave(file.path(figures_output_dir, paste0(file_info$stem, "_site_panels.png")), p, width = 14, height = 10, dpi = 300)
  }
}

# ============================================================
# 5.5 Combined NPP and GPP Plotting
# ============================================================
npp_list <- master_data_accumulator[grepl("\\.anpp$", names(master_data_accumulator))]
gpp_list <- master_data_accumulator[grepl("\\.agpp$", names(master_data_accumulator))]

if (length(npp_list) > 0 && length(gpp_list) > 0) {
  combined_df <- bind_rows(c(npp_list, gpp_list))
  
  p_combined <- ggplot(combined_df, aes(x = year, y = value, color = variable, group = variable)) +
    geom_line() + geom_point() + facet_wrap(~ site_name, ncol = 3) +
    scale_color_manual(values = c("anpp" = cb_palette[1], "agpp" = cb_palette[2])) +
    base_theme + theme(legend.position = "bottom")
  
  ggsave(file.path(figures_output_dir, "comparison_anpp_agpp_same_panel.png"), p_combined, width = 14, height = 10, dpi = 300)
}

# ============================================================
# 5.6 Combined Mortality Components
# ============================================================
mort_keys <- unlist(lapply(experimental_runs$input_folder, function(f) paste0(f, ".", c("mort", "mort_cav_day", "mort_min", "mort_greff"))))
mort_keys <- mort_keys[mort_keys %in% names(master_data_accumulator)]

if (length(mort_keys) > 0) {
  mort_df <- bind_rows(master_data_accumulator[mort_keys]) %>%
    mutate(variable = factor(variable, 
                             levels = c("mort", "mort_cav_day", "mort_min", "mort_greff"),
                             labels = c("Total", "Hydraulic", "Background", "Growth-Eff.")))
  
  p_mort <- ggplot(mort_df, aes(x = year, y = value, color = variable)) +
    geom_line(linewidth = 0.9) + 
    geom_point(size = 2.5, alpha = 0.7) + 
    facet_wrap(~ site_name, ncol = 3) +
    # Using the top 4 colors from your palette
    scale_color_manual(values = cb_palette[1:4]) +
    labs(title = "Annual Mortality Components", y = "Mortality Rate", x = "Year", color = "Type") +
    base_theme + 
    theme(legend.position = "bottom")
  
  ggsave(file.path(figures_output_dir, "combined_mortality_components_ANNUAL.png"), p_mort, width = 14, height = 10, dpi = 300)
  
  write.csv(mort_df, file.path(csv_dir, "mortality_components_ANNUAL_timeseries.csv"), row.names = FALSE)
}

# ============================================================
# 5.7 Final Summary
# ============================================================
all_data_combined <- bind_rows(master_data_accumulator, .id = "source")
summary_stats <- all_data_combined %>% group_by(variable, site_name) %>% 
  summarise(Mean = mean(value, na.rm = TRUE), .groups = "drop")
write.csv(summary_stats, file.path(csv_dir, "all_variables_summary_stats.csv"), row.names = FALSE)

message("\n=== PROCESSING COMPLETE ===")
message(paste("All figures saved to:", figures_output_dir))