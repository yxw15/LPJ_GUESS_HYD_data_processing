library(tidyverse)
library(lubridate)

# ============================================================
# 1. Global Settings & Directories
# ============================================================
base_lpj_dir <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD/results_lpj"
figures_output_dir <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD/Figures/Pine_Oak_mortality"
csv_dir <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD/results_lpj/results_csv"

dir.create(figures_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

experimental_runs <- tribble(
  ~input_folder, ~species, ~species_col,
  "Oak_mortality_visp", "Oak", "Que_pub",
  "Pine_mortality_visp", "Pine", "Pin_syl"
)

plot_year_min <- 1997 
plot_year_max <- 2004
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Create sequence for x-axis breaks (every 2 years starting from 1998)
x_axis_breaks <- seq(1998, plot_year_max, by = 2)  # Changed from 1997 to 1998
# Set x-axis limits to exactly from min to max
x_axis_limits <- c(plot_year_min, plot_year_max)

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
    panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    legend.position   = "bottom"
  )

# ============================================================
# 2. Variable Registry (updated processing_mode)
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
  "total_basal_area_fraction.out", "annual", "annual_raw", "Basal Area (m² ha⁻¹)", "Basal Area (m² ha⁻¹)", "total_basal_area_fraction",
  "landcover_frac.out", "annual", "annual_raw", "Forest Fraction", "Fraction", "forest_sum",
  "mort.out", "daily", "last_day_of_year", "Total Mortality", "Mortality (%)", "mort",
  "mort_cav.out", "daily", "last_day_of_year", "Hydraulic Mortality", "Mortality (%)", "mort_cav",
  "mort_min.out", "daily", "last_day_of_year", "Background Mortality", "Mortality (%)", "mort_min",
  "mort_greff.out", "daily", "last_day_of_year", "Growth Efficiency Mort.", "Mortality (%)", "mort_greff",
  "kappa_s_min.out", "daily", "last_day_of_year", "Fraction Cavitated Xylem", "Fraction", "kappa_s_min"
)

# ============================================================
# 3. Spatial Reference Loading
# ============================================================
sites_ready <- data.frame(
  short_name = "Visp",
  site_name = "Visp",
  Lon_round = 7.84,
  Lat_round = 46.30,
  stringsAsFactors = FALSE
)

# ============================================================
# 4. Enhanced Data Processor (handles both annual and daily)
# ============================================================
process_lpj_file <- function(filepath, species_col, species_name, processing_mode, file_type) {
  if (!file.exists(filepath)) {
    message(paste("File not found:", filepath))
    return(NULL)
  }
  
  # Read the file
  df <- read.table(filepath, header = TRUE, check.names = FALSE)
  
  # Check if required columns exist
  required_cols <- c("Lon", "Lat", species_col)
  if (file_type == "annual") {
    required_cols <- c(required_cols, "Year")
  } else if (file_type == "daily") {
    required_cols <- c(required_cols, "Year", "Day")
  }
  
  if (!all(required_cols %in% names(df))) {
    message(paste("Required columns not found in:", filepath))
    message(paste("Available columns:", paste(names(df), collapse = ", ")))
    return(NULL)
  }
  
  # Round coordinates
  df <- df %>%
    mutate(Lon_round = round(Lon, 2),
           Lat_round = round(Lat, 2))
  
  # Join with site information
  df_matched <- df %>%
    left_join(sites_ready, by = c("Lon_round", "Lat_round")) %>%
    filter(!is.na(short_name))
  
  if (nrow(df_matched) == 0) {
    message(paste("No matching coordinates found in:", filepath))
    return(NULL)
  }
  
  # Process based on file type and mode
  if (file_type == "annual") {
    # For annual files, just extract the annual values
    result <- df_matched %>%
      select(short_name, site_name, year = Year, value = all_of(species_col)) %>%
      mutate(species = species_name,
             year = as.numeric(year))
  } 
  else if (file_type == "daily" && processing_mode == "last_day_of_year") {
    # For daily files where we want the last day of each year
    # First, get the maximum day for each year
    result <- df_matched %>%
      group_by(short_name, site_name, Year) %>%
      filter(Day == max(Day, na.rm = TRUE)) %>%  # Get last day of year
      ungroup() %>%
      select(short_name, site_name, year = Year, value = all_of(species_col)) %>%
      mutate(species = species_name,
             year = as.numeric(year))
    
    message(paste("  Extracted last day values for", species_name, "-", nrow(result), "years"))
  }
  
  return(result)
}

# ============================================================
# 5. Collect Data for Both Species
# ============================================================
all_data <- list()

# Process each species and each variable
for (run_idx in seq_len(nrow(experimental_runs))) {
  run_info <- experimental_runs[run_idx, ]
  current_dir <- file.path(base_lpj_dir, run_info$input_folder)
  
  message(paste("\nProcessing", run_info$species, "from:", current_dir))
  
  for (i in seq_len(nrow(vars_to_plot))) {
    file_info <- vars_to_plot[i, ]
    file_path <- file.path(current_dir, file_info$file)
    
    # Process the file with appropriate method
    df_processed <- process_lpj_file(
      file_path, 
      run_info$species_col, 
      run_info$species,
      file_info$processing_mode,
      file_info$type
    )
    
    if (!is.null(df_processed)) {
      df_filtered <- df_processed %>% 
        filter(year >= plot_year_min, year <= plot_year_max) %>%
        mutate(variable = file_info$stem,
               var_label = file_info$var_label,
               ylab = file_info$ylab)
      
      if (nrow(df_filtered) > 0) {
        all_data[[paste0(run_info$species, "_", file_info$stem)]] <- df_filtered
        message(paste("  Loaded:", file_info$stem, "-", nrow(df_filtered), "rows"))
      }
    }
  }
}

# Combine all data
if (length(all_data) == 0) {
  stop("No data was loaded. Please check file paths and column names.")
}

combined_data <- bind_rows(all_data)
message(paste("\nTotal loaded:", nrow(combined_data), "rows"))

# Verify each variable has correct number of years per species
verification <- combined_data %>%
  group_by(variable, species) %>%
  summarise(years = n_distinct(year), .groups = "drop")
print(verification)

# ============================================================
# 6. Individual Variable Plots (Oak vs Pine side by side)
# ============================================================
for (var_name in unique(combined_data$variable)) {
  var_subset <- combined_data %>% filter(variable == var_name)
  
  if (nrow(var_subset) == 0) next
  
  # Get variable metadata
  var_info <- vars_to_plot %>% filter(stem == var_name)
  var_label <- ifelse(nrow(var_info) > 0, var_info$var_label[1], var_name)
  y_label <- ifelse(nrow(var_info) > 0, var_info$ylab[1], var_name)
  
  # Create comparison plot with x-axis breaks every 2 years (1998, 2000, 2002, 2004)
  p <- ggplot(var_subset, aes(x = year, y = value, color = species)) +
    geom_line(linewidth = 1.2) + 
    geom_point(size = 2, alpha = 0.7) +
    scale_x_continuous(breaks = x_axis_breaks, limits = x_axis_limits, expand = c(0.02, 0)) +
    scale_color_manual(values = c("Oak" = cb_palette[1], "Pine" = cb_palette[2])) +
    labs(title = paste(var_label, "- Visp Site"),
         subtitle = "Oak vs Pine Comparison (1997-2004)",
         y = y_label,
         x = "Year",
         color = "Species") +
    base_theme +
    theme(legend.position = "bottom")
  
  ggsave(file.path(figures_output_dir, paste0(var_name, "_Oak_vs_Pine_Visp.png")), 
         p, width = 10, height = 6, dpi = 300)
  
  message(paste("Saved plot for:", var_name))
}

# ============================================================
# 7. Mortality Components Plot 
# ============================================================
mort_vars <- c("mort", "mort_cav", "mort_min", "mort_greff")
mort_data <- combined_data %>% filter(variable %in% mort_vars)

if (nrow(mort_data) > 0) {
  # Check if we have annual data (one point per year)
  mort_check <- mort_data %>%
    group_by(variable, species, year) %>%
    summarise(n_records = n(), .groups = "drop")
  
  if (any(mort_check$n_records > 1)) {
    message("WARNING: Multiple records per year detected in mortality data!")
    print(mort_check %>% filter(n_records > 1))
  }
  
  mort_data <- mort_data %>%
    mutate(mort_type = factor(variable, 
                              levels = mort_vars,
                              labels = c("Total Mortality", "Hydraulic Mortality", 
                                         "Background Mortality", "Growth-Efficiency Mortality")))
  
  # Define colors for different mortality types
  mortality_colors <- c(
    "Total Mortality" = "#D55E00",
    "Hydraulic Mortality" = "#CC79A7",
    "Background Mortality" = "#009E73",
    "Growth-Efficiency Mortality" = "#F0E442"
  )
  
  # Define line types for species (Oak = solid, Pine = dashed)
  species_linetypes <- c("Oak" = "solid", "Pine" = "dashed")
  
  # Faceted by mortality type (one panel per mortality type)
  p_mort_facet <- ggplot(mort_data, aes(x = year, y = value, color = mort_type, linetype = species)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5, alpha = 0.8, aes(shape = species)) +
    facet_wrap(~ mort_type, ncol = 2, scales = "free_y") +
    scale_x_continuous(breaks = x_axis_breaks, limits = x_axis_limits, expand = c(0.02, 0)) +
    scale_color_manual(values = mortality_colors) +
    scale_linetype_manual(values = species_linetypes) +
    scale_shape_manual(values = c("Oak" = 16, "Pine" = 17)) +
    labs(title = "Mortality Components - Visp Site",
         subtitle = "Annual values (last day of year from daily output) - Oak vs Pine (1997-2004)",
         y = "Mortality (%)",
         x = "Year",
         color = "Mortality Type",
         linetype = "Species",
         shape = "Species") +
    base_theme +
    theme(strip.text = element_text(size = 11, face = "bold"))
  
  ggsave(file.path(figures_output_dir, "mortality_components_Oak_Pine_Visp.png"), 
         p_mort_facet, width = 12, height = 10, dpi = 300)
  
  # All mortality types in one plot (updated version with colors for types and linetypes for species)
  p_mort_combined <- ggplot(mort_data, aes(x = year, y = value, color = mort_type, linetype = species)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2, alpha = 0.7, aes(shape = species)) +
    scale_x_continuous(breaks = x_axis_breaks, limits = x_axis_limits, expand = c(0.02, 0)) +
    scale_color_manual(values = mortality_colors) +
    scale_linetype_manual(values = species_linetypes) +
    scale_shape_manual(values = c("Oak" = 16, "Pine" = 17)) +
    labs(title = "Mortality Components - Visp Site",
         subtitle = "Oak (solid) vs Pine (dashed) - Different colors represent mortality types (1997-2004)",
         y = "Mortality (%)",
         x = "Year",
         color = "Mortality Type",
         linetype = "Species",
         shape = "Species") +
    base_theme +
    theme(legend.box = "vertical",
          legend.title.align = 0.5)
  
  ggsave(file.path(figures_output_dir, "mortality_components_combined_Oak_Pine_Visp.png"), 
         p_mort_combined, width = 12, height = 7, dpi = 300)
  
  # Save mortality data
  write.csv(mort_data, file.path(csv_dir, "mortality_components_Visp_Oak_Pine.csv"), row.names = FALSE)
  
  # Print summary of mortality data
  message("\nMortality Data Summary:")
  print(mort_data %>%
          group_by(species, mort_type) %>%
          summarise(
            Mean_mortality = mean(value, na.rm = TRUE),
            Max_mortality = max(value, na.rm = TRUE),
            Year_of_max = year[which.max(value)],
            .groups = "drop"
          ))
}

# ============================================================
# 8. NPP vs GPP Comparison
# ============================================================
npp_data <- combined_data %>% filter(variable == "anpp")
gpp_data <- combined_data %>% filter(variable == "agpp")

if (nrow(npp_data) > 0 && nrow(gpp_data) > 0) {
  combined_flux <- bind_rows(
    npp_data %>% mutate(flux_type = "NPP"),
    gpp_data %>% mutate(flux_type = "GPP")
  )
  
  p_flux <- ggplot(combined_flux, aes(x = year, y = value, color = flux_type, linetype = species)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2, alpha = 0.7) +
    scale_x_continuous(breaks = x_axis_breaks, limits = x_axis_limits, expand = c(0.02, 0)) +
    scale_color_manual(values = c("NPP" = cb_palette[1], "GPP" = cb_palette[2])) +
    labs(title = "NPP vs GPP Comparison - Visp Site",
         subtitle = "Oak and Pine (1997-2004)",
         y = "Carbon Flux (kg C m-2 yr-1)",
         x = "Year",
         color = "Flux Type",
         linetype = "Species") +
    base_theme
  
  ggsave(file.path(figures_output_dir, "NPP_GPP_comparison_Oak_Pine_Visp.png"), 
         p_flux, width = 10, height = 6, dpi = 300)
}

# ============================================================
# 9. Additional Key Variables Comparison
# ============================================================
# Biomass comparison
biomass_data <- combined_data %>% filter(variable == "cmass")
if (nrow(biomass_data) > 0) {
  p_biomass <- ggplot(biomass_data, aes(x = year, y = value, color = species)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2, alpha = 0.7) +
    scale_x_continuous(breaks = x_axis_breaks, limits = x_axis_limits, expand = c(0.02, 0)) +
    scale_color_manual(values = c("Oak" = cb_palette[1], "Pine" = cb_palette[2])) +
    labs(title = "Carbon Biomass - Visp Site",
         subtitle = "Oak vs Pine (1997-2004)",
         y = "Biomass (kg C m-2)",
         x = "Year",
         color = "Species") +
    base_theme
  
  ggsave(file.path(figures_output_dir, "biomass_comparison_Oak_Pine_Visp.png"), 
         p_biomass, width = 10, height = 6, dpi = 300)
}

# Diameter/DBH comparison
diam_data <- combined_data %>% filter(variable %in% c("diam", "dbh"))
if (nrow(diam_data) > 0) {
  p_diam <- ggplot(diam_data, aes(x = year, y = value, color = species, linetype = variable)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2, alpha = 0.7) +
    scale_x_continuous(breaks = x_axis_breaks, limits = x_axis_limits, expand = c(0.02, 0)) +
    scale_color_manual(values = c("Oak" = cb_palette[1], "Pine" = cb_palette[2])) +
    labs(title = "Diameter Growth - Visp Site",
         subtitle = "Oak vs Pine (1997-2004)",
         y = "Diameter (cm)",
         x = "Year",
         color = "Species",
         linetype = "Variable") +
    base_theme
  
  ggsave(file.path(figures_output_dir, "diameter_comparison_Oak_Pine_Visp.png"), 
         p_diam, width = 10, height = 6, dpi = 300)
}

# Xylem cavitation
kappa_data <- combined_data %>% filter(variable == "kappa_s_min")
if (nrow(kappa_data) > 0) {
  p_kappa <- ggplot(kappa_data, aes(x = year, y = value, color = species)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2, alpha = 0.7) +
    scale_x_continuous(breaks = x_axis_breaks, limits = x_axis_limits, expand = c(0.02, 0)) +
    scale_color_manual(values = c("Oak" = cb_palette[1], "Pine" = cb_palette[2])) +
    labs(title = "Xylem Cavitation - Visp Site",
         subtitle = "Fraction of Cavitated Xylem (Oak vs Pine) (1997-2004)",
         y = "Fraction Cavitated",
         x = "Year",
         color = "Species") +
    base_theme
  
  ggsave(file.path(figures_output_dir, "xylem_cavitation_Oak_Pine_Visp.png"), 
         p_kappa, width = 10, height = 6, dpi = 300)
}

# ============================================================
# 10. Summary Statistics
# ============================================================
summary_stats <- combined_data %>% 
  group_by(variable, species) %>% 
  summarise(
    Mean = mean(value, na.rm = TRUE),
    SD = sd(value, na.rm = TRUE),
    Min = min(value, na.rm = TRUE),
    Max = max(value, na.rm = TRUE),
    CV = (SD/Mean)*100,
    .groups = "drop"
  )

write.csv(summary_stats, file.path(csv_dir, "summary_stats_Visp_Oak_Pine.csv"), row.names = FALSE)

# Save complete dataset
write.csv(combined_data, file.path(csv_dir, "complete_dataset_Visp_Oak_Pine.csv"), row.names = FALSE)

# ============================================================
# 11. Create a summary plot with all key variables
# ============================================================
key_vars_for_summary <- c("anpp", "cmass", "mort", "kappa_s_min")
summary_plot_data <- combined_data %>% filter(variable %in% key_vars_for_summary)

if (nrow(summary_plot_data) > 0) {
  # Get proper labels
  summary_plot_data <- summary_plot_data %>%
    mutate(var_label = case_when(
      variable == "anpp" ~ "NPP",
      variable == "cmass" ~ "Biomass",
      variable == "mort" ~ "Total Mortality",
      variable == "kappa_s_min" ~ "Xylem Cavitation",
      TRUE ~ variable
    ))
  
  p_summary <- ggplot(summary_plot_data, aes(x = year, y = value, color = species)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5, alpha = 0.6) +
    facet_wrap(~ var_label, scales = "free_y", ncol = 2) +
    scale_x_continuous(breaks = x_axis_breaks, limits = x_axis_limits, expand = c(0.02, 0)) +
    scale_color_manual(values = c("Oak" = cb_palette[1], "Pine" = cb_palette[2])) +
    labs(title = "Key Ecosystem Variables - Visp Site",
         subtitle = "Oak vs Pine Comparison (1997-2004)",
         y = "Value",
         x = "Year",
         color = "Species") +
    base_theme +
    theme(strip.text = element_text(size = 11, face = "bold"))
  
  ggsave(file.path(figures_output_dir, "summary_key_variables_Oak_Pine_Visp.png"), 
         p_summary, width = 12, height = 10, dpi = 300)
}

# ============================================================
# 12. Print completion message
# ============================================================
message("\n=== PROCESSING COMPLETE ===")
message(paste("All figures saved to:", figures_output_dir))
message(paste("Data saved to:", csv_dir))
message("\nVariables successfully loaded:")
for (var in unique(combined_data$variable)) {
  n_species <- length(unique(combined_data$species[combined_data$variable == var]))
  message(paste("  -", var, "(", n_species, "species)"))
}