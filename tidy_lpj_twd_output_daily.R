# =========================================================
# 0. SETUP & PATHS
# =========================================================
library(tidyverse)
library(lubridate)
library(plantecophys)

# Base directory
base_dir <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD"
setwd(base_dir)

# =========================================================
# 1. METEO DATA PROCESSING
# =========================================================
# Function to read meteo data for a specific treatment
read_meteo_for_treatment <- function(treatment) {
  # Define path based on treatment
  meteo_path <- ifelse(treatment == "control",
                       "MeteoSwiss/MeteoSwiss_station/all_stations_RUE_replaced_daytime_control.csv",
                       "MeteoSwiss/MeteoSwiss_station/all_stations_RUE_replaced_daytime_drought.csv"
  )
  
  read_csv(meteo_path) %>%
    filter(station_abbr == "RUE") %>%
    mutate(
      date   = as.Date(date),
      temp_C = temperature - 273.15,
      rh_100 = relative_humidity * 100,
      # Calculate Vapor Pressure Deficit
      vpd    = RHtoVPD(RH = rh_100, TdegC = temp_C)
    ) %>%
    select(date, temperature, temp_C, vpd, global_radiation, precipitation)
}

# Alternative: Read both meteo files into a named list for efficiency
meteo_data_list <- list(
  control = read_meteo_for_treatment("control"),
  drought = read_meteo_for_treatment("drought")
)

# =========================================================
# 2. HELPERS
# =========================================================
# Convert LPJ Day/Year to Date object
convert_lpj_date <- function(df) {
  df %>% mutate(date = as.Date(Day, origin = paste0(Year, "-01-01")))
}

# Generic reader to fetch specific columns from LPJ output files
read_lpj_var <- function(file_path, col_name, new_var_name) {
  if(!file.exists(file_path)) return(NULL)
  
  read.table(file_path, header = TRUE) %>%
    convert_lpj_date() %>%
    select(date, value = all_of(col_name)) %>%
    rename(!!new_var_name := value)
}

# =========================================================
# 3. DATA LOADING & SPECIES/TREATMENT MERGING
# =========================================================
species_map <- tribble(
  ~species, ~colname,
  "Beech",  "Fag_syl",
  "Oak",    "Que_rob",
  "Pine",   "Pin_syl",
  "Spruce", "Pic_abi"
)

# Define the variable mapping: File suffix -> Column name in final DF
var_files <- list(
  gc.out = "Gc",
  psi_leaf.out = "psi_leaf",
  psi_soil.out = "psi_soil",
  psi_xylem.out = "psi_xylem",
  stem_diameter.out = "stem_diameter",
  twd.out = "twd",
  stem_rwc.out = "stem_rwc"
)

# Vector of treatments to process
treatments <- c("control", "drought")

# Loop over treatments and species to build the unified dataset
final_daily <- map_dfr(treatments, function(trt) {
  
  # Get the meteo data for this treatment
  meteo_data <- meteo_data_list[[trt]]
  
  species_data <- species_map %>%
    group_by(species) %>%
    group_modify(~ {
      sp_name <- .y$species
      col_ref <- .x$colname
      
      # Dynamically switch folder between 'control' and 'drought'
      sp_path <- file.path("results_lpj/results_hoelstein_twd", trt, sp_name)
      
      # Read all physiological files for this species and join them
      map2(names(var_files), var_files, function(file, var_name) {
        read_lpj_var(file.path(sp_path, file), col_ref, var_name)
      }) %>%
        keep(~ !is.null(.x)) %>% # Safety check in case a file is missing
        reduce(full_join, by = "date")
    }) %>%
    ungroup() %>%
    mutate(treatment = trt) # Add tracking column
  
  # Merge species data with treatment-specific meteo data
  species_data %>%
    inner_join(meteo_data, by = "date")
  
}) %>%
  # Clean up records where core hydraulic parameters are missing
  drop_na(Gc, psi_leaf, psi_soil, psi_xylem) %>%
  # Reorder columns to put tracking info first
  select(treatment, species, date, everything())

# =========================================================
# 4. FILTERING & EXPORT WITH UPDATED FILE NAMES
# =========================================================

# 1. Define descriptive file names based on included contents
output_name_full <- "lpj_guess/lpj_guess_twd/lpj_control_drought_Gc_psiL_psiX_psiS_stem_diameter_twd_stem_rwc.csv"
output_name_filtered <- "lpj_guess/lpj_guess_twd/lpj_control_drought_Gc_psiL_psiX_psiS_stem_diameter_twd_stem_rwc_climate_filter.csv"

# Create directory if it doesn't exist
if(!dir.exists("lpj_guess")) {
  dir.create("lpj_guess", recursive = TRUE)
}

# 2. Export Full Dataset
write_csv(final_daily, output_name_full)

# 3. Apply Climate Filter (defined for summer dry/sunny days)
final_daily_filtered <- final_daily %>%
  filter(
    temperature > 287.15,           # >14°C
    precipitation < 1,              # Dry days
    global_radiation > 150,         # Sufficient light
    vpd > 0.3,                      # Evaporative demand
    month(date) %in% 6:9            # June to September
  )

# 4. Export Filtered Dataset
write_csv(final_daily_filtered, output_name_filtered)

# 5. Summary Printout
cat(glue::glue(
  "Processing Complete:
    - Total records (All Treatments): {nrow(final_daily)}
    - Filtered records (Summer/Dry): {nrow(final_daily_filtered)}
    - Full file saved to: {output_name_full}
    - Filtered file saved to: {output_name_filtered}\n"
))

# Optional: Print summary statistics by treatment and species
cat("\n\n=== Summary by Treatment and Species ===\n")
final_daily %>%
  group_by(treatment, species) %>%
  summarise(
    n_days = n(),
    mean_Gc = mean(Gc, na.rm = TRUE),
    mean_psi_leaf = mean(psi_leaf, na.rm = TRUE),
    mean_vpd = mean(vpd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(n = Inf)

