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

# Define the variable mapping:
#   file      = LPJ output filename
#   var_name  = column name in the final merged dataframe
#   col_name  = column name in the raw .out file
#               "species" = use species-specific column (e.g. Fag_syl, Que_rob)
#               otherwise = literal column name (e.g. "det_plant" for stand-level outputs)
var_files <- tribble(
  ~file,              ~var_name,        ~col_name,
  "det_plant.out",    "ET_plant",       "det_plant",
  "det_total.out",    "ET_total",       "det_total",
  "dgc.out",          "Gc",             "species",
  "dpsileaf.out",     "psi_leaf",       "species",
  "dpsisoil.out",     "psi_soil",       "species",
  "dpsixylem.out",    "psi_xylem",      "species",
  "hydraulic_lag.out","hydraulic_lag",  "species",
  "kappa_s_min.out",  "kappy_s_min",    "species",
  "mort.out",         "mort",           "species",
  "mort_cav.out",     "mort_cav",       "species",
  "mort_greff.out",   "mort_greff",     "species",
  "mort_min.out",     "mort_min",       "species",
  "stem_diameter.out","stem_diameter",  "species",
  "stem_rwc.out",     "stem_rwc",       "species",
  "twd.out",          "twd",            "species"
)

# Vector of treatments to process
treatments <- c("control", "drought")

# Loop over treatments and species to build the unified dataset
final_daily <- map_dfr(treatments, function(trt) {

  # Get the meteo data for this treatment
  meteo_data <- meteo_data_list[[trt]]

  # Loop over each species
  species_data <- map_dfr(seq_len(nrow(species_map)), function(i) {
    sp_name <- species_map$species[i]
    col_ref <- species_map$colname[i]

    # Dynamically switch folder between 'control' and 'drought'
    sp_path <- file.path("results_lpj/results_hoelstein_stem_storage", trt, sp_name)

    # Read all physiological files for this species and join them
    var_data_list <- map(seq_len(nrow(var_files)), function(j) {
      fpath <- file.path(sp_path, var_files$file[j])
      if (!file.exists(fpath)) return(NULL)

      # Resolve column: "species" = use species-specific column (e.g. Fag_syl),
      # otherwise use the literal column name from var_files (e.g. "det_plant")
      read_col <- if (var_files$col_name[j] == "species") col_ref else var_files$col_name[j]

      read.table(fpath, header = TRUE) %>%
        convert_lpj_date() %>%
        select(date, value = all_of(read_col)) %>%
        rename(!!var_files$var_name[j] := value)
    }) %>%
      keep(~ !is.null(.x))  # Safety check in case a file is missing

    # If no files were found for this species, skip it
    if (length(var_data_list) == 0) {
      cat(glue::glue("  WARNING: No output files found for {sp_name} ({trt}) in {sp_path}\n"))
      return(NULL)
    }

    # Join all variables on date
    reduce(var_data_list, full_join, by = "date") %>%
      mutate(species = sp_name)
  })

  # Merge species data with treatment-specific meteo data, then tag treatment
  if (nrow(species_data) == 0) {
    cat(glue::glue("  WARNING: No data for treatment '{trt}' — skipping\n"))
    return(NULL)
  }

  species_data %>%
    mutate(treatment = trt) %>%
    inner_join(meteo_data, by = "date") %>%
    select(treatment, species, date, everything())

}) 

# =========================================================
# 4. FILTERING & EXPORT WITH UPDATED FILE NAMES
# =========================================================

# 1. Define descriptive file names based on included contents (all var_files variables)
output_name_full <- paste0(
  "lpj_guess/lpj_guess_stem_storage/lpj_control_drought_",
  paste(var_files$var_name, collapse = "_"),
  ".csv"
)
output_name_filtered <- paste0(
  "lpj_guess/lpj_guess_stem_storage/lpj_control_drought_",
  paste(var_files$var_name, collapse = "_"),
  "_climate_filter.csv"
)

# Create directory if it doesn't exist
if(!dir.exists("lpj_guess/lpj_guess_stem_storage")) {
  dir.create("lpj_guess/lpj_guess_stem_storage", recursive = TRUE)
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

