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
# Process RUE station data once to use as a master climate reference
meteo_path <- "MeteoSwiss/MeteoSwiss_station/all_stations_RUE_replaced_daytime_drought.csv"

data_RUE <- read_csv(meteo_path) %>%
  filter(station_abbr == "RUE") %>%
  mutate(
    date   = as.Date(date),
    temp_C = temperature - 273.15,
    rh_100 = relative_humidity * 100,
    # Calculate Vapor Pressure Deficit
    vpd    = RHtoVPD(RH = rh_100, TdegC = temp_C)
  ) %>%
  select(date, temperature, temp_C, vpd, global_radiation, precipitation)

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
# 3. DATA LOADING & SPECIES MERGING
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
  dgcwater.out = "Gc",
  dpsileaf.out = "psiL",
  dpsisoil.out = "psiS",
  dpsixylem.out = "psiX"
)

final_daily <- species_map %>%
  group_by(species) %>%
  group_modify(~ {
    sp_name <- .y$species
    col_ref <- .x$colname
    sp_path <- file.path("results/hoelstein_drought", paste0(sp_name, "_hoelstein"))
    
    # Read all 4 physiological files for this species and join them
    map2(names(var_files), var_files, function(file, var_name) {
      read_lpj_var(file.path(sp_path, file), col_ref, var_name)
    }) %>%
      reduce(full_join, by = "date")
  }) %>%
  ungroup() %>%
  # Merge with Meteo data
  inner_join(data_RUE, by = "date") %>%
  # Clean up missing values across all hydraulic parameters
  drop_na(Gc, psiL, psiS, psiX)

# =========================================================
# 4. FILTERING & EXPORT
# =========================================================

# 1. Export Full Dataset
write_csv(final_daily, "lpj_guess/lpj_Gc_psiL_psiX_psiS_drought.csv")

# 2. Apply Climate Filter (defined for summer dry/sunny days)
final_daily_filtered <- final_daily %>%
  filter(
    temperature > 287.15,           # >14°C
    precipitation < 1,              # Dry days
    global_radiation > 150,         # Sufficient light
    vpd > 0.3,                      # Evaporative demand
    month(date) %in% 6:9            # June to September
  )

# 3. Export Filtered Dataset
write_csv(final_daily_filtered, "lpj_guess/lpj_Gc_psiL_psiX_psiS_drought_climate_filter.csv")

# 4. Summary Printout
cat(glue::glue(
  "Processing Complete:
   - Total records: {nrow(final_daily)}
   - Filtered records (Summer/Dry): {nrow(final_daily_filtered)}
   - Files saved to 'lpj_guess/' folder.\n"
))

