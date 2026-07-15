# =========================================================
# 0. SETUP & PATHS (AUGUST VALIDATION ONLY)
# =========================================================
library(tidyverse)
library(lubridate)
library(plantecophys)

VALIDATION_MONTH <- 8  # August only

base_dir <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD"
setwd(base_dir)

# =========================================================
# 1. METEO DATA PROCESSING
# =========================================================
read_meteo_for_treatment <- function(treatment) {
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
      vpd    = RHtoVPD(RH = rh_100, TdegC = temp_C)
    ) %>%
    select(date, temperature, temp_C, vpd, global_radiation, precipitation)
}

meteo_data_list <- list(
  control = read_meteo_for_treatment("control"),
  drought = read_meteo_for_treatment("drought")
)

# =========================================================
# 2. HELPERS
# =========================================================
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

var_files <- tribble(
  ~file,              ~var_name,        ~col_name,
  "det_plant.out",    "ET",             "det_plant",
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

treatments <- c("control", "drought")

final_daily <- map_dfr(treatments, function(trt) {
  meteo_data <- meteo_data_list[[trt]]

  species_data <- map_dfr(seq_len(nrow(species_map)), function(i) {
    sp_name <- species_map$species[i]
    col_ref <- species_map$colname[i]
    sp_path <- file.path("results_lpj/results_hoelstein_stem_storage", trt, sp_name)

    var_data_list <- map(seq_len(nrow(var_files)), function(j) {
      fpath <- file.path(sp_path, var_files$file[j])
      if (!file.exists(fpath)) return(NULL)
      read_col <- if (var_files$col_name[j] == "species") col_ref else var_files$col_name[j]
      read.table(fpath, header = TRUE) %>%
        convert_lpj_date() %>%
        select(date, value = all_of(read_col)) %>%
        rename(!!var_files$var_name[j] := value)
    }) %>% keep(~ !is.null(.x))

    if (length(var_data_list) == 0) {
      cat(glue::glue("  WARNING: No output files found for {sp_name} ({trt}) in {sp_path}\n"))
      return(NULL)
    }
    reduce(var_data_list, full_join, by = "date") %>%
      mutate(species = sp_name)
  })

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
# 4. FILTERING: AUGUST ONLY + CLIMATE FILTER
# =========================================================

output_name_full <- paste0(
  "lpj_guess/lpj_guess_stem_storage/lpj_control_drought_",
  paste(var_files$var_name, collapse = "_"),
  "_august.csv"
)
output_name_filtered <- paste0(
  "lpj_guess/lpj_guess_stem_storage/lpj_control_drought_",
  paste(var_files$var_name, collapse = "_"),
  "_august_climate_filter.csv"
)

if(!dir.exists("lpj_guess/lpj_guess_stem_storage")) {
  dir.create("lpj_guess/lpj_guess_stem_storage", recursive = TRUE)
}

# Export Full Dataset (August only)
final_daily_aug <- final_daily %>% filter(month(date) == VALIDATION_MONTH)
write_csv(final_daily_aug, output_name_full)

# Climate Filter (August only: dry, sunny days)
final_daily_filtered <- final_daily %>%
  filter(
    month(date) == VALIDATION_MONTH,   # AUGUST ONLY
    temperature > 287.15,
    precipitation < 1,
    global_radiation > 150,
    vpd > 0.3
  )

write_csv(final_daily_filtered, output_name_filtered)

cat(glue::glue(
  "August Validation Processing Complete:
    - August records (All): {nrow(final_daily_aug)}
    - August filtered (Dry/Sunny): {nrow(final_daily_filtered)}
    - Full file: {output_name_full}
    - Filtered file: {output_name_filtered}\n"
))

cat("\n=== August Summary by Treatment and Species ===\n")
final_daily_aug %>%
  group_by(treatment, species) %>%
  summarise(
    n_days = n(),
    mean_Gc = mean(Gc, na.rm = TRUE),
    mean_psi_leaf = mean(psi_leaf, na.rm = TRUE),
    mean_vpd = mean(vpd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(n = Inf)
