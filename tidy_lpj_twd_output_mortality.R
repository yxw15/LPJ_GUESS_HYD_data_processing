# ==============================================================================
# 0. SETUP & PATHS
# ==============================================================================
library(tidyverse)
library(lubridate)

# Base directory setup
base_dir <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD"
setwd(base_dir)

if(!dir.exists("lpj_guess")) {
  dir.create("lpj_guess", recursive = TRUE)
}

# ==============================================================================
# 1. METADATA REGISTRIES & HELPERS
# ==============================================================================
species_map <- tribble(
  ~species, ~colname,
  "Beech",  "Fag_syl",
  "Oak",    "Que_rob",
  "Pine",   "Pin_syl",
  "Spruce", "Pic_abi"
)

var_files <- list(
  "mort.out"           = "mort",
  "mort_cav.out"       = "mort_cav",
  "mort_min.out"       = "mort_min",
  "mort_greff.out"     = "mort_greff",
  "kappa_s_min.out"    = "kappa_s_min",
  "kappa_s_today.out"  = "kappa_s_today"
)

treatments <- c("control", "drought")

# Helper to map LPJ outputs to true Date objects
convert_lpj_date <- function(df) {
  df %>% 
    rename_with(~ c("year", "day"), 3:4) %>% 
    mutate(date = as.Date(day, origin = paste0(year, "-01-01")))
}

# Generic table reader: Removed year filters to allow full temporal range
read_lpj_var <- function(file_path, col_name, new_var_name) {
  if(!file.exists(file_path)) return(NULL)
  
  read.table(file_path, header = TRUE, check.names = FALSE) %>%
    convert_lpj_date() %>%
    select(date, year, value = all_of(col_name)) %>%
    rename(!!new_var_name := value)
}

# ==============================================================================
# 2. INGESTION & CO-REGISTRATION PIPELINE
# ==============================================================================
final_mortality_daily <- map_dfr(treatments, function(trt) {
  
  species_data <- species_map %>%
    group_by(species) %>%
    group_modify(~ {
      sp_path <- file.path("results_lpj/results_hoelstein_twd", tolower(trt), .y$species)
      
      file_data_list <- map2(names(var_files), var_files, function(file, var_name) {
        read_lpj_var(file.path(sp_path, file), .x$colname, var_name)
      }) %>%
        keep(~ !is.null(.x)) 
      
      if (length(file_data_list) == 0) return(tibble())
      
      # Join on both date and year to ensure alignment across the full timeseries
      file_data_list %>% reduce(full_join, by = c("date", "year"))
    }) %>%
    ungroup()
  
  if (nrow(species_data) > 0) {
    species_data <- species_data %>% mutate(treatment = trt)
  }
  
  return(species_data)
}) 

if (nrow(final_mortality_daily) == 0) {
  stop("Pipeline halted: No data files found. Check your path definitions.")
}

final_mortality_daily <- final_mortality_daily %>%
  select(treatment, species, date, year, everything())

# ==============================================================================
# 3. ANNUAL LAST-DAY CALCULATIONS & EXPORT
# ==============================================================================

# A. Save the compiled raw daily timeseries 
write_csv(final_mortality_daily, "lpj_guess/lpj_guess_twd/lpj_control_drought_mort_kappa_daily_full.csv")

# B. Aggregate to pull the last day of each year across the entire available range
final_mortality_annual_last <- final_mortality_daily %>%
  group_by(treatment, species, year) %>%
  arrange(date, .by_group = TRUE) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  select(treatment, species, year, date, everything())

write_csv(final_mortality_annual_last, "lpj_guess/lpj_guess_twd/lpj_control_drought_mort_kappa_annual_full.csv")

cat("Success: Full series processed from", min(final_mortality_annual_last$year), 
    "to", max(final_mortality_annual_last$year), "\n")