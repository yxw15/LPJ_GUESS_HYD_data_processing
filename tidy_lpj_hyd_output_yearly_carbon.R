# ==============================================================================
# 0. SETUP & PATHS
# ==============================================================================
library(tidyverse)

# Base directory setup
base_dir <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD"
setwd(base_dir)

if(!dir.exists("lpj_guess")) {
  dir.create("lpj_guess", recursive = TRUE)
}

# ==============================================================================
# 1. METADATA REGISTRIES & DYNAMIC READER
# ==============================================================================
species_map <- tribble(
  ~species, ~colname,
  "Beech",  "Fag_syl",
  "Oak",    "Que_rob",
  "Pine",   "Pin_syl",
  "Spruce", "Pic_abi"
)

yearly_files <- list(
  "cmass.out"      = "cmass",
  "cmass_mort.out" = "cmass_mort",
  "agpp.out"       = "agpp",
  "anpp.out"       = "anpp"
)

treatments <- c("control", "drought")

# UPDATED: Removed hard-coded year filters to process the full file content
read_lpj_yearly_var <- function(file_path, col_name, new_var_name) {
  if(!file.exists(file_path)) return(NULL)
  
  read.table(file_path, header = TRUE, check.names = FALSE) %>%
    rename_with(~ "year", matches("^Year$", ignore.case = TRUE)) %>%
    select(year, value = all_of(col_name)) %>%
    rename(!!new_var_name := value)
}

# ==============================================================================
# 2. INGESTION & CO-REGISTRATION PIPELINE
# ==============================================================================
final_yearly_carbon <- map_dfr(treatments, function(trt) {
  
  species_data <- species_map %>%
    group_by(species) %>%
    group_modify(~ {
      sp_path <- file.path("results_lpj/results_hoelstein_hyd", tolower(trt), .y$species)
      
      file_data_list <- map2(names(yearly_files), yearly_files, function(file, var_name) {
        read_lpj_yearly_var(file.path(sp_path, file), .x$colname, var_name)
      }) %>% keep(~ !is.null(.x))
      
      if (length(file_data_list) == 0) return(tibble())
      
      file_data_list %>% reduce(full_join, by = "year")
    }) %>% ungroup()
  
  if (nrow(species_data) > 0) species_data <- species_data %>% mutate(treatment = trt)
  return(species_data)
})

# ==============================================================================
# 3. EXPORT & DYNAMIC SUMMARY
# ==============================================================================
output_csv_path <- "lpj_guess/lpj_guess_hyd/lpj_control_drought_yearly_carbon_productivity_full.csv"
write_csv(final_yearly_carbon, output_csv_path)

# Dynamically calculate the range for the summary message
yr_range <- range(final_yearly_carbon$year, na.rm = TRUE)

cat(glue::glue(
  "\u2713 Success: Yearly data combined for range {yr_range[1]}-{yr_range[2]}!
    - Total consolidated rows: {nrow(final_yearly_carbon)}
    - Output saved to:         {output_csv_path}\n\n"
))

# Dynamically summarized print
cat("=== Mean Productivity & Biomass Stocks (Full Range) ===\n")
final_yearly_carbon %>%
  group_by(treatment, species) %>%
  summarise(across(c(cmass, cmass_mort, agpp, anpp), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  print(n = Inf)