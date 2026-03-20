setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

# ============================================================
# Settings
# ============================================================

base_dir <- "results"
out_dir  <- file.path(base_dir, "Figures", "times_series")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

day_starts_at_zero <- TRUE

year_min <- 1990
year_max <- 2022

plot_year_min <- 1990
plot_year_max <- 2022

stats_to_make <- c("mean", "min", "max")

species_map <- tribble(
  ~species,   ~colname,
  "Beech",    "Fag_syl",
  "Oak_pub",  "Que_pub",
  "Oak_rob",  "Que_rob",
  "Pine",     "Pin_syl",
  "Spruce",   "Pic_abi"
)

cb_palette <- c(
  Oak_pub  = "darkorange",
  Oak_rob  = "#F0E442",
  Beech    = "dodgerblue",
  Spruce   = "green4",
  Pine     = "purple1"
)

# ============================================================
# Variables
# ============================================================

daily_vars <- tribble(
  ~file,               ~var_label,                         ~ylab,                           ~stem,
  "dgc.out",           "Leaf conductance (gc)",            "Leaf conductance (mm/s)",       "dgc",
  "dgcwater.out",      "Canopy conductance (gcwater)",     "Canopy conductance (m/s)",      "dgcwater",
  "dpsisoil.out",      "Soil water potential",             "Soil water potential (MPa)",    "dpsisoil",
  "dpsileaf.out",      "Leaf water potential",             "Leaf water potential (MPa)",    "dpsileaf",
  "dpsixylem.out",     "Xylem water potential",            "Xylem water potential (MPa)",   "dpsixylem",
  "et_total.out",      "Transpiration",                    "Transpiration (mm/s)",          "et_total",
  "mort_cav_day.out",  "Hydraulic-failure mortality",      "Hydraulic mortality",           "mort_cav_day",
  "mort.out",          "Total mortality",                  "Total mortality",               "mort",
  "kappa_s_min.out",   "kappa_s_min",                      "Cavitation fraction",           "kappa_s_min",
  "kappa_s_today.out", "kappa_s_today",                    "Daily cavitation fraction",     "kappa_s_today"
)

monthly_vars <- tribble(
  ~file,       ~var_label,                     ~ylab,                          ~stem,
  "mevap.out", "Evapotranspiration",           "Evapotranspiration (mm/month)", "mevap",
  "mgpp.out",  "Gross Primary Production",     "GPP (kgC/m²/month)",            "mgpp",
  "mnpp.out",  "Net Primary Production",       "NPP (kgC/m²/month)",            "mnpp",
  "mlai.out",  "Leaf Area Index",              "LAI (m²/m²)",                  "mlai"
)

yearly_vars <- tribble(
  ~file,        ~var_label,                    ~ylab,                         ~stem,
  "agpp.out",   "Annual Gross Primary Production", "GPP (kgC/m²/year)",      "agpp",
  "anpp.out",   "Annual Net Primary Production",   "NPP (kgC/m²/year)",      "anpp",
  "cmass.out",  "Carbon Mass",                     "Carbon mass (kgC/m²)",  "cmass"
)

source("R_scripts/lpj_guess_hyd_functions.R")
source("R_scripts/run_lpj_guess_hyd_functions.R")

# ============================================================
# AVAILABLE VARIABLES TO PLOT
# ============================================================

# ---- DAILY variables (use with daily_plots$...) ----
# dgc          = Leaf conductance (gc)
# dgcwater     = Canopy conductance (gcwater)
# dpsisoil     = Soil water potential
# dpsileaf     = Leaf water potential
# dpsixylem    = Xylem water potential
# et_total     = Transpiration
# mort_cav_day = Hydraulic-failure mortality
# mort         = Total mortality
# mort_min     = Background mortality
# mort_greff   = Growth-efficiency mortality
# mort_cav     = Hydraulic mortality
# kappa_s_min  = kappa_s_min
# kappa_s_today= kappa_s_today

# ---- MONTHLY variables (use with monthly_plots$...) ----
# mevap = Evapotranspiration (mevap)
# mgpp  = GPP (mgpp)
# mlai  = LAI (mlai)
# mnpp  = NPP (mnpp)

# ---- YEARLY variables (use with yearly_plots$...) ----
# agpp  = Annual GPP (agpp)
# anpp  = Annual NPP (anpp)
# cmass = Carbon (cmass)

# ============================================================
# AVAILABLE SPECIES
# ============================================================
# ALL
# Beech
# Oak_pub
# Oak_rob
# Pine
# Spruce

# ============================================================
# AVAILABLE DAILY PLOT TYPES
# from DAILY source data
# ============================================================
# daily_raw_<species>
# monthly_<stat>_<species>
# yearly_<stat>_<species>
#
# where <stat> is:
# mean / min / max
#
# examples:
# daily_plots$dpsisoil$daily_raw_ALL
# daily_plots$dpsileaf$yearly_min_Beech
# daily_plots$et_total$monthly_max_Pine

# ============================================================
# AVAILABLE MONTHLY PLOT TYPES
# from MONTHLY source data
# ============================================================
# monthly_raw_<species>
# yearly_<stat>_<species>
#
# examples:
# monthly_plots$mlai$monthly_raw_ALL
# monthly_plots$mgpp$yearly_max_Oak_pub

# ============================================================
# AVAILABLE YEARLY PLOT TYPES
# from YEARLY source data
# ============================================================
# yearly_raw_<species>
#
# examples:
# yearly_plots$agpp$yearly_raw_ALL
# yearly_plots$cmass$yearly_raw_Beech

# ============================================================
# Show results
# ============================================================

# ---- daily source: raw daily + aggregated monthly/yearly ----
# daily_plots$dpsisoil$daily_raw_ALL
# daily_plots$dpsisoil$monthly_mean_ALL
# daily_plots$dpsisoil$yearly_mean_ALL
# daily_plots$dpsisoil$yearly_min_ALL
# 
# daily_plots$dpsixylem$daily_raw_ALL
# daily_plots$dpsixylem$monthly_min_ALL
# daily_plots$dpsixylem$yearly_min_ALL
# 
# daily_plots$dpsileaf$daily_raw_ALL
# daily_plots$dpsileaf$monthly_min_ALL
# daily_plots$dpsileaf$yearly_min_ALL
# 
# daily_plots$mort$yearly_mean_ALL
# daily_plots$mort_cav$yearly_mean_ALL
# daily_plots$mort_greff$yearly_mean_ALL
# daily_plots$mort_min$yearly_mean_ALL
# 
# daily_plots$kappa_s_min$daily_raw_ALL
# daily_plots$kappa_s_min$monthly_mean_ALL
# daily_plots$kappa_s_min$yearly_mean_ALL
# 
# daily_plots$kappa_s_today$daily_raw_ALL
# daily_plots$kappa_s_today$monthly_mean_ALL
# daily_plots$kappa_s_today$yearly_mean_ALL
# 
# # ---- monthly source: raw monthly + aggregated yearly ----
# monthly_plots$mlai$monthly_raw_ALL
# monthly_plots$mlai$yearly_mean_ALL
# 
# # ---- yearly source: raw yearly only ----
# yearly_plots$agpp$yearly_raw_ALL
# yearly_plots$cmass$yearly_raw_ALL
# 
# 
