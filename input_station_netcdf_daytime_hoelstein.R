setwd("~/Documents/Manuscript3")

source("MeteoSwiss_data_processing/download_nearest_station_ten.R")
source("MeteoSwiss_data_processing/download_nearest_station_one_minute_hour_daily.R")
source("LPJ_GUESS_HYD_data_processing/check_MeteoSwiss_Hoelstein_data.R")
source("LPJ_GUESS_HYD_data_processing/input_replace_daytime_hoelstein.R")
source("LPJ_GUESS_HYD_data_processing/input_station_to_netcdf_treatment.R")

