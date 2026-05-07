setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

# 10 minutes
radiation_mean_10minute <- read.csv("MeteoSwiss_variable/global_radiation__ten_minutes_mean_10min.csv")

# hourly
temp_mean_hourly <- read.csv("MeteoSwiss_variable/air_temperature_2_m_above_ground__hourly_mean_hourly.csv")
relative_humidity_mean_hourly <- read.csv("MeteoSwiss_variable/relative_air_humidity_2_m_above_ground__hourly_mean_hourly.csv")
radiation_mean_hourly <- read.csv("MeteoSwiss_variable/global_radiation__hourly_mean_hourly.csv")
wind_speed_hourly <- read.csv("MeteoSwiss_variable/wind_speed_scalar__hourly_mean_in_m_s_hourly.csv")

# daily
temp_mean_daily <- read.csv("MeteoSwiss_variable/air_temperature_2_m_above_ground__daily_mean_daily.csv")
relative_humidity_mean_daily <- read.csv("MeteoSwiss_variable/relative_air_humidity_2_m_above_ground__hourly_mean_hourly.csv")
wind_speed_mean_daily <-read.csv("MeteoSwiss_variable/wind_speed_scalar__daily_mean_in_m_s_daily.csv")
global_radiation_mean_daily <- read.csv("MeteoSwiss_variable/global_radiation__daily_mean_daily.csv")
precip_total_daily <- read.csv("MeteoSwiss_variable/precipitation__daily_total_0_utc_-_0_utc_daily.csv")

