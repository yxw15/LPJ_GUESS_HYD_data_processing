library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(plantecophys)
library(stringr)

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

# ============================================================
# BASE THEME
# ============================================================

base_theme <- theme_minimal() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(color = "black", size = 14),
    legend.position = "bottom",
    plot.title  = element_text(hjust = 0.5, size = 18),
    axis.title  = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    panel.border = element_blank(),
    strip.text = element_text(size = 14)
  )

# ============================================================
# LOAD DATA
# ============================================================

sap_flux_gc <- read.csv("SCCII/sap_flux_gc.csv")

meteo_temp   <- read.csv("MeteoSwiss_variable/air_temperature_2_m_above_ground__hourly_mean_hourly.csv")
meteo_precip <- read.csv("MeteoSwiss_variable/precipitation__daily_total_0_utc_-_0_utc_daily.csv")
meteo_rad    <- read.csv("MeteoSwiss_variable/global_radiation__ten_minutes_mean_10min.csv")
meteo_rh     <- read.csv("MeteoSwiss_variable/relative_air_humidity_2_m_above_ground__hourly_mean_hourly.csv")
meteo_wind   <- read.csv("MeteoSwiss_variable/wind_speed_scalar__hourly_mean_in_m_s_hourly.csv")

# ============================================================
# FILTER STATION
# ============================================================

meteo_temp   <- meteo_temp   %>% filter(station_abbr == "RUE")
meteo_precip <- meteo_precip %>% filter(station_abbr == "RUE")
meteo_rad    <- meteo_rad    %>% filter(station_abbr == "RUE")
meteo_rh     <- meteo_rh     %>% filter(station_abbr == "RUE")
meteo_wind   <- meteo_wind   %>% filter(station_abbr == "RUE")

# ============================================================
# PARSE METEO TIME
# ============================================================

meteo_temp <- meteo_temp %>%
  mutate(time_hour = floor_date(ymd_hms(datetime, tz="UTC"), "hour"))

meteo_rh <- meteo_rh %>%
  mutate(time_hour = floor_date(ymd_hms(datetime, tz="UTC"), "hour"))

meteo_wind <- meteo_wind %>%
  mutate(time_hour = floor_date(ymd_hms(datetime, tz="UTC"), "hour"))

meteo_rad <- meteo_rad %>%
  mutate(time_10min = floor_date(ymd_hms(datetime, tz="UTC"), "10 minutes"))

meteo_precip <- meteo_precip %>%
  mutate(date_day = as.Date(ymd_hms(datetime, tz="UTC")))

# ============================================================
# PREPARE SAP FLOW DATA
# ============================================================

sap_flux_gc <- sap_flux_gc %>%
  mutate(
    time_hour  = floor_date(ymd_hms(timestamp, tz="UTC"), "hour"),
    time_10min = floor_date(ymd_hms(timestamp, tz="UTC"), "10 minutes"),
    date_day   = as.Date(ymd_hms(timestamp, tz="UTC"))
  )

# ============================================================
# MERGE METEO DATA → FINAL DATASET NAME
# ============================================================

sap_flux_gc_hoelstein_meteo <- sap_flux_gc %>%
  # hourly meteo
  left_join(
    meteo_temp %>% select(time_hour, temp = value),
    by = "time_hour"
  ) %>%
  left_join(
    meteo_rh %>% select(time_hour, relhum = value),
    by = "time_hour"
  ) %>%
  left_join(
    meteo_wind %>% select(time_hour, wind = value),
    by = "time_hour"
  ) %>%
  
  # 10-min radiation
  left_join(
    meteo_rad %>% select(time_10min, radiation = value),
    by = "time_10min"
  ) %>%
  
  # daily precipitation
  left_join(
    meteo_precip %>% select(date_day, precip = value),
    by = "date_day"
  ) %>%
  
  # VPD
  mutate(
    vpd_meteo = RHtoVPD(relhum, temp)
  )

# ============================================================
# CHECK
# ============================================================

head(sap_flux_gc_hoelstein_meteo)
write.csv(sap_flux_gc_hoelstein_meteo, file = "SCCII/sap_flux_gc_hoelstein_meteo.csv", row.names = FALSE)

