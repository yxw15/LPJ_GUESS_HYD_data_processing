# =========================================================
# 0. SETUP
# =========================================================

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)

climate_folder <- "SCCII/Climate_Hoelstein"


# =========================================================
# 1. LOAD METEOSWISS DATA
# =========================================================

meteo <- read.csv(
  "MeteoSwiss/MeteoSwiss_station_to_netcdf/all_filtered_19910101_to_20251231.csv"
) %>%
  mutate(date = as.Date(date_day)) %>%
  select(
    station_abbr,
    date,
    mean_temperature = Air.temperature.2.m.above.ground..daily.mean,
    radiation = Global.radiation..daily.mean,
    relative_humidity = Relative.air.humidity.2.m.above.ground..daily.mean,
    wind_speed = Wind.speed.scalar..daily.mean.in.m.s,
    precipitation = Precipitation..daily.total.0.UTC...0.UTC,
    lat = primary_station_lat,
    lon = primary_station_lon
  )


# =========================================================
# 2. LOAD + MERGE CLIMATE FILES
# =========================================================

files <- list.files(
  climate_folder,
  pattern = "Climate_.*_archive\\.txt$",
  full.names = TRUE
)

climate_all <- files %>%
  lapply(function(f) {
    read.table(
      f,
      header = TRUE,
      sep = "",
      quote = "\"",
      na.strings = "NA",
      stringsAsFactors = FALSE
    )
  }) %>%
  bind_rows() %>%
  mutate(
    timestamp = gsub("[()]", "", timestamp_UTC),
    timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S"),
    date = as.Date(timestamp)
  )


# =========================================================
# 3. AGGREGATE CLIMATE TO DAILY
# =========================================================

climate_daily <- climate_all %>%
  mutate(
    timestamp = as.POSIXct(timestamp, tz = "UTC"),
    date_only = as.Date(timestamp),
    time_num = hour(timestamp) + minute(timestamp) / 60
  ) %>%
  group_by(date_only) %>%
  summarise(
    
    # daytime means (07:30–19:30)
    mean_temperature = mean(Temp_degC_crane[time_num >= 7.5 & time_num <= 19.5], na.rm = TRUE) + 273.15,
    radiation = mean(Solar_Wm.2_crane[time_num >= 7.5 & time_num <= 19.5], na.rm = TRUE),
    relative_humidity = mean(Humid_percent_crane[time_num >= 7.5 & time_num <= 19.5], na.rm = TRUE),
    
    # full-day precipitation (sum!)
    precipitation = ifelse(
      all(is.na(dRain_mm_cranegap)),
      NA,
      sum(dRain_mm_cranegap, na.rm = TRUE)
    )
    
    # wind_speed intentionally omitted (kept commented if needed)
    # wind_speed = ifelse(all(is.na(Wind_ms.1_crane)), NA,
    #                     mean(Wind_ms.1_crane, na.rm = TRUE))
    
  ) %>%
  ungroup() %>%
  mutate(station_abbr = "RUE")

# =========================================================
# 4. SPLIT METEO DATA
# =========================================================

meteo_rue <- meteo %>% filter(station_abbr == "RUE")
meteo_other <- meteo %>% filter(station_abbr != "RUE")


# =========================================================
# 5. REPLACE RUE WITH CLIMATE
# =========================================================

meteo_rue_replaced <- meteo_rue %>%
  left_join(climate_daily, by = c("date", "station_abbr")) %>%
  mutate(
    mean_temperature = coalesce(mean_temperature.y, mean_temperature.x),
    radiation = coalesce(radiation.y, radiation.x),
    relative_humidity = coalesce(relative_humidity.y, relative_humidity.x),
    # wind_speed = coalesce(wind_speed.y, wind_speed.x),
    precipitation = coalesce(precipitation.y, precipitation.x)
  ) %>%
  select(
    station_abbr, date,
    mean_temperature,
    radiation,
    relative_humidity,
    wind_speed,
    precipitation,
    lat, lon
  )


# =========================================================
# 6. REBUILD FULL DATASET
# =========================================================

meteo_final <- bind_rows(meteo_other, meteo_rue_replaced) %>%
  arrange(station_abbr, date)


# =========================================================
# 7. QUICK QA CHECK
# =========================================================

meteo_final %>%
  filter(station_abbr == "RUE") %>%
  summarise(
    n = n(),
    temp_na = sum(is.na(mean_temperature)),
    rad_na = sum(is.na(radiation))
  )


# =========================================================
# 8. PLOT RUE: FINAL VS CLIMATE
# =========================================================

meteo_rue_final <- meteo_final %>%
  filter(station_abbr == "RUE")

comparison_final <- meteo_rue_final %>%
  inner_join(climate_daily, by = c("date", "station_abbr"),
             suffix = c("_meteo_final", "_climate"))

ggplot(comparison_final, aes(date)) +
  geom_line(aes(y = mean_temperature_meteo_final, color = "Meteo final")) +
  geom_line(aes(y = mean_temperature_climate, color = "Climate")) +
  theme_minimal() +
  labs(
    title = "RUE Temperature: Meteo Final vs Climate",
    y = "Temperature (K)",
    color = "Dataset"
  )


# =========================================================
# 9. FULL TIME SERIES (ALL VARIABLES)
# =========================================================

df_long <- meteo_final %>%
  filter(station_abbr == "RUE") %>%
  select(date,
         mean_temperature,
         radiation,
         relative_humidity,
         wind_speed,
         precipitation) %>%
  pivot_longer(-date, names_to = "variable", values_to = "value")


ggplot(df_long, aes(date, value)) +
  geom_line(alpha = 0.6, color = "dodgerblue") +
  facet_wrap(~variable, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(title = "Meteo Final Time Series (RUE)")

df_monthly <- meteo_final %>%
  filter(station_abbr == "RUE") %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(
    mean_temperature = mean(mean_temperature, na.rm = TRUE),
    radiation = mean(radiation, na.rm = TRUE),
    relative_humidity = mean(relative_humidity, na.rm = TRUE),
    wind_speed = mean(wind_speed, na.rm = TRUE),
    precipitation = mean(precipitation, na.rm = TRUE)
  ) %>%
  ungroup()

df_monthly_long <- df_monthly %>%
  pivot_longer(-month, names_to = "variable", values_to = "value")

df_yearly <- meteo_final %>%
  filter(station_abbr == "RUE") %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarise(
    mean_temperature = mean(mean_temperature, na.rm = TRUE),
    radiation = mean(radiation, na.rm = TRUE),
    relative_humidity = mean(relative_humidity, na.rm = TRUE),
    wind_speed = mean(wind_speed, na.rm = TRUE),
    precipitation = mean(precipitation, na.rm = TRUE)
  ) %>%
  ungroup()

df_yearly_long <- df_yearly %>%
  pivot_longer(-year, names_to = "variable", values_to = "value")


ggplot(df_monthly_long, aes(month, value)) +
  geom_line(color = "dodgerblue", alpha = 0.8) +
  facet_wrap(~variable, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(title = "Monthly Mean (RUE)")

ggplot(df_yearly_long, aes(year, value)) +
  geom_line(color = "dodgerblue", linewidth = 1) +
  # geom_point(color = "dodgerblue") +
  facet_wrap(~variable, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(title = "Yearly Mean (RUE)")

# =========================================================
# 10. SAVE FINAL DATASET
# =========================================================

write.csv(
  meteo_final,
  file.path(
    "MeteoSwiss/Hoelstein_station_to_netcdf/all_filled_hoelstein_19910101_to_20251231.csv"
  ),
  row.names = FALSE
)
