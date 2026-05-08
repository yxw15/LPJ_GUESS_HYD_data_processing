# ============================================================
# 0) Setup
# ============================================================

library(dplyr)
library(lubridate)
library(tidyr)
library(readr)
library(ggplot2)

target_station <- "RUE"

# ============================================================
# 1) Load MeteoSwiss data
# ============================================================

meteo_daily_wholeday <- read.csv(
  "Data/MeteoSwiss_station/all_selected_19910101_to_20251231.csv"
)

hourly_wind_speed_mean <- read.csv(
  "Data/MeteoSwiss_variable/wind_speed_scalar__hourly_mean_in_m_s_hourly.csv"
)

hourly_global_radiation_mean <- read.csv(
  "Data/MeteoSwiss_variable/global_radiation__hourly_mean_hourly.csv"
)

hourly_relative_humidity_mean <- read.csv(
  "Data/MeteoSwiss_variable/relative_air_humidity_2_m_above_ground__hourly_mean_hourly.csv"
)

hourly_temperature_mean <- read.csv(
  "Data/MeteoSwiss_variable/air_temperature_2_m_above_ground__hourly_mean_hourly.csv"
)

daily_precipitation_total <- read.csv(
  "Data/MeteoSwiss_variable/precipitation__daily_total_0_utc_-_0_utc_daily.csv"
)

# ============================================================
# 2) Preprocess hourly MeteoSwiss variables
# ============================================================

hourly_relative_humidity_mean <- hourly_relative_humidity_mean %>%
  mutate(value = value / 100)

wind <- hourly_wind_speed_mean %>%
  rename(wind_speed = value)

rad <- hourly_global_radiation_mean %>%
  rename(global_radiation = value)

rh <- hourly_relative_humidity_mean %>%
  rename(relative_humidity = value)

temp <- hourly_temperature_mean %>%
  rename(temperature = value) %>%
  mutate(temperature = temperature + 273.15)

# ============================================================
# 3) Merge hourly data (FILTER RUE EARLY)
# ============================================================

hourly_all <- wind %>%
  full_join(rad, by = c("datetime", "station_abbr")) %>%
  full_join(rh, by = c("datetime", "station_abbr")) %>%
  full_join(temp, by = c("datetime", "station_abbr")) %>%
  filter(station_abbr == target_station) %>%   # ✅ IMPORTANT
  mutate(
    datetime = ymd_hms(datetime, tz = "UTC"),
    date = as.Date(datetime),
    time = format(datetime, "%H:%M")
  ) %>%
  filter(time >= "07:30" & time <= "19:30")

# ============================================================
# 4) MeteoSwiss DAYTIME daily aggregation
# ============================================================

meteo_daily_daytime <- hourly_all %>%
  group_by(station_abbr, date) %>%
  summarise(
    wind_speed = mean(wind_speed, na.rm = TRUE),
    global_radiation = mean(global_radiation, na.rm = TRUE),
    relative_humidity = mean(relative_humidity, na.rm = TRUE),
    temperature = mean(temperature, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 5) MeteoSwiss WHOLE DAY dataset (FILTER RUE)
# ============================================================

meteo_daily_wholeday <- meteo_daily_wholeday %>% 
  filter(station_abbr == target_station) %>%
  rename(
    wind_speed = wind_speed_daily_mean,
    global_radiation = global_radiation_daily_mean,
    precipitation = precipitation_daily_total,
    temperature = temperature_daily_mean,
    relative_humidity = relative_humidity_daily_mean
  ) %>%
  mutate(date = as.Date(date))

# ============================================================
# 6) Precipitation (FILTER RUE BEFORE JOIN)
# ============================================================

precipitation <- daily_precipitation_total %>%
  mutate(date = as.Date(ymd_hms(datetime, tz = "UTC"))) %>%
  filter(station_abbr == target_station) %>%   # ✅ IMPORTANT
  select(station_abbr, date, precipitation = value)

meteo_daily_daytime <- meteo_daily_daytime %>%
  left_join(precipitation, by = c("station_abbr", "date"))

# ============================================================
# 7) Hoelstein DAYTIME dataset
# ============================================================

hoelstein_folder <- "Data/SCCII/Climate_Hoelstein"

files <- list.files(
  hoelstein_folder,
  pattern = "Climate_.*_archive\\.txt$",
  full.names = TRUE
)

climate_all <- files %>%
  lapply(function(f) {
    read.table(f,
               header = TRUE,
               sep = "",
               quote = "\"",
               na.strings = "NA")
  }) %>%
  bind_rows() %>%
  mutate(
    timestamp = gsub("[()]", "", timestamp_UTC),
    timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S"),
    date = as.Date(timestamp),
    time_num = hour(timestamp) + minute(timestamp) / 60
  )

hoelstein_daily_daytime <- climate_all %>%
  group_by(date) %>%
  summarise(
    wind_speed = NA_real_,
    global_radiation = mean(Solar_Wm.2_crane[time_num >= 7.5 & time_num <= 19.5], na.rm = TRUE),
    relative_humidity = mean(Humid_percent_crane[time_num >= 7.5 & time_num <= 19.5], na.rm = TRUE) / 100,
    temperature = mean(Temp_degC_crane[time_num >= 7.5 & time_num <= 19.5], na.rm = TRUE) + 273.15,
    precipitation = sum(dRain_mm_cranegap, na.rm = TRUE),
    station_abbr = target_station
  ) %>%
  filter(date > as.Date("2016-12-31"))

# ============================================================
# 8) Combine datasets
# ============================================================

df_all <- bind_rows(
  meteo_daily_daytime %>% mutate(source = "MeteoSwiss_daytime"),
  meteo_daily_wholeday %>% mutate(source = "MeteoSwiss_whole_day"),
  hoelstein_daily_daytime %>% mutate(source = "Hoelstein_daytime")
) %>%
  select(-any_of(c("station", "lat", "lon", "date_day")))

df_long <- df_all %>%
  pivot_longer(
    cols = c(wind_speed,
             global_radiation,
             relative_humidity,
             temperature,
             precipitation),
    names_to = "variable",
    values_to = "value"
  )

# ============================================================
# 9) Plot
# ============================================================

ggplot(df_long, aes(x = date, y = value, color = source)) +
  geom_line(alpha = 0.6) +
  facet_wrap(~variable, scales = "free_y", ncol = 1) +
  scale_color_manual(
    values = c(
      "MeteoSwiss_daytime" = "orange",
      "MeteoSwiss_whole_day" = "dodgerblue",
      "Hoelstein_daytime" = "red"
    )
  ) +
  labs(
    title = paste("Meteorological comparison at station", target_station),
    x = "Date",
    y = "Value",
    color = "Dataset"
  ) +
  theme_minimal()

write.csv(
  df_all,
  "Data/MeteoSwiss_station/RUE_meteo_hoelstein.csv",
  row.names = FALSE
)

