# ============================================================
# 0) Setup
# ============================================================

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD/MeteoSwiss")

library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)

# ============================================================
# 1) Load data
# ============================================================

meteo_daily_wholeday <- read.csv(
  "MeteoSwiss_station/all_selected_19910101_to_20251231.csv"
)

data_daytime <- read.csv(
  "MeteoSwiss_station/RUE_meteo_hoelstein.csv"
)

# ============================================================
# 2) Clean WHOLE-DAY MeteoSwiss data
# ============================================================

meteo_clean <- meteo_daily_wholeday %>%
  select(-date_day) %>%
  rename(
    wind_speed = wind_speed_daily_mean,
    global_radiation = global_radiation_daily_mean,
    precipitation = precipitation_daily_total,
    temperature = temperature_daily_mean,
    relative_humidity = relative_humidity_daily_mean
  ) %>%
  mutate(date = as.Date(date))

# ============================================================
# 3) Split stations
# ============================================================

meteo_RUE <- meteo_clean %>%
  filter(station_abbr == "RUE") %>%
  mutate(date = as.Date(date))

meteo_other <- meteo_clean %>%
  filter(station_abbr != "RUE") %>%
  mutate(date = as.Date(date))

# ============================================================
# 4) Prepare replacement datasets (FIXED DATE TYPE)
# ============================================================

meteo_daytime <- data_daytime %>%
  filter(source == "MeteoSwiss_daytime") %>%
  select(-source) %>%
  mutate(date = as.Date(date))

hoelstein_daytime <- data_daytime %>%
  filter(source == "Hoelstein_daytime") %>%
  select(-source) %>%
  mutate(date = as.Date(date))

# ============================================================
# 5) Replace ONLY RUE values (Hoelstein > MeteoSwiss daytime)
# ============================================================

meteo_RUE_updated <- meteo_RUE %>%
  left_join(meteo_daytime,
            by = c("station_abbr", "date"),
            suffix = c("", "_day")) %>%
  left_join(hoelstein_daytime,
            by = c("station_abbr", "date"),
            suffix = c("", "_hoel")) %>%
  mutate(
    wind_speed = coalesce(wind_speed_hoel, wind_speed_day, wind_speed),
    global_radiation = coalesce(global_radiation_hoel, global_radiation_day, global_radiation),
    relative_humidity = coalesce(relative_humidity_hoel, relative_humidity_day, relative_humidity),
    temperature = coalesce(temperature_hoel, temperature_day, temperature),
    precipitation = coalesce(precipitation_hoel, precipitation_day, precipitation)
  ) %>%
  select(
    station_abbr, station, lat, lon, date,
    wind_speed, global_radiation,
    relative_humidity, temperature, precipitation
  )

# ============================================================
# 6) Recombine ALL stations
# ============================================================

meteo_final <- bind_rows(
  meteo_RUE_updated,
  meteo_other
)

# ============================================================
# 7) Station check
# ============================================================

station_counts <- meteo_final %>%
  count(station_abbr)

print(station_counts)

# ============================================================
# 8) RUE time series plot
# ============================================================

rue_plot_data <- meteo_final %>%
  filter(station_abbr == "RUE") %>%
  pivot_longer(
    cols = c(
      wind_speed,
      global_radiation,
      relative_humidity,
      temperature,
      precipitation
    ),
    names_to = "variable",
    values_to = "value"
  )

ggplot(rue_plot_data, aes(x = date, y = value)) +
  geom_line(color = "dodgerblue", alpha = 0.7) +
  facet_wrap(~variable, scales = "free_y", ncol = 1) +
  labs(
    title = "RUE Station (Updated: Hoelstein + MeteoSwiss daytime)",
    x = "Date",
    y = "Value"
  ) +
  theme_minimal()

# ============================================================
# 9) Save final dataset
# ============================================================

out_file <- "MeteoSwiss_station/all_stations_RUE_replaced_daytime.csv"

write.csv(meteo_final, out_file, row.names = FALSE)

cat("Saved final dataset to:", out_file, "\n")