setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD/MeteoSwiss")
# setwd("~/Documents/Manuscript3")

# ============================================================
# 0) Setup
# ============================================================

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

# ============================================================
# 10) Drought experiments
# ============================================================
# Define your specific start/end dates for each year
drought_periods <- tribble(
  ~year, ~start_date,  ~end_date,
  2023,  "2023-04-05", "2023-10-16",
  2024,  "2024-04-04", "2024-10-24",
  2025,  "2025-03-31", "2025-10-28" 
) %>%
  mutate(
    start_date = as.Date(start_date),
    end_date = as.Date(end_date)
  )

# Ensure your date column is in Date format
meteo_final <- meteo_final %>%
  mutate(date = as.Date(date),
         year_val = year(date))

# Join with our period definitions and apply logic
meteo_drought <- meteo_final %>%
  left_join(drought_periods, by = c("year_val" = "year")) %>%
  mutate(
    precipitation = if_else(
      !is.na(start_date) & date >= start_date & date <= end_date,
      precipitation * 0.5, # Reduce by 50%
      precipitation        # Keep original otherwise
    )
  ) %>%
  # Clean up the helper columns used for the calculation
  select(-start_date, -end_date, -year_val)

head(meteo_drought %>% filter(date == "2023-05-01"))

drought_out_file <- "MeteoSwiss_station/all_stations_RUE_replaced_daytime_drought.csv"

write.csv(meteo_drought, drought_out_file, row.names = FALSE)

cat("Success! Drought-adjusted dataset saved to:", drought_out_file, "\n")