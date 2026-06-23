setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD/MeteoSwiss/MeteoSwiss_Visp")

library(dplyr)
library(tidyr)
library(readr)
library(ncdf4)
library(lubridate)

# ============================================================
# 1) Load data
# ============================================================
daily_filtered <- read.csv("Data/MeteoSwiss_station/all_VIS_selected_meteo.csv")

daily_filtered <- daily_filtered %>%
  mutate(
    date = as.Date(date),
    date = as.POSIXct(date, tz = "UTC")
  )

# ============================================================
# 2) Build station table (0-based index, primary station cols)
# ============================================================
stations <- daily_filtered %>%
  group_by(station_abbr) %>%
  summarise(
    lat = as.numeric(sprintf("%.6f", first(primary_station_lat))),
    lon = as.numeric(sprintf("%.6f", first(primary_station_lon))),
    name = first(primary_station_name),
    .groups = "drop"
  ) %>%
  arrange(station_abbr) %>%
  mutate(landid = row_number() - 1) # Starts at 0

# ============================================================
# 3) Join landid & Rename variables
# ============================================================
daily_filtered <- daily_filtered %>%
  left_join(
    stations %>% select(station_abbr, landid),
    by = "station_abbr"
  ) %>%
  rename(
    temperature       = `Air.temperature.2.m.above.ground..daily.mean`,
    global_radiation  = `Global.radiation..daily.mean`,
    relative_humidity = `Relative.air.humidity.2.m.above.ground..daily.mean`,
    wind_speed        = `Wind.speed.scalar..daily.mean.in.m.s`,
    precipitation     = `Precipitation..daily.total.0.UTC...0.UTC`
  ) %>%
  mutate(relative_humidity = relative_humidity / 100) # Scaling relative humidity

# ============================================================
# 4) Time axis
# ============================================================
time_origin <- as.POSIXct("1991-01-01 00:00:00", tz = "UTC")
all_dates <- sort(unique(daily_filtered$date))
time_num <- as.numeric(difftime(as.POSIXct(all_dates, tz = "UTC"), time_origin, units = "days"))

ntime <- length(all_dates)
nland <- nrow(stations)

# ============================================================
# 5) Metadata with CF Standard Names
# ============================================================
var_info <- data.frame(
  original_col = c("temperature", "global_radiation", "relative_humidity", "wind_speed", "precipitation"),
  standard_name = c("air_temperature", "surface_downwelling_shortwave_flux", "relative_humidity", "wind_speed", "precipitation_amount"),
  long_name = c("Air temperature daily mean", "Global radiation daily mean", "Relative humidity daily mean", "Wind speed daily mean", "Precipitation daily total"),
  units = c("K", "W m-2", "1", "m s-1", "kg m-2"),
  stringsAsFactors = FALSE
)

# ============================================================
# 6) Matrix builder
# ============================================================
make_matrix <- function(col_name) {
  full_grid <- tidyr::crossing(landid = stations$landid, date = all_dates)
  
  tmp <- full_grid %>%
    left_join(daily_filtered, by = c("landid", "date")) %>%
    select(landid, date, value = all_of(col_name)) %>%
    arrange(landid, date)
  
  mat <- matrix(tmp$value, nrow = nland, ncol = ntime, byrow = TRUE)
  storage.mode(mat) <- "numeric"
  return(mat)
}

# ============================================================
# 7) NetCDF writer
# ============================================================
write_nc <- function(mat, varname, std_name, long_name, units, outfile) {
  fillvalue <- -9999
  if (file.exists(outfile)) file.remove(outfile)
  
  dim_land <- ncdim_def("station", "index", stations$landid)
  dim_time <- ncdim_def("time", "days since 1991-01-01 00:00:00", as.integer(time_num), calendar = "proleptic_gregorian")
  
  var_lon <- ncvar_def("lon", "degrees_east", list(dim_land), fillvalue, "longitude", prec="double")
  var_lat <- ncvar_def("lat", "degrees_north", list(dim_land), fillvalue, "latitude", prec="double")
  var_data <- ncvar_def(varname, units, list(dim_land, dim_time), fillvalue, long_name, prec="float")
  
  nc <- nc_create(outfile, list(var_lon, var_lat, var_data))
  ncvar_put(nc, "lon", stations$lon)
  ncvar_put(nc, "lat", stations$lat)
  
  mat[is.na(mat)] <- fillvalue
  ncvar_put(nc, var_data, mat)
  
  ncatt_put(nc, "time", "calendar", "proleptic_gregorian")
  ncatt_put(nc, "lon", "standard_name", "longitude")
  ncatt_put(nc, "lat", "standard_name", "latitude")
  ncatt_put(nc, varname, "standard_name", std_name)
  ncatt_put(nc, varname, "coordinates", "lon lat")
  
  nc_close(nc)
}

# ============================================================
# 8) Processing Loop
# ============================================================
out_dir <- "Data/MeteoSwiss_station_to_netcdf_Visp"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(var_info))) {
  cat("\nProcessing variable:", var_info$original_col[i], "\n")
  mat <- make_matrix(var_info$original_col[i])
  outfile <- file.path(out_dir, paste0(var_info$original_col[i], ".nc"))
  
  write_nc(mat, var_info$standard_name[i], var_info$standard_name[i], var_info$long_name[i], var_info$units[i], outfile)
}

# ============================================================
# 9) Create gridlist files
# ============================================================
write.table(stations %>% select(landid, station_abbr), file.path(out_dir, "gridlist.txt"), sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(stations %>% filter(station_abbr == "VIS") %>% select(landid, station_abbr), file.path(out_dir, "gridlist_VIS.txt"), sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ============================================================
# 10) Create soil_SCCII.dat (0-based)
# ============================================================
soil_data <- stations %>% transmute(lon_str = sprintf("%.6f", lon), lat_str = sprintf("%.6f", lat), val = landid)
write.table(soil_data, file.path(out_dir, "soil_list.dat"), sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("\n✅ SUCCESS: NetCDF/Soil sync complete with 0-based indexing and humidity scaled.\n")
