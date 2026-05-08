library(dplyr)
library(tidyr)
library(readr)
library(ncdf4)
library(lubridate)

# ============================================================
# 1) Core Processing Function
# ============================================================
process_meteo_to_netcdf <- function(input_csv, out_dir) {
  
  cat("\n============================================================\n")
  cat("Starting processing for:\n", input_csv, "\n")
  cat("============================================================\n")
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ============================================================
  # 2) Load data
  # ============================================================
  daily_filtered <- read.csv(input_csv)
  
  daily_filtered <- daily_filtered %>%
    mutate(
      date = as.Date(date),
      date = as.POSIXct(date, tz = "UTC")
    )
  
  # ============================================================
  # 3) Build station table (FORCED STRING PRECISION)
  # ============================================================
  stations <- daily_filtered %>%
    group_by(station_abbr) %>%
    summarise(
      lat = as.numeric(sprintf("%.6f", first(lat))),
      lon = as.numeric(sprintf("%.6f", first(lon))),
      .groups = "drop"
    ) %>%
    arrange(station_abbr) %>%
    mutate(landid = row_number())
  
  # ============================================================
  # 4) Join landid
  # ============================================================
  daily_filtered <- daily_filtered %>%
    left_join(
      stations %>% select(station_abbr, landid),
      by = "station_abbr"
    )
  
  # ============================================================
  # 5) Time axis
  # ============================================================
  time_origin <- as.POSIXct("1991-01-01 00:00:00", tz = "UTC")
  
  all_dates <- sort(unique(daily_filtered$date))
  
  time_num <- as.numeric(
    difftime(
      as.POSIXct(all_dates, tz = "UTC"),
      time_origin,
      units = "days"
    )
  )
  
  ntime <- length(all_dates)
  nland <- nrow(stations)
  
  # ============================================================
  # 6) Metadata with CF Standard Names
  # ============================================================
  var_info <- data.frame(
    original_col = c(
      "temperature",
      "global_radiation",
      "relative_humidity",
      "wind_speed",
      "precipitation"
    ),
    standard_name = c(
      "air_temperature",
      "surface_downwelling_shortwave_flux",
      "relative_humidity",
      "wind_speed",
      "precipitation_amount"
    ),
    long_name = c(
      "Air temperature daytime (7:30-19:30) mean",
      "Global radiation daytime (7:30-19:30) mean",
      "Relative humidity daytime (7:30-19:30) mean",
      "Wind speed daytime (7:30-19:30) mean",
      "Precipitation daily total"
    ),
    units = c(
      "K",
      "W m-2",
      "1",
      "m s-1",
      "kg m-2"
    ),
    stringsAsFactors = FALSE
  )
  
  # ============================================================
  # 7) Matrix builder
  # ============================================================
  make_matrix <- function(col_name) {
    
    full_grid <- tidyr::crossing(
      landid = stations$landid,
      date = all_dates
    )
    
    tmp <- full_grid %>%
      left_join(
        daily_filtered,
        by = c("landid", "date")
      ) %>%
      select(
        landid,
        date,
        value = all_of(col_name)
      ) %>%
      arrange(landid, date)
    
    mat <- matrix(
      tmp$value,
      nrow = nland,
      ncol = ntime,
      byrow = TRUE
    )
    
    storage.mode(mat) <- "numeric"
    
    return(mat)
  }
  
  # ============================================================
  # 8) NetCDF writer
  # ============================================================
  write_nc <- function(
    mat,
    varname,
    std_name,
    long_name,
    units,
    outfile
  ) {
    
    fillvalue <- -9999
    
    if (file.exists(outfile)) {
      file.remove(outfile)
    }
    
    dim_land <- ncdim_def(
      "station",
      "index",
      stations$landid
    )
    
    dim_time <- ncdim_def(
      "time",
      "days since 1991-01-01 00:00:00",
      as.integer(time_num),
      calendar = "proleptic_gregorian"
    )
    
    # Use DOUBLE precision for coordinates
    var_lon <- ncvar_def(
      "lon",
      "degrees_east",
      list(dim_land),
      fillvalue,
      "longitude",
      prec = "double"
    )
    
    var_lat <- ncvar_def(
      "lat",
      "degrees_north",
      list(dim_land),
      fillvalue,
      "latitude",
      prec = "double"
    )
    
    # Climate data stored as float
    var_data <- ncvar_def(
      varname,
      units,
      list(dim_land, dim_time),
      fillvalue,
      long_name,
      prec = "float"
    )
    
    nc <- nc_create(
      outfile,
      list(var_lon, var_lat, var_data)
    )
    
    # Write coordinates
    ncvar_put(nc, "lon", stations$lon)
    ncvar_put(nc, "lat", stations$lat)
    
    # Replace NA
    mat[is.na(mat)] <- fillvalue
    
    # Write climate data
    ncvar_put(nc, var_data, mat)
    
    # ============================================================
    # CF Attributes
    # ============================================================
    ncatt_put(nc, "time", "calendar", "proleptic_gregorian")
    
    ncatt_put(nc, "lon", "standard_name", "longitude")
    ncatt_put(nc, "lat", "standard_name", "latitude")
    
    ncatt_put(nc, varname, "standard_name", std_name)
    ncatt_put(nc, varname, "coordinates", "lon lat")
    
    nc_close(nc)
  }
  
  # ============================================================
  # 9) Processing Loop
  # ============================================================
  for (i in seq_len(nrow(var_info))) {
    
    cat("\nProcessing variable:",
        var_info$original_col[i],
        "\n")
    
    mat <- make_matrix(var_info$original_col[i])
    
    outfile <- file.path(
      out_dir,
      paste0(var_info$original_col[i], ".nc")
    )
    
    write_nc(
      mat = mat,
      varname = var_info$standard_name[i],
      std_name = var_info$standard_name[i],
      long_name = var_info$long_name[i],
      units = var_info$units[i],
      outfile = outfile
    )
  }
  
  # ============================================================
  # 10) Create gridlist_SCCII.txt
  # ============================================================
  write.table(
    stations %>% select(landid, station_abbr),
    file.path(out_dir, "gridlist_SCCII.txt"),
    sep = " ",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  
  # ============================================================
  # 11) Create soil_SCCII.dat
  # ============================================================
  soil_data <- stations %>%
    transmute(
      lon_str = sprintf("%.6f", lon),
      lat_str = sprintf("%.6f", lat),
      val = 6
    )
  
  write.table(
    soil_data,
    file.path(out_dir, "soil_SCCII.dat"),
    sep = " ",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  
  cat("\n✅ SUCCESS:",
      out_dir,
      "created with perfectly synchronized coordinates at 6 decimal places.\n")
}

# ============================================================
# 12) Run BOTH datasets
# ============================================================

# ------------------------------------------------------------
# A) NORMAL DATASET
# ------------------------------------------------------------
process_meteo_to_netcdf(
  input_csv =
    "Data/MeteoSwiss_station/all_stations_RUE_replaced_daytime.csv",
  
  out_dir =
    "Data/MeteoSwiss_station_to_netcdf_daytime"
)

# ------------------------------------------------------------
# B) DROUGHT DATASET
# ------------------------------------------------------------
process_meteo_to_netcdf(
  input_csv =
    "Data/MeteoSwiss_station/all_stations_RUE_replaced_daytime_drought.csv",
  
  out_dir =
    "Data/MeteoSwiss_station_to_netcdf_daytime_drought"
)

cat("\n============================================================\n")
cat("ALL PROCESSING FINISHED SUCCESSFULLY\n")
cat("Generated BOTH:\n")
cat("  1. Standard climate NetCDF dataset\n")
cat("  2. Drought climate NetCDF dataset\n")
cat("============================================================\n")