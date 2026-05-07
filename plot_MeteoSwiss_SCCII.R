setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

# ============================================================
# 1. VPD CALCULATION & NETCDF CREATION
# ============================================================
suppressPackageStartupMessages({
  library(ncdf4)
  library(dplyr)
  library(lubridate)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

data_folder <- "MeteoSwiss/MeteoSwiss_station_to_netcdf"
output_dir  <- "Figures/climate"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------ THEME ------------------
base_theme <- theme_minimal() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(color = "black", size = 14),
    legend.position = "bottom",
    plot.title  = element_text(hjust = 0.5, size = 18, color = "black"),
    axis.title  = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    strip.text = element_text(size = 14, face = "bold")
  )

# ------------------ LOAD DATA ------------------
temp_file <- file.path(data_folder, "mean_temperature.nc")
rh_file   <- file.path(data_folder, "relative_humidity.nc")

nc_temp <- nc_open(temp_file)
nc_rh   <- nc_open(rh_file)

temp_k   <- ncvar_get(nc_temp, "mean_temperature")
rh_frac  <- ncvar_get(nc_rh, "relative_humidity")
stations <- ncvar_get(nc_temp, "station")
time_vals <- ncvar_get(nc_temp, "time")
lons <- ncvar_get(nc_temp, "lon")
lats <- ncvar_get(nc_temp, "lat")

nc_close(nc_temp)
nc_close(nc_rh)

# ------------------ CLEAN + CONVERT ------------------
temp_k[temp_k == -9999] <- NA
rh_frac[rh_frac == -9999] <- NA

temp_c  <- temp_k - 273.15
rh_perc <- rh_frac * 100

# ------------------ VPD FUNCTION ------------------
calculate_vpd <- function(T, RH) {
  es <- 0.6108 * exp(17.27 * T / (T + 237.3))
  ea <- es * (RH / 100)
  es - ea
}

# ------------------ FAST VPD (VECTORIZED) ------------------
cat("Calculating VPD (vectorized)...\n")

vpd_matrix <- calculate_vpd(temp_c, rh_perc)

vpd_matrix[!is.finite(vpd_matrix)] <- NA
vpd_matrix[is.na(vpd_matrix)] <- -9999

# ------------------ WRITE NETCDF ------------------
dim_station <- ncdim_def("station", "index", 1:length(stations))
dim_time    <- ncdim_def("time", "days since 1991-01-01", as.integer(time_vals))

var_lon <- ncvar_def("lon", "degrees_east", list(dim_station))
var_lat <- ncvar_def("lat", "degrees_north", list(dim_station))
var_vpd <- ncvar_def("vpd", "kPa", list(dim_station, dim_time), -9999,
                     "Vapor Pressure Deficit")

vpd_out_path <- file.path(data_folder, "vpd.nc")

nc_out <- nc_create(vpd_out_path, list(var_lon, var_lat, var_vpd))

ncvar_put(nc_out, var_lon, lons)
ncvar_put(nc_out, var_lat, lats)
ncvar_put(nc_out, var_vpd, vpd_matrix)

ncatt_put(nc_out, "vpd", "units", "kPa")

nc_close(nc_out)

cat("✅ VPD NetCDF created\n")

# ============================================================
# 2. EXTRACTION & PLOTTING
# ============================================================

target_station_idx <- 6

extract_station_data <- function(file_path) {
  nc <- nc_open(file_path)
  
  var_name <- setdiff(names(nc$var), c("lon","lat","station"))[1]
  
  unit <- ncatt_get(nc, var_name, "units")$value
  if (is.null(unit)) unit <- "unknown"
  
  time_vals <- ncvar_get(nc, "time")
  origin <- as.Date(strsplit(ncatt_get(nc,"time","units")$value," ")[[1]][3])
  dates <- origin + time_vals
  
  values <- ncvar_get(nc, var_name)[target_station_idx, ]
  
  fill <- ncatt_get(nc, var_name, "_FillValue")$value
  if (!is.null(fill)) values[values == fill] <- NA
  
  nc_close(nc)
  
  tibble(date = dates, variable = var_name, unit = unit, value = values)
}

files <- list.files(data_folder, "\\.nc$", full.names = TRUE)

all_data_long <- map_dfr(files, extract_station_data) %>%
  filter(year(date) <= 2025)

# ------------------ AGGREGATION ------------------
monthly_df <- all_data_long %>%
  group_by(month = floor_date(date, "month"), variable, unit) %>%
  summarise(
    value = if (first(variable) == "precipitation") sum(value, na.rm=TRUE)
    else mean(value, na.rm=TRUE),
    .groups="drop"
  )

yearly_df <- all_data_long %>%
  group_by(year = floor_date(date, "year"), variable, unit) %>%
  summarise(
    value = if (first(variable) == "precipitation") sum(value, na.rm=TRUE)
    else mean(value, na.rm=TRUE),
    .groups="drop"
  )

# ------------------ PLOT FUNCTION ------------------
plot_all <- function(df, time_col, title_suffix) {
  
  df <- df %>%
    filter(!is.na(value)) %>%
    mutate(label = paste0(gsub("_"," ",variable)," (",unit,")"))
  
  year_min <- year(min(df[[time_col]]))
  year_max <- year(max(df[[time_col]]))
  
  breaks <- seq(floor(year_min/5)*5, ceiling(year_max/5)*5, 5)
  
  ggplot(df, aes(x=.data[[time_col]], y=value)) +
    geom_line(color="steelblue") +
    facet_wrap(~label, scales="free_y", ncol=2) +
    scale_x_date(
      breaks = as.Date(paste0(breaks,"-01-01")),
      labels = breaks
    ) +
    base_theme +
    labs(
      title = paste("MeteoSwiss Station", target_station_idx, "-", title_suffix),
      x = "Year",
      y = ""
    )
}

# ------------------ PLOT + SAVE ------------------
if (nrow(monthly_df) > 0) {
  p1 <- plot_all(monthly_df, "month", "Monthly (1991–2025)")
  print(p1)
  ggsave(file.path(output_dir, paste0("station_",target_station_idx,"_monthly.png")), p1, width=12, height=7, dpi=300)
}

if (nrow(yearly_df) > 0) {
  p2 <- plot_all(yearly_df, "year", "Yearly (1991–2025)")
  print(p2)
  ggsave(file.path(output_dir, paste0("station_",target_station_idx,"_yearly.png")), p2, width=12, height=7, dpi=300)
}

cat("✅ SUCCESS: Everything completed\n")

# ------------------ SUMMARY ------------------
all_data_long %>%
  group_by(variable) %>%
  summarise(
    n=sum(!is.na(value)),
    mean=mean(value,na.rm=TRUE),
    sd=sd(value,na.rm=TRUE),
    min=min(value,na.rm=TRUE),
    max=max(value,na.rm=TRUE)
  ) %>%
  print()