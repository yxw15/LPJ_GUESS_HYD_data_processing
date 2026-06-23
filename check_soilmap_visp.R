setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD/Switzerland_soil_maps")
library(terra)

sand_folder <- "sand_maps"
clay_folder <- "clay_maps"
organic_folder <- "soc_maps"

sand_files <- list.files(sand_folder, pattern = "\\.tif$", full.names = TRUE)
clay_files <- list.files(clay_folder, pattern = "\\.tif$", full.names = TRUE)
organic_files <- list.files(organic_folder, pattern = "\\.tif$", full.names = TRUE)

sand <- rast(sand_files)
clay <- rast(clay_files)
organic <- rast(organic_files)

sand.mean <- app(sand, mean)
clay.mean <- app(clay, mean)
organic.mean <- app(organic, mean)

##### SCCII site #####
# target_lat <- 47.439
# target_lon <- 7.776

##### VIS site #####
target_lat <- 46.3029
target_lon <- 7.842958

# Create point (WGS84)
p <- vect(data.frame(x = target_lon, y = target_lat),
          geom = c("x", "y"),
          crs = "EPSG:4326")

# Stack rasters
r_all <- c(sand.mean, clay.mean, organic.mean)
names(r_all) <- c("sand", "clay", "organic")

# Project point to raster CRS
p_proj <- project(p, crs(r_all))

# Extract values (nearest cell)
result <- extract(r_all, p_proj, cells = TRUE, xy = TRUE)

result
