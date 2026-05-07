vpd_func <- function(tmean, rh) {
  # saturation vapor pressure (hPa)
  es <- 6.11 * exp((2.5e6 / 461) * (1/273 - 1/(273 + tmean)))
  
  # vapor pressure deficit
  vpd <- ((100 - rh) / 100) * es
  
  return(vpd)
}