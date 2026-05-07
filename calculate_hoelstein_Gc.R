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
    plot.title  = element_text(hjust = 0.5, size = 18, color = "black"),
    axis.title  = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    panel.border = element_blank(),
    strip.text = element_text(size = 14)
  )

# ============================================================
# SPECIES ORDER + COLORS
# ============================================================

species_levels <- c("Oak", "Beech", "Spruce", "Pine")

species_palette <- c(
  Oak="#E69F00",
  Beech="#0072B2",
  Spruce="#009E73",
  Pine="#F0E442"
)

# ============================================================
# LOAD DATA
# ============================================================

sap_flux_raw <- readRDS("SCCII/sap_flux_density/sap_flux_density_yixuan.rds")

climate_folder <- "SCCII/Climate_Hoelstein"

files <- list.files(
  climate_folder,
  pattern = "Climate_.*_archive\\.txt$",
  full.names = TRUE
)

climate_all <- files %>%
  lapply(function(f) {
    read.table(f, header = TRUE, sep = "", quote = "\"",
               na.strings = "NA", stringsAsFactors = FALSE)
  }) %>%
  bind_rows() %>%
  mutate(
    timestamp = gsub("[()]", "", timestamp_UTC),
    timestamp = trimws(timestamp),
    timestamp = ifelse(timestamp %in% c("", "NA"), NA, timestamp),
    timestamp = parse_date_time(timestamp, orders = c("ymd HMS", "ymd HM")),
    timestamp = floor_date(timestamp, "10 minutes")
  ) %>%
  filter(!is.na(timestamp)) %>%
  select(timestamp, Humid_percent_crane, Temp_degC_crane, Solar_Wm.2_crane, dRain_mm_cranegap)

# ============================================================
# CLEAN SAP DATA
# ============================================================

sap_flux <- sap_flux_raw %>%
  mutate(
    timestamp_clean = str_squish(as.character(timestamp)),
    timestamp = parse_date_time(timestamp_clean, orders = c("dmy HM", "dmy HMS")),
    timestamp = floor_date(timestamp, "10 minutes")
  )

# ============================================================
# MERGE + STANDARDISE
# ============================================================

sap_flux_merged <- sap_flux %>%
  left_join(climate_all, by = "timestamp") %>%
  mutate(
    species = recode(species,
                     "Fs" = "Beech",
                     "Pa" = "Spruce",
                     "Ps" = "Pine",
                     "Qs" = "Oak"),
    species = factor(species, levels = species_levels)
  ) %>%
  filter(!is.na(sfd)) %>%
  mutate(
    vpd = RHtoVPD(Humid_percent_crane, Temp_degC_crane),
    year = year(timestamp),
    date = as.Date(timestamp)
  )

# ============================================================
# CANOPY CONDUCTANCE
# ============================================================

eta <- 44.6
T0 <- 273
h <- 500

sap_flux_gc <- sap_flux_merged %>%
  mutate(
    sfd_kg = sfd / 3600,
    G_asw = ((115.8 + 0.4236 * Temp_degC_crane) * sfd_kg / vpd) *
      eta * (T0 / (T0 + Temp_degC_crane)) *
      exp(-0.00012 * h),
    G_ms = G_asw / eta
  )

# ============================================================
# Q90 FUNCTION
# ============================================================

q90_daily <- function(df, value_col) {
  df %>%
    group_by(species, date) %>%
    filter(.data[[value_col]] >= quantile(.data[[value_col]], 0.9, na.rm = TRUE)) %>%
    summarise(value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")
}

# ============================================================
# SFD Q90
# ============================================================

sfd_q90 <- q90_daily(sap_flux_gc, "sfd")

p_sfd_q90 <- ggplot(sfd_q90, aes(date, value, color = species)) +
  geom_point(size = 0.6) +
  scale_color_manual(values = species_palette) +
  facet_wrap(~species, ncol = 2) +
  base_theme +
  labs(title = "SFD Q90", y = "SFD")

ggsave(
  filename = "Figures/Hoelstein/p_sfd_q90.png",
  plot = p_sfd_q90,
  width = 10,
  height = 6,
  dpi = 300
)

# ============================================================
# G_asw Q90
# ============================================================

g_asw_q90 <- q90_daily(sap_flux_gc, "G_asw")

p_gasw_q90 <- ggplot(g_asw_q90, aes(date, value, color = species)) +
  geom_point(size = 0.6) +
  scale_color_manual(values = species_palette) +
  facet_wrap(~species, ncol = 2) +
  base_theme +
  labs(title = "G_asw Q90", y = "G_asw")

ggsave(
  filename = "Figures/Hoelstein/p_Gasw_q90.png",
  plot = p_gasw_q90,
  width = 10,
  height = 6,
  dpi = 300
)

# ============================================================
# G_ms Q90
# ============================================================

g_ms_q90 <- q90_daily(sap_flux_gc, "G_ms")

p_gms_q90 <- ggplot(g_ms_q90, aes(date, value, color = species)) +
  geom_point(size = 0.6) +
  scale_color_manual(values = species_palette) +
  facet_wrap(~species, ncol = 2) +
  base_theme +
  labs(title = "G_ms Q90", y = "G_ms")

ggsave(
  filename = "Figures/Hoelstein/p_Gms_q90.png",
  plot = p_gms_q90,
  width = 10,
  height = 6,
  dpi = 300
)


# ============================================================
# YEARLY PLOTS FUNCTION
# ============================================================

plot_yearly <- function(df, value_col, name) {
  
  years <- unique(df$year)
  
  for (y in years) {
    
    df_y <- df %>%
      filter(year == y)
    
    p <- ggplot(df_y, aes(x = date, y = .data[[value_col]], color = species)) +
      geom_point(size = 0.6) +
      scale_color_manual(values = species_palette) +
      
      # 👉 MONTHLY X-AXIS
      scale_x_date(
        date_breaks = "1 month",
        date_labels = "%b"
      ) +
      
      facet_wrap(~species, ncol = 2) +
      base_theme +
      labs(
        title = paste(name, y),
        x = "Month",
        y = name
      )
    
    ggsave(
      filename = paste0("Figures/Hoelstein/", name, "_", y, ".png"),
      plot = p,
      width = 10,
      height = 6,
      dpi = 300
    )
  }
}

plot_yearly_q90 <- function(df, value_col, name) {
  
  # 1) compute daily Q90 from raw data
  df_q90_daily <- df %>%
    group_by(species, year, date) %>%
    mutate(
      q90 = quantile(.data[[value_col]], 0.9, na.rm = TRUE)
    ) %>%
    filter(.data[[value_col]] >= q90) %>%
    summarise(
      value = mean(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  years <- unique(df_q90_daily$year)
  
  # 2) plot
  for (y in years) {
    
    df_y <- df_q90_daily %>%
      filter(year == y)
    
    p <- ggplot(df_y, aes(x = date, y = value, color = species)) +
      geom_point(size = 0.6) +
      scale_color_manual(values = species_palette) +
      scale_x_date(date_breaks = "1 month", date_labels = "%b") +
      facet_wrap(~species, ncol = 2) +
      base_theme +
      labs(
        title = paste0(name, " (daily Q90) - ", y),
        x = "Month",
        y = paste0(name, " (top 10% within day)")
      )
    
    ggsave(
      filename = paste0("Figures/Hoelstein/", name, "_Q90_", y, ".png"),
      plot = p,
      width = 10,
      height = 6,
      dpi = 300
    )
  }
}

# ============================================================
# YEARLY OUTPUTS
# ============================================================

plot_yearly(sap_flux_gc, "sfd", "SFD")

plot_yearly(sap_flux_gc, "G_asw", "G_asw")

plot_yearly(sap_flux_gc, "G_ms", "G_ms")


plot_yearly_q90(sap_flux_gc, "sfd", "SFD")

plot_yearly_q90(sap_flux_gc, "G_asw", "G_asw")

plot_yearly_q90(sap_flux_gc, "G_ms", "G_ms")

# ============================================================
# YEARLY PLOTS FUNCTION (TREATMENT)
# ============================================================
plot_yearly_treatment <- function(df, value_col, name) {
  
  years <- unique(df$year)
  treatments <- unique(df$treatment)
  
  for (t in treatments) {
    
    for (y in years) {
      
      df_y <- df %>%
        filter(year == y, treatment == t)
      
      p <- ggplot(df_y,
                  aes(x = date, y = .data[[value_col]], color = species)) +
        geom_point(size = 0.6) +
        scale_color_manual(values = species_palette) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b") +
        facet_wrap(~species, ncol = 2) +
        base_theme +
        labs(
          title = paste(name, y, "-", t),
          subtitle = "Control vs Treatment separated",
          x = "Month",
          y = name
        )
      
      ggsave(
        filename = paste0("Figures/Hoelstein/", name, "_", y, "_", t, ".png"),
        plot = p,
        width = 10,
        height = 6,
        dpi = 300
      )
    }
  }
}

plot_yearly_q90_treatment <- function(df, value_col, name) {
  
  # 1) DAILY TOP 10% (within each day!)
  df_q90_daily <- df %>%
    group_by(species, treatment, year, date) %>%
    mutate(
      q90 = quantile(.data[[value_col]], 0.9, na.rm = TRUE)
    ) %>%
    filter(.data[[value_col]] >= q90) %>%
    summarise(
      value = mean(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  years <- unique(df_q90_daily$year)
  treatments <- unique(df_q90_daily$treatment)
  
  # 2) plot
  for (t in treatments) {
    for (y in years) {
      
      df_y <- df_q90_daily %>%
        filter(year == y, treatment == t)
      
      p <- ggplot(df_y,
                  aes(x = date, y = value, color = species)) +
        geom_point(size = 0.6) +
        scale_color_manual(values = species_palette) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b") +
        facet_wrap(~species, ncol = 2) +
        base_theme +
        labs(
          title = paste0(name, " (daily top 10%) - ", y, " - ", t),
          x = "Month",
          y = paste0(name, " (top 10% within day)")
        )
      
      ggsave(
        filename = paste0("Figures/Hoelstein/", name, "_Q90_", y, "_", t, ".png"),
        plot = p,
        width = 10,
        height = 6,
        dpi = 300
      )
    }
  }
}

# ============================================================
# YEARLY OUTPUTS
# ============================================================

plot_yearly_treatment(sap_flux_gc, "sfd", "SFD")

plot_yearly_treatment(sap_flux_gc, "G_asw", "G_asw")

plot_yearly_treatment(sap_flux_gc, "G_ms", "G_ms")


plot_yearly_q90_treatment(sap_flux_gc, "sfd", "SFD")

plot_yearly_q90_treatment(sap_flux_gc, "G_asw", "G_asw")

plot_yearly_q90_treatment(sap_flux_gc, "G_ms", "G_ms")


# ============================================================
# SAVE DATA
# ============================================================

sap_flux_gc <- sap_flux_gc %>%
  mutate(
    timestamp = parse_date_time(timestamp,
                                orders = c("ymd HMS", "ymd"),
                                tz = "UTC"),
    timestamp = format(timestamp, "%Y-%m-%d %H:%M:%S"))

write.csv(sap_flux_gc, file = "SCCII/sap_flux_gc.csv", row.names = FALSE)

head(sap_flux_gc)
