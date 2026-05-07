setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(purrr)
library(ggplot)

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


species_levels <- c("Oak", "Beech", "Spruce", "Pine")

cb_palette <- c(
  Oak      = "#E69F00",
  Beech    = "#0072B2",
  Spruce   = "#009E73",
  Pine     = "#F0E442"
)

sap_flux_gc_hoelstein_meteo <- read.csv("SCCII/sap_flux_gc_hoelstein_meteo.csv")

psi_L_folder <- "SCCII/Leaf_Water_Potential"

files <- list.files(
  psi_L_folder,
  pattern = "Water_potentials_.*_archive\\.txt$",
  full.names = TRUE
)

psi_L <- read.delim(files[1], sep = "\t", header = TRUE)
head(psi_L)

# day/month/year: 2018, 2019, 2020, 2021, 2022, 2023
# month/day/year: 2024, 2025

psiL_hoelstein <- map_dfr(files, function(file) {
  
  year <- str_extract(file, "\\d{4}") %>% as.numeric()
  
  df <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    
    select(date, species_id, tree_id, tree_nr, md_wp_av, pd_wp_av) %>%
    
    mutate(
      tree_nr = as.character(tree_nr),
      tree_id = as.character(tree_id),
      species_id = as.character(species_id)
    ) %>%
    
    filter(species_id %in% c("Fs", "Pa", "Ps", "Qs")) %>%
    
    mutate(
      date = if (year <= 2023) {
        dmy(date)
      } else {
        mdy(date)
      },
      
      species_name = recode(species_id,
                            "Fs" = "Beech",
                            "Pa" = "Spruce",
                            "Ps" = "Pine",
                            "Qs" = "Oak")
    )
  
  return(df)
})

write.csv(psiL_hoelstein, file = "SCCII/psiL_hoelstein.csv", row.names = FALSE)

psiL_hoelstein <- read.csv("SCCII/psiL_hoelstein.csv")

sap_daily <- sap_flux_gc_hoelstein_meteo %>%
  
  # =========================
  # FIX TIME + CREATE DAILY DATE
  # =========================
  mutate(
    timestamp = as.POSIXct(timestamp),
    date = as.Date(timestamp), 
    
    hour_min = hour(timestamp) + minute(timestamp) / 60,
    is_daytime = hour_min >= 7.5 & hour_min <= 19.5
  ) %>%
  
  group_by(date, species, tree_id, tree_nr, treatment) %>%
  
  summarise(
    
    # =========================
    # 🌤️ DAYTIME MEAN (07:30–19:30)
    # =========================
    relhum_crane     = mean(Humid_percent_crane[is_daytime], na.rm = TRUE),
    temp_crane       = mean(Temp_degC_crane[is_daytime], na.rm = TRUE),
    radiation_crane  = mean(Solar_Wm.2_crane[is_daytime], na.rm = TRUE),
    vpd_crane        = mean(vpd[is_daytime], na.rm = TRUE),
    
    temp_meteo      = mean(temp[is_daytime], na.rm = TRUE),
    relhum_meteo    = mean(relhum[is_daytime], na.rm = TRUE),
    wind_meteo      = mean(wind[is_daytime], na.rm = TRUE),
    radiation_meteo = mean(radiation[is_daytime], na.rm = TRUE),
    vpd_meteo       = mean(vpd_meteo[is_daytime], na.rm = TRUE),
    
    # =========================
    # 🌧️ DAILY TOTAL / MAX
    # =========================
    precip_crane   = sum(dRain_mm_cranegap, na.rm = TRUE),
    precip_meteo   = max(precip, na.rm = TRUE),
    
    # =========================
    # 🌳 TOP 10% (high activity)
    # =========================
    sfd = mean(sfd[sfd >= quantile(sfd, 0.9, na.rm = TRUE)], na.rm = TRUE),
    G_asw = mean(G_asw[G_asw >= quantile(G_asw, 0.9, na.rm = TRUE)], na.rm = TRUE),
    G_ms  = mean(G_ms[G_ms >= quantile(G_ms, 0.9, na.rm = TRUE)], na.rm = TRUE),
    
    .groups = "drop"
  )

write.csv(sap_daily, file = "SCCII/sap_daily.csv", row.names = FALSE)

climate_long <- sap_daily %>%
  select(date,
         relhum_crane, relhum_meteo,
         temp_crane, temp_meteo,
         radiation_crane, radiation_meteo,
         vpd_crane, vpd_meteo,
         precip_crane, precip_meteo) %>%
  
  pivot_longer(
    cols = -date,
    names_to = c("variable", "source"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(
    date = as.Date(date),
    year = format(date, "%Y")
  )

p <- ggplot(climate_long, aes(
  x = date,
  y = value,
  color = source,
  group = interaction(source, variable, year)
)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 0.6, size = 0.8) +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("crane" = "dodgerblue", "meteo" = "orange")) +
  theme_minimal() +
  base_theme

ggsave(
  filename = "Figures/Hoelstein/sap_daily/plot_climate_crane_meteo.png",
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
)

##### Remove based climate data #####
### VPD < 0.3 kPa 
### radiation < 150 
### Temp < 14 
### precip > 1

sap_daily_filtered <- sap_daily %>%
  filter(
    # VPD condition (either source must be OK)
    (vpd_crane >= 0.3 & vpd_meteo >= 0.3) &
      
      # radiation condition
      (radiation_crane >= 150 & radiation_meteo >= 150) &
      
      # temperature condition
      (temp_crane >= 15 & temp_meteo >= 15) &
      
      # precipitation condition (exclude heavy rain days)
      (precip_crane <= 1)
  )

climate_long <- sap_daily_filtered %>%
  select(date,
         relhum_crane, relhum_meteo,
         temp_crane, temp_meteo,
         radiation_crane, radiation_meteo,
         vpd_crane, vpd_meteo,
         precip_crane, precip_meteo) %>%
  
  pivot_longer(
    cols = -date,
    names_to = c("variable", "source"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(
    date = as.Date(date),
    year = format(date, "%Y")
  )

p <- ggplot(climate_long, aes(
  x = date,
  y = value,
  color = source,
  group = interaction(source, variable, year)
)) +
  geom_line(alpha = 0.8) +
  geom_point(alpha = 0.6, size = 0.8) +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("crane" = "dodgerblue", "meteo" = "orange")) +
  theme_minimal() +
  base_theme

ggsave(
  filename = "Figures/Hoelstein/sap_daily/plot_climate_crane_meteo_filtered.png",
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
)

##### Plot timeseries Gc #####
sap_long <- sap_daily_filtered %>%
  mutate(
    species = factor(species, levels = species_levels),
    year = format(date, "%Y")
  ) %>%
  select(date, year, species, treatment, sfd, G_asw, G_ms) %>%
  pivot_longer(
    cols = c(sfd, G_asw, G_ms),
    names_to = "variable",
    values_to = "value"
  ) 

ylim_df <- sap_long %>%
  group_by(variable) %>%
  summarise(
    ymin = min(value, na.rm = TRUE),
    ymax = max(value, na.rm = TRUE),
    .groups = "drop"
  )

plot_one <- function(df, yr, var, ylim_df) {
  
  df_sub <- df %>%
    filter(year == yr, variable == var)
  
  lims <- ylim_df %>% filter(variable == var)
  
  ggplot(df_sub, aes(
    x = date,
    y = value,
    color = treatment,
    group = treatment
  )) +
    
    geom_line(alpha = 0.8) +
    geom_point(alpha = 0.5, size = 0.6) +
    
    facet_wrap(~ species, ncol = 2) +
    
    scale_color_manual(values = c(
      control = "dodgerblue",
      treatment = "orange"
    )) +
    
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%b"
    ) +
    
    coord_cartesian(
      ylim = c(lims$ymin, lims$ymax)
    ) +
    
    labs(
      title = paste(var, "-", yr),
      x = "Month",
      y = var,
      color = "Treatment"
    ) +
    
    base_theme +
    
    theme(
      legend.position = "top",
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

years <- unique(sap_long$year)
vars  <- unique(sap_long$variable)

for (yr in years) {
  for (v in vars) {
    
    p <- plot_one(sap_long, yr, v, ylim_df)
    
    ggsave(
      filename = paste0("Figures/Hoelstein/sap_daily/", v, "_", yr, ".png"),
      plot = p,
      width = 10,
      height = 8,
      dpi = 300
    )
  }
}

##### combined sap psiL filtered ####
psiL_hoelstein <- psiL_hoelstein %>%
  mutate(
    date = as.Date(date),
    tree_id = as.character(tree_id),
    species = species_name   # unify name
  )

sap_daily_filtered <- sap_daily_filtered %>%
  mutate(
    date = as.Date(date),
    tree_id = as.character(tree_id),
    species = as.character(species)
  )

sap_psiL_daily_filtered <- sap_daily_filtered %>%
  left_join(
    psiL_hoelstein %>%
      select(date, tree_id, species, md_wp_av, pd_wp_av),
    by = c("tree_id", "species", "date")
  )

psi_plot <- sap_psiL_daily_filtered %>%
  mutate(
    species = factor(species, levels = species_levels)
  ) %>%
  select(date, species, treatment, md_wp_av, pd_wp_av) %>%
  pivot_longer(
    cols = c(md_wp_av, pd_wp_av),
    names_to = "variable",
    values_to = "value"
  )

plot_psi <- function(df, var_name) {
  
  df_sub <- df %>%
    filter(variable == var_name)
  
  ggplot(df_sub, aes(
    x = date,
    y = value,
    color = treatment,
    group = treatment
  )) +
    
    geom_line(alpha = 0.8) +
    geom_point(alpha = 0.6, size = 0.8) +
    
    facet_wrap(~ species, ncol = 2) +
    
    scale_color_manual(values = c(
      control = "dodgerblue",
      treatment = "orange"
    )) +
    
    labs(
      title = paste("Leaf Water Potential:", var_name),
      x = "Date",
      y = "ψL (MPa)",
      color = "Treatment"
    ) +
    
    base_theme +
    
    theme(
      legend.position = "top",
      strip.text = element_text(face = "bold")
    )
}

p_md <- plot_psi(psi_plot, "md_wp_av")
p_pd <- plot_psi(psi_plot, "pd_wp_av")

ggsave(
  "Figures/Hoelstein/sap_daily/psiL_md_wp_av.png",
  p_md,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  "Figures/Hoelstein/sap_daily/psiL_pd_wp_av.png",
  p_pd,
  width = 10,
  height = 7,
  dpi = 300
)

##### plot Gc/Gcmax - midday psiL #####
plot_df <- sap_psiL_daily_filtered %>%
  filter(
    !is.na(md_wp_av),
    !is.na(G_ms)
  ) %>%
  mutate(
    species = factor(species, levels = species_levels)
  )


p_Gms_vs_psiLmd <- ggplot(plot_df, aes(
  x = md_wp_av,
  y = G_ms,
  color = species
)) +
  
  geom_point(alpha = 0.7, size = 1.2) +
  
  facet_wrap(~ species, ncol = 2) +
  
  scale_color_manual(values = species_palette, drop = FALSE) +
  
  labs(
    title = "Stomatal Conductance vs Midday Water Potential",
    x = expression(Psi[L]~"(MPa)"),
    y = expression(G[ms])
  ) +
  
  base_theme +
  
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

print(p_Gms_vs_psiLmd)

ggsave(
  filename = "Figures/Hoelstein/sap_daily/Gms_vs_psiLmd.png",
  plot = p_Gms_vs_psiLmd,
  width = 10,
  height = 7,
  dpi = 300
)

p_Gms_vs_psiLmd_treatment <- ggplot(plot_df, aes(
  x = md_wp_av,
  y = G_ms,
  color = treatment
)) +
  
  geom_point(alpha = 0.7, size = 1.2) +
  
  facet_wrap(~ species, ncol = 2) +
  
  scale_color_manual(values = c(
    control = "dodgerblue",
    treatment = "orange"
  )) +
  
  labs(
    title = "G_ms vs Midday Water Potential",
    x = expression(Psi[L]~"(MPa)"),
    y = expression(G[ms]),
    color = "Treatment"
  ) +
  
  base_theme +
  
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

print(p_Gms_vs_psiLmd_treatment)
ggsave(
  filename = "Figures/Hoelstein/sap_daily/p_Gms_vs_psiLmd_treatment.png",
  plot = p_Gms_vs_psiLmd,
  width = 10,
  height = 7,
  dpi = 300
)

##### G_mx normalize #####
plot_df_norm <- plot_df %>%
  group_by(species) %>%
  mutate(
    G_ms_norm = G_ms / max(G_ms, na.rm = TRUE)
  ) %>%
  ungroup()

plot_df_norm_treat <- plot_df %>%
  group_by(species, treatment) %>%
  mutate(
    G_ms_norm = G_ms / max(G_ms, na.rm = TRUE)
  ) %>%
  ungroup()

p_norm_species <- ggplot(plot_df_norm, aes(
  x = md_wp_av,
  y = G_ms_norm,
  color = species
)) +
  
  geom_point(alpha = 0.7, size = 1.2) +
  
  facet_wrap(~ species, ncol = 2) +
  
  scale_color_manual(values = species_palette, drop = FALSE) +
  
  labs(
    title = "Normalized G_ms vs Midday Water Potential",
    x = expression(Psi[L]~"(MPa)"),
    y = expression(G[ms]~"/"~G[max])
  ) +
  
  base_theme +
  
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

print(p_norm_species)

ggsave(
  filename = "Figures/Hoelstein/sap_daily/Gms_vs_psiLmd_species_norm.png",
  plot = p_norm_species,
  width = 10,
  height = 7,
  dpi = 300
)

p_norm_treat <- ggplot(plot_df_norm_treat, aes(
  x = md_wp_av,
  y = G_ms_norm,
  color = treatment
)) +
  
  geom_point(alpha = 0.7, size = 1.2) +
  
  facet_wrap(~ species, ncol = 2) +
  
  scale_color_manual(values = c(
    control = "dodgerblue",
    treatment = "orange"
  )) +
  
  labs(
    title = "Normalized G_ms vs Midday Water Potential",
    x = expression(Psi[L]~"(MPa)"),
    y = expression(G[ms]~"/"~G[max]),
    color = "Treatment"
  ) +
  
  base_theme +
  
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

print(p_norm_treat)

ggsave(
  filename = "Figures/Hoelstein/sap_daily/Gms_vs_psiLmd_treatment_norm.png",
  plot = p_norm_treat,
  width = 10,
  height = 7,
  dpi = 300
)


##### plot one figure #####
plot_raw <- sap_psiL_daily_filtered %>%
  filter(!is.na(md_wp_av), !is.na(G_ms)) %>%
  mutate(species = factor(species, levels = species_levels))

plot_norm <- sap_psiL_daily_filtered %>%
  filter(!is.na(md_wp_av), !is.na(G_ms)) %>%
  group_by(species) %>%
  mutate(
    G_ms_norm = G_ms / max(G_ms, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(species = factor(species, levels = species_levels)) 
  # filter(G_ms_norm > 0.3)

p_raw <- ggplot(plot_raw, aes(
  x = md_wp_av,
  y = G_ms,
  color = species
)) +
  
  geom_point(alpha = 0.7, size = 1.2) +
  
  scale_color_manual(values = species_palette, drop = FALSE) +
  
  labs(
    title = "G_ms vs Midday Water Potential (Raw)",
    x = expression(Psi[L]~"(MPa)"),
    y = expression(G[ms])
  ) +
  
  base_theme +
  
  theme(
    legend.position = "top"
  )

print(p_raw)
ggsave(
  "Figures/Hoelstein/sap_daily/Gms_vs_psiL_raw_all_species.png",
  p_raw,
  width = 8,
  height = 6,
  dpi = 300
)

p_norm <- ggplot(plot_norm, aes(
  x = md_wp_av,
  y = G_ms_norm,
  color = species
)) +
  
  geom_point(alpha = 0.7, size = 1.2) +
  # geom_smooth() +
  scale_color_manual(values = species_palette, drop = FALSE) +
  
  labs(
    title = "Normalized G_ms vs Midday Water Potential",
    x = expression(Psi[L]~"(MPa)"),
    y = expression(G[ms]~"/"~G[max])
  ) +
  
  base_theme +
  
  theme(
    legend.position = "top"
  )

print(p_norm)

ggsave(
  "Figures/Hoelstein/sap_daily/Gms_vs_psiL_normalized_all_species.png",
  p_norm,
  width = 8,
  height = 6,
  dpi = 300
)


##### plot Gc/Gcmax - predawn psiL #####
plot_df <- sap_psiL_daily_filtered %>%
  filter(
    !is.na(pd_wp_av),
    !is.na(G_ms)
  ) %>%
  mutate(
    species = factor(species, levels = species_levels)
  )

p_raw <- ggplot(plot_df, aes(
  x = pd_wp_av,
  y = G_ms,
  color = species
)) +
  
  geom_point(alpha = 0.7, size = 1.2) +
  
  scale_color_manual(values = species_palette, drop = FALSE) +
  
  labs(
    title = "G_ms vs Predawn Water Potential (Raw)",
    x = expression(Psi[L]~"(predawn, MPa)"),
    y = expression(G[ms])
  ) +
  
  base_theme +
  
  theme(
    legend.position = "top"
  )

print(p_raw)

ggsave(
  "Figures/Hoelstein/sap_daily/Gms_vs_psiL_predawn_raw.png",
  p_raw,
  width = 8,
  height = 6,
  dpi = 300
)

p_treat <- ggplot(plot_df, aes(
  x = pd_wp_av,
  y = G_ms,
  color = treatment
)) +
  
  geom_point(alpha = 0.7, size = 1.2) +
  
  facet_wrap(~ species, ncol = 2) +
  
  scale_color_manual(values = c(
    control = "dodgerblue",
    treatment = "orange"
  )) +
  
  labs(
    title = "G_ms vs Predawn Water Potential",
    x = expression(Psi[L]~"(predawn, MPa)"),
    y = expression(G[ms]),
    color = "Treatment"
  ) +
  
  base_theme +
  
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

print(p_treat)

ggsave(
  "Figures/Hoelstein/sap_daily/Gms_vs_psiL_predawn_treatment.png",
  p_treat,
  width = 10,
  height = 7,
  dpi = 300
)

plot_df_norm <- plot_df %>%
  group_by(species) %>%
  mutate(
    G_ms_norm = G_ms / max(G_ms, na.rm = TRUE)
  ) %>%
  ungroup()

p_norm_species <- ggplot(plot_df_norm, aes(
  x = pd_wp_av,
  y = G_ms_norm,
  color = species
)) +
  
  geom_point(alpha = 0.7, size = 1.2) +
  
  facet_wrap(~ species, ncol = 2) +
  
  scale_color_manual(values = species_palette, drop = FALSE) +
  
  labs(
    title = "Normalized G_ms vs Predawn Water Potential",
    x = expression(Psi[L]~"(predawn, MPa)"),
    y = expression(G[ms]~"/"~G[max])
  ) +
  
  base_theme +
  
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

print(p_norm_species)

ggsave(
  "Figures/Hoelstein/sap_daily/Gms_vs_psiL_predawn_species_norm.png",
  p_norm_species,
  width = 10,
  height = 7,
  dpi = 300
)

plot_df_norm_treat <- plot_df %>%
  group_by(species, treatment) %>%
  mutate(
    G_ms_norm = G_ms / max(G_ms, na.rm = TRUE)
  ) %>%
  ungroup()


p_norm_treat <- ggplot(plot_df_norm_treat, aes(
  x = pd_wp_av,
  y = G_ms_norm,
  color = treatment
)) +
  
  geom_point(alpha = 0.7, size = 1.2) +
  
  facet_wrap(~ species, ncol = 2) +
  
  scale_color_manual(values = c(
    control = "dodgerblue",
    treatment = "orange"
  )) +
  
  labs(
    title = "Normalized G_ms vs Predawn Water Potential",
    x = expression(Psi[L]~"(predawn, MPa)"),
    y = expression(G[ms]~"/"~G[max]),
    color = "Treatment"
  ) +
  
  base_theme +
  
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

print(p_norm_treat)

ggsave(
  "Figures/Hoelstein/sap_daily/Gms_vs_psiL_predawn_treatment_norm.png",
  p_norm_treat,
  width = 10,
  height = 7,
  dpi = 300
)