# ==============================================================================
# 0. SETUP & PATHS
# ==============================================================================
library(tidyverse)

base_dir <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD"
setwd(base_dir)

out_dir  <- "Figures/lpj_guess_hyd"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. LOAD AND MERGE DATA
# ==============================================================================
mort_df   <- read_csv("lpj_guess/lpj_guess_hyd/lpj_control_drought_mort_kappa_annual_full.csv")
carbon_df <- read_csv("lpj_guess/lpj_guess_hyd/lpj_control_drought_yearly_carbon_productivity_full.csv")

# Filter and Merge: Ensure years < 1991 are removed immediately
master_annual_wide <- mort_df %>%
  select(-any_of("date")) %>%
  filter(year >= 1991) %>% 
  full_join(carbon_df, by = c("treatment", "species", "year")) %>%
  filter(year >= 1991) # Double-check filter after join

# Define range based on filtered data
min_yr <- 1991
max_yr <- max(master_annual_wide$year, na.rm = TRUE)

# ==============================================================================
# 2. METADATA & THEME
# ==============================================================================
species_levels <- c("oak", "beech", "spruce", "pine")
species_colors <- c("oak" = "#E69F00", "beech" = "#0072B2", "spruce" = "#009E73", "pine" = "#F0E442")
treatments     <- c("control", "drought")

var_registry <- tribble(
  ~var_name,       ~plot_title,                       ~ylab_text,
  "mort",          "total mortality rate",            "mortality rate (yr^-1)",
  "mort_min",      "background mortality rate",       "mortality rate (yr^-1)",
  "mort_greff",    "growth-efficiency mortality rate", "mortality rate (yr^-1)",
  "mort_cav",      "hydraulic-failure mortality rate", "mortality rate (yr^-1)",
  "kappa_s_min",   "maximum value of cavitated xylems", "cavitation fraction",
  "kappa_s_today", "fraction of cavitated xylems",     "cavitation fraction",
  "cmass",         "vegetation carbon biomass (cmass)", "carbon mass (kg c m^-2)",
  "cmass_mort",    "carbon mass lost to mortality",    "carbon mass loss (kg c m^-2 yr^-1)",
  "agpp",          "gross primary productivity (agpp)", "gpp (kg c m^-2 yr^-1)",
  "anpp",          "net primary productivity (anpp)",   "npp (kg c m^-2 yr^-1)"
)

base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.position   = "bottom",
    plot.title        = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title        = element_text(size = 12),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    strip.text        = element_text(size = 11, face = "bold")
  )

# ==============================================================================
# 3. PLOTTING LOOP
# ==============================================================================
plot_ready_df <- master_annual_wide %>%
  mutate(
    species   = factor(tolower(species), levels = species_levels),
    treatment = factor(treatment, levels = treatments)
  )

for(i in seq_len(nrow(var_registry))) {
  target_var  <- var_registry$var_name[i]
  title_text  <- str_to_title(var_registry$plot_title[i])
  y_axis_text <- str_to_sentence(var_registry$ylab_text[i])
  
  if (!target_var %in% names(plot_ready_df)) next
  
  message("Generating: ", target_var)
  
  p <- ggplot(plot_ready_df, aes(x = year, y = .data[[target_var]], color = species)) +
    geom_line(linewidth = 0.9, alpha = 0.8) +
    geom_point(size = 1.6) +
    facet_grid(species ~ treatment, scales = "free_y") + 
    scale_color_manual(values = species_colors) +
    # Force axis to start at 1991 and use a fixed number of breaks
    scale_x_continuous(
      limits = c(min_yr, max_yr), 
      breaks = scales::pretty_breaks(n = 6) 
    ) +
    labs(title = paste0(title_text, " (", min_yr, " - ", max_yr, ")"),
         x = "Year", y = y_axis_text) +
    base_theme +
    theme(legend.position = "none")
  
  suppressWarnings(
    ggsave(file.path(out_dir, paste0("ts_matrix_", target_var, ".png")), 
           plot = p, width = 11, height = 10, dpi = 300)
  )
}

