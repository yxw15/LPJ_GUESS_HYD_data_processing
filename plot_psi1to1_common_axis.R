# ==========================================================================
# COMMON AXIS LIMITS FOR BOTH SCATTER PLOTS
# Calculate once from ALL data (both control and drought)
# ==========================================================================

# Load both datasets (assuming you have them)
lpj_control   <- read.csv("lpj_guess/lpj_Gc_psiL_psiX_psiS.csv") %>% mutate(date = as.Date(date))
lpj_drought   <- read.csv("lpj_guess/lpj_Gc_psiL_psiX_psiS_drought.csv") %>% mutate(date = as.Date(date))

psiL_control <- read.csv("SCCII/psiL_hoelstein_control.csv") %>% mutate(date = as.Date(date))
psiL_drought <- read.csv("SCCII/psiL_hoelstein_drought.csv") %>% mutate(date = as.Date(date))

psiS_control <- read.csv("SCCII/psiS_hoelstein_control.csv") %>% 
  mutate(date = as.Date(date), psiS_mean = psiS_mean / 1000)
psiS_drought <- read.csv("SCCII/psiS_hoelstein_drought.csv") %>% 
  mutate(date = as.Date(date), psiS_mean = psiS_mean / 1000)

# Process both datasets through the same preparation steps to get all values
prepare_scatter_data <- function(lpj_data, obs_psiL, obs_psiS) {
  # Aggregate observations
  obs_psiL_daily <- obs_psiL %>%
    group_by(date, species_name) %>%
    summarise(
      psiLmd_obs = mean(md_wp_av, na.rm = TRUE),
      psiLpd_obs = mean(pd_wp_av, na.rm = TRUE),
      .groups = "drop"
    ) %>% rename(species = species_name)
  
  obs_psiS_daily <- obs_psiS %>%
    group_by(date) %>%
    summarise(psiS_obs = mean(psiS_mean, na.rm = TRUE), .groups = "drop")
  
  data_longterm <- lpj_data %>%
    filter(date >= "2018-01-01" & date <= "2025-12-31") %>%
    left_join(obs_psiL_daily, by = c("date", "species")) %>%
    left_join(obs_psiS_daily, by = "date")
  
  # Create comparison pairs
  plot_comparison <- data_longterm %>%
    filter(!is.na(psiLmd_obs) | !is.na(psiLpd_obs) | !is.na(psiS_obs)) %>%
    pivot_longer(
      cols = c(psiLmd_obs, psiLpd_obs, psiS_obs),
      names_to = "obs_type",
      values_to = "obs_value"
    ) %>%
    mutate(model_value = case_when(
      obs_type == "psiLmd_obs" ~ psiL,
      obs_type == "psiLpd_obs" ~ psiL,
      obs_type == "psiS_obs"   ~ psiS,
      TRUE ~ NA_real_
    )) %>%
    bind_rows(
      data_longterm %>%
        filter(!is.na(psiLpd_obs)) %>%
        transmute(species, date, obs_type = "psiX_vs_pd_obs", 
                  obs_value = psiLpd_obs, model_value = psiX),
      data_longterm %>%
        filter(!is.na(psiLmd_obs)) %>%
        transmute(species, date, obs_type = "psiX_vs_md_obs", 
                  obs_value = psiLmd_obs, model_value = psiX)
    ) %>%
    filter(!is.na(obs_value) & !is.na(model_value))
  
  return(plot_comparison)
}

# Prepare both datasets
control_comparison <- prepare_scatter_data(lpj_control, psiL_control, psiS_control)
drought_comparison <- prepare_scatter_data(lpj_drought, psiL_drought, psiS_drought)

# Combine all values to find common range
all_control_values <- c(control_comparison$obs_value, control_comparison$model_value)
all_drought_values <- c(drought_comparison$obs_value, drought_comparison$model_value)
all_values_combined <- c(all_control_values, all_drought_values)

# Calculate COMMON limits for both plots
plot_range <- range(all_values_combined, na.rm = TRUE)
buffer <- diff(plot_range) * 0.05
COMMON_LIMITS <- c(plot_range[1] - buffer, plot_range[2] + buffer)

# Also ensure ylim is consistent with these limits
Y_LIMITS <- c(min(COMMON_LIMITS[1], -4), max(COMMON_LIMITS[2], 0))

# ==========================================================================
# MODIFIED SCATTER PLOT FUNCTION WITH FIXED LIMITS
# ==========================================================================

create_scatter_plot <- function(plot_comparison, title_subtitle, filename) {
  
  # Ensure factor levels are correct
  plot_comparison <- plot_comparison %>%
    mutate(
      species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")),
      obs_type = factor(obs_type, levels = c("psiLmd_obs", "psiLpd_obs", "psiS_obs", 
                                             "psiX_vs_md_obs", "psiX_vs_pd_obs"))
    )
  
  p <- ggplot(plot_comparison, aes(x = obs_value, y = model_value, 
                                   color = obs_type, shape = obs_type)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_point(size = 2.5, alpha = 0.7) +
    facet_wrap(~species, ncol = 4) +
    
    # CRITICAL: Use EXACT same limits for both plots
    coord_equal(xlim = COMMON_LIMITS, ylim = COMMON_LIMITS) +
    
    scale_color_manual(
      values = c(
        "psiLmd_obs"     = "green4", 
        "psiLpd_obs"     = "purple", 
        "psiS_obs"       = "orange", 
        "psiX_vs_md_obs" = "dodgerblue", 
        "psiX_vs_pd_obs" = "pink"
      ),
      labels = c(
        "psiLmd_obs"     = expression(Psi["   L"]~vs~Obs~MD),
        "psiLpd_obs"     = expression(Psi["   L"]~vs~Obs~PD),
        "psiS_obs"       = expression(Psi["   S"]~vs~Obs~Soil),
        "psiX_vs_md_obs" = expression(Psi["   X"]~vs~Obs~MD),
        "psiX_vs_pd_obs" = expression(Psi["   X"]~vs~Obs~PD)
      )
    ) +
    scale_shape_manual(
      values = c(17, 16, 15, 2, 1), 
      labels = c(
        "psiLmd_obs"     = expression(Psi["   L"]~vs~Obs~MD),
        "psiLpd_obs"     = expression(Psi["   L"]~vs~Obs~PD),
        "psiS_obs"       = expression(Psi["   S"]~vs~Obs~Soil),
        "psiX_vs_md_obs" = expression(Psi["   X"]~vs~Obs~MD),
        "psiX_vs_pd_obs" = expression(Psi["   X"]~vs~Obs~PD)
      )
    ) +
    
    labs(
      title = title_subtitle,
      subtitle = "dashed = 1:1 reference line | identical axes between control and drought plots",
      x = expression(Observed~Psi~"(MPa)"),
      y = expression(Simulated~LPJ-GUESS~Psi~"(MPa)"),
      color = "Comparison Pair",
      shape = "Comparison Pair"
    ) +
    
    base_theme +
    theme(
      legend.position = "bottom", 
      legend.text.align = 0,
      legend.box = "vertical",
      panel.grid.minor = element_blank(),
      aspect.ratio = 1,
      plot.margin = margin(20, 20, 20, 20)
    )
  
  print(p)
  ggsave(filename, p, width = 16, height = 6, bg = "white")
  
  return(p)
}

# ==========================================================================
# CREATE BOTH PLOTS WITH IDENTICAL AXES
# ==========================================================================

# Print the common limits to verify
cat("Common axis limits for both plots:", round(COMMON_LIMITS, 2), "\n")

# Create control plot
p_control <- create_scatter_plot(
  control_comparison,
  "water potential model vs observation: 1:1 scatter plot (control)",
  "Figures/Hoelstein/psi_1to1.png"
)

# Create drought plot
p_drought <- create_scatter_plot(
  drought_comparison,
  "water potential model vs observation: 1:1 scatter plot (drought experiment since 2023)",
  "Figures/Hoelstein/psi_1to1_drought.png"
)
