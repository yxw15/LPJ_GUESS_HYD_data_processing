# ==========================================================================
# 1. SETUP, THEME, & PATHS
# ==========================================================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(stringr)
library(purrr)
library(scales)
library(hydroGOF)

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")
dir.create("Figures/lpj_guess_stem_storage/TWD", recursive = TRUE, showWarnings = FALSE)

species_levels <- c("oak", "beech", "spruce", "pine")
species_colors <- c("oak" = "#E69F00", "beech" = "#0072B2", "spruce" = "#009E73", "pine" = "#F0E442")

base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text       = element_text(color = "black", size = 11),
    legend.position   = "bottom",
    legend.box        = "vertical",
    plot.title        = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5, size = 10, color = "grey30"),
    axis.title        = element_text(size = 12),
    axis.text.x       = element_text(angle = 0, hjust = 0.5, size = 9),
    strip.text        = element_text(size = 11, face = "bold"),
    panel.grid.major  = element_line(color = "grey92"),
    panel.grid.minor  = element_blank()
  )

# ==========================================================================
# 2. DATA PREPARATION
# ==========================================================================
lpj_raw_twd <- read.csv("lpj_guess/lpj_guess_stem_storage/lpj_control_drought_ET_Gc_psi_leaf_psi_soil_psi_xylem_hydraulic_lag_kappy_s_min_mort_mort_cav_mort_greff_mort_min_stem_diameter_stem_rwc_twd.csv") %>%
  mutate(date = as.Date(date), species = factor(tolower(species), levels = species_levels), twd_model = twd * 1e6) %>%
  filter(!is.na(species)) %>% select(date, species, treatment, twd_model)

dendro_obs <- list.files(path = "SCCII/point_dendro", pattern = "^Point_dendrometers_.*_archive\\.txt$", full.names = TRUE) %>%
  map_dfr(~ read.delim(.x), .id = "source_file")
tree_info <- read.csv("SCCII/tree_info.csv") %>% mutate(treatment = ifelse(treatment == "treatment", "drought", treatment))

obs_twd_daily <- dendro_obs %>% inner_join(tree_info, by = "tree_id") %>%
  mutate(date = as.Date(str_extract(timestamp_UTC, "\\d{4}-\\d{2}-\\d{2}")), species = factor(tolower(species), levels = species_levels), treatment = tolower(treatment)) %>%
  group_by(date, species, treatment) %>% summarise(twd_mean_obs = mean(twd_micron_treenetproc, na.rm = TRUE), .groups = "drop")

combined_twd <- obs_twd_daily %>% inner_join(lpj_raw_twd, by = c("date", "species", "treatment"))

# Standardize and Calculate KGE
combined_std <- combined_twd %>% group_by(species, treatment) %>%
  mutate(std_twd_mean = (twd_mean_obs - mean(twd_mean_obs, na.rm=T))/sd(twd_mean_obs, na.rm=T),
         std_twd_model = (twd_model - mean(twd_model, na.rm=T))/sd(twd_model, na.rm=T)) %>%
  ungroup()

stats_summary <- combined_std %>% filter(!is.na(std_twd_mean), !is.na(std_twd_model)) %>%
  group_by(species, treatment) %>%
  summarise(
    kge = KGE(std_twd_model, std_twd_mean, na.rm = TRUE),
    r2 = cor(std_twd_mean, std_twd_model)^2,
    .groups = "drop"
  ) %>% mutate(label = paste0("KGE: ", round(kge, 2), "\nR²: ", round(r2, 2)))

# ==========================================================================
# 3. PLOTTING
# ==========================================================================

# Time Series
ts_plot <- ggplot(combined_std, aes(x = date)) +
  geom_line(aes(y = std_twd_model, color = species), alpha = 0.7) +
  geom_line(aes(y = std_twd_mean), color = "black", linetype = "dashed", alpha = 0.4) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = species_colors) +
  labs(title = "Standardized TWD Time Series", y = "Z-score", x = NULL) +
  base_theme

# Scatter
scatter_plot <- ggplot(combined_std, aes(x = std_twd_mean, y = std_twd_model, color = species)) +
  geom_point(alpha = 0.3, size = 1) +
  # 1:1 reference line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_smooth(method = "lm", color = "black", se = FALSE, linewidth = 0.5) +
  geom_text(data = stats_summary, aes(x = -Inf, y = Inf, label = label), 
            hjust = -0.1, vjust = 1.1, size = 3, color = "black", inherit.aes = FALSE) +
  facet_grid(treatment ~ species) +
  # coord_fixed ensures 1:1 ratio and keeps axes square
  coord_fixed(ratio = 1, xlim = c(-3, 3), ylim = c(-3, 3)) +
  scale_color_manual(values = species_colors) +
  labs(title = "Standardized TWD Comparison", x = "Observed (Z)", y = "Simulated (Z)") +
  base_theme

ggsave("Figures/lpj_guess_stem_storage/TWD/ts_standardized_all.png", ts_plot, width = 10, height = 6, dpi = 300)
ggsave("Figures/lpj_guess_stem_storage/TWD/scatter_standardized_all.png", scatter_plot, width = 10, height = 6, dpi = 300)

# ==========================================================================
# 4. MONTHLY MEAN AGGREGATION & PLOTS
# ==========================================================================

# Aggregate to monthly means
combined_monthly <- combined_twd %>%
  mutate(year_mon = floor_date(date, "month")) %>%
  group_by(species, treatment, year_mon) %>%
  summarise(
    twd_mean_obs = mean(twd_mean_obs, na.rm = TRUE),
    twd_model    = mean(twd_model, na.rm = TRUE),
    n_days       = n(),
    .groups      = "drop"
  ) %>%
  filter(n_days >= 5)  # require at least 5 days per month

# Re-standardize at monthly level
combined_monthly_std <- combined_monthly %>%
  group_by(species, treatment) %>%
  mutate(
    std_twd_mean  = (twd_mean_obs - mean(twd_mean_obs, na.rm = TRUE)) /
                     sd(twd_mean_obs, na.rm = TRUE),
    std_twd_model = (twd_model - mean(twd_model, na.rm = TRUE)) /
                     sd(twd_model, na.rm = TRUE)
  ) %>%
  ungroup()

# Monthly-level KGE and R² statistics
stats_monthly <- combined_monthly_std %>%
  filter(!is.na(std_twd_mean), !is.na(std_twd_model)) %>%
  group_by(species, treatment) %>%
  summarise(
    kge = KGE(std_twd_model, std_twd_mean, na.rm = TRUE),
    r2  = cor(std_twd_mean, std_twd_model)^2,
    n   = n(),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("KGE: ", round(kge, 2), "\nR²: ", round(r2, 2), "\nn: ", n, " months"))

# Monthly Time Series
ts_monthly_plot <- ggplot(combined_monthly_std, aes(x = year_mon)) +
  geom_line(aes(y = std_twd_model, color = species), alpha = 0.7) +
  geom_line(aes(y = std_twd_mean), color = "black", linetype = "dashed", alpha = 0.4) +
  facet_grid(treatment ~ species) +
  scale_color_manual(values = species_colors) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  labs(
    title    = "Standardized TWD Time Series (Monthly Means)",
    subtitle = "colored = simulated | black dashed = observed",
    y        = "Z-score",
    x        = NULL
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# Monthly Scatter
scatter_monthly_plot <- ggplot(combined_monthly_std,
                                aes(x = std_twd_mean, y = std_twd_model, color = species)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.6) +
  geom_smooth(method = "lm", color = "black", se = FALSE, linewidth = 0.5) +
  geom_text(data = stats_monthly,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.1, size = 3, color = "black",
            inherit.aes = FALSE) +
  facet_grid(treatment ~ species) +
  coord_fixed(ratio = 1, xlim = c(-3, 3), ylim = c(-3, 3)) +
  scale_color_manual(values = species_colors) +
  labs(
    title    = "Standardized TWD Comparison (Monthly Means)",
    subtitle = "dashed = 1:1 line | solid = linear regression",
    x        = "Observed (Z)",
    y        = "Simulated (Z)"
  ) +
  base_theme

# Save monthly plots
ggsave("Figures/lpj_guess_stem_storage/TWD/ts_standardized_monthly.png",
       ts_monthly_plot, width = 10, height = 6, dpi = 300)
ggsave("Figures/lpj_guess_stem_storage/TWD/scatter_standardized_monthly.png",
       scatter_monthly_plot, width = 10, height = 6, dpi = 300)

# Display monthly plots
print(ts_monthly_plot)
print(scatter_monthly_plot)

# ==========================================================================
# 5. STL DECOMPOSITION, THEIL-SEN TREND, MANN-KENDALL & KENDALL-TAU
# ==========================================================================

pkg_needed <- c("Kendall", "mblm", "zoo")
for (pkg in pkg_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", lib = user_lib)
  }
}
library(Kendall)   # Kendall() tau-b correlation, MannKendall()
library(mblm)      # mblm() Theil-Sen linear model for scatter plots
library(zoo)       # na.approx() for gap-filling
library(gridExtra) # arrangeGrob for multi-panel layouts
library(grid)      # textGrob for annotation

# 5b. Helper functions --------------------------------------------------------

# Safely compute STL decomposition on a monthly time series (freq = 6, semi-annual)
compute_stl <- function(x, freq = 3) {
  if (length(na.omit(x)) < 2 * freq) return(NULL)
  if (sd(x, na.rm = TRUE) < 1e-12)     return(NULL)
  tryCatch({
    stl(ts(x, frequency = freq), s.window = "periodic", robust = TRUE)
  }, error = function(e) NULL)
}

# Theil-Sen slope (units / month) and Mann-Kendall p-value
# Uses manual Theil-Sen (median of pairwise slopes) + Kendall::MannKendall
calc_theilsen_mk <- function(x) {
  x_clean <- na.omit(x)
  n <- length(x_clean)
  if (n < 3) return(c(slope = NA_real_, p_value = NA_real_))
  tryCatch({
    # Theil-Sen: median of all pairwise slopes
    v <- seq_len(n)
    pairs <- expand.grid(i = v, j = v)
    pairs <- pairs[pairs$i < pairs$j, ]
    slopes <- (x_clean[pairs$j] - x_clean[pairs$i]) /
              (pairs$j - pairs$i)
    ts_slope <- median(slopes, na.rm = TRUE)

    # Mann-Kendall test via Kendall package
    mk <- Kendall::MannKendall(x_clean)
    c(slope = ts_slope, p_value = mk$sl)
  }, error = function(e) c(slope = NA_real_, p_value = NA_real_))
}

# STL variance decomposition: proportion of var(trend), var(seasonal), var(remainder)
stl_var_props <- function(decomp) {
  if (is.null(decomp)) return(c(var_trend = NA_real_, var_seasonal = NA_real_, var_remainder = NA_real_))
  ts_mat <- decomp$time.series
  vars  <- apply(ts_mat, 2, var)
  props <- vars / sum(vars)
  # match by column name to ensure correct assignment
  out <- c(var_trend = NA_real_, var_seasonal = NA_real_, var_remainder = NA_real_)
  for (nm in names(props)) {
    out[paste0("var_", nm)] <- props[nm]
  }
  return(out)
}

# Significance star from p-value
sig_star <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*",
          ifelse(p < 0.1, ".", "")))))
}

# 5c. Prepare monthly time series per species x treatment ----------------------

# Build complete monthly date grid for each group, fill small gaps with
# interpolation so STL receives a regular series.
# Use standardized monthly TWD (Z-score within species × treatment)
combined_monthly_filled <- combined_monthly_std %>%
  group_by(species, treatment) %>%
  summarise(
    date_start = min(year_mon), date_end = max(year_mon),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    date_seq = list(seq.Date(date_start, date_end, by = "month"))
  ) %>%
  unnest(date_seq) %>%
  rename(year_mon = date_seq) %>%
  left_join(combined_monthly_std, by = c("species", "treatment", "year_mon")) %>%
  group_by(species, treatment) %>%
  mutate(
    std_twd_mean  = na.approx(std_twd_mean,  maxgap = 3, na.rm = FALSE),
    std_twd_model = na.approx(std_twd_model, maxgap = 3, na.rm = FALSE),
    # fall back to carry-forward for remaining edges
    std_twd_mean  = na.locf(std_twd_mean,  na.rm = FALSE),
    std_twd_model = na.locf(std_twd_model, na.rm = FALSE)
  ) %>%
  ungroup() %>%
  filter(!is.na(std_twd_mean), !is.na(std_twd_model))

# Nest for per-group processing
nested_data <- combined_monthly_filled %>%
  group_by(species, treatment) %>%
  summarise(
    n_months = n(),
    year_mon = list(year_mon),
    obs_raw  = list(std_twd_mean),
    sim_raw  = list(std_twd_model),
    .groups  = "drop"
  ) %>%
  filter(n_months >= 24)  # require at least 2 years

# 5d. Run STL, Theil-Sen, Mann-Kendall, and Kendall-Tau per group -------------

stl_plot_data  <- list()   # holds long-format data frames for plotting
all_stats      <- list()   # holds one-row summary stats per group

for (i in seq_len(nrow(nested_data))) {
  sp     <- as.character(nested_data$species[i])
  trt    <- as.character(nested_data$treatment[i])
  grp_id <- paste(sp, trt, sep = "_")

  dates  <- nested_data$year_mon[[i]]
  obs    <- nested_data$obs_raw[[i]]
  sim    <- nested_data$sim_raw[[i]]

  # --- STL decomposition ---
  stl_obs <- compute_stl(obs)
  stl_sim <- compute_stl(sim)

  # Extract components into data frames
  if (!is.null(stl_obs)) {
    df_obs <- as.data.frame(stl_obs$time.series)
    df_obs$date      <- dates
    df_obs$source    <- "observed"
    df_obs$data      <- obs
  }
  if (!is.null(stl_sim)) {
    df_sim <- as.data.frame(stl_sim$time.series)
    df_sim$date      <- dates
    df_sim$source    <- "simulated"
    df_sim$data      <- sim
  }

  if (!is.null(stl_obs) && !is.null(stl_sim)) {
    stl_long <- bind_rows(df_obs, df_sim) %>%
      pivot_longer(
        cols      = c(trend, seasonal, remainder, data),
        names_to  = "component",
        values_to = "value"
      ) %>%
      mutate(
        component = factor(component,
          levels = c("data", "trend", "seasonal", "remainder"),
          labels = c("Data", "Trend", "Seasonal", "Remainder")),
        species   = sp,
        treatment = trt
      )
    stl_plot_data[[grp_id]] <- stl_long
  }

  # --- Theil-Sen & Mann-Kendall on raw monthly data ---
  ts_obs <- calc_theilsen_mk(obs)
  ts_sim <- calc_theilsen_mk(sim)

  # --- Theil-Sen & Mann-Kendall on STL trend component ---
  ts_trend_obs <- if (!is.null(stl_obs)) {
    calc_theilsen_mk(stl_obs$time.series[, "trend"])
  } else c(slope = NA_real_, p_value = NA_real_)
  ts_trend_sim <- if (!is.null(stl_sim)) {
    calc_theilsen_mk(stl_sim$time.series[, "trend"])
  } else c(slope = NA_real_, p_value = NA_real_)

  # --- Kendall-Tau cross-correlation (obs vs sim) ---
  kt <- tryCatch({
    kr <- Kendall(obs, sim)
    c(tau = kr$tau, p_value = kr$sl)
  }, error = function(e) c(tau = NA_real_, p_value = NA_real_))

  # --- STL variance proportions ---
  vp_obs <- stl_var_props(stl_obs)
  vp_sim <- stl_var_props(stl_sim)

  # --- Assemble stats row ---
  all_stats[[grp_id]] <- data.frame(
    species            = sp,
    treatment          = trt,
    n_months           = nested_data$n_months[i],
    # raw series Theil-Sen (units/month) & MK p
    ts_slope_obs       = ts_obs["slope"],
    ts_p_obs           = ts_obs["p_value"],
    ts_slope_sim       = ts_sim["slope"],
    ts_p_sim           = ts_sim["p_value"],
    # STL trend-component Theil-Sen (units/month) & MK p
    ts_trend_slope_obs = ts_trend_obs["slope"],
    ts_trend_p_obs     = ts_trend_obs["p_value"],
    ts_trend_slope_sim = ts_trend_sim["slope"],
    ts_trend_p_sim     = ts_trend_sim["p_value"],
    # STL variance proportions (obs)
    stl_var_trend_obs     = vp_obs["var_trend"],
    stl_var_seasonal_obs  = vp_obs["var_seasonal"],
    stl_var_remainder_obs = vp_obs["var_remainder"],
    # STL variance proportions (sim)
    stl_var_trend_sim     = vp_sim["var_trend"],
    stl_var_seasonal_sim  = vp_sim["var_seasonal"],
    stl_var_remainder_sim = vp_sim["var_remainder"],
    # Kendall tau
    kendall_tau        = kt["tau"],
    kendall_p          = kt["p_value"],
    stringsAsFactors   = FALSE
  )
}

trend_stats <- bind_rows(all_stats) %>%
  mutate(species = factor(species, levels = species_levels))

# 5e. STL component plots (one per species × treatment) -----------------------

stl_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.position   = "bottom",
    legend.text       = element_text(size = 10),
    plot.title        = element_text(hjust = 0.5, size = 13, face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5, size = 9, color = "grey30"),
    strip.text        = element_text(size = 11, face = "bold"),
    axis.title        = element_text(size = 10),
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.major  = element_line(color = "grey92"),
    panel.grid.minor  = element_blank()
  )

for (grp_id in names(stl_plot_data)) {
  df_plot <- stl_plot_data[[grp_id]]
  sp      <- df_plot$species[1]
  trt     <- df_plot$treatment[1]

  # Stats for annotation
  st_row <- trend_stats %>% filter(species == sp, treatment == trt)

  # Build annotation string
  anno_text <- paste0(
    "Theil-Sen slope (obs): ",      round(st_row$ts_slope_obs, 4), " Z/mo",
    sig_star(st_row$ts_p_obs), "\n",
    "Theil-Sen slope (sim): ",      round(st_row$ts_slope_sim, 4), " Z/mo",
    sig_star(st_row$ts_p_sim), "\n",
    "Kendall τ = ",           round(st_row$kendall_tau, 3),
    ifelse(is.na(st_row$kendall_p), "",
           paste0(" (p = ", format.pval(st_row$kendall_p, digits = 2), ")"))
  )

  p <- ggplot(df_plot, aes(x = date, y = value, color = source, linetype = source)) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~ component, scales = "free_y", ncol = 1) +
    scale_color_manual(
      values = c("observed" = "grey50", "simulated" = unname(species_colors[sp])),
      name   = NULL,
      labels = c("observed" = "Observed", "simulated" = "Simulated")
    ) +
    scale_linetype_manual(
      values = c("observed" = "dashed", "simulated" = "solid"),
      name   = NULL
    ) +
    labs(
      title    = paste0("STL Decomposition — ", str_to_title(sp), " (", trt, ")"),
      subtitle = anno_text,
      x        = NULL, y = "TWD (Z-score)"
    ) +
    stl_theme +
    theme(plot.subtitle = element_text(hjust = 0, size = 8, color = "grey20",
                                        face = "plain", family = "mono"))

  fname <- paste0("Figures/lpj_guess_stem_storage/TWD/stl_", grp_id, ".png")
  ggsave(fname, p, width = 10, height = 8, dpi = 300)
}

cat("\nSaved", length(stl_plot_data), "STL component figures.\n")

# 5e2. Combined TREND-only figure (all species, one panel per treatment) --------

# Combine trend component across all groups
trend_combined <- bind_rows(lapply(names(stl_plot_data), function(grp_id) {
  stl_plot_data[[grp_id]] %>% filter(component == "Trend")
})) %>%
  mutate(species = factor(species, levels = species_levels))

# Color group: observed always grey, simulated uses species_colors
trend_color_map <- c("observed" = "grey50", species_colors)
trend_combined <- trend_combined %>%
  mutate(color_group = ifelse(source == "observed", "observed", as.character(species)))

# Per-panel annotation from trend_stats (trend-component Theil-Sen + Kendall τ)
trend_anno <- trend_stats %>%
  mutate(
    anno = paste0(
      "TS trend (obs): ", round(ts_trend_slope_obs, 4), " Z/mo",
      sig_star(ts_trend_p_obs), "  |  ",
      "TS trend (sim): ", round(ts_trend_slope_sim, 4), " Z/mo",
      sig_star(ts_trend_p_sim), "\n",
      "Kendall τ = ", round(kendall_tau, 3),
      ifelse(is.na(kendall_p), "",
             paste0(" (p = ", format.pval(kendall_p, digits = 2), ")"))
    )
  )

# Plot: species = panels, observed = grey, simulated = species color,
#       treatment = linetype (control = dashed, drought = solid)
p_trend_combined <- ggplot(trend_combined,
  aes(x = date, y = value, color = color_group, linetype = treatment)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ species, ncol = 2,
    labeller = labeller(species = c("oak" = "Oak", "beech" = "Beech",
                                     "spruce" = "Spruce", "pine" = "Pine"))) +
  scale_color_manual(
    values = trend_color_map,
    name   = NULL,
    breaks = c("observed", names(species_colors)),
    labels = c("observed" = "Observed", "oak" = "Oak", "beech" = "Beech",
               "spruce" = "Spruce", "pine" = "Pine")
  ) +
  scale_linetype_manual(
    values = c("control" = "dashed", "drought" = "solid"),
    name   = NULL,
    labels = c("control" = "Control", "drought" = "Drought")
  ) +
  labs(
    title    = "STL Trend Component — Standardized TWD",
    subtitle = "grey = observed  |  colored = simulated  |  dashed = control  |  solid = drought",
    x        = NULL,
    y        = "TWD Trend (Z-score)"
  ) +
  stl_theme +
  theme(
    legend.position = "bottom",
    legend.box     = "vertical"
  )

ggsave("Figures/lpj_guess_stem_storage/TWD/stl_trend_combined.png",
       p_trend_combined, width = 12, height = 8, dpi = 300)
print(p_trend_combined)

# 5f. Trend summary bar chart (STL trend-component Theil-Sen slopes) ------------

# Reshape trend_stats for side-by-side obs vs sim bars
trend_long <- trend_stats %>%
  select(species, treatment, ts_trend_slope_obs, ts_trend_slope_sim,
         ts_trend_p_obs, ts_trend_p_sim) %>%
  pivot_longer(
    cols      = c(ts_trend_slope_obs, ts_trend_slope_sim),
    names_to  = "source",
    values_to = "slope"
  ) %>%
  mutate(
    source = ifelse(grepl("obs", source), "observed", "simulated"),
    p_val  = ifelse(source == "observed", ts_trend_p_obs, ts_trend_p_sim),
    sig    = sig_star(p_val)
  )

# Compute label offset from the data (avoid reaching outside via trend_long$)
y_range <- max(abs(trend_long$slope), na.rm = TRUE)
y_off   <- if (is.finite(y_range)) y_range * 0.08 else 0.01

p_trend <- ggplot(trend_long, aes(x = species, y = slope, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge(0.9),
           color = "white", width = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_text(
    aes(label = sig, y = slope + sign(slope) * y_off),
    position = position_dodge(0.9), size = 4, vjust = 0.5
  ) +
  facet_wrap(~ source, ncol = 2,
    labeller = labeller(source = c("observed" = "Observed TWD", "simulated" = "Simulated TWD"))) +
  scale_fill_manual(
    values = c("control" = "#1f77b4", "drought" = "#d62728"),
    name   = "Treatment"
  ) +
  labs(
    title    = "STL Trend-Component Slope (Theil-Sen)",
    subtitle = "* p<0.05  ** p<0.01  *** p<0.001 (Mann-Kendall on STL trend)",
    x        = "Species",
    y        = "Slope (Z-score / month)"
  ) +
  base_theme +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, size = 11))

ggsave("Figures/lpj_guess_stem_storage/TWD/trend_summary.png",
       p_trend, width = 10, height = 6, dpi = 300)
print(p_trend)

# 5g. Kendall-Tau correlation heatmap ------------------------------------------

p_kt <- ggplot(trend_stats, aes(x = "Kendall τ", y = species, fill = kendall_tau)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(
    aes(label = paste0(round(kendall_tau, 3), sig_star(kendall_p))),
    size = 4.5, color = "black"
  ) +
  facet_wrap(~ treatment, ncol = 1) +
  scale_fill_gradient2(
    low  = "#d73027", mid = "white", high = "#4575b4",
    midpoint = 0, limits = c(-1, 1), na.value = "grey80",
    name = "Kendall τ"
  ) +
  labs(
    title    = "Kendall-τ Cross-Correlation: Observed vs Simulated TWD",
    subtitle = expression(paste("Monthly means | * p<0.05  ** p<0.01  *** p<0.001")),
    x = NULL, y = "Species"
  ) +
  base_theme +
  theme(
    legend.position   = "right",
    axis.text.x       = element_text(angle = 0, size = 11),
    axis.text.y       = element_text(size = 12),
    legend.key.height = unit(1.5, "cm")
  )

ggsave("Figures/lpj_guess_stem_storage/TWD/kendall_correlation.png",
       p_kt, width = 8, height = 6, dpi = 300)
print(p_kt)

# 5h. Robust Theil-Sen scatter plots (replace OLS with TS regression) ----------

# Add Theil-Sen fit lines to the monthly scatter (Section 4 scatter_monthly_plot)
# Compute TS slopes per species×treatment for the scatter
ts_scatter_stats <- combined_monthly_std %>%
  filter(!is.na(std_twd_mean), !is.na(std_twd_model)) %>%
  group_by(species, treatment) %>%
  summarise(
    ts_slope     = if (n() >= 3) {
                     x <- std_twd_mean; y <- std_twd_model
                     unname(coef(mblm(y ~ x))[2])
                   } else NA_real_,
    ts_intercept = if (n() >= 3) {
                     x <- std_twd_mean; y <- std_twd_model
                     unname(coef(mblm(y ~ x))[1])
                   } else NA_real_,
    .groups      = "drop"
  ) %>%
  mutate(
    label = paste0("TS slope: ", round(ts_slope, 3))
  )

scatter_monthly_ts <- ggplot(combined_monthly_std,
                              aes(x = std_twd_mean, y = std_twd_model, color = species)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.6) +
  geom_smooth(method = "lm", color = "grey60", se = FALSE, linewidth = 0.4) +
  # Theil-Sen robust regression line (thicker, black)
  geom_abline(data = ts_scatter_stats,
              aes(slope = ts_slope, intercept = ts_intercept),
              color = "black", linewidth = 0.7) +
  geom_text(data = stats_monthly %>%
              left_join(ts_scatter_stats %>% select(species, treatment, label_ts = label),
                        by = c("species", "treatment")) %>%
              mutate(label = paste0(label, "\n", label_ts)),
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.05, size = 2.5, color = "black",
            inherit.aes = FALSE) +
  facet_grid(treatment ~ species) +
  coord_fixed(ratio = 1, xlim = c(-3, 3), ylim = c(-3, 3)) +
  scale_color_manual(values = species_colors) +
  labs(
    title    = "Standardized TWD Comparison (Monthly Means)",
    subtitle = "dashed = 1:1 line | thin grey = OLS | thick black = Theil-Sen",
    x        = "Observed (Z)",
    y        = "Simulated (Z)"
  ) +
  base_theme

ggsave("Figures/lpj_guess_stem_storage/TWD/scatter_monthly_theilsen.png",
       scatter_monthly_ts, width = 10, height = 6, dpi = 300)
print(scatter_monthly_ts)

# 5i. Export statistics CSV ----------------------------------------------------
write.csv(trend_stats,
  "Figures/lpj_guess_stem_storage/TWD/trend_correlation_stats.csv",
  row.names = FALSE)

cat("\n=== Trend & Correlation Statistics ===\n")
print(trend_stats[, c("species", "treatment", "n_months",
                       "ts_slope_obs", "ts_p_obs",
                       "ts_slope_sim", "ts_p_sim",
                       "ts_trend_slope_obs", "ts_trend_p_obs",
                       "ts_trend_slope_sim", "ts_trend_p_sim",
                       "kendall_tau", "kendall_p")],
      row.names = FALSE)

cat("\nAll STL decomposition, trend, and correlation outputs saved to:\n",
    "  Figures/lpj_guess_stem_storage/TWD/\n")



