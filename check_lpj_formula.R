setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LPJ_GUESS_HYD")

library(tidyverse)
library(deSolve)

# --- 0. Define the base_theme ---
base_theme <- theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text       = element_text(color = "black", size = 9), # Smaller for many items
    legend.title      = element_text(size = 12, face = "bold"),
    legend.position   = "right", 
    plot.title        = element_text(hjust = 0.5, size = 18, color = "black", face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5, size = 14),
    axis.title        = element_text(size = 16),
    axis.text.x       = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y       = element_text(angle = 0, hjust = 0.5, size = 12),
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_line(color = "grey92", linewidth = 0.25)
  )

# --- 1. Setup Directories & Global Params ---
dir.create("Figures/lpj_formula", recursive = TRUE, showWarnings = FALSE)
times <- seq(0, 30, by = 0.1)
psi_soil_0 <- -0.2
get_psi_soil <- function(t) { psi_soil_0 - 0.08 * t }

# Helper: Create a rainbow palette for any number of lines
get_rainbow_palette <- function(n) {
  colorRampPalette(c("#ff7f7f", "#ffbf7f", "#ffff7f", "#7fff7f", "#7fffff", "#7f7fff", "#ff7fff"))(n)
}

# --- 2. Sensitivity Analysis: Lambda ---
lambda_vals <- seq(-0.5, 1, by = 0.1)
delta_psi_const <- 0.5 

hydro_ode_lambda <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    psi_s <- get_psi_soil(t)
    dPsiL <- psi_s * (1 - lambda) - psi_L - delta_psi_const
    return(list(dPsiL))
  })
}

results_lambda <- map_df(lambda_vals, function(l) {
  out <- ode(y = c(psi_L = -0.8), times = times, func = hydro_ode_lambda, parms = c(lambda = l))
  as.data.frame(out) %>% mutate(lambda = l, psi_soil = get_psi_soil(time))
}) %>% 
  mutate(lambda_f = factor(round(lambda, 1))) # Factor for individual legend lines

# Plots 1 & 2
lambda_palette <- get_rainbow_palette(length(unique(results_lambda$lambda_f)))

p1 <- ggplot(results_lambda, aes(x = time, y = psi_L, color = lambda_f)) +
  geom_line(aes(y = psi_soil), color = "black", linetype = "dashed", linewidth = 1.2) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = lambda_palette) +
  labs(
    title = "soil (dashed black) vs leaf (colored)",
    subtitle = expression(paste("colored lines vary ", lambda, " from -0.5 to 1.0")),
    x = "time", y = "water potential (MPa)", color = expression(lambda)
  ) + base_theme

p2 <- ggplot(results_lambda, aes(x = psi_soil, y = psi_L, color = lambda_f)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = lambda_palette) +
  labs(
    title = "leaf sensitivity to soil potential",
    subtitle = expression(paste("curves ordered from ", lambda, " = -0.5 to ", lambda, " = 1")),
    x = expression(Psi["   soil"] ~ "(MPa)"), y = expression(Psi["   leaf"] ~ "(MPa)"),
    color = expression(lambda)
  ) + base_theme

# --- 3. Sensitivity Analysis: DeltaPsiMax ---
lambda_fixed <- 0.4 
delta_psi_vals <- seq(0.2, 4.5, by = 0.1)

hydro_ode_delta <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    psi_s <- get_psi_soil(t)
    dPsiL <- psi_s * (1 - lambda_fixed) - psi_L - d_max
    return(list(dPsiL))
  })
}

results_delta <- map_df(delta_psi_vals, function(d) {
  out <- ode(y = c(psi_L = -0.8), times = times, func = hydro_ode_delta, parms = c(d_max = d))
  as.data.frame(out) %>% mutate(d_max = d, psi_soil = get_psi_soil(time))
}) %>% 
  mutate(d_max_f = factor(round(d_max, 1)))

# Plots 3 & 4
delta_palette <- get_rainbow_palette(length(unique(results_delta$d_max_f)))
leg_title_delta <- expression(Delta * " " * Psi["   max"]) # Added " " for spacing

p3 <- ggplot(results_delta, aes(x = time, y = psi_L, color = d_max_f)) +
  geom_line(aes(y = psi_soil), color = "black", linetype = "dashed", linewidth = 1.2) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = delta_palette) +
  guides(color = guide_legend(ncol = 3)) + # Multi-column for many lines
  labs(
    title = paste0("soil (dashed) vs leaf (solid), λ = ", lambda_fixed),
    subtitle = expression(paste("colored lines vary ", Delta, " ", Psi["   max"], " from 0.2 to 4.5 MPa")),
    x = "time", y = "water potential (MPa)", color = leg_title_delta
  ) + base_theme

p4 <- ggplot(results_delta, aes(x = psi_soil, y = psi_L, color = d_max_f)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = delta_palette) +
  guides(color = guide_legend(ncol = 3)) +
  labs(
    title = expression(paste("sensitivity: ", Psi["   leaf"], " vs ", Psi["   soil"], " (", lambda, " = 0.4)")),
    subtitle = expression(paste("curves ordered by ", Delta, " ", Psi["   max"], " (0.2 -> 4.5 MPa)")),
    x = expression(Psi["   soil"] ~ "(MPa)"), y = expression(Psi["   leaf"] ~ "(MPa)"),
    color = leg_title_delta
  ) + base_theme

# --- 4. Save Figures ---
ggsave("Figures/lpj_formula/psi_time_lambda.png", p1, width = 10, height = 8, dpi = 300)
ggsave("Figures/lpj_formula/psi_sensitivity_lambda.png", p2, width = 10, height = 8, dpi = 300)
ggsave("Figures/lpj_formula/psi_time_delta.png", p3, width = 12, height = 8, dpi = 300)
ggsave("Figures/lpj_formula/psi_sensitivity_delta.png", p4, width = 12, height = 8, dpi = 300)

print("Process complete. All figures updated with rainbow color and fixed symbols.")